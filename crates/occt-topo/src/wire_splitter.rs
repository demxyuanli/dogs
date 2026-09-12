//! Split a set of edges into closed wires (loops).
//!
//! Ports `BOPAlgo_WireSplitter` (TKBO). `perform()` runs MakeWire on a
//! regular block (even vertex degree after geometric dedup) and
//! `SplitBlock` / `Path` (`wire_splitter_block`) on an irregular block.
//!
//! Edges are matched by the coordinates of their endpoints (order
//! independent for identity, order aware for the chaining), taken from
//! `topo_tools_full::{edge_vertices, vertex_position}`.

use std::collections::{HashMap, HashSet};

use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools_full::{edge_vertices, edges_of, vertex_position};

/// Vertex identity: 3D coordinates quantized onto a `1e-6` grid.
type VKey = (i64, i64, i64);
/// Ordered endpoint pair of an edge (first child vertex -> second).
type EndKey = (VKey, VKey);

/// Grid size used to quantize vertex coordinates for identity tests.
const VERTEX_TOL: f64 = 1e-6;

/// Quantized identity key of a vertex (by its registered position).
fn vertex_key(v: &Vertex) -> VKey {
    let p = vertex_position(v);
    (
        (p.x() / VERTEX_TOL).round() as i64,
        (p.y() / VERTEX_TOL).round() as i64,
        (p.z() / VERTEX_TOL).round() as i64,
    )
}

/// Ordered endpoint keys of an edge (child order: first -> second).
fn edge_end_keys(e: &Edge) -> EndKey {
    let (a, b) = edge_vertices(e);
    let (Some(a), Some(b)) = (a, b) else {
        return ((0, 0, 0), (0, 0, 0));
    };
    (vertex_key(&a), vertex_key(&b))
}

/// The far end of edge `i` (given we stand at `current`).
fn other_endpoint(ends: &[(VKey, VKey)], i: usize, current: VKey) -> VKey {
    let (a, b) = ends[i];
    if a == current {
        b
    } else {
        a
    }
}

/// Choose the next edge to take when standing at `current` with several
/// unused candidates.
///
/// Preference order (the "select the loop back" rule):
/// 1. the candidate that immediately closes the loop back to `target`;
/// 2. a candidate whose far end is a vertex already on the current path
///    (the earliest such visit — the tightest loop);
/// 3. otherwise the candidate whose far end is closest to `target` (a greedy
///    heuristic that keeps the walk hugging the loop boundary).
fn select_next_edge(
    cands: &[usize],
    ends: &[(VKey, VKey)],
    current: VKey,
    target: VKey,
    path_verts: &[VKey],
) -> usize {
    for &i in cands {
        if other_endpoint(ends, i, current) == target {
            return i;
        }
    }
    let mut best: Option<(usize, usize)> = None;
    for &i in cands {
        let o = other_endpoint(ends, i, current);
        if let Some(k) = path_verts.iter().position(|&v| v == o) {
            if best.map_or(true, |(_, bk)| k < bk) {
                best = Some((i, k));
            }
        }
    }
    if let Some((i, _)) = best {
        return i;
    }
    let dist = |i: usize| -> f64 {
        let o = other_endpoint(ends, i, current);
        let dx = (o.0 - target.0) as f64;
        let dy = (o.1 - target.1) as f64;
        let dz = (o.2 - target.2) as f64;
        dx * dx + dy * dy + dz * dz
    };
    let mut it = cands.iter().cloned();
    let mut sel = it.next().unwrap();
    let mut sd = dist(sel);
    for i in it {
        let d = dist(i);
        if d < sd {
            sd = d;
            sel = i;
        }
    }
    sel
}

/// Deduplicate the edges by geometric endpoint pair: the on-face edges arrive
/// FORWARD and REVERSED, which are the same segment and would otherwise be seen
/// by the greedy walk as an immediate loop-back (every loop cut into single
/// edges). One representative per segment is kept.
fn dedup_edges(edges: &[Edge]) -> Vec<Edge> {
    let mut dedup: Vec<Edge> = Vec::new();
    let mut seen: HashSet<(VKey, VKey)> = HashSet::new();
    for e in edges {
        let (a, b) = edge_vertices(e);
        let (Some(av), Some(bv)) = (a, b) else { continue };
        let key = (vertex_key(&av).min(vertex_key(&bv)), vertex_key(&av).max(vertex_key(&bv)));
        if seen.insert(key) {
            dedup.push(e.clone());
        }
    }
    dedup
}

/// Whether every vertex has even degree when each geometric segment is counted
/// once (OCCT `IsRegular`). A regular block is a union of closed loops (a face
/// whose only on-face edges form closed section rings); an irregular block has
/// odd-degree vertices (a face cut by open section segments).
fn is_regular_block(edges: &[Edge]) -> bool {
    let mut deg: HashMap<VKey, usize> = HashMap::new();
    let mut seen: HashSet<(VKey, VKey)> = HashSet::new();
    for e in edges {
        let (a, b) = edge_vertices(e);
        let (Some(av), Some(bv)) = (a, b) else { continue };
        let key = (vertex_key(&av).min(vertex_key(&bv)), vertex_key(&av).max(vertex_key(&bv)));
        if !seen.insert(key) {
            continue;
        }
        *deg.entry(vertex_key(&av)).or_default() += 1;
        *deg.entry(vertex_key(&bv)).or_default() += 1;
    }
    deg.values().all(|&d| d % 2 == 0)
}

/// A set of edges together with the face they belong to.
///
/// Mirrors `BOPAlgo_WireEdgeSet`: the edges are the (intersection) result of
/// one face's splitting, `face` records the owner face for reference.
#[derive(Debug, Clone, Default)]
pub struct WireEdgeSet {
    face: Option<Face>,
    edges: Vec<Edge>,
}

impl WireEdgeSet {
    /// Empty edge set.
    pub fn new() -> Self {
        Self {
            face: None,
            edges: Vec::new(),
        }
    }

    /// Set the owning face.
    pub fn set_face(&mut self, face: Face) {
        self.face = Some(face);
    }

    /// The owning face, if any.
    pub fn face(&self) -> Option<&Face> {
        self.face.as_ref()
    }

    /// Add an edge.
    pub fn add_edge(&mut self, e: Edge) {
        self.edges.push(e);
    }

    /// Add several edges.
    pub fn add_edges(&mut self, es: &[Edge]) {
        self.edges.extend(es.iter().cloned());
    }

    /// The edges of the set.
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// Whether the set holds no edges.
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }

    /// Remove all edges.
    pub fn clear(&mut self) {
        self.edges.clear();
    }
}

/// Builds closed wires (loops) from a set of edges.
///
/// Mirrors `BOPAlgo_WireSplitter`. `perform()` consumes the edges of the
/// [`WireEdgeSet`], chaining each extracted loop into a closed `Wire`; the
/// resulting wires are available through [`WireSplitter::wires`].
pub struct WireSplitter {
    wes: WireEdgeSet,
    wires: Vec<TopoShape>,
}

impl Default for WireSplitter {
    fn default() -> Self {
        Self::new()
    }
}

impl WireSplitter {
    /// Empty splitter.
    pub fn new() -> Self {
        Self {
            wes: WireEdgeSet::new(),
            wires: Vec::new(),
        }
    }

    /// Set the wire edge set to process.
    pub fn set_wes(&mut self, wes: WireEdgeSet) {
        self.wes = wes;
    }

    /// The wire edge set.
    pub fn wes(&self) -> &WireEdgeSet {
        &self.wes
    }

    /// Mutable access to the wire edge set.
    pub fn wes_mut(&mut self) -> &mut WireEdgeSet {
        &mut self.wes
    }

    /// The closed wires built by `perform`.
    pub fn wires(&self) -> &[TopoShape] {
        &self.wires
    }

    /// Split `edges` into closed wires, consuming one loop per call to
    /// [`WireSplitter::make_wire`].
    fn split_edges(mut pool: Vec<Edge>) -> Result<Vec<TopoShape>, String> {
        let mut wires = Vec::new();
        while !pool.is_empty() {
            wires.push(Self::make_wire(&mut pool)?);
        }
        Ok(wires)
    }

    /// Run the splitting: chain the edges of the [`WireEdgeSet`] into closed
    /// wires.
    ///
    /// The on-face edges arrive FORWARD and REVERSED. A *regular* block — every
    /// geometric segment counted once gives each vertex even degree, i.e. a
    /// face whose only on-face edges form closed section rings — is chained
    /// after deduplicating the reversed copies: the greedy walk would otherwise
    /// see the reversed copy as an immediate loop-back and cut every ring into
    /// single edges. An irregular block (odd-degree vertices from open section
    /// segments) is chained directly; on failure a planar face falls back to
    /// the `SplitBlock` angle walk.
    pub fn perform(&mut self) -> Result<(), String> {
        self.wires.clear();
        if self.wes.edges.is_empty() {
            return Err("BOPAlgo_WireSplitter::perform: no input edges".to_string());
        }
        let regular = is_regular_block(&self.wes.edges);
        // OCCT MakeWires: a regular connexity block is MakeWire; an irregular
        // block (odd-degree vertices) always runs SplitBlock when a face is set.
        if !regular {
            if let Some(face) = self.wes.face().cloned() {
                self.wires = crate::wire_splitter_block::split_block(&face, &self.wes.edges)?;
                return Ok(());
            }
        }
        let edges = if regular {
            dedup_edges(&self.wes.edges)
        } else {
            self.wes.edges.clone()
        };
        match Self::split_edges(edges) {
            Ok(w) => {
                self.wires = w;
                Ok(())
            }
            Err(greedy_err) => {
                if let Some(face) = self.wes.face().cloned() {
                    self.wires = crate::wire_splitter_block::split_block(&face, &self.wes.edges)?;
                    return Ok(());
                }
                Err(greedy_err)
            }
        }
    }

    /// Chain the given edges by vertex adjacency into one closed wire.
    ///
    /// The edges of `edges` that form the first closed loop are consumed
    /// (removed from the slice); any remaining edges stay for a subsequent
    /// call. At a branch vertex the edge that leads back into the loop is
    /// preferred (see [`select_next_edge`]).
    ///
    /// Returns an error (never panics) when the edges cannot form a closed
    /// loop — e.g. a single open edge, or an open chain.
    pub fn make_wire(edges: &mut Vec<Edge>) -> Result<TopoShape, String> {
        if edges.is_empty() {
            return Err("BOPAlgo_WireSplitter::make_wire: no edges".to_string());
        }
        // A single edge: a closed edge (both endpoints coincide) is a valid
        // one-edge wire; an open edge cannot form a closed wire.
        if edges.len() == 1 {
            let e = edges[0].clone();
            let (a, b) = edge_vertices(&e);
            let (Some(a), Some(b)) = (a, b) else {
                return Err("make_wire: edge has no endpoints".to_string());
            };
            if vertex_key(&a) == vertex_key(&b) {
                edges.clear();
                let w = TopoBuilder::new().make_wire(&[e]);
                w.0.set_closed(true);
                return Ok(w.0);
            }
            return Err("make_wire: single open edge cannot form a closed wire".to_string());
        }
        // Ordered endpoint keys and vertex -> incident edge adjacency.
        let ends: Vec<EndKey> = edges.iter().map(edge_end_keys).collect();
        let mut adj: HashMap<VKey, Vec<usize>> = HashMap::new();
        for (i, &(a, b)) in ends.iter().enumerate() {
            adj.entry(a).or_default().push(i);
            if a != b {
                adj.entry(b).or_default().push(i);
            }
        }
        // Greedy walk starting from the first edge.
        let (s0, e0) = ends[0];
        let target = s0;
        let mut current = e0;
        let mut order: Vec<usize> = vec![0];
        let mut used: HashSet<usize> = HashSet::new();
        used.insert(0);
        let mut path_verts: Vec<VKey> = vec![s0];
        let mut prev_edge = 0usize;
        loop {
            if current == target {
                break;
            }
            let mut cands: Vec<usize> = adj
                .get(&current)
                .map(|v| {
                    v.iter()
                        .cloned()
                        .filter(|&i| i != prev_edge && !used.contains(&i))
                        .collect()
                })
                .unwrap_or_default();
            cands.sort_unstable();
            cands.dedup();
            if cands.is_empty() {
                return Err("make_wire: open chain, no outgoing edge at a vertex".to_string());
            }
            let next = select_next_edge(&cands, &ends, current, target, &path_verts);
            let o = other_endpoint(&ends, next, current);
            // Loop-back: the next edge lands on a vertex already on the path
            // (other than the start) — close the tight sub-loop as the wire.
            if o != target {
                if let Some(k) = path_verts.iter().position(|&v| v == o) {
                    let mut loop_edges: Vec<Edge> =
                        order[k..].iter().map(|&i| edges[i].clone()).collect();
                    loop_edges.push(edges[next].clone());
                    // Only the tight sub-loop (order[k..] + next) is consumed;
                    // the leading path (order[0..k]) belongs to the *outer*
                    // loop that is still being walked and must stay in the pool.
                    let mut consumed: HashSet<usize> = order[k..].iter().cloned().collect();
                    consumed.insert(next);
                    let remaining: Vec<Edge> = edges
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| !consumed.contains(i))
                        .map(|(_, e)| e.clone())
                        .collect();
                    *edges = remaining;
                    let w = TopoBuilder::new().make_wire(&loop_edges);
                    w.0.set_closed(true);
                    return Ok(w.0);
                }
            }
            order.push(next);
            used.insert(next);
            path_verts.push(current);
            current = o;
            prev_edge = next;
        }
        // Closed back at the start vertex: `order` is the whole loop.
        let wire_edges: Vec<Edge> = order.iter().map(|&i| edges[i].clone()).collect();
        let used_set: HashSet<usize> = order.iter().cloned().collect();
        let remaining: Vec<Edge> = edges
            .iter()
            .enumerate()
            .filter(|(i, _)| !used_set.contains(i))
            .map(|(_, e)| e.clone())
            .collect();
        *edges = remaining;
        let w = TopoBuilder::new().make_wire(&wire_edges);
        w.0.set_closed(true);
        Ok(w.0)
    }

    /// `BOPAlgo_WireSplitter::SplitBlock`.
    pub fn split_block(face: &Face, edges: Vec<Edge>) -> Result<Vec<TopoShape>, String> {
        crate::wire_splitter_block::split_block(face, &edges)
    }

    /// Build the closed wires of `face` — the outer loop plus, for a face
    /// with holes, one loop per hole.
    pub fn make_wires(face: &Face) -> Result<Vec<TopoShape>, String> {
        let mut wes = WireEdgeSet::new();
        wes.set_face(face.clone());
        wes.add_edges(&edges_of(&face.0));
        let mut ws = Self::new();
        ws.set_wes(wes);
        ws.perform()?;
        Ok(ws.wires().to_vec())
    }
}

/// Verify that a wire is topologically closed: every vertex of its edge set
/// has even degree and the vertices are all connected. (Equivalent to the
/// existence of a closed Eulerian walk over the wire edges.)
fn wire_chains_closed(wire: &TopoShape) -> bool {
    let edges = edges_of(wire);
    if edges.is_empty() {
        return false;
    }
    let mut deg: HashMap<VKey, usize> = HashMap::new();
    let mut adj: HashMap<VKey, Vec<VKey>> = HashMap::new();
    for e in &edges {
        let (a, b) = edge_end_keys(e);
        *deg.entry(a).or_default() += 1;
        if a != b {
            *deg.entry(b).or_default() += 1;
            adj.entry(a).or_default().push(b);
            adj.entry(b).or_default().push(a);
        }
    }
    if !deg.values().all(|&d| d % 2 == 0) {
        return false;
    }
    // Connectivity over the vertices.
    let start = *deg.keys().next().unwrap();
    let mut seen = HashSet::new();
    let mut stack = vec![start];
    while let Some(v) = stack.pop() {
        if !seen.insert(v) {
            continue;
        }
        if let Some(ns) = adj.get(&v) {
            for n in ns {
                if !seen.contains(n) {
                    stack.push(*n);
                }
            }
        }
    }
    seen.len() == deg.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abs::Orientation;
    use crate::brep_extrema::test_box::unit_box;
    use crate::shape::Wire;
    use crate::topo_tools_full::edges_of;
    use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt};
    use std::sync::Arc;
    use occt_geom::{GeomPlane, Surface};

    fn square_edges() -> Vec<Edge> {
        let b = TopoBuilder::new();
        let p = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let e1 = b.make_edge_segment(&p[0], &p[1]);
        let e2 = b.make_edge_segment(&p[1], &p[2]);
        let e3 = b.make_edge_segment(&p[2], &p[3]);
        let e4 = b.make_edge_segment(&p[3], &p[0]);
        vec![e1, e2, e3, e4]
    }

    fn planar_face(wires: &[Wire]) -> Face {
        let b = TopoBuilder::new();
        let pln = GpPln::new(GpAx3::new(
            GpPnt::zero(),
            GpDir::new(0.0, 0.0, 1.0).unwrap(),
            &GpDir::new(1.0, 0.0, 0.0).unwrap(),
        ).unwrap());
        let surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln));
        b.make_face(surface, wires)
    }

    #[test]
    fn make_wire_chains_a_square() {
        let mut edges = square_edges();
        let w = WireSplitter::make_wire(&mut edges).expect("square closes");
        assert!(w.is_wire());
        assert!(w.closed());
        assert!(edges.is_empty(), "all 4 edges consumed");
        assert!(wire_chains_closed(&w), "edges chain end-to-end and close");
    }

    #[test]
    fn single_open_edge_is_err_not_panic() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(1.0, 0.0, 0.0));
        let mut edges = vec![e];
        let r = WireSplitter::make_wire(&mut edges);
        assert!(r.is_err(), "single open edge cannot close");
        assert_eq!(edges.len(), 1, "edge not consumed on error");
    }

    #[test]
    fn single_closed_edge_is_a_valid_wire() {
        let b = TopoBuilder::new();
        // Both endpoint vertices coincide -> a closed (degenerate) edge.
        let p = GpPnt::new(1.0, 0.0, 0.0);
        let dir = GpDir::new(1.0, 0.0, 0.0).unwrap();
        let mut e = b.make_edge(
            Arc::new(occt_geom::GeomLine::new(
                occt_core::gp::GpLin::from_pnt_dir(p, dir),
            )),
            0.0,
            1.0,
        );
        let v = b.make_vertex(p, 0.0);
        b.add(&mut e.0, &v.0);
        b.add(&mut e.0, &v.0); // same vertex twice -> endpoints coincide
        let mut edges = vec![e];
        let w = WireSplitter::make_wire(&mut edges).expect("closed edge makes a wire");
        assert!(w.is_wire());
        assert!(edges.is_empty());
    }

    #[test]
    fn empty_make_wire_is_err() {
        let mut edges: Vec<Edge> = Vec::new();
        assert!(WireSplitter::make_wire(&mut edges).is_err());
    }

    #[test]
    fn box_faces_rechain_into_four_edge_loops() {
        let boxed = unit_box();
        for f in &boxed.faces {
            assert_eq!(edges_of(&f.0).len(), 4, "each box face has 4 edges");
            let wires = WireSplitter::make_wires(f).expect("face wires build");
            assert_eq!(wires.len(), 1, "one closed loop per box face");
            assert!(wires[0].is_wire());
            assert!(wires[0].closed());
            assert_eq!(edges_of(&wires[0]).len(), 4, "loop has the 4 edges");
            assert!(wire_chains_closed(&wires[0]));
        }
    }

    #[test]
    fn face_with_hole_yields_outer_and_inner_loop() {
        let b = TopoBuilder::new();
        // Outer square 0..1, inner square 0.25..0.75 (hole).
        let outer_p = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let inner_p = [
            GpPnt::new(0.25, 0.25, 0.0),
            GpPnt::new(0.75, 0.25, 0.0),
            GpPnt::new(0.75, 0.75, 0.0),
            GpPnt::new(0.25, 0.75, 0.0),
        ];
        let mut outer_edges = Vec::new();
        for i in 0..4 {
            outer_edges.push(b.make_edge_segment(&outer_p[i], &outer_p[(i + 1) % 4]));
        }
        let mut inner_edges = Vec::new();
        for i in 0..4 {
            inner_edges.push(b.make_edge_segment(&inner_p[i], &inner_p[(i + 1) % 4]));
        }
        let outer_wire = b.make_wire(&outer_edges);
        let inner_wire = b.make_wire(&inner_edges);
        let face = planar_face(&[outer_wire, inner_wire]);

        let wires = WireSplitter::make_wires(&face).expect("outer + hole loops");
        assert_eq!(wires.len(), 2, "outer loop and hole loop");
        for w in &wires {
            assert_eq!(edges_of(w).len(), 4);
            assert!(w.closed());
            assert!(wire_chains_closed(w));
        }

        // Covering-style WES (`BuildSplitFaces` 1.1 + 1.2): outer boundary
        // plus IN hole edges FORWARD and REVERSED. A regular block must
        // still walk the annulus (outer + inner) after dedup.
        let mut wes = WireEdgeSet::new();
        wes.set_face(face.clone());
        wes.add_edges(&outer_edges);
        for e in &inner_edges {
            wes.add_edges(&[e.clone()]);
            let mut rev = e.clone();
            rev.0.set_orientation(Orientation::Reversed);
            wes.add_edges(&[rev]);
        }
        let mut ws = WireSplitter::new();
        ws.set_wes(wes);
        ws.perform().expect("IN + boundary walks the annulus");
        assert_eq!(ws.wires().len(), 2, "outer loop and hole loop from IN+boundary");
        for w in ws.wires() {
            assert_eq!(edges_of(w).len(), 4);
            assert!(w.closed());
            assert!(wire_chains_closed(w));
        }
    }

    #[test]
    fn perform_via_wire_edge_set() {
        let boxed = unit_box();
        let f = &boxed.faces[0];
        let mut wes = WireEdgeSet::new();
        wes.set_face(f.clone());
        wes.add_edges(&edges_of(&f.0));
        assert_eq!(wes.edges().len(), 4);
        let mut ws = WireSplitter::new();
        ws.set_wes(wes);
        ws.perform().expect("split succeeds");
        assert_eq!(ws.wires().len(), 1);
        assert_eq!(edges_of(&ws.wires()[0]).len(), 4);
    }

    #[test]
    fn split_block_splits_square_by_vertical_line() {
        let b = TopoBuilder::new();
        let p = [
            b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 0.0),
            b.make_vertex(GpPnt::new(0.5, 0.0, 0.0), 0.0),
            b.make_vertex(GpPnt::new(1.0, 0.0, 0.0), 0.0),
            b.make_vertex(GpPnt::new(1.0, 1.0, 0.0), 0.0),
            b.make_vertex(GpPnt::new(0.5, 1.0, 0.0), 0.0),
            b.make_vertex(GpPnt::new(0.0, 1.0, 0.0), 0.0),
        ];
        let gp = |i: usize| crate::topo_tools_full::vertex_position(&p[i]);
        let e = |i: usize, j: usize| b.make_edge_segment_with_vertices(&gp(i), &gp(j), &p[i], &p[j]);
        let boundary = [e(0, 1), e(1, 2), e(2, 3), e(3, 4), e(4, 5), e(5, 0)];
        let face = planar_face(&[b.make_wire(&boundary)]);
        let section = e(1, 4);
        let mut all = edges_of(&face.0);
        all.push(section.clone());
        let mut rev = section;
        rev.0.set_orientation(Orientation::Reversed);
        all.push(rev);
        let wires = WireSplitter::split_block(&face, all).expect("split block");
        assert_eq!(wires.len(), 2, "the square splits into two pieces, got {}", wires.len());
        for w in &wires {
            assert!(w.closed());
            assert!(wire_chains_closed(w));
        }
    }

    #[test]
    fn perform_with_empty_set_errors() {
        let mut ws = WireSplitter::new();
        ws.set_wes(WireEdgeSet::new());
        assert!(ws.perform().is_err());
        assert!(ws.wires().is_empty());
    }
}
