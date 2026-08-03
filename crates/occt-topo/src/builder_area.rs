//! Area construction — a Rust port of `BOPAlgo_BuilderArea`.
//!
//! The root class of the algorithms that build closed regions ("areas") from a
//! set of boundary shapes. The pipeline is the four-phase virtual sequence of
//! OCCT:
//!
//! ```text
//!   ShapesToAvoid → Loops → Areas → InternalShapes
//! ```
//!
//! * `PerformShapesToAvoid` strips the shapes that cannot belong to a closed
//!   boundary (isolated vertices, dangling/leaf edges).
//! * `PerformLoops` orders the surviving edges into closed loops (wires).
//! * `PerformAreas` groups the loops into region faces.
//! * `PerformInternalShapes` classifies the shapes that sit *inside* an area
//!   (e.g. a full circular edge floating in the interior of a face).
//!
//! The loop-closing step here is a self-contained graph cycle tracer over
//! vertex adjacency (OCCT delegates it to `BOPAlgo_WireSplitter`); it does not
//! depend on the parallel `wire_splitter` module.
//!
//! Vertex identity is by *position* (quantized 3-D coordinates), so edges that
//! share a point but carry distinct `Vertex` TShapes (as produced by
//! `TopoBuilder::make_edge_segment`) are still recognized as touching.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt, GpVec};
use occt_geom::GeomPlane;

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Wire};
use crate::topo_tools_full::{edge_vertices, edges_of_wire, vertex_position};

/// Quantized 3-D vertex key: two vertices whose positions agree to within
/// ~1e-9 collapse to the same key.
type VKey = (i64, i64, i64);

/// Sentinel key for an edge that does not expose two endpoints.
const BAD_KEY: VKey = (i64::MIN, i64::MIN, i64::MIN);

fn vertex_key(p: &GpPnt) -> VKey {
    const K: f64 = 1e9;
    ((p.x() * K).round() as i64, (p.y() * K).round() as i64, (p.z() * K).round() as i64)
}

/// Endpoint keys of an edge, or `BAD_KEY` when the edge has no two vertices.
fn edge_end_keys(e: &Edge) -> (VKey, VKey) {
    let (a, b) = edge_vertices(e);
    match (a, b) {
        (Some(a), Some(b)) => {
            (vertex_key(&vertex_position(&a)), vertex_key(&vertex_position(&b)))
        }
        _ => (BAD_KEY, BAD_KEY),
    }
}

/// The four-phase area-construction protocol (`BOPAlgo_BuilderArea`).
///
/// `perform` runs the canonical OCCT sequence and is the normal entry point.
pub trait AreaBuilder {
    /// Set the input shapes (edges/faces, possibly vertices).
    fn set_shapes(&mut self, shapes: &[TopoShape]);

    /// Phase 1 — decide which shapes cannot participate in a closed area.
    fn perform_shapes_to_avoid(&mut self) -> Result<(), String>;

    /// Phase 2 — close the remaining edges into loops (wires).
    fn perform_loops(&mut self) -> Result<(), String>;

    /// Phase 3 — turn the loops into area faces.
    fn perform_areas(&mut self) -> Result<(), String>;

    /// Phase 4 — classify internal shapes.
    fn perform_internal_shapes(&mut self) -> Result<(), String>;

    /// Run the full `ShapesToAvoid → Loops → Areas → InternalShapes` sequence.
    fn perform(&mut self) -> Result<(), String> {
        self.perform_shapes_to_avoid()?;
        self.perform_loops()?;
        self.perform_areas()?;
        self.perform_internal_shapes()
    }

    /// The area faces produced by the pipeline.
    fn areas(&self) -> &[TopoShape];
}

/// Shared state of an area builder (`BOPAlgo_BuilderArea` members).
///
/// `shapes` is the input set; `loops` holds the closed boundary wires found by
/// `PerformLoops`; `avoid` holds the shapes excluded from the boundary;
/// `internal` holds the wires reconstructed from the shapes that lie inside an
/// area; `areas` is the final face list.
#[derive(Debug, Default)]
pub struct AreaBuilderBase {
    /// Input shapes (edges, and possibly vertices).
    pub shapes: Vec<TopoShape>,
    /// Closed boundary wires found by `PerformLoops`.
    pub loops: Vec<TopoShape>,
    /// Area faces produced by `PerformAreas`.
    pub areas: Vec<TopoShape>,
    /// Shapes excluded from the boundary by `PerformShapesToAvoid`.
    pub avoid: Vec<TopoShape>,
    /// Internal wires reconstructed from interior shapes.
    pub internal: Vec<TopoShape>,
    /// Non-fatal diagnostics.
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
    /// When true, internal shapes are *not* merged back into the result
    /// (`BOPAlgo_BuilderArea::SetAvoidInternalShapes`).
    pub avoid_internal_shapes: bool,
}

impl AreaBuilderBase {
    /// Whether `s` is currently marked to be avoided (same `TShape` identity).
    pub fn is_avoided(&self, s: &TopoShape) -> bool {
        self.avoid.iter().any(|a| a.same_tshape(s))
    }

    /// The input edges (shapes filtered down to `EDGE` type).
    fn input_edges(&self) -> Vec<Edge> {
        self.shapes.iter().filter(|s| s.is_edge()).map(|s| Edge(s.clone())).collect()
    }
}

impl AreaBuilder for AreaBuilderBase {
    fn set_shapes(&mut self, shapes: &[TopoShape]) {
        self.shapes = shapes.to_vec();
    }

    fn perform_shapes_to_avoid(&mut self) -> Result<(), String> {
        self.avoid.clear();
        // Non-edge shapes (vertices, faces, ...) cannot form a closed boundary.
        // In OCCT `myShapes` holds only edges for `BuilderFace`; the generic
        // base additionally drops stray vertices here.
        let mut edges = Vec::new();
        for s in &self.shapes {
            if s.is_edge() {
                edges.push(Edge(s.clone()));
            } else {
                self.avoid.push(s.clone());
            }
        }
        self.avoid.extend(shapes_to_avoid(&edges));
        Ok(())
    }

    fn perform_loops(&mut self) -> Result<(), String> {
        self.loops.clear();
        let usable: Vec<Edge> =
            self.input_edges().into_iter().filter(|e| !self.is_avoided(&e.0)).collect();
        let loops = collect_closed_loops(&usable);
        // Edges that could not be closed into a loop are set aside, exactly as
        // OCCT appends unprocessed edges to `myShapesToAvoid` after the split.
        let mut used: HashSet<usize> = HashSet::new();
        for l in &loops {
            for e in l {
                used.insert(Arc::as_ptr(&e.0.tshape) as usize);
            }
        }
        for e in &usable {
            if !used.contains(&(Arc::as_ptr(&e.0.tshape) as usize)) && !self.is_avoided(&e.0) {
                self.avoid.push(e.0.clone());
            }
        }
        for l in loops {
            self.loops.push(make_wire_from_edges(&l).0);
        }
        Ok(())
    }

    fn perform_areas(&mut self) -> Result<(), String> {
        self.areas.clear();
        let loops = self.loops.clone();
        let b = TopoBuilder::new();
        for l in loops {
            let wire = Wire(l);
            let edges = edges_of_wire(&wire);
            let pln = plane_from_loop(&edges)?;
            let face = b.make_face(Arc::new(GeomPlane::new(pln)), &[wire]);
            self.areas.push(face.0);
        }
        Ok(())
    }

    fn perform_internal_shapes(&mut self) -> Result<(), String> {
        self.internal.clear();
        // Interior shapes (edges that were set aside by `PerformShapesToAvoid`
        // but that still form closed loops, e.g. a full circular edge inside a
        // face) are grouped into internal wires. Dangling edges that could not
        // close are simply left unclassified.
        let avoid_edges: Vec<Edge> =
            self.avoid.iter().filter(|s| s.is_edge()).map(|s| Edge(s.clone())).collect();
        for l in collect_closed_loops(&avoid_edges) {
            self.internal.push(make_wire_from_edges(&l).0);
        }
        Ok(())
    }

    fn areas(&self) -> &[TopoShape] {
        &self.areas
    }
}

/// Edges that cannot belong to any closed loop — the `PerformShapesToAvoid`
/// computation.
///
/// Mirrors `BOPAlgo_BuilderFace::PerformShapesToAvoid`: repeatedly strip edges
/// whose endpoint vertex is touched by a single (non-degenerated) edge, and
/// the "lasso" case where a vertex is touched twice by the same edge with
/// distinct endpoints. The result is the set of edges that cannot be part of a
/// closed boundary.
///
/// ponytail: the OCCT branch that keeps edges attached to an `INTERNAL`
/// vertex is not ported (the flat model does not retain the vertex-child
/// orientation after `make_edge_segment`); all vertices are treated as
/// non-internal. Add back when per-edge vertex orientation is representable.
pub fn shapes_to_avoid(edges: &[Edge]) -> Vec<TopoShape> {
    let n = edges.len();
    let mut avoid = vec![false; n];
    let ends: Vec<(VKey, VKey)> = edges.iter().map(edge_end_keys).collect();

    loop {
        // Incidence: which (non-avoided) edges touch each vertex.
        let mut inc: HashMap<VKey, Vec<usize>> = HashMap::new();
        for i in 0..n {
            if avoid[i] {
                continue;
            }
            let (a, b) = ends[i];
            inc.entry(a).or_default().push(i);
            if a != b {
                inc.entry(b).or_default().push(i);
            }
        }

        let mut found = false;
        for i in 0..n {
            if avoid[i] {
                continue;
            }
            let (a, b) = ends[i];
            let mut dangling = false;
            for k in [a, b] {
                if let Some(list) = inc.get(&k) {
                    if list.len() == 1 {
                        // Only this edge touches the vertex → the edge hangs
                        // off the graph and can never close a loop.
                        if !BRepTool::is_degenerated(&edges[i]) {
                            dangling = true;
                            break;
                        }
                    } else if list.len() == 2 && list[0] == list[1] {
                        // The same edge touches the vertex twice (both ends at
                        // the same point) and is not a true closed loop.
                        if a != b {
                            dangling = true;
                            break;
                        }
                    }
                }
            }
            if dangling {
                avoid[i] = true;
                found = true;
            }
        }
        if !found {
            break;
        }
    }

    edges.iter().enumerate().filter(|&(i, _)| avoid[i]).map(|(_, e)| e.0.clone()).collect()
}

/// Ordered edge lists forming closed loops.
///
/// Self-contained cycle tracer over vertex adjacency: each loop starts at an
/// unused edge and walks from vertex to vertex, always taking the next unused
/// edge incident at the current vertex. A chain whose end does not return to
/// its start is an open chain and is reported as an error.
///
/// ponytail: the tracer assumes the vertex graph has degree ≤ 2 (a face
/// boundary or a single loop). For a general shape (e.g. the full edge set of
/// a box, where corners have degree 3) OCCT classifies edges with a 2-D face
/// classifier (`BOPAlgo_WireSplitter`); this port intentionally covers the
/// face-reconstruction and disjoint-region cases.
pub fn build_loops(edges: &[Edge]) -> Result<Vec<Vec<Edge>>, String> {
    let chains = trace_chains(edges);
    let mut out = Vec::new();
    for (order, closed) in chains {
        if !closed {
            return Err(format!(
                "build_loops: edges form an open chain starting at edge {}",
                order.first().copied().unwrap_or(0)
            ));
        }
        out.push(order.into_iter().map(|i| edges[i].clone()).collect());
    }
    Ok(out)
}

/// Like [`build_loops`] but tolerant: open chains are silently dropped instead
/// of reported. Used by `PerformLoops`/`PerformInternalShapes`, which move the
/// unclosable edges to the "avoid" set rather than fail.
pub fn collect_closed_loops(edges: &[Edge]) -> Vec<Vec<Edge>> {
    trace_chains(edges)
        .into_iter()
        .filter(|(_, closed)| *closed)
        .map(|(order, _)| order.into_iter().map(|i| edges[i].clone()).collect())
        .collect()
}

/// Trace vertex-adjacency chains. Returns `(edge indices, closed?)` per chain.
fn trace_chains(edges: &[Edge]) -> Vec<(Vec<usize>, bool)> {
    let n = edges.len();
    if n == 0 {
        return Vec::new();
    }
    let ends: Vec<(VKey, VKey)> = edges.iter().map(edge_end_keys).collect();

    let mut adj: HashMap<VKey, Vec<(usize, VKey)>> = HashMap::new();
    for i in 0..n {
        let (a, b) = ends[i];
        adj.entry(a).or_default().push((i, b));
        if a != b {
            adj.entry(b).or_default().push((i, a));
        }
    }

    let mut used = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if used[start] {
            continue;
        }
        let mut order = vec![start];
        used[start] = true;
        let start_v = ends[start].0;
        let mut cur_e = start;
        let mut cur_v = ends[start].1;
        loop {
            let next = adj
                .get(&cur_v)
                .and_then(|list| list.iter().find(|&&(ei, _)| !used[ei] && ei != cur_e).copied());
            match next {
                Some((ei, other)) => {
                    order.push(ei);
                    used[ei] = true;
                    cur_e = ei;
                    cur_v = other;
                }
                None => break,
            }
        }
        out.push((order, cur_v == start_v));
    }
    out
}

/// Whether the edges form a single closed loop: every vertex has degree
/// exactly 2 (boundary of one polygon) and a connectivity trace visits all
/// edges.
pub fn edges_form_closed_loop(edges: &[Edge]) -> bool {
    if edges.len() < 3 {
        return false;
    }
    let n = edges.len();
    let mut inc: HashMap<VKey, usize> = HashMap::new();
    let mut ends: Vec<(VKey, VKey)> = Vec::with_capacity(n);
    for e in edges {
        let (a, b) = edge_end_keys(e);
        if a == BAD_KEY || b == BAD_KEY {
            return false;
        }
        ends.push((a, b));
        *inc.entry(a).or_insert(0) += 1;
        if a != b {
            *inc.entry(b).or_insert(0) += 1;
        }
    }
    if inc.values().any(|&c| c != 2) {
        return false;
    }
    // Single connected cycle: BFS from edge 0 over shared vertices.
    let mut used = vec![false; n];
    let mut stack = vec![0usize];
    used[0] = true;
    let mut count = 1;
    while let Some(i) = stack.pop() {
        let (ia, ib) = ends[i];
        for k in [ia, ib] {
            for j in 0..n {
                if used[j] {
                    continue;
                }
                let (ja, jb) = ends[j];
                if ja == k || jb == k {
                    used[j] = true;
                    count += 1;
                    stack.push(j);
                }
            }
        }
    }
    count == n
}

/// Build a wire from the ordered edges, flagged closed.
pub fn make_wire_from_edges(edges: &[Edge]) -> Wire {
    let b = TopoBuilder::new();
    let wire = b.make_wire(edges);
    wire.set_closed(true);
    wire
}

/// A plane through the loop's vertices (`BRep_Tool::Surface` of a planar
/// face). Uses the first three non-collinear vertex positions.
pub fn plane_from_loop(edges: &[Edge]) -> Result<GpPln, String> {
    let mut pts = Vec::new();
    for e in edges {
        let (a, b) = edge_vertices(e);
        for v in [a, b].into_iter().flatten() {
            pts.push(vertex_position(&v));
        }
    }
    plane_from_points(&pts)
}

/// A plane through the first three non-collinear points in `pts`.
pub fn plane_from_points(pts: &[GpPnt]) -> Result<GpPln, String> {
    let n = pts.len();
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                let v1 = GpVec::from_pnts(&pts[i], &pts[j]);
                let v2 = GpVec::from_pnts(&pts[i], &pts[k]);
                let nrm = v1.xyz().crossed(v2.xyz());
                if nrm.modulus() > 1e-12 {
                    let d = GpDir::from_xyz(&nrm).map_err(|_| "plane_from_points: degenerate normal")?;
                    let z_axis = GpDir::new(0.0, 0.0, 1.0).unwrap();
                    let x_dir =
                        if d.is_normal(&z_axis) { z_axis } else { GpDir::new(1.0, 0.0, 0.0).unwrap() };
                    let ax3 = GpAx3::new(pts[i], d, &x_dir)
                        .map_err(|e| format!("plane_from_points: {e}"))?;
                    return Ok(GpPln::new(ax3));
                }
            }
        }
    }
    Err("plane_from_points: collinear or coincident vertices".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::TopoBuilder;
    use occt_core::gp::GpPnt;

    fn square_edges(b: &TopoBuilder, pts: &[GpPnt; 4]) -> Vec<Edge> {
        (0..4).map(|i| b.make_edge_segment(&pts[i], &pts[(i + 1) % 4])).collect()
    }

    #[test]
    fn square_edges_close_into_one_loop() {
        let b = TopoBuilder::new();
        let sq = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let edges = square_edges(&b, &sq);
        let loops = build_loops(&edges).expect("square closes");
        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].len(), 4);
        assert!(edges_form_closed_loop(&loops[0]));
    }

    #[test]
    fn disjoint_squares_close_separately() {
        let b = TopoBuilder::new();
        let mut all = Vec::new();
        for z in 0..3 {
            let sq = [
                GpPnt::new(0., 0., z as f64),
                GpPnt::new(1., 0., z as f64),
                GpPnt::new(1., 1., z as f64),
                GpPnt::new(0., 1., z as f64),
            ];
            all.extend(square_edges(&b, &sq));
        }
        let loops = build_loops(&all).expect("three disjoint loops");
        assert_eq!(loops.len(), 3);
    }

    #[test]
    fn open_chain_is_rejected() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0., 0., 0.), &GpPnt::new(1., 0., 0.));
        let e2 = b.make_edge_segment(&GpPnt::new(1., 0., 0.), &GpPnt::new(1., 1., 0.));
        assert!(build_loops(&[e1, e2]).is_err());
    }

    #[test]
    fn shapes_to_avoid_strips_dangling_and_vertices() {
        let b = TopoBuilder::new();
        // A closed square: every vertex has degree 2 → nothing stripped.
        let sq = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let square = square_edges(&b, &sq);
        assert!(shapes_to_avoid(&square).is_empty(), "closed square has no dangling edge");

        // A dangling segment attached to a square corner at z = 1 is stripped.
        let dangling = b.make_edge_segment(&GpPnt::new(0., 0., 1.), &GpPnt::new(1., 0., 1.));
        let mut mixed = square.clone();
        mixed.push(dangling.clone());
        let avoid = shapes_to_avoid(&mixed);
        assert_eq!(avoid.len(), 1);
        assert!(avoid[0].same_tshape(&dangling.0));
    }

    #[test]
    fn base_pipeline_builds_areas_for_disjoint_quads() {
        let b = TopoBuilder::new();
        let mut shapes = Vec::new();
        for z in 0..6 {
            let sq = [
                GpPnt::new(0., 0., z as f64),
                GpPnt::new(1., 0., z as f64),
                GpPnt::new(1., 1., z as f64),
                GpPnt::new(0., 1., z as f64),
            ];
            for e in square_edges(&b, &sq) {
                shapes.push(e.0);
            }
        }
        let mut builder = AreaBuilderBase::default();
        builder.set_shapes(&shapes);
        builder.perform().expect("area build");
        assert_eq!(builder.loops.len(), 6);
        assert_eq!(builder.areas.len(), 6);
        assert!(builder.avoid.is_empty(), "no dangling edges in six disjoint quads");
        assert!(builder.internal.is_empty());
    }
}
