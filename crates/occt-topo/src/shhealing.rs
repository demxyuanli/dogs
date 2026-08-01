//! Shape healing — repair free edges, remove small edges, close non-closed
//! wires, weld coincident vertices and relocate vertices.
//!
//! Port of the OCCT `ShapeHealing` toolkit: `ShapeFix_Wire` (FixSmall,
//! FixClosed, FixVertexTolerance), `ShapeFix_Shape`, `ShapeFix_Edge`,
//! `ShapeFix_Vertex`, and `ShapeAnalysis_FreeBounds`.
//!
//! Every fix follows the same pipeline: collect the real topology from the
//! `TShape` child tree, build maps from the old sub-shapes to their healed
//! replacements, then rebuild a fresh `TopoShape` tree that reuses the shared
//! `TShape` handles (preserving vertex/edge sharing between faces) and
//! re-registers geometry in the `GeometryRegistry` side-table for any newly
//! created shapes. Healing never mutates the input shape — it always returns a
//! new shape, so the caller can keep the original intact.
//!
//! Geometry (`BRep_TVertex`/`BRep_TEdge`/`BRep_TFace`) lives in the
//! process-wide `GeometryRegistry` keyed by `TShape` address. Rebuilding an
//! edge or vertex therefore creates a *new* `TShape` and registers its geometry
//! afresh; edges and vertices that survive unchanged are reused by handle, so
//! their existing side-table entries remain valid.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_geom::GeomLine;

use crate::abs::ShapeType;
use crate::brep_measure::edge_length;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{
    edges_of, edges_of_wire, edge_vertices, is_same, vertex_position, vertices_of, wires_of,
};

/// Summary of the fixes applied by [`heal_shape`].
///
/// Mirrors the aggregate status `ShapeFix_Shape` exposes after running its
/// sequence of repair operations. Each counter reports how many objects of a
/// given kind were repaired; `modified` is true when at least one fix fired.
#[derive(Debug, Default, Clone)]
pub struct HealReport {
    /// Number of free boundary edges that ceased to be free after healing.
    pub free_edges_fixed: usize,
    /// Number of edges shorter than the minimum length that were removed.
    pub small_edges_removed: usize,
    /// Number of wires that received a closing edge.
    pub wires_closed: usize,
    /// Number of vertex pairs merged because they lay within tolerance.
    pub vertices_welded: usize,
    /// True when any repair operation changed the shape.
    pub modified: bool,
}

/// The address of the `TShape` inside the shared `RwLock`, used as a stable
/// identity key for a `TopoShape` (matches `GeometryRegistry`'s keying and
/// `TopExp`-style uniqueness by `Handle(TShape)`).
#[inline]
fn ptr(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

/// Direct child `TopoShape`s of `s` (the `TShape::children` list).
fn children(s: &TopoShape) -> Vec<TopoShape> {
    s.tshape
        .read()
        .unwrap()
        .children
        .iter()
        .map(|h| TopoShape::from_handle(h.clone()))
        .collect()
}

/// Position-grouping index: the index of an existing representative point
/// within `tol` of `p`, or the freshly appended index when none matches.
///
/// This is the spatial hash used to decide whether two geometric points
/// (from two distinct `Vertex` instances) should be treated as one location.
/// OCCT achieves the same effect with tolerance-based `BRep_Tool` queries;
/// here an O(n) scan over the already-grouped representatives suffices for the
/// shape sizes this port targets.  (`ponytail:` O(n²) worst case — a spatial
/// grid buckets the points when heal performance ever matters.)
fn canonical_idx(pts: &mut Vec<GpPnt>, p: &GpPnt, tol: f64) -> usize {
    for (i, q) in pts.iter().enumerate() {
        if q.distance(p) <= tol {
            return i;
        }
    }
    pts.push(*p);
    pts.len() - 1
}

/// Build a straight-line edge from `p1` to `p2`, creating its own endpoint
/// vertices (`BRepBuilderAPI_MakeEdge(P1, P2)`). Unlike
/// `TopoBuilder::make_edge_segment` this never panics on coincident endpoints:
/// a degenerate zero-length edge is returned instead (the caller decides
/// whether to drop it). Used for closing edges and for rebuilt geometry where
/// the two endpoint points may coincide within tolerance.
fn make_segment_safe(builder: &TopoBuilder, p1: &GpPnt, p2: &GpPnt) -> Edge {
    let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2))
        .unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).unwrap());
    let lin = GpLin::from_pnt_dir(*p1, dir);
    let mut e = builder.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
    let v1 = builder.make_vertex(*p1, 0.0);
    let v2 = builder.make_vertex(*p2, 0.0);
    builder.add(&mut e.0, &v1.0);
    builder.add(&mut e.0, &v2.0);
    e
}

/// Build a straight-line edge between `p1` and `p2` whose endpoint vertex
/// children are the existing `v1`/`v2` handles rather than freshly created
/// vertices. This is the key primitive for healing: it lets a rebuilt edge
/// reference the canonical (merged) vertex instances so that edges shared
/// between faces keep pointing at the same vertex.
fn make_edge_with_vertices(
    builder: &TopoBuilder,
    p1: &GpPnt,
    p2: &GpPnt,
    v1: &Vertex,
    v2: &Vertex,
) -> Edge {
    let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2))
        .unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).unwrap());
    let lin = GpLin::from_pnt_dir(*p1, dir);
    let mut e = builder.make_edge(Arc::new(GeomLine::new(lin)), 0.0, p1.distance(p2));
    builder.add(&mut e.0, &v1.0);
    builder.add(&mut e.0, &v2.0);
    e
}

/// Rebuild an edge that keeps its existing curve geometry but whose endpoint
/// vertex children are replaced by `v1`/`v2`. Used for closed-loop edges (e.g.
/// a seam vertex that got welded) where the curve does not need to change but
/// the attached vertices do.
fn edge_with_canonical_vertices(old: &Edge, v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    let reg = GeometryRegistry::global();
    let curve = reg.edge_curve(&old.0)?;
    let (first, last) = reg.edge_parameters(&old.0);
    let builder = TopoBuilder::new();
    let mut ne = builder.make_edge(curve, first, last);
    builder.add(&mut ne.0, &v1.0);
    builder.add(&mut ne.0, &v2.0);
    Some(ne)
}

/// Disjoint-set (union-find) used to collapse groups of vertices connected by
/// small edges into a single canonical vertex. `union(a, b)` keeps `a`'s root
/// as the group root, so the canonical of a group is the *first* endpoint of
/// the first small edge that mentioned it — exactly the "collapse to a point"
/// behaviour `ShapeFix_Wire::FixSmall` needs for a tiny edge at the start of a
/// wire.
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    fn find(&mut self, i: usize) -> usize {
        let mut root = i;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        // Path compression: flatten the chain so later lookups are O(1).
        let mut cur = i;
        while self.parent[cur] != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            self.parent[rb] = ra;
        }
    }
}

/// The set of maps a healing operation feeds into [`rebuild_shape`].
///
/// - `vertex_map`: old vertex `TShape` address → the vertex to use in its
///   place (the canonical representative after welding, or a relocated copy).
/// - `edge_map`: old edge `TShape` address → the edge to use in its place
///   (kept, rebuilt with canonical vertices, or replaced by a closing edge's
///   sibling). Edges shared by several faces all resolve to the same entry, so
///   sharing survives the rebuild.
/// - `edge_removed`: old edge addresses that must be dropped from their wires
///   (degenerate after collapsing, or shorter than the minimum length).
/// - `add_closing`: old wire address → a straight edge to append at the end of
///   that wire (`ShapeFix_Wire::FixClosed`).
struct HealCtx {
    builder: TopoBuilder,
    vertex_map: HashMap<usize, Vertex>,
    edge_map: HashMap<usize, Edge>,
    edge_removed: HashSet<usize>,
    add_closing: HashMap<usize, Edge>,
}

impl HealCtx {
    fn new(
        builder: TopoBuilder,
        vertex_map: HashMap<usize, Vertex>,
        edge_map: HashMap<usize, Edge>,
        edge_removed: HashSet<usize>,
        add_closing: HashMap<usize, Edge>,
    ) -> Self {
        Self {
            builder,
            vertex_map,
            edge_map,
            edge_removed,
            add_closing,
        }
    }
}

/// Rebuild a wire: keep the surviving edges (in order), drop removed ones, and
/// append the closing edge when the fixer requested one.
fn rebuild_wire(s: &TopoShape, ctx: &HealCtx) -> TopoShape {
    let wire_ptr = ptr(s);
    let mut new_wire = Wire::new();
    let mut closed_flag = s.closed();
    if ctx.add_closing.contains_key(&wire_ptr) {
        closed_flag = true;
    }
    new_wire.set_closed(closed_flag);
    for child in children(s) {
        if !child.is_edge() {
            continue;
        }
        let edge_ptr = ptr(&child);
        if ctx.edge_removed.contains(&edge_ptr) {
            continue;
        }
        let edge = ctx
            .edge_map
            .get(&edge_ptr)
            .cloned()
            .unwrap_or(Edge(child));
        ctx.builder.add_edge(&mut new_wire, &edge);
    }
    if let Some(closer) = ctx.add_closing.get(&wire_ptr) {
        ctx.builder.add_edge(&mut new_wire, closer);
    }
    new_wire.0
}

/// Rebuild a face: re-register its surface geometry, rebuild each boundary
/// wire, and drop wires that collapsed to zero edges (`ShapeFix_Shape` removes
/// a wire whose edges were all consumed).
fn rebuild_face(s: &TopoShape, ctx: &HealCtx) -> TopoShape {
    let reg = GeometryRegistry::global();
    let mut new_face = Face::new();
    if let Some(g) = reg.face_geom(s) {
        reg.set_face(&new_face.0, g);
    }
    for child in children(s) {
        if !child.is_wire() {
            continue;
        }
        let rebuilt = rebuild_wire(&child, ctx);
        let w = Wire(rebuilt);
        if !edges_of_wire(&w).is_empty() {
            ctx.builder.add_wire(&mut new_face, &w);
        }
    }
    new_face.0
}

/// Rebuild an entire shape tree from the healing maps. Every level re-creates
/// the container `TShape` (wire/face/shell/solid/compound) while sub-shapes
/// come from the maps, so shared edges and vertices stay shared.
fn rebuild_shape(s: &TopoShape, ctx: &HealCtx) -> TopoShape {
    match s.shape_type() {
        ShapeType::Vertex => ctx
            .vertex_map
            .get(&ptr(s))
            .cloned()
            .map(|v| v.0)
            .unwrap_or_else(|| s.clone()),
        ShapeType::Edge => {
            if ctx.edge_removed.contains(&ptr(s)) {
                Edge::new().0
            } else {
                ctx.edge_map
                    .get(&ptr(s))
                    .cloned()
                    .map(|e| e.0)
                    .unwrap_or_else(|| s.clone())
            }
        }
        ShapeType::Wire => rebuild_wire(s, ctx),
        ShapeType::Face => rebuild_face(s, ctx),
        ShapeType::Shell => {
            let mut shell = Shell::new();
            for child in children(s) {
                let rebuilt = rebuild_shape(&child, ctx);
                if rebuilt.is_face() {
                    ctx.builder.add_face(&mut shell, &Face(rebuilt));
                }
            }
            shell.0
        }
        ShapeType::Solid | ShapeType::CompSolid => {
            let mut solid = Solid::new();
            for child in children(s) {
                let rebuilt = rebuild_shape(&child, ctx);
                match rebuilt.shape_type() {
                    ShapeType::Shell => {
                        ctx.builder.add_shell(&mut solid, &Shell(rebuilt));
                    }
                    ShapeType::Face => {
                        // A solid may rarely hold a bare face; wrap it so the
                        // rebuilt solid keeps a well-formed shell level.
                        let shell = ctx.builder.make_shell(&[Face(rebuilt)]);
                        ctx.builder.add_shell(&mut solid, &shell);
                    }
                    _ => {}
                }
            }
            solid.0
        }
        ShapeType::Compound | ShapeType::Shape => {
            let mut comp = crate::shape::Compound::new();
            for child in children(s) {
                if child.is_edge() && ctx.edge_removed.contains(&ptr(&child)) {
                    continue;
                }
                let rebuilt = rebuild_shape(&child, ctx);
                ctx.builder.add_compound(&mut comp, &rebuilt);
            }
            comp.0
        }
    }
}

/// Whether a wire already forms a closed loop, judged with a tight (~1e-9)
/// position epsilon. This deliberately does *not* treat a wire whose endpoints
/// merely lie within a heal tolerance as closed: such a wire has a gap and is
/// exactly the case `ShapeFix_Wire::FixClosed` repairs by adding an edge.
/// A lone closed loop (one edge whose two ends are the same vertex or coincide)
/// counts as closed. Mirrors the endpoint-parity test of `BRepCheck_Shell`.
fn wire_is_closed_loop(wire: &Wire) -> bool {
    const TIGHT: f64 = 1e-9;
    let edges = edges_of_wire(wire);
    if edges.is_empty() {
        return true;
    }
    if edges.len() == 1 {
        // A lone loop: one edge whose endpoints coincide.
        let (a, b) = edge_vertices(&edges[0]);
        return match (a, b) {
            (Some(a), Some(b)) => {
                is_same(&a.0, &b.0) || vertex_position(&a).distance(&vertex_position(&b)) <= TIGHT
            }
            _ => false,
        };
    }
    let mut pts: Vec<GpPnt> = Vec::new();
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for e in &edges {
        let (a, b) = edge_vertices(e);
        let (Some(a), Some(b)) = (a, b) else {
            return false;
        };
        let ia = canonical_idx(&mut pts, &vertex_position(&a), TIGHT);
        let ib = canonical_idx(&mut pts, &vertex_position(&b), TIGHT);
        *counts.entry(ia).or_insert(0) += 1;
        *counts.entry(ib).or_insert(0) += 1;
    }
    counts.values().all(|&c| c % 2 == 0)
}

/// Whether a wire is geometrically closed: every endpoint position is used by
/// an even number of edge ends (normally exactly twice). A lone closed loop —
/// a single edge whose two ends are the same vertex or coincide within `tol` —
/// counts as closed, matching `ShapeAnalysis_FreeBounds`' handling of seam
/// loops. Positions within `tol` are treated as one location.
fn wire_geometrically_closed_with_tol(wire: &Wire, tol: f64) -> bool {
    let edges = edges_of_wire(wire);
    if edges.is_empty() {
        return true;
    }
    if edges.len() == 1 {
        // A lone loop: one edge whose endpoints coincide.
        let (a, b) = edge_vertices(&edges[0]);
        return match (a, b) {
            (Some(a), Some(b)) => {
                is_same(&a.0, &b.0) || vertex_position(&a).distance(&vertex_position(&b)) <= tol
            }
            _ => false,
        };
    }
    let mut pts: Vec<GpPnt> = Vec::new();
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for e in &edges {
        let (a, b) = edge_vertices(e);
        let (Some(a), Some(b)) = (a, b) else {
            return false;
        };
        let ia = canonical_idx(&mut pts, &vertex_position(&a), tol);
        let ib = canonical_idx(&mut pts, &vertex_position(&b), tol);
        *counts.entry(ia).or_insert(0) += 1;
        *counts.entry(ib).or_insert(0) += 1;
    }
    counts.values().all(|&c| c % 2 == 0)
}

/// Weld coincident vertices (`ShapeFix_Wire::FixVertexTolerance`).
///
/// Every pair of `Vertex` instances whose registered points lie within `tol`
/// of each other is merged into a single canonical vertex; the shape is rebuilt
/// so all references to the merged vertices point at the canonical instance.
/// An edge whose two (distinct) endpoints collapse to the same canonical is
/// degenerate and is removed from its wire; a closed loop edge whose seam
/// vertex was welded keeps its curve and simply re-points both ends.
///
/// Returns the welded shape and the number of vertex merges performed.
pub fn weld_coincident_vertices(shape: &TopoShape, tol: f64) -> (TopoShape, usize) {
    let builder = TopoBuilder::new();
    let reg = GeometryRegistry::global();
    let verts = vertices_of(shape);

    // First pass: group vertices by position within `tol`. The first vertex of
    // each group becomes the canonical representative.
    let mut vertex_map: HashMap<usize, Vertex> = HashMap::new();
    let mut reps: Vec<Vertex> = Vec::new();
    let mut welded = 0usize;
    for v in &verts {
        let p = reg.vertex_point(&v.0);
        let mut found = None;
        for (ri, r) in reps.iter().enumerate() {
            if reg.vertex_point(&r.0).distance(&p) <= tol {
                found = Some(ri);
                break;
            }
        }
        match found {
            Some(ri) => {
                vertex_map.insert(ptr(&v.0), reps[ri].clone());
                welded += 1;
            }
            None => {
                vertex_map.insert(ptr(&v.0), v.clone());
                reps.push(v.clone());
            }
        }
    }

    if welded == 0 {
        return (shape.clone(), 0);
    }

    // Second pass: rebuild every edge against the canonical endpoints.
    let edges = edges_of(shape);
    let mut edge_map: HashMap<usize, Edge> = HashMap::new();
    let mut edge_removed: HashSet<usize> = HashSet::new();
    for e in &edges {
        let eptr = ptr(&e.0);
        let (a, b) = edge_vertices(e);
        match (a, b) {
            (Some(va), Some(vb)) => {
                let ia_ptr = ptr(&va.0);
                let ib_ptr = ptr(&vb.0);
                let ca = vertex_map.get(&ia_ptr).unwrap_or(&va).clone();
                let cb = vertex_map.get(&ib_ptr).unwrap_or(&vb).clone();
                if is_same(&ca.0, &va.0) && is_same(&cb.0, &vb.0) {
                    // Neither endpoint moved: keep the edge untouched.
                    edge_map.insert(eptr, e.clone());
                } else if is_same(&ca.0, &cb.0) {
                    if ia_ptr == ib_ptr {
                        // A closed loop whose seam vertex merged: keep the
                        // curve, re-point both ends at the canonical vertex.
                        if let Some(ne) = edge_with_canonical_vertices(e, &ca, &cb) {
                            edge_map.insert(eptr, ne);
                        } else {
                            edge_map.insert(eptr, e.clone());
                        }
                    } else {
                        // Two distinct endpoints collapsed to one point: the
                        // edge is degenerate and must disappear.
                        edge_removed.insert(eptr);
                    }
                } else {
                    let pa = reg.vertex_point(&ca.0);
                    let pb = reg.vertex_point(&cb.0);
                    let ne = make_edge_with_vertices(&builder, &pa, &pb, &ca, &cb);
                    edge_map.insert(eptr, ne);
                }
            }
            _ => {
                edge_map.insert(eptr, e.clone());
            }
        }
    }

    let healed = rebuild_shape(
        shape,
        &HealCtx::new(builder, vertex_map, edge_map, edge_removed, HashMap::new()),
    );
    (healed, welded)
}

/// Remove edges shorter than `min_len` (`ShapeFix_Wire::FixSmall`).
///
/// Each short edge is collapsed to a point: its two endpoint vertices are
/// unioned (transitively through chains of short edges) into one canonical
/// vertex, and the short edge is dropped from its wire. A surviving edge whose
/// endpoint was collapsed is rebuilt as a straight segment between the
/// canonical endpoints, so a two-edge corner `[tiny, long]` becomes the single
/// `long` edge spanning the collapsed corner (skipping the tiny edge). A wire
/// whose edges are all removed is dropped from its face.
///
/// Returns the healed shape and the number of removed edges.
pub fn remove_small_edges(shape: &TopoShape, min_len: f64) -> (TopoShape, usize) {
    let builder = TopoBuilder::new();
    let reg = GeometryRegistry::global();
    let edges = edges_of(shape);

    // Identify the short edges.
    let mut small: HashSet<usize> = HashSet::new();
    for e in &edges {
        if edge_length(e, 8) < min_len {
            small.insert(ptr(&e.0));
        }
    }
    if small.is_empty() {
        return (shape.clone(), 0);
    }

    // Union the endpoints of every short edge so they collapse to one vertex.
    let verts = vertices_of(shape);
    let mut idx_of_ptr: HashMap<usize, usize> = HashMap::new();
    for (i, v) in verts.iter().enumerate() {
        idx_of_ptr.insert(ptr(&v.0), i);
    }
    let mut uf = UnionFind::new(verts.len());
    for e in &edges {
        if !small.contains(&ptr(&e.0)) {
            continue;
        }
        let (a, b) = edge_vertices(e);
        if let (Some(va), Some(vb)) = (a, b) {
            let ia = idx_of_ptr[&ptr(&va.0)];
            let ib = idx_of_ptr[&ptr(&vb.0)];
            uf.union(ia, ib);
        }
    }

    // Canonical vertex per union group: the group root (first endpoint seen).
    let mut vertex_map: HashMap<usize, Vertex> = HashMap::new();
    for (i, v) in verts.iter().enumerate() {
        let root = uf.find(i);
        vertex_map.insert(ptr(&v.0), verts[root].clone());
    }

    // Rebuild surviving edges; any edge that touched a collapsed vertex is
    // rebuilt as a segment between the canonical endpoints.
    let mut edge_map: HashMap<usize, Edge> = HashMap::new();
    let mut edge_removed: HashSet<usize> = small.clone();
    for e in &edges {
        let eptr = ptr(&e.0);
        if small.contains(&eptr) {
            continue;
        }
        let (a, b) = edge_vertices(e);
        match (a, b) {
            (Some(va), Some(vb)) => {
                let ca = vertex_map
                    .get(&ptr(&va.0))
                    .cloned()
                    .unwrap_or_else(|| va.clone());
                let cb = vertex_map
                    .get(&ptr(&vb.0))
                    .cloned()
                    .unwrap_or_else(|| vb.clone());
                if is_same(&ca.0, &va.0) && is_same(&cb.0, &vb.0) {
                    edge_map.insert(eptr, e.clone());
                } else if is_same(&ca.0, &cb.0) {
                    // The whole span collapsed — nothing left to keep.
                    edge_removed.insert(eptr);
                } else {
                    let pa = reg.vertex_point(&ca.0);
                    let pb = reg.vertex_point(&cb.0);
                    let ne = make_edge_with_vertices(&builder, &pa, &pb, &ca, &cb);
                    edge_map.insert(eptr, ne);
                }
            }
            _ => {
                edge_map.insert(eptr, e.clone());
            }
        }
    }

    let healed = rebuild_shape(
        shape,
        &HealCtx::new(builder, vertex_map, edge_map, edge_removed, HashMap::new()),
    );
    (healed, small.len())
}

/// Close non-closed wires whose endpoints lie within `tol` of each other
/// (`ShapeFix_Wire::FixClosed`).
///
/// For every wire that is not already geometrically closed but whose last
/// edge's end point is within `tol` of its first edge's start point, a new
/// straight closing edge is appended from the last vertex back to the first
/// vertex. Wires that are already closed (every vertex used exactly twice, or
/// a lone closed loop) are left alone.
///
/// Returns the healed shape and the number of wires that received a closing
/// edge.
pub fn close_open_wires(shape: &TopoShape, tol: f64) -> (TopoShape, usize) {
    let builder = TopoBuilder::new();
    let mut add_closing: HashMap<usize, Edge> = HashMap::new();
    let mut count = 0usize;
    for w in wires_of(shape) {
        if wire_is_closed_loop(&w) {
            continue;
        }
        let edges = edges_of_wire(&w);
        if edges.is_empty() {
            continue;
        }
        let (first_a, _) = edge_vertices(&edges[0]);
        let (_, last_b) = edge_vertices(&edges[edges.len() - 1]);
        let (Some(va), Some(vb)) = (first_a, last_b) else {
            continue;
        };
        let pa = vertex_position(&va);
        let pb = vertex_position(&vb);
        if pa.distance(&pb) <= tol {
            let closer = make_segment_safe(&builder, &pb, &pa);
            add_closing.insert(ptr(&w.0), closer);
            count += 1;
        }
    }
    if count == 0 {
        return (shape.clone(), 0);
    }
    let healed = rebuild_shape(
        shape,
        &HealCtx::new(builder, HashMap::new(), HashMap::new(), HashSet::new(), add_closing),
    );
    (healed, count)
}

/// Apply every fix in sequence (`ShapeFix_Shape`): weld coincident vertices,
/// then remove small edges, then close non-closed wires. The report records
/// how many of each repair fired, and `modified` is set when anything changed.
pub fn heal_shape(shape: &TopoShape, tol: f64, min_edge_len: f64) -> (TopoShape, HealReport) {
    let mut report = HealReport::default();
    let free_before = free_edges(shape, tol).len();

    let (s1, w) = weld_coincident_vertices(shape, tol);
    report.vertices_welded = w;

    let (s2, r) = remove_small_edges(&s1, min_edge_len);
    report.small_edges_removed = r;

    let (s3, c) = close_open_wires(&s2, tol);
    report.wires_closed = c;

    let free_after = free_edges(&s3, tol).len();
    report.free_edges_fixed = free_before.saturating_sub(free_after);

    report.modified = w > 0 || r > 0 || c > 0 || report.free_edges_fixed > 0;
    (s3, report)
}

/// Free boundary edges of a shape (`ShapeAnalysis_FreeBounds`).
///
/// A wire is closed when every vertex position is used by exactly two edge
/// ends (a lone closed loop — one edge whose ends coincide — is also closed).
/// For every open wire this returns the edges that carry a free endpoint: an
/// endpoint used by exactly one edge in that wire. Each returned pair is the
/// `(first, last)` endpoint `Vertex` of a free boundary edge. A closed box's
/// wires are all closed, so the list is empty; a single open wire's edges all
/// show up.
pub fn free_edges(shape: &TopoShape, tol: f64) -> Vec<(Vertex, Vertex)> {
    let mut out = Vec::new();
    for w in wires_of(shape) {
        let edges = edges_of_wire(&w);
        if edges.is_empty() {
            continue;
        }
        // Group endpoint positions; record the canonical index of each end.
        let mut pts: Vec<GpPnt> = Vec::new();
        let mut edge_ends: Vec<(usize, usize)> = Vec::new();
        let mut complete = true;
        for e in &edges {
            let (a, b) = edge_vertices(e);
            let (Some(a), Some(b)) = (a, b) else {
                complete = false;
                break;
            };
            let ia = canonical_idx(&mut pts, &vertex_position(&a), tol);
            let ib = canonical_idx(&mut pts, &vertex_position(&b), tol);
            edge_ends.push((ia, ib));
        }
        if !complete {
            continue;
        }
        // Count how many edge ends use each canonical position.
        let mut usage = vec![0usize; pts.len()];
        for &(ia, ib) in &edge_ends {
            usage[ia] += 1;
            usage[ib] += 1;
        }
        // A closed wire (all ends paired) has no free edges.
        if usage.iter().all(|&c| c == 2) {
            continue;
        }
        for (k, e) in edges.iter().enumerate() {
            let (ia, ib) = edge_ends[k];
            if usage[ia] == 1 || usage[ib] == 1 {
                if let (Some(a), Some(b)) = edge_vertices(e) {
                    out.push((a, b));
                }
            }
        }
    }
    out
}

/// Whether every wire of the shape is geometrically closed after
/// [`heal_shape`] runs with the given tolerance (used for both the welding and
/// the closing tolerance, and as the minimum-edge length).
pub fn wire_is_closed_after_heal(shape: &TopoShape, tol: f64) -> bool {
    let (healed, _) = heal_shape(shape, tol, tol);
    wires_of(&healed)
        .iter()
        .all(|w| wire_geometrically_closed_with_tol(w, tol))
}

/// Relocate a vertex to `to` and update its incident edges' geometry
/// (`ShapeFix_Vertex`).
///
/// The moved vertex is replaced by a fresh `Vertex` instance registered at the
/// target point (the input shape is untouched). Every edge incident to it is
/// rebuilt as a straight segment between the moved endpoint and its other
/// endpoint, so the edge curves pass through the new location. All other
/// vertices and edges are preserved. Errors when `vertex` is not found in
/// `shape`.
pub fn move_vertex(shape: &TopoShape, vertex: &Vertex, to: &GpPnt) -> Result<TopoShape, String> {
    let target_ptr = ptr(&vertex.0);
    let verts = vertices_of(shape);
    if !verts.iter().any(|v| ptr(&v.0) == target_ptr) {
        return Err("move_vertex: vertex not found in shape".into());
    }

    let builder = TopoBuilder::new();
    let reg = GeometryRegistry::global();

    // The moved vertex becomes a new instance at the target point.
    let mut vertex_map: HashMap<usize, Vertex> = HashMap::new();
    for v in &verts {
        let vptr = ptr(&v.0);
        if vptr == target_ptr {
            let nv = builder.make_vertex(*to, reg.vertex_tolerance(v));
            vertex_map.insert(vptr, nv);
        } else {
            vertex_map.insert(vptr, v.clone());
        }
    }

    // Rebuild every incident edge so its curve runs through the new point.
    let edges = edges_of(shape);
    let mut edge_map: HashMap<usize, Edge> = HashMap::new();
    for e in &edges {
        let eptr = ptr(&e.0);
        let (a, b) = edge_vertices(e);
        match (a, b) {
            (Some(va), Some(vb)) => {
                let ia_ptr = ptr(&va.0);
                let ib_ptr = ptr(&vb.0);
                let incident = ia_ptr == target_ptr || ib_ptr == target_ptr;
                if !incident {
                    edge_map.insert(eptr, e.clone());
                    continue;
                }
                let ca = vertex_map.get(&ia_ptr).unwrap().clone();
                let cb = vertex_map.get(&ib_ptr).unwrap().clone();
                let pa = reg.vertex_point(&ca.0);
                let pb = reg.vertex_point(&cb.0);
                let ne = make_edge_with_vertices(&builder, &pa, &pb, &ca, &cb);
                edge_map.insert(eptr, ne);
            }
            _ => {
                edge_map.insert(eptr, e.clone());
            }
        }
    }

    let healed = rebuild_shape(
        shape,
        &HealCtx::new(builder, vertex_map, edge_map, HashSet::new(), HashMap::new()),
    );
    Ok(healed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_builder_api::{make_face_from_polygon, make_wire_from_points};
    use crate::brep_tool::BRepTool;
    use crate::primitives::BRepPrimBox;
    use crate::topo_tools_full::{faces_of, wire_is_closed};

    #[test]
    fn weld_two_boxes_glue_vertices() {
        // Two unit boxes sharing the x = 1 wall.
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 1.0, 1.0));
        let comp = TopoBuilder::new().make_compound_of(&[a.solid.0, b.solid.0]);
        assert_eq!(vertices_of(&comp.0).len(), 16);

        let (healed, welded) = weld_coincident_vertices(&comp.0, 1e-6);
        assert_eq!(welded, 4, "the four shared-wall corner pairs merge");
        assert_eq!(vertices_of(&healed).len(), 12);
    }

    #[test]
    fn weld_within_tolerance() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 0.0);
        let v2 = b.make_vertex(GpPnt::new(1e-4, 0.0, 0.0), 0.0);
        let comp = b.make_compound_of(&[v1.0.clone(), v2.0.clone()]);

        let (merged, n) = weld_coincident_vertices(&comp.0, 1e-3);
        assert_eq!(n, 1, "1e-4 apart merges at tol 1e-3");
        assert_eq!(vertices_of(&merged).len(), 1);

        let (_, n2) = weld_coincident_vertices(&comp.0, 1e-6);
        assert_eq!(n2, 0, "1e-4 apart does not merge at tol 1e-6");
    }

    #[test]
    fn remove_small_edge_collapses() {
        // A tiny first edge followed by a long one: the corner collapses so the
        // wire becomes a single edge from (0,0) to (1,0).
        let w = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.0001, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ])
        .unwrap();
        assert_eq!(edges_of_wire(&w).len(), 2);

        let (healed, n) = remove_small_edges(&w.0, 0.01);
        assert_eq!(n, 1);
        let healed_wire = Wire(healed);
        let he = edges_of_wire(&healed_wire);
        assert_eq!(he.len(), 1, "the tiny edge is collapsed away");
        let (a, z) = edge_vertices(&he[0]);
        assert!(
            vertex_position(&a.unwrap()).distance(&GpPnt::new(0.0, 0.0, 0.0)) < 1e-9,
            "collapsed start stays at the origin"
        );
        assert!(vertex_position(&z.unwrap()).distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9);
    }

    #[test]
    fn close_open_wire_adds_edge() {
        // A 3-edge wire missing the closing edge; endpoints ~1.4e-6 apart.
        let w = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1e-6, 1e-6, 0.0),
        ])
        .unwrap();

        let (healed, n) = close_open_wires(&w.0, 1e-3);
        assert_eq!(n, 1);
        let healed_wire = Wire(healed);
        assert_eq!(edges_of_wire(&healed_wire).len(), 4);
        assert!(wire_is_closed(&healed_wire), "closing edge completes the loop");
    }

    #[test]
    fn heal_shape_full_fixes_all() {
        // A wire with a tiny edge (2e-3: longer than the weld tol so it is the
        // small-edge fix that removes it), plus a near-closed wire whose ends
        // are 1.4e-6 apart (within the weld tol).
        let w1 = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.002, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ])
        .unwrap();
        let w2 = make_wire_from_points(&[
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1.0, 2.0, 0.0),
            GpPnt::new(1e-6, 1.0 + 1e-6, 0.0),
        ])
        .unwrap();
        let comp = TopoBuilder::new().make_compound_of(&[w1.0, w2.0]);

        let (healed, report) = heal_shape(&comp.0, 1e-3, 0.01);
        assert!(report.modified);
        assert!(report.small_edges_removed > 0, "tiny edge removed by FixSmall");
        assert!(report.vertices_welded > 0, "near-closed wire's ends welded");
        assert!(
            edges_of(&healed).iter().all(|e| edge_length(e, 8) >= 0.01 - 1e-9),
            "no edge shorter than the minimum remains"
        );
        // The near-closed wire is now a closed loop inside the healed compound.
        assert!(
            wires_of(&healed).iter().any(|w| wire_is_closed(w)),
            "the near-closed wire closes after healing"
        );
    }

    #[test]
    fn free_edges_shared_vs_free() {
        // A box's wires are all closed → no free edges.
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(free_edges(&boxy.solid.0, 1e-6).is_empty());

        // An open 2-edge wire: both edges carry a free endpoint.
        let w = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
        ])
        .unwrap();
        let free = free_edges(&w.0, 1e-6);
        assert_eq!(free.len(), 2, "both edges of the open wire are free");
    }

    #[test]
    fn wire_closed_after_heal() {
        let w = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1e-6, 1e-6, 0.0),
        ])
        .unwrap();
        assert!(!wire_is_closed(&w), "the wire starts out open");
        assert!(wire_is_closed_after_heal(&w.0, 1e-3));
    }

    #[test]
    fn move_vertex_relocates() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let verts = vertices_of(&b.solid.0);
        let origin = verts
            .iter()
            .find(|v| vertex_position(v).is_equal(&GpPnt::zero()))
            .expect("corner at the origin")
            .clone();

        let healed = move_vertex(&b.solid.0, &origin, &GpPnt::new(0.1, 0.0, 0.0)).unwrap();
        let hv = vertices_of(&healed);
        let moved = hv
            .iter()
            .find(|v| vertex_position(v).distance(&GpPnt::new(0.1, 0.0, 0.0)) < 1e-9)
            .expect("moved vertex present at the new point");

        // An incident edge's curve now starts near the moved point.
        let edges = edges_of(&healed);
        let incident = edges
            .iter()
            .find(|e| {
                let (a, z) = edge_vertices(e);
                let (Some(a), Some(z)) = (a, z) else { return false };
                is_same(&a.0, &moved.0) || is_same(&z.0, &moved.0)
            })
            .expect("incident edge");
        let curve = BRepTool::edge_curve(incident).expect("edge curve");
        let (f0, f1) = BRepTool::edge_parameters(incident);
        let touches_moved = curve.d0(f0).distance(&GpPnt::new(0.1, 0.0, 0.0)) < 1e-9
            || curve.d0(f1).distance(&GpPnt::new(0.1, 0.0, 0.0)) < 1e-9;
        assert!(touches_moved, "an incident edge's curve passes through the moved point");
    }

    #[test]
    fn heal_identity_on_good_shape() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (healed, report) = heal_shape(&b.solid.0, 1e-3, 0.01);
        assert!(!report.modified, "a valid box needs no healing");
        assert_eq!(report.small_edges_removed, 0);
        assert_eq!(report.wires_closed, 0);
        assert_eq!(report.vertices_welded, 0);
        assert_eq!(report.free_edges_fixed, 0);
        assert_eq!(vertices_of(&healed).len(), 8);
        assert_eq!(edges_of(&healed).len(), 12);
    }

    #[test]
    fn weld_merge_deduplicates_vertices() {
        // Two vertex instances at the exact same point merge into one.
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let v2 = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let comp = b.make_compound_of(&[v1.0, v2.0]);

        let (healed, n) = weld_coincident_vertices(&comp.0, 1e-9);
        assert_eq!(n, 1);
        assert_eq!(vertices_of(&healed).len(), 1);
    }

    #[test]
    fn remove_small_isolated_loop() {
        // A tiny standalone triangle loop is removed entirely.
        let face = make_face_from_polygon(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1e-4, 0.0, 0.0),
            GpPnt::new(0.0, 1e-4, 0.0),
        ])
        .unwrap();
        let (healed, n) = remove_small_edges(&face.0, 0.01);
        assert_eq!(n, 3, "all three tiny edges are removed");
        assert!(edges_of(&healed).is_empty(), "the loop is gone");
        assert!(wires_of(&healed).is_empty(), "the empty wire is dropped");
    }

    #[test]
    fn heal_preserves_solid_closure() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (healed, _) = heal_shape(&b.solid.0, 1e-6, 0.01);
        let shells = crate::topo_tools_full::shapes_of(&healed, ShapeType::Shell);
        assert_eq!(shells.len(), 1);
        let shell = Shell(shells[0].clone());
        assert!(crate::shell_check::shell_is_closed(&shell), "shell stays closed");
        let (v, e, f) = (
            vertices_of(&healed).len(),
            edges_of(&healed).len(),
            faces_of(&healed).len(),
        );
        assert_eq!((v, e, f), (8, 12, 6));
    }
}
