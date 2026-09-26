//! Repair options and edge-overlap repair. Split from `bop_builder.rs`.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;


use occt_core::gp::{GpPln, GpPnt, GpVec};
use occt_geom::{GeomPlane};



use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;

use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::shell_check::{shell_is_closed};

use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, shapes_of, vertex_position, vertices_of,
    wires_of_face,
};



use crate::bop_builder_dispatch::{expand_compound};
use crate::bop_builder_repair::{
    analyze_boundary, boolean_repaired, decompose_multi_result, face_connectivity_groups,
    repair_self_intersections, BoundaryAnalysisReport, MultiResult, RepairResult,
};

use crate::bop_builder_planar::*;

// ---------------------------------------------------------------------------
// Repair options
// ---------------------------------------------------------------------------

/// Options that tune [`repair_self_intersections_opts`].
#[derive(Debug, Clone)]
pub struct RepairOptions {
    /// Face-intersection tolerance (points closer than this are coincident).
    pub tolerance: f64,
    /// Maximum number of repair passes (a split can reveal a new crossing).
    pub max_passes: usize,
    /// Drop overlapping coplanar duplicate faces.
    pub remove_coplanar_overlaps: bool,
    /// Emit a warning when the repaired boundary is not closed.
    pub report_open_boundaries: bool,
}

impl Default for RepairOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            max_passes: 1,
            remove_coplanar_overlaps: true,
            report_open_boundaries: true,
        }
    }
}

/// Repair self-intersections with explicit [`RepairOptions`].
///
/// This is the configurable front-end over [`repair_self_intersections`]:
/// `max_passes` runs the repair to a fixpoint (see
/// [`repair_self_intersections_loop`]) and `report_open_boundaries` controls
/// the open-shell warning. `remove_coplanar_overlaps` is advisory — the base
/// repair always removes overlapping coplanar duplicates, so disabling it only
/// affects the documentation of the intent.  (`ponytail:` threaded flag is a
/// no-op; upgrade the pair loop if per-call control ever matters.)
pub fn repair_self_intersections_opts(shape: &TopoShape, opts: &RepairOptions) -> Result<RepairResult, String> {
    let mut result = repair_self_intersections(shape, opts.tolerance)?;
    for _ in 1..opts.max_passes.max(1) {
        if result.fixed_faces == 0 && result.removed_faces == 0 {
            break;
        }
        let next = repair_self_intersections(&result.repaired, opts.tolerance)?;
        result.fixed_faces += next.fixed_faces;
        result.removed_faces += next.removed_faces;
        result.warnings.extend(next.warnings);
        result.repaired = next.repaired;
        if next.fixed_faces == 0 && next.removed_faces == 0 {
            break;
        }
    }
    if !opts.report_open_boundaries {
        result.warnings.retain(|w| !w.contains("not closed"));
    }
    Ok(result)
}

// ---------------------------------------------------------------------------
// Edge-connectivity report (TopOpeBRep)
// ---------------------------------------------------------------------------

/// Per-component edge-connectivity summary of a boundary.
///
/// Each entry corresponds to one edge-connected component of a shape's face
/// set — the same grouping [`decompose_multi_result`] uses to split a boolean
/// result — and reports the face indices and the component's vertex/edge/face
/// counts, from which the Euler characteristic follows.
#[derive(Debug, Clone)]
pub struct ComponentEdgeReport {
    /// Face indices (into `faces_of`) of this component.
    pub face_indices: Vec<usize>,
    /// Distinct vertices in the component (edge endpoints, deduplicated).
    pub vertex_count: usize,
    /// Distinct edges in the component.
    pub edge_count: usize,
    /// Number of faces in the component.
    pub face_count: usize,
}

/// Per-component connectivity reports for a shape's boundary.
///
/// A closed box is one component (6 faces, 12 edges, 8 vertices → Euler 2). A
/// shell built from two crossing faces splits into one component per face (each
/// is a 1-face open component). Components are returned sorted by smallest face
/// index.
pub fn component_edge_reports(shape: &TopoShape, tol: f64) -> Vec<ComponentEdgeReport> {
    let faces = faces_of(shape);
    let groups = face_connectivity_groups(&faces);
    let mut reports: Vec<ComponentEdgeReport> = Vec::with_capacity(groups.len());
    for idx in groups {
        let mut edges: HashSet<usize> = HashSet::new();
        let mut verts: HashSet<usize> = HashSet::new();
        for &i in &idx {
            for e in edges_of(&faces[i].0) {
                edges.insert(Arc::as_ptr(&e.0.tshape) as usize);
                let (a, b) = crate::topo_tools_full::edge_vertices(&e);
                if let Some(a) = a {
                    verts.insert(Arc::as_ptr(&a.0.tshape) as usize);
                }
                if let Some(b) = b {
                    verts.insert(Arc::as_ptr(&b.0.tshape) as usize);
                }
            }
        }
        reports.push(ComponentEdgeReport {
            face_count: idx.len(),
            face_indices: idx,
            vertex_count: verts.len(),
            edge_count: edges.len(),
        });
    }
    let _ = tol;
    reports
}

/// The Euler characteristic V − E + F of a boundary component.
///
/// A closed manifold component has χ = 2 (a sphere-like boundary); an open
/// component has χ = 1 or less. This is the `TopExp`-style sanity check for a
/// decomposed result piece.
pub fn component_euler_characteristic(r: &ComponentEdgeReport) -> i32 {
    (r.vertex_count as i32) - (r.edge_count as i32) + (r.face_count as i32)
}

// ---------------------------------------------------------------------------
// Edge-overlap repair (BOPAlgo / TopOpeBRep depth)
// ---------------------------------------------------------------------------
//
// Ports of the edge-level repair and classification passes of `BOPAlgo_Builder`
// and `TopOpeBRep_BuildTool`: collinear overlapping-edge repair, face splitting
// along surface intersections, boolean-result edge classification and tolerance
// healing. All functions follow the crate's `Result<_, String>` convention and
// rebuild shapes through the weld/edge-map machinery the exact boolean uses, so
// coincident vertices and edges stay shared across the rebuilt boundary.

/// Is the edge's curve geometrically a straight line?
///
/// Samples 8 points along the curve and checks that they are collinear with the
/// first–last chord. A `GeomLine` (even trimmed) passes; a circle/arc fails.
fn edge_is_line_like(e: &Edge) -> bool {
    let Some(c) = BRepTool::edge_curve(e) else { return false };
    let (a, b) = BRepTool::edge_parameters(e);
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return false;
    }
    let n = 8;
    let p0 = c.d0(a);
    let pl = c.d0(b);
    let size = p0.distance(&pl);
    if size <= 1e-30 {
        return false;
    }
    let chord = GpVec::from_pnts(&p0, &pl);
    let tol = 1e-6 * size;
    for i in 1..n {
        let p = c.d0(a + (b - a) * i as f64 / n as f64);
        if GpVec::from_pnts(&p0, &p).cross_magnitude(&chord) > tol * size {
            return false;
        }
    }
    true
}

/// Report of [`repair_edge_overlaps`].
#[derive(Debug, Clone)]
pub struct RepairReport {
    /// The repaired shape (unchanged when nothing needed repairing).
    pub repaired: TopoShape,
    /// Number of edges that were split to resolve a partial overlap.
    pub repaired_edges: usize,
    /// Number of edges removed (zero-length edges and coincident duplicates).
    pub removed_edges: usize,
    /// Number of vertex pairs merged because they lay within `tol`.
    pub welded_vertices: usize,
    /// Non-fatal diagnostics (open rebuilt boundary, dropped faces, …).
    pub warnings: Vec<String>,
}

/// Detect and repair collinear overlapping edges of a boolean result/compound.
///
/// Looks at every pair of line-like edges of the shape and repairs three kinds
/// of overlap:
///
/// * **zero-length edges** — an edge whose two endpoints coincide (within
///   `tol`) is removed;
/// * **coincident duplicates** — two collinear edges spanning the same
///   interval are merged into one (the duplicate is counted as removed);
/// * **partial overlaps** — two collinear edges that overlap over a sub-segment
///   are split at the overlap boundaries so the shared sub-segment becomes a
///   single edge (each split edge counts as repaired).
///
/// The boundary is rebuilt through the same `Weld`/`EdgeMap` machinery the exact
/// boolean uses, so every vertex closer than `tol` is welded to one instance and
/// coincident sub-segments resolve to the same `Edge`. Non-line edges (arcs,
/// full circles) and non-planar faces are kept untouched. A `Compound` is
/// repaired component-wise and reassembled.
pub fn repair_edge_overlaps(shape: &TopoShape, tol: f64) -> Result<RepairReport, String> {
    let tol = tol.max(1e-9);
    if shape.is_compound() {
        let children = expand_compound(shape);
        let mut repaired_edges = 0usize;
        let mut removed_edges = 0usize;
        let mut welded_vertices = 0usize;
        let mut warnings: Vec<String> = Vec::new();
        let mut out: Vec<TopoShape> = Vec::with_capacity(children.len());
        for c in &children {
            let r = repair_edge_overlaps(c, tol)?;
            repaired_edges += r.repaired_edges;
            removed_edges += r.removed_edges;
            welded_vertices += r.welded_vertices;
            warnings.extend(r.warnings);
            out.push(r.repaired);
        }
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&out);
        return Ok(RepairReport { repaired: comp.0, repaired_edges, removed_edges, welded_vertices, warnings });
    }

    let faces = faces_of(shape);
    let edges = edges_of(shape);
    let n = edges.len();

    // Per-edge geometry: endpoints (from the curve), line-likeness and length.
    let mut pts: Vec<Option<(GpPnt, GpPnt)>> = Vec::with_capacity(n);
    let mut is_line: Vec<bool> = Vec::with_capacity(n);
    let mut lens: Vec<f64> = Vec::with_capacity(n);
    for e in &edges {
        let ep = BRepTool::edge_vertices(e);
        let len = ep.as_ref().map(|(a, b)| a.distance(b)).unwrap_or(0.0);
        let line = ep.as_ref().map(|_| edge_is_line_like(e)).unwrap_or(false);
        pts.push(ep);
        is_line.push(line);
        lens.push(len);
    }

    // Removed: zero-length edges.
    let mut removed: HashSet<usize> = HashSet::new();
    for i in 0..n {
        if lens[i] < tol {
            removed.insert(i);
        }
    }

    // Pass 1: coincident-duplicate detection (same span on the same line).
    let mut merge_with: Vec<Option<usize>> = vec![None; n];
    let mut merged_count = 0usize;
    for i in 0..n {
        if removed.contains(&i) || !is_line[i] {
            continue;
        }
        let Some((a1, a2)) = pts[i] else { continue };
        let li = lens[i];
        let di = GpVec::from_pnts(&a1, &a2).normalized();
        for j in (i + 1)..n {
            if removed.contains(&j) || !is_line[j] || merge_with[j].is_some() {
                continue;
            }
            let Some((b1, b2)) = pts[j] else { continue };
            let db = GpVec::from_pnts(&b1, &b2).normalized();
            if di.cross_magnitude(&db) > tol {
                continue;
            }
            if GpVec::from_pnts(&a1, &b1).cross_magnitude(&di) > tol {
                continue;
            }
            let t = |p: &GpPnt| GpVec::from_pnts(&a1, p).dot(&di);
            let (lo, hi) = (t(&b1).min(t(&b2)), t(&b1).max(t(&b2)));
            if (lo - 0.0).abs() <= tol && (hi - li).abs() <= tol {
                merge_with[j] = Some(i);
                merged_count += 1;
            }
        }
    }

    // Pass 2: collect split parameters (in edge-local arc-length) for every
    // edge whose span is cut by a collinear overlapping edge.
    let mut split_ts: Vec<Vec<f64>> = vec![Vec::new(); n];
    for i in 0..n {
        if removed.contains(&i) || !is_line[i] {
            continue;
        }
        let Some((a1, a2)) = pts[i] else { continue };
        let li = lens[i];
        let di = GpVec::from_pnts(&a1, &a2).normalized();
        let t = |p: &GpPnt| GpVec::from_pnts(&a1, p).dot(&di);
        let mut ts = vec![0.0, li];
        for j in 0..n {
            if i == j || removed.contains(&j) || !is_line[j] {
                continue;
            }
            let Some((b1, b2)) = pts[j] else { continue };
            let db = GpVec::from_pnts(&b1, &b2).normalized();
            if di.cross_magnitude(&db) > tol {
                continue;
            }
            if GpVec::from_pnts(&a1, &b1).cross_magnitude(&di) > tol {
                continue;
            }
            let (tj1, tj2) = (t(&b1), t(&b2));
            let (lo, hi) = (tj1.min(tj2), tj1.max(tj2));
            if hi.min(li) - lo.max(0.0) <= tol {
                continue; // disjoint or merely touching at an endpoint
            }
            ts.push(lo.max(0.0).min(li));
            ts.push(hi.max(0.0).min(li));
        }
        ts.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        ts.dedup_by(|x, y| (*x - *y).abs() <= tol.max(1e-12));
        split_ts[i] = ts;
    }
    let repaired_count = split_ts.iter().filter(|ts| ts.len() > 2).count();

    // Vertex merges: distinct vertex instances that collapse into one position.
    let mut welded_vertices = 0usize;
    let mut positions: Vec<GpPnt> = Vec::new();
    for v in vertices_of(shape) {
        let p = vertex_position(&v);
        if positions.iter().any(|q| q.distance(&p) <= tol) {
            welded_vertices += 1;
        } else {
            positions.push(p);
        }
    }

    // Edges whose identity must change: removed, split, or merged with a twin.
    let mut affected: HashSet<usize> = HashSet::new();
    for i in 0..n {
        if removed.contains(&i) || split_ts[i].len() > 2 {
            affected.insert(i);
        }
    }
    for (j, k) in merge_with.iter().enumerate() {
        if let Some(i) = k {
            affected.insert(*i);
            affected.insert(j);
        }
    }
    if affected.is_empty() {
        return Ok(RepairReport {
            repaired: shape.clone(),
            repaired_edges: 0,
            removed_edges: removed.len() + merged_count,
            welded_vertices,
            warnings: vec![],
        });
    }

    // Rebuild the boundary: only faces carrying an affected edge are rebuilt.
    let bld = TopoBuilder::new();
    let mut edge_idx: HashMap<usize, usize> = HashMap::new();
    for (i, e) in edges.iter().enumerate() {
        edge_idx.insert(Arc::as_ptr(&e.0.tshape) as usize, i);
    }
    let mut weld = Weld::new(tol.max(1e-7));
    let mut edge_map = EdgeMap::default();
    let mut result_faces: Vec<Face> = Vec::new();
    let mut warnings: Vec<String> = Vec::new();

    for f in &faces {
        let needs = edges_of(&f.0).iter().any(|e| {
            edge_idx
                .get(&(Arc::as_ptr(&e.0.tshape) as usize))
                .map_or(false, |idx| affected.contains(idx))
        });
        if !needs {
            result_faces.push(f.clone());
            continue;
        }
        let Some(pln) = face_plane_local(f) else {
            result_faces.push(f.clone());
            continue;
        };
        match rebuild_face_repaired(&bld, f, &pln, &edge_idx, &is_line, &pts, &split_ts, &removed, &mut weld, &mut edge_map)
        {
            Some(nf) => result_faces.push(nf),
            None => warnings.push("repair_edge_overlaps: a face collapsed to nothing and was dropped".into()),
        }
    }

    if result_faces.is_empty() {
        let empty = bld.make_compound_of(&[]);
        return Ok(RepairReport {
            repaired: empty.0,
            repaired_edges: repaired_count,
            removed_edges: removed.len() + merged_count,
            welded_vertices,
            warnings,
        });
    }

    let shell = bld.make_shell(&result_faces);
    let closed = shell_is_closed(&shell);
    if !closed {
        warnings.push("repair_edge_overlaps: rebuilt boundary is not closed".into());
    }
    let repaired: TopoShape = if closed { bld.make_solid(&[shell]).0 } else { shell.0 };
    Ok(RepairReport {
        repaired,
        repaired_edges: repaired_count,
        removed_edges: removed.len() + merged_count,
        welded_vertices,
        warnings,
    })
}

/// Rebuild a single face's wires with repaired edges.
///
/// Removed edges are dropped; split edges become their ordered sub-segments
/// (reversed when the wire traverses the edge backwards); unsplit line edges are
/// re-registered through the shared weld/edge-map so coincident edges across
/// faces resolve to one `Edge`; non-line edges are kept untouched. Returns
/// `None` when every wire collapsed.
#[allow(clippy::too_many_arguments)]
fn rebuild_face_repaired(
    bld: &TopoBuilder,
    f: &Face,
    pln: &GpPln,
    edge_idx: &HashMap<usize, usize>,
    is_line: &[bool],
    pts: &[Option<(GpPnt, GpPnt)>],
    split_ts: &[Vec<f64>],
    removed: &HashSet<usize>,
    weld: &mut Weld,
    edge_map: &mut EdgeMap,
) -> Option<Face> {
    let mut new_wires: Vec<Wire> = Vec::new();
    for w in wires_of_face(f) {
        let mut new_edges: Vec<Edge> = Vec::new();
        for we in edges_of_wire(&w) {
            let eptr = Arc::as_ptr(&we.0.tshape) as usize;
            let Some(&idx) = edge_idx.get(&eptr) else {
                new_edges.push(we.clone());
                continue;
            };
            if removed.contains(&idx) {
                continue;
            }
            let Some((p1, p2)) = pts[idx] else {
                new_edges.push(we.clone());
                continue;
            };
            let forward = !we.0.orientation().is_reversed();
            if !is_line[idx] {
                new_edges.push(we.clone());
                continue;
            }
            let ts = &split_ts[idx];
            if ts.len() >= 2 {
                let di = GpVec::from_pnts(&p1, &p2).normalized();
                let point_at = |t: f64| p1.translated_vec(&di.multiplied_scalar(t));
                let segs: Vec<(GpPnt, GpPnt)> = if forward {
                    ts.windows(2).map(|w| (point_at(w[0]), point_at(w[1]))).collect()
                } else {
                    ts.windows(2).rev().map(|w| (point_at(w[1]), point_at(w[0]))).collect()
                };
                for (a, b) in segs {
                    let ia = weld.weld(&a);
                    let ib = weld.weld(&b);
                    if ia != ib {
                        new_edges.push(edge_map.edge(bld, ia, ib, &weld.points));
                    }
                }
            } else {
                let (a, b) = if forward { (p1, p2) } else { (p2, p1) };
                let ia = weld.weld(&a);
                let ib = weld.weld(&b);
                if ia != ib {
                    new_edges.push(edge_map.edge(bld, ia, ib, &weld.points));
                }
            }
        }
        if !new_edges.is_empty() {
            new_wires.push(bld.make_wire(&new_edges));
        }
    }
    if new_wires.is_empty() {
        return None;
    }
    let mut nf = bld.make_face(Arc::new(GeomPlane::new(pln.clone())), &new_wires);
    nf.0.set_orientation(f.0.orientation());
    Some(nf)
}
