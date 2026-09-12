//! Self-intersection repair wrappers. Split from `bop_builder.rs`.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::geom::polygon_ops::{point_in_polygon2d, polygon_area2d};
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::{GeomPlane, Surface};

use crate::abs::ShapeType;
use crate::brep_extrema::is_inside;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::inttools::{edge_edge_intersections, edge_face_intersections};
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::shell_check::{shell_invariants, shell_is_closed};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, shapes_of, vertex_position, vertices_of,
    wires_of_face,
};

use crate::bop_builder::boolean;
use crate::bop_builder_core::{empty_result, single_shape_result, BoolOp, BooleanResult};
use crate::bop_builder_dispatch::{
    boolean_with_check, decompose_compound, detect_self_intersections, expand_compound,
    faces_polygon_overlap,
};
use crate::bop_builder_planar::*;

// ---------------------------------------------------------------------------
// Self-intersection repair
// ---------------------------------------------------------------------------
//
// A `BOPAlgo_ArgumentAnalyzer`-lite repair pass. `BOPAlgo_ArgumentAnalyzer`
// *detects* self-intersections (non-adjacent faces whose surfaces cross or
// overlap); the "repair" done here is the resolution step: split every
// crossing face along the intersection polyline and weld the split edges so
// the crossing faces become adjacent, and drop degenerate slivers and
// overlapping coplanar duplicates. This is the same split-and-weld machinery
// the exact boolean uses to produce manifold boundaries.

/// Result of [`repair_self_intersections`].
pub struct RepairResult {
    /// The repaired shape (unchanged when no self-intersection was found).
    pub repaired: TopoShape,
    /// Number of faces that were split to resolve a transversal crossing.
    pub fixed_faces: usize,
    /// Number of faces removed (overlapping coplanar duplicates, or every
    /// sub-face degenerating to zero area).
    pub removed_faces: usize,
    /// Non-fatal diagnostics (unrepairable non-planar crossings, open result).
    pub warnings: Vec<String>,
}

/// Repair a self-intersecting boundary (a `BOPAlgo_ArgumentAnalyzer`-lite port).
///
/// * Runs [`detect_self_intersections`]; a clean shape is returned unchanged
///   with `fixed_faces == 0`.
/// * For every flagged transversal crossing pair, the intersection polyline is
///   computed with [`crate::intpatch::surface_surface_intersection`], clipped
///   to both face polygons, and both faces are split along it through the same
///   weld/edge-map the boolean uses. The split edges coincide, so the crossing
///   faces become adjacent and the defect is resolved.
/// * Overlapping *coplanar* non-adjacent faces are resolved by removing the
///   duplicate (the overlap sliver); degenerate (zero-area) sub-faces are
///   dropped by the splitter.
/// * The surviving faces are rebuilt into a shell — or a solid when the shell
///   is closed — and an open result is reported as a warning.
/// * A `Compound` is repaired component-wise and reassembled.
pub fn repair_self_intersections(shape: &TopoShape, tol: f64) -> Result<RepairResult, String> {
    let tol = tol.max(1e-9);
    if shape.is_compound() {
        let mut fixed = 0usize;
        let mut removed = 0usize;
        let mut warnings: Vec<String> = Vec::new();
        let children = expand_compound(shape);
        let mut out: Vec<TopoShape> = Vec::with_capacity(children.len());
        for c in &children {
            let r = repair_self_intersections(c, tol)?;
            fixed += r.fixed_faces;
            removed += r.removed_faces;
            warnings.extend(r.warnings);
            out.push(r.repaired);
        }
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&out);
        return Ok(RepairResult { repaired: comp.0, fixed_faces: fixed, removed_faces: removed, warnings });
    }

    let report = detect_self_intersections(shape, tol);
    if !report.found {
        return Ok(RepairResult { repaired: shape.clone(), fixed_faces: 0, removed_faces: 0, warnings: vec![] });
    }

    let faces = faces_of(shape);
    let planes: Vec<Option<GpPln>> = faces.iter().map(face_plane_local).collect();
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();

    let mut segs: Vec<Vec<(GpPnt, GpPnt)>> = vec![Vec::new(); faces.len()];
    let mut drop: Vec<bool> = vec![false; faces.len()];
    let mut removed_count = 0usize;
    let mut warnings: Vec<String> = Vec::new();

    // Collect the split segments per face and the faces to drop.
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if face_edges[i].iter().any(|e| face_edges[j].contains(e)) {
                continue; // adjacent faces legitimately meet along an edge
            }
            let (Some(sa), Some(sb)) = (
                GeometryRegistry::global().face_surface(&faces[i].0),
                GeometryRegistry::global().face_surface(&faces[j].0),
            ) else {
                continue;
            };
            match crate::intpatch::surface_surface_intersection(&*sa, &*sb, tol) {
                crate::intpatch::SurfaceIntersection::Curves(curves) => {
                    if curves.is_empty() {
                        continue;
                    }
                    if planes[i].is_none() || planes[j].is_none() {
                        warnings.push(format!(
                            "self-intersection between faces {i} and {j} involves a non-planar face; left unrepaired"
                        ));
                        continue;
                    }
                    // The transversal crossing: the clipped intersection
                    // polyline is the cutting line on both faces.
                    let seg = face_face_segments_local(&faces[i], &faces[j], tol);
                    if seg.is_empty() {
                        continue;
                    }
                    for s in &seg {
                        segs[i].push(*s);
                        segs[j].push(*s);
                    }
                }
                crate::intpatch::SurfaceIntersection::Coincident => {
                    if faces_polygon_overlap(&faces[i], &faces[j], tol) {
                        if !drop[j] {
                            drop[j] = true;
                            removed_count += 1;
                            warnings.push(format!("removed overlapping coplanar face {j}"));
                        }
                    }
                }
                crate::intpatch::SurfaceIntersection::None => {}
            }
        }
    }

    // Count faces that actually split (an interior cutting segment). A segment
    // running exactly along a face boundary leaves the face whole, so it does
    // not count as "fixed".
    let mut fixed_count = 0usize;
    for i in 0..faces.len() {
        if drop[i] || segs[i].is_empty() {
            continue;
        }
        let Some(pln) = planes[i].clone() else { continue };
        let Some(poly) = face_polygon_local(&faces[i], &pln) else { continue };
        let segs2d: Vec<(GpPnt2d, GpPnt2d)> = segs[i]
            .iter()
            .filter_map(|(a, b)| {
                let a2 = project_point_to_plane(&pln, a);
                let b2 = project_point_to_plane(&pln, b);
                if a2.distance(&b2) < 1e-12 {
                    None
                } else {
                    Some((a2, b2))
                }
            })
            .collect();
        if trace_planar_regions(&poly, &segs2d).len() >= 2 {
            fixed_count += 1;
        }
    }

    // Rebuild the boundary: split every kept planar face along its cutting
    // segments, keeping non-planar faces whole.
    let bld = TopoBuilder::new();
    let mut kept_faces: Vec<Face> = Vec::new();
    let mut kept_planes: Vec<Option<GpPln>> = Vec::new();
    let mut kept_segs: Vec<Vec<(GpPnt, GpPnt)>> = Vec::new();
    let mut non_planar: Vec<Face> = Vec::new();
    for (i, f) in faces.iter().enumerate() {
        if drop[i] {
            continue;
        }
        if planes[i].is_some() {
            kept_faces.push(f.clone());
            kept_planes.push(planes[i].clone());
            kept_segs.push(segs[i].clone());
        } else {
            non_planar.push(f.clone());
        }
    }

    let mut weld = Weld::new(tol.max(1e-7));
    let mut edge_map = EdgeMap::default();
    let subs = split_faces(&bld, &kept_faces, &kept_planes, &kept_segs, &mut weld, &mut edge_map);
    let mut result_faces: Vec<Face> = non_planar;
    result_faces.extend(subs.into_iter().map(|sf| sf.face));

    if result_faces.is_empty() {
        warnings.push("repair removed every face; the result is empty".into());
        let empty = bld.make_compound_of(&[]);
        return Ok(RepairResult {
            repaired: empty.0,
            fixed_faces: fixed_count,
            removed_faces: removed_count,
            warnings,
        });
    }

    let shell = bld.make_shell(&result_faces);
    let closed = shell_is_closed(&shell);
    if !closed {
        warnings.push("repaired boundary is not closed (open shell)".into());
    }
    let repaired: TopoShape = if closed { bld.make_solid(&[shell]).0 } else { shell.0 };

    Ok(RepairResult { repaired, fixed_faces: fixed_count, removed_faces: removed_count, warnings })
}

/// Run a boolean and repair a self-intersecting result.
///
/// Runs [`crate::bop_curved::curved_boolean_full`], then
/// [`repair_self_intersections`] on the result. When the repair changed the
/// shape, the returned [`BooleanResult`] carries the repaired shape and a
/// warning mentioning the number of faces fixed/removed; a clean result is
/// returned unchanged.
pub fn boolean_repaired(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let r = crate::bop_curved::curved_boolean_full(a, b, op, tol)?;
    Ok(repair_boolean_result(r, tol))
}

/// Apply [`repair_self_intersections`] to an already-computed boolean result.
///
/// Shared by [`boolean_repaired`] and [`boolean_repaired_with_check`]. A clean
/// result is returned unchanged; otherwise the repaired shape replaces the
/// result shape (recomputing `solid`/`shells`/`faces`) and a warning naming the
/// number of fixed/removed faces is appended.
pub(crate) fn repair_boolean_result(r: BooleanResult, tol: f64) -> BooleanResult {
    let rep = match repair_self_intersections(&r.shape, tol) {
        Ok(rep) => rep,
        Err(e) => {
            let mut rr = r;
            rr.warnings.push(format!("self-intersection repair failed: {e}"));
            return rr;
        }
    };
    if rep.fixed_faces == 0 && rep.removed_faces == 0 {
        return r;
    }
    let mut rr = single_shape_result(&rep.repaired);
    rr.warnings.extend(r.warnings);
    rr.warnings.extend(rep.warnings);
    rr.warnings.push(format!(
        "boolean result had self-intersections; repair fixed {} face(s), removed {}",
        rep.fixed_faces, rep.removed_faces
    ));
    rr
}

// ---------------------------------------------------------------------------
// Multi-result decomposition
// ---------------------------------------------------------------------------

/// A boolean result decomposed into its disconnected pieces.
///
/// Mirrors the assembly step of `BOPAlgo_Builder` / `TopOpeBRep`: a boolean
/// result is a compound (or a multi-shell solid) holding several disconnected
/// results; this struct splits them apart and buckets them by shape type.
pub struct MultiResult {
    /// Every top-level result shape (compounds flattened, multi-shell solids
    /// split into one shape per connected boundary component).
    pub shapes: Vec<TopoShape>,
    /// The shapes of [`MultiResult::shapes`] that are solids.
    pub solids: Vec<TopoShape>,
    /// The shapes of [`MultiResult::shapes`] that are (open) shells.
    pub shells: Vec<TopoShape>,
    /// The shapes of [`MultiResult::shapes`] that are still compounds.
    pub compounds: Vec<TopoShape>,
}

/// Split a non-compound boundary into its edge-connected components, preserving
/// solidity: a closed component is returned as a one-shell solid, an open one
/// as a bare shell.
pub(crate) fn split_connected_boundaries(shape: &TopoShape, tol: f64) -> Vec<TopoShape> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return vec![shape.clone()];
    }
    if faces.len() == 1 {
        let bld = TopoBuilder::new();
        return vec![bld.make_shell(&faces).0];
    }
    let groups = face_connectivity_groups(&faces);
    let bld = TopoBuilder::new();
    let mut comps: Vec<TopoShape> = Vec::new();
    for idx in groups {
        let fs: Vec<Face> = idx.iter().map(|&i| faces[i].clone()).collect();
        let shell = bld.make_shell(&fs);
        if shell_is_closed(&shell) {
            comps.push(bld.make_solid(&[shell]).0);
        } else {
            comps.push(shell.0);
        }
    }
    comps
}

/// Group face indices by edge connectivity.
///
/// Two faces belong to the same group when they are connected through a chain
/// of shared boundary edges (identical `TShape` edge references). The groups
/// are returned sorted by their smallest face index, so the ordering is stable
/// for a given shape.
pub(crate) fn face_connectivity_groups(faces: &[Face]) -> Vec<Vec<usize>> {
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if face_edges[i].iter().any(|e| face_edges[j].contains(e)) {
                uf_unite(&mut parent, i, j);
            }
        }
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..faces.len() {
        groups.entry(uf_find(&mut parent, i)).or_default().push(i);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort_by_key(|g| g[0]);
    out
}

/// Union-find find with path compression.
fn uf_find(parent: &mut [usize], x: usize) -> usize {
    let mut r = x;
    while parent[r] != r {
        parent[r] = parent[parent[r]];
        r = parent[r];
    }
    r
}

/// Union-find union (root of `a` keeps the group).
fn uf_unite(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (uf_find(parent, a), uf_find(parent, b));
    if ra != rb {
        parent[ra] = rb;
    }
}

/// Decompose a shape into its edge-connected boundary components.
///
/// A public, shape-level form of the per-result decomposition used by
/// [`decompose_multi_result`]: compounds are flattened and every solid/shell is
/// split into one shape per connected boundary component (closed components
/// come back as one-shell solids, open ones as shells).
pub fn shape_boundary_components(shape: &TopoShape, tol: f64) -> Vec<TopoShape> {
    if shape.is_compound() {
        let mut out: Vec<TopoShape> = Vec::new();
        for s in expand_compound(shape) {
            out.extend(split_connected_boundaries(&s, tol));
        }
        out
    } else {
        split_connected_boundaries(shape, tol)
    }
}

/// Decompose a boolean result into its disconnected pieces.
///
/// * the top-level shape list comes from [`decompose_compound`] (compounds are
///   flattened);
/// * every remaining solid/shell is further split into its edge-connected
///   boundary components, so a multi-shell Cut result yields one shape per
///   piece;
/// * the pieces are bucketed into solids / shells / compounds.
pub fn decompose_multi_result(r: &BooleanResult) -> MultiResult {
    let out = shape_boundary_components(&r.shape, 1e-6);
    let solids = out.iter().filter(|s| s.is_solid()).cloned().collect();
    let shells = out.iter().filter(|s| s.is_shell()).cloned().collect();
    let compounds = out.iter().filter(|s| s.is_compound()).cloned().collect();
    MultiResult { shapes: out, solids, shells, compounds }
}

/// Run a boolean and decompose the result into its disconnected pieces.
///
/// A boolean that produces several disjoint results — a Fuse of disjoint
/// inputs (a compound), a Cut that leaves two or more pieces (a multi-shell
/// solid) — yields a [`MultiResult`] with one shape per piece.
pub fn boolean_split_result(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<MultiResult, String> {
    let r = crate::bop_curved::curved_boolean_full(a, b, op, tol)?;
    Ok(decompose_multi_result(&r))
}

/// True when the result decomposes into more than one non-empty shape.
///
/// A disjoint Fuse (a compound) or a Cut into several pieces is "multiple"; a
/// single overlapping Fuse/Common result is not.
pub fn result_is_multiple(r: &BooleanResult) -> bool {
    decompose_multi_result(r).shapes.len() > 1
}

/// Remove topological garbage from a shape.
///
/// Weld near-coincident vertices (merging duplicate vertex instances and
/// dropping edges whose two endpoints collapse to one point), then remove edges
/// shorter than `tol`. This is the `ShapeFix_Shape`-style hygiene pass that
/// keeps rebuilt boundaries minimal. The input shape is left untouched.
pub fn topologically_clean(shape: &TopoShape, tol: f64) -> TopoShape {
    let tol = tol.max(1e-9);
    let (weld, _) = crate::shhealing::weld_coincident_vertices(shape, tol);
    let (clean, _) = crate::shhealing::remove_small_edges(&weld, tol);
    clean
}

// ---------------------------------------------------------------------------
// Self-intersection analysis (BOPAlgo_ArgumentAnalyzer)
// ---------------------------------------------------------------------------
//
// `BOPAlgo_ArgumentAnalyzer` does not just report "self-intersection found" —
// it enumerates the offending face pairs, classifies the defect (transversal
// crossing vs. coplanar overlap), and exposes the intersection data. The items
// below are that detailed analysis, on top of the boolean summary
// [`detect_self_intersections`].

/// Classification of a self-intersection defect.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelfIntersectionKind {
    /// Two non-adjacent faces cross transversally: their surfaces intersect
    /// along a curve that passes through both face patches.
    Crossing,
    /// Two non-adjacent coplanar faces overlap in area.
    CoplanarOverlap,
}

impl SelfIntersectionKind {
    /// Human-readable label of the defect kind.
    pub fn label(&self) -> &'static str {
        match self {
            SelfIntersectionKind::Crossing => "transversal crossing",
            SelfIntersectionKind::CoplanarOverlap => "coplanar overlap",
        }
    }
}

/// One self-intersection defect found by [`analyze_self_intersections`].
///
/// Mirrors one entry of `BOPAlgo_ArgumentAnalyzer`'s self-interference list:
/// the two non-adjacent faces involved, the defect kind and the geometric
/// evidence (intersection sample points, and the segment clipped to both face
/// patches when it is computable).
#[derive(Debug, Clone)]
pub struct SelfIntersectionIssue {
    /// Index (into `faces_of`) of the first face of the pair.
    pub face_a: usize,
    /// Index (into `faces_of`) of the second face of the pair.
    pub face_b: usize,
    /// What kind of defect the pair exhibits.
    pub kind: SelfIntersectionKind,
    /// Sample points on the surface intersection (empty for a pure coplanar
    /// overlap, which has no curve).
    pub points: Vec<GpPnt>,
    /// The intersection segment clipped to both face polygons, when the two
    /// faces are planar and the line actually crosses both patches.
    pub segment: Option<(GpPnt, GpPnt)>,
}

/// Enumerate every self-intersection defect of a shape.
///
/// This is the per-pair detail behind [`detect_self_intersections`]: it walks
/// the same non-adjacent face pairs, but reports each defect with its kind and
/// geometry instead of collapsing them into a single boolean flag. A valid
/// closed box returns an empty list; a shell built from two crossing faces
/// returns one [`SelfIntersectionIssue`] of kind `Crossing`.
pub fn analyze_self_intersections(shape: &TopoShape, tol: f64) -> Vec<SelfIntersectionIssue> {
    let tol = tol.max(1e-9);
    let faces = faces_of(shape);
    if faces.len() < 2 {
        return Vec::new();
    }
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();
    let mut issues: Vec<SelfIntersectionIssue> = Vec::new();
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if face_edges[i].iter().any(|e| face_edges[j].contains(e)) {
                continue; // adjacent faces legitimately meet along an edge
            }
            let (Some(sa), Some(sb)) = (
                GeometryRegistry::global().face_surface(&faces[i].0),
                GeometryRegistry::global().face_surface(&faces[j].0),
            ) else {
                continue;
            };
            match crate::intpatch::surface_surface_intersection(&*sa, &*sb, tol) {
                crate::intpatch::SurfaceIntersection::Curves(curves) => {
                    let mut pts: Vec<GpPnt> = Vec::new();
                    for c in &curves {
                        pts.extend(c.points.iter().cloned());
                    }
                    if pts.is_empty() {
                        continue;
                    }
                    let clipped = face_face_segments_local(&faces[i], &faces[j], tol);
                    let segment = if clipped.is_empty() {
                        None
                    } else {
                        Some((clipped[0].0, clipped[0].1))
                    };
                    issues.push(SelfIntersectionIssue {
                        face_a: i,
                        face_b: j,
                        kind: SelfIntersectionKind::Crossing,
                        points: pts,
                        segment,
                    });
                }
                crate::intpatch::SurfaceIntersection::Coincident => {
                    if faces_polygon_overlap(&faces[i], &faces[j], tol) {
                        issues.push(SelfIntersectionIssue {
                            face_a: i,
                            face_b: j,
                            kind: SelfIntersectionKind::CoplanarOverlap,
                            points: Vec::new(),
                            segment: None,
                        });
                    }
                }
                crate::intpatch::SurfaceIntersection::None => {}
            }
        }
    }
    issues
}

/// One-line summary of a list of [`SelfIntersectionIssue`]s.
pub fn self_intersection_issues_summary(issues: &[SelfIntersectionIssue]) -> String {
    if issues.is_empty() {
        return "no self-intersections".to_string();
    }
    let crossings = issues
        .iter()
        .filter(|i| i.kind == SelfIntersectionKind::Crossing)
        .count();
    let overlaps = issues.len() - crossings;
    format!(
        "{} self-intersection issue(s): {crossings} crossing, {overlaps} coplanar overlap",
        issues.len()
    )
}

// ---------------------------------------------------------------------------
// Boundary structural analysis (BOPAlgo_ArgumentAnalyzer checks)
// ---------------------------------------------------------------------------

/// Is a face degenerate — a boundary polygon with (near-)zero area?
///
/// A face whose registered surface is not planar, or whose boundary cannot be
/// extracted, is conservatively reported as *not* degenerate (the check only
/// fires on faces it can measure).
pub fn face_is_degenerate(face: &Face, tol: f64) -> bool {
    let Some(pln) = face_plane_local(face) else {
        return false;
    };
    match face_polygon_local(face, &pln) {
        Some(poly) => polygon_area2d(&poly).abs() < tol.max(1e-12),
        None => false,
    }
}

/// Number of free (open-boundary) edges of a shape, within `tol`.
///
/// A free edge is an edge whose endpoint is used by exactly one edge end in its
/// wire — the signature of an open boundary. A closed box has none; a single
/// open face's wire edges all show up.
pub fn shape_free_edge_count(shape: &TopoShape, tol: f64) -> usize {
    crate::shhealing::free_edges(shape, tol).len()
}

/// Structured boundary analysis of a shape.
///
/// One report aggregates every check `BOPAlgo_ArgumentAnalyzer` runs on an
/// input/result before the boolean: shell count and closedness, free edges,
/// self-intersections, degenerate faces and small edges. A valid closed box
/// reports one closed shell and zeros everywhere else.
#[derive(Debug, Default, Clone)]
pub struct BoundaryAnalysisReport {
    /// Number of shells in the shape.
    pub shell_count: usize,
    /// Number of shells whose every boundary edge is used by exactly two faces.
    pub closed_shells: usize,
    /// Number of shells that are not closed.
    pub open_shells: usize,
    /// Number of free (open-boundary) edges.
    pub free_edges: usize,
    /// Number of self-intersecting non-adjacent face pairs.
    pub self_intersections: usize,
    /// Number of degenerate (near-zero-area) faces.
    pub degenerate_faces: usize,
    /// Number of edges shorter than `tol`.
    pub small_edges: usize,
}

/// Run the boundary analysis checks on `shape`.
pub fn analyze_boundary(shape: &TopoShape, tol: f64) -> BoundaryAnalysisReport {
    let tol = tol.max(1e-9);
    let shells: Vec<Shell> = shapes_of(shape, ShapeType::Shell).into_iter().map(Shell).collect();
    let closed_shells = shells.iter().filter(|s| shell_is_closed(s)).count();
    BoundaryAnalysisReport {
        shell_count: shells.len(),
        closed_shells,
        open_shells: shells.len() - closed_shells,
        free_edges: crate::shhealing::free_edges(shape, tol).len(),
        self_intersections: detect_self_intersections(shape, tol).edge_count,
        degenerate_faces: faces_of(shape).iter().filter(|f| face_is_degenerate(f, tol)).count(),
        small_edges: edges_of(shape)
            .iter()
            .filter(|e| crate::brep_measure::edge_length(e, 8) < tol)
            .count(),
    }
}

/// One-line diagnostic of a [`BoundaryAnalysisReport`].
pub fn boundary_analysis_summary(r: &BoundaryAnalysisReport) -> String {
    format!(
        "{} shell(s) ({} closed, {} open), {} free edge(s), {} self-intersection(s), {} degenerate face(s), {} small edge(s)",
        r.shell_count,
        r.closed_shells,
        r.open_shells,
        r.free_edges,
        r.self_intersections,
        r.degenerate_faces,
        r.small_edges
    )
}
