//! Boolean result analysis helpers. Split from `bop_builder.rs`.








use crate::abs::ShapeType;


use crate::builder::TopoBuilder;

use crate::shape::{Edge, Face, Solid, TopoShape};
use crate::shell_check::{shell_is_closed};

use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, shapes_of, vertex_position, vertices_of,
    wires_of_face,
};


use crate::bop_builder_core::{single_shape_result, BoolOp, BooleanResult};
use crate::bop_builder_dispatch::{
    boolean_multi, boolean_with_check, connected_components, decompose_compound,
    detect_self_intersections, expand_compound, validate_boolean_result,
};
use crate::bop_builder_repair::{
    analyze_boundary, boolean_repaired, boolean_split_result, decompose_multi_result,
    face_connectivity_groups, face_is_degenerate, repair_boolean_result,
    repair_self_intersections, shape_boundary_components, split_connected_boundaries,
    topologically_clean, BoundaryAnalysisReport, MultiResult, RepairResult, SelfIntersectionIssue,
    SelfIntersectionKind,
};

pub fn boolean_result_validate(r: &BooleanResult, tol: f64) -> Vec<String> {
    validate_boolean_result(r, tol)
}

/// Convenience wrapper: fuse every shape in `shapes` (see [`boolean_multi`]).
///
/// Kept as a named entry point so callers do not have to spell out
/// `BoolOp::Fuse`; behaves identically to `boolean_multi(shapes, BoolOp::Fuse,
/// tol)`.
pub fn boolean_fuse_all(shapes: &[TopoShape], tol: f64) -> Result<BooleanResult, String> {
    boolean_multi(shapes, BoolOp::Fuse, tol)
}

/// Bounding box of a boolean result, or `None` when the result is empty.
///
/// Delegates to [`crate::bbox_from_geometry::shape_bbox`] on the result shape,
/// so compounds (multi-solid results) are covered by the compound bbox.
pub fn boolean_result_bbox(r: &BooleanResult) -> Option<occt_core::bnd::BndBox> {
    let bb = crate::bbox_from_geometry::shape_bbox(&r.shape);
    if bb.is_void() {
        None
    } else {
        Some(bb)
    }
}

/// Normalize a boolean result into its "assembled" form.
///
/// When a result carries more than one shell (a multi-solid Cut/Common result
/// that was not merged into a compound), wrap each closed shell into its own
/// solid and return a compound of them; a single-shell result is returned
/// unchanged. This mirrors the final assembly step of `BOPAlgo_Builder`, which
/// hands the caller a compound of the disconnected result solids.
pub fn normalize_result(r: &BooleanResult) -> BooleanResult {
    if r.shape.is_compound() || r.shells.len() <= 1 {
        return r.clone();
    }
    let bld = TopoBuilder::new();
    let solids: Vec<TopoShape> = r
        .shells
        .iter()
        .filter(|sh| shell_is_closed(sh))
        .map(|sh| bld.make_solid(&[sh.clone()]).0)
        .collect();
    if solids.is_empty() {
        return r.clone();
    }
    let shape = if solids.len() == 1 {
        solids[0].clone()
    } else {
        bld.make_compound_of(&solids).0
    };
    let solid = Solid::wrap(shape.clone());
    BooleanResult {
        shape,
        solid,
        shells: r.shells.clone(),
        faces: r.faces.clone(),
        warnings: r.warnings.clone(),
    }
}
pub fn repair_result_summary(r: &RepairResult) -> String {
    let mut out = format!("fixed {} face(s), removed {}", r.fixed_faces, r.removed_faces);
    if !r.warnings.is_empty() {
        out.push_str(&format!("; {} warning(s): {}", r.warnings.len(), r.warnings.join("; ")));
    }
    out
}

// ---------------------------------------------------------------------------
// Multi-result helpers
// ---------------------------------------------------------------------------

/// One-line diagnostic of a [`MultiResult`].
pub fn multi_result_summary(m: &MultiResult) -> String {
    format!(
        "{} shape(s): {} solid, {} shell, {} compound",
        m.shapes.len(),
        m.solids.len(),
        m.shells.len(),
        m.compounds.len()
    )
}

/// Total volume of every component of a [`MultiResult`].
///
/// Volumes are measured from each component's tessellation (closed solids only
/// contribute real volume; open shells contribute ~0).
pub fn multi_result_total_volume(m: &MultiResult) -> f64 {
    m.shapes
        .iter()
        .map(|s| crate::brep_gprop::volume(s, 0.02))
        .sum()
}

/// Volume of every component of a [`MultiResult`], in component order.
pub fn component_volumes(m: &MultiResult) -> Vec<f64> {
    m.shapes.iter().map(|s| crate::brep_gprop::volume(s, 0.02)).collect()
}

/// The component with the largest volume, or `None` for an empty result.
pub fn largest_component(m: &MultiResult) -> Option<TopoShape> {
    let mut best: Option<(f64, TopoShape)> = None;
    for s in &m.shapes {
        let v = crate::brep_gprop::volume(s, 0.02);
        if best.as_ref().map_or(true, |(bv, _)| v > *bv) {
            best = Some((v, s.clone()));
        }
    }
    best.map(|(_, s)| s)
}

/// The component with the smallest volume, or `None` for an empty result.
pub fn smallest_component(m: &MultiResult) -> Option<TopoShape> {
    let mut best: Option<(f64, TopoShape)> = None;
    for s in &m.shapes {
        let v = crate::brep_gprop::volume(s, 0.02);
        if best.as_ref().map_or(true, |(bv, _)| v < *bv) {
            best = Some((v, s.clone()));
        }
    }
    best.map(|(_, s)| s)
}

/// Bucket sizes of a [`MultiResult`] as `(solids, shells, compounds)`.
pub fn component_counts_by_type(m: &MultiResult) -> (usize, usize, usize) {
    (m.solids.len(), m.shells.len(), m.compounds.len())
}

// ---------------------------------------------------------------------------
// Composite boolean pipelines
// ---------------------------------------------------------------------------

/// Run a boolean and return the disconnected result pieces as a plain vector.
///
/// Convenience wrapper over [`boolean_split_result`] that drops the
/// per-type buckets and returns just the shapes.
pub fn boolean_components(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<Vec<TopoShape>, String> {
    Ok(boolean_split_result(a, b, op, tol)?.shapes)
}

/// Run a boolean and return the volume of every disconnected result piece.
///
/// Useful for sanity-checking a Cut that should leave several pieces: the
/// returned volumes can be summed and compared against the analytic result.
pub fn boolean_component_volumes(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<Vec<f64>, String> {
    Ok(component_volumes(&boolean_split_result(a, b, op, tol)?))
}

/// Run a boolean, repair self-intersections, and decompose the repaired result.
pub fn boolean_repaired_multi(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<MultiResult, String> {
    let r = boolean_repaired(a, b, op, tol)?;
    Ok(decompose_multi_result(&r))
}

/// Run a boolean, repair self-intersections and topologically clean the result.
///
/// Returns the single repaired+cleaned result shape (a compound when the result
/// is disconnected).
pub fn boolean_repair_and_clean(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<TopoShape, String> {
    let r = boolean_repaired(a, b, op, tol)?;
    Ok(topologically_clean(&r.shape, tol))
}

// ---------------------------------------------------------------------------
// Repair iteration & single-defect tooling
// ---------------------------------------------------------------------------

/// Iterate [`repair_self_intersections`] until the boundary is clean.
///
/// A repair pass can introduce a *new* self-intersection when a face is split
/// near another crossing face, so `BOPAlgo`-style repair is applied to
/// fixpoint. Runs at most `max_iter` passes; the counters and warnings are
/// accumulated across passes.
pub fn repair_self_intersections_loop(shape: &TopoShape, tol: f64, max_iter: usize) -> RepairResult {
    let mut current = shape.clone();
    let mut fixed_total = 0usize;
    let mut removed_total = 0usize;
    let mut warnings: Vec<String> = Vec::new();
    for _ in 0..max_iter.max(1) {
        let rep = match repair_self_intersections(&current, tol) {
            Ok(r) => r,
            Err(e) => {
                warnings.push(format!("repair pass failed: {e}"));
                break;
            }
        };
        fixed_total += rep.fixed_faces;
        removed_total += rep.removed_faces;
        warnings.extend(rep.warnings);
        current = rep.repaired;
        if rep.fixed_faces == 0 && rep.removed_faces == 0 {
            break;
        }
    }
    RepairResult {
        repaired: current,
        fixed_faces: fixed_total,
        removed_faces: removed_total,
        warnings,
    }
}

/// List the degenerate (near-zero-area) faces of a shape.
pub fn shape_degenerate_faces(shape: &TopoShape, tol: f64) -> Vec<Face> {
    faces_of(shape)
        .into_iter()
        .filter(|f| face_is_degenerate(f, tol))
        .collect()
}

/// List the edges of a shape shorter than `tol`.
pub fn shape_small_edges(shape: &TopoShape, tol: f64) -> Vec<Edge> {
    let tol = tol.max(1e-9);
    edges_of(shape)
        .into_iter()
        .filter(|e| crate::brep_measure::edge_length(e, 8) < tol)
        .collect()
}

/// The face-index connectivity groups of a shape's boundary.
///
/// Returns the groups of [`face_connectivity_groups`] for every face under
/// `shape`, so callers can see how the boundary splits without rebuilding it.
pub fn boundary_connectivity(shape: &TopoShape, tol: f64) -> Vec<Vec<usize>> {
    let _ = tol;
    face_connectivity_groups(&faces_of(shape))
}

/// Structural validity issues of a shape, as human-readable strings.
///
/// A port of the `BOPAlgo_ArgumentAnalyzer` validity pass: invalid topology
/// structure, open shells, free edges, self-intersections and degenerate
/// faces. An empty list means the shape passed every check.
pub fn check_shape_validity(shape: &TopoShape, tol: f64) -> Vec<String> {
    let mut issues: Vec<String> = Vec::new();
    if !crate::topo_tools_full::structure_is_valid(shape) {
        issues.push("invalid topology structure (a sub-shape has an illegal parent type)".into());
    }
    let analysis = analyze_boundary(shape, tol);
    if analysis.open_shells > 0 {
        issues.push(format!("{} open shell(s)", analysis.open_shells));
    }
    if analysis.free_edges > 0 {
        issues.push(format!("{} free edge(s)", analysis.free_edges));
    }
    if analysis.self_intersections > 0 {
        issues.push(format!("{} self-intersection(s)", analysis.self_intersections));
    }
    if analysis.degenerate_faces > 0 {
        issues.push(format!("{} degenerate face(s)", analysis.degenerate_faces));
    }
    issues
}

// ---------------------------------------------------------------------------
// Component classification & selection
// ---------------------------------------------------------------------------

/// Broad shape kind of a result component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComponentKind {
    /// A solid (closed boundary).
    Solid,
    /// An open shell.
    Shell,
    /// A compound (still grouped result).
    Compound,
    /// A bare face.
    Face,
    /// A bare wire.
    Wire,
    /// Any other shape type.
    Other,
}

/// Classify a shape into a [`ComponentKind`].
pub fn component_kind(shape: &TopoShape) -> ComponentKind {
    match shape.shape_type() {
        ShapeType::Solid => ComponentKind::Solid,
        ShapeType::Shell => ComponentKind::Shell,
        ShapeType::Compound => ComponentKind::Compound,
        ShapeType::Face => ComponentKind::Face,
        ShapeType::Wire => ComponentKind::Wire,
        _ => ComponentKind::Other,
    }
}

/// Human-readable label of a [`ComponentKind`].
pub fn component_kind_label(kind: ComponentKind) -> &'static str {
    match kind {
        ComponentKind::Solid => "solid",
        ComponentKind::Shell => "shell",
        ComponentKind::Compound => "compound",
        ComponentKind::Face => "face",
        ComponentKind::Wire => "wire",
        ComponentKind::Other => "other",
    }
}

/// The components of a [`MultiResult`] that match `kind`.
pub fn multi_result_filter_by_kind(m: &MultiResult, kind: ComponentKind) -> Vec<TopoShape> {
    m.shapes
        .iter()
        .filter(|s| component_kind(s) == kind)
        .cloned()
        .collect()
}

/// The components of a [`MultiResult`], sorted by volume descending.
pub fn multi_result_sorted_by_volume(m: &MultiResult) -> Vec<TopoShape> {
    let mut sorted = m.shapes.clone();
    sorted.sort_by(|a, b| {
        crate::brep_gprop::volume(b, 0.02)
            .partial_cmp(&crate::brep_gprop::volume(a, 0.02))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    sorted
}

// ---------------------------------------------------------------------------
// Composite pipelines (check + repair, fuse-all, area)
// ---------------------------------------------------------------------------

/// Run a boolean with validation/tolerance-escalation, then repair the result.
///
/// Combines [`boolean_with_check`] (which re-runs the boolean at larger
/// tolerances when structural validation fails) with
/// [`repair_self_intersections`]. A repaired result carries a warning naming
/// the number of faces fixed/removed.
pub fn boolean_repaired_with_check(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let r = boolean_with_check(a, b, op, tol)?;
    Ok(repair_boolean_result(r, tol))
}

/// Fuse every shape in `shapes` and decompose the result into its pieces.
///
/// Convenience wrapper over [`boolean_multi`] + [`decompose_multi_result`]: a
/// Fuse of several disjoint boxes yields one [`MultiResult`] shape per box.
pub fn boolean_fuse_all_components(shapes: &[TopoShape], tol: f64) -> Result<MultiResult, String> {
    let r = boolean_multi(shapes, BoolOp::Fuse, tol)?;
    Ok(decompose_multi_result(&r))
}

/// Run a boolean, repair the result and return the repair report alongside it.
///
/// Like [`boolean_repaired`], but also exposes the underlying [`RepairResult`]
/// so callers can inspect exactly which faces were fixed/removed.
pub fn boolean_repaired_report(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<(BooleanResult, RepairResult), String> {
    let tol = tol.max(1e-9);
    // T-41: faithful `BOPAlgo_BOP` dispatch (see `boolean_repaired`).
    let r = crate::bop_builder_dispatch::boolean_dispatch(a, b, op, tol)?;
    let rep = repair_self_intersections(&r.shape, tol)?;
    if rep.fixed_faces == 0 && rep.removed_faces == 0 {
        return Ok((r, rep));
    }
    let mut rr = single_shape_result(&rep.repaired);
    rr.warnings.extend(r.warnings);
    rr.warnings.extend(rep.warnings.clone());
    rr.warnings.push(format!(
        "boolean result had self-intersections; repair fixed {} face(s), removed {}",
        rep.fixed_faces, rep.removed_faces
    ));
    Ok((rr, rep))
}

/// Surface area of a shape, from its tessellation (`BRepGProp::SurfaceProperties`).
pub fn boundary_area(shape: &TopoShape, deflection: f64) -> f64 {
    crate::brep_gprop::surface_area(shape, deflection.max(0.02))
}

/// Total surface area of every component of a [`MultiResult`].
pub fn multi_result_total_area(m: &MultiResult) -> f64 {
    m.shapes.iter().map(|s| crate::brep_gprop::surface_area(s, 0.02)).sum()
}

/// Volume of a single component (closed solid), or ~0 for open shells.
pub fn component_volume(shape: &TopoShape, deflection: f64) -> f64 {
    crate::brep_gprop::volume(shape, deflection)
}

// ---------------------------------------------------------------------------
// Boolean-result analysis helpers
// ---------------------------------------------------------------------------

/// The structural boundary analysis of a boolean result's shape.
pub fn boolean_result_boundary_report(r: &BooleanResult, tol: f64) -> BoundaryAnalysisReport {
    analyze_boundary(&r.shape, tol)
}

/// Number of disconnected pieces of a boolean result.
pub fn boolean_result_component_count(r: &BooleanResult) -> usize {
    decompose_multi_result(r).shapes.len()
}

/// Absolute difference between a boolean result's volume and an expectation.
///
/// `expected` is normally the analytic volume (e.g. `vol_a + vol_b` for a
/// Fuse); the delta is the mesh-measured deviation.
pub fn boolean_result_volume_delta(r: &BooleanResult, expected: f64) -> f64 {
    (crate::brep_gprop::volume(&r.shape, 0.02) - expected).abs()
}

/// Are two boolean results' volumes equal within a relative tolerance?
///
/// Uses the larger of the two volumes as the scale, so a 1% deviation on a
/// unit-volume result and on a 100-volume result are both "close".
pub fn boolean_results_close(a: &BooleanResult, b: &BooleanResult, rel_tol: f64) -> bool {
    let va = crate::brep_gprop::volume(&a.shape, 0.02);
    let vb = crate::brep_gprop::volume(&b.shape, 0.02);
    let scale = va.max(vb).max(1e-9);
    (va - vb).abs() <= rel_tol * scale
}

/// The broad shape kind of a boolean result's shape.
pub fn boolean_result_shape_kind(r: &BooleanResult) -> ComponentKind {
    component_kind(&r.shape)
}

/// Does a [`MultiResult`] contain at least one component of `kind`?
pub fn multi_result_contains_kind(m: &MultiResult, kind: ComponentKind) -> bool {
    m.shapes.iter().any(|s| component_kind(s) == kind)
}

/// Total volume of a boolean result's connected components.
pub fn boolean_result_total_volume(r: &BooleanResult) -> f64 {
    multi_result_total_volume(&decompose_multi_result(r))
}
