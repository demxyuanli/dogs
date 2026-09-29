//! N-ary / compound / degenerate boolean wrappers around [`crate::bop_builder::boolean`].
//! Split from `bop_builder.rs` without changing control flow.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::geom::polygon_ops::{point_in_polygon2d};
use occt_core::gp::{GpPnt};


use crate::abs::ShapeType;

use crate::builder::TopoBuilder;

use crate::shape::{Face, Shell, TopoShape};
use crate::shell_check::{shell_is_closed};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of, faces_of, shapes_of};


use crate::bop_builder_core::{empty_result, single_shape_result, BooleanResult, BoolOp};
use crate::bop_builder_planar::*;

/// Recursively flatten a compound into its non-compound sub-shapes.
///
/// A compound whose children are themselves compounds is fully flattened, so
/// nested results (e.g. a compound built by folding a multi-fuse) decompose
/// into their atomic solids/shells/faces.
pub(crate) fn expand_compound(s: &TopoShape) -> Vec<TopoShape> {
    if !s.is_compound() {
        return vec![s.clone()];
    }
    let mut out = Vec::new();
    for k in s.tshape.read().unwrap().children.clone() {
        out.extend(expand_compound(&k));
    }
    out
}

/// Split a `Compound` into its top-level sub-shapes (its direct children).
///
/// A non-compound shape is returned as a single-element slice. Nested
/// compounds are kept as-is at the top level (use [`boolean_multi`] or
/// [`expand_compound`] when full flattening is needed).
pub fn decompose_compound(shape: &TopoShape) -> Vec<TopoShape> {
    if !shape.is_compound() {
        return vec![shape.clone()];
    }
    shape
        .tshape
        .read()
        .unwrap()
        .children
        .clone()
}

/// Decompose a shape into its connected boundary components.
///
/// * a `Solid` → its shells, each wrapped back into a one-shell solid;
/// * a `Compound` → its top-level sub-shapes (recursively flattened);
/// * anything else → the shape itself.
///
/// This is the "result decomposition" step of `BOPAlgo_Builder`: after a
/// boolean, a result may hold several disconnected solids; this splits them
/// apart so callers can reason about each component.
pub fn shape_components(shape: &TopoShape) -> Vec<TopoShape> {
    match shape.shape_type() {
        ShapeType::Solid => {
            let bld = TopoBuilder::new();
            let mut comps = Vec::new();
            for sh in shapes_of(shape, ShapeType::Shell) {
                comps.push(bld.make_solid(&[Shell(sh)]).0);
            }
            if comps.is_empty() {
                vec![shape.clone()]
            } else {
                comps
            }
        }
        ShapeType::Compound => expand_compound(shape),
        _ => vec![shape.clone()],
    }
}

/// Does `s` carry any boundary faces? Empty compounds, null wires and other
/// degenerate inputs have none and are treated as the "empty" shape.
fn has_content(s: &TopoShape) -> bool {
    !faces_of(s).is_empty()
}

/// Is `s` usable as a solid operand? A `Solid` or a *closed* `Shell`; empty
/// shells/solids are rejected so the operands go to [`boolean_non_solid`]
/// (`BOPAlgo_BOP::CheckData` then decides, as in OCCT).
fn is_solid_input(s: &TopoShape) -> bool {
    if !has_content(s) {
        return false;
    }
    if s.is_solid() {
        return true;
    }
    if s.is_shell() {
        return shell_is_closed(&Shell(s.clone()));
    }
    false
}

/// Whether the bounding boxes of `a` and `b` overlap (touch counts).
fn bounding_boxes_overlap(a: &TopoShape, b: &TopoShape) -> bool {
    let ba = crate::bbox_from_geometry::shape_bbox(a);
    let bb = crate::bbox_from_geometry::shape_bbox(b);
    !ba.is_void() && !bb.is_void() && !ba.is_out_box(&bb)
}

/// Merge several result shapes into a single [`BooleanResult`].
///
/// Empty shapes are dropped. A single survivor becomes the result shape (and
/// its solid, when closed); multiple survivors are assembled into a compound.
/// This mirrors `BOPAlgo_Builder` returning a compound of the disconnected
/// result solids.
fn merge_shapes_result(shapes: Vec<TopoShape>, warnings: Vec<String>) -> BooleanResult {
    let kept: Vec<TopoShape> = shapes.into_iter().filter(has_content).collect();
    if kept.is_empty() {
        let mut r = empty_result(BoolOp::Fuse);
        r.warnings.extend(warnings);
        return r;
    }
    if kept.len() == 1 {
        let mut r = single_shape_result(&kept[0]);
        r.warnings.extend(warnings);
        return r;
    }
    let bld = TopoBuilder::new();
    let comp = bld.make_compound_of(&kept);
    let mut shells = Vec::new();
    let mut faces = Vec::new();
    for s in &kept {
        shells.extend(shapes_of(s, ShapeType::Shell).into_iter().map(Shell));
        faces.extend(faces_of(s));
    }
    BooleanResult { shape: comp.0, solid: None, shells, faces, warnings }
}

/// Fuse a list of shapes into a single result.
///
/// Overlapping arguments go through one `BOPAlgo_BOP` (`objects` = first,
/// `tools` = rest, one PaveFiller) — `BOPAlgo_BOP.cxx:106-134` requires both
/// groups, then `Perform` + `BuildShape`. Pairwise bbox folding is not in OCCT
/// and loses the N-way split. Mutually disjoint arguments stay a compound of
/// the unsplit inputs (no interferences, same as GF argument images).
fn fuse_components(shapes: &[TopoShape], tol: f64) -> Result<BooleanResult, String> {
    let components: Vec<TopoShape> = shapes.iter().filter(|s| has_content(s)).cloned().collect();
    if components.is_empty() {
        return Ok(empty_result(BoolOp::Fuse));
    }
    if components.len() == 1 {
        return Ok(single_shape_result(&components[0]));
    }
    let any_overlap = components.iter().enumerate().any(|(i, a)| {
        components[i + 1..]
            .iter()
            .any(|b| bounding_boxes_overlap(a, b))
    });
    if !any_overlap {
        return Ok(merge_shapes_result(components, vec![]));
    }
    let shape = crate::bop_builder2::builder_bop_with_fuzzy(
        std::slice::from_ref(&components[0]),
        &components[1..],
        crate::bop_builder2::BoolOp2::Fuse,
        tol,
    )?;
    let shape = if shape.is_compound() {
        let solids = shapes_of(&shape, ShapeType::Solid);
        if solids.len() == 1 {
            solids[0].clone()
        } else {
            shape
        }
    } else {
        shape
    };
    Ok(single_shape_result(&shape))
}

/// Top-level dispatch between the exact boolean paths.
///
/// * either operand is a compound → [`boolean_compound`] (expand, per-part);
/// * either operand is not a solid → [`boolean_non_solid`] (`BOPAlgo_BOP`'s
///   `CheckData` decides legality — OCCT has no best-effort layer);
/// * otherwise the full curved/planar boolean dispatcher.
pub(crate) fn boolean_dispatch(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    // `BOPAlgo_BOP::Perform` order: `CheckData` raises `BOPAlgo_AlertEmptyShape`
    // (`BOPAlgo_BOP.cxx:162-167`) and `TreatEmptyShape` (`:214-322`) then fixes the
    // result **before** anything else happens — in particular before any compound
    // handling. The port used to test `is_compound()` first, so an operand that is
    // an empty compound bypassed `TreatEmptyShape` and fell into the per-part path
    // (with no warning and, for two empty operands, a `too few arguments` error
    // instead of OCCT's empty result).
    let empty_objects = crate::bop_bop::is_empty_shape(a);
    let empty_tools = crate::bop_bop::is_empty_shape(b);
    if empty_objects || empty_tools {
        if let Some(chosen) = crate::bop_bop::treat_empty_shape(
            std::slice::from_ref(a),
            std::slice::from_ref(b),
            to_bool_op2(op),
        ) {
            let mut warnings: Vec<String> = Vec::new();
            if empty_objects {
                warnings.push("BOPAlgo_AlertEmptyShape (objects)".into());
            }
            if empty_tools {
                warnings.push("BOPAlgo_AlertEmptyShape (tools)".into());
            }
            let bld = TopoBuilder::new();
            let shape: TopoShape = bld.make_compound_of(&chosen).into();
            let mut r = single_shape_result(&shape);
            r.warnings.extend(warnings);
            return Ok(r);
        }
    }
    if a.is_compound() || b.is_compound() {
        return boolean_compound(a, b, op, tol);
    }
    if !is_solid_input(a) || !is_solid_input(b) {
        return boolean_non_solid(a, b, op, tol);
    }
    // T-41: OCCT has no `bop_curved` mesh/voxel engine; two non-compound
    // solids go through `BOPAlgo_BOP` like every other operand pair. The old
    // `curved_boolean_full` routing (UV-sampling / voxel boolean) had no OCCT
    // counterpart (audit A5 / T-41) and is removed here; its module-level tests
    // that call it directly are unaffected.
    crate::bop_builder::boolean_via_bopalgo(a, b, op, tol)
}

/// `BOPAlgo_BOP` on operands that are not both usable solids (faces, wires,
/// open shells, empty shapes).
///
/// OCCT has no best-effort "degenerate" branch: `BOPAlgo_BOP::CheckData`
/// (`BOPAlgo_BOP.cxx:140-210`) is the only gate. It skips an empty argument
/// with the `BOPAlgo_AlertEmptyShape` warning (`cxx:162-167`), lets the empty
/// group take the other group's dimension (`cxx:203-209`), and rejects a
/// `FUSE` of unequal dimensions (`cxx:181-186`, `:193-195`) or a `CUT` whose
/// objects are of larger dimension than the tools (`cxx:194`) with
/// `BOPAlgo_AlertBOPNotAllowed`; `COMMON` accepts any dimensions. The port's
/// `bop_bop::check_data` is that check and
/// `bop_builder2::builder_bop_with_fuzzy` is `BOPAlgo_BOP::Perform` +
/// `BuildShape`, so the operands are handed over unchanged.
fn boolean_non_solid(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let mut warnings: Vec<String> = Vec::new();
    if crate::bop_bop::is_empty_shape(a) {
        warnings.push("BOPAlgo_AlertEmptyShape (objects)".into());
    }
    if crate::bop_bop::is_empty_shape(b) {
        warnings.push("BOPAlgo_AlertEmptyShape (tools)".into());
    }
    let shape = crate::bop_builder2::builder_bop_with_fuzzy(
        std::slice::from_ref(a),
        std::slice::from_ref(b),
        to_bool_op2(op),
        tol,
    )?;
    let mut r = single_shape_result(&shape);
    r.warnings.extend(warnings);
    Ok(r)
}

/// The port's two operation enums (`bop_builder_core::BoolOp` and
/// `bop_builder2::BoolOp2`, both `Fuse`/`Cut`/`Common`).
pub(crate) fn to_bool_op2(op: BoolOp) -> crate::bop_builder2::BoolOp2 {
    match op {
        BoolOp::Fuse => crate::bop_builder2::BoolOp2::Fuse,
        BoolOp::Cut => crate::bop_builder2::BoolOp2::Cut,
        BoolOp::Common => crate::bop_builder2::BoolOp2::Common,
    }
}

/// Apply a boolean operation to a *list* of shapes (N-ary boolean).
///
/// `BOPAlgo_Builder`'s `Build`-with-many-arguments entry point. Handles:
///
/// * an empty list → an empty result (an empty compound);
/// * a single shape → the shape itself (identity);
/// * `Fuse` of mutually disjoint shapes → a `Compound` keeping every input
///   (not an error);
/// * `Fuse` of overlapping shapes → one `BOPAlgo_BOP` on all arguments;
/// * `Cut` → fold `a0 − a1 − a2 − …`;
/// * `Common` → fold `a0 ∩ a1 ∩ a2 ∩ …`.
///
/// Compounds among the inputs are flattened first (a compound is a set of
/// components, so the operation distributes over it).
pub fn boolean_multi(shapes: &[TopoShape], op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    if shapes.is_empty() {
        return Ok(empty_result(op));
    }
    let flat: Vec<TopoShape> = shapes.iter().flat_map(expand_compound).collect();
    if flat.is_empty() {
        return Ok(empty_result(op));
    }
    if flat.len() == 1 {
        return Ok(single_shape_result(&flat[0]));
    }
    match op {
        BoolOp::Fuse => fuse_components(&flat, tol),
        BoolOp::Cut | BoolOp::Common => {
            // Fold left. Empty and non-solid intermediates are handled by
            // `BOPAlgo_BOP` itself: an empty group takes the other group's
            // dimension (`BOPAlgo_BOP.cxx:203-209`), so a cut against an empty
            // shape leaves the accumulator unchanged and an empty common stays
            // empty — no special case is needed here.
            let mut acc = flat[0].clone();
            let mut warnings: Vec<String> = Vec::new();
            for s in &flat[1..] {
                let r = boolean_dispatch(&acc, s, op, tol)?;
                warnings.extend(r.warnings);
                acc = r.shape;
            }
            let mut res = single_shape_result(&acc);
            res.warnings.extend(warnings);
            Ok(res)
        }
    }
}

/// Boolean between shapes where either operand is a `Compound`.
///
/// Expands the compound(s) and applies the operation per component:
///
/// * `Fuse` — the union of *all* components of `a` and `b` (compounds are
///   unions of their parts, so `(A1 ∪ A2) ∪ (B1 ∪ B2)` is the full union);
/// * `Cut` — every component of `b` is cut out of every component of `a`;
/// * `Common` — the pairwise intersections `ai ∩ bj` are merged.
///
/// Solid operands are delegated to the exact boolean
/// (`boolean`/`curved_boolean_full`).
pub fn boolean_compound(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let subs_a = expand_compound(a);
    let subs_b = expand_compound(b);
    let mut warnings: Vec<String> = Vec::new();
    match op {
        BoolOp::Fuse => {
            let mut all = subs_a;
            all.extend(subs_b);
            fuse_components(&all, tol)
        }
        BoolOp::Cut => {
            let mut results: Vec<TopoShape> = Vec::new();
            for sa in &subs_a {
                let mut cur = sa.clone();
                for sb in &subs_b {
                    let r = boolean_dispatch(&cur, sb, BoolOp::Cut, tol)?;
                    warnings.extend(r.warnings);
                    cur = r.shape;
                }
                results.push(cur);
            }
            Ok(merge_shapes_result(results, warnings))
        }
        BoolOp::Common => {
            let mut results: Vec<TopoShape> = Vec::new();
            for sa in &subs_a {
                for sb in &subs_b {
                    let r = boolean_dispatch(sa, sb, BoolOp::Common, tol)?;
                    warnings.extend(r.warnings);
                    if has_content(&r.shape) {
                        results.push(r.shape);
                    }
                }
            }
            Ok(merge_shapes_result(results, warnings))
        }
    }
}

/// Report produced by [`detect_self_intersections`].
pub struct SelfIntersectionReport {
    /// Whether any pair of non-adjacent faces was found to intersect.
    pub found: bool,
    /// Number of intersecting face pairs (each pair counts as one
    /// self-intersection "edge" of the defect).
    pub edge_count: usize,
    /// Sample points on the found intersection curves.
    pub points: Vec<GpPnt>,
}

/// Do the two coplanar faces' boundary polygons overlap in area?
pub(crate) fn faces_polygon_overlap(f1: &Face, f2: &Face, tol: f64) -> bool {
    let pln = match face_plane_local(f1) {
        Some(p) => p,
        None => return false,
    };
    let (Some(poly1), Some(poly2)) = (face_polygon_local(f1, &pln), face_polygon_local(f2, &pln)) else {
        return false;
    };
    for v in &poly1 {
        if point_in_polygon2d(&poly2, v) {
            return true;
        }
    }
    for v in &poly2 {
        if point_in_polygon2d(&poly1, v) {
            return true;
        }
    }
    let _ = tol;
    false
}

/// Check a shape for self-intersecting faces.
///
/// A self-intersection is a pair of *non-adjacent* faces (faces that do not
/// share a boundary edge) whose underlying surfaces cross or overlap:
///
/// * transversal crossing → the sampled `SurfaceIntersection::Curves`;
/// * coplanar overlap → `SurfaceIntersection::Coincident` with overlapping
///   face polygons.
///
/// Adjacent faces are skipped because sharing a boundary edge is the normal
/// (and legal) way two faces of a solid meet. A valid closed box therefore
/// reports `found == false` (its only non-adjacent pairs are parallel faces),
/// while a shell built from two crossing faces reports `found == true`.
pub fn detect_self_intersections(shape: &TopoShape, tol: f64) -> SelfIntersectionReport {
    let tol = tol.max(1e-9);
    let faces = faces_of(shape);
    let mut report = SelfIntersectionReport { found: false, edge_count: 0, points: Vec::new() };
    if faces.len() < 2 {
        return report;
    }
    // Per-face boundary-edge identities, for the adjacency test. Two faces are
    // adjacent iff they share the same `TShape` edge.
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();
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
                    if !pts.is_empty() {
                        report.found = true;
                        report.edge_count += 1;
                        report.points.extend(pts);
                    }
                }
                crate::intpatch::SurfaceIntersection::Coincident => {
                    if faces_polygon_overlap(&faces[i], &faces[j], tol) {
                        report.found = true;
                        report.edge_count += 1;
                    }
                }
                crate::intpatch::SurfaceIntersection::None => {}
            }
        }
    }
    report
}

/// Structural validation of a boolean result.
///
/// Returns a list of human-readable issue strings (empty when the result is
/// clean): every shell must be closed (every boundary edge used by exactly two
/// faces — `shell_is_closed`) and the result must not be self-intersecting.
pub(crate) fn validate_boolean_result(r: &BooleanResult, tol: f64) -> Vec<String> {
    let mut issues: Vec<String> = Vec::new();
    if r.shells.is_empty() {
        issues.push("result has no shell".into());
    }
    for (i, sh) in r.shells.iter().enumerate() {
        if !shell_is_closed(sh) {
            issues.push(format!("result shell {i} is not closed"));
        }
    }
    let si = detect_self_intersections(&r.shape, tol);
    if si.found {
        issues.push(format!("result self-intersects ({} face pairs)", si.edge_count));
    }
    issues
}

/// Run a boolean and validate the result, re-computing at a larger tolerance
/// when validation finds issues.
///
/// The exact boolean is tolerance-sensitive near coincident faces; when the
/// first attempt fails validation (an open shell, a self-intersection), the
/// operation is re-run at `tol × 10` and `tol × 100`. The best (fewest issues)
/// result is returned, with the validation issues attached as warnings.
pub fn boolean_with_check(a: &TopoShape, b: &TopoShape, op: BoolOp, tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let mut result = boolean_dispatch(a, b, op, tol)?;
    let mut issues = validate_boolean_result(&result, tol);
    if !issues.is_empty() {
        for &scale in &[10.0f64, 100.0] {
            let nt = tol * scale;
            let candidate = boolean_dispatch(a, b, op, nt)?;
            let cand_issues = validate_boolean_result(&candidate, nt);
            if cand_issues.len() < issues.len() {
                result = candidate;
                issues = cand_issues;
                result.warnings.push(format!("recomputed at tol {nt:.1e} after validation (original {tol:.1e})"));
            }
            if issues.is_empty() {
                break;
            }
        }
    }
    result.warnings.extend(issues);
    Ok(result)
}


/// One-line diagnostic of a [`BooleanResult`].
///
/// Reports the result shape type, whether it is a closed solid, its face
/// count, and how many warnings it carries.
pub fn boolean_result_summary(r: &BooleanResult) -> String {
    let st = r.shape.shape_type().to_str();
    let kind = if r.solid.is_some() { "closed solid" } else { "open shell/compound" };
    let faces = faces_of(&r.shape).len();
    let warns = if r.warnings.is_empty() {
        "no warnings".to_string()
    } else {
        format!("{} warning(s): {}", r.warnings.len(), r.warnings.join("; "))
    };
    format!("{st} ({kind}) {faces} faces, {warns}")
}

/// Apply a *sequence* of boolean operations to a list of shapes.
///
/// The operations are applied left to right, each between the running
/// accumulator and the next shape: `((a0 op0 a1) op1 a2) op2 a3 …`. Unlike
/// [`boolean_multi`], every step may use a different operation, so a CSG tree
/// flattened into the alternating form `shape, op, shape, op, shape, …` can be
/// evaluated directly. Extra shapes (beyond `ops.len() + 1`) are ignored; a
/// missing operation for a remaining shape stops the fold.
pub fn boolean_fold(shapes: &[TopoShape], ops: &[BoolOp], tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    if shapes.is_empty() {
        return Ok(empty_result(BoolOp::Fuse));
    }
    let mut acc = shapes[0].clone();
    let mut warnings: Vec<String> = Vec::new();
    let steps = shapes.len().saturating_sub(1).min(ops.len());
    for i in 0..steps {
        let r = boolean_dispatch(&acc, &shapes[i + 1], ops[i], tol)?;
        warnings.extend(r.warnings);
        acc = r.shape;
    }
    let mut res = single_shape_result(&acc);
    res.warnings.extend(warnings);
    Ok(res)
}

/// Subtract every shape in `cuts` from `a`, one after the other (folded Cut).
///
/// Equivalent to `a − c1 − c2 − … − cn`. Each step routes through
/// [`boolean_dispatch`], so compounds and degenerate operands are handled.
/// The result is the final solid (or an empty compound when `a` is fully
/// removed).
pub fn boolean_cut_many(a: &TopoShape, cuts: &[TopoShape], tol: f64) -> Result<BooleanResult, String> {
    let tol = tol.max(1e-9);
    let mut acc = a.clone();
    let mut warnings: Vec<String> = Vec::new();
    for c in cuts {
        let r = boolean_dispatch(&acc, c, BoolOp::Cut, tol)?;
        warnings.extend(r.warnings);
        acc = r.shape;
    }
    let mut res = single_shape_result(&acc);
    res.warnings.extend(warnings);
    Ok(res)
}

/// Intersect a list of shapes (folded Common).
///
/// Equivalent to `a0 ∩ a1 ∩ a2 ∩ …`. Degenerate intermediates (an empty
/// common) stay empty through the fold, so the final result is empty as soon
/// as any pair is disjoint.
pub fn boolean_common_many(shapes: &[TopoShape], tol: f64) -> Result<BooleanResult, String> {
    boolean_multi(shapes, BoolOp::Common, tol)
}

/// Pairs of face indices (into [`faces_of`]) that share a boundary edge.
///
/// Two faces are *adjacent* when they reference the same `TShape` edge. This
/// is the raw face-adjacency graph of a boundary, useful for connectivity
/// analysis and for understanding where a self-intersection check skips.
pub fn face_adjacency(shape: &TopoShape) -> Vec<(usize, usize)> {
    let faces = faces_of(shape);
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if face_edges[i].iter().any(|e| face_edges[j].contains(e)) {
                pairs.push((i, j));
            }
        }
    }
    pairs
}

/// Split a shape into its edge-connected boundary components.
///
/// Two faces belong to the same component when they are connected through a
/// chain of shared boundary edges (the face-adjacency graph). A single closed
/// solid has exactly one component; a compound of disjoint solids yields one
/// component per solid. Each component is rebuilt as a one-shell solid, so the
/// result is the "decomposed" form of a boolean result (`BOPAlgo_Builder`
/// returns such disconnected solids inside a compound).
///
/// Faces that are topologically isolated (an open face with no shared edges)
/// each become their own component.
pub fn connected_components(shape: &TopoShape, tol: f64) -> Vec<TopoShape> {
    let _ = tol;
    let faces = faces_of(shape);
    if faces.len() <= 1 {
        return vec![shape.clone()];
    }
    let face_edges: Vec<HashSet<usize>> = faces
        .iter()
        .map(|f| edges_of(&f.0).into_iter().map(|e| Arc::as_ptr(&e.0.tshape) as usize).collect())
        .collect();

    // Union-find over faces (edge-sharing ⇒ same component).
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            parent[r] = parent[parent[r]];
            r = parent[r];
        }
        r
    }
    fn unite(parent: &mut Vec<usize>, a: usize, b: usize) {
        let (ra, rb) = (find(parent, a), find(parent, b));
        if ra != rb {
            parent[ra] = rb;
        }
    }

    let mut parent: Vec<usize> = (0..faces.len()).collect();
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if face_edges[i].iter().any(|e| face_edges[j].contains(e)) {
                unite(&mut parent, i, j);
            }
        }
    }

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..faces.len() {
        groups.entry(find(&mut parent, i)).or_default().push(i);
    }

    let bld = TopoBuilder::new();
    let mut comps: Vec<TopoShape> = Vec::new();
    for (_, idx) in groups {
        let fs: Vec<Face> = idx.iter().map(|&i| faces[i].clone()).collect();
        let shell = bld.make_shell(&fs);
        comps.push(bld.make_solid(&[shell]).0);
    }
    comps
}
