//! `XSAlgo_ShapeProcessor::ProcessShape` for the IGES default operator
//! `FixShape` (`IGESToBRep_Actor.cxx:204-228`, `IGESData.cxx:109`).
//!
//! `ShapeFix_Shape::Perform` (`ShapeFix_Shape.cxx:83-302`) walks the shape
//! with one `ShapeBuild_ReShape` and finishes with `Context()->Apply`
//! (`cxx:257`). The face arm calls `ShapeFix_Face::Perform`
//! (`cxx:193-210`). This file runs the steps of that face `Perform` that
//! are already ported, in the same order the STEP reader uses.
//!
//! A free wire has a null face, so `IsReady` is false and the face fixes
//! inside `ShapeFix_Wire::Perform` return immediately.

use crate::abs::{Orientation, ShapeType};
use crate::brep_class3d::SolidClassifier;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, Shell, TopoShape, Wire};
use crate::shape_fix_compose_shell::{ReShape, SharedReShape};
use crate::shhealing::{
    align_shell_closed_to_free_bounds, clear_shell_closed_if_analysis_free,
    fix_connected_at, fix_connected_closed, fix_reorder_wire_3d, fix_shell_face_orientation,
    fix_small_all, fix_add_curve3d_or_drop,
    fix_vertex_tolerance_edge, grown_vertex_tolerances_on_face, perform_wire_first_round,
    perform_wire_second_round, shape_fix_same_parameter, shell_open_by_free_edge,
    split_face_after_removed_segment, split_loop_wires, ShapeFixEdge, ShapeFixFace,
};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, vertices_of, wires_of_face};

use super::iges_data_iges_model::IgesModel;

/// `read.maxprecision.mode` default is 0 (`Interface_StaticStandards.cxx:45`).
/// `TrimTolerances` runs `LimitTolerance` only when the mode is 1
/// (`IGESToBRep_Actor.cxx:74-79`).
const MAX_PRECISION_MODE: i32 = 0;

/// `read.maxprecision.val` default (`Interface_StaticStandards.cxx:47`).
const MAX_PRECISION_VAL: f64 = 1.0;

/// `read.encoderegularity.angle` default (`Interface_StaticStandards.cxx:51`).
const ENCODE_REGULARITY_ANGLE: f64 = 0.01;

pub(in crate::iges_read) fn process_iges_shape(shape: &TopoShape, model: &IgesModel) -> TopoShape {
    let theeps = used_tolerance(model);
    let max_tol = (super::iges_to_brep_basic_curve::eps_geom(model.resolution) * model.unit_factor)
        .max(MAX_PRECISION_VAL);
    let mut ctx = SharedReShape::new();
    process_node(shape, &mut ctx, theeps, max_tol);
    // `ShapeFix_Shape.cxx:257`.
    let mut result = ctx.apply(shape);
    // `ShapeFix_Shape.cxx:259-261`. `NeedFix(myFixSameParameterMode)` is true for
    // the IGES default `-1` (`IGESData.cxx:116`). `ShapeFix_Shape::SameParameter`
    // (`ShapeFix_Shape.cxx:307-311`) is `ShapeFix::SameParameter(sh, enforce,
    // 0.0)` - `enforce` is false and `preci` is 0.
    shape_fix_same_parameter(&result, false, 0.0);
    // `cxx:267-291`. `NeedFix(myFixVertexTolMode)` is true (`IGESData.cxx:164`).
    fix_vertex_tolerances_after_same_parameter(&result, &mut ctx);
    // `cxx:294`.
    result = ctx.apply(&result);
    // `Actor.cxx:221-228`. A compound of one kind is not `TopAbs_SHAPE`.
    if shape_kind(&result) != ShapeType::Shape {
        encode_regul(&result);
        trim_tolerances(&result, theeps);
    }
    result
}

/// `ShapeFix_Shape.cxx:267-291`. Runs only when the result has more than one
/// face. `FixEdgeTool()` already carries the pass context, because each face
/// `Perform` forwarded it (`ShapeFix_Face.cxx:348`, `ShapeFix_Wire.cxx:307`),
/// so the update is `CopyVertex` (`ShapeFix_Edge.cxx:674-677`).
fn fix_vertex_tolerances_after_same_parameter(shape: &TopoShape, ctx: &mut SharedReShape) {
    let faces = faces_of(shape);
    if faces.len() <= 1 {
        return;
    }
    for face in &faces {
        for wire in wires_of_face(face) {
            for edge in edges_of_wire(&wire) {
                fix_vertex_tolerance_with_context(&edge, face, ctx);
            }
        }
    }
}

/// `ShapeFix_Edge::FixVertexTolerance(edge, face)` with a live context
/// (`ShapeFix_Edge.cxx:642-684`).
fn fix_vertex_tolerance_with_context(edge: &Edge, face: &Face, ctx: &mut SharedReShape) {
    let applied = ctx.apply(&edge.0);
    if !applied.is_edge() {
        return;
    }
    let Some(grown) = grown_vertex_tolerances_on_face(&Edge(applied), face) else {
        return;
    };
    ctx.copy_vertex(&grown.v1.0, grown.toler1);
    ctx.copy_vertex(&grown.v2.0, grown.toler2);
}

/// `IGESToBRep_Actor.cxx:186-190` over the ctor default `theeps = 0.0001`
/// (`cxx:88`).
fn used_tolerance(model: &IgesModel) -> f64 {
    if model.resolution > 1.0e-8 {
        model.resolution * model.unit_factor
    } else {
        1.0e-4
    }
}

fn process_node(shape: &TopoShape, ctx: &mut SharedReShape, preci: f64, max_tol: f64) {
    // `ShapeFix_Shape.cxx:117`: `S = Context()->Apply(myShape)`.
    let applied = ctx.apply(shape);
    match applied.shape_type() {
        ShapeType::Compound | ShapeType::CompSolid => {
            // `cxx:132-149`. `FixSameParameter` is cleared for the recursion
            // and restored after it, so only the outermost call runs it.
            let children = applied.tshape.read().unwrap().children.clone();
            for child in &children {
                process_node(child, ctx, preci, max_tol);
            }
        }
        ShapeType::Face => {
            // `cxx:193-207`. `NeedFix(myFixFaceMode)` is true (`IGESData.cxx:122`).
            perform_face(&Face(applied), ctx, preci, max_tol);
        }
        ShapeType::Shell => {
            // `NeedFix(myFixShellMode)` is true (`IGESData.cxx:120`).
            perform_iges_shell(&applied, ctx, preci, max_tol);
        }
        ShapeType::Solid => {
            // `NeedFix(myFixSolidMode)` is true (`IGESData.cxx:117`).
            perform_iges_solid(&applied, ctx, preci, max_tol);
        }
        ShapeType::Wire => {
            // `NeedFix(myFixWireMode)` is true (`IGESData.cxx:124`).
            perform_free_wire(&applied, ctx, preci, max_tol);
        }
        ShapeType::Edge => {
            // `cxx:238-246`. Context is set on this arm, and
            // `FixVertexTolerance(edge)` grows endpoint tolerances from every
            // pcurve (`checkAll=true`).
            fix_vertex_tolerance_edge(&crate::shape::Edge(applied));
        }
        _ => {}
    }
}

/// `ShapeFix_Shell::Perform` (`ShapeFix_Shell.cxx:101-173`).
/// `NeedFix(myFixFaceMode)` and `NeedFix(myFixOrientationMode)` are both true
/// (tool defaults `-1`). `CreateSolids` is not this function.
fn perform_iges_shell(shell: &TopoShape, ctx: &mut SharedReShape, preci: f64, max_tol: f64) {
    // `cxx:112`.
    let shell = ctx.apply(shell);
    if shell.shape_type() != ShapeType::Shell {
        return;
    }
    // `cxx:120-130`: each iterator face goes through `ShapeFix_Face::Perform`.
    for face in shell_faces_composed(&shell) {
        perform_face(&face, ctx, preci, max_tol);
    }
    // `cxx:139-142`, then the closed-flag write inside `FixFaceOrientation`
    // (`cxx:1474-1479`) before `GetShells`.
    let newsh = ctx.apply(&shell);
    if newsh.shape_type() == ShapeType::Shell {
        align_shell_closed_to_free_bounds(&newsh);
        let fixed = fix_shell_face_orientation(newsh.clone());
        if !fixed.same_tshape(&newsh) {
            ctx.replace(&newsh, &fixed);
        }
    }
    // `cxx:145-163`.
    let after = ctx.apply(&shell);
    if after.shape_type() == ShapeType::Shell {
        clear_shell_closed_if_analysis_free(&after);
    }
}

/// `ShapeFix_Solid::Perform` (`ShapeFix_Solid.cxx:460-643`).
/// `CreateOpenSolidMode` is 1 (`IGESData.cxx:119`). Several shells take
/// `CreateSolids` (`cxx:577-640`) when the shell count is not 1.
fn perform_iges_solid(solid: &TopoShape, ctx: &mut SharedReShape, preci: f64, max_tol: f64) {
    // `cxx:483-500`. `NeedFix(myFixShellMode)` is true (`IGESData.cxx:120`).
    let shells = shells_in(solid);
    for shell in &shells {
        perform_iges_shell(shell, ctx, preci, max_tol);
    }
    // `cxx:514`. `NeedFix(myFixShellOrientationMode)` is true (`IGESData.cxx:118`).
    let tmp = ctx.apply(solid);
    let shells = shells_in(&tmp);
    // The shortened `FixFaceOrientation` does not split a shell, so one input
    // shell is the `NbShells == 1` arm (`cxx:520-575`).
    if shells.len() != 1 {
        super::create_solids::replace_multi_shell(&tmp, ctx);
        return;
    }
    solid_from_one_shell(&tmp, &shells[0], ctx);
}

/// `ShapeFix_Wire::Perform` on a free wire (`ShapeFix_Shape.cxx:212-237`).
/// `SetFace` is null, so `IsReady` is false (`ShapeFix_Wire.lxx:68-70`).
/// Face fixes return immediately. `ModifyTopologyMode` is forced true, and
/// `ClosedWireMode` follows `S.Closed()`.
fn perform_free_wire(shape: &TopoShape, ctx: &mut SharedReShape, preci: f64, max_tol: f64) {
    let mut wire = Wire(shape.clone());
    let closed = shape.closed();
    let (reorder_ok, status) = fix_reorder_wire_3d(&wire);
    let mut fixed = status != crate::meshing::wire_order::WireOrderStatus::Same;
    let bare = Face(TopoShape::new(ShapeType::Face));
    fixed |= fix_small_all(&mut wire, &bare, preci, occt_core::precision::CONFUSION, !reorder_ok);
    if reorder_ok {
        fixed |= fix_connected_closed(&mut wire, preci, max_tol, closed);
    }
    fixed |= add_curve3d_or_drop(&mut wire, &bare, preci);
    copy_vertex_tolerances(&wire, ctx);
    if closed {
        let _ = fix_reorder_wire_3d(&wire);
    }
    copy_vertex_tolerances(&wire, ctx);
    if fixed {
        ctx.replace(shape, &wire.0);
    }
}

/// `FixAddCurve3d` (`ShapeFix_Wire.cxx:881-918`), shared with
/// `ShapeFix_Wire::perform_wire_first_round`'s `FixEdgeCurves` arm. The free
/// wire has no face, so the bare face is passed and the face-dependent
/// sub-calls of `FixClosed` (`cxx:909`) return false.
fn add_curve3d_or_drop(wire: &mut Wire, face: &Face, precision: f64) -> bool {
    fix_add_curve3d_or_drop(wire, face, precision)
}

/// `FixVertexTolerance(edge, null face)` with a context (`ShapeFix_Edge.cxx:674-677`).
fn copy_vertex_tolerances(wire: &Wire, ctx: &mut SharedReShape) {
    let bare = Face(TopoShape::new(ShapeType::Face));
    for edge in edges_of_wire(wire) {
        let Some(grown) = grown_vertex_tolerances_on_face(&edge, &bare) else {
            continue;
        };
        ctx.copy_vertex(&grown.v1.0, grown.toler1);
        ctx.copy_vertex(&grown.v2.0, grown.toler2);
    }
}

/// `ShapeFix_Solid::Perform` when `NbShells == 1` (`cxx:520-565`) and
/// `SolidFromShell` (`cxx:655-702`). `CreateOpenSolidMode` is 1, so an open
/// shell still enters `SolidFromShell`. `Context()->Replace` runs only when
/// the infinite point classifies `TopAbs_IN` (`ShapeExtend_DONE2`).
fn solid_from_one_shell(solid: &TopoShape, shell: &TopoShape, ctx: &mut SharedReShape) {
    let is_closed = !shell_open_by_free_edge(shell);
    shell.set_closed(is_closed);
    let mut sh = shell.clone();
    if !sh.free() {
        sh.set_free(true);
    }
    let builder = TopoBuilder::new();
    let trial = builder.make_solid(&[Shell(sh)]);
    let mut sc = SolidClassifier::new();
    sc.load(trial.0);
    sc.perform_infinite_point(occt_core::precision::CONFUSION);
    if sc.state() != FaceState::In {
        return;
    }
    let mut rev = shell.clone();
    if !rev.free() {
        rev.set_free(true);
    }
    rev.reverse();
    let soli2 = builder.make_solid(&[Shell(rev)]);
    ctx.replace(solid, &soli2.0);
}

/// `TopoDS_Iterator` faces with orientation composed onto the shell
/// (`ShapeFix_Shell.cxx:120`).
fn shell_faces_composed(shell: &TopoShape) -> Vec<Face> {
    let parent = shell.orientation();
    let stored = shell
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    stored
        .into_iter()
        .filter(|child| child.shape_type() == ShapeType::Face)
        .map(|mut child| {
            child.set_orientation(Orientation::compose(parent, child.orientation()));
            Face(child)
        })
        .collect()
}

/// `TopExp_Explorer(shape, TopAbs_SHELL)` (`ShapeFix_Solid.cxx:475`).
fn shells_in(shape: &TopoShape) -> Vec<TopoShape> {
    let mut out = Vec::new();
    collect_shells(shape, &mut out);
    out
}

fn collect_shells(shape: &TopoShape, out: &mut Vec<TopoShape>) {
    if shape.shape_type() == ShapeType::Shell {
        out.push(shape.clone());
        return;
    }
    let children = shape.tshape.read().expect("poisoned TShape lock").children.clone();
    for child in children {
        collect_shells(&child, out);
    }
}

/// Ported tail of `ShapeFix_Face::Perform` (`ShapeFix_Face.cxx:482-716`),
/// the same order `read_topology.rs` uses after the face is built.
fn perform_face(face: &Face, ctx: &mut SharedReShape, preci: f64, max_tol: f64) {
    let applied = ctx.apply(&face.0);
    let face = if applied.is_face() {
        Face(applied)
    } else {
        face.clone()
    };
    // `ShapeFix_Face.cxx:365-480`: wire `Perform` before the seam steps,
    // with lacking and self-intersection forced off.
    for mut wire in wires_stored_on_face(&face) {
        perform_wire_first_round(&mut wire, &face, preci, ctx);
    }
    let mut sff = ShapeFixFace::with_face_and_context(&face, ctx.clone());
    sff.max_tol = max_tol;
    sff.result = Some(face.0.clone());
    // `cxx:486-498`.
    sff.fix_periodic_degenerated();
    let seam_ret = sff.fix_missing_seam();
    if seam_ret {
        if let Some(res) = sff.result.clone() {
            let sub_faces = sub_faces_of(&res);
            if !sub_faces.is_empty() {
                let mut removed_segment = false;
                for rf in &sub_faces {
                    for mut wire in wires_stored_on_face(rf) {
                        if manifold_edge_count(&wire) == 0 {
                            continue;
                        }
                        if perform_wire_second_round(&mut wire, rf, preci) {
                            removed_segment = true;
                        }
                    }
                    drop_empty_manifold_wires(rf);
                }
                for rf in &sub_faces {
                    split_loop_wires(rf);
                }
                if removed_segment {
                    for rf in &sub_faces {
                        split_face_after_removed_segment(rf, ctx);
                    }
                }
                for rf in &sub_faces {
                    fix_face_second_part(&mut sff, rf);
                }
                return;
            }
        }
    }
    // `cxx:508-627` on the unsplit face, then `cxx:629-674` when a segment
    // was removed, then `cxx:676-716`.
    let mut removed_segment = false;
    for mut wire in wires_stored_on_face(&face) {
        if manifold_edge_count(&wire) == 0 {
            continue;
        }
        if perform_wire_second_round(&mut wire, &face, preci) {
            removed_segment = true;
        }
    }
    drop_empty_manifold_wires(&face);
    split_loop_wires(&face);
    if removed_segment {
        split_face_after_removed_segment(&face, ctx);
    }
    fix_face_second_part(&mut sff, &face);
}

fn sub_faces_of(shape: &TopoShape) -> Vec<Face> {
    if shape.shape_type() == ShapeType::Face {
        vec![Face(shape.clone())]
    } else if shape.shape_type() == ShapeType::Shell {
        faces_of(shape)
    } else {
        Vec::new()
    }
}

/// `TopoDS_Iterator(face, false)` (`ShapeFix_Face.cxx:400`).
fn wires_stored_on_face(face: &Face) -> Vec<Wire> {
    face.0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .iter()
        .filter(|child| child.shape_type() == ShapeType::Wire)
        .cloned()
        .map(Wire)
        .collect()
}

/// `ShapeExtend_WireData::NbEdges` (`ShapeFix_Face.cxx:541`): forward and
/// reversed edges only.
fn manifold_edge_count(wire: &Wire) -> i32 {
    edges_of_wire(wire)
        .iter()
        .filter(|edge| {
            matches!(
                edge.0.orientation(),
                Orientation::Forward | Orientation::Reversed
            )
        })
        .count() as i32
}

fn nonmanifold_edge_count(wire: &Wire) -> i32 {
    edges_of_wire(wire)
        .iter()
        .filter(|edge| {
            matches!(
                edge.0.orientation(),
                Orientation::Internal | Orientation::External
            )
        })
        .count() as i32
}

/// A wire with no edges is not added back (`cxx:547-551`). A wire that has
/// only non-manifold edges stays (`cxx:543-545`).
fn drop_empty_manifold_wires(face: &Face) {
    let children = face.0.tshape.read().expect("lock").children.clone();
    let kept: Vec<TopoShape> = children
        .into_iter()
        .filter(|child| {
            if child.shape_type() != ShapeType::Wire {
                return true;
            }
            let wire = Wire(child.clone());
            manifold_edge_count(&wire) > 0 || nonmanifold_edge_count(&wire) > 0
        })
        .collect();
    let mut ts = face.0.tshape.write().expect("lock");
    if kept.len() != ts.children.len() {
        ts.children = kept;
    }
}

/// `ShapeFix_Face::Perform` second-part loop (`ShapeFix_Face.cxx:676-716`).
fn fix_face_second_part(sff: &mut ShapeFixFace, face: &Face) -> TopoShape {
    let face = fix_face_wire_pairs(sff, face);
    let face = fix_face_orientation(sff, &face);
    // `BRepTools::Update(myFace)` (`cxx:700`) only drops the cached
    // triangulation. None is stored at transfer time.
    let mut need_split = true;
    sff.face = Some(face.clone());
    if sff.fix_add_natural_bound() {
        need_split = false;
    }
    let face = sff.face.clone().unwrap_or(face);
    if sff.fix_split_face_mode && need_split && sff.map_wires.len() > 1 {
        sff.face = Some(face.clone());
        if sff.fix_split_face() {
            if let Some(res) = sff.result.clone() {
                return res;
            }
        }
    }
    face.0
}

fn fix_face_wire_pairs(sff: &mut ShapeFixFace, face: &Face) -> Face {
    sff.face = Some(face.clone());
    sff.fix_wires_two_coinc_edges();
    if sff.fix_intersecting_wires_mode {
        sff.fix_intersecting_wires();
    }
    sff.face.clone().unwrap_or_else(|| face.clone())
}

fn fix_face_orientation(sff: &mut ShapeFixFace, face: &Face) -> Face {
    if !sff.fix_orientation_mode {
        return face.clone();
    }
    sff.face = Some(face.clone());
    sff.fix_orientation();
    sff.face.clone().unwrap_or_else(|| face.clone())
}

/// `ShapeExtend_Explorer::ShapeType(shape, true)` reduced to the compound
/// case `Actor.cxx:221` cares about: one child kind, or `Shape` when mixed.
fn shape_kind(shape: &TopoShape) -> ShapeType {
    match shape.shape_type() {
        ShapeType::Compound | ShapeType::CompSolid => {
            let children = shape.tshape.read().unwrap().children.clone();
            let mut kind = ShapeType::Shape;
            for child in &children {
                let next = shape_kind(child);
                if next == ShapeType::Shape {
                    return ShapeType::Shape;
                }
                if kind == ShapeType::Shape {
                    kind = next;
                } else if kind != next {
                    return ShapeType::Shape;
                }
            }
            kind
        }
        other => other,
    }
}

/// `EncodeRegul` (`IGESToBRep_Actor.cxx:44-66`).
fn encode_regul(shape: &TopoShape) {
    if ENCODE_REGULARITY_ANGLE <= 0.0 {
        return;
    }
    crate::brep_lib_continuity::encode_regularity(shape, ENCODE_REGULARITY_ANGLE);
}

/// `TrimTolerances` (`IGESToBRep_Actor.cxx:72-80`).
fn trim_tolerances(shape: &TopoShape, the_tolerance: f64) {
    if MAX_PRECISION_MODE == 1 {
        let tmax = the_tolerance.max(MAX_PRECISION_VAL);
        limit_tolerance(shape, 0.0, tmax);
    }
}

/// `ShapeFix_ShapeTolerance::LimitTolerance` (`ShapeFix_ShapeTolerance.cxx:35-137`)
/// for `styp == TopAbs_SHAPE`: vertices, edges, and faces.
fn limit_tolerance(shape: &TopoShape, tmin: f64, tmax: f64) {
    if tmin < 0.0 {
        return;
    }
    let iamax = tmax >= tmin;
    let reg = GeometryRegistry::global();
    for vertex in vertices_of(shape) {
        clamp_vertex(&reg, &vertex.0, tmin, tmax, iamax);
    }
    for edge in edges_of(shape) {
        let prec = crate::brep_tool::BRepTool::edge_tolerance(&edge);
        if let Some(tol) = clamped(prec, tmin, tmax, iamax) {
            reg.set_edge_tolerance(&edge.0, tol);
        }
    }
    for face in faces_of(shape) {
        let prec = crate::brep_tool::BRepTool::face_tolerance(&face);
        if let Some(tol) = clamped(prec, tmin, tmax, iamax) {
            reg.set_face_tolerance(&face.0, tol);
        }
    }
}

fn clamp_vertex(reg: &GeometryRegistry, shape: &TopoShape, tmin: f64, tmax: f64, iamax: bool) {
    let Some(mut geom) = reg.vertex_geom(shape) else {
        return;
    };
    if let Some(tol) = clamped(geom.tolerance, tmin, tmax, iamax) {
        geom.tolerance = tol;
        reg.set_vertex(shape, geom);
    }
}

fn clamped(prec: f64, tmin: f64, tmax: f64, iamax: bool) -> Option<f64> {
    if iamax && prec > tmax {
        Some(tmax)
    } else if prec < tmin {
        Some(tmin)
    } else {
        None
    }
}
