//! `ShapeFix::SameParameter` (`ShapeFix.cxx:70-275`).
//!
//! `ShapeFix_Shape::SameParameter` (`ShapeFix_Shape.cxx:307-311`) forwards to
//! this function with `preci = 0.0` and a null message registrator, and
//! `ShapeFix_Shape::Perform` calls it at `ShapeFix_Shape.cxx:259-261` whenever
//! `FixSameParameterMode` is not `NotFix`. Both read paths configure that key
//! as `-1` (`ShapeProcess_OperLibrary.cxx:813` from
//! `IGESData.cxx:116` / `STEPControl_Controller.cxx:211` and
//! `DE_ShapeFixParameters.hxx:39` `FixMode::FixOrNot`), and `NeedFix(-1)` is
//! true (`ShapeFix_Root.lxx:101-104`), so both operators run this function at
//! the tail of `FixShape`. It is *not* `BRepLib::SameParameter`
//! (`BRepLib.cxx:1016-1140`, `:1251-1739`): the per-edge work is
//! `ShapeFix_Edge::FixSameParameter` over the ancestor faces of
//! `TopExp::MapShapesAndAncestors`, plus a planar tolerance stage
//! (`ShapeFix.cxx:182-258`) that `BRepLib` has no equivalent of.

use super::prelude::*;

use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::iterator::cumulated_children;
use crate::shape_fix_compose_shell::{ReShape, SharedReShape};
use crate::topo_tools_full::{
    edges_of_wire, explorer_shapes, faces_of, map_shapes_and_ancestors, wires_of_face,
};

use super::{fix_same_parameter, fix_same_parameter_no_face, grown_vertex_tolerances_on_face};

/// `ShapeFix_Shape::Perform`'s FixVertexTolerance stage
/// (`ShapeFix_Shape.cxx:267-291`). Runs only when the result has more than one
/// face. `FixEdgeTool()` already carries the pass context, because each face
/// `Perform` forwarded it (`ShapeFix_Face.cxx:348`, `ShapeFix_Wire.cxx:307`),
/// so the update is `CopyVertex` (`ShapeFix_Edge.cxx:674-677`).
pub fn fix_vertex_tolerances_after_same_parameter(shape: &TopoShape, ctx: &mut SharedReShape) {
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

/// `ShapeFix::SameParameter(shape, enforce, preci, theProgress, theMsgReg)`
/// (`ShapeFix.cxx:70-275`).
///
/// `preci` reaches the function through `iatol`/`tol` only (`cxx:99-100`,
/// `cxx:125-128`) and `tol` is never read again, so it cannot change the
/// result; it is kept for the control-flow shape. `theMsgReg` is null on this
/// call path (`ShapeFix_Shape.cxx:250-255` passes the four-argument overload,
/// whose `theMsgReg` defaults to null), so the `else if` at `cxx:166-168` never
/// sends a message. Returns the `status` flag, which is false when any edge is
/// still not same-parameter after the fix.
pub fn shape_fix_same_parameter(shape: &TopoShape, enforce: bool, preci: f64) -> bool {
    let reg = GeometryRegistry::global();
    let iatol = preci > 0.0;
    // `cxx:90-91`: the ancestor map is built once, for the whole shape.
    let ancestors = map_shapes_and_ancestors(shape, ShapeType::Edge, ShapeType::Face);
    let mut status = true;
    // `cxx:102`, `cxx:111-112`: `TopExp_Explorer(shape, TopAbs_EDGE)` yields
    // occurrences, not distinct edges, so a shared edge is fixed once per
    // occurrence.
    for edge_shape in explorer_shapes(shape, ShapeType::Edge) {
        let edge = Edge(edge_shape);
        if !iatol {
            // `cxx:125-128`: `tol = BRep_Tool::Tolerance(E)`. `tol` is not read
            // on any later path, so the read is kept only as the control-flow
            // marker it is in OCCT.
            let _ = BRepTool::edge_tolerance(&edge);
        }
        if enforce {
            // `cxx:129-133`.
            reg.set_same_range(&edge.0, false);
            reg.set_same_parameter(&edge.0, false);
        }
        // `cxx:135-148`: `aMapEF.FindFromKey(E, aListOfFaces)`.
        let key = GeometryRegistry::shape_key(&edge.0);
        match ancestors.get(&key) {
            Some(faces) if !faces.is_empty() => {
                for face in faces {
                    fix_same_parameter(&edge, &Face(face.clone()));
                }
            }
            // `cxx:146-148` `K2-SEP97`: no ancestor face, empty-face overload.
            _ => fix_same_parameter_no_face(&edge),
        }
        // `cxx:151-166`.
        if !BRepTool::same_parameter(&edge) {
            status = false;
            reg.set_same_range(&edge.0, false);
            reg.set_same_parameter(&edge.0, false);
        }
    }
    // `cxx:182-266`: the planar tolerance stage.
    update_tolerances_on_planes(shape);
    status
}

/// `ShapeFix.cxx:182-266`: for every edge of every planar face, compare the 3D
/// curve against its projection on the plane at 23 uniform parameters and, when
/// the deviation exceeds the edge tolerance, raise the edge and both endpoint
/// vertex tolerances to `1.00005 * sqrt(max square distance)`.
fn update_tolerances_on_planes(shape: &TopoShape) {
    let reg = GeometryRegistry::global();
    // `cxx:186-188`: `aNbFaces` sizes the progress scope only.
    for face_shape in explorer_shapes(shape, ShapeType::Face) {
        let face = Face(face_shape);
        // `cxx:191`.
        let Some(surf) = BRepTool::face_surface(&face) else {
            continue;
        };
        // `cxx:192-206`: only `Geom_Plane` (directly, or as the basis of one
        // `Geom_RectangularTrimmedSurface`) is processed; a null `Geom_Plane`
        // continues the loop.
        let Some(plane) = geom_plane_of(&surf) else {
            continue;
        };
        // `cxx:209`: `TopExp_Explorer(face, TopAbs_EDGE)`.
        for edge_shape in explorer_shapes(&face.0, ShapeType::Edge) {
            let edge = Edge(edge_shape);
            // `cxx:213-218`: `BRep_Tool::Curve(edge, f, l)`, null -> skip.
            let Some(crv) = BRepTool::edge_curve(&edge) else {
                continue;
            };
            // `cxx:219-223`: `BRep_Tool::CurveOnSurface(edge, face, f, l)`,
            // null -> skip. The call overwrites `f`/`l` with the pcurve range,
            // which is what the sample parameters are built from.
            let Some((c2d, f, l)) = curve_on_surface(&edge, &face) else {
                continue;
            };
            // `cxx:224-225`: `Adaptor3d_CurveOnSurface(GHPC, AS)` with
            // `GEOM2DADAPTOR(GHPC, f, l)` and `AS = GeomAdaptor_Surface(plane)`.
            let acs = crate::meshing::edge_discret::CurveOnSurface::new(c2d, plane.clone(), f, l);
            let tol0 = BRepTool::edge_tolerance(&edge);
            let mut tol2 = tol0 * tol0;
            let mut is_changed = false;
            // `cxx:230-242`.
            const NCONTROL: i32 = 23;
            for i in 0..NCONTROL {
                let par = (f * f64::from(NCONTROL - 1 - i) + l * f64::from(i))
                    / f64::from(NCONTROL - 1);
                let pnt = crv.d0(par);
                let prj = acs.d0(par);
                let dist = pnt.square_distance(&prj);
                if tol2 < dist {
                    tol2 = dist;
                    is_changed = true;
                }
            }
            if !is_changed {
                continue;
            }
            // `cxx:244-255`: `tol = 1.00005 * sqrt(tol2)` (coefficient from
            // trj3_pm1-ct-203.stp #19681, edge 10) and the write happens only
            // when the new value does not shrink the tolerance.
            let tol = 1.00005 * tol2.sqrt();
            if tol >= tol0 {
                reg.set_edge_tolerance(&edge.0, tol);
                for vertex in cumulated_children(&edge.0) {
                    if vertex.shape_type() == ShapeType::Vertex {
                        Vertex(vertex).set_tolerance(tol);
                    }
                }
            }
        }
    }
}

/// `down_cast<Geom_Plane>(Surf)`, else `down_cast<Geom_Plane>(GRTS->BasisSurface())`
/// (`ShapeFix.cxx:193-205`). `Geom_RectangularTrimmedSurface::gp_pln` reports
/// its *basis* type for the adaptor (`GeomAdaptor_Surface.cxx:423-425`), so a
/// second trim level is rejected here the way a null `occ::down_cast` rejects
/// it in OCCT.
fn geom_plane_of(surf: &Arc<dyn Surface>) -> Option<Arc<dyn Surface>> {
    match surf.rectangular_trimmed_basis() {
        Some(basis) => {
            if basis.rectangular_trimmed_basis().is_none() && basis.gp_pln().is_some() {
                Some(basis)
            } else {
                None
            }
        }
        None => surf.gp_pln().is_some().then(|| surf.clone()),
    }
}

/// `BRep_Tool::CurveOnSurface(edge, face, First, Last)`
/// (`BRep_Tool.cxx:301-315` then `:327-372`).
///
/// Returns the edge's stored `CurveOnSurface` representation on the face's
/// surface - `PCurve2` when the representation is on a closed surface and the
/// face-composed local edge is `REVERSED` (`cxx:310-313`, `:354-357`) - with the
/// representation's range (`GC->Range`, `cxx:353`). With no representation on
/// that surface the call falls through to `CurveOnPlane` (`cxx:367-372`), which
/// returns a null handle unless the surface is a plane and the edge carries a
/// 3D curve, and reports the 3D range (`cxx:418-419`).
fn curve_on_surface(edge: &Edge, face: &Face) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    let pcurves = reg.edge_pcurves(&edge.0, face_key);
    if !pcurves.is_empty() {
        let mut reversed = edge.0.orientation().is_reversed();
        if face.0.orientation() == Orientation::Reversed {
            reversed = !reversed;
        }
        let pcurve = if pcurves.len() >= 2 && reversed {
            pcurves[1].clone()
        } else {
            pcurves[0].clone()
        };
        let (f, l) = reg
            .pcurve_range(&edge.0, face_key)
            .unwrap_or_else(|| reg.edge_parameters(&edge.0));
        return Some((pcurve, f, l));
    }
    let surf = reg.face_surface(&face.0)?;
    let curve = reg.edge_curve(&edge.0)?;
    let pcurve = crate::pcurve::curve_on_plane(curve.as_ref(), surf.as_ref())?;
    let (f, l) = reg.edge_parameters(&edge.0);
    Some((pcurve, f, l))
}
