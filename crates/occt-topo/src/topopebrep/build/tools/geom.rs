//! `TopOpeBRepBuild_Tools` — the geometric queries.
//!
//! Sources: `TopOpeBRepBuild_Tools.cxx:461-484` (`GetNormalToFaceOnEdge`),
//! `:486-541` (`GetNormalInNearestPoint`), `:543-607` (`GetTangentToEdgeEdge`),
//! `:609-627` (`GetTangentToEdge`), `:629-660` (`GetAdjacentFace`),
//! `:838-853` (`IsDegEdgesTheSame`) — OCCT 8.0.0.

use std::collections::HashMap;

use occt_core::gp::{GpPnt2d, GpVec, GpVec2d};
use occt_geom2d::curve::Curve2d;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{is_equal, is_same, shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};
use crate::topopebrep::tool::curvature::min_duv;

/// `const double PAR_T = 0.43213918` (`cxx:57`): `10 * e^(-PI)`.
pub const PAR_T: f64 = 0.43213918;

/// `FC2D_CurveOnSurface(E, F, f2, l2, tolpc, true)`: the p-curve of `E` on `F`
/// plus its parametric range.
pub(crate) fn fc2d_curve_on_surface(
    e: &Edge,
    f: &Face,
) -> Option<(std::sync::Arc<dyn Curve2d>, f64, f64)> {
    crate::boptools_2d::build_pcurve_for_edge_on_face(e, f).ok()?;
    crate::boptools_2d::curve_on_surface_range(e, f)
}

/// `GetNormalToFaceOnEdge(aFObj, anEdgeObj, aDirNormal)` (`cxx:461-484`):
/// the unnormalized surface normal `D1U ^ D1V` at the p-curve parameter
/// `f * PAR_T + (1 - PAR_T) * l`.
pub fn get_normal_to_face_on_edge(f_obj: &Face, e_obj: &Edge) -> Option<GpVec> {
    let (c2d, _f2, _l2) = fc2d_curve_on_surface(e_obj, f_obj)?;
    let a_ca = CurveAdaptor::of_edge(e_obj);
    let f = a_ca.first_parameter();
    let l = a_ca.last_parameter();
    let par = f * PAR_T + (1.0 - PAR_T) * l;
    let a_uv1 = c2d.d0(par);

    let a_sa = SurfaceAdaptor::of_face(f_obj);
    let (_p, a_tg1, a_tg2) = a_sa.d1(a_uv1.x(), a_uv1.y())?;
    Some(a_tg1.crossed(&a_tg2))
}

/// `GetNormalInNearestPoint(F, E, aNormal)` (`cxx:486-541`): the normal a
/// small step away from the p-curve, used where the normal *on* the curve
/// vanishes; the step direction is flipped when it leaves the face.
pub fn get_normal_in_nearest_point(f: &Face, e: &Edge) -> Option<GpVec> {
    let (c2d, f2, l2) = fc2d_curve_on_surface(e, f)?;
    let par = f2 * PAR_T + (1.0 - PAR_T) * l2;
    let (a_p, a_tangent) = c2d.d1(par);

    let x_norm = -a_tangent.y();
    let y_norm = a_tangent.x();

    let mut step = min_duv(f);
    step *= 1e-2;

    let a_pv = GpVec2d::new(a_p.x(), a_p.y());
    let mut a_step_v = occt_core::gp::GpDir2d::new(x_norm, y_norm).ok()?;
    let mut a_norm2d = a_pv.added(&GpVec2d::from_dir2d(&a_step_v).multiplied_scalar(step));

    let mut new_u = a_norm2d.x();
    let mut new_v = a_norm2d.y();

    let bs = SurfaceAdaptor::of_face(f);
    let (_p1, mut a_tg1, mut a_tg2) = bs.d1(new_u, new_v)?;

    // `BRepTopAdaptor_FClass2d FC(F, Precision::PConfusion())` (`cxx:521-523`).
    let a_p2d = GpPnt2d::new(new_u, new_v);
    let fc = crate::fclass2d::FClass2d::new(f, occt_core::precision::PCONFUSION).ok()?;
    let a_state = fc.perform(a_p2d);

    // Point out of the face: go the other way (`cxx:526-539`).
    if a_state == FaceState::Out {
        a_step_v.reverse();
        a_norm2d = a_pv.added(&GpVec2d::from_dir2d(&a_step_v).multiplied_scalar(step));
        new_u = a_norm2d.x();
        new_v = a_norm2d.y();
        let (p1, tg1, tg2) = bs.d1(new_u, new_v)?;
        let _ = p1;
        a_tg1 = tg1;
        a_tg2 = tg2;
    }

    Some(a_tg1.crossed(&a_tg2))
}

/// `GetTangentToEdgeEdge(aFObj, anEdgeObj, aOriEObj, aTangent)`
/// (`cxx:543-607`): the tangent of the piece edge, oriented like the original
/// edge. The returned `bool` is OCCT's return value: `true` when the tangent
/// had to be reversed because the two curves run opposite.
pub fn get_tangent_to_edge_edge(
    an_edge_obj: &Edge,
    a_ori_e_obj: &Edge,
) -> Option<(GpVec, bool)> {
    if BRepTool::is_degenerated(a_ori_e_obj) || BRepTool::is_degenerated(an_edge_obj) {
        return Some((get_tangent_to_edge(an_edge_obj)?, true));
    }

    let a_ca = CurveAdaptor::of_edge(an_edge_obj);
    let a_ca_ori = CurveAdaptor::of_edge(a_ori_e_obj);

    let f = a_ca.first_parameter();
    let l = a_ca.last_parameter();
    let par = f * PAR_T + (1.0 - PAR_T) * l;

    let (_p, mut a_tangent) = a_ca.d1(par)?;

    // `GeomAPI_ProjectPointOnCurve(aP, aCopyCurve, First, Last)` (`cxx:585-592`):
    // the adaptor carries the edge location and the same parameter range as the
    // transformed copy OCCT projects onto.
    let p = a_ca.d0(par)?;
    let proj = occt_geom::geom_api::project_point_on_curve(
        a_ca_ori.curve()?.as_ref(),
        &p,
        occt_core::precision::CONFUSION,
    )?;
    let par_ori = proj.parameter;

    let (_p_ori, a_tg_ori) = a_ca_ori.d1(par_ori)?;

    if an_edge_obj.0.orientation() == crate::abs::Orientation::Reversed {
        a_tangent.reverse();
    }

    if a_tg_ori.dot(&a_tangent) < 0.0 {
        a_tangent.reverse();
        return Some((a_tangent, true));
    }
    Some((a_tangent, false))
}

/// `GetTangentToEdge(anEdgeObj, aTangent)` (`cxx:609-627`): always succeeds in
/// OCCT (it returns `true` unconditionally).
pub fn get_tangent_to_edge(an_edge_obj: &Edge) -> Option<GpVec> {
    let a_ca = CurveAdaptor::of_edge(an_edge_obj);
    let f = a_ca.first_parameter();
    let l = a_ca.last_parameter();
    let par = f * PAR_T + (1.0 - PAR_T) * l;
    let (_p, a_tangent) = a_ca.d1(par)?;
    Some(a_tangent)
}

/// `GetAdjacentFace(aFaceObj, anEObj, anEdgeFaceMap, anAdjFaceObj)`
/// (`cxx:629-660`): the first face of `anEObj`'s ancestor list that is not
/// `aFaceObj` (`operator!=` is `!IsEqual`, orientation included).
pub fn get_adjacent_face(
    a_face_obj: &TopoShape,
    an_e_obj: &TopoShape,
    an_edge_face_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> Option<TopoShape> {
    let key = shape_map_key(an_e_obj);
    let list = an_edge_face_map.get(&key)?;
    list.iter()
        .find(|f| !is_equal(f, a_face_obj))
        .cloned()
}

/// `IsDegEdgesTheSame(anE1, anE2)` (`cxx:838-853`): both edges have a vertex
/// and their first vertices are `IsSame`.
pub fn is_deg_edges_the_same(an_e1: &TopoShape, an_e2: &TopoShape) -> bool {
    let map1 = super::state::map_shapes_of(an_e1, ShapeType::Vertex);
    let map2 = super::state::map_shapes_of(an_e2, ShapeType::Vertex);

    if map1.extent() == 0 || map2.extent() == 0 {
        return false;
    }
    match (map1.find_key(1), map2.find_key(1)) {
        (Some(a), Some(b)) => is_same(a, b),
        _ => false,
    }
}
