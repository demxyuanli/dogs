//! `TopOpeBRepTool_TOOL::Matter` — the oriented angle between two directions.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRepTool`):
//! - `Matter(d1, dR2, Ref)` (`TopOpeBRepTool_TOOL.cxx:1338-1361`)
//! - `Matter(d1, dR2)` (`TopOpeBRepTool_TOOL.cxx:1365-1373`), the 2D arm, which
//!   lifts both vectors into the plane `Z = 0` with `Ref = (0, 0, 1)`
//!
//! - `Matter(xx1, nt1, xx2, nt2, tola, ang)` (`cxx:1377-1427`)
//! - `Matter(f1, f2, e, par, tola, ang)` (`cxx:1661-1708`)
//! - `MatterKPtg(f1, f2, e, ang)` (`cxx:1710-1771`)

use occt_core::gp::{GpDir, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::ANGULAR;

use crate::shape::TopoShape;
use crate::topopebrep::tool::curvature;
use crate::topopebrep::tool::project;
use crate::topopebrep::tool::project_curve;
use crate::topopebrep::tool::topology;
use crate::topopebrep::tool::uv_app;

/// `TopOpeBRepTool_TOOL::Matter(d1, dR2, Ref)` (`cxx:1338-1361`): the angle
/// from `d1` to `-dR2`, measured around `Ref` in `[0, 2*Pi)`.
///
/// `0` and `Pi` are special-cased on `Precision::Angular()`; every other pair
/// goes through `gp_Vec::AngleWithRef`, whose negative answer is shifted by
/// `2*Pi`.
pub fn matter(d1: &GpVec, d_r2: &GpVec, v_ref: &GpVec) -> f64 {
    let d2 = d_r2.reversed();
    let tola = ANGULAR;
    let ang = d1.angle(&d2);
    if ang < tola {
        return 0.0;
    }
    if (std::f64::consts::PI - ang) < tola {
        return std::f64::consts::PI;
    }
    let ang = d1.angle_with_ref(&d2, v_ref);
    if ang < 0.0 {
        2.0 * std::f64::consts::PI + ang
    } else {
        ang
    }
}

/// `TopOpeBRepTool_TOOL::Matter(const gp_Vec2d&, const gp_Vec2d&)`
/// (`cxx:1365-1373`): the 3D arm with `Z = 0` and `Ref = (0, 0, 1)`.
pub fn matter_2d(d1: &GpVec2d, d_r2: &GpVec2d) -> f64 {
    let v1 = GpVec::new(d1.x(), d1.y(), 0.0);
    let vr2 = GpVec::new(d_r2.x(), d_r2.y(), 0.0);
    let v_ref = GpVec::new(0.0, 0.0, 1.0);
    matter(&v1, &vr2, &v_ref)
}

/// `TopOpeBRepTool_TOOL::Matter(xx1, nt1, xx2, nt2, tola, ang)`
/// (`cxx:1377-1427`): the matter angle of `f2` seen from `f1` around the
/// tangent `z = xx ^ nt`.
///
/// `None` is the C++ `false`: the two `z` axes are not opposite (`cxx:1405-1408`).
/// `cxx:1418-1420` pins a near-zero answer to `0.`; a negative
/// `AngleWithRef` is shifted by `2*Pi`.
pub fn matter_dirs(xx1: &GpDir, nt1: &GpDir, xx2: &GpDir, nt2: &GpDir, tola: f64) -> Option<f64> {
    let z1 = xx1.crossed(nt1).ok()?;
    let z2 = xx2.crossed(nt2).ok()?;
    if z2.dot(&z1) >= 0.0 {
        return None;
    }
    // `cxx:1414`: `-nti` points towards the 3d matter of `fi`, so `zi` is
    // reversed before the angle is measured.
    let z1 = z1.reversed();
    let mut ang = xx1.angle_with_ref(xx2, &z1);
    if ang.abs() < tola {
        return Some(0.0);
    }
    if ang < 0.0 {
        ang = 2.0 * std::f64::consts::PI + ang;
    }
    Some(ang)
}

/// `TopOpeBRepTool_TOOL::Matter(f1, f2, e, par, tola, ang)` (`cxx:1661-1708`):
/// the normal and in-face direction of both faces at `par` on `e`, then
/// [`matter_dirs`].
///
/// The face tolerances are scaled by `1.e2` (`f1`) and `2.e2` (`f2`), the two
/// `// nyitolxpu` provisional factors of OCCT.
pub fn matter_faces(
    f1: &TopoShape,
    f2: &TopoShape,
    e: &TopoShape,
    par: f64,
    tola: f64,
) -> Option<f64> {
    let face1 = crate::shape::Face(f1.clone());
    let face2 = crate::shape::Face(f2.clone());
    let tol1 = crate::brep_tool::BRepTool::face_tolerance(&face1) * 1.0e2;
    let uv1 = topology::fun_tool_par_on_e_f_tol(e, par, f1, tol1)?;
    let nt1 = curvature::nt(&uv1, &face1)?;
    let xx1 = curvature::xx(&uv1, f1, par, e)?;

    let tol2 = crate::brep_tool::BRepTool::face_tolerance(&face2) * 2.0e2;
    let uv2 = topology::fun_tool_par_on_e_f_tol(e, par, f2, tol2)?;
    let nt2 = curvature::nt(&uv2, &face2)?;
    let xx2 = curvature::xx(&uv2, f2, par, e)?;

    matter_dirs(&xx1, &nt1, &xx2, &nt2, tola)
}

/// `TopOpeBRepTool_TOOL::MatterKPtg(f1, f2, e, ang)` (`cxx:1710-1771`): the
/// coarse fallback used when the matter angle of the regular computation is
/// below its tolerance.
///
/// `cxx:1724-1726` takes `par = (1 - x) * first + x * last` with `x = 0.45678`
/// and the `0.123` ``UV`` step of `uvApp`; the only two answers are `0.` and
/// `2*Pi`, from the sign of `(proj(p_uvapp, f2) - p_uvapp) . Nt(f1)`.
///
/// `cxx:1728` ignores the `FUN_tool_paronEF` result, so a failed projection
/// leaves `uv1` at the default `(0, 0)` point.
pub fn matter_k_ptg(f1: &TopoShape, f2: &TopoShape, e: &TopoShape) -> Option<f64> {
    let (first, last) = project_curve::fun_tool_bounds(e);
    let x = 0.45678;
    let par = (1.0 - x) * first + x * last;
    let eps = 0.123; // NYIXPU190199
    let face1 = crate::shape::Face(f1.clone());
    let face2 = crate::shape::Face(f2.clone());

    let uv1 = topology::fun_tool_par_on_e_f(e, par, f1).unwrap_or(GpPnt2d::new(0.0, 0.0));
    let nt1 = curvature::nt(&uv1, &face1)?;
    let uv_app1 = uv_app::uv_app(f1, e, par, eps)?;
    let pf1 = topology::fun_tool_value_uv(&uv_app1, &face1)?;

    let (u2, v2, _d) = project::fun_tool_proj_pon_f(&pf1, &face2)?;
    let uv2 = GpPnt2d::new(u2, v2);
    let pf2 = topology::fun_tool_value_uv(&uv2, &face2)?;

    let v12 = GpDir::from_vec(&GpVec::new(
        pf2.x() - pf1.x(),
        pf2.y() - pf1.y(),
        pf2.z() - pf1.z(),
    ))
    .ok()?;
    let dot = v12.dot(&nt1);
    Some(if dot < 0.0 { 0.0 } else { 2.0 * std::f64::consts::PI })
}
