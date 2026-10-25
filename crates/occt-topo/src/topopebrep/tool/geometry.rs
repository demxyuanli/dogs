//! Geometry predicates and normal construction from
//! `TopOpeBRepTool_GEOMETRY.cxx`.
//!
//! Sources: `TKBool/TopOpeBRepTool/TopOpeBRepTool_GEOMETRY.cxx:86-358`.
//!
//! `BASISCURVE2D` (`cxx:35-50`) is folded into the callers: it unwraps
//! `Geom2d_OffsetCurve` and `Geom2d_TrimmedCurve` to their basis before the
//! type dispatch, which the 2-D adaptor port does itself via
//! `curve2d_kind`.

use occt_core::gp::{GpDir, GpPnt2d, GpVec, GpVec2d};
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_geom::surface::Surface;

use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::topopebrep::tool::c2d::C2dF;

/// `GeomAbs_CurveType` as far as the quadric predicate cares
/// (`TopOpeBRepTool_TOOL.cxx:205-228`, `GEOMETRY.cxx:246-278`).
fn curve_kind_is_quad(kind: crate::geom_abs_curve_type::GeomAbsCurveType) -> bool {
    use crate::geom_abs_curve_type::GeomAbsCurveType as T;
    matches!(
        kind,
        T::Line | T::Circle | T::Ellipse | T::Hyperbola | T::Parabola
    )
}

/// `FUN_quadCT(CT)` (`TopOpeBRepTool_GEOMETRY.cxx:205-228`): line, circle,
/// ellipse, hyperbola and parabola are "quadric".
pub fn fun_quad_ct(kind: crate::geom_abs_curve_type::GeomAbsCurveType) -> bool {
    curve_kind_is_quad(kind)
}

/// `FUN_tool_line(const Handle(Geom_Curve)&)` (`GEOMETRY.cxx:196-202`): the
/// basis curve's adaptor type is `GeomAbs_Line`.
pub fn fun_tool_line(curve: &dyn occt_geom::curve::Curve) -> bool {
    crate::geom_abs_curve_type::curve_kind(curve)
        == crate::geom_abs_curve_type::GeomAbsCurveType::Line
}

/// `FUN_tool_quad(const Handle(Geom_Curve)&)` (`GEOMETRY.cxx:232-243`).
pub fn fun_tool_quad(curve: &dyn occt_geom::curve::Curve) -> bool {
    curve_kind_is_quad(crate::geom_abs_curve_type::curve_kind(curve))
}

/// `FUN_tool_line(const Handle(Geom2d_Curve)&)` (`GEOMETRY.cxx:280-291`): the
/// basis p-curve's adaptor type is `GeomAbs_Line`.
pub fn fun_tool_line_2d(pc: &dyn occt_geom2d::curve::Curve2d) -> bool {
    use crate::geom_abs_curve_type::GeomAbsCurveType;
    let base = super::iso::basis_curve2d(pc);
    crate::geom_abs_curve_type::curve2d_kind(base) == GeomAbsCurveType::Line
}

/// `FUN_tool_quad(const Handle(Geom2d_Curve)&)` (`GEOMETRY.cxx:246-278`): line,
/// circle, ellipse, hyperbola and parabola p-curves are "quadric".
pub fn fun_tool_quad_2d(pc: &dyn occt_geom2d::curve::Curve2d) -> bool {
    let base = super::iso::basis_curve2d(pc);
    curve_kind_is_quad(crate::geom_abs_curve_type::curve2d_kind(base))
}

/// `TopOpeBRepTool_TOOL::Tg2d(iv, E, C2DF)` (`TOOL.cxx:718-730`): the unit
/// tangent of the p-curve of `e` at the parameter of its `iv` vertex.
pub fn tg_2d(iv: i32, e: &crate::shape::TopoShape, c2df: &C2dF) -> GpVec2d {
    let Some(pc) = c2df.curve.clone() else {
        return GpVec2d::new(0.0, 0.0);
    };
    let Some(par) = super::closed::par_e(iv, e) else {
        return GpVec2d::new(0.0, 0.0);
    };
    let (_p, tg2d) = pc.d1(par);
    tg2d.normalized().unwrap_or(tg2d)
}

/// `TopOpeBRepTool_TOOL::Tg2dApp(iv, E, C2DF, factor)` (`TOOL.cxx:734-757`).
///
/// The tangent is taken at `(1 - factor) * par(iv) + factor * par(iOOv)`, so a
/// quadric p-curve whose end tangent is degenerate still gives a direction.
/// `cxx:752-755` removed the `M_REVERSED` reversal on purpose: the caller wants
/// the original tangent.
pub fn tg_2d_app(
    iv: i32,
    e: &crate::shape::TopoShape,
    c2df: &C2dF,
    factor: f64,
) -> GpVec2d {
    let Some(pc) = c2df.curve.clone() else {
        return GpVec2d::new(0.0, 0.0);
    };
    let i_oo_v = if iv == 1 { 2 } else { 1 };
    let (Some(par), Some(oo_par)) = (
        super::closed::par_e(iv, e),
        super::closed::par_e(i_oo_v, e),
    ) else {
        return GpVec2d::new(0.0, 0.0);
    };
    let par_e = (1.0 - factor) * par + factor * oo_par;
    let (_p, tg2d) = pc.d1(par_e);
    tg2d.normalized().unwrap_or(tg2d)
}

/// `TopOpeBRepTool_TOOL::tryTg2dApp(iv, E, C2DF, factor)` (`TOOL.cxx:761-775`):
/// `Tg2dApp` only for a quadric non-line p-curve, `Tg2d` otherwise.
pub fn try_tg_2d_app(
    iv: i32,
    e: &crate::shape::TopoShape,
    c2df: &C2dF,
    factor: f64,
) -> GpVec2d {
    let Some(pc) = c2df.curve.clone() else {
        return GpVec2d::new(0.0, 0.0);
    };
    let is_quad = fun_tool_quad_2d(pc.as_ref());
    let line = fun_tool_line_2d(pc.as_ref());
    if !is_quad || line {
        return tg_2d(iv, e, c2df);
    }
    tg_2d_app(iv, e, c2df, factor)
}

/// `FUN_tool_quad(const Handle(Geom_Surface)&)` (`GEOMETRY.cxx:294-325`): plane,
/// cylinder, cone, sphere and torus are "quadric".
pub fn fun_tool_quad_surface(s: &dyn Surface) -> bool {
    matches!(
        surface_kind(s),
        GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Cone
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Torus
    )
}

/// `FUN_tool_plane(F)` (`TopOpeBRepTool_TOPOLOGY.cxx:681-686`): the basis
/// surface's adaptor type is `GeomAbs_Plane`.
pub fn fun_tool_plane(f: &crate::shape::Face) -> bool {
    let Some(s) = crate::brep_tool::BRepTool::face_surface_world(f) else {
        return false;
    };
    let s = super::closed::basis_surface(&s);
    surface_kind(s.as_ref()) == GeomAbsSurfaceType::Plane
}

/// `FUN_tool_cylinder(F)` (`TopOpeBRepTool_TOPOLOGY.cxx:689-694`).
pub fn fun_tool_cylinder(f: &crate::shape::Face) -> bool {
    let Some(s) = crate::brep_tool::BRepTool::face_surface_world(f) else {
        return false;
    };
    let s = super::closed::basis_surface(&s);
    surface_kind(s.as_ref()) == GeomAbsSurfaceType::Cylinder
}

/// `FUN_tool_closedS(F, uclosed, uperiod, vclosed, vperiod)`
/// (`TopOpeBRepTool_TOPOLOGY.cxx:697-711`): [`fun_tool_closed`] on the basis
/// surface of `F`. A face without a surface answers
/// `(false, false, 0, false, 0)` at the `S.IsNull()` gate (`cxx:705-708`).
pub fn fun_tool_closed_s(f: &crate::shape::Face) -> (bool, bool, f64, bool, f64) {
    let Some(s) = crate::brep_tool::BRepTool::face_surface_world(f) else {
        return (false, false, 0.0, false, 0.0);
    };
    let s = super::closed::basis_surface(&s);
    fun_tool_closed(s.as_ref())
}

/// `FUN_tool_closedS(F)` (`TOPOLOGY.cxx:714-720`): does `F` close in `U` or
/// `V`?
pub fn fun_tool_is_closed_s(f: &crate::shape::Face) -> bool {
    fun_tool_closed_s(f).0
}

/// `TopOpeBRepTool_TOOL::TolUV(F, tol3d)` (`TOOL.cxx:1925-1931`): the larger
/// of `BRepAdaptor_Surface::UResolution` and `VResolution` at `tol3d`.
pub fn tol_uv(f: &crate::shape::Face, tol3d: f64) -> f64 {
    let Some(s) = crate::brep_tool::BRepTool::face_surface_world(f) else {
        return 0.0;
    };
    let tolu = u_resolution(s.as_ref(), tol3d);
    let tolv = v_resolution(s.as_ref(), tol3d);
    tolu.max(tolv)
}

/// `FUN_tool_closed(S, uclosed, uperiod, vclosed, vperiod)`
/// (`GEOMETRY.cxx:327-358`).
///
/// Returns `(closed, uclosed, uperiod, vclosed, vperiod)`. A closed direction
/// that is not periodic is reported as open (`cxx:341`, `:350`).
pub fn fun_tool_closed(s: &dyn Surface) -> (bool, bool, f64, bool, f64) {
    let mut uperiod = 0.0;
    let mut vperiod = 0.0;
    let mut uclosed = s.is_u_closed();
    if uclosed {
        uclosed = s.is_u_periodic();
    }
    if uclosed {
        uperiod = s.u_period();
    }
    let mut vclosed = s.is_v_closed();
    if vclosed {
        vclosed = s.is_v_periodic();
    }
    if vclosed {
        vperiod = s.v_period();
    }
    let closed = uclosed || vclosed;
    (closed, uclosed, uperiod, vclosed, vperiod)
}

/// `FUN_tool_onapex(p2d, S)` (`GEOMETRY.cxx:86-110`).
///
/// The cone arm compares the 3-D distance from `S(p2d)` to `Cone().Apex()` with
/// `Precision::Confusion()`; the sphere arm tests `|v|` against `pi/2` within
/// `1.e-8`.
pub fn fun_tool_onapex(p2d: &GpPnt2d, s: &dyn Surface) -> bool {
    let tol = occt_core::precision::CONFUSION;
    let toluv = 1.0e-8;
    match surface_kind(s) {
        GeomAbsSurfaceType::Cone => match s.gp_cone() {
            Some(cone) => {
                let apex = cone.apex();
                let pnt = s.d0(p2d.x(), p2d.y());
                pnt.distance(&apex) < tol
            }
            None => false,
        },
        GeomAbsSurfaceType::Sphere => {
            let pisur2 = std::f64::consts::PI * 0.5;
            let v = p2d.y();
            (v - pisur2).abs() < toluv || (v + pisur2).abs() < toluv
        }
        _ => false,
    }
}

/// `FUN_tool_ngS(p2d, S)` (`GEOMETRY.cxx:113-192`).
///
/// The degenerate first-derivative branches (`du < tol` or `dv < tol`) are
/// translated for the cone and sphere apexes; the remaining `kpart` cases keep
/// the C++ `return gp_Dir(gp_Dir::D::Z)` (`cxx:185`), whose `NYI` comment
/// (`cxx:116-118`, `:183-185`) says all non-direct geometries are unhandled.
pub fn fun_tool_ng_s(p2d: &GpPnt2d, s: &dyn Surface) -> GpDir {
    let (p, d1u, d1v) = s.d1(p2d.x(), p2d.y());
    let _ = p;
    let du = d1u.magnitude();
    let dv = d1v.magnitude();
    let tol = occt_core::precision::CONFUSION;
    let kpart = du < tol || dv < tol;
    let z = GpDir::new(0.0, 0.0, 1.0).expect("unit Z");

    if kpart {
        let toluv = 1.0e-8;
        match surface_kind(s) {
            GeomAbsSurfaceType::Cone => {
                let Some(cone) = s.gp_cone() else {
                    return z;
                };
                let nullx = p2d.x().abs() < toluv;
                let apex = nullx && p2d.y().abs() < toluv;
                if apex {
                    let ax1 = cone.axis();
                    let axis = ax1.direction();
                    let mut ng = GpVec::new(axis.x(), axis.y(), axis.z());
                    ng.reverse();
                    return GpDir::from_vec(&ng).unwrap_or(z);
                }
                if du < tol {
                    let (vf, _) = s.v_range();
                    let onvf = (p2d.y() - vf).abs() < toluv;
                    let x = p2d.x();
                    let y = if onvf { p2d.y() + 1.0 } else { p2d.y() - 1.0 };
                    let (_, d1u, d1v) = s.d1(x, y);
                    let ng = d1u.crossed(&d1v);
                    return GpDir::from_vec(&ng).unwrap_or(z);
                }
            }
            GeomAbsSurfaceType::Sphere => {
                let pisur2 = std::f64::consts::PI * 0.5;
                let (u, v) = (p2d.x(), p2d.y());
                let vpisur2 = (v - pisur2).abs() < toluv;
                let vmoinspisur2 = (v + pisur2).abs() < toluv;
                if vpisur2 || vmoinspisur2 {
                    let Some(sphere) = s.gp_sphere() else {
                        return z;
                    };
                    let center = sphere.location();
                    let value = s.d0(u, v);
                    let ng = GpVec::from_pnts(&center, &value);
                    return GpDir::from_vec(&ng).unwrap_or(z);
                }
            }
            _ => {}
        }
        return z;
    }

    let udir = match GpDir::from_vec(&d1u) {
        Ok(d) => d,
        Err(_) => return z,
    };
    let vdir = match GpDir::from_vec(&d1v) {
        Ok(d) => d,
        Err(_) => return z,
    };
    match udir.crossed(&vdir) {
        Ok(d) => d,
        Err(_) => z,
    }
}
