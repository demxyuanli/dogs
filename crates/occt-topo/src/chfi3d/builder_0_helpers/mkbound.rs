//! `ChFi3d_mkbound` / `ChFi3d_BuildPCurve` — the `GeomFill_Boundary` builders
//! of the twisted-surface repair.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:1750-1967`
//! (OCCT 8.0.0):
//!
//! - `ChFi3d_mkbound(Fac, curv, sens1, pfac1, vfac1, sens2, pfac2, vfac2,
//!   t3d, ta)` (`cxx:1750-1774`),
//! - `ChFi3d_mkbound(Surf, curv, sens1, p1, v1, sens2, p2, v2, t3d, ta)`
//!   (`cxx:1777-1798`),
//! - `ChFi3d_mkbound(Geom_Surface, p1, p2, t3d, ta, isfreeboundary)`
//!   (`cxx:1802-1812`),
//! - `ChFi3d_mkbound(Adaptor3d_Surface, p1, p2, t3d, ta, isfreeboundary)`
//!   (`cxx:1815-1827`),
//! - `ChFi3d_mkbound(Adaptor3d_Surface, curv, t3d, ta, isfreeboundary)`
//!   (`cxx:1831-1846`), the one every other overload funnels into,
//! - `ChFi3d_mkbound(Fac, curv, p1, p2, t3d, ta, isfreeboundary)`
//!   (`cxx:1849-1862`),
//! - `ChFi3d_BuildPCurve(p1, d1, p2, d2, redresse)` (`cxx:1866-1897`),
//! - `ChFi3d_BuildPCurve(Surf, p1, v1, p2, v2, redresse)` (`cxx:1902-1935`),
//! - `ChFi3d_BuildPCurve(Surf, p1, v1, p2, v2, redresse)` on 3D vectors
//!   (`cxx:1937-1966`).
//!
//! Translation notes:
//! - Rust has no overloading, so the six `mkbound` overloads are
//!   [`mkbound_pfac2d`], [`mkbound_pfac3d`], [`mkbound_segment`],
//!   [`mkbound_segment_adaptor`], [`mkbound_curve`] and
//!   [`mkbound_curve_segment`], named after the argument pair they take,
//! - OCCT's `curv` out-parameter is part of the return value: the p-curve the
//!   function built is what the caller keeps (`PerformSetOfKGen` passes `PC1`
//!   to `CompleteData` afterwards),
//! - `occ::handle<Adaptor3d_Surface>` becomes `&SurfaceAdaptor` and
//!   `occ::handle<Geom_Surface>` becomes `&Arc<dyn Surface>`; the
//!   `GeomAdaptor_Surface(s)` wrapper of `cxx:1809` exists only to reach
//!   `Adaptor3d_Surface`, and its `UResolution` / `VResolution` are never read
//!   by the `mkbound` body, so the surface goes straight through,
//! - `Adaptor3d_CurveOnSurface(HC, HS)` is
//!   [`occt_geom::curve_on_surface::CurveOnSurface`], which the ported
//!   `GeomFill_BoundWithSurf` takes apart into `(curve, pcurve, surface)`,
//! - `gp_Dir2d` construction can fail (a null vector) where OCCT throws
//!   `Standard_ConstructionError`; that is the `Err`.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec, GpVec2d};
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_fill::boundary::{Boundary, SimpleBound};
use occt_geom::geom_fill::BoundWithSurf;
use occt_geom::surface::Surface;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::curve::Curve2d;

use crate::brep_approx::surface_tool::SurfaceTool;
use crate::chfi3d::coefficient::chfi3d_coefficient;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFi3d_BuildPCurve(p1, d1, p2, d2, redresse)` (`cxx:1866-1897`).
///
/// The "cubique a la mords moi le noeud": a Bezier with four poles whose two
/// inner poles are pushed along `d1` / `d2` by half the chord times the
/// cosine between the two directions.
///
/// `d1` / `d2` are `gp_Dir2d&` in OCCT and can be reversed in place when
/// `redresse` is set; both callers pass locals, so the port takes them by
/// value.
pub fn build_p_curve_dirs(
    p1: &GpPnt2d,
    mut d1: GpDir2d,
    p2: &GpPnt2d,
    mut d2: GpDir2d,
    redresse: bool,
) -> Result<Geom2dBezierCurve, String> {
    // `gp_Vec2d vref(p1, p2); gp_Dir2d dref(vref); double mref = vref.Magnitude();`
    let vref = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
    let dref = GpDir2d::new(vref.x(), vref.y())
        .map_err(|_| "gp_Dir2d: null vector".to_string())?;
    let mref = vref.magnitude();
    if redresse {
        // `cxx:1875-1886`.
        if d1.dot(&dref) < 0.0 {
            d1.reverse();
        }
        if d2.dot(&dref) > 0.0 {
            d2.reverse();
        }
    }
    // `cxx:1889-1896`.
    let mut lambda1 = d2.dot(&d1).abs().max(dref.dot(&d1).abs());
    lambda1 = (0.5 * mref * lambda1).max(1.0e-5);
    let pol2 = p1.translated_vec(&GpVec2d::new(lambda1 * d1.x(), lambda1 * d1.y()));
    let mut lambda2 = d1.dot(&d2).abs().max(dref.dot(&d2).abs());
    lambda2 = (0.5 * mref * lambda2).max(1.0e-5);
    let pol3 = p2.translated_vec(&GpVec2d::new(lambda2 * d2.x(), lambda2 * d2.y()));
    Geom2dBezierCurve::new(vec![*p1, pol2, pol3, *p2]).map_err(|e| e.to_string())
}

/// `ChFi3d_BuildPCurve(Surf, p1, v1, p2, v2, redresse)` on `gp_Vec2d`
/// (`cxx:1902-1935`): the same cubic in the resolution-normalized frame of
/// `Surf`, then scaled back pole by pole.
pub fn build_p_curve_vec2d(
    surf: &SurfaceAdaptor,
    p1: &GpPnt2d,
    v1: &GpVec2d,
    p2: &GpPnt2d,
    v2: &GpVec2d,
    redresse: bool,
) -> Result<Geom2dBezierCurve, String> {
    // `cxx:1910-1922`.
    let ures = SurfaceTool::u_resolution(surf, 1.0);
    let vres = SurfaceTool::v_resolution(surf, 1.0);
    let invures = 1.0 / ures;
    let invvres = 1.0 / vres;
    let pp1 = GpPnt2d::new(invures * p1.x(), invvres * p1.y());
    let pp2 = GpPnt2d::new(invures * p2.x(), invvres * p2.y());
    let vv1 = GpVec2d::new(invures * v1.x(), invvres * v1.y());
    let vv2 = GpVec2d::new(invures * v2.x(), invvres * v2.y());
    let d1 = GpDir2d::new(vv1.x(), vv1.y())
        .map_err(|_| "gp_Dir2d: null vector".to_string())?;
    let d2 = GpDir2d::new(vv2.x(), vv2.y())
        .map_err(|_| "gp_Dir2d: null vector".to_string())?;
    let mut pc = build_p_curve_dirs(&pp1, d1, &pp2, d2, redresse)?;
    // `cxx:1927-1933`: `pc->Pole(ip)` / `pc->SetPole(ip, pol)`, 1-based.
    for i in 0..pc.nb_poles() {
        let pol = pc.pole(i);
        pc.poles[i] = GpPnt2d::new(ures * pol.x(), vres * pol.y());
    }
    Ok(pc)
}

/// `ChFi3d_BuildPCurve(Surf, p1, v1, p2, v2, redresse)` on `gp_Vec`
/// (`cxx:1937-1966`): the 3D tangents become `(du, dv)` through
/// `ChFi3d_Coefficient`, an optional redressing against the chord, and then the
/// `gp_Vec2d` overload.
pub fn build_p_curve_vec3d(
    surf: &SurfaceAdaptor,
    p1: &GpPnt2d,
    v1: &GpVec,
    p2: &GpPnt2d,
    v2: &GpVec,
    redresse: bool,
) -> Result<Geom2dBezierCurve, String> {
    // `cxx:1945-1956`: `Surf->D1(p1.X(), p1.Y(), PP1, D1u, D1v)` then
    // `Surf->D1(p2.X(), p2.Y(), PP2, D1u, D1v)`.
    let (pp1, d1u, d1v) = surf
        .d1(p1.x(), p1.y())
        .ok_or_else(|| "ChFi3d_BuildPCurve: null surface adaptor".to_string())?;
    let (du, dv) = chfi3d_coefficient(v1, &d1u, &d1v);
    let mut vv1 = GpVec2d::new(du, dv);
    let (pp2, d1u, d1v) = surf
        .d1(p2.x(), p2.y())
        .ok_or_else(|| "ChFi3d_BuildPCurve: null surface adaptor".to_string())?;
    let (du, dv) = chfi3d_coefficient(v2, &d1u, &d1v);
    let mut vv2 = GpVec2d::new(du, dv);
    // `cxx:1957-1965`: `gp_Vec Vref(PP1, PP2)`.
    let vref = GpVec::new(pp2.x() - pp1.x(), pp2.y() - pp1.y(), pp2.z() - pp1.z());
    if redresse {
        if vref.dot(v1) < 0.0 {
            vv1.reverse();
        }
        if vref.dot(v2) > 0.0 {
            vv2.reverse();
        }
    }
    build_p_curve_vec2d(surf, p1, &vv1, p2, &vv2, false)
}

/// The surface a `CurveOnSurface` boundary is evaluated on
/// (`Adaptor3d_Surface::Surface()`); a null adaptor is the OCCT null
/// dereference.
fn adaptor_surface(surf: &SurfaceAdaptor) -> Result<Arc<dyn Surface>, String> {
    Ok(Arc::clone(surf.surface().ok_or_else(|| {
        "ChFi3d_mkbound: null Adaptor3d_Surface".to_string()
    })?))
}

/// `ChFi3d_mkbound(HS, curv, t3d, ta, isfreeboundary)` (`cxx:1831-1846`): the
/// overload every other one funnels into.
pub fn mkbound_curve(
    surface: &Arc<dyn Surface>,
    curv: &Arc<dyn Curve2d>,
    t3d: f64,
    ta: f64,
    isfreeboundary: bool,
) -> Arc<dyn Boundary> {
    // `Adaptor3d_CurveOnSurface COnS(HC, HS)` (`cxx:1837-1838`).
    let cons: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(
        Arc::clone(curv),
        Arc::clone(surface),
    ));
    if isfreeboundary {
        // `new GeomFill_SimpleBound(HCOnS, t3d, ta)` (`cxx:1842`).
        Arc::new(SimpleBound::new(cons, t3d, ta))
    } else {
        // `new GeomFill_BoundWithSurf(COnS, t3d, ta)` (`cxx:1844`).
        Arc::new(BoundWithSurf::new(
            cons,
            Arc::clone(curv),
            Arc::clone(surface),
            t3d,
            ta,
        ))
    }
}

/// `ChFi3d_mkbound(Adaptor3d_Surface, p1, p2, t3d, ta, isfreeboundary)`
/// (`cxx:1815-1827`): a two-pole Bezier segment between the two parameters.
pub fn mkbound_segment_adaptor(
    surf: &SurfaceAdaptor,
    p1: &GpPnt2d,
    p2: &GpPnt2d,
    t3d: f64,
    ta: f64,
    isfreeboundary: bool,
) -> Result<Arc<dyn Boundary>, String> {
    let surface = adaptor_surface(surf)?;
    let curv: Arc<dyn Curve2d> = Arc::new(
        Geom2dBezierCurve::new(vec![*p1, *p2]).map_err(|e| e.to_string())?,
    );
    Ok(mkbound_curve(&surface, &curv, t3d, ta, isfreeboundary))
}

/// `ChFi3d_mkbound(Geom_Surface, p1, p2, t3d, ta, isfreeboundary)`
/// (`cxx:1802-1812`).
pub fn mkbound_segment(
    surface: &Arc<dyn Surface>,
    p1: &GpPnt2d,
    p2: &GpPnt2d,
    t3d: f64,
    ta: f64,
    isfreeboundary: bool,
) -> Result<Arc<dyn Boundary>, String> {
    let curv: Arc<dyn Curve2d> = Arc::new(
        Geom2dBezierCurve::new(vec![*p1, *p2]).map_err(|e| e.to_string())?,
    );
    Ok(mkbound_curve(surface, &curv, t3d, ta, isfreeboundary))
}

/// `ChFi3d_mkbound(Fac, curv, p1, p2, t3d, ta, isfreeboundary)`
/// (`cxx:1849-1862`).
pub fn mkbound_curve_segment(
    surf: &SurfaceAdaptor,
    p1: &GpPnt2d,
    p2: &GpPnt2d,
    t3d: f64,
    ta: f64,
    isfreeboundary: bool,
) -> Result<(Arc<dyn Curve2d>, Arc<dyn Boundary>), String> {
    let curv: Arc<dyn Curve2d> = Arc::new(
        Geom2dBezierCurve::new(vec![*p1, *p2]).map_err(|e| e.to_string())?,
    );
    let bound = mkbound_segment_adaptor(surf, p1, p2, t3d, ta, isfreeboundary)?;
    Ok((curv, bound))
}

/// `ChFi3d_mkbound(Fac, curv, sens1, pfac1, vfac1, sens2, pfac2, vfac2, t3d,
/// ta)` (`cxx:1750-1774`): the 2D-tangent overload. Answers the p-curve it
/// built together with the boundary, because OCCT writes it back into `curv`.
#[allow(clippy::too_many_arguments)]
pub fn mkbound_pfac2d(
    surf: &SurfaceAdaptor,
    sens1: i32,
    pfac1: &GpPnt2d,
    vfac1: &GpVec2d,
    sens2: i32,
    pfac2: &GpPnt2d,
    vfac2: &GpVec2d,
    t3d: f64,
    ta: f64,
) -> Result<(Arc<dyn Curve2d>, Arc<dyn Boundary>), String> {
    // `cxx:1760-1770`.
    let mut v1 = GpDir2d::new(vfac1.x(), vfac1.y())
        .map_err(|_| "gp_Dir2d: null vector".to_string())?;
    if sens1 == 1 {
        v1.reverse();
    }
    let mut v2 = GpDir2d::new(vfac2.x(), vfac2.y())
        .map_err(|_| "gp_Dir2d: null vector".to_string())?;
    if sens2 == 1 {
        v2.reverse();
    }
    let pc = build_p_curve_dirs(pfac1, v1, pfac2, v2, false)?;
    let curv: Arc<dyn Curve2d> = Arc::new(pc);
    let bound = mkbound_curve(&adaptor_surface(surf)?, &curv, t3d, ta, false);
    Ok((curv, bound))
}

/// `ChFi3d_mkbound(Surf, curv, sens1, p1, v1, sens2, p2, v2, t3d, ta)`
/// (`cxx:1777-1798`): the 3D-tangent overload.
#[allow(clippy::too_many_arguments)]
pub fn mkbound_pfac3d(
    surf: &SurfaceAdaptor,
    sens1: i32,
    p1: &GpPnt2d,
    v1: &GpVec,
    sens2: i32,
    p2: &GpPnt2d,
    v2: &GpVec,
    t3d: f64,
    ta: f64,
) -> Result<(Arc<dyn Curve2d>, Arc<dyn Boundary>), String> {
    // `cxx:1788-1795`: `gp_Vec` is taken by value because `Reverse` is applied
    // to the caller's copy.
    let mut v1 = *v1;
    if sens1 == 1 {
        v1.reverse();
    }
    let mut v2 = *v2;
    if sens2 == 1 {
        v2.reverse();
    }
    let pc = build_p_curve_vec3d(surf, p1, &v1, p2, &v2, false)?;
    let curv: Arc<dyn Curve2d> = Arc::new(pc);
    let bound = mkbound_curve(&adaptor_surface(surf)?, &curv, t3d, ta, false);
    Ok((curv, bound))
}
