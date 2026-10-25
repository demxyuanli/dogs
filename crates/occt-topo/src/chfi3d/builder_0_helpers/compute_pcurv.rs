//! `ChFi3d_ComputePCurv` — the straight p-curve between two UV points, made
//! `SameParameter` with the 3D curve it belongs to.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx` (OCCT
//! 8.0.0):
//!
//! - `ChFi3d_ComputePCurv(UV1, UV2, Pcurv, Pardeb, Parfin, reverse)`
//!   (`cxx:1690-1747`), the four degenerate-line arms and the two-pole
//!   B-spline fallback,
//! - `ChFi3d_ComputePCurv(C3d, UV1, UV2, Pcurv, S, Pardeb, Parfin, tol3d,
//!   tolreached, reverse)` (`cxx:1647-1661`), the `Adaptor3d_Curve` /
//!   `Adaptor3d_Surface` overload,
//! - `ChFi3d_ComputePCurv(C3d, ..., S, ...)` (`cxx:1666-1685`), the
//!   `Geom_Curve` / `Geom_Surface` overload, which only wraps its arguments
//!   into the two adaptors above.
//!
//! Translation notes:
//! - the three overloads become [`ch_fi3d_compute_pcurv_uv`] (the 2D builder)
//!   and [`ch_fi3d_compute_pcurv`] (build, then `SameParameter`); the
//!   `Geom_Surface` overload differs from the adaptor one only by a
//!   `GeomAdaptor_Surface` wrapper that no longer exists in the port, and the
//!   `GeomAdaptor_Curve(C3d, Pardeb, Parfin)` wrapper is
//!   [`crate::topopebrep::tool::adaptor::GeomAdaptorCurve`], which the caller
//!   passes as the `c3d` argument,
//! - `occ::handle<Geom2d_Curve>& Pcurv` is in/out in OCCT because
//!   `SameParameter` may replace the curve; the port returns the tolerance
//!   reached and writes the curve back through `&mut Option<..>`,
//! - the boolean answer of `ChFi3d_SameParameter` is ignored by `cxx:1660`
//!   exactly like here.

use std::sync::Arc;

use occt_core::gp::{GpAx2d, GpDir2d, GpPnt2d};
use occt_core::precision::PCONFUSION;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::chfi3d::builder_helpers::same_parameter::same_parameter;

/// `ChFi3d_ComputePCurv(UV1, UV2, Pcurv, Pardeb, Parfin, reverse)`
/// (`cxx:1690-1747`).
///
/// The four `if` / `else if` arms cover the p-curves that a fillet boundary
/// has when it follows an iso-line of the surface with the very parameter
/// range of the 3D curve; anything else becomes a degree-1 B-spline with the
/// two poles and the double knots `[Pardeb, Pardeb, Parfin, Parfin]`. The
/// result is always a `Geom2d_TrimmedCurve` over `[Pardeb, Parfin]`
/// (`cxx:1747`).
pub fn ch_fi3d_compute_pcurv_uv(
    uv1: &GpPnt2d,
    uv2: &GpPnt2d,
    pardeb: f64,
    parfin: f64,
    reverse: bool,
) -> Result<Arc<dyn Curve2d>, String> {
    let tol = PCONFUSION;
    // `gp_Pnt2d p1, p2; if (!reverse) {p1=UV1;p2=UV2;} else {p1=UV2;p2=UV1;}`
    // (`cxx:1697-1708`).
    let (p1, p2) = if !reverse { (uv1, uv2) } else { (uv2, uv1) };
    let dx = p1.x() - p2.x();
    let dy = p1.y() - p2.y();
    let span = parfin - pardeb;

    let x2d = || GpDir2d { x: 1.0, y: 0.0 };
    let y2d = || GpDir2d { x: 0.0, y: 1.0 };

    // `cxx:1710-1714`.
    if dx.abs() <= tol && ((p2.y() - p1.y()) - span).abs() <= tol {
        let ppp = GpPnt2d::new(p1.x(), p1.y() - pardeb);
        let line = Geom2dLine::new(GpAx2d::new(ppp, y2d()));
        return Ok(Arc::new(Geom2dTrimmedCurve::new(
            Arc::new(line),
            pardeb,
            parfin,
        )));
    }
    // `cxx:1715-1720`.
    if dx.abs() <= tol && ((p1.y() - p2.y()) - span).abs() <= tol {
        let ppp = GpPnt2d::new(p1.x(), p1.y() + pardeb);
        let line = Geom2dLine::new(GpAx2d::new(ppp, y2d().reversed()));
        return Ok(Arc::new(Geom2dTrimmedCurve::new(
            Arc::new(line),
            pardeb,
            parfin,
        )));
    }
    // `cxx:1721-1726`.
    if dy.abs() <= tol && ((p2.x() - p1.x()) - span).abs() <= tol {
        let ppp = GpPnt2d::new(p1.x() - pardeb, p1.y());
        let line = Geom2dLine::new(GpAx2d::new(ppp, x2d()));
        return Ok(Arc::new(Geom2dTrimmedCurve::new(
            Arc::new(line),
            pardeb,
            parfin,
        )));
    }
    // `cxx:1727-1732`.
    if dy.abs() <= tol && ((p1.x() - p2.x()) - span).abs() <= tol {
        let ppp = GpPnt2d::new(p1.x() + pardeb, p1.y());
        let line = Geom2dLine::new(GpAx2d::new(ppp, x2d().reversed()));
        return Ok(Arc::new(Geom2dTrimmedCurve::new(
            Arc::new(line),
            pardeb,
            parfin,
        )));
    }

    // `NCollection_Array1<gp_Pnt2d> p(1, 2); k(1) = Pardeb; k(2) = Parfin;
    // m.Init(2); new Geom2d_BSplineCurve(p, k, m, 1)` (`cxx:1734-1744`).
    let bs = Geom2dBSplineCurve::new(
        vec![p1.x(), p2.x()],
        vec![p1.y(), p2.y()],
        vec![pardeb, pardeb, parfin, parfin],
        1,
    )
    .map_err(|e| format!("ChFi3d_ComputePCurv: Geom2d_BSplineCurve: {e}"))?;
    Ok(Arc::new(Geom2dTrimmedCurve::new(
        Arc::new(bs),
        pardeb,
        parfin,
    )))
}

/// `ChFi3d_ComputePCurv(C3d, UV1, UV2, Pcurv, S, Pardeb, Parfin, tol3d,
/// tolreached, reverse)` (`cxx:1647-1661`).
///
/// `c3d` is the `GeomAdaptor_Curve` the two adaptor / `Geom_Curve` overloads
/// build over `[Pardeb, Parfin]` — the p-curve is trimmed to the same basis
/// parameters, so the adaptor must keep them (see
/// [`crate::topopebrep::tool::adaptor::GeomAdaptorCurve`]). Answers
/// `tolreached`.
pub fn ch_fi3d_compute_pcurv(
    c3d: &dyn Curve,
    uv1: &GpPnt2d,
    uv2: &GpPnt2d,
    pcurv: &mut Option<Arc<dyn Curve2d>>,
    s: &dyn Surface,
    pardeb: f64,
    parfin: f64,
    tol3d: f64,
    reverse: bool,
) -> Result<f64, String> {
    // `ChFi3d_ComputePCurv(UV1, UV2, Pcurv, Pardeb, Parfin, reverse);`
    // (`cxx:1659`).
    let curve = ch_fi3d_compute_pcurv_uv(uv1, uv2, pardeb, parfin, reverse)?;
    // `ChFi3d_SameParameter(C3d, Pcurv, S, tol3d, tolreached);` (`cxx:1660`).
    let mut built = Some(curve);
    let (_, tolreached) = same_parameter(c3d, &mut built, s, tol3d);
    *pcurv = built;
    Ok(tolreached)
}
