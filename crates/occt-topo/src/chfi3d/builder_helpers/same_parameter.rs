//! `ChFi3d_CheckSameParameter` / `ChFi3d_SameParameter` — the two
//! `SameParameter` wrappers of `ChFi3d_Builder_0.cxx`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:1560-1642`
//! (OCCT 8.0.0).
//!
//! Divergences that are not heuristics:
//! - `occ::handle<Geom2d_Curve>& Pcurv` is an in/out parameter in OCCT because
//!   `SameParameter` replaces the p-curve when the approximation rebuilt it;
//!   the port takes `&mut Option<Arc<dyn Curve2d>>` and rewrites it from
//!   [`ApproxSameParameter::curve2d`], the same assignment as `cxx:1614`,
//! - `Adaptor3d_Curve` / `Adaptor3d_Surface` become the `Curve` / `Surface`
//!   traits the port's adaptors already implement.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;
use occt_geom::approx_same_parameter::ApproxSameParameter;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_geom::curve::Curve;

/// `ChFi3d_CheckSameParameter(C3d, Pcurv, S, tol3d, tolreached)`
/// (`cxx:1565-1595`): samples `Pcurv` against `C3d` at 45 parameters and
/// answers whether the largest 3d gap stays under `tol3d`.
///
/// `tolreached` is doubled on both exits, exactly like the source.
pub fn check_same_parameter(
    c3d: &dyn Curve,
    pcurv: &dyn Curve2d,
    s: &dyn Surface,
    tol3d: f64,
) -> (bool, f64) {
    let mut tolreached: f64 = 0.0;
    let f = c3d.first_parameter();
    let l = c3d.last_parameter();
    let nbp = 45;
    let step = 1.0 / (nbp - 1) as f64;
    for i in 0..nbp {
        let t = step * i as f64;
        let t = (1.0 - t) * f + t * l;
        let uv = pcurv.value(t);
        let ps: GpPnt = s.value(uv.x(), uv.y());
        let pc: GpPnt = c3d.value(t);
        let d2 = ps.square_distance(&pc);
        tolreached = tolreached.max(d2);
    }
    tolreached = tolreached.sqrt();
    if tolreached > tol3d {
        tolreached *= 2.0;
        return (false, tolreached);
    }
    tolreached *= 2.0;
    tolreached = tolreached.max(CONFUSION);
    (true, tolreached)
}

/// `ChFi3d_SameParameter(C3d, Pcurv, S, tol3d, tolreached)`
/// (`cxx:1602-1623`).
pub fn same_parameter(
    c3d: &dyn Curve,
    pcurv: &mut Option<Arc<dyn Curve2d>>,
    s: &dyn Surface,
    tol3d: f64,
) -> (bool, f64) {
    let Some(curve) = pcurv.as_ref() else {
        // OCCT dereferences the null handle in `CheckSameParameter`; a null
        // p-curve never reaches this wrapper there.
        return (false, 0.0);
    };
    let (ok, tolreached) = check_same_parameter(c3d, curve.as_ref(), s, tol3d);
    if ok {
        return (true, tolreached);
    }
    let sp = ApproxSameParameter::new(
        c3d,
        curve.as_ref(),
        s,
        c3d.first_parameter(),
        c3d.last_parameter(),
        tol3d,
    );
    if sp.done && !sp.same_parameter {
        *pcurv = sp.curve2d.clone();
    } else if !sp.done && !sp.same_parameter {
        return (false, sp.tol_reached);
    }
    (true, sp.tol_reached)
}

/// `ChFi3d_SameParameter(C3d, Pcurv, S, Pardeb, Parfin, tol3d, tolreached)`
/// (`cxx:1629-1640`): the `Geom_Curve` / `Geom_Surface` overload loads
/// `GeomAdaptor_Curve(C3d, Pardeb, Parfin)`, so the port trims the curve to
/// `[pardeb, parfin]` and forwards to [`same_parameter`].
pub fn same_parameter_range(
    c3d: &dyn Curve,
    pcurv: &mut Option<Arc<dyn Curve2d>>,
    s: &dyn Surface,
    pardeb: f64,
    parfin: f64,
    tol3d: f64,
) -> (bool, f64) {
    let trimmed = occt_geom::trimmed::GeomTrimmedCurve::new(
        Arc::from(c3d.clone_dyn()),
        pardeb,
        parfin,
    );
    same_parameter(&trimmed, pcurv, s, tol3d)
}
