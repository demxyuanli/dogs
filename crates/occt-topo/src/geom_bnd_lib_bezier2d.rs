//! Bounding box of a 2D Bezier curve from its control polygon.
//!
//! Source: `GeomBndLib_BezierCurve2d.cxx` `Box`. When `[U1,U2]` is a proper
//! sub-range of `[0,1]`, OCCT copies the curve and `Segment`s it; this port
//! evaluates the poles of the current curve when the range covers the whole
//! parameter interval (the pcurve case) and otherwise falls back to
//! `GeomBndLib_OtherCurve2d::Box` sampling, which is the same path
//! `BndLib_Add2dCurve` takes for a non-`Geom2dAdaptor` adaptor.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpPnt2d;
use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dBezierCurve;

use crate::geom_bnd_lib_other_curve2d::box_other_2d;

/// `GeomBndLib_BezierCurve2d::Box(theU1, theU2, theTol)` on an owned Bezier.
pub fn box_bezier_poles(curve: &Geom2dBezierCurve, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    let a_u1 = 0.0;
    let a_u2 = 1.0;
    let trim1 = u1.max(a_u1);
    let trim2 = u2.min(a_u2);
    if (a_u1 - trim1).abs() > PCONFUSION || (a_u2 - trim2).abs() > PCONFUSION {
        return box_other_2d(curve, trim1, trim2, tol);
    }
    let mut bnd_box = BndBox2d::new();
    for i in 0..curve.nb_poles() {
        bnd_box.add_point(&curve.pole(i));
    }
    bnd_box.enlarge(tol);
    bnd_box
}

/// Pole hull from a raw pole slice (when the caller already segmented).
pub fn box_bezier_pole_slice(poles: &[GpPnt2d], tol: f64) -> BndBox2d {
    let mut bnd_box = BndBox2d::new();
    for p in poles {
        bnd_box.add_point(p);
    }
    bnd_box.enlarge(tol);
    bnd_box
}

/// Bezier box from the `Curve2d` trait.
///
/// `GeomBndLib_BezierCurve2d::Box(theU1, theU2, theTol)`
/// (`GeomBndLib_BezierCurve2d.cxx:30-52`): when `[theU1, theU2]` is a proper
/// sub-range of the curve's own range, OCCT copies the curve and `Segment`s it
/// (`:38-43`); this port has no 2D Bezier `Segment`, so that arm falls back to
/// `GeomBndLib_OtherCurve2d::Box` sampling, and otherwise the control-polygon
/// hull is the box.
pub fn box_bezier_as_curve(curve: &dyn Curve2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    let a_u1 = curve.first_parameter();
    let a_u2 = curve.last_parameter();
    let trim1 = u1.max(a_u1);
    let trim2 = u2.min(a_u2);
    if (a_u1 - trim1).abs() > PCONFUSION || (a_u2 - trim2).abs() > PCONFUSION {
        return box_other_2d(curve, trim1, trim2, tol);
    }
    match curve.poles2d() {
        Some(poles) => box_bezier_pole_slice(&poles, tol),
        None => box_other_2d(curve, trim1, trim2, tol),
    }
}
