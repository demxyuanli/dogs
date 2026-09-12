//! Bounding box of a 2D Bezier curve from its control polygon.
//!
//! Source: `GeomBndLib_BezierCurve2d.cxx` `Box`. When `[U1,U2]` is a proper
//! sub-range of `[0,1]`, OCCT copies the curve and `Segment`s it; this port
//! evaluates the poles of the current curve when the range covers the whole
//! parameter interval (the pcurve case) and otherwise falls back to
//! `GeomBndLib_OtherCurve2d::Box` sampling, which is the same path
//! `BndLib_Add2dCurve` takes for a non-`Geom2dAdaptor` adaptor.

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpPnt2d;
use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dBezierCurve;

use crate::geom_bnd_lib_other2d::box_other;

/// `GeomBndLib_BezierCurve2d::Box(theU1, theU2, theTol)` on an owned Bezier.
pub fn box_bezier_poles(curve: &Geom2dBezierCurve, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    let a_u1 = 0.0;
    let a_u2 = 1.0;
    let a_trim1 = the_u1.max(a_u1);
    let a_trim2 = the_u2.min(a_u2);
    if (a_u1 - a_trim1).abs() > PCONFUSION || (a_u2 - a_trim2).abs() > PCONFUSION {
        return box_other(curve, a_trim1, a_trim2, the_tol);
    }
    let mut a_box = BndBox2d::new();
    for i in 0..curve.nb_poles() {
        a_box.add_point(&curve.pole(i));
    }
    a_box.enlarge(the_tol);
    a_box
}

/// Pole hull from a raw pole slice (when the caller already segmented).
pub fn box_bezier_pole_slice(poles: &[GpPnt2d], the_tol: f64) -> BndBox2d {
    let mut a_box = BndBox2d::new();
    for p in poles {
        a_box.add_point(p);
    }
    a_box.enlarge(the_tol);
    a_box
}

/// Bezier box from the `Curve2d` trait: sample poles via `d0` at Greville-like
/// nodes is not OCCT; without a downcast the OtherCurve2d `Box` path is used.
pub fn box_bezier_as_curve(curve: &dyn Curve2d, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    box_other(curve, the_u1, the_u2, the_tol)
}
