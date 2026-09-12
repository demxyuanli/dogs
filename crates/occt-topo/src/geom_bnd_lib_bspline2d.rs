//! Bounding box of a 2D B-spline from its control polygon.
//!
//! Source: `GeomBndLib_BSplineCurve2d.cxx` `Box`. OCCT copies and `Segment`s
//! when `[U1,U2]` is a proper sub-range, then adds every pole. This port adds
//! the live poles when the range covers the knot span; otherwise it uses
//! `GeomBndLib_OtherCurve2d::Box` (the adaptor fallback).

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpPnt2d;
use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dBSplineCurve;

use crate::geom_bnd_lib_other2d::box_other;

/// `GeomBndLib_BSplineCurve2d::Box(theU1, theU2, theTol)` on an owned B-spline.
pub fn box_bspline_poles(
    curve: &Geom2dBSplineCurve,
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> BndBox2d {
    let a_u1 = curve.first_parameter();
    let a_u2 = curve.last_parameter();
    let mut a_trim1 = the_u1.max(a_u1);
    let mut a_trim2 = the_u2.min(a_u2);
    if a_trim2 < a_trim1 {
        a_trim1 = a_u1;
        a_trim2 = a_u2;
    }
    if (a_u1 - a_trim1).abs() > PCONFUSION || (a_u2 - a_trim2).abs() > PCONFUSION {
        return box_other(curve, a_trim1, a_trim2, the_tol);
    }
    let mut a_box = BndBox2d::new();
    for i in 0..curve.nb_poles() {
        a_box.add_point(&GpPnt2d::new(curve.xs[i], curve.ys[i]));
    }
    a_box.enlarge(the_tol);
    a_box
}

/// Pole hull from separate x/y arrays.
pub fn box_bspline_xy(xs: &[f64], ys: &[f64], the_tol: f64) -> BndBox2d {
    let mut a_box = BndBox2d::new();
    let n = xs.len().min(ys.len());
    for i in 0..n {
        a_box.add_point(&GpPnt2d::new(xs[i], ys[i]));
    }
    a_box.enlarge(the_tol);
    a_box
}

/// B-spline box from the `Curve2d` trait without a downcast: OtherCurve2d Box.
pub fn box_bspline_as_curve(curve: &dyn Curve2d, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    box_other(curve, the_u1, the_u2, the_tol)
}
