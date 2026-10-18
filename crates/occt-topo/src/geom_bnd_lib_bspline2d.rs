//! Bounding box of a 2D B-spline from its control polygon.
//!
//! Source: `GeomBndLib_BSplineCurve2d.cxx` `Box`. OCCT copies and `Segment`s
//! when `[U1,U2]` is a proper sub-range, then adds every pole of the trimmed
//! control net; otherwise it adds the live poles.
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
use occt_geom2d::Geom2dBSplineCurve;

use crate::geom_bnd_lib_other_curve2d::box_other_2d;

/// `GeomBndLib_BSplineCurve2d::Box(theU1, theU2, theTol)` on an owned B-spline.
///
/// `aCurve->FirstParameter()` / `LastParameter()` (`:33-34`), the trim window
/// (`:35-42`), then `Copy()` + `Segment(aTrim1, aTrim2)` when the window is a
/// proper sub-range (`:44-50`), then the pole hull of the result (`:52-57`).
/// OCCT's `Segment` default tolerance is `Precision::PConfusion()`; the port
/// passes it explicitly (no default arguments in Rust).
pub fn box_bspline_poles(
    curve: &Geom2dBSplineCurve,
    u1: f64,
    u2: f64,
    tol: f64,
) -> BndBox2d {
    let mut curve_val = curve.clone();
    let u1_val = curve_val.first_parameter();
    let u2_val = curve_val.last_parameter();
    let mut trim1 = u1.max(u1_val);
    let mut trim2 = u2.min(u2_val);
    if trim2 < trim1 {
        trim1 = u1_val;
        trim2 = u2_val;
    }
    if (u1_val - trim1).abs() > PCONFUSION || (u2_val - trim2).abs() > PCONFUSION {
        if curve_val.segment(trim1, trim2, PCONFUSION).is_err() {
            return box_other_2d(curve, trim1, trim2, tol);
        }
    }
    let mut bnd_box = BndBox2d::new();
    for i in 0..curve_val.nb_poles() {
        bnd_box.add_point(&GpPnt2d::new(curve_val.xs[i], curve_val.ys[i]));
    }
    bnd_box.enlarge(tol);
    bnd_box
}

/// Pole hull from separate x/y arrays.
pub fn box_bspline_xy(xs: &[f64], ys: &[f64], tol: f64) -> BndBox2d {
    let mut bnd_box = BndBox2d::new();
    let n = xs.len().min(ys.len());
    for i in 0..n {
        bnd_box.add_point(&GpPnt2d::new(xs[i], ys[i]));
    }
    bnd_box.enlarge(tol);
    bnd_box
}

/// `GeomBndLib_BSplineCurve2d::Box(theU1, theU2, theTol)` from the `Curve2d`
/// trait without a downcast.
///
/// `GeomBndLib_Curve2d.cxx:150-157` feeds the evaluator the handle itself, or
/// `Geom2d_TrimmedCurve::BasisCurve()` when the handle is trimmed; the
/// evaluator is `GeomBndLib_BSplineCurve2d` only when that curve is a
/// `Geom2d_BSplineCurve` (`:137-143`). [`Curve2d::bspline_copy2d`] exposes
/// exactly that curve; when it is `None`, OCCT picked a different evaluator
/// (Line / Circle / ... / OtherCurve2d) and the port keeps the
/// `GeomBndLib_OtherCurve2d::Box` fallback.
pub fn box_bspline_as_curve(curve: &dyn Curve2d, u1: f64, u2: f64, tol: f64) -> BndBox2d {
    match curve.bspline_copy2d() {
        Some(curve_val) => box_bspline_poles(&curve_val, u1, u2, tol),
        None => box_other_2d(curve, u1, u2, tol),
    }
}
