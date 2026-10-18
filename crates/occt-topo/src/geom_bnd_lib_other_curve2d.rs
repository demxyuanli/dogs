//! Bounding box of a general 2D curve (`GeomAbs_OtherCurve`).
//!
//! Source: `GeomBndLib_OtherCurve2d.cxx` `Box` (not `BoxOptimal`).
//! `BndLib_Add2dCurve::Add(Geom2d_Curve, T1, T2, Tol, Box)` calls
//! `GeomBndLib_Curve2d::Add` which uses `Box(U1,U2,Tol)` — PerformAreas UV
//! bounds therefore take this sampling path for non-analytic pcurves:
//!
//! ```text
//! weakness = 1.5
//! N        = 33
//! tol      = FillBox2d(B1, curve, U1, U2, N)
//! B1.Enlarge(weakness * tol)
//! B1.Get(x, y, X, Y)
//! Box.Update(x, y, X, Y)
//! Box.Enlarge(theTol)
//! ```
//!
//! `FillBox2d` samples `2N` interior points (two steps of `dp = (Last-First)/(2N)`
//! per iteration) plus the start, and returns the max chordal deflection of the
//! mid-sample against the chord of the two outer samples. `BoxOptimal` (PSO +
//! Brent) is intentionally not ported here: PerformAreas never calls it.
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

/// `FillBox2d` from `GeomBndLib_OtherCurve2d.cxx`.
///
/// Walks `[theFirst, theLast]` with `dp = (Last-First) / (2 N)`. Each of the
/// `N` iterations adds two points and updates the max deflection of the odd
/// sample against the chord of the even pair. Returns that max deflection.
pub fn fill_box2d(
    the_box: &mut BndBox2d,
    the_curve: &dyn Curve2d,
    the_first: f64,
    the_last: f64,
    the_n: i32,
) -> f64 {
    let mut a_p1 = the_curve.d0(the_first);
    the_box.add_point(&a_p1);
    let mut p = the_first;
    let mut dp = the_last - the_first;
    let mut tol = 0.0_f64;
    if dp.abs() > PCONFUSION {
        dp /= 2.0 * the_n as f64;
        for _i in 1..=the_n {
            p += dp;
            let a_p2 = the_curve.d0(p);
            the_box.add_point(&a_p2);
            p += dp;
            let a_p3 = the_curve.d0(p);
            the_box.add_point(&a_p3);
            let a_pc = GpPnt2d::new(0.5 * (a_p1.x() + a_p3.x()), 0.5 * (a_p1.y() + a_p3.y()));
            tol = tol.max(a_pc.distance(&a_p2));
            a_p1 = a_p3;
        }
    } else {
        let a_p1b = the_curve.d0(the_first);
        the_box.add_point(&a_p1b);
        let a_p3 = the_curve.d0(the_last);
        the_box.add_point(&a_p3);
    }
    tol
}

const WEAKNESS: f64 = 1.5;
const N_SAMPLES: i32 = 33;

/// `GeomBndLib_OtherCurve2d::Box(theU1, theU2, theTol)`.
pub fn box_other(the_curve: &dyn Curve2d, the_u1: f64, the_u2: f64, the_tol: f64) -> BndBox2d {
    let mut a_box = BndBox2d::new();
    let mut a_b1 = BndBox2d::new();
    let tol = fill_box2d(&mut a_b1, the_curve, the_u1, the_u2, N_SAMPLES);
    a_b1.enlarge(WEAKNESS * tol);
    if let Some((x, y, xx, yy)) = a_b1.get() {
        a_box.update(x, y, xx, yy);
    }
    a_box.enlarge(the_tol);
    a_box
}

/// `GeomBndLib_OtherCurve2d::Box(theTol)` — full parameter range.
pub fn box_other_full(the_curve: &dyn Curve2d, the_tol: f64) -> BndBox2d {
    box_other(
        the_curve,
        the_curve.first_parameter(),
        the_curve.last_parameter(),
        the_tol,
    )
}

/// Adaptor-style 33-point box used when `BndLib_Add2dCurve::Add` cannot
/// downcast to `Geom2dAdaptor_Curve` (`BndLib_Add2dCurve.cxx:38-56`):
/// `N = 33`, `DU = (U2-U1)/(N-1)`, add `D0` at each node, then `Enlarge(Tol)`.
pub fn box_adaptor_sample(the_curve: &dyn Curve2d, a_u1: f64, a_u2: f64, a_tol: f64) -> BndBox2d {
    let mut a_box2d = BndBox2d::new();
    let n = 33;
    let mut u = a_u1;
    let du = if n > 1 { (a_u2 - a_u1) / (n - 1) as f64 } else { 0.0 };
    for _j in 1..n {
        a_box2d.add_point(&the_curve.d0(u));
        u += du;
    }
    a_box2d.add_point(&the_curve.d0(a_u2));
    a_box2d.enlarge(a_tol);
    a_box2d
}
