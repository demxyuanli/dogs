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
    bnd_box: &mut BndBox2d,
    curve: &dyn Curve2d,
    first: f64,
    last: f64,
    n: i32,
) -> f64 {
    let mut p1 = curve.d0(first);
    bnd_box.add_point(&p1);
    let mut p = first;
    let mut dp = last - first;
    let mut tol = 0.0_f64;
    if dp.abs() > PCONFUSION {
        dp /= 2.0 * n as f64;
        for _i in 1..=n {
            p += dp;
            let p2 = curve.d0(p);
            bnd_box.add_point(&p2);
            p += dp;
            let p3 = curve.d0(p);
            bnd_box.add_point(&p3);
            let pc = GpPnt2d::new(0.5 * (p1.x() + p3.x()), 0.5 * (p1.y() + p3.y()));
            tol = tol.max(pc.distance(&p2));
            p1 = p3;
        }
    } else {
        let p1b = curve.d0(first);
        bnd_box.add_point(&p1b);
        let p3 = curve.d0(last);
        bnd_box.add_point(&p3);
    }
    tol
}

const WEAKNESS: f64 = 1.5;
const N_SAMPLES: i32 = 33;

/// `GeomBndLib_OtherCurve2d::Box(theU1, theU2, theTol)`.
pub fn box_other_2d(curve: &dyn Curve2d, u1: f64, u2: f64, the_tol: f64) -> BndBox2d {
    let mut a_box = BndBox2d::new();
    let mut b1 = BndBox2d::new();
    let tol = fill_box2d(&mut b1, curve, u1, u2, N_SAMPLES);
    b1.enlarge(WEAKNESS * tol);
    if let Some((x, y, xx, yy)) = b1.get() {
        a_box.update(x, y, xx, yy);
    }
    a_box.enlarge(the_tol);
    a_box
}

/// `GeomBndLib_OtherCurve2d::Box(theTol)` — full parameter range.
pub fn box_other_full(curve: &dyn Curve2d, the_tol: f64) -> BndBox2d {
    box_other_2d(
        curve,
        curve.first_parameter(),
        curve.last_parameter(),
        the_tol,
    )
}

/// Adaptor-style 33-point box used when `BndLib_Add2dCurve::Add` cannot
/// downcast to `Geom2dAdaptor_Curve` (`BndLib_Add2dCurve.cxx:38-56`):
/// `N = 33`, `DU = (U2-U1)/(N-1)`, add `D0` at each node, then `Enlarge(Tol)`.
pub fn box_adaptor_sample(curve: &dyn Curve2d, a_u1: f64, a_u2: f64, a_tol: f64) -> BndBox2d {
    let mut box2d = BndBox2d::new();
    let n = 33;
    let mut u = a_u1;
    let du = if n > 1 { (a_u2 - a_u1) / (n - 1) as f64 } else { 0.0 };
    for _j in 1..n {
        box2d.add_point(&curve.d0(u));
        u += du;
    }
    box2d.add_point(&curve.d0(a_u2));
    box2d.enlarge(a_tol);
    box2d
}
