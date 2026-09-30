//! Add a 2D curve into a `Bnd_Box2d`.
//!
//! Source: `BndLib_Add2dCurve.cxx`. Four overloads:
//! * `Add(Adaptor2d, Tol, Box)` → `Add(C, First, Last, Tol, Box)`;
//! * `Add(Adaptor2d, U1, U2, Tol, Box)` — if the adaptor is not a
//!   `Geom2dAdaptor_Curve`, sample `N = 33` and enlarge; otherwise dispatch
//!   to the Geom2d handle overload;
//! * `Add(Geom2d_Curve, Tol, Box)` → first/last parameters;
//! * `Add(Geom2d_Curve, T1, T2, Tol, Box)` → `GeomBndLib_Curve2d(C).Add(...)`.
//!
//! PerformAreas UV bounds call `BndLib_Add2dCurve::Add(aC2D, aT1, aT2, 0., aBoxC)`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;
use occt_geom2d::curve::Curve2d;

use crate::geom_bnd_lib_curve2d::add_curve2d;
use crate::geom_bnd_lib_other2d::box_adaptor_sample;

/// `BndLib_Add2dCurve::Add(Adaptor2d_Curve2d, Tol, Box)`.
pub fn add_adaptor(curve: &dyn Curve2d, a_tol: f64, a_box2d: &mut BndBox2d) {
    add_adaptor_range(
        curve,
        curve.first_parameter(),
        curve.last_parameter(),
        a_tol,
        a_box2d,
    );
}

/// `BndLib_Add2dCurve::Add(Adaptor2d_Curve2d, U1, U2, Tol, Box)`.
///
/// The port always has a `Curve2d` (the Geom2d handle), so it takes the
/// Geom2d overload rather than the 33-point adaptor fallback. Call
/// [`add_adaptor_sample`] explicitly when the curve is not a Geom2d handle.
pub fn add_adaptor_range(
    curve: &dyn Curve2d,
    a_u1: f64,
    a_u2: f64,
    a_tol: f64,
    a_box2d: &mut BndBox2d,
) {
    add_geom2d_range(curve, a_u1, a_u2, a_tol, a_box2d);
}

/// The 33-point adaptor fallback (`BndLib_Add2dCurve.cxx:38-56`).
pub fn add_adaptor_sample(
    curve: &dyn Curve2d,
    a_u1: f64,
    a_u2: f64,
    a_tol: f64,
    a_box2d: &mut BndBox2d,
) {
    a_box2d.add_box(&box_adaptor_sample(curve, a_u1, a_u2, a_tol));
}

/// `BndLib_Add2dCurve::Add(Geom2d_Curve, Tol, Box)`.
pub fn add_geom2d(curve: &dyn Curve2d, a_tol: f64, a_box2d: &mut BndBox2d) {
    add_geom2d_range(
        curve,
        curve.first_parameter(),
        curve.last_parameter(),
        a_tol,
        a_box2d,
    );
}

/// `BndLib_Add2dCurve::Add(Geom2d_Curve, T1, T2, Tol, Box)`.
pub fn add_geom2d_range(
    curve: &dyn Curve2d,
    a_t1: f64,
    a_t2: f64,
    a_tol: f64,
    a_box2d: &mut BndBox2d,
) {
    add_curve2d(curve, a_t1, a_t2, a_tol, a_box2d);
}

/// `BndLib_Add2dCurve::AddOptimal` is not used by PerformAreas (`Add(..., 0.)`
/// is `Box()`, not `BoxOptimal`). Kept as an alias of [`add_geom2d_range`] so
/// a later caller matching the OCCT `AddOptimal` name does not invent a
/// different sampling budget.
pub fn add_geom2d_optimal(
    curve: &dyn Curve2d,
    a_t1: f64,
    a_t2: f64,
    a_tol: f64,
    a_box2d: &mut BndBox2d,
) {
    add_geom2d_range(curve, a_t1, a_t2, a_tol, a_box2d);
}
