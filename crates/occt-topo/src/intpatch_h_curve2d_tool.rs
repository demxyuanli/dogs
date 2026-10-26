//! `IntPatch_HCurve2dTool` - the `Adaptor2d_Curve2d` services the `IntPatch`
//! curve/surface intersection codes call.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_HCurve2dTool.hxx`,
//! `.lxx` and `.cxx`.
//!
//! Ported here (OCCT 8.0.0):
//! * `FirstParameter` / `LastParameter` (`lxx:38-45`).
//! * `Continuity` (`lxx:47-50`), `NbIntervals` (`lxx:52-57`), `Intervals`
//!   (`lxx:59-65`).
//! * `IsClosed` (`lxx:67-71`), `IsPeriodic` (`lxx:73-77`), `Period`
//!   (`lxx:79-83`).
//! * `Value` / `D0` / `D1` / `D2` / `D3` / `DN` (`lxx:85-134`).
//! * `Resolution` (`lxx:136-141`).
//! * `GetType` (`lxx:143-147`) plus the analytic accessors `Line` / `Circle` /
//!   `Ellipse` / `Parabola` / `Hyperbola` (`lxx:149-179`), which are the
//!   `gp_*2d` accessors of [`occt_geom2d::Curve2d`].
//! * `NbSamples(C, U0, U1)` (`cxx:20-53`).
//!
//! UNPORTED in this module:
//! * `Bezier` (`lxx:181-187`) and `BSpline` (`lxx:189-195`): the port's
//!   `Curve2d` trait carries the B-spline data (`bspline_copy2d`,
//!   `bspline_distinct_knots_mults`) instead of returning an OCCT handle.
//!
//! NOTE: [`crate::geom_abs_curve_type::curve2d_kind`] never reports
//! `GeomAbs_BezierCurve` (the `Curve2d` trait has no Bezier query, see
//! `occt-geom2d/src/curve.rs`), so the `GeomAbs_BezierCurve` arms below are the
//! faithful OCCT branches but stay unreachable until that query lands.

use occt_geom2d::Curve2d;

use crate::geom_abs_curve_type::{curve2d_kind, GeomAbsCurveType};

/// `FirstParameter(C)` (`lxx:38-41`).
pub fn first_parameter(c: &dyn Curve2d) -> f64 {
    c.first_parameter()
}

/// `LastParameter(C)` (`lxx:42-45`).
pub fn last_parameter(c: &dyn Curve2d) -> f64 {
    c.last_parameter()
}

/// `Continuity(C)` (`lxx:47-50`): the `GeomAbs_Shape` value of the adaptor.
pub fn continuity(c: &dyn Curve2d) -> u8 {
    c.continuity()
}

/// `NbIntervals(C, Sh)` (`lxx:52-57`).
pub fn nb_intervals(c: &dyn Curve2d, sh: u8) -> i32 {
    c.nb_intervals(sh)
}

/// `Intervals(C, Tab, Sh)` (`lxx:59-65`): `C->Intervals(Tab, Sh)`.
pub fn intervals(c: &dyn Curve2d, sh: u8) -> Vec<f64> {
    c.parameter_intervals(sh)
}

/// `IsClosed(C)` (`lxx:67-71`).
///
/// `Adaptor2d_Curve2d::IsClosed` is virtual; the default answers
/// `FirstParameter() == LastParameter()` (`Adaptor2d_Curve2d.hxx:62`). The port
/// has no per-class `IsClosed` override on `Curve2d`, so this is the base arm.
pub fn is_closed(c: &dyn Curve2d) -> bool {
    c.first_parameter() == c.last_parameter()
}

/// `IsPeriodic(C)` (`lxx:73-77`).
pub fn is_periodic(c: &dyn Curve2d) -> bool {
    c.is_periodic()
}

/// `Period(C)` (`lxx:79-83`).
pub fn period(c: &dyn Curve2d) -> f64 {
    c.period()
}

/// `Value(C, U)` (`lxx:85-90`).
pub fn value(c: &dyn Curve2d, u: f64) -> occt_core::gp::GpPnt2d {
    c.d0(u)
}

/// `D0(C, U, P)` (`lxx:92-98`).
pub fn d0(c: &dyn Curve2d, u: f64) -> occt_core::gp::GpPnt2d {
    c.d0(u)
}

/// `D1(C, U, P, T)` (`lxx:100-108`).
pub fn d1(c: &dyn Curve2d, u: f64) -> (occt_core::gp::GpPnt2d, occt_core::gp::GpVec2d) {
    c.d1(u)
}

/// `D2(C, U, P, T, N)` (`lxx:110-118`).
pub fn d2(
    c: &dyn Curve2d,
    u: f64,
) -> (
    occt_core::gp::GpPnt2d,
    occt_core::gp::GpVec2d,
    occt_core::gp::GpVec2d,
) {
    c.d2(u)
}

/// `D3(C, U, P, V1, V2, V3)` (`lxx:120-128`).
pub fn d3(
    c: &dyn Curve2d,
    u: f64,
) -> (
    occt_core::gp::GpPnt2d,
    occt_core::gp::GpVec2d,
    occt_core::gp::GpVec2d,
    occt_core::gp::GpVec2d,
) {
    c.d3(u)
}

/// `DN(C, U, N)` (`lxx:130-134`).
pub fn dn(c: &dyn Curve2d, u: f64, n: i32) -> occt_core::gp::GpVec2d {
    c.eval_dn(u, n)
}

/// `Resolution(C, R3d)` (`lxx:136-141`).
pub fn resolution(c: &dyn Curve2d, r3d: f64) -> f64 {
    c.resolution(r3d)
}

/// `GetType(C)` (`lxx:143-147`).
pub fn get_type(c: &dyn Curve2d) -> GeomAbsCurveType {
    curve2d_kind(c)
}

/// `Line(C)` (`lxx:149-155`).
pub fn line(c: &dyn Curve2d) -> Option<occt_core::gp::GpLin2d> {
    c.gp_lin2d()
}

/// `Circle(C)` (`lxx:157-162`).
pub fn circle(c: &dyn Curve2d) -> Option<occt_core::gp::GpCirc2d> {
    c.gp_circ2d()
}

/// `Ellipse(C)` (`lxx:164-169`).
pub fn ellipse(c: &dyn Curve2d) -> Option<occt_core::gp::GpElips2d> {
    c.gp_elips2d()
}

/// `Parabola(C)` (`lxx:171-175`).
pub fn parabola(c: &dyn Curve2d) -> Option<occt_core::gp::GpParab2d> {
    c.gp_parab2d()
}

/// `Hyperbola(C)` (`lxx:177-181`).
pub fn hyperbola(c: &dyn Curve2d) -> Option<occt_core::gp::GpHypr2d> {
    c.gp_hypr2d()
}

/// `NbSamples(C, U0, U1)` (`cxx:20-53`).
pub fn nb_samples(c: &dyn Curve2d, u0: f64, u1: f64) -> i32 {
    let nbs: f64;
    match curve2d_kind(c) {
        // `case GeomAbs_Line: return 2;`
        GeomAbsCurveType::Line => return 2,
        // `case GeomAbs_BezierCurve: nbs = (3 + C->NbPoles());`
        GeomAbsCurveType::BezierCurve => {
            nbs = 3.0 + c.bezier_nb_poles().unwrap_or(0) as f64;
        }
        // `case GeomAbs_BSplineCurve:` (`cxx:34-45`)
        GeomAbsCurveType::BSplineCurve => {
            let mut n = c.bspline_nb_knots().unwrap_or(0) as f64;
            n *= c.bspline_degree().unwrap_or(0) as f64;
            n *= u1 - u0;
            n /= c.last_parameter() - c.first_parameter();
            if n < 2.0 {
                n = 2.0;
            }
            nbs = n;
        }
        // `default: return 25;`
        _ => return 25,
    }
    if nbs > 50.0 {
        // `if (nbs > 50.0) return 50;`
        return 50;
    }
    // `return ((int)nbs);`
    nbs as i32
}
