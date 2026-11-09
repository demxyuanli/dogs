//! `HLRBRep_CurveTool` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_CurveTool.hxx:39-149`,
//! `.lxx:31-257`, `.cxx:22-83`).
//!
//! The `TheCurveTool` of the `IntCurve` instantiation `HLRBRep_CInter`
//! (`HLRBRep_CInter.hxx:26-34` = `IntCurve_IntCurveCurveGen<HLRBRep_CurveTool,
//! Geom2dIntGInter>`). Every function is a one-line forward to
//! `HLRBRep_Curve` (`HLRBRep_Curve.hxx:44-213`), except:
//!
//! * `GetInterval` (`.lxx:77-84`) reads two consecutive entries of the
//!   `Intervals` table,
//! * `Bezier` / `BSpline` (`.lxx:220-238`) are `nullptr` in OCCT itself - the
//!   two bodies are commented out upstream, so the port answers `None`,
//! * `EpsX` (`.lxx:242-245`) is the constant `1e-10`,
//! * `NbSamples` (`.cxx:22-49`) and `NbSamples(C, U0, U1)` (`.cxx:53-83`) are
//!   the only real bodies.
//!
//! `HLRBRep_Curve`'s own `Line` / `Circle` / `Ellipse` accessors answer
//! `Option` in the port (they are derived from the cached projected type,
//! `hlrbrep_curve/projection.rs:198-268`) where OCCT returns the value
//! unconditionally, and `Intervals` answers a `Vec<f64>` where OCCT fills an
//! `NCollection_Array1`. Both differences are already part of the ported
//! `HLRBRep_Curve`.
//!
//! UNPORTED: nothing. `Continuity` / `NbIntervals` / `Intervals` use
//! `GeomAbs_C1` because `IS_C2_COMPOSITE` is `0` (`.lxx:29`, `.lxx:57-72`).

use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d};
use occt_geom::curve::Curve;

use occt_topo::geom_abs_curve_type::GeomAbsCurveType;
use crate::hlrbrep_curve::HlrBRepCurve;

/// `GeomAbs_Shape::GeomAbs_C1` (`GeomAbs_Shape.hxx`): the continuity
/// `NbIntervals` / `Intervals` ask for (`.lxx:57-72`, `IS_C2_COMPOSITE == 0`).
const GEOMABS_C1: u8 = 2;

/// `GeomAbs_Shape::GeomAbs_CN` (`GeomAbs_Shape.hxx`), the continuity
/// `NbSamples(C, U0, U1)` asks the restricted 3D adaptor for (`.cxx:71`).
const GEOMABS_CN: u8 = 7;

/// `HLRBRep_CurveTool::NbSamples`'s `static double nbsOther = 10.0`
/// (`.cxx:26` / `.cxx:58`).
const NBS_OTHER: f64 = 10.0;

/// `HLRBRep_CurveTool::NbSamples`'s `if (nbs > 50) nbs = 50;`
/// (`.cxx:44-47` / `.cxx:78-81`).
const NBS_MAX: f64 = 50.0;

/// `HLRBRep_CurveTool::FirstParameter` (`.lxx:34-37`).
pub fn first_parameter(c: &HlrBRepCurve) -> f64 {
    c.first_parameter()
}

/// `HLRBRep_CurveTool::LastParameter` (`.lxx:42-45`).
pub fn last_parameter(c: &HlrBRepCurve) -> f64 {
    c.last_parameter()
}

/// `HLRBRep_CurveTool::Continuity` (`.lxx:49-52`).
pub fn continuity(c: &HlrBRepCurve) -> u8 {
    c.continuity()
}

/// `HLRBRep_CurveTool::NbIntervals` (`.lxx:56-63`): `C.NbIntervals(GeomAbs_C1)`.
pub fn nb_intervals(c: &HlrBRepCurve) -> i32 {
    c.nb_intervals(GEOMABS_C1)
}

/// `HLRBRep_CurveTool::Intervals` (`.lxx:67-73`): `C.Intervals(Tab, GeomAbs_C1)`.
pub fn intervals(c: &HlrBRepCurve) -> Vec<f64> {
    c.intervals(GEOMABS_C1)
}

/// `HLRBRep_CurveTool::GetInterval` (`.lxx:77-84`): `a = Tab.Value(i)`,
/// `b = Tab.Value(i + 1)` with the 1-based `i` OCCT passes.
pub fn get_interval(tab: &[f64], i: usize) -> (f64, f64) {
    (tab[i - 1], tab[i])
}

/// `HLRBRep_CurveTool::IsClosed` (`.lxx:88-91`).
pub fn is_closed(c: &HlrBRepCurve) -> bool {
    c.is_closed()
}

/// `HLRBRep_CurveTool::IsPeriodic` (`.lxx:95-98`).
pub fn is_periodic(c: &HlrBRepCurve) -> bool {
    c.is_periodic()
}

/// `HLRBRep_CurveTool::Period` (`.lxx:102-105`).
pub fn period(c: &HlrBRepCurve) -> f64 {
    c.period()
}

/// `HLRBRep_CurveTool::Value` (`.lxx:110-113`).
pub fn value(c: &HlrBRepCurve, u: f64) -> GpPnt2d {
    c.value(u)
}

/// `HLRBRep_CurveTool::D0` (`.lxx:117-120`).
pub fn d0(c: &HlrBRepCurve, u: f64) -> GpPnt2d {
    c.d0(u)
}

/// `HLRBRep_CurveTool::D1` (`.lxx:124-131`).
pub fn d1(c: &HlrBRepCurve, u: f64) -> (GpPnt2d, GpVec2d) {
    c.d1(u)
}

/// `HLRBRep_CurveTool::D2` (`.lxx:135-144`).
pub fn d2(c: &HlrBRepCurve, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
    c.d2(u)
}

/// `HLRBRep_CurveTool::D3` (`.lxx:148-155`).
pub fn d3(c: &HlrBRepCurve, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    c.d3(u)
}

/// `HLRBRep_CurveTool::DN` (`.lxx:159-162`).
pub fn dn(c: &HlrBRepCurve, u: f64, n: i32) -> GpVec2d {
    c.dn(u, n)
}

/// `HLRBRep_CurveTool::Resolution` (`.lxx:166-169`).
pub fn resolution(c: &HlrBRepCurve, r3d: f64) -> f64 {
    c.resolution(r3d)
}

/// `HLRBRep_CurveTool::GetType` (`.lxx:173-176`) and `TheType` (`.lxx:180-183`):
/// both answer `C.GetType()`, the cached projected type.
pub fn get_type(c: &HlrBRepCurve) -> GeomAbsCurveType {
    c.get_type()
}

/// `HLRBRep_CurveTool::Line` (`.lxx:186-189`).
pub fn line(c: &HlrBRepCurve) -> Option<GpLin2d> {
    c.line()
}

/// `HLRBRep_CurveTool::Circle` (`.lxx:192-195`).
pub fn circle(c: &HlrBRepCurve) -> Option<GpCirc2d> {
    c.circle()
}

/// `HLRBRep_CurveTool::Ellipse` (`.lxx:198-201`).
pub fn ellipse(c: &HlrBRepCurve) -> Option<GpElips2d> {
    c.ellipse()
}

/// `HLRBRep_CurveTool::Parabola` (`.lxx:204-207`).
pub fn parabola(c: &HlrBRepCurve) -> GpParab2d {
    c.parabola()
}

/// `HLRBRep_CurveTool::Hyperbola` (`.lxx:211-214`).
pub fn hyperbola(c: &HlrBRepCurve) -> GpHypr2d {
    c.hyperbola()
}

/// `HLRBRep_CurveTool::Bezier` (`.lxx:220-228`): OCCT returns `nullptr`
/// unconditionally (the `C->Bezier()` forward is commented out upstream).
pub fn bezier(_c: &HlrBRepCurve) -> Option<()> {
    None
}

/// `HLRBRep_CurveTool::BSpline` (`.lxx:230-238`): `nullptr` as above.
pub fn bspline(_c: &HlrBRepCurve) -> Option<()> {
    None
}

/// `HLRBRep_CurveTool::EpsX` (`.lxx:242-245`): the constant `1e-10`.
pub fn eps_x(_c: &HlrBRepCurve) -> f64 {
    1.0e-10
}

/// `HLRBRep_CurveTool::NbSamples(C)` (`.cxx:22-49`).
pub fn nb_samples(c: &HlrBRepCurve) -> i32 {
    let mut nbs = NBS_OTHER;
    match c.get_type() {
        GeomAbsCurveType::Line => nbs = 2.0,
        GeomAbsCurveType::BezierCurve => nbs = 3.0 + c.nb_poles() as f64,
        GeomAbsCurveType::BSplineCurve => {
            nbs = c.nb_knots() as f64;
            nbs *= c.degree() as f64;
            if nbs < 2.0 {
                nbs = 2.0;
            }
        }
        _ => {}
    }
    if nbs > NBS_MAX {
        nbs = NBS_MAX;
    }
    nbs as i32
}

/// `HLRBRep_CurveTool::NbSamples(C, u1, u2)` (`.cxx:53-83`): the `BSplineCurve`
/// arm builds `GeomAdaptor_Curve(Curve, u1, u2)` and reads
/// `NbIntervals(GeomAbs_CN)`, i.e. `BSplCLib::Intervals` over the restricted
/// range (`GeomAdaptor_Curve.cxx:371-414`). The port reaches the same routine
/// through `occt_core::bspl::adaptor_intervals` with `u1`/`u2` as the range, and
/// falls back to `1` for the `S == GeomAbs_C0` / `S <= Continuity()` short
/// circuit the OCCT adaptor answers before it samples anything.
pub fn nb_samples_range(c: &HlrBRepCurve, u1: f64, u2: f64) -> i32 {
    let mut nbs = NBS_OTHER;
    match c.get_type() {
        GeomAbsCurveType::Line => nbs = 2.0,
        GeomAbsCurveType::BezierCurve => nbs = 3.0 + c.nb_poles() as f64,
        GeomAbsCurveType::BSplineCurve => {
            nbs = restricted_nb_intervals(c, u1, u2) as f64 + 1.0;
            nbs *= c.degree() as f64;
            if nbs < 2.0 {
                nbs = 2.0;
            }
        }
        _ => {}
    }
    if nbs > NBS_MAX {
        nbs = NBS_MAX;
    }
    nbs as i32
}

/// `GeomAdaptor_Curve::NbIntervals(GeomAbs_CN)` (`GeomAdaptor_Curve.cxx:371-414`)
/// over the restricted range `[u1, u2]`.
fn restricted_nb_intervals(c: &HlrBRepCurve, u1: f64, u2: f64) -> i32 {
    let Some(curve) = c.curve() else {
        return 1;
    };
    let curve: &dyn Curve = curve.as_ref();
    // `GeomAdaptor_Curve.cxx:376-379`: a non-periodic B-spline whose continuity
    // already meets the request, or `GeomAbs_C0`, is a single interval.
    if !curve.is_periodic() && GEOMABS_CN <= curve.continuity() {
        return 1;
    }
    let Some(knots) = curve.bspline_knots() else {
        return curve.nb_intervals(GEOMABS_CN);
    };
    let Some(degree) = curve.nurbs_degree() else {
        return curve.nb_intervals(GEOMABS_CN);
    };
    let eps = curve
        .resolution(occt_core::precision::CONFUSION)
        .min(occt_core::precision::PCONFUSION);
    let tab = occt_core::bspl::adaptor_intervals(
        knots,
        degree,
        curve.is_periodic(),
        GEOMABS_CN,
        u1,
        u2,
        eps,
    );
    tab.len().saturating_sub(1).max(1) as i32
}

/// `HLRBRep_CurveTool::Degree` (`.lxx:249-252`).
pub fn degree(c: &HlrBRepCurve) -> i32 {
    c.degree()
}
