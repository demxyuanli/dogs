//! `HLRBRep_LineTool` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_LineTool.hxx:47-179`,
//! `.lxx:29-229`).
//!
//! "The LineTool class provides class methods to access the methods of the
//! Line" (`hxx:46-47`). Every member takes the `gp_Lin` itself, so the port is a
//! set of free functions over `GpLin`, the shape the `.lxx` gives them.
//!
//! The tool is what `HLRBRep_InterCSurf` and `HLRBRep_Intersector` use for the
//! curve side of their `gp_Lin` instantiations (`HLRBRep_InterCSurf.cxx:44-45`:
//! `using TheCurveTool = HLRBRep_LineTool`), i.e. the curve half of
//! `IntCurveSurface_Inter`'s `CurveTool` parameter. Two of its answers feed the
//! sampling:
//!
//! * `NbSamples(C, U0, U1)` is a constant `3` (`lxx:201-204`), where
//!   `IntCurveSurface_TheHCurveTool::NbSamples` (`...cxx:31-70`) answers `2` for
//!   a line;
//! * `SamplePars(C, U0, U1, Defl, NbMin)` answers `{U0, (U0+U1)/2, U1}`
//!   (`lxx:206-218`), where the `IntCurveSurface` tool answers `{U0, U1}` for a
//!   line. Both matter because the point count of the curve polygon changes the
//!   section points the interference engine reports.
//!
//! The `Circle` / `Ellipse` / `Hyperbola` / `Parabola` / `Bezier` / `BSpline`
//! accessors answer a null/empty object and `GetType` answers `GeomAbs_Line`
//! unconditionally (`lxx:124-167`), so the four non-line arms of
//! `IntCurveSurface_InterImpl::PerformBounds` (`IntCurveSurface_Inter.pxx:126-136`)
//! are unreachable in this instantiation.

use occt_core::elib::clib;
use occt_core::gp::{GpLin, GpPnt, GpVec};

/// `FirstParameter` (`lxx:29-32`): `RealFirst()`.
pub fn first_parameter(_c: &GpLin) -> f64 {
    -f64::MAX
}

/// `LastParameter` (`lxx:32-35`): `RealLast()`.
pub fn last_parameter(_c: &GpLin) -> f64 {
    f64::MAX
}

/// `Continuity` (`lxx:38-41`): `GeomAbs_CN`.
pub fn continuity(_c: &GpLin) -> u8 {
    6
}

/// `NbIntervals` (`lxx:43-46`): one interval.
pub fn nb_intervals(_c: &GpLin, _s: u8) -> usize {
    1
}

/// `Intervals` (`lxx:48-53`): fills nothing.
pub fn intervals(_c: &GpLin, _t: &mut Vec<f64>, _sh: u8) {}

/// `IntervalFirst` (`lxx:55-58`).
pub fn interval_first(_c: &GpLin) -> f64 {
    -f64::MAX
}

/// `IntervalLast` (`lxx:60-63`).
pub fn interval_last(_c: &GpLin) -> f64 {
    f64::MAX
}

/// `IntervalContinuity` (`lxx:65-68`): `GeomAbs_CN`.
pub fn interval_continuity(_c: &GpLin) -> u8 {
    6
}

/// `IsClosed` (`lxx:70-73`).
pub fn is_closed(_c: &GpLin) -> bool {
    false
}

/// `IsPeriodic` (`lxx:75-78`).
pub fn is_periodic(_c: &GpLin) -> bool {
    false
}

/// `Period` (`lxx:80-83`): `0.`.
pub fn period(_c: &GpLin) -> f64 {
    0.0
}

/// `Value` (`lxx:85-88`) = `ElCLib::Value`.
pub fn value(c: &GpLin, u: f64) -> GpPnt {
    clib::line_value(c, u)
}

/// `D0` (`lxx:90-93`).
pub fn d0(c: &GpLin, u: f64) -> GpPnt {
    clib::line_value(c, u)
}

/// `D1` (`lxx:95-98`) = `ElCLib::D1`.
pub fn d1(c: &GpLin, u: f64) -> (GpPnt, GpVec) {
    clib::line_d1(c, u)
}

/// `D2` (`lxx:100-105`): `ElCLib::D1` and a zero second derivative.
pub fn d2(c: &GpLin, u: f64) -> (GpPnt, GpVec, GpVec) {
    let (p, v1) = clib::line_d1(c, u);
    (p, v1, GpVec::new(0.0, 0.0, 0.0))
}

/// `D3` (`lxx:107-115`): `ElCLib::D1` and zero higher derivatives.
pub fn d3(c: &GpLin, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
    let (p, v1) = clib::line_d1(c, u);
    let z = GpVec::new(0.0, 0.0, 0.0);
    (p, v1, z, GpVec::new(0.0, 0.0, 0.0))
}

/// `DN` (`lxx:117-120`) = `ElCLib::DN`.
pub fn dn(c: &GpLin, _u: f64, n: i32) -> GpVec {
    clib::line_dn(c, n)
}

/// `Resolution` (`lxx:122-125`): a line's parameter is arc length.
pub fn resolution(_c: &GpLin, r3d: f64) -> f64 {
    r3d
}

/// `GetType` (`lxx:127-130`): always `GeomAbs_Line`.
pub fn get_type(_c: &GpLin) -> CurveType {
    CurveType::Line
}

/// `GeomAbs_CurveType` (`GeomAbs_CurveType.hxx:19-30`), restricted to the
/// members `HLRBRep_LineTool` can answer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveType {
    Line,
}

/// `Line` (`lxx:132-135`): the curve itself.
pub fn line(c: &GpLin) -> GpLin {
    c.clone()
}

/// `Degree` (`lxx:169-172`).
pub fn degree(_c: &GpLin) -> i32 {
    1
}

/// `NbPoles` (`lxx:174-177`): `0`.
pub fn nb_poles(_c: &GpLin) -> i32 {
    0
}

/// `IsRational` (`lxx:181-184`).
pub fn is_rational(_c: &GpLin) -> bool {
    false
}

/// `NbKnots` (`lxx:191-194`): `0`.
pub fn nb_knots(_c: &GpLin) -> i32 {
    0
}

/// `NbSamples(C, U0, U1)` (`lxx:201-204`): `3`.
pub fn nb_samples(_c: &GpLin, _u0: f64, _u1: f64) -> usize {
    3
}

/// `SamplePars(C, U0, U1, Defl, NbMin)` (`lxx:206-218`): three parameters,
/// `{U0, (U0+U1)/2, U1}`, ignoring `Defl` and `NbMin` for a line.
pub fn sample_pars(_c: &GpLin, u0: f64, u1: f64, _defl: f64, _nb_min: usize) -> Vec<f64> {
    vec![u0, (u0 + u1) * 0.5, u1]
}
