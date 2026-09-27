//! Port of `Geom2dInt_Geom2dCurveTool`
//! (`Geom2dInt_Geom2dCurveTool.{hxx,lxx,cxx}`): the `TheCurveTool` of the
//! `Geom2dInt_GInter` / `IntCurve` instantiation. It is a thin adaptor over the
//! port's `Curve2d` (= `Adaptor2d_Curve2d`).

use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d};

use crate::curve::Curve2d;
use super::curve_sampling::{nb_samples, nb_samples_range};

/// `GeomAbs_CurveType` (`GeomAbs_CurveType.hxx`), the values
/// `Geom2dAdaptor_Curve::GetType()` can return.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GeomAbsCurveType {
    Line,
    Circle,
    Ellipse,
    Hyperbola,
    Parabola,
    BezierCurve,
    BSplineCurve,
    OffsetCurve,
    OtherCurve,
}

/// `GeomAbs_Shape::GeomAbs_C1` (`GeomAbs_Shape.hxx`).
const GEOMABS_C1: u8 = 2;

/// `Geom2dInt_Geom2dCurveTool::GetType` (`lxx:32-35`) =
/// `Geom2dAdaptor_Curve::GetType()`.
pub fn get_type(c: &dyn Curve2d) -> GeomAbsCurveType {
    if c.is_line() {
        GeomAbsCurveType::Line
    } else if c.gp_circ2d().is_some() {
        GeomAbsCurveType::Circle
    } else if c.gp_elips2d().is_some() {
        GeomAbsCurveType::Ellipse
    } else if c.gp_hypr2d().is_some() {
        GeomAbsCurveType::Hyperbola
    } else if c.gp_parab2d().is_some() {
        GeomAbsCurveType::Parabola
    } else if c.bezier_nb_poles().is_some() {
        GeomAbsCurveType::BezierCurve
    } else if c.bspline_nb_knots().is_some() {
        GeomAbsCurveType::BSplineCurve
    } else if c.offset_basis().is_some() {
        GeomAbsCurveType::OffsetCurve
    } else {
        GeomAbsCurveType::OtherCurve
    }
}

/// `IsComposite` (`hxx:46`). The port carries no composite 2D curves.
pub fn is_composite(_c: &dyn Curve2d) -> bool {
    false
}

/// `Line` (`lxx:38-41`).
pub fn line(c: &dyn Curve2d) -> Option<GpLin2d> {
    c.gp_lin2d()
}

/// `Circle` (`lxx:44-47`).
pub fn circle(c: &dyn Curve2d) -> Option<GpCirc2d> {
    c.gp_circ2d()
}

/// `Ellipse` (`lxx:50-53`).
pub fn ellipse(c: &dyn Curve2d) -> Option<GpElips2d> {
    c.gp_elips2d()
}

/// `Parabola` (`lxx:56-59`).
pub fn parabola(c: &dyn Curve2d) -> Option<GpParab2d> {
    c.gp_parab2d()
}

/// `Hyperbola` (`lxx:62-65`).
pub fn hyperbola(c: &dyn Curve2d) -> Option<GpHypr2d> {
    c.gp_hypr2d()
}

/// `Value` (`lxx:68-71`).
pub fn value(c: &dyn Curve2d, u: f64) -> GpPnt2d {
    c.d0(u)
}

/// `D0` (`lxx:74-77`).
pub fn d0(c: &dyn Curve2d, u: f64) -> GpPnt2d {
    c.d0(u)
}

/// `D1` (`lxx:80-86`).
pub fn d1(c: &dyn Curve2d, u: f64) -> (GpPnt2d, GpVec2d) {
    c.d1(u)
}

/// `D2` (`lxx:89-97`).
pub fn d2(c: &dyn Curve2d, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
    c.d2(u)
}

/// `D3` (`lxx:100-109`).
pub fn d3(c: &dyn Curve2d, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
    c.d3(u)
}

/// `DN` (`lxx:112-117`).
pub fn dn(c: &dyn Curve2d, u: f64, n: i32) -> GpVec2d {
    c.eval_dn(u, n)
}

/// `FirstParameter` (`lxx:120-123`).
pub fn first_parameter(c: &dyn Curve2d) -> f64 {
    c.first_parameter()
}

/// `LastParameter` (`lxx:126-129`).
pub fn last_parameter(c: &dyn Curve2d) -> f64 {
    c.last_parameter()
}

/// `EpsX(C)` (`lxx:134-137`): the mathematical tolerance, constant 1.0e-10.
pub fn eps_x() -> f64 {
    1.0e-10
}

/// `EpsX(C, Eps_XYZ)` = `C.Resolution(Eps_XYZ)` (`lxx:140-143`). UNPORTED:
/// `Curve2d` has no `Resolution`. Not referenced by the `Geom2dInt` / `IntCurve`
/// kernels (a whole-tree grep finds only the one-argument `EpsX`).
pub fn eps_x_with(_c: &dyn Curve2d, _eps_xyz: f64) -> f64 {
    1.0e-10
}

/// `NbSamples(C)` (`cxx:73-91`).
pub fn nb_samples_curve(c: &dyn Curve2d) -> usize {
    nb_samples(c, c.first_parameter(), c.last_parameter())
}

/// `NbSamples(C, U0, U1)` (`cxx:23-70`).
pub fn nb_samples_curve_range(c: &dyn Curve2d, u0: f64, u1: f64) -> usize {
    nb_samples_range(c, u0, u1)
}

/// `NbIntervals` (`lxx:169-178`): `C.NbIntervals(GeomAbs_C1)`
/// (`IS_C2_COMPOSITE` is 0 at `lxx:29`).
pub fn nb_intervals(c: &dyn Curve2d) -> i32 {
    c.nb_intervals(GEOMABS_C1)
}

/// `Intervals` (`lxx:146-154`): `C.Intervals(Tab, GeomAbs_C1)`.
pub fn intervals(c: &dyn Curve2d) -> Vec<f64> {
    c.parameter_intervals(GEOMABS_C1)
}

/// `GetInterval` (`lxx:158-166`): `a = Tab(i)`, `b = Tab(i+1)` (1-based i).
pub fn get_interval(tab: &[f64], i: usize) -> (f64, f64) {
    (tab[i - 1], tab[i])
}

/// `Degree` (`lxx:182-185`) = `Adaptor2d_Curve2d::Degree`. UNPORTED: the
/// underlying `Geom2dAdaptor_Curve::Degree` is not carried by `Curve2d`; a
/// whole-tree grep shows no kernel calls this (`lxx:182` is its only mention).
pub fn degree(_c: &dyn Curve2d) -> i32 {
    1
}
