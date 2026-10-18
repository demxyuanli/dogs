//! Port of `Extrema_CurveTool` (`Extrema_CurveTool.hxx:38-142`): the
//! `TheCurveTool` of the `Adaptor3d_Curve` instantiation
//! (`Extrema_PCLocFOfLocEPCOfLocateExtPC.hxx:29-35`). A thin adaptor over the
//! port's `Curve` (= `Adaptor3d_Curve`).
//!
//! The 2D twin is `occt_geom2d::geom2d_int::curve_tool`
//! (`= Geom2dInt_Geom2dCurveTool`).

use occt_core::gp::{GpPnt, GpVec};

use crate::curve::Curve;

/// `GeomAbs_CurveType` (`GeomAbs_CurveType.hxx`), the values
/// `GeomAdaptor_Curve::GetType()` can return.
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

/// `GeomAdaptor_Curve::GetType()` for a port curve.
///
/// The order is `GeomAdaptor_Curve::load` (`GeomAdaptor_Curve.cxx:252-311`):
/// trimmed curves forward every query to their basis, then Circle, Line,
/// Ellipse, Parabola, Hyperbola, Bezier, BSpline, Offset, Other. Only the
/// concrete class decides -- a geometrically straight B-spline is
/// `GeomAbs_BSplineCurve` in OCCT, so `Curve::is_line` (a geometric test) is
/// deliberately **not** used here.
pub fn extrema_curve_type(c: &dyn Curve) -> GeomAbsCurveType {
    if c.gp_circ().is_some() {
        GeomAbsCurveType::Circle
    } else if c.gp_line().is_some() {
        GeomAbsCurveType::Line
    } else if c.gp_ellipse().is_some() {
        GeomAbsCurveType::Ellipse
    } else if c.gp_parabola().is_some() {
        GeomAbsCurveType::Parabola
    } else if c.gp_hyperbola().is_some() {
        GeomAbsCurveType::Hyperbola
    } else if c.bezier_poles().is_some() {
        GeomAbsCurveType::BezierCurve
    } else if c.bspline_knots().is_some() || c.nurbs_degree().is_some() {
        GeomAbsCurveType::BSplineCurve
    } else if c.offset_curve().is_some() {
        GeomAbsCurveType::OffsetCurve
    } else {
        GeomAbsCurveType::OtherCurve
    }
}

/// `D0` (`hxx:84-87`).
pub fn extrema_d0(c: &dyn Curve, u: f64) -> GpPnt {
    c.d0(u)
}

/// `D1` (`hxx:89-93`).
pub fn extrema_d1(c: &dyn Curve, u: f64) -> (GpPnt, GpVec) {
    c.d1(u)
}

/// `D2` (`hxx:95-104`).
pub fn extrema_d2(c: &dyn Curve, u: f64) -> (GpPnt, GpVec, GpVec) {
    c.d2(u)
}

/// `DN` (`hxx:112-115`).
pub fn extrema_dn(c: &dyn Curve, u: f64, n: i32) -> GpVec {
    c.eval_dn(u, n)
}

/// `FirstParameter` (`hxx:42`).
pub fn extrema_first_parameter(c: &dyn Curve) -> f64 {
    c.first_parameter()
}

/// `LastParameter` (`hxx:44`).
pub fn extrema_last_parameter(c: &dyn Curve) -> f64 {
    c.last_parameter()
}

/// `Resolution` (`hxx:74-77`).
pub fn resolution(c: &dyn Curve, r3d: f64) -> f64 {
    c.resolution(r3d)
}

/// `NbIntervals` (`hxx:50-53`).
pub fn extrema_nb_intervals(c: &dyn Curve, s: u8) -> i32 {
    c.nb_intervals(s)
}

/// `Intervals` (`hxx:56-63`).
pub fn extrema_intervals(c: &dyn Curve, s: u8) -> Vec<f64> {
    c.parameter_intervals(s)
}
