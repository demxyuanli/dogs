//! `IntCurve_PConicTool` (TKGeomAlgo, `IntCurve/IntCurve_PConicTool.hxx`,
//! `.cxx`). Parameter-space evaluation of an [`IntCurvePConic`].
use crate::elib::clib2d;
use crate::gp::{GpPnt2d, GpVec2d};
use crate::kernel::geomabs::CurveType;

use super::pconic::IntCurvePConic;

/// `IntCurve_PConicTool::Value` (`cxx:26-42`). The `default:` arm is
/// `GeomAbs_Hyperbola`, so an unexpected type falls through to the hyperbola
/// evaluation rather than the degenerate hyperbola of a "no match".
pub fn value(c: &IntCurvePConic, x: f64) -> GpPnt2d {
    match c.type_curve() {
        CurveType::Line => clib2d::line_value_ax2d(x, &c.axis2().x_axis()),
        CurveType::Circle => clib2d::circle_value_ax22d(x, c.axis2(), c.param1()),
        CurveType::Ellipse => {
            clib2d::ellipse_value_ax22d(x, c.axis2(), c.param1(), c.param2())
        }
        CurveType::Parabola => clib2d::parabola_value_ax22d(x, c.axis2(), c.param1()),
        _ => clib2d::hyperbola_value_ax22d(x, c.axis2(), c.param1(), c.param2()),
    }
}

/// `IntCurve_PConicTool::D1` (`cxx:48-78`). Unknown types leave `Pt`/`Tan`
/// untouched in OCCT; the Rust port returns the degenerate zero pair.
pub fn d1(c: &IntCurvePConic, u: f64) -> (GpPnt2d, GpVec2d) {
    match c.type_curve() {
        CurveType::Line => clib2d::line_d1_ax2d(u, &c.axis2().x_axis()),
        CurveType::Circle => clib2d::circle_d1_ax22d(u, c.axis2(), c.param1()),
        CurveType::Ellipse => clib2d::ellipse_d1_ax22d(u, c.axis2(), c.param1(), c.param2()),
        CurveType::Parabola => clib2d::parabola_d1_ax22d(u, c.axis2(), c.param1()),
        CurveType::Hyperbola => {
            clib2d::hyperbola_d1_ax22d(u, c.axis2(), c.param1(), c.param2())
        }
        _ => (GpPnt2d::new(0.0, 0.0), GpVec2d::new(0.0, 0.0)),
    }
}

/// `IntCurve_PConicTool::D2` (`cxx:80-114`). The line arm writes the tangent
/// twice — `Tan.SetCoord(0,0)` then `ElCLib::LineD1` overwrites it — so the
/// second derivative is zero, as in [`super::iconic_tool`].
pub fn d2(c: &IntCurvePConic, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
    match c.type_curve() {
        CurveType::Line => {
            let (p, tan) = clib2d::line_d1_ax2d(u, &c.axis2().x_axis());
            (p, tan, GpVec2d::new(0.0, 0.0))
        }
        CurveType::Circle => clib2d::circle_d2_ax22d(u, c.axis2(), c.param1()),
        CurveType::Ellipse => clib2d::ellipse_d2_ax22d(u, c.axis2(), c.param1(), c.param2()),
        CurveType::Parabola => clib2d::parabola_d2_ax22d(u, c.axis2(), c.param1()),
        CurveType::Hyperbola => {
            clib2d::hyperbola_d2_ax22d(u, c.axis2(), c.param1(), c.param2())
        }
        _ => (
            GpPnt2d::new(0.0, 0.0),
            GpVec2d::new(0.0, 0.0),
            GpVec2d::new(0.0, 0.0),
        ),
    }
}

/// `IntCurve_PConicTool::EpsX` (`cxx:116-119`).
pub fn eps_x(c: &IntCurvePConic) -> f64 {
    c.eps_x()
}

/// `IntCurve_PConicTool::NbSamples(const IntCurve_PConic&)` (`cxx:121-124`).
pub fn nb_samples(c: &IntCurvePConic) -> i32 {
    c.accuracy()
}

/// `IntCurve_PConicTool::NbSamples(const IntCurve_PConic&, U0, U1)`
/// (`cxx:127-131`). Both parameters are unnamed in OCCT and ignored here.
pub fn nb_samples_in_range(c: &IntCurvePConic, _u0: f64, _u1: f64) -> i32 {
    c.accuracy()
}
