//! `IntCurve_ProjectOnPConicTool` (TKGeomAlgo,
//! `IntCurve/IntCurve_ProjectOnPConicTool.hxx`, `.cxx`).
//!
//! Projects a point onto the *parametric* `IntCurve_PConic` and returns the
//! corresponding parameter, using the same correspondence as
//! `IntCurve_IConicTool::FindParameter`: the parameter of the point on the
//! implicit conic, wrapped into `[0, 2*PI]` for the closed conics.
use crate::elib::clib2d;
use crate::gp::GpPnt2d;
use crate::kernel::geomabs::CurveType;

use super::pconic::IntCurvePConic;

/// `IntCurve_ProjectOnPConicTool::FindParameter(C, Pnt, Low, High, Tol)`
/// (`cxx:22-90`).
///
/// `LowParameter`/`HighParameter` bound the search interval, and the result is
/// clamped into it. OCCT swaps the two when they arrive inverted (`cxx:30-39`)
/// but only clamps when they differ (`cxx:78`), which is reproduced here. The
/// tolerance is unnamed in OCCT and unused.
pub fn find_parameter_bounded(
    the_p_conic: &IntCurvePConic,
    p: &GpPnt2d,
    low_parameter: f64,
    high_parameter: f64,
    _tol: f64,
) -> f64 {
    let (param_sup, param_inf) = if low_parameter > high_parameter {
        (low_parameter, high_parameter)
    } else {
        (high_parameter, low_parameter)
    };

    let param = raw_parameter(the_p_conic, p);

    if param_inf != param_sup {
        if param < param_inf {
            return param_inf;
        }
        if param > param_sup {
            return param_sup;
        }
    }
    param
}

/// `IntCurve_ProjectOnPConicTool::FindParameter(C, Pnt, Tol)` (`cxx:92-140`):
/// the unbounded overload. The research runs on the natural parametric domain
/// of the curve. The tolerance is unnamed in OCCT and unused.
pub fn find_parameter(the_p_conic: &IntCurvePConic, p: &GpPnt2d, _tol: f64) -> f64 {
    raw_parameter(the_p_conic, p)
}

/// The shared body of both overloads (`cxx:40-75` / `cxx:101-135`).
fn raw_parameter(the_p_conic: &IntCurvePConic, p: &GpPnt2d) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    match the_p_conic.type_curve() {
        CurveType::Line => clib2d::line_parameter_ax2d(&the_p_conic.axis2().x_axis(), p),
        CurveType::Circle => {
            let mut param = clib2d::circle_parameter_ax22d(the_p_conic.axis2(), p);
            if param < 0.0 {
                param += two_pi;
            }
            param
        }
        CurveType::Ellipse => {
            let mut param = clib2d::ellipse_parameter_ax22d(
                the_p_conic.axis2(),
                the_p_conic.param1(),
                the_p_conic.param2(),
                p,
            );
            if param < 0.0 {
                param += two_pi;
            }
            param
        }
        CurveType::Parabola => clib2d::parabola_parameter_ax22d(the_p_conic.axis2(), p),
        CurveType::Hyperbola => clib2d::hyperbola_parameter_ax22d(
            the_p_conic.axis2(),
            the_p_conic.param1(),
            the_p_conic.param2(),
            p,
        ),
        _ => 0.0,
    }
}
