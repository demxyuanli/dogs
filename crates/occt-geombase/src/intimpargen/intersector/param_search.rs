use crate::intcurve::iconic_tool::IntCurveIConicTool;
use crate::intcurve::pconic::IntCurvePConic;
use crate::intcurve::{pconic_tool, project_on_pconic_tool};
use occt_gp::gp::GpPnt2d;
use occt_gp::intres2d::IntRes2dDomain;

use super::super::gen;
use super::structure::IntImpParGenIntersector;

impl IntImpParGenIntersector {
    /// `FindU` (`gxx:781-788`).
    pub(super) fn find_u(
        parameter: f64,
        point: &mut GpPnt2d,
        par_curve: &IntCurvePConic,
        imp_tool: &IntCurveIConicTool,
    ) -> f64 {
        *point = pconic_tool::pconic_value(par_curve, parameter);
        imp_tool.find_parameter(point)
    }

    /// `FindV` (`gxx:790-822`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn find_v(
        parameter: f64,
        point: &mut GpPnt2d,
        imp_tool: &IntCurveIConicTool,
        par_curve: &IntCurvePConic,
        par_curve_domain: &IntRes2dDomain,
        v0: f64,
        v1: f64,
        tolerance: f64,
    ) -> f64 {
        *point = imp_tool.value(parameter);
        if par_curve_domain.is_closed() {
            let v = project_on_pconic_tool::pconic_find_parameter(par_curve, point, tolerance);
            gen::normalize_on_domain(v, par_curve_domain)
        } else {
            let mut vv0 = v0;
            let mut vv1 = v1;
            if v1 < v0 {
                vv0 = v1;
                vv1 = v0;
            }
            let mut x = project_on_pconic_tool::find_parameter_bounded(
                par_curve,
                point,
                vv0,
                vv1,
                tolerance,
            );
            if x > vv1 {
                x = vv1;
            } else if x < vv0 {
                x = vv0;
            }
            x
        }
    }}
