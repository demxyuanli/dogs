use crate::intcurve::iconic_tool::IntCurveIConicTool;
use crate::intcurve::pconic::IntCurvePConic;
use crate::intcurve::pconic_tool;
use occt_gp::math_fn::{MathFunction, MathFunctionWithDerivative};
/// `IntCurve_MyImpParToolOfIntImpConicParConic`
/// (`IntCurve_MyImpParToolOfIntImpConicParConic.hxx:31-57`, `_0.cxx:26-59`):
/// the `math_FunctionWithDerivative` handed to `math_FunctionAllRoots`. It
/// evaluates the signed distance from a point of the parametric curve to the
/// implicit curve.
pub struct MyImpParTool<'a> {
    imp_tool: &'a IntCurveIConicTool,
    par_curve: &'a IntCurvePConic,
}

impl<'a> MyImpParTool<'a> {
    /// `IntCurve_MyImpParToolOfIntImpConicParConic(ITool, PC)` (`_0.cxx:26-31`).
    pub fn new(imp_tool: &'a IntCurveIConicTool, par_curve: &'a IntCurvePConic) -> Self {
        Self {
            imp_tool,
            par_curve,
        }
    }
}

impl MathFunction for MyImpParTool<'_> {
    /// `Value(Param, ApproxDistance)` (`_0.cxx:33-39`).
    fn value(&mut self, param: f64, approx_distance: &mut f64) -> bool {
        *approx_distance = self
            .imp_tool
            .distance(&pconic_tool::pconic_value(self.par_curve, param));
        true
    }
}

impl MathFunctionWithDerivative for MyImpParTool<'_> {
    /// `Derivative(Param, D)` (`_0.cxx:41-50`).
    fn derivative(&mut self, param: f64, d: &mut f64) -> bool {
        let grad = self
            .imp_tool
            .grad_distance(&pconic_tool::pconic_value(self.par_curve, param));
        let (_pt, tan) = pconic_tool::pconic_d1(self.par_curve, param);
        *d = grad.dot(&tan);
        true
    }

    /// `Values(Param, F, D)` (`_0.cxx:52-59`).
    fn values(&mut self, param: f64, approx_distance: &mut f64, deriv: &mut f64) -> bool {
        self.value(param, approx_distance);
        self.derivative(param, deriv);
        true
    }
}
