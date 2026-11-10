use crate::intcurve::iconic_tool::IntCurveIConicTool;
use crate::intcurve::pconic::IntCurvePConic;
use occt_gp::intres2d::{IntRes2dDomain, IntRes2dIntersection};

/// `IntImpParGen_Intersector` (`IntImpParGen_Intersector.gxx`). The OCCT base
/// `IntRes2d_Intersection` is held in [`Self::base`].
pub struct IntImpParGenIntersector {
    /// The `IntRes2d_Intersection` base.
    pub base: IntRes2dIntersection,
}

impl Default for IntImpParGenIntersector {
    fn default() -> Self {
        Self::new()
    }
}

impl IntImpParGenIntersector {
    /// `IntImpParGen_Intersector()` (`gxx:224-228`).
    pub fn new() -> Self {
        let mut base = IntRes2dIntersection::new();
        base.done = false;
        Self { base }
    }

    /// `IntImpParGen_Intersector(TheImpTool, TheImpCurveDomain, TheParCurve,
    /// TheParCurveDomain, TolConf, Tol)` (`gxx:233-242`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_params(
        imp_tool: &IntCurveIConicTool,
        imp_curve_domain: &IntRes2dDomain,
        par_curve: &IntCurvePConic,
        par_curve_domain: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut s = Self::new();
        s.perform(
            imp_tool,
            imp_curve_domain,
            par_curve,
            par_curve_domain,
            tol_conf,
            tol,
        );
        s
    }
}
