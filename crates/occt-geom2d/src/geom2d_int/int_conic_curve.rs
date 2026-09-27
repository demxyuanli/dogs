//! Port of `IntCurve_IntConicCurveGen` / `Geom2dInt_TheIntConicCurveOfGInter`
//! (`IntCurve_IntConicCurveGen.gxx:24-92`, `.lxx:30-130`; instantiations at
//! `Geom2dInt_TheIntConicCurveOfGInter_0.cxx:32-52`).

use std::f64::consts::PI;

use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d};
use occt_core::intcurve::IntCurveIConicTool;
use occt_core::intres2d::{IntRes2dDomain, IntRes2dIntersection};

use crate::curve::Curve2d;
use super::conic_curve::Geom2dIntConicCurve;

/// `IntCurve_IntConicCurveGen` (`IntCurve_IntConicCurveGen.hxx`), the
/// `Geom2dInt` instantiation (`Geom2dInt_TheIntConicCurveOfGInter`).
#[derive(Clone, Debug, Default)]
pub struct IntConicCurveGen {
    base: IntRes2dIntersection,
}

impl IntConicCurveGen {
    pub fn new() -> Self {
        Self { base: IntRes2dIntersection::new() }
    }

    /// `IntRes2d_Intersection::SetReversedParameters` (inherited).
    pub fn set_reversed_parameters(&mut self, flag: bool) {
        self.base.set_reversed_parameters(flag);
    }

    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    pub fn result(&self) -> &IntRes2dIntersection {
        &self.base
    }

    /// `Perform(ICurve, D1, PCurve, D2, TolConf, Tol)` (`lxx:119-130`).
    pub fn perform_iconic_tool(
        &mut self,
        i_curve: &IntCurveIConicTool,
        d1: &IntRes2dDomain,
        p_curve: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let mut my_intersection = Geom2dIntConicCurve::new();
        my_intersection
            .base
            .set_reversed_parameters(self.base.reversed_parameters());
        my_intersection.perform(i_curve, d1, p_curve, d2, tol_conf, tol);
        self.base.set_values(&my_intersection.base);
    }

    /// `Perform(const gp_Lin2d&, ...)` (`lxx:44-53`).
    pub fn perform_line(
        &mut self,
        l: &GpLin2d,
        d1: &IntRes2dDomain,
        p_curve: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let tool = IntCurveIConicTool::from_lin2d(l);
        self.perform_iconic_tool(&tool, d1, p_curve, d2, tol_conf, tol);
    }

    /// `Perform(const gp_Circ2d&, ...)` (`lxx:56-74`).
    pub fn perform_circle(
        &mut self,
        c: &GpCirc2d,
        d1: &IntRes2dDomain,
        p_curve: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let tool = IntCurveIConicTool::from_circ2d(c);
        if !d1.is_closed() {
            let mut d = *d1;
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI + PI);
            self.perform_iconic_tool(&tool, &d, p_curve, d2, tol_conf, tol);
        } else {
            self.perform_iconic_tool(&tool, d1, p_curve, d2, tol_conf, tol);
        }
    }

    /// `Perform(const gp_Elips2d&, ...)` (`lxx:77-94`).
    pub fn perform_ellipse(
        &mut self,
        e: &GpElips2d,
        d1: &IntRes2dDomain,
        p_curve: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let tool = IntCurveIConicTool::from_elips2d(e);
        if !d1.is_closed() {
            let mut d = *d1;
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI + PI);
            self.perform_iconic_tool(&tool, &d, p_curve, d2, tol_conf, tol);
        } else {
            self.perform_iconic_tool(&tool, d1, p_curve, d2, tol_conf, tol);
        }
    }

    /// `Perform(const gp_Parab2d&, ...)` (`lxx:97-105`).
    pub fn perform_parabola(
        &mut self,
        p: &GpParab2d,
        d1: &IntRes2dDomain,
        p_curve: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let tool = IntCurveIConicTool::from_parab2d(p);
        self.perform_iconic_tool(&tool, d1, p_curve, d2, tol_conf, tol);
    }

    /// `Perform(const gp_Hypr2d&, ...)` (`lxx:108-116`).
    pub fn perform_hyperbola(
        &mut self,
        h: &GpHypr2d,
        d1: &IntRes2dDomain,
        p_curve: &dyn Curve2d,
        d2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let tool = IntCurveIConicTool::from_hypr2d(h);
        self.perform_iconic_tool(&tool, d1, p_curve, d2, tol_conf, tol);
    }
}
