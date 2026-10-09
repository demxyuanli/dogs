//! Port of the `IntCurve_IntConicConic::Perform` overloads of
//! `IntCurve_IntConicConic.cxx` whose first curve is a `gp_Circ2d` and whose
//! second curve is not a circle or a line:
//! * `Perform(const gp_Circ2d&, ..., const gp_Parab2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:337-435`);
//! * `Perform(const gp_Circ2d&, ..., const gp_Elips2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:438-483`);
//! * `Perform(const gp_Circ2d&, ..., const gp_Hypr2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:486-581`).
//!
//! All three repeat the closed-domain fix-up of the circle
//! (`SetEquivalentParameters(FirstParameter, FirstParameter + 2*PI)` when
//! `!DC.IsClosed()`); the parabola and hyperbola variants additionally
//! offset-probe the open conic with a radius-shifted copy of the circle.

use std::f64::consts::PI;

use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpParab2d, GpPnt2d, GpVec2d};
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use occt_core::intcurve::{IntCurveIConicTool, IntCurvePConic};
use occt_core::intimpargen::IntImpParGenIntersector;
use occt_core::intres2d::IntRes2dDomain;

use super::int_conic_conic::IntCurveIntConicConic;
use super::int_conic_conic_ana_bounds::{
    infinite_window, probe_domain_tail, set_binf_bsup_from_int_ana2d_hypr,
    set_binf_bsup_from_int_ana2d_parab, DomainOutcome, PARAM_MAX_ON_HYPERBOLA,
    PARAM_MAX_ON_PARABOLA,
};

/// `IntRes2d_Domain D(DC); if (!DC.IsClosed()) D.SetEquivalentParameters(
/// DC.FirstParameter(), DC.FirstParameter() + M_PI + M_PI);`
fn closed_circle_domain(dc: &IntRes2dDomain) -> IntRes2dDomain {
    let mut d = *dc;
    if !dc.is_closed() {
        d.set_equivalent_parameters(dc.first_parameter(), dc.first_parameter() + PI + PI);
    }
    d
}

impl IntCurveIntConicConic {
    /// `Perform(const gp_Circ2d& C, const IntRes2d_Domain& DC,
    /// const gp_Parab2d& P, const IntRes2d_Domain& DP, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:337-435`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_circle_parabola(
        &mut self,
        c: &GpCirc2d,
        dc: &IntRes2dDomain,
        p: &GpParab2d,
        dp: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_circ2d(c);
        let mut p_curve = IntCurvePConic::from_parab2d(p);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        let d = closed_circle_domain(dc);

        let (mut binf, mut bsup) = infinite_window();
        let mut tolinf = 0.0;
        let mut tolsup = 0.0;

        let mut pntinf = GpPnt2d::new(0.0, 0.0);
        let mut pntsup = GpPnt2d::new(0.0, 0.0);

        let maxtol = c.radius / 10.0;
        let mut cp = *c;
        cp.set_radius(c.radius + maxtol);
        let mut the_int_ana2d = IntAna2dAnaIntersection::new();
        the_int_ana2d.perform_parab_conic(p, &IntAna2dConic::from_circ2d(&cp));
        set_binf_bsup_from_int_ana2d_parab(
            &the_int_ana2d,
            &mut binf,
            &mut pntinf,
            &mut bsup,
            &mut pntsup,
            p,
            maxtol,
            PARAM_MAX_ON_PARABOLA,
        );

        if c.radius > maxtol {
            cp.set_radius(c.radius - maxtol);
            the_int_ana2d.perform_parab_conic(p, &IntAna2dConic::from_circ2d(&cp));
            set_binf_bsup_from_int_ana2d_parab(
                &the_int_ana2d,
                &mut binf,
                &mut pntinf,
                &mut bsup,
                &mut pntsup,
                p,
                maxtol,
                PARAM_MAX_ON_PARABOLA,
            );
        }

        match probe_domain_tail(
            dp,
            &mut binf,
            &mut tolinf,
            &mut pntinf,
            &mut bsup,
            &mut tolsup,
            &mut pntsup,
            true,
            false,
        ) {
            DomainOutcome::Run(dp_modif) => {
                inter.perform(&i_tool, &d, &p_curve, &dp_modif, tol_conf, tol);
                self.result.set_values(&inter.base);
            }
            DomainOutcome::NoWindow | DomainOutcome::Stop => {
                self.result.done = true;
            }
        }
    }

    /// `Perform(const gp_Circ2d& C, const IntRes2d_Domain& DC,
    /// const gp_Elips2d& E, const IntRes2d_Domain& DE, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:438-483`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_circle_ellipse(
        &mut self,
        c: &GpCirc2d,
        dc: &IntRes2dDomain,
        e: &GpElips2d,
        de: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_circ2d(c);
        let mut p_curve = IntCurvePConic::from_elips2d(e);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        if !dc.is_closed() {
            let mut d1 = *dc;
            d1.set_equivalent_parameters(dc.first_parameter(), dc.first_parameter() + PI + PI);
            if !de.is_closed() {
                let mut d2 = *de;
                d2.set_equivalent_parameters(de.first_parameter(), de.first_parameter() + PI + PI);
                inter.perform(&i_tool, &d1, &p_curve, &d2, tol_conf, tol);
            } else {
                inter.perform(&i_tool, &d1, &p_curve, de, tol_conf, tol);
            }
        } else if !de.is_closed() {
            let mut d2 = *de;
            d2.set_equivalent_parameters(de.first_parameter(), de.first_parameter() + PI + PI);
            inter.perform(&i_tool, dc, &p_curve, &d2, tol_conf, tol);
        } else {
            inter.perform(&i_tool, dc, &p_curve, de, tol_conf, tol);
        }

        self.result.set_values(&inter.base);
    }

    /// `Perform(const gp_Circ2d& C, const IntRes2d_Domain& DC,
    /// const gp_Hypr2d& H, const IntRes2d_Domain& DH, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:486-581`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_circle_hyperbola(
        &mut self,
        c: &GpCirc2d,
        dc: &IntRes2dDomain,
        h: &GpHypr2d,
        dh: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_circ2d(c);
        let mut p_curve = IntCurvePConic::from_hypr2d(h);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        let d = closed_circle_domain(dc);

        let (mut binf, mut bsup) = infinite_window();
        let mut tolinf = 0.0;
        let mut tolsup = 0.0;

        let mut pntinf = GpPnt2d::new(0.0, 0.0);
        let mut pntsup = GpPnt2d::new(0.0, 0.0);

        let maxtol = c.radius / 10.0;
        let x_axis = h.x_axis();
        let offset = GpVec2d::new(maxtol * x_axis.vdir.x, maxtol * x_axis.vdir.y);
        let mut the_int_ana2d = IntAna2dAnaIntersection::new();

        let hp = h.translated_vec(&offset);
        the_int_ana2d.perform_hypr_conic(&hp, &IntAna2dConic::from_circ2d(c));
        set_binf_bsup_from_int_ana2d_hypr(
            &the_int_ana2d,
            &mut binf,
            &mut pntinf,
            &mut bsup,
            &mut pntsup,
            h,
            maxtol,
            PARAM_MAX_ON_HYPERBOLA,
        );

        let offset = offset.reversed();
        let hm = h.translated_vec(&offset);
        the_int_ana2d.perform_hypr_conic(&hm, &IntAna2dConic::from_circ2d(c));
        set_binf_bsup_from_int_ana2d_hypr(
            &the_int_ana2d,
            &mut binf,
            &mut pntinf,
            &mut bsup,
            &mut pntsup,
            h,
            maxtol,
            PARAM_MAX_ON_HYPERBOLA,
        );

        match probe_domain_tail(
            dh,
            &mut binf,
            &mut tolinf,
            &mut pntinf,
            &mut bsup,
            &mut tolsup,
            &mut pntsup,
            true,
            false,
        ) {
            DomainOutcome::Run(dh_modif) => {
                inter.perform(&i_tool, &d, &p_curve, &dh_modif, tol_conf, tol);
                self.result.set_values(&inter.base);
            }
            DomainOutcome::NoWindow | DomainOutcome::Stop => {
                self.result.done = true;
            }
        }
    }
}
