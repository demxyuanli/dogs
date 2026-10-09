//! Port of the `IntCurve_IntConicConic::Perform` overloads of
//! `IntCurve_IntConicConic.cxx` whose two curves are both closed or both open
//! insofar as neither is a `gp_Lin2d`:
//! * `Perform(const gp_Elips2d&, ..., const gp_Elips2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:913-957`);
//! * `Perform(const gp_Elips2d&, ..., const gp_Hypr2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:959-1062`);
//! * `Perform(const gp_Hypr2d&, ..., const gp_Hypr2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:1064-1166`).
//!
//! The two hyperbola-carrying overloads install the window with the extra
//! `if (binf >= bsup)` early return (`:1049-1053`, `:1152-1156`), and the
//! `Hypr2d`/`Hypr2d` one is the only `Perform` of the file whose `binf`/`bsup`
//! clamps have no nested early return (`:1139-1150`).

use std::f64::consts::PI;

use occt_core::gp::{GpElips2d, GpHypr2d, GpPnt2d, GpVec2d};
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use occt_core::intcurve::{IntCurveIConicTool, IntCurvePConic};
use occt_core::intimpargen::IntImpParGenIntersector;
use occt_core::intres2d::IntRes2dDomain;

use super::int_conic_conic::IntCurveIntConicConic;
use super::int_conic_conic_ana_bounds::{
    infinite_window, offset_tolerance, probe_domain_tail, set_binf_bsup_from_int_ana2d_hypr,
    DomainOutcome, PARAM_MAX_ON_HYPERBOLA,
};

impl IntCurveIntConicConic {
    /// `Perform(const gp_Elips2d& E1, const IntRes2d_Domain& DE1,
    /// const gp_Elips2d& E2, const IntRes2d_Domain& DE2, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:913-957`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_ellipse_ellipse(
        &mut self,
        e1: &GpElips2d,
        de1: &IntRes2dDomain,
        e2: &GpElips2d,
        de2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_elips2d(e1);
        let mut p_curve = IntCurvePConic::from_elips2d(e2);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        if !de1.is_closed() {
            let mut d1 = *de1;
            d1.set_equivalent_parameters(de1.first_parameter(), de1.first_parameter() + PI + PI);
            if !de2.is_closed() {
                let mut d2 = *de2;
                d2.set_equivalent_parameters(de2.first_parameter(), de2.first_parameter() + PI + PI);
                inter.perform(&i_tool, &d1, &p_curve, &d2, tol_conf, tol);
            } else {
                inter.perform(&i_tool, &d1, &p_curve, de2, tol_conf, tol);
            }
        } else if !de2.is_closed() {
            let mut d2 = *de2;
            d2.set_equivalent_parameters(de2.first_parameter(), de2.first_parameter() + PI + PI);
            inter.perform(&i_tool, de1, &p_curve, &d2, tol_conf, tol);
        } else {
            inter.perform(&i_tool, de1, &p_curve, de2, tol_conf, tol);
        }

        self.result.set_values(&inter.base);
    }

    /// `Perform(const gp_Elips2d& E, const IntRes2d_Domain& DE,
    /// const gp_Hypr2d& H, const IntRes2d_Domain& DH, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:959-1062`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_ellipse_hyperbola(
        &mut self,
        e: &GpElips2d,
        de: &IntRes2dDomain,
        h: &GpHypr2d,
        dh: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_elips2d(e);
        let mut p_curve = IntCurvePConic::from_hypr2d(h);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        let mut de_modif = *de;
        if !de.is_closed() {
            de_modif.set_equivalent_parameters(de.first_parameter(), de.first_parameter() + PI + PI);
        }

        let (mut binf, mut bsup) = infinite_window();
        let mut tolinf = 0.0;
        let mut tolsup = 0.0;

        let mut pntinf = GpPnt2d::new(0.0, 0.0);
        let mut pntsup = GpPnt2d::new(0.0, 0.0);

        let maxtol = e.minor_radius / 10.0;
        let x_axis = h.x_axis();
        let offset = GpVec2d::new(maxtol * x_axis.vdir.x, maxtol * x_axis.vdir.y);
        let mut the_int_ana2d = IntAna2dAnaIntersection::new();

        let hp = h.translated_vec(&offset);
        the_int_ana2d.perform_hypr_conic(&hp, &IntAna2dConic::from_elips2d(e));
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
        the_int_ana2d.perform_hypr_conic(&hm, &IntAna2dConic::from_elips2d(e));
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
            true,
        ) {
            DomainOutcome::Run(dh_modif) => {
                inter.perform(&i_tool, &de_modif, &p_curve, &dh_modif, tol_conf, tol);
                self.result.set_values(&inter.base);
            }
            DomainOutcome::NoWindow | DomainOutcome::Stop => {
                self.result.done = true;
            }
        }
    }

    /// `Perform(const gp_Hypr2d& H1, const IntRes2d_Domain& DH1,
    /// const gp_Hypr2d& H2, const IntRes2d_Domain& DH2, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:1064-1166`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_hyperbola_hyperbola(
        &mut self,
        h1: &GpHypr2d,
        dh1: &IntRes2dDomain,
        h2: &GpHypr2d,
        dh2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_hypr2d(h1);
        let mut p_curve = IntCurvePConic::from_hypr2d(h2);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        let (mut binf, mut bsup) = infinite_window();
        let mut tolinf = 0.0;
        let mut tolsup = 0.0;

        let mut pntinf = GpPnt2d::new(0.0, 0.0);
        let mut pntsup = GpPnt2d::new(0.0, 0.0);

        let maxtol = offset_tolerance(tol_conf, tol);
        let x_axis = h2.x_axis();
        let offset = GpVec2d::new(maxtol * x_axis.vdir.x, maxtol * x_axis.vdir.y);
        let mut the_int_ana2d = IntAna2dAnaIntersection::new();

        let hp = h2.translated_vec(&offset);
        the_int_ana2d.perform_hypr_conic(&hp, &IntAna2dConic::from_hypr2d(h1));
        set_binf_bsup_from_int_ana2d_hypr(
            &the_int_ana2d,
            &mut binf,
            &mut pntinf,
            &mut bsup,
            &mut pntsup,
            h2,
            maxtol,
            PARAM_MAX_ON_HYPERBOLA,
        );

        let offset = offset.reversed();
        let hm = h2.translated_vec(&offset);
        the_int_ana2d.perform_hypr_conic(&hm, &IntAna2dConic::from_hypr2d(h1));
        set_binf_bsup_from_int_ana2d_hypr(
            &the_int_ana2d,
            &mut binf,
            &mut pntinf,
            &mut bsup,
            &mut pntsup,
            h2,
            maxtol,
            PARAM_MAX_ON_HYPERBOLA,
        );

        match probe_domain_tail(
            dh2,
            &mut binf,
            &mut tolinf,
            &mut pntinf,
            &mut bsup,
            &mut tolsup,
            &mut pntsup,
            false,
            true,
        ) {
            DomainOutcome::Run(dh_modif) => {
                inter.perform(&i_tool, dh1, &p_curve, &dh_modif, tol_conf, tol);
                self.result.set_values(&inter.base);
            }
            DomainOutcome::NoWindow | DomainOutcome::Stop => {
                self.result.done = true;
            }
        }
    }
}
