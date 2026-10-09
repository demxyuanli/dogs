//! Port of the `IntCurve_IntConicConic::Perform` overloads of
//! `IntCurve_IntConicConic.cxx` in which a `gp_Parab2d` is the first curve:
//! * `Perform(const gp_Parab2d&, ..., const gp_Parab2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:584-688`);
//! * `Perform(const gp_Elips2d&, ..., const gp_Parab2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:690-807`);
//! * `Perform(const gp_Parab2d&, ..., const gp_Hypr2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:809-911`).

use std::f64::consts::PI;

use occt_core::gp::{GpElips2d, GpHypr2d, GpParab2d, GpPnt2d, GpVec2d};
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use occt_core::intcurve::{IntCurveIConicTool, IntCurvePConic};
use occt_core::intimpargen::IntImpParGenIntersector;
use occt_core::intres2d::IntRes2dDomain;

use super::int_conic_conic::IntCurveIntConicConic;
use super::int_conic_conic_ana_bounds::{
    infinite_window, offset_tolerance, probe_domain_tail,
    set_binf_bsup_from_int_ana2d_hypr, set_binf_bsup_from_int_ana2d_parab, DomainOutcome,
    PARAM_MAX_ON_HYPERBOLA, PARAM_MAX_ON_PARABOLA,
};

impl IntCurveIntConicConic {
    /// `Perform(const gp_Parab2d& P1, const IntRes2d_Domain& DP1,
    /// const gp_Parab2d& P2, const IntRes2d_Domain& DP2, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:584-688`).
    ///
    /// The probe offsets the *second* parabola (`P2`) along its `MirrorAxis` and
    /// intersects it with `P1`, but the window is measured on `P2`.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_parabola_parabola(
        &mut self,
        p1: &GpParab2d,
        dp1: &IntRes2dDomain,
        p2: &GpParab2d,
        dp2: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_parab2d(p1);
        let mut p_curve = IntCurvePConic::from_parab2d(p2);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        let (mut binf, mut bsup) = infinite_window();
        let mut tolinf = 0.0;
        let mut tolsup = 0.0;

        let mut pntinf = GpPnt2d::new(0.0, 0.0);
        let mut pntsup = GpPnt2d::new(0.0, 0.0);

        let maxtol = offset_tolerance(tol_conf, tol);
        let mirror_axis = p2.mirror_axis();
        let offset = GpVec2d::new(maxtol * mirror_axis.vdir.x, maxtol * mirror_axis.vdir.y);
        let mut the_int_ana2d = IntAna2dAnaIntersection::new();

        let pp = p2.translated_vec(&offset);
        the_int_ana2d.perform_parab_conic(&pp, &IntAna2dConic::from_parab2d(p1));
        set_binf_bsup_from_int_ana2d_parab(
            &the_int_ana2d,
            &mut binf,
            &mut pntinf,
            &mut bsup,
            &mut pntsup,
            p2,
            maxtol,
            PARAM_MAX_ON_PARABOLA,
        );

        let offset = offset.reversed();
        let pm = p2.translated_vec(&offset);
        the_int_ana2d.perform_parab_conic(&pm, &IntAna2dConic::from_parab2d(p1));
        set_binf_bsup_from_int_ana2d_parab(
            &the_int_ana2d,
            &mut binf,
            &mut pntinf,
            &mut bsup,
            &mut pntsup,
            p2,
            maxtol,
            PARAM_MAX_ON_PARABOLA,
        );

        match probe_domain_tail(
            dp2,
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
                inter.perform(&i_tool, dp1, &p_curve, &dp_modif, tol_conf, tol);
                self.result.set_values(&inter.base);
            }
            DomainOutcome::NoWindow | DomainOutcome::Stop => {
                self.result.done = true;
            }
        }
    }

    /// `Perform(const gp_Elips2d& E, const IntRes2d_Domain& DE,
    /// const gp_Parab2d& P, const IntRes2d_Domain& DP, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:690-807`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_ellipse_parabola(
        &mut self,
        e: &GpElips2d,
        de: &IntRes2dDomain,
        p: &GpParab2d,
        dp: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_elips2d(e);
        let mut p_curve = IntCurvePConic::from_parab2d(p);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        let mut d = *de;
        if !de.is_closed() {
            d.set_equivalent_parameters(de.first_parameter(), de.first_parameter() + PI + PI);
        }

        let (mut binf, mut bsup) = infinite_window();
        let mut tolinf = 0.0;
        let mut tolsup = 0.0;

        let mut pntinf = GpPnt2d::new(0.0, 0.0);
        let mut pntsup = GpPnt2d::new(0.0, 0.0);

        let maxtol = e.minor_radius / 10.0;
        let mut ep = *e;
        ep.set_major_radius(e.major_radius + maxtol);
        ep.set_minor_radius(e.minor_radius + maxtol);
        let mut the_int_ana2d = IntAna2dAnaIntersection::new();
        the_int_ana2d.perform_parab_conic(p, &IntAna2dConic::from_elips2d(&ep));
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

        if e.minor_radius > maxtol {
            ep.set_minor_radius(e.minor_radius - maxtol);
            ep.set_major_radius(e.major_radius - maxtol);
            the_int_ana2d.perform_parab_conic(p, &IntAna2dConic::from_elips2d(&ep));
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

    /// `Perform(const gp_Parab2d& P, const IntRes2d_Domain& DP,
    /// const gp_Hypr2d& H, const IntRes2d_Domain& DH, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:809-911`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_parabola_hyperbola(
        &mut self,
        p: &GpParab2d,
        dp: &IntRes2dDomain,
        h: &GpHypr2d,
        dh: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_parab2d(p);
        let mut p_curve = IntCurvePConic::from_hypr2d(h);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        let (mut binf, mut bsup) = infinite_window();
        let mut tolinf = 0.0;
        let mut tolsup = 0.0;

        let mut pntinf = GpPnt2d::new(0.0, 0.0);
        let mut pntsup = GpPnt2d::new(0.0, 0.0);

        let maxtol = offset_tolerance(tol_conf, tol);
        let x_axis = h.x_axis();
        let offset = GpVec2d::new(maxtol * x_axis.vdir.x, maxtol * x_axis.vdir.y);
        let mut the_int_ana2d = IntAna2dAnaIntersection::new();

        let hp = h.translated_vec(&offset);
        the_int_ana2d.perform_hypr_conic(&hp, &IntAna2dConic::from_parab2d(p));
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
        the_int_ana2d.perform_hypr_conic(&hm, &IntAna2dConic::from_parab2d(p));
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
                inter.perform(&i_tool, dp, &p_curve, &dh_modif, tol_conf, tol);
                self.result.set_values(&inter.base);
            }
            DomainOutcome::NoWindow | DomainOutcome::Stop => {
                self.result.done = true;
            }
        }
    }
}
