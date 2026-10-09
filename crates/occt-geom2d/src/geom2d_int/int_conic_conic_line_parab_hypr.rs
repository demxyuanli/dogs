//! Port of the `IntCurve_IntConicConic::Perform` overloads of
//! `IntCurve_IntConicConic.cxx` whose first curve is a `gp_Lin2d` and whose
//! second curve is one of the two open conics:
//! * `Perform(const gp_Lin2d&, ..., const gp_Parab2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:109-225`);
//! * `Perform(const gp_Lin2d&, ..., const gp_Hypr2d&, ...)`
//!   (`IntCurve_IntConicConic.cxx:229-333`).
//!
//! Both offset-probe the open conic to build its parametric window
//! (see `int_conic_conic_ana_bounds.rs`) and then run the generic
//! `IntImpParGen` intersector `Inter`.

use occt_core::gp::{GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d};
use occt_core::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use occt_core::intcurve::{IntCurveIConicTool, IntCurvePConic};
use occt_core::intimpargen::IntImpParGenIntersector;
use occt_core::intres2d::IntRes2dDomain;

use super::int_conic_conic::IntCurveIntConicConic;
use super::int_conic_conic_ana_bounds::{
    bounded_domain, infinite_window, offset_tolerance, offset_tolerance_line_parabola,
    probe_domain_tail, set_binf_bsup_from_int_ana2d_hypr, set_binf_bsup_from_int_ana2d_parab,
    set_bounded_domain, DomainOutcome, PARAM_MAX_ON_HYPERBOLA, PARAM_MAX_ON_PARABOLA,
    TOL_EXACT_INTER,
};

impl IntCurveIntConicConic {
    /// `Perform(const gp_Lin2d& L, const IntRes2d_Domain& DL,
    /// const gp_Parab2d& P, const IntRes2d_Domain& DP, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:109-225`).
    ///
    /// The only overload of the file with the `TOL_EXACT_INTER` retry at
    /// `:205-213` (`wasSet`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_line_parabola(
        &mut self,
        l: &GpLin2d,
        dl: &IntRes2dDomain,
        p: &GpParab2d,
        dp: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_lin2d(l);
        let mut p_curve = IntCurvePConic::from_parab2d(p);
        p_curve.set_accuracy(20);

        let mut inter = IntImpParGenIntersector::new();
        inter.base.set_reversed_parameters(self.result.reversed_parameters());

        let (mut binf, mut bsup) = infinite_window();
        let maxtol = offset_tolerance_line_parabola(tol_conf, tol);
        let mut was_set = false;

        let mut pntinf = GpPnt2d::new(0.0, 0.0);
        let mut pntsup = GpPnt2d::new(0.0, 0.0);
        let mut the_int_ana2d = IntAna2dAnaIntersection::new();

        // `gp_Vec2d Offset(maxtol * L.Direction().Y(), maxtol * L.Direction().X())`
        // (`cxx:138`): the two components are swapped, as in OCCT.
        let offset = GpVec2d::new(maxtol * l.pos.vdir.y, maxtol * l.pos.vdir.x);
        let lp = l.translated_vec(&offset);
        the_int_ana2d.perform_parab_conic(p, &IntAna2dConic::from_lin2d(&lp));
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

        let offset = offset.reversed();
        let lm = l.translated_vec(&offset);
        the_int_ana2d.perform_parab_conic(p, &IntAna2dConic::from_lin2d(&lm));
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

        if binf <= bsup {
            if !bounded_domain(dp) {
                let mut tolinf = 0.0;
                let mut tolsup = 0.0;
                if set_bounded_domain(
                    dp,
                    &mut binf,
                    &mut tolinf,
                    &mut pntinf,
                    &mut bsup,
                    &mut tolsup,
                    &mut pntsup,
                ) {
                    let dp_modif =
                        IntRes2dDomain::bounded(&pntinf, binf, tolinf, &pntsup, bsup, tolsup);
                    inter.perform(&i_tool, dl, &p_curve, &dp_modif, tol_conf, tol);
                } else {
                    self.result.done = true;
                    return;
                }
            } else {
                let mut ft = 0.0;
                let mut lt = 0.0;
                if binf < dp.first_parameter() {
                    binf = dp.first_parameter();
                    pntinf = *dp.first_point();
                    ft = dp.first_tolerance();
                    if bsup < dp.first_parameter() {
                        self.result.done = true;
                        return;
                    }
                }
                if bsup > dp.last_parameter() {
                    bsup = dp.last_parameter();
                    pntsup = *dp.last_point();
                    lt = dp.last_tolerance();
                    if binf > dp.last_parameter() {
                        self.result.done = true;
                        return;
                    }
                }
                let dp_modif = IntRes2dDomain::bounded(&pntinf, binf, ft, &pntsup, bsup, lt);
                inter.perform(&i_tool, dl, &p_curve, &dp_modif, TOL_EXACT_INTER, TOL_EXACT_INTER);
                self.result.set_values(&inter.base);
                was_set = true;
                if self.result.is_done() && self.result.nb_points() == 0 {
                    self.result.reset_fields();
                    inter.perform(&i_tool, dl, &p_curve, &dp_modif, tol_conf, tol);
                    was_set = false;
                }
            }
            if !was_set {
                self.result.set_values(&inter.base);
            }
        } else {
            self.result.done = true;
        }
    }

    /// `Perform(const gp_Lin2d& L, const IntRes2d_Domain& DL,
    /// const gp_Hypr2d& H, const IntRes2d_Domain& DH, TolConf, Tol)`
    /// (`IntCurve_IntConicConic.cxx:229-333`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_line_hyperbola(
        &mut self,
        l: &GpLin2d,
        dl: &IntRes2dDomain,
        h: &GpHypr2d,
        dh: &IntRes2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.result.reset_fields();

        let i_tool = IntCurveIConicTool::from_lin2d(l);
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
        the_int_ana2d.perform_hypr_conic(&hp, &IntAna2dConic::from_lin2d(l));
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
        the_int_ana2d.perform_hypr_conic(&hm, &IntAna2dConic::from_lin2d(l));
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
                inter.perform(&i_tool, dl, &p_curve, &dh_modif, tol_conf, tol);
                self.result.set_values(&inter.base);
            }
            DomainOutcome::NoWindow | DomainOutcome::Stop => {
                self.result.done = true;
            }
        }
    }
}
