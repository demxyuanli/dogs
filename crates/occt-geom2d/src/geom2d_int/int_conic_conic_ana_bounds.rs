//! File-static helpers of `IntCurve_IntConicConic.cxx`
//! (`IntCurve_IntConicConic.cxx:48-51 :53-87 :1170-1266`): the two parameter
//! limits `PARAM_MAX_ON_PARABOLA` / `PARAM_MAX_ON_HYPERBOLA`, the exact-tolerance
//! constant `TOL_EXACT_INTER`, `BOUNDED_DOMAIN`, `SET_BOUNDED_DOMAIN` and the two
//! `SetBinfBsupFromIntAna2d` overloads.
//!
//! The `Perform` overloads of that file first probe the open conic with a
//! slightly offset copy of the other element (`IntAna2d_AnaIntersection`) to
//! build a parametric window on the open conic, then hand that window to the
//! generic `IntImpParGen` intersector (`Inter.Perform`). These helpers are the
//! probe-to-window half of that scheme.

use occt_core::elib::clib2d;
use occt_core::gp::{GpHypr2d, GpParab2d, GpPnt2d};
use occt_core::intana2d::IntAna2dAnaIntersection;
use occt_core::intres2d::IntRes2dDomain;

/// `PARAM_MAX_ON_PARABOLA` (`IntCurve_IntConicConic.cxx:48`).
pub(crate) const PARAM_MAX_ON_PARABOLA: f64 = 100000000.0;

/// `PARAM_MAX_ON_HYPERBOLA` (`IntCurve_IntConicConic.cxx:49`).
pub(crate) const PARAM_MAX_ON_HYPERBOLA: f64 = 10000.0;

/// `TOL_EXACT_INTER` (`IntCurve_IntConicConic.cxx:50`).
pub(crate) const TOL_EXACT_INTER: f64 = 1.0e-7;

/// `BOUNDED_DOMAIN` (`IntCurve_IntConicConic.cxx:53-56`).
pub(crate) fn bounded_domain(domain: &IntRes2dDomain) -> bool {
    domain.has_first_point() && domain.has_last_point()
}

/// `SET_BOUNDED_DOMAIN` (`IntCurve_IntConicConic.cxx:58-87`).
///
/// Faithful to the OCCT body, including the `bsup > domain.FirstParameter()`
/// test in the `HasLastPoint()` branch (`:74`, which reads `FirstParameter`
/// where the symmetric branch reads `LastParameter`) and the commented-out early
/// return, so the only result is `bsup > binf`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn set_bounded_domain(
    domain: &IntRes2dDomain,
    binf: &mut f64,
    tolinf: &mut f64,
    pntinf: &mut GpPnt2d,
    bsup: &mut f64,
    tolsup: &mut f64,
    pntsup: &mut GpPnt2d,
) -> bool {
    if domain.has_first_point() && *binf < domain.first_parameter() {
        *pntinf = *domain.first_point();
        *binf = domain.first_parameter();
        *tolinf = domain.first_tolerance();
    }
    if domain.has_last_point() && *bsup > domain.first_parameter() {
        *pntsup = *domain.last_point();
        *bsup = domain.last_parameter();
        *tolsup = domain.last_tolerance();
    }
    *bsup > *binf
}

/// `SetBinfBsupFromIntAna2d(..., const gp_Parab2d&, Maxtol, LIMITE)`
/// (`IntCurve_IntConicConic.cxx:1170-1216`).
pub(crate) fn set_binf_bsup_from_int_ana2d_parab(
    the_int_ana2d: &IntAna2dAnaIntersection,
    binf: &mut f64,
    pntinf: &mut GpPnt2d,
    bsup: &mut f64,
    pntsup: &mut GpPnt2d,
    pr: &GpParab2d,
    maxtol: f64,
    limite: f64,
) {
    if !the_int_ana2d.is_done() || the_int_ana2d.is_empty() {
        return;
    }
    for p in 1..=the_int_ana2d.nb_points() {
        let mut param = the_int_ana2d.point(p).param_on_first();
        if param.abs() >= limite {
            continue;
        }
        let (_pt, v) = clib2d::parabola_d1_ax22d(param, &pr.pos, pr.focal);
        let norme_d1 = v.magnitude();
        let mut dparam = 100.0 * maxtol / norme_d1;
        if dparam < 1e-3 {
            dparam = 1e-3;
        }
        param -= dparam;
        if param < *binf {
            *binf = param;
            *pntinf = clib2d::parabola_value_ax22d(param, &pr.pos, pr.focal);
        }
        param += dparam + dparam;
        if param > *bsup {
            *bsup = param;
            *pntsup = clib2d::parabola_value_ax22d(param, &pr.pos, pr.focal);
        }
    }
}

/// `SetBinfBsupFromIntAna2d(..., const gp_Hypr2d&, Maxtol, LIMITE)`
/// (`IntCurve_IntConicConic.cxx:1218-1266`).
pub(crate) fn set_binf_bsup_from_int_ana2d_hypr(
    the_int_ana2d: &IntAna2dAnaIntersection,
    binf: &mut f64,
    pntinf: &mut GpPnt2d,
    bsup: &mut f64,
    pntsup: &mut GpPnt2d,
    h: &GpHypr2d,
    maxtol: f64,
    limite: f64,
) {
    if !the_int_ana2d.is_done() || the_int_ana2d.is_empty() {
        return;
    }
    for p in 1..=the_int_ana2d.nb_points() {
        let mut param = the_int_ana2d.point(p).param_on_first();
        if param.abs() >= limite {
            continue;
        }
        let (_pt, v) = clib2d::hyperbola_d1_ax22d(param, &h.pos, h.major_radius, h.minor_radius);
        let norme_d1 = v.magnitude();
        let mut dparam = 100.0 * maxtol / norme_d1;
        if dparam < 1e-3 {
            dparam = 1e-3;
        }
        param -= dparam;
        if param < *binf {
            *binf = param;
            *pntinf = clib2d::hyperbola_value_ax22d(param, &h.pos, h.major_radius, h.minor_radius);
        }
        param += dparam + dparam;
        if param > *bsup {
            *bsup = param;
            *pntsup = clib2d::hyperbola_value_ax22d(param, &h.pos, h.major_radius, h.minor_radius);
        }
    }
}

/// The `double binf = Precision::Infinite(), bsup = -Precision::Infinite()`
/// initialisation every `Perform` overload of `IntCurve_IntConicConic.cxx`
/// opens its domain probe with.
pub(crate) fn infinite_window() -> (f64, f64) {
    (occt_core::precision::INFINITE, -occt_core::precision::INFINITE)
}

/// The `maxtol` of `Perform(const gp_Lin2d&, ..., const gp_Parab2d&, ...)`
/// (`IntCurve_IntConicConic.cxx:125-136`): `max(Tol, TolConf)`, floored at
/// `1.e-7`, then scaled by 100. The floor is applied *before* the scaling here,
/// unlike [`offset_tolerance`], which is why the two are not merged.
pub(crate) fn offset_tolerance_line_parabola(tol_conf: f64, tol: f64) -> f64 {
    let maxtol = if tol > tol_conf { tol } else { tol_conf };
    let maxtol = if maxtol < 1.0e-7 { 1.0e-7 } else { maxtol };
    maxtol * 100.0
}

/// The `maxtol` shared by the `Perform` overloads that offset a hyperbola
/// (`IntCurve_IntConicConic.cxx:246-257`, `:602-611`, `:824-835`, `:1085-1096`):
/// `max(Tol, TolConf)`, scaled by 100, then floored at `0.000001`.
pub(crate) fn offset_tolerance(tol_conf: f64, tol: f64) -> f64 {
    let maxtol = if tol > tol_conf { tol } else { tol_conf };
    let maxtol = maxtol * 100.0;
    if maxtol < 0.000001 {
        0.000001
    } else {
        maxtol
    }
}

/// Outcome of the `if (binf <= bsup) { ... } else { done = true; }` tail shared
/// by every `Perform` overload that probes an open conic.
pub(crate) enum DomainOutcome {
    /// `binf > bsup`: the caller sets `done = true` and returns.
    NoWindow,
    /// An inner early return of the tail (`SET_BOUNDED_DOMAIN` false, or one of
    /// the `binf`/`bsup` window tests): the caller sets `done = true` and
    /// returns *without* `SetValues`.
    Stop,
    /// `Inter.Perform(ITool, D1, PCurve, domain, TolConf, Tol)` must be run,
    /// after which the caller runs `SetValues(Inter)`.
    Run(IntRes2dDomain),
}

/// The domain-window tail shared by the `Perform` overloads of
/// `IntCurve_IntConicConic.cxx` that bound an open conic
/// (`:162-221`, `:283-326`, `:385-428`, `:531-575`, `:636-683`, `:752-801`,
/// `:861-905`, `:1006-1058`, `:1118-1162`). The two `Perform` overloads that end
/// in `:1049-1053` and `:1152-1156` add an extra `if (binf >= bsup)` test
/// (`reversed_window_check`), and the `Hypr2d`/`Hypr2d` overload at `:1139-1150`
/// is the only one whose `binf`/`bsup` clamps have no nested early return
/// (`nested_early_out = false`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn probe_domain_tail(
    domain: &IntRes2dDomain,
    binf: &mut f64,
    tolinf: &mut f64,
    pntinf: &mut GpPnt2d,
    bsup: &mut f64,
    tolsup: &mut f64,
    pntsup: &mut GpPnt2d,
    nested_early_out: bool,
    reversed_window_check: bool,
) -> DomainOutcome {
    if *binf > *bsup {
        return DomainOutcome::NoWindow;
    }
    if !bounded_domain(domain) {
        if set_bounded_domain(domain, binf, tolinf, pntinf, bsup, tolsup, pntsup) {
            DomainOutcome::Run(IntRes2dDomain::bounded(
                pntinf, *binf, *tolinf, pntsup, *bsup, *tolsup,
            ))
        } else {
            DomainOutcome::Stop
        }
    } else {
        if *binf < domain.first_parameter() {
            *binf = domain.first_parameter();
            *pntinf = *domain.first_point();
            *tolinf = domain.first_tolerance();
            if nested_early_out && *bsup < domain.first_parameter() {
                return DomainOutcome::Stop;
            }
        }
        if *bsup > domain.last_parameter() {
            *bsup = domain.last_parameter();
            *pntsup = *domain.last_point();
            *tolsup = domain.last_tolerance();
            if nested_early_out && *binf > domain.last_parameter() {
                return DomainOutcome::Stop;
            }
        }
        if reversed_window_check && *binf >= *bsup {
            return DomainOutcome::Stop;
        }
        DomainOutcome::Run(IntRes2dDomain::bounded(
            pntinf, *binf, *tolinf, pntsup, *bsup, *tolsup,
        ))
    }
}
