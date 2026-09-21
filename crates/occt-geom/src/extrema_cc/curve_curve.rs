use super::prelude::*;
use super::*;

/// All local extrema of |C1−C2| via grid → bracket sign changes → Newton.

/// Deduplicate by parameter proximity and sort ascending by distance.
pub(super) fn dedupe_sort(v: Vec<ExtremaPair>) -> Vec<ExtremaPair> {
    let mut out: Vec<ExtremaPair> = Vec::new();
    for e in v {
        if let Some(o) = out.iter_mut().find(|o| (o.u1 - e.u1).abs() < 1e-5 && (o.u2 - e.u2).abs() < 1e-5) {
            if e.distance < o.distance {
                *o = e;
            }
        } else {
            out.push(e);
        }
    }
    out.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
    out
}

pub(super) fn param_in(u: f64, lo: f64, hi: f64) -> bool {
    let lo_ok = !lo.is_finite() || u >= lo - CONFUSION;
    let hi_ok = !hi.is_finite() || u <= hi + CONFUSION;
    lo_ok && hi_ok
}

/// Map `u` into `[lo, lo+period)` for a periodic curve (OCCT `ElCLib::InPeriod`).
pub(super) fn in_period(u: f64, lo: f64, period: f64) -> f64 {
    if period <= 0.0 || !lo.is_finite() {
        u
    } else {
        lo + (u - lo).rem_euclid(period)
    }
}

// ---------------------------------------------------------------------------
// Public dispatch.
// ---------------------------------------------------------------------------

/// `Extrema_ExtCC` over the curves' own parameter ranges.
pub fn curve_curve_extrema_all(c1: &dyn Curve, c2: &dyn Curve) -> Vec<ExtremaPair> {
    let (a1, b1) = (c1.first_parameter(), c1.last_parameter());
    let (a2, b2) = (c2.first_parameter(), c2.last_parameter());
    curve_curve_extrema_all_range(c1, c2, a1, b1, a2, b2)
}

/// `Extrema_ExtCC(C1, C2, U1, U2, V1, V2)` — the extrema over explicit parameter
/// ranges (`Extrema_ExtCC.cxx:150-190` installs them with `Initialize`/`SetParams`
/// and hands them to `Extrema_ECC`, i.e. `Extrema_GGenExtCC`, as
/// `myLowBorder`/`myUppBorder`, `cxx:180`).
///
/// This is what a caller that owns a *bounded* piece of an unbounded curve (an
/// edge whose curve is an infinite line) must use: feeding the adaptor's own
/// `±Precision::Infinite()` bounds to the optimizer cannot resolve extrema near
/// small parameters.
///
/// Lines/circles plus elementary partners are solved analytically
/// (`Extrema_ExtElC`, `Extrema_ExtCC.cxx:251-305`); everything else goes through
/// [`GGenExtCC`]. Deduplicated and sorted by distance.
pub fn curve_curve_extrema_all_range(
    c1: &dyn Curve,
    c2: &dyn Curve,
    a1: f64,
    b1: f64,
    a2: f64,
    b2: f64,
) -> Vec<ExtremaPair> {

    // `Extrema_ExtCC::Perform` dispatches on `Adaptor3d_Curve::GetType()`
    // (`Extrema_ExtCC.cxx:251-305`): a line plus an elementary curve, or two
    // circles, go to the analytic `Extrema_ExtElC`; everything else to
    // `Extrema_ECC` (= `Extrema_GGenExtCC`). The port reads the type from the
    // curve's own queries (`gp_line`/`gp_circ`), not from a sampling
    // classifier — the previous `is_line`/`classify_circle` reconstruction is
    // gone with the invented seed set.
    let l1 = c1.gp_line();
    let l2 = c2.gp_line();
    let g1 = c1.gp_circ();
    let g2 = c2.gp_circ();

    let analytic: Vec<ExtremaPair> = match (l1, l2, g1, g2) {
        (Some(l1), Some(l2), _, _) => line_line_extrema(&l1, &l2),
        (Some(l), None, _, Some(c)) => line_circle_extrema(&l, &c),
        (None, Some(l), Some(c), _) => {
            // Swap: line is c2, circle is c1.
            let mut out = line_circle_extrema(&l, &c);
            for e in out.iter_mut() {
                std::mem::swap(&mut e.p1, &mut e.p2);
                std::mem::swap(&mut e.u1, &mut e.u2);
                std::mem::swap(&mut e.v1, &mut e.v2);
            }
            out
        }
        (None, None, Some(g1), Some(g2)) => circle_circle_extrema(&g1, &g2),
        _ => Vec::new(),
    };

    if !analytic.is_empty() {
        // Map periodic parameters into the curve range before the range check.
        let p1 = if c1.is_periodic() { c1.period() } else { 0.0 };
        let p2 = if c2.is_periodic() { c2.period() } else { 0.0 };
        let analytic: Vec<ExtremaPair> = analytic
            .into_iter()
            .map(|mut e| {
                e.u1 = in_period(e.u1, a1, p1);
                e.u2 = in_period(e.u2, a2, p2);
                e
            })
            .collect();
        let filtered: Vec<ExtremaPair> = analytic
            .into_iter()
            .filter(|e| param_in(e.u1, a1, b1) && param_in(e.u2, a2, b2))
            .collect();
        if !filtered.is_empty() {
            return dedupe_sort(filtered);
        }
    }

    // `Extrema_ExtCC` falls back to `Extrema_ECC`, which *is*
    // `Extrema_GGenExtCC` (`Extrema_ECC.hxx:23-28`), for every pair that is not
    // an elementary `Extrema_ExtElC` case (`Extrema_ExtCC.cxx:248-305`).
    let mut g = GGenExtCC::new(c1, c2, (a1, a2), (b1, b2));
    if g.perform().is_err() || !g.is_done() {
        return Vec::new();
    }
    let mut out: Vec<ExtremaPair> = Vec::new();
    for n in 1..=g.nb_ext() {
        let (u, p1, v, p2) = g.points(n);
        out.push(pair_cc(p1, u, p2, v));
    }
    dedupe_sort(out)
}

/// Minimum distance between two curves (with the closest points and
/// parameters). The first entry of `curve_curve_extrema_all` (sorted ascending).
///
/// `None` means OCCT's `myDone = false` — `Extrema_GGenExtCC::Perform`
/// (`Extrema_GGenExtCC.hxx:691-695`) returns without solutions when the search
/// finds none, and no fabricated pair is produced. The previous body fell back to
/// a 16×16 parameter grid and finally to a zero-distance pair at the origin
/// (audit A15).
pub fn curve_curve_extrema(c1: &dyn Curve, c2: &dyn Curve) -> Option<ExtremaPair> {
    curve_curve_extrema_all(c1, c2).into_iter().next()
}

/// Local extremum of the distance between two curves from seed parameters
/// (Newton refinement). Port of `Extrema_LocateExtCC` / `Extrema_LocECC`.
pub fn locate_extcc(c1: &dyn Curve, c2: &dyn Curve, u0: f64, v0: f64) -> Option<ExtremaPair> {
    let (a1, b1) = (c1.first_parameter(), c1.last_parameter());
    let (a2, b2) = (c2.first_parameter(), c2.last_parameter());
    let (u, v) = refine_curve_curve(c1, c2, u0, v0, a1, b1, a2, b2);
    let p1 = c1.d0(u);
    let p2 = c2.d0(v);
    if p1.x().is_finite() && p2.x().is_finite() {
        Some(pair_cc(p1, u, p2, v))
    } else {
        None
    }
}
