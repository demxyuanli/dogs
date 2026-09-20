use super::prelude::*;
use super::*;

/// All local extrema of |C1−C2| via grid → bracket sign changes → Newton.
pub(super) fn newton_curve_curve_all(c1: &dyn Curve, c2: &dyn Curve) -> Vec<ExtremaPair> {
    let (a1, b1) = (c1.first_parameter(), c1.last_parameter());
    let (a2, b2) = (c2.first_parameter(), c2.last_parameter());
    let us = build_samples(c1, a1, b1);
    let vs = build_samples(c2, a2, b2);
    let n1 = us.len();
    let n2 = vs.len();

    // Sample F1, F2 and the squared distance on the grid.
    let mut f1g = vec![vec![f64::NAN; n2]; n1];
    let mut f2g = vec![vec![f64::NAN; n2]; n1];
    let mut dg = vec![vec![f64::NAN; n2]; n1];
    for i in 0..n1 {
        for j in 0..n2 {
            let u = us[i];
            let v = vs[j];
            let p1 = c1.d0(u);
            let p2 = c2.d0(v);
            let (_, du) = c1.d1(u);
            let (_, dv) = c2.d1(v);
            if !p1.x().is_finite() || !p2.x().is_finite() {
                continue;
            }
            let d = GpVec::from_pnts(&p2, &p1);
            let ndu = du.magnitude();
            let ndv = dv.magnitude();
            if ndu > 1e-12 {
                f1g[i][j] = d.dot(&du) / ndu;
            }
            if ndv > 1e-12 {
                f2g[i][j] = d.dot(&dv) / ndv;
            }
            dg[i][j] = p1.square_distance(&p2);
        }
    }

    let mut seeds: Vec<(f64, f64)> = Vec::new();

    // Sign changes of F1 (along u) and F2 (along v) bracket a root.
    for i in 0..n1.saturating_sub(1) {
        for j in 0..n2.saturating_sub(1) {
            let a = f1g[i][j];
            let b = f1g[i + 1][j];
            let c = f2g[i][j];
            let d = f2g[i][j + 1];
            if !(a.is_finite() && b.is_finite() && c.is_finite() && d.is_finite()) {
                continue;
            }
            let s1 = (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0) || a.abs() < 1e-14 || b.abs() < 1e-14;
            let s2 = (c < 0.0 && d > 0.0) || (c > 0.0 && d < 0.0) || c.abs() < 1e-14 || d.abs() < 1e-14;
            if s1 && s2 {
                seeds.push((0.5 * (us[i] + us[i + 1]), 0.5 * (vs[j] + vs[j + 1])));
            }
        }
    }

    // Local min/max suppression on the sampled squared distance (catches
    // tangency extrema missed by sign changes).
    for i in 1..n1.saturating_sub(1) {
        for j in 1..n2.saturating_sub(1) {
            let d = dg[i][j];
            if !d.is_finite() {
                continue;
            }
            let neighbors = [dg[i - 1][j], dg[i + 1][j], dg[i][j - 1], dg[i][j + 1]];
            if !neighbors.iter().all(|x| x.is_finite()) {
                continue;
            }
            let is_min = d <= neighbors[0] && d <= neighbors[1] && d <= neighbors[2] && d <= neighbors[3];
            let is_max = d >= neighbors[0] && d >= neighbors[1] && d >= neighbors[2] && d >= neighbors[3];
            if is_min || is_max {
                seeds.push((us[i], vs[j]));
            }
        }
    }

    // Boundary seeds (best grid point on each edge) for trimmed curves.
    if n1 >= 2 && n2 >= 2 {
        let mut edges: [Option<(f64, f64)>; 4] = [None, None, None, None];
        for i in 0..n1 {
            if dg[i][0].is_finite() && edges[0].map_or(true, |(_, d)| dg[i][0] < d) {
                edges[0] = Some((us[i], vs[0]));
            }
            if dg[i][n2 - 1].is_finite() && edges[1].map_or(true, |(_, d)| dg[i][n2 - 1] < d) {
                edges[1] = Some((us[i], vs[n2 - 1]));
            }
        }
        for j in 0..n2 {
            if dg[0][j].is_finite() && edges[2].map_or(true, |(_, d)| dg[0][j] < d) {
                edges[2] = Some((us[0], vs[j]));
            }
            if dg[n1 - 1][j].is_finite() && edges[3].map_or(true, |(_, d)| dg[n1 - 1][j] < d) {
                edges[3] = Some((us[n1 - 1], vs[j]));
            }
        }
        for e in edges.into_iter().flatten() {
            seeds.push(e);
        }
    }

    let mut out: Vec<ExtremaPair> = Vec::new();
    for (u0, v0) in seeds {
        let (u, v) = refine_curve_curve(c1, c2, u0, v0, a1, b1, a2, b2);
        let p1 = c1.d0(u);
        let p2 = c2.d0(v);
        if p1.x().is_finite() && p2.x().is_finite() {
            out.push(pair_cc(p1, u, p2, v));
        }
    }
    dedupe_sort(out)
}

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

/// All local extrema (minima and maxima) of the distance between two curves,
/// deduplicated and sorted by distance.
///
/// Lines and circles are classified and solved analytically (exact); every
/// other pair goes through the grid + Newton path. Unbounded curves use
/// expanding-window grids.
pub fn curve_curve_extrema_all(c1: &dyn Curve, c2: &dyn Curve) -> Vec<ExtremaPair> {
    let (a1, b1) = (c1.first_parameter(), c1.last_parameter());
    let (a2, b2) = (c2.first_parameter(), c2.last_parameter());

    let l1 = is_line(c1).then(|| line_of_curve(c1)).flatten();
    let l2 = is_line(c2).then(|| line_of_curve(c2)).flatten();
    let g1 = classify_circle(c1);
    let g2 = classify_circle(c2);

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

    dedupe_sort(newton_curve_curve_all(c1, c2))
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
