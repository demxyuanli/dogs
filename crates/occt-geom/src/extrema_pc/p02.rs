use super::prelude::*;
use super::*;


/// Deduplicate by point proximity and sort ascending by distance.
pub(super) fn dedupe_sort(v: Vec<ExtremaPair>) -> Vec<ExtremaPair> {
    let mut out: Vec<ExtremaPair> = Vec::new();
    for e in v {
        if let Some(o) = out.iter_mut().find(|o| o.p2.distance(&e.p2) < 1e-6) {
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
/// Refine the original-curve parameter that hits the analytic point `q`,
/// seeded near the affine-mapped circle parameter.
pub(super) fn param_for_point(c: &dyn Curve, q: &GpPnt, seed: f64, a: f64, b: f64) -> f64 {
    let mut u = seed.clamp(a, b);
    for _ in 0..16 {
        let (cq, d1, d2) = c.d2(u);
        let v = GpVec::from_pnts(q, &cq);
        let fp = d1.square_magnitude() + v.dot(&d2);
        if fp.abs() < 1e-300 {
            break;
        }
        let un = u - v.dot(&d1) / fp;
        if !un.is_finite() {
            break;
        }
        if (un - u).abs() < 1e-12 * (1.0 + u.abs()) {
            u = un;
            break;
        }
        u = un.clamp(a, b);
    }
    u
}

// ---------------------------------------------------------------------------
// Public dispatch.
// ---------------------------------------------------------------------------

/// All local extrema (minima and maxima) of the point–curve distance,
/// deduplicated and sorted by distance.
///
/// Lines and circles are classified and solved analytically (exact); every
/// other curve goes through the grid + Newton path.
pub fn point_curve_extrema_all(c: &dyn Curve, p: &GpPnt) -> Vec<ExtremaPair> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let mut pairs: Vec<ExtremaPair> = Vec::new();

    if is_line(c) {
        // Analytic line. d0(u) = C0 + u·T0 with constant tangent T0, so the
        // projection parameter is u* = (P-C0)·T0 / |T0|².
        let u0 = if a.is_finite() { a } else { 0.0 };
        let c0 = c.d0(u0);
        let t0 = c.d1(u0).1;
        let l2 = t0.square_magnitude();
        if l2 > 0.0 {
            let up = ((p.x() - c0.x()) * t0.x() + (p.y() - c0.y()) * t0.y() + (p.z() - c0.z()) * t0.z()) / l2;
            if up.is_finite() && up >= a - CONFUSION && up <= b + CONFUSION {
                pairs.push(pair(c, p, up));
            }
        }
        if a.is_finite() {
            pairs.push(pair(c, p, a));
        }
        if b.is_finite() {
            pairs.push(pair(c, p, b));
        }
        return dedupe_sort(pairs);
    }

    if let Some(gc) = classify_circle(c) {
        // Analytic circle solutions within the curve's range, then map each
        // back to the original parameterization (Newton on the original curve,
        // verified against the analytic point) so parameter offsets and curve
        // orientation cannot introduce a mismatch.
        let period = 2.0 * std::f64::consts::PI;
        for s in circle_all(&gc, p, a, b) {
            let seed_a = s.u1;
            let seed_b = 2.0 * a - s.u1;
            let mut best = s.u1;
            let mut best_err = f64::INFINITY;
            for seed in [seed_a, seed_b] {
                let ws = a + (seed - a).rem_euclid(period);
                if ws < a - CONFUSION || ws > b + CONFUSION {
                    continue;
                }
                let u = param_for_point(c, &s.p2, ws, a, b);
                let err = c.d0(u).distance(&s.p2);
                if err < best_err {
                    best_err = err;
                    best = u;
                }
            }
            pairs.push(pair(c, p, best));
        }
        if a.is_finite() {
            pairs.push(pair(c, p, a));
        }
        if b.is_finite() {
            pairs.push(pair(c, p, b));
        }
        return dedupe_sort(pairs);
    }

    dedupe_sort(newton_point_curve_all(c, p))
}

/// Minimum distance from `p` to `c` (with the closest point and parameter).
pub fn point_curve_extrema(c: &dyn Curve, p: &GpPnt) -> ExtremaPair {
    match point_curve_extrema_all(c, p).into_iter().next() {
        Some(e) => e,
        None => {
            // Degenerate curve: fall back to the golden-section refine.
            let (u, q) =
                crate::extrema::refine_curve_point(c, p, c.first_parameter(), c.last_parameter());
            ExtremaPair { p1: *p, p2: q, distance: p.distance(&q), u1: u, v1: None, u2: u, v2: None }
        }
    }
}

/// Maximum distance from `p` to `c` (farthest local extremum).
pub fn point_curve_max_extrema(c: &dyn Curve, p: &GpPnt) -> ExtremaPair {
    match point_curve_extrema_all(c, p).into_iter().last() {
        Some(e) => e,
        None => point_curve_extrema(c, p),
    }
}
