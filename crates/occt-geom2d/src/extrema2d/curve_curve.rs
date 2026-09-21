use super::prelude::*;
use super::*;

pub(super) fn param_for_point2d(c: &dyn Curve2d, q: &GpPnt2d, seed: f64, a: f64, b: f64) -> f64 {
    let mut u = seed.clamp(a, b);
    for _ in 0..16 {
        let (cq, d1, d2) = c.d2(u);
        let v = GpVec2d::new(cq.x() - q.x(), cq.y() - q.y());
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

/// All local extrema (minima and maxima) of the point–curve distance in the
/// plane, deduplicated and sorted by distance. Lines and circles are
/// classified and solved analytically; everything else uses grid + Newton.
pub fn point_curve_extrema2d_all(c: &dyn Curve2d, p: &GpPnt2d) -> Vec<Extrema2d> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let mut pairs: Vec<Extrema2d> = Vec::new();

    if is_line2d(c) {
        let u0 = if a.is_finite() { a } else { 0.0 };
        let c0 = c.d0(u0);
        let t0 = c.d1(u0).1;
        let l2 = t0.square_magnitude();
        if l2 > 0.0 {
            let up = ((p.x() - c0.x()) * t0.x() + (p.y() - c0.y()) * t0.y()) / l2;
            if up.is_finite() && up >= a - CONFUSION && up <= b + CONFUSION {
                pairs.push(pair2d_c(c, p, up));
            }
        }
        if a.is_finite() {
            pairs.push(pair2d_c(c, p, a));
        }
        if b.is_finite() {
            pairs.push(pair2d_c(c, p, b));
        }
        return dedupe_sort2d(pairs);
    }

    if let Some(gc) = classify_circle2d(c) {
        let period = 2.0 * std::f64::consts::PI;
        for s in circle_all2d(&gc, p, a, b) {
            let mut best = s.u1;
            let mut best_err = f64::INFINITY;
            for seed in [s.u1, 2.0 * a - s.u1] {
                let ws = a + (seed - a).rem_euclid(period);
                if ws < a - CONFUSION || ws > b + CONFUSION {
                    continue;
                }
                let u = param_for_point2d(c, &s.p2, ws, a, b);
                let err = c.d0(u).distance(&s.p2);
                if err < best_err {
                    best_err = err;
                    best = u;
                }
            }
            pairs.push(pair2d_c(c, p, best));
        }
        if a.is_finite() {
            pairs.push(pair2d_c(c, p, a));
        }
        if b.is_finite() {
            pairs.push(pair2d_c(c, p, b));
        }
        return dedupe_sort2d(pairs);
    }

    dedupe_sort2d(newton_point_curve_all2d(c, p))
}

/// Minimum distance from point `p` to curve `c` (with closest point).
/// Exact analytic dispatch for lines/circles; Newton refinement otherwise
/// (see `point_curve_extrema2d_all`).
pub fn point_curve_extrema2d(c: &dyn Curve2d, p: &GpPnt2d) -> Extrema2d {
    match point_curve_extrema2d_all(c, p).into_iter().next() {
        Some(e) => e,
        None => {
            // Degenerate curve: golden-section fallback.
            let f = |u: f64| {
                let q = c.d0(u);
                (q.x() - p.x()).powi(2) + (q.y() - p.y()).powi(2)
            };
            let (a0, b0) = (c.first_parameter(), c.last_parameter());
            let (a, b) = if a0.is_finite() && b0.is_finite() && b0 > a0 {
                (a0, b0)
            } else {
                let mut w = 1.0;
                let mut fm = f(0.0);
                for _ in 0..10 {
                    if f(-w) >= fm && f(w) >= fm {
                        break;
                    }
                    w *= 8.0;
                    fm = f(0.0);
                }
                (-w, w)
            };
            let phi = (5.0f64.sqrt() - 1.0) * 0.5;
            let (mut lo, mut hi) = (a, b);
            for _ in 0..80 {
                let x1 = hi - phi * (hi - lo);
                let x2 = lo + phi * (hi - lo);
                if f(x1) < f(x2) {
                    hi = x2;
                } else {
                    lo = x1;
                }
            }
            let u = 0.5 * (lo + hi);
            let q = c.d0(u);
            Extrema2d { p1: *p, p2: q, distance: p.distance(&q), u1: u, u2: u }
        }
    }
}

/// Maximum distance from point `p` to curve `c` (sampling-based).
pub fn point_curve_max_extrema2d(c: &dyn Curve2d, p: &GpPnt2d, samples: usize) -> Extrema2d {
    let (a, b) = bound(c);
    let n = samples.max(2);
    let mut best: Option<Extrema2d> = None;
    for i in 0..=n {
        let u = a + (b - a) * i as f64 / n as f64;
        let q = c.d0(u);
        let d = p.distance(&q);
        if best.as_ref().map_or(true, |e: &Extrema2d| d > e.distance) {
            best = Some(Extrema2d { p1: *p, p2: q, distance: d, u1: u, u2: u });
        }
    }
    best.unwrap()
}

/// Minimum distance between two curves (all extrema, first = closest).
///
/// Delegates to `curve_curve_extrema2d_all` (exact analytic for classified
/// lines/circles, grid + Newton otherwise). `samples` is kept for signature
/// compatibility but no longer drives the computation.
pub fn curve_curve_extrema2d(c1: &dyn Curve2d, c2: &dyn Curve2d, _samples: usize) -> Vec<Extrema2d> {
    curve_curve_extrema2d_all(c1, c2)
}

// ---------------------------------------------------------------------------
// Curve–curve extrema (port of `Extrema_ExtCC2d`, `Extrema_ExtElC2d`).
// ---------------------------------------------------------------------------

pub(super) fn pair2d_cc(p1: GpPnt2d, u1: f64, p2: GpPnt2d, u2: f64) -> Extrema2d {
    Extrema2d { p1, p2, distance: p1.distance(&p2), u1, u2 }
}

pub(super) fn line2d_value(l: &GpLin2d, u: f64) -> GpPnt2d {
    GpPnt2d::new(l.pos.loc.x() + u * l.pos.vdir.x, l.pos.loc.y() + u * l.pos.vdir.y)
}

pub(super) fn circle_param2d(c: &GpCirc2d, p: &GpPnt2d) -> f64 {
    let v = GpVec2d::new(p.x() - c.location().x(), p.y() - c.location().y());
    let xd = c.pos.vxdir;
    let yd = c.pos.vydir;
    // C(u) = O + r·cos·xd + r·sin·yd → u = atan2(v·yd, v·xd).
    v.dot(&GpVec2d::new(yd.x, yd.y)).atan2(v.dot(&GpVec2d::new(xd.x, xd.y)))
}

/// Exact extrema between two 2D lines (port of `Extrema_ExtElC2d`).
/// Parallel lines have constant distance; otherwise they intersect (distance 0).
pub fn line_line_extrema2d(l1: &GpLin2d, l2: &GpLin2d) -> Vec<Extrema2d> {
    let d1 = GpVec2d::from_dir2d(&l1.pos.vdir);
    let d2 = GpVec2d::from_dir2d(&l2.pos.vdir);
    if d1.is_parallel(&d2, ANGULAR) {
        let p1 = l1.pos.loc;
        let u2 = (p1.x() - l2.pos.loc.x()) * d2.x() + (p1.y() - l2.pos.loc.y()) * d2.y();
        let p2 = line2d_value(l2, u2);
        return vec![pair2d_cc(p1, 0.0, p2, u2)];
    }
    // Cramer's rule for the intersection.
    let a_p1p2 = GpVec2d::new(l2.pos.loc.x() - l1.pos.loc.x(), l2.pos.loc.y() - l1.pos.loc.y());
    let delim = 1.0 / d1.crossed(&d2);
    let param1 = a_p1p2.crossed(&d2) * delim;
    let param2 = -(d1.crossed(&a_p1p2)) * delim;
    let p1 = line2d_value(l1, param1);
    let p2 = line2d_value(l2, param2);
    vec![pair2d_cc(p1, param1, p2, param2)]
}

/// Exact extrema between a 2D line and a 2D circle (port of
/// `Extrema_ExtElC2d(const gp_Lin2d&, const gp_Circ2d&)`), plus the
/// zero-distance intersections.
pub fn line_circle_extrema2d(l: &GpLin2d, c: &GpCirc2d) -> Vec<Extrema2d> {
    let d = l.pos.vdir;
    let x2 = c.pos.vxdir;
    let y2 = c.pos.vydir;
    let dx = d.dot(&x2);
    let dy = d.dot(&y2);
    let o1 = l.pos.loc;

    // Radius perpendicular to the line direction.
    let mut tetas: Vec<f64> = Vec::new();
    if dy.abs() <= 1e-15 {
        tetas.push(PI / 2.0);
    } else {
        tetas.push((-dx / dy).atan());
    }
    tetas.push(tetas[0] + PI);
    if tetas[0] < 0.0 {
        tetas[0] += 2.0 * PI;
    }

    let mut out = Vec::new();
    for teta in tetas {
        let p2 = circle2d_val(c, teta);
        let u1 = (p2.x() - o1.x()) * d.x + (p2.y() - o1.y()) * d.y;
        let p1 = line2d_value(l, u1);
        out.push(pair2d_cc(p1, u1, p2, teta));
    }

    // Intersections: |O1 + t·D − O2|² = r².
    let o2 = c.location();
    let r = c.radius();
    let v12 = GpVec2d::new(o1.x() - o2.x(), o1.y() - o2.y());
    let a = d.x * d.x + d.y * d.y;
    let b = 2.0 * (v12.x() * d.x + v12.y() * d.y);
    let c0 = v12.square_magnitude() - r * r;
    let disc = b * b - 4.0 * a * c0;
    if disc >= 0.0 && a > 1e-300 {
        let sq = disc.sqrt();
        for t in [(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
            let p1 = line2d_value(l, t);
            let u2 = circle_param2d(c, &p1);
            out.push(pair2d_cc(p1, t, p1, u2));
        }
    }
    dedupe_sort_cc2d(out)
}

/// Exact extrema between two 2D circles (port of
/// `Extrema_ExtElC2d(const gp_Circ2d&, const gp_Circ2d&)`), plus the
/// zero-distance intersections for crossing circles.
pub fn circle_circle_extrema2d(c1: &GpCirc2d, c2: &GpCirc2d) -> Vec<Extrema2d> {
    let o1 = c1.location();
    let o2 = c2.location();
    let r1 = c1.radius();
    let r2 = c2.radius();

    // Concentric: constant distance |r1 − r2|.
    let d12 = o1.distance(&o2);
    if d12 < CONFUSION {
        let u = GpVec2d::new(c1.pos.vxdir.x, c1.pos.vxdir.y);
        let p1 = o1.translated_vec(&u.multiplied_scalar(r1));
        let p2 = o2.translated_vec(&u.multiplied_scalar(r2));
        return vec![pair2d_cc(p1, circle_param2d(c1, &p1), p2, circle_param2d(c2, &p2))];
    }

    let dir = GpVec2d::new(o2.x() - o1.x(), o2.y() - o1.y()).multiplied_scalar(1.0 / d12);
    let mut out = Vec::with_capacity(6);

    // Four collinear critical pairs.
    let p11 = o1.translated_vec(&dir.multiplied_scalar(-r1));
    let p12 = o1.translated_vec(&dir.multiplied_scalar(r1));
    let p21 = o2.translated_vec(&dir.multiplied_scalar(-r2));
    let p22 = o2.translated_vec(&dir.multiplied_scalar(r2));
    for (pa, pb) in [(p11, p21), (p11, p22), (p12, p21), (p12, p22)] {
        out.push(pair2d_cc(pa, circle_param2d(c1, &pa), pb, circle_param2d(c2, &pb)));
    }

    // Intersections (zero distance).
    let b_out = d12 > (r1 + r2 + CONFUSION);
    let b_in = d12 < (r1 - r2).abs() - CONFUSION;
    if !b_out && !b_in {
        let alpha = 0.5 * (r1 * r1 - r2 * r2 + d12 * d12) / d12;
        let val = r1 * r1 - alpha * alpha;
        let beta = val.abs().sqrt();
        let pt = o1.translated_vec(&dir.multiplied_scalar(alpha));
        // Perpendicular direction in the plane.
        let dlt = GpVec2d::new(-dir.y(), dir.x());
        let pl1 = pt.translated_vec(&dlt.multiplied_scalar(beta));
        let pl2 = pt.translated_vec(&dlt.multiplied_scalar(-beta));
        if pl1.square_distance(&pl2) > CONFUSION * CONFUSION {
            out.push(pair2d_cc(pl1, circle_param2d(c1, &pl1), pl1, circle_param2d(c2, &pl1)));
            out.push(pair2d_cc(pl2, circle_param2d(c1, &pl2), pl2, circle_param2d(c2, &pl2)));
        } else {
            out.push(pair2d_cc(pl1, circle_param2d(c1, &pl1), pl1, circle_param2d(c2, &pl1)));
        }
    }
    dedupe_sort_cc2d(out)
}

pub(super) fn line_of_curve2d(c: &dyn Curve2d) -> Option<GpLin2d> {
    let u0 = if c.first_parameter().is_finite() { c.first_parameter() } else { 0.0 };
    let p = c.d0(u0);
    let d = c.d1(u0).1;
    let m = d.magnitude();
    if m < 1e-12 {
        return None;
    }
    GpDir2d::from_vec2d(&d).ok().map(|dir| GpLin2d::from_pnt_dir(p, dir))
}

/// F1 = (C1−C2)·C1′/|C1′| and F2 = (C1−C2)·C2′/|C2′| plus the Jacobian (2D).
pub(super) fn f_and_jac2d(
    c1: &dyn Curve2d,
    c2: &dyn Curve2d,
    u: f64,
    v: f64,
) -> Option<(f64, f64, f64, f64, f64, f64)> {
    let p1 = c1.d0(u);
    let p2 = c2.d0(v);
    let (_, du, duu) = c1.d2(u);
    let (_, dv, dvv) = c2.d2(v);
    if !p1.x().is_finite() || !p2.x().is_finite() {
        return None;
    }
    let ndu = du.magnitude();
    let ndv = dv.magnitude();
    if ndu < 1e-12 || ndv < 1e-12 {
        return None;
    }
    let d = GpVec2d::new(p1.x() - p2.x(), p1.y() - p2.y());
    let f1 = d.dot(&du) / ndu;
    let f2 = d.dot(&dv) / ndv;
    let j11 = ndu + d.dot(&duu) / ndu - f1 * du.dot(&duu) / (ndu * ndu);
    let j12 = -dv.dot(&du) / ndu;
    let j21 = du.dot(&dv) / ndv;
    let j22 = -ndv + d.dot(&dvv) / ndv - f2 * dv.dot(&dvv) / (ndv * ndv);
    Some((f1, f2, j11, j12, j21, j22))
}

/// Newton refine (clamped to bounds) from a seed (2D).
pub(super) fn refine_curve_curve2d_newton(
    c1: &dyn Curve2d,
    c2: &dyn Curve2d,
    mut u: f64,
    mut v: f64,
    lo1: f64,
    hi1: f64,
    lo2: f64,
    hi2: f64,
) -> (f64, f64) {
    for _ in 0..32 {
        match f_and_jac2d(c1, c2, u, v) {
            Some((f1, f2, j11, j12, j21, j22)) => {
                let det = j11 * j22 - j12 * j21;
                if det.abs() < 1e-300 {
                    break;
                }
                let du = (f2 * j12 - f1 * j22) / det;
                let dv = (f1 * j21 - j11 * f2) / det;
                let un = (u + du).clamp(lo1.min(hi1), lo1.max(hi1));
                let vn = (v + dv).clamp(lo2.min(hi2), lo2.max(hi2));
                let converged = (un - u).abs() < 1e-12 * (1.0 + u.abs())
                    && (vn - v).abs() < 1e-12 * (1.0 + v.abs());
                u = un;
                v = vn;
                if converged {
                    break;
                }
            }
            None => break,
        }
    }
    (u, v)
}

/// All local extrema of |C1−C2| in the plane via grid → Newton (2D).
pub(super) fn newton_curve_curve_all2d(c1: &dyn Curve2d, c2: &dyn Curve2d) -> Vec<Extrema2d> {
    let (a1, b1) = (c1.first_parameter(), c1.last_parameter());
    let (a2, b2) = (c2.first_parameter(), c2.last_parameter());
    let us = build_samples2d(c1, a1, b1);
    let vs = build_samples2d(c2, a2, b2);
    let n1 = us.len();
    let n2 = vs.len();

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
            let d = GpVec2d::new(p1.x() - p2.x(), p1.y() - p2.y());
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

    let mut out: Vec<Extrema2d> = Vec::new();
    for (u0, v0) in seeds {
        let (u, v) = refine_curve_curve2d_newton(c1, c2, u0, v0, a1, b1, a2, b2);
        let p1 = c1.d0(u);
        let p2 = c2.d0(v);
        if p1.x().is_finite() && p2.x().is_finite() {
            out.push(pair2d_cc(p1, u, p2, v));
        }
    }
    dedupe_sort_cc2d(out)
}

/// Deduplicate curve–curve extrema by parameter proximity, sorted by distance.
pub(super) fn dedupe_sort_cc2d(v: Vec<Extrema2d>) -> Vec<Extrema2d> {
    let mut out: Vec<Extrema2d> = Vec::new();
    for e in v {
        if let Some(o) = out
            .iter_mut()
            .find(|o| (o.u1 - e.u1).abs() < 1e-5 && (o.u2 - e.u2).abs() < 1e-5)
        {
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

pub(super) fn param_in2d(u: f64, lo: f64, hi: f64) -> bool {
    let lo_ok = !lo.is_finite() || u >= lo - CONFUSION;
    let hi_ok = !hi.is_finite() || u <= hi + CONFUSION;
    lo_ok && hi_ok
}

/// Map `u` into `[lo, lo+period)` for a periodic curve (OCCT `ElCLib::InPeriod`).
pub(super) fn in_period2d(u: f64, lo: f64, period: f64) -> f64 {
    if period <= 0.0 || !lo.is_finite() {
        u
    } else {
        lo + (u - lo).rem_euclid(period)
    }
}

/// All local extrema (minima and maxima) of the distance between two planar
/// curves, deduplicated and sorted by distance. Lines and circles are solved
/// analytically; everything else uses grid + Newton.
pub fn curve_curve_extrema2d_all(c1: &dyn Curve2d, c2: &dyn Curve2d) -> Vec<Extrema2d> {
    let (a1, b1) = (c1.first_parameter(), c1.last_parameter());
    let (a2, b2) = (c2.first_parameter(), c2.last_parameter());

    let l1 = is_line2d(c1).then(|| line_of_curve2d(c1)).flatten();
    let l2 = is_line2d(c2).then(|| line_of_curve2d(c2)).flatten();
    let g1 = classify_circle2d(c1);
    let g2 = classify_circle2d(c2);

    let analytic: Vec<Extrema2d> = match (l1, l2, g1, g2) {
        (Some(l1), Some(l2), _, _) => line_line_extrema2d(&l1, &l2),
        (Some(l), None, _, Some(c)) => line_circle_extrema2d(&l, &c),
        (None, Some(l), Some(c), _) => {
            let mut out = line_circle_extrema2d(&l, &c);
            for e in out.iter_mut() {
                std::mem::swap(&mut e.p1, &mut e.p2);
                std::mem::swap(&mut e.u1, &mut e.u2);
            }
            out
        }
        (None, None, Some(g1), Some(g2)) => circle_circle_extrema2d(&g1, &g2),
        _ => Vec::new(),
    };

    if !analytic.is_empty() {
        // Map periodic parameters into the curve range before the range check.
        let p1 = if c1.is_periodic() { c1.period() } else { 0.0 };
        let p2 = if c2.is_periodic() { c2.period() } else { 0.0 };
        let analytic: Vec<Extrema2d> = analytic
            .into_iter()
            .map(|mut e| {
                e.u1 = in_period2d(e.u1, a1, p1);
                e.u2 = in_period2d(e.u2, a2, p2);
                e
            })
            .collect();
        let filtered: Vec<Extrema2d> = analytic
            .into_iter()
            .filter(|e| param_in2d(e.u1, a1, b1) && param_in2d(e.u2, a2, b2))
            .collect();
        if !filtered.is_empty() {
            return dedupe_sort_cc2d(filtered);
        }
    }
    dedupe_sort_cc2d(newton_curve_curve_all2d(c1, c2))
}

/// Local extremum of the distance between two planar curves from seed
/// parameters (Newton refinement). Port of `Extrema_LocateExtCC2d`.
pub fn locate_extcc2d(c1: &dyn Curve2d, c2: &dyn Curve2d, u0: f64, v0: f64) -> Option<Extrema2d> {
    let (a1, b1) = (c1.first_parameter(), c1.last_parameter());
    let (a2, b2) = (c2.first_parameter(), c2.last_parameter());
    let (u, v) = refine_curve_curve2d_newton(c1, c2, u0, v0, a1, b1, a2, b2);
    let p1 = c1.d0(u);
    let p2 = c2.d0(v);
    if p1.x().is_finite() && p2.x().is_finite() {
        Some(pair2d_cc(p1, u, p2, v))
    } else {
        None
    }
}

/// All curve–curve intersection points (zero-distance extrema).
pub fn curve_curve_intersections2d(c1: &dyn Curve2d, c2: &dyn Curve2d, tol: f64) -> Vec<GpPnt2d> {
    let (a1, b1) = bound(c1);
    let (a2, b2) = bound(c2);
    let n = 32;
    let mut pts: Vec<GpPnt2d> = Vec::new();
    for i in 0..=n {
        let u = a1 + (b1 - a1) * i as f64 / n as f64;
        for j in 0..=n {
            let v = a2 + (b2 - a2) * j as f64 / n as f64;
            let p = c1.d0(u);
            let q = c2.d0(v);
            if p.distance(&q) < tol {
                if pts.iter().all(|x| x.distance(&p) > 1e-6) {
                    pts.push(p);
                }
            }
        }
    }
    pts
}

/// Distance from a point to the closest point on a polyline.
pub fn point_polyline_extrema2d(poly: &[GpPnt2d], p: &GpPnt2d) -> Extrema2d {
    let mut best: Option<Extrema2d> = None;
    for w in poly.windows(2) {
        let (a, b) = (w[0], w[1]);
        let abx = b.x() - a.x();
        let aby = b.y() - a.y();
        let len2 = abx * abx + aby * aby;
        let t = if len2 < 1e-30 {
            0.0
        } else {
            (((p.x() - a.x()) * abx + (p.y() - a.y()) * aby) / len2).clamp(0.0, 1.0)
        };
        let q = GpPnt2d::new(a.x() + t * abx, a.y() + t * aby);
        let d = p.distance(&q);
        if best.as_ref().map_or(true, |e: &Extrema2d| d < e.distance) {
            best = Some(Extrema2d { p1: *p, p2: q, distance: d, u1: t, u2: t });
        }
    }
    best.unwrap()
}

/// Tangent direction of a 2D curve at `u` (normalized).
pub fn tangent2d(c: &dyn Curve2d, u: f64) -> GpVec2d {
    let (_, d) = c.d1(u);
    let m = d.magnitude();
    if m > 1e-12 {
        GpVec2d::new(d.x() / m, d.y() / m)
    } else {
        GpVec2d::new(1.0, 0.0)
    }
}
