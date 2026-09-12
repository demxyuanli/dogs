use super::prelude::*;
use super::*;

pub(super) const PI: f64 = std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Polynomial root helpers.
//
// Local copies of the kernels in `extrema_pc` (those are private). Ports of
// `math_DirectPolynomialRoots` / `math_TrigonometricFunctionRoots`.
// ---------------------------------------------------------------------------

pub(super) fn poly4(a: f64, b: f64, c: f64, d: f64, e: f64, x: f64) -> f64 {
    ((a * x + b) * x + c) * x * x + d * x + e
}

pub(super) fn quadratic_roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() < 1e-300 {
        if b.abs() < 1e-300 {
            return Vec::new();
        }
        return vec![-c / b];
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return Vec::new();
    }
    let sq = disc.sqrt();
    let q = -0.5 * (b + b.signum() * sq);
    let r1 = q / a;
    let r2 = if q.abs() > 1e-300 { c / q } else { (-b - sq) / (2.0 * a) };
    let mut v = vec![r1, r2];
    v.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    v
}

pub(super) fn cubic_roots(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a.abs() < 1e-300 {
        return quadratic_roots(b, c, d);
    }
    let b = b / a;
    let c = c / a;
    let d = d / a;
    let p = c - b * b / 3.0;
    let q = 2.0 * b * b * b / 27.0 - b * c / 3.0 + d;
    let disc = q * q / 4.0 + p * p * p / 27.0;
    let shift = b / 3.0;
    let mut roots = Vec::new();
    if disc > 1e-14 {
        let sq = disc.sqrt();
        let u = (-q / 2.0 + sq).cbrt();
        let v = (-q / 2.0 - sq).cbrt();
        roots.push(u + v - shift);
    } else if disc >= -1e-14 {
        if q.abs() < 1e-14 {
            roots.push(-shift);
        } else {
            let u = (-q / 2.0).cbrt();
            roots.push(2.0 * u - shift);
            roots.push(-u - shift);
        }
    } else {
        let r = 2.0 * (-p / 3.0).sqrt();
        let theta = ((3.0 * q / (2.0 * p)) * (-3.0 / p).sqrt()).acos() / 3.0;
        let two_pi = 2.0 * PI;
        roots.push(r * theta.cos() - shift);
        roots.push(r * (theta - two_pi / 3.0).cos() - shift);
        roots.push(r * (theta + two_pi / 3.0).cos() - shift);
    }
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    let mut out: Vec<f64> = Vec::new();
    for x in roots {
        if out.last().map_or(true, |&l| (x - l).abs() > 1e-9 * (1.0 + x.abs())) {
            out.push(x);
        }
    }
    out
}

pub(super) fn polish4(a: f64, b: f64, c: f64, d: f64, e: f64, mut x: f64) -> f64 {
    for _ in 0..10 {
        let fp = (4.0 * a * x + 3.0 * b) * x * x + 2.0 * c * x + d;
        if fp.abs() < 1e-300 {
            break;
        }
        let xn = x - poly4(a, b, c, d, e, x) / fp;
        if !xn.is_finite() || (xn - x).abs() > 1.0 + x.abs() {
            break;
        }
        if (xn - x).abs() < 1e-13 * (1.0 + x.abs()) {
            return xn;
        }
        x = xn;
    }
    x
}

/// Real roots of a·x⁴ + b·x³ + c·x² + d·x + e = 0 (Ferrari method).
pub(super) fn quartic_roots(a4: f64, a3: f64, a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    if a4.abs() < 1e-300 {
        return cubic_roots(a3, a2, a1, a0);
    }
    let a = a3 / a4;
    let b = a2 / a4;
    let c = a1 / a4;
    let d = a0 / a4;
    let p = b - 3.0 * a * a / 8.0;
    let q = c - a * b / 2.0 + a * a * a / 8.0;
    let r = d - a * c / 4.0 + a * a * b / 16.0 - 3.0 * a * a * a * a / 256.0;
    let ms = cubic_roots(8.0, -4.0 * p, -8.0 * r, 4.0 * p * r - q * q);
    let mut m: Option<f64> = None;
    for mm in ms {
        if 2.0 * mm - p >= 0.0 && m.map_or(true, |best| mm > best) {
            m = Some(mm);
        }
    }
    let m = match m {
        Some(m) => m,
        None => return Vec::new(),
    };
    let two_m_p = 2.0 * m - p;
    let mut yroots: Vec<f64> = Vec::new();
    if two_m_p < 1e-12 {
        let rad = m * m - r;
        if rad < 0.0 {
            return Vec::new();
        }
        let sr = rad.sqrt();
        for ym in [sr, -sr] {
            let val = -m + ym;
            if val >= 0.0 {
                let y = val.sqrt();
                yroots.push(y);
                if y > 1e-12 {
                    yroots.push(-y);
                }
            }
        }
    } else {
        let s = two_m_p.sqrt();
        let t = q / (2.0 * s);
        yroots.extend(quadratic_roots(1.0, -s, m + t));
        yroots.extend(quadratic_roots(1.0, s, m - t));
    }
    let mut out: Vec<f64> = Vec::new();
    for y in yroots {
        let x = y - a / 4.0;
        if poly4(1.0, a, b, c, d, x).abs() < 1e-6 {
            out.push(polish4(1.0, a, b, c, d, x));
        }
    }
    out.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    let mut dedup: Vec<f64> = Vec::new();
    for x in out {
        if dedup.last().map_or(true, |&l| (x - l).abs() > 1e-8 * (1.0 + x.abs())) {
            dedup.push(x);
        }
    }
    dedup
}

/// Push `u + 2πk` for every integer k that lands in [uinf, usup].
pub(super) fn push_periodic(roots: &mut Vec<f64>, u: f64, uinf: f64, usup: f64) {
    if !(uinf.is_finite() && usup.is_finite()) {
        return;
    }
    let period = 2.0 * PI;
    let k0 = ((uinf - u) / period).ceil() as i64;
    let k1 = ((usup - u) / period).floor() as i64;
    for k in k0..=k1 {
        roots.push(u + k as f64 * period);
    }
}

/// Evaluate `cc·cos² + 2·sc·cos·sin + c·cos + s·sin + cte` and its derivative.
pub(super) fn trig_poly(cc: f64, sc: f64, c: f64, s: f64, cte: f64, u: f64) -> (f64, f64) {
    let (su, cu) = u.sin_cos();
    let g = cc * cu * cu + 2.0 * sc * cu * su + c * cu + s * su + cte;
    let gp = -2.0 * cc * cu * su + 2.0 * sc * (cu * cu - su * su) - c * su + s * cu;
    (g, gp)
}

/// Newton-polish a root of the full trigonometric polynomial.
pub(super) fn polish_trig_full(cc: f64, sc: f64, c: f64, s: f64, cte: f64, mut u: f64) -> f64 {
    for _ in 0..8 {
        let (g, gp) = trig_poly(cc, sc, c, s, cte, u);
        if gp.abs() < 1e-300 {
            break;
        }
        let un = u - g / gp;
        if !un.is_finite() || (un - u).abs() < 1e-13 * (1.0 + u.abs()) {
            break;
        }
        u = un;
    }
    u
}

/// All roots of `cc·cos² + 2·sc·cos·sin + c·cos + s·sin + cte = 0` in
/// [uinf, usup].  Substituting t = tan(u/2) yields the quartic
/// (cc−c+cte)·t⁴ + (−4sc+2s)·t³ + (−2cc+2cte)·t² + (4sc+2s)·t + (cc+c+cte) = 0;
/// u = π (t = ∞) is a root iff cc−c+cte = 0.
pub(super) fn trig_roots_full(cc: f64, sc: f64, c: f64, s: f64, cte: f64, uinf: f64, usup: f64) -> Vec<f64> {
    let eps = 1e-13;
    let a = cc - c + cte;
    let b = -4.0 * sc + 2.0 * s;
    let c2 = -2.0 * cc + 2.0 * cte;
    let d = 4.0 * sc + 2.0 * s;
    let e = cc + c + cte;
    let mut roots: Vec<f64> = Vec::new();
    for t in quartic_roots(a, b, c2, d, e) {
        push_periodic(&mut roots, 2.0 * t.atan(), uinf, usup);
    }
    if a.abs() < eps {
        push_periodic(&mut roots, PI, uinf, usup);
    }
    for r in roots.iter_mut() {
        *r = polish_trig_full(cc, sc, c, s, cte, *r);
    }
    let amax = cc.abs().max(sc.abs()).max(c.abs()).max(s.abs()).max(cte.abs());
    let tol = 1e-8 * (1.0 + amax);
    roots.retain(|u| *u >= uinf - 1e-8 && *u <= usup + 1e-8 && trig_poly(cc, sc, c, s, cte, *u).0.abs() < tol);
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    roots.dedup_by(|x, y| (*x - *y).abs() < 1e-8);
    roots
}

// ---------------------------------------------------------------------------
// Analytic solvers (port of `Extrema_ExtElC`).
// ---------------------------------------------------------------------------

pub(super) fn pair_cc(p1: GpPnt, u1: f64, p2: GpPnt, u2: f64) -> ExtremaPair {
    ExtremaPair { p1, p2, distance: p1.distance(&p2), u1, v1: None, u2, v2: None }
}

pub(super) fn circle_param(c: &GpCirc, p: &GpPnt) -> f64 {
    let o = c.location();
    let v = GpVec::from_pnts(&o, p);
    let xd = *c.x_axis().direction();
    let yd = *c.y_axis().direction();
    // C(u) = O + r·cos·xd + r·sin·yd → u = atan2(v·yd, v·xd).
    v.dot(&GpVec::from_xyz(yd.xyz())).atan2(v.dot(&GpVec::from_xyz(xd.xyz())))
}

pub(super) fn line_value(l: &GpLin, u: f64) -> GpPnt {
    clib::line_value(l, u)
}

pub(super) fn circle_value(c: &GpCirc, u: f64) -> GpPnt {
    clib::circle_value(c, u)
}

/// Exact extrema between two (unbounded) straight lines.
///
/// Port of `Extrema_ExtElC(const gp_Lin&, const gp_Lin&, AngTol)`. Skew lines
/// have a single closest pair; parallel lines have constant distance.
pub fn line_line_extrema(l1: &GpLin, l2: &GpLin) -> Vec<ExtremaPair> {
    let d1 = l1.direction();
    let d2 = l2.direction();
    let cosa = d1.dot(&d2);
    let sqsina = 1.0 - cosa * cosa;
    if sqsina < RESOLUTION || d1.is_parallel(&d2) {
        // Parallel: distance is the distance from any point of C1 to C2.
        let p1 = l1.location();
        let u2 = GpVec::from_pnts(&l2.location(), &p1).dot(&GpVec::from_xyz(d2.xyz()));
        let p2 = line_value(l2, u2);
        return vec![pair_cc(p1, 0.0, p2, u2)];
    }
    let l1l2 = GpVec::from_pnts(&l1.location(), &l2.location());
    let d1l = GpVec::from_xyz(d1.xyz()).dot(&l1l2);
    let d2l = GpVec::from_xyz(d2.xyz()).dot(&l1l2);
    let u1 = (d1l - cosa * d2l) / sqsina;
    let u2 = (cosa * d1l - d2l) / sqsina;
    let p1 = line_value(l1, u1);
    let p2 = line_value(l2, u2);
    vec![pair_cc(p1, u1, p2, u2)]
}

/// Extrema between a straight line and a circle.
///
/// Port of `Extrema_ExtElC(const gp_Lin&, const gp_Circ&, Tol)`. Lines parallel
/// to the circle plane are solved in the plane (2D line–circle plus
/// intersections); otherwise the full 3D trigonometric condition
/// `A1·cos² + 2·A2·cos·sin + A3·cos + A4·sin + A5 = 0` is solved.
pub fn line_circle_extrema(l: &GpLin, c: &GpCirc) -> Vec<ExtremaPair> {
    let dir_c = *c.axis().direction();
    let dir_l = l.direction();
    let o1 = l.location();
    let o2 = c.location();
    let r = c.radius();

    // Planar (line parallel to the circle plane): solve in the circle frame.
    if dir_c.dot(&dir_l).abs() <= ANGULAR {
        return planar_line_circle(l, c);
    }

    // Full 3D case. Work in the circle's reference frame.
    let x2 = *c.x_axis().direction();
    let y2 = *c.y_axis().direction();
    let dx = dir_l.dot(&x2);
    let dy = dir_l.dot(&y2);
    let _dz = dir_l.dot(&dir_c);

    // O2O1 in the circle frame.
    let o2o1 = GpVec::from_pnts(&o2, &o1);
    let vx = o2o1.dot(&GpVec::from_xyz(x2.xyz()));
    let vy = o2o1.dot(&GpVec::from_xyz(y2.xyz()));

    // V = (D·O2O1)·D − O2O1 in the circle frame.
    let do2o1 = dx * vx + dy * vy;
    let vx = do2o1 * dx - vx;
    let vy = do2o1 * dy - vy;

    // Coefficients of the trigonometric equation (divided by R for A3/A4).
    let a5 = r * dx * dy;
    let a1 = -2.0 * a5;
    let a2 = 0.5 * r * (dx * dx - dy * dy);
    let a3 = vy;
    let a4 = -vx;

    // Infinite solutions: line is the circle axis (constant distance R).
    if (a1.abs() + a2.abs() + a3.abs() + a4.abs()) < 1e-10 && a5.abs() < 1e-10 {
        let p2 = circle_value(c, 0.0);
        let u1 = GpVec::from_pnts(&o1, &p2).dot(&GpVec::from_xyz(dir_l.xyz()));
        let p1 = line_value(l, u1);
        return vec![pair_cc(p1, u1, p2, 0.0)];
    }

    let us = trig_roots_full(a1, a2, a3, a4, a5, 0.0, 2.0 * PI);
    let mut out = Vec::with_capacity(us.len());
    for u2 in us {
        let p2 = circle_value(c, u2);
        let u1 = GpVec::from_pnts(&o1, &p2).dot(&GpVec::from_xyz(dir_l.xyz()));
        let p1 = line_value(l, u1);
        out.push(pair_cc(p1, u1, p2, u2));
    }
    dedupe_sort(out)
}

/// Planar line–circle (port of `Extrema_ExtElC::PlanarLineCircleExtrema`).
///
/// The line is projected onto the circle plane; the 2D tangent extrema (radius
/// perpendicular to the projected direction) and the 2D intersections are both
/// mapped back to the actual 3D line/circle parameters (the "intersections" of
/// the projection are candidate closest pairs in 3D, not literal crossings).
pub(super) fn planar_line_circle(l: &GpLin, c: &GpCirc) -> Vec<ExtremaPair> {
    let o1 = l.location();
    let o2 = c.location();
    let r = c.radius();
    let x2 = *c.x_axis().direction();
    let y2 = *c.y_axis().direction();

    // Project the line onto the circle plane.
    let dl = GpVec::from_xyz(l.direction().xyz());
    let p2dx = GpVec::from_xyz(x2.xyz());
    let p2dy = GpVec::from_xyz(y2.xyz());
    let ldx = dl.dot(&p2dx);
    let ldy = dl.dot(&p2dy);
    let o2o1 = GpVec::from_pnts(&o2, &o1);
    let loc2d_x = o2o1.dot(&p2dx);
    let loc2d_y = o2o1.dot(&p2dy);

    let mut out = Vec::new();

    // 2D tangent extrema: radius perpendicular to the projected line direction.
    let mut tetas: Vec<f64> = Vec::new();
    if ldy.abs() <= 1e-15 {
        tetas.push(PI / 2.0);
    } else {
        tetas.push((-ldx / ldy).atan());
    }
    tetas.push(tetas[0] + PI);
    if tetas[0] < 0.0 {
        tetas[0] += 2.0 * PI;
    }
    for teta in tetas {
        let p2 = circle_value(c, teta);
        let u1 = GpVec::from_pnts(&o1, &p2).dot(&dl);
        let p1 = line_value(l, u1);
        out.push(pair_cc(p1, u1, p2, teta));
    }

    // 2D projected intersections: |(loc2d_x, loc2d_y) + t·(ldx, ldy)|² = r².
    let a = ldx * ldx + ldy * ldy;
    let b = 2.0 * (loc2d_x * ldx + loc2d_y * ldy);
    let c0 = loc2d_x * loc2d_x + loc2d_y * loc2d_y - r * r;
    let disc = b * b - 4.0 * a * c0;
    if disc >= 0.0 && a > 1e-300 {
        let sq = disc.sqrt();
        for t in [(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
            let p1 = line_value(l, t);
            let x2d = loc2d_x + t * ldx;
            let y2d = loc2d_y + t * ldy;
            let u2 = y2d.atan2(x2d);
            let p2 = circle_value(c, u2);
            out.push(pair_cc(p1, t, p2, u2));
        }
    }
    dedupe_sort(out)
}

/// Extrema between two circles.
///
/// Port of `Extrema_ExtElC(const gp_Circ&, const gp_Circ&)`. Only same-plane
/// circles are handled exactly (coaxial constant-distance, the four collinear
/// critical pairs, and intersections). Skew circles return an empty vector and
/// must fall back to the general Newton path.
pub fn circle_circle_extrema(c1: &GpCirc, c2: &GpCirc) -> Vec<ExtremaPair> {
    let oc1 = c1.location();
    let oc2 = c2.location();
    let dc1 = *c1.axis().direction();
    let dc2 = *c2.axis().direction();
    let tol_d2 = CONFUSION * CONFUSION;

    // Same-plane check.
    let plane_dist2 = {
        let v = GpVec::from_pnts(&oc1, &oc2);
        let n = GpVec::from_xyz(dc1.xyz());
        let d = v.dot(&n);
        d * d
    };
    if !(dc1.is_parallel(&dc2) && plane_dist2 < tol_d2) {
        return Vec::new(); // skew circles
    }

    let r1 = c1.radius();
    let r2 = c2.radius();

    // Coaxial (concentric): constant distance |r1 − r2|.
    let d12 = oc1.distance(&oc2);
    if d12 < CONFUSION {
        let u = GpVec::from_xyz(c1.x_axis().direction().xyz());
        let p1 = oc1.translated_vec(&u.multiplied_scalar(r1));
        let p2 = oc2.translated_vec(&u.multiplied_scalar(r2));
        return vec![pair_cc(p1, circle_param(c1, &p1), p2, circle_param(c2, &p2))];
    }

    let dir12 = GpVec::from_pnts(&oc1, &oc2).multiplied_scalar(1.0 / d12);
    let mut out = Vec::with_capacity(6);

    // Four collinear critical pairs.
    let p11 = oc1.translated_vec(&dir12.multiplied_scalar(-r1));
    let p12 = oc1.translated_vec(&dir12.multiplied_scalar(r1));
    let p21 = oc2.translated_vec(&dir12.multiplied_scalar(-r2));
    let p22 = oc2.translated_vec(&dir12.multiplied_scalar(r2));
    for (pa, pb) in [(p11, p21), (p11, p22), (p12, p21), (p12, p22)] {
        out.push(pair_cc(pa, circle_param(c1, &pa), pb, circle_param(c2, &pb)));
    }

    // Intersections (zero distance) when the circles cross or touch.
    let b_out = d12 > (r1 + r2 + CONFUSION);
    let b_in = d12 < (r1 - r2).abs() - CONFUSION;
    if !b_out && !b_in {
        let alpha = 0.5 * (r1 * r1 - r2 * r2 + d12 * d12) / d12;
        let val = r1 * r1 - alpha * alpha;
        let beta = val.abs().sqrt();
        let pt = oc1.translated_vec(&dir12.multiplied_scalar(alpha));
        let dlt = GpVec::from_xyz(dc1.xyz()).crossed(&dir12);
        let pl1 = pt.translated_vec(&dlt.multiplied_scalar(beta));
        let pl2 = pt.translated_vec(&dlt.multiplied_scalar(-beta));
        if pl1.square_distance(&pl2) > tol_d2 {
            out.push(pair_cc(pl1, circle_param(c1, &pl1), pl1, circle_param(c2, &pl1)));
            out.push(pair_cc(pl2, circle_param(c1, &pl2), pl2, circle_param(c2, &pl2)));
        } else {
            out.push(pair_cc(pl1, circle_param(c1, &pl1), pl1, circle_param(c2, &pl1)));
        }
    }
    dedupe_sort(out)
}

// ---------------------------------------------------------------------------
// Classification of a `dyn Curve` into an analytic type.
// ---------------------------------------------------------------------------

/// Constant-tangent-direction detection ⇒ line (same as `extrema_pc`).
pub(super) fn is_line(c: &dyn Curve) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let us: Vec<f64> = if a.is_finite() && b.is_finite() && b > a {
        (0..=5).map(|i| a + (b - a) * i as f64 / 5.0).collect()
    } else {
        vec![-2.0, -1.0, 0.0, 1.0, 2.0]
    };
    let t0 = c.d1(us[0]).1;
    let m0 = t0.magnitude();
    if m0 < 1e-9 {
        return false;
    }
    for u in &us[1..] {
        let t = c.d1(*u).1;
        let m = t.magnitude();
        if m < 1e-9 {
            return false;
        }
        if t0.cross_magnitude(&t) > 1e-7 * m0 * m {
            return false;
        }
    }
    true
}

pub(super) fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let ab = GpVec::from_pnts(a, b);
    let ac = GpVec::from_pnts(a, c);
    let n = ab.crossed(&ac);
    let n2 = n.square_magnitude();
    if n2 < 1e-20 {
        return None;
    }
    let ab2 = ab.square_magnitude();
    let ac2 = ac.square_magnitude();
    let n_ac = n.crossed(&ac);
    let ab_n = ab.crossed(&n);
    let coef = 0.5 / n2;
    Some(GpPnt::new(
        a.x() + coef * (ab2 * n_ac.x() + ac2 * ab_n.x()),
        a.y() + coef * (ab2 * n_ac.y() + ac2 * ab_n.y()),
        a.z() + coef * (ab2 * n_ac.z() + ac2 * ab_n.z()),
    ))
}

/// Reconstruct a circle from a bounded curve whose samples are coplanar and
/// equidistant from a single center (same as `extrema_pc`).
pub(super) fn classify_circle(c: &dyn Curve) -> Option<GpCirc> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite() && b > a) {
        return None;
    }
    let n = 8;
    let pts: Vec<GpPnt> = (0..n).map(|i| c.d0(a + (b - a) * i as f64 / (n - 1) as f64)).collect();
    let e0 = GpVec::from_pnts(&pts[0], &pts[1]);
    let e1 = GpVec::from_pnts(&pts[0], &pts[2]);
    let nrm = e0.crossed(&e1);
    let nmag = nrm.magnitude();
    if nmag < 1e-9 {
        return None;
    }
    for i in 1..n - 2 {
        let ei = GpVec::from_pnts(&pts[i], &pts[i + 1]);
        let ei2 = GpVec::from_pnts(&pts[i], &pts[i + 2]);
        let ni = ei.crossed(&ei2);
        if ni.magnitude() < 1e-9 {
            continue;
        }
        if nrm.cross_magnitude(&ni) > 1e-4 * nmag * ni.magnitude() {
            return None;
        }
    }
    let o = circumcenter(&pts[0], &pts[1], &pts[2])?;
    let r = pts[0].distance(&o);
    if r < 1e-12 {
        return None;
    }
    for p in &pts {
        if (p.distance(&o) - r).abs() > 1e-4 * r {
            return None;
        }
    }
    let normal = GpDir::from_vec(&nrm).ok()?;
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&o, &pts[0])).ok()?;
    let ax = GpAx2::new(o, normal, xdir).ok()?;
    Some(GpCirc { pos: ax, radius: r })
}

pub(super) fn line_of_curve(c: &dyn Curve) -> Option<GpLin> {
    let u0 = if c.first_parameter().is_finite() { c.first_parameter() } else { 0.0 };
    let p = c.d0(u0);
    let d = c.d1(u0).1;
    let m = d.magnitude();
    if m < 1e-12 {
        return None;
    }
    GpDir::from_vec(&d).ok().map(|dir| GpLin::from_pnt_dir(p, dir))
}

// ---------------------------------------------------------------------------
// General Newton path (port of `Extrema_GGenExtCC` + `Extrema_GFuncExtCC`).
// ---------------------------------------------------------------------------

/// F1 = (C1−C2)·C1′/|C1′| and F2 = (C1−C2)·C2′/|C2′| plus the Jacobian.
pub(super) fn f_and_jac(
    c1: &dyn Curve,
    c2: &dyn Curve,
    u: f64,
    v: f64,
) -> Option<(f64, f64, f64, f64, f64, f64)> {
    // Note: `d2` on `GeomLine` returns the point at parameter 0 (clib quirk),
    // so the point is always fetched via `d0`.
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
    let d = GpVec::from_pnts(&p2, &p1); // C1 − C2
    let f1 = d.dot(&du) / ndu;
    let f2 = d.dot(&dv) / ndv;
    let j11 = ndu + d.dot(&duu) / ndu - f1 * du.dot(&duu) / (ndu * ndu);
    let j12 = -dv.dot(&du) / ndu;
    let j21 = du.dot(&dv) / ndv;
    let j22 = -ndv + d.dot(&dvv) / ndv - f2 * dv.dot(&dvv) / (ndv * ndv);
    Some((f1, f2, j11, j12, j21, j22))
}

/// Newton refine (clamped to bounds) from a seed; ports `Extrema_LocECC`.
pub(super) fn refine_curve_curve(
    c1: &dyn Curve,
    c2: &dyn Curve,
    mut u: f64,
    mut v: f64,
    lo1: f64,
    hi1: f64,
    lo2: f64,
    hi2: f64,
) -> (f64, f64) {
    for _ in 0..32 {
        match f_and_jac(c1, c2, u, v) {
            Some((f1, f2, j11, j12, j21, j22)) => {
                let det = j11 * j22 - j12 * j21;
                if det.abs() < 1e-300 {
                    break;
                }
                let du = (f2 * j12 - f1 * j22) / det;
                let dv = (f1 * j21 - j11 * f2) / det;
                let un = (u + du).clamp(lo1.min(hi1), lo1.max(hi1));
                let vn = (v + dv).clamp(lo2.min(hi2), lo2.max(hi2));
                let converged =
                    (un - u).abs() < 1e-12 * (1.0 + u.abs()) && (vn - v).abs() < 1e-12 * (1.0 + v.abs());
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

/// Build a sorted parameter sample grid: uniformly over a bounded range, or
/// expanding windows around 0 for unbounded curves.
pub(super) fn build_samples(_c: &dyn Curve, a: f64, b: f64) -> Vec<f64> {
    let mut s: Vec<f64> = Vec::new();
    if a.is_finite() && b.is_finite() && b > a {
        let span = b - a;
        let n = ((span / 0.1).ceil() as usize).clamp(24, 256);
        for i in 0..=n {
            s.push(a + span * i as f64 / n as f64);
        }
    } else {
        let n = 24;
        for k in 0..9 {
            let w = 2.0f64.powi(k as i32);
            for i in 0..=n {
                s.push(-w + 2.0 * w * i as f64 / n as f64);
            }
        }
    }
    s.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    s.dedup_by(|x, y| (*x - *y).abs() < 1e-12);
    s
}
