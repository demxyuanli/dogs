use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Polynomial root helpers.
//
// `occt-geom` does not depend on `occt-math`, so the closed-form quadratic /
// cubic and the Ferrari quartic are kept local. Ports `math_DirectPolynomial
// Roots` and the trigonometric reduction used by `math_TrigonometricFunction
// Roots`.
// ---------------------------------------------------------------------------

/// Horner evaluation of a·x⁴ + b·x³ + c·x² + d·x + e.

pub(super) fn poly4(a: f64, b: f64, c: f64, d: f64, e: f64, x: f64) -> f64 {
    ((a * x + b) * x + c) * x * x + d * x + e
}

/// Real roots of a·x² + b·x + c = 0 (up to two, sorted ascending).
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

/// Real roots of a·x³ + b·x² + c·x + d = 0 (Cardano / trigonometric form).
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
            roots.push(-shift); // triple root
        } else {
            let u = (-q / 2.0).cbrt();
            roots.push(2.0 * u - shift);
            roots.push(-u - shift);
        }
    } else {
        let r = 2.0 * (-p / 3.0).sqrt();
        let theta = ((3.0 * q / (2.0 * p)) * (-3.0 / p).sqrt()).acos() / 3.0;
        let two_pi = 2.0 * std::f64::consts::PI;
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

/// Newton-polish a root of a·x⁴ + b·x³ + c·x² + d·x + e.
pub(super) fn polish4(a: f64, b: f64, c: f64, d: f64, e: f64, mut x: f64) -> f64 {
    for _ in 0..10 {
        // f′ = 4a·x³ + 3b·x² + 2c·x + d
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
    // Normalize to monic x⁴ + A x³ + B x² + C x + D.
    let a = a3 / a4;
    let b = a2 / a4;
    let c = a1 / a4;
    let d = a0 / a4;
    // Depress: x = y - A/4 → y⁴ + p y² + q y + r = 0.
    let p = b - 3.0 * a * a / 8.0;
    let q = c - a * b / 2.0 + a * a * a / 8.0;
    let r = d - a * c / 4.0 + a * a * b / 16.0 - 3.0 * a * a * a * a / 256.0;
    // Resolvent cubic: 8m³ - 4p m² - 8r m + (4pr - q²) = 0.
    let ms = cubic_roots(8.0, -4.0 * p, -8.0 * r, 4.0 * p * r - q * q);
    // Pick the largest root with 2m - p >= 0 (guaranteed to exist for a
    // quartic with real roots).
    let mut m: Option<f64> = None;
    for mm in ms {
        if 2.0 * mm - p >= 0.0 && m.map_or(true, |best| mm > best) {
            m = Some(mm);
        }
    }
    let m = match m {
        Some(m) => m,
        None => return Vec::new(), // no real roots
    };
    let two_m_p = 2.0 * m - p;
    let mut yroots: Vec<f64> = Vec::new();
    if two_m_p < 1e-12 {
        // s = √(2m-p) ≈ 0 forces q ≈ 0: (y²+m)² = m²-r, so y² = -m ± √(m²-r).
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
        // y² + m = ±(s·y - t) → two quadratics.
        yroots.extend(quadratic_roots(1.0, -s, m + t));
        yroots.extend(quadratic_roots(1.0, s, m - t));
    }
    let mut out: Vec<f64> = Vec::new();
    for y in yroots {
        let x = y - a / 4.0;
        // Filter: only keep candidates that actually satisfy the polynomial.
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
    let period = 2.0 * std::f64::consts::PI;
    let k0 = ((uinf - u) / period).ceil() as i64;
    let k1 = ((usup - u) / period).floor() as i64;
    for k in k0..=k1 {
        roots.push(u + k as f64 * period);
    }
}

/// Newton-polish a root of `d·cos·sin + e·cos + f·sin = 0`.
pub(super) fn polish_trig(d: f64, e: f64, f: f64, mut u: f64) -> f64 {
    for _ in 0..8 {
        let (s, c) = u.sin_cos();
        let g = d * s * c + e * c + f * s;
        let gp = d * (c * c - s * s) - e * s + f * c;
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

/// All roots of `d·cos(u)·sin(u) + e·cos(u) + f·sin(u) = 0` in [uinf, usup].
///
/// Substituting t = tan(u/2) yields the quartic
/// `-e·t⁴ + 2(f-d)·t³ + 2(f+d)·t + e = 0`; u = π (t = ∞) is a root iff e = 0.
pub(super) fn trig_roots_sincos(d: f64, e: f64, f: f64, uinf: f64, usup: f64) -> Vec<f64> {
    let eps = 1e-13;
    let mut roots: Vec<f64> = Vec::new();
    let ts = quartic_roots(-e, 2.0 * (f - d), 0.0, 2.0 * (f + d), e);
    for t in ts {
        push_periodic(&mut roots, 2.0 * t.atan(), uinf, usup);
    }
    if e.abs() < eps {
        push_periodic(&mut roots, std::f64::consts::PI, uinf, usup);
    }
    for r in roots.iter_mut() {
        *r = polish_trig(d, e, f, *r);
    }
    roots.retain(|u| *u >= uinf - 1e-8 && *u <= usup + 1e-8);
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    roots.dedup_by(|x, y| (*x - *y).abs() < 1e-8);
    roots
}

/// Map `u` into [uinf, uinf+period) (OCCT `ElCLib::AdjustPeriodic`); a value
/// landing exactly at the upper bound is wrapped back to the lower bound.
pub(super) fn adjust_periodic(uinf: f64, period: f64, tol: f64, u: f64) -> f64 {
    let mut r = uinf + (u - uinf).rem_euclid(period);
    if (uinf + period - r) < tol {
        r = uinf;
    }
    r
}

// ---------------------------------------------------------------------------
// Analytic solvers (port of `Extrema_ExtPElC`).
// ---------------------------------------------------------------------------

pub(super) fn pair_from(p: &GpPnt, u: f64, q: GpPnt) -> ExtremaPair {
    ExtremaPair { p1: *p, p2: q, distance: p.distance(&q), u1: u, v1: None, u2: u, v2: None }
}

/// Point to line: the unique closest point, `u = (P-O)·dir`.
pub(super) fn line_all(l: &GpLin, p: &GpPnt, uinf: f64, usup: f64) -> Vec<ExtremaPair> {
    let dir = GpVec::from_xyz(l.direction().xyz());
    let loc = l.location();
    let v = GpVec::from_pnts(&loc, p);
    let u = v.dot(&dir);
    if u >= uinf - CONFUSION && u <= usup + CONFUSION {
        vec![pair_from(p, u, loc.translated_vec(&dir.multiplied_scalar(u)))]
    } else {
        Vec::new()
    }
}

/// Point to circle: near point (minimum) at angle of the projected point and
/// the opposite (maximum). Port of `Extrema_ExtPElC::Perform(gp_Circ)`.
pub(super) fn circle_all(c: &GpCirc, p: &GpPnt, uinf: f64, usup: f64) -> Vec<ExtremaPair> {
    let o = c.location();
    let axe = *c.axis().direction();
    let op = GpVec::from_pnts(&o, p);
    let pp = p.translated_vec(&GpVec::from_xyz(axe.xyz()).multiplied_scalar(-op.dot(&GpVec::from_xyz(axe.xyz()))));
    let opp = GpVec::from_pnts(&o, &pp);
    if opp.magnitude() < CONFUSION {
        return Vec::new(); // point on the circle axis: infinite solutions
    }
    let xdir = *c.x_axis().direction();
    let opp_dir = match GpDir::from_vec(&opp) {
        Ok(d) => d,
        Err(_) => return Vec::new(),
    };
    let mut u0 = xdir.angle_with_ref(&opp_dir, &axe); // in (-π, π]
    if u0 + std::f64::consts::PI < ANGULAR {
        u0 = -std::f64::consts::PI;
    } else if u0 - std::f64::consts::PI > -ANGULAR {
        u0 = std::f64::consts::PI;
    }
    let u1 = u0 + std::f64::consts::PI;
    let r = c.radius();
    let tolu = if r > RESOLUTION { CONFUSION / r } else { f64::INFINITY };
    let period = 2.0 * std::f64::consts::PI;
    let mut out = Vec::new();
    for us in [u0, u1] {
        let u = adjust_periodic(uinf, period, tolu, us);
        if u >= uinf - tolu && u <= usup + tolu {
            out.push(pair_from(p, u, clib::circle_value(c, u)));
        }
    }
    out
}

/// Point to ellipse. Port of `Extrema_ExtPElC::Perform(gp_Elips)`: solve
/// (B²-A²)·cos·sin + A·X·sin + B·Y·cos = 0 for u. The sign of the Y·cos term
/// matches this crate's ellipse parameterization C(u) = O + A·cos·X − B·sin·Y.
pub(super) fn ellipse_all(e: &GpElips, p: &GpPnt, uinf: f64, usup: f64) -> Vec<ExtremaPair> {
    let o = e.location();
    let axe = *e.axis().direction();
    let op = GpVec::from_pnts(&o, p);
    let pp = p.translated_vec(&GpVec::from_xyz(axe.xyz()).multiplied_scalar(-op.dot(&GpVec::from_xyz(axe.xyz()))));
    let a = e.major_radius();
    let b = e.minor_radius();
    let opp = GpVec::from_pnts(&o, &pp);
    if opp.magnitude() < CONFUSION && (a - b).abs() < CONFUSION {
        return Vec::new(); // concentric point on a circle: infinite solutions
    }
    let xdir = *e.x_axis().direction();
    let ydir = *e.y_axis().direction();
    let x = opp.dot(&GpVec::from_xyz(xdir.xyz()));
    let y = opp.dot(&GpVec::from_xyz(ydir.xyz()));
    let us = trig_roots_sincos(b * b - a * a, b * y, a * x, uinf, usup);
    us.into_iter().map(|u| pair_from(p, u, clib::ellipse_value(e, u))).collect()
}

/// Point to hyperbola. Port of `Extrema_ExtPElC::Perform(gp_Hypr)`: quartic
/// in v = e^u, keep roots v > 0.
pub(super) fn hyperbola_all(h: &GpHypr, p: &GpPnt, uinf: f64, usup: f64) -> Vec<ExtremaPair> {
    let o = h.location();
    let axe = *h.axis().direction();
    let op = GpVec::from_pnts(&o, p);
    let pp = p.translated_vec(&GpVec::from_xyz(axe.xyz()).multiplied_scalar(-op.dot(&GpVec::from_xyz(axe.xyz()))));
    let r = h.major_radius;
    let s = h.minor_radius;
    let opp = GpVec::from_pnts(&o, &pp);
    let xdir = *h.x_axis().direction();
    let ydir = *h.y_axis().direction();
    let x = opp.dot(&GpVec::from_xyz(xdir.xyz()));
    let y = opp.dot(&GpVec::from_xyz(ydir.xyz()));
    let c1 = (r * r + s * s) / 4.0;
    let vs = quartic_roots(c1, -(x * r + y * s) / 2.0, 0.0, (x * r - y * s) / 2.0, -c1);
    let tol2 = CONFUSION * CONFUSION;
    let mut out = Vec::new();
    let mut seen: Vec<GpPnt> = Vec::new();
    for v in vs {
        if v > 0.0 {
            let u = v.ln();
            if u >= uinf && u <= usup {
                let q = clib::hyperbola_value(h, u);
                if seen.iter().all(|t| t.square_distance(&q) >= tol2) {
                    seen.push(q);
                    out.push(pair_from(p, u, q));
                }
            }
        }
    }
    out
}

/// Point to parabola. Port of `Extrema_ExtPElC::Perform(gp_Parab)`: cubic
/// (1/(4F))·u³ + (2F-X)·u − 2F·Y = 0.
pub(super) fn parabola_all(pa: &GpParab, p: &GpPnt, uinf: f64, usup: f64) -> Vec<ExtremaPair> {
    let o = pa.location();
    let axe = *pa.axis().direction();
    let op = GpVec::from_pnts(&o, p);
    let pp = p.translated_vec(&GpVec::from_xyz(axe.xyz()).multiplied_scalar(-op.dot(&GpVec::from_xyz(axe.xyz()))));
    let f = pa.focal;
    if f.abs() < RESOLUTION {
        return Vec::new();
    }
    let opp = GpVec::from_pnts(&o, &pp);
    let xdir = *pa.x_axis().direction();
    let ydir = *pa.y_axis().direction();
    let x = opp.dot(&GpVec::from_xyz(xdir.xyz()));
    let y = opp.dot(&GpVec::from_xyz(ydir.xyz()));
    let us = cubic_roots(1.0 / (4.0 * f), 0.0, 2.0 * f - x, -2.0 * f * y);
    let tol2 = CONFUSION * CONFUSION;
    let mut out = Vec::new();
    let mut seen: Vec<GpPnt> = Vec::new();
    for u in us {
        if u >= uinf && u <= usup {
            let q = clib::parabola_value(pa, u);
            if seen.iter().all(|t| t.square_distance(&q) >= tol2) {
                seen.push(q);
                out.push(pair_from(p, u, q));
            }
        }
    }
    out
}

/// Public: closest point on the line.
pub fn point_line_extrema(l: &GpLin, p: &GpPnt) -> ExtremaPair {
    let inf = f64::INFINITY;
    line_all(l, p, -inf, inf).first().copied().unwrap_or_else(|| pair_from(p, 0.0, clib::line_value(l, 0.0)))
}

/// Public: closest point on the circle (full range [0, 2π]).
pub fn point_circle_extrema(c: &GpCirc, p: &GpPnt) -> ExtremaPair {
    let period = 2.0 * std::f64::consts::PI;
    circle_all(c, p, 0.0, period)
        .into_iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal))
        .unwrap_or_else(|| pair_from(p, 0.0, clib::circle_value(c, 0.0)))
}

/// Public: closest point on the ellipse (full range [0, 2π]).
pub fn point_ellipse_extrema(e: &GpElips, p: &GpPnt) -> ExtremaPair {
    let period = 2.0 * std::f64::consts::PI;
    ellipse_all(e, p, 0.0, period)
        .into_iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal))
        .unwrap_or_else(|| pair_from(p, 0.0, clib::ellipse_value(e, 0.0)))
}

/// Public: closest point on the hyperbola.
pub fn point_hyperbola_extrema(h: &GpHypr, p: &GpPnt) -> ExtremaPair {
    let inf = f64::INFINITY;
    hyperbola_all(h, p, -inf, inf)
        .into_iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal))
        .unwrap_or_else(|| pair_from(p, 0.0, clib::hyperbola_value(h, 0.0)))
}

/// Public: closest point on the parabola.
pub fn point_parabola_extrema(pa: &GpParab, p: &GpPnt) -> ExtremaPair {
    let inf = f64::INFINITY;
    parabola_all(pa, p, -inf, inf)
        .into_iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal))
        .unwrap_or_else(|| pair_from(p, 0.0, clib::parabola_value(pa, 0.0)))
}
