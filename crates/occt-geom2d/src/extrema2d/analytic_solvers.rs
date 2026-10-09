use super::prelude::*;

pub(super) const PI: f64 = std::f64::consts::PI;

/// A solved extremum between two 2D objects.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Extrema2d {
    pub p1: GpPnt2d,
    pub p2: GpPnt2d,
    pub distance: f64,
    pub u1: f64,
    pub u2: f64,
}

pub(super) fn bound(c: &dyn Curve2d) -> (f64, f64) {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-1.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Polynomial root helpers (same kernels as `occt_geom::extrema_pc`; kept local
// because occt-geom2d does not depend on occt-geom). Ports of
// `math_DirectPolynomialRoots` / `math_TrigonometricFunctionRoots`.
// ---------------------------------------------------------------------------

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

pub(super) fn poly4(a: f64, b: f64, c: f64, d: f64, e: f64, x: f64) -> f64 {
    ((a * x + b) * x + c) * x * x + d * x + e
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

pub(super) fn adjust_periodic(uinf: f64, period: f64, tol: f64, u: f64) -> f64 {
    let mut r = uinf + (u - uinf).rem_euclid(period);
    if (uinf + period - r) < tol {
        r = uinf;
    }
    r
}

// ---------------------------------------------------------------------------
// 2D analytic solvers (port of `Extrema_ExtPElC2d`).
// ---------------------------------------------------------------------------

/// Point on a 2D circle: O + r·(cos·xd + sin·yd).
pub(super) fn circle2d_val(c: &GpCirc2d, u: f64) -> GpPnt2d {
    let r = c.radius;
    let o = c.location();
    let xd = c.pos.vxdir;
    let yd = c.pos.vydir;
    GpPnt2d::new(o.x() + r * u.cos() * xd.x + r * u.sin() * yd.x, o.y() + r * u.cos() * xd.y + r * u.sin() * yd.y)
}

/// Point on a 2D ellipse: `O + a·cos·xd + b·sin·yd`
/// (`ElCLib::EllipseValue`, `ElCLib.cxx:543-555`, and
/// `Geom2d_Ellipse::EvalD0`, `Geom2d_Ellipse.cxx:256-259`).
pub(super) fn ellipse2d_val(e: &GpElips2d, u: f64) -> GpPnt2d {
    let a = e.major_radius;
    let b = e.minor_radius;
    let o = e.pos.point;
    let xd = e.pos.vxdir;
    let yd = e.pos.vydir;
    GpPnt2d::new(o.x() + a * u.cos() * xd.x + b * u.sin() * yd.x, o.y() + a * u.cos() * xd.y + b * u.sin() * yd.y)
}

/// Point on a 2D hyperbola: O + a·cosh·xd + b·sinh·yd.
pub(super) fn hyperbola2d_val(h: &GpHypr2d, u: f64) -> GpPnt2d {
    let a = h.major_radius;
    let b = h.minor_radius;
    let o = h.pos.point;
    let xd = h.pos.vxdir;
    let yd = h.pos.vydir;
    GpPnt2d::new(o.x() + a * u.cosh() * xd.x + b * u.sinh() * yd.x, o.y() + a * u.cosh() * xd.y + b * u.sinh() * yd.y)
}

/// Point on a 2D parabola: O + (u²/4f)·xd + u·yd.
pub(super) fn parabola2d_val(pa: &GpParab2d, u: f64) -> GpPnt2d {
    let f = pa.focal;
    let o = pa.pos.point;
    let xd = pa.pos.vxdir;
    let yd = pa.pos.vydir;
    let x = u * u / (4.0 * f);
    GpPnt2d::new(o.x() + x * xd.x + u * yd.x, o.y() + x * xd.y + u * yd.y)
}

pub(super) fn pair2d(p: &GpPnt2d, u: f64, q: GpPnt2d) -> Extrema2d {
    Extrema2d { p1: *p, p2: q, distance: p.distance(&q), u1: u, u2: u }
}

pub(super) fn line_all2d(l: &GpLin2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
    let loc = l.pos.loc;
    let dirx = l.pos.vdir.x;
    let diry = l.pos.vdir.y;
    let u = (p.x() - loc.x()) * dirx + (p.y() - loc.y()) * diry;
    if u >= uinf - CONFUSION && u <= usup + CONFUSION {
        vec![pair2d(p, u, GpPnt2d::new(loc.x() + u * dirx, loc.y() + u * diry))]
    } else {
        Vec::new()
    }
}

pub(super) fn circle_all2d(c: &GpCirc2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
    let o = c.location();
    let r = c.radius();
    let vx = o.x() - p.x();
    let vy = o.y() - p.y();
    let m = (vx * vx + vy * vy).sqrt();
    if m < CONFUSION {
        return Vec::new(); // P == center: infinite solutions
    }
    let ux = vx / m;
    let uy = vy / m;
    let xd = c.pos.vxdir;
    let yd = c.pos.vydir;
    let (p1x, p1y) = (o.x() + r * ux, o.y() + r * uy);
    let (p2x, p2y) = (o.x() - r * ux, o.y() - r * uy);
    // Angle of (P1-O) in the circle frame.
    let (ddx, ddy) = (p1x - o.x(), p1y - o.y());
    let u1 = (ddx * yd.x + ddy * yd.y).atan2(ddx * xd.x + ddy * xd.y);
    let u2 = u1 + std::f64::consts::PI;
    let period = 2.0 * std::f64::consts::PI;
    let tolu = if r > RESOLUTION { CONFUSION / r } else { f64::INFINITY };
    let mut out = Vec::new();
    for (us, qx, qy) in [(u1, p1x, p1y), (u2, p2x, p2y)] {
        let u = adjust_periodic(uinf, period, tolu, us);
        if u >= uinf - tolu && u <= usup + tolu {
            out.push(pair2d(p, u, GpPnt2d::new(qx, qy)));
        }
    }
    out
}

pub(super) fn ellipse_all2d(e: &GpElips2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
    let o = e.pos.point;
    let a = e.major_radius;
    let b = e.minor_radius;
    let xd = e.pos.vxdir;
    let yd = e.pos.vydir;
    let v = GpVec2d::new(p.x() - o.x(), p.y() - o.y());
    let x = v.dot(&GpVec2d::new(xd.x, xd.y));
    let y = v.dot(&GpVec2d::new(yd.x, yd.y));
    if v.magnitude() < CONFUSION && (a - b).abs() < CONFUSION {
        return Vec::new();
    }
    // `F(u) = (P - E(u))·E'(u) = (b² - a²)·cos·sin - b·y·cos + a·x·sin = 0`
    // (`Extrema_ExtPElC2d.cxx:192`, `math_TrigonometricFunctionRoots(0,
    // (b² - a²)/2, -b·y, a·x, 0, Uinf, Usup)`), in the
    // `d·cos·sin + e·cos + f·sin` form `trig_roots_sincos` solves.
    let us = trig_roots_sincos(b * b - a * a, -b * y, a * x, uinf, usup);
    us.into_iter().map(|u| pair2d(p, u, ellipse2d_val(e, u))).collect()
}

pub(super) fn hyperbola_all2d(h: &GpHypr2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
    let o = h.pos.point;
    let r = h.major_radius;
    let s = h.minor_radius;
    let xd = h.pos.vxdir;
    let yd = h.pos.vydir;
    let v = GpVec2d::new(p.x() - o.x(), p.y() - o.y());
    let x = v.dot(&GpVec2d::new(xd.x, xd.y));
    let y = v.dot(&GpVec2d::new(yd.x, yd.y));
    let c1 = (r * r + s * s) / 4.0;
    let vs = quartic_roots(c1, -(x * r + y * s) / 2.0, 0.0, (x * r - y * s) / 2.0, -c1);
    let tol2 = CONFUSION * CONFUSION;
    let mut out = Vec::new();
    let mut seen: Vec<GpPnt2d> = Vec::new();
    for vv in vs {
        if vv > 0.0 {
            let u = vv.ln();
            if u >= uinf && u <= usup {
                let q = hyperbola2d_val(h, u);
                if seen.iter().all(|t| t.square_distance(&q) >= tol2) {
                    seen.push(q);
                    out.push(pair2d(p, u, q));
                }
            }
        }
    }
    out
}

pub(super) fn parabola_all2d(pa: &GpParab2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
    let o = pa.pos.point;
    let f = pa.focal;
    if f.abs() < RESOLUTION {
        return Vec::new();
    }
    let xd = pa.pos.vxdir;
    let yd = pa.pos.vydir;
    let v = GpVec2d::new(p.x() - o.x(), p.y() - o.y());
    let x = v.dot(&GpVec2d::new(xd.x, xd.y));
    let y = v.dot(&GpVec2d::new(yd.x, yd.y));
    let us = cubic_roots(1.0 / (4.0 * f), 0.0, 2.0 * f - x, -2.0 * f * y);
    let tol2 = CONFUSION * CONFUSION;
    let mut out = Vec::new();
    let mut seen: Vec<GpPnt2d> = Vec::new();
    for u in us {
        if u >= uinf && u <= usup {
            let q = parabola2d_val(pa, u);
            if seen.iter().all(|t| t.square_distance(&q) >= tol2) {
                seen.push(q);
                out.push(pair2d(p, u, q));
            }
        }
    }
    out
}

/// Public: closest point on the 2D line.
pub fn point_line_extrema2d(l: &GpLin2d, p: &GpPnt2d) -> Extrema2d {
    let inf = f64::INFINITY;
    line_all2d(l, p, -inf, inf).first().copied().unwrap_or_else(|| {
        Extrema2d { p1: *p, p2: l.pos.loc, distance: p.distance(&l.pos.loc), u1: 0.0, u2: 0.0 }
    })
}

/// Public: closest point on the 2D circle (full range [0, 2π]).
pub fn point_circle_extrema2d(c: &GpCirc2d, p: &GpPnt2d) -> Extrema2d {
    let period = 2.0 * std::f64::consts::PI;
    circle_all2d(c, p, 0.0, period)
        .into_iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal))
        .unwrap_or_else(|| {
            let q = circle2d_val(c, 0.0);
            Extrema2d { p1: *p, p2: q, distance: p.distance(&q), u1: 0.0, u2: 0.0 }
        })
}

/// Public: closest point on the 2D ellipse (full range [0, 2π]).
pub fn point_ellipse_extrema2d(e: &GpElips2d, p: &GpPnt2d) -> Extrema2d {
    let period = 2.0 * std::f64::consts::PI;
    ellipse_all2d(e, p, 0.0, period)
        .into_iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal))
        .unwrap_or_else(|| {
            let q = ellipse2d_val(e, 0.0);
            Extrema2d { p1: *p, p2: q, distance: p.distance(&q), u1: 0.0, u2: 0.0 }
        })
}

/// Public: closest point on the 2D hyperbola.
pub fn point_hyperbola_extrema2d(h: &GpHypr2d, p: &GpPnt2d) -> Extrema2d {
    let inf = f64::INFINITY;
    hyperbola_all2d(h, p, -inf, inf)
        .into_iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal))
        .unwrap_or_else(|| {
            let q = hyperbola2d_val(h, 0.0);
            Extrema2d { p1: *p, p2: q, distance: p.distance(&q), u1: 0.0, u2: 0.0 }
        })
}

/// Public: closest point on the 2D parabola.
pub fn point_parabola_extrema2d(pa: &GpParab2d, p: &GpPnt2d) -> Extrema2d {
    let inf = f64::INFINITY;
    parabola_all2d(pa, p, -inf, inf)
        .into_iter()
        .min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal))
        .unwrap_or_else(|| {
            let q = parabola2d_val(pa, 0.0);
            Extrema2d { p1: *p, p2: q, distance: p.distance(&q), u1: 0.0, u2: 0.0 }
        })
}

// ---------------------------------------------------------------------------
// Classification + general Newton path for a `dyn Curve2d`.
// ---------------------------------------------------------------------------

pub(super) fn is_line2d(c: &dyn Curve2d) -> bool {
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
        if t0.crossed(&t).abs() > 1e-7 * m0 * m {
            return false;
        }
    }
    true
}

pub(super) fn circumcenter2d(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d) -> Option<GpPnt2d> {
    let d = 2.0 * (a.x() * (b.y() - c.y()) + b.x() * (c.y() - a.y()) + c.x() * (a.y() - b.y()));
    if d.abs() < 1e-20 {
        return None;
    }
    let a1 = a.x() * a.x() + a.y() * a.y();
    let b1 = b.x() * b.x() + b.y() * b.y();
    let c1 = c.x() * c.x() + c.y() * c.y();
    let ux = (a1 * (b.y() - c.y()) + b1 * (c.y() - a.y()) + c1 * (a.y() - b.y())) / d;
    let uy = (a1 * (c.x() - b.x()) + b1 * (a.x() - c.x()) + c1 * (b.x() - a.x())) / d;
    Some(GpPnt2d::new(ux, uy))
}

pub(super) fn classify_circle2d(c: &dyn Curve2d) -> Option<GpCirc2d> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite() && b > a) {
        return None;
    }
    let n = 8;
    let pts: Vec<GpPnt2d> = (0..n).map(|i| c.d0(a + (b - a) * i as f64 / (n - 1) as f64)).collect();
    let o = circumcenter2d(&pts[0], &pts[1], &pts[2])?;
    let r = pts[0].distance(&o);
    if r < 1e-12 {
        return None;
    }
    for p in &pts {
        if (p.distance(&o) - r).abs() > 1e-4 * r {
            return None;
        }
    }
    let xd = GpDir2d::from_vec2d(&GpVec2d::new(pts[0].x() - o.x(), pts[0].y() - o.y())).ok()?;
    Some(GpCirc2d { pos: GpAx22d::from_xdir(o, xd), radius: r })
}

pub(super) fn fval2d(c: &dyn Curve2d, p: &GpPnt2d, u: f64) -> Option<f64> {
    let (q, d1) = c.d1(u);
    let m = d1.magnitude();
    if m < 1e-12 || !q.x().is_finite() {
        return None;
    }
    Some(GpVec2d::new(q.x() - p.x(), q.y() - p.y()).dot(&d1))
}

pub(super) fn fprime2d(c: &dyn Curve2d, p: &GpPnt2d, u: f64) -> Option<f64> {
    let (q, d1, d2) = c.d2(u);
    let m = d1.magnitude();
    if m < 1e-12 {
        return None;
    }
    let v = GpVec2d::new(q.x() - p.x(), q.y() - p.y());
    Some(d1.square_magnitude() + v.dot(&d2))
}

pub(super) fn dist22d(c: &dyn Curve2d, p: &GpPnt2d, u: f64) -> Option<f64> {
    let q = c.d0(u);
    if !q.x().is_finite() {
        return None;
    }
    Some(q.square_distance(p))
}

pub(super) fn solve_f_zero2d(c: &dyn Curve2d, p: &GpPnt2d, mut lo: f64, mut hi: f64) -> f64 {
    let mut u = 0.5 * (lo + hi);
    for _ in 0..60 {
        let fl = match (fval2d(c, p, lo), fval2d(c, p, hi)) {
            (Some(fl), Some(_)) => fl,
            _ => break,
        };
        let un = match (fval2d(c, p, u), fprime2d(c, p, u)) {
            (Some(f), Some(fp)) if fp.abs() > 1e-300 => {
                let s = f / fp;
                if s.is_finite() && u - s > lo && u - s < hi { u - s } else { 0.5 * (lo + hi) }
            }
            _ => 0.5 * (lo + hi),
        };
        let fu = match fval2d(c, p, un) {
            Some(v) => v,
            None => break,
        };
        if fl * fu <= 0.0 {
            hi = un;
        } else {
            lo = un;
        }
        if (hi - lo).abs() < 1e-13 * (1.0 + hi.abs()) {
            break;
        }
        u = un;
    }
    0.5 * (lo + hi)
}

pub(super) fn refine_seed2d(c: &dyn Curve2d, p: &GpPnt2d, mut u: f64, lo: f64, hi: f64) -> f64 {
    for _ in 0..24 {
        match (fval2d(c, p, u), fprime2d(c, p, u)) {
            (Some(f), Some(fp)) if fp.abs() > 1e-300 => {
                let un = u - f / fp;
                if !un.is_finite() {
                    break;
                }
                if (un - u).abs() < 1e-12 * (1.0 + u.abs()) {
                    u = un;
                    break;
                }
                u = un.clamp(lo, hi);
            }
            _ => break,
        }
    }
    u
}

pub(super) fn build_samples2d(_c: &dyn Curve2d, a: f64, b: f64) -> Vec<f64> {
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

pub(super) fn pair2d_c(c: &dyn Curve2d, p: &GpPnt2d, u: f64) -> Extrema2d {
    pair2d(p, u, c.d0(u))
}

/// All local extrema of |C(u)-P| via grid → bracket sign changes → Newton.
pub(super) fn newton_point_curve_all2d(c: &dyn Curve2d, p: &GpPnt2d) -> Vec<Extrema2d> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let samples = build_samples2d(c, a, b);
    let mut pairs: Vec<Extrema2d> = Vec::new();
    for w in samples.windows(2) {
        let (u0, u1) = (w[0], w[1]);
        if (u1 - u0).abs() < 1e-14 {
            continue;
        }
        if let (Some(f0), Some(f1)) = (fval2d(c, p, u0), fval2d(c, p, u1)) {
            if f0 * f1 <= 0.0 {
                pairs.push(pair2d_c(c, p, solve_f_zero2d(c, p, u0, u1)));
            }
        }
    }
    for i in 1..samples.len().saturating_sub(1) {
        let (u0, u1, u2) = (samples[i - 1], samples[i], samples[i + 1]);
        let (d0, d1, d2) = match (dist22d(c, p, u0), dist22d(c, p, u1), dist22d(c, p, u2)) {
            (Some(a), Some(b), Some(cc)) => (a, b, cc),
            _ => continue,
        };
        if (d1 <= d0 && d1 <= d2) || (d1 >= d0 && d1 >= d2) {
            pairs.push(pair2d_c(c, p, refine_seed2d(c, p, u1, u0, u2)));
        }
    }
    if a.is_finite() {
        pairs.push(pair2d_c(c, p, a));
    }
    if b.is_finite() {
        pairs.push(pair2d_c(c, p, b));
    }
    pairs
}

pub(super) fn dedupe_sort2d(v: Vec<Extrema2d>) -> Vec<Extrema2d> {
    let mut out: Vec<Extrema2d> = Vec::new();
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
