//! 2D extrema — point/curve minimum and maximum distances in the plane.
//!
//! Port of `Extrema_ExtPC2d`, `Extrema_ExtCC2d` (discretization + local
//! refinement). Source: `Extrema` (TKGeomBase).

use std::cmp::Ordering;

use occt_core::gp::{
    GpAx22d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d,
};
use occt_core::precision::{ANGULAR, CONFUSION, RESOLUTION};

use crate::curve::Curve2d;

const PI: f64 = std::f64::consts::PI;

/// A solved extremum between two 2D objects.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Extrema2d {
    pub p1: GpPnt2d,
    pub p2: GpPnt2d,
    pub distance: f64,
    pub u1: f64,
    pub u2: f64,
}

fn bound(c: &dyn Curve2d) -> (f64, f64) {
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

fn quadratic_roots(a: f64, b: f64, c: f64) -> Vec<f64> {
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

fn cubic_roots(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
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

fn poly4(a: f64, b: f64, c: f64, d: f64, e: f64, x: f64) -> f64 {
    ((a * x + b) * x + c) * x * x + d * x + e
}

fn polish4(a: f64, b: f64, c: f64, d: f64, e: f64, mut x: f64) -> f64 {
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
fn quartic_roots(a4: f64, a3: f64, a2: f64, a1: f64, a0: f64) -> Vec<f64> {
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

fn push_periodic(roots: &mut Vec<f64>, u: f64, uinf: f64, usup: f64) {
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

fn polish_trig(d: f64, e: f64, f: f64, mut u: f64) -> f64 {
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
fn trig_roots_sincos(d: f64, e: f64, f: f64, uinf: f64, usup: f64) -> Vec<f64> {
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

fn adjust_periodic(uinf: f64, period: f64, tol: f64, u: f64) -> f64 {
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
fn circle2d_val(c: &GpCirc2d, u: f64) -> GpPnt2d {
    let r = c.radius;
    let o = c.location();
    let xd = c.pos.vxdir;
    let yd = c.pos.vydir;
    GpPnt2d::new(o.x() + r * u.cos() * xd.x + r * u.sin() * yd.x, o.y() + r * u.cos() * xd.y + r * u.sin() * yd.y)
}

/// Point on a 2D ellipse: O + a·cos·xd − b·sin·yd (matches `Geom2d_Ellipse`).
fn ellipse2d_val(e: &GpElips2d, u: f64) -> GpPnt2d {
    let a = e.major_radius;
    let b = e.minor_radius;
    let o = e.pos.point;
    let xd = e.pos.vxdir;
    let yd = e.pos.vydir;
    GpPnt2d::new(o.x() + a * u.cos() * xd.x - b * u.sin() * yd.x, o.y() + a * u.cos() * xd.y - b * u.sin() * yd.y)
}

/// Point on a 2D hyperbola: O + a·cosh·xd + b·sinh·yd.
fn hyperbola2d_val(h: &GpHypr2d, u: f64) -> GpPnt2d {
    let a = h.major_radius;
    let b = h.minor_radius;
    let o = h.pos.point;
    let xd = h.pos.vxdir;
    let yd = h.pos.vydir;
    GpPnt2d::new(o.x() + a * u.cosh() * xd.x + b * u.sinh() * yd.x, o.y() + a * u.cosh() * xd.y + b * u.sinh() * yd.y)
}

/// Point on a 2D parabola: O + (u²/4f)·xd + u·yd.
fn parabola2d_val(pa: &GpParab2d, u: f64) -> GpPnt2d {
    let f = pa.focal;
    let o = pa.pos.point;
    let xd = pa.pos.vxdir;
    let yd = pa.pos.vydir;
    let x = u * u / (4.0 * f);
    GpPnt2d::new(o.x() + x * xd.x + u * yd.x, o.y() + x * xd.y + u * yd.y)
}

fn pair2d(p: &GpPnt2d, u: f64, q: GpPnt2d) -> Extrema2d {
    Extrema2d { p1: *p, p2: q, distance: p.distance(&q), u1: u, u2: u }
}

fn line_all2d(l: &GpLin2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
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

fn circle_all2d(c: &GpCirc2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
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

fn ellipse_all2d(e: &GpElips2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
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
    let us = trig_roots_sincos(b * b - a * a, b * y, a * x, uinf, usup);
    us.into_iter().map(|u| pair2d(p, u, ellipse2d_val(e, u))).collect()
}

fn hyperbola_all2d(h: &GpHypr2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
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

fn parabola_all2d(pa: &GpParab2d, p: &GpPnt2d, uinf: f64, usup: f64) -> Vec<Extrema2d> {
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

fn is_line2d(c: &dyn Curve2d) -> bool {
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

fn circumcenter2d(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d) -> Option<GpPnt2d> {
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

fn classify_circle2d(c: &dyn Curve2d) -> Option<GpCirc2d> {
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

fn fval2d(c: &dyn Curve2d, p: &GpPnt2d, u: f64) -> Option<f64> {
    let (q, d1) = c.d1(u);
    let m = d1.magnitude();
    if m < 1e-12 || !q.x().is_finite() {
        return None;
    }
    Some(GpVec2d::new(q.x() - p.x(), q.y() - p.y()).dot(&d1))
}

fn fprime2d(c: &dyn Curve2d, p: &GpPnt2d, u: f64) -> Option<f64> {
    let (q, d1, d2) = c.d2(u);
    let m = d1.magnitude();
    if m < 1e-12 {
        return None;
    }
    let v = GpVec2d::new(q.x() - p.x(), q.y() - p.y());
    Some(d1.square_magnitude() + v.dot(&d2))
}

fn dist22d(c: &dyn Curve2d, p: &GpPnt2d, u: f64) -> Option<f64> {
    let q = c.d0(u);
    if !q.x().is_finite() {
        return None;
    }
    Some(q.square_distance(p))
}

fn solve_f_zero2d(c: &dyn Curve2d, p: &GpPnt2d, mut lo: f64, mut hi: f64) -> f64 {
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

fn refine_seed2d(c: &dyn Curve2d, p: &GpPnt2d, mut u: f64, lo: f64, hi: f64) -> f64 {
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

fn build_samples2d(_c: &dyn Curve2d, a: f64, b: f64) -> Vec<f64> {
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

fn pair2d_c(c: &dyn Curve2d, p: &GpPnt2d, u: f64) -> Extrema2d {
    pair2d(p, u, c.d0(u))
}

/// All local extrema of |C(u)-P| via grid → bracket sign changes → Newton.
fn newton_point_curve_all2d(c: &dyn Curve2d, p: &GpPnt2d) -> Vec<Extrema2d> {
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

fn dedupe_sort2d(v: Vec<Extrema2d>) -> Vec<Extrema2d> {
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

fn param_for_point2d(c: &dyn Curve2d, q: &GpPnt2d, seed: f64, a: f64, b: f64) -> f64 {
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

fn pair2d_cc(p1: GpPnt2d, u1: f64, p2: GpPnt2d, u2: f64) -> Extrema2d {
    Extrema2d { p1, p2, distance: p1.distance(&p2), u1, u2 }
}

fn line2d_value(l: &GpLin2d, u: f64) -> GpPnt2d {
    GpPnt2d::new(l.pos.loc.x() + u * l.pos.vdir.x, l.pos.loc.y() + u * l.pos.vdir.y)
}

fn circle_param2d(c: &GpCirc2d, p: &GpPnt2d) -> f64 {
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

fn line_of_curve2d(c: &dyn Curve2d) -> Option<GpLin2d> {
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
fn f_and_jac2d(
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
fn refine_curve_curve2d_newton(
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
fn newton_curve_curve_all2d(c1: &dyn Curve2d, c2: &dyn Curve2d) -> Vec<Extrema2d> {
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
fn dedupe_sort_cc2d(v: Vec<Extrema2d>) -> Vec<Extrema2d> {
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

fn param_in2d(u: f64, lo: f64, hi: f64) -> bool {
    let lo_ok = !lo.is_finite() || u >= lo - CONFUSION;
    let hi_ok = !hi.is_finite() || u <= hi + CONFUSION;
    lo_ok && hi_ok
}

/// Map `u` into `[lo, lo+period)` for a periodic curve (OCCT `ElCLib::InPeriod`).
fn in_period2d(u: f64, lo: f64, period: f64) -> f64 {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Geom2dCircle, Geom2dEllipse, Geom2dLine};
    use occt_core::gp::{GpAx22d, GpAx2d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d};

    fn line2d() -> Geom2dLine {
        Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap()))
    }

    fn circle2d() -> Geom2dCircle {
        Geom2dCircle::new(GpCirc2d::new(occt_core::gp::GpAx22d::standard(), 1.0))
    }

    #[test]
    fn point_line_min_distance() {
        let c = line2d();
        let e = point_curve_extrema2d(&c, &GpPnt2d::new(3.0, 4.0));
        assert!((e.distance - 4.0).abs() < 1e-6, "dist {}", e.distance);
        assert!((e.p2.x() - 3.0).abs() < 1e-5, "closest x {}", e.p2.x());
    }

    #[test]
    fn point_circle_min_and_max() {
        let c = circle2d();
        let p = GpPnt2d::new(3.0, 0.0);
        let e = point_curve_extrema2d(&c, &p);
        assert!((e.distance - 2.0).abs() < 1e-6, "min {}", e.distance);
        let m = point_curve_max_extrema2d(&c, &p, 64);
        assert!((m.distance - 4.0).abs() < 1e-6, "max {}", m.distance);
    }

    #[test]
    fn curve_curve_circles_min() {
        // Two unit circles centers 3 apart → min distance 1.
        let c1 = circle2d();
        let c2 = Geom2dCircle::new(GpCirc2d::new(occt_core::gp::GpAx22d::standard(), 1.0).translated_vec(&GpVec2d::new(3.0, 0.0)));
        let es = curve_curve_extrema2d(&c1, &c2, 24);
        assert!(!es.is_empty());
        assert!((es[0].distance - 1.0).abs() < 1e-5, "min {}", es[0].distance);
    }

    #[test]
    fn curve_curve_intersections_found() {
        // Circle at origin and horizontal line y=0 → 2 intersection points.
        let c1 = circle2d();
        let l = line2d();
        let pts = curve_curve_intersections2d(&c1, &l, 1e-3);
        assert_eq!(pts.len(), 2, "circle∩x-axis points {pts:?}");
        for p in &pts {
            assert!((p.x().abs() - 1.0).abs() < 0.05, "on circle x {}", p.x());
        }
    }

    #[test]
    fn point_polyline_closest() {
        let poly = vec![GpPnt2d::new(0.0, 0.0), GpPnt2d::new(2.0, 0.0)];
        let e = point_polyline_extrema2d(&poly, &GpPnt2d::new(1.0, 3.0));
        assert!((e.distance - 3.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.p2.x() - 1.0).abs() < 1e-9 && e.p2.y().abs() < 1e-9);
    }

    #[test]
    fn point_polyline_endpoint_t() {
        let poly = vec![GpPnt2d::new(0.0, 0.0), GpPnt2d::new(2.0, 0.0)];
        // Point beyond the segment end clamps to the endpoint.
        let e = point_polyline_extrema2d(&poly, &GpPnt2d::new(5.0, 1.0));
        assert!((e.p2.x() - 2.0).abs() < 1e-9, "clamp x {}", e.p2.x());
    }

    #[test]
    fn tangent_horizontal_line() {
        let c = line2d();
        let t = tangent2d(&c, 0.5);
        assert!((t.x() - 1.0).abs() < 1e-9 && t.y().abs() < 1e-9, "tangent {t:?}");
    }

    #[test]
    fn skew_lines_min_distance() {
        // Horizontal line and vertical line offset — skew in the plane means
        // they cross; min distance 0. Use parallel lines instead: y=0 and y=3.
        let l1 = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap()));
        let l2 = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap()));
        let es = curve_curve_extrema2d(&l1, &l2, 8);
        assert!((es[0].distance - 3.0).abs() < 1e-5, "min {}", es[0].distance);
    }

    const PI: f64 = std::f64::consts::PI;

    #[test]
    fn analytic_line_extrema2d() {
        let l = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap());
        let e = point_line_extrema2d(&l, &GpPnt2d::new(3.0, 4.0));
        assert!((e.distance - 4.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.u1 - 3.0).abs() < 1e-9, "u {}", e.u1);
        assert!((e.p2.x() - 3.0).abs() < 1e-9 && e.p2.y().abs() < 1e-9);
    }

    #[test]
    fn analytic_circle_extrema2d_min_max() {
        let c = GpCirc2d::new(GpAx22d::standard(), 1.0);
        let p = GpPnt2d::new(3.0, 0.0);
        let all = circle_all2d(&c, &p, 0.0, 2.0 * PI);
        assert_eq!(all.len(), 2, "min+max {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.u1).abs() < 1e-9, "min u {}", min.u1);
        assert!((max.distance - 4.0).abs() < 1e-9, "max {}", max.distance);
        assert!((max.u1 - PI).abs() < 1e-9, "max u {}", max.u1);
    }

    #[test]
    fn analytic_ellipse_extrema2d() {
        // a=2, b=1, point (3,0): closest (2,0) at u=0 (dist 1), farthest (-2,0) at u=π.
        let e = GpElips2d::new(GpAx22d::standard(), 2.0, 1.0);
        let p = GpPnt2d::new(3.0, 0.0);
        let em = point_ellipse_extrema2d(&e, &p);
        assert!((em.distance - 1.0).abs() < 1e-7, "min {}", em.distance);
        assert!((em.p2.x() - 2.0).abs() < 1e-6 && em.p2.y().abs() < 1e-6, "closest {em:?}");
        let all = ellipse_all2d(&e, &p, 0.0, 2.0 * PI);
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((max.distance - 5.0).abs() < 1e-7, "max {}", max.distance);
    }

    #[test]
    fn analytic_hyperbola_extrema2d_origin() {
        let h = GpHypr2d::new(GpAx22d::standard(), 1.0, 1.0);
        let e = point_hyperbola_extrema2d(&h, &GpPnt2d::new(0.0, 0.0));
        assert!((e.distance - 1.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.u1).abs() < 1e-9, "u {}", e.u1);
        assert!((e.p2.x() - 1.0).abs() < 1e-9, "closest {e:?}");
    }

    #[test]
    fn analytic_parabola_extrema2d() {
        // F=1, C(u) = (u²/4, u). Point (4,0): min √12 at u = ±2√2.
        let pa = GpParab2d::new(GpAx22d::standard(), 1.0);
        let e = point_parabola_extrema2d(&pa, &GpPnt2d::new(4.0, 0.0));
        assert!((e.distance - (12.0f64).sqrt()).abs() < 1e-7, "dist {}", e.distance);
        assert!((e.u1.abs() - 2.0 * (2.0f64).sqrt()).abs() < 1e-6, "u {}", e.u1);
    }

    #[test]
    fn point_curve_extrema2d_all_line_and_circle() {
        let line = line2d();
        let all = point_curve_extrema2d_all(&line, &GpPnt2d::new(3.0, 4.0));
        assert!(!all.is_empty());
        assert!((all[0].distance - 4.0).abs() < 1e-9, "line min {}", all[0].distance);

        let circle = circle2d();
        let p = GpPnt2d::new(3.0, 0.0);
        let all = point_curve_extrema2d_all(&circle, &p);
        assert!((all[0].distance - 2.0).abs() < 1e-6, "circle min {}", all[0].distance);
        assert!((all[all.len() - 1].distance - 4.0).abs() < 1e-6, "circle max {}", all[all.len() - 1].distance);
    }

    #[test]
    fn newton_path_2d_bspline_min() {
        // Convex degree-2 arc: the Newton path must return a genuine minimum.
        let c = crate::bspline_curve::Geom2dBSplineCurve::new(
            vec![0.0, 2.0, 4.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let p = GpPnt2d::new(2.0, 4.0);
        let all = point_curve_extrema2d_all(&c, &p);
        assert!(!all.is_empty());
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let (_, d1) = c.d1(min.u1);
        let v = GpVec2d::new(min.p2.x() - p.x(), min.p2.y() - p.y());
        assert!(v.dot(&d1).abs() < 1e-6, "dF={}", v.dot(&d1));
        // The wrapper must agree.
        let e = point_curve_extrema2d(&c, &p);
        assert!((e.distance - min.distance).abs() < 1e-6, "wrapper {}", e.distance);
    }

    #[test]
    fn conic_via_dyn_curve_matches_analytic() {
        // A Geom2dEllipse driven through the dyn path: classification falls to
        // Newton (ellipse is not classified) and must still find the min.
        let e = Geom2dEllipse::from_axes(GpAx22d::standard(), 2.0, 1.0);
        let p = GpPnt2d::new(3.0, 0.0);
        let em = point_curve_extrema2d(&e, &p);
        assert!((em.distance - 1.0).abs() < 1e-6, "ellipse min {}", em.distance);
        assert!((em.p2.x() - 2.0).abs() < 1e-5, "closest {em:?}");
    }

    // --- Curve-curve extrema (Phase 13 Wave 2) ---

    #[test]
    fn line_line_extrema2d_parallel_constant_distance() {
        let l1 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap());
        let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap());
        let all = line_line_extrema2d(&l1, &l2);
        assert_eq!(all.len(), 1);
        assert!((all[0].distance - 3.0).abs() < 1e-9, "dist {}", all[0].distance);
    }

    #[test]
    fn line_line_extrema2d_intersecting_zero() {
        // x-axis and y-axis cross at the origin.
        let l1 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap());
        let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(0.0, 1.0).unwrap());
        let all = line_line_extrema2d(&l1, &l2);
        assert_eq!(all.len(), 1);
        assert!(all[0].distance < 1e-9, "dist {}", all[0].distance);
        assert!((all[0].p1.x()).abs() < 1e-9 && (all[0].p1.y()).abs() < 1e-9, "{:?}", all[0]);
    }

    #[test]
    fn line_circle_extrema2d_min_max() {
        // Line y=3 and unit circle at origin: min 2, max 4.
        let l = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap());
        let c = GpCirc2d::new(GpAx22d::standard(), 1.0);
        let all = line_circle_extrema2d(&l, &c);
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((max.distance - 4.0).abs() < 1e-9, "max {}", max.distance);
    }

    #[test]
    fn line_circle_extrema2d_intersecting_zero() {
        // Line y=0 through the unit circle: two zero-distance intersection pairs.
        let l = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap());
        let c = GpCirc2d::new(GpAx22d::standard(), 1.0);
        let all = line_circle_extrema2d(&l, &c);
        assert!(all.iter().any(|e| e.distance < 1e-9), "no intersection {all:?}");
    }

    #[test]
    fn circle_circle_extrema2d_min_and_max() {
        // Unit circles at (0,0) and (3,0): min 1, max 5.
        let c1 = GpCirc2d::new(GpAx22d::standard(), 1.0);
        let c2 = GpCirc2d::new(GpAx22d::standard(), 1.0).translated_vec(&GpVec2d::new(3.0, 0.0));
        let all = circle_circle_extrema2d(&c1, &c2);
        assert!(all.iter().any(|e| (e.distance - 1.0).abs() < 1e-7), "min missing {all:?}");
        assert!(all.iter().any(|e| (e.distance - 5.0).abs() < 1e-7), "max missing {all:?}");
    }

    #[test]
    fn curve_curve_extrema2d_all_circles_min_max() {
        let c1 = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let c2 = Geom2dCircle::new(
            GpCirc2d::new(GpAx22d::standard(), 1.0).translated_vec(&GpVec2d::new(3.0, 0.0)),
        );
        let all = curve_curve_extrema2d_all(&c1, &c2);
        assert!(!all.is_empty());
        assert!((all[0].distance - 1.0).abs() < 1e-7, "min {}", all[0].distance);
        assert!((all[all.len() - 1].distance - 5.0).abs() < 1e-7, "max {}", all[all.len() - 1].distance);
        // The samples-taking wrapper delegates to the same set.
        let es = curve_curve_extrema2d(&c1, &c2, 16);
        assert!((es[0].distance - 1.0).abs() < 1e-7, "wrapper min {}", es[0].distance);
    }

    #[test]
    fn newton_path_2d_bspline_curve_curve() {
        // Convex degree-2 arc (max y = 1 at u = 0.5) vs line y=3: min distance 2.
        let line = Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap());
        let bs = crate::bspline_curve::Geom2dBSplineCurve::new(
            vec![0.0, 2.0, 4.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let all = curve_curve_extrema2d_all(&line, &bs);
        assert!(!all.is_empty(), "no extrema");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-5, "min {}", min.distance);
        assert!((min.u2 - 0.5).abs() < 1e-3, "bspline param {}", min.u2);
        // Extremum condition: (C1−C2)·C′ ≈ 0 on both curves.
        let v = GpVec2d::new(min.p1.x() - min.p2.x(), min.p1.y() - min.p2.y());
        let (_, d1l) = line.d1(min.u1);
        let (_, d1b) = bs.d1(min.u2);
        assert!(v.dot(&d1l).abs() < 1e-4, "dF1 {}", v.dot(&d1l));
        assert!(v.dot(&d1b).abs() < 1e-4, "dF2 {}", v.dot(&d1b));
    }

    #[test]
    fn locate_extcc2d_refines_seed() {
        let line = Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap());
        let bs = crate::bspline_curve::Geom2dBSplineCurve::new(
            vec![0.0, 2.0, 4.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let e = locate_extcc2d(&line, &bs, 2.0, 0.5).expect("locate failed");
        assert!((e.distance - 2.0).abs() < 1e-5, "dist {}", e.distance);
        assert!((e.u2 - 0.5).abs() < 1e-3, "param {}", e.u2);
    }
}
