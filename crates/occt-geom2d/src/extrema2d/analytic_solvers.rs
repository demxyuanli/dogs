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

// ---------------------------------------------------------------------------
// Polynomial root helpers.
//
// These are the ported `math_DirectPolynomialRoots`
// (`occt-core/src/math_direct_poly_roots.rs`) consumed as the OCCT
// `Value(1 .. NbSolutions)` list. The local closed-form quadratic/cubic, the
// Ferrari quartic and their `poly4`/`polish4` helpers were a fourth copy of the
// same invented solver and are gone (`occt-geom2d` depends on `occt-core`,
// where the port lives, so nothing had to be kept local).
// ---------------------------------------------------------------------------

pub(super) use occt_core::math_direct_poly_roots::{direct_cubic_roots, quartic_roots};

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
    let us = direct_cubic_roots(1.0 / (4.0 * f), 0.0, 2.0 * f - x, -2.0 * f * y);
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

// The self-invented `span/0.1` seed grid, its bracket bisection and its Newton
// refiners used to live here (`build_samples2d`, `newton_point_curve_all2d`,
// `fval2d`, `fprime2d`, `dist22d`, `solve_f_zero2d`, `refine_seed2d`,
// `pair2d_c`, `dedupe_sort2d`). They stood in for `Extrema_ExtPC2d` =
// `Extrema_GGExtPC<Adaptor2d_Curve2d, ...>` and are deleted; the faithful
// engine now lives in `extrema2d::ext_pc2d` / `extrema2d::gfunc_ext_pc2d`.
