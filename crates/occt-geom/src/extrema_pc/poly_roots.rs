use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Polynomial root helpers.
//
// These are the ported `math_DirectPolynomialRoots`
// (`occt-core/src/math_direct_poly_roots.rs`) consumed as the OCCT
// `Value(1 .. NbSolutions)` list. The closed-form quadratic/cubic, the Ferrari
// quartic and their `poly4`/`polish4` helpers that used to live here were a
// second, invented solver (with `1e-14`/`1e-12`/`1e-9`/`1e-8`/`1e-6`/`1e-300`
// thresholds of their own); the header note claiming `occt-geom` has no
// `occt-math` dependency was stale -- it depends on `occt-math` and
// `occt-core`.
// ---------------------------------------------------------------------------

pub(super) use occt_core::math_direct_poly_roots::{direct_cubic_roots as cubic_roots, direct_quadratic_roots as quadratic_roots, quartic_roots};

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
/// (B²-A²)·cos·sin − B·Y·cos + A·X·sin = 0 for u (`Extrema_ExtPElC.cxx:253-260`
/// with `ko3 = -B*Y`), i.e. `(C(u)-Pp)·C'(u) = 0` for
/// `C(u) = O + A·cos·X + B·sin·Y` (`ElCLib::EllipseValue`, `cxx:176-189`).
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
    let us = trig_roots_sincos(b * b - a * a, -b * y, a * x, uinf, usup);
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
