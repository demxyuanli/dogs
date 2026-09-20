//! Curve/surface measurement utilities — lengths, areas, curvatures.
//!
//! **UNPORTED (audit A10 §8 / A15)**: OCCT has no such functions in
//! `ElCLib`/`ElSLib`. The self-contained Simpson integration, chord-length
//! accumulation (`curve_arc_length`, `point_at_arc_length`) and midpoint
//! rectangle areas (`surface_patch_area`) below are port-local; OCCT's
//! equivalents are `GCPnts_AbscissaPoint` (arc length via Newton on the
//! integral, with a tolerance) and `GProp_*` (areas/volumes).
use crate::gp::{GpPnt, GpVec, GpLin, GpCirc, GpElips, GpHypr, GpParab, GpPln};
use crate::elib::{clib, slib};

/// Exact line length between parameters (infinite line — length unbounded).
pub fn line_length(_l: &GpLin, a: f64, b: f64) -> f64 {
    // |b - a| * unit direction = parameter span
    (b - a).abs()
}

/// Exact circle arc length over [a, b] (in radians).
pub fn circle_arc_length(c: &GpCirc, a: f64, b: f64) -> f64 {
    c.radius * (b - a).abs()
}

/// Exact circle circumference.
pub fn circle_circumference(c: &GpCirc) -> f64 { 2.0 * std::f64::consts::PI * c.radius }

/// Ellipse arc length — no closed form; Gauss-Legendre quadrature.
pub fn ellipse_arc_length(e: &GpElips, a: f64, b: f64, n: usize) -> f64 {
    let major = e.major_radius;
    let minor = e.minor_radius;
    // Speed: ds/du = sqrt((major*sin u)² + (minor*cos u)²)
    let integrand = |u: f64| {
        let su = u.sin(); let cu = u.cos();
        (major*major*su*su + minor*minor*cu*cu).sqrt()
    };
    // Simpson's rule (robust, no dependency)
    simpson(&integrand, a, b, n)
}

/// Parabola arc length over [a, b]. y = u²/(4f) along X, u along Y.
pub fn parabola_arc_length(p: &GpParab, a: f64, b: f64, n: usize) -> f64 {
    let f = p.focal;
    let integrand = |u: f64| {
        (1.0 + (u / (2.0 * f)) * (u / (2.0 * f))).sqrt()
    };
    simpson(&integrand, a, b, n)
}

/// Hyperbola arc length via quadrature.
pub fn hyperbola_arc_length(h: &GpHypr, a: f64, b: f64, n: usize) -> f64 {
    let major = h.major_radius; let minor = h.minor_radius;
    let integrand = |u: f64| {
        let su = u.sinh(); let cu = u.cosh();
        (major*major*cu*cu + minor*minor*su*su).sqrt()
    };
    simpson(&integrand, a, b, n)
}

/// Simpson's rule numerical integration.
pub fn simpson<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> f64 {
    let n = (n / 2 * 2).max(2); // even
    let h = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for i in 1..n {
        let x = a + i as f64 * h;
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(x);
    }
    sum * h / 3.0
}

/// Numeric arc length of any parametric curve.
pub fn curve_arc_length<F: Fn(f64) -> GpPnt>(pt: &F, a: f64, b: f64, n: usize) -> f64 {
    let n = n.max(2);
    let h = (b - a) / n as f64;
    let mut prev = pt(a);
    let mut total = 0.0;
    for i in 1..=n {
        let x = a + i as f64 * h;
        let cur = pt(x);
        total += cur.coord.subtracted(&prev.coord).modulus();
        prev = cur;
    }
    total
}

/// Plane area (unbounded — returns infinity).
pub fn plane_area(_p: &GpPln) -> f64 { f64::INFINITY }

/// Circle disc area.
pub fn circle_disc_area(c: &GpCirc) -> f64 { std::f64::consts::PI * c.radius * c.radius }

/// Ellipse area.
pub fn ellipse_area(e: &GpElips) -> f64 { std::f64::consts::PI * e.major_radius * e.minor_radius }

/// Circle curvature (constant = 1/r).
pub fn circle_curvature(c: &GpCirc) -> f64 {
    if c.radius.abs() > 1e-30 { 1.0 / c.radius.abs() } else { f64::INFINITY }
}

/// Ellipse curvature at parameter u.
pub fn ellipse_curvature(e: &GpElips, u: f64) -> f64 {
    let a = e.major_radius; let b = e.minor_radius;
    let su = u.sin(); let cu = u.cos();
    let denom = (a*a*su*su + b*b*cu*cu).powf(1.5);
    if denom.abs() < 1e-30 { return 0.0; }
    a * b / denom
}

/// Line curvature = 0.
pub fn line_curvature(_l: &GpLin) -> f64 { 0.0 }

/// Normalize parameter to circle's [0, 2π) period.
pub fn normalize_circle_param(u: f64) -> f64 {
    u.rem_euclid(2.0 * std::f64::consts::PI)
}

/// Evaluate point at a given arc-length offset along a sampled curve.
/// Returns interpolated point.
pub fn point_at_arc_length<F: Fn(f64) -> GpPnt>(pt: &F, a: f64, b: f64, target_len: f64, n: usize) -> GpPnt {
    // Walk samples accumulating length until >= target
    let n = n.max(4);
    let h = (b - a) / n as f64;
    let mut prev = pt(a);
    let mut acc = 0.0;
    for i in 1..=n {
        let x = a + i as f64 * h;
        let cur = pt(x);
        let seg = cur.coord.subtracted(&prev.coord).modulus();
        if acc + seg >= target_len {
            let frac = if seg > 1e-30 { (target_len - acc) / seg } else { 0.0 };
            let dx = cur.coord.subtracted(&prev.coord).multiplied(frac);
            return GpPnt::from_xyz(&prev.coord.added(&dx));
        }
        acc += seg;
        prev = cur;
    }
    prev
}

/// Surface patch area via numerical integration of |du × dv|.
pub fn surface_patch_area<F: Fn(f64, f64) -> (GpPnt, GpVec, GpVec)>(
    d1: &F, u0: f64, u1: f64, v0: f64, v1: f64, nu: usize, nv: usize) -> f64 {
    let (nu, nv) = (nu.max(2), nv.max(2));
    let hu = (u1 - u0) / nu as f64;
    let hv = (v1 - v0) / nv as f64;
    let mut area = 0.0;
    for i in 0..nu {
        for j in 0..nv {
            let u = u0 + (i as f64 + 0.5) * hu;
            let v = v0 + (j as f64 + 0.5) * hv;
            let (_, du, dv) = d1(u, v);
            let n = du.xyz().crossed(dv.xyz());
            area += n.modulus() * hu * hv;
        }
    }
    area
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::{GpDir, GpAx2};

    #[test]
    fn circle_arc() {
        let c = GpCirc::new(GpAx2::standard(), 2.0);
        let len = circle_arc_length(&c, 0.0, std::f64::consts::FRAC_PI_2);
        assert!((len - std::f64::consts::PI).abs() < 1e-12); // r*θ = 2*π/2 = π
    }

    #[test]
    fn simpson_sin() {
        let v = simpson(&|x| x.sin(), 0.0, std::f64::consts::PI, 100);
        assert!((v - 2.0).abs() < 1e-6, "v={v}");
    }

    #[test]
    fn ellipse_area_basic() {
        let e = GpElips::new(GpAx2::standard(), 2.0, 1.0);
        assert!((ellipse_area(&e) - 2.0 * std::f64::consts::PI).abs() < 1e-12);
    }

    #[test]
    fn ellipse_curvature_extrema() {
        // At u=0 (major end): curvature = a/b² = 2/1 = 2
        let e = GpElips::new(GpAx2::standard(), 2.0, 1.0);
        let k = ellipse_curvature(&e, 0.0);
        assert!((k - 2.0).abs() < 1e-10, "k={k}");
    }

    #[test]
    fn test_circle_curvature() {
        let c = GpCirc::new(GpAx2::standard(), 5.0);
        assert!((circle_curvature(&c) - 0.2).abs() < 1e-12);
    }
}
