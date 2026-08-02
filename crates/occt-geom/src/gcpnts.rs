//! Curve discretization points. Port of the `GCPnts` package (TKGeomBase):
//! `GCPnts_AbscissaPoint`, `GCPnts_UniformAbscissa`, `GCPnts_QuasiUniformAbscissa`,
//! plus parameter accessors for the deflection samplers already implemented in
//! `occt-core` (`GCPnts_UniformDeflection`, `GCPnts_TangentialDeflection`).
//!
//! The `Curve` trait has no arc-length method, so lengths are integrated
//! adaptively (adaptive Simpson) over `|C′(u)|` and the abscissa inversion is
//! a safeguarded Newton solve, mirroring `CPnts_AbscissaPoint::Perform`.

use occt_core::gcpnts::{CurveDeriv, CurveSample, TangentialDeflection, UniformDeflection};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::CONFUSION;

use crate::curve::Curve;

const MAX_DEPTH: usize = 20;

fn speed(c: &dyn Curve, u: f64) -> f64 {
    let (_, d1) = c.d1(u);
    d1.magnitude()
}

/// Adaptive Simpson integration of `|C′(u)|` over `[a, b]` to absolute
/// tolerance `tol`.
fn adaptive_length(c: &dyn Curve, a: f64, b: f64, tol: f64) -> f64 {
    fn simpson(a: f64, b: f64, fa: f64, fm: f64, fb: f64) -> f64 {
        (b - a) / 6.0 * (fa + 4.0 * fm + fb)
    }
    fn rec(
        c: &dyn Curve,
        a: f64,
        b: f64,
        fa: f64,
        fm: f64,
        fb: f64,
        whole: f64,
        tol: f64,
        depth: usize,
    ) -> f64 {
        let mid = 0.5 * (a + b);
        let lm = 0.5 * (a + mid);
        let rm = 0.5 * (mid + b);
        let flm = speed(c, lm);
        let frm = speed(c, rm);
        let left = simpson(a, mid, fa, flm, fm);
        let right = simpson(mid, b, fm, frm, fb);
        if depth >= MAX_DEPTH || (left + right - whole).abs() <= 15.0 * tol {
            left + right + (left + right - whole) / 15.0
        } else {
            rec(c, a, mid, fa, flm, fm, left, tol * 0.5, depth + 1)
                + rec(c, mid, b, fm, frm, fb, right, tol * 0.5, depth + 1)
        }
    }
    let (fa, fm, fb) = (speed(c, a), speed(c, 0.5 * (a + b)), speed(c, b));
    let whole = simpson(a, b, fa, fm, fb);
    rec(c, a, b, fa, fm, fb, whole, tol, 0)
}

/// Arc length of `c` over `[a, b]` (absolute tolerance `tol`).
pub fn curve_length_range(c: &dyn Curve, a: f64, b: f64, tol: f64) -> f64 {
    if !(a.is_finite() && b.is_finite()) {
        return f64::INFINITY;
    }
    if b <= a {
        return 0.0;
    }
    let span = b - a;
    // Scale the tolerance to the span so absolute small curves still resolve.
    let tol = tol.max(span * 1e-10);
    adaptive_length(c, a, b, tol)
}

/// Total arc length of `c` over its parameter range.
pub fn curve_length(c: &dyn Curve) -> f64 {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() {
        curve_length_range(c, a, b, CONFUSION * 0.1)
    } else {
        f64::INFINITY
    }
}

/// Parameter `u` at arc length `abscissa` from the point of parameter `from`.
///
/// Positive `abscissa` walks forward, negative backward. Mirrors
/// `GCPnts_AbscissaPoint`; errors if the requested distance exceeds the
/// available curve length.
pub fn abscissa_point(c: &dyn Curve, abscissa: f64, from: f64) -> Result<f64, String> {
    if abscissa.abs() <= CONFUSION {
        return Ok(from);
    }
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return Err("abscissa_point: unbounded curve".to_string());
    }
    if from < a - 1e-12 || from > b + 1e-12 {
        return Err(format!("abscissa_point: from {from} outside [{a}, {b}]"));
    }
    let from = from.clamp(a, b);
    let tol = CONFUSION * 0.1;

    let (lo, hi): (f64, f64) = if abscissa > 0.0 {
        (from, b)
    } else {
        (a, from)
    };
    // g(u) is monotone increasing with the signed arc length from `from`.
    // For abscissa < 0 the walk is backward: g(u) = |abscissa| − arc_len(u, from).
    let g = |u: f64| -> f64 {
        if abscissa > 0.0 {
            curve_length_range(c, from, u, tol) - abscissa
        } else {
            -abscissa - curve_length_range(c, u, from, tol)
        }
    };
    let sp = |u: f64| speed(c, u);

    let glo = g(lo);
    let ghi = g(hi);
    if glo > 0.0 || ghi < 0.0 {
        return Err("abscissa_point: abscissa beyond curve length".to_string());
    }
    if (ghi - glo).abs() < 1e-15 {
        return Ok(0.5 * (lo + hi));
    }

    let mut lo = lo;
    let mut hi = hi;
    let mut u = 0.5 * (lo + hi);
    for _ in 0..80 {
        let gu = g(u);
        if gu.abs() < 1e-9 * (1.0 + u.abs()) {
            return Ok(u);
        }
        let s = sp(u);
        let step = if s > 1e-14 { gu / s } else { 0.0 };
        let un = if step.is_finite() && step.abs() <= hi - lo {
            u - step
        } else {
            0.5 * (lo + hi)
        };
        let un = un.clamp(lo, hi);
        if g(un) <= 0.0 {
            lo = un;
        } else {
            hi = un;
        }
        u = un;
        if hi - lo < 1e-12 * (1.0 + hi.abs()) {
            break;
        }
    }
    Ok(u)
}

/// `n + 1` parameters at equal arc-length spacing across the curve
/// (`GCPnts_UniformAbscissa`). A zero-length curve returns `[a; n + 1]`.
pub fn uniform_abscissa(c: &dyn Curve, n: usize) -> Result<Vec<f64>, String> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return Err("uniform_abscissa: unbounded curve".to_string());
    }
    let n = n.max(1);
    let total = curve_length_range(c, a, b, CONFUSION * 0.1);
    if total <= CONFUSION {
        // Degenerate: every parameter maps to the same point.
        return Ok(vec![a; n + 1]);
    }
    let step = total / n as f64;
    let mut params = Vec::with_capacity(n + 1);
    params.push(a);
    let mut prev = a;
    for _ in 1..n {
        // Find u with arc_length(prev, u) = step, in [prev, b].
        let tol = CONFUSION * 0.1;
        let g = |u: f64| curve_length_range(c, prev, u, tol) - step;
        let sp = |u: f64| speed(c, u);
        let mut lo = prev;
        let mut hi = b;
        let mut u = 0.5 * (lo + hi);
        let mut result = u;
        for _ in 0..80 {
            let gu = g(u);
            if gu.abs() < 1e-9 * (1.0 + u.abs()) {
                result = u;
                break;
            }
            let s = sp(u);
            let step_n = if s > 1e-14 { gu / s } else { 0.0 };
            let un = if step_n.is_finite() && step_n.abs() <= hi - lo {
                u - step_n
            } else {
                0.5 * (lo + hi)
            };
            let un = un.clamp(lo, hi);
            if g(un) <= 0.0 {
                lo = un;
            } else {
                hi = un;
            }
            u = un;
            result = u;
            if hi - lo < 1e-12 * (1.0 + hi.abs()) {
                break;
            }
        }
        let u = result.min(b);
        params.push(u);
        prev = u;
    }
    params.push(b);
    Ok(params)
}

/// `n` parameters distributed at equal *chord* intervals over a polyline
/// approximation of the curve (port of `GCPnts_QuasiUniformAbscissa`, which
/// uses a 2n-sample chord-length table).
pub fn quasi_uniform_abscissa(c: &dyn Curve, n: usize) -> Result<Vec<f64>, String> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return Err("quasi_uniform_abscissa: unbounded curve".to_string());
    }
    let n = n.max(2);
    let total = curve_length_range(c, a, b, CONFUSION * 0.1);
    if total <= CONFUSION {
        return Ok(vec![a; n]);
    }
    // Chord-length vs parameter table over 2n samples (last lands on b).
    let samples = 2 * n;
    let du = (b - a) / (samples - 1) as f64;
    let mut cum = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    cum.push(0.0);
    params.push(a);
    let mut prev = c.d0(a);
    for k in 1..samples {
        let u = a + k as f64 * du;
        let p = c.d0(u);
        let l = cum[k - 1] + prev.distance(&p);
        cum.push(l);
        params.push(u);
        prev = p;
    }
    let total = *cum.last().unwrap();
    let dcorde = total / (n - 1) as f64;
    let mut out = Vec::with_capacity(n);
    out.push(a);
    let mut idx = 1usize;
    for i in 1..n - 1 {
        let target = dcorde * i as f64;
        while idx < samples - 1 && cum[idx] < target {
            idx += 1;
        }
        let denom = cum[idx] - cum[idx - 1];
        let alpha = if denom.abs() > 1e-30 {
            (target - cum[idx - 1]) / denom
        } else {
            0.0
        };
        out.push(params[idx - 1] + alpha * (params[idx] - params[idx - 1]));
    }
    out.push(b);
    Ok(out)
}

/// Adapter from `&dyn Curve` to the `occt-core` sampler traits.
struct Adapter<'a>(&'a dyn Curve);

impl CurveSample for Adapter<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.0.d0(u)
    }
}

impl CurveDeriv for Adapter<'_> {
    fn tangent(&self, u: f64) -> GpVec {
        self.0.d1(u).1
    }
}

/// Parameters such that every chord deviates from the curve by at most `tol`
/// (port of `GCPnts_UniformDeflection`, parameter accessor).
pub fn uniform_deflection(c: &dyn Curve, tol: f64) -> Vec<f64> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return vec![a, b];
    }
    UniformDeflection::from_curve_with_deflection(&Adapter(c), a, b, tol).params
}

/// Parameters refined by both chord deviation `tol` and tangent-angle change
/// `angle_tol` (port of `GCPnts_TangentialDeflection`, parameter accessor).
pub fn tangential_deflection(c: &dyn Curve, tol: f64, angle_tol: f64) -> Vec<f64> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite()) {
        return vec![a, b];
    }
    TangentialDeflection::from_curve_with_deriv(&Adapter(c), a, b, tol, angle_tol).params
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GeomCircle, GeomLine, GeomTrimmedCurve};
    use occt_core::elib::clib;
    use occt_core::gp::{GpAx2, GpCirc, GpDir, GpPnt};
    use std::sync::Arc;

    const PI: f64 = std::f64::consts::PI;

    fn unit_circle() -> GeomCircle {
        GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))
    }

    fn unit_circle_arc() -> GeomTrimmedCurve {
        GeomTrimmedCurve::new(Arc::new(unit_circle()), 0.0, PI / 2.0)
    }

    #[test]
    fn length_unit_circle() {
        let c = unit_circle();
        let l = curve_length(&c);
        assert!((l - 2.0 * PI).abs() < 1e-7, "length {l}");
    }

    #[test]
    fn length_quarter_arc() {
        let c = unit_circle_arc();
        let l = curve_length(&c);
        assert!((l - PI / 2.0).abs() < 1e-7, "length {l}");
    }

    #[test]
    fn abscissa_point_half_circle() {
        let c = unit_circle();
        // Arc length π from u=0 lands at u=π.
        let u = abscissa_point(&c, PI, 0.0).unwrap();
        assert!((u - PI).abs() < 1e-6, "u {u}");
        let p = c.d0(u);
        assert!(p.distance(&GpPnt::new(-1.0, 0.0, 0.0)) < 1e-6, "point {p:?}");
    }

    #[test]
    fn abscissa_point_total_length_returns_last_parameter() {
        let c = unit_circle();
        let u = abscissa_point(&c, 2.0 * PI, 0.0).unwrap();
        assert!((u - 2.0 * PI).abs() < 1e-6, "u {u}");
        // Beyond total length errors.
        assert!(abscissa_point(&c, 3.0 * PI, 0.0).is_err());
    }

    #[test]
    fn abscissa_point_negative_goes_backward() {
        let c = unit_circle();
        // Backward by π/2 from u=π lands at u=π/2 (the point (0,1,0)).
        let u = abscissa_point(&c, -PI / 2.0, PI).unwrap();
        assert!((u - PI / 2.0).abs() < 1e-6, "u {u}");
        let p = c.d0(u);
        assert!(p.distance(&GpPnt::new(0.0, 1.0, 0.0)) < 1e-6, "point {p:?} at u {u}");
    }

    #[test]
    fn uniform_abscissa_unit_circle_quarter_steps() {
        let c = unit_circle();
        let params = uniform_abscissa(&c, 4).unwrap();
        assert_eq!(params.len(), 5);
        let expected = [0.0, PI / 2.0, PI, 3.0 * PI / 2.0, 2.0 * PI];
        for (g, e) in params.iter().zip(expected.iter()) {
            assert!((g - e).abs() < 1e-6, "got {g} expected {e}: {params:?}");
        }
        // Points are equally spaced on the circle.
        let pts: Vec<GpPnt> = params.iter().map(|&u| c.d0(u)).collect();
        for w in pts.windows(2) {
            let chord = w[0].distance(&w[1]);
            assert!((chord - (2.0f64).sqrt()).abs() < 1e-5, "chord {chord}");
        }
    }

    #[test]
    fn uniform_abscissa_line_is_uniform() {
        // A 10-unit trimmed line has parameter range [0, 1]; equal arc-length
        // spacing of 10/4 = 2.5 lands at params 0, 0.25, 0.5, 0.75, 1.
        let line = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let c = GeomTrimmedCurve::new(Arc::new(line), 0.0, 10.0);
        let params = uniform_abscissa(&c, 4).unwrap();
        for (i, u) in params.iter().enumerate() {
            assert!((u - 0.25 * i as f64).abs() < 1e-9, "params {params:?}");
        }
    }

    #[test]
    fn uniform_abscissa_degenerate_zero_length() {
        let line = GeomLine::from_pnt_dir(GpPnt::new(1., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let c = GeomTrimmedCurve::new(Arc::new(line), 0.0, 0.0);
        let params = uniform_abscissa(&c, 3).unwrap();
        assert_eq!(params, vec![0.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn quasi_uniform_abscissa_line() {
        // 10-unit trimmed line over [0, 1]: 5 points at equal chord intervals
        // 2.5 → params 0, 0.25, 0.5, 0.75, 1.
        let line = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let c = GeomTrimmedCurve::new(Arc::new(line), 0.0, 10.0);
        let params = quasi_uniform_abscissa(&c, 5).unwrap();
        assert_eq!(params.len(), 5);
        assert!((params[0] - 0.0).abs() < 1e-9);
        assert!((params[4] - 1.0).abs() < 1e-9);
        assert!((params[1] - 0.25).abs() < 1e-6, "params {params:?}");
        assert!((params[2] - 0.5).abs() < 1e-6, "params {params:?}");
        assert!((params[3] - 0.75).abs() < 1e-6, "params {params:?}");
    }

    #[test]
    fn uniform_deflection_circle_within_tolerance() {
        let c = unit_circle();
        let tol = 0.01;
        let params = uniform_deflection(&c, tol);
        assert!(params.len() >= 9, "only {} params", params.len());
        // Chord error over every span is within tol.
        for w in params.windows(2) {
            let (ua, ub) = (w[0], w[1]);
            let (pa, pb) = (c.d0(ua), c.d0(ub));
            for k in 1..9 {
                let p = c.d0(ua + (ub - ua) * k as f64 / 9.0);
                let d = point_seg_dist(&p, &pa, &pb);
                assert!(d <= tol + 1e-9, "deviation {d} > tol {tol} in [{ua}, {ub}]");
            }
        }
        // First and last parameters are the range bounds.
        assert!((params[0] - 0.0).abs() < 1e-12);
        assert!((params[params.len() - 1] - 2.0 * PI).abs() < 1e-12);
    }

    #[test]
    fn tangential_deflection_refines() {
        let c = unit_circle();
        let params = tangential_deflection(&c, 0.5, 0.1);
        assert!(params.len() > 8, "only {} params", params.len());
        assert!((params[0] - 0.0).abs() < 1e-12);
        assert!((params[params.len() - 1] - 2.0 * PI).abs() < 1e-12);
    }

    fn point_seg_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
        let ab = GpVec::from_pnts(a, b);
        let l2 = ab.square_magnitude();
        if l2 <= f64::EPSILON {
            return p.distance(a);
        }
        let t = ((p.x() - a.x()) * ab.x() + (p.y() - a.y()) * ab.y() + (p.z() - a.z()) * ab.z()) / l2;
        let t = t.clamp(0.0, 1.0);
        let q = GpPnt::new(a.x() + t * ab.x(), a.y() + t * ab.y(), a.z() + t * ab.z());
        p.distance(&q)
    }

    #[test]
    fn sampled_circle_points_consistent() {
        let c = unit_circle();
        let params = uniform_abscissa(&c, 8).unwrap();
        for &u in &params {
            let p = c.d0(u);
            let q = clib::circle_value(c.circ(), u);
            assert!(p.distance(&q) < 1e-12);
        }
    }
}
