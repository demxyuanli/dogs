//! Polyline approximation and parametrization utilities.
//! Used for tessellation and curve discretization.
use crate::gp::{GpPnt, GpVec};

/// Adaptive polyline sampling of a parametric curve.
/// curve: evaluates point at parameter u.
/// a, b: parameter interval. tol: max chord deviation from curve.
/// Returns (points, parameters).
pub fn adaptive_sample<F: Fn(f64) -> GpPnt>(curve: &F, a: f64, b: f64, tol: f64, max_depth: usize) -> (Vec<GpPnt>, Vec<f64>) {
    let mut pts = Vec::new();
    let mut params = Vec::new();
    let mut stack = vec![(a, b)];
    let mut depth = 0usize;

    while let Some((u0, u1)) = stack.pop() {
        if depth > max_depth { break; }
        let p0 = curve(u0);
        let p1 = curve(u1);
        let um = 0.5 * (u0 + u1);
        let pm = curve(um);
        let chord = p0.coord.subtracted(&p1.coord).modulus();
        if chord < 1e-30 {
            pts.push(p0); params.push(u0);
            continue;
        }
        // Distance from midpoint to chord
        let d = point_segment_distance(&pm, &p0, &p1);
        if d < tol && (u1 - u0) < (b - a) * 0.5 {
            if pts.is_empty() { pts.push(p0); params.push(u0); }
            pts.push(pm); params.push(um);
            pts.push(p1); params.push(u1);
        } else {
            stack.push((um, u1));
            stack.push((u0, um));
            depth += 1;
        }
    }
    (pts, params)
}

/// Uniform polyline sampling (fixed count).
pub fn uniform_sample<F: Fn(f64) -> GpPnt>(curve: &F, a: f64, b: f64, n: usize) -> (Vec<GpPnt>, Vec<f64>) {
    let mut pts = Vec::with_capacity(n);
    let mut params = Vec::with_capacity(n);
    if n == 0 { return (pts, params); }
    if n == 1 { let p = curve(a); pts.push(p); params.push(a); return (pts, params); }
    for i in 0..n {
        let u = a + (b - a) * i as f64 / (n - 1) as f64;
        pts.push(curve(u));
        params.push(u);
    }
    (pts, params)
}

/// Chord-length parameterization of sampled points → cumulative parameter values.
pub fn chord_length_params(pts: &[GpPnt]) -> Vec<f64> {
    let n = pts.len();
    if n < 2 { return vec![0.0]; }
    let mut params = vec![0.0f64; n];
    for i in 1..n {
        params[i] = params[i-1] + pts[i].coord.subtracted(&pts[i-1].coord).modulus();
    }
    let total = params[n-1];
    if total > 1e-30 { for p in &mut params { *p /= total; } }
    params
}

/// Distance from point to segment.
pub fn point_segment_distance(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let ap = p.coord.subtracted(&a.coord);
    let len2 = ab.square_modulus();
    if len2 < 1e-30 { return ap.modulus(); }
    let t = (ap.dot(&ab) / len2).clamp(0.0, 1.0);
    let proj = ab.multiplied(t);
    ap.subtracted(&proj).modulus()
}

/// Approximate total arc length of a curve by sampling.
pub fn arc_length<F: Fn(f64) -> GpPnt>(curve: &F, a: f64, b: f64, n: usize) -> f64 {
    let (pts, _) = uniform_sample(curve, a, b, n);
    let mut total = 0.0;
    for w in pts.windows(2) {
        total += w[1].coord.subtracted(&w[0].coord).modulus();
    }
    total
}

/// Tangent direction at a curve point via finite differences.
pub fn tangent_at<F: Fn(f64) -> GpPnt>(curve: &F, u: f64, h: f64) -> GpVec {
    let p1 = curve(u + h);
    let p2 = curve(u - h);
    let v = p1.coord.subtracted(&p2.coord).divided(2.0 * h);
    let m = v.modulus();
    if m > 1e-30 { GpVec::from_xyz(&v.divided(m)) } else { GpVec::zero() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_line() {
        let f = |u: f64| GpPnt::new(u, 0.0, 0.0);
        let (pts, params) = uniform_sample(&f, 0.0, 1.0, 3);
        assert_eq!(pts.len(), 3);
        assert!((pts[0].x() - 0.0).abs() < 1e-14);
        assert!((pts[2].x() - 1.0).abs() < 1e-14);
        assert!((params[1] - 0.5).abs() < 1e-14);
    }

    #[test]
    fn arc_length_line() {
        let f = |u: f64| GpPnt::new(u, u, 0.0);
        let len = arc_length(&f, 0.0, 1.0, 100);
        assert!((len - 2.0f64.sqrt()).abs() < 1e-6);
    }

    #[test]
    fn point_seg_dist() {
        let a = GpPnt::new(0., 0., 0.);
        let b = GpPnt::new(1., 0., 0.);
        let d = point_segment_distance(&GpPnt::new(0.5, 1., 0.), &a, &b);
        assert!((d - 1.0).abs() < 1e-14);
    }
}
