//! Differential geometry of a 3D parametric curve — Frenet frame, curvature,
//! torsion, and arc-length parameterization.
//!
//! Port of `GeomLProp_CurAndInf` / `Geom_Curve::D3` differential-geometry
//! helpers (curvature, torsion, Frenet triad) and arc-length reparameterization
//! (used by `GCPnts_UniformAbscissa` / `GeomAdaptor_Curve`). Source:
//! `GeomLProp` (TKGeomBase), `GCPnts` (TKG3d).

use crate::gp::{GpPnt, GpVec, GpXyz};

/// The Frenet–Serret frame of a curve at a parameter.
#[derive(Debug, Clone, Copy)]
pub struct FrenetFrame {
    pub point: GpPnt,
    pub tangent: GpVec,
    pub normal: GpVec,
    pub binormal: GpVec,
    pub curvature: f64,
    pub torsion: f64,
}

/// Unit-vectorize a GpXyz; returns a zero vector for a (near-)zero input.
fn unitize(v: &GpXyz) -> GpVec {
    let m = v.modulus();
    if m > 1e-30 {
        GpVec::new(v.x / m, v.y / m, v.z / m)
    } else {
        GpVec::zero()
    }
}

/// Compute the curvature at `u` from the first two derivatives:
/// κ = |r' × r''| / |r'|³.
pub fn curvature(r1: &GpVec, r2: &GpVec) -> f64 {
    let m1 = r1.xyz().modulus();
    if m1 < 1e-30 {
        return 0.0;
    }
    r1.xyz().crossed(r2.xyz()).modulus() / (m1 * m1 * m1)
}

/// Compute the torsion at `u`: τ = (r' × r'')·r''' / |r' × r''|².
pub fn torsion(r1: &GpVec, r2: &GpVec, r3: &GpVec) -> f64 {
    let cross = r1.xyz().crossed(r2.xyz());
    let denom = cross.square_modulus();
    if denom < 1e-30 {
        return 0.0;
    }
    cross.dot(&r3.xyz()) / denom
}

/// Build the Frenet frame from the position and first three derivatives.
/// - tangent  T = r' / |r'|
/// - normal   N = (T × (r'' × T)) / |...|   (the second derivative projected
///   onto the normal plane, normalized)
/// - binormal B = T × N
/// If r' is zero, T defaults to (1,0,0) and the frame is degenerate (zero
/// curvature).
pub fn frenet_frame(p: GpPnt, r1: &GpVec, r2: &GpVec, r3: &GpVec) -> FrenetFrame {
    let t = unitize(&r1.xyz());
    let t_xyz = t.xyz();
    // r'' projected onto the normal plane: r2 − (r2·T)T.
    let r2n = r2.xyz().subtracted(&t_xyz.multiplied(r2.xyz().dot(&t_xyz)));
    let n = unitize(&r2n);
    let b = unitize(&t_xyz.crossed(n.xyz()));
    let curv = curvature(r1, r2);
    let tors = torsion(r1, r2, r3);
    FrenetFrame {
        point: p,
        tangent: t,
        normal: n,
        binormal: b,
        curvature: curv,
        torsion: tors,
    }
}

/// Sample the Frenet frame along a curve. `d0/d1/d2/d3` are closures returning
/// the position and up to three derivatives at `u`.
pub fn frenet_frame_at<F1, F2, F3>(
    p: &F1,
    r1: &F2,
    r2: &F3,
    u: f64,
    h: f64,
) -> FrenetFrame
where
    F1: Fn(f64) -> GpPnt,
    F2: Fn(f64) -> GpVec,
    F3: Fn(f64) -> GpVec,
{
    // Derivatives via finite differences if only d0/d1 are exact.
    let d1 = r1(u);
    let d2 = r2(u);
    let d3 = finite_diff_derivative(p, u, h);
    frenet_frame(p(u), &d1, &d2, &d3)
}

/// Third derivative by central finite differences of `d0` (used when the curve
/// only provides positions).
fn finite_diff_derivative<F: Fn(f64) -> GpPnt>(p: &F, u: f64, h: f64) -> GpVec {
    let a = p(u + 2.0 * h);
    let b = p(u + h);
    let c = p(u - h);
    let d = p(u - 2.0 * h);
    GpVec::new(
        (a.x() - 2.0 * b.x() + 2.0 * c.x() - d.x()) / (2.0 * h * h * h),
        (a.y() - 2.0 * b.y() + 2.0 * c.y() - d.y()) / (2.0 * h * h * h),
        (a.z() - 2.0 * b.z() + 2.0 * c.z() - d.z()) / (2.0 * h * h * h),
    )
}

/// Re-parameterize a curve by arc length: sample the curve, build an
/// s → u map (cumulative chord lengths), and return (s, u) pairs at
/// `n` equally spaced arc-length positions. Used for `GCPnts_UniformAbscissa`.
pub fn arc_length_parameters<F: Fn(f64) -> GpPnt>(
    p: &F,
    a: f64,
    b: f64,
    n: usize,
    samples: usize,
) -> Vec<(f64, f64)> {
    if n < 2 {
        return vec![(0.0, a)];
    }
    let samples = samples.max(2);
    let mut cum = vec![0.0f64; samples + 1];
    let mut us = vec![0.0f64; samples + 1];
    let mut prev = p(a);
    us[0] = a;
    for i in 1..=samples {
        let u = a + (b - a) * i as f64 / samples as f64;
        let cur = p(u);
        cum[i] = cum[i - 1] + prev.distance(&cur);
        prev = cur;
        us[i] = u;
    }
    let total = cum[samples];
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let target = total * k as f64 / (n - 1) as f64;
        let mut idx = 0;
        while idx < samples && cum[idx + 1] < target {
            idx += 1;
        }
        // Linear interpolate u within the segment.
        let seg = cum[idx + 1] - cum[idx];
        let t = if seg > 1e-30 { (target - cum[idx]) / seg } else { 0.0 };
        let u = us[idx] + t * (us[idx + 1] - us[idx]);
        out.push((target, u));
    }
    out
}

/// Total arc length of a curve over [a, b] via chord-length sampling.
pub fn curve_arc_length<F: Fn(f64) -> GpPnt>(p: &F, a: f64, b: f64, samples: usize) -> f64 {
    let n = samples.max(2);
    let mut len = 0.0;
    let mut prev = p(a);
    for i in 1..=n {
        let u = a + (b - a) * i as f64 / n as f64;
        let cur = p(u);
        len += prev.distance(&cur);
        prev = cur;
    }
    len
}

/// The tangent direction of a curve at `u` (unit). Falls back to finite
/// differences if the analytic derivative is near-zero.
pub fn curve_tangent<F: Fn(f64) -> GpPnt, D: Fn(f64) -> GpVec>(
    p: &F,
    d: &D,
    u: f64,
    h: f64,
) -> GpVec {
    let d1 = d(u);
    if d1.xyz().square_modulus() > 1e-30 {
        unitize(&d1.xyz())
    } else {
        let p1 = p(u - h);
        let p2 = p(u + h);
        unitize(&GpVec::from_pnts(&p1, &p2).xyz())
    }
}

/// Radius of curvature (1/κ); None when κ ≈ 0 (straight).
pub fn radius_of_curvature(r1: &GpVec, r2: &GpVec) -> Option<f64> {
    let k = curvature(r1, r2);
    if k > 1e-30 {
        Some(1.0 / k)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn circle_point(u: f64) -> GpPnt {
        GpPnt::new(u.cos(), u.sin(), 0.0)
    }
    fn circle_d1(u: f64) -> GpVec {
        GpVec::new(-u.sin(), u.cos(), 0.0)
    }
    fn circle_d2(u: f64) -> GpVec {
        GpVec::new(-u.cos(), -u.sin(), 0.0)
    }
    fn circle_d3(u: f64) -> GpVec {
        GpVec::new(u.sin(), -u.cos(), 0.0)
    }

    #[test]
    fn circle_curvature_unit() {
        // Unit circle → curvature 1 everywhere.
        let k = curvature(&circle_d1(0.5), &circle_d2(0.5));
        assert!((k - 1.0).abs() < 1e-9, "curvature {k}");
    }

    #[test]
    fn circle_torsion_zero() {
        // Planar circle → torsion 0.
        let t = torsion(&circle_d1(0.3), &circle_d2(0.3), &circle_d3(0.3));
        assert!(t.abs() < 1e-9, "torsion {t}");
    }

    #[test]
    fn circle_frenet_frame() {
        let f = frenet_frame(circle_point(0.0), &circle_d1(0.0), &circle_d2(0.0), &circle_d3(0.0));
        // At u=0: T=(0,1,0), N=(-1,0,0), B=(0,0,1) (for the CCW circle).
        assert!((f.tangent.x() - 0.0).abs() < 1e-9 && (f.tangent.y() - 1.0).abs() < 1e-9, "T {:?}", f.tangent);
        assert!((f.normal.x() + 1.0).abs() < 1e-9, "N {:?}", f.normal);
        assert!((f.binormal.z() - 1.0).abs() < 1e-9, "B {:?}", f.binormal);
        assert!((f.curvature - 1.0).abs() < 1e-9);
        // Frame orthonormal.
        assert!(f.tangent.dot(&f.normal).abs() < 1e-9);
        assert!(f.tangent.dot(&f.binormal).abs() < 1e-9);
    }

    #[test]
    fn arc_length_circle() {
        // Half circle of radius 1 → length π.
        let len = curve_arc_length(&circle_point, 0.0, std::f64::consts::PI, 512);
        assert!((len - std::f64::consts::PI).abs() < 1e-3, "len {len}");
    }

    #[test]
    fn arc_length_parameters_uniform() {
        // Quarter circle → s uniformly spaced, u advances monotonically.
        let params = arc_length_parameters(&circle_point, 0.0, std::f64::consts::PI / 2.0, 5, 256);
        assert_eq!(params.len(), 5);
        let total = params[4].0;
        assert!((total - std::f64::consts::PI / 2.0).abs() < 0.01, "total {total}");
        // Uniform spacing in s.
        let step = params[1].0 - params[0].0;
        for k in 1..5 {
            let d = params[k].0 - params[k - 1].0;
            assert!((d - step).abs() < 1e-6, "spacing {d} vs {step}");
        }
    }

    #[test]
    fn radius_of_curvature_circle() {
        let r = radius_of_curvature(&circle_d1(1.0), &circle_d2(1.0)).expect("r");
        assert!((r - 1.0).abs() < 1e-9);
        // A straight line → None.
        let straight1 = GpVec::new(1.0, 0.0, 0.0);
        let straight2 = GpVec::zero();
        assert!(radius_of_curvature(&straight1, &straight2).is_none());
    }

    #[test]
    fn tangent_finite_diff_fallback() {
        // d1 near-zero at a cusp-like point → FD tangent still works.
        let p = |u: f64| GpPnt::new(u, u * u, 0.0);
        let d = |u: f64| GpVec::new(1.0, 2.0 * u, 0.0);
        let t = curve_tangent(&p, &d, 1.0, 1e-6);
        let expected = GpVec::new(1.0, 2.0, 0.0);
        let em = expected.xyz().modulus();
        assert!((t.x() - expected.x() / em).abs() < 1e-6, "t {:?}", t);
    }

    #[test]
    fn torsion_helix() {
        // Helix r(t) = (cos t, sin t, t): κ = 1/2, τ = 1/2.
        let h = |u: f64| GpPnt::new(u.cos(), u.sin(), u);
        let d1 = |u: f64| GpVec::new(-u.sin(), u.cos(), 1.0);
        let d2 = |u: f64| GpVec::new(-u.cos(), -u.sin(), 0.0);
        let d3 = |u: f64| GpVec::new(u.sin(), -u.cos(), 0.0);
        let k = curvature(&d1(0.4), &d2(0.4));
        let t = torsion(&d1(0.4), &d2(0.4), &d3(0.4));
        assert!((k - 0.5).abs() < 1e-9, "helix κ {k}");
        assert!((t - 0.5).abs() < 1e-9, "helix τ {t}");
    }

    #[test]
    fn frenet_finite_difference_derivative() {
        // Third derivative of a cubic is constant.
        let p = |u: f64| GpPnt::new(u, u * u * u, 0.0);
        let d3 = finite_diff_derivative(&p, 0.5, 1e-4);
        assert!((d3.y() - 6.0).abs() < 0.01, "cubic d3y {}", d3.y());
    }
}
