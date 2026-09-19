//! Point generators for curves (`GCPnts` package).
//!
//! **Ported (faithful)**: the `GCPnts_TangentialDeflection` engine lives in
//! `gcpnts_perform.rs` / `gcpnts_estim.rs` — `PerformLinear` (`cxx:321-348`),
//! `PerformCurve` (`cxx:522-916`) and `EstimDefl` (`cxx:960-1010`), re-exported
//! below as `perform_linear` / `perform_tangential_curve`.
//!
//! **UNPORTED** (do not fake these — no OCCT branch exists here yet):
//! - `GCPnts_UniformDeflection` (`GCPnts_UniformDeflection.cxx`): Linear /
//!   Circular / Curved / Composite dispatch plus the `Controle` end-point
//!   correction. Callers that need a chord-deflection polyline must currently
//!   use the tangential engine (what `BRepMesh_CurveTessellator` does in OCCT:
//!   `BRepMesh_CurveTessellator.cxx` -> `GCPnts_TangentialDeflection`).
//! - `GCPnts_QuasiUniformDeflection` (`GCPnts_QuasiUniformDeflection.cxx`).
//! - `GCPnts_UniformAbscissa` (`GCPnts_UniformAbscissa.cxx`).
//!
//! The remaining helpers in this file (`UniformPoints`, `sample_line`,
//! `sample_circle`, `polyline_length`, `reparametrize_by_chord_length`,
//! `evaluate_at`, `total_chord_error`) are **port-internal utilities, not OCCT
//! translations**; they carry no OCCT line references and must not be used as
//! evidence of alignment.

use crate::gp::GpPnt;

#[path = "gcpnts_perform.rs"]
mod perform;
pub use perform::{perform_linear, perform_tangential_curve, CurveSecondDeriv};

/// Minimal curve abstraction for parameter-uniform sampling.
pub trait CurveSample {
    /// Point of the curve at parameter `u`.
    fn point(&self, u: f64) -> GpPnt;
}

/// Uniformly spaced points over a parameter interval.
///
/// Port-internal utility (not an OCCT class; `GCPnts_UniformAbscissa` is
/// arc-length based and is UNPORTED, see the module docs).
pub struct UniformPoints {
    pub params: Vec<f64>,
    pub points: Vec<GpPnt>,
}

impl UniformPoints {
    /// `n` points spread uniformly over `[a, b]` (n-1 parameter steps).
    pub fn from_curve<C: CurveSample>(c: &C, a: f64, b: f64, n: usize) -> Self {
        let mut params = Vec::with_capacity(n);
        let mut points = Vec::with_capacity(n);
        if n == 0 {
            return Self { params, points };
        }
        if n == 1 {
            params.push(a);
            points.push(c.point(a));
            return Self { params, points };
        }
        let step = (b - a) / (n - 1) as f64;
        for i in 0..n {
            let u = a + step * i as f64;
            params.push(u);
            points.push(c.point(u));
        }
        Self { params, points }
    }

    /// Wrap already-sampled parameter/point pairs.
    pub fn from_sampled(params: Vec<f64>, points: Vec<GpPnt>) -> Self {
        Self { params, points }
    }

    /// Number of sampled points.
    pub fn nb_points(&self) -> usize {
        self.points.len()
    }

    /// `i`-th sampled point (0-based).
    pub fn point(&self, i: usize) -> GpPnt {
        self.points[i]
    }

    /// Parameter of the `i`-th sample (0-based).
    pub fn parameter(&self, i: usize) -> f64 {
        self.params[i]
    }
}

/// Total length of the polyline through `pts`.
pub fn polyline_length(pts: &[GpPnt]) -> f64 {
    pts.windows(2).map(|w| w[0].distance(&w[1])).sum()
}

/// Normalized cumulative chord-length parameters in `[0, 1]` (last is `1`).
pub fn reparametrize_by_chord_length(points: &[GpPnt]) -> Vec<f64> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    let mut cum = Vec::with_capacity(n);
    cum.push(0.0);
    for i in 1..n {
        cum.push(cum[i - 1] + points[i].distance(&points[i - 1]));
    }
    let total = *cum.last().unwrap();
    if total <= f64::EPSILON {
        return vec![0.0; n];
    }
    for v in cum.iter_mut() {
        *v /= total;
    }
    cum
}

/// Linearly interpolate the sampled polyline at parameter `u`.
pub fn evaluate_at(points: &[GpPnt], params: &[f64], u: f64) -> GpPnt {
    let n = points.len();
    if n == 0 {
        return GpPnt::new(0.0, 0.0, 0.0);
    }
    if n == 1 || u <= params[0] {
        return points[0];
    }
    if u >= params[n - 1] {
        return points[n - 1];
    }
    // Find the bracket i with params[i] <= u < params[i+1].
    let mut i = 0;
    while i + 1 < n && params[i + 1] < u {
        i += 1;
    }
    let (ua, ub) = (params[i], params[i + 1]);
    let t = if ub > ua { (u - ua) / (ub - ua) } else { 0.0 };
    lerp(&points[i], &points[i + 1], t)
}

/// `n` points uniformly spaced along the line segment `a..b` (inclusive).
pub fn sample_line(a: &GpPnt, b: &GpPnt, n: usize) -> Vec<GpPnt> {
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = if n == 1 { 0.0 } else { i as f64 / (n - 1) as f64 };
        pts.push(lerp(a, b, t));
    }
    pts
}

/// `n` points evenly spaced around a circle of `radius` centered at `center`.
///
/// The circle lies in the XY plane by default; pass `z_normal = true` to lay
/// it in the XZ plane instead.
pub fn sample_circle(center: &GpPnt, radius: f64, n: usize, z_normal: bool) -> Vec<GpPnt> {
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let th = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let (cx, cy) = (radius * th.cos(), radius * th.sin());
        let (dx, dy, dz) = if z_normal { (cx, 0.0, cy) } else { (cx, cy, 0.0) };
        pts.push(GpPnt::new(center.x() + dx, center.y() + dy, center.z() + dz));
    }
    pts
}

/// Maximum deviation of the sampled polyline from the true curve.
///
/// For each span the true curve is sampled at 8 interior parameters and the
/// furthest of those points from the chord is taken; the result is the
/// maximum over all spans. Port-internal measurement helper (not OCCT).
pub fn total_chord_error(points: &[GpPnt], orig: &dyn Fn(f64) -> GpPnt, params: &[f64]) -> f64 {
    let mut max_dev = 0.0f64;
    let n = params.len().min(points.len());
    for i in 0..n.saturating_sub(1) {
        let pa = &points[i];
        let pb = &points[i + 1];
        let (ua, ub) = (params[i], params[i + 1]);
        for k in 1..9 {
            let u = ua + (ub - ua) * k as f64 / 9.0;
            let d = point_segment_dist(&orig(u), pa, pb);
            if d > max_dev {
                max_dev = d;
            }
        }
    }
    max_dev
}

/// Linear interpolation between two points.
fn lerp(a: &GpPnt, b: &GpPnt, t: f64) -> GpPnt {
    GpPnt::new(
        a.x() + t * (b.x() - a.x()),
        a.y() + t * (b.y() - a.y()),
        a.z() + t * (b.z() - a.z()),
    )
}

/// Distance from point `p` to the line segment `a..b`.
pub(crate) fn point_segment_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let abx = b.x() - a.x();
    let aby = b.y() - a.y();
    let abz = b.z() - a.z();
    let len2 = abx * abx + aby * aby + abz * abz;
    if len2 <= f64::EPSILON {
        return p.distance(a);
    }
    let t = ((p.x() - a.x()) * abx + (p.y() - a.y()) * aby + (p.z() - a.z()) * abz) / len2;
    let t = t.clamp(0.0, 1.0);
    let (cx, cy, cz) = (a.x() + t * abx, a.y() + t * aby, a.z() + t * abz);
    let (dx, dy, dz) = (p.x() - cx, p.y() - cy, p.z() - cz);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    struct LinearX;
    impl CurveSample for LinearX {
        fn point(&self, u: f64) -> GpPnt {
            GpPnt::new(u, 0.0, 0.0)
        }
    }

    #[test]
    fn uniform_points_linear() {
        let up = UniformPoints::from_curve(&LinearX, 0.0, 10.0, 11);
        assert_eq!(up.nb_points(), 11);
        let p5 = up.point(5);
        assert!((p5.x() - 5.0).abs() < 1e-12, "p5.x={}", p5.x());
        assert!(p5.y().abs() < 1e-12 && p5.z().abs() < 1e-12);
        assert!((up.parameter(10) - 10.0).abs() < 1e-12);
    }

    #[test]
    fn sample_line_and_circle() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(10.0, 0.0, 0.0);
        let line = sample_line(&a, &b, 5);
        assert_eq!(line.len(), 5);
        assert!(line[0].distance(&a) < 1e-12);
        assert!(line[4].distance(&b) < 1e-12);
        assert!((line[2].x() - 5.0).abs() < 1e-12);

        let center = GpPnt::new(1.0, 2.0, 3.0);
        let circ = sample_circle(&center, 2.0, 8, false);
        assert_eq!(circ.len(), 8);
        for p in &circ {
            assert!((p.distance(&center) - 2.0).abs() < 1e-12);
            assert!((p.z() - 3.0).abs() < 1e-12);
        }
        let circ_xz = sample_circle(&center, 2.0, 8, true);
        for p in &circ_xz {
            assert!((p.distance(&center) - 2.0).abs() < 1e-12);
            assert!((p.y() - 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn chord_parametrization_and_length() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1.0, 1.0, 1.0),
        ];
        assert!((polyline_length(&pts) - 3.0).abs() < 1e-12);
        let p = reparametrize_by_chord_length(&pts);
        assert_eq!(p.len(), 4);
        assert!(p[0].abs() < 1e-12);
        assert!((p[1] - 1.0 / 3.0).abs() < 1e-12);
        assert!((p[2] - 2.0 / 3.0).abs() < 1e-12);
        assert!((p[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn evaluate_at_interpolates() {
        let pts = vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(10.0, 0.0, 0.0)];
        let params = vec![0.0, 10.0];
        let mid = evaluate_at(&pts, &params, 5.0);
        assert!((mid.x() - 5.0).abs() < 1e-12);
        assert!((evaluate_at(&pts, &params, -1.0).x() - 0.0).abs() < 1e-12);
        assert!((evaluate_at(&pts, &params, 99.0).x() - 10.0).abs() < 1e-12);
    }
}
