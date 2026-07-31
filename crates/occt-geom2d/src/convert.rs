//! Curve conversion utilities for the `occt-geom2d` crate.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec2d};

use crate::curve::Curve2d;

/// Sample a 2D curve into an open polyline of `n` segments (`n + 1` points).
///
/// Points are evaluated with [`Curve2d::d0`] at uniformly spaced parameters in
/// `[a, b]`.
pub fn sample_curve2d(curve: &Arc<dyn Curve2d>, a: f64, b: f64, n: usize) -> Vec<GpPnt2d> {
    if n == 0 {
        return vec![curve.d0(a)];
    }
    (0..=n)
        .map(|i| curve.d0(a + (b - a) * i as f64 / n as f64))
        .collect()
}

/// Axis-aligned bounding box of a sampled curve, as `(xmin, xmax, ymin, ymax)`.
///
/// Computed from the polyline produced by [`sample_curve2d`].
pub fn curve2d_bbox(
    curve: &Arc<dyn Curve2d>,
    a: f64,
    b: f64,
    n: usize,
) -> (f64, f64, f64, f64) {
    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for p in sample_curve2d(curve, a, b, n) {
        xmin = xmin.min(p.x());
        xmax = xmax.max(p.x());
        ymin = ymin.min(p.y());
        ymax = ymax.max(p.y());
    }
    (xmin, xmax, ymin, ymax)
}

/// Approximate 2D curve length by summing chord lengths of a sampled polyline.
pub fn curve2d_length(curve: &Arc<dyn Curve2d>, a: f64, b: f64, n: usize) -> f64 {
    sample_curve2d(curve, a, b, n)
        .windows(2)
        .map(|w| (w[1] - w[0]).magnitude())
        .sum()
}

/// Unit tangent vector at parameter `u`, computed from the normalized first
/// derivative. Degenerate derivatives yield the +X axis.
pub fn tangent2d(curve: &Arc<dyn Curve2d>, u: f64) -> GpVec2d {
    let (_, v) = curve.d1(u);
    let m = v.magnitude();
    if m > 1e-15 {
        GpVec2d::new(v.x() / m, v.y() / m)
    } else {
        GpVec2d::new(1.0, 0.0)
    }
}

/// Point on the curve at parameter `u`.
pub fn point_at(curve: &Arc<dyn Curve2d>, u: f64) -> GpPnt2d {
    curve.d0(u)
}

/// Whether the curve is closed: distance between the endpoint samples ≤ `tol`.
pub fn is_closed2d(curve: &Arc<dyn Curve2d>, tol: f64) -> bool {
    let p0 = curve.d0(curve.first_param());
    let p1 = curve.d0(curve.last_param());
    (p1 - p0).magnitude() <= tol
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::line::Geom2dLine;
    use occt_core::gp::{GpDir2d, GpLin2d};

    fn unit_x_line() -> Arc<dyn Curve2d> {
        Arc::new(Geom2dLine::new(GpLin2d::from_pnt_dir(
            GpPnt2d::new(0.0, 0.0),
            GpDir2d::new(1.0, 0.0),
        )))
    }

    #[test]
    fn tangent_at_midpoint_has_unit_length() {
        let curve = unit_x_line();
        let t = tangent2d(&curve, 0.5);
        assert!((t.magnitude() - 1.0).abs() < 1e-9);
    }
}
