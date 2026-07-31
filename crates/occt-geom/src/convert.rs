//! Curve/surface conversion utilities for the `occt-geom` crate.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpLin, GpPnt};

use crate::curve::Curve;
use crate::line::GeomLine;
use crate::surface::Surface;

/// Wrap a concrete [`Curve`] value in a trait-object [`Arc`].
pub fn curve_to_arc<C: Curve + 'static>(c: C) -> Arc<dyn Curve> {
    Arc::new(c)
}

/// Wrap a concrete [`Surface`] value in a trait-object [`Arc`].
pub fn surface_to_arc<S: Surface + 'static>(s: S) -> Arc<dyn Surface> {
    Arc::new(s)
}

/// Recover the underlying [`GpLin`] when `curve` is backed by a [`GeomLine`].
///
/// Stub: real type-erasure would require `fn as_any(&self) -> &dyn std::any::Any`
/// on the [`Curve`] trait before the `Arc<dyn Curve>` could be downcast. Once
/// that exists this can return `Some(...)` for the concrete geometry types and
/// `None` for everything else.
pub fn try_line(_curve: &Arc<dyn Curve>) -> Option<GpLin> {
    None
}

/// Sample a curve into an open polyline of `n` segments (`n + 1` points).
///
/// Points are evaluated with [`Curve::d0`] at uniformly spaced parameters in
/// `[a, b]`.
pub fn sample_curve_to_polyline(
    curve: &Arc<dyn Curve>,
    a: f64,
    b: f64,
    n: usize,
) -> Vec<GpPnt> {
    if n == 0 {
        return vec![curve.d0(a)];
    }
    (0..=n)
        .map(|i| curve.d0(a + (b - a) * i as f64 / n as f64))
        .collect()
}

/// Axis-aligned bounding box of a sampled curve.
///
/// The box is the min/max of the polyline produced by
/// [`sample_curve_to_polyline`], so it only contains sampled points.
pub fn curve_bbox(curve: &Arc<dyn Curve>, a: f64, b: f64, n: usize) -> occt_core::bnd::BndBox {
    let mut bb = occt_core::bnd::BndBox::new();
    for p in sample_curve_to_polyline(curve, a, b, n) {
        bb.add(p);
    }
    bb
}

/// Approximate curve length by summing the chord lengths of a sampled polyline.
pub fn approximate_curve_length(curve: &Arc<dyn Curve>, a: f64, b: f64, n: usize) -> f64 {
    sample_curve_to_polyline(curve, a, b, n)
        .windows(2)
        .map(|w| (w[1] - w[0]).magnitude())
        .sum()
}

/// Whether the curve reports itself as periodic.
pub fn is_curve_periodic(curve: &Arc<dyn Curve>) -> bool {
    curve.is_periodic()
}

/// Continuity of the curve (`0` = C0, `1` = C1, ...).
pub fn curve_continuity(curve: &Arc<dyn Curve>) -> u8 {
    curve.continuity()
}

/// Build the unique line through `p1` and `p2`.
///
/// The direction is `p2 - p1` normalized; when the two points coincide the
/// default direction is used.
pub fn make_line_from_points(p1: &GpPnt, p2: &GpPnt) -> GeomLine {
    let v = *p2 - *p1;
    let dir = GpDir::from_xyz(&v.xyz()).unwrap_or_default();
    GeomLine::new(GpLin::from_pnt_dir(*p1, dir))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_x_line() -> Arc<dyn Curve> {
        Arc::new(make_line_from_points(
            &GpPnt::new(0.0, 0.0, 0.0),
            &GpPnt::new(1.0, 0.0, 0.0),
        ))
    }

    #[test]
    fn sampled_line_points_are_collinear() {
        let curve = unit_x_line();
        let pts = sample_curve_to_polyline(&curve, 0.0, 1.0, 8);
        assert_eq!(pts.len(), 9);
        for p in &pts {
            // All sampled points lie on the X axis.
            assert!(p.y().abs() < 1e-12);
            assert!(p.z().abs() < 1e-12);
        }
    }

    #[test]
    fn unit_line_length_is_approx_one() {
        let curve = unit_x_line();
        let len = approximate_curve_length(&curve, 0.0, 1.0, 16);
        assert!((len - 1.0).abs() < 1e-9, "len = {len}");
    }
}
