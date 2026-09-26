//! Curve conversion utilities for the `occt-geom2d` crate.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec2d};

use crate::curve::Curve2d;

// The former polyline substitutes `sample_curve2d` / `curve2d_bbox` /
// `curve2d_length(curve, a, b, n)` were removed here: they had no consumer
// left once the faithful ports landed (2-D length is
// `crate::curve_ops::curve2d_length` = `GCPnts_AbscissaPoint::Length`; the
// bounding box used in production is `crate::geom2d_api::curve2d_bbox`).
// Removing them closes the "curve2d_length = Simpson/chord" leftovers of
// board card R2-18.

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
