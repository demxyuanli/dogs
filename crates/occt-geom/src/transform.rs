//! Transform helpers for curves and surfaces.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpPnt, GpTrsf, GpVec};

use crate::curve::Curve;
use crate::surface::Surface;

/// Apply an arbitrary transformation to a curve.
pub fn transform_curve(curve: &Arc<dyn Curve>, t: &GpTrsf) -> Arc<dyn Curve> {
    curve.transformed(t)
}

/// Apply an arbitrary transformation to a surface.
pub fn transform_surface(surf: &Arc<dyn Surface>, t: &GpTrsf) -> Arc<dyn Surface> {
    surf.transformed(t)
}

/// Translate a curve by vector `v`.
pub fn translate_curve(curve: &Arc<dyn Curve>, v: &GpVec) -> Arc<dyn Curve> {
    let mut t = GpTrsf::default();
    t.set_translation_vec(v);
    transform_curve(curve, &t)
}

/// Translate a surface by vector `v`.
pub fn translate_surface(surf: &Arc<dyn Surface>, v: &GpVec) -> Arc<dyn Surface> {
    let mut t = GpTrsf::default();
    t.set_translation_vec(v);
    transform_surface(surf, &t)
}

/// Mirror a curve about a point.
pub fn mirror_curve(curve: &Arc<dyn Curve>, p: &GpPnt) -> Arc<dyn Curve> {
    let mut t = GpTrsf::default();
    t.set_mirror_pnt(p);
    transform_curve(curve, &t)
}

/// Rotate a curve about an axis by `angle` radians.
pub fn rotate_curve(curve: &Arc<dyn Curve>, axis: &GpAx1, angle: f64) -> Arc<dyn Curve> {
    let mut t = GpTrsf::default();
    t.set_rotation_ax1(axis, angle);
    transform_curve(curve, &t)
}

/// Scale a curve about point `p` by factor `s`.
pub fn scale_curve(curve: &Arc<dyn Curve>, p: &GpPnt, s: f64) -> Arc<dyn Curve> {
    let mut t = GpTrsf::default();
    t.set_scale(p, s);
    transform_curve(curve, &t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::convert::make_line_from_points;

    #[test]
    fn translate_moves_start_point() {
        let line = make_line_from_points(
            &GpPnt::new(0.0, 0.0, 0.0),
            &GpPnt::new(1.0, 0.0, 0.0),
        );
        let curve: Arc<dyn Curve> = Arc::new(line);
        let moved = translate_curve(&curve, &GpVec::new(10.0, 0.0, 0.0));
        let p = moved.d0(0.0);
        assert!((p.x() - 10.0).abs() < 1e-9);
        assert!(p.y().abs() < 1e-9);
        assert!(p.z().abs() < 1e-9);
    }
}
