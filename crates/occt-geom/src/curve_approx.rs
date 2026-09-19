//! Curve approximation to polylines.
//!
//! **Provenance**: the previous header cited `GCPnts_UniformDeflection.hxx` and
//! `GeomConvert_CurveToPolyline.hxx`; the latter class **does not exist** in
//! OCCT 8.0.0 (no `*CurveToPolyline*` file in the source tree) and
//! `GCPnts_UniformDeflection` is **UNPORTED** (`occt-core/src/gcpnts.rs` module
//! docs). The polyline is therefore produced by the faithful
//! `GCPnts_TangentialDeflection` engine with its angular term disabled
//! (`angular_deflection = PI`, the same convention as
//! `meshing::edge_discret::CurveTessellator::from_range`), which is what OCCT's
//! mesh pipeline uses to tessellate curves.

use crate::curve::Curve;
use occt_core::gcpnts::{perform_tangential_curve, polyline_length, CurveSample, CurveSecondDeriv, UniformPoints};
use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};

/// Adapts `&dyn Curve` to the object-safe `CurveSample` sampler trait.
struct CurveAdapter<'a>(&'a dyn Curve);

impl CurveSample for CurveAdapter<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.0.value(u)
    }
}

/// Adapts `&dyn Curve` to the `GCPnts_TangentialDeflection` engine.
struct TdCurve<'a>(&'a dyn Curve);

impl CurveSecondDeriv for TdCurve<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.0.value(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.0.d2(u)
    }
}

/// Adaptive polyline approximation of `c` over its parameter range so that
/// every chord deviates from the curve by at most `tol`. Curves with a
/// non-finite parameter range (e.g. lines) fall back to uniform sampling of
/// `[0, 1]` of the natural parameter.
pub fn curve_to_polyline(c: &dyn Curve, tol: f64) -> Vec<GpPnt> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() {
        return curve_to_polyline_uniform(c, 64);
    }
    let mut intervals = c.parameter_intervals(6);
    if intervals.len() < 2 {
        intervals = vec![a, b];
    }
    let degree_min_nb = c.nurbs_degree().map(|d| (d + 1).max(2)).unwrap_or(2);
    let (_, points) = perform_tangential_curve(
        &TdCurve(c),
        a,
        b,
        std::f64::consts::PI,
        tol,
        2,
        PCONFUSION,
        CONFUSION,
        &intervals,
        degree_min_nb,
    );
    if points.len() >= 2 {
        points
    } else {
        vec![c.value(a), c.value(b)]
    }
}

/// `n` points sampled uniformly across the curve's parameter range
/// (`[0, 1]` for curves with a non-finite range).
pub fn curve_to_polyline_uniform(c: &dyn Curve, n: usize) -> Vec<GpPnt> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let (a, b) = finite_param(a, b);
    UniformPoints::from_curve(&CurveAdapter(c), a, b, n.max(2)).points
}

/// Distance of `mid_pt` from the chord `[c(a), c(b)]`.
pub fn adaptive_chord_error(c: &dyn Curve, a: f64, b: f64, mid_pt: &GpPnt) -> f64 {
    let pa = c.value(a);
    let pb = c.value(b);
    point_segment_dist(mid_pt, &pa, &pb)
}

/// Length of the polyline produced by [`curve_to_polyline`] at `tol`.
pub fn curve_approx_length(c: &dyn Curve, tol: f64) -> f64 {
    polyline_length(&curve_to_polyline(c, tol))
}

fn finite_param(a: f64, b: f64) -> (f64, f64) {
    if a.is_finite() && b.is_finite() {
        (a, b)
    } else {
        (0.0, 1.0)
    }
}

/// Distance from point `p` to the line segment `a..b`.
fn point_segment_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
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
    use crate::{GeomCircle, GeomLine};
    use crate::trimmed::GeomTrimmedCurve;
    use occt_core::gp::{GpAx2, GpCirc, GpDir, GpPnt};
    use std::sync::Arc;

    #[test]
    fn line_endpoints_and_length() {
        let line = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let trimmed = GeomTrimmedCurve::new(Arc::new(line), 0.0, 3.0);
        let pl = curve_to_polyline(&trimmed, 1e-6);
        assert_eq!(pl.first().unwrap().distance(&GpPnt::new(0., 0., 0.)), 0.0);
        assert_eq!(pl.last().unwrap().distance(&GpPnt::new(3., 0., 0.)), 0.0);
        let len = curve_approx_length(&trimmed, 1e-6);
        assert!((len - 3.0).abs() < 1e-9, "len={len}");
    }

    #[test]
    fn circle_semicircle() {
        let circ = GeomCircle::new(GpCirc::new(GpAx2::standard(), 2.0));
        let pi = std::f64::consts::PI;
        let trimmed = GeomTrimmedCurve::new(Arc::new(circ), 0.0, pi);
        let pl = curve_to_polyline(&trimmed, 1e-3);
        assert!(pl.len() >= 2);
        let first = pl.first().unwrap();
        let last = pl.last().unwrap();
        let on_p = |p: &GpPnt| {
            p.distance(&GpPnt::new(2., 0., 0.)) < 1e-6 || p.distance(&GpPnt::new(-2., 0., 0.)) < 1e-6
        };
        assert!(on_p(first), "first endpoint {first:?}");
        assert!(on_p(last), "last endpoint {last:?}");
        // Arc length of a radius-2 semicircle is 2*pi.
        let len = curve_approx_length(&trimmed, 1e-3);
        assert!((len - 2.0 * pi).abs() < 0.05, "len={len}");
    }

    #[test]
    fn uniform_sampling() {
        let line = GeomLine::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let trimmed = GeomTrimmedCurve::new(Arc::new(line), 0.0, 2.0);
        let pl = curve_to_polyline_uniform(&trimmed, 11);
        assert_eq!(pl.len(), 11);
        assert_eq!(pl[5].distance(&GpPnt::new(1., 0., 0.)), 0.0);
    }
}
