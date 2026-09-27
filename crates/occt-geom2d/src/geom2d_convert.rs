//! `Geom2dConvert` subset used by `ShapeBuild_Edge::TransformPCurve`
//! (`Geom2dConvert.cxx`).
//!
//! `CurveToBSplineCurve` is exact for a `Geom2d_BezierCurve` (same poles with
//! the clamped knot vector `[0 x (degree+1), 1 x (degree+1)]`) and the
//! identity for a `Geom2d_BSplineCurve`.
//!
//! UNPORTED: the `Geom2d_Conic` arm (`Convert_QuasiAngular`) and
//! `Geom2dConvert_ApproxCurve` (`Geom2dConvert_ApproxCurve.cxx`); callers that
//! hit them must mark the branch UNPORTED.

use crate::bspline_curve::Geom2dBSplineCurve;
use crate::curve::Curve2d;

/// `Geom2dConvert::CurveToBSplineCurve(C)` (`Geom2dConvert.cxx`), the exact
/// non-conic cases: Bezier -> BSpline with the identical poles, BSpline ->
/// itself. Returns `None` for a conic or any other curve type (UNPORTED).
pub fn curve_to_bspline_curve(c: &dyn Curve2d) -> Option<Box<dyn Curve2d>> {
    if c.is_bspline2d() {
        return Some(c.clone_dyn());
    }
    if c.is_bezier2d() {
        let poles = c.poles2d()?;
        let degree = poles.len().checked_sub(1)?;
        let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
        let mut knots = vec![0.0f64; degree + 1];
        knots.extend(std::iter::repeat(1.0f64).take(degree + 1));
        let bs = Geom2dBSplineCurve::new(xs, ys, knots, degree).ok()?;
        return Some(Box::new(bs));
    }
    None
}
