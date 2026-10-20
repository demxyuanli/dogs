//! BSpline conversion depth. Port of the `GeomConvert` gap-fillers
//! (TKGeomBase): `GeomConvert_BSplineCurveKnotSplitting`,
//! `GeomConvert_BSplineSurfaceKnotSplitting` and
//! `GeomConvert_BSplineSurfaceToBezierSurface`.
//!
//! Curve→Bezier is already covered by `crate::bspline_to_bezier` and is not
//! duplicated here; `GeomConvert_ApproxCurve` (adaptive spline approximation)
//! is ported in `crate::convert_approx_curve` on top of
//! `crate::adv_approx::ApproxAFunction3d` (the `AdvApprox_ApproxAFunction`
//! port) and `AdvApprox_PrefAndRec` cutting.
//!
//! `GeomConvert_CompCurveToBSplineCurve` is ported (constructor + `Add`,
//! `GeomConvert_CompCurveToBSplineCurve.cxx:32-273`); it needs
//! `Geom_BSplineCurve::IncreaseDegree` and `Geom_BSplineCurve::RemoveKnot`
//! (both ported on `bspline_curve.rs`).
//!
//! Both `Geom_OffsetCurve` arms of `GeomConvert::CurveToBSplineCurve`
//! (`GeomConvert.cxx:340-354`, `:436-450`) call `GeomConvert_ApproxCurve(C,
//! 1e-4, C2, 16, 14)`; they are wired to `crate::convert_approx_curve`.
//!
//! `GeomConvert::CurveToBSplineCurve` itself is ported for the trimmed
//! line/conic arms (including the `RationalC1` + `U2-U1>=6` split-stitch arm)
//! and the Bezier / B-spline copy arms (batch 92).

mod bezier_surface;
mod comp_curve;
mod curve_convert;
mod knot_splitting;

pub use self::bezier_surface::{bspline_surface_to_bezier_surface, BezierSurfacePatch};
pub use self::comp_curve::CompCurveToBSplineCurve;
pub use self::curve_convert::{bspline_copy, curve_to_bspline_curve, is_bspline};
pub use self::knot_splitting::{
    knot_splitting_curve, knot_splitting_curve_with_continuity, knot_splitting_surface,
    knot_splitting_surface_with_continuity,
};

#[cfg(test)]
mod tests;
