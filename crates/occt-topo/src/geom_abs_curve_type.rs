//! `GeomAbs_CurveType` and the `GeomAdaptor_Curve::load` dispatch.
//!
//! Source: `TKG3d/GeomAdaptor/GeomAdaptor_Curve.cxx:238-313` (the
//! `DynamicType()` chain that sets `myTypeCurve`), `GeomAbs_CurveType.hxx`.
//!
//! OCCT keeps the type in the adaptor (`myTypeCurve`) and `GetType()` returns
//! it; the port's curves carry no such field, so the type is recomputed from
//! the same concrete-class queries in the same order.

use occt_geom::curve::Curve;
use occt_geom2d::curve::Curve2d;

/// `GeomAbs_CurveType` as far as the dispatchers of this crate read it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsCurveType {
    Line,
    Circle,
    Ellipse,
    Hyperbola,
    Parabola,
    BezierCurve,
    BSplineCurve,
    OffsetCurve,
    OtherCurve,
}

/// `GeomAdaptor_Curve::load`'s type chain (`GeomAdaptor_Curve.cxx:238-313`):
/// `Geom_TrimmedCurve` recurses into its basis curve, then Circle, Line,
/// Ellipse, Parabola, Hyperbola, Bezier, BSpline, Offset, else `OtherCurve`.
pub fn curve_kind(c: &dyn Curve) -> GeomAbsCurveType {
    if let Some((basis, _, _)) = c.untrimmed_basis() {
        return curve_kind(basis.as_ref());
    }
    if c.gp_circ().is_some() {
        GeomAbsCurveType::Circle
    } else if c.gp_line().is_some() {
        GeomAbsCurveType::Line
    } else if c.gp_ellipse().is_some() {
        GeomAbsCurveType::Ellipse
    } else if c.gp_parabola().is_some() {
        GeomAbsCurveType::Parabola
    } else if c.gp_hyperbola().is_some() {
        GeomAbsCurveType::Hyperbola
    } else if c.bezier_poles().is_some() {
        GeomAbsCurveType::BezierCurve
    } else if c.bspline_poles().is_some() {
        GeomAbsCurveType::BSplineCurve
    } else if c.offset_curve().is_some() {
        GeomAbsCurveType::OffsetCurve
    } else {
        GeomAbsCurveType::OtherCurve
    }
}

/// `Geom2dAdaptor_Curve::load`'s type chain (`Geom2dAdaptor_Curve.cxx:1069-1143`):
/// `Geom2d_TrimmedCurve` recurses into its basis curve, then Circle, Line,
/// Ellipse, Parabola, Hyperbola, BSpline, Offset, else `OtherCurve`.
///
/// `Geom2d_BezierCurve` is the one link the port cannot test: `Curve2d` carries
/// no Bezier query (`crates/occt-geom2d/src/curve.rs:236` exposes
/// `bspline_poles2d` only), so a Bezier p-curve reports `BSplineCurve`.
/// `FUN_tool_quad` / `FUN_tool_line` answer `false` for that type either way.
pub fn curve2d_kind(c: &dyn Curve2d) -> GeomAbsCurveType {
    if let Some(basis) = c.trimmed_basis() {
        return curve2d_kind(basis);
    }
    if c.gp_circ2d().is_some() {
        GeomAbsCurveType::Circle
    } else if c.gp_lin2d().is_some() {
        GeomAbsCurveType::Line
    } else if c.gp_elips2d().is_some() {
        GeomAbsCurveType::Ellipse
    } else if c.gp_parab2d().is_some() {
        GeomAbsCurveType::Parabola
    } else if c.gp_hypr2d().is_some() {
        GeomAbsCurveType::Hyperbola
    } else if c.bspline_poles2d().is_some() {
        GeomAbsCurveType::BSplineCurve
    } else if c.offset_basis().is_some() {
        GeomAbsCurveType::OffsetCurve
    } else {
        GeomAbsCurveType::OtherCurve
    }
}
