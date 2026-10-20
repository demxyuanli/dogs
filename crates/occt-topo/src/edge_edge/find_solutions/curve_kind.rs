use occt_geom::Curve;

// ---------------------------------------------------------------------------
// Curve type (`GeomAdaptor_Curve::GetType`)
// ---------------------------------------------------------------------------

/// `GeomAbs_CurveType` — the arms `IntTools_EdgeEdge` distinguishes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum GeomCurveKind {
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

/// `GeomAdaptor_Curve::GetType()` for a port curve.
///
/// The order is `GeomAdaptor_Curve::load` (`GeomAdaptor_Curve.cxx:252-311`)
/// verbatim: trimmed curves are unwrapped to their basis (the port's trimmed
/// types forward these queries), then Circle, Line, Ellipse, Parabola,
/// Hyperbola, Bezier, BSpline, Offset, Other. Only the concrete class decides —
/// a geometrically straight B-spline is `GeomAbs_BSplineCurve` in OCCT, so
/// `Curve::is_line` (a geometric test) is deliberately **not** used here.
pub(crate) fn geom_curve_kind(c: &dyn Curve) -> GeomCurveKind {
    if c.gp_circ().is_some() {
        GeomCurveKind::Circle
    } else if c.gp_line().is_some() {
        GeomCurveKind::Line
    } else if c.gp_ellipse().is_some() {
        GeomCurveKind::Ellipse
    } else if c.gp_parabola().is_some() {
        GeomCurveKind::Parabola
    } else if c.gp_hyperbola().is_some() {
        GeomCurveKind::Hyperbola
    } else if c.bezier_poles().is_some() {
        GeomCurveKind::BezierCurve
    } else if c.bspline_knots().is_some() || c.nurbs_degree().is_some() {
        GeomCurveKind::BSplineCurve
    } else if c.offset_curve().is_some() {
        GeomCurveKind::OffsetCurve
    } else {
        GeomCurveKind::OtherCurve
    }
}

/// `TypeToInteger` (`IntTools_EdgeEdge.cxx:1456-1482`).
pub(crate) fn curve_kind_to_integer(kind: GeomCurveKind) -> i32 {
    match kind {
        GeomCurveKind::Line => 0,
        GeomCurveKind::Hyperbola | GeomCurveKind::Parabola => 1,
        GeomCurveKind::Circle | GeomCurveKind::Ellipse => 2,
        GeomCurveKind::BezierCurve | GeomCurveKind::BSplineCurve => 3,
        _ => 4,
    }
}
