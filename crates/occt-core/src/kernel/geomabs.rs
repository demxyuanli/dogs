//! Geometry abstraction enums. Source: `GeomAbs/`
//! These enums classify curves, surfaces, and geometric properties.

/// Type of parametric curve. Source: `GeomAbs_CurveType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveType { Line, Circle, Ellipse, Hyperbola, Parabola, BezierCurve, BSplineCurve, OffsetCurve, OtherCurve }

/// Type of parametric surface. Source: `GeomAbs_SurfaceType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceType { Plane, Cylinder, Cone, Sphere, Torus, BezierSurface, BSplineSurface, SurfaceOfRevolution, SurfaceOfExtrusion, OffsetSurface, OtherSurface }

/// Continuity class. Source: `GeomAbs_Shape`
///
/// Numeric order matches OCCT (`GeomAbs_C0`=0 … `GeomAbs_CN`=6).
///
/// **Known deviation (task T-64)**: OCCT returns `GeomAbs_CN` for every
/// analytic curve/surface (`Geom_Conic.cxx:32-35`, `Geom_Line.cxx:128-131`,
/// `Geom_ElementarySurface.cxx:25` for plane/cylinder/cone/sphere/torus,
/// `Geom2d_Line.cxx:172`, `Geom2d_Conic.cxx:48`), while this port's analytic
/// classes still return `G2` (3). Raising them alone empties the Sphere mesh in
/// `step_to_obj` (a consumer treats the reported continuity as a span bound),
/// so the change must be done together with that consumer's OCCT semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Shape { C0, G1, C1, G2, C2, C3, CN }

/// Join type for curve connections. Source: `GeomAbs_JoinType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JoinType { Arc, Tangent, Intersection }

/// Curve form (open or closed). Source: `GeomAbs_CurveForm`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveForm { Open, Closed, Periodic }

/// B-spline knot distribution. Source: `GeomAbs_BSplKnotDistribution`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BSplKnotDistribution { NonUniform, Uniform, QuasiUniform, PiecewiseBezier }

/// Iso-line type. Source: `GeomAbs_IsoType`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsoType { IsoU, IsoV, None }
