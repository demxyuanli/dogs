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
/// Analytic curves and surfaces report `GeomAbs_CN` (6) like OCCT does
/// (`Geom_Conic.cxx:32-35`, `Geom_Line.cxx:128-131`,
/// `Geom_ElementarySurface.cxx:25` for plane/cylinder/cone/sphere/torus,
/// `Geom2d_Line.cxx:172`, `Geom2d_Conic.cxx:48`). Task T-64 had to hold them at
/// `G2` (3) because a mesh consumer treated the reported continuity as a span
/// bound and the Sphere mesh came out empty; that consumer now derives its spans
/// the OCCT way (`GeomAdaptor_Curve::NbIntervals` returns a single interval for
/// every non-B-spline, non-offset curve), so the values are faithful again.
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
