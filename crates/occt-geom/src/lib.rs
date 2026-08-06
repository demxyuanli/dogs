//! OCCT 3D geometry (TKG3d). Parametric curves + surfaces via traits.
//! Replaces OCCT's Handle(Geom_Curve) and Handle(Geom_Surface) with Arc<dyn Trait>.

use std::sync::Arc;

pub mod curve;
pub mod bspline_curve;
pub mod interp_curve;
pub mod surface;
pub mod line;
pub mod circle;
pub mod ellipse;
pub mod hyperbola;
pub mod parabola;
pub mod trimmed;
pub mod revolved;
pub mod offset;
pub mod plane;
pub mod cylinder;
pub mod cone;
pub mod sphere;
pub mod torus;

pub use curve::Curve;
pub use bspline_curve::GeomBSplineCurve;
pub use surface::Surface;
// Phase 3 modules (curve/surface conversion).
pub mod bspline_to_bezier;
pub mod bezier_to_polyline;
pub mod curve_approx;
pub mod surface_to_grid;
pub mod convert_geom;
// Phase 4 modules.
pub mod geom_api;
pub mod curve_reparam;
pub mod surface_fit;
pub mod surface_ops;
pub use line::GeomLine;
pub use circle::GeomCircle;
pub use ellipse::GeomEllipse;
pub use hyperbola::GeomHyperbola;
pub use parabola::GeomParabola;
pub use plane::GeomPlane;
pub use cylinder::GeomCylinder;
pub use cone::GeomCone;
pub use sphere::GeomSphere;
pub use torus::GeomTorus;
pub use trimmed::GeomTrimmedCurve;
pub use revolved::GeomRevolvedSurface;
pub use surface_of_revolution::GeomSurfaceOfRevolution;
pub use offset::GeomOffsetCurve;
pub mod offset_surface;
pub use offset_surface::GeomOffsetSurface;
pub mod extrema;
pub mod extrema_pc;
pub mod bspline_surface;
// Phase 13 (TKGeomAlgo depth): wave 2-4 placeholders.
pub mod extrema_cc;
pub mod extrema_surf;
pub mod extrema_ss;
pub mod intana;
pub mod gcpnts;
pub mod convert_bspl;
pub mod surface_of_revolution;

pub type HandleCurve = Arc<dyn Curve>;
pub type HandleSurface = Arc<dyn Surface>;
