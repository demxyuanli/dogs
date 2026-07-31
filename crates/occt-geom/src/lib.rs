//! OCCT 3D geometry (TKG3d). Parametric curves + surfaces via traits.
//! Replaces OCCT's Handle(Geom_Curve) and Handle(Geom_Surface) with Arc<dyn Trait>.

use std::sync::Arc;

pub mod curve;
pub mod surface;
pub mod line;
pub mod circle;
pub mod ellipse;
pub mod hyperbola;
pub mod parabola;
pub mod plane;
pub mod cylinder;
pub mod cone;
pub mod sphere;
pub mod torus;

pub use curve::Curve;
pub use surface::Surface;
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

pub type HandleCurve = Arc<dyn Curve>;
pub type HandleSurface = Arc<dyn Surface>;
