//! OCCT 2D geometry (TKG2d). Parametric curves via trait-based polymorphism.
//! Replaces OCCT's Handle(Geom2d_Curve) pattern with Arc<dyn Curve2d>.
//!
//! # Architecture
//! ```text
//! curve       (Curve2d trait: D0, D1, D2, Continuity, Transform)
//! line        (Geom2dLine: gp_Lin2d wrapper)
//! circle      (Geom2dCircle: gp_Circ2d wrapper)
//! ellipse     (Geom2dEllipse: gp_Elips2d wrapper)
//! hyperbola   (Geom2dHyperbola: gp_Hypr2d wrapper)
//! parabola    (Geom2dParabola: gp_Parab2d wrapper)
//! ```

use std::sync::Arc;

pub mod curve;
pub mod line;
pub mod circle;
pub mod ellipse;
pub mod hyperbola;
pub mod parabola;
pub mod trimmed;
pub mod offset;

pub use curve::Curve2d;
pub use line::Geom2dLine;
pub use circle::Geom2dCircle;
pub use ellipse::Geom2dEllipse;
pub use hyperbola::Geom2dHyperbola;
pub use parabola::Geom2dParabola;

/// Type alias for OCCT's Handle(Geom2d_Curve) — Arc is the Rust equivalent.
pub type HandleCurve2d = Arc<dyn Curve2d>;
