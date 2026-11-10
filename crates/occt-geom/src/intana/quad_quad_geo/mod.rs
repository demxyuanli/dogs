//! `IntAna_QuadQuadGeo` analytic arms (cylinder/cylinder, cylinder/sphere,
//! sphere/cone, cone/cone, cylinder/cone), the shared `AxeOperator` helpers and
//! `IntAna_IntLinTorus`.
//!
//! Source: `ModelingAlgorithms/TKGeomBase/IntAna/IntAna_QuadQuadGeo.cxx` and
//! `IntAna_IntLinTorus.cxx` (OCCT 8.0.0). The torus-pair arms live beside this
//! module in `intana_torus.rs`; the plane/plane, plane/sphere, plane/cylinder
//! and plane/cone arms live in `intana/analytic_intersections.rs`.

use super::*;

mod axes;
mod cone_pairs;
mod cylinder_pairs;
mod dispatch;
mod int_lin_torus;
mod sphere_cone;

pub use axes::*;
pub use cone_pairs::*;
pub use cylinder_pairs::*;
pub use dispatch::*;
pub use int_lin_torus::*;
pub use sphere_cone::*;
