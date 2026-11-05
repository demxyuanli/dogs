//! `GeomPlate` package — the plate solving front end.
//!
//! Source: `src/ModelingAlgorithms/TKGeomAlgo/GeomPlate/` (OCCT 8.0.0).
//!
//! Layout:
//! - [`Aij`]: `GeomPlate_Aij` (`GeomPlate_Aij.cxx:22-32`),
//! - [`PointConstraint`]: `GeomPlate_PointConstraint`
//!   (`GeomPlate_PointConstraint.cxx:29-215`),
//! - [`CurveConstraint`]: `GeomPlate_CurveConstraint`
//!   (`GeomPlate_CurveConstraint.cxx:39-434`),
//! - [`PlateSurface`]: `GeomPlate_Surface`
//!   (`GeomPlate_Surface.cxx:42-309`),
//! - [`BuildAveragePlane`]: `GeomPlate_BuildAveragePlane`
//!   (`GeomPlate_BuildAveragePlane.cxx:47-750`),
//! - [`PlateG0Criterion`] / [`PlateG1Criterion`]:
//!   `GeomPlate_PlateG0Criterion` (`GeomPlate_PlateG0Criterion.cxx:33-122`)
//!   and `GeomPlate_PlateG1Criterion`
//!   (`GeomPlate_PlateG1Criterion.cxx:31-133`).

mod aij;
mod build_average_plane;
mod build_plate_surface;
mod criteria;
mod curve_constraint;
mod make_approx;
mod point_constraint;
mod surface;

pub use aij::Aij;
pub use build_average_plane::{half_space, BuildAveragePlane};
pub use build_plate_surface::{BuildPlateSurface, PerformOutcome};
pub use criteria::{PlateG0Criterion, PlateG1Criterion};
pub use curve_constraint::{BoundaryCurve, CurveConstraint};
pub use make_approx::MakeApprox;
pub use point_constraint::PointConstraint;
pub use surface::PlateSurface;
