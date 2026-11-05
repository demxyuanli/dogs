//! `Plate` package (OCCT `TKGeomAlgo/Plate`): variational spline surface.
//!
//! The package solves for a two-variable function that satisfies a set of
//! point / derivative / continuity constraints while minimising an energy-like
//! criterion. It is the numeric core behind `GeomPlate_BuildPlateSurface` and
//! `BRepFill_Filling`.
//!
//! Layout (one module per OCCT class / responsibility):
//! - [`d1`] [`d2`] [`d3`]: `Plate_D1` / `Plate_D2` / `Plate_D3` derivative value
//!   containers,
//! - [`pinpoint_constraint`]: `Plate_PinpointConstraint`,
//! - [`linear_xyz_constraint`] / [`linear_scalar_constraint`]: linear
//!   combinations of pinpoint constraints,
//! - [`line_constraint`] / [`plane_constraint`] / [`sampled_curve_constraint`] /
//!   [`global_translation_constraint`]: geometric constraints wrapping the
//!   linear ones,
//! - [`gto_c_constraint`] / [`free_gto_c_constraint`]: `G1` / `G2` / `G3`
//!   continuity constraints,
//! - [`plate`]: `Plate_Plate`, the solver.

mod d1;
mod d2;
mod d3;
mod free_gto_c_constraint;
mod global_translation_constraint;
mod gto_c_constraint;
mod line_constraint;
mod linear_scalar_constraint;
mod linear_xyz_constraint;
mod pinpoint_constraint;
mod plane_constraint;
mod plate;
mod sampled_curve_constraint;

pub use d1::PlateD1;
pub use d2::PlateD2;
pub use d3::PlateD3;
pub use free_gto_c_constraint::FreeGtoCConstraint;
pub use global_translation_constraint::GlobalTranslationConstraint;
pub use gto_c_constraint::GtoCConstraint;
pub use line_constraint::LineConstraint;
pub use linear_scalar_constraint::LinearScalarConstraint;
pub use linear_xyz_constraint::LinearXyzConstraint;
pub use pinpoint_constraint::PinpointConstraint;
pub use plane_constraint::PlaneConstraint;
pub use plate::Plate;
pub use sampled_curve_constraint::SampledCurveConstraint;
