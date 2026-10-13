//! Edge to face UV curve (pcurve) projection onto analytic and general surfaces.
//!
//! This is the former monolithic `surface_projector.rs`, split by OCCT
//! control-flow region. Every item keeps the OCCT citation it already carried:
//!
//! - [`classify`] - surface-type classification plus the sampled analytic frames
//!   (`cylinder_params` / `cone_params` / `torus_params`) and their helpers.
//! - [`projector`] - [`Projector`], the analytic per-type UV projection of a 3D
//!   point or curve, and the [`project_point_on_surface`] entry point.
//! - [`analytic`] - the `generateCurvePoints` sample count and the basis-surface
//!   unwrap behind the `ProjLib` projection route.
//! - [`newton_uv`] - `ShapeAnalysis_Surface::ValueOfUV` with the `SurfaceNewton`
//!   polish and `NextValueOfUV`.
//! - [`closed_state`] - `RestrictBounds` and the `IsUClosed` / `IsVClosed` family.
//! - [`iso_param`] - `isAnIsoparametric` and `ShapeAnalysis_Curve::Project`.
//! - [`period_jump`] - the seam / degeneracy jump resolution on a closed surface
//!   (`needResolveUJump` / `needResolveVJump`, `AdjustOverDegen`).
//! - [`curve_projection`] - `ShapeConstruct_ProjectCurveOnSurface::Perform` and
//!   the iso / analytic p-curve builders it dispatches to.
//! - [`reparam`] - [`ReparamCurve2d`] and the helpers that align a pcurve range
//!   with the edge 3D range.
//!
//! Source: `ShapeAnalysis_Surface` and `ShapeConstruct_ProjectCurveOnSurface`.

use super::prelude::*;

mod analytic;
mod classify;
mod closed_state;
mod curve_projection;
mod iso_param;
mod newton_uv;
mod period_jump;
mod projector;
mod reparam;

pub(crate) use analytic::*;
pub use classify::*;
pub(crate) use closed_state::*;
pub use curve_projection::*;
pub(crate) use iso_param::*;
pub(crate) use newton_uv::*;
use period_jump::*;
pub use projector::*;
pub use reparam::*;
