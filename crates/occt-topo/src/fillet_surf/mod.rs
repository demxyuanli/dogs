//! `FilletSurf` — geometric information about a constant-radius fillet.
//!
//! Sources (OCCT 8.0.0), `ModelingAlgorithms/TKFillet/FilletSurf/`:
//! - `FilletSurf_Builder.cxx:1-387` — the public API ([`FilletSurfBuilder`]),
//! - `FilletSurf_InternalBuilder.cxx:1-879` — the private worker
//!   ([`FilletSurfInternalBuilder`]),
//! - `FilletSurf_StatusDone.hxx`, `FilletSurf_ErrorTypeStatus.hxx`,
//!   `FilletSurf_StatusType.hxx` — the three status enumerations.
//!
//! The package is the only consumer of `ChFi3d_FilBuilder` that drives the
//! builder one `PerformSetOfKGen` at a time instead of going through
//! `ChFi3d_Builder::Compute`, and the only one that overrides the face / face
//! `PerformSurf` to keep the walked line for its section queries. It closes the
//! `BlendFunc_ConstRad` → `ChFi3d_FilBuilder` → `FilletSurf` reference loop.

mod builder;
mod internal_builder;
mod status;

pub use builder::FilletSurfBuilder;
pub use internal_builder::edge_parameter::{
    compute_edge_parameter, int_plan_edge, unbounded_plane_window, UNBOUNDED_PLANE_WINDOW,
};
pub use internal_builder::ops::FilletSurfOps;
pub use internal_builder::FilletSurfInternalBuilder;
pub use status::{ErrorTypeStatus, StatusDone, StatusType};
