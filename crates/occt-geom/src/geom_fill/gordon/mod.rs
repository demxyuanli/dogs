//! `GeomFill_Gordon` family (OCCT `TKGeomAlgo/GeomFill`).
//!
//! - [`gordon`] — `GeomFill_Gordon`, the high-level driver
//!   (`GeomFill_Gordon.cxx`).
//! - [`reparam_curve`] — the anonymous-namespace helpers of that file
//!   (`GeomFill_Gordon.cxx:33-429`).
//! - [`builder`] — `GeomFill_GordonBuilder`, the surface builder
//!   (`GeomFill_GordonBuilder.cxx`).

pub mod builder;
pub mod gordon;
pub mod reparam_curve;

pub use builder::GordonBuilder;
pub use gordon::Gordon;
pub use reparam_curve::{invert_reparam_bisection, kink_parameters, linspace_with_breaks, reparam_curve};
