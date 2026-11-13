//! `GeomFill_Gordon` family (OCCT `TKGeomAlgo/GeomFill`, OCCT 8.0.1).
//!
//! - [`gordon`] — `GeomFill_Gordon`, the high-level driver
//!   (`GeomFill_Gordon.cxx`).
//! - [`network_surface`] — `GeomFill_NetworkSurface`
//!   (`GeomFill_NetworkSurface.cxx`).
//! - [`utilities`] — `GordonUtilities` (`GeomFill_GordonUtilities.pxx`).

pub mod gordon;
pub mod network_surface;
pub mod utilities;

pub use gordon::{
    ApproximationMode, BuildReport, BuildStage, Gordon, ResultStatus,
};
pub use network_surface::NetworkSurface;
