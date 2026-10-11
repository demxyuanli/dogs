//! ShapeUpgrade tools used by the IGES reader.
//!
//! `C0ShapeToC1Shape` (`ShapeAlgo_AlgoContainer.cxx:390-399`).

mod divide_continuity;
mod split_surface;

pub use divide_continuity::correct_c0_basis;
pub use split_surface::C0Basis;
