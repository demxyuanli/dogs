//! 2D `ElCLib` evaluators used by the `GeomBndLib_*` boxes.
//!
//! The port moved to `occt_geom::geom_bnd_lib::elclib2d`; this module re-exports
//! it so the crate's existing call sites stay unchanged.
pub use occt_geom::geom_bnd_lib::elclib2d::*;
