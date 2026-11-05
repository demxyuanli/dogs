//! Revolution sweep used when an IGES 120 generatrix is not a single curve.
//!
//! `BRepPrimAPI_MakeRevol` -> `BRepSweep_Revol` -> `BRepSweep_Rotation`.

mod continuity;
mod grid;
mod quilt;
mod num_shape;
mod prism;
mod revol;
mod revolution_surface;
mod rotation;
mod translation;

pub use prism::Prism;
pub use quilt::Quilt;
pub use revol::{make_revol, Revol};
pub use translation::sweep_make_prism;

pub(crate) use revolution_surface::revolution_ax3;

