//! Revolution sweep used when an IGES 120 generatrix is not a single curve.
//!
//! `BRepPrimAPI_MakeRevol` -> `BRepSweep_Revol` -> `BRepSweep_Rotation`.

mod continuity;
mod grid;
mod quilt;
mod num_shape;
mod revol;
mod revolution_surface;
mod rotation;
mod translation;

pub use revol::make_revol;
pub use translation::sweep_make_prism;
