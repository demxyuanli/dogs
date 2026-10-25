//! `GeomBndLib_Curve.cxx` for 3D curves (`Box` only).
//!
//! The port moved to `occt_geom::geom_bnd_lib::curve3d` (TKGeomBase belongs to
//! `occt-geom`; see that module for the OCCT source references). This module
//! re-exports it so the crate's existing `crate::geom_bnd_lib_curve3d::box_curve`
//! call sites stay unchanged.
pub use occt_geom::geom_bnd_lib::curve3d::*;
