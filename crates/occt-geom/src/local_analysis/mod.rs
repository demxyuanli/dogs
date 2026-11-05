//! `TKGeomAlgo/LocalAnalysis` — local continuity analysis between two
//! curves / surfaces.
//!
//! Source: `src/ModelingAlgorithms/TKGeomAlgo/LocalAnalysis/`.
//!
//! The package belongs to TKGeomAlgo, so it lives in `occt-geom`:
//! `GeomPlate_BuildPlateSurface::EcartContraintesMil`
//! (`GeomPlate_BuildPlateSurface.cxx:833-863`) and
//! `GeomPlate_BuildPlateSurface::VerifPoints` (`:2771-2779`) both call
//! `LocalAnalysis_SurfaceContinuity`.
//!
//! | OCCT | here |
//! |---|---|
//! | `LocalAnalysis_SurfaceContinuity` (`SurfaceContinuity.hxx:33-241`) | [`surface_continuity::SurfaceContinuity`] |
//! | `LocalAnalysis_StatusErrorType` (`StatusErrorType.hxx:21-29`) | [`surface_continuity::StatusErrorType`] |

pub mod surface_continuity;

pub use surface_continuity::{StatusErrorType, SurfaceContinuity};
