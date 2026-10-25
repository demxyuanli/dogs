//! `TKGeomAlgo/LocalAnalysis` — local continuity analysis between two
//! curves / surfaces.
//!
//! Source: `src/ModelingAlgorithms/TKGeomAlgo/LocalAnalysis/`.
//!
//! | OCCT | here |
//! |---|---|
//! | `LocalAnalysis_SurfaceContinuity` (`SurfaceContinuity.hxx:33-241`) | [`surface_continuity::SurfaceContinuity`] |
//! | `LocalAnalysis_StatusErrorType` (`StatusErrorType.hxx:21-29`) | [`surface_continuity::StatusErrorType`] |

pub mod surface_continuity;

pub use surface_continuity::{StatusErrorType, SurfaceContinuity};
