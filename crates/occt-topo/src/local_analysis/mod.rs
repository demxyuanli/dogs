//! `TKGeomAlgo/LocalAnalysis` — the code lives in `occt-geom`
//! (`occt-geom/src/local_analysis/`), because `LocalAnalysis` is a TKGeomAlgo
//! package and `GeomPlate_BuildPlateSurface` is its only consumer
//! (`GeomPlate_BuildPlateSurface.cxx:833-863`, `:2771-2779`).
//!
//! The re-export keeps the `crate::local_analysis::...` paths of the
//! `occt-topo` callers (`chfi3d/tangent_faces.rs`,
//! `shape_upgrade/divide_continuity.rs`) unchanged.

pub use occt_geom::local_analysis::surface_continuity;

pub use surface_continuity::{StatusErrorType, SurfaceContinuity};
