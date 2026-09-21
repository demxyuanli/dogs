//! IntTools data classes — Phase 16 Wave A1.
//!
//! Ports the `IntTools_*` data containers from TKBO
//! (`ModelingAlgorithms/TKBO/IntTools`) into self-contained Rust value types:
//!
//! - [`IntRange`] / [`IntToolsRange`] — a 1D range `[first, last]`
//!   (`IntTools_Range`).
//! - [`CommonPrt`] / [`CommonPartType`] — a common part between shapes
//!   (`IntTools_CommonPrt`).
//! - [`IntRoot`] / [`RootType`] — a root of a 1D function, kept as the
//!   interval the root was found in (`IntTools_Root`).
//! - [`PntOnFace`] / [`PntOn2Faces`] — a 3D point plus UV parameters on one /
//!   two faces (`IntTools_PntOnFace`, `IntTools_PntOn2Faces`).
//! - [`IntCurve`] / [`CurveKind`] — an intersection-curve container
//!   (`IntTools_Curve`).
//! - [`MarkedRangeSet`] — a flagged interval set with union / intersection /
//!   difference (`IntTools_MarkedRangeSet`).
//! - [`LocalizeData`] / [`LocalizeData2`] — local edge/face intersection
//!   bookkeeping (face range, edge ranges, UV point lists).
//!
//! These are pure data classes: they carry geometry-adjacent numbers and shape
//! references but perform no surface/curve evaluation themselves. The module is
//! self-contained — it depends only on `GpPnt` (coordinates) and `TopoShape`
//! (shape handles), never on sibling Phase 16 modules.
mod prelude {

pub(crate) use crate::shape::TopoShape;
pub(crate) use occt_core::gp::GpPnt;

}


mod types;
mod localize_data;
pub use types::*;
pub use localize_data::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
