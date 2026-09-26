//! IntTools interval sampling and topology tools — Phase 16b.
//!
//! Ports the TKBO `IntTools` sample/range-localization classes into
//! self-contained Rust value types:
//!
//! - [`BaseRangeSample`] — base class for range-index management
//!   (`IntTools_BaseRangeSample`).
//! - [`CurveRangeSample`] — a curve parameter sub-range addressed by
//!   `(depth, index)`; `get_range` materialises the `[first, last]` bounds
//!   (`IntTools_CurveRangeSample`).
//! - [`SurfaceRangeSample`] — the 2D analogue: a `(U, V)` cell addressed by
//!   `(depth_u, index_u, depth_v, index_v)` (`IntTools_SurfaceRangeSample`).
//! - [`CurveRangeLocalizeData`] / [`SurfaceRangeLocalizeData`] — split a
//!   curve/UV domain by sample points and map each cell to a curve/surface
//!   index, tracking which cells are already known "out"
//!   (`IntTools_CurveRangeLocalizeData`, `IntTools_SurfaceRangeLocalizeData`).
//! - [`TopolTool`] — the sample-point generator for intersection algorithms:
//!   computes a uniform `U × V` grid (`ComputeSamplePoints`), answers
//!   `sample_point` (UV + evaluated 3D point), and builds an
//!   deflection-adaptive grid for BSpline surfaces (`SamplePnts`)
//!   (`IntTools_TopolTool`).
//!
//! The range samples reuse [`IntRange`](crate::inttools_data::IntRange) from
//! the Phase 16a data port. The module depends only on `occt_geom::Surface`
//! (via `dyn Surface`) and the surface-type classifier reused from the BRepMesh
//! range-splitter port — never on sibling Phase 16 modules.
//!
//! `ponytail:` `IntTools_TopolTool` in OCCT reads the concrete `Geom_*` surface
//! type, radii and pole/knot counts through `Adaptor3d_Surface` accessors. The
//! Rust port only sees `dyn Surface`, so the analytic type is recovered with
//! [`classify_surface`], radii are measured by sampling diametrically opposite
//! iso-parameter points, and BSpline/Bezier pole counts (unreachable through the
//! trait object) fall back to a fixed 10×10 base grid that `sample_pnts` then
//! refines adaptively.
mod prelude {

pub(crate) use occt_core::gp::{GpPnt, GpPnt2d};
pub(crate) use occt_core::precision::ANGULAR;
pub(crate) use occt_geom::Surface;

pub(crate) use crate::inttools_data::IntRange;
pub(crate) use crate::meshing::range_splitter::{classify_surface, SurfaceType};

}


mod samples;
mod surface_sample;
pub use samples::*;


#[cfg(test)]
#[path = "tests.rs"]
mod tests;
