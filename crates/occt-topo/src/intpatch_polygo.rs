//! `IntPatch_Polygo` - the abstract base of the `IntPatch` 2D polygons.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_Polygo.hxx`,
//! `.lxx` and `.cxx`.
//!
//! `IntPatch_Polygo` derives from `Intf_Polygon2d` and only adds the stored
//! `myError` (`hxx:55`) plus the `NbSegments` / `Segment` defaults built on the
//! still-virtual `NbPoints` / `Point` pair. The port keeps that shape as
//! [`IntPatchPolygo`]: `Error` / `NbPoints` / `Point` live here, and the three
//! `Intf_Polygon2d` members are implemented per class from them - OCCT's
//! `DeflectionOverEstimation` (`lxx:24-27`) returns `myError`, `NbSegments`
//! (`lxx:31-34`) returns `NbPoints() - 1`, and `Segment(i)`
//! (`lxx:38-43`) returns `(Point(i), Point(i + 1))`.
//!
//! `Dump` (`cxx:28-47`) is an OCCT debug dump with no caller in the port and is
//! not translated.

use occt_core::gp::GpPnt2d;
use occt_core::intf::IntfPolygon2d;

/// `IntPatch_Polygo` (`hxx:27-57`).
pub trait IntPatchPolygo: IntfPolygon2d {
    /// `myError` (`hxx:55`) - `Error()` (`lxx:18-21`).
    fn error(&self) -> f64;

    /// `NbPoints()` (`hxx:34`), 1-based indices.
    fn nb_points(&self) -> usize;

    /// `Point(Index)` (`hxx:36`), 1-based.
    fn point(&self, index: usize) -> GpPnt2d;
}
