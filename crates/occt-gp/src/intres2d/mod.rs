//! Port of the `IntRes2d` package
//! (`src/ModelingAlgorithms/TKGeomAlgo/IntRes2d/FILES.cmake`): the 2D
//! intersection result and domain types shared by every 2D intersector.
//!
//! `IntRes2d_SequenceOfIntersectionPoint` and
//! `IntRes2d_SequenceOfIntersectionSegment` are `NCollection_Sequence` aliases
//! (`IntRes2d_SequenceOfIntersectionPoint.hxx` /
//! `IntRes2d_SequenceOfIntersectionSegment.hxx`); the port uses `Vec` for them,
//! which is what [`IntRes2dIntersection::lpnt`] / [`IntRes2dIntersection::lseg`]
//! already are.

pub mod domain;
pub mod intersection;
pub mod transition;

pub use domain::IntRes2dDomain;
pub use intersection::{
    IntRes2dIntersection, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment,
};
pub use transition::{
    IntRes2dPosition, IntRes2dSituation, IntRes2dTransition, IntRes2dTypeTrans,
};
