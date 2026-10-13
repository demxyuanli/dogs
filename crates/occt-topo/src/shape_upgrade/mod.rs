//! ShapeUpgrade tools used by the IGES reader.
//!
//! `C0ShapeToC1Shape` (`ShapeAlgo_AlgoContainer.cxx:390-399`).
//! `C0BSplineToSequenceOfC1BSplineCurve` (`ShapeAlgo_AlgoContainer.cxx:636-650`).

mod c0_bspline_sequence;
mod divide_continuity;
mod split_surface;

pub use c0_bspline_sequence::{
    c0_bspline2d_to_sequence_of_c1_bspline_curve, c0_bspline_to_sequence_of_c1_bspline_curve,
    GEOMABS_C1,
};
pub use divide_continuity::correct_c0_basis;
pub use split_surface::C0Basis;
