//! `Approx_BSplineApproxInterp` — constrained least-squares B-spline
//! approximation with exact interpolation (and optional C0 kink) constraints.
//!
//! Source: `ModelingData/TKGeomBase/Approx/Approx_BSplineApproxInterp.cxx:1-800`
//! and `Approx_BSplineApproxInterp.hxx:28-236` (OCCT 8.0.0).
//!
//! The solver builds the KKT saddle-point system of `hxx:36-48`
//!
//! ```text
//! | A^T*A  C^T  L^T |   | x |   | A^T*b |
//! | C       0    0  | * | l | = |   d   |
//! | L       0    0  |   | m |   |   0   |
//! ```
//!
//! and factorises it once with `MathLin::LU` (`MathLin_Gauss.hxx:44-155`), then
//! solves for X, Y and Z with the shared forward/back substitution of
//! `cxx:592-625`.
//!
//! Translation notes:
//! - `NCollection_Array1<gp_Pnt>` (1..N) and `NCollection_Array1<double>`
//!   (1..N) are `Vec`s indexed `0..N`; every `Value(theParams.Lower() + i)`
//!   becomes `params[i]`. `myPoints` keeps OCCT's `Array1(1, max(Length, 1))`
//!   size by padding an empty input with one default point.
//! - `NCollection_DynamicArray<int>` (`myInterpolated`, `myApproximated`,
//!   `myKinks`) is a `Vec<usize>` of 0-based point indices, exactly what OCCT
//!   appends.
//! - `GeomAdaptor_Curve` is the curve itself: `EvalD0` / `EvalD2` are
//!   [`Curve::d0`] / [`Curve::d2`]. The adaptor cache (`cxx:205`,
//!   `cxx:652-668`) therefore has no counterpart and is dropped.
//! - `MathLin::LU(aLHS, myMinPivot)` is `MathLin_Gauss.hxx:51-155` verbatim,
//!   down to the scaled implicit pivoting and the `myMinPivot` singularity
//!   test, so [`BSplineApproxInterp::set_min_pivot`] keeps its meaning.
//! - `BSplCLib::EvalBsplineBasis` is
//!   [`occt_core::bspl::eval_basis::eval_bspline_basis`];
//!   `BSplCLib::KnotSequenceLength` / `KnotSequence` are
//!   [`occt_core::bspl::knots::knot_sequence_length`] and
//!   [`occt_core::bspl::banded_interp::knot_sequence`].

mod constants;
mod linear_algebra;
mod matrix;
mod solver;

pub use solver::BSplineApproxInterp;
