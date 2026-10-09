//! 2D extrema — point/curve minimum and maximum distances in the plane.
//!
//! Port of `Extrema_ExtPC2d`, `Extrema_ExtCC2d` / `Extrema_ECC2d` and the
//! general engine behind them. Source: `Extrema` (TKGeomBase).
//!
//! `curve2d_intersections` now runs the faithful general route:
//! `Extrema_ExtCC2d` (`ext_cc2d`) -> `Extrema_ECC2d` = `Extrema_GGenExtCC`
//! (`general_extrema`) -> `Extrema_GlobOptFuncCCC2` (`glob_opt_func`) ->
//! `math_GlobOptMin` (`occt_math::globoptmin`). See `curve_ops.rs`.

mod prelude {

pub(crate) use std::cmp::Ordering;

pub(crate) use occt_core::gp::{
    GpAx22d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d,
};
pub(crate) use occt_core::precision::{
    Precision, ANGULAR, CONFUSION, INFINITE, PCONFUSION, RESOLUTION, SQUARE_CONFUSION,
};

pub(crate) use crate::curve::Curve2d;

}


mod analytic_solvers;
mod curve_curve;
mod curve2d_tool;
mod general_extrema;
mod glob_opt_func;
mod ext_cc2d;
mod ext_el_c2d;
pub use analytic_solvers::*;
pub use curve_curve::*;
pub use ext_cc2d::*;
pub use ext_el_c2d::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
