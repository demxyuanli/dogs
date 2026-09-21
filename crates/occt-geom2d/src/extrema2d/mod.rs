//! 2D extrema — point/curve minimum and maximum distances in the plane.
//!
//! Port of `Extrema_ExtPC2d`, `Extrema_ExtCC2d` (discretization + local
//! refinement). Source: `Extrema` (TKGeomBase).
mod prelude {

pub(crate) use std::cmp::Ordering;

pub(crate) use occt_core::gp::{
    GpAx22d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpVec2d,
};
pub(crate) use occt_core::precision::{ANGULAR, CONFUSION, RESOLUTION};

pub(crate) use crate::curve::Curve2d;

}


mod analytic_solvers;
mod curve_curve;
pub use analytic_solvers::*;
pub use curve_curve::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
