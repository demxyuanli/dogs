//! Curve–curve extrema (3D).
//!
//! Port of `Extrema_ExtCC`, `Extrema_ExtElC`, `Extrema_ECC` and
//! `Extrema_LocateExtCC` (TKGeomBase). Analytic exact solvers for the cheap
//! elementary pairs (skew-line distance, line–circle, coaxial/parallel
//! circle–circle) plus a general grid → bracket-sign-change → Newton path for
//! arbitrary curves (F1(u,v) = (C1−C2)·C1′/|C1′| = 0,
//! F2(u,v) = (C1−C2)·C2′/|C2′| = 0).
//!
//! `dyn Curve` cannot be downcast, so analytic dispatch classifies by geometric
//! invariants (mirroring `extrema_pc`): only lines (constant tangent) and
//! circles (coplanar + equidistant samples) are classified; anything else goes
//! through the Newton path.
mod prelude {

pub(crate) use std::cmp::Ordering;

pub(crate) use occt_core::elib::clib;
pub(crate) use occt_core::gp::{GpAx2, GpCirc, GpDir, GpLin, GpPnt, GpVec};
pub(crate) use occt_core::precision::{ANGULAR, CONFUSION, RESOLUTION};

pub(crate) use crate::curve::Curve;
pub(crate) use crate::extrema::ExtremaPair;

}


mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
