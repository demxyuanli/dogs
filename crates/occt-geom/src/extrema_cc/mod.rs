//! Curve–curve extrema (3D).
//!
//! Port of `Extrema_ExtCC`, `Extrema_ExtElC`, `Extrema_ECC` and
//! `Extrema_LocateExtCC` (TKGeomBase). Analytic exact solvers for the cheap
//! elementary pairs (skew-line distance, line–circle, coaxial/parallel
//! circle–circle) plus a Newton refinement of F1(u,v) = (C1−C2)·C1′/|C1′| = 0,
//! F2(u,v) = (C1−C2)·C2′/|C2′| = 0 (`poly_roots::f_and_jac` / `refine_curve_curve`).
//!
//! **UNPORTED (audit A7 / task T-66)** — the **seed set** for the general path
//! is invented: `curve_curve` builds a uniform `n1 × n2` grid, keeps local min/max of
//! the sampled squared distance plus the best point of each grid boundary, then
//! Newton-refines each seed. OCCT `Extrema_GenExtCC::Perform` instead runs
//! `math_GlobOptMin` over the 2D parameter box (and `Extrema_ECC` combines the
//! analytic `Extrema_ExtElC` arms with `math_FunctionSetRoot` from defined
//! starts). The in-repo `occt-math/globoptmin.rs` already ports
//! `math_GlobOptMin`, so the faithful replacement is available; until then this
//! path must not be treated as aligned.
//!
//! `dyn Curve` cannot be downcast, so analytic dispatch classifies by geometric
//! invariants (mirroring `extrema_pc`): only lines (constant tangent) and
//! circles (coplanar + equidistant samples) are classified; anything else goes
//! through the seed + Newton path.
mod prelude {

pub(crate) use std::cmp::Ordering;

pub(crate) use occt_core::elib::clib;
pub(crate) use occt_core::gp::{GpCirc, GpLin, GpPnt, GpVec};
pub(crate) use occt_core::precision::{ANGULAR, CONFUSION, RESOLUTION};

pub(crate) use crate::curve::Curve;
pub(crate) use crate::extrema::ExtremaPair;

}


mod poly_roots;
mod curve_curve;
mod glob_opt_func;
mod general_extrema;
pub use poly_roots::*;
pub use curve_curve::*;
pub use glob_opt_func::*;
pub use general_extrema::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
