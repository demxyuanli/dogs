//! Curve–curve extrema (3D).
//!
//! Port of `Extrema_ExtCC`, `Extrema_ExtElC`, `Extrema_ECC` and
//! `Extrema_LocateExtCC` (TKGeomBase). Analytic exact solvers for the cheap
//! elementary pairs (skew-line distance, line–circle, coaxial/parallel
//! circle–circle) plus a Newton refinement of F1(u,v) = (C1−C2)·C1′/|C1′| = 0,
//! F2(u,v) = (C1−C2)·C2′/|C2′| = 0 (`p01::f_and_jac` / `refine_curve_curve`).
//!
//! **UNPORTED (audit A7 / task T-66)** — the **seed set** for the general path
//! is invented: `p02` builds a uniform `n1 × n2` grid, keeps local min/max of
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
pub(crate) use occt_core::gp::{GpAx2, GpCirc, GpDir, GpLin, GpPnt, GpVec};
pub(crate) use occt_core::precision::{ANGULAR, CONFUSION, RESOLUTION};

pub(crate) use crate::curve::Curve;
pub(crate) use crate::extrema::ExtremaPair;

}


mod p01;
mod p02;
mod p03;
pub use p01::*;
pub use p02::*;
pub use p03::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
