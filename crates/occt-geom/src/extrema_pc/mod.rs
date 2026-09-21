//! Exact point–curve extrema.
//!
//! Replaces the sampling-based point–curve path with the exact analytic
//! solvers of `Extrema_ExtPElC` (line, circle, ellipse, hyperbola, parabola)
//! plus the Newton-refined general path of `Extrema_GGExtPC` /
//! `Extrema_GFuncExtPC` (bracket sign changes of (C-P)·C′ = 0, refine each
//! root with d1/d2). Source: `Extrema` (TKGeomBase).
//!
//! `dyn Curve` cannot be downcast, so analytic dispatch classifies by
//! geometric invariants (mirroring `brep_surface::classify_surface`): only
//! lines (constant tangent) and circles (coplanar + equidistant samples) are
//! classified; anything else goes through the general Newton path.
mod prelude {

pub(crate) use std::cmp::Ordering;

pub(crate) use occt_core::elib::clib;
pub(crate) use occt_core::gp::{GpCirc, GpDir, GpElips, GpHypr, GpLin, GpParab, GpPnt, GpVec};
pub(crate) use occt_core::precision::{ANGULAR, CONFUSION, RESOLUTION};

pub(crate) use crate::curve::Curve;
pub(crate) use crate::extrema::ExtremaPair;

}


mod poly_roots;
mod point_curve;
mod general_extrema_pc;
pub use poly_roots::*;
pub use point_curve::*;
pub use general_extrema_pc::{extrema_ext_pc_min_in_range, extrema_ext_pc_range, ExtPcSolution};

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
