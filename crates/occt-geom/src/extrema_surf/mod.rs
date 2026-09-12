//! Point-surface / curve-surface extrema. Port of `Extrema_ExtPS`,
//! `Extrema_ExtPElS`, `Extrema_ExtPRevS`, `Extrema_ExtPExtS`,
//! `Extrema_ExtCS`, `Extrema_ExtElCS`, `Extrema_FuncExtCS` (TKGeomBase).
//!
//! The analytic point-to-surface solvers (`Extrema_ExtPElS`: plane, sphere,
//! cylinder, cone, torus) are ported exactly from the OCCT `.cxx`. Revolved
//! and extruded surfaces (`ExtPRevS` / `ExtPExtS`) are heavy (they reduce to a
//! point-curve extrema on the generating curve) and are routed through the
//! general Newton path, exactly as OCCT does for the non-analytically-
//! computable cases.
//!
//! `dyn Surface` cannot be downcast, so analytic dispatch classifies by
//! geometric invariants (mirroring `brep_surface::classify_surface`): only
//! planes (constant normal) and spheres (equidistant samples from a solved
//! center) are classified; everything else goes through the Newton path.
//!
//! The Newton systems are solved with a NUMERIC Jacobian: the `Surface` trait
//! exposes only `d0`/`d1` (no second derivatives), so the Jacobian of the
//! orthogonality conditions `F = ((S-P)·Su, (S-P)·Sv)` is computed by central
//! finite differences of `d1` with `eps ≈ 1e-6` relative to the parameter
//! range. // ponytail: numeric Jacobian, Surface trait lacks d2
mod prelude {

pub(crate) use std::cmp::Ordering;

pub(crate) use occt_core::elib::{clib, slib};
pub(crate) use occt_core::gp::{GpAx3, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpVec};
pub(crate) use occt_core::precision::CONFUSION;

pub(crate) use crate::curve::Curve;
pub(crate) use crate::extrema::ExtremaPair;
pub(crate) use crate::surface::Surface;

}


mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
