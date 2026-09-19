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
//! **UNPORTED (audit A1 / task T-67)** — the general point-surface path
//! (`point_surface_newton_all*`) is a **substitute, not a port**: it seeds a
//! 24×24 parameter grid (local extrema + global min/max + corners) and
//! Newton-refines with a **numeric** Jacobian. OCCT has no such grid: its
//! general path is `Extrema_GenExtPS::Perform` (`Extrema_GenExtPS.cxx`,
//! 1056 lines) — per-C2-interval sampling with `mySample`, then
//! `math_FunctionSetRoot` on `Extrema_GFuncExtPS`'s analytic system; the type
//! dispatch and iso-degenerate handling live in `Extrema_ExtPS`
//! (`Extrema_ExtPS.cxx`, 376 lines) with the analytic arms in
//! `Extrema_ExtPElS` (454 lines). Porting those three files is the faithful
//! fix; until then this module must not be cited as an `Extrema_ExtPS` port.
//!
//! `Surface::d2` **does** exist (`surface.rs:12`, analytic surfaces override
//! it; the trait default is a central difference, see A15/T-51), so the numeric
//! Jacobian here is not a trait limitation — it is part of the substitute.
//!
//! `dyn Surface` cannot be downcast, so analytic dispatch classifies by
//! geometric invariants (mirroring `brep_surface::classify_surface`): only
//! planes (constant normal) and spheres (equidistant samples from a solved
//! center) are classified; everything else goes through the Newton path.
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
