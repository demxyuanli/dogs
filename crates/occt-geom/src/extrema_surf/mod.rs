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
//! **UNPORTED (audit A1 / task T-67)** — the general point-surface *engine*
//! (`point_surface_newton_all*`) is a **substitute, not a port**: it seeds a
//! 24×24 parameter grid (local extrema + global min/max + corners) and
//! Newton-refines with a **numeric** Jacobian. OCCT has no such grid: its
//! general path is `Extrema_GenExtPS::Perform` (`Extrema_GenExtPS.cxx:968`,
//! 1195 lines in 8.0.0) — per-C2-interval sampling with `GetGridPoints` /
//! `BuildGrid` over `GeomGridEval_Surface`, a `Bnd_Sphere` UBTree
//! (`BuildTree`), then `math_FunctionSetRoot` on `Extrema_FuncPSNorm`'s
//! analytic system. Porting that engine is the remaining faithful fix.
//!
//! The **layer above it is ported**: [`ExtPs`] (`p03.rs`) is `Extrema_ExtPS` —
//! the type dispatch to the analytic `Extrema_ExtPElS` arms, the ±1e10 window
//! clamp, the `nbU/nbV` sampling counts with the 300 sample `IsoIsDeg` rule,
//! `TreatSolution`'s periodic normalization and window test, and the
//! `IsDone`/`NbExt` semantics. Two of its three engines are still unported and
//! are recorded at their call sites: `Extrema_ExtPExtS`/`Extrema_ExtPRevS`
//! (extrusion/revolution, `Extrema_ExtPS.cxx:292-343`) and `Extrema_GenExtPS`
//! (`Extrema_ExtPS.cxx:346`).
//!
//! `Surface::d2` **does** exist (`surface.rs:12`, analytic surfaces override
//! it; the trait default is a central difference, see A15/T-51), so the numeric
//! Jacobian here is not a trait limitation — it is part of the substitute.
//!
//! `dyn Surface` cannot be downcast, so analytic dispatch classifies with the
//! `GetType()`-equivalent trait queries (`gp_pln` / `gp_cylinder` / `gp_cone` /
//! `gp_sphere` / `gp_torus` / `is_surface_of_linear_extrusion` /
//! `is_surface_of_revolution` / `is_bspline_surface` / `is_bezier_surface`,
//! see [`ext_ps_surface_type`]), not from a sampling classifier.
mod prelude {

pub(crate) use std::cmp::Ordering;

pub(crate) use occt_core::elib::{clib, slib};
pub(crate) use occt_core::gp::{GpAx3, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpVec};
pub(crate) use occt_core::precision::{CONFUSION, PCONFUSION};

pub(crate) use crate::curve::Curve;
pub(crate) use crate::extrema::ExtremaPair;
pub(crate) use crate::surface::Surface;

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
