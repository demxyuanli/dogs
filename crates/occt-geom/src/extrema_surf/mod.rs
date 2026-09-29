//! Point-surface / curve-surface extrema. Port of `Extrema_ExtPS`,
//! `Extrema_ExtPElS`, `Extrema_ExtPRevS`, `Extrema_ExtPExtS`,
//! `Extrema_ExtCS`, `Extrema_ExtElCS`, `Extrema_FuncExtCS` (TKGeomBase).
//!
//! The analytic point-to-surface solvers (`Extrema_ExtPElS`: plane, sphere,
//! cylinder, cone, torus) are ported exactly from the OCCT `.cxx`. Revolved
//! and extruded surfaces (`ExtPRevS` / `ExtPExtS`) are **ported** too
//! (`revolution_point_extrema.rs` / `extrusion_point_extrema.rs`): they reduce
//! to a point-curve extrema on the generating curve and fall back to the
//! general `Extrema_GenExtPS` engine exactly where OCCT does
//! (`IsCaseAnalyticallyComputable` false).
//!
//! The general point-surface engine **is ported**: [`gen_ext_ps`] is
//! `Extrema_GenExtPS` (8.0.0, 1195 lines): `Initialize`, the per-geometry
//! `GetGridPoints` / `fillParams` sampling, `BuildGrid` (node / UIsoEdge /
//! VIsoEdge / Face square distances), `BuildTree`'s B-spline sample raise and
//! `CorrectNbSamples`, `FindSolution` over `math_FunctionSetRoot` with the
//! analytic `Extrema_FuncPSNorm` system, and the MIN / MAX scans of
//! `Perform`. See `gen_ext_ps.rs` for the exact `cxx` lines.
//!
//! The **layer above it is ported**: [`ExtPs`] (`point_surface_extrema.rs`) is `Extrema_ExtPS` —
//! the type dispatch to the analytic `Extrema_ExtPElS` arms, the ±1e10 window
//! clamp, the `nbU/nbV` sampling counts with the 300 sample `IsoIsDeg` rule,
//! `TreatSolution`'s periodic normalization and window test, and the
//! `IsDone`/`NbExt` semantics. Its two reduced engines
//! `Extrema_ExtPExtS` / `Extrema_ExtPRevS` (extrusion/revolution,
//! `Extrema_ExtPS.cxx:292-343`) are ported and dispatched to; the shared
//! `Extrema_ExtPElC` reduction lives in `elementary_curve_extrema.rs`.
//!
//! **UNPORTED substitute (audit A15, closed by T-67)** — `point_surface_newton_all*`
//! (`numeric_extrema.rs`) is a **substitute, not a port** (24x24 grid seeding plus a
//! numeric-Jacobian Newton). It is no longer on any library path: the `Extrema_ExtPS`
//! dispatch runs [`gen_ext_ps`], and `point_surface_extrema_box` now reports
//! "no solution" on an empty window so its caller takes the OCCT `UVFromIso` branch
//! (`ShapeAnalysis_Surface.cxx:1449-1459`, ported as `pcurve_full::uv_from_iso`).
//! It survives only for the `newton_path_bspline_paraboloid_min` regression test.
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
pub(crate) use occt_core::gp::{GpAx1, GpAx2, GpAx3, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpTrsf, GpVec};
pub(crate) use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION, RESOLUTION, SQUARE_CONFUSION};

pub(crate) use crate::curve::Curve;
pub(crate) use crate::extrema::ExtremaPair;
pub(crate) use crate::surface::Surface;

}


mod analytic_solvers;
mod elementary_curve_extrema;
mod extrusion_point_extrema;
mod gen_ext_ps;
mod numeric_extrema;
mod point_surface_extrema;
mod revolution_point_extrema;
pub use analytic_solvers::*;
pub use extrusion_point_extrema::*;
pub use gen_ext_ps::*;
pub use numeric_extrema::*;
pub use point_surface_extrema::*;
pub use revolution_point_extrema::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
