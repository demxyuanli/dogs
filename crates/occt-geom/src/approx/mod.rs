//! `Approx_*` — the general approximation drivers of `TKGeomBase`.
//!
//! Ported: `Approx_MCurvesToBSpCurve`, `Approx_ComputeLine` (the
//! `AppDef_Compute` instantiation) together with its `Parameters`,
//! `Perform`, `CheckMultiCurve`, `Compute`, `ComputeCurve` and the four
//! tangency helpers, `Approx_BSplComputeLine` (the `AppDef_BSplineCompute`
//! instantiation) with its `Perform`, `Parameters`, `Compute`, `Interpol` and
//! tangency helpers, and `Approx_BSplineApproxInterp` (used by
//! `GeomFill_Gordon`).
//!
//! Not ported yet: `Approx_Curve3d`, `Approx_SweepApproximation`
//! (`Approx_SweepApproximation.cxx`).

pub mod bspline_approx_interp;
pub mod bspl_compute_line;
pub mod compute_line;
pub mod line_tool;
pub mod m_curves_to_bsp_curve;
