//! `AppBlend_Debug` — the two process-wide switches `AppBlend_AppSurf` reads
//! while it builds the approximation input.
//!
//! Source: `TKGeomAlgo/AppBlend/AppBlend_Debug.cxx:22-45` (OCCT 8.0.0).
//!
//! Translation notes:
//! - OCCT's file-scope `static bool` become [`AtomicBool`], because the port is
//!   used from several threads; both start at `false`, exactly like `cxx:22` /
//!   `cxx:36`.
//! - `AppBlend_GetContextSplineApprox` has no caller inside
//!   `AppBlend_AppSurf.gxx`; it is kept because the OCCT API exposes it.

use std::sync::atomic::{AtomicBool, Ordering};

/// `static bool AppBlend_ContextSplineApprox` (`cxx:22`).
static CONTEXT_SPLINE_APPROX: AtomicBool = AtomicBool::new(false);

/// `static bool AppBlend_ContextApproxWithNoTgt` (`cxx:36`).
static CONTEXT_APPROX_WITH_NO_TGT: AtomicBool = AtomicBool::new(false);

/// `AppBlend_SetContextSplineApprox(const bool)` (`cxx:24-27`).
pub fn app_blend_set_context_spline_approx(b: bool) {
    CONTEXT_SPLINE_APPROX.store(b, Ordering::Relaxed);
}

/// `AppBlend_GetContextSplineApprox()` (`cxx:29-32`).
pub fn app_blend_get_context_spline_approx() -> bool {
    CONTEXT_SPLINE_APPROX.load(Ordering::Relaxed)
}

/// `AppBlend_SetContextApproxWithNoTgt(const bool)` (`cxx:38-41`).
pub fn app_blend_set_context_approx_with_no_tgt(b: bool) {
    CONTEXT_APPROX_WITH_NO_TGT.store(b, Ordering::Relaxed);
}

/// `AppBlend_GetContextApproxWithNoTgt()` (`cxx:43-46`).
pub fn app_blend_get_context_approx_with_no_tgt() -> bool {
    CONTEXT_APPROX_WITH_NO_TGT.load(Ordering::Relaxed)
}
