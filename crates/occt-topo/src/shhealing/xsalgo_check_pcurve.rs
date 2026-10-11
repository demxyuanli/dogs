//! `XSAlgo_ShapeProcessor::CheckPCurve` — public entry used by the STEP
//! `TranslateEdgeLoop` advanced check.
//!
//! Full control flow lives in [`super::pcurve_ranges::check_pcurve`]
//! (`XSAlgo_ShapeProcessor.cxx:344-505`), including the `:403-505` deviation /
//! re-projection tail. Call sites that do not know seam status pass
//! `is_seam = false` (OCCT's `CheckPCurves` uses `sbwd->IsSeam(i)`).

use super::prelude::*;

/// Returns `true` when the edge keeps a pcurve on `face` after the check.
///
/// Mirrors `XSAlgo_ShapeProcessor::CheckPCurve` with `theIsSeam = false`.
pub fn xsalgo_check_pcurve(edge: &Edge, face: &Face, preci: f64) -> bool {
    xsalgo_check_pcurve_seam(edge, face, preci, false)
}

/// `XSAlgo_ShapeProcessor::CheckPCurve` with an explicit seam flag
/// (`XSAlgo_ShapeProcessor.cxx:344`, `theIsSeam`).
pub fn xsalgo_check_pcurve_seam(edge: &Edge, face: &Face, preci: f64, is_seam: bool) -> bool {
    super::pcurve_ranges::check_pcurve(edge, face, preci, is_seam)
}
