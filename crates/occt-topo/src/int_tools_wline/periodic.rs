//! Geometric U / V periods of a surface, used by the walking-line helpers.
//!
//! Source: sub-functions of `IntTools_WLineTool::DecompositionOfWLine`
//! (`ModelingAlgorithms/TKTopAlgo/IntTools/IntTools_WLineTool.cxx`).

use occt_geom::Surface;

/// Geometric U period, or `None` when the surface is not U-periodic.
pub fn u_period(s: &dyn Surface) -> Option<f64> {
    if !s.is_u_periodic() {
        return None;
    }
    let (a, b) = s.u_range();
    if a.is_finite() && b.is_finite() && (b - a).abs() > 1e-30 {
        Some((b - a).abs())
    } else {
        Some(2.0 * std::f64::consts::PI)
    }
}

/// Geometric V period, or `None` when the surface is not V-periodic.
pub fn v_period(s: &dyn Surface) -> Option<f64> {
    if !s.is_v_periodic() {
        return None;
    }
    let (a, b) = s.v_range();
    if a.is_finite() && b.is_finite() && (b - a).abs() > 1e-30 {
        Some((b - a).abs())
    } else {
        Some(2.0 * std::f64::consts::PI)
    }
}
