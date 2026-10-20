//! `BRepMesh_EdgeDeflection` — absolute deflection and consistency checks.

// ---------------------------------------------------------------------------
// EdgeDeflection
// ---------------------------------------------------------------------------

use super::prelude::*;
use super::*;

/// EdgeDeflection computation helpers (`BRepMesh_Deflection`).
pub struct EdgeDeflection;

impl EdgeDeflection {
    /// Absolute deflection for a shape of `max_shape_size` given a relative
    /// deflection. Mirrors `ComputeAbsoluteDeflection` (clamped 0.5×..2×).
    pub fn compute_absolute_deflection(relative: f64, max_shape_size: f64) -> f64 {
        if relative <= 0.0 {
            return 0.0;
        }
        if max_shape_size <= 0.0 {
            return relative;
        }
        let shape_size = max_shape_size.max(relative);
        let coeff = (max_shape_size / (2.0 * shape_size)).clamp(0.5, 2.0);
        coeff * shape_size * relative
    }

    /// Whether an existing polygon deflection fits the required one.
    /// Mirrors `IsConsistent`: current must not exceed required by `ratio`, and
    /// (when quality decrease is allowed) must not undercut it either.
    pub fn is_consistent(
        current: f64,
        required: f64,
        allow_quality_decrease: bool,
        ratio: f64,
    ) -> bool {
        current < (1.0 + ratio) * required
            && (!allow_quality_decrease || current > (1.0 - ratio) * required)
    }
}
