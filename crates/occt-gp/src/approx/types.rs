//! `Approx_Status` and `Approx_ParametrizationType` — the two plain enums the
//! `Approx` / `AppDef` / `AppBlend` drivers pass around.
//!
//! Source: `TKGeomBase/Approx/Approx_Status.hxx:20-26` and
//! `TKGeomBase/Approx/Approx_ParametrizationType.hxx:20-29` (OCCT 8.0.0).
//!
//! Both are unscoped C++ enums; the ports are Rust enums with the same
//! variants and the same declaration order.

/// `Approx_Status` (`Approx_Status.hxx:20-26`).
///
/// Value order is kept: `Approx_PointsAdded`, `Approx_NoPointsAdded`,
/// `Approx_NoApproximation`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproxStatus {
    /// `Approx_PointsAdded`.
    PointsAdded,
    /// `Approx_NoPointsAdded`.
    NoPointsAdded,
    /// `Approx_NoApproximation`.
    NoApproximation,
}

/// `Approx_ParametrizationType` (`Approx_ParametrizationType.hxx:20-29`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproxParametrizationType {
    /// `Approx_ChordLength`: parameters proportional to point distances.
    ChordLength,
    /// `Approx_Centripetal`: parameters proportional to the square roots of
    /// point distances.
    Centripetal,
    /// `Approx_IsoParametric`: parameters distributed uniformly.
    IsoParametric,
}

impl Default for ApproxParametrizationType {
    /// `AppBlend_AppSurf` initialises `paramtype` to `Approx_ChordLength`
    /// (`AppBlend_AppSurf.gxx:58`); the same default is used here.
    fn default() -> Self {
        Self::ChordLength
    }
}
