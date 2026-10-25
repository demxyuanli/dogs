//! `GeomFill_Trihedron` — the trihedron modes a sweep can pick.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Trihedron.hxx:20-31` (OCCT 8.0.0).

/// `enum GeomFill_Trihedron` (`GeomFill_Trihedron.hxx:20-31`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomFillTrihedron {
    /// `GeomFill_IsCorrectedFrenet`.
    CorrectedFrenet,
    /// `GeomFill_IsFixed`.
    Fixed,
    /// `GeomFill_IsFrenet`.
    Frenet,
    /// `GeomFill_IsConstantNormal`.
    ConstantNormal,
    /// `GeomFill_IsDarboux`.
    Darboux,
    /// `GeomFill_IsGuideAC`.
    GuideAC,
    /// `GeomFill_IsGuidePlan`.
    GuidePlan,
    /// `GeomFill_IsGuideACWithContact`.
    GuideACWithContact,
    /// `GeomFill_IsGuidePlanWithContact`.
    GuidePlanWithContact,
    /// `GeomFill_IsDiscreteTrihedron`.
    DiscreteTrihedron,
}
