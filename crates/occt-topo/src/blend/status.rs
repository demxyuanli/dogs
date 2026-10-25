//! `Blend_Status` and `Blend_DecrochStatus` — the outcome of one walking step
//! and the "has the solution left a restriction" flag.
//!
//! Source: `ModelingAlgorithms/TKFillet/Blend/Blend_Status.hxx:19-32`,
//! `Blend_DecrochStatus.hxx:19-28` (OCCT 8.0.0). Both are plain unscoped C++
//! enums without explicit values, so the declarations keep the OCCT order.

/// `enum Blend_Status` (`Blend_Status.hxx:19-32`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BlendStatus {
    /// `Blend_StepTooLarge`
    StepTooLarge,
    /// `Blend_StepTooSmall`
    StepTooSmall,
    /// `Blend_Backward`
    Backward,
    /// `Blend_SamePoints`
    SamePoints,
    /// `Blend_OnRst1`
    OnRst1,
    /// `Blend_OnRst2`
    OnRst2,
    /// `Blend_OnRst12`
    OnRst12,
    /// `Blend_OK`
    #[default]
    Ok,
}

/// `enum Blend_DecrochStatus` (`Blend_DecrochStatus.hxx:19-28`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BlendDecrochStatus {
    /// `Blend_NoDecroch`
    #[default]
    NoDecroch,
    /// `Blend_DecrochRst1`
    DecrochRst1,
    /// `Blend_DecrochRst2`
    DecrochRst2,
    /// `Blend_DecrochBoth`
    DecrochBoth,
}
