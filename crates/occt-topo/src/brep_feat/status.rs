//! `BRepFeat_Status` — the placement result reported by the hole feature.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Status.hxx:21-27`
//! (the whole file is the enum). This is not `BRepFeat_StatusError`
//! ([`super::status_error::StatusError`]): that one is the `myStatusError`
//! code of `BRepFeat_RibSlot`, while this one is the `myStatus` field of
//! `BRepFeat_MakeCylindricalHole`.

/// `enum BRepFeat_Status` (`hxx:21-27`), same order and names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// `BRepFeat_NoError`
    NoError,
    /// `BRepFeat_InvalidPlacement`
    InvalidPlacement,
    /// `BRepFeat_HoleTooLong`
    HoleTooLong,
}

impl Default for Status {
    /// The hole constructor's initial value (`.lxx:21-27`).
    fn default() -> Self {
        Self::NoError
    }
}
