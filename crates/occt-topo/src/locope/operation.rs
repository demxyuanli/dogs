//! `LocOpe_Operation` (`LocOpe_Operation.hxx:20-25`): the operation a
//! `LocOpe_Gluer` performs, decided while it binds the glued faces.
//!
//! `LocOpe_Gluer::Bind` starts in `LocOpe_INVALID` (`.lxx:26`) and the first
//! bound face decides between `FUSE` and `CUT` from the relative orientation of
//! the two faces; any later bind that disagrees puts it back to `INVALID`
//! (`LocOpe_Gluer.cxx:118-121`) and `Perform` then refuses to run
//! (`cxx:163-166`).

/// `enum LocOpe_Operation` (`LocOpe_Operation.hxx:20-25`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LocOpeOperation {
    /// `LocOpe_INVALID` (`hxx:24`, the initial value).
    #[default]
    Invalid,
    /// `LocOpe_FUSE` (`hxx:22`).
    Fuse,
    /// `LocOpe_CUT` (`hxx:23`).
    Cut,
}

impl LocOpeOperation {
    /// `LocOpe_INVALID` (`hxx:24`): no bind has decided the operation yet, or
    /// the binds disagreed.
    pub fn is_invalid(&self) -> bool {
        *self == LocOpeOperation::Invalid
    }
}
