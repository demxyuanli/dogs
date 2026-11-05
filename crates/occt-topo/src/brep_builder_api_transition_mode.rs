//! `BRepBuilderAPI_TransitionMode` — the corner mode of a pipe shell.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepBuilderAPI/BRepBuilderAPI_TransitionMode.hxx:21-27`
//! (OCCT 8.0.0). A plain enumerator with no `.cxx`.
//!
//! `BRepOffsetAPI_MakePipeShell::SetTransitionMode`
//! (`BRepOffsetAPI_MakePipeShell.cxx:198-201`) forwards it with a C cast to
//! `BRepFill_TransitionStyle`, so the two enums line up position by position.
//! [`TransitionMode::to_transition_style`] is that cast.

use crate::brep_fill_sweep::TransitionStyle;

/// `BRepBuilderAPI_TransitionMode` (`BRepBuilderAPI_TransitionMode.hxx:21-27`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitionMode {
    /// `BRepBuilderAPI_Transformed`.
    Transformed,
    /// `BRepBuilderAPI_RightCorner`.
    RightCorner,
    /// `BRepBuilderAPI_RoundCorner`.
    RoundCorner,
}

impl TransitionMode {
    /// The `(BRepFill_TransitionStyle)Mode` cast of
    /// `BRepOffsetAPI_MakePipeShell.cxx:199`. `BRepFill_TransitionStyle`
    /// (`BRepFill_TransitionStyle.hxx:20-25`) holds only `Modified`, `Right`
    /// and `Round`, so the identifiers agree index by index; in particular the
    /// API's `Transformed` — the documented default — really selects
    /// `BRepFill_Modified`, exactly as the OCCT cast does.
    pub fn to_transition_style(self) -> TransitionStyle {
        match self {
            TransitionMode::Transformed => TransitionStyle::Modified,
            TransitionMode::RightCorner => TransitionStyle::Right,
            TransitionMode::RoundCorner => TransitionStyle::Round,
        }
    }
}
