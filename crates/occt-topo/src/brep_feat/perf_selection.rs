//! `BRepFeat_PerfSelection` — the selection semantics of a local `Perform`.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_PerfSelection.hxx:30-36`
//! (the whole file is the enum).

/// `enum BRepFeat_PerfSelection` (`hxx:30-36`), same order and names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerfSelection {
    /// `BRepFeat_NoSelection`
    NoSelection,
    /// `BRepFeat_SelectionFU` — a face up to which the operation runs.
    SelectionFU,
    /// `BRepFeat_SelectionU` — a point up to which the operation runs.
    SelectionU,
    /// `BRepFeat_SelectionSh` — a shape the operation applies on.
    SelectionSh,
    /// `BRepFeat_SelectionShU` — a shape up to which the operation runs.
    SelectionShU,
}

impl Default for PerfSelection {
    /// The `BRepFeat_Form` constructor's initial value (`lxx:22-35`).
    fn default() -> Self {
        Self::NoSelection
    }
}
