//! `TopOpeBRepDS_ListOfShapeOn1State` — a list of shapes with a split bit.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_ListOfShapeOn1State.hxx:26-50`,
//! `TopOpeBRepDS_ListOfShapeOn1State.cxx:22-66`.

use crate::shape::TopoShape;

/// Source: `TopOpeBRepDS_ListOfShapeOn1State.hxx:46-48`
#[derive(Debug, Clone, Default)]
pub struct ListOfShapeOn1State {
    /// `myList`
    list: Vec<TopoShape>,
    /// `mySplits`
    splits: i32,
}

impl ListOfShapeOn1State {
    /// `TopOpeBRepDS_ListOfShapeOn1State()` (`cxx:22-25`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `ListOnState()` (`cxx:29-32`).
    pub fn list_on_state(&self) -> &[TopoShape] {
        &self.list
    }

    /// `ChangeListOnState()` (`cxx:36-39`).
    pub fn change_list_on_state(&mut self) -> &mut Vec<TopoShape> {
        &mut self.list
    }

    /// `Clear()` (`cxx:42-46`): empties the list and resets the split bits.
    pub fn clear(&mut self) {
        self.list.clear();
        self.splits = 0;
    }

    /// `IsSplit()` (`cxx:50-57`): bit 0 of `mySplits`.
    pub fn is_split(&self) -> bool {
        self.splits & 1 != 0
    }

    /// `Split(B = true)` (`cxx:61-72`): sets or clears bit 0 of `mySplits`.
    pub fn split(&mut self, b: bool) {
        let mask = 1;
        if b {
            self.splits |= mask;
        } else {
            self.splits &= !mask;
        }
    }
}
