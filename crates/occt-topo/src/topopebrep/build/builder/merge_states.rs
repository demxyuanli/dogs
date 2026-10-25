//! The `Merge.cxx` state predicates of `TopOpeBRepBuild_Builder`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Merge.cxx`
//! (OCCT 8.0.0): `Opec12` (`853-857`), `Opec21` (`861-865`), `Opecom`
//! (`869-873`), `Opefus` (`877-881`), `Classify` (`341-346`) and
//! `ChangeClassify` (`348-354`). They only read `myState1` / `myState2` and the
//! classify pair, which is why they live outside `builder/mod.rs`.

use crate::abs::State;
use crate::shape::TopoShape;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;

use super::Builder;

impl Builder {
    /// `Opec12()` (`.cxx:853-857`).
    pub fn opec12(&self) -> bool {
        self.my_state1 == State::Out && self.my_state2 == State::In
    }

    /// `Opec21()` (`.cxx:861-865`).
    pub fn opec21(&self) -> bool {
        self.my_state1 == State::In && self.my_state2 == State::Out
    }

    /// `Opecom()` (`.cxx:869-873`).
    pub fn opecom(&self) -> bool {
        self.my_state1 == State::In && self.my_state2 == State::In
    }

    /// `Opefus()` (`.cxx:877-881`).
    pub fn opefus(&self) -> bool {
        self.my_state1 == State::Out && self.my_state2 == State::Out
    }

    /// `Classify()` (`Merge.cxx:341-346`).
    pub fn classify(&self) -> bool {
        if !self.my_classify_def {
            return true;
        }
        self.my_classify_val
    }

    /// `ChangeClassify(classify)` (`Merge.cxx:348-354`).
    pub fn change_classify(&mut self, classify: bool) {
        self.my_classify_val = classify;
        self.my_classify_def = true;
    }
}

/// `TopExp::MapShapes(S, M)` (`TopExp.cxx:31-45`): `S` and every sub-shape,
/// `S` first, duplicates skipped through `Contains`.
pub(crate) fn map_shapes_into(s: &TopoShape, m: &mut IndexedShapeSet) {
    if m.find_index(s) != 0 {
        return;
    }
    m.add(s);
    let children = match s.tshape.read() {
        Ok(t) => t.children.clone(),
        Err(_) => return,
    };
    for c in children {
        map_shapes_into(&c, m);
    }
}
