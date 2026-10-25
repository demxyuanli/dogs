//! `TopOpeBRepDS_ShapeWithState` — a shape split into its `IN` / `OUT` / `ON`
//! parts, plus the state of the whole shape.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_ShapeWithState.hxx:24-58`,
//! `TopOpeBRepDS_ShapeWithState.cxx:20-126`.

use crate::abs::State;
use crate::shape::TopoShape;

/// Source: `TopOpeBRepDS_ShapeWithState.hxx:51-56`
#[derive(Debug, Clone)]
pub struct ShapeWithState {
    /// `myPartIn`
    part_in: Vec<TopoShape>,
    /// `myPartOut`
    part_out: Vec<TopoShape>,
    /// `myPartOn`
    part_on: Vec<TopoShape>,
    /// `myState`
    state: State,
    /// `myIsSplitted`
    is_splitted: bool,
}

impl Default for ShapeWithState {
    /// `TopOpeBRepDS_ShapeWithState()` (`cxx:22-26`).
    fn default() -> Self {
        Self {
            part_in: Vec::new(),
            part_out: Vec::new(),
            part_on: Vec::new(),
            state: State::Unknown,
            is_splitted: false,
        }
    }
}

impl ShapeWithState {
    /// `TopOpeBRepDS_ShapeWithState()` (`cxx:22-26`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Part(aState)` (`cxx:30-45`): an empty list for `UNKNOWN`.
    pub fn part(&self, state: State) -> &[TopoShape] {
        match state {
            State::In => &self.part_in,
            State::Out => &self.part_out,
            State::On => &self.part_on,
            State::Unknown => &[],
        }
    }

    /// `AddPart(aShape, aState)` (`cxx:49-62`): `UNKNOWN` appends nothing.
    pub fn add_part(&mut self, shape: TopoShape, state: State) {
        match state {
            State::In => self.part_in.push(shape),
            State::Out => self.part_out.push(shape),
            State::On => self.part_on.push(shape),
            State::Unknown => {}
        }
    }

    /// `AddParts(aListOfShape, aState)` (`cxx:66-93`).
    pub fn add_parts(&mut self, shapes: &[TopoShape], state: State) {
        let target = match state {
            State::In => &mut self.part_in,
            State::Out => &mut self.part_out,
            State::On => &mut self.part_on,
            State::Unknown => return,
        };
        for s in shapes {
            target.push(s.clone());
        }
    }

    /// `SetState(aState)` (`cxx:97-100`).
    pub fn set_state(&mut self, state: State) {
        self.state = state;
    }

    /// `State()` (`cxx:104-107`).
    pub fn state(&self) -> State {
        self.state
    }

    /// `SetIsSplitted(anIsSplitted)` (`cxx:111-114`).
    pub fn set_is_splitted(&mut self, is_splitted: bool) {
        self.is_splitted = is_splitted;
    }

    /// `IsSplitted()` (`cxx:118-121`).
    pub fn is_splitted(&self) -> bool {
        self.is_splitted
    }
}
