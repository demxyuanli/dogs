//! `HatchGen_IntersectionPoint` and its two enums.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/HatchGen/HatchGen_IntersectionPoint.hxx:30-94`,
//! `HatchGen_IntersectionPoint.cxx:20-172`, `HatchGen_IntersectionType.hxx:20-28`,
//! `HatchGen_ErrorStatus.hxx:21-28` (OCCT 8.0.0).
//!
//! `HatchGen_IntersectionPoint` is the abstract base of
//! `HatchGen_PointOnElement` and `HatchGen_PointOnHatching`; OCCT reaches the
//! fields through inheritance, the port embeds one [`IntersectionPoint`] per
//! derived type and forwards the accessors.

use crate::abs::{Orientation, State};

/// `RealLast()` (`Standard_Real.hxx:179-182`).
const REAL_LAST: f64 = f64::MAX;

/// `HatchGen_IntersectionType` (`hxx:20-28`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntersectionType {
    /// `HatchGen_TRUE`
    True,
    /// `HatchGen_TOUCH`
    Touch,
    /// `HatchGen_TANGENT`
    Tangent,
    /// `HatchGen_UNDETERMINED`
    Undetermined,
}

/// `HatchGen_ErrorStatus` (`hxx:21-28`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ErrorStatus {
    /// `HatchGen_NoProblem`
    NoProblem,
    /// `HatchGen_TrimFailure`
    TrimFailure,
    /// `HatchGen_TransitionFailure`
    TransitionFailure,
    /// `HatchGen_IncoherentParity`
    IncoherentParity,
    /// `HatchGen_IncompatibleStates`
    IncompatibleStates,
}

/// `HatchGen_IntersectionPoint` (`hxx:30-94`).
#[derive(Clone, Debug)]
pub struct IntersectionPoint {
    /// `myIndex`
    my_index: i32,
    /// `myParam`
    my_param: f64,
    /// `myPosit`
    my_posit: Orientation,
    /// `myBefore`
    my_before: State,
    /// `myAfter`
    my_after: State,
    /// `mySegBeg`
    my_seg_beg: bool,
    /// `mySegEnd`
    my_seg_end: bool,
}

impl Default for IntersectionPoint {
    fn default() -> Self {
        Self::new()
    }
}

impl IntersectionPoint {
    /// `HatchGen_IntersectionPoint()` (`cxx:22-31`).
    ///
    /// OCCT makes the constructor `protected` (only the two derived classes
    /// build it); the port keeps it public because the embedding structs need
    /// it.
    pub fn new() -> Self {
        Self {
            my_index: 0,
            my_param: REAL_LAST,
            my_posit: Orientation::Internal,
            my_before: State::Unknown,
            my_after: State::Unknown,
            my_seg_beg: false,
            my_seg_end: false,
        }
    }

    /// `SetIndex(Index)` (`cxx:38-40`).
    pub fn set_index(&mut self, index: i32) {
        self.my_index = index;
    }

    /// `Index()` (`cxx:47-49`).
    pub fn index(&self) -> i32 {
        self.my_index
    }

    /// `SetParameter(Parameter)` (`cxx:57-59`).
    pub fn set_parameter(&mut self, parameter: f64) {
        self.my_param = parameter;
    }

    /// `Parameter()` (`cxx:66-68`).
    pub fn parameter(&self) -> f64 {
        self.my_param
    }

    /// `SetPosition(Position)` (`cxx:76-78`).
    pub fn set_position(&mut self, position: Orientation) {
        self.my_posit = position;
    }

    /// `Position()` (`cxx:85-87`).
    pub fn position(&self) -> Orientation {
        self.my_posit
    }

    /// `SetStateBefore(State)` (`cxx:95-97`).
    pub fn set_state_before(&mut self, state: State) {
        self.my_before = state;
    }

    /// `StateBefore()` (`cxx:104-106`).
    pub fn state_before(&self) -> State {
        self.my_before
    }

    /// `SetStateAfter(State)` (`cxx:115-117`).
    pub fn set_state_after(&mut self, state: State) {
        self.my_after = state;
    }

    /// `StateAfter()` (`cxx:124-126`).
    pub fn state_after(&self) -> State {
        self.my_after
    }

    /// `SetSegmentBeginning(State = true)` (`cxx:134-136`).
    pub fn set_segment_beginning(&mut self, state: bool) {
        self.my_seg_beg = state;
    }

    /// `SegmentBeginning()` (`cxx:144-152`).
    pub fn segment_beginning(&self) -> bool {
        self.my_seg_beg
    }

    /// `SetSegmentEnd(State = true)` (`cxx:159-161`).
    pub fn set_segment_end(&mut self, state: bool) {
        self.my_seg_end = state;
    }

    /// `SegmentEnd()` (`cxx:168-170`).
    pub fn segment_end(&self) -> bool {
        self.my_seg_end
    }
}
