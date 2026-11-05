//! `HatchGen_PointOnElement`.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/HatchGen/HatchGen_PointOnElement.hxx:30-70`,
//! `.lxx:18-33`, `.cxx:25-292` (OCCT 8.0.0).

use occt_core::intres2d::{
    IntRes2dIntersectionPoint, IntRes2dPosition, IntRes2dSituation, IntRes2dTypeTrans,
};

use crate::abs::{Orientation, State};

use super::intersection_point::{IntersectionPoint, IntersectionType};

/// `HatchGen_PointOnElement` (`hxx:30-70`).
#[derive(Clone, Debug)]
pub struct PointOnElement {
    /// The inherited `HatchGen_IntersectionPoint` part.
    base: IntersectionPoint,
    /// `myType`
    my_type: IntersectionType,
}

impl Default for PointOnElement {
    fn default() -> Self {
        Self::new()
    }
}

impl PointOnElement {
    /// `HatchGen_PointOnElement()` (`cxx:27-30`).
    pub fn new() -> Self {
        Self {
            base: IntersectionPoint::new(),
            my_type: IntersectionType::Undetermined,
        }
    }

    /// `HatchGen_PointOnElement(const IntRes2d_IntersectionPoint&)` (`cxx:32-177`).
    ///
    /// `myPosit` comes from the second transition, the before/after states from
    /// the first one. The `IntRes2d_Touch` arm splits on `Situation()` and on
    /// `TransitionOfSecond().IsOpposite()`; `myIndex` stays 0 (`cxx:39`) and is
    /// set afterwards by the caller, and the two segment flags are reset at the
    /// end (`cxx:174-175`).
    pub fn from_intersection_point(point: &IntRes2dIntersectionPoint) -> Self {
        let trs_h = point.transition_of_first();
        let trs_e = point.transition_of_second();
        let mut p = Self::new();

        p.base.set_parameter(point.param_on_second());

        match trs_e.position_on_curve() {
            IntRes2dPosition::Head => p.base.set_position(Orientation::Forward), // `cxx:45-47`
            IntRes2dPosition::Middle => p.base.set_position(Orientation::Internal), // `cxx:48-50`
            IntRes2dPosition::End => p.base.set_position(Orientation::Reversed), // `cxx:51-53`
        }

        match trs_h.transition_type() {
            IntRes2dTypeTrans::In => {
                // `cxx:58-63`
                p.base.set_state_before(State::Out);
                p.base.set_state_after(State::In);
                p.my_type = if p.base.position() == Orientation::Internal {
                    IntersectionType::True
                } else {
                    IntersectionType::Touch
                };
            }
            IntRes2dTypeTrans::Out => {
                // `cxx:64-69`
                p.base.set_state_before(State::In);
                p.base.set_state_after(State::Out);
                p.my_type = if p.base.position() == Orientation::Internal {
                    IntersectionType::True
                } else {
                    IntersectionType::Touch
                };
            }
            IntRes2dTypeTrans::Touch => {
                // `cxx:71-162`
                match trs_h.situation() {
                    IntRes2dSituation::Inside => {
                        p.my_type = IntersectionType::Tangent;
                        match p.base.position() {
                            Orientation::Forward => {
                                if trs_e.is_opposite() {
                                    p.base.set_state_before(State::In);
                                    p.base.set_state_after(State::Out);
                                } else {
                                    p.base.set_state_before(State::Out);
                                    p.base.set_state_after(State::In);
                                }
                            }
                            Orientation::Internal => {
                                p.base.set_state_before(State::In);
                                p.base.set_state_after(State::In);
                            }
                            Orientation::Reversed => {
                                if trs_e.is_opposite() {
                                    p.base.set_state_before(State::Out);
                                    p.base.set_state_after(State::In);
                                } else {
                                    p.base.set_state_before(State::In);
                                    p.base.set_state_after(State::Out);
                                }
                            }
                            Orientation::External => {}
                        }
                    }
                    IntRes2dSituation::Outside => {
                        p.my_type = IntersectionType::Tangent;
                        match p.base.position() {
                            Orientation::Forward => {
                                if trs_e.is_opposite() {
                                    p.base.set_state_before(State::Out);
                                    p.base.set_state_after(State::In);
                                } else {
                                    p.base.set_state_before(State::In);
                                    p.base.set_state_after(State::Out);
                                }
                            }
                            Orientation::Internal => {
                                p.base.set_state_before(State::Out);
                                p.base.set_state_after(State::Out);
                            }
                            Orientation::Reversed => {
                                if trs_e.is_opposite() {
                                    p.base.set_state_before(State::In);
                                    p.base.set_state_after(State::Out);
                                } else {
                                    p.base.set_state_before(State::Out);
                                    p.base.set_state_after(State::In);
                                }
                            }
                            Orientation::External => {}
                        }
                    }
                    IntRes2dSituation::Unknown => {
                        p.base.set_state_before(State::Unknown);
                        p.base.set_state_after(State::Unknown);
                        p.my_type = IntersectionType::Tangent;
                    }
                }
            }
            IntRes2dTypeTrans::Undecided => {
                // `cxx:164-169`
                p.base.set_state_before(State::Unknown);
                p.base.set_state_after(State::Unknown);
                p.my_type = IntersectionType::Undetermined;
            }
        }

        p.base.set_segment_beginning(false);
        p.base.set_segment_end(false);
        p
    }

    /// `SetIntersectionType(Type)` (`.lxx:20-23`).
    pub fn set_intersection_type(&mut self, the_type: IntersectionType) {
        self.my_type = the_type;
    }

    /// `IntersectionType()` (`.lxx:29-32`).
    pub fn intersection_type(&self) -> IntersectionType {
        self.my_type
    }

    /// `IsIdentical(Point, Confusion)` (`cxx:181-188`).
    pub fn is_identical(&self, point: &PointOnElement, confusion: f64) -> bool {
        let delta = (self.base.parameter() - point.base.parameter()).abs();
        delta <= confusion
            && self.base.index() == point.base.index()
            && self.base.position() == point.base.position()
            && self.my_type == point.my_type
            && self.base.state_before() == point.base.state_before()
            && self.base.state_after() == point.base.state_after()
            && self.base.segment_beginning() == point.base.segment_beginning()
            && self.base.segment_end() == point.base.segment_end()
    }

    /// `IsDifferent(Point, Confusion)` (`cxx:194-202`).
    pub fn is_different(&self, point: &PointOnElement, confusion: f64) -> bool {
        !self.is_identical(point, confusion)
    }

    /// `SetIndex(Index)` (`cxx:38-40`).
    pub fn set_index(&mut self, index: i32) {
        self.base.set_index(index);
    }

    /// `Index()` (`cxx:47-49`).
    pub fn index(&self) -> i32 {
        self.base.index()
    }

    /// `SetParameter(Parameter)` (`cxx:57-59`).
    pub fn set_parameter(&mut self, parameter: f64) {
        self.base.set_parameter(parameter);
    }

    /// `Parameter()` (`cxx:66-68`).
    pub fn parameter(&self) -> f64 {
        self.base.parameter()
    }

    /// `SetPosition(Position)` (`cxx:76-78`).
    pub fn set_position(&mut self, position: Orientation) {
        self.base.set_position(position);
    }

    /// `Position()` (`cxx:85-87`).
    pub fn position(&self) -> Orientation {
        self.base.position()
    }

    /// `SetStateBefore(State)` (`cxx:95-97`).
    pub fn set_state_before(&mut self, state: State) {
        self.base.set_state_before(state);
    }

    /// `StateBefore()` (`cxx:104-106`).
    pub fn state_before(&self) -> State {
        self.base.state_before()
    }

    /// `SetStateAfter(State)` (`cxx:115-117`).
    pub fn set_state_after(&mut self, state: State) {
        self.base.set_state_after(state);
    }

    /// `StateAfter()` (`cxx:124-126`).
    pub fn state_after(&self) -> State {
        self.base.state_after()
    }

    /// `SetSegmentBeginning(State = true)` (`cxx:134-136`).
    pub fn set_segment_beginning(&mut self, state: bool) {
        self.base.set_segment_beginning(state);
    }

    /// `SegmentBeginning()` (`cxx:144-152`).
    pub fn segment_beginning(&self) -> bool {
        self.base.segment_beginning()
    }

    /// `SetSegmentEnd(State = true)` (`cxx:159-161`).
    pub fn set_segment_end(&mut self, state: bool) {
        self.base.set_segment_end(state);
    }

    /// `SegmentEnd()` (`cxx:168-170`).
    pub fn segment_end(&self) -> bool {
        self.base.segment_end()
    }
}
