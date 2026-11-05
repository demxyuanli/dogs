//! `HatchGen_PointOnHatching`.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/HatchGen/HatchGen_PointOnHatching.hxx:33-88`,
//! `.cxx:25-237` (OCCT 8.0.0).

use occt_core::intres2d::{IntRes2dIntersectionPoint, IntRes2dPosition};

use crate::abs::{Orientation, State};

use super::intersection_point::IntersectionPoint;
use super::point_on_element::PointOnElement;

/// `HatchGen_PointOnHatching` (`hxx:33-88`).
#[derive(Clone, Debug)]
pub struct PointOnHatching {
    /// The inherited `HatchGen_IntersectionPoint` part.
    base: IntersectionPoint,
    /// `myPoints`
    my_points: Vec<PointOnElement>,
}

impl Default for PointOnHatching {
    fn default() -> Self {
        Self::new()
    }
}

impl PointOnHatching {
    /// `HatchGen_PointOnHatching()` (`cxx:27-29`, `= default`).
    pub fn new() -> Self {
        Self {
            base: IntersectionPoint::new(),
            my_points: Vec::new(),
        }
    }

    /// `HatchGen_PointOnHatching(const IntRes2d_IntersectionPoint&)` (`cxx:33-53`).
    ///
    /// `myParam` is the parameter on the *first* curve (the hatching), the
    /// position comes from the first transition. `myIndex` is 0 and both states
    /// are `UNKNOWN`, as in the base constructor.
    pub fn from_intersection_point(point: &IntRes2dIntersectionPoint) -> Self {
        let mut p = Self::new();
        p.base.set_parameter(point.param_on_first());
        match point.transition_of_first().position_on_curve() {
            IntRes2dPosition::Head => p.base.set_position(Orientation::Forward), // `cxx:38-40`
            IntRes2dPosition::Middle => p.base.set_position(Orientation::Internal), // `cxx:41-43`
            IntRes2dPosition::End => p.base.set_position(Orientation::Reversed), // `cxx:44-46`
        }
        p.base.set_state_before(State::Unknown);
        p.base.set_state_after(State::Unknown);
        p.base.set_segment_beginning(false);
        p.base.set_segment_end(false);
        p.my_points.clear();
        p
    }

    /// `AddPoint(Point, Confusion)` (`cxx:60-73`): the point is appended only
    /// when every stored point on element is different from it.
    pub fn add_point(&mut self, point: &PointOnElement, confusion: f64) {
        let nb_pnt = self.my_points.len();
        let mut ipnt = 0usize;
        while ipnt < nb_pnt && self.my_points[ipnt].is_different(point, confusion) {
            ipnt += 1;
        }
        if ipnt >= nb_pnt {
            self.my_points.push(point.clone());
        }
    }

    /// `NbPoints()` (`cxx:82-85`).
    pub fn nb_points(&self) -> usize {
        self.my_points.len()
    }

    /// `Point(Index)` (`cxx:91-94`), 1-based; OCCT raises `OutOfRange` for
    /// `Index > NbPoints`, the port panics through the indexing operation.
    pub fn point(&self, index: usize) -> &PointOnElement {
        &self.my_points[index - 1]
    }

    /// `RemPoint(Index)` (`cxx:102-105`), 1-based.
    pub fn rem_point(&mut self, index: usize) {
        self.my_points.remove(index - 1);
    }

    /// `ClrPoints()` (`cxx:112-115`).
    pub fn clr_points(&mut self) {
        self.my_points.clear();
    }

    /// `IsLower(Point, Confusion)` (`cxx:122-126`).
    pub fn is_lower(&self, point: &PointOnHatching, confusion: f64) -> bool {
        point.base.parameter() - self.base.parameter() > confusion
    }

    /// `IsEqual(Point, Confusion)` (`cxx:132-136`).
    pub fn is_equal(&self, point: &PointOnHatching, confusion: f64) -> bool {
        (point.base.parameter() - self.base.parameter()).abs() <= confusion
    }

    /// `IsGreater(Point, Confusion)` (`cxx:142-146`).
    pub fn is_greater(&self, point: &PointOnHatching, confusion: f64) -> bool {
        self.base.parameter() - point.base.parameter() > confusion
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
