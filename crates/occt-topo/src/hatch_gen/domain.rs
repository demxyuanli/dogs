//! `HatchGen_Domain`: a parametric interval of a hatching, bounded or not.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/HatchGen/HatchGen_Domain.hxx:28-90`,
//! `.lxx:18-126`, `.cxx:25-95` (OCCT 8.0.0).

use super::point_on_hatching::PointOnHatching;

/// `HatchGen_Domain` (`hxx:28-90`).
#[derive(Clone, Debug)]
pub struct Domain {
    /// `myHasFirstPoint`
    my_has_first_point: bool,
    /// `myFirstPoint`
    my_first_point: PointOnHatching,
    /// `myHasSecondPoint`
    my_has_second_point: bool,
    /// `mySecondPoint`
    my_second_point: PointOnHatching,
}

impl Default for Domain {
    fn default() -> Self {
        Self::new()
    }
}

impl Domain {
    /// `HatchGen_Domain()` (`cxx:27-31`): an infinite domain, both points are
    /// absent.
    pub fn new() -> Self {
        Self {
            my_has_first_point: false,
            my_first_point: PointOnHatching::new(),
            my_has_second_point: false,
            my_second_point: PointOnHatching::new(),
        }
    }

    /// `HatchGen_Domain(P1, P2)` (`cxx:35-41`).
    pub fn from_two_points(p1: &PointOnHatching, p2: &PointOnHatching) -> Self {
        Self {
            my_has_first_point: true,
            my_first_point: p1.clone(),
            my_has_second_point: true,
            my_second_point: p2.clone(),
        }
    }

    /// `HatchGen_Domain(P, First)` (`cxx:45-60`): a semi-infinite domain whose
    /// single point is the first one when `first` holds, the second one
    /// otherwise.
    pub fn from_point(point: &PointOnHatching, first: bool) -> Self {
        if first {
            Self {
                my_has_first_point: true,
                my_has_second_point: false,
                my_first_point: point.clone(),
                my_second_point: PointOnHatching::new(),
            }
        } else {
            Self {
                my_has_first_point: false,
                my_has_second_point: true,
                my_first_point: PointOnHatching::new(),
                my_second_point: point.clone(),
            }
        }
    }

    /// `SetPoints(P1, P2)` (`.lxx:23-32`).
    pub fn set_points(&mut self, p1: &PointOnHatching, p2: &PointOnHatching) {
        self.my_has_first_point = true;
        self.my_first_point = p1.clone();
        self.my_has_second_point = true;
        self.my_second_point = p2.clone();
    }

    /// `SetPoints()` (`.lxx:39-43`): both points at the infinite.
    pub fn set_points_infinite(&mut self) {
        self.my_has_first_point = false;
        self.my_has_second_point = false;
    }

    /// `SetFirstPoint(P)` (`.lxx:50-54`).
    pub fn set_first_point(&mut self, point: &PointOnHatching) {
        self.my_has_first_point = true;
        self.my_first_point = point.clone();
    }

    /// `SetFirstPoint()` (`.lxx:62-65`): the first point at the infinite.
    pub fn set_first_point_infinite(&mut self) {
        self.my_has_first_point = false;
    }

    /// `SetSecondPoint(P)` (`.lxx:72-76`).
    pub fn set_second_point(&mut self, point: &PointOnHatching) {
        self.my_has_second_point = true;
        self.my_second_point = point.clone();
    }

    /// `SetSecondPoint()` (`.lxx:82-85`): the second point at the infinite.
    pub fn set_second_point_infinite(&mut self) {
        self.my_has_second_point = false;
    }

    /// `HasFirstPoint()` (`.lxx:92-95`).
    pub fn has_first_point(&self) -> bool {
        self.my_has_first_point
    }

    /// `FirstPoint()` (`.lxx:102-106`). OCCT raises `Standard_DomainError`
    /// when the first point is absent; the port panics with that message.
    pub fn first_point(&self) -> &PointOnHatching {
        assert!(self.my_has_first_point, "HatchGen_Domain::FirstPoint");
        &self.my_first_point
    }

    /// `HasSecondPoint()` (`.lxx:112-115`).
    pub fn has_second_point(&self) -> bool {
        self.my_has_second_point
    }

    /// `SecondPoint()` (`.lxx:122-126`). OCCT raises `Standard_DomainError`
    /// when the second point is absent; the port panics with that message.
    pub fn second_point(&self) -> &PointOnHatching {
        assert!(self.my_has_second_point, "HatchGen_Domain::SecondPoint");
        &self.my_second_point
    }
}
