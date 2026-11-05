//! `MAT2d_Connexion` (`MAT2d_Connexion.cxx:23-232`, `.hxx:...`): one link
//! between two elements of the `MAT2d` graph, carrying the parameters and
//! points of both sides plus their distance.
//!
//! `IsAfter` (`cxx:204-231`) orders two connexions that start on the same line
//! by item index, then by parameter, and finally by the orientation of the two
//! chords, which is what `MAT2d_MiniPath` needs to walk a `MAT_Arc`'s
//! connexions in order.

use occt_core::gp::{GpPnt2d, GpVec2d};

/// `class MAT2d_Connexion : public Standard_Transient`
/// (`MAT2d_Connexion.hxx`).
#[derive(Debug, Clone, PartialEq)]
pub struct Mat2dConnexion {
    /// `lineA` (`.hxx`).
    line_a: i32,
    /// `lineB`.
    line_b: i32,
    /// `itemA`.
    item_a: i32,
    /// `itemB`.
    item_b: i32,
    /// `distance`.
    distance: f64,
    /// `parameterOnA`.
    parameter_on_a: f64,
    /// `parameterOnB`.
    parameter_on_b: f64,
    /// `pointA`.
    point_a: GpPnt2d,
    /// `pointB`.
    point_b: GpPnt2d,
}

impl Default for Mat2dConnexion {
    /// `MAT2d_Connexion()` (`cxx:25-34`).
    fn default() -> Self {
        Self::new()
    }
}

impl Mat2dConnexion {
    /// `MAT2d_Connexion()` (`cxx:25-34`): every field zeroed.
    pub fn new() -> Self {
        Self {
            line_a: 0,
            line_b: 0,
            item_a: 0,
            item_b: 0,
            distance: 0.0,
            parameter_on_a: 0.0,
            parameter_on_b: 0.0,
            point_a: GpPnt2d::new(0.0, 0.0),
            point_b: GpPnt2d::new(0.0, 0.0),
        }
    }

    /// `MAT2d_Connexion(LineA, LineB, ItemA, ItemB, Distance, ParameterOnA,
    /// ParameterOnB, PointA, PointB)` (`cxx:37-58`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_values(
        line_a: i32,
        line_b: i32,
        item_a: i32,
        item_b: i32,
        distance: f64,
        parameter_on_a: f64,
        parameter_on_b: f64,
        point_a: GpPnt2d,
        point_b: GpPnt2d,
    ) -> Self {
        Self {
            line_a,
            line_b,
            item_a,
            item_b,
            distance,
            parameter_on_a,
            parameter_on_b,
            point_a,
            point_b,
        }
    }

    /// `IndexFirstLine()` (`cxx:62-65`).
    pub fn index_first_line(&self) -> i32 {
        self.line_a
    }

    /// `IndexSecondLine()` (`cxx:69-72`).
    pub fn index_second_line(&self) -> i32 {
        self.line_b
    }

    /// `IndexItemOnFirst()` (`cxx:76-79`).
    pub fn index_item_on_first(&self) -> i32 {
        self.item_a
    }

    /// `IndexItemOnSecond()` (`cxx:82-85`).
    pub fn index_item_on_second(&self) -> i32 {
        self.item_b
    }

    /// `ParameterOnFirst()` (`cxx:90-93`).
    pub fn parameter_on_first(&self) -> f64 {
        self.parameter_on_a
    }

    /// `ParameterOnSecond()` (`cxx:97-100`).
    pub fn parameter_on_second(&self) -> f64 {
        self.parameter_on_b
    }

    /// `PointOnFirst()` (`cxx:104-107`).
    pub fn point_on_first(&self) -> GpPnt2d {
        self.point_a
    }

    /// `PointOnSecond()` (`cxx:112-115`).
    pub fn point_on_second(&self) -> GpPnt2d {
        self.point_b
    }

    /// `Distance()` (`cxx:119-122`).
    pub fn distance(&self) -> f64 {
        self.distance
    }

    /// `IndexFirstLine(anIndex)` (`cxx:126-129`).
    pub fn set_index_first_line(&mut self, index: i32) {
        self.line_a = index;
    }

    /// `IndexSecondLine(anIndex)` (`cxx:133-136`).
    pub fn set_index_second_line(&mut self, index: i32) {
        self.line_b = index;
    }

    /// `IndexItemOnFirst(anIndex)` (`cxx:140-143`).
    pub fn set_index_item_on_first(&mut self, index: i32) {
        self.item_a = index;
    }

    /// `IndexItemOnSecond(anIndex)` (`cxx:147-150`).
    pub fn set_index_item_on_second(&mut self, index: i32) {
        self.item_b = index;
    }

    /// `ParameterOnFirst(aParameter)` (`cxx:154-157`).
    pub fn set_parameter_on_first(&mut self, parameter: f64) {
        self.parameter_on_a = parameter;
    }

    /// `ParameterOnSecond(aParameter)` (`cxx:161-164`).
    pub fn set_parameter_on_second(&mut self, parameter: f64) {
        self.parameter_on_b = parameter;
    }

    /// `PointOnFirst(aPoint)` (`cxx:168-171`).
    pub fn set_point_on_first(&mut self, point: GpPnt2d) {
        self.point_a = point;
    }

    /// `PointOnSecond(aPoint)` (`cxx:175-178`).
    pub fn set_point_on_second(&mut self, point: GpPnt2d) {
        self.point_b = point;
    }

    /// `Distance(d)` (`cxx:182-185`).
    pub fn set_distance(&mut self, d: f64) {
        self.distance = d;
    }

    /// `Reverse()` (`cxx:189-199`): the same connexion seen from the other
    /// side.
    pub fn reversed(&self) -> Self {
        Self::with_values(
            self.line_b,
            self.line_a,
            self.item_b,
            self.item_a,
            self.distance,
            self.parameter_on_b,
            self.parameter_on_a,
            self.point_b,
            self.point_a,
        )
    }

    /// `IsAfter(C2, Sense)` (`cxx:204-231`).
    pub fn is_after(&self, c2: &Self, sense: f64) -> bool {
        if self.line_a != c2.index_first_line() {
            return false;
        }
        if self.item_a > c2.index_item_on_first() {
            return true;
        } else if self.item_a == c2.index_item_on_first() {
            if self.parameter_on_a > c2.parameter_on_first() {
                return true;
            } else if self.parameter_on_a == c2.parameter_on_first() {
                let v1 = vec_between(&c2.point_on_first(), &c2.point_on_second());
                let v2 = vec_between(&self.point_a, &self.point_b);
                if v1.crossed(&v2) * sense > 0.0 {
                    return true;
                }
            }
        }
        false
    }
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.hxx`): the vector `P2 - P1`.
fn vec_between(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}
