//! `TopOpeBRepDS_PointData` — a `TopOpeBRepDS_Point` plus the two shapes it
//! separates.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_PointData.hxx:24-51`,
//! `TopOpeBRepDS_PointData.cxx:22-65`.

use super::geometry_data::GeometryData;
use super::point::Point;

/// Source: `TopOpeBRepDS_PointData.hxx:46-49`
#[derive(Debug, Clone, Default)]
pub struct PointData {
    /// The `TopOpeBRepDS_GeometryData` base.
    pub(super) geometry: GeometryData,
    /// `myPoint`
    pub(super) point: Point,
    /// `myS1`
    pub(super) s1: i32,
    /// `myS2`
    pub(super) s2: i32,
}

impl PointData {
    /// `TopOpeBRepDS_PointData()` (`cxx:25-29`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepDS_PointData(const TopOpeBRepDS_Point& P)` (`cxx:33-38`).
    pub fn with_point(p: Point) -> Self {
        Self {
            geometry: GeometryData::new(),
            point: p,
            s1: 0,
            s2: 0,
        }
    }

    /// `TopOpeBRepDS_PointData(const Point& P, const int I1, const int I2)`
    /// (`cxx:42-49`).
    pub fn with_point_and_shapes(p: Point, i1: i32, i2: i32) -> Self {
        Self {
            geometry: GeometryData::new(),
            point: p,
            s1: i1,
            s2: i2,
        }
    }

    /// The inherited `TopOpeBRepDS_GeometryData` part.
    pub fn geometry(&self) -> &GeometryData {
        &self.geometry
    }

    /// The inherited `TopOpeBRepDS_GeometryData` part, mutably.
    pub fn geometry_mut(&mut self) -> &mut GeometryData {
        &mut self.geometry
    }

    /// `SetShapes(I1, I2)` (`cxx:53-57`).
    pub fn set_shapes(&mut self, i1: i32, i2: i32) {
        self.s1 = i1;
        self.s2 = i2;
    }

    /// `GetShapes(I1, I2)` (`cxx:61-65`).
    pub fn get_shapes(&self) -> (i32, i32) {
        (self.s1, self.s2)
    }
}
