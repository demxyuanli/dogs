//! `TopOpeBRepDS_PointIterator` — the interference walker that yields point
//! and vertex interferences, plus the parameter / orientation / config queries
//! the vertex builders read.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_PointIterator.hxx:24-72`,
//! `TopOpeBRepDS_PointIterator.cxx:22-129`.

use std::sync::RwLockReadGuard;

use super::interference::{Interference, InterferenceHandle};
use super::interference_iterator::{DerivedMatch, InterferenceIterator, InterferenceList};
use super::kind::{Config, Kind};
use crate::abs::{Orientation, State};

/// Source: `TopOpeBRepDS_PointIterator.hxx:24-72`.
#[derive(Debug, Clone)]
pub struct PointIterator {
    it: InterferenceIterator,
}

impl PointIterator {
    /// `TopOpeBRepDS_PointIterator(L)` (`cxx:22-30`): the base ctor matches
    /// unfiltered, then the body calls `Match()` again with the derived
    /// predicate in effect.
    pub fn new(l: InterferenceList) -> Self {
        let mut it = InterferenceIterator::of_list(l);
        it.set_derived_match(DerivedMatch::Point);
        it.match_filter();
        Self { it }
    }

    /// `More()`.
    pub fn more(&self) -> bool {
        self.it.more()
    }

    /// `Next()`.
    pub fn next(&mut self) {
        self.it.next();
    }

    /// `Value()`.
    pub fn value(&self) -> &InterferenceHandle {
        self.it.value()
    }

    /// Read guard on `Value()`.
    pub fn current(&self) -> RwLockReadGuard<'_, Interference> {
        self.it.current()
    }

    /// `Current()` (`cxx:42-46`): `Value()->Geometry()`.
    pub fn geometry(&self) -> i32 {
        self.it.current().geometry()
    }

    /// `Orientation(S)` (`cxx:50-56`).
    pub fn orientation(&self, s: State) -> Orientation {
        self.it.current().transition().orientation(s, crate::abs::ShapeType::Face)
    }

    /// `Parameter()` (`cxx:60-76`): the parameter of a
    /// `TopOpeBRepDS_CurvePointInterference` or
    /// `TopOpeBRepDS_EdgeVertexInterference`; every other dynamic type raises
    /// `Standard_ProgramError`.
    pub fn parameter(&self) -> Result<f64, String> {
        let i = self.it.current();
        match i.dynamic_type() {
            "TopOpeBRepDS_CurvePointInterference" | "TopOpeBRepDS_EdgeVertexInterference" => i
                .parameter()
                .ok_or_else(|| "TopOpeBRepDS_PointIterator::Parameter()".to_string()),
            _ => Err("TopOpeBRepDS_PointIterator::Parameter()".into()),
        }
    }

    /// `IsVertex()` (`cxx:80-83`).
    pub fn is_vertex(&self) -> bool {
        self.it.current().geometry_type() == Kind::Vertex
    }

    /// `IsPoint()` (`cxx:85-90`).
    pub fn is_point(&self) -> bool {
        self.it.current().geometry_type() == Kind::Point
    }

    /// `DiffOriented()` (`cxx:94-107`): the `Config` of an
    /// `TopOpeBRepDS_EdgeVertexInterference`; other types raise
    /// `Standard_ProgramError`.
    pub fn diff_oriented(&self) -> Result<bool, String> {
        let i = self.it.current();
        if i.dynamic_type() != "TopOpeBRepDS_EdgeVertexInterference" {
            return Err("TopOpeBRepDS_PointIterator::DiffOriented()".into());
        }
        Ok(i.config() == Some(Config::DiffOriented))
    }

    /// `SameOriented()` (`cxx:111-124`).
    pub fn same_oriented(&self) -> Result<bool, String> {
        let i = self.it.current();
        if i.dynamic_type() != "TopOpeBRepDS_EdgeVertexInterference" {
            return Err("TopOpeBRepDS_PointIterator::SameOriented()".into());
        }
        Ok(i.config() == Some(Config::SameOriented))
    }

    /// `Support()` (`cxx:127-129`).
    pub fn support(&self) -> i32 {
        self.it.current().support()
    }
}
