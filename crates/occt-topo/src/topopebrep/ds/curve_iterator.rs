//! `TopOpeBRepDS_CurveIterator` — the interference walker that yields curve
//! interferences.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_CurveIterator.hxx:22-58`,
//! `TopOpeBRepDS_CurveIterator.cxx:21-65`.

use std::sync::RwLockReadGuard;

use occt_geom2d::curve::Curve2d;

use super::interference::{Interference, InterferenceHandle};
use super::interference_iterator::{DerivedMatch, InterferenceIterator, InterferenceList};
use crate::abs::{Orientation, State};

/// Source: `TopOpeBRepDS_CurveIterator.hxx:22-58`.
#[derive(Debug, Clone)]
pub struct CurveIterator {
    it: InterferenceIterator,
}

impl CurveIterator {
    /// `TopOpeBRepDS_CurveIterator(L)` (`cxx:21-30`): the base ctor runs
    /// `Init(L)` (interference-iterator.cxx:43-48, unfiltered `Match()`), then
    /// the body calls `Match()` again — now with the derived
    /// `MatchInterference` in effect, so the filter does apply from the first
    /// element.
    pub fn new(l: InterferenceList) -> Self {
        let mut it = InterferenceIterator::of_list(l);
        it.set_derived_match(DerivedMatch::Curve);
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

    /// `Current()` (`cxx:43-48`): `Value()->Geometry()`.
    pub fn geometry(&self) -> i32 {
        self.it.current().geometry()
    }

    /// `Orientation(S)` (`cxx:52-58`):
    /// `Value()->Transition().Orientation(S)`; the C++ default `T` is
    /// `TopAbs_FACE` (`hxx:104-105`).
    pub fn orientation(&self, s: State) -> Orientation {
        self.it.current().transition().orientation(s, crate::abs::ShapeType::Face)
    }

    /// `PCurve()` (`cxx:62-65`): the 2d curve of the
    /// `TopOpeBRepDS_SurfaceCurveInterference`.
    ///
    /// OCCT's `down_cast` there is unchecked; `None` here means the
    /// interference is not a surface/curve one.
    pub fn p_curve(&self) -> Option<Curve2dRef> {
        self.it.current().p_curve().cloned()
    }
}

/// A cloned 2d-curve handle (`occ::handle<Geom2d_Curve>`).
pub type Curve2dRef = std::sync::Arc<dyn Curve2d>;
