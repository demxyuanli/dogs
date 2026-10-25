//! `TopOpeBRepDS_SurfaceIterator` — the interference walker that yields
//! surface interferences.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_SurfaceIterator.hxx:22-48`,
//! `TopOpeBRepDS_SurfaceIterator.cxx:21-45`.

use std::sync::RwLockReadGuard;

use super::interference::{Interference, InterferenceHandle};
use super::interference_iterator::{InterferenceIterator, InterferenceList};
use super::kind::Kind;
use crate::abs::{Orientation, State};

/// Source: `TopOpeBRepDS_SurfaceIterator.hxx:22-48`.
#[derive(Debug, Clone)]
pub struct SurfaceIterator {
    it: InterferenceIterator,
}

impl SurfaceIterator {
    /// `TopOpeBRepDS_SurfaceIterator(L)` (`cxx:21-27`).
    ///
    /// The ctor only calls `InterferenceIterator::GeometryKind(SURFACE)` and
    /// **does not** re-run `Match()` afterwards, while the base ctor's `Init`
    /// already matched unfiltered (the cursor stays on element 1). The filter
    /// therefore takes effect only from the second element on — an OCCT quirk
    /// this port keeps.
    pub fn new(l: InterferenceList) -> Self {
        let mut it = InterferenceIterator::of_list(l);
        it.geometry_kind(Kind::Surface);
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

    /// `Current()` (`cxx:31-36`): `Value()->Geometry()`.
    pub fn geometry(&self) -> i32 {
        self.it.current().geometry()
    }

    /// `Orientation(S)` (`cxx:40-45`).
    pub fn orientation(&self, s: State) -> Orientation {
        self.it.current().transition().orientation(s, crate::abs::ShapeType::Face)
    }
}
