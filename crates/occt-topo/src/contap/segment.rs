//! `Contap_TheSegmentOfTheSearch` - a solution span on a restriction arc.
//!
//! Sources: `Contap_TheSegmentOfTheSearch.hxx:29-111`, `_0.cxx:24-45`.

use std::sync::Arc;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;

use super::path_point::ContapThePathPointOfTheSearch;

/// `Contap_TheSegmentOfTheSearch` (`hxx:29-66`).
#[derive(Clone)]
pub struct ContapTheSegmentOfTheSearch {
    /// `arc` (`hxx:61`).
    arc: Arc<Curve2dAdaptor>,
    /// `hasfp` (`hxx:62`).
    has_fp: bool,
    /// `thefp` (`hxx:63`).
    the_fp: ContapThePathPointOfTheSearch,
    /// `haslp` (`hxx:64`).
    has_lp: bool,
    /// `thelp` (`hxx:65`).
    the_lp: ContapThePathPointOfTheSearch,
}

impl ContapTheSegmentOfTheSearch {
    /// `Contap_TheSegmentOfTheSearch()` (`_0.cxx:24-28`).
    pub fn new() -> Self {
        Self {
            arc: Arc::new(Curve2dAdaptor::new()),
            has_fp: false,
            the_fp: ContapThePathPointOfTheSearch::new(),
            has_lp: false,
            the_lp: ContapThePathPointOfTheSearch::new(),
        }
    }

    /// `SetValue(A)` (`hxx:73-78`): clears both limits.
    pub fn set_value(&mut self, a: Arc<Curve2dAdaptor>) {
        self.has_fp = false;
        self.has_lp = false;
        self.arc = a;
    }

    /// `SetLimitPoint(V, First)` (`_0.cxx:30-45`).
    pub fn set_limit_point(&mut self, v: &ContapThePathPointOfTheSearch, first: bool) {
        if first {
            self.has_fp = true;
            self.the_fp = v.clone();
        } else {
            self.has_lp = true;
            self.the_lp = v.clone();
        }
    }

    /// `Curve()` (`hxx:79-82`).
    pub fn curve(&self) -> &Arc<Curve2dAdaptor> {
        &self.arc
    }

    /// `HasFirstPoint()` (`hxx:84-87`).
    pub fn has_first_point(&self) -> bool {
        self.has_fp
    }

    /// `FirstPoint()` (`hxx:89-95`): `None` where OCCT raises
    /// `Standard_DomainError`.
    pub fn first_point(&self) -> Option<&ContapThePathPointOfTheSearch> {
        if self.has_fp {
            Some(&self.the_fp)
        } else {
            None
        }
    }

    /// `HasLastPoint()` (`hxx:97-100`).
    pub fn has_last_point(&self) -> bool {
        self.has_lp
    }

    /// `LastPoint()` (`hxx:102-110`): `None` where OCCT raises
    /// `Standard_DomainError`.
    pub fn last_point(&self) -> Option<&ContapThePathPointOfTheSearch> {
        if self.has_lp {
            Some(&self.the_lp)
        } else {
            None
        }
    }
}

impl Default for ContapTheSegmentOfTheSearch {
    fn default() -> Self {
        Self::new()
    }
}
