//! Interference between a projected edge and a face. Source:
//! `HLRAlgo_Interference`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_Interference.hxx:30-160`,
//! `HLRAlgo_Interference_0.cxx:24-39`).

use crate::abs::Orientation;
use crate::hlr_algo::coincidence::HlrAlgoCoincidence;
use crate::hlr_algo::intersection::HlrAlgoIntersection;

/// An interference ties an [`HlrAlgoIntersection`] to the
/// [`HlrAlgoCoincidence`] data of the hiding edge and to the three
/// orientations of the interference. Source: `HLRAlgo_Interference`.
///
/// The OCCT default constructor (`_0.cxx:24`) is `= default`, leaving the
/// three orientations uninitialized; the Rust defaults are `FORWARD`,
/// matching zero-initialized memory.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HlrAlgoInterference {
    intersection: HlrAlgoIntersection,
    boundary: HlrAlgoCoincidence,
    orientation: Orientation,
    transition: Orientation,
    b_transition: Orientation,
}

impl Default for HlrAlgoInterference {
    /// `HLRAlgo_Interference() = default` (`_0.cxx:24`).
    fn default() -> Self {
        Self {
            intersection: HlrAlgoIntersection::default(),
            boundary: HlrAlgoCoincidence::default(),
            orientation: Orientation::Forward,
            transition: Orientation::Forward,
            b_transition: Orientation::Forward,
        }
    }
}

impl HlrAlgoInterference {
    /// `HLRAlgo_Interference()` (`_0.cxx:24`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `HLRAlgo_Interference(Inters, Bound, Orient, Trans, BTrans)`
    /// (`_0.cxx:27-39`).
    pub fn from_parts(
        inters: &HlrAlgoIntersection,
        bound: &HlrAlgoCoincidence,
        orient: Orientation,
        trans: Orientation,
        b_trans: Orientation,
    ) -> Self {
        Self {
            intersection: *inters,
            boundary: *bound,
            orientation: orient,
            transition: trans,
            b_transition: b_trans,
        }
    }

    /// `Intersection(I)` (`hxx:78-81`).
    pub fn set_intersection(&mut self, i: &HlrAlgoIntersection) {
        self.intersection = *i;
    }

    /// `Boundary(B)` (`hxx:85-88`).
    pub fn set_boundary(&mut self, b: &HlrAlgoCoincidence) {
        self.boundary = *b;
    }

    /// `Orientation(O)` (`hxx:92-95`).
    pub fn set_orientation(&mut self, o: Orientation) {
        self.orientation = o;
    }

    /// `Transition(Tr)` (`hxx:99-102`).
    pub fn set_transition(&mut self, tr: Orientation) {
        self.transition = tr;
    }

    /// `BoundaryTransition(BTr)` (`hxx:106-109`).
    pub fn set_boundary_transition(&mut self, b_tr: Orientation) {
        self.b_transition = b_tr;
    }

    /// `Intersection()` (`hxx:112-116`).
    pub fn intersection(&self) -> &HlrAlgoIntersection {
        &self.intersection
    }

    /// `ChangeIntersection()` (`hxx:120-124`).
    pub fn change_intersection(&mut self) -> &mut HlrAlgoIntersection {
        &mut self.intersection
    }

    /// `Boundary()` (`hxx:126-130`).
    pub fn boundary(&self) -> &HlrAlgoCoincidence {
        &self.boundary
    }

    /// `ChangeBoundary()` (`hxx:132-136`).
    pub fn change_boundary(&mut self) -> &mut HlrAlgoCoincidence {
        &mut self.boundary
    }

    /// `Orientation()` (`hxx:140-144`).
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// `Transition()` (`hxx:146-150`).
    pub fn transition(&self) -> Orientation {
        self.transition
    }

    /// `BoundaryTransition()` (`hxx:152-157`).
    pub fn boundary_transition(&self) -> Orientation {
        self.b_transition
    }
}
