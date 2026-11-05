//! An intersection on an edge to hide. Source: `HLRAlgo_Intersection`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_Intersection.hxx`, `.lxx`,
//! `HLRAlgo_Intersection.cxx:22-50`).

use crate::abs::{Orientation, State};

/// Describes an intersection on an edge to hide: a parameter and a state
/// (`ON` = on the face, `OUT` = above the face, `IN` = under the face).
/// Source: `HLRAlgo_Intersection`.
///
/// The OCCT default constructor (`cxx:22-29`) leaves `myOrien` and `myState`
/// uninitialized; the Rust defaults are the zero-valued enum members
/// (`FORWARD` / `IN`), matching a zero-initialized allocation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HlrAlgoIntersection {
    orien: Orientation,
    seg_index: i32,
    index: i32,
    level: i32,
    param: f64,
    toler: f32,
    state: State,
}

impl Default for HlrAlgoIntersection {
    /// `HLRAlgo_Intersection()` (`cxx:22-29`).
    fn default() -> Self {
        Self {
            orien: Orientation::Forward,
            seg_index: 0,
            index: 0,
            level: 0,
            param: 0.0,
            toler: 0.0,
            state: State::In,
        }
    }
}

impl HlrAlgoIntersection {
    /// `HLRAlgo_Intersection()` (`cxx:22-29`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `HLRAlgo_Intersection(Ori, Lev, SegInd, Ind, P, Tol, S)`
    /// (`cxx:32-50`).
    pub fn from_parts(
        ori: Orientation,
        lev: i32,
        seg_ind: i32,
        ind: i32,
        p: f64,
        tol: f32,
        s: State,
    ) -> Self {
        Self {
            orien: ori,
            seg_index: seg_ind,
            index: ind,
            level: lev,
            param: p,
            toler: tol,
            state: s,
        }
    }

    /// `Orientation(Ori)` (`lxx:19-22`).
    pub fn set_orientation(&mut self, ori: Orientation) {
        self.orien = ori;
    }

    /// `Orientation()` (`lxx:26-29`).
    pub fn orientation(&self) -> Orientation {
        self.orien
    }

    /// `Level(Lev)` (`lxx:35-38`).
    pub fn set_level(&mut self, lev: i32) {
        self.level = lev;
    }

    /// `Level()` (`lxx:42-45`).
    pub fn level(&self) -> i32 {
        self.level
    }

    /// `SegIndex(SegInd)` (`lxx:49-52`).
    pub fn set_seg_index(&mut self, seg_ind: i32) {
        self.seg_index = seg_ind;
    }

    /// `SegIndex()` (`lxx:56-59`).
    pub fn seg_index(&self) -> i32 {
        self.seg_index
    }

    /// `Index(Ind)` (`lxx:63-66`).
    pub fn set_index(&mut self, ind: i32) {
        self.index = ind;
    }

    /// `Index()` (`lxx:70-73`).
    pub fn index(&self) -> i32 {
        self.index
    }

    /// `Parameter(P)` (`lxx:77-80`).
    pub fn set_parameter(&mut self, p: f64) {
        self.param = p;
    }

    /// `Parameter()` (`lxx:84-87`).
    pub fn parameter(&self) -> f64 {
        self.param
    }

    /// `Tolerance(T)` (`lxx:91-94`).
    pub fn set_tolerance(&mut self, t: f32) {
        self.toler = t;
    }

    /// `Tolerance()` (`lxx:98-101`).
    pub fn tolerance(&self) -> f32 {
        self.toler
    }

    /// `State(S)` (`lxx:105-108`).
    pub fn set_state(&mut self, s: State) {
        self.state = s;
    }

    /// `State()` (`lxx:112-115`).
    pub fn state(&self) -> State {
        self.state
    }
}
