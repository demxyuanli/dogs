//! `HLRBRep_AreaLimit`
//! (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_AreaLimit.{hxx,cxx}`,
//! `hxx:26-89`, `cxx:24-153`): one vertex of a hiding edge together with the
//! four states [`crate::hlrbrep_edge_builder`] needs to walk the areas between
//! two consecutive vertices.
//!
//! `myStateBefore` / `myStateAfter` describe the visible/hidden state of the
//! **edge** on either side of the vertex (i.e. the interference transition),
//! `myEdgeBefore` / `myEdgeAfter` the state of the boundary of the hiding face.
//!
//! Deviation: OCCT links the limits with `handle<HLRBRep_AreaLimit>` pointers
//! (`hxx:86-87`), so the list can be circular (`cxx` `EdgeBuilder:147-151`).
//! The port keeps the nodes in an arena owned by the builder and stores the
//! links as `Option<usize>` arena indices, which keeps the same circular shape
//! without reference counting. `Clear()` (`cxx:148-153`) drops both links.

use crate::abs::State;
use crate::hlr_algo::HlrAlgoIntersection;

/// `HLRBRep_AreaLimit` (`hxx:26-89`).
#[derive(Debug, Clone)]
pub struct HlrBRepAreaLimit {
    /// `HLRAlgo_Intersection myVertex` (`hxx:78`).
    my_vertex: HlrAlgoIntersection,
    /// `bool myBoundary` (`hxx:79`).
    my_boundary: bool,
    /// `bool myInterference` (`hxx:80`).
    my_interference: bool,
    /// `TopAbs_State myStateBefore` (`hxx:81`).
    my_state_before: State,
    /// `TopAbs_State myStateAfter` (`hxx:82`).
    my_state_after: State,
    /// `TopAbs_State myEdgeBefore` (`hxx:83`).
    my_edge_before: State,
    /// `TopAbs_State myEdgeAfter` (`hxx:84`).
    my_edge_after: State,
    /// `handle<HLRBRep_AreaLimit> myPrevious` (`hxx:86`) as an arena index.
    my_previous: Option<usize>,
    /// `handle<HLRBRep_AreaLimit> myNext` (`hxx:87`) as an arena index.
    my_next: Option<usize>,
}

impl HlrBRepAreaLimit {
    /// `HLRBRep_AreaLimit(V, Boundary, Interference, StateBefore, StateAfter,
    /// EdgeBefore, EdgeAfter)` (`cxx:24-40`): both links null.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        v: &HlrAlgoIntersection,
        boundary: bool,
        interference: bool,
        state_before: State,
        state_after: State,
        edge_before: State,
        edge_after: State,
    ) -> Self {
        Self {
            my_vertex: *v,
            my_boundary: boundary,
            my_interference: interference,
            my_state_before: state_before,
            my_state_after: state_after,
            my_edge_before: edge_before,
            my_edge_after: edge_after,
            my_previous: None,
            my_next: None,
        }
    }

    /// `StateBefore(St)` (`cxx:44-47`).
    pub fn set_state_before(&mut self, st: State) {
        self.my_state_before = st;
    }

    /// `StateAfter(St)` (`cxx:51-54`).
    pub fn set_state_after(&mut self, st: State) {
        self.my_state_after = st;
    }

    /// `EdgeBefore(St)` (`cxx:58-61`).
    pub fn set_edge_before(&mut self, st: State) {
        self.my_edge_before = st;
    }

    /// `EdgeAfter(St)` (`cxx:65-68`).
    pub fn set_edge_after(&mut self, st: State) {
        self.my_edge_after = st;
    }

    /// `Previous(P)` (`cxx:72-75`), as an arena index.
    pub fn set_previous(&mut self, p: Option<usize>) {
        self.my_previous = p;
    }

    /// `Next(N)` (`cxx:78-81`), as an arena index.
    pub fn set_next(&mut self, n: Option<usize>) {
        self.my_next = n;
    }

    /// `Vertex()` (`cxx:86-89`).
    pub fn vertex(&self) -> &HlrAlgoIntersection {
        &self.my_vertex
    }

    /// `IsBoundary()` (`cxx:93-96`).
    pub fn is_boundary(&self) -> bool {
        self.my_boundary
    }

    /// `IsInterference()` (`cxx:100-103`).
    pub fn is_interference(&self) -> bool {
        self.my_interference
    }

    /// `StateBefore()` (`cxx:107-110`).
    pub fn state_before(&self) -> State {
        self.my_state_before
    }

    /// `StateAfter()` (`cxx:114-117`).
    pub fn state_after(&self) -> State {
        self.my_state_after
    }

    /// `EdgeBefore()` (`cxx:121-124`).
    pub fn edge_before(&self) -> State {
        self.my_edge_before
    }

    /// `EdgeAfter()` (`cxx:128-131`).
    pub fn edge_after(&self) -> State {
        self.my_edge_after
    }

    /// `Previous()` (`cxx:136-139`), as an arena index.
    pub fn previous(&self) -> Option<usize> {
        self.my_previous
    }

    /// `Next()` (`cxx:143-146`), as an arena index.
    pub fn next(&self) -> Option<usize> {
        self.my_next
    }

    /// `Clear()` (`cxx:148-153`): both links nullified.
    pub fn clear(&mut self) {
        self.my_previous = None;
        self.my_next = None;
    }
}
