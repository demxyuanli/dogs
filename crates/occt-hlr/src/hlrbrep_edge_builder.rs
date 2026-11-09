//! `HLRBRep_EdgeBuilder`
//! (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_EdgeBuilder.{hxx,cxx}`,
//! `hxx:31-115`, `cxx:30-521`): splits one hiding edge into the alternating
//! areas "state `toBuild`" / "state not `toBuild`" delimited by the vertices of
//! a [`HlrBRepVertexList`], and yields the boundary vertices of each area as a
//! `(start, end)` pair.
//!
//! Construction (`cxx:30-230`) turns every vertex into an
//! [`HlrBRepAreaLimit`] carrying the four states of `cxx:50-121`, links them
//! into a (circular for a periodic edge) list, resolves the `UNKNOWN` states by
//! walking the list once (`cxx:156-225`) and then calls `Builds(TopAbs_IN)`.
//!
//! The area cursor is OCCT's `left`/`right` handle pair:
//! * `InitAreas` puts `right` on the head so the first area is the one that
//!   wraps the seam of a circular list (`cxx:235-238`);
//! * `HasArea` stops when both ends are exhausted or `right` came back to the
//!   head (`cxx:264-278`);
//! * `AreaState` / `AreaEdgeState` read `(left->After, right->Before)`
//!   (`cxx:282-310`);
//! * the vertex cursor `current` counts `1` = left limit, `2` = right limit,
//!   `3` = no vertex (`cxx:339-347`, `:399-421`).
//!
//! Deviation: the OCCT list is `handle<HLRBRep_AreaLimit>` with `Next()` /
//! `Previous()` pointers; the port owns all nodes in `limits` and stores the
//! links as `Option<usize>` indices (see [`HlrBRepAreaLimit`]). The three
//! `Standard_NoSuchObject` / `Standard_DomainError` raises become `None` / `Err`.

use occt_topo::abs::{Orientation, State};
use crate::hlr_algo::HlrAlgoIntersection;
use crate::hlrbrep_area_limit::HlrBRepAreaLimit;
use crate::hlrbrep_vertex_list::HlrBRepVertexList;

/// `HLRBRep_EdgeBuilder` (`hxx:31-115`).
#[derive(Debug, Clone)]
pub struct HlrBRepEdgeBuilder {
    /// `TopAbs_State toBuild` (`hxx:108`).
    to_build: State,
    /// `handle<HLRBRep_AreaLimit> myLimits` (`hxx:109`), as an arena index.
    my_limits: usize,
    /// The arena behind `myLimits` / `left` / `right`.
    limits: Vec<HlrBRepAreaLimit>,
    /// `handle<HLRBRep_AreaLimit> left` (`hxx:110`).
    left: Option<usize>,
    /// `handle<HLRBRep_AreaLimit> right` (`hxx:111`).
    right: Option<usize>,
    /// `int current` (`hxx:112`): `1` left, `2` right, `3` none.
    current: i32,
}

impl HlrBRepEdgeBuilder {
    /// `HLRBRep_EdgeBuilder(VList)` (`cxx:30-230`).
    ///
    /// Returns `Err` where OCCT raises `Standard_DomainError`: an empty vertex
    /// list (`cxx:39`) or no interference at all (`cxx:186`).
    pub fn new(vlist: &mut HlrBRepVertexList) -> Result<Self, &'static str> {
        if !vlist.more() {
            return Err("HLRBRep_EdgeBuilder: empty vertex list");
        }

        let mut limits: Vec<HlrBRepAreaLimit> = Vec::new();
        let mut last: Option<usize> = None;

        while vlist.more() {
            let mut before = State::Unknown;
            let mut after = State::Unknown;
            let mut ebefore = State::Unknown;
            let mut eafter = State::Unknown;

            if vlist.is_boundary() {
                match vlist.orientation().ok_or("EdgeBuilder: no boundary orientation")? {
                    Orientation::Forward => {
                        ebefore = State::Out;
                        eafter = State::In;
                    }
                    Orientation::Reversed => {
                        ebefore = State::In;
                        eafter = State::Out;
                    }
                    Orientation::Internal => {
                        ebefore = State::In;
                        eafter = State::In;
                    }
                    Orientation::External => {
                        ebefore = State::Out;
                        eafter = State::Out;
                    }
                }
            }

            if vlist.is_interference() {
                match vlist.transition().ok_or("EdgeBuilder: no transition")? {
                    Orientation::Forward => {
                        before = State::Out;
                        after = State::In;
                    }
                    Orientation::Reversed => {
                        before = State::In;
                        after = State::Out;
                    }
                    Orientation::Internal => {
                        before = State::In;
                        after = State::In;
                    }
                    Orientation::External => {
                        before = State::Out;
                        after = State::Out;
                    }
                }

                match vlist
                    .boundary_transition()
                    .ok_or("EdgeBuilder: no boundary transition")?
                {
                    Orientation::Forward => after = State::On,
                    Orientation::Reversed => before = State::On,
                    Orientation::Internal => {
                        before = State::On;
                        after = State::On;
                    }
                    Orientation::External => {}
                }
            }

            let v = *vlist.current().ok_or("EdgeBuilder: no current vertex")?;
            let cur = limits.len();
            limits.push(HlrBRepAreaLimit::new(
                &v,
                vlist.is_boundary(),
                vlist.is_interference(),
                before,
                after,
                ebefore,
                eafter,
            ));
            if let Some(l) = last {
                limits[l].set_next(Some(cur));
                limits[cur].set_previous(Some(l));
            }
            last = Some(cur);
            vlist.next();
        }

        let first = 0usize;
        // periodicity: close the list (`cxx:146-151`)
        if vlist.is_periodic() {
            let l = last.ok_or("EdgeBuilder: empty vertex list")?;
            limits[l].set_next(Some(first));
            limits[first].set_previous(Some(l));
        }

        // process UNKNOWN areas (`cxx:153-183`)
        let mut stat = State::Unknown;
        let mut estat = State::Unknown;
        let mut cur = Some(first);
        while let Some(c) = cur {
            if stat == State::Unknown {
                stat = limits[c].state_before();
                if stat == State::Unknown {
                    stat = limits[c].state_after();
                }
            }
            if estat == State::Unknown {
                estat = limits[c].edge_before();
                if estat == State::Unknown {
                    estat = limits[c].edge_after();
                }
            }
            cur = limits[c].next();
            if cur == Some(first) {
                break;
            }
        }

        // error if no interferences (`cxx:186`)
        if stat == State::Unknown {
            return Err("HLRBRep_EdgeBuilder: no interferences");
        }
        // if no boundary the edge covers the whole curve (`cxx:187-190`)
        if estat == State::Unknown {
            estat = State::In;
        }

        // propagate states (`cxx:192-226`)
        let mut cur = Some(first);
        while let Some(c) = cur {
            if limits[c].state_before() == State::Unknown {
                limits[c].set_state_before(stat);
            } else {
                stat = limits[c].state_after();
            }
            if limits[c].state_after() == State::Unknown {
                limits[c].set_state_after(stat);
            }
            if limits[c].edge_before() == State::Unknown {
                limits[c].set_edge_before(estat);
            } else {
                estat = limits[c].edge_after();
            }
            if limits[c].edge_after() == State::Unknown {
                limits[c].set_edge_after(estat);
            }

            cur = limits[c].next();
            if cur == Some(first) {
                break;
            }
        }

        let mut eb = Self {
            to_build: State::Unknown,
            my_limits: first,
            limits,
            left: None,
            right: None,
            current: 0,
        };
        // initialise with IN parts (`cxx:229`)
        eb.builds(State::In);
        Ok(eb)
    }

    /// `InitAreas()` (`cxx:233-238`).
    pub fn init_areas(&mut self) {
        self.left = self.limits[self.my_limits].previous();
        self.right = Some(self.my_limits);
    }

    /// `NextArea()` (`cxx:242-248`).
    pub fn next_area(&mut self) {
        self.left = self.right;
        if let Some(r) = self.right {
            self.right = self.limits[r].next();
        }
    }

    /// `PreviousArea()` (`cxx:252-260`).
    pub fn previous_area(&mut self) {
        self.right = self.left;
        if let Some(l) = self.left {
            self.left = self.limits[l].previous();
        }
    }

    /// `HasArea()` (`cxx:264-278`).
    pub fn has_area(&self) -> bool {
        if self.left.is_none() && self.right.is_none() {
            return false;
        }
        if self.right == Some(self.my_limits) {
            return false;
        }
        true
    }

    /// `AreaState()` (`cxx:282-294`).
    pub fn area_state(&self) -> State {
        let mut stat = State::Unknown;
        if let Some(l) = self.left {
            stat = self.limits[l].state_after();
        }
        if let Some(r) = self.right {
            stat = self.limits[r].state_before();
        }
        stat
    }

    /// `AreaEdgeState()` (`cxx:298-310`).
    pub fn area_edge_state(&self) -> State {
        let mut stat = State::Unknown;
        if let Some(l) = self.left {
            stat = self.limits[l].edge_after();
        }
        if let Some(r) = self.right {
            stat = self.limits[r].edge_before();
        }
        stat
    }

    /// `LeftLimit()` (`cxx:314-317`).
    pub fn left_limit(&self) -> Option<usize> {
        self.left
    }

    /// `RightLimit()` (`cxx:321-324`).
    pub fn right_limit(&self) -> Option<usize> {
        self.right
    }

    /// `Builds(ToBuild)` (`cxx:328-349`).
    pub fn builds(&mut self, to_build: State) {
        self.to_build = to_build;
        self.init_areas();
        loop {
            if self.area_state() == to_build && self.area_edge_state() == State::In {
                self.current = if self.left.is_none() { 2 } else { 1 };
                return;
            }
            self.next_area();
            if !self.has_area() {
                break;
            }
        }
        self.current = 3;
    }

    /// `MoreEdges()` (`cxx:353-356`).
    pub fn more_edges(&self) -> bool {
        self.has_area()
    }

    /// `NextEdge()` (`cxx:360-386`).
    pub fn next_edge(&mut self) {
        // clean the current edge
        while self.area_state() == self.to_build {
            self.next_area();
        }
        // go to the next edge
        while self.has_area() {
            if self.area_state() == self.to_build && self.area_edge_state() == State::In {
                self.current = if self.left.is_none() { 2 } else { 1 };
                return;
            }
            self.next_area();
        }
    }

    /// `MoreVertices()` (`cxx:390-393`).
    pub fn more_vertices(&self) -> bool {
        self.current < 3
    }

    /// `NextVertex()` (`cxx:397-422`).
    pub fn next_vertex(&mut self) -> Option<()> {
        if self.current == 1 {
            self.current = 2;
            if self.right.is_none() {
                self.current = 3;
            }
        } else if self.current == 2 {
            self.next_area();
            if self.area_state() == self.to_build && self.area_edge_state() == State::In {
                self.current = 2;
            } else {
                self.current = 3;
            }
        } else {
            // `Standard_NoSuchObject("EdgeBuilder::NextVertex ...")`
            return None;
        }
        Some(())
    }

    /// `Current()` (`cxx:426-438`).
    pub fn current(&self) -> Option<&HlrAlgoIntersection> {
        if self.current == 1 {
            self.left.map(|l| self.limits[l].vertex())
        } else if self.current == 2 {
            self.right.map(|r| self.limits[r].vertex())
        } else {
            None
        }
    }

    /// `IsBoundary()` (`cxx:442-457`).
    pub fn is_boundary(&self) -> Option<bool> {
        if self.current == 1 {
            self.left.map(|l| self.limits[l].is_boundary())
        } else if self.current == 2 {
            self.right.map(|r| self.limits[r].is_boundary())
        } else {
            None
        }
    }

    /// `IsInterference()` (`cxx:461-474`).
    pub fn is_interference(&self) -> Option<bool> {
        if self.current == 1 {
            self.left.map(|l| self.limits[l].is_interference())
        } else if self.current == 2 {
            self.right.map(|r| self.limits[r].is_interference())
        } else {
            None
        }
    }

    /// `Orientation()` (`cxx:478-505`): `INTERNAL` when the vertex does not
    /// change the state, `FORWARD` on the left limit and `REVERSED` on the
    /// right one, `EXTERNAL` when there is no current vertex.
    pub fn orientation(&self) -> Orientation {
        if self.current == 1 {
            if let Some(l) = self.left {
                let lim = &self.limits[l];
                if lim.state_before() == lim.state_after()
                    && lim.edge_before() == lim.edge_after()
                {
                    return Orientation::Internal;
                }
            }
            Orientation::Forward
        } else if self.current == 2 {
            if let Some(r) = self.right {
                let lim = &self.limits[r];
                if lim.state_before() == lim.state_after()
                    && lim.edge_before() == lim.edge_after()
                {
                    return Orientation::Internal;
                }
            }
            Orientation::Reversed
        } else {
            Orientation::External
        }
    }

    /// `Destroy()` (`cxx:509-521`): unlink every node and reset the cursor.
    pub fn destroy(&mut self) {
        for lim in self.limits.iter_mut() {
            lim.clear();
        }
        self.left = None;
        self.right = None;
        self.limits.clear();
        self.my_limits = 0;
    }
}
