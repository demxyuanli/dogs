//! `HLRBRep_VertexList`
//! (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_VertexList.{hxx,cxx}`,
//! `hxx:32-78`, `cxx:27-155`): the merge of the two vertices of the edge being
//! hidden with the sorted interference list, in parameter order.
//!
//! The list starts on the first of the two sources (`cxx:36-37`, the
//! constructor's `InitVertices()` + `Next()`), and `Next()` (`cxx:57-83`)
//! advances both cursors then keeps only the earlier one; when the tool says
//! the interference sits on the current vertex
//! ([`HlrBRepEdgeInterferenceTool::same_vertex_and_interference`]) both stay
//! current and `Current()` reports the boundary vertex, exactly as OCCT does.
//!
//! Deviations:
//! * `myIterator` (`hxx:74`) is an `NCollection_List` iterator, i.e. a live
//!   cursor into the caller's list. The port takes a snapshot of the
//!   interference list (`HlrBRep_Hider` never mutates it after the list is
//!   handed over) and keeps an index into it.
//! * `Current()` / `Orientation()` / `Transition()` / `BoundaryTransition()`
//!   raise `Standard_NoSuchObject` / `Standard_DomainError` in OCCT; the port
//!   returns `None`.
//! * `myTool` (`hxx:75`) is a *copy* of the interference tool in OCCT too, so
//!   the port stores a clone.

use crate::abs::Orientation;
use crate::hlr_algo::{HlrAlgoInterference, HlrAlgoIntersection};
use crate::hlrbrep_edge_interference_tool::HlrBRepEdgeInterferenceTool;

/// `HLRBRep_VertexList` (`hxx:32-78`).
#[derive(Debug, Clone)]
pub struct HlrBRepVertexList {
    /// `NCollection_List<HLRAlgo_Interference>::Iterator myIterator`
    /// (`hxx:74`), as an index into [`Self::interferences`].
    my_iterator: usize,
    /// The snapshot of the list the OCCT iterator walks.
    interferences: Vec<HlrAlgoInterference>,
    /// `HLRBRep_EdgeInterferenceTool myTool` (`hxx:75`).
    my_tool: HlrBRepEdgeInterferenceTool,
    /// `bool fromEdge` (`hxx:76`).
    from_edge: bool,
    /// `bool fromInterf` (`hxx:77`).
    from_interf: bool,
}

impl HlrBRepVertexList {
    /// `HLRBRep_VertexList(T, I)` (`cxx:27-38`): copies the tool and the
    /// iterator, resets the tool's vertex cursor and advances once.
    pub fn new(t: &HlrBRepEdgeInterferenceTool, interferences: &[HlrAlgoInterference]) -> Self {
        let mut vlist = Self {
            my_iterator: 0,
            interferences: interferences.to_vec(),
            my_tool: t.clone(),
            from_edge: false,
            from_interf: false,
        };
        vlist.my_tool.init_vertices();
        vlist.next();
        vlist
    }

    /// `IsPeriodic()` (`cxx:42-45`).
    pub fn is_periodic(&self) -> bool {
        self.my_tool.is_periodic()
    }

    /// `More()` (`cxx:49-52`).
    pub fn more(&self) -> bool {
        self.from_edge || self.from_interf
    }

    /// `Next()` (`cxx:57-83`).
    pub fn next(&mut self) {
        if self.from_interf {
            self.my_iterator += 1;
        }
        if self.from_edge {
            self.my_tool.next_vertex();
        }
        self.from_interf = self.my_iterator < self.interferences.len();
        self.from_edge = self.my_tool.more_vertices();
        if self.from_edge && self.from_interf {
            let intf = self.interferences[self.my_iterator];
            if !self.my_tool.same_vertex_and_interference(&intf) {
                if self.my_tool.current_parameter()
                    < self.my_tool.parameter_of_interference(&intf)
                {
                    self.from_interf = false;
                } else {
                    self.from_edge = false;
                }
            }
        }
    }

    /// `Current()` (`cxx:87-99`).
    pub fn current(&self) -> Option<&HlrAlgoIntersection> {
        if self.from_edge {
            Some(self.my_tool.current_vertex())
        } else if self.from_interf {
            Some(self.interferences[self.my_iterator].intersection())
        } else {
            None
        }
    }

    /// `IsBoundary()` (`cxx:105-108`).
    pub fn is_boundary(&self) -> bool {
        self.from_edge
    }

    /// `IsInterference()` (`cxx:112-115`).
    pub fn is_interference(&self) -> bool {
        self.from_interf
    }

    /// `Orientation()` (`cxx:119-127`).
    pub fn orientation(&self) -> Option<Orientation> {
        if self.from_edge {
            Some(self.my_tool.current_orientation())
        } else {
            None
        }
    }

    /// `Transition()` (`cxx:133-143`).
    pub fn transition(&self) -> Option<Orientation> {
        if self.from_interf {
            Some(self.interferences[self.my_iterator].transition())
        } else {
            None
        }
    }

    /// `BoundaryTransition()` (`cxx:147-155`).
    pub fn boundary_transition(&self) -> Option<Orientation> {
        if self.from_interf {
            Some(self.interferences[self.my_iterator].boundary_transition())
        } else {
            None
        }
    }
}
