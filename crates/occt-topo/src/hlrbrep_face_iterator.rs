//! `HLRBRep_FaceIterator` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_FaceIterator.hxx:27-90`,
//! `HLRBRep_FaceIterator.cxx:25-60`, `HLRBRep_FaceIterator.lxx:21-91`): walks the
//! `HLRAlgo_WiresBlock` / `HLRAlgo_EdgesBlock` structure of an `HLRBRep_FaceData`
//! in wire-major, edge-minor order.
//!
//! The OCCT iterator keeps raw handles to the blocks; the port keeps `Arc`
//! clones of the same blocks (the wire and edge blocks are immutable for the
//! duration of an exploration).

use std::sync::Arc;

use crate::abs::Orientation;
use crate::hlr_algo::{HlrAlgoEdgesBlock, HlrAlgoWiresBlock};
use crate::hlrbrep_face_data::HlrBRepFaceData;

/// `HLRBRep_FaceIterator`.
#[derive(Debug, Clone, Default)]
pub struct HlrBRepFaceIterator {
    /// `int iWire` (`hxx:77`).
    i_wire: i32,
    /// `int nbWires` (`hxx:78`).
    nb_wires: i32,
    /// `int iEdge` (`hxx:79`).
    i_edge: i32,
    /// `int nbEdges` (`hxx:80`).
    nb_edges: i32,
    /// `occ::handle<HLRAlgo_WiresBlock> myWires` (`hxx:81`).
    my_wires: Option<Arc<HlrAlgoWiresBlock>>,
    /// `occ::handle<HLRAlgo_EdgesBlock> myEdges` (`hxx:82`).
    my_edges: Option<Arc<HlrAlgoEdgesBlock>>,
}

impl HlrBRepFaceIterator {
    /// `HLRBRep_FaceIterator()` (`cxx:25`): `= default`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `InitEdge(fd)` (`cxx:29-40`).
    pub fn init_edge(&mut self, fd: &HlrBRepFaceData) {
        self.i_wire = 0;
        self.my_wires = fd.wires().cloned();
        self.nb_wires = self
            .my_wires
            .as_ref()
            .map(|w| w.nb_wires() as i32)
            .unwrap_or(0);

        self.i_edge = 0;
        self.nb_edges = 0;
        self.next_edge();
    }

    /// `MoreEdge()` (`lxx:21-24`).
    pub fn more_edge(&self) -> bool {
        self.i_wire <= self.nb_wires
    }

    /// `NextEdge()` (`cxx:44-56`).
    pub fn next_edge(&mut self) {
        self.i_edge += 1;
        if self.i_edge > self.nb_edges {
            self.i_wire += 1;
            if self.i_wire <= self.nb_wires && self.i_wire >= 1 {
                self.i_edge = 1;
                self.my_edges = self
                    .my_wires
                    .as_ref()
                    .and_then(|w| w.wire(self.i_wire as usize))
                    .cloned();
                self.nb_edges = self
                    .my_edges
                    .as_ref()
                    .map(|e| e.nb_edges() as i32)
                    .unwrap_or(0);
            }
        }
    }

    /// `BeginningOfWire()` (`lxx:29-32`).
    pub fn beginning_of_wire(&self) -> bool {
        self.i_edge == 1
    }

    /// `EndOfWire()` (`lxx:37-40`).
    pub fn end_of_wire(&self) -> bool {
        self.i_edge == self.nb_edges
    }

    /// `SkipWire()` (`lxx:45-49`).
    pub fn skip_wire(&mut self) {
        self.i_edge = self.nb_edges;
        self.next_edge();
    }

    /// `Wire()` (`lxx:54-57`).
    pub fn wire(&self) -> Option<Arc<HlrAlgoEdgesBlock>> {
        self.my_wires
            .as_ref()
            .and_then(|w| w.wire(self.i_wire as usize))
            .cloned()
    }

    /// `Edge()` (`lxx:62-65`).
    pub fn edge(&self) -> i32 {
        self.my_edges
            .as_ref()
            .map(|e| e.edge(self.i_edge as usize))
            .unwrap_or(0)
    }

    /// `Orientation()` (`lxx:70-73`).
    pub fn orientation(&self) -> Orientation {
        self.my_edges
            .as_ref()
            .map(|e| e.orientation(self.i_edge as usize))
            .unwrap_or(Orientation::Forward)
    }

    /// `OutLine()` (`lxx:78-81`).
    pub fn out_line(&self) -> bool {
        self.my_edges
            .as_ref()
            .map(|e| e.out_line(self.i_edge as usize))
            .unwrap_or(false)
    }

    /// `Internal()` (`lxx:86-89`).
    pub fn internal(&self) -> bool {
        self.my_edges
            .as_ref()
            .map(|e| e.internal(self.i_edge as usize))
            .unwrap_or(false)
    }

    /// `Double()` (`lxx:94-97`).
    pub fn double(&self) -> bool {
        self.my_edges
            .as_ref()
            .map(|e| e.is_double(self.i_edge as usize))
            .unwrap_or(false)
    }

    /// `IsoLine()` (`lxx:102-105`).
    pub fn iso_line(&self) -> bool {
        self.my_edges
            .as_ref()
            .map(|e| e.iso_line(self.i_edge as usize))
            .unwrap_or(false)
    }
}
