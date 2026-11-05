//! Set of edge blocks of one wire block. Source: `HLRAlgo_WiresBlock`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_WiresBlock.hxx:27-66`,
//! `HLRAlgo_WiresBlock.cxx:19-23`).

use std::sync::Arc;

use super::edges_block::{HlrAlgoEdgesBlock, MinMaxIndices};

/// A set of [`HlrAlgoEdgesBlock`]s. Source: `HLRAlgo_WiresBlock`.
///
/// OCCT stores `NCollection_Array1<occ::handle<Standard_Transient>>` and
/// reinterprets each slot as `occ::handle<HLRAlgo_EdgesBlock>` (`hxx:49-52`);
/// the port stores the strongly typed handle directly.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct HlrAlgoWiresBlock {
    wires: Vec<Option<Arc<HlrAlgoEdgesBlock>>>,
    min_max: MinMaxIndices,
}

impl HlrAlgoWiresBlock {
    /// `HLRAlgo_WiresBlock(const int NbWires)` (`hxx:39-43`).
    pub fn new(nb_wires: usize) -> Self {
        Self {
            wires: vec![None; nb_wires],
            min_max: MinMaxIndices::default(),
        }
    }

    /// `NbWires()` (`hxx:45`): `myWires.Upper()`.
    pub fn nb_wires(&self) -> usize {
        self.wires.len()
    }

    /// `Set(I, W)` (`hxx:47`).
    pub fn set(&mut self, i: usize, w: Arc<HlrAlgoEdgesBlock>) {
        self.wires[i - 1] = Some(w);
    }

    /// `Wire(I)` (`hxx:49-52`).
    pub fn wire(&self, i: usize) -> Option<&Arc<HlrAlgoEdgesBlock>> {
        self.wires[i - 1].as_ref()
    }

    /// `UpdateMinMax(theMinMaxes)` (`hxx:54-57`).
    pub fn update_min_max(&mut self, min_max: &MinMaxIndices) {
        self.min_max = *min_max;
    }

    /// `MinMax()` (`hxx:59`).
    pub fn min_max(&self) -> &MinMaxIndices {
        &self.min_max
    }

    /// `MinMax()` mutable form (the OCCT non-const accessor).
    pub fn change_min_max(&mut self) -> &mut MinMaxIndices {
        &mut self.min_max
    }
}
