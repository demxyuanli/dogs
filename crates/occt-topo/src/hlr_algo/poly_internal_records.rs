//! The node / segment records of the outline update. Sources:
//! * `HLRAlgo_PolyInternalNode` (`HLRAlgo_PolyInternalNode.hxx:29-70`,
//!   `HLRAlgo_PolyInternalNode.cxx:18`).
//! * `HLRAlgo_PolyInternalSegment` (`HLRAlgo_PolyInternalSegment.hxx:25-33`).

use occt_core::gp::{GpXY, GpXyz};

/// `HLRAlgo_PolyInternalNode::NodeIndices` (`hxx:32-35`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PolyInternalNodeIndices {
    pub nd_sg: i32,
    pub flag: i32,
    pub edg1: i32,
    pub edg2: i32,
}

/// `HLRAlgo_PolyInternalNode::NodeData` (`hxx:37-51`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolyInternalNodeData {
    pub point: GpXyz,
    pub normal: GpXyz,
    pub uv: GpXY,
    pub pcu1: f64,
    pub pcu2: f64,
    pub scal: f64,
}

impl Default for PolyInternalNodeData {
    /// `NodeData()` (`hxx:38-44`): only `PCu1`, `PCu2`, `Scal` are zeroed;
    /// `Point`, `Normal` and `UV` are default constructed.
    fn default() -> Self {
        Self {
            point: GpXyz::new(0.0, 0.0, 0.0),
            normal: GpXyz::new(0.0, 0.0, 0.0),
            uv: GpXY::new(0.0, 0.0),
            pcu1: 0.0,
            pcu2: 0.0,
            scal: 0.0,
        }
    }
}

/// Node of the outline update. Source: `HLRAlgo_PolyInternalNode`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct HlrAlgoPolyInternalNode {
    indices: PolyInternalNodeIndices,
    data: PolyInternalNodeData,
}

impl HlrAlgoPolyInternalNode {
    /// `HLRAlgo_PolyInternalNode()` (`hxx:51-57`): zeroes the four indices.
    pub fn new() -> Self {
        Self {
            indices: PolyInternalNodeIndices::default(),
            data: PolyInternalNodeData::default(),
        }
    }

    /// `Indices()` (`hxx:60`).
    pub fn indices(&self) -> &PolyInternalNodeIndices {
        &self.indices
    }

    /// `Indices()` mutable form.
    pub fn change_indices(&mut self) -> &mut PolyInternalNodeIndices {
        &mut self.indices
    }

    /// `Data()` (`hxx:62`).
    pub fn data(&self) -> &PolyInternalNodeData {
        &self.data
    }

    /// `Data()` mutable form.
    pub fn change_data(&mut self) -> &mut PolyInternalNodeData {
        &mut self.data
    }
}

/// `HLRAlgo_PolyInternalSegment` (`PolyInternalSegment.hxx:25-33`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct HlrAlgoPolyInternalSegment {
    pub lst_sg1: i32,
    pub lst_sg2: i32,
    pub nxt_sg1: i32,
    pub nxt_sg2: i32,
    pub conex1: i32,
    pub conex2: i32,
}
