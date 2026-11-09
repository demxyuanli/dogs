//! Set of edges of one wire, with per-edge flags. Source: `HLRAlgo_EdgesBlock`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_EdgesBlock.hxx:27-162`,
//! `HLRAlgo_EdgesBlock.cxx:27-31`).

use occt_topo::abs::Orientation;

/// Flag masks of `HLRAlgo_EdgesBlock::EMskFlags` (`hxx:148-156`).
const EMASK_ORIENT: i32 = 15;
const EMASK_OUT_LINE: i32 = 16;
const EMASK_INTERNAL: i32 = 32;
const EMASK_DOUBLE: i32 = 64;
const EMASK_ISO_LINE: i32 = 128;

/// `(int)TopAbs_Orientation` of the four enum members
/// (`TopAbs_Orientation.hxx:22-28`).
fn orientation_to_bits(o: Orientation) -> i32 {
    match o {
        Orientation::Forward => 0,
        Orientation::Reversed => 1,
        Orientation::Internal => 2,
        Orientation::External => 3,
    }
}

/// `(TopAbs_Orientation)(myFlags(I) & EMaskOrient)` (`hxx:97-100`); the four
/// values above are the only ones `Orientation(I, Or)` stores.
fn bits_to_orientation(bits: i32) -> Orientation {
    match bits {
        1 => Orientation::Reversed,
        2 => Orientation::Internal,
        3 => Orientation::External,
        _ => Orientation::Forward,
    }
}

/// 16-entry min / max arrays used by the edge block pruning. Source:
/// `HLRAlgo_EdgesBlock::MinMaxIndices` (`hxx:45-80`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinMaxIndices {
    pub min: [i32; 8],
    pub max: [i32; 8],
}

impl Default for MinMaxIndices {
    /// The OCCT struct has no constructor; zero-initialized arrays.
    fn default() -> Self {
        Self {
            min: [0; 8],
            max: [0; 8],
        }
    }
}

impl MinMaxIndices {
    /// `Minimize` (`hxx:48-63`): element-wise `min`. Note the OCCT body min's
    /// **both** `Min` and `Max` arrays.
    pub fn minimize(&mut self, other: &Self) -> &mut Self {
        for i in 0..8 {
            if self.min[i] > other.min[i] {
                self.min[i] = other.min[i];
            }
            if self.max[i] > other.max[i] {
                self.max[i] = other.max[i];
            }
        }
        self
    }

    /// `Maximize` (`hxx:65-79`): element-wise `max` on both arrays.
    pub fn maximize(&mut self, other: &Self) -> &mut Self {
        for i in 0..8 {
            if self.min[i] < other.min[i] {
                self.min[i] = other.min[i];
            }
            if self.max[i] < other.max[i] {
                self.max[i] = other.max[i];
            }
        }
        self
    }
}

/// A set of edges of one wire plus its per-edge flags. Source:
/// `HLRAlgo_EdgesBlock`.
///
/// The OCCT arrays are 1-based (`myEdges(1..NbEdges)`); the accessors here
/// take the same 1-based `i` so that the call sites of `HLRBRep_Data` port
/// without index arithmetic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HlrAlgoEdgesBlock {
    edges: Vec<i32>,
    flags: Vec<i32>,
    min_max: MinMaxIndices,
}

impl HlrAlgoEdgesBlock {
    /// `HLRAlgo_EdgesBlock(const int NbEdges)` (`cxx:27-31`).
    pub fn new(nb_edges: usize) -> Self {
        Self {
            edges: vec![0; nb_edges],
            flags: vec![0; nb_edges],
            min_max: MinMaxIndices::default(),
        }
    }

    /// `NbEdges()` (`hxx:85`): `myEdges.Upper()`.
    pub fn nb_edges(&self) -> usize {
        self.edges.len()
    }

    /// `Edge(I, EI)` (`hxx:87`).
    pub fn set_edge(&mut self, i: usize, ei: i32) {
        self.edges[i - 1] = ei;
    }

    /// `Edge(I)` (`hxx:89`).
    pub fn edge(&self, i: usize) -> i32 {
        self.edges[i - 1]
    }

    /// `Orientation(I, Or)` (`hxx:91-95`): clear the 4-bit orientation field,
    /// then OR the new value in.
    pub fn set_orientation(&mut self, i: usize, or: Orientation) {
        let f = &mut self.flags[i - 1];
        *f &= !EMASK_ORIENT;
        *f |= orientation_to_bits(or) & EMASK_ORIENT;
    }

    /// `Orientation(I)` (`hxx:97-100`).
    pub fn orientation(&self, i: usize) -> Orientation {
        bits_to_orientation(self.flags[i - 1] & EMASK_ORIENT)
    }

    /// `OutLine(I)` (`hxx:102`).
    pub fn out_line(&self, i: usize) -> bool {
        self.flags[i - 1] & EMASK_OUT_LINE != 0
    }

    /// `OutLine(I, B)` (`hxx:104-110`).
    pub fn set_out_line(&mut self, i: usize, b: bool) {
        set_flag(&mut self.flags[i - 1], EMASK_OUT_LINE, b);
    }

    /// `Internal(I)` (`hxx:112`).
    pub fn internal(&self, i: usize) -> bool {
        self.flags[i - 1] & EMASK_INTERNAL != 0
    }

    /// `Internal(I, B)` (`hxx:114-120`).
    pub fn set_internal(&mut self, i: usize, b: bool) {
        set_flag(&mut self.flags[i - 1], EMASK_INTERNAL, b);
    }

    /// `Double(I)` (`hxx:122`).
    pub fn is_double(&self, i: usize) -> bool {
        self.flags[i - 1] & EMASK_DOUBLE != 0
    }

    /// `Double(I, B)` (`hxx:124-130`).
    pub fn set_double(&mut self, i: usize, b: bool) {
        set_flag(&mut self.flags[i - 1], EMASK_DOUBLE, b);
    }

    /// `IsoLine(I)` (`hxx:132`).
    pub fn iso_line(&self, i: usize) -> bool {
        self.flags[i - 1] & EMASK_ISO_LINE != 0
    }

    /// `IsoLine(I, B)` (`hxx:134-140`).
    pub fn set_iso_line(&mut self, i: usize, b: bool) {
        set_flag(&mut self.flags[i - 1], EMASK_ISO_LINE, b);
    }

    /// `UpdateMinMax(TotMinMax)` (`hxx:142`).
    pub fn update_min_max(&mut self, tot_min_max: &MinMaxIndices) {
        self.min_max = *tot_min_max;
    }

    /// `MinMax()` (`hxx:144`).
    pub fn min_max(&self) -> &MinMaxIndices {
        &self.min_max
    }

    /// `MinMax()` mutable form (the OCCT non-const accessor).
    pub fn change_min_max(&mut self) -> &mut MinMaxIndices {
        &mut self.min_max
    }
}

fn set_flag(flag: &mut i32, mask: i32, b: bool) {
    if b {
        *flag |= mask;
    } else {
        *flag &= !mask;
    }
}
