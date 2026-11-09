//! The small polygon records of `TKHLR/HLRAlgo`: the triangle record, the
//! hiding-triangle record, and the face / edge flag masks.
//!
//! Sources:
//! * `HLRAlgo_PolyMask.hxx:20-33` — `HLRAlgo_PolyMask`.
//! * `HLRAlgo_TriangleData.hxx:25-31` — `HLRAlgo_TriangleData`.
//! * `HLRAlgo_PolyHidingData.hxx:25-72` — `HLRAlgo_PolyHidingData`.

use occt_core::gp::GpXyz;

/// `HLRAlgo_PolyMask` (`PolyMask.hxx:20-33`): the three outline / three
/// regular-edge bits and the eight face bits.
pub mod poly_mask {
    pub const EMSK_OUT_LIN1: i32 = 1;
    pub const EMSK_OUT_LIN2: i32 = 2;
    pub const EMSK_OUT_LIN3: i32 = 4;
    pub const EMSK_GR_A_LIN1: i32 = 8;
    pub const EMSK_GR_A_LIN2: i32 = 16;
    pub const EMSK_GR_A_LIN3: i32 = 32;
    pub const FMSK_BACK: i32 = 64;
    pub const FMSK_SIDE: i32 = 128;
    pub const FMSK_HIDING: i32 = 256;
    pub const FMSK_FLAT: i32 = 512;
    pub const FMSK_ON_OUT_L: i32 = 1024;
    pub const FMSK_OR_BACK: i32 = 2048;
    pub const FMSK_FR_BACK: i32 = 4096;
}

/// `HLRAlgo_TriangleData` (`TriangleData.hxx:25-31`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct HlrAlgoTriangleData {
    pub node1: i32,
    pub node2: i32,
    pub node3: i32,
    pub flags: i32,
}

/// `HLRAlgo_PolyHidingData::TriangleIndices` (`hxx:31-34`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TriangleIndices {
    pub index: i32,
    pub min: i32,
    pub max: i32,
}

/// `HLRAlgo_PolyHidingData::PlaneT` (`hxx:36-45`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolyHidingPlane {
    pub normal: GpXyz,
    pub d: f64,
}

impl Default for PolyHidingPlane {
    /// `PlaneT()` (`hxx:38-42`): `D = 0`, `Normal` default constructed.
    fn default() -> Self {
        Self {
            normal: GpXyz::new(0.0, 0.0, 0.0),
            d: 0.0,
        }
    }
}

/// Data structure of a set of hiding triangles. Source:
/// `HLRAlgo_PolyHidingData`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct HlrAlgoPolyHidingData {
    indices: TriangleIndices,
    plane: PolyHidingPlane,
}

impl HlrAlgoPolyHidingData {
    /// `HLRAlgo_PolyHidingData() = default` (`hxx:47`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Set(Index, Minim, Maxim, A, B, C, D)` (`hxx:49-63`).
    pub fn set(&mut self, index: i32, minim: i32, maxim: i32, a: f64, b: f64, c: f64, d: f64) {
        self.indices.index = index;
        self.indices.min = minim;
        self.indices.max = maxim;
        self.plane.normal = GpXyz::new(a, b, c);
        self.plane.d = d;
    }

    /// `Indices()` (`hxx:65`).
    pub fn indices(&self) -> &TriangleIndices {
        &self.indices
    }

    /// `Plane()` (`hxx:67`).
    pub fn plane(&self) -> &PolyHidingPlane {
        &self.plane
    }
}
