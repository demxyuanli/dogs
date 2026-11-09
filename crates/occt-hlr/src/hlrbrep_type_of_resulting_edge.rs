//! `HLRBRep_TypeOfResultingEdge` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_TypeOfResultingEdge.hxx:22-39`):
//! the `typ` selector `HLRBRep_HLRToShape::InternalCompound` filters with
//! (`HLRBRep_HLRToShape.cxx:171-197`).
//!
//! The enumerators keep the OCCT values, because `InternalCompound` compares
//! `typ` against the integer literals `1`..`5` and takes `int` in the header
//! (`HLRBRep_HLRToShape.hxx:152`).

/// `HLRBRep_TypeOfResultingEdge`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum HlrBRepTypeOfResultingEdge {
    /// `HLRBRep_Undefined` (`hxx:24`): falls into the default sharp/silhouette
    /// arm of `cxx:196-198`.
    Undefined = 0,
    /// `HLRBRep_IsoLine` (`hxx:27`): isoparametric lines.
    IsoLine = 1,
    /// `HLRBRep_OutLine` (`hxx:30`): outline ("silhouette") edges.
    OutLine = 2,
    /// `HLRBRep_Rg1Line` (`hxx:33`): G1-continuity edges between two surfaces.
    Rg1Line = 3,
    /// `HLRBRep_RgNLine` (`hxx:36`): CN-continuity edges on one surface.
    RgNLine = 4,
    /// `HLRBRep_Sharp` (`hxx:39`): C0-continuity edges.
    Sharp = 5,
}

impl HlrBRepTypeOfResultingEdge {
    /// The enumerator as the `int typ` the private methods take.
    pub fn as_int(self) -> i32 {
        self as i32
    }
}
