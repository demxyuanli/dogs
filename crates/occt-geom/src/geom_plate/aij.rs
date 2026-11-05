//! `GeomPlate_Aij` — the index pair plus cross product container.
//!
//! Source: `TKGeomAlgo/GeomPlate/GeomPlate_Aij.hxx`,
//! `GeomPlate_Aij.cxx:22-32` (OCCT 8.0.0).
//!
//! `GeomPlate_BuildAveragePlane` builds one of these per `(i, j)` normal pair
//! and only reads the three fields back, so the layout is the class.

use occt_core::gp::GpVec;

/// `struct GeomPlate_Aij` (`hxx:29-43`): "a structure containing indexes of
/// two normals and its cross product".
#[derive(Debug, Clone, Copy)]
pub struct Aij {
    /// `Ind1` (`hxx:40`).
    pub ind1: i32,
    /// `Ind2` (`hxx:41`).
    pub ind2: i32,
    /// `Vec` (`hxx:42`).
    pub vec: GpVec,
}

impl Aij {
    /// `GeomPlate_Aij()` (`cxx:22-26`): `Ind1 = Ind2 = 0`, `Vec` left at the
    /// `gp_Vec` default (the null vector).
    pub fn new() -> Self {
        Self {
            ind1: 0,
            ind2: 0,
            vec: GpVec::zero(),
        }
    }

    /// `GeomPlate_Aij(anInd1, anInd2, aVec)` (`cxx:28-32`).
    pub fn from_parts(ind1: i32, ind2: i32, vec: GpVec) -> Self {
        Self { ind1, ind2, vec }
    }
}
