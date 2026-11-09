//! `HLRTopoBRep_FaceData` (`ModelingAlgorithms/TKHLR/HLRTopoBRep/HLRTopoBRep_FaceData.hxx:29-55`,
//! `HLRTopoBRep_FaceData.cxx:21`, `HLRTopoBRep_FaceData.lxx:19-57`): the three
//! shape lists of one face - internal outlines, outlines on restriction and
//! iso lines.
//!
//! OCCT stores `NCollection_List<TopoDS_Shape>`; the port stores `Vec` with the
//! same append-only behaviour (the only mutator is the non-const `Add*L`,
//! `lxx:40-57`).

use crate::shape::TopoShape;

/// `HLRTopoBRep_FaceData`.
#[derive(Debug, Clone, Default)]
pub struct HlrTopoBRepFaceData {
    /// `NCollection_List<TopoDS_Shape> myIntL` (`hxx:49`).
    my_int_l: Vec<TopoShape>,
    /// `NCollection_List<TopoDS_Shape> myOutL` (`hxx:50`).
    my_out_l: Vec<TopoShape>,
    /// `NCollection_List<TopoDS_Shape> myIsoL` (`hxx:51`).
    my_iso_l: Vec<TopoShape>,
}

impl HlrTopoBRepFaceData {
    /// `HLRTopoBRep_FaceData()` (`cxx:21`): `= default`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `FaceIntL()` (`lxx:19-22`).
    pub fn face_int_l(&self) -> &[TopoShape] {
        &self.my_int_l
    }

    /// `FaceOutL()` (`lxx:26-29`).
    pub fn face_out_l(&self) -> &[TopoShape] {
        &self.my_out_l
    }

    /// `FaceIsoL()` (`lxx:33-36`).
    pub fn face_iso_l(&self) -> &[TopoShape] {
        &self.my_iso_l
    }

    /// `AddIntL()` (`lxx:40-43`).
    pub fn add_int_l(&mut self) -> &mut Vec<TopoShape> {
        &mut self.my_int_l
    }

    /// `AddOutL()` (`lxx:47-50`).
    pub fn add_out_l(&mut self) -> &mut Vec<TopoShape> {
        &mut self.my_out_l
    }

    /// `AddIsoL()` (`lxx:54-57`).
    pub fn add_iso_l(&mut self) -> &mut Vec<TopoShape> {
        &mut self.my_iso_l
    }
}
