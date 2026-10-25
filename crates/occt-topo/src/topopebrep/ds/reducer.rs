//! `TopOpeBRepDS_Reducer` — the front end that runs the edge (`EIR`) and the
//! face (`FIR`) interference reducers.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Reducer.hxx:26-44`,
//! `TopOpeBRepDS_Reducer.cxx:22-48` (OCCT 8.0.0).
//!
//! `FUN_ds_FEIGb1TO0(myHDS, M)` (`Reducer.cxx:46`) is commented out in OCCT and
//! stays commented out here.

use super::eir::Eir;
use super::fir::Fir;
use super::h_data_structure::HDataStructure;
use super::tool::DShapeOn1StateMap;

/// `TopOpeBRepDS_Reducer` (`hxx:26-44`).
pub struct Reducer<'a> {
    /// `myHDS`
    hds: &'a mut HDataStructure,
}

impl<'a> Reducer<'a> {
    /// `TopOpeBRepDS_Reducer(HDS)` (`cxx:22-25`).
    pub fn new(hds: &'a mut HDataStructure) -> Self {
        Self { hds }
    }

    /// `ProcessEdgeInterferences()` (`cxx:30-34`).
    pub fn process_edge_interferences(&mut self) -> Result<(), String> {
        let mut eir = Eir::new(self.hds);
        eir.process_all_edge_interferences()
    }

    /// `ProcessFaceInterferences(M)` (`cxx:38-48`).
    pub fn process_face_interferences(
        &mut self,
        m: &DShapeOn1StateMap,
    ) -> Result<(), String> {
        // modified by NIZHNY-MZV  Tue Nov 16 16:12:15 1999
        // FUN_ds_FEIGb1TO0(myHDS,M); //xpu250199

        let mut fir = Fir::new(self.hds);
        fir.process_all_face_interferences(m)
    }
}
