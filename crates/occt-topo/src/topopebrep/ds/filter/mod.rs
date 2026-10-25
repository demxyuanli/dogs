//! `TopOpeBRepDS_Filter` — the per-shape interference clean-up run by
//! `TopOpeBRep_DSFiller::Filter`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Filter.hxx:26-84`,
//! `TopOpeBRepDS_Filter.cxx:22-93`,
//! `TopOpeBRepDS_FilterEdgeInterferences.cxx:46-89`,
//! `TopOpeBRepDS_FilterFaceInterferences.cxx:37-68`,
//! `TopOpeBRepDS_FilterCurveInterferences.cxx:21-64` (OCCT 8.0.0).
//!
//! Every `Process*Interferences(SIX)` works on the live
//! `BDS.ChangeShapeInterferences(SIX)` list. The port snapshots that list,
//! runs the same `FUN_*` helpers on the snapshot, and writes the result back;
//! the helpers themselves stay one-to-one with the `.cxx` functions.

pub mod edge_tools;
pub mod face_tools;
pub mod funk;

use std::cell::RefCell;
use std::rc::Rc;

use crate::abs::{ShapeType, State};

use super::h_data_structure::HDataStructure;
use super::interference_iterator::InterferenceList;
use super::kind::Kind;
use super::process_interferences_tool::{
    fun_reducedoublons, fun_select_gk_interference, fun_select_sk_interference,
    fun_select_trasha_interference, fun_select_traunk_interference,
};
use super::tool::DShapeOn1StateMap;
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;

/// `TopOpeBRepDS_Filter` (`hxx:26-84`).
pub struct Filter<'a> {
    /// `myHDS`
    hds: &'a mut HDataStructure,
    /// `myPShapeClassif`
    p_shape_classif: Option<Rc<RefCell<ShapeClassifier>>>,
}

impl<'a> Filter<'a> {
    /// `TopOpeBRepDS_Filter(HDS, pClassif)` (`cxx:22-28`).
    pub fn new(
        hds: &'a mut HDataStructure,
        p_shape_classif: Option<Rc<RefCell<ShapeClassifier>>>,
    ) -> Self {
        Self {
            hds,
            p_shape_classif,
        }
    }

    /// `ProcessInterferences()` (`cxx:32-36`): edges, then curves.
    ///
    /// `connex` is `TopOpeBRepDS_connex.cxx`'s global connexity, which this port
    /// threads explicitly (see `Connexity::prepare`).
    pub fn process_interferences(&mut self, connex: &super::connex::Connexity) {
        self.process_edge_interferences(connex);
        self.process_curve_interferences();
    }

    /// `ProcessEdgeInterferences()` (`cxx:38-57`): every `EDGE` shape of the DS.
    pub fn process_edge_interferences(&mut self, connex: &super::connex::Connexity) {
        let nshape = self.hds.ds().nb_shapes();
        for i in 1..=nshape {
            let Some(s) = self.hds.ds().shape(i, true) else {
                continue;
            };
            if s.shape_type() == ShapeType::Edge {
                self.process_edge_interferences_at(i, connex);
            }
        }
    }

    /// `ProcessEdgeInterferences(const int SIX)` (`FilterEdgeInterferences.cxx:64-90`).
    pub fn process_edge_interferences_at(
        &mut self,
        six: i32,
        connex: &super::connex::Connexity,
    ) {
        let Some(s) = self.hds.ds().shape(six, true) else {
            return;
        };
        let mut li: InterferenceList = self.hds.ds().shape_interferences(&s, true).to_vec();

        let mut lall = li.clone();
        // xpu : 270398 : delete the faulty interferences attached to the section
        // edge EIX.
        edge_tools::fun_purge_ds_on_se(&*self.hds, connex, six, &mut lall);
        let mut l_uu = InterferenceList::new();
        fun_select_traunk_interference(&mut lall, &mut l_uu);
        edge_tools::fun_resolve_e_unknown(&mut l_uu, self.hds.ds(), six);

        let mut lw = InterferenceList::new();
        lw.extend(lall);
        lw.extend(l_uu);
        let mut l_f = InterferenceList::new();
        let mut l_e = InterferenceList::new();
        fun_select_trasha_interference(&mut lw, ShapeType::Face, &mut l_f);
        fun_select_trasha_interference(&mut lw, ShapeType::Edge, &mut l_e);
        fun_filter_edge(&mut l_f, &*self.hds, six);
        // FUN_unkeepEsymetrictransitions(lF, BDS, SIX); // commented in OCCT
        fun_filter_edge(&mut l_e, &*self.hds, six);

        li.clear();
        li.extend(l_f);
        li.extend(l_e);
        *self.hds.change_ds().change_shape_interferences(&s) = li;
    }

    /// `ProcessFaceInterferences(MEsp)` (`cxx:59-81`): every `FACE` shape.
    pub fn process_face_interferences(&mut self, me_sp: &DShapeOn1StateMap) {
        let nshape = self.hds.ds().nb_shapes();
        for i in 1..=nshape {
            let Some(s) = self.hds.ds().shape(i, true) else {
                continue;
            };
            if s.shape_type() == ShapeType::Face {
                self.process_face_interferences_at(i, me_sp);
            }
        }
    }

    /// `ProcessFaceInterferences(const int SIX, MEsp)`
    /// (`FilterFaceInterferences.cxx:36-70`).
    pub fn process_face_interferences_at(&mut self, six: i32, me_sp: &DShapeOn1StateMap) {
        let Some(s) = self.hds.ds().shape(six, true) else {
            return;
        };
        let ds = self.hds.ds();
        let mut li: InterferenceList = ds.shape_interferences(&s, true).to_vec();
        fun_reducedoublons(&mut li, ds, six);

        let mut lall = li.clone();
        let mut l_uu = InterferenceList::new();
        fun_select_traunk_interference(&mut lall, &mut l_uu);
        face_tools::fun_resolve_f_unknown(
            &mut l_uu,
            ds,
            six,
            me_sp,
            self.p_shape_classif.clone(),
        );

        let mut lw = InterferenceList::new();
        lw.extend(lall);
        lw.extend(l_uu);
        let mut l_f = InterferenceList::new();
        let mut l_fe = InterferenceList::new();
        let mut l_fef = InterferenceList::new();
        let mut l_e = InterferenceList::new();
        fun_select_trasha_interference(&mut lw, ShapeType::Face, &mut l_f);
        fun_select_gk_interference(&mut l_f, Kind::Edge, &mut l_fe);
        fun_select_sk_interference(&mut l_fe, Kind::Face, &mut l_fef);
        fun_select_trasha_interference(&mut lw, ShapeType::Edge, &mut l_e);

        li.clear();
        li.extend(l_f);
        li.extend(l_fe);
        li.extend(l_fef);
        li.extend(l_e);
        *self.hds.change_ds().change_shape_interferences(&s) = li;
    }

    /// `ProcessCurveInterferences()` (`cxx:83-93`): every curve of the DS.
    pub fn process_curve_interferences(&mut self) {
        let ncurve = self.hds.ds().nb_curves();
        for i in 1..=ncurve {
            self.process_curve_interferences_at(i);
        }
    }

    /// `ProcessCurveInterferences(const int CIX)`
    /// (`FilterCurveInterferences.cxx:22-65`): drop the `VERTEX` interferences
    /// that repeat a later one with the same orientation.
    pub fn process_curve_interferences_at(&mut self, cix: i32) {
        let mut li: InterferenceList = self.hds.ds().curve_interferences(cix).to_vec();
        let mut i = 0;
        while i < li.len() {
            let (g1, gt1, o1) = {
                let h = li[i].read().expect("interference lock");
                let o = h.transition().orientation(State::In, ShapeType::Face);
                (h.geometry(), h.geometry_type(), o)
            };
            if gt1 == Kind::Vertex {
                let mut j = i + 1;
                while j < li.len() {
                    let (g2, gt2, o2) = {
                        let h = li[j].read().expect("interference lock");
                        let o = h.transition().orientation(State::In, ShapeType::Face);
                        (h.geometry(), h.geometry_type(), o)
                    };
                    // xpu140898 : USA60111 : CPI(FORWARD,v10,C1) + CPIREV(REVERSED,v10,C1)
                    // do NOT delete CPIREV.
                    let remove = gt2 == gt1 && g2 == g1 && o1 == o2;
                    if remove {
                        li.remove(j);
                    } else {
                        j += 1;
                    }
                }
            }
            i += 1;
        }
        *self.hds.change_ds().change_curve_interferences(cix) = li;
    }
}

/// `FUN_FilterEdge` (`FilterEdgeInterferences.cxx:44-59`).
fn fun_filter_edge(li: &mut InterferenceList, hds: &HDataStructure, six: i32) {
    edge_tools::fun_order_ff_same_domain(li, hds, six);
    edge_tools::fun_order_state_trans(li, hds, six);
    edge_tools::fun_unkeep_e_interferences(li, hds.ds(), six);
    fun_reducedoublons(li, hds.ds(), six);
}
