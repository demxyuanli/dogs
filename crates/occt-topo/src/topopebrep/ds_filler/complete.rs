//! `TopOpeBRep_DSFiller` complete pass.
//!
//! Source: `TKBool/TopOpeBRep/TopOpeBRep_DSFiller.cxx:721-888`,
//! `TopOpeBRep_DSFiller.cxx:1001-1024` (OCCT 8.0.0).
//!
//! `Complete` (`cxx:723-738`) runs five steps in order: `GapFiller`,
//! `CompleteDS`, `Filter`, `Reducer`, `RemoveUnsharedGeometry`, `Checker`.
//! `Insert1d` (`cxx:1129-1130`) calls `CompleteDS` alone and `Insert2d`
//! (`cxx:898`) calls `CompleteDS2d` alone.
//!
//! UNPORTED steps, with the OCCT source of each:
//!
//! * `GapFiller` (`cxx:742-750`) — `TopOpeBRepDS_GapFiller::Perform`
//!   (`TopOpeBRepDS_GapFiller.cxx:52-84`), run against the edge connexity of
//!   `FDSCNX_Prepare` (`cxx:127-220`).
//! * `CompleteDS` (`cxx:754-788`) — the twelve `FUN_ds_*` completion passes
//!   live in `TopOpeBRepDS_EXPORT.cxx` and are ported in
//!   [`crate::topopebrep::ds::complete_ds`], which [`DsFiller::complete_ds`]
//!   calls in the OCCT order.
//! * `Filter` (`cxx:792-802`) — `TopOpeBRepDS_Filter::ProcessEdgeInterferences`
//!   / `ProcessCurveInterferences` (`TopOpeBRepDS_Filter.cxx`). [`filter`]
//!   is a no-op until that class is translated.
//! * `Checker` (`cxx:882-888`) — commented out in OCCT itself, so there is
//!   nothing to call.

use crate::abs::ShapeType;
use crate::topopebrep::ds::export::fun_ds_sdm;
use crate::topopebrep::ds::filter::Filter;
use crate::topopebrep::ds::gap_filler::GapFiller;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::reducer::Reducer;

use super::ds_tools::{brep_unfill_same_domain, fun_share_nog};
use super::DsFiller;

impl DsFiller {
    /// `Complete(HDS)` (`cxx:723-738`).
    pub fn complete(&mut self, hds: &mut HDataStructure) {
        self.gap_filler(hds);
        self.complete_ds(hds);
        self.filter(hds);
        self.reducer(hds);
        self.remove_unshared_geometry(hds);
        self.checker(hds);
    }

    /// `GapFiller(HDS)` (`cxx:742-750`):
    /// `TopOpeBRepDS_GapFiller(myHDS).Perform()`.
    pub fn gap_filler(&mut self, hds: &mut HDataStructure) {
        let mut gf = {
            let imm: &HDataStructure = hds;
            GapFiller::new(imm)
        };
        gf.perform(hds, &self.connex);
    }

    /// `CompleteDS(HDS)` (`cxx:754-788`).
    pub fn complete_ds(&mut self, hds: &mut HDataStructure) {
        // `cxx:759-765`: both intersected shapes must be present.
        let s1 = self.shape_intersector.shape(1);
        let s2 = self.shape_intersector.shape(2);
        if s1.is_none() || s2.is_none() {
            return;
        }

        // Each line below is one `FUN_ds_*` call of `cxx:770-784`, in order;
        // `FUN_ds_completeforSE7` is commented out at `cxx:780` and
        // `FUN_ds_unkeepEVIonGb1` (`cxx:769`) is commented out in OCCT too.
        crate::topopebrep::ds::complete_ds::run_complete_ds(hds, &self.connex);
    }

    /// `Filter(HDS)` (`cxx:792-802`).
    ///
    /// OCCT builds `TopOpeBRepDS_Filter F(HDS)`, i.e. the constructor whose
    /// `pClassif` defaults to null (`Filter.hxx`), and runs the edge pass then
    /// the curve pass.
    pub fn filter(&mut self, hds: &mut HDataStructure) {
        let mut f = Filter::new(hds, None);
        f.process_edge_interferences(&self.connex);
        f.process_curve_interferences();
    }

    /// `Reducer(HDS)` (`cxx:806-813`).
    pub fn reducer(&mut self, hds: &mut HDataStructure) {
        let mut r = Reducer::new(hds);
        // `ProcessEdgeInterferences` raises `Standard_ProgramError` on a
        // missing edge; the `Result` is dropped on purpose here because
        // `Complete` is `void` in OCCT (`cxx:806-813`).
        let _ = r.process_edge_interferences();
    }

    /// `RemoveUnsharedGeometry(HDS)` (`cxx:817-878`).
    pub fn remove_unshared_geometry(&mut self, hds: &mut HDataStructure) {
        let nbs = hds.ds().nb_shapes();
        for i in 1..=nbs {
            // `cxx:836-842`
            let Some(s) = hds.ds().shape(i, true) else {
                continue;
            };
            if s.shape_type() != ShapeType::Face {
                continue;
            }
            // `cxx:842-846`
            if !hds.has_same_domain(&s, true) {
                continue;
            }
            // `cxx:847-851`
            if hds.ds().ancestor_rank(i) != 1 {
                continue;
            }

            // `cxx:852-871`
            let l_s_sd = hds.ds().shape_same_domain(&s).to_vec();
            for s_sd in l_s_sd {
                // `cxx:858-862`
                let rk_s_sd = hds.ds().shape_index(&s_sd, true);
                if rk_s_sd > 0 && hds.ds().ancestor_rank(rk_s_sd) == 1 {
                    continue;
                }

                let mut unfill = fun_share_nog(hds, &s, &s_sd);
                // `cxx:863`: `unfill && FUN_ds_sdm(BDS,S,Ssd) && FUN_ds_sdm(BDS,Ssd,S)`
                unfill = unfill
                    && fun_ds_sdm(hds.ds(), &s, &s_sd)
                    && fun_ds_sdm(hds.ds(), &s_sd, &s);
                if unfill {
                    let psc = self.ensure_p_shape_classifier();
                    let mut sc = psc.borrow_mut();
                    let _ = brep_unfill_same_domain(&s, &s_sd, hds, &mut sc);
                }
            }
        }
    }

    /// `Checker(HDS)` (`cxx:882-888`): the body is commented out in OCCT
    /// (`TopOpeBRepDS_Checker` is NYI), so the port has nothing to run.
    pub fn checker(&mut self, _hds: &mut HDataStructure) {}

    /// `CompleteDS2d(HDS)` (`cxx:1001-1024`).
    pub fn complete_ds_2d(&mut self, hds: &mut HDataStructure) {
        // `cxx:1004-1011`
        let s1 = self.shape_intersector2d.shape(1);
        let s2 = self.shape_intersector2d.shape(2);
        if s1.is_none() || s2.is_none() {
            return;
        }
        let (Some(s1), Some(s2)) = (s1, s2) else {
            return;
        };
        hds.add_ancestors(&s1);
        hds.add_ancestors(&s2);

        // `cxx:1013-1017`: Filter then Reducer.
        let mut f = Filter::new(hds, None);
        f.process_edge_interferences(&self.connex);
        f.process_curve_interferences();
        let mut r = Reducer::new(hds);
        let _ = r.process_edge_interferences();
        // `cxx:1018`: `TopOpeBRepDS_Checker` is NYI in OCCT.
    }
}
