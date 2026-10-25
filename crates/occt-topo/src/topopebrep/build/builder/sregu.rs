//! `TopOpeBRepBuild_Builder::RegularizeSolids` / `RegularizeSolid`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_SREGU.cxx:39-333`
//! (OCCT 8.0.0).
//!
//! `RegularizeSolid` regularizes the shells of one new solid
//! (`TopOpeBRepTool::RegularizeShells`), rebuilds solids out of the new shells
//! (`TopOpeBRepBuild_ShellToSolid`), and then folds the second-level face
//! splits back into the builder's split map. `RegularizeSolids` runs that over
//! a list of solids and then re-expands, for the same-domain solids of `SO`,
//! every split face that the previous step recorded in `myMemoSplit`.
//!
//! Release-visible behaviour only: the `OCCT_DEBUG` context switches
//! `GetcontextNOREGUSO` / `GetcontextREGUXPU` (`cxx:143-146`) are debug-build
//! only, so `toregu` and `usestos` are the constant `true` here.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};
use crate::topopebrep::build::ShellToSolid;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;
use crate::topopebrep::tool::regularize_shells;

use super::Builder;

impl Builder {
    /// `RegularizeSolids(SO, lnewSolid, LOSO)` (`cxx:41-126`).
    pub fn regularize_solids(
        &mut self,
        so: &TopoShape,
        lnew_solid: &[TopoShape],
        loso: &mut Vec<TopoShape>,
    ) {
        loso.clear();
        self.my_memo_split.clear();

        // `cxx:48-58`.
        for new_solid in lnew_solid {
            let mut new_solid_loso: Vec<TopoShape> = Vec::new();
            self.regularize_solid(so, new_solid, &mut new_solid_loso);
            loso.extend(new_solid_loso);
        }

        // `cxx:63-67`.
        if self.my_memo_split.extent() == 0 {
            return;
        }

        // `cxx:69-73`: lsosdSO = the same-domain solids of SO.
        let mut l1: Vec<TopoShape> = Vec::new();
        let mut l2: Vec<TopoShape> = Vec::new();
        self.find_same_domain_from(so, &mut l1, &mut l2);
        l1.extend(l2);

        // `cxx:75-125`.
        for sosd_so in l1 {
            for f in explorer_shapes(&sosd_so, ShapeType::Face) {
                let rankf = self.shape_rank(&f);
                let staf = if rankf == 1 {
                    self.my_state1
                } else {
                    self.my_state2
                };

                // `cxx:87-90`.
                if !self.is_split(&f, staf) {
                    continue;
                }

                // `cxx:93-95`: snapshot `lspf` (`ChangeSplit` aliases the map).
                let lspf = self.splits(&f, staf).to_vec();
                let mut new_lspf: Vec<TopoShape> = Vec::new();
                for fsp in lspf {
                    // `cxx:100-108`.
                    if !self.my_memo_split.contains(&fsp) {
                        new_lspf.push(fsp);
                    } else {
                        let lspfsp = self.splits(&fsp, staf).to_vec();
                        new_lspf.extend(lspfsp);
                    }
                }

                // `cxx:110-111`: `lspf.Clear(); GCopyList(newlspf, lspf);`
                let lspf = self.change_split(&f, staf);
                lspf.clear();
                lspf.extend(new_lspf);
            }
        }
    }

    /// `RegularizeSolid(SS, anewSolid, LOSO)` (`cxx:129-333`).
    pub fn regularize_solid(
        &mut self,
        ss: &TopoShape,
        anew_solid: &TopoShape,
        loso: &mut Vec<TopoShape>,
    ) {
        loso.clear();
        let new_solid = anew_solid.clone();

        // `cxx:136-137`, `cxx:150-155`: `toregu` is always true outside
        // `OCCT_DEBUG`.
        let toregu = true;
        if !toregu {
            loso.push(new_solid);
            return;
        }

        // `cxx:157-161`: OldShells --> NewShells, then the ancestor-face splits.
        let mut osns: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut f_splits: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let rw = regularize_shells(&new_solid, &mut osns, &mut f_splits);
        // OCCT passes `myFSplits` itself, so a partial answer is kept.
        self.my_f_splits = f_splits;

        // `cxx:164-168`.
        if !rw {
            loso.push(new_solid);
            return;
        }

        // `cxx:170-211`: `usestos` is always true outside `OCCT_DEBUG`.
        let mut stos = ShellToSolid::new();
        // `cxx:176-207`: walk the map (OCCT's `DataMap` order is a bucket
        // order; the shell order of `newSolid` is the deterministic choice).
        for sh in explorer_shapes(&new_solid, ShapeType::Shell) {
            let Some(lns) = osns.get(&shape_map_key(&sh)) else {
                continue;
            };
            for w in lns.clone() {
                stos.add_shell(&w);
            }
        }
        let mut new_solids: Vec<TopoShape> = Vec::new();
        stos.make_solids(&new_solid, &mut new_solids);
        let rf = !new_solids.is_empty();

        // `cxx:214-220`.
        if !rf {
            loso.push(new_solid);
            return;
        }

        // `cxx:230-236`.
        loso.extend(new_solids);

        // `cxx:241-248`: mfns = the faces of newSolid.
        let mut mfns = IndexedShapeSet::new();
        for f in explorer_shapes(&new_solid, ShapeType::Face) {
            mfns.add(&f);
        }

        // `cxx:250-254`: lssdSS = the same-domain solids of SS.
        let mut l1: Vec<TopoShape> = Vec::new();
        let mut l2: Vec<TopoShape> = Vec::new();
        self.find_same_domain_from(ss, &mut l1, &mut l2);
        let mut lssd_ss = l1;
        lssd_ss.extend(l2);

        // `cxx:257-331`.
        for ssd_ss in lssd_ss {
            let rankssd_ss = self.shape_rank(&ssd_ss);
            let stassd_ss = if rankssd_ss == 1 {
                self.my_state1
            } else {
                self.my_state2
            };

            for ssd_ssf in explorer_shapes(&ssd_ss, ShapeType::Face) {
                let stassd_ssf = stassd_ss;
                // `cxx:289`: snapshot, `ChangeSplit` aliases the map.
                let lspssd_ssf = self.change_split(&ssd_ssf, stassd_ssf).clone();

                for fspssd_ssf in lspssd_ssf {
                    // `cxx:301-305`.
                    if !mfns.contains(&fspssd_ssf) {
                        continue;
                    }

                    // `cxx:310-319`.
                    let Some(lrfsplit) = self
                        .my_f_splits
                        .get(&shape_map_key(&fspssd_ssf))
                        .cloned()
                    else {
                        continue;
                    };

                    // `cxx:322`: record that fspssdSSf was re-split.
                    self.my_memo_split.add(&fspssd_ssf);

                    // `cxx:325-326`: `ChangeSplit(fspssdSSf, stassdSSf)` then
                    // `GCopyList(lrfsplit, lsp)`.
                    let lsp = self.change_split(&fspssd_ssf, stassd_ssf);
                    lsp.extend(lrfsplit);
                }
            }
        }
    }
}
