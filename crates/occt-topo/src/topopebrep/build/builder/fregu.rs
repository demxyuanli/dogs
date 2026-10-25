//! `TopOpeBRepBuild_Builder::RegularizeFaces` / `RegularizeFace`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_FREGU.cxx`
//! (OCCT 8.0.0, 373 lines): `RegularizeFaces` (`:39-131`) and `RegularizeFace`
//! (`:160-372`).
//!
//! UNPORTED:
//! - `FUN_setInternal` (`cxx:134-158`) is inside a `/* ... */` comment, and so
//!   is its call at `cxx:184`.
//! - the `OCCT_DEBUG` branches (`GetcontextNOREGUFA` / `GetcontextREGUXPU` /
//!   `GtraceSPS` / `debregufa`): they only pick `toregu` / `usewtof`, which are
//!   both `true` without a debug context (`cxx:165-166`).

use crate::abs::{ShapeType, State};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, shape_map_key};
use crate::topopebrep::ds::data_structure::IndexedShapeSet;
use crate::topopebrep::tool::regularize_w::{
    regularize_face as tool_regularize_face, regularize_wires,
};

use crate::topopebrep::build::wire_to_face::WireToFace;
use super::Builder;

impl Builder {
    /// `RegularizeFaces(FF, lnewFace, LOF)` (`cxx:39-131`).
    pub fn regularize_faces(
        &mut self,
        ff: &TopoShape,
        lnew_face: &[TopoShape],
        lof: &mut Vec<TopoShape>,
    ) {
        lof.clear();
        self.my_memo_split.clear();

        for new_face in lnew_face {
            let mut new_face_lof: Vec<TopoShape> = Vec::new();
            self.regularize_face(ff, new_face, &mut new_face_lof);
            lof.extend(new_face_lof);
        }

        let nr = self.my_memo_split.extent();
        if nr == 0 {
            return;
        }

        // lfsdFF = faces SameDomain de FF (`cxx:68-72`)
        let mut lfsd1: Vec<TopoShape> = vec![ff.clone()];
        let mut lfsd2: Vec<TopoShape> = Vec::new();
        self.find_same_domain(&mut lfsd1, &mut lfsd2);
        let mut lfsd: Vec<TopoShape> = Vec::new();
        lfsd.extend(lfsd1);
        lfsd.extend(lfsd2);

        for fsd_ff in lfsd {
            for e in explorer_shapes(&fsd_ff, ShapeType::Edge) {
                let ranke = self.shape_rank(&e);
                let staeope = if ranke == 1 {
                    self.my_state1
                } else {
                    self.my_state2
                };

                for iiista in 1..=2 {
                    let stae = if iiista == 2 { State::On } else { staeope };

                    if !self.is_split(&e, stae) {
                        continue;
                    }

                    // `cxx:107-127`: replace every memorised split by its own
                    // splits, keep the others.
                    let lspe_snapshot = self.splits(&e, stae).to_vec();
                    let mut new_lspe: Vec<TopoShape> = Vec::new();
                    for esp in &lspe_snapshot {
                        let espmemo = self.my_memo_split.contains(esp);
                        if !espmemo {
                            new_lspe.push(esp.clone());
                        } else {
                            let lspesp = self.splits(esp, stae);
                            Self::copy_list(lspesp, &mut new_lspe);
                        }
                    }
                    {
                        let lspe = self.change_split(&e, stae);
                        lspe.clear();
                        Self::copy_list(&new_lspe, lspe);
                    }
                }
            }
        }
    }

    /// `RegularizeFace(FF, anewFace, LOF)` (`cxx:160-372`).
    pub fn regularize_face(
        &mut self,
        ff: &TopoShape,
        anew_face: &TopoShape,
        lof: &mut Vec<TopoShape>,
    ) {
        lof.clear();
        let new_face = anew_face.clone();

        // `cxx:168-169`: both are only turned off by the debug contexts.
        let toregu = true;
        let usewtof = true;

        if !toregu {
            lof.push(new_face);
            return;
        }

        // `cxx:193-198`
        let mut ownw: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
        self.my_e_splits.clear();

        let rw = regularize_wires(&new_face, &mut ownw, &mut self.my_e_splits);
        if !rw {
            lof.push(new_face);
            return;
        }

        let mut newfaces: Vec<TopoShape> = Vec::new();
        let rf;
        if usewtof {
            // `cxx:206-235`
            let mut wtof = WireToFace::new();
            for (ow, lw) in &ownw {
                // if lnw is empty, then ow is kept unchanged
                if lw.is_empty() {
                    wtof.add_wire(ow);
                }
                for w in lw {
                    wtof.add_wire(w);
                }
            }
            // `cxx:232`: OCCT runs this on a throwaway `Builder` created by
            // `WireToFace`; the port passes the builder it is called on.
            wtof.make_faces(&new_face, self, &mut newfaces);
            rf = !newfaces.is_empty();
        } else {
            rf = tool_regularize_face(&new_face, &ownw, &mut newfaces);
        }

        if !rf {
            lof.push(new_face);
            return;
        }

        // LOF = new regularized faces of newFace (`cxx:258-262`)
        lof.extend(newfaces);

        // update of the split edges (`cxx:264-370`)
        let mut menf = IndexedShapeSet::new();
        for e in explorer_shapes(&new_face, ShapeType::Edge) {
            menf.add(&e);
        }

        let mut lfsd1: Vec<TopoShape> = vec![ff.clone()];
        let mut lfsd2: Vec<TopoShape> = Vec::new();
        self.find_same_domain(&mut lfsd1, &mut lfsd2);
        let mut lfsd: Vec<TopoShape> = Vec::new();
        lfsd.extend(lfsd1);
        lfsd.extend(lfsd2);

        for fsd_ff in lfsd {
            let rankfsd_ff = self.shape_rank(&fsd_ff);
            let stafsd_ff = if rankfsd_ff == 1 {
                self.my_state1
            } else {
                self.my_state2
            };

            for fsd_ffe in explorer_shapes(&fsd_ff, ShapeType::Edge) {
                for iiista in 1..=2 {
                    let stafsd_ffe = if iiista == 2 { State::On } else { stafsd_ff };

                    let lsp_snapshot = self.splits(&fsd_ffe, stafsd_ffe).to_vec();
                    for espfsd_ffe in lsp_snapshot {
                        if !menf.contains(&espfsd_ffe) {
                            continue;
                        }

                        // espfsdFFe is a split edge of fsdFFe and belongs to
                        // newFace; RegularizeWires has re-split it.
                        let resplitloc = self
                            .my_e_splits
                            .get(&shape_map_key(&espfsd_ffe))
                            .cloned();
                        if let Some(lresplit) = resplitloc {
                            self.my_memo_split.add(&espfsd_ffe);
                            let lsp = self.change_split(&espfsd_ffe, stafsd_ffe);
                            lsp.clear();
                            Self::copy_list(&lresplit, lsp);
                        }
                    }
                }
            }
        }
    }
}
