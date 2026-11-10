//! `TopOpeBRepBuild_ffsfs.cxx` — filling a `ShellFaceSet` with the faces of a
//! same-domain family.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_ffsfs.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `FUN_setmotherope` / `FUN_unsetmotherope` / `FUN_ismotheropedef` / `FUN_motherope` | `55-77` (no reader; see `grid_ss_globals`) |
//! | `FUN_getAncestorFsp(B, SC, LF, fsp, p3ddef, p3d)` | `111-163` |
//! | `FUN_getAncestorFsp(B, SC, LF1, LF2, fsp)` | `165-231` |
//! | `FUN_getAncestorFsp(B, SC, LF1, LF2, FOR, SplitAnc)` | `233-282` |
//! | `GFillFaceSFS(FOR, LSO2, Gin, SFS)` | `293-691` |
//!
//! The `GLOBAL_*` state these three share with `GridSS.cxx` lives in
//! [`super::grid_ss_globals`].
//!
//! UNPORTED: the `OCCT_DEBUG` blocks (`cxx:301-307`, `:313-322`, `:348-353`,
//! `:406-419`, `:476-491`, `:564-579`, `:606-612`, `:619-621`, `:628-631`,
//! `:652-664`, `:683-689`) — trace and dump hooks only.

use crate::abs::{Orientation, State};
use crate::brep_class3d::find_a_point_in_the_face;
use crate::fclass2d::{FaceState, FClass2d};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::build::builder::grid_ff_edges::set_global_classify_split_edge;
use crate::topopebrep::build::fuse_face::FuseFace;
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::shell_face_set::ShellFaceSet;
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::tool::fun_tool_proj_pon_f;

use super::grid_ss_globals as gss;
use super::Builder;

impl Builder {
    /// `GFillFaceSFS(FOR, LSO2, Gin, SFS)` (`cxx:293-691`).
    pub fn g_fill_face_sfs(
        &mut self,
        for_: &TopoShape,
        lso2: &[TopoShape],
        gin: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        let (tb1, tb2) = gin.states_on();

        let tosplit = self.g_to_split(for_, tb1);
        let tomerge = self.to_merge(for_);
        let rk_for = self.ancestor_rank_of(for_);

        let mut ff = for_.clone();
        ff.set_orientation(Orientation::Forward);
        let hsd = self.has_same_domain(for_);
        gss::set_lfrtoprocess(false);

        if tosplit && tomerge {
            // `cxx:335-339`: merge the same-domain faces only from the
            // reference face of the family.
            let fref = match self.data_structure() {
                Some(hds) => {
                    let iref = hds.same_domain_reference(for_);
                    hds.shape(iref, false)
                }
                None => None,
            };
            let isref = fref.as_ref().map(|f| is_same(for_, f)).unwrap_or(false);

            // `cxx:341-343`: `xpu280199 PRO16958`, `GToMerge(fref)` may be
            // false while the face is still merged.
            let makemergeref = fref.as_ref().map(|f| self.to_merge(f)).unwrap_or(false);
            let mut makemerge = isref;
            makemerge = makemerge || (!makemergeref && rk_for == 1);

            if makemerge {
                // `cxx:357-359`: `performfufa` is true outside `OCCT_DEBUG`.
                let performfufa = true;

                // `cxx:360-370`.
                gss::set_lfrtoprocess(true);
                if gss::lfrtoprocess() {
                    gss::lfr1_clear();
                }

                // `cxx:372-381`.
                gss::split_anc_prepare();

                // `cxx:384-391`.
                let (mut lfso, mut lfdo) = (Vec::new(), Vec::new());
                self.find_same_domain_same_orientation_from(&ff, &mut lfso, &mut lfdo);
                let rank_f = self.shape_rank(&ff);
                let rank_x = if rank_f != 0 {
                    if rank_f == 1 {
                        2
                    } else {
                        1
                    }
                } else {
                    0
                };
                let mut lfso1: Vec<TopoShape> = Vec::new();
                let mut lfdo1: Vec<TopoShape> = Vec::new();
                let mut lfso2: Vec<TopoShape> = Vec::new();
                let mut lfdo2: Vec<TopoShape> = Vec::new();
                self.find_same_rank(&lfso, rank_f, &mut lfso1);
                self.find_same_rank(&lfdo, rank_f, &mut lfdo1);
                self.find_same_rank(&lfso, rank_x, &mut lfso2);
                self.find_same_rank(&lfdo, rank_x, &mut lfdo2);

                // `cxx:401-408`: WES of the faces with the same topological
                // orientation.
                let mut lf1 = lfso1.clone();
                let mut lf2 = lfso2.clone();
                let mut dodo = !lf1.is_empty() && !lf2.is_empty();

                // `cxx:410`: `FUN_unsetmotherope();` — the port passes the
                // `GTopo` into `g_fill_faces_wes_make_faces` directly, so the
                // static keeps no `GTopo` reader.

                let mut gm = gin.clone();
                gm.change_config(Config::SameOriented, Config::SameOriented);
                if dodo {
                    set_global_classify_split_edge(true);
                    self.g_fill_faces_wes_make_faces(&lf1, &lf2, &gm);
                    set_global_classify_split_edge(false);

                    gss::set_revownsplfacori(true);
                    gss::build_ancestor_rank_prepare(
                        self,
                        &lf1,
                        &lf2,
                        Config::SameOriented,
                        Config::SameOriented,
                    );
                    if hsd {
                        self.fun_get_ancestor_fsp_split_anc(&lf1, &lf2, for_);
                    }

                    // `cxx:441-451`.
                    gss::lfr1_clear();
                    self.g_split_face_sfs(for_, lso2, &gm, sfs);
                    gss::set_lfrtoprocess(false);
                    gss::set_revownsplfacori(false);
                }

                // `cxx:457-471`: WES `FOR` + the faces of the opposite
                // topological orientation.
                lf1 = lfso1.clone();
                lf2 = lfdo2.clone();
                dodo = !lf1.is_empty() && !lf2.is_empty();

                gm = gin.clone();
                let tb = tb2;
                let ntb = State::In;
                gm.change_value_state(tb, State::On, false);
                gm.change_value_state(ntb, State::On, true);
                gm.change_config(Config::SameOriented, Config::DiffOriented);
                // `cxx:481`: `FUN_setmotherope(GM);` — see the note above.
                if dodo {
                    set_global_classify_split_edge(true);
                    self.g_fill_faces_wes_make_faces(&lf1, &lf2, &gm);
                    set_global_classify_split_edge(false);

                    gss::set_revownsplfacori(true);
                    gss::build_ancestor_rank_prepare(
                        self,
                        &lf1,
                        &lf2,
                        Config::SameOriented,
                        Config::DiffOriented,
                    );
                    if hsd {
                        self.fun_get_ancestor_fsp_split_anc(&lf1, &lf2, for_);
                    }

                    // `cxx:512-545`: for a COMMON the split of the reference
                    // face is kept, the splits of the family members are reset.
                    if self.opecom() {
                        let issplit_in = self.is_split(for_, State::In);
                        if issplit_in {
                            let sp_for = self.splits(for_, State::In).to_vec();
                            for f in &lf1 {
                                if self.is_split(f, State::In) {
                                    self.change_split(f, State::In).clear();
                                }
                            }
                            for f in &lf2 {
                                if self.is_split(f, State::In) {
                                    self.change_split(f, State::In).clear();
                                }
                            }
                            let l = self.change_split(for_, State::In);
                            l.extend(sp_for);
                        }
                    }

                    self.g_split_face_sfs(for_, lso2, &gm, sfs);
                    gss::set_revownsplfacori(false);
                }

                // `cxx:551-603`: the FUSE grid.
                if !self.opecom() {
                    gm = gin.clone();
                    let tb = tb1;
                    let ntb = State::In;
                    gm.change_value_state(State::On, tb, false);
                    gm.change_value_state(State::On, ntb, true);
                    gm.change_config(Config::SameOriented, Config::DiffOriented);
                    // `cxx:556`: `FUN_setmotherope(GM);` — see the note above.
                    if dodo {
                        set_global_classify_split_edge(true);
                        self.g_fill_faces_wes_make_faces(&lf1, &lf2, &gm);
                        set_global_classify_split_edge(false);

                        gss::set_revownsplfacori(true);
                        gss::build_ancestor_rank_prepare(
                            self,
                            &lf1,
                            &lf2,
                            Config::SameOriented,
                            Config::DiffOriented,
                        );
                        if hsd {
                            self.fun_get_ancestor_fsp_split_anc(&lf1, &lf2, for_);
                        }
                        self.g_split_face_sfs(for_, lso2, &gm, sfs);
                        gss::set_revownsplfacori(false);
                    }
                }

                // `cxx:608`: `FUN_unsetmotherope();` — see the note above.

                // `cxx:613`: `GLOBAL_SplitAnc->Clear()`.
                gss::split_anc_clear();

                // `cxx:616`: `SFS.ChangeStartShapes().Extent();` — a read with
                // no effect on the set.

                // `cxx:618-667`: `TopOpeBRepBuild_FuseFace` on `GLOBAL_lfr1`.
                // `addinternal = 1` (`cxx:636`), i.e. `CXM == 1` ->
                // `myInternal = false`, and `ldum` is the empty `LIF`.
                if performfufa {
                    let mut fufa = FuseFace::new();
                    fufa.init(&[], &gss::lfr1_snapshot(), 1);
                    fufa.perform_face();
                    // `cxx:639-642`.
                    if !fufa.is_done() {
                        return;
                    }
                    // `cxx:650-667`: the fused faces are start elements.
                    for f in fufa.l_fuse_face().to_vec() {
                        sfs.add_start_element(&f);
                    }
                }
            }
        } else if tosplit && !tomerge {
            // `cxx:672-676`.
            self.g_split_face(for_, gin, lso2);
            self.g_split_face_sfs(for_, lso2, gin, sfs);
        } else if !tosplit && !tomerge {
            // `cxx:677-680`.
            self.g_split_face_sfs(for_, lso2, gin, sfs);
        }

        // `cxx:682`.
        self.my_edge_avoid.clear();
    }

    /// `FUN_getAncestorFsp(B, SC, LF, fsp, p3ddef, p3d)` (`cxx:111-163`): the
    /// index of the face of `LF` that holds the IN part of `fsp`.
    fn fun_get_ancestor_fsp_in_list(
        &mut self,
        lf: &[TopoShape],
        fsp: &TopoShape,
        p3ddef: &mut bool,
        p3d: &mut occt_core::gp::GpPnt,
    ) -> i32 {
        for f in lf {
            let st = self
                .my_shape_classifier
                .borrow_mut()
                .state_shape_shape(fsp, f, 1);
            let mut st = match st {
                Ok(st) => st,
                // `StateShapeShape` throws on a null reference; the port keeps
                // OCCT's `continue` on an unusable answer.
                Err(_) => continue,
            };
            if st == State::Unknown || st == State::Out {
                continue;
            }
            if st == State::On {
                if !*p3ddef {
                    let Some(face) = Face::wrap(fsp.clone()) else {
                        return 0;
                    };
                    if !find_a_point_in_the_face(&face, p3d) {
                        return 0;
                    }
                    *p3ddef = true;
                }
                let Some(face) = Face::wrap(f.clone()) else {
                    return 0;
                };
                let Some((pf, pl, dd)) = fun_tool_proj_pon_f(p3d, &face) else {
                    return 0;
                };
                // `cxx:147`: `tolf = BRep_Tool::Tolerance(f) * 1.e1`.
                let tolf = crate::brep_tool::BRepTool::face_tolerance(&face) * 1.0e1;
                if dd > tolf {
                    return 0;
                }
                // `cxx:151-152`: `TolClass = 1e-8` then `FClass2d.Perform(p2d)`.
                let st_class = match FClass2d::new(&face, 1e-8) {
                    Ok(cl) => cl.perform_tab_orien(occt_core::gp::GpPnt2d::new(pf, pl)),
                    Err(_) => FaceState::Unknown,
                };
                st = match st_class {
                    FaceState::In => State::In,
                    FaceState::Out => State::Out,
                    FaceState::On => State::On,
                    FaceState::Unknown => State::Unknown,
                };
            }
            if st == State::In {
                return self.shape_index_of(f);
            }
        }
        0
    }

    /// `FUN_getAncestorFsp(B, SC, LF1, LF2, fsp)` (`cxx:165-231`).
    fn fun_get_ancestor_fsp2(&mut self, lf1: &[TopoShape], lf2: &[TopoShape], fsp: &TopoShape) -> i32 {
        let (of1, of2) = gss::build_ancestor_rank_get(fsp);
        let rk_fsp = if of1 && !of2 {
            1
        } else if of2 && !of1 {
            2
        } else {
            0
        };
        let unk = rk_fsp == 0;

        let (rk_f1, rk_f2) = match self.data_structure() {
            Some(hds) => (
                lf1.first()
                    .map(|f| hds.ds().ancestor_rank_of(f))
                    .unwrap_or(0),
                lf2.first()
                    .map(|f| hds.ds().ancestor_rank_of(f))
                    .unwrap_or(0),
            ),
            None => (0, 0),
        };
        let mut p3ddef = false;
        let mut p3d = occt_core::gp::GpPnt::new(0.0, 0.0, 0.0);

        let ison1 = rk_f1 == rk_fsp;
        let ison2 = rk_f2 == rk_fsp;

        let mut ianc1 = 0;
        let mut ianc2 = 0;
        if ison1 || unk {
            ianc1 = self.fun_get_ancestor_fsp_in_list(lf1, fsp, &mut p3ddef, &mut p3d);
        }
        if ison1 {
            return ianc1;
        }
        if ison2 || unk {
            ianc2 = self.fun_get_ancestor_fsp_in_list(lf2, fsp, &mut p3ddef, &mut p3d);
        }
        if ison2 {
            return ianc2;
        }
        if ianc1 + ianc2 > 0 {
            if ianc1 == 0 {
                return ianc2;
            } else if ianc2 == 0 {
                return ianc1;
            }
            // `fsp` has two ancestor faces.
            return 0;
        }
        0
    }

    /// `FUN_getAncestorFsp(B, SC, LF1, LF2, FOR, SplitAnc)` (`cxx:233-282`).
    fn fun_get_ancestor_fsp_split_anc(
        &mut self,
        lf1: &[TopoShape],
        lf2: &[TopoShape],
        for_: &TopoShape,
    ) {
        if !gss::split_anc_defined() {
            return;
        }

        let issplit_in = self.is_split(for_, State::In);
        let issplit_out = self.is_split(for_, State::Out);
        if !issplit_in && !issplit_out {
            return;
        }

        let mut sp_for: Vec<TopoShape> = Vec::new();
        if issplit_in {
            sp_for.extend_from_slice(self.splits(for_, State::In));
        }
        if issplit_out {
            sp_for.extend_from_slice(self.splits(for_, State::Out));
        }

        for fsp in &sp_for {
            if gss::split_anc_is_bound(fsp) {
                continue;
            }
            let ianc = self.fun_get_ancestor_fsp2(lf1, lf2, fsp);
            if ianc != 0 {
                gss::split_anc_bind(fsp, ianc);
            }
        }
    }

    /// `BDS.Shape(F)` (`HDataStructure::Shape(S)`): the DS index of `s`; the
    /// port answers `-1` where OCCT would answer `0`.
    fn shape_index_of(&self, s: &TopoShape) -> i32 {
        match self.data_structure() {
            Some(hds) => hds.shape_index(s, false),
            None => 0,
        }
    }

    /// `BDS.AncestorRank(S)` for a shape.
    fn ancestor_rank_of(&self, s: &TopoShape) -> i32 {
        match self.data_structure() {
            Some(hds) => hds.ds().ancestor_rank_of(s),
            None => 0,
        }
    }

    /// `myDataStructure->HasSameDomain(S)` (`hxx:`, default `FindKeep` true).
    pub(crate) fn has_same_domain(&self, s: &TopoShape) -> bool {
        match self.data_structure() {
            Some(hds) => hds.has_same_domain(s, true),
            None => false,
        }
    }
}
