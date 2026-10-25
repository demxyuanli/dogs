//! `TopOpeBRepBuild_GridSS.cxx` — the face half of the solid/solid grid.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GridSS.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GSplitFaceSFS(FOR, LSclass, G1, SFS)` | `832-1028` |
//! | `GMergeFaceSFS(FOR, G1, SFS)` | `1031-1182` |
//!
//! `GMergeFaceSFS` is documented by OCCT as unused (`cxx:1029`); it is ported
//! for completeness because `Builder::GMergeFaceSFS` is part of the public
//! surface of the class.
//!
//! UNPORTED: the `OCCT_DEBUG` blocks (`cxx:844-858`, `:891-899`, `:971-984`,
//! `:1008-1019`, `:1033-1044`, `:1080-1089`, `:1118-1125`, `:1158-1170`,
//! `:1175-1181`) — trace and dump hooks only.

use crate::abs::{Orientation, State};
use crate::shape::TopoShape;
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::g_tool;
use crate::topopebrep::build::shell_face_set::ShellFaceSet;
use crate::topopebrep::ds::kind::Config;

use super::grid_ss_globals as gss;
use super::Builder;

impl Builder {
    /// `GSplitFaceSFS(FOR, LSclass, G1, SFS)` (`cxx:832-1028`).
    pub fn g_split_face_sfs(
        &mut self,
        for_: &TopoShape,
        ls_class: &[TopoShape],
        g1: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        let (tb1, tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();
        let ori_f = for_.orientation();
        let mut new_ori_f = Builder::orient(ori_f, rev_ori1);

        let issplit = self.is_split(for_, tb1);
        if issplit {
            // `cxx:869-871`: the split faces all carry the same topological
            // orientation; `FUNBUILD_ORIENTLOFS` re-orients them from the edge
            // ancestors when the grid asked for it.
            if gss::revownsplfacori() {
                let lsf = self.change_split(for_, tb1);
                gss::build_orient_lofs(tb1, tb2, lsf);
            }

            for new_f in self.splits(for_, tb1).to_vec() {
                if gss::split_anc_defined() {
                    let hasoridef = gss::split_anc_is_bound(&new_f);
                    let ope_fus = self.opefus();
                    let opec12 = self.opec12();
                    let opec21 = self.opec21();
                    let ope_cut = opec12 || opec21;
                    let ope_com = self.opecom();

                    if hasoridef {
                        let i_anc = gss::split_anc_find(&new_f).unwrap_or(0);

                        let (rk_anc, o_anc, iref_anc, c_anc, oref_anc) =
                            match self.data_structure() {
                                Some(hds) => {
                                    let iref_anc = hds.ds().same_domain_ref(i_anc);
                                    (
                                        hds.ds().ancestor_rank(i_anc),
                                        hds.shape(i_anc, false)
                                            .map(|s| s.orientation())
                                            .unwrap_or(Orientation::Forward),
                                        iref_anc,
                                        hds.ds().same_domain_ori(i_anc),
                                        hds.shape(iref_anc, false)
                                            .map(|s| s.orientation())
                                            .unwrap_or(Orientation::Forward),
                                    )
                                }
                                None => (
                                    0,
                                    Orientation::Forward,
                                    0,
                                    Config::UnshGeometry,
                                    Orientation::Forward,
                                ),
                            };

                        if ope_com {
                            // `cxx:900-904`: `orifspIN = orifanc`.
                            new_ori_f = o_anc;
                        } else if ope_cut {
                            // `cxx:905-928`.
                            let mut tb_anc = State::Unknown;
                            if opec12 {
                                tb_anc = if rk_anc == 1 { State::Out } else { State::In };
                            }
                            if opec21 {
                                tb_anc = if rk_anc == 2 { State::Out } else { State::In };
                            }
                            if tb_anc == State::Out {
                                new_ori_f = o_anc;
                            } else {
                                new_ori_f = gss::complement(o_anc);
                            }
                        } else if ope_fus {
                            // `cxx:929-932`.
                            new_ori_f = o_anc;
                        }

                        // `cxx:934-959`: the new face is built on the geometry
                        // of the reference ancestor.
                        let mut reverse = false;
                        if iref_anc != i_anc {
                            let mut samegeom = c_anc == Config::SameOriented;
                            if c_anc == Config::DiffOriented {
                                samegeom = false;
                            }
                            if o_anc != oref_anc {
                                samegeom = !samegeom;
                            }
                            reverse = !samegeom;
                        }
                        if reverse {
                            new_ori_f = gss::complement(new_ori_f);
                        }
                    }
                }

                let mut new_f = new_f;
                new_f.set_orientation(new_ori_f);

                // `cxx:966-984`.
                if gss::lfrtoprocess() {
                    gss::lfr1_append(&new_f);
                } else {
                    sfs.add_start_element(&new_f);
                }
            }
        } else {
            // `cxx:990-1002`: `FOR` has no split result on `TB1`, keep it when
            // it is located `TB1` against `LSclass`.
            let mut add = true;

            let (hs, hg) = match self.data_structure() {
                Some(hds) => (hds.has_shape(for_, true), hds.has_geometry(for_)),
                None => (false, false),
            };
            let testkeep = hs && !hg;

            if testkeep {
                add = self.g_keep_shape(for_, ls_class, tb1);
            }
            if add {
                let mut f = for_.clone();
                f.set_orientation(new_ori_f);
                sfs.add_element(&f);
            }
        }
    }

    /// `GMergeFaceSFS(FOR, G1, SFS)` (`cxx:1031-1182`).
    pub fn g_merge_face_sfs(&mut self, for_: &TopoShape, g1: &GTopo, sfs: &mut ShellFaceSet) {
        let tomerge = self.to_merge(for_);
        if !tomerge {
            return;
        }

        // `cxx:1056-1057`: `StatesON` is called for its `ON`-state guard.
        let _ = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();

        let ori_f = for_.orientation();
        let new_ori_f = Builder::orient(ori_f, rev_ori1);

        let mut ff = for_.clone();
        ff.set_orientation(Orientation::Forward);

        let makecomsam = Builder::take_common_of_same(g1);
        let makecomdif = Builder::take_common_of_diff(g1);
        if !makecomsam && !makecomdif {
            return;
        }

        // `cxx:1072-1078`: same-domain families of the two parent shapes.
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

        // `cxx:1092-1116`.
        let mut performcom = false;
        let mut lf1: Vec<TopoShape> = Vec::new();
        let mut lf2: Vec<TopoShape> = Vec::new();
        if makecomsam {
            performcom = !lfso1.is_empty() && !lfso2.is_empty();
            if performcom {
                lf1 = lfso1.clone();
                lf2 = lfso2.clone();
            }
        } else if makecomdif {
            performcom = !lfso1.is_empty() && !lfdo2.is_empty();
            if performcom {
                lf1 = lfso1.clone();
                lf2 = lfdo2.clone();
            }
        }

        if !performcom {
            return;
        }

        // `cxx:1127-1140`: the merge grid on faces.
        let mut gf = g_tool::GTool::g_com_unsh(crate::abs::ShapeType::Face, crate::abs::ShapeType::Face);
        if makecomsam {
            gf.change_config(Config::SameOriented, Config::SameOriented);
        } else if makecomdif {
            gf.change_config(Config::SameOriented, Config::DiffOriented);
        }

        self.g_merge_faces(&lf1, &lf2, &gf);

        // `cxx:1142-1172`: keep the merge result only from the reference face
        // of the same-domain family.
        let addmerge = match self.data_structure() {
            Some(hds) => {
                let i_fref = hds.same_domain_reference(for_);
                match hds.shape(i_fref, false) {
                    Some(fref) => crate::topo_tools::is_same(for_, &fref),
                    None => false,
                }
            }
            None => false,
        };

        if addmerge {
            for new_f in self.merged(for_, State::In).to_vec() {
                let mut new_f = new_f;
                new_f.set_orientation(new_ori_f);
                sfs.add_start_element(&new_f);
            }
        }
    }
}
