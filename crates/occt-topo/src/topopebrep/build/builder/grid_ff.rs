//! `TopOpeBRepBuild_GridFF` — the face/face half of the `GTopo` grid.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GridFF.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `FUN_computeLIFfaces2d(BU, F, E, pDS2d)` | `144-308` |
//! | `FUN_computeLIFfaces2d(BU, F, pDS2d)` | `312-327` |
//! | `GMergeFaces(LF1, LF2, G1)` | `335-463` |
//! | `FUN_validF1edge(F)` | `474-507` |
//! | `GFillFacesWESMakeFaces(LLF1, LF2, LSO, GM)` | `511-635` |
//! | `GFillFaceWES(FOR1, LFclass, G1, WES)` | `639-706` |
//! | `GFillWireWES(W, LSclass, G1, WES)` | `709-757` |
//! | `GFillEdgeWES(EOR, LSclass, G1, WES)` | `761-800` |
//!
//! `GLOBAL_DS2d` (`cxx:332`) is a process-wide `TopOpeBRepDS_PDataStructure`
//! that `GMergeFaces` re-`Init`s and then fills with the 2D interferences of
//! the merged face pair. `Perform2d` reads it while `GLOBAL_faces2d` is set, so
//! the port keeps it in the thread-local cell of `builder/on/perform_2d.rs`
//! rather than in a local (`T-358`).
//!
//! `GLOBAL_faces2d` (`cxx:88`) is a plain `bool` that OCCT sets around the
//! `GFillFacesWESK` pass; the port keeps the same member field.
//!
//! UNPORTED: the `OCCT_DEBUG` blocks (`cxx:353-363`, `:549-558`, `:614-617`,
//! `:649-653`, `:700-703`, `:717-735`, `:741-752`, `:767-795`) — trace and
//! dump only.

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::build::face_builder::FaceBuilder;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::data_structure::{DataStructure, IndexedShapeSet};
use crate::topopebrep::ds::export::{fds_data, fds_idata, fun_ds_has_fei, fun_ds_sdm};
use crate::topopebrep::ds::kind::{Config, Kind};
use crate::topopebrep::ds::process_interferences_tool::fun_reducedoublons;
use crate::topopebrep::tool::{fun_tool_find_ancestor, fun_tool_par_e};

use super::super::g_topo::GTopo;
use super::grid_ee::set_global_faces2d;
use super::on::{init_global_ds2d, with_global_ds2d_mut};
use super::Builder;
impl Builder {
    /// `GMergeFaces(LF1, LF2, G1)` (`cxx:335-463`).
    pub fn g_merge_faces(&mut self, lf1: &[TopoShape], lf2: &[TopoShape], g1: &GTopo) {
        // `cxx:339-343`
        if lf1.is_empty() {
            return;
        }

        // `cxx:345-348`: `GLOBAL_DS2d` is created on first use and `Init`ed on
        // every call; the port keeps it in the thread-local cell of
        // `builder/on/perform_2d.rs`, which is where `Perform2d` reads it.
        let (tb1, tb2) = g1.states_on();

        let Some(f1) = lf1.first().cloned() else {
            return;
        };

        // `cxx:369-404`: the 2D interferences of both face lists.
        init_global_ds2d();
        with_global_ds2d_mut(|ds2d| {
            // `cxx:369-380`
            for f in lf1 {
                ds2d.add_shape_with_ancestor(f, 1);
            }
            for f in lf2 {
                ds2d.add_shape_with_ancestor(f, 2);
            }

            // `cxx:382-394`
            for f in lf1 {
                self.fun_compute_lif_faces2d(f, ds2d);
            }
            for f in lf2 {
                self.fun_compute_lif_faces2d(f, ds2d);
            }

            // `cxx:397-404`
            for ii in 1..=ds2d.nb_shapes() {
                let mut li = ds2d.change_shape_interferences_at(ii).clone();
                fun_reducedoublons(&mut li, ds2d, ii);
                *ds2d.change_shape_interferences_at(ii) = li;
            }
        });

        // `cxx:406-407`
        self.my_face_reference = Some(f1.clone());
        let mut wes = WireEdgeSet::new(&f1);

        // `cxx:409-413`
        set_global_faces2d(true);
        self.g_fill_faces_wes_k(lf1, lf2, g1, &mut wes, 1);
        self.g_fill_faces_wes_k(lf1, lf2, g1, &mut wes, 3);
        set_global_faces2d(false);

        // `cxx:415-422`
        let mut f1f = lf1.first().cloned().unwrap_or_else(|| f1.clone());
        f1f.set_orientation(Orientation::Forward);
        let force_class = true;
        let mut fabu = FaceBuilder::new();
        fabu.init_face_builder(&mut wes, &f1f, force_class);

        // `cxx:424-432`
        let mut lfm: Vec<TopoShape> = Vec::new();
        let mut m_wis_old = std::collections::HashMap::new();
        self.g_fabumake_faces(&f1f, &mut fabu, &mut lfm, &mut m_wis_old);

        // `cxx:434-438`
        let mut new_lfm: Vec<TopoShape> = Vec::new();
        self.regularize_faces(&f1f, &lfm, &mut new_lfm);
        lfm = new_lfm;

        // `cxx:440-461`: the new faces are the merged faces of both lists.
        for f1x in lf1 {
            let tomerge = !self.is_merged(f1x, tb1);
            if tomerge {
                *self.change_merged(f1x, tb1) = lfm.clone();
            }
        }
        for f2 in lf2 {
            let tomerge = !self.is_merged(f2, tb2);
            if tomerge {
                *self.change_merged(f2, tb2) = lfm.clone();
            }
        }
    }

    /// `FUN_computeLIFfaces2d(BU, F, pDS2d)` (`cxx:312-327`).
    fn fun_compute_lif_faces2d(&self, f: &TopoShape, p_ds2d: &mut DataStructure) {
        for e in explorer_shapes(f, ShapeType::Edge) {
            if !self.fun_compute_lif_faces2d_edge(f, &e, p_ds2d) {
                return;
            }
        }
    }

    /// `FUN_computeLIFfaces2d(BU, F, E, pDS2d)` (`cxx:144-308`).
    ///
    /// Computes the new face/face interferences `F FTRA` for a same-domain
    /// partner of `F`: `{I = (T(F), ES, FTRA)}` where `ES` interferes with the
    /// split edge `E` of `F`.
    fn fun_compute_lif_faces2d_edge(
        &self,
        f: &TopoShape,
        e: &TopoShape,
        p_ds2d: &mut DataStructure,
    ) -> bool {
        let Some(hds) = self.my_data_structure.as_ref() else {
            return true;
        };
        let bds = hds.ds();

        // `cxx:152-157`
        let li = bds.shape_interferences(e, true).to_vec();
        let ie = bds.shape_index(e, true);
        let if_ = bds.shape_index(f, true);
        let rk_f = bds.ancestor_rank_of(f);

        // `cxx:157-161`
        let mut hassp_e = self.is_split(e, State::On);
        if hassp_e {
            hassp_e = !self.splits(e, State::On).is_empty();
        }

        // `cxx:162-176`
        let mut ff_found = IndexedShapeSet::new();
        let mut fsdm: Vec<TopoShape> = Vec::new();
        for sd in bds.shape_same_domain(f) {
            if bds.ancestor_rank_of(sd) == rk_f {
                continue;
            }
            fsdm.push(sd.clone());
        }

        // `cxx:178-306`
        for i in li {
            let Ok(guard) = i.read() else {
                continue;
            };
            // `cxx:183-189`
            let (sb, ib, _sa, _ia, _gt, _g, st, s) = fds_idata(&guard);
            if st != Kind::Edge {
                return false;
            }

            // `cxx:191-207`
            let mut itra = ib;
            let ftra = if sb == ShapeType::Face {
                bds.shape(ib, true)
            } else if sb == ShapeType::Edge {
                let Some(edge_s) = bds.shape(s, true) else {
                    continue;
                };
                match fun_tool_find_ancestor(&fsdm, &edge_s) {
                    Some(ftra) => {
                        itra = bds.shape_index(&ftra, true);
                        Some(ftra)
                    }
                    None => return false,
                }
            } else {
                None
            };
            let Some(ftra) = ftra else {
                continue;
            };
            let found = ff_found.contains(&ftra);

            // `cxx:217-221`
            let Some(par_e) = guard.parameter() else {
                continue;
            };
            let Some(es) = bds.shape(s, true) else {
                continue;
            };
            let mut hassp_es = self.is_split(&es, State::On);
            if hassp_es {
                hassp_es = !self.splits(&es, State::On).is_empty();
            }
            // `cxx:219`: OCCT assigns to `hasspE` here (a typo that widens the
            // face-level flag); the port keeps the assignment as written.
            if hassp_es {
                hassp_e = !self.splits(&es, State::On).is_empty();
            }

            // `cxx:225-229`
            let sdm = fun_ds_sdm(bds, e, &es);
            let mut mk_ton_esdm = sdm && hassp_e && !found;
            let has_fei_f_e_ftra = fun_ds_has_fei(p_ds2d, f, ie, itra);
            mk_ton_esdm = mk_ton_esdm && !has_fei_f_e_ftra;

            // `cxx:230-244`
            if mk_ton_esdm {
                ff_found.add(&ftra);
                let dummy = Edge::new().0;
                if let Some(mut new_t) =
                    crate::topopebrep::ds::complete_ds::fun_ds_mk_ton_fsdm(
                        hds, if_, itra, s, ie, par_e, &dummy, true,
                    )
                {
                    new_t.set_index(itra);
                    let new_i = crate::topopebrep::ds::interference_tool::make_face_edge_interference(
                        &new_t,
                        itra,
                        ie,
                        true,
                        Config::SameOriented,
                    );
                    p_ds2d.add_shape_interference(f, new_i);
                }
            }

            // `cxx:245-277`
            let mut mk_ton_es_sdm = sdm && hassp_es;
            let has_fei_ftra_e_f = fun_ds_has_fei(p_ds2d, &ftra, ie, if_);
            mk_ton_es_sdm = mk_ton_es_sdm && !has_fei_ftra_e_f;
            if mk_ton_es_sdm {
                let dummy = Edge::new().0;
                let Some(par_es) = fun_tool_par_e(e, par_e, &es) else {
                    continue;
                };
                if let Some(mut new_t) =
                    crate::topopebrep::ds::complete_ds::fun_ds_mk_ton_fsdm(
                        hds, itra, if_, ie, s, par_es, &dummy, true,
                    )
                {
                    new_t.set_index(if_);
                    let new_i = crate::topopebrep::ds::interference_tool::make_face_edge_interference(
                        &new_t,
                        if_,
                        ie,
                        false,
                        Config::SameOriented,
                    );
                    p_ds2d.add_shape_interference(&ftra, new_i);
                }

                if let Some(mut new_t) =
                    crate::topopebrep::ds::complete_ds::fun_ds_mk_ton_fsdm(
                        hds, itra, if_, ie, ie, par_e, &dummy, true,
                    )
                {
                    new_t.set_index(if_);
                    let new_i = crate::topopebrep::ds::interference_tool::make_face_edge_interference(
                        &new_t,
                        if_,
                        s,
                        true,
                        Config::SameOriented,
                    );
                    p_ds2d.add_shape_interference(&ftra, new_i);
                }
            }

            // `cxx:279-305`
            let mut mk_ton_es = hassp_es;
            let has_fei_f_s_ftra = fun_ds_has_fei(p_ds2d, f, s, itra);
            mk_ton_es = mk_ton_es && !has_fei_f_s_ftra;
            if mk_ton_es {
                let Some(par_es) = fun_tool_par_e(e, par_e, &es) else {
                    continue;
                };
                let dummy = Edge::new().0;
                if let Some(mut new_t) =
                    crate::topopebrep::ds::complete_ds::fun_ds_mk_ton_fsdm(
                        hds, if_, itra, s, s, par_es, &dummy, true,
                    )
                {
                    new_t.set_index(itra);
                    let new_i = crate::topopebrep::ds::interference_tool::make_face_edge_interference(
                        &new_t,
                        itra,
                        s,
                        false,
                        Config::SameOriented,
                    );
                    p_ds2d.add_shape_interference(f, new_i);
                }
            }
        }

        true
    }

    /// `GFillFacesWESMakeFaces(LLF1, LF2, LSO, GM)` (`cxx:511-635`).
    pub fn g_fill_faces_wes_make_faces(
        &mut self,
        llf1: &[TopoShape],
        lf2: &[TopoShape],
        gm: &GTopo,
    ) {
        let (tb1, tb2) = gm.states_on();
        if llf1.is_empty() {
            return;
        }

        // `cxx:524-541`: the reference face of the same-domain family first.
        let mut lf1: Vec<TopoShape> = Vec::new();
        for fcur in llf1 {
            let (icur, iref) = match self.my_data_structure.as_ref() {
                Some(hds) => (
                    hds.shape_index(fcur, true),
                    hds.same_domain_reference(fcur),
                ),
                None => (0, 0),
            };
            if icur == iref {
                lf1.insert(0, fcur.clone());
            } else {
                lf1.push(fcur.clone());
            }
        }

        // `cxx:545-548`
        let mut ff = lf1.first().cloned().unwrap_or_else(|| TopoShape::new(ShapeType::Face));
        ff.set_orientation(Orientation::Forward);
        let iff = match self.my_data_structure.as_ref() {
            Some(hds) => hds.shape_index(&ff, true),
            None => 0,
        };

        // `cxx:549`
        let mut wes = WireEdgeSet::new(&ff);

        // `cxx:561-573`
        let mut n1;
        set_global_faces2d(true);
        self.g_fill_faces_wes_k(&lf1, lf2, gm, &mut wes, 1);
        set_global_faces2d(false);
        n1 = wes.start_elements().len() as i32;

        self.g_fill_faces_wes_k(&lf1, lf2, gm, &mut wes, 2);
        n1 = wes.start_elements().len() as i32;

        self.g_fill_faces_wes_k(&lf1, lf2, gm, &mut wes, 3);
        n1 = wes.start_elements().len() as i32;

        // `cxx:577-578`: the start edges added by `GFillCurveTopologyWES` are
        // the ones to avoid from now on.
        let n2 = wes.start_elements().len() as i32;
        self.my_edge_avoid.clear();
        let start_elements = wes.start_elements().to_vec();
        Self::copy_list_range(&start_elements, n1 + 1, n2, &mut self.my_edge_avoid);

        // `cxx:579-580`
        let mut lof: Vec<TopoShape> = Vec::new();
        self.gwes_make_faces(&ff, &mut wes, &mut lof);

        // `cxx:585-597`: drop the faces whose only edge is not a valid bound.
        let mut keep: Vec<TopoShape> = Vec::new();
        for face in lof {
            if fun_valid_f1edge(&face) {
                keep.push(face);
            }
        }
        let lof = keep;

        // `cxx:600-608`
        let mut lofs: Vec<TopoShape> = Vec::new();
        self.g_keep_shapes(&ff, &[], tb1, &lof, &mut lofs);

        // `cxx:606-608`: the faces take the original orientation of `FF`.
        let ods_ff = match self.my_data_structure.as_ref() {
            Some(hds) => hds
                .shape(iff, true)
                .map(|s| s.orientation())
                .unwrap_or(Orientation::Forward),
            None => Orientation::Forward,
        };
        for f in lofs.iter_mut() {
            f.set_orientation(ods_ff);
        }

        // `cxx:610-636`: OCCT appends to the split list without clearing it.
        for s in &lf1 {
            self.mark_split(s, tb1, true);
            let ls1 = self.change_split(s, tb1);
            Self::copy_list(&lofs, ls1);
        }
        for s in lf2 {
            self.mark_split(s, tb2, true);
            let ls2 = self.change_split(s, tb2);
            Self::copy_list(&lofs, ls2);
        }
    }

    /// `GFillFaceWES(FOR1, LFclass, G1, WES)` (`cxx:639-706`).
    pub fn g_fill_face_wes(
        &mut self,
        for1: &TopoShape,
        lf_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        // `cxx:644-646`
        let (tb1, _tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();

        // `cxx:657-666`
        let ope_cut = self.opec12() || self.opec21();
        let com_of_cut = ope_cut && tb1 == _tb2 && tb1 == State::In;
        let hsdm = self
            .my_data_structure
            .as_ref()
            .map(|h| h.has_same_domain(for1, true))
            .unwrap_or(false);
        if hsdm && com_of_cut {
            return;
        }

        // `cxx:668-669`
        let mut ff = for1.clone();
        ff.set_orientation(Orientation::Forward);
        self.my_face_to_fill = Some(ff.clone());

        // `cxx:671-698`
        for w in explorer_shapes(&ff, ShapeType::Wire) {
            let hasshape = self
                .my_data_structure
                .as_ref()
                .map(|h| h.has_shape(&w, true))
                .unwrap_or(false);

            if !hasshape {
                // `cxx:678-692`: the wire is not in the DS, so classify it.
                let mut pos = State::Unknown;
                let keep = self.g_keep_shape1(&w, lf_class, tb1, &mut pos);
                if keep {
                    let ori_w = w.orientation();
                    let new_ori_w = Builder::orient(ori_w, rev_ori1);
                    let mut w2 = w.clone();
                    w2.set_orientation(new_ori_w);
                    wes.add_shape(&w2);
                } else if self.my_process_on && pos == State::On {
                    self.my_on_elem_map.add(&w);
                }
            } else {
                // `cxx:696`: the wire has geometric edges, split them.
                self.g_fill_wire_wes(&w, lf_class, g1, wes);
            }
        }
    }

    /// `GFillWireWES(W, LSclass, G1, WES)` (`cxx:709-757`).
    pub fn g_fill_wire_wes(
        &mut self,
        w: &TopoShape,
        ls_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        for eor in explorer_shapes(w, ShapeType::Edge) {
            self.g_fill_edge_wes(&eor, ls_class, g1, wes);
        }
    }
}

/// `FUN_validF1edge(F)` (`cxx:474-507`): `true` when `F` has more than one
/// distinct edge, or exactly one that closes on itself.
fn fun_valid_f1edge(f: &TopoShape) -> bool {
    let mut m_et = IndexedShapeSet::new();
    let mut n_e = 0;
    for e in explorer_shapes(f, ShapeType::Edge) {
        if m_et.contains(&e) {
            continue;
        }
        m_et.add(&e);
        n_e += 1;
        if n_e > 2 {
            break;
        }
    }
    if n_e > 1 {
        return true;
    }
    if n_e == 1 {
        let edges = explorer_shapes(f, ShapeType::Edge);
        let Some(e) = edges.first() else {
            return false;
        };
        return crate::topopebrep::tool::closed_e(e).is_some();
    }
    false
}
