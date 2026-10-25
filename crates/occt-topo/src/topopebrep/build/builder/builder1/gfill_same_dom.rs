//! `TopOpeBRepBuild_Builder1` — the same-domain face family.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GFillFaceSameDomSFS(FOR, LSO2, Gin, SFS)` | `533-613` |
//! | `GFillFaceSameDomWES(FOR1, LFclass, G1, WES)` | `614-746` |
//! | `GFillWireSameDomWES(W, LSclass, G1, WES)` | `747-817` |
//! | `GFillEdgeSameDomWES(EOR, LSclass, G1, WES)` | `818-1150` |
//! | `GWESMakeFaces(FF, WES, LOF)` (virtual, `hxx:80-82`) | `1373-1412` |
//!
//! `GFillFaceSameDomWES` walks the same-domain closure of `FOR1` while it
//! grows (`cxx:628-641`): every face reached adds its own same-domain faces to
//! the working map and its edges to `mySDEdgeMap`, which `PerformONParts`
//! reads. Faces of the tool operand (`AncestorRank == 2`) fill with the
//! permuted grid, which is why the local `GFTW` is copied per face.
//!
//! `GWESMakeFaces` is the third virtual override: it runs the base member and
//! then normalizes (`NormalizeFace`) and 2d-corrects (`CorrectFace2d`) every
//! built face. It is reached through `Builder::gwes_make_faces`.
//!
//! UNPORTED: none. `TopOpeBRepBuild_Tools::NormalizeFace` and
//! `CorrectFace2d` are both ported (`tools::normalize_face`,
//! `tools::correct_face_2d`).

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key};
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::shell_face_set::ShellFaceSet;
use crate::topopebrep::build::tools::{
    correct_face_2d, get_normal_in_nearest_point, get_normal_to_face_on_edge,
    get_tangent_to_edge, normalize_face, update_edge_on_face, update_pcurves,
};
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;

use super::fill_2d::classify_edge_to_face_by_one_point;
use super::global_use_new_builder;
use super::super::Builder;

/// `TopOpeBRepBuild_Tools::GetTangentToEdge` / `GetNormalToFaceOnEdge`
/// answer `false` and leave the vector untouched; the port passes a zero
/// vector in that case (`fill_2d.rs` does the same).
fn zero_vec() -> occt_core::gp::GpVec {
    occt_core::gp::GpVec::new(0.0, 0.0, 0.0)
}

impl Builder {
    /// `GFillFaceSameDomSFS(FOR, LSO2, Gin, SFS)` (`cxx:533-613`).
    pub(crate) fn b1_g_fill_face_same_dom_sfs(
        &mut self,
        for_: &TopoShape,
        lso2: &[TopoShape],
        gin: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        // `cxx:538-542`
        self.b1.my_processed_parts_out_2d.clear();
        self.b1.my_processed_parts_on_2d.clear();
        self.b1.my_map_of_edge_with_face_state.clear();
        super::sd_edge_map_clear();
        self.b1.my_source_shapes.clear();

        // `cxx:544-547`: only the object operand drives the same-domain pass.
        let rank = self
            .data_structure()
            .map(|hds| hds.ds().ancestor_rank_of(for_))
            .unwrap_or(0);
        if rank != 1 {
            return;
        }
        let g1 = gin;

        let (tb1, _tb2) = g1.states_on();

        // `cxx:556-557`
        let mut ff = for_.clone();
        ff.set_orientation(Orientation::Forward);

        // `cxx:560`
        let mut wes = WireEdgeSet::new(&ff);

        // `cxx:563`
        self.b1_g_fill_face_same_dom_wes(for_, lso2, g1, &mut wes);

        // `cxx:565-568`
        self.my_edge_avoid.clear();
        self.mark_split(&ff, tb1, true);

        // `cxx:571-572`
        let mut lof: Vec<TopoShape> = Vec::new();
        self.gwes_make_faces(&ff, &mut wes, &mut lof);

        // `cxx:575-607`: orient every new face and hand it to the set.
        let orig_rev = for_.orientation() != Orientation::Forward;
        let mut ori_lof: Vec<TopoShape> = Vec::new();
        for a_face in lof.into_iter() {
            let mut a_face = a_face;
            let mut rev = false;
            for an_edge in explorer_shapes(&a_face, ShapeType::Edge) {
                // `myMapOfEdgeWithFaceState.Find(anEdge, rev)` stops at the
                // first edge the map knows.
                if let Some(v) = self
                    .b1
                    .my_map_of_edge_with_face_state
                    .get(&shape_map_key(&an_edge))
                {
                    rev = *v;
                    break;
                }
            }
            if orig_rev {
                a_face.reverse();
            }
            if rev {
                a_face.reverse();
            }
            ori_lof.push(a_face.clone());
            sfs.add_start_element(&a_face);
        }

        // `cxx:609-610`
        let mut lofs: Vec<TopoShape> = Vec::new();
        self.base_g_keep_shapes_same_dom(&ff, tb1, &ori_lof, &mut lofs);
        let lo = self.change_split(&ff, tb1);
        lo.clear();
        lo.extend(lofs);
    }

    /// `myEmptyShapeList` is the `GKeepShapes` class list of
    /// `GFillFaceSameDomSFS` (`cxx:609`).
    fn base_g_keep_shapes_same_dom(
        &mut self,
        ff: &TopoShape,
        tb1: State,
        ori_lof: &[TopoShape],
        lofs: &mut Vec<TopoShape>,
    ) {
        self.g_keep_shapes(ff, &[], tb1, ori_lof, lofs);
    }

    /// `GFillFaceSameDomWES(FOR1, LFclass, G1, WES)` (`cxx:614-746`).
    pub(crate) fn b1_g_fill_face_same_dom_wes(
        &mut self,
        for1: &TopoShape,
        _lf_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        let (tb1, tb2) = g1.states_on();

        // `cxx:622`
        self.b1.my_base_face_to_fill = Some(for1.clone());

        // `cxx:625-628`
        let mut cur_same_dom: Vec<TopoShape> = vec![for1.clone()];

        // `cxx:629-641`: the closure grows while it is walked.
        let mut i = 0usize;
        while i < cur_same_dom.len() {
            let cur = cur_same_dom[i].clone();
            let same_dom = self
                .data_structure()
                .map(|hds| hds.same_domain(&cur).to_vec())
                .unwrap_or_default();
            for sdf in same_dom {
                if cur_same_dom.iter().any(|s| is_same(s, &sdf)) {
                    continue;
                }
                cur_same_dom.push(sdf.clone());
                self.b1.my_same_dom_map.add(&sdf);
                super::sd_edge_map_add_shapes(&sdf);
            }
            i += 1;
        }

        // `cxx:643-645`
        let n_f = cur_same_dom.len();
        for cur_f in cur_same_dom.iter().take(n_f) {
            let mut cur_ff = cur_f.clone();
            cur_ff.set_orientation(Orientation::Forward);

            // `cxx:654`
            self.b1.my_sd_face_to_fill = Some(cur_f.clone());
            let iref = self
                .data_structure()
                .map(|hds| hds.ds().ancestor_rank_of(&cur_ff))
                .unwrap_or(0);

            // `cxx:657-669`: the grid follows the ancestor rank.
            let (tb, rev_ori) = if iref == 1 {
                (tb1, g1.is_to_reverse1())
            } else {
                (tb2, g1.is_to_reverse2())
            };
            let mut gftw = g1.clone();
            if iref != 1 && rev_ori {
                gftw = g1.copy_permuted();
            }

            // `cxx:672`
            self
                .g_fill_curve_topology_wes(&cur_ff, &gftw, wes);

            // `cxx:675`: the `SDFaces` argument is unused by OCCT.
            self.perform_on_parts(&cur_ff, g1, wes);

            // `cxx:677`
            let lsf: Vec<TopoShape> = self
                .data_structure()
                .map(|hds| hds.ds().shape_same_domain(&cur_ff).to_vec())
                .unwrap_or_default();

            // `cxx:679-742`
            for w in explorer_shapes(&cur_ff, ShapeType::Wire) {
                let hasshape = self
                    .data_structure()
                    .map(|hds| hds.has_shape(&w, true))
                    .unwrap_or(false);
                let sh_st = self
                    .data_structure()
                    .and_then(|hds| hds.ds().get_shape_with_state(&w))
                    .map(|sws| sws.state())
                    .unwrap_or(State::Unknown);

                if !hasshape && sh_st != State::On {
                    // `cxx:688-735`
                    let keep = sh_st == tb;
                    if !keep {
                        continue;
                    }
                    let mut w = w;
                    let newori_w = Builder::orient(w.orientation(), rev_ori);

                    let base_face = self.b1.my_base_face_to_fill.clone();
                    let sd_face = self.b1.my_sd_face_to_fill.clone();
                    let same = match (base_face.as_ref(), sd_face.as_ref()) {
                        (Some(a), Some(b)) => is_same(a, b),
                        _ => true,
                    };
                    if !same {
                        if let (Some(b), Some(s)) = (base_face.clone(), sd_face.clone()) {
                            if let (Some(bf), Some(sf)) = (Face::wrap(b), Face::wrap(s)) {
                                update_pcurves(&w, &sf, &bf);
                            }
                        }
                    } else {
                        self.b1.my_source_shapes.add(&w);
                    }

                    let mut use_edges = false;
                    for eor in explorer_shapes(&w, ShapeType::Edge) {
                        let oldori = eor.orientation();
                        let mut eor = eor;
                        let state_of_face_ori = self.b1_orientate_edge_on_face(
                            &mut eor,
                            g1,
                        );
                        if eor.orientation() != oldori {
                            use_edges = true;
                            wes.add_start_element(&eor);
                        }
                        self
                            .b1
                            .my_map_of_edge_with_face_state
                            .insert(shape_map_key(&eor), state_of_face_ori);
                    }

                    if !use_edges {
                        w.set_orientation(newori_w);
                        wes.add_shape(&w);
                    }
                } else {
                    // `cxx:737-740`
                    self.b1_g_fill_wire_same_dom_wes(&w, &lsf, g1, wes);
                }
            }
        }
    }

    /// `OrientateEdgeOnFace(E, myBaseFaceToFill, mySDFaceToFill, G1, s)`
    /// (`cxx:2227-2290`) with the two faces the fill state holds.
    fn b1_orientate_edge_on_face(&self, edge: &mut TopoShape, g1: &GTopo) -> bool {
        let base = self.b1.my_base_face_to_fill.clone();
        let sd = self.b1.my_sd_face_to_fill.clone();
        let (Some(base), Some(sd)) = (base, sd) else {
            return false;
        };
        let (Some(base_f), Some(sd_f)) = (Face::wrap(base), Face::wrap(sd)) else {
            return false;
        };
        self.orientate_edge_on_face(edge, &base_f, &sd_f, g1)
    }

    /// `GFillWireSameDomWES(W, LSclass, G1, WES)` (`cxx:747-817`).
    pub(crate) fn b1_g_fill_wire_same_dom_wes(
        &mut self,
        w: &TopoShape,
        ls_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        let (tb1, tb2) = g1.states_on();

        // `cxx:757-770`: the rank picks the state and the reverse flag.
        let iref = self
            .data_structure()
            .map(|hds| hds.ds().ancestor_rank_of(w))
            .unwrap_or(0);
        let (tb, rev_ori) = if iref == 1 {
            (tb1, g1.is_to_reverse1())
        } else {
            (tb2, g1.is_to_reverse2())
        };

        // `cxx:772-814`
        for eor in explorer_shapes(w, ShapeType::Edge) {
            let hasshape = self
                .data_structure()
                .map(|hds| hds.has_shape(&eor, true))
                .unwrap_or(false);

            if !hasshape {
                let sh_st = self
                    .data_structure()
                    .and_then(|hds| hds.ds().get_shape_with_state(&eor))
                    .map(|sws| sws.state())
                    .unwrap_or(State::Unknown);
                let keep = sh_st == tb;
                if !keep {
                    continue;
                }
                // `cxx:786`: `Orient(oriE, RevOri)` is called and its result
                // dropped; kept for the side-effect-free fidelity of the
                // translation.
                let _ = Builder::orient(eor.orientation(), rev_ori);

                let base_face = self.b1.my_base_face_to_fill.clone();
                let sd_face = self.b1.my_sd_face_to_fill.clone();
                let same = match (base_face.as_ref(), sd_face.as_ref()) {
                    (Some(a), Some(b)) => is_same(a, b),
                    _ => true,
                };
                if !same {
                    if let (Some(b), Some(s)) = (base_face.clone(), sd_face.clone()) {
                        if let (Some(bf), Some(sf)) = (Face::wrap(b), Face::wrap(s)) {
                            if let Some(e) = Edge::wrap(eor.clone()) {
                                update_edge_on_face(&e, &sf, &bf);
                            }
                        }
                    }
                } else {
                    self.b1.my_source_shapes.add(&eor);
                }

                let mut eor = eor;
                let state_of_face_ori = self.b1_orientate_edge_on_face(&mut eor, g1);
                self
                    .b1
                    .my_map_of_edge_with_face_state
                    .insert(shape_map_key(&eor), state_of_face_ori);
                wes.add_element(&eor);
            } else {
                // `cxx:811-813`
                self.b1_g_fill_edge_same_dom_wes(&eor, ls_class, g1, wes);
            }
        }
    }

    /// `GFillEdgeSameDomWES(EOR, LSclass, G1, WES)` (`cxx:818-1150`).
    pub(crate) fn b1_g_fill_edge_same_dom_wes(
        &mut self,
        eor: &TopoShape,
        ls_class: &[TopoShape],
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        let (tb1, tb2) = g1.states_on();

        // `cxx:826-838`
        let iref = self
            .data_structure()
            .map(|hds| hds.ds().ancestor_rank_of(eor))
            .unwrap_or(0);
        let (tb, rev_ori) = if iref == 1 {
            (tb1, g1.is_to_reverse1())
        } else {
            (tb2, g1.is_to_reverse2())
        };

        // `cxx:840-841`: the result of `Orient` is dropped, as in OCCT.
        let ori_e = eor.orientation();
        let _ = Builder::orient(ori_e, rev_ori);

        // 1) the split parts of state `TB` (`cxx:844-872`).
        for new_e in self.parts_in_state(eor, tb) {
            let mut new_e = new_e;
            new_e.set_orientation(ori_e);
            self.b1_update_on_fill_faces(&mut new_e);
            let state_of_face_ori = self.b1_orientate_edge_on_face(&mut new_e, g1);
            self
                .b1
                .my_map_of_edge_with_face_state
                .insert(shape_map_key(&new_e), state_of_face_ori);
            wes.add_start_element(&new_e);
        }

        // 2) the ON parts (`cxx:875-1145`).
        for a_split_part in self.parts_in_state(eor, State::On) {
            let mut a_split_part = a_split_part;
            a_split_part.set_orientation(eor.orientation());
            self.b1_fill_edge_on_part(eor, &a_split_part, ls_class, g1, tb, ori_e, wes);
        }
    }

    /// `UpdateEdgeOnFace` or `mySourceShapes.Add`, as the base/target face
    /// pair decides (`cxx:851-866`, `:1090-1098`, `:1138-1141`).
    fn b1_update_on_fill_faces(&mut self, edge: &mut TopoShape) {
        let base_face = self.b1.my_base_face_to_fill.clone();
        let sd_face = self.b1.my_sd_face_to_fill.clone();
        let same = match (base_face.as_ref(), sd_face.as_ref()) {
            (Some(a), Some(b)) => is_same(a, b),
            _ => true,
        };
        if !same {
            if let (Some(b), Some(s)) = (base_face, sd_face) {
                if let (Some(bf), Some(sf), Some(e)) =
                    (Face::wrap(b), Face::wrap(s), Edge::wrap(edge.clone()))
                {
                    update_edge_on_face(&e, &sf, &bf);
                }
            }
        } else {
            self.b1.my_source_shapes.add(edge);
        }
    }
}

impl Builder {
    /// `GFillEdgeSameDomWES` — one ON part of the edge (`cxx:878-1145`).
    ///
    /// Split out of [`Builder1::b1_g_fill_edge_same_dom_wes`] only to keep the
    /// member readable; this is exactly the body of the `it` loop at
    /// `cxx:878-1145`.
    #[allow(clippy::too_many_arguments)]
    fn b1_fill_edge_on_part(
        &mut self,
        eor: &TopoShape,
        a_split_part: &TopoShape,
        ls_class: &[TopoShape],
        g1: &GTopo,
        tb: State,
        ori_e: Orientation,
        wes: &mut WireEdgeSet,
    ) {
        // ---- ON2D (`cxx:884-947`) --------------------------------------
        let has_sd = self.has_same_domain(a_split_part) || self.has_same_domain(eor);
        if has_sd {
            // `cxx:888-895`: a REVERSED `mySDFaceToFill` flips both shapes.
            let mut e_on = a_split_part.clone();
            let mut n_eor = eor.clone();
            let sd_face_ori = self
                .b1
                .my_sd_face_to_fill
                .as_ref()
                .map(|f| f.orientation());
            if sd_face_ori == Some(Orientation::Reversed) {
                e_on.reverse();
                n_eor.reverse();
            }

            // `cxx:898-901`: the OUT2D pieces are computed once, for the object.
            if self.b1.my_processed_parts_out_2d.contains(&e_on) {
                return;
            }

            let Some(sd_face_shape) = self.b1.my_sd_face_to_fill.clone() else {
                return;
            };
            let mut a_pieces: Vec<TopoShape> = Vec::new();
            let mut a_faces: Vec<TopoShape> = Vec::new();
            let mut a_pieces_out_2d: Vec<TopoShape> = Vec::new();
            let flag = self.perform_piece_on_2d(
                &e_on,
                &sd_face_shape,
                &n_eor,
                &mut a_pieces,
                &mut a_faces,
                &mut a_pieces_out_2d,
            );

            // `cxx:904-914`
            for a_p in a_pieces_out_2d {
                let mut a_fp = a_p.clone();
                a_fp.reverse();
                wes.add_start_element(&a_fp);
                wes.add_start_element(&a_p);
                self.b1.my_processed_parts_out_2d.add(&a_fp);
            }

            // `cxx:915-941`: only the pieces lying on the face being filled
            // are kept, oriented like the original edge.
            for (a_piece, a_piece_face) in a_pieces.iter().zip(a_faces.iter()) {
                if !is_same(a_piece_face, &sd_face_shape) {
                    continue;
                }
                let mut a_piece_to_keep = a_piece.clone();
                let is_rev = a_piece_to_keep.orientation() == n_eor.orientation();
                a_piece_to_keep.set_orientation(ori_e);
                let state_of_face_ori =
                    self.b1_orientate_edge_on_face(&mut a_piece_to_keep, g1);
                // `cxx:930-934`: undo the original orientation mismatch.
                if !is_rev {
                    a_piece_to_keep.reverse();
                }
                self
                    .b1
                    .my_map_of_edge_with_face_state
                    .insert(shape_map_key(&a_piece_to_keep), state_of_face_ori);
                wes.add_start_element(&a_piece_to_keep);
            }

            // `cxx:941-944`: `flag != 0` means the 2d pass decided.
            if flag != 0 {
                return;
            }
        }

        // ---- IN2D classification (`cxx:949-978`) ------------------------
        let mut a_sd_shape: Option<TopoShape> = None;
        let mut a_state = State::Unknown;
        if ls_class.len() == 1 {
            a_sd_shape = Some(ls_class[0].clone());
            if let (Some(e), Some(f)) = (Edge::wrap(a_split_part.clone()), Face::wrap(ls_class[0].clone()))
            {
                a_state = classify_edge_to_face_by_one_point(&e, &f);
            }
        } else {
            // A face with several same-domain twins: the first one that puts
            // the part IN or ON wins (`cxx:958-972`).
            for cur_sd in ls_class.iter() {
                let Some(f) = Face::wrap(cur_sd.clone()) else {
                    continue;
                };
                let Some(e) = Edge::wrap(a_split_part.clone()) else {
                    continue;
                };
                a_state = classify_edge_to_face_by_one_point(&e, &f);
                if a_state == State::In || a_state == State::On {
                    a_sd_shape = Some(cur_sd.clone());
                    break;
                }
            }
        }

        // `cxx:976-979`: the ON and UNKNOWN parts are handled above.
        if a_state == State::On || a_state == State::Unknown {
            return;
        }

        // ---- OUT2D (`cxx:982-1112`) ------------------------------------
        if a_state == State::Out || a_sd_shape.is_none() {
            let keep = self.b1_edge_out_2d_keep(eor, a_split_part, tb);

            // `cxx:1089-1110`
            if keep {
                let mut a_split_part = a_split_part.clone();
                self.b1_update_on_fill_faces(&mut a_split_part);
                let state_of_face_ori =
                    self.b1_orientate_edge_on_face(&mut a_split_part, g1);
                self
                    .b1
                    .my_map_of_edge_with_face_state
                    .insert(shape_map_key(&a_split_part), state_of_face_ori);
                wes.add_start_element(&a_split_part);
            }
            return;
        }

        // ---- IN2D (`cxx:1116-1144`) ------------------------------------
        let Some(a_sd) = a_sd_shape else {
            return;
        };
        let Some(sd_face_shape) = self.b1.my_sd_face_to_fill.clone() else {
            return;
        };
        let (Some(a_piece), Some(a_orig), Some(edge_face), Some(to_face)) = (
            Edge::wrap(a_split_part.clone()),
            Edge::wrap(eor.clone()),
            Face::wrap(sd_face_shape.clone()),
            Face::wrap(a_sd.clone()),
        ) else {
            return;
        };
        let keep = self.perform_piece_in_2d(&a_piece, &a_orig, &edge_face, &to_face, g1);
        if !keep {
            return;
        }

        // `cxx:1131-1141`
        self.b1.my_splits_on_to_keep.add(a_split_part);

        let mut a_split = a_split_part.clone();
        let base_face = self.b1.my_base_face_to_fill.clone();
        let state_of_face_ori = match base_face {
            Some(b) => match (Face::wrap(b), Face::wrap(a_sd.clone())) {
                (Some(b), Some(t)) => self.orientate_edge_on_face(&mut a_split, &b, &t, g1),
                _ => false,
            },
            None => false,
        };
        self
            .b1
            .my_map_of_edge_with_face_state
            .insert(shape_map_key(&a_split), state_of_face_ori);

        if is_same(
            self
                .b1
                .my_base_face_to_fill
                .as_ref()
                .unwrap_or(&sd_face_shape),
            &sd_face_shape,
        ) {
            self.b1.my_source_shapes.add(&a_split);
        }
        wes.add_start_element(&a_split);
    }

    /// The OUT2D keep decision of `GFillEdgeSameDomWES` (`cxx:985-1087`).
    fn b1_edge_out_2d_keep(
        &self,
        eor: &TopoShape,
        a_split_part: &TopoShape,
        tb: State,
    ) -> bool {
        // `cxx:996-1017`: the faces of the other operand under `EOR`, at most
        // two of them, the first same-domain one is the adjacent face.
        let a_fel = self
            .b1
            .my_map_of_edge_faces
            .find_from_key(eor)
            .cloned()
            .unwrap_or_default();
        let mut a_adj_sd_face: Option<TopoShape> = None;
        if a_fel.len() <= 2 {
            for a_adj in a_fel.iter() {
                if let Some(sdf) = self.b1.my_sd_face_to_fill.as_ref() {
                    if is_same(sdf, a_adj) {
                        continue;
                    }
                }
                if self.has_same_domain(a_adj) {
                    a_adj_sd_face = Some(a_adj.clone());
                    break;
                }
            }
        }

        let Some(a_adj_sd_face) = a_adj_sd_face else {
            // `cxx:1085-1087`: no adjacent face, a plain IN/OUT test.
            return tb == State::In;
        };

        // `cxx:1021-1030`: the copy of `EOR` under the adjacent face and the
        // first same-domain face of that adjacent face.
        let adj_eor = explorer_shapes(&a_adj_sd_face, ShapeType::Edge)
            .into_iter()
            .find(|e| is_same(e, eor));
        let a_sd_to_adj = self
            .data_structure()
            .and_then(|hds| hds.same_domain(&a_adj_sd_face).first().cloned());
        let (Some(adj_eor), Some(a_sd_to_adj)) = (adj_eor, a_sd_to_adj) else {
            // `aEAdjMap.FindIndex(EOR)` returning 0 (`cxx:1023`) or an empty
            // `SameDomain` list (`cxx:1028`) makes OCCT read out of range; the
            // port keeps the "no adjacent face" answer instead.
            return tb == State::In;
        };

        let mut a_split_p = a_split_part.clone();
        a_split_p.set_orientation(adj_eor.orientation());
        let Some(a_split_p_e) = Edge::wrap(a_split_p.clone()) else {
            return tb == State::In;
        };

        // `cxx:1032-1045`
        let mut a_n2 = match Face::wrap(a_sd_to_adj.clone()) {
            Some(f) => get_normal_to_face_on_edge(&f, &a_split_p_e).unwrap_or_else(zero_vec),
            None => zero_vec(),
        };
        if a_sd_to_adj.orientation() == Orientation::Reversed {
            a_n2.reverse();
        }
        let mut a_n3 = match Face::wrap(a_adj_sd_face.clone()) {
            Some(f) => get_normal_to_face_on_edge(&f, &a_split_p_e).unwrap_or_else(zero_vec),
            None => zero_vec(),
        };
        if a_adj_sd_face.orientation() == Orientation::Reversed {
            a_n3.reverse();
        }
        // OCCT computes `aN3` and never reads it (`cxx:1040-1045`).
        let _ = a_n3;

        // `cxx:1047-1051`
        let mut a_tg = get_tangent_to_edge(&a_split_p_e).unwrap_or_else(zero_vec);
        if a_split_p.orientation() == Orientation::Reversed {
            a_tg.reverse();
        }
        let a_bi_n = a_tg.crossed(&a_n2);

        // `cxx:1055-1078`
        let Some(sd_face_shape) = self.b1.my_sd_face_to_fill.clone() else {
            return tb == State::In;
        };
        let Some(sd_face) = Face::wrap(sd_face_shape.clone()) else {
            return tb == State::In;
        };
        let mut a_n1 = get_normal_to_face_on_edge(&sd_face, &a_split_p_e).unwrap_or_else(zero_vec);
        if sd_face_shape.orientation() == Orientation::Reversed {
            a_n1.reverse();
        }
        let mut scalar_pr = a_bi_n.dot(&a_n1);
        if scalar_pr.abs() <= 1.0e-10 {
            // `cxx:1064-1076`: step inside the face to get a usable normal.
            if let Some(n) = get_normal_in_nearest_point(&sd_face, &a_split_p_e) {
                a_n1 = n;
            }
            if sd_face_shape.orientation() == Orientation::Reversed {
                a_n1.reverse();
            }
            scalar_pr = a_bi_n.dot(&a_n1);
            // `cxx:1074-1077` assigns `keep` here and `cxx:1079-1081`
            // overwrites it unconditionally, so only the later test matters.
        }
        let a_part_state = if scalar_pr > 0.0 {
            State::In
        } else {
            State::Out
        };
        a_part_state == tb
    }
}

