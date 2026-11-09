//! The rebuild driver `IntersectFaces(const Map& theVertsToAvoid, const
//! Message_ProgressRange&)` (`cxx:5924-6565`).
//!
//! The driver judges one invalidity at a time and intersects only the faces
//! that neighbour it:
//!
//! 1. `cxx:5937-5957`: the vertices of the invalid edges, split into the
//!    "all" set and the "still invalid" set. `bLookVertToAvoid` follows the
//!    second set's emptiness.
//! 2. `cxx:5959-5986`: `PrepareFacesForIntersection` builds the lookup tables
//!    and `FindVerticesToAvoid` extends the caller's avoid set.
//! 3. `cxx:5997-6162`: the blocks of artificially invalid faces. Per block a
//!    face belongs to only one offset face (`aMFArt.Extent() == 1`, `:6034-6039`)
//!    makes it `bAlone`; the invalid edges of the block give `aMEInvOnArt` and,
//!    for the alone case, the images of their origins give `aMEAvoid`
//!    (`:6056-6079`). The remaining edges are mapped vertex -> edge
//!    (`:6082-6086`), and an edge whose vertex belongs to the invalidity and
//!    that has no second edge there becomes an "alone edge" (`:6090-6110`).
//!    When such an edge adds nothing to the intersection table
//!    (`aDMSF(aE).Extent() <= 1`, `:6113-6116`) and none of its vertices sits on
//!    more than two faces (`:6118-6133`), the other images of its origins that
//!    were already mapped are added too (`:6135-6160`).
//! 4. `cxx:6164-6187`: `aMAllInvs`, the invalidities of all faces, used to
//!    avoid reusing a common edge.
//! 5. `cxx:6205-6551`: per invalid face, per connexity block of its splits
//!    (`:6219-6240`), per connexity block of the block's invalid edges
//!    (`:6283-6288`), the faces to intersect are looked up
//!    (`FindFacesForIntersection`, `:6311-6319`), pairs closer than three
//!    faces are skipped (`:6320-6325`), and every remaining pair either
//!    reuses the common edges (`ProcessCommonEdges`, `:6436-6465`), reuses a
//!    previously computed intersection line (`UpdateIntersectedFaces`,
//!    `:6474-6499`) or intersects the two faces (`IntersectFaces`,
//!    `:6505-6518`). Afterwards the intersection edges of the block are
//!    trimmed (`IntersectAndTrimEdges`, `:6522-6532`) and recorded per
//!    invalidity edge (`:6535-6550`).
//!
//! The reference ends with `UpdateValidEdges(aFLE, aDMOENEdges, aMVBounds,
//! aMEInvOnArt, aMECheckExt, theVertsToAvoid, aEImages, aDMEETrim, ...)`
//! (`cxx:6556-6564`), which the port calls in the same place with the same
//! maps.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::Dmvv;
use crate::shape::TopoShape;
use crate::topo_tools::{
    children_keep_orientation, explorer_shapes, is_same, shape_map_key, ShapeMapKey,
};

use super::invalid_edges_filter::make_compound_from;
use super::rebuild_prep::{map_shapes, map_vertices_to_edges};
use super::statics::{append_to_list, find_common_parts};
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::IntersectFaces(Map& theVertsToAvoid,
    /// const Message_ProgressRange& theRange)` (`cxx:5924-6565`).
    ///
    /// Returns early where the reference does when nothing has to be rebuilt
    /// (`cxx:5928-5932`); otherwise it ends with `UpdateValidEdges`.
    pub fn intersect_faces_driver(&mut self, the_verts_to_avoid: &mut ShapeIndexedMap) {
        // `cxx:5928-5932`.
        if self.my_faces_to_rebuild.extent() == 0 {
            return;
        }

        // `cxx:5936-5957`: the vertices of the invalid edges. An edge that is
        // also valid contributes to `aMVInvAll` only.
        let mut a_mv_inv = ShapeIndexedMap::new();
        let mut a_mv_inv_all = ShapeIndexedMap::new();
        let a_nb_inv = self.my_invalid_edges.extent();
        for i in 1..=a_nb_inv {
            let Some(a_e_inv) = self.my_invalid_edges.value(i).cloned() else {
                continue;
            };
            let b_valid = self.my_valid_edges.contains(&a_e_inv);
            for a_v in explorer_shapes(&a_e_inv, ShapeType::Vertex) {
                a_mv_inv_all.add(&a_v);
                if !b_valid {
                    a_mv_inv.add(&a_v);
                }
            }
        }

        // `cxx:5957`.
        let b_look_vert_to_avoid = a_mv_inv.extent() > 0;

        // `cxx:5959-5975`.
        let mut a_fle = Dmvv::new();
        let mut a_m_done: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut a_dmsf: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut a_me_inf_e_trim: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut a_dmve_full: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut a_dmef_inv = Dmvv::new();
        self.prepare_faces_for_intersection(
            b_look_vert_to_avoid,
            &mut a_fle,
            &mut a_m_done,
            &mut a_dmsf,
            &mut a_me_inf_e_trim,
            &mut a_dmve_full,
            &mut a_dmef_inv,
        );

        // `cxx:5982-5983`.
        let mut a_mvr_inv = the_verts_to_avoid.clone();
        self.find_vertices_to_avoid(&a_dmef_inv, &a_dmve_full, &mut a_mvr_inv);

        // `cxx:5999-6011`: split -> the artificial invalid face it belongs to,
        // and one compound of all the artificial splits.
        let mut a_dmf_imf: HashMap<ShapeMapKey, TopoShape> = HashMap::new();
        let mut a_cf_art_splits: Vec<TopoShape> = Vec::new();
        // `NCollection_DataMap::Iterator` (`cxx:6002-6011`) walks the bound
        // faces; the port keeps them in bind order.
        let a_art_keys: Vec<TopoShape> = self.my_art_invalid_face_order.clone();
        for a_f in &a_art_keys {
            let a_lf_inv: Vec<TopoShape> = self
                .my_invalid_faces
                .seek(a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for a_f_im in a_lf_inv {
                a_cf_art_splits.push(a_f_im.clone());
                a_dmf_imf.insert(shape_map_key(&a_f_im), a_f.clone());
            }
        }
        let a_cf_art = make_compound_from(a_cf_art_splits.iter());
        // `cxx:6014-6015`.
        let a_lcb_art =
            AlgoTools::make_connexity_blocks_in(&a_cf_art, ShapeType::Vertex, ShapeType::Face);

        // `cxx:6018-6019`.
        let mut a_me_alone = ShapeIndexedMap::new();
        let mut a_me_inv_on_art = ShapeIndexedMap::new();

        // `cxx:6021-6162`.
        for a_cb in &a_lcb_art {
            // `cxx:6031-6039`: the block is "alone" when its splits come from a
            // single offset face.
            let mut a_mf_art = ShapeIndexedMap::new();
            for a_exp_f in explorer_shapes(a_cb, ShapeType::Face) {
                if let Some(a_f) = a_dmf_imf.get(&shape_map_key(&a_exp_f)) {
                    a_mf_art.add(a_f);
                }
            }
            let b_alone = a_mf_art.extent() == 1;

            // `cxx:6042-6050`.
            let mut a_mve_inv = ShapeIndexedMap::new();
            let mut a_mfence = ShapeIndexedMap::new();
            let mut a_me_avoid = ShapeIndexedMap::new();
            let mut a_dmve_val = Dmvv::new();

            // `cxx:6053-6087`.
            for a_e in explorer_shapes(a_cb, ShapeType::Edge) {
                if self.my_invalid_edges.contains(&a_e) {
                    a_me_inv_on_art.add(&a_e);
                    for a_v in children_keep_orientation(&a_e) {
                        a_mve_inv.add(&a_v);
                    }
                    // `cxx:6064-6078`: the images of the origins of an alone
                    // edge must not become alone edges themselves.
                    if b_alone {
                        if let Some(p_le_or) = self.my_oe_origins.get(&shape_map_key(&a_e)).cloned()
                        {
                            for a_e_or in &p_le_or {
                                // `myOEImages.Find` raises on an unbound edge
                                // in the reference; the port skips it.
                                let Some(a_le_im) =
                                    self.my_oe_images.seek(a_e_or).map(|l| l.to_vec())
                                else {
                                    continue;
                                };
                                for a_e_im in &a_le_im {
                                    a_me_avoid.add(a_e_im);
                                }
                            }
                        }
                    }
                    continue;
                }
                // `cxx:6082-6086`.
                if a_mfence.add(&a_e) == a_mfence.extent() {
                    map_vertices_to_edges(&a_e, &mut a_dmve_val);
                }
            }

            // `cxx:6090-6161`: the edges with a free vertex.
            let a_nb_v = a_dmve_val.extent();
            for i in 1..=a_nb_v {
                let Some(a_v) = a_dmve_val.keys().get(i - 1).cloned() else {
                    continue;
                };
                if !a_mve_inv.contains(&a_v) {
                    continue;
                }
                let a_lev: Vec<TopoShape> = a_dmve_val.seek(&a_v).map(|l| l.to_vec()).unwrap_or_default();
                // `cxx:6100-6103`.
                if a_lev.len() > 1 {
                    continue;
                }
                let Some(a_e) = a_lev.first().cloned() else {
                    continue;
                };
                // `cxx:6107-6110`.
                if a_me_avoid.contains(&a_e) {
                    continue;
                }
                a_me_alone.add(&a_e);

                // `cxx:6113-6116`: the edge adds nothing to the intersection
                // table.
                let a_dmsf_e = a_dmsf
                    .get(&shape_map_key(&a_e))
                    .map(|l| l.len())
                    .unwrap_or(0);
                if a_dmsf_e > 1 {
                    continue;
                }
                // `cxx:6118-6133`: none of the vertices may sit on more than
                // two faces.
                let mut b_many_faces = false;
                for a_ve in children_keep_orientation(&a_e) {
                    let a_extent = a_dmsf
                        .get(&shape_map_key(&a_ve))
                        .map(|l| l.len())
                        .unwrap_or(0);
                    if a_extent > 2 {
                        b_many_faces = true;
                        break;
                    }
                }
                if b_many_faces {
                    continue;
                }
                // `cxx:6135-6160`: look for the other images of the same
                // origins.
                let Some(p_le_or) = self.my_oe_origins.get(&shape_map_key(&a_e)).cloned() else {
                    continue;
                };
                for a_e_or in &p_le_or {
                    // `myOEImages.Find` raises on an unbound edge in the
                    // reference; the port skips it.
                    let Some(a_le_im) = self.my_oe_images.seek(a_e_or).map(|l| l.to_vec()) else {
                        continue;
                    };
                    for a_e_im in &a_le_im {
                        if a_mfence.contains(a_e_im) {
                            a_me_alone.add(a_e_im);
                        }
                    }
                }
            }
        }

        // `cxx:6164-6187`: the invalidities of all faces.
        let mut a_m_all_invs = ShapeIndexedMap::new();
        let a_nb_inv = self.my_invalid_faces.extent();
        for k in 1..=a_nb_inv {
            let Some(a_f) = self.my_invalid_faces.keys().get(k - 1).cloned() else {
                continue;
            };
            let a_lf: Vec<TopoShape> = self
                .my_invalid_faces
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for a_f_im in &a_lf {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    if self.my_invalid_edges.contains(&a_e) || a_me_alone.contains(&a_e) {
                        a_m_all_invs.add(&a_e);
                        for a_v in children_keep_orientation(&a_e) {
                            a_m_all_invs.add(&a_v);
                        }
                    }
                }
            }
        }

        // `cxx:6190-6203`.
        let mut a_mv_bounds = ShapeIndexedMap::new();
        let mut a_me_check_ext = ShapeIndexedMap::new();
        let mut a_dmeet_trim: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut a_e_images: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut a_dmoen_edges = Dmvv::new();

        // `cxx:6205-6552`.
        let a_nb_inv = self.my_invalid_faces.extent();
        for k in 1..=a_nb_inv {
            let Some(a_f_inv) = self.my_invalid_faces.keys().get(k - 1).cloned() else {
                continue;
            };
            // `cxx:6214-6216`.
            let b_self_reb_avoid = self.my_f_self_reb_avoid.contains(&a_f_inv);
            let a_lf_inv: Vec<TopoShape> = self
                .my_invalid_faces
                .seek(&a_f_inv)
                .map(|l| l.to_vec())
                .unwrap_or_default();

            // `cxx:6219-6240`: the blocks of the invalid splits.
            let a_lcb: Vec<TopoShape> = if a_lf_inv.len() > 1 {
                let a_cf_inv = make_compound_from(a_lf_inv.iter());
                AlgoTools::make_connexity_blocks_in(&a_cf_inv, ShapeType::Edge, ShapeType::Face)
            } else {
                a_lf_inv.clone()
            };

            // `cxx:6243`.
            let b_artificial = self.my_art_invalid_faces.contains_key(&shape_map_key(&a_f_inv));

            // `cxx:6245-6551`.
            for a_cb_inv in &a_lcb {
                // `cxx:6252-6259`.
                let mut a_me_fence = ShapeIndexedMap::new();
                let mut a_cbe_edges: Vec<TopoShape> = Vec::new();
                let mut an_inside_edges = ShapeIndexedMap::new();
                let mut an_inside_vertices = ShapeIndexedMap::new();

                // `cxx:6260-6280`.
                for a_e in explorer_shapes(a_cb_inv, ShapeType::Edge) {
                    if self.my_invalid_edges.contains(&a_e)
                        || (b_artificial && a_me_alone.contains(&a_e))
                    {
                        if a_me_fence.add(&a_e) == a_me_fence.extent() {
                            a_cbe_edges.push(a_e.clone());
                            if !self.my_edges_to_avoid.contains(&a_e)
                                && self.my_invalid_edges.contains(&a_e)
                            {
                                an_inside_edges.add(&a_e);
                                for a_v in children_keep_orientation(&a_e) {
                                    an_inside_vertices.add(&a_v);
                                }
                            }
                        }
                    }
                }

                // `cxx:6283-6288`.
                let a_cbe = make_compound_from(a_cbe_edges.iter());
                let a_lcbe =
                    AlgoTools::make_connexity_blocks_in(&a_cbe, ShapeType::Vertex, ShapeType::Edge);

                // `cxx:6289-6550`.
                for a_cbe_loc in &a_lcbe {
                    // `cxx:6294-6303`: the edges and vertices of the block;
                    // `aMECV` starts as the edge set.
                    let mut a_me = ShapeIndexedMap::new();
                    for a_e in explorer_shapes(a_cbe_loc, ShapeType::Edge) {
                        a_me.add(&a_e);
                    }
                    let mut a_mecv = a_me.clone();
                    for a_v in explorer_shapes(a_cbe_loc, ShapeType::Vertex) {
                        a_me.add(&a_v);
                    }

                    // `cxx:6304-6309`.
                    let mut a_mf_int = ShapeIndexedMap::new();
                    let mut a_mf_int_ext = ShapeIndexedMap::new();
                    let mut a_lf_int: Vec<TopoShape> = Vec::new();
                    let mut a_mf_avoid = ShapeIndexedMap::new();

                    // `cxx:6311-6319`.
                    self.find_faces_for_intersection(
                        &a_f_inv,
                        &a_me,
                        &a_dmsf,
                        &a_mv_inv_all,
                        b_artificial,
                        &mut a_mf_avoid,
                        &mut a_mf_int,
                        &mut a_mf_int_ext,
                        &mut a_lf_int,
                    );

                    // `cxx:6320-6325`: nothing to intersect.
                    if a_mf_int.extent() < 3 {
                        continue;
                    }

                    // `cxx:6327-6331`.
                    let p_mf_inter = self
                        .my_intersection_pairs
                        .get(&shape_map_key(&a_f_inv))
                        .cloned();

                    // `cxx:6334`.
                    let mut a_me_to_int = ShapeIndexedMap::new();
                    let a_nb = a_mf_int.extent();

                    // `cxx:6336-6519`.
                    for i in 1..=a_nb {
                        let Some(a_fi) = a_mf_int.value(i).cloned() else {
                            continue;
                        };
                        // `cxx:6342-6346`.
                        if b_self_reb_avoid && is_same(&a_fi, &a_f_inv) {
                            continue;
                        }
                        // `cxx:6348-6352`.
                        let Some(a_lf_imi) = self.my_of_images.seek(&a_fi).map(|l| l.to_vec()) else {
                            continue;
                        };
                        // `cxx:6354-6358`.
                        if !a_fle.contains(&a_fi) {
                            continue;
                        }
                        let mut a_lfei: Vec<TopoShape> = a_fle
                            .seek(&a_fi)
                            .map(|l| l.to_vec())
                            .unwrap_or_default();
                        // `cxx:6360`.
                        let mut a_lf_done: Vec<TopoShape> =
                            a_m_done.get(&shape_map_key(&a_fi)).cloned().unwrap_or_default();

                        // `cxx:6362-6367`.
                        let p_inter_fi: Option<ShapeIndexedMap> = match &p_mf_inter {
                            None => None,
                            Some(a_map) => a_map.get(&shape_map_key(&a_fi)).cloned(),
                        };
                        if p_mf_inter.is_some() && p_inter_fi.is_none() {
                            continue;
                        }

                        // `cxx:6369-6376`: the edges and vertices of the
                        // first face's splits.
                        let mut a_mev_im = ShapeIndexedMap::new();
                        map_shapes(a_lf_imi.iter(), ShapeType::Edge, &mut a_mev_im);
                        map_shapes(a_lf_imi.iter(), ShapeType::Vertex, &mut a_mev_im);
                        let is_i_contains_e = has_intersection(&a_mev_im, &an_inside_edges);
                        let is_i_contains_v = has_intersection(&a_mev_im, &an_inside_vertices);

                        // `cxx:6378-6518`.
                        for j in (i + 1)..=a_nb {
                            let Some(a_fj) = a_mf_int.value(j).cloned() else {
                                continue;
                            };
                            // `cxx:6379-6383`.
                            if b_self_reb_avoid && is_same(&a_fj, &a_f_inv) {
                                continue;
                            }
                            // `cxx:6385-6388`.
                            if let Some(a_inter) = &p_inter_fi {
                                if !a_inter.contains(&a_fj) {
                                    continue;
                                }
                            }
                            // `cxx:6390-6395`.
                            let Some(a_lf_imj) = self.my_of_images.seek(&a_fj).map(|l| l.to_vec())
                            else {
                                continue;
                            };
                            // `cxx:6397-6401`.
                            if !a_fle.contains(&a_fj) {
                                continue;
                            }
                            let mut a_lfej: Vec<TopoShape> = a_fle
                                .seek(&a_fj)
                                .map(|l| l.to_vec())
                                .unwrap_or_default();

                            // `cxx:6404-6410`.
                            a_mev_im.clear();
                            map_shapes(a_lf_imj.iter(), ShapeType::Edge, &mut a_mev_im);
                            map_shapes(a_lf_imj.iter(), ShapeType::Vertex, &mut a_mev_im);
                            let is_j_contains_e = has_intersection(&a_mev_im, &an_inside_edges);
                            let is_j_contains_v = has_intersection(&a_mev_im, &an_inside_vertices);

                            // `cxx:6412-6424`: a face connected to an inside
                            // edge needs the other one connected to an inside
                            // vertex, and the other way round; when that fails
                            // the splits must already share a vertex.
                            if (is_i_contains_e && !is_j_contains_v)
                                || (is_j_contains_e && !is_i_contains_v)
                            {
                                let mut a_lvc: Vec<TopoShape> = Vec::new();
                                find_common_parts(
                                    &a_lf_imi,
                                    &a_lf_imj,
                                    &mut a_lvc,
                                    ShapeType::Vertex,
                                );
                                if a_lvc.is_empty() {
                                    continue;
                                }
                            }

                            // `cxx:6426-6465`: reuse the common edges instead of
                            // intersecting again.
                            let mut a_lec: Vec<TopoShape> = Vec::new();
                            find_common_parts(
                                &a_lf_imi,
                                &a_lf_imj,
                                &mut a_lec,
                                ShapeType::Edge,
                            );
                            if !a_lec.is_empty() {
                                let b_force_use = a_mf_int_ext.contains(&a_fi)
                                    || a_mf_int_ext.contains(&a_fj);
                                self.process_common_edges(
                                    &a_lec,
                                    &a_me,
                                    &a_me_inf_e_trim,
                                    &a_m_all_invs,
                                    b_force_use,
                                    &mut a_mecv,
                                    &mut a_me_check_ext,
                                    &mut a_dmeet_trim,
                                    &mut a_lfei,
                                    &mut a_lfej,
                                    &mut a_me_to_int,
                                );

                                // `cxx:6447-6464`: the common vertices that do
                                // not belong to the common edges are used for
                                // trimming.
                                let mut a_mv_on_ce = ShapeIndexedMap::new();
                                for it_e in &a_lec {
                                    for a_v in explorer_shapes(it_e, ShapeType::Vertex) {
                                        a_mv_on_ce.add(&a_v);
                                    }
                                }
                                let mut a_lev: Vec<TopoShape> = Vec::new();
                                find_common_parts(
                                    &a_lf_imi,
                                    &a_lf_imj,
                                    &mut a_lev,
                                    ShapeType::Vertex,
                                );
                                for it_v in &a_lev {
                                    if !a_mv_on_ce.contains(it_v) {
                                        a_mecv.add(it_v);
                                    }
                                }
                                // `cxx:6463-6465`: the local copies go back into
                                // the maps `ProcessCommonEdges` grew.
                                a_fle.bind_list(&a_fi, &a_lfei);
                                a_fle.bind_list(&a_fj, &a_lfej);
                                continue;
                            }

                            // `cxx:6467-6471`: two non-artificial invalid
                            // faces sharing edges are left alone.
                            if self.my_invalid_faces.contains(&a_fi)
                                && self.my_invalid_faces.contains(&a_fj)
                                && !self
                                    .my_art_invalid_faces
                                    .contains_key(&shape_map_key(&a_fi))
                                && !self
                                    .my_art_invalid_faces
                                    .contains_key(&shape_map_key(&a_fj))
                            {
                                continue;
                            }

                            // `cxx:6473-6499`: the pair was already treated -
                            // reuse its intersection line and refine the
                            // origins.
                            let b_done = a_lf_done.iter().any(|a_f| is_same(a_f, &a_fj));
                            if b_done {
                                self.update_intersected_faces(
                                    &a_f_inv,
                                    &a_fi,
                                    &a_fj,
                                    &a_lf_inv,
                                    &a_lf_imi,
                                    &a_lf_imj,
                                    &a_lfei,
                                    &a_lfej,
                                    &mut a_me_to_int,
                                );
                                a_fle.bind_list(&a_fi, &a_lfei);
                                a_fle.bind_list(&a_fj, &a_lfej);
                                continue;
                            }

                            // `cxx:6500-6504`.
                            if a_mf_avoid.contains(&a_fi) || a_mf_avoid.contains(&a_fj) {
                                continue;
                            }

                            // `cxx:6506-6507`: the pair is marked before the
                            // intersection so the second direction is not
                            // intersected again.
                            append_to_list(&mut a_lf_done, &a_fj);
                            a_m_done.insert(shape_map_key(&a_fi), a_lf_done.clone());
                            a_m_done
                                .entry(shape_map_key(&a_fj))
                                .or_default()
                                .push(a_fi.clone());

                            // `cxx:6509-6518`.
                            self.intersect_faces_of_pair(
                                &a_f_inv,
                                &a_fi,
                                &a_fj,
                                &a_lf_inv,
                                &a_lf_imi,
                                &a_lf_imj,
                                &mut a_lfei,
                                &mut a_lfej,
                                &mut a_mecv,
                                &mut a_me_to_int,
                            );
                            a_fle.bind_list(&a_fi, &a_lfei);
                            a_fle.bind_list(&a_fj, &a_lfej);
                        }

                        // `cxx:6521-6532`: trim the intersection edges of this
                        // chain.
                        let a_ss_interfs = if b_artificial {
                            Some(&self.my_ss_interfs_art)
                        } else {
                            None
                        };
                        self.intersect_and_trim_edges(
                            &a_mf_int,
                            &a_me_to_int,
                            &a_dmeet_trim,
                            &a_me,
                            &a_mecv,
                            &a_mv_inv,
                            &a_mvr_inv,
                            &a_me_check_ext,
                            a_ss_interfs,
                            &mut a_mv_bounds,
                            &mut a_e_images,
                        );

                        // `cxx:6534-6550`: every edge of the block is linked to
                        // the edges to intersect.
                        let a_nb_e_to_int = a_me_to_int.extent();
                        for i_e in 1..=a_nb_e_to_int {
                            let Some(a_e_int) = a_me_to_int.value(i_e).cloned() else {
                                continue;
                            };
                            for a_e in explorer_shapes(a_cbe_loc, ShapeType::Edge) {
                                if !a_dmoen_edges.contains(&a_e) {
                                    a_dmoen_edges.bind_list(&a_e, &[]);
                                }
                                a_dmoen_edges.add(&a_e, &a_e_int);
                            }
                        }
                    }
                }
            }
        }

        // `cxx:6556-6564`: filter the obtained edges.
        self.update_valid_edges(
            a_fle,
            a_dmoen_edges,
            a_mv_bounds,
            a_me_inv_on_art,
            a_me_check_ext,
            the_verts_to_avoid,
            a_e_images,
            a_dmeet_trim,
        );
    }

}

/// `NCollection_MapAlgo::HasIntersection(theMap1, theMap2)` (`:6374`,
/// `:6376`, `:6406`, `:6408`): whether the two sets share an element. The
/// reference's second set is an `NCollection_Map`, the first one is an
/// `IndexedMap`; both test `IsSame`.
fn has_intersection(the_map1: &ShapeIndexedMap, the_map2: &ShapeIndexedMap) -> bool {
    the_map1.iter().any(|s| the_map2.contains(s))
}
