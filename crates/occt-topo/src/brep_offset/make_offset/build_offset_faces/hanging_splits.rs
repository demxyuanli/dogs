//! `RemoveHangingParts` (`cxx:4962-5164`), `RemoveValidSplits`
//! (`cxx:5165-5223`) and `RemoveInvalidSplits` (`cxx:5225-5315`).
//!
//! `RemoveHangingParts` finds the face blocks of the volume builder that no
//! solid uses and drops the ones that are isolated invalid material
//! (`cxx:5005-5147`). `RemoveValidSplits` / `RemoveInvalidSplits` then empty
//! the `myOFImages` / `myInvalidFaces` lists of the splits marked for removal,
//! collecting their edges in `theMERemoved` so `RemoveInsideFaces` can
//! classify them (`cxx:4619-4627`).

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::bop_maker_volume::BopMakerVolume;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, map_shapes_and_ancestors, shape_map_key, ShapeMapKey};

use super::invalid_edges_filter::make_compound_from;
use super::statics::take_modified;
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::RemoveHangingParts(const
    /// BOPAlgo_MakerVolume& theMV, const DataMap<Shape, Shape>& theDMFImF,
    /// const IndexedMap& theMFInv, Map& theMFToRem)` (`cxx:4962-5164`).
    ///
    /// A block is removed when it contains an inverted edge that no solid uses
    /// (`cxx:5046-5061`), or when it holds an invalid split and is either not
    /// connected to an invalid part of the solids or comes from a single
    /// offset face (`cxx:5143-5148`).
    pub fn remove_hanging_parts(
        &self,
        the_mv: &BopMakerVolume,
        the_dmf_im_f: &HashMap<ShapeMapKey, TopoShape>,
        the_mf_inv: &ShapeIndexedMap,
        the_mf_to_rem: &mut ShapeIndexedMap,
    ) {
        // `cxx:4976-4978`: the faces of the result solids, used as a fence.
        let mut a_mf_s = ShapeIndexedMap::new();
        for a_f in explorer_shapes(the_mv.shape(), ShapeType::Face) {
            a_mf_s.add(&a_f);
        }

        let mut a_cf_hangs_splits: Vec<TopoShape> = Vec::new();

        // `cxx:4984-4986`.
        let a_mv_ims = the_mv.images();

        // `cxx:4988-4997`: every split of the arguments that is not on a solid.
        for a_arg in the_mv.arguments() {
            for a_f in explorer_shapes(a_arg, ShapeType::Face) {
                take_modified(&a_f, a_mv_ims, &mut a_cf_hangs_splits, Some(&mut a_mf_s));
            }
        }

        // `cxx:4999-5005`.
        let a_cf_hangs = make_compound_from(a_cf_hangs_splits.iter());
        let a_lcb_hangs =
            AlgoTools::make_connexity_blocks_in(&a_cf_hangs, ShapeType::Edge, ShapeType::Face);
        if a_lcb_hangs.is_empty() {
            return;
        }

        // `cxx:5012-5018`: the result edges/vertices and the faces around them.
        let a_dm_ef = map_shapes_and_ancestors(the_mv.shape(), ShapeType::Edge, ShapeType::Face);
        let a_dm_ve = map_shapes_and_ancestors(the_mv.shape(), ShapeType::Vertex, ShapeType::Edge);

        // `cxx:5020-5033`: the invalid and inverted edges with their images.
        let mut a_me_inv = ShapeIndexedMap::new();
        let a_nb_e = self.my_invalid_edges.extent();
        for i in 1..=a_nb_e {
            let Some(a_e) = self.my_invalid_edges.value(i).cloned() else {
                continue;
            };
            take_modified(&a_e, a_mv_ims, &mut a_me_inv, None);
        }
        let mut a_me_inverted = ShapeIndexedMap::new();
        let a_nb_inv = self.my_inverted_edges.extent();
        for i in 1..=a_nb_inv {
            let Some(a_e) = self.my_inverted_edges.value(i).cloned() else {
                continue;
            };
            take_modified(&a_e, a_mv_ims, &mut a_me_inverted, None);
        }

        // `cxx:5036-5038`.
        let a_mv_ors = the_mv.origins();

        // `cxx:5040-5148`.
        let mut a_blocks_to_remove: Vec<TopoShape> = Vec::new();

        for a_cbh in &a_lcb_hangs {
            // `cxx:5046-5061`: a block with an inverted edge no solid uses.
            let mut b_has_inverted = false;
            for a_e in explorer_shapes(a_cbh, ShapeType::Edge) {
                b_has_inverted =
                    !a_dm_ef.contains_key(&shape_map_key(&a_e)) && a_me_inverted.contains(&a_e);
                if b_has_inverted {
                    break;
                }
            }
            if b_has_inverted {
                a_blocks_to_remove.push(a_cbh.clone());
                continue;
            }

            // `cxx:5063-5068`.
            let mut b_has_invalid_face = false;
            let mut b_is_connected = false;
            let mut a_block_me = ShapeIndexedMap::new();
            for a_e in explorer_shapes(a_cbh, ShapeType::Edge) {
                a_block_me.add(&a_e);
            }
            let mut a_m_offset_f = ShapeIndexedMap::new();

            // `cxx:5071-5141`.
            for a_f in explorer_shapes(a_cbh, ShapeType::Face) {
                // `cxx:5076-5079`.
                if !b_has_invalid_face {
                    b_has_invalid_face = the_mf_inv.contains(&a_f);
                }

                // `cxx:5081-5117`.
                if !b_is_connected {
                    // `cxx:5086-5097`: through the edges of the split.
                    for a_e in explorer_shapes(&a_f, ShapeType::Edge) {
                        if let Some(p_lf) = a_dm_ef.get(&shape_map_key(&a_e)) {
                            for a_f_around in p_lf {
                                if the_mf_inv.contains(a_f_around) {
                                    b_is_connected = true;
                                    break;
                                }
                            }
                        }
                        if b_is_connected {
                            break;
                        }
                    }
                    // `cxx:5100-5116`: through the vertices.
                    if !b_is_connected {
                        for a_v in explorer_shapes(&a_f, ShapeType::Vertex) {
                            if let Some(p_le) = a_dm_ve.get(&shape_map_key(&a_v)) {
                                for a_e_around in p_le {
                                    if !a_block_me.contains(a_e_around)
                                        && a_me_inv.contains(a_e_around)
                                    {
                                        b_is_connected = true;
                                        break;
                                    }
                                }
                            }
                            if b_is_connected {
                                break;
                            }
                        }
                    }
                }

                // `cxx:5119-5141`: the offset faces the split comes from.
                match a_mv_ors.get(&shape_map_key(&a_f)) {
                    Some(p_lf_or) => {
                        for a_or in p_lf_or.clone() {
                            if let Some(p_f_offset) = the_dmf_im_f.get(&shape_map_key(&a_or)) {
                                a_m_offset_f.add(p_f_offset);
                            }
                        }
                    }
                    None => {
                        if let Some(p_f_offset) = the_dmf_im_f.get(&shape_map_key(&a_f)) {
                            a_m_offset_f.add(p_f_offset);
                        }
                    }
                }
            }

            // `cxx:5143-5148`.
            let b_remove = b_has_invalid_face && (!b_is_connected || a_m_offset_f.extent() == 1);
            if b_remove {
                a_blocks_to_remove.push(a_cbh.clone());
            }
        }

        // `cxx:5150-5162`.
        for a_cbh in &a_blocks_to_remove {
            for a_f in explorer_shapes(a_cbh, ShapeType::Face) {
                the_mf_to_rem.add(&a_f);
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::RemoveValidSplits(const Map&
    /// theSpRem, BOPAlgo_Builder& theGF, IndexedMap& theMERemoved)`
    /// (`cxx:5165-5223`).
    ///
    /// A valid split goes when it is marked for removal, or when every one of
    /// its images in `theGF` is (`cxx:5190-5216`). Both arms record the edges of
    /// the removed split.
    pub fn remove_valid_splits(
        &mut self,
        the_sp_rem: &ShapeIndexedMap,
        the_gf: &BopMakerVolume,
        the_me_removed: &mut ShapeIndexedMap,
    ) {
        let a_nb = self.my_of_images.extent();
        // `cxx:5175-5178`.
        if a_nb == 0 {
            return;
        }

        for i in 1..=a_nb {
            let Some(a_f) = self.my_of_images.keys().get(i - 1).cloned() else {
                continue;
            };
            let Some(a_ls_im) = self.my_of_images.change_seek(&a_f) else {
                continue;
            };
            let mut j = 0usize;
            while j < a_ls_im.len() {
                let a_s_im = a_ls_im[j].clone();

                // `cxx:5186-5194`.
                if the_sp_rem.contains(&a_s_im) {
                    for a_e in explorer_shapes(&a_s_im, ShapeType::Edge) {
                        the_me_removed.add(&a_e);
                    }
                    a_ls_im.remove(j);
                    continue;
                }

                // `cxx:5196-5217`: the images of the split.
                let a_ls_im_im = the_gf.modified(&a_s_im);
                if !a_ls_im_im.is_empty() {
                    let mut b_all_rem = true;
                    for a_s_im_im in &a_ls_im_im {
                        if the_sp_rem.contains(a_s_im_im) {
                            for a_e in explorer_shapes(a_s_im_im, ShapeType::Edge) {
                                the_me_removed.add(&a_e);
                            }
                        } else {
                            b_all_rem = false;
                        }
                    }
                    // `cxx:5211-5217`.
                    if b_all_rem {
                        for a_e in explorer_shapes(&a_s_im, ShapeType::Edge) {
                            the_me_removed.add(&a_e);
                        }
                        a_ls_im.remove(j);
                        continue;
                    }
                }

                j += 1;
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::RemoveInvalidSplits(const Map&
    /// theSpRem, BOPAlgo_Builder& theGF, IndexedMap& theMERemoved)`
    /// (`cxx:5225-5315`).
    ///
    /// Unlike [`Self::remove_valid_splits`], a split removed because it is
    /// marked does *not* record its edges in `theMERemoved` (`cxx:5249-5254`);
    /// only the "every image removed" (`cxx:5280-5284`) and "every invalid
    /// edge removed" (`cxx:5306-5308`) arms do.
    pub fn remove_invalid_splits(
        &mut self,
        the_sp_rem: &ShapeIndexedMap,
        the_gf: &BopMakerVolume,
        the_me_removed: &mut ShapeIndexedMap,
    ) {
        let a_nb = self.my_invalid_faces.extent();
        // `cxx:5233-5236`.
        if a_nb == 0 {
            return;
        }

        for i in 1..=a_nb {
            let Some(a_s) = self.my_invalid_faces.keys().get(i - 1).cloned() else {
                continue;
            };
            // `cxx:5241`.
            let b_art = self
                .my_art_invalid_faces
                .contains_key(&shape_map_key(&a_s));

            let Some(a_ls_im) = self.my_invalid_faces.change_seek(&a_s) else {
                continue;
            };
            let mut j = 0usize;
            while j < a_ls_im.len() {
                let a_s_im = a_ls_im[j].clone();

                // `cxx:5248-5255`.
                if the_sp_rem.contains(&a_s_im) {
                    a_ls_im.remove(j);
                    continue;
                }

                // `cxx:5257-5263`.
                let a_ls_im_im = the_gf.modified(&a_s_im);
                if a_ls_im_im.is_empty() {
                    j += 1;
                    continue;
                }

                // `cxx:5265-5284`.
                let mut b_all_rem = true;
                let mut a_me_removed = ShapeIndexedMap::new();
                for a_s_im_im in &a_ls_im_im {
                    if the_sp_rem.contains(a_s_im_im) {
                        for a_e in explorer_shapes(a_s_im_im, ShapeType::Edge) {
                            a_me_removed.add(&a_e);
                        }
                    } else {
                        b_all_rem = false;
                    }
                }
                if b_all_rem {
                    a_ls_im.remove(j);
                    continue;
                }

                // `cxx:5287-5291`.
                if b_art {
                    j += 1;
                    continue;
                }

                // `cxx:5293-5313`: all invalid edges of the split are removed.
                let mut b_has_invalid_left = false;
                for a_e_inv in explorer_shapes(&a_s_im, ShapeType::Edge) {
                    if self.my_invalid_edges.contains(&a_e_inv)
                        && !a_me_removed.contains(&a_e_inv)
                    {
                        b_has_invalid_left = true;
                        break;
                    }
                }
                if !b_has_invalid_left {
                    for a_e in explorer_shapes(&a_s_im, ShapeType::Edge) {
                        the_me_removed.add(&a_e);
                    }
                    a_ls_im.remove(j);
                } else {
                    j += 1;
                }
            }
        }
    }
}
