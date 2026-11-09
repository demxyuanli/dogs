//! Face-level validity bookkeeping of `BRepOffset_BuildOffsetFaces` -
//! `FindInvalidFaces` (`cxx:2856-3085`) and `FilterInvalidFaces`
//! (`cxx:5369-5549`).
//!
//! `FindInvalidFaces` is the reading of the edge classification produced by the
//! two `FindInvalidEdges` passes, per split of one offset face. The reference's
//! comment at `cxx:2870-2880` states the contract:
//!
//! * a face is invalid when it has been reverted (at least two unconnected
//!   edges changed orientation) or when every checked edge is invalid;
//! * a face is removed from the splits when all its checked edges are invalid
//!   for this face but valid for another split;
//! * a face is kept when some of its edges are valid for it.
//!
//! `FilterInvalidFaces` is the second-thought pass run after the invalid splits
//! have been dropped: it re-validates a face whose images are gone and refuses
//! to remove a face whose removal would leave a free edge.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::Dmvv;
use crate::brep_tools::outer_wire;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{
    explorer_shapes, is_same, map_shapes_and_ancestors_into, shape_map_key, ShapeMapKey,
};

use super::invalid_edges_filter::make_compound_from;
use super::statics::append_to_list;
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::FindInvalidFaces(List& theLFImages,
    /// const DataMap& theDMFMVE, const BRepOffset_DataMapOfShapeIndexedMap&
    /// theDMFMIE, const Map& theMENeutral, const Map&
    /// theEdgesInvalidByVertex, const Map& theEdgesValidByVertex, const Map&
    /// theMFHoles, IndexedMap& theMFInvInHole, List& theInvFaces, List&
    /// theInvertedFaces)` (`cxx:2856-3085`).
    ///
    /// `the_lf_images` is the reference's `NCollection_List` of splits: the
    /// method both walks it and removes entries from it (`cxx:2995`, `:3003`,
    /// `:3015`), which the `Vec` does by index so the removal keeps the
    /// reference's iterator position (the element after the removed one).
    pub fn find_invalid_faces(
        &mut self,
        the_lf_images: &mut Vec<TopoShape>,
        the_dmfmve: &HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_dmfmie: &HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_me_neutral: &ShapeIndexedMap,
        the_edges_invalid_by_vertex: &ShapeIndexedMap,
        the_edges_valid_by_vertex: &ShapeIndexedMap,
        the_mf_holes: &ShapeIndexedMap,
        the_mf_inv_in_hole: &mut ShapeIndexedMap,
        the_inv_faces: &mut Vec<TopoShape>,
        the_inverted_faces: &mut Vec<TopoShape>,
    ) {
        // `cxx:2885`: a single split is treated as "inverted is invalid".
        let b_treat_inverted_as_invalid = the_lf_images.len() == 1;

        // `cxx:2888`: neutral edges to remove.
        let mut a_men_rem = ShapeIndexedMap::new();
        // `cxx:2891`: faces for the post treat.
        let mut a_lfpt: Vec<TopoShape> = Vec::new();

        // `cxx:2894-2900`: edge -> split over the images.
        //
        // The reference fills `aDMEF` with
        // `MapShapesAndAncestors(aFIm, EDGE, FACE, aDMEF)` and then appends
        // through `AppendToList` (`cxx:3726-3736`), which keeps one entry per
        // `IsSame` split. Walking the edges of each split gives the same key
        // set and the same list lengths.
        let mut a_dmef = Dmvv::new();
        for a_f_im in the_lf_images.iter() {
            for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                if a_dmef.seek(&a_e).is_none() {
                    a_dmef.bind_list(&a_e, &[]);
                }
                if let Some(a_lf) = a_dmef.change_seek(&a_e) {
                    append_to_list(a_lf, a_f_im);
                }
            }
        }

        // `cxx:2902-3036`: the first pass.
        let mut i = 0usize;
        while i < the_lf_images.len() {
            let Some(a_f_im) = Face::wrap(the_lf_images[i].clone()) else {
                i += 1;
                continue;
            };

            // `cxx:2908-2911`: the valid / invalid edges recorded for the split.
            let p_mve = the_dmfmve.get(&shape_map_key(&a_f_im.0));
            let p_mie = the_dmfmie.get(&shape_map_key(&a_f_im.0));

            let mut b_has_valid = false;
            let mut b_all_valid = true;
            let mut b_all_invalid = true;
            let mut b_has_really_invalid = false;
            let mut b_all_inv_neutral = true;
            let mut b_is_invalid_by_inverted = true;
            let mut b_has_inverted = false;
            let mut a_nb_checked = 0usize;

            // `cxx:2922-2973`: the outer wire of the split.
            let a_w_im = outer_wire(&a_f_im);
            let a_w_edges = match &a_w_im {
                Some(w) => crate::topo_tools::edges_of_wire(w)
                    .into_iter()
                    .map(|e| e.0)
                    .collect(),
                None => Vec::new(),
            };
            for a_e_im in &a_w_edges {
                let b_valid = self.my_valid_edges.contains(a_e_im);
                let mut b_invalid = self.my_invalid_edges.contains(a_e_im);
                let b_neutral = the_me_neutral.contains(a_e_im);

                // `cxx:2932-2936`: an edge nobody checked.
                if !b_valid && !b_invalid && !b_neutral {
                    continue;
                }

                // `cxx:2939-2943`: skip non-boundary edges originating from a
                // vertex.
                if (the_edges_invalid_by_vertex.contains(a_e_im)
                    || the_edges_valid_by_vertex.contains(a_e_im))
                    && a_dmef.seek(a_e_im).map(|l| l.len()).unwrap_or(0) != 1
                {
                    continue;
                }

                a_nb_checked += 1;

                // `cxx:2947-2951`.
                let b_invalid_loc = p_mie.map(|m| m.contains(a_e_im)).unwrap_or(false);
                b_has_really_invalid = b_invalid
                    && b_invalid_loc
                    && !b_valid
                    && !the_edges_invalid_by_vertex.contains(a_e_im);
                if b_has_really_invalid {
                    break;
                }

                // `cxx:2955-2961`.
                let b_valid_loc = p_mve.map(|m| m.contains(a_e_im)).unwrap_or(false);
                let b_inverted = self.my_inverted_edges.contains(a_e_im);
                if !b_invalid && !b_invalid_loc && b_treat_inverted_as_invalid {
                    b_invalid = b_inverted;
                }

                // `cxx:2963-2971`.
                if b_valid_loc && b_neutral {
                    b_has_valid = true;
                }
                b_all_valid &= b_valid_loc;
                b_all_invalid &= b_invalid || b_invalid_loc;
                b_all_inv_neutral &= b_all_invalid && b_neutral;
                b_is_invalid_by_inverted &= b_invalid_loc || b_inverted;
                b_has_inverted |= b_inverted;
            }

            // `cxx:2975-2979`.
            if a_nb_checked == 0 {
                i += 1;
                continue;
            }

            // `cxx:2980-2997`: the part filled by overlapping splits.
            if !b_has_really_invalid && (b_all_inv_neutral && !b_has_valid) && a_nb_checked > 1 {
                if b_has_inverted {
                    i += 1;
                    continue;
                }
                map_shapes_into(&a_f_im.0, ShapeType::Edge, &mut a_men_rem);
                the_lf_images.remove(i);
                continue;
            }

            // `cxx:2999-3006`: an invalid split.
            if b_has_really_invalid
                || (b_all_invalid
                    && !(b_has_valid || b_all_valid)
                    && (!b_all_inv_neutral || a_nb_checked != 1))
            {
                the_inv_faces.push(a_f_im.0.clone());
                if the_mf_holes.contains(&a_f_im.0) {
                    the_mf_inv_in_hole.add(&a_f_im.0);
                }
                i += 1;
                continue;
            }

            // `cxx:3008-3016`: a split inside a hole.
            if the_mf_holes.contains(&a_f_im.0) {
                map_shapes_into(&a_f_im.0, ShapeType::Edge, &mut a_men_rem);
                the_lf_images.remove(i);
                continue;
            }

            // `cxx:3018-3022`.
            if b_is_invalid_by_inverted && !(b_has_valid || b_all_valid) {
                the_inverted_faces.push(a_f_im.0.clone());
            }

            // `cxx:3024-3033`.
            if !b_all_inv_neutral {
                a_lfpt.push(a_f_im.0.clone());
            } else {
                map_shapes_into(&a_f_im.0, ShapeType::Edge, &mut a_men_rem);
            }
            i += 1;
        }

        // `cxx:3038-3041`.
        if a_lfpt.is_empty() || a_men_rem.extent() == 0 {
            return;
        }

        // `cxx:3043-3082`: re-check the splits with the removed neutrals gone.
        for a_f_shape in &a_lfpt {
            let Some(a_f_im) = Face::wrap(a_f_shape.clone()) else {
                continue;
            };
            let p_mve = the_dmfmve.get(&shape_map_key(&a_f_im.0));

            let mut b_has_valid = false;
            let mut b_all_valid = true;
            let mut b_all_invalid = true;

            let a_w_im = outer_wire(&a_f_im);
            let a_w_edges = match &a_w_im {
                Some(w) => crate::topo_tools::edges_of_wire(w)
                    .into_iter()
                    .map(|e| e.0)
                    .collect(),
                None => Vec::new(),
            };
            for a_e_im in &a_w_edges {
                // `cxx:3059`: the reference computes `bValid` here and never
                // reads it (only `bValidLoc` feeds the flags below).
                let _b_valid = self.my_valid_edges.contains(a_e_im);
                let mut b_invalid = self.my_invalid_edges.contains(a_e_im);
                let b_neutral = the_me_neutral.contains(a_e_im) && !a_men_rem.contains(a_e_im);
                let b_valid_loc = p_mve.map(|m| m.contains(a_e_im)).unwrap_or(false);

                // `cxx:3064-3067`.
                if !b_invalid && b_treat_inverted_as_invalid {
                    b_invalid = self.my_inverted_edges.contains(a_e_im);
                }

                // `cxx:3069-3073`.
                if b_valid_loc && b_neutral {
                    b_has_valid = true;
                }
                b_all_valid &= b_valid_loc;
                b_all_invalid &= b_invalid;
            }

            // `cxx:3078-3081`.
            if b_all_invalid && !b_has_valid && !b_all_valid {
                the_inv_faces.push(a_f_im.0.clone());
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::FilterInvalidFaces(const
    /// IndexedDataMap<Shape, List>& theDMEF, const IndexedMap&
    /// theMERemoved)` (`cxx:5369-5549`).
    ///
    /// Rebuilds `myInvalidFaces` from the splits that are still invalid once
    /// the removed edges are known: a face is kept as invalid only when its own
    /// image list is gone, or when removing it would leave a free edge, or when
    /// one of its invalid edges is shared with a split outside the face.
    pub fn filter_invalid_faces(
        &mut self,
        the_dmef: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_me_removed: &ShapeIndexedMap,
    ) {
        // `cxx:5380`: the rebuilt invalid-face map.
        let mut a_really_inv_faces = Dmvv::new();
        // `cxx:5382-5383`: edge -> face over both invalid and valid splits.
        let mut a_dmef_all: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        let a_nb = self.my_invalid_faces.extent();

        for i in 1..=a_nb {
            let Some(a_f) = self.my_invalid_faces.keys().get(i - 1).cloned() else {
                continue;
            };
            let a_lf_inv: Vec<TopoShape> = self
                .my_invalid_faces
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();

            // `cxx:5394-5407`: an artificially invalid face keeps its entry
            // only while it still had invalid splits.
            if self.my_art_invalid_faces.contains_key(&shape_map_key(&a_f)) {
                if a_lf_inv.is_empty() {
                    self.unbind_art_invalid_face(&a_f);
                } else {
                    a_really_inv_faces.bind_list(&a_f, &a_lf_inv);
                }
                continue;
            }

            // `cxx:5409-5412`.
            if a_lf_inv.is_empty() {
                continue;
            }

            // `cxx:5414`: the current splits of the face.
            let a_lf_im: Vec<TopoShape> = self
                .my_of_images
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            let mut b_invalid = a_lf_im.is_empty();

            // `cxx:5417-5440`: do the invalid splits intersect the current
            // ones?
            if !b_invalid {
                let mut b_found = false;
                for a_f_inv in &a_lf_inv {
                    for a_f_im in &a_lf_im {
                        if is_same(a_f_inv, a_f_im) {
                            b_found = true;
                            break;
                        }
                    }
                    if b_found {
                        break;
                    }
                }
                b_invalid = b_found;
            }

            // `cxx:5442-5476`: a free edge on either list keeps the face.
            if !b_invalid {
                for j in 0..2 {
                    let a_li: &[TopoShape] = if j == 0 { &a_lf_im } else { &a_lf_inv };
                    let mut b_free = false;
                    for a_f_im in a_li {
                        let mut b_this_free = false;
                        for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                            if !the_me_removed.contains(&a_e) {
                                if let Some(p_lef) = the_dmef.get(&shape_map_key(&a_e)) {
                                    if p_lef.len() == 1 {
                                        b_this_free = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if b_this_free {
                            b_free = true;
                            break;
                        }
                    }
                    if b_free {
                        b_invalid = true;
                        break;
                    }
                }
            }

            // `cxx:5477-5543`.
            if b_invalid {
                if a_dmef_all.is_empty() {
                    // `cxx:5479-5488`.
                    a_dmef_all = the_dmef.clone();
                    for i_f in 1..=a_nb {
                        let Some(a_key) = self.my_invalid_faces.keys().get(i_f - 1) else {
                            continue;
                        };
                        let a_lf: Vec<TopoShape> = self
                            .my_invalid_faces
                            .seek(a_key)
                            .map(|l| l.to_vec())
                            .unwrap_or_default();
                        for a_f_inv in &a_lf {
                            map_shapes_and_ancestors_into(
                                a_f_inv,
                                ShapeType::Edge,
                                ShapeType::Face,
                                &mut a_dmef_all,
                            );
                        }
                    }
                }

                // `cxx:5490-5497`.
                let mut a_local_splits = ShapeIndexedMap::new();
                for j in 0..2 {
                    let a_li: &[TopoShape] = if j == 0 { &a_lf_im } else { &a_lf_inv };
                    for a_f_im in a_li {
                        a_local_splits.add(a_f_im);
                    }
                }

                // `cxx:5499-5531`: every invalid edge of the face must stay
                // inside the face's own splits.
                let mut b_shared = false;
                for a_f_im in &a_lf_inv {
                    let mut b_this_shared = false;
                    for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                        if self.my_invalid_edges.contains(&a_e) && !the_me_removed.contains(&a_e) {
                            if let Some(a_lf) = a_dmef_all.get(&shape_map_key(&a_e)) {
                                if a_lf.iter().any(|f| !a_local_splits.contains(f)) {
                                    b_this_shared = true;
                                    break;
                                }
                            }
                        }
                    }
                    if b_this_shared {
                        b_shared = true;
                        break;
                    }
                }
                b_invalid = b_shared;

                // `cxx:5532-5539`: the face is valid after all - hand its
                // splits back to the face.
                if !b_invalid {
                    let mut a_lf_im = a_lf_im.clone();
                    for a_f_inv in &a_lf_inv {
                        append_to_list(&mut a_lf_im, a_f_inv);
                    }
                    self.my_of_images.bind_list(&a_f, &a_lf_im);
                }
            }

            // `cxx:5546-5549`.
            if b_invalid {
                a_really_inv_faces.bind_list(&a_f, &a_lf_inv);
            }
        }

        self.my_invalid_faces = a_really_inv_faces;
    }

    /// `void BRepOffset_BuildOffsetFaces::RemoveInvalidSplitsFromValid(const
    /// DataMap& theDMFMVIE)` (`cxx:3893-4070`).
    ///
    /// The invalid splits are grouped into connexity blocks over the original
    /// faces; a block is dropped when every one of its splits belongs to a face
    /// that has no valid image left (`cxx:3946-3974`), or when the split
    /// connected by inverted edges may go (`cxx:3976-4045`). Artificial invalid
    /// faces never take part (`cxx:3916-3920`).
    ///
    /// `the_dmfmvie` is the reference's `DataMap<Shape, Map<Shape>>` of the
    /// valid inverted edges per split (`theDMFMVIE`).
    pub fn remove_invalid_splits_from_valid(
        &mut self,
        the_dmfmvie: &HashMap<ShapeMapKey, ShapeIndexedMap>,
    ) {
        // `cxx:3904`: the fence and the splits to remove.
        let mut a_m_fence = ShapeIndexedMap::new();
        let mut a_mf_to_rem = ShapeIndexedMap::new();

        // `cxx:3911-3932`: the compound of every still-invalid split, plus
        // split -> original face.
        let mut a_ce_inv_edges: Vec<TopoShape> = Vec::new();
        let mut a_dmifof: HashMap<ShapeMapKey, TopoShape> = HashMap::new();

        let a_nb = self.my_invalid_faces.extent();
        for i in 1..=a_nb {
            let Some(a_f) = self.my_invalid_faces.keys().get(i - 1).cloned() else {
                continue;
            };
            // `cxx:3917-3920`: artificial invalid faces are not removed.
            if self.my_art_invalid_faces.contains_key(&shape_map_key(&a_f)) {
                continue;
            }
            let a_lf_inv: Vec<TopoShape> = self
                .my_invalid_faces
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for a_f_im in &a_lf_inv {
                if a_m_fence.add(a_f_im) == a_m_fence.extent() {
                    a_ce_inv_edges.push(a_f_im.clone());
                    a_dmifof.insert(shape_map_key(a_f_im), a_f.clone());
                }
            }
        }

        // `cxx:3935-3937`: the connexity blocks of the invalid splits.
        let a_cf_inv = make_compound_from(a_ce_inv_edges.iter());
        let a_lcb_inv = AlgoTools::make_connexity_blocks_in(
            &a_cf_inv,
            ShapeType::Edge,
            ShapeType::Face,
        );

        // `cxx:3939-4046`: one block at a time.
        for a_cb in &a_lcb_inv {
            let a_cb_faces = explorer_shapes(a_cb, ShapeType::Face);

            // `cxx:3944-3947`: a single-split block is removable; otherwise
            // look for a split whose face has no valid image left.
            let mut b_removable = a_cb_faces.len() < 2;
            if !b_removable {
                // `cxx:3950-3962`.
                for a_f_im in &a_cb_faces {
                    let Some(a_f) = a_dmifof.get(&shape_map_key(a_f_im)) else {
                        continue;
                    };
                    let a_lf_im_len = self.my_of_images.seek(a_f).map(|l| l.len()).unwrap_or(0);
                    let a_lf_inv_len = self
                        .my_invalid_faces
                        .seek(a_f)
                        .map(|l| l.len())
                        .unwrap_or(0);
                    if a_lf_im_len == a_lf_inv_len {
                        b_removable = true;
                        break;
                    }
                }
            }

            // `cxx:3966-3974`.
            if b_removable {
                for a_f in &a_cb_faces {
                    a_mf_to_rem.add(a_f);
                }
                continue;
            }

            // `cxx:3977-3997`: edge -> split over the block, and original face
            // -> the block's splits of it.
            let mut a_dmef: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
            map_shapes_and_ancestors_into(
                a_cb,
                ShapeType::Edge,
                ShapeType::Face,
                &mut a_dmef,
            );
            let mut a_dmff: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
            let mut a_dmff_order: Vec<TopoShape> = Vec::new();
            for a_f_cb in &a_cb_faces {
                let Some(a_f) = a_dmifof.get(&shape_map_key(a_f_cb)) else {
                    continue;
                };
                let a_key = shape_map_key(a_f);
                if !a_dmff.contains_key(&a_key) {
                    a_dmff_order.push(a_f.clone());
                }
                let a_lf = a_dmff.entry(a_key).or_default();
                a_lf.push(a_f_cb.clone());
            }

            // `cxx:3999-4044`: either remove all the splits of one original
            // face or none of them.
            for a_f in &a_dmff_order {
                let a_lf_cb = match a_dmff.get(&shape_map_key(a_f)) {
                    Some(l) => l.clone(),
                    None => continue,
                };
                let p_valid_inverted = the_dmfmvie.get(&shape_map_key(a_f));

                let mut b_remove_group = false;
                for a_f_cb in &a_lf_cb {
                    let mut b_face_passed = true;
                    for a_e_cb in explorer_shapes(a_f_cb, ShapeType::Edge) {
                        // `cxx:4019-4022`: a valid inverted edge keeps the split.
                        if let Some(a_mvie) = p_valid_inverted {
                            if a_mvie.contains(&a_e_cb) {
                                b_face_passed = false;
                                break;
                            }
                        }
                        // `cxx:4023-4029`: an edge shared with another split has
                        // to be inverted itself.
                        let a_neighbours = a_dmef
                            .get(&shape_map_key(&a_e_cb))
                            .map(|l| l.len())
                            .unwrap_or(0);
                        if a_neighbours > 1 && !self.my_inverted_edges.contains(&a_e_cb) {
                            b_face_passed = false;
                            break;
                        }
                    }
                    // `cxx:4031-4034`: the first split that passes marks the
                    // whole group.
                    if b_face_passed {
                        b_remove_group = true;
                        break;
                    }
                }
                // `cxx:4035-4043`.
                if b_remove_group {
                    for a_f_cb in &a_lf_cb {
                        a_mf_to_rem.add(a_f_cb);
                    }
                }
            }
        }

        // `cxx:4048-4069`: drop the removed splits from the faces' images.
        if a_mf_to_rem.extent() == 0 {
            return;
        }
        let a_nb = self.my_invalid_faces.extent();
        for i in 1..=a_nb {
            let Some(a_f) = self.my_invalid_faces.keys().get(i - 1).cloned() else {
                continue;
            };
            let Some(a_lf_images) = self.my_of_images.change_seek(&a_f) else {
                continue;
            };
            let mut j = 0usize;
            while j < a_lf_images.len() {
                if a_mf_to_rem.contains(&a_lf_images[j]) {
                    a_lf_images.remove(j);
                } else {
                    j += 1;
                }
            }
        }
    }
}

/// `TopExp::MapShapes(theS, theType, theMap)` accumulated into `map`
/// (`TopExp.cxx:27-46`).
fn map_shapes_into(the_s: &TopoShape, the_type: ShapeType, map: &mut ShapeIndexedMap) {
    for a_s in explorer_shapes(the_s, the_type) {
        map.add(&a_s);
    }
}
