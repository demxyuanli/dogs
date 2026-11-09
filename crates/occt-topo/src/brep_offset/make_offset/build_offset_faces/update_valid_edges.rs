//! `UpdateValidEdges` (`cxx:7679-8236`): the rebuild epilogue that turns the
//! intersection results into the class's final edge maps.
//!
//! Stages, in the reference's order:
//!
//! 1. `cxx:7703-7733`: build `aLE` (the edges to update) and `aMELF` (edge ->
//!    faces) from `theFLE`. An edge under extended check - or already skipped
//!    once - without images is dropped from `theMECheckExt` and remembered in
//!    `aMETmp` (`cxx:7715-7720`).
//! 2. `cxx:7748-7757`: `TrimNewIntersectionEdges` trims the splits of `aLE` by
//!    the bounding vertices, filling `aMEB` (bounding splits), `aMVOld` (their
//!    vertices), `aMENew` (new edges) and `aDMEOr` (split -> its untrimmed
//!    edge).
//! 3. `cxx:7759-7766`: with no images left, only `UpdateNewIntersectionEdges`
//!    runs.
//! 4. `cxx:7771-7881`: `theOENEdges` becomes a compound split into connexity
//!    blocks (`VERTEX`/`EDGE`); per block the new edges are fenced from the
//!    block's edge images, their splits are listed, and more than one split is
//!    intersected among themselves by `IntersectEdges` (`cxx:7859-7878`).
//! 5. `cxx:7887-7966`: the first filter stage - per block, the bounds come from
//!    the splits of the block's faces only (`GetBounds`), the invalid splits
//!    from `GetInvalidEdgesByBounds`, and the surviving edges become the block's
//!    valid compound. `FilterSplits(aLE, aMEVal, false, ...)` then drops the
//!    invalid images.
//! 6. `cxx:7981-8057`: the second filter stage - the previously removed new
//!    edges are added back as markers (`cxx:7986-7997`), the surviving blocks
//!    are intersected together when there is more than one (`cxx:7999-8012`),
//!    the bounds now come from *all* `myOFImages` keys (`cxx:8025-8031`), and
//!    `FilterSplits(aLE, aMEInv, true, ...)` keeps only the valid images.
//! 7. `cxx:8061-8117`: the faces whose splits the valid splits touch. A face of
//!    `myOFImages` is updated when it is a key of `theFLE`, or when one of its
//!    splits carries a vertex of `aSplits`. `GetBoundsToUpdate` then answers its
//!    bounding edges (`aLABounds`), the offset edges to re-image (`aLAValid`) and
//!    the bounds compound (`aBounds`).
//! 8. `cxx:8110-8129`: `BOPAlgo_Builder` between the bounds and the splits, then
//!    `UpdateImages` / `UpdateNewIntersectionEdges` / `UpdateOrigins` /
//!    `UpdateIntersectedEdges` refresh the maps.
//! 9. `cxx:8130-8235`: rebuild `myEdgesToAvoid`. An avoided edge is kept when
//!    the builder left it untouched, or when its single image *is* the edge and
//!    its pave block has no non-trivial common block (`cxx:8161-8186`).
//!    Otherwise its images are cut against `myInsideEdges` (`cxx:8208-8228`),
//!    and every surviving piece is kept.
//!
//! The two `HasErrors` arms of stages 3-8 return without touching the maps, as
//! the reference's early `return`s do.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::bop_builder2::{builder_bop_algo, BoolOp2, BopBuilder};
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::Dmvv;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::invalid_edges_filter::make_compound_from;
use super::statics::{append_to_list, update_images, update_images_dmvv, update_origins};
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void UpdateValidEdges(const IndexedDataMap& theFLE, const
    /// IndexedDataMap& theOENEdges, const Map& theMVBounds, const Map&
    /// theMEInvOnArt, Map& theMECheckExt, Map& theVertsToAvoid, DataMap&
    /// theEImages, DataMap& theEETrim, const Message_ProgressRange&
    /// theRange)` (`cxx:7556-7564`, `:7679-8236`).
    #[allow(clippy::too_many_arguments)]
    pub fn update_valid_edges(
        &mut self,
        the_fle: Dmvv,
        the_oen_edges: Dmvv,
        the_mv_bounds: ShapeIndexedMap,
        the_me_inv_on_art: ShapeIndexedMap,
        mut the_me_check_ext: ShapeIndexedMap,
        the_verts_to_avoid: &mut ShapeIndexedMap,
        mut the_e_images: HashMap<ShapeMapKey, Vec<TopoShape>>,
        mut the_eetrim: HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) {
        // `cxx:7703-7706`.
        let mut a_le: Vec<TopoShape> = Vec::new();
        let mut a_melf: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut a_me_tmp = ShapeIndexedMap::new();

        // `cxx:7708-7732`.
        for (i, a_f) in the_fle.keys().iter().enumerate() {
            let Some(a_le_int) = the_fle.seek(a_f) else {
                continue;
            };
            let _ = i;
            for a_e in a_le_int {
                // `cxx:7715-7720`: an image-less edge under check is dropped
                // from the check set and remembered.
                let b_checked = the_me_check_ext.contains(a_e) || a_me_tmp.contains(a_e);
                if b_checked && !the_e_images.contains_key(&shape_map_key(a_e)) {
                    the_me_check_ext.remove(a_e);
                    a_me_tmp.add(a_e);
                    continue;
                }
                // `cxx:7720-7729`: first use of the edge binds the key.
                let a_key = shape_map_key(a_e);
                if !a_melf.contains_key(&a_key) {
                    a_melf.insert(a_key.clone(), Vec::new());
                    a_le.push(a_e.clone());
                }
                if let Some(a_lf) = a_melf.get_mut(&a_key) {
                    a_lf.push(a_f.clone());
                }
            }
        }

        // `cxx:7731-7734`.
        if a_le.is_empty() {
            return;
        }

        // `cxx:7737-7746`.
        let mut a_meb = ShapeIndexedMap::new();
        let mut a_menew = ShapeIndexedMap::new();
        let mut a_mv_old = ShapeIndexedMap::new();
        let mut a_dmeor: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        // `cxx:7748-7757`.
        self.trim_new_intersection_edges(
            &a_le,
            &the_eetrim,
            &the_mv_bounds,
            &mut the_me_check_ext,
            &mut the_e_images,
            &mut a_meb,
            &mut a_mv_old,
            &mut a_menew,
            &mut a_dmeor,
            &mut a_melf,
        );

        // `cxx:7759-7766`.
        if the_e_images.is_empty() {
            self.update_new_intersection_edges(&a_le, &a_melf, &the_e_images, &mut the_eetrim);
            return;
        }

        // `cxx:7777-7787`: the compound of the invalidity edges.
        let a_ce_all = make_compound_from(the_oen_edges.keys().iter());
        // `cxx:7789-7790`.
        let a_lblocks = AlgoTools::make_connexity_blocks_in(
            &a_ce_all,
            ShapeType::Vertex,
            ShapeType::Edge,
        );

        // `cxx:7797-7880`: `aMBlocksSp`, an `IndexedDataMap` (key order kept).
        let mut a_mblocks_sp: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
        for a_block in &a_lblocks {
            // `cxx:7809-7826`: the new edges of the block, fenced.
            let mut a_block_lenew: Vec<TopoShape> = Vec::new();
            let mut a_me_fence = ShapeIndexedMap::new();
            for a_e in explorer_shapes(a_block, ShapeType::Edge) {
                let Some(a_le_int) = the_oen_edges.seek(&a_e) else {
                    continue;
                };
                for a_e_int in a_le_int {
                    if a_me_fence.add(a_e_int) == a_me_fence.extent() {
                        a_block_lenew.push(a_e_int.clone());
                    }
                }
            }

            // `cxx:7828-7831`.
            if a_block_lenew.is_empty() {
                continue;
            }

            // `cxx:7835-7857`: the splits of those edges.
            let mut a_lsplits: Vec<TopoShape> = Vec::new();
            for a_e in &a_block_lenew {
                let Some(a_le_im) = the_e_images.get(&shape_map_key(a_e)) else {
                    continue;
                };
                for a_e_im in a_le_im {
                    a_lsplits.push(a_e_im.clone());
                }
            }

            // `cxx:7851-7854`.
            if a_lsplits.is_empty() {
                continue;
            }

            // `cxx:7858-7877`.
            let a_ce = if a_lsplits.len() > 1 {
                self.intersect_edges(
                    &a_lsplits,
                    &a_block_lenew,
                    &the_mv_bounds,
                    the_verts_to_avoid,
                    &mut a_menew,
                    &mut the_me_check_ext,
                    &mut the_e_images,
                    &mut a_dmeor,
                    &mut a_melf,
                )
            } else {
                a_lsplits[0].clone()
            };

            // `cxx:7879`: `aMBlocksSp.Add(aCE, aBlockLENew)` - an existing key
            // keeps its position and takes the new list.
            match a_mblocks_sp
                .iter_mut()
                .find(|(a_key, _)| is_same(a_key, &a_ce))
            {
                Some(a_slot) => a_slot.1 = a_block_lenew,
                None => a_mblocks_sp.push((a_ce, a_block_lenew)),
            }
        }

        // `cxx:7889-7966`: the first filter stage, per block.
        let mut a_me_val = ShapeIndexedMap::new();
        let mut a_lval_blocks: Vec<TopoShape> = Vec::new();
        for (a_ce, a_block_lenew) in &a_mblocks_sp {
            // `cxx:7904-7921`: the faces behind the block's edges.
            let mut a_lfaces: Vec<TopoShape> = Vec::new();
            for a_e in a_block_lenew {
                let Some(a_lf) = a_melf.get(&shape_map_key(a_e)) else {
                    continue;
                };
                for a_f in a_lf {
                    append_to_list(&mut a_lfaces, a_f);
                }
            }
            // `cxx:7925`.
            let a_filter_bounds = self.get_bounds(&a_lfaces, &a_meb);
            // `cxx:7928-7939`.
            let mut a_me_inv_loc = ShapeIndexedMap::new();
            self.get_invalid_edges_by_bounds(
                a_ce,
                &a_filter_bounds,
                &a_mv_old,
                &a_menew,
                &a_dmeor,
                &a_melf,
                &the_e_images,
                &the_me_check_ext,
                &the_me_inv_on_art,
                the_verts_to_avoid,
                &mut a_me_inv_loc,
            );

            // `cxx:7943-7965`: keep the block's valid edges.
            let mut a_ce_val_edges: Vec<TopoShape> = Vec::new();
            for a_e_sp in explorer_shapes(a_ce, ShapeType::Edge) {
                if !a_me_inv_loc.contains(&a_e_sp) && a_me_val.add(&a_e_sp) == a_me_val.extent() {
                    a_ce_val_edges.push(a_e_sp);
                }
            }
            if !a_ce_val_edges.is_empty() {
                a_lval_blocks.push(make_compound_from(a_ce_val_edges.iter()));
            }
        }

        // `cxx:7968-7969`: filter the images after the first stage. The
        // answered compound is not read - stage 2 overwrites `aSplits1` in both
        // of its arms (`cxx:8006`, `:8013`) - but the filtering of
        // `theEImages` is the point of the call.
        let _a_splits1 = self.filter_splits(&a_le, &a_me_val, false, &mut the_e_images);

        // `cxx:7968-7974`: nothing survived - only the intersection edges are
        // updated.
        if a_lval_blocks.is_empty() {
            self.update_new_intersection_edges(&a_le, &a_melf, &the_e_images, &mut the_eetrim);
            return;
        }

        // `cxx:7981-7998`: the second stage adds the removed new edges as
        // markers.
        for (a_ce, _) in &a_mblocks_sp {
            for a_e_im in explorer_shapes(a_ce, ShapeType::Edge) {
                if a_menew.contains(&a_e_im) && !a_me_val.contains(&a_e_im) {
                    a_lval_blocks.push(a_e_im);
                }
            }
        }

        // `cxx:7999-8014`.
        let a_splits1 = if a_lval_blocks.len() > 1 {
            self.intersect_edges(
                &a_lval_blocks,
                &a_le,
                &the_mv_bounds,
                the_verts_to_avoid,
                &mut a_menew,
                &mut the_me_check_ext,
                &mut the_e_images,
                &mut a_dmeor,
                &mut a_melf,
            )
        } else {
            a_lval_blocks[0].clone()
        };

        // `cxx:8025-8031`: the bounds of every offset face's splits.
        let a_lfaces_all: Vec<TopoShape> = self.my_of_images.keys().to_vec();
        let a_filter_bounds = self.get_bounds(&a_lfaces_all, &a_meb);

        // `cxx:8034-8048`.
        let mut a_me_inv = ShapeIndexedMap::new();
        self.get_invalid_edges_by_bounds(
            &a_splits1,
            &a_filter_bounds,
            &a_mv_old,
            &a_menew,
            &a_dmeor,
            &a_melf,
            &the_e_images,
            &the_me_check_ext,
            &the_me_inv_on_art,
            the_verts_to_avoid,
            &mut a_me_inv,
        );

        // `cxx:8050-8054`: keep only the valid images.
        let a_splits = self.filter_splits(&a_le, &a_me_inv, true, &mut the_e_images);

        // `cxx:8061-8108`: the faces whose splits the valid splits touch.
        let mut a_mv_sp = ShapeIndexedMap::new();
        for a_v in explorer_shapes(&a_splits, ShapeType::Vertex) {
            a_mv_sp.add(&a_v);
        }
        let mut a_lf: Vec<TopoShape> = Vec::new();
        for a_f in self.my_of_images.keys().to_vec() {
            // `cxx:8071-8075`: a face already in `theFLE` is always updated.
            if the_fle.contains(&a_f) {
                a_lf.push(a_f);
                continue;
            }
            // `cxx:8078-8101`: otherwise some split must carry a new vertex.
            let a_lf_im = self.my_of_images.seek(&a_f).unwrap_or(&[]).to_vec();
            let mut b_found = false;
            for a_f_im in &a_lf_im {
                if explorer_shapes(a_f_im, ShapeType::Vertex)
                    .iter()
                    .any(|a_v| a_mv_sp.contains(a_v))
                {
                    b_found = true;
                    break;
                }
            }
            if b_found {
                a_lf.push(a_f);
            }
        }

        // `cxx:8109-8117`.
        let (a_labounds, a_lavalid, a_bounds) = self.get_bounds_to_update(&a_lf, &a_meb);

        // `cxx:8118-8127`: the general fuse of the bounds and the valid splits.
        let mut a_gf = BopBuilder::new();
        a_gf.add_argument(&a_bounds);
        a_gf.add_argument(&a_splits);
        a_gf.perform().ok();

        // `cxx:8117-8129`.
        update_images(&a_le, &mut the_e_images, &a_gf, &mut self.my_modified_edges);
        self.update_new_intersection_edges(&a_le, &a_melf, &the_e_images, &mut the_eetrim);
        update_images_dmvv(&a_lavalid, &mut self.my_oe_images, &a_gf, &mut self.my_modified_edges);
        update_origins(&a_labounds, &mut self.my_oe_origins, &a_gf);
        update_origins(&a_labounds, &mut self.my_edges_origins, &a_gf);
        self.update_intersected_edges(&a_labounds, &a_gf);

        // `cxx:8131-8139`: the edges the fuse created from `aSplits`.
        let mut a_new_edges = ShapeIndexedMap::new();
        let a_splits_im: Vec<TopoShape> = a_gf
            .history()
            .images()
            .get(&shape_map_key(&a_splits))
            .cloned()
            .unwrap_or_default();
        for a_s in &a_splits_im {
            for a_e in explorer_shapes(a_s, ShapeType::Edge) {
                a_new_edges.add(&a_e);
            }
        }

        // `cxx:8141-8146`: the compound of the inside edges.
        let an_inside_edges =
            make_compound_from(self.my_inside_edges.iter());

        // `cxx:8149-8202`: rebuild the avoid map.
        let mut a_me_avoid = ShapeIndexedMap::new();
        let mut a_ce_avoid_children: Vec<TopoShape> = Vec::new();
        let a_nb_e = self.my_edges_to_avoid.extent();
        for i in 1..=a_nb_e {
            let Some(a_e) = self.my_edges_to_avoid.value(i).cloned() else {
                continue;
            };
            // `cxx:8160`.
            let a_le_im: Vec<TopoShape> = a_gf
                .history()
                .modified(&a_e)
                .map(|l| l.to_vec())
                .unwrap_or_default();

            // `cxx:8163-8186`: only untouched and fully coinciding edges are
            // kept in the avoid map.
            let mut b_keep = a_le_im.is_empty();
            if a_le_im.len() == 1 && is_same(&a_e, &a_le_im[0]) {
                // `cxx:8167-8185`.
                if let Some(an_index) = a_gf.ds().index(&a_e) {
                    let a_lpb = a_gf.ds().pave_blocks(an_index);
                    if a_lpb.len() == 1 {
                        if let Some(a_cb) = a_gf.ds().common_block(&a_lpb[0]) {
                            // A pave block of the common block whose original
                            // edge has more than one pave block breaks the
                            // keep.
                            let mut b_break = false;
                            for a_pb_cb in a_cb.pave_blocks() {
                                if a_gf
                                    .ds()
                                    .pave_blocks(a_pb_cb.original_edge())
                                    .len()
                                    > 1
                                {
                                    b_break = true;
                                    break;
                                }
                            }
                            b_keep = !b_break;
                        }
                    }
                }
            }

            // `cxx:8188-8192`.
            if b_keep {
                a_me_avoid.add(&a_e);
                continue;
            }

            // `cxx:8194-8201`.
            for a_e_im in &a_le_im {
                if !a_new_edges.contains(a_e_im) {
                    if !a_ce_avoid_children.iter().any(|a| is_same(a, a_e_im)) {
                        a_ce_avoid_children.push(a_e_im.clone());
                    }
                }
            }
        }

        // `cxx:8204-8228`: cut the avoid edges against the inside edges.
        let mut b_is_cut = false;
        if !a_ce_avoid_children.is_empty() {
            let a_ce_avoid = make_compound_from(a_ce_avoid_children.iter());
            let mut a_bop = builder_bop_algo(&[a_ce_avoid], &[an_inside_edges], BoolOp2::Cut);
            b_is_cut = !a_bop.has_errors();
            let _ = &mut a_bop;
            if b_is_cut {
                for a_e in &a_ce_avoid_children {
                    if !a_bop.history().is_deleted(a_e) {
                        a_me_avoid.add(a_e);
                    }
                }
            }
        }

        // `cxx:8229-8235`.
        if !b_is_cut {
            let a_ce_avoid = make_compound_from(a_ce_avoid_children.iter());
            for a_e in explorer_shapes(&a_ce_avoid, ShapeType::Edge) {
                a_me_avoid.add(&a_e);
            }
        }

        // `cxx:8236`.
        self.my_edges_to_avoid = a_me_avoid;
    }
}

