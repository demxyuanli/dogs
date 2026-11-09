//! Removal of the invalid splits created by inverted edges -
//! `CheckInvertedBlock` (`cxx:3545-3681`), its file-local helper
//! `GetVerticesOnEdges` (`cxx:3522-3542`) and
//! `RemoveInvalidSplitsByInvertedEdges` (`cxx:3685-3886`).
//!
//! The reference refuses to remove a split whose removal would leave a hole in
//! the offset shell. Three conditions make a connexity block of splits
//! removable (`cxx:3559-3678`):
//!
//! 1. the block holds more than one split;
//! 2. it holds at least two vertex-connected inverted edges with different
//!    origins - two images of one edge do not count;
//! 3. none of the vertices of its inverted edges is shared with another block.
//!
//! `RemoveInvalidSplitsByInvertedEdges` then treats an invalid split as
//! removable when at least two of its inverted edges are not connected to each
//! other (`cxx:3821-3841`) and one of their vertices is used by more than three
//! splits (`cxx:3786-3799`), which is what makes the removal leave no hole.

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
use super::statics::append_to_list;
use super::BRepOffsetBuildOffsetFaces;

/// `static void GetVerticesOnEdges(const TopoDS_Shape& theCB, const
/// IndexedMap& theEdges, Map& theVerticesOnEdges, Map& theAllVertices)`
/// (`cxx:3522-3542`).
///
/// Every vertex of the block goes into `the_all_vertices`; the vertices of the
/// edges present in `the_edges` also go into `the_vertices_on_edges`.
fn get_vertices_on_edges(
    the_cb: &TopoShape,
    the_edges: &ShapeIndexedMap,
    the_vertices_on_edges: &mut ShapeIndexedMap,
    the_all_vertices: &mut ShapeIndexedMap,
) {
    for a_e in explorer_shapes(the_cb, ShapeType::Edge) {
        // `cxx:3532`.
        let is_on_given_edges = the_edges.contains(&a_e);
        // `cxx:3533`: `TopoDS_Iterator aItV(aE)`.
        for a_v in children_keep_orientation(&a_e) {
            the_all_vertices.add(&a_v);
            if is_on_given_edges {
                the_vertices_on_edges.add(&a_v);
            }
        }
    }
}

impl BRepOffsetBuildOffsetFaces {
    /// `bool BRepOffset_BuildOffsetFaces::CheckInvertedBlock(const
    /// TopoDS_Shape& theCB, const List& theLCBF, DataMap&
    /// theDMCBVInverted, DataMap& theDMCBVAll)` (`cxx:3545-3681`).
    ///
    /// The two `DataMap<Shape, Map<Shape>>` caches are the reference's
    /// per-block vertex sets, kept across the blocks of one face so each block
    /// is walked once. Both are carried by
    /// `HashMap<ShapeMapKey, ShapeIndexedMap>`.
    pub fn check_inverted_block(
        &mut self,
        the_cb: &TopoShape,
        the_lcb_f: &[TopoShape],
        the_dm_cb_v_inverted: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_dm_cb_v_all: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
    ) -> bool {
        // 1. `cxx:3561-3564`: more than one split in the block.
        if children_keep_orientation(the_cb).len() < 2 {
            return false;
        }

        // 2. `cxx:3568-3588`: the inverted edges of the block.
        let mut a_me_cb_inv = ShapeIndexedMap::new();
        let mut a_ce_cb_inv_edges: Vec<TopoShape> = Vec::new();
        for a_e in explorer_shapes(the_cb, ShapeType::Edge) {
            if self.my_inverted_edges.contains(&a_e) {
                if a_me_cb_inv.add(&a_e) == a_me_cb_inv.extent() {
                    a_ce_cb_inv_edges.push(a_e);
                }
            }
        }
        if a_me_cb_inv.extent() < 2 {
            return false;
        }

        // `cxx:3591-3632`: at least one vertex-connected group of inverted edges
        // with two different origins.
        let a_ce_cb_inv = make_compound_from(a_ce_cb_inv_edges.iter());
        let a_lcbe = AlgoTools::make_connexity_blocks_in(
            &a_ce_cb_inv,
            ShapeType::Vertex,
            ShapeType::Edge,
        );
        let mut b_has_pair = false;
        for a_cbe in &a_lcbe {
            let mut a_nb_unique = 0usize;
            let mut a_me_origins = ShapeIndexedMap::new();
            // `cxx:3602`: `TopoDS_Iterator aItE(aCBE)`.
            for a_e in children_keep_orientation(a_cbe) {
                match self.my_oe_origins.get(&shape_map_key(&a_e)) {
                    // `cxx:3607-3611`: no origins recorded - the edge is its
                    // own origin.
                    None => {
                        a_me_origins.add(&a_e);
                        a_nb_unique += 1;
                    }
                    Some(a_le_or) => {
                        let a_le_or = a_le_or.clone();
                        for a_e_or in &a_le_or {
                            if a_me_origins.add(a_e_or) == a_me_origins.extent() {
                                a_nb_unique += 1;
                            }
                        }
                    }
                }
            }
            // `cxx:3623-3626`.
            if a_nb_unique >= 2 {
                b_has_pair = true;
                break;
            }
        }
        if !b_has_pair {
            return false;
        }

        // 3. `cxx:3635-3649`: the vertex sets of this block.
        let a_cb_key = shape_map_key(the_cb);
        if !the_dm_cb_v_inverted.contains_key(&a_cb_key) {
            let mut a_mv_inverted = ShapeIndexedMap::new();
            let mut a_mv_all = ShapeIndexedMap::new();
            get_vertices_on_edges(the_cb, &self.my_inverted_edges, &mut a_mv_inverted, &mut a_mv_all);
            the_dm_cb_v_inverted.insert(a_cb_key.clone(), a_mv_inverted);
            the_dm_cb_v_all.insert(a_cb_key.clone(), a_mv_all);
        }

        // `cxx:3651-3676`.
        for a_cb1 in the_lcb_f {
            if is_same(a_cb1, the_cb) {
                continue;
            }
            let a_cb1_key = shape_map_key(a_cb1);
            if !the_dm_cb_v_inverted.contains_key(&a_cb1_key) {
                let mut a_mv_inverted1 = ShapeIndexedMap::new();
                let mut a_mv_all1 = ShapeIndexedMap::new();
                get_vertices_on_edges(
                    a_cb1,
                    &self.my_inverted_edges,
                    &mut a_mv_inverted1,
                    &mut a_mv_all1,
                );
                the_dm_cb_v_inverted.insert(a_cb1_key.clone(), a_mv_inverted1);
                the_dm_cb_v_all.insert(a_cb1_key.clone(), a_mv_all1);
            }

            // `cxx:3675-3678`: `NCollection_MapAlgo::HasIntersection`.
            let a_intersects = {
                let a_mv_inverted = the_dm_cb_v_inverted.get(&a_cb_key);
                let a_mv_all1 = the_dm_cb_v_all.get(&a_cb1_key);
                match (a_mv_inverted, a_mv_all1) {
                    (Some(a), Some(b)) => a.iter().any(|v| b.contains(v)),
                    _ => false,
                }
            };
            if a_intersects {
                return false;
            }
        }

        true
    }

    /// `void BRepOffset_BuildOffsetFaces::RemoveInvalidSplitsByInvertedEdges(
    /// IndexedMap& theMERemoved)` (`cxx:3685-3886`).
    pub fn remove_invalid_splits_by_inverted_edges(&mut self, the_me_removed: &mut ShapeIndexedMap) {
        // `cxx:3690-3693`.
        if self.my_inverted_edges.extent() == 0 {
            return;
        }

        // `cxx:3701-3702`: edges that must be kept, and vertex -> splits.
        let mut a_me_avoid = ShapeIndexedMap::new();
        let mut a_dm_vf: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        // `cxx:3703-3779`: per offset face, judge the regularity of its splits.
        let a_of_keys: Vec<TopoShape> = self.my_of_images.keys().to_vec();
        for a_f_offset in a_of_keys {
            let Some(a_lf_im) = self.my_of_images.seek(&a_f_offset) else {
                continue;
            };
            let a_lf_im: Vec<TopoShape> = a_lf_im.to_vec();

            // `cxx:3707-3712`.
            let a_cf_im = make_compound_from(a_lf_im.iter());
            let mut a_dmef: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

            // `cxx:3713-3750`.
            for a_f in &a_lf_im {
                // `cxx:3719-3736`: only the outer edges may be removed.
                for a_e in explorer_shapes(a_f, ShapeType::Edge) {
                    let a_key = shape_map_key(&a_e);
                    if a_dmef.contains_key(&a_key) {
                        a_me_avoid.add(&a_e);
                    }
                    let a_lf = a_dmef.entry(a_key).or_default();
                    append_to_list(a_lf, a_f);
                }
                // `cxx:3738-3749`.
                for a_v in explorer_shapes(a_f, ShapeType::Vertex) {
                    let a_lf = a_dm_vf.entry(shape_map_key(&a_v)).or_default();
                    append_to_list(a_lf, a_f);
                }
            }

            // `cxx:3753-3758`: the splits have to form a single block.
            let a_lcb_f = AlgoTools::make_connexity_blocks_in(
                &a_cf_im,
                ShapeType::Edge,
                ShapeType::Face,
            );
            if a_lcb_f.len() == 1 {
                continue;
            }

            // `cxx:3761-3779`.
            let mut a_dm_cb_v_inverted: HashMap<ShapeMapKey, ShapeIndexedMap> = HashMap::new();
            let mut a_dm_cb_v_all: HashMap<ShapeMapKey, ShapeIndexedMap> = HashMap::new();
            for a_cb in &a_lcb_f {
                if !self.check_inverted_block(
                    a_cb,
                    &a_lcb_f,
                    &mut a_dm_cb_v_inverted,
                    &mut a_dm_cb_v_all,
                ) {
                    // `cxx:3775-3777`: none of the edges of this block may go.
                    for a_e in explorer_shapes(a_cb, ShapeType::Edge) {
                        a_me_avoid.add(&a_e);
                    }
                }
            }
        }

        // `cxx:3782-3802`: an inverted edge that no block protects and whose
        // vertex is shared by more than three splits may be removed.
        let mut a_me_rem = ShapeIndexedMap::new();
        let a_nb_inverted = self.my_inverted_edges.extent();
        for i_inverted in 1..=a_nb_inverted {
            let Some(a_e) = self.my_inverted_edges.value(i_inverted).cloned() else {
                continue;
            };
            if a_me_avoid.contains(&a_e) {
                continue;
            }
            // `cxx:3788`: `TopoDS_Iterator aIt(aE)`.
            for a_v in children_keep_orientation(&a_e) {
                if let Some(a_lf) = a_dm_vf.get(&shape_map_key(&a_v)) {
                    if a_lf.len() > 3 {
                        a_me_rem.add(&a_e);
                        break;
                    }
                }
            }
        }
        if a_me_rem.extent() == 0 {
            return;
        }

        // `cxx:3807-3857`: an invalid split whose removed inverted edges form
        // two or more unconnected groups can be dropped.
        let mut a_inv_faces = Dmvv::new();
        let mut a_mf_rem = ShapeIndexedMap::new();
        let mut a_mf_to_update = ShapeIndexedMap::new();

        let a_inv_keys: Vec<TopoShape> = self.my_invalid_faces.keys().to_vec();
        for a_f in a_inv_keys {
            let Some(a_lf_im) = self.my_invalid_faces.seek(&a_f) else {
                continue;
            };
            let mut a_lf_im: Vec<TopoShape> = a_lf_im.to_vec();

            // `cxx:3819-3852`: remove the splits in place.
            let mut i = 0usize;
            while i < a_lf_im.len() {
                let a_f_im = a_lf_im[i].clone();
                // `cxx:3823-3837`: the removed inverted edges of the split.
                let mut a_ce_inv_edges: Vec<TopoShape> = Vec::new();
                for a_e in explorer_shapes(&a_f_im, ShapeType::Edge) {
                    if a_me_rem.contains(&a_e) {
                        a_ce_inv_edges.push(a_e);
                    }
                }
                let a_ce_inv = make_compound_from(a_ce_inv_edges.iter());
                let a_lcbe = AlgoTools::make_connexity_blocks_in(
                    &a_ce_inv,
                    ShapeType::Vertex,
                    ShapeType::Edge,
                );
                // `cxx:3841-3851`.
                if a_lcbe.len() >= 2 {
                    a_mf_to_update.add(&a_f);
                    a_mf_rem.add(&a_f_im);
                    a_lf_im.remove(i);
                } else {
                    i += 1;
                }
            }

            // `cxx:3854-3857`.
            if !a_lf_im.is_empty() {
                a_inv_faces.bind_list(&a_f, &a_lf_im);
            }
        }

        if a_mf_rem.extent() == 0 {
            return;
        }

        // `cxx:3865`.
        self.my_invalid_faces = a_inv_faces;

        // `cxx:3866-3886`: drop the removed splits from the faces' images.
        let a_nb_update = a_mf_to_update.extent();
        for i in 1..=a_nb_update {
            let Some(a_f) = a_mf_to_update.value(i).cloned() else {
                continue;
            };
            let Some(a_lf_im) = self.my_of_images.change_seek(&a_f) else {
                continue;
            };
            let mut j = 0usize;
            while j < a_lf_im.len() {
                let a_f_im = a_lf_im[j].clone();
                if a_mf_rem.contains(&a_f_im) {
                    for a_e in explorer_shapes(&a_f_im, ShapeType::Edge) {
                        the_me_removed.add(&a_e);
                    }
                    a_lf_im.remove(j);
                } else {
                    j += 1;
                }
            }
        }
    }
}
