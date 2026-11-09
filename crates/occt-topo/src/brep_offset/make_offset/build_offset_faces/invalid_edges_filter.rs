//! Edge-level filtering of `BRepOffset_BuildOffsetFaces`:
//!
//! * `MakeInvertedEdgesInvalid` (`cxx:2778-2853`) - an inverted edge locked
//!   inside a closed loop of invalid edges is invalid as well;
//! * `CheckEdgesCreatedByVertex` (`cxx:5555-5601`) - an unchecked edge whose
//!   origins are all vertices is invalid;
//! * `FilterInvalidEdges` (`cxx:5606-5750`) - drops the invalid edges that
//!   belong only to removed blocks and collects the ones a rebuild may use;
//! * `FilterEdgesImages` (`cxx:5320-5363`) - keeps only the edge images that
//!   survived the rebuild and rebuilds `myOEOrigins` from them.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::shape::TopoShape;
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_same, map_shapes_and_ancestors, shape_map_key, ShapeMapKey,
};

use super::statics::append_to_list;
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::MakeInvertedEdgesInvalid(const List&
    /// theLFOffset)` (`cxx:2778-2853`).
    ///
    /// The loops are the connexity blocks of the invalid edges connected
    /// through vertices (`cxx:2793`); only closed loops (every vertex used by
    /// exactly two edges, `cxx:2809-2818`) take part.
    pub fn make_inverted_edges_invalid(&mut self, the_lf_offset: &[TopoShape]) {
        // `cxx:2779-2782`.
        if self.my_invalid_edges.extent() == 0 || self.my_inverted_edges.extent() == 0 {
            return;
        }

        // `cxx:2785-2794`: the compound of all invalid edges and its loops.
        let a_cb_e_inv = make_compound_from(self.my_invalid_edges.iter());
        let a_lcb = crate::algo_tools::AlgoTools::make_connexity_blocks_in(
            &a_cb_e_inv,
            ShapeType::Vertex,
            ShapeType::Edge,
        );

        // `cxx:2797`: vertex -> loop of the closed loops.
        let mut a_dm_vcb: HashMap<ShapeMapKey, TopoShape> = HashMap::new();

        // `cxx:2799-2827`.
        for a_cb in &a_lcb {
            let a_dmve = map_shapes_and_ancestors(a_cb, ShapeType::Vertex, ShapeType::Edge);
            let mut is_closed = true;
            for a_lf in a_dmve.values() {
                if a_lf.len() != 2 {
                    is_closed = false;
                    break;
                }
            }
            if !is_closed {
                continue;
            }
            // `cxx:2822-2825`.
            for a_v_key in a_dmve.keys() {
                a_dm_vcb.insert(a_v_key.clone(), a_cb.clone());
            }
        }

        // `cxx:2830-2849`.
        for a_f_offset in the_lf_offset {
            let Some(a_lf_im) = self.my_of_images.seek(a_f_offset) else {
                continue;
            };
            let a_lf_im: Vec<TopoShape> = a_lf_im.to_vec();
            for a_f_im in &a_lf_im {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    if self.my_invalid_edges.contains(&a_e)
                        || !self.my_inverted_edges.contains(&a_e)
                    {
                        continue;
                    }
                    // `cxx:2841-2844`: both endpoints must sit in the same
                    // closed loop.
                    let Some(a_edge) = crate::shape::Edge::wrap(a_e.clone()) else {
                        continue;
                    };
                    let (a_v1, a_v2) = edge_vertices(&a_edge);
                    let (Some(a_v1), Some(a_v2)) = (a_v1, a_v2) else {
                        continue;
                    };
                    let p_cb1 = a_dm_vcb.get(&shape_map_key(&a_v1.0));
                    let p_cb2 = a_dm_vcb.get(&shape_map_key(&a_v2.0));
                    if let (Some(p_cb1), Some(p_cb2)) = (p_cb1, p_cb2) {
                        if is_same(p_cb1, p_cb2) {
                            self.my_invalid_edges.add(&a_e);
                        }
                    }
                }
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::CheckEdgesCreatedByVertex()`
    /// (`cxx:5555-5601`).
    pub fn check_edges_created_by_vertex(&mut self) {
        let a_nb_f = self.my_invalid_faces.extent();
        for i in 1..=a_nb_f {
            let Some(a_f) = self.my_invalid_faces.keys().get(i - 1).cloned() else {
                continue;
            };
            // `cxx:5562-5565`.
            if self.my_art_invalid_faces.contains_key(&shape_map_key(&a_f)) {
                continue;
            }
            let a_lf_im: Vec<TopoShape> = self
                .my_invalid_faces
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for a_f_im in &a_lf_im {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    if self.my_invalid_edges.contains(&a_e) || self.my_valid_edges.contains(&a_e) {
                        continue;
                    }
                    // `cxx:5580-5595`: the edge is created from a vertex when
                    // *all* of its origins are vertices.
                    let Some(a_le_or) = self.my_edges_origins.get(&shape_map_key(&a_e)) else {
                        continue;
                    };
                    if a_le_or
                        .iter()
                        .all(|o| o.shape_type() == ShapeType::Vertex)
                    {
                        self.my_invalid_edges.add(&a_e);
                    }
                }
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::FilterInvalidEdges(const
    /// BRepOffset_DataMapOfShapeIndexedMap& theDMFMIE, const IndexedMap&
    /// theMERemoved, const IndexedMap& theMEInside, Map& theMEUseInRebuild)`
    /// (`cxx:5606-5750`).
    pub fn filter_invalid_edges(
        &mut self,
        the_dmfmie: &HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_me_removed: &ShapeIndexedMap,
        the_me_inside: &ShapeIndexedMap,
        the_me_use_in_rebuild: &mut ShapeIndexedMap,
    ) {
        // `cxx:5612-5638`: the invalid edges of the invalid faces.
        let mut a_me_inv = ShapeIndexedMap::new();
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
            for a_f_im in &a_lf_inv {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    a_me_inv.add(&a_e);
                }
            }
        }
        let a_ce_inv = make_compound_from(self.my_invalid_edges.iter());

        // `cxx:5640-5662`: a loop whose edges are all removed is avoided.
        let mut a_me_inv_to_avoid = ShapeIndexedMap::new();
        let a_lcbe = crate::algo_tools::AlgoTools::make_connexity_blocks_in(
            &a_ce_inv,
            ShapeType::Vertex,
            ShapeType::Edge,
        );
        for a_cbe in &a_lcbe {
            let a_this_edges = explorer_shapes(a_cbe, ShapeType::Edge);
            let mut b_all_removed = true;
            for a_e in &a_this_edges {
                if !the_me_removed.contains(a_e) {
                    b_all_removed = false;
                    break;
                }
            }
            if b_all_removed {
                for a_e in &a_this_edges {
                    a_me_inv_to_avoid.add(a_e);
                }
            }
        }

        // `cxx:5664-5707`: the invalid edges that are still meaningful.
        let mut a_really_inv_edges = ShapeIndexedMap::new();
        let a_nb = self.my_invalid_faces.extent();
        for i in 1..=a_nb {
            let Some(a_f) = self.my_invalid_faces.keys().get(i - 1).cloned() else {
                continue;
            };
            if self.my_art_invalid_faces.contains_key(&shape_map_key(&a_f)) {
                // `cxx:5674-5686`.
                let Some(a_mie) = the_dmfmie.get(&shape_map_key(&a_f)) else {
                    continue;
                };
                for i_e in 1..=a_mie.extent() {
                    let Some(a_e) = a_mie.value(i_e) else {
                        continue;
                    };
                    if a_me_inv.contains(a_e) && !a_me_inv_to_avoid.contains(a_e) {
                        a_really_inv_edges.add(a_e);
                    }
                }
            } else {
                // `cxx:5688-5703`.
                let a_lf_inv: Vec<TopoShape> = self
                    .my_invalid_faces
                    .seek(&a_f)
                    .map(|l| l.to_vec())
                    .unwrap_or_default();
                for a_f_im in &a_lf_inv {
                    for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                        if self.my_invalid_edges.contains(&a_e)
                            && !a_me_inv_to_avoid.contains(&a_e)
                        {
                            a_really_inv_edges.add(&a_e);
                        }
                    }
                }
            }
        }
        self.my_invalid_edges = a_really_inv_edges;

        // `cxx:5709-5749`: an inside invalid edge shared with a valid one is
        // usable for a rebuild when none of its images is an outside invalid
        // edge of the same origin.
        let a_nb = self.my_invalid_edges.extent();
        for i in 1..=a_nb {
            let Some(a_e) = self.my_invalid_edges.value(i).cloned() else {
                continue;
            };
            if !the_me_inside.contains(&a_e) || !self.my_valid_edges.contains(&a_e) {
                continue;
            }
            let Some(a_e_origins) = self.my_oe_origins.get(&shape_map_key(&a_e)) else {
                continue;
            };
            let a_e_origins = a_e_origins.clone();
            let mut b_has_inv_outside = false;
            for a_e_or in &a_e_origins {
                let Some(a_e_ims) = self.my_oe_images.seek(a_e_or) else {
                    continue;
                };
                let a_e_ims: Vec<TopoShape> = a_e_ims.to_vec();
                for a_e_im in &a_e_ims {
                    b_has_inv_outside = self.my_invalid_edges.contains(a_e_im)
                        && !the_me_inside.contains(a_e_im);
                    if b_has_inv_outside {
                        break;
                    }
                }
                if b_has_inv_outside {
                    break;
                }
            }
            if !b_has_inv_outside {
                the_me_use_in_rebuild.add(&a_e);
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::FilterEdgesImages(const
    /// TopoDS_Shape& theS)` (`cxx:5320-5363`).
    ///
    /// The edges with no image left stay in `myOEImages` as keys with an empty
    /// list (`cxx:5341-5344`), which is what keeps them from being used when
    /// the splits are rebuilt.
    pub fn filter_edges_images(&mut self, the_s: &TopoShape) {
        // `cxx:5322-5323`: the surviving edges.
        let mut a_me = ShapeIndexedMap::new();
        for a_e in explorer_shapes(the_s, ShapeType::Edge) {
            a_me.add(&a_e);
        }

        // `cxx:5325`.
        self.my_oe_origins.clear();

        // `cxx:5326-5362`.
        let a_keys: Vec<TopoShape> = self.my_oe_images.keys().to_vec();
        for a_e in a_keys {
            let Some(a_le_im) = self.my_oe_images.change_seek(&a_e) else {
                continue;
            };
            // `cxx:5334-5346`: drop the images that did not survive, keeping
            // the order of the rest.
            let mut a_kept: Vec<TopoShape> = Vec::new();
            for a_e_im in a_le_im.iter() {
                if !a_me.contains(a_e_im) {
                    continue;
                }
                a_kept.push(a_e_im.clone());
            }
            *a_le_im = a_kept;

            // `cxx:5348-5358`: the origins of every surviving image.
            for a_e_im in a_le_im.iter() {
                let a_l_or = self
                    .my_oe_origins
                    .entry(shape_map_key(a_e_im))
                    .or_default();
                append_to_list(a_l_or, &a_e);
            }
        }
    }
}

/// `TopoDS_Compound` plus `BRep_Builder::Add` over `shapes`
/// (`cxx:2785-2789`, `cxx:5616-5638`).
pub(super) fn make_compound_from<'a>(shapes: impl Iterator<Item = &'a TopoShape>) -> TopoShape {
    let a_bb = crate::builder::TopoBuilder::new();
    let mut a_compound = crate::shape::Compound::new();
    for a_s in shapes {
        a_bb.add(&mut a_compound.0, a_s);
    }
    a_compound.0
}
