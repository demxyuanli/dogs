//! `IntersectAndTrimEdges` (`cxx:7319-7604`): the intersection edges of one
//! invalidity are intersected among themselves and, afterwards, the common
//! intersection edges are trimmed by the newly created ones.
//!
//! The first half (`cxx:7338-7450`) gathers the arguments of the edge
//! intersection: the vertices of the splits of the faces to rebuild that are
//! not to be avoided and that belong to the invalidity (`cxx:7351-7372`), the
//! vertices carried by the artificial-face intersection records
//! (`cxx:7379-7392`), the vertices reachable through the invalidity shapes
//! (`cxx:7394-7412`), and every edge to intersect together with its bounding
//! vertices (`cxx:7424-7450`). `BOPAlgo_Builder` then intersects them
//! (`cxx:7452-7458`).
//!
//! The bounding vertices of the result (`cxx:7461-7477`) feed `theMVBounds`,
//! `GetInvalidEdges` marks the splits meeting an invalid vertex
//! (`cxx:7479-7480`), and the surviving splits of the edges to intersect are
//! collected in the `aCEIm` compound and in `theEImages` (`cxx:7486-7526`).
//!
//! The second half (`cxx:7527-7601`) intersects the common intersection edges
//! with `aCEIm`, skips the ones that coincide with a common block whose other
//! original edges are not in `theMECheckExt` (`cxx:7546-7582`), and otherwise
//! saves their images and the bounding vertices of the trimmed parts
//! (`cxx:7585-7599`).

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::bop_builder2::BopBuilder;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::shape::TopoShape;
use crate::topo_tools::{
    explorer_shapes, map_shapes_and_ancestors_into, shape_map_key, ShapeMapKey,
};

use super::statics::append_to_list;
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::IntersectAndTrimEdges(const
    /// IndexedMap& theMFInt, const IndexedMap& theMEInt, const DataMap&
    /// theDMEETrim, const IndexedMap& theMSInv, const IndexedMap& theMVE,
    /// const Map& theVertsToAvoid, const Map& theNewVertsToAvoid, const Map&
    /// theMECheckExt, const DataMap* theSSInterfs, Map& theMVBounds, DataMap&
    /// theEImages)` (`cxx:7319-7604`).
    #[allow(clippy::too_many_arguments)]
    pub fn intersect_and_trim_edges(
        &self,
        the_mf_int: &ShapeIndexedMap,
        the_me_int: &ShapeIndexedMap,
        the_dmeet_trim: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_ms_inv: &ShapeIndexedMap,
        the_mve: &ShapeIndexedMap,
        the_verts_to_avoid: &ShapeIndexedMap,
        the_new_verts_to_avoid: &ShapeIndexedMap,
        the_me_check_ext: &ShapeIndexedMap,
        the_ss_interfs: Option<&HashMap<ShapeMapKey, Vec<TopoShape>>>,
        the_mv_bounds: &mut ShapeIndexedMap,
        the_e_images: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) {
        // `cxx:7337-7340`.
        if the_me_int.extent() == 0 {
            return;
        }

        // `cxx:7341-7345`.
        let mut a_l_args: Vec<TopoShape> = Vec::new();
        let mut a_mfence = ShapeIndexedMap::new();
        let mut a_dmve: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        // `cxx:7351-7372`: the vertices of the splits close to the invalidity.
        let a_nb = the_mf_int.extent();
        for i in 1..=a_nb {
            let Some(a_f) = the_mf_int.value(i).cloned() else {
                continue;
            };
            // `myFacesToRebuild.FindFromKey(aF)` raises on an unbound face in
            // the reference; the port reads an empty list.
            let a_le: Vec<TopoShape> = self
                .my_faces_to_rebuild
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for a_e in &a_le {
                map_shapes_and_ancestors_into(a_e, ShapeType::Vertex, ShapeType::Edge, &mut a_dmve);
                for a_v1 in explorer_shapes(a_e, ShapeType::Vertex) {
                    if !the_verts_to_avoid.contains(&a_v1)
                        && the_mve.contains(&a_v1)
                        && a_mfence.add(&a_v1) == a_mfence.extent()
                    {
                        append_to_list(&mut a_l_args, &a_v1);
                    }
                }
            }
        }

        // `cxx:7374-7413`: the vertices of the invalidity shapes.
        let a_nb = the_ms_inv.extent();
        for i in 1..=a_nb {
            let Some(a_s) = the_ms_inv.value(i).cloned() else {
                continue;
            };
            // `cxx:7379-7392`: the vertices of the artificial-face section
            // records. They are appended without the `aMFence` check.
            if let Some(a_map) = the_ss_interfs {
                if let Some(a_lv) = a_map.get(&shape_map_key(&a_s)) {
                    for a_value in a_lv {
                        if a_value.shape_type() == ShapeType::Vertex {
                            append_to_list(&mut a_l_args, a_value);
                        }
                    }
                }
            }
            // `cxx:7394-7412`: the vertex case.
            if let Some(a_lve) = a_dmve.get(&shape_map_key(&a_s)).cloned() {
                for a_e in &a_lve {
                    for a_v1 in explorer_shapes(a_e, ShapeType::Vertex) {
                        if !the_verts_to_avoid.contains(&a_v1)
                            && a_mfence.add(&a_v1) == a_mfence.extent()
                        {
                            append_to_list(&mut a_l_args, &a_v1);
                        }
                    }
                }
            }
        }

        // `cxx:7415-7450`.
        let mut a_lv_bounds: Vec<TopoShape> = Vec::new();
        let mut a_le_new: Vec<TopoShape> = Vec::new();
        let mut a_le_int: Vec<TopoShape> = Vec::new();
        let mut a_lce: Vec<TopoShape> = Vec::new();
        let a_nb = the_me_int.extent();
        for i in 1..=a_nb {
            let Some(a_e) = the_me_int.value(i).cloned() else {
                continue;
            };
            // `cxx:7429-7434`: the additional common edges are trimmed
            // separately.
            if the_me_check_ext.contains(&a_e) {
                append_to_list(&mut a_lce, &a_e);
                continue;
            }
            // `cxx:7436-7439`: the edges whose trimming failed. The reference
            // fills this list and never reads it again, so the port keeps it
            // filled for the same control flow only.
            if !the_dmeet_trim.contains_key(&shape_map_key(&a_e)) {
                append_to_list(&mut a_le_new, &a_e);
            }
            append_to_list(&mut a_le_int, &a_e);
            append_to_list(&mut a_l_args, &a_e);
            // `cxx:7445-7449`.
            for a_v in explorer_shapes(&a_e, ShapeType::Vertex) {
                append_to_list(&mut a_lv_bounds, &a_v);
            }
        }
        let _ = a_le_new;

        // `cxx:7452-7458`.
        let mut a_gf = BopBuilder::new();
        a_gf.set_arguments(&a_l_args);
        let _ = a_gf.perform();
        if a_gf.has_errors() {
            return;
        }

        // `cxx:7460-7477`: the bounding vertices of the trimmed result.
        for a_v in &a_lv_bounds {
            let a_lv_im = a_gf.history().modified(a_v).unwrap_or_default();
            match a_lv_im.first() {
                Some(a_v_im) => {
                    the_mv_bounds.add(a_v_im);
                }
                None => {
                    the_mv_bounds.add(a_v);
                }
            }
        }

        // `cxx:7479-7480`.
        let mut a_me_inv = ShapeIndexedMap::new();
        self.get_invalid_edges(the_new_verts_to_avoid, the_mv_bounds, &a_gf, &mut a_me_inv);

        // `cxx:7482-7526`: drop the splits that touch an invalid vertex and
        // collect the rest.
        let mut a_ce_im: Vec<TopoShape> = Vec::new();
        for a_e in &a_le_int {
            let a_le_im = a_gf.history().modified(a_e).unwrap_or_default();
            if a_le_im.is_empty() {
                continue;
            }
            // `cxx:7498-7512`: the invalid splits are removed from the local
            // copy only; the valid ones go into the compound.
            let mut a_valid: Vec<TopoShape> = Vec::new();
            for a_e_im in a_le_im {
                if a_me_inv.contains(a_e_im) {
                    continue;
                }
                append_to_list(&mut a_ce_im, a_e_im);
                a_valid.push(a_e_im.clone());
            }
            // `cxx:7516-7524`.
            if !a_valid.is_empty() {
                let entry = the_e_images.entry(shape_map_key(a_e)).or_default();
                entry.extend(a_valid);
            }
        }

        // `cxx:7526-7529`.
        if a_lce.is_empty() {
            return;
        }

        // `cxx:7530-7540`: trim the common edges by the other intersection
        // edges. `aCEIm` is added as an argument even when it is empty
        // (`cxx:7533`).
        let a_ce_im_compound = super::invalid_edges_filter::make_compound_from(a_ce_im.iter());
        let mut a_gfce = BopBuilder::new();
        a_gfce.set_arguments(&a_lce);
        a_gfce.add_argument(&a_ce_im_compound);
        let _ = a_gfce.perform();
        if a_gfce.has_errors() {
            return;
        }

        // `cxx:7542-7599`.
        let p_ds = a_gfce.ds();
        for a_e in &a_lce {
            let a_le_im = a_gfce.history().modified(a_e).unwrap_or_default();
            if a_le_im.is_empty() {
                continue;
            }
            // `cxx:7552-7582`: skip a split coinciding with an intersection
            // edge. The outer loop stops as soon as one common block holds an
            // original edge outside `theMECheckExt`.
            let mut b_skip = false;
            if let Some(n_e) = p_ds.index(a_e) {
                for pb in p_ds.pave_blocks(n_e) {
                    if !p_ds.is_common_block(pb) {
                        continue;
                    }
                    let mut b_outside = false;
                    if let Some(cb) = p_ds.common_block(pb) {
                        for pbc in &cb.pave_blocks {
                            let b_contains = p_ds
                                .shape(pbc.original_edge())
                                .map(|s| the_me_check_ext.contains(s))
                                .unwrap_or(false);
                            if !b_contains {
                                b_outside = true;
                                break;
                            }
                        }
                    }
                    if b_outside {
                        b_skip = true;
                        break;
                    }
                }
            }
            if b_skip {
                continue;
            }

            // `cxx:7585-7593`: save the images.
            let entry = the_e_images.entry(shape_map_key(a_e)).or_default();
            entry.extend(a_le_im.iter().cloned());

            // `cxx:7595-7599`: save the bounding vertices.
            for a_v in explorer_shapes(a_e, ShapeType::Vertex) {
                let a_lv_im = a_gfce.history().modified(&a_v).unwrap_or_default();
                match a_lv_im.first() {
                    Some(a_v_im) => {
                        the_mv_bounds.add(a_v_im);
                    }
                    None => {
                        the_mv_bounds.add(&a_v);
                    }
                }
            }
        }
    }
}
