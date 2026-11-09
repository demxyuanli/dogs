//! The edge-trimming chain of `UpdateValidEdges`: `TrimNewIntersectionEdges`
//! (`cxx:8239-8401`), `IntersectEdges` (`cxx:8408-8545`) and
//! `UpdateNewIntersectionEdges` (`cxx:9048-9292`).
//!
//! `TrimNewIntersectionEdges` walks the new edges and, when an edge has more
//! than one split, fuses the splits together with the two bounding vertices of
//! the untrimmed edge (`cxx:8297-8338`). Every split containing a bounding
//! vertex is dropped (`cxx:8347-8367`); the survivors are rebound under the new
//! edge together with the faces of the old one (`cxx:8374-8398`), and the
//! `theMECheckExt` / `theMENew` flags migrate to them.
//!
//! `IntersectEdges` intersects a list of splits among themselves to remove
//! self-intersections (`cxx:8423-8427`), updates the images of the trimmed
//! edges, re-images the newly created edges in `theMENew` / `theMECheckExt`
//! (`cxx:8461-8522`), drops the invalid splits found by `GetInvalidEdges`
//! (`cxx:8525-8542`) and finally refreshes the origins (`cxx:8545`).
//!
//! `UpdateNewIntersectionEdges` folds the new splits back into the maps of the
//! class: `myETrimEInf` / `myModifiedEdges` (`cxx:9095-9104`), `myAsDes`
//! (`cxx:9109-9115`), `myOEImages` / `myOEOrigins` (`cxx:9118-9257`) and
//! `myEdgesOrigins` (`cxx:9140-9166`, `cxx:9253-9288`). Two arms differ:
//! a *new* edge binds its images and origins outright, while a *trimmed* edge
//! replaces its old splits in place (`cxx:9194-9225`).

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::bop_builder2::BopBuilder;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::shape::TopoShape;
use crate::topo_tools::{children_keep_orientation, explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::invalid_edges_filter::make_compound_from;
use super::statics::{append_to_list, update_images, update_origins};
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void TrimNewIntersectionEdges(const List& theLE, const DataMap&
    /// theEETrim, const Map& theMVBounds, Map& theMECheckExt, DataMap&
    /// theEImages, Map& theMEB, Map& theMVOld, Map& theMENew, DataMap&
    /// theDMEOr, DataMap& theMELF)` (`cxx:8239-8401`).
    pub fn trim_new_intersection_edges(
        &self,
        the_le: &[TopoShape],
        the_eetrim: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_mv_bounds: &ShapeIndexedMap,
        the_me_check_ext: &mut ShapeIndexedMap,
        the_e_images: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_meb: &mut ShapeIndexedMap,
        the_mv_old: &mut ShapeIndexedMap,
        the_menew: &mut ShapeIndexedMap,
        the_dmeor: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_melf: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) {
        for a_e in the_le {
            // `cxx:8261`: the flag also drops the edge from the set.
            let b_check_ext = the_me_check_ext.remove(a_e);

            // `cxx:8263-8279`: the trimmed parts of the edge and their
            // vertices are the bounding ones.
            let b_old = the_eetrim.contains_key(&shape_map_key(a_e));
            if b_old {
                let a_let = the_eetrim
                    .get(&shape_map_key(a_e))
                    .cloned()
                    .unwrap_or_default();
                for a_et in &a_let {
                    the_meb.add(a_et);
                    for a_v in explorer_shapes(a_et, ShapeType::Vertex) {
                        the_mv_old.add(&a_v);
                    }
                }
            }

            // `cxx:8281-8291`.
            let Some(a_le_im) = the_e_images.get(&shape_map_key(a_e)).cloned() else {
                continue;
            };
            if a_le_im.is_empty() {
                the_e_images.remove(&shape_map_key(a_e));
                continue;
            }

            // `cxx:8293-8344`: fuse the splits of a multi-split edge together
            // with the two bounding vertices.
            let mut a_mev_bounds = ShapeIndexedMap::new();
            let a_ce_im: TopoShape;
            let mut a_mv = ShapeIndexedMap::new();
            let mut a_v1: Option<TopoShape> = None;
            let mut a_v2: Option<TopoShape> = None;
            if a_le_im.len() > 1 {
                let mut a_gfe = BopBuilder::new();
                for a_e_im in &a_le_im {
                    a_gfe.add_argument(a_e_im);
                    for a_v in explorer_shapes(a_e_im, ShapeType::Vertex) {
                        a_mv.add(&a_v);
                    }
                }
                // `cxx:8310-8317`: the two bounding vertices take part in the
                // operation. A degenerate edge answers the same vertex twice
                // in the reference, so the padding repeats it here.
                let a_verts = explorer_shapes(a_e, ShapeType::Vertex);
                if let Some(a_first) = a_verts.first().cloned() {
                    let a_second = a_verts.get(1).cloned().unwrap_or_else(|| a_first.clone());
                    a_v1 = Some(a_first.clone());
                    a_v2 = Some(a_second.clone());
                    a_gfe.add_argument(&a_first);
                    a_gfe.add_argument(&a_second);
                    a_mv.add(&a_first);
                    a_mv.add(&a_second);
                }
                a_gfe.perform().ok();
                if !a_gfe.has_errors() {
                    // `cxx:8326-8337`: the bounding vertices' first image.
                    let a_nb_v = a_mv.extent();
                    for i_v in 1..=a_nb_v {
                        let Some(a_v) = a_mv.value(i_v).cloned() else {
                            continue;
                        };
                        let b_is_own = a_v1.as_ref().is_some_and(|v| is_same(v, &a_v))
                            || a_v2.as_ref().is_some_and(|v| is_same(v, &a_v));
                        if the_mv_bounds.contains(&a_v) || b_is_own {
                            let a_lv_im = a_gfe.history().modified(&a_v).unwrap_or(&[]);
                            let a_v_bound = a_lv_im.first().cloned().unwrap_or_else(|| a_v.clone());
                            a_mev_bounds.add(&a_v_bound);
                        }
                    }
                    a_ce_im = a_gfe.result().clone();
                } else {
                    // The reference leaves `aCEIm` untouched on the error arm
                    // (`cxx:8319`), so the edge walk below sees an empty shape.
                    a_ce_im = make_compound_from(std::iter::empty());
                }
            } else {
                a_ce_im = a_le_im[0].clone();
            }

            // `cxx:8347-8367`: keep only the splits free of bounding vertices.
            let mut a_kept: Vec<TopoShape> = Vec::new();
            for a_e_im in explorer_shapes(&a_ce_im, ShapeType::Edge) {
                let mut b_bound = false;
                for a_v in children_keep_orientation(&a_e_im) {
                    if a_mev_bounds.contains(&a_v) || the_mv_bounds.contains(&a_v) {
                        b_bound = true;
                        break;
                    }
                }
                if b_bound {
                    continue;
                }
                a_kept.push(a_e_im.clone());
                // `cxx:8365`: split -> the untrimmed edge it comes from.
                the_dmeor
                    .entry(shape_map_key(&a_e_im))
                    .or_default()
                    .push(a_e.clone());
            }

            // `cxx:8369-8399`.
            if a_kept.is_empty() {
                the_e_images.remove(&shape_map_key(a_e));
                continue;
            }
            let a_lfe = the_melf
                .get(&shape_map_key(a_e))
                .cloned()
                .unwrap_or_default();
            for a_e_im in &a_kept {
                let a_lfe_im = the_melf
                    .entry(shape_map_key(a_e_im))
                    .or_default();
                for a_f in &a_lfe {
                    append_to_list(a_lfe_im, a_f);
                }
                // `cxx:8393-8398`.
                if b_check_ext {
                    the_me_check_ext.add(a_e_im);
                } else if !b_old {
                    the_menew.add(a_e_im);
                }
            }
            the_e_images.insert(shape_map_key(a_e), a_kept);
        }
    }

    /// `void IntersectEdges(const List& theLA, const List& theLE, const Map&
    /// theMVBounds, const Map& theVertsToAvoid, Map& theMENew, Map&
    /// theMECheckExt, DataMap& theEImages, DataMap& theDMEOr, DataMap& theMELF,
    /// Shape& theSplits)` (`cxx:8408-8545`).
    pub fn intersect_edges(
        &mut self,
        the_la: &[TopoShape],
        the_le: &[TopoShape],
        the_mv_bounds: &ShapeIndexedMap,
        the_verts_to_avoid: &ShapeIndexedMap,
        the_menew: &mut ShapeIndexedMap,
        the_me_check_ext: &mut ShapeIndexedMap,
        the_e_images: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_dmeor: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_melf: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) -> TopoShape {
        // `cxx:8423-8427`.
        let mut a_gfa = BopBuilder::new();
        a_gfa.set_arguments(the_la);
        a_gfa.perform().ok();
        if a_gfa.has_errors() {
            // `cxx:8428-8439`: the input is the result.
            return make_compound_from(the_la.iter());
        }

        // `cxx:8442`.
        update_images(the_le, the_e_images, &a_gfa, &mut self.my_modified_edges);
        // `cxx:8446`.
        let a_splits = a_gfa.result().clone();

        // `cxx:8451-8459`: the edges of the arguments.
        let mut a_le_input: Vec<TopoShape> = Vec::new();
        for a_s in the_la {
            for a_e in explorer_shapes(a_s, ShapeType::Edge) {
                a_le_input.push(a_e);
            }
        }

        // `cxx:8462-8484`: re-image the newly created edges.
        for a_e in &a_le_input {
            if !the_menew.contains(a_e) {
                continue;
            }
            let a_le_im: Vec<TopoShape> = match a_gfa.history().modified(a_e) {
                Some(a_list) if !a_list.is_empty() => a_list.to_vec(),
                _ => continue,
            };
            the_menew.remove(a_e);
            for a_e_im in &a_le_im {
                the_menew.add(a_e_im);
            }
        }

        // `cxx:8486-8522`: the edges under extended check migrate their flag
        // and their faces to the new splits.
        for a_e in &a_le_input {
            let a_le_im: Vec<TopoShape> = match a_gfa.history().modified(a_e) {
                Some(a_list) if !a_list.is_empty() => a_list.to_vec(),
                _ => continue,
            };
            if the_me_check_ext.contains(a_e) {
                for a_e_im in &a_le_im {
                    the_me_check_ext.add(a_e_im);
                }
                the_me_check_ext.remove(a_e);
            }
            let a_lfe = the_melf
                .get(&shape_map_key(a_e))
                .cloned()
                .unwrap_or_default();
            for a_e_im in &a_le_im {
                let a_lfe_im = the_melf.entry(shape_map_key(a_e_im)).or_default();
                for a_f in &a_lfe {
                    append_to_list(a_lfe_im, a_f);
                }
            }
        }

        // `cxx:8524-8543`: drop the invalid splits.
        let mut a_me_inv = ShapeIndexedMap::new();
        self.get_invalid_edges(the_verts_to_avoid, the_mv_bounds, &a_gfa, &mut a_me_inv);
        let a_splits = if a_me_inv.extent() > 0 {
            let a_kept: Vec<TopoShape> = explorer_shapes(&a_splits, ShapeType::Edge)
                .into_iter()
                .filter(|a_e| !a_me_inv.contains(a_e))
                .collect();
            make_compound_from(a_kept.iter())
        } else {
            a_splits
        };

        // `cxx:8545`.
        update_origins(&a_le_input, the_dmeor, &a_gfa);

        a_splits
    }

    /// `void UpdateNewIntersectionEdges(const List& theLE, const DataMap&
    /// theMELF, const DataMap& theEImages, DataMap& theEETrim)`
    /// (`cxx:9048-9292`).
    pub fn update_new_intersection_edges(
        &mut self,
        the_le: &[TopoShape],
        the_melf: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_e_images: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_eetrim: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) {
        for a_e in the_le {
            let a_key = shape_map_key(a_e);

            // `cxx:9065-9091`: an edge without images keeps only the trimmed
            // parts that are still invalid or inverted.
            if !the_e_images.contains_key(&a_key) {
                let Some(a_let) = the_eetrim.get(&a_key).cloned() else {
                    continue;
                };
                let a_kept: Vec<TopoShape> = a_let
                    .into_iter()
                    .filter(|a_et| {
                        self.my_invalid_edges.contains(a_et) || self.my_inverted_edges.contains(a_et)
                    })
                    .collect();
                if a_kept.is_empty() {
                    continue;
                }
                the_eetrim.insert(a_key.clone(), a_kept);
            }

            // `cxx:9093-9094`: the images of the edge, empty when unbound.
            let a_lenew: Vec<TopoShape> = the_e_images
                .get(&a_key)
                .cloned()
                .unwrap_or_default();

            // `cxx:9097-9104`.
            for a_et in &a_lenew {
                self.my_e_trim_e_inf.insert(shape_map_key(a_et), a_e.clone());
                self.my_modified_edges.add(a_et);
            }

            // `cxx:9106-9167`: a new edge.
            if !the_eetrim.contains_key(&a_key) {
                let a_lf = the_melf.get(&a_key).cloned().unwrap_or_default();
                // `cxx:9112-9114`: add the edge to AsDes under every face.
                for a_f in &a_lf {
                    self.my_as_des.add(a_f, a_e);
                }
                // `cxx:9118-9120`.
                self.my_oe_images.bind_list(a_e, &a_lenew);
                self.my_modified_edges.add(a_e);
                // `cxx:9123-9138`.
                for a_e_new in &a_lenew {
                    match self.my_oe_origins.get_mut(&shape_map_key(a_e_new)) {
                        Some(a_origins) => append_to_list(a_origins, a_e),
                        None => {
                            self.my_oe_origins
                                .insert(shape_map_key(a_e_new), vec![a_e.clone()]);
                        }
                    }
                }
                // `cxx:9141-9166`: carry the initial origins over.
                if let Some(a_le_or_init) = self.my_edges_origins.get(&a_key).cloned() {
                    for a_e_new in &a_lenew {
                        match self.my_edges_origins.get_mut(&shape_map_key(a_e_new)) {
                            Some(a_lenew_or) => {
                                for a_e_or in &a_le_or_init {
                                    append_to_list(a_lenew_or, a_e_or);
                                }
                            }
                            None => {
                                self.my_edges_origins
                                    .insert(shape_map_key(a_e_new), a_le_or_init.clone());
                            }
                        }
                    }
                }
                continue;
            }

            // `cxx:9169-9289`: a trimmed edge - the old splits give way to the
            // new ones.
            let a_le_old = the_eetrim.get(&a_key).cloned().unwrap_or_default();
            let mut an_init_origins: Vec<TopoShape> = Vec::new();
            for a_e_old in &a_le_old {
                match self.my_oe_origins.get(&shape_map_key(a_e_old)).cloned() {
                    Some(a_e_origins) => {
                        for a_e_or in &a_e_origins {
                            // `cxx:9198`.
                            self.my_modified_edges.add(a_e_or);
                            // `cxx:9200-9208`: drop the old edge from the
                            // images; the key must be bound.
                            let Some(a_e_images) =
                                self.my_oe_images.change_seek(a_e_or)
                            else {
                                continue;
                            };
                            a_e_images.retain(|a_e_im| !is_same(a_e_im, a_e_old));
                            // `cxx:9210-9226`: append the new images and their
                            // origins.
                            for a_e_new in &a_lenew {
                                append_to_list(a_e_images, a_e_new);
                            }
                        }
                        for a_e_new in &a_lenew {
                            for a_e_or in &a_e_origins {
                                match self.my_oe_origins.get_mut(&shape_map_key(a_e_new)) {
                                    Some(a_new_origins) => {
                                        append_to_list(a_new_origins, a_e_or)
                                    }
                                    None => {
                                        self.my_oe_origins
                                            .insert(shape_map_key(a_e_new), vec![a_e_or.clone()]);
                                    }
                                }
                            }
                        }
                    }
                    None => {
                        // `cxx:9230-9252`.
                        self.my_oe_images.bind_list(a_e_old, &a_lenew);
                        self.my_modified_edges.add(a_e_old);
                        for a_e_new in &a_lenew {
                            match self.my_oe_origins.get_mut(&shape_map_key(a_e_new)) {
                                Some(a_origins) => append_to_list(a_origins, a_e_old),
                                None => {
                                    self.my_oe_origins
                                        .insert(shape_map_key(a_e_new), vec![a_e_old.clone()]);
                                }
                            }
                        }
                    }
                }
                // `cxx:9256-9265`.
                if let Some(a_le_or_init) = self.my_edges_origins.get(&shape_map_key(a_e_old)) {
                    for a_e_or_init in a_le_or_init {
                        append_to_list(&mut an_init_origins, a_e_or_init);
                    }
                }
            }

            // `cxx:9268-9289`.
            if !an_init_origins.is_empty() {
                for a_e_new in &a_lenew {
                    match self.my_edges_origins.get_mut(&shape_map_key(a_e_new)) {
                        Some(a_lenew_or) => {
                            for a_e_or in &an_init_origins {
                                append_to_list(a_lenew_or, a_e_or);
                            }
                        }
                        None => {
                            self.my_edges_origins
                                .insert(shape_map_key(a_e_new), an_init_origins.clone());
                        }
                    }
                }
            }
        }
    }
}
