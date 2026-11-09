//! Edge validity bookkeeping of `BRepOffset_BuildOffsetFaces` - both
//! `FindInvalidEdges` overloads and the file-local `addAsNeutral`.
//!
//! * [`BRepOffsetBuildOffsetFaces::find_invalid_edges`] -
//!   `cxx:2057-2528`. Walks the splits of one offset face and classifies every
//!   edge of every split as valid / invalid / neutral by comparing it with the
//!   original edges it came from.
//! * `addAsNeutral` (`cxx:2534-2561`) - [`add_as_neutral`].
//! * [`BRepOffsetBuildOffsetFaces::find_unclassified_invalid_edges`] -
//!   `cxx:2569-2774`, the second overload, run once over all processed faces to
//!   find the edges left unclassified by the per-face pass.
//!
//! Classification contract of the first overload (the reference's comment at
//! `cxx:2073-2078`):
//!
//! 1. the edge orientation on the split differs from the original edge/face;
//! 2. the vertices of the edge swapped places compared to the original;
//! 3. edges created from vertices (intersection of two faces connected only by
//!    a vertex) get their orientation from the original face's edges adjacent
//!    to that vertex.

use std::collections::HashMap;

use occt_core::gp::GpDir;
use occt_core::precision::{ANGULAR, CONFUSION};

use crate::abs::{Orientation, ShapeType};
use crate::bop_tools_set::BopToolsSet;
use crate::builder::TopoBuilder;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::Dmvv;
use crate::shape::{Compound, Edge, Face, TopoShape};
use crate::topo_tools::{
    children_keep_orientation, explorer_shapes, is_same, map_shapes_and_ancestors, shape_map_key,
    AncestorMap, ShapeMapKey,
};
use crate::topopebrep::build::tools::get_normal_to_face_on_edge;

use super::statics::{append_to_list, find_shape, get_average_tangent, nb_points, AddToContainer};
use super::BRepOffsetBuildOffsetFaces;

/// `static void addAsNeutral(const TopoDS_Shape& theE, const TopoDS_Shape&
/// theFInv, const TopoDS_Shape& theFVal, BRepOffset_DataMapOfShapeIndexedMap&
/// theLocInvEdges, DataMap<Shape, Map>& theLocValidEdges)` (`cxx:2534-2561`).
///
/// Appends `the_e` to `the_f_inv`'s invalid map and to `the_f_val`'s valid map,
/// binding an empty container first when the face is not there yet.
pub fn add_as_neutral(
    the_e: &TopoShape,
    the_f_inv: &TopoShape,
    the_f_val: &TopoShape,
    the_loc_inv_edges: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
    the_loc_valid_edges: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
) {
    the_loc_inv_edges
        .entry(shape_map_key(the_f_inv))
        .or_default()
        .add(the_e);
    the_loc_valid_edges
        .entry(shape_map_key(the_f_val))
        .or_default()
        .add(the_e);
}

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::FindInvalidEdges(const TopoDS_Face&
    /// theF, const List& theLFImages, DataMap& theDMFMVE, DataMap& theDMFMNE,
    /// BRepOffset_DataMapOfShapeIndexedMapOfShape& theDMFMIE, DataMap&
    /// theDMFMVIE, DataMap<Shape, List>& theDMEOrLEIm, Map&
    /// theEdgesInvalidByVertex, Map& theEdgesValidByVertex, const
    /// Message_ProgressRange&)` (`cxx:2057-2528`).
    ///
    /// `the_dmfmve` / `the_dmfmne` / `the_dmfmvie` are the reference's
    /// `DataMap<Shape, Map<Shape>>` (valid, neutral, valid-inverted edges per
    /// split) and `the_dmfmie` its `DataMap<Shape, IndexedMap<Shape>>` (invalid
    /// edges per split); both are carried by
    /// `HashMap<ShapeMapKey, ShapeIndexedMap>` because the per-key container is
    /// the only thing the consumers read.
    pub fn find_invalid_edges(
        &mut self,
        the_f: &Face,
        the_lf_images: &[TopoShape],
        the_dmfmve: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_dmfmne: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_dmfmie: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_dmfmvie: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_dme_or_le_im: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_edges_invalid_by_vertex: &mut ShapeIndexedMap,
        the_edges_valid_by_vertex: &mut ShapeIndexedMap,
    ) {
        // `cxx:2085-2086`: the original face of the offset face.
        let a_f_or = match self.my_faces_origins.get(&shape_map_key(&the_f.0)) {
            Some(a_f) => match Face::wrap(a_f.clone()) {
                Some(a_f) => a_f,
                None => return,
            },
            // `myFacesOrigins->Find` asserts on an absent key in OCCT.
            None => return,
        };

        // `cxx:2088-2098`.
        let mut a_me_inv = ShapeIndexedMap::new();
        let mut a_me_val = ShapeIndexedMap::new();
        let mut a_me_int = ShapeIndexedMap::new();
        let mut a_dmve = Dmvv::new();
        let mut a_dmef = Dmvv::new();
        let mut a_m_edges = ShapeIndexedMap::new();
        let mut an_images: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        // `cxx:2102-2161`: first pass - the split topology maps.
        for a_f_im_shape in the_lf_images {
            let Some(a_f_im) = Face::wrap(a_f_im_shape.clone()) else {
                continue;
            };
            for a_e in explorer_shapes(&a_f_im.0, ShapeType::Edge) {
                // `cxx:2116`.
                a_m_edges.add(&a_e);

                // `cxx:2118-2122`: edge -> split (AppendToList de-duplicates).
                match a_dmef.change_seek(&a_e) {
                    Some(p_lf) => append_to_list(p_lf, &a_f_im.0),
                    None => a_dmef.add(&a_e, &a_f_im.0),
                }

                // `cxx:2124-2138`: `TopoDS_Iterator(aE)` - the edge's vertices,
                // one entry per vertex occurrence.
                for a_v in children_keep_orientation(&a_e) {
                    match a_dmve.change_seek(&a_v) {
                        Some(p_le) => append_to_list(p_le, &a_e),
                        None => a_dmve.add(&a_v, &a_e),
                    }
                }

                // `cxx:2140-2160`: back map, original edge in the original face
                // -> its offset images.
                let Some(p_l_or) = self.my_edges_origins.get(&shape_map_key(&a_e)).cloned() else {
                    continue;
                };
                for a_s_or in &p_l_or {
                    let mut a_s_in_f = TopoShape::new(a_s_or.shape_type());
                    if !find_shape(a_s_or, &a_f_or.0, self.my_analyzer.as_ref(), &mut a_s_in_f) {
                        continue;
                    }
                    let p_images = an_images.entry(shape_map_key(&a_s_in_f)).or_default();
                    append_to_list(p_images, &a_e);
                }
            }
        }

        // `cxx:2163-2168`: the original face's vertex -> edge map, filled on
        // first necessity.
        let mut a_dmvef_or: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        // `cxx:2169-2466`: second pass - classify the edges of every split.
        for a_f_im_shape in the_lf_images {
            let Some(a_f_im) = Face::wrap(a_f_im_shape.clone()) else {
                continue;
            };
            // `cxx:2178-2180`.
            let mut a_mve = ShapeIndexedMap::new();
            let mut a_mie = ShapeIndexedMap::new();

            for a_e_im_shape in explorer_shapes(&a_f_im.0, ShapeType::Edge) {
                let Some(a_e_im) = Edge::wrap(a_e_im_shape) else {
                    continue;
                };
                // `cxx:2188-2192`.
                if a_e_im.orientation() == Orientation::Internal {
                    a_me_int.add(&a_e_im.0);
                    continue;
                }
                // `cxx:2194-2198`.
                let Some(p_le_or) = self
                    .my_edges_origins
                    .get(&shape_map_key(&a_e_im.0))
                    .cloned()
                else {
                    continue;
                };
                if p_le_or.is_empty() {
                    continue;
                }
                // `cxx:2200-2207`.
                let a_nb_v_or = p_le_or
                    .iter()
                    .filter(|s| s.shape_type() == ShapeType::Vertex)
                    .count();
                if a_nb_v_or > 1 && (p_le_or.len() - a_nb_v_or) > 1 {
                    continue;
                }

                // `cxx:2212-2220`.
                let mut a_me = ShapeIndexedMap::new();
                let mut a_mv = ShapeIndexedMap::new();
                let mut a_mf = ShapeIndexedMap::new();
                let mut b_invalid = false;
                let mut b_checked = false;
                let a_nb_p = nb_points(&a_e_im);
                let mut a_nb_inv = 0usize;
                let b_use_vertex = a_nb_v_or == 1
                    && a_dmef
                        .seek(&a_e_im.0)
                        .map(|l| l.len())
                        .unwrap_or(0)
                        == 1
                    && !self
                        .my_oe_origins
                        .contains_key(&shape_map_key(&a_e_im.0));

                // `cxx:2222-2378`: compare with every origin of the edge.
                for a_s_or in &p_le_or {
                    let b_vertex = a_s_or.shape_type() == ShapeType::Vertex;
                    // `aEOrF` is null until an original edge (or the compound of
                    // the adjacent original edges) is found.
                    let mut a_e_or_f = TopoShape::new(ShapeType::Edge);
                    let mut b_found = false;

                    if b_vertex {
                        // `cxx:2229-2233`.
                        if !b_use_vertex {
                            continue;
                        }
                        // `cxx:2234-2244`.
                        if a_dmvef_or.is_empty() {
                            a_dmvef_or = map_shapes_and_ancestors(
                                &a_f_or.0,
                                ShapeType::Vertex,
                                ShapeType::Edge,
                            );
                        }
                        if let Some(p_le_f_or) = a_dmvef_or.get(&shape_map_key(a_s_or)).cloned() {
                            // `cxx:2245-2266`.
                            let a_bb = TopoBuilder::new();
                            let mut a_ce_or = Compound::new();
                            let mut b_all_tgt = true;
                            let a_v_ref = get_average_tangent(&p_le_f_or[0], a_nb_p);
                            for a_e_or in &p_le_f_or {
                                a_bb.add(&mut a_ce_or.0, a_e_or);
                                let a_v_cur = get_average_tangent(a_e_or, a_nb_p);
                                if !a_v_ref.is_parallel_ang(&a_v_cur, ANGULAR) {
                                    b_all_tgt = false;
                                }
                            }
                            // `cxx:2263-2266`: a vertex between tangent edges
                            // stays unusable, so `aEOrF` remains null.
                            if !b_all_tgt {
                                a_e_or_f = a_ce_or.0;
                                b_found = true;
                            }
                        }
                    } else {
                        // `cxx:2269-2279`.
                        b_found = find_shape(
                            a_s_or,
                            &a_f_or.0,
                            self.my_analyzer.as_ref(),
                            &mut a_e_or_f,
                        );
                        let p_le_im = the_dme_or_le_im
                            .entry(shape_map_key(a_s_or))
                            .or_default();
                        append_to_list(p_le_im, &a_e_im.0);
                    }

                    // `cxx:2281-2285`.
                    if !b_found {
                        continue;
                    }

                    if b_vertex {
                        // `cxx:2287-2333`: the original edges sharing the vertex
                        // must keep sharing it.
                        let mut a_mv_total = ShapeIndexedMap::new();
                        let mut a_nb_checked = 0usize;
                        for a_e_or in children_keep_orientation(&a_e_or_f) {
                            let Some(a_l_im) =
                                an_images.get(&shape_map_key(&a_e_or)).cloned()
                            else {
                                continue;
                            };
                            a_nb_checked += 1;
                            // `cxx:2303-2308`: one map fed by every image, so a
                            // vertex met in two images counts twice.
                            let mut a_mv_loc: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
                            let mut a_mv_loc_index: HashMap<ShapeMapKey, usize> = HashMap::new();
                            for an_im in &a_l_im {
                                let a_loc = AncestorMap::map_shapes_and_ancestors(
                                    an_im,
                                    ShapeType::Vertex,
                                    ShapeType::Edge,
                                );
                                for j in 0..a_loc.length() {
                                    let a_v = a_loc.find_key(j);
                                    match a_mv_loc_index.get(&shape_map_key(a_v)) {
                                        Some(&k) => {
                                            a_mv_loc[k].1.extend(a_loc.ancestors(j).iter().cloned())
                                        }
                                        None => {
                                            a_mv_loc_index.insert(
                                                shape_map_key(a_v),
                                                a_mv_loc.len(),
                                            );
                                            a_mv_loc
                                                .push((a_v.clone(), a_loc.ancestors(j).to_vec()));
                                        }
                                    }
                                }
                            }
                            // `cxx:2309-2321`.
                            for (a_v, a_le) in &a_mv_loc {
                                if a_le.len() > 1 && !a_mv_total.add_to_container(a_v) {
                                    b_invalid = true;
                                    the_edges_invalid_by_vertex.add(&a_e_im.0);
                                    break;
                                }
                            }
                            if b_invalid {
                                break;
                            }
                        }
                        // `cxx:2322-2331`.
                        if !b_invalid && a_nb_checked < 2 {
                            continue;
                        }
                        the_edges_valid_by_vertex.add(&a_e_im.0);
                    } else {
                        // `cxx:2334-2375`: the same 3d direction means the same
                        // orientation; a reversed direction means inverted.
                        let mut a_v_sum1 = get_average_tangent(&a_e_im.0, a_nb_p);
                        let mut a_v_sum2 = get_average_tangent(&a_e_or_f, a_nb_p);
                        a_v_sum1.normalize();
                        a_v_sum2.normalize();
                        let a_cos = a_v_sum1.dot(&a_v_sum2);
                        if a_cos.abs() < 0.9999 {
                            continue;
                        }
                        // `cxx:2352-2360`.
                        a_me.add(&a_e_or_f);
                        for a_v in explorer_shapes(&a_e_or_f, ShapeType::Vertex) {
                            a_mv.add(&a_v);
                        }
                        // `cxx:2361-2369`.
                        if let Some(a_analyzer) = self.my_analyzer.as_ref() {
                            for a_f in a_analyzer.ancestors(&a_e_or_f).to_vec() {
                                a_mf.add(&a_f);
                            }
                        }
                        // `cxx:2371-2375`.
                        if a_cos < CONFUSION {
                            b_invalid = true;
                            a_nb_inv += 1;
                        }
                    }
                    b_checked = true;
                }

                // `cxx:2380-2383`.
                if !b_checked {
                    continue;
                }

                // `cxx:2385-2386`.
                let mut b_local_only = a_nb_v_or > 1 && (p_le_or.len() - a_nb_v_or) > 1;
                let a_nb_e = a_me.extent();
                let a_nb_v = a_mv.extent();

                // `cxx:2387-2428`.
                if a_nb_e > 1 && a_nb_v == 2 * a_nb_e {
                    let mut b_skip = true;
                    if a_mf.extent() > 2 && (a_nb_inv == 0 || a_nb_inv == a_nb_e) {
                        if the_lf_images.len() > 2 {
                            let mut b_broken = false;
                            for a_v in children_keep_orientation(&a_e_im.0) {
                                let Some(a_le_v) = a_dmve.seek(&a_v).map(|l| l.to_vec()) else {
                                    // `aDMVE.FindFromKey(itV.Value())` asserts on
                                    // an unbound vertex in OCCT.
                                    continue;
                                };
                                let mut b_edge_short = false;
                                for a_e in &a_le_v {
                                    let a_nb_f = a_dmef.seek(a_e).map(|l| l.len()).unwrap_or(0);
                                    if a_nb_f < 2 {
                                        b_edge_short = true;
                                        break;
                                    }
                                }
                                if b_edge_short {
                                    b_broken = true;
                                    break;
                                }
                            }
                            // `cxx:2418`: `bSkip = itV.More()`.
                            b_skip = b_broken;
                        }
                    }
                    // `cxx:2420-2427`.
                    if b_skip {
                        continue;
                    }
                    b_local_only = true;
                }

                // `cxx:2430-2439`.
                if b_invalid {
                    if !b_local_only {
                        self.my_invalid_edges.add(&a_e_im.0);
                    }
                    a_mie.add(&a_e_im.0);
                    a_me_inv.add(&a_e_im.0);
                    continue;
                }

                // `cxx:2441-2442`.
                let b_inverted = if a_nb_e == 0 || b_local_only {
                    false
                } else {
                    self.check_inverted(&a_e_im, &a_f_or, &a_dmve, &a_m_edges)
                };

                // `cxx:2444-2453`.
                if !b_inverted || a_nb_v_or == 0 {
                    if !b_local_only {
                        self.my_valid_edges.add(&a_e_im.0);
                    }
                    a_mve.add(&a_e_im.0);
                    a_me_val.add(&a_e_im.0);
                }
            }

            // `cxx:2456-2464`.
            if a_mve.extent() > 0 {
                the_dmfmve.insert(shape_map_key(&a_f_im.0), a_mve);
            }
            if a_mie.extent() > 0 {
                the_dmfmie.insert(shape_map_key(&a_f_im.0), a_mie);
            }
        }

        // `cxx:2468-2527`: the invalid edges that are also images of inverted
        // edges make their inverted images invalid too.
        let mut a_mvie = ShapeIndexedMap::new();
        let mut a_mne = ShapeIndexedMap::new();
        for a_nb_e_inv in 1..=a_me_inv.extent() {
            let Some(a_e_im) = a_me_inv.value(a_nb_e_inv).cloned() else {
                continue;
            };
            // `cxx:2480-2484`: neutral edges - valid in one split, invalid in
            // another.
            if a_me_val.contains(&a_e_im) {
                a_mne.add(&a_e_im);
                continue;
            }
            // `cxx:2487-2490`.
            if !self.my_inverted_edges.contains(&a_e_im) {
                continue;
            }
            // `cxx:2492-2496`.
            let Some(p_l_oe_or) = self.my_oe_origins.get(&shape_map_key(&a_e_im)).cloned() else {
                continue;
            };
            // `cxx:2498-2515`.
            for a_oe_or in &p_l_oe_or {
                let Some(a_le_im1) = self.my_oe_images.seek(a_oe_or).map(|l| l.to_vec()) else {
                    // `myOEImages.Find(aOEOr)` asserts on an unbound key.
                    continue;
                };
                for a_e_im1 in &a_le_im1 {
                    if a_m_edges.contains(a_e_im1)
                        && !a_me_inv.contains(a_e_im1)
                        && !a_me_int.contains(a_e_im1)
                        && self.my_inverted_edges.contains(a_e_im1)
                    {
                        self.my_invalid_edges.add(a_e_im1);
                        a_mvie.add(a_e_im1);
                    }
                }
            }
        }

        // `cxx:2516-2527`.
        if a_mne.extent() > 0 {
            the_dmfmne.insert(shape_map_key(&the_f.0), a_mne);
        }
        if a_mvie.extent() > 0 {
            the_dmfmvie.insert(shape_map_key(&the_f.0), a_mvie);
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::FindInvalidEdges(const List&
    /// theLFOffset, BRepOffset_DataMapOfShapeIndexedMapOfShape& theLocInvEdges,
    /// DataMap& theLocValidEdges, DataMap& theNeutralEdges)` (`cxx:2569-2774`).
    ///
    /// The additional pass over all processed faces:
    ///
    /// 1. find the edges classified as valid/invalid somewhere but not in the
    ///    current split (`cxx:2580`, `:2634-2650`);
    /// 2. for every such edge take the other splits where the *same* edge set
    ///    lives (`cxx:2686-2703`) and compare the face normal on the edge;
    /// 3. the split where the edge orientation flips marks it neutral, and
    ///    `addAsNeutral` moves it between the invalid/valid maps
    ///    (`cxx:2704-2758`).
    pub fn find_unclassified_invalid_edges(
        &mut self,
        the_lf_offset: &[TopoShape],
        the_loc_inv_edges: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_loc_valid_edges: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
        the_neutral_edges: &mut HashMap<ShapeMapKey, ShapeIndexedMap>,
    ) {
        // `cxx:2583-2588`: edge -> splits where it was not classified, and
        // split -> the offset face it belongs to.
        let mut a_me_unclassified = Dmvv::new();
        let mut a_f_split_f_offset: HashMap<ShapeMapKey, TopoShape> = HashMap::new();

        // `cxx:2589-2608`: the offset images of the analyzer's new faces are
        // skipped (they are artificial).
        let mut a_new_faces = ShapeIndexedMap::new();
        if let Some(a_analyzer) = self.my_analyzer.as_ref() {
            let mut a_map_new_tmp = ShapeIndexedMap::new();
            for a_f in a_analyzer.new_faces() {
                a_map_new_tmp.add(a_f);
            }
            for a_f_offset in the_lf_offset {
                let Some(a_f_origin) = self.my_faces_origins.get(&shape_map_key(a_f_offset)) else {
                    // `myFacesOrigins->Find` asserts on an absent key.
                    continue;
                };
                if a_map_new_tmp.contains(a_f_origin) {
                    a_new_faces.add(a_f_offset);
                }
            }
        }

        // `cxx:2610-2659`: the first pass.
        let mut an_ef_map = Dmvv::new();
        for a_f in the_lf_offset {
            if a_new_faces.contains(a_f) {
                continue;
            }
            let Some(a_lf_images) = self.my_of_images.seek(a_f).map(|l| l.to_vec()) else {
                // `myOFImages.FindFromKey(aF)` asserts on an unbound key.
                continue;
            };
            for a_f_im in &a_lf_images {
                // `cxx:2625`: `MapShapesAndAncestors(aFIm, EDGE, FACE, anEFMap)`
                // accumulates into the map shared by all splits.
                for a_face in explorer_shapes(a_f_im, ShapeType::Face) {
                    for a_e in explorer_shapes(&a_face, ShapeType::Edge) {
                        match an_ef_map.change_seek(&a_e) {
                            Some(p_lf) => append_to_list(p_lf, &a_face),
                            None => an_ef_map.add(&a_e, &a_face),
                        }
                    }
                }

                // `cxx:2627-2630`: the per-split classification made above.
                let p_me_invalid = the_loc_inv_edges.get(&shape_map_key(a_f_im));
                let p_me_valid = the_loc_valid_edges.get(&shape_map_key(a_f_im));

                // `cxx:2632-2652`.
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    if self.my_invalid_edges.contains(&a_e) == self.my_valid_edges.contains(&a_e) {
                        continue;
                    }
                    let b_in_invalid = p_me_invalid.map_or(false, |m| m.contains(&a_e));
                    let b_in_valid = p_me_valid.map_or(false, |m| m.contains(&a_e));
                    if b_in_invalid || b_in_valid {
                        continue;
                    }
                    match a_me_unclassified.change_seek(&a_e) {
                        Some(p_map) => append_to_list(p_map, a_f_im),
                        None => a_me_unclassified.add(&a_e, a_f_im),
                    }
                    a_f_split_f_offset.insert(shape_map_key(a_f_im), a_f.clone());
                }
            }
        }

        // `cxx:2658-2661`.
        if a_me_unclassified.extent() == 0 {
            return;
        }

        // `cxx:2663-2767`: analyze the unclassified edges.
        let a_nb_e = a_me_unclassified.extent();
        for i_e in 1..=a_nb_e {
            let a_e = a_me_unclassified.keys()[i_e - 1].clone();
            let Some(a_mf_unclassified) =
                a_me_unclassified.seek(&a_e).map(|l| l.to_vec())
            else {
                continue;
            };
            // `cxx:2670-2671`.
            let Some(a_lf) = an_ef_map.seek(&a_e).map(|l| l.to_vec()) else {
                // `anEFMap.FindFromKey(aE)` asserts on an unbound key.
                continue;
            };

            // `cxx:2673-2763`.
            for a_f_classified in &a_lf {
                if a_mf_unclassified
                    .iter()
                    .any(|a_f| is_same(a_f, a_f_classified))
                {
                    continue;
                }

                // `cxx:2682-2684`: the edge set of the classified split.
                let an_edge_set_class = BopToolsSet::from_shape(a_f_classified, ShapeType::Edge);

                // `cxx:2686-2688`.
                let mut a_e_classified = TopoShape::new(ShapeType::Edge);
                find_shape(&a_e, a_f_classified, None, &mut a_e_classified);
                let an_ori_class = a_e_classified.orientation();

                // `cxx:2689-2691`: `TopoDS::Edge` / `TopoDS::Face` casts. OCCT
                // throws on a wrong type here; the port skips the split.
                let (Some(a_e_classified_edge), Some(a_f_classified_face)) = (
                    Edge::wrap(a_e_classified.clone()),
                    Face::wrap(a_f_classified.clone()),
                ) else {
                    continue;
                };
                let Some(a_dn_class) =
                    get_normal_to_face_on_edge(&a_f_classified_face, &a_e_classified_edge)
                else {
                    // `BRepGProp_Face` would throw in OCCT; without the normal
                    // no orientation can be compared.
                    continue;
                };

                // `cxx:2692-2694`.
                let p_me_invalid = the_loc_inv_edges.get(&shape_map_key(a_f_classified));
                let b_invalid = p_me_invalid.map_or(false, |m| m.contains(&a_e));

                // `cxx:2696-2759`.
                for a_f_unclassified in &a_mf_unclassified {
                    let an_edge_set_unclass =
                        BopToolsSet::from_shape(a_f_unclassified, ShapeType::Edge);
                    // `cxx:2706`.
                    if !(an_edge_set_class == an_edge_set_unclass) {
                        continue;
                    }

                    // `cxx:2708-2712`: the normals' agreement decides whether
                    // the orientations are directly comparable.
                    let Some(a_f_unclassified_face) = Face::wrap(a_f_unclassified.clone()) else {
                        continue;
                    };
                    let Some(a_dn_unclass) =
                        get_normal_to_face_on_edge(&a_f_unclassified_face, &a_e_classified_edge)
                    else {
                        continue;
                    };
                    let b_same_ori = match (GpDir::from_vec(&a_dn_class), GpDir::from_vec(&a_dn_unclass))
                    {
                        (Ok(a_dn_class), Ok(a_dn_unclass)) => {
                            a_dn_class.is_equal_tol(&a_dn_unclass, ANGULAR)
                        }
                        _ => false,
                    };

                    // `cxx:2715-2717`.
                    let Some(a_f_offset) =
                        a_f_split_f_offset.get(&shape_map_key(a_f_unclassified)).cloned()
                    else {
                        continue;
                    };
                    let Some(a_lf_splits) = self.my_of_images.seek(&a_f_offset).map(|l| l.to_vec())
                    else {
                        // `myOFImages.FindFromKey(aFOffset)` asserts.
                        continue;
                    };

                    // `cxx:2718-2760`: the other splits of the same offset face.
                    // The reference's trailing `if (itLFSp.More()) break;`
                    // (`cxx:2760`) cannot fire - the loop has no `break`.
                    for a_f_sp in &a_lf_splits {
                        if is_same(a_f_sp, a_f_unclassified)
                            || !a_mf_unclassified.iter().any(|a_f| is_same(a_f, a_f_sp))
                        {
                            continue;
                        }
                        let mut a_e_unclassified = TopoShape::new(ShapeType::Edge);
                        find_shape(&a_e, a_f_sp, None, &mut a_e_unclassified);
                        let mut an_ori_unclass = a_e_unclassified.orientation();
                        if !b_same_ori {
                            an_ori_unclass = an_ori_unclass.reversed();
                        }
                        // `cxx:2735-2759`.
                        if an_ori_class != an_ori_unclass {
                            the_neutral_edges
                                .entry(shape_map_key(&a_f_offset))
                                .or_default()
                                .add(&a_e);
                            if b_invalid && b_same_ori {
                                // invalid in aFUnclassified, valid in aFSp.
                                add_as_neutral(
                                    &a_e,
                                    a_f_classified,
                                    a_f_sp,
                                    the_loc_inv_edges,
                                    the_loc_valid_edges,
                                );
                            } else {
                                // invalid in aFSp, valid in aFUnclassified.
                                add_as_neutral(
                                    &a_e,
                                    a_f_sp,
                                    a_f_classified,
                                    the_loc_inv_edges,
                                    the_loc_valid_edges,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
