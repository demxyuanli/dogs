//! `BuildSplitsOfFaces` (`cxx:1342-1833`): the driver that splits every offset
//! face, classifies the resulting edges as valid / invalid / neutral, drops the
//! invalid splits and finally leaves `myInvalidFaces` / `myInvalidEdges` /
//! `myEdgesToAvoid` in the state `BuildSplitsOfExtendedFaces` retries on.
//!
//! Stages, in the reference's order:
//!
//! 1. `cxx:1377-1548`: per face of `myFaces`, `GetEdges` collects the splitting
//!    edges; the invalid ones (`aMapEInv`) make the first split attempt
//!    suspicious, so the face is split again over `aCE + aMapEInv`
//!    (`cxx:1423-1432`) and the pieces made only of those edges are dropped
//!    (`cxx:1435-1476`). When the face stays artificial it is recorded as such
//!    and its images are all taken as invalid (`cxx:1490-1521`). Otherwise
//!    `FindInvalidEdges` classifies every edge of every split
//!    (`cxx:1525-1534`) and the splits are stored under the offset face.
//! 2. `cxx:1550-1561`: with nothing invalid found, the driver stops. Otherwise
//!    the second `FindInvalidEdges` pass looks at the edges left unclassified by
//!    the first one (`cxx:1557`), and `MakeInvertedEdgesInvalid` marks the
//!    inverted edges sitting inside loops of invalid edges (`cxx:1561`).
//! 3. `cxx:1600-1691`: an `(EDGE, FACE)` map over all splits feeds
//!    `FindFacesInsideHoleWires` / `FindInvalidFaces`, which fill
//!    `myInvalidFaces`, the hole-face compound and the inverted-face list.
//! 4. `cxx:1693-1735`: nothing invalid -> clear `myInvalidEdges` and stop. Then
//!    the invalid splits built from inverted edges are removed
//!    (`RemoveInvalidSplitsByInvertedEdges`), the invalid splits are removed
//!    from the valid lists (`RemoveInvalidSplitsFromValid`) and the faces inside
//!    the volume are removed (`RemoveInsideFaces`).
//! 5. `cxx:1738-1794`: the surviving splits become one compound; its
//!    `(EDGE, FACE)` map drives `FilterEdgesImages`, `FilterInvalidFaces` and
//!    `FilterInvalidEdges`, and `CheckEdgesCreatedByVertex` re-checks the edges
//!    created from a vertex.
//! 6. `cxx:1808-1832`: `myLastInvEdges` takes every invalid edge, an invalid
//!    edge the rebuild does not use also goes to `myEdgesToAvoid`, and the
//!    per-face rebuild counter `myAlreadyInvFaces` is bumped.
//!
//! `Message_ProgressScope` (`cxx:1377`, `:1379`, `:1624`) is the only use of
//! `theRange`; the port has no progress reporting, so the three
//! `if (!aPS*.More()) return;` arms (`cxx:1384`, `:1631`) are absent. Every
//! other branch is kept.
//!
//! Two reference details are reproduced explicitly:
//!
//! * the loop at `cxx:1435-1476` removes from `aLFImages1` while iterating and
//!   only advances the iterator in the keep arm, so the port walks the `Vec` by
//!   index and does not increment after a removal;
//! * `checkConnectionsOfFace` (`cxx:1461`) reads `aLFImages1` *in its current,
//!   partly filtered state*, so the port hands it the live list.
//!
//! `pLFIm` (`cxx:1389`) can be null in the reference; the two writes that follow
//! (`cxx:1495`, `cxx:1539-1545`) then dereference it. The port writes through
//! `Dmvv::bind_list`, which binds the key when absent, so the unbound case is
//! safe here instead of crashing.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{
    explorer_shapes, map_shapes_and_ancestors, map_shapes_and_ancestors_into, shape_map_key,
    ShapeMapKey,
};

use super::invalid_edges_filter::make_compound_from;
use super::splits::build_splits_of_face;
use super::statics::{check_connections_of_face, has_data};
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::BuildSplitsOfFaces(const
    /// Message_ProgressRange&)` (`cxx:1342-1833`).
    pub fn build_splits_of_faces(&mut self) {
        // `cxx:1348-1374`: the maps the classification fills.
        let mut a_lf_done: Vec<TopoShape> = Vec::new();
        // Extended face -> neutral edges (valid in one split, invalid in the
        // other).
        let mut a_dmfmne: HashMap<ShapeMapKey, ShapeIndexedMap> = HashMap::new();
        // Offset face -> its valid edges.
        let mut a_dmfmve: HashMap<ShapeMapKey, ShapeIndexedMap> = HashMap::new();
        // Offset face -> its invalid edges.
        let mut a_dmfmie: HashMap<ShapeMapKey, ShapeIndexedMap> = HashMap::new();
        // Offset face -> its valid inverted edges.
        let mut a_dmfmvie: HashMap<ShapeMapKey, ShapeIndexedMap> = HashMap::new();
        // Splits to check for internals.
        let mut a_mf_to_check_int = ShapeIndexedMap::new();
        // Edges created from a vertex and marked invalid / valid.
        let mut a_me_edge_invalid_by_vertex = ShapeIndexedMap::new();
        let mut a_me_edge_valid_by_vertex = ShapeIndexedMap::new();
        // Original edge -> its splits.
        let mut a_dme_or_le_im: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        // `cxx:1380-1548`: build the splits of every face.
        let a_faces = self.my_faces.clone();
        for a_f_shape in a_faces {
            let Some(a_f) = Face::wrap(a_f_shape) else {
                continue;
            };

            // `cxx:1389-1393`: an already bound key with an empty list is done.
            let b_bound = self.my_of_images.contains(&a_f.0);
            let b_empty = self
                .my_of_images
                .seek(&a_f.0)
                .map(|l| l.is_empty())
                .unwrap_or(false);
            if b_bound && b_empty {
                continue;
            }

            // `cxx:1395-1401`: the edges by which the face must be split.
            let mut a_map_e_inv = ShapeIndexedMap::new();
            let (a_ce, b_found) = self.get_edges(&a_f, Some(&mut a_map_e_inv));
            if !b_found {
                continue;
            }

            // `cxx:1412-1413`.
            let mut a_lf_images: Vec<TopoShape> = Vec::new();
            build_splits_of_face(&a_f, &a_ce, &mut self.my_faces_origins, &mut a_lf_images);

            // `cxx:1415-1522`: the invalid splitting edges are present.
            if a_map_e_inv.extent() != 0 {
                // `cxx:1418-1420`: no splits at all, or the splits do not cover
                // every possible piece.
                let mut a_men_inv = ShapeIndexedMap::new();
                let mut b_artificial_case = a_lf_images.is_empty()
                    || self.check_if_artificial(
                        &a_f.0,
                        &a_lf_images,
                        &a_ce,
                        &a_map_e_inv,
                        &mut a_men_inv,
                    );

                // `cxx:1423-1432`: split again over `aCE` plus the invalid edges.
                let mut a_ce1_children: Vec<TopoShape> = Vec::new();
                a_ce1_children.push(a_ce.clone());
                for i in 1..=a_map_e_inv.extent() {
                    if let Some(a_e) = a_map_e_inv.value(i) {
                        a_ce1_children.push(a_e.clone());
                    }
                }
                let a_ce1 = make_compound_from(a_ce1_children.iter());
                let mut a_lf_images1: Vec<TopoShape> = Vec::new();
                build_splits_of_face(
                    &a_f,
                    &a_ce1,
                    &mut self.my_faces_origins,
                    &mut a_lf_images1,
                );

                // `cxx:1435-1476`: drop the pieces made only of invalid edges;
                // the ones that still connect to a neighbour stay.
                let mut i = 0;
                while i < a_lf_images1.len() {
                    let a_f_im = a_lf_images1[i].clone();
                    let mut b_all_inv = true;
                    // `cxx:1441`: the re-check starts from the artificial case.
                    let mut a_to_recheck_face = b_artificial_case;
                    // `!aExpE.More()` of the reference: true only when no arm
                    // broke out early.
                    let mut b_explored_all = true;
                    for a_e in explorer_shapes(&a_f_im, ShapeType::Edge) {
                        if !a_map_e_inv.contains(&a_e) {
                            b_all_inv = false;
                            if !a_men_inv.contains(&a_e) {
                                a_to_recheck_face = false;
                                b_explored_all = false;
                                break;
                            }
                        }
                    }
                    // `cxx:1459-1462`: the number of connections decides.
                    if b_explored_all && a_to_recheck_face {
                        a_to_recheck_face = check_connections_of_face(&a_f_im, &a_lf_images1);
                    }
                    // `cxx:1464-1475`.
                    if b_explored_all && !a_to_recheck_face {
                        if b_all_inv {
                            a_mf_to_check_int.add(&a_f_im);
                        }
                        // The reference's `Remove(it)` leaves the iterator on
                        // the next element, so the index is not advanced.
                        a_lf_images1.remove(i);
                    } else {
                        i += 1;
                    }
                }

                // `cxx:1478-1488`.
                if b_artificial_case {
                    if a_lf_images.len() == a_lf_images1.len() {
                        b_artificial_case = false;
                    } else {
                        a_lf_images = a_lf_images1;
                    }
                }

                // `cxx:1490-1521`: the face stays artificial - every edge of
                // the invalid set makes both the edge and the face invalid.
                if b_artificial_case {
                    let mut a_me_inv = ShapeIndexedMap::new();
                    // `*pLFIm = aLFImages` (`cxx:1495`).
                    self.my_of_images.bind_list(&a_f.0, &a_lf_images);
                    for a_f_im in &a_lf_images {
                        for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                            if a_map_e_inv.contains(&a_e) {
                                self.my_invalid_edges.add(&a_e);
                                a_me_inv.add(&a_e);
                            } else {
                                self.my_valid_edges.add(&a_e);
                            }
                        }
                    }
                    // `cxx:1516-1518`.
                    self.bind_art_invalid_face(&a_f.0, a_me_inv.clone());
                    a_dmfmie.insert(shape_map_key(&a_f.0), a_me_inv);
                    a_lf_done.push(a_f.0.clone());
                    continue;
                }
            }

            // `cxx:1525-1534`: classify the edges of every split.
            self.find_invalid_edges(
                &a_f,
                &a_lf_images,
                &mut a_dmfmve,
                &mut a_dmfmne,
                &mut a_dmfmie,
                &mut a_dmfmvie,
                &mut a_dme_or_le_im,
                &mut a_me_edge_invalid_by_vertex,
                &mut a_me_edge_valid_by_vertex,
            );

            // `cxx:1537-1547`: store the splits under the offset face.
            self.my_of_images.bind_list(&a_f.0, &a_lf_images);
            a_lf_done.push(a_f.0.clone());
        }

        // `cxx:1550-1553`: nothing invalid - nothing to rebuild.
        if self.my_invalid_edges.extent() == 0
            && self.my_art_invalid_faces.is_empty()
            && a_dmfmie.is_empty()
        {
            return;
        }

        // `cxx:1557`: the edges left unclassified by the first pass.
        self.find_unclassified_invalid_edges(
            &a_lf_done,
            &mut a_dmfmie,
            &mut a_dmfmve,
            &mut a_dmfmne,
        );

        // `cxx:1561`: the inverted edges inside loops of invalid edges.
        self.make_inverted_edges_invalid(&a_lf_done);

        // `cxx:1600-1613`: `(EDGE, FACE)` over every split of every offset face.
        let mut an_ef_map: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let a_nb_f = self.my_of_images.extent();
        for i in 1..=a_nb_f {
            let Some(a_f) = self.my_of_images.keys().get(i - 1).cloned() else {
                continue;
            };
            let a_lf_im: Vec<TopoShape> = self
                .my_of_images
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for a_f_im in &a_lf_im {
                map_shapes_and_ancestors_into(
                    a_f_im,
                    ShapeType::Edge,
                    ShapeType::Face,
                    &mut an_ef_map,
                );
            }
        }

        // `cxx:1615-1622`.
        let mut a_mf_inv_in_hole = ShapeIndexedMap::new();
        let mut a_f_holes_children: Vec<TopoShape> = Vec::new();
        let mut an_inverted_faces: Vec<TopoShape> = Vec::new();

        // `cxx:1627-1691`: find the invalid faces.
        let a_lf_done_snapshot = a_lf_done.clone();
        for a_f_done in a_lf_done_snapshot {
            let Some(a_f) = Face::wrap(a_f_done) else {
                continue;
            };
            let a_key = shape_map_key(&a_f.0);
            let mut a_lf_images: Vec<TopoShape> = self
                .my_of_images
                .seek(&a_f.0)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            let mut a_lf_inv: Vec<TopoShape> = Vec::new();
            // `cxx:1638-1642`.
            let b_artificial_case = self.my_art_invalid_faces.contains_key(&a_key);
            if b_artificial_case {
                a_lf_inv = a_lf_images.clone();
            } else {
                // `cxx:1646-1650`: the neutral edges, empty when unbound.
                let an_empty_map = ShapeIndexedMap::new();
                let p_mne = a_dmfmne.get(&a_key).unwrap_or(&an_empty_map);

                // `cxx:1653-1654`: the hole wires of the original face.
                let Some(a_f_or) = self
                    .my_faces_origins
                    .get(&a_key)
                    .cloned()
                    .and_then(Face::wrap)
                else {
                    // `myFacesOrigins->Find` raises on an absent key.
                    continue;
                };
                let mut a_mf_holes = ShapeIndexedMap::new();
                self.find_faces_inside_hole_wires(
                    &a_f_or,
                    &a_f,
                    &a_lf_images,
                    &a_dme_or_le_im,
                    &an_ef_map,
                    &mut a_mf_holes,
                );

                // `cxx:1656-1660`: the hole faces.
                for a_f_hole in a_mf_holes.iter() {
                    a_f_holes_children.push(a_f_hole.clone());
                }

                // `cxx:1663-1672`: the invalid faces of this offset face.
                self.find_invalid_faces(
                    &mut a_lf_images,
                    &a_dmfmve,
                    &a_dmfmie,
                    p_mne,
                    &a_me_edge_invalid_by_vertex,
                    &a_me_edge_valid_by_vertex,
                    &a_mf_holes,
                    &mut a_mf_inv_in_hole,
                    &mut a_lf_inv,
                    &mut an_inverted_faces,
                );
            }

            // `cxx:1675-1690`.
            if !a_lf_inv.is_empty() {
                if let Some(&a_nb) = self.my_already_inv_faces.get(&a_key) {
                    if a_nb > 2 {
                        // `cxx:1681-1684`: everything is invalid - drop all.
                        if a_lf_inv.len() == a_lf_images.len() && !b_artificial_case {
                            a_lf_images.clear();
                        }
                        a_lf_inv.clear();
                    }
                }
                self.my_invalid_faces.bind_list(&a_f.0, &a_lf_inv);
            }

            // `aLFImages` is a reference into `myOFImages` (`cxx:1635`).
            self.my_of_images.bind_list(&a_f.0, &a_lf_images);
        }

        // `cxx:1693-1697`.
        if self.my_invalid_faces.extent() == 0 {
            self.my_invalid_edges.clear();
            return;
        }

        // `cxx:1716-1723`: drop the invalid splits built from inverted edges.
        let mut a_meremoved = ShapeIndexedMap::new();
        self.remove_invalid_splits_by_inverted_edges(&mut a_meremoved);
        if self.my_invalid_faces.extent() == 0 {
            self.my_invalid_edges.clear();
            return;
        }

        // `cxx:1726`.
        self.remove_invalid_splits_from_valid(&a_dmfmvie);

        // `cxx:1729-1735`: drop the faces inside the volume.
        let a_f_holes = make_compound_from(a_f_holes_children.iter());
        let mut a_me_inside = std::mem::take(&mut self.my_inside_edges);
        self.remove_inside_faces(
            &an_inverted_faces,
            &a_mf_to_check_int,
            &a_mf_inv_in_hole,
            &a_f_holes,
            &mut a_meremoved,
            &mut a_me_inside,
        );
        self.my_inside_edges = a_me_inside;

        // `cxx:1738-1751`: the compound of the valid splits.
        let mut a_cf_im_children: Vec<TopoShape> = Vec::new();
        let a_nb = self.my_of_images.extent();
        for i in 1..=a_nb {
            let Some(a_f) = self.my_of_images.keys().get(i - 1).cloned() else {
                continue;
            };
            let a_lf_im: Vec<TopoShape> = self
                .my_of_images
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for a_f_im in a_lf_im {
                a_cf_im_children.push(a_f_im);
            }
        }
        let a_cf_im = make_compound_from(a_cf_im_children.iter());

        // `cxx:1753-1755`.
        let a_dmef = map_shapes_and_ancestors(&a_cf_im, ShapeType::Edge, ShapeType::Face);

        // `cxx:1758`: filter the image / origin maps.
        self.filter_edges_images(&a_cf_im);

        // `cxx:1761`: filter the invalid faces; the ternary picks the inside
        // edges when something was removed.
        if a_meremoved.extent() != 0 {
            let a_me_inside = std::mem::take(&mut self.my_inside_edges);
            self.filter_invalid_faces(&a_dmef, &a_me_inside);
            self.my_inside_edges = a_me_inside;
        } else {
            self.filter_invalid_faces(&a_dmef, &a_meremoved);
        }

        // `cxx:1762-1767`.
        if self.my_invalid_faces.extent() == 0 {
            self.my_invalid_edges.clear();
            return;
        }

        // `cxx:1787-1791`: filter the invalid edges, same ternary.
        let mut a_me_use_in_rebuild = ShapeIndexedMap::new();
        if a_meremoved.extent() != 0 {
            let a_me_inside = std::mem::take(&mut self.my_inside_edges);
            self.filter_invalid_edges(
                &a_dmfmie,
                &a_meremoved,
                &a_me_inside,
                &mut a_me_use_in_rebuild,
            );
            self.my_inside_edges = a_me_inside;
        } else {
            self.filter_invalid_edges(
                &a_dmfmie,
                &a_meremoved,
                &a_meremoved,
                &mut a_me_use_in_rebuild,
            );
        }

        // `cxx:1794`: re-check the edges created from a vertex.
        self.check_edges_created_by_vertex();

        // `cxx:1808-1818`.
        self.my_last_inv_edges.clear();
        let a_nb = self.my_invalid_edges.extent();
        for i in 1..=a_nb {
            let Some(a_e) = self.my_invalid_edges.value(i).cloned() else {
                continue;
            };
            self.my_last_inv_edges.add(&a_e);
            if !a_me_use_in_rebuild.contains(&a_e) {
                self.my_edges_to_avoid.add(&a_e);
            }
        }

        // `cxx:1820-1832`: bump the rebuild counter of every invalid face.
        let a_nb = self.my_invalid_faces.extent();
        for i in 1..=a_nb {
            let Some(a_f) = self.my_invalid_faces.keys().get(i - 1).cloned() else {
                continue;
            };
            let a_key = shape_map_key(&a_f);
            match self.my_already_inv_faces.get_mut(&a_key) {
                Some(a_nb) => *a_nb += 1,
                None => {
                    self.my_already_inv_faces.insert(a_key, 1);
                }
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::BuildSplitsOfExtendedFaces(const
    /// Message_ProgressRange& theRange)` (`cxx:1088-1169`).
    ///
    /// The retry loop around [`Self::build_splits_of_faces`]: every attempt
    /// starts from clean invalidity maps, splits the faces again, and - when
    /// invalid faces are left - rebuilds them through
    /// [`Self::intersect_faces_driver`]. The loop stops as soon as the splits
    /// come out clean or nothing is left to rebuild, and is capped at ten
    /// attempts because each invalid face can only be rebuilt three times
    /// (`cxx:1113-1115`).
    ///
    /// `myImage` is passed in as `the_image` (see the module note of
    /// [`super::BRepOffsetBuildOffsetFaces`]); the history is filled at the end
    /// (`cxx:1168`), as the reference does. `Message_ProgressScope`
    /// (`cxx:1097`, `:1142`) and the `aPart` halving it feeds
    /// (`cxx:1102`, `:1120-1121`) are the only use of `theRange`; the port has
    /// no progress reporting, so the two `if (!aPS.More()) return;` arms
    /// (`cxx:1106`, `:1123`) are absent. Every other branch is kept.
    pub fn build_splits_of_extended_faces(&mut self, the_image: &mut BRepAlgoImage) {
        // `cxx:1091-1095`.
        if !has_data(&self.my_faces)
            || !has_data(&self.my_edges_origins)
            || !has_data(&self.my_faces_origins)
            || !has_data(&self.my_e_trim_e_inf)
        {
            return;
        }

        // `cxx:1105`: fuse the trimmed offset edges to avoid self-intersecting
        // splits.
        self.intersect_trimmed_edges();

        // `cxx:1111`.
        let mut a_verts_to_avoid = ShapeIndexedMap::new();

        // `cxx:1115`.
        const A_NB_MAX_ATTEMPTS: usize = 10;

        // `cxx:1121-1161`.
        for _i_count in 1..=A_NB_MAX_ATTEMPTS {
            // `cxx:1128-1137`: every attempt starts from clean invalidity data.
            self.my_invalid_faces = crate::brep_offset::Dmvv::new();
            self.my_art_invalid_faces.clear();
            self.my_art_invalid_face_order.clear();
            self.my_invalid_edges.clear();
            self.my_inverted_edges.clear();
            self.my_ss_interfs.clear();
            self.my_ss_interfs_art.clear();
            self.my_intersection_pairs.clear();
            self.my_solids = None;
            self.my_faces_to_rebuild = crate::brep_offset::Dmvv::new();
            self.my_f_self_reb_avoid.clear();

            // `cxx:1145`.
            self.build_splits_of_faces();
            // `cxx:1146-1149`.
            if self.my_invalid_faces.extent() == 0 {
                break;
            }

            // `cxx:1152-1156`.
            self.find_faces_to_rebuild();
            if self.my_faces_to_rebuild.extent() == 0 {
                break;
            }

            // `cxx:1159-1160`.
            self.my_modified_edges.clear();
            self.intersect_faces_driver(&mut a_verts_to_avoid);
        }

        // `cxx:1165`: fill the gaps the splits still leave.
        self.fill_gaps();

        // `cxx:1168`.
        self.fill_history(the_image);
    }
}
