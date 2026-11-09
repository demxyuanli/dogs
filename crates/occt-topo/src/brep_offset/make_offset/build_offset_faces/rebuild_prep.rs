//! The preparation half of the rebuild: `FindFacesToRebuild`
//! (`cxx:5758-5897`), the file-local `mapShapes` (`cxx:5904-5917`),
//! `PrepareFacesForIntersection` (`cxx:6571-6674`), `FindVerticesToAvoid`
//! (`cxx:6680-6782`) and `FindFacesForIntersection` (`cxx:6788-6969`).
//!
//! `FindFacesToRebuild` decides which offset faces must be re-split: a face
//! goes in `myFacesToRebuild` when one of its splits carries an edge or a
//! vertex of an invalid face (`cxx:5784-5806`, `:5857-5879`), when it is
//! connected to another face through `mySSInterfs` / `mySSInterfsArt`
//! (`cxx:5813-5821`), or when it has splits as well as invalid splits
//! (`cxx:5884-5890`, recorded in `myFSelfRebAvoid`).
//!
//! `PrepareFacesForIntersection` builds the four lookup tables the
//! intersection pass needs: split edge -> offset face and infinite edge ->
//! trimmed edges (`theDMSF`, `theMEInfETrim`, `cxx:6599-6619`), plus the
//! vertex -> face and vertex -> edge tables used only when vertices have to be
//! avoided (`theDMVEFull`, `cxx:6621-6645`), and invalid edge -> offset faces
//! of the non-artificial invalid splits (`theDMEFInv`, `cxx:6648-6672`).
//!
//! `FindVerticesToAvoid` collects the vertices an invalid edge ends at when all
//! the edges around them are invalid, or when at least two of them are
//! inverted (`cxx:6747-6780`); those vertices are then not used as split ends.
//!
//! `FindFacesForIntersection` walks the shapes of one invalidity
//! (`theME`) through `theDMSF` to find the faces sharing them
//! (`cxx:6803-6847`), then adds the faces recorded in `mySSInterfs` /
//! `mySSInterfsArt` for the same shapes that are connected to what was already
//! gathered (`cxx:6873-6926`) - and keeps the ones forming a block of a single
//! offset face (`cxx:6931-6968`).

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::algo_tools::AlgoTools;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::Dmvv;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::invalid_edges_filter::make_compound_from;
use super::statics::append_to_list;
use super::BRepOffsetBuildOffsetFaces;

/// `template <class Container> static void mapShapes(const Container&
/// theVecShapes, const TopAbs_ShapeEnum theType, Map& theMap)`
/// (`cxx:5904-5917`): every shape of the type under every entry of the
/// container into the map.
pub(super) fn map_shapes<'a, I>(the_vec_shapes: I, the_type: ShapeType, the_map: &mut ShapeIndexedMap)
where
    I: IntoIterator<Item = &'a TopoShape>,
{
    for a_shape in the_vec_shapes {
        for a_sub in explorer_shapes(a_shape, the_type) {
            the_map.add(&a_sub);
        }
    }
}

/// `theDMSF.ChangeSeek(aKey)` + `Bound` + `AppendToList` (`cxx:6606-6619`).
fn append_under(the_map: &mut HashMap<ShapeMapKey, Vec<TopoShape>>, key: &TopoShape, value: &TopoShape) {
    append_to_list(the_map.entry(shape_map_key(key)).or_default(), value);
}

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::FindFacesToRebuild()`
    /// (`cxx:5758-5897`).
    pub fn find_faces_to_rebuild(&mut self) {
        let a_nb = self.my_of_images.extent();
        // `cxx:5761-5764`.
        if a_nb == 0 {
            return;
        }

        // `cxx:5768-5771`: the edges and vertices of the invalid faces, the
        // faces to rebuild, and the vertices to avoid per face.
        let mut a_me_reb = ShapeIndexedMap::new();
        let mut a_mf_reb = ShapeIndexedMap::new();
        let mut a_dmflv: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        // `cxx:5774-5821`: the invalid faces.
        let a_nb = self.my_invalid_faces.extent();
        for i in 1..=a_nb {
            let Some(a_f) = self.my_invalid_faces.keys().get(i - 1).cloned() else {
                continue;
            };
            // `cxx:5778-5779`: cleared per face.
            let mut a_mfence = ShapeIndexedMap::new();
            let mut a_lv_avoid: Vec<TopoShape> = Vec::new();

            let a_lf_im: Vec<TopoShape> = self
                .my_invalid_faces
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();
            for a_f_im in &a_lf_im {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    a_me_reb.add(&a_e);
                    // `cxx:5790-5802`.
                    if self.my_invalid_edges.contains(&a_e) {
                        for a_v in explorer_shapes(&a_e, ShapeType::Vertex) {
                            if a_mfence.add(&a_v) == a_mfence.extent() {
                                append_to_list(&mut a_lv_avoid, &a_v);
                                a_me_reb.add(&a_v);
                            }
                        }
                    }
                }
            }

            // `cxx:5805-5808`.
            if !a_lv_avoid.is_empty() {
                a_dmflv.insert(shape_map_key(&a_f), a_lv_avoid);
            }

            // `cxx:5810-5821`.
            let p_lf = if !self.my_art_invalid_faces.contains_key(&shape_map_key(&a_f)) {
                self.my_ss_interfs.get(&shape_map_key(&a_f)).cloned()
            } else {
                self.my_ss_interfs_art.get(&shape_map_key(&a_f)).cloned()
            };
            if let Some(a_lf) = p_lf {
                for a_fe in &a_lf {
                    a_mf_reb.add(a_fe);
                }
            }
        }

        // `cxx:5824-5895`: one offset face at a time.
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

            // `cxx:5830-5841`.
            let mut a_mv_avoid = ShapeIndexedMap::new();
            if let Some(a_lv_avoid) = a_dmflv.get(&shape_map_key(&a_f)) {
                for a_v in a_lv_avoid.clone() {
                    a_mv_avoid.add(&a_v);
                }
            }

            // `cxx:5843-5845`.
            let mut b_rebuild = a_mf_reb.contains(&a_f);
            let mut a_le_valid: Vec<TopoShape> = Vec::new();
            let mut a_mfence = ShapeIndexedMap::new();

            // `cxx:5847-5880`.
            for a_f_im in &a_lf_im {
                for a_e_im in explorer_shapes(a_f_im, ShapeType::Edge) {
                    // `cxx:5854-5861`.
                    if !self.my_invalid_edges.contains(&a_e_im) {
                        if a_mfence.add(&a_e_im) == a_mfence.extent() {
                            append_to_list(&mut a_le_valid, &a_e_im);
                        }
                    }
                    // `cxx:5863-5866`.
                    if !b_rebuild {
                        b_rebuild = a_me_reb.contains(&a_e_im);
                    }
                    // `cxx:5868-5879`.
                    if !b_rebuild {
                        for a_v in explorer_shapes(&a_e_im, ShapeType::Vertex) {
                            if !a_mv_avoid.contains(&a_v) {
                                b_rebuild = a_me_reb.contains(&a_v);
                            }
                            if b_rebuild {
                                break;
                            }
                        }
                    }
                }
            }

            // `cxx:5883-5890`: splits together with invalid splits.
            if !b_rebuild {
                b_rebuild = !a_lf_im.is_empty() && self.my_invalid_faces.contains(&a_f);
                if b_rebuild {
                    self.my_f_self_reb_avoid.add(&a_f);
                }
            }

            // `cxx:5892-5895`: `myFacesToRebuild.Add(aF, aLEValid)`.
            if b_rebuild {
                let mut a_cur = self
                    .my_faces_to_rebuild
                    .seek(&a_f)
                    .map(|l| l.to_vec())
                    .unwrap_or_default();
                a_cur.extend_from_slice(&a_le_valid);
                self.my_faces_to_rebuild.bind_list(&a_f, &a_cur);
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::PrepareFacesForIntersection(const
    /// bool theLookVertToAvoid, IndexedDataMap& theFLE, DataMap& theMDone,
    /// DataMap& theDMSF, DataMap& theMEInfETrim, DataMap& theDMVEFull,
    /// IndexedDataMap& theDMEFInv)` (`cxx:6571-6674`).
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_faces_for_intersection(
        &self,
        the_look_vert_to_avoid: bool,
        the_fle: &mut Dmvv,
        the_m_done: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_dmsf: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_me_inf_e_trim: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_dmve_full: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_dmef_inv: &mut Dmvv,
    ) {
        let a_nb = self.my_faces_to_rebuild.extent();
        for i in 1..=a_nb {
            let Some(a_f) = self.my_faces_to_rebuild.keys().get(i - 1).cloned() else {
                continue;
            };
            // `cxx:6592-6594`: bind an empty edge list for the face.
            the_fle.bind_list(&a_f, &[]);
            the_m_done.insert(shape_map_key(&a_f), Vec::new());

            // `cxx:6596`.
            let a_lf_im: Vec<TopoShape> = self
                .my_of_images
                .seek(&a_f)
                .map(|l| l.to_vec())
                .unwrap_or_default();

            for a_f_im in &a_lf_im {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    // `cxx:6604-6610`: save the connection to the untrimmed
                    // face.
                    append_under(the_dmsf, &a_e, &a_f);

                    // `cxx:6612-6619`: save the connection to the untrimmed
                    // edge. `myETrimEInf->Find` raises on an unbound edge in
                    // the reference; the port skips it.
                    let Some(a_e_inf) = self.my_e_trim_e_inf.get(&shape_map_key(&a_e)).cloned()
                    else {
                        continue;
                    };
                    append_under(the_me_inf_e_trim, &a_e_inf, &a_e);

                    // `cxx:6621-6644`.
                    for a_v in explorer_shapes(&a_e, ShapeType::Vertex) {
                        // Save the connection to the face.
                        append_under(the_dmsf, &a_v, &a_f);
                        if the_look_vert_to_avoid {
                            // Save the connection to the edges.
                            append_under(the_dmve_full, &a_v, &a_e);
                        }
                    }
                }
            }

            // `cxx:6647-6672`: the edges of the invalid splits, stacked under
            // their edge with the offset face.
            if the_look_vert_to_avoid {
                let p_lf_inv = self.my_invalid_faces.seek(&a_f).map(|l| l.to_vec());
                if p_lf_inv.is_none()
                    || self.my_art_invalid_faces.contains_key(&shape_map_key(&a_f))
                {
                    continue;
                }
                for a_f_inv in &p_lf_inv.unwrap_or_default() {
                    for a_e in explorer_shapes(a_f_inv, ShapeType::Edge) {
                        if !the_dmef_inv.contains(&a_e) {
                            the_dmef_inv.bind_list(&a_e, &[]);
                        }
                        the_dmef_inv.add(&a_e, &a_f);
                    }
                }
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::FindVerticesToAvoid(const
    /// IndexedDataMap& theDMEFInv, const DataMap& theDMVEFull, Map&
    /// theMVRInv)` (`cxx:6680-6782`).
    pub fn find_vertices_to_avoid(
        &self,
        the_dmef_inv: &Dmvv,
        the_dmve_full: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_mvr_inv: &mut ShapeIndexedMap,
    ) {
        // `cxx:6688`.
        let mut a_mfence = ShapeIndexedMap::new();

        let a_nb = the_dmef_inv.extent();
        for i in 1..=a_nb {
            // `cxx:6692-6697`.
            let Some(a_e) = the_dmef_inv.keys().get(i - 1).cloned() else {
                continue;
            };
            let a_lf_inv: Vec<TopoShape> = the_dmef_inv.seek(&a_e).map(|l| l.to_vec()).unwrap_or_default();
            if a_lf_inv.len() == 1 {
                continue;
            }

            // `cxx:6699-6707`.
            if !self.my_invalid_edges.contains(&a_e) || self.my_valid_edges.contains(&a_e) {
                continue;
            }
            if a_mfence.add(&a_e) != a_mfence.extent() {
                continue;
            }

            // `cxx:6709-6733`: the edges around the vertices of the edge, taken
            // from the images of its origins when they are known.
            let mut a_mve_edges = Dmvv::new();
            match self.my_oe_origins.get(&shape_map_key(&a_e)).cloned() {
                Some(p_le_or) => {
                    for a_e_or in &p_le_or {
                        // `myOEImages.Find` raises on an unbound edge in the
                        // reference; the port skips it.
                        let Some(a_le_im) = self.my_oe_images.seek(a_e_or).map(|l| l.to_vec())
                        else {
                            continue;
                        };
                        for a_e_im in &a_le_im {
                            a_mfence.add(a_e_im);
                            map_vertices_to_edges(a_e_im, &mut a_mve_edges);
                        }
                    }
                }
                None => {
                    map_vertices_to_edges(&a_e, &mut a_mve_edges);
                }
            }

            // `cxx:6735-6781`.
            let a_nb_v = a_mve_edges.extent();
            for j in 1..=a_nb_v {
                let Some(a_v) = a_mve_edges.keys().get(j - 1).cloned() else {
                    continue;
                };
                let a_le: Vec<TopoShape> = a_mve_edges.seek(&a_v).map(|l| l.to_vec()).unwrap_or_default();
                // `cxx:6737-6740`: only the vertices with a single edge.
                if a_le.len() != 1 {
                    continue;
                }
                // `cxx:6742-6746`.
                if a_mfence.add(&a_v) != a_mfence.extent() {
                    continue;
                }
                // `cxx:6747-6754`.
                let Some(p_le) = the_dmve_full.get(&shape_map_key(&a_v)) else {
                    // Isolated vertex.
                    the_mvr_inv.add(&a_v);
                    continue;
                };

                // `cxx:6758-6781`.
                let mut i_nb_e_inverted = 0usize;
                let mut b_all_edges_inv = true;
                for a_ev in p_le {
                    if self.my_inverted_edges.contains(a_ev) {
                        i_nb_e_inverted += 1;
                    }
                    if b_all_edges_inv {
                        b_all_edges_inv = self.my_invalid_edges.contains(a_ev);
                    }
                }
                if i_nb_e_inverted > 1 || b_all_edges_inv {
                    the_mvr_inv.add(&a_v);
                }
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::FindFacesForIntersection(const
    /// TopoDS_Shape& theFInv, const IndexedMap& theME, const DataMap&
    /// theDMSF, const Map& theMVInvAll, const bool theArtCase, IndexedMap&
    /// theMFAvoid, IndexedMap& theMFInt, IndexedMap& theMFIntExt, List&
    /// theLFImInt)` (`cxx:6788-6969`).
    #[allow(clippy::too_many_arguments)]
    pub fn find_faces_for_intersection(
        &self,
        the_f_inv: &TopoShape,
        the_me: &ShapeIndexedMap,
        the_dmsf: &HashMap<ShapeMapKey, Vec<TopoShape>>,
        the_mv_inv_all: &ShapeIndexedMap,
        the_art_case: bool,
        the_mf_avoid: &mut ShapeIndexedMap,
        the_mf_int: &mut ShapeIndexedMap,
        the_mf_int_ext: &mut ShapeIndexedMap,
        the_lf_im_int: &mut Vec<TopoShape>,
    ) {
        let a_nb_e = the_me.extent();
        let mut a_m_shapes = ShapeIndexedMap::new();

        // `cxx:6803-6847`: the faces sharing the shapes of the invalidity.
        for i in 1..=a_nb_e {
            let Some(a_s) = the_me.value(i).cloned() else {
                continue;
            };
            let Some(a_lf) = the_dmsf.get(&shape_map_key(&a_s)).cloned() else {
                continue;
            };

            // `cxx:6812-6814`: in the artificial case only the vertices that
            // are not on an invalid edge are avoided, and only artificial
            // invalid faces are marked.
            let b_avoid = the_art_case
                && a_s.shape_type() == ShapeType::Vertex
                && !the_mv_inv_all.contains(&a_s);

            for a_f in &a_lf {
                // `cxx:6821-6824`.
                if the_mf_int.contains(a_f) {
                    continue;
                }
                // `cxx:6827-6830`.
                if b_avoid && self.my_art_invalid_faces.contains_key(&shape_map_key(a_f)) {
                    the_mf_avoid.add(a_f);
                }
                // `cxx:6832`.
                the_mf_int.add(a_f);
                // `cxx:6834`.
                let b_use = !is_same(a_f, the_f_inv);

                // `cxx:6836-6846`.
                let a_lf_im: Vec<TopoShape> = self
                    .my_of_images
                    .seek(a_f)
                    .map(|l| l.to_vec())
                    .unwrap_or_default();
                for a_f_im in &a_lf_im {
                    append_to_list(the_lf_im_int, a_f_im);
                    if b_use {
                        map_shapes_into_edges(a_f_im, &mut a_m_shapes);
                    }
                }
            }
        }

        // `cxx:6849-6852`.
        let a_ss_interfs_map = if the_art_case {
            &self.my_ss_interfs_art
        } else {
            &self.my_ss_interfs
        };
        let Some(p_lf_inv) = a_ss_interfs_map.get(&shape_map_key(the_f_inv)).cloned() else {
            return;
        };

        // `cxx:6856-6862`.
        let mut a_mf = ShapeIndexedMap::new();
        for a_f in &p_lf_inv {
            a_mf.add(a_f);
        }

        // `cxx:6868-6871`: the splits to add and their original faces.
        let mut a_mf_to_add = ShapeIndexedMap::new();
        let mut a_dmf_or: HashMap<ShapeMapKey, TopoShape> = HashMap::new();
        let mut a_cf_splits: Vec<TopoShape> = Vec::new();

        // `cxx:6873-6926`.
        for i in 1..=a_nb_e {
            let Some(a_s) = the_me.value(i).cloned() else {
                continue;
            };
            let Some(p_lf) = a_ss_interfs_map.get(&shape_map_key(&a_s)).cloned() else {
                continue;
            };

            for a_f in &p_lf {
                // `cxx:6886-6889`.
                if the_mf_int.contains(a_f) || a_mf_to_add.contains(a_f) || !a_mf.contains(a_f) {
                    continue;
                }

                let a_lf_im: Vec<TopoShape> = self
                    .my_of_images
                    .seek(a_f)
                    .map(|l| l.to_vec())
                    .unwrap_or_default();

                // `cxx:6891-6915`: the face has to touch what was already
                // gathered.
                if !the_art_case {
                    let mut b_connected = false;
                    for a_f_im in &a_lf_im {
                        for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                            if a_m_shapes.contains(&a_e) {
                                b_connected = true;
                                break;
                            }
                        }
                        if b_connected {
                            break;
                        }
                    }
                    if !b_connected {
                        continue;
                    }
                }

                // `cxx:6918-6926`.
                a_mf_to_add.add(a_f);
                for a_f_im in &a_lf_im {
                    a_dmf_or.insert(shape_map_key(a_f_im), a_f.clone());
                    a_cf_splits.push(a_f_im.clone());
                }
            }
        }

        // `cxx:6928-6931`.
        if a_mf_to_add.extent() == 0 {
            return;
        }

        // `cxx:6933-6938`.
        let a_cf = make_compound_from(a_cf_splits.iter());
        let a_lcb = AlgoTools::make_connexity_blocks_in(&a_cf, ShapeType::Edge, ShapeType::Face);
        if a_lcb.len() == 1 && a_mf_to_add.extent() > 1 {
            return;
        }

        // `cxx:6941-6968`.
        for a_cb in &a_lcb {
            a_mf_to_add.clear();
            for a_f_im in explorer_shapes(a_cb, ShapeType::Face) {
                if let Some(a_f) = a_dmf_or.get(&shape_map_key(&a_f_im)) {
                    a_mf_to_add.add(a_f);
                }
            }

            // `cxx:6952-6968`.
            if a_mf_to_add.extent() == 1 {
                let Some(a_f) = a_mf_to_add.value(1).cloned() else {
                    continue;
                };
                the_mf_int.add(&a_f);
                the_mf_int_ext.add(&a_f);

                let a_lf_im: Vec<TopoShape> = self
                    .my_of_images
                    .seek(&a_f)
                    .map(|l| l.to_vec())
                    .unwrap_or_default();
                for a_f_im in &a_lf_im {
                    append_to_list(the_lf_im_int, a_f_im);
                }
            }
        }
    }
}

/// `TopExp::MapShapes(aFIm, TopAbs_EDGE, aMShapes)` (`cxx:6843`), accumulated.
fn map_shapes_into_edges(the_shape: &TopoShape, the_map: &mut ShapeIndexedMap) {
    for a_e in explorer_shapes(the_shape, ShapeType::Edge) {
        the_map.add(&a_e);
    }
}

/// `TopExp::MapShapesAndAncestors(aShape, TopAbs_VERTEX, TopAbs_EDGE, theMap)`
/// (`cxx:6725`, `:6731`), accumulated into a [`Dmvv`] so the vertex keys stay
/// readable (`FindKey(j)` in the reference). Every vertex of `the_shape` is
/// paired with `the_shape` itself, which is the reference's ancestor.
pub(super) fn map_vertices_to_edges(the_shape: &TopoShape, the_map: &mut Dmvv) {
    for a_v in explorer_shapes(the_shape, ShapeType::Vertex) {
        if !the_map.contains(&a_v) {
            the_map.bind_list(&a_v, &[]);
        }
        the_map.add(&a_v, the_shape);
    }
}
