//! Validity tests of the offset face splits used by
//! `BRepOffset_BuildOffsetFaces::BuildSplitsOfFaces`
//! (`BRepOffset_MakeOffset_1.cxx:1342-1838`).
//!
//! * `CheckIfArtificial` (`cxx:1960-2056`) - tells the "artificial" case: the
//!   splits built without the invalid edges already use every boundary edge,
//!   so they have to be rebuilt with the invalid edges added.
//! * `CheckInverted` (`cxx:3316-3516`) - compares the first->last vertex
//!   direction of an offset edge with the same direction of the original edge
//!   and records the edges whose vertices swapped places in `myInvertedEdges`.

use std::collections::HashMap;

use occt_core::gp::GpVec;

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::Dmvv;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_same, shape_map_key, vertex_position, AncestorMap,
    ShapeMapKey,
};

use super::statics::{append_to_list, find_shape, AddToContainer};
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `bool BRepOffset_BuildOffsetFaces::CheckIfArtificial(const
    /// TopoDS_Shape& theF, const List& theLFImages, const TopoDS_Shape& theCE,
    /// const IndexedMap& theMapEInv, Map& theMENInv)` (`cxx:1960-2056`).
    ///
    /// `theMENInv` receives the neutral edges: the edges of `theCE` that share
    /// a vertex with an invalid edge and are not used by the splits. Answers
    /// `true` when the walk over `myAsDes->Descendant(theF)` stopped early -
    /// the reference's trailing `return aItLE.More()`.
    pub fn check_if_artificial(
        &self,
        the_f: &TopoShape,
        the_lf_images: &[TopoShape],
        the_ce: &TopoShape,
        the_map_e_inv: &ShapeIndexedMap,
        the_men_inv: &mut ShapeIndexedMap,
    ) -> bool {
        // `cxx:1967-1976`: every edge and vertex used by the splits.
        let mut a_me_used = ShapeIndexedMap::new();
        for a_f_im in the_lf_images {
            for a_s in explorer_shapes(a_f_im, ShapeType::Edge) {
                a_me_used.add(&a_s);
            }
            for a_s in explorer_shapes(a_f_im, ShapeType::Vertex) {
                a_me_used.add(&a_s);
            }
        }

        // `cxx:1978-1979`: the vertices of the split edges -> those edges.
        let a_mve = AncestorMap::map_shapes_and_ancestors(the_ce, ShapeType::Vertex, ShapeType::Edge);

        // `cxx:1981-2003`: an edge of `theCE` touching an invalid edge that no
        // split uses is neutral.
        for i in 1..=the_map_e_inv.extent() {
            let Some(a_e_inv) = the_map_e_inv.value(i) else {
                continue;
            };
            let Some(a_e_inv) = Edge::wrap(a_e_inv.clone()) else {
                continue;
            };
            for a_v in explorer_shapes(&a_e_inv.0, ShapeType::Vertex) {
                let Some(p_le_n_inv) = a_mve.find_from_key(&a_v) else {
                    continue;
                };
                for a_en_inv in p_le_n_inv {
                    if !a_me_used.contains(a_en_inv) {
                        the_men_inv.add(a_en_inv);
                    }
                }
            }
        }

        // `cxx:2005-2008`.
        if the_men_inv.extent() == 0 {
            return false;
        }

        // `cxx:2010-2012`.
        let mut a_me_found = ShapeIndexedMap::new();
        for a_e in explorer_shapes(the_ce, ShapeType::Edge) {
            a_me_found.add(&a_e);
        }

        // `cxx:2014-2052`: the first descendant edge (or all of its images)
        // that is a split edge, not neutral and not used by the splits makes
        // the case artificial.
        let a_le = self.my_as_des.descendant(the_f).to_vec();
        for a_e_shape in &a_le {
            let Some(a_e) = Edge::wrap(a_e_shape.clone()) else {
                continue;
            };
            match self.my_oe_images.seek(&a_e.0).map(|l| l.to_vec()) {
                // `cxx:2018-2041`.
                Some(p_le_im) => {
                    let mut b_checked = false;
                    let mut b_used = false;
                    for a_e_im_shape in &p_le_im {
                        let Some(a_e_im) = Edge::wrap(a_e_im_shape.clone()) else {
                            continue;
                        };
                        // `cxx:2027-2030`.
                        if !a_me_found.contains(&a_e_im.0) || the_men_inv.contains(&a_e_im.0) {
                            continue;
                        }
                        b_checked = true;
                        // `cxx:2033-2036`.
                        if a_me_used.contains(&a_e_im.0) {
                            b_used = true;
                            break;
                        }
                    }
                    // `cxx:2039-2041`: every image was checked and none is used.
                    if b_checked && !b_used {
                        return true;
                    }
                }
                // `cxx:2043-2049`.
                None => {
                    if a_me_found.contains(&a_e.0)
                        && !the_men_inv.contains(&a_e.0)
                        && !a_me_used.contains(&a_e.0)
                    {
                        return true;
                    }
                }
            }
        }

        // `cxx:2052`: the walk ran to the end.
        false
    }

    /// `bool BRepOffset_BuildOffsetFaces::CheckInverted(const TopoDS_Edge&
    /// theEIm, const TopoDS_Face& theFOr, const IndexedDataMap& theDMVE,
    /// const IndexedMap& theMEdges)` (`cxx:3316-3516`).
    ///
    /// `theDMVE` is the vertex -> split edge map `FindInvalidEdges` builds and
    /// `theMEdges` the set of all split edges. When the offset edge is in
    /// `myOEOrigins`, the compared vertices are the two loose ends of its
    /// images; otherwise they are its own vertices.
    pub fn check_inverted(
        &mut self,
        the_e_im: &Edge,
        the_f_or: &Face,
        the_dmve: &Dmvv,
        the_m_edges: &ShapeIndexedMap,
    ) -> bool {
        // `cxx:3329-3335`: `TopExp::Vertices(theEIm, aVI1, aVI2)`. OCCT
        // dereferences both vertices unconditionally; an edge missing one
        // cannot be compared, so the port answers `false`.
        let (a_vi1, a_vi2) = edge_vertices(the_e_im);
        let (Some(mut a_vi1), Some(mut a_vi2)) = (a_vi1, a_vi2) else {
            return false;
        };

        // `cxx:3337-3412`.
        let mut a_le_images: Vec<TopoShape> = Vec::new();
        match self
            .my_oe_origins
            .get(&shape_map_key(&the_e_im.0))
            .cloned()
        {
            Some(a_l_offset_or) => {
                // `cxx:3340-3361`: the images of the origins that are split
                // edges, once each.
                let mut an_images_edges: Vec<Edge> = Vec::new();
                let mut a_m_im_fence = ShapeIndexedMap::new();
                for a_e_offset_or in &a_l_offset_or {
                    let Some(a_l_images) =
                        self.my_oe_images.seek(a_e_offset_or).map(|l| l.to_vec())
                    else {
                        // `myOEImages.Find` asserts on an unbound key in OCCT.
                        continue;
                    };
                    for an_im_shape in &a_l_images {
                        let Some(an_im) = Edge::wrap(an_im_shape.clone()) else {
                            continue;
                        };
                        if the_m_edges.contains(&an_im.0)
                            && a_m_im_fence.add_to_container(&an_im.0)
                        {
                            an_images_edges.push(an_im.clone());
                            a_le_images.push(an_im.0.clone());
                        }
                    }
                }

                // `cxx:3364-3380`: the vertices belonging to a single image
                // edge are the loose ends of the images wire.
                let an_images = TopoBuilder::new().make_wire(&an_images_edges);
                let a_dm_im_ve = AncestorMap::map_shapes_and_ancestors(
                    &an_images.0,
                    ShapeType::Vertex,
                    ShapeType::Edge,
                );
                let mut a_lv_alone: Vec<TopoShape> = Vec::new();
                for i in 0..a_dm_im_ve.length() {
                    if a_dm_im_ve.ancestors(i).len() == 1 {
                        a_lv_alone.push(a_dm_im_ve.find_key(i).clone());
                    }
                }

                // `cxx:3381-3403`: the loose end closer to the first vertex of
                // the compared edge is that vertex.
                if a_lv_alone.len() > 1 {
                    let (Some(a_vw1), Some(a_vw2)) = (
                        Vertex::wrap(a_lv_alone[0].clone()),
                        Vertex::wrap(a_lv_alone[a_lv_alone.len() - 1].clone()),
                    ) else {
                        return false;
                    };
                    let a_p_i1 = vertex_position(&a_vi1);
                    let a_p_w1 = vertex_position(&a_vw1);
                    let a_p_w2 = vertex_position(&a_vw2);
                    if a_p_i1.square_distance(&a_p_w1) < a_p_i1.square_distance(&a_p_w2) {
                        a_vi1 = a_vw1;
                        a_vi2 = a_vw2;
                    } else {
                        a_vi1 = a_vw2;
                        a_vi2 = a_vw1;
                    }
                }
            }
            None => {
                // `cxx:3408-3411`.
                a_le_images.push(the_e_im.0.clone());
            }
        }

        // `cxx:3414-3415`: the split edges incident to those two vertices.
        let a_lie1: Vec<TopoShape> = the_dmve.seek(&a_vi1.0).unwrap_or(&[]).to_vec();
        let a_lie2: Vec<TopoShape> = the_dmve.seek(&a_vi2.0).unwrap_or(&[]).to_vec();

        // `cxx:3419-3451`: the original edges of both lists that belong to the
        // original face.
        let mut a_loe1: Vec<TopoShape> = Vec::new();
        let mut a_loe2: Vec<TopoShape> = Vec::new();
        for i in 0..2 {
            let a_lie = if i == 0 { &a_lie1 } else { &a_lie2 };
            let a_loe = if i == 0 { &mut a_loe1 } else { &mut a_loe2 };
            let mut a_m_fence = ShapeIndexedMap::new();
            for a_ei in a_lie {
                let Some(a_le_origins) = self.my_edges_origins.get(&shape_map_key(a_ei)).cloned()
                else {
                    continue;
                };
                for a_eo in &a_le_origins {
                    // `cxx:3439`.
                    if a_eo.shape_type() != ShapeType::Edge || !a_m_fence.add_to_container(a_eo) {
                        continue;
                    }
                    // `FindShape(aEO, theFOr, nullptr, aEOin)`.
                    let mut a_eo_in = TopoShape::new(ShapeType::Edge);
                    if find_shape(a_eo, &the_f_or.0, None, &mut a_eo_in) {
                        append_to_list(a_loe, a_eo);
                    }
                }
            }
        }

        // `cxx:3453-3456`.
        if a_loe1.len() < 2 || a_loe2.len() < 2 {
            return false;
        }

        // `cxx:3458-3486`: the vertex shared by the largest number of edges,
        // mapped over all edges of the list into one `MapShapesAndAncestors`.
        let mut a_vo: [Option<TopoShape>; 2] = [None, None];
        for i in 0..2 {
            let a_loe = if i == 0 { &a_loe1 } else { &a_loe2 };
            let mut a_dmve_loc: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
            let mut a_index: HashMap<ShapeMapKey, usize> = HashMap::new();
            for a_e in a_loe {
                let a_loc = AncestorMap::map_shapes_and_ancestors(
                    a_e,
                    ShapeType::Vertex,
                    ShapeType::Edge,
                );
                for j in 0..a_loc.length() {
                    let a_v = a_loc.find_key(j);
                    match a_index.get(&shape_map_key(a_v)) {
                        Some(&k) => a_dmve_loc[k].1.extend(a_loc.ancestors(j).iter().cloned()),
                        None => {
                            a_index.insert(shape_map_key(a_v), a_dmve_loc.len());
                            a_dmve_loc.push((a_v.clone(), a_loc.ancestors(j).to_vec()));
                        }
                    }
                }
            }

            // `cxx:3471-3481`.
            let mut a_nb_e_max = 0;
            let mut a_v_found: Option<TopoShape> = None;
            for (a_v, a_le) in &a_dmve_loc {
                let a_nb_e = a_le.len();
                if a_nb_e > 1 && a_nb_e > a_nb_e_max {
                    a_v_found = Some(a_v.clone());
                    a_nb_e_max = a_nb_e;
                }
            }
            // `cxx:3482-3485`.
            match a_v_found {
                Some(a_v) => a_vo[i] = Some(a_v),
                None => return false,
            }
        }

        // `cxx:3488-3491`.
        let (Some(a_vo1), Some(a_vo2)) = (a_vo[0].clone(), a_vo[1].clone()) else {
            return false;
        };
        if is_same(&a_vo1, &a_vo2) {
            return false;
        }

        // `cxx:3493-3516`: the directions of the two vertex pairs.
        let (Some(a_vo1), Some(a_vo2)) = (Vertex::wrap(a_vo1), Vertex::wrap(a_vo2)) else {
            return false;
        };
        let a_vi = GpVec::from_pnts(&vertex_position(&a_vi1), &vertex_position(&a_vi2));
        let a_vo = GpVec::from_pnts(&vertex_position(&a_vo1), &vertex_position(&a_vo2));
        let an_angle = a_vi.angle(&a_vo);
        let b_inverted = (an_angle - std::f64::consts::PI).abs() < 1.0e-4;
        if b_inverted {
            for a_e_invr in &a_le_images {
                self.my_inverted_edges.add(a_e_invr);
            }
        }
        b_inverted
    }
}
