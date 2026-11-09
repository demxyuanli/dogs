//! Edge side of `BRepOffset_BuildOffsetFaces`: the edge selection
//! (`GetEdges`, `_1.cxx:1839-1953`) and the two edge passes
//! (`IntersectTrimmedEdges` `_1.cxx:1207-1292`,
//! `UpdateIntersectedEdges` `_1.cxx:1172-1203`).

use crate::abs::ShapeType;
use crate::bop_builder2::BopBuilder;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::builder::TopoBuilder;
use crate::shape::{Compound, Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, shape_map_key};

use super::statics::{append_to_list, process_micro_edge, AddToContainer};
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `bool BRepOffset_BuildOffsetFaces::GetEdges(const TopoDS_Face&,
    /// TopoDS_Shape&, IndexedMap* theInvMap = nullptr)` (`_1.cxx:1839-1953`).
    ///
    /// Builds the compound of offset edges that split `the_face` and answers
    /// `bFound && bUpdate`:
    ///
    /// * `bFound` - at least one such edge is not a boundary edge of the face
    ///   after the `myOEImages` substitution (`_1.cxx:1914-1917`, `:1947`);
    /// * `bUpdate` - at least one of the face's descendants or of their images
    ///   is in `myModifiedEdges` (`_1.cxx:1878-1881`, `:1919-1922`,
    ///   `:1928-1931`).
    ///
    /// `the_inv_map` receives the edges skipped because they are in
    /// `myEdgesToAvoid` (`_1.cxx:1896-1903`, `:1928-1936`).
    pub fn get_edges(
        &self,
        the_face: &crate::shape::Face,
        mut the_inv_map: Option<&mut ShapeIndexedMap>,
    ) -> (TopoShape, bool) {
        // `_1.cxx:1844-1861`: the boundary edges, with the splits from
        // `myOEImages` substituted for the originals.
        let mut a_m_bounds = ShapeIndexedMap::new();
        for a_e in explorer_shapes(&the_face.0, ShapeType::Edge) {
            match self.my_oe_images.seek(&a_e) {
                Some(p_le_im) => {
                    for a_e_im in p_le_im {
                        a_m_bounds.add(a_e_im);
                    }
                }
                None => {
                    a_m_bounds.add(&a_e);
                }
            }
        }

        // `_1.cxx:1863-1869`.
        let a_bb = TopoBuilder::new();
        let mut b_found = false;
        let mut b_update = false;
        // `TopoDS_Compound anEdges; aBB.MakeCompound(anEdges);`
        let mut an_edges = Compound::new();
        let mut a_m_fence = ShapeIndexedMap::new();

        // `_1.cxx:1871-1950`: the edges by which the offset face is split.
        let a_le = self.my_as_des.descendant(&the_face.0).to_vec();
        for a_e_shape in &a_le {
            // `TopoDS::Edge(aItLE.Value())` (`_1.cxx:1876`).
            let Some(a_e) = Edge::wrap(a_e_shape.clone()) else {
                continue;
            };
            // `_1.cxx:1878-1881`.
            if !b_update {
                b_update = self.my_modified_edges.contains(&a_e.0);
            }

            match self.my_oe_images.seek(&a_e.0).map(|l| l.to_vec()) {
                // `_1.cxx:1883-1925`.
                Some(p_le_im) => {
                    for a_e_im_shape in &p_le_im {
                        let Some(a_e_im) = Edge::wrap(a_e_im_shape.clone()) else {
                            continue;
                        };
                        // `_1.cxx:1890-1893`.
                        if !a_m_fence.add_to_container(&a_e_im.0) {
                            continue;
                        }
                        // `_1.cxx:1895-1907`.
                        if self.my_edges_to_avoid.contains(&a_e_im.0) {
                            if let Some(inv) = the_inv_map.as_deref_mut() {
                                inv.add(&a_e_im.0);
                            }
                            if !b_update {
                                b_update = self.my_last_inv_edges.contains(&a_e_im.0);
                            }
                            continue;
                        }
                        // `_1.cxx:1908-1912`.
                        if process_micro_edge(&a_e_im) {
                            continue;
                        }
                        // `_1.cxx:1913-1923`.
                        a_bb.add(&mut an_edges.0, &a_e_im.0);
                        if !b_found {
                            b_found = !a_m_bounds.contains(&a_e_im.0);
                        }
                        if !b_update {
                            b_update = self.my_modified_edges.contains(&a_e_im.0);
                        }
                    }
                }
                // `_1.cxx:1926-1949`: no images - the reference repeats the
                // same block for the edge itself (and skips the fence).
                None => {
                    if self.my_edges_to_avoid.contains(&a_e.0) {
                        if let Some(inv) = the_inv_map.as_deref_mut() {
                            inv.add(&a_e.0);
                        }
                        if !b_update {
                            b_update = self.my_last_inv_edges.contains(&a_e.0);
                        }
                        continue;
                    }
                    if process_micro_edge(&a_e) {
                        continue;
                    }
                    a_bb.add(&mut an_edges.0, &a_e.0);
                    if !b_found {
                        b_found = !a_m_bounds.contains(&a_e.0);
                    }
                }
            }
        }

        // `_1.cxx:1951-1953`.
        (an_edges.0, b_found && b_update)
    }

    /// `void BRepOffset_BuildOffsetFaces::UpdateIntersectedEdges(const
    /// NCollection_List<TopoDS_Shape>&, BOPAlgo_Builder&)`
    /// (`_1.cxx:1172-1203`).
    ///
    /// Binds every image of an intersected edge to the infinite edge the source
    /// edge maps to in `myETrimEInf`. `the_gf.Modified(aS)` maps onto
    /// [`crate::bop_hist::BopHistory::modified`].
    pub fn update_intersected_edges(&mut self, the_la: &[TopoShape], the_gf: &BopBuilder) {
        for a_s in the_la {
            // `_1.cxx:1182-1186`.
            let Some(p_e_inf) = self.my_e_trim_e_inf.get(&shape_map_key(a_s)).cloned() else {
                continue;
            };
            // `_1.cxx:1188-1192`.
            let a_ls_im: Vec<TopoShape> = match the_gf.history().modified(a_s) {
                Some(a_list) if !a_list.is_empty() => a_list.to_vec(),
                _ => continue,
            };
            // `_1.cxx:1194-1201`.
            for a_e_im in &a_ls_im {
                if !self.my_e_trim_e_inf.contains_key(&shape_map_key(a_e_im)) {
                    self.my_e_trim_e_inf.insert(shape_map_key(a_e_im), p_e_inf.clone());
                }
            }
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::IntersectTrimmedEdges(const
    /// Message_ProgressRange&)` (`_1.cxx:1207-1292`).
    ///
    /// Collects the descendants of the offset faces that are not micro edges
    /// and not already in `myModifiedEdges`, fuses them with a general fuse,
    /// then records the splits in `myOEImages` / `myOEOrigins` and refreshes
    /// `myEdgesOrigins` / `myETrimEInf`.
    ///
    /// `BOPAlgo_Builder aGFE; aGFE.SetArguments(aLS); aGFE.Perform()`
    /// (`_1.cxx:1243-1246`) is exactly `BopBuilder::set_arguments` +
    /// `perform` - the general-fuse algorithm without the object/tool filter
    /// (`BOPAlgo_Builder::Perform`). `aGFE.HasErrors()` maps onto
    /// [`BopBuilder::has_errors`].
    pub fn intersect_trimmed_edges(&mut self) {
        // `_1.cxx:1209-1237`: the edges to intersect.
        let mut a_ls: Vec<TopoShape> = Vec::new();
        for a_f_shape in self.my_faces.clone() {
            let Some(a_f) = crate::shape::Face::wrap(a_f_shape) else {
                continue;
            };
            // `_1.cxx:1222`: `myAsDes->Descendant(aF)`.
            let a_le = self.my_as_des.descendant(&a_f.0).to_vec();
            for a_e_shape in &a_le {
                let Some(a_e) = Edge::wrap(a_e_shape.clone()) else {
                    continue;
                };
                // `_1.cxx:1229-1232`.
                if process_micro_edge(&a_e) {
                    continue;
                }
                // `_1.cxx:1233-1236`: `myModifiedEdges.Add(aE)` answers whether
                // the edge was new.
                if self.my_modified_edges.add_to_container(&a_e.0) {
                    a_ls.push(a_e.0.clone());
                }
            }
        }

        // `_1.cxx:1240-1243`: nothing to intersect.
        if a_ls.len() < 2 {
            return;
        }

        // `_1.cxx:1244-1254`: the fuse of the collected edges.
        let mut a_gfe = BopBuilder::new();
        a_gfe.set_arguments(&a_ls);
        if a_gfe.perform().is_err() || a_gfe.has_errors() {
            return;
        }

        // `_1.cxx:1256-1288`: the images and origins of every split edge.
        let mut a_la: Vec<TopoShape> = Vec::new();
        for a_e in &a_ls {
            let a_le_im: Vec<TopoShape> = match a_gfe.history().modified(a_e) {
                Some(a_list) if !a_list.is_empty() => a_list.to_vec(),
                _ => continue,
            };
            // `_1.cxx:1271-1273`.
            a_la.push(a_e.clone());
            self.my_oe_images.bind_list(a_e, &a_le_im);
            // `_1.cxx:1276-1285`: `myOEOrigins` with `AppendToList` de-duplication.
            for a_e_im in &a_le_im {
                let key = shape_map_key(a_e_im);
                match self.my_oe_origins.get_mut(&key) {
                    Some(p_le_or) => append_to_list(p_le_or, a_e),
                    None => {
                        self.my_oe_origins.insert(key, vec![a_e.clone()]);
                    }
                }
            }
        }

        // `_1.cxx:1290-1291`.
        super::statics::update_origins(&a_la, &mut self.my_edges_origins, &a_gfe);
        self.update_intersected_edges(&a_la, &a_gfe);
    }
}
