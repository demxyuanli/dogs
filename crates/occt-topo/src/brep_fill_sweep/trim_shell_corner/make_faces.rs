//! `MakeFacesNonSec` and `MakeFacesSec` of `BRepFill_TrimShellCorner`.
//!
//! Moved here verbatim from the single-file port; both methods keep the
//! `BRepFill_TrimShellCorner.cxx` branch structure.

use super::*;
use super::helpers::*;
use super::choose_section::*;
use super::from_edges::*;

impl TrimShellCorner {
    /// `MakeFacesNonSec` (`cxx:443-783`).
    ///
    /// The bounds of the corner row do not share a section curve: the two face
    /// columns are joined at the common vertex found on their bound edges
    /// (`FindCommonVertex`) and/or between their `u`-iso edge pairs, then the
    /// substituted faces are recorded in `myHistMap`.
    pub(super) fn make_faces_non_sec(
        &mut self,
        the_index: usize,
        the_ds: &BopdsDS,
        the_face_index1: usize,
        the_face_index2: usize,
    ) -> Result<bool, String> {
        // `cxx:448-452`.
        let mut has_new_edge = false;
        let mut new_edge: Option<Edge> = None;

        let Some(bounds) = self.bounds.clone() else {
            return Ok(false);
        };
        let Some(u_edges) = self.u_edges.clone() else {
            return Ok(false);
        };
        let (Some(a_e1), Some(a_e2)) = (
            bounds.value(the_index, 1).cloned(),
            bounds.value(the_index, 2).cloned(),
        ) else {
            return Ok(false);
        };

        // `cxx:455-464`: the common vertex between the two bounds.
        let an_index1 = the_ds.index(&a_e1).unwrap_or(usize::MAX);
        let an_index2 = the_ds.index(&a_e2).unwrap_or(usize::MAX);
        let mut common_vertex: Option<TopoShape> = None;
        let mut apar1 = 0.0f64;
        let mut apar2 = 0.0f64;
        let vertex_found = find_common_vertex(
            the_ds,
            an_index1,
            an_index2,
            &self.int_point_cross_dir,
            &mut common_vertex,
            &mut apar1,
            &mut apar2,
        );

        let mut substitutor = MapReShape::new();

        // `cxx:467-517`: the common vertices between the `u`-iso edge pairs.
        let mut common_vertices: Vec<TopoShape> = Vec::new();
        let mut common_flag = 0i32;
        let low_col = 1usize;
        let up_col = u_edges.nb_cols();
        for (ueit, eindex) in [(1usize, the_index), (2usize, the_index + 1)] {
            let (Some(a_shape1), Some(a_shape2)) = (
                u_edges.value(eindex, low_col).cloned(),
                u_edges.value(eindex, up_col).cloned(),
            ) else {
                continue;
            };
            // `cxx:479-497`: an already-substituted `u`-iso edge wins.
            let a_ue1 = hist_first(&self.hist_map, &a_shape1).unwrap_or_else(|| a_shape1.clone());
            let a_ue2 = hist_first(&self.hist_map, &a_shape2).unwrap_or_else(|| a_shape2.clone());
            // `cxx:499-506`.
            if !is_same(&a_shape1, &a_ue1) {
                substitutor.replace(&forward(&a_shape1), &forward(&a_ue1));
            }
            if !is_same(&a_shape2, &a_ue2) {
                substitutor.replace(&forward(&a_shape2), &forward(&a_ue2));
            }
            // `cxx:509-517`: `LastVertex(aUE1)` / `FirstVertex(aUE2)`.
            let v1 = last_vertex(&Edge(a_ue1));
            let v2 = first_vertex(&Edge(a_ue2));
            if let (Some(v1), Some(v2)) = (v1, v2) {
                if is_same(&v1.0, &v2.0) {
                    common_flag = if common_flag == 0 { ueit as i32 } else { 3 };
                    common_vertices.push(v1.0);
                }
            }
        }

        // `cxx:520-537`.
        if vertex_found {
            if common_vertices.len() != 1 {
                return Ok(false);
            }
            let Some(cv) = common_vertex.clone() else {
                return Ok(false);
            };
            let cv = Vertex(cv);
            let first_cv = Vertex(common_vertices[0].clone());
            new_edge = if common_flag == 1 {
                make_edge_vv(&first_cv, &cv)
            } else {
                make_edge_vv(&cv, &first_cv)
            };
            has_new_edge = true;
        }

        // `cxx:539-546`.
        if common_vertices.len() == 2 {
            new_edge = make_edge_vv(
                &Vertex(common_vertices[0].clone()),
                &Vertex(common_vertices[1].clone()),
            );
            has_new_edge = true;
        }

        let builder = TopoBuilder::new();
        for fit in 1..=2usize {
            // `cxx:548-577`: assemble the compound whose free vertices are the
            // unmatched ends of the two `u`-iso edges.
            let mut comp_children: Vec<TopoShape> = Vec::new();
            let mut map_v: Vec<TopoShape> = Vec::new();
            for eindex in [the_index, the_index + 1] {
                let Some(a_shape) = u_edges.value(eindex, low_col + fit - 1).cloned() else {
                    continue;
                };
                let a_ue =
                    hist_first(&self.hist_map, &a_shape).unwrap_or_else(|| a_shape.clone());
                let a_v = if fit == 1 {
                    first_vertex(&Edge(a_ue.clone()))
                } else {
                    last_vertex(&Edge(a_ue.clone()))
                };
                if let Some(v) = a_v {
                    map_v.push(v.0);
                }
                comp_children.push(a_ue);
            }
            if has_new_edge {
                if let Some(ne) = new_edge.clone() {
                    comp_children.push(ne.0.clone());
                    store_vedge_in_hist_map(&self.v_edges, the_index, &ne.0, &mut self.hist_map);
                }
            }
            let comp = builder.make_compound_of(&comp_children);

            // `cxx:579-584`.
            let mut alone_vertices: Vec<TopoShape> = Vec::new();
            find_free_vertices(&comp.0, &map_v, &mut alone_vertices);
            if !alone_vertices.is_empty() && alone_vertices.len() != 2 {
                return Ok(false);
            }

            // `cxx:586-591`.
            let a_face_index = if fit == 1 {
                the_face_index1
            } else {
                the_face_index2
            };
            let Some(mut a_face) = the_ds.shape(a_face_index).cloned() else {
                return Ok(false);
            };
            let a_face_ori = a_face.orientation();
            a_face.set_orientation(Orientation::Forward);

            // `cxx:593-618`.
            if has_new_edge {
                if let Some(mut ne) = new_edge.clone() {
                    ne.0.set_orientation(Orientation::Forward);
                    // `cxx:596-612`: reverse the new edge when both of its ends
                    // sit closer to the *other* end of the reference bound edge.
                    let a_e_ref = Edge(if fit == 1 {
                        a_e1.clone()
                    } else {
                        a_e2.clone()
                    });
                    let (v1, v2) = (
                        first_vertex(&ne),
                        last_vertex(&ne),
                    );
                    let (r1, r2) = (
                        first_vertex(&a_e_ref),
                        last_vertex(&a_e_ref),
                    );
                    if let (Some(v1), Some(v2), Some(r1), Some(r2)) = (v1, v2, r1, r2) {
                        let p1 = BRepTool::vertex_point(&v1);
                        let p2 = BRepTool::vertex_point(&v2);
                        let q1 = BRepTool::vertex_point(&r1);
                        let q2 = BRepTool::vertex_point(&r2);
                        let p11 = p1.distance(&q1);
                        let p22 = p2.distance(&q2);
                        let p12 = p1.distance(&q2);
                        let p21 = p2.distance(&q1);
                        if p11 > p12 && p22 > p21 {
                            ne.0.reverse();
                        }
                    }
                    // `cxx:614-618`: add the pcurve on a non-planar face.
                    ShapeFixEdge::new().fix_add_pcurve(&ne, &a_face);
                    new_edge = Some(ne);
                }
            }

            // `cxx:620-726`.
            let mut ordered: Vec<Edge> = Vec::new();
            if alone_vertices.is_empty() {
                if let Some(ne) = new_edge.clone() {
                    ordered.push(ne);
                }
            } else {
                let an_e_index = if fit == 1 { an_index1 } else { an_index2 };
                let mut found1 = false;
                let mut found2 = false;
                let mut param1 = 0.0f64;
                let mut param2 = 0.0f64;
                // `cxx:629-657`: the two paves of the bound edge that carry the
                // free vertices.
                for pave in the_ds.paves(an_e_index) {
                    let Some(a_v) = the_ds.shape(pave.index()).cloned() else {
                        continue;
                    };
                    if is_same(&a_v, &alone_vertices[0]) && !found1 {
                        param1 = pave.parameter();
                        found1 = true;
                    }
                    if is_same(&a_v, alone_vertices.last().unwrap()) && !found2 {
                        param2 = pave.parameter();
                        found2 = true;
                    }
                }
                if !(found1 && found2) {
                    return Ok(false); // `cxx:716-719`.
                }
                // `cxx:659-688`: re-cut the bound edge between the two paves.
                let src = if fit == 1 { &a_e1 } else { &a_e2 };
                let (v_first, v_last, lo, hi) = if param1 < param2 {
                    (
                        alone_vertices[0].clone(),
                        alone_vertices[alone_vertices.len() - 1].clone(),
                        param1,
                        param2,
                    )
                } else {
                    (
                        alone_vertices[alone_vertices.len() - 1].clone(),
                        alone_vertices[0].clone(),
                        param2,
                        param1,
                    )
                };
                let mut v1 = Vertex(v_first);
                let mut v2 = Vertex(v_last);
                v1.0.set_orientation(Orientation::Forward);
                v2.0.set_orientation(Orientation::Reversed);
                let mut new_bound = empty_copied_edge(src);
                new_bound.set_orientation(Orientation::Forward);
                builder.add(&mut new_bound, &v1.0);
                builder.add(&mut new_bound, &v2.0);
                GeometryRegistry::global().set_ranges_all(&new_bound, lo, hi);
                ordered.push(Edge(new_bound));
                // `cxx:692-714`: prepend the new edge when it starts at `aV2`.
                if let Some(ne) = new_edge.clone() {
                    let starts_at_v2 = explorer_shapes(
                        &ne.0,
                        ShapeType::Vertex,
                    )
                    .iter()
                    .any(|v| is_same(v, &v2.0));
                    if starts_at_v2 {
                        ordered.push(ne);
                    } else {
                        ordered.insert(0, ne);
                    }
                }
            }

            // `cxx:729-749`.
            if !ordered.is_empty() {
                let mut wire = Wire::new();
                for e in &ordered {
                    builder.add(&mut wire.0, &e.0);
                }
                let target = if fit == 1 { &a_e1 } else { &a_e2 };
                substitutor.replace(&forward(target), &wire.0);
            }

            // `cxx:749-752`.
            substitutor.apply(&a_face);
            let new_face = substitutor.value(&a_face).unwrap_or_else(|| a_face.clone());
            let mut new_face = new_face;
            new_face.set_orientation(a_face_ori);
            self.hist_map
                .insert(shape_map_key(&a_face), vec![new_face]);

            // `cxx:753-779`: every edge of the source face that the
            // substitution rewrote is re-recorded as its list of new edges.
            for e in explorer_shapes(&a_face, ShapeType::Edge) {
                let Some(a_new_value) = substitutor.value(&e) else {
                    continue;
                };
                if is_same(&a_new_value, &e) {
                    continue;
                }
                if self.hist_map.contains_key(&shape_map_key(&e)) {
                    continue;
                }
                let list =
                    explorer_shapes(&a_new_value, ShapeType::Edge);
                self.hist_map.insert(shape_map_key(&e), list);
            }
        }

        Ok(true) // `cxx:782`.
    }

    /// `MakeFacesSec` (`cxx:789-1019`).
    /// `MakeFacesSec` (`cxx:789-1018`).
    ///
    /// The two face columns share a real section curve, so the bound edges are
    /// replaced by the section-edges chain that runs from the common vertex of
    /// the left `u`-iso pair to the common vertex of the right pair
    /// (`FindFromUEdge` then `FindFromVEdge`).
    pub(super) fn make_faces_sec(
        &mut self,
        the_index: usize,
        the_ds: &BopdsDS,
        the_face_index1: usize,
        the_face_index2: usize,
        the_ss_interf_index: usize,
    ) -> Result<bool, String> {
        // `cxx:792-797`.
        let a_b_curves = match the_ds.interf_ff().get(the_ss_interf_index) {
            Some(ff) => ff.curves().to_vec(),
            None => return Ok(false),
        };

        // `cxx:799-805`: `FilterSectionEdges(aBCurves, aSecPlane, theDS, aSecEdges)`.
        // `aSecPlane` is declared but never assigned, so it stays a null face
        // and the reference's `IntTools_BeanFaceIntersector` arm is dead.
        let mut a_sec_edges = filter_section_edges(&a_b_curves, the_ds);

        // `cxx:806-808`.
        let Some(u_edges) = self.u_edges.clone() else {
            return Ok(false);
        };
        let Some(bounds) = self.bounds.clone() else {
            return Ok(false);
        };

        // `cxx:808-813`.
        let (Some(left_e1), Some(left_e2), Some(right_e1), Some(right_e2)) = (
            u_edges.value(the_index, 1).cloned(),
            u_edges.value(the_index, 2).cloned(),
            u_edges.value(the_index + 1, 1).cloned(),
            u_edges.value(the_index + 1, 2).cloned(),
        ) else {
            return Ok(false);
        };

        // `cxx:815-818`.
        let (Some(index_left_e1), Some(index_left_e2), Some(index_right_e1), Some(index_right_e2)) = (
            the_ds.index(&left_e1),
            the_ds.index(&left_e2),
            the_ds.index(&right_e1),
            the_ds.index(&right_e2),
        ) else {
            return Ok(false);
        };

        // `cxx:819-838`: the two common vertices of the `u`-iso pairs.
        let mut first_vertex: Option<TopoShape> = None;
        let mut last_vertex: Option<TopoShape> = None;
        let mut _param_left_e1 = 0.0f64;
        let mut _param_left_e2 = 0.0f64;
        let mut _param_right_e1 = 0.0f64;
        let mut _param_right_e2 = 0.0f64;
        find_common_vertex(
            the_ds,
            index_left_e1,
            index_left_e2,
            &self.int_point_cross_dir,
            &mut first_vertex,
            &mut _param_left_e1,
            &mut _param_left_e2,
        );
        find_common_vertex(
            the_ds,
            index_right_e1,
            index_right_e2,
            &self.int_point_cross_dir,
            &mut last_vertex,
            &mut _param_right_e1,
            &mut _param_right_e2,
        );

        // `cxx:840-857`: `ChooseSection`; a found wire replaces the section
        // compound by its own edges and is recorded as the `v`-iso image.
        if let Some((sec_wire, _sec_plane, _is_singular)) = choose_section(
            self.transition,
            &self.axe_of_bis_plane,
            &a_sec_edges,
            first_vertex.as_ref(),
            last_vertex.as_ref(),
        ) {
            a_sec_edges = explorer_shapes(&sec_wire, ShapeType::Edge);
            store_vedge_in_hist_map(&self.v_edges, the_index, &sec_wire, &mut self.hist_map);
        }

        // `cxx:860`.
        let mut substitutor = MapReShape::new();

        for fit in 0..2usize {
            // `cxx:862-869`.
            let a_face_index = if fit == 0 {
                the_face_index1
            } else {
                the_face_index2
            };
            let Some(a_face) = the_ds.shape(a_face_index).cloned() else {
                return Ok(false);
            };
            let a_face_ori = a_face.orientation();
            let mut a_face_f = a_face.clone();
            a_face_f.set_orientation(Orientation::Forward);
            let Some(a_bound_edge) = bounds.value(the_index, 1 + fit).cloned() else {
                return Ok(false);
            };
            // `cxx:870`.
            let Some(a_bound_edge_index) = the_ds.index(&a_bound_edge) else {
                return Ok(false);
            };

            // `cxx:871-878`.
            let Some((a_ue1, a_ue2)) = get_u_edges(
                the_index,
                fit,
                &u_edges,
                &a_bound_edge,
                &a_face_f,
            ) else {
                return Ok(false);
            };
            let a_ue1_old = a_ue1.clone();
            let a_ue2_old = a_ue2.clone();

            // `cxx:882-911`: substitute each `u`-iso edge through the map,
            // matching the reference's `lst.First().Oriented(edge.Orientation())`.
            let mut a_ue1 = a_ue1;
            let mut a_ue2 = a_ue2;
            for (slot, old) in [(&mut a_ue1, &a_ue1_old), (&mut a_ue2, &a_ue2_old)] {
                if let Some(an_edge) = hist_first(&self.hist_map, &old.0) {
                    let an_edge = oriented_to(&an_edge, old.0.orientation());
                    if !is_same(&old.0, &an_edge) {
                        substitutor.replace(&forward(&old.0), &forward(&an_edge));
                    }
                    *slot = Edge(an_edge);
                }
            }

            // `cxx:912-919`.
            let mut comp_of_sec_edges = a_sec_edges.clone();
            let mut list_of_wire_edges: Vec<TopoShape> = Vec::new();

            // `cxx:920-958`.
            let Some(u_result) = find_from_u_edge(
                &a_ue1_old.0,
                &a_ue2_old.0,
                &a_ue1.0,
                &a_ue2.0,
                &a_face_f,
                &a_sec_edges,
                fit,
                &a_bound_edge,
                a_bound_edge_index,
                the_ds,
                &self.hist_map,
                comp_of_sec_edges,
                list_of_wire_edges,
            ) else {
                return Ok(false);
            };
            comp_of_sec_edges = u_result.comp_of_sec_edges;
            list_of_wire_edges = u_result.list_of_wire_edges;

            // `cxx:933-957`.
            let Some(v_result) = find_from_v_edge(
                &u_result.found_pave,
                u_result.is_on_u_edge,
                &a_ue1_old.0,
                &a_ue2_old.0,
                &a_face_f,
                &comp_of_sec_edges,
                fit,
                &a_bound_edge,
                a_bound_edge_index,
                the_ds,
                &self.hist_map,
            ) else {
                return Ok(false);
            };
            if !v_result.is_section_found && list_of_wire_edges.is_empty() {
                return Ok(false);
            }
            list_of_wire_edges.extend(v_result.list_of_wire_edges);

            // `cxx:965-981`: a raw `BRep_Builder` wire (`MakeWire` + `Add`), so
            // the edges keep the order and orientation of the list.
            if !list_of_wire_edges.is_empty() {
                let mut edges: Vec<Edge> = Vec::new();
                for s in &list_of_wire_edges {
                    if !is_same(&a_bound_edge, s) {
                        edges.push(Edge(s.clone()));
                    }
                }
                let a_w = TopoBuilder::new().make_wire(&edges);
                substitutor.replace(&forward(&a_bound_edge), &a_w.0);
            }

            // `cxx:984-988`.
            substitutor.apply(&a_face);
            let Some(mut new_face) = substitutor.value(&a_face) else {
                return Ok(false);
            };
            new_face.set_orientation(a_face_ori);
            self.hist_map
                .insert(shape_map_key(&a_face), vec![new_face]);

            // `cxx:990-1015`.
            for e in explorer_shapes(&a_face, ShapeType::Edge) {
                let Some(a_new_value) = substitutor.value(&e) else {
                    continue; // `cxx:996-999`.
                };
                if is_same(&a_new_value, &e) {
                    continue;
                }
                if self.hist_map.contains_key(&shape_map_key(&e)) {
                    continue; // `cxx:1001-1004`.
                }
                let list = explorer_shapes(&a_new_value, ShapeType::Edge);
                self.hist_map.insert(shape_map_key(&e), list);
            }
        }

        Ok(true) // `cxx:1018`.
    }
}
