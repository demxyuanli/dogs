//! `BRepOffset_MakeLoops::BuildFaces` (`BRepOffset_MakeLoops.cxx:455-661`).

use std::collections::HashMap;

use crate::abs::Orientation;
use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_algo_loop::BRepAlgoLoop;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same, shape_map_key, ShapeMapKey};

use super::{substitute_vertices, BRepOffsetMakeLoops};

impl BRepOffsetMakeLoops {
    /// `BuildFaces(LF, AsDes, Image, theRange)` (`cxx:455-661`): rebuild only the
    /// faces that hold at least one already-cut edge, after mapping the old
    /// vertices of those edges onto the new ones.
    pub fn build_faces(
        &mut self,
        lf: &[TopoShape],
        as_des: &mut AsDes,
        image: &mut BRepAlgoImage,
    ) {
        let mut loops = BRepAlgoLoop::new();
        // `cxx:464`: `Loops.VerticesForSubstitute(myVerVerMap)`.
        loops.set_vertices_for_substitute(self.my_ver_ver_map.clone());

        for f in lf {
            let face = Face(f.clone());
            loops.init(&face);
            let mut to_rebuild = false;
            let mut added_edges: Vec<TopoShape> = Vec::new();

            // `cxx:481`: only an unimaged face is examined.
            if !image.has_image(&face.0) {
                let le: Vec<TopoShape> = as_des.descendant(&face.0).to_vec();

                // `NCollection_DataMap MONV` (`cxx:491`): old vertex -> new
                // vertex, `IsSame` identity.
                let mut monv: HashMap<ShapeMapKey, TopoShape> = HashMap::new();

                // `cxx:494-531`: the first pass, which decides `ToRebuild` and
                // fills `MONV`.
                for e_shape in &le {
                    let e = Edge(e_shape.clone());
                    if !image.has_image(&e.0) {
                        continue;
                    }
                    let lce: Vec<TopoShape> = image.image(&e.0);
                    if lce.len() == 1 && is_same(&lce[0], &e.0) {
                        // `cxx:501-508`: the edge maps onto itself.
                        let ce = lce[0].clone().oriented(e.0.orientation());
                        loops.add_const_edge(&Edge(ce));
                        continue;
                    }
                    // `cxx:509-530`: `F` has to be reconstructed.
                    to_rebuild = true;
                    for ce_shape in &lce {
                        let ce_shape = ce_shape.clone().oriented(e.0.orientation());
                        let ce = Edge(ce_shape);
                        // `TopExp::Vertices(E, OV1, OV2)` /
                        // `TopExp::Vertices(CE, NV1, NV2)` (`:518-519`), both
                        // with the default `CumOri = false`.
                        let (ov1, ov2) = edge_vertices(&e);
                        let (nv1, nv2) = edge_vertices(&ce);
                        // `cxx:520-527`: a null endpoint would make OCCT bind a
                        // null shape and later call `UpdateVertex` on it; the
                        // port only maps real pairs.
                        if let (Some(ov), Some(nv)) = (ov1.as_ref(), nv1.as_ref()) {
                            if !is_same(&ov.0, &nv.0) {
                                monv.insert(shape_map_key(&ov.0), nv.0.clone());
                            }
                        }
                        if let (Some(ov), Some(nv)) = (ov2.as_ref(), nv2.as_ref()) {
                            if !is_same(&ov.0, &nv.0) {
                                monv.insert(shape_map_key(&ov.0), nv.0.clone());
                            }
                        }
                        loops.add_const_edge(&ce);
                    }
                }

                if to_rebuild {
                    // `cxx:539-576`: the edges that were not reconstructed, with
                    // the vertices `MONV` substitutes.
                    for e_shape in &le {
                        let e = Edge(e_shape.clone());
                        // `BRep_Tool::Range(E, f, l)` (`:542-543`).
                        let (f_par, l_par) = BRepTool::edge_parameters(&e);
                        if image.has_image(&e.0) {
                            continue;
                        }
                        let (ov1, ov2) = edge_vertices(&e);
                        let mut lv: Vec<TopoShape> = Vec::new();
                        if let Some(vv) = ov1
                            .as_ref()
                            .and_then(|ov| monv.get(&shape_map_key(&ov.0)))
                        {
                            // `cxx:550-555`: FORWARD, then `UpdateVertex` at the
                            // edge's first parameter.
                            let vv = vv.clone().oriented(Orientation::Forward);
                            let tol = BRepTool::vertex_tolerance(&Vertex(vv.clone()));
                            let _ = crate::brep_builder_update_vertex::update_vertex_on_edge(
                                &vv.oriented(Orientation::Internal),
                                f_par,
                                &e.0,
                                tol,
                            );
                            lv.push(vv);
                        }
                        if let Some(vv) = ov2
                            .as_ref()
                            .and_then(|ov| monv.get(&shape_map_key(&ov.0)))
                        {
                            // `cxx:556-565`: REVERSED, at the last parameter.
                            let vv = vv.clone().oriented(Orientation::Reversed);
                            let tol = BRepTool::vertex_tolerance(&Vertex(vv.clone()));
                            let _ = crate::brep_builder_update_vertex::update_vertex_on_edge(
                                &vv.oriented(Orientation::Internal),
                                l_par,
                                &e.0,
                                tol,
                            );
                            lv.push(vv);
                        }
                        // `cxx:566-576`.
                        if lv.is_empty() {
                            loops.add_const_edge(&e);
                        } else {
                            loops.add_edge(&e, &lv);
                            added_edges.push(e.0.clone());
                        }
                    }
                }
            }

            // `cxx:579-611`: reconstruction.
            if to_rebuild {
                loops.perform();
                loops.wires_to_faces();

                let nf: Vec<TopoShape> = loops.new_faces().to_vec();
                let _ = image.bind_list(&face.0, &nf);

                for e_shape in &added_edges {
                    let e = Edge(e_shape.clone());
                    let Some(new_edges) = loops.new_edges(&e).cloned() else {
                        continue;
                    };
                    if image.has_image(&e.0) {
                        let _ = image.add_list(&e.0, &new_edges);
                    } else {
                        let _ = image.bind_list(&e.0, &new_edges);
                    }
                }
            }
        }

        // `cxx:613-616`.
        self.my_ver_ver_map = loops.vertices_for_substitute().clone();
        if self.my_ver_ver_map.is_empty() {
            return;
        }

        // `cxx:617-660`.
        let map = self.my_ver_ver_map.clone();
        substitute_vertices(lf, image, as_des, &map);
    }
}
