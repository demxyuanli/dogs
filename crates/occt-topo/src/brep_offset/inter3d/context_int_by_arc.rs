//! `ContextIntByArc` (`BRepOffset_Inter3d.cxx:1276-1494`): the intersections of
//! the context (un-offset) faces with the offset faces lying on the ancestors
//! of their edges and vertices.

use crate::abs::{Orientation, ShapeType};
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_offset::{Analyse, BRepOffsetTool};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key, ShapeMapKey};
use std::collections::HashSet;

use super::extent_edge::extent_edge;
use super::BRepOffsetInter3d;

impl BRepOffsetInter3d {
    /// `ContextIntByArc(ContextFaces, InSide, Analyse, InitOffsetFace,
    /// InitOffsetEdge, theRange)` (`cxx:1276-1494`).
    ///
    /// `context_faces` is the `NCollection_IndexedMap` in index order.
    pub fn context_int_by_arc(
        &mut self,
        context_faces: &[TopoShape],
        in_side: bool,
        analyse: &Analyse,
        init_offset_face: &BRepAlgoImage,
        init_offset_edge: &mut BRepAlgoImage,
    ) {
        let mut l_int1: Vec<TopoShape> = Vec::new();
        let mut l_int2: Vec<TopoShape> = Vec::new();
        // `cxx:1286`: `NCollection_Map` of the vertices already treated.
        let mut mv: HashSet<ShapeMapKey> = HashSet::new();

        // `cxx:1294-1298`.
        for cf in context_faces {
            self.touched.add(cf);
        }

        for cf_shape in context_faces {
            let cf = Face(cf_shape.clone());

            // `cxx:1308`: `TopExp_Explorer(CF.Oriented(FORWARD), TopAbs_EDGE)`.
            let cf_forward = cf.0.oriented(Orientation::Forward);
            for e_shape in explorer_shapes(&cf_forward, ShapeType::Edge) {
                let e = Edge(e_shape);

                // `cxx:1311-1346`: the edge of the context face is not in the
                // initial shape.
                if !analyse.has_ancestor(&e.0) {
                    if in_side {
                        self.as_des.add(&cf.0, &e.0);
                    } else if !init_offset_edge.has_image(&e.0) {
                        let (f, l) = BRepTool::edge_parameters(&e);
                        let tol = BRepTool::edge_tolerance(&e);
                        let mut ne = extent_edge(&e, &cf);
                        // `TopExp::Vertices(E, V1, V2)` (default `CumOri`).
                        let v = edge_vertices(&e);
                        // `NE.Orientation(FORWARD)` for the AsDes links.
                        ne.0.set_orientation(Orientation::Forward);
                        if let Some(v1) = v.0.as_ref() {
                            self.as_des
                                .add(&ne.0, &v1.0.oriented(Orientation::Reversed));
                        }
                        if let Some(v2) = v.1.as_ref() {
                            self.as_des
                                .add(&ne.0, &v2.0.oriented(Orientation::Forward));
                        }
                        if let Some(v1) = v.0.as_ref() {
                            let _ = crate::brep_builder_update_vertex::update_vertex_on_edge(
                                &v1.0.oriented(Orientation::Internal),
                                f,
                                &ne.0,
                                tol,
                            );
                        }
                        if let Some(v2) = v.1.as_ref() {
                            let _ = crate::brep_builder_update_vertex::update_vertex_on_edge(
                                &v2.0.oriented(Orientation::Internal),
                                l,
                                &ne.0,
                                tol,
                            );
                        }
                        ne.0.set_orientation(e.0.orientation());
                        self.as_des.add(&cf.0, &ne.0);
                        let _ = init_offset_edge.bind(&e.0, &ne.0);
                    } else {
                        // `NE = InitOffsetEdge.Image(E).First()`.
                        if let Some(ne) = init_offset_edge.image(&e.0).first() {
                            self.as_des
                                .add(&cf.0, &ne.oriented(e.0.orientation()));
                        }
                    }
                    continue;
                }

                // `cxx:1350-1357`: `OF1` is the offset face of the first
                // ancestor of `E`, `OE` the offset edge of `E`.
                let Some(si) = analyse.ancestors(&e.0).first().cloned() else {
                    continue;
                };
                let Some(of1) = first_image(init_offset_face, &si) else {
                    continue;
                };
                let Some(oe) = first_edge_image(init_offset_edge, &e.0) else {
                    continue;
                };

                // `cxx:1359-1368`: `OE` must carry a pcurve on both faces.
                if !has_pcurve(&oe, &cf) || !has_pcurve(&oe, &of1) {
                    continue;
                }

                // `cxx:1379-1387`: `OE` itself is the intersection, oriented by
                // `OrientSection` with the first orientation reversed.
                l_int1.clear();
                l_int1.push(oe.0.clone());
                l_int2.clear();
                let Some((o1, _o2)) = BRepOffsetTool::orient_section(&oe, &cf, &of1) else {
                    continue;
                };
                l_int1[0].set_orientation(o1.reversed());
                self.store(&cf, &of1, &l_int1, &l_int2);

                // `cxx:1389-1447`: the offsets on the ancestors of the two
                // vertices of `E`.
                let v = edge_vertices(&e);
                for vi in [v.0, v.1].into_iter().flatten() {
                    if !mv.insert(shape_map_key(&vi.0)) {
                        continue;
                    }
                    let le: Vec<TopoShape> = analyse.ancestors(&vi.0).to_vec();
                    for ev_shape in &le {
                        if ev_shape.shape_type() != ShapeType::Edge {
                            continue;
                        }
                        if !init_offset_face.has_image(ev_shape) {
                            continue;
                        }
                        let Some(new_of1) = first_image(init_offset_face, ev_shape) else {
                            continue;
                        };
                        let Some(new_oe) = first_edge_image(init_offset_edge, &vi.0) else {
                            continue;
                        };
                        if !has_pcurve(&new_oe, &cf) || !has_pcurve(&new_oe, &new_of1) {
                            continue;
                        }
                        l_int1.clear();
                        l_int1.push(new_oe.0.clone());
                        l_int2.clear();
                        let Some((no1, _no2)) =
                            BRepOffsetTool::orient_section(&new_oe, &cf, &new_of1)
                        else {
                            continue;
                        };
                        l_int1[0].set_orientation(no1.reversed());
                        self.store(&cf, &new_of1, &l_int1, &l_int2);
                    }
                }
            }

            // `cxx:1450-1493`: the offset faces of the faces that carry an
            // edge ancestor of a vertex of `CF`.
            let cf_forward = cf.0.oriented(Orientation::Forward);
            for v_shape in explorer_shapes(&cf_forward, ShapeType::Vertex) {
                if !analyse.has_ancestor(&v_shape) {
                    continue;
                }
                let le: Vec<TopoShape> = analyse.ancestors(&v_shape).to_vec();
                for ev in &le {
                    if ev.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    let lf: Vec<TopoShape> = analyse.ancestors(ev).to_vec();
                    for fev in &lf {
                        if fev.shape_type() != ShapeType::Face {
                            continue;
                        }
                        let Some(of1) = first_image(init_offset_face, fev) else {
                            continue;
                        };
                        if self.is_done(&of1, &cf) {
                            continue;
                        }
                        // `cxx:1475-1481`: every edge of `OF1`.
                        let of1_forward = of1.0.oriented(Orientation::Forward);
                        let loe: Vec<Edge> = explorer_shapes(&of1_forward, ShapeType::Edge)
                            .into_iter()
                            .map(Edge)
                            .collect();
                        // `TryProject(CF, OF1, LOE, LInt1, LInt2, mySide, myTol)`
                        // overwrites `LInt1`/`LInt2`, so the port's narrower
                        // `Vec<Edge>` output is funnelled back into them.
                        let mut lp1: Vec<Edge> = Vec::new();
                        let mut lp2: Vec<Edge> = Vec::new();
                        let projected = BRepOffsetTool::try_project(
                            &cf, &of1, &loe, &mut lp1, &mut lp2, self.side, self.tol,
                        );
                        l_int1 = lp1.into_iter().map(|e| e.0).collect();
                        l_int2 = lp2.into_iter().map(|e| e.0).collect();
                        if !projected || l_int1.is_empty() {
                            BRepOffsetTool::inter3d(
                                &cf, &of1, &mut l_int1, &mut l_int2, self.side, None, None, None,
                            );
                        }
                        self.store(&cf, &of1, &l_int1, &l_int2);
                    }
                }
            }
        }
    }
}

/// `TopoDS::Face(Image(S).First())`.
fn first_image(image: &BRepAlgoImage, s: &TopoShape) -> Option<Face> {
    image.image(s).first().map(|i| Face(i.clone()))
}

/// `TopoDS::Edge(Image(S).First())`.
fn first_edge_image(image: &BRepAlgoImage, s: &TopoShape) -> Option<Edge> {
    image.image(s).first().map(|i| Edge(i.clone()))
}

/// `BRep_Tool::CurveOnSurface(E, F, f, l).IsNull()`; the port's helper answers
/// `None` for a missing pcurve.
fn has_pcurve(e: &Edge, f: &Face) -> bool {
    crate::boptools_2d::curve_on_surface_range(e, f).is_some()
}
