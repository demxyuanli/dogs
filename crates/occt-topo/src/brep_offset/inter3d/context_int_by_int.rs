//! `ContextIntByInt` (`BRepOffset_Inter3d.cxx:1045-1272`): the intersections of
//! the context (un-offset) faces with the offset faces lying on the ancestors
//! of their edges and vertices, in integer (non-tangent) mode.
//!
//! Control flow, in the reference's order:
//!
//! 1. every context face becomes touched (`cxx:1069-1075`); with
//!    `ExtentContext`, `EnLargeFace(CF, NF, false, false)` caches its enlarged
//!    version in `MES`;
//! 2. per context face (`cxx:1080-1271`): the edges (and, in planar mode, the
//!    vertices) of `CF.Oriented(FORWARD)` fill `VEmap`;
//! 3. an edge with no ancestor in the initial shape (`cxx:1117-1165`) either
//!    joins the `AsDes` as it is, or is stretched by the file-local
//!    `ExtentEdge` (`:1137`), its two vertices are linked to the new edge, the
//!    vertex parameters are updated, and the new edge is cached in `MES` under
//!    the context edge;
//! 4. a vertex (`cxx:1167-1209`) contributes the ancestors of its non-generated
//!    ancestor edges that are not themselves in `VEmap`, whose faces carry no
//!    edge of `VEmap`;
//! 5. every ancestor face `F` of that list (`cxx:1212-1269`) takes
//!    `OF = MapSF(F).Face()` (enlarged once, cached in `MES`), then
//!    `Inter3D(WCF, NF, ..., OUT, E, CF, F)` runs unless `IsDone(NF, CF)`:
//!    a non-empty `LInt1` is stored in the `AsDes` and bound to `aS` in
//!    `Build` (the first edge alone, otherwise a compound together with the
//!    edges already bound), an empty one sends `aS` to `Failed`.
//!
//! `OE = MapSF(F).Generated(E)` (`cxx:1217-1218`) is read into a `LOE` list the
//! function never consumes, so the port does not call `Generated` at all.
//!
//! Port guards, each replacing an OCCT exception or a null-handle dereference,
//! never a branch:
//!
//! * `MapSF(F)` raising `Standard_NoSuchObject` (`cxx:1215`) and an offset face
//!   that was never built skip that ancestor face.
//! * `EnLargeFace` leaving `NF` null because `BRep_Tool::Surface` is null
//!   (`BRepOffset_Tool.cxx:3461`) skips that ancestor face.
//! * `MES(E)` on an unbindable key (`cxx:1156`) answers an empty shape, which
//!   links nothing instead of dereferencing.
//! * the `BRep_Builder::UpdateVertex` result (`cxx:1145-1148`) is discarded as
//!   in the reference; a vertex that cannot be updated is left as it was.

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_offset::{Analyse, BRepOffsetTool, EnLargeFaceOptions, ShapeMap};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key};

use super::extent_edge::extent_edge;
use super::{BRepOffsetInter3d, MapSf, ShapeIndexedMap};

impl BRepOffsetInter3d {
    /// `ContextIntByInt(ContextFaces, ExtentContext, MapSF, Analyse, MES,
    /// Build, Failed, theRange, bIsPlanar)` (`cxx:1045-1272`).
    pub fn context_int_by_int(
        &mut self,
        context_faces: &[TopoShape],
        extent_context: bool,
        map_sf: &MapSf,
        analyse: &Analyse,
        mes: &mut ShapeMap,
        build: &mut ShapeMap,
        failed: &mut Vec<TopoShape>,
        is_planar: bool,
    ) {
        let builder = TopoBuilder::new();

        // `cxx:1069-1075`: the touched faces and, with `ExtentContext`, their
        // enlarged versions.
        for cf in context_faces {
            self.touched.add(cf);
            if extent_context {
                if let Some((nf, _)) =
                    BRepOffsetTool::enlarge_face(&Face(cf.clone()), &EnLargeFaceOptions::default())
                {
                    mes.insert(shape_map_key(cf), nf.0.clone());
                }
            }
        }

        // `cxx:1077`: `TopAbs_State Side = TopAbs_OUT`.
        let side = State::Out;

        // `cxx:1080-1271`: one context face at a time.
        for cf_shape in context_faces {
            let cf = Face(cf_shape.clone());
            // `cxx:1087-1093`: `WCF` is the enlarged face when it exists.
            let wcf = if extent_context {
                match mes.get(&shape_map_key(cf_shape)) {
                    Some(nf) => Face(nf.clone()),
                    None => cf.clone(),
                }
            } else {
                cf.clone()
            };

            // `TopExp::MapShapes(CF.Oriented(FORWARD), ...)` (`cxx:1096-1103`).
            let cf_forward = cf.0.oriented(Orientation::Forward);
            let mut vemap = ShapeIndexedMap::new();
            for s in explorer_shapes(&cf_forward, ShapeType::Edge) {
                vemap.add(&s);
            }
            if is_planar {
                for s in explorer_shapes(&cf_forward, ShapeType::Vertex) {
                    vemap.add(&s);
                }
            }

            let n_ve = vemap.extent();
            for j in 1..=n_ve {
                let Some(a_s) = vemap.value(j).cloned() else {
                    continue;
                };
                let b_edge = a_s.shape_type() == ShapeType::Edge;
                // `E` is only assigned in the edge arm (`cxx:1114`); the vertex
                // arm passes the null edge to `Inter3D`.
                let mut e_ref: Option<Edge> = None;
                let mut anc: Vec<TopoShape> = Vec::new();

                if b_edge {
                    let e = Edge(a_s.clone());
                    if !analyse.has_ancestor(&e.0) {
                        // `cxx:1119-1165`: an edge of the context face that is
                        // not in the initial shape.
                        if !extent_context {
                            self.as_des.add(&cf.0, &e.0);
                            self.new_edges.add(&e.0);
                        } else if !mes.contains_key(&shape_map_key(&e.0)) {
                            let (f, l) = BRepTool::edge_parameters(&e);
                            let tol = BRepTool::edge_tolerance(&e);
                            let mut ne = extent_edge(&e, &cf);
                            let v = edge_vertices(&e);
                            ne.0.set_orientation(Orientation::Forward);
                            if let Some(v1) = v.0.as_ref() {
                                self.as_des.add(&ne.0, &v1.0.oriented(Orientation::Reversed));
                            }
                            if let Some(v2) = v.1.as_ref() {
                                self.as_des.add(&ne.0, &v2.0.oriented(Orientation::Forward));
                            }
                            // `B.UpdateVertex(V1, f, NE, Tol)` / `(V2, l, ...)`.
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
                            self.new_edges.add(&ne.0);
                            mes.insert(shape_map_key(&e.0), ne.0.clone());
                        } else if let Some(ne) = mes.get(&shape_map_key(&e.0)).cloned() {
                            self.as_des.add(&cf.0, &ne.oriented(e.0.orientation()));
                        }
                        continue;
                    }
                    e_ref = Some(e);
                    anc = analyse.ancestors(&a_s).to_vec();
                } else {
                    // `cxx:1167-1209`: the faces reached through the vertex.
                    if !analyse.has_ancestor(&a_s) {
                        continue;
                    }
                    for edge_shape in analyse.ancestors(&a_s).to_vec() {
                        if edge_shape.shape_type() != ShapeType::Edge {
                            continue;
                        }
                        let a_e = Edge(edge_shape);
                        if BRepTool::is_degenerated(&a_e) {
                            continue;
                        }
                        if vemap.contains(&a_e.0) {
                            continue;
                        }
                        for a_f in analyse.ancestors(&a_e.0).to_vec() {
                            let mut b_add = true;
                            for a_ef in explorer_shapes(&a_f, ShapeType::Edge) {
                                if vemap.contains(&a_ef) {
                                    b_add = false;
                                    break;
                                }
                            }
                            if b_add {
                                anc.push(a_f);
                            }
                        }
                    }
                }

                // `cxx:1212-1269`: every ancestor face of the list.
                for f_shape in anc {
                    let f = Face(f_shape);
                    let Some(of) = map_sf
                        .get(&shape_map_key(&f.0))
                        .and_then(|o| o.face())
                        .cloned()
                    else {
                        continue;
                    };
                    // `cxx:1220-1228`: `EnLargeFace(OF, NF, true, true)`, cached.
                    let nf = if let Some(nf) = mes.get(&shape_map_key(&of.0)) {
                        Face(nf.clone())
                    } else {
                        let opts = EnLargeFaceOptions {
                            can_extent_surface: true,
                            update_pcurve: true,
                            ..EnLargeFaceOptions::default()
                        };
                        let Some((nf, _)) = BRepOffsetTool::enlarge_face(&of, &opts) else {
                            continue;
                        };
                        mes.insert(shape_map_key(&of.0), nf.0.clone());
                        nf
                    };

                    if self.is_done(&nf, &cf) {
                        continue;
                    }
                    // `cxx:1231-1236`: `Inter3D(WCF, NF, ..., Side, E, CF, F)`.
                    let mut l_int1: Vec<TopoShape> = Vec::new();
                    let mut l_int2: Vec<TopoShape> = Vec::new();
                    BRepOffsetTool::inter3d(
                        &wcf,
                        &nf,
                        &mut l_int1,
                        &mut l_int2,
                        side,
                        e_ref.as_ref(),
                        Some(&cf),
                        Some(&f),
                    );
                    self.set_done(&nf, &cf);
                    if !l_int1.is_empty() {
                        self.store(&cf, &nf, &l_int1, &l_int2);
                        // `cxx:1239-1264`.
                        if l_int1.len() == 1 && !build.contains_key(&shape_map_key(&a_s)) {
                            build.insert(shape_map_key(&a_s), l_int1[0].clone());
                        } else {
                            let mut shapes: Vec<TopoShape> = Vec::new();
                            if let Some(se) = build.get(&shape_map_key(&a_s)) {
                                shapes.extend(explorer_shapes(se, ShapeType::Edge));
                            }
                            shapes.extend(l_int1.iter().cloned());
                            build.insert(shape_map_key(&a_s), builder.make_compound_of(&shapes).0);
                        }
                    } else {
                        failed.push(a_s.clone());
                    }
                }
            }
        }
    }
}
