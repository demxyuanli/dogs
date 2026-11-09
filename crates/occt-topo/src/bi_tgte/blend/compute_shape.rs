//! `BiTgte_Blend::ComputeShape` (`BiTgte_Blend.cxx:2221-2507`): assemble the
//! result shape - the untouched faces, the tubes, and the reconstructions of
//! the faces the unwinding cut.
//!
//! Control flow, all from the `.cxx`:
//!
//! * `:2243-2292` invariant edges: for a face that is in `myFaces`, bound in
//!   `myMapSF` and whose offset face survived the unwinding, an offset edge
//!   that unwinds onto itself (`Image(CurOE).First().IsSame(CurOE)`) is
//!   recorded in `myCreated(CurOE)(CurF) = {CurE}`.
//! * `:2296-2305` a face outside `myFaces` is copied to the result untouched.
//! * `:2306-2478` otherwise every piece `LastImage(CurOF)` gives is rebuilt on
//!   the *source* surface of `CurF`: `:2333-2337` `MakeFace` + `UpdateFace`
//!   with `CurF`'s surface and tolerance, `:2343-2377` the vertex images of the
//!   cut edges are recorded in `myCreated`, `:2380-2471` every wire of the
//!   piece is rebuilt edge by edge (`UpdateEdge` with the two p-curves of a
//!   cut edge, or with the recorded / newly made edge and vertices), then
//!   `:2470-2475` the piece orientation is copied and the face added.
//! * `:2480-2505` the result is sewn again, `SameParameter` runs on every edge
//!   of the sewed shape, and `myResult` becomes the sewed shape.
//!
//! Deviations / guards, none of them branches of their own:
//! - `myCreated(E)` bound but `(CurF)` not (`cxx:2416-2419` / `:2425-2428` /
//!   `:2446-2449`): OCCT leaves `OE` / `OV1` / `OV2` null and dereferences
//!   them; the port builds them, i.e. runs the `else` arm.
//! - a null `BRep_Tool::CurveOnSurface` (`cxx:2391`) throws in the reference;
//!   the port drops that edge from the rebuilt wire.
//! - `BRep_Tool::Surface(CurF, L)` answering null (`cxx:2335`) throws; the port
//!   skips the piece.
//! - `BRepTools::Update` (`cxx:2473`) has no counterpart (no cached UV points).
//! - the local `Created` map of `cxx:2233` is never read in the reference, so
//!   the port keeps no counterpart.

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::brep_builder_api_sewing::Sewing;
use crate::brep_lib_build_curves3d::build_curves3d;
use crate::brep_lib_same_parameter::same_parameter_shape;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::pcurve_on_face;
use super::Blend;

impl Blend {
    /// `void BiTgte_Blend::ComputeShape()` (`cxx:2221-2507`).
    pub(super) fn compute_shape(&mut self) {
        let builder = TopoBuilder::new();
        let Some(shape) = self.my_shape.clone() else {
            return;
        };

        // ---------------------------------------------------------------
        // `cxx:2243-2292`: invariant edges stay where they are.
        // ---------------------------------------------------------------
        for cur_f in explorer_shapes(&shape, ShapeType::Face) {
            if !self.my_faces.contains(&cur_f) {
                continue; // the face is not touched
            }
            if !self.my_map_sf.is_bound(&cur_f) {
                continue; // inverted or degenerated
            }
            let Some(offset) = self.my_map_sf.get(&cur_f).cloned() else {
                continue;
            };
            let Some(cur_of) = offset.face().cloned() else {
                continue;
            };
            if !self.my_image_offset.has_image(&cur_of.0) {
                continue; // the face disappears in the unwinding
            }

            for cur_e in explorer_shapes(&cur_f, ShapeType::Edge) {
                let Some(cur_oe) = offset.generated(&cur_e) else {
                    continue;
                };
                if !self.my_image_offset.has_image(&cur_oe) {
                    continue;
                }
                let image = self.my_image_offset.image(&cur_oe);
                let Some(im_e) = image.first() else {
                    continue;
                };
                if is_same(im_e, &cur_oe) {
                    self.my_created.bind_empty(&cur_oe);
                    self.my_created.inner_bind_empty(&cur_oe, &cur_f);
                    self.my_created.inner_append(&cur_oe, &cur_f, &cur_e);
                }
            }
        }

        // ---------------------------------------------------------------
        // `cxx:2296-2478`: untouched faces, then the rebuilt pieces.
        // ---------------------------------------------------------------
        for cur_f in explorer_shapes(&shape, ShapeType::Face) {
            if !self.my_faces.contains(&cur_f) {
                if let Some(result) = self.my_result.as_mut() {
                    builder.add(result, &cur_f);
                }
                continue;
            }

            if !self.my_map_sf.is_bound(&cur_f) {
                continue; // inverted or degenerated
            }
            let Some(cur_of) = self.my_map_sf.get(&cur_f).and_then(|of| of.face().cloned()) else {
                continue;
            };
            if !self.my_image_offset.has_image(&cur_of.0) {
                continue; // the face disappears in the unwinding
            }

            let mut lim: Vec<TopoShape> = Vec::new();
            self.my_image_offset.last_image(&cur_of.0, &mut lim);

            for debouc in lim {
                let cur_f_shape = Face(cur_f.clone());
                let Some(s) = BRepTool::face_surface(&cur_f_shape) else {
                    continue;
                };
                // `cxx:2333-2337`.
                let mut new_f = builder.make_face_with_tolerance(
                    s.clone(),
                    &[],
                    BRepTool::face_tolerance(&cur_f_shape),
                );

                // `cxx:2340-2342`.
                let face = Face(debouc.oriented(Orientation::Forward));

                // `cxx:2343-2377`: record the vertex images of the cut edges.
                for e in explorer_shapes(&face.0, ShapeType::Edge) {
                    if !self.my_created.is_bound(&e) || !self.my_created.inner_is_bound(&e, &cur_f) {
                        continue;
                    }
                    let Some(oe) = self.my_created.inner_first(&e, &cur_f).cloned() else {
                        continue;
                    };
                    let (v1, v2) = edge_vertices(&Edge(e.clone()));
                    let (ov1, ov2) = edge_vertices(&Edge(oe));
                    for (v, ov) in [(v1, ov1), (v2, ov2)] {
                        let (Some(v), Some(ov)) = (v, ov) else {
                            continue;
                        };
                        if !self.my_created.is_bound(&v.0) {
                            self.my_created.bind_empty(&v.0);
                        }
                        if !self.my_created.inner_is_bound(&v.0, &cur_f) {
                            self.my_created.inner_bind_empty(&v.0, &cur_f);
                            self.my_created.inner_append(&v.0, &cur_f, &ov.0);
                        }
                    }
                }

                // `cxx:2380-2471`: rebuild every wire.
                let mut map_ss: HashMap<ShapeMapKey, TopoShape> = HashMap::new();
                let wires = explorer_shapes(&face.0, ShapeType::Wire);
                for w in wires {
                    let mut new_edges: Vec<Edge> = Vec::new();
                    for e in explorer_shapes(&w.oriented(Orientation::Forward), ShapeType::Edge) {
                        let e = Edge(e);
                        let Some((c2d, _surf, first, last)) = pcurve_on_face(&e, &face) else {
                            continue;
                        };
                        let tol_e = BRepTool::edge_tolerance(&e);

                        let existing = map_ss.get(&shape_map_key(&e.0)).cloned();
                        let oe = if let Some(existing) = existing {
                            // `cxx:2394-2411`: a cut edge.
                            let oe = Edge(existing);
                            let rev = Edge(e.0.oriented(Orientation::Reversed));
                            let Some((c2d_1, _, _, _)) = pcurve_on_face(&rev, &face) else {
                                continue;
                            };
                            let curves = if e.0.orientation() == Orientation::Forward {
                                vec![c2d.clone(), c2d_1]
                            } else {
                                vec![c2d_1, c2d.clone()]
                            };
                            GeometryRegistry::global().set_edge_pcurves_on_surface(
                                &oe.0,
                                &s,
                                curves,
                                tol_e,
                            );
                            GeometryRegistry::global().set_ranges_all(&oe.0, first, last);
                            oe
                        } else {
                            // `cxx:2413-2465`.
                            let recorded = if self.my_created.is_bound(&e.0) {
                                self.my_created.inner_first(&e.0, &cur_f).cloned()
                            } else {
                                None
                            };
                            let oe = if let Some(rec) = recorded {
                                Edge(rec)
                            } else {
                                let mut oe = Edge::new();
                                // `cxx:2423-2459`: the endpoint vertices.
                                let (v1, v2) = edge_vertices(&e);
                                let mut out_vertices: Vec<Option<TopoShape>> = Vec::new();
                                for v in [v1.clone(), v2.clone()].into_iter().flatten() {
                                    let recorded_v = if self.my_created.is_bound(&v.0)
                                        && self.my_created.inner_is_bound(&v.0, &cur_f)
                                    {
                                        self.my_created.inner_first(&v.0, &cur_f).cloned()
                                    } else {
                                        None
                                    };
                                    if let Some(rv) = recorded_v {
                                        out_vertices.push(Some(rv));
                                        continue;
                                    }
                                    // `cxx:2435-2442` / `:2452-2459`.
                                    let u = BRepTool::parameter_of_vertex_on_face(&v.0, &e, &face);
                                    let p2d = c2d.d0(u);
                                    let p = s.d0(p2d.x(), p2d.y()).transformed(
                                        &cur_f_shape.0.location().transformation(),
                                    );
                                    let nv = builder.make_vertex(p, BRepTool::vertex_tolerance(&v));
                                    self.my_created.bind_empty(&v.0);
                                    self.my_created.inner_bind_empty(&v.0, &cur_f);
                                    self.my_created.inner_append(&v.0, &cur_f, &nv.0);
                                    out_vertices.push(Some(nv.0));
                                }
                                let (ov1, ov2) = (
                                    out_vertices.first().cloned().flatten(),
                                    out_vertices.get(1).cloned().flatten(),
                                );
                                if let Some(v1) = v1.as_ref() {
                                    if let Some(ov1) = ov1.as_ref() {
                                        builder
                                            .add(&mut oe.0, &ov1.oriented(v1.0.orientation()));
                                    }
                                }
                                if let Some(v2) = v2.as_ref() {
                                    if let Some(ov2) = ov2.as_ref() {
                                        builder
                                            .add(&mut oe.0, &ov2.oriented(v2.0.orientation()));
                                    }
                                }
                                oe
                            };
                            GeometryRegistry::global().set_edge_pcurve(
                                &oe.0,
                                GeometryRegistry::shape_key(&new_f.0),
                                c2d,
                            );
                            GeometryRegistry::global().set_ranges_all(&oe.0, first, last);
                            map_ss.insert(shape_map_key(&e.0), oe.0.clone());
                            oe
                        };
                        new_edges.push(Edge(oe.0.oriented(e.0.orientation())));
                    }
                    // `cxx:2469`.
                    let new_wire = builder.make_wire(&new_edges);
                    let oriented_wire = new_wire.0.oriented(w.orientation());
                    builder.add(&mut new_f.0, &oriented_wire);
                }

                // `cxx:2472-2475`.
                new_f.0.set_orientation(debouc.orientation());
                if let Some(result) = self.my_result.as_mut() {
                    builder.add(result, &new_f.0);
                }
            }
        }

        // ---------------------------------------------------------------
        // `cxx:2481-2505`: sew the result once more.
        // ---------------------------------------------------------------
        let mut sew = Sewing::new(self.my_tol, true, true, true, false);
        let Some(result) = self.my_result.clone() else {
            return;
        };
        let _ = build_curves3d(&result);
        for f in explorer_shapes(&result, ShapeType::Face) {
            sew.add(&f);
        }
        if sew.perform().is_err() {
            return;
        }
        let Some(sewed) = sew.sewed_shape().cloned() else {
            return;
        };
        for e in explorer_shapes(&sewed, ShapeType::Edge) {
            same_parameter_shape(&e, BRepTool::edge_tolerance(&Edge(e.clone())), false);
        }
        self.my_result = Some(sewed);
    }
}




