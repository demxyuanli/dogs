//! `BiTgte_Blend::ComputeSurfaces` (`BiTgte_Blend.cxx:1761-2217`): build one
//! tube per unwound centre edge and one sphere per vertex where three tubes
//! meet.
//!
//! Control flow, all from the `.cxx`:
//!
//! * `:1764` `myFaces.Clear()`; `:1780` `myResult` becomes the compound the
//!   tubes are added to; `:1792` `Co` collects the centre lines.
//! * `:1802-1815` for every edge of `myEdges` whose two ascendants are faces
//!   `F1` / `F2`; `:1824-1825` their orientation in the two faces (`Orientation`
//!   against the descendant lists); `:1832-1833` `Or1` / `Or2` from
//!   `myInitOffsetFace`.
//! * `:1836-1884` the support shape of each side is either the rolled edge
//!   (`OF1isEdge`: its 3D curve, location applied) or the offset face (its
//!   surface); a vertex support stops the whole edge (`:1862-1865`,
//!   `:1884-1887`).
//! * `:1891-1907` the cube pieces: `myImageOffset.LastImage(CurE)`, or `CurE`
//!   itself when both supports are edges (a tube-tube intersection line with no
//!   image is kept as is).
//! * `:1911-2038` for each piece: its p-curves on `F1` / `F2`, the two
//!   `FindCreatedEdge` calls and their vertex sets, then `E1` / `E2` built
//!   either from `MakeCurve` (edge support) or from the p-curve and surface
//!   (`BRepLib_MakeEdge` + `KPartCurve3d`).
//! * `:2041-2050` `myCreated` records `E1` / `E2` under `CurCutE`.
//! * `:2057-2094` when `FindCreatedEdge` missed a tangent edge, the vertex set
//!   is scanned for an edge joining the two found vertices.
//! * `:2096-2100` the orientations, then
//!   `BRepOffset_Offset(CurCutE, E1, E2, -myRadius, E1f, E1l, myNubs, myTol,
//!   GeomAbs_C2)` into `myMapSF`, `myCenters` and the result.
//! * `:2108-2178` with `myBuildShape`: the tube's first / second edge vertices
//!   are projected onto the support curve to register their parameters on
//!   `Or1` / `Or2` (`BRep_Builder::UpdateVertex`) and appended to `myCutEdges`.
//! * `:2185-2216` the spheres: every vertex of `Co` with exactly three
//!   incident centre edges gets `BRepOffset_Offset(V, LOE, -myRadius, ...)`
//!   with the three generated edges reversed.
//!
//! Deviations / guards, none of them branches of their own:
//! - `BRep_Tool::Surface(OF)` raising on a null surface and
//!   `TopoDS::Face(L.First())` raising on a vertex support map to `continue`
//!   (`cxx:1815`, `:1832`, `:1862`), which is what the reference does for the
//!   filtered cases around them.
//! - `BRepLib_MakeEdge` failing (`cxx:1983`, `:2033`) answers a null edge in
//!   the reference and `KPartCurve3d` then dereferences it; the port drops that
//!   cube piece instead.
//! - `GeomAPI_ProjectPointOnCurve` failing (`cxx:2120`, `:2123`, `:2151`,
//!   `:2155`) throws in the reference; the port skips that vertex update.

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;
use occt_core::precision::CONFUSION;
use occt_geom::geom_api::project_point_on_curve;
use occt_geom::{Curve, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_offset::BRepOffsetOffset;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same, shape_map_key, AncestorMap};

use super::make_curve::make_curve;
use super::statics::{
    curve_on_edge, find_created_edge, find_vertex, kpart_curve3d, make_degenerated_edge,
    orientation, EdgeSet,
};
use super::{make_edge_2d, make_edge_3d, pcurve_on_face, Blend};

impl Blend {
    /// `void BiTgte_Blend::ComputeSurfaces()` (`cxx:1761-2217`).
    pub(super) fn compute_surfaces(&mut self) {
        // `cxx:1764`.
        self.my_faces = crate::brep_offset::ShapeIndexedMap::new();

        // `cxx:1776-1780`.
        let tol_angle = 2.0 * (self.my_tol / (self.my_radius * 0.5).abs()).asin();
        let center_analyse = match self.my_result.clone() {
            Some(result) => crate::brep_offset::Analyse::with_shape(&result, tol_angle).ok(),
            None => None,
        };
        let Some(center_analyse) = center_analyse else {
            return;
        };

        // `cxx:1785-1792`.
        let builder = TopoBuilder::new();
        self.my_result = Some(builder.make_shape(ShapeType::Compound));
        let mut co = builder.make_shape(ShapeType::Compound);

        // `cxx:1802-1815`.
        let nb_edges = self.my_edges.extent();
        for i in 0..nb_edges {
            let cur_e = Edge(self.my_edges.value(i).expect("myEdges(i)").clone());

            let l = self.my_as_des.ascendant(&cur_e.0).to_vec();
            if l.len() != 2 {
                continue;
            }
            let f1 = Face(l[0].clone());
            let f2 = Face(l[1].clone());

            // `cxx:1824-1825`.
            let ld1 = self.my_as_des.descendant(&f1.0).to_vec();
            let ld2 = self.my_as_des.descendant(&f2.0).to_vec();
            let orien1 = orientation(&cur_e, &f1, &ld1);
            let orien2 = orientation(&cur_e, &f2, &ld2);

            // `cxx:1832-1833`.
            let (Ok(or1), Ok(or2)) = (
                self.my_init_offset_face.image_from(&f1.0),
                self.my_init_offset_face.image_from(&f2.0),
            ) else {
                continue;
            };

            // `cxx:1836-1837`.
            self.my_faces.add(&or1);
            self.my_faces.add(&or2);

            // `cxx:1839-1884`.
            let mut of1_is_edge = false;
            let mut of2_is_edge = false;
            let mut oe1: Option<Edge> = None;
            let mut oe2: Option<Edge> = None;
            let mut gc1: Option<Arc<dyn Curve>> = None;
            let mut gc2: Option<Arc<dyn Curve>> = None;
            let mut of1: Option<Face> = None;
            let mut of2: Option<Face> = None;
            let mut gs1: Option<Arc<dyn Surface>> = None;
            let mut gs2: Option<Arc<dyn Surface>> = None;

            if or1.shape_type() == ShapeType::Edge {
                of1_is_edge = true;
                let e = Edge(or1.clone());
                // `cxx:1848-1849`: the curve with `Loc` applied.
                gc1 = BRepTool::edge_curve_world(&e);
                if gc1.is_none() {
                    continue;
                }
                oe1 = Some(e);
            } else if or1.shape_type() == ShapeType::Face {
                let f = Face(or1.clone());
                gs1 = BRepTool::face_surface(&f);
                if gs1.is_none() {
                    continue;
                }
                of1 = Some(f);
            }

            // `cxx:1862-1865`.
            if of1.is_none() && oe1.is_none() {
                continue;
            }

            if or2.shape_type() == ShapeType::Edge {
                of2_is_edge = true;
                let e = Edge(or2.clone());
                gc2 = BRepTool::edge_curve_world(&e);
                if gc2.is_none() {
                    continue;
                }
                oe2 = Some(e);
            } else if or2.shape_type() == ShapeType::Face {
                let f = Face(or2.clone());
                gs2 = BRepTool::face_surface(&f);
                if gs2.is_none() {
                    continue;
                }
                of2 = Some(f);
            }

            // `cxx:1884-1887`.
            if of2.is_none() && oe2.is_none() {
                continue;
            }

            // `cxx:1889-1907`.
            let mut cur_l: Vec<TopoShape> = Vec::new();
            if !self.my_image_offset.has_image(&cur_e.0) {
                if of1_is_edge && of2_is_edge {
                    cur_l.push(cur_e.0.clone());
                } else {
                    continue;
                }
            } else {
                self.my_image_offset.last_image(&cur_e.0, &mut cur_l);
            }

            // `cxx:1910-1912`.
            for item in cur_l {
                let cur_cut_e = Edge(item);

                // `cxx:1914-1923`.
                let (Some((pc1, _, f1p, l1p)), Some((pc2, _, f2p, l2p))) = (
                    pcurve_on_face(&cur_cut_e, &f1),
                    pcurve_on_face(&cur_cut_e, &f2),
                ) else {
                    continue;
                };

                // `cxx:1925-1937`.
                let mut map_on_v1f = EdgeSet::new();
                let mut map_on_v1l = EdgeSet::new();
                let (v1f, v1l) = edge_vertices(&cur_cut_e);
                let mut e1f: Option<Edge> = None;
                let mut e1l: Option<Edge> = None;
                if let Some(v1f) = v1f.as_ref() {
                    e1f = find_created_edge(
                        v1f,
                        &cur_cut_e,
                        &self.my_map_sf,
                        &mut map_on_v1f,
                        &center_analyse,
                        self.my_radius,
                        self.my_tol,
                    );
                }
                if let Some(v1l) = v1l.as_ref() {
                    e1l = find_created_edge(
                        v1l,
                        &cur_cut_e,
                        &self.my_map_sf,
                        &mut map_on_v1l,
                        &center_analyse,
                        self.my_radius,
                        self.my_tol,
                    );
                }

                // `cxx:1940-1988`: the first support edge.
                let (mut e1, vf_on_e1, vl_on_e1): (Edge, Option<Vertex>, Option<Vertex>) = if of1_is_edge {
                    let Some(oe1) = oe1.as_ref() else {
                        continue;
                    };
                    let cone = curve_on_edge(&cur_cut_e, oe1);
                    let c = make_curve(&cone);
                    let p1 = c.d0(c.first_parameter());
                    let p2 = c.d0(c.last_parameter());
                    let vf = find_vertex(&p1, &map_on_v1f, self.my_tol)
                        .or_else(|| find_vertex(&p1, &map_on_v1l, self.my_tol));
                    let vl = find_vertex(&p2, &map_on_v1l, self.my_tol)
                        .or_else(|| find_vertex(&p2, &map_on_v1f, self.my_tol));
                    let e = if p1.square_distance(&p2) < self.my_tol * self.my_tol {
                        make_degenerated_edge(&c, vf.as_ref())
                    } else {
                        let Some(edge) = make_edge_3d(&c, vf.as_ref(), vl.as_ref()) else {
                            continue;
                        };
                        edge
                    };
                    (e, vf, vl)
                } else {
                    let (Some(gs1), Some(of1)) = (gs1.as_ref(), of1.as_ref()) else {
                        continue;
                    };
                    let p2d = pc1.d0(f1p);
                    let p1 = gs1.d0(p2d.x(), p2d.y());
                    let p2d = pc1.d0(l1p);
                    let p2 = gs1.d0(p2d.x(), p2d.y());
                    let vf = find_vertex(&p1, &map_on_v1f, self.my_tol);
                    let vl = find_vertex(&p2, &map_on_v1l, self.my_tol);
                    let Some(e) = make_edge_2d(
                        &pc1,
                        gs1,
                        of1,
                        vf.as_ref(),
                        vl.as_ref(),
                        f1p,
                        l1p,
                    ) else {
                        continue;
                    };
                    kpart_curve3d(&e, &pc1, gs1);
                    (e, vf, vl)
                };

                // `cxx:1991-2038`: the second support edge.
                let (mut e2, vf_on_e2, vl_on_e2): (Edge, Option<Vertex>, Option<Vertex>) = if of2_is_edge {
                    let Some(oe2) = oe2.as_ref() else {
                        continue;
                    };
                    let cone = curve_on_edge(&cur_cut_e, oe2);
                    let c = make_curve(&cone);
                    let p1 = c.d0(c.first_parameter());
                    let p2 = c.d0(c.last_parameter());
                    let vf = find_vertex(&p1, &map_on_v1f, self.my_tol)
                        .or_else(|| find_vertex(&p1, &map_on_v1l, self.my_tol));
                    let vl = find_vertex(&p2, &map_on_v1l, self.my_tol)
                        .or_else(|| find_vertex(&p2, &map_on_v1f, self.my_tol));
                    let e = if p1.square_distance(&p2) < self.my_tol * self.my_tol {
                        make_degenerated_edge(&c, vf.as_ref())
                    } else {
                        let Some(edge) = make_edge_3d(&c, vf.as_ref(), vl.as_ref()) else {
                            continue;
                        };
                        edge
                    };
                    (e, vf, vl)
                } else {
                    let (Some(gs2), Some(of2)) = (gs2.as_ref(), of2.as_ref()) else {
                        continue;
                    };
                    let p2d = pc2.d0(f2p);
                    let p1 = gs2.d0(p2d.x(), p2d.y());
                    let p2d = pc2.d0(l2p);
                    let p2 = gs2.d0(p2d.x(), p2d.y());
                    let vf = find_vertex(&p1, &map_on_v1f, self.my_tol);
                    let vl = find_vertex(&p2, &map_on_v1l, self.my_tol);
                    let Some(e) = make_edge_2d(
                        &pc2,
                        gs2,
                        of2,
                        vf.as_ref(),
                        vl.as_ref(),
                        f2p,
                        l2p,
                    ) else {
                        continue;
                    };
                    kpart_curve3d(&e, &pc2, gs2);
                    (e, vf, vl)
                };

                // `cxx:2054-2094`: E1f / E1l fallback on the vertex sets.
                if e1f.is_none() {
                    if let (Some(vf1), Some(vf2)) = (vf_on_e1.as_ref(), vf_on_e2.as_ref()) {
                        e1f = join_on_vertices(&map_on_v1f, vf1, vf2);
                    }
                }
                if e1l.is_none() {
                    if let (Some(vl1), Some(vl2)) = (vl_on_e1.as_ref(), vl_on_e2.as_ref()) {
                        e1l = join_on_vertices(&map_on_v1l, vl1, vl2);
                    }
                }

                // `cxx:2096-2097`.
                e1.0.set_orientation(orien1);
                e2.0.set_orientation(orien2);

                // `cxx:2098-2099`: `BRepOffset_Offset(CurCutE, E1, E2, -myRadius,
                // E1f, E1l, myNubs, myTol, GeomAbs_C2)`. The two edges are
                // rewritten in place; `myCreated` records them afterwards so
                // that it sees that final state.
                let mut an_offset = BRepOffsetOffset::new();
                if an_offset
                    .init_pipe_first_last(
                        &cur_cut_e,
                        &mut e1,
                        &mut e2,
                        -self.my_radius,
                        e1f.as_ref(),
                        e1l.as_ref(),
                        self.my_nubs,
                        self.my_tol,
                        Shape::C2,
                    )
                    .is_err()
                {
                    continue;
                }

                // `cxx:2042-2049`.
                if self.my_build_shape {
                    self.my_created.bind_empty(&cur_cut_e.0);
                    self.my_created.inner_bind_empty(&cur_cut_e.0, &or1);
                    self.my_created.inner_append(&cur_cut_e.0, &or1, &e1.0);
                    self.my_created.inner_bind_empty(&cur_cut_e.0, &or2);
                    self.my_created.inner_append(&cur_cut_e.0, &or2, &e2.0);
                }

                // `cxx:2100-2103`.
                let Some(tuyo) = an_offset.face().cloned() else {
                    continue;
                };
                self.my_map_sf.bind(&cur_cut_e.0, an_offset);
                self.my_centers.add(&cur_cut_e.0);
                builder.add(&mut co, &cur_cut_e.0);
                if let Some(result) = self.my_result.as_mut() {
                    builder.add(result, &tuyo.0);
                }

                // `cxx:2108-2178`.
                if self.my_build_shape {
                    let tube_edges = explorer_shapes(&tuyo.0, ShapeType::Edge);
                    if of1_is_edge {
                        if let (Some(gc1), Some(eon_f1)) = (gc1.as_ref(), tube_edges.first()) {
                            update_cut_edges(
                                &Edge(eon_f1.clone()),
                                gc1,
                                &or1,
                                self.my_tol,
                                &mut self.my_cut_edges,
                            );
                        }
                    }
                    if of2_is_edge {
                        if let (Some(gc2), Some(eon_f2)) = (gc2.as_ref(), tube_edges.get(1)) {
                            update_cut_edges(
                                &Edge(eon_f2.clone()),
                                gc2,
                                &or2,
                                self.my_tol,
                                &mut self.my_cut_edges,
                            );
                        }
                    }
                }
            }
        }

        // `cxx:2184-2216`: spheres on vertices with three incident centre
        // edges.
        let map = AncestorMap::map_shapes_and_ancestors(&co, ShapeType::Vertex, ShapeType::Edge);
        for v in explorer_shapes(&co, ShapeType::Vertex) {
            let Some(anc) = map.find_from_key(&v) else {
                continue;
            };
            if anc.len() != 3 {
                continue;
            }
            let mut loe: Vec<Edge> = Vec::new();
            for ae in anc {
                let Some(of) = self.my_map_sf.get(ae) else {
                    continue;
                };
                let Some(gen) = of.generated(&v) else {
                    continue;
                };
                // `cxx:2201-2206`: `Reverse` is the constant `true`, so the
                // generated edge is reversed.
                loe.push(Edge(gen.oriented(Orientation::Reversed)));
            }
            if loe.len() != 3 {
                continue;
            }

            let vertex = Vertex(v.clone());
            let mut oft = BRepOffsetOffset::new();
            if oft
                .init_vertex(
                    &vertex,
                    &loe,
                    -self.my_radius,
                    self.my_nubs,
                    self.my_tol,
                    Shape::C2,
                )
                .is_err()
            {
                continue;
            }
            let Some(face) = oft.face().cloned() else {
                continue;
            };
            self.my_map_sf.bind(&v, oft);
            self.my_centers.add(&v);
            if let Some(result) = self.my_result.as_mut() {
                builder.add(result, &face.0);
            }
        }
    }

}

/// `cxx:2110-2172`: project the tube's end vertices onto the support curve and
/// register their parameters on the support edge, then collect the vertices.
fn update_cut_edges(
    eon_f: &Edge,
    support: &Arc<dyn Curve>,
    or: &TopoShape,
    tol: f64,
    cut_edges: &mut std::collections::HashMap<crate::topo_tools::ShapeMapKey, Vec<TopoShape>>,
) {
    let (v1, v2) = edge_vertices(eon_f);
    for v in [v1.clone(), v2.clone()].into_iter().flatten() {
        let p = BRepTool::vertex_point(&v);
        let Some(proj) = project_point_on_curve(support.as_ref(), &p, CONFUSION) else {
            continue;
        };
        let _ = update_vertex_on_edge(
            &v.0.oriented(Orientation::Internal),
            proj.parameter,
            or,
            tol,
        );
    }
    cut_edges
        .entry(shape_map_key(or))
        .or_default()
        .extend([v1, v2].into_iter().flatten().map(|v| v.0));
}

/// `cxx:2060-2072` / `:2078-2092`: an edge of `map` whose endpoints are the
/// two given vertices.
fn join_on_vertices(map: &EdgeSet, v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    for e in map.iter() {
        let (a, b) = edge_vertices(e);
        let (Some(a), Some(b)) = (a, b) else {
            continue;
        };
        if (is_same(&a.0, &v1.0) && is_same(&b.0, &v2.0))
            || (is_same(&b.0, &v1.0) && is_same(&a.0, &v2.0))
        {
            return Some(e.clone());
        }
    }
    None
}



