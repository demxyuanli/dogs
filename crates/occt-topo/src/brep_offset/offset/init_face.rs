//! `BRepOffset_Offset::Init(Face, Offset, Created, OffsetOutside, JoinType)`
//! (`BRepOffset_Offset.cxx:466-1094`), the general face-offset arm.
//!
//! The surface half is `singularity.rs` (`:500-807`); this file is
//! `:811-1093`: build `myFace` on the offset surface, seed `MapSS` with the
//! vertex images supplied through `Created` (`:822-870`), then rebuild every
//! wire (`:872-1089`).
//!
//! Every branch below is the `.cxx` branch with the same predicate: nothing
//! here is special-cased for a given STEP.

use std::collections::HashSet;
use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::CONFUSION;
use occt_geom::Surface;
use occt_geom2d::{Curve2d, Geom2dLine};

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_offset::brep_offset_surface;
use crate::brep_tool::BRepTool;
use crate::brep_tools;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    edge_vertices, edges_of_wire, explorer_shapes, shape_map_key, wires_of_face, ShapeMapKey,
};

use super::compute_curve3d::{compute_curve3d, update_edge_pcurve, update_edge_pcurves};
use super::singularity::process_singularities;
use super::{BRepOffsetOffset, MapSs};

impl BRepOffsetOffset {
    /// `BRepOffset_Offset::Init(Face, Offset, Created, OffsetOutside,
    /// JoinType)` (`BRepOffset_Offset.cxx:466-1094`).
    ///
    /// `created` is OCCT's `Created` map: the images of the source edges and
    /// vertices that an adjacent, G1-continuous face already produced, so the
    /// shared boundary is built once (`cxx:843-869`, `:959`).
    pub fn init_face(
        &mut self,
        face: &Face,
        offset: f64,
        created: &MapSs,
        offset_outside: bool,
        join_type: JoinType,
    ) -> Result<(), String> {
        let reg = GeometryRegistry::global();
        let builder = TopoBuilder::new();

        // `cxx:473-478`.
        self.my_shape = Some(face.0.clone());
        let mut my_offset = offset;
        if face.0.orientation() == Orientation::Reversed {
            my_offset = -my_offset;
        }

        // `cxx:480-499`: strip the trimming, then bake the location into the
        // surface for the types whose `Copy` + `Transform` OCCT trusts.
        let loc = face.0.location().clone();
        let mut s =
            BRepTool::face_surface(face).ok_or("BRepOffset_Offset::Init: face has no surface")?;
        if let Some(basis) = s.rectangular_trimmed_basis() {
            s = basis;
        }
        let mut is_transformed = false;
        if (s.is_bspline_surface()
            || s.is_surface_of_linear_extrusion()
            || s.is_surface_of_revolution()
            || s.is_offset_surface())
            && !loc.is_identity()
        {
            let t = loc.transformation();
            s = Arc::from(s.transformed(&t));
            is_transformed = true;
        }

        // `cxx:500-522`: the cone particular case. `vv1`/`vv2` are
        // `BRepTools::UVBounds(Face, ...)` (`cxx:509`), not the surface
        // bounds used further down.
        if let Some(cone) = s.gp_cone() {
            let apex = cone.apex();
            let (_, vc) =
                slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), &apex);
            let (_, _, vv1, vv2) = brep_tools::uv_bounds(face);
            if vv2 < vc && cone.semi_angle() > 0.0 {
                my_offset = -my_offset;
            } else if vv1 > vc && cone.semi_angle() < 0.0 {
                my_offset = -my_offset;
            }
            if !cone.position().is_direct() {
                my_offset = -my_offset;
            }
        }

        // `cxx:524`.
        let (the_box, status) = brep_offset_surface(s.as_ref(), my_offset, false)?;
        self.my_status = status;
        let mut the_surf: Arc<dyn Surface> = Arc::from(the_box);

        // `cxx:526-807`.
        let singular = process_singularities(face, &s, &mut the_surf, offset_outside, join_type)?;

        // `cxx:811-820`: `MakeFace` + `UpdateFace(myFace, TheSurf, L, tol)`.
        let mut my_face =
            builder.make_face_with_tolerance(the_surf.clone(), &[], BRepTool::face_tolerance(face));
        if !is_transformed {
            my_face.0.set_location(&loc);
        }

        // `cxx:822-831`.
        let cur_face = Face(face.0.oriented(Orientation::Forward));
        let (_, _, vf1, _) = brep_tools::uv_bounds(face);
        let (_, _, sur_v1, sur_v2) = (0.0, 0.0, the_surf.v_range().0, the_surf.v_range().1);

        // `cxx:829` `VonDegen`.
        let mut von_degen: HashSet<ShapeMapKey> = HashSet::new();
        let mut map_ss: MapSs = MapSs::new();

        // `cxx:833-870`: seed `MapSS` with the images `Created` already holds.
        for shape in explorer_shapes(&cur_face.0, ShapeType::Edge) {
            let e = Edge(shape);
            let (v1, v2) = edge_vertices(&e);
            if singular.has && BRepTool::is_degenerated(&e) {
                if let Some(v1) = &v1 {
                    von_degen.insert(shape_map_key(&v1.0));
                }
            }
            if let Some(oe) = created.get(&shape_map_key(&e.0)) {
                let oe = Edge(oe.clone());
                let (ov1, ov2) = edge_vertices(&oe);
                if let (Some(v1), Some(ov1)) = (&v1, &ov1) {
                    map_ss.entry(shape_map_key(&v1.0)).or_insert(ov1.0.clone());
                }
                if let (Some(v2), Some(ov2)) = (&v2, &ov2) {
                    map_ss.entry(shape_map_key(&v2.0)).or_insert(ov2.0.clone());
                }
            }
            if let Some(v1) = &v1 {
                if let Some(o) = created.get(&shape_map_key(&v1.0)) {
                    map_ss.entry(shape_map_key(&v1.0)).or_insert(o.clone());
                }
            }
            if let Some(v2) = &v2 {
                if let Some(o) = created.get(&shape_map_key(&v2.0)) {
                    map_ss.entry(shape_map_key(&v2.0)).or_insert(o.clone());
                }
            }
        }

        // `cxx:872-1089`.
        for w in wires_of_face(&cur_face) {
            let wire_fwd = Wire(w.0.oriented(Orientation::Forward));
            let mut new_edges: Vec<Edge> = Vec::with_capacity(4);
            for e in edges_of_wire(&wire_fwd) {
                let Some((v1, v2)) = endpoints(&e) else {
                    continue;
                };
                // `cxx:888`: OCCT dereferences the null handle when the edge
                // has no p-curve on this face.
                let Some((c2d, read_f, read_l)) = boptools_2d::curve_on_surface_range(&e, &cur_face)
                else {
                    continue;
                };
                let mut f = read_f;
                let mut l = read_l;
                let tol_e = BRepTool::edge_tolerance(&e);
                let k_e = shape_map_key(&e.0);
                let k1 = shape_map_key(&v1.0);
                let k2 = shape_map_key(&v2.0);

                // `cxx:890-908`: an edge already in `MapSS` whose endpoints
                // are not singular is a seam occurrence.
                let oe: Edge = if map_ss.contains_key(&k_e)
                    && !von_degen.contains(&k1)
                    && !von_degen.contains(&k2)
                {
                    let base = Edge(map_ss[&k_e].clone());
                    let e_rev = Edge(e.0.oriented(Orientation::Reversed));
                    let Some((c2d_1, rf, rl)) =
                        boptools_2d::curve_on_surface_range(&e_rev, &cur_face)
                    else {
                        continue;
                    };
                    f = rf;
                    l = rl;
                    if e.orientation() == Orientation::Forward {
                        update_edge_pcurves(&base, &c2d, &c2d_1, &my_face, tol_e);
                    } else {
                        update_edge_pcurves(&base, &c2d_1, &c2d, &my_face, tol_e);
                    }
                    reg.set_edge_range(&base.0, f, l);
                    base
                } else {
                    // `cxx:910-913`.
                    let e_fwd = Edge(e.0.oriented(Orientation::Forward));
                    let Some(t1) =
                        BRepTool::parameter_of_vertex_on_face_checked(&v1.0, &e_fwd, &cur_face)
                    else {
                        continue;
                    };
                    let Some(t2) =
                        BRepTool::parameter_of_vertex_on_face_checked(&v2.0, &e_fwd, &cur_face)
                    else {
                        continue;
                    };
                    let p2d1 = c2d.value(t1);
                    let p2d2 = c2d.value(t2);

                    // `cxx:914-957`.
                    let (p1, v_start) = if von_degen.contains(&k1) {
                        if (p2d1.y() - vf1).abs() <= CONFUSION {
                            (singular.min_apex, sur_v1)
                        } else {
                            (singular.max_apex, sur_v2)
                        }
                    } else {
                        let mut p1 = the_surf.d0(p2d1.x(), p2d1.y());
                        if !loc.is_identity() && !is_transformed {
                            p1.transform(&loc.transformation());
                        }
                        (p1, p2d1.y())
                    };
                    let (p2, v_end) = if von_degen.contains(&k2) {
                        if (p2d2.y() - vf1).abs() <= CONFUSION {
                            (singular.min_apex, sur_v1)
                        } else {
                            (singular.max_apex, sur_v2)
                        }
                    } else {
                        let mut p2 = the_surf.d0(p2d2.x(), p2d2.y());
                        if !loc.is_identity() && !is_transformed {
                            p2.transform(&loc.transformation());
                        }
                        (p2, p2d2.y())
                    };

                    // `cxx:958-1009`.
                    let base = if let Some(oe) = created.get(&k_e) {
                        Edge(oe.clone())
                    } else if let Some(oe) = map_ss.get(&k_e) {
                        Edge(oe.clone())
                    } else {
                        let mut oe = Edge(TopoShape::new(ShapeType::Edge));
                        if !is_transformed {
                            oe.0.set_location(&loc);
                        }
                        let ov1 = match map_ss.get(&k1) {
                            Some(ov) => Vertex(ov.clone()),
                            None => {
                                let ov =
                                    builder.make_vertex(p1, BRepTool::vertex_tolerance(&v1));
                                map_ss.insert(k1.clone(), ov.0.clone());
                                ov
                            }
                        };
                        let ov2 = match map_ss.get(&k2) {
                            Some(ov) => Vertex(ov.clone()),
                            None => {
                                let ov =
                                    builder.make_vertex(p2, BRepTool::vertex_tolerance(&v2));
                                map_ss.insert(k2.clone(), ov.0.clone());
                                ov
                            }
                        };
                        builder.add_edge_vertices(&mut oe, &ov1, &ov2);
                        if BRepTool::is_degenerated(&e) {
                            reg.set_degenerated(&oe.0, true);
                        }
                        oe
                    };

                    // `cxx:1010-1079`.
                    let mut c2d_final = c2d.clone();
                    if von_degen.contains(&k1) || von_degen.contains(&k2) {
                        let mut q1 = p2d1;
                        let mut q2 = p2d2;
                        if von_degen.contains(&k1) {
                            q1.set_y(v_start);
                        }
                        if von_degen.contains(&k2) {
                            q2.set_y(v_end);
                        }
                        c2d_final = Arc::new(two_point_line(q1, q2)?);
                        f = 0.0;
                        l = q1.distance(&q2);
                        if map_ss.contains_key(&k_e) {
                            let Some((existing, ef, el)) =
                                boptools_2d::curve_on_surface_range(&base, &my_face)
                            else {
                                continue;
                            };
                            f = ef;
                            l = el;
                            if e.orientation() == Orientation::Forward {
                                update_edge_pcurves(&base, &c2d_final, &existing, &my_face, tol_e);
                            } else {
                                update_edge_pcurves(&base, &existing, &c2d_final, &my_face, tol_e);
                            }
                        } else {
                            update_edge_pcurve(&base, &c2d_final, &my_face, tol_e);
                        }
                        // `cxx:1040`: `Range(OE, myFace, f, l)`.
                        let face_key = GeometryRegistry::shape_key(&my_face.0);
                        reg.set_pcurve_range(&base.0, face_key, f, l);

                        // `cxx:1041-1065`.
                        if !BRepTool::is_degenerated(&e) && the_surf.is_u_closed() {
                            let e_rev = Edge(e.0.oriented(Orientation::Reversed));
                            let Some((c2d_1, _, _)) =
                                boptools_2d::curve_on_surface_range(&e_rev, &cur_face)
                            else {
                                continue;
                            };
                            let Some(r1) =
                                BRepTool::parameter_of_vertex_on_face_checked(&v1.0, &e, &cur_face)
                            else {
                                continue;
                            };
                            let Some(r2) =
                                BRepTool::parameter_of_vertex_on_face_checked(&v2.0, &e, &cur_face)
                            else {
                                continue;
                            };
                            let mut q1 = c2d_1.value(r1);
                            let mut q2 = c2d_1.value(r2);
                            if von_degen.contains(&k1) {
                                q1.set_y(v_start);
                            }
                            if von_degen.contains(&k2) {
                                q2.set_y(v_end);
                            }
                            let second: Arc<dyn Curve2d> = Arc::new(two_point_line(q1, q2)?);
                            if e.orientation() == Orientation::Forward {
                                update_edge_pcurves(&base, &c2d_final, &second, &my_face, tol_e);
                            } else {
                                update_edge_pcurves(&base, &second, &c2d_final, &my_face, tol_e);
                            }
                        }
                    } else {
                        // `cxx:1075-1079`.
                        update_edge_pcurve(&base, &c2d_final, &my_face, tol_e);
                        reg.set_edge_range(&base.0, f, l);
                    }

                    // `cxx:1080-1083`.
                    if !BRepTool::is_degenerated(&base) {
                        compute_curve3d(&base, &c2d_final, &the_surf, tol_e);
                    }
                    map_ss.insert(k_e, base.0.clone());
                    base
                };

                // `cxx:1086`.
                new_edges.push(Edge(oe.0.oriented(e.orientation())));
            }
            // `cxx:1088`.
            let ow = builder.make_wire(&new_edges);
            builder.add_wire(&mut my_face, &Wire(ow.0.oriented(w.orientation())));
        }

        // `cxx:1091-1093`.
        my_face.0.set_orientation(face.0.orientation());
        self.my_face = Some(my_face);
        Ok(())
    }
}

/// `TopExp::Vertices(E, V1, V2)` (`cxx:884`) with both endpoints present.
fn endpoints(e: &Edge) -> Option<(Vertex, Vertex)> {
    let (v1, v2) = edge_vertices(e);
    Some((v1?, v2?))
}

/// `new Geom2d_Line(P2d1, gp_Vec2d(P2d1, P2d2))` (`cxx:1020`, `:1056`).
fn two_point_line(p1: GpPnt2d, p2: GpPnt2d) -> Result<Geom2dLine, String> {
    let v = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
    let dir = GpDir2d::from_vec2d(&v)
        .map_err(|_| "BRepOffset_Offset: null 2d line direction".to_string())?;
    Ok(Geom2dLine::from_pnt_dir(p1, dir))
}
