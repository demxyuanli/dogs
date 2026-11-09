//! `BRepOffset_MakeOffset::CorrectConicalFaces`
//! (`BRepOffset_MakeOffset.cxx:2396-2853`) and the file-local `GetEdgePoints`
//! (`:515-530`).
//!
//! The method repairs the conical faces the arc-join offset produced around an
//! apex: every degenerated edge that is **not** really degenerated (its three
//! sample points are not within `TolApex = 1e-5`, `:2435-2441`) gives the apex
//! vertex of its root edge (`:2445-2446`), and all the faces meeting at one
//! apex are collected in `FacesOfCone` (`:2447-2463`). Nothing else happens when
//! that map is empty (`:2757-2760`), which is the state of every offset that
//! produced no conical apex.
//!
//! For each apex the reference then rebuilds one spherical cap (`:2467-2760`):
//!
//! * `:2485-2527` walks the faces of the apex, picks the degenerated edge whose
//!   root starts at the apex (`:2489-2501`), turns it into a normal edge
//!   (`:2503-2506`) and, on the first face, builds the `gp_Sphere` centred at
//!   the apex with radius `myOffset` from the frame spanned by the edge's
//!   `fPnt / mPnt / lPnt` (`:2511-2527`).
//! * `:2531-2546` re-hangs the edge on the chain of vertices the cap shares
//!   (the first face's `fPnt` closes the loop) and `:2547-2573` puts the
//!   corresponding iso-circle of the sphere on it as its 3D curve plus a
//!   `Geom2d_Line` pcurve.
//! * `:2576-2622` moves the faces' other edges from the old apex vertex to the
//!   new chain.
//! * `:2625-2746` orders the wire, cuts the spherical face
//!   (`BRepLib_MakeFace(aSphSurf, Ufirst, Ulast, Vfirst, Vlast, Confusion)`) and
//!   substitutes the chain into it with `BRepTools_Substitution`.
//! * `:2748-2756` hangs the new face on the (single) shell of `myOffsetShape`.
//!
//! `:2762-2851` then repeats the solid/compound re-assembly of `MakeSolid`
//! (`:3818-3893`) whenever the loop modified the shape, gated by
//! `myShape.ShapeType() == TopAbs_SOLID || myThickening` (`:2766`).
//!
//! Port notes:
//!
//! * `TopExp::MapShapes(aFace, TopAbs_EDGE, Emap)` (`:2424`, `:2581`) answers an
//!   `IsSame`-unique map; [`explorer_shapes`] repeats a seam edge met twice, so
//!   the port de-duplicates by [`shape_map_key`] at both sites.
//! * `BRepAlgo_Image::Root(anEdge)` (`:2445`, `:2498`) throws
//!   `Standard_NoSuchObject` when the edge is not bound; the port skips the
//!   edge.
//! * `BRepAdaptor_Curve(anEdge, aFace)` (`:2434`) uses the pcurve of the edge on
//!   the face when there is one and the 3D curve otherwise; both are available
//!   here, so [`curve_value_on_face`] follows `BRepAdaptor_Curve.cxx:40-64`.
//! * `TopoDS::Vertex/Face/Edge` casts (`:2418`, `:2479`, `:2487`, `:2587`,
//!   `:2600`, `:2633`, `:2662`, `:2679`, `:2710`, `:2719`) throw
//!   `Standard_TypeMismatch` on a wrong type; the port skips the shape.
//! * `BRep_Tool::CurveOnSurface(E, S, L, f, l)` on a bare surface (`:2637`,
//!   `:2666`) has no face to key the lookup in OCCT; the port keys the pcurve
//!   map by the surface pointer, which is the equality
//!   `BRep_CurveRepresentation::IsCurveOnSurface(S, l)` tests anyway.
//! * `gp_Sphere::SetRadius(aRadius)` raises `Standard_ConstructionError` for a
//!   negative radius (`gp_Sphere.cxx:88-96`); the port skips the apex group.
//! * `RealLast()` / `RealFirst()` (`:2626`) are `f64::MAX` / `f64::MIN` here.
//! * `theShell.Closed(true)` after `BRep_Tool::IsClosed(theShell)`
//!   (`:2751-2755`); `BRep_Tool::IsClosed(Shell)` is
//!   `bop_build_solids::geometrically_open` negated.
//! * `Message_ProgressScope` does not appear in this method, so nothing is
//!   dropped for progress.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpAx2, GpDir, GpPnt, GpPnt2d, GpSphere, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::sphere::GeomSphere;
use occt_geom::surface::Surface;
use occt_geom::trimmed::GeomTrimmedCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_lib_make_face::make_face_uv;
use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_offset::BRepOffsetMakeOffset;
use crate::brep_tool::BRepTool;
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key, ShapeMapKey};

/// `TolApex` (`cxx:2404`).
const TOL_APEX: f64 = 1.0e-5;

/// `GetEdgePoints(anEdge, aFace, fPnt, mPnt, lPnt)` (`cxx:515-530`): the 3D
/// points of the edge's pcurve at its first, middle and last parameter.
fn get_edge_points(edge: &Edge, face: &Face) -> Option<(GpPnt, GpPnt, GpPnt)> {
    let (c2d, f, l) = curve_on_surface_range(edge, face)?;
    let surf = BRepTool::face_surface(face)?;
    let f_pnt2d = c2d.d0(f);
    let l_pnt2d = c2d.d0(l);
    let m_pnt2d = c2d.d0(0.5 * (f + l));
    Some((
        surf.d0(f_pnt2d.x(), f_pnt2d.y()),
        surf.d0(m_pnt2d.x(), m_pnt2d.y()),
        surf.d0(l_pnt2d.x(), l_pnt2d.y()),
    ))
}

/// `BRepAdaptor_Curve(E, F).Value(theParam)` (`BRepAdaptor_Curve.cxx:40-64`):
/// the pcurve of `E` on `F` lifted through the surface, the 3D curve when the
/// edge carries no pcurve on that face.
fn curve_value_on_face(edge: &Edge, face: &Face, param: f64) -> Option<GpPnt> {
    if let Some((c2d, _, _)) = curve_on_surface_range(edge, face) {
        let surf = BRepTool::face_surface(face)?;
        let uv = c2d.d0(param);
        return Some(surf.d0(uv.x(), uv.y()));
    }
    BRepTool::edge_curve(edge).map(|c| c.d0(param))
}

/// `BRep_Tool::CurveOnSurface(E, S, L, f, l)` (`BRep_Tool.cxx:167-180`) for a
/// bare `Geom_Surface`: the representation is keyed by the surface handle.
fn pcurve_on_surface(
    edge: &Edge,
    surf: &Arc<dyn Surface>,
) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let ptr = Arc::as_ptr(surf) as *const () as usize;
    GeometryRegistry::global()
        .edge_pcurve_reps(&edge.0)
        .into_iter()
        .find(|(s, _, _, _)| Arc::as_ptr(s) as *const () as usize == ptr)
        .map(|(_, c, f, l)| (c, f, l))
}

/// `BRep_Builder::Remove(theS, V)` for a vertex listed by `TopExp::Vertices`.
///
/// `TopoDS_Builder::Remove` matches `(TShape, orientation, location)`
/// (`TopoDS_Builder.cxx:107-133`) while `TopExp::Vertices` hands back the
/// vertex occurrence with its stored orientation; the port therefore removes
/// every vertex child sharing the `TShape`, which is the set the reference
/// call removes.
fn remove_vertex_occurrences(edge: &mut TopoShape, vertex: &TopoShape) {
    let builder = TopoBuilder::new();
    let children = edge.tshape.read().expect("poisoned TShape lock").children.clone();
    for c in children {
        if c.shape_type() == ShapeType::Vertex && c.same_tshape(vertex) {
            builder.remove(edge, &c);
        }
    }
}

/// `IsSame`-unique edge list of a face, the `TopExp::MapShapes(F, EDGE, M)`
/// walk of `cxx:2581`.
fn unique_edges_of_face(face: &Face) -> Vec<TopoShape> {
    let mut seen: Vec<ShapeMapKey> = Vec::new();
    let mut out: Vec<TopoShape> = Vec::new();
    for e in explorer_shapes(&face.0, ShapeType::Edge) {
        let k = shape_map_key(&e);
        if !seen.contains(&k) {
            seen.push(k);
            out.push(e);
        }
    }
    out
}

impl BRepOffsetMakeOffset {
    /// `CorrectConicalFaces()` (`cxx:2396-2853`).
    pub fn correct_conical_faces(&mut self) {
        let Some(offset_shape) = self.my_offset_shape.clone() else {
            return;
        };
        let builder = TopoBuilder::new();
        let reg = GeometryRegistry::global();

        // `cxx:2404`: `double TolApex = 1.e-5`.

        // `cxx:2410-2411`: `FacesOfCone`, an `IndexedDataMap` keyed by the apex
        // vertex. The port keeps the first-seen order instead of the reference's
        // hash order (only the sequence of the rebuilt caps can differ).
        let mut faces_of_cone: Vec<(ShapeMapKey, TopoShape, Vec<TopoShape>)> = Vec::new();

        // `cxx:2415-2466`.
        if self.my_join == occt_core::kernel::geomabs::JoinType::Arc {
            for face_shape in explorer_shapes(&offset_shape, ShapeType::Face) {
                let Some(a_face) = Face::wrap(face_shape) else {
                    continue;
                };
                for edge_shape in unique_edges_of_face(&a_face) {
                    let Some(an_edge) = Edge::wrap(edge_shape) else {
                        continue;
                    };
                    if !BRepTool::is_degenerated(&an_edge) {
                        continue;
                    }
                    // `cxx:2433-2441`: a really degenerated edge is skipped.
                    let Some((p_first, p_mid, p_last)) = get_edge_points(&an_edge, &a_face)
                    else {
                        continue;
                    };
                    if p_first.distance(&p_last) <= TOL_APEX
                        && p_first.distance(&p_mid) <= TOL_APEX
                    {
                        continue;
                    }
                    // `cxx:2445-2446`: the apex is the first vertex of the root
                    // edge of the offset edge.
                    let Ok(or_edge_shape) = self.my_init_offset_edge.root(&an_edge.0) else {
                        continue;
                    };
                    let Some(or_edge) = Edge::wrap(or_edge_shape) else {
                        continue;
                    };
                    let (v_f, _) = edge_vertices(&or_edge);
                    let Some(v_f) = v_f else {
                        continue;
                    };
                    // `cxx:2447-2463`.
                    let key = shape_map_key(&v_f.0);
                    if let Some(entry) = faces_of_cone.iter_mut().find(|e| e.0 == key) {
                        entry.2.push(a_face.0.clone());
                    } else {
                        faces_of_cone.push((key, v_f.0.clone(), vec![a_face.0.clone()]));
                    }
                }
            }
        }

        // `cxx:2757-2760`: nothing collected, nothing modified.
        if faces_of_cone.is_empty() {
            return;
        }

        // `cxx:2468-2473`.
        let mut is_modified = false;

        // `cxx:2474-2759`: one spherical cap per apex.
        for (apex_key, apex_shape, cone_faces) in &faces_of_cone {
            let Some(an_apex) = Vertex::wrap(apex_shape.clone()) else {
                continue;
            };
            let apex_pnt = BRepTool::vertex_point(&an_apex);

            // `cxx:2475-2479`.
            let mut the_sphere: Option<GpSphere> = None;
            let mut a_sph_surf: Option<Arc<dyn Surface>> = None;
            let mut sphere_wire: Vec<TopoShape> = Vec::new();
            let mut is_first_face = true;
            let mut first_point: Option<GpPnt> = None;
            let mut the_first_vertex: Option<TopoShape> = None;
            let mut cur_first_vertex: Option<TopoShape> = None;

            // `cxx:2483-2624`.
            for face_shape in cone_faces {
                let Some(a_face) = Face::wrap(face_shape.clone()) else {
                    continue;
                };
                // `cxx:2487-2502`: the degenerated edge of this face whose root
                // starts at the apex.
                let mut deg_edge: Option<TopoShape> = None;
                for e_shape in explorer_shapes(&a_face.0, ShapeType::Edge) {
                    let Some(e) = Edge::wrap(e_shape) else {
                        continue;
                    };
                    if !BRepTool::is_degenerated(&e) {
                        continue;
                    }
                    let Ok(root) = self.my_init_offset_edge.root(&e.0) else {
                        continue;
                    };
                    let Some(root_e) = Edge::wrap(root) else {
                        continue;
                    };
                    let (v_f, _) = edge_vertices(&root_e);
                    if v_f.as_ref().map(|v| shape_map_key(&v.0)) == Some(apex_key.clone()) {
                        deg_edge = Some(e.0);
                        break;
                    }
                }
                let Some(deg_edge) = deg_edge else {
                    continue;
                };

                // `cxx:2503-2506`.
                let mut cur_edge = deg_edge.oriented(Orientation::Forward);
                reg.set_degenerated(&cur_edge, false);
                reg.set_same_range(&cur_edge, false);
                reg.set_same_parameter(&cur_edge, false);

                // `cxx:2508-2510`.
                let Some((f_pnt, m_pnt, l_pnt)) =
                    get_edge_points(&Edge(cur_edge.clone()), &a_face)
                else {
                    continue;
                };
                let (_f_par, _l_par) = BRepTool::edge_parameters(&Edge(cur_edge.clone()));

                // `cxx:2511-2528`: the sphere, built on the first face only.
                if is_first_face {
                    let a_vec1 = GpVec::from_pnts(&f_pnt, &m_pnt);
                    let a_vec2 = GpVec::from_pnts(&f_pnt, &l_pnt);
                    let a_norm = a_vec1.crossed(&a_vec2);
                    let apex_to_fpnt = GpVec::from_pnts(&apex_pnt, &f_pnt);
                    let y_dir = a_norm.crossed(&apex_to_fpnt);
                    let x_dir = y_dir.crossed(&a_norm);
                    let (Ok(nz), Ok(nx)) = (GpDir::from_vec(&a_norm), GpDir::from_vec(&x_dir))
                    else {
                        // `gp_Dir` raises `Standard_ConstructionError` on a
                        // null vector.
                        continue;
                    };
                    let Ok(an_ax2) = GpAx2::new(apex_pnt, nz, nx) else {
                        continue;
                    };
                    let Ok(base) = GpSphere::new(an_ax2.to_ax3(), 1.0) else {
                        continue;
                    };
                    let mut sphere = base;
                    // `gp_Sphere::SetRadius` raises on a negative radius.
                    if sphere.set_radius(self.my_offset).is_err() {
                        continue;
                    }
                    a_sph_surf = Some(Arc::new(GeomSphere::new(sphere)));
                    the_sphere = Some(sphere);
                    first_point = Some(f_pnt);
                    let v = builder.make_vertex(f_pnt, CONFUSION);
                    the_first_vertex = Some(v.0.clone());
                    cur_first_vertex = Some(v.0);
                }
                let (Some(sphere), Some(surf)) = (the_sphere, a_sph_surf.clone()) else {
                    continue;
                };
                let (Some(first_pnt), Some(first_vert), Some(cur_first)) =
                    (first_point, the_first_vertex.clone(), cur_first_vertex.clone())
                else {
                    continue;
                };

                // `cxx:2531-2546`: re-hang the edge on the cap vertex chain.
                let (v1, v2) = edge_vertices(&Edge(cur_edge.clone()));
                let end_vert = if l_pnt.distance(&first_pnt) <= CONFUSION {
                    first_vert.clone()
                } else {
                    builder.make_vertex(l_pnt, CONFUSION).0
                };
                cur_edge.set_free(true);
                if let Some(v1) = &v1 {
                    remove_vertex_occurrences(&mut cur_edge, &v1.0);
                }
                if let Some(v2) = &v2 {
                    remove_vertex_occurrences(&mut cur_edge, &v2.0);
                }
                builder.add(&mut cur_edge, &cur_first.oriented(Orientation::Forward));
                builder.add(&mut cur_edge, &end_vert.oriented(Orientation::Reversed));

                // `cxx:2547-2573`: the iso-circle of the sphere as the 3D curve.
                let (u_f, v_f) = slib::sphere_parameters(sphere.position(), &f_pnt);
                let (mut u_l, _v_l) = slib::sphere_parameters(sphere.position(), &l_pnt);
                if u_l.abs() <= CONFUSION {
                    u_l = 2.0 * std::f64::consts::PI;
                }
                let Some(basis) = surf.v_iso_curve(v_f) else {
                    continue;
                };
                let trim_curv = GeomTrimmedCurve::new(basis, u_f, u_l);
                reg.set_edge_curve3d(&cur_edge, Arc::new(trim_curv), CONFUSION);
                reg.set_edge_range(&cur_edge, u_f, u_l);
                let lin2d = Geom2dLine::from_pnt_dir(
                    GpPnt2d::new(0.0, v_f),
                    occt_core::gp::dir2d::GpDir2d::new(1.0, 0.0).expect("DX2d"),
                );
                let trim_lin2d = Geom2dTrimmedCurve::new(Arc::new(lin2d), u_f, u_l);
                reg.update_edge_pcurve_on_surface(
                    &cur_edge,
                    &surf,
                    Some(Arc::new(trim_lin2d)),
                    CONFUSION,
                );
                reg.set_pcurve_range_on_surface(&cur_edge, &surf, u_f, u_l);
                // `BRepLib::SameParameter(theEdge)` (`BRepLib.cxx:1237-1246`)
                // with the default `Precision::Confusion()`.
                same_parameter_edge_inplace(&Edge(cur_edge.clone()), CONFUSION);
                sphere_wire.push(cur_edge.clone());

                // `cxx:2576-2579`: the 2D ends of the edge on the face.
                let Some((ba_c2d, ba_f, ba_l)) =
                    curve_on_surface_range(&Edge(cur_edge.clone()), &a_face)
                else {
                    continue;
                };
                let f_pnt2d = ba_c2d.d0(ba_f);
                let l_pnt2d = ba_c2d.d0(ba_l);
                let _ = &l_pnt2d;

                // `cxx:2581-2622`: move the face's edges off the old vertex.
                let mut ee: Vec<TopoShape> = Vec::new();
                for a_shape in unique_edges_of_face(&a_face) {
                    let Some(an_edge) = Edge::wrap(a_shape) else {
                        continue;
                    };
                    if BRepTool::is_degenerated(&an_edge) {
                        continue;
                    }
                    let (v1e, v2e) = edge_vertices(&an_edge);
                    let hits = |v: &Option<Vertex>| {
                        v.as_ref()
                            .map(|v| v1.as_ref().map(|a| a.0.same_tshape(&v.0)).unwrap_or(false))
                            .unwrap_or(false)
                    };
                    if hits(&v1e) || hits(&v2e) {
                        ee.push(an_edge.0);
                    }
                }
                for e_shape in ee {
                    let mut e_fwd = e_shape.oriented(Orientation::Forward);
                    e_fwd.set_free(true);
                    let (v1e, v2e) = edge_vertices(&Edge(e_fwd.clone()));
                    let Some((ee_c2d, ee_f, ee_l)) =
                        curve_on_surface_range(&Edge(e_fwd.clone()), &a_face)
                    else {
                        continue;
                    };
                    let p2d1 = ee_c2d.d0(ee_f);
                    let p2d2 = ee_c2d.d0(ee_l);
                    let is_first = v1e
                        .as_ref()
                        .map(|v| v1.as_ref().map(|a| a.0.same_tshape(&v.0)).unwrap_or(false))
                        .unwrap_or(false);
                    if is_first {
                        let new_v = if p2d1.distance(&f_pnt2d) <= CONFUSION {
                            &first_vert
                        } else {
                            &end_vert
                        };
                        if let Some(v1e) = &v1e {
                            remove_vertex_occurrences(&mut e_fwd, &v1e.0);
                        }
                        builder.add(&mut e_fwd, &new_v.oriented(Orientation::Forward));
                    } else {
                        let new_v = if p2d2.distance(&f_pnt2d) <= CONFUSION {
                            &first_vert
                        } else {
                            &end_vert
                        };
                        if let Some(v2e) = &v2e {
                            remove_vertex_occurrences(&mut e_fwd, &v2e.0);
                        }
                        builder.add(&mut e_fwd, &new_v.oriented(Orientation::Reversed));
                    }
                }

                // `cxx:2623-2624`.
                is_first_face = false;
                cur_first_vertex = Some(end_vert);
            }

            // `cxx:2625-2656`: the U extent of the cap.
            let Some(surf) = a_sph_surf.clone() else {
                continue;
            };
            let mut u_first = f64::MAX;
            let mut u_last = f64::MIN;
            for e_shape in &sphere_wire {
                let Some(e) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                let Some((c2d, f, l)) = pcurve_on_surface(&e, &surf) else {
                    continue;
                };
                let p2d1 = c2d.d0(f);
                let p2d2 = c2d.d0(l);
                for x in [p2d1.x(), p2d2.x()] {
                    if x < u_first {
                        u_first = x;
                    }
                    if x > u_last {
                        u_last = x;
                    }
                }
            }
            if u_first > u_last {
                continue;
            }

            // `cxx:2657-2692`: the wire, starting at the edge whose first 2D
            // point is `u_first`.
            let mut edges_of_wire: Vec<TopoShape> = sphere_wire.clone();
            let mut first_edge: Option<TopoShape> = None;
            let mut p2d1_of_first = GpPnt2d::new(0.0, 0.0);
            let mut i = 0;
            while i < edges_of_wire.len() {
                let Some(e) = Edge::wrap(edges_of_wire[i].clone()) else {
                    i += 1;
                    continue;
                };
                let Some((c2d, f, _l)) = pcurve_on_surface(&e, &surf) else {
                    i += 1;
                    continue;
                };
                let p2d1 = c2d.d0(f);
                first_edge = Some(e.0);
                p2d1_of_first = p2d1;
                if (p2d1.x() - u_first).abs() <= CONFUSION {
                    edges_of_wire.remove(i);
                    break;
                }
                i += 1;
            }
            let Some(first_edge) = first_edge else {
                continue;
            };
            let mut new_edges: Vec<TopoShape> = Vec::new();
            new_edges.push(first_edge.oriented(Orientation::Forward));
            let (vf1, cur_vertex) = edge_vertices(&Edge(first_edge.clone()));
            let Some(vf1) = vf1 else {
                continue;
            };
            let mut cur_vertex = match cur_vertex {
                Some(v) => v,
                None => continue,
            };
            let mut i = 0;
            while i < edges_of_wire.len() {
                let Some(an_edge) = Edge::wrap(edges_of_wire[i].clone()) else {
                    i += 1;
                    continue;
                };
                let (v1e, v2e) = edge_vertices(&an_edge);
                let hit1 = v1e
                    .as_ref()
                    .map(|v| v.0.same_tshape(&cur_vertex.0))
                    .unwrap_or(false);
                let hit2 = v2e
                    .as_ref()
                    .map(|v| v.0.same_tshape(&cur_vertex.0))
                    .unwrap_or(false);
                if hit1 || hit2 {
                    new_edges.push(an_edge.0.oriented(Orientation::Forward));
                    let next = if hit1 { v2e } else { v1e };
                    let Some(next) = next else {
                        break;
                    };
                    cur_vertex = next;
                    edges_of_wire.remove(i);
                } else {
                    i += 1;
                }
            }

            // `cxx:2694-2707`: the V extent and the cap face.
            let (v_first, v_last) = if p2d1_of_first.y() > 0.0 {
                (p2d1_of_first.y(), std::f64::consts::FRAC_PI_2)
            } else {
                (-std::f64::consts::FRAC_PI_2, p2d1_of_first.y())
            };
            let mut new_spherical_face =
                make_face_uv(surf.clone(), u_first, u_last, v_first, v_last, CONFUSION);

            // `cxx:2708-2728`: the degenerated vertex of the cap and the edge
            // that does not touch it.
            let cap_edges = explorer_shapes(&new_spherical_face.0, ShapeType::Edge);
            let mut deg_vertex: Option<TopoShape> = None;
            for e_shape in &cap_edges {
                let Some(e) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                if BRepTool::is_degenerated(&e) {
                    let (v, _) = edge_vertices(&e);
                    deg_vertex = v.map(|v| v.0);
                    break;
                }
            }
            let Some(deg_vertex) = deg_vertex else {
                continue;
            };
            let mut old_edge: Option<TopoShape> = None;
            for e_shape in &cap_edges {
                let Some(e) = Edge::wrap(e_shape.clone()) else {
                    continue;
                };
                let (v1e, v2e) = edge_vertices(&e);
                let touches = |v: &Option<Vertex>| {
                    v.as_ref()
                        .map(|v| v.0.same_tshape(&deg_vertex))
                        .unwrap_or(false)
                };
                if !touches(&v1e) && !touches(&v2e) {
                    old_edge = Some(e.0);
                    break;
                }
            }
            let Some(old_edge) = old_edge else {
                continue;
            };
            let Some(old_edge_w) = Edge::wrap(old_edge.clone()) else {
                continue;
            };
            let (ov1, ov2) = edge_vertices(&old_edge_w);
            let (Some(ov1), Some(ov2)) = (ov1, ov2) else {
                continue;
            };

            // `cxx:2729-2746`: substitution of the chain into the cap.
            let lv1: Vec<TopoShape> = vec![vf1.oriented(Orientation::Forward)];
            let lv2: Vec<TopoShape> = vec![cur_vertex.0.oriented(Orientation::Forward)];
            let mut substitutor = BRepToolsSubstitution::new();
            substitutor.substitute(&ov1.0.oriented(Orientation::Forward), &lv1);
            if !ov1.0.same_tshape(&ov2.0) {
                substitutor.substitute(&ov2.0.oriented(Orientation::Forward), &lv2);
            }
            substitutor.substitute(&old_edge.oriented(Orientation::Forward), &new_edges);
            substitutor.build(&new_spherical_face.0);
            if substitutor.is_copied(&new_spherical_face.0) {
                if let Some(list) = substitutor.copy(&new_spherical_face.0) {
                    if let Some(first) = list.first() {
                        if let Some(f) = Face::wrap(first.clone()) {
                            new_spherical_face = f;
                        }
                    }
                }
            }

            // `cxx:2748-2756`: hang the cap on the shell.
            let shells = explorer_shapes(&offset_shape, ShapeType::Shell);
            let Some(mut the_shell) = shells.first().and_then(|s| Shell::wrap(s.clone())) else {
                continue;
            };
            the_shell.0.set_free(true);
            builder.add(&mut the_shell.0, &new_spherical_face.0);
            is_modified = true;
            if !the_shell.0.is_closed()
                && !crate::bop_build_solids::geometrically_open(&the_shell.0)
            {
                the_shell.0.set_closed(true);
            }
        }

        // `cxx:2761-2764`.
        if !is_modified {
            return;
        }

        // `cxx:2766-2851`: the solid / compound re-assembly, the same tail as
        // `MakeSolid` (`:3818-3893`).
        let is_solid = self
            .my_shape
            .as_ref()
            .map(|s| s.shape_type() == ShapeType::Solid)
            .unwrap_or(false);
        if !is_solid && !self.my_thickening {
            return;
        }
        let Some(offset_shape) = self.my_offset_shape.clone() else {
            return;
        };
        let mut nb_shell = 0;
        let mut nc: Vec<TopoShape> = Vec::new();
        let mut s1: Option<TopoShape> = None;
        let mut sol_shells: Vec<Shell> = Vec::new();
        for sh_shape in explorer_shapes(&offset_shape, ShapeType::Shell) {
            nb_shell += 1;
            if sh_shape.is_closed() {
                sol_shells.push(Shell(sh_shape));
            } else {
                nc.push(sh_shape.clone());
                if nb_shell == 1 {
                    s1 = Some(sh_shape);
                }
            }
        }
        let nbs = sol_shells.len();
        let sol = builder.make_solid(&sol_shells);
        sol.0.set_closed(true);
        let mut sol_is_null = nbs == 0;
        // `cxx:2806-2821`: the port's `BRepCheck_Analyzer` is a stub (see
        // `solids.rs`), so this arm never fires; it is kept for the shape of the
        // control flow.
        if nbs > 1 && !crate::brep_check::is_valid(&sol.0) {
            sol_is_null = true;
        }
        let nc_is_null = nc.is_empty();
        // `cxx:2825-2850`.
        if !sol_is_null && !nc_is_null {
            nc.push(sol.0);
            let comp = builder.make_compound_of(&nc);
            self.my_offset_shape = Some(comp.0);
        } else if sol_is_null && !nc_is_null {
            self.my_offset_shape = if nb_shell == 1 {
                s1
            } else {
                let comp = builder.make_compound_of(&nc);
                Some(comp.0)
            };
        } else if !sol_is_null && nc_is_null {
            self.my_offset_shape = Some(sol.0);
        } else {
            let comp = builder.make_compound_of(&nc);
            self.my_offset_shape = Some(comp.0);
        }
    }
}
