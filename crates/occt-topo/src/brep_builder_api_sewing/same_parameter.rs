//! `SameRange` / `SameParameter` and both `SameParameterEdge` overloads of
//! `BRepBuilderAPI_Sewing` (`cxx:109-139`, `:304-321`, `:335-456`,
//! `:600-1063`), plus the two static helpers they share:
//! `findNMVertices` (`cxx:517-574`) and `ComputeToleranceVertex`
//! (`cxx:576-660`).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};
use occt_geom::extrema_pc::extrema_ext_pc_range_done;
use occt_geom::geom_lib_same_range_2d;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

use super::closed_surface::pcurve_on_surface;
use super::shape_maps::{ShapeKeyedMap, ShapeSet};
use super::Sewing;

/// `Geom2d_Curve::Copy()`.
fn copy_curve2d(c: &Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
    Arc::from(c.clone_dyn())
}

/// `BRep_Tool::Range(edge)`.
fn edge_range(e: &TopoShape) -> (f64, f64) {
    GeometryRegistry::global().edge_parameters(e)
}

/// `TopExp::Vertices(edge, v1, v2)`.
fn edge_vertices(e: &TopoShape) -> (Option<TopoShape>, Option<TopoShape>) {
    match crate::topo_tools::edge_vertices(&Edge(e.clone())) {
        (Some(a), Some(b)) => (Some(a.0), Some(b.0)),
        (a, b) => (a.map(|v| v.0), b.map(|v| v.0)),
    }
}

/// `BRep_Tool::CurveOnSurface(edge, surface, loc, first, last)`.
fn pcurve_range_on_surface(
    edge: &TopoShape,
    surf: &Arc<dyn occt_geom::Surface>,
) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    pcurve_on_surface(edge, surf)
}

/// `findNMVertices(edge, seqVert, seqPars)` (`cxx:517-574`): the INTERNAL /
/// EXTERNAL vertices of the edge with their parameter on its curve.
pub(crate) fn find_nm_vertices(edge: &TopoShape) -> (Vec<TopoShape>, Vec<f64>) {
    let mut verts: Vec<TopoShape> = Vec::new();
    let children = edge.tshape.read().expect("poisoned TShape lock").children.clone();
    for child in &children {
        if matches!(
            child.orientation(),
            Orientation::Internal | Orientation::External
        ) {
            verts.push(child.clone());
        }
    }
    if verts.is_empty() {
        return (verts, Vec::new());
    }
    let Some(c3d) = GeometryRegistry::global().edge_curve(edge) else {
        return (verts, Vec::new());
    };
    let (first, last) = edge_range(edge);
    let pfirst = c3d.d0(first);
    let plast = c3d.d0(last);
    let mut pars = Vec::with_capacity(verts.len());
    for v in &verts {
        let pt = BRepTool::vertex_point_world(&Vertex(v.clone()));
        let dist_f2 = pfirst.square_distance(&pt);
        let dist_l2 = plast.square_distance(&pt);
        let mut apar = if dist_f2 > dist_l2 { last } else { first };
        let (sols, done) = extrema_ext_pc_range_done(c3d.as_ref(), &pt, first, last);
        if done && !sols.is_empty() {
            let mut dist2_min = dist_f2.min(dist_l2);
            let mut ind_min = 0usize;
            for (ind, s) in sols.iter().enumerate() {
                if s.sq_dist < dist2_min {
                    ind_min = ind;
                    dist2_min = s.sq_dist;
                }
            }
            if ind_min > 0 {
                apar = sols[ind_min].u;
            }
            pars.push(apar);
        }
    }
    (verts, pars)
}

/// `ComputeToleranceVertex(V1, V2, newV)` (`cxx:515-560`).
fn compute_tolerance_vertex2(v1: &TopoShape, v2: &TopoShape) -> Vertex {
    let a_p = [
        BRepTool::vertex_point_world(&Vertex(v1.clone())),
        BRepTool::vertex_point_world(&Vertex(v2.clone())),
    ];
    let a_r = [
        BRepTool::vertex_tolerance(&Vertex(v1.clone())),
        BRepTool::vertex_tolerance(&Vertex(v2.clone())),
    ];
    let (m, n) = if a_r[0] < a_r[1] { (1, 0) } else { (0, 1) };
    let d_r = a_r[m] - a_r[n];
    let a_vd = GpVec::from_pnts(&a_p[m], &a_p[n]);
    let a_d = a_vd.magnitude();
    let builder = TopoBuilder::new();
    if a_d <= d_r || a_d < Sewing::EPS {
        return builder.make_vertex(a_p[m], a_r[m]);
    }
    let a_rr = 0.5 * (a_r[m] + a_r[n] + a_d);
    let xyz_m = a_p[m].xyz();
    let xyz_n = a_p[n].xyz();
    let vd = a_vd.xyz();
    let k = d_r / a_d;
    let x = 0.5 * (xyz_m.x + xyz_n.x - vd.x * k);
    let y = 0.5 * (xyz_m.y + xyz_n.y - vd.y * k);
    let z = 0.5 * (xyz_m.z + xyz_n.z - vd.z * k);
    builder.make_vertex(GpPnt::new(x, y, z), a_rr)
}

/// `ComputeToleranceVertex(V1, V2, V3, newV)` (`cxx:562-598`).
fn compute_tolerance_vertex3(v1: &TopoShape, v2: &TopoShape, v3: &TopoShape) -> Vertex {
    let vv = [v1, v2, v3];
    let mut a_p = [GpPnt::zero(); 3];
    let mut a_r = [0.0f64; 3];
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for i in 0..3 {
        a_p[i] = BRepTool::vertex_point_world(&Vertex(vv[i].clone()));
        a_r[i] = BRepTool::vertex_tolerance(&Vertex(vv[i].clone()));
        x += a_p[i].x();
        y += a_p[i].y();
        z += a_p[i].z();
    }
    let center = GpPnt::new(x / 3.0, y / 3.0, z / 3.0);
    let mut a_dmax = -1.0f64;
    for i in 0..3 {
        let a_di = center.distance(&a_p[i]) + a_r[i];
        if a_di > a_dmax {
            a_dmax = a_di;
        }
    }
    TopoBuilder::new().make_vertex(center, a_dmax)
}

impl Sewing {
    /// `SameRange(CurvePtr, FirstOnCurve, LastOnCurve, RequestedFirst,
    /// RequestedLast)` (`cxx:109-139`): `GeomLib::SameRange` wrapped in the
    /// `try`/`catch` arm that turns a failure into a null curve.
    pub(crate) fn same_range(
        &self,
        curve: &Arc<dyn Curve2d>,
        first_on_curve: f64,
        last_on_curve: f64,
        requested_first: f64,
        requested_last: f64,
    ) -> Option<Arc<dyn Curve2d>> {
        geom_lib_same_range_2d::same_range_2d(
            PCONFUSION,
            curve,
            first_on_curve,
            last_on_curve,
            requested_first,
            requested_last,
        )
        .ok()
    }

    /// `SameParameter(edge)` (`cxx:307-321`): `BRepLib::SameParameter(edge)`
    /// with the exception arm folded into "leave the edge alone".
    pub(crate) fn same_parameter_edge_geom(&self, edge: &Edge) {
        let tol = BRepTool::edge_tolerance(edge);
        crate::brep_lib_same_parameter::same_parameter_edge_inplace(edge, tol);
    }

    /// `SameParameterEdge(edge, seqEdges, seqForward, mapMerged, locReShape)`
    /// (`cxx:357-513`): merge the sequence of sections onto one edge, keeping
    /// the 3D curve / first pcurve / parametrisation of the reference section.
    pub(crate) fn same_parameter_edge_seq(
        &mut self,
        edge: &TopoShape,
        seq_edges: &[TopoShape],
        seq_forward: &[bool],
        merged_faces: &mut ShapeSet,
        map_merged: &mut ShapeSet,
        loc_reshape: &mut crate::shape_fix_compose_shell::MapReShape,
        loc_is_main: bool,
    ) -> Option<Edge> {
        // Retrieve reference section (`cxx:341-346`).
        let mut edge1 = self.reshape.apply_option(edge)?;
        if !loc_is_main {
            // `cxx:343-346`: the local reshaper is applied only when it is a
            // different context, otherwise `Edge1` keeps the first result.
            edge1 = loc_reshape.apply_option(&edge1)?;
        }
        let mut edge1 = Edge(edge1).clone();
        let mut is_forward = true;
        let mut list_faces1: Vec<TopoShape> = Vec::new();

        if self.sewing {
            // Fill listFaces1 with the faces of Edge1 (`cxx:352-368`).
            let mut bnd1 = edge.clone();
            if let Some(b) = self.section_bound.get(edge) {
                bnd1 = b;
            }
            if let Some(faces) = self.bound_faces.get(&bnd1) {
                for f in faces {
                    if merged_faces.add(&f) {
                        list_faces1.push(f);
                    }
                }
            }
        } else {
            // Create presentation edge (`cxx:369-389`).
            let (v1, v2) = edge_vertices(&edge1.0);
            let v1 = v1
                .map(|v| self.vertex_node.get(&v).unwrap_or(v))
                .unwrap_or_else(|| TopoShape::new(ShapeType::Vertex));
            let v2 = v2
                .map(|v| self.vertex_node.get(&v).unwrap_or(v))
                .unwrap_or_else(|| TopoShape::new(ShapeType::Vertex));
            let new_edge = crate::shape_fix_compose_shell::empty_copied_edge(&edge1.0);
            let builder = TopoBuilder::new();
            let mut an_edge = new_edge.oriented(Orientation::Forward);
            builder.add(&mut an_edge, &v1.oriented(Orientation::Forward));
            builder.add(&mut an_edge, &v2.oriented(Orientation::Reversed));
            edge1 = Edge(new_edge);
        }

        let mut is_done = false;
        for (i, oedge2) in seq_edges.iter().enumerate() {
            if self.sewing {
                let mut edge2 = Edge(self.reshape.apply_option(oedge2)?);
                if !loc_is_main {
                    edge2 = Edge(loc_reshape.apply_option(&edge2.0)?);
                }
                // Calculate relative orientation (`cxx:409-411`).
                let mut orientation = seq_forward.get(i).copied().unwrap_or(true);
                if !is_forward {
                    orientation = !orientation;
                }
                // Retrieve the faces of the second edge (`cxx:413-420`).
                let mut bnd2 = oedge2.clone();
                if let Some(b) = self.section_bound.get(oedge2) {
                    bnd2 = b;
                }
                let Some(list_faces2) = self.bound_faces.get(&bnd2) else {
                    continue; // Skip floating edge.
                };
                let mut which_sec = 1;
                let Some(new_edge) = self.same_parameter_edge_pair(
                    &edge1.0,
                    &edge2.0,
                    &list_faces1,
                    &list_faces2,
                    orientation,
                    &mut which_sec,
                    true,
                ) else {
                    continue;
                };
                // Record the faces of the temporary merged edge (`cxx:428-432`).
                for f in &list_faces2 {
                    if merged_faces.add(f) {
                        list_faces1.push(f.clone());
                    }
                }
                if !orientation && which_sec != 1 {
                    is_forward = !is_forward;
                }
                edge1 = new_edge;
            }
            map_merged.add(oedge2);
            is_done = true;
            if !self.nonmanifold {
                break;
            }
        }

        if is_done {
            edge1
                .0
                .set_orientation(if is_forward {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                });
            Some(edge1)
        } else {
            None
        }
    }

    /// `SameParameterEdge(edgeFirst, edgeLast, listFacesFirst, listFacesLast,
    /// secForward, whichSec, firstCall)` (`cxx:661-1186`).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn same_parameter_edge_pair(
        &mut self,
        edge_first: &TopoShape,
        edge_last: &TopoShape,
        list_faces_first: &[TopoShape],
        list_faces_last: &[TopoShape],
        sec_forward: bool,
        which_sec: &mut i32,
        first_call: bool,
    ) -> Option<Edge> {
        // Do not process floating edges (`cxx:608-610`).
        if list_faces_first.is_empty() || list_faces_last.is_empty() {
            return None;
        }
        let (e1, e2);
        if first_call {
            // Take the longest edge as first (`cxx:610-635`).
            let (f, l) = edge_range(edge_first);
            let len1 = match BRepTool::edge_curve_adaptor_world(&Edge(edge_first.clone())) {
                Some(c) => occt_geom::gcpnts::curve_length_range(c.as_ref(), f, l, CONFUSION),
                None => 0.0,
            };
            let (f, l) = edge_range(edge_last);
            let len2 = match BRepTool::edge_curve_adaptor_world(&Edge(edge_last.clone())) {
                Some(c) => occt_geom::gcpnts::curve_length_range(c.as_ref(), f, l, CONFUSION),
                None => 0.0,
            };
            if len1 < len2 {
                e1 = edge_last.clone();
                e2 = edge_first.clone();
                *which_sec = 2;
            } else {
                e1 = edge_first.clone();
                e2 = edge_last.clone();
                *which_sec = 1;
            }
        } else if *which_sec == 1 {
            e1 = edge_last.clone();
            e2 = edge_first.clone();
            *which_sec = 2;
        } else {
            e1 = edge_first.clone();
            e2 = edge_last.clone();
            *which_sec = 1;
        }

        let (first, last) = edge_range(&e1);
        let reg = GeometryRegistry::global();

        // To keep NM vertices on the edge (`cxx:657-661`).
        let (mut nm_verts, mut nm_pars) = find_nm_vertices(&e1);
        let (nm2_verts, nm2_pars) = find_nm_vertices(&e2);
        nm_verts.extend(nm2_verts);
        nm_pars.extend(nm2_pars);

        // Create the new edge (`cxx:664-680`).
        let mut edge = TopoShape::new(ShapeType::Edge);
        edge.set_orientation(e1.orientation());
        // `BRep_Tool::Curve(edge1, loc3d, first3d, last3d)` with the location
        // composed in, which is what the code right after does by hand
        // (`cxx:671-675`).
        let c3d_world = BRepTool::edge_curve_adaptor_world(&Edge(e1.clone()));
        let (first3d, last3d) = edge_range(&e1);
        let c3d = c3d_world?;
        let tol1 = BRepTool::edge_tolerance(&Edge(e1.clone()));
        reg.set_edge_curve3d(&edge, c3d.clone(), tol1);
        reg.set_edge_range(&edge, first, last);
        // `aBuilder.SameRange(edge, Standard_False)` / `SameParameter(...)`
        // (`cxx:679-680`).
        reg.set_same_range(&edge, false);
        reg.set_same_parameter(&edge, false);

        // Create and add the new vertices (`cxx:681-783`).
        let (v11, v12) = edge_vertices(&e1);
        let (v21, v22) = edge_vertices(&e2);
        let (v11, v12, v21, v22) = (v11?, v12?, v21?, v22?);
        let is_closed1 = crate::topo_tools::is_same(&v11, &v12);
        let is_closed2 = crate::topo_tools::is_same(&v21, &v22);
        if !is_closed1 && !is_closed2 {
            if sec_forward {
                if crate::topo_tools::is_same(&v11, &v22) || crate::topo_tools::is_same(&v12, &v21) {
                    return None;
                }
            } else if crate::topo_tools::is_same(&v11, &v21)
                || crate::topo_tools::is_same(&v12, &v22)
            {
                return None;
            }
        }
        let (v1_new, v2_new);
        if is_closed1 || is_closed2 {
            let v = if is_closed1 && is_closed2 {
                compute_tolerance_vertex2(&v11, &v21)
            } else if is_closed1 {
                compute_tolerance_vertex3(&v22, &v21, &v11)
            } else {
                compute_tolerance_vertex3(&v11, &v12, &v21)
            };
            v1_new = v.0.clone();
            v2_new = v.0;
        } else {
            let is_old_first = if sec_forward {
                crate::topo_tools::is_same(&v11, &v21)
            } else {
                crate::topo_tools::is_same(&v11, &v22)
            };
            let is_old_last = if sec_forward {
                crate::topo_tools::is_same(&v12, &v22)
            } else {
                crate::topo_tools::is_same(&v12, &v21)
            };
            let mut a = None;
            let mut b = None;
            if sec_forward {
                if !is_old_first {
                    a = Some(compute_tolerance_vertex2(&v11, &v21).0);
                }
                if !is_old_last {
                    b = Some(compute_tolerance_vertex2(&v12, &v22).0);
                }
            } else {
                if !is_old_first {
                    a = Some(compute_tolerance_vertex2(&v11, &v22).0);
                }
                if !is_old_last {
                    b = Some(compute_tolerance_vertex2(&v12, &v21).0);
                }
            }
            v1_new = if is_old_first { v11 } else { a? };
            v2_new = if is_old_last { v12 } else { b? };
        }
        {
            let builder = TopoBuilder::new();
            let mut an_edge = edge.oriented(Orientation::Forward);
            builder.add(&mut an_edge, &v1_new.oriented(Orientation::Forward));
            builder.add(&mut an_edge, &v2_new.oriented(Orientation::Reversed));
            for v in &nm_verts {
                builder.add(&mut an_edge, v);
            }
        }

        // Retrieve the second pcurves (`cxx:785-955`).
        let mut is_res_edge = false;
        let faces2: &[TopoShape] = if *which_sec == 1 {
            list_faces_last
        } else {
            list_faces_first
        };
        // `BRep_Tool::Surface(fac2, loc2)` survives the loop (`cxx:783-786`).
        let mut last_surf2: Option<Arc<dyn occt_geom::Surface>> = None;
        let mut last_loc2 = occt_core::toploc::TopLocLocation::identity();
        for fac2_shape in faces2 {
            let fac2 = Face(fac2_shape.clone());
            let Some(surf2) = BRepTool::face_surface(&fac2) else {
                continue;
            };
            last_surf2 = Some(surf2.clone());
            last_loc2 = fac2.0.location().clone();
            // `cxx:805-816`.
            let is_seam2 = (self.is_u_closed_surface(&surf2, &e2)
                || self.is_v_closed_surface(&surf2, &e2))
                && BRepTool::is_closed_edge_face(&Edge(e2.clone()), &fac2);
            let mut c2d21: Option<Arc<dyn Curve2d>> = None;
            let mut c2d21_range = (0.0f64, 0.0f64);
            if is_seam2 {
                if !self.nonmanifold {
                    return None;
                }
                let rev = e2.oriented(Orientation::Reversed);
                match pcurve_range_on_surface(&rev, &surf2) {
                    Some((c, f, l)) => {
                        c2d21 = Some(c);
                        c2d21_range = (f, l);
                    }
                    None => continue,
                }
            }
            let c2d2_rep = pcurve_range_on_surface(&e2, &surf2);
            if c2d2_rep.is_none() && c2d21.is_none() {
                continue;
            }
            let mut first_old;
            let mut last_old;
            if let Some(c21_in) = c2d21.take() {
                let mut c21 = copy_curve2d(&c21_in);
                first_old = c2d21_range.0;
                last_old = c2d21_range.1;
                if !sec_forward {
                    if c21.is_geom2d_line() {
                        c21 = Arc::new(Geom2dTrimmedCurve::new(c21, first_old, last_old));
                    }
                    let first2d = first_old;
                    let last2d = last_old;
                    first_old = c21.reversed_parameter(last2d);
                    last_old = c21.reversed_parameter(first2d);
                    c21 = Arc::from(c21.reversed());
                }
                c2d21 = self.same_range(&c21, first_old, last_old, first, last);
            }
            let Some((c2d2_raw, f2, l2)) = c2d2_rep else {
                continue;
            };
            first_old = f2;
            last_old = l2;
            let mut c2d2 = copy_curve2d(&c2d2_raw);
            if !sec_forward {
                if c2d2.is_geom2d_line() {
                    c2d2 = Arc::new(Geom2dTrimmedCurve::new(c2d2, first_old, last_old));
                }
                let first2d = first_old;
                let last2d = last_old;
                first_old = c2d2.reversed_parameter(last2d);
                last_old = c2d2.reversed_parameter(first2d);
                c2d2 = Arc::from(c2d2.reversed());
            }
            let Some(c2d2) = self.same_range(&c2d2, first_old, last_old, first, last) else {
                continue;
            };

            // Add the second pcurve (`cxx:852-955`).
            let mut is_seam = false;
            let mut ori: Orientation;
            let faces1: &[TopoShape] = if *which_sec == 1 {
                list_faces_first
            } else {
                list_faces_last
            };
            for fac1_shape in faces1 {
                let fac1 = Face(fac1_shape.clone());
                let Some(surf1) = BRepTool::face_surface(&fac1) else {
                    continue;
                };
                let is_seam1 = (self.is_u_closed_surface(&surf1, &e1)
                    || self.is_v_closed_surface(&surf1, &e1))
                    && BRepTool::is_closed_edge_face(&Edge(e1.clone()), &fac1);
                let c2d1_rep = pcurve_range_on_surface(&e1, &surf1);
                ori = e1.orientation();
                if fac1.0.orientation() == Orientation::Reversed {
                    ori = ori.reversed();
                }
                let mut c2d11: Option<Arc<dyn Curve2d>> = None;
                if is_seam1 {
                    if !self.nonmanifold {
                        return None;
                    }
                    let rev = e1.oriented(Orientation::Reversed);
                    c2d11 = pcurve_range_on_surface(&rev, &surf1).map(|(c, _, _)| c);
                    let fk = GeometryRegistry::shape_key(&fac1.0);
                    match (&c2d1_rep, &c2d11) {
                        (Some((c1, _, _)), Some(c11)) => {
                            if ori == Orientation::Forward {
                                reg.set_edge_pcurves(&edge, fk, vec![c1.clone(), c11.clone()]);
                            } else {
                                reg.set_edge_pcurves(&edge, fk, vec![c11.clone(), c1.clone()]);
                            }
                        }
                        _ => {}
                    }
                } else if let Some((c1, _, _)) = &c2d1_rep {
                    reg.set_edge_pcurve(&edge, GeometryRegistry::shape_key(&fac1.0), c1.clone());
                }
                let (c2d1, _c2d11) = match (&c2d1_rep, &c2d11) {
                    (None, None) => continue,
                    (a, b) => (
                        a.as_ref().map(|(c, _, _)| c.clone()),
                        b.clone(),
                    ),
                };

                if Arc::ptr_eq(&surf2, &surf1) && fac2.0.location().is_equal(&fac1.0.location()) {
                    let uclosed = self.is_u_closed_surface(&surf2, &e2);
                    let vclosed = self.is_v_closed_surface(&surf2, &e2);
                    if uclosed || vclosed {
                        // OCCT dereferences `c2d1` here; a null one keeps the
                        // non-seam path instead of panicking.
                        if let Some(c1) = c2d1.as_ref() {
                            let pf = c1.first_parameter();
                            let p1n = c1.d0(first.max(pf));
                            let p21n = c2d2.d0(first.max(c2d2.first_parameter()));
                            let p22n = c2d2.d0(last.min(c2d2.last_parameter()));
                            let a_dist = p1n.distance(&p21n).min(p1n.distance(&p22n));
                            let (u1, u2) = surf2.u_range();
                            let (v1r, v2r) = surf2.v_range();
                            is_seam = (uclosed && a_dist > 0.75 * (u2 - u1).abs())
                                || (vclosed && a_dist > 0.75 * (v2r - v1r).abs());
                            if !is_seam
                                && BRepTool::is_closed_edge_face(&Edge(edge.clone()), &fac1)
                            {
                                continue;
                            }
                        }
                    }
                }

                is_res_edge = true;
                let fk = GeometryRegistry::shape_key(&fac1.0);
                if is_seam {
                    let c1 = c2d1.as_ref().unwrap();
                    if ori == Orientation::Forward {
                        reg.set_edge_pcurves(&edge, fk, vec![c1.clone(), c2d2.clone()]);
                    } else {
                        reg.set_edge_pcurves(&edge, fk, vec![c2d2.clone(), c1.clone()]);
                    }
                } else if let (true, Some(c21)) = (is_seam2, c2d21.clone()) {
                    let mut init_ori = e2.orientation();
                    let sec_ori = edge.orientation();
                    let mut sec_ori = sec_ori;
                    if fac2.0.orientation() == Orientation::Reversed {
                        init_ori = init_ori.reversed();
                        sec_ori = sec_ori.reversed();
                    }
                    if !sec_forward {
                        init_ori = init_ori.reversed();
                    }
                    let _ = sec_ori;
                    if init_ori == Orientation::Forward {
                        reg.set_edge_pcurves(&edge, fk, vec![c2d2.clone(), c21]);
                    } else {
                        reg.set_edge_pcurves(&edge, fk, vec![c21, c2d2.clone()]);
                    }
                } else {
                    reg.set_edge_pcurve(&edge, fk, c2d2.clone());
                }
                if is_seam {
                    // `cxx:862`: the face loop runs while `!isSeam`.
                    break;
                }
            }
        }

        // `cxx:956-968`.
        let mut tol_reached = INFINITE;
        let mut is_same_par = false;
        if is_res_edge {
            self.same_parameter_edge_geom(&Edge(edge.clone()));
        }
        if BRepTool::same_parameter(&Edge(edge.clone())) {
            is_same_par = true;
            tol_reached = BRepTool::edge_tolerance(&Edge(edge.clone()));
        }

        if first_call && (!is_res_edge || !is_same_par || tol_reached > self.tolerance) {
            let mut which_secn = *which_sec;
            // Try to merge on the second section (`cxx:975-997`).
            let s_edge = self.same_parameter_edge_pair(
                edge_first,
                edge_last,
                list_faces_first,
                list_faces_last,
                sec_forward,
                &mut which_secn,
                false,
            );
            let mut second_ok = false;
            if let Some(s_edge) = s_edge {
                let tol_reached_2 = BRepTool::edge_tolerance(&s_edge);
                second_ok =
                    BRepTool::same_parameter(&s_edge) && tol_reached_2 < tol_reached;
                if second_ok {
                    edge = s_edge.0;
                    *which_sec = which_secn;
                    tol_reached = tol_reached_2;
                }
            }
            if !second_ok {
                // `cxx:1000-1056`: recompute the tolerance from the pcurves.
                let nbp = 23;
                let delta_t = (last3d - first3d) / (nbp - 1) as f64;
                let mut c3dpnt = Vec::with_capacity(nbp);
                for i in 1..=nbp {
                    c3dpnt.push(c3d.d0(first3d + (i - 1) as f64 * delta_t));
                }
                let mut max_tol = -1.0f64;
                // `cxx:1015-1042`: the pcurve index loop.
                let Some(surf2) = last_surf2.clone() else {
                    return Some(Edge(edge));
                };
                let reps = pcurves_on_surface(&edge, &surf2);
                let surf2_world: Arc<dyn occt_geom::Surface> = if !last_loc2.is_identity() {
                    Arc::from(surf2.transformed(&last_loc2.transformation()))
                } else {
                    surf2.clone()
                };
                for (c2d2, _pf, _pl) in reps {
                    let delta = (last - first) / (nbp - 1) as f64;
                    let mut dist2 = 0.0f64;
                    for i in 1..=nbp {
                        let p2d = c2d2.d0(first + (i - 1) as f64 * delta);
                        let p2 = surf2_world.d0(p2d.x(), p2d.y());
                        let dist = p2.square_distance(&c3dpnt[i - 1]);
                        if dist > dist2 {
                            dist2 = dist;
                        }
                    }
                    max_tol = (dist2.sqrt() * (1.0 + 1e-7)).max(CONFUSION);
                }
                if max_tol >= 0.0 && max_tol < tol_reached {
                    if tol_reached > self.max_tolerance() {
                        // Set the tolerance directly (`cxx:1044-1048`).
                        GeometryRegistry::global()
                            .set_edge_tolerance(&edge, max_tol);
                    } else {
                        // Just update the tolerance with the computed distance.
                        GeometryRegistry::global()
                            .raise_edge_tolerance(&edge, max_tol);
                    }
                    GeometryRegistry::global().set_same_parameter(&edge, true);
                }
            }
        }

        let tol_edge1 = BRepTool::edge_tolerance(&Edge(edge.clone()));
        if tol_edge1 > self.max_tolerance() {
            return None;
        }
        Some(Edge(edge))
    }

    /// The `(edge, surface)` list of pcurves used by the tolerance fallback:
    /// `BRep_Tool::CurveOnSurface(edge, c2d, surf, loc, f, l, j)` in index
    /// order (`cxx:1015-1018`).
    #[allow(dead_code)]
    pub(crate) fn pcurves_of(&self, edge: &TopoShape) -> Vec<Arc<dyn Curve2d>> {
        GeometryRegistry::global()
            .edge_pcurve_reps(edge)
            .into_iter()
            .map(|(_, c, _, _)| c)
            .collect()
    }
}

/// `BRep_Tool::CurveOnSurface(edge, c2d, surf, loc, first, last, index)` loop
/// of the tolerance fallback: every pcurve representation of `edge` on the
/// `surf` handle.
fn pcurves_on_surface(
    edge: &TopoShape,
    surf: &Arc<dyn occt_geom::Surface>,
) -> Vec<(Arc<dyn Curve2d>, f64, f64)> {
    let want = Arc::as_ptr(surf) as *const () as usize;
    GeometryRegistry::global()
        .edge_pcurve_reps(edge)
        .into_iter()
        .filter(|(s, _, _, _)| Arc::as_ptr(s) as *const () as usize == want)
        .map(|(_, c, f, l)| (c, f, l))
        .collect()
}

/// `TColStd_SequenceOfReal` / `TColStd_SequenceOfShape` pair used by the
/// merging stage; re-exported for the callers that build the sequences.


#[allow(dead_code)]
fn _unused(_: &ShapeKeyedMap<TopoShape>) {}




