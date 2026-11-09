//! `BRepFeat_RibSlot::ExtremeFaces`
//! (`ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx:747-1319`).
//!
//! Finds the two extreme faces / edges / vertices of the profile wire on the
//! basis: the support faces the tool starts from and ends on, plus the flags
//! telling whether the profile endpoint already lies on a vertex of the basis.
//!
//! A one-edge wire projects its two endpoints with `PtOnEdgeVertex` and falls
//! back to `LocOpe_CSIntersector` on the whole profile when a face is still
//! missing (`cxx:793-1074`). A multi-edge wire walks the basis faces and the
//! intersections of each profile edge with them (`cxx:1078-1319`).

use std::sync::Arc;

use occt_core::gp::{GpPln, GpPnt};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::trimmed::GeomTrimmedCurveBasis;

use crate::abs::ShapeType;
use crate::brep_feat::rib_slot::{
    choice_of_faces, int_par, pt_on_edge_vertex, RibSlot,
};
use crate::brep_int_curve_surface_inter::BRepIntCurveSurfaceInter;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::ShapeSet;
use crate::locope::CSIntersector;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::shhealing::{
    first_vertex as edge_first_vertex, last_vertex as edge_last_vertex,
};
use crate::topo_tools::{explorer_shapes, is_same, vertex_position};

/// The outputs of `ExtremeFaces` (`cxx:750-762`).
#[derive(Clone, Debug)]
pub struct ExtremeFaces {
    /// `FirstEdge`.
    pub first_edge: Option<Edge>,
    /// `LastEdge`.
    pub last_edge: Option<Edge>,
    /// `FirstFace`.
    pub first_face: Option<Face>,
    /// `LastFace`.
    pub last_face: Option<Face>,
    /// `FirstVertex` (in-out: the reference reads it back in the multi-edge
    /// branch at `cxx:1119`).
    pub first_vertex: Option<Vertex>,
    /// `LastVertex` (in-out, same as `FirstVertex`).
    pub last_vertex: Option<Vertex>,
    /// `OnFirstFace`.
    pub on_first_face: bool,
    /// `OnLastFace`.
    pub on_last_face: bool,
    /// `PtOnFirstEdge`.
    pub pt_on_first_edge: bool,
    /// `PtOnLastEdge`.
    pub pt_on_last_edge: bool,
    /// `OnFirstEdge`.
    pub on_first_edge: Option<Edge>,
    /// `OnLastEdge`.
    pub on_last_edge: Option<Edge>,
}

/// The basis faces that use `on` as a sub-shape, in `TopExp_Explorer` order.
///
/// `dedupe_faces` / `stop_on_first` mirror the two loops of each branch:
/// `cxx:845-867` runs a `ShapeMapHasher` map and scans every edge of a candidate;
/// `cxx:888-901` stops at the first matching vertex; `cxx:1164-1182` and
/// `cxx:1233-1249` run without a map and skip `exclude` by `IsSame`.
fn faces_containing_edge(
    sbase: &TopoShape,
    on: &Edge,
    exclude: Option<&Face>,
    dedupe_faces: bool,
    stop_on_first: bool,
) -> Vec<TopoShape> {
    let mut faces = Vec::new();
    let mut map = ShapeSet::new();
    for f in explorer_shapes(sbase, ShapeType::Face) {
        if dedupe_faces && !map.add(&f) {
            continue;
        }
        if let Some(ex) = exclude {
            if is_same(&f, &ex.0) {
                continue;
            }
        }
        for e in explorer_shapes(&f, ShapeType::Edge) {
            if is_same(&e, &on.0) {
                faces.push(f.clone());
                if stop_on_first {
                    break;
                }
            }
        }
    }
    faces
}

/// The vertex twin of [`faces_containing_edge`].
fn faces_containing_vertex(
    sbase: &TopoShape,
    on: &Vertex,
    exclude: Option<&Face>,
    dedupe_faces: bool,
    stop_on_first: bool,
) -> Vec<TopoShape> {
    let mut faces = Vec::new();
    let mut map = ShapeSet::new();
    for f in explorer_shapes(sbase, ShapeType::Face) {
        if dedupe_faces && !map.add(&f) {
            continue;
        }
        if let Some(ex) = exclude {
            if is_same(&f, &ex.0) {
                continue;
            }
        }
        for v in explorer_shapes(&f, ShapeType::Vertex) {
            if is_same(&v, &on.0) {
                faces.push(f.clone());
                if stop_on_first {
                    break;
                }
            }
        }
    }
    faces
}

impl RibSlot {
    /// `ExtremeFaces(RevolRib, bnd, Pln, FirstEdge, LastEdge, FirstFace,
    /// LastFace, FirstVertex, LastVertex, OnFirstFace, OnLastFace,
    /// PtOnFirstEdge, PtOnLastEdge, OnFirstEdge, OnLastEdge)` (`cxx:747-1319`).
    ///
    /// `Ok(false, _)` is the reference's `Data = false`: no pair of extreme
    /// faces could be found. `Err` is a raised `Standard_*` failure
    /// (`Standard_OutOfRange` of `LocOpe_CSIntersector::Point` at `cxx:1032`,
    /// `Standard_ConstructionError` of its `Perform`).
    pub fn extreme_faces(
        &self,
        revol_rib: bool,
        bnd: f64,
        pln: &GpPln,
        first_vertex: Option<Vertex>,
        last_vertex: Option<Vertex>,
    ) -> Result<(bool, ExtremeFaces), String> {
        // `cxx:770-777`.
        let mut out = ExtremeFaces {
            first_edge: None,
            last_edge: None,
            first_face: None,
            last_face: None,
            first_vertex,
            last_vertex,
            on_first_face: false,
            on_last_face: false,
            pt_on_first_edge: false,
            pt_on_last_edge: false,
            on_first_edge: None,
            on_last_edge: None,
        };

        // A null `mySbase` cannot host a support face; the reference would
        // simply explore it and find nothing.
        let Some(sbase) = self.my_sbase.clone() else {
            return Ok((false, out));
        };
        let b = TopoBuilder::new();
        let edges: Vec<TopoShape> = match &self.my_wire {
            Some(w) => explorer_shapes(&w.0, ShapeType::Edge),
            None => Vec::new(),
        };

        if edges.len() == 1 {
            // ---the wire includes only one edge (`cxx:793-1074`)
            let e = Edge::wrap(edges[0].clone()).ok_or_else(|| "Standard_DomainError".to_string())?;
            let (first_par, last_par) = BRepTool::edge_parameters(&e);
            let cc =
                BRepTool::edge_curve_world(&e).ok_or_else(|| "Standard_NullObject".to_string())?;
            let p1 = edge_first_vertex(&e)
                .map(|v| vertex_position(&v))
                .ok_or_else(|| "Standard_NullObject: TopExp::FirstVertex".to_string())?;
            let p2 = edge_last_vertex(&e)
                .map(|v| vertex_position(&v))
                .ok_or_else(|| "Standard_NullObject: TopExp::LastVertex".to_string())?;

            // `cxx:817-836`.
            let r1 = pt_on_edge_vertex(revol_rib, &sbase, &p1);
            let r2 = pt_on_edge_vertex(revol_rib, &sbase, &p2);

            // `cxx:839-917`.
            if r1.pt_on_edge {
                // `cxx:842-851` / `cxx:853-917` both replace the whole list.
                let faces = if !r1.pt_on_vertex {
                    let on_edge = r1
                        .on_edge
                        .clone()
                        .ok_or_else(|| "Standard_NullObject".to_string())?;
                    faces_containing_edge(&sbase, &on_edge, None, true, false)
                } else {
                    let on_vertex = r1
                        .on_vertex
                        .clone()
                        .ok_or_else(|| "Standard_NullObject".to_string())?;
                    faces_containing_vertex(&sbase, &on_vertex, None, true, true)
                };
                if !faces.is_empty() {
                    if let Some(fff) =
                        choice_of_faces(&faces, &cc, first_par + bnd / 50.0, bnd / 50.0, pln)?
                    {
                        out.first_face = Some(fff);
                    }
                }
                out.first_edge = Some(e.clone());
                // `BRepLib_MakeVertex v(p1); FirstVertex = v;`
                out.first_vertex = Some(b.make_vertex(p1, CONFUSION));
                out.on_first_face = true;
            }

            // `cxx:919-997`.
            if r2.pt_on_edge {
                // `cxx:922-931` / `cxx:933-997` both replace the whole list.
                let faces = if !r2.pt_on_vertex {
                    let on_edge = r2
                        .on_edge
                        .clone()
                        .ok_or_else(|| "Standard_NullObject".to_string())?;
                    faces_containing_edge(&sbase, &on_edge, None, true, false)
                } else {
                    let on_vertex = r2
                        .on_vertex
                        .clone()
                        .ok_or_else(|| "Standard_NullObject".to_string())?;
                    faces_containing_vertex(&sbase, &on_vertex, None, true, true)
                };
                if !faces.is_empty() {
                    if let Some(fff) =
                        choice_of_faces(&faces, &cc, last_par - bnd / 50.0, bnd / 50.0, pln)?
                    {
                        out.last_face = Some(fff);
                    }
                }
                out.last_edge = Some(e.clone());
                out.last_vertex = Some(b.make_vertex(p2, CONFUSION));
                out.on_last_face = true;
            }

            // `cxx:999-1003`.
            if out.first_face.is_some() && out.last_face.is_some() {
                return Ok((true, out));
            }

            // `cxx:1008-1059`: the profile itself against the basis.
            let mut asi = CSIntersector::with_shape(&sbase);
            asi.perform_curves(&[Some(cc.clone())])?;
            let npts = asi.nb_points(1).unwrap_or(0);
            if !asi.is_done() || npts < 2 {
                // `Data = false; return Data;`
                return Ok((false, out));
            }
            let mut lastpar = asi.point(1, npts)?.parameter();
            let mut lastindex = npts;
            if lastpar > last_par {
                // `for (int jj = NbPoints - 1; jj >= 1; jj--)` (`cxx:1020`).
                for jj in (1..npts).rev() {
                    let par = asi.point(1, jj)?.parameter();
                    if par <= last_par {
                        lastpar = par;
                        lastindex = jj;
                        break;
                    }
                }
            }
            let firstindex = lastindex - 1;
            let firstpar = asi.point(1, firstindex)?.parameter();

            let mut firstpoint: Option<GpPnt> = None;
            let mut lastpoint: Option<GpPnt> = None;
            if out.first_face.is_none() {
                out.first_face = Some(asi.point(1, firstindex)?.face().clone());
                let pp = cc.d0(firstpar);
                firstpoint = Some(pp);
                out.first_vertex = Some(b.make_vertex(pp, CONFUSION));
                out.first_edge = Some(e.clone());
            }
            if out.last_face.is_none() {
                out.last_face = Some(asi.point(1, lastindex)?.face().clone());
                let pp = cc.d0(lastpar);
                lastpoint = Some(pp);
                out.last_vertex = Some(b.make_vertex(pp, CONFUSION));
                out.last_edge = Some(e.clone());
            }

            // `cxx:1063-1071`: the point read only exists when the face was
            // missing, which is exactly when the flag is still false.
            if !out.on_first_face {
                out.on_first_face = firstpoint
                    .map(|fp| p1.distance(&fp) <= CONFUSION)
                    .unwrap_or(false);
            }
            if !out.on_last_face {
                out.on_last_face = lastpoint
                    .map(|lp| p2.distance(&lp) <= CONFUSION)
                    .unwrap_or(false);
            }

            let data = out.first_face.is_some() && out.last_face.is_some();
            return Ok((data, out));
        }

        // ---The wire consists of several edges (`cxx:1078-1319`)
        let Some(wire) = self.my_wire.clone() else {
            return Ok((false, out));
        };
        let mut first_ok = false;
        let mut last_ok = false;
        let mut ex = WireExplorer::of_wire(&wire);
        while ex.more() {
            let Some(e) = ex.current().cloned() else {
                ex.next();
                continue;
            };
            let (mut f, mut l) = BRepTool::edge_parameters(&e);
            let cur =
                BRepTool::edge_curve_world(&e).ok_or_else(|| "Standard_NullObject".to_string())?;
            // `f = f - bnd/10000; l = l + bnd/10000;`
            f -= bnd / 10000.0;
            l += bnd / 10000.0;
            let curve: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(cur, f, l));
            let p2 = edge_last_vertex(&e)
                .map(|v| vertex_position(&v))
                .ok_or_else(|| "Standard_NullObject: TopExp::LastVertex".to_string())?;

            for face_sh in explorer_shapes(&sbase, ShapeType::Face) {
                let Some(a_cur_face) = Face::wrap(face_sh.clone()) else {
                    continue;
                };
                let face_tol = BRepTool::face_tolerance(&a_cur_face);
                let mut inter = BRepIntCurveSurfaceInter::new();
                inter.init_shape_curve(&a_cur_face.0, curve.clone(), face_tol);
                if !inter.more() {
                    continue;
                }
                while inter.more() {
                    let the_point = inter.pnt()?;
                    // `cxx:1119-1127`.
                    if let Some(fv) = &out.first_vertex {
                        let point = vertex_position(fv);
                        if point.distance(&the_point) <= face_tol {
                            inter.next();
                            continue;
                        }
                    }
                    let intpar = int_par(Some(curve.as_ref()), &the_point);
                    let the_edge = Some(e.clone());
                    let the_face = Some(a_cur_face.clone());
                    let the_vertex = Some(b.make_vertex(the_point, CONFUSION));
                    // `cxx:1132-1137`.
                    if !first_ok && the_point.distance(&p2) <= CONFUSION {
                        inter.next();
                        continue;
                    }
                    // `cxx:1140-1149`.
                    let r = pt_on_edge_vertex(revol_rib, &a_cur_face.0, &the_point);

                    // `cxx:1152-1220`.
                    if out.first_edge.is_none()
                        && the_edge.is_some()
                        && the_face.is_some()
                        && the_vertex.is_some()
                    {
                        out.first_edge = the_edge.clone();
                        out.first_face = the_face.clone();
                        out.first_vertex = the_vertex.clone();
                        out.pt_on_first_edge = r.pt_on_edge;
                        out.on_first_edge = r.on_edge.clone();
                        if r.pt_on_edge && !r.pt_on_vertex {
                            let on_edge = r
                                .on_edge
                                .clone()
                                .ok_or_else(|| "Standard_NullObject".to_string())?;
                            let mut faces = vec![out
                                .first_face
                                .clone()
                                .ok_or_else(|| "Standard_NullObject".to_string())?
                                .0
                                .clone()];
                            faces.extend(faces_containing_edge(
                                &sbase,
                                &on_edge,
                                out.first_face.as_ref(),
                                false,
                                false,
                            ));
                            if let Some(fff) = choice_of_faces(
                                &faces,
                                &curve,
                                intpar + bnd / 10.0,
                                bnd / 10.0,
                                pln,
                            )? {
                                out.first_face = Some(fff);
                            }
                        } else if r.pt_on_edge && r.pt_on_vertex {
                            let on_vertex = r
                                .on_vertex
                                .clone()
                                .ok_or_else(|| "Standard_NullObject".to_string())?;
                            let mut faces = vec![out
                                .first_face
                                .clone()
                                .ok_or_else(|| "Standard_NullObject".to_string())?
                                .0
                                .clone()];
                            faces.extend(faces_containing_vertex(
                                &sbase,
                                &on_vertex,
                                out.first_face.as_ref(),
                                false,
                                false,
                            ));
                            if let Some(fff) = choice_of_faces(
                                &faces,
                                &curve,
                                intpar + bnd / 10.0,
                                bnd / 10.0,
                                pln,
                            )? {
                                out.first_face = Some(fff);
                            }
                        }
                        if out.first_edge.is_some()
                            && out.first_face.is_some()
                            && out.first_vertex.is_some()
                        {
                            first_ok = true;
                        }
                    }

                    // `cxx:1221-1289`.
                    if out.last_edge.is_none()
                        && the_edge.is_some()
                        && the_face.is_some()
                        && the_vertex.is_some()
                        && out.first_edge.is_some()
                    {
                        out.last_edge = the_edge.clone();
                        out.last_face = the_face.clone();
                        out.last_vertex = the_vertex.clone();
                        out.pt_on_last_edge = r.pt_on_edge;
                        out.on_last_edge = r.on_edge.clone();
                        if r.pt_on_edge && !r.pt_on_vertex {
                            let on_edge = r
                                .on_edge
                                .clone()
                                .ok_or_else(|| "Standard_NullObject".to_string())?;
                            let mut faces = vec![out
                                .last_face
                                .clone()
                                .ok_or_else(|| "Standard_NullObject".to_string())?
                                .0
                                .clone()];
                            faces.extend(faces_containing_edge(
                                &sbase,
                                &on_edge,
                                out.last_face.as_ref(),
                                false,
                                false,
                            ));
                            if let Some(fff) = choice_of_faces(
                                &faces,
                                &curve,
                                intpar - bnd / 10.0,
                                bnd / 10.0,
                                pln,
                            )? {
                                out.last_face = Some(fff);
                            }
                        } else if r.pt_on_edge && r.pt_on_vertex {
                            let on_vertex = r
                                .on_vertex
                                .clone()
                                .ok_or_else(|| "Standard_NullObject".to_string())?;
                            let mut faces = vec![out
                                .last_face
                                .clone()
                                .ok_or_else(|| "Standard_NullObject".to_string())?
                                .0
                                .clone()];
                            faces.extend(faces_containing_vertex(
                                &sbase,
                                &on_vertex,
                                out.last_face.as_ref(),
                                false,
                                false,
                            ));
                            if let Some(fff) = choice_of_faces(
                                &faces,
                                &curve,
                                intpar - bnd / 10.0,
                                bnd / 10.0,
                                pln,
                            )? {
                                out.last_face = Some(fff);
                            }
                        }
                        if out.last_edge.is_some()
                            && out.last_face.is_some()
                            && out.last_vertex.is_some()
                        {
                            last_ok = true;
                        }
                        break;
                    }
                    inter.next();
                }
            }
            ex.next();
        }

        // `cxx:1292-1318`.
        if first_ok && last_ok {
            let first_edge = out
                .first_edge
                .clone()
                .ok_or_else(|| "Standard_NullObject".to_string())?;
            let last_edge = out
                .last_edge
                .clone()
                .ok_or_else(|| "Standard_NullObject".to_string())?;
            let pp1 = edge_first_vertex(&first_edge)
                .map(|v| vertex_position(&v))
                .ok_or_else(|| "Standard_NullObject: TopExp::FirstVertex".to_string())?;
            let pp2 = edge_last_vertex(&last_edge)
                .map(|v| vertex_position(&v))
                .ok_or_else(|| "Standard_NullObject: TopExp::LastVertex".to_string())?;
            let p1 = out
                .first_vertex
                .as_ref()
                .map(vertex_position)
                .ok_or_else(|| "Standard_NullObject".to_string())?;
            let p2 = out
                .last_vertex
                .as_ref()
                .map(vertex_position)
                .ok_or_else(|| "Standard_NullObject".to_string())?;
            if p1.distance(&pp1)
                <= BRepTool::face_tolerance(
                    out.first_face
                        .as_ref()
                        .ok_or_else(|| "Standard_NullObject".to_string())?,
                )
            {
                out.on_first_face = true;
            }
            if p2.distance(&pp2)
                <= BRepTool::face_tolerance(
                    out.last_face
                        .as_ref()
                        .ok_or_else(|| "Standard_NullObject".to_string())?,
                )
            {
                out.on_last_face = true;
            }
            return Ok((true, out));
        }
        Ok((false, out))
    }
}
