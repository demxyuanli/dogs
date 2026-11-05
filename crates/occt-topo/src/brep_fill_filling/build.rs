//! `BRepFill_Filling::Build` — the driver of the N-side filling.
//!
//! Source: `TKBool/BRepFill/BRepFill_Filling.cxx:571-777` (OCCT 8.0.0).
//!
//! # Translation notes
//!
//! * The plate solve runs through
//!   [`BuildPlateSurface`](occt_geom::geom_plate::BuildPlateSurface); the
//!   approximation through [`MakeApprox`](occt_geom::geom_plate::MakeApprox)
//!   with the reference's `Continuity = GeomAbs_C1` and
//!   `EnlargeCoeff = 1.1` defaults (`GeomPlate_MakeApprox.hxx:45-57`).
//! * `GeomAPI_ProjectPointOnSurf` is
//!   [`project_point_on_surface`](occt_geom::geom_api::project_point_on_surface):
//!   the reference keeps one `Projector` and re-`Init`s it per vertex; the port
//!   builds the engine per call, which is the same `Extrema_ExtPS` window
//!   (`GeomAPI_ProjectPointOnSurf.cxx:88-96`). A failed projection is
//!   `LowerDistance() > Precision::Confusion()`: OCCT raises `StdFail_NotDone`
//!   when the extrema engine found nothing, the port takes the same `continue`
//!   as the "too far" arm.
//! * `BRepLib_MakeEdge(Line2d, CurSurface, FirstVtx, LastVtx)` (`cxx:657`) is
//!   [`make_edge_on_surface`], the `Init(C2d, S, V1, V2)` path of
//!   `BRepLib_MakeEdge.cxx:1045-1215`.
//! * `BRep_Builder::UpdateEdge(E, C2d, S, L, Tol)` is
//!   [`GeometryRegistry::update_edge_pcurve_on_surface`]; `BRep_Builder::Add(E,
//!   V)` is [`TopoBuilder::add`] on the edge shape with the vertex orientation
//!   the reference writes (`FORWARD` then `REVERSED`).
//! * `BRepLib::SameParameter(NewEdge, dmax, true)` (`cxx:770`) is
//!   [`same_parameter_edge_inplace`] — the reference reaches
//!   `InternalSameParameter` through the single-shape overload (`BRepLib.cxx:906-910`);
//!   for one edge the two entry points do the same work.
//! * `BRepLib_MakeFace(Surface, FinalWire)` (`cxx:777`) is
//!   [`make_face_from_surface_and_wire`](crate::brep_lib_make_face::make_face_from_surface_and_wire).
//! * `GeomPlate_MakeApprox` never throws: a failed approximation leaves
//!   `mySurface` null (`GeomPlate_MakeApprox.cxx:325`), and
//!   `BRepLib_MakeFace(null, W)` then leaves `myShape` null while `myIsDone`
//!   stays `true` (`cxx:680-688` runs before the face is built). The port
//!   reproduces that observable state: `Face()` answers `None`, `IsDone()` stays
//!   true.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpXY, GpXyz};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{Precision, CONFUSION, PCONFUSION};
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::geom_plate::{BuildPlateSurface, MakeApprox, PlateG0Criterion};
use occt_geom::surface::Surface;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::brep_fill_filling::{
    build_wires, find_extremities_of_holes, wire_from_list, EdgeFaceAndOrder, Filling,
};
use crate::brep_lib_make_face::make_face_from_surface_and_wire;
use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::shape_fix_compose_shell::empty_copied_edge;
use crate::tgeometry::{EdgeGeom, GeometryRegistry};
use crate::topo_tools::{edge_vertices, shape_map_key};

impl Filling {
    /// `void BRepFill_Filling::Build()` (`cxx:571-777`).
    ///
    /// `Err` carries the two `Standard_Failure` messages the reference throws
    /// from `AddConstraints` (`cxx:338`) and the `Standard_OutOfRange` of
    /// `VerSeq(i + 1)` on an odd-length extremity sequence.
    pub fn build(&mut self) -> Result<(), String> {
        // `cxx:574-582`.
        self.builder = Some(BuildPlateSurface::new(
            self.degree,
            self.nb_pts_on_cur,
            self.nb_iter,
            self.tol2d,
            self.tol3d,
            self.tol_ang,
            self.anisotropie,
        )?);

        // `cxx:584-588`.
        if self.boundary.is_empty() {
            self.is_done = false;
            return Ok(());
        }

        // `cxx:594`.
        let mut ver_seq: Vec<TopoShape> = Vec::new();

        // `cxx:596-601`.
        let mut edge_list: Vec<TopoShape> = self.boundary.iter().map(|b| b.edge.0.clone()).collect();

        // `cxx:603-604`.
        let wire_list = build_wires(&mut edge_list, &mut self.old_new_map);
        find_extremities_of_holes(&wire_list, &mut ver_seq);

        // `cxx:606-660`: build the missing bounds out of the free constraints.
        for j in 0..self.free_constraints.len() {
            let cur_face = self.free_constraints[j].face.clone();
            let order = self.free_constraints[j].order;
            // `cxx:612-617`: `BRepAdaptor_Surface(F)` -> `BRep_Tool::Surface`.
            let Some(cur_surface) = BRepTool::face_surface_world(&cur_face) else {
                continue;
            };

            let mut i = 0usize;
            while i < ver_seq.len() {
                let first_vtx = Vertex::wrap(ver_seq[i].clone()).ok_or_else(|| {
                    "Standard_OutOfRange: BRepFill_Filling::Build: VerSeq entry is not a vertex"
                        .to_string()
                })?;
                let last_vtx = ver_seq.get(i + 1).cloned().ok_or_else(|| {
                    "Standard_OutOfRange: BRepFill_Filling::Build: VerSeq(i + 1)".to_string()
                })?;
                let last_vtx = Vertex::wrap(last_vtx).ok_or_else(|| {
                    "Standard_OutOfRange: BRepFill_Filling::Build: VerSeq entry is not a vertex"
                        .to_string()
                })?;

                // `cxx:621-629`.
                let first_pnt = BRepTool::vertex_point_world(&first_vtx);
                let Some(p1) = project_point_on_surface(cur_surface.as_ref(), &first_pnt, CONFUSION)
                else {
                    i += 2;
                    continue;
                };
                if p1.distance > CONFUSION {
                    i += 2;
                    continue;
                }
                let (u1, v1) = (p1.u, p1.v);

                // `cxx:636-645`.
                let last_pnt = BRepTool::vertex_point_world(&last_vtx);
                let Some(p2) = project_point_on_surface(cur_surface.as_ref(), &last_pnt, CONFUSION)
                else {
                    i += 2;
                    continue;
                };
                if p2.distance > CONFUSION {
                    i += 2;
                    continue;
                }
                let (u2, v2) = (p2.u, p2.v);

                // `cxx:650-658`.
                let points = vec![GpPnt2d::new(u1, v1), GpPnt2d::new(u2, v2)];
                let line2d = Geom2dBezierCurve::new(points)
                    .map_err(|e| format!("Standard_ConstructionError: Geom2d_BezierCurve: {e}"))?;
                let line2d: Arc<dyn Curve2d> = Arc::new(line2d);
                let e = make_edge_on_surface(&line2d, &cur_surface, &first_vtx, &last_vtx)?;
                self.add_with_support(&e, &cur_face, order, true);
                ver_seq.remove(i);
                ver_seq.remove(i);
                break;
            }
        }

        // `cxx:663-668`: load the initial surface when one was given.
        if self.is_init_face_given {
            if let Some(init_face) = self.init_face.clone() {
                if let Some(surf) = BRepTool::face_surface_world(&init_face) {
                    if let Some(builder) = self.builder.as_mut() {
                        builder.load_init_surface(surf);
                    }
                }
            }
        }

        // `cxx:670-678`.
        let boundary: Vec<EdgeFaceAndOrder> = self.boundary.clone();
        let constraints: Vec<EdgeFaceAndOrder> = self.constraints.clone();
        let init_face = if self.is_init_face_given {
            self.init_face.clone()
        } else {
            None
        };
        let builder = self.builder.as_mut().expect("builder set above");
        Filling::add_constraints(
            builder,
            &boundary,
            self.nb_pts_on_cur,
            self.tol3d,
            self.tol_ang,
            self.tol_curv,
            init_face.as_ref(),
        )?;
        builder.set_nb_bounds(boundary.len() as i32);
        Filling::add_constraints(
            builder,
            &constraints,
            self.nb_pts_on_cur,
            self.tol3d,
            self.tol_ang,
            self.tol_curv,
            init_face.as_ref(),
        )?;
        // `cxx:671-677`: `myBuilder->Add(myPoints(i))` on the reference's shared
        // handles; the port clones the constraints into the builder and keeps
        // `myPoints` populated for a later `Build`.
        for i in 0..self.points.len() {
            builder.add_point_constraint(self.points[i].clone());
        }

        // `cxx:679-688`.
        builder.perform()?;
        if !builder.is_done() {
            self.is_done = false;
            return Ok(());
        }
        self.is_done = true;

        // `cxx:690-714`: approximation of the plate surface.
        let Some(g_plate) = builder.surface().cloned() else {
            // `myBuilder->Surface()` on a done builder always answers a
            // handle; a null `GeomPlate_Surface` here would have thrown in
            // OCCT's `MakeApprox` constructor.
            return Ok(());
        };
        let dmax = 1.1 * builder.g0_error();
        let approx: Option<GeomBSplineSurface> = if !self.is_init_face_given {
            let mut s2d: Vec<GpXY> = Vec::new();
            builder.disc2d_contour(&mut s2d)?;
            let mut s3d: Vec<GpXyz> = Vec::new();
            builder.disc3d_contour(0, &mut s3d)?;
            let seuil = self.tol3d.max(10.0 * builder.g0_error());
            let criterion = PlateG0Criterion::with_defaults(
                s2d.iter().map(|p| (p.x, p.y)).collect(),
                s3d,
                seuil,
            );
            MakeApprox::from_criterion(
                g_plate,
                &criterion,
                self.tol3d,
                self.max_segments,
                self.max_deg,
                Shape::C1,
                1.1,
            )
            .ok()
            .and_then(|a| a.into_surface())
        } else {
            // `GeomPlate_MakeApprox(GPlate, myTol3d, myMaxSegments, myMaxDeg,
            // dmax, 0)` (`cxx:711`).
            MakeApprox::from_crit_order(
                g_plate,
                self.tol3d,
                self.max_segments,
                self.max_deg,
                dmax,
                0,
                Shape::C1,
                1.1,
            )
            .ok()
            .and_then(|a| a.into_surface())
        };

        let Some(bspline) = approx else {
            // `Surface` stays null; see the module note.
            self.face = None;
            return Ok(());
        };
        let surface: Arc<dyn Surface> = Arc::new(bspline);

        // `cxx:717-719`.
        let curves_on_plate = builder.curves2d()?;

        // `cxx:721-775`.
        let b = TopoBuilder::new();
        let reg = GeometryRegistry::global();
        let mut final_edges: Vec<TopoShape> = Vec::new();
        let boundary = self.boundary.clone();
        for (i, record) in boundary.iter().enumerate() {
            let init_edge = record.edge.clone();
            let mut an_edge = init_edge.clone();
            an_edge.0.set_orientation(Orientation::Forward);
            if let Some(mapped) = self.old_new_map.get(&shape_map_key(&an_edge.0)).cloned() {
                an_edge = Edge::wrap(mapped)
                    .ok_or_else(|| "BRepFill_Filling::Build: mapped shape is not an edge".to_string())?;
            }

            let a_curve_on_plate = curves_on_plate.get(i).cloned().flatten();

            let mut new_edge = empty_copied_edge(&an_edge.0);
            let (v1, v2) = edge_vertices(&an_edge);

            let new_v1 = Self::new_or_mapped_vertex(&mut self.old_new_map, v1.as_ref(), dmax)?;
            let new_v2 = Self::new_or_mapped_vertex(&mut self.old_new_map, v2.as_ref(), dmax)?;

            if let Some(v) = new_v1.as_ref() {
                b.add(&mut new_edge, &v.0.oriented(Orientation::Forward));
            }
            if let Some(v) = new_v2.as_ref() {
                b.add(&mut new_edge, &v.0.oriented(Orientation::Reversed));
            }
            // `BB.UpdateEdge(NewEdge, aCurveOnPlate, Surface, Loc, dmax)`.
            reg.update_edge_pcurve_on_surface(&new_edge, &surface, a_curve_on_plate, dmax);
            // `BRepLib::SameParameter(NewEdge, dmax, true)`.
            let edge = Edge(new_edge.clone());
            same_parameter_edge_inplace(&edge, dmax);

            final_edges.push(new_edge);
            self.old_new_map
                .insert(shape_map_key(&init_edge.0), edge.0.oriented(Orientation::Forward));
        }

        // `cxx:775-777`.
        let final_wire: Wire = wire_from_list(&final_edges);
        self.face = Some(make_face_from_surface_and_wire(surface, &final_wire));
        Ok(())
    }

    /// The vertex half of the `cxx:737-762` block: reuse the mapped vertex when
    /// `myOldNewMap` has one, otherwise `BRepLib_MakeVertex(BRep_Tool::Pnt(V))`
    /// `+ BB.UpdateVertex(NewV, dmax)` and remember the mapping.
    fn new_or_mapped_vertex(
        old_new_map: &mut crate::brep_fill_filling::OldNewMap,
        v: Option<&Vertex>,
        dmax: f64,
    ) -> Result<Option<Vertex>, String> {
        let Some(v) = v else {
            return Ok(None);
        };
        if let Some(mapped) = old_new_map.get(&shape_map_key(&v.0)).cloned() {
            return Vertex::wrap(mapped)
                .map(Some)
                .ok_or_else(|| "BRepFill_Filling::Build: mapped shape is not a vertex".to_string());
        }
        let a_pnt = BRepTool::vertex_point_world(v);
        let b = TopoBuilder::new();
        let new_v = b.make_vertex(a_pnt, CONFUSION);
        // `BB.UpdateVertex(NewV1, dmax)` (`cxx:747` / `:759`): raise-only.
        GeometryRegistry::global().raise_vertex_tolerance(&new_v.0, dmax);
        old_new_map.insert(shape_map_key(&v.0), new_v.0.clone());
        Ok(Some(new_v))
    }
}

/// `BRepLib_MakeEdge(C2d, S, V1, V2)` (`BRepLib_MakeEdge.cxx:1105-1112`): the
/// two vertices are projected onto the curve-on-surface and the edge is built
/// with the given handles.
///
/// The `while (CT = down_cast<Geom2d_TrimmedCurve>` strip loop (`cxx:1109-1115`)
/// is vacuous on this call site (the curve is the fresh two-pole
/// `Geom2d_BezierCurve` of `cxx:655`) and the port's `dyn Curve2d` has no
/// downcast, so the basis-curve strip is not written.
fn make_edge_on_surface(
    pc: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    vv1: &Vertex,
    vv2: &Vertex,
) -> Result<Edge, String> {
    let p1 = project_vertex_on_curve_on_surface(pc, surf, vv1)?;
    let p2 = project_vertex_on_curve_on_surface(pc, surf, vv2)?;
    init_2d_with_vertices(pc, surf, vv1, vv2, p1, p2)
}

/// `static bool Project(C2d, S, V, p)` (`BRepLib_MakeEdge.cxx:115-176`): the
/// parameter of `V` on the curve-on-surface.
///
/// Only the two end-distance tests are written. The `Extrema_ExtPC` fallback
/// (`cxx:151-172`) has no port, and it is unreachable here: the pcurve is the
/// two-pole Bezier whose poles are exactly the projections of `vv1` / `vv2`
/// (`BRepFill_Filling.cxx:650-657`), so one of the two ends is always within the
/// vertex tolerance.
fn project_vertex_on_curve_on_surface(
    pc: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    v: &Vertex,
) -> Result<f64, String> {
    let acos = CurveOnSurface::new(pc.clone(), surf.clone());
    let p = BRepTool::vertex_point_world(v);
    let tol = BRepTool::vertex_tolerance(v);
    let eps2 = tol * tol;

    let first = acos.first_parameter();
    let last = acos.last_parameter();
    let p_first = acos.d0(first);
    let p_last = acos.d0(last);
    let d1 = p_first.square_distance(&p);
    let d2 = p_last.square_distance(&p);
    if d1 < d2 && d1 <= eps2 {
        return Ok(first);
    }
    if d2 < d1 && d2 <= eps2 {
        return Ok(last);
    }
    Err(
        "BRepLib_PointProjectionFailed: BRepLib_MakeEdge::Project (UNPORTED Extrema_ExtPC arm)"
            .to_string(),
    )
}

/// `BRepLib_MakeEdge::Init(C2d, S, V1, V2, p1, p2)` (`BRepLib_MakeEdge.cxx:1045-1215`).
///
/// `Err` carries the `BRepLib_EdgeError` the reference stores when it bails out
/// before `Done()`; OCCT then hands back a null edge, the port reports the same
/// state as an error.
fn init_2d_with_vertices(
    pc: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    vv1: &Vertex,
    vv2: &Vertex,
    pp1: f64,
    pp2: f64,
) -> Result<Edge, String> {
    let epsilon = PCONFUSION;
    let cf = pc.first_parameter();
    let cl = pc.last_parameter();
    let mut p1 = pp1;
    let mut p2 = pp2;
    let v1;
    let v2;
    let mut reverse = false;
    if pc.is_periodic() {
        // `cxx:1094-1100`.
        occt_core::elib::clib2d::adjust_periodic(cf, cl, epsilon, &mut p1, &mut p2);
        v1 = vv1.clone();
        v2 = vv2.clone();
    } else {
        // `cxx:1102-1127`.
        if p1 < p2 {
            v1 = vv1.clone();
            v2 = vv2.clone();
        } else {
            v2 = vv1.clone();
            v1 = vv2.clone();
            std::mem::swap(&mut p1, &mut p2);
            reverse = true;
        }
        if cf - p1 > epsilon || p2 - cl > epsilon {
            return Err("BRepLib_ParameterOutOfRange".to_string());
        }
    }

    let p1_inf = Precision::is_negative_infinite(p1);
    let p2_inf = Precision::is_positive_infinite(p2);
    // `cxx:1132-1146`: the 3D points of the two ends.
    let point_at = |u: f64| {
        let p2d = pc.d0(u);
        surf.d0(p2d.x(), p2d.y())
    };
    let q1 = if p1_inf { None } else { Some(point_at(p1)) };
    let q2 = if p2_inf { None } else { Some(point_at(p2)) };

    // `cxx:1148`: `preci = BRepLib::Precision()`.
    let preci = CONFUSION;
    // `cxx:1152-1155`.
    let closed = match (&q1, &q2) {
        (Some(a), Some(b)) => a.distance(b) <= preci,
        _ => false,
    };

    // `cxx:1157-1268`: the vertex checks.
    let mut v1 = Some(v1);
    let mut v2 = Some(v2);
    if closed {
        match (v1.clone(), v2.clone()) {
            (None, None) => {
                let p = q1.clone().unwrap_or_else(|| point_at(p1));
                let b = TopoBuilder::new();
                let made = b.make_vertex(p, preci);
                v1 = Some(made.clone());
                v2 = Some(made);
            }
            (None, Some(other)) => v1 = Some(other),
            (Some(other), None) => v2 = Some(other),
            (Some(a), Some(b)) => {
                if !crate::topo_tools::is_same(&a.0, &b.0) {
                    return Err("BRepLib_DifferentPointsOnClosedCurve".to_string());
                }
                if let Some(q) = q1.as_ref() {
                    let tol = preci.max(BRepTool::vertex_tolerance(&a));
                    if q.distance(&BRepTool::vertex_point_world(&a)) > tol {
                        return Err("BRepLib_DifferentPointsOnClosedCurve".to_string());
                    }
                }
            }
        }
    } else {
        // `cxx:1229-1268`.
        if p1_inf {
            if v1.is_some() {
                return Err("BRepLib_PointWithInfiniteParameter".to_string());
            }
        } else if let Some(a) = v1.as_ref() {
            let tol = preci.max(BRepTool::vertex_tolerance(a));
            let dist = q1.as_ref().map_or(0.0, |q| q.distance(&BRepTool::vertex_point_world(a)));
            if dist > tol {
                return Err("BRepLib_DifferentsPointAndParameter".to_string());
            }
        } else {
            let b = TopoBuilder::new();
            v1 = Some(b.make_vertex(q1.clone().unwrap_or_else(|| point_at(p1)), preci));
        }

        if p2_inf {
            if v2.is_some() {
                return Err("BRepLib_PointWithInfiniteParameter".to_string());
            }
        } else if let Some(a) = v2.as_ref() {
            let tol = preci.max(BRepTool::vertex_tolerance(a));
            let dist = q2.as_ref().map_or(0.0, |q| q.distance(&BRepTool::vertex_point_world(a)));
            if dist > tol {
                return Err("BRepLib_DifferentsPointAndParameter".to_string());
            }
        } else {
            let b = TopoBuilder::new();
            v2 = Some(b.make_vertex(q2.clone().unwrap_or_else(|| point_at(p2)), preci));
        }
    }

    // `cxx:1270-1300`: `B.MakeEdge(E)`, `B.UpdateEdge(E, C, S, Location(), preci)`,
    // the two vertices, `B.Range(E, p1, p2)`.
    let reg = GeometryRegistry::global();
    let b = TopoBuilder::new();
    let mut e = Edge::new();
    let cos = CurveOnSurface::new(pc.clone(), surf.clone());
    let mut geom = EdgeGeom::new(Arc::new(cos), p1, p2);
    geom.tolerance = preci;
    reg.set_edge(&e.0, geom);
    reg.update_edge_pcurve_on_surface(&e.0, surf, Some(pc.clone()), preci);
    if let Some(v) = v1.as_ref() {
        b.add(&mut e.0, &v.0.oriented(Orientation::Forward));
    }
    if let Some(v) = v2.as_ref() {
        b.add(&mut e.0, &v.0.oriented(Orientation::Reversed));
    }
    // `cxx:1302-1305`.
    if reverse {
        e.0.set_orientation(Orientation::Reversed);
    }
    Ok(e)
}

