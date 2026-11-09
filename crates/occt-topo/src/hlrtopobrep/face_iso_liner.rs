//! `HLRTopoBRep_FaceIsoLiner` (`ModelingAlgorithms/TKHLR/HLRTopoBRep/HLRTopoBRep_FaceIsoLiner.hxx:36-91`,
//! `.cxx:52-489`): the isoparametric lines a face contributes when the caller
//! asks for them (`nbIso != 0`).
//!
//! `Perform` (`cxx:55-414`) hatches the face's parameter rectangle with the two
//! families of iso lines through `Geom2dHatch_Hatcher`, keeps the domains that
//! lie inside the face and stores one edge per domain through `MakeIsoLine`
//! (`cxx:467-489`). `MakeVertex` (`cxx:419-462`) resolves the domain endpoints
//! onto the existing edges of the face.
//!
//! Deviations:
//! * `Hatcher.AddElement(Geom2dAdaptor_Curve, Ori)` and
//!   `Hatcher.AddHatching(Geom2dAdaptor_Curve)` are
//!   [`crate::geom2d_hatch::Hatcher`]'s `Arc<dyn Curve2d>` overloads; the
//!   `Geom2d_TrimmedCurve` arm is [`Geom2dTrimmedCurve`].
//! * `BRep_Tool::CurveOnSurface(E, F, First, Last)` is
//!   [`crate::boptools_2d::curve_on_surface_range`].
//! * `TopoDS_Edge E; B.MakeEdge(E);` creates an edge with only
//!   `myRange = (-RealLast(), RealLast())`; the port's `Edge::new()` carries the
//!   same default (`tshape.rs:404-405`), and `BRep_TEdge::UpdateTolerance`
//!   becomes `set_edge_tolerance`.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::precision::{self, Precision};
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::Curve2d;
use occt_geom2d::Geom2dLine;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom2d_hatch::{Hatcher, Intersector};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::abs::ShapeType;
use crate::topo_tools::{edge_vertices, explorer_shapes};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::data::HlrTopoBRepData;

/// `IntersectorConfusion` (`cxx:43`).
const INTERSECTOR_CONFUSION: f64 = 1.0e-10;
/// `IntersectorTangency` (`cxx:44`).
const INTERSECTOR_TANGENCY: f64 = 1.0e-10;
/// `HatcherConfusion2d` (`cxx:45`).
const HATCHER_CONFUSION_2D: f64 = 1.0e-8;
/// `HatcherConfusion3d` (`cxx:46`).
const HATCHER_CONFUSION_3D: f64 = 1.0e-8;
/// `Infinite` (`cxx:47`).
const INFINITE: f64 = 100.;

/// `HLRTopoBRep_FaceIsoLiner::Perform(FI, F, DS, nbIsos)` (`cxx:55-414`).
pub fn perform(fi: i32, f: &Face, ds: &mut HlrTopoBRepData, nb_isos: i32) {
    let _ = fi; // `(void)FI;` (`cxx:59`)

    // `TopoDS_Face TF = F; TF.Orientation(TopAbs_FORWARD);` (`cxx:66-67`).
    let mut tf_shape = f.shape().clone();
    tf_shape.set_orientation(Orientation::Forward);
    let Some(tf) = Face::wrap(tf_shape) else {
        return;
    };

    let mut intersector = Intersector::new(INTERSECTOR_CONFUSION, INTERSECTOR_TANGENCY);
    let mut hatcher = Hatcher::new(
        intersector.clone(),
        HATCHER_CONFUSION_2D,
        HATCHER_CONFUSION_3D,
        true,
        false,
    );
    let _ = &mut intersector;

    // `BRepTools::UVBounds(TF, UMin, UMax, VMin, VMax)` (`cxx:74`).
    let (mut u_min, mut u_max, mut v_min, mut v_max) = BRepTool::uv_bounds(&tf);

    // `cxx:75-105`: replace the infinite bounds by a finite window.
    let inf_u_min = Precision::is_negative_infinite(u_min);
    let inf_u_max = Precision::is_positive_infinite(u_max);
    let inf_v_min = Precision::is_negative_infinite(v_min);
    let inf_v_max = Precision::is_positive_infinite(v_max);

    if inf_u_min && inf_u_max {
        u_min = -INFINITE;
        u_max = INFINITE;
    } else if inf_u_min {
        u_min = u_max - INFINITE;
    } else if inf_u_max {
        u_max = u_min + INFINITE;
    }
    if inf_v_min && inf_v_max {
        v_min = -INFINITE;
        v_max = INFINITE;
    } else if inf_v_min {
        v_min = v_max - INFINITE;
    } else if inf_v_max {
        v_max = v_min + INFINITE;
    }

    // `for (ExpEdges.Init(TF, TopAbs_EDGE); ...) ne++;` (`cxx:108-113`).
    let face_edges = explorer_shapes(tf.shape(), ShapeType::Edge);
    let mut ne = face_edges.len();
    // `if (DS.FaceHasIntL(TF)) { for (...) ne++; }` (`cxx:115-125`).
    let int_l: Vec<TopoShape> = if ds.face_has_int_l(&tf) {
        ds.face_int_l(&tf).to_vec()
    } else {
        Vec::new()
    };
    ne += int_l.len();

    // `NCollection_Array1<TopoDS_Shape> SH(1, ne); NCollection_Array1<bool> IL(1, ne);`
    // (`cxx:127-128`). 1-based slots, `None` where OCCT would have a null shape.
    let mut sh: Vec<Option<TopoShape>> = vec![None; ne + 1];
    let mut il: Vec<bool> = vec![false; ne + 1];

    // `cxx:130-153`: the real edges of the face.
    for e_shape in &face_edges {
        let Some(new_e) = Edge::wrap(e_shape.clone()) else {
            continue;
        };
        let Some((pc, u1, u2)) = crate::boptools_2d::curve_on_surface_range(&new_e, &tf) else {
            continue;
        };
        let ind_e = if (pc.first_parameter() - u1).abs() <= precision::PCONFUSION
            && (pc.last_parameter() - u2).abs() <= precision::PCONFUSION
        {
            // `IndE = Hatcher.AddElement(PC, newE.Orientation());` (`cxx:137`).
            hatcher.add_element(pc, new_e.shape().orientation())
        } else {
            // `Geom2d_TrimmedCurve TPC(PC, U1, U2); IndE =
            // Hatcher.AddElement(TPC, newE.Orientation());` (`cxx:141-143`).
            let tpc: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(pc, u1, u2));
            hatcher.add_element(tpc, new_e.shape().orientation())
        };
        if ind_e < 1 || ind_e as usize > ne {
            continue;
        }
        sh[ind_e as usize] = Some(new_e.0.clone());
        // `IL(IndE) = DS.IsOutLFaceEdge(TF, newE);` (`cxx:145-151`).
        il[ind_e as usize] = ds.is_out_l_face_edge(&tf, &new_e);
    }

    // `cxx:155-177`: the internal outlines already built on the face.
    for e_shape in &int_l {
        let Some(new_e) = Edge::wrap(e_shape.clone()) else {
            continue;
        };
        let Some((pc, u1, u2)) = crate::boptools_2d::curve_on_surface_range(&new_e, &tf) else {
            continue;
        };
        let ind_e = if (pc.first_parameter() - u1).abs() <= precision::PCONFUSION
            && (pc.last_parameter() - u2).abs() <= precision::PCONFUSION
        {
            hatcher.add_element(pc, Orientation::Internal)
        } else {
            let tpc: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(pc, u1, u2));
            hatcher.add_element(tpc, Orientation::Internal)
        };
        if ind_e < 1 || ind_e as usize > ne {
            continue;
        }
        sh[ind_e as usize] = Some(new_e.0.clone());
        il[ind_e as usize] = true; // `IL(IndE) = true;` (`cxx:175`)
    }

    // `BRepAdaptor_Surface Surface(TF); double Tolerance = Tolerance(TF);`
    // (`cxx:184-185`).
    let surface = SurfaceAdaptor::of_face_restricted(&tf);
    let tolerance = BRepTool::face_tolerance(&tf);

    // `double DeltaU = |UMax - UMin|; DeltaV = |VMax - VMin|;
    // Confusion = min(DeltaU, DeltaV) * HatcherConfusion3d;
    // Hatcher.Confusion3d(Confusion);` (`cxx:188-192`).
    let delta_u = (u_max - u_min).abs();
    let delta_v = (v_max - v_min).abs();
    let confusion = delta_u.min(delta_v) * HATCHER_CONFUSION_3D;
    hatcher.set_confusion3d(confusion);

    // `cxx:197-302`: the U iso curves.
    let step_u = delta_u / (nb_isos as f64);
    if step_u > confusion {
        let mut u_prm = u_min + step_u / 2.;
        let dir = GpDir2d::new(0.0, 1.0).expect("unit gp_Dir2d");
        for _ in 1..=nb_isos {
            let iso = Geom2dLine::from_pnt_dir(GpPnt2d::new(u_prm, 0.), dir);
            let iso_arc: Arc<dyn Curve2d> = Arc::new(iso.clone());
            let ind_h = hatcher.add_hatching(Arc::clone(&iso_arc));
            hatcher.trim_ind(ind_h);
            if hatcher.trim_done(ind_h) && !hatcher.trim_failed(ind_h) {
                hatcher.compute_domains_ind(ind_h);
            }
            if !hatcher.is_done_ind(ind_h) {
                // `cxx:215-240`: the hatcher rejected the line.
                hatcher.rem_hatching(ind_h);
                u_prm += step_u;
                continue;
            }
            let nb_dom = hatcher.nb_domains(ind_h);
            for idom in 1..=nb_dom {
                let domains = read_domain(&hatcher, ind_h, idom as i32, v_min, v_max);
                iso_u_domain(&tf, ds, &surface, &iso, tolerance, domains, &sh, &il);
            }
            hatcher.rem_hatching(ind_h);
            u_prm += step_u;
        }
    }

    // `cxx:307-414`: the V iso curves.
    let step_v = delta_v / (nb_isos as f64);
    if step_v > confusion {
        let mut v_prm = v_min + step_v / 2.;
        let dir = GpDir2d::new(1.0, 0.0).expect("unit gp_Dir2d");
        for _ in 1..=nb_isos {
            let iso = Geom2dLine::from_pnt_dir(GpPnt2d::new(0., v_prm), dir);
            let iso_arc: Arc<dyn Curve2d> = Arc::new(iso.clone());
            let ind_h = hatcher.add_hatching(Arc::clone(&iso_arc));
            hatcher.trim_ind(ind_h);
            if hatcher.trim_done(ind_h) && !hatcher.trim_failed(ind_h) {
                hatcher.compute_domains_ind(ind_h);
            }
            if !hatcher.is_done_ind(ind_h) {
                hatcher.rem_hatching(ind_h);
                v_prm += step_v;
                continue;
            }
            let nb_dom = hatcher.nb_domains(ind_h);
            for idom in 1..=nb_dom {
                let domains = read_domain(&hatcher, ind_h, idom as i32, v_min, v_max);
                iso_v_domain(&tf, ds, &surface, &iso, tolerance, domains, &sh, &il);
            }
            hatcher.rem_hatching(ind_h);
            v_prm += step_v;
        }
    }
}

/// The `HatchGen_Domain` data the two iso-pass arms read, copied out of the
/// hatcher so the borrow can end before `RemHatching`.
struct DomainInfo {
    has_first: bool,
    first_param: f64,
    first_points: Vec<(i32, f64)>,
    has_second: bool,
    second_param: f64,
    second_points: Vec<(i32, f64)>,
}

/// `Hatcher.Domain(IndH, IDom)` (`cxx:249`, `:359`).
fn read_domain(hatcher: &Hatcher, ind_h: i32, idom: i32, v_min: f64, v_max: f64) -> DomainInfo {
    let dom = hatcher.domain(ind_h, idom);
    let has_first = dom.has_first_point();
    let first_param = if has_first {
        dom.first_point().parameter()
    } else {
        v_min - INFINITE
    };
    let mut first_points = Vec::new();
    if has_first {
        let pnt_h = dom.first_point();
        for ipnt in 1..=pnt_h.nb_points() {
            let pnt_e = pnt_h.point(ipnt);
            first_points.push((pnt_e.index(), pnt_e.parameter()));
        }
    }
    let has_second = dom.has_second_point();
    let second_param = if has_second {
        dom.second_point().parameter()
    } else {
        v_max + INFINITE
    };
    let mut second_points = Vec::new();
    if has_second {
        let pnt_h = dom.second_point();
        for ipnt in 1..=pnt_h.nb_points() {
            let pnt_e = pnt_h.point(ipnt);
            second_points.push((pnt_e.index(), pnt_e.parameter()));
        }
    }
    DomainInfo {
        has_first,
        first_param,
        first_points,
        has_second,
        second_param,
        second_points,
    }
}

/// The body of the `IDom` loop of the U iso pass (`cxx:246-298`).
#[allow(clippy::too_many_arguments)]
fn iso_u_domain(
    tf: &Face,
    ds: &mut HlrTopoBRepData,
    surface: &SurfaceAdaptor,
    iso: &Geom2dLine,
    tolerance: f64,
    dom: DomainInfo,
    sh: &[Option<TopoShape>],
    il: &[bool],
) {
    let u11 = dom.first_param;
    let u21 = dom.second_param;
    let p1 = iso.d0(u11);
    let p1 = surface
        .d0(p1.x(), p1.y())
        .unwrap_or_else(|| occt_core::gp::GpPnt::new(0.0, 0.0, 0.0));
    let p2 = iso.d0(u21);
    let p2 = surface
        .d0(p2.x(), p2.y())
        .unwrap_or_else(|| occt_core::gp::GpPnt::new(0.0, 0.0, 0.0));

    let mut v1u: Option<Vertex> = None;
    let mut v2u: Option<Vertex> = None;
    if dom.has_first {
        for (index, parameter) in &dom.first_points {
            if let Some(v) = make_vertex_at(sh, *index, &p1, *parameter, tolerance, ds) {
                if il.get(*index as usize).copied().unwrap_or(false) {
                    ds.add_out_v(&v);
                }
                v1u = Some(v);
            }
        }
    }
    if dom.has_second {
        for (index, parameter) in &dom.second_points {
            if let Some(v) = make_vertex_at(sh, *index, &p2, *parameter, tolerance, ds) {
                if il.get(*index as usize).copied().unwrap_or(false) {
                    ds.add_out_v(&v);
                }
                v2u = Some(v);
            }
        }
    }
    if let (Some(v1), Some(v2)) = (v1u, v2u) {
        make_iso_line(tf, iso, &v1, &v2, u11, u21, tolerance, ds);
    }
}

/// The body of the `IDom` loop of the V iso pass (`cxx:356-408`).
#[allow(clippy::too_many_arguments)]
fn iso_v_domain(
    tf: &Face,
    ds: &mut HlrTopoBRepData,
    surface: &SurfaceAdaptor,
    iso: &Geom2dLine,
    tolerance: f64,
    dom: DomainInfo,
    sh: &[Option<TopoShape>],
    il: &[bool],
) {
    let u12 = dom.first_param;
    let u22 = dom.second_param;
    let p1 = iso.d0(u12);
    let p1 = surface
        .d0(p1.x(), p1.y())
        .unwrap_or_else(|| occt_core::gp::GpPnt::new(0.0, 0.0, 0.0));
    let p2 = iso.d0(u22);
    let p2 = surface
        .d0(p2.x(), p2.y())
        .unwrap_or_else(|| occt_core::gp::GpPnt::new(0.0, 0.0, 0.0));

    let mut v1v: Option<Vertex> = None;
    let mut v2v: Option<Vertex> = None;
    if dom.has_first {
        for (index, parameter) in &dom.first_points {
            if let Some(v) = make_vertex_at(sh, *index, &p1, *parameter, tolerance, ds) {
                if il.get(*index as usize).copied().unwrap_or(false) {
                    ds.add_out_v(&v);
                }
                v1v = Some(v);
            }
        }
    }
    if dom.has_second {
        for (index, parameter) in &dom.second_points {
            if let Some(v) = make_vertex_at(sh, *index, &p2, *parameter, tolerance, ds) {
                if il.get(*index as usize).copied().unwrap_or(false) {
                    ds.add_out_v(&v);
                }
                v2v = Some(v);
            }
        }
    }
    if let (Some(v1), Some(v2)) = (v1v, v2v) {
        make_iso_line(tf, iso, &v1, &v2, u12, u22, tolerance, ds);
    }
}

/// `HLRTopoBRep_FaceIsoLiner::MakeVertex(TopoDS::Edge(SH(Index)), P, Par, Tol,
/// DS)` (`cxx:265-269`, `:281-285`, `:373-377`, `:393-397`).
fn make_vertex_at(
    sh: &[Option<TopoShape>],
    index: i32,
    p: &occt_core::gp::GpPnt,
    parameter: f64,
    tolerance: f64,
    ds: &mut HlrTopoBRepData,
) -> Option<Vertex> {
    if index < 1 {
        return None;
    }
    let shape = sh.get(index as usize)?.as_ref()?;
    let e = Edge::wrap(shape.clone())?;
    Some(make_vertex(&e, p, parameter, tolerance, ds))
}

/// `HLRTopoBRep_FaceIsoLiner::MakeVertex(E, P, Par, Tol, DS)` (`cxx:419-462`).
fn make_vertex(
    e: &Edge,
    p: &occt_core::gp::GpPnt,
    par: f64,
    tol: f64,
    ds: &mut HlrTopoBRepData,
) -> Vertex {
    let maker = TopoBuilder;
    // `TopExp::Vertices(E, VF, VL);` (`cxx:427`).
    let (vf, vl) = edge_vertices(e);
    if let Some(vf) = &vf {
        // `if (P.IsEqual(Pnt(VF), Tolerance(VF))) return VF;` (`cxx:428-431`).
        if p.is_equal_with_tolerance(&BRepTool::vertex_point(vf), BRepTool::vertex_tolerance(vf)) {
            return vf.clone();
        }
    }
    if let Some(vl) = &vl {
        // `if (P.IsEqual(Pnt(VL), Tolerance(VL))) return VL;` (`cxx:432-435`).
        if p.is_equal_with_tolerance(&BRepTool::vertex_point(vl), BRepTool::vertex_tolerance(vl)) {
            return vl.clone();
        }
    }

    // `for (DS.InitVertex(E); DS.MoreVertex(); DS.NextVertex())` (`cxx:437`).
    let mut found: Option<Vertex> = None;
    let mut inserted = false;
    ds.init_vertex(e);
    while ds.more_vertex() {
        let Some(cur_v) = ds.vertex().cloned() else {
            break;
        };
        let cur_p = ds.parameter();
        let Some(cur_v) = Vertex::wrap(cur_v) else {
            break;
        };
        if p.is_equal_with_tolerance(&BRepTool::vertex_point(&cur_v), BRepTool::vertex_tolerance(&cur_v))
        {
            found = Some(cur_v);
            break;
        } else if par < cur_p {
            // `B.MakeVertex(V, P, Tol); V.Orientation(INTERNAL); DS.InsertBefore(V, Par);`
            // (`cxx:446-451`).
            let mut v = maker.make_vertex(*p, tol);
            v.0.set_orientation(Orientation::Internal);
            ds.insert_before(&v.0, par);
            found = Some(v);
            inserted = true;
            break;
        }
        ds.next_vertex();
    }
    if !inserted && !ds.more_vertex() {
        // `B.MakeVertex(V, P, Tol); V.Orientation(INTERNAL); DS.Append(V, Par);`
        // (`cxx:455-459`).
        let mut v = maker.make_vertex(*p, tol);
        v.0.set_orientation(Orientation::Internal);
        ds.append(&v.0, par);
        found = Some(v);
    }
    found.expect("HLRTopoBRep_FaceIsoLiner::MakeVertex resolved a vertex")
}

/// `HLRTopoBRep_FaceIsoLiner::MakeIsoLine(F, Iso, V1, V2, U1, U2, Tol, DS)`
/// (`cxx:467-489`).
#[allow(clippy::too_many_arguments)]
fn make_iso_line(
    f: &Face,
    iso: &Geom2dLine,
    v1: &Vertex,
    v2: &Vertex,
    u1: f64,
    u2: f64,
    tol: f64,
    ds: &mut HlrTopoBRepData,
) {
    let maker = TopoBuilder;
    // `TopoDS_Edge E; E.Orientation(TopAbs_INTERNAL);` (`cxx:478-479`).
    let mut e = Edge::new();
    e.0.set_orientation(Orientation::Internal);
    // `V1.Orientation(FORWARD); V2.Orientation(REVERSED);` (`cxx:480-481`).
    let mut v1 = v1.clone();
    let mut v2 = v2.clone();
    v1.0.set_orientation(Orientation::Forward);
    v2.0.set_orientation(Orientation::Reversed);

    // `B.MakeEdge(E); B.UpdateEdge(E, Iso, F, Tol);` (`cxx:482-483`):
    // `BRep_Builder::UpdateEdge(E, C, F, Tol)` binds the pcurve to
    // `BRep_Tool::Surface(F)` and only raises the tolerance.
    if let Some(surface) = BRepTool::face_surface(f) {
        GeometryRegistry::global().update_edge_pcurve_on_surface(
            &e.0,
            &surface,
            Some(Arc::new(iso.clone()) as Arc<dyn Curve2d>),
            tol,
        );
    }

    // `B.Add(E, V1); B.UpdateVertex(V1, U1, E, Tol);` (`cxx:484-485`).
    maker.add(&mut e.0, &v1.0);
    GeometryRegistry::global().set_vertex_parameter_on_edge(&e.0, &v1.0, u1);
    // `B.Add(E, V2); B.UpdateVertex(V2, U2, E, Tol);` (`cxx:486-487`).
    maker.add(&mut e.0, &v2.0);
    GeometryRegistry::global().set_vertex_parameter_on_edge(&e.0, &v2.0, u2);

    // `DS.AddIsoL(F).Append(E);` (`cxx:488`).
    ds.add_iso_l(f).push(e.0);
}
