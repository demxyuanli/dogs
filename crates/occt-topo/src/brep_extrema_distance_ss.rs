//! `BRepExtrema_DistanceSS` — the distance between one pair of sub-shapes.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepExtrema/BRepExtrema_DistanceSS.cxx`
//! (1206 lines) and `.hxx:32-133` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(S1, S2, B1, B2, DstRef, Deflection, Flag, Algo)` | `hxx:50-64` | [`DistanceSS::new`] |
//! | `Perform(S1, S2, B1, B2)` type switch | `cxx:656-775` | [`DistanceSS::perform`] |
//! | `Perform(V, V)` | `cxx:780-802` | [`DistanceSS::perform_vertex_vertex`] |
//! | `Perform(V, E)` | `cxx:806-864` | [`DistanceSS::perform_vertex_edge`] |
//! | `Perform(V, F)` | `cxx:868-925` | [`DistanceSS::perform_vertex_face`] |
//! | `Perform(E, E)` | `cxx:929-1023` | [`DistanceSS::perform_edge_edge`] |
//! | `Perform(E, F)` | `cxx:1026-1187` | [`DistanceSS::perform_edge_face`] |
//! | `Perform(F, F)` | `cxx:1191-1206` | [`DistanceSS::perform_face_face`] |
//! | `TRI_SOLUTION` | `cxx:60-75` | [`tri_solution`] |
//! | `MIN_SOLUTION` | `cxx:79-97` | [`min_solution`] |
//! | `TRIM_INFINIT_EDGE` | `cxx:102-245` | [`trim_infinit_edge`] |
//! | `TRIM_INFINIT_FACE` | `cxx:249-484` | [`trim_infinit_face`] |
//! | `PERFORM_C0` | `cxx:488-630` | [`perform_c0`] |
//! | `isOnBoundary` | `cxx:637-653` | [`is_on_boundary`] |
//!
//! # Translation notes
//!
//! * `BRepClass_FaceClassifier` is [`crate::brep_class::FaceClassifier`]; the
//!   accepted states are `IN` only in the vertex-face (`cxx:907`),
//!   edge-face (`cxx:1080`) and face-face arms, exactly as the `.cxx` reads.
//! * `Extrema_ExtPC` is [`occt_geom::extrema_pc::extrema_ext_pc_range`] plus
//!   its `IsDone()`-carrying twin; the wrapper class lives in
//!   [`super::brep_extrema_ext_pc`].
//! * `GeomAdaptor_Curve::NbIntervals(GeomAbs_C1)` / `Intervals(...)` are
//!   [`occt_geom::extrema_pc::extrema_nb_intervals`] /
//!   [`occt_geom::extrema_pc::extrema_intervals`] (`s = 1` is `GeomAbs_C1`).
//! * `BRep_Tool::Curve(E, F, L)` is [`BRepTool::edge_curve`] +
//!   [`BRepTool::edge_parameters`].
//! * `BRepBuilderAPI_MakeVertex(P)` is
//!   [`crate::builder::TopoBuilder::make_vertex`] with
//!   `Precision::Confusion()` (the constructor's default tolerance).
//! * `TRIM_INFINIT_EDGE` rebuilds an edge with `BRepBuilderAPI_MakeEdge(new
//!   Geom_TrimmedCurve(pCurv, Umin, Umax))` (`cxx:240-241`); the port uses
//!   [`crate::builder::TopoBuilder::make_edge`] over a
//!   [`occt_geom::trimmed::GeomTrimmedCurve`].
//! * `TRIM_INFINIT_FACE` rebuilds a face with `BRepBuilderAPI_MakeFace(new
//!   Geom_RectangularTrimmedSurface(...), Precision::Confusion())`
//!   (`cxx:471-472`); the port uses
//!   [`crate::builder::TopoBuilder::make_face_with_tolerance`].
//! * The `NCollection_Sequence` pairs stay two parallel `Vec`s, like OCCT's
//!   two sequences.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::Precision;
use occt_geom::curve::Curve;
use occt_geom::extrema_surf::ExtremaPOnSurf;
use occt_geom::{extrema_pc, Surface};

use crate::brep_class::FaceClassifier;
use crate::fclass2d::FaceState;
use crate::brep_bnd_lib;
use crate::brep_tool::BRepTool;
use crate::brep_extrema_ext_cc::ExtCC;
use crate::brep_extrema_ext_cf::ExtCF;
use crate::brep_extrema_ext_ff::ExtFF;
use crate::brep_extrema_ext_pc::ExtPC;
use crate::brep_extrema_ext_pf::ExtPF;
use crate::brep_extrema_solution_elem::{SolutionElem, SupportType};
use crate::builder::TopoBuilder;
use crate::abs::ShapeType;
use crate::shape::{Edge, Face, TopoShape, Vertex};

use occt_geom::extrema_surf::{ExtremaExtAlgo, ExtremaExtFlag};

/// `TRI_SOLUTION` (`BRepExtrema_DistanceSS.cxx:60-75`): true when no already
/// recorded solution is within `Precision::Confusion()` of `pt`.
pub fn tri_solution(seq_sol: &[SolutionElem], pt: &GpPnt) -> bool {
    for s in seq_sol {
        if s.point().distance(pt) <= Precision::CONFUSION {
            return false;
        }
    }
    true
}

/// `MIN_SOLUTION` (`cxx:79-97`): keep the pairs whose distance equals `dst_ref`
/// within `eps`.
fn min_solution(
    seq_sol1: &[SolutionElem],
    seq_sol2: &[SolutionElem],
    dst_ref: f64,
    eps: f64,
) -> (Vec<SolutionElem>, Vec<SolutionElem>) {
    let mut out1 = Vec::new();
    let mut out2 = Vec::new();
    for (s1, s2) in seq_sol1.iter().zip(seq_sol2.iter()) {
        if (s1.dist() - dst_ref).abs() < eps {
            out1.push(s1.clone());
            out2.push(s2.clone());
        }
    }
    (out1, out2)
}

/// `isOnBoundary` (`cxx:637-653`): true when `(the_sol, the_param)` is the
/// edge's vertex at that parameter, within the vertex tolerance.
pub fn is_on_boundary(the_edge: &Edge, the_sol: &GpPnt, the_param: f64, the_ptol: f64) -> bool {
    let kids = the_edge.0.tshape.read().unwrap().children.clone();
    for k in kids {
        if k.shape_type() != ShapeType::Vertex {
            continue;
        }
        let v = match Vertex::wrap(k.clone()) {
            Some(v) => v,
            None => continue,
        };
        let v_param = BRepTool::parameter_of_vertex(&the_edge.0, &k);
        if (v_param - the_param).abs() < the_ptol
            && BRepTool::vertex_point(&v).distance(the_sol) < BRepTool::vertex_tolerance(&v)
        {
            return true;
        }
    }
    false
}

/// `TRIM_INFINIT_EDGE` (`cxx:102-245`).
///
/// Returns `(res_edge, is_trim1, is_trim2)`; `res_edge` defaults to `s2`
/// (`cxx:113`).
pub fn trim_infinit_edge(
    s1: &Edge,
    s2: &Edge,
) -> (Option<Edge>, bool, bool) {
    let mut b_is_trim1 = false;
    let mut b_is_trim2 = false;
    // `cxx:108-111`.
    if BRepTool::is_degenerated(s1) || BRepTool::is_degenerated(s2) {
        return (None, b_is_trim1, b_is_trim2);
    }
    // `cxx:114-115`.
    let (a_first1, a_last1) = BRepTool::edge_parameters(s1);
    let (a_first2, a_last2) = BRepTool::edge_parameters(s2);
    let p_curv1 = BRepTool::edge_curve(s1);
    let p_curv2 = BRepTool::edge_curve(s2);

    // `cxx:117-121`.
    if Precision::is_infinite(a_first1)
        && Precision::is_infinite(a_last1)
        && Precision::is_infinite(a_first2)
        && Precision::is_infinite(a_last2)
    {
        return (None, b_is_trim1, b_is_trim2);
    }

    // `cxx:123-126`.
    let mut u_min = 0.0;
    let mut u_max = 0.0;
    let mut b_umin = false;
    let mut b_umax = false;

    // `cxx:128-158`: pick the infinite curve and its known-bound side.
    let mut p_curv: Option<Arc<dyn Curve>> = None;
    if let Some(c) = &p_curv1 {
        if Precision::is_infinite(a_first1) || Precision::is_infinite(a_last1) {
            p_curv = Some(c.clone());
            b_is_trim1 = true;
            if !Precision::is_infinite(a_first1) {
                b_umin = true;
                u_min = a_first1;
            } else if !Precision::is_infinite(a_last1) {
                b_umax = true;
                u_max = a_last1;
            }
        }
    }
    if p_curv.is_none() {
        if let Some(c) = &p_curv2 {
            if Precision::is_infinite(a_first2) || Precision::is_infinite(a_last2) {
                p_curv = Some(c.clone());
                b_is_trim2 = true;
                if !Precision::is_infinite(a_first2) {
                    b_umin = true;
                    u_min = a_first2;
                } else if !Precision::is_infinite(a_last2) {
                    b_umax = true;
                    u_max = a_last2;
                }
            }
        }
    }

    let p_curv = match p_curv {
        Some(c) => c,
        // `cxx:159`: neither edge was infinite.
        None => return (None, b_is_trim1, b_is_trim2),
    };

    // `cxx:161-181`: the other edge's bounding box, then its 8 corners.
    let other = if b_is_trim1 { s2 } else { s1 };
    let mut a_edge_box = BndBox::new();
    brep_bnd_lib::add(&other.0, &mut a_edge_box, false);
    let (xmin, ymin, zmin, xmax, ymax, zmax) = match a_edge_box.get() {
        Some(v) => v,
        None => return (None, b_is_trim1, b_is_trim2),
    };
    let corners = [
        GpPnt::new(xmin, ymin, zmin),
        GpPnt::new(xmin, ymax, zmin),
        GpPnt::new(xmin, ymax, zmax),
        GpPnt::new(xmin, ymin, zmax),
        GpPnt::new(xmax, ymax, zmin),
        GpPnt::new(xmax, ymax, zmax),
        GpPnt::new(xmax, ymin, zmax),
        GpPnt::new(xmax, ymin, zmin),
    ];

    // `cxx:183-198`: `GeomAPI_ProjectPointOnCurve`'s `LowerDistanceParameter`
    // for each corner.
    let mut arr_u = [0.0f64; 8];
    for (i, c) in corners.iter().enumerate() {
        arr_u[i] = occt_geom::geom_api::project_point_on_curve(p_curv.as_ref(), c, Precision::CONFUSION)
            .map(|p| p.parameter)
            .unwrap_or(0.0);
    }

    // `cxx:200-207`.
    if !b_umin {
        u_min = arr_u[0];
    }
    if !b_umax {
        u_max = arr_u[0];
    }
    // `cxx:209-223`.
    for a_u in arr_u {
        if a_u < u_min {
            u_min = a_u;
        } else if a_u > u_max {
            u_max = a_u;
        }
    }

    // `cxx:225-233`.
    let mut tol = Precision::CONFUSION;
    if b_is_trim1 {
        tol = BRepTool::edge_tolerance(s1);
    } else if b_is_trim2 {
        tol = BRepTool::edge_tolerance(s2);
    }

    // `cxx:235-239`.
    let eps_u = p_curv.resolution(3.0 * tol);
    if (u_min - u_max).abs() < eps_u {
        u_min -= eps_u;
        u_max += eps_u;
    }

    // `cxx:241-242`: `new Geom_TrimmedCurve(pCurv, Umin, Umax)` then
    // `BRepBuilderAPI_MakeEdge(result)`.
    let trimmed = occt_geom::trimmed::GeomTrimmedCurve::new(p_curv, u_min, u_max);
    let b = TopoBuilder::new();
    let edge = b.make_edge(Arc::new(trimmed), u_min, u_max);
    (Some(edge), b_is_trim1, b_is_trim2)
}

/// `TRIM_INFINIT_FACE` (`cxx:249-484`).
///
/// Returns `(res_face, is_infinite)`; `res_face` is the second face argument
/// when the pair has no infinite face (`cxx:284`).
pub fn trim_infinit_face(s1: &TopoShape, s2: &TopoShape) -> (Option<Face>, bool) {
    // `cxx:252`.
    let mut b_is_infinit = false;

    // `cxx:254-283`: orient the (edge, face) pair.
    let (a_e, a_f) = match (s1.shape_type(), s2.shape_type()) {
        (ShapeType::Edge, ShapeType::Face) => {
            let e = match Edge::wrap(s1.clone()) {
                Some(e) => e,
                None => return (None, b_is_infinit),
            };
            // `cxx:264-267`.
            if BRepTool::is_degenerated(&e) {
                return (None, b_is_infinit);
            }
            match Face::wrap(s2.clone()) {
                Some(f) => (e, f),
                None => return (None, b_is_infinit),
            }
        }
        (ShapeType::Face, ShapeType::Edge) => {
            let e = match Edge::wrap(s2.clone()) {
                Some(e) => e,
                None => return (None, b_is_infinit),
            };
            // `cxx:272-275`.
            if BRepTool::is_degenerated(&e) {
                return (None, b_is_infinit);
            }
            match Face::wrap(s1.clone()) {
                Some(f) => (e, f),
                None => return (None, b_is_infinit),
            }
        }
        // `cxx:278-283`.
        _ => return (None, false),
    };

    // `cxx:285-287`.
    let p_surf = match BRepTool::face_surface(&a_f) {
        Some(s) => s,
        None => return (None, b_is_infinit),
    };
    // `cxx:288`.
    let b_restrict = BRepTool::natural_restriction(&a_f);

    // `cxx:290-293`.
    let mut u_min = f64::MAX;
    let mut u_max = f64::MIN;
    let mut v_min = f64::MAX;
    let mut v_max = f64::MIN;
    let (mut b_umin, mut b_umax, mut b_vmin, mut b_vmax) = (false, false, false, false);
    let mut b_is_trim = false;

    // `cxx:295-337`.
    if b_restrict {
        let (u1, u2) = p_surf.u_range();
        let (v1, v2) = p_surf.v_range();
        if Precision::is_infinite(u1) {
            b_is_trim = true;
        } else {
            u_min = u1;
            b_umin = true;
        }
        if Precision::is_infinite(u2) {
            b_is_trim = true;
        } else {
            u_max = u2;
            b_umax = true;
        }
        if Precision::is_infinite(v1) {
            b_is_trim = true;
        } else {
            v_min = v1;
            b_vmin = true;
        }
        if Precision::is_infinite(v2) {
            b_is_trim = true;
        } else {
            v_max = v2;
            b_vmax = true;
        }
    } else {
        // `cxx:340-347`.
        let (u1, u2, v1, v2) = BRepTool::uv_bounds(&a_f);
        if Precision::is_infinite(u1)
            && Precision::is_infinite(u2)
            && Precision::is_infinite(v1)
            && Precision::is_infinite(v2)
        {
            b_is_trim = true;
        }
    }

    if !b_is_trim {
        // `cxx:478-482`.
        return (Some(a_f), false);
    }

    // `cxx:351-355`.
    let mut a_edge_box = BndBox::new();
    brep_bnd_lib::add(&a_e.0, &mut a_edge_box, false);
    if a_edge_box.is_whole() {
        return (None, b_is_infinit);
    }
    let (xmin, ymin, zmin, xmax, ymax, zmax) = match a_edge_box.get() {
        Some(v) => v,
        None => return (None, b_is_infinit),
    };
    let corners = [
        GpPnt::new(xmin, ymin, zmin),
        GpPnt::new(xmin, ymax, zmin),
        GpPnt::new(xmin, ymax, zmax),
        GpPnt::new(xmin, ymin, zmax),
        GpPnt::new(xmax, ymax, zmin),
        GpPnt::new(xmax, ymax, zmax),
        GpPnt::new(xmax, ymin, zmax),
        GpPnt::new(xmax, ymin, zmin),
    ];

    // `cxx:372-410`: `GeomAPI_ProjectPointOnSurf`'s
    // `LowerDistanceParameters`; failed projections leave `(0, 0)`.
    let mut arr_u = [0.0f64; 8];
    let mut arr_v = [0.0f64; 8];
    for (i, c) in corners.iter().enumerate() {
        if let Some(p) =
            occt_geom::geom_api::project_point_on_surface(p_surf.as_ref(), c, Precision::CONFUSION)
        {
            arr_u[i] = p.u;
            arr_v[i] = p.v;
        }
    }

    // `cxx:412-426`.
    if !b_umin {
        u_min = arr_u[0];
    }
    if !b_umax {
        u_max = arr_u[0];
    }
    if !b_vmin {
        v_min = arr_v[0];
    }
    if !b_vmax {
        v_max = arr_v[0];
    }
    // `cxx:428-453`.
    for i in 0..8 {
        let a_u = arr_u[i];
        if a_u < u_min {
            u_min = a_u;
        } else if a_u > u_max {
            u_max = a_u;
        }
        let a_v = arr_v[i];
        if a_v < v_min {
            v_min = a_v;
        } else if a_v > v_max {
            v_max = a_v;
        }
    }

    // `cxx:455-470`.
    let tol = BRepTool::face_tolerance(&a_f);
    let eps_u = occt_geom::approx_same_parameter::u_resolution(p_surf.as_ref(), 3.0 * tol);
    if (u_min - u_max).abs() < eps_u {
        u_min -= eps_u;
        u_max += eps_u;
    }
    let eps_v = occt_geom::approx_same_parameter::v_resolution(p_surf.as_ref(), 3.0 * tol);
    if (v_min - v_max).abs() < eps_v {
        v_min -= eps_v;
        v_max += eps_v;
    }

    // `cxx:471-475`.
    let trimmed = match occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface::uv(
        p_surf, u_min, u_max, v_min, v_max,
    ) {
        Ok(t) => t,
        Err(_) => return (None, b_is_infinit),
    };
    let b = TopoBuilder::new();
    let face = b.make_face_with_tolerance(Arc::new(trimmed), &[], Precision::CONFUSION);
    b_is_infinit = true;
    (Some(face), b_is_infinit)
}

/// `PERFORM_C0` (`cxx:488-630`).
#[allow(clippy::too_many_arguments)]
pub fn perform_c0(
    s1: &Edge,
    s2: &Edge,
    seq_sol1: &mut Vec<SolutionElem>,
    seq_sol2: &mut Vec<SolutionElem>,
    dst_ref: f64,
    m_dst_ref: &mut f64,
    eps: f64,
) {
    // `cxx:495-498`.
    if BRepTool::is_degenerated(s1) || BRepTool::is_degenerated(s2) {
        return;
    }

    // `cxx:500-513`.
    for i_e in 0..2 {
        let (e, e_other) = if i_e == 0 { (s1, s2) } else { (s2, s1) };
        // `cxx:515-518`.
        let (a_first, a_last) = BRepTool::edge_parameters(e);
        let p_curv = match BRepTool::edge_curve(e) {
            Some(c) => c,
            None => continue,
        };
        let (a_f_other, a_l_other) = BRepTool::edge_parameters(e_other);
        let p_curv_other = match BRepTool::edge_curve(e_other) {
            Some(c) => c,
            None => continue,
        };

        // `cxx:520`.
        if p_curv.continuity() != 0 {
            continue;
        }

        // `cxx:522-536`.
        let eps_p = Precision::PCONFUSION;
        let arr_inter = extrema_pc::extrema_intervals(p_curv.as_ref(), 1);
        let arr_inter_other = extrema_pc::extrema_intervals(p_curv_other.as_ref(), 1);
        // `cxx:537-538`.
        let (u_deb, u_fin) = (a_f_other, a_l_other);

        // `cxx:541-628`.
        for &a_parameter in &arr_inter {
            let a_pnt = p_curv.d0(a_parameter);
            // `cxx:544`.
            let b = TopoBuilder::new();
            let v1 = b.make_vertex(a_pnt, Precision::CONFUSION);

            // `cxx:547-548`.
            let ext = ExtPC::new(&v1, e_other);
            let nb_extrema = if ext.is_done() { ext.nb_ext() } else { 0 };
            if nb_extrema > 0 {
                // `cxx:553-560`.
                let mut dstmin = ext.square_distance(1);
                for ii in 2..=nb_extrema {
                    let s_dst = ext.square_distance(ii);
                    if s_dst < dstmin {
                        dstmin = s_dst;
                    }
                }
                dstmin = dstmin.sqrt();

                // `cxx:563-604`.
                if dstmin < dst_ref - eps || (dstmin - dst_ref).abs() < eps {
                    for ii in 1..=nb_extrema {
                        if (dstmin - ext.square_distance(ii).sqrt()).abs() < eps {
                            let pt = ext.point(ii);
                            // `cxx:570-578`: the recorded point lies on
                            // `Eother`, so the already-known solutions are the
                            // mirror sequence.
                            let tri = if i_e == 0 {
                                tri_solution(seq_sol2, &pt)
                            } else {
                                tri_solution(seq_sol1, &pt)
                            };
                            if tri {
                                // `cxx:581-582`.
                                let t = ext.parameter(ii);
                                if (t - u_deb).abs() >= eps_p && (t - u_fin).abs() > eps_p {
                                    if *m_dst_ref > dstmin {
                                        *m_dst_ref = dstmin;
                                    }
                                    // `cxx:588-599`.
                                    let sol1 = SolutionElem::of_edge(
                                        dstmin,
                                        &a_pnt,
                                        SupportType::IsOnEdge,
                                        e,
                                        a_parameter,
                                    );
                                    let sol2 = SolutionElem::of_edge(
                                        dstmin,
                                        &pt,
                                        SupportType::IsOnEdge,
                                        e_other,
                                        t,
                                    );
                                    if i_e == 0 {
                                        seq_sol1.push(sol1);
                                        seq_sol2.push(sol2);
                                    } else {
                                        seq_sol1.push(sol2);
                                        seq_sol2.push(sol1);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // `cxx:606-627`.
            for &a_parameter_other in &arr_inter_other {
                let a_pnt_other = p_curv_other.d0(a_parameter_other);
                let dst = a_pnt.distance(&a_pnt_other);
                if dst < dst_ref - eps || (dst - dst_ref).abs() < eps {
                    if *m_dst_ref > dst {
                        *m_dst_ref = dst;
                    }
                    let sol1 = SolutionElem::of_edge(
                        dst,
                        &a_pnt,
                        SupportType::IsOnEdge,
                        e,
                        a_parameter,
                    );
                    let sol2 = SolutionElem::of_edge(
                        dst,
                        &a_pnt_other,
                        SupportType::IsOnEdge,
                        e_other,
                        a_parameter_other,
                    );
                    if i_e == 0 {
                        seq_sol1.push(sol1);
                        seq_sol2.push(sol2);
                    } else {
                        seq_sol1.push(sol2);
                        seq_sol2.push(sol1);
                    }
                }
            }
            let _ = a_first;
            let _ = a_last;
            let _ = a_f_other;
            let _ = a_l_other;
        }
    }
}

/// `BRepExtrema_DistanceSS` (`BRepExtrema_DistanceSS.hxx:32-133`).
pub struct DistanceSS {
    /// `myDstRef` (`hxx:126`).
    my_dst_ref: f64,
    /// `myModif` (`hxx:127`) — `IsDone()`.
    my_modif: bool,
    /// `myEps` (`hxx:128`).
    my_eps: f64,
    /// `mySeqSolShape1` (`hxx:130`).
    my_seq_sol_shape1: Vec<SolutionElem>,
    /// `mySeqSolShape2` (`hxx:131`).
    my_seq_sol_shape2: Vec<SolutionElem>,
    /// `myFlag` / `myAlgo` (`hxx:132-133`).
    my_flag: ExtremaExtFlag,
    my_algo: ExtremaExtAlgo,
}

impl DistanceSS {
    /// `BRepExtrema_DistanceSS(S1, S2, Box1, Box2, DstRef, Deflection, Flag,
    /// Algo)` (`hxx:50-64`): stores the reference distance and runs `Perform`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        the_s1: &TopoShape,
        the_s2: &TopoShape,
        the_box1: &BndBox,
        the_box2: &BndBox,
        the_dst_ref: f64,
        the_deflection: f64,
        the_ext_flag: ExtremaExtFlag,
        the_ext_algo: ExtremaExtAlgo,
    ) -> Self {
        let mut this = Self {
            my_dst_ref: the_dst_ref,
            my_modif: false,
            my_eps: the_deflection,
            my_seq_sol_shape1: Vec::new(),
            my_seq_sol_shape2: Vec::new(),
            my_flag: the_ext_flag,
            my_algo: the_ext_algo,
        };
        this.perform(the_s1, the_s2, the_box1, the_box2);
        this
    }

    /// `Perform(S1, S2, Box1, Box2)` (`cxx:656-775`).
    pub fn perform(
        &mut self,
        the_s1: &TopoShape,
        the_s2: &TopoShape,
        the_box1: &BndBox,
        the_box2: &BndBox,
    ) {
        // `cxx:660-662`.
        self.my_seq_sol_shape1.clear();
        self.my_seq_sol_shape2.clear();
        self.my_modif = false;

        // `cxx:664-671`: the box–box gap is the early exit. `Bnd_Box::Distance`
        // and `brep_extrema::shape_bbox_distance` compute the same axis-gap
        // norm for finite boxes.
        let a_bb_dist = crate::brep_extrema::shape_bbox_distance(the_box1, the_box2);
        if a_bb_dist - self.my_dst_ref > self.my_eps {
            return;
        }

        // `cxx:673-772`: the type switch. Every arm borrows the two sequences
        // mutably, so the arms are dispatched explicitly instead of through a
        // closure.
        match (the_s1.shape_type(), the_s2.shape_type()) {
            (ShapeType::Vertex, ShapeType::Vertex) => {
                let (v1, v2) = match (Vertex::wrap(the_s1.clone()), Vertex::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                self.perform_vertex_vertex(&v1, &v2);
            }
            (ShapeType::Vertex, ShapeType::Edge) => {
                let (v1, e2) = match (Vertex::wrap(the_s1.clone()), Edge::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                self.perform_vertex_edge(&v1, &e2);
            }
            (ShapeType::Vertex, ShapeType::Face) => {
                let (v1, f2) = match (Vertex::wrap(the_s1.clone()), Face::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                self.perform_vertex_face(&v1, &f2);
            }
            (ShapeType::Edge, ShapeType::Vertex) => {
                let (e1, v2) = match (Edge::wrap(the_s1.clone()), Vertex::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
                let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
                let mut dst = self.my_dst_ref;
                let mut modif = self.my_modif;
                vertex_edge_impl(
                    &v2,
                    &e1,
                    &mut seq2,
                    &mut seq1,
                    &mut dst,
                    self.my_eps,
                    &mut modif,
                );
                self.my_seq_sol_shape1 = seq1;
                self.my_seq_sol_shape2 = seq2;
                self.my_dst_ref = dst;
                self.my_modif = modif;
            }
            (ShapeType::Edge, ShapeType::Edge) => {
                let (e1, e2) = match (Edge::wrap(the_s1.clone()), Edge::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                // `cxx:710-723`: `TRIM_INFINIT_EDGE`.
                let (trim, trim1, trim2) = trim_infinit_edge(&e1, &e2);
                let (t1, t2) = match trim {
                    Some(t) => (
                        if trim1 { t.clone() } else { e1.clone() },
                        if trim2 { t } else { e2.clone() },
                    ),
                    None => (e1.clone(), e2.clone()),
                };
                self.perform_edge_edge(&t1, &t2);
            }
            (ShapeType::Edge, ShapeType::Face) => {
                let (e1, f2) = match (Edge::wrap(the_s1.clone()), Face::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                // `cxx:725-734`: `TRIM_INFINIT_FACE(aE1, aF2, ...)`.
                let (trim, is_inf) = trim_infinit_face(&e1.0, &f2.0);
                let f2 = if is_inf {
                    match trim {
                        Some(f) => f,
                        None => f2.clone(),
                    }
                } else {
                    f2.clone()
                };
                self.perform_edge_face(&e1, &f2);
            }
            (ShapeType::Face, ShapeType::Vertex) => {
                let (f1, v2) = match (Face::wrap(the_s1.clone()), Vertex::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
                let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
                let mut dst = self.my_dst_ref;
                let mut modif = self.my_modif;
                vertex_face_impl(
                    &v2,
                    &f1,
                    &mut seq2,
                    &mut seq1,
                    &mut dst,
                    self.my_eps,
                    self.my_flag,
                    self.my_algo,
                    &mut modif,
                );
                self.my_seq_sol_shape1 = seq1;
                self.my_seq_sol_shape2 = seq2;
                self.my_dst_ref = dst;
                self.my_modif = modif;
            }
            (ShapeType::Face, ShapeType::Edge) => {
                let (f1, e2) = match (Face::wrap(the_s1.clone()), Edge::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                // `cxx:753-759`: `TRIM_INFINIT_FACE(aF1, aE2, ...)`.
                let (trim, is_inf) = trim_infinit_face(&f1.0, &e2.0);
                let f1 = if is_inf {
                    match trim {
                        Some(f) => f,
                        None => f1.clone(),
                    }
                } else {
                    f1.clone()
                };
                let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
                let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
                let mut dst = self.my_dst_ref;
                let mut modif = self.my_modif;
                edge_face_impl(
                    &e2,
                    &f1,
                    &mut seq2,
                    &mut seq1,
                    &mut dst,
                    self.my_eps,
                    &mut modif,
                );
                self.my_seq_sol_shape1 = seq1;
                self.my_seq_sol_shape2 = seq2;
                self.my_dst_ref = dst;
                self.my_modif = modif;
            }
            (ShapeType::Face, ShapeType::Face) => {
                let (f1, f2) = match (Face::wrap(the_s1.clone()), Face::wrap(the_s2.clone())) {
                    (Some(a), Some(b)) => (a, b),
                    _ => return,
                };
                self.perform_face_face(&f1, &f2);
            }
            // `cxx:771-772`.
            _ => {}
        }
    }

    /// `IsDone()` (`hxx:69`): `myModif`.
    pub fn is_done(&self) -> bool {
        self.my_modif
    }

    /// `DistValue()` (`hxx:72`).
    pub fn dist_value(&self) -> f64 {
        self.my_dst_ref
    }

    /// `Seq1Value()` (`hxx:75`).
    pub fn seq1_value(&self) -> &[SolutionElem] {
        &self.my_seq_sol_shape1
    }

    /// `Seq2Value()` (`hxx:78`).
    pub fn seq2_value(&self) -> &[SolutionElem] {
        &self.my_seq_sol_shape2
    }

    /// `Perform(V1, V2, Seq1, Seq2)` (`cxx:780-802`).
    pub fn perform_vertex_vertex(&mut self, s1: &Vertex, s2: &Vertex) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        vertex_vertex_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }

    /// `Perform(V, E, Seq1, Seq2)` (`cxx:806-864`).
    pub fn perform_vertex_edge(&mut self, s1: &Vertex, s2: &Edge) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        vertex_edge_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }

    /// `Perform(V, F, Seq1, Seq2)` (`cxx:868-925`).
    pub fn perform_vertex_face(&mut self, s1: &Vertex, s2: &Face) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        vertex_face_impl(
            s1,
            s2,
            &mut seq1,
            &mut seq2,
            &mut dst,
            self.my_eps,
            self.my_flag,
            self.my_algo,
            &mut modif,
        );
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }

    /// `Perform(E1, E2, Seq1, Seq2)` (`cxx:929-1023`).
    pub fn perform_edge_edge(&mut self, s1: &Edge, s2: &Edge) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        edge_edge_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }

    /// `Perform(E, F, Seq1, Seq2)` (`cxx:1026-1187`).
    pub fn perform_edge_face(&mut self, s1: &Edge, s2: &Face) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        edge_face_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }

    /// `Perform(F1, F2, Seq1, Seq2)` (`cxx:1191-1206`).
    pub fn perform_face_face(&mut self, s1: &Face, s2: &Face) {
        let mut seq1 = std::mem::take(&mut self.my_seq_sol_shape1);
        let mut seq2 = std::mem::take(&mut self.my_seq_sol_shape2);
        let mut dst = self.my_dst_ref;
        let mut modif = self.my_modif;
        face_face_impl(s1, s2, &mut seq1, &mut seq2, &mut dst, self.my_eps, &mut modif);
        self.my_seq_sol_shape1 = seq1;
        self.my_seq_sol_shape2 = seq2;
        self.my_dst_ref = dst;
        self.my_modif = modif;
    }
}

/// `Perform(V1, V2, ...)` (`cxx:780-802`).
fn vertex_vertex_impl(
    the_s1: &Vertex,
    the_s2: &Vertex,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:787-788`.
    let p1 = BRepTool::vertex_point(the_s1);
    let p2 = BRepTool::vertex_point(the_s2);
    // `cxx:789-801`.
    let dst = p1.distance(&p2);
    if dst < *dst_ref - eps || (dst - *dst_ref).abs() < eps {
        if *dst_ref > dst {
            *dst_ref = dst;
        }
        *modif = true;
        seq1.push(SolutionElem::of_vertex(dst, &p1, SupportType::IsVertex, the_s1));
        seq2.push(SolutionElem::of_vertex(dst, &p2, SupportType::IsVertex, the_s2));
    }
}

/// `Perform(V, E, ...)` (`cxx:806-864`).
fn vertex_edge_impl(
    the_s1: &Vertex,
    the_s2: &Edge,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:812-815`.
    if BRepTool::is_degenerated(the_s2) {
        return;
    }
    // `cxx:817-820`.
    let ext = ExtPC::new(the_s1, the_s2);
    let nb_extrema = if ext.is_done() { ext.nb_ext() } else { 0 };
    if nb_extrema == 0 {
        return;
    }
    // `cxx:822-831`.
    let mut dstmin = ext.square_distance(1);
    for i in 2..=nb_extrema {
        let s = ext.square_distance(i);
        if s < dstmin {
            dstmin = s;
        }
    }
    dstmin = dstmin.sqrt();
    // `cxx:832-861`.
    if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
        let p1 = BRepTool::vertex_point(the_s1);
        let eps_p = Precision::PCONFUSION;
        for i in 1..=nb_extrema {
            if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                let pt = ext.point(i);
                // `cxx:844`.
                if tri_solution(seq2, &pt) {
                    let t = ext.parameter(i);
                    // `cxx:847`.
                    if !is_on_boundary(the_s2, &pt, t, eps_p) {
                        if *dst_ref > dstmin {
                            *dst_ref = dstmin;
                        }
                        *modif = true;
                        seq1.push(SolutionElem::of_vertex(
                            dstmin,
                            &p1,
                            SupportType::IsVertex,
                            the_s1,
                        ));
                        seq2.push(SolutionElem::of_edge(
                            dstmin,
                            &pt,
                            SupportType::IsOnEdge,
                            the_s2,
                            t,
                        ));
                    }
                }
            }
        }
    }
}

/// `Perform(V, F, ...)` (`cxx:868-925`).
#[allow(clippy::too_many_arguments)]
fn vertex_face_impl(
    the_s1: &Vertex,
    the_s2: &Face,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    flag: ExtremaExtFlag,
    algo: ExtremaExtAlgo,
    modif: &mut bool,
) {
    // `cxx:875-876`.
    let ext = ExtPF::new(the_s1, the_s2, flag, algo);
    let nb_extrema = if ext.is_done() { ext.nb_ext() } else { 0 };
    if nb_extrema == 0 {
        return;
    }
    // `cxx:878-888`.
    let mut dstmin = ext.square_distance(1);
    for i in 2..=nb_extrema {
        let s = ext.square_distance(i);
        if s < dstmin {
            dstmin = s;
        }
    }
    dstmin = dstmin.sqrt();
    // `cxx:890-922`.
    if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
        let p1 = BRepTool::vertex_point(the_s1);
        let mut classifier = FaceClassifier::new();
        let tol = BRepTool::face_tolerance(the_s2);
        for i in 1..=nb_extrema {
            if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                let pt = ext.point(i).2;
                // `cxx:902`.
                if tri_solution(seq2, &pt) {
                    let (u, v) = ext.parameter(i);
                    classifier.perform(the_s2, GpPnt2d::new(u, v), tol);
                    // `cxx:908`.
                    if classifier.state() == FaceState::In {
                        if *dst_ref > dstmin {
                            *dst_ref = dstmin;
                        }
                        *modif = true;
                        seq1.push(SolutionElem::of_vertex(
                            dstmin,
                            &p1,
                            SupportType::IsVertex,
                            the_s1,
                        ));
                        seq2.push(SolutionElem::of_face(
                            dstmin,
                            &pt,
                            SupportType::IsInFace,
                            the_s2,
                            u,
                            v,
                        ));
                    }
                }
            }
        }
    }
}

/// `Perform(E1, E2, ...)` (`cxx:929-1023`).
fn edge_edge_impl(
    the_s1: &Edge,
    the_s2: &Edge,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:935-938`.
    if BRepTool::is_degenerated(the_s1) || BRepTool::is_degenerated(the_s2) {
        return;
    }
    // `cxx:940`: the reference captured before the extremum search.
    let dst_ref0 = *dst_ref;

    // `cxx:942-943`.
    let ext = ExtCC::new(the_s1, the_s2);
    let nb_extrema = if ext.is_done() {
        if ext.is_parallel() {
            0
        } else {
            ext.nb_ext()
        }
    } else {
        0
    };
    // `cxx:944-989`.
    if nb_extrema > 0 {
        let mut dstmin = ext.square_distance(1);
        for i in 2..=nb_extrema {
            let s = ext.square_distance(i);
            if s < dstmin {
                dstmin = s;
            }
        }
        dstmin = dstmin.sqrt();
        if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
            let eps_p = Precision::PCONFUSION;
            for i in 1..=nb_extrema {
                if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                    let pt1 = ext.point_on_e1(i);
                    let pt2 = ext.point_on_e2(i);
                    // `cxx:969`.
                    if tri_solution(seq1, &pt1) || tri_solution(seq2, &pt2) {
                        let t1 = ext.parameter_on_e1(i);
                        let t2 = ext.parameter_on_e2(i);
                        // `cxx:973`.
                        if !is_on_boundary(the_s1, &pt1, t1, eps_p)
                            && !is_on_boundary(the_s2, &pt2, t2, eps_p)
                        {
                            if *dst_ref > dstmin {
                                *dst_ref = dstmin;
                            }
                            *modif = true;
                            seq1.push(SolutionElem::of_edge(
                                dstmin,
                                &pt1,
                                SupportType::IsOnEdge,
                                the_s1,
                                t1,
                            ));
                            seq2.push(SolutionElem::of_edge(
                                dstmin,
                                &pt2,
                                SupportType::IsOnEdge,
                                the_s2,
                                t2,
                            ));
                        }
                    }
                }
            }
        }
    }

    // `cxx:992-1022`.
    let mut seq_solution1 = Vec::new();
    let mut seq_solution2 = Vec::new();
    perform_c0(
        the_s1,
        the_s2,
        &mut seq_solution1,
        &mut seq_solution2,
        dst_ref0,
        dst_ref,
        eps,
    );
    let (sol1, sol2) = min_solution(&seq_solution1, &seq_solution2, *dst_ref, eps);
    for (s1, s2) in sol1.iter().zip(sol2.iter()) {
        let pt1 = *s1.point();
        let pt2 = *s2.point();
        // `cxx:1014`.
        if tri_solution(seq1, &pt1) || tri_solution(seq2, &pt2) {
            seq1.push(s1.clone());
            seq2.push(s2.clone());
            *modif = true;
        }
    }
}

/// `Perform(E, F, ...)` (`cxx:1026-1187`).
fn edge_face_impl(
    the_s1: &Edge,
    the_s2: &Face,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:1032-1035`.
    if BRepTool::is_degenerated(the_s1) {
        return;
    }

    let mut classifier = FaceClassifier::new();

    // `cxx:1039-1040`.
    let ext = ExtCF::new(the_s1, the_s2);
    let nb_extrema = if ext.is_done() {
        if ext.is_parallel() {
            0
        } else {
            ext.nb_ext()
        }
    } else {
        0
    };
    // `cxx:1041-1095`.
    if nb_extrema > 0 {
        let mut dstmin = ext.square_distance(1);
        for i in 2..=nb_extrema {
            let s = ext.square_distance(i);
            if s < dstmin {
                dstmin = s;
            }
        }
        dstmin = dstmin.sqrt();
        if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
            let tol = BRepTool::face_tolerance(the_s2);
            let eps_p = Precision::PCONFUSION;
            for i in 1..=nb_extrema {
                if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                    let pt1 = ext.point_on_edge(i);
                    let pt2 = ext.point_on_face(i);
                    // `cxx:1069`.
                    if tri_solution(seq1, &pt1) || tri_solution(seq2, &pt2) {
                        let t1 = ext.parameter_on_edge(i);
                        // `cxx:1072`.
                        if !is_on_boundary(the_s1, &pt1, t1, eps_p) {
                            let (u, v) = ext.parameter_on_face(i);
                            classifier.perform(the_s2, GpPnt2d::new(u, v), tol);
                            // `cxx:1079`.
                            if classifier.state() == FaceState::In {
                                if *dst_ref > dstmin {
                                    *dst_ref = dstmin;
                                }
                                *modif = true;
                                seq1.push(SolutionElem::of_edge(
                                    dstmin,
                                    &pt1,
                                    SupportType::IsOnEdge,
                                    the_s1,
                                    t1,
                                ));
                                seq2.push(SolutionElem::of_face(
                                    dstmin,
                                    &pt2,
                                    SupportType::IsInFace,
                                    the_s2,
                                    u,
                                    v,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    // `cxx:1097-1186`: the `GeomAbs_C0` continuation.
    let (a_first, a_last) = BRepTool::edge_parameters(the_s1);
    let p_curv = match BRepTool::edge_curve(the_s1) {
        Some(c) => c,
        None => return,
    };
    // `cxx:1099`.
    if p_curv.continuity() != 0 {
        return;
    }

    // `cxx:1100-1110`.
    let arr_inter = extrema_pc::extrema_intervals(p_curv.as_ref(), 1);
    let tol = BRepTool::face_tolerance(the_s2);

    let mut seq_solution1 = Vec::new();
    let mut seq_solution2 = Vec::new();
    let b = TopoBuilder::new();
    for &a_parameter in &arr_inter {
        // `cxx:1117-1119`.
        let a_pnt = p_curv.d0(a_parameter);
        let v1 = b.make_vertex(a_pnt, Precision::CONFUSION);
        // `cxx:1121`.
        let ext_pf = ExtPF::new(
            &v1,
            the_s2,
            ExtremaExtFlag::MinMax,
            ExtremaExtAlgo::Grad,
        );
        let nb_extrema_pf = if ext_pf.is_done() { ext_pf.nb_ext() } else { 0 };
        if nb_extrema_pf == 0 {
            continue;
        }
        // `cxx:1125-1134`.
        let mut dstmin = ext_pf.square_distance(1);
        for ii in 2..=nb_extrema_pf {
            let s = ext_pf.square_distance(ii);
            if s < dstmin {
                dstmin = s;
            }
        }
        dstmin = dstmin.sqrt();
        // `cxx:1136-1171`.
        if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
            for ii in 1..=nb_extrema_pf {
                if (dstmin - ext_pf.square_distance(ii).sqrt()).abs() < eps {
                    let (u, v) = ext_pf.parameter(ii);
                    classifier.perform(the_s2, GpPnt2d::new(u, v), tol);
                    // `cxx:1146`.
                    if classifier.state() == FaceState::In {
                        if *dst_ref > dstmin {
                            *dst_ref = dstmin;
                        }
                        *modif = true;
                        seq_solution1.push(SolutionElem::of_edge(
                            dstmin,
                            &a_pnt,
                            SupportType::IsOnEdge,
                            the_s1,
                            a_parameter,
                        ));
                        seq_solution2.push(SolutionElem::of_face(
                            dstmin,
                            &ext_pf.point(ii).2,
                            SupportType::IsInFace,
                            the_s2,
                            u,
                            v,
                        ));
                    }
                }
            }
        }
    }

    // `cxx:1173-1186`.
    let (sol1, sol2) = min_solution(&seq_solution1, &seq_solution2, *dst_ref, eps);
    for (s1, s2) in sol1.iter().zip(sol2.iter()) {
        seq1.push(s1.clone());
        seq2.push(s2.clone());
    }
    let _ = (a_first, a_last);
}

/// `Perform(F1, F2, ...)` (`cxx:1191-1206`).
fn face_face_impl(
    the_s1: &Face,
    the_s2: &Face,
    seq1: &mut Vec<SolutionElem>,
    seq2: &mut Vec<SolutionElem>,
    dst_ref: &mut f64,
    eps: f64,
    modif: &mut bool,
) {
    // `cxx:1197-1198`.
    let ext = ExtFF::new(the_s1, the_s2);
    let nb_extrema = if ext.is_done() {
        if ext.is_parallel() {
            0
        } else {
            ext.nb_ext()
        }
    } else {
        0
    };
    if nb_extrema == 0 {
        return;
    }
    // `cxx:1201-1210`.
    let mut dstmin = ext.square_distance(1);
    for i in 2..=nb_extrema {
        let s = ext.square_distance(i);
        if s < dstmin {
            dstmin = s;
        }
    }
    dstmin = dstmin.sqrt();
    // `cxx:1211-1261`.
    if dstmin < *dst_ref - eps || (dstmin - *dst_ref).abs() < eps {
        let tol1 = BRepTool::face_tolerance(the_s1);
        let tol2 = BRepTool::face_tolerance(the_s2);
        let mut classifier = FaceClassifier::new();
        for i in 1..=nb_extrema {
            if (dstmin - ext.square_distance(i).sqrt()).abs() < eps {
                let pt1 = ext.point_on_face1(i);
                let pt2 = ext.point_on_face2(i);
                // `cxx:1229`.
                if tri_solution(seq1, &pt1) || tri_solution(seq2, &pt2) {
                    let (u1, v1) = ext.parameter_on_face1(i);
                    classifier.perform(the_s1, GpPnt2d::new(u1, v1), tol1);
                    // `cxx:1234`.
                    if classifier.state() == FaceState::In {
                        let (u2, v2) = ext.parameter_on_face2(i);
                        classifier.perform(the_s2, GpPnt2d::new(u2, v2), tol2);
                        // `cxx:1240`.
                        if classifier.state() == FaceState::In {
                            if *dst_ref > dstmin {
                                *dst_ref = dstmin;
                            }
                            *modif = true;
                            seq1.push(SolutionElem::of_face(
                                dstmin,
                                &pt1,
                                SupportType::IsInFace,
                                the_s1,
                                u1,
                                v1,
                            ));
                            seq2.push(SolutionElem::of_face(
                                dstmin,
                                &pt2,
                                SupportType::IsInFace,
                                the_s2,
                                u2,
                                v2,
                            ));
                        }
                    }
                }
            }
        }
    }
}
