//! Remaining `IntTools_Tools` helpers used by `MakeCurve` / `PrepareLines3D`.
//!
//! Source: `IntTools_Tools.cxx`:
//! - `HasInternalEdge` (59)
//! - `IsClosed` (78)
//! - `RejectLines` (106)
//! - `IsDirsCoinside` (164, 177)
//! - `SplitCurve` (191)
//! - `IntermediatePoint` (254)
//! - `IsVertex` (263, 284, 309, 326)
//! - `ComputeVV` (354)
//! - `MakeFaceFromWireAndFace` (376)
//! - `ClassifyPointByFace` (390)
//! - `IsMiddlePointsEqual` (404)
//! - `CurveTolerance` (430)
//!
//! `CheckCurve` / `IsOnPave` / `SegPln` / `ComputeIntRange` live in
//! `int_tools_curve_box` and `int_tools_segpln`.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, SQUARE_CONFUSION};
use occt_geom::{Curve, GeomTrimmedCurve};
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::intermediate_point;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::{FClass2d, FaceState};
use crate::int_face_face::FaceFaceCurve;
use crate::inttools_data::{CurveKind, IntRange};
use crate::shape::{Edge, Face, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::edges_of_wire;
use crate::topexp::Explorer;

pub use crate::boptools_2d::intermediate_point as int_tools_intermediate_point;

/// `IntTools_Tools::HasInternalEdge`.
pub fn has_internal_edge(wire: &Wire) -> bool {
    for e in edges_of_wire(wire) {
        if e.orientation() == Orientation::Internal {
            return true;
        }
    }
    false
}

/// `IntTools_Tools::IsClosed` on a 3D curve.
pub fn is_closed_curve(a_c3d: &dyn Curve) -> bool {
    let a_f = a_c3d.first_parameter();
    let a_l = a_c3d.last_parameter();
    if !a_f.is_finite() || !a_l.is_finite() {
        return false;
    }
    let a_p1 = a_c3d.d0(a_f);
    let a_p2 = a_c3d.d0(a_l);
    a_p1.square_distance(&a_p2) < SQUARE_CONFUSION
}

/// `IntTools_Tools::IsDirsCoinside(D1, D2)` with default `dLim = 0.0002`.
pub fn is_dirs_coinside(d1: &GpDir, d2: &GpDir) -> bool {
    is_dirs_coinside_lim(d1, d2, 0.0002)
}

/// `IntTools_Tools::IsDirsCoinside(D1, D2, dLim)`.
pub fn is_dirs_coinside_lim(d1: &GpDir, d2: &GpDir, d_lim: f64) -> bool {
    let p1 = GpPnt::new(d1.x(), d1.y(), d1.z());
    let p2 = GpPnt::new(d2.x(), d2.y(), d2.z());
    let d = p1.distance(&p2);
    d < d_lim || (2.0 - d).abs() < d_lim
}

fn curve_direction_at_mid(c: &dyn Curve) -> Option<GpDir> {
    let a = c.first_parameter();
    let b = c.last_parameter();
    let t = if a.is_finite() && b.is_finite() {
        0.5 * (a + b)
    } else {
        0.0
    };
    let (_, tau) = c.d1(t);
    GpDir::from_vec(&tau).ok()
}

/// `IntTools_Tools::RejectLines`.
///
/// Line-ness is taken from [`FaceFaceCurve::kind`] (no `dyn Curve` downcast).
pub fn reject_lines(a_s_in: &[FaceFaceCurve]) -> Vec<FaceFaceCurve> {
    let a_nb = a_s_in.len();
    let mut a_s_out = Vec::new();
    let mut a_d1: Option<GpDir> = None;
    for (i, ic) in a_s_in.iter().enumerate() {
        if ic.kind != CurveKind::Line {
            return a_s_in.to_vec();
        }
        let Some(a_d2) = curve_direction_at_mid(ic.curve.as_ref()) else {
            return a_s_in.to_vec();
        };
        if i == 0 {
            a_s_out.push(ic.clone());
            a_d1 = Some(a_d2);
            continue;
        }
        let Some(d1) = a_d1 else {
            return a_s_in.to_vec();
        };
        if !is_dirs_coinside(&d1, &a_d2) {
            a_s_out.push(ic.clone());
            return a_s_out;
        }
        let _ = a_nb;
    }
    a_s_out
}

/// `IntTools_Tools::SplitCurve` — split a closed curve at the intermediate
/// parameter into two trimmed halves.
pub fn split_curve(ic: &FaceFaceCurve) -> Vec<FaceFaceCurve> {
    if !is_closed_curve(ic.curve.as_ref()) {
        return Vec::new();
    }
    let a_f = ic.curve.first_parameter();
    let a_l = ic.curve.last_parameter();
    let mut a_mid = 0.5 * (a_f + a_l);
    if ic.kind == CurveKind::BSpline || ic.kind == CurveKind::Other {
        a_mid = intermediate_point(a_f, a_l);
    }
    let c3d_f: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(ic.curve.clone(), a_f, a_mid));
    let c3d_l: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(ic.curve.clone(), a_mid, a_l));
    let pc1_f = ic.pcurve1.as_ref().map(|pc| {
        Arc::new(Geom2dTrimmedCurve::new(pc.clone(), a_f, a_mid)) as Arc<dyn occt_geom2d::curve::Curve2d>
    });
    let pc1_l = ic.pcurve1.as_ref().map(|pc| {
        Arc::new(Geom2dTrimmedCurve::new(pc.clone(), a_mid, a_l)) as Arc<dyn occt_geom2d::curve::Curve2d>
    });
    let pc2_f = ic.pcurve2.as_ref().map(|pc| {
        Arc::new(Geom2dTrimmedCurve::new(pc.clone(), a_f, a_mid)) as Arc<dyn occt_geom2d::curve::Curve2d>
    });
    let pc2_l = ic.pcurve2.as_ref().map(|pc| {
        Arc::new(Geom2dTrimmedCurve::new(pc.clone(), a_mid, a_l)) as Arc<dyn occt_geom2d::curve::Curve2d>
    });
    let mk = |curve: Arc<dyn Curve>, p1, p2| {
        let range = IntRange::new_unchecked(curve.first_parameter(), curve.last_parameter());
        FaceFaceCurve {
            kind: ic.kind,
            curve,
            range,
            face1_idx: ic.face1_idx,
            face2_idx: ic.face2_idx,
            pcurve1: p1,
            pcurve2: p2,
            tolerance: ic.tolerance,
            tangential_tolerance: ic.tangential_tolerance,
        }
    };
    vec![mk(c3d_f, pc1_f, pc2_f), mk(c3d_l, pc1_l, pc2_l)]
}

/// `IntTools_Tools::IsVertex(P, TolPV, V)`.
pub fn is_vertex_pnt(a_p: &GpPnt, a_tol_pv: f64, a_v: &Vertex) -> bool {
    let a_tol_v = BRepTool::vertex_tolerance(a_v) + a_tol_pv + CONFUSION;
    let a_pv = BRepTool::vertex_point(a_v);
    a_pv.square_distance(a_p) <= a_tol_v * a_tol_v
}

/// `IntTools_Tools::IsVertex(E, V, t)`.
pub fn is_vertex_edge_at(a_e: &Edge, a_v: &Vertex, t: f64) -> bool {
    let Some(c) = GeometryRegistry::global().edge_curve(&a_e.0) else {
        return false;
    };
    let a_pt = c.d0(t);
    let a_tol_v = BRepTool::vertex_tolerance(a_v);
    let a_pv = BRepTool::vertex_point(a_v);
    a_pv.square_distance(&a_pt) < a_tol_v * a_tol_v
}

/// `IntTools_Tools::IsVertex(E, t)` — OCCT overwrites the vertex-tolerance
/// square with `1.e-12` before the distance test.
pub fn is_vertex_on_edge(a_e: &Edge, t: f64) -> bool {
    let Some(c) = GeometryRegistry::global().edge_curve(&a_e.0) else {
        return false;
    };
    let a_pt = c.d0(t);
    let mut exp = Explorer::new(&a_e.0, ShapeType::Vertex);
    while exp.more() {
        let a_v = Vertex(exp.current().clone());
        let a_pv = BRepTool::vertex_point(&a_v);
        if a_pv.square_distance(&a_pt) < 1.0e-12 {
            return true;
        }
        exp.next();
    }
    false
}

/// `IntTools_Tools::ComputeVV`. Returns `0` when the vertices coincide.
pub fn compute_vv(a_v1: &Vertex, a_v2: &Vertex) -> i32 {
    let a_tol_sum = BRepTool::vertex_tolerance(a_v1) + BRepTool::vertex_tolerance(a_v2);
    let a_p1 = BRepTool::vertex_point(a_v1);
    let a_p2 = BRepTool::vertex_point(a_v2);
    if a_p1.square_distance(&a_p2) < a_tol_sum * a_tol_sum {
        0
    } else {
        -1
    }
}

/// `IntTools_Tools::MakeFaceFromWireAndFace`.
pub fn make_face_from_wire_and_face(a_w: &Wire, a_f: &Face) -> Face {
    match BRepTool::face_surface(a_f) {
        Some(s) => TopoBuilder::new().make_face(s, std::slice::from_ref(a_w)),
        None => a_f.clone(),
    }
}

/// `IntTools_Tools::ClassifyPointByFace`.
pub fn classify_point_by_face(a_f: &Face, a_p2d: GpPnt2d) -> FaceState {
    let a_face_tolerance = BRepTool::face_tolerance(a_f);
    match FClass2d::new(a_f, a_face_tolerance) {
        Ok(c) => c.perform(a_p2d),
        Err(_) => FaceState::Unknown,
    }
}

/// `IntTools_Tools::IsMiddlePointsEqual`.
pub fn is_middle_points_equal(a_e1: &Edge, a_e2: &Edge) -> bool {
    let Some(c1) = GeometryRegistry::global().edge_curve(&a_e1.0) else {
        return false;
    };
    let Some(c2) = GeometryRegistry::global().edge_curve(&a_e2.0) else {
        return false;
    };
    let (f1, l1) = GeometryRegistry::global().edge_parameters(&a_e1.0);
    let (f2, l2) = GeometryRegistry::global().edge_parameters(&a_e2.0);
    let m1 = 0.5 * (f1 + l1);
    let m2 = 0.5 * (f2 + l2);
    let a_p1 = c1.d0(m1);
    let a_p2 = c2.d0(m2);
    let a_sum_tol = BRepTool::edge_tolerance(a_e1) + BRepTool::edge_tolerance(a_e2);
    a_p1.square_distance(&a_p2) < a_sum_tol * a_sum_tol
}

/// `IntTools_Tools::CurveTolerance`. Without a trimmed-parabola downcast the
/// base tolerance is returned (OCCT also returns `aTolBase` for non-parabolas).
pub fn curve_tolerance(_a_c3d: &dyn Curve, a_tol_base: f64) -> f64 {
    a_tol_base
}
