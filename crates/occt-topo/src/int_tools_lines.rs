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
pub fn is_closed_curve(c3d: &dyn Curve) -> bool {
    let f = c3d.first_parameter();
    let l = c3d.last_parameter();
    if !f.is_finite() || !l.is_finite() {
        return false;
    }
    let p1 = c3d.d0(f);
    let p2 = c3d.d0(l);
    p1.square_distance(&p2) < SQUARE_CONFUSION
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
pub fn reject_lines(s_in: &[FaceFaceCurve]) -> Vec<FaceFaceCurve> {
    let nb = s_in.len();
    let mut s_out = Vec::new();
    let mut a_d1: Option<GpDir> = None;
    for (i, ic) in s_in.iter().enumerate() {
        if ic.kind != CurveKind::Line {
            return s_in.to_vec();
        }
        let Some(d2) = curve_direction_at_mid(ic.curve.as_ref()) else {
            return s_in.to_vec();
        };
        if i == 0 {
            s_out.push(ic.clone());
            a_d1 = Some(d2);
            continue;
        }
        let Some(d1) = a_d1 else {
            return s_in.to_vec();
        };
        if !is_dirs_coinside(&d1, &d2) {
            s_out.push(ic.clone());
            return s_out;
        }
        let _ = nb;
    }
    s_out
}

/// `IntTools_Tools::SplitCurve` — split a closed curve at the intermediate
/// parameter into two trimmed halves.
pub fn split_curve(ic: &FaceFaceCurve) -> Vec<FaceFaceCurve> {
    if !is_closed_curve(ic.curve.as_ref()) {
        return Vec::new();
    }
    let f = ic.curve.first_parameter();
    let l = ic.curve.last_parameter();
    let mut mid = 0.5 * (f + l);
    if ic.kind == CurveKind::BSpline || ic.kind == CurveKind::Other {
        mid = intermediate_point(f, l);
    }
    let c3d_f: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(ic.curve.clone(), f, mid));
    let c3d_l: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(ic.curve.clone(), mid, l));
    let pc1_f = ic.pcurve1.as_ref().map(|pc| {
        Arc::new(Geom2dTrimmedCurve::new(pc.clone(), f, mid)) as Arc<dyn occt_geom2d::curve::Curve2d>
    });
    let pc1_l = ic.pcurve1.as_ref().map(|pc| {
        Arc::new(Geom2dTrimmedCurve::new(pc.clone(), mid, l)) as Arc<dyn occt_geom2d::curve::Curve2d>
    });
    let pc2_f = ic.pcurve2.as_ref().map(|pc| {
        Arc::new(Geom2dTrimmedCurve::new(pc.clone(), f, mid)) as Arc<dyn occt_geom2d::curve::Curve2d>
    });
    let pc2_l = ic.pcurve2.as_ref().map(|pc| {
        Arc::new(Geom2dTrimmedCurve::new(pc.clone(), mid, l)) as Arc<dyn occt_geom2d::curve::Curve2d>
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
pub fn is_vertex_pnt(p: &GpPnt, tol_pv: f64, v: &Vertex) -> bool {
    let tol_v = BRepTool::vertex_tolerance(v) + tol_pv + CONFUSION;
    let pv = BRepTool::vertex_point(v);
    pv.square_distance(p) <= tol_v * tol_v
}

/// `IntTools_Tools::IsVertex(E, V, t)`.
pub fn is_vertex_edge_at(a_e: &Edge, v: &Vertex, t: f64) -> bool {
    let Some(c) = GeometryRegistry::global().edge_curve(&a_e.0) else {
        return false;
    };
    let pt = c.d0(t);
    let tol_v = BRepTool::vertex_tolerance(v);
    let pv = BRepTool::vertex_point(v);
    pv.square_distance(&pt) < tol_v * tol_v
}

/// `IntTools_Tools::IsVertex(E, t)` — OCCT overwrites the vertex-tolerance
/// square with `1.e-12` before the distance test.
pub fn is_vertex_on_edge(a_e: &Edge, t: f64) -> bool {
    let Some(c) = GeometryRegistry::global().edge_curve(&a_e.0) else {
        return false;
    };
    let pt = c.d0(t);
    let mut exp = Explorer::new(&a_e.0, ShapeType::Vertex);
    while exp.more() {
        let v = Vertex(exp.current().clone());
        let pv = BRepTool::vertex_point(&v);
        if pv.square_distance(&pt) < 1.0e-12 {
            return true;
        }
        exp.next();
    }
    false
}

/// `IntTools_Tools::ComputeVV`. Returns `0` when the vertices coincide.
pub fn compute_vv(v1: &Vertex, v2: &Vertex) -> i32 {
    let tol_sum = BRepTool::vertex_tolerance(v1) + BRepTool::vertex_tolerance(v2);
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    if p1.square_distance(&p2) < tol_sum * tol_sum {
        0
    } else {
        -1
    }
}

/// `IntTools_Tools::MakeFaceFromWireAndFace`.
pub fn make_face_from_wire_and_face(w: &Wire, f: &Face) -> Face {
    match BRepTool::face_surface(f) {
        Some(s) => TopoBuilder::new().make_face(s, std::slice::from_ref(w)),
        None => f.clone(),
    }
}

/// `IntTools_Tools::ClassifyPointByFace`.
pub fn classify_point_by_face(f: &Face, p2d: GpPnt2d) -> FaceState {
    let a_face_tolerance = BRepTool::face_tolerance(f);
    match FClass2d::new(f, a_face_tolerance) {
        Ok(c) => c.perform(p2d),
        Err(_) => FaceState::Unknown,
    }
}

/// `IntTools_Tools::IsMiddlePointsEqual`.
pub fn is_middle_points_equal(e1: &Edge, e2: &Edge) -> bool {
    let Some(c1) = GeometryRegistry::global().edge_curve(&e1.0) else {
        return false;
    };
    let Some(c2) = GeometryRegistry::global().edge_curve(&e2.0) else {
        return false;
    };
    let (f1, l1) = GeometryRegistry::global().edge_parameters(&e1.0);
    let (f2, l2) = GeometryRegistry::global().edge_parameters(&e2.0);
    let m1 = 0.5 * (f1 + l1);
    let m2 = 0.5 * (f2 + l2);
    let p1 = c1.d0(m1);
    let p2 = c2.d0(m2);
    let sum_tol = BRepTool::edge_tolerance(e1) + BRepTool::edge_tolerance(e2);
    p1.square_distance(&p2) < sum_tol * sum_tol
}

/// `IntTools_Tools::CurveTolerance`. Without a trimmed-parabola downcast the
/// base tolerance is returned (OCCT also returns `aTolBase` for non-parabolas).
pub fn curve_tolerance(_a_c3d: &dyn Curve, tol_base: f64) -> f64 {
    tol_base
}
