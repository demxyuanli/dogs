//! `TopOpeBRepDS_InterferenceTool` — the interference factory functions.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_InterferenceTool.hxx:30-90`,
//! `TopOpeBRepDS_InterferenceTool.cxx:1-144`.
//!
//! OCCT returns `handle<TopOpeBRepDS_Interference>`; the port's handle is
//! [`InterferenceHandle`] = `Arc<RwLock<Interference>>`, so every factory wraps
//! the value it constructs in a fresh handle, exactly as `new` does.
//!
//! `Parameter(I)` / `Parameter(I, Par)` (`cxx:131-143`) are
//! `down_cast<CurvePointInterference>`, so the port takes the `Parameter()`
//! that [`Interference::parameter`] already exposes.

use std::sync::Arc;
use std::sync::RwLock;

use super::interference::{Interference, InterferenceHandle};
use super::kind::{Config, Kind};
use super::transition::Transition;

/// Wraps a constructed interference in a fresh `handle`.
pub fn make_handle(i: Interference) -> InterferenceHandle {
    Arc::new(RwLock::new(i))
}

/// `MakeEdgeInterference(T, SK, SI, GK, GI, P)` (`cxx:29-38`): a
/// `TopOpeBRepDS_CurvePointInterference`.
pub fn make_edge_interference(
    t: &Transition,
    sk: Kind,
    si: i32,
    gk: Kind,
    gi: i32,
    p: f64,
) -> InterferenceHandle {
    make_handle(Interference::curve_point(t.clone(), sk, si, gk, gi, p))
}

/// `MakeCurveInterference(T, SK, SI, GK, GI, P)` (`cxx:42-52`): identical to
/// [`make_edge_interference`] in OCCT.
pub fn make_curve_interference(
    t: &Transition,
    sk: Kind,
    si: i32,
    gk: Kind,
    gi: i32,
    p: f64,
) -> InterferenceHandle {
    make_handle(Interference::curve_point(t.clone(), sk, si, gk, gi, p))
}

/// `DuplicateCurvePointInterference(I)` (`cxx:56-64`): the transition is
/// complemented, the four indices and the parameter are copied.
pub fn duplicate_curve_point_interference(i: &Interference) -> InterferenceHandle {
    make_handle(Interference::curve_point(
        i.transition().complement(),
        i.support_type(),
        i.support(),
        i.geometry_type(),
        i.geometry(),
        // cxx:62 calls Parameter(I), a down_cast to CurvePointInterference.
        i.parameter().unwrap_or(0.0),
    ))
}

/// `MakeFaceCurveInterference(T, FaceI, CurveI, PC)` (`cxx:68-80`).
pub fn make_face_curve_interference(
    transition: &Transition,
    face_i: i32,
    curve_i: i32,
    pc: Option<Arc<dyn occt_geom2d::curve::Curve2d>>,
) -> InterferenceHandle {
    make_handle(Interference::surface_curve(
        transition.clone(),
        Kind::Face,
        face_i,
        Kind::Curve,
        curve_i,
        pc,
    ))
}

/// `MakeSolidSurfaceInterference(T, SolidI, SurfaceI)` (`cxx:84-94`).
pub fn make_solid_surface_interference(
    transition: &Transition,
    solid_i: i32,
    surface_i: i32,
) -> InterferenceHandle {
    make_handle(Interference::solid_surface(
        transition.clone(),
        Kind::Solid,
        solid_i,
        Kind::Surface,
        surface_i,
    ))
}

/// `MakeEdgeVertexInterference(T, EdgeI, VertexI, VertexIsBound, C, param)`
/// (`cxx:98-114`).
pub fn make_edge_vertex_interference(
    transition: &Transition,
    edge_i: i32,
    vertex_i: i32,
    vertex_is_bound: bool,
    c: Config,
    param: f64,
) -> InterferenceHandle {
    make_handle(Interference::edge_vertex_on_edge(
        transition.clone(),
        edge_i,
        vertex_i,
        vertex_is_bound,
        c,
        param,
    ))
}

/// `MakeFaceEdgeInterference(T, FaceI, EdgeI, EdgeIsBound, C)` (`cxx:118-127`).
pub fn make_face_edge_interference(
    transition: &Transition,
    face_i: i32,
    edge_i: i32,
    edge_is_bound: bool,
    c: Config,
) -> InterferenceHandle {
    make_handle(Interference::face_edge(
        transition.clone(),
        face_i,
        edge_i,
        edge_is_bound,
        c,
    ))
}

/// `InterferenceTool::Parameter(I)` (`cxx:131-134`).
pub fn parameter(i: &Interference) -> f64 {
    i.parameter().unwrap_or(0.0)
}

/// `InterferenceTool::Parameter(I, Par)` (`cxx:138-143`).
pub fn set_parameter(i: &mut Interference, par: f64) {
    i.set_parameter(par);
}
