use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::extrema::point_surface_extrema;
use occt_geom::extrema_pc::point_curve_extrema_all;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topexp::Explorer;

use super::*;

/// The nearest vertex of `shape` to `p`.
pub fn closest_vertex(shape: &TopoShape, p: &GpPnt) -> Option<Vertex> {
    let mut best: Option<Vertex> = None;
    let mut best_d = f64::INFINITY;
    let mut seen = std::collections::HashSet::new();
    let mut ex = Explorer::new(shape, ShapeType::Vertex);
    while ex.more() {
        let vs = ex.current().clone();
        ex.next();
        if !seen.insert(Arc::as_ptr(&vs.tshape) as usize) {
            continue;
        }
        if let Some(v) = Vertex::wrap(vs) {
            let d = BRepTool::vertex_point(&v).distance(p);
            if d < best_d {
                best_d = d;
                best = Some(v);
            }
        }
    }
    best
}

/// Closest point on an edge to `p`, returned as `(parameter, point)`.
///
/// Solves the point–curve extrema with `occt_geom`'s analytic/Newton solver
/// ([`point_curve_extrema_all`]), then clamps the winning parameter to the
/// edge's bounded parameter range (matching the previous sampling behaviour for
/// unbounded curves: a generous `±1e6` window). Every returned extremum is
/// clamped and evaluated so the global minimum of the *bounded* edge is found
/// even when the unbounded curve's minimum lies outside the edge. For an edge
/// without a curve the parameter is `NaN`.
pub fn closest_point_on_edge(edge: &Edge, p: &GpPnt, _samples: usize) -> (f64, GpPnt) {
    let curve = match BRepTool::edge_curve(edge) {
        Some(c) => c,
        None => return (f64::NAN, GpPnt::zero()),
    };
    let (mut f, mut l) = BRepTool::edge_parameters(edge);
    if !f.is_finite() || !l.is_finite() {
        // Unbounded curve (infinite line): clamp to a generous window.
        f = -1e6;
        l = 1e6;
    }
    let (f, l) = if f <= l { (f, l) } else { (l, f) };
    let extrema = point_curve_extrema_all(&*curve, p);
    if extrema.is_empty() {
        // Degenerate curve (e.g. zero-length): golden-section fallback.
        let (u, q) = occt_geom::extrema::refine_curve_point(&*curve, p, f, l);
        return (u, q);
    }
    let mut best_u = f;
    let mut best_d = f64::INFINITY;
    let mut consider = |u: f64| {
        let u = u.clamp(f, l);
        let d = curve.d0(u).square_distance(p);
        if d < best_d {
            best_d = d;
            best_u = u;
        }
    };
    for e in &extrema {
        consider(e.u1);
    }
    // The bounded-edge candidate must also include both endpoints; the solver
    // reports them only when the underlying curve's range is finite.
    consider(f);
    consider(l);
    (best_u, curve.d0(best_u))
}

/// Closest point on a face to `p`, returned as `((u, v), point)`.
///
/// Solves the point–surface extrema with `occt_geom`'s analytic/Newton solver
/// ([`point_surface_extrema`]), clamps the parameter to the face's UV window,
/// and then also checks each boundary edge (a bounded face's closest point can
/// lie on an edge; the previous UV-grid sampler found those implicitly). When
/// an edge wins, its `(u, v)` is recovered by projecting the 3D point back onto
/// the surface.
pub fn closest_point_on_face(face: &Face, p: &GpPnt, _nu: usize, _nv: usize) -> ((f64, f64), GpPnt) {
    let surface = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return ((0.0, 0.0), GpPnt::zero()),
    };
    let (u1, u2, v1, v2) = match face_uv_window(face) {
        Some(w) => w,
        None => return ((0.0, 0.0), GpPnt::zero()),
    };
    // `Extrema_ExtPS` on the surface; `None` is `!IsDone()` (no extremum at
    // all), which leaves only the boundary-edge candidates below.
    let e = point_surface_extrema(&*surface, p);
    // The analytic solver reports its own parameter frame (a plane is
    // reconstructed with arbitrary X/Y axes), so take its frame-independent
    // closest POINT and re-project it into the face's own (u, v) before
    // clamping to the face's UV window.
    let mut best_u = u1;
    let mut best_v = v1;
    let mut best_p = surface.d0(best_u, best_v);
    let mut best_d = best_p.square_distance(p);
    if let Some(e) = e {
        let (pu, pv) = project_to_surface(&surface, &e.p2);
        let pu = pu.clamp(u1, u2);
        let pv = pv.clamp(v1, v2);
        let q = surface.d0(pu, pv);
        let d = q.square_distance(p);
        if d < best_d {
            best_d = d;
            best_u = pu;
            best_v = pv;
            best_p = q;
        }
    }
    for edge in face_boundary_edges(face) {
        let (_, q) = closest_point_on_edge(&edge, p, 0);
        let d = q.square_distance(p);
        if d < best_d {
            let (pu, pv) = project_to_surface(&surface, &q);
            best_d = d;
            best_u = pu;
            best_v = pv;
            best_p = q;
        }
    }
    ((best_u, best_v), best_p)
}

