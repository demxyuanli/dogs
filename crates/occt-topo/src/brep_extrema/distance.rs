
use occt_core::gp::GpPnt;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topexp::Explorer;

use super::*;

/// Minimum distance from a point to any part of a shape.
///
/// Exact for vertices; exact analytic / Newton for edges and faces (via
/// [`closest_point_on_edge`] / [`closest_point_on_face`]).
pub fn point_shape_distance(p: &GpPnt, shape: &TopoShape) -> f64 {
    let mut best = f64::INFINITY;

    let mut ex = Explorer::new(shape, ShapeType::Vertex);
    while ex.more() {
        let vs = ex.current().clone();
        ex.next();
        if let Some(v) = Vertex::wrap(vs) {
            best = best.min(BRepTool::vertex_point(&v).distance(p));
        }
    }

    let mut ex = Explorer::new(shape, ShapeType::Edge);
    while ex.more() {
        let es = ex.current().clone();
        ex.next();
        if let Some(e) = Edge::wrap(es) {
            let (_, q) = closest_point_on_edge(&e, p, 32);
            best = best.min(q.distance(p));
        }
    }

    let mut ex = Explorer::new(shape, ShapeType::Face);
    while ex.more() {
        let fs = ex.current().clone();
        ex.next();
        if let Some(f) = Face::wrap(fs) {
            let (_, q) = closest_point_on_face(&f, p, 16, 16);
            best = best.min(q.distance(p));
        }
    }

    best
}

/// Approximate minimum distance between two shapes.
///
/// The distance is minimized over vertex–vertex, vertex–edge, vertex–face,
/// edge–vertex and edge–face sample pairs. Edge–edge pairs are covered
/// indirectly through the sampled edge points vs. faces. The result is a
/// coarse upper bound on the true distance; it is not an exact extremum
/// computation.
pub fn shape_distance(a: &TopoShape, b: &TopoShape) -> f64 {
    let va = collect_vertices(a);
    let vb = collect_vertices(b);
    let ea = collect_edges(a);
    let eb = collect_edges(b);
    let _fa = collect_faces(a);
    let fb = collect_faces(b);

    let mut best = f64::INFINITY;

    for av in &va {
        for bv in &vb {
            best = best.min(BRepTool::vertex_point(av).distance(&BRepTool::vertex_point(bv)));
        }
    }
    for av in &va {
        let p = BRepTool::vertex_point(av);
        for be in &eb {
            let (_, q) = closest_point_on_edge(be, &p, 16);
            best = best.min(q.distance(&p));
        }
    }
    for av in &va {
        let p = BRepTool::vertex_point(av);
        for bf in &fb {
            let (_, q) = closest_point_on_face(bf, &p, 8, 8);
            best = best.min(q.distance(&p));
        }
    }
    for bv in &vb {
        let p = BRepTool::vertex_point(bv);
        for ae in &ea {
            let (_, q) = closest_point_on_edge(ae, &p, 16);
            best = best.min(q.distance(&p));
        }
    }
    for ae in &ea {
        let (f, l) = BRepTool::edge_parameters(ae);
        if !f.is_finite() || !l.is_finite() {
            continue;
        }
        let curve = BRepTool::edge_curve(ae);
        if curve.is_none() {
            continue;
        }
        let curve = curve.unwrap();
        for i in 0..=8 {
            let u = f + (l - f) * (i as f64 / 8.0);
            let q = curve.d0(u);
            for bf in &fb {
                let (_, fq) = closest_point_on_face(bf, &q, 8, 8);
                best = best.min(fq.distance(&q));
            }
        }
    }

    best
}

