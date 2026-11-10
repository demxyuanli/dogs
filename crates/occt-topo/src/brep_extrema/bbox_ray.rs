
use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpXyz};

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{TopoShape, Vertex};
use crate::topexp::Explorer;


/// Axis-aligned bounding box of a shape's vertices.
pub(super) fn shape_bbox(shape: &TopoShape) -> Option<BndBox> {
    let mut b = BndBox::new();
    let mut any = false;
    let mut ex = Explorer::new(shape, ShapeType::Vertex);
    while ex.more() {
        let vs = ex.current().clone();
        ex.next();
        if let Some(v) = Vertex::wrap(vs) {
            b.add_point(&BRepTool::vertex_point(&v));
            any = true;
        }
    }
    if any {
        Some(b)
    } else {
        None
    }
}

/// Möller–Trumbore ray/triangle test for a ray from `origin` in the +X
/// direction.
/// Möller–Trumbore ray (+X) / triangle intersection; returns the hit parameter
/// `t` when the ray crosses the triangle.
pub(super) fn ray_plus_x_hits(a: &GpPnt, b: &GpPnt, c: &GpPnt, origin: &GpPnt) -> Option<f64> {
    const EPS: f64 = 1e-9;
    let dir = GpXyz::new(1.0, 0.0, 0.0);
    let edge1 = b.coord.subtracted(&a.coord);
    let edge2 = c.coord.subtracted(&a.coord);
    let h = dir.crossed(&edge2);
    let det = edge1.dot(&h);
    if det.abs() < EPS {
        return None;
    }
    let inv = 1.0 / det;
    let s = origin.coord.subtracted(&a.coord);
    let u = s.dot(&h) * inv;
    if u < -EPS || u > 1.0 + EPS {
        return None;
    }
    let q = s.crossed(&edge1);
    let v = dir.dot(&q) * inv;
    if v < -EPS || u + v > 1.0 + EPS {
        return None;
    }
    let t = edge2.dot(&q) * inv;
    if t > EPS { Some(t) } else { None }
}
