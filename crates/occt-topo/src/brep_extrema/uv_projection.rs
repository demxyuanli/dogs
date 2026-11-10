use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::Surface;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};


/// A finite parametric window `(u_min, u_max, v_min, v_max)` for a face.
///
/// Uses the surface's own range when finite; otherwise projects the boundary
/// vertices of the face's wires onto the surface (Newton solve) and takes the
/// bounding box of the resulting `(u, v)` points.
pub(super) fn face_uv_window(face: &Face) -> Option<(f64, f64, f64, f64)> {
    let surface = BRepTool::face_surface(face)?;
    let (u1, u2, v1, v2) = BRepTool::uv_bounds(face);
    if u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite() {
        return Some((u1, u2, v1, v2));
    }
    let pts = face_boundary_points(face);
    if pts.is_empty() {
        return None;
    }
    let mut us: Vec<f64> = Vec::new();
    let mut vs: Vec<f64> = Vec::new();
    for p in &pts {
        let (u, v) = project_to_surface(&surface, p);
        us.push(u);
        vs.push(v);
    }
    let ua = us.iter().copied().fold(f64::INFINITY, f64::min);
    let ub = us.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let va = vs.iter().copied().fold(f64::INFINITY, f64::min);
    let vb = vs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let (ua, ub) = if ub - ua < 1e-12 { (ua - 1.0, ub + 1.0) } else { (ua, ub) };
    let (va, vb) = if vb - va < 1e-12 { (va - 1.0, vb + 1.0) } else { (va, vb) };
    Some((ua, ub, va, vb))
}

/// Project a 3D point onto a surface by Newton iteration on the `d1` partials.
///
/// Converges in a single step for planes (affine mapping); a few steps for
/// other analytic surfaces.
pub(super) fn project_to_surface(surface: &Arc<dyn Surface>, p: &GpPnt) -> (f64, f64) {
    let mut u = 0.0;
    let mut v = 0.0;
    for _ in 0..20 {
        let (pt, du, dv) = surface.d1(u, v);
        let rx = p.x() - pt.x();
        let ry = p.y() - pt.y();
        let rz = p.z() - pt.z();
        let duu = du.coord.square_modulus();
        let dvv = dv.coord.square_modulus();
        let duv = du.coord.dot(&dv.coord);
        let b1 = rx * du.x() + ry * du.y() + rz * du.z();
        let b2 = rx * dv.x() + ry * dv.y() + rz * dv.z();
        let det = duu * dvv - duv * duv;
        if det.abs() < 1e-30 {
            break;
        }
        let d_u = (b1 * dvv - b2 * duv) / det;
        let d_v = (b2 * duu - b1 * duv) / det;
        u += d_u;
        v += d_v;
        if (d_u * d_u + d_v * d_v).sqrt() < 1e-12 {
            break;
        }
    }
    (u, v)
}

/// All edges of a face's wires (boundary edges).
pub(super) fn face_boundary_edges(face: &Face) -> Vec<Edge> {
    let mut out: Vec<Edge> = Vec::new();
    let face_kids = face.0.tshape.read().unwrap().children.clone();
    for wire in face_kids {
        if wire.shape_type() != ShapeType::Wire {
            continue;
        }
        let wire_kids = wire.tshape.read().unwrap().children.clone();
        for e in wire_kids {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            if let Some(edge) = Edge::wrap(e) {
                out.push(edge);
            }
        }
    }
    out
}

/// Distinct 3D endpoint points of all edges in a face's wires.
pub(super) fn face_boundary_points(face: &Face) -> Vec<GpPnt> {
    let mut out: Vec<GpPnt> = Vec::new();
    let face_kids = face.0.tshape.read().unwrap().children.clone();
    for wire in face_kids {
        if wire.shape_type() != ShapeType::Wire {
            continue;
        }
        let wire_kids = wire.tshape.read().unwrap().children.clone();
        for e in wire_kids {
            if e.shape_type() != ShapeType::Edge {
                continue;
            }
            if let Some(edge) = Edge::wrap(e) {
                if let Some((a, b)) = BRepTool::edge_vertices(&edge) {
                    if !out.iter().any(|q| q.distance(&a) < 1e-9) {
                        out.push(a);
                    }
                    if !out.iter().any(|q| q.distance(&b) < 1e-9) {
                        out.push(b);
                    }
                }
            }
        }
    }
    out
}
