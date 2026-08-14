//! Real bounding box computed from geometry (not from a pre-built mesh).
//!
//! Dispatch on shape type: vertices contribute their point, edges sample their
//! curve, faces sample their surface on a UV grid (and recurse their wires),
//! and containers recurse. Missing geometry is skipped gracefully (void box).

use occt_core::bnd::BndBox;

use crate::abs::ShapeType;
use crate::shape::{Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::wireframe::face_uv_bounds;

fn reg() -> &'static GeometryRegistry {
    GeometryRegistry::global()
}

/// Axis-aligned bounding box of `shape` derived from its registered geometry.
pub fn shape_bbox(shape: &TopoShape) -> BndBox {
    let mut b = BndBox::new();
    add_to_bbox(shape, &mut b);
    b
}

/// Bounding box enlarged by `margin` in every direction.
pub fn shape_bbox_with_margin(shape: &TopoShape, margin: f64) -> BndBox {
    let mut b = shape_bbox(shape);
    b.enlarge(margin);
    b
}

/// Quick reject: `true` if `shape`'s bounding box does not overlap `other`.
pub fn is_shape_outside(shape: &TopoShape, other: &BndBox) -> bool {
    shape_bbox(shape).is_out_box(other)
}

fn add_to_bbox(s: &TopoShape, b: &mut BndBox) {
    match s.shape_type() {
        ShapeType::Vertex => {
            b.add_point(&reg().vertex_point(s));
        }
        ShapeType::Edge => {
            if let Some(curve) = reg().edge_curve(s) {
                let (a0, a1) = reg().edge_parameters(s);
                if a0.is_finite() && a1.is_finite() && a1 > a0 {
                    for i in 0..=64 {
                        let u = a0 + (a1 - a0) * i as f64 / 64.0;
                        b.add_point(&curve.d0(u));
                    }
                }
            }
            add_children(s, b);
        }
        ShapeType::Face => {
            if let Some(surf) = reg().face_surface(s) {
                let f = Face(s.clone());
                let (u0, u1, v0, v1) = face_uv_bounds(&f, surf.as_ref());
                if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()
                    && u1 > u0 && v1 > v0
                {
                    for i in 0..=16 {
                        for j in 0..=16 {
                            let u = u0 + (u1 - u0) * i as f64 / 16.0;
                            let v = v0 + (v1 - v0) * j as f64 / 16.0;
                            b.add_point(&surf.d0(u, v));
                        }
                    }
                }
            }
            add_children(s, b);
        }
        // Wire, Shell, Solid, Compound, CompSolid: recurse children.
        _ => add_children(s, b),
    }
}

fn add_children(s: &TopoShape, b: &mut BndBox) {
    let kids = s.tshape.read().unwrap().children.clone();
    for child in kids {
        add_to_bbox(&child, b);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpPnt;
    use crate::shape::Vertex;
    use crate::tgeometry::VertexGeom;
    use crate::wireframe::tests::{segment_edge, square_face};

    #[test]
    fn vertex_bbox_contains_point() {
        let v = Vertex::new();
        reg().set_vertex(&v.0, VertexGeom { point: GpPnt::new(1.0, 2.0, 3.0), tolerance: 0.0 });
        let b = shape_bbox(&v.0);
        assert!(!b.is_out(&GpPnt::new(1.0, 2.0, 3.0)));
        assert!(b.is_out(&GpPnt::new(5.0, 2.0, 3.0)));
    }

    #[test]
    fn edge_bbox_from_curve() {
        let e = segment_edge(0.0, 0.0, 1.0, 0.0);
        let b = shape_bbox(&e.0);
        let (x0, x1, y0, y1, z0, z1) = b.get().unwrap();
        assert!(x0.abs() < 1e-9 && (x1 - 1.0).abs() < 1e-9);
        assert!(y0.abs() < 1e-9 && y1.abs() < 1e-9);
        assert!(z0.abs() < 1e-9 && z1.abs() < 1e-9);
    }

    #[test]
    fn face_bbox_from_surface_and_wires() {
        let f = square_face();
        let b = shape_bbox(&f.0);
        let (x0, x1, y0, y1, z0, z1) = b.get().unwrap();
        assert!(x0.abs() < 1e-9 && (x1 - 1.0).abs() < 1e-9);
        assert!(y0.abs() < 1e-9 && (y1 - 1.0).abs() < 1e-9);
        assert!(z0.abs() < 1e-9 && z1.abs() < 1e-9);
    }

    #[test]
    fn bbox_with_margin_enlarges() {
        let v = Vertex::new();
        reg().set_vertex(&v.0, VertexGeom { point: GpPnt::zero(), tolerance: 0.0 });
        let b = shape_bbox_with_margin(&v.0, 0.5);
        let (x0, x1, y0, _y1, z0, _z1) = b.get().unwrap();
        assert!((x0 + 0.5).abs() < 1e-12 && (x1 - 0.5).abs() < 1e-12);
        assert!((y0 + 0.5).abs() < 1e-12 && (z0 + 0.5).abs() < 1e-12);
    }

    #[test]
    fn is_shape_outside_quick_reject() {
        let v = Vertex::new();
        reg().set_vertex(&v.0, VertexGeom { point: GpPnt::zero(), tolerance: 0.0 });
        let far = BndBox::from_corners(&GpPnt::new(10.0, 10.0, 10.0), &GpPnt::new(11.0, 11.0, 11.0));
        assert!(is_shape_outside(&v.0, &far));
        let near = BndBox::from_corners(&GpPnt::new(-1.0, -1.0, -1.0), &GpPnt::new(1.0, 1.0, 1.0));
        assert!(!is_shape_outside(&v.0, &near));
    }
}
