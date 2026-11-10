use super::*;

use occt_core::gp::GpPnt;

use crate::brep_tool::BRepTool;
use occt_core::bnd::BndBox;

fn box_solid() -> test_box::UnitBox {
    test_box::unit_box()
}

#[test]
fn point_to_box_distance() {
    let b = box_solid();
    let d = point_shape_distance(&GpPnt::new(2.0, 0.0, 0.0), &b.solid.0);
    assert!((d - 1.0).abs() < 1e-6, "distance {d}");
}

#[test]
fn closest_vertex_is_corner() {
    let b = box_solid();
    let v = closest_vertex(&b.solid.0, &GpPnt::new(2.0, 0.0, 0.0)).expect("a vertex");
    let p = BRepTool::vertex_point(&v);
    assert!(p.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9, "vertex at {p:?}");
}

#[test]
fn point_in_box() {
    let b = box_solid();
    assert!(is_inside(&b.solid.0, &GpPnt::new(0.5, 0.5, 0.5)));
    assert!(!is_inside(&b.solid.0, &GpPnt::new(2.0, 0.0, 0.0)));
    assert!(!is_inside(&b.solid.0, &GpPnt::new(-0.5, 0.5, 0.5)));
}

#[test]
fn closest_on_edge_and_face() {
    let b = box_solid();
    let edge = &b.edges[0]; // (0,0,0)-(1,0,0)
    let (u, q) = closest_point_on_edge(edge, &GpPnt::new(0.5, 1.0, 0.0), 32);
    assert!((u - 0.5).abs() < 1e-3, "u {u}");
    assert!(q.distance(&GpPnt::new(0.5, 0.0, 0.0)) < 1e-6, "q {q:?}");

    let face = &b.faces[2]; // front (y = 0)
    let ((u, v), q) = closest_point_on_face(face, &GpPnt::new(0.5, 1.0, 0.5), 8, 8);
    assert!((u - 0.5).abs() < 1e-2, "u {u}");
    assert!((v - 0.5).abs() < 1e-2, "v {v}");
    assert!(q.distance(&GpPnt::new(0.5, 0.0, 0.5)) < 1e-6, "q {q:?}");
}

#[test]
fn shape_self_distance_zero() {
    let b = box_solid();
    let d = shape_distance(&b.solid.0, &b.solid.0);
    assert!(d < 1e-6, "distance {d}");
}

#[test]
fn bbox_separation() {
    let a = BndBox::from_corners(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
    let b = BndBox::from_corners(&GpPnt::new(3.0, 0.0, 0.0), &GpPnt::new(4.0, 1.0, 1.0));
    assert!((shape_bbox_distance(&a, &b) - 2.0).abs() < 1e-12);
    let c = BndBox::from_corners(&GpPnt::new(0.5, 0.5, 0.5), &GpPnt::new(2.0, 2.0, 2.0));
    assert_eq!(shape_bbox_distance(&a, &c), 0.0);
}
