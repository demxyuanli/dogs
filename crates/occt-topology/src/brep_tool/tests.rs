use super::*;
use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use occt_core::gp::{GpAx3, GpPln, GpPnt};


fn edge() -> Edge {
    Edge(TopoShape::new(ShapeType::Edge))
}

/// Release registry entries for a shape tree so tests don't leave stale
/// geometry keyed by a freed Arc address in the process-wide side-table.
fn clear_tree(s: &TopoShape) {
    GeometryRegistry::global().clear_shape(s);
    let children = s.tshape.read().unwrap().children.clone();
    for c in children {
        clear_tree(&c);
    }
}

#[test]
fn edge_parameters_are_unbounded() {
    let (first, last) = BRepTool::edge_parameters(&edge());
    assert_eq!(first, f64::NEG_INFINITY);
    assert_eq!(last, f64::INFINITY);
}

#[test]
fn unbounded_edge_is_not_closed() {
    assert!(!BRepTool::is_closed_edge(&edge()));
}

#[test]
fn shape_type_str_names_vertex() {
    assert_eq!(BRepTool::shape_type_str(ShapeType::Vertex), "Vertex");
}

#[test]
fn empty_shape_list_is_integral() {
    let shapes: Vec<TopoShape> = Vec::new();
    assert!(BRepTool::check_shape_integrity(&shapes));
}

#[test]
fn make_vertex_stores_point() {
    let b = TopoBuilder::new();
    let v = b.make_vertex(GpPnt::new(1.0, 2.0, 3.0), 1e-3);
    assert!(BRepTool::vertex_point(&v).is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
    assert_eq!(BRepTool::vertex_tolerance(&v), 1e-3);
    clear_tree(&v.0);
}

#[test]
fn make_edge_stores_curve() {
    let b = TopoBuilder::new();
    let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(10.0, 0.0, 0.0));
    assert_eq!(BRepTool::edge_parameters(&e), (0.0, 10.0));
    let c = BRepTool::edge_curve(&e).expect("curve present");
    let mid = c.d0(5.0);
    assert!(mid.is_equal(&GpPnt::new(5.0, 0.0, 0.0)));
    assert!(BRepTool::same_parameter(&e));
    assert!(!BRepTool::is_degenerated(&e));
    assert_eq!(BRepTool::parameter_on_edge(&e, 0), 0.0);
    assert_eq!(BRepTool::parameter_on_edge(&e, 1), 10.0);
    clear_tree(&e.0);
}

#[test]
fn segment_has_two_vertices_and_endpoints() {
    let b = TopoBuilder::new();
    let p1 = GpPnt::new(0.0, 0.0, 0.0);
    let p2 = GpPnt::new(3.0, 4.0, 0.0);
    let e = b.make_edge_segment(&p1, &p2);
    assert_eq!(e.0.tshape.read().unwrap().children.len(), 2);
    let (a, z) = BRepTool::edge_vertices(&e).expect("endpoints");
    assert!(a.is_equal(&p1));
    assert!(z.is_equal(&p2));
    clear_tree(&e.0);
}

#[test]
fn make_face_plane_registers_surface() {
    let b = TopoBuilder::new();
    let f = b.make_face_plane(&GpPln::new(GpAx3::standard()));
    assert!(BRepTool::face_surface(&f).is_some());
    let (u1, u2, v1, v2) = BRepTool::uv_bounds(&f);
    assert_eq!(u1, f64::NEG_INFINITY);
    assert_eq!(u2, f64::INFINITY);
    assert_eq!(v1, f64::NEG_INFINITY);
    assert_eq!(v2, f64::INFINITY);
    assert!(BRepTool::natural_restriction(&f));
    clear_tree(&f.0);
}

#[test]
fn child_trees_build() {
    let b = TopoBuilder::new();
    let e1 = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(1.0, 0.0, 0.0));
    let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
    let wire = b.make_wire(&[e1, e2]);
    let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
    let shell = b.make_shell(&[face.clone()]);
    let solid = b.make_solid(&[shell.clone()]);
    assert_eq!(wire.0.tshape.read().unwrap().children.len(), 2);
    assert_eq!(face.0.tshape.read().unwrap().children.len(), 0);
    assert_eq!(shell.0.tshape.read().unwrap().children.len(), 1);
    assert_eq!(solid.0.tshape.read().unwrap().children.len(), 1);
    clear_tree(&solid.0);
    clear_tree(&wire.0);
}

#[test]
fn unregistered_shapes_fall_back() {
    let v = Vertex::new();
    assert!(BRepTool::vertex_point(&v).is_equal(&GpPnt::zero()));
    assert_eq!(BRepTool::vertex_tolerance(&v), occt_core::precision::CONFUSION);
    assert!(BRepTool::edge_curve(&edge()).is_none());
    assert!(BRepTool::edge_vertices(&edge()).is_none());
    assert!(BRepTool::face_surface(&Face::new()).is_none());
    assert_eq!(
        BRepTool::uv_bounds(&Face::new()),
        (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY)
    );
    assert!(BRepTool::natural_restriction(&Face::new()));
}
