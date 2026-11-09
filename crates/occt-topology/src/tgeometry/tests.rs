use super::*;

use crate::abs::ShapeType;
use occt_geom::GeomLine;

#[test]
fn vertex_roundtrip() {
    let reg = GeometryRegistry::global();
    let v = TopoShape::new(ShapeType::Vertex);
    reg.set_vertex(&v, VertexGeom { point: GpPnt::new(1.0, 2.0, 3.0), tolerance: 1e-3 });
    let g = reg.vertex_geom(&v).expect("registered");
    assert!(g.point.is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
    assert_eq!(g.tolerance, 1e-3);
    assert!(reg.vertex_point(&v).is_equal(&GpPnt::new(1.0, 2.0, 3.0)));
    reg.clear_shape(&v);
    assert!(reg.vertex_geom(&v).is_none());
}

#[test]
fn edge_roundtrip_with_curve() {
    use occt_core::gp::{GpAx1, GpDir, GpLin};
    let reg = GeometryRegistry::global();
    let e = TopoShape::new(ShapeType::Edge);
    let lin = GpLin::new(GpAx1::new(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap()));
    let curve = Arc::new(GeomLine::new(lin));
    reg.set_edge(&e, EdgeGeom::new(curve, 0.0, 5.0));
    assert_eq!(reg.edge_parameters(&e), (0.0, 5.0));
    let c = reg.edge_curve(&e).expect("curve present");
    let p = c.d0(2.5);
    assert!((p.x() - 2.5).abs() < 1e-12);
    reg.clear_shape(&e);
    assert_eq!(reg.edge_parameters(&e), (f64::NEG_INFINITY, f64::INFINITY));
}

#[test]
fn unregistered_shapes_fall_back() {
    let reg = GeometryRegistry::global();
    let e = TopoShape::new(ShapeType::Edge);
    // Default: unbounded range, no curve, origin point for vertex.
    assert_eq!(reg.edge_parameters(&e), (f64::NEG_INFINITY, f64::INFINITY));
    assert!(reg.edge_curve(&e).is_none());
    let v = TopoShape::new(ShapeType::Vertex);
    assert!(reg.vertex_point(&v).is_equal(&GpPnt::zero()));
}
