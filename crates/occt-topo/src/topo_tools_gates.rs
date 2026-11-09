//! Unit tests for `occt_topology::topo_tools`, kept in `occt-topo` because the
//! reference solid needs `crate::primitives::BRepPrimBox`.

#[cfg(test)]
mod topo_tools_gates {
    use occt_topology::builder::TopoBuilder;
    use occt_topology::topo_tools::*;
    use crate::primitives::BRepPrimBox;
    use occt_core::gp::GpPnt;
    use occt_topology::abs::ShapeType;

    #[test]
    fn box_counts_and_types() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let c = shape_counts(&b.solid.0);
        assert_eq!(c[&ShapeType::Vertex], 8);
        assert_eq!(c[&ShapeType::Edge], 12);
        assert_eq!(c[&ShapeType::Face], 6);
        assert_eq!(c[&ShapeType::Wire], 6);
        assert_eq!(c[&ShapeType::Shell], 1);
        assert_eq!(c[&ShapeType::Solid], 1);
    }

    #[test]
    fn map_shapes_and_specific_collectors() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let map = map_shapes(&b.solid.0, &[ShapeType::Vertex, ShapeType::Face]);
        assert_eq!(map[&ShapeType::Vertex].len(), 8);
        assert_eq!(map[&ShapeType::Face].len(), 6);
        assert!(vertices_of(&b.solid.0).len() == 8);
        assert!(edges_of(&b.solid.0).len() == 12);
        assert!(faces_of(&b.solid.0).len() == 6);
    }

    #[test]
    fn wire_closed_and_edge_vertices() {
        let b = TopoBuilder::new();
        // A closed square wire.
        let p = [GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.),
                 GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.)];
        let e1 = b.make_edge_segment(&p[0], &p[1]);
        let e2 = b.make_edge_segment(&p[1], &p[2]);
        let e3 = b.make_edge_segment(&p[2], &p[3]);
        let e4 = b.make_edge_segment(&p[3], &p[0]);
        let wire = b.make_wire(&[e1.clone(), e2, e3, e4]);
        assert!(wire_is_closed(&wire));

        let (a, z) = edge_vertices(&e1);
        assert!(vertex_position(&a.unwrap()).is_equal(&p[0]));
        assert!(vertex_position(&z.unwrap()).is_equal(&p[1]));
    }

    #[test]
    fn structure_valid_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(structure_is_valid(&b.solid.0));
    }

    #[test]
    fn is_same_and_wires_of_face() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&b.solid.0);
        assert!(!faces.is_empty());
        let w = wires_of_face(&faces[0]);
        assert_eq!(w.len(), 1);
        let ws = edges_of_wire(&w[0]);
        assert_eq!(ws.len(), 4);
        assert!(is_same(&faces[0].0, &faces[0].0));
    }
}
