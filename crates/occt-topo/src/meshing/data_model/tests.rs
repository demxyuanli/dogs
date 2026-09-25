use super::prelude::*;
use super::*;
    use crate::builder::TopoBuilder;
    use crate::tgeometry::GeometryRegistry;
    use occt_core::gp::{GpAx3, GpPln, GpPnt, GpPnt2d};

    /// Release registry entries for a shape tree so tests don't leave stale
    /// geometry keyed by a freed Arc address in the process-wide side-table.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn sample_edge(b: &TopoBuilder) -> Edge {
        b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0))
    }

    fn sample_wire(b: &TopoBuilder) -> Wire {
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        b.make_wire(&[e1, e2])
    }

    fn sample_face(b: &TopoBuilder) -> Face {
        b.make_face_plane(&GpPln::new(GpAx3::standard()))
    }

    #[test]
    fn mesh_status_bit_flags() {
        let mut s = MeshStatus::NO_ERROR;
        assert!(s.is_empty());
        assert_eq!(s.status_mask(), 0x0);
        s.set_status(MeshStatus::OPEN_WIRE);
        assert!(s.is_set(MeshStatus::OPEN_WIRE));
        assert!(!s.is_set(MeshStatus::FAILURE));
        s.set_status(MeshStatus::FAILURE);
        assert_eq!(s.status_mask(), 0x1 | 0x4);
        assert_eq!(s.union(MeshStatus::REMESH).status_mask(), 0x1 | 0x4 | 0x8);
        s.unset_status(MeshStatus::OPEN_WIRE);
        assert!(!s.is_set(MeshStatus::OPEN_WIRE));
        assert_eq!(s.status_mask(), 0x4);
    }

    #[test]
    fn model_add_remove_query() {
        let b = TopoBuilder::new();
        let mut model = MeshModel::new(TopoShape::new(crate::abs::ShapeType::Compound));
        assert_eq!(model.faces_nb(), 0);
        assert_eq!(model.edges_nb(), 0);
        assert_eq!(model.wires_nb(), 0);

        let e1 = sample_edge(&b);
        let e2 = sample_edge(&b);
        let ei0 = model.add_edge(e1);
        let ei1 = model.add_edge(e2);
        assert_eq!(model.edges_nb(), 2);

        let fi0 = model.add_face(sample_face(&b));
        assert_eq!(model.faces_nb(), 1);

        let wi0 = model.add_wire(sample_wire(&b));
        assert_eq!(model.wires_nb(), 1);

        // Wire references both edges; face references the wire.
        model.wire_mut(wi0).unwrap().add_edge(ei0, Orientation::Forward);
        model.wire_mut(wi0).unwrap().add_edge(ei1, Orientation::Reversed);
        model.face_mut(fi0).unwrap().add_wire(wi0);
        assert_eq!(model.wire(wi0).unwrap().edges_nb(), 2);
        assert_eq!(model.face(fi0).unwrap().wires_nb(), 1);
        assert_eq!(model.wire(wi0).unwrap().edge(0).unwrap(), ei0);
        assert_eq!(
            model.wire(wi0).unwrap().edge_orientation(1).unwrap(),
            Orientation::Reversed
        );

        // Query with out-of-range index fails.
        assert!(model.edge(99).is_err());
        assert!(model.face(99).is_err());

        // Removal returns the removed item and shrinks the collection.
        let removed = model.remove_edge(ei0).unwrap();
        assert_eq!(removed.pcurves_nb(), 0);
        assert_eq!(model.edges_nb(), 1);

        let removed_face = model.remove_face(fi0).unwrap();
        assert_eq!(removed_face.wires_nb(), 1);
        assert_eq!(model.faces_nb(), 0);

        clear_tree(&model.shape().unwrap().clone());
    }

    #[test]
    fn model_status_aggregation() {
        let b = TopoBuilder::new();
        let mut model = MeshModel::new(TopoShape::new(crate::abs::ShapeType::Compound));
        assert_eq!(model.status_mask(), 0);

        let fi = model.add_face(sample_face(&b));
        let ei = model.add_edge(sample_edge(&b));
        model.face_mut(fi).unwrap().set_status(MeshStatus::FAILURE);
        model.edge_mut(ei).unwrap().set_status(MeshStatus::OPEN_WIRE);
        model.edge_mut(ei).unwrap().set_status(MeshStatus::UNORIENTED_WIRE);

        assert_eq!(model.status_mask(), 0x4 | 0x1 | 0x10);
        assert!(model.has_status(MeshStatus::FAILURE));
        assert!(model.has_status(MeshStatus::OPEN_WIRE));
        assert!(!model.has_status(MeshStatus::REUSED));

        // Unsetting a flag on one entity is reflected in the aggregation.
        model.edge_mut(ei).unwrap().unset_status(MeshStatus::OPEN_WIRE);
        assert!(!model.has_status(MeshStatus::OPEN_WIRE));

        // A face that is only ReMesh / UnorientedWire stays valid.
        model.face_mut(fi).unwrap().unset_status(MeshStatus::FAILURE);
        model.face_mut(fi).unwrap().set_status(MeshStatus::REMESH);
        assert!(model.face(fi).unwrap().is_valid());

        clear_tree(&model.shape().unwrap().clone());
    }

    #[test]
    fn curve_and_pcurve_creation() {
        let mut c = MeshCurve::new();
        assert_eq!(c.parameters_nb(), 0);
        c.add_point(GpPnt::new(0.0, 0.0, 0.0), 0.0);
        c.add_point(GpPnt::new(1.0, 0.0, 0.0), 1.0);
        c.add_point(GpPnt::new(2.0, 0.0, 0.0), 2.0);
        c.insert_point(1, GpPnt::new(0.5, 0.0, 0.0), 0.5).unwrap();
        assert_eq!(c.parameters_nb(), 4);
        assert_eq!(c.get_parameter(1).unwrap(), 0.5);
        assert!(c.get_point(1).unwrap().is_equal(&GpPnt::new(0.5, 0.0, 0.0)));
        assert_eq!(c.points().len(), c.parameters().len());

        // Out-of-range access is an error.
        assert!(c.get_point(99).is_err());

        // Clear keeps only the end points when requested.
        c.clear(true);
        assert_eq!(c.parameters_nb(), 2);
        assert_eq!(c.get_parameter(0).unwrap(), 0.0);
        assert_eq!(c.get_parameter(1).unwrap(), 2.0);
        c.clear(false);
        assert_eq!(c.parameters_nb(), 0);

        // Deflection attribute.
        c.set_deflection(0.01);
        assert_eq!(c.deflection(), 0.01);

        let mut pc = MeshPCurve::new(3, Orientation::Reversed);
        pc.add_point(GpPnt2d::new(0.0, 0.0), 0.0);
        pc.add_point(GpPnt2d::new(1.0, 1.0), 1.0);
        assert_eq!(pc.parameters_nb(), 2);
        assert_eq!(pc.get_point(1).unwrap(), GpPnt2d::new(1.0, 1.0));
        assert_eq!(pc.get_index(0).unwrap(), 0);
        *pc.get_index_mut(1).unwrap() = 7;
        assert_eq!(pc.get_index(1).unwrap(), 7);
        assert!(!pc.is_forward());
        assert!(!pc.is_internal());
        assert_eq!(pc.face(), 3);
        assert_eq!(pc.orientation(), Orientation::Reversed);
        assert_eq!(pc.indices().len(), pc.points().len());

        pc.set_deflection(0.02);
        assert_eq!(pc.deflection(), 0.02);

        let fwd = MeshPCurve::new(0, Orientation::Forward);
        assert!(fwd.is_forward());
    }

    #[test]
    fn mesh_edge_reads_topological_flags_and_curve() {
        let b = TopoBuilder::new();
        let e = sample_edge(&b);
        let mut me = MeshEdge::new(e);
        assert!(me.curve().is_some());
        assert_eq!(me.first_parameter(), 0.0);
        assert_eq!(me.last_parameter(), 1.0);
        assert!(me.same_param());
        assert!(me.same_range());
        assert!(!me.degenerated());
        assert!(me.is_free());
        assert_eq!(me.angular_deflection(), f64::MAX);

        // Add pcurves and look them up per face / orientation.
        let fi0 = me.add_pcurve(0, Orientation::Forward);
        me.add_pcurve(0, Orientation::Reversed);
        me.add_pcurve(1, Orientation::Forward);
        assert_eq!(me.pcurves_nb(), 3);
        assert_eq!(me.pcurve(fi0).unwrap().orientation(), Orientation::Forward);
        assert_eq!(me.pcurve_for(0, Orientation::Reversed).unwrap().orientation(), Orientation::Reversed);
        assert_eq!(me.pcurve_for(1, Orientation::Forward).unwrap().face(), 1);
        assert_eq!(me.pcurves_for(0), vec![0, 1]);
        assert!(!me.is_free());

        clear_tree(&me.edge().0);
    }
