//! Unit tests of the node-insertion algorithms.

    use super::*;
    use crate::meshing::model_builder::ModelBuilder;
    use crate::primitives::BRepPrimBox;
    use crate::tgeometry::GeometryRegistry;

    /// The 3x3 lattice of UV points.
    fn grid_3x3() -> Vec<GpPnt2d> {
        (0..3)
            .flat_map(|i| (0..3).map(move |j| GpPnt2d::new(i as f64, j as f64)))
            .collect()
    }

    /// Octagonal boundary (counter-clockwise, positive area) of the 3x3 grid's
    /// convex hull plus its closing loop.
    fn octagon() -> (Vec<GpPnt2d>, Vec<(i32, i32)>) {
        let pts = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(2.0, 0.0),
            GpPnt2d::new(2.0, 1.0),
            GpPnt2d::new(2.0, 2.0),
            GpPnt2d::new(1.0, 2.0),
            GpPnt2d::new(0.0, 2.0),
            GpPnt2d::new(0.0, 1.0),
        ];
        let edges = (0..8).map(|i| (i as i32, ((i + 1) % 8) as i32)).collect();
        (pts, edges)
    }

    #[test]
    fn grid_3x3_uv_points_inserted_into_delaun_produce_triangles() {
        let pts = grid_3x3();
        let mut algo = DelaunayNodeInsertionMeshAlgo::new();
        assert_eq!(algo.insert_internal_nodes(&pts).unwrap(), 9);

        let res = algo.triangulate().expect("triangulate");
        assert!(!res.triangles.is_empty());
        // OCCT BRepMesh_Delaun.cxx:703 calls ProcessConstraints() unconditionally;
        // frontierAdjust() ends with cleanupMesh() (cxx:1028) which prunes boundary
        // triangles whose neighbor touches the super-triangle. These nodes carry only
        // Free links, so the mesh keeps 6 triangles over 7 nodes instead of 8 over 9.
        assert_eq!(res.triangles.len(), 6);
        assert_eq!(res.nodes.len(), 7);
        // Every triangle references live node indices only.
        for t in &res.triangles {
            for &v in &t.vertex_indices {
                assert!((1..=9).contains(&v), "triangle references node {v}");
            }
        }
    }

    #[test]
    fn boundary_nodes_registered_as_frontier_with_constraint_links() {
        let pts = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(0.0, 1.0),
        ];
        let edges = vec![(0i32, 1), (1, 2), (2, 3), (3, 0)];

        let mut algo = DelaunayNodeInsertionMeshAlgo::new();
        assert_eq!(algo.insert_boundary_nodes(&pts, &vec![GpPnt::zero(); 4], &edges).unwrap(), 4);
        assert_eq!(algo.boundary_indices(), &[1, 2, 3, 4][..]);

        let ds = algo.structure();
        assert_eq!(ds.nb_nodes(), 4);
        for i in 1..=4 {
            assert_eq!(ds.node(i).state, VertexState::Frontier);
        }
        // The four constraint edges are registered as Frontier links.
        assert_eq!(ds.links_of_domain().len(), 4);
        for &l in ds.links_of_domain() {
            assert_eq!(ds.link_movability(l), VertexState::Frontier);
        }
    }

    #[test]
    fn internal_insertion_yields_2n_minus_2_minus_h() {
        let (boundary, edges) = octagon();
        let mut algo = DelaunayNodeInsertionMeshAlgo::new();
        assert_eq!(algo.insert_boundary_nodes(&boundary, &vec![GpPnt::zero(); 8], &edges).unwrap(), 8);
        assert_eq!(algo.insert_internal_nodes(&[GpPnt2d::new(1.0, 1.0)]).unwrap(), 1);

        let res = algo.triangulate().expect("triangulate");
        let n = 9; // 8 boundary + 1 internal point
        let h = 8; // octagonal hull
        assert_eq!(res.triangles.len() as i32, 2 * n - 2 - h);
        assert_eq!(res.triangles.len(), 8);
    }

    #[test]
    fn perform_errors_on_face_without_boundary_pcurves() {
        let shape = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0.clone();
        let params = MeshParameters::default();
        let mut model = ModelBuilder::build_model(&shape, &params).expect("model built");
        let mut algo = DelaunayNodeInsertionMeshAlgo::new();
        // Pcurves are empty until the edge discretizers run, so no boundary UV
        // points can be collected and the face must be rejected.
        let err = algo.perform(&mut model, 0, &params);
        assert!(err.is_err(), "unpopulated pcurves yield no boundary UV points");
        GeometryRegistry::global().clear_shape(&shape);
    }
