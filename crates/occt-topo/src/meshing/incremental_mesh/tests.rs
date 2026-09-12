use super::prelude::*;
use super::*;
    use crate::mesh::mesh_surface_area;
    use crate::primitives::BRepPrimBox;

    fn unit_box() -> TopoShape {
        BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0.clone()
    }

    #[test]
    fn perform_produces_triangles_for_box() {
        let shape = unit_box();
        let mut inc = IncrementalMesh::new();
        inc.set_shape(&shape);
        inc.change_parameters().deflection = 0.05;
        let mesh = inc.perform().expect("mesh produced");
        assert!(inc.is_done(), "IsDone set after perform");
        assert!(inc.is_modified(), "modified after perform");
        assert!(mesh.triangles.len() > 0, "triangles {}", mesh.triangles.len());
        assert!(mesh.vertices.len() >= 8, "vertices {}", mesh.vertices.len());
        // Unit box: surface area ≈ 6.
        let area = mesh_surface_area(&mesh);
        assert!((area - 6.0).abs() < 0.5, "box area {area}");
    }

    #[test]
    fn status_flags_are_no_error_for_clean_box() {
        let shape = unit_box();
        let params = MeshParameters { deflection: 0.05, ..MeshParameters::default() };
        let inc = IncrementalMesh::from_parameters(&shape, params);
        assert!(inc.mesh().is_some(), "mesh cached");
        assert_eq!(inc.get_status_flags(), 0, "flags {}", inc.get_status_flags());
    }

    #[test]
    fn constructor_performs_automatically() {
        let shape = unit_box();
        let inc = IncrementalMesh::from_deflection(&shape, 0.05, false, 0.5);
        assert!(inc.is_done(), "constructor performs");
        assert!(inc.mesh().is_some(), "constructor caches mesh");
        let m = inc.mesh().unwrap();
        assert!(m.triangles.len() > 0);
    }

    #[test]
    fn perform_without_shape_fails() {
        let mut inc = IncrementalMesh::new();
        assert!(inc.perform().is_err());
        assert!(!inc.is_done());
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let shape = unit_box();
        let mut inc = IncrementalMesh::new();
        inc.set_shape(&shape);
        inc.change_parameters().deflection = 0.0; // below Precision::Confusion()
        assert!(inc.perform().is_err());
        assert!(!inc.is_done());
    }

    #[test]
    fn incremental_mesh_to_shape_mesh_produces_triangles() {
        let shape = unit_box();
        let mesh = incremental_mesh_to_shape_mesh(&shape, 0.05).expect("mesh produced");
        assert!(mesh.triangles.len() > 0, "triangles {}", mesh.triangles.len());
        assert!(!mesh.vertices.is_empty(), "vertices non-empty");
        // The 6 box faces tile the unit surface: area ≈ 6.
        let area = mesh_surface_area(&mesh);
        assert!((area - 6.0).abs() < 0.5, "box area {area}");
    }

    #[test]
    fn pipeline_preserves_face_count() {
        let shape = unit_box();
        let mut inc = IncrementalMesh::new();
        inc.set_shape(&shape);
        inc.change_parameters().deflection = 0.05;
        // Build + pre-process the model and run the per-face pipeline directly.
        let mut model = ModelBuilder::build_model(&shape, inc.parameters()).expect("model built");
        ModelPreProcessor::perform(&mut model, inc.parameters());
        let tris = inc.triangulate_model_faces(&mut model).expect("pipeline ran");
        assert_eq!(tris.len(), 6, "box has 6 faces");
        for ft in &tris {
            assert!(ft.triangles.len() > 0, "face {} has triangles", ft.face_index);
        }
    }
