use super::prelude::*;
use super::*;
    use std::sync::Arc;

    use occt_core::gp::{GpDir, GpLin};
    use occt_geom::GeomLine;

    use crate::brep_extrema::test_box::unit_box;
    use crate::builder::TopoBuilder;
    use crate::shape::TopoShape;
    use crate::tgeometry::GeometryRegistry;

    /// Drop the global geometry entries owned by a shape subtree.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("unit direction")
    }

    // -----------------------------------------------------------------------
    // Classifier cache
    // -----------------------------------------------------------------------

    #[test]
    fn classifier_cache_reused_per_face() {
        let bx = unit_box();
        let mut ctx = IntToolsContext::new();
        assert_eq!(ctx.fclass2d_cache_len(), 0);

        // First query on the bottom face builds the classifier.
        assert!(ctx
            .is_point_in_face(&bx.faces[0], &GpPnt::new(0.25, 0.75, 0.0), None, 1e-6)
            .unwrap());
        assert_eq!(ctx.fclass2d_cache_len(), 1);

        // A second query on the same face reuses it.
        assert!(ctx
            .is_point_in_face(&bx.faces[0], &GpPnt::new(0.5, 0.5, 0.0), None, 1e-6)
            .unwrap());
        assert_eq!(ctx.fclass2d_cache_len(), 1);

        // A different face builds a second classifier.
        let _ = ctx
            .is_point_in_face(&bx.faces[1], &GpPnt::new(0.25, 0.25, 1.0), None, 1e-6)
            .unwrap();
        assert_eq!(ctx.fclass2d_cache_len(), 2);

        ctx.clear_cached();
        assert_eq!(ctx.fclass2d_cache_len(), 0);
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Point-in-face
    // -----------------------------------------------------------------------

    #[test]
    fn is_point_in_face_inside_outside_boundary() {
        let bx = unit_box();
        let face = &bx.faces[0]; // bottom (z = 0)
        let mut ctx = IntToolsContext::new();
        let tol = 1e-6;

        // Face centre -> strictly In.
        assert!(ctx.is_point_in_face(face, &GpPnt::new(0.5, 0.5, 0.0), None, tol).unwrap());
        // In the plane but outside the boundary -> false.
        assert!(!ctx.is_point_in_face(face, &GpPnt::new(-0.5, 0.5, 0.0), None, tol).unwrap());
        // Above the surface (distance > tol) -> false.
        assert!(!ctx.is_point_in_face(face, &GpPnt::new(0.5, 0.5, 1.0), None, tol).unwrap());

        // A boundary point: strict In is false, in-or-on is true, state is On.
        let bnd = GpPnt::new(0.5, 0.0, 0.0);
        assert!(!ctx.is_point_in_face(face, &bnd, None, tol).unwrap());
        assert!(ctx.is_point_in_on_face(face, &bnd, None, tol).unwrap());
        let (u, v) = ctx.project_point_on_face(face, &bnd).unwrap();
        assert_eq!(ctx.state_point_face(face, (u, v), tol).unwrap(), FaceState::On);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn is_point_in_face_explicit_uv_skips_projection() {
        let bx = unit_box();
        let face = &bx.faces[0];
        let mut ctx = IntToolsContext::new();
        // Supplying the exact UV (0.5, 0.5) classifies In without projecting.
        let p = GpPnt::new(0.5, 0.5, 0.0);
        assert!(ctx.is_point_in_face(face, &p, Some((0.5, 0.5)), 1e-9).unwrap());
        // A supplied UV that is far from the 3D point fails the distance check.
        assert!(!ctx.is_point_in_face(face, &p, Some((0.5, -0.5)), 1e-9).unwrap());
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn state_point_face_in_on_out() {
        let bx = unit_box();
        let face = &bx.faces[0];
        let mut ctx = IntToolsContext::new();
        let tol = 1e-6;
        assert_eq!(ctx.state_point_face(face, (0.5, 0.5), tol).unwrap(), FaceState::In);
        assert_eq!(ctx.state_point_face(face, (0.5, 0.0), tol).unwrap(), FaceState::On);
        assert_eq!(ctx.state_point_face(face, (-0.5, 0.5), tol).unwrap(), FaceState::Out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn is_valid_point_for_face_in_on_out() {
        let bx = unit_box();
        let face = &bx.faces[0];
        let mut ctx = IntToolsContext::new();
        assert!(ctx.is_valid_point_for_face((0.5, 0.5), face).unwrap());
        // Boundary counts as valid (in-or-on).
        assert!(ctx.is_valid_point_for_face((0.5, 0.0), face).unwrap());
        assert!(!ctx.is_valid_point_for_face((-0.5, 0.5), face).unwrap());
        assert!(ctx
            .is_valid_point_for_faces((0.5, 0.5), &bx.faces[0], &bx.faces[1])
            .unwrap());
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Blocks
    // -----------------------------------------------------------------------

    #[test]
    fn is_valid_block_for_face_interior_and_exterior() {
        let bx = unit_box();
        let face = &bx.faces[0];
        let mut ctx = IntToolsContext::new();

        // A u-range inside the face's u-domain, sampled at v = 0.5.
        let inside = IntRange::new(0.2, 0.8).unwrap();
        assert!(ctx.is_valid_block_for_face(inside, face).unwrap());
        // The whole u-domain (endpoints On the boundary) is still valid.
        let whole = IntRange::new(0.0, 1.0).unwrap();
        assert!(ctx.is_valid_block_for_face(whole, face).unwrap());
        // A u-range that leaves the domain is rejected.
        let outside = IntRange::new(1.5, 2.5).unwrap();
        assert!(!ctx.is_valid_block_for_face(outside, face).unwrap());

        // Both faces: valid in both, invalid when one rejects.
        assert!(ctx
            .is_valid_block_for_faces(inside, &bx.faces[0], &bx.faces[1])
            .unwrap());
        assert!(!ctx
            .is_valid_block_for_faces(outside, &bx.faces[0], &bx.faces[1])
            .unwrap());
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Projection
    // -----------------------------------------------------------------------

    #[test]
    fn project_point_on_edge_returns_parameter() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        let edge = &bx.edges[0]; // (0,0,0) -> (1,0,0), params [0, 1]
        let t = ctx.project_point_on_edge(edge, &GpPnt::new(0.5, 0.0, 0.0)).expect("projection");
        assert!((t - 0.5).abs() < 1e-6, "param {t}");
        // A point off the line still projects to its closest parameter.
        let t2 = ctx.project_point_on_edge(edge, &GpPnt::new(0.25, 3.0, 0.0)).expect("projection");
        assert!((t2 - 0.25).abs() < 1e-6, "param {t2}");
        // An edge without geometry has no projection.
        assert!(ctx.project_point_on_edge(&Edge::new(), &GpPnt::zero()).is_none());
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn project_point_on_face_returns_uv() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        let face = &bx.faces[0];
        let (u, v) = ctx.project_point_on_face(face, &GpPnt::new(0.25, 0.75, 0.0)).unwrap();
        let surf = BRepTool::face_surface(face).expect("surface");
        let q = surf.d0(u, v);
        assert!(q.distance(&GpPnt::new(0.25, 0.75, 0.0)) < 1e-6, "uv {u},{v}");
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Geometric classification
    // -----------------------------------------------------------------------

    #[test]
    fn compute_pe_classifies_vertex_on_edge() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        // Vertex 0 (0,0,0) lies on edge 0 ((0,0,0) -> (1,0,0)).
        assert_eq!(ctx.compute_pe(&bx.vertices[0], &bx.edges[0], 1e-7), 0);
        // Vertex 6 (1,1,1) is farther than the tolerance from edge 0.
        assert_eq!(ctx.compute_pe(&bx.vertices[6], &bx.edges[0], 1e-7), -4);
        // A bare edge (no geometry) -> -2.
        assert_eq!(ctx.compute_pe(&bx.vertices[0], &Edge::new(), 1e-7), -2);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn compute_ve_classifies_edges() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        // edge 0 and edge 1 share the vertex (1,0,0) -> touch.
        assert_eq!(ctx.compute_ve(&bx.edges[0], &bx.edges[1], 1e-7), 0);
        // edge 0 and edge 4 are parallel disjoint (z=0 vs z=1) -> separated.
        assert_eq!(ctx.compute_ve(&bx.edges[0], &bx.edges[4], 1e-7), -4);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn compute_vf_classifies_vertex_on_face() {
        let bx = unit_box();
        let b = TopoBuilder::new();
        let interior = b.make_vertex(GpPnt::new(0.5, 0.5, 0.0), 1e-7);
        let above = b.make_vertex(GpPnt::new(0.5, 0.5, 2.0), 1e-7);
        let mut ctx = IntToolsContext::new();

        // Projection strictly inside the bottom face -> on.
        assert_eq!(ctx.compute_vf(&interior, &bx.faces[0], 1e-7), 0);
        // Distance too large -> -2.
        assert_eq!(ctx.compute_vf(&above, &bx.faces[0], 1e-7), -2);
        // Box corner vertex projects onto the face boundary -> -3.
        assert_eq!(ctx.compute_vf(&bx.vertices[0], &bx.faces[0], 1e-7), -3);
        clear_tree(&bx.solid.0);
        clear_tree(&interior.0);
        clear_tree(&above.0);
    }

    #[test]
    fn is_vertex_on_line_detects_hits() {
        let bx = unit_box();
        // A line along edge 0: through (0,0,0) direction (1,0,0).
        let lin = GpLin::from_pnt_dir(GpPnt::zero(), dir(1.0, 0.0, 0.0));
        let curve: Arc<dyn Curve> = Arc::new(GeomLine::new(lin));
        let ctx = IntToolsContext::new();
        // Vertex 0 is on the line.
        assert!(ctx.is_vertex_on_line(&bx.vertices[0], curve.as_ref(), 1e-7));
        // Vertex 6 (1,1,1) is off the line.
        assert!(!ctx.is_vertex_on_line(&bx.vertices[6], curve.as_ref(), 1e-7));
        clear_tree(&bx.solid.0);
    }

    // -----------------------------------------------------------------------
    // Misc
    // -----------------------------------------------------------------------

    #[test]
    fn uv_bounds_delegates_to_surface() {
        let bx = unit_box();
        let ctx = IntToolsContext::new();
        let (u0, u1, v0, v1) = ctx.uv_bounds(&bx.faces[0]);
        // The bottom face is an (unbounded) plane.
        assert!(!u0.is_finite() && !u1.is_finite() && !v0.is_finite() && !v1.is_finite());
        let (a0, a1, b0, b1) = crate::brep_surface::face_uv_bounds(&bx.faces[0]);
        assert_eq!((u0, u1, v0, v1), (a0, a1, b0, b1));
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn projection_tolerance_setter() {
        let mut ctx = IntToolsContext::new();
        assert_eq!(ctx.pon_s_projection_tolerance(), 1e-12);
        ctx.set_pon_s_projection_tolerance(1e-8);
        assert_eq!(ctx.pon_s_projection_tolerance(), 1e-8);
    }
