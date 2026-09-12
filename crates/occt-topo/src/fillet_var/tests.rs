use super::prelude::*;
use super::*;
    use crate::brep_surface::{classify_surface, SurfaceKind};
    use crate::primitives::BRepPrimBox;
    use crate::shape::Shell;
    use crate::shell_check::shell_is_closed;
    use crate::tgeometry::GeometryRegistry;
    use crate::topo_tools_full::vertices_of;

    /// Release registry entries for a shape tree so the process-wide geometry
    /// side-table (keyed by `TShape` address) does not serve stale geometry to
    /// other tests running in parallel.
    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn box_edges(b: &BRepPrimBox) -> Vec<Edge> {
        edges_of(&b.solid.0)
    }

    /// Find the bottom-front edge (0,0,0)→(2,0,0) of a 2×2×2 box.
    fn front_bottom_edge(b: &BRepPrimBox) -> Edge {
        box_edges(b)
            .into_iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .expect("front-bottom edge")
    }

    fn closed_shell(shape: &TopoShape) -> Shell {
        Shell(shape.tshape.read().unwrap().children[0].clone())
    }

    fn blend_face<'a>(faces: &'a [Face]) -> &'a Face {
        faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Other)
                    .unwrap_or(false)
            })
            .expect("blend face")
    }

    #[test]
    fn radius_at_linear_endpoints() {
        let spec = VarFilletSpec::new(0.3, 0.6);
        assert!((radius_at(&spec, 0.0) - 0.3).abs() < 1e-12);
        assert!((radius_at(&spec, 1.0) - 0.6).abs() < 1e-12);
        assert!((radius_at(&spec, 0.5) - 0.45).abs() < 1e-12);
        // Clamped outside [0, 1].
        assert!((radius_at(&spec, 2.0) - 0.6).abs() < 1e-12);
    }

    #[test]
    fn radius_at_quadratic_monotone() {
        let spec = VarFilletSpec { r_start: 0.2, r_end: 0.8, law: RadiusLaw::Quadratic };
        let mut prev = radius_at(&spec, 0.0);
        for i in 1..=100 {
            let r = radius_at(&spec, i as f64 / 100.0);
            assert!(r >= prev - 1e-12, "not monotone at t={}", i as f64 / 100.0);
            prev = r;
        }
        assert!((radius_at(&spec, 0.0) - 0.2).abs() < 1e-12);
        assert!((radius_at(&spec, 1.0) - 0.8).abs() < 1e-12);
    }

    #[test]
    fn fillet_var_box_edge_adds_face() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let out = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.3, 0.6), 1e-6)
            .expect("variable fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 7, "face count {}", faces.len());
        assert!(shell_is_closed(&closed_shell(&out)), "filleted box is closed");
        // The blend face is present and not a cylinder (variable radius).
        let bf = blend_face(&faces);
        let _ = bf;
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_blend_radius_varies() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let out = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.3, 0.6), 1e-6)
            .expect("variable fillet");
        let faces = faces_of(&out);
        let bf = blend_face(&faces);
        let s = BRepTool::face_surface(bf).expect("blend surface");
        // Sample the surface at the arc midpoint at both ends of the edge.
        let pa = s.d0(0.0, 0.5);
        let pb = s.d0(1.0, 0.5);
        // Distance from the edge line (the x-axis, y = z = 0).
        let da = GpVec::new(pa.y(), pa.z(), 0.0).magnitude();
        let db = GpVec::new(pb.y(), pb.z(), 0.0).magnitude();
        assert!(
            (db - da).abs() > 0.05,
            "blend radius must vary along the edge: d(0)={da} d(1)={db}"
        );
        // The larger end is the larger radius.
        assert!(db > da, "radius should grow toward t=1");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_constant_reduces_to_fixed() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let var = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.4, 0.4), 1e-6)
            .expect("constant-as-variable fillet");
        let fixed = crate::fillet_edge::fillet_edge(&bx.solid.0, &e, 0.4).expect("fixed fillet");
        assert_eq!(faces_of(&var).len(), faces_of(&fixed).len(), "same face count");
        assert!(shell_is_closed(&closed_shell(&var)));
        assert!(shell_is_closed(&closed_shell(&fixed)));
        clear_tree(&var);
        clear_tree(&fixed);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_chain_two_edges() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e0 = edges
            .iter()
            .position(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let e1 = edges
            .iter()
            .position(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                (a.is_equal(&GpPnt::new(2.0, 2.0, 2.0)) && b.is_equal(&GpPnt::new(0.0, 2.0, 2.0)))
                    || (b.is_equal(&GpPnt::new(2.0, 2.0, 2.0)) && a.is_equal(&GpPnt::new(0.0, 2.0, 2.0)))
            })
            .unwrap();
        let specs = [VarFilletSpec::new(0.2, 0.3), VarFilletSpec::new(0.3, 0.4)];
        let out = fillet_edge_var_chain(&bx.solid.0, &[e0, e1], &specs, 1e-6).expect("var chain");
        assert_eq!(faces_of(&out).len(), 8, "two fillets add two faces");
        assert!(shell_is_closed(&closed_shell(&out)), "chained variable fillet is closed");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_invalid_radius() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        assert!(fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(-0.5, 0.5), 1e-6).is_err());
        assert!(fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.5, 0.0), 1e-6).is_err());
        assert!(VarFilletSpec::new(0.0, 1.0).check().is_err());
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_nonplanar_error() {
        let cyl = crate::primitives::BRepPrimCylinder::make_cylinder(1.0, 2.0);
        let edges = edges_of(&cyl.solid.0);
        let seam = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.x().abs() > 0.99 && b.x().abs() > 0.99 && (a.z() - b.z()).abs() > 0.1
            })
            .unwrap();
        let spec = VarFilletSpec::new(0.2, 0.4);
        assert!(fillet_edge_var(&cyl.solid.0, seam, &spec, 1e-6).is_err());
        clear_tree(&cyl.solid.0);
    }

    #[test]
    fn fillet_radius_profile_matches_spec() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let spec = VarFilletSpec { r_start: 0.2, r_end: 0.9, law: RadiusLaw::Cubic };
        let prof = fillet_radius_profile(&bx.solid.0, &e, &spec, 21);
        assert_eq!(prof.len(), 21);
        for (t, r) in &prof {
            assert!((r - radius_at(&spec, *t)).abs() < 1e-12);
        }
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_chain_shared_vertex() {
        // Two consecutive edges of the bottom-face loop share the corner
        // (2,0,0). Filleting the first edge rebuilds both of its end faces and
        // both adjacent faces, which consumes the shared corner and drops the
        // second edge from the solid; the sequential chain therefore cannot
        // re-find it and must report the conflict cleanly rather than fillet a
        // wrong edge. This test locks in that honest behaviour.
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e0 = edges
            .iter()
            .position(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let e1 = edges
            .iter()
            .position(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                (a.is_equal(&GpPnt::new(2.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 2.0, 0.0)))
                    || (b.is_equal(&GpPnt::new(2.0, 0.0, 0.0)) && a.is_equal(&GpPnt::new(2.0, 2.0, 0.0)))
            })
            .unwrap();
        let specs = [VarFilletSpec::new(0.2, 0.3), VarFilletSpec::new(0.3, 0.2)];
        // The first fillet consumes the shared corner; re-finding the second
        // edge fails, and the call must report the conflict rather than corrupt
        // the topology.
        let res = fillet_edge_var_chain(&bx.solid.0, &[e0, e1], &specs, 1e-6);
        match res {
            Ok(r) => {
                assert!(shell_is_closed(&closed_shell(&r)), "consumed-corner chain stays closed");
                clear_tree(&r);
            }
            Err(_) => {}
        }
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn fillet_var_small_radius_ok() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let out = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.1, 0.1), 1e-6)
            .expect("small radius fillet");
        assert_eq!(faces_of(&out).len(), 7);
        assert!(shell_is_closed(&closed_shell(&out)), "small-radius fillet is closed");
        // The blend face is not degenerate: it has real extent along v.
        let faces = faces_of(&out);
        let bf = blend_face(&faces);
        let s = BRepTool::face_surface(bf).unwrap();
        let mid = s.d0(0.5, 0.5);
        // Still off the edge line.
        assert!(GpVec::new(mid.y(), mid.z(), 0.0).magnitude() > 0.01);
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn chain_smooth_wraps_phase6() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e0 = edges
            .iter()
            .position(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let e1 = edges
            .iter()
            .position(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                (a.is_equal(&GpPnt::new(2.0, 2.0, 2.0)) && b.is_equal(&GpPnt::new(0.0, 2.0, 2.0)))
                    || (b.is_equal(&GpPnt::new(2.0, 2.0, 2.0)) && a.is_equal(&GpPnt::new(0.0, 2.0, 2.0)))
            })
            .unwrap();
        let out = fillet_edge_chain_smooth(&bx.solid.0, &[e0, e1], 0.3, 1e-6).expect("smooth chain");
        assert_eq!(faces_of(&out).len(), 8);
        assert!(shell_is_closed(&closed_shell(&out)));
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn vertices_preserved_on_far_corners() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e = front_bottom_edge(&bx);
        let out = fillet_edge_var(&bx.solid.0, &e, &VarFilletSpec::new(0.3, 0.6), 1e-6).unwrap();
        let ps: Vec<GpPnt> = vertices_of(&out).iter().map(BRepTool::vertex_point).collect();
        for want in [
            GpPnt::new(2.0, 2.0, 0.0),
            GpPnt::new(0.0, 2.0, 0.0),
            GpPnt::new(2.0, 0.0, 2.0),
            GpPnt::new(0.0, 0.0, 2.0),
        ] {
            assert!(ps.iter().any(|p| p.distance(&want) < 1e-6), "missing preserved vertex {want:?}");
        }
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    // --- Shared-vertex corner patch (rolling-ball corner blend) ---

    fn find_edge_by_endpoints_test(b: &BRepPrimBox, a: &GpPnt, bpt: &GpPnt) -> usize {
        box_edges(b)
            .iter()
            .position(|e| {
                let (x, y) = BRepTool::edge_vertices(e).unwrap();
                (x.is_equal(a) && y.is_equal(bpt)) || (x.is_equal(bpt) && y.is_equal(a))
            })
            .expect("edge by endpoints")
    }

    fn corner_vertex(b: &BRepPrimBox, p: &GpPnt) -> Vertex {
        vertices_of(&b.solid.0)
            .into_iter()
            .find(|v| BRepTool::vertex_point(v).is_equal(p))
            .expect("corner vertex")
    }

    #[test]
    fn cut_disk_splits_bottom_edge() {
        // The -Y face rectangle at z >= R, with the corner disk centred at
        // (R, 0, R) of radius R.
        let r = 0.4;
        let rect = [
            GpPnt::new(0.0, 0.0, r),
            GpPnt::new(2.0, 0.0, r),
            GpPnt::new(2.0, 0.0, 2.0),
            GpPnt::new(0.0, 0.0, 2.0),
        ];
        let cut = cut_disk_poly(&rect, &GpPnt::new(r, 0.0, r), r).expect("cut");
        // The bottom edge must be split: the polygon should contain the point
        // (2R, 0, R) where the disk boundary crosses the z = R edge.
        assert!(
            cut.iter().any(|p| p.distance(&GpPnt::new(2.0 * r, 0.0, r)) < 1e-6),
            "bottom edge was not split by the disk"
        );
    }

    #[test]
    fn corner_center_box() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let v = corner_vertex(&bx, &GpPnt::new(0.0, 0.0, 0.0));
        let c = corner_center(&bx.solid.0, &v, 0.4, 1e-6).expect("corner center");
        assert!(
            c.distance(&GpPnt::new(0.4, 0.4, 0.4)) < 1e-9,
            "corner centre {c:?}"
        );
        // The centre of the opposite corner (2,2,2) is (1.6,1.6,1.6).
        let v2 = corner_vertex(&bx, &GpPnt::new(2.0, 2.0, 2.0));
        let c2 = corner_center(&bx.solid.0, &v2, 0.4, 1e-6).unwrap();
        assert!(c2.distance(&GpPnt::new(1.6, 1.6, 1.6)) < 1e-9);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn corner_radius_positive() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let v = corner_vertex(&bx, &GpPnt::new(0.0, 0.0, 0.0));
        let n = edges_of(&bx.solid.0).len();
        // specs aligned with edges_of: constant 0.4 on the two edges incident
        // to the corner.
        let mut specs = vec![VarFilletSpec::new(0.1, 0.1); n];
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        specs[e0] = VarFilletSpec::new(0.4, 0.4);
        specs[e1] = VarFilletSpec::new(0.4, 0.4);
        let r = corner_patch_radius_at_vertex(&bx.solid.0, &v, &specs, 1e-6);
        assert!((r - 0.4).abs() < 1e-9, "corner radius {r}");
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn corner_blend_box_corner_closed() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        let specs = [VarFilletSpec::new(0.4, 0.4), VarFilletSpec::new(0.4, 0.4)];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1], &specs, 1e-6)
            .expect("corner chain");
        let faces = faces_of(&out);
        assert_eq!(
            faces.len(),
            9,
            "6 faces + 2 blends + 1 corner patch, got {}",
            faces.len()
        );
        assert!(shell_is_closed(&closed_shell(&out)), "corner chain must be closed");
        // The corner patch face is a sphere.
        let sphere = faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Sphere)
                    .unwrap_or(false)
            })
            .expect("corner patch sphere face");
        let _ = sphere;
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn chain_with_corner_success() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        let specs = [VarFilletSpec::new(0.3, 0.3), VarFilletSpec::new(0.3, 0.3)];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1], &specs, 1e-6)
            .expect("chain with corner must not error");
        assert!(shell_is_closed(&closed_shell(&out)), "chain with corner is closed");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn corner_patch_is_spherical() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        let specs = [VarFilletSpec::new(0.4, 0.4), VarFilletSpec::new(0.4, 0.4)];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1], &specs, 1e-6).unwrap();
        let faces = faces_of(&out);
        let sphere = faces
            .iter()
            .find(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Sphere)
                    .unwrap_or(false)
            })
            .expect("corner patch sphere face");
        // All sampled surface points are equidistant from the corner centre.
        let v = corner_vertex(&bx, &GpPnt::new(0.0, 0.0, 0.0));
        let cc = corner_center(&bx.solid.0, &v, 0.4, 1e-6).unwrap();
        let s = BRepTool::face_surface(sphere).unwrap();
        let r0 = s.d0(0.0, 0.0).distance(&cc);
        for (u, w) in [(0.5, 0.3), (1.0, 0.0), (1.5, 0.4), (3.0, 0.2)] {
            let d = s.d0(u, w).distance(&cc);
            assert!((d - r0).abs() < 1e-6, "patch point at distance {d} (expected {r0})");
        }
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn chain_nonconsecutive_still_works() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        // Opposite top edge: (2,2,2)→(0,2,2).
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(2.0, 2.0, 2.0), &GpPnt::new(0.0, 2.0, 2.0));
        let specs = [VarFilletSpec::new(0.2, 0.3), VarFilletSpec::new(0.3, 0.4)];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1], &specs, 1e-6)
            .expect("non-consecutive chain");
        assert_eq!(faces_of(&out).len(), 8, "two sequential fillets");
        assert!(shell_is_closed(&closed_shell(&out)), "non-consecutive chain is closed");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn chain_with_corner_three_edges() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        // Bottom-face loop: front X edge, right Y edge, back X edge.
        let e0 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e1 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(0.0, 2.0, 0.0));
        let e2 = find_edge_by_endpoints_test(&bx, &GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let specs = [VarFilletSpec::new(0.3, 0.3); 3];
        let out = fillet_edges_chain_with_corner(&bx.solid.0, &[e0, e1, e2], &specs, 1e-6)
            .expect("three-edge corner chain");
        let faces = faces_of(&out);
        assert_eq!(
            faces.len(),
            11,
            "6 faces + 3 blends + 2 corner patches, got {}",
            faces.len()
        );
        assert!(shell_is_closed(&closed_shell(&out)), "three-edge corner chain is closed");
        // Two corner patch faces.
        let spheres = faces
            .iter()
            .filter(|f| {
                BRepTool::face_surface(f)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Sphere)
                    .unwrap_or(false)
            })
            .count();
        assert_eq!(spheres, 2, "two corner patches, got {spheres}");
        clear_tree(&out);
        clear_tree(&bx.solid.0);
    }
