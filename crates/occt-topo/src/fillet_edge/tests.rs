use super::prelude::*;
use super::*;
    use crate::brep_surface::classify_surface;
    use crate::primitives::BRepPrimBox;
    use crate::shape::Shell;
    use crate::shell_check::shell_is_closed;
    use crate::topo_tools_full::vertices_of;

    fn box_edges(b: &BRepPrimBox) -> Vec<Edge> {
        edges_of(&b.solid.0)
    }

    fn find_face_by_surface<'a>(
        faces: &'a [Face],
        f: impl Fn(&dyn Surface) -> bool,
    ) -> Option<&'a Face> {
        faces
            .iter()
            .find(|fa| BRepTool::face_surface(fa).map(|s| f(s.as_ref())).unwrap_or(false))
    }

    #[test]
    fn corner_wire_fillet_alias() {
        // The 2D corner fillet module and this 3D edge fillet coexist without
        // a name clash: `fillet::fillet_corner` (wire) vs `fillet_edge` here.
        assert!(FilletSpec { radius: 0.5 }.check().is_ok());
        assert!(FilletSpec { radius: 0.0 }.check().is_err());
    }

    #[test]
    fn fillet_box_edge_90deg() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        // Edge from (0,0,0) → (2,0,0): the X edge at y=0, z=0.
        let e = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .expect("find X edge");
        let out = fillet_edge(&bx.solid.0, e, 0.4).expect("fillet");
        let faces = faces_of(&out);
        // 6 → 7 faces (blend cylinder added).
        assert_eq!(faces.len(), 7, "face count {}", faces.len());
        // Blend face classifies as a cylinder.
        let blend = find_face_by_surface(&faces, |s| classify_surface_full(s) == SurfaceKind::Cylinder)
            .expect("blend face");
        let r = cylinder_radius(BRepTool::face_surface(blend).unwrap().as_ref()).unwrap();
        assert!((r - 0.4).abs() < 1e-6, "blend radius {r}");
        // Closed shell.
        let shell = Shell(out.tshape.read().unwrap().children[0].clone());
        assert!(shell_is_closed(&shell), "filleted box is closed");
        // Every edge endpoint of the result lies on the original box surface
        // (within tolerance). The box spans [0,2]³; the surface is the union
        // of the six coordinate planes at 0 and 2.
        let on_box_surface = |p: &GpPnt| {
            let d = [
                p.x(), p.x() - 2.0, p.y(), p.y() - 2.0, p.z(), p.z() - 2.0,
            ];
            d.iter().map(|v| v.abs()).fold(f64::INFINITY, f64::min) < 1e-6
        };
        for edge in edges_of(&out) {
            let (a, b) = BRepTool::edge_vertices(&edge).unwrap();
            assert!(on_box_surface(&a), "endpoint {a:?} off box surface");
            assert!(on_box_surface(&b), "endpoint {b:?} off box surface");
        }
    }

    #[test]
    fn fillet_increases_face_count() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let out = fillet_edge(&bx.solid.0, e, 0.3).expect("fillet");
        assert_eq!(faces_of(&out).len(), 7);
    }

    #[test]
    fn fillet_blend_surface_is_cylinder() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let out = fillet_edge(&bx.solid.0, e, 0.4).unwrap();
        let faces = faces_of(&out);
        let blend = find_face_by_surface(&faces, |s| classify_surface_full(s) == SurfaceKind::Cylinder)
            .expect("blend face");
        // Both classifiers see the cylinder: the analytic one reports the exact
        // `GetType()` kind (`Geom_CylindricalSurface` → `Cylinder`) and the
        // extended wrapper agrees. (The analytic classifier used to be blind to
        // cylinders and reported `Other`.)
        let blend_surf = BRepTool::face_surface(blend).unwrap();
        assert_eq!(classify_surface(blend_surf.as_ref()), SurfaceKind::Cylinder);
        assert_eq!(classify_surface_full(blend_surf.as_ref()), SurfaceKind::Cylinder);
    }

    #[test]
    fn fillet_chain_two_edges() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        // Edge 0: (0,0,0)→(2,0,0). Edge at index of the opposite top edge
        // (2,2,2)→(0,2,2): find by endpoints.
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
        let out = fillet_edge_chain(&bx.solid.0, &[e0, e1], 0.3).expect("chain");
        assert_eq!(faces_of(&out).len(), 8, "two fillets add two faces");
        let shell = Shell(out.tshape.read().unwrap().children[0].clone());
        assert!(shell_is_closed(&shell), "chained fillet is closed");
    }

    #[test]
    fn fillet_corner_solid_sphere_blend() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let verts = vertices_of(&bx.solid.0);
        let corner = verts
            .iter()
            .find(|v| BRepTool::vertex_point(v).is_equal(&GpPnt::new(0.0, 0.0, 0.0)))
            .expect("corner vertex");
        let out = fillet_corner_solid(&bx.solid.0, corner, 0.4).expect("corner fillet");
        let faces = faces_of(&out);
        assert_eq!(faces.len(), 7, "corner fillet adds one face");
        let sphere = find_face_by_surface(&faces, |s| classify_surface(s) == SurfaceKind::Sphere)
            .expect("sphere blend face");
        let _ = sphere;
        let shell = Shell(out.tshape.read().unwrap().children[0].clone());
        assert!(shell_is_closed(&shell), "corner fillet is closed");
    }

    #[test]
    fn nonplanar_edge_errors() {
        let cyl = crate::primitives::BRepPrimCylinder::make_cylinder(1.0, 2.0);
        let edges = edges_of(&cyl.solid.0);
        // The cap circles are adjacent to a planar cap and the lateral face;
        // the seam is adjacent to the lateral (non-planar) face on both sides.
        // Filleting the seam must fail because both adjacent faces are the
        // same non-planar lateral surface.
        let seam = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.x().abs() > 0.99 && b.x().abs() > 0.99 && (a.z() - b.z()).abs() > 0.1
            })
            .unwrap();
        assert!(fillet_edge(&cyl.solid.0, seam, 0.2).is_err());
    }

    #[test]
    fn radius_positive_required() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e = &edges[0];
        assert!(fillet_edge(&bx.solid.0, e, 0.0).is_err());
        assert!(fillet_edge(&bx.solid.0, e, -1.0).is_err());
    }

    #[test]
    fn vertex_endpoint_geometry_preserved() {
        let bx = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let edges = box_edges(&bx);
        let e = edges
            .iter()
            .find(|e| {
                let (a, b) = BRepTool::edge_vertices(e).unwrap();
                a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
            })
            .unwrap();
        let out = fillet_edge(&bx.solid.0, e, 0.4).unwrap();
        // The two far endpoints of the filleted edge are the +Y/+Z corners on
        // the original surface; the tangency vertices lie on the adjacent
        // planes. Check the far corner (2,2,0) and (0,2,0) remain present.
        let verts = vertices_of(&out);
        let ps: Vec<GpPnt> = verts.iter().map(BRepTool::vertex_point).collect();
        for want in [
            GpPnt::new(2.0, 2.0, 0.0),
            GpPnt::new(0.0, 2.0, 0.0),
            GpPnt::new(2.0, 0.0, 2.0),
            GpPnt::new(0.0, 0.0, 2.0),
        ] {
            assert!(
                ps.iter().any(|p| p.distance(&want) < 1e-6),
                "missing preserved vertex {want:?}"
            );
        }
        // The original sharp-edge endpoints are replaced by tangency points on
        // the adjacent planes (z=0 and y=0 for this edge).
        assert!(
            ps.iter().any(|p| (p.distance(&GpPnt::new(0.4, 0.0, 0.0)) < 1e-6)
                || (p.distance(&GpPnt::new(0.0, 0.4, 0.0)) < 1e-6)),
            "tangency vertices present"
        );
    }
