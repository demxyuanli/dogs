use super::prelude::*;
use super::*;
    use std::f64::consts::PI;

    use crate::primitives::{BRepPrimBox, BRepPrimCylinder, BRepPrimSphere};
    use crate::shape::Shell;
    use crate::shell_check::shell_is_closed;
    use crate::tgeometry::GeometryRegistry;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-6 * b.abs().max(1.0)
    }

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn z0_plane() -> GpPln {
        GpPln::new(GpAx3::standard())
    }

    #[test]
    fn offset_polygon_square_outward() {
        let plane = z0_plane();
        let pts = [
            GpPnt::new(-0.5, -0.5, 0.0),
            GpPnt::new(0.5, -0.5, 0.0),
            GpPnt::new(0.5, 0.5, 0.0),
            GpPnt::new(-0.5, 0.5, 0.0),
        ];
        let out = offset_polygon(&pts, &plane, 0.5).expect("offset");
        assert_eq!(out.len(), 4);
        for p in &out {
            assert!(approx(p.x().abs(), 1.0), "x {:?}", p);
            assert!(approx(p.y().abs(), 1.0), "y {:?}", p);
            assert!(approx(p.z(), 0.0));
        }
    }

    #[test]
    fn offset_polygon_square_inward() {
        let plane = z0_plane();
        let pts = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(1.0, -1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(-1.0, 1.0, 0.0),
        ];
        let out = offset_polygon(&pts, &plane, -0.5).expect("offset");
        assert_eq!(out.len(), 4);
        for p in &out {
            assert!(approx(p.x().abs(), 0.5), "x {:?}", p);
            assert!(approx(p.y().abs(), 0.5), "y {:?}", p);
        }
    }

    #[test]
    fn offset_circle_radius() {
        let plane = z0_plane();
        let c = offset_circle(GpPnt::zero(), 1.0, &plane, 0.5).expect("offset +");
        assert!(approx(c.radius(), 1.5));
        let c2 = offset_circle(GpPnt::zero(), 1.0, &plane, -0.5).expect("offset -");
        assert!(approx(c2.radius(), 0.5));
    }

    #[test]
    fn offset_wire_2d_closed_square() {
        let b = TopoBuilder::new();
        let pts = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(1.0, -1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(-1.0, 1.0, 0.0),
        ];
        let edges = [
            b.make_edge_segment(&pts[0], &pts[1]),
            b.make_edge_segment(&pts[1], &pts[2]),
            b.make_edge_segment(&pts[2], &pts[3]),
            b.make_edge_segment(&pts[3], &pts[0]),
        ];
        let wire = b.make_wire(&edges);
        let out = offset_wire_2d(&wire, &z0_plane(), 0.3).expect("offset");
        for v in vertices_of(&out.0) {
            let p = BRepTool::vertex_point(&v);
            assert!(approx(p.x().abs(), 1.3), "x {}", p.x());
            assert!(approx(p.y().abs(), 1.3), "y {}", p.y());
        }
        clear_tree(&out.0);
    }

    #[test]
    fn offset_face_plane_translates() {
        let pts = [
            GpPnt::new(-1.0, -1.0, 0.0),
            GpPnt::new(1.0, -1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(-1.0, 1.0, 0.0),
        ];
        let face = crate::brep_builder_api::make_face_from_polygon(&pts).expect("square face");
        let off = offset_face(&face, 1.0).expect("offset");
        let s = BRepTool::face_surface(&off).expect("surface");
        let p = s.d0(0.0, 0.0);
        assert!(approx(p.z(), 1.0), "z {}", p.z());
        clear_tree(&off.0);
    }

    #[test]
    fn offset_face_sphere_grows() {
        let s = BRepPrimSphere::make_sphere(1.0);
        let face = faces_of(&s.solid.0)[0].clone();
        let off = offset_face(&face, 0.5).expect("offset");
        let surf = BRepTool::face_surface(&off).expect("surface");
        let c = sphere_center(surf.as_ref()).expect("center");
        assert!(approx(c.distance(&GpPnt::zero()), 0.0), "center {:?}", c);
        let r = surf.d0(0.0, 0.0).distance(&c);
        assert!(approx(r, 1.5), "radius {}", r);
        clear_tree(&off.0);
    }

    #[test]
    fn offset_face_cylinder_radius() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 3.0);
        let lateral = faces_of(&c.solid.0)
            .into_iter()
            .find(|f| !face_is_planar(f))
            .expect("lateral face");
        let off = offset_face(&lateral, 0.2).expect("offset");
        let surf = BRepTool::face_surface(&off).expect("surface");
        let (_, _, r) = extract_cylinder(surf.as_ref()).expect("cylinder params");
        assert!(approx(r, 1.2), "radius {}", r);
        clear_tree(&off.0);
    }

    #[test]
    fn offset_solid_box_grows() {
        let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let result = offset_shell(&b.solid.0, 0.5).expect("offset");
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        let mut zs = Vec::new();
        for v in vertices_of(&result) {
            let p = BRepTool::vertex_point(&v);
            xs.push(p.x());
            ys.push(p.y());
            zs.push(p.z());
        }
        let (min, max) = (
            |v: &Vec<f64>| v.iter().cloned().fold(f64::INFINITY, f64::min),
            |v: &Vec<f64>| v.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        );
        assert!(approx(min(&xs), -0.5) && approx(max(&xs), 2.5), "x [{},{}]", min(&xs), max(&xs));
        assert!(approx(min(&ys), -0.5) && approx(max(&ys), 2.5));
        assert!(approx(min(&zs), -0.5) && approx(max(&zs), 2.5));
        // Side length 3.0.
        assert!(approx(max(&xs) - min(&xs), 3.0));
        // Closed manifold shell.
        let shell = Shell(result.tshape.read().unwrap().children[0].clone());
        assert!(shell_is_closed(&shell), "offset box must be a closed shell");
        clear_tree(&result);
    }

    #[test]
    fn offset_solid_box_shrinks() {
        let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let result = offset_shell(&b.solid.0, -0.4).expect("offset");
        let xs: Vec<f64> = vertices_of(&result)
            .iter()
            .map(|v| BRepTool::vertex_point(v).x())
            .collect();
        let ys: Vec<f64> = vertices_of(&result)
            .iter()
            .map(|v| BRepTool::vertex_point(v).y())
            .collect();
        let mn = |v: &Vec<f64>| v.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = |v: &Vec<f64>| v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(approx(mn(&xs), 0.4) && approx(mx(&xs), 1.6), "x [{},{}]", mn(&xs), mx(&xs));
        assert!(approx(mn(&ys), 0.4) && approx(mx(&ys), 1.6));
        assert!(approx(mx(&xs) - mn(&xs), 1.2), "side length");
        clear_tree(&result);
    }

    #[test]
    fn offset_negative_shrinks_sphere() {
        let s = BRepPrimSphere::make_sphere(1.0);
        let result = offset_shell(&s.solid.0, -0.3).expect("offset");
        let face = faces_of(&result)[0].clone();
        let surf = BRepTool::face_surface(&face).expect("surface");
        let c = sphere_center(surf.as_ref()).expect("center");
        let r = surf.d0(0.0, 0.0).distance(&c);
        assert!(approx(r, 0.7), "radius {}", r);
        clear_tree(&result);
    }

    #[test]
    fn plane_plane_plane_corner() {
        let px = GpPln::new(
            GpAx3::new(GpPnt::new(1.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap(), &GpDir::new(0.0, 1.0, 0.0).unwrap())
                .expect("x plane"),
        );
        let py = GpPln::new(
            GpAx3::new(GpPnt::new(0.0, 1.0, 0.0), GpDir::new(0.0, 1.0, 0.0).unwrap(), &GpDir::new(0.0, 0.0, 1.0).unwrap())
                .expect("y plane"),
        );
        let pz = GpPln::new(
            GpAx3::new(GpPnt::new(0.0, 0.0, 1.0), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
                .expect("z plane"),
        );
        let pt = plane_plane_plane_intersection(&px, &py, &pz).expect("corner");
        assert!(pt.distance(&GpPnt::new(1.0, 1.0, 1.0)) < 1e-9, "{:?}", pt);
        // Two parallel planes (x=1 and x=2) with z=1: no common point.
        let px2 = GpPln::new(
            GpAx3::new(GpPnt::new(2.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap(), &GpDir::new(0.0, 1.0, 0.0).unwrap())
                .expect("x=2 plane"),
        );
        assert!(plane_plane_plane_intersection(&px, &px2, &pz).is_none());
    }

    #[test]
    fn offset_wire_3d_circle() {
        let b = TopoBuilder::new();
        let ax3 = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
            .expect("circle frame");
        let curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(GpCirc::new(ax3.ax2(), 1.0)));
        let mut e = b.make_edge(curve.clone(), 0.0, 2.0 * PI);
        let seam = b.make_vertex(curve.d0(0.0), 0.0);
        b.add(&mut e.0, &seam.0);
        b.add(&mut e.0, &seam.0);
        let wire = b.make_wire(&[e]);
        let out = offset_wire_3d(&wire, 0.2).expect("offset");
        for v in vertices_of(&out.0) {
            let p = BRepTool::vertex_point(&v);
            assert!(approx(p.distance(&GpPnt::zero()), 1.2), "radius {}", p.distance(&GpPnt::zero()));
        }
        clear_tree(&out.0);
    }
