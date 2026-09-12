use super::prelude::*;
use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::geom::polygon_ops::polygon_area2d;
    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpDir};
    use occt_geom::{GeomCircle, GeomPlane};

    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::shape::{Compound, TopoShape};
    use crate::tgeometry::GeometryRegistry;
    use crate::topo_tools_full::faces_of;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&c);
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    /// Planar square face `origin + [0,size]·u_dir + [0,size]·v_dir`.
    fn make_square_face(b: &TopoBuilder, origin: GpPnt, u_dir: GpDir, v_dir: GpDir, size: f64) -> Face {
        let uv = GpVec::from_xyz(u_dir.xyz()).multiplied_scalar(size);
        let vv = GpVec::from_xyz(v_dir.xyz()).multiplied_scalar(size);
        let c0 = origin;
        let c1 = c0.translated_vec(&uv);
        let c2 = c1.translated_vec(&vv);
        let c3 = c0.translated_vec(&vv);
        let e1 = b.make_edge_segment(&c0, &c1);
        let e2 = b.make_edge_segment(&c1, &c2);
        let e3 = b.make_edge_segment(&c2, &c3);
        let e4 = b.make_edge_segment(&c3, &c0);
        let wire = b.make_wire(&[e1, e2, e3, e4]);
        let normal = GpVec::from_xyz(u_dir.xyz()).crossed(&GpVec::from_xyz(v_dir.xyz()));
        let n = GpDir::from_vec(&normal).expect("normal");
        let ax3 = GpAx3::new(origin, n, &u_dir).expect("frame");
        b.make_face(Arc::new(GeomPlane::new(GpPln::new(ax3))), &[wire])
    }

    #[test]
    fn crossing_line_edges_single_hit() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let hits = edge_edge_intersections(&e1, &e2, 1e-9);
        assert_eq!(hits.len(), 1, "hits: {hits:?}");
        assert!(hits[0].point.distance(&GpPnt::new(1.0, 1.0, 0.0)) < 1e-6);
        let d = 2.0f64.sqrt();
        assert!((hits[0].u1 - d).abs() < 1e-6, "u1={}", hits[0].u1);
        assert!((hits[0].u2 - d).abs() < 1e-6, "u2={}", hits[0].u2);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_edge_through_planar_face() {
        let b = TopoBuilder::new();
        let face = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let e = b.make_edge_segment(&GpPnt::new(0.5, 0.5, -1.0), &GpPnt::new(0.5, 0.5, 1.0));
        let hits = edge_face_intersections(&e, &face, 1e-9);
        assert_eq!(hits.len(), 1, "hits: {hits:?}");
        assert!(hits[0].1.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6);
        assert!((hits[0].0 - 1.0).abs() < 1e-6, "u={}", hits[0].0);
        clear_tree(&face.0);
        clear_tree(&e.0);
    }

    #[test]
    fn line_edge_missing_face() {
        let b = TopoBuilder::new();
        let face = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let e = b.make_edge_segment(&GpPnt::new(5.0, 5.0, -1.0), &GpPnt::new(5.0, 5.0, 1.0));
        assert!(edge_face_intersections(&e, &face, 1e-9).is_empty());
        clear_tree(&face.0);
        clear_tree(&e.0);
    }

    #[test]
    fn perpendicular_faces_intersection_segment() {
        let b = TopoBuilder::new();
        let f1 = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let f2 = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), 1.0);
        let segs = face_face_intersection_segments(&f1, &f2, 1e-9);
        assert_eq!(segs.len(), 1, "segs: {segs:?}");
        let (p, q) = segs[0];
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let c = GpPnt::new(1.0, 0.0, 0.0);
        assert!((p.distance(&a) < 1e-6 && q.distance(&c) < 1e-6) || (p.distance(&c) < 1e-6 && q.distance(&a) < 1e-6));
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn parallel_faces_no_segment() {
        let b = TopoBuilder::new();
        let f1 = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let f2 = make_square_face(&b, GpPnt::new(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        assert!(face_face_intersection_segments(&f1, &f2, 1e-9).is_empty());
        clear_tree(&f1.0);
        clear_tree(&f2.0);
    }

    #[test]
    fn point_on_face_inside_outside() {
        let b = TopoBuilder::new();
        let face = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        assert!(point_on_face(&face, &GpPnt::new(0.5, 0.5, 0.0), 1e-9));
        assert!(!point_on_face(&face, &GpPnt::new(2.0, 0.5, 0.0), 1e-9));
        assert!(!point_on_face(&face, &GpPnt::new(0.5, 0.5, 1.0), 1e-9));
        clear_tree(&face.0);
    }

    #[test]
    fn face_polygon_2d_unit_square() {
        let b = TopoBuilder::new();
        let face = make_square_face(&b, GpPnt::zero(), dir(1.0, 0.0, 0.0), dir(0.0, 1.0, 0.0), 1.0);
        let poly = face_polygon_2d(&face).expect("polygon");
        assert_eq!(poly.len(), 4);
        assert!((polygon_area2d(&poly).abs() - 1.0).abs() < 1e-9, "area={}", polygon_area2d(&poly));
        clear_tree(&face.0);
    }

    #[test]
    fn split_edge_at_params_wrapper() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(10.0, 0.0, 0.0));
        let parts = split_edges_at_params(&e, &[3.0, 7.0]);
        assert_eq!(parts.len(), 3);
        let ranges: Vec<(f64, f64)> = parts.iter().map(|p| BRepTool::edge_parameters(p)).collect();
        assert_eq!(ranges, vec![(0.0, 3.0), (3.0, 7.0), (7.0, 10.0)]);
        for p in &parts {
            clear_tree(&p.0);
        }
        clear_tree(&e.0);
    }

    #[test]
    fn line_circle_edges_coplanar_intersect() {
        let b = TopoBuilder::new();
        let ce = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let le = b.make_edge_segment(&GpPnt::new(-1.5, 0.0, 0.0), &GpPnt::new(1.5, 0.0, 0.0));
        let hits = edge_edge_intersections(&le, &ce, 1e-6);
        assert_eq!(hits.len(), 2, "hits: {hits:?}");
        for h in &hits {
            assert!((h.point.coord.modulus() - 1.0).abs() < 1e-6, "point {:?}", h.point);
        }
        clear_tree(&ce.0);
        clear_tree(&le.0);
    }

    #[test]
    fn circle_circle_edges_coplanar_intersect() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(1.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let hits = edge_edge_intersections(&c1, &c2, 1e-6);
        assert_eq!(hits.len(), 2, "hits: {hits:?}");
        for h in &hits {
            assert!((h.point.x() - 0.5).abs() < 1e-6, "x={}", h.point.x());
            assert!((h.point.y().abs() - 0.75f64.sqrt()).abs() < 1e-4, "y={}", h.point.y());
        }
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    fn find_face_by_normal(bx: &crate::primitives::BRepPrimBox, n: (f64, f64, f64)) -> Face {
        faces_of(&bx.solid.0)
            .into_iter()
            .find(|f| {
                face_plane_from_face(f).map_or(false, |pln| {
                    let d = GpVec::from_xyz(pln.axis().direction().xyz());
                    (d.x() - n.0).abs() < 1e-9 && (d.y() - n.1).abs() < 1e-9 && (d.z() - n.2).abs() < 1e-9
                })
            })
            .expect("box face with normal")
    }

    #[test]
    fn box_face_polygon_chains_shared_edges() {
        let bx = crate::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&bx.solid.0);
        assert_eq!(faces.len(), 6);
        for f in faces {
            let poly = face_polygon_2d(&f).expect("box face polygon");
            assert_eq!(poly.len(), 4, "poly: {poly:?}");
            assert!((polygon_area2d(&poly).abs() - 1.0).abs() < 1e-9, "area={}", polygon_area2d(&poly));
        }
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn adjacent_box_faces_intersect_along_edge() {
        let bx = crate::primitives::BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let fx = find_face_by_normal(&bx, (1.0, 0.0, 0.0));
        let fy = find_face_by_normal(&bx, (0.0, 1.0, 0.0));
        let segs = face_face_intersection_segments(&fx, &fy, 1e-9);
        assert_eq!(segs.len(), 1, "segs: {segs:?}");
        let (p, q) = segs[0];
        assert!((p.x() - 1.0).abs() < 1e-6 && (p.y() - 1.0).abs() < 1e-6);
        assert!((q.x() - 1.0).abs() < 1e-6 && (q.y() - 1.0).abs() < 1e-6);
        let zs = [p.z(), q.z()];
        assert!(zs.contains(&0.0) && zs.contains(&1.0), "z endpoints {zs:?}");
        clear_tree(&bx.solid.0);
    }

    #[test]
    fn collect_edge_intersection_points_across_shapes() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let mut comp = Compound::new();
        b.add_compound(&mut comp, &e2.0);
        let res = collect_edge_intersection_points(&[e1.clone()], &[comp.0.clone()], 1e-9);
        assert_eq!(res.len(), 1, "res: {res:?}");
        assert_eq!(res[0].0, 0);
        assert!(res[0].2.distance(&GpPnt::new(1.0, 1.0, 0.0)) < 1e-6);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
        clear_tree(&comp.0);
    }
