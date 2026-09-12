use super::prelude::*;
use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::builder_face::build_face_with_holes;
    use crate::topo_tools_full::wire_is_closed;
    use occt_core::gp::{GpAx3, GpDir, GpLin, GpPln, GpPnt};
    use occt_geom::GeomLine;

    #[test]
    fn make_new_vertex_at_box_corner() {
        let b = unit_box();
        let v = AlgoTools::make_new_vertex(&b.corners[0], 1e-7).expect("vertex");
        assert!(v.is_vertex());
        let vtx = Vertex(v);
        assert!(BRepTool::vertex_point(&vtx).is_equal(&b.corners[0]), "point must match corner");
        assert_eq!(BRepTool::vertex_tolerance(&vtx), 1e-7);
    }

    #[test]
    fn make_edge_from_box_edge_curve_connects_endpoints() {
        let b = unit_box();
        let curve = BRepTool::edge_curve(&b.edges[0]).expect("box edge curve");
        let e = AlgoTools::make_edge(
            curve,
            Some(&b.vertices[0].0),
            0.0,
            Some(&b.vertices[1].0),
            1.0,
            1e-7,
        )
        .expect("edge");
        assert!(e.is_edge());
        assert_eq!(e.tshape.read().unwrap().children.len(), 2);

        let p0 = AlgoTools::point_on_edge(&Edge(e.clone()), 0.0).expect("start point");
        let p1 = AlgoTools::point_on_edge(&Edge(e.clone()), 1.0).expect("end point");
        assert!(p0.is_equal(&b.corners[0]), "start {:?}", p0);
        assert!(p1.is_equal(&b.corners[1]), "end {:?}", p1);
        assert_eq!(BRepTool::edge_tolerance(&Edge(e)), 1e-7);

        // Vertex tolerances were bumped to cover the curve ends.
        assert!(
            BRepTool::vertex_tolerance(&b.vertices[0]) >= 1e-7,
            "v0 tolerance grew"
        );
    }

    #[test]
    fn make_pcurve_agrees_with_pcurve_full() {
        let b = unit_box();
        let edge = &b.edges[0];
        let face = &b.faces[0];
        let a = AlgoTools::make_pcurve(edge, face).expect("make_pcurve");
        let full = pcurve_full::make_pcurve_full(edge, face).expect("pcurve_full");
        for i in 0..=8 {
            let t = i as f64 / 8.0;
            let pa = a.d0(t);
            let pf = full.d0(t);
            assert!(
                (pa.x() - pf.x()).abs() < 1e-9 && (pa.y() - pf.y()).abs() < 1e-9,
                "sample t={t}"
            );
        }
    }

    #[test]
    fn compute_vv_coincident_returns_one_separated_returns_zero() {
        let b = unit_box();
        // The box corner vertex coincides with its own point → 1.
        assert_eq!(AlgoTools::compute_vv(&b.vertices[0].0, &b.corners[0], 1e-7), 1);
        // A far-away point is separated → 0.
        assert_eq!(
            AlgoTools::compute_vv(&b.vertices[0].0, &GpPnt::new(10.0, 10.0, 10.0), 1e-7),
            0
        );
    }

    #[test]
    fn point_on_edge_at_endpoints() {
        let b = unit_box();
        let p0 = AlgoTools::point_on_edge(&b.edges[0], 0.0).expect("t=0");
        let p1 = AlgoTools::point_on_edge(&b.edges[0], 1.0).expect("t=1");
        assert!(p0.is_equal(&b.corners[0]), "t=0 {:?}", p0);
        assert!(p1.is_equal(&b.corners[1]), "t=1 {:?}", p1);
    }

    #[test]
    fn update_vertex_grows_tolerance_to_cover_edge_point() {
        let b = unit_box();
        // Vertex near (0,0,0.5): edge 0 runs (0,0,0)→(1,0,0), so the point at
        // t=0.5 is (0.5,0,0) — distance ≈ 0.707 → tolerance must grow above it.
        let mut v = TopoBuilder::new().make_vertex(GpPnt::new(0.0, 0.0, 0.5), 1e-7);
        AlgoTools::update_vertex(&b.edges[0], 0.5, &mut v.0);
        let tol = BRepTool::vertex_tolerance(&v);
        assert!(tol > 0.7, "tolerance {tol} must cover the edge point");

        // A vertex already covering the edge point is left untouched.
        let mut v2 = TopoBuilder::new().make_vertex(GpPnt::new(0.5, 0.0, 0.0), 0.1);
        AlgoTools::update_vertex(&b.edges[0], 0.5, &mut v2.0);
        assert_eq!(BRepTool::vertex_tolerance(&v2), 0.1);
    }

    #[test]
    fn get_normal_to_surface_plane_is_z() {
        let b = TopoBuilder::new();
        let face = b.make_face_plane(&GpPln::new(GpAx3::standard()));
        let s = BRepTool::face_surface(&face).expect("plane surface");
        let n = AlgoTools::get_normal_to_surface(s.as_ref(), 0.5, 0.5).expect("normal");
        // Standard plane: unit normal along ±Z.
        assert!(n.xyz().z.abs() > 0.99, "normal {:?}", n);
        assert!(n.square_magnitude() - 1.0 < 1e-6, "unit normal");
    }

    #[test]
    fn edge_to_face_and_adjust_pcurve_on_surf() {
        let b = unit_box();
        let edge = &b.edges[0];
        let face = &b.faces[0];
        let pc = AlgoTools2D::edge_to_face(edge, face).expect("edge to face");
        // Edge 0 (0,0,0)→(1,0,0) on the bottom face maps to a (0, v) isoline.
        let q0 = pc.d0(0.0);
        assert!((q0.x() - 0.0).abs() < 1e-6 && (q0.y() - 0.0).abs() < 1e-6, "start {:?}", q0);

        // Adjusting keeps the pcurve inside the face UV bounds.
        let adj = AlgoTools2D::adjust_pcurve_on_surf(&pc, face, 1e-7).expect("adjusted");
        let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(face);
        for i in 0..=8 {
            let q = adj.d0(i as f64 / 8.0);
            assert!(q.x() >= umin - 1e-6 && q.x() <= umax + 1e-6, "u {} at i={}", q.x(), i);
            assert!(q.y() >= vmin - 1e-6 && q.y() <= vmax + 1e-6, "v {} at i={}", q.y(), i);
        }
    }

    #[test]
    fn shape_list_dedupes_and_type_count() {
        let b = unit_box();
        let faces: Vec<TopoShape> = b.faces.iter().map(|f| f.0.clone()).collect();
        assert_eq!(BOPToolsSet::type_count(&faces, ShapeType::Face), 6);
        assert_eq!(BOPToolsSet::type_count(&faces, ShapeType::Edge), 0);

        // Duplicated faces collapse to the 6 distinct ones.
        let mut dup = faces.clone();
        dup.extend(faces.iter().cloned());
        assert_eq!(BOPToolsSet::shape_list(&dup).len(), 6);

        // Repeated identical edges collapse to one.
        let e0 = b.edges[0].0.clone();
        let lst = vec![e0.clone(), e0.clone(), b.edges[1].0.clone()];
        assert_eq!(BOPToolsSet::shape_list(&lst).len(), 2);
    }

    // ---- Phase 18b tests ----------------------------------------------------

    #[test]
    fn compute_state_face_in_out_on() {
        let b = unit_box();
        let face = &b.faces[0]; // bottom (z = 0)
        // Interior point of the face → In.
        assert_eq!(
            AlgoTools::compute_state(&face.0, &GpPnt::new(0.5, 0.5, 0.0), 1e-6).unwrap(),
            FaceState::In
        );
        // Point above the face plane → Out.
        assert_eq!(
            AlgoTools::compute_state(&face.0, &GpPnt::new(0.5, 0.5, 0.1), 1e-6).unwrap(),
            FaceState::Out
        );
        // Point on the face boundary edge → On.
        assert_eq!(
            AlgoTools::compute_state(&face.0, &GpPnt::new(0.5, 0.0, 0.0), 1e-6).unwrap(),
            FaceState::On
        );
    }

    #[test]
    fn compute_state_solid_in_out_on() {
        let b = unit_box();
        assert_eq!(
            AlgoTools::compute_state(&b.solid.0, &GpPnt::new(0.5, 0.5, 0.5), 1e-6).unwrap(),
            FaceState::In
        );
        assert_eq!(
            AlgoTools::compute_state(&b.solid.0, &GpPnt::new(2.0, 0.0, 0.0), 1e-6).unwrap(),
            FaceState::Out
        );
        // On the front face (y = 0).
        assert_eq!(
            AlgoTools::compute_state(&b.solid.0, &GpPnt::new(0.5, 0.0, 0.5), 1e-6).unwrap(),
            FaceState::On
        );
    }

    #[test]
    fn compute_state_edge_and_vertex() {
        let b = unit_box();
        // Edge 0 runs (0,0,0)→(1,0,0).
        assert_eq!(
            AlgoTools::compute_state(&b.edges[0].0, &GpPnt::new(0.5, 0.0, 0.0), 1e-6).unwrap(),
            FaceState::On
        );
        assert_eq!(
            AlgoTools::compute_state(&b.edges[0].0, &GpPnt::new(0.5, 1.0, 0.0), 1e-6).unwrap(),
            FaceState::Out
        );
        // Vertex at a corner: coincident → On, far away → Out.
        assert_eq!(
            AlgoTools::compute_state(&b.vertices[0].0, &b.corners[0], 1e-6).unwrap(),
            FaceState::On
        );
        assert_eq!(
            AlgoTools::compute_state(&b.vertices[0].0, &GpPnt::new(5.0, 5.0, 5.0), 1e-6).unwrap(),
            FaceState::Out
        );
    }

    #[test]
    fn make_connexity_block_box_faces_single() {
        let b = unit_box();
        let faces: Vec<TopoShape> = b.faces.iter().map(|f| f.0.clone()).collect();
        let block = AlgoTools::make_connexity_block(&faces);
        assert_eq!(block.shapes().len(), 6);
        assert!(block.is_regular());
    }

    #[test]
    fn make_connexity_blocks_box_one_block() {
        let b = unit_box();
        let faces: Vec<TopoShape> = b.faces.iter().map(|f| f.0.clone()).collect();
        let blocks = AlgoTools::make_connexity_blocks(&faces);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].shapes().len(), 6);
    }

    #[test]
    fn make_connexity_blocks_separates_disjoint_faces() {
        let tb = TopoBuilder::new();
        let sq = |tb: &TopoBuilder, off: f64| -> TopoShape {
            let pts = [
                GpPnt::new(off, 0.0, 0.0),
                GpPnt::new(off + 1.0, 0.0, 0.0),
                GpPnt::new(off + 1.0, 1.0, 0.0),
                GpPnt::new(off, 1.0, 0.0),
            ];
            let edges: Vec<Edge> = (0..4)
                .map(|i| tb.make_edge_segment(&pts[i], &pts[(i + 1) % 4]))
                .collect();
            build_face_with_holes(&edges, &[]).expect("square face").0
        };
        let f1 = sq(&tb, 0.0);
        let f2 = sq(&tb, 10.0);
        let blocks = AlgoTools::make_connexity_blocks(&[f1, f2]);
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].shapes().len(), 1);
        assert_eq!(blocks[1].shapes().len(), 1);
    }

    #[test]
    fn make_blocks_connects_transitive_chain() {
        // 0~1, 1~2 (no direct 0~2): MakeBlocks joins the chain into one block.
        let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
        adj.entry(0).or_default().extend([1]);
        adj.entry(1).or_default().extend([0, 2]);
        adj.entry(2).or_default().extend([1]);
        let blocks = AlgoTools::make_blocks(&adj);
        assert_eq!(blocks.len(), 1, "transitive chain must collapse into one block");
        assert_eq!(blocks[0].len(), 3);
        assert!(blocks[0].contains(&0) && blocks[0].contains(&1) && blocks[0].contains(&2));
    }

    #[test]
    fn make_blocks_separates_disconnected() {
        // Two disconnected pairs and one isolated key.
        let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
        adj.entry(0).or_default().extend([1]);
        adj.entry(1).or_default().extend([0]);
        adj.entry(2).or_default().extend([3]);
        adj.entry(3).or_default().extend([2]);
        adj.entry(4).or_default(); // isolated: present as a key, no neighbours
        let blocks = AlgoTools::make_blocks(&adj);
        assert_eq!(blocks.len(), 3);
        let mut sizes: Vec<usize> = blocks.iter().map(|b| b.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![1, 2, 2]);
    }

    #[test]
    fn orient_edges_on_wire_chains_shuffled_square() {
        let tb = TopoBuilder::new();
        let p = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let e1 = tb.make_edge_segment(&p[0], &p[1]);
        let e2 = tb.make_edge_segment(&p[1], &p[2]);
        let e3 = tb.make_edge_segment(&p[2], &p[3]);
        let e4 = tb.make_edge_segment(&p[3], &p[0]);
        let mut wire = tb.make_wire(&[e2, e4, e1, e3]);
        assert!(!wire_is_closed(&wire));
        AlgoTools::orient_edges_on_wire(&mut wire.0);
        assert!(wire_is_closed(&wire));
    }

    #[test]
    fn orient_faces_on_shell_dedupes_and_orders() {
        let b = unit_box();
        let shell = TopoBuilder::new().make_shell(&b.faces);
        let mut sh = shell.0;
        AlgoTools::orient_faces_on_shell(&mut sh);
        assert_eq!(sh.tshape.read().unwrap().children.len(), 6);
        for h in sh.tshape.read().unwrap().children.iter() {
            assert_eq!(h.shape_type(), ShapeType::Face);
        }
    }

    #[test]
    fn copy_edge_preserves_geometry() {
        let b = unit_box();
        let copy = AlgoTools::copy_edge(&b.edges[0]).expect("copy edge");
        assert!(!copy.0.same_tshape(&b.edges[0].0));
        assert_eq!(BRepTool::edge_parameters(&copy), BRepTool::edge_parameters(&b.edges[0]));
        let (a1, b1) = BRepTool::edge_vertices(&b.edges[0]).expect("orig endpoints");
        let (a2, b2) = BRepTool::edge_vertices(&copy).expect("copy endpoints");
        assert!(a1.is_equal(&a2) && b1.is_equal(&b2));
        assert_eq!(copy.0.tshape.read().unwrap().children.len(), 2);
    }

    #[test]
    fn make_split_edge_bounds_at_params() {
        let b = unit_box();
        let p1 = AlgoTools::point_on_edge(&b.edges[0], 0.25).expect("p1");
        let p2 = AlgoTools::point_on_edge(&b.edges[0], 0.75).expect("p2");
        let v1 = AlgoTools::make_new_vertex(&p1, 1e-7).expect("v1");
        let v2 = AlgoTools::make_new_vertex(&p2, 1e-7).expect("v2");
        let e = AlgoTools::make_split_edge(&b.edges[0], Some(&v1), 0.25, Some(&v2), 0.75)
            .expect("split edge");
        assert_eq!(BRepTool::edge_parameters(&e), (0.25, 0.75));
        let q1 = AlgoTools::point_on_edge(&e, 0.25).expect("start");
        let q2 = AlgoTools::point_on_edge(&e, 0.75).expect("end");
        assert!(q1.distance(&GpPnt::new(0.25, 0.0, 0.0)) < 1e-9, "start {q1:?}");
        assert!(q2.distance(&GpPnt::new(0.75, 0.0, 0.0)) < 1e-9, "end {q2:?}");
    }

    #[test]
    fn is_micro_edge_short_true_long_false() {
        let b = unit_box();
        assert!(!AlgoTools::is_micro_edge(&b.edges[0], 0.5));
        let tb = TopoBuilder::new();
        let micro = tb.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1e-6, 0.0, 0.0));
        assert!(AlgoTools::is_micro_edge(&micro, 0.5));
    }

    #[test]
    fn is_open_shell_closed_false_open_true() {
        let b = unit_box();
        let shell = b.solid.0.tshape.read().unwrap().children[0].clone();
        assert!(!AlgoTools::is_open_shell(&shell));

        // Five faces of the box form an open shell (the right face is missing).
        let open = TopoBuilder::new().make_shell(&b.faces[0..5]);
        assert!(AlgoTools::is_open_shell(&open.0));
    }

    #[test]
    fn is_inverted_solid_normal_box_false() {
        let b = unit_box();
        assert!(!AlgoTools::is_inverted_solid(&b.solid.0));
    }

    #[test]
    fn sense_perpendicular_box_faces_is_zero() {
        let b = unit_box();
        // Bottom (normal −Z) and front (normal −Y) share edge 0; the normals
        // are perpendicular so the sense is 0.
        assert_eq!(AlgoTools::sense(&b.faces[0], &b.faces[2]), 0);
        // Two faces that do not share an edge also report 0.
        assert_eq!(AlgoTools::sense(&b.faces[0], &b.faces[1]), 0);
    }

    #[test]
    fn is_hole_detects_inner_wire() {
        let tb = TopoBuilder::new();
        let outer = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let hole = [
            GpPnt::new(0.75, 0.75, 0.0),
            GpPnt::new(0.75, 0.25, 0.0),
            GpPnt::new(0.25, 0.25, 0.0),
            GpPnt::new(0.25, 0.75, 0.0),
        ];
        let square = |pts: &[GpPnt; 4]| -> Vec<Edge> {
            (0..4).map(|i| tb.make_edge_segment(&pts[i], &pts[(i + 1) % 4])).collect()
        };
        let face = build_face_with_holes(&square(&outer), &[square(&hole)]).expect("face with hole");
        let wires = wires_of_face(&face);
        assert_eq!(wires.len(), 2);
        assert!(!AlgoTools::is_hole(&wires[0].0, &face), "outer wire is not a hole");
        assert!(AlgoTools::is_hole(&wires[1].0, &face), "inner wire is a hole");
    }

    #[test]
    fn dimension_by_shape_type() {
        let b = unit_box();
        assert_eq!(AlgoTools::dimension(&b.vertices[0].0), 0);
        assert_eq!(AlgoTools::dimension(&b.edges[0].0), 1);
        assert_eq!(AlgoTools::dimension(&b.faces[0].0), 2);
        assert_eq!(AlgoTools::dimension(&b.solid.0), 3);
        let wire = TopoBuilder::new().make_wire(&b.edges[0..4]);
        assert_eq!(AlgoTools::dimension(&wire.0), 1);
    }

    #[test]
    fn correct_tolerances_grows_vertex_tolerance() {
        let tb = TopoBuilder::new();
        let v = tb.make_vertex(GpPnt::new(0.0, 0.0, 0.5), 1e-7);
        let lin = GpLin::from_pnt_dir(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let mut e = tb.make_edge(Arc::new(GeomLine::new(lin)), 0.0, 1.0);
        tb.add(&mut e.0, &v.0);
        // The vertex sits 0.5 above the curve start; correct_tolerances must
        // grow its tolerance to cover the endpoint (capped by the 1.0 budget).
        AlgoTools::correct_tolerances(&e.0, 1.0);
        assert!(
            BRepTool::vertex_tolerance(&v) >= 0.5,
            "tol {}",
            BRepTool::vertex_tolerance(&v)
        );
    }

    #[test]
    fn get_edge_off_none_for_edge_on_face() {
        let b = unit_box();
        // Edge 0 lies on the bottom face → no off-face parameter.
        assert!(AlgoTools::get_edge_off(&b.edges[0], &b.faces[0]).is_none());
    }
