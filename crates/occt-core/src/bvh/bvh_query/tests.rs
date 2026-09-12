use super::prelude::*;
use super::*;

    /// A unit cube centred at the origin: 8 vertices, 12 triangles.
    ///
    /// The winding matches `geom::mesh_analysis::box_mesh` (outward CCW) so
    /// the divergence-theorem volume comes out positive.
    fn cube_triangles() -> Vec<(GpPnt, GpPnt, GpPnt)> {
        let v: Vec<GpPnt> = vec![
            GpPnt::new(-0.5, -0.5, -0.5), GpPnt::new(0.5, -0.5, -0.5),
            GpPnt::new(0.5, 0.5, -0.5), GpPnt::new(-0.5, 0.5, -0.5),
            GpPnt::new(-0.5, -0.5, 0.5), GpPnt::new(0.5, -0.5, 0.5),
            GpPnt::new(0.5, 0.5, 0.5), GpPnt::new(-0.5, 0.5, 0.5),
        ];
        let t: Vec<(usize, usize, usize)> = vec![
            (0, 3, 2), (0, 2, 1), // -Z
            (4, 5, 6), (4, 6, 7), // +Z
            (0, 1, 5), (0, 5, 4), // -Y
            (3, 7, 6), (3, 6, 2), // +Y
            (0, 4, 7), (0, 7, 3), // -X
            (1, 2, 6), (1, 6, 5), // +X
        ];
        t.into_iter()
            .map(|(i, j, k)| (v[i], v[j], v[k]))
            .collect()
    }

    /// Cube corner-index mesh (parallel to `cube_triangles`).
    fn cube_index_mesh() -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
        let v: Vec<GpPnt> = vec![
            GpPnt::new(-0.5, -0.5, -0.5), GpPnt::new(0.5, -0.5, -0.5),
            GpPnt::new(0.5, 0.5, -0.5), GpPnt::new(-0.5, 0.5, -0.5),
            GpPnt::new(-0.5, -0.5, 0.5), GpPnt::new(0.5, -0.5, 0.5),
            GpPnt::new(0.5, 0.5, 0.5), GpPnt::new(-0.5, 0.5, 0.5),
        ];
        let t: Vec<(usize, usize, usize)> = vec![
            (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
            (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
        ];
        (v, t)
    }

    #[test]
    fn ray_hits_triangle_nearest() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // Ray from below the cube along +Z: front face at z = -0.5, so t = 4.5.
        // Offset in xy keeps the hit strictly inside a face triangle (off the
        // shared diagonal).
        let hit = ray_cast_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, -5.0), GpVec::new(0.0, 0.0, 1.0), 100.0)
            .expect("ray hits the cube");
        assert!((hit.t - 4.5).abs() < 1e-9, "t = {}", hit.t);
        assert!((hit.point.z() + 0.5).abs() < 1e-9, "z = {}", hit.point.z());
        assert!(hit.triangle_index < tris.len());
    }

    #[test]
    fn ray_misses_mesh() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // Pointing away from the cube (downward, cube is above origin).
        assert!(ray_cast_mesh(&bvh, &tris, GpPnt::new(0.0, 0.0, -5.0), GpVec::new(0.0, 0.0, -1.0), 100.0).is_none());
        // Direction is zero.
        assert!(ray_cast_mesh(&bvh, &tris, GpPnt::new(0.0, 0.0, -5.0), GpVec::zero(), 100.0).is_none());
    }

    #[test]
    fn segment_clamped_hit() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // Short segment stopping before the front face (z = -1) -> miss.
        let short = segment_query_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, -5.0), GpPnt::new(0.2, 0.1, -1.0));
        assert!(short.is_none(), "short segment must not reach the cube");
        // Long segment through the cube -> hit at distance 4.5 from a.
        let long = segment_query_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, -5.0), GpPnt::new(0.2, 0.1, 5.0))
            .expect("segment through the cube hits");
        assert!((long.t - 4.5).abs() < 1e-9, "t = {}", long.t);
    }

    #[test]
    fn box_query_overlap() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // A box fully inside the cube overlaps every face's triangle bbox.
        let mut q = BndBox::from_corners(&GpPnt::new(-0.2, -0.2, -0.2), &GpPnt::new(0.2, 0.2, 0.2));
        let hits = box_query_union(&bvh, &q);
        assert!(!hits.is_empty(), "interior box must overlap triangles");
        q = BndBox::from_corners(&GpPnt::new(10.0, 10.0, 10.0), &GpPnt::new(11.0, 11.0, 11.0));
        assert!(box_query_union(&bvh, &q).is_empty(), "distant box must miss");
    }

    #[test]
    fn point_inside_closed() {
        // Tetrahedron (0,0,0),(1,0,0),(0,1,0),(0,0,1) — 4 triangles, closed.
        let o = GpPnt::new(0.0, 0.0, 0.0);
        let x = GpPnt::new(1.0, 0.0, 0.0);
        let y = GpPnt::new(0.0, 1.0, 0.0);
        let z = GpPnt::new(0.0, 0.0, 1.0);
        let tris = vec![(o, x, y), (o, z, x), (o, y, z), (x, z, y)];
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        assert!(point_inside_box(&bvh, &tris, GpPnt::new(0.1, 0.1, 0.1)), "interior point");
        assert!(!point_inside_box(&bvh, &tris, GpPnt::new(10.0, 10.0, 10.0)), "exterior point");
    }

    #[test]
    fn closest_point_on_mesh() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        // Point above the +Z face (z = 0.5): distance should be 1.5.
        let (pt, d) = closest_point_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, 2.0));
        assert!((d - 1.5).abs() < 1e-9, "distance = {d}");
        assert!((pt.z() - 0.5).abs() < 1e-9, "closest z = {}", pt.z());
        // Point exactly on a face: distance ~ 0.
        let (_, d0) = closest_point_mesh(&bvh, &tris, GpPnt::new(0.2, 0.1, 0.5));
        assert!(d0 < 1e-9, "on-surface distance = {d0}");
    }

    #[test]
    fn mesh_volume_box() {
        let tris = cube_triangles();
        let bvh = crate::bvh::builder_tri::build_tri_bvh(&tris, 2);
        let vol = mesh_volume(&bvh, &tris);
        assert!((vol - 1.0).abs() < 0.05, "unit cube volume = {vol}");
    }

    #[test]
    fn surface_area_unit_square() {
        // Two triangles on a unit square in the z=0 plane.
        let tris = vec![
            (GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0)),
            (GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)),
        ];
        assert!((mesh_surface_area(&tris) - 1.0).abs() < 1e-12, "area = {}", mesh_surface_area(&tris));
    }

    #[test]
    fn edge_topology_boundary_and_manifold() {
        // Single triangle: 3 boundary edges, 0 non-manifold.
        assert_eq!(mesh_edge_topology(&[(0, 1, 2)]), (3, 0));
        // Two triangles sharing edge (0,2): 4 boundary edges, 1 shared edge.
        assert_eq!(mesh_edge_topology(&[(0, 1, 2), (0, 2, 3)]), (4, 0));
        // Closed cube: 0 boundary.
        let (_, t) = cube_index_mesh();
        assert_eq!(mesh_edge_topology(&t), (0, 0));
    }

    #[test]
    fn manifold_closed_cube() {
        let (_, t) = cube_index_mesh();
        assert!(mesh_is_manifold_closed(&t), "cube is closed and manifold");
        // A single triangle is not closed (boundary edges).
        assert!(!mesh_is_manifold_closed(&[(0, 1, 2)]));
    }

    #[test]
    fn connected_components_two_boxes() {
        // Two disjoint cubes: vertices 0..7 and 8..15.
        let v: Vec<GpPnt> = (0..16)
            .map(|i| {
                if i < 8 {
                    GpPnt::new(i as f64, 0.0, 0.0)
                } else {
                    GpPnt::new(100.0 + i as f64, 0.0, 0.0)
                }
            })
            .collect();
        let cube: Vec<(usize, usize, usize)> = vec![
            (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
            (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
        ];
        let shifted: Vec<(usize, usize, usize)> =
            cube.iter().map(|&(a, b, c)| (a + 8, b + 8, c + 8)).collect();
        let mut tris = cube;
        tris.extend(shifted);
        let comps = mesh_connected_components(&tris);
        assert_eq!(comps.len(), 2, "two disjoint cubes -> two components");
        // Sorted by smallest vertex: first component is 0..8.
        assert_eq!(comps[0].len(), 8);
        assert_eq!(comps[1].len(), 8);
        assert!(comps.iter().all(|c| !c.is_empty()));
        let _ = v;
    }

    #[test]
    fn boundary_loops_square() {
        // Square: (0,0)-(1,0)-(1,1)-(0,1), two triangles.
        let tris = vec![(0, 1, 2), (0, 2, 3)];
        let loops = mesh_boundary_loops(&tris);
        assert_eq!(loops.len(), 1, "a square has one boundary loop");
        // Closed loop: first == last, four unique vertices -> 5 entries.
        let l = &loops[0];
        assert_eq!(l.len(), 5, "loop = {l:?}");
        assert_eq!(l[0], l[4], "loop is closed");
        // Closed cube: no boundary loops at all.
        let (_, t) = cube_index_mesh();
        assert!(mesh_boundary_loops(&t).is_empty());
    }

    #[test]
    fn ear_clipping_square() {
        let sq = [
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(0.0, 1.0),
        ];
        let tris = polygon_triangulate_ear(&sq).expect("square triangulates");
        assert_eq!(tris.len(), 2);
        let mut area = 0.0;
        for (a, b, c) in &tris {
            let (pa, pb, pc) = (sq[*a], sq[*b], sq[*c]);
            area += 0.5
                * ((pb.x() - pa.x()) * (pc.y() - pa.y()) - (pb.y() - pa.y()) * (pc.x() - pa.x()))
                    .abs();
        }
        assert!((area - 1.0).abs() < 1e-12, "triangulated area = {area}");
    }

    #[test]
    fn incident_triangle_lists() {
        let (v, t) = cube_index_mesh();
        let incident = mesh_to_polygon_indices(&v, &t);
        assert_eq!(incident.len(), 8);
        for (i, list) in incident.iter().enumerate() {
            assert!(!list.is_empty(), "vertex {i} has no incident triangles");
            // Every triangle in the list references vertex i.
            for &ti in list {
                let (a, b, c) = t[ti];
                assert!(a == i || b == i || c == i);
            }
        }
    }
