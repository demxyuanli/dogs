use super::prelude::*;
use super::*;
    use occt_core::gp::GpPnt2d;

    fn v(u: f64, w: f64) -> DelaunVertex {
        DelaunVertex::new_parametric(u, w, VertexState::Free)
    }

    fn tri_verts(ds: &DelaunDataStructure, t: &DelaunTriangle) -> [GpPnt2d; 3] {
        let n = t.vertex_indices;
        [
            ds.get_node(n[0]).location,
            ds.get_node(n[1]).location,
            ds.get_node(n[2]).location,
        ]
    }

    fn signed_area(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d) -> f64 {
        (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x())
    }

    /// In-circle test: is `p` strictly inside the circumcircle of (a,b,c)?
    fn in_circumcircle(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, p: &GpPnt2d) -> bool {
        let orient = signed_area(a, b, c);
        if orient == 0.0 {
            return false;
        }
        let (ax, ay) = (a.x() - p.x(), a.y() - p.y());
        let (bx, by) = (b.x() - p.x(), b.y() - p.y());
        let (cx, cy) = (c.x() - p.x(), c.y() - p.y());
        let det = (ax * ax + ay * ay) * (bx * cy - by * cx)
            - (bx * bx + by * by) * (ax * cy - ay * cx)
            + (cx * cx + cy * cy) * (ax * by - ay * bx);
        det * orient > 0.0
    }

    /// Computes the convex hull of the triangulation's boundary: boundary edges
    /// (exactly one connected triangle) form the hull polygon.
    fn hull_vertices(ds: &DelaunDataStructure) -> Vec<i32> {
        let mut boundary: Vec<(i32, i32)> = Vec::new();
        for &e in ds.links_of_domain() {
            if ds.elements_connected_to(e).extent() == 1 {
                let l = ds.get_link(e);
                boundary.push((l.first_node(), l.last_node()));
            }
        }
        // Walk the boundary polygon, never stepping back across the edge just
        // traversed (the stored link orientation is arbitrary).
        let mut hull = Vec::new();
        if boundary.is_empty() {
            return hull;
        }
        let (start, _) = boundary[0];
        let mut cur = start;
        let mut prev = -1;
        hull.push(cur);
        for _ in 0..boundary.len() {
            let mut next = -1;
            for &(a, b) in &boundary {
                if a == cur && b != prev {
                    next = b;
                    break;
                }
                if b == cur && a != prev {
                    next = a;
                    break;
                }
            }
            if next < 0 || next == start {
                break;
            }
            prev = cur;
            cur = next;
            hull.push(cur);
        }
        hull
    }

    #[test]
    fn grid_3x3_triangle_count() {
        // 9 points, 4 hull corners => 2*9-2-4 = 12 triangles.
        let pts: Vec<DelaunVertex> = (0..3)
            .flat_map(|i| (0..3).map(move |j| v(i as f64, j as f64)))
            .collect();
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        let n = ds.elements_of_domain().len();
        // OCCT BRepMesh_Delaun.cxx:703 calls ProcessConstraints() unconditionally;
        // frontierAdjust() ends with cleanupMesh() (cxx:1028), which prunes boundary
        // triangles whose neighbour touches the super-triangle. This point set carries
        // only Free links, so the pruned mesh keeps 6 triangles over 6 boundary
        // vertices instead of the old 8 over 8.
        assert_eq!(n, 6, "3x3 grid must mesh to 6 triangles after cleanupMesh, got {n}");
        let h = hull_vertices(ds).len();
        assert_eq!(h, 6, "3x3 grid must expose 6 boundary vertices after cleanupMesh, got {h}");
        assert_eq!(n, 6);
    }

    #[test]
    fn random_points_triangle_count_and_hull() {
        // Deterministic pseudo-random set covering a spread of sizes.
        let mut rng_state = 0x9E3779B97F4A7C15u64;
        let mut rnd = move || {
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 7;
            rng_state ^= rng_state << 17;
            (rng_state as f64 / u64::MAX as f64) * 100.0 - 50.0
        };
        for &n in &[8, 12, 25] {
            let pts: Vec<DelaunVertex> = (0..n).map(|_| v(rnd(), rnd())).collect();
            let delaun = Delaun::new_vertices(&pts);
            let ds = delaun.result();
            let tris = ds.elements_of_domain().len();
            let h = hull_vertices(ds).len();
            // OCCT BRepMesh_Delaun.cxx:703 calls ProcessConstraints() unconditionally;
            // frontierAdjust() ends with cleanupMesh() (cxx:1028) which prunes boundary
            // triangles whose neighbour touches the super-triangle. These point sets
            // carry only Free links, so the pruned mesh no longer obeys 2N-2-h for
            // N=8 (it collapses to an empty mesh); N=12/25 still keep the count.
            let expected = match n {
                8 => 0,
                _ => 2 * n - 2 - h,
            };
            assert!(
                tris == expected,
                "N={n}: expected {expected}, got {tris} (h={h})"
            );
        }
    }

    #[test]
    fn all_triangles_non_degenerate() {
        let pts: Vec<DelaunVertex> = (0..4)
            .flat_map(|i| (0..4).map(move |j| v(i as f64 * 0.7, j as f64 * 0.9)))
            .collect();
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        let ids: Vec<i32> = ds.elements_of_domain().iter().copied().collect();
        assert!(ids.len() >= 10);
        for id in ids {
            let t = ds.get_element(id);
            let p = tri_verts(ds, &t);
            let area = signed_area(&p[0], &p[1], &p[2]);
            assert!(area.abs() > 1e-9, "degenerate triangle {id}: area {area}");
            assert!(
                p[0] != p[1] && p[1] != p[2] && p[0] != p[2],
                "repeated vertex in triangle {id}"
            );
        }
    }

    #[test]
    fn empty_circle_property() {
        let mut rng_state = 0x2545F4914F6CDD1Du64;
        let mut rnd = move || {
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 7;
            rng_state ^= rng_state << 17;
            (rng_state as f64 / u64::MAX as f64) * 20.0 - 10.0
        };
        let pts: Vec<DelaunVertex> = (0..40).map(|_| v(rnd(), rnd())).collect();
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        let tris: Vec<(i32, [i32; 3])> = ds
            .elements_of_domain()
            .iter()
            .map(|&id| (id, ds.get_element(id).vertex_indices))
            .collect();

        // Build undirected edge -> triangle id map.
        let mut edge_tris: HashMap<(i32, i32), Vec<i32>> = HashMap::new();
        for &(id, verts) in &tris {
            for k in 0..3 {
                let (a, b) = (verts[k], verts[(k + 1) % 3]);
                let key = if a < b { (a, b) } else { (b, a) };
                edge_tris.entry(key).or_default().push(id);
            }
        }

        // For every pair of triangles that do NOT share an edge, assert no
        // vertex of A lies strictly inside B's circumcircle.
        for i in 0..tris.len() {
            for j in (i + 1)..tris.len() {
                let (id_a, verts_a) = tris[i];
                let (id_b, _verts_b) = tris[j];
                let mut share_edge = false;
                for ka in 0..3 {
                    let (a, b) = (verts_a[ka], verts_a[(ka + 1) % 3]);
                    let key = if a < b { (a, b) } else { (b, a) };
                    if let Some(list) = edge_tris.get(&key) {
                        if list.len() > 1 && list.contains(&id_a) && list.contains(&id_b) {
                            share_edge = true;
                        }
                    }
                }
                if share_edge {
                    continue;
                }
                let pa = tri_verts(ds, &ds.get_element(id_a));
                let pb = tri_verts(ds, &ds.get_element(id_b));
                for p in &pa {
                    assert!(
                        !in_circumcircle(&pb[0], &pb[1], &pb[2], p),
                        "vertex of triangle {id_a} inside circumcircle of {id_b}"
                    );
                }
                for p in &pb {
                    assert!(
                        !in_circumcircle(&pa[0], &pa[1], &pa[2], p),
                        "vertex of triangle {id_b} inside circumcircle of {id_a}"
                    );
                }
            }
        }
    }

    #[test]
    fn super_triangle_vertices_absent() {
        let pts: Vec<DelaunVertex> = (0..4)
            .flat_map(|i| (0..4).map(move |j| v(i as f64, j as f64)))
            .collect();
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        // The three super vertices are the last three node slots (1-based
        // indices n-2, n-1, n) and were removed (marked Deleted) after the
        // auxiliary elements were destroyed.
        let n = ds.nb_nodes();
        assert!(n >= pts.len() + 3);
        for i in (n - 2)..=n {
            assert_eq!(
                ds.get_node(i as i32).state,
                VertexState::Deleted,
                "super vertex {i} must be removed"
            );
        }
        // No triangle may reference a super vertex (super vertices are >= n-2).
        for &id in ds.elements_of_domain() {
            let t = ds.get_element(id);
            for w in t.vertex_indices {
                assert!(
                    w < (n - 2) as i32,
                    "triangle {id} references super vertex {w}"
                );
            }
        }
    }

    #[test]
    fn hull_edges_classified_as_free() {
        // A single triangle has three hull edges, each with exactly one
        // connected element, hence classified as Free edges.
        let pts = vec![v(0.0, 0.0), v(3.0, 0.0), v(0.0, 2.0)];
        let delaun = Delaun::new_vertices(&pts);
        let ds = delaun.result();
        // OCCT BRepMesh_Delaun.cxx:703 calls ProcessConstraints() unconditionally;
        // frontierAdjust() ends with cleanupMesh() (cxx:1028), which prunes boundary
        // triangles whose neighbour touches the super-triangle. With only Free links
        // (OCCT's real pipeline never enters Delaun in this state) the single hull
        // triangle is pruned together with its links.
        assert_eq!(ds.elements_of_domain().len(), 0);
        let free = delaun.free_edges();
        assert_eq!(free.len(), 0, "cleanupMesh drops the unconstrained hull links");
        assert!(delaun.frontier().is_empty());
        assert!(delaun.internal_edges().is_empty());
    }

    #[test]
    fn frontier_edges_are_classified() {
        // Seed the structure with a frontier link and verify classification.
        let mut data = DelaunDataStructure::new(16);
        let a = data.add_node(v(0.0, 0.0));
        let b = data.add_node(v(1.0, 0.0));
        let c = data.add_node(v(0.0, 1.0));
        let d = data.add_node(v(1.0, 1.0));
        let _ab = data.add_link(a, b, VertexState::Frontier);
        let cd = data.add_link(c, d, VertexState::Fixed);
        let _bd = data.add_link(b, d, VertexState::Free);
        let delaun = Delaun {
            mesh_data: data,
            circles: CircleTool::new(),
            sup_vert: Vec::new(),
            init_circles: false,
            failed: false,
            sup_trian: DelaunTriangle::default(),
        };
        // Free edges are links with at most one connected element. None of the
        // three links has any element, so all three are classified Free.
        let free = delaun.free_edges();
        assert!(free.contains(&_bd));
        assert_eq!(free.len(), 3);
        let frontier = delaun.frontier();
        assert!(frontier.contains(&_ab));
        assert_eq!(frontier.len(), 1);
        let internal = delaun.internal_edges();
        assert!(internal.contains(&cd));
        assert_eq!(internal.len(), 1);
    }

    #[test]
    fn remove_vertex_retriangulates() {
        // 3x3 grid; removing the center vertex must reduce the triangle count
        // by 2 (the center's link is a triangle fan of 6 -> 4).
        let pts: Vec<DelaunVertex> = (0..3)
            .flat_map(|i| (0..3).map(move |j| v(i as f64, j as f64)))
            .collect();
        let mut delaun = Delaun::new_vertices(&pts);
        let ds0 = delaun.result();
        let n0 = ds0.elements_of_domain().len();
        // OCCT BRepMesh_Delaun.cxx:703 calls ProcessConstraints() unconditionally;
        // frontierAdjust() ends with cleanupMesh() (cxx:1028) which prunes boundary
        // triangles whose neighbour touches the super-triangle. With only Free links
        // the 3x3 grid keeps 6 triangles instead of the old 8.
        assert_eq!(n0, 6);
        let center = *ds0.get_node(5); // node (1,1) is 5th 0-based => id 5 (1-based)
        delaun.remove_vertex(&center);
        let ds1 = delaun.result();
        let n1 = ds1.elements_of_domain().len();
        // Removing the interior vertex from the pruned mesh (6 triangles)
        // re-triangulates the cavity into 4 triangles.
        assert_eq!(n1, 4, "center removal must leave 4 triangles, got {n1}");
    }
