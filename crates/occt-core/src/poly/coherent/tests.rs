use super::prelude::*;
use super::*;

    fn unit_cube_triangulation() -> PolyTriangulation {
        // 8 vertices of the unit cube, 12 triangles (2 per face).
        let nodes = vec![
            GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.), GpPnt::new(0., 1., 0.),
            GpPnt::new(0., 0., 1.), GpPnt::new(1., 0., 1.), GpPnt::new(1., 1., 1.), GpPnt::new(0., 1., 1.),
        ];
        let tris = vec![
            // bottom
            super::super::triangulation::Triangle::new(0, 1, 2),
            super::super::triangulation::Triangle::new(0, 2, 3),
            // top
            super::super::triangulation::Triangle::new(4, 6, 5),
            super::super::triangulation::Triangle::new(4, 7, 6),
            // front
            super::super::triangulation::Triangle::new(0, 5, 1),
            super::super::triangulation::Triangle::new(0, 4, 5),
            // back
            super::super::triangulation::Triangle::new(3, 2, 6),
            super::super::triangulation::Triangle::new(3, 6, 7),
            // left
            super::super::triangulation::Triangle::new(0, 3, 7),
            super::super::triangulation::Triangle::new(0, 7, 4),
            // right
            super::super::triangulation::Triangle::new(1, 5, 6),
            super::super::triangulation::Triangle::new(1, 6, 2),
        ];
        PolyTriangulation::from_parts(nodes, tris)
    }

    #[test]
    fn build_cube_connectivity() {
        let tri = unit_cube_triangulation();
        let ct = PolyCoherentTriangulation::from_poly_triangulation(&tri).unwrap();
        assert_eq!(ct.n_nodes(), 8);
        assert_eq!(ct.n_triangles(), 12);
        // Per-node incident-triangle counts for this cube triangulation.
        let expected = [6usize, 4, 4, 4, 4, 4, 6, 4];
        for i in 0..8 {
            assert_eq!(ct.nodes[i].triangles.len(), expected[i], "node {i}");
        }
        // A closed cube has no boundary edges.
        assert!(ct.boundary_loops().is_empty());
    }

    #[test]
    fn boundary_loops_open_square() {
        // A flat square split into two triangles: one boundary loop of 4 nodes.
        let mut ct = PolyCoherentTriangulation::new();
        ct.set_node(&GpXyz::new(0., 0., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(1., 0., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(1., 1., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(0., 1., 0.), -1).unwrap();
        ct.add_triangle(0, 1, 2).unwrap();
        ct.add_triangle(0, 2, 3).unwrap();
        let loops = ct.boundary_loops();
        assert_eq!(loops.len(), 1);
        let l = &loops[0];
        // 4 distinct nodes plus the repeated closure node (OCCT convention).
        assert_eq!(l.len(), 5);
        assert_eq!(l[0], l[l.len() - 1]);
        let mut distinct = l.clone();
        distinct.pop();
        distinct.sort_unstable();
        assert_eq!(distinct, vec![0, 1, 2, 3]);
    }

    #[test]
    fn triangle_adjacency_queries() {
        let tri = unit_cube_triangulation();
        let mut ct = PolyCoherentTriangulation::from_poly_triangulation(&tri).unwrap();
        ct.compute_links().unwrap();
        assert!(ct.n_links() > 0);
        // Every triangle in the cube has 3 neighbours.
        for t in 0..ct.triangles.len() {
            let tr = &ct.triangles[t];
            assert_eq!(tr.n_connections(), 3, "triangle {t} has 3 neighbours");
            let neigh = tr.neighbours();
            assert!(neigh.iter().all(|n| n.is_some()));
        }
        // find_triangle on the interior cube edge (0,1): two triangles share it.
        let link = PolyCoherentLink::new_nodes(0, 1);
        let [left, right] = ct.find_triangle(&link);
        assert!(left.is_some());
        assert!(right.is_some());
        assert_ne!(left, right);
    }

    #[test]
    fn find_triangle_boundary_edge_one_side() {
        // Open square: the edge (0,1) belongs to only one triangle.
        let mut ct = PolyCoherentTriangulation::new();
        ct.set_node(&GpXyz::new(0., 0., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(1., 0., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(1., 1., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(0., 1., 0.), -1).unwrap();
        ct.add_triangle(0, 1, 2).unwrap();
        ct.add_triangle(0, 2, 3).unwrap();
        let link = PolyCoherentLink::new_nodes(0, 1);
        let [left, right] = ct.find_triangle(&link);
        assert_eq!(left.and(right).is_some(), false); // exactly one triangle
        assert!(left.is_some() || right.is_some());
    }

    #[test]
    fn link_construction_orders_nodes() {
        // Triangle (2, 0, 1): side 2 is opposite node index 2 (node value 1),
        // the shared edge is between node values 2 and 0.
        let tri = PolyCoherentTriangle::new3(2, 0, 1);
        let link = PolyCoherentLink::from_triangle(&tri, 2);
        // Edge endpoints ordered ascending: 0 then 2.
        assert_eq!(link.node[0], 0);
        assert_eq!(link.node[1], 2);
        // Opposite node on the left (this triangle) is node index 2 (value 1).
        assert_eq!(link.opposite_node[1], 1);
        assert_eq!(link.opposite_node[0], -1); // no neighbour yet
    }

    #[test]
    fn node_free_and_removal() {
        let mut ct = PolyCoherentTriangulation::new();
        ct.set_node(&GpXyz::new(0., 0., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(1., 0., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(0., 1., 0.), -1).unwrap();
        assert!(ct.nodes[0].is_free_node());
        let t0 = ct.add_triangle(0, 1, 2).unwrap();
        assert!(!ct.nodes[0].is_free_node());
        assert_eq!(ct.n_triangles(), 1);
        ct.remove_triangle(t0).unwrap();
        assert!(ct.nodes[0].is_free_node());
        assert_eq!(ct.n_triangles(), 0);
        assert_eq!(ct.get_free_nodes(), vec![0, 1, 2]);
    }

    #[test]
    fn remove_degenerated_sliver() {
        // A nearly-degenerate triangle 0-2-1 (edge 1-2 is tiny) plus a healthy one.
        let mut ct = PolyCoherentTriangulation::new();
        ct.set_node(&GpXyz::new(0., 0., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(1., 0., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(1.0000001, 0., 0.), -1).unwrap();
        ct.set_node(&GpXyz::new(0., 1., 0.), -1).unwrap();
        let _t0 = ct.add_triangle(0, 2, 1).unwrap(); // sliver
        let _t1 = ct.add_triangle(0, 3, 2).unwrap();
        let pairs = ct.remove_degenerated(1e-3).unwrap();
        assert_eq!(pairs.len(), 1);
        // Sliver removed; one triangle remains.
        assert_eq!(ct.n_triangles(), 1);
    }

    #[test]
    fn export_roundtrip_compacts_free_nodes() {
        let tri = unit_cube_triangulation();
        let mut ct = PolyCoherentTriangulation::from_poly_triangulation(&tri).unwrap();
        // Add an extra free node.
        ct.set_node(&GpXyz::new(9., 9., 9.), -1).unwrap();
        let out = ct.get_poly_triangulation().unwrap();
        assert_eq!(out.nb_nodes(), 8);
        assert_eq!(out.nb_triangles(), 12);
    }
