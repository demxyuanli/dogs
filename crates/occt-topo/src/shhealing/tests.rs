use super::prelude::*;
use super::*;
    use crate::brep_builder_api::{make_face_from_polygon, make_wire_from_points};
    use crate::brep_tool::BRepTool;
    use crate::primitives::BRepPrimBox;
    use crate::topo_tools_full::{faces_of, wire_is_closed};

    #[test]
    fn weld_two_boxes_glue_vertices() {
        // Two unit boxes sharing the x = 1 wall.
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = BRepPrimBox::make_box_corner(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 1.0, 1.0));
        let comp = TopoBuilder::new().make_compound_of(&[a.solid.0, b.solid.0]);
        assert_eq!(vertices_of(&comp.0).len(), 16);

        let (healed, welded) = weld_coincident_vertices(&comp.0, 1e-6);
        assert_eq!(welded, 4, "the four shared-wall corner pairs merge");
        assert_eq!(vertices_of(&healed).len(), 12);
    }

    #[test]
    fn weld_within_tolerance() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 0.0);
        let v2 = b.make_vertex(GpPnt::new(1e-4, 0.0, 0.0), 0.0);
        let comp = b.make_compound_of(&[v1.0.clone(), v2.0.clone()]);

        let (merged, n) = weld_coincident_vertices(&comp.0, 1e-3);
        assert_eq!(n, 1, "1e-4 apart merges at tol 1e-3");
        assert_eq!(vertices_of(&merged).len(), 1);

        let (_, n2) = weld_coincident_vertices(&comp.0, 1e-6);
        assert_eq!(n2, 0, "1e-4 apart does not merge at tol 1e-6");
    }

    #[test]
    fn remove_small_edge_collapses() {
        // A tiny first edge followed by a long one. `ShapeAnalysis_Wire::
        // CheckSmall` flags the tiny edge (its endpoints and its curve midpoint
        // are all within `precsmall`), so `ShapeFix_Wire::FixSmall` drops it.
        // The `FixConnected(n, precsmall)` follow-up only merges the
        // neighbouring vertices when *they* are within `precsmall` - here they
        // are the ends of the one surviving edge, ~1.0 apart, so nothing is
        // merged and the survivor keeps its own endpoints
        // (`ShapeFix_Wire.cxx:1459-1468`, `ShapeAnalysis_Wire.cxx:693-762`).
        let w = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.0001, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ])
        .unwrap();
        assert_eq!(edges_of_wire(&w).len(), 2);

        let (healed, n) = remove_small_edges(&w.0, 0.01);
        assert_eq!(n, 1);
        let healed_wire = Wire(healed);
        let he = edges_of_wire(&healed_wire);
        assert_eq!(he.len(), 1, "the tiny edge is removed from the wire");
        let (a, z) = edge_vertices(&he[0]);
        assert!(
            vertex_position(&a.unwrap()).distance(&GpPnt::new(0.0001, 0.0, 0.0)) < 1e-9,
            "the survivor keeps the vertex it already had"
        );
        assert!(vertex_position(&z.unwrap()).distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9);
    }

    #[test]
    fn close_open_wire_adds_edge() {
        // A 3-edge wire missing the closing edge; endpoints ~1.4e-6 apart.
        let w = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1e-6, 1e-6, 0.0),
        ])
        .unwrap();

        let (healed, n) = close_open_wires(&w.0, 1e-3);
        assert_eq!(n, 1);
        let healed_wire = Wire(healed);
        assert_eq!(edges_of_wire(&healed_wire).len(), 4);
        assert!(wire_is_closed(&healed_wire), "closing edge completes the loop");
    }

    #[test]
    fn heal_shape_full_fixes_all() {
        // A wire with a tiny edge (2e-3: longer than the weld tol so it is the
        // small-edge fix that removes it), plus a near-closed wire whose ends
        // are 1.4e-6 apart (within the weld tol).
        let w1 = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.002, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ])
        .unwrap();
        let w2 = make_wire_from_points(&[
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1.0, 2.0, 0.0),
            GpPnt::new(1e-6, 1.0 + 1e-6, 0.0),
        ])
        .unwrap();
        let comp = TopoBuilder::new().make_compound_of(&[w1.0, w2.0]);

        let (healed, report) = heal_shape(&comp.0, 1e-3, 0.01);
        assert!(report.modified);
        assert!(report.small_edges_removed > 0, "tiny edge removed by FixSmall");
        assert!(report.vertices_welded > 0, "near-closed wire's ends welded");
        assert!(
            edges_of(&healed).iter().all(|e| edge_length(e, 8) >= 0.01 - 1e-9),
            "no edge shorter than the minimum remains"
        );
        // The near-closed wire is now a closed loop inside the healed compound.
        assert!(
            wires_of(&healed).iter().any(|w| wire_is_closed(w)),
            "the near-closed wire closes after healing"
        );
    }

    #[test]
    fn free_edges_shared_vs_free() {
        // A box's wires are all closed → no free edges.
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(free_edges(&boxy.solid.0, 1e-6).is_empty());

        // An open 2-edge wire: both edges carry a free endpoint.
        let w = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
        ])
        .unwrap();
        let free = free_edges(&w.0, 1e-6);
        assert_eq!(free.len(), 2, "both edges of the open wire are free");
    }

    #[test]
    fn wire_closed_after_heal() {
        let w = make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1e-6, 1e-6, 0.0),
        ])
        .unwrap();
        assert!(!wire_is_closed(&w), "the wire starts out open");
        assert!(wire_is_closed_after_heal(&w.0, 1e-3));
    }

    #[test]
    fn move_vertex_relocates() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let verts = vertices_of(&b.solid.0);
        let origin = verts
            .iter()
            .find(|v| vertex_position(v).is_equal(&GpPnt::zero()))
            .expect("corner at the origin")
            .clone();

        let healed = move_vertex(&b.solid.0, &origin, &GpPnt::new(0.1, 0.0, 0.0)).unwrap();
        let hv = vertices_of(&healed);
        let moved = hv
            .iter()
            .find(|v| vertex_position(v).distance(&GpPnt::new(0.1, 0.0, 0.0)) < 1e-9)
            .expect("moved vertex present at the new point");

        // An incident edge's curve now starts near the moved point.
        let edges = edges_of(&healed);
        let incident = edges
            .iter()
            .find(|e| {
                let (a, z) = edge_vertices(e);
                let (Some(a), Some(z)) = (a, z) else { return false };
                is_same(&a.0, &moved.0) || is_same(&z.0, &moved.0)
            })
            .expect("incident edge");
        let curve = BRepTool::edge_curve(incident).expect("edge curve");
        let (f0, f1) = BRepTool::edge_parameters(incident);
        let touches_moved = curve.d0(f0).distance(&GpPnt::new(0.1, 0.0, 0.0)) < 1e-9
            || curve.d0(f1).distance(&GpPnt::new(0.1, 0.0, 0.0)) < 1e-9;
        assert!(touches_moved, "an incident edge's curve passes through the moved point");
    }

    #[test]
    fn heal_identity_on_good_shape() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (healed, report) = heal_shape(&b.solid.0, 1e-3, 0.01);
        assert!(!report.modified, "a valid box needs no healing");
        assert_eq!(report.small_edges_removed, 0);
        assert_eq!(report.wires_closed, 0);
        assert_eq!(report.vertices_welded, 0);
        assert_eq!(report.free_edges_fixed, 0);
        assert_eq!(vertices_of(&healed).len(), 8);
        assert_eq!(edges_of(&healed).len(), 12);
    }

    #[test]
    fn weld_merge_deduplicates_vertices() {
        // Two vertex instances at the exact same point merge into one.
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let v2 = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let comp = b.make_compound_of(&[v1.0, v2.0]);

        let (healed, n) = weld_coincident_vertices(&comp.0, 1e-9);
        assert_eq!(n, 1);
        assert_eq!(vertices_of(&healed).len(), 1);
    }

    #[test]
    fn remove_small_isolated_loop() {
        // A tiny standalone triangle loop. Every edge is small, but
        // `ShapeAnalysis_Wire::CheckSmall` bails out once the wire holds a
        // single edge (`ShapeAnalysis_Wire.cxx:768-771`), so `FixSmall` removes
        // two of the three and the loop survives with one edge. A face whose
        // wires all vanish is what `ShapeFix_Face::FixSmallAreaWire`
        // (`ShapeFix_Face.cxx:2331-2394`) drops, not this pass.
        let face = make_face_from_polygon(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1e-4, 0.0, 0.0),
            GpPnt::new(0.0, 1e-4, 0.0),
        ])
        .unwrap();
        let (healed, n) = remove_small_edges(&face.0, 0.01);
        assert_eq!(n, 2, "CheckSmall stops at one remaining edge");
        assert_eq!(edges_of(&healed).len(), 1);
        assert_eq!(wires_of(&healed).len(), 1, "the wire is not dropped");
    }

    #[test]
    fn heal_preserves_solid_closure() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (healed, _) = heal_shape(&b.solid.0, 1e-6, 0.01);
        let shells = crate::topo_tools_full::shapes_of(&healed, ShapeType::Shell);
        assert_eq!(shells.len(), 1);
        let shell = Shell(shells[0].clone());
        assert!(crate::shell_check::shell_is_closed(&shell), "shell stays closed");
        let (v, e, f) = (
            vertices_of(&healed).len(),
            edges_of(&healed).len(),
            faces_of(&healed).len(),
        );
        assert_eq!((v, e, f), (8, 12, 6));
    }
