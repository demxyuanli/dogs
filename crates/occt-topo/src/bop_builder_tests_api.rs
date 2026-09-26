use super::*;
use super::tests::{box_vol, clear_tree, disjoint_box_shapes, overlapping_boxes, test_box_at};
use crate::builder::TopoBuilder;
use crate::primitives::BRepPrimBox;
use crate::shape::{Shell};

use occt_core::gp::GpPnt;

    #[test]
    fn boolean_multi_single_shape_identity() {
        let shapes = disjoint_box_shapes();
        let r = boolean_multi(&shapes[..1], BoolOp::Fuse, 1e-6).expect("single fuse ok");
        assert_eq!(r.shape.shape_type(), ShapeType::Solid, "one shape → itself");
        assert!(r.solid.is_some());
        assert!((box_vol(&r.shape) - 1.0).abs() < 0.05);
        clear_tree(&r.shape);
        clear_tree(&shapes[0]);
    }

    #[test]
    fn boolean_multi_overlap_folds() {
        // Three boxes overlapping along x (each overlaps the next).
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0));
        let b3 = test_box_at(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 1.0, 1.0));
        let shapes = vec![b1.0.clone(), b2.0.clone(), b3.0.clone()];
        let r = boolean_multi(&shapes, BoolOp::Fuse, 1e-6).expect("multi fuse ok");
        assert!(r.solid.is_some(), "overlapping fuse produces a solid");
        assert!(shell_is_closed(&r.shells[0]), "overlapping fuse shell closed");
        let v = box_vol(&r.shape);
        // Analytic union: b1∪b2 = 1.5; ∪b3 = [0,2]×[0,1]×[0,1] = 2.0.
        assert!((v - 2.0).abs() < 0.2, "multi fuse volume {v} (expected ~2.0)");
        clear_tree(&r.shape);
        for s in &shapes {
            clear_tree(s);
        }
    }

    #[test]
    fn boolean_compound_cut() {
        let big = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 3.0, 3.0));
        let s1 = test_box_at(&GpPnt::new(0.5, 0.5, 0.5), &GpPnt::new(1.5, 1.5, 1.5));
        let s2 = test_box_at(&GpPnt::new(1.7, 0.5, 0.5), &GpPnt::new(2.7, 1.5, 1.5));
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[s1.0.clone(), s2.0.clone()]);
        let r = boolean_compound(&big.0, &comp.0, BoolOp::Cut, 1e-6).expect("compound cut ok");
        let v = box_vol(&r.shape);
        // 3³ − 1 − 1 = 25.
        assert!((v - 25.0).abs() < 0.15, "compound cut volume {v} (expected ~25.0)");
        assert!(r.solid.is_some(), "one big box minus two internal boxes → one solid");
        clear_tree(&r.shape);
        clear_tree(&big.0);
        clear_tree(&s1.0);
        clear_tree(&s2.0);
        clear_tree(&comp.0);
    }

    #[test]
    fn boolean_compound_common() {
        let a1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let a2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let b = test_box_at(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(2.5, 1.0, 1.0));
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&[a1.0.clone(), a2.0.clone()]);
        let r = boolean_compound(&comp.0, &b.0, BoolOp::Common, 1e-6).expect("compound common ok");
        let v = box_vol(&r.shape);
        // a1∩b = 0.5 volume, a2∩b = 0.5 volume.
        assert!((v - 1.0).abs() < 0.15, "compound common volume {v} (expected ~1.0)");
        clear_tree(&r.shape);
        clear_tree(&comp.0);
        clear_tree(&a1.0);
        clear_tree(&a2.0);
        clear_tree(&b.0);
    }

    #[test]
    fn self_intersection_detected() {
        use occt_core::gp::{GpAx3, GpDir, GpPln};
        use occt_geom::GeomPlane;
        use std::sync::Arc;

        // A valid box has no self-intersections.
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = detect_self_intersections(&boxed.solid.0, 1e-6);
        assert!(!rep.found, "a box must not self-intersect");

        // A shell with a horizontal face (z=0) crossed by a vertical face
        // (y=0): the two surfaces intersect along a line through both faces.
        let bld = TopoBuilder::new();
        let h_wire = bld.make_wire(&[
            bld.make_edge_segment(&GpPnt::new(-1.0, -1.0, 0.0), &GpPnt::new(1.0, -1.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, -1.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, 1.0, 0.0), &GpPnt::new(-1.0, 1.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(-1.0, 1.0, 0.0), &GpPnt::new(-1.0, -1.0, 0.0)),
        ]);
        let f_h = bld.make_face(Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard()))), &[h_wire]);
        let v_wire = bld.make_wire(&[
            bld.make_edge_segment(&GpPnt::new(-1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 1.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, 0.0, 1.0), &GpPnt::new(-1.0, 0.0, 1.0)),
            bld.make_edge_segment(&GpPnt::new(-1.0, 0.0, 1.0), &GpPnt::new(-1.0, 0.0, 0.0)),
        ]);
        let pln_y = GpPln::new(
            GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 1.0, 0.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
                .unwrap(),
        );
        let f_v = bld.make_face(Arc::new(GeomPlane::new(pln_y)), &[v_wire]);
        let shell = bld.make_shell(&[f_h, f_v]);
        let rep2 = detect_self_intersections(&shell.0, 1e-6);
        assert!(rep2.found, "crossing faces must be detected as self-intersecting");
        assert!(rep2.edge_count >= 1, "at least one intersecting face pair");
        clear_tree(&shell.0);
        clear_tree(&boxed.solid.0);
    }

    #[test]
    fn boolean_with_check_repairs() {
        let (a, b) = overlapping_boxes();
        let r = boolean_with_check(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("with-check fuse ok");
        assert!(r.solid.is_some(), "fuse produces a solid");
        assert!(!r.shells.is_empty(), "result has a shell");
        assert!(r.shells.iter().all(shell_is_closed), "result shells are closed");
        // Hard failures (open shell / repair failure) must not occur; benign
        // diagnostics such as a volume cross-check warning are acceptable.
        let hard = r.warnings.iter().any(|w| w.contains("not closed") || w.contains("no shell"));
        assert!(!hard, "no hard failures: {:?}", r.warnings);
        let v = box_vol(&r.shape);
        assert!((v - 1.5).abs() < 0.05, "with-check fuse volume {v} (expected 1.5)");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn decompose_compound_top_level() {
        let bld = TopoBuilder::new();
        let shapes = disjoint_box_shapes();
        let comp = bld.make_compound_of(&shapes);
        let subs = decompose_compound(&comp.0);
        assert_eq!(subs.len(), 3, "compound of 3 shapes decomposes into 3");
        for s in &subs {
            assert!(!s.is_compound(), "sub-shapes are atomic");
        }
        // Non-compound input → itself.
        let single = decompose_compound(&shapes[0]);
        assert_eq!(single.len(), 1);
        // shape_components splits a solid into its shells.
        let comps = shape_components(&shapes[0]);
        assert!(!comps.is_empty());
        assert!(comps[0].is_solid());
        clear_tree(&comp.0);
        for s in &shapes {
            clear_tree(s);
        }
    }

    #[test]
    fn boolean_result_summary_string() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let s = boolean_result_summary(&r);
        assert!(s.contains(r.shape.shape_type().to_str()), "summary has the shape type: {s}");
        assert!(s.contains("faces"), "summary mentions faces: {s}");
        assert!(s.contains("solid") || s.contains("compound"), "summary describes the solidity: {s}");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn connected_components_splits_compound() {
        let shapes = disjoint_box_shapes();
        let bld = TopoBuilder::new();
        let comp = bld.make_compound_of(&shapes);
        let comps = connected_components(&comp.0, 1e-6);
        assert_eq!(comps.len(), 3, "compound of 3 disjoint boxes → 3 components");
        for c in &comps {
            assert!(c.is_solid(), "each component is a solid");
        }
        // A single solid has one component.
        let one = connected_components(&shapes[0], 1e-6);
        assert_eq!(one.len(), 1, "a single solid → one component");
        clear_tree(&comp.0);
        for s in &shapes {
            clear_tree(s);
        }
    }

    #[test]
    fn boolean_fold_sequence() {
        // ((a ∪ b) − c): fuse two boxes, then cut a third.
        let a = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b = test_box_at(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0));
        let c = test_box_at(&GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 0.5, 1.0));
        let shapes = vec![a.0.clone(), b.0.clone(), c.0.clone()];
        let ops = vec![BoolOp::Fuse, BoolOp::Cut];
        let r = boolean_fold(&shapes, &ops, 1e-6).expect("fold ok");
        // a∪b = [0,1.5]×[0,1]×[0,1] (vol 1.5); cut c = [0.5,1.5]×[0,0.5]×[0,1]
        // (vol 0.5) → 1.5 − 0.5 = 1.0.
        let v = box_vol(&r.shape);
        assert!((v - 1.0).abs() < 0.2, "fold volume {v} (expected ~1.0)");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
        clear_tree(&c.0);
    }

    #[test]
    fn boolean_cut_many_works() {
        let big = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 3.0, 3.0));
        let c1 = test_box_at(&GpPnt::new(0.5, 0.5, 0.5), &GpPnt::new(1.5, 1.5, 1.5));
        let c2 = test_box_at(&GpPnt::new(1.7, 0.5, 0.5), &GpPnt::new(2.7, 1.5, 1.5));
        let r = boolean_cut_many(&big.0, &[c1.0.clone(), c2.0.clone()], 1e-6).expect("cut many ok");
        let v = box_vol(&r.shape);
        assert!((v - 25.0).abs() < 0.2, "cut-many volume {v} (expected ~25.0)");
        clear_tree(&r.shape);
        clear_tree(&big.0);
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn boolean_result_validate_clean() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        // A closed, non-self-intersecting fuse has no hard structural issues.
        let issues = boolean_result_validate(&r, 1e-6);
        let hard: Vec<&String> = issues.iter().filter(|s| s.contains("not closed") || s.contains("no shell")).collect();
        assert!(hard.is_empty(), "no open-shell issues: {issues:?}");
        assert!(r.shells.iter().all(shell_is_closed), "fuse shell is closed");
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    // ------------------------------------------------------------------
    // Self-intersection repair + multi-result decomposition tests
    // ------------------------------------------------------------------

    /// A shell with a horizontal face (z = 0) crossed by a vertical face
    /// (y = 0): the two faces intersect along the x-axis segment.
    pub(super) fn crossing_shell() -> Shell {
        use occt_core::gp::{GpAx3, GpDir, GpPln};
        use occt_geom::GeomPlane;
        use std::sync::Arc;
        let bld = TopoBuilder::new();
        let h_wire = bld.make_wire(&[
            bld.make_edge_segment(&GpPnt::new(-1.0, -1.0, 0.0), &GpPnt::new(1.0, -1.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, -1.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, 1.0, 0.0), &GpPnt::new(-1.0, 1.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(-1.0, 1.0, 0.0), &GpPnt::new(-1.0, -1.0, 0.0)),
        ]);
        let f_h = bld.make_face(Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard()))), &[h_wire]);
        let v_wire = bld.make_wire(&[
            bld.make_edge_segment(&GpPnt::new(-1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 1.0)),
            bld.make_edge_segment(&GpPnt::new(1.0, 0.0, 1.0), &GpPnt::new(-1.0, 0.0, 1.0)),
            bld.make_edge_segment(&GpPnt::new(-1.0, 0.0, 1.0), &GpPnt::new(-1.0, 0.0, 0.0)),
        ]);
        let pln_y = GpPln::new(
            GpAx3::new(
                GpPnt::zero(),
                GpDir::new(0.0, 1.0, 0.0).unwrap(),
                &GpDir::new(1.0, 0.0, 0.0).unwrap(),
            )
            .unwrap(),
        );
        let f_v = bld.make_face(Arc::new(GeomPlane::new(pln_y)), &[v_wire]);
        bld.make_shell(&[f_h, f_v])
    }

    #[test]
    fn repair_noop_on_good_box() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = repair_self_intersections(&boxy.solid.0, 1e-6).expect("repair ok");
        assert_eq!(rep.fixed_faces, 0, "a valid box needs no face fixes");
        assert_eq!(rep.removed_faces, 0);
        assert!(rep.warnings.is_empty(), "no warnings: {:?}", rep.warnings);
        let v = box_vol(&rep.repaired);
        assert!((v - 1.0).abs() < 0.05, "box volume preserved {v}");
        clear_tree(&rep.repaired);
        clear_tree(&boxy.solid.0);
    }

    #[test]
    fn repair_detects_crossing() {
        let shell = crossing_shell();
        let rep = detect_self_intersections(&shell.0, 1e-6);
        assert!(rep.found, "crossing shell is flagged");
        let out = repair_self_intersections(&shell.0, 1e-6).expect("repair ok");
        assert!(out.fixed_faces >= 1, "at least one face fixed, got {}", out.fixed_faces);
        assert!(
            out.repaired.is_shell() || out.repaired.is_solid(),
            "repaired shape is a boundary (shell or solid)"
        );
        let no_longer_crossing = detect_self_intersections(&out.repaired, 1e-6);
        assert!(!no_longer_crossing.found, "crossing resolved after repair");
        clear_tree(&out.repaired);
        clear_tree(&shell.0);
    }

    #[test]
    fn repair_preserves_volume() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = repair_self_intersections(&boxy.solid.0, 1e-6).expect("repair ok");
        let v = box_vol(&rep.repaired);
        assert!((v - 1.0).abs() < 0.05, "repaired box volume {v} (expected ~1.0)");
        clear_tree(&rep.repaired);
        clear_tree(&boxy.solid.0);
    }

    #[test]
    fn boolean_repaired_warns_on_fix() {
        // A self-intersecting shell fused with a disjoint box keeps the
        // crossing inside the compound; the repaired wrapper flags it.
        let shell = crossing_shell();
        let boxy = BRepPrimBox::make_box_corner(&GpPnt::new(5.0, 5.0, 5.0), &GpPnt::new(6.0, 6.0, 6.0));
        let r = boolean_repaired(&shell.0, &boxy.solid.0, BoolOp::Fuse, 1e-6).expect("repaired fuse ok");
        assert!(
            r.warnings.iter().any(|w| w.contains("fixed")),
            "warning mentions fixed faces: {:?}",
            r.warnings
        );

        // A clean boolean has no repair warning.
        let (a, b) = overlapping_boxes();
        let r2 = boolean_repaired(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("clean fuse ok");
        assert!(
            !r2.warnings.iter().any(|w| w.contains("fixed")),
            "clean result has no repair warning: {:?}",
            r2.warnings
        );
        clear_tree(&r.shape);
        clear_tree(&r2.shape);
        clear_tree(&shell.0);
        clear_tree(&boxy.solid.0);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn boolean_repaired_matches_plain() {
        let (a, b) = overlapping_boxes();
        let plain = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("plain fuse");
        let repaired = boolean_repaired(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("repaired fuse");
        let vp = box_vol(&plain.shape);
        let vr = box_vol(&repaired.shape);
        assert!((vp - vr).abs() < 0.05, "repaired {vr} vs plain {vp}");
        clear_tree(&plain.shape);
        clear_tree(&repaired.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn decompose_multi_result_three_boxes() {
        let shapes = disjoint_box_shapes();
        let r = boolean_multi(&shapes, BoolOp::Fuse, 1e-6).expect("multi fuse ok");
        let m = decompose_multi_result(&r);
        assert_eq!(m.shapes.len(), 3, "three disjoint boxes decompose into 3, got {}", m.shapes.len());
        assert!(m.solids.len() >= 1, "at least one solid: {}", m.solids.len());
        assert_eq!(m.compounds.len(), 0, "no nested compounds remain");
        clear_tree(&r.shape);
        for s in &shapes {
            clear_tree(s);
        }
    }

    #[test]
    fn boolean_split_result_fuse_disjoint() {
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let m = boolean_split_result(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("split fuse ok");
        assert!(m.shapes.len() >= 2, "disjoint fuse produces multiple shapes: {}", m.shapes.len());
        assert!(m.solids.len() >= 2, "both pieces are solids: {}", m.solids.len());
        clear_tree(&b1.0);
        clear_tree(&b2.0);
    }

    #[test]
    fn boolean_split_result_cut_two_pieces() {
        // A 3×1×1 box cut by a thin full-width slab into two unit pieces.
        let box_a = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let slab = test_box_at(&GpPnt::new(1.0, -0.1, -0.1), &GpPnt::new(2.0, 1.1, 1.1));
        let m = boolean_split_result(&box_a.0, &slab.0, BoolOp::Cut, 1e-6).expect("split cut ok");
        assert!(m.shapes.len() >= 2, "cut into two pieces yields >=2 shapes, got {}", m.shapes.len());
        assert!(m.solids.len() >= 2, "two pieces are solids: {}", m.solids.len());
        let total: f64 = m.solids.iter().map(|s| box_vol(s)).sum();
        assert!((total - 2.0).abs() < 0.2, "total volume {total} (expected ~2.0)");
        clear_tree(&box_a.0);
        clear_tree(&slab.0);
    }

    #[test]
    fn result_is_multiple_disjoint_vs_overlap() {
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let r_disjoint = boolean(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("disjoint fuse");
        assert!(result_is_multiple(&r_disjoint), "disjoint fuse is multiple");

        let (a, b) = overlapping_boxes();
        let r_overlap = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("overlap fuse");
        assert!(!result_is_multiple(&r_overlap), "overlapping fuse is a single result");
        clear_tree(&r_disjoint.shape);
        clear_tree(&r_overlap.shape);
        clear_tree(&b1.0);
        clear_tree(&b2.0);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn multi_result_solids_and_shells() {
        use crate::brep_builder_api::make_face_from_polygon;
        let boxy = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let face = make_face_from_polygon(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ])
        .expect("face ok");
        let bld = TopoBuilder::new();
        let open_shell = bld.make_shell(&[face]);
        let comp = bld.make_compound_of(&[boxy.0.clone(), open_shell.0.clone()]);
        let r = BooleanResult {
            shape: comp.0.clone(),
            solid: None,
            shells: vec![open_shell],
            faces: vec![],
            warnings: vec![],
        };
        let m = decompose_multi_result(&r);
        assert!(m.shapes.len() >= 2, "solid + open shell decomposes into 2: {}", m.shapes.len());
        assert!(!m.solids.is_empty(), "solids populated");
        assert!(!m.shells.is_empty(), "shells populated");
        clear_tree(&comp.0);
        clear_tree(&boxy.0);
    }

    #[test]
    fn topologically_clean_removes_duplicates() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let v2 = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let comp = b.make_compound_of(&[v1.0.clone(), v2.0.clone()]);
        let before = vertices_of(&comp.0).len();
        let cleaned = topologically_clean(&comp.0, 1e-6);
        let after = vertices_of(&cleaned).len();
        assert!(after < before, "clean reduces vertex count: {before} -> {after}");
        assert_eq!(after, 1, "duplicate vertices merge into one");
        clear_tree(&comp.0);
        clear_tree(&cleaned);
    }

    #[test]
    fn analyze_self_intersections_reports_crossing() {
        let shell = crossing_shell();
        let issues = analyze_self_intersections(&shell.0, 1e-6);
        assert_eq!(issues.len(), 1, "one crossing pair, got {}", issues.len());
        assert_eq!(issues[0].kind, SelfIntersectionKind::Crossing);
        assert!(!issues[0].points.is_empty(), "crossing carries sample points");
        assert!(issues[0].segment.is_some(), "crossing has a clipped segment");
        let summary = self_intersection_issues_summary(&issues);
        assert!(summary.contains("1 crossing"), "summary: {summary}");

        // A valid box has no issues.
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(analyze_self_intersections(&boxy.solid.0, 1e-6).is_empty());
        clear_tree(&shell.0);
        clear_tree(&boxy.solid.0);
    }
