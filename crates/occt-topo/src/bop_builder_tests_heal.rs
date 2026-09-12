use super::*;
use super::tests::{box_vol, clear_tree, disjoint_box_shapes, overlapping_boxes, test_box_at};
use super::tests_api::crossing_shell;
use crate::builder::TopoBuilder;
use crate::primitives::BRepPrimBox;
use crate::shape::{Face, Shell};
use crate::tgeometry::GeometryRegistry;
use occt_core::gp::GpPnt;

    #[test]
    fn analyze_boundary_box_vs_crossing() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = analyze_boundary(&boxy.solid.0, 1e-6);
        assert_eq!(rep.shell_count, 1);
        assert_eq!(rep.closed_shells, 1);
        assert_eq!(rep.open_shells, 0);
        assert_eq!(rep.free_edges, 0);
        assert_eq!(rep.self_intersections, 0);
        assert_eq!(rep.degenerate_faces, 0);
        assert_eq!(rep.small_edges, 0);
        let summary = boundary_analysis_summary(&rep);
        assert!(summary.contains("1 closed"), "summary: {summary}");

        let shell = crossing_shell();
        let rep2 = analyze_boundary(&shell.0, 1e-6);
        assert!(rep2.self_intersections >= 1, "crossing shell reports a self-intersection");
        assert!(rep2.open_shells >= 1, "crossing shell is open");
        clear_tree(&boxy.solid.0);
        clear_tree(&shell.0);
    }

    #[test]
    fn multi_result_helpers() {
        let b1 = test_box_at(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 1.0));
        let b2 = test_box_at(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 1.0, 1.0));
        let m = boolean_split_result(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("split fuse ok");
        assert_eq!(m.shapes.len(), 2);
        assert_eq!(component_counts_by_type(&m), (2, 0, 0));
        let vols = component_volumes(&m);
        assert!((vols.iter().sum::<f64>() - 2.0).abs() < 0.1, "volumes {vols:?}");
        assert!((multi_result_total_volume(&m) - 2.0).abs() < 0.1);
        let big = largest_component(&m).expect("largest");
        assert!(box_vol(&big) > 0.9, "largest component has volume");
        let summary = multi_result_summary(&m);
        assert!(summary.contains("2 shape(s)"), "summary: {summary}");

        // Pipeline wrappers agree.
        let comps = boolean_components(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("components");
        assert_eq!(comps.len(), 2);
        let cvols = boolean_component_volumes(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("volumes");
        assert_eq!(cvols.len(), 2);
        let single = boolean_repair_and_clean(&b1.0, &b2.0, BoolOp::Fuse, 1e-6).expect("clean");
        assert!(single.is_compound(), "two disjoint boxes fuse to a compound");
        clear_tree(&b1.0);
        clear_tree(&b2.0);
    }

    #[test]
    fn repair_loop_converges() {
        let shell = crossing_shell();
        let out = repair_self_intersections_loop(&shell.0, 1e-6, 4);
        assert!(out.fixed_faces >= 1, "loop fixed at least one face");
        assert!(
            !detect_self_intersections(&out.repaired, 1e-6).found,
            "loop converges to a clean boundary"
        );
        // The repaired boundary is still a valid (possibly open) shell.
        assert!(out.repaired.is_shell() || out.repaired.is_solid());
        let summary = repair_result_summary(&out);
        assert!(summary.contains("fixed"), "summary: {summary}");
        clear_tree(&out.repaired);
        clear_tree(&shell.0);
    }

    #[test]
    fn validity_checks_box_clean() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert!(check_shape_validity(&boxy.solid.0, 1e-6).is_empty(), "box passes every check");
        assert!(shape_degenerate_faces(&boxy.solid.0, 1e-6).is_empty());
        assert!(shape_small_edges(&boxy.solid.0, 1e-6).is_empty());
        let conn = boundary_connectivity(&boxy.solid.0, 1e-6);
        assert_eq!(conn.len(), 1, "a box is one connected boundary");
        assert_eq!(conn[0].len(), 6, "all six faces connected");

        let shell = crossing_shell();
        let issues = check_shape_validity(&shell.0, 1e-6);
        assert!(
            issues.iter().any(|s| s.contains("self-intersection")),
            "crossing shell flags a self-intersection: {issues:?}"
        );
        clear_tree(&boxy.solid.0);
        clear_tree(&shell.0);
    }

    #[test]
    fn component_kind_classification() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        assert_eq!(component_kind(&boxy.solid.0), ComponentKind::Solid);
        assert_eq!(component_kind_label(ComponentKind::Solid), "solid");
        let face = crate::brep_builder_api::make_face_from_polygon(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ])
        .expect("face ok");
        assert_eq!(component_kind(&face.0), ComponentKind::Face);
        clear_tree(&boxy.solid.0);
        clear_tree(&face.0);
    }

    #[test]
    fn boolean_repaired_with_check_clean() {
        let (a, b) = overlapping_boxes();
        let r = boolean_repaired_with_check(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("with-check repaired fuse");
        assert!(r.solid.is_some(), "fuse produces a solid");
        let v = box_vol(&r.shape);
        assert!((v - 1.5).abs() < 0.05, "with-check repaired volume {v} (expected 1.5)");
        let (report_r, rep) = boolean_repaired_report(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("report fuse");
        assert!(report_r.solid.is_some());
        assert_eq!(rep.fixed_faces, 0, "clean boolean needs no repair");
        clear_tree(&r.shape);
        clear_tree(&report_r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn boolean_result_analysis_helpers() {
        let (a, b) = overlapping_boxes();
        let r = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        let report = boolean_result_boundary_report(&r, 1e-6);
        assert_eq!(report.closed_shells, 1, "fused box is one closed shell");
        assert_eq!(boolean_result_component_count(&r), 1, "single solid result");
        assert_eq!(boolean_result_shape_kind(&r), ComponentKind::Solid);
        let delta = boolean_result_volume_delta(&r, 1.5);
        assert!(delta < 0.05, "volume delta {delta}");
        assert!(boolean_results_close(&r, &r, 0.01), "a result is close to itself");
        assert!((boolean_result_total_volume(&r) - 1.5).abs() < 0.05);
        clear_tree(&r.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    #[test]
    fn component_edge_reports_box_and_crossing() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let reports = component_edge_reports(&boxy.solid.0, 1e-6);
        assert_eq!(reports.len(), 1, "a box is one component");
        assert_eq!(reports[0].face_count, 6);
        assert_eq!(reports[0].edge_count, 12);
        assert_eq!(reports[0].vertex_count, 8);
        assert_eq!(component_euler_characteristic(&reports[0]), 2, "closed box has χ = 2");

        let shell = crossing_shell();
        let reports2 = component_edge_reports(&shell.0, 1e-6);
        assert_eq!(reports2.len(), 2, "two crossing faces are two components");
        assert!(reports2.iter().all(|r| r.face_count == 1));
        clear_tree(&boxy.solid.0);
        clear_tree(&shell.0);
    }

    #[test]
    fn repair_opts_multi_pass() {
        let shell = crossing_shell();
        let opts = RepairOptions { tolerance: 1e-6, max_passes: 3, remove_coplanar_overlaps: true, report_open_boundaries: true };
        let out = repair_self_intersections_opts(&shell.0, &opts).expect("opts repair ok");
        assert!(out.fixed_faces >= 1, "multi-pass repair fixes the crossing");
        assert!(!detect_self_intersections(&out.repaired, 1e-6).found, "crossing resolved");

        // Default options: a single pass still fixes the crossing.
        let out2 = repair_self_intersections_opts(&shell.0, &RepairOptions::default()).expect("default repair ok");
        assert!(out2.fixed_faces >= 1);
        clear_tree(&out.repaired);
        clear_tree(&out2.repaired);
        clear_tree(&shell.0);
    }

    // ------------------------------------------------------------------
    // Phase 12 — BOPAlgo/TopOpeBRep depth: edge-overlap repair, face
    // splitting, edge classification and tolerance healing.
    // ------------------------------------------------------------------

    /// Does every shell under `shape` form a closed manifold boundary?
    fn closed_of(shape: &TopoShape) -> bool {
        let shells: Vec<Shell> = shapes_of(shape, ShapeType::Shell).into_iter().map(Shell).collect();
        !shells.is_empty() && shells.iter().all(shell_is_closed)
    }

    #[test]
    fn repair_edge_overlaps_noop_on_box() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let rep = repair_edge_overlaps(&boxy.solid.0, 1e-6).expect("repair ok");
        assert_eq!(rep.repaired_edges, 0, "no overlapping edges in a box");
        assert_eq!(rep.removed_edges, 0, "no edges removed from a box");
        assert_eq!(rep.welded_vertices, 0, "a box has no near-coincident vertex pairs");
        assert!(rep.warnings.is_empty(), "no warnings: {:?}", rep.warnings);
        let v = box_vol(&rep.repaired);
        assert!((v - 1.0).abs() < 0.05, "box volume preserved {v}");
        assert!(closed_of(&rep.repaired), "repaired box stays closed");
        clear_tree(&rep.repaired);
        clear_tree(&boxy.solid.0);
    }

    #[test]
    fn repair_edge_overlaps_splits_partial_overlap() {
        // Two coplanar triangles whose base edges partially overlap along the
        // x-axis: [0,1] and [0.5,1.5]. Both edges must split at 0.5 / 1.0.
        let f1 = crate::brep_builder_api::make_face_from_polygon(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ])
        .expect("face1");
        let f2 = crate::brep_builder_api::make_face_from_polygon(&[
            GpPnt::new(0.5, 0.0, 0.0),
            GpPnt::new(1.5, 0.0, 0.0),
            GpPnt::new(0.5, 1.0, 0.0),
        ])
        .expect("face2");
        let bld = TopoBuilder::new();
        let shell = bld.make_shell(&[f1, f2]);

        let rep = repair_edge_overlaps(&shell.0, 1e-6).expect("repair ok");
        assert!(rep.repaired_edges >= 2, "both overlapping edges are split, got {}", rep.repaired_edges);
        assert!(!closed_of(&rep.repaired), "two open faces form an open shell");
        // The overlapping sub-segment [0.5,1] is shared: distinct edge count rises
        // from 6 (two triangles) to 7 (the overlap resolved into three x-axis
        // segments plus four triangle side edges).
        let ec = edges_of(&rep.repaired).len();
        assert_eq!(ec, 7, "distinct edges after repair: {ec}");
        clear_tree(&rep.repaired);
        clear_tree(&shell.0);
    }

    #[test]
    fn repair_edge_overlaps_merges_coincident() {
        // Two identical coplanar square faces: every boundary edge has a
        // coincident twin. The duplicates merge and the corner vertices weld.
        let square = |p: &GpPnt| {
            crate::brep_builder_api::make_face_from_polygon(&[
                *p,
                GpPnt::new(p.x() + 1.0, p.y(), p.z()),
                GpPnt::new(p.x() + 1.0, p.y() + 1.0, p.z()),
                GpPnt::new(p.x(), p.y() + 1.0, p.z()),
            ])
            .expect("square face")
        };
        let f1 = square(&GpPnt::new(0.0, 0.0, 0.0));
        let f2 = square(&GpPnt::new(0.0, 0.0, 0.0));
        let bld = TopoBuilder::new();
        let shell = bld.make_shell(&[f1, f2]);

        let rep = repair_edge_overlaps(&shell.0, 1e-6).expect("repair ok");
        assert!(rep.removed_edges >= 4, "four coincident edges removed, got {}", rep.removed_edges);
        assert!(rep.welded_vertices >= 4, "four corner pairs welded, got {}", rep.welded_vertices);
        clear_tree(&rep.repaired);
        clear_tree(&shell.0);
    }

    #[test]
    fn split_faces_along_intersections_crossing() {
        // A horizontal face crossed by a vertical face: splitting along their
        // intersection turns the horizontal face into two sub-faces.
        let shell = crossing_shell();
        let faces_before = faces_of(&shell.0).len();
        let pairs = vec![(0usize, 1usize)];
        let split = split_faces_along_intersections(&shell.0, &pairs, 1e-6).expect("split ok");
        let faces_after = faces_of(&split).len();
        assert!(faces_after > faces_before, "splitting grows the face count: {faces_before} -> {faces_after}");
        assert!(faces_after >= 3, "crossing faces split into at least 3, got {faces_after}");
        clear_tree(&split);
        clear_tree(&shell.0);
    }

    #[test]
    fn classify_boolean_edges_box_all_external() {
        let boxy = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let classes = classify_boolean_edges(&boxy.solid.0, BoolOp::Fuse);
        assert_eq!(classes.len(), 12, "a box has 12 distinct edges");
        let counts = edge_class_counts(&classes);
        assert_eq!(counts.external, 12, "every box edge is external: {counts:?}");
        assert_eq!(counts.internal, 0);
        assert_eq!(counts.shared, 0);
        assert_eq!(counts.on_face, 0);
        let summary = edge_class_counts_summary(&counts);
        assert!(summary.contains("12 edge(s)"), "summary: {summary}");
        clear_tree(&boxy.solid.0);
    }

    #[test]
    fn classify_boolean_edges_shared_three_faces_on_one_edge() {
        // A fan of three triangle faces around a single shared edge: the shared
        // edge is referenced by 3 faces -> Shared; each triangle's two unique
        // edges are referenced by 1 face -> Internal.
        let bld = TopoBuilder::new();
        let e0 = bld.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let pln = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));
        let face_on = |apex: GpPnt| {
            let e1 = bld.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &apex);
            let e2 = bld.make_edge_segment(&apex, &GpPnt::new(0.0, 0.0, 0.0));
            let wire = bld.make_wire(&[e0.clone(), e1, e2]);
            bld.make_face(pln.clone(), &[wire])
        };
        let shell = bld.make_shell(&[
            face_on(GpPnt::new(0.0, 1.0, 0.0)),
            face_on(GpPnt::new(0.0, -1.0, 0.0)),
            face_on(GpPnt::new(0.0, 0.0, 1.0)),
        ]);
        let classes = classify_boolean_edges(&shell.0, BoolOp::Fuse);
        let counts = edge_class_counts(&classes);
        assert_eq!(counts.shared, 1, "the shared base edge is used by three faces: {counts:?}");
        assert_eq!(counts.internal, 6, "six unique triangle edges are single-face edges");
        assert_eq!(counts.shared + counts.on_face + counts.internal + counts.external, classes.len());
        clear_tree(&shell.0);
    }

    #[test]
    fn heal_tolerance_welds_and_removes() {
        // A near-closed wire: the last vertex is within tol of the first, so
        // welding merges them.
        let w = crate::brep_builder_api::make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1e-6, 1e-6, 0.0),
        ])
        .expect("wire");
        let rep = heal_tolerance_report(&w.0, 1e-3).expect("heal ok");
        assert!(rep.welded_vertices >= 1, "near-coincident vertices weld, got {}", rep.welded_vertices);
        // A tiny edge is removed by the small-edge pass.
        let w2 = crate::brep_builder_api::make_wire_from_points(&[
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(5e-4, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ])
        .expect("wire2");
        let rep2 = heal_tolerance_report(&w2.0, 1e-4).expect("heal2 ok");
        assert!(rep2.removed_edges >= 1, "tiny edge removed, got {}", rep2.removed_edges);
        // The TopoShape entry point returns just the healed shape.
        let healed = heal_tolerance(&w.0, 1e-3).expect("heal shape ok");
        assert!(healed.is_wire() || healed.is_compound(), "healed shape is a boundary");
        clear_tree(&rep.healed);
        clear_tree(&rep2.healed);
        clear_tree(&healed);
        clear_tree(&w.0);
        clear_tree(&w2.0);
    }

    #[test]
    fn phase12_integration_fuse_cut_split_repair_classify() {
        let (a, b) = overlapping_boxes();
        let fuse = boolean(&a.0, &b.0, BoolOp::Fuse, 1e-6).expect("fuse ok");
        assert!(fuse.solid.is_some(), "fuse is a solid");
        assert!(shell_is_closed(&fuse.shells[0]), "fuse shell is closed");
        let fuse_vol = box_vol(&fuse.shape);
        assert!((fuse_vol - 1.5).abs() < 0.1, "fuse volume {fuse_vol} (expected 1.5)");

        // Remaining face intersections in the clean fuse result (likely none).
        let issues = analyze_self_intersections(&fuse.shape, 1e-6);
        let pairs: Vec<(usize, usize)> = issues.iter().map(|i| (i.face_a, i.face_b)).collect();

        // Split faces along those intersections, then repair and heal.
        let split = split_faces_along_intersections(&fuse.shape, &pairs, 1e-6).expect("split ok");
        assert!(closed_of(&split), "split result stays closed");
        let rep = repair_edge_overlaps(&split, 1e-6).expect("repair ok");
        assert!(closed_of(&rep.repaired), "repaired result stays closed");
        let healed = heal_tolerance(&rep.repaired, 1e-6).expect("heal ok");
        assert!(closed_of(&healed), "healed result stays closed");

        // Volume is preserved through the whole pipeline.
        let v = box_vol(&healed);
        assert!((v - fuse_vol).abs() < 0.1, "volume preserved {v} vs {fuse_vol}");

        // Classify the repaired result: closed manifold -> no internal edges.
        let classes = classify_boolean_edges(&healed, BoolOp::Fuse);
        assert_eq!(classes.len(), edges_of(&healed).len(), "one class per edge");
        let counts = edge_class_counts(&classes);
        assert_eq!(counts.internal, 0, "closed result has no internal edges: {counts:?}");
        assert_eq!(counts.shared + counts.on_face + counts.internal + counts.external, classes.len());
        assert!(counts.external > 0, "outer boundary edges present");

        // The operand-aware classifier also reports a sane total.
        let op_classes = classify_edges_with_operands(&a.0, &b.0, &healed, BoolOp::Fuse, 1e-6);
        assert_eq!(op_classes.len(), classes.len());

        // Cut path: A − B is a closed solid too.
        let cut = boolean(&a.0, &b.0, BoolOp::Cut, 1e-6).expect("cut ok");
        assert!(cut.solid.is_some(), "cut produces a solid");
        let cut_vol = box_vol(&cut.shape);
        assert!((cut_vol - 0.5).abs() < 0.1, "cut volume {cut_vol} (expected 0.5)");
        assert!(closed_of(&cut.shape), "cut result stays closed");

        clear_tree(&split);
        clear_tree(&rep.repaired);
        clear_tree(&healed);
        clear_tree(&fuse.shape);
        clear_tree(&cut.shape);
        clear_tree(&a.0);
        clear_tree(&b.0);
    }

    // ------------------------------------------------------------------
    // Subtask 3 gate: hole-aware (FClass2d) ring-face classification
    // ------------------------------------------------------------------

    /// A ring face: outer square `0..2`, square hole `0.5..1.5`, in `z=0`.
    fn ring_face(b: &TopoBuilder) -> Face {
        use crate::builder_face::build_face_with_holes;
        let mk = |pts: &[GpPnt]| {
            (0..4).map(|i| b.make_edge_segment(&pts[i], &pts[(i + 1) % 4])).collect::<Vec<Edge>>()
        };
        let outer = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
            GpPnt::new(2.0, 2.0, 0.0),
            GpPnt::new(0.0, 2.0, 0.0),
        ];
        let hole = [
            GpPnt::new(0.5, 0.5, 0.0),
            GpPnt::new(0.5, 1.5, 0.0),
            GpPnt::new(1.5, 1.5, 0.0),
            GpPnt::new(1.5, 0.5, 0.0),
        ];
        build_face_with_holes(&mk(&outer), &[mk(&hole)]).expect("ring face")
    }

    /// The ported `IntTools_FClass2d` classifies a ring face's hole vs its
    /// material correctly: hole-interior points are `Out`, ring-material points
    /// are `In`, both boundary rings are `On`, and the face is a growth (not a
    /// hole).
    #[test]
    fn fclass2d_classifies_ring_face_hole_vs_material() {
        use crate::fclass2d::{FaceState, FClass2d};
        let b = TopoBuilder::new();
        let ring = ring_face(&b);
        let pln = face_plane_local(&ring).expect("plane");
        let cl = FClass2d::new(&ring, 1e-7).expect("classifier");
        assert!(!cl.is_hole(), "ring face is a growth, not a hole");
        let p = |x: f64, y: f64| project_point_to_plane(&pln, &GpPnt::new(x, y, 0.0));
        // Ring material points → In.
        assert_eq!(cl.perform(p(0.1, 0.1)), FaceState::In, "ring corner material");
        assert_eq!(cl.perform(p(1.9, 1.9)), FaceState::In, "ring far corner material");
        // Hole points → Out (a hole-interior point is not part of the face).
        assert_eq!(cl.perform(p(1.0, 1.0)), FaceState::Out, "hole center");
        assert_eq!(cl.perform(p(0.75, 0.75)), FaceState::Out, "hole interior");
        // Boundaries → On.
        assert_eq!(cl.perform(p(0.0, 1.0)), FaceState::On, "outer boundary");
        assert_eq!(cl.perform(p(0.5, 1.0)), FaceState::On, "hole boundary");
        // Outside → Out.
        assert_eq!(cl.perform(p(3.0, 3.0)), FaceState::Out, "outside");
        clear_tree(&ring);
    }

    /// The boolean's own hole-aware point-in-face test (`point_in_face_holes`,
    /// which mirrors `FClass2d`'s region semantics in the plane frame) agrees:
    /// a hole-interior point is not "on" the face, a ring-material point is.
    #[test]
    fn boolean_classify_ring_face_hole_vs_material() {
        let b = TopoBuilder::new();
        let ring = ring_face(&b);
        let pln = face_plane_local(&ring).expect("plane");
        for (label, x, y, expect) in [
            ("ring material (0.1,0.1)", 0.1, 0.1, true),
            ("ring material (1.9,1.9)", 1.9, 1.9, true),
            ("hole center (1,1)", 1.0, 1.0, false),
            ("hole interior (0.75,0.75)", 0.75, 0.75, false),
            ("outside (3,3)", 3.0, 3.0, false),
        ] {
            let on = point_in_face_holes(&ring, &pln, &GpPnt::new(x, y, 0.0), 1e-6);
            assert_eq!(on, expect, "{label}");
        }
        clear_tree(&ring);
    }

    /// `PerformShapesToAvoid` translation: the arrangement strips dangling
    /// section segments (degree-1 endpoints) that can never close a loop, so a
    /// slit or a floating segment leaves the face's regions intact instead of
    /// stalling the half-edge traversal.
    #[test]
    fn trace_regions_strips_dangling_section_segment() {
        let sq = [
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(0.0, 1.0),
        ];
        // A segment fully floating inside the face: both ends dangling → stripped.
        let regions = trace_planar_regions(&sq, &[(GpPnt2d::new(0.3, 0.3), GpPnt2d::new(0.7, 0.7))]);
        assert_eq!(regions.len(), 1, "floating segment does not split the face");
        // A slit from the bottom boundary inward: the inner end is dangling.
        let regions = trace_planar_regions(&sq, &[(GpPnt2d::new(0.5, 0.0), GpPnt2d::new(0.5, 0.4))]);
        assert_eq!(regions.len(), 1, "boundary slit does not split the face");
        // A closed loop still splits the face (dangling strips are isolated).
        let loop_pts = [
            GpPnt2d::new(0.25, 0.25),
            GpPnt2d::new(0.75, 0.25),
            GpPnt2d::new(0.75, 0.75),
            GpPnt2d::new(0.25, 0.75),
        ];
        let mut closed: Vec<(GpPnt2d, GpPnt2d)> = Vec::new();
        for i in 0..4 {
            closed.push((loop_pts[i], loop_pts[(i + 1) % 4]));
        }
        closed.push((GpPnt2d::new(0.1, 0.1), GpPnt2d::new(0.2, 0.9))); // dangling diagonal
        let regions = trace_planar_regions(&sq, &closed);
        assert_eq!(regions.len(), 2, "closed loop splits into ring + disk");
    }

