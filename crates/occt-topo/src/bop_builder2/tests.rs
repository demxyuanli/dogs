use super::prelude::*;
use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::primitives::BRepPrimBox;
    use occt_core::gp::GpPnt;

    /// A second unit box at `x ∈ [3, 4]` — disjoint from the `[0, 1]³` box.
    fn far_box() -> TopoShape {
        BRepPrimBox::make_box_corner(&GpPnt::new(3.0, 0.0, 0.0), &GpPnt::new(4.0, 1.0, 1.0)).solid.0
    }

    /// Counts the direct solid children of a compound.
    fn count_solids(shape: &TopoShape) -> usize {
        shape
            .tshape
            .read()
            .unwrap()
            .children
            .iter()
            .filter(|h| h.shape_type() == ShapeType::Solid)
            .count()
    }

    #[test]
    fn defaults() {
        let b = BopBuilder::new();
        assert!(b.arguments().is_empty());
        assert!(!b.has_errors());
        assert!(b.errors().is_empty());
        assert!(b.warnings().is_empty());
        assert!(b.history().is_empty());
        assert!(b.origins().is_empty());
        assert!(b.shapes_sd().is_empty());
        assert_eq!(b.ds().nb_shapes(), 0);
        assert!(b.result().is_compound());
    }

    #[test]
    fn set_arguments_dedupes_by_tshape() {
        let a = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), a.solid.0.clone()]);
        assert_eq!(b.arguments().len(), 1);
        b.add_argument(&a.solid.0);
        assert_eq!(b.arguments().len(), 1);
        // A freshly built box is a distinct TShape and is added.
        let c = unit_box();
        b.add_argument(&c.solid.0);
        assert_eq!(b.arguments().len(), 2);
    }

    #[test]
    fn too_few_arguments_fails() {
        let mut b = BopBuilder::new();
        let err = b.perform().unwrap_err();
        assert!(err.contains("too few"), "err: {err}");
        assert!(b.has_errors());
    }

    #[test]
    fn perform_on_two_overlapping_boxes_succeeds() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        let result = b.perform().unwrap();
        assert!(!b.has_errors(), "errors: {:?}", b.errors());
        assert!(result.is_compound() || result.is_solid(), "type {:?}", result.shape_type());
        // The GF result is a compound of the two (unmodified) solids while the
        // face/solid stages are stubbed; the history still recorded the
        // vertex/edge images produced by the intersection.
        assert!(!b.history().is_empty(), "history should record split images");
        assert!(!b.origins().is_empty(), "origins back-map should be populated");
    }

    #[test]
    fn prepare_history_fills_modified_generated_removed() {
        // Fuse of two identical boxes: the coincident vertex/edge/face images
        // become modified relations, all kept in the result (`PrepareHistory`,
        // `BOPAlgo_Builder_4.cxx`).
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.perform().unwrap();
        let m = b.history().modified_map();
        assert!(!m.is_empty(), "identical boxes must record modified relations");
        let result_sub: Vec<TopoShape> = all_subshapes(b.result());
        for splits in m.values() {
            for sp in splits {
                assert!(
                    result_sub.iter().any(|r| r.same_tshape(sp)),
                    "a modified split must be kept in the result"
                );
            }
        }
        // Any generated relation also references a result shape.
        for gens in b.history().generated_map().values() {
            for gg in gens {
                assert!(result_sub.iter().any(|r| r.same_tshape(gg)));
            }
        }
        // Cut of a disjoint tool: the tool solid has no trace in the result and
        // no surviving splits, so it is marked removed.
        let a = unit_box();
        let bx = far_box();
        let mut b2 = BopBuilder::new();
        b2.set_arguments(&[a.solid.0.clone(), bx.clone()]);
        let (os, ts) = BoolOp2::Cut.states();
        b2.perform_internal(&[a.solid.0.clone()], os, &[bx.clone()], ts).unwrap();
        assert!(b2.history().is_deleted(&bx), "the cut-away tool solid is removed");
    }

    #[test]
    fn fill_images_vertices_maps_same_domain_corners() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        // The 8 coincident corners of the two identical boxes are same-domain.
        assert!(!b.ds().shapes_sd().is_empty(), "coincident corners must be same-domain");
        b.fill_images_vertices().unwrap();
        let sd = b.ds().shapes_sd().clone();
        for (&n_v, &n_vsd) in &sd {
            let v = b.ds().shape(n_v).expect("vertex shape");
            let vsd = b.ds().shape(n_vsd).expect("sd vertex shape");
            let img = b.history().image(v).expect("vertex has an image");
            assert!(!img.is_empty());
            assert!(img[0].same_tshape(vsd), "image must be the same-domain vertex");
            let ors = b.origins().get(&shape_key(vsd)).expect("origin recorded");
            assert!(ors.iter().any(|o| o.same_tshape(v)), "origins back-map must point at the source");
        }
    }

    #[test]
    fn fill_images_edges_records_split_edges() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        b.fill_images_edges().unwrap();
        // Every source edge that carries pave blocks has a recorded image
        // (its split piece — the whole edge when no split was needed).
        let n = b.ds().nb_source_shapes();
        let mut n_edge_images = 0;
        for i in 0..n {
            let Some(si) = b.ds().shape_info(i) else { continue };
            if si.shape_type() != ShapeType::Edge {
                continue;
            }
            if !b.ds().has_pave_blocks(i) {
                continue;
            }
            let e = si.shape();
            let img = b.history().image(e).expect("edge has an image");
            assert!(!img.is_empty());
            n_edge_images += 1;
        }
        assert!(n_edge_images >= 24, "12 box edges × 2 boxes = 24 source edges");
    }

    #[test]
    fn builder_bop_fuse_separate_boxes_two_bodies() {
        let a = unit_box();
        let bx = far_box();
        let result = builder_bop(&[a.solid.0.clone()], &[bx], BoolOp2::Fuse).unwrap();
        assert!(result.is_compound() || result.is_solid());
        assert_eq!(count_solids(&result), 2, "fuse of two disjoint boxes keeps both solids");
    }

    #[test]
    fn builder_bop_cut_and_common_run_without_error() {
        let a = unit_box();
        let bx = far_box();
        let cut = builder_bop(&[a.solid.0.clone()], &[bx.clone()], BoolOp2::Cut).unwrap();
        assert!(!cut.is_null());
        let cmn = builder_bop(&[a.solid.0.clone()], &[bx], BoolOp2::Common).unwrap();
        assert!(!cmn.is_null());
    }

    #[test]
    fn history_and_origins_are_accessible() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.perform().unwrap();
        assert!(!b.history().images().is_empty());
        // history_mut allows the sibling Phase-20 modules to add images.
        b.history_mut().add_image(&a.solid.0, c.solid.0.clone());
        assert!(b.history().has_image(&a.solid.0));
    }

    #[test]
    fn trait_support_methods_are_available() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        // Options are forwarded from the PaveFiller.
        assert_eq!(b.fuzzy_value(), 1e-7);
        assert!(!b.non_destructive());
        // origins_mut exposes the back-map for the sibling modules.
        b.origins_mut().entry(shape_key(&c.solid.0)).or_default().push(a.solid.0.clone());
        assert!(b.origins().get(&shape_key(&c.solid.0)).is_some());
        // Same-domain binding round-trips through the DS.
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        let sd = b.ds().shapes_sd().clone();
        if let Some((&n, &m)) = sd.iter().next() {
            let s = b.ds().shape(n).unwrap().clone();
            let rep = b.ds().shape(m).unwrap().clone();
            b.bind_shapes_sd(s.clone(), rep.clone());
            let found = b.seek_shapes_sd(&s).expect("same-domain representative");
            assert!(found.same_tshape(&rep));
        }
    }

    #[test]
    fn bool_op_states_follow_occt_conversion() {
        assert_eq!(BoolOp2::Fuse.states(), (FaceState::Out, FaceState::Out));
        assert_eq!(BoolOp2::Cut.states(), (FaceState::Out, FaceState::In));
        assert_eq!(BoolOp2::Common.states(), (FaceState::In, FaceState::In));
    }

    #[test]
    fn build_bop_rejects_unknown_shapes_and_bad_states() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        // A shape that is not an argument of the operation is rejected.
        let err = b
            .build_bop(&[c.solid.0.clone()], FaceState::Out, &[], FaceState::Out)
            .unwrap_err();
        assert!(err.contains("unknown shape"), "err: {err}");

        // An invalid state is rejected. Use a fresh builder — the failure above
        // left the report dirty, and `build_bop` short-circuits on a dirty
        // report before validating the states.
        let mut b2 = BopBuilder::new();
        b2.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b2.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b2.filler_mut().perform().unwrap();
        let err2 = b2
            .build_bop(&[a.solid.0.clone()], FaceState::On, &[], FaceState::Out)
            .unwrap_err();
        assert!(err2.contains("invalid state"), "err: {err2}");
    }

    #[test]
    fn clear_resets_state_keeps_nothing() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.perform().unwrap();
        assert!(!b.history().is_empty());
        b.clear();
        assert!(b.arguments().is_empty());
        assert!(b.history().is_empty());
        assert!(!b.has_errors());
    }

    #[test]
    fn bop_builder_implements_build_ops_and_like() {
        // The two Phase-20 host traits must be implemented and reachable
        // through the generic entry points.
        fn via_ops<B: BopBuildOps>(f: &mut B) -> usize {
            f.ds().nb_shapes() + f.fuzzy_value() as usize
        }
        fn via_like<B: BopBuilderLike>(f: &mut B) -> bool {
            f.has_errors()
        }
        let mut b = BopBuilder::new();
        assert_eq!(via_ops(&mut b), 0, "empty DS, fuzzy 1e-7 floors to 0");
        assert!(!via_like(&mut b));

        // The mutable accessors of both traits delegate to the builder fields.
        let a = unit_box();
        b.history_mut().add_image(&a.solid.0, a.solid.0.clone());
        assert!(b.history().has_image(&a.solid.0));
        b.origins_mut().entry(shape_key(&a.solid.0)).or_default().push(a.solid.0.clone());
        assert!(b.origins().get(&shape_key(&a.solid.0)).is_some());
    }

    #[test]
    fn fill_images_faces_callable_through_trait() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        // The face stage runs directly through the trait surface.
        crate::bop_build_faces::fill_images_faces(&mut b).unwrap();
        assert!(!b.has_errors(), "errors: {:?}", b.errors());
    }

    #[test]
    fn fill_images_solids_callable_through_trait() {
        let a = unit_box();
        let c = unit_box();
        let mut b = BopBuilder::new();
        b.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        b.filler_mut().perform().unwrap();
        // The solid stage runs directly through the trait surface.
        crate::bop_build_solids::fill_images_solids(&mut b).unwrap();
        assert!(!b.has_errors(), "errors: {:?}", b.errors());
    }

    #[test]
    fn perform_on_two_make_box_boxes_stays_valid() {
        // The validation case: two overlapping boxes (built with
        // BRepPrimBox::make_box) through the full General Fuse pipeline, with
        // the face and solid stages wired. Two coincident unit boxes: the
        // split-solid stage rebuilds at least one solid image whose union
        // covers the coincident pair.
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let c = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut br = BopBuilder::new();
        br.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        let result = br.perform().unwrap();
        assert!(!br.has_errors(), "errors: {:?}", br.errors());
        assert!(result.is_compound() || result.is_solid(), "type {:?}", result.shape_type());
        assert!(!result.is_null());
        // At least one argument solid was rebuilt into an image by the
        // split-solid stage (coincident faces produce section edges), and the
        // pipeline completes without errors.
        let any_image = br
            .arguments()
            .iter()
            .any(|arg| arg.is_solid() && br.history().has_image(arg));
        assert!(any_image, "a solid image was recorded");
    }
