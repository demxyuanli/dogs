use super::prelude::*;
use std::collections::HashMap;
    use std::sync::Arc;

    use occt_core::gp::{GpAx3, GpPln, GpPnt};
    use occt_geom::{GeomPlane, Surface};

    use super::*;
    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::primitives::BRepPrimBox;
    use crate::shape::{Edge, Face, Vertex};
    use crate::topo_tools_full::{faces_of, shapes_of};

    // -----------------------------------------------------------------------
    // Stub host
    // -----------------------------------------------------------------------

    struct StubBuilder {
        ds: BopdsDS,
        history: BopHistory,
        fuzzy: f64,
        args: Vec<TopoShape>,
        origins: HashMap<usize, Vec<TopoShape>>,
    }

    impl BopBuildOps for StubBuilder {
        fn ds(&self) -> &BopdsDS {
            &self.ds
        }
        fn history(&self) -> &BopHistory {
            &self.history
        }
        fn history_mut(&mut self) -> &mut BopHistory {
            &mut self.history
        }
        fn fuzzy_value(&self) -> f64 {
            self.fuzzy
        }
        fn arguments(&self) -> &[TopoShape] {
            &self.args
        }
        fn origins_mut(&mut self) -> &mut HashMap<usize, Vec<TopoShape>> {
            &mut self.origins
        }
    }

    fn stub(ds: BopdsDS, history: BopHistory, args: Vec<TopoShape>) -> StubBuilder {
        StubBuilder { ds, history, fuzzy: 1e-7, args, origins: HashMap::new() }
    }

    /// The face of `faces` whose boundary-vertex mean lies at height `z`.
    fn face_at_z(faces: &[Face], z: f64) -> Face {
        faces
            .iter()
            .find(|f| {
                let vs = crate::topo_tools_full::vertices_of(&f.0);
                if vs.is_empty() {
                    return false;
                }
                let zavg =
                    vs.iter().map(|v| BRepTool::vertex_point(v).z()).sum::<f64>() / vs.len() as f64;
                (zavg - z).abs() < 1e-9
            })
            .cloned()
            .expect("face at height")
    }

    /// A fresh face on the same surface and boundary edges as `src` (a new
    /// TShape, so it is a distinct split image).
    fn re_face(src: &Face) -> Face {
        let b = TopoBuilder::new();
        let edges = crate::topo_tools_full::edges_of(&src.0);
        let wire = b.make_wire(&edges);
        let surf = BRepTool::face_surface(src).expect("face surface");
        b.make_face(surf, &[wire])
    }

    // -----------------------------------------------------------------------
    // fill_images_containers
    // -----------------------------------------------------------------------

    #[test]
    fn empty_builder_leaves_history_untouched() {
        let ds = BopdsDS::new();
        let history = BopHistory::new();
        let mut b = stub(ds.clone(), history, Vec::new());
        fill_images_containers(&mut b, ShapeType::Shell).unwrap();
        assert!(!b.history().has_any_images());
    }

    #[test]
    fn box_shell_with_split_top_face_reassembles_closed_shell() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shell = shapes_of(&boxed.solid.0, ShapeType::Shell)[0].clone();
        let faces = faces_of(&shell);
        assert_eq!(faces.len(), 6);

        // Split the top face: register a fresh, geometrically identical face as
        // its only image.
        let top = face_at_z(&faces, 1.0);
        let new_top = re_face(&top);
        assert!(!new_top.0.same_tshape(&top.0));

        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top.0, new_top.0.clone());
        let mut b = stub(ds, history, vec![boxed.solid.0.clone()]);

        fill_images_containers(&mut b, ShapeType::Shell).unwrap();

        let imgs = b.history().image(&shell).expect("shell has an image");
        assert_eq!(imgs.len(), 1, "one closed shell image");
        let img = &imgs[0];
        assert!(img.is_shell(), "image is a shell");
        assert!(img.closed(), "rebuilt shell is closed");
        assert_eq!(faces_of(img).len(), 6, "all six split faces present");
    }

    #[test]
    fn unmodified_box_shell_gets_no_image() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shell = shapes_of(&boxed.solid.0, ShapeType::Shell)[0].clone();
        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut b = stub(ds, BopHistory::new(), vec![boxed.solid.0.clone()]);
        fill_images_containers(&mut b, ShapeType::Shell).unwrap();
        assert!(
            b.history().image(&shell).is_none(),
            "no face was split -> no container image"
        );
    }

    // -----------------------------------------------------------------------
    // fill_images_compounds
    // -----------------------------------------------------------------------

    #[test]
    fn compound_with_split_child_is_rebuilt() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e2_new = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let comp = b.make_compound_of(&[e1.0.clone(), e2.0.clone()]);

        let mut ds = BopdsDS::new();
        ds.init(&[comp.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&e2.0, e2_new.0.clone());
        let mut stub = stub(ds, history, vec![comp.0.clone()]);

        fill_images_compounds(&mut stub).unwrap();

        let imgs = stub.history().image(&comp.0).expect("compound has an image");
        assert_eq!(imgs.len(), 1);
        let img = &imgs[0];
        assert!(img.is_compound());
        let kids = direct_children(img);
        assert!(kids.iter().any(|k| k.same_tshape(&e1.0)), "unmodified child kept");
        assert!(kids.iter().any(|k| k.same_tshape(&e2_new.0)), "split child replaced");
        assert!(
            !kids.iter().any(|k| k.same_tshape(&e2.0)),
            "original split child not in the image"
        );
    }

    #[test]
    fn unmodified_compound_gets_no_image() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 0.0);
        let v2 = b.make_vertex(GpPnt::new(1.0, 0.0, 0.0), 0.0);
        let comp = b.make_compound_of(&[v1.0, v2.0]);
        let mut ds = BopdsDS::new();
        ds.init(&[comp.0.clone()]);
        let mut stub = stub(ds, BopHistory::new(), vec![comp.0.clone()]);
        fill_images_compounds(&mut stub).unwrap();
        assert!(stub.history().image(&comp.0).is_none());
    }

    // -----------------------------------------------------------------------
    // fill_internal_vertices
    // -----------------------------------------------------------------------

    fn square_face_with_alone_vertices() -> (Face, Face, Vertex, Vertex) {
        let b = TopoBuilder::new();
        let p = [GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0), GpPnt::new(1.0, 1.0, 0.0), GpPnt::new(0.0, 1.0, 0.0)];
        let edges: Vec<Edge> = (0..4)
            .map(|i| b.make_edge_segment(&p[i], &p[(i + 1) % 4]))
            .collect();
        let wire = b.make_wire(&edges);
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(GpPln::new(GpAx3::standard())));
        let mut face = b.make_face(surf.clone(), &[wire]);

        // A fresh split face with the same boundary.
        let wire2 = b.make_wire(&edges);
        let face_im = b.make_face(surf, &[wire2]);

        // Two alone vertices: one inside the square, one outside.
        let v_in = b.make_vertex(GpPnt::new(0.5, 0.5, 0.0), 0.0);
        let v_out = b.make_vertex(GpPnt::new(2.0, 2.0, 0.0), 0.0);
        b.add(&mut face.0, &v_in.0);
        b.add(&mut face.0, &v_out.0);
        (face, face_im, v_in, v_out)
    }

    #[test]
    fn alone_vertices_inside_split_are_added_as_internal() {
        let (face, face_im, v_in, v_out) = square_face_with_alone_vertices();

        let mut ds = BopdsDS::new();
        ds.init(&[face.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&face.0, face_im.0.clone());
        let mut b = stub(ds, history, vec![face.0.clone()]);

        fill_internal_vertices(&mut b).unwrap();

        // The flat model stores children without a per-child orientation marker,
        // so the INTERNAL annotation set by the port is not observable; what
        // matters is that the inside vertex became a child of the face image
        // while the outside one did not.
        let kids = direct_children(&face_im.0);
        assert!(
            kids.iter().any(|k| k.same_tshape(&v_in.0)),
            "inside alone vertex added as a child of the split face"
        );
        assert!(
            !kids.iter().any(|k| k.same_tshape(&v_out.0)),
            "outside alone vertex dropped"
        );
    }

    #[test]
    fn no_alone_vertices_adds_nothing() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = faces_of(&boxed.solid.0);
        let top = face_at_z(&faces, 1.0);
        let new_top = re_face(&top);
        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top.0, new_top.0.clone());
        let mut b = stub(ds, history, vec![boxed.solid.0.clone()]);
        fill_internal_vertices(&mut b).unwrap();
        // The box faces have no alone vertices, so no INTERNAL vertex is added
        // to the split top face (its wire child remains the only child).
        let kids = direct_children(&new_top.0);
        assert!(
            !kids.iter().any(|k| k.is_vertex() && k.orientation() == Orientation::Internal),
            "no internal vertex added to the split face"
        );
    }

    // -----------------------------------------------------------------------
    // fill_internal_shapes
    // -----------------------------------------------------------------------

    #[test]
    fn inside_vertex_settles_into_original_solid_as_copy() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = TopoBuilder::new();
        let v_in = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 0.0);
        let mut solid = boxed.solid.0.clone();
        b.add(&mut solid, &v_in.0); // internal vertex child of the solid

        let mut ds = BopdsDS::new();
        ds.init(&[solid.clone()]);
        let mut stub = stub(ds, BopHistory::new(), vec![solid.clone()]);

        fill_internal_shapes(&mut stub).unwrap();

        let imgs = stub.history().image(&solid).expect("solid gains a copy image");
        assert_eq!(imgs.len(), 1);
        let img = &imgs[0];
        assert!(img.is_solid());
        // The flat model does not persist a per-child orientation marker, so the
        // assertion is on presence of the vertex child (the INTERNAL annotation
        // set by the port is not observable through the child list).
        let kids = direct_children(img);
        assert!(
            kids.iter().any(|k| k.same_tshape(&v_in.0)),
            "internal vertex present in the solid copy"
        );
        // The original solid is preserved (only one shell, no internal vertex added).
        assert!(direct_children(&solid).iter().filter(|c| c.is_shell()).count() >= 1);
    }

    #[test]
    fn outside_vertex_is_not_settled() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let b = TopoBuilder::new();
        let v_out = b.make_vertex(GpPnt::new(5.0, 5.0, 5.0), 0.0);
        let mut solid = boxed.solid.0.clone();
        b.add(&mut solid, &v_out.0);

        let mut ds = BopdsDS::new();
        ds.init(&[solid.clone()]);
        let mut stub = stub(ds, BopHistory::new(), vec![solid.clone()]);

        fill_internal_shapes(&mut stub).unwrap();

        assert!(
            stub.history().image(&solid).is_none(),
            "no vertex lies inside the solid -> no split solid"
        );
    }

    // -----------------------------------------------------------------------
    // build_draft_solid
    // -----------------------------------------------------------------------

    #[test]
    fn box_shell_wraps_into_unit_solid() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shell = shapes_of(&boxed.solid.0, ShapeType::Shell)[0].clone();
        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut b = stub(ds, BopHistory::new(), vec![boxed.solid.0.clone()]);

        let solid = build_draft_solid(&mut b, &shell).unwrap();
        assert!(solid.is_solid());
        let v = crate::brep_gprop::volume(&solid, 0.02);
        assert!((v - 1.0).abs() < 0.01, "unit box volume, got {v}");
    }

    #[test]
    fn split_box_shell_rebuilds_closed_shell_into_solid() {
        let boxed = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shell = shapes_of(&boxed.solid.0, ShapeType::Shell)[0].clone();
        let faces = faces_of(&shell);
        let top = face_at_z(&faces, 1.0);
        let new_top = re_face(&top);

        let mut ds = BopdsDS::new();
        ds.init(&[boxed.solid.0.clone()]);
        let mut history = BopHistory::new();
        history.add_image(&top.0, new_top.0.clone());
        let mut b = stub(ds, history, vec![boxed.solid.0.clone()]);

        // First reassemble the shell image.
        fill_images_containers(&mut b, ShapeType::Shell).unwrap();
        let shell_imgs = b.history().image(&shell).unwrap().to_vec();
        assert_eq!(shell_imgs.len(), 1);

        // Then build a solid from the rebuilt (closed) shell.
        let solid = build_draft_solid(&mut b, &shell_imgs[0]).unwrap();
        assert!(solid.is_solid());
        let v = crate::brep_gprop::volume(&solid, 0.02);
        assert!((v - 1.0).abs() < 0.01, "split box still has unit volume, got {v}");
    }

    #[test]
    fn coincident_boxes_draft_solid_volume_is_one() {
        // Regression: two fully coincident boxes (same make_box parameters)
        // went through FillSameDomainFaces + BuildDraftSolid with their volume
        // inflated to ~17. The coincident faces were over-split into degenerate
        // sliver pieces and the draft-solid assembly counted every piece; both
        // the degenerate-piece skip and the same-domain self-image guard bring
        // each rebuilt box back to volume 1.
        let a = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let c = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut br = crate::bop_builder2::BopBuilder::new();
        br.set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        br.filler_mut().set_arguments(&[a.solid.0.clone(), c.solid.0.clone()]);
        br.filler_mut().perform().unwrap();
        crate::bop_build_faces::fill_images_faces(&mut br).unwrap();
        assert!(!br.has_errors(), "errors: {:?}", br.errors());

        for solid in [a.solid.0.clone(), c.solid.0.clone()] {
            let shell = shapes_of(&solid, ShapeType::Shell)[0].clone();
            let draft = build_draft_solid(&mut br, &shell).unwrap();
            let v = crate::brep_gprop::volume(&draft, 0.02);
            assert!((v - 1.0).abs() < 0.01, "coincident box draft volume, got {v}");
        }
    }
