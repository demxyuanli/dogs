use super::prelude::*;
use super::*;
    use crate::abs::ShapeType;
    use crate::brep_tool::BRepTool;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::topo_tools_full::{faces_of, vertices_of};
    use occt_core::gp::GpVec;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static TMP_SEQ: AtomicUsize = AtomicUsize::new(0);

    fn write_temp_step(content: &str) -> String {
        let n = TMP_SEQ.fetch_add(1, Ordering::SeqCst);
        let path = { let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target").join("test-out"); let _ = std::fs::create_dir_all(&d); d }.join(format!(
            "occt_xcaf_doc_{}_{n}.step",
            std::process::id()
        ));
        std::fs::write(&path, content).expect("write temp step");
        path.to_str().unwrap().to_string()
    }

    fn sample_doc() -> XcafDocument {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut leaf = XcafDocNode::with_shape(b.solid.0);
        leaf.attrs = XcafAttrs {
            name: Some("Leaf".into()),
            color: Some((0.0, 1.0, 0.0)),
            layer: Some("L1".into()),
            material: Some("steel".into()),
        };
        leaf.instance_name = Some("leaf-01".into());
        let mut mid = XcafDocNode::new();
        mid.attrs = XcafAttrs::named("Mid");
        mid.children.push(leaf);
        let mut root = XcafDocNode::new();
        root.attrs = XcafAttrs::named("Root");
        root.children.push(mid);
        XcafDocument {
            root,
            version: "1.0".into(),
            views: Vec::new(),
            annotations: Vec::new(),
        }
    }

    #[test]
    fn xcaf_to_bincaf_roundtrip() {
        let doc = sample_doc();
        let bin = xcaf_to_bincaf(&doc);
        let bytes = crate::bincaf::serialize_bincaf(&bin).expect("serialize");
        let bin2 = crate::bincaf::deserialize_bincaf(&bytes).expect("deserialize");
        let got = bincaf_to_xcaf(&bin2);
        assert_eq!(got.root.attrs.name.as_deref(), Some("Root"));
        assert_eq!(got.root.children[0].attrs.name.as_deref(), Some("Mid"));
        let leaf = &got.root.children[0].children[0];
        assert_eq!(leaf.attrs.name.as_deref(), Some("Leaf"));
        assert_eq!(leaf.attrs.color, Some((0.0, 1.0, 0.0)));
        assert_eq!(leaf.attrs.layer.as_deref(), Some("L1"));
        assert_eq!(leaf.attrs.material.as_deref(), Some("steel"));
        assert_eq!(leaf.instance_name.as_deref(), Some("leaf-01"));
        assert_eq!(
            leaf.shape.as_ref().map(|s| s.shape_type()),
            Some(ShapeType::Solid)
        );
    }

    #[test]
    fn attrs_strings_roundtrip() {
        let a = XcafAttrs {
            name: Some("Widget".into()),
            color: Some((1.0, 0.0, 0.0)),
            layer: Some("L1".into()),
            material: Some("steel".into()),
        };
        let s = attrs_to_strings(&a);
        assert_eq!(s.len(), 4);
        let b = strings_to_attrs(&s);
        assert_eq!(a, b);
        // kind:value text form.
        let p = parse_attr_string("color:1,0,0").expect("parse attr");
        assert_eq!((p.kind.as_str(), p.value.as_str()), ("color", "1,0,0"));
        assert_eq!(format_attr_string(&p), "color:1,0,0");
    }

    #[test]
    fn xcaf_xml_roundtrip() {
        let doc = sample_doc();
        let xml = xcaf_to_xml(&doc);
        let got = xml_to_xcaf(&xml).expect("parse xml");
        assert_eq!(got.version, "1.0");
        assert_eq!(got.root.attrs.name.as_deref(), Some("Root"));
        assert_eq!(got.root.children.len(), 1);
        let leaf = &got.root.children[0].children[0];
        assert_eq!(leaf.attrs.name.as_deref(), Some("Leaf"));
        assert_eq!(leaf.attrs.material.as_deref(), Some("steel"));
        assert_eq!(leaf.instance_name.as_deref(), Some("leaf-01"));
        assert_eq!(
            leaf.shape.as_ref().map(|s| s.shape_type()),
            Some(ShapeType::Solid)
        );
    }

    #[test]
    fn find_by_name_nested() {
        let doc = sample_doc();
        assert!(find_by_name(&doc, "Root").is_some());
        assert!(find_by_name(&doc, "Mid").is_some());
        let leaf = find_by_name(&doc, "Leaf").expect("leaf found");
        assert_eq!(leaf.attrs.name.as_deref(), Some("Leaf"));
        assert!(find_by_name(&doc, "Nope").is_none());
        assert_eq!(count_nodes(&doc.root), 3);
        let lines = document_tree_lines(&doc);
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn apply_attrs_flattens() {
        let doc = sample_doc();
        let flat = apply_attrs_to_shapes(&doc);
        assert_eq!(flat.len(), 1, "only the leaf carries a shape");
        assert_eq!(flat[0].1.name.as_deref(), Some("Leaf"));
        assert_eq!(flat[0].1.color, Some((0.0, 1.0, 0.0)));
        assert_eq!(flat[0].0.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&flat[0].0).len(), 6);
    }

    #[test]
    fn read_step_points_lines() {
        let step = "ISO-10303-21;\nDATA;\n\
#1=CARTESIAN_POINT('P1',(1.,2.,3.));\n\
#2=DIRECTION('',(0.,0.,1.));\n\
#3=VECTOR('',#2,5.);\n\
#4=LINE('LineA',#1,#3);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let path = write_temp_step(step);
        let entities = read_step_entities(&path).expect("read entities");
        let lengths = step_curve_lengths(&path).expect("curve lengths");
        std::fs::remove_file(&path).ok();

        let mut n_points = 0;
        let mut n_lines = 0;
        for e in entities {
            match e {
                StepEntity::Point((id, p)) => {
                    assert_eq!(id, "1");
                    assert_eq!(p, [1., 2., 3.]);
                    n_points += 1;
                }
                StepEntity::Line((id, name, origin, dir)) => {
                    assert_eq!(id, "4");
                    assert_eq!(name, "LineA");
                    assert_eq!(origin, [1., 2., 3.]);
                    assert_eq!(dir, [0., 0., 5.]);
                    n_lines += 1;
                }
                _ => {}
            }
        }
        assert_eq!(n_points, 1);
        assert_eq!(n_lines, 1);
        // LINE length = VECTOR magnitude.
        assert_eq!(lengths, vec![("4".to_string(), 5.0)]);
    }

    #[test]
    fn read_step_circle_length() {
        let step = "ISO-10303-21;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n\
#2=DIRECTION('',(0.,0.,1.));\n\
#3=DIRECTION('',(1.,0.,0.));\n\
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5=CIRCLE('',#4,2.);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let path = write_temp_step(step);
        let lengths = step_curve_lengths(&path).expect("lengths");
        std::fs::remove_file(&path).ok();
        assert_eq!(lengths.len(), 1);
        assert_eq!(lengths[0].0, "5");
        assert!(
            (lengths[0].1 - 2.0 * PI * 2.0).abs() < 1e-9,
            "length {}",
            lengths[0].1
        );
    }

    #[test]
    fn read_step_solid_present() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let step = crate::step::write_shape_step(&b.solid.0);
        let path = write_temp_step(&step);
        let entities = read_step_entities(&path).expect("read entities");
        std::fs::remove_file(&path).ok();
        let mut solids = 0;
        for e in &entities {
            if let StepEntity::Solid((_, shape)) = e {
                assert_eq!(shape.shape_type(), ShapeType::Solid);
                assert_eq!(faces_of(shape).len(), 6);
                solids += 1;
            }
        }
        assert_eq!(solids, 1, "one solid expected");
    }

    #[test]
    fn read_step_summary_types() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let step = crate::step::write_shape_step(&b.solid.0);
        let path = write_temp_step(&step);
        let summary = read_step_summary(&path).expect("summary");
        std::fs::remove_file(&path).ok();
        assert!(summary.contains("MANIFOLD_SOLID_BREP: 1"), "{summary}");
        assert!(summary.contains("ADVANCED_FACE"), "{summary}");
        assert!(summary.contains("CARTESIAN_POINT"), "{summary}");
    }

    #[test]
    fn document_to_step_writes() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let s = BRepPrimSphere::make_sphere(0.5);
        let mut part = XcafDocNode::with_shape(b.solid.0);
        part.attrs = XcafAttrs {
            name: Some("CubePart".into()),
            layer: Some("L2".into()),
            ..Default::default()
        };
        part.instance_name = Some("cube-01".into());
        let mut root = XcafDocNode::with_shape(s.solid.0);
        root.attrs = XcafAttrs::named("AssemblyRoot");
        root.children.push(part);
        let doc = XcafDocument {
            root,
            version: "1.0".into(),
            views: Vec::new(),
            annotations: Vec::new(),
        };
        let step = document_to_step(&doc);
        assert!(step.contains("CubePart"), "product name:\n{step}");
        assert!(step.contains("AssemblyRoot"), "assembly name");
        assert!(step.contains("L2"), "layer");
        assert!(step.contains("cube-01"), "instance name");
        assert!(step.contains("XCAFDOC-BEGIN"), "tree block");
    }

    #[test]
    fn document_to_step_roundtrip() {
        let b1 = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let b2 = BRepPrimSphere::make_sphere(1.0);
        let mut part = XcafDocNode::with_shape(b1.solid.0);
        part.attrs = XcafAttrs::named("BoxPart");
        let mut root = XcafDocNode::with_shape(b2.solid.0);
        root.attrs = XcafAttrs::named("Assembly");
        root.children.push(part);
        let doc = XcafDocument {
            root,
            version: "1.0".into(),
            views: Vec::new(),
            annotations: Vec::new(),
        };
        let step = document_to_step(&doc);
        assert!(
            step.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"),
            "assembly:\n{step}"
        );
        let model = crate::step::read_step(&step).expect("read back");
        assert_eq!(model.len(), 2, "two solids");
        let names: Vec<&str> = model.names();
        assert!(names.contains(&"Assembly"), "names: {names:?}");
        assert!(names.contains(&"BoxPart"), "names: {names:?}");
        for ms in &model.shapes {
            assert_eq!(ms.shape.shape_type(), ShapeType::Solid);
        }
    }

    #[test]
    fn empty_doc_ok() {
        let doc = XcafDocument::new();
        let bin = xcaf_to_bincaf(&doc);
        let got = bincaf_to_xcaf(&bin);
        assert!(got.root.attrs.name.is_none());
        assert!(got.root.children.is_empty());
        let xml = xcaf_to_xml(&doc);
        let got2 = xml_to_xcaf(&xml).expect("xml parse");
        assert!(got2.root.children.is_empty());
        let step = document_to_step(&doc);
        assert!(step.contains("DATA;"), "{step}");
        assert!(strings_to_attrs(&[]).is_empty());
    }

    fn view(name: &str, proj: XcafProjection) -> XcafView {
        XcafView {
            name: name.into(),
            eye: GpPnt::new(0.0, -10.0, 5.0),
            target: GpPnt::zero(),
            up: GpDir::new(0.0, 0.0, 1.0).unwrap(),
            projection: proj,
        }
    }

    fn annotation(id: &str, kind: XcafAnnotationKind) -> XcafAnnotation {
        XcafAnnotation {
            id: id.into(),
            kind,
            anchor: GpPnt::new(1.0, 2.0, 3.0),
            text: "width".into(),
            value: 5.0,
        }
    }

    #[test]
    fn views_annotations_roundtrip() {
        let mut doc = XcafDocument::new();
        xcaf_add_view(&mut doc, view("Front", XcafProjection::Perspective));
        xcaf_add_view(&mut doc, view("Top", XcafProjection::Orthographic));
        xcaf_add_annotation(&mut doc, annotation("d1", XcafAnnotationKind::Dimension));
        xcaf_add_annotation(&mut doc, annotation("n1", XcafAnnotationKind::Note));
        assert_eq!(xcaf_list_views(&doc).len(), 2);
        assert_eq!(xcaf_annotations(&doc).len(), 2);
        assert!(xcaf_find_annotation(&doc, "d1").is_some());
        assert!(xcaf_find_annotation(&doc, "missing").is_none());

        // Binary round-trip.
        let bytes = crate::bincaf::serialize_bincaf(&xcaf_to_bincaf(&doc)).expect("serialize");
        let got = bincaf_to_xcaf(&crate::bincaf::deserialize_bincaf(&bytes).expect("deserialize"));
        assert_eq!(got.views.len(), 2);
        assert_eq!(got.views[0].name, "Front");
        assert_eq!(got.views[0].projection, XcafProjection::Perspective);
        assert_eq!(got.views[1].projection, XcafProjection::Orthographic);
        assert_eq!(got.views[0].up, GpDir::new(0.0, 0.0, 1.0).unwrap());
        assert_eq!(got.annotations.len(), 2);
        assert_eq!(got.annotations[0].id, "d1");
        assert_eq!(got.annotations[0].kind, XcafAnnotationKind::Dimension);
        assert_eq!(got.annotations[0].text, "width");
        assert!((got.annotations[0].value - 5.0).abs() < 1e-9);
        assert_eq!(got.annotations[1].kind, XcafAnnotationKind::Note);

        // XML round-trip.
        let got2 = xml_to_xcaf(&xcaf_to_xml(&doc)).expect("xml parse");
        assert_eq!(got2.views.len(), 2);
        assert_eq!(got2.views[1].name, "Top");
        assert_eq!(got2.views[0].eye.distance(&GpPnt::new(0.0, -10.0, 5.0)), 0.0);
        assert_eq!(got2.annotations.len(), 2);
        assert_eq!(got2.annotations[0].id, "d1");

        // Add / remove / find mutation.
        assert!(xcaf_remove_annotation(&mut doc, "d1"));
        assert!(xcaf_find_annotation(&doc, "d1").is_none());
        assert_eq!(xcaf_annotations(&doc).len(), 1);
        assert!(!xcaf_remove_annotation(&mut doc, "missing"));
    }

    #[test]
    fn instance_placement_expansion() {
        // Single placed occurrence through the API.
        let mut doc = XcafDocument::new();
        let mut part = XcafDocNode::with_shape(BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0);
        part.attrs = XcafAttrs::named("Part");
        part.instance_name = Some("part-01".into());
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(0.0, 0.0, 7.0));
        let name = xcaf_add_instance(&mut doc, part, t);
        assert_eq!(name, "part-01");
        assert_eq!(xcaf_instance_count(&doc), 1);
        let inst = xcaf_find_instance(&doc, "part-01").expect("instance found");
        assert_eq!(inst.name, "part-01");
        assert!((inst.transform.translation_part().z - 7.0).abs() < 1e-9);

        // Nested placements: sub-assembly at (0,3,0), leaf at (5,0,0) → (5,3,0).
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut leaf = XcafDocNode::with_shape(b.solid.0);
        leaf.attrs = XcafAttrs::named("Leaf");
        leaf.instance_name = Some("leaf-01".into());
        let mut t1 = GpTrsf::identity();
        t1.set_translation_vec(&GpVec::new(5.0, 0.0, 0.0));
        leaf.location = Some(t1);

        let mut sub = XcafDocNode::new();
        sub.attrs = XcafAttrs::named("Sub");
        sub.instance_name = Some("sub-01".into());
        let mut t2 = GpTrsf::identity();
        t2.set_translation_vec(&GpVec::new(0.0, 3.0, 0.0));
        sub.location = Some(t2);
        sub.children.push(leaf);

        let mut doc2 = XcafDocument::new();
        doc2.root.children.push(sub);
        assert_eq!(xcaf_instance_count(&doc2), 2);

        let sub_inst = xcaf_find_instance(&doc2, "sub-01").expect("sub found");
        assert_eq!(sub_inst.node.children.len(), 1);
        assert!((sub_inst.transform.translation_part().y - 3.0).abs() < 1e-9);

        let expanded = xcaf_expand_instances(&doc2);
        assert_eq!(expanded.len(), 1, "only the leaf carries a shape");
        assert_eq!(expanded[0].0, "leaf-01");
        let shape = &expanded[0].1;
        let mut xs: Vec<f64> = vertices_of(shape)
            .iter()
            .map(|v| BRepTool::vertex_point(v).x())
            .collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((xs[0] - 5.0).abs() < 1e-9, "min x {xs:?}");
        assert!((xs[7] - 6.0).abs() < 1e-9, "max x {xs:?}");
        let ys: Vec<f64> = vertices_of(shape)
            .iter()
            .map(|v| BRepTool::vertex_point(v).y())
            .collect();
        let zs: Vec<f64> = vertices_of(shape)
            .iter()
            .map(|v| BRepTool::vertex_point(v).z())
            .collect();
        let minmax = |v: &[f64]| {
            (
                v.iter().cloned().fold(f64::INFINITY, f64::min),
                v.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            )
        };
        let (min_y, max_y) = minmax(&ys);
        let (min_z, max_z) = minmax(&zs);
        assert!((min_y - 3.0).abs() < 1e-9 && (max_y - 4.0).abs() < 1e-9, "y {ys:?}");
        assert!(min_z.abs() < 1e-9 && (max_z - 1.0).abs() < 1e-9, "z {zs:?}");

        // Locations survive the binary round-trip.
        let bytes = crate::bincaf::serialize_bincaf(&xcaf_to_bincaf(&doc2)).expect("serialize");
        let got = bincaf_to_xcaf(&crate::bincaf::deserialize_bincaf(&bytes).expect("deserialize"));
        let leaf_back = find_by_name(&got, "leaf-01").expect("leaf back");
        let loc = leaf_back.location.as_ref().expect("location preserved");
        assert!((loc.translation_part().x - 5.0).abs() < 1e-9);
        let xml = xml_to_xcaf(&xcaf_to_xml(&doc2)).expect("xml parse");
        let leaf_xml = find_by_name(&xml, "leaf-01").expect("leaf in xml");
        let loc2 = leaf_xml.location.as_ref().expect("xml location preserved");
        assert!((loc2.translation_part().x - 5.0).abs() < 1e-9);
        assert!((loc2.translation_part().y - 0.0).abs() < 1e-9);
    }

    #[test]
    fn step_assembly_with_placements() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let s = BRepPrimSphere::make_sphere(0.5);
        let mut part = XcafDocNode::with_shape(b.solid.0);
        part.attrs = XcafAttrs::named("BoxPart");
        part.instance_name = Some("box-01".into());
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(10.0, 0.0, 0.0));
        part.location = Some(t);
        let mut root = XcafDocNode::with_shape(s.solid.0);
        root.attrs = XcafAttrs::named("Assy");
        root.children.push(part);
        let mut doc = XcafDocument::new();
        doc.root = root;
        xcaf_add_view(&mut doc, view("Front", XcafProjection::Perspective));

        let step = write_step_assembly_with_placements(&doc);
        assert!(step.contains("AXIS2_PLACEMENT_3D"), "placement keyword:\n{step}");
        assert!(
            step.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"),
            "assembly records:\n{step}"
        );
        assert!(step.contains("XCAFDOC-BEGIN"), "metadata block");
        assert!(step.contains("VIEW|"), "view metadata:\n{step}");

        // The box was translated +10 in x: its x=1 corner is at world x=11.
        assert!(step.contains("11."), "transformed coordinate:\n{step}");

        let model = crate::step::read_step(&step).expect("parse back");
        assert_eq!(model.len(), 2, "two solids: {step}");
        let names = model.names();
        assert!(names.contains(&"Assy"), "names: {names:?}");
        assert!(names.contains(&"box-01"), "names: {names:?}");
        let box_shape = model
            .shapes
            .iter()
            .find(|ms| ms.name == "box-01")
            .expect("box shape")
            .shape
            .clone();
        let mut xs: Vec<f64> = vertices_of(&box_shape)
            .iter()
            .map(|v| BRepTool::vertex_point(v).x())
            .collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((xs[0] - 10.0).abs() < 1e-6, "min x {xs:?}");
        assert!((xs[7] - 11.0).abs() < 1e-6, "max x {xs:?}");
    }
