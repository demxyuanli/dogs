use super::prelude::*;
use super::*;
    use crate::brep_surface::SurfaceKind;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::topo_tools_full::{faces_of, vertices_of};

    #[test]
    fn box_roundtrip_counts_points_and_faces() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let step = write_shape_step(&b.solid.0);
        let model = read_step(&step).expect("parse step");
        assert_eq!(model.len(), 1);

        let (nv, ne, nf) = step_roundtrip_counts(&b.solid.0).expect("roundtrip counts");
        assert_eq!((nv, ne, nf), (8, 12, 6));

        let shape = &model.shapes[0].shape;
        let verts = vertices_of(shape);
        assert_eq!(verts.len(), 8);
        let pts: Vec<GpPnt> = verts
            .iter()
            .map(|v| GeometryRegistry::global().vertex_point(&v.0))
            .collect();
        for &(x, y, z) in &[
            (0., 0., 0.),
            (2., 0., 0.),
            (2., 3., 0.),
            (0., 3., 0.),
            (0., 0., 4.),
            (2., 0., 4.),
            (2., 3., 4.),
            (0., 3., 4.),
        ] {
            assert!(
                pts.iter().any(|p| p.is_equal(&GpPnt::new(x, y, z))),
                "missing corner ({x},{y},{z})"
            );
        }

        let faces = faces_of(shape);
        assert_eq!(faces.len(), 6);
        for f in &faces {
            let surf = GeometryRegistry::global().face_surface(&f.0);
            assert!(surf.is_some(), "face has a surface");
            if let Some(s) = surf {
                assert_eq!(classify_surface(s.as_ref()), SurfaceKind::Plane);
            }
        }
    }

    #[test]
    fn sphere_roundtrip_surface() {
        let s = BRepPrimSphere::make_sphere(2.5);
        let step = write_shape_step(&s.solid.0);
        let model = read_step(&step).expect("parse sphere step");
        assert_eq!(model.len(), 1);
        let shape = &model.shapes[0].shape;
        let faces = faces_of(shape);
        assert_eq!(faces.len(), 1);
        let surf = GeometryRegistry::global()
            .face_surface(&faces[0].0)
            .expect("sphere face has surface");
        assert_eq!(classify_surface(surf.as_ref()), SurfaceKind::Sphere);
    }

    #[test]
    fn model_two_shapes_names_preserved() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let s = BRepPrimSphere::make_sphere(1.0);
        let mut model = BRepModel::new();
        model.add("Box", b.solid.0.clone());
        model.add("Sphere", s.solid.0.clone());
        let step = write_step(&model);
        let got = read_step(&step).expect("parse model step");
        assert_eq!(got.len(), 2);
        let mut names: Vec<&str> = got.names();
        names.sort();
        assert_eq!(names, vec!["Box", "Sphere"]);
    }

    #[test]
    fn malformed_missing_end_iso_errors() {
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n#1=DIRECTION('',(1.,0.,0.));\nENDSEC;\n";
        assert!(read_step(s).is_err());
    }

    #[test]
    fn reader_skips_unknown_entity_with_warning() {
        // A minimal file whose representation references an unsupported entity;
        // the reader must collect a warning and skip the shape gracefully.
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n\
#1=APPLICATION_CONTEXT('AUTOMOTIVE_DESIGN');\n\
#2=WIDGET_FROBNICATOR('',42.);\n\
#3=GEOMETRIC_REPRESENTATION_CONTEXT('','',3);\n\
#4=ADVANCED_BREP_SHAPE_REPRESENTATION('widget',(#2),#3);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let (model, warnings) = read_step_with_warnings(s).expect("parse");
        assert!(model.is_empty(), "unsupported shape is skipped");
        assert!(
            warnings.iter().any(|w| w.contains("WIDGET_FROBNICATOR")),
            "expected a warning about the unknown entity, got {warnings:?}"
        );
    }

    #[test]
    fn step_real_formats_decimal_point() {
        assert_eq!(step_real(0.0), "0.0");
        assert_eq!(step_real(2.5), "2.5");
        assert_eq!(step_real(-3.0), "-3.0");
        assert!(step_real(1e20).contains('.'));
    }

    #[test]
    fn split_top_handles_nested_lists() {
        let args = split_top("'',(1.,0.,0.),(#2,#3)");
        assert_eq!(args.len(), 3);
        assert_eq!(args[0], "''");
        assert_eq!(args[1], "(1.,0.,0.)");
        assert_eq!(args[2], "(#2,#3)");
    }

    #[test]
    fn cylinder_roundtrip() {
        use crate::primitives::BRepPrimCylinder;
        let c = BRepPrimCylinder::make_cylinder(2.0, 10.0);
        let (nv, ne, nf) = step_roundtrip_counts(&c.solid.0).expect("counts");
        assert_eq!((nv, ne, nf), (2, 3, 3));
        let s = write_shape_step(&c.solid.0);
        let m = read_step(&s).expect("read cylinder");
        assert_eq!(m.len(), 1);
        let fs = faces_of(&m.shapes[0].shape);
        assert_eq!(fs.len(), 3);
        // Two cap faces should classify as planes.
        let planar = fs
            .iter()
            .filter(|f| {
                GeometryRegistry::global()
                    .face_surface(&f.0)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Plane)
                    .unwrap_or(false)
            })
            .count();
        assert_eq!(planar, 2);
    }

    #[test]
    fn step_file_io_roundtrip() {
        use crate::primitives::BRepPrimSphere;
        let s = BRepPrimSphere::make_sphere(1.5);
        let mut model = BRepModel::new();
        model.add("Ball", s.solid.0.clone());
        let path = { let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target").join("test-out"); let _ = std::fs::create_dir_all(&d); d }.join("occt_step_test.step");
        let p = path.to_str().unwrap();
        write_step_file(p, &model).expect("write file");
        let got = read_step_file(p).expect("read file");
        assert_eq!(got.len(), 1);
        assert_eq!(got.names(), vec!["Ball"]);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn bspline_curve_written() {
        let c = GeomBSplineCurve::new(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(2.0, 1.0, 0.0),
                GpPnt::new(3.0, 0.0, 0.0),
            ],
            vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let mut w = StepWriter::new();
        let id = write_bspline_curve(&mut w, &c).unwrap();
        let line = &w.lines[id - 1];
        assert!(line.contains("B_SPLINE_CURVE_WITH_KNOTS"), "{line}");
        // Degree 2, four control-point refs, polynomial (SELF weights).
        assert!(line.contains("('',2,(#1,#2,#3,#4)"), "{line}");
        assert!(line.contains("SELF"), "{line}");
        // Knots 0/0.5/1 with multiplicities 3/1/3, and the UNSPECIFIED defaults.
        assert!(line.contains("(0.0,0.5,1.0)"), "{line}");
        assert!(line.contains("(3,1,3)"), "{line}");
        assert!(line.contains("UNSPECIFIED"), "{line}");
    }

    #[test]
    fn rational_bspline_curve_written() {
        let c = GeomBSplineCurve::rational(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(2.0, 0.0, 0.0),
            ],
            vec![1.0, 0.5, 1.0],
            vec![0.0, 0.0, 0.5, 1.0, 1.0],
            1,
        )
        .unwrap();
        let mut w = StepWriter::new();
        let id = write_bspline_curve(&mut w, &c).unwrap();
        let line = &w.lines[id - 1];
        assert!(line.contains("(1.0,0.5,1.0)"), "{line}");
        assert!(!line.contains("SELF"), "{line}");
    }

    #[test]
    fn bspline_surface_written() {
        let (ku, kv) = occt_geom::bspline_surface::bspline_surface_uniform_knots(3, 3, 2, 2);
        let poles = vec![
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
                GpPnt::new(0.0, 2.0, 0.0),
            ],
            vec![
                GpPnt::new(1.0, 0.0, 1.0),
                GpPnt::new(1.0, 1.0, 1.0),
                GpPnt::new(1.0, 2.0, 1.0),
            ],
            vec![
                GpPnt::new(2.0, 0.0, 0.0),
                GpPnt::new(2.0, 1.0, 0.0),
                GpPnt::new(2.0, 2.0, 0.0),
            ],
        ];
        let s = GeomBSplineSurface::new(poles, ku, kv, 2, 2).unwrap();
        let mut w = StepWriter::new();
        let id = write_bspline_surface(&mut w, &s).unwrap();
        let line = &w.lines[id - 1];
        assert!(line.contains("B_SPLINE_SURFACE_WITH_KNOTS"), "{line}");
        // u_degree, v_degree both 2; 3×3 pole grid; polynomial.
        assert!(line.contains("('',2,2,((#1,#2,#3),(#4,#5,#6),(#7,#8,#9))"), "{line}");
        // `GeomToStep_MakeBSplineSurfaceWithKnots.cxx:150-163`: the attribute
        // order is surface form, U/V closed, self-intersect, U/V multiplicities,
        // U/V knots, knot spec; the multiplicities precede the knot values, and
        // a polynomial surface has no weights attribute (the grid appears only
        // in the rational complex form, `...AndRational....cxx:170-176`). Same
        // record layout as OCCT's own output, `data/occ/T0M.stp:901`.
        assert!(
            line.contains(",.UNSPECIFIED.,.F.,.F.,.F.,(3,3),(3,3),(0.0,1.0),(0.0,1.0),.UNSPECIFIED."),
            "{line}"
        );
        assert!(!line.contains("SELF"), "{line}");
    }

    #[test]
    fn trimmed_curve_written() {
        let mut w = StepWriter::new();
        let id = write_trimmed_curve(&mut w, 42, 1.5, 0.5).unwrap();
        let line = &w.lines[id - 1];
        // Bounds are normalised to ascending order.
        assert!(
            line.contains("TRIMMED_CURVE('',#42,1,0.5,1.5,PARAMETER)"),
            "{line}"
        );
    }

    #[test]
    fn offset_curve_written() {
        let mut w = StepWriter::new();
        let id = write_offset_curve(&mut w, 7, 2.5, 9).unwrap();
        let line = &w.lines[id - 1];
        assert!(
            line.contains("OFFSET_CURVE_3D('',#7,#9,2.5,UNSPECIFIED,UNSPECIFIED)"),
            "{line}"
        );
    }

    #[test]
    fn circle_ellipse_params() {
        let ax2 = GpAx2::new(GpPnt::zero(), dir_z(), dir_x()).unwrap();
        let c = GeomCircle::new(GpCirc::new(ax2.clone(), 2.0));
        let mut w = StepWriter::new();
        let id = write_conic_params(&mut w, &c).unwrap().expect("circle is a conic");
        let line = &w.lines[id - 1];
        assert!(line.contains("CIRCLE"), "{line}");
        assert!(line.contains(",2.0)"), "{line}");

        let e = GeomEllipse::new(GpElips::new(ax2, 3.0, 1.5));
        let id = write_conic_params(&mut w, &e).unwrap().expect("ellipse is a conic");
        let line = &w.lines[id - 1];
        assert!(line.contains("ELLIPSE"), "{line}");
        assert!(line.contains("3.0,1.5)"), "{line}");
    }

    /// A small curved B-spline surface grid (z = u² + v³) reused by the
    /// spline round-trip tests.
    fn bspline_test_surface() -> GeomBSplineSurface {
        let (nu, nv) = (5, 5);
        let points: Vec<Vec<GpPnt>> = (0..nu)
            .map(|i| {
                (0..nv)
                    .map(|j| {
                        let u = i as f64 / (nu - 1) as f64;
                        let v = j as f64 / (nv - 1) as f64;
                        GpPnt::new(u, v, u * u + v * v * v)
                    })
                    .collect()
            })
            .collect();
        occt_geom::bspline_surface::fit_surface_grid(&points, 3, 3).unwrap()
    }

    #[test]
    fn step_bspline_roundtrip() {
        let surf = bspline_test_surface();
        let seen = surf.poles.len();
        let face = crate::brep_builder_full::BRepBuilderFace::from_surface(Arc::new(surf));
        let step = write_step_with_splines(&face).unwrap();
        let m = read_step(&step).expect("read bspline step");
        assert_eq!(m.len(), 1);
        assert!(m.shapes[0].shape.is_face(), "shape type preserved");
        let fs = faces_of(&m.shapes[0].shape);
        assert!(!fs.is_empty(), "reconstructed shape has faces");
        let back = GeometryRegistry::global()
            .face_surface(&fs[0].0)
            .and_then(|s| s.osculating_bspline())
            .expect("surface reads back as a B-spline");
        assert_eq!(back.poles.len(), seen, "pole grid preserved");
        assert!(!back.is_rational(), "polynomial surface carries no weights");
        // Rational arm: the merged complex entity must read back with weights.
        let r = bspline_test_surface();
        let w: Vec<Vec<f64>> = r
            .poles
            .iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .map(|(j, _)| 1.0 + 0.1 * (i + j) as f64)
                    .collect()
            })
            .collect();
        let rr = occt_geom::bspline_surface::GeomBSplineSurface::rational(
            r.poles.clone(),
            w,
            r.knots_u.clone(),
            r.knots_v.clone(),
            r.deg_u,
            r.deg_v,
        )
        .unwrap();
        let rface = crate::brep_builder_full::BRepBuilderFace::from_surface(Arc::new(rr));
        let rstep = write_step_with_splines(&rface).unwrap();
        assert!(
            rstep.contains("RATIONAL_B_SPLINE_SURFACE"),
            "rational surface uses the complex form"
        );
        let rm = read_step(&rstep).expect("read rational bspline step");
        let rfs = faces_of(&rm.shapes[0].shape);
        assert!(!rfs.is_empty(), "rational shape has faces");
        let rback = GeometryRegistry::global()
            .face_surface(&rfs[0].0)
            .and_then(|s| s.osculating_bspline())
            .expect("rational surface reads back as a B-spline");
        assert!(rback.is_rational(), "weights preserved through the complex form");
        assert_eq!(
            rback.weights.as_ref().map_or(0, |w| w.len()),
            rback.poles.len(),
            "weight grid matches the pole grid"
        );
    }

    #[test]
    fn step_spline_reader_handles() {
        let surf = bspline_test_surface();
        let face = crate::brep_builder_full::BRepBuilderFace::from_surface(Arc::new(surf));
        let step = write_step_with_splines(&face).unwrap();
        let path = { let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target").join("test-out"); let _ = std::fs::create_dir_all(&d); d }.join("occt_step_spline_reader.step");
        let p = path.to_str().unwrap();
        std::fs::write(p, &step).expect("write spline step");
        let m = read_step_file(p).expect("read spline step file");
        assert_eq!(m.len(), 1);
        let fs = faces_of(&m.shapes[0].shape);
        assert!(!fs.is_empty());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn box_roundtrip_splines() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let step = write_step_with_splines(&b.solid.0).unwrap();
        let m = read_step(&step).expect("read spline box step");
        assert_eq!(m.len(), 1);
        let c = crate::topo_tools_full::shape_counts(&m.shapes[0].shape);
        assert_eq!(c.get(&ShapeType::Vertex).copied().unwrap_or(0), 8);
        assert_eq!(c.get(&ShapeType::Edge).copied().unwrap_or(0), 12);
        assert_eq!(c.get(&ShapeType::Face).copied().unwrap_or(0), 6);
    }

    #[test]
    fn conic_circle_roundtrip() {
        let b = TopoBuilder::new();
        let ax2 = GpAx2::new(GpPnt::zero(), dir_z(), dir_x()).unwrap();
        let circ = GpCirc::new(ax2, 2.0);
        let mut e = b.make_edge(Arc::new(GeomCircle::new(circ)), 0.0, 2.0 * PI);
        let v = b.make_vertex(GpPnt::new(2.0, 0.0, 0.0), 0.0);
        b.add_edge_vertices(&mut e, &v, &v);
        let w = b.make_wire(&[Edge(e.0)]);
        let ax3 = GpAx3::new(GpPnt::zero(), dir_z(), &dir_x()).unwrap();
        let face = crate::brep_builder_full::BRepBuilderFace::from_wire(&w, &GpPln::new(ax3)).unwrap();
        let step = write_step_with_splines(&face).unwrap();
        let m = read_step(&step).expect("read circle face step");
        assert_eq!(m.len(), 1);
        let fs = faces_of(&m.shapes[0].shape);
        assert!(!fs.is_empty(), "reconstructed circle face");
    }

    #[test]
    fn step_assembly_written() {
        let a = StepAssembly {
            name: "Assy".into(),
            products: vec![
                ("PartA".into(), BRepPrimBox::make_box(1.0, 2.0, 3.0).solid.0),
                ("PartB".into(), BRepPrimBox::make_box(2.0, 1.0, 1.0).solid.0),
            ],
            children: vec![("Assy".into(), vec!["PartA".into(), "PartB".into()])],
        };
        let out = write_step_assembly(&a).unwrap();
        assert!(out.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"), "{out}");
        // The assembly root plus the two parts.
        assert!(out.matches("PRODUCT('").count() >= 2, "{out}");
        // The usage records reference the correct product definitions.
        assert!(out.contains("PartA"), "{out}");
        assert!(out.contains("PartB"), "{out}");
    }

    #[test]
    fn step_assembly_unknown_child_errors() {
        let a = StepAssembly {
            name: "Assy".into(),
            products: vec![("PartA".into(), BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0)],
            children: vec![("Assy".into(), vec!["Ghost".into()])],
        };
        let err = write_step_assembly(&a).unwrap_err();
        assert!(err.contains("Ghost"), "{err}");
    }

    #[test]
    fn step_color_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let out = write_step_with_color(&b.solid.0, (0.8, 0.2, 0.1)).unwrap();
        assert!(out.contains("COLOUR_RGB"), "{out}");
        assert!(out.contains("0.8"), "{out}");
        assert!(out.contains("0.2"), "{out}");
        assert!(out.contains("0.1"), "{out}");
        // The style chain is attached to the shape representation.
        assert!(out.contains("SURFACE_STYLE_FILL_AREA"), "{out}");
        assert!(out.contains("SURFACE_STYLE_USAGE"), "{out}");
        // The file still reads back as the shape.
        let m = read_step(&out).expect("read colored step");
        assert_eq!(m.len(), 1);
    }

    #[test]
    fn step_name_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let out = write_step_with_name(&b.solid.0, "MyBox").unwrap();
        assert!(out.contains("MyBox"), "{out}");
        // The PRODUCT name survives a read (representation name is derived
        // from the product name by the writer).
        let m = read_step(&out).expect("read named step");
        assert_eq!(m.names(), vec!["MyBox"]);
    }

    #[test]
    fn step_units_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let units = StepUnits {
            length_unit_m: 0.001,
            angle_unit_rad: 1.0,
        };
        let out = write_step_with_units(&b.solid.0, "Metric", &units).unwrap();
        assert!(out.contains("Metric"), "{out}");
        assert!(out.contains("SI_UNIT"), "{out}");
        assert!(out.contains("DIMENSIONAL_EXPONENTS"), "{out}");
        // The file still reads back as the shape (unit records are skipped).
        let m = read_step(&out).expect("read units step");
        assert_eq!(m.len(), 1);
    }

    #[test]
    fn step_name_color_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let out = write_step_with_name_and_color(&b.solid.0, "RedBox", (1.0, 0.0, 0.0)).unwrap();
        assert!(out.contains("RedBox"), "{out}");
        assert!(out.contains("COLOUR_RGB"), "{out}");
        assert!(out.contains("SURFACE_STYLE_USAGE"), "{out}");
        let m = read_step(&out).expect("read named+color step");
        assert_eq!(m.names(), vec!["RedBox"]);
    }

    #[test]
    fn step_assembly_roundtrip() {
        let a = StepAssembly {
            name: "Rig".into(),
            products: vec![
                ("Leg".into(), BRepPrimBox::make_box(0.5, 0.5, 2.0).solid.0),
                ("Foot".into(), BRepPrimBox::make_box(0.8, 0.3, 0.2).solid.0),
            ],
            children: vec![
                ("Rig".into(), vec!["Leg".into()]),
                ("Leg".into(), vec!["Foot".into()]),
            ],
        };
        let out = write_step_assembly(&a).unwrap();
        let got = read_step_assembly(&out).expect("read assembly step");
        // The tree is reconstructed: root "Rig", two parts, two usage edges.
        assert_eq!(got.name, "Rig");
        assert_eq!(got.products.len(), 2);
        let names: Vec<String> = got.products.iter().map(|(n, _)| n.clone()).collect();
        assert!(names.contains(&"Leg".into()), "{names:?}");
        assert!(names.contains(&"Foot".into()), "{names:?}");
        let flat: Vec<(String, String)> = got
            .children
            .iter()
            .flat_map(|(p, kids)| kids.iter().map(move |k| (p.clone(), k.clone())))
            .collect();
        assert!(flat.contains(&("Rig".into(), "Leg".into())), "{flat:?}");
        assert!(flat.contains(&("Leg".into(), "Foot".into())), "{flat:?}");
    }

    #[test]
    fn step_spline_reader_plain_bspline() {
        // A hand-written STEP file using the plain (knotless) B_SPLINE_CURVE
        // and B_SPLINE_SURFACE entities must read back into a valid face.
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n\
#1=APPLICATION_CONTEXT('AUTOMOTIVE_DESIGN');\n\
#2=CARTESIAN_POINT('',(0.,0.,0.));\n\
#3=CARTESIAN_POINT('',(0.,1.,0.));\n\
#4=CARTESIAN_POINT('',(1.,0.,1.));\n\
#5=CARTESIAN_POINT('',(1.,1.,1.));\n\
#6=CARTESIAN_POINT('',(2.,0.,0.));\n\
#7=CARTESIAN_POINT('',(2.,1.,0.));\n\
#8=CARTESIAN_POINT('',(3.,0.,1.));\n\
#9=CARTESIAN_POINT('',(3.,1.,1.));\n\
#10=B_SPLINE_SURFACE('',1,1,((#2,#3),(#4,#5),(#6,#7),(#8,#9)),UNSPECIFIED,.F.,.F.,.F.);\n\
#11=GEOMETRIC_REPRESENTATION_CONTEXT('','',3);\n\
#12=ADVANCED_FACE('',#10,(),.T.);\n\
#13=ADVANCED_BREP_SHAPE_REPRESENTATION('plain',(#12),#11);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let m = read_step(s).expect("read plain bspline step");
        assert_eq!(m.len(), 1);
        let fs = faces_of(&m.shapes[0].shape);
        assert!(!fs.is_empty(), "plain B_SPLINE_SURFACE yields a face");
    }

    #[test]
    fn polyline_written_and_read() {
        let pts = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(2.0, 1.0, 0.0),
        ];
        let mut w = StepWriter::new();
        let id = write_polyline(&mut w, &pts).unwrap();
        assert!(w.lines[id - 1].contains("POLYLINE('',(#1,#2,#3))"), "{}", w.lines[id - 1]);

        // A wire built from a polyline edge round-trips through a hand-written
        // STEP file that references the POLYLINE.
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n\
#1=APPLICATION_CONTEXT('AUTOMOTIVE_DESIGN');\n\
#2=CARTESIAN_POINT('',(0.,0.,0.));\n\
#3=CARTESIAN_POINT('',(1.,0.,0.));\n\
#4=CARTESIAN_POINT('',(2.,1.,0.));\n\
#5=POLYLINE('',(#2,#3,#4));\n\
#6=DIRECTION('',(1.,0.,0.));\n\
#7=VECTOR('',#6,1.);\n\
#8=LINE('',#2,#7);\n\
#9=VERTEX_POINT('',#2);\n\
#10=VERTEX_POINT('',#4);\n\
#11=EDGE_CURVE('',#9,#10,#5,.T.);\n\
#12=ORIENTED_EDGE('',#9,#10,#11,.T.);\n\
#13=EDGE_LOOP('',(#12));\n\
#14=GEOMETRIC_REPRESENTATION_CONTEXT('','',3);\n\
#15=ADVANCED_BREP_SHAPE_REPRESENTATION('poly',(#11),#14);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let m = read_step(s).expect("read polyline step");
        assert_eq!(m.len(), 1);
        let e = crate::topo_tools_full::edges_of(&m.shapes[0].shape);
        assert_eq!(e.len(), 1);
    }

    #[test]
    fn step_options_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let opts = StepWriteOptions {
            name: "OptBox".into(),
            color: Some((0.25, 0.5, 0.75)),
            units: StepUnits {
                length_unit_m: 0.001,
                angle_unit_rad: 1.0,
            },
            splines: true,
        };
        let out = write_step_with_options(&b.solid.0, &opts).unwrap();
        assert!(out.contains("OptBox"), "{out}");
        assert!(out.contains("COLOUR_RGB"), "{out}");
        assert!(out.contains("0.25"), "{out}");
        assert!(out.contains("SI_UNIT"), "{out}");
        let m = read_step(&out).expect("read options step");
        assert_eq!(m.names(), vec!["OptBox"]);
    }

    #[test]
    fn step_header_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let h = StepHeader {
            description: "Test part".into(),
            name: "part.step".into(),
            timestamp: "2030-01-01T00:00:00".into(),
            author: "alice".into(),
            organization: "acme".into(),
            preprocessor: "occt-topo".into(),
            originator: "bob".into(),
            schema: "AP242".into(),
        };
        let out = write_step_with_header(&b.solid.0, "Part", &h).unwrap();
        assert!(out.contains("Test part"), "{out}");
        assert!(out.contains("part.step"), "{out}");
        assert!(out.contains("2030-01-01T00:00:00"), "{out}");
        assert!(out.contains("AP242"), "{out}");
        let m = read_step(&out).expect("read header step");
        assert_eq!(m.names(), vec!["Part"]);
    }

    #[test]
    fn step_assembly_with_splines_written() {
        let a = StepAssembly {
            name: "Assy".into(),
            products: vec![("PartA".into(), BRepPrimBox::make_box(1.0, 2.0, 3.0).solid.0)],
            children: vec![("Assy".into(), vec!["PartA".into()])],
        };
        let out = write_step_assembly_with_splines(&a).unwrap();
        assert!(out.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"), "{out}");
        let got = read_step_assembly(&out).expect("read assembly step");
        assert_eq!(got.name, "Assy");
        assert_eq!(got.products.len(), 1);
    }
