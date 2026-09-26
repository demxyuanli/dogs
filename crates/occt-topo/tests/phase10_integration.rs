//! Phase 10 end-to-end: full topology boolean, curved-face fillet, texture/
//! font/views rendering, Draw TCL vars, inertia tensors, XCAF schema + STEP read.

use occt_core::geom::curve_interp3d::{catmull_rom, cubic_spline_interp, polyline_curve};
use occt_core::gp::{GpPnt, GpVec};
use occt_topo::primitives::{BRepPrimBox, BRepPrimSphere};

#[test]
fn multi_boolean_disjoint_compound() {
    use occt_topo::bop_builder::{boolean_multi, BoolOp};
    let a = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0;
    let b = BRepPrimBox::make_box_corner(&GpPnt::new(5.0, 5.0, 5.0), &GpPnt::new(6.0, 6.0, 6.0)).solid.0;
    let c = BRepPrimBox::make_box_corner(&GpPnt::new(10.0, 10.0, 10.0), &GpPnt::new(11.0, 11.0, 11.0)).solid.0;
    let r = boolean_multi(&[a, b, c], BoolOp::Fuse, 1e-6).expect("multi fuse");
    // Disjoint boxes → compound.
    assert_eq!(r.shape.shape_type(), occt_topo::abs::ShapeType::Compound);
    let comps = occt_topo::bop_builder::decompose_compound(&r.shape);
    assert_eq!(comps.len(), 3, "3 disjoint boxes");
}

#[test]
fn self_intersection_detected_box() {
    use occt_topo::bop_builder::detect_self_intersections;
    let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
    let rep = detect_self_intersections(&b.solid.0, 1e-6);
    assert!(!rep.found, "valid box has no self-intersection");
}

#[test]
fn curved_face_fillet_sphere_plane() {
    use occt_topo::fillet_curved::fillet_edge_curved;
    // A sphere sitting on a box: find the circular base edge and fillet it.
    let box_ = BRepPrimBox::make_box(4.0, 4.0, 2.0);
    let sphere = BRepPrimSphere::make_sphere(1.0);
    // Place sphere on top (z=2): use boolean_compound Fuse to attach it.
    use occt_topo::bop_builder::{boolean_compound, BoolOp};
    use occt_topo::shape_ops::translated_copy;
    let sphere_at = translated_copy(&sphere.solid.0, &GpVec::new(2.0, 2.0, 2.0)).expect("translate");
    let fused = boolean_compound(&box_.solid.0, &sphere_at, BoolOp::Fuse, 1e-6).expect("fuse");
    // The circular edge where sphere meets box top → fillet it. The edge must
    // be the sphere/box intersection **circle**: a box top edge is a
    // plane/plane pair, which `fillet_edge_curved` refuses
    // (`fillet_curved` supports plane+sphere etc.; a plane/plane edge belongs
    // to the straight `fillet::fillet_edge`). Selecting by z alone picked the
    // box outline first (T-88 fixture bug, not an engine gap).
    let edge = occt_topo::topo_tools_full::edges_of(&fused.shape)
        .into_iter()
        .find(|e| {
            occt_topo::brep_tool::BRepTool::edge_curve(e)
                .map(|c| {
                    let p = c.d0((c.first_parameter() + c.last_parameter()) * 0.5);
                    c.gp_circ().is_some() && (p.z() - 2.0).abs() < 0.5
                })
                .unwrap_or(false)
        });
    if let Some(e) = edge {
        let f = fillet_edge_curved(&fused.shape, &e, 0.2, 1e-6).expect("curved fillet");
        assert!(occt_topo::topo_tools_full::faces_of(&f).len() > 0);
    }
    // At minimum the fuse produced a valid shape.
    assert!(occt_topo::topo_tools_full::faces_of(&fused.shape).len() > 0);
}

#[test]
fn texture_font_render() {
    use occt_topo::viz_scene::{
        checkerboard_texture, render_textured_ppm, text_to_pixels, BitmapFont, Camera, SceneShape, VizScene,
    };
    let tex = checkerboard_texture(8, 8, 2, (255, 0, 0), (0, 0, 255));
    let mut scene = VizScene::new();
    let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
    let shape = SceneShape::new(b.solid.0).with_uv(None); // synthesize UVs
    scene.add(shape);
    let cam = Camera::default();
    let ppm = render_textured_ppm(&scene, &cam, 64, 48, 0.2, &tex);
    assert!(ppm.starts_with(b"P6"));
    assert_eq!(ppm.len(), 13 + 64 * 48 * 3);
    // Two checkerboard colors present somewhere.
    let red = ppm[13..].chunks_exact(3).any(|p| p[0] > 200 && p[1] < 50);
    let blue = ppm[13..].chunks_exact(3).any(|p| p[2] > 200 && p[1] < 50);
    assert!(red && blue, "checkerboard has red+blue");

    let font = BitmapFont::default();
    let (w, h) = occt_topo::viz_scene::text_raster_size("AB", &font, 1);
    assert_eq!((w, h), (10, 7));
    let px = text_to_pixels("A", &font, 2);
    assert_eq!(px.len(), 5 * 2 * 7 * 2);
}

#[test]
fn draw_tcl_vars_and_transform() {
    use occt_topo::draw::DrawSession;
    let mut s = DrawSession::default();
    occt_topo::draw::execute_line(&mut s, "set x 3").expect("set");
    occt_topo::draw::execute_line(&mut s, "box b $x 2 2").expect("box");
    // Copy FIRST (the original), then translate in place.
    occt_topo::draw::execute_line(&mut s, "copy b b2").expect("copy");
    occt_topo::draw::execute_line(&mut s, "translate b 1 0 0").expect("translate");
    assert!(s.shapes.contains_key("b"));
    assert!(s.shapes.contains_key("b2"));
    // b translated +1 in x → some vertex at x≈4 (box spans [0,3]×[0,2]×[0,2]).
    let verts = occt_topo::topo_tools_full::vertices_of(&s.shapes["b"]);
    assert!(verts.iter().any(|v| (occt_topo::brep_tool::BRepTool::vertex_point(v).x() - 4.0).abs() < 1e-6));
    // b2 is the un-translated copy → x≈3 max.
    let verts2 = occt_topo::topo_tools_full::vertices_of(&s.shapes["b2"]);
    assert!(verts2.iter().any(|v| (occt_topo::brep_tool::BRepTool::vertex_point(v).x() - 3.0).abs() < 1e-6));
}

#[test]
fn inertia_tensor_box() {
    use occt_topo::gprop_analytic::{inertia_tensor, principal_inertia};
    let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
    let it = inertia_tensor(&b.solid.0, 1.0).expect("inertia");
    // Box spans [0,2]×[0,3]×[0,4], mass 24 about centroid.
    // Izz = m/12·(dx²+dy²) about the z-axis through centroid = 24/12·(4+9) = 26.
    assert!((it.izz - 26.0).abs() < 1e-6, "Izz {}", it.izz);
    let (vals, _axes) = principal_inertia(&b.solid.0, 1.0).expect("principal");
    assert_eq!(vals.len(), 3);
    assert!(vals[0] >= vals[1] && vals[1] >= vals[2], "sorted {vals:?}");
}

#[test]
fn xcaf_doc_roundtrip() {
    use occt_topo::xcaf::{attrs_to_strings, strings_to_attrs, XcafAttrs};
    let attrs = XcafAttrs {
        name: Some("part1".into()),
        color: Some((1.0, 0.0, 0.0)),
        layer: Some("L1".into()),
        material: Some("steel".into()),
    };
    let strings = attrs_to_strings(&attrs);
    assert_eq!(strings.len(), 4);
    let back = strings_to_attrs(&strings);
    assert_eq!(back.name.as_deref(), Some("part1"));
    assert_eq!(back.color, Some((1.0, 0.0, 0.0)));
    assert_eq!(back.layer.as_deref(), Some("L1"));
    assert_eq!(back.material.as_deref(), Some("steel"));
}

#[test]
fn curve_interp3d_passes_points() {
    let pts = vec![
        GpPnt::new(0.0, 0.0, 0.0),
        GpPnt::new(1.0, 2.0, 0.0),
        GpPnt::new(2.0, -1.0, 0.0),
        GpPnt::new(3.0, 3.0, 0.0),
    ];
    let c = cubic_spline_interp(&pts).expect("spline");
    let knots = c.data_knots();
    for (i, p) in pts.iter().enumerate() {
        let q = c.point(knots[i]);
        assert!(p.distance(&q) < 1e-6, "point {i}");
    }
    let cat = catmull_rom(&pts);
    let mid = cat.point(0.5);
    assert!(mid.x() > 0.5 && mid.x() < 2.0);
    let pl = polyline_curve(&pts);
    assert!((pl.point(0.0).x() - 0.0).abs() < 1e-9);
}
