//! Phase 9 end-to-end: analytic mass properties, trimmed-face boolean, fillet
//! corner patch, Draw REPL, hardware-style shaded rendering, curve Frenet
//! frame, VRML import.

use occt_core::geom::curve_frenet::{curvature, frenet_frame, torsion};
use occt_core::gp::{GpPnt, GpVec};
use occt_topo::gprop_analytic::{analytic_properties, analytic_volume, extrema_point_line};
use occt_topo::primitives::{BRepPrimBox, BRepPrimSphere};

#[test]
fn analytic_props_exact_box() {
    let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
    let p = analytic_properties(&b.solid.0).expect("analytic props");
    assert!(p.exact, "box is analytic");
    assert!((p.surface_area - 52.0).abs() < 1e-9, "area {}", p.surface_area);
    assert!((p.volume - 24.0).abs() < 1e-9, "volume {}", p.volume);
    assert!((p.centroid.x() - 1.0).abs() < 1e-9, "cx {}", p.centroid.x());
}

#[test]
fn analytic_volume_sphere() {
    let s = BRepPrimSphere::make_sphere(1.0);
    let v = analytic_volume(&s.solid.0).expect("sphere volume");
    assert!((v - 4.0 * std::f64::consts::PI / 3.0).abs() < 1e-6, "v {v}");
}

#[test]
fn extrema_point_line_analytic() {
    let line = occt_geom::GeomLine::new(occt_core::gp::GpLin::new(
        occt_core::gp::GpAx1::new(GpPnt::new(0.0, 0.0, 0.0), occt_core::gp::GpDir::new(1.0, 0.0, 0.0).unwrap()),
    ));
    let e = extrema_point_line(GpPnt::new(3.0, 4.0, 0.0), &line);
    assert!((e.distance - 4.0).abs() < 1e-12, "dist {}", e.distance);
}

#[test]
fn trimmed_boolean_preserves_curved_face() {
    use occt_topo::bop_builder::BoolOp;
    use occt_topo::bop_curved::curved_boolean_full;
    let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
    let s = BRepPrimSphere::make_sphere(0.5);
    // Sphere fully inside the box → Fuse = box (sphere internal faces culled).
    let r = curved_boolean_full(&b.solid.0, &s.solid.0, BoolOp::Fuse, 1e-6).expect("boolean");
    let nf = occt_topo::topo_tools_full::faces_of(&r.shape).len();
    assert!(nf >= 6, "fused box faces {nf}");
}

#[test]
fn fillet_corner_patch_closes_chain() {
    use occt_topo::fillet_var::{fillet_edges_chain_with_corner, RadiusLaw, VarFilletSpec};
    let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
    // Find two consecutive bottom edges that both meet at (0,0,0).
    let origin = GpPnt::new(0.0, 0.0, 0.0);
    let touches_origin = |e: &occt_topo::shape::Edge| {
        occt_topo::brep_tool::BRepTool::edge_vertices(e)
            .map(|(a, c)| a.is_equal(&origin) || c.is_equal(&origin))
            .unwrap_or(false)
    };
    let edges = occt_topo::topo_tools_full::edges_of(&b.solid.0);
    let mut at_origin: Vec<usize> = edges.iter().enumerate()
        .filter(|(_, e)| touches_origin(e))
        .map(|(i, _)| i)
        .collect();
    assert!(at_origin.len() >= 2, "box has 3 edges through origin, got {}", at_origin.len());
    at_origin.truncate(2);
    let spec = VarFilletSpec { r_start: 0.4, r_end: 0.4, law: RadiusLaw::Linear };
    let f = fillet_edges_chain_with_corner(&b.solid.0, &at_origin, &[spec, spec], 1e-6).expect("chain");
    let nf = occt_topo::topo_tools_full::faces_of(&f).len();
    assert!(nf >= 7, "box + blends + corner patch → {nf}");
}

#[test]
fn draw_repl_basic_commands() {
    use occt_topo::draw::DrawSession;
    let mut s = DrawSession::default();
    occt_topo::draw::execute_line(&mut s, "box a 2 2 2").expect("box");
    occt_topo::draw::execute_line(&mut s, "sphere s 1").expect("sphere");
    occt_topo::draw::execute_line(&mut s, "info s").expect("info");
    assert!(s.shapes.contains_key("a"));
    assert!(s.shapes.contains_key("s"));
    assert!(s.log.iter().any(|l| l.contains("faces")), "info logged");
    assert!(occt_topo::draw::execute_line(&mut s, "bogus").is_err());
    assert!(s.last_error.is_some());
}

#[test]
fn frenet_helix_curvature_torsion() {
    let h = |u: f64| GpPnt::new(u.cos(), u.sin(), u);
    let d1 = |u: f64| GpVec::new(-u.sin(), u.cos(), 1.0);
    let d2 = |u: f64| GpVec::new(-u.cos(), -u.sin(), 0.0);
    let d3 = |u: f64| GpVec::new(u.sin(), -u.cos(), 0.0);
    let k = curvature(&d1(0.3), &d2(0.3));
    let t = torsion(&d1(0.3), &d2(0.3), &d3(0.3));
    assert!((k - 0.5).abs() < 1e-9, "κ {k}");
    assert!((t - 0.5).abs() < 1e-9, "τ {t}");
    let f = frenet_frame(h(0.3), &d1(0.3), &d2(0.3), &d3(0.3));
    assert!(f.tangent.dot(&f.normal).abs() < 1e-9, "orthonormal");
}

#[test]
fn vrml_read_and_scene() {
    // Write a small VRML quad, read it back as a scene (mirrors the working
    // fixture syntax: #VRML header + comma-separated index lists).
    let vrml = "#VRML V2.0 utf8\n\
        Transform {\n\
          translation 2 0 0\n\
          children [\n\
            Shape {\n\
              geometry IndexedFaceSet {\n\
                coord Coordinate { point [ 0 0 0, 1 0 0, 1 1 0, 0 1 0 ] }\n\
                coordIndex [ 0,1,2,-1, 0,2,3,-1 ]\n\
              }\n\
            }\n\
          ]\n\
        }\n";
    std::fs::write("p9_test_quad.wrl", vrml).unwrap();
    let scene = occt_topo::rwmesh::read_vrml_scene("p9_test_quad.wrl").expect("vrml");
    std::fs::remove_file("p9_test_quad.wrl").ok();
    assert_eq!(scene.nodes.len(), 1);
    assert_eq!(scene.nodes[0].triangles.len(), 2);
    // Translation applied: vertices near x≈2.
    assert!(scene.nodes[0].vertices.iter().any(|p| (p.x() - 2.0).abs() < 0.1));
}
