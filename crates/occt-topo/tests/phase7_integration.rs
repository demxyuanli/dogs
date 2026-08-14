//! Phase 7 end-to-end: BRepFeat features (draft/groove/neck), variable-radius
//! fillet, shape healing, XML XCAF, and math depth (RK45).

use occt_core::gp::{GpDir, GpPnt};
use occt_math::ode_rk45::rk45_solve;
use occt_topo::brepfeat::draft;
use occt_topo::fillet_var::{fillet_edge_var, radius_at, RadiusLaw, VarFilletSpec};
use occt_topo::primitives::BRepPrimBox;
use occt_topo::shhealing::heal_shape;
use occt_topo::xmlcaf::{from_xml, to_xml, XmlAttribute, XmlEntry, XmlXcafDoc};

fn approx(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

#[test]
fn brepfeat_draft() {
    let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
    // Draft the top face (z=2 plane) by 10°, pivoting about a vertical plane
    // (x=0) so the hinge line is not parallel to the face.
    let top_pts = vec![GpPnt::new(1.0, 1.0, 2.0)];
    let pivot = occt_core::gp::GpPln::new(
        occt_core::gp::GpAx3::new(
            GpPnt::new(0.0, 0.0, 0.0),
            occt_core::gp::GpDir::new(1.0, 0.0, 0.0).unwrap(),
            &occt_core::gp::GpDir::new(0.0, 0.0, 1.0).unwrap(),
        )
        .expect("pivot ax3"),
    );
    let d = draft(&b.solid, &top_pts, 10.0, &pivot, 1e-6).expect("draft");
    assert!(d.volume > 0.0);
}

// NOTE: brepfeat neck/groove (revolve + boolean) is heavy in debug builds and
// is exercised in depth by brepfeat.rs's own 17 module tests; the integration
// suite keeps the fast draft feature here only.

#[test]
fn fillet_var_linear_profile() {
    let spec = VarFilletSpec { r_start: 0.3, r_end: 0.6, law: RadiusLaw::Linear };
    assert!(approx(radius_at(&spec, 0.0), 0.3, 1e-9));
    assert!(approx(radius_at(&spec, 1.0), 0.6, 1e-9));
    assert!(approx(radius_at(&spec, 0.5), 0.45, 1e-9));

    let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
    // Fillet the front-bottom edge (0,0,0)→(2,0,0) with variable radius.
    let edge = occt_topo::topo_tools_full::edges_of(&b.solid.0)
        .into_iter()
        .find(|e| {
            occt_topo::brep_tool::BRepTool::edge_vertices(e)
                .map(|(a, b2)| {
                    a.is_equal(&GpPnt::new(0.0, 0.0, 0.0)) && b2.is_equal(&GpPnt::new(2.0, 0.0, 0.0))
                })
                .unwrap_or(false)
        })
        .expect("front-bottom edge (0,0,0)→(2,0,0)");
    let f = fillet_edge_var(&b.solid.0, &edge, &spec, 1e-6).expect("var fillet");
    let nf = occt_topo::topo_tools_full::faces_of(&f).len();
    assert!(nf >= 7, "box (6 faces) + blend → {nf}");
}

#[test]
fn shhealing_removes_small_edge() {
    let b = occt_topo::builder::TopoBuilder::new();
    let p0 = GpPnt::new(0.0, 0.0, 0.0);
    let p1 = GpPnt::new(1e-4, 0.0, 0.0);
    let p2 = GpPnt::new(1.0, 0.0, 0.0);
    let p3 = GpPnt::new(1.0, 1.0, 0.0);
    let e01 = b.make_edge_segment(&p0, &p1);
    let e12 = b.make_edge_segment(&p1, &p2);
    let e23 = b.make_edge_segment(&p2, &p3);
    let wire = b.make_wire(&[e01, e12, e23]);
    let (healed, report) = heal_shape(&wire, 1e-6, 0.01);
    assert!(report.small_edges_removed >= 1, "tiny edge removed, report {report:?}");
    let ne = occt_topo::topo_tools_full::edges_of(&healed).len();
    assert!(ne <= 3, "edge count {ne}");
}

#[test]
fn xmlcaf_roundtrip_nested() {
    let mut doc = XmlXcafDoc::default();
    doc.root = XmlEntry {
        name: "assembly".into(),
        shape: None,
        attributes: vec![XmlAttribute { kind: "color".into(), value: "1,0,0".into() }],
        children: vec![XmlEntry {
            name: "part".into(),
            shape: Some(BRepPrimBox::make_box(1.0, 2.0, 3.0).solid.0),
            attributes: vec![],
            children: vec![],
        }],
    };
    let xml = to_xml(&doc).expect("to_xml");
    assert!(xml.starts_with("<?xml"));
    let back = from_xml(&xml).expect("from_xml");
    assert_eq!(back.root.name, "assembly");
    assert_eq!(back.root.attributes[0].value, "1,0,0");
    assert_eq!(back.root.children.len(), 1);
    let nf = occt_topo::topo_tools_full::faces_of(back.root.children[0].shape.as_ref().unwrap()).len();
    assert_eq!(nf, 6, "box faces after XML round-trip");
}

#[test]
fn math_depth_rk45() {
    // RK45: harmonic oscillator y'' = −y → sin(t).
    let (t_end, y) = rk45_solve(
        &|t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[1];
            dy[1] = -y[0];
            let _ = t;
        },
        0.0,
        &[0.0, 1.0],
        6.0,
        1e-8,
        1e-8,
    )
    .expect("rk45");
    assert!(approx(t_end, 6.0, 1e-9));
    assert!(approx(y[0], (6.0f64).sin(), 1e-4), "sin(6) ≈ {}", y[0]);
}
