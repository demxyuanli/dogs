//! Phase 4 end-to-end: STEP exchange, fillet, revolve, projection, HLR/VRML,
//! geometry transforms, and shape identification — the ModelingAlgorithms +
//! DataExchange milestone.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::gp::quaternion_ext::{from_axis_angle, rotate_point};
use occt_geom::geom_api::project_point_on_curve;
use occt_geom2d::geom2d_api::intersect_circle_line;
use occt_topo::brep_exchange::{brep_to_obj, brep_to_stl_ascii};
use occt_topo::brep_tool::BRepTool;
use occt_topo::brep_builder_api::{make_cylinder_full, make_edge_arc, make_face_from_polygon};
use occt_topo::fillet::{fillet_corner, tangent_points};
use occt_topo::primitives::BRepPrimBox;
use occt_topo::shape_naming::{find_vertex_at, geometric_hash};
use occt_topo::shape_ops::translated_copy;
use occt_topo::sweep_revolve::revolve_rectangle;
use occt_topo::step::{read_step, write_step, write_shape_step};

fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-3 * b.abs().max(1.0)
}

#[test]
fn step_box_roundtrip_pipeline() {
    // Create a real-geometry box → write STEP → read back → verify.
    let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
    let step_text = write_shape_step(&b.solid.0);
    assert!(step_text.starts_with("ISO-10303-21;"));
    assert!(step_text.contains("CARTESIAN_POINT"));
    assert!(step_text.contains("ADVANCED_FACE"));
    assert!(step_text.contains("MANIFOLD_SOLID_BREP"));

    let model = read_step(&step_text).expect("STEP text parses back");
    assert_eq!(model.len(), 1);
    let shape = model.get(0).expect("one shape").shape.clone();
    let counts = occt_topo::topo_tools_full::shape_counts(&shape);
    assert_eq!(counts[&occt_topo::ShapeType::Vertex], 8);
    assert_eq!(counts[&occt_topo::ShapeType::Edge], 12);
    assert_eq!(counts[&occt_topo::ShapeType::Face], 6);
}

#[test]
fn fillet_and_measure() {
    // Fillet a 90° corner.
    let p1 = GpPnt::new(2.0, 0.0, 0.0);
    let p2 = GpPnt::new(0.0, 0.0, 0.0);
    let p3 = GpPnt::new(0.0, 2.0, 0.0);
    let (t1, t2) = tangent_points(&p1, &p2, &p3, 0.5).expect("tangency");
    assert!(approx(t1.distance(&GpPnt::new(0.5, 0.0, 0.0)), 0.0));
    assert!(approx(t2.distance(&GpPnt::new(0.0, 0.5, 0.0)), 0.0));

    let wire = fillet_corner(&p1, &p2, &p3, 0.5).expect("fillet wire");
    let edges = occt_topo::topo_tools_full::edges_of_wire(&wire);
    assert_eq!(edges.len(), 3);
    // The arc edge is a quarter-circle of radius 0.5 centered at (0.5,0.5).
    let arc = BRepTool::edge_curve(&edges[1]).unwrap();
    let mid = arc.d0(0.5 * (arc.first_parameter() + arc.last_parameter()));
    assert!(approx(mid.distance(&GpPnt::new(0.5, 0.5, 0.0)), 0.5));
}

#[test]
fn revolve_and_prism_volume() {
    let cyl = revolve_rectangle(1.0, 3.0, 24).expect("revolved cylinder");
    assert!(approx(occt_topo::sweep_revolve::revolved_volume(&cyl), std::f64::consts::PI * 3.0));
    // Lateral area: use the open (single-segment) profile → 2πrh.
    let open = occt_topo::sweep_revolve::revolve_polyline_around_z(
        &[occt_core::gp::GpPnt2d::new(1.0, 0.0), occt_core::gp::GpPnt2d::new(1.0, 3.0)], 24,
    ).expect("open cylinder");
    assert!(approx(occt_topo::sweep_revolve::revolved_lateral_area(&open), 2.0 * std::f64::consts::PI * 3.0));
}

#[test]
fn projection_and_intersection() {
    // Project (0,2,0) onto the unit circle → distance 1.
    let circle = occt_geom::GeomCircle::new(occt_core::gp::GpCirc::new(
        occt_core::gp::GpAx2::standard(), 1.0));
    let pr = project_point_on_curve(&circle, &GpPnt::new(0.0, 2.0, 0.0), 1e-6).expect("project");
    assert!(approx(pr.distance, 1.0));

    // 2D circle/line exact intersection.
    let ax22 = occt_core::gp::GpAx22d::new(
        occt_core::gp::GpPnt2d::new(0.0, 0.0),
        occt_core::gp::GpDir2d::new(1.0, 0.0).unwrap(),
        occt_core::gp::GpDir2d::new(0.0, 1.0).unwrap(),
    ).unwrap();
    let c2 = occt_core::gp::GpCirc2d::new(ax22, 2.0);
    let l2 = occt_core::gp::GpLin2d::from_pnt_dir(
        occt_core::gp::GpPnt2d::new(0.0, 0.0),
        occt_core::gp::GpDir2d::new(1.0, 0.0).unwrap(),
    );
    let pts = intersect_circle_line(&c2, &l2);
    assert_eq!(pts.len(), 2);
    assert!(approx(pts[0].x().abs(), 2.0) && approx(pts[1].x().abs(), 2.0));
}

#[test]
fn geometry_transform_and_hash() {
    let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let moved = translated_copy(&b.solid.0, &GpVec::new(4.0, 0.0, 0.0)).unwrap();
    let v = find_vertex_at(&moved, &GpPnt::new(4.0, 0.0, 0.0), 1e-9).expect("moved origin vertex");
    assert!(BRepTool::vertex_point(&v).distance(&GpPnt::new(4.0, 0.0, 0.0)) < 1e-9);
    assert_ne!(geometric_hash(&b.solid.0), geometric_hash(&moved));
}

#[test]
fn quaternion_rotation() {
    let q = from_axis_angle(&occt_core::gp::GpDir::new(0.0, 0.0, 1.0).unwrap(), std::f64::consts::FRAC_PI_2);
    let p = rotate_point(&q, &GpPnt::new(1.0, 0.0, 0.0));
    assert!(approx(p.x(), 0.0) && approx(p.y(), 1.0));
}

#[test]
fn exchange_exporters_still_work() {
    let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let obj = brep_to_obj(&b.solid.0, 0.1);
    assert!(obj.contains("v "));
    let stl = brep_to_stl_ascii(&b.solid.0, 0.1);
    assert!(stl.contains("solid"));
}

#[test]
fn arc_edge_and_cylinder_full() {
    let arc = make_edge_arc(
        &GpPnt::new(1.0, 0.0, 0.0),
        &GpPnt::new(0.0, 1.0, 0.0),
        &GpPnt::new(-1.0, 0.0, 0.0),
    ).expect("3-point arc");
    let c = BRepTool::edge_curve(&arc).unwrap();
    let mid = c.d0(0.5 * (c.first_parameter() + c.last_parameter()));
    assert!(approx(mid.distance(&GpPnt::new(0.0, 1.0, 0.0)), 0.0));

    let cyl = make_cylinder_full(1.0, 2.0).expect("cylinder full");
    assert_eq!(occt_topo::topo_tools_full::faces_of(&cyl.solid.0).len(), 3);
    assert!(approx(std::f64::consts::PI * cyl.radius * cyl.radius * cyl.height, std::f64::consts::PI * 2.0));
}

#[test]
fn polygon_and_face_building() {
    let sq = [
        GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.),
    ];
    let face = make_face_from_polygon(&sq).expect("square face");
    assert!(occt_topo::brep_surface::face_is_planar(&face));
    let poly = occt_core::geom::polygon_ops::polygon_area3d(&sq);
    assert!(approx(poly, 1.0));
}
