//! Phase 5 end-to-end: exact boolean, real BRepMesh, IGES, SVG render,
//! BVH queries, offset2d and the extended math/numerics.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::geom::curve_ops3d::{curve_length3d, curvature3d};
use occt_topo::bop_builder::{boolean, BoolOp};
use occt_topo::brep_exchange::brep_to_obj;
use occt_topo::brepmesh::incremental_mesh;
use occt_topo::iges::write_shape_iges;
use occt_topo::primitives::BRepPrimBox;
use occt_topo::render_svg::{render_svg, SvgRenderOptions};

fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() < 0.05 * b.abs().max(1.0)
}

fn box_at() -> (occt_topo::shape::TopoShape, occt_topo::shape::TopoShape) {
    let a = BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0;
    let b = BRepPrimBox::make_box_corner(
        &GpPnt::new(0.5, 0.0, 0.0),
        &GpPnt::new(1.5, 1.0, 1.0),
    ).solid.0;
    (a, b)
}

fn volume_of(shape: &occt_topo::shape::TopoShape) -> f64 {
    occt_topo::brep_gprop::volume(shape, 0.05)
}

#[test]
fn exact_boolean_fuse_cut_common() {
    let (a, b) = box_at();
    // Fuse: two unit boxes with 0.5 overlap → volume 1.5.
    let fused = boolean(&a, &b, BoolOp::Fuse, 1e-6).expect("fuse");
    assert!(fused.solid.is_some(), "fuse produces a solid");
    let v = volume_of(&fused.shape);
    assert!(approx(v, 1.5), "fuse volume {v} (expect 1.5)");

    // Cut: A − B → 0.5.
    let cut = boolean(&a, &b, BoolOp::Cut, 1e-6).expect("cut");
    let v = volume_of(&cut.shape);
    assert!(approx(v, 0.5), "cut volume {v} (expect 0.5)");

    // Common → 0.5.
    let common = boolean(&a, &b, BoolOp::Common, 1e-6).expect("common");
    let v = volume_of(&common.shape);
    assert!(approx(v, 0.5), "common volume {v} (expect 0.5)");
}

#[test]
fn real_brepmesh_and_iges() {
    // Adaptive mesh of a sphere → area near 4πr².
    let sphere = occt_topo::primitives::BRepPrimSphere::make_sphere(2.0);
    let inc = incremental_mesh(&sphere.solid.0, 0.05).expect("incremental mesh");
    let area = occt_topo::mesh::mesh_surface_area(&inc.mesh);
    assert!(approx(area, 4.0 * std::f64::consts::PI * 4.0), "sphere area {area}");

    // IGES export of a box.
    let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let iges = write_shape_iges(&b.solid.0);
    assert!(iges.starts_with('S'), "IGES starts with S record");
    assert!(iges.lines().all(|l| l.chars().count() == 80), "all IGES lines 80 chars");
    assert!(iges.contains("186"), "has MANIFOLD SOLID BREP");
    assert!(iges.contains("510"), "has FACE");
}

#[test]
fn svg_render_and_obj_export() {
    let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let opts = SvgRenderOptions::default();
    let svg = render_svg(&b.solid.0, &opts, 0.1);
    assert!(svg.contains("<svg"));
    assert!(svg.contains("</svg>"));
    assert!(occt_topo::render_svg::svg_polygon_count(&svg) > 0);

    let obj = brep_to_obj(&b.solid.0, 0.1);
    assert!(obj.contains("v "));
}

#[test]
fn bvh_ray_and_point_queries() {
    use occt_core::bvh::builder_tri::build_tri_bvh;
    use occt_core::bvh::bvh_ops::{bvh_point_in_mesh, bvh_ray_cast};
    let tris = vec![
        ((GpPnt::new(1.,0.,0.)), (GpPnt::new(1.,1.,0.)), (GpPnt::new(1.,1.,1.))),
        ((GpPnt::new(1.,0.,0.)), (GpPnt::new(1.,1.,1.)), (GpPnt::new(1.,0.,1.))),
        ((GpPnt::new(0.,0.,0.)), (GpPnt::new(0.,1.,0.)), (GpPnt::new(0.,1.,1.))),
        ((GpPnt::new(0.,0.,0.)), (GpPnt::new(0.,1.,1.)), (GpPnt::new(0.,0.,1.))),
    ];
    let bvh = build_tri_bvh(&tris, 2);
    assert!(bvh_point_in_mesh(&bvh, &tris, &GpPnt::new(0.5, 0.5, 0.5)));
    let hit = bvh_ray_cast(&bvh, &tris, &GpPnt::new(-1.0, 0.5, 0.5), &GpVec::new(1.0, 0.0, 0.0));
    assert!(hit.is_some());
    assert!((hit.unwrap().1 - 1.0).abs() < 1e-9);
}

#[test]
fn curve_ops_3d_measuring() {
    let f = |t: f64| GpPnt::new(2.0 * t.cos(), 2.0 * t.sin(), 0.0);
    let len = curve_length3d(&f, 0.0, std::f64::consts::PI, 64);
    assert!(approx(len, 2.0 * std::f64::consts::PI));
    let k = curvature3d(&f, 1.0, 1e-3);
    assert!(approx(k, 0.5), "curvature {k}");
}

#[test]
fn geom2d_offset_and_math_stats() {
    use occt_core::gp::GpPnt2d;
    let sq = [
        GpPnt2d::new(0.,0.), GpPnt2d::new(1.,0.), GpPnt2d::new(1.,1.), GpPnt2d::new(0.,1.),
    ];
    let off = occt_geom2d::offset2d::offset_polygon(&sq, 0.1).expect("offset");
    let area = occt_geom2d::offset2d::polygon_offset_area(&off, 0.0);
    assert!(approx(area, 1.44), "offset square area {area}");

    // Statistics + sparse CG.
    let xs: Vec<f64> = (0..10).map(|i| i as f64).collect();
    let ys: Vec<f64> = xs.iter().map(|x| 2.0 * x + 1.0).collect();
    let (b, a) = occt_math::statistics_full::linear_regression(&xs, &ys);
    assert!((b - 2.0).abs() < 1e-9 && (a - 1.0).abs() < 1e-9);

    let mut sm = occt_math::matrix_sparse::SparseMatrix::new(3usize);
    for i in 0usize..3 {
        sm.set(i, i, 2.0);
        if i + 1 < 3 {
            sm.set(i, i + 1, 1.0);
            sm.set(i + 1, i, 1.0);
        }
    }
    let bvec = vec![4.0, 8.0, 8.0]; // A·(1,2,3) = (4,8,8)
    let x = occt_math::matrix_sparse::conjugate_gradient(&sm, &bvec, &[0.0, 0.0, 0.0], 100, 1e-8).unwrap();
    assert!((x[0] - 1.0).abs() < 1e-6 && (x[1] - 2.0).abs() < 1e-6 && (x[2] - 3.0).abs() < 1e-6);
}

#[test]
fn shell_orientation_and_shape_checks() {
    let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let shells = occt_topo::topo_tools_full::shapes_of(&b.solid.0, occt_topo::ShapeType::Shell);
    let shell = occt_topo::shape::Shell(shells[0].clone());
    assert!(occt_topo::brep_shell::shell_is_valid(&shell));
    let (errors, _) = occt_topo::shape_checks::analyse(&b.solid.0);
    assert!(errors.is_empty(), "box has no shape-check errors: {errors:?}");
}
