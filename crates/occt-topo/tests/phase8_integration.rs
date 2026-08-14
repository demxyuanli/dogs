//! Phase 8 end-to-end: RWMesh import, BRepBuilderAPI full, visualization-lite,
//! BVH queries, and math depth.

use occt_core::bvh::builder_tri::{build_tri_bvh, TriBvh};
use occt_core::bvh::bvh_query::{ray_cast_mesh, segment_query_mesh};
use occt_core::gp::{GpPnt, GpVec};
use occt_topo::brep_builder_full::BRepBuilderSolid;
use occt_topo::primitives::BRepPrimBox;
use occt_topo::rwmesh::{mesh_from_shape, scene_bounds, scene_to_compound};
use occt_topo::viz_scene::{camera_view_matrix, project_point, Camera, CameraProjection, SceneShape, VizScene};

fn box_mesh() -> (Vec<GpPnt>, Vec<(usize, usize, usize)>) {
    let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
    let m = occt_topo::shape_mesh::mesh_shape(&b.solid.0, 0.1);
    (m.vertices, m.triangles.iter().map(|t| (t.n0, t.n1, t.n2)).collect())
}

#[test]
fn rwmesh_scene_roundtrip() {
    let scene = mesh_from_shape(&BRepPrimBox::make_box(2.0, 3.0, 4.0).solid.0, 0.2);
    assert_eq!(scene.nodes.len(), 1);
    assert!(scene.nodes[0].triangles.len() > 0);
    // Scene bounds contain the box corners.
    let (mn, mx) = scene_bounds(&scene).expect("bounds");
    assert!(mn.x().abs() < 1e-9 && mn.y().abs() < 1e-9 && mn.z().abs() < 1e-9);
    assert!((mx.x() - 2.0).abs() < 0.1 && (mx.y() - 3.0).abs() < 0.1 && (mx.z() - 4.0).abs() < 0.1);
    // Compound from the scene has faces.
    let comp = scene_to_compound(&scene, 1e-6).expect("compound");
    let nf = occt_topo::topo_tools_full::faces_of(&comp).len();
    assert!(nf > 0, "faces {nf}");
}

#[test]
fn brep_builder_solid_box() {
    let solid = BRepBuilderSolid::box_corners(GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(2.0, 2.0, 2.0))
        .expect("box solid");
    let nv = occt_topo::topo_tools_full::vertices_of(&solid.0).len();
    assert_eq!(nv, 8, "box vertices {nv}");
    let nf = occt_topo::topo_tools_full::faces_of(&solid.0).len();
    assert_eq!(nf, 6, "box faces {nf}");
}

#[test]
fn viz_camera_projection() {
    let cam = Camera {
        eye: GpPnt::new(0.0, 0.0, 5.0),
        target: GpPnt::new(0.0, 0.0, 0.0),
        up: GpVec::new(0.0, 1.0, 0.0),
        fov_deg: 45.0,
        projection: CameraProjection::Perspective,
        near: 0.1,
        far: 100.0,
    };
    let w = 200;
    let h = 200;
    // Origin maps to screen center.
    let (cx, cy, _) = project_point(&cam, GpPnt::new(0.0, 0.0, 0.0), w, h).expect("origin");
    assert!((cx - w as f64 / 2.0).abs() < 2.0, "cx {cx}");
    assert!((cy - h as f64 / 2.0).abs() < 2.0, "cy {cy}");
    // +X maps to the right half.
    let (px, _, _) = project_point(&cam, GpPnt::new(1.0, 0.0, 0.0), w, h).expect("+x");
    assert!(px > w as f64 / 2.0, "px {px}");
    // Behind the camera → None.
    assert!(project_point(&cam, GpPnt::new(0.0, 0.0, 6.0), w, h).is_none());
    // View matrix has a non-trivial rotation (forward row non-zero).
    let m = camera_view_matrix(&cam);
    let fwd = occt_core::gp::GpVec::new(m.m[2][0], m.m[2][1], m.m[2][2]);
    assert!(fwd.xyz().modulus() > 0.5, "forward row non-zero");
}

#[test]
fn bvh_ray_and_segment_queries() {
    let (v, t) = box_mesh();
    let tris: Vec<(GpPnt, GpPnt, GpPnt)> = t.iter().map(|&(a, b, c)| (v[a], v[b], v[c])).collect();
    let bvh: TriBvh = build_tri_bvh(&tris, 4);
    // Ray from below the box along +Z hits it.
    let hit = ray_cast_mesh(&bvh, &tris, GpPnt::new(1.0, 1.0, -5.0), GpVec::new(0.0, 0.0, 1.0), 100.0);
    assert!(hit.is_some(), "ray hits box");
    // A short segment that stops before the box → None.
    let seg = segment_query_mesh(&bvh, &tris, GpPnt::new(1.0, 1.0, -5.0), GpPnt::new(1.0, 1.0, -0.5));
    assert!(seg.is_none(), "segment stops short");
}

#[test]
fn viz_scene_add_and_render() {
    let mut scene = VizScene::new();
    let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    scene.add(SceneShape::new(b.solid.0));
    assert_eq!(scene.len(), 1);
    let cam = Camera::default();
    let svg = occt_topo::viz_scene::render_scene_svg(&scene, &cam, 128, 128, 0.15);
    assert!(svg.contains("<polygon") || svg.contains("polygon"), "svg has polygons");
}
