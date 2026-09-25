//! Phase 3 end-to-end integration: create real-geometry BRep → query →
//! measure → mesh → export → re-import.
//!
//! Exercises the milestone chain A1→A4 (geometry kernel + children tree),
//! A2/A3 (builder + BRepTool), A5-A7 (wireframe/shape_mesh/bbox), B
//! (curve/surface conversion via shape_mesh internals) and D (exchange).

use occt_core::gp::{GpPnt, GpVec};
use occt_core::io::obj::ObjMesh;
use occt_core::poly::Triangulation;
use occt_topo::brep_exchange::{brep_to_obj, brep_to_ply, brep_to_stl_ascii, brep_to_stl_binary};
use occt_topo::brep_surface::face_is_planar;
use occt_topo::brep_tool::BRepTool;
use occt_topo::bbox_from_geometry::shape_bbox;
use occt_topo::brep_extrema::is_inside;
use occt_topo::primitives::{BRepPrimBox, BRepPrimCylinder, BRepPrimSphere};
use occt_topo::shape_mesh::{mesh_shape, shape_surface_area, shape_volume};
use occt_topo::sweep::prism_from_polygon;
use occt_topo::topo_tools_full::{faces_of, vertices_of, structure_is_valid};

fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-3 * b.abs().max(1.0)
}

#[test]
fn real_box_full_pipeline() {
    // 1. Create a real-geometry box.
    let box_ = BRepPrimBox::make_box(2.0, 3.0, 4.0);
    let (nv, ne, nf) = box_.counts();
    assert_eq!((nv, ne, nf), (8, 12, 6));
    assert!(structure_is_valid(&box_.solid.0));

    // 2. Query geometry back.
    assert!(faces_of(&box_.solid.0).iter().all(|f| face_is_planar(f)));
    let verts = vertices_of(&box_.solid.0);
    let mut has_origin = false;
    let mut has_max = false;
    for v in &verts {
        let p = BRepTool::vertex_point(v);
        if p.is_equal(&GpPnt::zero()) {
            has_origin = true;
        }
        if p.is_equal(&GpPnt::new(2.0, 3.0, 4.0)) {
            has_max = true;
        }
    }
    assert!(has_origin && has_max);

    // 3. Measure (analytic volume + mesh area).
    assert!(approx(box_.volume(), 24.0));
    let mesh = mesh_shape(&box_.solid.0, 0.1);
    assert!(approx(shape_surface_area(&box_.solid.0, 0.1), 52.0));
    assert!(!mesh.triangles.is_empty());

    // 4. Real bounding box from geometry.
    let bb = shape_bbox(&box_.solid.0);
    let (x0, x1, _y0, y1, _z0, z1) = bb.get().expect("non-empty bbox");
    assert!(approx(x0, 0.0) && approx(x1, 2.0) && approx(y1, 3.0) && approx(z1, 4.0));

    // 5. Export OBJ/STL/PLY.
    let obj = brep_to_obj(&box_.solid.0, 0.2);
    assert!(obj.contains("v "));
    assert!(obj.contains("f "));
    let stl_ascii = brep_to_stl_ascii(&box_.solid.0, 0.2);
    assert!(stl_ascii.contains("solid"));
    let stl_bin = brep_to_stl_binary(&box_.solid.0, 0.2);
    assert!(stl_bin.len() > 80);
    let ply = brep_to_ply(&box_.solid.0, 0.2);
    assert!(ply.contains("ply"));
}

#[test]
fn primitives_measure_correctly() {
    let sphere = BRepPrimSphere::make_sphere(2.0);
    // Mesh-based volume of a closed sphere ~ (4/3)π·8 ≈ 33.5.
    assert!(approx(sphere.volume(), 33.510321638291124));
    let v = shape_volume(&sphere.solid.0, 0.15);
    assert!(v > 30.0 && v < 37.0, "sphere mesh volume {v}");
    assert!(is_inside(&sphere.solid.0, &GpPnt::new(0.5, 0.5, 0.5)));

    let cyl = BRepPrimCylinder::make_cylinder(1.0, 3.0);
    assert!(approx(cyl.volume(), 3.0 * std::f64::consts::PI));
    let faces = faces_of(&cyl.solid.0);
    assert_eq!(faces.len(), 3);
}

#[test]
fn sweep_prism_and_export_roundtrip() {
    // Extrude a right triangle → volume 0.5·base·height.
    let tri = [GpPnt::new(0.,0.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(0.,2.,0.)];
    let prism = prism_from_polygon(&tri, &GpVec::new(0.0, 0.0, 1.0), 5.0);
    assert!(approx(occt_topo::sweep::prism_volume(&prism), 10.0));
    assert_eq!((prism.vertices.len(), prism.edges.len()), (6, 9));

    // Export → re-import via OBJ parser.
    let obj = brep_to_obj(&prism.solid.0, 0.1);
    let mesh: Result<ObjMesh, String> = ObjMesh::parse(&obj);
    let mesh = mesh.expect("OBJ text parses");
    assert!(!mesh.vertices.is_empty());
    assert!(!mesh.faces.is_empty());
    let tri: Triangulation = mesh.to_triangulation();
    assert!(tri.triangles.len() > 0);
}

#[test]
fn bbox_and_mesh_agree() {
    let box_ = BRepPrimBox::make_box(1.0, 1.0, 1.0);
    let mesh_bb = occt_topo::mesh::mesh_bbox(&mesh_shape(&box_.solid.0, 0.1));
    let geom_bb = shape_bbox(&box_.solid.0);
    let (_, a1, _, b1, _, c1) = mesh_bb.get().unwrap();
    let (_, x1, _, y1, _, z1) = geom_bb.get().unwrap();
    assert!(a1 > 0.99 && a1 <= 1.01);
    assert!(x1 > 0.99 && x1 <= 1.01 && b1 == y1 && c1 == z1);
}
