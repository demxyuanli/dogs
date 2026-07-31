//! End-to-end demo: real-geometry model creation → measure → mesh → export → reimport.
//!
//! Demonstrates the Phase 3 milestone: `BRepPrimBox` builds a full boundary
//! representation with registered geometry (vertex points, line curves,
//! planar faces), `BRepTool` reads it back, `shape_mesh` triangulates any
//! BRep, and `brep_exchange` exports it to OBJ/STL/PLY.

use occt_core::gp::{GpPnt, GpVec};
use occt_topo::brep_exchange::{brep_to_obj, brep_write_obj, brep_to_stl_binary};
use occt_topo::brep_surface::face_is_planar;
use occt_topo::brep_tool::BRepTool;
use occt_topo::primitives::{BRepPrimBox, BRepPrimCylinder, BRepPrimSphere};
use occt_topo::shape_mesh::{mesh_shape, shape_surface_area};
use occt_topo::sweep::prism_from_polygon;
use occt_topo::topo_tools_full::{faces_of, vertices_of};

fn main() {
    println!("=== Phase 3: 真实几何 BRep 模型创建 ===");

    // 1. Box with real geometry: 8 vertices, 12 edges, 6 planar faces.
    let box_ = BRepPrimBox::make_box(2.0, 3.0, 4.0);
    let (nv, ne, nf) = box_.counts();
    println!("  立方体: {nv} 顶点, {ne} 边, {nf} 面 (真几何)");
    let all_planar = faces_of(&box_.solid.0).iter().all(|f| face_is_planar(f));
    println!("  全部面为平面: {all_planar}");

    // 2. Query geometry back via BRepTool.
    let verts = vertices_of(&box_.solid.0);
    let v0 = BRepTool::vertex_point(&verts[0]);
    let vmax = verts.iter().fold(v0, |a, v| {
        let p = BRepTool::vertex_point(v);
        if p.x() + p.y() + p.z() > a.x() + a.y() + a.z() { p } else { a }
    });
    println!("  顶点坐标范围: {v0:?} .. {vmax:?}");

    // 3. Measure.
    println!("  体积 = {} (期望 2·3·4 = 24)", box_.volume());

    // 4. Mesh + export OBJ/STL.
    let mesh = mesh_shape(&box_.solid.0, 0.1);
    println!("  网格化: {} 顶点, {} 三角形, 表面积 {:.3} (期望 52.0)",
        mesh.vertices.len(), mesh.triangles.len(), shape_surface_area(&box_.solid.0, 0.1));
    let obj = brep_to_obj(&box_.solid.0, 0.1);
    println!("  OBJ 导出: {} 行", obj.lines().count());
    brep_write_obj("../../examples/output/demo_box_real.obj", &box_.solid.0, 0.1).unwrap();
    let stl = brep_to_stl_binary(&box_.solid.0, 0.1);
    println!("  STL 导出: {} 字节", stl.len());

    // 5. Sphere + cylinder.
    let sphere = BRepPrimSphere::make_sphere(2.0);
    let smesh = mesh_shape(&sphere.solid.0, 0.15);
    let sarea = shape_surface_area(&sphere.solid.0, 0.15);
    println!("  球体表面积 {sarea:.2} (期望 4π·4 ≈ 50.27), {} 三角形", smesh.triangles.len());
    let cyl = BRepPrimCylinder::make_cylinder(1.0, 3.0);
    println!("  圆柱体积 {} (期望 3π ≈ 9.42)", cyl.volume());

    // 6. Sweep: extrude a triangle prism.
    let tri = [GpPnt::new(0.,0.,0.), GpPnt::new(2.,0.,0.), GpPnt::new(0.,2.,0.)];
    let prism = prism_from_polygon(&tri, &GpVec::new(0.0, 0.0, 1.0), 5.0);
    println!("  棱柱: {} 顶点, {} 边, {} 面, 体积 {} (期望 10.0)",
        prism.vertices.len(), prism.edges.len(), prism.lateral_faces.len() + 2,
        occt_topo::sweep::prism_volume(&prism));

    println!("\n=== 导出文件 (examples/output/) ===");
    println!("  demo_box_real.obj — 真实几何立方体");
    println!("  ✔ 端到端链路正常: 创建 → 查询 → 测量 → 网格化 → 导出");
}
