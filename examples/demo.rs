//! End-to-end demo: real-geometry BRep → measure → mesh → STEP/VRML/OBJ/STL
//! export, fillet, revolve, and HLR projection.
//!
//! Demonstrates the Phase 3+4 milestones: real-geometry primitives, arbitrary
//! BRep meshing, STEP (ISO 10303-21) exchange, VRML/OBJ/STL export, corner
//! fillet, lathe revolve, and hidden-line-removal projection.

use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_geom::geom_api::project_point_on_curve;
use occt_topo::brep_exchange::{brep_to_obj, brep_to_stl_binary};
use occt_topo::brep_surface::face_is_planar;
use occt_topo::brep_tool::BRepTool;
use occt_topo::fillet::fillet_corner;
use occt_topo::hlr::wireframe_projection;
use occt_topo::primitives::BRepPrimBox;
use occt_topo::shape_mesh::{mesh_shape, shape_surface_area};
use occt_topo::step::write_shape_step;
use occt_topo::sweep_revolve::revolve_rectangle;
use occt_topo::topo_tools_full::{faces_of, vertices_of};
use occt_topo::vrml::write_vrml_shape;

fn main() {
    println!("=== Phase 3+4: 真实几何 BRep 模型全链路 ===");

    // 1. Real-geometry box.
    let box_ = BRepPrimBox::make_box(2.0, 3.0, 4.0);
    let (nv, ne, nf) = box_.counts();
    let all_planar = faces_of(&box_.solid.0).iter().all(|f| face_is_planar(f));
    println!("  立方体: {nv} 顶点, {ne} 边, {nf} 面, 全平面: {all_planar}, 体积 {}", box_.volume());

    // 2. STEP (ISO 10303-21) export.
    let step = write_shape_step(&box_.solid.0);
    println!("  STEP 导出: {} 字节 ({} 实体行)", step.len(), step.lines().count());
    std::fs::write("../../examples/output/demo_box.step", &step).unwrap();

    // 3. Mesh + OBJ/STL.
    let mesh = mesh_shape(&box_.solid.0, 0.1);
    println!("  网格化: {} 三角形, 表面积 {:.3} (期望 52.0)",
        mesh.triangles.len(), shape_surface_area(&box_.solid.0, 0.1));
    brep_to_obj(&box_.solid.0, 0.2);
    std::fs::write("../../examples/output/demo_box.obj", brep_to_obj(&box_.solid.0, 0.2)).unwrap();
    std::fs::write("../../examples/output/demo_box.stl", brep_to_stl_binary(&box_.solid.0, 0.2)).unwrap();
    println!("  OBJ/STL 导出完成");

    // 4. VRML export.
    let vrml = write_vrml_shape(&box_.solid.0, "box", 0.15);
    std::fs::write("../../examples/output/demo_box.wrl", &vrml).unwrap();
    println!("  VRML 导出: {} 字节", vrml.len());

    // 5. Fillet a corner.
    let p1 = GpPnt::new(2.0, 0.0, 0.0);
    let p2 = GpPnt::new(0.0, 0.0, 0.0);
    let p3 = GpPnt::new(0.0, 2.0, 0.0);
    let wire = fillet_corner(&p1, &p2, &p3, 0.4).expect("fillet");
    let edges = occt_topo::topo_tools_full::edges_of_wire(&wire);
    println!("  圆角: 90° 角 → {} 边 (2 直线 + 1 圆弧, r=0.4)", edges.len());

    // 6. Revolve → cylinder.
    let cyl = revolve_rectangle(1.0, 3.0, 24).expect("revolve");
    println!("  旋转体: 圆柱 半径1×高3 → 体积 {:.3} (期望 {:.3})",
        occt_topo::sweep_revolve::revolved_volume(&cyl),
        std::f64::consts::PI * 3.0);

    // 7. Projection query.
    let circle = occt_geom::GeomCircle::new(occt_core::gp::GpCirc::new(
        occt_core::gp::GpAx2::standard(), 1.0));
    if let Some(pr) = project_point_on_curve(&circle, &GpPnt::new(0.0, 2.0, 0.0), 1e-6) {
        println!("  投影: (0,2,0) → 单位圆, 距离 {:.3}", pr.distance);
    }

    // 8. Hidden-line removal.
    let vis = wireframe_projection(&box_.solid.0, &GpVec::new(0.0, 0.0, 1.0), 0.15);
    println!("  消隐: +Z 视角 → {} 条可见线段", vis.len());

    // 9. Geometry queries.
    let verts = vertices_of(&box_.solid.0);
    let v0 = BRepTool::vertex_point(&verts[0]);
    println!("  顶点[0] = {v0:?}");
    println!("\n=== 导出文件 (examples/output/) ===");
    println!("  demo_box.step (STEP) · demo_box.obj · demo_box.stl · demo_box.wrl (VRML)");
    println!("  ✔ Phase 4 链路正常: 创建 → STEP → 网格 → 圆角 → 旋转体 → 消隐 → 导出");
}
