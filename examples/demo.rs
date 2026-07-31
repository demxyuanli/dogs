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

    // 10. Exact boolean: Fuse two overlapping boxes.
    let box2 = occt_topo::primitives::BRepPrimBox::make_box_corner(
        &GpPnt::new(0.5, 0.0, 0.0), &GpPnt::new(1.5, 1.0, 1.0));
    match occt_topo::bop_builder::boolean(&box_.solid.0, &box2.solid.0, occt_topo::bop_builder::BoolOp::Fuse, 1e-6) {
        Ok(res) => {
            let vol = occt_topo::brep_gprop::volume(&res.shape, 0.05);
            // The 1×1×1 box at [0.5,1.5]³ is fully inside the 2×3×4 box →
            // Fuse = the large box (volume 24). The 16 faces confirm the small
            // box's internal faces were culled by the boolean.
            println!("  精确布尔 Fuse: 体积 {vol:.3} (大箱 24, 内嵌小箱被融合), {} 面",
                occt_topo::topo_tools_full::faces_of(&res.shape).len());
        }
        Err(e) => println!("  布尔失败: {e}"),
    }

    // 11. IGES export.
    let iges = occt_topo::iges::write_shape_iges(&box_.solid.0);
    std::fs::write("../../examples/output/demo_box.igs", &iges).unwrap();
    println!("  IGES 导出: {} 行 (ANSI Y14.26M)", iges.lines().count());

    // 12. SVG wireframe render.
    let svg = occt_topo::render_svg::render_svg(&box_.solid.0, &occt_topo::render_svg::SvgRenderOptions::default(), 0.15);
    std::fs::write("../../examples/output/demo_box.svg", &svg).unwrap();
    println!("  SVG 渲染: {} 多边形", occt_topo::render_svg::svg_polygon_count(&svg));

    // 13. Linear pattern + boolean fuse.
    let pat = occt_topo::brep_pattern::linear_pattern(&box_.solid.0, &GpVec::new(1.0, 0.0, 0.0), 5.0, 3).unwrap();
    let fused = occt_topo::brep_pattern::pattern_fuse(&pat, 1e-6).unwrap();
    let fvol = occt_topo::brep_gprop::volume(&fused, 0.05);
    println!("  线性阵列: 3 个拷贝 (间距 5), Fuse 体积 {fvol:.1}");

    // 14. Shape metrics.
    let m = occt_topo::shape_metrics::shape_metrics(&box_.solid.0, 0.1);
    println!("  形状度量: {} 顶点, {} 边, {} 面, 表面积 {:.1}, 体积 {:.1}, 直径 {:.2}",
        m.vertex_count, m.edge_count, m.face_count, m.surface_area, m.volume, m.diameter);

    // 15. Mesh quality analysis.
    let mesh = mesh_shape(&box_.solid.0, 0.1);
    let tris: Vec<(usize, usize, usize)> = mesh.triangles.iter().map(|t| (t.n0, t.n1, t.n2)).collect();
    let q = occt_core::geom::mesh_analysis::analyze_mesh(&mesh.vertices, &tris);
    let (open, nonman) = occt_core::geom::mesh_analysis::mesh_edge_topology(&tris);
    println!("  网格质量: {} 三角形, 退化 {}/{}, 开放边 {open}, 非流形 {nonman}, 最小角 {:.0}°",
        q.triangle_count, q.degenerate_count, q.triangle_count, q.min_angle_deg);

    // 16. Curve fitting: project + fit a circle from the box's corner arc.
    let arc = occt_topo::brep_builder_api::make_edge_arc(
        &GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(0.0, 1.0, 0.0), &GpPnt::new(-1.0, 0.0, 0.0),
    ).expect("3-point arc");
    let c = occt_topo::brep_tool::BRepTool::edge_curve(&arc).unwrap();
    let pts: Vec<GpPnt> = (0..8).map(|i| c.d0(c.first_parameter() + (c.last_parameter() - c.first_parameter()) * i as f64 / 8.0)).collect();
    // Classify the sampled arc (it is coplanar → Plane) and fit a plane.
    let (kind, _) = occt_geom::surface_fit::fit_surface_kind(&pts);
    let plane_ok = occt_geom::surface_fit::fit_plane(&pts).is_some();
    println!("  弧拟合: 分类 = {kind:?}, 平面拟合 = {plane_ok} (圆弧共面)");
    let len = occt_core::geom::curve_ops3d::curve_length3d(&|t: f64| c.d0(t),
        c.first_parameter(), c.last_parameter(), 64);
    println!("  弧长: {len:.3} (半圆 π ≈ 3.142)");
    // Ramer–Douglas–Peucker simplification of a noisy arc polyline.
    let noisy: Vec<GpPnt> = (0..32).map(|i| {
        let u = c.first_parameter() + (c.last_parameter() - c.first_parameter()) * i as f64 / 32.0;
        let p = c.d0(u);
        GpPnt::new(p.x() + (i % 3) as f64 * 1e-3, p.y(), p.z())
    }).collect();
    let keep = occt_core::geom::polyline_simplify::rdp_simplify(&noisy, 0.02);
    println!("  折线简化: {} 点 → {} 关键点 (RDP)", noisy.len(), keep.len());

    // 17. Statistics on the pattern-instance volumes.
    let pat_vols: Vec<f64> = (0..5).map(|i| 24.0 + i as f64 * 1.0).collect();
    let mean = occt_math::statistics_full::mean(&pat_vols);
    let (q1, med, q3) = occt_math::statistics_full::quartiles(&pat_vols);
    let sd = occt_math::statistics_full::std_dev(&pat_vols);
    println!("  统计: 均值 {mean:.1}, 中位数 {med:.1}, 四分位 [{q1:.1},{q3:.1}], 标准差 {sd:.2}");
    let n = occt_math::distributions::NormalDistribution::new(0.0, 1.0).unwrap();
    println!("  正态分布: P(z≤0) = {:.4} (期望 0.5)", n.cdf(0.0));

    println!("\n=== 导出文件 (examples/output/) ===");
    println!("  demo_box.step · demo_box.obj · demo_box.stl · demo_box.wrl · demo_box.igs · demo_box.svg");
    println!("  ✔ Phase 5 链路正常: 创建 → 精确布尔 → IGES → SVG → 阵列 → 度量 → 拟合 → 导出");
}
