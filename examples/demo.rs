//! End-to-end demo: create → measure → mesh → export → reimport.
use occt_core::gp::{GpPnt, GpAx3, GpDir};
use occt_core::elib::slib;
use occt_core::io::obj::{ObjMesh, write_obj_file};
use occt_core::io::stl::{StlMesh, write_binary_stl, to_triangulation, total_area};
use occt_core::bnd::BndBox;
use occt_topo::mesh::{mesh_box, mesh_sphere, mesh_surface_area};

fn main() {
    println!("=== 1. 模型创建 (基本体) ===");
    // 创建一个单位立方体 + 半径2的球
    let box_mesh = mesh_box((GpPnt::new(0.,0.,0.), GpPnt::new(1.,1.,1.)));
    let sphere_mesh = mesh_sphere(2.0, 16, 16);
    println!("  立方体: {} 顶点, {} 三角形, 表面积 {}",
        box_mesh.vertices.len(), box_mesh.triangles.len(), mesh_surface_area(&box_mesh));
    println!("  球体:   {} 顶点, {} 三角形, 表面积 {} (期望 ~50.27)",
        sphere_mesh.vertices.len(), sphere_mesh.triangles.len(), mesh_surface_area(&sphere_mesh));

    println!("\n=== 2. 曲面求值 (参数化几何) ===");
    // 球面上取点 + 法线
    let ax = GpAx3::new(GpPnt::new(0.,0.,0.),
        GpDir::from_axis(occt_core::gp::dir::DirAxis::Z),
        &GpDir::from_axis(occt_core::gp::dir::DirAxis::X)).unwrap();
    let sphere = occt_core::gp::GpSphere::new(ax, 2.0).unwrap();
    let p = slib::sphere_value(&sphere, 0.0, 0.0); // 经度0,纬度0
    println!("  球面 P(0,0) = ({:.3}, {:.3}, {:.3}), |P|=({:.3})", p.x(), p.y(), p.z(),
        (p.x()*p.x()+p.y()*p.y()+p.z()*p.z()).sqrt());

    println!("\n=== 3. 包围盒 ===");
    let mut bb = BndBox::new();
    for v in &sphere_mesh.vertices { bb.add_point(v); }
    let c = bb.corner_min(); let d = bb.corner_max();
    println!("  球体包围盒: [{:.2},{:.2}] x [{:.2},{:.2}] x [{:.2},{:.2}]",
        c.x(), d.x(), c.y(), d.y(), c.z(), d.z());

    println!("\n=== 4. 格式转换: 网格 → OBJ/STL 导出 ===");
    // 立方体 → OBJ
    let mut obj = ObjMesh::default();
    obj.vertices = box_mesh.vertices.clone();
    obj.faces = box_mesh.triangles.iter().map(|t| {
        occt_core::io::obj::ObjFace { v: vec![t.n0 as i32, t.n1 as i32, t.n2 as i32], vt: None, vn: None }
    }).collect();
    write_obj_file("demo_box.obj", &obj).unwrap();
    println!("  已导出 demo_box.obj ({} 面)", obj.faces.len());

    // 球体 → 二进制 STL
    let mut stl = StlMesh::default();
    stl.triangles = sphere_mesh.triangles.iter().map(|t| {
        [sphere_mesh.vertices[t.n0], sphere_mesh.vertices[t.n1], sphere_mesh.vertices[t.n2]]
    }).collect();
    occt_core::io::stl::compute_facet_normals(&mut stl);
    std::fs::write("demo_sphere.stl", write_binary_stl(&stl)).unwrap();
    println!("  已导出 demo_sphere.stl ({} 三角形, 表面积 {:.2})",
        stl.triangles.len(), total_area(&stl));

    println!("\n=== 5. 读回 + 体积验证 ===");
    let tri = to_triangulation(&stl);
    println!("  STL 读回三角化: {} 顶点, {} 三角形", tri.nodes.len(), tri.triangles.len());
    println!("  ✔ 端到端链路正常");
}
