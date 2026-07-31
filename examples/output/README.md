# 测试数据 / Test data

由 `examples/demo.rs` 生成（`cargo run --example demo -p occt-topo`）。

保留这些文件用于验证导出管线，**不要删除**。

| 文件 | 内容 | 用途 |
|---|---|---|
| `demo_box.obj` | 单位立方体 (8 顶点, 12 三角形面) | OBJ 导出验证，可用 Blender/MeshLab 打开 |
| `demo_sphere.stl` | 半径 2 球体 (512 三角形, 表面积 ≈ 49.38 ≈ 4πr²) | 二进制 STL 导出验证，可用 3D 查看器打开 |
