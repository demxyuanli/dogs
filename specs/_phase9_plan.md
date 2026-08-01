# Phase 9 Plan — 最终缺口收尾（≥5,000 行）

> 覆盖 `_coverage.md` ⬜ 清单的 6 项。目标：本轮新增 Rust ≥ 5,000 行（实际力争 ~8,000）。

## 缺口 → 方案

| # | 缺口 | 方案 | 模块 | 预估 |
|---|---|---|---|---|
| 1 | NURBS 曲面布尔完整拓扑 | 一般曲面布尔升级为修剪 B-Rep 面（NURBS 表面 + 求交曲线边界线框） | bop_curved 扩展 | 1,400 |
| 2 | TKFillet 链式共享顶点角块 | 连续两边的球面角块补丁，链式不报错 | fillet_var 扩展 | 900 |
| 3 | RWMesh 纹理/VRML | OBJ/glTF 纹理坐标+材质、VRML 读入 | rwmesh 扩展 | 900 |
| 4 | 硬件 Visualization | 材质/光照/深度缓冲 PPM 光栅 + 相机轨道/平移/缩放 | viz_scene 扩展 | 1,200 |
| 5 | Draw/TEST | Draw_Interpretor-lite REPL + 脚本执行 | draw.rs | 800 |
| 6 | BRepGProp/Extrema 解析深度 | 解析面积/体积/质心 + 解析最近点 | gprop_analytic + extrema 扩展 | 1,100 |
| 7 | 我的并行模块 | — | — | 600 |

## 依赖
- 各 Agent 只写自己的文件；已有模块（bop_curved/fillet_var/rwmesh/viz_scene/extrema/brep_gprop）为扩展，读全再改。
- draw.rs 依赖多个模块做命令，但自身新文件。
