# Phase 8 Plan（Trellis 收尾）— 规划剩余缺口，完整迁移翻译

> 依据 `specs/_coverage.md` 未移植清单 + `docs/trellis/plan.md` 分层。目标：本轮新增 Rust ≥ 10,000 行（最低 5,000），逐条消除剩余缺口。

## 剩余缺口全图（trellis 式：缺口 → 影响 → 方案）

| # | 缺口（OCCT 包） | 现状 | Phase 8 方案 | 预估 |
|---|---|---|---|---|
| 1 | **RWMesh 全量** | 仅 OBJ/PLY/STL 写 + OBJ/PLY 读 | glTF 读入 + STL 读 + 格式间转换 + 网格场景 | 1,300 |
| 2 | **BRepBuilderAPI 完整**（30 类中 ~10 类） | 仅 MakeEdge/Polygon/Arc 等 | MakeEdge 多曲线、MakeWire/Face/Shell/Solid/Vertex 统一包装 | 1,500 |
| 3 | **BVH / Poly 深度** | BVH 三角 + 遍历；Poly 基础 | BVH 线/面/盒/射线查询、Poly_Triangulation 工具、网格壳/边拓扑 | 1,400 |
| 4 | **Visualization-lite**（TKV3d 最小子集） | 仅 SVG/HLR | AIS_Shape 场景图、相机投影、正交/透视 PPM 渲染 | 1,400 |
| 5 | **math 深度** | 60+ 模块 | Adams-Bashforth 多步 ODE、罚函数约束优化、插值曲线构造 | 1,400 |
| 6 | **TKXMesh/TKMeshVS** | 无 | 网格细分细化、网格修复、网格到 TShell | 1,000 |
| 7 | **我的并行模块** | — | mesh_pipeline（拓扑链修复）+ geom2d 样条工具 | 900 |

## 里程碑

- **M1 (1+2)**: 交换读入闭环（glTF/STL read）+ 构造 API 完整
- **M2 (3)**: 空间查询与网格算法
- **M3 (4)**: 最小渲染器（场景→相机→投影→像素）
- **M4 (5+6+7)**: 数值/网格/样条收尾

## 依赖

- 各模块独立，Agent 只写自己的文件；lib.rs 已预注册。
- P8-A 依赖 occt-core io/（已有 OBJ/PLY 解析可复用）。
- P8-D 依赖 render_svg（投影思想）。
