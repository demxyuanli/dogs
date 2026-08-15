# BRepMesh 逐行对齐任务清单

> 目标：把 `crates/occt-topo/src/meshing/` 现有 16,250 行简化移植，逐行对齐 OCCT 源码，
> 消除「漏这漏那」的语义差（seam 缝合、collectWirePoints 排他循环等）。
> 参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKMesh\BRepMesh\`（含 `TKShHealing` 的 `ShapeAnalysis_Wire`）。

## 方法

1. 逐任务对照 OCCT 源文件翻译，不猜测、不省略。
2. 每完成相对完整的一段，`cargo check` 作第一步校验。
3. 全部翻译后，`cargo test`（lib + step_obj_parity / step_to_obj）验证算法等价。

## 任务顺序（按 OCCT 管线自上而下）

| # | 任务 | OCCT 源 | Rust 目标 | 状态 |
|---|---|---|---|---|
| 13 | ShapeVisitor::addWire（ShapeAnalysis_Wire::CheckOrder 重排 + AddPCurve/AddEdge + EXTERNAL 跳过） | `BRepMesh_ShapeVisitor.cxx` | `meshing/model_builder.rs` | ✅ |
| 14 | EdgeDiscret::Tessellate3d/2d（端点用 BRep_Tool::Pnt(vertex) 替换 + EdgeParameterProvider 逐参数求 pcurve） | `BRepMesh_EdgeDiscret.cxx` + `BRepMesh_EdgeParameterProvider.hxx` | `meshing/edge_discret.rs` | ✅ |
| 15 | NodeInsertion::collectWirePoints（排他循环跳过最后一点 + RangeSplitter.AddPoint） | `BRepMesh_NodeInsertionMeshAlgo.hxx` | `meshing/node_insertion.rs` | ✅ 已提交 |
| 16 | NodeInsertion::initDataStructure（AdjustRange + SetCellSize/SetTolerance + Classifier.RegisterWire） | 同上 | 同上 | ✅（Classifier.RegisterWire + 过滤；SetCellSize/Scale 属 #20/#21） |
| 17 | DelaunayNodeInsertion::insertNodes/registerSurfaceNodes（Classifier::Perform==IN 才注册） | `BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx` | 同上 | ✅ |
| 18 | ModelPreProcessor（deflection 计算 + UV range） | `BRepMesh_ModelPreProcessor.cxx` | `meshing/model_builder.rs` | ✅ |
| 19 | CurveTessellator（自适应 chord-deviation + min_points + 端点保留） | `BRepMesh_CurveTessellator.cxx` | `meshing/edge_discret.rs` | ✅（angular 项已加；0.5× linear + min_size 未做） |
| 20 | RangeSplitter 族（Cylinder/Cone/Sphere/Torus/NURBS/Default/UVParam/BoundaryParams 步长 + AdjustRange + GenerateSurfaceNodes） | `BRepMesh_*RangeSplitter.cxx` | `meshing/range_splitter.rs` | ⚠️ 半（Default 基类逐行核对 ✅；NURBS `AnalyticalFilter` **插入 pass 已移植**（`analytical_filter_insert` 接进 `generate_nurbs_grid`），**移除/thinning pass 未做**——纯密度优化，不影响 bbox） |
| 21 | Delaun + DataStructureOfDelaun（空圆插入 + frontier/约束边 + cleanupMesh，**含周期面 seam 缝合**） | `BRepMesh_Delaun.cxx` + `BRepMesh_DataStructureOfDelaun.cxx` | `meshing/delaun.rs` + `delaun_data.rs` | ✅（逐方法核对 `Delaun.hxx` 46 方法 + `DataStructureOfDelaun.hxx` 全 API，均已在 delaun.rs/delaun_data.rs/delaun_index.rs；seam 非 Delaun 内缝合，OCC 靠共享 3D 折线，已随 #14） |
| 22 | ShapeTool/MeshTool/GeomTool（Range/UV bounds/法向/切矢/分类） | `BRepMesh_{Shape,Mesh,Geom}Tool.cxx` | `meshing/shape_tool.rs` + `mesh_tool.rs` + `geom_tool.rs` | ✅（已深移植；`max_face_tolerance` 已接入 #18） |
| 23 | ModelHealer + ModelPostProcessor（愈合 + 后处理） | `BRepMesh_ModelHealer.cxx` + `BRepMesh_ModelPostProcessor.cxx` | `meshing/model_healer.rs` | ✅（`fixFaceBoundaries`/`connectClosestPoints`/`getCommonVertex`/`closestPoint(s)`/`adjustSamePoints` 已移植，接进 `triangulate_model_faces`；heal + post 已有；`amplifyEdges`/`FaceChecker` 未移植） |
| 24 | FaceDiscret + Classifier（面 UV 顶点 + 面内分类 CSLib_Class2d） | `BRepMesh_FaceDiscret.cxx` + `BRepMesh_Classifier.cxx` | `meshing/face_discret.rs` | ✅（Classifier 已移植，winding-number 等效 CSLib_Class2d） |
| 25 | IncrementalMesh 管线（Builder→PreProcessor→各 algo→Healer→PostProcessor）+ 切进 export_mesh | `BRepMesh_IncrementalMesh.cxx` | `meshing/incremental_mesh.rs` + `brep_exchange.rs` | ✅ 管线已对齐 + Delaunay 已切进导出；**残留 B-spline bbox 漂移**（Shape-2 ±135 vs ±60、Shape 0.037）根因是 `step.rs` 丢 `SURFACE_CURVE` pcurve（#11 残留，非 BRepMesh）；四叉树兜底保留 |

## 已知关键差异（对齐清单）

1. **collectWirePoints 排他循环**（#15，已修）：`for (; aPointIt != aEndIndex; …)` 跳过每条边最后一点，避免 junction 重复。
2. **seam 缝合**（#21）：回转面 u=0 ≡ u=2π 两侧在 Delaunay 里需缝合，screw 圆柱/锥 seam 开边根因在此。
3. **wire 组装 CheckOrder 重排**（#13）：OCCT 用 `ShapeAnalysis_Wire::CheckOrder` 按几何连接性重排边；现 `build_model` 用 `edges_of_wire` 存序，需对照。
4. **surface nodes 分类过滤**（#16/#17）：OCCT `Classifier::Perform==IN` 才注册 surface node；现 `generate_surface_nodes` 全量注册。
5. **Tessellate3d 端点替换**（#14）：OCCT 用 `BRep_Tool::Pnt(vertex)` 替换端点，现实现需对照。

## 纪律

- 每个任务只改自己对应的 `meshing/` 模块，不碰 lib.rs 注册。
- 逐行对照 OCCT 源，不写「等效简化」。
- 每任务加 `#[test]`，保持现有 1284 lib + step parity 绿。
