# OCCT ↔ Rust 对齐与覆盖矩阵

> 审查日期：2026-07-31。
> Rust：5 crate · 49,237 行 · 709 测试（core 195 / math 121 / geom 44 / geom2d 38 / topo 311）
> 对照：`D:\source\occt-src`（~13,395 文件；已移植范围内 ~73 万行源码）

## 架构对齐（toolkit → crate）

| OCCT toolkit | Rust crate | 对齐度 |
|---|---|---|
| TKernel（24 包） | occt-core/kernel | 核心类型对齐；异常→`Result<T,OCCError>` |
| TKMath/gp+bnd+poly+bvh | occt-core/{gp,bnd,poly,bvh} | gp 34/42 类 (81%) |
| TKMath/math（54 类） | occt-math（36 模块） | ~67% |
| TKG3d（Geom 39 类） | occt-geom（29 模块） | ~64% 常用曲线曲面 |
| TKG2d（Geom2d 22 类） | occt-geom2d（15 模块） | ~64% |
| TKBRep | occt-topo（47 模块） | BRep/TopExp/TopTools 主干对齐 |

设计保真：`TopoShape{TShape,Location,Orientation}` ✓ · `Arc<dyn Curve/Surface>` 替代 Handle ✓ ·
`OCCError` 替代 `Standard_Failure` ✓ · math 1-based 索引 ✓。
**偏差**：几何存 `GeometryRegistry` 侧表（trait 对象不可 Clone/Debug）；布尔/求交/网格为采样或体素近似。

## 覆盖矩阵

### ✅ 已完成（≥80%）

| OCCT 包 | Rust | 覆盖 |
|---|---|---|
| gp（42 类） | gp/（34 类型） | 81% |
| Bnd（9） | bnd/（9 模块） | 100% |
| TopExp（2） | topexp + topo_tools_full | 100% |
| BRep_Tool / BRep_Builder | brep_tool / builder | 100% 语义 |
| ElCLib / ElSLib | elib/ | ~90% |
| BSplCLib / BSplSLib | bspl/（14 模块） | ~70% |
| STEPControl | step.rs（1,706 行，26 实体） | 读写闭环 |
| OBJ/STL/PLY/VRML (DE) | io/ + vrml.rs | 100% |

### ◐ 部分覆盖（25–70%）

| OCCT 包 | Rust | 覆盖 | 说明 |
|---|---|---|---|
| math（54） | occt-math 36 模块 | 67% | 缺分布采样等 |
| Geom（39） | 25 类型 + geom_api | 64% | |
| Geom2d（22） | 14 类型 + geom2d_api | 64% | |
| BRepPrimAPI（12） | primitives + brep_builder_api | 58% | Box/Cyl/Sphere/Cone/Torus/Prism/Revol/Pipe |
| BRepBuilderAPI（30） | brep_builder_api | 33% | MakeEdge/Vertex/Wire/Face/Polygon/Arc |
| BRepCheck（10） | shape_analysis + brep_measure | 50% | |
| BRepExtrema（18） | brep_extrema | 33% | |
| BRepGProp（13） | brep_gprop + gprop | 23% | |
| BRepAlgoAPI（10） | boolean_ops + solid_union | 30% | **体素近似** |
| BRepOffsetAPI（14） | loft + pipe | 36% | |
| BRepSweep（9） | sweep + sweep_revolve + pipe | 33% | |
| BRepFilletAPI（4） | fillet | 25% | |
| BRepMesh（59） | wireframe + shape_mesh | **3%** | 采样近似 |
| GeomAPI（11） | geom_api | ~80%（常用） | |
| Geom2dAPI（5） | geom2d_api | ~80% | |
| GCPnts（11） | gcpnts | ~45% | |
| Extrema（72） | geom_api + brep_extrema | <10% | 采样近似 |
| IntAna/IntCurvesFace（10） | face_face + geom_api | ~40% | |
| Poly（18） | poly/ | 22% | |
| BVH（34） | bvh/ | 12% | |
| ShapeAnalysis（19） | shape_analysis | ~25% | |
| GProp（10） | gprop + brep_gprop | ~40% | |

### ◔ 新增（Phase 5，2026-08-01）
- **精确布尔（平面多面体）**：inttools + bop_builder（Fuse/Cut/Common 体积验证）✅ 平面；NURBS 曲面仍采样近似
- **IGES 写器**（ANSI Y14.26M，POINT/LINE/ARC/NURBS/FACE/SHELL/MSB）✅
- **真实 BRepMesh**（Deflection 自适应四叉树）✅ 曲面细分；BRepMesh 完整 59 类仍缺
- **SVG 渲染**（投影+画家算法）· **特征建模**（boss/hole/protrusion/pocket）· **XCAF-lite**
- **数学深度**：分布/二维插值/稀疏矩阵 CG/FFT2D/特征值/统计/优化器/多边形/Delaunay/BVH/曲面拟合

### ⬜ 未移植（0%）

- **精确布尔（NURBS 曲面）**：IntPatch（40 类）、BOPAlgo 完整、TopOpeBRep（37 类）
- **特征深度**：TKFeat（BRepFeat、LocOpe）
- **精确边倒圆**：TKFillet（95,710 行源码）
- **TKShHealing · TKXMesh · TKMeshVS**
- **DataExchange 其余**：BinXCAF、XmlXCAF、GLTF、RWMesh
- **整个 Visualization**（TKOpenGl/TKV3d/TKService/TKIVtk）
- **Draw / TEST**

## 总体结论

| 指标 | 数值 |
|---|---|
| Rust 移植量 | 39,231 行（含测试） |
| OCCT 已移植范围源码 | ~73 万行（20 toolkit） |
| 行数比 | ~5.4% |
| 核心几何主干覆盖 | ~85% |
| 全量 OCCT 类覆盖（估算） | ~15–20% |

强项：底层几何内核 + 拓扑数据结构 + 交换（STEP/OBJ/STL/PLY/VRML）已形成可用闭环。
最大缺口：精确布尔 → 真实 BRepMesh → IGES/BinXCAF → Visualization。

## 语义偏差风险（与 OCCT"对齐"的差异点）

1. **体素布尔 / 采样求交结果不是精确 TopoDS** —— 离散网格结果与 OCCT 精确布尔在拓扑上不等价（已标注）。
2. `is_inside`/`point_in_mesh` 曾因射线双计数误判（已修复：重合命中参数去重）。
3. `make_box_corner` 曾只改 bbox 不改几何（已修复：经 prism 构建真实位置几何）。
4. 曲线/曲面 trait 对象无法下转型 —— 类型分类靠几何不变量采样（平面/球面可判，柱/锥近似）。
