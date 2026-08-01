# OCCT ↔ Rust 对齐与覆盖矩阵

> 审查日期：2026-07-31。
> Rust：5 crate · 89,209 行 · 1,259 测试（core 247 / math 202 / geom 68 / geom2d 56 / topo 686）
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
| BRepAlgoAPI（10） | bop_builder + boolean_ops | 45% | **精确（平面）+ 体素（曲面）** |
| BRepOffsetAPI（14） | loft + pipe | 36% | |
| BRepSweep（9） | sweep + sweep_revolve + pipe | 33% | |
| BRepFilletAPI（4） | fillet | 25% | |
| BRepMesh（59） | brepmesh + wireframe + shape_mesh | 25% | Deflection 自适应（平面精确/曲面细分） |
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

### ◔ 新增（Phase 6，2026-08-01）
- **IntPatch 曲面求交**（intpatch）：平面-球/球-球解析圆交线 + 平面-柱/锥/torus + 一般曲面 trace 回退 ✅
- **曲面布尔**（bop_curved）：球/柱面实体的 Fuse/Cut/Common（分类 + 重建）✅；NURBS 一般曲面仍近似
- **实体边倒圆**（fillet_edge）：滚动球圆柱/球面混合（盒边/盒角）✅
- **BRepOffset 深度**（brep_offset）：多边形 miter 偏移/面偏移/壳偏移（盒 ±d 精确）✅
- **BinXCAF**（bincaf）：二进制装配容器（拓扑+属性+子级）读写闭环 ✅
- **glTF 2.0**（gltf）：JSON+BIN 网格导出（positions/normals/indices，data-URI 可选）✅
- **多边形布尔**（polygon_boolean）：2D 交/并/差（凸 Sutherland-Hodgman + 一般 Weiler-Atherton）✅
- **几何极值**（geom extrema）：点-曲线/曲线-曲线/曲线-曲面/曲面-曲面最近距离 ✅
- **修复 gp 变换 bug**：circ/cone/cylinder/pln 方向向量被平移污染 → transforms_xyz_dir ✅

### ⬜ 未移植（0%）

- **BOPAlgo 完整（自交修复/多结果拓扑）/ TopOpeBRep 全拓扑**
- **TKFillet 全（曲面-曲面一般混合、多法向复杂角）**
- **TKOpenGl 硬件渲染 / TKService 字体纹理全量 / TKV3d 交互选择高亮 / TKIVtk**
- **Draw 完整 TCL（表达式解析/过程/命令表）**
- **BRepGProp 更全（体积惯性、张量高级）**
- **XCAF 完整 schema（视图/标注/实例树）/ STEP 全实体写**

## 总体结论

| 指标 | 数值 |
|---|---|
| Rust 移植量 | 89,209 行（含测试） |
| OCCT 已移植范围源码 | ~73 万行（20 toolkit） |
| 行数比 | ~12.2% |
| 核心几何主干覆盖 | ~85% |
| 全量 OCCT 类覆盖（估算） | ~26–30% |

### 覆盖矩阵行更新（2026-08-01，Phase 6）
| OCCT 包 | Rust | 覆盖 | 说明 |
|---|---|---|---|
| IntPatch / BOPAlgo（曲面） | intpatch + bop_curved | 15%→35% | 解析曲面求交 + 曲面实体布尔 |
| BRepFilletAPI / ChFi3d | fillet_edge | 25%→40% | 滚动球边倒圆（圆柱/球面混合） |
| BRepOffsetAPI / BRepOffset | brep_offset | 36%→50% | 多边形/面/壳偏移 |
| BinXCAF / RWMesh | bincaf + gltf | 0%→30% | 二进制装配容器 + glTF 2.0 |
| Extrema（点-曲线/曲面） | occt-geom extrema | <10%→35% | 4 类极值 |
| 平面布尔 2D | polygon_boolean | — | 交/并/差（凸+一般） |

### 覆盖矩阵行更新（2026-08-01，Phase 7）
| OCCT 包 | Rust | 覆盖 | 说明 |
|---|---|---|---|
| IntPatch 一般曲面 / BOPAlgo | intpatch + bop_curved 扩展 | 35%→50% | B样条曲面 trace 求交 + 一般曲面布尔 |
| BRepFeat / LocOpe（TKFeat） | brepfeat | 0%→40% | 拔模/凹槽/颈部/加强筋/通孔 |
| BRepFilletAPI 变半径 | fillet_var | 40%→55% | 沿边半径线性/二次/三次变化 + 链式 |
| ShapeHealing / ShapeFix | shhealing | 0%→45% | 焊点/小边移除/闭合线框/顶点移动 |
| XmlXCAF | xmlcaf | 0%→35% | XML 装配容器（bincaf 文本版） |
| Extrema2d | geom2d extrema2d | — | 2D 点/曲线极值 + 交点 |
| 网格处理 | mesh_ops | — | Laplacian 平滑/细分/减面 |
| 数值深度 | ode_rk45 + lsq_nonlinear | — | DOPRI5 自适应 ODE + Gauss-Newton |

### 覆盖矩阵行更新（2026-08-01，Phase 8）
| OCCT 包 | Rust | 覆盖 | 说明 |
|---|---|---|---|
| RWMesh | rwmesh | 0%→45% | OBJ/PLY/STL/glTF 读入 + 场景合成 + 格式转换 |
| BRepBuilderAPI | brep_builder_full | 33%→70% | MakeEdge/Wire/Face/Shell/Solid/Polygon 统一包装 |
| BVH / Poly | bvh_query | 12%→40% | 射线/线段/盒查询 + 网格拓扑/体积/分量 |
| Visualization-lite | viz_scene | 0%→25% | 场景图/相机/透视/正交渲染（SVG+PPM） |
| math 深度 | ode_multistep + constraint_opt | 67%→75% | Adams 多步 ODE + 罚函数/增广拉格朗日约束优化 |
| 曲线插值 | geom interp_curve | — | 全局 B样条/带切矢插值 |
| TKXMesh/TKMeshVS | mesh_pipeline | 0%→45% | 细分/平滑/减面/修复/一致性/到 BRep |
| 2D 曲线工具 | geom2d curve_tools2d | — | 3 点圆弧/2D 插值/均匀取点/圆角 |

### 覆盖矩阵行更新（2026-08-01，Phase 9）
| OCCT 包 | Rust | 覆盖 | 说明 |
|---|---|---|---|
| NURBS 曲面布尔完整拓扑 | bop_curved 扩展 | 50%→70% | 修剪 B-Rep 面（保持 NURBS/解析曲面 + 求交曲线边界线框） |
| TKFillet 链式共享顶点角块 | fillet_var 扩展 | 55%→70% | 滚动球角块球面补丁，连续边链式闭合 |
| RWMesh 纹理/VRML | rwmesh 扩展 | 45%→60% | OBJ .mtl 材质+UV、glTF 材质/纹理、VRML 2.0 读入 |
| 硬件 Visualization | viz_scene 扩展 | 25%→45% | Phong/材质/光照/深度缓冲光栅 + 相机轨道/平移/缩放/拾取 |
| Draw/TEST | draw | 0%→40% | Draw_Interpretor-lite REPL + 脚本 + 批量测试驱动 |
| BRepGProp/Extrema 解析深度 | gprop_analytic | 23%→45% | 解析面积/体积/质心（盒精确 24）+ 解析最近点 |
| 曲线微分几何 | core curve_frenet | — | Frenet 帧/曲率/挠率/弧长参数化 |

### 覆盖矩阵行更新（2026-08-01，Phase 10）
| OCCT 包 | Rust | 覆盖 | 说明 |
|---|---|---|---|
| 拓扑布尔完整（多实体/自交/退化） | bop_builder 扩展 | 45%→65% | boolean_multi/compound、自交检测、退化、结果分解 |
| TKFillet 曲面面混合 | fillet_curved | 55%→70% | 平面-球/柱/球-球滚动球混合（torus 混合面） |
| TKService 纹理/字体/多视图 | viz_scene 扩展 | 45%→60% | UV 纹理映射、5×7 位图字体标签、多视图网格渲染 |
| Draw 完整交互 | draw 扩展 | 40%→65% | TCL 变量/表达式、for/if、视图命令、形状变换 |
| BRepGProp 高阶 | gprop_analytic 扩展 | 45%→60% | 惯性张量、主惯量/主轴、解析曲线长 |
| XCAF 全 schema + STEP 读 | xcaf 扩展 | 35%→55% | 名称/颜色/图层/材质全属性 + STEP 实体解析 |
| 3D 曲线插值 | core curve_interp3d | — | Catmull-Rom/三次样条/折线 + 弧长重参数化 |

强项：底层几何内核 + 拓扑数据结构 + 交换（STEP/OBJ/STL/PLY/VRML）已形成可用闭环。
最大缺口：精确布尔 → 真实 BRepMesh → IGES/BinXCAF → Visualization。

## 语义偏差风险（与 OCCT"对齐"的差异点）

1. **体素布尔 / 采样求交结果不是精确 TopoDS** —— 离散网格结果与 OCCT 精确布尔在拓扑上不等价（已标注）。
2. `is_inside`/`point_in_mesh` 曾因射线双计数误判（已修复：重合命中参数去重）。
3. `make_box_corner` 曾只改 bbox 不改几何（已修复：经 prism 构建真实位置几何）。
4. 曲线/曲面 trait 对象无法下转型 —— 类型分类靠几何不变量采样（平面/球面可判，柱/锥近似）。
