> ⚠️ **仅历史（T-29 登记）**：本文件是 2026-09-09 的覆盖量测快照，行数/百分比与"未移植件"清单均已过期。
> 现行有效基线与任务状态以 `specs/_board.md` §2 / §3 为准；本文件不参与任何门禁。
>
> 审查日期：**2026-09-09**（对照 `D:\source\OCCT-src` OCC 8.0.0；刷新量测）。
> 前一版：2026-08-22。
> Rust：5 crate · **~216,620** 行源码（core 22,572 / math 12,446 / geom 14,358 / geom2d 4,568 / topo **162,676**）+ topo tests 2,144 · **741** 个 `.rs`。
> OCCT 8.0.0：`src/` hxx+cxx+lxx **~2,775,909** 行（~13,769 文件）；含 gxx 等 **~3,097,008**；`.hxx` ~7,072。
> CAD 内核（Foundation+ModelingData+ModelingAlgorithms）**~1,814,140** 行。
> 行数比：全库 **~7.9%**；CAD 内核 **~12.1%**。类/API 加权（CAD 内核）约 **48–52%**（较 8-22 的 ~45% 上浮，主因 IntPatch/PaveFiller/BOP 加深）。
>
> topo 簇行数（约）：int 28k · bop 27k · meshing 18k · pave 16k · brep 12k · fillet 8k · step/viz/draw/xcaf 各 2–4k。
> 源码标记粗扫：`Source: OCCT…` ~1.5k 处；`not ported`/`stub`/`unimplemented` 类注释 ~167；`approx`/`fallback`/`voxel`/`heuristic` ~357。

## 2026-09-09 再审查（对照 8.0.0）

| 层 | OCCT loc（现测） | 对齐 | 最大缺口 / 缺陷 |
|---|---|---|---|
| Foundation (TKernel+TKMath) | ~341k | ~70% | gp 长尾；**严重缺陷**：`GeomCylinder/Sphere/Torus::d1` 返回零切矢 |
| ModelingData | ~562k | ~55–60% | Adaptor/RectangularTrimmed/LinearExtrusion；几何仍在 `GeometryRegistry` 侧表 |
| ModelingAlgorithms | ~911k | ~42–48% | ImpPrm RLine/端点未完；`GetFaceOff` 未移植；ChFi3d/TopOpe 浅 |
| DataExchange | ~651k | ~14–18% | STEP 实体子集 vs 42 包；IGES/XCAF 非全 schema |
| Viz / Draw / OCAF | ~597k | ~8–20% | lite 替代（viz_scene / draw / xcaf），不算 toolkit 完成 |

**绿路径（已接线）**
- 布尔：`boolean` → `bop_builder2::builder_bop_with_fuzzy` → PaveFiller（含 RepeatIntersection / ForceInterfEE/EF / ProcessDE）→ FillImages* → BuildBOP。平面 `boolean_planar_legacy` **不再**作回退；无面 operand 仍 **voxel_fallback**。
- 面交：`PatchIntersection` = ImpImp →（Fail）PrmPrm；单侧解析二次曲面 → ImpPrm（SearchInside+IWalking）；否则 PrmPrm。IntPatch 相关 ~**8.3k** 行。
- **整体测试导出（STEP→OBJ）**：`rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj`。读 `data/*.step`，偏转 0.1，写仓库根 `output/<stem>.obj`。对照基线仍是 `tests/step_obj_parity` / `data/occ-*.obj`。

**语义偏差 / 缺陷（现行代码）**
1. **`cylinder.rs` / `sphere.rs` / `torus.rs` 的 `d1`** — 已改走 `surface_eval` / ElSLib（2026-09-09）。
2. **`GeometryRegistry` 侧表** — 边/面几何不在 `TShape` 内；`transform_shape` / `transformed_copy` 仅在有 `vertex_geom` 时写回，不再对未注册顶点写原点。
3. **`edge_vertices`** — `TopExp::FirstVertex/LastVertex`（存储 FORWARD/REVERSED）；`MakeEdge` / 盒/柱/锥原语按 `BRepLib_MakeEdge` 写向。其它 `add` 位点仍可能两个都是 FORWARD（回退插入顺序）。
4. **`GetFaceOff` / 角法向** — live 路径 `algo_tools_face::get_face_off` / `is_internal_face`。leftover `is_covering_face` 仅死路径 + 合成测试，未接到 `BopBuilder`。
5. **ImpPrm** — SearchInside+IWalking；`fleche` 已接入偏转；IWLine 端点 Destination/Recadre/MakeTransition 已接。HVertex 合并（cxx 329–465）仍标未移植。
6. **体素残留** — 无面 operand 仍 `voxel_fallback`。`validate` 的 voxel 体积对拍已去掉。绿路径仍是 `bop_builder2`。
7. **BRepMesh** — `FaceChecker` + `amplifyEdges` 接到 `IncrementalMesh::heal_self_intersecting_wires`（锥缝 `amplify_cone_seams` 仍在 `ModelPreProcessor`）。

## 2026-08-22 再审查（历史快照）

| 层 | OCCT loc | 对齐 | 最大缺口 |
|---|---|---|---|
| Foundation (TKernel+TKMath) | 223k | ~68% | gp 缺 Euler/GTrsf2d/NLerp；TKernel 无 Storage/Plugin |
| ModelingData | 324k | ~55% | Geom Bezier 已导出；仍缺 RectangularTrimmed / LinearExtrusion / Adaptor |
| ModelingAlgorithms | 862k | ~38% | PaveFiller 序列不全；两条布尔管线；TKGeomAlgo/ChFi3d/TopOpe |
| DataExchange | 603k | ~14% | STEP 42 包 vs `step.rs` lite |
| Viz / Draw / OCAF | 556k | ~8–18% | 有意的 lite 替代，不算 toolkit 完成 |

**语义偏差（8-22）**：几何侧表 `GeometryRegistry`；长尾仍直接读 `TShape.children`；`edge_vertices` 非 FirstVertex；当时记 PaveFiller 缺 RepeatIntersection / ForceInterf / ProcessDE — **2026-09 已接线进 `perform_internal`**。

---

# OCCT ↔ Rust 对齐与覆盖矩阵

> 审查日期：2026-08-04。
> Rust：5 crate · 158,239 行 · 2,052 测试（core 302 / math 215 / geom 139 / geom2d 72 / topo 1340）
> STEP→OBJ 对拍门禁（`tests/step_obj_parity.rs`，oracle=data/occ-*.obj）：Cube/Cone 精确、Sphere/Torus 采样差、Cylinder/rev 宽容差（缺口，Phase 21 移植）
> BRepMesh 59 类迁移完成（TKMesh，4 波，~17,400 行新增：框架/数据/Delaunay/细化/分割器/愈合/工厂）
> Phase 13（TKGeomAlgo 深度）：Extrema 解析/牛顿化 + IntAna + GCPnts + GeomConvert + IntCurvesFace
> Phase 14（基础补全）：Poly + GProp + BRepGProp + math GlobOptMin/FRPR（+7,286 行）
> Phase 15（精确 NURBS 布尔·机械波）：BOPDS + Area/Shell/Wire 构建器 + Options/历史（+5,340 行）
> Phase 16（精确布尔·波 A）：IntTools 机械层 + pcurve 基建（+5,118 行）
> Phase 17（精确布尔·波 C1）：曲线-曲面求交 + FClass2d + BeanFaceIntersector + EdgeFace（+5,692 行）
> Phase 18（精确布尔·波 C2a）：EdgeEdge 一般曲线 + FaceFace 精确 + AlgoTools 族 + Context/2D 完整（+5,130 行）
> Phase 19（精确布尔·波 C2b-1）：PaveFiller 全量求交管线（+4,997 行）
> Phase 20（精确布尔·波 C2b-2）：BOPAlgo_Builder 重建收官——交叠盒 Fuse/Cut/Common 体积精确（+4,415 行）
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
| math（54） | occt-math 36 模块 | 85% | GlobOptMin/FRPR 补全，残余碎片化 |
| Geom（39） | 25 类型 + geom_api | 64% | |
| Geom2d（22） | 14 类型 + geom2d_api | 64% | |
| BRepPrimAPI（12） | primitives + brep_builder_api | 58% | Box/Cyl/Sphere/Cone/Torus/Prism/Revol/Pipe |
| BRepBuilderAPI（30） | brep_builder_api | 33% | MakeEdge/Vertex/Wire/Face/Polygon/Arc |
| BRepCheck（10） | shape_analysis + brep_measure | 50% | |
| BRepExtrema（18） | brep_extrema | 33% | |
| BRepGProp（13） | brep_gprop + brep_gprop_full | 80% | 精确 Gauss 积分（Linear/Surface/Volume/GK） |
| BRepAlgoAPI（10） | bop_builder + boolean_ops | 45% | **绿路径走 BOPAlgo**（`boolean` → `bop_builder2`）；平面 2D 排列仅作未闭合/报错回退；NURBS 曲面仍采样近似 |
| BOPAlgo/BOPDS/BOPTools | bopds + builder_area/face + shell/wire_splitter | 35%→55% | Phase 15 机械波：DS 信息中枢 + 闭合面/壳/线构建 + Options/历史 |
| IntTools（24 类） | inttools_data/range/sample/roots + pcurve + intcurvesurface + fclass2d + bean_face + edge_face + edge_edge + int_face_face + int_tools_full + int_curve | 0%→85% | Phase 16-18：求交核全（解析 <1e-9，一般 <1e-6）；PaveFiller 接入波 C2b-2 |
| BOPAlgo_PaveFiller（12 文件 10,668 行） | pave_filler/intersect/blocks/common（4,997 行） | 0%→70% | Phase 19：VV/VE/EE/VF/EF/FF + MakeBlocks(CommonBlock) + MakePCurves + MakeSplitEdges + ShrunkData/自交检测 |
| BOPAlgo_Builder（5,372 行） | bop_builder2/build_faces/build_common/build_solids（4,415 行） | 0%→70% | Phase 20：FillImages* + BuildSplitFaces/Solids + BOP Fuse/Cut/Common——**交叠盒布尔体积精确**（fuse 1.5/cut 0.5/common 0.5） |
| BRepOffsetAPI（14） | loft + pipe | 36% | |
| BRepSweep（9） | sweep + sweep_revolve + pipe | 33% | |
| BRepFilletAPI（4） | fillet | 25% | |
| BRepMesh（59） | brepmesh + wireframe + shape_mesh | 25% | Deflection 自适应（平面精确/曲面细分） |
| GeomAPI（11） | geom_api | ~80%（常用） | |
| Geom2dAPI（5） | geom2d_api | ~80% | |
| GCPnts（11） | gcpnts | ~45% | |
| Extrema（72） | geom_api + brep_extrema | <10% | 采样近似 |
| IntAna/IntCurvesFace（10） | face_face + geom_api | ~40% | |
| Poly（18） | poly/（9 模块） | 85% | Triangulation/Coherent/Connect/MakeLoops/MergeNodes 全 |
| BVH（34） | bvh/ | 12% | 类保真留待（TriBvh 已覆盖实际查询） |
| ShapeAnalysis（19） | shape_analysis | ~25% | |
| GProp（10） | gprop/（gprops.rs 框架） | 85% | GProps 累加 + PrincipalProps + PG/Sel/Vel/Cel + PEquation |
| BRepGProp（13） | brep_gprop_full（2,085 行） | 80% | 精确 Gauss 积分：Linear/Surface/Volume/GK，盒/球/柱对拍 1e-6 |

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

### 覆盖矩阵行更新（2026-08-02，Phase 12）
| OCCT 包 | Rust | 覆盖 | 说明 |
|---|---|---|---|
| BOPAlgo/TopOpeBRep 深度 | bop_builder 扩展 | 65%→75% | 边重叠修复、面沿交线分割、边/顶点分类（Shared/OnFace/Internal/External）、容差愈合 |
| TKFillet 多法向复杂角 | fillet_curved 扩展 | 70%→80% | 三边曲面角块确定性补丁（torus 扇区+球带+球角块），链首尾闭合 |
| TKV3d 交互选择/高亮 | viz_scene 扩展 | 60%→70% | 射线拾取 pick_shape/pick_point、选中高亮渲染、双字体、屏幕文字叠加 |
| Draw 过程/函数 + 命令表 | draw 扩展 | 65%→80% | proc/def/call/return、递归调用栈、数组/lappend/concat、集中命令注册表 |
| XCAF 视图/标注/实例树 | xcaf 扩展 | 55%→70% | 命名视图、尺寸/注释标注、带放置实例树 + 展开、STEP 装配写 |

### 覆盖矩阵行更新（2026-08-02，BRepMesh 59 类迁移完成）
| OCCT 包 | Rust | 覆盖 | 说明 |
|---|---|---|---|
| BRepMesh（59 类，TKMesh） | meshing/（24 模块，~17,400 行） | 25%→90% | 现代增量管线逐类移植：IMeshTools/IMeshData 框架、BRepMeshData 模型、IncrementalMesh、GeomTool/Deflection、Edge/Face 离散化、UV Delaunay 三角化（Delaun + DataStructure）、deflection 控制细化（DelaunayDeflectionControlMeshAlgo）、RangeSplitter 族（解析+周期 seam）、ModelHealer/PostProcessor、工厂、Delabella 备选算法、Triangulator/FastDiscret，接入 shape_mesh |

### 覆盖矩阵行更新（2026-08-02，Phase 13 TKGeomAlgo 深度）
| OCCT 包 | Rust | 覆盖 | 说明 |
|---|---|---|---|
| Extrema（72 类，TKGeomBase） | extrema_pc + extrema_cc + extrema_surf + extrema_ss（~5,000 行） | <10%→75% | 采样近似→解析/牛顿：ExtPElC 点-直线/圆/椭/双曲/抛物精确解、ExtPC/ExtPC2d 网格+牛顿（d1/d2）、ExtCC/ExtElC（斜交线/线-圆/圆-圆解析配对）、ExtPElS 点-平面/球/柱/锥/torus 精确解、ExtPS/ExtCS/ExtSS 数值 Jacobian 牛顿（Surface trait 无 d2,ponytail 标记）、ExtPRevS/ExtPExtS 走通用路径；`extrema.rs` 改为纯门面,4 个采样函数全部切新实现（签名保留） |
| IntAna（8 类）+ IntCurvesFace（2） | intana.rs（957 行）+ brep_face_intersect.rs | 40%→85% | Int3Pln 三平面、IntQuadQuad/QuadQuadGeo 解析（平面-平面/平面-球/球-球/平面-柱/平面-锥含二次曲线）、IntLinTorus;IntCurvesFace 薄模块复用 inttools/intpatch/face_face + intana 平面∩球/球∩球精确圆 |
| GCPnts（11） | gcpnts.rs（466 行） | 45%→85% | AbscissaPoint（自适应积分+牛顿）、UniformAbscissa、QuasiUniform、UniformDeflection、TangentialDeflection |
| GeomConvert（16） | convert_bspl.rs（510 行） | 部分→70% | KnotSplitting（曲线+曲面）、BSplineSurfaceToBezierSurface（全重数插入,rational 感知）、CompCurveToBSplineCurve;BSplineCurveToBezier 复用既有模块;ApproxCurve 暂由 curve_approx 覆盖 |
| BRepExtrema（形状级） | brep_extrema.rs | 33%→60% | closest_point_on_edge/face 改走解析/牛顿 extrema,带边界钳制（含逆向边、边界边遍历）、平面求解器参数框架修正 |
| 备注 | intana vs intpatch | — | 平面∩球/球∩球/平面∩柱两处重叠（不同 API 层,互补非重复）,待决定是否合并 |

强项：底层几何内核 + 拓扑数据结构 + 交换（STEP/OBJ/STL/PLY/VRML）已形成可用闭环。
最大缺口：精确布尔 → 真实 BRepMesh → IGES/BinXCAF → Visualization。

## 语义偏差风险（与 OCCT"对齐"的差异点）

1. **体素布尔 / 采样求交结果不是精确 TopoDS** —— 离散网格结果与 OCCT 精确布尔在拓扑上不等价（已标注）。
2. `is_inside`/`point_in_mesh` 曾因射线双计数误判（已修复：重合命中参数去重）。
3. `make_box_corner` 曾只改 bbox 不改几何（已修复：经 prism 构建真实位置几何）。
4. 曲线/曲面 trait 对象无法下转型 —— 类型分类靠几何不变量采样（平面/球面可判，柱/锥近似）。
