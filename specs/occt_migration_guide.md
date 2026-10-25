# OCCT → Rust 移植指导说明（忠实度审计 + 下一轮路线图）

> 本文是画板「OCCT 移植忠实度审计与下一轮 10 万行路线图」的可执行版本：把画板上的审计结论、
> lite 替代件缺口、7 波次路线图、门禁与禁令整理成一份可直接照着做的迁移指导。
>
> - 画板：`~/.cursor/projects/d-source-repos-dogs/canvases/occt-migration-audit-and-100k-roadmap.canvas.tsx`
> - 审计条目台账：`specs/_audit/_index.md`（A0–A33）、分区报告 `specs/_audit/{occt-core,occt-geom,bop-boolean,topo-rest,mesh-exchange}.md`
> - 任务台账：`specs/_tasks.md`（T-204、T-209–T-213 等）
> - 板真相：`docs/goals/drift-and-gap-autofix/state.yaml`（`goal.oracle` 四道基线）
> - 门禁基线冻结记录：`docs/goals/drift-and-gap-autofix/notes/T002-gate-baseline.md`
>
> **本文所有数字为 2026-10-22 在本机实测**（PowerShell + `cargo`），不是引用画板旧数；
> 与画板不一致处在 §7 逐条订正。

---

## 0. 量测口径（先读这一节，后面每个数字都按此口径）

| 对象 | 口径 | 工具 |
|---|---|---|
| OCCT 侧源码根 | `D:\source\OCCT-src\src`（OCCT 8.0.0） | `Get-ChildItem -Recurse` |
| OCCT 行数 | `.cxx` + `.hxx`（计划口径）；另标 `.lxx` / `.gxx` / `GTests` | `Get-Content \| Measure-Object -Line` |
| Rust 侧 | `crates/{occt-core,occt-math,occt-geom,occt-geom2d,occt-topo}`，**排除 `target/`** | 同上 |
| Rust 行数 | 逐文件求和，**含 `tests/` 与 `examples/`** | 同上 |
| 类数 | 各包 `.hxx` 里 `class`/`struct` 的**去重类名**计数 | `Select-String` + `Sort-Object -Unique` |
| 忠实度标记 | 源码中 `UNPORTED` 出现次数与分布文件数 | `rg -c UNPORTED` |
| 门禁 | `cargo check` / `--lib` / `--tests` / `step_obj_gates` / `export_data_obj` | `cargo --offline` |

**行号口径**：本文引用源码用**物理行号**（`grep`/`rg` 的行号），不用 `read` 工具的折行显示行号。

**排除项**：OCCT 各包的 `GTests/`（C++ 测试脚手架）不是移植目标；`FILES.cmake` 也不计入类/文件体量。

---

## 1. 全景：Rust 侧体量、测试与门禁现状

### 1.1 五个 crate 的实测体量（2026-10-24）

| crate | `.rs` 文件 | 行数 | lib 单测 | `UNPORTED` 出现 | 分布文件 |
|---|---:|---:|---:|---:|---:|
| `occt-core` | 227 | 47,757 | 283 passed / 0 failed | 27 | 20 |
| `occt-math` | 68 | 16,553 | 210 passed / 0 failed | 4 | 3 |
| `occt-geom` | 323 | 87,844 | 143 passed / 0 failed | 77 | 38 |
| `occt-geom2d` | 75 | 20,442 | 72 passed / 0 failed | 10 | 7 |
| `occt-topo` | 1,130 | 323,455 | 1254 passed / 0 failed | 236 | 164 |
| **合计** | **1,823** | **496,051** | **1,962** | **354** | **232** |

- 全部 `UNPORTED` 标记**均带 OCCT 文件与行号出处**，无「假出处」。
- OCCT 全树（`.cxx/.hxx/.lxx/.gxx`）：**13,803 文件 / 2,439,487 行**。
- 本仓 Rust 行数 / OCCT 全树 ≈ **20.3%**；下一轮 7 波次再加 390,414 行 OCCT 对应实现。
- `UNPORTED` 较 2026-10-22（168 处 / 118 文件）增 185 处 / 113 文件，增量集中在
  `occt-topo/src/topopebrep/`（303 文件 / 73,399 行）、`approx_int/`（7 文件 / 2,212 行）、
  `brep_approx/`（10 文件 / 2,168 行）三处新落地的模块——属「先按分支标出来再补」，不是回退。

### 1.2 门禁基线（本轮实跑，2026-10-24）

| 门禁 | 命令 | 读数 |
|---|---|---|
| 编译 | `cargo check --manifest-path crates/<crate>/Cargo.toml --offline --all-targets`（5 个 crate） | **0 errors**（警告为既有基线数） |
| 单元测试 | `cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --lib` | **1254 passed / 0 failed**（27.9s） |
| 其余四 crate | `cargo test --manifest-path crates/<crate>/Cargo.toml --offline --lib` | core 283 / math 210 / geom 143 / geom2d 72，**全 0 failed** |
| 门禁套件 | `cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --tests` | **13 个套件、69 passed / 0 failed**（含 lib 共 14 套件 / 1323 passed） |
| STEP 对拍 | `cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --test step_obj_gates` | **4 passed**（23 个模型对 `data/occ-*.obj`，203.5s） |
| 整体导出 | `cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj` | **15 ok / 0 err** |

**13 个门禁套件**（`crates/occt-topo/tests/`，套件名沿用不改）：

| 套件 | 用例数 | 套件 | 用例数 |
|---|---:|---|---:|
| `boolean_mesh_exchange_gates` | 7 | `modeling_exchange_gates` | 9 |
| `bop_builder2_boss` | 2 | `modeling_gates` | 3 |
| `bopalgo_builder_gates` | 5 | `pave_filler_gates` | 5 |
| `curved_boolean_xcaf_gates` | 8 | `rwmesh_builder_viz_gates` | 5 |
| `features_healing_math_gates` | 5 | `step_obj_gates` | 4 |
| `geometry_pipeline_gates` | 4 | `intersection_offset_exchange_gates` | 5 |
| `massprops_draw_render_gates` | 7 | | |

**`export_data_obj` 与 OCCT 参考网格的逐模型顶点/面数对照**（`data/occ-*.obj` = OCCT 8.0.0 参考）：

| 模型 | OCCT v/f | 本仓 v/f | 模型 | OCCT v/f | 本仓 v/f |
|---|---|---|---|---|---|
| Cube | 24/12 | 24/12 | Shape | 6150/11372 | 6150/11372 |
| Cone | 195/301 | 195/301 | Shape-1 | 3343/4336 | 3343/4336 |
| Cylinder | 146/140 | 146/140 | Shape-2 | 3098/4778 | 3098/4778 |
| Extrusion | 24/12 | 24/12 | Sphere | 642/1244 | 642/1244 |
| HoledPlate | 180/128 | 180/128 | Torus | 1369/2592 | 1369/2592 |
| Offset | 712/892 | 712/892 | linkrods | 3494/5078 | 3494/5078 |
| OffsetPlaneHoleEdge | 48/32 | 48/32 | rev | 104/92 | 104/92 |
| screw | 600/790 | 600/790 | | | |

→ **15/15 逐位一致**（顶点/面数）；网格顶点包围盒抽查 `screw`、`Shape-1` 亦逐位一致。

`data/occ/` 下 8 个 OCCT 原始测试模型（`ATU01038 / bottom / motoc / top / T0M / TDB / acs10 / a3n00`）
与 OCCT 参考**存在已知残差**（顶点/面数差 ≤ 0.15%，`motoc` 完全一致），由 `step_obj_gates` 的
逐模型 bbox/面积容差覆盖；这些差值是**已登记的 STEP 读入残差**，不是本轮回归。

### 1.3 S 波次类覆盖快照（2026-10-24，口径：OCCT 类全名在 Rust 源码出现）

| 波次 | 包 | 类数 | 已落地 | 未落地行数 | 说明 |
|---|---|---:|---:|---:|---|
| S1 | GeomFill | 61 | 61 | 0 | 全收口 |
| S1 | Law | 11 | 11 | 0 | 全收口 |
| S2 | TopOpeBRep | 37 | 30 | 306 | 余 `define` + 5 个 `P*` + `traceSIFF`（死码） |
| S2 | TopOpeBRepBuild | 44 | 40 | 123 | 余 `define` + 3 个 `P*` |
| S2 | TopOpeBRepDS | 56 | 55 | 72 | 余 `define` 宏头 |
| S2 | TopOpeBRepTool | 37 | 31 | 350 | 余 `defineG`/`Plos`/`PSoClassif`/`AncestorsTool`/`KRO`+`STATE`（死码） |
| **S2 合计** | | **174** | **156** | **851** | 2026-10-23 时为 86/171、缺 38,451 行 |
| S2-IntPatch | IntPatch | 40 | 25 | 4,983 | 余 `RstInt` / `Polyhedron` 族 |

- 「未落地行数」= 该类 `.hxx` + `.cxx` 物理行合计。
- **S2 的算法体已全部落地**（T-278 之后）。剩余 851 行分类：`define`/`defineG` 宏头 195、
  `P*` 句柄别名 180、`OCCT_DEBUG` 死码 381（`traceSIFF` 176 + `STATE` 143 + `KRO` 62）、
  `AncestorsTool` 的 `TopExp` 单行转发 74、`Plos` 的 `typedef List*` 21。
  这几类在 Rust 侧由 `const` / `Arc` / `Vec<TopoShape>` / `topo_tools::map_shapes_and_ancestors`
  天然覆盖，不单独立类。

---

## 2. 忠实度判定口径（照抄仓库规则，不新增标准）

| 档 | 定义 | 处置 |
|---|---|---|
| **忠实移植** | 与 OCCT `.cxx` 逐行对应（分支、常量、控制流顺序一致） | 通过 |
| **等价替换** | 用本仓已有原语替代 OCCT 工具类；语义等价 | 必须在注释登记**翻译边界 + OCCT 出处** |
| **自创 / 未移植** | ① OCCT 无该规则/阈值/分支；② 以采样/网格/迭代近似替代 OCCT 解析或控制流**且改变结果**；③ 静默默认值掩盖失败；④ 产物根本不是从 OCCT 翻译来的（含**假出处**：文件头写 `Source: X` 而 X 无此物） | 必须摘除、改成忠实件，或**标 `UNPORTED` 并写清 OCCT 文件与行号** |

**忠实的判定标准是「对着 `.cxx` 逐分支核」，不是「跑出来数字对」。**
以前某次 Rust 跑出来的数字不是规格，OCCT 才是。

---

## 3. 已完成移植的审计现状（A0–A33）

条目编号沿用 `specs/_audit/_index.md`。下表是 2026-10-22 对**工作树现状**（而非文档结项时间）的复核。

| 编号 | 条目 | 现状 | 证据 / 落点 |
|---|---|---|---|
| A0 | Offset 曲线 `d0/d1/d2` 公式错误（平移 vs 法向） | ✅ 已修 | `occt-geom/src/offset.rs` 按 `Geom_OffsetCurveUtils.pxx` 重写；残留仅 4/5 阶 `EvalDN`（非初等/非 B 样条基曲线），已标 `UNPORTED` |
| A1 | `surface_closest_params` 网格扫描替代点-面极值 | ✅ 忠实件已接 | `occt-geom/src/extrema_surf/gen_ext_ps/`（`Extrema_GenExtPS` 逐段移植）；仅 `numeric_extrema.rs` 保留为**已下线的替代件**并标 `UNPORTED`（该目录现共 2 处 `UNPORTED`） |
| A2 | `is_inside` 7×7 采样网格 + 抖动射线 | ✅ 已修 | 实体走 `SolidClassifier::classify`；旧实现更名 `is_inside_mesh` 并标 `UNPORTED`（仅非实体形状） |
| A3 | `step::classify_curve` 6 点采样 + 0.02 阈值 | ✅ 已修 | 分类器整体删除，写侧改走 `GeomToStep_MakeCurve.cxx:50-104` 的 `IsKind` 顺序 |
| A4 | `AlgoTools::compute_state` 网格投影 | ✅ 已修 | F 分支走 `BRepClass_FaceClassifier::Perform` |
| A5 | 曲面布尔体素/网格装配兜底 | ✅ 已摘除 | `bop_curved/` 已整目录删除、`voxel_fallback` 已删；曲面分支走 `boolean_via_bopalgo`（**残留见 §7.3**：`solid_union.rs` 仍留体素布尔） |
| A6 | `gcpnts` 伪 `UniformDeflection` 等（递归二分 + 深度截断） | ✅ 已修 | 改走 `CurveTessellator` / `GCPnts_TangentialDeflection`；三个伪结构体已删 |
| A7 | 点-曲线/曲线-曲线网格+牛顿种子集 | ✅ 已修 | `Extrema_ExtPC` / `Extrema_ExtCC` / `GGenExtCC` 逐段移植，自创种子集与采样分类器全删 |
| A8 | `CompCurveToBSplineCurve` 采样重拟 | ◐ 主体已修 | `GeomConvert::CurveToBSplineCurve` 已接；残留：修剪 Bezier 的 rational 子支（`convert_bspl/curve_convert.rs:247-248`） |
| A9 | 7 处**假出处**（文件头 `Source:` 指向 OCCT 中不存在的物） | ✅ 已修（但仍有 1 处同类，见 §7.3） | `convert/`、`cslib/`、`gprop/`、`bnd/obb_pca.rs`、`bnd/intersect.rs`、`polyline_simplify.rs`、`int/curve_curve.rs` 文件头均已订正 |
| A10 | 静默默认值（球/环面零二阶导、缺「被包含」分支、凭空造交点、`choose_left_way` 只取首候选） | ✅ 已修 | `ElSLib::SphereD2/TorusD2`、`Bnd_Sphere.cxx:73-101`、`Poly_MakeLoops::chooseLeftWay` 均已补 |
| A11 | 布尔退化分支（质心判定 / 采样解并入精确解） | ✅ 已修 | `boolean_degenerate` 整套删除；`edge_edge` 采样解并入段删除 |
| A12 | `brepfeat` 解析覆盖体积 / 网格夹具冒充 `MakeCylinder` / 车削网格重建 / `clamp(16,64)` | ✅ 已修 | `boss_thru_all`/`boss` 改走解析 `translated_cylinder`（`BRepPrimAPI_MakeCylinder`），`revolve_profile_about` 改走 `brep_sweep::make_revol`（`BRepSweep_Revol`）；`mesh_cylinder` / `resolution_for` 已不存在（`rg` 零命中） |
| A13 | `face_to_triangles` 未裁剪 UV 窗口建网格 | ✅ 已修 | 改委托 `IncrementalMesh::discretize_face`（`BRepMesh_FaceDiscret` + `BRepMesh_Delaun`） |
| A14 | 整包 lite 替代件（借 OCCT 包名） | ◐ 已登记未移植 | `hlr.rs`(254) / `viz_scene/`(3,727) / `draw/`(2,562) / `xcaf/`(2,020) / `render_svg.rs`(355) / `validate.rs`(191) 模块头均已加 `UNPORTED` + 真实出处；**下一轮 S7 直接覆盖 `hlr`** |
| A15 | `IsDone` 被兜底吞掉 + 自创阈值（默认 `d2` 前向差分等） | ◐ 主体已修 | `occt-geom/src/surface.rs:25` 仍有 `let h = 1e-6;` 的默认差分（5 个初等面 + B 样条/修剪面/回转面/拉伸面均已覆写解析 `D2`，其余 3 个实现者继承近似，已逐个列名）；`gcpnts` 的 `UniformDeflection`/`QuasiUniformDeflection` 仍 `UNPORTED`；`adv_app2var` 若干臂 `UNPORTED` |
| A16 | `geom_api` 256×256 采样求交/投影器 | ✅ 已修 | 采样族整体删除；曲线–曲线走 `IntTools_EdgeEdge`，非平面边–面走 `IntCurveSurface_HInter` |
| A17 | 网格出口偏转形参失效 + 三级静默回退 | ✅ 已修 | 两级回退删除；形参改名 `maximal_chordial_deviation` |
| A18 | 平面面片耳切 + 最近点对桥洞 | ✅ 已修 | `brepmesh.rs` 整体删除（−602 行），`wireframe.rs` 平面机器删除（−288 行），改走约束 Delaunay |
| A19 | 失败率阈值整形状换算法（`WIREFRAME_FALLBACK_RATIO_MAX = 0.10`） | ✅ 已修 | 阈值与两级兜底全删（`rg WIREFRAME_FALLBACK_RATIO / MAX_FACE_DEPTH` 零命中）；逐面失败只置 `IMeshData_Failure` |
| A20 | STEP 读入边参数用距离代替 `IsSame` 拓扑判据 | ✅ 已修 | 两次 `Project`（`StepToTopoDS_TranslateEdge.cxx:442-444`） |
| A21 | `adjustSamePoints` 退化支路左右接反 | ✅ 已修 | `meshing/model_healer.rs` 按 `BRepMesh_ModelHealer.cxx:491-512` 两支逐行对齐 |
| A22 | `ProjectAct` 缺椭圆/抛物/双曲精确臂 | ✅ 已修 | `ShapeAnalysis_Curve.cxx:355-477` 分派齐全 |
| A23 | `face_uv_bounds` 9 点 pcurve 采样 | ✅ 已修 | 委托 `brep_tools::uv_bounds`（`BRepTools::UVBounds` → `AddUVBounds`） |
| A24 | 网格面型靠周期标志 + 半径采样猜 | ✅ 已修 | `classify_surface` 改为 `GeomAdaptor_Surface::Load`（`cxx:422-513`）精确类判定 |
| A25 | `Delaun::decomposeSimplePolygon` 自造删三角形补丁 | ✅ 已修 | 补丁删除；失败改为置 `IMeshData_Failure` |
| A26 | 交换/写侧偏离（门禁外） | ◐ 主体已修 | PLY 焊接、STL 退化法向、IGES 曲线族、VRML `solid FALSE` 已修；**残留**：STEP 写侧修剪 Bezier 的 `Segment` 臂（`step/write_context.rs:366-371`）、`RECTANGULAR_TRIMMED_SURFACE` 包装（`write_context.rs:51-54`）、pcurve 的 B 样条拟合臂（`write_context.rs:246-249`）仍 `UNPORTED` |
| A27 | 解析曲线/曲面 `continuity()` 一律返回 3 | ✅ 已修 | 15 处改为 `GeomAbs_CN`(6) |
| A28 | `BRepGProp_Face` 边界弧用 3D 曲线 + UV 反演 | ✅ 已修 | `BoundaryArc` 改 `ArcGeom::Pcurve`，优先取已存 pcurve |
| A29 | 直线解析臂建在 `d0(uinf)` 上导致参数整体平移 | ✅ 已修 | 改取 `Curve::gp_line()` |
| A30 | `brep_surface::classify_surface` 8×8 采样 | ✅ 已修 | 改精确类判定（顺带修好柱/锥 pcurve 单位换算被整段跳过） |
| A31 | 伪造回转面类 `GeomRevolvedSurface` | ✅ 已修 | `occt-geom/src/revolved.rs` 已整文件删除 |
| A32 | OBJ/STL/PLY 写出漏 `TriangleOriented`（REVERSED 面反绕） | ✅ 已修 | `face_occurrence_reversed` / `oriented_triangle` 复刻 `RWMesh_FaceIterator::TriangleOriented` |
| A33 | 共享 `ShapeBuild_ReShape` 被环境开关默认关闭 | ✅ 已修 | 开关删除，共享 `ReShape` 成为唯一路径 |

### 3.1 三个专项条目（非 A 编号）

| 条目 | 现状 | 证据 |
|---|---|---|
| **`BRepCheck_Edge`** | ◐ 部分移植 | `Minimum` 已移植（`brep_check.rs:154` `edge_check`）；`InContext` / `Tolerance` / `CheckPolygonOnTriangulation` 未移植（`brep_check.rs:153` 自承 `InContext` not ported）。`brep_check.rs` **无生产消费者、不被任何门禁覆盖** ⇒ 按规则不做无法验证的移植 |
| **`ShapeFix` 模式体系** | ✅ 已收口，◐ 残留 | T-209–T-213 已把 `ShapeFix_Wire` / `Face` / `Shape` / `Solid` / `Shell` 的模式门字段化（新增 `shhealing/wire_fix/modes.rs`、`shhealing/shape_fix_face/modes.rs`、`shhealing/shape_fix_face/small_area_wire.rs`、`shhealing/shape_fix_shape_modes.rs`），两轮 wire 循环、`AutoCorrectPrecision`、`FixLoopWire`、`FixSmallAreaWire`、`CreateOpenSolidMode` 均已接线。**残留**：`XSAlgo_ShapeProcessor.cxx:696-710` 的 `DE_ShapeFixParameters` 灌入层未移植（三个预设仍硬编码）；STEP 侧 `ShapeFix_Solid` 多壳 `CreateSolids` 臂（`cxx:577-640`）未接线（仅 IGES 路径可达）；`ShapeFix::FixVertexPosition` 未移植（两读入路径该 key 均为 `0`，不触发） |
| **T-204 STEP 写侧基线** | ⛔ BLOCKED（非待办） | 读侧基线齐全，**写侧无 OCCT 对照产物**：`specs/occt_probe/` 无 `STEPControl_Writer`/`STEPCAFControl_Writer` 探针，`data/` 下无 OCCT 写出的 STEP。只能自对自的往返断言被规则禁止。解封条件：① 先做 OCCT 侧写侧探针；或 ② 另开允许往返 harness 的批次 |

---

## 4. 真正的缺口：被 lite 替代件占位的 OCCT 包

下一轮不追新特性，而是把**已被 lite 替代件占位**的 OCCT 包按 `.cxx` 逐分支补齐。
左边是 OCCT 真实体量，右边是本仓现有替代件体量 —— 差距就是本轮工作量。

| OCCT 包 | OCCT 行数 | 本仓替代件 | 替代件行数 | 缺口性质 |
|---|---:|---|---:|---|
| TKFillet（全部 9 子包） | 85,017 | `fillet.rs` + `fillet_edge/` + `fillet_curved/` + `fillet_var/` | 7,972 | 仅平面-平面柱面 + 平面/球/柱环面特例；**无 ChFi3d 骨架、无可变半径律** |
| TKBool/BRepFill | 29,942 | `brep_fill.rs` + `loft.rs` + `brep_pipe.rs` | 889 | 直纹/折线扫掠；**无 SectionLaw / Sweep / ApproxSeewing** |
| TKOffset（BRepOffset + API + Draft + BiTgte） | 42,160 | `brep_offset/` | 1,166 | 面/壳偏移的简化臂；`MakeOffsetShape` 一般臂仍 `UNPORTED` |
| TKFeat（LocOpe + BRepFeat） | 28,455 | `brepfeat/` | 877 | 棱柱/旋转体特征；**无 LocOpe 局部布尔与 RibSlot** |
| TKHLR（HLRBRep + HLRAlgo + HLRTopoBRep + Contap） | 33,556 | `hlr.rs` | 254 | **画家算法替代，不是 HLRBRep 隐藏线**；模块头已自承 `UNPORTED` |
| TKGeomAlgo/GeomFill + Law | 31,812 | `occt-geom/src/geom_fill.rs` | 314 | 只有 `Profiler` + `Generator` 两个类；**无 SectionLaw / Sweep / 约束填充** |
| **合计** | **250,942** | | **11,472** | 替代件仅覆盖 OCCT 真实体量的 **4.6%** |

> 这 6 个包的 OCCT 行数与 §5 的 7 波次 390,414 行互有重叠：S1 = GeomFill + Law；
> S3 = TKFillet；S4 = BRepFill + LocOpe + BRepFeat；S6 = TKOffset；S7 = TKHLR。

---

## 5. 下一轮路线图：7 波次 / 390,414 行 OCCT

### 5.1 波次定义与依赖顺序

| 波次 | 范围 | 包 | `.cxx` | `.hxx` | 合计 | 累计 | OCCT 类 | `.cxx+.hxx` 文件 |
|---|---|---|---:|---:|---:|---:|---:|---:|
| **S1** | GeomFill + Law 曲面填充基座 | GeomFill, Law | 24,831 | 6,981 | **31,812** | 31,812 | 90 | 140 |
| **S2** | TopOpeBRep 经典布尔引擎 | TopOpeBRep, TopOpeBRepBuild, TopOpeBRepDS, TopOpeBRepTool | 87,094 | 12,051 | **99,145** | 130,957 | 213 | 375 |
| **S3** | TKFillet 完整圆角与变半径 | ChFi3d, BRepBlend, BlendFunc, ChFi2d, ChFiKPart, ChFiDS, Blend, FilletSurf, BRepFilletAPI | 72,127 | 12,890 | **85,017** | 215,974 | 195 | 243 |
| **S4** | BRepFill / LocOpe / BRepFeat 拓扑放样与局部特征 | BRepFill, LocOpe, BRepFeat | 52,616 | 5,781 | **58,397** | 274,371 | 125 | 139 |
| **S5** | Geom2dGcc + GccAna + GeomPlate + FairCurve | Geom2dGcc, GccAna, GeomPlate, FairCurve | 33,228 | 7,099 | **40,327** | 314,698 | 100 | 156 |
| **S6** | TKOffset 完整三维/二维偏移与 BiTgte | BRepOffset, BRepOffsetAPI, Draft, BiTgte | 38,089 | 4,071 | **42,160** | 356,858 | 88 | 72 |
| **S7** | TKHLR 隐藏线消除管线 | HLRBRep, HLRAlgo, HLRTopoBRep, Contap | 24,299 | 9,257 | **33,556** | **390,414** | 200 | 181 |
| | | **合计** | **332,284** | **58,130** | **390,414** | | **1,011** | **1,306** |

- **10 万行下限**：S2 结束时累计 **130,957 行**，越过门槛（3.9× 于下限）。
- 文件数口径：仅 `.cxx` + `.hxx`；各包另有 `.lxx`（98 个文件）与 `GTests`（不移植）。

### 5.2 依赖 DAG

```text
        S1 (GeomFill+Law) ─────────┬──> S3 (TKFillet) <── S2 (TopOpeBRep)
                                  │
                                  └──> S4 (BRepFill/LocOpe/BRepFeat)

        S5 (Gcc+Plate+FairCurve)   S6 (TKOffset)   S7 (TKHLR)   ← 无波内前置，可并行
```

- **S1 与 S2 是仅有的两个前置基座**；S5 / S6 / S7 波内无前置，可与 S1/S2 并行开工。
- S3 同时依赖 S1（`ChFi3d` 用 `GeomFill_ConstrainedFilling`）与 S2（`ChFi3d` 用 `TopOpeBRepDS` 干涉），
  因此 **S3 必须排在 S1、S2 之后**，不允许「先做简化圆角再补骨架」。

### 5.3 每波的前置与落点（开工前先确认落点目录是否已存在）

| 波次 | 前置 | 移植落点（现状：多数为空目录，需新建） |
|---|---|---|
| S1 | 无（仓内已有 `GeomFill_Generator`/`Profiler` 雏形 = `occt-geom/src/geom_fill.rs`） | `occt-geom/src/geom_fill/`：`GeomFill_SectionLaw` / `Sweep` / `ShapeLaw` / `ConstrainedFilling`；`occt-math/src/law/`：`Law_Function` 族 |
| S2 | 无（`BOPAlgo` 已在仓内；本波是其经典引擎底座） | `occt-topo/src/topopebrep/`：DS 数据结构 / 干涉 / 布尔构建 / 工具层 |
| S3 | **S1 + S2** | `occt-topo/src/chfi3d/` + `occt-topo/src/blend/`：`ChFi3d` 骨架 + `BRepBlend` 曲面族 + `BlendFunc` 段形 + 变半径律 |
| S4 | S1 | `occt-topo/src/brep_fill/`、`occt-topo/src/locope/`、`occt-topo/src/brepfeat/`：脊线放样 / 壳转角修剪 / LocOpe 局部布尔 / BRepFeat 棱肋槽 |
| S5 | 无（`GccEnt`/`GccInt` 为包内轻件） | `occt-geom2d/src/gcc/`：圆-圆/圆-线通用裁剪；`occt-geom/src/plate/` + `occt-geom/src/fair_curve/`：板面拟合与光顺曲线 |
| S6 | 无（`BRepTools_Quilt`/`Modification` 已在 `brep_sweep/` 内） | `occt-topo/src/brep_offset/` 一般臂、`occt-topo/src/draft/`、`occt-topo/src/bi_tgte/` |
| S7 | 无（`Intf`/`Intf_Polygon2d` 已在 `occt-core/src/intf`） | `occt-topo/src/hlr/`：替换 `hlr.rs` 画家算法，落地 `HLRBRep_Algo` + `HLRToShape` |

**已核实**：S2/S3/S4/S5/S6/S7 的落点目录（`topopebrep/`、`chfi3d/`、`blend/`、`locope/`、
`gcc/`、`plate/`、`fair_curve/`、`draft/`、`bi_tgte/`、`hlr/`）**当前都不存在**；
S1 的 `occt-geom/src/law/` 存在但**未接线**（见 §7.3）。

### 5.4 逐包体量（29 个包）

| 波次 | OCCT 包 | `.cxx` | `.hxx` | 合计 | 类数 |
|---|---|---:|---:|---:|---:|
| S1 | TKGeomAlgo/GeomFill | 21,316 | 5,776 | 27,092 | 78 |
| S1 | TKGeomAlgo/Law | 3,515 | 1,205 | 4,720 | 12 |
| S2 | TKBool/TopOpeBRep | 17,507 | 2,461 | 19,968 | 55 |
| S2 | TKBool/TopOpeBRepBuild | 33,788 | 3,445 | 37,233 | 60 |
| S2 | TKBool/TopOpeBRepDS | 18,382 | 3,526 | 21,908 | 55 |
| S2 | TKBool/TopOpeBRepTool | 17,417 | 2,619 | 20,036 | 43 |
| S3 | TKFillet/ChFi3d | 34,411 | 2,075 | 36,486 | 33 |
| S3 | TKFillet/BRepBlend | 12,767 | 3,280 | 16,047 | 52 |
| S3 | TKFillet/BlendFunc | 10,065 | 2,239 | 12,304 | 31 |
| S3 | TKFillet/ChFi2d | 4,233 | 786 | 5,019 | 10 |
| S3 | TKFillet/ChFiKPart | 4,938 | 542 | 5,480 | 5 |
| S3 | TKFillet/ChFiDS | 2,926 | 1,216 | 4,142 | 28 |
| S3 | TKFillet/Blend | 764 | 1,416 | 2,180 | 15 |
| S3 | TKFillet/FilletSurf | 1,172 | 407 | 1,579 | 10 |
| S3 | TKFillet/BRepFilletAPI | 851 | 929 | 1,780 | 11 |
| S4 | TKBool/BRepFill | 27,264 | 2,678 | 29,942 | 62 |
| S4 | TKFeat/LocOpe | 12,077 | 1,463 | 13,540 | 34 |
| S4 | TKFeat/BRepFeat | 13,275 | 1,640 | 14,915 | 29 |
| S5 | TKGeomAlgo/Geom2dGcc | 11,712 | 3,268 | 14,980 | 38 |
| S5 | TKGeomAlgo/GccAna | 13,832 | 1,873 | 15,705 | 21 |
| S5 | TKGeomAlgo/GeomPlate | 4,943 | 1,009 | 5,952 | 25 |
| S5 | TKGeomAlgo/FairCurve | 2,741 | 949 | 3,690 | 16 |
| S6 | TKOffset/BRepOffset | 27,068 | 1,483 | 28,551 | 24 |
| S6 | TKOffset/BRepOffsetAPI | 5,147 | 1,765 | 6,912 | 28 |
| S6 | TKOffset/Draft | 3,119 | 407 | 3,526 | 15 |
| S6 | TKOffset/BiTgte | 2,755 | 416 | 3,171 | 21 |
| S7 | TKHLR/HLRBRep | 14,565 | 5,648 | 20,213 | 111 |
| S7 | TKHLR/HLRAlgo | 3,814 | 1,503 | 5,317 | 33 |
| S7 | TKHLR/HLRTopoBRep | 1,746 | 359 | 2,105 | 15 |
| S7 | TKHLR/Contap | 4,174 | 1,747 | 5,921 | 41 |

---

## 6. 每波必须过的门禁与禁令

### 6.1 门禁（沿用现有基线，**不新写测试**）

| 阶段 | 动作 | 通过判据 |
|---|---|---|
| 改完先编译 | `cargo check --manifest-path crates/<crate>/Cargo.toml --offline --all-targets` | 0 errors（警告数不增） |
| 改完跑单测 | `cargo test --manifest-path crates/<crate>/Cargo.toml --offline --lib` | 0 failed，测试数不减少 |
| 动了共享管线 | `cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --tests` | 13 套件 69 passed / 0 failed |
| 动了 BOP / 网格 / STEP 读入 | `cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --test step_obj_gates` | 4 passed（23 模型） |
| 整体网格复查 | `cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj` | 15 ok / 0 err，v/f 与 `data/occ-*.obj` 逐位一致 |

> `step_obj_gates` 单次约 3 分钟，**不要每个小包都跑**；按 `state.yaml` 的 cadence
> 「每 3–5 个包或每次触碰共享管线时跑一次」。
> 全量 `cargo test` 只在提交前跑一次；改单一 crate 时只跑该 crate。

### 6.2 移植禁令（仓库规则，零例外）

1. **不得为某一个 STEP / 特征 / 特例调参。**
   网格管线是共享的，特例补丁会带偏已对齐的模型。
   缺分支时先定位缺了哪段 OCCT 控制流，按同一套 `.cxx` 分支修复；确实没有同等分支，
   就标 `UNPORTED` 并在注释写清 **OCCT 文件与行号**。
2. **不得引入 OCCT 里不存在的规则。**
   不加谓词、启发式、面积比、长度滤边、体积门这类 OCCT 里不存在的规则。
3. **不得用假输入证明已对齐。**
   不用 `StubBopBuilder` 一类假输入，不用自对自往返断言当基线。
4. **不得为对齐新写单元测试。**
   门禁只能引用既有基线（既有测试名、harness / parity 输出）。
5. **不得静默降级。**
   失败要可观测（`Option`/`Result`/`IsDone` 语义），不换算法兜底。
6. **完成即停。**
   当前任务门禁满足 → 报告结果与发现的旁支问题 → 停，等下一步指令。
   旁支问题写进报告，**不要动手修**，不要自创下一步。

### 6.3 代码组织约束

- 单文件超 800 行必须拆成模块目录：`<module>/mod.rs` + `mod <主题>;` + `pub use <主题>::*;`，
  调用方看到的仍是 `<module>::<item>`。**拆分为 move-only，不得夹带行为改动**。
- 禁止 `pNN.rs` 这类无语义文件名；文件名取该文件承载的 OCCT 类/职责。
- 注释与字符串一律 ASCII 英文；文档可中文。

---

## 7. 本次整理对画板数据的订正（2026-10-22 实测）

### 7.1 数字订正

| 项 | 画板值 | 实测值 | 说明 |
|---|---:|---:|---|
| Rust 文件数 | 1,333 | **1,294** | 画板含 `target/` 或按旧口径统计 |
| Rust 行数 | 360,914 | **361,344** | +430 行（工作树未提交改动） |
| `occt-topo` 文件 / 行 | 796 / 237,200 | **778 / 237,101** | 同上 |
| `occt-topo` 测试 | 1323 | **1254（lib）+ 69（门禁）= 1323** | 画板口径一致 ✓ |
| `UNPORTED` 总数 | 169 处 / 119 文件 | **168 处 / 118 文件** | 按 `rg -c` 精确计数 |
| `occt-topo` `UNPORTED` | 99 | 99 处 / **71 文件**（画板写 99 处未给文件数） | ✓ |
| S3 `.cxx+.hxx` 文件数 | 247 | **243** | 画板口径含 `.lxx`/`FILES.cmake`；本文按 `.cxx+.hxx` 计 |
| 各波次 `.cxx`/`.hxx`/类数 | 画板值 | **逐包实测全部一致** ✓ | 29 个包的行数与类数零偏差 |
| 门禁 `--lib` | 1254 ✓ / 2031 合计 ✓ | 全部复现 ✓ | |

### 7.2 审计条目状态订正

| 条目 | 画板描述 | 2026-10-22 实测 |
|---|---|---|
| A1 | 「部分移植；`extrema_surf` 仍留 2 处 UNPORTED」 | **2 处仍在**（`numeric_extrema.rs`、`mod.rs`），但那是**已下线的替代件**；忠实件 `gen_ext_ps/`（`Extrema_GenExtPS`）已在 live 路径上 ⇒ 实质已修 |
| A13-余 | 「`bop_builder_report::component_volume` 仍用网格测量体积」 | **已过时**。`bop_builder_report.rs:415` → `brep_gprop::volume` → `brep_gprop::shape_mass_properties`，**解析 Gauss 积分优先**，仅在解析积分失败时回退网格。残留的是**回退粒度差异**：OCCT 按面做 `NoSurf`/`NoTri`/surface 三分（`BRepGProp.cxx:344-360`），本仓按形状回退（`brep_gprop.rs:9-15` 已自承） |
| A26 | 「PLY 焊接 / STL 阈值 / IGES 采样 / VRML TRUE / OBJ 空 vn」 | 前四项**已修**（T-62/批 16/批 18）。残留：STEP 写侧修剪 Bezier `Segment`、`RECTANGULAR_TRIMMED_SURFACE` 包装、pcurve B 样条拟合臂 |
| BRepCheck | 「`InContext` / `Tolerance` / `CheckPolygonOnTriangulation`」 | 三项仍未移植；`Minimum` 已移植。`brep_check.rs` 无消费者、不被门禁覆盖 ⇒ 按规则不做无法验证的移植 |
| ShapeFix | 「`CreateSolids` 臂、`DE_ShapeFixParameters` 灌入层未接线」 | **仍成立**，且 T-209–T-213 已把模式字段化落地；残留面缩小为：`DE_ShapeFixParameters` 灌入层 + STEP 侧多壳 `CreateSolids` |

### 7.3 新发现的旁支问题（**本轮只记录，不动手修**）

1. **S1 的 Law + GeomFill 已进编译，但还没有消费者。**
   `crates/occt-geom/src/law/` 已从 T-196 时的 2 个文件补到 7 个（`function` / `bspline` /
   `interpolate` / `knot_split` / `composite` / `services`），`lib.rs:91` 已挂 `pub mod law;`。
   `crates/occt-geom/src/geom_fill/` 现在覆盖 `Filling` / `Stretch` / `Coons` / `Curved` /
   `BSplineCurves` / `BezierCurves` / `Profiler` / `Generator` / `Boundary`（`SimpleBound`、
   `DegeneratedBound`、`BoundWithSurf`）/ `CoonsAlgPatch` / `CornerState` / `TgtField`
   （`TgtOnCoons`）/ `ConstrainedFilling`，`lib.rs:88` 已挂 `pub mod geom_fill;`，
   `occt-core` 侧补了 `AdvApprox_ApproxAFunction` 的 `Num3DSS = N` 形态
   （`ApproxAFunction3dMulti`）。⇒ 这一波只做「翻译 + 编译」，**尚无调用方**：
   `BRepFill_Filling` / `ChFi3d`（S3）接线前，几何基线（`step_obj_parity`、
   `fuse_box_cylinder_is_closed_solid`）都不会经过它。
   `GeomFill_ConstrainedFilling` 的 4 个 `Check*` 诊断臂（`cxx:1465-1710`，只往 `std::cout`
   打印并调用 `Law_draw1dcurve`）按「无绘图」登记为未移植。
2. **`solid_union.rs` 的文件头是 A9 同类的「假出处」。**
   `crates/occt-topo/src/solid_union.rs:1-10` 声明 `Source: BRepAlgoAPI_Fuse / BRepAlgoAPI_BooleanOperation`，
   但文件实现的是三角汤并集 + **体素布尔近似**（`voxel_union_grid` / `voxel_union`），OCCT
   `BRepAlgoAPI_Fuse` 无对应分支（该文件自承 `ponytail:` 近似）。按 A9 口径应改为
   「非 OCCT 翻译」或标 `UNPORTED` + 出处。
3. **`intcurvesurface/types.rs:6` 的 `ANALYTIC_TOL = 1e-9` 与 OCCT 差一个数量级。**
   OCCT 同用途常量为 `THE_TOLTANGENCY = 1e-8`（`IntCurveSurface_InterUtils.pxx:56`）。
   T-212 / T-213 两次登记为旁支，仍未改。

---

## 8. 风险与处置

| 风险 | 影响 | 处置 |
|---|---|---|
| S3 的 `ChFi3d` 依赖 `TopOpeBRepDS` 干涉与 `GeomFill_ConstrainedFilling` | 高 | 已排在 S1、S2 之后；**不允许先做「简化圆角」再补骨架** |
| 移植后输出与冻结基线不一致（面数/顶点/面积） | 高 | 每波必须跑 `step_obj_gates` + `export_data_obj`；**失配即停**并报告差异，不得调基线 |
| 为凑某个 STEP 的结果加特例补丁 | 高 | 规则明令禁止：只走 `.cxx` 同等分支；缺分支就标 `UNPORTED` 并写明 OCCT 文件与行号 |
| 单文件超 800 / 1000 行 | 中 | 按 OCCT 类职责拆模块目录（`<module>/mod.rs` + `pub use`），**move-only，不改逻辑** |
| S7 隐藏线管线体量大且无 OCCT 几何基线 | 中 | 先移植 `HLRBRep_ShapeToHLR` + 轮廓提取；隐藏面判定对照现有 `specs/occt_probe/`；`hlr.rs` 画家算法在替换完成前保持 `UNPORTED` 状态不动 |
| IGES 曲线族采样与整球自造回转面等写侧偏离（A26） | 中 | 属 S4/S6 之外的独立收口项，按 `GeomToIGES_*` 实体逐条对齐，**不与本轮混做** |
| `DE_ShapeFixParameters` 灌入层缺失导致模式预设硬编码 | 中 | S1–S7 都不触碰 ShapeHealing；如需接线，先做 `XSAlgo_ShapeProcessor.cxx:696-710` 单点移植 |
| `intcurvesurface` 容差与 OCCT 差一个数量级 | 中 | 属数值级旁支，与本轮 7 波次无关；单独立项、单独对拍 |

---

## 9. 开工清单（把本文当作业指导用）

按波次逐个推进，每个包一个提交，每包一次窄门禁。

**S1（GeomFill + Law，31,812 行）**

- [x] `crates/occt-geom/src/law/` 补齐 `bspline/`（`bsp_func.rs` / `eval.rs` / `modify.rs` / `query.rs`）/ `composite.rs` / `interpolate.rs` / `knot_split.rs` / `services.rs` / `function.rs`
- [x] `crates/occt-geom/src/lib.rs` 加 `pub mod law;`
- [x] `occt-geom/src/geom_fill.rs` → `geom_fill/` 目录，扩 `SectionLaw` / `Sweep` / `ShapeLaw` / `ConstrainedFilling`
- [x] `Law_Function` 族落在 `crates/occt-geom/src/law/function.rs`（**未**放 `occt-math`：`Law_BSpline` 依赖 `Geom_BSplineCurve` 的求值栈，放在 `occt-math` 会造成反向依赖）
- [x] `AppDef` / `AppParCurves` 基础层：`app_def/`（`MultiPointConstraint` / `MultiLine` / `MyLineTool`）、`app_par_curves/`（`Constraint` / `BernsteinMatrix` / `SplineFunction`）、`occt-core/src/approx/`（`Approx_Status` / `Approx_ParametrizationType`）
- [x] `GeomFill_GordonBuilder`（`geom_fill/gordon/builder.rs`）及其依赖：`PLib::EvalLagrange`、`GeomAPI_Interpolate`（`geom_api_interpolate/`）、`Geom_BSplineSurface::IncreaseDegree` / `InsertUKnots` / `InsertVKnots` / `ExchangeUV` / `SetPole` / `SetWeight`、`Geom_BSplineCurve::IsEqual`
- [x] `math_Householder`（`occt-math/src/householder.rs`）+ `math_Recipes::DACTCL_Decompose`/`DACTCL_Solve`（`math_recipes.rs`）+ `AppParCurves_LeastSquare`（`app_par_curves/least_square/`，`AppParCurves_LeastSquare.gxx:1-1860`）
- [x] `math_Crout`（`occt-math/src/math_crout.rs`，`math_Crout.cxx:31-208`）+ `math_Uzawa`（`occt-math/src/uzawa.rs`，`math_Uzawa.cxx:39-285`）+ `AppParCurves_ResolConstraint`（`app_par_curves/resol_constraint/`，`AppParCurves_ResolConstraint.gxx:40-946`）
- [x] `AppParCurves_Function`（`app_par_curves/function/`，`AppParCurves_Function.gxx:35-695`）+ `AppParCurves_BSpFunction`（`app_par_curves/bsp_function.rs`，`AppParCurves_BSpFunction.gxx:27-338`）
- [x] `AppParCurves_Gradient`（`app_par_curves/gradient.rs`，`AppParCurves_Gradient.gxx:43-217`）+ `AppParCurves_BSpGradient`（`app_par_curves/bsp_gradient.rs`，`AppParCurves_BSpGradient.gxx:82-410`）；两者依赖的 `BSplCLib::PolesCoefficients` / `CoefsD1` 已在 `occt-core/src/bspl/coefs.rs`
- [x] `Approx_BSplComputeLine`（`approx/bspl_compute_line/`，`Approx_BSplComputeLine.gxx:1-1458`）→ `AppBlend_AppSurf::InternalPerform` 的 `SpApprox` 臂（`app_surf/perform.rs`，`AppBlend_AppSurf.gxx:420-471`）→ 收口 `NSections::ComputeSurface` 的无参考面分支
- [ ] `AppBlend_AppSurf::PerformSmoothing`（`AppDef_Variational.cxx:1-1413`，`AppBlend_AppSurf.gxx:472-541`）仍 UNPORTED；`GeomFill_Sweep.cxx:288-327`（`GeomConvert_ApproxSurface`）仍 UNPORTED。T-227 已补 `PLib_HermitJacobi`（`occt-math/src/plib_hermit_jacobi.rs`），`FEmTool` 链（`FEmTool_Curve` / `ProfileMatrix` / `LinearTension|Flexion|Jerk` / `Assembly`）→ `AppDef_LinearCriteria` → `AppDef_Variational` 是收口这条的下一段
- [x] `GeomFill_Gordon`（`geom_fill/gordon/`，`GeomFill_Gordon.cxx` + `GeomFill_GordonBuilder.cxx`）及其依赖：`Approx_BSplineApproxInterp`（`approx/bspline_approx_interp.rs`）、`GeomAPI_Interpolate`（`geom_api_interpolate/`）、`extrema_cc::curve_curve_extrema_all`、`GeomFill_Profiler`
- [x] `GeomFill_NSections`（`geom_fill/n_sections/`，`GeomFill_NSections.cxx:201-1047`）+ `pipe::init` 的 `Init(Path, FirstSect, LastSect)` / `Init(Path, NSections)` 接线（`GeomFill_Pipe.cxx:606-728`）
- [x] 门禁：`cargo check` + `occt-geom --lib`（143 passed）+ `occt-math --lib`（210 passed）
- [x] **S1 类覆盖收口（2026-10-23 复核）**：`GeomFill` 60/60 + `Law` 10/10 = **70/70 类**（口径：OCCT 类全名，下划线可作分隔符，在 Rust 源码中出现即已落地）。`geom_fill/` 87 文件 / 23,955 行 + `law/` 11 文件 / 3,300 行。余 `AppBlend_AppSurf::PerformSmoothing` 与 `GeomFill_Sweep` 的 `GeomConvert_ApproxSurface` 臂仍 UNPORTED（见上）。

**S2（TopOpeBRep，99,145 行）** — 类覆盖 **94/171**（2026-10-24 复核）

- [x] 新建 `crates/occt-topo/src/topopebrep/`（`ds/`、`tool/` 两个子目录已落地；`interf/`、`build/` 待建）
- [x] `TopOpeBRepDS` 基座（T-218）：`kind` / `transition` / `geometry_data` / `interference` / `shape_data` / `point(_data)` / `surface(_data)` / `curve(_data)` / `shape_with_state` / `list_of_shape_on1_state` / `data_structure/`（四文件，`cxx:33-1355`）/ `curve_explorer` / `surface_explorer` / `point_explorer`
- [x] `TopOpeBRepTool` 定向族（T-219）：`ShapesSameOriented` / `FacesSameOriented` / `SurfacesSameOriented` / `EdgesSameOriented` / `CurvesSameOriented` + `EdgeData`（`TopOpeBRepTool_ShapeTool.cxx:394-653`）→ `tool/shape_tool.rs`；被读到的 `BRepAdaptor_Surface` / `BRepAdaptor_Curve` 切片 → `tool/adaptor.rs`；`LProp_CurveUtils` → `tool/lprop.rs`；`FUN_tool_projPonS` → `tool/project.rs`；`GeomAbs_CurveType` → `geom_abs_curve_type.rs`
- [x] `TopOpeBRepDS` 的 `FillShapesSameDomain` 两个重载（T-219，`TopOpeBRepDS_DataStructure.cxx:828-991`）
- [x] `TopOpeBRepDS_TKI`（T-220，`TopOpeBRepDS_TKI.cxx:1-432`）→ `ds/tki.rs`；`TopOpeBRepDS_HDataStructure`（`HDataStructure.cxx:40-813`）→ `ds/h_data_structure/`（查询 + `MakeCurve`/`SortOnParameter`/`ScanInterfList` 与 `StoreInterference*` 分两文件）
- [x] `TopOpeBRepDS_InterferenceTool`（T-220，`InterferenceTool.cxx:29-143`）→ `ds/interference_tool.rs`；`TopOpeBRepDS_ProcessInterferencesTool` 的纯列表过滤族（`ProcessInterferencesTool.cxx:38-391`、`501-528`）→ `ds/process_interferences_tool/mod.rs`
- [x] `TopOpeBRepDS_Marker` / `TopOpeBRepDS_Association`（T-220，`Marker.cxx:24-93`、`Association.cxx:22-108`）→ `ds/marker.rs` / `ds/association.rs`；`TopOpeBRepDS_EXPORT` 纯访问器与 same-domain 查询（`EXPORT.cxx:61-223`、`226-274`、`464-495`、`647-712`）→ `ds/export.rs`
- [x] `TopOpeBRepTool_TOOL` 闭合/定向族（T-221，`TOOL.cxx:142-472`）→ `tool/closed.rs`；切线族 `TggeomE` / `TgINSIDE` / `FUN_tool_dirC`（`TOOL.cxx:619-717`、`:1039-1048`）→ `tool/tangent.rs`；`FUN_tool_bounds` / `projPonC` / `projPonE` / `value` / `parE`（`PROJECT.cxx:33-291`、`TOPOLOGY.cxx:313-340`、`:932-1048`）→ `tool/project_curve.rs`；`FUN_tool_maxtol` / `nbshapes`（`TOPOLOGY.cxx:1148-1190`）→ `tool/measures.rs`；`FUN_tool_ngS` / `onapex` / `quadCT` / `quad` / `closed`（`GEOMETRY.cxx:86-358`）→ `tool/geometry.rs`
- [x] `TopOpeBRepDS_TOOL`（T-221，`TOOL.cxx:82-404`）→ `ds/tool.rs`；`FUN_reducedoublons`（`ProcessInterferencesTool.cxx:394-497`）→ `ds/process_interferences_tool/mod.rs`；`BRep_Tool::UVPoints(E,F)`（`BRep_Tool.cxx:1005-1106`）→ `brep_tool.rs`；`FDS_aresamdom` 两个重载（`EXPORT.cxx:282-357`）→ `ds/export.rs`
- [x] `TopOpeBRepDS_connex`（T-222，`connex.cxx:41-296`）→ `ds/connex.rs`（`thread_local` 全局改为显式 `Connexity` 结构）；`FDS_EdgeIsConnexToSameDomainFaces`（`EXPORT.cxx:357-391`）→ `ds/export.rs`
- [x] `TopOpeBRepDS_FIR`（T-223，`FIR.cxx:36-571`，含 `FUN_isPonF` / `FUN_findPonF` / `FUN_reduceEDGEgeometry1` / `FUN_GmapS` / `FUN_stateedgeface` / `FUN_reduceEDGEgeometry`）→ `ds/fir.rs`；`TopTrans_CurveTransition`（`TopTrans_CurveTransition.cxx:26-441`）→ `top_trans_curve_transition.rs`；`TopOpeBRepDS_EdgeInterferenceTool`（`EdgeInterferenceTool.cxx:30-174`）→ `ds/edge_interference_tool.rs`；`TopOpeBRepDS_FaceInterferenceTool` + `GeomLProp_SLProps` + `TopTrans_SurfaceTransition`（同批）→ `ds/face_interference_tool.rs` / `occt-geom/src/geom_lprop/sl_props.rs` / `top_trans_surface_transition.rs`
- [x] `TopOpeBRepDS_EIR` + `TopOpeBRepDS_Reducer` + `ProcessInterferencesTool` 1d/2d/3d 选择族（T-225，`EIR.cxx:1-1221`、`Reducer.cxx:22-48`、`ProcessInterferencesTool.cxx:531-1185`）→ `ds/eir/` / `ds/reducer.rs` / `ds/process_interferences_tool/{select_2d,select_3d}.rs`；`FUN_tool_Eshared`（`TOPOLOGY.cxx:843-883`）→ `tool/topology.rs`
- [x] `TopOpeBRepDS_Check` + `HDataStructure::ChkIntg` 接线；`TopOpeBRepTool_OutCurveType` / `GeomTool`（T-226，`Check.cxx:20-656`、`HDataStructure.cxx:133-139`、`OutCurveType.hxx:20-25`、`GeomTool.cxx:22-141`）→ `ds/check.rs` / `tool/out_curve_type.rs` / `tool/geom_tool.rs`
- [x] `TopOpeBRep` `interf` 基础值类型（T-227，`P2Dstatus.hxx:22-28`、`Bipoint.cxx:20-43`、`WPointInter.cxx:20-82`、`Hctxff2d.cxx:28-245`）→ `interf/{p2d_status,bipoint,w_point_inter,hctxff2d}.rs`；`BRepAdaptor_Surface::Face()`（`BRepAdaptor_Surface.hxx:44`）→ `tool/adaptor.rs` 的 `SurfaceAdaptor::face`
- [x] `GeomLib_CheckBSplineCurve` / `GeomLib_Check2dBSplineCurve`（T-227，`CheckBSplineCurve.cxx:24-188`、`Check2dBSplineCurve.cxx:26-187`）→ `occt-geom/src/geom_lib_check_bspline.rs`；`PLib_HermitJacobi`（T-227，`PLib_HermitJacobi.cxx:24-372`）→ `occt-math/src/plib_hermit_jacobi.rs`（`FEmTool` 前置）
- [x] `TopOpeBRep` `interf` 值类型补全（T-228）：`Point2d`（`Point2d.hxx:40-144` / `.lxx` / `.cxx:26-176`）、`Hctxee2d`（`Hctxee2d.cxx:44-251`）、`PointGeomTool`（`PointGeomTool.cxx:29-98`）、`PointClassifier`（`PointClassifier.cxx:26-73`）→ `interf/{point2d,hctxee2d,point_geom_tool,point_classifier}.rs`；`IntSurf_Transition` / `TypeTrans` / `Situation`（`IntSurf_Transition.cxx:17-44`、`.lxx:16-74`）→ `intsurf_transition.rs`；`BRepAdaptor_Surface::Initialize(F)`（`BRepAdaptor_Surface.cxx:56-82`）→ `SurfaceAdaptor::initialize`
- [x] `TopOpeBRep` `interf` 值类型补全（T-229）：完整 `IntPatch_Point`（`IntPatch_Point.hxx:33-263` / `.lxx:16-303` / `.cxx:26-261`）→ `intpatch_point.rs`；`BRepAdaptor_Curve2d`（`BRepAdaptor_Curve2d.cxx:19-77`）→ `brep_adaptor_curve2d.rs`；`BRepTopAdaptor_HVertex`（`BRepTopAdaptor_HVertex.cxx:17-200`）→ `brep_top_adaptor_hvertex.rs`；`TopOpeBRep_VPointInter`（`VPointInter.hxx:32-194` / `.lxx:21-303` / `.cxx:36-483`）→ `interf/v_point_inter.rs`；`BRep_Tool::Parameter(V,E,F)`（`BRep_Tool.cxx:1519-1624`）→ `BRepTool::parameter_of_vertex_on_face`；`TopOpeBRep_TypeLineCurve`（`TypeLineCurve.hxx:19-30`）→ `topopebrep/type_line_curve.rs`；`TopOpeBRep_LineInter`（`LineInter.hxx:42-204` / `.lxx:19-303` / `.cxx:19-700`）→ `topopebrep/line_inter.rs`
- [x] `TopOpeBRep` `interf` 余项（T-230）：`VPointInterIterator`（`VPointInterIterator.cxx:23-137`）/ `WPointInterIterator`（`WPointInterIterator.cxx:23-82`）→ `topopebrep/{v_point_inter_iterator,w_point_inter_iterator}.rs`；`TopOpeBRep_FFTransitionTool`（`FFTransitionTool.cxx:31-435`）→ `topopebrep/ff_transition_tool.rs`；`TopOpeBRep_GeomTool`（`GeomTool.cxx:31-188`）→ `topopebrep/geom_tool.rs`；`TopOpeBRepTool_CurveTool::MakeBSpline1fromPnt(2d)`（`CurveTool.cxx:898-972`）→ `topopebrep/tool/curve_tool.rs`；`Extrema_ExtPC` 定区间最小解 → `occt_geom::extrema_pc::point_curve_extrema_range`
- [x] `TopOpeBRep` 扫描 + 面-边相交（T-231）：`TopOpeBRepTool_HBoxTool`（`HBoxTool.cxx:51-260`）→ `tool/h_box_tool.rs`；`TopOpeBRepTool_BoxSort`（`BoxSort.cxx:53-356`）→ `tool/box_sort.rs`；`TopOpeBRepTool_ShapeExplorer`（`ShapeExplorer.hxx:24-77`）→ `tool/shape_explorer.rs`；`TopOpeBRep_ShapeScanner`（`ShapeScanner.cxx:37-151`）→ `topopebrep/shape_scanner.rs`；`Bnd_Box::IsOut(gp_Lin)`（`Bnd_Box.cxx:761-887`）→ `BndBox::is_out_line`；`BRepIntCurveSurface_Inter`（`BRepIntCurveSurface_Inter.cxx:30-262`）→ `brep_int_curve_surface_inter.rs`；`TopOpeBRep_FaceEdgeIntersector`（`FaceEdgeIntersector.cxx:26-352`）→ `topopebrep/face_edge_intersector.rs`
- [x] `TopOpeBRep` 面分类前置（T-232）：`BRepClass_FaceExplorer`（`FaceExplorer.cxx:26-307`）→ `brep_class/face_explorer.rs`；`BRepClass_FClassifier` + `TopClass_FaceClassifier`（`FClassifier.cxx:25-100`、`TopClass_FaceClassifier.pxx:28-176`）→ `brep_class/f_classifier.rs`；`TopOpeBRep_VPointInterClassifier`（`VPointInterClassifier.cxx:26-290`）→ `interf/v_point_inter_classifier.rs`；`TopOpeBRepTool_tol`（`tol.cxx:41-311`）→ `tool/tol.rs`
- [x] `BRepTopAdaptor_TopolTool`（T-233，`BRepTopAdaptor_TopolTool.cxx:43-653`）→ `brep_top_adaptor_topol_tool.rs`
- [x] `IntPatch_Line` 基态 + `IntPatch_LineConstructor` 前缀（T-234）：`IntPatch_Line`（`IntPatch_Line.hxx:35-137` / `.lxx:17-94` / `.cxx:21-65`）→ `intpatch_line.rs`；`GeomIntLine` 的 transition/situation 入口与 `TopOpeBRep_LineInter::{TransitionOnS1,TransitionOnS2,SituationS1,SituationS2,ComputeFaceFaceTransition,Bounds}`（`LineInter.cxx:404-429`、`:592-598`）→ `topopebrep/line_inter.rs`；`IntPatch_LineConstructor` 静态函数族（`LineConstructor.cxx:44-589`、`:896-946`）+ `AddLine`（`:592-880`）→ `intpatch_line_ctor/{helpers,add_line}.rs`
- [x] 四类交线顶点参数归一（T-235）：`IntPatch_WLine::ComputeVertexParameters` + `RecadreMemePeriode` / `CompareVertexAndPoint` / `SameVtxRst` / `CompareVerticesOnS1/2`（`IntPatch_WLine.cxx:88-321`、`:323-1211`）→ `intpatch_vparams/{common,wline}.rs`；`IntPatch_RLine`（`IntPatch_RLine.cxx:143-433`）→ `intpatch_vparams/rline.rs`；`IntPatch_GLine`（`IntPatch_GLine.cxx:414-1102`）→ `intpatch_vparams/gline.rs`；`IntPatch_ALine`（`IntPatch_ALine.cxx:70-676`）→ `intpatch_vparams/aline.rs`；`IntPatch_Point` 的 `SetArc` / `SetVertex` / `IsVertexOnS1/2` / `ParameterOnArc1/2` → `int_tools_wline/patch_point.rs`；`IntPatch_ImpImpIntersection.cxx:2977-2995` 的按类型分派 → `intpatch_impimp_put.rs`；`LineConstructor.cxx:688` / `:737` 接线 → `intpatch_line_ctor/add_line.rs`
- [x] `TopOpeBRep_FacesFiller` 的退化边族（T-240）：`FUN_VPIndex`（`TopOpeBRep_vpr.cxx:350-469`）→ `faces_filler/vp_index.rs`；`FUN_FillVof12` / `FUN_addmapve` / `FUN_GetdgData`（`TopOpeBRep_vprdeg.cxx:105-379`）→ `faces_filler/degen_data.rs`；`FUN_putInterfonDegenEd` / `ProcessVPondgE` / `local_FindVertex` / `local_FindTreatedEdgeOnVertex` / `local_ReduceMapOfTreatedVertices`（`vprdeg.cxx:381-862`）→ `faces_filler/vprdeg_interf.rs`；`TopOpeBRepTool_mkTondgE`（`mkTondgE.cxx:56-530`：`Initialize` `:58-101`、`SetclE` `:105-108`、`IsT2d` `:113-116`、`SetRest` `:121-132`、`MkTonE` `:272-329` / `:334-530`、`FUN_getEc` `:192-215`、`FUN_MkTonE` `:217-267`）→ `tool/mk_tondg_e.rs`；`TopOpeBRepTool_TOOL::Matter`（`TOOL.cxx:1338-1373`）→ `tool/matter.rs`；`UVISO` / `ParISO`（`TOOL.cxx:1191-1243`、`:544-576`）→ `tool/iso.rs`；`TolP`（`TOOL.cxx:1935-1939`）与 `FUN_tool_inS`（`TOPOLOGY.cxx:832-840`）→ `tool/topology.rs`
- [x] `TopOpeBRep_FacesFiller::ProcessVPonR` + `FUNBREP_topokpart` + `FUN_newtransEdge`（T-241）：`ProcessVPonR`（`TopOpeBRep_vpr.cxx:877-1732`）→ `faces_filler/vponr.rs` + `faces_filler/vponr_ds.rs`；`FUN_onedge`（`vpr.cxx:848-866`）、`FUN_brep_ONfirstP`（`_1.cxx:634-644`）→ `vponr.rs`；`FUN_LineRestF` / `FUN_newtransEdge` / `FUN_ScanInterfList` / `FUN_selectTRAISHAinterference` / `FUN_selectGinterference` / `FUN_sameGsameS`（`vpr.cxx:471-780`）→ `faces_filler/vpr_trans_edge.rs`；`FUN_processCPI`（`vpr.cxx:782-853`）→ `faces_filler/vpr_cpi.rs`；`FUNBREP_Periodize` / `FUNBREP_HasSameGPoint` / `FUNBREP_SameUV` / `FUNBREP_topowalki` / `FUNBREP_topogline` / `TopoParameter` / `FUNBREP_topokpart`（`TopOpeBRep_kpart.cxx:49-916`）→ `topopebrep/kpart.rs`；`FDS_getupperlower`（`TopOpeBRepDS_EXPORT.cxx:617-643`）→ `ds/export.rs`
- [x] `TopOpeBRep_FacesFiller::ProcessRLine`（T-242，`TopOpeBRep_FacesFiller_1.cxx:647-1125`）→ `faces_filler/rline.rs`；`ProcessVPondgE` 的 `PVKind`/`PVIndex` 改为 in/out 入参（`TopOpeBRep_vprdeg.cxx:636-658`）→ `faces_filler/vprdeg_interf.rs`
- [x] `TopOpeBRep_EdgesFiller`（T-243，`EdgesFiller.hxx:39-110` / `.cxx:60-648`：构造 `:60-66`、`Insert` `:72-363`、`SetShapeTransition` `:372-425`、`GetGeometry` `:430-440`、`MakeGeometry` `:444-473`、`Face(ISI,F)` `:477-491`、`Face(ISI)` `:495-509`、`StorePI` `:514-528`、`StoreVI` `:551-565`、`ToRecompute` `:581-590`、`StoreRecompute` `:594-603`、`RecomputeInterferences` `:608-648`）→ `topopebrep/edges_filler.rs`
- [x] `TopOpeBRepDS_BuildTool` 的 `Make*` 族与边更新族（T-256）：`MakeEdge` 四个重载（`BuildTool.cxx:157-263`）、`MakeFace`（`:272-277`）→ `ds/build_tool/make_shapes.rs`；`GetOrientedEdgeVertices`（`:320-339`）、`UpdateEdgeCurveTol`（`:343-468`）、`UpdateEdge`（`:1006-1028`）、`Project`（`:1031-1055`）、`Parameter(C,E,V)`（`:1058-1065`）、`PCurve(F,E,CDS,PC)`（`:1243-1291`）→ `ds/build_tool/edge_update.rs`；`ds/build_tool.rs` 拆成目录模块（`mod.rs` 436 行），`BRep_TVertex::UpdateTolerance`（`BRep_TVertex.cxx:29-35`）→ `GeometryRegistry::raise_vertex_tolerance`，详见 `_tasks.md` T-256
- [ ] `TopOpeBRepDS` 余项：`TopOpeBRepDS_BuildTool`（`BuildTool.cxx` 行数按批统计，`Parameter` 已由 T-250 落地、`Make*`/边更新族已由 T-256 落地、`GetGeomTool`/`ChangeGeomTool`/`Approximation` + `ComputePCurves` `:654-789` / `PutPCurves` `:793-813` / `RecomputeCurves` `:816-890` 由 T-275 落地、`ApproxCurves` `:471-653` 由 T-277 第 5 批落地；`RecomputeCurve` `:1344-1376` / `RecomputeBSpline1Curve` `:1377-1453` / `RecomputeCurveOnCone` `:1454-1580` 在 OCCT 8.0.0 里整体处于 `cxx:1342` 起的 `/* ... */` 注释块内，属死码不译 —— 余项只有 `EXPORT` 的 `FDS_SIisGIofIofSBAofTofI` / `FDS_getupperlower`（缺 `TOOL::NgApp` / `Nt` / `XX` / `Getduv`））；`CurveTool` 侧 `MakeBSpline1fromPnt`/`MakeBSpline1fromPnt2d` 由 T-230、`IsProjectable`/`MakePCurveOnFace`/`MakePCurve` 与 `ProjLib_ProjectedCurve` 的解析臂由 T-274 落地，`MakePCurveOnFace` 的 `BRepTools::UVBounds` 口径与 `ShapeTool::Resolution3d` 由 T-275 落地，`BRepApprox` 的评价层（`ApproxInt_SvSurfaces` / `ApproxInt_MultiLine` / `BRepApprox_ApproxLine` / `BRepApprox_SurfaceTool` / `BRepApprox_ThePrmPrmSvSurfacesOfApprox`）由 T-277 第 1-4 批落地，`CurveTool::MakeCurves` 与 `BuildTool::ApproxCurves` 由第 5 批落地（第 1 批评价层 `SvSurfaces`/`MultiLine`/`ApproxLine`/`SurfaceTool`/`PrmPrmSvSurfaces`；第 2 批 `ApproxInt_MultiLineTool` 落成 `LineTool` 并把 `Approx_ComputeLine`/`Approx_BSplComputeLine`/`AppParCurves` 全部泛化到 `&L: LineTool`、`ApproxInt_KnotTools`；第 3 批 `ApproxInt_Approx` 驱动 + `BRepApprox_Approx` 的面型分发；第 4 批 `ApproxInt_ImpPrmSvSurfaces` + `IntImp_ZerImpFunc` + `SurfaceTool` 的 `Plane`/`Cylinder`/`Sphere`/`Cone`/`Torus`/`Is*Periodic`/`*Period`/`*Resolution` + `gxx:246-307` quadric 臂；第 5 批 `TopOpeBRepTool_CurveTool::MakeCurves` `CurveTool.cxx:300-895` + `CurveImprovement` 抽点 + `BuildTool::ApproxCurves` + `RecomputeCurves` 的 `iswalk && approx` 臂；UNPORTED 清单见 `_tasks.md` T-274 (3)、T-275 (3)、T-277 (3)(7)(8)(9)(10)）
- [ ] `TopOpeBRepTool` 余项（CORRISO / REGUW / RegularizeW / CLASSI / connexity / face / C2DF / PURGE 已由 T-264 落地）：`FuseEdges`（`FuseEdges.cxx:882`）、`makeTransition`（`:797`）、`AncestorsTool`（`:28`）、`Plos`、`PShapeClassifier`、`PSoClassif`、`PurgeInternalEdges`、`STATE`
- [x] `TopOpeBRep_FaceEdgeFiller`（T-244，`FaceEdgeFiller.hxx:30-70` / `.cxx:36-316`：构造 `:36`、`Insert` `:40-169`、`ScanInterfList` `:187-203`、`GetGeometry(IT,FEINT,G,BDS)` `:206-220`、`StoreInterference` `:248-286`、`MakeGeometry` `:289-316`）→ `topopebrep/face_edge_filler.rs`；`TopOpeBRep_PointGeomTool::MakePoint(const TopOpeBRep_FaceEdgeIntersector&)`（`PointGeomTool.cxx:78-82`）→ `interf/point_geom_tool.rs`
- [x] `TopOpeBRep_DSFiller`（T-246，`DSFiller.hxx:44-160` / `.cxx:1215`：构造 `:137-141`、`PShapeClassifier` `:155-158`、`BREP_correctgbound` `:174-278`、`BREP_UnfillSameDomain` `:282-300`、`FUN_shareNOG` `:305-424`、`Insert` `:427-439`、`InsertIntersection` `:444-719`、`Complete` `:723-738`、`CompleteDS` `:754-788`、`Filter` `:792-802`、`Reducer` `:806-813`、`RemoveUnsharedGeometry` `:817-878`、`Checker` `:882-888`、`Insert2d` `:892-905`、`InsertIntersection2d` `:907-996`、`CompleteDS2d` `:1001-1024`、`IsMadeOf1d` `:1026-1066`、`IsContext1d` `:1072-1079`、`Insert1d` `:1084-1133`、`ClearShapeSameDomain` `:1145-1179`）→ `topopebrep/ds_filler/{mod,insert,complete,ds_tools}.rs`；`BREP_sortonparameter`（`TopOpeBRep_sort.cxx:26-104`）与 `BREP_mergePDS`（`TopOpeBRep_mergePDS.cxx:29-196`）→ `ds_filler/ds_tools.rs`；`Interference::is_curve_point`（`CurvePointInterference.hxx:28-45` 的 RTTI 判据）→ `ds/interference.rs`；`myPShapeClassifier` 改为 `Rc<RefCell<ShapeClassifier>>` 共享（`FacesFiller::set_pshape_classifier` 同步）
- [x] `TopOpeBRep_ShapeIntersector2d`（T-245，`ShapeIntersector2d.hxx:47-120` / `.cxx:28-411`：构造 `:28-32`、`Reset` `:38-48`、`Init` `:52-57`、`SetIntersectionDone` `:62-65`、`CurrentGeomShape` `:69-97`、`InitIntersection` `:102-106`、`MoreIntersection` `:109-125`、`NextIntersection` `:187-217`、`Init/Find FF` `:226-280`、`NextFFCouple` `:285-298`、`Init/Find EEFF` `:306-367`、`NextEEFFCouple` `:381-388`、`Shape` `:392-405`、`ChangeEdgesIntersector` `:409-411`）→ `topopebrep/shape_intersector2d.rs`
- [ ] `TopOpeBRep` 余项（29/36 类已落地）：`TopOpeBRep_FacesFiller` / `FacesFiller_1`、`PEdgesIntersector` / `PFacesFiller` / `PFacesIntersector`、`PIntRes2d_IntersectionPoint`、`PPntOn2S`、`FFDumper` / `traceSIFF`
- [x] `TopOpeBRepDS_GapTool` + `TopOpeBRepDS_GapFiller`（T-247，`GapTool.hxx:34-96` / `.cxx:275`：`Init` `:52-84`、`Curve` `:89-115`、`EdgeSupport` `:144-175`、`FacesSupport` `:179-191`、`ParameterOnEdge` `:193-227`、`SetParameterOnEdge` `:230-259`、`SetPoint` `:263-275`；`GapFiller.hxx:30-106` / `.cxx:537`：`Perform` `:52-84`、`FindAssociatedPoints` `:103-199`、`FilterByFace` `:239-258`、`FilterByEdge` `:287-303`、`Normal` `:330-394`、`FilterByIncidentDistance` `:398-462`、`ReBuildGeom` `:465-514`、`BuildNewGeometries` `:517-537`）→ `ds/gap_tool.rs` + `ds/gap_filler.rs`；`DsFiller` 持 `Connexity` 并接通 `FDSCNX_Prepare` / `GapFiller`
- [x] `TopOpeBRepDS_Filter`（T-248，`Filter.cxx:22-93`：构造 `:22-28`、`ProcessInterferences` `:32-36`、`ProcessEdgeInterferences()` `:39-57`、`ProcessFaceInterferences(MEsp)` `:61-80`、`ProcessCurveInterferences()` `:85-93`；`FilterEdgeInterferences.cxx:46-89`；`FilterFaceInterferences.cxx:37-68`；`FilterCurveInterferences.cxx:21-64`）→ `ds/filter/mod.rs`；`TopOpeBRepDS_ProcessEdgeInterferences.cxx`（685 行）的 `FUN_vertexofedge` `:35-53`、`FUN_keepEinterference` `:55-131`、`FUN_unkeepEinterferences` `:134-167`、`FUN_unkeepEsymetrictransitions` `:245-350`、`FUN_orderFFsamedomain` `:352-405`、`FUN_orderSTATETRANSonG` `:408-453`、`FUN_orderSTATETRANS` `:455-486`、`FUN_resolveEUNKNOWN` `:488-571`、`FUN_purgeDSonSE` `:575-685` 与 `TopOpeBRepDS_ProcessFaceInterferences.cxx`（434 行）的 `FUN_mkTonF` `:69-143`、`FUN_edgeofface` `:146-163`、`FUN_keepFinterference` `:165-186`、`FUN_unkeepFdoubleGBoundinterferences` `:189-277`、`FUN_resolveFUNKNOWN` `:280-431` → `ds/filter/{edge_tools,face_tools}.rs`；`TopOpeBRepDS_funk.cxx:44-279` → `ds/filter/funk.rs`；前置 `BRepTools::UVBounds` / `AddUVBounds`（`BRepTools.cxx:64-367`）→ `brep_tools_uv.rs`；`BRepExtrema_ExtPF`（`BRepExtrema_ExtPF.cxx:31-100`）→ `brep_extrema_pf.rs`；`FUN_tool_projPonboundedF`（`TopOpeBRepTool_PROJECT.cxx:304-347`）→ `tool/project.rs`；`DsFiller::filter` / `complete_ds_2d` 接通
- [x] `TopOpeBRepDS_EXPORT` 的 `CompleteDS` 十三段 + `TopOpeBRepDS_redu`（T-249，`TopOpeBRepDS_EXPORT.cxx`：`FUN_ds_hasSDMancestorfaces` `:715-758`、`FUN_ds_PURGEforE9` `:761-859`、`FUN_ds_completeforSE1` `:862-1003`、`FUN_ds_completeforSE2` `:1005-1350`、`FUN_ds_completeforSE3` `:1352-1421`、`FUN_ds_shareG` `:1423-1522`、`FUN_ds_mkTonFsdm` `:1523-1668`、`FUN_ds_oriEinF` `:1676-1791`、`FUN_ds_hasI2d` `:1997-2051`、`FUN_ds_completeforSE4` `:2055-2153`、`FUN_ds_completeforSE5` `:2157-2333`、`FUN_ds_completeforSE6` `:2336-2503`、`FUN_ds_completeforE7` `:2504-2682`、`FUN_ds_completeforSE8` `:2684-2778`、`FUN_ds_completeforSE9` `:2779-2866`、`FUN_ds_PointToVertex` `:2868-3039`、`FUN_redusamshaonE` `:3041-3107`、`FUN_ds_redusamsha` `:3109-3277`、`FUN_ds_complete1dForSESDM` `:3777-3952`；`TopOpeBRepDS_redu.cxx`：`FUN_scanloi` `:35-54`、`FUN_ds_redu2d1d` `:56-179`、`FUN_ds_GetTr` `:181-353`）→ `ds/complete_ds/`（14 个文件）+ `ds/redu.rs`；`DsFiller::complete_ds` 按 `DSFiller.cxx:770-784` 的调用序接通 `run_complete_ds`
- [ ] `TopOpeBRepDS_EXPORT` 余项（依赖 `TopOpeBRepBuild` 包，非本波次）：`FUN_ds_hasFEI`（`:3279-3306`，调用方 `TopOpeBRepBuild_GridFF.cxx:227/246/281`；`FUN_ds_ONesd` `:3308-3345` 已随 T-257 落地）；`FUN_ds_unkeepEVIonGb1`（`:3364-3425`）与 `FUN_ds_FEIGb1TO0`（`:3647-3774`）在 OCCT 内调用点本身被注释（`TopOpeBRep_DSFiller.cxx:768`、`TopOpeBRepDS_Reducer.cxx:44`），属死码不译
- [x] `TopOpeBRepBuild` 首批 10/43 类（T-250）：`LoopEnum`（`LoopEnum.hxx:18-27`）、`AreaBuilder`（`AreaBuilder.cxx:30-430`）、`Area1d/2d/3dBuilder`、`FaceAreaBuilder`、`SolidAreaBuilder`、`FaceBuilder`（`FaceBuilder.cxx:58-661`）、`EdgeBuilder`（`.cxx:22-82`）、`SolidBuilder`（`.cxx:22-147`）、`VertexInfo`（`.cxx:16-153`）→ `topopebrep/build/`；`TopOpeBRepDS_BuildTool::Parameter` + `BRep_Builder::UpdateVertex`（`BRep_Builder.cxx:1220-1311`）→ `ds/build_tool.rs`；`TopOpeBRepBuild_Builder` 簿记层（`Builder.cxx:77-924`、`Merge.cxx:174-758`、`Section.cxx:248-254`）→ `build/builder/mod.rs`，详见 `_tasks.md` T-250
- [x] `TopOpeBRepBuild_Builder` 的 split / fill 族（T-251）：`AddIntersectionEdges`（`Builder.cxx:150-179`）、`SplitEdge`/`SplitEdge1`/`SplitEdge2`（`:925-1153`）、`SplitFace`/`SplitFace1`/`SplitFace2`（`:1155-1533`）、`SplitSolid`（`:1535-1711`）、`SplitShapes` + `FUN_touched`（`:1713-1874`）→ `build/builder/splits.rs`；`FillShape`/`FillFace`/`FillSolid`（`:1876-1989`）、`FillVertexSet`/`FillVertexSetOnValue`（`:1991-2057`）→ `build/builder/fill.rs`，详见 `_tasks.md` T-251
- [x] `TopOpeBRepBuild_Builder` 的 `Perform` / `MergeShapes` / `Section` 三面（T-252）：`Perform(HDS)`（`Builder.cxx:116-133`）、`Perform(HDS,S1,S2)`（`:138-145`）→ `build/builder/perform.rs`；`BuildVertices`（`BuildVertices.cxx:25-34`）→ 同文件；`IsKPart`（`KPart.cxx:169-172`）、`MergeShapes`（`Merge.cxx:174-336`）、`MergeSolids`（`:356-361`）、`MergeSolid`（`:368-371`）→ `build/builder/merge.rs`；`Section` / `SectionCurves` / `SectionEdges` / `FillSecEdgeAncestorMap`（`Section.cxx:1234-1392`）→ `build/builder/section.rs`；`Builder::myShapeClassifier` 改为共享单元以对上 `Filter` 的裸指针语义（`Builder.cxx:131`），详见 `_tasks.md` T-252
- [x] `TopOpeBRepBuild_HBuilder` 公开 API（T-253）：`TopOpeBRepBuild_HBuilder`（`HBuilder.cxx:39-738`：`Perform` `:50-62`、`Clear` `:67-71`、`Merge*` `:97-119`、`IsSplit`/`Splits`/`IsMerged`/`Merged` `:124-155`、`Section` 游标 `:203-296`、`MakeEdgeAncestorMap` `:301-416`、`GetDSEdgeFromSectEdge` `:420-446`、`GetDSFaceFromDSEdge` `:449-477`、`MakeCurveAncestorMap` `:480-504`、`GetDSFaceFromDSCurve` `:526-549`、`GetDSPointFromNewVertex` `:553-577`、`EdgeCurveAncestors` `:580-607`、`EdgeSectionAncestors` `:612-714`、`IsKPart`/`MergeKPart`/`ChangeBuilder` `:718-738`）→ `build/h_builder.rs`；`TopOpeBRepTool_TOOL::SplitE`（`TopOpeBRepTool_TOOL.cxx:473-528`）+ `FUN_tool_sortVonE`（`:436-469`）+ `TrslUVModifE`（`:1319-1335`）+ `FUN_ds_CopyEdge`/`FUN_ds_Parameter`（`TopOpeBRepTool_TOPOLOGY.cxx:1573-1623`）→ `tool/split_e.rs`；`Builder::shape_of_key` + `SectionCurves(LSE)` / `SectionEdges(LSE)` 的外部列表形态（`Section.cxx:1255-1268`、`:1271-1351`），详见 `_tasks.md` T-253
- [x] `TopOpeBRepBuild` 的 `G*` 驱动族（T-254）：`GToMerge`（`Grid.cxx:95-120`）、`GTakeCommonOfSame`（`:125-147`）、`GTakeCommonOfDiff`（`:152-174`）、`GFindSamDom` 两个重载（`:182-190`、`:197-256`）、`GFindSamDomSODO` 两个重载（`:260-268`、`:272-340`）、`GCopyList` 两个重载（`:441-464`）、`GParamOnReference`（`:1055-1116`）→ `build/builder/grid.rs`；`GFillEdgePVS`（`GridEE.cxx:198-214`）、`FUN_MoreSHAINT`（`:218-246`）、`FUN_isonbound`（`:305-319`）、`GFillPointTopologyPVS` 两个重载（`:331-776`、`:778-1051`）、`GLOBAL_issp`/`GLOBAL_faces2d`/`GLOBAL_hassd` → `build/builder/grid_ee.rs`；`GEDBUMakeEdges`（`makeedges.cxx:79-228`）→ `build/builder/make_edges.rs`；`TopOpeBRepBuild_FUN_aresamegeom`（`GridFF.cxx:132-141`）→ `build/builder/grid.rs`，详见 `_tasks.md` T-254
- [x] `TopOpeBRepBuild_Builder` 的 section-edge 拆分（T-255）：`SplitSectionEdges`（`Section.cxx:304-915`）+ `FUN_periodicS`（`:85-114`）、`FUN_periodic`（`:116-123`）、`FUN_onboundsper` 两个重载（`:125-176`）、`FUN_PinC`/`FUN_PinE`/`FUN_EstaEE`（`:178-244`）、`FUN_selLEE`（`:257-277`）、`FUN_FindEinSLOS`（`:282-300`）→ `build/builder/section_edge.rs`；`SplitSectionEdge`（`:928-1230`）→ `build/builder/section_edge_split.rs`；接线 `Perform` `Builder.cxx:130`、`MergeShapes` `Merge.cxx:200`、`SectionEdges` `Section.cxx:1294` 三处，详见 `_tasks.md` T-255
- [x] `TopOpeBRepBuild_Builder::BuildFaces`（T-256）：`BuildFaces(iS, HDS)`（`BuildFaces.cxx:40-96`）+ `BuildFaces(HDS)`（`:99-107`）→ `build/builder/build_faces.rs`，接线 `Perform` `Builder.cxx:127`（依赖 `BuildTool::MakeFace` 与四参 `PCurve`，均由 T-256 落地），详见 `_tasks.md` T-256
- [x] `TopOpeBRepBuild_BuilderON` 前段（T-257）：`GFillONPartsWES`（`on.cxx:28-49`）、`Perform`（`BuilderON.cxx:158-188`）、`GFillONCheckI`（`:190-249`）、`GFillONPartsWES1`（`:251-290`）、`GFillONPartsWES2` 的 `yap00`/`yap0`/`yapc1`/`yapc2`/`yapc3` 五臂（`:554-1603`）→ `build/builder/on/{mod,helpers,wes2,wes2_yap0,wes2_yapc}.rs`；`FUN_tool_staPinE`（`topology.cxx:340-367`）→ `tool/topology.rs`；`FUN_ds_ONesd`（`EXPORT.cxx:3308-3345`）→ `ds/export.rs`，详见 `_tasks.md` T-257
- [x] `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` 余段（T-271）：`FCX` 块（`BuilderON.cxx:1605-1787`）与 `yap1`（`:1790-2071`）/`yap2`（`:2074-2439`）/`yap6`（`:2442-2808`）/`yap3`（`:2811-2914`）/`yap5`（`:2917-3010`）/`yap4`（`:3013-3043`）/`yap7`（`:3046-3121`）→ `build/builder/on/{wes2_fcx,wes2_yap2,wes2_yap6,wes2_yap3,wes2_yap5,wes2_yap4,wes2_yap7}.rs`，`Wes2Ctx` 增 `factor`（`cxx:652`），详见 `_tasks.md` T-271
- [x] `TopOpeBRepBuild_Builder` 的 KPart 五分支（T-258）：`FindIsKPart`（`KPart.cxx:111-166`）、`MergeKPart`（`:176-208`）、`KPisdisj`/`KPisfafa`/`KPissoso`/`KPiskoletge` 与 `*sh` 子谓词（`:1195-2194`）、`KPisdisjanalyse`/`KPiskoletgeanalyse`（`:2010-2194`）、`MergeKPartisdisj`（`:840-1030`）+ `BuildNewSolid`（`:2486-2655`）+ `disjPerform*`（`:2660-2978`）、`MergeKPartisfafa`/`issoso`（`:1034-1192`）、`MergeKPartiskoletge`（`:578-836`）+ `FUN_makefaces`/`FUN_rebuildfc`/`FUN_comparekoletgesh`/`FUN_changev`/`FUN_updatev`/`FUN_addf`（`:268-570`）、`KPlhg`/`KPlhsd`/`KPls`/`KPContains`/`KPclasSS`/`KPclassF`/`KPclassFF`（`:2196-2367`）→ `build/builder/k_part/{mod,lists,predicates,merge_same,merge_disj,merge_koletge}.rs`；`FUN_tool_comparebndkole`（`TOPOLOGY.cxx:1218-1282`）、`FUN_tool_pcurveonF`（`:1347-1367`）、`FUN_tool_inS`（`:832-840`）、`FUN_tool_SameOri`（`:1285-1336`）→ `tool/topology.rs` + `k_part/merge_koletge.rs`；接线 `Perform(HDS,S1,S2)` `Builder.cxx:144` 与 `MergeShapes` `Merge.cxx:203-208`（`GMergeFaces`/`GMergeSolids` 仍未译），详见 `_tasks.md` T-258
- [x] `TopOpeBRepBuild_Builder` 的 `iskole` 分支与 `End()`（T-259）：`KPiskole`（`kpkole.cxx:576-732`）、`KPiskoleanalyse`（`:734-843`）、`KPiskolesh`（`:1085-1157`）→ `build/builder/k_part/iskole.rs`；`MergeKPartiskole`（`:56-572`）、`KPmakeface`（`:945-1031`）、`FUNKP_KPmakefaces`（`:845-941`）→ `k_part/merge_iskole.rs`；`TopOpeBRepBuild_WireToFace`（`WireToFace.cxx:23-61`）→ `build/wire_to_face.rs`；`End()` 块 1 连续性重编码（`End.cxx:47-193`）+ 块 2 顶点容差更新（`:195-365`）→ `build/builder/end.rs`；接线 `FindIsKPart` 的 `iskole` arm（`KPart.cxx:139-145`）、`MergeKPart` 五臂 + `End()`（`:185-208`）、`KPClearMaps`（`:1453-1456`）（`End` 块 3 的 `TopOpeBRepTool_FuseEdges` `:366-410` 与 `UpdateSplitAndMerged` `:413-520` 由 T-272 落地），详见 `_tasks.md` T-259
- [x] `TopOpeBRepBuild_Tools`（T-260）：状态传播族 `FindState`（`Tools.cxx:60-92`）、`PropagateState`（`:94-252`）、`FindState2`（`:254-290`）、`FindState1`（`:292-333`）、`FindStateThroughVertex`（`:335-387`）、`SpreadStateToChild`（`:389-410`）、`PropagateStateForWires`（`:412-457`）→ `build/tools/state.rs`；`GetNormalToFaceOnEdge`（`:461-484`）、`GetNormalInNearestPoint`（`:486-541`）、`GetTangentToEdgeEdge`（`:543-607`）、`GetTangentToEdge`（`:609-627`）、`GetAdjacentFace`（`:629-660`）、`IsDegEdgesTheSame`（`:838-853`）→ `build/tools/geom.rs`；`UpdatePCurves`（`:664-677`）、`UpdateEdgeOnFace`（`:679-714`）、`UpdateEdgeOnPeriodicalFace`（`:716-836`）、`NormalizeFace`（`:855-906`）→ `build/tools/pcurve.rs`；`ShapeClassifier::StateShapeReference(S, null)`（`TopOpeBRepTool_ShapeClassifier.cxx:190-198`）→ `tool/shape_classifier.rs`（`CorrectFace2d` `:908-924` 见 T-261），详见 `_tasks.md` T-260
- [x] `TopOpeBRepBuild_CorrectFace2d` + `TopOpeBRepBuild_Tools2d`（T-261）：`Perform`（`CorrectFace2d.cxx:587-642`）、`CheckFace`（`:226-336`）、`BuildCopyData`（`:116-222`）、`ConnectWire`（`:340-585`）、`MakeRightWire`（`:646-683`）、`GetP2dFL`（`:687-702`）、`MakeHeadList`（`:706-744`）、`CheckList`（`:748-847`）、`TranslateCurve2d`（`:851-862`）、`OuterWire`（`:866-890`）、`BndBoxWire`（`:894-907`）、`MoveWire2d`（`:911-977`）、`MoveWires2d`（`:982-1030`）、两个 `UpdateEdge` 重载（`:1033-1117`）→ `build/correct_face_2d/{mod,copy_data,wires_2d}.rs`；`Tools2d::Path`（`Tools2d.cxx:38-66`）、`BuildPath`（`:70-158`）、`MakeMapOfShapeVertexInfo`（`:162-210`）、`DumpMapOfShapeVertexInfo`（`:213-243`）→ `build/tools2d.rs`；`TopOpeBRepBuild_Tools::CorrectFace2d` wrapper（`Tools.cxx:908-924`）→ `build/tools/pcurve.rs` 的 `correct_face_2d`（`CORRECT_FACE2D_UNPORTED` 占位已替换），详见 `_tasks.md` T-261
- [x] `TopOpeBRepBuild_Builder1` 族（T-262）：`GLOBAL_USE_NEW_BUILDER`（`Builder1.cxx:61`）、构造（`:65-73`）、`Clear`（`:86-97`）、`Perform` 两个重载（`:101-104`、`:108-141`）、`MergeKPart` 两个重载（`:145-148`、`:152-182`）→ `build/builder/builder1/mod.rs`；`PerformShapeWithStates` 两个重载 + `PerformFacesWithStates`（`Builder1_1.cxx:183-552`）→ `builder1/shape_states.rs`；`StatusEdgesToSplit`（`:556-855`）、`SplitEdge`（`:859-1115`）、`ClassifyEdgeToSolidByOnePoint`（`:1117-1144`）→ `builder1/status_edges.rs`；`CorrectResult2d`（`Builder1_2.cxx:34-224`）→ `builder1/correct_result_2d.rs`；`HBuilder` 成员改为 `Builder1`（`HBuilder.hxx:170`）（`GFill*` 覆盖族与 `PerformONParts`/`PerformPieceIn2D` 由 T-269 补齐），详见 `_tasks.md` T-262
- [x] `TopOpeBRepBuild_fctwes` + `GFABUMakeFaces` + `FUN_tool_ClosedW` + `TopOpeBRepTool_2d`（T-263）：`GFillCurveTopologyWES` 两个重载（`fctwes.cxx:60-136`、`:138-268`）→ `build/builder/curve_topology_wes.rs`；`FUN_purgeFon1nonoriE`（`makefaces.cxx:248-276`）、`FUN_ReOrientIntExtEdge`（`:284-333`）、`FUN_CheckORI`（`:337-360`）、`GFABUMakeFaces`（`:365-755`）→ `build/builder/make_faces.rs`；`FUN_getv`（`PURGE.cxx:75-101`）、`FUN_addtomap`（`:58-73`）、`FUN_tool_ClosedW`（`:103-296`）→ `tool/purge.rs`；`FC2D_*` 全族（`2d.cxx:96-406`）+ `TopOpeBRepTool_C2DF` → `tool/c2d.rs`（`GWESMakeFaces` 与 `CORRISO`/`PurgeClosingEdges`/`MakeFaces`/`CorrectONUVISO`/`RegularizeFaces` 仍待，见 T-264），详见 `_tasks.md` T-263
- [x] `TopOpeBRepTool` 的 CORRISO / REGUW / RegularizeW / CLASSI / connexity / face 与 `TopOpeBRepBuild` 的 FREGU / GWESMakeFaces（T-264）：`TopOpeBRepTool_face.cxx:26-129` → `tool/face.rs`；`TopOpeBRepTool_CLASSI.cxx:47-615` → `tool/classi.rs`；`TopOpeBRepTool_connexity.cxx:47-243` → `tool/connexity.rs`；`TopOpeBRepTool_REGUW.cxx:119-1325` → `tool/reguw/{mod,map_s,block,regu}.rs`；`TopOpeBRepTool_RegularizeW.cxx:86-887` → `tool/regularize_w.rs`；`TopOpeBRepTool_TOOL::WireToFace`（`TOOL.cxx:1943-1974`）→ `tool/wire_to_faces.rs`；`TopOpeBRepBuild_FREGU.cxx:39-372` → `build/builder/fregu.rs`；`TopOpeBRepBuild_makefaces.cxx:85-246` → `build/builder/gwes_make_faces.rs`（`GridFF` 的 `GMergeFaces`/`GSplitEdgeWES`/`GMergeEdgeWES`/`GSplitEdge` 与 `GridSS` 全族仍待，见 T-264(f)），详见 `_tasks.md` T-264
- [x] `TopOpeBRepTool_REGUS` 与 matter 辅助（T-265）：`FUN_tool_isobounds`（`TOPOLOGY.cxx:155-220`）、`FUN_nearestISO`（`:585-643`）、`FUN_tool_nC2dINSIDES`（`matter.cxx:45-62`）、`FUN_tool_getdxx`（`:150-230`）、`FUN_tool_getxx`（`:106-147`）、`tryOriEinF`（`TOOL.cxx:779-814`）、`NgApp`（`:818-918`）、`tryNgApp`（`:922-944`）、`uvApp`（`:1491-1512`）、`outUVbounds`（`:1915-1921`）、`MkShell`（`:1819-1827`）、`Remove`（`:1831-1848`）→ `tool/{uv_app,ng_app,mk_shell}.rs`；`TopOpeBRepTool_REGUS.cxx:88-1039` 全类 → `tool/regus/{mod,map_s,split_faces,regu,block,nearest_f}.rs`；`TopOpeBRepTool_RegularizeS.cxx:39-74` 的 `RegularizeShells` → `tool/regularize_s.rs`；`TopOpeBRepBuild_ShellToSolid.cxx:25-61` 与 `TopOpeBRepBuild_SREGU.cxx:41-333` 见 T-266，详见 `_tasks.md` T-265
- [x] `TopOpeBRepBuild_ShellToSolid` + `Builder::RegularizeSolids` / `RegularizeSolid`（T-266）：`ShellToSolid.cxx:25-61`（`Init`/`AddShell`/`MakeSolids`）→ `build/shell_to_solid.rs`；`SREGU.cxx:41-126`（`RegularizeSolids`：`myMemoSplit` 门 + 同域实体 split 表回填）、`:129-333`（`RegularizeSolid`：`RegularizeShells` → `ShellToSolid::MakeSolids` → `myFSplits` 二级切分写回 `ChangeSplit`）→ `build/builder/sregu.rs`（`myFSplits` 从此有写入点，`myMemoSplit` 首次使用），详见 `_tasks.md` T-266
- [x] `TopOpeBRepBuild_GridFF` 全部 + `ffwesk` + `keep` + `GToSplit`（T-267）：`keep.cxx:41-202` → `build/builder/keep.rs`；`Grid.cxx:30-91` 的 `GToSplit` → `build/builder/grid.rs::g_to_split`；`ffwesk.cxx:47-258` → `build/builder/ffwesk.rs`；`GridFF.cxx:144-327`（`FUN_computeLIFfaces2d` 两个重载）、`:335-463`（`GMergeFaces`）、`:474-507`（`FUN_validF1edge`）、`:511-635`（`GFillFacesWESMakeFaces`）、`:639-706`（`GFillFaceWES`）、`:709-757`（`GFillWireWES`）→ `build/builder/grid_ff.rs`；`:761-800`（`GFillEdgeWES`）、`:803-814`（`FUN_samgeomori`）、`:826-1207`（`GSplitEdgeWES`）、`:1217-1347`（`GMergeEdgeWES`）、`:1350-1485`（`GSplitEdge`）、`:1493-1521`（`ClassifyEdgeToSolidByOnePoint`）→ `build/builder/grid_ff_edges.rs`；`makeedges.cxx:41-76` 的 `GPVSMakeEdges` → 同上；`EXPORT.cxx:3279-3306` 的 `FUN_ds_hasFEI` → `ds/export.rs`；`TOPOLOGY.cxx:1548-1564` 的 `FUN_tool_findAncestor` → `tool/topology.rs`（`GLOBAL_classifysplitedge` 的写入方 `ffsfs.cxx` 未译，恒 `false`），详见 `_tasks.md` T-267
- [x] `TopOpeBRepBuild_GridSS` 的 SFS 族 + `ffsfs` + `makesolids`（T-268）：`makesolids.cxx:35-172` → `build/builder/make_solids.rs`；`GridSS.cxx:282-350`（`GMergeSolids`）、`:353-428`（`GFillSolidsSFS`）、`:433-494`（`GFillSolidSFS`）、`:500-590`（`GFillSurfaceTopologySFS` 两个重载）、`:595-645`（`GFillShellSFS`）→ `build/builder/grid_ss.rs`；`:832-1028`（`GSplitFaceSFS`）、`:1031-1182`（`GMergeFaceSFS`）→ `build/builder/grid_ss_faces.rs`；`:1348-1453`（`GSplitFace`）→ `build/builder/grid_ss_split.rs`；`:1454-1592`（`AddONPatchesSFS`）→ `build/builder/grid_ss_on.rs`；`:648-803`（`FUNBUILD_MAPSUBSHAPES`/`FUNBUILD_MAPANCSPLSHAPES`/`FUNBUILD_ANCESTORRANKPREPARE`/`FUNBUILD_ANCESTORRANKGET`/`FUNBUILD_ORIENTLOFS` 与 `GLOBAL_*` 全局）→ `build/builder/grid_ss_globals.rs`；`ffsfs.cxx:55-77`（`FUN_setmotherope` 族）、`:111-282`（`FUN_getAncestorFsp` 三负载）、`:293-691`（`GFillFaceSFS`）→ `build/builder/ffsfs.rs`（`FillOnPatches`/`AreFacesCoincideInArea`/`FindFacesTouchingEdge`/`SplitEvisoONperiodicF` 由 T-269 补齐），详见 `_tasks.md` T-268
- [x] `TopOpeBRepBuild_Builder1` 的 `GFill*` 虚函数族 + `GWESMakeFaces` 覆写 + `SplitEvisoONperiodicF`（T-269）：`Builder1.cxx:186-532`（`GFillSolidSFS`/`GFillShellSFS` 两个虚覆写 + `GFillFaceNotSameDomSFS`/`GFillFaceNotSameDomWES`/`GFillWireNotSameDomWES`/`GFillEdgeNotSameDomWES`）→ `builder1/gfill_not_same_dom.rs`；`:533-1150`（`GFillFaceSameDomSFS`/`GFillFaceSameDomWES`/`GFillWireSameDomWES`/`GFillEdgeSameDomWES`）→ `builder1/gfill_same_dom.rs`；`:1373-1412`（`GWESMakeFaces` 覆写：`KP_ISSOSO` 走 `NormalizeFace`、随后全体 `CorrectFace2d`）→ `builder1/gwes_make_faces.rs`；子类字段 `hxx:192-203` → `builder1/state.rs` 的 `Builder1State`（挂 `Builder.b1`，配 `Builder.is_builder1` 做虚分派）；`GridSS.cxx:89-116`（`FUN_getG`）、`:138-190`（`FUN_EPIforEvisoONperiodicF`）、`:237-279`（`FUN_getEPI`/`FUN_getEPIonEds`）、`:1184-1294`（`FUN_SplitEvisoONperiodicF`）、`:1314-1346`（`SplitEvisoONperiodicF`）→ `build/builder/split_eviso_on_periodic.rs`，`Perform` 的 `cxx:125` 接上；`GridSS.cxx:1441` 的 `FillOnPatches` 调用点接回 `on_patches.rs`；`Merge.cxx` 的 `Opec12`/`Opec21`/`Opecom`/`Opefus`/`Classify`/`ChangeClassify` 与 `TopExp::MapShapes` 迁到 `build/builder/merge_states.rs`（`builder/mod.rs` 755 行），详见 `_tasks.md` T-269
- [x] `TopOpeBRepBuild_FuseFace` 全类 + `GFillFaceSFS` 的 fuse 尾段（T-270）：`TopOpeBRepBuild_FuseFace.cxx:49-1398` 全类 → `build/fuse_face/{mod,perform_face,perform_edge,group_shape,group_edge,make_edge,same_support,shape_list_map}.rs`（`Init` `:73-118`、`PerformFace` `:126-582`、`PerformEdge` `:590-681`、`ClearEdge` `:688-822`、`ClearVertex` `:829-841`、`GroupShape` `:845-971`、`GroupEdge` `:975-1074`、`MakeEdge` `:1078-1177`、`SameSupport` `:1181-1398`）；`ffsfs.cxx:618-667` 的 `Init(ldum, *GLOBAL_lfr1, addinternal=1)` + `PerformFace` + `IsDone` + `LFuseFace` → `SFS.AddStartElement` 接回 `build/builder/ffsfs.rs`（`MakeEdge` 的 `BRepLib_MakeEdge(curve, V1, V2)` 重建由 T-272 落地）
- [x] `BRepLib_MakeEdge::Init(C,V1,V2)` / `Project` + `BRepTools_Substitution` + `TopOpeBRepTool_FuseEdges` + `Builder::End` 第三块（T-272）：`BRepLib_MakeEdge.cxx:55-109`（`Project`）、`:543-572`（`Init(C,V1,V2)`）、`:603-796`（`Init(C,V1,V2,p1,p2)` 顶点臂）→ `brep_lib_make_edge.rs`（`adaptor_curve` 复刻 `GeomAdaptor_Curve` 的「剥 trim、区间留基参数」语义）；`BRepTools_Substitution.cxx:35-171` 全类 → `brep_tools_substitution.rs`；`TopOpeBRepTool_FuseEdges.cxx` 全 1,075 行 → `topopebrep/tool/fuse_edges/{mod,connex,result_edges,update_pcurve,perform,same_support}.rs`；`End.cxx:366-410`（`FuseEdges` 块）+ `:413-520`（`UpdateSplitAndMerged`）→ `build/builder/end_fuse_edges.rs`，`Builder::end()` 接第三块；`FuseFace.cxx:1078-1177` 的 `MakeEdge` 重建臂接 `init_3d_projecting`，详见 `_tasks.md` T-272
- [x] `ProjLib` 的 projector 族（T-273）：`ProjLib_Projector.cxx:34-251`（`IsDone`/`Done`/`GetType`/`SetType`/`IsPeriodic`/`SetPeriodic`/`Line`/`Circle`/`Ellipse`/`Hyperbola`/`Parabola`/五个 `Project` 基体 `:145-176`/`UFrame` `:179-199`/`VFrame` `:204-224`/`Set(Get)Bezier`/`Set(Get)BSpline`）→ `crates/occt-geom/src/projlib/projector.rs`；`ProjLib_Plane.cxx:24-169`（`Init` `:76-83`、`Project(gp_Lin)` `:102-107`、`Project(gp_Circ)` `:111-122`、`Project(gp_Elips)` `:126-139`、`Project(gp_Parab)` `:143-155`、`Project(gp_Hypr)` `:159-169`）→ `projlib/plane.rs`；`ProjLib_Cylinder.cxx:28-175`（`Project(gp_Lin)` `:96-118`、`Project(gp_Circ)` `:122-157`、`Project(gp_Elips)` `:159-164` 空体、`:167-176` 基体）→ `projlib/cylinder.rs`；`ProjLib_Cone.cxx:30-180`（`Project(gp_Lin)` `:69-104`、`Project(gp_Circ)` `:108-165`、`:167-180` 基体）→ `projlib/cone.rs`；`ProjLib_Sphere.cxx:33-249`（`EvalPnt2d` `:66-92`、`Project(gp_Circ)` `:98-179`、`SetInBounds` `:206-249`、`:182-200` 基体）→ `projlib/sphere.rs`；`ProjLib_Torus.cxx:28-193`（`EvalPnt2d` `:59-77`、`Project(gp_Circ)` `:82-174`、`:176-193` 基体）→ `projlib/torus.rs`；`projlib.rs` 改为 `projlib/mod.rs` + 六个类文件并再导出 `Projector`/`Plane`/`Cylinder`/`Cone`/`Sphere`/`Torus`；`gp_Dir::IsNormal(other, tol)`（`gp_Dir.hxx:171-174`）补进 `occt-core/src/gp/dir.rs`，详见 `_tasks.md` T-273
- [x] `ProjLib_ProjectedCurve` + `TopOpeBRepTool_CurveTool` 的 pcurve 构造（T-274）：`ProjLib_ProjectedCurve.cxx:25-1178` 全类（`ComputeTolU` `:53-61`、`ComputeTolV` `:64-74`、`IsoIsDeg` `:80-126`、`TrimC3d` `:129-159`、`Project` `:241-270` 双分发、四个构造 `:274-327`、`Perform` `:369-886` 的解析臂与球面极点裁剪 `:466-493`、周期收尾 `:787-886`、`Set*`/`GetTolerance`/`GetType`）→ `crates/occt-geom/src/projlib/projected_curve.rs`（`projected_curve.rs:1-25` 逐条记录三个 UNPORTED 近似器）；`CurveTool.cxx:111-139`（`MakePCurve`）、`:975-1037`（`IsProjectable`）、`:1040-1156`（`MakePCurveOnFace` 的 trim/`MakePCurve`/球面准周期镜像 `:1082-1116`/U·V 周期平移 `:1118-1153`）→ `crates/occt-topo/src/topopebrep/tool/curve_tool.rs`（`is_projectable`/`make_pcurve_on_face`/`make_pcurve`/`degree1_knots`）；配套 `Geom2dBSplineCurve::knot(index)`（1 基去重结）、`GeomTrimmedCurveBasis` 的 `gp_line`/`gp_ellipse`/`gp_hyperbola`/`gp_parabola` 转发，详见 `_tasks.md` T-274
- [x] `TopOpeBRepDS_BuildTool` 的曲线重算族 + `ShapeTool::Resolution3d`（T-275）：`BuildTool.cxx` 的 `GetGeomTool` `:134-138` / `ChangeGeomTool` `:142-147` / `Approximation` `:1315-1318`（`TypeC3D() != BSPLINE1`）与 `geom_tool` 成员（`hxx:210`，默认 `GeomTool(APPROX)`）→ `ds/build_tool/mod.rs`；`FUN_getUV` `:557-573`（`GeomAPI_ProjectPointOnSurf` + `LowerDistanceParameters`）、`FUN_makeUisoLineOnSphe` `:602-651`（采样 `0.243234`/`0.543345`、`vdir = ±Y`、`SetLin2d` 重建）、`ComputePCurves` `:654-789`、`PutPCurves` `:793-813`、`RecomputeCurves` `:816-890` → `ds/build_tool/recompute_curves.rs`；`TopOpeBRepTool_ShapeTool::Resolution3dU/V/Resolution3d(SU)/(F)`（`ShapeTool.cxx:655-694`）→ `topopebrep/tool/shape_tool.rs`；`MakePCurveOnFace` 补 `uv_bounds` 参数，把「球面镜像读自然区间（`BAS` 以 `Restriction = false` 构造）」与「`du`/`dv` 读 `BRepTools::UVBounds`」两个来源分开（`CurveTool.cxx:1055`、`:1073-1074`、`:1085-1106`），详见 `_tasks.md` T-275
- [x] `TopOpeBRepBuild_Builder::BuildEdges` 两个重载（T-276）：`BuildEdges.cxx:38-100`（`BuildEdges(iC, HDS)`：`MakeEdge(anEdge, curC, DS)` → `PaveSet`/`PaveClassifier`/`EdgeBuilder` → `MakeEdges` → 逐边 `RecomputeCurves` → `RemoveCurve(iC)` 或 `UpdateEdge`）、`:104-218`（`BuildEdges(HDS)`：`ick` 定位 → keep 曲线 `RemoveCurve` + `ChangeNbCurves(ick-1)` → 对 `Mother == 0` 的曲线走单曲线重载 → `tp` 点数引用 → `RemovePoint`）→ `crates/occt-topo/src/topopebrep/build/builder/build_edges.rs`，`perform.rs` 在 `Builder.cxx:126` 调用 `build_edges_all()`；`FindKeep` 默认值按 `CurveExplorer.hxx:37` / `DataStructure.hxx:254`、`:149` 取 `true`，两处前序遍历显式 `false`，详见 `_tasks.md` T-276
- [ ] `TopOpeBRepBuild` 余项（36/43 类）：`BuilderON`（`BuilderON.cxx:2903`，仅 `BuilderON2d` 的 `GLOBAL_faces2d` 调试分支未译，见 T-257/T-271）、`GridSS`（`:1903`，SFS 族与 `FillOnPatches`/`AreFacesCoincideInArea`/`FindFacesTouchingEdge`/`SplitEvisoONperiodicF` 已由 T-268/T-269 落地）、`Builder`（`:1813`，`Perform` 的三个步骤 `BuildEdges`/`BuildFaces`/`SplitEvisoONperiodicF` 与 `Perform`/`Merge*`/`Section`/`BuildVertices`/`SplitSectionEdge*`/`KPart*`/`End` 块 1-3（块 3 见 T-272）与全部 `G*` 已落地，见 T-252、T-254、T-255、T-256、T-258、T-259、T-269、T-272、T-276）、`Section`（`:1263`）与 `Tools` / `makefaces`（`GWESMakeFaces` 已由 T-264 落地） / `Merge` / `ffsfs`（`GFillFaceSFS` 主体与 `FuseFace` 尾段已由 T-268/T-270 落地） / `Griddump` / `FaceBuilder` / `WireEdgeSet` / `WireEdgeClassifier` / `PaveSet` / `PaveClassifier` / `AreaBuilder` 族 / `BlockBuilder` / `ShellFaceSet` / `ShellFaceClassifier` / `SolidBuilder` / `CompositeClassifier` / `LoopSet` / `PBuilder` / `PGTopo` / `GTopo` / `GTool` / `GIter` / `VertexInfo` / `ShapeSet` / `ShapeListOfShape` / `PWireEdgeSet` / `EdgeBuilder` / `FaceAreaBuilder` / `SolidAreaBuilder` / `Loop` / `LoopClassifier` / `LoopEnum` / `Pave` / `kpresu`
- [ ] `IntPatch` 余项（25/39 类已落地，缺 4,983 行）：`RstInt`（`RstInt.cxx:1266`）、`InterferencePolyhedron`（`:1091`）、`Polyhedron`（`:631`）、`TheIWLineOfTheIWalking`（`:257`）、`PolyArc`（`:256`）、`PolyLine`（`:247`）、`PolyhedronBVH`（`:147`）、`CSFunction`（`:125`）、`BVHTraversal`（`:93`）、`HCurve2dTool`（`:52`）、`CurvIntSurf`（`:90`）、`Polygo`（`:82`）、`PolyhedronTool`（`:69`）、`SpecPntType` / `SearchPnt`
- [x] `IntPatch_LineConstructor` 尾段（T-237）：`TestWLineIsARLine`（`LineConstructor.cxx:949-1005`）、`TestIfWLineIsRestriction`（`:1008-1071`）、`ProjectOnArc`（`:1073-1094`）、`TestWLineToRLine`（`:1095-1432`）→ `intpatch_line_ctor/wline_to_rline.rs`；`Perform`（`:1433-2085`：Analytic `:1451-1514`、Walking `:1515-1689`、Geometric `:1690-1867`、Restriction `:1868-2084`）→ `intpatch_line_ctor/perform.rs`；`IntPatch_RLine::Add(LineOn2S)` 的 `curv` 点序列 → `geom_int_types.rs` 的 `RLine::points`；`IntPatch_Point::ReverseTransition`（`IntPatch_Point.cxx:87-137`）→ `PatchPoint::reverse_transition`
- [x] `TopOpeBRep_FacesIntersector`（T-238，`FacesIntersector.hxx:35-140` / `.cxx:155-1752`：`Perform` `:166-266`、`PrepareLines` `:376-485`、`TestWLineAlongRestriction` `:650-702`、`BuildRLineBasedOnWLine` `:704-868`、`BuildRLine` `:874-1058`、`TestWLinesToAnArc` `:1064-1210`、`MergeWLinesIfAllSegmentsAlongRestriction` `:1226-1344`、`GetArc` `:1357-1623`、`IsPointOK` `:1631-1652`、`GetPointOn2S` `:1658-1685`、`GetMergedWLineOnRestriction` `:1692-1752`）→ `topopebrep/faces_intersector/{mod,rline_build,merged_wline,wlines_to_arc}.rs`
- [x] `TopOpeBRep_EdgesIntersector`（T-239，`EdgesIntersector.hxx:30-206` / `.cxx:1-1199` / `_1.cxx:1-464`：`SetFaces` `.cxx:98-176`、`ForceTolerances` `:177-185`、`TransitionEqualAndExtremity` `:186-217`、`IsTangentSegment` `:219-251`、`EdgesIntersector_checkT1D` `:253-308`、`Perform` `:316-689`、`Dimension` `:691-706`、`ComputeSameDomain` `:708-771`、`SetSameDomain` `:772-779`、`MakePoints2d` `:780-821`、`ReduceSegment` `:822-917`、`ReduceSegments` `:918-950`、点/边访问器 `:951-1164`；`_1.cxx:32-464` 的 `Segment1` / `IsOpposite1` / `InitPoint1` / `MorePoint1` / `NextPoint1` / `Point1` / `Status1` / `IsPointOfSegment1` / `Index1` / `EdgesConfig1` / `Transition1` / `Parameter1` / `IsVertex1` / `Vertex1` / `Value1`）→ `topopebrep/edges_intersector/{mod,perform,point_ops}.rs`
- [x] `TopOpeBRep_ShapeIntersector`（T-239，`ShapeIntersector.hxx:57-180` / `.cxx:72-1341`：`InitIntersection` `:373-408`、`NextIntersection` `:550-622`、六个偶对族 `Init/Find/More/Next` `:624-1032`、`CurrentGeomShape` `:119-181`、`RejectedFaces` `:187-366`、`OneShapeIsHalfSpace` `:1085-1236`、`GetNewSolid` `:1244-1341`）→ `topopebrep/shape_intersector/{mod,couples,rejected_faces}.rs`
- [ ] 按 `TopOpeBRepDS.cxx` → `TopOpeBRep.cxx` → `TopOpeBRepBuild.cxx` → `TopOpeBRepTool.cxx` 的依赖顺序
- [ ] 门禁：编译 + `occt-topo --lib` + `bopalgo_builder_gates` / `bop_builder2_boss`

**S3（TKFillet，85,017 行）** — 前置 S1、S2 已收口，本波已开工

- [x] `ChFiDS` 值层（T-279）：`ChFiDS_State` / `TypeOfConcavity` / `ChamfMethod` / `ChamfMode` / `ErrorStatus` → `occt-topo/src/chfi_ds/enums.rs`；`ChFiDS_CircSection` + `NCollection_HArray1<ChFiDS_CircSection>`、`ChFiDS_CommonPoint`、`ChFiDS_FaceInterference`、`ChFiDS_Regul`、`ChFiDS_Map`、`ChFiDS_SurfData` → `chfi_ds/{circ_section,common_point,face_interference,regul,map,surf_data}.rs`
- [x] `ChFiDS` 引导线与条带（T-280）：`ChFiDS_ElSpine`、`ChFiDS_Spine`（+ `SpineHandle`）、`ChFiDS_Stripe`（+ `StripeHandle`）、`ChFiDS_StripeMap` → `chfi_ds/{el_spine,spine/,stripe,stripe_map}.rs`；配套 `GCPnts_AbscissaPoint` 带初值反演（`occt-geom/src/gcpnts/abscissa.rs`）与 `CurveAdaptor` 的 `Curve` 实现（`topopebrep/tool/adaptor.rs`）
- [x] `ChFiDS` 两类专用引导线（T-281）：`ChFiDS_FilSpine`（半径律，`parandrad` + 每条基本脊线的 `Law_Composite`）→ `chfi_ds/fil_spine/{mod,law}.rs`；`ChFiDS_ChamfSpine`（截面距离）→ `chfi_ds/chamf_spine.rs`；子类状态按 `SpineClass` 折进 `Spine`
- [x] `ChFi3d` 自由函数层（T-282）：`DefineConnectType` / `IsTangentFaces` / `ConcaveSide` / `NextSide`（两个重载）/ `SameSide` / `ChFi3d_Coefficient` / `ChFi3d_FilletShape` / `ChFi3d_SearchSing` → `chfi3d/{connect_type,tangent_faces,side,coefficient,fillet_shape,search_sing}.rs`；前置依赖 `LocalAnalysis_SurfaceContinuity` + `LocalAnalysis_StatusErrorType` → `local_analysis/surface_continuity.rs`，`BRepTools::OriEdgeInFace` → `brep_tools.rs`，`GeomLProp_SLProps::D1U/D1V/D2U/D2V` → `occt-geom/src/geom_lprop/sl_props.rs`，`Adaptor3d_Surface::D2` → `topopebrep/tool/adaptor.rs`
- [x] `ChFi3d_Builder_0.cxx` 辅助自由函数（T-283）：`ChFi3d_InPeriod` / `ChFi3d_SearchPivot` / 两个 `ChFi3d_Boite` / `ChFi3d_conexfaces` / `ChFi3d_EdgeState` / `ChFi3d_evalconti` / `ChFi3d_cherche_{face1,element,edge,vertex}` / `ChFi3d_nbface` / `ChFi3d_edge_common_faces` / `ChFi3d_AngleEdge` / `ChFi3d_ChercheBordsLibres` / `ChFi3d_Nb{,NotDegenerated}Sharp?Edges` / `ChFi3d_Couture{,OnVertex}` / `ChFi3d_IsPseudoSeam` / `ChFi3d_CircularSpine` / `ChFi3d_Spine` / `ChFi3d_Parameters` / `ChFi3d_IsSmooth` → `chfi3d/builder_helpers/{period,box2d,faces,vertex_edges,topology_stats,seam,spine_geom,surface_params,smooth}.rs`；前置依赖 `BRepTools::IsReallyClosed` → `brep_tools.rs`，`GeomLProp_CLProps::SetParameter` / `CentreOfCurvature` → `occt-geom/src/geom_lprop/cl_props.rs`，`IntAna_QuadQuadGeo` 平面/平面线臂 `plane_plane_line` 提升为 `pub` → `occt-geom/src/intana/analytic_intersections.rs`。仍 UNPORTED：`ChFi3d_BoundSrfImpl`/`BoundFac`/`BoundSrf`（`:158-240`、`:691-715`）、`ChFi3d_KParticular`（`:604-683`）、`ChFi3d_InterPlaneEdge`（`:718-767`）、`ChFi3d_ApproxByC2`（`:5835-5851`）
- [x] `ChFi3d_Builder` 数据层（T-284）：构造器 + 五个 `ChFiDS_Map` 采集 + `SetParams` / `SetContinuity` + 结果/故障访问器（`IsDone` / `Shape` / `NbFaultyContours` / `FaultyContour` / `NbComputedSurfaces` / `ComputedSurface` / `NbFaultyVertices` / `FaultyVertex` / `HasResult` / `BadShape` / `StripeStatus`）+ 轮廓访问器（`Remove` / `Value` / `NbElements` / 两个 `Contains` / `Length` / `FirstVertex` / `LastVertex` / `Abscissa` / `RelativeAbscissa` / `Closed` / `ClosedAndTangent`）+ `FaceTangency` + `Reset` + `Generated` + `UpdateTolesp` → `chfi3d/builder/mod.rs`；前置依赖 `ChFiDS_Spine::GetTolesp`（`ChFiDS_Spine.lxx:167-170`）→ `chfi_ds/spine/mod.rs`。仍 UNPORTED（等 `BRepBlend` / `BlendFunc`）：`Compute` / `ExtentAnalyse` / `PerformSingularCorner` / `PerformFilletOnVertex` / `PerformExtremity` / `PerformElement` / `SimulKPart` / 三个 `SimulSurf` 重载 / `PerformTwoCornerbyInter` / `SplitKPart` / `ComputeData` / `CompleteData` / `SimulData` / `StoreData` / `PerformSetOfSurfOnElSpine` / `ConexFaces` / `PerformOneCorner` / `IntersectMoreCorner` / `PerformMoreSurfdata` / `PerformIntersectionAtEnd` / `FindFace` / `MoreSurfdata` / `StartSol`（其中 `ComputeData` / `CompleteData` / `SimulData` / `StoreData` 已在 T-294 落地）
- [x] `TKFillet` 混合值层（T-285）：`Blend_Point`（9 个构造 + 10 个 `SetValue` 重载 + 20 个访问器，`Blend_Point.hxx:31-308` / `.lxx:19-200` / `.cxx:21-483`）→ `occt-topo/src/blend/point.rs`；`BRepBlend_PointOnRst`（`hxx:34-75` / `.lxx:16-35` / `.cxx:22-47`）→ `brep_blend/point_on_rst.rs`；`BRepBlend_Extremity`（`hxx:42-155` / `.lxx:16-87` / `.cxx:26-142`）→ `brep_blend/extremity.rs`；`BRepBlend_Line`（`hxx:37-186` / `.lxx:16-130` / `.cxx:25-58`）→ `brep_blend/line.rs`。仍 UNPORTED：`Blend_Function` / `Blend_FuncInv` / `Blend_SurfRstFunction` / `Blend_SurfPointFuncInv` / `Blend_SurfCurvFuncInv` / `Blend_RstRstFunction` / `Blend_CurvPointFuncInv` / `Blend_CurvPointRadInv` / `Blend_RstRstEval`（`blend/mod.rs` 头注释）；`SurfCurv*Inv` / `SurfPoint*Inv`（`brep_blend/mod.rs` 头注释；`BRepBlend_AppFunc*` / `AppSurf` / `AppSurface` 已在 T-293 落地，`Blend_Function` 族 / `BlendTool` / `HCurve*Tool` 已在 T-286 落地，`Walking` / `CSWalking` / `SurfRstLineBuilder` / `RstRstLineBuilder` 已在 T-288..T-291 落地，四个 `*Rad` 曲面函数已在 T-292 落地）
- [x] `Blend` 函数接口族 + `BRepBlend` 几何查询（T-286）：`Blend_AppFunction`（`hxx:39-192` / `.cxx:20-23`）→ `blend/app_function.rs`；`Blend_Function`（`hxx:44-121` / `.cxx:18-58`）→ `blend/function.rs`；`Blend_CSFunction`（`hxx:44-199` / `.cxx:18-52`）→ `blend/cs_function.rs`；`Blend_SurfRstFunction`（`hxx:38-223` / `.cxx:19-35`）→ `blend/surf_rst_function.rs`；`Blend_RstRstFunction`（`hxx:42-233` / `.cxx:19-35`）→ `blend/rst_rst_function.rs`；`Blend_FuncInv`（`hxx:37-94` / `.cxx:18-21`）→ `blend/func_inv.rs`；`Blend_SurfPointFuncInv` / `Blend_SurfCurvFuncInv` / `Blend_CurvPointFuncInv`（`hxx:31-88` / `.cxx:21-25`）→ `blend/inv_functions.rs`；`Blend_Status` / `Blend_DecrochStatus`（`Blend_Status.hxx:19-32` / `Blend_DecrochStatus.hxx:19-28`）→ `blend/status.rs`；`BRepBlend_BlendTool`（`hxx:35-90` / `.cxx:17-142` / `.lxx:17-58`）→ `brep_blend/blend_tool.rs`；`BRepBlend_HCurveTool` / `HCurve2dTool` 的 `NbSamples`（两个 `.cxx:25-57`）→ `brep_blend/curve_tool.rs`；前置 `Extrema_EPCOfExtPC2d` = `Extrema_GGenExtPC`（`Extrema_GGenExtPC.hxx:36-232`，`Perform` `:160-173`）→ `occt-geom2d/src/extrema2d/epc_of_ext_pc2d.rs`。仍 UNPORTED：`Blend_CurvPointRadInv` / `Blend_RstRstEval` / `Blend_AppSurface`、`BRepBlend_AppFunc*` / `AppSurf` / `AppSurface`（等 `AppBlend_AppSurf`）、四个行走器与七个 `BRepBlend_*Inv` / 四个 `BRepBlend_*Rad`（后三项已在 T-288..T-292 落地；`BlendFunc` 段形包已在 T-292 落地，`ComputeNormal` / `ComputeDNormal` 除外；`BRepBlend_AppFunc*` / `AppSurf` / `AppSurface` 已在 T-293 落地）
- [x] `BRepBlend` 逆函数五类（T-287）：`BRepBlend_SurfPointConstRadInv`（`hxx:32-84` / `.cxx:17-276`）→ `brep_blend/inv/surf_point_const_rad_inv.rs`；`BRepBlend_SurfPointEvolRadInv`（`hxx:32-84` / `.cxx:17-269`）→ `brep_blend/inv/surf_point_evol_rad_inv.rs`；`BRepBlend_SurfCurvConstRadInv`（`hxx:33-90` / `.cxx:17-284`）→ `brep_blend/inv/surf_curv_const_rad_inv.rs`；`BRepBlend_SurfCurvEvolRadInv`（`hxx:33-90` / `.cxx:17-303`）→ `brep_blend/inv/surf_curv_evol_rad_inv.rs`；`BRepBlend_CurvPointRadInv`（`hxx:31-80` / `.cxx:17-130`）→ `brep_blend/inv/curv_point_rad_inv.rs`。四个行走器 `BRepBlend_Walking` / `CSWalking` / `SurfRstLineBuilder` / `RstRstLineBuilder` 已在 T-288..T-291 全部落地，四个 `*Rad` 曲面函数已在 T-292 落地，`BRepBlend_AppFunc*` / `AppSurf` / `AppSurface` 已在 T-293 落地
- [x] `BRepBlend_Walking` 行走器（T-288）：`BRepBlend_Walking`（`hxx:41-186` / `.cxx:95-2769`）→ `brep_blend/walking/`。构造 + `SetDomainsToRecadre` + `AddSingularPoint` + `Perform` 与 `ClassificationOnS1/S2` / `Check` / `Check2d` / `IsDone` / `Line` / `TwistOnS1/S2` → `walking/mod.rs`；`PerformFirstSection` / `PerformFirstSectionRecadre` → `walking/first_section.rs`；`TestArret` 的两种重载 → `walking/test_arret.rs`；`Recadre` / `ArcToRecadre` → `walking/recadre.rs`；`Transition` / `MakeExtremity` / `MakeSingularExtremity` / `CorrectExtremityOnOneRst` / `RecadreIfPeriodic` → `walking/extremity.rs`；`Continu` 两个重载 + `Complete` → `walking/continu.rs`；`InternalPerform` 与 `OnBothRestrictions` / `SetBothExtremities` / `RecordPoint` / `EvalPInit` → `walking/internal_perform.rs`；前置 `IntSurf::MakeTransition` → `intsurf_transition.rs`，`Adaptor3d_TopolTool::Identical` → `brep_top_adaptor_topol_tool.rs`，`gce_MakePln(P1,P2,P3)` → `occt-geom/src/intana/line_torus.rs` 的 `plane_through_three`
- [x] `BRepBlend_CSWalking` 行走器（T-289）：`BRepBlend_CSWalking`（`hxx:46-117` / `.cxx:88-1118`）→ `brep_blend/cs_walking/`。构造 + `Perform` + `Complete` + `IsDone` / `Line` → `cs_walking/mod.rs`；`TestArret` / `CheckDeflectionOnSurf` / `CheckDeflectionOnCurv` → `cs_walking/test_arret.rs`；`Transition` / `MakeExtremity` → `cs_walking/extremity.rs`；`InternalPerform` → `cs_walking/internal_perform.rs`；复用 `brep_blend::walking::{solve_root_set, arc_curve}`。仍 UNPORTED：`BRepBlend_CSWalking::Recadre`（`cxx:589-738`，OCCT 整段位于 `/* */` 内）
- [x] `BRepBlend_SurfRstLineBuilder` 行走器（T-290）：`BRepBlend_SurfRstLineBuilder`（`hxx:64-176` / `.lxx:19-50` / `.cxx:137-1777`）→ `brep_blend/surf_rst_line_builder/`。构造 + `ArcToRecadre` + `Perform` + `IsDone` / `Line` / `DecrochStart` / `DecrochEnd` → `surf_rst_line_builder/mod.rs`；`PerformFirstSection` + `Complete` → `surf_rst_line_builder/first_section.rs`；`InternalPerform` → `surf_rst_line_builder/internal_perform.rs`；三个 `Recadre` 重载 → `surf_rst_line_builder/recadre.rs`；`Transition` + `MakeExtremity` → `surf_rst_line_builder/extremity.rs`；`CheckDeflectionOnSurf` / `CheckDeflectionOnRst` / `TestArret` / `CheckInside` / `ConvOrToTra` → `surf_rst_line_builder/test_arret.rs`；复用 `brep_blend::walking::{solve_root_set, arc_curve}`。四个 `*Rad` 曲面函数已在 T-292 落地，`BRepBlend_AppFunc*` / `AppSurf` / `AppSurface` 已在 T-293 落地
- [x] `BRepBlend_RstRstLineBuilder` 行走器（T-291）：`BRepBlend_RstRstLineBuilder`（`hxx:64-205` / `.lxx:19-61` / `.cxx:151-1964`）→ `brep_blend/rst_rst_line_builder/`。构造 + `Perform` + `IsDone` / `Line` / `Decroch1Start` / `Decroch1End` / `Decroch2Start` / `Decroch2End` → `rst_rst_line_builder/mod.rs`；`PerformFirstSection` + `Complete` → `rst_rst_line_builder/first_section.rs`；`InternalPerform` → `rst_rst_line_builder/internal_perform.rs`；四个 `Recadre` 重载 → `rst_rst_line_builder/recadre.rs`；`Transition` + `MakeExtremity` → `rst_rst_line_builder/extremity.rs`；`CheckDeflectionOnRst1` / `CheckDeflectionOnRst2` / `TestArret` / `CheckInside` / `ConvOrToTra` → `rst_rst_line_builder/test_arret.rs`；复用 `brep_blend::walking::{solve_root_set, arc_curve}`。四个 `BRepBlend` 行走器至此全部落地。四个 `*Rad` 曲面函数已在 T-292 落地，`BRepBlend_AppFunc*` / `AppSurf` / `AppSurface` 已在 T-293 落地
- [x] `BRepBlend` 四个 `*Rad` 曲面函数（T-292）：`BRepBlend_SurfRstConstRad`（`hxx:41-244` / `.cxx:51-1045`）→ `brep_blend/surf_rst_const_rad/`：`Value` / `Derivatives` / `Values` + `t3dto2d` + `CenterCircle` 族 → `equations.rs`；`IsSolution` + `Decroch` → `solution.rs`；三个 `Section` + `Resolution` + 查询族 → `section.rs`；构造 / `Set` / 三个 trait impl → `mod.rs`。`BRepBlend_SurfRstEvolRad`（`hxx:41-248` / `.cxx:113-1166`）→ `surf_rst_evol_rad/`：三方程 + `t3dto2d` → `equations.rs`；`IsSolution` + `Decroch` → `solution.rs`；`FusionneIntervalles` + 查询族 → `queries.rs`；`Section` 族 → `section.rs`；`mod.rs` 带 `tevol` / `fevol` 与 `Set(Choix)` 的 `sg1`。`BRepBlend_RstRstConstRad`（`hxx:41-265` / `.cxx:53-981`）→ `rst_rst_const_rad/`。`BRepBlend_RstRstEvolRad`（`hxx:41-268` / `.cxx:117-1092`）→ `rst_rst_evol_rad/`：两方程 + `CenterCircleRst1Rst2` → `equations.rs`；`IsSolution` + `Decroch` → `solution.rs`；`FusionneIntervalles` → `queries.rs`；`Section` 族 → `section.rs`。四个类都实现在 `Blend_RstRstFunction` / `Blend_SurfRstFunction` + `Blend_AppFunction` + `math_FunctionSetWithDerivatives` 上，全部经 `brep_blend/mod.rs` 的 `pub use` 导出；四个类共用的 `BlendFunc` 段形包（`GetShape` / `GetMinimalWeights` / `NextShape` + `BlendFunc_SectionShape`，`BlendFunc.hxx:31-71` / `.cxx:40-152`）→ `blend/blend_func.rs`。仍 UNPORTED：`BRepBlend_AppFunc*` / `AppSurf` / `AppSurface`；`BlendFunc::ComputeNormal` / `ComputeDNormal`（`.cxx:155-279`，依赖未移植的 `CSLib`）
- [x] `BRepBlend` 近似层（T-293）：`BRepBlend_AppFuncRoot`（`hxx:41-185` / `.cxx:26-432`）→ `brep_blend/app_func/`：构造 + `D0` / `D1` / `D2` + `Nb2dCurves` / `SectionShape` / `Knots` / `Mults` / `IsRational` / `NbIntervals` / `Intervals` / `SetInterval` / `Resolution` / `GetTolerance` / `SetTolerance` / `BarycentreOfSurf` / `MaximalSection` / `GetMinimalWeight` → `app_func/mod.rs`；`SearchPoint`（`:300-376`）+ `SearchLocation`（`:387-432`）+ `FuncSetRef` 适配 → `app_func/search.rs`；`Point` / `Vec` 三个子类的实现（`BRepBlend_AppFunc.cxx:33-44`、`AppFuncRst.cxx:32-44`、`AppFuncRstRst.cxx:33-44`）→ `app_func/contact.rs`。`BRepBlend_AppSurface`（`hxx:33-122` / `.cxx:25-137` / `.lxx:26-94`）→ `brep_blend/app_surface.rs`。`BRepBlend_AppSurf`（`hxx:44-193` / `BRepBlend_AppSurf_0.cxx:17-30`，即 `AppBlend_AppSurf.gxx` 以 `TheSectionGenerator = Blend_AppFunction`、`TheLine = BRepBlend_Line` 实例化）→ `brep_blend/app_surf.rs`：复用 `occt_geom::geom_fill::AppSurf` 与 `AppBlendSectionGenerator`，`LineSectionGenerator` 把 `Line` 的点序喂给 `Blend_AppFunction`；`perform` / `perform_smoothing` / `perform_nb_max_p` 三个入口。配套改动：`occt-core/src/approx/sweep_function.rs` 去掉 `ApproxSweepFunction` 的 `Send + Sync` 约束（OCCT 是普通继承 + `void*` 别名，port 用 `Rc<RefCell<..>>`）；`occt-geom/src/approx_sweep_approximation.rs` 补 `Dump`（`cxx:846-872`）。消费方 `ChFi3d_Builder_6.cxx:528-610`（`CompleteData` 的三个重载）已在 T-294 落地，`AppSurf` 实例在 OCCT 8.0.0 里无调用方。仍 UNPORTED：`AppBlend_AppSurf::PerformSmoothing` 的 `AppDef_Variational` 臂（`AppBlend_AppSurf.gxx:472-541`）；`ChFi3d_Builder::CompleteData` 本体已在 T-294 落地。
- [x] `ChFi3d_Builder_6.cxx` 算法层（T-294）：`CompParam`（`_6.cxx:209-302`）→ `chfi3d/builder_helpers/comp_param.rs`；`UpdateLine`（`:344-387`）→ `chfi3d/builder/update_line.rs`；四个 `CompleteData`（`:393-624`：Geom/Iso、`Blend_Function`+`AppFunc`+`AppSurface`、`Blend_SurfRstFunction`、`Blend_RstRstFunction`）→ `chfi3d/builder/complete_data.rs`；`StoreData`（`:629-1005`）→ `chfi3d/builder/store_data.rs`；两个 `ComputeData`（SurfRst `:1009-1171` / RstRst `:1175-1342`）→ `chfi3d/builder/compute_data.rs`；两个 `SimulData`（SurfRst `:1346-1501` / RstRst `:1506-1669`）→ `chfi3d/builder/simul_data.rs`；`ChFi3d_CheckSameParameter`（`ChFi3d_Builder_0.cxx:1565-1595`）+ 两个 `ChFi3d_SameParameter`（`:1602-1640`）→ `chfi3d/builder_helpers/same_parameter.rs`；`ChFi3d_TrsfTrans`（`:2174-2186`）+ `ChFi3d_FilCommonPoint`（`:2193-2296`）→ `chfi3d/builder_helpers/fil_common_point.rs`。前置：`BSplCLib::TangExtendToConstraint`（`BSplCLib.cxx:3794-4238`）+ `PLib::CoefficientsPoles`（`PLib.cxx:1522-1603`）→ `occt-geom/src/bspl_tang_extend.rs`；`BSplCLib::Eval` 的任意维重载（`:3640-3790`）→ `occt-core/src/bspl/eval_extrap.rs` 的 `eval_dim`；`GeomLib::ExtendKPart` / `ExtendSurfByLength`（`GeomLib.cxx:1416-1975`）→ `occt-geom/src/geom_lib_extend_surf.rs`；`Geom_RectangularTrimmedSurface::IsSurfaceOfLinearExtrusion` 转发 basis（`GeomAdaptor_Surface.cxx:423-425` 的同一解包）→ `occt-geom/src/rectangular_trimmed.rs`；`Geom2dInt_GInter` 的 `ComputeDomain` 显式区间重载（`Geom2dInt_GInter.gxx:120-177`）→ `occt-geom2d/src/geom2d_int/ginter.rs` 的 `domain_over_range`；`BRepBlend_PointOnRst::arc` 改存 `Curve2dAdaptor`（`ChFi3d_FilCommonPoint` 按 `_0.cxx:2238-2239` 下转 `BRepAdaptor_Curve2d` 读 `Edge()`），四个行走器的 `add_arc` 同步改签名；`Geom_BSplineSurface` 的 `Rational(Weights)` 判定（`Geom_BSplineSurface.cxx:110-138`）→ `occt-geom/src/bspline_surface/mod.rs` 的 `is_really_rational`（`StoreData` 据此决定是否保留权重，替换掉原来按 `PConfusion` 偏离 1 的 port-only 判据）。仍 UNPORTED：`ChFi3d_Builder_6.cxx` 的面/面重载 `ComputeData(..., Spine, ..., Blend_Function&, Blend_FuncInv&, ...)`（`:1674-2508`）与 `SimulData(..., AdditionalHGuide, ...)`（`:2509-2703`）**已在 T-295 落地**；`ChFi3d_Builder::Compute` / `ExtentAnalyse` / `Perform*` 族与 `ChFi3d_ChBuilder` / `ChFi3d_FilBuilder`
- [x] `ChFi3d_Builder` 面/面行走层（T-295）：`ChFi3d_Builder::ComputeData` 面/面重载（`ChFi3d_Builder_6.cxx:1676-2507`）→ `chfi3d/builder/compute_data_face_face/`（`mod.rs` 驱动 + `extension.rs` 的 `apply_extensions`，`:2336-2505`）；`SimulData` 面/面重载（`:2511-2703`）→ `chfi3d/builder/simul_data_face_face.rs`；`ChFi3d_Builder::SearchFace`（`ChFi3d_Builder_2.cxx:1523-1627`）→ `chfi3d/builder/search_face.rs`；其三个文件静态 `IsInput`（`_2.cxx:399-471`）/ `IsG1`（`:475-521`）/ `SearchFaceOnV`（`:532-589`）→ `chfi3d/builder_2_helpers/search_face.rs`；`ChFi3d_Builder_6.cxx` 的文件静态 `SearchIndex`（`:94-103`）/ `IsVois`+`IsObst`+`nbedconnex`（`:107-205`）/ `CompBlendPoint`（`:304-339`）→ `chfi3d/builder_6_helpers/`；前置 `ChFiDS_ElSpine::SetOrigin`（`ChFiDS_ElSpine.cxx:230-244`）→ `chfi_ds/el_spine.rs` 的 `ElSpine::set_origin`（经 `Curve::as_bspline_curve` 下转 + `GeomBSplineCurve::set_origin_u`，`Geom_BSplineCurve.cxx:913-970`）。`ChFi3d_Builder_6.cxx` 的四个 `ComputeData` / 四个 `SimulData` 至此全部落地。`PerformExtremity` / `PerformElement` 已在 T-296 落地；仍 UNPORTED：`ChFi3d_Builder::Compute` / `ExtentAnalyse` / `PerformTwoCorner*` / `PerformThreeCorner` / `PerformFourCorner` / `SplitKPart` 族，以及 `ChFi3d_ChBuilder` / `ChFi3d_FilBuilder`
- [x] `ChFi3d_Builder` 脊线行走层（T-296）：`ChFi3d_Builder_1.cxx` 的五个文件静态 `ReorderFaces`（`:65-123`）/ `ConcatCurves`（`:125-147`）/ `MakeOffsetEdge`（`:148-304`）/ `TangentExtremity`（`:625-711`）/ `TangentOnVertex`（`:713-728`）→ `chfi3d/builder/perform_helpers.rs`；`ChFi3d_Builder::PerformExtremity`（`:730-885`）→ `chfi3d/builder/perform_extremity.rs`，其 `ConexFaces`（`ChFi3d_Builder_2.cxx:859-880`）→ 同文件的 `ChFi3dBuilder::cone_x_faces`；`ChFi3d_Builder::PerformElement`（`:887-1179`）→ `chfi3d/builder/perform_element.rs`（下游 + 上游两段切向行走、`SetEdges` / `PutInFirst` 与两条 offset 链）。仍 UNPORTED：`ChFi3d_Builder::Compute` / `ExtentAnalyse` / `PerformSingularCorner` / `PerformFilletOnVertex` / `PerformTwoCorner*` / `PerformThreeCorner` / `PerformFourCorner` / `SplitKPart` 族，以及 `ChFi3d_ChBuilder` / `ChFi3d_FilBuilder`
- [x] `ChFi3d_Builder` 曲面驱动层（T-297）：`ChFi3d_Builder::CallPerformSurf`（`ChFi3d_Builder_2.cxx:635-818`）→ `chfi3d/builder/call_perform_surf.rs`；`ChFi3d_Builder::StripeOrientations`（`:826-857`）→ 同文件；OCCT 的两个 protected 纯虚面/面虚函数 `PerformSurf`（`ChFi3d_Builder.hxx:445-466`）/ `SimulSurf`（`:203-221`）→ `chfi3d/builder/surf_ops.rs` 的 `ChFi3dSurfOps`。六组边接触 `PerformSurf` / `SimulSurf` 重载已在 T-299 以基类 `throw` 的默认方法落地，`PerformSetOfKPart` 已在 T-301 落地；仍 UNPORTED：`PerformSetOfKGen`（`:3298-`）、`PerformSetOfSurf`（`:3882-`）、`SimulKPart` 的实装（`hxx:201`；基类默认体已进 `ChFi3dSurfOps`）、`Compute` 家族与两个子类构造器
- [x] `ChFi3d_Builder` 首解搜索层（T-298）：`ChFi3d_Builder::StartSol`（`ChFi3d_Builder_2.cxx:900-1108`）→ `chfi3d/builder/start_sol.rs`；其两个文件静态 `BonVoisin`（`:213-318`）/ `Projection`（`:320-370`）→ `chfi3d/builder_2_helpers/start_sol.rs`；`ChFi3d_InterPlaneEdge`（`ChFi3d_Builder_0.cxx:718-767`）→ `chfi3d/builder_0_helpers/inter_plane_edge.rs`；`PerformFirstSection`（`ChFi3d_Builder.hxx:577-588`）→ `ChFi3dSurfOps`。第二个 `StartSol` 重载与 `ChFi3d_BuildPlane` 已在 T-299 落地，`PerformSetOfSurfOnElSpine` 已在 T-300 落地，`PerformSetOfKPart` 已在 T-301 落地；仍 UNPORTED：`SimulKPart` 的实装（`hxx:201`）、`PerformSetOfKGen`（`:3298-`）、`PerformSetOfSurf`（`:3882-`）
- [x] `ChFi3d_Builder` 障碍侧首解层（T-299）：第二个 `StartSol` 重载（`ChFi3d_Builder_2.cxx:1162-1520`，`BRepAdaptor_Curve2d` 障碍点版本）→ `chfi3d/builder/start_sol_obstacle.rs`；`ChFi3d_BuildPlane`（`:1113-1144`）→ `chfi3d/builder_2_helpers/build_plane.rs`；`isTangentToArc`（`:192-209`）→ `chfi3d/builder_2_helpers/is_tangent_to_arc.rs`；`BRep_Tool::Parameters(V, F)`（`BRep_Tool.cxx:1661-1703`）→ `brep_tool.rs` 的 `parameters_of_vertex_on_face`；六组边接触 `PerformSurf` / `SimulSurf`（`hxx:223-302`、`:470-556`，基类体 `ChFi3d_Builder_NotImp.cxx:29-220` 全部 `throw`）→ `ChFi3dSurfOps` 的六个默认方法；`ChFi3d_Purge`（`_2.cxx:2047-2101`）/ `InsertAfter`（`:2108-2132`）/ `RemoveSD`（`:2136-2172`）/ `InsertBefore`（`:2179-2203`）→ `chfi3d/builder/stripe_surgery.rs`。`PerformSetOfSurfOnElSpine` 已在 T-300 落地，`PerformSetOfKPart` 已在 T-301 落地；仍 UNPORTED：`PerformSetOfKGen`（`_2.cxx:3298-`）、`PerformSetOfSurf`（`:3882-`）、`SimulKPart` 的实装（`hxx:201`）、`Compute` 家族与两个子类构造器
- [x] `ChFi3d_Builder` 脊线行走主体（T-300）：`PerformSetOfSurfOnElSpine`（`ChFi3d_Builder_2.cxx:2207-3003`）→ `chfi3d/builder/perform_set_of_surf/`：`mod.rs`（驱动 + `Walk` 状态 + `OffsetHGuide` 下标回写 `flush_offset`，`:2207-2230` / `:2233-2257` / `:2989-2999`）、`start.rs`（`:2283-2412`：`StartSol(Stripe, ...)`、周期段收缩、双向重启、`Vref`、`nbed` 重写、`bidf`-`bidl`、`MaxStep`）、`step.rs`（`:2416-2988`：延长端 `intf`/`intl`、`ChFi3d_BuildPlane`、`CallPerformSurf`、`done` 两臂、`Inside` 收尾、周期同弧 `return`、反向重启与终止判定）、`obstacles.rs`（`:2567-2809`：both / S1 / S2 三支的 `Choix` 与 `Soldep`）、`splice.rs`（`:2870-2928`：`InsertAfter` / `InsertBefore` + `ChFi3d_Purge`）。前置：`BRepAdaptor_Surface()` 空构造（`BRepAdaptor_Surface.cxx:28`）→ `topopebrep/tool/adaptor.rs` 的 `SurfaceAdaptor::new()`；`ChFiDS_Spine::ChangeOffsetElSpines()` 读口（`ChFiDS_Spine.hxx:151-154`）→ `chfi_ds/spine/mod.rs`。三障碍支调用的边接触 `PerformSurf` / `SimulSurf` 仍是 `ChFi3d_Builder_NotImp.cxx:29-220` 的 `throw` 默认体，因此在当前基线恒 `Err`，待 `ChFi3d_FilBuilder` / `ChFi3d_ChBuilder` 落地后才可执行。仍 UNPORTED：`PerformSetOfKGen`（`:3298-`）、`PerformSetOfSurf`（`:3882-`）、`SimulKPart` 的实装（`hxx:201`）、`Compute` 家族与两个子类构造器
- [x] `ChFi3d_Builder` 特例角行走层（T-301）：`ChFi3d_Builder::PerformSetOfKPart`（`ChFi3d_Builder_2.cxx:3004-3281`）→ `chfi3d/builder/perform_set_of_k_part.rs`（边扫描 + `SurfData` 序列 + 周期乱序重排 `:3082-3121` + `ElSpine` 配对 `:3170-3213` + 收尾区间 `:3222-3265` + 两个 `PerformElSpine` 循环 `:3267-3278`）；前置：`ChFi3d_KParticular`（`ChFi3d_Builder_0.cxx:604-683`）→ `chfi3d/builder_0_helpers/k_particular.rs`；`ChFi3d_SolidIndex`（`_0.cxx:2300-2321`）→ `chfi3d/builder_0_helpers/solid_index.rs`；文件静态 `TgtKP`（`_2.cxx:373-398`）→ `chfi3d/builder_2_helpers/tgt_kp.rs`；`ChFi3d_BoxDiag`（`_2.cxx:3282-3297`）→ `chfi3d/builder_2_helpers/box_diag.rs`；`ChFi3dSurfOps::SimulKPart`（`hxx:201`，基类默认体 `ChFi3d_Builder_NotImp.cxx` 的 `throw`）→ `chfi3d/builder/surf_ops.rs` 默认方法。三个惰性挂点：`ChFiKPart_ComputeData::Compute`（`ChFiKPart` 包）→ `chfi3d/chfi_k_part.rs`、`ChFi3d_Builder::SplitKPart`（`ChFi3d_Builder_SpKP.cxx:749-1160`，前置 `Geom2dHatch` + `Tri`/`FillSD`/`CpSD`）→ `chfi3d/builder/split_k_part.rs`、`ChFi3d_PerformElSpine`（`_0.cxx:4797-5352`，前置 `ChFi3d_ApproxByC2` + `CurveCleaner`）→ `chfi3d/builder_0_helpers/perform_el_spine.rs`，三者均答 `Ok(false)` / `Err` 并写清 `.cxx` 行号。`ChFi3d_Builder_2.cxx` 的三条行走至此全部落地（T-300 的 `PerformSetOfSurfOnElSpine`、T-301 的 `PerformSetOfKPart`）；仍 UNPORTED：`PerformSetOfKGen`（`:3298-3881`）、`PerformSetOfSurf`（`:3882-3900`）、上面三个惰性挂点的主体、`ChFi3d_Builder::Compute` / `ExtentAnalyse` / `PerformTwoCorner*` / `PerformThreeCorner` / `PerformFourCorner`、`ChFi3d_ChBuilder` / `ChFi3d_FilBuilder`
- [ ] `ChFi3d` 构造器家族（38,891 行）：`ChFi3d_Builder`（`Builder.cxx` 916 / `Builder_0.cxx` 5,655 / `Builder_1.cxx` 1,257 / `Builder_2.cxx` 3,722 / `Builder_6.cxx` 2,575 / `Builder_C1.cxx` 5,071 / `Builder_C2.cxx` 778 / `Builder_CnCrn.cxx` 3,773 / `Builder_SpKP.cxx` 1,451）→ `ChFi3d_ChBuilder`（3,190）→ `ChFi3d_FilBuilder`（5,173）→ `BRepBlend` 曲面族（19,287 行）→ `BlendFunc` 段形（14,652 行）→ 变半径律 → `ChFi2d` / `ChFiKPart` / `ChFiDS` / `Blend` / `FilletSurf` / `BRepFilletAPI`
- [ ] 门禁：编译 + `occt-topo --lib` + `step_obj_gates` + `export_data_obj`

**S4 / S5 / S6 / S7** — 见 §5.3 落点表；S5/S6/S7 可与 S1/S2 并行

---

## 10. 收尾与文档同步

每波收尾时同步更新：

- `specs/_tasks.md` —— 任务台账（新增 T 编号、记录收口证据与门禁读数）
- `specs/_audit/_index.md` —— 对应 A 条目的「处置」行与状态
- `docs/goals/drift-and-gap-autofix/state.yaml` —— `goal.oracle` 读数与包收据
- `CLAUDE.md` / `README.md` —— 若模块落点、命名或门禁命令有变

**引用规范**：报告里引用源码用**真实文件名 + 物理行号**（`step/read_geometry.rs:476`），
不要引用已废弃的 `pNN`，不要用 `read` 工具的折行显示行号。

---

### 附：本文的对照基线

| 基线 | 位置 |
|---|---|
| OCCT 参考网格 | `data/occ-*.obj`（15 个）、`data/occ/occ-*.obj`（8 个） |
| 本仓产物 | `data/output/*.obj`（`export_data_obj`）、`data/output/gates/*.obj`（`step_obj_gates`） |
| OCCT 侧探针 | `specs/occt_probe/`（C++，需本机 OCCT 构建） |
| 审计条目 | `specs/_audit/_index.md`（A0–A33） |
| 门禁基线冻结记录 | `docs/goals/drift-and-gap-autofix/notes/T002-gate-baseline.md` |
