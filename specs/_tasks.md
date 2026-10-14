# _tasks — OCCT → Rust migration map (current)

> Last updated: 2026-10-14. Implementation `.rs` files are split to <=800 lines.
> Phase 3 detail: [specs/_phase3_plan.md](_phase3_plan.md)

## 已完整 (✅)

### occt-core (208 tests)
| 域 | 模块 |
|---|---|
| gp 几何 (37 类型) | gp/ 全部 — 点/线/圆/椭圆/双曲/抛物/平面/柱/锥/球/环/2D/变换 |
| 求值 | elib/ — ElCLib + ElSLib + intersect + surface_eval + measure |
| B样条/Bezier | bspl/ — knots/eval/bezier/plib/surface/poles/rational/曲线工具 (rational eval 修复) |
| 包围体 | bnd/ — Box/Box2d/Sphere/OBB(PCA)/SortBox |
| 空间结构 | bvh/ — BVH + 遍历 + 三角BVH |
| 网格 | poly/ — Triangulation/Polygon2D/Polygon3D |
| 交换 | io/ — OBJ + STL + PLY 读写 |
| 建模基础 | csg/ (体素布尔) · hull/ · geom/ (采样/三角剖分/拟合/折线/多边形布尔) |
| 质量属性 | gprop/ — 质心/惯性/面积/体积 |
| 内核 | kernel/ — OCCError/Handle/容器/TCollection/TColStd/GeomAbs/Units/OSD/Resource/字符串 |
| 其他 | numeric/ · quantity/ · message/ · toploc/ · convert/ · cslib/ · gcpnts/ · precision |

### occt-math (160 tests)
vector/matrix/intvec/status + SVD/Crout/Jacobi/Householder/Newton/BFGS/Powell/
Gauss-Legendre/Kronrod/Eigen/Trig/LeastSquares/MultiInt + Cholesky/LU/QR-full/
cubic-spline/polyfit/Levenberg/RK4/bisection/secant/golden/brent/polynomial
roots + matrix_ext + stats/fft/rng + **P6**: spline_surface (张量积 B样条曲面插值) ·
eig_qr (Hessenberg+QR 一般矩阵特征) · nurbs_fit (最小二乘 B样条拟合) ·
distrib_extra (χ²/t/F/Γ/Beta/Poisson 分布族)

### occt-geom (TKG3d, 52 tests)
Curve/Surface traits + Line/Circle/Ellipse/Hyperbola/Parabola/Plane/Cylinder/Cone/
Sphere/Torus + B样条/Bezier/Trimmed/Offset/Revolved + **B1-B5 转换** +
**P6 extrema**: 点-曲线/曲线-曲线/曲线-曲面/曲面-曲面最近距离
(Extrema_ExtPC/CC/PS/CS/SS)

### occt-geom2d (TKG2d, 38 tests)
Curve2d trait + Line/Circle/Ellipse/Hyperbola/Parabola/Trimmed/Offset/BSpline +
Geom2dBezierCurve · curve_ops(长度/交点/最近点) · bspline2d_to_bezier

### occt-topo (TKBRep, `--lib` 1276 tests)
- **几何内核**: tgeometry(GeometryRegistry 侧表) · 真实 children 树 ·
  builder(全部 make_*) · brep_tool(BRep_Tool) · primitives(真几何 box/cyl/sphere/cone/torus)
- **网格化**: wireframe · shape_mesh(任意 BRep 网格化+法线+焊接) · incremental_mesh(逐行移植 BRepMesh_FaceDiscret+Delaun)
- **拓扑分析**: face_face · edge_split · shell_check · solid_union ·
  boolean_ops · brep_extrema · shape_analysis · brep_gprop
- **精确建模 (P5+P6)**: inttools + bop_builder (平面精确布尔) ·
  **intpatch (曲面解析求交) + bop_curved (曲面实体布尔)** ·
  **fillet_edge (滚动球边/角倒圆)** · **brep_offset (多边形/面/壳偏移)**
- **交换**: STEP(ISO 10303-21) · IGES(ANSI Y14.26M) · OBJ/STL/PLY/VRML ·
  **bincaf (二进制装配) · gltf (glTF 2.0)**
- **其他**: sweep · loft · pipe · sweep_revolve · fillet · hlr · render_svg ·
  feature(凸台/凹槽) · patterns · shape_metrics · brep_projection · xcaf

## 里程碑状态

| 里程碑 | 状态 |
|---|---|
| M1 (A1-A4): BRep 可存几何、可查询 | ✅ make_box 真几何, BRepTool 取回坐标/曲面 |
| M2 (A5-A7 + B): 任意形状可网格化/可视化/导出 | ✅ 任意 BRep 三角化 → OBJ/STL/PLY/IGES/glTF |
| M3 (C+D): 网格导入 + 基础布尔 | ✅ 体素布尔 + mesh→BRep + 交换 |
| M4 (E, P5): 平面精确布尔 + Deflection 网格 + SVG | ✅ 精确平面 Fuse/Cut/Common |
| M5 (F, P6): 曲面求交/布尔 + 边倒圆 + 偏移 | ✅ 球/柱面精确布尔 + 滚动球倒圆 |

## BOPAlgo remaining vs OCCT (large)

Do not grow the dual-pipeline / planar / volume-gate scaffold. Do not invent FillIn3D IN counts or spoke filters. Do not change `IsHole`. Small wiring (BuildSplitFaces edge images, `SetContext`, history `myImages` keys) is deferred.

| Id | Gap | OCCT | Port | Notes |
|---|---|---|---|---|
| B1 | Dual boolean entry + `bop_builder.rs` size | `BOPAlgo_BOP` only | `boolean()` unchanged; entry ~200 lines; API/tests/planar split under 800 | Gate not shrunk. Split 2026-08-29 |
| B2 | `BuilderSolid::PerformAreas` BoxTree / `myBoxes` | `BOPAlgo_BuilderSolid.cxx:397` | Hole `AabbTree` + `IsInside` owner + `SetWhole` unbound; `IsHole` still containment | Ported 2026-08-29 |
| B3 | `BuilderSolid::PerformInternalShapes` | `ClassifyFaces` at `:673` | `classify_faces_occt`; uses `myBoxes` when present | Ported 2026-08-29 |
| B4 | `PostTreat` avoid-map | `BOPAlgo_Builder.cxx:450` `aMA` | `correct_tolerances_avoid`; filled when `NonDestructive` | Ported 2026-08-29 |
| B5 | CompSolid rebuild | `BOPAlgo_BOP.cxx:1267` `CollectContainers` | Collect + connexity wrap; SOLID args still early-return | Ported 2026-08-29 |
| B6 | Leftover split-solids | `BuildSplitSolids` + `aMST` | Live: `bop_split_solids_occt` + `bop_bop::build_solid` | Isolated in private `bop_build_solids` leftover submodule 2026-08-29 |
| B7 | `IntTools_Context::ComputeVE` (edge/edge) | `IntTools_EdgeEdge` | Sampled min-distance; no live PaveFiller caller | Deferred — VE is vertex/edge `pave_ve`; EE is `edge_edge.rs` |
| B8 | PaveFiller BoxTree | `IsExistingPaveBlock` / `IntersectVertices` | `AabbTree` + `bnd2bvh3d`; same `IsOut` | Ported 2026-08-29 |

## 远期 (⬜)
- BOPAlgo 完整（自交修复/多结果拓扑）/ TopOpeBRep 全拓扑
- TKFillet 曲面-曲面一般混合、多法向复杂角
- TKOpenGl 硬件渲染 / TKService 字体纹理全量 / TKV3d 交互选择高亮 / TKIVtk
- Draw 完整 TCL（表达式解析/过程/命令表）
- BRepGProp 更全（体积惯性、张量高级）
- XCAF 完整 schema（视图/标注/实例树）/ STEP 全实体写

## 已知自创残余（须按 OCCT 控制流收口，勿用启发式替代）

| Id | 残余 | 现状 |
|---|---|---|
| T-125 | `brepfeat::solid_volume` 用线框端点扇形和 + Newell 法向（`brepfeat/features.rs`，本轮前 `:506/553/567/583`），OCCT 无此规则；同族自创件为网格夹具冒充解析件（`features.rs:377 translated_mesh_cylinder`、`:399 mesh_cylinder(..,24)`、`:325 revolve_profile_about` 车削网格、`:465 resolution_for` 的 `clamp(16,64)`） | **2026-10-13 分开实测，结论一分为二**：<br>(1) **体积积分已换忠实件**：`solid_volume` 删除自创扇形和/Newell 法向，改调 `brep_gprop_full::volume_properties`（`BRepGProp::VolumeProperties` → `BRepGProp_Vinert`/`BRepGProp_Gauss::Vinert`，`BRepGProp.cxx:555-586`、`BRepGProp_Vinert.cxx:197-204`），签名去掉已死的 `_nu/_nv` UV 网格参数，删掉随之无调用者的 `polygon_fan_volume`/`polygon_normal`/`wire_polygon_3d`。`brepfeat --lib` **17/17 绿**，`occt-topo --lib` **1276/1276 绿**；同一 shape 上对拍，gprop 比自创扇形和更准：rib 4.090000（精确）vs 4.150000、neck 4.489835（= 16 边形精确）vs 4.816391、boss 4.155291（= 24 边形精确）vs 4.284701。**旧记载的"gprop 对结果面积分有缺口"不成立**——`volume_properties` 在 `brepfeat` 结果面上是正确的。<br>(2) **解析件替换仍被 BOP 阻塞（未做）**：3 项失配全部来自 `bop_builder::boolean`（`BOPAlgo_BOP`）对解析旋转体/贯通圆柱工具的处理，与体积无关。逐项隔离实测（`brepfeat` 探针，2026-10-13）：box ∪ 解析 neck（z=1 贴合）→`res_solid=true` 但工具丢失、体积 = box 4.000000（应 4.489835）；box ∪ 解析贯通 boss → `res_solid=false`、整根圆柱未裁剪地加上去 2.706858（应 2.565487）；96 面圆柱 − 解析 groove 环（`Cut`）→ `res_solid=false` 且刃口**增加**体积 9.337869 > 9.317486。换成车削网格夹具三者全部回到基线 ⇒ 缺口属 **A5/T-79（`bop_curved`/`BOPAlgo_BOP` 曲面分支）**，不属网格器（`GeomCylinder::d1` 返回真实偏导，`occt-geom/src/cylinder.rs:20-22`）。网格夹具因此保留并就地标为 UNPORTED（`features.rs` `revolve_profile_about`/`mesh_cylinder`/`brepfeat/mod.rs` 模块头），工作理由由原来那句错的"解析面 `d1` 为零"改为上述实测数字 |
| T-79 | `bop_curved` 网格级布尔（曲面 operand，`bop_curved/region_trim.rs:51`、`face_meshing.rs:418`） | **2026-10-14 复核：分派已收口**——`boolean_dispatch` 的曲面分支改调 `bop_builder::boolean_via_bopalgo`（`bop_builder_dispatch.rs:227-232`），`curved_boolean_full`/`_ext`/`_volume` 已无生产调用点（仅测试），T-80 已闭环、旧阻塞解除。**仍未修**：忠实路径自身失配（`zz_bop_curved_probe`：`FUSE` faces=11/vol=10.72000 vs OCCT 8.50265/faces=8），根因是盒面 × 柱面求交曲线的容差 `3.734e0`（应 1e-7 量级），致盒顶面未切开、柱侧面塌成 0.048 碎面。见 `_audit/bop-boolean.md` 发现 3 处置、`_audit/_index.md` A5 |

## A26 IGES 写侧余项（2026-10-13 收口）

| Id | 余项 | 状态 |
|---|---|---|
| T-78a | 周期面反周期化：`GeomToIGES_GeomSurface.cxx:303-343` 的 `SetUOrigin`/`SetVOrigin` 重定原点（issue 26138）+ `SetUNotPeriodic` | ✅ 已落地（`iges/emit_surface.rs` 的 `u_reorigin`/`v_reorigin` 块（`:316-355`）→ `GeomBSplineSurface::set_u_origin`/`set_v_origin` → `bspl::unperiodize_direction`）；旧"未移植"注释已清除 |
| T-78b | 全周椭圆支：`GeomToIGES_GeomCurve.cxx:620-645`（104→126）；圈整周支 `cxx:565-570` | ✅ 已落地（`iges/emit_conic.rs::emit_whole_period_ellipse`（`:17`）、`iges/emit_curve.rs::emit_circle_transfer`（`:126`）的 `plast = pfirst`） |
| T-78c | 对拍基线 | `iges_check` 14/14 `ok`；参数段记录多重集 vs `_igesdump/oracle/*.iges`：10/14 逐字节一致（Shape-2 的 28 条 128 记录逐字节一致），3 个模型浮点噪声/读侧极点密度差异，`OffsetPlaneHoleEdge` 仅缺 314/406（颜色/名字） |
| T-78d | 仍 UNPORTED（有 OCCT 出处，非自创） | Bezier 窄区间 `Geom_BezierCurve::Segment`；`cxx:476-502` 的 BRep-mode `analyticMode` 臂（默认 Faces 模式不可达） |

## 文件拆分记录（>800 行按 OCCT 职责切分，控制流不变）

命名规则见 `CLAUDE.md`「文件与模块命名」：`<module>/mod.rs` + `mod <主题>;` + `pub use <主题>::*;`，文件名取该文件承载的 OCCT 类/职责。拆分为纯搬迁：`pub use` 重导出保持原调用路径，**未改任何分支**。

| 日期 | 原文件 | 行数 | 新布局（各行数） |
|---|---|---|---|
| 2026-10-13 | `occt-topo/src/iges.rs` | 4361 | `iges/{mod 90, entry 48, writer 356, card 312, entity 188, curve_kind 403, frame_geometry 130, wire_helpers 239, pcurve_transfer 412, emit_conic 345, emit_curve 424, emit_surface 618, emit_face 604, tests 77}.rs`；`iges::write_iges` 调用路径不变 |
| 2026-10-13 | `occt-topo/src/pcurve_full/surface_projector.rs` | 3223 | `pcurve_full/surface_projector/{mod 42, classify 194, projector 170, analytic 78, newton_uv 485, closed_state 403, iso_param 643, period_jump 276, curve_projection 633, reparam 271}.rs`；`super::surface_projector::*` 重导出不变 |
| 2026-10-09 | `occt-topo/src/shhealing/intersection_tool.rs` | 1465（intake 快照数，仓库无该单体版本的提交） | `shhealing/intersection_tool.rs` 66（模块根）+ `shhealing/intersection_tool/{helpers 80, edge_split 302, vertex_union 184, self_intersect 759, fix_intersecting_wires 901}.rs`（**非空行数**；总行数依次 92/330/194/786/934）；`ShapeFix_IntersectionTool` 调用路径不变 |

**引用规则**：规格/报告/注释里旧写的 `iges.rs:NNNN`、`surface_projector.rs:NNNN` 已按新文件名就近订正（如 `iges_curve_kind` → `iges/curve_kind.rs:36`、`ReparamCurve2d` → `surface_projector/reparam.rs:10`）。`_igesdump/` 与 `docs/goals/**` 里的历史记录保持原样，不改写已完成目标的留档。2026-10-14 例外：`docs/goals/**` 中指向已拆模块的**路径 token** 已订正（`iges.rs` → `iges/`、`surface_projector.rs` → `surface_projector/`），并给已消失的 `canvases/iges-read-geometry.canvas.tsx` 加了删除标注；各 receipt 的历史数字、命令与结论未改。

**仍待拆（>1000 行，按体量降序）**：`step/read_topology.rs` 2341、`meshing/node_insertion.rs` 1887、`occt-geom/bspline_surface.rs` 1842、`occt-math/function_set_root.rs` 1686、`meshing/incremental_mesh/discret_root.rs` 1533、`occt-geom2d/geom2d_int/int_poly_poly_gen.rs` 1496、`meshing/edge_discret.rs` 1398、`iges_read/iges_to_brep_topo_curve.rs` 1263、`occt-geom/adv_app2var/gauss_gss_xj.rs` 1257、`shhealing/shape_fix_face.rs` 1230。
