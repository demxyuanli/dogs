# occt-topo（BOP/网格/交换以外）移植忠实度审查（2026-09-20）

范围：`crates/occt-topo/src/`。只读，未改源码。
统计：自创 14 条、等价替换未登记 4 条（共 18 条，按严重度排序）。

**结论摘要**：最高危三条 —— ①`brepfeat/features.rs:109` 用解析公式**覆盖布尔结果的体积**；②`wireframe.rs:423+431`（2026-10-06 复核；原 `:399+407` 已漂移）对曲面按**未裁剪自然 UV 域**建网格，`shape_volume`/`shape_surface_area`/`brep_gprop` 全链失真且被 `bop_builder_report` 当修复判据（**✅ 已修 2026-10-14**：A13 换掉建网格的栅格采样，本轮再把 `shape_volume`/`shape_surface_area`/`brep_gprop` 的网格路径整段换成 `BRepGProp` 精确积分，见发现 6/11）；③`brepfeat/features.rs:394`+`boolean_ops.rs:173`（`feature.rs` 已删，订正 2026-10-06） 体素布尔在精确 `bop_builder::boolean`（BOPAlgo_BOP）已存在时仍是 feature 默认路径，注释仍称"exact BRep boolean 不可用"（过期）。另：`hlr.rs`/`viz_scene/`/`draw/`/`xcaf/`/`render_svg.rs` 借用了 OCCT 包名，但不是移植。所有条目均未进入 `export_data_obj`/`step_obj_parity` 主门禁，除第 11、13 条（经 `bop_builder_report`/`shape_mesh` 间接进入；第 6 条已随本轮修复退出该链；第 11、13 条已于 2026-10-14 修复，见各自处置）。

## 发现

### 1. `boss_thru_all` 用解析公式覆盖布尔结果体积
- 判定：自创 ｜ live：**是**（公有 `brepfeat::boss_thru_all`，`brepfeat/tests.rs:171`）
- 证据：`crates/occt-topo/src/brepfeat/features.rs:109`
```rust
result.volume = (v0 + PI * radius * radius * height - overlap).max(0.0);
```
- OCCT：`TKFeat/BRepFeat/BRepFeat_MakePrism.cxx`（`BRepAlgoAPI_Fuse`，体积事后由 `BRepGProp` 积分）；`BRepFeat_MakeRevol.cxx:310-315` 无此步。
- 影响/建议：返回体积与 `result.shape` 实际体积不一致 → `volume_delta`/`ratio`/`kind`（同文件 `:161-179`）全失真，`tests.rs:172` 因覆盖恒真。删覆盖，改用 `brep_gprop_full`；否则标未移植。
- **订正（2026-10-10，已修）**：解析覆盖已删（`1ab71402`）。体积改由 `brepfeat/features.rs` 的 `solid_volume` 按 wire 的存储朝向累加（`c3b3abf9`，OCCT `BRepGProp` 语义），`boss_thru_all_pierces` 实测体积 2.5590491，与 box + tool − overlap 的解析值一致。`tests.rs:172` 仍只断言 `volume > before`，未断言精确体积。

### 2. 网格夹具 `mesh_cylinder`（24 切片）冒充解析圆柱
- 判定：自创 ｜ live：**是**（`translated_mesh_cylinder` ← `boss_thru_all:98`、`boss:123`；`feature.rs:155 z_cylinder_mesh` 为副本）
- 证据：`crates/occt-topo/src/brepfeat/features.rs:419`
```rust
let cyl = mesh_cylinder(radius, height, 24);
```
- OCCT：`TKPrim/BRepPrimAPI_MakeCylinder`（解析圆柱面）；`slices=24` 在 OCCT 中不存在。
- 影响/建议：boss/hole 结果为 24 面棱柱，半径内切误差 `1-cos(π/24)≈0.86%`；`tests.rs:185-187` 的 `r∈(0.1,0.6)` （`:184-188` 弧段）断言正为掩盖棱柱化。改 `primitives::BRepPrimCylinder` 或 `brep_builder_api::make_cylinder_full:263`；不能网格化就标未移植。

### 3. `revolve_profile_about` 用车削网格重建旋转体
- 判定：自创 ｜ live：**是**（`groove:49`、`neck:57` 唯一工具构造路径）
- 证据：`crates/occt-topo/src/brepfeat/features.rs:361`
```rust
pub(super) fn revolve_profile_about(profile: &[GpPnt2d], steps: usize, axis: &GpAx1) -> Result<TopoShape, String> {
```
- OCCT：`BRepFeat_MakeRevol.cxx:215 theRevol.Perform` → `LocOpe_Revol` 产出解析旋转面，无三角网格。
- 影响/建议：`steps`（调用方传 16）决定 groove/neck 几何；`tests.rs:96-109` 的 `<0.5` 容差即网格化后果。按 `TKFeat/LocOpe/LocOpe_Revol.cxx` 移植；不可移植则标未移植。

### 4. `resolution_for` 的 `[16,64]` 夹取
- 判定：自创 ｜ live：**是**（`brepfeat::boolean_feature:334`、`feature.rs:101`；`feature.rs:197` 为同逻辑副本）
- 证据：`crates/occt-topo/src/brepfeat/features.rs:508`
```rust
r.clamp(16, 64)
```
- OCCT：未找到。OCCT 离散化由 deflection 决定，无体素分辨率夹取。
- 影响/建议：`tol` 被吞（`tol<span/64` 恒 64 格，`tol>span/16` 恒 16 格），16/64 为裸常量。标未移植；面分类应走 `BRepClass3d_SolidClassifier`，不做体素。

### 5. `draft` 只旋转面所在平面、不动边界线，且写全局注册表
- 判定：自创 ｜ live：**是**（`phase7_integration.rs:30 brepfeat_draft` 门禁、`brepfeat/tests.rs:56`）
- 证据：`crates/occt-topo/src/brepfeat/features.rs:37`
```rust
GeometryRegistry::global().set_face(
```
- OCCT：`BRepFeat_MakeDPrism::Perform` —— 锥形棱柱 Tool + `BRepAlgoAPI_Fuse/Cut`，侧面由扫掠新建、顶点随之移动；无"改面平面"步骤。
- 影响/建议：面平面被转、顶点仍在原位 → 顶点不在其面上，面/线框自洽性破坏；`tests.rs:80` 只断言面数仍 6 以绕开；全局写入污染共享 `TShape` 的其他视图。按 `MakeDPrism.cxx` 走 Tool+布尔，否则标未移植。

### 6. `face_to_triangles` 对曲面按未裁剪自然 UV 域建网格
- 判定：自创 ｜ live：**是**（`shape_mesh::mesh_shape:29` → `shape_volume`/`shape_surface_area`（`shape_mesh.rs:109,117`）、`brep_gprop.rs:39`、`render_svg.rs:56`、`draw`；`bop_builder_report.rs:405,410` 用作修复判据）
- 证据：`crates/occt-topo/src/wireframe.rs:399`
```rust
let (u0, u1, v0, v1) = face_uv_bounds(f, surface.as_ref());
let nu = ((du / def).ceil() as usize + 1).clamp(3, 64);
```
- OCCT：`BRepMesh` 沿面裁剪边界离散；`BRepGProp_Gauss::Compute`（`BRepGProp.cxx:277`、`BRepGProp_Sinert.cxx:80-109`）为精确积分。
- 影响/建议：90° 圆柱面片段会被当整圈 2π 曲面计面积/体积（`brep_measure.rs:37` 自承 "samples the full natural surface"）。按 `BRepGProp_Gauss` 落地（`brep_gprop_full` 已有基础）。
- **处置（✅ 已修：2026-10-14）**：`face_to_triangles` 本身已由 A13（2026-10-12）改委托 `IncrementalMesh::discretize_face`（`BRepMesh` 忠实管线）；本轮把**测量链**从网格路径摘干净，不再让任何体积/面积经过网格：(1) `shape_mesh::shape_volume`（自创 `(1/6)|Σ a·(b×c)|` 散度和，`shape_mesh.rs:117`）与 `shape_mesh::shape_surface_area` 删除，调用方改 `brep_gprop::{volume,surface_area}`（`brepfeat`/`brep_pattern`/`shape_metrics`/`draw`/`bop_curved`/`phase3_integration`）；(2) `brep_gprop`（见发现 11）改为精确积分；(3) `shape_metrics` 的面积/体积改由 `brep_gprop` 提供。**对拍**（`zz_bop_curved_probe`，deflection 0.02，解析值 = `vol_gauss`）：`cyl r=0.4 h=1` 1.43993 → **0.50265**；`r=0.4 h=2` 2.87986 → **1.00531**；`r=0.4 h=4` 5.75972 → **2.01062**；`r=1 h=1` 1.99957 → **3.14159**；box 8.00000 不变。门禁：`step_obj_gates` 4/4、`export_data_obj` 15/15、`phase3..10/19/20_integration` 与 `bop_builder2_boss` 全绿。

### 7. ~~`feature.rs` 是 `brepfeat` 的重复体素实现 + 过期注释~~（订正 2026-10-06：文件已删）
- 判定：自创 ｜ live：**是**（公有 `protrusion`/`pocket`/`boss`/`hole`；仓库内仅 `feature.rs:217` 自测）
- 证据：`crates/occt-topo/src/feature.rs:101`
```rust
let resolution = resolution_for(solid, tool, tol);
```
- OCCT：`BRepFeat_MakePrism`/`BRepAlgoAPI_Fuse|Cut`。`feature.rs:6-7` 称 "exact BRep boolean 尚不可用"，但 `bop_builder::boolean` 已存在并被 `brepfeat` 使用。
- 影响/建议：同操作产出两条不同几何（体素 vs 平面精确），`warnings` 为假信息。转发 `bop_builder::boolean`（与 `brepfeat/features.rs:319-327` 同构）或并入 `brepfeat` 并标未移植。

### 8. `gprop_analytic` 顶点质量惯量张量兜底
- 判定：自创 ｜ live：**是**（公有 `inertia_tensor`/`full_analytic_properties`；消费方仅 `tests.rs`、`phase10_integration.rs:106`）
- 证据：`crates/occt-topo/src/gprop_analytic/mass_properties.rs:405`
```rust
let m_i = m / verts.len() as f64;
```
- OCCT：`BRepGProp_Vinert`/`VinertGK`（Gauss–Kronrod，任意形状精确），无顶点均摊兜底。
- 影响/建议：非 box/sphere/cyl/cone/torus 实体的惯量继承粗糙估计，API 签名不体现降级。标未移植（注释提了 coarse，但 API 未登记）。

### 9. `shape_metrics` 曲率对角采样 `(i*7 % samples)`
- 判定：自创 ｜ live：**是**（公有 `shape_metrics::shape_metrics`）
- 证据：`crates/occt-topo/src/shape_metrics.rs:146`
```rust
let v = v0 + (v1 - v0) * (i * 7 % samples) as f64 / samples.max(1) as f64;
```
- OCCT：未找到。`BRepGProp`/`BRepLProp` 无此接口。
- 影响/建议：`7` 为裸常量；min/max/avg 只反映对角线 8 点，对整面无代表性（`:148` 又自承有限差分符号有噪声故取绝对值）。标未移植，或改 Gauss 面积 + 真实曲率工具。
- **范围裁定（2026-10-16）**：本项在 T001 清单里登记为自创件残余（D 段），且无对应 OCCT 类可移植（`BRepGProp`/`BRepLProp` 无该接口），故不立包；保持登记、不改代码。同条记录见 `specs/occt-scope-drift.canvas.tsx` 2.4。

### 10. `validate.rs` 不是 `BRepCheck_Analyzer` 的移植
- 判定：自创 ｜ live：**是**（`pub mod validate`，`phase3/5` 集成测试使用；`brep_check.rs` 才是 `BRepCheck` 的移植）
- 证据：`crates/occt-topo/src/validate.rs:69`
```rust
r.messages.push(format!("shell has only {faces} faces (a closed shell needs ≥4)"));
pub fn is_valid(&self) -> bool { true }   // validate.rs:131
```
- OCCT：`BRepCheck_Shell::CheckClosed` 由边沿面使用次数判定闭壳，无"面数≥4"规则；`IsValid` 非常量 true。
- 影响/建议：`Analyzer::is_valid()` 恒真使相关门禁失效；4 面阈值把合法小面片壳判 `Undefined`。删除启发式，转调 `brep_check.rs`；保留 API 标未移植。
- **范围裁定（2026-10-16）**：本项缺的是 `BRepCheck_Analyzer` 这一整层（T001 清单 C-B：需先移植未移植层），`brep_check.rs` 是另一条 `BRepCheck_Edge` 移植线；本 tranche 未立包，保持登记、不改代码。使用约束照旧：`validate.rs` 不得作门禁（清单 F）。同条记录见 `specs/occt-scope-drift.canvas.tsx` 2.4。

### 11. `brep_gprop.rs` 是网格采样，不是 `BRepGProp`
- 判定：自创（模块头自承 "Approximate port"）｜ live：**是**（`bop_builder_report.rs:405,410` 修复前后面积判据、`brep_pipe.rs:77`、`loft.rs:386`、`gprop_analytic/analytic_props.rs:418+`）
- 证据：`crates/occt-topo/src/brep_gprop.rs:98`
```rust
let n = ((1.0 / deflection.max(1e-4)).ceil() as usize).clamp(1, 32);
```
- OCCT：`BRepGProp::SurfaceProperties`/`VolumeProperties` → `BRepGProp_Gauss`。
- 影响/建议：`1/deflection` 夹 `[1,32]` 是裸规则，面积/体积/质心随 deflection 漂移，而 OCCT 有 `myEpsilon` 误差估计。转发 `brep_gprop_full`。
- **处置（✅ 已修：2026-10-14）**：本地网格 `mesh_shape`（`brep_extrema::mesh_faces(shape, n, n)` 的 `1/deflection` 夹 `[1,32]` 栅格）与 `compute_inertia` 散度和整段删除；`shape_mass_properties` 改走 `brep_gprop_full::{volume_properties, surface_properties}`（`BRepGProp_Vinert`/`Sinert` + `BRepGProp_Gauss`，精确，无 deflection）。唯一保留的网格分支是 OCCT `BRepGProp.cxx:344-360` 的 `NoSurf` 臂：解析积分无法建立的面（`FaceGauss::new` 失败/无面）改用 `brep_gprop_full::mesh_props`（`BRepGProp_MeshProps`）在 `shape_mesh::mesh_shape` 的 `BRepMesh` 三角化上积分，`deflection` 只在该回退里生效。公开签名（`volume`/`surface_area`/`centroid`/`solid_density_mass`/`shape_mass_properties`）未变，调用方零改动。

### 12. `offset_by_sampling` 把一般曲面采成 16×16 双线性面（✅ 2026-10-14 已修）
- 判定：自创 ｜ live：原为**是**（`offset_face_surface:644` 的 cone/torus/其他分支 → `offset_face` → `offset_shell`，即 `draw` 的 offset 命令）
- 证据：`crates/occt-topo/src/brep_offset/curve_face_offset.rs`（原 `:588`）
```rust
let (nu, nv) = (16usize, 16usize);
```
- OCCT：`BRepOffset_Offset`/`BRepOffsetAPI_MakeOffsetShape` 构造精确 `Geom_OffsetSurface`。
- 影响/建议：cone/torus 偏移面退化为 16×16 分段双线性面；原 `curve_face_offset.rs:648-649` 还自承非平面面"边不投影到偏移面"。改用 `Geom_OffsetSurface`，否则标未移植，不返回近似面。
- **处置（✅ 已修：2026-10-14）**：
  - 新增 `brep_offset/brep_offset_surface.rs::brep_offset_surface(s, offset, allow_c0) -> (Box<dyn Surface>, BRepOffsetStatus)`，逐分支移植 `BRepOffset::Surface`（`BRepOffset.cxx:44-200`）：`Geom_RectangularTrimmedSurface` 递归基面并保留 `(u1,u2,v1,v2)` 窗口（`:180-191`）；Plane 沿 `x_dir ^ y_dir * Offset` 平移（`:54-61`）；Cylinder 半径 `±Offset`（`r <= -Confusion` → 轴旋转 π、status `Reversed`；`|r| < Confusion` → `Degenerated`，`:62-90`）；Cone `Radius = RefRadius + Offset*cos(Alpha)` 且参考面沿轴移 `-Offset*sin(Alpha)`（`:92-111`）；Sphere（`:113-141`）；Torus 仅 spindle（`MinorRadius < MajorRadius`）有解析偏移（`:143-171`）。其余类型（`SurfaceOfRevolution`/`SurfaceOfLinearExtrusion`/`BSplineSurface`/`BezierSurface`/`OffsetSurface`）落到精确 `Geom_OffsetSurface`（`:193-197`）。`BRepOffset_Status` 按 `BRepOffset_Status.hxx:24-31` 定序。无任何采样拟合。
  - `offset_face`（`curve_face_offset.rs:513`）删除 `offset_by_sampling`/`offset_plane_surface`/`extract_cylinder`/`is_cylindrical`/`is_spherical`，改直接调 `brep_offset_surface`（与 `BRepOffset_Offset::Init` 同一调用）。
  - **未移植（已在函数头登记）**：OCCT 一般臂按源 pcurve 重建每条边（`BRepOffset_Offset.cxx:911-958`、写回 pcurve `:1076-1078`、`ComputeCurve3d` `:1080-1083`）尚未移植，故无解析偏移的面（planar B-spline、cone 等）仍保留源 wire。`allowC0` 仅 OCCT 传给 `Geom_OffsetSurface` 构造（`:196`），端口 `GeomOffsetSurface::new` 恒建 osculating 探针（等价 `allowC0 = false`），故只接受 `false`。
  - **对拍**：OCCT 8.0.0 探针 `occt_probe --offsetsurf <dist>`（固定 `u∈{0.4,1.1,2.7}`、`v∈{-0.5,0.3,1.2}` 网格，先按 `BRepOffset_Offset::Init` 拆 trim）对照端口 `examples/zz_offset_surf.rs`：`Cone`/`Cylinder`/`Sphere`/`Torus`/`rev`/`Extrusion`/`Offset`/`Cube`/`Shape-1`/`Shape-2`/`screw`/`OffsetPlaneHoleEdge` × `{0.5,-0.3,2.0}` **36/36 逐位一致**（含类型、`BRepOffset_Status`、采样点）。`HoledPlate` 3/3 仅面枚举次序不同（按 `f` 去掉序号后结果集合完全相同）。一般面（`OffsetSurface` 基面如 `Shape` f=9；负半径 torus 如 `linkrods` -0.3）两侧同为 `Geom_OffsetSurface`、status 相同，取值差异属 `occt-geom` 的 `GeomOffsetSurface::Value` 子系统（已记 `specs/_audit/occt-geom.md:180`）。`linkrods` 0.5/2.0 探针在 f=13 求值越界参数（u=2.7 远超该面 UV 窗口）时 OCCT 侧崩（`0xC0000409`），无 f≥13 的 OCCT 真值。
  - **门禁**：`cargo test -p occt-topo --lib brep_offset` **12 passed / 0 failed**；`step_obj_gates` **4/4**。

### 13. `edge_length` 是等距折线采样（✅ 2026-10-14 已修）
- 判定：自创 ｜ live：**是**（`bop_builder_report.rs:253`、`bop_builder_repair.rs:634`、`algo_tools/construct.rs:546`、`shhealing/common.rs:492` 作短边判据；`shape_metrics.rs:40`）
- 证据：`crates/occt-topo/src/brep_measure.rs:12`
```rust
pub fn edge_length(e: &Edge, samples: usize) -> f64 {
```
- OCCT：`BRepGProp::LinearProperties`（`BRepGProp.cxx:121`）→ `BRepGProp_EdgeTool` 积分，带误差估计。
- 影响/建议：采样数由调用方给（8/32），圆弧短边在 8 段折线下系统性低估 → 短边判据结果被替换规则左右。改调 `BRepGProp` 线性积分（`brep_gprop_full`）。
- **处置（✅ 已修：2026-10-14）**：
  - 新增 `brep_gprop_full::edge_linear_properties(e)`（`brep_gprop_full/properties.rs:107`）：单边走 `BRepGProp_Cinert::Perform`（`cinert_perform`），积分阶数取 `curve_integration_order` = `BRepGProp_EdgeTool::IntegrationOrder`（直线 2、其余 10），朝反向边按 `TopoDS_Shape::Orientation` 反转曲线与参数区间——与 `linear_properties` 的逐边循环同一段代码。`GProps::dim` 即边长。
  - `brep_measure::edge_length(e, samples)` → `edge_length(e)`（`brep_measure.rs:15`）：删除 `samples` 形参与折线循环，改为委托 `edge_linear_properties`；同类聚合 `wire_length`/`shape_perimeter`/`max_edge_length`/`min_edge_length` 一并去掉 `samples`。调用方同步为 `bop_builder_report.rs:253`、`bop_builder_repair.rs:638`、`algo_tools/construct.rs:577`（`is_degenerated_edge`）、`shape_metrics.rs:40,252,256`、`brep_faces.rs:207`（`face_wire_length`）、`gprop_analytic/edge_length.rs`（`CurveKind::Other`/`circle_radius=None` 的兜底由折线改为同一 Gauss 积分）与三处测试。
  - **对拍**：`cargo test -p occt-topo` **1317 passed / 0 failed**（1 ignored）；`step_obj_gates` **4/4**；`export_data_obj` **15/15** 且 obj 顶点/面数与基线逐位一致（本次改动改变的是短边判据的数值，未改变任何门禁几何）。

### 14. 用采样识别二次曲面代替 downcast（✅ 已全部收口：2026-10-15）
- 判定：等价替换未登记（后果改变）｜ live：**否**（原 `offset_face_surface:613,619,632` 的 `is_cylindrical`/`is_spherical` 分派已随 §12 删除；`offset_face` 的 `is_planar(new_surf, 8, 8, 1e-6)` 与 `sphere_center` 估判已改为 `gp_pln()`/`gp_sphere()` 精确判据；`brep_surface::is_planar`/`face_is_planar` 采样版及全部调用点已删除；`edge_face_kind` 的曲线/曲面族与 `IntTools_EdgeFace` 的 `IsEqDistance`/`IsCoplanar`/`IsRadius` 亦已换成 `GetType()`/`gp_*` 同等 downcast）
- 证据（原 `crates/occt-topo/src/brep_offset/curve_face_offset.rs:551`，已删除）
```rust
pub(super) fn is_cylindrical(s: &dyn Surface, loc: &GpPnt, dir: &GpDir, r: f64) -> bool {
```
- OCCT：`BRepOffset_Offset` 用 `Geom_Cylinder`/`IsKind(STANDARD_TYPE(...))` 类型判据，无采样。
- 影响/建议：小曲率 B 样条等可被误判为圆柱 → 走错偏移分支。`tgeometry.rs` 已声明 `dyn Surface` 无法 downcast，但本模块未登记。
- **现状（2026-10-14）**：§12 删除了 `offset_face` 里采样判圆柱/球的分派与 `offset_by_sampling`；`brep_surface::classify_surface`（`brep_surface.rs:89`，`GeomAdaptor_Surface::Load` 等价）已存在并经审核 A30 落地。`offset_face`（`curve_face_offset.rs:520`）选分支仍用 `is_planar(new_surf, 8, 8, 1e-6)` 而非 `classify_surface(...) == Plane`，`sphere_center`（`brep_surface.rs:174`）亦非 `Geom_Sphere` downcast。残留待改，未登记为已对齐。
- **处置（✅ 已修：2026-10-15）**：
  - `brep_surface::sphere_center`（`brep_surface.rs:136`）改为精确 downcast：`s.gp_sphere().map(|sp| sp.pos.location())`，即 `GeomAdaptor_Surface::Sphere()`（`GeomAdaptor_Surface.cxx:1997-2006`）的 `Geom_SphericalSurface` 判据；`gp_sphere` 已按 `GeomAdaptor_Surface.cxx:423-425` 递归 `Geom_RectangularTrimmedSurface` 基面（`occt-geom/src/rectangular_trimmed.rs:609-611`），`Geom_OffsetSurface` 不作答（`GeomAdaptor_Surface.cxx:496-523` 给 `GeomAbs_OffsetSurface`，非 `GeomAbs_Sphere`）。原实现取三个采样点的外接圆圆心：对真球精确，但对柱/锥/环/**任意 B 样条片**同样返回 `Some`，调用方据 `sphere_center(s).is_some()` 把 B 样条判成球（`edge_face_kind.rs:207` 的 `surface_kind` 是最后兜底臂，`topo-rest.md` 各处球分支同理）。`det3`/`solve3`/`circumcenter`（原 `brep_surface.rs:121-167`）随之成为死代码，整段删除。
  - `offset_face`（`curve_face_offset.rs:517`）两条分支改为精确面型判据：平面臂 `if let Some(p) = new_surf.gp_pln()`（`:524`），法向直接取 `gp_Pln::Position().Direction()`，与 `BRepOffset::Surface` 的平移向量 `XDirection ^ YDirection * Offset`（`BRepOffset.cxx:54-61`）同源；球臂改用 `(surf.gp_sphere(), new_surf.gp_sphere())`（`:539`），半径取 `gp_Sphere::Radius()` 而非 `d0(0,0).distance(center)`（原写法在裁剪球上会求值到 UV 窗口外）。原 `is_planar(new_surf, 8, 8, 1e-6)` 会把**几何上平的 B 样条面**判成平面臂（OCCT 给 `GeomAbs_BSplineSurface`，落一般臂），而 `sphere_center` 的 `Some` 兜底会把**锥/柱面**拖进球臂做一次以伪圆心为中心的 `SetScale`——面边界被错误缩放。修后只有真平面/真球进入这两臂，其余落 `keep source wires` 臂（该臂本身仍是 UNPORTED 的一般臂，见 `curve_face_offset.rs:508-516` 函数头）。
  - `brep_offset/mod.rs:15` 的 prelude 去掉已无调用方的 `is_planar`/`surface_normal`，`sphere_center` 改为 testing-only 直接引入（`brep_offset/tests.rs:9`）。
  - **批次 2（✅ 已修：2026-10-15）**：BOP 交叉与 EF 家族 9 处调用点按 `.cxx` 判据换成精确面型判据，新增 `brep_surface::is_plane_type`（`brep_surface.rs:107`，即 `GeomAdaptor_Surface::GetType() == GeomAbs_Plane`）与面级包装 `brep_surface::face_is_plane_type`（`brep_surface.rs:112`）：
    - `pave_ff_perform.rs:401`（原 `:372`）与 `pave_intersect/face_face.rs:86`（原 `:83`）：`BOPAlgo_PaveFiller_6.cxx:380` 的 `aBAS1.GetType() == GeomAbs_Plane && aBAS2.GetType() == GeomAbs_Plane`，它门控 `CheckPlanes`。原采样版会把**几何上平的 B 样条面**当平面，从而在两张 B 样条平面片之间错误地走 `Append InterfFF + Init(0,0)`（跳过 FaceFace），OCCT 会真做相交。
    - `pave_ff_perform.rs:471-472`（原 `:442-443`）与 `pave_intersect/face_face.rs:165`（原 `:159`）：切向容差的平面臂是 `IntTools_FaceFace::PerformPlanes`（`IntTools_FaceFace.cxx:2541-2555` 的 `sqrt(aDt^2 + TolF1^2)`），由 `Perform` 在 `aType1 == GeomAbs_Plane && aType2 == GeomAbs_Plane`（`cxx:395`）时进入。`int_tools_curve_box::curve_tangential_tolerance` 的文档注释原写 `ComputeTolReached3d`，实为 `PerformPlanes`，一并纠正（`int_tools_curve_box.rs:182-191`）。
    - `pave_ef_perform.rs:229`：`bLinePlane`，`BOPAlgo_PaveFiller_5.cxx:392-393` 的 `aBAC.GetType() == GeomAbs_Line && myContext->SurfaceAdaptor(aF).GetType() == GeomAbs_Plane`。
    - `pave_force_ef.rs:251`：`BOPAlgo_PaveFiller_5.cxx:1038` 的 `aSurfAdaptor.GetType() != GeomAbs_Plane || aBAC.GetType() != GeomAbs_Line`（`bUseAddTol` 前置）。
    - `pave_pcurves.rs:201` `is_based_on_plane`：`BOPAlgo_PaveFiller_7.cxx:936` 的 `IsBasedOnPlane`——`down_cast<Geom_RectangularTrimmedSurface>` 后取基面再 `down_cast<Geom_Plane>`，与 `GetType() == GeomAbs_Plane` 同解（`classify_surface` 递归裁剪基面，`gp_pln` 亦递归）。
    - `edge_face.rs:407`：`IntTools_EdgeFace::DistanceFunction` 走 `GeomAPI_ProjectPointOnSurf` → `Extrema_ExtPS::Perform` 的 `case GeomAbs_Plane`（`Extrema_ExtPS.cxx:279-281`，解析 `Extrema_ExtPElS`）；采样版把平的 B 样条面也送进解析投影。
  - `pave_ff_perform::is_plane_ff`（`pave_ff_perform.rs:84`）原注释写「Geom_Plane, or plane under offset / trim」但实现只有 `classify_surface == Plane`（不含 offset）。本批按 `BOPAlgo_PaveFiller_6.cxx:84-102` 重写为真正的 `IsPlaneFF`（平面 / 裁剪基面为平面 / 偏移基面为平面），并让 `shift_closed_faces` 恢复 OCCT 的双判据：`:404` 的跳过用 `GetType()`，`:410,416` 传给 `IsClosedFF` 的 `anIsPlane1/2` 用 `IsPlaneFF`（`pave_ff_perform.rs:241-256`）。
  - **批次 3（✅ 已修：2026-10-15，`is_planar`/`face_is_planar` 采样版整段删除）**：批次 2 末项列出的余 19 处调用点全部按各自 `.cxx` 判据换成 `is_plane_type`/`face_is_plane_type`，随后 `brep_surface::is_planar`（原 `:61`）与 `brep_surface::face_is_planar`（原 `:188`）及 8×8 法向采样体整段删除，仓内已无 `is_planar` 符号（仅注释保留迁移说明）：
    - `brep_surface.rs:153`（`face_plane`）、`pcurve.rs:249`（`make_pcurve_on_face` 的 `ProjLib::IsAnaSurf` 平面臂，`ProjLib.cxx:217-232`、`BOPTools_AlgoTools2D.cxx:501-560`；`pcurve.rs:391` 的 `face_projector` 走 `classify_surface == Plane`）。
    - `fillet_curved/surface_info.rs:139,220`（`is_curved_surface` / `plane_from_surface`：`ChFi3d` 用 `GeomAbs_Plane`）、`fillet_curved/torus_blend.rs:603`（`face_is_planar_face`）。
    - `gprop_analytic/analytic_props.rs:378,587`、`gprop_analytic/mass_properties.rs:17,61,325,443`：质量属性积分路径本身不含采样语义，原先只借用采样版判「全体平面」，现直接用 `GetType()`。
    - `brep_offset/shell_offset.rs:122`、`brep_faces.rs:203`（`face_is_planar_face` → `face_is_plane_type`）。
    - 测试夹具同批改口径：`brep_builder_full/mod.rs:247`、`brep_sketch.rs:234`、`brep_builder_api.rs:372,405-406`、`fclass2d/tests.rs:160`、`brep_offset/tests.rs:132`、`fillet_curved/tests/helpers.rs:693`、`brep_surface.rs:288,304`。
  - **批次 3 同族（✅ 已修：2026-10-15）：`edge_face_kind.rs` 从采样族改为 `GetType()` 同等 downcast**，全部按 `IntTools_EdgeFace` 的分支落位（模块头注明 `IntTools_EdgeFace.cxx:74,247,532,620,650,798,827`）：
    - `curve_kind`（`:59`）＝ `Adaptor3d_Curve::GetType()` 的 Line/Circle/BSpline 三臂：`gp_line()`（`GeomAdaptor_Curve.cxx:252-311` 的 `GeomAbs_Line`）、`gp_circ()`（`GeomAbs_Circle`）、`nurbs_degree()`（`GeomAbs_BezierCurve`/`GeomAbs_BSplineCurve`），其余走 `default`（`IntTools_EdgeFace` 视为既非线、非圆、亦非 B 样条）。原先采 8 个点做共线/等距判定，会把**几何上是直线或圆的 B 样条**报成 `Line`/`Circle` 并送进解析臂。
    - `circle_geometry`（`:34`）＝ `CurveTool::Circle`（`IntTools_EdgeFace.cxx:800,829`）的 `gp_circ()` downcast；原先用 8 个采样点拟合圆，B 样条圆同样作答。
    - `surface_kind`（`:207`）直接委托 `classify_surface`（`Adaptor3d_Surface::GetType()`）；原先分两段采样（从 `d0` 反推二次曲面参数、再用 8×8 栅格与自造 `1e-4` 容差复核），会把形似柱/锥/环的 B 样条片报成解析型，也会因 UV 窗口漂移而否掉真二次曲面。
    - `is_eq_distance`（`:217`）＝ `IntTools_EdgeFace::IsEqDistance`（`cxx:240-299`）：分支取 `GetType()`，frame 依次取 `gp_cylinder()`/`gp_cone()`/`gp_torus()`（`SurfaceTool::Cylinder/Cone/Torus`），容差用调用方传入的字面 `1e-7`；锥顶按 `gp_Cone::Apex()`（`Location - V * (RefRadius / tan(SemiAngle))`）解析求出，不再靠采样点列。
    - `is_coplanar`（`:271`）/`is_radius`（`:284`）＝ `IsCoplanar`/`IsRadius`（`cxx:787-811`、`cxx:816-840`）：`aCType == GeomAbs_Circle && aSType == GeomAbs_Plane` 双门控，轴向比用 `gp_Pln::Axis().Direction()`（`SurfaceTool::Plane`），`aPln.Distance(aCenter)` 取自平面自身 placement 而非采样 `d0`。
    - 随之删除的采样/拟合件：`is_line_like`、`is_circle_like`、`first_three_spanning`、`det3`、`solve3`、`circumcenter`、`cylinder_matches`、`cone_matches`、`torus_matches`、`quadric_u0`、`midpoint`（均在原 `edge_face_kind.rs` 内）。
  - **同族残留（✅ 已修：2026-10-15，批次 4）**：`edge_face_kind::curve_resolution`（原 `:73`）当时仍是 64 段差分估计最大速率再取 `tol / max_speed`，而 OCCT 的调用点（`IntTools_EdgeFace.cxx:326,327,373`）走的是 `myC.Resolution(myCriteria)`，即 `BRepAdaptor_Curve::Resolution` → `GeomAdaptor_Curve::Resolution`（`GeomAdaptor_Curve.cxx:1116-1148`）的**解析闭式**：Line→`R3D`、Circle→`2*asin(R3D/(2R))`、Ellipse→`R3D/MajorRadius`、Bezier/BSpline→`Geom_BezierCurve::Resolution`/`Geom_BSplineCurve::Resolution`（→`BSplCLib::Resolution`）、default→`Precision::Parametric(R3D)`。该闭式在 `occt-geom` 已作为 `Curve::resolution(r3d)` 逐型落地（`curve.rs:135` 的 default、`line.rs`、`circle.rs:37`、`ellipse.rs`、`bezier_curve.rs:636`、`bspline_curve.rs:984`）。批次 4 把该族全部换成 `Curve::resolution` / `Curve2d::resolution`，明细见下方「批次 4」。

  #### 批次 4（✅ 已修：2026-10-15）：`Resolution` 采样族全量换成 `Adaptor::Resolution`

  - **`occt-geom2d` 新增 `Curve2d::resolution`**，即 `Geom2dAdaptor_Curve::Resolution`（`Geom2dAdaptor_Curve.cxx:1186-1219`）：trait default＝`Precision::Parametric(Ruv)`（该函数 `default:` 臂，落地 `OffsetCurve`/`Parabola`/`Hyperbola`/`OtherCurve`），逐型覆盖 `line.rs`（`Ruv`）、`circle.rs`（`R > Ruv/2` 时 `2*asin(Ruv/(2R))`，否则 `2*pi`）、`ellipse.rs`（`Ruv/MajorRadius`）、`bezier_curve.rs`（`Geom2d_BezierCurve::Resolution` `cxx:743-760`，poles 抬到 `z=0` 后调共享的 `BSplCLib::Resolution`，knot 序列为 `deg+1` 个 0 接 `deg+1` 个 1）、`bspline_curve.rs`（转调既有 inherent `Geom2dBSplineCurve::resolution`，`bspline_curve.rs:724`）、`trimmed.rs`（转发基曲线，对应 `Load` 对 `Geom2d_TrimmedCurve` 的递归 `cxx:284-288`；本移植的 2D 裁剪保留基参数空间，故无需再缩放）。原先全仓 `Curve2d` 无此查询，调用方各自用 `tol/|D1|` 或中点差分估计。
  - **调用点替换**（每处都回各自 `.cxx` 确认过）：
    - `edge_face_kind::curve_resolution`（`:73`）→ `c.resolution(tol)`，调用方 `edge_face_type.rs:26,59` 的 `myC.Resolution(myCriteria)`（`IntTools_EdgeFace.cxx:326-327` `MakeType`、`:373` `CheckTouch`）。
    - `inttools_range`（原 `curve_resolution` 三处调用）→ `curve.resolution(..)`：`IntTools_ShrunkRange.cxx:162` 的 `aPTolE = aBAC.Resolution(aTolE)`（后面那条 `min(aPTolE, range/100)` 保留）、`BRepLib_1.cxx:202` 的 `anEps = max(theCurve.Resolution(theTolE)*0.1, ...)`、`BRepLib_1.cxx:61` 的 `aStep = theCurve.Resolution(theTol)*1.01`。原 8 点中心差分求最大速率的函数体整段删除（`find_valid_range` 里对 `[first,last]` 的依赖随之消失，`Resolution` 与范围无关）。
    - `algo_tools_range::correct_range_ee` / `correct_range_ef` → `BOPTools_AlgoTools_2.cxx:299-341` / `:385-415`：早退改 `c.is_line()`（`aCT == GeomAbs_Line`），freeform 判据改 `aCT == GeomAbs_BezierCurve || BSplineCurve || OffsetCurve || OtherCurve`（`nurbs_degree()` / `offset_curve()` / 其余落 `OtherCurve`），`else` 臂（Circle/Ellipse/Parabola/Hyperbola）改 `c.resolution(a_res)`。原 `curve_is_line` 用 3 点共线采样、`curve_is_freeform = !curve_is_line` 把四类二次曲线也拖进 D1 臂——OCCT 对这些走 `aBC.Resolution`。
    - `bean_face::resolution`（`:206`）→ `self.curve().resolution(tol)`，即 `myCurveResolution = myCurve.Resolution(myCriteria)`（`IntTools_BeanFaceIntersector.cxx:178,200,223,583` 写、`:998,1189,1198,1897` 读）。
    - `wire_splitter_block::angle_2d`（`:489`）→ `dt = c.resolution(tol2d).max(PCONFUSION)`，非直线臂改 `!c.is_line()`（`aGAC2D.GetType() != GeomAbs_Line`，`BOPAlgo_WireSplitter_1.cxx:792,795`），替掉原先的 `pc_curve_kind(..) != CurveKind::Line` 周期/无界启发式；`curve_resolution` 保留为一行转调。该 `pc_curve_kind` 分类器本身随后由 T-132（2026-10-15）整段删除，全仓改走 `occt_geom2d::adaptor_curve::get_type2d`（`Geom2dAdaptor_Curve::GetType` 同等件）。
    - `int_tools_vertex_line::adaptor_resolution`（`:415`）→ `curve.resolution(tol3d)`，即 `aPTol = aGAC.Resolution(max(aTolR3D, aTolV))`（`BOPAlgo_PaveFiller_6.cxx:3001-3002`）；原先对无界曲线短路返回 `tol3d`（OCCT 的 `Resolution` 只看 `GetType()`，无界直线照样返回 `Ruv`）。
    - **装饰器补 `resolution` 以免退回 trait default**：`pcurve_full/surface_projector/reparam.rs` 的 `ReparamCurve2d`（`GeomLib::SameRange` 的线性重参数化替身；直线/圆臂只平移/旋转故 `scale=±1`，其余臂 OCCT 转 BSpline 并 `BSplCLib::Reparametrize`，参数增量按 `scale` 折算）；`pcurve.rs` 的 `PlaneKeepParam2d`（`GeomProjLib::ProjectOnPlane` KeepParam=true，参数空间与源曲线一致，按该 wrapper 既有 kind 委托惯例取源曲线答案）。加之前这两类 wrapper 会落到 `Precision::Parametric`，把步长缩小若干量级。
  - **门禁（批次 4）**：`cargo check` 五个 crate 与 `engine/occt/runner` 无 error（`engine/occt/driver` 的 build script 因缺 C++ 工具链失败，与本次改动无关）；`cargo test -p occt-geom2d` **72 passed / 0 failed**；`cargo test -p occt-topo --lib` **1251 passed / 0 failed**；`--tests` **1317 passed / 0 failed**（`step_obj_gates` **4/4**）；`export_data_obj` **15/15**，obj 顶点/面数与包围盒与基线逐位一致（`Cone` v=195 f=301、`Cube` v=24 f=12、`Cylinder` v=146 f=140、`Extrusion` v=24 f=12、`HoledPlate` v=180 f=128、`Offset` v=712 f=892、`OffsetPlaneHoleEdge` v=48 f=32、`Shape-1` v=3335 f=4328、`Shape-2` v=3081 f=4744、`Shape` v=6150 f=11372、`Sphere` v=642 f=1244、`Torus` v=1369 f=2592、`linkrods` v=3494 f=5078、`rev` v=104 f=92、`screw` v=601 f=792）。
  - **同批发现、未改**：`inttools_range::find_nearest_valid_point` 仍缺 `BRepLib_1.cxx:69-83` 的 `aD1Mag` 局部奇点加速分支（`GetType()` 判 Bezier/BSpline、Offset 则先下到基曲线，`aD1Mag = 1/Resolution(1.)*0.01` 再平方）；该缺口在原注释里已标为有意省略（`ponytail:`）。**已于 T-133（2026-10-15）收口**：`aD1Mag` 阈值（`cxx:69-81`）与 `D1` 加速倍增循环（`cxx:108-141`）全部落地，判型复用 `edge_edge::curve_kind`＝`GeomAdaptor_Curve::GetType()`（`GeomAdaptor_Curve.cxx:252-311`），`ponytail:` 注记删除；门禁 `--lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 且 v/f 逐位一致。

  - **门禁（批次 2）**：`cargo test -p occt-topo --lib` **1251 passed / 0 failed**；`cargo test -p occt-topo --tests` **1317 passed / 0 failed**；`step_obj_gates` **4/4**，23 模型字号与面积比全 `1.0000`（`our/occ`）；隔离复算（`git stash` 掉本批改动）得到逐位相同的 `v/f` 与面积，证明本批对网格无位移。
  - **门禁（批次 3）**：`cargo test -p occt-topo --lib` **1251 passed / 0 failed**；`cargo test -p occt-topo --tests` 全部 `ok`（`step_obj_gates` **4/4**、`phase3..10/19/20_integration` 8/5/5/4/9/7/5/5/5/8、`bop_builder2_boss` 2）；`export_data_obj` **15/15**，obj 顶点/面数与基线逐位一致（`Sphere.step` v=642 f=1244、`Torus.step` v=1369 f=2592、`Offset.step` v=712 f=892、`Shape.step` v=6150 f=11372、`Cylinder.step` v=146 f=140、`rev.step` v=104 f=92 等均未变）。

### 15. `hlr.rs` 是画家算法近似，不是 `HLRBRep`
- 判定：自创 ｜ live：**是**（`render_svg.rs:57`、`draw` 的 view 命令）
- 证据：`crates/occt-topo/src/hlr.rs:126`
```rust
/// ponytail: conservative segment test — a partial occlusion (an edge crosses
```
- OCCT：`TKHLR/HLRBRep/**`（`HLRBRep_Algo`/`Hider`/`EdgeBuilder`）基于精确边面求交与轮廓线。
- 影响/建议：凹形、自遮挡、轮廓线的可见边集合与 OCCT 不同。标未移植（模块名与 OCCT 包同名，易误认为已移植，应在 `lib.rs:59` 注明）。
- **范围裁定（2026-10-16）**：本项归 T001 清单 C-C（lite 替代，`hlr`），与 `viz_scene`/`draw`/`xcaf`/`render_svg` 同类；本 tranche 不立包，保持登记、不改代码。同条记录见 `specs/occt-scope-drift.canvas.tsx` 2.4。

### 16. `viz_scene/`、`draw/`、`xcaf/`、`render_svg.rs` 借用 OCCT 包名但实现自写
- 判定：自创 ｜ live：**是**（均 `pub mod`；`draw::run_script` 有 doctest）
- 证据：`crates/occt-topo/src/viz_scene/mod.rs:3`、`draw/mod.rs:269`、`xcaf/mod.rs:6`
```rust
//! `V3d_View` ... for offline rendering.
//! * Boolean operations route through [`curved_boolean_ext`] ...
//! XCAF ... embedded as a comment block before the `DATA` section
```
- OCCT：`TKV3d` 走 OpenGL 驱动；`Draw`（`src/Draw/TKTopTest`）是 Tcl 解释器 + 命令注册；`XCAFDoc` 需 STEP 214 AP 的 XCAF 结构。
- 影响/建议：`draw` 语法（`set`/`for`/`if`）与 OCCT Draw 命令不兼容；`xcaf` 元数据塞在 STEP 注释块，OCCT 读不到。三处在模块头统一写"未移植：本模块不是 OCCT X 的移植，是替代实现"。

### 17. `pipe_along_polyline` 只平移剖面、不建移动标架
- 判定：自创 ｜ live：**否**（仅 `brep_pipe.rs:122+` 自测；`draw` 未暴露 pipe 命令）
- 证据：`crates/occt-topo/src/brep_pipe.rs:50`
```rust
let sections: Vec<Wire> = path
```
- OCCT：`BRepOffsetAPI_MakePipe` → `BRepFill_Pipe::Perform`（`TKBool/BRepFill/BRepFill_Pipe.cxx:193`）用 `GeomFill_LocationLaw`（Frenet/CorrectedFrenet/DiscreteTrihedron）旋转剖面。
- 影响/建议：非直线路径上剖面朝向错误，而 `brep_pipe.rs:1-5` 写 "Port of `BRepOffsetAPI_MakePipe`"。按 `BRepFill_Pipe.cxx:193-260` 建标架，否则标未移植。
- **范围裁定（2026-10-16）**：本项确定对应件是 `BRepFill_Pipe::Perform`（未移植层），且 `live：**否**`（仅自测，`draw` 未暴露 pipe 命令），故本 tranche 不立包，保持登记、不改代码。同条记录见 `specs/occt-scope-drift.canvas.tsx` 2.4。

### 18. `solid_union::shape_bbox` 9×9 采样兜底 + `brep_class3d` 确定性探针
- 判定：等价替换未登记 ｜ live：**部分**（`solid_union` 本体仅 `bop_builder_tests.rs:174` 对拍；`bbox_union` 被 `boolean_ops.rs:184` 用；`brep_class3d.rs:96 find_a_point_in_the_face` 供分类器使用）
- 证据：`crates/occt-topo/src/solid_union.rs:103`
```rust
for i in 0..=8 {
```
- OCCT：`BRepBndLib::Add`/`BndLib` 解析极值；`BRepClass3d_SolidExplorer::FindAPointInTheFace`（`BRepClass3d_SolidExplorer.cxx:512-560`）逐面取参数点并带 `aTestInvert` 反转重试。
- 影响/建议：9×9 格点会漏极值 → 体素网格范围偏小；探针方向固定，不等价于 OCCT 多方向重试。包围盒改 `geom_bnd_lib_*` 解析路径；`OtherSegment` 补 `aTestInvert` 分支并登记简化。

### 19. `bbox_from_geometry::shape_bbox` 用 64 点边折线 + 16×16 UV 栅格代替解析极值（✅ 2026-10-15 已修，T-134）
- 判定：自创（采样代替解析极值/闭式） ｜ live：**原为是**（`bop_builder.rs:86`、`bop_builder_dispatch.rs:109`、`bop_classify_occt.rs:78`、`builder_solid.rs:128` 的 BOP 预筛与网格自适应偏转）
- 证据：`crates/occt-topo/src/bbox_from_geometry.rs:56-74`（已删除）
```rust
for i in 0..=64 {
    let u = a0 + (a1 - a0) * i as f64 / 64.0;
    b.add_point(&curve.d0(u));
}
```
- OCCT：`BRepBndLib::Add(S, B, /*useTriangulation*/ false)`（`BRepBndLib.cxx:81-215`）：面走 `BndLib_AddSurface::Add(BS, Tolerance(F), B)`（`cxx:109,118`，仅平面面再加 `BndLib_Add3dCurve::Add` 于其边界边 `:128` 并 `B.Enlarge(Tolerance(F))` `:131`），面上/外的边与顶点分别走 `:146-149`、`:212`，全部取 `geom_bnd_lib_*` 的解析或控制网极值。
- 影响/建议：64/16×16 栅格会漏极值 → BOP 预筛盒偏小（可能把真相交的形状判为不相交），且该盒同时是网格自适应偏转的尺度基准；且本项无任何审计文档登记（原由 T001 记为 DOCDRIFT）。
- **处置（✅ 已修：2026-10-15）**：
  - `bbox_from_geometry::shape_bbox` 整体委托既有忠实件 `brep_bnd_lib::shape_bnd_box`（= `BRepBndLib::Add(S, B, false)`），采样版 `add_to_bbox` 递归体（顶点/边/面三分支 + 子形状递归）与它唯一的 `face_uv_bounds` 依赖一并删除；三个入口 `shape_bbox` / `shape_bbox_with_margin` / `is_shape_outside` 签名不变，BOP/网格调用点零改动。
  - **可观察差异（OCCT 语义，非回归）**：`Bnd_Box::Get()` 现在带 `Precision::Confusion()` 的两侧 gap —— 每条边/顶点按 `BRep_Tool::Tolerance` 取容差，而 `BRep_Tool.cxx:881-895` 把低于 `Precision::Confusion()` 的容差抬到该值；平面面再加一次 `B.Enlarge(Tolerance(F))`。故单位盒 `Get()` 为 `[-1e-7, 1+1e-7]`，此前采样版给 `[0, 1]`。这正是 `BRepBndLib::Add` 的行为（`Bnd_Box` 的 gap 是 OCCT 本身的设计），旧读数才是假的。
  - 随之按 cxx 盒订正 5 处既有断言（未新增测试、未放宽门禁容差）：`bbox_from_geometry::tests::{edge_bbox_from_curve, face_bbox_from_surface_and_wires}` 改判 `[-t, 1+t]`（`t = Precision::Confusion()`），`draw::tests::bbox_and_verts_commands` 的打印串改 `-0.0000001 / 2.0000001`，`meshing::model_builder::tests::build_model_sets_max_size` 的相对档尺度改 `1 + 2t`，`viz_scene::tests::scene_camera::scene_bounds_box` 改判 `[-t, 1+t]`。
  - **门禁**：`cargo check` 五 crate 0 error；`occt-topo --lib` **1251 passed / 0 failed**；`step_obj_gates` **4/4**（23 模型 f-ratio 1.00，`Shape-2` 0.99；面积比 1.0000，`Shape-2` 0.9999 在其登记容差内，与 `notes/T002-gate-baseline.md` 基线一致）；`export_data_obj` **15/15 ok**，15 个模型 v/f 与包围盒与基线逐位一致（`Cone` 195/301、`Cylinder` 146/140、`Sphere` 642/1244、`Torus` 1369/2592、`Shape` 6150/11372、`linkrods` 3494/5078、`screw` 601/792 全未变）。
  - **旁支（未修）**：`solid_union::shape_bbox`（本文件 §18）仍有独立的 9×9 采样兜底，属另一条路径，未随本项收口。

### 20. 三处自创夹取 Newton 冒充 `Extrema_LocateExtPC`（✅ 2026-10-16 已修，T-164）
- 判定：自创 ｜ live：**原为是**（`int_tools_vertex_line.rs` ← `IntTools_Context::IsVertexOnLine`；`shhealing/wire_fix/same_parameter.rs` ← `BRepLib_ValidateEdge::processApprox`；`brep_lib_continuity.rs` ← `BRepLib::ContinuityOfFaces`）
- 证据（原实现，均非 OCCT 控制流）：
```rust
// int_tools_vertex_line.rs
const LOCATE_ITERS: usize = 12;          // fixed iteration count
if g1.square_norm() <= EXTREMA_EPS * EXTREMA_EPS { break; }   // invented speed gate
if (u1 - u0).abs() <= EXTREMA_EPS { break; }                  // displacement test
// shhealing/wire_fix/same_parameter.rs
for _ in 0..12 { /* Newton */ }          // 1e-20 speed gate + 1e-10 displacement test
// brep_lib_continuity.rs
for _ in 0..100 { /* Newton */ }         // plus a has_analytic_d2() precondition
```
- OCCT 对应：`Extrema_GLocateExtPC::Perform`（`Extrema_GLocateExtPC.hxx:125-307`）不是「固定次数迭代」，也没有位移收敛判据 —— 自由曲线臂在 `NbIntervals(GeomAbs_C2)` 里定位含 `U0` 的区间并 `Extrema_GenLocateExtPC::Initialize` + `Perform`（`hxx:167-176`），未命中则向两侧邻区间走，以区间端点切矢点积 `<= RealEpsilon()`（`hxx:206`/`:243`）判端点极值、否则对该区间再做一次 `Perform`（`hxx:217-219`/`:252-254`）；`Extrema_GenLocateExtPC::Perform`（`Extrema_GenLocateExtPC.hxx:108-127`）本体是 `math_FunctionRoot(myF, theU0, mytolU, myumin, myusup)`，done 后再用 `|ff| >= 1.e-07` 复核并清 `IsDone`（`hxx:116-126`）；初等类型臂跑整段 `Extrema_ExtPC` 并按 `|u - U0|` 取最近解（`hxx:282-303`）。三处消费点的 cxx 出处：`IntTools_Context.cxx:842-871`、`BRepLib_ValidateEdge.cxx:168-227`（每曲线一个持久实例、循环外 `Initialize`）、`BRepLib.cxx:2341-2355`（`aHC2 = new BRepAdaptor_Curve(anEdgeInFace2, theFace2)`，`cxx:2344-2346`）。
- 影响：`brep_lib_continuity.rs` 的旧 `has_analytic_d2` 前置条件使曲面二阶导非解析的边整条静默放弃标记（不是「算错」而是「不算」）；另两处固定迭代次数/位移判据与 cxx 的区间行走语义不同。
- **处置（✅ 已修：2026-10-16，包 P034，commit `c272c08a`）**：新增四个忠实件 `extrema_pc/{curve_tool,gfunc_ext_pc,gen_locate_ext_pc,locate_ext_pc}.rs`，三处调用点改接 `LocateExtPC`，旧 Newton 与其私有助手（`newton_step`/`clamp_param`/`cos_d1`/`cos_d2`/`extrema_f`/`extrema_df`/`locate_on_cos`/`has_analytic_d2`）全部删除；`general_extrema_pc.rs` 增 `extrema_ext_pc_range_done` 以导出初等臂所需的 `IsDone`（`Extrema_GGExtPC.hxx:462-469`：解表可非空而 `mydone == false`）。门禁：五 crate 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15（后两者与 P033 逐字一致）。
- **旁支（未修，已在源码就地标 UNPORTED）**：`CurveOnSurface` 缺 `Adaptor3d_CurveOnSurface::Resolution` 覆写（`Adaptor3d_CurveOnSurface.cxx:1364-1370`：`myCurve->Resolution(min(UResolution, VResolution))`，本仓退回 `Curve` 的 `Precision::Parametric(Confusion)`，`occt-geom/src/curve.rs:135-137`）；`CurveOnSurface` 不暴露 `gp_line`，故 `EvalKPart` 判为 `GeomAbs_Line` 的 cos（`cxx:1565`/`:1687`/`:1701`）在本仓型别为 `GeomAbs_OtherCurve`（两条路都进自由曲线臂，无可观察差异）。覆写落在共享离散适配器（`occt-topo/src/meshing/edge_discret.rs`）上，blast radius 覆盖网格管线，故只登记。
- **诚实边界**：`BRepLib_ValidateEdge` 类本身仍未移植，`shhealing/wire_fix/same_parameter.rs` 只承载 `processApprox`（`cxx:103-229`）的控制流；`brep_check.rs` 的 `BRepCheck_Edge::InContext`（`BRepCheck_Edge.cxx:263-556`）仍 UNPORTED，其依赖是那个未移植的类，与本次 `LocateExtPC` 改动无关（见「补记」末条）。

## 补记（2026-10-06）：T-104 / T-105 / T-106 落地的忠实件

本节补记 §9.588–§9.623 中落在本报告范围内的新移植件。这些都属于"原本缺失、现已按 OCCT 补齐"的控制流，故不改变上文 14 条自创的判定，但修正了若干"未移植"的隐含前提与多处行号漂移。

- **`BRepLib::BuildCurve3d` 族（§9.588）** —— `crates/occt-topo/src/brep_lib_same_range.rs`（新模块）：`rep_ranges:42`、`check_same_range:86`（= `BRepLib::CheckSameRange`）、`same_range:106`（= `BRepLib::SameRange`；文件内 `:35` 登记一处 **UNPORTED (ordering)**：OCCT 走插入有序表，端口未对齐顺序；`:103` 注明 `Tolerance` 是 `GeomLib::SameRange` 的第一实参）。**订正**：旧端口用 `CurveOnSurface` 适配器伪造 3D 曲线，现 `shhealing/shape_build_edge.rs::ShapeBuildEdge::build_curve3d`（`:110`）走真实 `BRepLib::BuildCurve3d`（OCCT `BRepLib.cxx:301-455`；同族 `BRepLib::CheckSameRange` = `:149-183`、`BRepLib::SameRange` = `:187-263`、`evaluateMaxSegment`（`BRepLib.cxx` 文件内静态，`:273-297`）），平面臂用 `GeomLib::To3d`（`occt-geom/src/geom_lib.rs::to_3d:63`，`GeomLib.cxx:559-679`；`geom_lib.rs:403` 注明平面臂上 `GeomLib::To3d` 抛 `Standard_NotImplemented` 的情形），一般臂用 `GeomLib::build_curve3d`（`geom_lib.rs:434`，`GeomLib.cxx:1051-1165`，内部 `is_iso_line:185` / `build_c3d_on_iso_line:258` = `GeomLib::isIsoLine`（`GeomLib.cxx:2991-3077`）/ `GeomLib::buildC3dOnIsoLine`（`GeomLib.cxx:3079-3221`），经 `AdvApprox`）。相配套 `shhealing/shape_fix_edge.rs::{temp_same_range:30, ShapeFixEdge::fix_add_curve3d:84}`（`ShapeFix_Edge.cxx:335-464` / `:618-638`）与 `BRepLib.cxx:149-183`（`CheckSameRange`）/ `:1237`（`SameParameter` 本体）入口。
- **周期锥 loop 修复（§9.590–§9.591）** —— `crates/occt-topo/src/shhealing/shape_fix_face.rs::is_periodic_conical_loop`（`:109`）与 `ShapeFixFace::fix_periodic_degenerated`（`:277`）；接线顺序为 `perform_fix_missing_seam`（`:253`）先 `FixPeriodicDegenerated` 再 `FixMissingSeam`（`:416`）。OCCT：`ShapeFix_Face.cxx:3018-3098`（`IsPeriodicConicalLoop`）、`:3101-3259`（`FixPeriodicDegenerated`）、`:486-489`（先调 `FixPeriodicDegenerated`）/`:492-498`（后调 `FixMissingSeam`）、`:1739`（`FixMissingSeam` 头部 `Context()->Apply(myFace)`）、`ShapeFix_Root.lxx:101`、`ShapeProcess_OperLibrary.cxx:830`。实证（§9.591）：T0M face 1764 由"1 边 0 跨度、`invalid discrete range`"恢复为全 V 段并出网格。
- **`GeomLib::SameRange` 的 2D 臂与 `TransformPCurve`（§9.601 / §9.614–§9.623）** —— `crates/occt-topo/src/shhealing/wire_fix/pcurve_range.rs::geom_lib_same_range`（`:25`，包 P026 后 `shhealing/wire_fix.rs` 已拆为 `wire_fix/` 目录）：等跨 `Geom2d_Line` 平移臂（`GeomLib.cxx:864-871`）、等跨 `Geom2d_Circle` 旋转臂（`cxx:872-889`）、等跨 `Geom2d_TrimmedCurve` 递归臂（`cxx:890-901`）、等跨 `other` 臂（`cxx:908-922`）、不等跨臂（`cxx:924-969`，含 `aCCheck` 周期判据 `:929-932`、`Udeb/Ufin` 夹取 `:948-952`、退化回落 `:954-958`）；`shhealing/shape_build_edge.rs::ShapeBuildEdge::transform_pcurve`（`:280`，`ShapeBuild_Edge.cxx:596-700`，圆锥臂 + 兜底臂 + Bezier 臂 + BSpline 臂）；`crates/occt-topo/src/geom_bnd_lib_bspline2d.rs`（`GeomBndLib_BSplineCurve2d.cxx:34-57`、`GeomBndLib_Curve2d.cxx:137-157` / `:150-157`、`GeomBndLib_BSplineSurface.cxx:30` / `:103` / `:118`；文件头 `:3-7` 登记"忠实但 OCCT 侧消费方尚未全部接线"）；`crates/occt-topo/src/pcurve_full/surface_projector/reparam.rs` 的 `ReparamCurve2d`（`:23`）/ `same_range_window`（`:160`）/ `trimmed_reparam_bspline`（`:201`）/ `same_range_bspline2d`（`:226`）/ `reparam_curve2d`（`:260`）。**2026-10-16 订正（包 P026）**：两条转换尾原本各有一份副本 —— `reparam.rs` 与 `pcurve_range.rs` 各写一套「`Geom2d_TrimmedCurve` → `Geom2dConvert::CurveToBSplineCurve` → `BSplCLib::Reparametrize` → `SetKnots`」（`GeomLib.cxx:912-921` 与 `:962-968`），且 `pcurve_range.rs` 的副本在 cxx 不做判型的不等跨臂上按曲线类型分派（`line_trim_to_reparam_bspline` / `unequal_trim_window` → `same_range_unequal_line` / `same_range_unequal_other`）。现两份合并为唯一忠实件 `trimmed_reparam_bspline`（`cxx:912-921` + `:962-968`：`Geom2dTrimmedCurve(CurvePtr, U1, U2)` → `Geom2dConvert::CurveToBSplineCurve` → 对 `BS->Knots()` 做 `BSplCLib::Reparametrize(RequestedFirst, RequestedLast, Knots)` → `BS->SetKnots(Knots)`）与 `same_range_window`（`cxx:934-959`），上述分派件与 `ReparamCurve2d::bspline_copy2d` 整体删除（wrapper 只在 `Geom2dConvert::CurveToBSplineCurve` 拒绝时存在，`Geom2dConvert.cxx:442-444`，故 trait 默认 `None` 才是正确答案）。`occt-geom2d` 侧的 `CompCurveToBSplineCurve` / `ApproxCurve` / `Geom2d_BSplineCurve` 操作族见 `occt-geom.md` 的 2026-10-06 补记；网格侧 `add_wire` 整链反转删除与 `pcurve_and_range` 见 `mesh-exchange.md` 发现 5c。
- **订正：上文条目的行号与前提（2026-10-06 逐条 grep 复核）**
  - 发现 1/2/3/4/5（`brepfeat/features.rs`）：文件仍在，物理行号已漂移 —— `result.volume = ...` 现 `:109`（未变）；`revolve_profile_about` 现 `:336`（原文写 `:361`）；`mesh_cylinder(radius, height, 24)` 现 `:394`（原文写 `:419`）；`r.clamp(16, 64)` 现 `:483`（原文写 `:508`）；`GeometryRegistry::global().set_face` 现 `:37`（未变）。
  - 发现 6（`wireframe.rs:399` / `:407-408`）：现为 `:423`（`face_uv_bounds(f, surface.as_ref())`）与 `:431-432`（`clamp(3, 64)`）；`face_uv_bounds` 函数本体现 `:486`。**仍未修**（被 A13/A18/T-69 阻塞）。
  - 发现 7（`feature.rs` 顶层重复体素实现）：**该文件已不存在**（`crates/occt-topo/src/feature.rs` 已删除），仅 `brepfeat/features.rs` 保留 —— 本条据此订正为"重复体素实现已收敛到 `brepfeat`"。
  - 发现 8/9/10/11 的行号经复核仍成立：`validate.rs:150` 的 `is_valid()` 恒真仍在（`:3` 已加 **UNPORTED**）；`brep_gprop.rs:98`、`shape_metrics.rs:146`、`solid_union.rs:103` 均在。**订正（2026-10-14）**：发现 11 与 13 均已修（发现 11 见本条处置，发现 13 的 `brep_measure.rs:12` 折线采样已换成 `BRepGProp::LinearProperties` 委托、形参 `samples` 删除），故第 13 条不再是未整改项。
  - `补充说明` 的 `brep_check.rs edge_check`：`BRepCheck_Edge::Minimum` **已移植**（2026-10-13，`cxx:74-258`）——`!SameRange && SameParameter` → `InvalidSameParameterFlag`、无 3D 表示 → `No3DCurve`、`Last <= First` 与基曲线周期/范围比对 → `InvalidRange`、空曲线表示 + `Degenerated` → `InvalidDegeneratedFlag`，`Status` 按 `BRepCheck_Status.hxx:20-59` 重排补齐；`InContext`（`:263-556`，依赖未移植的 `BRepLib_ValidateEdge`）、`Tolerance`（`:598-707`）、`CheckPolygonOnTriangulation`（`:711-839`）与表示的 `TopLoc_Location`（`:163-171`）仍 UNPORTED，已在函数头登记。`BRepLib::CheckSameRange` 已随 §9.588 落地（见上），两者不是同一函数，勿混同。

## 补充说明
- `brep_offset/shell_offset.rs:130` 对非凸/曲面直接 `Err`（保守拒绝而非自创几何），但未登记"未移植 `BRepOffset_MakeOffset` 迭代求交"。
- `shape_checks.rs:128 orientation_consistency_check` 以面心-壳质心判外向（自承 heuristic），被 `brep_check.rs:106` 复用为 `NotClosed` 判据；OCCT `BRepCheck_Shell` 用朝向一致性而非质心。
- `fillet.rs:1-8`、`fillet_edge/chain.rs:29-37`（球八分之一角片代替 `ChFi3d` 滚动球角）、`fillet_curved/`、`fillet_var/mod.rs:1-14` 均自承 "simplified port" 且无门禁覆盖（`fillet_curved` 仅 `cone_dbg`/`phase10_integration.rs:31`），建议按组登记而非逐条计自创。
- `brep_check.rs edge_check` 已是 `BRepCheck_Edge::Minimum`（`cxx:74-258`）的移植（2026-10-13）；`InContext`/`Tolerance`/`CheckPolygonOnTriangulation` 仍 UNPORTED（见上）。
