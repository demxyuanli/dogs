# occt-topo（BOP/网格/交换以外）移植忠实度审查（2026-09-20）

范围：`crates/occt-topo/src/`。只读，未改源码。
统计：自创 14 条、等价替换未登记 4 条（共 18 条，按严重度排序）。

**结论摘要**：最高危三条 —— ①`brepfeat/features.rs:109` 用解析公式**覆盖布尔结果的体积**；②`wireframe.rs:399+407` 对曲面按**未裁剪自然 UV 域**建网格，`shape_volume`/`shape_surface_area`/`brep_gprop` 全链失真且被 `bop_builder_report` 当修复判据；③`feature.rs:100`+`brepfeat/features.rs:318`+`boolean_ops.rs:173` 体素布尔在精确 `bop_builder::boolean`（BOPAlgo_BOP）已存在时仍是 feature 默认路径，注释仍称"exact BRep boolean 不可用"（过期）。另：`hlr.rs`/`viz_scene/`/`draw/`/`xcaf/`/`render_svg.rs` 借用了 OCCT 包名，但不是移植。所有条目均未进入 `export_data_obj`/`step_obj_parity` 主门禁，除第 6、11、13 条（经 `bop_builder_report`/`shape_mesh` 间接进入）。

## 发现

### 1. `boss_thru_all` 用解析公式覆盖布尔结果体积
- 判定：自创 ｜ live：**是**（公有 `brepfeat::boss_thru_all`，`brepfeat/tests.rs:171`）
- 证据：`crates/occt-topo/src/brepfeat/features.rs:109`
```rust
result.volume = (v0 + PI * radius * radius * height - overlap).max(0.0);
```
- OCCT：`TKFeat/BRepFeat/BRepFeat_MakePrism.cxx`（`BRepAlgoAPI_Fuse`，体积事后由 `BRepGProp` 积分）；`BRepFeat_MakeRevol.cxx:310-315` 无此步。
- 影响/建议：返回体积与 `result.shape` 实际体积不一致 → `volume_delta`/`ratio`/`kind`（同文件 `:161-179`）全失真，`tests.rs:172` 因覆盖恒真。删覆盖，改用 `brep_gprop_full`；否则标未移植。

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

### 7. `feature.rs` 是 `brepfeat` 的重复体素实现 + 过期注释
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

### 10. `validate.rs` 不是 `BRepCheck_Analyzer` 的移植
- 判定：自创 ｜ live：**是**（`pub mod validate`，`phase3/5` 集成测试使用；`brep_check.rs` 才是 `BRepCheck` 的移植）
- 证据：`crates/occt-topo/src/validate.rs:69`
```rust
r.messages.push(format!("shell has only {faces} faces (a closed shell needs ≥4)"));
pub fn is_valid(&self) -> bool { true }   // validate.rs:131
```
- OCCT：`BRepCheck_Shell::CheckClosed` 由边沿面使用次数判定闭壳，无"面数≥4"规则；`IsValid` 非常量 true。
- 影响/建议：`Analyzer::is_valid()` 恒真使相关门禁失效；4 面阈值把合法小面片壳判 `Undefined`。删除启发式，转调 `brep_check.rs`；保留 API 标未移植。

### 11. `brep_gprop.rs` 是网格采样，不是 `BRepGProp`
- 判定：自创（模块头自承 "Approximate port"）｜ live：**是**（`bop_builder_report.rs:405,410` 修复前后面积判据、`brep_pipe.rs:77`、`loft.rs:386`、`gprop_analytic/analytic_props.rs:418+`）
- 证据：`crates/occt-topo/src/brep_gprop.rs:98`
```rust
let n = ((1.0 / deflection.max(1e-4)).ceil() as usize).clamp(1, 32);
```
- OCCT：`BRepGProp::SurfaceProperties`/`VolumeProperties` → `BRepGProp_Gauss`。
- 影响/建议：`1/deflection` 夹 `[1,32]` 是裸规则，面积/体积/质心随 deflection 漂移，而 OCCT 有 `myEpsilon` 误差估计。转发 `brep_gprop_full`。

### 12. `offset_by_sampling` 把一般曲面采成 16×16 双线性面
- 判定：自创 ｜ live：**是**（`offset_face_surface:644` 的 cone/torus/其他分支 → `offset_face:650` → `offset_shell`，即 `draw` 的 offset 命令）
- 证据：`crates/occt-topo/src/brep_offset/curve_face_offset.rs:588`
```rust
let (nu, nv) = (16usize, 16usize);
```
- OCCT：`BRepOffset_Offset`/`BRepOffsetAPI_MakeOffsetShape` 构造精确 `Geom_OffsetSurface`。
- 影响/建议：cone/torus 偏移面退化为 16×16 分段双线性面；`brep_offset/curve_face_offset.rs:648-649` 还自承非平面面"边不投影到偏移面"。改用 `Geom_OffsetSurface`，否则标未移植，不返回近似面。

### 13. `edge_length` 是等距折线采样
- 判定：自创 ｜ live：**是**（`bop_builder_report.rs:253`、`bop_builder_repair.rs:634`、`algo_tools/construct.rs:546`、`shhealing/common.rs:492` 作短边判据；`shape_metrics.rs:40`）
- 证据：`crates/occt-topo/src/brep_measure.rs:12`
```rust
pub fn edge_length(e: &Edge, samples: usize) -> f64 {
```
- OCCT：`BRepGProp::LinearProperties`（`BRepGProp.cxx:121`）→ `BRepGProp_EdgeTool` 积分，带误差估计。
- 影响/建议：采样数由调用方给（8/32），圆弧短边在 8 段折线下系统性低估 → 短边判据结果被替换规则左右。改调 `BRepGProp` 线性积分（`brep_gprop_full`）。

### 14. 用采样识别二次曲面代替 downcast
- 判定：等价替换未登记（后果改变）｜ live：**是**（`offset_face_surface:613,619,632` 的分派依据，连同 `brep_surface::is_planar` 采样版）
- 证据：`crates/occt-topo/src/brep_offset/curve_face_offset.rs:551`
```rust
pub(super) fn is_cylindrical(s: &dyn Surface, loc: &GpPnt, dir: &GpDir, r: f64) -> bool {
```
- OCCT：`BRepOffset_Offset` 用 `Geom_Cylinder`/`IsKind(STANDARD_TYPE(...))` 类型判据，无采样。
- 影响/建议：小曲率 B 样条等可被误判为圆柱 → 走错偏移分支。`tgeometry.rs` 已声明 `dyn Surface` 无法 downcast，但本模块未登记。在模块头登记该替换与误差条件，或给 `occt_geom::Surface` 加类型标签。

### 15. `hlr.rs` 是画家算法近似，不是 `HLRBRep`
- 判定：自创 ｜ live：**是**（`render_svg.rs:57`、`draw` 的 view 命令）
- 证据：`crates/occt-topo/src/hlr.rs:126`
```rust
/// ponytail: conservative segment test — a partial occlusion (an edge crosses
```
- OCCT：`TKHLR/HLRBRep/**`（`HLRBRep_Algo`/`Hider`/`EdgeBuilder`）基于精确边面求交与轮廓线。
- 影响/建议：凹形、自遮挡、轮廓线的可见边集合与 OCCT 不同。标未移植（模块名与 OCCT 包同名，易误认为已移植，应在 `lib.rs:59` 注明）。

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

### 18. `solid_union::shape_bbox` 9×9 采样兜底 + `brep_class3d` 确定性探针
- 判定：等价替换未登记 ｜ live：**部分**（`solid_union` 本体仅 `bop_builder_tests.rs:174` 对拍；`bbox_union` 被 `boolean_ops.rs:184` 用；`brep_class3d.rs:96 find_a_point_in_the_face` 供分类器使用）
- 证据：`crates/occt-topo/src/solid_union.rs:103`
```rust
for i in 0..=8 {
```
- OCCT：`BRepBndLib::Add`/`BndLib` 解析极值；`BRepClass3d_SolidExplorer::FindAPointInTheFace`（`BRepClass3d_SolidExplorer.cxx:512-560`）逐面取参数点并带 `aTestInvert` 反转重试。
- 影响/建议：9×9 格点会漏极值 → 体素网格范围偏小；探针方向固定，不等价于 OCCT 多方向重试。包围盒改 `geom_bnd_lib_*` 解析路径；`OtherSegment` 补 `aTestInvert` 分支并登记简化。

## 补充说明
- `brep_offset/shell_offset.rs:130` 对非凸/曲面直接 `Err`（保守拒绝而非自创几何），但未登记"未移植 `BRepOffset_MakeOffset` 迭代求交"。
- `shape_checks.rs:128 orientation_consistency_check` 以面心-壳质心判外向（自承 heuristic），被 `brep_check.rs:106` 复用为 `NotClosed` 判据；OCCT `BRepCheck_Shell` 用朝向一致性而非质心。
- `fillet.rs:1-8`、`fillet_edge/chain.rs:29-37`（球八分之一角片代替 `ChFi3d` 滚动球角）、`fillet_curved/`、`fillet_var/mod.rs:1-14` 均自承 "simplified port" 且无门禁覆盖（`fillet_curved` 仅 `cone_dbg`/`phase10_integration.rs:31`），建议按组登记而非逐条计自创。
- `brep_check.rs:56-72 edge_check` 只查曲线存在与范围有限，`BRepCheck_Edge` 的 SameParameter/SameRange/重合边检查未移植。
