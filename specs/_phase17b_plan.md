# Phase 17b Plan — 波 C1 补足：BeanFaceIntersector + EdgeFace（+4,000 行，总量达 5,000）

> Phase 17 交付 HInter（1,660）+ FClass2d（813）= 2,473 行，不达 5,000。补足：波 C 主线数值核 **BeanFaceIntersector**（EdgeFace 的数值心脏）+ **EdgeFace**（完整求交编排）。
> 关键：BeanFaceIntersector 全部依赖已就绪——ExtCS/GenExtCS（Phase 13 `extrema_surf.rs`）、IntAna（`intana.rs`）、HInter（`intcurvesurface.rs`）、FClass2d（`fclass2d.rs`）、Context/CurveRangeSample/SurfaceRangeSample（波 A）。
> 参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKBO\IntTools\IntTools_BeanFaceIntersector.cxx`（2,639）+ `IntTools_EdgeFace.cxx`（875）。
> 基线：5 crate · 138,072 行 · 1,785+27=1,812 测试（topo 1112）。

## 方案（各 agent 独立文件，无冲突）

| # | 缺口 | 方案 | 模块（新建） | 预估 |
|---|---|---|---|---|
| 1 | 边-面求交数值核（无） | `IntTools_BeanFaceIntersector`：Perform → 边落在面内的参数区间集（`Result`）。路径：FastComputeAnalytic（解析对快速判定/重合）→ ComputeLinePlane（直线×平面）→ ComputeAroundExactIntersection（HInter 精确交点附近区间）→ ComputeUsingExtremum（ExtCS 极值距离判定）→ ComputeLocalized（CurveRange/SurfaceRange 局部化）。用 HInter 交点 + ExtCS 距离 + 符号判定确定边穿越面的区间 | `occt-topo/src/bean_face.rs` | 2,500 |
| 2 | 边-面完整求交（无） | `IntTools_EdgeFace`：Perform → CommonPrt 序列（边落在面上的公共区间）。CheckData/IsProjectable/DistanceFunction/MakeType/CheckTouch/CheckTouchVertex/IsCoincident。核心：BeanFaceIntersector 结果区间 → 每区间 MakeType（点/边公共部分类型）→ 分类（FClass2d 或 3D 判定） | `occt-topo/src/edge_face.rs` | 1,500 |

**合计 ≈ 4,000 行** → 阶段总量 ≈ 6,500。

## OCCT 类 → Rust 映射要点

- **BeanFaceIntersector**：`new(curve, surface, tol_e, tol_f)`、`set_bean_parameters(first,last)`、`set_surface_parameters(umin,umax,vmin,vmax)`、`set_context(ctx)`、`perform() -> Result<(),String>`、`is_done()`、`result() -> &[IntRange]`、`minimal_square_distance()`。
  - `FastComputeAnalytic`：解析对（直线×平面/球/柱/锥/环）用 intana/解析式判定是否重合或可快速求交；成功则直接给区间。
  - `ComputeAroundExactIntersection`：HInter 交点处向两侧扩展参数区间直到离开面（d0 距离 > tol + 符号变化）。
  - `ComputeUsingExtremum`：用 `crate::extrema_surf`（Phase 13 的 ExtCS，先查函数名）算曲线-曲面极值距离，若全距 > tol → 无交。
  - `ComputeLocalized`：`CurveRangeSample`/`SurfaceRangeSample` 细分参数域局部化。
  - 输出 `Vec<IntRange>`（边落在面上的参数区间）。
- **EdgeFace**：`new(edge, face)`、`perform()`、`common_parts() -> &[CommonPrt]`、`is_done()`。
  - `CommonPrt` 用 `crate::inttools_data::CommonPrt`（Phase 16 已移植）。
  - `MakeType`：区间类型（点=IsRoot / 边=Edge / 未知）。`CheckTouch`/`CheckTouchVertex`：相切判定（用 FClass2d 或距离符号）。
  - 分类：区间中点 3D 点 → `fclass2d::FClass2d::perform`（UV 判定）或距离符号。
  - `IsCoincident`：边曲线完全在面上（每采样点距离 < tol）。
- **接线**：不接 bop_builder（波 C2 统一接）；仅模块内部自洽。

## 验证（每 agent 自测）

- **BeanFaceIntersector**：直线穿过盒面 → 1 个区间（交段）；直线在面外 → 0；直线贴面（距离 < tol）→ 重合区间；BSpline 曲线 vs 球面 → 穿越区间正确。
- **EdgeFace**：盒面 + 穿过边 → CommonPrt 含 1 个 Edge 型区间；边在面外 → 空；边端点贴面 → IsRoot 型；与 `crate::inttools::edge_face_intersections` 结果一致（现有近似版对拍，参数接近）。
- 用 `crate::brep_extrema::test_box::unit_box`（先确认路径）的 face + 构造边。

## 依赖与纪律

- 每个 Agent 只写自己的新模块文件；`lib.rs` 由 orchestrator 注册。
- 匹配现有风格（`Result<T,String>`、`mod tests`、`Arc<dyn Curve/Surface>`、GpPnt、`crate::inttools_data::IntRange/CommonPrt`）。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target`。
- 不碰 bop_builder/bop_curved/boolean_ops 等既有模块。
