# Phase 18 Plan — 精确 NURBS 布尔 · 波 C2a：EdgeEdge 一般曲线 + FaceFace 精确 + AlgoTools（≥5,000 行）

> 精确布尔收官第二波第一子波。侦察结论：整个 C2（EdgeEdge 1,659 + FaceFace 3,111 + BOPTools_AlgoTools 3,300 + PaveFiller 10,668 C++）≈ 13,000-15,000 Rust 行，需拆 C2a（本阶段）+ C2b（PaveFiller 独立阶段）。
> **C2a**：EdgeEdge 一般曲线求交 + FaceFace 精确求交（薄编排器）+ BOPTools_AlgoTools 族（PaveFiller 构建依赖）。
> 参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKBO\IntTools\{EdgeEdge,FaceFace}.cxx` + `BOPTools\{AlgoTools,AlgoTools2D,Set,AlgoTools3D}.cxx`。
> 基线：5 crate · 143,770 行 · 1,830 测试（topo 1130）。纯增量新模块。

## 缺口 → 方案（各 agent 独立文件，无冲突）

| # | 缺口 | 方案 | 模块（新建） | 预估 |
|---|---|---|---|---|
| 1 | 边-边一般曲线求交（现仅线/圆解析） | `IntTools_EdgeEdge`：Prepare（类型+容差）→ Perform → FindSolutions（extrema_cc 极值附近区间 + 符号二分）→ FindParameters（交点参数对）→ MergeSolutions。线线/圆圆解析保留（复用现有），一般曲线用 `crate::extrema_cc::curve_curve_extrema_all` 找极值 + 区间收缩 | `occt-topo/src/edge_edge.rs` | 1,800 |
| 2 | 面-面精确求交（现为网格采样） | `IntTools_FaceFace`：Perform → 交线集合。**薄编排器**（跳过 IntPatch 行走器）：解析对（平面×平面/球/柱/锥/环、球×球等）用 `crate::intpatch` 解析（intersect_plane_sphere/intersect_sphere_sphere 等）→ `IntersectionCurve`；一般曲面回退 intpatch 采样/网格路径。输出交线 + 每面 pcurve | `occt-topo/src/face_face.rs` | 1,500 |
| 3 | PaveFiller 构建依赖（无） | `BOPTools_AlgoTools`：MakeNewVertex（交点→新顶点）、MakeEdge（交线→边）、MakePCurve（边→面 pcurve，复用 `pcurve_full`）、ComputeVV（顶点合并）、PointOnEdge/UpdateVertex、GetNormalToSurface（法向，d1 或数值差分）+ AlgoTools2D 的 EdgeToFace/AdjustPCurve + Set 工具 | `occt-topo/src/algo_tools.rs` | 2,200 |

**合计 ≈ 5,500 行**。

## OCCT 类 → Rust 映射要点

- **IntTools_EdgeEdge**：`new()`、`set_edge1/edge2`、`set_range1/range2`、`set_fuzzy_value`、`perform() -> Result<(),String>`、`is_done()`、`common_parts() -> &[CommonPrt]`、`points() -> &[PntOn2Faces]`。
  - 一般曲线：`curve_curve_extrema_all` → 距离 < tol 的极值对 → 每对参数做区间二分（距离符号）→ `FindParameters`（用 `geom2d_api::project_point_on_curve` 或 extrema 精确参）。
  - `ComputeLineLine`/`ComputeCircleCircle`（解析）保留——复用现有 `inttools::edge_edge_intersections` 或重写解析式。
  - `MergeSolutions`：参数去重排序（复用 `inttools_roots::remove_identical_roots`）。
- **IntTools_FaceFace**：`new()`、`set_face1/face2`、`set_tolerance`、`perform() -> Result<(),String>`、`is_done()`、`nb_curves()/curve(i) -> IntersectionCurve`、`nb_paves(i)/pave(i)`。
  - `IntersectionCurve { kind: CurveKind, curve: Arc<dyn Curve>, range: IntRange, face1_idx, face2_idx, pcurve1: Option<Arc<dyn Curve2d>>, pcurve2 }`。
  - 解析对：面分类（`brep_surface::classify_surface`）→ `intpatch` 解析函数（先查其 `IntersectionCurve` 类型与全部解析对函数）→ 一般曲面回退 `intpatch::surface_surface_intersection`（现有采样，先查签名）。
  - 输出交线 + 每面 pcurve（`pcurve_full::make_pcurve_full`）。
- **BOPTools_AlgoTools**：
  - `MakeNewVertex(p, tol) -> Vertex`、`MakeEdge(c, v1, t1, v2, t2, tol) -> Edge`（用 `crate::builder::TopoBuilder`，先查 API）、`MakePCurve(edge, face) -> Arc<dyn Curve2d>`（委托 `pcurve_full`）。
  - `ComputeVV(v, p, tol) -> i32`（顶点重合判定）、`PointOnEdge(edge, t) -> GpPnt`、`UpdateVertex(ic, t, v)`、`GetNormalToSurface(s, u, v) -> GpVec`（d1 或数值差分，`brep_surface::surface_normal` 已有）。

## 验证（每 agent 自测）

- **EdgeEdge**：两直线相交 → 1 交点参数对；两圆相交 → 2 交点；直线×BSpline → 极值+交点；分离曲线 → 空；与 `inttools::edge_edge_intersections` 对拍一致。用 `crate::brep_extrema::test_box::unit_box`（先确认路径）的边。
- **FaceFace**：两平面相交 → 直线交线（解析精确）；平面×球 → 圆交线；球×球 → 圆交线；一般曲面×球 → 回退路径交线合理；与 `intpatch::surface_surface_intersection` 对拍。
- **AlgoTools**：MakeNewVertex 在交点处建顶点、MakeEdge 建边（顶点连接）、MakePCurve 与 `pcurve_full` 一致、ComputeVV 重合/分离判定、GetNormalToSurface 平面法向正确。

## 依赖与纪律

- 每个 Agent 只写自己的新模块文件；`lib.rs` 由 orchestrator 注册。
- 匹配现有风格（`Result<T,String>`、`mod tests`、`Arc<dyn Curve>`、`dyn Surface`、GpPnt、`crate::inttools_data::IntRange/CommonPrt`）。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target` 避免并行文件锁与磁盘满。
- 不碰 bop_builder/bop_curved/boolean_ops 等既有模块。
