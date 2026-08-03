# Phase 17 Plan — 精确 NURBS 布尔 · 波 C1：曲线-曲面求交 + 2D 面分类（≥5,000 行）

> 精确布尔数值内核核心波 C 的第一子波。侦察结论：波 C 全量（BeanFaceIntersector 2,639 + EdgeEdge 1,659 + FaceFace 3,111 + IntCurveSurface 2,132）≈ 7,000-9,000 Rust 行高风险；切 C1（本阶段，独立可验证）+ C2（后续）。
> **C1 内容**：`IntCurveSurface_HInter`（曲线-曲面精确求交，解析+一般）+ `IntTools_FClass2d`（面 2D 分类，EdgeFace 前置依赖）。
> 参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKGeomAlgo\IntCurveSurface\`（HInter 581 + InterImpl.lxx + 辅助类，2,132 行）+ `TKBO\IntTools\IntTools_FClass2d.hxx/.cxx`（950 行）。
> 基线：5 crate · 138,072 行 · 1,785 测试（topo 1085）。纯增量新模块 + 可选小改。

## 缺口 → 方案（各 agent 独立文件，无冲突）

| # | 缺口 | 方案 | 模块（新建） | 预估 |
|---|---|---|---|---|
| 1 | 曲线-曲面精确求交（无） | `IntCurveSurface_HInter`：`Perform(曲线, 曲面)` → 交点集（`IntersectionPoint`：UV+参数+3D+状态）+ 段集。解析路径（直线/圆/椭圆/抛物/双曲 × 平面/球/柱/锥/环，`PerformConicSurf` 用现有 intana/解析式）；一般路径（InterImpl 模板：多面体细分 + BoundSortBox 裁剪 + 区间精确/牛顿） | `occt-topo/src/intcurvesurface.rs` | 3,200 |
| 2 | 面 2D 分类（无） | `IntTools_FClass2d`：`Perform(Puv) -> State{In,Out,On,Unknown}` + `PerformInfinitePoint` + `TestOnRestriction` + `IsHole`。用面 UV 边界（`brep_surface::edge_pcurve_on_face`/`pcurve_full`）构造 2D 区域，`occt_core::geom::polygon_ops::point_in_polygon2d`（凸）+ 环分解做 in/out；弯曲边用 pcurve 采样多边形。BRepClass_FaceExplorer 简化（面的 UV 环→多边形集合） | `occt-topo/src/fclass2d.rs` | 1,800 |

**合计 ≈ 5,000 行**，达目标。

## OCCT 类 → Rust 映射要点

- **IntCurveSurface_HInter**：
  - `IntersectionPoint { param: f64, u: f64, v: f64, pnt: GpPnt, state: TopAbsState }`；`IntersectionSegment { first: IntersectionPoint, last: IntersectionPoint }`；`State { In, Out, On, Unknown }`（映射 `TopAbs_State`）。
  - `Perform(curve: &dyn Curve, surface: &dyn Surface, u/v 域, curve 域)` → `Result<HInterResult, String>`（`nb_points/point(i)/nb_segments/segment(i)`）。
  - **解析路径**：`PerformConicSurf`（直线/圆/椭圆/抛物/双曲 × 解析曲面）——对平面：直线→代入平面方程求根；对球/柱/锥/环：参数化代入 → 一元方程（用 `occt_math` roots）。复用现有 `occt_geom` 的解析曲面类型（GeomPlane/GeomCylinder/GeomSphere/...）。
  - **一般路径**（InterImpl 简化）：曲线采样成折线（`occt_core::gcpnts`）→ 每段与曲面包围盒（`BndBox`）粗筛 → 细分区间，每区间做 d0 距离 + 符号判断 → 区间收缩/牛顿精确根（`occt_math::newton`）。目标鲁棒（平面/球/柱精确 <1e-9，一般曲面 <1e-6）。
- **FClass2d**：
  - `Init(face, tol)`：从面取 UV 边界环（外环 + 孔环），各边用 pcurve（`pcurve_full::make_pcurve_full` 或 `brep_surface::edge_pcurve_on_face` 采样）转 2D 多边形。
  - `Perform(Puv)`：2D 点在面区域判定——外环包含且不在孔内 → In；在边界上 → On；否则 Out。弯曲面用多边形近似（pcurve 采样）。
  - `PerformInfinitePoint`：UV 域外左下角点 → 期望 Out（开放面）或 In（闭合周期面）。
  - 周期性处理：U/V 周期面 `AdjustPeriodic` 把 UV 折叠进 [Umin,Umax]×[Vmin,Vmax]。
- **接线**：不接 bop_builder（波 C2 统一接）；但 `inttools.rs` 的 `edge_face_intersections` 可选加 FClass2d 增强——本阶段不做，避免破坏现有路径。

## 验证（每 agent 自测）

- **HInter**：直线穿过球面 → 2 个交点（解析精确，距离 <1e-9）；直线与平面平行 → 0 交点；圆边 vs 平面 → 2 交点或相切 1；一般 BSpline 曲线 vs 球面 → 交点与网格求交一致（<1e-6）。`occt_geom` 有 GeomSphere 等构造（先查 `crates/occt-geom/src/lib.rs`）。
- **FClass2d**：单位盒面（UV 正方形）中心 → In、角外 → Out、边界中点 → On；带孔面 → 孔内 Out；圆柱侧面（周期 U）→ 周期折叠后正确。用 `crate::brep_extrema::test_box::unit_box`（先确认路径）的 face + `pcurve_full` 采样。

## 依赖与纪律

- 每个 Agent 只写自己的新模块文件；`lib.rs` 由 orchestrator 注册。
- 匹配现有风格（`Result<T,String>`、`mod tests`、`Arc<dyn Curve/Surface>`、GpPnt/GpPnt2d）。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target`。
- 不碰 bop_builder/bop_curved/boolean_ops 等既有模块；不接现有 inttools（防回归）。
