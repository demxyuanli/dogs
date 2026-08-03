# Phase 16 Plan — 精确 NURBS 布尔 · 波 A：IntTools 机械层 + pcurve 基建（≥5,000 行）

> 精确布尔数值内核第一阶段。侦察结论：原「波 A(IntTools 求交)+波 B(分类)」不可行——IntTools FaceFace/EdgeFace 会拖进 BeanFaceIntersector(2.6k)+IntCurveSurface_HInter(2.1k)+IntPatch(35.9k)，单波 15k+ 全高风险。
> **重切分后的波 A**：只做 IntTools 机械层（数据类/Tools/ShrunkRange/Context/CommonPrt）+ pcurve 基建（tgeometry 字段 + MakePCurveOnFace）。纯移植、低风险、不碰现有布尔 pipeline。
> 参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKBO\IntTools\`（24 文件 15,109 行）+ `BOPTools\AlgoTools2D`（pcurve）。
> 基线：5 crate · 132,905 行 · 1,725 测试（topo 997）。纯增量新模块 + 1 处小改动（tgeometry EdgeGeom）。

## 缺口 → 方案（各 agent 独立文件，无冲突）

| # | 缺口 | 方案 | 模块（新建） | 预估 |
|---|---|---|---|---|
| 1 | IntTools 数据类（无） | `CommonPrt`（公共部分）+ `Root`/`Range` + `PntOnFace`/`PntOn2Faces` + `Curve` + `MarkedRangeSet` + `LocalizeData`×2 + `SequenceOf*` 别名 | `occt-topo/src/inttools_data.rs` | 1,500 |
| 2 | IntTools 工具/收缩/上下文（无） | `IntTools_ShrunkRange`（GCPnts_AbscissaPoint + 边收缩）+ `IntTools_Tools`（静态工具）+ `IntTools_Context`（FClass2d/SolidClassifier/Hatcher 缓存壳） | `occt-topo/src/inttools_range.rs` | 1,800 |
| 3 | pcurve 基建（无存储） | `MakePCurveOnFace`（ge2d 曲线构造：line/circle/trimmed + 一般曲线采样拟）+ tgeometry `EdgeGeom` 加 `pcurves: HashMap<face_key, Arc<dyn Curve2d>>` 字段 + registry 存取器 | `occt-topo/src/pcurve.rs` + 小改 `tgeometry.rs` | 1,800 |

**合计 ≈ 5,100 行**，达 5,000 目标。

## OCCT 类 → Rust 映射要点

- **数据类**：`IntTools_CommonPrt`（type/range/face/vertices）；`IntTools_Root`（iRoot/typeRoot/range/conflict）与 `IntTools_Range`（first/last）；`PntOnFace`（face+UV+3D）、`PntOn2Faces`（两面上的点+连接）；`IntTools_Curve`（curve type/surface/range）；`MarkedRangeSet`（区间并/交集合，排序）。
- **ShrunkRange**：`SetShrunkRange(edge, face, tol)`——用 `GCPnts_AbscissaPoint`（`crate::gcpnts` Phase 13 已有）沿边做弧长收缩，输出内缩 [first,last]；`IsDone`。
- **Tools**：静态函数（`ComputeTolerance`/`IsInRange`/`IsVertex`/`MiddlePoint` 等）。
- **Context**：`SetEdge(edge, face)` 缓存 `BRepAdaptor_Curve`/`Surface` + `FClass2d`（波 B 才填，本波留壳）；`IsPointOnFace`/`ProjectPointOnFace` 用 `brep_surface::surface_closest_params`。
- **MakePCurveOnFace**：edge→face 的 UV 曲线。线/圆/平面直接用 `occt-geom2d` 的 `Geom2dLine/Circle/TrimmedCurve` 构造；一般曲线沿面采样（`brep_surface::edge_pcurve_on_face` 已有采样逻辑）再拟 `bspline_curve`。返回 `Arc<dyn Curve2d>`。
- **tgeometry 改动**（唯一既有文件改动，~90 行）：`EdgeGeom` 加 `pcurves: HashMap<usize, Arc<dyn Curve2d>>`（face key→curve）+ `set_pcurve/get_pcurve`；`GeometryRegistry` 加 `edge_pcurve` 存取。向后兼容（新字段默认空 map），不破坏现有 997 topo 测试。

## 验证（每 agent 自测）

- **数据类**：CommonPrt 构造/访问、Range 包含/排序、MarkedRangeSet 并/交/差正确。
- **ShrunkRange**：单位盒一条边长 1 → 收缩后 first>0 且 last<1（弧长 >0）；退化边 Err。
- **pcurve**：盒面直线边 → pcurve 是直线（端点 UV 与面投影一致）；圆边 → 圆弧 pcurve 端点正确；一般曲线面 → 采样 pcurve 两端点命中面边界。

## 依赖与纪律

- Agent 1/2 各写独立文件；Agent 3 写 pcurve.rs + **只改 tgeometry.rs 一处**（加字段/存取器，不碰其它）。
- 匹配现有风格（`Result<T,String>`、`mod tests`、`Arc<dyn Curve/Surface>`、`GpPnt2d`、`occt_geom2d` 的 Curve2d trait）。ShrunkRange 复用 `occt_core::gcpnts`。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target` 避免并行文件锁。
- 不碰 bop_builder/bop_curved/boolean_ops/bopds 等既有模块。
