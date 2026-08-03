# Phase 18b Plan — 波 C2a 补足：AlgoTools 完整 + IntTools_Context 完整 + AlgoTools2D（+4,700 行）

> Phase 18 交付 edge_edge(841) + int_face_face(1004) + algo_tools(385) = 2,230 行，不达 5,000。补足波 C2a 真实剩余：BOPTools_AlgoTools 剩余机械方法 + IntTools_Context 完整版 + BOPTools_AlgoTools2D 完整版。
> 参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKBO\BOPTools\{AlgoTools,AlgoTools2D}.cxx` + `IntTools\IntTools_Context.cxx`。
> 基线：5 crate · 143,770 行 · 1,830 测试（topo 1130，Phase 18 后 1163）。

## 缺口 → 方案（各 agent 独立文件，无冲突）

| # | 缺口 | 方案 | 模块（扩展/新建） | 预估 |
|---|---|---|---|---|
| 1 | AlgoTools 剩余机械方法（~40 个未移植） | `ComputeState`（点/边/面对形状 in/out，用 `crate::fclass2d` 2D + 距离 3D）、`MakeConnexityBlock(s)`（复用 `crate::connexity_block`）、`OrientEdgesOnWire`/`OrientFacesOnShell`、`CopyEdge`/`MakeSplitEdge`、`IsMicroEdge`、`IsOpenShell`、`Sense`、`IsHole`、`Dimension`、`CorrectTolerances`、`IsInvertedSolid`、`GetEdgeOff`/`GetFaceOff` | `occt-topo/src/algo_tools.rs`（扩展追加方法） | 2,000 |
| 2 | IntTools_Context 完整版（波 A 只做壳） | `IsPointInFace/IsPointInOnFace`（FClass2d 缓存）、`StatePointFace`、`IsValidBlockForFace(s)`、`IsValidPointForFace(s)`、`IsVertexOnLine`、`ProjectPointOnEdge`、`ComputePE/VE/VF`（点-边/边-边/边-面分类）、`UVBounds`、`clear_cached`。FClass2d/SolidClassifier 缓存用 `crate::fclass2d` | `occt-topo/src/int_tools_full.rs`（新建） | 1,200 |
| 3 | BOPTools_AlgoTools2D 完整版（仅 2 个已移植） | `Make2D`（边→面 pcurve，委托 `pcurve_full`）、`CurveOnSurface`、`EdgeTangent`、`IsEdgeIsoline`、`IntermediatePoint`、`AttachExistingPCurve`、`BuildPCurveForEdgeOnFace`、`PointOnSurface`、`HasCurveOnSurface` | `occt-topo/src/boptools_2d.rs`（新建） | 1,500 |

**合计 ≈ 4,700 行** → 阶段总量 ≈ 6,900。

## OCCT 类 → Rust 映射要点

- **AlgoTools 扩展**（追加到现有 `impl AlgoTools` 块，不破坏已有方法）：
  - `compute_state(shape, p, tol) -> TopAbsState`：面→`fclass2d` 2D 判定（投影到面 UV）；边→距离符号；体→`brep_extrema::is_inside`。
  - `make_connexity_block(shapes) -> ConnexityBlock`（复用 `crate::connexity_block::ConnexityBlock`，先读 API）。
  - `orient_edges_on_wire(wire)`/`orient_faces_on_shell(shell)`：按邻接重定向（法向一致性）。
  - `copy_edge(edge) -> Edge`、`make_split_edge(edge, v1, t1, v2, t2)`（用 `TopoBuilder`）。
  - `is_micro_edge(edge, tol)`（长度 < tol）、`is_open_shell(shell)`（边界边非空）。
  - `sense(f1, f2) -> i32`（同向/反向）、`is_hole(wire, face)`、`dimension(shape) -> usize`。
  - `correct_tolerances(shape, tol)`：遍历顶点/边修正容差。
- **IntTools_Context 完整版**（`int_tools_full.rs`）：
  - `IntToolsContext`：缓存 `HashMap<ShapeKey, FClass2d>`（用 `crate::shape_naming::ShapeId`）。
  - `is_point_in_face(face, p, uv, tol)`（FClass2d::perform）、`state_point_face`、`is_valid_block_for_face(range, face)`。
  - `project_point_on_edge(edge, p) -> Option<f64>`（`brep_extrema` 或距离）。
  - `compute_pe/ve/vf`：点-边/边-边/边-面几何分类（delegating `edge_edge`/`edge_face`，先查签名）。
- **BOPTools_AlgoTools2D**（`boptools_2d.rs`）：
  - `make_2d(edge, face)` = `pcurve_full::make_pcurve_full` 别名；`curve_on_surface(edge, face) -> Option<Arc<dyn Curve2d>>`（查 `tgeometry` pcurves）。
  - `edge_tangent(edge, t) -> GpVec2d`（pcurve 导数）、`is_edge_isoline(edge, face)`、`intermediate_point(t1, t2)`（OCCT 常数 0.43213918，同 `inttools_roots::middle_point`）。
  - `attach_existing_pcurve`（注册到 tgeometry）、`build_pcurve_for_edge_on_face`、`point_on_surface(face, u, v) -> GpPnt`、`has_curve_on_surface`。

## 验证（每 agent 自测）

- **AlgoTools 扩展**：盒边 `is_micro_edge` 对大边 false、微边 true；`compute_state` 对盒面内/外点 In/Out；`make_connexity_block` 对盒 6 面 1 块；`orient_edges_on_wire` 后环闭合；`copy_edge`/`make_split_edge` 端点正确。
- **IntTools_Context**：`is_point_in_face` 与 `fclass2d` 一致；`is_valid_block_for_face` 对范围内/外判定；`project_point_on_edge` 返回正确参数；`compute_pe/ve/vf` 与 `edge_edge`/`edge_face` 对拍。
- **AlgoTools2D**：`make_2d` 与 `pcurve_full` 一致；`edge_tangent` 对盒直线边返回常数方向；`is_edge_isoline` 对盒边 true；`intermediate_point` = 0.43213918。

## 依赖与纪律

- Agent 1 只扩展 algo_tools.rs（追加方法，不改已有）；Agent 2/3 写新文件。`lib.rs` 由 orchestrator 注册新模块。
- 匹配现有风格（`Result<T,String>`、`mod tests`、`TopoShape`、GpPnt）。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target` 避免磁盘满。
- 不碰 bop_builder/bop_curved/boolean_ops 等既有模块。
