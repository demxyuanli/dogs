# Phase 19 Plan — 精确 NURBS 布尔 · 波 C2b-1：PaveFiller 核心（≥7,000 行）

> 精确布尔收官第三波。PaveFiller 全量 10,668 C++ 行（~15,000 Rust），拆 C2b-1（本阶段：主类+求交+分块）+ C2b-2（后续：FF 接入 + bop_builder pipeline）。
> 所有前置依赖已就绪：BOPDS（bopds.rs）、求交核（edge_edge/edge_face/bean_face/int_face_face/intcurvesurface/fclass2d）、AlgoTools（algo_tools.rs）、Context（int_tools_full.rs）、pcurve（pcurve_full.rs）、boptools_2d.rs。
> 参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKBO\BOPAlgo\BOPAlgo_PaveFiller{,_1.._11}.cxx`（12 文件 10,668 行）。
> 基线：5 crate · 148,909 行 · 1,886 测试（topo 1215）。纯增量新模块。

## 缺口 → 方案（各 agent 独立文件，无冲突）

| # | 缺口 | 方案 | 模块（新建） | 预估 |
|---|---|---|---|---|
| 1 | PaveFiller 主类（无） | `PaveFiller`：`Init`（BOPDS 初始化 + 参数数组）+ `Perform`/`PerformInternal`（管线编排：VV→VE→EE→VF→EF→FF + 每步 HasErrors 检查 + UpdatePaveBlocks）+ `Clear` + `SetGlue`/`RepeatIntersection` + error/warning 累积 | `occt-topo/src/pave_filler.rs` | 600 |
| 2 | 求交执行（无） | `PerformVV`（顶点合并，AlgoTools::ComputeVV）+ `PerformVE`/`IntersectVE`（边-顶点，edge_face 定位）+ `PerformEE`（边-边，edge_edge 交点 + MakeNewVertex + pave 放置）+ `PerformVF`/`TreatVerticesEE`（顶点-面，fclass2d 判定）+ `SplitPaveBlocks` + `TreatNewVertices`（新顶点插入 BOPDS） | `occt-topo/src/pave_intersect.rs` | 3,200 |
| 3 | pave 分块与边切分（无） | `MakeBlocks`（PaveBlock 区间分组，BOPDS UpdatePaveBlocks）+ `FilterPavesOnCurves` + `MakeSplitEdges`/`MakeSplitEdge`（AlgoTools::MakeEdge + 新边写入 BOPDS）+ `MakePCurves`（每面 pcurve，pcurve_full）+ `FillPaves`/`FindPaveBlocks` + `ProcessDE` | `occt-topo/src/pave_blocks.rs` | 2,500 |
| 4 | 收缩数据与更新（无） | `FillShrunkData`/`AnalyzeShrunkData`（ShrunkRange 填充）+ `UpdatePaveBlocksWithSDVertices`/`UpdateEdgeTolerance`/`UpdateInterfsWithSDVertices`/`CorrectToleranceOfSE` + `CheckSelfInterference` + `SetNonDestructive` | `occt-topo/src/pave_common.rs` | 1,200 |

**合计 ≈ 7,500 行**。

## OCCT 类 → Rust 映射要点

- **PaveFiller 主类**：持 `ds: BopdsDS`（`crate::bopds::BopdsDS`，先读其 API：append/shape_info/pave_blocks/change_pave_blocks_mut/update_pave_blocks/index 等）、`context: IntToolsContext`（`crate::int_tools_full`）、`fuzzy_value`、`errors/warnings`。
  - `PerformInternal`：`prepare` → `perform_vv` → `perform_ve` → update → `perform_ee` → update → `perform_vf` → update → `perform_ef` → update → `perform_ff`（本阶段 FF 可先空/基础）→ `make_blocks` → `make_pcurves` → `make_split_edges`。每步 `has_errors()` 短路。
- **求交执行**：
  - `PerformVV`：两顶点重合 → AlgoTools::compute_vv + MakeNewVertex → BOPDS append。
  - `PerformVE`：边-顶点 → edge_face 或距离定位参数 t。
  - `PerformEE`：`crate::edge_edge::EdgeEdge`（先读 API）→ 交点 → `algo_tools::make_new_vertex` + `make_edge` → BOPDS PaveBlock 切分。
  - `PerformVF`：顶点-面 → `fclass2d::FClass2d::perform`（in/out/on）。
  - `SplitPaveBlocks`/`TreatNewVertices`：新顶点插入 `ds.change_pave_blocks_mut(edge_idx)`。
- **pave 分块**：`MakeBlocks` 用 PaveBlock 的区间端点排序分组（`inttools_data::IntRange`）；`MakeSplitEdges` 用 `algo_tools::make_edge` 建新边；`MakePCurves` 用 `pcurve_full::make_pcurve_full` 存 `tgeometry`。
- **收缩/更新**：`FillShrunkData` 用 `crate::inttools_range::ShrunkRange`（先读 API）；`CheckSelfInterference` 用 `bopds::BopdsIteratorSI`（自交候选）。

## 验证（每 agent 自测）

- **主类**：`Init` 对两立方体参数 → BOPDS 含 2 体全子形状；`Perform` 后 has_errors false。
- **求交**：两盒交叠 → `PerformEE` 得交点 + 新顶点入 BOPDS；`PerformVV` 盒角重合合并。
- **分块**：`MakeBlocks` 对盒边多个 pave → PaveBlock 分组正确；`MakeSplitEdges` 交点处切边 → 新边入 BOPDS。
- **收缩**：`FillShrunkData` 盒边收缩区间有效；`CheckSelfInterference` 对自交形状报告。
- 用 `crate::brep_extrema::test_box::unit_box`（先确认路径）构造交叠场景。

## 依赖与纪律

- 每个 Agent 只写自己的新模块文件；`lib.rs` 由 orchestrator 注册。
- 匹配现有风格（`Result<T,String>`、`mod tests`、`TopoShape/Face/Edge`、`crate::bopds::BopdsDS`、`IntRange`）。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target` 避免磁盘满。
- 不碰 bop_builder/bop_curved/boolean_ops 等既有模块；不接 bop_builder（C2b-2 统一接）。
