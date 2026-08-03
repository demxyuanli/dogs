# Phase 15 Plan — 精确 NURBS 布尔 · 机械波：BOPDS 数据结构 + 面/壳/线构建（≥6,000 行）

> 精确布尔收官第一波：BOPDS 信息中枢 + BOPAlgo_BuilderArea/BuilderFace/ShellSplitter/WireSplitter 构建机械层。
> **不碰数值内核**（IntPatch SS 求交 / FClass2d / PaveFiller 全量留后续波）——本波先建结构与图算法，求交源复用现有 `intpatch`/`inttools`。
> 参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKBO\`（BOPDS 28 文件 6,242 行 / BOPAlgo 相关 / BOPTools 17 文件 7,972 行）。
> 基线：5 crate · 127,556 行 · 1,656 测试。纯增量新模块，不动既有 bop_builder/bop_curved/boolean_ops。

## 缺口 → 方案（各 agent 独立文件，无冲突）

| # | 缺口 | 方案 | 模块（新建） | 预估 |
|---|---|---|---|---|
| 1 | BOPDS 信息中枢（无等价物） | `BOPDS_DS`（参数数组/形状索引/pave block/face info/rank/range）+ `PaveBlock` + `Pave` + `CommonBlock` + `FaceInfo` + `ShapeInfo` + `IndexRange` + `Interf` + `Iterator`/`SubIterator`/`IteratorSI` | `occt-topo/src/bopds.rs` | 3,500 |
| 2 | 闭合面构建（无） | `BOPAlgo_BuilderArea`（基类：shapes/loops/areas）+ `BOPAlgo_BuilderFace`（从 split faces 重建闭合 face，含线框→面） | `occt-topo/src/builder_area.rs` + `builder_face.rs` | 1,800 |
| 3 | 闭合 shell/wire 图回路 | `BOPAlgo_ShellSplitter`（连通块→shell）+ `BOPTools_ConnexityBlock`（连通性块）+ `BOPAlgo_WireSplitter`（线框→闭合 wire，`SplitBlock`/`MakeWire`） | `occt-topo/src/shell_splitter.rs` + `wire_splitter.rs` + `connexity_block.rs` | 1,800 |
| 4 | 选项/告警/历史（无） | `BOPAlgo_Options`（fuzzy/parallel/report）+ `BOPAlgo_Alerts`（错误/警告）+ 历史/命名侧表（shape_naming 复用） | `occt-topo/src/bopalgo_options.rs` + `bop_hist.rs` | 1,500 |

**合计 ≈ 8,600 行**，远超 5,000 目标。

## OCCT 类 → Rust 映射要点

- **BOPDS_DS**：核心是 `NCollection_IndexedMapOfShape` → Rust `IndexedMap`（`occt-core/src/kernel/containers.rs` 已有）+ `BOPDS_ShapeInfo` 数组（DynamicArray → `Vec`）+ `PaveBlock` 列表（每边一组）。`Append(shape) -> index`、`ShapeInfo(i)`、`Index(shape)`、`Rank(i)`、`HasPaveBlocks/ChangePaveBlocks`、`UpdatePaveBlocks`、`UpdateCommonBlock`、`FaceInfoPool`。形状索引用 `shape_naming.rs` 的 ShapeId/geometric_hash。
- **PaveBlock**：`PaveBlock(edge_index, first, last)` + 分块；`Pave`：边上的参数点（index, param）；`CommonBlock`：多边共享块（同一几何区间的若干边）；`FaceInfo`：面的边信息（pave blocks 所属面）。
- **BuilderArea**（虚基类）：`SetShapes/PerformShapesToAvoid/PerformLoops/PerformAreas/PerformInternalShapes` → Rust 用 trait `AreaBuilder` + 具体 `FaceBuilder`/`ShellBuilder`。
- **ShellSplitter**：`AddStartElement`/`Perform`/`Shells` + `SplitBlock`（把 ConnexityBlock 按连通性拆成闭合 shell）；`WireSplitter`：`SetWES(线框边集)/Perform/MakeWires/MakeWire`（边图回路闭合）。
- **Options**：`fuzzy_value/run_parallel/report/errors/warnings`；Alerts 是字符串集合。
- **接线**：本波产出独立结构 + 构建器，`bop_builder.rs` **不接**（留待数值波统一接入）；但 `bopds.rs` 提供 `IndexedMap` 复用 kernel/containers.rs。

## 差分验证（每波必做）

- **BOPDS**：单位盒 8 顶 12 边 6 面 1 壳 1 体 → DS 索引完整（Append 计数、Index 往返、Rank、Range）；插入/移除 pave block 后 UpdatePaveBlocks 一致。
- **ShellSplitter/WireSplitter**：6 个四边形边集合重组成 6 个闭合环；单位盒 shell 拆连通块 → 1 个闭合 shell（连通性验证）；分离两个立方体 → 2 个 shell。
- **BuilderFace**：带 1 内孔的平面 face（环+子环）重建 → 闭合 face，面面积与输入一致。
- **Options**：fuzzy 设置/清空、错误告警累积。

## 依赖与纪律

- 每个 Agent 只写自己的新模块文件（各不同文件，无冲突）；`lib.rs` 由 orchestrator 统一注册。
- 匹配现有风格（`Result<T,String>`、`mod tests`、TopoShape/Face/Edge/Wire/Shell/Solid、`occt_core::kernel::containers::IndexedMap`、`shape_naming::ShapeId`）。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target` 避免并行文件锁。
- 不碰 bop_builder/bop_curved/boolean_ops/meshing 等既有模块。
