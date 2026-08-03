# Phase 20 Plan — 精确 NURBS 布尔 · 波 C2b-2a：BOPAlgo_Builder 重建（≥5,000 行）

> 精确布尔收官第四波。PaveFiller（Phase 19）完成求交后，BOPAlgo_Builder 用 DS 的 split 结果重建面/壳/体——替换现有 bop_builder 采样 pipeline。
> 所有依赖已就绪：PaveFiller（pave_*）、BOPDS（bopds.rs 含 ShapesSD/Images）、构建器（builder_area/face/shell_splitter/wire_splitter）、历史（bop_hist.rs）、分类（fclass2d/bean_face）。
> 参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKBO\BOPAlgo\BOPAlgo_Builder{,_1.._4}.cxx`（5,372 行）+ `BOPAlgo_BOP.cxx`（1,784 行）。
> 基线：5 crate · 153,809 行 · 1,942 测试（topo 1276）。纯增量新模块。

## 缺口 → 方案（各 agent 独立文件，无冲突）

| # | 缺口 | 方案 | 模块（新建） | 预估 |
|---|---|---|---|---|
| 1 | Builder 主类（无） | `BopBuilder`：`Perform`（PaveFiller 求交 → BuildBOP → BuildResult → PostTreat）+ `BuildBOP`（Objects/Tools 状态筛选）+ `FillImagesVertices`（ShapesSD 同域 → Images/Origins 三表，用 `crate::bop_hist::BopHistory`）+ `FillImagesEdges` + `BuildResult` + `PostTreat` + `Bop` 包装（Fuse/Cut/Common 状态选择） | `occt-topo/src/bop_builder2.rs` | 2,000 |
| 2 | 面重建（无） | `FillImagesFaces` + `BuildSplitFaces`（每面用 `builder_face::FaceBuilder` 从 split 环重建闭合面）+ `FillSameDomainFaces`（同域面合并） | `occt-topo/src/bop_build_faces.rs` | 2,000 |
| 3 | 容器/内部形状（无） | `FillImagesContainers`/`FillImagesCompounds`（wire/shell/compound 组装）+ `FillInternalVertices`/`FillInternalShapes`（内部孤形状）+ `BuildDraftSolid`（草稿实体从壳组装） | `occt-topo/src/bop_build_common.rs` | 1,500 |

**合计 ≈ 5,500 行**。

## OCCT 类 → Rust 映射要点

- **BopBuilder 主类**：持 `filler: PaveFiller`（`crate::pave_filler`）、`ds: BopdsDS`（`filler.ds()`）、`images/history: BopHistory`（`crate::bop_hist`，先读其 API：add_image/add_modified/images 等）、`origins: HashMap<ShapeKey, Vec<TopoShape>>`、`shapes_sd: HashMap<usize,usize>`（来自 ds）。
  - `Perform`：`PaveFiller::perform` → `BuildBOP(objects, obj_state, tools, tools_state)` → `BuildResult` → `PostTreat`。每步 `has_errors()` 短路。
  - `FillImagesVertices`：遍历 `ds.shapes_sd()`（先读 bopds 是否暴露 ShapesSD）→ images.add_image(v, vsd) + origins 反向表。
  - `FillImagesEdges`：遍历每边 → 其 pave blocks 切分 → 每块 `algo_tools::make_split_edge` → images。
  - `BuildResult`：按结果类型（solid/shell/compound）组装 `images` 里选中的形状。
  - `PostTreat`：容差修正/移除微形状。
- **面重建**：`BuildSplitFaces` 每面 → 从 split 边环（wire_splitter/make_wire）→ `builder_face::FaceBuilder` 重建闭合面（先读其 API：SetFace/perform/areas）。`FillSameDomainFaces`：共享几何的面对 → 保留一个。
- **容器/内部**：`FillImagesContainers` 把 split 面组装成 wire/shell；`FillInternalShapes` 把内部孤边/点分类保留。
- **Bop 包装**：`bop_builder2.rs` 提供 `Fuse/Cut/Common` 函数：`builder_bop(objects, tools, op) -> Result<TopoShape, String>`（状态选择：Fuse=In∪In、Cut=In-Out、Common=In∩In）。

## 验证（每 agent 自测）

- **主类**：两交叠盒 → `Perform` 后无错误 + history 有 split 记录；`FillImagesVertices` 盒角同域映射正确；`BuildResult` 产出有效 shape。
- **面重建**：被切分的盒面 → `BuildSplitFaces` 重建闭合面、面积正确；同域面合并。
- **容器**：`FillImagesContainers` 组装 shell 闭合；内部孤形状分类。
- 用 `crate::brep_extrema::test_box::unit_box`（先确认路径）或 `primitives::BRepPrimBox` 构造交叠场景。最终 `Fuse` 两盒 → 体积 = 1+1-重叠。

## 依赖与纪律

- 每个 Agent 只写自己的新模块文件；`lib.rs` 由 orchestrator 注册。
- 匹配现有风格（`Result<T,String>`、`mod tests`、`TopoShape/Face/Edge`、`crate::bopds::BopdsDS`、`crate::pave_filler::PaveFiller`）。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target` 避免磁盘满。
- 不碰 bop_builder/bop_curved/boolean_ops 等既有模块（新 Builder 独立，后续统一接）。
