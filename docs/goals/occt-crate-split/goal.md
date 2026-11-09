# OCCT Rust Crate 分拆重构（SRS 落地）

## Objective

把 `docs/goals/OCCT-Rust-Crate-Architecture-SRS.md` 描述的分层 crate 架构落到现有代码上。
现状是 5 个 crate（`occt-core` 49.5k 行、`occt-math` 20.7k、`occt-geom2d` 54.7k、`occt-geom` 117.4k、`occt-topo` 612.5k），
其中 `occt-topo` 一个 crate 就同时承载 TKBRep、TKTopAlgo、TKBO、TKFillet、TKOffset、TKPrim、TKShHealing、TKMesh、
TKBool、Data Exchange 与可视化，正是 SRS §1.2 警告的巨型 crate 形态。

本 tranche 的边界：

- 包含：Rust crate 层的物理拆分（新 crate、模块迁移、Cargo 依赖与 feature 对齐、通过 re-export 保持既有路径可用）。
- 不包含：任何算法控制流改动、任何与 OCCT `.cxx` 对齐的修正、任何为拆分而调参或加启发式的行为变化。
- 不包含：新增测试、新增说明性 Markdown（本目录 `notes/` 与 charter 除外）。

## T001 架构合理性评估结论

### 结论

SRS 的**分层意图成立**，但**规格里写的 14-crate 依赖图不能按原文实现**，必须修正三点后才能落地。

支持 SRS 的证据：

1. `occt-topo` 612,451 行 / 2,269 个 `.rs`，是次大 crate（`occt-geom` 117,408 行）的 5.2 倍，占全仓 855k 行的 72%。
   一个 crate 覆盖 9 个 OCCT toolkit，改 `meshing` 或 `brep_offset` 任何一行都要重编这个 612k 行的 crate，SRS G1/G2 的诉求在现状下完全不成立。
2. `occt-core` 内部 38 个模块**无环**（SCC 计算：0 个强连通分量），`foundation(0 出边) -> gp -> geombase` 三层切分线干净，
   跨层 `crate::` 引用只有 74 + 111 + 20 = 205 处，`pub(crate)` 仅 11 处。底层拆分成本低、风险低，可直接按 SRS P0/P1 执行。
3. `occt-geometry` 不依赖拓扑（`occt-geom`/`occt-geom2d` 均不依赖 `occt-topo`），符合 SRS §4.1.4。

与 SRS 冲突、必须修正的三点：

1. **§4.1.3「算法 crate 之间默认不得互相依赖」不可实现。**
   实测模块级依赖存在，且与 OCCT 自身 toolkit 依赖一致：
   `brep_offset`（TKOffset）→ `chfi3d`（TKFillet）、`bop_builder2`/`topopebrep`（TKBool）、`brep_fill_sweep`（TKOffset/TKPrim）、`meshing`；
   `chfi3d`（TKFillet）→ `topopebrep`、`blend`/`brep_blend`、`intpatch`。
   即真实顺序是 `brep -> algo-core -> 经典布尔(TopOpeBRep) -> prim -> fillet -> offset`，不是 SRS 画的平级兄弟。
   规格需改为：**算法 crate 之间允许 DAG 依赖，禁止环**，并接受 `prim -> fillet -> offset` 这类真实顺序。
2. **§3.3.7 把 `occt-mesh` 放在 `occt-brep` 之上的平级位置不可实现。**
   `fclass2d`、`brep_extrema`、`brep_lib_build_curves3d`、`brep_lib_continuity`、`brep_lib_same_parameter`、`brep_offset`、
   `brep_tools`、`wireframe`、`shhealing` 都直接 `use crate::meshing`。
   mesh 的**数据结构**必须下沉到 brep 层（Poly/Triangulation 已在 `occt-core::poly`），只有 `BRepMesh_IncrementalMesh`
   这类**算法**才能放在上层。
3. **§3.4.1 把整个 Data Exchange 放在算法层之上不可实现。**
   `topopebrep` 依赖 `brep_check`，`brep_offset_api_make_pipe` 依赖 `brep_fill_pipe`，`step` 依赖 `shhealing`，
   `iges_read` 依赖 `brep_fill`。即 Exchange 的**几何构造部分**（BRepTools_ShapeSet、STEP/IGES 的 BRep 构建）
   天然横跨算法层。可拆的是**文件格式读写**（`step`/`iges`/`iges_read`/`xmlcaf`/`bincaf`/`rwmesh`），
   它们不反向被建模层依赖（实测 204 模块 SCC 无任何出边指向这些模块）。

### 另一个决定性事实

对 `occt-topo` 做模块级 SCC 计算，得到**一个 204 个模块的强连通分量**，它同时包含：

- TKBRep 层：`shape`、`tshape`、`abs`、`iterator`、`transform`、`topexp`、`tgeometry`、`brep_tool`、`brep_tools*`
- TKTopAlgo 层：`builder`、`brep_builder_api*`、`brep_lib_*`、`brep_class*`、`brep_extrema*`、`brep_check`、`fclass2d`
- 算法层：`inttools*`、`intpatch*`、`bop*`、`bopds*`、`pave*`、`chfi*`、`brep_offset`、`shhealing`、`meshing`、`primitives`

Rust 的 crate 依赖必须无环，因此**这 204 个模块今天不能拆进两个 crate**。
拆开它们的唯一合法办法是先打断具体回边（例如 `tshape -> brep_tool`、`builder -> brep_lib_make_wire`、
`brep_extrema -> shape_mesh`、`fclass2d -> meshing`），每打断一条都是一次真实的依赖重排，不是重命名。

### 因此采纳的执行顺序

自底向上、只做**能保住 API 的拆分**：被拆出的下层 crate 由原 crate 用 `pub use` 重新导出，
`occt_core::gp::...` 这类既有路径全部继续可用，因此下游 5 个 crate 与 13 个门禁测试零改动。

反面（自顶向下）不成立：把 `step` 移出 `occt-topo` 后，`occt-topo` 无法再 `pub use` 一个依赖自己的 crate，
13 个门禁测试与 27 个 example 的 `occt_topo::step::...` 路径会全部失效。

## Goal Oracle

全部引用仓库内既有基线，不新增测试：

1. `rtk cargo check --manifest-path crates/occt-topo/Cargo.toml --all-targets --offline` 0 error
   （基线：本轮起点实测 0 error / 646 warnings）。
2. `rtk cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --test step_obj_gates` 4 passed，
   且 `rtk cargo run ... --example export_data_obj` 的 v/f 计数与 `data/` 的 OCCT 基线一致。
3. 拆分前后 `occt-core` 的公开符号集合不变：迁移完成后对 `crates/` 全量 `cargo check --all-targets` 0 error，
   且 `crates/occt-geom`、`crates/occt-geom2d`、`crates/occt-topo`、`crates/occt-math` 的 `use occt_core::...` 语句零改动。
4. 每个新 crate 独立 `cargo check -p <crate>` 通过，且 `cargo tree` 中不出现回到上层 crate 的边（无环）。

## 非协商约束

- 移植门禁（`.cursor/rules/occt-port-gates.mdc`）：本 goal 不改控制流、不改算法、不对齐特例；拆分只做位置与可见性调整。
- 不新增测试；不新增说明性 Markdown（`notes/` 除外）；代码注释与字符串使用 ASCII 英文。
- 不改动工作区中既有的未提交改动（`brep_sweep/revol.rs`、`brepfeat/*`、`specs/*`）。
- 命令走 Windows PowerShell，`rtk` 前缀。
- 单文件超 800 行需分拆；拆出的 crate 内出现超 1000 行文件时记录为旁支问题，不在本 tranche 修。
- 完成即停：每个任务写 receipt，旁支问题写进 receipt，不自行扩展范围。

## Current Tranche

1. T001 Scout（已完成）：评估 SRS 合理性、测出 crate/模块依赖与 SCC。
2. T002 PM（已完成）：输出 charter 与画板。
3. T003 Worker：建 workspace 根 `Cargo.toml`；把 `occt-core` 拆成 `occt-foundation` + `occt-gp` + `occt-geombase`，
   `occt-core` 退化为 re-export 门面。
4. T004 Judge：阶段评审，确认 oracle 1-4 与下一刀切点。
5. T005 Worker：拆出 `occt-exchange`（格式读写 + XCAF），重接 consumers。
6. T006 Worker：拆出 `occt-visual`（HLR/场景/rwmesh）与 `occt` 门面（`draw` + 门禁测试 + examples）。
7. T007 Scout：为 TKBRep 抽出列出现存的 204-SCC 回边清单与每条的最小打断方案。
8. T008 Judge：裁定回边打断批次与 `occt-topo`/`occt-brep` 的最终切面。
9. T009 Worker（已完成）：按批准批次打断回边并抽出 `occt-topo`（TKBRep）——批 B1b，17 模块。
10. T010 Worker（已完成）：批 B2，打断 7 条可逆回边后把 `brep_tools{,_uv,_modifier,_wire_explorer}` / `brep_uv_bounds` /
    `brep_bnd_lib` / `brep_surface` 并入 `occt-topology`；`brep_faces` / `brep_measure` 留在 `occt-topo`
    （`brep_measure -> brep_gprop_full` 是 TKTopAlgo 的高斯积分，不是错位的数据层 helper）。
11. T999 Judge（已完成）：逐条复核 Goal Oracle 1-4 全部通过，记录 `full_outcome_complete: true`，
    未拆的 `occt-brep` / 算法层 crate 与 204-SCC 回边清单见 `state.yaml` 的 T999 receipt。

## Run Command

```text
/goal Follow docs/goals/occt-crate-split/goal.md.
```

## Canonical Board

机器真相在 `docs/goals/occt-crate-split/state.yaml`。若本文件与 `state.yaml` 冲突，以 `state.yaml` 为准。
