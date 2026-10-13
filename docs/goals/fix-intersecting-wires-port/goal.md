# FixIntersectingWires 端口对齐审计与收敛

## Objective

对照 OCCT `ShapeFix_IntersectionTool::FixIntersectingWires`（`D:\source\OCCT-src\src\ModelingAlgorithms\TKShHealing\ShapeFix\ShapeFix_IntersectionTool.cxx:1834-2530`，共 697 行）审计 Rust 端口 `crates/occt-topo/src/shhealing/intersection_tool/fix_intersecting_wires.rs`（916 行）。逐段核实控制流，仅修复被证实的偏差（且修复必须走 OCCT 同一 `.cxx` 分支），并保持 STEP→OBJ 门禁通过。

本 tranche 的边界：

- 包含：`FixIntersectingWires` 本体（cxx:1834-2530）、它直接调用的 `intersection_tool.rs` 辅助函数、`ShapeFix_Face::FixIntersectingWires` 包装（`ShapeFix_Face.cxx:2821-2825`）、`read_topology.rs` 的 `fix_face_wire_pairs` 接线（`cxx:676-685`）。
- 不包含：`FixSelfIntersectWire`（cxx:1029-1833）、`ShapeFix_Face` 中 FixMissingSeam / ComposeShell 相关的未移植段、网格管线。

## Original Request

根据FixIntersectingWires ported 内容做一个goal任务

## Intake Summary

- Input shape: `specific`
- Audience: 本项目 owner（OCCT 到 Rust 的 STEP 读取与网格管线移植）
- Authority: `requested`（审计与收敛偏差；不提交、不推送、不回滚他人改动）
- Proof type: `metric`（既有门禁数字）+ `artifact`（逐段对照表）
- Completion proof: 逐段对照表中每个 cxx 段都有 Rust 对应位置或明确"未移植"标注；被证实的偏差已按同分支修复；`step_obj_gates` 4/4 通过；`export_data_obj` 计数与 OCCT 基线一致。
- Goal oracle: 见下方 Goal Oracle。
- Likely misfire:
  1. 把规格 §9.633 的"逐段对照"当作已验证结论直接采信，而不重新核对。
  2. 为凑某个 STEP 的顶点/面数，加入面积比、长度滤边、体积门或按模型的特例分支。
  3. 只跑 `cargo check` 就宣称对齐。
  4. 把 `shape_fix_face.rs:1002`、`:1223` 的 UNPORTED 注释（引用的是 `ShapeFix_Face.cxx` 的 FixMissingSeam / ComposeShell 段）误算为本 tranche 缺口。
- Blind spots considered:
  - §9.633 的逐段对照是文档自述，不是独立证据，需要 Scout 重新逐段核对。
  - `intersection_tool.rs` 为 1465 行，超过项目"单文件超过 800 行需分拆"规则，属于本 tranche 相邻的结构债。
  - `step_obj_gates` 单次约 11-12 分钟；全量 `cargo test` 不属于本 tranche 的常规门禁（项目规则：仅提交前跑一次）。
  - 工作区存在大量未提交改动（含上一轮 `is_outer_bound` 修复，见 `crates/occt-topo/src/shape_analysis.rs`）。Worker 只能改 `allowed_files`，不得回滚或格式化无关文件。
  - `cargo test --lib` 中 26 项因写入 `%TEMP%` 被拒（`os error 5`）失败，属环境问题，不作为本 tranche 的回归判据。
  - 项目规则要求代码注释与字符串使用 ASCII 英文；Worker 新增内容须遵守。
- Existing plan facts:
  - 规格 `specs/_a3n00_gap_analysis.md` §9.633：11 个 cxx 行段的逐段对照表，逐段对照表另见 canvas `intersecting-wires-port`。
  - §9.633 C 节保真要点 1-6（`Context()->Apply` 只在 `FixWiresTwoCoincEdges` 里做；过期的 `aBox1`；`ChangeValue` 别名写回；`IsModified` 死写；`E1 == E2` 的等价判据；重建面的位置先于子形状）。

## Goal Oracle

本 goal 的 oracle 为以下四条，全部基于仓库内既有基线，不新增测试：

1. `rtk cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --test step_obj_gates` 结果为 4 passed（既有测试 `step_obj_parity_bboxes_match_occt`、`step_obj_area_matches_occt`、validity 与 export 三项，容差见 `crates/occt-topo/tests/common/mod.rs`）。
2. `rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj` 对 `data/*.step` 全部输出 `ok`，且计数与 OCCT 基线一致。已知 OCCT 基线（来自 `step_obj_parity` 的 `occ v=… f=…` 输出）：Cube 24/12，Cone 195/301，Cylinder 146/140，Sphere 642/1244，Torus 1369/2592，HoledPlate 180/128，OffsetPlaneHoleEdge 48/32，rev 104/92，Offset 712/892，Shape 6150/11372，Shape-1 3343/4336（格式为 v/f）。
3. 逐段对照表：`ShapeFix_IntersectionTool.cxx:1834-2530` 的每个段（至少覆盖 1842-1862、1864-1872、1874-1888、1889-1909、1910-2045、2047-2482、2483-2505、2508-2529）都标注为 faithful / divergent（附修复） / unported（附 cxx 行号与理由）。
4. `rtk cargo check --manifest-path crates/occt-topo/Cargo.toml --all-targets --offline` 0 error。

门禁判据只能引用上述既有基线。出现失配即停，报告差异，不得用特例或启发式掩盖。

## Goal Kind

`specific`

## Current Tranche

按顺序推进：

1. T001 Scout：逐段审计，产出对照表与偏差清单（只读）。
2. T002 Judge：从偏差清单选出最大安全切片，给出 T003 的 exact `allowed_files` / `verify` / `stop_if`，并裁定 `intersection_tool.rs` 拆分是否纳入 T004。
3. T003 Worker：执行 T002 批准的偏差修复。若 T001 未发现已证实偏差，则以 receipt 记录"无偏差"并跳过。
4. T004 Worker：仅当 T002 裁定拆分时执行，行为保持的文件拆分。
5. T999 Judge：最终审计，逐条映射 Goal Oracle，记录 `full_outcome_complete`。

## Non-Negotiable Constraints

- 移植门禁（`.cursor/rules/occt-port-gates.mdc`）：改完先 `cargo check`；对着 `.cxx` 审控制流再改；有同等分支才改，没有就标"未移植"并注释写清 OCCT 文件与行号；不加 OCCT 里不存在的谓词、启发式、面积比、长度滤边、体积门；不为某个 STEP 调参或写特例。
- 不为对齐新写单元测试；不用 `StubBopBuilder` 一类假输入证明对齐。
- 不自动生成测试程序；不新增说明性 Markdown（本目录 `notes/` 下的 receipt 除外）。
- 代码注释与字符串使用 ASCII 英文，不使用中文注释。
- 命令按 Windows PowerShell 执行，并使用 `rtk` 前缀；C++ 相关命令不适用于本 tranche。
- Worker 只能写 `allowed_files`；不得回滚、格式化或提交无关改动；不 commit、不 push。
- 单文件超过 800 行需分拆，超过 1000 行提示分拆。
- 完成即停：每个任务完成后写 receipt，旁支问题写入 receipt，不自行扩展范围。

## Stop Rule

只有 T999 的最终审计 receipt 同时满足 Goal Oracle 四条并记录 `full_outcome_complete: true`，才可结束。

以下情况不得停止整个 goal：

- 某个偏差需要人工判断（改为 blocked 并写 receipt，继续其它本地安全切片）。
- 门禁因环境失败（先判断是否为环境问题，并记录；不得把环境失败当作对齐成功）。

门禁失配（step_obj_gates 失败或计数偏离基线）时，必须回到 Scout/Judge 定位，不得直接改特例通过。

## Slice Sizing

- T001 只读，覆盖整个 cxx 段范围，不拆成逐段小任务。
- T003 的切片由 T002 决定；应取同一 `.cxx` 分支的完整修复，不做零碎的单点改动。

## Canonical Board

机器真相位于 `docs/goals/fix-intersecting-wires-port/state.yaml`。若本文件与 `state.yaml` 冲突，以 `state.yaml` 为准。

## Run Command

```text
/goal Follow docs/goals/fix-intersecting-wires-port/goal.md.
```

## PM Loop

每次 `/goal` 续跑时：

1. 读取本 charter。
2. 读取 `state.yaml`。
3. 只处理当前 active 任务。
4. 按任务类型分派 Scout / Judge / Worker / PM。
5. 写入紧凑 receipt，更新 `state.yaml`。
6. 若仍有安全的本地工作，选择下一个最大可逆 Worker 切片继续。
7. 仅在阶段、风险、门禁失败或最终审计时安排 Judge。
8. 只有 T999 的 Judge 审计 receipt 可将 goal 标记为 done。
