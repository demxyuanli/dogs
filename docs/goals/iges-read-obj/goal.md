# IGES 读侧：读入 IGES 并输出 OBJ

## Objective

把 OCCT IGES 读侧移植到 Rust：读入 IGES，重建形状，再走现有网格路径写出 OBJ。第一批证明只用写侧已经与 `_igesdump/oracle` 参数段一致的模型，对照这些模型既有的 STEP 网格基线。

## Original Request

现在开始实现读侧，读iges然后输出obj。下一步先给出任务看板，再启动长程来完成。

## Intake Summary

- Input shape: `specific`
- Audience: 本项目 owner（OCCT 到 Rust 的几何移植）
- Authority: `requested`（实现读侧并出 OBJ；不提交、不推送）
- Proof type: `artifact`（OBJ）+ `metric`（既有 STEP 网格基线的包围盒与面积）
- Completion proof: 读 `_igesdump/oracle` 中 Cube、Sphere、Cone、Extrusion、HoledPlate、rev、Offset 的 IGES，用与 `export_data_obj` 相同的偏转 0.1 网格路径写出 OBJ；包围盒与面积落在 `step_obj_parity` 已有容差内；`cargo check` 0 error。
- Goal oracle: 见下方 Goal Oracle。
- Likely misfire:
  1. 用端口自己写出的、仍有几何差的 IGES 当规格，反过来改读侧去迁就错误写侧。
  2. 为凑某个模型的顶点数加入面积比、长度滤边、体积门或按模型特例。
  3. 只跑 `cargo check` 或只解析文件头就宣称读侧完成。
  4. 新写单元测试代替对拍既有 OBJ 基线。
  5. 把写侧未闭合的 Shape / Shape-1 / Shape-2 / screw / Torus / Cylinder 差带进本 tranche。
- Blind spots considered:
  - 回转面经 IGES 120 重建后，网格顶点数可能与 STEP 直读不同；判据用既有包围盒和面积容差，顶点数只报告、不用特例掰齐。
  - oracle 是 OCCT 写出的 IGES，读侧规格是 `IGESControl_Reader` / `IGESToBRep`，不是端口写侧。
  - 工作区有未提交的网格/愈合改动和大量 `_igesdump` 临时文件。Worker 不得回滚或收进这些文件。
  - `iges.rs` 已超过 1000 行。读侧不要再堆进该文件，新代码按 OCCT 类拆文件。
  - 全量 `cargo test` 只在提交前跑。本 tranche 不提交。
- Existing plan facts:
  - 写侧已推送。Faces 模式是 120+124，402 form 为 1，全局段版本 11、精度 308/15。
  - 与 oracle 参数段一致：Cube、Sphere、Cone、Extrusion、HoledPlate、rev、Offset。
  - 既有 STEP 网格计数（`step_obj_parity`）：Cube 24/12，Cone 195/301，Sphere 642/1244，HoledPlate 180/128，rev 104/92，Offset 712/892。Extrusion 的计数以当次 `export_data_obj` 对 `data/Extrusion.step` 的输出为准，不预写数字。

## Goal Oracle

本 goal 的 oracle 为以下三条，全部基于仓库内既有基线，不新增测试：

1. `rtk cargo check --manifest-path crates/occt-topo/Cargo.toml --offline` 0 error。
2. 读 `_igesdump/oracle/<name>.iges`（Cube、Sphere、Cone、Extrusion、HoledPlate、rev、Offset），经读侧重建形状后，用与 `export_data_obj` 相同的偏转 0.1 写出 OBJ。
3. 每个 OBJ 的包围盒与表面积对照该模型既有 STEP 网格（`data/output/<stem>.obj` 或 `step_obj_parity` 已记录的基线），容差只用 `crates/occt-topo/tests/common/mod.rs` 里已有的容差。失配即停，报告差异。

门禁判据只能引用上述既有基线。

## Goal Kind

`specific`

## Current Tranche

按顺序推进，直到 oracle 三条同时成立：

1. T001 Scout：只读，定位 OCCT 读侧入口、这 7 个 oracle 文件实际用到的实体，以及现有 OBJ 导出入口。
2. T002 Judge：选出最大的一块可逆实现，写明 `allowed_files` / `verify` / `stop_if`。
3. T003 Worker：实现这块。之后若 oracle 仍未满足，PM 追加下一块 Worker，而不是停在“骨架已好”。
4. T999 Judge：把三条 oracle 映射到 receipt，记录 `full_outcome_complete`。

## Non-Negotiable Constraints

- 移植门禁：改完先 `cargo check`；对着 `.cxx` 审控制流再改；有同等分支才改，没有就标未移植并写清 OCCT 文件与行号；不加 OCCT 里不存在的谓词、启发式、面积比、长度滤边、体积门。
- 不为对齐新写单元测试。不新增说明性 Markdown（本目录 `notes/` 下的 receipt 除外）。
- 代码注释与字符串使用 ASCII 英文。
- 命令按 Windows 执行，并使用 `rtk` 前缀。PowerShell 会吃掉 `$`，环境变量用 `cmd /c "set VAR=...&& rtk cargo ..."`。
- Worker 只能写 `allowed_files`。不回滚、不提交、不推送无关改动。
- 读侧新文件按 OCCT 类命名，禁止 `pNN.rs`。单文件超过 800 行要拆，超过 1000 行必须拆。
- 改完代码后运行 `graphify update .`。不要运行 `graphify label`。
- 完成即停：旁支写入 receipt，不自行扩展到写侧残留差。

## Stop Rule

只有 T999 的最终审计 receipt 同时满足 Goal Oracle 三条并记录 `full_outcome_complete: true`，才可结束。

以下情况不得停止整个 goal：

- 只完成了勘察或判定，还没有读入并写出 OBJ。
- 只读通了 Cube，其余 6 个模型仍安全可做。
- 某一实体类型在 OCCT 中没有对应分支。该类型标未移植，其余类型继续。

## Slice Sizing

一块任务应当是一条能读进文件并朝 OBJ 推进的垂直切片，而不是单独一个解析辅助函数。

## Canonical Board

`docs/goals/iges-read-obj/state.yaml`

章程与 `state.yaml` 不一致时，以 `state.yaml` 的任务状态、活动任务和 receipt 为准。

## Run Command

```text
/goal Follow docs/goals/iges-read-obj/goal.md.
```
