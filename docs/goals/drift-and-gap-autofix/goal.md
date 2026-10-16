# 偏移与缺口自动收敛

## Objective

对全仓残留的「范围/语义偏移」（`~/.cursor/projects/d-source-repos-dogs/canvases/occt-code-review-status.canvas.tsx`、`specs/_audit/*.md` 的 SEMANTIC
族条目）与「未移植缺口」（代码里的 UNPORTED 标注、`specs/_audit/_index.md` §9.3、
`specs/_a3n00_gap_analysis.md`）做长程自动修正：逐条回 OCCT `.cxx` 核实，只在 `.cxx` 有同等分支时
改稿，改完跑既有基线复核，通过即提交，然后立即挑下一条，不停在「等指令」。

本 tranche 的边界：

- 包含：`crates/occt-topo`、`crates/occt-geom`、`crates/occt-geom2d`、`crates/occt-core`、
  `crates/occt-math` 内被审计文档登记为 live 的偏移条目与 UNPORTED 条目；对应的审计文档/看板状态订正。
- 不包含：网格偏转/密度类几何产出对齐（那是 `A19` 族，已有独立读数）、IGES/STEP 写侧新特性开发、
  Visualization/Draw/ApplicationFramework 的 lite 替代扩写、`bop_curved`（已删除）。

## Original Request

提交，然后启动goal长程任务对偏移和缺口进行自动修正动作，不再逐步推进

## Intake Summary

- Input shape: `audit`
- Audience: 本项目 owner（OCCT 到 Rust 的几何内核移植）
- Authority: `approved`（可改源码、可本地提交；**不推送、不打 tag、不 rebase、不回滚他人改动**）
- Proof type: `metric`（既有门禁数字）+ `artifact`（审计文档逐条订正）
- Completion proof: live 偏移/缺口清单每条都有结论（已按同分支改完并过门禁，或写明 `.cxx` 无同等分支
  而永久标 UNPORTED）；审计文档与看板状态与代码一致；四条门禁（编译 / `--lib` / `step_obj_gates` /
  `export_data_obj`）保持基线不退。
- Goal oracle: 见下方 Goal Oracle（四条既有基线命令，全部不新增测试）。
- Likely misfire:
  1. 为凑某个模型的顶点/面数/体积，加入面积比、长度滤边、体积门、按 STEP 的特例分支
     （`.cursor/rules/occt-port-gates.mdc` 明令禁止）。
  2. 把审计文档的自述当作已验证结论直接采信，不重新对照 `.cxx`。
  3. 只跑 `cargo check` 就宣称对齐。
  4. 为了让看板"变干净"而把 live 条目直接改写成「已移植」而不改代码。
  5. 把"删掉自创件"当成"已对齐"——只有 `.cxx` 同等分支落地才算闭环。
- Blind spots considered:
  - `~/.cursor/projects/d-source-repos-dogs/canvases/occt-code-review-status.canvas.tsx` 的 `OPEN` / `DOCDRIFT` 条目 state 多为自述，须逐条重新核实真伪
    （本轮已发现两处标题本身有误的先例：`A19 余项`、`GCPnts_AbsComposite`）。
  - `step_obj_gates` 单次约 11-12 分钟，`export_data_obj` 约 15 个模型；两者是本 tranche 的慢门禁，
    不要每个小包都全跑，按包的影响面选跑。
  - 全量 `cargo test` 只在提交前跑一次（项目规则）。
  - 多数 `--lib` / `%TEMP%` 写入类失败属环境问题，须先区分环境失败与真回归。
  - 单文件超 800 行需分拆（项目规则）：`fix_intersecting_wires.rs` 等已知超限文件属结构债，
    可作独立包处理，但不得借分拆之机改行为。
  - 代码注释与字符串必须 ASCII 英文（项目规则）；文档可用中文。
  - 工作树当前干净（HEAD `e77be18e`），每个 Worker 包结束即提交，保持干净。

## Goal Oracle

四条既有基线，全部不新增测试：

1. `rtk cargo check --manifest-path crates/occt-topo/Cargo.toml --offline`（及 geom/geom2d/core/math 四个 crate）
   无 error。
2. `rtk cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --lib` ≥ 1251 passed / 0 failed
   （基线 `1251/0`，只允许单调不减；环境类失败须先证明与改动无关）。
3. `rtk cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --test step_obj_gates` 4/4，
   且 `--nocapture` 下 23 个模型面积比 1.0000。
4. `rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj`
   15/15 `ok`（`-- data/occ data` 口径 23/23），OBJ 顶点/面数与 `data/occ/occ-*.obj` 逐位一致。

任何一条退坡即视为该包未通过：回滚该包的改动（局部），改走 `.cxx` 同等分支重做，而不是放宽门禁。

## Constraints

- 遵 `.cursor/rules/occt-port-gates.mdc` 五道门禁：先编译；对着 `.cxx` 审控制流再改；对齐验证跑真实几何
  对既有基线；完成即停（但本 goal 是「下一条继续」，不停）；任务卡只引用既有基线。
- 禁止谓词/启发式/面积比/长度滤边/体积门这类 OCCT 里不存在的规则。
- 每包一个提交，提交信息写清 `.cxx` 文件与行号。
- 不推送。

## Goal Oracle Cadence

- 每个 Worker 包结束：跑该包影响面内最窄的门禁（编译 + 相关单测 / `--lib`）。
- 每 3-5 个包或触及 BOP/网格/读 STEP 共享管线时：跑 `step_obj_gates` 4/4 + `export_data_obj` 15/15。
- 最终审计：四条门禁全跑一次 + 清单逐条对照。
