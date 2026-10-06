# dogs — OCCT 的 Rust 移植（STEP → OBJ 网格对齐）

本仓库把 **OpenCT 8.0.0**（下称 OCCT）按类、按控制流逐步移植到 Rust，目标是让
**STEP 读入 → 整形（ShapeFix）→ 网格化 → OBJ 导出**这条链在**真实几何**上与 OCCT 对齐：
逐模型的网格面积比、顶点/三角计数、以及逐面的结构（wire 数、每环边数）都要能对拍。

判断对齐“靠对拍现有基线”，不靠堆断言；`.cxx` 是规格，移植代码必须能指到对应的 OCCT 分支。

---

## 目录结构

```
crates/                 Rust 移植（5 个 crate，均 v0.1.0，**没有根 workspace**）
  occt-core/            基础类型与常量（gp_*、Precision/Confusion、TopoShape…）
  occt-math/            数值（求解、极值、插值…）
  occt-geom/            3D 几何（曲线/曲面）
  occt-geom2d/          2D 几何（pcurve）
  occt-topo/            拓扑 + STEP 读写 + ShapeFix/ShapeAnalysis + 网格管线（主战场）
data/                   模型与参考基线
  occ/*.stp             输入模型；occ-*.obj 为 **OCCT 参考网格（GT）**
  iges/                 IGES 样例
specs/                  规格与对拍记录
  board.canvas.tsx      任务看板（卡面/进度/nextAction，人机共编）
  _a3n00_gap_analysis.md  a3n00 缺口分析（§9.x 逐轮记录，含每次实测数字与口径）
  occt_probe/           **OCCT 侧**只读探针（C++，需本机 OCCT 构建）
tools/, engine/, docs/  辅助脚本 / 引擎侧代码 / 文档
output/                 OBJ 导出目录（`export_data_obj` 的产物）
CLAUDE.md               **本仓库的工作纪律**（门禁、命名、TEMP 插桩、提交要求）——改代码前先读
.target-gate/           gitignored 的临时产物（门禁输出、对拍脚本、探针 dump）
```

---

## 环境要求

- **Rust**（stable；无 `rust-toolchain.toml`，用本机默认工具链）
- 本仓库的文档与历史命令使用 `rtk cargo …` 包装器；**没有 `rtk` 时直接用 `cargo` 即可**，参数完全等价
- 依赖已 vendored/缓存，日常用 `--offline` 构建
- 只有需要跑 **OCCT 侧探针**（`specs/occt_probe/`）时才需要本机 OCCT 构建；
  仓库里的示例路径是 `D:\source\occt-8.0.0`（含 `env.bat`、`win64/vc14/{bin,lib}`），按你的机器改

---

## 构建与测试

仓库没有根 workspace，所以**必须用 `--manifest-path`**：

```bash
# 1) 先编译（第一条门禁）
cargo check --manifest-path crates/occt-topo/Cargo.toml --offline

# 2) 单元测试（提交前跑一次；当前基线 1281 passed / 0 failed）
cargo test  --manifest-path crates/occt-topo/Cargo.toml --offline --lib

# 3) 逐模型网格对齐门禁（area ratio / f-ratio，当前 5 passed / 0 failed）
cargo test  --manifest-path crates/occt-topo/Cargo.toml --offline --test step_obj_gates -- --nocapture

# 4) 导出 OBJ（data/*.step → data/output/<stem>.obj，偏转 0.1）
cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj
```

其余集成测试见 `crates/occt-topo/tests/`（`phase*_integration.rs`、`bop_builder2_boss.rs` 等）。

---

## 对齐工作流（摘自 `CLAUDE.md`，务必遵守）

1. **改完先编译。** 编译不过就先修编译，再谈对齐。
2. **对着 `.cxx` 审控制流再改。** 有缺口时先定位“缺了哪段 OCCT 控制流”，修复走**同一套 `.cxx` 分支**；
   不为某个 STEP / 特征 / 特例调参（网格管线是共享的，特例补丁会带偏已对齐的模型）。
   没有同等分支就标成 **UNPORTED**（注释写清 OCCT 文件与行号），不加 OCCT 里不存在的谓词/启发式。
3. **对齐验证 = 跑当前实例的真实几何，对照现有基线。** 失配即停并报告差异；
   **不要**为对齐新写单元测试，**不要**改断言/基线/`area_tol`（任务卡的门禁只能引用现有基线）。
4. **完成即停。** 门禁满足 → 报告结果与旁支问题 → 停，等下一步指令。
5. **TEMP 插桩必须可撤除**：一律用 `edit` 反向替换（**禁止 `git checkout`** 甩掉工作树），
   提交前用 `git status --porcelain` / `git diff --stat` 证明库代码已还原。
6. **命名**：禁止 `pNN.rs` 这类无语义文件名；拆模块用 `<module>/mod.rs` + `mod <主题>;`；
   引用源码用真实文件名与行号（如 `step/read_topology.rs:476`）。

---

## 仪器（只读探针）

**端口侧**（`crates/occt-topo/examples/`，`cargo run --example <名字>`）：

| 探针 | 用途 |
|---|---|
| `export_data_obj` | 批量导出 OBJ（`data/*.step` → `data/output/<stem>.obj`） |
| `zz_probe_a3n00` | a3n00 逐面读数：`--fstats`（每面 mv/mt/bbox）、`--fdump`、`--fixms`、`--ecensus`（逐面 wire/边数） |
| `zz_uv_feed` | 面/模型结构：`--model <f>`、`--ids`（逐面 wire 数 + 直方图） |
| `zz_seam_fix` | 单面 `FixMissingSeam` 前后结构（`<stp> <face>`：BEFORE/result/HEALED） |
| `zz_share_probe` / `zz_pcurve_key_probe` / `cone_dbg` / `linkrods_dbg` / `iges_check` | 专题排查 |

**OCCT 侧**（`specs/occt_probe/`，需本机 OCCT）：

```bash
cmd /c "specs\occt_probe\run_dbg.bat data\occ\a3n00.stp"          # 逐面 census（整形后）
cmd /c "specs\occt_probe\run_dbg.bat data\occ\a3n00.stp noop"     # 关 ShapeProcess（未整形）
cmd /c "specs\occt_probe\run_dbg.bat data\occ\a3n00.stp wdump 113" # 按面序号 dump 边/端点
specs\occt_probe\probe.bat data\occ\a3n00.stp --facestats 1.076007 0.349066   # 同参数逐面 nodes/tris
```

> 对拍铁律：**跨侧配对必须按几何（bbox/顶点集合），不得按面号**；
> 逐面读数必须**在同一次运行内取**（跨探针会出现不一致）。

---

## 数据与产物

- 输入：`data/occ/*.stp`（a3n00、acs10、ATU01038、bottom、motoc、T0M、top、TDB…）
- 参考（GT）：`data/occ/occ-*.obj` —— 逐模型面积比与 f-ratio 的对照基准
- 产物：`data/output/<stem>.obj`（已 gitignore）—— STEP→OBJ 唯一导出目录（`export_data_obj` 与 `step_obj_gates` 共用）
- 临时：`.target-gate/`（已 gitignore）存放门禁输出、对拍脚本与探针 dump

---

## 文档索引

- `CLAUDE.md` —— 工作纪律（**改代码前必读**）
- `specs/board.canvas.tsx` —— 任务看板：卡面、状态、`nextAction`、活动流
- `specs/_a3n00_gap_analysis.md` —— a3n00 缺口分析：§9.x 逐轮记录**实测数字与其口径**、
  已排除项、以及每条结论对应的 `.cxx` 行号
- `docs/` —— 阶段性审查文档

---

## 现状（如实）

- **已对齐**：STEP 读入/拓扑组装、`ShapeFix_Face::FixMissingSeam`（导入期接回，
  a3n00 面积比 0.8627 → **0.8996**，带倒角法兰 16 面 2 wires → 1 wire/5 边）、
  网格侧 ModelBuilder/Healer/FaceChecker 主流程；逐模型门禁 `step_obj_gates` 5/5，
  其中 a3n00 0.8996、T0M 0.9987、acs10 0.9846。
- **已知最大缺口**：**没有 `ShapeProcess`/`FixShape` 驱动器**（端口在 reader 里按面直调子步）。
  实测关掉该算子时 OCCT 有 **131/226 面无法网格化**，说明它是承重环节；端口仅缺分派器与参数下发。
- **未结个案**：a3n00 的“螺母斜切面”（端口 F113）仍为 **4 wires / 空面**。
  已定位到：`fix_missing_seam` 在该面上把 28 条边分成了 **6 条 wire / 5 张面**（其中 1 张 2-wire 面 = 6e+7e，
  另 4 张单 wire 面 = 5/4/3/3），而 OCCT 是 **2 条 wire（22+6）/ 1 张面**——**总边数相同（28），这一步不丢边**；
  结果被包成 `Shell(5)`，而 reader 只接受 `Face` ⇒ 整包被丢弃。打包与剪枝语义均已证明与 `.cxx` 同构
  （更正记录见 §9.450：早期文档里的「5 块单 wire 补丁 6/5/4/3/3」是解析错误），分歧指向
  `SplitWires/BreakWires` 的 **wire 分组**。全过程与后续路径见 `specs/_a3n00_gap_analysis.md` §9.370–§9.450 与看板 T-101。
