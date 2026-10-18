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
   引用源码用真实文件名与行号（如 `step/read_topology/translate_tool.rs:38`）。

---

## 仪器（只读探针）

**端口侧**（`crates/occt-topo/examples/`，`cargo run --example <名字>`）：

| 探针 | 用途 |
|---|---|
| `export_data_obj` | 批量导出 OBJ（`data/*.step` → `data/output/<stem>.obj`） |
| `model_probe` | a3n00 逐面读数：`--fstats`（每面 mv/mt/bbox）、`--fdump`、`--fixms`、`--ecensus`（逐面 wire/边数） |
| `mesh_uv_feed_probe` | 面/模型结构：`--model <f>`、`--ids`（逐面 wire 数 + 直方图） |
| `seam_fix_probe` | 单面 `FixMissingSeam` 前后结构（`<stp> <face>`：BEFORE/result/HEALED） |
| `seam_pcurve_share_probe` / `pcurve_key_probe` / `cone_classify_check` / `mesh_diagnostics` / `iges_check` | 专题排查。`iges_check` 写 IGES 并做卡片结构检查：参数是 `data/<name>.step`，含 `/` 的参数按仓库根路径读 `.step` 或 `.stp`。`IGES_DUMP` 指向目录时写出 `<参数名>.iges` |

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

- 输入：`data/occ/*.stp` 与 `data/occ/*.step`（a3n00、acs10、ATU01038、bottom、motoc、T0M、top、TDB），以及 `data/*.step`（含 `linkrods.step`）
- 参考（GT）：`data/occ/occ-*.obj` —— 逐模型面积比与 f-ratio 的对照基准
- glTF 参考（GT）：`data/occ-{cube,cone,linkrods}.gltf` + `.bin` —— 由 `data/_occ_ref_gltf.tcl` 在 DRAWEXE 里重生成（`ReadStep` → `incmesh` 用 `bbox_max*0.001*4`、`-angular 20` → `WriteGltf`，默认非二进制、`systemCoordSys Zup`、`nodeNameFormat instOrProd`、`meshNameFormat product`；设 `OCC_FORCE_UV=1` 则加 `-forceUVExport`，用于逐节点 UV 对拍）。对拍用 `rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example gltf1_read_probe`：`model_to_gltf_caf` 的 Cube/Cone 输出与基线逐字节相等。linkrods 的 accessor/bufferView/mesh 布局一致，几何残差读数（按显式 bin 路径比对，**不要**让脚本按 `buffers[0].uri` 回退路径，否则会读到 OCCT 的 bin）：37 个 primitive 里位置 21 个逐位相同、10 个差 ≤1e-4、5 个差 1e-4..7.5e-4；索引 25 个不同但三角形按 1e-3 取整做集合比对只 5/5078 不匹配；法向 29 个不同。对拍脚本所在 `_tmp/` 已 gitignore
- glTF 读入：`rwmesh::read_gltf_scene_with_options` / `read_gltf_model`（`RWGltf_CafReader`，含 glTF 1.0/2.0、GLB、bufferView/data URI 贴图、late-data 两阶段）。文档式回读见 `gltf1_read_probe` 的 doc OBJ 段
- 产物：`data/output/<stem>.obj` + `<stem>.mtl`（已 gitignore）—— STEP→OBJ 导出目录，由 `export_data_obj` 独占写入（带颜色时同时写 `.mtl` sidecar）
- 贴图目录：`data/output/<stem>_textures/` —— 样式带 `BaseColorTexture` 时由 `RWMesh_MaterialMap::CreateTextureFolder` / `CopyTexture` 写出（纯图片文件按原名拷贝，块/缓冲区贴图按 `<n>.<ext>` 重命名），`.mtl` 里对应 `map_Kd <stem>_textures/<file>`。STEP 读侧当前不带贴图，故该目录只在实际带贴图的样式上出现
- 门禁产物：`data/output/gates/<stem>.obj` + `<stem>.mtl` —— `step_obj_gates` 写入，**与导出目录分离**，避免跑门禁时覆盖上面的 CAF 产物；门禁走同一 `model_to_obj_caf` 管线，同一模型两边文件逐字节相同
- 读回校验：`rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example rwmesh_obj_readback`（默认扫 `data/output/gates/`，也可传 OBJ 文件或目录）。用 `rwmesh::read_obj_scene` 把上面的产物读回，逐节点核对 `usemtl` 是否挂到节点、`Ka`/`Kd`/`Ks`/`Ns`/`Tr` 是否与 `.mtl` 写出值一致（写侧没写的字段必须落在 `RWObj_Material` 默认值上）、`map_Kd`/`map_Ks`/`map_Bump` 是否按 `processTexturePath` 解析到同一路径（含绝对路径的 `findRelativePath` 回退）、三角形下标是否越界；任一不符退出码 1
- IGES 对拍：`_igesdump/oracle/` 是 OCCT 写出的 14 个样例。端口导出用 `iges_check`，`IGES_DUMP` 只作人工对照，不入库
  - **降级 oracle：`_igesdump/oracle/Torus.iges` 只作写侧参数对照，不作"读回正确输出"基准。** 它复现的是 **OCCT 8.0.0 Faces 模式往返缺陷**（2026-10-12 定位，详见 `docs/goals/iges-read-full/state.yaml` T027）：Faces 模式（`IGESControl_Writer("MM", 0)`，`GeomToIGES_GeomSurface.cxx:576-598` 的 `myBRepMode && myAnalytic` 为假）把解析环面写成 `144 + 120` 回转面 + 相位不自洽的 `142`；OCCT 8.0.0 自己读回得到 **1 面 / 2 wires / v 跨 4π ≈ 302 节点 146 三角**（面积 67.24，真值 789.57），而 OCCT **7.8.1**（FreeCAD 1.1 内置）读同一文件是 **1 面 / 1 wire / v 跨 2π / 完整环面**（1639 点 3122 面，面积 786.66）。同一 IGES 若按 **BRep 模式**写出则是 `198,9,11,10.,2.,13`（解析环面），OCCT 8.0.0 读回 **1 面 / 1 wire / 1369 节点 2592 三角 / bbox ±12**（完整环面）。端口在此条上与规格 8.0.0 **逐点一致**（302/146 到 1e-7），故 Torus 的读侧 sweep 仍可作 parity 检查，但 **302/146 不得判为端口缺陷或"正确输出"**。
- 临时：`.target-gate/`（已 gitignore）存放门禁输出、对拍脚本与探针 dump

---

## 文档索引

- `CLAUDE.md` —— 工作纪律（**改代码前必读**）
- `~/.cursor/projects/d-source-repos-dogs/canvases/occt-code-review-status.canvas.tsx` —— Cursor 自带画布（IDE 内渲染，不在仓库内）：量测、门禁表（`GATES`）、逐模型对拍（`PARITY`）、仍开登记项（`OPEN`）、文档偏移（`DOCDRIFT`）、包时间线（`PACKAGES`）。仓库内自建板 `specs/occt-scope-drift.canvas.tsx`（import 从未存在的 `dsh/canvas`，Cursor 从不渲染）已于 2026-10-16 退役删除
- `docs/goals/drift-and-gap-autofix/state.yaml` —— 长跑自主板的板真相（`active_task` / 包收据 / 门禁读数）
- `specs/_a3n00_gap_analysis.md` —— a3n00 缺口分析：§9.x 逐轮记录**实测数字与其口径**、
  已排除项、以及每条结论对应的 `.cxx` 行号
- `docs/` —— 阶段性审查文档

---

## 现状（如实）

- **已对齐**：STEP 读入/拓扑组装、`ShapeFix_Face::FixMissingSeam`（导入期接回，
  a3n00 面积比 0.8627 → **0.8996**，带倒角法兰 16 面 2 wires → 1 wire/5 边）、
  网格侧 ModelBuilder/Healer/FaceChecker 主流程；逐模型门禁 `step_obj_gates` 4/4（2026-10-09 复跑通过，逐模型比值沿用 §9.633 记录），
  其中 a3n00 0.8996、T0M 0.9987、acs10 0.9846。
- **已知最大缺口**：**没有 `ShapeProcess`/`FixShape` 驱动器**（端口在 reader 里按面直调子步）。
  实测关掉该算子时 OCCT 有 **131/226 面无法网格化**，说明它是承重环节；端口仅缺分派器与参数下发。
- **IGES 写侧**（读侧未做）：`write_iges` 一条路径服务全部模型。全局段按 `IGESData.cxx` 写 IGES 5.3（版本 `11`，精度 `308,15,308,15`，单位 `MM`），序号 `%7.7d` 补零，壳组 402 的 form 为 `1`。初等面在 Faces 模式下是 120 回转面 + 124，不是 192/194/196/198（OCCT 8.0.0 在 Faces 模式下同样如此：`GeomToIGES_GeomSurface.cxx:576-598` 只在 BRep 模式 `myBRepMode && myAnalytic` 成立时才走 198/194/192/196，Faces 模式的对照文件 `_igesdump/oracle/*.iges` 也是 120 回转面）。`iges_check` 已通过：`linkrods`，以及 `data/occ` 的 4 个 `.step` 和 4 个 `.stp`。与 oracle 参数段一致的有 Cube、Sphere、Cone、Extrusion、HoledPlate、rev、Offset；Cone 的目录 form 与全局段对齐后，查看器里是完整圆锥。仍有几何差的是 Cylinder、Torus、screw、Shape、Shape-1、Shape-2，以及 OffsetPlaneHoleEdge 的颜色/名字尾。
- **未结个案**：a3n00 的“螺母斜切面”（端口 F113）仍为 **4 wires / 空面**。
  已定位到：`fix_missing_seam` 在该面上把 28 条边分成了 **6 条 wire / 5 张面**（其中 1 张 2-wire 面 = 6e+7e，
  另 4 张单 wire 面 = 5/4/3/3），而 OCCT 是 **2 条 wire（22+6）/ 1 张面**——**总边数相同（28），这一步不丢边**；
  结果被包成 `Shell(5)`，而 reader 只接受 `Face` ⇒ 整包被丢弃。打包与剪枝语义均已证明与 `.cxx` 同构
  （更正记录见 §9.450：早期文档里的「5 块单 wire 补丁 6/5/4/3/3」是解析错误），分歧指向
  `SplitWires/BreakWires` 的 **wire 分组**。全过程与后续路径见 `specs/_a3n00_gap_analysis.md` §9.370–§9.450 与看板 T-101。
