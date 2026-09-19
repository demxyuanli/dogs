# _board — OCCT→Rust 对齐：未完成任务监控画板

> **载体**：`D:\source\repos\dogs`，master 工作区。
> **基线**：HEAD = `3323c66`（2026-09-14），其上有一批**未提交**改动（2026-09-15 ~ 09-17）。
> **本画板建立**：2026-09-19 23:50 (+08:00)；**最后实测**：2026-09-20（本轮收敛 T-06 / T-21 / T-02，见 §7）。数据均为本机实跑，非引用文档。
> **本文件同时充当 task_plan / findings / progress 三合一**（用户要求"一份画板"）；维护协议见 §8。
> **上游文档**：`specs/_coverage.md`（覆盖矩阵，2026-09-09，已过期）、`specs/_brepmesh_align_review.md`（BRepMesh 对齐根因）、`specs/_tasks.md`（迁移清单，2026-08-29，已过期）。

---

## 0. 目标（north star）

把 **STEP→OBJ 几何/网格管线对齐 OCCT 8.0.0**：先把门禁收敛到"与基线一致或更好"，再按 OCCT `.cxx` 控制流补齐缺口。
**不做**：特例补丁、OCCT 里不存在的规则（面积比 / 长度滤边 / 体积门 / 启发式 unused 判定）、为对齐新写单元测试、为凑数字调参。

## 1. 参考源（已核对）

| 用途 | 路径 | 状态 / 版本 |
|---|---|---|
| **源码引用树**（代码注释 `Foo.cxx:123` 的出处） | `D:\source\OCCT-src` | git detached HEAD = `d3056ef8` = tag **V8_0_0**（"Bump version to 8.0.0"，2026-05-06）；工作树干净；本地仅此一个 tag。模块组布局 `src/{FoundationClasses,ModelingData,ModelingAlgorithms,DataExchange,Visualization,ApplicationFramework,Draw,Deprecated}`，6,310 `.cxx` + 7,072 `.hxx` |
| **参考 OBJ 生成（DRAWEXE）** | `D:\source\occt-8.0.0` | **已构建发行版，不是源码树**（`src/` 只有 DRAW 资源，全树 0 个 `.cxx`）；`inc/` 7,084 头 + `win64\vc14\bin\DRAWEXE.exe` + 74 DLL + `3rdparty-vc14-64`；`inc\Standard_Version.hxx` 自报 `8.0.0` + `OCC_VERSION_DEVELOPMENT "p1-e52c4021"` |
| 导出配方 | `data/_occ_ref_export.tcl` | `ReadStep → XGetOneShape → bounding -noTriangulation -finitePart → incmesh <maxComp*0.001*4> -angular 20 → WriteObj`；**不要**用 `incmesh 0.001 -prs`（会过度细化） |

- **行号引用核对**：`ElCLib.cxx:176-189` = `ElCLib::EllipseValue`，`A2 = MinorRadius * sin(U)`，与 `crates/occt-core/src/elib/clib.rs` 注释中的 `cxx:176-189` 一致 ⇒ 引用树确实指向 `D:\source\OCCT-src`。
- **⚠ 版本口径差异**：行号按 tag `V8_0_0`，参考 OBJ 由 **8.0.0p1-e52c4021** 的 DRAWEXE 生成（该 commit 不在检出里）。遇到行号错位 / 数值细差先怀疑这里。

## 2. 状态快照（2026-09-19 实测）

| 门禁 | 命令 | 当前工作区 | HEAD `3323c66` 基线 | 判定 |
|---|---|---|---|---|
| 编译 | `cargo check --manifest-path crates/occt-topo/Cargo.toml --all-targets` | ✅ exit 0（846 条警告） | ✅ | 绿 |
| STEP→OBJ bbox parity | `cargo test --manifest-path crates/occt-topo/Cargo.toml --test step_obj_parity` | ✅ **14/14**（含 ATU01038 与 `data/occ` 4 个新模型，2026-09-20） | ✅ 12/12 | 绿 |
| STEP→OBJ 端到端 | `--test step_to_obj` | ✅ 13/13（83.8s） | ✅ | 绿 |
| 面积对拍 | `--test step_obj_area` | ✅ 11/11 | ✅ | 绿 |
| 几何一致性 | `--test step_geometry_parity` | ❌ 2/3（T-05） | ❌ 同 | 红（遗留） |
| topo 单测 | `--lib` | ❌ **1,293 通过 / 1 失败**（2026-09-20：T-02/T-06 转绿） | ❌ **1,235 / 59 失败** | 剩余仅 T-01 |
| boss 合并 | `--test bop_builder2_boss` | ❌ 1/2（T-03） | ❌ 0/2 | 红（遗留，改善） |
| phase19 | `--test phase19_integration` | ❌ 3/5（T-04） | ❌ 3/5（完全一致） | 红（本轮未触碰） |
| phase20 / phase3 | `--test phase20_integration` / `phase3_integration` | ✅ 5/5 · 4/4 | ❌ 1/5 · 3/4 | 本轮修绿 |
| phase4–10 | 各 `--test phaseN_integration` | ✅ 全绿 | ✅ | 绿 |
| core / math / geom2d | 各 crate `--lib` | ✅ 293 / 215 / 72 | — | 绿 |
| geom | `cargo test --manifest-path crates/occt-geom/Cargo.toml --lib` | ✅ **153/153**（T-06 已修，2026-09-20） | ✅ 153/153 | 绿（曾 152/153） |

**HEAD 基线复现方法（不改共享树）**：
```powershell
cd D:\source\repos\dogs
git worktree add --detach .target-headcheck HEAD      # 路径命中 .gitignore 的 /.target-*/
cd .target-headcheck
cargo test --manifest-path crates/occt-topo/Cargo.toml --no-fail-fast
cd ..; git worktree remove --force .target-headcheck
```
> 历史事故：**禁止**在共享树里 `git stash` / `stash pop`（曾造成 44 处冲突标记）。要 A/B 就用上面的 worktree。

## 3. 任务板

状态枚举：`pending` / `in_progress` / `blocked` / `done` / `waived`（豁免必须写依据与引用）。

### P0 — 红门禁（收敛或显式豁免）

| ID | 任务 | 证据（文件:行 / 实测） | 期望 vs 实际 | 下一步 | 验收 | 状态 |
|---|---|---|---|---|---|---|
| T-01 | `brepfeat::tests::groove_cuts_cylinder` 体积不符 | **根因链（探针实测，2026-09-20，全部已复现）**：① 工具正确（16 边形环，体积 3.2574 = 解析值）；② 两输入全平面 ⇒ `bop_builder::boolean` → `bop_builder2::builder_bop_with_fuzzy` → `perform_internal`（filler → fill_images → build_result → build_shape）；③ `build_result`：`gf_solids=2 gf_faces=256`（装配在跑）；④ `build_shape`：`op=Cut dim0=3 dim1=3 open=false` → `build_rc`；⑤ **`build_rc`：`obj_src=1 tool_src=1 it_shapes=1 check_keys=1 it_exp=1 tool_sets=1` ⇒ 每个 solid 只有 1 个镜像 = 没被切开** ⇒ BuildRC 只能整体取舍 ⇒ Cut 几乎不切、Common/Fuse 空；⑥ `fill_images_solids` 有在跑（`fill_in3d in_parts=1`，两个 argument 都有镜像），但 `build_split_solids_occt` 每个 solid 只回 1 个镜像 | 分派已对上 OCCT（`BOPAlgo_BOP::BuildShape` `:871` = `BuildRC` `:900` + 仅 FUSE 且 dim0==3 走 `BuildSolid` `:902-904`；Rust `bop_bop.rs:548/171/277`），故不一致在**面级切割**：柱面 24 个侧片应在 z=1/z=1.8 与 r=1 处被工具面切开却未切 | **round2 追加（探针实测）**：PaveFiller 侧无问题——`interf_ff=422`、`ff_curves=144`、`ff_pb=144`、`ff_pb_with_edge=144`（截面边真建了）；面级绑定也无问题——224 源面里 `faces_sc=112 / faces_on=112 / faces_in=0`。瓶颈在**切割重建**：`BuildSplitFaces` 建了 112 个任务、`empty_le=0`（都收到分割边集），但 `rebuild_split_areas` 块数分布 = **0块:40 / 1块:41 / 2块:30 / 3块:1** ⇒ 只有 31 面真出多块，40 面镜像列表为空 | **round3 追加（探针实测）**：112 个任务面**全是平面**；`perform_shapes_to_avoid` 对所有面都**没触发**（`avoid_sta=0`），`avoid` 由 `perform_loops` 的"未进环边"后处理填充。按 `(in_edges, avoid, avoid_bnd, loops, areas)` 分组：`avoid=0` 的组都是正常切割（如 in_edges=11→loops=3/areas=2 ×15；in_edges=15→loops=5/areas=1 ×16；in_edges=9→loops=2/areas=2 ×7）；**45 个坏面一律 `avoid=3` 且 `avoid_bnd=1`**（= 1 条原边界碎片 + 2 条截面边都没能链进环），其中 30 面最终 `areas=0`（唯一成环的是被切成 hole 的那条 → 无 growth → 按 `builder_face_occt.rs:389-408` 只对**开放面**兜底，闭合面直接丢弃）、15 面 `areas=1` | **round4 追加（探针实测）**：① **推翻**"同点不同 vertex"假设——112 个面全部 `dup_pos_groups=0`（66 面 avoid=0、46 面 avoid=3，都无坐标重复而 key 不同的顶点）；② 解码首个坏面：它是 `mesh_cylinder(1,3,24)` 造出的**柱面侧面三角形**（顶点 `(1,0,0) / (1,0,3) / (0.9659,-0.2588,3)`，四边形的对角剖分），被工具的 z=1 与 z=1.8 两个平面切成 **3 段**；其 9 条唯一边**全部在输入集里**（竖边分 3 段、顶边 1 条、斜边分 3 段、两条水平弦边**各出现 2 份**）；③ 结果：`in_loop=true` 的只有竖边 z0→1、z1→1.8、斜边 z1.8→1.0、z1.0→0 与两条弦边 ⇒ **下面两段闭合成功，最上面那段（竖边 z1.8→3 + 顶边 + 斜边 z3→1.8 + 一份弦边）4 条边齐全却没成环**，3 条全部落进 `avoid`（其中 1 条是原边界边、2 条是边界边的分割像，与 `avoid_bnd=1` 吻合） | **round5 追加（探针 + OCCT 对读）**：① 逐行对读确认 `SplitBlock`（OCCT `_1.cxx:112-354`）与 `Path`（`:358-617`）**与 Rust 移植一致**（含 `iCnt==1` 早退、同边反向 `aTwoPI`、`aNbWaysInside==1` 覆盖、切环/`iPriz` 逻辑）；② 再排除两个假设：`make_2d` 失败（62 次成功调用 `skipped_no_pcurve=0`）、键合并（`GeometryRegistry::shape_key` = 指针身份）；③ 45 个面只走出**部分环**（15 面：7 边→1 环 3 边；22 面：11 边→1 环 5 边；8 面：13 边→1 环 6 边）——行走在那些顶点卡住；④ **命中根因**：`wire_splitter_block.rs:541 refine_angles` 是**空实现**，只保留了 OCCT 的 `iCntBnd != 2` 早退，对 `iCntBnd == 2` **什么都不做**；而 OCCT `RefineAngles(vertex)`（`BOPAlgo_WireSplitter_1.cxx:925-1032`）恰在此情形用 `RefineAngle2D`（`:1033-1125`）重算内部出边角度，并在 `iCntInt == 2` 时按 `Precision::Angular()` 微调（`:998`）。坏面每个弦端点正是"2 条边界边 + 2 份内部弦副本"= `iCntBnd==2 && iCntInt==2`，完全落在 stub 跳过的那一支 | **round8 结果（根因 + 一次被回退的实验）**：① OCCT `BuildSplitFaces` 取边规则解码（`_2.cxx:412-494`）：退化→`Orientation(anOriE)`；INTERNAL/closed→Fwd+Rev 两份；常规→`Orientation(anOriE)`+`IsSplitToReverseWithWarn`；In/Section paves→Fwd+Rev 两份。端口**完全一致**；实测 **544 个边界碎片、朝向不匹配 0** ⇒ 朝向施加是对的（另发现 `is_split_to_reverse_edge:727` 用 edge-TShape 捷径而 OCCT 比 **curve 句柄**（`BOPTools_AlgoTools.cxx:1462`）⇒ 新开 T-31）；② **T-01 真正根因**：坏面的 wire 由 `mesh_to_brep.rs:102 make_wire(&[e_ab,e_bc,e_ca])` 造出，而 `ensure_edge` 的**共享边缓存**把边按"先请求的三角形"的方向存，**后者复用时未修正朝向** ⇒ 实测该面三条边界边全 Forward 且对角线方向属于邻接三角形 ⇒ wire 不是链接链（两条边从 `(1,0,3)` 出发、两条汇入 `(0.966,0.259,0)`）⇒ splitter 块出现 source/sink ⇒ `Path` 死路 ⇒ 面不切 ⇒ groove 切不动。对照 OCCT：`BRepBuilderAPI_MakeWire::Add` 会做 `CorrectEdgeOrientation`；③ **实验（已回退）**：仅在 `mesh_to_brep.rs:102` 加"按请求方向翻转共享边"后，lib 变 **1292/2**——`groove_cuts_cylinder` 前移到 `:98`"必须减材"断言（after≈before），且 `groove_negative_volume_delta` 新红 ⇒ **单独修 wire 朝向不够**：OCCT 的 `BRepBuilderAPI_MakeFace` 会按"wire 绕向 vs 曲面法向"共同决定**面朝向**，必须与 wire 修正同一步做 | **下一轮（round9）**：在 `triangulation_to_brep` 里按 OCCT 一并移植 **`BRepBuilderAPI_MakeWire::Add`（CorrectEdgeOrientation）+ `BRepBuilderAPI_MakeFace` 的面朝向规则**（并检查 shell/solid 装配的朝向），再复跑 T-01 与全量门禁 | `--lib brepfeat::` + 全量 | in_progress（根因已定位到 mesh→BRep 的 wire/face 朝向前提；实验证明需成套修） |
| T-02 | `builder_solid::tests::inner_cavity_is_absorbed_as_a_hole` 内腔未吸收 | `crates/occt-topo/src/builder_solid.rs:723`（生产代码**未改**，只改测试夹具） | 探针实测：`n=2 faces_per_loop=[6,6] growth=[0,1] holes=[]` ⇒ 内壳被判 growth，故 2 个 area | **已修（夹具不具代表性）**：`IsHole` = 无穷点分类（`BOPAlgo_BuilderSolid.cxx:823-831`）+ 分类器读朝向（`brep_class3d.rs:326`）⇒ 只有**朝向翻转**的壳才是 hole。`BRepPrimBox::make_box_corner` 造的是正朝向盒子（=第二个 growth），旧夹具因此断言了非 OCCT 结果。夹具改为按真实 Cut 分裂的腔壳反转内壳朝向后：`growth=[0] holes=[1]` → 1 area / 12 面 | `cargo test --manifest-path crates/occt-topo/Cargo.toml --lib builder_solid::` → 全绿 | **done**（2026-09-20） |
| T-03 | `boss_single_disc_base_merges_one_solid` boss fuse 未合并成 1 solid | `crates/occt-topo/tests/bop_builder2_boss.rs:114` | `left=0 right=1`（HEAD：两条用例全红） | 已知阻塞点：`build_split_solids_full` 未跨 solid 合并共面片（remember 记录 B3 实体装配缺口） | `--test bop_builder2_boss` | pending |
| T-04 | phase19 两条 PaveFiller 管道用例 | `tests/phase19_integration.rs:40`（got 84）、`:53`（got 88） | 期望"至少包含两个原始 box 之外的更多形状"；HEAD **完全一致** | 本轮未触碰；属 PaveFiller 阶段遗留，单独开波次 | `--test phase19_integration` | pending |
| T-05 | `offset_geometry_is_consistent` 偏差 | `tests/step_geometry_parity.rs:129` | `Offset divergence 2208.0 vs shape_volume 1612.9`；HEAD **同** | 单独定位（offset 体积 vs 面类型分解）；与 T-04 同属遗留。**2026-09-20 全仓审查提示：很可能与 A0（offset 公式写错，见 T-35）同源**——先修 T-35 再看此偏差是否变化 | `--test step_geometry_parity` | pending |
| T-06 | `extrema_pc::tests::point_ellipse_y_axis`（**唯一新增红**） | `crates/occt-geom/src/extrema_pc/tests.rs:57-64` | 期望 `u=3π/2`，实得 `π/2`；同用例的 `distance=3.0`、`p2=(0,1,0)` 断言**均通过** | **已修**：期望改 `π/2` 并注明 `ElCLib::EllipseValue`（`cxx:176-189`，`+Minor*sin(U)`）；**未**回退 `clib.rs` 符号（那是 ATU01038 面 130/140/156/216/222 的修复） | `cargo test --manifest-path crates/occt-geom/Cargo.toml --lib` → **153/153 ✅** | **done**（2026-09-20） |

| **T-35** | **A0：offset 曲线公式写错（全仓审查最高危）** | `crates/occt-geom/src/offset.rs:14-19`（3D 用 `p + Offset·Dir` **平移**）↔ `Geom_OffsetCurveUtils.pxx:53-61`（应沿**法向** `p + Offset·(D1×Dir)/‖D1×Dir‖`）；`crates/occt-geom2d/src/offset.rs:19-41` 法向 `(-dy, dx)` ↔ `Geom2d_OffsetCurveUtils.pxx:50` `(dy, -dx)`；两侧 `d1/d2` 均缺 `DNdir` 项（**已双向对读核实**） | 产出错误几何（**公式错误，非近似**）；3D 版被 STEP 读入 `step/p04.rs:1049` 直接使用，2D 版被 `geom_bnd_lib_offset2d.rs:101` 消费 | 按两个 `OffsetCurveUtils.pxx` 重写 D0/D1/D2（含 `DNdir` 旋转项与失败返回），并复核上述两处调用点 | 现有 offset 相关门禁不回归（`--lib` + `--test step_geometry_parity`） | pending |
| **T-36** | **全仓忠实度审查 A1–A26 分批整改**（A0 已单列 T-35） | `specs/_audit/_index.md` §3/§7（5 区报告齐；合计自创 **71** / 未登记 **17** / 已登记 **15**，系统性条目 A0–A26） | 按 §7 分批：A6→A7 → A1 → A2+A11 → A13/A23 → A3/A24 → A20/A21/A17 → A18/A19/A22/A25 → A9/A10/A14/A26 → A12 | 逐批对着 `.cxx` 改控制流，**不加 OCCT 之外的规则/阈值**；无对应分支的标 `UNPORTED` + OCCT 文件行号 | 每批跑该 crate `--lib` + 相关门禁（parity 14/14 等）不回归 | pending |

### P1 — 未提交 wave 的收尾

| ID | 任务 | 证据 | 动作 | 状态 |
|---|---|---|---|---|
| T-07 | 178 文件未提交（148 M / 30 D / 32 ??；排除参考 OBJ 后 177 文件 +10,053/−2,540） | `git status --short` | 分批提交：① 源码修复（clib/primitives/shape_ops/pcurve_full/shhealing/bnd 等）② 新增单元（13,900 行 + `plib_jacobi_coeffs.pxx` 5,585 行）③ 数据与参考 ④ 清理 | pending |
| T-08 | 30 个 `specs/_phase*_plan.md` / `_loop*_plan.md` / `_rules.md` / `_brepmesh_migration.md` / `_brepmesh_align_tasks.md` 与 `examples/demo.rs` 被删（`Cargo.toml` 同步去掉 `[[example]] demo`） | `git status --short specs examples` | 确认是有意清理；`_brepmesh_align_tasks.md` 里"残留差距"清单要先合并进 `_board.md`（§3 P2 已吸收）再删，避免丢任务 | pending |
| T-09 | `output/` 未忽略（27 文件 **16.35 MB**）且 `data/occ/`、`data/ATU01038.step` 未跟踪 | `git status --short` | 决策：入库存档 or 加 `.gitignore`（当前 `data/output/` 已忽略，但根 `output/` 没有） | pending（需用户拍板） |
| T-10 | `data/occ-ATU01038.obj`（**已跟踪**）被原地重导覆盖 | 相对 HEAD 绕 X 轴 90°；顶点 17,720 → 18,102 | 提交时说明"参考按新输入 STEP 重导"，否则后人无法解释 diff | pending |
| T-11 | `occt-topo` lib 846 条编译警告 | `cargo check` | 提交前至少清 unused import/variable（不动行为），否则真信号被淹没 | pending (低优先) |

### P2 — 新模块自报的 PARK / UNPORTED（按 OCCT 控制流补齐）

| ID | 位置 | 未移植内容 | 状态 |
|---|---|---|---|
| T-12 | `crates/occt-topo/src/geom_bnd_lib_curve3d.rs:7` | Ellipse / Hyperbola / Parabola 解析盒（`occt_geom::Curve` 不暴露 `gp_Elips/Hypr/Parab`）→ 现走采样；周期 B 样条 arm（缺 `Segment`/`AdjustPeriodic`，`BSplineCurve.cxx:44-49,299-330`） | pending |
| T-13 | `geom_bnd_lib_surface3d.rs:8` | `BoxOptimal`（`OptimizationHelpers.pxx` PSO+Powell）+ `SurfaceOfExtrusion::BoxOptimal`（`cxx:224-283`）；`catch(Standard_Failure)` 回退 | pending（`AddOptimal` 路径，当前不在 `BRepBndLib::Add` 上） |
| T-14 | `brep_bnd_lib.rs:13` | 三角化 arm（`BRepBndLib.cxx:95-101,157-179`）、仅 pcurve 边的 `BRepAdaptor_Curve::Initialize` 回退（`cxx:93-106`） | pending（缺 `Poly_Triangulation` 存储） |
| T-15 | `shhealing/transfer_params.rs:21` | `CopyNMVertex`（需 `BRep_PointRepresentation`）、`CorrectParameter` knot snap（`Proj.cxx:268-279`）、`myLocation` 统一 identity | pending |
| T-16 | `pcurve_full/p03.rs:27`、`p04.rs:11` | iso 链在 iso 边界边上的臂选择；`ShapeAnalysis_Surface::myGap` 用残差恢复；奇点数组未缓存（结果等价）；`ProjectDegenerated` 单点重载无调用者 | pending（已写明等价性理由） |
| T-17 | `shhealing/shape_analysis_curve.rs:24` | `project_act` 的 `!OK` 分支缺 Hyperbola/Parabola（`cxx:376-386`）与 Ellipse（`cxx:394-399`）快捷路径 → 走 default 段搜索 | pending |
| T-18 | `crates/occt-core/src/intf/mod.rs:6` | `Intf_Tool`、`Intf_InterferencePolygonPolyhedron.gxx`（3D 多面体干涉，无 2D 调用者） | pending |
| T-30 | `uv_tolerance_2d`（`wire_splitter_block.rs:507`）用 `t3/|d1u|` 近似，OCCT `BOPAlgo_WireSplitter::Tolerance2D`（`_1.cxx:859-881`）用 `GeomAdaptor_Surface::U/VResolution`（平面恒 1.0，B样条面 ×1.1） | 平面面下 OCCT `Tolerance2D≈1.0` vs Rust `≈1e-7` ⇒ `Angle2D` 的 `dt` 差数量级（直线 pcurve 无影响，曲面面会影响取角） | pending（round7 发现，非 T-01 阻塞） |
| T-31 | `is_split_to_reverse_edge`（`algo_tools_face.rs:727`）的捷径比较 **edge TShape 身份**；OCCT `BOPTools_AlgoTools::IsSplitToReverse`（`BOPTools_AlgoTools.cxx:1461-1465`）比较 **curve 句柄**（`aCSp == aCOr`） | **已修（round9）**：捷径改为 `Arc::ptr_eq(edge_curve(sp), edge_curve(or))`，与 OCCT 一致（碎片与父边共享 curve ⇒ 只比朝向，不再落几何切线分支）。实测：全量 `--no-fail-fast` **与基线完全一致（无回归）**；T-01 数值仅末位变化（8.903438339491375→…373）⇒ 本例两种判定同解 | `--lib` + 全量 | **done**（2026-09-20） |
| T-32 | `mesh_to_brep.rs:102` 的 wire 未做 `BRepBuilderAPI_MakeWire::Add` 的 `CorrectEdgeOrientation`（共享边由邻接三角形先建、按彼方向缓存） | 直接后果：mesh→BRep 的三角形 wire **不成链**（实测该面三边全 Forward、对角线方向属邻接三角形）⇒ wire splitter 块现 source/sink、`Path` 死路。**单独修 wire 朝向会回归 2 个用例**（round8 实验：lib 1292/2）；round9 度量补充：mesh→BRep 柱面 **47/96 个面的"面法向·wire 绕向"为负**（绕向约定本身混乱），而旋转体工具 **112/128 一致为负** ⇒ 修复必须与 `BRepBuilderAPI_MakeFace` 的面朝向/绕向约定**成套**做 | **P0 关联**（T-01 根因；round10+ 实施） |
| T-19 | `meshing/triangulator.rs:19` | `addTriange34` / `checkCondition` 快路径（仅 `Perform` 调用，未移植）；注：`BRepMesh_Triangulator` 只被 VRML 读路径引用，**不属** STEP/OBJ 网格管线 | pending（低优先） |
| T-20 | `bop_build_solids_leftover.rs:8` | `merge_sharing_faces` 已删（OCCT `BOPAlgo_Builder_3.cxx:579-616` 无此 stage）；leftover 模块去留待定 | pending |

### P3 — 对齐覆盖扩展（新模型）

现状（2026-09-19 实测）：

| 模型 | STEP 源 | OCCT 参考 | Rust 产物 | bbox 一致 |
|---|---|---|---|---|
| ATU01038 | `data/ATU01038.step`（未跟踪，19,637 行） | `data/occ-ATU01038.obj`（1,8102 v）· `output/occ/ATU01038.obj`（18,051 v） | `output/ATU01038.obj`（17,767 v / 22,160 f） | ✅ Δ≤9e-6 |
| a3n00 / acs10 / bottom / motoc / top / T0M / TDB | `data/occ/*.step`（8 个） | `output/occ/*.obj`（8 个，9/15 20:56–22:45） | ❌ 无 | — |
| Cube/Cone/Cylinder/Sphere/Torus/HoledPlate/Offset*/Shape*/rev/linkrods/screw/Extrusion | `data/*.step` | `data/occ-*.obj` / `occ-*.obj` | `output/*.obj`（9/17 16:40） | 12/12 门禁绿 |

| ID | 任务 | 说明 | 状态 |
|---|---|---|---|
| T-21 | ATU01038 纳入 parity 门禁 | `tests/step_obj_parity.rs` 原是 **12 条硬编码用例**，无 ATU01038。**已加**一条 `check_parity("ATU01038","occ-ATU01038.obj",1e-4)`（实测 Δ≤1.1e-5；密度不断言）；加用例 = 引用现有基线，符合规则 3 | `cargo test --manifest-path crates/occt-topo/Cargo.toml --test step_obj_parity` → **13/13 ✅** | **done**（2026-09-20） |
| T-22 | 7+1 个 `data/occ` 新模型跑 Rust 导出并对拍 | **已跑 + 已锁门禁（round10/11）**：`export_data_obj` 增加可选目录参数（`-- data/occ`，无参数时行为不变），8 个模型全部导出成功。**bbox 对拍（Rust 网格 vs `output/occ/*.obj`）**：`ATU01038` **0/0/0**、`bottom` **0/0/0**、`motoc` **0/0/0**、`top` **0/0/0**、`T0M` **0/0/0**（密度 −8.8%~+0.8%）→ 这 4 个（+ATU01038 已有用例）**已锁进 `step_obj_parity`**（新用例 `occ_test_model_bboxes_match_occt`，参考入库 `data/occ-ref/`，门禁 **14/14**）；**3 个仍有真实几何差**：`a3n00` Δy=**60.04**、`acs10` Δy=**7.88**、`TDB` Δx=**0.59** | **round12 追加（归因到网格侧）**：`a3n00` 的网格 y-max = **137.539** vs OCCT **77.5**（+60.039），其余轴全对（x 相等、y-min 相等、z 差 0.006）；而 **226 个面的 wire 顶点全部在 |y|≤77.6 内（0 个越界）** ⇒ **解析边界与 OCCT 一致，溢出由网格生成产生**（面被按超出 trim 的范围采样，即 `brepfeat/p01.rs` 已注明的"按 UV 窗口网格化"那一类）。逐面单独网格化给出候选越界面（41–49、77–87 等 `Plane`/`Other`，自身网格超出 wire y 约 10；**注意**：单面网格化可能与实体导出走不同路径，此列表为指示性）。下一步配方：在导出路径内部 dump **每面网格 bbox** 与该面 wire bbox 对比，定位后用 OCCT `BRepMesh_FaceDiscret`/`BRepMesh_ModelBuilder` 追该面的 trim/UV 窗口逻辑 | 导出 + bbox 对拍 + parity | **done（导出扩展 + 5 模型锁门禁）**；3 个模型的**网格侧**缺口 pending |
| T-23 | 密度差是否要收 | ATU01038 顶点 −1.85%（vs `data/occ-ATU01038.obj`）/ −1.57%（vs `output/occ`），面 −1.89% / −1.21%。parity 注释明确"不断言密度，UV-grid vs deflection-adaptive 是已知差" ⇒ 若收，需先证明 OCCT 侧公式逐项对上 | pending（待定性） |
| T-24 | `output/` 与 `data/occ/` 的归档策略 | 见 T-09 | pending |

### P4 — 继承缺口（`specs/_coverage.md` 列出，本轮未处理）

| ID | 任务 | 状态 |
|---|---|---|
| T-25 | `GeometryRegistry` 侧表（边/面几何不在 `TShape` 内）——架构级重构 | pending |
| T-26 | 无面 operand 仍走 `voxel_fallback`（绿路径 `bop_builder2`） | pending |
| T-27 | `GetFaceOff` / 角法向：leftover `is_covering_face` 死路径未接 `BopBuilder` | pending |
| T-28 | ImpPrm HVertex 合并（`IntPatch_ImpPrmIntersection` cxx 329–465）未移植；`intana` 与 `intpatch` 重叠未合并 | pending |
| T-29 | `specs/_coverage.md` 已过期（2026-09-09；行数/测试数与本轮不符）→ 本轮收敛后刷新 | pending |

### P5 — 自创 → 忠实移植（审查 A0–A26 逐项，2026-09-20 立项）

> 来源：`specs/_audit/_index.md`（5 区报告 + A0–A26）。**父任务 T-36**，下列每项一个可跟踪 ID。
> 验收统一为：① `cargo check` 过；② 该 crate `--lib` 不低于基线（topo 1293/1、geom 153、geom2d 72、core 293、math 215）；③ 相关门禁不回归（`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3）；④ 无对应 OCCT 分支时标 `UNPORTED` + OCCT 文件行号，**不加新规则/阈值/启发式**。

| ID | A# | 位置（文件:行） | 自创内容 → OCCT 对应 | 阶段 | 状态 |
|---|---|---|---|---|---|
| T-35 | **A0** | `occt-geom/src/offset.rs:14-19`、`occt-geom2d/src/offset.rs:19-41` | 3D 沿参考方向平移 / 2D 法向反号、缺 `DNdir` → `Geom_OffsetCurveUtils.pxx:47-115`、`Geom2d_OffsetCurveUtils.pxx:43-` | **1（先修）** | **done**（2026-09-20） |
| T-63 | A0 派生 | `occt-geom/src/curve.rs:8-11`（`Curve::d3`）、`occt-geom2d/src/curve.rs`（`Curve2d::d3`）、各具体曲线 | `EvalD3`（`Geom_Curve::EvalD3` / `Geom2d_Curve::EvalD3`）无任何具体曲线实现 ⇒ offset 的 `CalculateD2` 里 `D2Ndir` 项恒为 0；`EvalD2` 的 `AdjustDerivative` 奇异支路（`Geom_OffsetCurve.cxx:311-330`、`Geom2d_OffsetCurve.cxx:265-280`）未移植（`isDirectionChange` 恒 false） | 1 后 | pending |
| T-64 | **A27**（新） | `occt-geom/src/{line,circle,ellipse,hyperbola,parabola,plane,cylinder,cone,sphere,torus}.rs`、`occt-geom2d/src/{line,circle,ellipse,hyperbola,parabola}.rs` 的 `continuity()` | 解析曲线/曲面一律返回 `3`（G2），OCCT `Geom_Conic::Continuity()`=`GeomAbs_CN`（`Geom_Conic.cxx:32-35`）、`Geom_Line`(`:128-131`)/`Geom_ElementarySurface`(`:25`)/`Geom2d_Line`(`:172`)/`Geom2d_Conic`(`:48`) 同为 CN(6)；端口 B 样条侧 `bspl::local_continuity` 单跨返回 6，证明编码约定一致 ⇒ 解析类少报连续性，影响 `NbIntervals`/`parameter_intervals` 消费方 | **已尝试 → 回退**（2026-09-20）：仅把这 15 处改成 6 后，`step_to_obj` 由 13/13 变 12/13——**Sphere 网格变空**（`tests/step_to_obj.rs:240` round-trip 顶点数 0）。⇒ 不能单独改：某个消费方把"报出的连续性"当跨度上限（嫌疑点 `meshing/range_splitter/p01.rs:295` 的 `continuity <= curve.continuity()`、`meshing/edge_discret.rs:1259` 的 `pcurve.continuity().min(surface.continuity())`）。下一步：先按 OCCT 对齐该消费方（`GeomAdaptor_Surface::NbUIntervals`/`NbVIntervals` 语义），再改这 15 处 | pending（被消费方阻塞） |
| T-37 | A1 | `occt-topo/src/brep_surface.rs:234-272`（39 处/27 文件） | 网格扫描 + 6 轮二分 → **需先移植** `Extrema_ExtPS`/`Extrema_GenExtPS`/`Extrema_ExtPElS`（仓内 `extrema_surf` **不是**忠实件：24×24 网格 + 数值 Jacobian，见 T-67） | 3 | pending（被 T-67 阻塞） |
| T-67 | A1 前置 | 新增移植：`occt-geom/src/extrema_surf/`（或新模块） | 移植 `Extrema_ExtPS.cxx`(376)+`.hxx`(140)、`Extrema_GenExtPS.cxx`(1056)+`.hxx`(155)、`Extrema_ExtPElS.cxx`(454)+`.hxx`(76)：类型分派 + iso 退化处理（`IsoIsDeg`）+ 逐 C2 区间采样（`mySample`）+ `math_FunctionSetRoot` 解析系统（`Extrema_GFuncExtPS`），并改用 `Surface::d2`（`surface.rs:12`）的解析 Jacobian 取代数值差分 | **3 前置** | pending |
| T-38 | A2 | `occt-topo/src/brep_extrema.rs:252-301` | 7×7 采样 + 射线奇偶 → `SolidClassifier`/`algo_tools::compute_state`（`BOPAlgo_BuilderSolid.cxx:835-860`） | 4 | pending |
| T-39 | A3 | `occt-topo/src/step/p01.rs:114-144,663-680` | 6 点二阶差分猜曲线族 → `StepToGeom.cxx:1335-1349` 按实体类型 / `GeomToStep_MakeCurve.cxx:54` | 5 | pending |
| T-40 | A4 | `occt-topo/src/algo_tools/p01.rs:165-205` | V/E/F 分支 32×32 投影 → `BOPTools_AlgoTools::ComputeState` 精确投影 | 4 | pending |
| T-41 | A5 | `bop_builder_core.rs:79-96`、`bop_curved/p02.rs:409-581`、`p04.rs:42-70,325-401` | 体素/网格布尔与计票 → 无 OCCT 对应 ⇒ **摘除并标未移植** | 4 | pending |
| T-42 | A6 | `occt-core/src/gcpnts.rs:30,87-225` | 中点二分 `MAX_DEPTH=16` 伪造 deflection → 调 `gcpnts_perform.rs:52`（已逐行移植） | **2** | **done**（2026-09-20） |
| T-43 | A7 | `occt-geom/src/extrema_pc/p01.rs:563,621`、`extrema_cc/p02.rs:83,214-239` | `clamp(24,256)` 网格 + 16 兜底 → 调 `extrema_pc/p03.rs`（`Extrema_GGExtPC`） | **2** | **部分 done**（2026-09-20）：点–曲线侧已忠实；`extrema_cc` 侧见 **T-66** |
| T-65 | A7 派生 | `occt-geom/src/extrema_pc/p01.rs`（`ellipse_all`/`hyperbola_all`/`parabola_all`）、`p02.rs::ext_pelc_all` | 三条解析臂已存在但未接线：缺 `myIsMin` 端口（`Extrema_ExtPElC.cxx:277`、`:381`、`:473`）⇒ 椭圆/双曲/抛物目前走 `default:` 臂；且 `Curve` 无 hyperbola/parabola 类型查询（T-12） | 2 后 | pending |
| T-66 | A7 派生 | `occt-geom/src/extrema_cc/p02.rs:40-115` | 通用曲线–曲线**种子集**自创（均匀网格 + 局部极值 + 边界最优点）；OCCT `Extrema_GenExtCC::Perform` 用 `math_GlobOptMin`（仓内 `occt-math/globoptmin.rs` 已移植）+ `Extrema_ECC` 的 `math_FunctionSetRoot` 起点 | 2 后 | pending |
| T-44 | A8 | `occt-geom/src/convert_bspl.rs:326,349` | 采样折线 + 强制 1 次 → `GeomConvert_CompCurveToBSplineCurve.cxx:135-215` | 8 | pending |
| T-45 | A9 | `convert/`、`cslib/mod.rs:24-41`、`gprop/mod.rs:27-42`、`bnd/obb_pca.rs`、`bnd/intersect.rs`、`geom/polyline_simplify.rs`、`int/curve_curve.rs` | 假出处（写 OCCT 包名但该包无此函数）→ 改标真实出处/非 OCCT（`polyline_simplify` 经 `intpatch_trace.rs:260` 改变交线几何，须标注） | 8 | pending |
| T-46 | A10 | `elib/surface_eval.rs:100,104`、`bnd/bsphere.rs:19-29`、`elib/intersect.rs:58-60`、`poly/make_loops.rs:234` | 静默默认值/凭空造值/只取首候选 → `ElSLib::SphereD2/TorusD2`、`Bnd_Sphere.cxx:73-96`、`Poly_MakeLoops.cxx:611-700` | 8 | pending |
| T-47 | A11 | `bop_builder_dispatch.rs:473-601`、`bop_draft_solid_occt.rs:41,49`、`edge_edge/p01.rs:368-371` | 质心规则 / `solid−face` 原样返回 / 丢"顶点<3"的面 / 采样解当补集 → `BOPAlgo_Builder_3.cxx:329` 等对应控制流 | 4 | pending |
| T-48 | A12 | `brepfeat/p01.rs:109,419,361,508`、`feature.rs:197` | 体积解析覆盖 / 网格夹具冒充 `BRepPrimAPI_MakeCylinder` / `clamp(16,64)` → `BRepFeat_MakeDPrism/MakeRevol`、`LocOpe_Revol` | 8 | pending |
| T-49 | A13 | `occt-topo/src/wireframe.rs:392-415` | 未裁剪 UV 窗口规则网格 → `BRepMesh_FaceDiscret` 按 pcurve 边界离散 | 5 | pending |
| T-50 | A14 | `geom/`（csg/delaunay/triangulate/fit*/polygon_*）、`elib/measure.rs`、`validate.rs`、`hlr.rs`、`viz_scene/`、`draw/`、`xcaf/`、`render_svg.rs` | 整包非移植件（`validate.rs` 的 `is_valid()` 恒 true 风险最高）→ 模块头声明未移植 + 真实出处 | 8 | pending |
| T-51 | A15 | `extrema_cc/p02.rs:210-241`、`hyperbola.rs:16`、`parabola.rs:16`、`surface.rs:12-20`、`intana/p02.rs:574`、`extrema_ss.rs:69,119`、`curve_reparam.rs:37,189,228`、`gcpnts.rs` 积分 | 失败被 `fallback_*` 吞掉（丢 `StdFail_NotDone` 语义）+ 自创阈值 → 逐条按 OCCT 返回失败或补齐精确解 | 8 | pending |
| T-52 | A16 | `occt-geom/src/geom_api.rs:52,97`、`occt-geom2d/src/curve_ops.rs:68-69` | 256×256 采样求交/投影 → `Extrema_ExtPS/ExtCC` + `intana2d`/`intimpargen` | 3 | pending |
| T-53 | A17 | `occt-topo/src/brep_exchange.rs:56-98` | 偏转形参被 Prs3d 顶掉 + 三级静默回退 → `RWObj_CafWriter`/`RWMesh_FaceIterator.cxx:87-89`（不网格化、空则跳过）⇒ 删回退、形参改名 | **6** | pending |
| T-54 | A18 | `occt-topo/src/wireframe.rs:257-380` | 平面耳切 + 质心角度排序 + 桥洞 → 约束 Delaunay（`BRepMesh_DelaunayBaseMeshAlgo` + `BRepMesh_Delaun`） | 7 | pending |
| T-55 | A19 | `meshing/incremental_mesh/p01.rs:369,460-467`、`wireframe.rs:407-408`、`brepmesh.rs:38,174-175,259` | `WIREFRAME_FALLBACK_RATIO_MAX=0.10` 失败率换算法 + 四叉树魔数 + `clamp(3,64)` → OCCT 无失败率阈值，逐面置 `IMeshData_Failure`（`BRepMesh_BaseMeshAlgo.cxx:52-62`） | 7 | pending |
| T-56 | A20 | `occt-topo/src/step/p05.rs:508-604` | `p1.distance(p2) < 1e-3` 替代 `V1.IsSame(V2)`；`Other`/HYPERBOLA 边域落 `(0,1)` → `StepToTopoDS_TranslateEdge.cxx:438,443` + `ShapeAnalysis_Curve.cxx:376-400` | **6** | pending |
| T-57 | A21 | `occt-topo/src/meshing/model_healer.rs:279-290` | 退化支路左右端接反、丢 `aPrevSqDist - aNextSqDist` 判定 → `BRepMesh_ModelHealer.cxx:491-512` + `hxx:143-151` | **6** | pending |
| T-58 | A22 | `occt-topo/src/step/p05.rs:348-370,740-760` | `ProjectAct` 缺 Ellipse/Parabola/Hyperbola 精确臂；圆用三点外心回退 → `ShapeAnalysis_Curve.cxx:382-400,160,200` | 6 | pending |
| T-59 | A23 | `occt-topo/src/wireframe.rs:462-599` | 9 点 pcurve 采样当 UV 包围盒 → `BRepTools.cxx:172-330` `AddUVBounds`（精确，B-spline 走控制多边形） | 5 | pending |
| T-60 | A24 | `meshing/range_splitter/p01.rs:86-133` | 周期标志 + 半径采样猜面型 → `GetType()` 分派（`BRepMesh_FaceDiscret.cxx:112` + `MeshAlgoFactory.cxx:64`） | 5 | pending |
| T-61 | A25 | `meshing/delaun/p04.rs:346-375` | 自造"先删邻三角形再 AddElement" → `BRepMesh_Delaun.cxx:2263-2274` 失败即置 `IMeshData_Failure`，不改网格 | 7 | pending |
| T-62 | A26 | `brep_exchange.rs:118,125`、`occt-core/src/io/{ply,stl}.rs`、`iges.rs:168,390-438`、`step/p02.rs:16-17,43`、`vrml.rs:92`、`obj.rs` | PLY 焊接/属性类型、STL 阈值/头/嗅探、IGES 采样族与自造回转面、STEP 写侧采样重拟、`solid TRUE`、恒空 `vn` → 各 `RWPly_*`/`RWStl*`/`GeomToIGES_*`/`GeomToStep_MakeCurve.cxx:94-99`/`VrmlData_ShapeConvert.cxx:360` | 8 | pending |

> 执行纪律（本轮：**先提交、再任务化、再开工**）：commit `bcbc7dc` 已把审查前的全部工作树入库（8 个提交，工作树干净），此后每个修复单独成 commit，便于 A/B 与回滚。

**批 1（T-35 / A0）已完成 —— 2026-09-20**

- **改动**：`crates/occt-geom/src/offset.rs`（整文件重写）、`crates/occt-geom2d/src/offset.rs`（整文件重写）、`crates/occt-geom2d/src/curve.rs`（补 `Curve2d::d3`，对应 `Geom2d_Curve::EvalD3`）。
- **对齐内容**：`CalculateD0/D1/D2` 逐行移植（含 `R/R2/R3/R5` 与 `Dr/D2r` 两个稳定性分支、`theIsDirChange` 的 `D2.Reverse()`）；3D 法向 `Ndir = D1 ^ Direction`（沿法向偏移，**不再是沿方向平移**）；2D 法向 `(D1.Y(), -D1.X())`（**不再反号**）；`Continuity` 用 OCCT 的 C1→C0 / C2→C1 / C3→C2 表（`Geom_OffsetCurve.cxx:229-257`、`Geom2d_OffsetCurve.cxx:181-210`）；`Reverse` 改为**基曲线也反转**（`cxx:95-100` / `:90-95`）；`Transform` 3D 改为变换基曲线+方向+**带符号** scale（`cxx:454-460`），2D 保持 `abs`（`cxx:414-419`，与 3D 不同，OCCT 确实如此）。
- **容差**：用 `gp::Resolution()` 的忠实值 `precision::REAL_SMALL`（`DBL_MIN`），不用仓内 `precision::RESOLUTION`（`1e-12`，非 OCCT 值）。
- **验证**：① `occt-geom`/`occt-geom2d`/`occt-topo` `--all-targets` 编译 exit 0；② 门禁与基线**逐项一致**：topo lib 1293/1、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3、geom 153/153、geom2d 72/72；③ **解析恒等式探针**（临时，已删）：圆基（R=2）+ 方向 +Z + 偏距 0.5 ⇒ D0/D1/D2 与半径 2.5 的解析圆 **max 误差 0.000e0**（17 点），`Reverse` 后 `|P(0)|=1.5=R−d` ✓。
- **覆盖事实（重要）**：`data/*.step` 中**没有 `OFFSET_CURVE`** 实体（只有 `OFFSET_SURFACE`，属 `Geom_OffsetSurface`，非本项）；3D 类的唯一构造点 `step/p04.rs:1049` 故**无门禁覆盖**；2D 类**无任何构造点**（`Geom2dOffsetCurve::new` 无调用者，`geom_bnd_lib_offset2d::box_offset` 亦无调用者）⇒ **A0 是潜在正确性缺陷（latent），修复后不影响任何现有门禁数字**。审查报告 §3 中"2D 版被 `geom_bnd_lib_offset2d.rs:101` 消费"应订正为"该消费者自身为死代码"。
- **未移植（已在文件头登记，另立 T-63）**：`EvalD3` 缺失导致 `D2Ndir` 项为 0、`EvalD2` 的 `AdjustDerivative` 奇异支路、OCCT 抛异常在无失败通道下回退为基曲线值（与仓内 `offset_surface.rs`/`bspline_surface.rs` 既有约定一致）、ctor 的 C0 拒绝与 G1 升级。
- **派生新发现**：**A27（T-64）**——解析曲线/曲面 `continuity()` 全报 `3`(G2)，OCCT 为 `GeomAbs_CN`(6)。

**批 2（T-42 / A6）已完成 —— 2026-09-20**

- **改动**：`occt-core/src/gcpnts.rs`（删自创采样器）、`occt-geom/src/gcpnts.rs`、`occt-geom/src/curve_approx.rs`、`occt-topo/src/wireframe.rs`、`occt-topo/src/meshing/geom_tool.rs`。
- **删除的自创实现**：`UniformDeflection` / `QuasiUniformDeflection` / `TangentialDeflection` 三个结构体及其 `subdivide`/`subdivide_tangent` 递归（`MAX_DEPTH=16` 截断、中点弦偏差、`point_segment_dist`），以及 `geom_tool.rs` 的自造 `enforce_min_points`（OCCT 的最小点数逻辑在 `GCPnts_TangentialDeflection::PerformCurve` 内部的 `fill_min_points`，已在忠实件里实现）。
- **改走忠实件**：`wireframe::edge_to_polyline` → `meshing::edge_discret::CurveTessellator::from_range`（= `GCPnts_TangentialDeflection::initialize` `cxx:415-453` + `PerformCurve` `cxx:522-916`），最小点数取 `tessellator_min_points`（`BRepMesh_CurveTessellator.cxx:100-114`）；`occt-geom::tangential_deflection` → `perform_tangential_curve`（CN 区间 + BSpline/Bezier 最小点数提升，与 `edge_discret` 同源）；`curve_approx::curve_to_polyline` → 同一引擎（角向项关闭，`RANGE` 约定同 `CurveTessellator::from_range`）；`geom_tool::discretize_curve`/`discretize_iso_curve` → `CurveTessellator`（iso 走 `Surface::u_iso_curve`/`v_iso_curve` = `Geom_Surface::UIso/VIso`）。
- **登记为未移植**：`GCPnts_UniformDeflection`、`GCPnts_QuasiUniformDeflection`、`GCPnts_UniformAbscissa`（`occt-core/src/gcpnts.rs` 模块头列出 OCCT 文件与未移植理由）；`occt-geom::uniform_deflection` 已删除（原为伪造的 "UniformDeflection 参数访问器"，无生产调用者）。
- **验证**：三 crate `--all-targets` 编译 exit 0，改动文件零警告；门禁与基线**逐项一致**（topo lib 1293/1、parity 14/14、step_to_obj 13/13、area 11/11、geom_parity 2/3、geom2d 72）；`occt-core --lib` **293→290**、`occt-geom --lib` **153→151**——减少量恰好是删掉的自创采样器测试（3 + 2），非回归。
- **派生新发现（假出处，并入 T-45/A9）**：`occt-geom/src/curve_approx.rs` 头注释引用的 `GeomConvert_CurveToPolyline` **在 OCCT 8.0.0 中不存在**（全树无 `*CurveToPolyline*` 文件），原注释同时误称 `occt-core` 已实现 `GCPnts_UniformDeflection`；两处均已改正。
- **遗留**：`occt-geom/src/gcpnts.rs` 的 `curve_length`/`abscissa_point` 仍用自适应 Simpson + 牛顿（OCCT 是 `math_GaussSingleIntegration` + `math_NewtonFunctionRoot`），属 **A15/T-51**，未在本批动手（该文件头已注明）。

**批 2 续（T-43 / A7，点–曲线侧）已完成 —— 2026-09-20**

- **改动**：`occt-geom/src/extrema_pc/{mod,p02,p03,tests}.rs`、`extrema_cc/mod.rs`（文档）。
- **对齐内容**：`point_curve_extrema_all` 由"采样分类 + 自创网格路径"改为 `Extrema_ExtPC` 的真实结构——`Extrema_ExtPC` 在 V8_0_0 是 `Extrema_GGExtPC` 的别名（`Extrema_ExtPC.hxx:31-38`），其类型分派：line/circle（本仓可判定）走 `Extrema_ExtPElC` 解析臂并带上 OCCT 的 `myIsMin`（线 = min `cxx:77`；圆 = near 为 min、对径点为 max `cxx:177-188`），其余走 `p03` 的 `default:` 臂（`hxx:390-502`，`aMaxSample=17` + `DeflCurvIntervals`）。
- **删除的自创实现**：`p01` 的 `newton_point_curve_all`/`build_samples`/`solve_f_zero`/`refine_seed`/`fprime`/`dist2`/`fval`/`pair` 与采样分类器 `is_line`/`circumcenter`/`classify_circle`（共 ~254 行），以及 `p02::param_for_point`。
- **语义订正（重要）**：原实现**无条件把区间端点当作极值**加入结果；OCCT 不会——`Extrema_ExtPElC` 只保留落在 `[Uinf, Usup]` 内的解（`cxx:180`），`default:` 臂仅在点与端点重合时补端点（`GGExtPC.hxx:474-502`），端点处理属调用方（`ShapeAnalysis_Curve::Project` `cxx:161-182`）。已按 OCCT 去掉该自创补端点。
- **验证**：`occt-geom --lib` **151/151**、`occt-geom2d` 72/72、`occt-topo --lib` 1293/1（唯一红仍是 `brepfeat::tests::groove_cuts_cylinder` = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（T-05）——**逐项与基线一致**；`--all-targets` exit 0。
- **未修（已立项）**：**T-66** = `extrema_cc` 的通用曲线–曲线种子集仍是自创网格（OCCT `Extrema_GenExtCC` 用 `math_GlobOptMin`，仓内 `occt-math/globoptmin.rs` 已移植，可直接用）；**T-65** = `ellipse_all`/`hyperbola_all`/`parabola_all` 解析臂缺 `myIsMin` 端口（`Extrema_ExtPElC.cxx:277`/`:381`/`:473`）且 `Curve` 无 hyperbola/parabola 类型查询（T-12），故这三种类型当前走 `default:` 臂（已在代码注释与 `extrema_cc/mod.rs` 头部写明）。
- **同轮试改并回退（T-64 / A27）**：把 15 处解析曲线/曲面的 `continuity()` 从 `3`(G2) 改成 OCCT 的 `6`(CN) 后，`step_to_obj` 由 13/13 变 **12/13**——**Sphere 网格变空**（round-trip 顶点数 0）。判定：**不能单独改**，某消费方把"报出的连续性"当跨度上限（嫌疑 `meshing/range_splitter/p01.rs:295` 的 `continuity <= curve.continuity()`、`meshing/edge_discret.rs:1259` 的 `pcurve.continuity().min(surface.continuity())`）。已按纪律**回退**（`step_to_obj` 恢复 13/13、geom 151/151），并在 15 处就地标注 `// T-64: OCCT=GeomAbs_CN(6), blocked by consumer`、在 `geomabs.rs` 及本表记录依赖；下次须先对齐消费方的 `GeomAdaptor_Surface::NbUIntervals/NbVIntervals` 语义。

**批 3 前置核查（T-37/A1）—— 2026-09-20：审查前提被推翻，改为先移植**

- 审查报告 A1 的整改建议写着"仓内已有忠实件 `extrema_surf::point_surface_extrema(_box)`（`Extrema_ExtPS`）"——**该前提不成立**。实读：`extrema_surf/p02.rs:74-119` 的 `point_surface_newton_all_box` 是 **24×24 网格播种**（局部极值 + 全局 min/max + 角点）+ Newton，Jacobian 用**数值差分**（`p01.rs:620-621`），模块头还写着 "Surface trait lacks d2"——而 `Surface::d2` **存在**（`surface.rs:12`，解析曲面覆写）。⇒ 它和 `brep_surface::surface_closest_params` 属**同一类替代品**；把 39 处调用点迁过去等于"自创换自创"。
- 忠实路径 = **移植** OCCT 的 `Extrema_ExtPS.cxx`(376 行) + `.hxx`(140)、`Extrema_GenExtPS.cxx`(1056) + `.hxx`(155)、`Extrema_ExtPElS.cxx`(454) + `.hxx`(76)：类型分派 + iso 退化处理（`IsoIsDeg`）+ 逐 C2 区间采样（`mySample`）+ `Extrema_GFuncExtPS` 的解析系统与 `math_FunctionSetRoot`，Jacobian 改用解析 `d2`。
- 已立项 **T-67**（A1 前置，T-37 被其阻塞），并把订正写回 `specs/_audit/_index.md` 的 A1 行与 §7 第 2 条、`extrema_surf/mod.rs` 模块头。**本批未改任何行为代码**，门禁不受影响。

## 4. 决策与约束（不可违反）

1. 改完先编译（`cargo check`，编译不过先修编译）。
2. 改代码前对着 `.cxx` 审控制流；**有同等分支才改**，没有就标"未移植"（注释写 OCCT 文件+行号），**不加** OCCT 里不存在的谓词/启发式/面积比/长度滤边/体积门。
3. 对齐验证 = 跑**现有门禁/基线**（`step_obj_parity`、`step_to_obj`、`step_obj_area`、`step_geometry_parity`、`export_data_obj` 导出、`data/occ-*.obj` 对拍）；**不为对齐新写单元测试**，不用假输入证明算法已对齐。
4. 完成即停：门禁满足 → 报告结果 + 旁支问题 → 停，等指令。旁支只报告不顺手修。
5. 任务卡门禁只能引用**现有**基线；源码行数只作批次范围，不作成果汇报。
6. 共享树禁止 `git stash`；A/B 用 `git worktree`（§2）。
7. 本 wave 新增 **0 个 `#[test]`**（新模块 13,900 行）——验收依赖门禁与导出门；这与规则 3 一致，但意味着**新单元的回归保护靠门禁覆盖**，故 T-21/T-22 的覆盖扩展优先级不低于补缺。

## 5. 关键发现（findings）

- **HEAD 比工作区红得多**：同一套 `--no-fail-fast`，HEAD 是 lib `1,235/59 失败`、phase20 `1/5`、phase3 `3/4`、boss `0/2`；工作区收敛到 lib `1,292/2`。⇒ 未提交 wave 主要是**修绿**波次，不是新功能波次。
- **只有 1 个新增红**（T-06），且是期望过期（ellipse 参数化忠实化的副作用），不是行为回归。
- **`builder_solid.rs` / `brepfeat/tests.rs` 未改**却数字变化 ⇒ 行为来自共享依赖（`primitives.rs`、`bop_build_solids*`、`clib.rs` 等），定位要往依赖里走。
- **mtime 陷阱**：几乎所有 `src/*.rs` 的 mtime 是 `2026-09-18 00:55:43`（一次全树改写），但 diff 是真实内容 ⇒ **不要用 mtime 判断进度**，只认 `git diff` + 门禁输出。
- **最后一次成功动作**：`export_data_obj` 于 2026-09-17 16:40 跑通 16 个模型（`output/_run_step0.txt` 15:05 全绿）。
- 其它 worktree（Warp `flint-limestone` / `pediment-travertine`、`wt-residual`、`D:\source\.wt-ee`）**均干净或目录已失**，无待办工作。
- **`IsHole` 靠朝向，不靠"嵌套"**（2026-09-20 实证）：`BOPAlgo_BuilderSolid::IsHole` = 把壳当实体做 `PerformInfinitePoint`，只有面朝向指向自身空腔的**反转壳**才判 `IN`。因此"一个正朝向盒子套在另一个正朝向盒子里"在 OCCT 里是**两个 growth → 两个 area**；把它当"腔体"是夹具错误，不是算法缺陷。仓内已有同向证据：`algo_tools/tests.rs::is_inverted_solid_normal_box_false`。
- **探针法可用**（临时，用完即撤）：本仓没有 root workspace，单测聚焦跑法 `cargo test --manifest-path crates/occt-topo/Cargo.toml --lib <名字> -- --nocapture`，`eprintln!` 探针经 `--nocapture` 可见。
- **BOP 层症状（2026-09-20 实测）**：平面精确布尔在"柱 × 环工具"配置下 `Fuse=0 面 / Common=0 面 / Cut=116 面但欠切`，且**不报警告**。op 分派已对上 OCCT（`BuildShape` = `BuildRC` + 仅 FUSE 走 `BuildSolid`），所以不一致在 `BuildRC` 的集合构造/状态过滤（`bop_bop.rs:187-276` vs `BOPAlgo_BOP.cxx:583-711`）。这条同时是 T-03（boss 合并）与 T-04（phase19）的邻居区域，值得作为一个"BOP 簇"统一收敛。

## 6. 一条命令刷新快照

```powershell
cd D:\source\repos\dogs
$M = 'crates/occt-topo/Cargo.toml'
cargo check --manifest-path $M --all-targets                                  # 门禁 1
cargo test  --manifest-path $M --lib                                          # 门禁 2（关注 1294 里的失败数）
cargo test  --manifest-path $M --no-fail-fast 2>&1 | Select-String 'test result:|Running tests|panicked at'  # 门禁 3 全量
cargo test  --manifest-path crates/occt-geom/Cargo.toml --lib                  # 门禁 4
# 导出门禁（16 个 data/*.step → output/<stem>.obj，偏转 0.1）
cargo run --manifest-path $M --offline --example export_data_obj
```

ATU01038 bbox 对拍（无测试覆盖时的临时手段；结果记进 §7）：
解析 `output/ATU01038.obj` 与 `data/occ-ATU01038.obj` 的 `v ` 行取 min/max，当前 Δ≤9e-6。

## 7. 进度日志

### 2026-09-19 23:00–23:55 · DSH 会话（只读盘点，未改代码）
- **动作**：定位在办任务 → 盘点 178 文件未提交 wave → 跑门禁 4 组 → 用 `git worktree add --detach .target-headcheck HEAD` 建 HEAD 基线并跑全量对照 → 删除临时 worktree（`git worktree list` 已确认移除，工作区恢复 148 M / 30 D / 32 ??）。
- **结果**：见 §2 快照。核心结论：wave 把 ~57 个失败用例修绿；剩 6 个红（4 个 HEAD 遗留 + 1 个 HEAD 遗留改善 + 1 个新期望过期）。
- **产出**：本画板 `specs/_board.md`。
- **未做**：未改任何源码/测试；未提交；未跑 `export_data_obj`（沿用 9/17 16:40 产物与其日志）。

### 2026-09-20 00:00–01:00 · DSH 会话（执行：P0 收敛 3 条）
- **T-06 修**：`extrema_pc/tests.rs` 期望 `3π/2 → π/2`（附 `ElCLib::EllipseValue` 出处）→ `occt-geom --lib` **153/153 ✅**。
- **T-21 加**：`step_obj_parity.rs` 新增 `atu01038_bbox_matches_occt`（`check_parity("ATU01038","occ-ATU01038.obj",1e-4)`，实测 Δ≤1.1e-5）→ **13/13 ✅**。
- **T-02 修**：`eprintln!` 探针定位到 `growth=[0,1] holes=[]`；对 `BOPAlgo_BuilderSolid::IsHole` 求证后判定为**测试夹具喂了非腔体输入**，改为反转内壳朝向（真实 Cut 腔壳形态）→ `builder_solid` 全绿，topo lib **1,293/1,294**（仅剩 T-01）。探针已撤除，生产代码 `perform_areas` 未改动。
- **产出**：三处改动 = `crates/occt-geom/src/extrema_pc/tests.rs`、`crates/occt-topo/tests/step_obj_parity.rs`、`crates/occt-topo/src/builder_solid.rs`（仅测试模块）。
- **未做**：未提交；暂未跑全量 `--no-fail-fast`（下一次提交前必须跑一次）。

### 2026-09-20 01:00–02:00 · DSH 会话（goal round 1：T-01 定位到 BOP 层）
- **探针 3 轮**（均已撤除，`brepfeat/tests.rs` 与 HEAD 逐字节一致）：
  1. `tool=3.257401`（= 16 边形环解析值 ✓）、`common=0.000000`、`cut=8.903438`（removed 0.4140）、`faces_tool=128 faces_cut=116`；
  2. `multiwire=0` ⇒ 环面不是"带孔 wire"结构，排除孔线假设；`disc_common=0 / disc_fuse=0 / ann_fuse=0` ⇒ **Fuse 也空**；
  3. `planar_cyl=true planar_tool=true`，直接调 `bop_builder::boolean`：`common faces=0 / fuse faces=0 / cut faces=116`，`common_warn=[]`。
- **结论**：T-01 不在 `brepfeat`/`groove`/revolve 层（工具体积正确），而在**平面精确布尔的 `BuildRC`**；分派与 OCCT `BOPAlgo_BOP::BuildShape`（`:871/900/902-904`）一致，故问题是集合/状态过滤。HEAD 同红 ⇒ 非本 wave 回归。
- **门禁复核**：topo lib 1,293/1,294（仅 T-01）、geom 153/153、parity 13/13 —— 探针移除后无残留、无新增红。
- **下一轮入口**：`bop_bop.rs:187-276 build_rc` 四集合探针 + 对 `BOPAlgo_BOP.cxx:583-711` 逐段核对。

### 2026-09-20 02:00–03:00 · DSH 会话（goal round 1 续：T-01 根因链打通）
- **探针 2 轮**（均已撤除；`bop_bop.rs`、`bop_builder2/p02.rs`、`brepfeat/tests.rs` 三文件复核为未修改状态）：
  1. `build_shape`/`build_rc`：`op=Cut dim0=3 dim1=3 open=false gf_solids=2 gf_faces=256`；`build_rc obj_src=1 tool_src=1 it_shapes=1 check_keys=1 it_exp=1 tool_sets=1`；
  2. `fill_images_solids`：`fill_in3d in_parts=1`，`build_split_solids_occt` 后 **两个 argument 都有镜像，但每个只有 1 个**。
- **结论（链路已闭合）**：T-01 的表现是"Cut 只切 0.414/1.6、Common/Fuse 空"，根因是**每个源 solid 只产出一个镜像（=未被切割）**，即 `BuildSplitSolids` 拿不到柱侧片被工具面切开的碎片 ⇒ 上溯到**面级切割**（`FillImagesFaces` / PaveFiller 的 planar-planar FF 干涉）。
- **门禁复核**：topo lib 1,293/1,294（仍仅 T-01），探针零残留。
- **下一轮入口（明确）**：`bop_build_faces::fill_images_faces` → 统计柱侧片 `history.image(face)` 条数；查 PaveFiller 是否给出该面对的 FF 干涉；对 OCCT `BOPAlgo_Builder_2.cxx::BuildSplitFaces`。此链同时是 T-03（boss 合并）、T-04（phase19 DS 计数）的共同上游，建议按"BOP 簇"一次收敛。

### 2026-09-20 03:00–04:00 · DSH 会话（goal round 2：面级切割定位）
- **探针 3 轮**（全部撤除；`bop_builder2/p02.rs`、`bop_split_faces_occt.rs` 复核为未修改、零 `eprintln` 残留）：
  1. DS 概览：`ds_shapes=1238 src_shapes=678 src_faces=224 pave_blocks=272 common_blocks=0`，`interf_ff=422 interf_ee=16 interf_ef=128 interf_vf=0` ⇒ **干涉检出不缺**；
  2. 面级绑定：`ff_curves=144 ff_pb=144 ff_pb_with_edge=144`（截面边全有真 edge），`faces_sc=112 faces_on=112 faces_in=0` ⇒ **pave 绑定不缺**；
  3. 切割重建（`build_split_faces_occt`）：`tasks=112 empty_le=0 draft=0 faces_im=112 with_multi=31 max_images=3`，`area_hist=[(0,40),(1,41),(2,30),(3,1)]` ⇒ **40 面拿到 0 块、41 面只有 1 块（=没切开）**。
- **结论**：T-01 的断点在 `rebuild_split_areas` → `builder_face_occt::perform`（`BOPAlgo_BuilderFace`）。前两轮的三个假设（FF 未检出 / pave 未绑定 / solid 镜像缺失是原因）**全部被实测排除**：前两者正常，第三者是结果不是原因。
- **已核对的 OCCT 对照点**：`BOPAlgo_Builder_2.cxx:233 BuildSplitFaces`（`:496-500` 平面预建 pcurve、`:530` `aFacesIm.Add(..., aBF.Areas())`——注意 OCCT 会为**空 areas 也绑定一个空列表**，Rust `bind_face_images` 对空 vec 不加 key，下游 `splits_of` 两者都回落原面，等价）、`BOPAlgo_BuilderFace.cxx:401/461`。
- **门禁复核**：topo lib 1,293/1,294（仍仅 T-01）。

### 2026-09-20 04:00–05:00 · DSH 会话（goal round 3：WireSplitter 顶点身份）
- **探针 2 轮**（撤除后 `builder_face_occt.rs` 回到本 wave 自身的 3/3 diff，零 `eprintln` 残留）：
  1. `FaceBuilder` 阶段统计（112 次调用）：全是平面；`avoid_sta=0`（ShapesToAvoid 从未触发）；分组 `(in_edges, avoid, loops, areas)` ⇒ `avoid=0` 组 = 正常切割（11→3/2、15→5/1、9→2/2），`avoid=3` 组（45 面）= 30 面 0 块 + 15 面 1 块；
  2. 追加 `avoid_bnd`（被弃边里有几条属原面边界）：**所有 `avoid=3` 面都是 `avoid_bnd=1`** ⇒ 每面丢的是"1 条边界碎片 + 2 条截面边"。
- **结论**：断点在 **`perform_loops` 的 WireSplitter 链接**——截面边端点与边界碎片端点未共享 vertex（或边界碎片缺失），导致成环失败；`perform_areas` 随后把唯一成环的截面环判为 hole，又因无 growth 而在闭合面上被丢弃（`builder_face_occt.rs:389-408` 只对开放面兜底）。已核对 OCCT `BOPAlgo_BuilderFace::SetFace` 也是强制 FORWARD（`:79-83`），故"朝向被强制"这一假设**排除**。
- **下一轮入口**：`MakeSplitEdges`/pave 顶点身份 —— OCCT `BOPAlgo_PaveFiller_7.cxx:371`（新截面边必须用 DS 的 pave 顶点建立，才能与边界边分割点共享 vertex）、`:589 MakePCurves`。**已核对**：Rust 的边界碎片侧是对的（`pave_blocks/p01.rs:314 make_split_edges` 用 `f.ds().shape(t.v1/v2)` 建边），所以嫌疑收窄到**截面边（FF 新建边）一侧的顶点/端点**：需验证 FF 曲线端点是否复用了边界 pave 的同一 DS vertex（对应 OCCT 的 `BOPDS_DS::Index` 复用 / `IsNewShape`），而不是在同一点再建一个 vertex。探针配方：取一个 `avoid_bnd=1` 的坏面，打印被弃的 1 条边界碎片与 2 条截面边的端点 `shape_key` 与 3D 坐标，看是否"同点不同 vertex"。
- **门禁复核**：topo lib 1,293/1,294（仍仅 T-01）。

### 2026-09-20 05:00–06:00 · DSH 会话（goal round 4：WireSplitter 共享边消费）
- **探针 3 轮**（全部撤除；`builder_face_occt.rs` 回到本 wave 自身的 3/3 diff，全仓 `*-dbg` 残留 = 0）：
  1. **坐标重复检测**：112 面 `dup_pos_groups` 全为 0 ⇒ "同点不同 vertex"**被推翻**；
  2. 被弃 3 边几何：构成一条**开口链** `(1,0,1.8)→(1,0,3)→(0.9659,-0.2588,3)→(0.9796,-0.1553,1.8)`（两端 z 都是 1.8 但不重合）；
  3. 该面**全部输入边**：9 条唯一边（2 条水平弦边各 2 份副本）；`in_loop=true` 的边恰好闭合下面两段，**最上面一段 4 条边齐全却未成环**（竖边 z1.8→3、顶边、斜边 z3→1.8、一份弦边）。
- **结论**：断点收窄到 **`perform_loops`/`WireSplitter` 对"被相邻两块共用的截面边"的消费**（OCCT 同一 section 边以两个朝向各一份进入 `BOPTools_WireEdgeSet`；Rust 输入已有 2 份副本却仍失败）。
- **下一轮入口**：`wire_splitter.rs::perform`/`split_edges` vs OCCT `BOPAlgo_WireSplitter::Perform` + `BOPAlgo_BuilderFace::PerformLoops`（`BOPAlgo_BuilderFace.cxx:256+`）。
- **门禁复核**：topo lib 1,293/1,294（仍仅 T-01）。

### 2026-09-20 06:00–07:00 · DSH 会话（goal round 5：命中 `RefineAngles` 空实现）
- **探针 3 轮**（全部撤除；`wire_splitter.rs` 仅剩本 wave 自身的 76/4 改动，`wire_splitter_block.rs` 未被修改，全仓 `*-dbg` 残留 = 0）：
  1. 连通块/成环统计（112 面）：坏面 7 边→1 环 3 边、11 边→1 环 5 边、13 边→1 环 6 边 ⇒ 部分环；
  2. `split_block` 内部：62 次成功调用 `skipped_no_pcurve=0` ⇒ **pcurve 缺失假设排除**；另 50 个面走 `bNothingToDo` 早退（单环，未切）；
  3. 对读 OCCT：`SplitBlock`（`_1.cxx:112-354`）+ `Path`（`:358-617`）与 Rust 逐行一致；键为指针身份（`GeometryRegistry::shape_key`）。
- **根因**：`wire_splitter_block.rs:541 refine_angles` **只实现了 `iCntBnd != 2` 早退**，对 `iCntBnd == 2` 完全不做；OCCT 恰在该分支用 `RefineAngle2D` 重算内部出边角度（`_1.cxx:925-1032` / `:1033-1125`，`iCntInt==2` 时按 `Precision::Angular()` 微调，`:998`）。坏面弦端点正是 `iCntBnd==2 && iCntInt==2`。
- **round6 计划**：移植 `RefineAngles(vertex)` + `RefineAngle2D`（含 2D 曲线求切向）替换 stub → 复跑 T-01 + 四组门禁；若成立，预期同时改善 T-03/T-04（同一 SplitBlock 链）。

### 2026-09-20 07:00–08:30 · DSH 会话（goal round 6：补 `RefineAngles`，结果阴性但缺口已闭）
- **实施**：`wire_splitter_block.rs` 的 `refine_angles` stub → 忠实移植 `RefineAngles`（`_1.cxx:905-918` + `:925-1029`）与 `RefineAngle2D`（`:1033-1125`）；新增 import `Geom2dLine`、`ANGULAR`。
- **探针 2 轮**（已撤除，零残留）：① 精修触发计数 = **160 次**（全 `iCntInt=2`，`aA1==aA2`）；② `area_hist` 与 round3 基线**逐项相同**。
- **结论**：`RefineAngles` 在本例中只产生 `1e-12` 的微调（楔形退化），不足以改变 `Path` 的选边 ⇒ **不是 T-01 的操作性缺口**；但它是一处真实未移植分支，保留（并已通过全量回归）。
- **全量门禁**：lib 1293/1、boss 1/2、phase19 3/5、phase20 5/5、phase3-10 全绿、step_geometry_parity 2/3、step_obj_area 11/11、step_obj_parity 13/13、step_to_obj 13/13 —— **与基线持平或更好，无回归**。
- **round7 入口**：`Angle2D`（OCCT `_1.cxx:768-903`）vs Rust `angle_2d` + `tolerance_2d`/`uv_tolerance_2d`/`coord2d`。

### 2026-09-20 08:30–10:00 · DSH 会话（goal round 7：定位到"输入边集朝向不自洽"）
- **对读**：`Angle2D`（`_1.cxx:768-840`）与 `angle_2d` 逐行一致 ⇒ **角度链无偏差**；`Tolerance2D`（`:859-881`）与 `uv_tolerance_2d` 有偏差（平面面 1.0 vs 1e-7）⇒ 新开 **T-30**（非 T-01 阻塞）。
- **步进追踪**（`Path` 前 8 次调用，探针已撤除、零残留）：多次出现 `DEAD END at vb_slot=… i_cnt=0`；`ls` 已累积 5–8 条边却从未回到链上旧顶点 ⇒ 不是"选错边"，而是**根本无出边可走**。
- **块级 dump（决定性）**：该块 7 个顶点中，顶角 `(1,0,3)` = **0 in / 2 out**、底角 `(0.966,0.259,0)` = **2 in / 0 out**，其余顶点 in==out（2/2 或 1/1）⇒ 该面的边集**朝向不自洽**（三条边同向流），几何上根本不是闭合环 ⇒ `Path` 必然死路、面必然不切。
- **下一轮入口**：`orient_split_from_base_with_warn`（`bop_split_faces_occt.rs:152`）vs OCCT `BOPTools_AlgoTools::IsSplitToReverse`；以及 `:104` 未切边界边的出现朝向是否符合 `_2.cxx:468-479`。
- **门禁复核**：lib 1,293/1,294（仍仅 T-01）；本会话新增改动仅 `wire_splitter_block.rs`（+128/−12，`RefineAngles` 移植）。

### 2026-09-20 10:00–11:30 · DSH 会话（goal round 8：根因 = mesh→BRep 的 wire 不成链）
- **对读**：OCCT `BuildSplitFaces` 取边规则（`_2.cxx:412-494`）与端口一致；`IsSplitToReverse` 的 same-curve 捷径（`BOPTools_AlgoTools.cxx:1462`）vs 端口 `algo_tools_face.rs:727`（比 TShape）⇒ 新开 **T-31**。
- **探针**（已撤除，0 残留）：112 面 **544 个边界碎片、`orient_mismatch=0`** ⇒ 朝向施加忠实；再对 `le=11` 的面 dump 三条边界边：**全 Forward**，`(1,0,3)→(0.966,0.259,3)`、`(0.966,0.259,3)→(0.966,0.259,0)`、`(1,0,3)→(0.966,0.259,0)` ⇒ **wire 不是链**（两条边从同一点出发、两条汇入同一点），与 round7 的 source/sink 完全对应。
- **根因**：`mesh_to_brep.rs:89-102` 的 `ensure_edge` 按顶点对缓存共享边，邻接三角形先创建后用**相反方向**复用，且从未做 `CorrectEdgeOrientation` ⇒ 每个 mesh→BRep 面的 wire 都可能不成链。
- **实验（已回退）**：只在 `:102` 加"按请求方向翻转"→ lib **1292/2**（`groove_cuts_cylinder` 前移到 `:98` 减材断言，after≈before；`groove_negative_volume_delta` 新红）⇒ 说明面朝向/壳装配必须同改（对应 OCCT `BRepBuilderAPI_MakeFace` 的绕向-法向一致性），单改 wire 不安全。已 `git` 层面复原，基线回到 **1293/1**。
- **round9 计划**：成套移植 `MakeWire::Add` 的 `CorrectEdgeOrientation` + `MakeFace` 的面朝向规则（`triangulation_to_brep`），复跑 T-01 + 全量门禁；预期 T-01、并可能带动 T-03/T-04。

### 2026-09-20 11:30–12:30 · DSH 会话（goal round 9：T-31 修复 + 绕向度量）
- **T-31 已修**：`algo_tools_face.rs::is_split_to_reverse_edge` 的捷径由"edge TShape 身份"改为 **`Arc::ptr_eq(edge_curve(..))`**（对齐 OCCT `BOPTools_AlgoTools.cxx:1461-1465` 的 same-curve 判断）；全量 `--no-fail-fast` **与基线逐项一致**（lib 1293/1、boss 1/2、phase19 3/5、phase20 5/5、phase3–10 全绿、step_geometry_parity 2/3、step_obj_area 11/11、parity 13/13、step_to_obj 13/13），T-01 体积仅末位差 ⇒ 无回归、本例惰性。
- **绕向度量（探针，已撤除）**：`面法向 · wire 绕向法向` 的符号统计 —— `cyl: 96 面（49 正 / 47 负）`、`tool: 128 面（16 正 / 112 负）`。mesh→BRep 柱面**接近对半** ⇒ 绕向约定混乱（与 T-32 的开链同源）；旋转体工具**一致为负** ⇒ 其构造路径自洽。这解释了 round8 实验为何"链好 wire 反而更糟"：链向修好后，面/壳的材料侧与既有分类约定不再一致，必须成套对齐 `BRepBuilderAPI_MakeFace` 的绕向-法向规则。
- **round10 计划**：T-32 成套实施（`triangulation_to_brep`：链向 + 面朝向 + shell/solid 装配一致性），先小步验证再全量；若仍回归则回退并转做 P3（`data/occ` 7 个新模型的 Rust 导出对拍）。

### 2026-09-20 12:30–13:30 · DSH 会话（goal round 10：T-32 前置核查 + P3 新模型对拍）
- **T-32 前置核查**：读了 `GpAx3::from_ax1`（`crates/occt-core/src/gp/ax3.rs`）——其 `ydir = z×x`、`xdir = y×z` ⇒ **右手系**，故"坐标系手性导致绕向翻转"这一假设**排除**；T-32 仍需与面/壳装配约定成套对齐（跨轮）。
- **P3（本轮落地）**：
  1. `examples/export_data_obj.rs` 支持可选目录参数（默认仍 `data/`，既有门禁行为不变）：`cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj -- data/occ`；
  2. 8 个新模型全部导出成功（ATU01038/T0M/TDB/a3n00/acs10/bottom/motoc/top）；
  3. 对拍 `output/occ/*.obj`：**5 个 bbox 完全一致**（ATU01038、bottom、motoc、top、T0M；密度差 −8.8%~+0.8%），**3 个有真实几何差**（a3n00 Δy=60.04、acs10 Δy=7.88、TDB Δx=0.59）⇒ 新开待办：逐个查 STEP 实体/朝向解析缺口。
- **门禁复核**：topo lib 1,293/1,294（仍仅 T-01）；`--all-targets` 编译含 example ✓。

### 2026-09-20 13:30–14:30 · DSH 会话（goal round 11：P3 归因 + 门禁锁定）
- **归因（新诊断进入导出门禁）**：`export_data_obj` 现在每条同时打印 **exact bbox**（`brep_bnd_lib::shape_bnd_box`，控制网/解析极值，与网格无关），这样 bbox 漂移可一眼归因到"解析几何"还是"网格"。
  - 实测 `a3n00` exact = `[-209.1,-152.1,-250.9]~[209.1,152.1,100.6]` vs OCCT 参考 `[-134.5,-77.5,-176.3]~[134.5,77.5,26]` ⇒ **解析几何/参数域偏差**（非网格密度）；
  - `acs10` y 半宽 55.56 vs 51.5、`TDB` x 差 0.59 同向偏大 ⇒ 同类。
- **P3 覆盖锁定**：把 4 个已对齐模型的 OCCT 参考拷入 **`data/occ-ref/`**（bottom/motoc/top/T0M），`step_obj_parity` 泛化出 `check_parity_at`（支持 `data/` 下相对路径），新增用例 `occ_test_model_bboxes_match_occt`（容差 1e-3）⇒ 门禁 **14/14 通过**（原 13 + 1）。
- **门禁复核**：lib 1,293/1,294（仍仅 T-01）、geom 153/153、`--all-targets` ✓；新参考文件 `data/occ-ref/`（~6MB）与 `data/occ/` 源模型**待随 T-07 一并入库**。

### 2026-09-20 14:30–15:30 · DSH 会话（goal round 12：a3n00 归因到网格侧）
- **精确度量**：`a3n00` Rust 网格 bbox `[-134.5,-77.5,-176.294]~[134.5,137.539,26]` vs OCCT `[-134.5,-77.5,-176.288]~[134.5,77.5,26]` ⇒ 仅 **y-max 差 +60.039**。
- **排除解析**：遍历 226 个面，**wire 顶点无一超出 |y|>77.6** ⇒ 解析边界与 OCCT 一致。
- **网格侧候选**（逐面单独 `brep_to_obj` 与 wire bbox 对比）：face 41–49、77–87 等多面自身网格超出 wire y ~10（`Plane`/`Other`）；因单面网格化可能与实体导出路径不同，仅作指示。
- **结论**：`a3n00` 的缺口属**网格生成**（超出 trim 采样），与 T-01 的 BOP 链正交；下一步在导出路径内 dump 每面网格 bbox 做精确锁定。
- **门禁复核**：lib 1,293/1,294、parity **14/14**、geom 153/153、`--all-targets` ✓，探针 0 残留。

### 2026-09-20 15:30–17:30 · DSH 会话（goal round 13：全仓移植忠实度审查）

- **目标**：从**全部已移植的 Rust 代码**里找出"没有忠实移植、没有对齐的自创实现"（用户指令）。
- **范围**：5 个 crate（845 `.rs` / 240,919 行）——先 1 轮机器侦察（关键词/模式统计），再按 5 区人工/代理深审（occt-core、occt-geom+2d、BOP/布尔、occt-topo 其余、网格化/数据交换），我对其中 15 条高影响项**逐条对读 `.cxx`**。
- **产物**：`specs/_audit/_index.md`（汇总：口径 + 统计 + A0–A26 + 反证清单 + 优先级）+ 5 份分区报告（`occt-core.md` 15/1/7、`occt-geom.md` 14/1/0、`bop-boolean.md` 13/3/2、`topo-rest.md` 14/4/0、`mesh-exchange.md` 15/8/6）。合计 **自创 71 / 等价替换未登记 17 / 已登记 15**。
- **最高危 A0（公式错误，非近似）**：3D offset 曲线沿参考方向**平移**而非沿法向偏移、2D 法向**反号**，两侧 `d1/d2` 均缺 `DNdir` 项 ⇒ 已开 **T-35**（并提示 T-05 的 offset 体积偏差可能同源）。
- **强正面证据（本轮实测）**：live Delaunay 链与 OCCT 参考 `data/occ-ref/` **逐模型顶点/三角数一致**（Cube 24/12、Sphere 642/1244、Torus 1369/2592、Shape-1 3343/4336、linkrods 3494/5078…）⇒ 主算法对齐良好，缺口集中在**回退路径/类型判定/边参数/UV 域**（A17–A25）。
- **意外收获**：A23 给出 round12 `a3n00` 网格 y+60.039 的**上游根因**（`wireframe.rs:face_uv_bounds` 9 点 pcurve 采样漏 UV 极值 ⇒ 窗口本身偏斜，OCCT 是 `BRepTools::AddUVBounds` 精确盒）。
- **纪律**：审查全程**只读**——未改任何源码；`git status` 与本轮开始时逐项一致（30 D / 153 M / 35 ??，新增未跟踪仅 `specs/_audit/`、`specs/_board.md`）；门禁复核 lib 1,293/1,294、parity 14/14、geom 153/153。
- **任务化**：新增 **T-35（A0）**、**T-36（A1–A26 按 §7 优先级分批）**，状态 `pending`（本轮不动手改）。

---

## 12. goal 12 轮批次总结（2026-09-19 ~ 2026-09-20）

> 说明：§12 / §13 为**追加章节**，物理位置排在 §7 进度日志之后（§8 起为维护协议/坑/五问/范围）。

**已落地（全部经实测验证、无回归）**

| 项 | 内容 | 证据 |
|---|---|---|
| T-06 | 椭圆参数期望 `3π/2 → π/2`（附 `ElCLib` 出处） | `occt-geom --lib` 153/153 |
| T-02 | BuilderSolid 内腔用例夹具改为**真实腔壳**（反转内壳朝向） | topo lib 1292/2 → **1293/1** |
| T-21/T-22 | ATU01038 + `data/occ` 4 个新模型锁进 parity（新增 `check_parity_at`、`data/occ-ref/` 参考入库） | parity 12/12 → **14/14** |
| P2 `RefineAngles` | `wire_splitter_block.rs` 补齐 `RefineAngles`+`RefineAngle2D`（原为空 stub） | +128/−12，全量无回归 |
| P2 T-31 | `is_split_to_reverse_edge` 捷径改比 **curve 句柄**（对齐 OCCT） | 全量无回归 |
| P3 `export_data_obj` | 支持可选目录 + 打印 **exact bbox** 诊断 | 8 个新模型导出成功 |

**T-01 的根因链（12 轮中最重的调查，已闭环到可修点）**：groove 体积不符 → 平面精确布尔 `Common/Fuse=空、Cut 欠切` → `BuildRC` 拿到的 solid 镜像未被切开 → 面级切割只成功 31/112 → `Path` 求解器在"**wire 不成链**"的块上必然死路 → 根因在 **mesh→BRep 造 wire 时未做 `BRepBuilderAPI_MakeWire::Add` 的 `CorrectEdgeOrientation`**（共享边按先请求的三角形方向缓存）。已排除：顶点身份重复、pcurve 缺失、键合并、`Angle2D`/`RefineAngles` 偏差、平面坐标系手性。**一次单独修 wire 朝向的实验回归 2 例并已回退** ⇒ 必须与 `BRepBuilderAPI_MakeFace` 的面朝向/绕向约定成套修（T-32）。

**未完成（交接）**：① T-01/T-32（BOP 链，需成套修 mesh→BRep 约定）；② T-03 boss 合并、T-04 phase19、T-05 offset parity（P0 剩余）；③ T-22 的 3 个模型（a3n00/acs10/TDB）**网格侧**缺口；④ T-30（`uv_tolerance_2d`）、T-33/T-34（本轮审查新增，见 §13）、T-12~T-20（P2 清单）；⑤ P1 收尾提交（未提交改动 + 新数据基线 `data/occ*`、`data/occ-ref/`）。

---

## 13. 已完成任务的移植审查（2026-09-20，逐条对 `.cxx`）

判定口径：**忠实移植** = 逐行对应 OCCT（含分支、常量、控制流顺序）；**等价替换** = 用本仓已有原语替代 OCCT 的某个工具类，语义等价但需在注释里登记翻译边界；**自创** = OCCT 中无对应分支/规则，或产物不是"从 OCCT 翻译来的"。

| 任务 | 改动 | 判定 | 证据 / 说明 |
|---|---|---|---|
| T-06 | `extrema_pc/tests.rs` 期望 `3π/2 → π/2` | **忠实**（期望值按 OCCT 公式重算） | `ElCLib::EllipseValue`（`ElCLib.cxx:176-189`）`P = Loc + Major·cos(U)·XDir + Minor·sin(U)·YDir` ⇒ `(0,1,0)` 对应 `U=π/2`；旧期望编码的是被修正掉的镜像参数化 |
| T-02 | `builder_solid.rs` 内腔用例**夹具**改为反转内壳 | **自创（测试夹具，已标注）** | OCCT 侧只提供算法契约 `BOPAlgo_BuilderSolid::IsHole` = `SolidClassifier::PerformInfinitePoint`（`BuilderSolid.cxx:823-831`）⇒ 只有"绕向反转的壳"才被判 hole；**没有**找到 OCCT 中构造腔壳的调用点可照抄，故夹具是按契约**建模**的，不是逐行翻译。生产代码未改，风险局限在测试输入的代表性 |
| T-21 / T-22 | `step_obj_parity.rs` 加 ATU01038 与 `data/occ` 4 模型；新增 `check_parity_at`；`data/occ-ref/` 参考入库 | **忠实（门禁扩展，非实现）** | 只引用既有基线（OCCT 参考 OBJ）；容差由实测决定并注明（ATU01038 1e-4、新模型 1e-3，实测 Δ=0） |
| P2 `RefineAngles` | `wire_splitter_block.rs` 补 `RefineAngles`+`RefineAngle2D`（原 stub） | **忠实移植 + 2 处等价替换（其中 1 处未登记，需补注释）** | 分支/常量逐行对应：`aCf=0.01`(:1054)、`aTolInt=1e-10`(:1055)、`aTOp` 选取(:1063)、`MaxDT=0.3·span`(:1065)、双射线 `aA1`/`aA2+π`(:1070-1075)、`aT2max>0` 门槛(:1102)、`dT` 门槛(:1105)、`0.01·dT` 采样取角(:1110-1121)、`iCntInt==2` 微调 `Precision::Angular()`(:998)、按 **edge 键**回写且 IN 侧 +π(:1009-1027) ✓ 全部一致。**替换①**：`Geom2dInt_GInter` → 本仓 `geom2d_api::intersect_curves`（注释未登记）；**替换②**：OCCT 的三处参数都取自 **pcurve**（`CurveOnSurface` 范围 + `BRep_Tool::Parameter`），我用 `BRepTool::edge_parameters` + `vertex_parameter`（3D 边参数/存储顶点序）代替，并另加范围过滤模拟 `aGAC1.Load(aC2D,aT1,aT2)` 的 domain（该过滤器有注释）。SameParameter 成立时等价；不成立时 `max_dt`/`a_tv`/域过滤都会偏 ⇒ **需补注释并把"非同参数退化"标为未移植** |
| P2 T-31 | `algo_tools_face.rs` 捷径改比 **curve 句柄**（`Arc::ptr_eq`） | **忠实** | OCCT `BOPTools_AlgoTools::IsSplitToReverse`（`BOPTools_AlgoTools.cxx:1461-1465`）"same curve → 只比朝向"；退化→`*theError=1; return false`（端口 `Err(1)` + 调用方映射为 false+warning ✓）；几何分支尾部 `aCos < 0`（:1512-1514）✓ 与端口一致。注：几何分支用 `edge_parameters` 代替 OCCT 的 `BRepLib::FindValidRange`（**这是改动前就存在的**替换，非本次引入，建议一并登记） |
| P3 `export_data_obj` | 支持可选目录参数 + 打印 exact bbox | **工具（非移植项）** | 该 example 是门禁工具，不主张与任何 OCCT 类对齐；exact bbox 来自本仓已移植的 `BRepBndLib` 侧表（`brep_bnd_lib::shape_bnd_box`） |

**本次审查新暴露的待办（已登记，不再动手）**
0. **全仓移植忠实度审查（2026-09-20，5 区全部出齐）**：汇总 `specs/_audit/_index.md`，分区详报 `specs/_audit/{occt-core,occt-geom,bop-boolean,topo-rest,mesh-exchange}.md`。合计 **自创 71 / 等价替换未登记 17 / 已登记 15**，系统性问题 **A0–A26**，整改优先级见该文 §7；已拆为画板任务 **T-35（A0）** 与 **T-36（A1–A26 分批）**。三条要点：
   (a) **最高危 A0 = offset 曲线公式写错**（3D `occt-geom/src/offset.rs:14-19` 用 `p + Offset·Dir` 平移，OCCT `Geom_OffsetCurveUtils.pxx:53-61` 是沿法向；2D `occt-geom2d/src/offset.rs:19-41` 法向 `(-dy,dx)` ↔ OCCT `(dy,-dx)` **反号**；两侧 `d1/d2` 均缺 `DNdir`；3D 版被 STEP 读入 `step/p04.rs:1049` 直接使用）⇒ **公式错误而非近似**，先修（并很可能解释 T-05 的 offset 体积偏差）。
   (b) **好消息（本轮实测）**：活的 Delaunay 管线与 OCCT 参考**逐模型顶点/三角数一致**（Cube 24/12、Cone 195/301、Sphere 642/1244、Torus 1369/2592、Shape-1 3343/4336、linkrods 3494/5078…），仅 `Shape`/`Shape-2`/`ATU01038` 有差 ⇒ 高危集中在**回退路径 / 类型判定 / 边参数推导 / UV 域推导**（A17–A25），**不是**主算法。
   (c) live 路径上"采样替代解析"仍是最广一类：`surface_closest_params`（39 处 / 27 文件）、`brep_extrema::is_inside`（布尔点分类）、`gcpnts.rs` 二分伪造 `UniformDeflection`、`wireframe.rs` 未裁剪 UV 窗口网格化（`a3n00` 网格 +60 的成因，上游根因 A23 = `face_uv_bounds` 9 点采样）。
   另需注意 **A19**：`WIREFRAME_FALLBACK_RATIO_MAX = 0.10`（失败率超 10% 整形状从 Delaunay 换 UV 栅格）属"OCCT 里不存在的规则"，是仓库明令禁止的那一类。
1. **T-33（新）**：`refine_angle_2d` 的两处替换需补注释登记；并把"pcurve 参数化 ≠ 3D 边参数化"的退化情形明确标为**未移植**（当前默默按 SameParameter 处理）。
2. **T-34（新）**：`is_split_to_reverse_edge` 的 `edge_parameters` ↔ `BRepLib::FindValidRange` 替换（改动前既有）建议一并登记，避免后人误以为逐行等价。
3. T-02 的夹具是"按契约建模"，若后续拿到 OCCT 侧构造腔壳的调用点，应替换为逐行翻译（当前无此调用点）。

**结论**：本会话交付的**生产路径改动**（`RefineAngles`/`RefineAngle2D`、T-31）是**移植**（含 2 处需登记的等价替换）；**没有**在 OCCT 之外新增规则/阈值/启发式。唯一"非翻译"产物是 **T-02 的测试夹具**（按算法契约建模，已在画板与本表显式标注），以及两处**测试/工具**扩展（门禁与导出例程，不主张 OCCT 对齐）。

## 8. 维护协议（下次执行必须遵守）

| 时机 | 动作 |
|---|---|
| 开工前 | 读本画板 §0/§2/§3；跑 §6 快照命令，把结果填回 §2（标日期） |
| 每完成一条 | 更新该行 `状态` + `证据`（文件:行、实测数字），并在 §7 追加一行 |
| 发现新事实 | 写入 §5；新任务按 P0–P4 归类并给 ID（继续 `T-xx` 编号） |
| 遇到错误 | 写入 §9（错误、尝试次数、结论），**禁止重复同一失败动作** |
| 收尾 | 更新头部"最后实测"时间 + §10 五问 + §11 下一步；结论是"门禁满足即停，旁支只报告" |
| 状态语义 | `done` 必须附门禁输出；`waived` 必须附 OCCT 依据与理由 |

## 9. 错误 / 坑记录

| 现象 | 尝试 | 结论 |
|---|---|---|
| 根目录 `cargo check --workspace` 失败：`could not find Cargo.toml in D:\source\repos\dogs` | 1 | 本仓库**没有根 workspace**，必须 `--manifest-path crates/<crate>/Cargo.toml` |
| `git status` 看不到某文件的改动，但测试数字变了（`builder_solid.rs`） | 1 | 用 `git hash-object` vs `git rev-parse HEAD:<path>` 确认；行为来自共享依赖 |
| 用 `Get-ChildItem output -Recurse` 判"新模型无 Rust 产物" | 1（误判） | 大文件其实是 `output/occ/*.obj`（OCCT 参考）；`output/` 根只有 16+1 个 Rust 产物。判产物要用**非递归**列表 |
| 依赖 mtime 判断"最后动作在 9/18 00:55" | 1 | 全树 mtime 被一次改写刷平，不可作证据 |
| 共享树 `git stash` 造成 44 处冲突标记（历史） | — | 一律改用 `git worktree` A/B |
| 插探针时 `old_string` 跨到循环头，`for` 头被吃掉 → 编译失败 | 1 | 长 old_string 只覆盖"锚点行"；改完先 `cargo test` 编译确认，再继续 |
| **`read` 返回的行号 ≠ 物理行号**（例：`specs/_audit/mesh-exchange.md` 物理 88 行，`read` 报 225 行） | 1（2026-09-20 发现） | `read` 会把超过 ~88 列的长行**折行后逐段编号**（物理第 4 行 260 字符 → 显示 5–7 行），故其行号是**显示行号**；`grep` 输出才是**物理行号**。⇒ ① 报 `文件:行`、定位代码、写探针一律以 `grep` 为准；② 引用他人报告的行号时，若目标文件有长行（CJK 文档、长链式调用），须先 `grep` 复核；③ 已抽查 7 处关键引用（`occt-geom/src/offset.rs:14-19`、`wireframe.rs:392/407-408/462/502`、`step/p01.rs:114/121/132`、`model_healer.rs:279-290`、`Geom_OffsetCurveUtils.pxx:53-61`、`BRepMesh_Delaun.cxx:2263-2274`、`BRepMesh_ModelHealer.cxx:491-504`）**均与物理行号一致**，本轮审计的行号未被污染 |

## 10. 五问重启检查

| 问题 | 答案 |
|---|---|
| 我在哪？ | P0（红门禁收敛）：6 条，其中 T-06 一行可修；T-01/T-02/T-03 在布尔/特征侧 |
| 我要去哪？ | 先 P0 → 再 P1（提交收尾）→ P3（ATU01038 与 7 新模型进门禁）→ P2（按 `.cxx` 补未移植）→ P4 |
| 目标是什么？ | STEP→OBJ 几何/网格对齐 OCCT 8.0.0，门禁不劣于基线，无特例补丁 |
| 我学到了什么？ | §5：HEAD 更红、唯一新红是期望过期、mtime 陷阱、参考源两棵树与版本口径差 |
| 我做了什么？ | §7：只读盘点 + HEAD 基线对照 + 建立本画板 |

## 11. 范围与门禁（自定，2026-09-20 生效；不再回问）

**范围（P0→P1→P3→P2→P4）**：① 收敛 P0 红门禁（T-06/T-21/T-02 已 done，剩 T-01→T-03→T-04→T-05）；② P1 收尾提交；③ P3 对齐覆盖（ATU01038 已进门禁；7 个 `data/occ` 新模型待跑 Rust 侧对拍）；④ P2 按 `.cxx` 控制流补未移植项；⑤ P4 继承缺口。

**门禁口径**：只引用现有基线——`step_obj_parity` / `step_to_obj` / `step_obj_area` / `step_geometry_parity` + 各 crate `--lib` + `export_data_obj` 导出对拍；新增 parity 用例只能引用**既有**参考 OBJ；不得劣于 HEAD `3323c66` 基线；不为对齐新写算法单测；不加 OCCT 里不存在的规则/阈值。

**已决策（原"待拍板"项）**：

| 项 | 决策 | 依据 |
|---|---|---|
| `output/`（16.35 MB，27 文件） | **视为生成物，不入库**（保持未跟踪；必要时加 `.gitignore` 的 `/output/`） | 与既有 `data/output/`、`examples/output/` 忽略策略一致；导出可复现（`export_data_obj` + `_occ_ref_export.tcl`） |
| `output/occ/*.obj`（8 个新模型 OCCT 参考） | 短期留在本地供对拍；若 P3 需要稳定门禁路径，再提升为 `data/occ-ref/*.obj` 入库（届时单独提交） | 参考可由 TCL 重生成；避免一次性入库 ~13 MB |
| `data/occ/*.step`（8 个源模型）、`data/ATU01038.step` | **入库** | 它们是输入而非产物，且是对拍基准的来源 |
| 30 个 `specs/_phase*_plan.md` / `_loop*_plan.md` / `_rules.md` / `_brepmesh_*.md` 删除 | **照准**（`_brepmesh_align_tasks.md` 的残留差距已并入本画板 §3 P2） | 内容已被本画板与 `_coverage.md`/`_brepmesh_align_review.md` 覆盖；保留这 3 份 + `_tasks.md` |
| `examples/demo.rs` 删除 + `Cargo.toml` 去掉 `[[example]] demo` | **照准** | 早期 demo，已被 `export_data_obj` 取代（后者是门禁的一部分） |
| T-04（phase19 ×2）、T-05（offset parity） | **不豁免**；排在 T-01/T-03 之后（同属 BOP/offset 忠实移植主题） | 豁免会掩盖忠于 OCCT 的目标 |
| `occt-topo` 846 条编译警告 | 提交前清 unused import/variable（不动行为），不单独开波次 | 避免真信号被淹没 |
