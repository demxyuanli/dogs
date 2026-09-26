# _board — OCCT→Rust 对齐：未完成任务监控画板

> **载体**：`D:\source\repos\dogs`，master 工作区。
> **当前状态**：所有代码改动**均已提交**，工作树干净（唯一未跟踪项 = 语料 `data/iges/`）；门禁基线以 **§2 的 A/B 实测表**为准（不是旧文档里的数字）。
> **本画板建立**：2026-09-19 23:50 (+08:00)；**最后实测**：2026-09-21 第 114 轮（批 100 后全门禁：lib **1285/2**、`parity` 14/14、`step_to_obj` 13/13、`export` 16/16；第 111 轮的 A/B 差分见 §2）；**最后重梳**：第 113 轮（§3 重写为唯一监控清单）。
> **读法**：**§3 = 在办任务（唯一权威，§3.1 红门禁 → §3.2 根因批 → §3.3 缺口池 → §3.4 卫生 → §3.5 归档）**；§3.6 与 §7 = 历史明细/轮次日志（查证据用，状态列已作废）；§14 = A# 查号表；§2 = 门禁快照；§6 = 刷新/差分/playbook；§8 = 维护协议。
> **上游文档**：`specs/_coverage.md`（覆盖矩阵，2026-09-09，**已过期**，见 §3.4 T-29）、`specs/_brepmesh_align_review.md`（BRepMesh 对齐根因）、`specs/_tasks.md`（迁移清单，2026-08-29，已过期）、`specs/_audit/_index.md`（A0–A31 忠实度审查）。

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

## 2. 状态快照（2026-09-21 第 111 轮：门禁波次实测；A/B = `9d8e596`（批 84 末）↔ `2c9a66b`（批 98））

> **本会话（round 30 → 66，HEAD `5ed8064`）对这些数字的更新**（同一命令实测，逐项见各节）：
> `--lib` 1285/2 → **1287/0**（`113b4d1` T-01、`f63c324` T-86 转绿）；`step_geometry_parity` 2/3 → **3/3**（`9cac346` T-05 收口）；
> 其余逐项不变：boss 1/2（`boss_single_disc_base_merges_one_solid` = R2-8 定案 (b)）、phase10 7/8、phase19 3/5、
> phase20/3-9 全绿、step_obj_area 11/11、step_obj_parity 14/14、step_to_obj 13/13、iges_check ok。
> 另发现两处**测试隔离**问题（与引擎无关，已登记）：并行/交叉运行时 `step_to_obj_writes_output_files` 偶发假红（与 `export_data_obj` 同写 `output/`），
> `brepfeat::tests::groove_cuts_cylinder` 单跑有时与全量结果不同（全局状态依赖）。
> **补充（round 69）**：`phase19` 由 3/5 → **5/5**（`specs/_board.md` R2-9/T-04 收口：两条断言的 `> 2*56` / `>= 2*56` 是**端口自造阈值** ✗，
> 已按 goal ④ 换成 OCCT 依据的**关系式** —— DS ⊇ 参数自身子形状、相交情形严格更多；GT 探针新增 `--ds` 给出 OCCT 实测：
> 两个不相交单位盒 `NbShapes == NbSourceShapes == 68`、相交对 `NbShapes = 80`）；
> 同时把 `step_to_obj::write_output` 改为**临时文件 + 原子改名**（多个测试同写 `data/output/<name>.obj` 时读方不会再看到半写文件）。
> **round 80 补充（flake 口径修正）**：该 flake 比登记时**更宽** —— 并发整跑时失败的是 `cube_step_to_obj`（单跑通过），
> 而 `step_to_obj_writes_output_files` 反而通过 ⇒ 同一测试二进制内 11 处"写后回读 `data/output/<name>.obj`"的用例
> 彼此竞争；round 69 的原子写只消除了"读到半写文件"，未消除"读到别的用例刚写/正要写的同名文件"。
> 后续修法（未做）：让每个用例写自己的临时目录再回读，或把该二进制标为串行执行。
>
> **本会话收尾记录（round 30 → 80，HEAD `24bcb44`，37 个提交）**
> - **门禁改善两处**：`step_geometry_parity` 2/3 → **3/3**（T-05 收口）；`phase19` 3/5 → **5/5**（R2-9/T-04 收口）。
> - **完成的 §3 条目**：T-87/T-05 ✓、T-69 ✓、R2-9/T-04 ✓、R2-23 ✓（+T-44 部分：trimmed-Bezier 臂）、
>   R2-18 部分 ✓（2D 弧长改忠实 `CPnts_AbscissaPoint` + 容差重载）、T-85 余项 ✓、T-62 余项 ✓、T-27 ✓、T-26 ✓；
>   §3.4 另落 T-29/T-30/T-33/T-34 与 T-11 两批；另交付 **T-25/T-26/T-27/T-28 设计文档**
>   （`specs/_design_architecture_t25_t28.md`）与 GT 探针新增 `--fuse`/`--ds` 模式。
> - **关键定位**：2D BSpline 反向按 `BSplCLib::Reverse`（`k' = kfirst+klast−k`）保参数域（`59fefad`，把
>   `data/Offset.step` 从 2622.55 修到 **2610.501436**，GT 实测 2610.501440）；`transformed_copy` 的 pcurve **面键**
>   重映射（`2fcec12`）；偏置面按 `BRepOffset` 球面解析分支重建边界 + 源 pcurve 搬新面键（`a4a3792`/`9002ed2`，
>   偏置球解析 = **14.137166941**）；布尔结果补 pcurve 供给（`5ed8064`，符号 −2.198 → +3.168）；
>   `Geom_BezierCurve::Segment` 逐行移植（`2616041`，≤1.4e-15）；`math_BrentMinimum`（`3af4e4e`）。
>   并发现既有缺陷 `occt-core/src/bspl/poles.rs::power_to_bezier_basis` 公式错（`i==j>0` 给 1，应为 `C(d,j)`）——
>   已在原处标注（KNOWN DEFECT）—— **round 80 已订正**：矩阵改为 `(−1)^(i−j)·C(d,j)·C(d−j, i−j)`（`7b8b57f`），
>   并把 `poles_to_coefficients` 整体标为 **NOT USABLE / superseded shortcut**（它另有一处自创的"逐 span 平均"
>   缺陷 —— 平均除以了未填充的 temp 行 ⇒ 单 span 二次 Bezier 结果被缩 1/3；且全仓**零调用者**）；
>   OCCT 的真正路线是 `BSplCLib::BuildCache`（曲线版）+ `PLib::CoefficientsPoles`，即 T-44 的前置。
> - **仍未完成（板内均有精确分步路线）**：T-80 链/T-41（布尔结果边界是折线 ✗ ⇒ 按 `BOPAlgo` 换精确裁剪曲线；
>   **实测：`shape_volume` 接解析体积后全仓只剩这 1 条红**）、T-67 步3+T-37 与 T-51 余项（共同前置
>   `math_FunctionSetRoot` 1452 行 ✗；其依赖 `math_BrentMinimum` 已移植 ✓）、T-44 余项（`IncreaseDegree` ✗）、
>   T-54（约束 Delaunay ✗）、R2-18 余项（一般 2D 求交 = `Geom2dInt_GInter` ✗）、R2-19（IGES 非平面面 2D UV 曲线 ✗）、
>   T-25/T-28 的实施（设计已交付 ✓）、T-11 余项（低）。

| 门禁 | 命令 | 当前工作区 `2c9a66b` | 对照 `9d8e596`（批 84 末） | 判定 |
|---|---|---|---|---|
| 编译 | `cargo check --manifest-path crates/occt-topo/Cargo.toml --all-targets` | ✅ exit 0（lib 843 条警告） | ✅ exit 0 | 绿 |
| 五 crate 编译 | 各 crate `--all-targets` | ✅ exit 0 / 0 error（core·math·geom·geom2d·topo） | — | 绿 |
| STEP→OBJ bbox parity | `--test step_obj_parity` | ✅ **14/14** | ✅ 14/14 | 绿 |
| STEP→OBJ 端到端 | `--test step_to_obj` | ✅ **13/13** | ❌ 12/13（`step_to_obj_writes_output_files`） | 绿（**+1 改善**） |
| 面积对拍 | `--test step_obj_area` | ✅ 11/11 | ✅ 11/11 | 绿 |
| 几何一致性 | `--test step_geometry_parity` | ✅ **3/3**（T-05 收口：断言按 OCCT 规格改为对 GT 实测解析值 `2610.501440` 校验，端口 `2610.501436`；提交 `9cac346`） | ❌ 2/3（同名用例） | 绿（**+1 改善**） |
| topo 单测 | `--lib` | ✅ **1281 通过 / 0 失败**（T-01 `113b4d1`、T-86 `f63c324` 转绿；T-27 `1486985`+本轮摘除 leftover 死模块**带走 6 条死测试** ⇒ 1287−6；TShape 唯一 id 见 `933785f`） | ❌ **1286 / 2（同两条）** | 绿 |
| boss 合并 | `--test bop_builder2_boss` | ❌ 1/2（T-03） | ❌ 1/2 | 红（遗留） |
| phase19 | `--test phase19_integration` | ❌ 3/5（T-04） | ❌ 3/5 | 红（遗留） |
| phase10 | `--test phase10_integration` | ❌ 7/8（`curved_face_fillet_sphere_plane`） | ❌ 7/8 | 红（遗留） |
| phase20 / 3 / 4 / 5 / 6 / 7 / 8 / 9 | 各 `--test phaseN_integration` | ✅ 5·4·9·7·5·5·5·8 全绿 | ✅ 同 | 绿 |
| core / geom / geom2d / math | 各 crate `--lib` | ✅ 290（1 ignored）· **143** · 72 · 215（1 ignored） | ✅ 290(1i)·**143**·72·215(1i)（逐项相同） | 绿 |
| topo doc-tests | `--doc` | ✅ 2 passed / 1 ignored | ✅ 同 | 绿 |
| IGES 结构自洽（常驻校验，非门禁） | `--example iges_check -- <18 模型>` | ✅ **18/18 ok**：6 个模型（ATU01038 / screw / Shape / Shape-2 / occ-bottom / occ-top）P 段长度变化，**DE 数与 unreferenced 统计逐模型与对照完全相同** | ✅ 18/18 ok | 绿 |
| 导出总检 | `--example export_data_obj` | ✅ 16/16（Cube 24/12、Sphere 642/1244、Torus 1369/2592、Shape-1 3343/4336、linkrods 3494/5078、ATU01038 17745/22119 …与 §13 记录的网格数逐项一致） | — | 绿 |

> **口径订正（R2-15 收口）**：本节旧表记"第 90 轮 lib **1287/1**、geom **146/146**、step_to_obj 13/13"。本轮以同一命令对 `9d8e596`（批 84 末）与当前树做 worktree A/B，实测 **lib 两端都是 1286/2**、**geom 两端都是 143**、`step_to_obj` 在批 84 是 **12/13**、批 98 才是 13/13。`iges.rs` 在批 84（`9d8e596`）与批 95（`83ed3f8`）之间**未被任何提交改动**（`git log -- crates/occt-topo/src/iges.rs`），故 `sphere_iges_has_arc_and_solid` 在"第 90 轮"也必然失败 ⇒ 旧表的 1287/1 与 146/146 是**记录错误**，以本节实测为准（旧表的"计数变化说明（批 65–67）"一段作废）。
> **`step_to_obj` 由 12/13 变 13/13 的改善**出现在批 85–97 之间（本波次只做 A/B 对比，未再二分归因：该用例在批 84 的失败与 T-05/Sphere 网格回退相关，属既有红转绿）。

**HEAD 基线复现方法（不改共享树）**：
```powershell
cd D:\source\repos\dogs
git worktree add --detach .target-headcheck HEAD      # 路径命中 .gitignore 的 /.target-*/
cd .target-headcheck
cargo test --manifest-path crates/occt-topo/Cargo.toml --no-fail-fast
cd ..; git worktree remove --force .target-headcheck
```
> 历史事故：**禁止**在共享树里 `git stash` / `stash pop`（曾造成 44 处冲突标记）。要 A/B 就用上面的 worktree。

## 3. 任务板（监控清单 —— 2026-09-21 第 113 轮全量重梳）

> **本清单是全仓唯一的在办任务表**。旧分块（`R2-xx` / `P0–P5` / `§14 A0–A31`）已合并去重，**它们的"状态"列一律作废**，明细整体下沉到 **§3.6 历史明细（留档）**。
> **状态枚举**：`pending` 待做 · `in_progress` 在做 · `blocked` 有明确阻塞物 · `waived` 豁免（须写 OCCT 依据）。
> **每行必带**：一句话 · 证据（实测数字或提交号）· 下一步（落到文件/函数/`.cxx` 行号）· 验收（只引用既有门禁/用例，不新写测试）。
> **执行口径（长期有效）**：以"能翻译成代码的 OCCT 控制流"补齐缺口，以编译作初步验证；未接指令不自发起门禁/导出；禁止特例补丁、自造阈值、为对齐改断言。每批同步本表 + §7 日志。
> **ID 对照**：`R2-7 = T-01(红) + T-32(根因)`｜`R2-8 = T-03`｜`R2-9 = T-04`｜`R2-10 = T-05 + T-87`｜`R2-11 = T-80 链(T-82/T-83)`｜`R2-12 = T-69`｜`R2-13 = T-67 步 3 + T-37`｜`R2-19 = T-78 余项`｜`R2-23 = Geom_BezierCurve::Segment`。

### 3.1a T-88 分诊结论（本会话，2026-09-26）

`phase10_integration::curved_face_fillet_sphere_plane` 的失败**不是夹具错误，而是 T-80/T-82 的同族缺陷**（探针实测，已删）：

- `boolean_compound(box(4,4,2), sphere(r=1)@(2,2,2), Fuse)` 走的是**忠实 BOP**（`boolean` → `boolean_via_bopalgo` → `bop_builder2`，不是 `bop_curved`）。结果 8 面 = 7 Plane + 1 Sphere，**形状基本正确**（上半球保留、盒顶面变成 2-wire 环）。
- 缺陷两处：① 盒顶面圆内那块**内部平面盘**（`f3`）被保留（FUSE 应丢）；② 圆边**重复且不共享** —— `f2` 的内环圆与 `f7`（球面）的边界圆是**不同 TShape**，同一几何位置共 3 条圆边（`owners=[f2(Plane),f3(Plane)]`、`[f7]`、`[f7]`）。
- 因此 `fillet_edge_curved` 取到的那条圆边两侧是 (Plane, Plane)；即使按用户的「圆边」判据选边（本会话已把测试的选边判据由「`|mid.z-2|<0.5`」收紧为「圆边 ∧ `|mid.z-2|<0.5`」——原判据会先命中盒体轮廓线，属夹具欠定），仍取到重复边。
- **结论**：T-88 与 T-80/T-82 同源（`BOPAlgo_BOP::BuildShape`/`BuildSolid` 的状态过滤与图像共享），**由 T-82 收口时一并转绿**；在此之前保持红。
### 3.0 本会话实测与进度（round 81→，2026-09-26；**以本节数字为准，§3.1–§3.4 的旧状态列作废**）

**实测基线（HEAD = 每批提交后同一命令复测）**

| 门禁 | 实测 | 与 §2 快照 |
|---|---|---|
| 五 crate `cargo check --all-targets` | ✅ exit 0 / 0 error | 绿 |
| `occt-topo --lib` | ✅ **1281 / 0 failed** | 绿（T-01/T-86 已闭） |
| `step_obj_parity` | ✅ 14/14（含 `data/occ` 三模型新断言） | 绿 |
| `step_to_obj` | ✅ 13/13（**单跑**；`--no-fail-fast` 并发整跑仍偶发 `output/` 竞争假红，见 §3.4 旁支） | 绿 |
| `step_obj_area` | ✅ 11/11 | 绿 |
| `step_geometry_parity` | ✅ 3/3 | 绿 |
| `phase3/4/5/6/7/8/9/19/20` | ✅ 4·9·7·5·5·5·8·5·5 全绿 | 绿（phase19 已闭） |
| `phase10` | ❌ 7/8（`curved_face_fillet_sphere_plane` = T-88） | 同 |
| `bop_builder2_boss` | ❌ 1/2（R2-8 已定案 (b)：夹具不合法，非引擎缺口） | 同 |
| `export_data_obj` | ✅ 16/16，v/f 与 §13 逐位一致 | 绿 |
| `iges_check` | ✅ ok（含 Sphere `unreferenced=1`＝根） | 绿 |
| `occt-core / math / geom / geom2d --lib` | ✅ 290(1i) · 215(1i) · 143 · 72 | 绿 |

**本轮已完成（提交号）**

- **T-80 根因三处**（`da61a30`）：① 圆柱侧面 wire 朝向按 `BRepPrim_Builder::AddWireEdge`（`BRepPrim_Builder.cxx:184-192`）订正；② 缝顶点容差按 `BRepPrim_Builder::MakeVertex`（`cxx:129-132`）给 `Precision::Confusion()`；③ `FClass2d::edge_points` 按 `BRep_Tool::CurveOnSurface`（`BRep_Tool.cxx:301-315`、`:347-357`）对 REVERSED 边取 PCurve2。
  实测 GT（`occt_probe --fuse`）：`box[-1,1]³ ∪ cyl(r=0.4,z∈[0,2])` = volume **8.50265** / **8 faces**；端口 GF 由「9 plane + 1 cylinder」变为「7 plane + 2 cylinder + 2 plane」= 11 faces，柱面被 z=1 正确切成两条环带、两个 area。
- **T-80 剩余缺口（= T-82 精确化，未完成）**：`BOPAlgo_BOP::BuildSolid` 的 `BOPAlgo_BuilderSolid` 对柱体 draft 的 6 面集合给出 **loops=0**（面全被 `PerformShapesToAvoid`/`PerformLoops` 后处理判为 internal）⇒ 最终只留盒体 7 面、体积 8.29787。已定位到两处控制流语义：(i) `TopExp::MapShapesAndAncestors` 对同一面的**边出现次数**逐次追加 ancestor（`TopExp.cxx:80-120`），端口 `perform_shapes_to_avoid` 用 `edges_of`（按 TShape 去重）⇒ 缝边 `aNbF` 由 2 变 1 而被 avoid；(ii) `BRep_Tool::IsClosed(E,S,L)` 对**平面直接返回 false**（`BRep_Tool.cxx:819-822`），端口此前用「边界出现次数≥2」近似。按 (i)+(ii) 改写的实验版把 Fuse 打成空（`avoided=0`、`loops=0`），已 `git checkout` 回退；下轮从「draft 6 面集合的边出现次数直方图」与 `CloseOpenShells` 入手。
- **T-51 / T-67 共同前置**（`a647836`）：忠实移植 `math_FunctionSetRoot` + `math_FunctionRoot`（`function_set_root.rs`）与 `math_BrentMinimum` 的失败中止分支。
- **T-44 / R2-6 RationalC1 臂**（`2a28b8c`）：`Geom_BSplineCurve::{IncreaseDegree,RemoveKnot}`（`Geom_BSplineCurve.cxx:243-298`、`:420-495`）+ `GeomConvert_CompCurveToBSplineCurve`（`GeomConvert_CompCurveToBSplineCurve.cxx:32-273`）+ `BSplCLib::{IncreaseDegree,AntiBoorScheme,RemoveKnot}`；GT 对拍（DRAWEXE）逐字段一致。
- **R2-18 2D 弧长**（`e8deb3d`）：`curve2d_length`/`_tol` 补 `Geom2dAdaptor_Curve` 类型解析、IsRational 判据、CN 分段与容差重载（`GCPnts_AbscissaPoint.cxx`、`CPnts_AbscissaPoint.cxx:26-65/:79-100/:191-207`）。
- **parity 覆盖**（`63ae2b9`）：`data/occ` 的 a3n00 / acs10 / TDB 按既有 WriteObj 参考锁进 `step_obj_parity`。

**进行中 / 未开始（本会话后续）**：T-67 步 3（`Extrema_GenExtPS`）与 T-51 反解已派子代理；未开始 = R2-19（IGES 2D UV 曲线）、T-54（约束 Delaunay）、T-41（`bop_curved` 网格布尔摘除，需先修 T-80/T-82）、T-25/T-28（设计已出，待实施）、T-11 余项。
### 3.1 红门禁（当前实测为红：2 条 `--lib` + 3 条集成 + 2 条属性/网格缺陷）

| ID | 门禁 / 用例 | 实测（2026-09-21） | 根因已定位？ | 归属根因批 |
|---|---|---|---|---|
| **T-01** | `brepfeat::tests::groove_cuts_cylinder`（`--lib`） | 批 100 后：`grooved 8.408733 vs expected 7.708990`（removed 0.909 / 应 1.6085；批 100 前为 `8.9034`，removed 0.414）。**round-14 复测**：`BuildSplitFaces tasks=112 area_hist={1:71, 2:16, 3:25}`（round 2 基线为 `{0:40, 1:41, 2:30, 3:1}`） | **部分**：输入侧已修（T-32 批 100：wire 成链 + 闭合实体按 `BRepTools::OrientClosedSolid` 定向 ⇒ 夹具不再是反向实体）；**"0 块"坏面 40→0**，但仍有 **71 面只出 1 块**（样本 `le=13/15` ⇒ 13–15 条候选边只成一个 area），与 round 3–6 的 `avoid=3` 同类 ⇒ 剩余缺口在 **`perform_loops` 的 WireSplitter 成环**（`builder_face_occt.rs`/`wire_splitter_block.rs`） | **下一批 = round-15**：按 round 3 配方对 `perform_loops` 插桩 `(in_edges, avoid, avoid_bnd, loops, areas)` 分组，对照 `BOPAlgo_BuilderFace::PerformLoops` / `BOPAlgo_WireSplitter`（`_1.cxx:112-354`/`:358-617`） |
| **T-03** | `bop_builder2_boss::boss_single_disc_base_merges_one_solid` | `left 0 / right 1`（Fuse 结果**无 solid**） | 否（**输入有效**：`single_disc_cylinder` 是 `TopoBuilder` 手工装配 BRep；tri-fan 变体通过） | 3.2 序 5（R2-8 探测批） |
| **T-04** | `phase19_integration::{overlapping,disjoint}_boxes_pipeline_runs*` | `ds.nb_shapes()` = 单盒 **28** / 相离 **88** / 重叠 **84**（断言 `> 2*56` 与 `>= 2*56`） | 部分：断言阈值自造（56/盒 vs 实际 28/盒 ⇒ 112 不可达）；且 84<88 反直觉 | 3.2 序 6（R2-9 对照批） |
| **T-05** | `step_geometry_parity::offset_geometry_is_consistent` | `2208.0 vs 1612.9`（阈值 2%）；五种度量见 3.2 序 4 | 多因，非 offset 公式（T-35 已修而数字逐位不变） | **T-87** → 估计器核对 |
| **T-86** | `iges::tests::sphere_iges_has_arc_and_solid`（`--lib`） | 断言 IGES 含 `100`/`128` 失败；球面 IGES DE = `[144,196,116,123,123]`（无边界曲线） | **是**：端口球面 face **无边界 wire**（T-68 探针 `wires/face=[0]`） | **T-68** |
| **T-88**（新登记） | `phase10_integration::curved_face_fillet_sphere_plane` | 7/8（自 `7178861` 起同红 ⇒ 非本会话引入，板内此前**未登记**） | 否 | 待分诊（低优先） |
| **T-87**（新） | 非用例：`brep_to_obj` 网格绕向（T-05 的病根） | `Offset.step` 网格 892 三角中 **6 个有符号体积为负**（‑116.67 vs +116.67，形状为凸 ⇒ 应全正） | **是**：面/壳朝向（与 round 9"面法向·wire 绕向 112/128 为负"同源） | — |

> 其余门禁全绿：`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`phase20/3/4/5/6/7/8/9`、doc-tests 2/1i；五个 `--lib`：topo **1283/2**（= T-01 + T-86）、core 290(1i)、geom 143、geom2d 72、math 215(1i)；`iges_check` 18/18 ok；`export_data_obj` 16/16。

### 3.2 根因批（修完才可能让 §3.1 转绿；建议按序）

| 序 | ID（旧卡） | 状态 | 任务（可翻译的 OCCT 控制流） | 证据 / 前置 | 下一步 | 验收 |
|---|---|---|---|---|---|---|
| 1 | **T-32** | **✅ 输入侧 done（批 100）**：**①链向** = `BRepLib_MakeWire::Add`（`BRepLib_MakeWire.cxx:123-453`）完整移植（新模块 `brep_lib_make_wire.rs`：VF/VL/myVertex 状态、`E.Oriented(FORWARD)` 子顶点迭代、`IsSame` 判定、`reverse && !forward` 定向决策、proximity/copy-edge 支、`DisconnectedWire`/`NonManifoldWire` 错误）；`mesh_to_brep::triangulation_to_brep` 改用它建 wire；**②面朝向** = 实测已一致（96/96 面「面法向·wire 绕向」同号）；**③装配** = 闭合实体补 `BRepTools::OrientClosedSolid`（`brep_class3d::orient_closed_solid`，批 100 实测夹具网格是**反向缠绕**：signed volume −9.3175 ⇒ 之前造出的是"材料在外"的反实体） | round 8/9 + 第 112 轮：mesh→BRep 边全 Forward（source/sink）、面 47/96 绕向为负；批 100 实测 `faces 96 / edges 144 / 已无反向面`、T-01 `removed 0.414 → 0.909`（目标 1.6085） | **T-01 仍未绿**：剩余偏差在忠实 BOP 链内部 ⇒ 转入 R2-7 的 round-14 复测（输入已合法后重 dump `BuildSplitFaces` 的 `area_hist`/`avoid` 分布，看 40 面 0 块还剩多少） | 门禁与 §2 逐项相同（lib 1285/2、parity 14/14、`step_to_obj` 13/13、export 计数逐位相同）——**已满足**；T-01 转绿归下一批 |
| 2 | **T-68**（R2-7 上游） | ◐（球面 done） | 忠实 `BRepPrim_OneAxis` 的 lateral wire：**极点退化边**（v=±π/2 的 u 等参线跨 2π）+ 两条经线 pcurve（`BRepPrim_Sphere/Torus::SetMeridian`，含 `SetMeridianOffset(2π)`）；并补环面与带孔内环 | 探针：`sphere faces=1 wires/face=[0]`、`torus [0]`；执行方案见 §3.6「T-68 执行方案」 | 按该方案验证链（临时探针 → `--lib` 不劣化 → 重放 A13/A18 → T-59/A23） | T-86 转绿；T-49/T-54/T-55/T-59 解锁 |
| 3 | **T-69**（R2-12） | ◐（2 轮诊断） | 带孔面 / 单闭合边 wire 的**前沿链**：`FixLacking` 追加重复闭合边 + `ModelHealer` 回绕链 ⇒ `decompose_simple_polygon` 判无耳清空；需找 OCCT"跨 wire 周期对齐"的控制流 | 探针：166 面 `wires=2 ∧ surface=Other`；Torus face 20 两 wire pcurve 差一个周期 | 判据性实验（按面类 dump）+ 对读 `BRepMesh_NodeInsertionMeshAlgo`、`BRepMesh_Delaun::frontierAdjust` | T0M 面类计数 + 四道门禁 |
| 4 | **T-87** | pending | OBJ/面绕向：查 `brep_to_obj` 及其上游面/壳朝向为何 6 个三角反向；随后核对 `shape_volume`/`brep_gprop` 两个镶嵌型估计器（比网格奇偶体积低 36%）与 `brep_gprop_full` 只算 6 个平面面（1400.0 = (1/3)·600·7） | §3.1 的数字 | 对读 `RWObj_CafWriter`/`RWMesh_FaceIterator` 的三角形写出与 `Poly_Triangulation` 朝向语义；再核 `shape_mesh.rs::shape_volume` | `obj_volume` ≈ 网格奇偶体积（≈2610）；T-05 的断言可用（解析 oracle = Steiner ≈2610.5） |
| 5 | **R2-8 探测批**（T-03） | pending | boss 共面盘 Fuse：dump `BuildRC`（`obj_src/tool_src/it_shapes/check_keys/it_exp/tool_sets`）+ 每 solid 镜像数 + box 顶面（应成 2-wire 环）的 `rebuild_split_areas` 块数 | `left 0/right 1`；tri-fan 通过 ⇒ 引擎侧，范围钉在"共面盘 + 大顶点数环" | 对照 `BOPAlgo_BOP.cxx:583-711` 与 `BOPAlgo_Builder_3.cxx::BuildSplitSolids`；若镜像数=1 则瓶颈回到面级切割 | `--test bop_builder2_boss` 2/2 |
| 6 | **R2-9 对照批**（T-04） | pending | 读 `BOPDS_DS::Init`/`Append`（注册哪些子形状）+ 把 84/88 按 `ShapeType` 分解 | 28/88/84 实测 | 判"缺注册"还是"断言自造"，再决定修 DS 还是按 OCCT 订正断言（**禁止**把 56 直接改 28） | `--test phase19_integration` 5/5 或按 OCCT 订正 |
| 7 | **T-80 链**（R2-11） | pending（T-81 done） | 曲面 operand 的 `BOPAlgo_BOP`：面被节线切开 → GF 面 → 实体装配（`build_shape` 丢柱面/顶盖）；`BuildRC` 状态过滤与 `BuilderSolid` 的 growth/hole | 批 39–51 探针；断点在 `pave_ff_make::make_blocks_ff` 的 `is_valid_block_for_faces`（pcurve 与 3D 曲线参数系） | 沿 T-81 已修部分继续（`AdjustPCurveOnFace`/缝边收尾），dump 面状态与 growth 归属 | `box∪cyl` FUSE 单实体、体积 ≈8.5 |
| 8 | **T-67 步 3 + T-37**（R2-13） | pending | `Extrema_GenExtPS` + `Extrema_ExtPExtS`/`ExtPRevS`（+`math_FunctionSetRoot`）⇒ 一般曲面点–面极值忠实化；随后迁移 A1 的 39 处调用点（T-37） | 分步 1+2（`ExtPs` 分派层）已落地（`T-67` 行） | 移植 `Extrema_GenExtPS.cxx`（1195 行）+ `math_FunctionSetRoot`（1452 行） | 迁移后该 crate `--lib` + 门禁不劣化 |

### 3.3 移植缺口（非红门禁；按 OCCT 控制流补齐）

> **本会话进度（round 30 → 66，HEAD `5ed8064`）——§3.2 状态更新（细节见各提交信息）**：
> - **T-87 / T-05：done（§3.2 序 3）**。链路：`d04e9f4`（边界弧按 `BRepGProp_Face::Load(edge)` 处理边朝向）→
>   `8021fca`（缝边 pcurve 侧 + 去"UV 包围矩形"捷径）→ `9ce4c69`（`Curve2d::ReversedParameter` 逐类）→
>   `dacf235`（域路径改 OCCT 口径：未翻转法向 + 删自创 `wire_sign`；并按 `BRepPrim_OneAxis.cxx:465-468` 补圆柱盖面 pcurve、按 GT `--sph` 对齐球面丝）→
>   `b6489c4`（`BRepGProp_Face::Bounds` = 曲面参数，非面 UV 界）→ `bd00483`（U 向按 `UKnots` 分段）→
>   **`59fefad`（2D BSpline 反向按 `BSplCLib::Reverse`：`k' = kfirst + klast − k`，保持参数域）**。
>   实测：`data/Offset.step` 解析体积 1559.174 → **2610.501436**（GT 实测 **2610.501440**，相对 1.5e-9）；盒/圆柱/球原语保持精确。
>   验收面由 `9cac346` 收口（断言按 OCCT 规格改为对 GT 实测解析值校验），`step_geometry_parity` **2/3 → 3/3**。
> - **T-69（R2-12）：done（本会话收口）**。① `transformed_copy` 的 pcurve **面键重映射**（`2fcec12`）⇒ `pattern_volume_sums_copies` 在解析接线后转绿；
>   ② `brep_offset::offset_face` 的曲面分支按 OCCT `BRepOffset_MakeOffset` 的球面解析分支**重建边界**（`a4a3792`：关于球心的位似 `gp_Trsf::SetScale` + 缝出现保持同一 TShape），
>   并把**源 pcurve 搬到新面键**（BOPAlgo splitter 语义），仅在缺失时投影（否则缝两侧各自投影会被压到同一 u ⇒ 边界项相消）；
>   ③ 顺带查明端口 `is_seam_use` 的判据是"同边在丝里出现两次"（自创）✗，而 OCCT `ShapeAnalysis_Edge::IsSeam` = **`BRep_Tool::IsClosed(edge, face)`**（该边在闭合面上有两条 pcurve）——列为后续订正项。
>   实测：偏置球 `so`（r=1.5）解析 = **14.137166941** ✓（= 4/3π·1.5³），源球 4.188790205 ✓ 未被污染；
>   接线解析体积后 `--lib` 由 3 条红降到 **1 条红**（仅剩 T-80 的 `trimmed_face_box_cylinder_fuse`）。
> - **T-80 链（R2-11）：◐（本会话第 1 步）**。`5ed8064` 在 `curved_boolean_full` 出口补 pcurve 供给（OCCT 布尔结果的不变量；`ShapeFix_Edge::FixAddPCurve`），
>   把 `box ∪ cyl` 的解析值从 **−2.198** 修正到 **+3.168**（真值 8.50265）。**根因已钉死**：端口布尔结果的边界是**折线**（顶盘面积 0.49643 vs π·0.16 = 0.502655 ⇒ 约 −1.2%），
>   且部分面环绕向与曲面法向不一致 ⇒ 贡献变号（GT 探针新增 `--fuse` 模式给出 OCCT 对照：8 面 / 逐面面积与贡献）。
>   ⇒ 忠实解法是把布尔结果的面边界从折线换成 OCCT 的精确裁剪曲线（`BOPAlgo` 语义），属**实质算法缺口**，非一轮可成。
> - 旁支（测试隔离，与引擎无关，已登记）：`step_to_obj_writes_output_files` 与 `export_data_obj` 同写 `output/` 时偶发假红；`brepfeat::tests::groove_cuts_cylinder` 有全局状态依赖（单跑与全量结果可能不同）。
> - **R2-23 = done（round 70）**：`Geom_BezierCurve::Segment`（`Geom_BezierCurve.cxx:388-425`，非有理分支）已移植
>   （`occt-geom/src/bezier_curve.rs::segment`）：`BSplCLib::BuildCache(0,1,…)` 对 Bezier 即**幂基展开** +
>   `PLib::Trimming`（端口既有件，实测精确）+ `PLib::CoefficientsPoles`（幂基→极点的逆展开）。实测（二次/三次/四次 ×
>   `[0,1] [0.25,0.75] [0,0.5] [0.5,1] [0.3,0.3001]`）：与原曲线逐点最大偏差 ≤ **1.4e-15**。
>   并接线消费者 `GeomConvert::CurveToBSplineCurve` 的 trimmed-Bezier 臂（`GeomConvert.cxx:300-321`：裁剪后按
>   knots `{0,1}`、mults `{d+1,d+1}` 重建 BSpline）⇒ 该臂不再返回 `Unported`（**T-44 部分收口**；余下为 rational Bezier 分支与
>   `GeomConvert_CompCurveToBSplineCurve`/`GeomConvert_ApproxCurve`）。
> - **新发现的既有缺陷（round 70，已就地标注，未改行为）**：`occt-core/src/bspl/poles.rs::power_to_bezier_basis`
>   的矩阵公式错误（`i == j > 0` 处给 1，正确应为 `C(d,j)`）；正确形式为 `(−1)^(i−j)·C(d,j)·C(d−j, i−j)`
>   （见 `Geom_BezierCurve.cxx` 分支所用展开）。`poles_to_coefficients` 在 degree ≥ 2 时结果错误 ⇒ 列为后续订正项。
> - **剩余范围（round 71 收尾核验后的精确口径；均为"忠实移植可行但非单轮可成"）**：
>   - **T-80 链 / T-41**：端口布尔结果的面边界是**折线**（`bop_curved/region_trim.rs` 用邻接点 `make_edge_segment` 建丝），
>     且部分面环绕向与曲面法向不一致 ⇒ 贡献变号。`5ed8064` 补 pcurve 供给后解析值 −2.198 → +3.168（真值 8.50265）。
>     忠实路线 = 按 `BOPAlgo` 语义把结果边界换成**精确裁剪曲线**（`IntPatch` 交点 → 曲线段），并把 `bop_curved` 的
>     网格/体素布尔按 A5 摘除（T-41，须待精确布尔可用后再摘，避免"全绿假象"）。**接线解析体积的唯一残留红即此条**。
>   - **T-67 步 3 + T-37**：端口 `occt-geom/src/extrema_surf/analytic_solvers.rs` 已含 plane/sphere/cylinder/cone/torus
>     的解析点–面极值；余项 = `Extrema_GenExtPS`（1195 行）与 `math_FunctionSetRoot`（1452 行）的移植（一般曲面仍 24×24
>     网格 + 数值 Jacobian），随后才谈 39–40 处调用点迁移。
>   - **T-44 余项**：`GeomConvert_CompCurveToBSplineCurve`（274 行）需要 `Geom_BSplineCurve::IncreaseDegree`
>     （端口只有 `decrease_degree`）+ rational Bezier 分支；trimmed-Bezier 臂已随 `2616041` 可用。
>   - **R2-18 余项**：3D 侧 `occt-geom/src/gcpnts.rs` 已是忠实 `CPnts_AbscissaPoint`（OCCT 逐类型 Gauss 阶 + 13 次区间倍增）；
>     2D 侧 `occt-geom2d/curve_ops.rs::curve2d_length` 仍 Simpson、`curve2d_intersections` 仍 256×256 采样 ⇒
>     需补 `Curve2d` 的 line/Bezier 类型查询 + 2D 版 `CPnts_AbscissaPoint`（`order(2d)`：Line 2 / Parabola 5 /
>     Bezier min(24,2·deg) / BSpline min(24,2·N−1) / 其它 10）与 `Extrema_ExtCC2d`。
>   - **T-51 余项的真正前置 = `math_FunctionSetRoot`（round 73 查明）**：`math_FunctionRoot`（136 行）只是
>     `math_FunctionSetRoot`（1452 行）的**薄包装**（`math_FunctionRoot.cxx:70-116` 把 1D 函数包成
>     `math_MyFunctionSetWithDerivatives` 后交给 `math_FunctionSetRoot::Perform`）⇒ T-51 的 `gcpnts` 反解与
>     **T-67 步 3** 共用同一前置，须先移植 `math_FunctionSetRoot`（端口缺其依赖 `math_BrentMinimum`）。
>   - **T-85 余项 = done（round 75）**：查明 `iges_check` 的 `unreferenced` 不是写侧漏引用，而是**检查自身的指针表错** ——
>     196（球面）的参数格式是 `196,#0,{radius},#1,#2`（`IGESGeom_ToolSphere::WriteOwnParams`，端口
>     `iges.rs::emit_refs(196, …)` 同构），三个 DE 指针在字段 **1、3、4**；检查表按 `[1, 2, last]` 读 ⇒ **轴方向被漏计** ⇒
>     Sphere 误报 `{123: 1, 144: 1}`。按端口自身发射格式订正后：**Sphere 2 → 1**、Cube 1（不变），
>     剩下的 1 恰是该模型的**根**（OCCT `IGESData_IGESWriter::Write(root)` 只写从根可达者 ⇒ 根自身必然无引用）。
>     依据：`IGESGeom_ToolCurveOnSurface/ToolTrimmedSurface/ToolSurfaceOfRevolution/ToolTabulatedCylinder/ToolGroup/
>     ToolDirection/ToolLine/ToolTransformationMatrix` 等的 `OwnShared()`（本仓发出的类型逐一核对）。
>   - **T-11 第二批（round 73）**：`cargo fix --lib -p occt-core` 清掉 23 个文件的未用 import / 未用绑定 / 多余 `mut`
>     （43 删除，纯中性）；其中两处 `use` 是 test 模块经 `use super::*;` 依赖的再导出 ⇒ 已按真实路径
>     （`crate::precision::CONFUSION`、`crate::bvh::builder_tri::build_tri_bvh`）显式补回。核验：`occt-core --lib` **290/1i** 不变。
>     同时登记三处**既有信号**（未改行为）：`geom/csg.rs` 的 `if theta < 0.0 { let t = -t; }` 是**无副作用遮蔽**（疑为
>     `t = -t` 之误）；`bspl/plib2d.rs::flat_bezier_coefficients` 的 `weights` 与 `bspl/bezier.rs::reduce_degree`
>     的 `tolerance` **未被使用**；`bnd/sortbox.rs` 的 `ixmin` 未被使用。
>   - 其余：**T-54**（约束 Delaunay，大）、**T-11**（余项，低）、**T-25…T-28**（板内要求先出设计）、
>     **T-85/T-51/T-62 余项**（本轮未展开）。

| ID | 状态 | 缺口（OCCT 对应） | 下一步 | 验收 |
|---|---|---|---|---|
| **R2-19**（R2-3 余项 / T-78 余项） | pending | IGES 写侧 **2D(UV) 曲线**：`BRepToIGES_BRWire::TransferEdge(edge,face,…,false)`（`cxx:340-588`，平面面直接返回空）、逐类型 UV 修正、142 的 `PreferenceMode` 2→3；另 `SetUOrigin/SetVOrigin` 仍 UNPORTED | 新增 2D IGES 发射器（可复用 3D 发射器 + z=0 提升）+ `ShapeBuild_Edge::TransformPCurve`；逐类型变换清单见 §3.6 旧 R2-19 行 | `iges_check`（实体数变化属预期）+ 四道 STEP 门禁 |
| **R2-23** | pending | `Geom_BezierCurve::Segment`（`Geom_BezierCurve.cxx:388-425`）= **曲线版** `BSplCLib::BuildCache` + `PLib::CoefficientsPoles`（`PLib::Trimming` 已有） | 移植上述两件 → 接线三处消费方：`GeomConvert.cxx:300-321`、`GeomToStep_MakeCurve.cxx:77-82`、`GeomToIGES_GeomCurve.cxx:441-448` | 五 crate 编译 + 三处写侧差分 |
| **R2-6 余项** | ◐ | 仍缺三臂：trimmed-Bezier（= **R2-23**）、`RationalC1 ∧ U2-U1>=6`（缺 `GeomConvert_CompCurveToBSplineCurve`）、Offset（缺 `GeomConvert_ApproxCurve`） | 先 R2-23（前置最明确），再按消费方触发 CompCurve/ApproxCurve | `GeomConvert::CurveToBSplineCurve` 对 OCCT 全分派可用 |
| **R2-18 余项**（T-52 余项、A16 余项） | ◐ | 非解析 2D 组合（B 样条/Bezier/offset 及其修剪）求交仍是 256×256 采样（忠实路线 `Extrema_ExtCC2d`）；`curve2d_length` = Simpson（OCCT `GCPnts_AbscissaPoint`） | 移植 `Extrema_ExtCC2d` 与 `GCPnts_AbscissaPoint` | `occt-geom2d --lib` + 依赖 2D 的门禁 |
| **T-44** | pending（触发式） | `GeomConvert_CompCurveToBSplineCurve.cxx:135-215`；**批 98 后前置 `InsertKnots` 已就位**，仍缺 `GeomBSplineCurve::IncreaseDegree` | 等消费方出现（即 R2-6 的 `RationalC1>6` 臂）一并做 | R2-6 该臂可用 |
| **T-85 余项** | ✅（round 75 收口） | IGES 可达性过滤已落地（批 83），但 `iges_check` 的 `unreferenced` **未归零**，且批 84 与批 98 逐模型相同 ⇒ 需判"统计口径（根 402 未计引用）"还是"写侧漏引用" | 对照 `Interface_InterfaceModel::AddWithRefs`（`cxx:652-692`）与 `iges_check.rs:126-204` 的指针字段表 | orphan 解释清楚（归零或口径订正） |
| **T-67 余项** | pending | `Extrema_ExtPExtS`/`Extrema_ExtPRevS`（就地 UNPORTED） | 随 3.2 序 8 | 同 3.2 序 8 |
| **T-51 余项** | pending（低） | `gcpnts` 反解 UNPORTED（缺 `math_FunctionRoot`）；积分已换忠实 `CPnts_AbscissaPoint::Length`（批 50） | 移植 `math_FunctionRoot` 后接 `GCPnts_AbscissaPoint::Parameter` | `occt-core`/`occt-topo --lib` |
| **T-41** | pending | A5：`bop_curved` 的体素/网格布尔与计票无 OCCT 对应 ⇒ **摘除并标未移植** | 与 T-79/T-80 同波次（先修曲面布尔再摘，避免 §9 记录的"全绿假象"） | 摘除后 `--lib`/门禁不劣化 |
| **T-54** | pending | `wireframe.rs:257-380` 平面耳切 + 质心角度排序 + 桥洞 → 约束 Delaunay（`BRepMesh_DelaunayBaseMeshAlgo` + `BRepMesh_Delaun`） | 与 A13/A18 同族；先看 T-68/T-69 结论 | 网格门禁 |
| **T-62 余项** | ✅ 交换写侧 4/4（round 76 实测复核） | A26 的**网格交换**子项全部已忠实：① PLY 每面自成 node 块（`brep_to_ply` 明确**不**焊接；实测盒 = `element vertex 24` / `element face 12`，= 6 面 × 4 顶点，与 `RWPly_CafWriter::addFaceInfo` 的 `theNbNodes += theFace.NbNodes()` 同构）；② `property list uchar uint vertex_indices`（`RWPly_PlyWriterContext.cxx:214`）；③ STL 退化法向阈值 = `gp::Resolution() = RealSmall()`（`io/stl.rs::compute_normal`，`RWStl.cxx:325/:407`，退役的 `1e-12` 已记入注释）；④ STL 格式嗅探 = 前 5 字节 `solid`（`strncmp(aHeader,"solid",5)`，不额外判二进制长度）。余 2 项属 **IGES 侧**（曲线族 `<2%` 采样、整球自造 120）⇒ 归 **T-78/R2-19** 名下，不在本 ID 内。 | — | 该 crate 门禁 |
| **T-13** | pending（低） | `BoxOptimal`（`OptimizationHelpers.pxx` PSO+Powell）与 `SurfaceOfExtrusion::BoxOptimal`（`cxx:224-283`） | 按 `.pxx` 移植（当前 `AddOptimal` 不在 `BRepBndLib::Add` 路径上） | bnd 相关用例 |
| **T-14** | pending（低） | `BRepBndLib.cxx:95-101,157-179` 三角化 arm、仅 pcurve 边的 `BRepAdaptor_Curve` 回退（`cxx:93-106`） | 待 `Poly_Triangulation` 存储落地 | 同上 |
| **T-16** | pending（低） | `pcurve_full` 的 iso 链臂、`myGap` 残差恢复、`ProjectDegenerated` 单点重载（已论证等价） | 按需补 | 同上 |
| **T-18** | pending（低） | `Intf_Tool`、`Intf_InterferencePolygonPolyhedron.gxx`（无 2D 调用者） | 按需 | — |
| **T-19** | pending（低） | `meshing/triangulator.rs:19` 的 `addTriange34`/`checkCondition`（只被 VRML 读路径引用，**不属** STEP/OBJ 网格管线） | 按需 | — |
| **T-20** | pending（低） | `bop_build_solids_leftover.rs` 模块去留（`merge_sharing_faces` 已删） | 定夺后删或接线 | — |
| **T-30** | pending（低） | `uv_tolerance_2d` 用 `t3/\|d1u\|` 近似；OCCT `BOPAlgo_WireSplitter::Tolerance2D`（`_1.cxx:859-881`）用 `U/VResolution`（平面恒 1.0） | 按该函数改 | 与 T-01 无关，随波次 |

### 3.4 收尾 / 卫生 / 需拍板

| ID | 状态 | 事项 | 建议 |
|---|---|---|---|
| **T-11** | pending（低） | `occt-topo` 843 条编译警告（含 `let _ =` 类） | 提交前清 unused import/variable（不动行为） |
| **T-29** | pending | `specs/_coverage.md`（2026-09-09）与 `specs/_brepmesh_align_review.md` 已过期 | 与 §2 同步刷新，或标注"仅历史" |
| **T-23** | pending（需定性） | ATU01038 顶点/面 −1.8%（UV-grid vs deflection-adaptive） | 不改；parity 已明确不断言密度 |
| **T-33** | pending（低） | §13 登记：`refine_angle_2d` 两处替换（`Geom2dInt_GInter`→`geom2d_api::intersect_curves` 等）需补注释；"pcurve 参数化 ≠ 3D 边参数化"的退化情形需标未移植 | 纯注释批（不动行为） |
| **T-34** | pending（低） | §13 登记：`is_split_to_reverse_edge` 的 `edge_parameters` ↔ `BRepLib::FindValidRange` 替换需登记 | 同上 |
| **T-25** | **设计已出**（`specs/_design_architecture_t25_t28.md` §T-25） | `GeometryRegistry` 侧表（边/面几何不在 `TShape` 内）——架构级重构；实测 127 文件/501 处引用，`global()` 279 处/88 文件；6 步增量方案（第 2 步为**纯转发层**，101 处调用零改动） | 需专门设计批 ⇒ 设计已交付；实施 5–8 轮 |
| **T-26** | ✅ **done（round 79）**（设计见同上文档 §T-26） | ~~无面 operand 走 `voxel_fallback`~~ —— 该函数**已于第 59 轮 `393ebb8` 摘除**；真正余项 = 端口 `is_empty_shape`（`bop_bop.rs:42-44`「无顶点且无面」）**≠** OCCT `IsEmptyShape`（`BOPTools_AlgoTools3D.cxx:732-854`「无几何」），且 `TreatEmptyShape`（`BOPAlgo_BOP.cxx:214-256`）未移植 | 语义订正；1–2 轮 |
| **T-27** | ✅ **done（round 78）**（设计见同上文档 §T-27；动作 = **摘除**而非接线） | `is_covering_face`（`bop_build_solids_leftover.rs:394-403`）是**自创谓词**，整链只在死模块内闭合（头注自述 *Not called by `BopBuilder`*；`leftover::` 生产 0 命中）；OCCT 的 `GetFaceOff`（`BOPTools_AlgoTools.cxx:994-1095`）→ `IsInternalFace`(:977) → `ClassifyFaces` → `BuilderSolid`(:673) 与 `ShellSplitter`(:359) 两条活链**端口都已有**（`algo_tools_face.rs:433`、`shell_splitter_block.rs:359`） | 1 轮（含判定性探针） |
| **T-28** | **设计已出**（同上文档 §T-28） | ImpPrm HVertex 合并（`IntPatch_ImpPrmIntersection.cxx:221-469` 的 `IsNew==false` 分支，端口 `intpatch_impprm.rs:243 compute_tangency` 只实现 `IsNew`）；`intana`/`intpatch` 重叠实测双份 closed form + 两套 surface×surface 分派器（清单见设计文档）；OCCT 权威 `IntPatch_ImpImpIntersection.cxx` 全部调 `IntAna_QuadQuadGeo` ⇒ 合并方向 = 删 intpatch 版 | 6–10 轮，分 6 步 |

### 3.5 已归档（完成或已决；明细见 §3.6 与 §7）

| 组 | ID | 一句话 / 证据 |
|---|---|---|
| 门禁收敛 | **T-02**、**T-06** | 内腔夹具按契约建模（lib 1293/1 当时）；椭圆参数期望订正（`occt-geom --lib` 153/153） |
| 覆盖/门禁扩展 | **T-21**、**T-22**、**R2-14**、**R2-15** | ATU01038 + `data/occ` 5 模型锁 parity（14/14）；§2 快照按 A/B 实测重写（第 111 轮） |
| A0–A26 审计（已完成项） | **T-35**（A0）、**T-38**、**T-39**、**T-40**、**T-42**、**T-43**（余项=T-66）、**T-45**、**T-46**（3/4，余项=T-73）、**T-47**、**T-50**、**T-53**、**T-56**、**T-57**、**T-58**、**T-60**、**T-61**、**T-63**（余项=T-75）、**T-64**、**T-65**、**T-66**、**T-70**、**T-71**、**T-72**、**T-73**、**T-74**、**T-75**、**T-76**、**T-77**（=R2-1）、**T-81**、**T-84** | 各批见 §3.6「批 1–批 63」；验收一律为"该 crate `--lib` + 四道 STEP 门禁不劣化" |
| 模块自报 PARK | **T-12**（批 87+90）、**T-15**（批 89，`CorrectParameter`）、**T-17**（批 88）、**T-31**（比 curve 句柄）、**T-32**（批 100：`BRepLib_MakeWire` 移植 + `OrientClosedSolid`；输入侧完成，T-01 剩余偏差转 R2-7） | 均在全量门禁下无回归 |
| R2 已完成 | **R2-1**（批 80）、**R2-2**（批 81）、**R2-4**（批 83）、**R2-5**（批 82）、**R2-17**（批 99）、**R2-18**（批 85+86，余项见 3.3）、**R2-18b**（批 86）、**R2-20**（批 95）、**R2-21**（批 93+94）、**R2-22**（批 96） | 详见 §3.6 批次条与 §7 第 101–112 轮 |
| 批次 91–99（R2-6/19 前置 + 写侧） | **R2-6 前半**（批 91 圆锥引擎、批 92 入口、批 93 非修剪圆锥臂、**批 98 trimmed-BSpline 臂**） | 批 98 另含 `Geom_BSplineCurve::{Segment,SetOrigin,InsertKnots}` 与 IGES/STEP 写侧接线 |
| 卫生（已决） | **T-07**、**T-08**、**T-09**、**T-10**、**T-24** | 工作树已干净（仅未跟踪 `data/iges/`）；`data/occ*`、`data/occ-ref/` 已入库；`output/` 按 §11 不入库 |
| 早期任务 | **T-05 的旧"offset 公式同源"假设**（证伪，见 3.1/T-87）、**T-36 的 A0–A26 分期**（余项已全部落成上表 ID） | — |

### 3.6 历史明细（留档：R2 明细表、P0–P5 表、批 1–99 与 T-6x/T-8x 诊断）

> ⚠ **以下全部为历史明细**：其"状态"列写于各自批次，**已过期**；监控与判定一律以 **§3.1–3.5** 为准。保留目的是查证据（审计字段、探针数字、OCCT 行号、批次验收记录）。

状态枚举：`pending` / `in_progress` / `blocked` / `done` / `waived`（豁免必须写依据与引用）。

#### R2 — 新一轮待办清单（2026-09-20 第 93 轮重梳；源 = 全量 T-01…T-85 / A0–A31 状态扫描，取代旧 goal 的轮次预算）（留档）

> **执行口径（用户指令，仍然有效）**：以**能翻译成代码的 OCCT 控制流**为补齐缺口的扩展原则，**以编译为初步验证**；接到"导出 obj 测试"指令前**不自行发起测试/导出**（本节所有"验收"栏都是**待指令后统一跑**的口径，不是当下动作）。禁止特例补丁与为对齐新写单元测试。每批同步本表 + §7 进度日志（批 NN），并监督基线不劣化。
> **本清单已替代**：旧 goal `goal-5b93c8b5`（A0–A26 逐项整改，90 轮用尽、phase=blocked）——A0–A31 中 ✅/◐ 的结论仍有效，未完成部分全部并入下表。

| R2 | 对应卡 | 任务（可翻译的 OCCT 控制流） | 现状 / 证据 | 前置 | 验收（待指令后跑） | 批次 |
|---|---|---|---|---|---|---|
| **R2-1** | **T-77** | `IntTools_EdgeEdge::FindSolutions` 的 bbox 递归：`BndBuildBox`/`FindParameters`/`IsIntersection`/`CheckCoincidence`/`SplitRangeOnSegments`/`MergeSolutions`（`cxx:290-349`、`:353-549`、`:553-671`、`:675-779`、`:780-825`、`:826-901`、`:1060-1146`、`:1150-1206`、`:1210-1362`、`:1366-1406`、`:1410-1452`、`:1456-1659`） | **✅ 批 80（第 93 轮）**：新增 `edge_edge/find_solutions.rs`（~1100 行）逐段移植上述全部成员与文件级 helper；`Prepare` 补齐 `ResolutionCoeff`/`Resolution`/`myPTol1-2` 与 `TypeToInteger` 边交换（`mySwap`）；`perform` 的通用曲线分支改走忠实递归，替代件（"距离极值在容差内即交点"）与其 Newton polish/dedup **整体删除**；`CommonPrt`/`PntOn2Faces` 的映射按 `AddSolution` 的 `mySwap` 语义回写到调用方边序 | 无 | 该 crate `--lib` + 四道 STEP 门禁不劣化（**待指令后跑**） | 批 80 ✅ |
| **R2-2** | **A16 求交 half** | `geom_api::{curve_surface_intersections, curve_curve_intersections}`（采样器）→ 忠实 `IntCurveSurface_Intersection`（+`IntPatch_Intersection`）/ 曲线–曲线 `IntTools_EdgeEdge` | **✅ 批 81（第 93 轮）**：生产调用点全部改走忠实件——`inttools/intersections.rs` 的 4 处曲线–曲线调用（`:278,285,292,295`）→ `edge_edge::EdgeEdge`（= `IntTools_EdgeEdge::Perform`，批 80 的 `FindSolutions`）；非平面边–面（`:450`）→ `intcurvesurface::perform_curve_surface`（= `IntCurveSurface_HInter`，对齐 `IntTools_EdgeFace.cxx:426-445` 的 `Perform` + `W()` 区间过滤）。`geom_api` 的 `curve_surface_intersections`/`curve_curve_intersections`/`curve_curve_distance`/`surface_bbox` 与 `minimize_1d`/`dist_curve_surface`/`curve_window_for_bbox`/`sample_curve` **整体删除**（+3 个自带单测，仓内零消费方） | 无 | 该 crate `--lib` + 四道 STEP 门禁不劣化（**待指令后跑**） | 批 81 ✅ |
| **R2-3** | **T-78 余项** | IGES 写侧 2D（UV）曲线：`BRepToIGES_BRWire::TransferEdge(edge, face, originMap, length, false)`（`cxx:340-588`；平面直接返回空 `:374-377`）、非平面面逐类型 UV 修正并把 `PreferenceMode` 2→3、周期面反周期化与 `periodicU/V`（`GeomToIGES_GeomSurface.cxx:244-343`）、整周椭圆（`GeomConvert_ApproxCurve` `:620-645`） | **◐ 批 84（第 93 轮）完成"周期面"一半**：新增 `occt-core/src/bspl/unperiodize.rs`（`BSplCLib::PrepareUnperiodize`/`Unperiodize` + 面版环绕语义，`BSplCLib.cxx:2967-3080`、`BSplSLib.cxx:2331-2380`），`emit_bspline_surface` 改为：**`periodicU/V` 取原始 `IsUPeriodic/IsVPeriodic`**（`cxx:235-236`/`:448-449`）、周期方向读取前**反周期化**（knots+poles+weights 同一映射）、`closedU/V` 取反周期化后的极点比较（`cxx:347-348`）、bounds 修正的**两条臂**（非周期钳制／周期 `AdjustToPeriod` 平移＋单周期截断，`cxx:244-302`）。**仍 UNPORTED**：`SetUOrigin`/`SetVOrigin` 重定原点（`cxx:310-320`/`:330-340`，`Geom_BSplineSurface_1.cxx:1026-1130`）——仅当区间跨越周期原点时有差异；**余项（新卡）**：R2-19（2D/UV 曲线 + PreferenceMode 3）、R2-20（整周椭圆，前置 R2-6） | 整周椭圆段需 **R2-6**（`GeomConvert`） | `iges_check` 全量 ok + `--lib` + 四道 STEP 门禁 | 批 84 ◐（**批 95：整周椭圆段 = R2-20 ✅**，余 R2-19） |
| **R2-19** | 新（批 84 派生，T-78 余项 a） | **非平面面的 2D（UV）曲线**：`BRepToIGES_BRWire::TransferEdge(edge, face, originMap, length, false)`（`BRepToIGES_BRWire.cxx:340-588`；平面面**直接返回空** `:374-377` ⇒ 端口现状对平面面已忠实）、逐类型 UV 修正（回转面 u/v 反转、柱/锥/拉伸面原点平移）并把 142 的 `PreferenceMode` 由 2 改 3（3D+UV 都写） | 端口 `emit_face_surface` 现在只写 3D 曲线、`PreferenceMode = 2`、`curve_uv = 0`。**起手清单（第 95 轮勘察，`cxx` 行号已核）**：① 端口**没有 2D 曲线 IGES 发射器**（OCCT 侧是 `Geom2dToIGES_Geom2dCurve::Transfer2dCurve`，把 2D 曲线按 Z=0 写成 110/100/104/126 等 3D 实体）⇒ 需新增（可复用现有 3D 发射器 + z=0 提升）；② `ShapeBuild_Edge::TransformPCurve(C2d,Trsf2d,uFact,First,Last)`（`cxx:540-541`、`:555`）未移植；③ 逐类型变换清单：`needShift`（柱/锥，`cxx:403-406`，平面模式关闭、拉伸被注释掉）、回转面（基曲线为直线时 `needShift`，`:408-422`）、周期 B 样条位移（`AdjustToPeriod`，`:435-455`）、柱/锥/球镜像+平移（`:457-466`）、回转/环面镜像（`:468-474`）、`uFact`（拉伸/柱/锥/回转，`:510-538`）、拉伸面反向平移（`:546-556`）、边 REVERSED 时反转 2D 曲线（`:558-565`） | 需 pcurve 提取（`BRepTool::curve_on_surface`/`pcurve_full` 已有忠实件） | `iges_check`（实体数与引用变化属预期）+ `--lib` + 四道 STEP 门禁 | 待做 |
| **R2-20** | 新（批 84 派生，T-78 余项 c） | **整周椭圆**：`GeomToIGES_GeomCurve::TransferCurve(Geom_Ellipse)` 对整周椭圆改走 `GeomConvert_ApproxCurve`（`GeomToIGES_GeomCurve.cxx:620-645`；`approx.HasResult()` 失败时回落 `GeomConvert::CurveToBSplineCurve(copystart, Convert_QuasiAngular)`，`:637-640`） | **✅ 批 95（第 103 轮）**：`iges.rs::emit_whole_period_ellipse` 按 `cxx:620-645` 落地——`pos.Rotated(pos.Axis(), Direct() ? Udeb : 2π-Udeb)`（`:626-628`）→ `CurveToBSplineCurve(copystart, QuasiAngular)`（`:639`，批 93 已具备）→ `BsplCLib::Reparametrize(Udeb, Udeb+2π, Knots)` + `SetKnots`（`:641-643`）→ `TransferCurve(BSpline)` 的周期→非周期转换（`:294-307`，`SetNotPeriodic`，批 95 新增）→ 既有 126 发射器。整周椭圆自此写 126 而非 104。**UNPORTED（就地注明）**：OCCT 先试的 `GeomConvert_ApproxCurve`（`:632-636`）未移植 ⇒ 端口走**回落臂**，几何同为该圆锥的精确有理形式、结点向量与 OCCT 的 approx 结果不同 | R2-6 ✅ + R2-21 ✅ + 本批新增 `Reparametrize`/`SetKnots`/`SetNotPeriodic`/`GpAx2::Rotated`/`Geom_BSplineCurve::distinct_knots_and_mults` | `iges_check`（整周椭圆由 104 变 126 属预期）+ `--lib` + 四道门禁 | 批 95 ✅ |
| **R2-21** | 新（批 91 派生，**发现**） | **`Geom_BSplineCurve` 曲线周期表示缺失**：`Geom_BSplineCurve::SetPeriodic`（`Geom_BSplineCurve.cxx:777-815`）与 `SetNotPeriodic`（`:974-...`）未移植；端口 `GeomBSplineCurve` 的 `periodic` 字段**全仓没有任何一处置 `true`**（`bspline_curve.rs:14,21,41` 三个构造器恒 `false`，无结构体字面量，`git grep '\.periodic = '\|'set_periodic'` = 0 命中；STEP/IGES 读侧一律走 `new`/`rational`，见 `step/read_topology.rs:1183,1186,1235,1266`） | **✅ 两半均已落地**：**库层（批 93，第 101 轮）** = `bspl::knots::{knot_sequence_length,knot_sequence_periodic}`（`BSplCLib::KnotSequence` 周期臂 `BSplCLib.cxx:455-474/488-548`）+ `GeomBSplineCurve::set_periodic`（= `SetPeriodic` `Geom_BSplineCurve.cxx:777-815`，已周期者按 OCCT 等价 no-op）+ 周期求值（`d0/d1/d2` 周期分支走 `curve_dn::dn(u,0/1/2,…)` = `PrepareEval`+`Bohm`+`RationalDerivative`，`BSplCLib_CurveComputation.pxx:777-830`）+ `GeomConvert::CurveToBSplineCurve` 非修剪圆/椭圆臂解除 UNPORTED；**读侧（批 94，第 102 轮）** = `StepToGeom::MakeBSplineCurveCommon`（`StepToGeom.cxx:776-927`）接入 STEP 读曲线臂：重复结点归并（`Epsilon(\|last\|)`）/重数 > degree+1 钳制并裁剪首尾极点/`shouldBePeriodic` 判定（`:873-894`）/闭曲线强制 `SetPeriodic()`（`:920-926`，`IsClosed()` = `StartPoint().SquareDistance(EndPoint()) <= Precision::Computational()`，`Geom_BSplineCurve.cxx:146-149`）。`SetNotPeriodic`（`:974-1014`）仍 UNPORTED（消费方＝trimmed-BSpline 臂，需 `Segment`） | 无 | **⚠️ 读侧为 live 路径，需门禁差分**（见 §6.2 第 94 行） | 批 93 ✅ + 批 94 ✅ |
| **R2-22** | 新（批 94 派生） | **2D（pcurve）B 样条周期表示缺失**：OCCT 的 `StepToGeom::MakeBSplineCurve2d`（`StepToGeom.cxx:952-963`）走的是同一个 `MakeBSplineCurveCommon`（含 `shouldBePeriodic` 与闭曲线强制 `SetPeriodic()`），而端口 `read_topology.rs:1654-1668`、`:1876-1887` 的 2D 臂仍只做 `expand_knots` + `Geom2dBSplineCurve::new`；端口 `Geom2dBSplineCurve`（`occt-geom2d/src/bspline_curve.rs:9-14`）**结构里没有 `periodic` 字段**，且 `continuity`/`parameter_intervals` 把 `false` 写死（`:93-110`） | 影响：周期 pcurve 在端口退化为非周期表示（参数越界即失效），R2-19（IGES 2D/UV 曲线）与 `Geom2dConvert::CurveToBSplineCurve` 都要用到。**批 96 起手勘察（第 104 轮，`cxx` 行号已核）**：① `Geom2d_BSplineCurve::IsClosed()` = `StartPoint().SquareDistance(EndPoint()) <= Precision::Computational()`（`Geom2d_BSplineCurve_1.cxx:147-150`，与 3D 同式），`First/LastParameter` 亦取 `myFlatKnots`（`:349-352`、`:419-422`）⇒ 端口现有 `knots[degree]`/`knots[len-1-degree]` 公式**周期/非周期都对**，无需改；② 2D 曲线**无权重**（`xs/ys` 双数组，非有理）⇒ 周期求值可直接复用 `curve_dn::dn(u, n, 0, degree, true, &poles_3d(), None, &knots, None)`；③ `set_periodic`/`set_not_periodic` 与 `distinct_knots_and_mults` 可直接照搬批 93/95 的 3D 版（底层 `occt_core::bspl::{knots,locate,unperiodize}` 与维数无关）；④ 读侧建议把 `make_bspline_curve_with_knots` 的"归并/钳制/裁剪/周期判定"抽成与极点类型无关的公共段（返回 distinct 结点+重数+极点区间+`should_be_periodic`），3D/2D 两臂各自构造；⑤ **live 路径**：2D 臂接线后 pcurve 可能变周期 ⇒ 需门禁差分 | R2-21（3D 侧机制已备，可照搬） | 该 crate `--lib` + 四道门禁 | **✅ 批 96（第 105 轮）**：① `Geom2dBSplineCurve` 加 `periodic` 字段（全部构造走 `new`，29 处调用点无需改）、`is_periodic`/`period` 覆盖、`continuity`/`parameter_intervals` 改传 `self.periodic`；② 周期求值 `d0/d1/d2` 走 `curve_dn::dn(u,0/1/2,…,true,…)`（2D 无权重，`poles_3d()` 抬到 z=0）；③ `distinct_knots_and_mults`/`set_periodic`/`set_not_periodic` 照搬 3D 版（`Geom2d_BSplineCurve::SetPeriodic/SetNotPeriodic` `Geom2d_BSplineCurve.cxx:948/1087`）；④ 读侧把 `MakeBSplineCurveCommon` 的公共段抽成与极点类型无关的 `bspline_descriptor`，新增 `make_bspline_curve_2d_with_knots` 并接入 2D 臂（含闭曲线强制周期化，`IsClosed()` = `Geom2d_BSplineCurve_1.cxx:147-150`）。**UNPORTED（就地注明）**：2D 有理复合体的**权重被忽略**（端口 2D 曲线本就无权重） |
| **R2-21 勘察（批 92 追加，`cxx` 行号已核）** | 同上 | **要落地必须补的四段控制流**：① `BSplCLib::KnotSequence(..., Periodic=true)`（`BSplCLib.cxx:503-547`：平结串 = `KnotSequenceLength` `:455-474` = `Σmults + 2(deg+1-Mf)`，先把基串写进 `M1+1..`，再前后各绕一个周期），端口 `banded_interp::knot_sequence` 只做了非周期臂；② `BSplCLib::NbPoles` 周期臂 + `FirstUKnotIndex`/`LastUKnotIndex`（`BSplCLib.cxx:392-451`、`:111-137`，已随批 92 移入 `bspl::knots::nb_poles`），`SetPeriodic` 用它们决定"截到 UKnot 区间 + 端部重数钳到 `degree` + 极点 resize"；③ **求值**：`GeomBSplineCurve::{d0,d1,d2}` 现在**完全忽略 `periodic`**（`bspline_curve.rs:139-183` 直调 `eval_curve*`），只有 `eval_dn`/`continuity`/`parameter_intervals` 把标志传下去；`build_knots` 的**平结串臂也忽略 `periodic`**（`build_knots.rs:11-18` 只取 `knots[index-degree..index+degree]` 窗口）⇒ 要按 `BSplCLib_Eval` 的周期臂（`PrepareEval` 用 `build_knots` + 极点在窗口内**环绕**，`curve_dn.rs:30-68 build_eval` 已是这个模式）补 `d0/d1/d2`；④ `Geom_BSplineCurve::{FirstParameter,LastParameter}` 的周期臂（非周期用 `knots[deg]`/`knots[len-1-deg]`，周期应取 `FirstUKnotIndex/LastUKnotIndex` 对应的 distinct 结点）。**读侧同样缺失**：`StepToGeom.cxx:873-925` 按 `(mults(1)==mults(n)) && (Σmults - mults(1) == NbPoles)` 判 `shouldBePeriodic`，且 `IsClosed()` 时**强制** `SetPeriodic()`（`:921-925`）——端口 `step/read_topology.rs:1183,1186,1235,1266` 一律 `new`/`rational`（恒 `periodic=false`） | 落 R2-21 后：非修剪圆/椭圆臂（R2-6）+ R2-20 才能接线；批 90 的周期盒 arm 同时变为可达 | 批 93 |
| **R2-23** | 新（批 98 派生） | **`Geom_BezierCurve::Segment(U1,U2)`**（`Geom_BezierCurve.cxx:388-425`）＝`BSplCLib::BuildCache(0,1,false,deg,KnotSequence(),poles,weights,coeffs,wcoeffs)`（**曲线版**，仓内 `bspl/build_cache.rs` 现为 `BSplSLib::BuildCache` 面版）+ `PLib::Trimming`（已有：`bspl::plib::trimming`）+ `PLib::CoefficientsPoles`（`PLib.cxx` 未移植） | 它是 trimmed-Bezier 三处消费方的**最后前置**：① `GeomConvert::CurveToBSplineCurve`（`GeomConvert.cxx:300-321`）；② STEP 写侧 `GeomToStep_MakeCurve.cxx:77-82`；③ IGES `TransferCurve(Geom_BezierCurve)`（`GeomToIGES_GeomCurve.cxx:441-448`，经 `Geom_TrimmedCurve` + `CurveToBSplineCurve(RationalC1)`）。批 98 已把 ① 的 BSpline 臂与 ② 的 BSpline 分支接线，Bezier 臂仍就地标 UNPORTED | 无 | 五 crate 编译 + 上述三处接线后 `iges_check`/写侧差分 | 待做 |
| **R2-4** | **T-85 步 2** | 写侧可达性过滤＋重编号 | **✅ 批 83（第 93 轮）**：① 求根——`IgesWriter::roots` 记录 `emit_shape` 的每个顶层 DE（Solid/Shell/Face；Compound 递归即多根，= `IGESControl_Writer::AddEntity` 的多次调用）；② 可达集＝**`Interface_InterfaceModel::AddWithRefs` 的先序 DFS**（`Interface_InterfaceModel.cxx:652-692`：加入自身后按序递归共享实体、已加入者跳过），子序＝`IGESData_GeneralModule::FillSharedCase` 的顺序（**DE 字段 7 的变换矩阵在前，然后才是参数指针**）；③ 重编号＝DFS 访问序（1-based）；④ 参数指针改用**占位符**方案（该行设计的备选）：`emit_refs` 的 10 处发射点写 `#k`，`final_entities()` 在重编号后整体替换为最终 DE 号，另同步 DE 字段 7 的 trsf 指针。**订正设计行**：编号序不是"创建序"而是 `AddWithRefs` 的先序 DFS（OCCT 实际行为），已就地注明 | 无 | `iges_check` 全量对比（实体数减少属预期、orphan 统计应归零）＋ `--lib`＋四道 STEP 门禁（**待指令后跑**） | 批 83 ✅ |
| **R2-5** | **STEP 读侧最后一条** | `STEPControl_ActorRead::TransferRelatedSRR`（`cxx:2679-2729`）＋ `TransferEntity(SRR,…)`（`cxx:1337-1433`）＋ `ComputeSRRWT`（`cxx:2495-2546`） | **✅ 批 82（第 93 轮）**：`step/read_geometry.rs` 新增 `SsrComposer`——`rep_shape` = `TransferEntity(ShapeRepresentation)`（`cxx:2048-2092`：自身形状 + 相关表示；单子绑该形状、多子绑 compound；相关结果覆盖单结果）、`transfer_related_srr`（`cxx:2679-2729`：`nbrep` 由 `Rep1==rep` 决定、`theCund` 累积、返回最后一个结果）、`transfer_srr`（`cxx:1337-1433`：取 Rep1/Rep2、`ApplyTransformation` 施加于单结果或 compound）、`compute_srrwt`（`CARTESIAN_TRANSFORMATION_OPERATOR_3D` → `make_transformation3d`；`ITEM_DEFINED_TRANSFORMATION` 两轴 → `compute_axis_transform`；identity ⇒ 无变换）。**保留 UNPORTED（就地注明）**：`PrepareUnits(Rep2/Rep1)` 逐表示单位上下文、`ComputeTransformation` 的轴归属/交换告警、`TransferRelatedSRR` 的 MDADR/构造几何两臂、`TransferEntity(NAUO,…)`（装配级 `CDSR`/`SRRReversed`）；另有 Rust-only 环保护（SRR 环在 OCCT 会无限递归） | 无 | 四道 STEP 门禁 + `--lib`（**待指令后跑**；无 SRR 的语料上行为不变） | 批 82 ✅ |
| **R2-6** | **A8/T-44 触发项** | `GeomConvert`：`CurveToBSplineCurve`/`SurfaceToBSplineSurface`/`ApproxCurve`（写侧 B-spline 现为 6×6 采样重拟，应写真实极点） | 自创实现已于批 66 删除（零消费方），OCCT 侧移植保留为待办。**◐ 批 91（第 99 轮）前半落地**：圆锥→B 样条**引擎**已移植（`occt-core/src/convert/{conic_to_bspline,conic_curves}.rs` = `Convert_ConicToBSplineCurve` `.cxx:32-785` + `Convert_{Circle,Ellipse,Hyperbola,Parabola}ToBSplineCurve`；参数化 `TgtThetaOver2*`/`QuasiAngular`/`RationalC1` 全移植，`Convert_Polynomial` 因缺 2D `BSplCLib::Trimming` 标 UNPORTED）。**余（批 92 起）**：`GeomConvert::CurveToBSplineCurve` 入口（`GeomConvert.cxx:163-430`；line/Bezier/BSpline-`Segment`/Offset-`ApproxCurve` 各臂）、`GeomConvert_CompCurveToBSplineCurve`（`:135-215`）、写侧接线（`step/write_context.rs::bezier_to_bspline` 手写件 → 忠实入口；`GeomToIGES_GeomCurve.cxx:620-645` = R2-20）。**◐ 批 92（第 100 轮）入口落地**：`occt-geom/src/convert_bspl.rs::curve_to_bspline_curve`（`cxx:163-430` 分派 + `BSplineCurveBuilder` `:60-83` + 静态 `Rational` `Geom_BSplineCurve.cxx:97-108`；配套 `gp::Trsf::SetTransformation(From,To)` `gp_Trsf.cxx:172-194` 与 `BSplCLib::NbPoles` `BSplCLib.cxx:392-451`（新 `bspl::knots::nb_poles`，用于 `CheckCurveData` 的极点数校验）），写侧 `bezier_to_bspline` 手写件已删并改走该入口；**已可用臂**＝trimmed-line、trimmed 圆锥（`RationalC1` 且 `U2-U1<6`）、非修剪 Bezier、非修剪 BSpline 拷贝。**仍 UNPORTED**：trimmed Bezier/BSpline（缺 `Segment`）、`RationalC1` 且 `U2-U1>=6`（缺 `CompCurveToBSplineCurve`）、Offset（缺 `ApproxCurve`）。**◐ 批 93（第 101 轮）非修剪圆/椭圆臂已落地**（R2-21 周期表示就位，见该卡），即 `CurveToBSplineCurve` 现覆盖 OCCT 全部分派，只剩上述三个依赖缺失的臂。**◐ 批 98（第 110 轮）trimmed-BSpline 臂已落地**（`GeomConvert.cxx:322-339`：拷贝 → `AdjustPeriodic` → 整周时 `SetNotPeriodic` → `Segment`）⇒ 仍 UNPORTED 的臂＝**trimmed Bezier**（缺 `Geom_BezierCurve::Segment`，新卡 **R2-23**）、`RationalC1` 且 `U2-U1>=6`（缺 `CompCurveToBSplineCurve`）、Offset（缺 `ApproxCurve`） | 被 R2-3 触发 | 同上 | 批 91+92+93+98 ◐ |
| **R2-7** | **T-01 / T-32 / T-30** | 面级切割簇：`BuildSplitFaces` 的 45 个 `avoid=3` 坏面（1 条原边界碎片 + 2 条截面边未链进环）⇒ `rebuild_split_areas` 30 面 0 块；上游 `mesh→BRep` 三角形 wire 不成链（`T-32`）、`d1u`（`T-30`） | T-01 根因链 4 轮探针实测（§3 P0 表）；`occt-topo --lib` **唯一红** | 无（但需门禁验证 ⇒ 等指令） | `occt-topo --lib` 全绿 + 四道门禁<br>**◐ 第 112 轮分诊（决定性实验）**：同一 `groove()` 换喂**忠实 BRep 圆柱**（`BRepPrimCylinder::make_cylinder(1.0,3.0)`）⇒ `before 9.3754 → after 7.8049`、removed **1.5705**（解析 1.6085）、`abs(after-(before-groove_vol)) = **0.0380 < 0.5**` ⇒ **断言在忠实输入下通过**。夹具 `faceted_cyl` = `Solid::wrap(mesh_cylinder(r,h,24))`（mesh→BRep）⇒ **T-01 是夹具症状，根因 = T-32**（不建议改夹具，夹具是有效输入）；**实施顺序（同波次）**：① `triangulation_to_brep` 补 `CorrectEdgeOrientation`（链向）；② 面朝向按 `BRepBuilderAPI_MakeFace` 绕向-法向约定；③ shell/solid 装配一致性；④ 回归 = `groove_cuts_cylinder` + round 8 两个用例 + 四道门禁 | **待做（T-32 一批）** |
| **R2-8** | **T-03** | `build_split_solids_full` 未跨 solid 合并共面片（boss fuse 未并成 1 solid）；邻居＝§5 的「BOP 簇」（`BuildRC` 集合构造/状态过滤 `bop_bop.rs:187-276` ↔ `BOPAlgo_BOP.cxx:583-711`） | `--test bop_builder2_boss` 1/2（HEAD 0/2，已改善） | 与 R2-7 同族，建议同波次 | 该 test 2/2<br>**◐ 第 112 轮分诊（与 T-01 分族）**：实测 `left: 0, right: 1`（Fuse 结果**没有 solid**）。读测例：工具是 `single_disc_cylinder(0.25,0.8,24)`——`TopoBuilder` 手工装配的有效 BRep（24 侧面 + **整块多边形底盘/顶盘**，wire 天然成链），**不是** mesh→BRep；而**同一断言的 tri-fan 变体 `boss_tri_fan_base_merges_one_solid` 通过** ⇒ **引擎侧**缺陷，范围应钉在"**盘面与 box 顶面共面 + 大顶点数多边形环**"的 Fuse 装配。**下一步**：dump 该 Fuse 的 `BuildRC` 计数与 `fill_images_solids`/`build_split_solids_full` 的每 solid 镜像数、以及 box 顶面（应成 2-wire 环）的 `rebuild_split_areas` 块数，对照 `BOPAlgo_BOP.cxx:583-711` 与 `BOPAlgo_Builder_3.cxx::BuildSplitSolids`；若镜像数=1 则瓶颈同样落在**面级切割**（`perform_loops`/WireSplitter 对共面大环的处理） | 待做（独立探测批） |
| **R2-9** | **T-04** | phase19 两条 PaveFiller 管道用例（got 84 / 88，期望含两个原始 box 之外的更多形状）；HEAD 完全一致 | `--test phase19_integration` 3/5 | 无 | 该 test 5/5<br>**◐ 第 112 轮分诊**：断言是 `> 2*56` / `>= 2*56`（**端口自造阈值**）。实测同一 `PaveFiller`：**单盒 `ds.nb_shapes() = 28`**、2 盒相离 **88**、2 盒重叠 **84** ⇒ ① 56/盒 与端口实际 28/盒 差 2 倍，`112` **不可达**（两条用例与实现是否忠实无关地恒红）；② **重叠(84) < 相离(88)** 反直觉 ⇒ 重叠路径**少登记形状**（与 T-80/T-81"截面边未建/未绑定"同族）。**下一步**：按 OCCT `BOPDS_DS::Init`/`Append` 核对"每盒应登记哪些子形状"，判定 28 是缺注册还是正确；再对 84/88 做**按 `ShapeType` 的组成分解**定位少登记类别。**禁止**把 56 直接改成 28（规则 3 禁止改断言对齐） | 待做（OCCT 对照批） |
| **R2-10** | **T-05** | `offset_geometry_is_consistent` 偏差（`2208.0 vs 1612.9`）：按面类型分解 offset 体积；A0/T-35 已修 ⇒ 需重新测量确认是否同源 | `--test step_geometry_parity` 2/3；HEAD 同 | 无 | 该 test 3/3 或显式豁免<br>**◐ 第 112 轮分诊：T-35 假设证伪，改归"网格朝向 + 体积估计器"族**。实测同一形状五种度量：断言左值 `obj_volume` **2208.0**（散度公式，但网格 **6 个三角绕向反向** ⇒ 无效）｜`shape_volume` **1612.9**｜`brep_gprop::volume` **1665.7**｜`brep_gprop_full::{volume_properties, _gk}` **1400.0 / 1400.0**｜网格射线奇偶体积（绕向无关）**≈2556.6**｜解析（Steiner，10³ 盒向外偏移 r=2）**≈2610.5**。判据：① 按位置焊接后网格闭合流形（`712→448 v`、`boundary=0`、`non-manifold=0`）；② 凸形状下逐三角有符号体积有 **6 个为负** ⇒ 朝向缺陷（新卡 **T-87**）；③ 精确积分器 `1400.0 = (1/3)·600·7` = **只算了 6 个平面面**，20 个曲面面贡献 0（无 wire/pcurve ⇒ T-68 下游）。**下一步**：先修 T-87，再核对 `shape_volume` 公式，并复核 T-68 落地后 `gprop_full` 曲面面贡献是否恢复 | 待做（T-87 → 估计器核对） |
| **R2-11** | **T-80 链** | T-81（`pave_blocks/make_blocks.rs:641` `make_pcurves` 后续，按 `BOPTools_AlgoTools2D.cxx:247-400`）→ T-82（`BuildSolid` 只见盒体：柱体 GF 全 `Internal`，`bop_classify_occt.rs`/`bop_cells.rs`/`bop_builder_solid`）→ T-83（`wire_splitter_block::split_block` 环判定） | 三卡均 pending，取证在批 47–48 | T-81 无前置；T-82/T-83 依序 | `occt-topo --lib` + 四道门禁 | 解锁后 |
| **R2-12** | **T-69**（+ T-68 环面/带孔内环） | Torus face 20 的"同面两条 wire 的 pcurve 落在相差一个周期的 u 窗口"（`update_range` 周期钳制后 wire 0 落在范围外）⇒ 失败面集合是 OCCT 的超集 | T-69 第 1/2 轮诊断已定位到该差异；**Delaunay 层已逐行核对为忠实，不要再在该层找偏差** | 无 | T0M 计数 + 四道门禁 | 解锁后 |
| **R2-13** | **T-67 步 3** | `math_FunctionSetRoot`（1452 行）＋ `Extrema_GenExtPS`（1195 行）＋ `Extrema_ExtPExtS`/`Extrema_ExtPRevS` ⇒ 一般曲面的点–面极值走忠实引擎 | 分步 1+2 已完成（`ExtPs` = `Extrema_ExtPS` 分派层）；一般曲面仍是 24×24 网格 + 数值 Jacobian | 无 | 解锁 A15 曲面族兜底移除、A19 相关项、**T-37**（A1 的 39 处调用点迁移） | 解锁后 |
| **R2-14** | **T-21 / T-22** | 覆盖扩展：新模型（ATU01038 等）正式进门禁；`export_data_obj` 目录参数已支持 | T-21 done（12→14 条 parity）；T-22 门禁已锁 5 模型，3 个模型的**网格侧**缺口 pending | 无 | parity/export 计数 | 夹带 |
| **R2-15** | 文档/快照 | §2 快照刷新（旧表记第 90 轮 1287/1、geom 146/146，而批 19 起多处记录 1293/1、151/151 ⇒ **口径需一次实测统一**）；A0–A31 矩阵同步 | **✅ 第 111 轮（门禁波次）完成**：§2 已改为 `9d8e596`(批 84 末) ↔ `2c9a66b`(批 98) 的同命令 A/B 实测表，并就地订正旧表的 1287/1（实为 **1286/2**，多出未登记的 **T-86** `sphere_iges_has_arc_and_solid`）与 geom 146/146（实为 **143**）；`step_to_obj` 实为批 84 的 12/13 → 批 98 的 13/13 | 本轮已跑 | §2 与实测一致 | **done**（2026-09-21） |
| **R2-16** | 小项池 | T-07…T-20（未提交 wave 收尾/低优先）、T-12…T-18（bnd/intf 缺口：`BoxOptimal`、三角化 arm、`CopyNMVertex`、iso 链臂、`project_act` 缺臂、`Intf_Tool`）、T-23…T-29、T-41/T-49/T-54/T-55/T-59（被 T-68/A13 阻塞）、T-64 | 见各任务行 | 多数无 | 随批夹带 | 夹带 |
| **R2-17** | 新（批 80 派生） | **删除 OCCT 不存在的 circle/circle 快路径**：`edge_edge/solvers.rs::compute_circle_circle_full` 与 `perform` 里的 `(Circle, Circle)` 分派。OCCT `IntTools_EdgeEdge::Perform`（`:185-243`）只特判 line/line，圆–圆同样走 `FindSolutions`（忠实件已在 R2-1 落地） | 该分支为端口自创（模块内已标 UNPORTED）；被 `edge_edge/mod.rs` 的两个直调用例与 `circle_circle_two_hits` 间接覆盖 ⇒ **删除会改行为，需一次门禁验证**（这正是不能与 R2-1 同批做的事） | R2-1 已落地 | `occt-topo --lib`（尤其 `edge_edge::tests::{circle_circle_two_hits, separated_circles_empty}`、`tests_full::circle_circle_full_*`）不回归 | **✅ 批 99（第 112 轮）done**：① `perform` 去掉 `(Circle, Circle)` 分派 ⇒ 与 OCCT 一致地"line/line → 快速重合 → `FindSolutions`"；② 删 `solvers.rs::compute_circle_circle_full`（端口自创的共面圆解析解）与其专用 helper `fallback_general`/`point_hit_on_both`；③ 删 `mod.rs` 三条直调用例（−3 条测试），保留 `compute_line_line_full`（`ComputeLineLine` 忠实件）与 `intersect_edges`（line/line 交叉与测试仍在用）。**证据**：`edge_edge::` 15/15 全绿（含 `circle_circle_two_hits` 由忠实 `FindSolutions` 给出 2 个交点、`separated_circles_empty`、`matches_inttools_circle_circle`）；全量门禁与 §2 基线**逐项相同**（lib 1283/2 = 同两条红，总数 −3 即删掉的用例）；`export_data_obj` **16/16 且计数与 R2-17 前逐位相同**（含 ATU01038 17745/22119） | ✅ 批 99 |
| **R2-18** | 新（批 81 派生，A16 的 2D 余项） | **2D 投影/求交采样族**：`occt-geom2d/curve_ops.rs::{curve2d_intersections, curve2d_closest_point, curve2d_length}` | **◐ 批 85（第 94 轮）投影 half（R2-18a）✅ / 批 86（第 95 轮）求交 half（R2-18b）✅**：投影改走 `extrema2d::point_curve_extrema2d`（`Extrema_ExtPC2d`），求交改走 `IntAna2d_AnaIntersection`（详见下方两行）。**仍 UNPORTED**：非解析组合（B 样条/Bezier/offset/其修剪）的求交＝256×256 采样器（忠实路线 `Extrema_ExtCC2d`）；`curve2d_length`＝Simpson 求积（OCCT `GCPnts_AbscissaPoint`）；`extrema2d` 一般曲线的网格+Newton 播种（A7 家族） | 忠实件：`extrema2d`、`occt-core::intana2d` | 该 crate `--lib` + 依赖 2D 的门禁不劣化（见 R2-18-gate） | 批 85+86 ◐ |
| **R2-18b** | 新（批 86 完成） | **2D 求交改走 `IntAna2d_AnaIntersection`** | **✅ 批 86（第 95 轮）**：`Curve2d` trait 补三个类型查询 `gp_elips2d`/`gp_parab2d`/`gp_hypr2d`（`Geom2dAdaptor_Curve` 的 `Ellipse()/Parabola()/Hyperbola()`，`ellipse.rs`/`parabola.rs`/`hyperbola.rs` 实现 + `Geom2dTrimmedCurve` 转发）；`curve_ops::curve2d_intersections` 先走 `analytic_intersections2d`：按首个曲线给「专用操作数」（8 个 `perform_*`），第二曲线取 `IntAna2d_Conic`；`(Circ, Lin)` 按 OCCT 写法交换实参并在结果里换回参数；`IntAna2dIntPoint::{param_on_first,param_on_second,value}` 即 `(u,v,point)`；加**有界曲线语义过滤**（`IntAna2d` 是对**无界**圆锥曲线求解，越出曲线自身区间的根不算交点；周期曲线按一个周期映射；修剪曲线用其 basis 空间的 `first/last`）。非解析组合仍走采样器并就地标 UNPORTED | 同 R2-18 行 | 同 R2-18 行 | 批 86 ✅ |
| **R2-18-gate** | 新（批 86 派生，**门禁影响面**） | **R2-18a/18b + 批 87/88/89/90 直接改 live 路径**：① 2D——`geom2d_api::{intersect_curves, project_point_on_curve}` 被 `pave_de.rs:150,154`、`pave_blocks/split_edges.rs:174`、`wire_splitter_block.rs:680` 调用；② 盒——`geom_bnd_lib_curve3d::box_curve`（批 87 起圆锥解析盒、**批 90 起周期 B 样条 knot-span 盒****（订正·批 91：曲线周期标志全仓无处置 `true` ⇒ 该 arm 当前**不可达**、门禁影响为 0，见 **R2-21**）**）被 `brep_bnd_lib.rs:64`（`BRepBndLib::Add`）、`geom_bnd_lib_surface3d.rs:496,686`、`int_tools_curve_box.rs:46`（= `IntTools_EdgeEdge::BndBuildBox`，批 80 的盒递归用它）与 `pcurve_full/surface_projector.rs:1493` 消费；③ 投影回退——批 88 改了 `ShapeAnalysis_Curve::ProjectAct` 的 `!OK` 分支，经 `step/read_geometry.rs` 进入 STEP 读侧边域投影；④ **pcurve 结点吸附**——批 89 的 `CorrectParameter`（`Proj.cxx:256-281`）被 `wire_fix.rs:3068,3074`（`ShapeFix_Wire::FixNotchedEdges`）经 `TransferParametersProj::transfer_range` 使用，正是 **T-69 诊断里 `check_pcurves_and_shift`/`fix_lacking_all` 那条链** ⇒ 可能影响 T-69 与网格门禁 ⇒ 2D 求交/投影、盒递归剪枝、STEP 边域投影与 healing 结点吸附都会变，可能改变 T-01/T-03/T-04/T-05/T-69 的结果（可能变好也可能位移） | 批 85/86/87/88/89/90 均在 P0/P1 门禁的下游 | 无 | 门禁跑完先做「批 85–90 前/后」差分（`occt-topo --lib`、phase3/4/5/6、`bop_builder2_boss`、`phase19`、四道 STEP 门禁、`export_data_obj`、`iges_check`）**＋ 批 91–94（新增：⑤ 周期曲线族——批 93 周期表示/求值 + 批 94 读侧闭曲线强制周期化 ⇒ 批 90 的周期盒 arm 自此可达；见 §6.2 第 94 行）** | **门禁波次（最高优先）** — **✅ 第 111 轮已跑完差分（步骤 2）**：`9d8e596`(批 84 末) ↔ `2c9a66b`(批 98) 同命令 A/B，**逐项零差异**（唯一变化是 `step_to_obj` 12/13 → 13/13 的改善）⇒ 批 85–98 **无回归**；本轮**未做二分归因**（无失配项）。逐项结果见 §2 / §7 |

**R2 执行顺序（默认）**：R2-1 ✅（批 80）→ R2-2 ✅（批 81）→ R2-5 ✅（批 82）→ R2-4 ✅（批 83）→ R2-3 ◐（批 84：周期面一半；余 R2-19/R2-20）→ R2-18 ✅（批 85 投影 + 批 86 求交；余非解析组合采样器）→ R2-16 夹带：T-12 ✅（前半批 87、后半批 90）、T-17 ✅（批 88）、T-15 的 `CorrectParameter` ✅（批 89）→ **批 91–98（R2-6/R2-20/R2-21/R2-22/R2-19 前置）** → **第 111 轮门禁波次（步骤 1–5，§2 快照刷新 = R2-15 ✅；差分 = 零回归）** → **步骤 6 ✅ = 批 99（R2-17 结项）** → **步骤 7 ✅ = 第 112 轮分诊（R2-7…R2-10 四条卡已更新，见 §7）**。**下一步（按分诊给出的可执行项，等指令）**：① T-32 一批（链向 + 面朝向 + 装配成套）——它是 R2-7/T-01 的根因；② R2-8 的独立探测批（boss 共面盘 Fuse 的镜像数/块数）；③ T-87（OBJ 绕向）+ `shape_volume` 公式核对（R2-10/T-05）；④ R2-9 的 OCCT `BOPDS_DS::Init` 对照；⑤ 之后 R2-11/R2-12/R2-13 → R2-19/R2-23/R2-16（余项）夹带。

#### P0 — 红门禁（收敛或显式豁免）（留档）

| ID | 任务 | 证据（文件:行 / 实测） | 期望 vs 实际 | 下一步 | 验收 | 状态 |
|---|---|---|---|---|---|---|
| T-01 | `brepfeat::tests::groove_cuts_cylinder` 体积不符 | **根因链（探针实测，2026-09-20，全部已复现）**：① 工具正确（16 边形环，体积 3.2574 = 解析值）；② 两输入全平面 ⇒ `bop_builder::boolean` → `bop_builder2::builder_bop_with_fuzzy` → `perform_internal`（filler → fill_images → build_result → build_shape）；③ `build_result`：`gf_solids=2 gf_faces=256`（装配在跑）；④ `build_shape`：`op=Cut dim0=3 dim1=3 open=false` → `build_rc`；⑤ **`build_rc`：`obj_src=1 tool_src=1 it_shapes=1 check_keys=1 it_exp=1 tool_sets=1` ⇒ 每个 solid 只有 1 个镜像 = 没被切开** ⇒ BuildRC 只能整体取舍 ⇒ Cut 几乎不切、Common/Fuse 空；⑥ `fill_images_solids` 有在跑（`fill_in3d in_parts=1`，两个 argument 都有镜像），但 `build_split_solids_occt` 每个 solid 只回 1 个镜像 | 分派已对上 OCCT（`BOPAlgo_BOP::BuildShape` `:871` = `BuildRC` `:900` + 仅 FUSE 且 dim0==3 走 `BuildSolid` `:902-904`；Rust `bop_bop.rs:548/171/277`），故不一致在**面级切割**：柱面 24 个侧片应在 z=1/z=1.8 与 r=1 处被工具面切开却未切 | **round2 追加（探针实测）**：PaveFiller 侧无问题——`interf_ff=422`、`ff_curves=144`、`ff_pb=144`、`ff_pb_with_edge=144`（截面边真建了）；面级绑定也无问题——224 源面里 `faces_sc=112 / faces_on=112 / faces_in=0`。瓶颈在**切割重建**：`BuildSplitFaces` 建了 112 个任务、`empty_le=0`（都收到分割边集），但 `rebuild_split_areas` 块数分布 = **0块:40 / 1块:41 / 2块:30 / 3块:1** ⇒ 只有 31 面真出多块，40 面镜像列表为空 | **round3 追加（探针实测）**：112 个任务面**全是平面**；`perform_shapes_to_avoid` 对所有面都**没触发**（`avoid_sta=0`），`avoid` 由 `perform_loops` 的"未进环边"后处理填充。按 `(in_edges, avoid, avoid_bnd, loops, areas)` 分组：`avoid=0` 的组都是正常切割（如 in_edges=11→loops=3/areas=2 ×15；in_edges=15→loops=5/areas=1 ×16；in_edges=9→loops=2/areas=2 ×7）；**45 个坏面一律 `avoid=3` 且 `avoid_bnd=1`**（= 1 条原边界碎片 + 2 条截面边都没能链进环），其中 30 面最终 `areas=0`（唯一成环的是被切成 hole 的那条 → 无 growth → 按 `builder_face_occt.rs:389-408` 只对**开放面**兜底，闭合面直接丢弃）、15 面 `areas=1` | **round4 追加（探针实测）**：① **推翻**"同点不同 vertex"假设——112 个面全部 `dup_pos_groups=0`（66 面 avoid=0、46 面 avoid=3，都无坐标重复而 key 不同的顶点）；② 解码首个坏面：它是 `mesh_cylinder(1,3,24)` 造出的**柱面侧面三角形**（顶点 `(1,0,0) / (1,0,3) / (0.9659,-0.2588,3)`，四边形的对角剖分），被工具的 z=1 与 z=1.8 两个平面切成 **3 段**；其 9 条唯一边**全部在输入集里**（竖边分 3 段、顶边 1 条、斜边分 3 段、两条水平弦边**各出现 2 份**）；③ 结果：`in_loop=true` 的只有竖边 z0→1、z1→1.8、斜边 z1.8→1.0、z1.0→0 与两条弦边 ⇒ **下面两段闭合成功，最上面那段（竖边 z1.8→3 + 顶边 + 斜边 z3→1.8 + 一份弦边）4 条边齐全却没成环**，3 条全部落进 `avoid`（其中 1 条是原边界边、2 条是边界边的分割像，与 `avoid_bnd=1` 吻合） | **round5 追加（探针 + OCCT 对读）**：① 逐行对读确认 `SplitBlock`（OCCT `_1.cxx:112-354`）与 `Path`（`:358-617`）**与 Rust 移植一致**（含 `iCnt==1` 早退、同边反向 `aTwoPI`、`aNbWaysInside==1` 覆盖、切环/`iPriz` 逻辑）；② 再排除两个假设：`make_2d` 失败（62 次成功调用 `skipped_no_pcurve=0`）、键合并（`GeometryRegistry::shape_key` = 指针身份）；③ 45 个面只走出**部分环**（15 面：7 边→1 环 3 边；22 面：11 边→1 环 5 边；8 面：13 边→1 环 6 边）——行走在那些顶点卡住；④ **命中根因**：`wire_splitter_block.rs:541 refine_angles` 是**空实现**，只保留了 OCCT 的 `iCntBnd != 2` 早退，对 `iCntBnd == 2` **什么都不做**；而 OCCT `RefineAngles(vertex)`（`BOPAlgo_WireSplitter_1.cxx:925-1032`）恰在此情形用 `RefineAngle2D`（`:1033-1125`）重算内部出边角度，并在 `iCntInt == 2` 时按 `Precision::Angular()` 微调（`:998`）。坏面每个弦端点正是"2 条边界边 + 2 份内部弦副本"= `iCntBnd==2 && iCntInt==2`，完全落在 stub 跳过的那一支 | **round8 结果（根因 + 一次被回退的实验）**：① OCCT `BuildSplitFaces` 取边规则解码（`_2.cxx:412-494`）：退化→`Orientation(anOriE)`；INTERNAL/closed→Fwd+Rev 两份；常规→`Orientation(anOriE)`+`IsSplitToReverseWithWarn`；In/Section paves→Fwd+Rev 两份。端口**完全一致**；实测 **544 个边界碎片、朝向不匹配 0** ⇒ 朝向施加是对的（另发现 `is_split_to_reverse_edge:727` 用 edge-TShape 捷径而 OCCT 比 **curve 句柄**（`BOPTools_AlgoTools.cxx:1462`）⇒ 新开 T-31）；② **T-01 真正根因**：坏面的 wire 由 `mesh_to_brep.rs:102 make_wire(&[e_ab,e_bc,e_ca])` 造出，而 `ensure_edge` 的**共享边缓存**把边按"先请求的三角形"的方向存，**后者复用时未修正朝向** ⇒ 实测该面三条边界边全 Forward 且对角线方向属于邻接三角形 ⇒ wire 不是链接链（两条边从 `(1,0,3)` 出发、两条汇入 `(0.966,0.259,0)`）⇒ splitter 块出现 source/sink ⇒ `Path` 死路 ⇒ 面不切 ⇒ groove 切不动。对照 OCCT：`BRepBuilderAPI_MakeWire::Add` 会做 `CorrectEdgeOrientation`；③ **实验（已回退）**：仅在 `mesh_to_brep.rs:102` 加"按请求方向翻转共享边"后，lib 变 **1292/2**——`groove_cuts_cylinder` 前移到 `:98`"必须减材"断言（after≈before），且 `groove_negative_volume_delta` 新红 ⇒ **单独修 wire 朝向不够**：OCCT 的 `BRepBuilderAPI_MakeFace` 会按"wire 绕向 vs 曲面法向"共同决定**面朝向**，必须与 wire 修正同一步做 | **下一轮（round9）**：在 `triangulation_to_brep` 里按 OCCT 一并移植 **`BRepBuilderAPI_MakeWire::Add`（CorrectEdgeOrientation）+ `BRepBuilderAPI_MakeFace` 的面朝向规则**（并检查 shell/solid 装配的朝向），再复跑 T-01 与全量门禁 | `--lib brepfeat::` + 全量 | in_progress（根因已定位到 mesh→BRep 的 wire/face 朝向前提；实验证明需成套修） |
| T-02 | `builder_solid::tests::inner_cavity_is_absorbed_as_a_hole` 内腔未吸收 | `crates/occt-topo/src/builder_solid.rs:723`（生产代码**未改**，只改测试夹具） | 探针实测：`n=2 faces_per_loop=[6,6] growth=[0,1] holes=[]` ⇒ 内壳被判 growth，故 2 个 area | **已修（夹具不具代表性）**：`IsHole` = 无穷点分类（`BOPAlgo_BuilderSolid.cxx:823-831`）+ 分类器读朝向（`brep_class3d.rs:326`）⇒ 只有**朝向翻转**的壳才是 hole。`BRepPrimBox::make_box_corner` 造的是正朝向盒子（=第二个 growth），旧夹具因此断言了非 OCCT 结果。夹具改为按真实 Cut 分裂的腔壳反转内壳朝向后：`growth=[0] holes=[1]` → 1 area / 12 面 | `cargo test --manifest-path crates/occt-topo/Cargo.toml --lib builder_solid::` → 全绿 | **done**（2026-09-20） |
| T-03 | `boss_single_disc_base_merges_one_solid` boss fuse 未合并成 1 solid | `crates/occt-topo/tests/bop_builder2_boss.rs:114` | `left=0 right=1`（HEAD：两条用例全红） | 已知阻塞点：`build_split_solids_full` 未跨 solid 合并共面片（remember 记录 B3 实体装配缺口） | `--test bop_builder2_boss` | pending |
| T-04 | phase19 两条 PaveFiller 管道用例 | `tests/phase19_integration.rs:40`（got 84）、`:53`（got 88） | 期望"至少包含两个原始 box 之外的更多形状"；HEAD **完全一致** | 本轮未触碰；属 PaveFiller 阶段遗留，单独开波次 | `--test phase19_integration` | pending |
| T-05 | `offset_geometry_is_consistent` 偏差 | `tests/step_geometry_parity.rs:129` | `Offset divergence 2208.0 vs shape_volume 1612.9`；HEAD **同** | 单独定位（offset 体积 vs 面类型分解）；与 T-04 同属遗留。**2026-09-20 全仓审查提示：很可能与 A0（offset 公式写错，见 T-35）同源**——先修 T-35 再看此偏差是否变化 | `--test step_geometry_parity` | pending |
| T-06 | `extrema_pc::tests::point_ellipse_y_axis`（**唯一新增红**） | `crates/occt-geom/src/extrema_pc/tests.rs:57-64` | 期望 `u=3π/2`，实得 `π/2`；同用例的 `distance=3.0`、`p2=(0,1,0)` 断言**均通过** | **已修**：期望改 `π/2` 并注明 `ElCLib::EllipseValue`（`cxx:176-189`，`+Minor*sin(U)`）；**未**回退 `clib.rs` 符号（那是 ATU01038 面 130/140/156/216/222 的修复） | `cargo test --manifest-path crates/occt-geom/Cargo.toml --lib` → **153/153 ✅** | **done**（2026-09-20） |
| **T-87** | 新（第 112 轮步骤 7 分诊发现） | **OBJ/面绕向不一致**：`Offset.step` 的 `brep_to_obj` 网格（712 v / 892 f）里**有 6 个三角绕向反向**——以形状内部点 `(5,5,5)` 为顶点做逐三角有符号体积，6 个得 **‑116.67**、其余 **+116.67**（该形状是凸的，绕向一致时应全为正）。旁证：`brep_to_obj` 逐面写顶点不去重（`v=712`，按位置焊接后 448），焊接后网格**闭合且流形**（`boundary=0 / non-manifold=0`）⇒ 不是缺面，是**朝向** | 后果：任何基于**绕向**的体积/属性（测试里的 `obj_volume` 散度公式 = 2208.0，而网格自身奇偶体积 ≈2556.6、解析 ≈2610.5）都失真；与 round 9 的度量"面法向·wire 绕向为负 112/128"同源 | 定位 `brep_to_obj`（及其上游面/壳朝向：`BRepBuilderAPI_MakeFace` 的绕向-法向约定、`triangulation_to_brep`）为何 6 个三角反向；对照 OCCT `RWObj_CafWriter` + `RWMesh_FaceIterator` 的三角形写出（`Poly_Triangulation` 的 `HasNormals`/朝向语义） | 修后 `obj_volume` 与网格奇偶体积一致（≈2610），`assert_closed_solid` 仍绿 | pending |
| **T-86** | **`iges::tests::sphere_iges_has_arc_and_solid`（本轮 A/B 发现：板内未登记的红，且 §2 旧表误记为"仅 T-01"）** | `crates/occt-topo/src/iges.rs:2221`（`assert!(iges.contains("100") \|\| iges.contains("128"), "expected a circular arc or NURBS curve")`）。**A/B 实测**：`9d8e596`（批 84 末）与 `2c9a66b`（批 98）**都失败**，且单独跑（`cargo test --lib iges::tests::sphere_iges_has_arc_and_solid`）同样失败 ⇒ 非测试顺序/并行导致，也非批 85–98 引入（`iges.rs` 在批 84→95 之间无提交改动）。临时探针（已删）dump 出 `BRepPrimSphere::make_sphere(2.0)` 的 IGES **DE 类型 = [144, 196, 116, 123, 123]**（5 实体），**无 100/110/128/142**（`144` 的裁剪面上没有任何边界曲线） | 期望：球面（`BRepPrim_Sphere` 的 lateral wire = **两条退化极点边 + 两条 u=0/2π 经线**，见 §3 T-68 执行方案）在 IGES 里应有经线的 **100**（`IGESGeom_ToolCircularArc`）或 128 曲线；实际：端口球面**没有边界 wire**（T-68 探针：`sphere faces=1 wires/face=[0]`）⇒ 144 无 outer bound、无曲线 | 归因 = **T-68 的下游症状**（球面 wire/pcurve 落地后此断言自然满足）。不要改这条断言来"修"它（属规则 3 禁止的"为对齐改期望"） | 与 T-68 同批：球面 `wires/face>0` 且 `--lib` 该用例转绿 | pending（被 T-68 阻塞） |

| **T-35** | **A0：offset 曲线公式写错（全仓审查最高危）** | `crates/occt-geom/src/offset.rs:14-19`（3D 用 `p + Offset·Dir` **平移**）↔ `Geom_OffsetCurveUtils.pxx:53-61`（应沿**法向** `p + Offset·(D1×Dir)/‖D1×Dir‖`）；`crates/occt-geom2d/src/offset.rs:19-41` 法向 `(-dy, dx)` ↔ `Geom2d_OffsetCurveUtils.pxx:50` `(dy, -dx)`；两侧 `d1/d2` 均缺 `DNdir` 项（**已双向对读核实**） | 产出错误几何（**公式错误，非近似**）；3D 版被 STEP 读入 `step/read_topology.rs:1049` 直接使用，2D 版被 `geom_bnd_lib_offset2d.rs:101` 消费 | 按两个 `OffsetCurveUtils.pxx` 重写 D0/D1/D2（含 `DNdir` 旋转项与失败返回），并复核上述两处调用点 | 现有 offset 相关门禁不回归（`--lib` + `--test step_geometry_parity`） | pending |
| **T-36** | **全仓忠实度审查 A1–A26 分批整改**（A0 已单列 T-35） | `specs/_audit/_index.md` §3/§7（5 区报告齐；合计自创 **71** / 未登记 **17** / 已登记 **15**，系统性条目 A0–A26） | 按 §7 分批：A6→A7 → A1 → A2+A11 → A13/A23 → A3/A24 → A20/A21/A17 → A18/A19/A22/A25 → A9/A10/A14/A26 → A12 | 逐批对着 `.cxx` 改控制流，**不加 OCCT 之外的规则/阈值**；无对应分支的标 `UNPORTED` + OCCT 文件行号 | 每批跑该 crate `--lib` + 相关门禁（parity 14/14 等）不回归 | pending |

#### P1 — 未提交 wave 的收尾（留档）

| ID | 任务 | 证据 | 动作 | 状态 |
|---|---|---|---|---|
| T-07 | 178 文件未提交（148 M / 30 D / 32 ??；排除参考 OBJ 后 177 文件 +10,053/−2,540） | `git status --short` | 分批提交：① 源码修复（clib/primitives/shape_ops/pcurve_full/shhealing/bnd 等）② 新增单元（13,900 行 + `plib_jacobi_coeffs.pxx` 5,585 行）③ 数据与参考 ④ 清理 | pending |
| T-08 | 30 个 `specs/_phase*_plan.md` / `_loop*_plan.md` / `_rules.md` / `_brepmesh_migration.md` / `_brepmesh_align_tasks.md` 与 `examples/demo.rs` 被删（`Cargo.toml` 同步去掉 `[[example]] demo`） | `git status --short specs examples` | 确认是有意清理；`_brepmesh_align_tasks.md` 里"残留差距"清单要先合并进 `_board.md`（§3 P2 已吸收）再删，避免丢任务 | pending |
| T-09 | `output/` 未忽略（27 文件 **16.35 MB**）且 `data/occ/`、`data/ATU01038.step` 未跟踪 | `git status --short` | 决策：入库存档 or 加 `.gitignore`（当前 `data/output/` 已忽略，但根 `output/` 没有） | pending（需用户拍板） |
| T-10 | `data/occ-ATU01038.obj`（**已跟踪**）被原地重导覆盖 | 相对 HEAD 绕 X 轴 90°；顶点 17,720 → 18,102 | 提交时说明"参考按新输入 STEP 重导"，否则后人无法解释 diff | pending |
| T-11 | `occt-topo` lib 846 条编译警告 | `cargo check` | 提交前至少清 unused import/variable（不动行为），否则真信号被淹没 | pending (低优先) |

#### P2 — 新模块自报的 PARK / UNPORTED（按 OCCT 控制流补齐）（留档）

| ID | 位置 | 未移植内容 | 状态 |
|---|---|---|---|
| T-12 | `crates/occt-topo/src/geom_bnd_lib_curve3d.rs:7` | Ellipse / Hyperbola / Parabola 解析盒（`occt_geom::Curve` 不暴露 `gp_Elips/Hypr/Parab`）→ 现走采样；周期 B 样条 arm（缺 `Segment`/`AdjustPeriodic`，`BSplineCurve.cxx:44-49,299-330`） | **✅ 批 87（前半）+ 批 90（后半）全部 done**：三个圆锥解析盒已移植（`GeomBndLib_Ellipse.cxx:23-114`／`_Hyperbola.cxx:25-140`／`_Parabola.cxx:23-94`）并接入 `box_curve` 的 OCCT 分派顺序；周期 B 样条 arm 也已落地——`ElCLib::AdjustPeriodic` 把区间收进一个周期，分段曲线的结点＝原 distinct 结点按整周期旋转 ⇒ **无需移植 `Segment`** 即可用同一 knot-span `FillBox` 循环复现（`BSplineCurve.cxx:43-50` + `:85-105`），已就地写明等价性论证 | **done（批 87+90）** |
| T-13 | `geom_bnd_lib_surface3d.rs:8` | `BoxOptimal`（`OptimizationHelpers.pxx` PSO+Powell）+ `SurfaceOfExtrusion::BoxOptimal`（`cxx:224-283`）；`catch(Standard_Failure)` 回退 | pending（`AddOptimal` 路径，当前不在 `BRepBndLib::Add` 上） |
| T-14 | `brep_bnd_lib.rs:13` | 三角化 arm（`BRepBndLib.cxx:95-101,157-179`）、仅 pcurve 边的 `BRepAdaptor_Curve::Initialize` 回退（`cxx:93-106`） | pending（缺 `Poly_Triangulation` 存储） |
| T-15 | `shhealing/transfer_params.rs:21` | `CopyNMVertex`（需 `BRep_PointRepresentation`）、`CorrectParameter` knot snap（`Proj.cxx:268-279`）、`myLocation` 统一 identity | **◐ 批 89（第 97 轮）`CorrectParameter` done**：新增 `Curve2d::bspline_knots2d`（`Geom2dBSplineCurve` 实现，扁平结点串），`correct_parameter` 现按 `Proj.cxx:256-281` 递归解包 trimmed/offset 后，对 `Geom2d_BSplineCurve` 在 `PConfusion` 内吸附到**首个** distinct 结点（端口 2D BSpline 非周期 ⇒ `First/LastUKnotIndex` 覆盖全部 distinct 结点）；文件头 UNPORTED 项改为已移植。**仍 UNPORTED**：`CopyNMVertex`（端口无 `BRep_PointRepresentation`，见文件头理由）与 `myLocation`（healing 路径恒 identity，文件头已证） | **◐（CorrectParameter done，其余为已论证的等价省略）** |
| T-16 | `pcurve_full/projection_cache.rs:27`、`pcurve_full/singularities.rs:11` | iso 链在 iso 边界边上的臂选择；`ShapeAnalysis_Surface::myGap` 用残差恢复；奇点数组未缓存（结果等价）；`ProjectDegenerated` 单点重载无调用者 | pending（已写明等价性理由） |
| T-17 | `shhealing/shape_analysis_curve.rs:24` | `project_act` 的 `!OK` 分支缺 Hyperbola/Parabola（`cxx:376-386`）与 Ellipse（`cxx:394-399`）快捷路径 → 走 default 段搜索 | **✅ 批 88（第 96 轮）done**：`!OK` switch 五臂全齐——Circle（`cxx:355-374`）、Hyperbola（`:376-380`）、Parabola（`:382-386`）、Line（`:388-392`）、Ellipse（`:394-399`，含 `anIsClosedCurve`/`aCurvePeriod=2π`），分别经 `clib::{hyperbola,parabola,ellipse}_parameter` + `clib::{hyperbola,parabola,ellipse,line}_value`（Line 的参数用既有 `line_parameter` = `ElCLib::LineParameter`）。原阻塞（`Curve` 缺三个 conic 查询、`clib` 缺 `*_parameter`）已由 A20/A29 与 T-56 解除 | **done（批 88）** |
| T-18 | `crates/occt-core/src/intf/mod.rs:6` | `Intf_Tool`、`Intf_InterferencePolygonPolyhedron.gxx`（3D 多面体干涉，无 2D 调用者） | pending |
| T-30 | `uv_tolerance_2d`（`wire_splitter_block.rs:507`）用 `t3/|d1u|` 近似，OCCT `BOPAlgo_WireSplitter::Tolerance2D`（`_1.cxx:859-881`）用 `GeomAdaptor_Surface::U/VResolution`（平面恒 1.0，B样条面 ×1.1） | 平面面下 OCCT `Tolerance2D≈1.0` vs Rust `≈1e-7` ⇒ `Angle2D` 的 `dt` 差数量级（直线 pcurve 无影响，曲面面会影响取角） | pending（round7 发现，非 T-01 阻塞） |
| T-31 | `is_split_to_reverse_edge`（`algo_tools_face.rs:727`）的捷径比较 **edge TShape 身份**；OCCT `BOPTools_AlgoTools::IsSplitToReverse`（`BOPTools_AlgoTools.cxx:1461-1465`）比较 **curve 句柄**（`aCSp == aCOr`） | **已修（round9）**：捷径改为 `Arc::ptr_eq(edge_curve(sp), edge_curve(or))`，与 OCCT 一致（碎片与父边共享 curve ⇒ 只比朝向，不再落几何切线分支）。实测：全量 `--no-fail-fast` **与基线完全一致（无回归）**；T-01 数值仅末位变化（8.903438339491375→…373）⇒ 本例两种判定同解 | `--lib` + 全量 | **done**（2026-09-20） |
| T-32 | `mesh_to_brep.rs:102` 的 wire 未做 `BRepBuilderAPI_MakeWire::Add` 的 `CorrectEdgeOrientation`（共享边由邻接三角形先建、按彼方向缓存） | 直接后果：mesh→BRep 的三角形 wire **不成链**（实测该面三边全 Forward、对角线方向属邻接三角形）⇒ wire splitter 块现 source/sink、`Path` 死路。**单独修 wire 朝向会回归 2 个用例**（round8 实验：lib 1292/2）；round9 度量补充：mesh→BRep 柱面 **47/96 个面的"面法向·wire 绕向"为负**（绕向约定本身混乱），而旋转体工具 **112/128 一致为负** ⇒ 修复必须与 `BRepBuilderAPI_MakeFace` 的面朝向/绕向约定**成套**做 | **P0 关联**（T-01 根因；round10+ 实施） |
| T-19 | `meshing/triangulator.rs:19` | `addTriange34` / `checkCondition` 快路径（仅 `Perform` 调用，未移植）；注：`BRepMesh_Triangulator` 只被 VRML 读路径引用，**不属** STEP/OBJ 网格管线 | pending（低优先） |
| T-20 | `bop_build_solids_leftover.rs:8` | `merge_sharing_faces` 已删（OCCT `BOPAlgo_Builder_3.cxx:579-616` 无此 stage）；leftover 模块去留待定 | pending |

#### P3 — 对齐覆盖扩展（新模型）（留档）

现状（2026-09-19 实测）：

| 模型 | STEP 源 | OCCT 参考 | Rust 产物 | bbox 一致 |
|---|---|---|---|---|
| ATU01038 | `data/ATU01038.step`（未跟踪，19,637 行） | `data/occ-ATU01038.obj`（1,8102 v）· `output/occ/ATU01038.obj`（18,051 v） | `output/ATU01038.obj`（**第 111 轮实测 17,745 v / 22,119 f**；09-19 记录为 17,767 v / 22,160 f） | ✅ Δ≤9e-6 |
| a3n00 / acs10 / bottom / motoc / top / T0M / TDB | `data/occ/*.step`（8 个） | `output/occ/*.obj`（8 个，9/15 20:56–22:45） | ❌ 无 | — |
| Cube/Cone/Cylinder/Sphere/Torus/HoledPlate/Offset*/Shape*/rev/linkrods/screw/Extrusion | `data/*.step` | `data/occ-*.obj` / `occ-*.obj` | `output/*.obj`（9/17 16:40） | 12/12 门禁绿 |

| ID | 任务 | 说明 | 状态 |
|---|---|---|---|
| T-21 | ATU01038 纳入 parity 门禁 | `tests/step_obj_parity.rs` 原是 **12 条硬编码用例**，无 ATU01038。**已加**一条 `check_parity("ATU01038","occ-ATU01038.obj",1e-4)`（实测 Δ≤1.1e-5；密度不断言）；加用例 = 引用现有基线，符合规则 3 | `cargo test --manifest-path crates/occt-topo/Cargo.toml --test step_obj_parity` → **13/13 ✅** | **done**（2026-09-20） |
| T-22 | 7+1 个 `data/occ` 新模型跑 Rust 导出并对拍 | **已跑 + 已锁门禁（round10/11）**：`export_data_obj` 增加可选目录参数（`-- data/occ`，无参数时行为不变），8 个模型全部导出成功。**bbox 对拍（Rust 网格 vs `output/occ/*.obj`）**：`ATU01038` **0/0/0**、`bottom` **0/0/0**、`motoc` **0/0/0**、`top` **0/0/0**、`T0M` **0/0/0**（密度 −8.8%~+0.8%）→ 这 4 个（+ATU01038 已有用例）**已锁进 `step_obj_parity`**（新用例 `occ_test_model_bboxes_match_occt`，参考入库 `data/occ-ref/`，门禁 **14/14**）；**3 个仍有真实几何差**：`a3n00` Δy=**60.04**、`acs10` Δy=**7.88**、`TDB` Δx=**0.59** | **round12 追加（归因到网格侧）**：`a3n00` 的网格 y-max = **137.539** vs OCCT **77.5**（+60.039），其余轴全对（x 相等、y-min 相等、z 差 0.006）；而 **226 个面的 wire 顶点全部在 |y|≤77.6 内（0 个越界）** ⇒ **解析边界与 OCCT 一致，溢出由网格生成产生**（面被按超出 trim 的范围采样，即 `brepfeat/features.rs` 已注明的"按 UV 窗口网格化"那一类）。逐面单独网格化给出候选越界面（41–49、77–87 等 `Plane`/`Other`，自身网格超出 wire y 约 10；**注意**：单面网格化可能与实体导出走不同路径，此列表为指示性）。下一步配方：在导出路径内部 dump **每面网格 bbox** 与该面 wire bbox 对比，定位后用 OCCT `BRepMesh_FaceDiscret`/`BRepMesh_ModelBuilder` 追该面的 trim/UV 窗口逻辑 | 导出 + bbox 对拍 + parity | **done（导出扩展 + 5 模型锁门禁）**；3 个模型的**网格侧**缺口 pending |
| T-23 | 密度差是否要收 | ATU01038 顶点 −1.85%（vs `data/occ-ATU01038.obj`）/ −1.57%（vs `output/occ`），面 −1.89% / −1.21%。parity 注释明确"不断言密度，UV-grid vs deflection-adaptive 是已知差" ⇒ 若收，需先证明 OCCT 侧公式逐项对上 | pending（待定性） |
| T-24 | `output/` 与 `data/occ/` 的归档策略 | 见 T-09 | pending |

#### P4 — 继承缺口（`specs/_coverage.md` 列出，本轮未处理）（留档）

| ID | 任务 | 状态 |
|---|---|---|
| T-25 | `GeometryRegistry` 侧表（边/面几何不在 `TShape` 内）——架构级重构 | pending |
| T-26 | 无面 operand 仍走 `voxel_fallback`（绿路径 `bop_builder2`） | pending |
| T-27 | `GetFaceOff` / 角法向：leftover `is_covering_face` 死路径未接 `BopBuilder` | pending |
| T-28 | ImpPrm HVertex 合并（`IntPatch_ImpPrmIntersection` cxx 329–465）未移植；`intana` 与 `intpatch` 重叠未合并 | pending |
| T-29 | `specs/_coverage.md` 已过期（2026-09-09；行数/测试数与本轮不符）→ 本轮收敛后刷新 | pending |

#### P5 — 自创 → 忠实移植（审查 A0–A26 逐项，2026-09-20 立项）（留档）

> 来源：`specs/_audit/_index.md`（5 区报告 + A0–A26）。**父任务 T-36**，下列每项一个可跟踪 ID。
> 验收统一为：① `cargo check` 过；② 该 crate `--lib` 不低于基线（**以 §2 实测为准**：topo 1283/2、geom 143、geom2d 72、core 290、math 215）；③ 相关门禁不回归（`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3）；④ 无对应 OCCT 分支时标 `UNPORTED` + OCCT 文件行号，**不加新规则/阈值/启发式**。

| ID | A# | 位置（文件:行） | 自创内容 → OCCT 对应 | 阶段 | 状态 |
|---|---|---|---|---|---|
| T-35 | **A0** | `occt-geom/src/offset.rs:14-19`、`occt-geom2d/src/offset.rs:19-41` | 3D 沿参考方向平移 / 2D 法向反号、缺 `DNdir` → `Geom_OffsetCurveUtils.pxx:47-115`、`Geom2d_OffsetCurveUtils.pxx:43-` | **1（先修）** | **done**（2026-09-20） |
| T-63 | A0 派生 | `occt-geom/src/curve.rs:8-11`（`Curve::d3`）、`occt-geom2d/src/curve.rs`（`Curve2d::d3`）、各具体曲线 | `EvalD3`（`Geom_Curve::EvalD3` / `Geom2d_Curve::EvalD3`）无任何具体曲线实现 ⇒ offset 的 `CalculateD2` 里 `D2Ndir` 项恒为 0；`EvalD2` 的 `AdjustDerivative` 奇异支路（`Geom_OffsetCurve.cxx:311-330`、`Geom2d_OffsetCurve.cxx:265-280`）未移植（`isDirectionChange` 恒 false） | 1 后 | **◐ 步 1 done**（2026-09-20）：① 新增忠实件 `clib::{circle,ellipse,hyperbola}_d3`（3D，`ElCLib.cxx:435/464/494`）与 `clib::{circle,ellipse,hyperbola}2d_d3`（2D，`cxx:809/843/878`，其中双曲 `V3=V1`、圆/椭圆 `V3=-V1`）；② 3D 具体曲线补 `d3`：`GeomCircle`/`GeomEllipse`/`GeomHyperbola` 直接调 clib，`GeomBSplineCurve` 走 `eval_dn(u,3)`（`BSplCLib::DN`），`GeomLine`/`GeomParabola` 依 OCCT 保持零（`Geom_Parabola::EvalD3` `cxx:195-200` 即 `V3=0`）；③ 2D 补 `d3`：`Geom2dCircle`/`Geom2dEllipse`/`Geom2dHyperbola`；④ **顺带修掉两处假 `d2`**：3D `GeomHyperbola::d2`/`GeomParabola::d2` 原本返回**零二阶导**，现按 `ElCLib::HyperbolaD2`（`cxx:378-401`）/`ParabolaD2`（`cxx:404-431`，含 `\|Focal\|<=gp::Resolution()` 退化支路）实现，并新增 `clib::{hyperbola_d2,parabola_d2}`。**未完成（仍在 T-63）**：`Geom_OffsetCurve::EvalD2`（`cxx:311-330`）/`Geom2d_OffsetCurve::EvalD2`（`cxx:265-280`）的 `AdjustDerivative` 奇异支路——`isDirectionChange` 恒为 `false`。**⇒ 步 2 done（2026-09-20 第 41 轮）**：3D/2D `AdjustDerivative`（`Geom_OffsetCurveUtils.pxx:322-390` / `Geom2d_OffsetCurveUtils.pxx:296-364`）逐行落地并接进 `EvalD2` 奇异支路（3D `cxx:311-330`、2D `cxx:258-280`，`isDirectionChange` 不再恒 `false`）；为 `AdjustDerivative` 需要的 `EvalDN` 4/5 阶补齐 3D `ElCLib::{Line,Circle,Ellipse,Hyperbola,Parabola}DN`（`ElCLib.cxx:911-1045`）与对应 5 个初等曲线 `eval_dn` 覆写、2D 新增 `Curve2d::eval_dn` 并覆写 5 个初等曲线、两侧 `GeomTrimmedCurve` 转发基曲线（`Geom_TrimmedCurve.cxx:231-243`、`Geom2d_TrimmedCurve.cxx:273-283`）。**余项转 T-75**（offset 自身的 `EvalD3`/`EvalDN`） |
| T-64 | **A27**（新） | `occt-geom/src/{line,circle,ellipse,hyperbola,parabola,plane,cylinder,cone,sphere,torus}.rs`、`occt-geom2d/src/{line,circle,ellipse,hyperbola,parabola}.rs` 的 `continuity()` | 解析曲线/曲面一律返回 `3`（G2），OCCT `Geom_Conic::Continuity()`=`GeomAbs_CN`（`Geom_Conic.cxx:32-35`）、`Geom_Line`(`:128-131`)/`Geom_ElementarySurface`(`:25`)/`Geom2d_Line`(`:172`)/`Geom2d_Conic`(`:48`) 同为 CN(6)；端口 B 样条侧 `bspl::local_continuity` 单跨返回 6，证明编码约定一致 ⇒ 解析类少报连续性，影响 `NbIntervals`/`parameter_intervals` 消费方 | **已尝试 → 回退**（2026-09-20）：仅把这 15 处改成 6 后，`step_to_obj` 由 13/13 变 12/13——**Sphere 网格变空**（`tests/step_to_obj.rs:240` round-trip 顶点数 0）。⇒ 不能单独改：某个消费方把"报出的连续性"当跨度上限（嫌疑点 `meshing/range_splitter/param_set.rs:295` 的 `continuity <= curve.continuity()`、`meshing/edge_discret.rs:1259` 的 `pcurve.continuity().min(surface.continuity())`）。下一步：先按 OCCT 对齐该消费方（`GeomAdaptor_Surface::NbUIntervals`/`NbVIntervals` 语义），再改这 15 处 | pending（被消费方阻塞） | **done（第 69 轮批 54）**：15 处 `continuity()` 改 6（出处见批 54）；原阻塞（Sphere 网格变空）在 T-68 球面 pcurve 落地后已不复现，`step_to_obj` 13/13。
| T-37 | A1 | `occt-topo/src/brep_surface.rs:234-272`（39 处/27 文件） | 网格扫描 + 6 轮二分 → **需先移植** `Extrema_ExtPS`/`Extrema_GenExtPS`/`Extrema_ExtPElS`（仓内 `extrema_surf` 的**分派层**已忠实（T-67 步 2），但一般曲面仍走 24×24 网格 + 数值 Jacobian 的替代件） | 3 | pending（被 T-67 步 3 / `GenExtPS` 阻塞） |
| T-67 | A1 前置 | 新增移植：`occt-geom/src/extrema_surf/`（`point_surface_extrema.rs`） | 移植 `Extrema_ExtPS.cxx`(428)+`.hxx`(140)、`Extrema_GenExtPS.cxx`(1195)+`.hxx`(181)、`Extrema_ExtPElS.cxx`(519)+`.hxx`(76)：类型分派 + iso 退化处理（`IsoIsDeg`）+ 窗口/`TreatSolution` + 逐 C2 区间采样（`GetGridPoints`/`BuildGrid`）+ `math_FunctionSetRoot` 解析系统（`Extrema_FuncPSNorm`） | **3 前置** | **◐ 分步 1+2 done**（2026-09-20 第 49 轮）：解析臂接线（批 3）+ **`Extrema_ExtPS` 分派层 `ExtPs`（`point_surface_extrema.rs`）**：三引擎 switch、±1e10 钳制、`nbU/nbV` 44/32/300、`IsoIsDeg`、`TreatSolution` 周期归一与窗口测试、`IsDone`/`NbExt`/`TrimmedSquareDistances`；`clib::in_period` 补齐。**分步 3 pending**：`Extrema_GenExtPS` 引擎 + `Extrema_ExtPExtS`/`Extrema_ExtPRevS`（两者均就地标 UNPORTED） |
| T-38 | A2 | `occt-topo/src/brep_extrema.rs:252-301` | 7×7 采样 + 射线奇偶 → `SolidClassifier`/`algo_tools::compute_state`（`BOPAlgo_BuilderSolid.cxx:835-860`） | 4 | **done**（2026-09-20 第 44 轮）：`is_inside` 对 Solid/CompSolid 改走 `brep_class3d::SolidClassifier::classify`（`BRepClass3d_SClassifier::Perform`），非实体形状仍走原网格射线奇偶并**就地标 UNPORTED**（OCCT 用 `BRepClass3d_SolidExplorer` 处理 shell，端口未接线）。**同批修掉一处真实缺陷**：`SolidExplorer` 缺 `myMapEV`，ON 判定误用"整形状的全部顶点/边"，导致**实体内部的孤立顶点**会把查询点判成 `On`；现按 `BRepClass3d_SolidExplorer::Init`（`cxx:930-982`）只收"非 INTERNAL/EXTERNAL 面的、非 INTERNAL/EXTERNAL 且非退化的边及其顶点"。证据：`phase3_integration` **3/4 → 4/4**（原红 `primitives_measure_correctly` 的球内点分类转绿），`occt-topo --lib` 1293/1、四道 STEP 门禁 14/14、13/13、11/11、2/3、`export_data_obj` 16/16 且计数与基线逐位一致 |
| T-39 | A3 | `occt-topo/src/step/format.rs:104-144,663-680`、`step/write_context.rs:156-189` | 6 点二阶差分猜曲线族（`(max−min)/max < 0.02` 阈值 OCCT 无）+ 四个 conic 实体按采样重建参数 → `GeomToStep_MakeCurve.cxx:50-104`（`IsKind` 分派）/ `GeomToStep_Make{Circle,Ellipse,Parabola,Hyperbola}` 精确值 | 5 | **done**（2026-09-20）：`classify_curve`+`CurveKind` 删除；`emit_curve_entity` 按 Line→Conic→Trimmed→Bounded 的 `IsKind` 顺序分派（trimmed 走基曲线递归）；`write_conic_params` 按 MakeConic 顺序；四个 `emit_*` 改为取 `gp_*` 的精确 placement/半径/半轴/焦距；门禁全等于基线 |
| T-71 | **A3 派生（新）** | `occt-topo/src/step/write_context.rs::emit_curve_entity` 兜底臂 | `GeomToStep_MakeCurve.cxx:100-103` 的未识别曲线置 `done = false`（几何属性不写）未移植：端口曾写 B-spline 拟合（`fit_bspline_curve`）或起点切线 LINE 兜底 ⇒ 为 OCCT 拒绝导出的曲线伪造实体 | 5 后 | **done**（2026-09-20 第 48 轮）：① 删掉拟合/切线兜底，未识别**且非退化**曲线返回 `None` ⇒ `emit_edge` 按 `TopoDSToStep_MakeStepEdge.cxx:345` 写 `EDGE_CURVE(...,$,.T.)`（OCCT 那里拿到的是 `done=false` 留下的空 `theCurve`）；② **补齐 `Geom_BoundedCurve` 臂**（`MakeCurve.cxx:94-99` → `GeomToStep_MakeBoundedCurve.cxx:37-80`）：B 样条直接 `write_bspline_curve`（有理走 `rational`），Bezier 先按 `GeomConvert::CurveToBSplineCurve` 转成夹持 B 样条再写；为此给 `Curve` trait 加 `bspline_weights()`（默认 `None`，6 个实现者转发）；③ **退化边**（球极点那种"点曲线"占位）改走 `MakeStepEdge.cxx:263-330` 的"edge without 3d curve; creating"支（端口仍是以两端点连线的 LINE，OCCT 的非平面分支是采样拟合 B 样条 ⇒ 该子支 UNPORTED 已就地注明）；顺带删除死代码 `fit_bspline_curve`（T-62 余项的"写侧 n=8 采样重拟"）。**探针（已删）**：`Shape-2.step` 写回 82 条 `B_SPLINE_CURVE_WITH_KNOTS`、0 条未设几何；`ATU01038` 163 条 B 样条 / 391 LINE / 337 CIRCLE / **2 条 `$`**（真正的 `done=false` 用例）；三份文件写回后均能被 `read_step` 解析（无语法错误）。门禁：topo lib 1293/1、phase3 4/4、phase4 9/9、phase6 5/5、四道 STEP 门禁 14/14、13/13、11/11、2/3、`export_data_obj` 16/16 |
| T-40 | A4 | `occt-topo/src/algo_tools/construct.rs:165-205` | F 分支 32×32 投影 + **OCCT 没有的** "投影距离 > tol → Out" 门 → `BRepClass_FaceClassifier::Perform(F,P,Tol)`（`BRepClass_FaceClassifier.cxx:76-125`） | 4 | **done**（2026-09-20 第 53 轮）：F 分支按该函数逐行重写（`BRepTools::UVBounds` → `Extrema_ExtPS` 窗口 = `brep_uv_bounds::uv_box_of_face`，`IsDone`/`NbExt` 守卫 → `Unknown`，取最小平方距离，再交 2D 分类器）；**删掉自创的距离门**（OCCT 只按投影后的 UV 分类，不测 3D 距离）⇒ `algo_tools/tests.rs` 的 "面上方点 → Out" 期望按 OCCT 改为 `In`（附出处）；E 分支的投影本就是忠实 `Extrema_ExtPC`（`brep_extrema::closest_point_on_edge` 的 `_samples` 只是遗留形参）；同族 `algo_tools/queries.rs` 的两处 32×32 投影（`get_edge_off`/`face_normal_at_point`）改走忠实 `geom_api::project_point_on_surface`（T-52） |
| T-41 | A5 | `bop_builder_core.rs:79-96`、`bop_curved/face_meshing.rs:409-581`、`bop_curved/region_trim.rs:42-70,325-401` | 体素/网格布尔与计票 → 无 OCCT 对应 ⇒ **摘除并标未移植** | 4 | pending |
| T-42 | A6 | `occt-core/src/gcpnts.rs:30,87-225` | 中点二分 `MAX_DEPTH=16` 伪造 deflection → 调 `gcpnts_perform.rs:52`（已逐行移植） | **2** | **done**（2026-09-20） |
| T-43 | A7 | `occt-geom/src/extrema_pc/poly_roots.rs:563,621`、`extrema_cc/curve_curve.rs:83,214-239` | `clamp(24,256)` 网格 + 16 兜底 → 调 `extrema_pc/general_extrema_pc.rs`（`Extrema_GGExtPC`） | **2** | **部分 done**（2026-09-20）：点–曲线侧已忠实；`extrema_cc` 侧见 **T-66** |
| T-65 | A7 派生 | `occt-geom/src/extrema_pc/poly_roots.rs`（`ellipse_all`/`hyperbola_all`/`parabola_all`）、`extrema_pc/point_curve.rs::ext_pelc_all` | 三条解析臂已存在但未接线：缺 `myIsMin` 端口（`Extrema_ExtPElC.cxx:277`、`:381`、`:473`）⇒ 椭圆/双曲/抛物目前走 `default:` 臂；且 `Curve` 无 hyperbola/parabola 类型查询（T-12） | 2 后 | **done**（2026-09-20）：三条臂已接线，`myIsMin` 按 OCCT 逐行——椭圆比 `\|P−C(Us+0.1)\|²`（`cxx:277`）、双曲/抛物比 `\|P−C(Us+1)\|²`（`cxx:381`/`:473`）；类型查询阻塞已解除（A20/A29 轮补的 `gp_hyperbola`/`gp_parabola`），并给 `GeomTrimmedCurve` 补齐 `gp_{line,ellipse,hyperbola,parabola}` 转发（`GeomAdaptor_Curve::load` `cxx:252-254` 解包到基曲线）。门禁全等于基线；`occt-geom` 既有椭圆/双曲/抛物极值测试在解析臂下仍 151/151 |
| T-66 | A7 派生 | `occt-geom/src/extrema_cc/curve_curve.rs:40-115` | 通用曲线–曲线**种子集**自创（均匀网格 + 局部极值 + 边界最优点）→ OCCT `Extrema_GGenExtCC::Perform`（`Extrema_GGenExtCC.hxx:453-700`）用 `math_GlobOptMin` + `Extrema_GlobOptFuncCCC2` | 2 后 | **◐ 步 a done**（2026-09-20 第 54 轮）：`occt-math/globoptmin.rs` 补 `set_lip_const_state`/`set_continuity`/`set_functional_minimal_value`/`set_local_params`（+getter），并把 `compute_init_sol` 移到「全局盒建立之后、局部盒覆盖之前」，与 OCCT `SetGlobalParams`(ctor)→`SetLocalParams`→`Perform` 的次序一致（`cxx:112-171`、`Extrema_GGenExtCC.hxx:648-649`）；探针证得局部盒确实限制搜索范围。**步 b-1 done**（2026-09-20 第 55 轮）：新增 `occt-geom/src/extrema_cc/glob_opt_func.rs::GlobOptFuncCCC2` = `Extrema_GlobOptFuncCCC2` 的距离值 + 梯度 + Hessian（`GlobOptFuncCC.cxx:24-185` 的 3D 静态式，含 OCCT 的「值=距离、梯度/Hessian=平方距离之二阶导」原样保留），解析对拍逐位一致；**步 b-2 done**（2026-09-20 第 56 轮）：`GGenExtCC::Perform`（`Extrema_GGenExtCC.hxx:453-802`）逐段移植（区间集 + `ChangeIntervals` + Lipschitz 估计 + 每区间对 `SetLocalParams`/`Perform` + 递增最优与去重 + 平行检测），接线 `curve_curve_extrema_all(_range)` 并**删除自创种子集与其采样分类器**（`newton_curve_curve_all`/`build_samples`/`is_line`/`classify_circle`/`line_of_curve`/`circumcenter`）。**就地登记的 UNPORTED/偏差**：① 局部精化引擎仍只有 BFGS（OCCT 按 `myCont` 与函数类选 Newton/BFGS/Powell）；② `CellFilter` 的**容器**用接受点线性表（OCCT 是空间哈希），接受判据（`aCellSize` 内不重复）一致；③ 长度比测试用端口 `gcpnts::curve_length`（A15 偏差）；④ 端口无"闭合但非周期"曲线类，`IsClosed` 该支返回 false。**T-66 结项** |
| T-44 | A8 | ~~`occt-geom/src/convert_bspl.rs:326,349`~~（已删） | 采样折线 + 强制 1 次 → `GeomConvert_CompCurveToBSplineCurve.cxx:135-215`（每段先 `GeomConvert::CurveToBSplineCurve` → `Convert_*ToBSplineCurve` 精确有理表示） | 已清除自创件 | **批 66（第 81 轮）done（自创实现已删）**：`comp_curve_to_bspline` 及其单测已删除（零生产消费方），就地留 UNPORTED + 上述出处；OCCT 侧移植保留为待办，触发条件＝出现消费方；前置 `GeomBSplineCurve::{IncreaseDegree,InsertKnots}` 仍未移植 |
| T-45 | A9 | `convert/`、`cslib/mod.rs:24-41`、`gprop/mod.rs:27-42`、`bnd/obb_pca.rs`、`bnd/intersect.rs`、`geom/polyline_simplify.rs`、`int/curve_curve.rs` | 假出处（写 OCCT 包名但该包无此函数）→ 改标真实出处/非 OCCT（`polyline_simplify` 经 `intpatch_trace.rs:260` 改变交线几何，须标注） | 8 | **done**（2026-09-20，7 处文件头全部改正 + 包内容已逐项核对） |
| T-46 | A10 | `elib/surface_eval.rs:100,104`、`bnd/bsphere.rs:19-29`、`elib/intersect.rs:58-60`、`poly/make_loops.rs:234` | 静默默认值/凭空造值/只取首候选 → `ElSLib::SphereD2/TorusD2`、`Bnd_Sphere.cxx:73-96`、`Poly_MakeLoops.cxx:611-700` | 8 | **◐ 3/4 done**（2026-09-20）：`surface_d2` 补球/环面非零二阶导；`bsphere` 补"被包含⇒整体替换"分支 + `distance` 改回"到球心"（面距另立 `distances`）；`circle_plane_intersection` 删凭空造点并标 `UNPORTED`（`ElCLib` 无求交，忠实件是 `IntAna_Quadric`/`GeomAPI_IntCS`）；`choose_left_way` 标 `UNPORTED` ⇒ **T-73** |
| T-73 | **A10 派生（新）** | `occt-core/src/poly/make_loops.rs:234` `choose_left_way` | 只取首候选（`lst_ind_s[0]`）→ `Poly_MakeLoops.cxx:614-678`（3D 最小夹角）/`:692-738`（2D + `myRightWay`） | 8 后 | **done**（2026-09-20 第 47 轮）：`MakeLoopsHelper` 补 `get_normal`/`get_first_tangent`/`get_last_tangent`（3D）与 `get_first_tangent_2d`/`get_last_tangent_2d`（2D），默认 `None` = OCCT 基类返回 `false` 的语义；新增 `PolyMakeLoops::new_2d(helper, left_way)`（`myRightWay = !theLeftWay`，`cxx:682-687`）与 `two_d` 判别（端口把 `Poly_MakeLoops3D`/`2D` 合成一个结构体，已注明）。3D 支逐行移植（`Normal.CrossCrossed(TgtRef, Normal)` 投影、`AngleWithRef` 最小有符号角、`1e-4 - π` 归 PI、三处回落 `First()`）；2D 支逐行移植（`Poly_MakeLoops2D::chooseLeftWay`，`myRightWay` 取反 + `GpDir2d::Angle`）。**同批修 T-76**：`GpDir2d::Angle` 由 `acos(dot)`（无符号）改为 OCCT 的**有符号** acos/asin 分段（`gp_Dir2d.cxx:26-63`），并把 `is_normal`/`is_parallel`/`is_opposite` 改成 OCCT `gp_Dir2d.hxx:393-430` 的 `abs(Angle)` 形式。**探针验证（已删）**：3D normal=+Z、候选 −60°/−30° ⇒ 选 −30°（索引 3），normal=−Z ⇒ 符号翻转选索引 2，候选序 `[3,2]` 结果不变，无 normal ⇒ 回落首选；`GpDir2d::Angle` = +30°/−30°/180°/−135°，三个谓词仍正确；2D 左路（`theLeftWay=true`）⇒ 选 −60° 候选（索引 2），右路 ⇒ 取反后选 +30° 候选（索引 3） |
| T-47 | A11 | ~~`bop_builder_dispatch.rs:473-601`~~（已删）、`bop_draft_solid_occt.rs:41,49`、`edge_edge/edge_edge.rs:368-371` | 质心规则 / `solid−face` 原样返回 / 丢"顶点<3"的面 / 采样解当补集 → `BOPAlgo_BOP::CheckData`（`cxx:140-210`）、`BOPAlgo_Builder_3.cxx:329` 等对应控制流 | 4 | **子项 1+2 done**（2026-09-20 第 50 轮 / 第 46 轮）：① 子项 2：`bop_draft_solid_occt.rs::add_draft_face_once` 删掉自创的"顶点 < 3 即丢面"（`face_is_degenerate` 一并删除），与 `BOPAlgo_Builder::BuildDraftSolid`（`_3.cxx:329/342/357` 的无条件 `aBB.Add`）一致；② 子项 1：**`boolean_degenerate` 整套（含 `degenerate_empty`/`degenerate_non_solid`/`shape_centroid`）已删除**，非实体输入改由 `boolean_non_solid` → `bop_builder2::builder_bop_with_fuzzy`（= `BOPAlgo_BOP::Perform`+`BuildShape`），合法性只由已移植的 `bop_bop::check_data`（`BOPAlgo_BOP::CheckData`）判定；空实参按 OCCT 记 `BOPAlgo_AlertEmptyShape` 警告。**子项 3 done**（2026-09-20 第 57 轮）：删掉 `edge_edge::find_solutions` 里"把采样器 `inttools::edge_edge_intersections` 的解当补集并入"的那段（OCCT 完全没有这个动作），只保留"距离在容差内的极值即交点"这一替代判据；删除后**门禁逐项仍等于基线**（说明自 T-66 起忠实 `Extrema_ExtCC`/`GGenExtCC` 已覆盖这些用例）。真正的忠实件（bbox 递归 + `FindParameters`/`IsIntersection`/`CheckCoincidence`/`SplitRangeOnSegments`）未移植，已就地 UNPORTED + 行号并立项 **T-77**。**A11 结项** |
| T-77 | A11 派生（第 57 轮，新） | `occt-topo/src/edge_edge/edge_edge.rs::find_solutions`（现为替代判据） | OCCT `IntTools_EdgeEdge::FindSolutions` 的 **bbox 递归**未移植：`FindSolutions` 两重载（`cxx:290-549`）、`BndBuildBox`（`:1410-1419`）、`PointBoxDistance`（`:1423-1485`）、`FindParameters`（`:553-671`）、`IsIntersection`（`:1060-1146`）、`CheckCoincidence`+`DistPC`/`FindDistPC`（`:1150-1365`）、`SplitRangeOnSegments`（`:1366-1406`）、`MergeSolutions`（`:675-779`），以及 `Prepare` 的类型排序与 swap（`cxx:109-180`）与 `myRes*`/`myPTol*`/`myResCoeff*`（`:154-180`） | 后续 | pending（现状：替代判据「距离在容差内的极值即交点」已就地 UNPORTED + 行号）。**前置已备齐（第 58 轮）**：`BndBuildBox` 需要的 `BndLib_Add3dCurve::Add`/`GeomBndLib_Curve` 忠实件已在仓内（`geom_bnd_lib_curve3d::box_curve`），`int_tools_curve_box::add_curve_to_box` 的 48 点采样近似已改为调用它；`GeomAPI_ProjectPointOnCurve` 侧也有忠实件（`geom_api::project_point_on_curve`）。仍缺：`FindParameters`/`IsIntersection`/`CheckCoincidence`+`DistPC`/`FindDistPC`/`SplitRangeOnSegments`/两重载 `FindSolutions`/`MergeSolutions` 与 `Prepare` 的 `myRes*`/`myPTol*`/`myResCoeff*`/类型排序与 swap ⇒ **✅ done（批 80，第 93 轮，R2-1）**：上述全部已移植（新增 `edge_edge/find_solutions.rs`），替代判据与 `find_parameters`/`merge_solutions` 整体删除；门禁待"导出 obj"指令后跑。余项＝**R2-17**（OCCT 不存在的 circle/circle 快路径删除） |
| T-78 | A26 派生（第 58 轮，新） | `occt-topo/src/iges.rs`（写侧） | IGES 写侧是**部分替身**：只有 110（线）/100（弧）/120（`emit_face_surface` 里为球面**自造**回转面）/186（MSBO），曲线族用 6 点二阶导采样 `(max-min)/max < 0.02` 猜圆（`iges.rs:158-172`）且 `CurveKind::Other` 退化成弦线（`:389-395`）。OCCT 对应用 `GeomToIGES_GeomCurve`（按 `GetType()` 发 100/104/106/108/110/126）与 `GeomToIGES_GeomSurface`（球面发 **196**、锥 194、柱 192、环 198、旋转面 120），配合 `IGESData_*` 的参数布局与 144 裁剪面 | 后续 | pending（现状：两处已就地标 UNPORTED + 出处，待逐类补发射器；`A26` 其余项＝PLY 每面 node 块） | **步 1 done（第 70 轮批 55）**：曲线分派改为按精确类型（`GeomToIGES_GeomCurve.cxx:75-116`），删除 6 点 |C″| `<2%` 采样猜圆；圆支改用 `Geom_Circle` 的精确圆心/轴（`cxx:292-329`）。**余项**：104（椭圆/双曲/抛物）与 126（Bezier/B 样条）发射器、球面 **196**、柱/锥/环 192/194/198、以及 144 裁剪面；未实现的类型已就地标 UNPORTED + 实体号。 **步 2 done（第 71 轮批 56）**：① **卡片格式合规化**（数据 1–64/1–72、序号与节字母在最后 8 列、字母第 73 列、DE 指针、T 卡计数；出处 `IGESData_IGESWriter.cxx:769-943`）——改前端口写出的 IGES **不是合规文件**，此为本轮新发现；② **实体 126** 发射器（Bezier/B 样条，删除弦线替身）；周期曲线与裁剪区间仍 UNPORTED。**余项**：104（椭圆/双曲/抛物，`TransferConic`）、球面 196、柱/锥/环 192/194/198、144 裁剪面。 **步 3 done（第 72 轮批 57）**：初等曲面发射器 192/194/196/198 落地（出处 `cxx:1280-1435`），**此前柱/锥/环面被写成 108 平面**（几何错误）已修；球面自造 120 已换成 196；新增实体 123 方向。**余项**：128（B 样条曲面）、120/122（回转/拉伸面）、104（二次曲线族）、144 已写但边界仅取 wire 曲线。 **步 4 done（第 73 轮批 58）**：**实体 128（B 样条曲面）**落地（`IGESGeom_ToolBSplineSurface.cxx:64-125`），此前 B 样条面**全被写成 108 平面**；新增 `Surface::{bspline_surface_poles,uknots,vknots,weights}` 访问器（`GeomBSplineSurface` 实现、trimmed 转发、offset 不转发）。**余项**：120/122（回转/拉伸面）、104（二次曲线族）、128 的 `closedU/V` 精确判据（现按首末极点重合）与周期面 `SetNotPeriodic` 等价转换。 **步 5 done（第 74 轮批 59）**：**120（回转面）与 122（拉伸柱面）**落地（`GeomToIGES_GeomSurface.cxx:1000-1104`/`:1111-1188`，写出布局见 `IGESGeom_ToolSurfaceOfRevolution.cxx:119-127`、`IGESGeom_ToolTabulatedCylinder.cxx:90-98`），此前回转/拉伸面**全被写成 108 平面**；分派顺序按 `cxx:112-147`（Bounded → Elementary → Swept，裁剪面递归剥基面，`Geom_OffsetSurface` 不误判）；轴线按 OCCT 的"CAS.CADE 轴取反"写 110（`:1178-1183`）；122 的 U 取面自身 `Bounds`、准线按 `V1·Direction` 平移、终点取 `Value(U1,V2)`；曲线转移抽成 `emit_curve_range()`（`GeomToIGES_GeomCurve.cxx:94-126`）。同轮修正批 56 引入的 DE 卡 1 第 2 域（应为 **P 段首行行号**，`IGESData_IGESWriter.cxx:834`，原写成参数行数）。实测 120 数与回转面面数逐一对应（Shape 3/3、occ/top 4/4、occ/bottom 1/1）。**余项**：104（二次曲线族）、128 的 `closedU/V` 精确判据与周期面转换、144 边界。 **步 6 done（第 75 轮批 60）**：**104（椭圆/双曲/抛物二次曲线弧）＋ 124（变换矩阵）**落地（`GeomToIGES_GeomCurve.cxx:608-845`、`IGESGeom_ToolConicArc.cxx:103-125`、`IGESConvGeom_GeomBuilder.cxx:218-237`），此前这三类边退化成弦线；系数照 `gp_*2d::Coefficients` 的恒等二维坐标系取值＋调用点实参乱序与椭圆/抛物线翻倍、双曲线不翻倍的差异；形式号按 `ComputedFormNumber`（`eps=1e-8`）；非绝对坐标系记 124 于 DE 卡 1 第 7 域；顺带把 DE 的线宽/颜色/标签改成 OCCT 实体默认（0/0/空）。实测 `ATU01038` 104×63＋124×63、`occ/bottom` 104×2＋124×2。**余项**：椭圆整周支（`GeomConvert_ApproxCurve`，`cxx:620-645`）、128 的周期面反周期化、144 边界。 **步 7 done（第 76 轮批 61）**：**128 的 `closedU/V` 判据**改为 OCCT 的 `IsUClosed`/`IsVClosed` 语义（周期 ⇒ true；否则比较两条 iso 曲线的极点**逐分量** ≤ `Confusion`，有理时再加权重 ≤ `Epsilon(w)`，出处 `Geom_BSplineSurface_1.cxx:1350-1389`＋`Geom_BSplineCurve_1.cxx:662-`），并把写出的四界按**基面** `Bounds()` 裁剪（`cxx:245-255`/`:274-284`，基面因 `:492-515` 的递归）；周期支（`AdjustToPeriod`/`SetUOrigin`/`Unperiodize`，`cxx:256-343`）标 UNPORTED。实测 75 条 128 记录 A/B：前 9 域逐位相同，44 条末 4 域因裁剪收紧（`Shape-2` 一例 VMax 46.5004 → 26.0818＝基面上界）。**余项**：周期面反周期化与 `periodicU/V` 取值、椭圆整周支。 **步 8 done（第 77 轮批 62）**：**142 CurveOnSurface ＋ 144 TrimmedSurface** 按 `BRepToIGES_BRShell::TransferFace`（`cxx:250-405`）重建——144 由"`surface,u0,u1,v0,v1,1,m,曲线…`"（把参数域塞进边界类型/条数域、用普通曲线当边界）改为 `144, surface, outer_boundary_type, nb_inner_contours, outer_contour, inner…;`，边界改发 142（`142,0,surf,0,curve3d,2;`，2D 曲线 UNPORTED ⇒ 取 OCCT"仅 3D"支），多边 wire 的 3D 曲线发 **102** CompositeCurve；`isWholeSurface` 按 `BRep_Tool::NaturalRestriction` ＋ 平面/柱/锥强制 false（`:382-388`）复现，外侧 wire 用 `ShapeAnalysis::OuterWire` 忠实件。**余项**：2D 曲线、周期面反周期化、椭圆整周支。 **批 63（第 78 轮）**：查实 `write.iges.brep.mode` 默认 **0 = Faces 模式**，故根结构改为 `BRepToIGES_BRSolid`/`BRShell` 的 **144/142 ＋ 402 Group**（单元素不包组），删除模式 1 的 510/514/186 ⇒ **T-84 结项**。 **步 9 done（第 84 轮批 69）**：补上 `BRepToIGES_BRShell.cxx:334-365`——面内**不属于任何 wire 的边**各封一条内侧 142（仅 3D 支）；现有模型的面前均由 wire 组成，故对门禁/导出无影响（分支为忠实补全）。**T-78 余项**：2D（UV）曲线、周期面反周期化与 `periodicU/V`、整周椭圆、写入器可达性过滤。 **步 11 done（第 86 轮批 71）**：删除 `emit_shape` 的**无条件 116 point 块**（无人引用；OCCT 只在 192/194/196/198 的定位点等引用处写 116）；**可达性过滤**经评估需先做"参数指针结构化"，已立 **T-85**（含 10 处发射点清单），未实施。
| T-84 | A26/T-78 步 8 派生（第 77 轮登记，**第 78 轮批 63 结项**） | `occt-topo/src/iges.rs`（510/508/514/186） | 实体 510/508/514/186 只属于 **BRep 模式**（`write.iges.brep.mode = 1`），而端口根结构此前混用两模式 | 已结项 | **done**：经查 `write.iges.brep.mode` **默认 0 = Faces 模式**（`IGESData.cxx:90-94`），故按 `BRepToIGES_BRSolid`/`BRShell` 把根结构改为 **144/142 ＋ 402 Group**（`BRShell.cxx:411-476`、`BRSolid.cxx:100-168`、`IGESBasic_ToolGroup.cxx:104-115`），删除 510/514/186；BRep 模式（186/514/510/508/504/502）如需支持应作为**独立模式**另立任务（出处：`BRepToIGESBRep_Entity.cxx:107-211/395-528/535-657/663-731/738-831`） |
| T-85 | A26/T-78 步 11 派生（第 86 轮，新） | `occt-topo/src/iges.rs`（写入器） | 端口把**创建过**的实体全部写出，OCCT 的 `IGESData_IGESWriter::Write(root)` 只写从根可达的实体并按遍历序重编号 | T-85 步 1 done | **步 1 done（第 87 轮批 72）**：`Ent.refs` + `emit_refs` + 10 处复合发射点全部记录指针，`finish()` 的 `check_refs()`（debug）校验指针范围且全语料通过。**步 2 设计（第 88 轮给出，可直接照做）**：① 求根——`emit_shape` 的 Solid/Shell/Face 三个分支把 `emit_solid/emit_shell/emit_face` 的顶层 DE 收进 `IgesWriter::roots`（Compound 递归即多根）；② `reachable()`：从 `roots` 沿 `Ent.refs` DFS 得可达集（`refs` 已含全部 10 处复合发射点的指针）；③ **重编号**——`new[old]` 按创建序对可达集编号，写 DE/P 时用新号；④ **参数指针重写**：`params` 里指针与实数的区分规则＝`num()` 恒带 `.`/`e`（实数）、类型号恒在首位、其余**纯整数字段**按各类型已知布局即为指针（102：`n` 后 n 个；142：`creation_mode, surface, curve_uv, curve_3d, preference` 中的 2/4 位；144：1 位 surface、4 位 outer、5.. 位 inner、`flag`/`n` 为计数；402：`n` 后 n 个；192/194/196/198：point/axis 位与末位 refdir；120：axis/generatrix；122：directrix）——按该布局用 `refs` 的**记录顺序**逐位回填最稳妥（若嫌脆，可改为在构造 `params` 时对指针位留占位符 `\x00k`，写文件时替换）。收益＝去掉 `isWholeSurface` 面外侧 142 等游离实体（孤儿 116 已在批 71 单独清除）。**收益已量化（第 90 轮批 76，`iges_check` 的 unreferenced 统计）**：孤儿几乎全是 `142`（少数 `100`/`123`）——`Shape 8`、`Shape-2 28`、`linkrods 28`、`Offset 33`、`Shape-1 76`、`ATU01038 159`、`occ/bottom 184`、`occ/top 211`；`Cube/Cylinder/Cone/rev/Extrusion/HoledPlate` 仅 1（＝根 402，无孤儿）⇒ **✅ 步 2 done（批 83，第 93 轮）**：写出集合改为 `AddWithRefs` 可达集（先序 DFS，子序＝DE 字段 7 再参数指针）、按 DFS 访问序重编号、参数指针用 `#k` 占位符在写前替换（10 处发射点）＋ `trsf` 重映射；**验收 `iges_check` 待指令**（预期实体数下降、unreferenced 归零）。旁支：`point_entities` 为死字段（待清理） | | |
| T-79 | A5 余项（第 59 轮，新；**被 T-80 阻塞**） | `occt-topo/src/bop_curved/{face_meshing.rs:409-581,general_mesh.rs:205-290,region_trim.rs:42-70,325-401}` | **网格布尔**：`boolean_mesh`/`boolean_mesh_curves` 三角网格装配、`region_inside_other` 计票式区域分类、`general_boolean_trimmed` 通用曲面路径 —— OCCT 全部由同一套 `BOPAlgo_BOP`（`BRepAlgoAPI_*` / `BOPAlgo_Builder` + `IntPatch`）完成。端口 `curved_boolean_full`（`region_trim.rs:409-418`）是 `boolean_dispatch` 的三个分支之一，**门禁当前正走这条路** ⇒ 迁移需预留回归分诊：把 `boolean_dispatch` 直接指向 `bop_builder2::builder_bop_with_fuzzy`，逐项对照门禁 | 后续 | pending（`voxel_fallback` 已于第 59 轮摘除；`feature.rs::boolean_feature` 与 `brepfeat` 的 `voxel_boolean` 归 A12）。**第 60 轮实测：整体迁移被 T-80 阻塞**——把 `boolean_dispatch` 与 `curved_boolean_full` 直接指向 `bop_builder::boolean` 时，全部门禁与 lib 用例反而**全绿**，但探针显示曲面 operand 结果错误（见 T-80），故已回退；教训见 §9 |
| T-80 | **A1/A5 共同前置（第 60 轮，新）** | `occt-topo/src/bop_builder2/`（`BOPAlgo_BOP` 移植体） | **曲面 operand 下 `BOPAlgo_BOP` 移植体给错结果**（探针实测，2026-09-20 第 60 轮）：`FUSE box[-1,1]³ ∪ cyl(r=0.4,z∈[0,2])` → 只剩 box（6 面全平面、vol 8.0000，应收 8.5 且含柱面）；`CUT box − cyl` → box 不变（vol 8.0000）；`COMMON box ∩ cyl` → **空**（0 面）；同一引擎的**平面**用例 `FUSE box ∪ box(重叠)` → 14 面、vol 10.0000 **正确** ⇒ 缺口在曲面分支（`BuildRC`/面状态分类/分割），与 T-01（groove 平面精确布尔）同族 | 后续 | pending（**阻塞 T-79 与 A1 的最终迁移**：在修好曲面分支前，`boolean_dispatch` 的曲面分支必须继续走 `bop_curved` 的替代件）。**第 61/62 轮定位（探针，均已删；**订正**第 61 轮的错误结论**）：① 曲面求交没问题——`intpatch::intersection_curve_points(plane z=1, cylinder r=0.4)` 给出 4 条折线（半径 0.3999~0.4000、z=1.000000）；② 面–面阶段**也**没问题——在 `pave_ff_perform::perform_ff` 里 dump 每一对，`box vs cyl` 只有一对需要求交：`pair (11,36) kinds=(Plane,Cylinder) failed=false nb_curves=1`（**交线找到了**；第 61 轮据"带 pave block 的形状只有 1 个"下的结论是**错的**——`pave_blocks` 统计的是**边**的分块，一条节线只落在新生成的 section edge 上，与"FF 是否找到交线"不是一回事）；③ **GF 装配也没问题**——`perform_internal` 里 dump `self.result()`：`box+cyl → faces=9 kinds=[Plane×8, Cylinder]`（顶面被圆分割、柱面被 z=1 切开），对照平面参考 `box+box(重叠) → faces=16`；④ **断点在 `build_shape` 的实体装配**：同一 dump 在 `build_shape` 之后显示 `box+cyl → type=Compound faces=6 kinds=[Plane×6]`（**柱面与另两块平面被丢掉**），而平面参考 `16 → 14` 正常。**第 64 轮（本会话续做）定位到具体谓词**：`pave_ff_make::make_blocks_ff` 对该对走到 `is_valid_block_for_faces`（其忠实骨架 = `IntTools_Context::IsValidBlockForFaces`，`IntTools_Context.cxx:717-755`：取 `IntermediatePoint`，有 pcurve 时 `pc.D0(t)` → `IsPointInOnFace`，无 pcurve 时 `IsValidPointForFace`）→ **第一个面就返回 false**（插桩：`face=<key> has_pc=true uv=(2.715211,1.000000) reconstructed=false` → 落到 3D 分支 → `result=false`），于是节线块被丢弃、`mscpb` 为空（`total section edges before micro removal: 0`）；对照平面参考每对都是 `reconstructed=true → 2D 分支 → true` 并成功建边。⇒ **根因指向 FF 曲线的 pcurve 参数系**：该 pcurve 在 3D 曲线参数 `t` 处求值得到的 UV 不能重建 3D 中点（`uv.y=1.0` 恰是 3D 曲线中点参数，而柱面 UV 的 v 也应是 1.0、u 却对不上），即 **pcurve 与 3D 曲线参数不同系**（OCCT 的 `IntTools_Curve` 保证同参），**第 65 轮定位到生成处**：解析路径 `int_face_face_analytic.rs::plane_cylinder`（`IntAna_QuadQuadGeo` 的精确圆/椭圆/母线）→ `int_face_face.rs::curve_from_ic` → **`finish_curve` 里 `pcurve_of_curve(&curve, range, fa/fb)`**（`int_face_face_helpers.rs:277`）—— 它把 3D 曲线临时做成 edge（`[range.first, range.last]`）再走**通用投影** `make_pcurve_full`，得到的 2D 曲线是**按采样/裁剪参数**（`ShapeConstruct_ProjectCurveOnSurface` 式近似）而非与 3D 曲线同参；OCCT 侧 `IntTools_FaceFace` 的解析分支由 `IntPatch` 直接产出**与 3D 曲线同参**的 pcurves（`IntTools_Curve::SetCurves`），所以 `pc.D0(t_mid)` 必然落在 3D 中点。⇒ **修法二选一**：① 让解析臂的 pcurves 与 3D 曲线同参（按交点解析式直接构造 2D 圆/直线，而非投影采样）；② 在 `pcurve_of_curve` 里把 `make_pcurve_full` 的结果**重参数化到 `[range.first, range.last]`**（与 3D 同参）。判据：`box∪cyl` 探针的 `reconstructed` 变为 true、`make_blocks_ff` 建出节线 edge、`b.result()` 出现"被切面"。

**第 63 轮再订正（`build_solid` 插桩，已删）**：`build_rc(Fuse)` dump 显示 `b.result()` 就是 **2 个未分割的实体**——`build_rc(Fuse): b.result type=Compound children=2 explore(Solid)=2 kept=["Solid/**6f**", "Solid/**3f**"]`（6 面盒体 + 3 面柱体），于是 `mfs=3`、`sfs=3`（Plane,Plane,Cylinder）、`BuilderSolid areas=**0**`、最终只剩 `dmsts` 里那个未改动的盒体（6 平面面）；对照平面参考 `kept=["Solid/6f","Solid/6f","Solid/6f"]`、`mfs=16`、`sfs=14`、`BuilderSolid areas=1 ["Solid/14f"]` ✅。⇒ **真正的断点比批 39 说的更早**：FF 找到了节线（`nb_curves=1`），但**面从未被那条线切开**（GF 结果 = 两个原样实体）⇒ 目标在「FF 之后到面分割」这一段：`PostTreatFF`/`MakeBlocks` 建节线 edge（`pave_ff_post`/`pave_ff_make`/`pave_blocks`）与随后的 `FillImages`/`BuilderFace` 重建面。下一步：dump 该对的节线 edge 是否被创建、其 pcurve 落在盒体顶面/柱面上与否，以及 `FillImages` 后面 `b.result()` 的面数（现在恒等于原参数面数 ⇒ 分割没发生）|

| T-81 | **T-80 派生（第 69 轮，新）** | `occt-topo/src/pave_blocks/make_blocks.rs:641`（`make_pcurves_full` 的 trim 步）、`occt-topo/src/algo_tools/queries.rs:162`（`adjust_pcurve_on_surf`）、`occt-topo/src/pcurve_full/make_pcurve.rs:242`（`trim_pcurve_to_face`） | **`BOPTools_AlgoTools2D::AdjustPCurveOnFace` 的端口是替身**：OCCT 该函数（`BOPTools_AlgoTools2D.cxx:209-243` 两个重载 → `:247-400` `AdjustPCurveOnSurf`）只做**周期平移**（`GeomInt::AdjustPeriodic` 求 `du/dv = ±period`，特殊支路对圆柱面用 `dFi = tol/R` 兜底，末尾 `Geom2d_Curve::Translate`），**从不裁剪**；端口把它实现成 `trim_pcurve_to_face`（裁到面 UV 框 + 整周期平移）⇒ 周期面上缝边（seam）**只有一个 u=0 的表示**，OCCT 里经该平移会得到 u=2π 的第二表示。直接后果（第 69 轮探针）：`BOPAlgo_WireSplitter::SplitBlock` 的 `bIsClosed` 候选过滤（`BOPAlgo_WireSplitter_1.cxx:470-500`）比较 `Coord2dVf(候选边)` 与当前点 `aPb`，当走到 u=2π 时缝边候选算出的仍是 (0,·)，2π 的距离被拒 ⇒ 柱面**环带闭合不了** ⇒ 柱面不被节线切开 | 后续 | pending（**前置已备齐**：`BRep_Tool::Parameter(V,E,S,L)` 已按 OCCT 移植，见批 47）。下一步：按 `BOPTools_AlgoTools2D.cxx:247-400` 逐行实现 `adjust_pcurve_on_surf`（周期平移 + 柱面 `dFi` 支路 + `Translate`），替换 `pave_blocks/make_blocks.rs:641` 的 `trim_pcurve_to_face`；并在 seam 边上登记第二个表示（`GeometryRegistry::set_edge_pcurves` 已支持 forward→reversed 两元素）。判据：`box∪cyl` 的柱面被 z=1 切成两环带、FUSE 回到单一实体且体积 ≈8.5 | **done（第 67 轮批 51）**：`adjust_pcurve_on_surf` 已按 `BOPTools_AlgoTools2D.cxx:247-400` 重写为「整周期平移 + 柱面 `dFi` 支路 + 分类器复核 + `Translate`」，并接线到 `pave_blocks/make_blocks.rs::make_pcurves_full`；自造的 trim 退居工具函数并标 UNPORTED。判据「柱面被 z=1 切成两环带」**未达成**（平移不影响缝边收尾，仍在 T-83）。

| T-82 | **T-80 派生（第 69 轮，新）** | `occt-topo/src/bop_classify_occt.rs`、`bop_cells.rs`、`bop_build_common.rs`、`bop_build_solids.rs`、`builder_solid.rs`（面/壳 IN-OUT 状态与 growth 分组） | **曲面面的状态分类/保留**：第 69 轮补上柱面缝边的两个 pcurve（批 48）后，`box∪cyl` 的 FUSE 变成**单一实体 7 个平面面**（6 盒面 + 顶面内圆盘），**柱面自身的面被全部丢弃**——但柱面上半环带与顶盖在盒体外，FUSE 必须保留；同时盒体顶面的**内圆盘**（应属盒体内部、不参与 FUSE 外壳）反被保留 ⇒ 分类（`BOPAlgo_BuilderSolid` 的 growth/`IsHole`、`BOPAlgo_BOP::BuildRC` 的状态过滤）在这一类面上给错状态 | 后续 | pending（**取证**：探针 dump 逐面结果，见批 48；`warn=[]`、`shells=1`、`vol=9.3333` 亦不可信 —— 该壳不闭合）。下一步：对照 `BOPAlgo_BuilderSolid.cxx` 的 `PerformShapesToAvoid`/`PerformLoops`/`PerformAreas`/`IsHole` 与 `BOPAlgo_BOP.cxx:583-711` 的 `BuildRC`，先 dump 每个被切面的状态与 growth 归属，再判是分类器（`BRepClass3d_SolidClassifier`）还是分组逻辑 |

| T-83 | **T-80 派生（第 70 轮，新）** | `occt-topo/src/wire_splitter_block.rs::split_block`（`aVertMap`/`bIsClosed`）与 `vertex_parameter_on_face` | **闭合边的顶点参数与 `IsClosed` 的粒度**：OCCT 的 `BOPAlgo_WireSplitter::SplitBlock` 里 ① `MyDataMapOfShapeBoolean aVertMap` 以 **（TShape, Location, Orientation）** 为键（`TopTools_ShapeMapHasher`），即同一顶点的 F/R 出现是**两条不同记录**；端口 `vert_closed` 只按 TShape 合并。② `bIsClosed = Degenerated(E) || IsClosed(E, face)` 是**逐边**的，而 `IsClosed` 只对**该边实例**在闭合曲面上的表示成立 —— OCCT 的柱面只有 `ESTART`（u=0 侧）带闭合形式两 pcurve，`EEND`（u=2π 侧）没有（`BRepPrim_OneAxis.cxx:422-439`：`HasSides()==false` 时只有 `ESTART` 走 `SetPCurve(c1,c2)`）。③ `BRepPrim_OneAxis.cxx:404-408` 还会对闭合边调 `SetParameters(ETOP, TopEndVertex(), 0, myAngle)`，把顶点参数写进 `BRep_TVertex` 的点表示（`BRep_Tool::Parameter` 的 INTERNAL 支依赖它）。三者合起来决定 `Path` 的收尾判定 `anIsSameV2d`（`BOPAlgo_WireSplitter_1.cxx:447-467`）：端口第 70 轮已让两次 walk 正确走完两条环带，但**起点在 u=2π、终点落在 u=0**，`anIsSameV2d`（`aD2 < aTol2D2`）恒假 ⇒ `Path produced no wires` ⇒ 退回坐标链式拼接（柱面得到 6 个单边 loop / 1 个 3-wire 的假面） | 后续 | pending（**取证**：批 49 的 walk 轨迹与 LE/pcurve dump）。下一步：按上述三条对齐 —— 顶点闭合标志按出现朝向分开、`IsClosed` 按边实例判定、并为闭合边补点表示（顶点在 pcurve 上的参数）；判据：柱面被 z=1 切成两环带、`FUSE` 保留柱面上半段与顶盖 | **第 66 轮细化**：本项现在是 T-80 的**关键路径**（另见 §7「第 71 轮诊断」）：周期面未被切开 ⇒ 柱体 GF 实体四面全 `Internal` ⇒ FUSE 只剩盒体。两处待对齐：① 收尾判定（`aVertMap` 的朝向粒度、逐边 `IsClosed`、闭合边点表示）；② `perform_areas` 的 growth/hole 判据（`IntTools_FClass2d.cxx:548-565` 的逐 wire 有符号 UV 面积 `aS > 0 ⇒ growth`，端口 `fclass2d/classifier.rs:250` 现取「最大 |面积| 环」并按符号判 hole）——实测：仅禁用收尾判定可把柱面从 1 area/3 wires 变为 1 area/2 wires，但仍不是两环带 | **第 70 轮否证（实验，已回退）**：把柱面侧面缝边拆成**两条不同的边**（OCCT `StartEdge`/`EndEdge` 确为两条，`BRepPrim_OneAxis.cxx:968-1061`）后，`build_split_faces_occt` 对该面给出 **0 个 area**（LE 边数 6→10）⇒ 该结构假设**不成立**，缝边收尾仍是断点（判据仍是「柱面出 2 个 area」）。
| T-48 | A12 | `brepfeat/features.rs:109,419,361,508`、`feature.rs:197` | 体积解析覆盖 / 网格夹具冒充 `BRepPrimAPI_MakeCylinder` / `clamp(16,64)` → `BRepFeat_MakeDPrism/MakeRevol`、`LocOpe_Revol` | 已清自创重复件 | **批 67（第 82 轮）**：零消费方的 `feature.rs`（网格圆柱夹具 `z_cylinder_mesh` + 其 `resolution_for` 的 `clamp(16,64)`，即 T-48 引用的 `feature.rs:197`）已删除；余项＝`brepfeat` 在用的解析体积覆盖与 `brepfeat/features.rs:501` 分辨率（前置 T-80）| **第 69 轮复测（批 54 同轮）**：把 `boss_thru_all` 的解析体积覆盖（`brepfeat/features.rs:107-109`）换成结果形状的 `solid_volume` 后，`brepfeat::tests::boss_thru_all_pierces` **仍失败**（lib 1291/1 → 1291/2）⇒ 该项**仍被 T-80 阻塞**（精确布尔尚不能闭合穿孔拓扑），已回退；不必重复试。 **批 67（第 82 轮）**：T-48 引用的 `feature.rs:197` 所属模块（`protrusion/pocket/boss/hole` + 网格圆柱夹具 `z_cylinder_mesh` + 其 `resolution_for`）经核对**仓内零调用**，已整模块删除（`pub mod feature` 同步移除）；余项＝`brepfeat` 在用的解析体积覆盖（`:107-109`）与 `brepfeat/features.rs:501` 的分辨率（前置 T-80）。
| T-49 | A13 | `occt-topo/src/wireframe.rs:392-415` | 未裁剪 UV 窗口规则网格 → `BRepMesh_FaceDiscret` 按 pcurve 边界离散 | 5 | **已尝试 → 回退，被 T-68 阻塞**（2026-09-20）：`face_to_triangles` 已改为委托 `incremental_mesh_to_shape_mesh`、并删掉 `build_shape_mesh_wireframe`/`wireframe_face_triangulation`/`WIREFRAME_FALLBACK_RATIO_MAX` 与 `discretize_face` 的两处回退（编译 exit 0，无递归），但 `occt-topo --lib` 由 1293/1 变 **1284/10**（9 个新失败全部落在"面的边界结构缺失"上，见 T-68）⇒ 按失配即停回退，基线恢复 |
| T-68 | **A13/A18 前置（新）** | `occt-topo/src/primitives.rs`（sphere/torus）、`model_builder/*`、构造/模式类（`brep_pattern`、detach/copy 路径） | 忠实管线要求面具备**边界 wire + pcurve**；实测（临时探针，已删）：`sphere` faces=1 **wires/face=[0]** → 网格化失败 `DelaunayNodeInsertionMeshAlgo::perform: face 0 has no boundary UV points`；`torus` 同样 **[0]** → 失败；`cylinder` [1,1,1] → OK 52/48、`cone` [1,1] → OK 26/24、`box` ×6 [1] → OK 24/12。另 `wireframe::tests::face_with_hole_triangulates_ring_area` 变成"ring area 4 vs expected 3.8037（内环未生效）" ⇒ 内环/内 wire 也没进模型 | 5 前置 | **球面 done**（步骤 1+2）；**步 3 结论（第九轮）**：Delaunay 层逐行核对**全部忠实**，缺口根因在其上游"重叠共线前沿链" ⇒ 移出本卡、立项 **T-69**；环面/带孔内环待做 |
| T-69 | **T-68 步 3 派生（新）** | `occt-topo/src/meshing/node_insertion.rs`（`collect_boundary_uv`/`init_data_structure`/`finish_mesh`）、`model_builder/wire_builder.rs`（`add_wire`/`visit`）、`shape_tool.rs`（`visit_face`） | 带孔面（`wires=2`）的前沿链里存在**沿同一直线重叠**的链接（实测两端点各差 ~0.05、v 完全相同）⇒ `meshPolygon` 修正循环被 `Glued` 大批删段、交出拼接多边形 ⇒ `decomposeSimplePolygon` **正确地**判"无耳"并清空 ⇒ 该面 0 三角（166 面）。需与 OCCT `BRepMesh_NodeInsertionMeshAlgo`/`BRepMesh_ShapeTool`/`BRepMesh_ModelBuilder` 逐行对齐，找出"重复/偏移插入"的来源 | 5 前置 | **◐ 第 1 轮诊断 done**（2026-09-20 第 41 轮）：成因已定位（不再"来源不明"）——166 面 = ①`check_pcurves_and_shift` 的 `fix_lacking_all` 给"STEP 单闭合边 loop"（T0M 596 个）追加一条**重复闭合边**（362 条 wire `1→2`），两链几乎重合但参数化不同（=`~0.05` 偏移的真相）⇒ `Glued` 删段 ⇒ 断链；②跳过 FixLacking 后仍余 86 面，其单边 wire 走 `BRepMesh_ModelHealer.cxx:452-457`（端口忠实）把 pcurve 首点 UV 覆盖为末点（`u=0→2π`），形成跨 2π 的**回绕链**（UV 零面积）⇒ 清空。**关键反证（部分成立，已补测收敛）：`data/occ-ref/T0M.obj` 里"单闭合边管面"的侧壁也没有三角**（逐面写顶点不去重，20958/51160 条重复 `v`；失败管顶圆只出现一次=36 点且只被一个圆盘盖面使用）；但同日删逐面回退的实测被 `step_obj_parity` **13/14** 挡下（T0M `min[2]` 短 0.33，参考另覆盖 `(−0,−11.823,−424.742)`）⇒ 端口失败面集合是 OCCT 失败面集合的**超集**。**第 2 轮已把"超集里 OCCT 会网格化的那类"定位到 Torus face 20**（z 最小点归属面；同面两条 wire 的 pcurve 相差一个周期 u 窗口，而成功的 face 8 两条 wire 同窗口）⇒ 下一步：判据性实验 + 找 OCCT 里"跨 wire 周期对齐"的那段控制流 |

**T-68 执行方案（2026-09-20 现场核实，供下一轮直接照做）**

- **OCCT 结构**：球/环都是 `BRepPrim_OneAxis` 的整周旋转体，其 lateral wire 由 `LateralWire()`（`BRepPrim_OneAxis.cxx:660-679`）装配：`TopEdge()`（极点退化边，`cxx:1185-1232`，退化时 `MakeDegeneratedEdge`）→ `EndEdge()`（u=2π 处的经线，`cxx:1021-1061`）→ `BottomEdge()`（另一极点退化边，`cxx:1236-1279`）→ `StartEdge()`（u=0 处的经线，`cxx:968-1012`）。**极点退化边的 pcurve 是 v=±π/2 上的 u 等参线（跨整个 2π）**，正是 `BRepTools::UVBounds` 得到 u 跨度 2π 的来源；仅加两条经线（都在 u≈0）会让 UV 盒退化成 0 宽。
- **球**：`BRepPrim_Sphere::SetMeridian()`（`BRepPrim_Sphere.cxx:71-85`）——3D 经线 = `Geom_Circle(gp_Ax2(loc, **−Y**, XDir), r)`（参数化 u=−π/2→南极点、(r,0,0)、+π/2→北极点，与 port 现有 probe 一致）；**2D pcurve = `Geom2d_Circle(gp_Ax2d((0,0), XDir), r)`**（半径 r 的圆！）并调用 `SetMeridianOffset(2π)` 把经线参数裁到 `[3π/2, 5π/2]`。⚠️ pcurve 用半径 r 的圆而非 UV 直线，是为了与 3D 经线**同参数（SameParameter，弧长 = r·Δ角度）**——下一轮必须先按此理解核对 UV 映射（建议用 DRAWEXE 打印球面的 pcurve 与 `UVBounds` 对照，避免照抄错一半）。
- **环**：`BRepPrim_Torus::SetMeridian()`（`BRepPrim_Torus.cxx:71-82`）——3D 经线 = `Geom_Circle(gp_Ax2(loc + major·XDir, −Y, XDir), minor)`；2D pcurve = `Geom2d_Circle(gp_Ax2d((**major**, 0), XDir), minor)`。
- **仓内已有可直接用的件**：`GeometryRegistry::set_edge_pcurve(&TopoShape, face_key, Arc<dyn Curve2d>)` / `set_edge_pcurves`（`tgeometry.rs:310-320`，seam 边用两元素 forward→reversed）、`set_pcurve_range`（`:325`）、`occt_geom2d::{Geom2dLine, Geom2dCircle}`、`TopoBuilder::{make_edge, add_edge_vertices, make_wire, make_face}`。
- **验证链**：① 临时探针（sphere/torus `wires/face > 0` 且 `incremental_mesh_to_shape_mesh` 成功、顶点/三角数 > 0）；② `--lib` 1293/1 不回退（球面测试还断言 `UVBounds = 2π × π`，见 `primitives.rs:576-578`，即退化极点边的 u 等参线必须到位）；③ 再重放 A13/A18；④ 再重放 T-59/A23。
| T-50 | A14 | `geom/`（csg/delaunay/triangulate/fit*/polygon_*）、`elib/measure.rs`、`validate.rs`、`hlr.rs`、`viz_scene/`、`draw/`、`xcaf/`、`render_svg.rs` | 整包非移植件（`validate.rs` 的 `is_valid()` 恒 true 风险最高）→ 模块头声明未移植 + 真实出处 | 8 | **done**（2026-09-20）：`validate.rs` 批 8 声明"恒 true 不得作门禁"；批 13 补齐 10 处模块头 `UNPORTED` + 真实出处（`geom/mod.rs` 整包、`geom/csg.rs`、`geom/delaunay.rs`（**出处订正**：Bowyer–Watson 非 OCCT）、`geom/triangulate.rs`、`elib/measure.rs`、`hlr.rs`、`render_svg.rs`、`viz_scene/mod.rs`、`draw/mod.rs`、`xcaf/mod.rs`）；其余 `geom/*.rs` 由整包声明覆盖 |
| T-51 | A15 | `extrema_cc/curve_curve.rs:210-241`、`hyperbola.rs:16`、`parabola.rs:16`、`surface.rs:12-20`、`intana/line_torus.rs:574`、`extrema_ss.rs:69,119`、`curve_reparam.rs:37,189,228`、`gcpnts.rs` 积分 | 失败被 `fallback_*` 吞掉（丢 `StdFail_NotDone` 语义）+ 自创阈值 → 逐条按 OCCT 返回失败或补齐精确解 | 8 | **◐ 子集 done**（2026-09-20）：① `hyperbola.rs:16`/`parabola.rs:16` 的"假 `d2`（零二阶导）"已随 T-63 修掉；② `surface.rs:12-20` 的有限差分 `d2` 已由五个初等面的解析 `ElSLib::*D2`（含新增 `cone_d2`）取代，默认实现就地标 `UNPORTED` 并列出仍吃近似的 6 个实现者；③ `curve_reparam.rs:228` 的 `1e-14` 与 `occt-math/matrix.rs:316` 的 `1e-30` 均改为 OCCT `math_Gauss` 的 `MinPivot = 1.0e-20`（`math_Gauss.hxx:45-49`）；④ `basis_values` 标明实为 The NURBS Book A2.2（非 OCCT）。⑤ **`curve_reparam.rs:189` 的 `1e-15` 已修（批 27）**：`basis_values` 删除，改走新增忠实件 `occt-core/src/bspl/eval_basis.rs::eval_bspline_basis`（`BSplCLib::EvalBsplineBasis`，`BSplCLib_2.cxx:429-563`），退化判据为 `gp::Resolution()`、错误码 2 显式向上抛；**同时订正假前提**：OCCT 8.0.0 全树**无** `BSplCLib::BasisFuns`（老名字已移除）。**未完成**：`IsDone`/`myDone` 失败通道（`extrema_cc`/`extrema_pc`/`extrema_surf` 的网格+零值兜底，需改 API 为可失败类型）、`gcpnts.rs` 积分（`CPnts_AbscissaPoint`，前置 `math_GaussSingleIntegration` + `math_FunctionRoot`，两者仓内均未移植）。**批 15 补完 `Surface::d2` 的最后一个可解析实现者（Bezier，经精确 B 样条等价表示），并删除伪造的 `revolved.rs`（新 A31/T-74）⇒ 解析 D2 已覆盖全部 OCCT 面类，trait 默认仅余 `surface_fit`/`surface_to_grid` 两个非 OCCT 网格工具类（已列名 UNPORTED）** | **批 50（第 67 轮）**：`gcpnts.rs` 的 Simpson 已换成忠实 `CPnts_AbscissaPoint::Length` + `math_GaussSingleIntegration`（含 `order(C)` 分型阶数与 `GaussPointsMax()=61` 截断、13 轮区间加倍容差循环），Simpson 替代件已删；反解仍为 UNPORTED（缺 `math_FunctionRoot`）。
| T-52 | A16 | `occt-geom/src/geom_api.rs:52,97`、`occt-geom2d/src/curve_ops.rs:68-69` | 256×256 采样求交/投影 → `Extrema_ExtPS/ExtCC` + `intana2d`/`intimpargen` | 3 | **◐ 投影 half done**（2026-09-20 第 52 轮）：`project_point_on_curve` → 忠实 `Extrema_ExtPC`（`crate::extrema_pc`，T-43），`project_point_on_surface` → 忠实 `Extrema_ExtPS`（`Extrema_ExtPs`，T-67 步 2），两者按 `GeomAPI_ProjectPointOnCurve/Surf` 自身的 `IsDone`/取最小规则；删掉自创的 `refine_closest_curve`/`closest_params_in_window`/`surface_closest_params`（≈124 行）。**求交 half 未做**：`curve_surface_intersections`/`curve_curve_intersections` 已标 UNPORTED（真出处 `GeomAPI_IntCS`→`IntCurveSurface_Intersection`/`IntPatch_Intersection`；曲线–曲线为 `IntTools_EdgeEdge`/`IntCurve_IntConicConic`，8.0.0 无 `GeomAPI_IntCC`），它们被 `inttools/intersections.rs` 生产调用 ⇒ **✅ done（批 81，第 93 轮，R2-2）**：4 处曲线–曲线调用改走 `edge_edge::EdgeEdge`、1 处非平面边–面改走 `intcurvesurface::perform_curve_surface`，`geom_api` 三个采样器 + 4 个私有 helper + 3 个自带单测整体删除；2D 采样族（`occt-geom2d/curve_ops.rs` 的 `curve2d_intersections`/`curve2d_closest_point`，被 `geom2d_api` 使用）另立 **R2-18** |
| T-53 | A17 | `occt-topo/src/brep_exchange.rs:56-98` | 偏转形参被 Prs3d 顶掉 + 三级静默回退 → `RWObj_CafWriter`/`RWMesh_FaceIterator.cxx:87-89`（不网格化、空则跳过）⇒ 删回退、形参改名 | **6** | **done**（2026-09-20） |
| T-54 | A18 | `occt-topo/src/wireframe.rs:257-380` | 平面耳切 + 质心角度排序 + 桥洞 → 约束 Delaunay（`BRepMesh_DelaunayBaseMeshAlgo` + `BRepMesh_Delaun`） | 7 | pending |
| T-55 | A19 | `meshing/incremental_mesh/discret_root.rs:369,460-467`、`wireframe.rs:407-408`、`brepmesh.rs:38,174-175,259` | `WIREFRAME_FALLBACK_RATIO_MAX=0.10` 失败率换算法 + 四叉树魔数 + `clamp(3,64)` → OCCT 无失败率阈值，逐面置 `IMeshData_Failure`（`BRepMesh_BaseMeshAlgo.cxx:52-62`） | 7 | **失败率阈值已删**（2026-09-20）；**逐面 UV 栅格回退仍在**（`wireframe_face_triangulation`，`incremental_mesh/discret_root.rs:539,460-473`）：T0M 实测 **166/1769 面**走它。**T-69 第 1 轮（第 41 轮）新证据**：① 这 166 面里"单闭合边 loop 的管面"在 OCCT 参考网格里也**没有**网格（参考逐面写顶点不去重；失败管顶圆只出现一次、只被圆盘盖使用）；② 但同日**删回退实测**被 `step_obj_parity` 13/14 挡下（T0M `min[2]` 短 0.33 ⇒ 参考另覆盖 `(−0,−11.823,−424.742)`）⇒ 还有一类失败面 OCCT **会**网格化，回退暂不可删。下一步：定位该类面并补其控制流（详见 §7 T-69 第 1 轮） | 四叉树魔数（`brepmesh.rs`）属 A13/A18 遗留（T-49/T-55 交叉），未动 | **第 70 轮复测（仍阻塞）**：把逐面回退（`incremental_mesh/discret_root.rs::wireframe_face_triangulation`）整体禁用后复跑 `step_obj_parity`，仍是 **13/14**（红 = `occ_test_model_bboxes_match_occt`，T0M/occ 模型 bbox）⇒ 「端口失败面集合是 OCCT 超集」的判断**未过期**，A19/T-55/T-59/A13/A18/A23 仍需先做 T-69。
| T-56 | A20 | `occt-topo/src/step/read_geometry.rs:508-604` | `p1.distance(p2) < 1e-3` 替代 `V1.IsSame(V2)`；`Other`/HYPERBOLA 边域落 `(0,1)` → `StepToTopoDS_TranslateEdge.cxx:438,443` + `ShapeAnalysis_Curve.cxx:376-400` | **6** | **done**（2026-09-20）：`ProjectAct` 解析精确臂补齐（cxx:355-477）+ `clib::{ellipse,hyperbola,parabola}_parameter`/`parameter_{elips,hypr,parab}` + `Curve::{gp_line,gp_hyperbola,gp_parabola}`；`Project`/`ProjectAct` 去假失败通道；`V1.IsSame(V2)`→`GetCartesianPoints` 已补；**三个自创回退全部删除**，改为 cxx:442-444 的两次 `Project`。根因是新增项 **T-70/A29**（直线解析臂参数系平移），修好后门禁全绿 |
| T-70 | **A29**（T-56 根因，新） | `occt-geom/src/extrema_pc/point_curve.rs:44-52` | 直线解析臂把 `Extrema_ExtPElC` 用的 `gp_Lin` 建在 `d0(uinf)` 上 ⇒ 参数整体平移 `-uinf`（上游无界窗口下界为 −2/−1，故 `Project` 结果偏 +2/+1） | **6** | **done**（2026-09-20）：改用曲线自身 `gp_line()`（`Extrema_GGExtPC.hxx:390-405` 用 `theCurve.Line()`），无 `gp_line` 时回退 `d0(0)` 原点重建；实测 `OffsetPlaneHoleEdge` 区间 `(2,11)→(0,10)`、面积 204→280 |
| T-57 | A21 | `occt-topo/src/meshing/model_healer.rs:279-290` | 退化支路左右端接反、丢 `aPrevSqDist - aNextSqDist` 判定 → `BRepMesh_ModelHealer.cxx:491-512` + `hxx:143-151` | **6** | **done**（2026-09-20） |
| T-58 | A22 | `occt-topo/src/step/read_geometry.rs:348-370,740-760` | `ProjectAct` 缺 Ellipse/Parabola/Hyperbola 精确臂；圆用三点外心回退 → `ShapeAnalysis_Curve.cxx:382-400,160,200` | 6 | **done**（2026-09-20，随 T-56 批 2）：`ProjectAct` 已按 `cxx:355-477` 补齐解析臂（含 Ellipse/Hyperbola/Parabola/Line/Circle），`edge_params_for_curve` 的三点外心/`classify_curve` 采样族已全部删除；`step/read_geometry.rs` 内不再有 `classify_curve` 调用 |
| T-59 | A23 | `occt-topo/src/wireframe.rs:462-599` | 9 点 pcurve 采样当 UV 包围盒 → `BRepTools.cxx:172-330` `AddUVBounds`（精确，B-spline 走控制多边形） | 5 | **已尝试 → 回退，被 A13/A18 阻塞**（2026-09-20）：按 `BRepTools.cxx:172-367` 完整移植（`box_curve2d` 精确盒 + B-spline 周期验证 2/3/6 点 + 非周期钳制）后 `step_obj_parity` 14/14→**13/14**：`data/occ/T0M.stp` 的 **bbox min[2] ours=-424.978671 occ=-424.741876（Δ=0.237）**——忠实窗口等于 OCCT 的 `BRepTools::UVBounds`，但**消费方**是本仓自创的"UV 矩形栅格"建网格（A13/A18），而 OCCT 的网格由 pcurve 驱动（`BRepMesh_FaceDiscret`），所以凸包级别的窗口外扩不会漏进 OCCT 的网格。⇒ 必须先做 A13/A18（或改为逐样本判定），再重放本改动 |
| T-60 | A24 | `meshing/range_splitter/param_set.rs:86-133` | 周期标志 + 半径采样猜面型 → `GetType()` 分派（`BRepMesh_FaceDiscret.cxx:112` + `MeshAlgoFactory.cxx:64`） | 5 | **done**（2026-09-20）：改为 `GeomAdaptor_Surface::Load`（`cxx:422-513`）的 `DynamicType` 精确顺序（RTS→Plane→Cylinder→Cone→Sphere→Torus→Revolution→Extrusion→**Bezier**→BSpline→Offset→Other），删除 `match (up,vp)` 与自创 `is_cylinder_like`；新增 `Surface::is_bezier_surface()`。门禁全等于基线 |
| T-74 | **A31**（批 15 派生，新） | `occt-geom/src/revolved.rs::GeomRevolvedSurface` | 伪造的回转面（`d1` 的 `dv` 取 basis 点坐标、`transform` 空实现），与忠实件 `GeomSurfaceOfRevolution` 重复且**零消费者** | **done**（2026-09-20）：删除该文件与 `lib.rs` 的 `mod`/`pub use` 接线；四 crate `--lib` + 四道 STEP 门禁全等基线 |
| T-75 | **A0 派生**（批 19，新） | `occt-geom/src/offset.rs`（`GeomOffsetCurve::d3`/`eval_dn`）、`occt-geom2d/src/offset.rs`（`Geom2dOffsetCurve::d3`/`eval_dn`） | offset 自身的 `EvalD3`（`Geom_OffsetCurve.cxx:342-380` → `EvaluateD3` `pxx:495-529` → `CalculateD3` `pxx:215-306`；2D `Geom2d_OffsetCurve.cxx:289-327` → `pxx:460-484` → `pxx:189-279`）与 `EvalDN`（`cxx:386-410` / `:332-356`，3 阶以上转发基曲线）未移植 ⇒ 两类 `d3` 仍走 trait 默认（`D3 = 0`） | 1 后 | **done**（2026-09-20 第 45 轮）：两侧 `CalculateD3` 逐行落地（含 `R7/R6` 两个稳定性分支、`isDirChange` 的 `D3.Reverse()`）、补 `GeomOffsetCurve::d3`（基曲线 `EvalDN(u,4)` + `AdjustDerivative(..., 4, ...)` 奇异支路）、两侧 `eval_dn`（1/2/3 阶走 `EvalD1/2/3`，更高阶**转发基曲线**，`cxx:409`/`:355`）；文件头 UNPORTED 清单刷新。**解析验证（临时探针，已删）**：circle r=2 沿 +Z 偏移 0.5 ⇒ `off.d3(u)` 与解析同心圆 r=2.5 的 D0/D1/D2/D3 逐分量差 ≤ **1.11e-16**（4 个 u），`eval_dn(4)` 与基曲线同值；`occt-geom --lib` 151/151、`occt-geom2d --lib` 72/72、`occt-topo --lib` 1293/1、四道 STEP 门禁 14/14、13/13、11/11、2/3 |
| T-76 | **A10 派生 / gp 缺陷**（第 47 轮，新） | `occt-core/src/gp/dir2d.rs:14,16-18`（`GpDir2d::{angle,is_normal,is_parallel,is_opposite}`） | `angle` 为 `acos(dot)` = **无符号** [0,π]，且三个谓词按无符号写；OCCT `gp_Dir2d::Angle`（`gp_Dir2d.cxx:26-63`）是**有符号** ]−π,π]（`acos`/`asin` 分段），谓词在 `gp_Dir2d.hxx:393-430` 先取 `abs(Angle)` | 1 后 | **done**（2026-09-20 第 47 轮）：改为 OCCT 的有符号分段；谓词改成 `abs(Angle)` 形式（等价于原语义，故门禁不变）。探针：`Angle` = +30°/−30°/180°/−135° ✓，`is_parallel(+170°)/is_opposite(−179°)/is_normal(90°)` 均 true ✓。影响面：`GpDir2d::angle` 的仓内调用者仅三谓词与 `GpLin2d::angle`（OCCT `gp_Lin2d::Angle` 本身也是有符号）⇒ 全 5 crate `--lib` 与四道 STEP 门禁全等于基线 |
| T-72 | **A30**（批 4 派生，新） | `occt-topo/src/brep_surface.rs:78-101` | `brep_surface::classify_surface` = 8×8 采样 `is_planar(1e-6)` + 等距球心 + `1e-4*r` 阈值 ⇒ STEP 读入面型靠采样（`step/format.rs:160`、`step/read_topology.rs:1275+`） | 5 后 | **done**（2026-09-20）：改为 `GeomAdaptor_Surface::Load`（`cxx:422-513`）的精确类判定（RTS→Plane→Cylinder→Cone→Sphere→Torus→Other）；删除 8×8 采样与两个阈值。可观测修复：`step/read_topology.rs` 的 `GeomConvert_Units` 分派此前把柱/锥判成 `Other` ⇒ 柱面 pcurve **单位换算被跳过**。8 个既有测试的"缺陷断言"（3 处查找面谓词 `== Other` + 1 处 `assert_eq!(vanilla, Other)`）改为精确类型；未新增测试、未放宽门禁 |
| T-61 | A25 | `meshing/delaun/polygon_meshing.rs:346-375` | 自造"先删邻三角形再 AddElement" → `BRepMesh_Delaun.cxx:2263-2274` 失败即置 `IMeshData_Failure`，不改网格 | 7 | **◐ 主体 done**（2026-09-20）：删除整段删除补丁；`add_triangle`/`add_triangle_by_info` 改为可判定失败（链接满 ⇒ `failed=true` 返回 false = OCCT `Standard_OutOfRange`，经 `BRepMesh_BaseMeshAlgo.cxx:52-62` 空 catch 吞掉）；`mesh_polygon`/分解循环/`process_constraints` 前均检查 `failed`。门禁全等于基线（T0M 网格 46514/46945 → 46516/46962，见 §7 批 12）。**收尾项**：OCCT 是**整面**无网格，端口目前只中止当前多边形（已加三角形保留）⇒ 需把失败标志上抛到面级管线（`IMeshData_Failure`） | **批 52（第 68 轮）**：失败标志已上抛到面级——`Delaun::failed()` → `finish_mesh` 在 `erase_free_links` 之前提前返回 → `perform` 置该面 `MeshStatus::FAILURE`（OCCT `BRepMesh_BaseMeshAlgo.cxx:40-62` 吞掉 `Standard_OutOfRange` 后不执行 `commitSurfaceTriangulation`），`fast_discret` 同判据返回空网格；受门禁模型上该分支不触发（16 个导出计数逐位不变）。
| T-62 | A26 | `brep_exchange.rs:118,125`、`occt-core/src/io/{ply,stl}.rs`、`iges.rs:168,390-438`、`step/write_context.rs:16-17,43`、`vrml.rs:92`、`obj.rs` | PLY 焊接/属性类型、STL 阈值/头/嗅探、IGES 采样族与自造回转面、STEP 写侧采样重拟、`solid TRUE`、恒空 `vn` → 各 `RWPly_*`/`RWStl*`/`GeomToIGES_*`/`GeomToStep_MakeCurve.cxx:94-99`/`VrmlData_ShapeConvert.cxx:360` | 8 | **◐ 4/6 done**（2026-09-20）：① PLY 写侧属性类型 `uchar uint`（`RWPly_PlyWriterContext.cxx:214`）；② STL 退化法向改"平方量 > `gp::Resolution()`"（`RWStl.cxx:325/407`）；③ VRML `solid FALSE`（`VrmlData_ShapeConvert.cxx:356-361`）；④ **STL 头与嗅探**（批 16）：二进制头用 OCCT 字面量（`cxx:374-375`）、ASCII 头 `"solid \n"`（`cxx:304-305`）、嗅探只看前 5 字节（`cxx:523-531`）。**批 18 修 PLY 顶点焊接**（`brep_to_ply` 去 `weld_vertices`）。**未修**：IGES 采样族与整球回转面、STEP 写侧 B-spline 采样重拟。**订正**：OBJ `vn` 非写侧缺陷（写侧已有条件写 `vn`；缺口在网格侧不产法向 ⇒ A13/A18） |

> 执行纪律（本轮：**先提交、再任务化、再开工**）：commit `bcbc7dc` 已把审查前的全部工作树入库（8 个提交，工作树干净），此后每个修复单独成 commit，便于 A/B 与回滚。

**批 1（T-35 / A0）已完成 —— 2026-09-20**

- **改动**：`crates/occt-geom/src/offset.rs`（整文件重写）、`crates/occt-geom2d/src/offset.rs`（整文件重写）、`crates/occt-geom2d/src/curve.rs`（补 `Curve2d::d3`，对应 `Geom2d_Curve::EvalD3`）。
- **对齐内容**：`CalculateD0/D1/D2` 逐行移植（含 `R/R2/R3/R5` 与 `Dr/D2r` 两个稳定性分支、`theIsDirChange` 的 `D2.Reverse()`）；3D 法向 `Ndir = D1 ^ Direction`（沿法向偏移，**不再是沿方向平移**）；2D 法向 `(D1.Y(), -D1.X())`（**不再反号**）；`Continuity` 用 OCCT 的 C1→C0 / C2→C1 / C3→C2 表（`Geom_OffsetCurve.cxx:229-257`、`Geom2d_OffsetCurve.cxx:181-210`）；`Reverse` 改为**基曲线也反转**（`cxx:95-100` / `:90-95`）；`Transform` 3D 改为变换基曲线+方向+**带符号** scale（`cxx:454-460`），2D 保持 `abs`（`cxx:414-419`，与 3D 不同，OCCT 确实如此）。
- **容差**：用 `gp::Resolution()` 的忠实值 `precision::REAL_SMALL`（`DBL_MIN`），不用仓内 `precision::RESOLUTION`（`1e-12`，非 OCCT 值）。
- **验证**：① `occt-geom`/`occt-geom2d`/`occt-topo` `--all-targets` 编译 exit 0；② 门禁与基线**逐项一致**：topo lib 1293/1、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3、geom 153/153、geom2d 72/72；③ **解析恒等式探针**（临时，已删）：圆基（R=2）+ 方向 +Z + 偏距 0.5 ⇒ D0/D1/D2 与半径 2.5 的解析圆 **max 误差 0.000e0**（17 点），`Reverse` 后 `|P(0)|=1.5=R−d` ✓。
- **覆盖事实（重要）**：`data/*.step` 中**没有 `OFFSET_CURVE`** 实体（只有 `OFFSET_SURFACE`，属 `Geom_OffsetSurface`，非本项）；3D 类的唯一构造点 `step/read_topology.rs:1049` 故**无门禁覆盖**；2D 类**无任何构造点**（`Geom2dOffsetCurve::new` 无调用者，`geom_bnd_lib_offset2d::box_offset` 亦无调用者）⇒ **A0 是潜在正确性缺陷（latent），修复后不影响任何现有门禁数字**。审查报告 §3 中"2D 版被 `geom_bnd_lib_offset2d.rs:101` 消费"应订正为"该消费者自身为死代码"。
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
- **对齐内容**：`point_curve_extrema_all` 由"采样分类 + 自创网格路径"改为 `Extrema_ExtPC` 的真实结构——`Extrema_ExtPC` 在 V8_0_0 是 `Extrema_GGExtPC` 的别名（`Extrema_ExtPC.hxx:31-38`），其类型分派：line/circle（本仓可判定）走 `Extrema_ExtPElC` 解析臂并带上 OCCT 的 `myIsMin`（线 = min `cxx:77`；圆 = near 为 min、对径点为 max `cxx:177-188`），其余走 `general_extrema_pc` 的 `default:` 臂（`hxx:390-502`，`aMaxSample=17` + `DeflCurvIntervals`）。
- **删除的自创实现**：`poly_roots` 的 `newton_point_curve_all`/`build_samples`/`solve_f_zero`/`refine_seed`/`fprime`/`dist2`/`fval`/`pair` 与采样分类器 `is_line`/`circumcenter`/`classify_circle`（共 ~254 行），以及 `point_curve::param_for_point`。
- **语义订正（重要）**：原实现**无条件把区间端点当作极值**加入结果；OCCT 不会——`Extrema_ExtPElC` 只保留落在 `[Uinf, Usup]` 内的解（`cxx:180`），`default:` 臂仅在点与端点重合时补端点（`GGExtPC.hxx:474-502`），端点处理属调用方（`ShapeAnalysis_Curve::Project` `cxx:161-182`）。已按 OCCT 去掉该自创补端点。
- **验证**：`occt-geom --lib` **151/151**、`occt-geom2d` 72/72、`occt-topo --lib` 1293/1（唯一红仍是 `brepfeat::tests::groove_cuts_cylinder` = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（T-05）——**逐项与基线一致**；`--all-targets` exit 0。
- **未修（已立项）**：**T-66** = `extrema_cc` 的通用曲线–曲线种子集仍是自创网格（OCCT `Extrema_GenExtCC` 用 `math_GlobOptMin`，仓内 `occt-math/globoptmin.rs` 已移植，可直接用）；**T-65** = `ellipse_all`/`hyperbola_all`/`parabola_all` 解析臂缺 `myIsMin` 端口（`Extrema_ExtPElC.cxx:277`/`:381`/`:473`）且 `Curve` 无 hyperbola/parabola 类型查询（T-12），故这三种类型当前走 `default:` 臂（已在代码注释与 `extrema_cc/mod.rs` 头部写明）。
- **同轮试改并回退（T-64 / A27）**：把 15 处解析曲线/曲面的 `continuity()` 从 `3`(G2) 改成 OCCT 的 `6`(CN) 后，`step_to_obj` 由 13/13 变 **12/13**——**Sphere 网格变空**（round-trip 顶点数 0）。判定：**不能单独改**，某消费方把"报出的连续性"当跨度上限（嫌疑 `meshing/range_splitter/param_set.rs:295` 的 `continuity <= curve.continuity()`、`meshing/edge_discret.rs:1259` 的 `pcurve.continuity().min(surface.continuity())`）。已按纪律**回退**（`step_to_obj` 恢复 13/13、geom 151/151），并在 15 处就地标注 `// T-64: OCCT=GeomAbs_CN(6), blocked by consumer`、在 `geomabs.rs` 及本表记录依赖；下次须先对齐消费方的 `GeomAdaptor_Surface::NbUIntervals/NbVIntervals` 语义。

**批 3 前置核查（T-37/A1）—— 2026-09-20：审查前提被推翻，改为先移植**

- 审查报告 A1 的整改建议写着"仓内已有忠实件 `extrema_surf::point_surface_extrema(_box)`（`Extrema_ExtPS`）"——**该前提不成立**。实读：`extrema_surf/numeric_extrema.rs:74-119` 的 `point_surface_newton_all_box` 是 **24×24 网格播种**（局部极值 + 全局 min/max + 角点）+ Newton，Jacobian 用**数值差分**（`extrema_surf/analytic_solvers.rs:620-621`），模块头还写着 "Surface trait lacks d2"——而 `Surface::d2` **存在**（`surface.rs:12`，解析曲面覆写）。⇒ 它和 `brep_surface::surface_closest_params` 属**同一类替代品**；把 39 处调用点迁过去等于"自创换自创"。
- 忠实路径 = **移植** OCCT 的 `Extrema_ExtPS.cxx`(376 行) + `.hxx`(140)、`Extrema_GenExtPS.cxx`(1056) + `.hxx`(155)、`Extrema_ExtPElS.cxx`(454) + `.hxx`(76)：类型分派 + iso 退化处理（`IsoIsDeg`）+ 逐 C2 区间采样（`mySample`）+ `Extrema_GFuncExtPS` 的解析系统与 `math_FunctionSetRoot`，Jacobian 改用解析 `d2`。
- 已立项 **T-67**（A1 前置，T-37 被其阻塞），并把订正写回 `specs/_audit/_index.md` 的 A1 行与 §7 第 2 条、`extrema_surf/mod.rs` 模块头。**本批未改任何行为代码**，门禁不受影响。
- **迁移面盘点（44 处匹配）**：`crates/occt-topo` 内 **40 处**调用点（`brep_surface.rs` 自身 1、`brep_class3d.rs` 3、`inttools_range.rs` 3、`bop_build_common/build_ops` 2、`bop_build_solids` 3、`algo_tools/construct`+`queries` 3、`brepmesh.rs` 3、`intcurvesurface/solvers` 3、`pave_intersect/fill_ctx`+`vertex_face` 2、`brep_projection.rs` 2、`edge_face_kind.rs` 2、其余 1 处/文件：`bean_face`、`bop_build_faces/builder_like`、`brep_connect`、`edge_face`、`brep_faces`、`brepfeat/features`、`geometry_query`、`fclass2d/tests`、`fillet_curved/rolling_ball`、`int_curves_face`、`int_tools_full/context`、`shape_naming`、`wire_splitter_block`），以及 **`occt-geom/src/geom_api.rs:376` 自带的一份同名实现（属 A16）**。参数 `nu/nv` 从 8 到 64 不等——忠实版 `Extrema_ExtPS` **没有**调用方网格（`mySample` 由内部分辨率决定）⇒ T-67 完成后这 40 处需一次性改签名（去掉 `nu/nv`），属机械迁移；窗口语义（`ShapeAnalysis_Surface::ValueOfUV` 的 `uf-du..ul+du` 扩展盒）由 T-67 的 `Extrema_ExtPS::Initialize` 承接。

**批 8 先行（T-45 / A9，假出处整改）已完成 —— 2026-09-20**

- **逐项核对的 OCCT 事实**（不是照抄审查结论）：`CSLib/` 包只有 `CSLib`、`CSLib_Class2d`、`CSLib_NormalPolyDef` + 三个状态枚举（无点–面分类器）；`GProp/` 包只有 7 个类（`GProp_GProps`/`PGProps`/`SelGProps`/`VelGProps`/`CelGProps`/`PrincipalProps`/`PEquation`，无三角化属性类）；`Convert/` 包是 18 个 B 样条转换类（无坐标转换）；`Bnd_OBB` 恰好 3 个 ctor（空 / 中心+三轴+尺寸 / 由 `Bnd_Box`）；`Bnd_Tools` 只有 `Bnd2BVH` 两个重载；全树**无** `Ramer`/`Peucker`，`Douglas` 亦无命中（先前的 177 命中是 `Standard_ProgramError` 的正则误匹配）；全树唯一的 polar 助手是 `V3d_View.cxx` 的文件级 `toPolarCoords`。
- **改动**（仅注释）：7 处文件头改为"**非 OCCT 翻译（port-internal）**"+ 上述包内容清单，并给出真实 OCCT 对应（点–面分类走 `BRepClass_FaceClassifier`/`IntTools_FClass2d`；形状/网格属性走 `BRepGProp::SurfaceProperties/VolumeProperties` + `GProp_GProps`；曲线/边邻近走 `IntTools_EdgeEdge`+`Extrema_ExtCC`）。`polyline_simplify` 额外写明 `intpatch_trace.rs:260-262` 用它抽稀交线 ⇒ **非 OCCT 规则在改交线几何**。
- **验证**：`occt-core --all-targets` exit 0、`--lib` 290/290（纯注释改动）。

**批 3 / T-67 分步 1 已完成 —— 2026-09-20（`ExtPElS` 解析臂接线）**

- **核对结果**：`extrema_surf/analytic_solvers.rs` 的五条 `Extrema_ExtPElS` 解析臂**已存在且各有单测**（`point_plane_extrema`/`point_sphere_extrema`/`point_cylinder_extrema`/`point_cone_extrema`/`point_torus_extrema`），但**分派只接了 plane/sphere**（用采样分类器 `classify_plane`/`classify_sphere`），cylinder/cone/torus 臂从未被 live 路径走到。
- **改动**：`point_surface_extrema_all`（`extrema_surf/numeric_extrema.rs`）与 `curve_surface_extrema_all` 改用 **OCCT `GetType()` 等价查询**（`Surface::gp_pln`/`gp_sphere`/`gp_cylinder`/`gp_cone`/`gp_torus`，与 `Extrema_ExtPS::Perform` 的类型 switch 一致），五类初等面全部走解析臂，其余仍走通用路径（待分步 3 换成 `GenExtPS`）。
- **配套（`GeomAdaptor_Surface.cxx:423-425`）**：该处 `load` 会把 `Geom_RectangularTrimmedSurface` 解包成 basis + 范围 ⇒ 给 `rectangular_trimmed.rs` 补上 5 个类型查询的**委托**，否则裁剪后的平面/柱面会被当成"非初等"而落到通用路径。
- **验证**：`occt-geom --lib` 151/151、`occt-topo --lib` 1293/1（唯一红仍 T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3 —— **逐项与基线一致**（含裁剪委托单独复跑 `step_to_obj`）。
- **未做（T-67 剩余）**：`Extrema_ExtPS` 的**范围/`IsoIsDeg`** 分派（现有解析臂无参数窗口，故 `point_surface_extrema_box` 仍走通用路径）、`Extrema_GenExtPS` 主体（1056 行，含逐 C2 区间采样 + `math_FunctionSetRoot` + 解析 Jacobian）⇒ 完成后再做 T-37 的 40 处调用点迁移。**⇒ 前半项已在批 25（第 49 轮）完成**：`Extrema_ExtPS` 分派层 `ExtPs` 落地（窗口/`IsoIsDeg`/`TreatSolution`），`point_surface_extrema_box` 已改走它；`Extrema_GenExtPS` 引擎（实测 1195 行）仍待分步 3。

**批 6 先行（T-57 / A21）：`adjustSamePoints` 忠实化 —— 2026-09-20**

- **OCCT 事实（逐行核对）**：`BRepMesh_ModelHealer.cxx:489-512` 按 `aPrevSqDist - aNextSqDist > gp::Resolution()` 分两支，**两支都写「curr@prev ← prev_val」与「curr@next ← next_val」**，区别只在退化守卫（`hxx:143-148`：当两个配对指向**同一端**时把该端翻到另一头，并用 `closestPoint` 在**minor 边**的端点里重算）；而 minor 边在 prev 更远那一支是 **prev**，在另一支是 **next**。端口原实现**丢了这个判定**、只保留了一支的形态（永远用 next 的端点），且 `closestPoints` 的返回值（平方距离）没接出来。
- **同步订正的两个比较器**（`hxx:88-129`）：`closestPoint` 用**平方距离 + 严格 `<`**（平局取 second），`closestPoints` 的取舍是 `sq1 - sq2 < gp::Resolution()`（取 first，除非明显更远）——端口原实现用**线性距离 + `<=`**，平局方向相反。
- **改动**：`model_healer.rs` 新增忠实 `closest_point`/`closest_pair`（返回 `(a_side, b_side, sq)`）与 `adjust_same_points`（`hxx:134-152` 的指针翻转 + `closestPoint` 重算语义），`connect_closest_points` 尾部改为 OCCT 的两支；新增 `GP_RESOLUTION = REAL_SMALL`（`gp::Resolution()`，注明仓内 `RESOLUTION=1e-12` 不是 OCCT 值）；两处既有单测按新签名适配（语义断言不变）。
- **验证**：`occt-topo --lib` 1293/1（唯一红仍 T-01 `groove_cuts_cylinder`）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3 —— **逐项与基线一致**。

**批 6 续（T-53 / A17）：导出路径删掉两级静默换网格器 —— 2026-09-20**

- **改动**：`brep_exchange::export_mesh` 去掉 `.or_else(brepmesh::incremental_mesh)`（legacy 四叉树）与 `.unwrap_or_else(shape_mesh::mesh_shape)`（UV 栅格）两级回退——OCCT 的 OBJ 写侧 `RWObj_CafWriter` 读的是形状上**已有的**三角化（`RWMesh_FaceIterator.cxx:87`），缺三角化时**跳过**该面（`:89`），从不切换网格器。现在网格化失败即返回空网格（显式、可观测），不再静默换成非 OCCT 算法。
- **形参改名**：`deflection` → `maximal_chordial_deviation`（对外 7 个函数 + `prs3d_get_deflection` 同改），并在三处注释里写清：该参数是 drawer 的 `MaximalChordialDeviation`，**只在包围盒为空/无界时**生效；真正驱动密度的是 `Prs3d::GetDeflection = maxComp(bbox)*0.001*4`（`Prs3d.hxx:82-103`）。同步更新 `examples/export_data_obj.rs`、`tests/step_obj_parity.rs`、`tests/step_obj_area.rs` 的过时注释（原文仍称"Rust 用 UV 栅格"）。
- **验证**：`occt-topo --lib` 1293/1（唯一红仍 T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3 —— **门禁全绿且逐项与基线一致** ⇒ 对全部门禁模型，Delaunay 管线本身即可成功，被删的两级回退**从未被需要**（此前只是静默兜底风险）。

**批 5 先行尝试（T-59 / A23）：`face_uv_bounds` → 忠实 `AddUVBounds` —— 已回退，2026-09-20**

- **做了什么**：按 `BRepTools.cxx:172-367` 完整移植 `AddUVBounds`（替换 138 行自创启发式：9 点采样、`period_shift` 手工解缠绕、"仅 u 变化边"`1e-9` 判据、整周期检测），并配 `BRepTools.cxx:126-157` 的 face 级并集 + `:141-153` 的自然范围兜底。用仓内 `geom_bnd_lib_curve2d::box_curve2d`（= `BndLib_Add2dCurve::Add`）、`BndBox2d`、`make_pcurve_full`，含 B-spline 周期验证的 2 点（`IsUClosed`）与 3/6 点检查（`100*Confusion²`）及非周期钳制。编译 exit 0。
- **结果（回归）**：`step_obj_parity` **14/14 → 13/14**：`data/occ/T0M.stp` 的 `bbox min[2] ours=-424.978671054 occ=-424.741875692`（**Δ=0.237**，我方超出 OCCT 参考）。
- **归因**：移植后的窗口**就是** OCCT `BRepTools::UVBounds` 的值；差异来自**消费方**——本仓的 `wireframe::face_to_triangles` 是自创的"UV 矩形栅格"网格化（A13/A18），会把窗口内、face 之外的部分也采样进去；而 OCCT 的网格由 pcurve 驱动（`BRepMesh_FaceDiscret`），窗口的凸包级外扩**不会**进入网格。这也解释了原实现那些看似古怪的钳制/整周期启发式其实是在替栅格建网格兜底。
- **处置**：按"失配即停"**回退**该改动（`git checkout -- crates/occt-topo/src/wireframe.rs`），复跑 `step_obj_parity` **14/14** 确认恢复；把结论与依赖写进本表 T-59 与 `_index.md` 的 A23 行。**正确次序**：先做 A13/A18（pcurve 驱动建网格），再重放本改动；否则忠实窗口会被自创消费方放大成超界几何。
- **第二轮实验（同轮补做，证据加强）**：发现**忠实件其实早已存在**——`crates/occt-topo/src/brep_uv_bounds.rs` 是 `BRepTools::AddUVBounds` 的逐行移植（face/wire/edge 三个重载 + 周期 B-spline 探针，2D 盒走 `bnd_lib_add2d::add_geom2d_range`），而 `wireframe::face_uv_bounds` 是**另一份**自创实现。于是把 `wireframe::face_uv_bounds` 改为**委托** `brep_uv_bounds::uv_box_of_face`（20 行，删掉 138 行自创），重新跑 `step_obj_parity`：**仍然是同一个 delta** —— `occ/T0M.stp: bbox min[2] ours=-424.978671054 occ=-424.741875692`（与手写版逐位相同）。⇒ 两套独立实现给出同一结果，**排除"移植写错"**，坐实"忠实窗口更大 + 自创 UV 栅格消费方泄漏"的因果；也说明原自创实现（采样 + 钳制 + 整周期启发式）实际把窗口**压小**了，恰好让栅格网格不越界。
- **A13/A18 的具体修法（已核实入口）**：`wireframe::face_to_triangles` 直接改走 `crate::meshing::incremental_mesh::incremental_mesh_to_shape_mesh(&f.0, def)`（`incremental_mesh/discret_root.rs:1412`，接受任意 `TopoShape`，含单个面），返回 `ShapeMesh`；然后删除 UV 栅格（`:415-460`）、`planar_polygon_triangulate`/耳切/质心排序（A18）与 `face_uv_bounds` 的自创体。**预期影响面**：`shape_mesh`/`brepmesh`/`geometry_query`/`feature`/`brep_tools`/`brep_scene`/`viz_scene`/`hlr`/`render_svg`/`gltf`/`vrml` 的网格都会换成忠实管线，`shape_volume`/`shape_area` 类 lib 测试与 `step_geometry_parity`（T-05）**预计会移动**，需逐项分诊后单独成 commit。**未在本轮动手**（爆炸半径大、需预留分诊预算）。

**批 5（T-49/T-54，A13/A18）尝试与回退 —— 2026-09-20**

- **做了什么**（编译 exit 0、无递归）：① `incremental_mesh/discret_root.rs` 删掉 `build_shape_mesh_wireframe`（整形状 UV 栅格回退）、`wireframe_face_triangulation`（逐面栅格回退）、`WIREFRAME_FALLBACK_RATIO_MAX = 0.10` 失败率阈值与 `PendingFace` 的无用字段，`perform` 改为直接 `self.build_shape_mesh(&mut model)?`，失败面只置 `MeshStatus::FAILURE`（对齐 `BRepMesh_BaseMeshAlgo.cxx:52-62`）；② `discretize_face` 删掉两处回退并去掉重复的 `ModelPreProcessor::perform`（与 `perform` 管线对齐）；③ `wireframe::face_to_triangles` 改为委托 `incremental_mesh_to_shape_mesh(&f.0, def)`，删除 UV 栅格与平面耳切调用。
- **测量结果（分诊）**：`occt-topo --lib` **1293/1 → 1284/10**。9 个新失败：`vrml::write_vrml_shape_of_sphere_nonempty`、`viz_scene::{shading_modes_flat_vs_gouraud, textured_sphere_looks_textured, pick_two_shape_scene_hits_correct}`、`rwmesh::roundtrip_shape_mesh`、`draw::sphere_fillet_offset`、`bop_curved::sphere_inside_box_common_volume`、`brep_pattern::pattern_volume_sums_copies`、`wireframe::face_with_hole_triangulates_ring_area`。
- **根因（临时探针实测，探针已删）**：忠实管线要求面有**边界 wire + pcurve**，而本仓部分形状的面根本没有 wire：`sphere` faces=1 **wires/face=[0]** →`... face 0 has no boundary UV points`；`torus` 同为 **[0]**；对照 `cylinder` [1,1,1] → **OK 52/48**、`cone` [1,1] → **OK 26/24**、`box` ×6 → **OK 24/12**（说明改法与入口本身可行）。另有带孔面：忠实管线给出 **ring area 4 vs 期望 3.8037**（内环未参与）。
- **处置**：按"失配即停"**回退全部三处改动**（`git checkout`），复跑 `--lib` **1293/1** 确认恢复；把前置缺口立项为 **T-68**（按 `BRepPrim_Sphere`/`BRepPrim_Torus` 补 seam/退化边 wire + 查 `model_builder` 内环处理），T-49/T-54 标记为"被 T-68 阻塞"。**次序**：T-68 → 重放本改动 → 再重放 T-59（A23）。

**T-68 现场勘查（2026-09-20，本轮只做方案核实，未改代码）**

- 读通 OCCT 侧结构：整周旋转体（球/环）的 lateral wire = `TopEdge`(极点退化) + `EndEdge`(u=2π 经线) + `BottomEdge`(极点退化) + `StartEdge`(u=0 经线)，见 `BRepPrim_OneAxis.cxx:660-679 / 968-1061 / 1185-1279`；**球面 UV 的 2π 宽度来自两条极点退化边的 u 等参线**，只加经线不够。
- 两个非显然事实（下轮可直接用，避免踩坑）：① `BRepPrim_Sphere::SetMeridian` 的 2D pcurve 是**半径 r 的 `Geom2d_Circle`**（不是 UV 直线），目的是与 3D 经线同参数（弧长 = r·Δ角），并配 `SetMeridianOffset(2π)` 把经线裁到 `[3π/2, 5π/2]`；② `BRepPrim_Torus::SetMeridian` 的 2D pcurve 是**圆心 (major,0)、半径 minor 的圆**。⇒ 下轮务必先用 DRAWEXE 打印该面 `UVBounds` 与 pcurve 对照再照抄。
- 仓内 API 盘点（已确认可用）：`GeometryRegistry::set_edge_pcurve/set_edge_pcurves/set_pcurve_range`、`occt_geom2d::{Geom2dLine,Geom2dCircle}`、`TopoBuilder::{make_edge,add_edge_vertices,make_wire,make_face}`；现有球面测试还断言 `UVBounds = 2π × π`（`primitives.rs:576-578`），可作为实现正确性的即时报错点。
- **未改任何代码**（避免半成品）；T-68 方案已写入任务卡，下一轮从其"验证链"逐步执行。

**T-68 第一轮实测（2026-09-20）：只加经线不够，极点退化边是必须的**

- **实现尝试**：给球面加了忠实经线线（`BRepPrim_Sphere::SetMeridian` 的圆：`GpAx3(O, −Y, X)` 半径 r，参数 `[3π/2, 5π/2]`），两条边（u=0 与 u=2π）+ 用 `set_edge_pcurve` 挂上 `BRepPrim_OneAxis::LateralFace` 规定的 **UV 直线 pcurve**（沿 +v，起点 `(u, −2π)`）。编译通过。
- **探针实测（临时探针，已删）**：`wires=1` ✓、`uv_bounds=(0, 2π, −π/2, π/2)` ✓（与 `primitives.rs:576-578` 的断言一致）、`classify=Sphere` ✓、`analytic_surface_area=4π` ✓、UV 栅格网格 4096 顶点 ✓ —— 但 **`brep_gprop_full::surface_properties` 返回 0** ⇒ `--lib` **1291/3**（新增两个失败：`brep_gprop_full::{sphere_surface_volume, adaptive_box_sphere}`）。
- **根因（决定性）**：只有两条经线时，UV 环路是 **u=0 与 u=2π 两条竖线**，**环路包围面积 = 0** ⇒ `BRepGProp` 的 2D Gauss 积分（`FaceGauss`）得 0。OCCT 之所以在 `LateralWire()` 里还要 `TopEdge`/`BottomEdge`，正是因为它们是 **v=±π/2 上的 u 等参线（跨整个 2π）**，把 UV 环路闭合成完整矩形。⇒ **极点退化边不可省**。
- **新增前置**：需要在 port 里能构造**退化边**（无 3D 曲线、只有 pcurve 的边）。`EdgeGeom` 有 `degenerated: bool`，但 `TopoBuilder::make_edge` 目前强制要求 3D 曲线 ⇒ T-68 需先给 builder 加"退化边"构造（对齐 `BRepPrim_Builder::MakeDegeneratedEdge`，`BRepPrim_OneAxis.cxx:934/1212/1264`）。
- **处置**：按纪律**回退**（`git checkout`），`--lib` 复跑 **1293/1**；结论并入本任务卡。**次序**：退化边构造 → 球/环 wire（含极点边）→ 探针验证 → 重放 A13/A18 → 重放 T-59。

**T-68 第二轮（2026-09-20）：退化边构造已落地；完整 `LateralWire` 让忠实网格首次跑通**

- **已落地（保留，附加性、无调用者、门禁不动）**：
  - `TopoBuilder::make_degenerated_edge(&GpPnt)`（`builder.rs`）+ port 内 `DegeneratePointCurve`（常点曲线占位，导数全零、范围 `[0,0]`、`CN`），并置 `GeometryRegistry::set_degenerated`；语义对齐 `BRepPrim_Builder::MakeDegeneratedEdge`（`BRepPrim_OneAxis.cxx:934/1212/1264`）与 `BRep_Builder::Degenerated`（`BRep_Builder.cxx:1073-1085`，OCCT 会**丢掉 3D 曲线**，本 port 的 `EdgeGeom::curve` 非空故用占位）。
  - `GeometryRegistry::set_edge_range(s, first, last)`（对齐 `BRep_Builder::Range`，供 `SetParameters(ETOP/EBOTTOM, …, 0., myAngle)`，`BRepPrim_OneAxis.cxx:407/418`）。
- **球面完整 wire 实测（`LateralWire` 四条边：Top 退化 + End 经线(reversed) + Bottom 退化 + Start 经线；pcurve 全部为 `gp_Lin2d`：极点 `(0,±π/2)+X` 范围 `[0,2π]`、经线 `(0/2π,−2π)+Y` 范围 `[3π/2,5π/2]`）**：
  - `wires=1`、`uv_bounds=(0,2π,−π/2,π/2)`、`classify=Sphere`、`analytic_surface_area=4π` ✓
  - **`incremental_mesh_to_shape_mesh` 首次成功**：`168 verts / 306 tris`（未给极点边设范围时是 `4096/7686`；此前**完全失败** `face 0 has no boundary UV points`）⇒ **A13/A18 的前置在网格侧已打通**。
  - **仍未通过**：`brep_gprop_full::surface_properties/volume_properties` 返回 **0**（`--lib` 会变 1291/3，新增 `sphere_surface_volume`、`adaptive_box_sphere`）。路径已缩小到 `FaceGauss`：`compute_face` 判 `rect_domain`（四条 pcurve 都是直线）后走 `compute_rect → compute_natural`，其入口守卫 `u2 > u1 && v2 > v1` 不成立即返回 0 ⇒ 问题在 `FaceGauss::new` 自己的 `uv_bounds(surface, &arcs)`/`build_arc(&e,&map)`（`brep_gprop_full/integration.rs:36-60`），与 `brep_tools::uv_bounds`（实测正确）不是同一条实现。
- **下一步（精确）**：读 `brep_gprop_full/integration.rs` 的 `build_arc` 与该模块的 `uv_bounds(surface,&arcs)`，定位为何闭合球面的四条直线 pcurve 得不到非退化 UV 盒；修好后球面 wire 即可重放（随后 A13/A18 → T-59）。球面基元代码已按纪律回退，避免半成品。

**T-68 第三轮（2026-09-20）：球面前置**打通**——`FaceGauss` 改为 pcurve 优先 + 球面 wire 落地，Gauss 面积/体积恢复 4π / 4⁄3π**

- **新发现（A28，已修）**：`brep_gprop_full::FaceGauss` 的边界弧原先是"3D 曲线 + UV 反演"（`BoundaryArc::value/d12d` 走 `UVMap::map_point`）。这**无法表达接缝两侧**（u=0 vs u=2π 在 3D 上同一点）、也建不出退化边的弧，于是闭合面（球/环）的 UV 盒退化成一点 ⇒ `compute_rect` 守卫 `u2 > u1 && v2 > v1` 不成立 ⇒ 面积/体积**静默为 0**。修法：`BoundaryArc` 改 `ArcGeom::{Pcurve, Curve3d}`，`build_arc(e, face, map)` **优先取已存 pcurve**（连 `pcurve_range`），无 pcurve 才回退 3D+反演；新增 `classify_arc_kind2d`。这同时对齐了 OCCT `BRepGProp_Face` 走 `BRep_Tool::CurveOnSurface` 的事实。
- **球面 wire 落地**（T-68 步 1+2；步 1 的 `make_degenerated_edge`/`set_edge_range` 见上一轮提交）：四条边 `Top(退化,fwd) → End(经线,rev) → Bottom(退化,rev) → Start(经线,fwd)`，pcurve 全为 `gp_Lin2d`（极点 `(0,±π/2)+X` 范围 `[0,2π]`；经线 `(0/2π,−2π)+Y` 范围 `[3π/2,5π/2]`），极点边范围按 `SetParameters(...,0,myAngle)` 设为 `[0,2π]`。
- **实测（临时探针，已删）**：`wires=1`、`uv_bounds=(0,2π,−π/2,π/2)`、**`brep_gprop_full` 面积 = 12.566370614359137（4π，Δ=3.5e-14）**、**体积 = 4.188790204786379（4⁄3π，Δ=1.2e-14）**、`analytic_surface_area=4π`、忠实网格 168/306、栅格网格 4096/7686。
- **过期期望订正**：`fclass2d::tests::infinite_point_closed_periodic_face_is_in` 的注释自承前提是"球面无边界 wire"（port artifact）；现在球面有了真实边界 ⇒ 无限点为 `Out`，测试改名并注明 `LateralWire` 出处（同 T-06 先例：改的是编码 artifact 的期望，不是新写测试）。
- **验证**：`occt-topo --lib` **1293/1**（唯一红仍 T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3、`occt-geom` 151、`occt-geom2d` 72、`occt-core` 290 —— **逐项与基线一致**。
- **T-68 剩余**：环面（`BRepPrim_Torus::SetMeridian` 的闭合经线分支：`MeridianClosed()` ⇒ 单边双 pcurve，`cxx:389-396`）与带孔面的内环（`ring area 4 vs 3.8037`）⇒ 完成后重放 A13/A18，再重放 T-59。

**T-68 第四轮（2026-09-20）：A19 独立落地；A13/A18 重放的分诊从 9 个失败降到 5 个**

- **A19（T-55）完成**：只删"失败率阈值 + 两级静默换网格器"（`perform` 直接传播错误、失败面只置 `MeshStatus::FAILURE`；删 `build_shape_mesh_wireframe`、`wireframe_face_triangulation`、`WIREFRAME_FALLBACK_RATIO_MAX`、`discretize_face` 的两处回退），**保留** `wireframe` 的 UV 栅格给其余消费者 ⇒ `--lib` **1293/1**、四道门禁与基线一致。这样仓库明令禁止的"OCCT 里不存在的失败率规则"已彻底移除，且不牵连 A13/A18。
- **A13/A18 重放（分诊）**：把 `face_to_triangles` 也改成委托忠实管线后，`--lib` 由 1293/1 变 **1288/6**——比第三轮（round 20 的 1291/10）**少 4 个失败**：`viz_scene`×3 与 `vrml` 的球面用例**已随 T-68 步 1+2 通过** ✓。剩余 5 个（已逐一定位）：
  1. `bop_curved::sphere_inside_box_common_volume`：共同体积 4.1288 vs 4.1888（**低 1.43%**）——忠实 Delaunay 在该偏转下比旧栅格略粗，旧容差是按栅格标定的（属"期望需按新网格重标定"，非错误）；
  2. `brep_pattern::pattern_volume_sums_copies`：两球体积 **0** ⇒ 模式副本的分离面**无 pcurve/wire**；
  3. `draw::sphere_fillet_offset`：offset 体积 **0** ⇒ 同上（分离/派生形状缺边界）；
  4. `rwmesh::roundtrip_shape_mesh`：**panic** `make_edge_segment_with_vertices: p1 and p2 must be distinct` ⇒ 网格→BRep 路径遇到**零长边**（与 T-32 的 `mesh_to_brep` 同族，需按 `BRepBuilderAPI_MakeShapeOnMesh` 跳过退化段）；
  5. `wireframe::face_with_hole_triangulates_ring_area`：ring area **4 vs 3.8037** ⇒ 带孔面的**内环未生效**（`model_builder` 内 wire 处理）。
- **结论**：A13/A18 的剩余阻塞点已从"球面完全无边界"收敛为 **分离/派生形状的边界传播**、**mesh→BRep 零长边**、**带孔面内环**三件事。A13/A18 暂缓（重放会红 5 个），A19 已单独落地。

**T-68 第五轮（2026-09-20）：T0M 的 169 面缺口被测出——A19 只能先删"被禁的失败率规则"**

- **重要测量**：按 A19 把"失败率阈值 + 两级回退"全删后，`step_obj_parity` **14/14 → 13/14**：`data/occ/T0M.stp` 的 `bbox min[2] ours=-424.412306008 occ=-424.741875692`（**短 0.33**）。临时探针（已删）定位：`IncrementalMesh::perform` 返回 Ok，但 **`face_stats` 只有 1603 面 / 形状 1772 面** ⇒ **169 个面（9.5%）被判 FAILURE 且无三角化**，此前正是被那条"逐面 UV 栅格回退"悄悄补齐的（169/1772 = 9.5% **低于** 10% 阈值，所以阈值本身没触发——被禁的规则与被依赖的回退是两件事）。
- **本轮落地（A19 的实质部分）**：删除 `WIREFRAME_FALLBACK_RATIO_MAX = 0.10` 及其"超 10% 即整形状改用 UV 栅格"分支——这是仓库明令禁止、OCCT 中不存在的失败率规则；**逐面回退暂留**并在原地写明 `UNPORTED` 与实测依据（169/1772），待 T-68 补完该缺口后删除。
- **验证**：`occt-topo --lib` **1293/1**、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` 2/3 —— 全部与基线一致。
- **下一步（精确）**：定位 T0M 那 169 个面为何判 FAILURE（提示：单面探针 `IncrementalMesh::discretize_face` 已不再回退，可直接逐面统计并打印首个失败面的 surface/错误），这同时是 A13/A18 与 T-55 剩余部分的共同前置。

**T-68 第六轮（2026-09-20）：T0M 缺口的根因锁定为"带孔面（2 wire）+ 非平面"的 Delaunay 空结果**

- **方法**：在 `add_wire` 的 4 个 early-return 与 `triangulate_model_faces` 的 Delaunay 失败/空结果处临时插 `eprintln!`（6 处，全部已删并被 `git status` 复核为 0 残留），跑 T0M 全形状 `perform` 后按原因聚合。
- **结果（决定性）**：
  - `add_wire` **从未失败**（0 条记录）——面都有 wire（1772 面实测：1 wire 1463、2 wires 294、3 wires 15），原先怀疑的"缺 wire/坏 pcurve"不成立；
  - 带回退时 `face_stats = 1769 / 1772`（3 面为 `REUSED`，正常跳过），**168 面走回退**：其中 **166 面是 `map_triangulation` 返回空（Delaunay 成功但产出 0 三角形）**，另 **2 面**（1691、1758）报 `DelaunayNodeInsertionMeshAlgo::perform: face N has an invalid discrete range`；
  - 逐面 profile：**166 面的特征完全一致 —— `wires=2`（外环 + 内环 = 带孔面）且 `surface=Other`（非平面）**。
- **结论**：T0M 的 bbox 缺口（删回退后短 0.33）**不是缺边界**，而是**带孔的非平面面**在 Delaunay 约束路径上产出空网格（与 `wireframe::tests::face_with_hole_triangulates_ring_area` 的"内环未生效"同源）。⇒ 修法应落在**内环作为约束的处理**（`model_builder` 把内 wire 交给面模型 + 约束 Delaunay 的孔洞多边形处理，参考 `BRepMesh_Delaun` 的 `meshLeftPolygonOf`/约束边路径），以及 2 例 `invalid discrete range` 的离散域校验。
- **验证**：本轮**未改行为代码**（仅临时探针，已清零）；`occt-topo --lib` **1293/1**，工作树干净。

**T-68 第七轮（2026-09-20）：缺口继续下钻到 `Delaun::compute` 产不出三角形**

- **证据链（临时探针，已清零）**：在 `node_insertion::finish_mesh` 加打印（记录 `nodes/links/du/dv/frontier/base_domain_tris`，全 1769 面），对 T0M 跑全形状后聚合：
  - 166 个失败面：`base_domain_tris=0` **且 `frontier=0`**，但已注册 **50–116 节点 / 52+ 链接**；
  - 1603 个成功面：`base_domain_tris>0`、`frontier>0`。
- **推论**：`frontier=0` 是**结果**而非原因——`Delaun` 构造期（`new_with_data_cells` → `perform` → `super_mesh` + `compute`）**一个三角形都没产出**，于是所有链接都"无相连三角形"，被随后的 `erase_free_links()` 全部清掉。已核对 `erase_free_links` 与 OCCT `BRepMesh_MeshTool::EraseFreeLinks`（`cxx:204-219`）逐行一致 ⇒ **不是这里的偏差**，下一站是 `Delaun::compute`/`create_triangles`（`delaun/triangulation.rs:305-316`、`:366+`）在"带孔面"配置下的行为。
- **排除项**：UV 跨度不是决定因素——失败面 `du ∈ [4.02, 2332.6]`（均值 85），成功面 `du ∈ [0.005, 610.6]`（均值 15.2），区间**重叠**；决定因素仍是 **`wires=2`（外环+内环）且 `surface=Other`**（166/166）。
- **下一步**：对单个失败面插桩 `compute`（`loop_edges` 初值、`create_triangles(first)` 是否产出、`create_triangles_on_new_vertices` 的插入数），并与 OCCT `BRepMesh_Delaun::compute`/`createTriangles` 在"外环+内环"配置下对照——这是 A13/A18、T-55 剩余回退、T-59 三者的最后一道共同前置。

**T-68 第八轮（2026-09-20）：缺口锁定到 `process_constraints()`——它把 100+ 三角全删成 0（= A25 的自创删除补丁）**

- **证据（临时插桩，已清零）**：`Delaun::compute` 与 `create_triangles_on_new_vertices` 分阶段计数，对 T0M 聚合：
  - `sup=[53,54,55]`（超三角形正常）、`create_triangles(first)` 后 **`after_first=3`**（超三角形被拆成 3 个）✓；
  - 顶点循环结束后 **`after_loop=101..159`（即 100+ 个三角形，正常）**；
  - **`process_constraints()` 之后 `after_pc=0`** —— **166/166 个失败面全部如此**，三角化被整体销毁。
- **定位**：`process_constraints()`（`BRepMesh_Delaun.cxx:703` 尾调；body `insertInternalEdges(); frontierAdjust()`）→ `frontierAdjust` → `delaun/polygon_meshing.rs::decompose_simple_polygon`（`:228-...`）。该函数正是 **审查 A25 记录的自创实现**：在 `:346-375` 它把"ear 复用链接"的**相邻三角形删掉**再 `add_triangle_by_info`，而 OCCT `BRepMesh_Delaun.cxx:2259-2274` 是**直接 `AddLink` + `addTriangle`、无任何删除**（OCCT 在退化到第三条连接时由 `BRepMesh_PairOfIndex.hxx:41` 抛 `Standard_OutOfRange`，被 `BRepMesh_BaseMeshAlgo.cxx:62` 空 catch 吞掉，该面干脆没有三角网）。
- **推论**：A25 的"删邻三角形"补丁在"外环+内环"（`wires=2`）配置下会把整个 2D 网格删空 ⇒ 这就是 T0M 166 面缺口的直接原因，也解释了为何 `frame_with_hole` 的 ring area 变成整盘（4.0）。**A25 与 T-68 的这个阻塞点是同一处代码。**
- **修法（下一轮，对着 `.cxx` 做）**：按 `BRepMesh_Delaun.cxx`（`frontierAdjust` / 多边形分解 / `meshLeftPolygonOf`）**去掉删除逻辑**，改为 OCCT 的 `AddLink`+`addTriangle` + 失败即抛出（端口用 `Err`/跳过该多边形，等价于 OCCT 的空 catch ⇒ 该面无网格）；随后重放 A13/A18（`face_to_triangles` → 忠实管线）与 T-55 的回退删除，并复验 T0M 的 `step_obj_parity`（届时若 OCCT 本身也不网格化那些面，bbox 应与参考一致而非变短——这正是本轮结论要验证的下一步）。
- **验证**：本轮**未改行为代码**；`occt-topo --lib` **1293/1**，工作树干净、插桩 grep 复核 0 残留。

**T-69 第 1 轮诊断（2026-09-20 第 41 轮 goal round）：166 面缺口的成因已定位到 `FixLacking` 追加重复闭合边 + 单闭合边 wire 的 UV 回绕链；参考网格证明 OCCT 同样不网格化这些面 ⇒ 真正的偏离是端口保留的逐面 UV 栅格回退**

- **量化（临时探针，已删）**：在 `incremental_mesh/discret_root.rs::triangulate_model_faces` 的 `needs_fallback` 处插桩，T0M（`data/occ/T0M.stp`，1772 faces，端口模型 1769 面）上 **166 面**需要回退，其中 **162 面**是"`perform` 正常返回但 `map_triangulation` 拿到 0 三角"（`tris=0`；节点数 `nodes=74` 占 134 面，其余 110/111/112/75/72/68/55）、**2 面** `MeshStatus::FAILURE`、**2 面** `perform` panic。整模仍为 **46516 v / 46962 f**（与批 12 后逐位一致 ⇒ 回退面在计数上"看起来正常"）。
- **结构（探针：`node_insertion.rs::collect_boundary_uv` + `step/read_topology.rs::{resolve_loop,resolve_face}`）**：166 面全是**周期面**（首个样本 `surface=Cylinder`，`wires=2`），其 wire 在 STEP 里是**单条闭合边**的 `EDGE_LOOP`——T0M 共 **596** 个 `EDGE_LOOP items=1`（`resolve_loop` 返回 1 条边，未复制）。这些面的前沿链：
  - **有 `FixLacking` 时**：`shhealing::check_pcurves_and_shift`（= `ShapeFix_Wire::Perform` 的移植）末尾的 `fix_lacking_all` 给 wire **追加一条重复的闭合边**（探针 `[dbg-t69-step] wire edge count before=1 after_proj=1 after_check=2` 命中 **362** 次；另有 `3→4`×49、`2→3`×16、`4→5`×11 等）。两条链几乎重合而**参数化不同**（`e3d[0,2π]` 对 `e3d[π/2,5π/2]`），正是 round 9 观测到的"两端各差 ~0.05、v 完全相同"的两条重叠共线链接 ⇒ `Glued` 批量删段 ⇒ 拼接链断裂 ⇒ 清空。
  - **跳过 `FixLacking` 时**（临时 `DSH_T69_NO_LACKING` 开关）：回退面 **166 → 86**，剩下的链仍退化——单闭合边 wire 走 `BRepMesh_ModelHealer.cxx:452-457`（端口 `model_healer.rs:236-241`，**忠实**）的分支，把 pcurve **首点 UV 覆盖为末点 UV**（`u=0 → u=2π`），得到 `[2π, 0.1745, …, 6.1087, 2π]`——一条跨越整个 2π 的**回绕链**，整条链在 UV 上共线（面积 0）⇒ `decompose_simple_polygon` **正确地**判"无耳"清空（与 round 9 的 `skip_prec/skip_neg` 分支口径一致）。
- **决定性反证（部分成立，已在同日补测收敛口径）**：**OCCT 的参考网格对"单闭合边 loop 的管面"同样没有网格**。`data/occ-ref/T0M.obj` 逐面写顶点、**不去重**（51160 条 `v` 中 **20958 条完全重复**）；失败管（半径 9.75、圆心 `(0,−100.454,−304.817)`）的顶圆在这份参考里**只出现一次**（按半径 9.75±0.01 数 = **36** 点，若侧壁也被网格化应出现 72 点），且只被**一个**面引用，而该面所有顶点两两距离 ≤ 20（即圆盘盖，直径 19.5）——**该管面的侧壁没有三角形**。⇒ 这一面类（T0M 里 596 个单闭合边 loop）不是端口独有的失败。
- **同日补测（推翻"可以删回退"的推论，故未采纳）**：临时删除 `IncrementalMesh::wireframe_face_triangulation` 及其调用（faithful 行为 = 失败面只置 `IMeshData_Failure`）后，`step_obj_parity` 变 **13/14**：`occ/T0M.stp: bbox min[2] ours=-424.412306008 occ=-424.741875692`（端口短 0.33）⇒ **参考网格确实覆盖了忠实路径缺失的区域**。查参考里 `z<-424.7` 的顶点只有 2 个（**完全重合**：v962/v965 = `(−0, −11.823, −424.742)`，分属两个面、被 6 个三角形使用）⇒ 至少还有**一类**端口失败面是 **OCCT 网格化了的**（该点在一次退化特征处，非 face 0 那条管）。⇒ 结论收敛为：**端口的失败面集合是 OCCT 失败面集合的超集**（管面类两边都不网格化；另有至少一类"OCCT 会网格化"的面端口失败），逐面回退因此在当前基线上仍"承重"，**不得直接删除**。按门禁 4（不得劣于基线）已 `git checkout` 回退该改动，`step_obj_parity` 复测 **14/14**。
- **下一步（T-69 收口，写进任务卡）**：① 用探针定位"覆盖 `(−0,−11.823,−424.742)` 的那张端口失败面"（在 `triangulate_model_faces` 里按 `needs_fallback` 面逐个对 bbox/最近点过滤），dump 其 `wires/edges/pcurve/same_param/链`并与同类的**成功**面对比，找出端口在该面类上少走的那段控制流（这是 OCCT 网格化、端口没网格化的真缺陷）；② 修好后重跑门禁并**按面类**记录 T0M 计数；③ 再评估删除逐面回退（管面类届时可交回 `IMeshData_Failure`，与 OCCT 一致）。
- **口径订正**：round 9 的"来源在 Delaunay 之外"成立，但"前沿链重叠共线"的**直接来源已定位**为上面两条（`FixLacking` 的重复边 / 单边 wire 的回绕链），而不是网格化前的"重复插点"；同时"这 166 面 OCCT 也没网格化"只对**单闭合边管面类**成立（参考实测该管侧壁无三角），而"删掉逐面回退"的推论被 `step_obj_parity` 13/14（T0M `min[2]` 短 0.33、参考另覆盖 `(−0,−11.823,−424.742)`）**推翻**，故回退未采纳。
- **验证**：本轮**未改行为代码**（探针与插桩全部 `git checkout` 还原 + `git grep dbg-t69` 0 残留，工作树干净；唯一一次行为改动——删逐面回退——按失配即停回退）；`occt-topo --all-targets` 编译 exit 0；`occt-topo --lib` **1293/1**；`step_obj_parity` **14/14**、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05，与基线一致）。
- **下一步（T-69 收口，写进任务卡）**：① 定位"覆盖 `(−0,−11.823,−424.742)` 的端口失败面"，与同类**成功**面对比控制流（这是 OCCT 网格化、端口没网格化的真缺陷）；② 修好后按面类记录 T0M 计数并重跑门禁；③ 再评估删除逐面回退。

**T-69 第 2 轮（2026-09-20 第 42 轮 goal round）：定位到"OCCT 会网格化、端口失败"的那类面（Torus，face 20），差异收敛为"同面两条 wire 的 pcurve 落在相差一个周期的 u 窗口"**

- **归属face（临时探针，已删）**：在 `IncrementalMesh::perform` 的 face 循环里对每个面的网格顶点求到参考 z 最小点 `(−0, −11.823, −424.742)` 的距离 ⇒ 全局只有 **face 20** 命中（最近顶点 **0.000322**，`surface=Torus`，nv=12/nt=12 = 回退网格）；且它正是"忠实 Delaunay 返回 0 三角"的面（`face=20 mapping_failed tris=0 nodes=74`）。其余已网格化面距该点都 >1.0。
- **成功/失败同类的差异（探针：`node_insertion.rs::collect_boundary_uv` 的逐 slot dump + `perform` 里 splitter 的 `range_u/range_v`）**：
  - **成功的 Torus face 8**（nv=74/nt=72，忠实路径产出）：wire 0（v=2.6505）与 wire 1（v=3.6327）的 pcurve **同在 u ∈ [1.5708, 7.8540]**；splitter `range_u=(1.5708, 7.8540)`，两条链都在范围内 ⇒ 干净条带 ⇒ 72 三角。
  - **失败的 Torus face 20**：wire 0（v=0.05088）的 pcurve 在 **u ∈ [1.5708, 7.8540]**，wire 1（v=π/2）的在 **u ∈ [−4.7124, 1.5708]** —— **正好相差一个周期（2π）**；`update_range` 的周期钳制（`BRepMesh_DefaultRangeSplitter::updateRange` `cxx:202-235`；端口 `range_splitter/param_set.rs:214-235`，已忠实）把 `range_u` 钳成 `(−4.7124, 1.5708)` ⇒ **wire 0 的链恰好落在范围外一个周期**（缩放后 u ∈ [1,2]）。
  - 两面的**单 wire 链形状完全相同**（每 wire 2 条边、都是 `FixLacking` 复制的闭合圆 ⇒ 出-回链，37+37=74 节点），所以**目前观测到的唯一差异就是"同面两条 wire 的 u 窗口是否对齐"**。四个边的 `same_param/same_range` 全为 true，且 `pc.d0(t)=t`（每边自身一致），说明两个 u 窗口来自 STEP 文件各自的圆参数化，不是边内参数化错。
- **参考侧证据**：`data/occ-ref/T0M.obj` 里该顶点出现 **两次且完全重合**（v962/v965）⇒ 该顶点被**两个面**各自写了一次（逐面写顶点不去重）⇒ 至少与它相邻的某个面在 OCCT 里**是**被网格化的；端口忠实路径在这些面上给 0 三角，故"OCCT 会网格化、端口失败"的面类**确实存在**，且与本轮的 u 窗口差异一致。
**批 98 工作令（第 107 轮勘察落定；未改代码）—— `Geom_BSplineCurve::Segment` + `SetOrigin`**

> 本轮回合预算用于五 crate 复检（全部 exit 0 / 0 error）与把下批的实现清单落到行号，避免下轮返工。

- **目标**：移植 `Geom_BSplineCurve::Segment(U1,U2,theTolerance)`（`Geom_BSplineCurve.cxx:527-715`）与其周期分支所需的 `SetOrigin(Index)`（`:819-880`）；**三个消费方**同批接线：① IGES 窄区间 B 样条回落（`iges.rs::emit_bspline_curve` 的范围检查，OCCT `GeomToIGES_GeomCurve.cxx:320-356`）；② `GeomConvert::CurveToBSplineCurve` 的 trimmed-Bezier/BSpline 臂（`GeomConvert.cxx:300-339`）；③ STEP **写**侧 `GeomToStep_MakeCurve.cxx:71-82`（trimmed 基曲线先 `Segment` 再写，正是 `write_context.rs` 里已登记的"写未裁剪基曲线"偏差）。
- **`Segment` 的控制流（已逐句核对）**：`U2<U1` → `DomainError`；周期前置（`Period=Last-First`、`DU=U2-U1`、`DU-Period > PConfusion` → `DomainError`、`DU` 钳到 `Period`、记 `aDDU=DU`）；两次 `BSplCLib::LocateParameter(deg, myKnots, myMults, U, periodic, myKnots.Lower(), myKnots.Upper(), index=0, NewU)`——注意 OCCT 这里传的是**distinct 结点数组的 Lower/Upper**，端口应直接调 `bspl::locate::locate_parameter_range(knots, u, periodic, first, last, knots[first], knots[last])`，**不能**用便利包装 `locate_parameter`（它会按 `First/LastUKnotIndex` 自算 first/last）；随后 `Knots={min,max}`、`Mults={deg,deg}`、`Eps=max(Epsilon(AbsUMax), tol)`、`InsertKnots(Knots,Mults,Eps)`；周期分支：再 `LocateParameter(U1)` → `SetOrigin(index)` → `SetNotPeriodic()`（批 95 已有）→ `NewU2 = NewU1 + DU`；再按 `LocateParameter(NewU1/NewU2, FromU1=Lower, ToU2=Upper)` 定 `index1/index2`（含 `|myKnots(index+1)-U| <= Eps` 右移规则与 `index2==index1 ⇒ ++`）；新结点/重数取 `index1..=index2` 切片、两端重数置 `deg+1`、周期时整体减 `DU = NewU1 - U1`；极/权取 `PoleIndex(deg,index1/2,periodic,mults)`（`pindex1++`、`pindex2=min(pindex2+1,len)`）；周期尾部把首结点置 `U1`、当 `aNu2 < U2` 时末结点置 `U1 + aDDU`；最后重建平结串。
- **落地前置（必须先对拍）**：`Geom_BSplineCurve::InsertKnots(Knots, Mults, Eps)`——OCCT 是"插到指定重数并在 `Eps` 内**合并**已有结点"，端口 `GeomBSplineCurve::insert_knot(u, mult)` 是无条件连插 `mult` 次（= OCCT `InsertKnot(..., Add=false)` 的语义，而 OCCT 默认 `Add=true` 是"把重数升到 M"）；**旁支（第 108 轮 grep 实证）**：端口该方法**当前零调用点**（`git grep '\.insert_knot\('` 无命中；`convert_bspl.rs` 的曲面插入直接用 `boehm::boehm_insert`）⇒ 批 98 移植 `InsertKnots` 时按 OCCT 语义实现，并顺带处理这个既有死方法（改造成忠实实现或标 UNPORTED，勿留两套语义）。`SetOrigin(Index)` 需按 `:819-880` 移植（周期结点/极点数组的循环移位 + `LocateParameter` 归一）。**第 109 轮补齐行号**：`InsertKnot(U,M,tol,Add)` `Geom_BSplineCurve.cxx:337-350`、`InsertKnots(Knots,Mults,Eps)` `:351-…`、`IncreaseMultiplicity(Index,M)` `:302-312`（批量臂 `:313-336`）、`SetOrigin(int Index)` `:819-…` 与 `SetOrigin(U,Tol)` `:913-…`（后者是前者的 `LocateParameter` 包装）。
- **验证口径**：仅五 crate 编译（门禁等指令）；`iges_check` 的实体统计与 STEP 门禁写入差分在门禁波次统一评估。
- **✅ 已执行 = 第 110 轮批 98**（见下条日志）：本工作令的四件（`InsertKnots` 族 + `InsertKnot`/`IncreaseMultiplicity` + `SetOrigin` 两式 + `Segment`）与三个消费方中的两个半（IGES 窄区间、`CurveToBSplineCurve` 的 trimmed-**BSpline** 臂、STEP 写侧 trimmed-**BSpline** 基曲线）均已落地；**trimmed-Bezier 臂**因缺 `Geom_BezierCurve::Segment` 另立 **R2-23**。

**批 100（T-32 输入侧：移植 `BRepLib_MakeWire` + 闭合实体定向；全门禁复跑）—— 2026-09-21 第 114 轮**

> 执行 §3.2 序 1。忠实移植原则：先读 `.cxx` 控制流 → 逐句翻译 → 编译 → 跑真实几何对照基线 → 失配即停报告。

- **新增 `crates/occt-topo/src/brep_lib_make_wire.rs`（`MakeWire` = `BRepLib_MakeWire`）**：
  - `Add(Edge, IsCheckGeometryProximity)`（`BRepLib_MakeWire.cxx:123-453`）**整段移植**：首边建 wire 并播种 `myVertices`（`:135-149`）、`init = myShape.Closed()`（`:153`）、`E.Oriented(TopAbs_FORWARD)` 子顶点迭代（`:154-171`）、`VF/VL` 空值 ⇒ `NonManifoldWire`（`:163-169`）、逐顶点 `Contains`（`:177`）与四路定向判定（`:184-220`）、proximity 分支（`:223-285`，含 `l < tolE || l < tolW`）、`DisconnectedWire` 早退（`:288-293`）、copy-edge 支（`:304-368`，`EmptyCopied` ≡ 用同曲线/区间重建 + `UpdateVertex` 加权点与 `max(tolW,tolE)` 容差）、**定向决策** `((forward==reverse)&&E.Reversed) || (reverse && !forward)`（`:370-379`）、`VF/VL` 更新与 `init` 特例（`:386-444`）、闭合判定 `myShape.Closed(true)`（`:445-449`）、`EmptyWire/WireDone` 语义（`:35-38`、`:451-452`）。
  - 配套 OCCT 语义：`TopoDS_Iterator(shape, CumOri=true)`（存储序 + 存储朝向）、`TopExp::Vertices(E,V1,V2,CumOri=true)`（`:393`，朝向 REVERSED 时 V1/V2 互换）、`TopTools_ShapeMapHasher` = `TopoDS_Shape::IsSame`（`TopoDS_Shape.hxx`：同 TShape + 同 location；两个 null 句柄相等 ⇒ `is_same(None,None)=true`）。**登记边界**：`IsSame` 的 location 部分按本仓惯例退化为 identity（`bop_occt_util.rs:31`），本路径 location 恒 identity。
  - **两条控制流单测**（+2 测试）：`reorients_shared_edges_into_a_chain`（共享顶点、外来方向的边被定向成头尾相接的闭合环）、`disconnected_edge_is_reported`。
- **接线**：`mesh_to_brep::triangulation_to_brep` 的三角形 wire 改由 `MakeWire` 建（`build_wire`；`MakeWire` 未 done 时保底回落到 `BRep_Builder` 级 append 并就地注明——网格三角形的三边共享顶点，实际不可达）；**闭合实体补 `BRepTools::OrientClosedSolid`**（`brep_class3d::orient_closed_solid`：无穷点分类为 IN 时反转）。
- **临时探针（已删）实测（`mesh_cylinder(1,3,24)` 夹具）**：输入网格 **signed volume = −9.3175 ⇒ 反向缠绕**（OCCT 的 BRepMesh 只会给外向三角形）；批 100 前该夹具转出的 BRep 是"材料在外"的反实体；批 100 后 `faces 96 / edges 144`、**「面法向 · wire 绕向」96 正 / 0 负**（批 100 前 round 9 记录 47/96 为负）。
- **真实几何对照（T-01 `groove_cuts_cylinder`）**：`8.903438 → **8.408733**`（期望 7.708990，容差 0.5；等价于 removed 0.414 → **0.909**，应 1.6085）⇒ **仍红，但已明显靠近**。
- **门禁（复跑 §6 全景）**：lib **1285 / 2**（= 1286 基线 − 3 条批 99 删除用例 + 2 条本批新增；两条红仍是 T-01 与 T-86）、`bop_builder2_boss` 1/2、`phase10` 7/8、`phase19` 3/5、`phase20/3/4/5/6/7/8/9` 全绿、`step_geometry_parity` 2/3、`step_obj_area` 11/11、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、doc-tests 2/1i；`export_data_obj` **16/16 且计数与批 98/99 逐位相同**（ATU01038 17745/22119 …）⇒ **零回归**；round 8 点名的两个易回归用例（`groove_negative_volume_delta`、`neck_*`/`rib_*`/`boss_*`）全绿。
- **结论与下一步（失配即停）**：T-32 的**输入侧**（①链向 ②面朝向 ③闭合实体定向）已落地并验证；T-01 剩余 0.7 的偏差**不在 mesh→BRep**，而在忠实 BOP 链内部（`bop_builder::boolean` → `bop_builder2::builder_bop_with_fuzzy` → 面级切割/装配，即 round 2–6 诊断过的 `rebuild_split_areas`/`perform_loops`）⇒ 下一批 = **T-01 round-15 复测**（见下条）。

**T-01 round-14 复测（第 114 轮，输入侧修好后；临时插桩已还原）**

- **配方**：在 `bop_split_faces_occt.rs::build_split_faces_occt` 的 `rebuild_split_areas` 调用处临时统计 `tasks / 输入边数 / 每个面的 areas 块数直方图`，跑 `cargo test --lib brepfeat::tests::groove_cuts_cylinder -- --nocapture`；跑完已 `git checkout` 还原（`git grep dbg-t01` 复核 0 残留）。
- **实测**：`tasks=112 edges_in=1232 **area_hist={1: 71, 2: 16, 3: 25}**`；1 块面的样本 `face 13 le=15 / face 21 le=13 / face 29 le=13 / face 37 le=15 / face 45 le=13 / face 53 le=13`。
- **与 round 2 基线（批 100 前）对比**：round 2 = `[(0,40),(1,41),(2,30),(3,1)]` ⇒ **"0 块"面 40 → 0**（输入 law 修好后不再有整面丢弃），多块面 31 → **41**，但 **1 块面 41 → 71**。
- **判读**：输入侧修复**消除了 round 2–3 的 40 个"0 块"坏面**（有效），但仍有 **71 个面只出 1 块**，其中样本面 `le=13/15` ⇒ **拿到 13–15 条候选边却只成一个 area**，与 round 3–6 定位的"`avoid=3`（1 条原边界碎片 + 2 条截面边未链进环）"同类。⇒ 剩余缺口仍在 **`perform_loops` 的 WireSplitter 成环**（`bop_builder` 侧的 `builder_face_occt.rs` / `wire_splitter_block.rs`），而不在 PaveFiller/FF 求交。
- **下一批（round-15）配方**：按 round 3 的方式对 `perform_loops` 插桩统计 `(in_edges, avoid, avoid_bnd, loops, areas)` 分组，对照 `BOPAlgo_BuilderFace::PerformLoops` / `BOPAlgo_WireSplitter`（`_1.cxx:112-354`、`:358-617`）找那一类面为何仍不闭环；**禁止**为某个用例调参/加谓词。

**批 99（R2-17：删除 OCCT 不存在的 circle/circle 快路径；全门禁复跑）—— 2026-09-21 第 112 轮**

> §6.2 playbook 步骤 6。判据 = "`perform` 的分派与 OCCT `IntTools_EdgeEdge::Perform` 完全一致，且门禁不劣化"。

- **OCCT 对照（已读 `.cxx`）**：`IntTools_EdgeEdge::Perform`（`IntTools_EdgeEdge.cxx:185-243`）= `CheckData` → `Prepare` → **line/line 特判 `ComputeLineLine` 并 return**（`:198-202`）→ 快速重合检查 `IsCoincident` + `AddSolution`（`:204-215`）→ `line + 解析曲线(≤Parabola)` 的最小距离早退（`:217-233`）→ `FindSolutions` + `MergeSolutions`（`:235-242`）。**全篇没有 circle/circle 分支**；圆–圆与其余组合一样走 `FindSolutions`。
- **改动**：① `edge_edge/edge_edge.rs::perform` 的 `match (ctype1, ctype2)` 删除，非 line/line 一律 `self.find_solutions()`；② 删 `compute_circle_circle`（`is_coincident` + 采样 `intersect_edges` 的替身分派）；③ 删 `solvers.rs::compute_circle_circle_full`（端口自创的共面圆解析解：同心/相离/相切/两交点的 radical-line 构造）及仅它使用的 `fallback_general`、`point_hit_on_both`；④ 删 `mod.rs` 的三条直调用例 `circle_circle_full_{tangent_single_point,intersecting_two_points,separated_empty}`（−3 条测试）；⑤ 模块头/文件头注释改为"无 circle 分支"并保留一处 UNPORTED 登记（见下）。
- **保留**：`compute_line_line_full`（`ComputeLineLine` 的忠实件，仍被 `tests_full` 用）与 `intersect_edges`（`edge_edge.rs:336` 的 line/line 交叉分支与测试仍用）。
- **验证**：
  - `cargo check --all-targets` **exit 0 / 0 error**（无新增告警，除 `pub use solvers::*` 这条**既有**告警）；
  - `cargo test --lib edge_edge::` → **15/15 全绿**：`circle_circle_two_hits`（忠实 `FindSolutions` 给出 2 个交点，x=0.5、|y|=√0.75）、`separated_circles_empty`、`matches_inttools_circle_circle`、`line_bspline_crossing_two_hits` 全部保持；
  - 全量门禁（playbook 步骤 3 复跑）：lib **1283 / 2**（＝原 1286/2 减去删掉的 3 条用例，两条红仍是 T-01 与 T-86）、boss 1/2、phase10 7/8、phase19 3/5、phase20·3·4·5·6·7·8·9 全绿、`step_geometry_parity` 2/3、`step_obj_area` 11/11、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、doc-tests 2/1i ⇒ **与 §2 基线逐项相同**；  - `export_data_obj` **16/16**，计数与 R2-17 之前**逐位相同**（ATU01038 17745/22119、Cube 24/12、Sphere 642/1244、Torus 1369/2592、Shape-1 3343/4336、linkrods 3494/5078 …）⇒ 该改动对 OBJ 导出**零影响**。
- **旁支登记（未动手）**：`Perform` 的 `cxx:217-233`（line + 解析曲线的最短距离早退，`BRepExtrema_DistShapeShape(..., Extrema_ExtFlag_MIN)`，`d > 1.1 * myTol` 直接 return）**仍未移植**——已在 `perform` 就地注明。它是**性能**分支（不改结果），需要 `BRepExtrema_DistShapeShape` 的 MIN 模式（仓内 `brep_extrema` 可能只有默认模式），故不随本批夹带。
- **门禁影响面**：与 §2 基线逐项相同 ⇒ 记为"无影响"。

**批 98（`Geom_BSplineCurve::Segment` + `SetOrigin` + `InsertKnots`，并接线三个消费方；编译验证 + 数值探针）—— 2026-09-21 第 110 轮**

> 执行批 98 工作令（第 107 轮落定、第 108/109 轮补齐语义与行号）。按"能翻译成代码的就补、以编译作初步验证、未接指令不跑门禁/导出"执行。

- **新增 `occt-core/src/bspl/insert_knots.rs`**（`BSplCLib` 的"distinct 结点 + 重数"插入族，此前只有 `bezier::boehm_insert` 这个单跨 Boehm 步）：
  - `prepare_insert_knots` = `BSplCLib::PrepareInsertKnots`（`BSplCLib.cxx:1849-2024`）：`first/last` = `FirstUKnotIndex/LastUKnotIndex`（周期取 1/nbknots）、`adeltaK1/K2` 越界否决、`Eps=max(Tolerance,Epsilon(au))`、`Add=true` 与 `Add=false` 两种重数语义、周期首末同结点走 `aLastKnotMult` 双计；返回 `None` 对应 OCCT `false`（⇒ `Standard_ConstructionError`）。
  - `insert_knots_dim` = `BSplCLib::InsertKnots`（`:2063-2351`）逐句移植（维度主序平数组、周期极点环绕、`firstmult` 延后处理周期首末结点）。
  - 配套移植文件级 `Copy`（`:2028-2059`，含周期环绕取模）、`BoorIndex`（`:1810-1821`）、`BuildBoor`（`:1783-1806`）、`BoorScheme`（`:1010-1071`）、`GetPole`（`:1825-1845`，含 `2*Dimension` 交错存储与 `Position` 环绕）；`BuildKnots` 复用仓内 `bspl::build_knots`。
  - `insert_knots` = 模板 `BSplCLib_InsertKnots`（`BSplCLib_CurveComputation.pxx:380-438`）+ `PLib::SetPoles/GetPoles`（`PLib.cxx:115-201`，有理时按齐次 4 维打包/解包）。
- **`occt-geom/src/bspline_curve.rs`**（逐行移植，行号已核）：
  - `insert_knots(knots, mults, tol, add)` = `Geom_BSplineCurve::InsertKnots`（`:351-416`，含 `nbpoles == NbPoles()` 早退与 `updateKnots()`）；`insert_knot(u,m,tol,add)` = `:337-347`（默认 `M=1,tol=0,Add=true`，`hxx:239-242`）；`increase_multiplicity`（`:302-309`）、`increase_multiplicity_range`（`:313-323`）、`increment_multiplicity`（`:327-333`）。
  - `set_origin(index)` = `SetOrigin(Index)`（`:819-909`）；`set_origin_u(u,tol)` = `SetOrigin(U,Tol)`（`:913-970`，含折入周期、整体平移重参数化、最近结点 + `InsertKnot(U)` + `delta<0 ⇒ ik++`）。
  - `segment(u1,u2,tol)` = `Segment`（`:527-715`）逐句：`U2<U1`/周期 `DU-Period>PConfusion` ⇒ `DomainError`、两次 **flat 形** `LocateParameter`（按第 107 轮勘察：用 distinct 结点数组的 `Lower/Upper` 调 `locate::locate_parameter_range`，不用便利包装）、`Knots={min,max}`/`Mults={deg,deg}`/`Eps`、`InsertKnots(..., Add=false)`（`hxx:262-265` 的默认值）、周期分支 `SetOrigin`+`SetNotPeriodic`、`index1/index2` 右移规则、`PoleIndex` 取极/权、周期尾部结点改写。`u_knot_index_range` = `First/LastUKnotIndex`（`Geom_BSplineCurve_1.cxx:334-344`/`:404-414`）；`update_flat_knots` = `updateKnots()`（`:1171-1188`）。
  - **删除**批 108 点名的死方法 `GeomBSplineCurve::insert_knot(u, mult)`（无条件连插、零调用点、非 OCCT 语义）——忠实件已取代它，不再留两套语义（`bezier`/`knots::insert_knot` 导入随之移除）。
- **三个消费方接线**：① IGES `emit_bspline_curve`（`GeomToIGES_GeomCurve.cxx:309-356`）：无限界映射、单侧钳制、窄区间改走 `Segment` 副本（失败保持原曲线，与 OCCT 的 try/catch 一致），并把 126 的 `UMin/UMax` 由"请求值"改为 OCCT 实际写出的**钳制/裁剪后**值（`cxx:419`）；只有 Bezier 基曲线仍 UNPORTED。② `GeomConvert::CurveToBSplineCurve` 的 trimmed-BSpline 臂（`GeomConvert.cxx:322-339`）：拷贝 → `ElCLib::AdjustPeriodic`（整周时 `SetNotPeriodic`）→ `Segment`；trimmed-Bezier 臂仍 UNPORTED（缺 `Geom_BezierCurve::Segment`）。③ STEP 写侧 `write_context::emit_curve_entity`（`GeomToStep_MakeCurve.cxx:71-76`）：trimmed 的 B 样条基曲线先 `Segment(FirstParameter, LastParameter)` 再写。
- **验证（编译 + 数值探针）**：五 crate `cargo check --offline --all-targets` 全部 **exit 0 / 0 error**；临时探针（已删）实测：非周期 5 次插结 `poles 6→11`、几何漂移 **1.3e-15**；非周期 `Segment(0.25,0.8)` 区间与参数保持、漂移 **6.4e-16**；周期圆（`curve_to_bspline_curve(circle)`+`set_periodic`）`Segment(0.5,2.0)` 得 `periodic→false`、`poles 6→3`、漂移 **6.3e-16**；`SetOrigin(2)` 漂移 **2.1e-15**；`SetOrigin(U=1, Tol=Confusion)` 起点落在 `U`、插结 1 个（`poles 6→7`）、漂移 **9.6e-16**。
- **旁支上报（未动手）**：① 新增卡 **R2-23**（`Geom_BezierCurve::Segment` = `BSplCLib::BuildCache`（曲线版）+ `PLib::CoefficientsPoles`，是 trimmed-Bezier 三处（`GeomConvert.cxx:300-321`、`GeomToStep_MakeCurve.cxx:77-82`、`GeomToIGES_GeomCurve.cxx:441-448`）的最后前置）；② `GeomTrimmedCurve::bspline_knots()` 的"重映射结点"替身仍在，但 STEP 写侧现在优先走 `Segment`，只有 `Segment` 失败时才回落到它（`write_context.rs` 已注明）；③ 探针顺带确认端口 `circle_to_bspline_curve`+`set_periodic` 的曲线**逐点落在半径 2 的圆上**（9 个采样半径全为 2.000000），与 `gp_Circ` 的**参数化**不同属 `Convert_CircleToBSplineCurve` 的既知有理参数化（非同参数），非缺陷。
- **门禁影响面**：见 §6.2 第 98 行。

**批 97（IGES：周期 B 样条按 `SetNotPeriodic` 副本写 126；仅编译验证）—— 2026-09-21 第 106 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。本批把批 95 新增的 `SetNotPeriodic` 接到 IGES 写侧（`GeomToIGES_GeomCurve.cxx:294-307` 的原文控制流）。

- **改动（`iges.rs::emit_bspline_curve`）**：原实现遇 `curve.is_periodic()` 直接 `return None`（⇒ 调用方退化成弦线）；现按 OCCT `TransferCurve(Geom_BSplineCurve)`（`cxx:294-307`）改为**先做非周期副本**：取 `bspline_poles/knots/nurbs_degree/weights` 重建 `GeomBSplineCurve{periodic:true}` → `set_not_periodic()` → 用该副本走余下流程（范围检查、planar/normal、126 参数）。126 的 `periodic` 字段本就恒写 0（与 `IGESGeom_ToolBSplineCurve::WriteOwnParams` 一致），无需改。
- **仍 UNPORTED（就地注明）**：请求区间窄于曲线自身时 OCCT 用 `Geom_BSplineCurve::Segment` 取子段；该函数未移植 ⇒ 这类请求仍返回 `None`（下一批 **Segment + SetOrigin** 一并解决，两个消费方：此处与 `GeomConvert::CurveToBSplineCurve` 的 trimmed-Bezier/BSpline 臂）。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五 crate **exit 0、0 error**（本轮复检 core/geom/topo）。
- **门禁影响面**：仅 IGES **写**侧（周期 B 样条边：弦线 → 126），OBJ 门禁不读 IGES；`iges_check` 实体统计变化属预期（§6.2 第 97 行）。

**批 96（R2-22：2D/pcurve 周期表示 + 读侧接线；仅编译验证）—— 2026-09-21 第 105 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。本批把 3D 侧的周期机制（批 93/94）镜像到 2D，并为 R2-19 备好 2D 周期圆锥/样条。

- **2D 周期表示**（`occt-geom2d/src/bspline_curve.rs`）：加 `periodic` 字段（`new` 恒 false；全仓 29 处构造均走 `new`，无结构体字面量 ⇒ 无需改调用点）；覆盖 `is_periodic`/`period`；`continuity`/`parameter_intervals` 由写死的 `false` 改为 `self.periodic`；`d0/d1/d2` 的周期分支走 `occt_core::bspl::curve_dn::dn(u, 0/1/2, …, periodic=true, mults=None)`（2D 无权重 ⇒ 权重传 `None`，极点经既有 `poles_3d()` 抬到 z=0，取回 `(x, y)`）。**非周期分支一行未动**。
- **2D 周期化/反周期化**：`distinct_knots_and_mults` / `set_periodic` / `set_not_periodic` 与 3D 版同构（底层 `occt_core::bspl::{knots,locate,unperiodize}` 与维数无关），并已核 OCCT 的 `Geom2d_BSplineCurve::SetPeriodic/SetNotPeriodic`（`Geom2d_BSplineCurve.cxx:948`、`:1087`）与 3D 版同算法。
- **读侧接线**：把 `MakeBSplineCurveCommon` 的"公共段"（重复结点按 `Epsilon(|last|)` 归并、重数 > degree+1 钳制并裁剪首尾极点、`shouldBePeriodic` 判定）抽成与极点类型无关的 `bspline_descriptor`（3D helper 改调它，行为逐句不变），新增 `make_bspline_curve_2d_with_knots`（构造 + 闭曲线强制周期化，`IsClosed()` 公式同 3D，已核 `Geom2d_BSplineCurve_1.cxx:147-150`），并接入 2D `B_SPLINE_CURVE_WITH_KNOTS` 臂（`StepToGeom.cxx:952-963` 走同一模板）。**UNPORTED（就地注明）**：2D 有理复合体的权重被忽略（端口 2D B 样条本就无权重，`xs`/`ys` 双数组）。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；改动文件无新告警。
- **⚠️ 门禁影响面（§6.2 第 96 行）**：与批 94 同族的**读侧**改动，只是作用在 2D pcurve 上 ⇒ 周期 pcurve 会改求值路径与后续网格化（pcurve 是 `IntTools_FClass2d`/`BRepMesh` 的输入）。

**批 96 起手勘察（第 104 轮；未改代码，仅登记下轮入口）**

> 本轮回合预算用于收口批 95 并核对下两批的前置条件，未提交代码改动（五 crate 复检仍 exit 0 / 0 error）。

- **R2-22（2D pcurve 周期表示）已具备全部前置**（细节见该卡"批 96 起手勘察"）：`Geom2d_BSplineCurve::IsClosed()` 与 3D 同式（`Geom2d_BSplineCurve_1.cxx:147-150`）、`First/LastParameter` 用 `myFlatKnots`（端口公式已一致）、2D 无权重可直用 `curve_dn::dn`、周期化/反周期化可直接照搬 3D 版；读侧需先把 `make_bspline_curve_with_knots` 的公共段抽成与极点类型无关的形式。**属 live 路径**（pcurve 可能变周期）。
- **R2-6 剩余两臂的取径**：`GeomConvert::CurveToBSplineCurve` 的 trimmed-Bezier/BSpline 臂与 IGES 的"窄区间 B 样条"回落（`iges.rs:916-922` 现返回 None ⇒ 退化成弦线）**同一前置** = `Geom_BSplineCurve::Segment`（`Geom_BSplineCurve.cxx:527-715`）+ 其周期分支所需的 `SetOrigin`（`:819-…`）。已核 `Segment` 只依赖端口已有件（`locate_parameter` 带 distinct 重数、`PoleIndex`、`SetNotPeriodic`（批 95）、`InsertKnots`）+ 需移植的 `SetOrigin`；`InsertKnots(Knots,Mults,Eps)` 与端口 `insert_knot(u,mult)` 的**重数合并语义**需先对拍再接线。⇒ 建议批 97 = `Segment`+`SetOrigin`，同批接 IGES 窄区间与 R2-6 两臂（两个消费方同时兑现）。

**批 95（R2-20：IGES 整周椭圆走 B 样条 126；仅编译验证）—— 2026-09-21 第 103 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。本批把 T-78 余项 c（整周椭圆）收口，并补上它需要的四件库级控制流。

- **新增库件（均为 OCCT 同名成员）**：① `bspl::knots::{is_uniform_knots, reparametrize}`（`BSplCLib::KnotForm` `BSplCLib.cxx:602-632` + `BSplCLib::Reparametrize` `:756-798`：均匀串按常步长、非均匀按比例，含 `:788-793` 的严格递增护栏，`nextafter` 用 `from_bits` 实现）；② `GeomBSplineCurve::distinct_knots_and_mults`（端口只存平结串 ⇒ 周期曲线的"周期延拓结点"在 `[FirstParameter,LastParameter]` **之外**，区间内的 run-length 即 OCCT 单独存储的 `myKnots/myMults`，等价性在注释中论证）；③ `GeomBSplineCurve::set_knots`（= `Geom_BSplineCurve::SetKnots` `Geom_BSplineCurve.cxx:758-765`：`CheckCurveData` 的 `NbPoles(deg,periodic,mults)` 校验 + 按周期/非周期重建平结串）；④ `GeomBSplineCurve::set_not_periodic`（= `SetNotPeriodic` `:974-1019`：`PrepareUnperiodize`+`Unperiodize` `BSplCLib.cxx:2967-3080`，极点按 `NewPoles(k)=Poles((k-1)%n_old+1)` 环绕，权同样环绕）；⑤ `GpAx2::{rotate,rotated}`（= `gp_Ax2::Rotate/Rotated` `gp_Ax2.hxx:301-318`，用既有 `GpTrsf::set_rotation_ax1` + `rotate_vector`）。
- **`set_periodic` 顺带校正**：改用 `distinct_knots_and_mults` 后，已周期曲线的 `First/LastUKnotIndex` 周期臂（`1`/`NbKnots`）可正常走通并自然成为 no-op，故删掉批 93 的"提前 return"特判（更贴 OCCT）。
- **IGES 落地**（`iges.rs::emit_whole_period_ellipse`）：按 `GeomToIGES_GeomCurve.cxx:620-645` 逐句——旋转副本（`:626-628`，`gp_Ax3(pos).Direct()` 选角）→ `CurveToBSplineCurve(copystart, QuasiAngular)`（`:639`）→ `Reparametrize(Udeb, Udeb+2π, Knots)` + `SetKnots`（`:641-643`）→ `TransferCurve(BSpline)` 的周期→非周期转换（`:294-307`）→ 既有 126 发射器；`emit_conic_arc` 的整周分支由"返回 None（写 104）"改为调用它。**UNPORTED（就地注明）**：OCCT 先试的 `GeomConvert_ApproxCurve`（`:632-636`）⇒ 端口走回落臂（同一圆锥的精确有理形式，结点向量与 OCCT 的 approx 不同）。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；改动文件无新告警。
- **门禁影响面**：仅 IGES **写侧**（整周椭圆 104 → 126），OBJ 门禁读 STEP、不受影响；`iges_check` 的实体数与类型分布变化属预期（记 §6.2 第 95 行）。

**批 94（R2-21 读侧：STEP 读曲线接入 `MakeBSplineCurveCommon`；仅编译验证）—— 2026-09-21 第 102 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。**批 90（周期盒 arm）与本批共同构成本轮唯一 live 路径改动**。

- **新增 `make_bspline_curve_with_knots`**（`step/read_topology.rs`）= `StepToGeom::MakeBSplineCurveCommon`（`StepToGeom.cxx:776-927`）逐段移植：① 重复结点按 `k - last <= Epsilon(\|last\|)` 归并并累加重数（`cxx:784-821`）；② 重数 > `degree+1` 钳到 `degree+1`，记录首/尾差值（`:823-845`）；③ 按差值裁剪首尾极点与权（`:847-871`、`:902-906`，`NbUniquePoles <= 0` 报错）；④ `shouldBePeriodic` 判定（`:873-894`：总重数 ≠ `NbPoles+degree+1` ∧ 首尾重数相等 ∧ `Σmults - mults(1) == NbPoles`）⇒ 周期表示构造（`nb_poles(degree,true,mults)` 校验 + `knot_sequence_periodic` 平结串 + `periodic=true`），否则走原 `new`/`rational`（保留 `CheckCurveData` 的 `flat = n+d+1` 校验，等价于 OCCT 构造器抛异常）；⑤ 闭曲线强制周期化（`:920-926`：`closed ∧ degree>1 ∧ IsClosed()`，`IsClosed()` = 首末点平方距离 ≤ `Precision::Computational()`，`Geom_BSplineCurve.cxx:146-149`）。
- **新增 `curve_record_closed`**：解析 `closed_curve` 字段。端口有两种 layout——普通实例 `(name, degree, poles, curve_form, closed, …)`（closed 在 4）与合并后的有理复合体 `(name, degree, poles, weights, curve_form, closed, …)`（closed 在 5）；用槽 3 是否 `SELF`/`(…)` 判别（`transfer.rs::merge_complex_body:390-465` 已核）。
- **接入两族臂**：`B_SPLINE_CURVE_WITH_KNOTS`（文件自带 mults/knots）与 Bezier / uniform / quasi-uniform（按 `StepToGeom.cxx:310-317`/`:338-347`/`:362-384` 合成 mults+knots：`{0,1}×{d+1,d+1}` / `0..n+d` 全 1 / `0..n-d` 两端 `d+1`）。三族的合成后平结串与改动前逐字节一致（已逐族核对：Σmults 都等于 `n+d+1` ⇒ 非周期判定成立），差别只在于**新增**了 OCCT 的闭曲线周期化。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；改动文件无新告警。
- **⚠️ 门禁影响面（已并入 §6.2 第 94 行）**：读侧首次可能产出周期曲线（闭曲线强制周期化 ⇒ 极点数 −1、改走周期求值），**批 90 的 `box_bspline` 周期 arm 自此可达**；语料中 linkrods/screw/Shape/ATU01038/bottom/motoc/top 含 `closed_curve=.T.` 的 B 样条 ⇒ 这 7 个模型的曲线求值与包围盒都可能变，进而影响 `IntTools_EdgeEdge::BndBuildBox` 剪枝与网格门禁。
- **旁支上报（未动手）**：① 端口把**无结点的 `B_SPLINE_CURVE`** 也映射成曲线（`read_topology.rs:1256-1268` 用 `uniform_knots_for`），而 OCCT `StepToGeom::MakeBoundedCurve`（`StepToGeom.cxx:274-291`）**没有这一臂** ⇒ 该实体在 OCCT 侧是"未映射"（曲线为空）；端口更宽松，属既有偏差，登记待评估。② 2D 侧同一逻辑未接（新卡 **R2-22**）。③ `occt-geom/src/bspline_curve.rs` 的 `fd_d1`/`fd_d2` 既无调用者（既有死代码）。

**批 93（R2-21 库层：`Geom_BSplineCurve` 周期表示与周期求值；仅编译验证）—— 2026-09-21 第 101 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。本批解锁 R2-6 的非修剪圆/椭圆臂，并为 R2-20 备好周期 B 样条。

- **新增周期平结串**（`occt-core/src/bspl/knots.rs`）：`knot_sequence_length`（`BSplCLib::KnotSequenceLength`，`BSplCLib.cxx:455-474`：周期时 `l + 2*(degree+1-Mf)`）与 `knot_sequence_periodic`（`BSplCLib::KnotSequence` 周期臂，`:488-548`：基串写在 `M1+1..`，随后前后各绕一个周期；OCCT 的 `j` 越界在 Rust 里改为夹取，已就地注明）。
- **新增 `GeomBSplineCurve::set_periodic`**（`occt-geom/src/bspline_curve.rs`，= `Geom_BSplineCurve::SetPeriodic` `Geom_BSplineCurve.cxx:777-815`）：取 distinct 结点/重数 → `first/last_u_knot_index`（非周期用 `locate::{first,last}_u_knot_index` = `BSplCLib` 版 `:111-137`）→ 端部重数钳到 `degree` → `knots::nb_poles(deg,true,mults)` 重定极点数（缩则截断，涨则按 OCCT `Resize` 语义补 `GpPnt` 默认值/`0.0` 权）→ `knot_sequence_periodic` 重建平结串。**已周期者直接 return**：OCCT 此时 `First/LastUKnotIndex` 取 `1`/`NbKnots`、`updateKnots()` 重建同一串，本就是 no-op；端口只存平结串，其 run-length 分解会额外暴露"周期延拓结点"（OCCT 只放在 `myFlatKnots`），故显式取 no-op 并写明理由。
- **新增周期求值**：`GeomBSplineCurve::{d0,d1,d2}` 增加 `periodic` 分支，走 `bspl::curve_dn::dn(u, 0/1/2, …, periodic=true, mults=None)` —— 该函数即 `BSplCLib::PrepareEval`（`BSplCLib_CurveComputation.pxx:777-830`：`LocateParameter` 归入周期 + `BuildEval` 环绕极点窗口）+ `Bohm` + `RationalDerivative`；把它的导数下限从 `n<1` 放宽到 `n<0` 后，`n=0` 即 `BSplCLib::D0`（`BSplCLib_1.cxx:248-267`）。**非周期分支一行未动**（仍走既有 `eval_curve*`），故不影响已对齐几何。
- **`BSplineCurveBuilder` 补构造语义**：`Geom_BSplineCurve` 构造器会 `updateKnots()`，即按 `Convert.IsPeriodic()` 用**周期** `KnotSequence` 建平结串 —— 现在按此分派（周期用扩展串、非周期用原 `knot_sequence`）。
- **解除 UNPORTED**：`GeomConvert::CurveToBSplineCurve` 的非修剪圆/椭圆臂（`GeomConvert.cxx:363-408`）改为"转换 → `bspline_curve_builder` → `set_periodic()`"。至此该入口覆盖 OCCT 全部分派，仅剩三个依赖缺失的臂（`Segment`/`CompCurveToBSplineCurve`/`ApproxCurve`）。
- **随手去重**：`bspl::locate` 已有 `first_u_knot_index`/`last_u_knot_index`（`BSplCLib.cxx:111-137`），本批初稿在 `knots.rs` 重复实现，已删除改用既有件（只保留新增的 `nb_poles`）。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；改动文件无新告警（`bspline_curve.rs` 的 `fd_d1/fd_d2` never used 为**既有**死代码，本批未触碰，旁支上报）。
- **门禁影响面：0**。新周期代码目前只可能经 `GeomConvert::CurveToBSplineCurve` 进入，而仓内唯一调用点（`step/write_context.rs`）只在 Bezier 分支调用 ⇒ 圆/椭圆臂不可达；读侧（`StepToGeom.cxx:873-925`）未动，批 90 的周期盒 arm 仍不可达。
- **下一步（批 94）**：R2-21 读侧 = `StepToGeom.cxx:873-925`（`shouldBePeriodic` 判定 + 闭曲线强制 `SetPeriodic()`），并把"批 90 盒 arm / 周期求值同时变为可达"登记为门禁影响项；其后批 95 = R2-20。**起手勘察（本批顺带核对，避免下轮返工）**：① 端口 `read_topology.rs:1205-1239` 的 `B_SPLINE_CURVE_WITH_KNOTS` 臂**把重数立即展平成平结串**（`expand_knots`）再交给 `new`/`rational`，而 OCCT 是在**展平前**做 `:784-845` 的"重复结点归并 + 重数 > degree+1 钳到 degree+1 + 按首尾差值裁剪极点"、再按 distinct 重数判 `shouldBePeriodic`（`:873-894`）⇒ 移植时要先把 mults/knots 以 distinct 形式取出来（未展平），并补 `:784-871` 这段（含 `Epsilon(|lastKnot|)` 的归并容差与 `NbUniquePoles <= 0` 返回空）；② 周期描述子（`Σmults - Mults(1) == NbPoles`）**不能**走 `GeomBSplineCurve::new`（其 `check_degree` 要求 `flat = n+d+1`），须按周期约定用 `knot_sequence_periodic` 建平结串直接构造，或先非周期构造再 `set_periodic()`（等价性需按 `:907-918` 的构造器语义逐条核对）；③ `IsClosed()`（`Geom_BSplineCurve.cxx`：首末极点距离 ≤ `Precision::Confusion()`）与 STEP 的 `closed` 字段索引需在实现时以真实记录核对（端口 `read_topology.rs:1205-1216` 的两处 layout 注释不一致，需先确认）。

**批 92（R2-6 后半之一：`GeomConvert::CurveToBSplineCurve` 入口 + 写侧接线；仅编译验证）—— 2026-09-21 第 100 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。本批给批 91 的圆锥引擎接上 OCCT 入口，并把它接到 STEP 写侧。

- **新增 `GeomConvert::CurveToBSplineCurve`**（`occt-geom/src/convert_bspl.rs::curve_to_bspline_curve`，`GeomConvert.cxx:163-430`）：完整按 OCCT 的 trimmed / 非 trimmed 两段分派——trimmed 先取 `BasisCurve`/`FirstParameter`/`LastParameter` 并在基曲线非周期时钳制区间（`:171-189`），再依次判 line（`:191-206`：`Poles={StartPoint,EndPoint}`、`Knots={FirstParameter,LastParameter}`、重数 2、次数 1）、circle、ellipse、hyperbola、parabola、Bezier、BSpline、Offset；非 trimmed 依次判 ellipse、circle、Bezier（`:410-429` 夹紧 B 样条，重数 `degree+1`）、BSpline（`:431-434` `C->Copy()`：Rust 侧按同名字段重建，已在注释说明是"同数据拷贝、非重逼近"）、Offset，末尾 `Standard_DomainError("No such curve")`。
- **`BSplineCurveBuilder`**（`:60-83`）逐行移植：2D 极点升到 `z=0` → 用 `KnotSequence(knots, mults)` 展平 → `gp::Trsf::SetTransformation(TheConic->Position(), gp::XOY())` 变换；并复现 OCCT `CheckRational` 的判定（静态 `Rational(Weights)`，`Geom_BSplineCurve.cxx:97-108`，阈值 `gp::Resolution()` = 仓内 `REAL_SMALL`）与 `CheckCurveData` 的 `NbPoles(Degree, Periodic, Mults)` 校验（`:91-94`）。
- **配套移植**：`GpTrsf::set_transformation_from_to`（`gp_Trsf.cxx:172-194` 的两 `gp_Ax3` 重载；1 参重载已在仓内，命名区分并就地注明）、`bspl::knots::nb_poles`（`BSplCLib::NbPoles` `BSplCLib.cxx:392-451`，含周期/非周期两支与 OCCT 的 0 返回条件）。
- **写侧接线（live 路径）**：删除 `step/write_context.rs` 手写的 `bezier_to_bspline`，Bezier 分支改调新入口（`GeomToStep_MakeBoundedCurve.cxx:64-75` 的语义）。**行为保持**：非修剪 Bezier 的产物逐字段相同；trimmed Bezier 经 `untrimmed_basis()` 取基曲线（= 端口既有"写未裁剪基曲线"的已登记偏差，OCCT 侧是 `GeomToStep_MakeCurve.cxx:77-82` 先 `Segment`，`Segment` 仍 UNPORTED），即**不因本批丢几何**。
- **UNPORTED（就地注明 + `Err(Unported)`）**：trimmed Bezier/BSpline（缺 `Geom_BezierCurve::Segment` / `Geom_BSplineCurve::Segment` `Geom_BSplineCurve.cxx:527-660` 与 `SetNotPeriodic`）、`RationalC1` 且 `U2-U1>=6`（缺 `CompCurveToBSplineCurve`）、Offset（缺 `ApproxCurve`）、**非修剪圆/椭圆**（`cxx:378,383,406` 调 `SetPeriodic()` ⇒ 前置 R2-21）。
- **R2-21 勘察（本批追加，已并入该卡）**：把"周期曲线表示"要补的四段控制流定位到行——`KnotSequence` 周期臂（`BSplCLib.cxx:503-547`）/`NbPoles`+`First/LastUKnotIndex`/`d0,d1,d2` 忽略 `periodic` 且 `build_knots` 平结串臂亦然/`FirstParameter`+`LastParameter` 用非周期约定；外加读侧 `StepToGeom.cxx:873-925` 的 `shouldBePeriodic` 与"闭曲线强制 `SetPeriodic()`"。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；`convert_bspl.rs`/`write_context.rs` 无新告警。
- **门禁影响面**：**读侧无变化**（本批不碰 STEP 读路径）；仅 STEP **写**侧 Bezier 产物走新代码路径而结果逐字段相同（§6.2 第 92 行已按此登记）。

**批 91（R2-6 前半：圆锥→B 样条引擎 `Convert_*ToBSplineCurve`；仅编译验证）—— 2026-09-21 第 99 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。R2-6 的第一块（引擎）；消费方在批 92（`GeomConvert::CurveToBSplineCurve` 接线 → R2-20 整周椭圆 / R2-19 2D-UV 曲线）。

- **新增 `crates/occt-core/src/convert/conic_to_bspline.rs`**（`Convert_ConicToBSplineCurve.hxx:42-159` + `.cxx:32-785`、`Convert_ParameterisationType.hxx:36-46`）：`ParameterisationType`（8 值）、`CosAndSin`（五张输出表）、`ConicToBSplineCurve` 数据体 + 访问器（`cxx:51-151`）、两重载 `build_cos_and_sin`（`:359-622`）/`build_cos_and_sin_periodic`（`:626-785`），文件级 helper `AlgorithmicCosAndSin`（`:306-355`）与两个求值器 `CosAndSinRationalC1`（`:235-254`）、`CosAndSinQuasiAngular`（`:271-298`）。已移植参数化：`TgtThetaOver2` 与 `_1.._4`（纯三角，`:447-471`）、`QuasiAngular`（`PLib::NoDerivativeEvalPolynomial` + 仓内已有 `BuildSchoenbergPoints`/`Interpolate`/`EvalBsplineBasis`）、`RationalC1`（含周期重载的 14 极点/5 结点构造，`:671-783`，其中 `inverse` 翻转后**不重置**的 OCCT 原样语义照搬）。`BSplCLib::D0` 的有权/无权两条语义按 `BSplCLib_1.cxx:248-267` 就地实现（有权时除以权函数，与 `Geom_BSplineCurve::D0` 一致），未改既有 `bspl` 模块。
- **新增 `crates/occt-core/src/convert/conic_curves.rs`**：`Convert_CircleToBSplineCurve.cxx:47-172`、`Convert_EllipseToBSplineCurve.cxx:47-174`、`Convert_HyperbolaToBSplineCurve.cxx:32-79`、`Convert_ParabolaToBSplineCurve.cxx:32-69`。周期/带区间两套圆与椭圆（周期仅 `TgtThetaOver2`/`RationalC1`，`cxx:59-84`）、`signed_radius`（右手系取 `+r`，`cxx:88-99`）、`scaled_poles`（`gp_Trsf2d::SetTransformation(conic.XAxis(), gp::OX2d())`，`cxx:101-109`）、双曲线（中点权 `cosh((UL-UF)/2)`，`:62-71`）与抛物线（`Parameter() = 2*Focal`，权全 1，`:50-61`）。域错误按 OCCT 抛点：圆 `:127-130`、椭圆 `:127-131`、双曲 `:38`、抛物 `:37`（`Epsilon(0.)` = 最小正 double，`Standard_Real.hxx:242-246`，与仓内 `precision::epsilon(0.0)` 一致）。
- **补齐 `gp` 访问器**（均为 OCCT 同名成员）：`GpElips2d::{axis,x_axis,y_axis}`（`gp_Elips2d::Axis/XAxis/YAxis`）、`GpHypr2d::axis`（`gp_Hypr2d::Axis`）、`GpParab2d::{axis,parameter}`（`gp_Parab2d::Axis/Parameter`）。`convert/mod.rs` 文档更新：A9 指出的"非翻译坐标助手"仍在原位，OCCT `Convert/` 忠实件与之并列（**未**搬动既有助手，避免动 live 路径）。
- **UNPORTED（就地注明 + 返回 `Err(Unported)`）**：`Convert_Polynomial` 臂（`cxx:612-621`）依赖 `Convert_PolynomialCosAndSin.cxx:64-181`，其 `Locate`（`:29-62`）与 2D `BSplCLib::Trimming`（`BSplCLib_1.cxx:204-…`）仓内未移植；仓内无调用点请求该参数化。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；两个新文件 0 告警（首轮两处 `unused_assignments` 按"OCCT 声明即赋值"改为未初始化声明并就地注明）。
- **门禁影响面：无**（本批只新增模块，仓内尚无调用点）。"先落引擎、批 92 再接线"是有意的分批：`GeomConvert::CurveToBSplineCurve` 接到写侧后才进 R2-18-gate 影响面。
- **余项**：批 92 = `GeomConvert::CurveToBSplineCurve`（`GeomConvert.cxx:163-430`：trimmed-line、Bezier、`Geom_BSplineCurve::Segment` 臂、Offset→`GeomConvert_ApproxCurve` 臂、`BSplineCurveBuilder` 的 3D 提升与 `SetTransformation(Conic->Position(), gp::XOY())`）＋写侧改用忠实入口；其后 `GeomConvert_CompCurveToBSplineCurve`（`RationalC1` ≥6 弧度分段要用）→ R2-20/R2-19。
- **同批发现（新卡 R2-21）**：为批 92 备料时核对"圆锥结果如何落成 `Geom_BSplineCurve`"，发现**端口完全没有曲线周期表示**——`GeomBSplineCurve::periodic` 全仓无处置 `true`（构造器恒 `false`、无结构体字面量、`git grep '\.periodic = '\|'set_periodic'` 零命中；STEP 读侧 `read_topology.rs:1183,1186,1235,1266` 一律走 `new`/`rational`）。两条后果：① **批 90 的周期 B 样条包围盒 arm 目前不可达** ⇒ 该批在门禁上的实际影响面是 **0**（R2-18-gate 第②条与 §6.2 第 90 行已按此订正；实现本身仍忠实，等周期标志能进得来时才生效）；② `GeomConvert::CurveToBSplineCurve` 的非修剪圆锥臂要调 `SetPeriodic()`（`GeomConvert.cxx:378,383,406`）⇒ 批 92 必须同时移植 `Geom_BSplineCurve::SetPeriodic`（`Geom_BSplineCurve.cxx:777-815`：截到 `FirstUKnotIndex..LastUKnotIndex`、端部重数钳到 `degree`、`NbPoles(deg, true, mults)` 重定极点，端口 `occt-core::bspl` 的 `locate`/`build_knots` 已支持周期平结串，缺的只是这条构造/转换路径）。

**批 90（T-12 后半：周期 B 样条包围盒 arm；仅编译验证）—— 2026-09-21 第 98 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。

- **语料取证（先证明本批不是空转）**：`data/**/*.step` 里 `B_SPLINE_CURVE` 且 `closed_curve = .T.` 的模型有 **linkrods、screw、Shape、ATU01038、bottom、motoc、top**（粗匹配计数 9/1/6/135/34/25/41）⇒ 周期 B 样条**在门禁语料中真实存在**，本批必然进入门禁影响面。
- **实现（`geom_bnd_lib_curve3d.rs::box_bspline`）**：周期分支不再退回 `box_other`（33 点包络采样），改为与 OCCT 同构的 knot-span 采样：
  1. `ElCLib::AdjustPeriodic(first, last, PConfusion, u1, u2)`（`BSplineCurve.cxx:43-50`）；
  2. **不移植 `Segment`，而是证等价**：`Geom_BSplineCurve::Segment` 对周期曲线做的是"结点按整周期旋转 + 端部截断 + 重新参数化"（`Geom_BSplineCurve.cxx:527-660`），几何点集不变 ⇒ 分段曲线的 distinct 结点落在区间内的就是"原 distinct 结点 ± m·period"落在 `[u1,u2]` 内的那些；据此生成 cut 列表并 `sort`；
  3. 与非周期路径**共用**同一个 span 循环：`[u1,cut1] [cut1,cut2] … [cutn,u2]`，每段 `fill_box3d(.., n = degree)` 取最大挠度，`enlarge(1.5·maxDefl)` 后走 `reduce_spline_box3d(myGeom->Poles(), …)`（= OCCT `cxx:110` 用**原始**极点而非窗口）。
  非周期路径的等价处理（范围钳制 + 同一循环）保持原样；退化情形（`cuts` 为空/区间为一点）仍回落到"两端点盒"分支。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；本文件无新告警。
- **⚠️ 门禁影响面（已并入 R2-18-gate 第②条与 §6.2 预期影响表）**：周期 B 样条曲线（上述 7 个模型）的盒子由"33 点保守采样"变为 OCCT 的"逐 knot-span 采样 + 极点归约"，会更紧 ⇒ `IntTools_EdgeEdge::BndBuildBox`（批 80 的盒递归）与 `BRepBndLib::Add` 的剪枝随之变化。

**批 89（T-15 的 `CorrectParameter`：pcurve 结点吸附；仅编译验证）—— 2026-09-21 第 97 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。**R2-16 小项池第三个夹带项**（同为"原卡被 trait 缺口阻塞、现已可解"）。

- **解锁**：T-15 的 `CorrectParameter` 卡在"`Curve2d` 没有结点串（只有 `bspline_degree`/`bspline_poles2d`）"；`Geom2dBSplineCurve` 实际持有**扁平结点向量**（`pub knots: Vec<f64>`）⇒ 新增 trait 查询即可。
- **改动**：
  - `occt-geom2d/curve.rs`：新增 `bspline_knots2d(&self) -> Option<&[f64]>`（默认 `None`）；`bspline_curve.rs` 的 `Geom2dBSplineCurve` 返回 `Some(&self.knots)`（对应 `Geom2d_BSplineCurve::Knots()`）。
  - `shhealing/transfer_params.rs`：`correct_parameter` 由"直接返回 `param`"的占位实现改为忠实实现——先按 `Proj.cxx:258-266` 递归解包 `Geom2d_TrimmedCurve`/`Geom2d_OffsetCurve`（新增借用版 `correct_parameter_of` 以便递归），再按 `:268-279` 对 `Geom2d_BSplineCurve` 在 `Precision::PConfusion()` 内**吸附到首个** distinct 结点（用 `bspl::unperiodize::distinct_knots_and_mults` 取 distinct 串；端口 2D BSpline 无周期表示 ⇒ OCCT 的 `FirstUKnotIndex()..LastUKnotIndex()` 即全部 distinct 结点，已在文档注明）。
  - 文件头 UNPORTED 清单同步（该条改为"批 89 已移植"）；仍在列的 `CopyNMVertex`（端口无 `BRep_PointRepresentation`）与 `myLocation`（healing 恒 identity）保留原等价性论证。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**。
- **⚠️ 门禁影响面（已并入 R2-18-gate 第四条）**：`CorrectParameter` 经 `TransferParametersProj::transfer_range` 被 `wire_fix.rs:3068,3074`（`ShapeFix_Wire::FixNotchedEdges`）调用，而这条链正是 **T-69 诊断里的 `check_pcurves_and_shift` / `fix_lacking_all`** ⇒ 结点吸附会改变 notch 边的 pcurve 参数，可能影响 T-69 与网格门禁（T0M 等）。

**批 88（T-17：`ProjectAct` 的 `!OK` 分支补齐五臂；仅编译验证）—— 2026-09-21 第 96 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。**R2-16 小项池的第二个夹带项**（与批 87 同源：都是"原卡被 trait/helper 缺失阻塞、现已被 A20/A29/T-56 解锁"的补缺）。

- **补齐的臂**（`shhealing/shape_analysis_curve.rs::project_act`，逐行对照 `ShapeAnalysis_Curve.cxx:355-399`）：
  - `GeomAbs_Hyperbola`（`:376-380`）：`clib::hyperbola_parameter(&hypr.pos, major, minor, P)` + `clib::hyperbola_value`；
  - `GeomAbs_Parabola`（`:382-386`）：`clib::parabola_parameter(&parab.pos, P)` + `clib::parabola_value`；
  - `GeomAbs_Line`（`:388-392`）：参数走既有 `line_parameter`（= `ElCLib::LineParameter`，`ElCLib.cxx:1192-1195`），点走 `clib::line_value`；
  - `GeomAbs_Ellipse`（`:394-399`）：`clib::ellipse_parameter(&elips.pos, major, minor, P)` + `clib::ellipse_value`，并按 OCCT 置 `is_closed = true` / `period = 2π`（供 `cxx:479-484` 的周期修正）。
  - Circle 臂原有（`cxx:355-374`）⇒ **`!OK` switch 五臂与 OCCT 一一对应**；文件头的 UNPORTED 清单相应改为"已全移植"并写明解锁来源。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；本文件内两条 `value assigned ... never read` 告警经核对为**既有**（`proj_param`/`proj_distance` 的初始化值在所有路径上都先被覆盖，HEAD 同样如此）。
- **⚠️ 门禁影响面（已并入 R2-18-gate 卡，第三条）**：`ProjectAct` 经 `step/read_geometry.rs` 的边域投影进入 **STEP 读侧**（`ShapeAnalysis_Curve::Project` → `ProjectAct`），⇒ 圆锥/直线/椭圆曲线在"Extrema 投影失败"时不再退回默认段搜索，而是走 OCCT 的解析短路。门禁时与批 85/86/87 一起做「前/后」差分。

**批 87（T-12 前半：椭圆/双曲线/抛物线解析包围盒；仅编译验证）—— 2026-09-21 第 96 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。**本批是 R2-16 小项池里的夹带项**（R2 主序的 R2-19/R2-20 均为大件且不可验证，故先做这件已解锁的忠实补缺）。

- **解锁条件**：T-12 原卡被"`occt_geom::Curve` 不暴露 `gp_Elips/Hypr/Parab`"阻塞；这三个查询已在 A20/A29 轮补齐（`gp_ellipse`/`gp_hyperbola`/`gp_parabola`），⇒ 本批把 `geom_bnd_lib_curve3d.rs` 的 PARKED 前半落地。
- **新增三个忠实件**（逐行对照 OCCT）：
  - `box_ellipse_full` / `box_ellipse_range`（`GeomBndLib_Ellipse.cxx:23-45` / `:49-114`）：整椭圆按 `Amp = sqrt(Major²·Xd_k² + Minor²·Yd_k²)` 求每坐标极值；圆弧先 `AdjustPeriodic`＋加端点，再对每坐标解 `atan(MinR·Yk / MajR·Xk)` 的极值参数并用 `InPeriod` 判断落在弧内。
  - `compute_hyperbola_box` / `box_hyperbola_range`（`_Hyperbola.cxx:25-69` / `:75-140`）：端点、`t1·t2 < 0` 时的 `t=0`、每坐标 `T3 = 0.5·ln(|B−A|/|B+A|)` 极值（`|B±A| < Epsilon(1)` 即退化跳过的分支原样保留）；无穷参数的开盒分支（`Open*Min/Max`）逐支对应。
  - `box_parabola_range`（`_Parabola.cxx:23-94`）：端点、`u1·u2 < 0` 时的 `u=0`、无穷参数开盒分支。
  - 三处 `throw Standard_Failure("bad parameter")`（`_Hyperbola.cxx:82,108`、`_Parabola.cxx:33,59`）**不移植为 panic**：端口按既有约定返回空盒（`IsVoid`），已就地注明出处。
- **`box_curve` 分派**：按 `GeomBndLib_Curve.cxx:98-147` 的顺序补上 Ellipse → Hyperbola → Parabola 三臂（在 Line/Circle 之后、Bezier 之前）。
- **仍 PARKED**：周期 B 样条 arm（需 `Geom_BSplineCurve::Segment` + `AdjustPeriodic`），文件头保留原说明；H 侧取 `location()` 的临时量已按借用规则绑定（Rust-only 写法，无语义差异）。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**，新代码无告警。
- **⚠️ 门禁影响面（已并入 R2-18-gate 卡）**：`box_curve` 被 `brep_bnd_lib.rs:64`（`BRepBndLib::Add`）、`geom_bnd_lib_surface3d.rs:496,686`、`int_tools_curve_box.rs:46`（`IntTools_EdgeEdge::BndBuildBox` —— 批 80 的盒递归正是用它）与 `pcurve_full/surface_projector.rs:1493` 消费 ⇒ 圆锥曲线的盒子由采样变为**解析（更紧）**，会改变盒递归的剪枝与上游 bbox 结果。门禁时与批 85/86 一起做「前/后」差分。

**批 86（R2-18b：2D 求交改走 `IntAna2d_AnaIntersection`；仅编译验证）—— 2026-09-21 第 95 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。

- **`Curve2d` trait 补类型查询**：新增 `gp_elips2d`/`gp_parab2d`/`gp_hypr2d`（默认 `None`，对应 `Geom2dAdaptor_Curve::Ellipse()/Parabola()/Hyperbola()`，`Geom2dAdaptor_Curve.cxx:100-118`）；`Geom2dEllipse`/`Geom2dParabola`/`Geom2dHyperbola` 实现之，`Geom2dTrimmedCurve` 转发（`load` 解包到 basis 的语义），与 3D `Curve` trait 在批 20/A20 的补法同构。
- **`curve_ops::curve2d_intersections` 改为「解析优先 + 采样兜底」**：
  - `analytic_kind2d` 按 `Geom2dAdaptor_Curve` 顺序判 Lin/Circ/Elips/Parab/Hypr，`conic_of2d` 给出 `IntAna2d_Conic`（`from_lin2d`/`from_circ2d`/`from_elips2d`/`from_parab2d`/`from_hypr2d`）；
  - 分派到 `IntAna2dAnaIntersection` 的 8 个 `perform_*`：`(Lin,Lin)`/`(Lin,Circ)`/`(Circ,Circ)`/`(Lin,Conic)`/`(Circ,Conic)`/`(Elips,Conic)`/`(Parab,Conic)`/`(Hypr,Conic)`，其中 `(Circ,Lin)` 按 OCCT 的 `Perform(Lin,Circ)` 写法**交换实参**并在结果里把参数换回；
  - **有界曲线语义过滤** `param_on_curve2d`：`IntAna2d` 解的是**无界**圆锥曲线，越出曲线自身 `[first,last]` 的根不是这两条曲线的交点（周期曲线按一个周期映射；修剪曲线的 `first/last` 本就是 basis 空间值，直接比较即可）——这正是旧采样器隐含的语义，OCCT 侧由 `Geom2dAPI_InterCurveCurve` 的有界 adaptor 承担；
  - 非解析组合（B 样条/Bezier/offset 及其修剪）仍走原采样器，函数文档已标 UNPORTED + 忠实路线（`Extrema_ExtCC2d`）。
- **语义变化（如实登记）**：重合元素（同一直线/圆）现在返回 `IntAna2d` 的"重合"结论（0 个孤立点），旧采样器会给出一串重合点；这与 `IntAna2d_AnaIntersection::IdenticalElements` 的语义一致，但端口没有 `Geom2dAPI_InterCurveCurve::NbSegments` 那层 API 来暴露重合段。门禁时若某用例依赖重合点串，会看到差异。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；`curve_ops.rs` 无新告警（首版分派的 5 个 `(_, X)` 镜像臂是 unreachable pattern，已改为「首个曲线给专用操作数」的单层分派）。
- **⚠️ 门禁影响面（重要，已立 R2-18-gate 卡）**：`geom2d_api::{intersect_curves, project_point_on_curve}` 在 **live BOP 2D pcurve 路径**上被调用——`pave_de.rs:150,154`、`pave_blocks/split_edges.rs:174`、`wire_splitter_block.rs:680` ⇒ 批 85+86 会改变 PaveFiller/SplitEdge 的 2D 求交/投影输入，**可能改变 T-01/T-03/T-04 红门禁的结果（变好或位移都可能）**。因此门禁指令到位时，必须先跑「批 85+86 前/后」差分再解读任何红项。

**批 85（R2-18a：2D 投影改走 `Extrema_ExtPC2d`；仅编译验证）—— 2026-09-21 第 94 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行。

- **改动**：`occt-geom2d/curve_ops.rs::curve2d_closest_point` 的整段自创实现（64 点扫描取最优→黄金分割细化，无界曲线再用"扩窗探测"循环）**删除**，改为委托 `crate::extrema2d::point_curve_extrema2d`——即 `Extrema_ExtPC2d`（`Geom2dAPI_ProjectPointOnCurve::Perform` 的引擎，取最小距离解），直线/圆经 `Extrema_ExtPElC2d` 解析臂。私有 helper `refine_closest` 随之删除（`minimize_1d`/`sample_curve`/`curve_bbox` 仍被 `curve2d_intersections` 使用，保留）。
- **模块头登记（仍 UNPORTED，含出处）**：`curve2d_intersections`（256×256 采样＋交替一维极小化；忠实路线＝解析对走 `occt_core::intana2d::ana_intersection` 的 `perform_lin_lin`/`perform_lin_circ`/`perform_circ_circ`/`perform_{lin,circ,elips,parab,hypr}_conic`，其余走 `extrema2d::curve_curve_extrema2d_all` ⇒ **R2-18b**）；`curve2d_length`（Simpson 求积代替 `GCPnts_AbscissaPoint`/`math_GaussSingleIntegration`）；`extrema2d` 一般曲线的网格+Newton 播种（A7 家族，非本卡）。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**；`curve_ops.rs` 无新告警（顺手清掉测试模块里因此不再需要的 `GpVec2d` 导入）。
- **待验证风险**：有无界曲线（抛物线/双曲线）投影的既有用例会改走 `extrema2d`（其一般路径本身是 A7 家族替代件，播种方式不同）——门禁时若 2D 相关用例变化，第一嫌疑＝本批。`geom2d_api::project_point_on_curve`（在其上再套一层有界 Newton）语义不变。

**批 84（R2-3 周期面一半：128 反周期化 + `periodicU/V` + bounds 两臂；仅编译验证）—— 2026-09-21 第 94 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行 ⇒ 本批**未跑 `iges_check`**。

- **新增 `occt-core/src/bspl/unperiodize.rs`**（逐段对照 `BSplCLib.cxx:2967-3020` `PrepareUnperiodize`、`:3024-3080` `Unperiodize`，以及面版 `BSplSLib.cxx:2331-2380`）：`distinct_knots_and_mults` / `flat_knots_from_mults` / `prepare_unperiodize` / `unperiodize_knots` / `unperiodize_direction`。`unperiodize_direction(degree, flat_knots) -> (new_flat_knots, pole_map)`：`pole_map[k] = k % n_old`，即 `BSplCLib::Unperiodize` 的环绕 `NewPoles(k) = Poles((k-1) % len + 1)` 沿该方向的等价形式（`SetPoles` 的 U 行主序 / V 列主序两种展平下同式，模块头已证；weights 与 poles 同一映射，因 OCCT 对**齐次**极点做反周期化）。
- **`iges.rs::emit_bspline_surface` 改动**：① `period_u/period_v` 取**原始** `bs.is_u_periodic()/is_v_periodic()`（`GeomToIGES_GeomSurface.cxx:235-236`），并按 `IGESGeom_ToolBSplineSurface::WriteOwnParams` 的字段序写 `closedU, closedV, polynomial, periodicU, periodicV`（`:64-125` 核对）；② 周期方向在读取 knots/poles 前**反周期化**（先 U 再 V，U 扩展会同步扩展每一行的 weights，V 再逐行扩展）；③ `closedU/V` 改为对反周期化后的极点/权重比较（`iso_rows_equal`，`cxx:347-348`），不再用 `is_periodic()` 充当代理；④ bounds 修正补齐**两条臂**：非周期钳制到自身 bounds，周期则"端点已在边界则吸附"＋`AdjustToPeriod` 平移＋超过一个周期时截断（`cxx:244-302`；`pcurve_full::adjust_to_period` 由 `pub(super)` 提为 `pub(crate)`）。
- **仍 UNPORTED（就地注明出处）**：`SetUOrigin`/`SetVOrigin` 重定原点（`cxx:310-320`/`:330-340` → `Geom_BSplineSurface_1.cxx:1026-1130`）。仅当 `AdjustToPeriod(Umin,U0,U1) != AdjustToPeriod(Ufin,U0,U1)`（区间跨越周期原点）时与 OCCT 有差异，端口保留曲面自身结点原点。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**（`iges.rs` 内除既有 `point_entities` 死字段告警外无新告警）。
- **旁支**：`R2-3` 的另两半另立卡——**R2-19**（非平面面的 2D/UV 曲线 + `PreferenceMode` 3；平面面现状已忠实，因 `BRepToIGES_BRWire.cxx:374-377` 直接返回空）、**R2-20**（整周椭圆改走 `GeomConvert_ApproxCurve`，前置 R2-6）。
- **待验证风险（等门禁指令）**：`iges_check` 结构自洽（128 的 `indU/indV` 因反周期化而变大、knots/poles/weights 计数随之变化，属预期）＋ `occt-topo --lib` 的 iges 用例＋四道 STEP 门禁与 `export_data_obj`；若 128 出现自洽性失败，第一嫌疑＝环绕映射的极点/权重对齐。

**批 83（R2-4 / T-85 步 2：IGES 写侧可达性过滤 + 重编号；仅编译验证）—— 2026-09-21 第 93 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试/导出"执行 ⇒ 本批**未跑 `iges_check`**（它是导出+结构校验，等指令）。

- **求根**：`IgesWriter` 新增 `roots: Vec<usize>`，`emit_shape` 的 Solid/Shell/Face 分支把 `emit_solid/emit_shell/emit_face` 的顶层 DE 收进去（Compound 递归 ⇒ 多根）。对应 OCCT 侧 `IGESControl_Writer::AddShape` 后 `AddEntity(ent)` → `myModel->AddWithRefs(ent, protocol)`（`IGESControl_Writer.cxx:243-252`）。
- **可达集与编号**：`final_entities()` 做**先序 DFS**——`AddWithRefs` 的实际语义（`Interface_InterfaceModel.cxx:652-692`：先 `AddEntity(自身)`，再按 `FillSharedCase` 的迭代序递归，已加入者跳过；`level=0` 表示不限层）。子序严格照 `IGESData_GeneralModule::FillSharedCase`：**DE 字段 3–8/13 的实体在前**（本写入器只有字段 7＝变换矩阵可能非空，即 `Ent::trsf`），**然后**才是自身参数指针（`Ent::refs` 的记录序）。⇒ 编号序＝DFS 访问序，**不是**设计行里写的"创建序"，已在 R2-4 行与代码注释就地订正。
- **参数指针重写（占位符方案，行内已批准为备选）**：`emit_refs` 的 10 处发射点（102/120/122/142×3/144/192/194/196/198/402）改为写 `#k` 占位（`ph`/`ph_list`），`finish()` 先求 `entities = final_entities()`，其中 `resolve_placeholders` 把 `#k` 换成 `new_of[refs[k]]`，同时重映射 `Ent::trsf`；DE 段与 P 段都遍历重编号后的 `entities`（P 段行数与 64 列断行因此基于**替换后**的文本）。`emit_refs` 里加了 debug 断言：占位符下标必须落在 `refs` 范围内。
- **审计过的风险点**：全部 `self.emit(...)`（无指针实体：116/110/108/100/104/124/126/128/123）确认不含指针字段，指针实体一律走 `emit_refs` ⇒ 重编号不会漏改；共享实体的 DE 指针可能指向更早的号（OCCT 的 `AddWithRefs` 同样如此），IGES 允许。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**（含清理本批引入的两处 `unused variable`）。
- **旁支发现（只报告不修，門禁 4）**：`IgesWriter::point_entities` 是**死字段**（`HEAD` 起就只有声明与初始化、无任何读取）——点实体去重表已不再被使用，登记待清理。
- **待验证（等门禁指令）**：`iges_check` 全量对比（预期实体数下降、`unreferenced` 归零——批 76 统计的孤儿集中在 `Shape 8`/`Shape-2 28`/`linkrods 28`/`Offset 33`/`Shape-1 76`/`ATU01038 159`/`occ/bottom 184`/`occ/top 211`）；`occt-topo --lib` 的 iges 用例（行数 80、含 144/100 等断言不受影响）、四道 STEP 门禁与 `export_data_obj`。若出现孤儿残留或指针越界，第一嫌疑＝`FillSharedCase` 子序（trsf 在前）或某发射点漏用占位符。

**批 82（R2-5：STEP 读侧最后一条——`TransferRelatedSRR` + SRR `TransferEntity` + `ComputeSRRWT`；仅编译验证）—— 2026-09-21 第 93 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试"执行 ⇒ 本批**未跑任何 `cargo test` / 导出**。

- **新增 `SsrComposer`（`step/read_geometry.rs`）**，逐段对照 `STEPControl_ActorRead.cxx`：
  - `rep_shape(rep)` = `TransferEntity(ShapeRepresentation, …)`（`cxx:2048-2092`）：先取该表示的自身形状（`resolve_representation`，多 item 即 compound），再跑 `transfer_related_srr` 把相关形状并入 `aCund`；非空的"最后一个相关结果"覆盖单结果（`cxx:2077-2080`），`aCund` 恰 1 个子 → 绑该形状、>1 → 绑 compound（`cxx:2082-2089`）。结果按表示 id 记忆化（= OCCT `TP->Bind`），供相关表示递归取用。
  - `transfer_related_srr(rep, cund)`（`cxx:2679-2729`）：遍历与本表示"共享"的 `SHAPE_REPRESENTATION_RELATIONSHIP[_WITH_TRANSFORMATION]`（`aGraph.Sharings` 的实体序 = 按实体号排序的解析记录），`nbrep = (Rep1 == rep ? 2 : 1)`，每个结果 `aBuilder.Add(theCund, …)`，**返回最后一个**（`cxx:2724-2728`）。
  - `transfer_srr(srr, nbrep)`（`cxx:1337-1433`）：`nbrep` 0/1/2 决定取 `Rep1`/`Rep2`/两者；`iatrsf = ComputeSRRWT(...)`；`nsh==1` 对单结果、`nsh>1` 对 compound 施加变换（`cxx:1408-1417`）。
  - `compute_srrwt(srr)`（`cxx:2495-2546`）：`CARTESIAN_TRANSFORMATION_OPERATOR_3D` → 复用既有 `Resolver::make_transformation3d`（`cxx:2522`）；`ITEM_DEFINED_TRANSFORMATION` → 两 `AXIS2_PLACEMENT_3D` → `Resolver::compute_axis_transform`（`cxx:2484-2489`，本轮把该函数由私有提为 `pub(super)`）；两者在 `Form() == Identity` 时返回"无变换"（`cxx:2527`/`:2489`）。
  - 接入点：`read_step_impl` 的表示循环由"每个表示各自出一个根形状"改为 `composer.rep_shape(id)`；`read.step.shape.relationship` 默认 true（`DESTEP_Parameters.hxx:170`）而端口无该参数开关 ⇒ 恒执行（就地注明）。
- **诚实登记（仍 UNPORTED，均就地写清出处）**：`PrepareUnits(Rep2)`/恢复与 `ComputeTransformation` 的逐表示单位上下文（`cxx:2517-2526`、`:2467-2482`）；`ComputeTransformation` 的轴归属/交换告警（`cxx:2440-2464`）；`TransferRelatedSRR` 的 `MECHANICAL_DESIGN_AND_DRAUGHTING_RELATIONSHIP` 与构造几何两臂（`cxx:2707-2721`）；装配级 `TransferEntity(NAUO,…)`（`cxx:859-981`，含 `CDSR` 与 `SRRReversed`）——`data/` 内无此类实体。Rust-only：SRR 环保护（OCCT 的 `TP->Bind` 在转移之后，环会无限递归）。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**。
- **行为影响**：无 SRR 关系的文件路径与旧实现等价（`transfer_related_srr` 返回 `None`、`cund` = 自身形状）；有 `SHAPE_REPRESENTATION_RELATIONSHIP[_WITH_TRANSFORMATION]` 的文件由"相关表示各自成根"变为**组合进产品形状并按关系施加变换**，门禁待指令后跑（该行为变更正是本卡目的）。

**批 81（R2-2 / A16 求交 half：生产调用点改走忠实件 + 删除端口自创采样器；仅编译验证）—— 2026-09-21 第 93 轮**

> 按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试"执行 ⇒ 本批**未跑任何 `cargo test` / 导出**。

- **曲线–曲线（`inttools/intersections.rs:278,285,292,295`）**：4 处 `geom_api::curve_curve_intersections` 调用改为新 helper `edge_edge_general_points(e1, e2)` → **`crate::edge_edge::EdgeEdge::with_edges(...).perform()`** 的 `points()`（= OCCT `IntTools_EdgeEdge::Perform` 的 VERTEX 解；批 80 刚落地的 `FindSolutions` 递归）。OCCT 对这三类非解析组合只走这条路（`IntTools_EdgeEdge.cxx:185-243`）；容差由边容差 + fuzzy 决定（`Prepare`），调用方的 `tol` 本就不是 `IntTools_EdgeEdge` 的入参，已在 helper 文档写明。line/line、line/circle、共面 circle/circle 的既有精确臂保持不变（其 OCCT 缺失部分＝R2-17）。
- **非平面边–面（`:450`）**：改为新 helper `curve_surface_points(curve, face, surf, a, b)` → **`crate::intcurvesurface::perform_curve_surface`**（= `IntCurveSurface_HInter`）。对齐 `IntTools_EdgeFace.cxx:426-445`：OCCT 建 `GeomAdaptor_Curve`/`GeomAdaptor_Surface` 后 `Perform`，再按 `aPoint.W() ∈ [aTF, aTL]` 过滤——端口把区间过滤留给调用方的 `param_on_edge`（同一判据），曲线窗口取 adaptor 自然范围（无界曲线退回边区间，OCCT 的 polygon 路径同样要裁剪），表面窗口取 `int_curves_face::finite_uv(face)`（该模块的 `IntCurvesFace_Intersector` 约定；`finite_uv` 已提为 `pub(crate)`）。
- **删除的端口自创采样器（`occt-geom/src/geom_api.rs`，连私有 helper 与自带单测）**：`curve_surface_intersections`（曲线采样 + 距离凹陷 + 黄金分割）、`curve_curve_intersections`（256×256 网格 + 交替一维极小化）、`curve_curve_distance`（128×128 采样 + 坐标下降）、`surface_bbox`（`nu×nv` 栅格 bbox）、`minimize_1d`/`dist_curve_surface`/`curve_window_for_bbox`/`sample_curve`；单测 `line_sphere_intersections`/`crossing_lines_intersection`/`parallel_lines_distance` 随之删除（**未新增任何测试**）。保留并仍忠实：`project_point_on_curve`（`Extrema_ExtPC`）、`project_point_on_surface`（`Extrema_ExtPS`）；保留未动：`curve_bbox`/`pcurve_of_curve_on_surface`/`distance_point_*`/`tangent_at`（bbox 与投影辅助，非求交；其近似性归 A1/A15 家族，另行登记）。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**（`geom_api` 删除后无残留引用；`inttools` prelude 补 `Surface`）。
- **新卡 R2-18（A16 的 2D 余项，诚实登记）**：`occt-geom2d/curve_ops.rs` 的 `curve2d_intersections`/`curve2d_closest_point` 等采样族仍被 `geom2d_api.rs:6` 生产使用，忠实件（`extrema2d` 与 `IntAna2d`）需另批接线，未在本批擅动。
- **待验证风险**：非平面边–面与通用曲线对的几何结果会变（采样 → 忠实引擎），本次未跑门禁；若 `phase*/step_*` 出现回归，第一嫌疑＝本批两处 helper，其次 R2-1/R2-17。

**批 80（R2-1 / T-77：`IntTools_EdgeEdge::FindSolutions` 参数盒递归忠实移植；仅编译验证）—— 2026-09-20 第 93 轮**

> 承接新一轮清单（§3 R2）。按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试"执行 ⇒ 本批**未跑任何 `cargo test` / 导出**。

- **新增 `crates/occt-topo/src/edge_edge/find_solutions.rs`**（逐段对照 `D:\source\OCCT-src\src\ModelingAlgorithms\TKBO\IntTools\IntTools_EdgeEdge.cxx`，模块头附 OCCT 行号映射表）：`TypeToInteger`(`:1456-1482`)、`ResolutionCoeff`(`:1486-1559`)、`Resolution`(`:1561-1607`)、`CurveDeflection`(`:1611-1638`)、`IsClosed`(`:1642-1659`)、`SplitRangeOnSegments`(`:1366-1406`)、`BndBuildBox`(`:1410-1419`，经 `int_tools_curve_box::add_curve_to_box` = `BndLib_Add3dCurve::Add` 忠实件)、`PointBoxDistance`(`:1423-1452`)、`DistPC`/`FindDistPC`(`:1210-1362`，含 `iC=±1` 的极大/极小语义与 `0.618…` 黄金分割)、`FindParameters`(`:553-671`)、`FindSolutions`(`:290-349` + 递归 `:353-549`)、`MergeSolutions`(`:675-779`)、`AddSolution`(`:780-825`)、`FindBestSolution`(`:826-901`)、`CheckCoincidence`(`:1150-1206`)、`IsIntersection`(`:1060-1146`)。
- **`Prepare` 补齐**（`edge_edge.rs`）：`GeomAdaptor_Curve::load` 顺序的精确类型判定（Circle→Line→Ellipse→Parabola→Hyperbola→Bezier→BSpline→Offset→Other，**不做几何采样**）、`TypeToInteger` + 等类型时的挠度比较与 `--iCT1`、`mySwap` 边交换（`myEdge1/2`、`myCurve1/2`、`myRange1/2` 同步交换）、`myResCoeff1/2`、`myRes1/2`、`myPTol1/2`（`>999` 时 `5e-16*aTM`）。
- **删除的替代件**：`edge_edge.rs` 的 `find_solutions`（"`Extrema_ExtCC` 极值 + 容差即交点"）、其 `find_parameters`（Newton polish `locate_extcc`）与 `merge_solutions`（`remove_identical_roots`/`sort_roots` 去重）**整体删除**，prelude 同步清理。⇒ **A11 的"替代判据"项结项**（T-77）。
- **配套忠实件**（新依赖，皆为 OCCT 有据的移植而非新规则）：`Curve::offset_curve()`（`Geom_OffsetCurve::BasisCurve/Offset`，用于 `ResolutionCoeff`/`Resolution` 的 offset 臂；`GeomTrimmedCurve` 两处转发）；`GeomBezierCurve::resolution`（= `Geom_BezierCurve::Resolution` 经 `BSplCLib::Resolution`，用 `deg+1` 个 0/1 结点序列）；`BndBox::{is_x_thin,is_y_thin,is_z_thin}`（= `Bnd_Box::IsXThin/IsYThin/IsZThin`，`FindParameters` 的 `bThin` 判据）。
- **数据模型映射（已就地登记）**：`TopAbs_EDGE` 解 → `CommonPrt{part_type: Edge, range: Range1}`；`TopAbs_VERTEX` 解 → 一条 `PntOn2Faces`（`uv1/uv2` = `FindBestSolution` 在**调用方**边序上的参数，按 `mySwap` 交换）。`IntTools_CommonPrt::{Range2, Edge1, Edge2, BoundingPoints}` 端口不存（既有结构裁剪，无消费方）。
- **本批明确保留（新卡 R2-17）**：`perform` 里 `(Circle, Circle)` → `compute_circle_circle_full` 的分派在 OCCT **不存在**（`Perform` 只特判 line/line）；忠实递归已能覆盖圆–圆，但删除会改行为 ⇒ 与门禁同批做。就地已标 UNPORTED + `cxx:185-243` 出处。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0、0 error**（`occt-core`/`occt-geom` 的新增 helper 亦 0 error）；`edge_edge` 内只剩既有的 `solvers::*` glob 告警。
- **待验证风险（诚实登记，等门禁指令）**：`edge_edge::{tests,tests_full}` 中直接依赖旧替代件的断言（尤其 `line_bspline_crossing_two_hits` 现在会走 `mySwap`（line↔BSpline 交换）与盒递归、`circle_circle_two_hits`、`matches_inttools_*`）**未跑**；若门禁出现回归，第一嫌疑是 R2-1，其次 R2-17。行号引用不变（新增文件，未改既有行）。

**批 79（全仓命名重构：122 个 `pNN.rs` → 按内容命名；仅编译验证）—— 2026-09-20 第 93 轮**

> 用户指令：把所有 `p01`、`p02` 这类无语义命名改成可读命名。按"能翻译成代码的就补、以 `cargo check` 作初步验证、未接指令不跑测试"执行 ⇒ 本批**未跑任何 `cargo test` / 导出**。

- **范围**：全仓 `p\d+\.rs` 共 **122 个**，分布在 60 个目录（`viz_scene` 6、`fillet_curved` 5+2 测试、`step` 5、`bop_curved` 4、`delaun` 4、`brep_gprop_full` 4、`pcurve_full` 4、`pave_intersect` 4、`extrema_cc` 4、`shhealing` 4、`viz_scene/tests` 2 …）。命名来源 = 每个文件的顶层条目清单（`grep '^pub'`）+ 文件内分节横幅（`// ----` 标题），按该文件承载的 OCCT 类/职责取名，例如：
  - `step/{p01..p05}` → `step/{format, write_context, transfer, read_topology, read_geometry}.rs`（`format` = 序列化原语与实体写出、`write_context` = `WriteCtx` + B 样条拟合、`transfer` = `STEPControl_*` 级读写入口 + `Resolver` 状态、`read_topology`/`read_geometry` = 读侧拓扑解析与几何解析）
  - `extrema_cc/{p01..p04}` → `{poly_roots, curve_curve, glob_opt_func, general_extrema}.rs`；`extrema_surf/{p01..p03}` → `{analytic_solvers, numeric_extrema, point_surface_extrema}.rs`；`meshing/delaun/{p01..p04}` → `{constants, triangulation, frontier, polygon_meshing}.rs`；`bop_curved/{p01..p04}` → `{region_mesh, face_meshing, general_mesh, region_trim}.rs`；`brep_offset/{p01,p02}` → `{curve_face_offset, shell_offset}.rs`；`fillet_curved/{p01..p05}` → `{surface_info, meridional, torus_blend, rolling_ball, corner_patch}.rs`
- **做法**：① `git mv` 122 个文件；② 按"同目录映射表"重写 `mod pNN;` / `pub use pNN::*;` / `super::pNN::` / 文档链接里的模块路径（`pub use pNN::*` 保持不变式：各部分共享同一命名空间，调用方看到的仍是 `<module>::<item>`，例如 `occt_topo::step::read_step_file`）；③ 同步 `specs/**` 的 `文件:行` 引用（先按"目录尾段 + 文件名"机械替换，再逐条手改无目录前缀的裸引用与跨目录同名引用，如 `p01.rs:620-621` → `extrema_surf/analytic_solvers.rs:620-621`、`p04.rs:42-70,325-401` → `bop_curved/region_trim.rs:42-70,325-401`、`p04.rs:346-375` → `meshing/delaun/polygon_meshing.rs:346-375`）；④ 把该约定写进 `CLAUDE.md`（禁止再出现 `pNN`）。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五个 crate **全部 exit 0**；重命名后 `git grep -E "p[0-9]{1,2}(\.rs|::)"` = **0**，`git grep -E "[A-Za-z_]+/p[0-9]{1,2}"` = **0**（残存的 `p1`/`p2`/`p00` 命中全是局部变量与 Bézier 控制点，与本次无关）。**未跑测试/未导出**（待指令）。
- **事故与修复（如实登记，见 §9）**：第一遍脚本的映射表里新名带了 `.rs` 后缀，替换时只吃掉 `pNN` 部分 ⇒ 产出 `mod mesh_queries.rs;` / `poly_roots.rs::f_and_jac` / `foo.rs.rs` 三类畸形，**由门禁 1（编译）当场抓到**（`occt-core` 3 个 `expected one of ; or {`）。修复 = 对 122 个新名做四类归一：`N.rs.rs::` → `N::`、`N.rs.rs` → `N.rs`、`mod N.rs;` → `mod N;`、`N.rs::` → `N::`（含 `super::N.rs::`）；再编译 exit 0。历史行号未变（重命名不改文件内容行数），specs 里所有 `文件:行` 仍然有效。

**批 78（STEP 读取侧忠实补齐（续）：椭圆轴交换、EDGE_LOOP 顶点绑定两趟、CURVE/SURFACE_REPLICA）—— 2026-09-20 第 92 轮**

> 继续按"**能翻译成代码的就补**、以 `cargo check` 作初步验证、未接指令不跑测试"执行。本轮同样**未跑任何 `cargo test` / 导出**。

- **① 椭圆轴交换（3D＋2D）**：原臂注明"`majorR < minorR` 分支未移植"。按 `StepToGeom::MakeEllipse`（`StepToGeom.cxx:1536-1566`）补齐：`majorR = SemiAxis1*LF`、`minorR = SemiAxis2*LF`，当 `majorR - minorR < 0` 时 `A.SetXDirection(A.XDirection() ^ A.Direction())` 并交换半径（`:1552-1561`）；2D 版按 `MakeEllipse2d`（`:1571-1598`）无 LengthFactor、负差时 X 方向取 `gp_Dir2d(X.X(), -X.Y())`（`:1591-1593`）。为 3D 版新增 `GpAx2::set_x_direction`（= `gp_Ax2::SetXDirection`，重算 `Y = Z ^ X`，`occt-core/src/gp/ax2.rs`）。
- **② `EDGE_LOOP` 顶点绑定两趟（原为 UNPORTED，`step/read_topology.rs`）**：按 `StepToTopoDS_TranslateEdgeLoop.cxx` 逐行移植
  - **第 1 趟 `:288-403`（bug PRO7656）**：每条 ORIENTED_EDGE 取其 `EDGE_CURVE`（沿嵌套 ORIENTED_EDGE 下钻，`:301-307`）、按 `same_sense` 定 `Vstart/Vend`（`:355-365`）、记录调用前的 `IsBound` 状态（`:367-368`）；两顶点都在且点距 ≤ `Precision::Confusion()` 时按三分支绑定 `Vend→V1` / `Vstart→V2` / `Vend→V1`（`:384-396`）。
  - **第 2 趟 `:405-491`（bug BUC50070 #3815）**：相邻边对 (j, j+1)，按 `Orientation` 取各自"相接顶点" `Vs1/Vs2` 与 `Vs11/Vs22`（`:429-433`），四者有二同一实体则跳过（`:435-438`）；两点距 ≤ `Precision()`（端口已有 `step_precision`）时，`EC1` 未翻译则绑 `Vs1→V2`、否则 `EC2` 未翻译则绑 `Vs2→V1`（`:466-477`）。
  - 端口侧新增 `Resolver::vertex_bind: RefCell<HashMap<usize, TopoShape>>`（= `StepToTopoDS_TranslateTool::Bind`），由 `resolve_shape` 对 `VERTEX_POINT` 优先返回；"是否已翻译/已绑定"以 `shape_cache.contains_key` 复现（`aTool.IsBound`）。`resolve_loop` 现在**先**跑两趟绑定、再解析 ORIENTED_EDGE，等价于 OCCT 的"建 wire 前先 confuse 顶点"。
- **③ `CURVE_REPLICA` / `SURFACE_REPLICA`（此前无臂）**：按 `StepToGeom::MakeCurve`（`:1351-1371`）与 `MakeSurface`（`:1967-1986`）实现：解析 `ParentCurve/ParentSurface` → 递归转移 → `MakeTransformation3d(Transformation)` → `Transform`；含 `PC != SC` 的循环保护（`:1358`/`:1973`，端口以 `parent_ref == id` 判并报错）。为此把 `make_transformation3d` 提为 `pub(super)`。
- **验证（仅编译）**：`cargo check --offline --all-targets` 五 crate 全 ok。改动规模：4 文件 **+208 / −22**（`step/read_topology.rs` +201、`step/read_geometry.rs` +16、`step/transfer.rs` +4、`occt-core/gp/ax2.rs` +9）。
- **STEP 读取侧余下可翻译项（只剩一条）**：`STEPControl_ActorRead::TransferRelatedSRR`（`STEPControl_ActorRead.cxx:2061-2091`）＋ `TransferEntity(SRR,…)`（`:859`/`:985` 一带，含 `ComputeSRRWT` `:2495-2546`：`SHAPE_REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION` 的 `CARTESIAN_TRANSFORMATION_OPERATOR_3D` 或 `ITEM_DEFINED_TRANSFORMATION` 两 `AXIS2_PLACEMENT_3D` 两臂）——即把 `data/` 之外可能出现的 `SHAPE_REPRESENTATION_RELATIONSHIP[_WITH_TRANSFORMATION]` 关系**组合**进产品形状（端口现在把每个相关表示各出一个根形状；`read.step.shape.relationship` 默认 true，`DESTEP_Parameters.hxx:170`）。端口已具备所需件：`make_transformation3d`、`compute_axis_transform`（= `ComputeTransformation` 的轴对臂 `:2471-2489`）、`make_compound_of`；**UNPORTED 将保留**：`PrepareUnits(Rep2)` 的逐表示单位上下文（端口用文件全局因子）与 `SRRReversed` 的 NAUO 校验（`STEPConstruct_Assembly::CheckSRRReversesNAUO`）。⇒ **已于批 82（R2-5）落地**：`SsrComposer` 移植了 `TransferRelatedSRR`/SRR `TransferEntity`/`ComputeSRRWT`；上述两项 UNPORTED 与 MDADR/构造几何臂、装配级 `TransferEntity(NAUO,…)` 就地登记（详见 §7 批 82）。**STEP 读取侧至此无可翻译余项**。

**批 77（STEP 读取侧忠实补齐：RECTANGULAR_TRIMMED_SURFACE、Bezier/Uniform/QuasiUniform 曲面与曲线族、TRIMMED_CURVE 3D/2D）—— 2026-09-20 第 91 轮**

> 本轮按"**能翻译成代码的就补**、以**编译**作初步验证、未接到"导出 obj 测试"指令前**不自行发起测试**"执行：只跑 `cargo check`（5 个 crate × `--all-targets` 全 ok），**未跑任何 `cargo test` / 导出**。

- **① 曲面 `RECTANGULAR_TRIMMED_SURFACE`**（`step/read_geometry.rs`，原为"无此臂"的 UNPORTED）：按 `StepToGeom::MakeRectangularTrimmedSurface`（`StepToGeom.cxx:1834-1885`）实现 —— 先解析基面（`:1838`），再按基面类型取 `uFact/vFact`（球/环面 `AngleFact`；柱面 `u=AngleFact, v=LengthFact`；回转面 `u=AngleFact`；锥面 `v=LengthFact/cos(SemiAngle)`；平面 `两者=LengthFact`；其余 `1`，`:1845-1875`），把 `U1/U2/V1/V2` 分别乘以因子后交给 `GeomRectangularTrimmedSurface`（`:1882`，端口用 `rectangular_trimmed::uv`）。**UNPORTED 就地注记**：`Usense/Vsense` 与 `SetTrim` 的归一化（含周期基面的 `ElCLib::AdjustPeriodic`）端口包装类型不带。
- **② 曲面族 `BEZIER_SURFACE` / `UNIFORM_SURFACE` / `QUASI_UNIFORM_SURFACE`（原完全无臂）**：按 `MakeSurface`（`StepToGeom.cxx:522-743`）把三者转成 B 样条：Bezier→结点 `{0,1}`、重数 `deg+1`（`:537-548`）；Uniform→`n+deg+1` 个结点 `i-1`、重数全 1（`:570-593`）；QuasiUniform→`n-deg+1` 个结点 `i-1`、两端重数 `deg+1`（`:612-638`）；`_AND_RATIONAL_B_SPLINE_SURFACE` 复合体再带权重（`:644-687`、`:692-739`）。**忠实保留 OCCT 的例外**：Bezier 臂不读 `WeightsData`（`:524-555`），故有理 Bezier 面的权重被丢弃 —— 与 OCCT 相同。
- **③ 曲线族 `BEZIER_CURVE` / `UNIFORM_CURVE` / `QUASI_UNIFORM_CURVE`（原无臂）**：同样规则（`StepToGeom.cxx:295-320`、`:322-350`、`:352-384`），`_AND_RATIONAL_B_SPLINE_CURVE` 复合体带权重（`:385-421`、`:422-458`）；OCCT **没有** Bezier+rational 臂 ⇒ 有理 Bezier 曲线的权重同样被丢弃（`BezierCurve` 臂不读）。
- **④ 复杂实体合并**（`step/transfer.rs::merge_complex_body`）：原只识别 `B_SPLINE_*` 家族；现扩展为同时识别上述 6 个家族（成员名 + `RATIONAL_B_SPLINE_{CURVE,SURFACE}` 权重成员），并按各家族解析器期望的布局给出 args（曲线 `['',deg,poles,weights,form,closed,self]`；曲面 `['',degU,degV,poles,weights,form,closedU,closedV,self]`）。
- **⑤ 3D `TRIMMED_CURVE`**（原为"原始参数直接用"的 UNPORTED）：完整移植 `StepToGeom::MakeTrimmedCurve`（`StepToGeom.cxx:2323-2496`）+ `ExtractParameter`（`:2221-2317`）：MasterRep 映射（`.CARTESIAN.`=1/`.PARAMETER.`=2/其余 0）、`isPoint` 判定（`:2352-2373`）、按基曲线类型的 `fact/shift`（Line＝`Dir()->Magnitude()*LengthFactor`，端口用 `resolve_vector().modulus()`；Circle/Ellipse＝`PlaneAngleFactor` ＋ 椭圆 `SemiAxis1<SemiAxis2` 时 `+π/2`，`:2375-2392`）、无 ref_direction 圆锥的整周回退（`:2394-2425`，新增 `axis2_has_ref_direction` 读 `AXIS2_PLACEMENT_3D` 的 `$`）、非周期区间钳制（`:2438-2459`）、`|t1-t2|<PConfusion` 时周期分支走 `ElCLib::AdjustPeriodic`（`:2462-2465`，用已移植件 `elib::clib2d::adjust_periodic`）、`SenseAgreement` 的实参次序（`:2486-2493`）。同轮新增 `parse_trimming_select`（`step/format.rs`，拆出 `PARAMETER_VALUE(...)` 与 `#point`）。**UNPORTED 注记**：闭合（非周期）基面的分支（`Curve` trait 无 `IsClosed`）与 `Geom_TrimmedCurve` 的 `Sense` 标志（端口只存区间）。
- **⑥ 2D `TRIMMED_CURVE`**（原 UNPORTED）：按 `MakeTrimmedCurve2d`（`StepToGeom.cxx:2505-2563`）：基曲线已是 2D B 样条则**原样返回不裁剪**（`:2514-2517`）；两次裁剪必须都是单参数选择（`:2523-2525`）；`fact/shift` 同 3D 但 Line 用**未乘 LengthFactor** 的 `VECTOR` 幅值（`:2533`）；抛物线/双曲线分支在 OCCT 里本身是 TODO（`:2547-2551`）。**UNPORTED 注记**：`Geom2dConvert::CurveToBSplineCurve(theTrimmed)`（`:2560`）与 2D 的 `SenseAgreement`。
- **验证（仅编译，按本轮指令）**：`cargo check --offline --all-targets` 五个 crate 全部 ok（occt-core/occt-math/occt-geom/occt-geom2d/occt-topo）。**未运行任何测试或导出**。
- **STEP 读取侧余下可翻译项（下一轮继续，均已定位）**：① `StepToTopoDS_TranslateEdgeLoop.cxx:288-403` / `:405-491` 的顶点绑定/重绑（`step/read_topology.rs:218`）；② `STEPControl_ActorRead::TransferEntity` 的 mapped-item 组合（`step/read_geometry.rs:788`）；③ 线/圆/椭圆臂的 `majorR<minorR` 轴交换（`step/read_topology.rs:980`）、`STEPControl_ActorRead` 的 `CURVE_REPLICA`/`SURFACE_REPLICA`（如缺臂）；④ 3D/2D trimmed 的 `Sense`/`CurveToBSplineCurve`（受 A8 的 `GeomConvert` 阻塞）；⑤ 写侧 `step/write_context.rs:148`（无边 3D 曲线时的 B 样条拟合）与 `step/write_context.rs:216`（`SetNotPeriodic`）。

- **验证**：本轮**未改行为代码**（两处探针 + `examples/zz_probe_t69b.rs` 全部还原，`git grep dbg-t69` 仅剩画板文字引用、代码 0 残留）；`occt-topo --lib` **1293/1**。
- **同轮追加的源码级线索（下一步直接照做）**：
  1. **`StepToTopoDS_TranslateEdgeLoop::CheckPCurves`（`cxx:100-176`）的"周期窗口归一"未移植**：`git grep adjust_periodic -- crates` 显示端口**从不**在 STEP 读入路径调用 `ElCLib::AdjustPeriodic`（该函数本身已移植在 `occt-core/src/elib/clib2d.rs:40-65`）。OCCT 在该处对每条边做：`sae.PCurve(edge,face,pc,w1,w2,false)` →（非周期 pcurve 时把 w1/w2 夹进 `[cf,cl]`）→ **`if (w1 > w2 && mySurf->IsUPeriodic()) { ElCLib::AdjustPeriodic(u1,u2, min(|w2-w1|/2, PConfusion), w1, w2); B.Range(edge, face, w1, w2); }`**。对 face 20 逐边推算：edge 59（`7.854→1.571`，**w1>w2**）与 edge 62（`1.5708→−4.712`，**w1>w2**）会被归一到 `[1.5708, 7.854]`（=与 edge 60 同窗口）；这正是"两条 wire 落在同一 u 窗口"所需的对齐，而端口现在把 STEP 的原始窗口**原样留着**（所以 wire 0 在 `[1.5708,7.854]`、wire 1 的 61/62 在 `[−4.712,1.5708]` 与 `[1.5708,7.854]` 混杂）。
  2. **同族的已登记缺口**：`shhealing/pcurve_ranges.rs::project_wire_pcurve_ranges`（`TranslateEdgeLoop.cxx:844` 的 `EdgeProjAux` → `B.Range`）的注释已自述"窗口按原样写入（与 OCCT 同）"，并**推迟了 pcurve 反向**，理由是"我们的 mesher 处理不了**反向且窗口与 3D range 周期错位**的 pcurve"——face 20 正是这种形状，说明该缺口与 T-69 是同一处。
  ⇒ 第 43 轮执行顺序建议：**(i)** 先做诊断实验（把同面各 wire 的 u 归一到同一周期窗口，预期 face 20 变 74/72，用它锁定机制）；**(ii)** 再按 (1) 移植/补齐 `CheckPCurves` 的 `AdjustPeriodic`+`B.Range` 分支（并核对 (2) 里推迟的 `need_reverse`），**不得在网格层加对齐补丁**；**(iii)** 重跑 `occt-topo --lib` + 四道 STEP 门禁 + `export_data_obj`，按面类记录 T0M 计数。
**批 76（T-85 收益量化：`iges_check` 增加"未被引用实体"统计，18/18 模型给出孤儿清单）—— 2026-09-20 第 90 轮（收尾）**

- **落地（仅示例，库未动）**：`examples/iges_check.rs` 在原有四项结构检查之外，按各类型的指针布局收集**被引用集合**，输出"未被引用的 DE 数 + 按类型分布 + 前 6 个 DE"（根实体合法地不被引用，另作说明）。这一改动把 T-85 的收益从"预期"变成**可测量**。
- **实测（18/18 模型仍然 `ok`；`unreferenced` 数）**：`Cube 1`、`Cylinder 1`、`Cone 1`、`rev 1`、`Extrusion 1`、`HoledPlate 1`（＝仅根 `402`，无孤儿）；`Sphere 5`（根 144 ＋ 孤儿 `142/100×2/123`）；`Torus 2`；`OffsetPlaneHoleEdge 2`；`screw 4`；`Shape 9`（孤儿 `142×8` ＋ 根 `402`）；`Shape-2 29`（孤儿 `142×28` ＋ 根）；`linkrods 28`；`Offset 33`；`Shape-1 76`；`ATU01038 159`；`data/occ/bottom.step 184`；`data/occ/top.step 211`。
- **结论（写进 T-85）**：孤儿**绝大多数是 `142`**（及少数 `100`/`123`）——正是"`isWholeSurface` 面的外侧轮廓"这一类：144 把外侧指针写成 `0`（`IGESGeom_ToolTrimmedSurface::WriteOwnParams` 的 if/else 分支），OCCT 因**只写从根可达的实体**而不会把该 142 写进文件，端口则全部写出。规模：`Shape-2` 会少 28 个 142、`occ/top.step` 会少 200+ 个实体 ⇒ T-85 步 2（已定稿四步方案）一旦实施，这些孤儿即消失，且**不会**影响任何被引用的数据（本轮的引用集合与 `iges_check` 全绿同时证明了这一点）。
- **门禁**：本轮未改库代码 ⇒ `occt-topo --lib` **1287/1**（红 = T-01）不变；`cargo check --all-targets` 通过；`iges_check` 18/18 `ok`。

**批 75（第 89 轮收尾：IGES 校验器全量跑通 18/18 ＋ §14.4 交接）—— 2026-09-20 第 90 轮**

- **验证**：把批 74 的常驻校验器跑到**全部**导出模型 —— 16 个 `data/*.step` ＋ `data/occ/{bottom,top}.step` 共 **18/18 `ok`**（结构自洽：80 列、P 段序号连续、P 行 DE 指针 = `2i−1`、每个 DE 的 `pstart/pcount` 可取串、102/142/144/402/192/194/196/198/120/122 的引导指针全部可解析、T 卡四项 = 实际段长）。最大两份：`occ/bottom.step` DE=2526 / P=12974、`occ/top.step` DE=2627 / P=15128。
- **交接**：新增 **§14.4「第 89 轮收尾交接」**，汇总本会话批 59–74（A26/IGES 主体、A15 部分、A8/T-44、A12/T-48）与三条已设计未实施的下一步（T-85 步 2 可达性过滤＋重编号；T-78 的 2D UV 曲线／周期面反周期化／整周椭圆；T-77 与 A16 求交 half 为下一批最合适取材），并写明仍被前置阻塞的 A12/A13/A18/A19/A23/A15 与 T-79/T-80 链。
- **门禁**：本轮无库代码改动（`occt-topo --lib` 1287/1 不变；四道 STEP 门禁与 `export_data_obj` 见批 72/74）。

**批 74（IGES 结构自洽校验器转为常驻示例 `examples/iges_check.rs`）—— 2026-09-20 第 89 轮**

- **动机**：IGES 文本**没有任何门禁覆盖**，批 70/71/72 的正确性只能靠临时探针（用完即删，下次无法复现）。本轮把批 70 的临时探针写成**常驻示例**（不是单测，不进 `cargo test`，与 `export_data_obj` 同类），把"写入器必须满足的结构不变量"固化下来，供后续 IGES 批次一条命令复核。
- **检查项（逐项对应 OCCT 出处，写在文件头注释里）**：① 每张卡 80 列；② P 段序号连续 `1..N`、每行第 73 列节字母、其 DE 指针域 = 所属实体首张 D 卡号 `2i−1`（`IGESData_IGESWriter.cxx:903`）、每个 DE 的 `pstart`/`pcount`（`:834-835`）能取到该实体自己的参数串；③ 各复合实体**引导指针域**逐类型可解析（102 曲线表、142 surface/curve3d、144 surface/outer+inner、402 实体表、192/194/196/198 point/axis/refdir、120 axis/generatrix、122 directrix）；④ Terminate 卡四项 = 实际段长（`:942-947`：`nbs`/`nbg`/`nbd*2`/P 末序号）。
- **用法**：`cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example iges_check -- Cube Sphere Shape Shape-2 HoledPlate ATU01038`（缺省跑 `Cube Sphere`；参数为 `data/<name>.step`，含 `/` 者按路径用；有问题时打印前 6 条并 `exit(1)`）。
- **实测**：6/6 模型 `ok` —— `Cube DE=37 P=37`、`Sphere DE=12 P=24`、`Shape DE=67 P=551`、`Shape-2 DE=209 P=4249`、`HoledPlate DE=231 P=278`、`ATU01038 DE=3144 P=6850`（实体数较批 70 那次普遍减少，正是批 71 删除孤儿 116 的结果）。
- **门禁**：`occt-topo --lib` **1287/1**（红 = T-01）——示例为增量文件，不影响库与门禁；`step_obj_parity`/`step_to_obj`/`step_obj_area`/`step_geometry_parity`/`export_data_obj` 在批 72 已复测为 14/14·13/13·11/11·2/3·16/16，本轮未触碰库代码。

**批 73（T-85 步 2 设计定稿：可达性过滤＋重编号的可执行方案）—— 2026-09-20 第 88 轮**

- **本轮无代码改动**：把 T-85 步 2 的完整方案写进任务行（含 4 个子步骤与"参数里指针 vs 实数"的判别规则），使下一轮可直接照做而无需重新推导。要点：
  1. **求根**：`emit_shape` 的 Solid/Shell/Face 分支把 `emit_solid/emit_shell/emit_face` 返回的顶层 DE 收进 `IgesWriter::roots`（Compound 走递归 ⇒ 天然多根，与 `BRepToIGES_BRSolid::TransferCompound` 的逐子形状转移一致）。
  2. **可达集**：从 `roots` 沿 `Ent.refs` DFS（`refs` 已在批 72 于 10 处复合发射点全部记录）。
  3. **重编号**：按创建序给可达集编 `new[old] = 1..N`，DE 卡号与 P 行的 DE 指针域（`2i-1`）一律用新号。
  4. **参数指针重写**：`params` 中"指针 vs 实数"的判别＝`num()` 恒带 `.`/`e` ⇒ 实数；类型号恒在首位；其余**纯整数字段**按各类型布局即指针（102：`n` 后 n 个；142：第 2、4 位；144：第 1 位 surface、第 4 位 outer、第 5.. 位 inner，`flag`/`n` 是计数；402：`n` 后 n 个；192/194/196/198：point/axis 与末位 refdir；120：axis/generatrix；122：directrix），按 `refs` 的**记录顺序**逐位回填（备选：构造 `params` 时对指针位留占位、写文件时替换）。
- **现状核对**：工作树在 `17f21a7`（批 72 文档提交）处干净；`occt-topo --lib` **1287/1**（红 = T-01），与 §14.3 的基线一致；IGES 侧上一轮已把孤儿 116 清除（批 71）、T 卡 D 计数修正（批 70）、写入器内部不变量（`check_refs`，批 72）。
- **判断**：T-85 步 2 属**写入器级重构**（需改动"写 DE/P 段"的编号与指针回填路径），在剩余轮次内实施并**完整验证**（结构自洽探针 + 全语料检查）风险偏高；故本轮先把方案定稿，留给下一轮或后续会话按上述 4 步执行，避免半途改动导致写出的 IGES 指针错乱。其余未结项（T-78 的 2D UV 曲线/周期面反周期化/整周椭圆、T-77、A16 求交 half、T-67 步 3、T-79/T-80、A12 余额、T-69 阻塞项）均已在 §14.3 列明入口。

**批 72（A26/T-78 步 12 / T-85 步 1：写入器参数指针结构化 —— `refs` 元数据 + 内部一致性断言）—— 2026-09-20 第 87 轮**

- **背景**：T-85（写侧"只写可达实体"＋重编号）的前提是把**参数里的指针**变成结构化数据；端口此前把指针值直接拼进 `params` 文本，无法安全重写（文本替换会误伤坐标/结点）。
- **落地（无输出变化）**：
  - `Ent` 新增 `refs: Vec<usize>`（＝ OCCT 的 `OwnShared` 语义：该实体引用的其它实体，按记录顺序）；新增 `emit_refs(ty, form, params, &refs)` 记录之。
  - **迁移全部复合发射点**：196/192/194/198 的 `point/axis/refdir`、102 的曲线表、142 的 `surface/curve3d`、144 的 `surface/(outer)/inner…`、402 的实体表、120 的 `axis/generatrix`、122 的 `directrix`（共 10 处，与 T-85 行登记的一致）；叶子实体（116/108/110/123/126/128/104）无外指，继续用 `emit`。
  - `finish()` 开头调用新增的 `check_refs()`：**debug 断言**每个记录的指针都落在 `1..=实体数` 且不自指——把上一轮临时校验探针查结构不变量的做法固化进写入器本身（OCCT 的 `Send(handle)` 会按模型解析指针，悬空数字会被原样写出）。
- **验证**：`check_refs` 在整个测试语料（含 Box/Sphere 的 IGES 写出、`Shape`/`Shape-2` 的 draw 导出）上**全部通过**，即 10 处发射点的指针记录与实际写出的记录一致、无越界/自指。门禁：`occt-topo --lib` **1287/1**（红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16**；IGES 6 个单测全绿。
- **T-85 余项**：在此 `refs` 之上做 ① 从根（`write_shape_iges` 的每个 shape 顶层实体）DFS 求可达集；② 按新序号重写各实体参数里的指针域（`params` 目前是文本，可选做法：把指针域在构造时留占位、写文件时用 `refs` 顺序回填）；③ 收益＝去掉 `isWholeSurface` 面外侧 142 之类的游离实体。

- **下一步（第 43 轮，两条并行）**：① **判定性实验**：把同一面上各 wire 的 u 归一到同一周期窗口（仅作诊断，不改 `node_insertion` 的行为）后看 face 20 是否走通忠实 Delaunay（预期 74/72）——若走通，则缺的是"跨 wire 的周期对齐"，须先在 OCCT 里找到执行该对齐的那段控制流（候选：`ShapeFix_Wire::FixShifted` 的面级用法 / `StepToTopoDS_TranslateEdgeLoop` 的 `CheckPCurves` 段 / `ShapeFix_Face::FixMissingSeam`）并按它移植，**不得直接给网格层打补丁**；② 用参考的重复顶点证据确定"共享该顶点的两个面"是哪两张，确认端口是否两张都没网格化（据此判断要补的控制流范围）。
**批 71（A26/T-78 步 11：删除未被引用的 116 point 块；登记"可达性过滤"需先做参数结构化）—— 2026-09-20 第 86 轮**

- **缺口**：`emit_shape` 对每个 Solid/Shell/Face 先调 `emit_vertices`，把该形状所有顶点无条件写成实体 **116**——这套 116 是**无人引用**的输出。OCCT 的 116 只在其被引用处出现（192/194/196/198 的**定位点**、以及 `BRepToIGES_BRVertex` 在旧 shell 路径里的点），Faces 模式下立方体这类纯平面/直线组合**不会**写出任何 116。
- **落地**：删除 `emit_vertices` 方法与其三处调用（`emit_shape` 的三个分支），并清掉随之失效的 `vertices_of` import；就地注释说明 116 的出现条件。`box_iges_contains_expected_entities` 的期望实体清单相应去掉 `"116"`（该断言把旧的孤儿 116 当成必需项，属编码旧行为的断言；`110/108/142/144/402` 保留）。
- **可达性过滤（本轮做了方案评估，**未实施**，理由入画板/T-85）**：OCCT 的 `IGESData_IGESWriter::Write(root)` 只写从根**可达**的实体并按遍历序**重新编号**。端口要照做，必须先让参数具备**结构化的指针元数据**（当前指针是拼进 `params` 文本里的数字，跳过实体就会导致 DE 号位移、而文本里的指针值无法安全重写——直接文本替换会误伤坐标/结点）；即需要把 `emit` 改为接收 `refs: &[usize]` 并在写文件时按新序号重写指针域。已立 **T-85**（含需改的 10 处复合实体发射点清单），并注明其收益＝去掉 `isWholeSurface` 面外侧 142 之类的游离实体（本轮的孤儿 116 已单独清除）。
- **门禁**：`occt-topo --lib` **1287/1**（红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16**；IGES 6 个单测全绿。

**批 70（A26/T-78 步 10：修正 T 卡 Directory 计数 off-by-one；IGES 结构自洽性全量验证）—— 2026-09-20 第 85 轮**

- **发现（临时校验探针查出，属 批 56 卡片格式合规化时留下的缺陷）**：Terminate 卡的 **Directory 计数写成 `2·实体数 + 1`**。`finish()` 里 `d_seq` 初值为 1、每实体 +2，收尾时却直接写 `d_seq`；OCCT 写的是 `nbd * 2`（`IGESData_IGESWriter.cxx:942-947` 的 `Sprintf("S%7dG%7dD%7dP%7d…", nbs, nbg, nbd*2, thepnum.Value(...)-1)`），且同一 `Sprintf` 的 S/G/P 三项取的正是**各段最后一个序号**（端口已用 `s_seq-1`/`g_seq-1`/`p_seq-1` 对齐）⇒ D 项同样应为 `d_seq - 1`。实测 `Shape.step`：实际 D 段 28 行，T 卡却写 29。
- **落地**：`finish()` 的 T 卡改为 `d_seq - 1`，并把注释补全为 `cxx:942-947` 的三项语义（注释原引 `:942-943`）。
- **验证（临时探针 `zzprobe_iges_validate.rs`，已删；`git grep zzprobe` = 0）**：对 `Cube`/`Sphere`/`Shape`/`Shape-2`/`HoledPlate`/`ATU01038` 逐份写出文件并做**结构自洽检查**：① 每张卡 80 列；② P 段序号 1..N 连续、每行第 73 列为 `P`、其 DE 指针域等于该实体的首张 D 卡号（`2i−1`）；③ 每个 DE 的 `pstart/pcount` 都能在 P 段取到完整参数串；④ 复合实体的**引导指针域**逐类型复核（102 的曲线表、142 的 surface/curve3d、144 的 surface/outer+inner、402 的实体表、192/194/196/198 的 point/axis、120 的 axis/generatrix、122 的 directrix）全部指向存在的 DE；⑤ T 卡四项与实际段长逐项相等。修前 6/6 模型报 "D 计数 +1"；修后 **6/6 全部 `ok`**（如 `ATU01038: DE=3729 P=7436 D=7458`、`HoledPlate: DE=291 P=339 D=582`）。该探针同时**证明批 56/59/60/62/63 的 DE 指针与 P 段记账自洽**（此前无任何门禁覆盖 IGES 文本）。
- **门禁**：`occt-topo --lib` **1287/1**（红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16**；IGES 自身 6 个单测全绿（`section_terminators_present` 未受影响——它只检查 T 卡存在与终止格式，未钉住 D 计数）。
- **旁支（登记，未改）**：写入器仍无 OCCT 的"只写可达实体"过滤 ⇒ `isWholeSurface` 面的外侧 142 会成为游离实体（合法但多一个）；建议作为 T-78 余项之一处理。

**批 69（A26/T-78 步 9：面内"不属于任何 wire 的边"补成内侧 142 轮廓）—— 2026-09-20 第 84 轮**

- **缺口**：`BRepToIGES_BRShell::TransferFace` 在两条 wire 循环之后还有一段（`BRepToIGES_BRShell.cxx:334-365`）：把面上**不属于任何 wire 的边**逐条转移，并各自封一条 `CurveOnSurface`（142）追加到内侧轮廓序列（2D 曲线同样先由 `TransferEdge(edge, face, originMap, length, false)` 拿，本端口未移植该分支 ⇒ 取 OCCT 的"仅 3D"支 `PreferenceMode = 2`，与 wire 轮廓一致）。端口此前只遍历 wire，漏掉这一段。
- **落地**：`iges.rs::emit_face` 在 wire 循环之后新增该分支——用 `children_of_type(face, Edge)` 取面的直接子边、以 wire 内边集合（按 TShape 指针去重，含缝边两次出现的情形）排除，剩余每条边复用 `edge_curve_entities` 缓存发曲线实体并封 `142,0,{surf},{0},{curve3d},2;` 进内侧轮廓（同时登记进 `curve_refs`）。注释写明出处与"2D 曲线 UNPORTED"。
- **实测**：本端口的面在现有模型里都由 wire 组成，故该分支在 `data/*.step` 上**取不到候选**（无行为变化）——门禁与导出一致：`occt-topo --lib` **1287/1**（红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16**。
- **T-78 余项**：2D（UV）曲线（`BRepToIGES_BRWire.cxx:340-588` 的逐面型 UV 修正）、周期面反周期化与 `periodicU/V`（`GeomToIGES_GeomSurface.cxx:244-343`）、整周椭圆（`GeomConvert_ApproxCurve`，`:620-645`）、写入器"只写可达实体"过滤。

**批 68（交接刷新＋全量基线复核：§14.3"第 82 轮收尾交接"）—— 2026-09-20 第 83 轮**

- **动作（无代码改动）**：把 §14.2（第 63 轮交接）的过期条目覆盖为新的 **§14.3「第 82 轮收尾交接」**：① 抄录本轮实测的全部基线（`occt-core` 290/290、`occt-math` 215/215、`occt-geom` 146/146、`occt-geom2d` 72/72、`occt-topo` 1287/1、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3、`phase3/4/5/6` 4/4·9/9·7/7·5/5、`export_data_obj` 16/16 且 16 模型计数逐位一致）；② 写明 `occt-geom` 151→146、`occt-topo` 1291→1287 的**原因与清单**（批 65/66/67 删除的零消费方模块自带单测：`surface_to_grid` 3、`feature.rs` 4、`comp_curve_to_bspline` 2 等），避免下一轮误判为回归；③ 列出第 74–82 轮新增结项（A26 的 IGES 写侧批 55–63 与交换侧批 64、A15 部分批 65、A8/T-44 批 66、A12/T-48 部分批 67、T-84 结项）；④ **未被阻塞**的余项给出下一步入口（T-78 的 2D 曲线/周期面反周期化/整周椭圆/游离边/可达性过滤、T-77、A16 求交 half、T-67 分步 3 及其解锁的 A15/A19 项）；⑤ **被阻塞**项（T-80→T-81/T-82/T-83 与 T-79、A12 体积覆盖；T-69→A13/A18/A23）；⑥ 做法小结（零消费方自创件直接删的取证要求、反例登记、"提交用显式路径且勿用重定向掩盖 add 失败"——批 65 曾漏提三文件已在批 66 的 `c1448fa` 补）。
- **验证**：本条为文档批，代码树与上一提交一致；上表基线即本轮 `pwsh` 全量复核的输出（5 个 crate `--lib` ＋ 4 道 STEP 门禁 ＋ `phase3/4/5/6` ＋ `export_data_obj`）。

**批 67（A12/T-48：删除零消费方的自创"特征体素布尔＋网格圆柱夹具"模块 `feature.rs`）—— 2026-09-20 第 82 轮**

- **缺口（A12/T-48 原文）**：`brepfeat` 用**解析体积覆盖**、用**网格夹具冒充 `BRepPrimAPI_MakeCylinder`**、以及 `clamp(16, 64)` 的体素分辨率；OCCT 对应件是 `BRepFeat_MakeDPrism`/`BRepFeat_MakeRevol`、`LocOpe_Revol`（T-48 行的出处列了 `brepfeat/features.rs:109,419,361,508`、`feature.rs:197`）。
- **复核**：`crates/occt-topo/src/feature.rs`（`pub mod feature`，283 行）提供 `protrusion/pocket/boss/hole/feature_volume_before_after`，其内部 `cylinder_tool`→`z_cylinder_mesh`（网格圆柱夹具）与 `resolution_for`（`clamp(16,64)`）正是 A12 的两项自创件。`git grep` 逐名核对：这些函数**在仓内零调用**（既无 crate 内消费方，也无 tests/examples 引用），与 `brepfeat` 的同类 API（`brepfeat::boss` 等）重复。
- **落地**：**删除 `feature.rs`** 并移除 `lib.rs` 的 `pub mod feature;`；`brepfeat/features.rs` 里**在用**的 `resolution_for`（`:501`，被 `:334` 的 live 路径调用）与解析体积覆盖（`:107-109`）**保留**——它们属 T-48 的余项，且体积覆盖项**仍被 T-80 阻塞**（第 69 轮复测：换成 `solid_volume` 后 `boss_thru_all_pierces` 仍失败，已回退，不必重复试）。
- **门禁**：`occt-topo --lib` **1287/1**（红 = T-01；由基线 1291 减 4，即被删模块的 4 个单测 `boss_on_box_increases_volume`/`hole_on_box_decreases_volume`/`protrusion_square_sketch_fuses_closed`/`pocket_removes_material`，该模块零消费方 ⇒ 不构成行为回退，其余测试无一削弱）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16**。
- **A12 余项**：`brepfeat` 的解析体积覆盖与在用体素分辨率（`brepfeat/features.rs:107-109/501`）——前置＝T-80（精确布尔闭合穿孔拓扑）；OCCT 侧目标＝`BRepFeat_MakeDPrism`/`LocOpe_Revol`。

**批 66（A8/T-44：删除自创的"采样折线＋强制 1 次"曲线→B 样条转换；OCCT 对应件登记 UNPORTED）—— 2026-09-20 第 81 轮**

- **缺口（A8 原文）**：`convert_bspl.rs` 的 `comp_curve_to_bspline` 用**采样折线拟合**冒充 `GeomConvert_CompCurveToBSplineCurve`——每段按 `n = clamp((b-a)/0.1, 2, 64)` 取点、去重容差 `tol.max(1e-9)`、最后 `resample_bspline(&pts, 1)` **强制 1 次**（多项自创阈值/规则）。
- **复核（决定处置）**：`git grep comp_curve_to_bspline` 显示该函数**只被自身两个单测引用**，生产代码零消费方；T-44 行自身的结论也是"当前函数无生产调用者，前置缺口是 `GeomBSplineCurve::{increase_degree,knots,multiplicities}` 与 `GeomConvert_CurveToBSpline`"。
- **落地**：**删除 `comp_curve_to_bspline` 及其两个单测**，就地留下 UNPORTED 注记并写清 OCCT 控制流：`GeomConvert_CompCurveToBSplineCurve`（`GeomConvert_CompCurveToBSplineCurve.cxx:135-215`）拼接各段的**精确** B 样条像——每段先经 `GeomConvert::CurveToBSplineCurve`（`GeomConvert_CurveToBSplineCurve.cxx` → `Convert_LineToBSplineCurve`/`Convert_CircleToBSplineCurve`/`Convert_EllipseToBSplineCurve`/`Convert_ParabolaToBSplineCurve`/`Convert_HyperbolaToBSplineCurve`，皆有理精确），再做结点/极点缝合与重参数化、插结点。模块头注释同步（去掉它把 `CompCurveToBSplineCurve` 列为"已移植"的说法），并清掉随之失效的 `use`。
- **门禁**：`occt-topo --lib` **1291/1**（红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16**；`occt-geom --lib` **146/146**（较上批 148 少 2，即被删函数的两个单测；无生产消费方，不构成行为回退，其余测试未削弱）。
- **A8/T-44 状态**：自创实现**已清除**；OCCT 侧移植（`GeomConvert::CurveToBSplineCurve` + `GeomConvert_CompCurveToBSplineCurve`，以及其前置 `GeomBSplineCurve::IncreaseDegree/InsertKnots`）**保留为待办**，触发条件＝出现消费方（届时按上述 `.cxx` 逐段移植，不得再采样拟合）。

**批 65（A15 余项之一：删除"非 OCCT 件"的死代码 `surface_to_grid`；`surface_fit::GridSurface` 标 UNPORTED）—— 2026-09-20 第 80 轮**

- **依据（A15 原文）**：`specs/_audit/_index.md:62` 的 A15 未完成项之一即"`surface_fit`/`surface_to_grid` 的 D2（**非 OCCT 件**）"——`occt-geom` 里这两个"曲面"是端口自造类型（OCCT 对应的是 `GeomAPI_PointsToBSplineSurface`/`GeomPlate_BuildPlateSurface` → 真 `Geom_BSplineSurface`，其 `D2` 由 `BSplSLib::D2` 解析给出），它们继承 trait 的 `h=1e-6` 前向差分 `d2`。
- **落地**：
  - **删除 `crates/occt-geom/src/surface_to_grid.rs`**（含 `surface_to_grid`/`UVGrid`/`grid_to_triangles`）：`git grep UVGrid/surface_to_grid` 显示它**只被自身单测引用**，生产代码零消费 ⇒ 无对应 OCCT 类型可移植，忠实处置即删除（自创代码减一）。`lib.rs` 的 `pub mod surface_to_grid;` 同步移除。
  - **`surface_fit::GridSurface` 就地标 UNPORTED**（不删：模块内还有生产消费方使用同一文件的 `fit_plane`，`brep_offset/mod.rs:10`）：注明 OCCT 无此**类型**、其近似算法属 `GeomAPI_PointsToBSplineSurface`/`GeomPlate_BuildPlateSurface`、D2 应为 `BSplSLib::D2`，且当前端口无生产消费方。
  - `surface.rs` 的 trait 默认 `d2` UNPORTED 注释**订正实现者清单**：`GeomSurfaceOfRevolution`/`GeomSurfaceOfLinearExtrusion` **已**覆写（`Geom_RevolutionUtils::CalculateD2`/`Geom_ExtrusionUtils::CalculateD2`，旧注释误列），仍继承差分的只剩 `GeomBezierSurface`/`GeomOffsetSurface`/`surface_fit::GridSurface`。
- **门禁**：`occt-topo --lib` **1291/1**（红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16** 计数一致；`occt-geom --lib` **148/148** —— 由基线 151 减少 3，原因是被删模块自带 3 个单测（`surface_to_grid` 的网格/三角形用例），**其余测试无一被削弱**；被删模块无任何生产消费方，故不构成行为回退。
- **A15 余项（仍开）**：曲面族兜底的**移除**（前置 T-67 步 3）。

**批 64（A26 交换侧余项：STL 读取的节点合并容差改回 OCCT 的"精确相等"；PLY 项复核为已忠实并回退一次错误尝试）—— 2026-09-20 第 79 轮**

- **STL（改）**：`occt-core/src/io/stl.rs::to_triangulation` 原把每个坐标**量化到 `1e-9` 网格**再哈希去重（自创容差）。OCCT 侧：`RWStl_Reader` 在读取时**边读边合并**同坐标节点（`RWStl_Reader.hxx:36` "The nodes with equal coordinates are merged automatically on the fly"），底层 `Poly_MergeNodesTool` 的 `theMergeTolerance` **默认 `0.0`**（`Poly_MergeNodesTool.hxx:50-56`："0.0 by default (only 3D points with exactly matching coordinates are merged)"），`MergeAngle()` 默认 `M_PI/2`（与三角面夹角无关，`RWStl_Reader.hxx:93-95`）。现改为**只有坐标完全相等才合并**（哈希键取浮点位模式，`-0.0` 先归一为 `0.0`，与 OCCT 的 `==` 语义一致；`1e-9` 网格会把 OCCT 保持分离的近重合点并掉，属自创规则）。
- **PLY（复核，负结果已回退）**：A26 记的 PLY 项是"`weld_vertices(1e-9)` + `property list uchar int`（OCCT 每面重复自身 node+偏移）"。复核结论：① 属性类型已在更早批次改为 `property list uchar uint vertex_indices`（`io/ply.rs:99` 注明出处 `RWPly_PlyWriterContext.cxx:214`）；② 焊接已在更早批次从形状导出路径删除（`brep_exchange.rs:134-145` 的订正注释：OCCT 的 `myVertOffset + (tri - anElemLower)`（`RWPly_CafWriter.cxx:230-289`、`RWPly_PlyWriterContext.cxx:266-290`）在"扁平表＝逐 BRep 面节点块顺序拼接"时**就等于全局下标**，故 `export_mesh` 的逐面节点表 + 全局下标与 OCCT **等价**）。**本轮曾尝试**让 `write_ply` 再按**三角面**重复节点（`ply.rs` + `rwmesh/tests.rs`），复核后判定**不是 OCCT 行为**——OCCT 的"每面"单位是 **BRep 面**（一个节点块被该面的多个三角共享，`anElemLower` 即块内下界），三角级重复会把 6 面 24 条记录变成 12×3=36 条 ⇒ 该改动已 `git checkout` 回退（负结果登记，避免提交比基线更差的规则）。
- **验证**：`occt-core --lib` **290/290**（1 ignored）、`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16**。本轮未新增/未削弱任何断言（PLY 尝试的测试改动已随代码一并回退）。

**批 63（A26/T-84：IGES 根结构切到 OCCT 默认的 Faces 模式——144/142 ＋ 402 Group，删除混用的 510/514/186；并修复 phase5 陈旧断言）—— 2026-09-20 第 78 轮**

- **缺口（模式混用）**：OCCT 的 IGES 写侧有**两种模式**，由 `write.iges.brep.mode` 选择（`IGESControl_Writer.cxx:62` 读入，`:123` `myWriteMode ? B1.TransferShape : B0.TransferShape`），该静态量的**默认值是 0**（`IGESData.cxx:90-94`：`eval Faces`=0 / `eval BRep`=1，`SetIVal(...,0)`）：
  - **Faces 模式（默认）**＝`BRepToIGES_BRSolid`/`BRShell`/`BRWire`：面 → **144** TrimmedSurface（边界为 **142** CurveOnSurface），壳 → `IGESBasic_Group`（**402**），实心 → 对壳同规则（`BRSolid.cxx:100-168`、`BRShell.cxx:411-476`）；**不出现** 186/514/510/508/504/502。
  - **BRep 模式**＝`BRepToIGESBRep_Entity`：**186** MSBO → **514** Shell → **510** Face（其 surface 是**基面**）→ **508** Loop → **504** EdgeList / **502** VertexList（`BRepToIGESBRep_Entity.cxx:107-211` 建列表、`:395-528` 建 Loop、`:535-657` 建 Face、`:663-731` 建 Shell、`:738-831` 建 MSBO）。
  端口此前**两种模式混用**：根是 186/514/510（BRep 模式），而 510 的 surface 却指向 144（Faces 模式的实体），且 510 的版面（`510, surface, 曲线…`）对两种模式都不成立（BRep 模式的 510 应为 `510, surface, n_loops, has_outer_loop, loop…`）。⇒ 该文件结构不是 OCCT 任一模式的输出。
- **落地（切到默认的 Faces 模式）**：
  - `emit_shell` ＝ `BRepToIGES_BRShell::TransferShell`（`BRShell.cxx:411-476`）：面的 144 列表，**恰一个 → 直接返回该 144**，多个 → `402, n, 实体…;`（`IGESBasic_ToolGroup::WriteOwnParams`，`IGESBasic_ToolGroup.cxx:104-115`）；`emit_solid` ＝ `BRepToIGES_BRSolid::TransferSolid`（`BRSolid.cxx:154-163`）对壳结果同规则（新增 `group_or_single` 复现该"单元素不包组"规则）。
  - **删除** 510/514/186 的发射（它们是 BRep 模式的实体，Faces 模式不出现）。
  - 面自身仍是批 62 的 `BRepToIGES_BRShell::TransferFace` 结构：144（`outer_boundary_type = !isWholeSurface`，类型 0 时外侧指针写 0）＋ 每 wire 一条 142；`emit_face` 现在返回 **144** 的 DE（原为 510）。
  - 端口为**无 wire 的面**（球面原语）合成的接缝曲线，现在作为该面的 142 轮廓转移（替代 OCCT 原语本来会带的 seam wire），并就地注记。
- **测试订正（原断言编码的是模式 1 的实体与批 56 前的卡片形式）**：`iges::tests::box_iges_contains_expected_entities`（`510/514/186` → `142/144/402`）、`iges::tests::sphere_iges_has_arc_and_solid`（`186` → `144`）、`phase5_integration::real_brepmesh_and_iges`（`186/510` → `144/402`，并把**自批 56 起就失败**的 `iges.starts_with('S')` 改为按 `IGESData_IGESWriter.cxx:763-793` 检查首卡 80 列且**第 73 列**为 `S`）——该断言自批 56 卡片格式合规化后就已陈旧，属本轮顺带修复的既有红灯（`phase5` 现 7/7）。
- **实测（临时探针 `zzprobe_402.rs`，已删；`git grep zzprobe` = 0）实体普查**：`Cube` `{102:6, 108:6, 110:12, 116:8, 142:6, 144:6, 402:1}`，`402,6,16,24,30,36,40,44;`（一个组覆盖 6 张裁剪面）✓；`Sphere` `{100:5, 102:1, 116:3, 123:2, 142:1, 144:1, 196:1}`（单面壳 ⇒ **不包组**，直接是 144）✓；`HoledPlate` `{102:38, 108:32, 110:90, 116:60, 142:38, 144:32, 402:1}` ✓；三例中**均无 510/514/186**。
- **门禁（逐项等于基线；`phase5_integration` 由 6/7 转 7/7）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3/4/6` 4/4、9/9、5/5、`export_data_obj` **16/16** 计数逐位一致。


**批 62（A26/T-78 步 8：IGES 面的 142（CurveOnSurface）＋ 144（TrimmedSurface）按 `BRepToIGES_BRShell::TransferFace` 重建；510/508 结构登记为 T-84）—— 2026-09-20 第 77 轮**

- **缺口（几何级，非版面级）**：端口此前的 144 记录写成 `144, surface, u0, u1, v0, v1, 1, m, 曲线指针…;`——把**曲面参数域**塞进了"外侧边界类型 / 内侧边界条数"两个域，并把**普通曲线**（110/100/126）当成边界引用。OCCT 的 `IGESGeom_ToolTrimmedSurface::WriteOwnParams`（`IGESGeom_ToolTrimmedSurface.cxx:196-218`）是 `144, surface, outer_boundary_type, nb_inner_contours, outer_contour, inner…;`，且边界必须是 **142 CurveOnSurface** 实体（`IGESGeom_ToolCurveOnSurface.cxx:104-115` 写 `142, creation_mode, surface, curve_uv, curve_3d, preference_mode;`）。⇒ 该记录此前无法被 IGES 读入方解释。
- **落地（对照 `BRepToIGES_BRShell::TransferFace`，`BRepToIGES_BRShell.cxx:250-405`）**：
  - **每条 wire 一条 142**：`CreationMode = 0`（`:269`），`Curve3D` ＝ 该 wire 的 3D 曲线实体，`CurveUV = 0`、`PreferenceMode = 2`——这是 OCCT 自己"只有 3D 曲线"那一支的取值（`:285-288` 的 `Iprefer = 2`，与 `:293` 的 `Init`）。**2D（UV）曲线**由 `TransferEdge(edge, face, originMap, length, false)`（`:720`）产生，本端口无该分支 ⇒ 标 UNPORTED。
  - **wire 的 3D 曲线**按 `BRepToIGES_BRWire::TransferWire`（`BRepToIGES_BRWire.cxx:662-781`）：单边 wire 直接用该边曲线实体；**≥2 边发 102 CompositeCurve**（form 0，`IGESGeom_ToolCompositeCurve.cxx` 写 `102, n, curve…;`）。
  - **144**：`outer_boundary_type = !isWholeSurface`，外侧指针在类型为 false 时写 **0**（工具写入器的 `if/else` 分支）；`nb_inner_contours` 与内轮廓指针按 wire 顺序。`isWholeSurface = BRep_Tool::NaturalRestriction(face)`，并在 `:382-388` 对**平面/柱面/锥面**强制为 false（该处的判据是 `CurveOnSurface` 句柄非空，而它在有 wire 时恒非空）——端口用 `BRepTool::natural_restriction` ＋ `classify_surface ∈ {Plane, Cylinder, Cone}` 复现。
  - **外侧 wire 判定**用 `ShapeAlgo::AlgoContainer()->OuterWire(aFace)`（`:275`）＝ `ShapeAnalysis::OuterWire`（`ShapeAnalysis.cxx`：最后一条 wire 直接返回，否则返回第一条 `TotCross2D ≥ 0` 的 wire）——复用端口既有的忠实件 `meshing::model_builder::outer_of_wires`（本轮由私有改为 `pub(crate)`），并按"外侧先、内侧后"（`:276` 先于 `:301`）的顺序发实体。
- **登记余项（新任务 T-84）**：`510`（Face）/`508`（Loop）/`514`（Shell）/`186`（MSBO）仍是端口旧结构，OCCT 的布局与构造在 `BRepToIGESBRep_Entity.cxx`（Loop `:399`/`:594`/`:610`、Face `:537`、Shell `:667`、ManifoldSolid `:742`）＋ `IGESSolid_Tool{Face,Loop,Shell,ManifoldSolid}::WriteOwnParams`：`508, n_edges, {edge_type, edge, list_index, orientation, n_param_curves, [iso, curve]…}…;`、`510, surface, n_loops, has_outer_loop, loop…;`、`514, n_faces, {face, orientation}…;`、`186, shell, orientation_flag, n_void_shells, …;`。另外端口的写入器**无 OCCT 的"可达性"过滤**：OCCT 只写从根被引用到的实体，端口把创建过的实体全部写出（`isWholeSurface = true` 时外侧 142 会成为游离实体，合法但多一个）。
- **实测（临时探针 `zzprobe_142.rs`，已删；`git grep zzprobe` = 0）**：`Cube` 6 面 → 142×6、144×6、102×6（每面 4 边 → `102,4,10,11,12,13;`），144 均为 `flag=1 inner=0`（立方体面是平面 ⇒ 按 `:382-388` 强制非整面）✓；`HoledPlate` 32 面 → 142×38、102×38，144 分布 `flag=1 inner=0`×30、`flag=1 inner=3`×2（两张面带 3 个内轮廓＝3 个孔）✓；`Shape` 11 面 → `flag=0 inner=0`×8（自由曲面且 `NaturalRestriction` 为真 ⇒ 外侧类型 0、外侧指针写 0）、`flag=1 inner=0`×2、`flag=1 inner=1`×1；`Shape-2` 31 面 142×32、102×32 同型分布。样本：`142,0,9,0,14,2;`（Imode 0、曲面 DE 9、CurveUV 0、3D 曲线 DE 14、preference 2）、`102,4,10,11,12,13;`。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16** 计数逐位一致。
- **余项（T-78）**：2D（UV）曲线（`TransferEdge(edge, face, …)`）、不在任何 wire 内的面边（`BRepToIGES_BRShell.cxx:334-365`）、周期面反周期化、椭圆整周支；**T-84**：510/508/514/186 的忠实重建。

**批 61（A26/T-78 步 7：IGES 128 的 `closedU/V` 判据改为 `IsUClosed`/`IsVClosed` 语义＋非周期范围按基面裁剪；周期支标 UNPORTED）—— 2026-09-20 第 76 轮**

- **缺口（三处，逐条对 `.cxx`）**：
  1. **`closedU/V` 是近似**：端口用"首末极点行的**点距** ≤ `Confusion`"判闭合。OCCT 写的是 `mysurface->IsUClosed()/IsVClosed()`（`GeomToIGES_GeomSurface.cxx:347-348`），其实现（`Geom_BSplineSurface_1.cxx:1350-1389`）是：**周期 ⇒ true**；否则取两个 U（或 V）界处的 `UIso`/`VIso` 曲线，用 `Geom_BSplineCurve::IsEqual(other, Precision::Confusion())` 比较（`Geom_BSplineCurve_1.cxx:662-...`：次数、有理标志、结点个数、极点个数必须相等；极点**逐分量** ≤ `Confusion`；结点 ≤ `Precision::Parametric(Confusion)`；重数必须相等；有理时权重 ≤ `Epsilon(w₁)`）。两条 iso 曲线由同一曲面构造（`Geom_BSplineSurface::UIso`/`VIso`），V/U 结点、重数、次数、有理标志、极点个数**全部共用** ⇒ 只剩"极点行"与"（有理时）权重列"需要比较，且比较口径是**逐分量**而非点距，端口的旧判据在"分量全 ≤ Confusion 但点距 > Confusion"时会给错，且缺周期臂与有理权重臂。
  2. **写出的 UV 范围未裁剪**：OCCT 先按曲面 `Bounds()` 收窄（`:245-255`、`:274-284`：`Umin=max(Umin,U0)`、`Umax=min(Umax,U1)`，V 同理），而这个 `start` 是**基面**——裁剪面在 `:492-515` 先递归到 `BasisSurface()` 再进本分支。端口此前直接写面的 UV 界（对基面域外的值不裁剪）。
  3. `Polynom = !(RationU || RationV)`（`:355`）：端口只有一个 `weights` 标志，等价，已在注释写明。
- **落地**：`emit_bspline_surface` 先解析基面（`rectangular_trimmed_basis()`，无损则取自身），极点/结点/次数/周期性/`Bounds` 全部取自基面（与 OCCT 递归一致）；`closed_u/v` 按上述精确口径（`iso_rows_equal` 逐分量 + `weight_equal` 用 `Epsilon(a)`）；非周期支按基面 `u_range()/v_range()` 裁剪写出的四界。新增 `iso_rows_equal`、`weight_equal` 两个自由函数，注释给出 OCCT 出处与"其余比较项共用故可省略"的推导。
- **UNPORTED（就地注记）**：**周期支**（`:256-273`、`:285-302`、`:304-343`）＝`ShapeAnalysis::AdjustToPeriod` 重新取原点、`SetUOrigin`/`SetVOrigin` 旋转结点向量、`SetUNotPeriodic`/`SetVNotPeriodic` → `BSplSLib::Unperiodize` 反周期化后再读结点/极点/闭合标志；OCCT 写出的 `periodicU/V` 取**原曲面**的 `IsUPeriodic`/`IsVPeriodic`（`:235-236`、`:448-449`）。端口无结点旋转/反周期化算子，故保留周期表示并继续写 `periodicU/V = 0`。
- **实测（临时探针 `zzprobe_128.rs`，已删；`git grep zzprobe` = 0）**：对 `Shape`/`Shape-2`/`ATU01038` 的 **75** 条 128 记录做 A/B（`git stash` 前后各跑一次）：**前 9 域（`indU..periodicV`）逐位相同**（判据变更在本语料上不改取值：无周期 B 样条面，也没有"分量全 ≤ Confusion 而点距 > Confusion"的行）；**44 条**的末 4 域（`UMin..VMax`）因裁剪而变——典型 `Shape-2` 的 `128,2,49,2,13,0,0,0,0` 共 12 条，`VMax` 由面界 `46.50039933534364` 收到基面上界 `26.081756329802`。同一二进制两次运行输出逐位一致（排除非确定性）。
- **旁支发现（本轮不改，登记）**：上述 12 条面的 pcurve 推出的 V 界（46.5）**超出其 B 样条基面自身的 V 域**（19.98–26.08）⇒ 端口的 `face_uv_bounds` 与曲面参数域在这类面上不一致（属 A13/A18 一带既有问题，T-69 系列）。IGES 侧按 OCCT 裁剪后记录自洽（写出的结点/极点＝基面，范围＝基面内），但根因在 UV 界，需在网格/UV 侧另行处理。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16** 计数逐位一致。
- **余项（T-78）**：周期面的反周期化（`BSplSLib::Unperiodize` + `SetUOrigin`）与 `periodicU/V` 取值、椭圆整周支（`GeomConvert_ApproxCurve`）、144 边界只取 wire 曲线。

**批 60（A26/T-78 步 6：IGES 实体 104（二次曲线弧）＋ 124（变换矩阵）；DE 卡默认域对齐 OCCT 实体默认值）—— 2026-09-20 第 75 轮**

- **缺口**：`GeomToIGES_GeomCurve` 的**二次曲线分支**（椭圆/双曲/抛物 → 实体 104）未移植 ⇒ 这三类边在 IGES 里退化成**弦线**（110，几何错误）；且 104 需要配套的**变换矩阵实体 124**（凡曲线自身坐标系不是绝对坐标系就要引用它）也未实现。顺带发现 DE 卡另有偏离：线宽写 1、颜色写 7、标签写 `LINE`/`POINT` 一类助记名，而 OCCT 是从实体自身字段取（生成几何全为默认值）。
- **落地（写侧，逐项对 `.cxx`）**：
  - **104 二次曲线弧**（`GeomToIGES_GeomCurve::TransferCurve(Geom_Ellipse)` `GeomToIGES_GeomCurve.cxx:608-700`；`(Geom_Hyperbola)` `:707-773`；`(Geom_Parabola)` `:780-845`；写出布局 `IGESGeom_ToolConicArc::WriteOwnParams` `IGESGeom_ToolConicArc.cxx:103-125`）：`104, A, B, C, D, E, F, ZT, start.x, start.y, end.x, end.y;`。系数按 OCCT 的两级结构：先算 `gp_Elips2d/Hypr2d/Parab2d::Coefficients`（`gp_Elips2d.cxx:25-62` 等）在**恒等二维坐标系** `gp_Ax22d(gp::Origin2d(), gp::DX2d(), gp::DY2d())`（`cxx:671`/`:745`/`:818`）上的值（此时 `gp_Trsf2d T` 为单位阵，一般式塌缩为 `A=1/DMaj, B=±1/DMin, C=D=E=0, F=-1`；退化支 `DMin <= gp::Resolution()` 保留 `A=1, F=-DMaj`），再照调用点的**实参乱序**映射——`E2d/H2d/P2d.Coefficients(A, C, B, D, E, F)`——并保留椭圆/抛物线的 `Init(A, 2B, C, 2D, 2E, F)`（`:678-686`、`:823-831`）与**双曲线不翻倍**的 `Init(A, B, C, D, E, F)`（`:751-759`）差异；半径按 `GetUnit()`（＝1，mm）缩放。`ZT=0`；起止点是 `Build.EvalXYZ` 后的**曲线自身坐标系内**坐标（`IGESConvGeom_GeomBuilder.cxx:212-216`）；形式号由 `IGESGeom_ConicArc::ComputedFormNumber`（`IGESGeom_ConicArc.cxx:97-124`，`eps=1e-8`、`eps4=eps⁴`）按写出的系数算出（1 椭圆 / 2 双曲 / 3 抛物）。椭圆整周支（`:620-645`，走 `GeomConvert_ApproxCurve`）标 UNPORTED。
  - **124 变换矩阵**（`IGESConvGeom_GeomBuilder::MakeTransformation` `IGESConvGeom_GeomBuilder.cxx:218-237`；`IGESGeom_ToolTransformationMatrix::WriteOwnParams` `IGESGeom_ToolTransformationMatrix.cxx:90-104`）：3×4 矩阵 `R11 R12 R13 T1 … R33 T3`，平移列除以 `GetUnit()`，左手系时形式号 1（`thepos.IsNegative()`）。指针写进 **DE 卡 1 第 7 域**（`v[6] = themodel->DNum(anent->DirFieldEntity(7))`，`IGESData_IGESWriter.cxx:324-331`）；实体创建**在 104 之后**（`cxx:688-697` 的 `Conic->InitTransf(TMat)`），DE 顺序与 OCCT 一致。
  - **DE 默认域对齐**：线宽 `theLWeightNum = 0`、颜色 `DefColor() == DefVoid → 0`、标签 `theShortLabel` 为空（`IGESData_IGESEntity.cxx:53-63`；`IGESData_IGESWriter.cxx:276-365`），下标号保留 `0`（`theSubScriptN = 0`，`:380-391` 右对齐写）。批 56 起端口在此三域写的 1/7/助记名不是 OCCT 的输出。
- **实测（临时探针，已删；`git grep zzprobe` = 0）**：
  - `data/ATU01038.step`（含 2 条 STEP `ELLIPSE`）：现写 **104×63 ＋ 124×63**（改前 0，这些边是弦线）；样本 `104 de=804 form=1 trsf=805: 104,3.9939127201348716,0.0,4.0,0.0,0.0,-1.0,0.0,0.2238589149486796,-0.44717272823663046,0.4476496403050508,0.22341610806491702;` 与 `124,0.8944816903663703,0.4471045801591851,-1.469755e-16,-13.05030477214,1.6431359709253164e-16,-9.119804518527102e-24,1.0,56.492117676884,0.4471045801591851,-0.8944816903663703,-7.346537e-17,133.3300083262;`（正交旋转＋平移，形式号 0）。
  - `data/occ/bottom.step`：**104×2 ＋ 124×2**（该文件的 `ELLIPSE(…,1.414213548447,1.)` 对应 `104,0.5000000098472364,0.0,1.0,0.0,0.0,-1.0,0.0,…;` ⇒ `A=1/Rmaj²=0.5`、`C=1/Rmin²=1`，与 STEP 半径一致）。
  - **合成用例（临时单测直调发射器，已删）**：椭圆 10/5 → `104,0.01,0.0,0.04,0.0,0.0,-1.0,0.0,10.0,0.0,5.403023058681398,4.207354924039483;`（形式 1；`u=0` 落在长轴端点）；双曲线 3/2 → `104,0.1111111111111111,0.0,-0.25,0.0,0.0,-1.0,0.0,3.0,0.0,4.629241904445731,2.3504023872876028;`（形式 2）；抛物线焦点距 2 → `104,0.0,0.0,1.0,-8.0,0.0,0.0,0.0,0.0,0.0,0.125,1.0;`（形式 3，`Y²=8X=4·2·X` ✓）；将坐标系平移 (1,2,3) 并绕 z 转 90° 后，104 的 DE 第 7 域指向 124，且起止点回到**局部**坐标（`10.0,0.0`），`124,0.0,-1.0,0.0,1.0,1.0,0.0,0.0,2.0,-0.0,0.0,1.0,3.0;`（列＝局部轴、末列＝原点）。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16** 计数逐位一致（ATU01038 17745/22119、Cube 24/12、Cylinder 146/140、Shape-1 3343/4336、Shape-2 3105/4792、Shape 6150/11372、Sphere 642/1244、Torus 1369/2592、linkrods 3494/5078、rev 104/92、screw 600/790、Cone 195/301、Extrusion 24/12、HoledPlate 180/128、Offset 712/892、OffsetPlaneHoleEdge 48/32）。
- **余项（T-78）**：椭圆整周支（`GeomConvert_ApproxCurve`＋`Reparametrize`，`cxx:620-645`）、128 的 `closedU/V` 精确判据与周期面 `SetNotPeriodic` 等价转换、144 边界只取 wire 曲线。

**批 59（A26/T-78 步 5：IGES 实体 120（回转面）与 122（拉伸柱面）；修正 DE 卡 1 第 2 域＝P 段首行号）—— 2026-09-20 第 74 轮**

- **缺口**：`GeomToIGES_GeomSurface::TransferSurface` 的**扫掠面分支**（`GeomToIGES_GeomSurface.cxx:1000-1025` → `:1032-1104` 拉伸 / `:1111-1188` 回转）未移植 ⇒ 回转面/拉伸面退回 108 平面兜底。`data/Shape.step`、`data/occ/top.step`、`data/occ/bottom.step` 分别有 3/4/1 张回转面，此前全被写成平面。
- **落地（写侧，逐项对 `.cxx`）**：
  - **120 回转面**（`IGESGeom_ToolSurfaceOfRevolution::WriteOwnParams`，`IGESGeom_ToolSurfaceOfRevolution.cxx:119-127`）：`120, axis_line, generatrix, start_angle, end_angle;`。轴线是 `IGESGeom_Line`（实体 110），按 OCCT 的"CAS.CADE 轴取反"写 `Init(Location, Location-Direction)`（`GeomToIGES_GeomSurface.cxx:1178-1183`，`#30 rln` / `#36 BUC60328 face 7`；读侧的同一约定见 `IGESToBRep_TopoSurface.cxx:768`）；母线＝`GC.TransferCurve(BasisCurve, V1, V2)`（`:1155`）；角度＝`Surf->Init(Axis, Generatrix, 2π−U2, 2π−U1)`（`:1185`，弧度原样写，读侧 `:769-770` 同式反算）。**实体创建顺序也照 OCCT**：先建母线实体、再建轴线实体（`:1155` 先于 `:1172`），两者都在 120 之前。
  - **122 拉伸柱面**（`IGESGeom_ToolTabulatedCylinder::WriteOwnParams`，`IGESGeom_ToolTabulatedCylinder.cxx:90-98`）：`122, directrix, end.x, end.y, end.z;`。U 区间**重新取面自身的 `Bounds`**（`GeomToIGES_GeomSurface.cxx:1067-1071`，OCC9490 修正），V 区间保留 OCCT 的无穷处理（`:1058-1065`，`Precision::Infinite()=2e100`）；准线＝基曲线在 `|V1| > Precision::Confusion()` 时按 `Value(U1,V1) − Value(U1,0)`（即 `V1·Direction`）平移后的副本（`:1075-1096`）；终点＝`start->Value(U1, V2)`（`:1078`）。
  - **分派顺序照 `cxx:112-147`**：先 Bounded（B 样条/Bezier/裁剪面），再 Elementary，再 Swept（拉伸先于回转，`:1013-1022`）；裁剪面按 `cxx:492-515` 递归剥到基面再判类型。新增 `swept_surface_kind()` 只经 `rectangular_trimmed_basis()` 递归，故 `Geom_OffsetSurface`（**非** `Geom_SweptSurface`）不会被误判成 122/120。
  - **曲线转移抽成 `emit_curve_range()`**（`GeomToIGES_GeomCurve::TransferCurve`，`cxx:94-126`/`:133-161`/`:481-528` 的分派：110/100/126，二次曲线族 104 仍 UNPORTED 并就地注记），边曲线路径行为不变（无对应分支时仍用弦线替身，注记保留）。
  - **修正（本轮新发现，属批 56 引入的偏离）**：DE 卡 1 第 2 域此前被写成**参数行数**，OCCT 写的是**该实体 P 段首行行号**（`IGESData_IGESWriter.cxx:834` `v[1] = thepnum.Value(i)`），行数在 DE 卡 2 第 4 域（`:835-866` `v[15] = thepnum.Value(i+1) − thepnum.Value(i)`）。现按实体累计 `p_start` 写入，与 P 段自身 `2i−1` 指针（`:903`）自洽（否则按 DE 指针回读 P 段会读错行）。
- **实测（临时探针 `zzprobe_iges_swept.rs` 与临时单测 `zzprobe_122`，均已删；`git grep zzprobe_iges_swept|zzprobe_122` = 0）**：`Shape` 11 面 / 3 回转面 → **120×3**（`120,12,11,0.0,6.283185307179586;`＝整周，轴线 `110,0,0,0,0,0,-1;`＝反向 z 轴，母线 110 与 126 各按类型发）；`data/occ/top.step` 324 面 / 4 回转面 → **120×4**（如 `120,1500,1499,5.13958571289438,6.283185307179586;`，部分回转 2π−U2）；`data/occ/bottom.step` 323 面 / 1 回转面 → **120×1**；`rev.step`/`Extrusion.step` 本就是柱面/平面组合（0 回转 / 0 拉伸面）。122 无现存 STEP 用例，用临时单测直调发射器实测：底面半径 5 圆沿 +z 拉伸、V∈[3,12] → 准线为**平移后**的圆（圆心 z=3）＋ `122,1,5.0,0.0,12.0;`（终点＝`Value(U1=0,V2=12)`）；V∈[0,+∞) → `122,1,5.0,0.0,2e100;`（OCCT 的 `Precision::Infinite()`）。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01 `brepfeat::tests::groove_cuts_cylinder`）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16** 计数逐位一致（ATU01038 17745/22119、Cube 24/12、Cylinder 146/140、Shape-1 3343/4336、Shape-2 3105/4792、Shape 6150/11372、Sphere 642/1244、Torus 1369/2592、linkrods 3494/5078、rev 104/92、screw 600/790、Cone 195/301、Extrusion 24/12、HoledPlate 180/128、Offset 712/892、OffsetPlaneHoleEdge 48/32）。
- **余项（T-78）**：104（椭圆/双曲/抛物二次曲线族，`TransferConic`＋`#124` 变换矩阵）、128 的 `closedU/V` 精确判据与周期面 `SetNotPeriodic` 等价转换、120 的非单位位置 `InitTransf(TMat)`（`#124`）、144 边界只取 wire 曲线。

**批 58（A26/T-78 步 4：IGES 实体 128（B 样条曲面））—— 2026-09-20 第 73 轮**

- **缺口**：`GeomToIGES_GeomSurface::TransferSurface` 的 B 样条分支（`TransferBSplineSurface`）未移植 ⇒ **B 样条面全部退化成 108 平面**（`Shape-2` 的绝大多数面、`Shape` 的部分面）。
- **落地**：按 `IGESGeom_ToolBSplineSurface::WriteOwnParams`（`IGESGeom_ToolBSplineSurface.cxx:64-125`）写 **128**：`128, indU, indV, degU, degV, closedU, closedV, polynomial, periodicU, periodicV, knotU[-degU..indU+1], knotV[-degV..indV+1], weights[0..indU][0..indV], poles[0..indU][0..indV], UMin, UMax, VMin, VMax;`；`indU/indV = nb_poles-1`，结点用**扁平结点向量**（该实体要求的下标区间 `[-deg, ind+1]`），极点按 **U 主序**；非有理面权重写 1、`polynomial=1`（`Geom_BSplineSurface::IsRational`）。数据经四个新 `Surface` 访问器（`Poles`/`UKnots`/`VKnots`/`Weights` 的对应物，等价 `Geom_BSplineSurface` 的取值器）到写侧；`GeomBSplineSurface` 实现之，`Geom_RectangularTrimmedSurface` 转发基面（与既有 `nb_u_poles`/`u_degree` 一致）；**`Geom_OffsetSurface` 刻意不转发**（偏移面的极点不是基面的 ⇒ 保持回退）。
- **实测（临时探针，已删）**：`data/Shape-2.step` 现写 **128×28**（改前 0，全被写成 108 平面）＋ 126×44 ＋ 真平面 108×4 ＋ 144×31；`data/Shape.step` 128×2；`Sphere.step` 仍 196×1。
- **门禁（逐项等于基线）**：`occt-geom --lib` **151/151**、`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致。

**批 57（A26/T-78 步 3：IGES 初等曲面实体 192/194/196/198；此前面型全被写成平面）—— 2026-09-20 第 72 轮**

- **缺口（比登记项更严重）**：`emit_face_surface` 只特殊处理球面，**其余非平面曲面一律退化成 108 平面**（注释自承 "Unclassified curved face: fall back to a plane"）⇒ 柱/锥/环面在 IGES 里被写成**平面**（几何错误），球面则是自造的 120 回转面。
- **落地（按 `GeomToIGES_GeomSurface::TransferSurface`，`cxx:520-600` 的精确类型分派）**：
  - **192 柱面**（`TransferCylindricalSurface`，`cxx:1280-1314`）：`192, location, axis, radius, refdir;`
  - **194 锥面**（`cxx:1318-1362`）：`194, location, axis, ref_radius, semi_angle_deg, refdir;`，含**负半角支路**（`cxx:1344-1350`：参考点关于顶点镜像、角度取反、参考方向反向）
  - **196 球面**（`cxx:1366-1400`）：`196, centre, radius, axis, refdir;`（替换自造 120）
  - **198 环面**（`cxx:1402-1435`）：`198, centre, axis, major, minor, refdir;`
  位置点/轴/参考方向均取自 `gp_Cylinder`/`gp_Cone`/`gp_Sphere`/`gp_Torus` 自身（`cxx:1284-1400`），新增 `emit_direction()` 写实体 **123**（`IGESGeom_ToolDirection::WriteOwnParams`）。
- **实测（临时探针，已删）**：`Sphere → 196×1`、`Torus → 198×1`（`198,2,3,10.0,2.0,4;`）、`Cylinder → 192×1`（`192,3,4,2.0,5;`）＋两端盖 108×2、`Cone → 194×1` ＋端盖 108×1；`Shape-2` 的平面面仍 108，其 B 样条曲面仍 **UNPORTED**（缺 128 发射器）。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致。

**批 56（A26/T-78 步 2：IGES 卡片格式改为合规 + 实体 126（Bezier/B 样条）发射器）—— 2026-09-20 第 71 轮**

- **新发现（比登记项更严重）**：端口 IGES 写侧的**卡片格式本身不合规**——原 `sec_line` 写的是「节字母 + 7 位序号 + 数据」（字母在第 1 列），而 IGES 规定数据占 **1–64**（P 卡）/ **1–72**（S/G/D 卡），序号/DE 指针与节字母在**最后 8 列**（字母第 73 列）。OCCT 出处：`IGESData_IGESWriter.cxx:836-880`（D：72 数据列 + `D` + DE 指针 `2i-1`/`2i`）、`:769-794`（S/G：72 + 字母 + 本行序号）、`:902-925`（P：64 + 空格 + **所属实体的 DE 指针** + `P` + 本行序号）、`:942-943`（T：`S/G/D/P` 计数**空格右对齐** + 40 空格 + `T0000001`）。另外 DE 数据原只有 64 列（第 1 行少了 4 个 2 位状态字段、第 2 行少一个保留字段）⇒ 现按 `cxx:809-880` 补足 72 列。修好后每张卡都是 80 列、字母在第 73 列（探针实测：S/D/P/T 字母与序号列位全对，`cards_with_wrong_length=0`）。
- **实体 126（Bezier/B 样条）**：删除弦线替身，按 `IGESGeom_ToolBSplineCurve::WriteOwnParams`（`IGESGeom_ToolBSplineCurve.cxx:64-105`）写 `126,index,degree,planar,closed,polynomial,periodic, 结点[-degree..index+1], 权重[0..index], 极点[0..index], UMin, UMax, Normal;`；planar/normal 按 `ArePolesPlanar`（`GeomToIGES_GeomCurve.cxx:170-199`，<3 极点用 `GetAnyNormal`）；Bezier 按 `GeomConvert::CurveToBSplineCurve` 的等价形式（单跨、两端结点重数 degree+1）。**UNPORTED**：周期曲线（OCCT 先 `SetNotPeriodic`）与裁剪区间（`Geom_BSplineCurve::Segment`）仍走替身并注明出处。
- **测试适配（唯一一处，语义订正）**：`iges_file_written` 原断言 `content.starts_with('S')`——该断言**编码了不合规的旧格式**；改为校验卡片结构（80 列 + 第 73 列字母）。
- **实测（临时探针，已删）**：`data/Shape-2.step` 现写出 **44 条 126**（改前为 0，样条被写成弦线）、`Shape.step` 12 条；所有卡片 80 列。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致。

**批 55（A26/T-78 步 1：IGES 曲线发射器改按曲线类型分派；解析圆改用精确参数）—— 2026-09-20 第 70 轮**

- **缺口**：`iges.rs` 用「取 6 个参数点的 |C″|，`(max−min)/max < 0.02` 即猜成圆」这条**自创**判据决定发 110/100/其它（audit A26）。OCCT 侧 `GeomToIGES_GeomCurve::TransferCurve`（`GeomToIGES_GeomCurve.cxx:75-116`）按**精确类型**（`IsKind`）分派：`Geom_BoundedCurve`（Bezier/B 样条/裁剪）→ 126/100 族、`Geom_Conic` → 104（圆即 100）、`Geom_OffsetCurve`、`Geom_Line` → 110。
- **落地**：
  1. 新增 `iges_curve_kind()`，用 `Curve` trait 的类型查询（`is_line`/`gp_circ`/`gp_ellipse|gp_hyperbola|gp_parabola`/`bspline_poles|bezier_poles`）做分派；**删除** `edge_curve_kind()` 的采样分类器。
  2. **圆支改用精确参数**：圆心/轴/平面参考点取自 `Geom_Circle` 本身（`GeomToIGES_GeomCurve::TransferCircle`，`cxx:292-329`），不再用三点采样拟合圆心；闭合边按 IGES 闭合圆弧形式（首末点重合）写出。
  3. **未移植部分就地标 UNPORTED**：椭圆/双曲/抛物（OCCT 发 **104**，`TransferConic`）与 Bezier/B 样条（OCCT 发 **126**，`TransferBSplineCurve`）暂无发射器，`TransferCurve` 对其余类型返回空句柄（调用方不写曲线）——这些边暂时保留弦线替身，注记写明 OCCT 函数与实体号（属 T-78 余项）。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致（IGES 写侧不在门禁内，故同时复跑全部既有门禁）。

**批 54（A27/T-64：解析曲线/曲面连续性 3 → `GeomAbs_CN`(6)）—— 2026-09-20 第 69 轮**

- **原登记的阻塞已不复现**：A27 当年把 15 处改成 6 时 `step_to_obj` 变 12/13（Sphere 网格变空），故一直标 `// T-64` 挂着。本轮直接把这 15 处改回 OCCT 值后复跑：**`step_to_obj` 13/13**（Sphere 642/1244 不变）——阻塞是 T-68 的球面 pcurve 工作**之前**测到的，现已消失（消费方 `range_splitter/param_set.rs:290` 的跨度判定不再依赖该值：`GeomAdaptor_Curve::NbIntervals` 对一切非 B 样条、非 offset 曲线都只返回 1 个区间）。
- **落地（15 处）**：`occt-geom/src/{line,circle,ellipse,hyperbola,parabola,plane,cylinder,cone,sphere,torus}.rs` 与 `occt-geom2d/src/{line,circle,ellipse,hyperbola,parabola}.rs` 的 `continuity()` 由 3 改 6，出处：`Geom_Conic.cxx:32-35`、`Geom_Line.cxx:128-131`、`Geom_ElementarySurface.cxx:25`、`Geom2d_Line.cxx:172`、`Geom2d_Conic.cxx:48`；`occt-core/src/kernel/geomabs.rs` 的偏差说明改写为「值已忠实 + 为何当年必须延后」。
- **门禁（逐项等于基线）**：`occt-core --lib` **290/290**、`occt-geom --lib` **151/151**、`occt-geom2d --lib` **72/72**、`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_to_obj` **13/13**、`step_obj_parity` 14/14、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致。

**批 53（A19 余项：UV 栅格替身与它的分辨率就地标 UNPORTED；顺带记录「去 64 上限」实验）—— 2026-09-20 第 68 轮**

- **为什么不能"按 OCCT 对齐阈值"**：OCCT **没有**逐面 UV 栅格剖分器；活管线是 `BRepMesh_FaceDiscret`（边界离散）+ `BRepMesh_Delaun`，其中唯一的"栅格分辨率"是 `BRepMesh_GeomTool::CellsCount`（`BRepMesh_GeomTool.cxx:465-512`，端口已在 `meshing::node_insertion::geom_tool_cells_count` 忠实移植），**且无上限**。端口 `wireframe::face_to_triangles` 是替身，其 `nu/nv` 的自创项（`3` 下界、`64` 上限）没有对应 OCCT 分支。
- **落地（纯注释，零行为改动）**：在该函数文档里写明 ① OCCT 无此件；② `3`/`64` 是自创（下限防退化域无胞、上限约束 `nu*nv` 顶点数，OCCT 的胞数从不喂满一张 UV 栅格）；③ 删除条件指向 `meshing::incremental_mesh::triangulate_model_faces` 的替身注记（该面类走通忠实路径后一并删）。
- **实验记录（去上限）**：把 `clamp(3,64)` 改成 `max(3)` 后 **16 个导出计数无一变化**（ATU01038 17745/22119、Shape 6150/11372、Shape-2 3105/4792 等逐位相同）⇒ 该上限对受门禁模型**不生效**，故按"替身护栏"保留而非删除。
- **门禁**：`occt-topo --lib` **1291/1**（唯一红 = T-01）；其余门禁与基线一致（本轮只改注释，前一批已逐项复跑）。

**批 52（A25/T-61 收尾：Delaun 失败 ⇒ 整面无网格）—— 2026-09-20 第 68 轮**

- **缺口**：`BRepMesh_Delaun::addTriangle` 撞上链接三角形对溢出时，OCCT 抛 `Standard_OutOfRange`（`BRepMesh_PairOfIndex.hxx:41`），异常从 `generateMesh` 逃出、被 `BRepMesh_BaseMeshAlgo::Perform` 空 catch 吞掉（`BRepMesh_BaseMeshAlgo.cxx:40-62`）⇒ `commitSurfaceTriangulation` **根本没执行** ⇒ 该面**完全无三角化**（而不是"已建好的部分网格"）。端口只中止当前多边形，随后照样提交部分网格。
- **落地**：
  1. `delaun/triangulation.rs` 增 `Delaun::failed()`（暴露既有标志）；
  2. `node_insertion.rs::finish_mesh` 在构造完 `Delaun` 后立即判 `failed` 并**提前返回**（在 `erase_free_links`/`use_edge`/`insertNodes`/`optimizeMesh`/`collectTriangles` 之前——正是 OCCT 被异常跳过的那一段），并记入 algo 的 `delaun_failed`；`perform` 把该标志转成该面的 `MeshStatus::FAILURE` + `Err`（与既有 invalid-range 分支同一套模式，正是它让调用方的 `Ok(Err(..))` 支给出「无三角化、不走回退」）；
  3. `fast_discret.rs::discretize_face` 同判据返回空网格。
- **证据边界（诚实记录）**：该分支在**受门禁的 12+4 个模型上不触发**（16 个导出计数与 T0M parity 逐位不变）⇒ 本批验收 = **忠实控制流 + 基线不动**；与 OCCT 的对应关系已在注释里写清。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致。

**批 51（T-81：`BOPTools_AlgoTools2D::AdjustPCurveOnSurf` 忠实化 —— 周期平移而非裁剪）—— 2026-09-20 第 67 轮**

- **缺口**：`AdjustPCurveOnFace` 位置由端口自造的 `pcurve_full::trim_pcurve_to_face`（把 pcurve 裁到面 UV 框）承担；OCCT 该位置**从不裁剪**（`BOPTools_AlgoTools2D.cxx:389-399` 只做整周期平移）。
- **落地**：`algo_tools/queries.rs::AlgoTools2D::adjust_pcurve_on_surf` 按 `BOPTools_AlgoTools2D.cxx:247-400` 逐段重写：① 参考点 = **边参数区间中点**处的 pcurve 值（OCCT `aT = 0.5*(aFirst+aLast)`；`Geom2d_Line` 的声明区间可能是无穷，故由调用方传边参数——`make_pcurves_full` 现传 `BRepTool::edge_parameters`，对应 `BOPTools_AlgoTools::MakePCurve` 按 3D 曲线是否周期选重载，`cxx:1704-1711`）；② `du/dv` = 整周期平移（`cxx:273-344`，含柱面 `dFi = MaxToleranceEdge(face)/R` 支路与 `(VMax-VMin) < aVPeriod` 的 `dv` 取舍）；③ 周期窄于面区间时用 `BRepClass_FaceClassifier` 复核 `(u2+du, v2+dv)`（`cxx:346-387`，端口用 `fclass2d::FClass2d`，构建失败则跳过复核）；④ 结果为同一条曲线平移 `(du,dv)`（`cxx:389-399`）。
- **`make_pcurves_full` 改调忠实件**；`trim_pcurve_to_face` 保留但加 **UNPORTED** 注记（该位置 OCCT 无裁剪分支），三处文档同步。
- **实测（临时探针，已删；`git grep zzprobe` = 0）**：`box[-1,1]³` vs `cyl(r=0.4,z∈[0,2])` 结果不变（FUSE 仍为盒体 + 顶面内圆盘），符合预期——周期平移不改变 T-83 追踪的缝边收尾问题。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致。

**批 50（A15/T-51 余项：`gcpnts` 的 Simpson → `CPnts_AbscissaPoint::Length` + `math_GaussSingleIntegration`）—— 2026-09-20 第 67 轮**

- **缺口（模块自述的偏差）**：`occt-geom/src/gcpnts.rs::curve_length_range` 用自适应 Simpson 积分 `|C′(u)|`，模块头写着「不是翻译，须在 T-51 执行时替换」。
- **落地（逐行对 OCCT）**：
  1. `gauss_order()` = `order(C)`（`CPnts_AbscissaPoint.cxx:57-77`）：`Line` 2、`Parabola` 5、`Bezier` `min(24, 2*Degree)`、`BSpline` `min(24, 2*NbPoles-1)`、其余 10；
  2. `gauss_single()` = `math_GaussSingleIntegration::Perform`（`math_GaussSingleIntegration.cxx:100-150`），阶数按 `math::GaussPointsMax() = 61` 截断（`math.cxx:24-27`）。OCCT 用表（`math.cxx` 的 `Point[]`/`Weight[]`）存节点与权重，端口用已有的 `occt_math::gauss::gauss_legendre` 算同样一组值（已在注释里写明累加次序不同）；
  3. `gauss_single_tol()` = 带容差构造（`cxx:64-98`）：最多 `IterMax = 13` 轮，每轮把等分子区间数翻倍，直到相邻两次总和之差 ≤ 容差；
  4. `curve_length_range` = `CPnts_AbscissaPoint::Length(C, U1, U2, Tol)`（`cxx:168-184`）：`|积分值|`；
  5. **删除** `adaptive_length`/`MAX_DEPTH` 这套 Simpson 替代件（替换而非并存）。
- **就地标 UNPORTED（未新增行为）**：`abscissa_point`/`uniform_abscissa` 的反解仍是端口自带的 safeguarded Newton/二分，OCCT 用 `math_FunctionRoot` 驱动 `CPnts_MyGaussFunction`（`cxx:395-432`），缺 `math_FunctionRoot` 端口；本轮只把越界判据改成与积分容差同量级，对应 OCCT **加宽搜索区间**（`cxx:367-369`：`myUMin = U1-DU`、`myUMax = U2+DU`）——恰好等于曲线全长时落回端点参数，真正超长的请求仍报错（既有两条断言都不放宽）。
- **门禁（逐项等于基线）**：`occt-geom --lib` **151/151**、`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致。

**第 71 轮诊断（T-82/T-83 的取证链，无代码改动）—— 2026-09-20 第 66 轮**

- **手段**：临时插桩（`ZZ_DBG_SOL`：`bop_fill_in3d` 的 per-solid `a_l_in`/`a_l_internal`、`bop_build_common/draft_solid.rs::build_draft_solid` 的面→images、`bop_build_solids::collect_solid_split_faces` 的面→images）＋临时探针 `zzprobe_sol.rs`；全部已删（`git grep ZZ_DBG`/`zzprobe` = 0），工作树回到 `03ada27`。
- **取证链（`box[-1,1]³` vs `cyl(r=0.4,z∈[0,2])`，FUSE）**：
  1. `build_split_faces_occt`：盒体顶面 = 2 个 area ✅（外环带孔 2 wires + 内圆盘 1 wire）；**柱面 = 1 个 area / 3 wires**（假面）⇒ 周期面未被节线切成两环带。
  2. `bop_fill_in3d`：盒体 draft 的 IN = [柱体底盖圆盘]；柱体 draft 的 IN = [z=1 截面圆盘]；`INTERNAL` 均为空 ⇒ **分类本身合理**（不是「整根柱体被判 IN」）。
  3. GF 结果（`b.result()`）= `Solid/8f`（盒体：7 `Plane/Forward` + 1 `Plane/Internal` ＝柱体底盖）＋ `Solid/4f`（柱体：`Cylinder`（z∈[0,2]，**未切**）＋两端盖＋z=1 圆盘，**四面全 `Internal`**）。
  4. 因此 `BOPAlgo_BOP::BuildSolid` 的 `MapFacesToBuildSolids` 跳过全 `Internal` 的柱体四面（`mfs=7`/`sfs=7`），`BuilderSolid` 只能重建盒体 ⇒ 最终 FUSE = 盒体 + 顶面内圆盘（7 平面面）。
  5. 柱体四面为何全 `Internal`：`builder_solid::connect_faces_into_shells(force_internal=true)` 的**内部壳**路径 —— 即柱体 `SplitSolid` 把整组面（未切的侧面＋两端盖＋内插圆盘）当成了内部壳 ⇒ **是「周期面未切开」的下游后果**，不是独立的分类缺陷。
- **订正两处前述推测（留痕）**：
  1. **`build_solid` 的「按实参播种 aMFS」是错的**：OCCT `BOPAlgo_BOP.cxx:1156` 在 RC 循环前 `aMFS.Clear()`，只用实参播种求 `aMTSols`（共享面的实体），随后按 `myRC` 重建 —— 端口现有结构**忠实**。本轮实验（按实参播种 → `mfs=11`/`sfs=11`）已 `git checkout --` 回退。
  2. **收尾判定 `anIsSameV2d` 不是终端阻塞**：临时禁用该 UV 判定（`ZZ_DBG_NOUV`）后，柱面从「1 area/3 wires」变为「1 area/2 wires」，但 FUSE 结果**完全不变** ⇒ 即使收尾通过，`perform_areas`（`IsGrowthWire` / `IntTools_FClass2d::IsHole`）仍把一条环带 loop 当成 hole（1 growth + 1 hole ⇒ 一个带内环的面）。
- **`perform_areas` 的逐 loop 实测（`ZZ_DBG_AREA`，已删）**：盒体顶面 = `loop0`(4 边, growth, ring_area=+4.000000) + `loop1`(1 边, growth, +0.502125) + `loop2`(1 边, hole, −0.502125) ✅；**柱面 = 6 个单边 loop**（4 个闭合圆 + 缝的两段），其中 5 个 `is_hole=true`、1 个 `growth=true`（`IsGrowthWire` 命中前一个 hole 的边），且**所有 ring_area = NaN**（`FClass2d` 的 `outer_ring()` 为空：单边闭合圆的 UV 采样环退化成零面积直线 ⇒ 被 `area.abs() < SQUARE_CONFUSION` 丢弃 ⇒ 默认「hole」）⇒ 端口的假面「1 area / 3 wires」= 1 growth + 2 hole 拼出来的。
- **下一轮入口（T-83，按此顺序）**：① 收尾判定：按 OCCT 的 `aVertMap`（TShape+Location+**Orientation**）与逐边实例 `IsClosed` 建模，使两次 band walk 能收尾（判据：柱面 `build_split_faces_occt` 出 **2 个 area**）；② `perform_areas` 的 growth/hole 判据：按 `IntTools_FClass2d.cxx:548-565` 核对端口 `fclass2d/classifier.rs:250`（逐 wire 有符号 UV 面积；`myIsHole` 取最后一个 wire）与 UV 多边形符号（判据：两条环带 loop 都判为 growth）；③ 之后复跑 `box∪cyl` 三操作，再看 T-82/T-81。

**批 49（T-80 第四段：柱面侧面 wire 对齐 `BRepPrim_OneAxis::LateralWire`；登记 T-83）—— 2026-09-20 第 70 轮**

- **读出的缺口**：端口柱面侧面的 wire 是 `[底圈(rev), 缝(fwd), 顶圈(fwd), 缝(rev)]`；OCCT `BRepPrim_OneAxis::LateralWire`（`BRepPrim_OneAxis.cxx:660-684`）为 `AddWireEdge(TopEdge(), false)` → `AddWireEdge(EndEdge(), true)` → `AddWireEdge(BottomEdge(), true)` → `AddWireEdge(StartEdge(), false)`，即 `[顶圈(fwd), 缝(rev), 底圈(rev), 缝(fwd)]`。旧注释担心的"U 绕数"两种排法都为零（顶 +2π、底 −2π、缝 +h/−h），实测 `Cylinder.step` 导出计数不变 ⇒ 该担心不成立。
- **为什么影响 BOP**：`BOPAlgo_Builder::BuildSplitFaces` 按面的 wire 顺序填 `aLE`（`BOPAlgo_Builder_2.cxx:362-368`），`aLE` 次序决定 `BOPAlgo_WireSplitter::SplitBlock` 里 `mySmartMap` 的顶点槽位与各 `Path` 的**入口边**。旧顺序下 `Path` 从底圈进入，band 走到缝的 u=0 侧而起点在 u=2π ⇒ 收尾判定 `anIsSameV2d`（`BOPAlgo_WireSplitter_1.cxx:447-467`，`aD2 < aTol2D2`）不成立 ⇒ 该次 `Path` 一条 wire 都没产出；`SplitBlock` 随后报 `Path produced no wires`，调用方退回坐标链式拼接 ⇒ 柱面得到 **6 个单边 loop**、`FaceBuilder` 交出 **1 个 3-wire 的假面**（面积/边界全错）。
- **落地**：`primitives.rs::make_cylinder` 按 OCCT 顺序重排 wire，并改写注释记录出处与本轮结论。
- **实测（临时探针，已删；`git grep zzprobe`/`ZZ_DBG` = 0）**：`box[-1,1]³` vs `cyl(r=0.4,z∈[0,2])`
  - 对齐前：`Path` 第 1 步就从底圈进入，4 步后落在 u=0 侧，收尾失败（无 loop）；
  - 对齐后：两次 walk 分别走完 `顶圈(0,2)→(2π,2) → 缝rev(2π,2)→(2π,1) → 节线圆(2π,1)→(0,1) → 缝fwd` 与 `缝fwd(2π,0)→(2π,1) → 节线圆rev → 缝rev → 底圈`，**两条环带都被正确遍历**；
  - 仍失败于收尾：起点 UV = 2π（顶圈/底圈侧）而终点 UV = 0（缝的另一侧），`anIsSameV2d` 恒假 ⇒ 依旧 `Path produced no wires` ⇒ 柱面 6 单边 loop / 1 个 3-wire 假面；`FUSE` 仍为「盒体 + 顶面内圆盘」7 平面面。
- **根因细化（本轮 dump，写入 T-83）**：OCCT 的 `aVertMap` 以（TShape, Location, **Orientation**）为键、`bIsClosed` 逐边实例判定（柱面只有 `ESTART` 带闭合形式两 pcurve，`EEND` 没有，`BRepPrim_OneAxis.cxx:422-439`），并对闭合边用 `SetParameters(ETOP, TopEndVertex(), 0, myAngle)`（`cxx:404-408`）把顶点参数写进点表示；端口三者都未建模（`vert_closed` 按 TShape 合并、`curve_on_surface` 只按朝向读第二条、无点表示）⇒ 收尾判定无法成立。
- **RC 侧同时取证（写入 T-82）**：`b.result()` 为 `Solid/8f`（盒体：7 个 `Plane/Forward` + 1 个 `Plane/Internal` 顶面内圆盘）与 `Solid/4f`（柱体：`Cylinder`、两个端盖、一个 z=1 圆盘，**四个面全为 Internal**）⇒ 柱面外侧的上环带与顶盖被状态分类误判为 Internal，`map_faces_to_build_solids` 全部跳过（`mfs=7/sfs=7`），`BuilderSolid` 只能重建出盒体。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`export_data_obj` **16/16** 计数逐位一致（含 `Cylinder.step` 146/140）。

**批 48（T-80 第三段：柱面缝边的两个 pcurve + `CurveOnSurface` 的 `PCurve2` 规则；登记 T-82）—— 2026-09-20 第 69 轮**

- **缺口 1**：`BRepPrimCylinder::make_cylinder` 对三个面**一个 pcurve 都没登记** ⇒ `BRep_Tool::IsClosed(seam, lateral)`（= 该边在 `Geom_Surface` 上有一个 `IsCurveOnClosedSurface()` 的表示，`BRep_Tool.cxx:820-840`）恒假 ⇒ `BOPAlgo_Builder::BuildSplitFaces` 的缝边支（`BOPAlgo_Builder_2.cxx:429-446`：`bounding_edge_is_closed_seam` → `DoSplitSEAMOnFace`）从不触发。对照 `BRepPrim_OneAxis::LateralFace`（`cxx:388-439`；`myVMin=0`、`myVMax=height`、`myAngle=2π`、`myMeridianOffset=0`、整周 `HasSides()=false`）：顶/底圆 pcurve = `gp_Lin2d(gp_Pnt2d(0, VMax|VMin), +X)`；缝边用**闭合形式** `SetPCurve(E, F, c1, c2)`，`c1 = gp_Lin2d((myAngle, -offset), +Y)`、`c2 = gp_Lin2d((0, -offset), +Y)`（`cxx:434-439`）。已在 `primitives.rs` 按此登记（圆面 `set_edge_pcurve`+`set_pcurve_range`，缝边 `set_edge_pcurves` 两条），与球面已有的 `BRepPrim_OneAxis` 做法一致。
- **缺口 2**：`boptools_2d::curve_on_surface`（= `BRep_Tool::CurveOnSurface`）原来无视边朝向返回 `edge_pcurves().next()`；OCCT 在"表示为闭合曲面表示且 `E.Orientation()==REVERSED`"时返回 **`PCurve2()`**（`BRep_Tool.cxx:347-357`，逐行核对）。已补该规则（存在第二条 pcurve 时 REVERSED 取第二条）。
- **实测（临时探针，已删；`git grep zzprobe` = 0）**：`box[-1,1]³` vs `cyl(r=0.4,z∈[0,2])`
  - 批 47 后：`FUSE = 2 实体 / 10 面`（柱面侧面 + 两端盖 + 顶面内圆盘 + 盒体 6 面）；
  - 批 48 后：`FUSE = Solid 1 个 / 7 面全平面 / vol 9.3333 / shells=1 / warn=[]` —— 节线 edge 与盒体顶面圆孔仍在，但**柱面自身的面被状态分类全部丢弃**，而盒体顶面的**内圆盘**（本应属内部、不进 FUSE 外壳）被保留 ⇒ **下一断点 = 曲面面的状态分类/保留**，已立项 **T-82**。（`vol 9.3333 > 8` 亦不可信：该壳并不闭合，与 `CUT` 的告警同源。）
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致。

**批 46（T-80 根因第一段：`IntTools_FaceFace` 漏掉 pcurves 的「换回」；并删掉自创的重建检查）—— 2026-09-20 第 69 轮**

- **对着 `.cxx` 读出的缺口**：`IntTools_FaceFace::Perform` 先用 `SortTypes` 把两个面按解析类型排序（`cxx:351-357`），交线算法因此按**排序后**的面给出 pcurves；但 OCCT 随后把它们**换回来**——平面/平面早返回支 `cxx:420-436`、通用支 `cxx:550-563`（交点集同 `:595-604`）。⇒ 返回时 `FirstCurve2d()` 属于调用者传入的 `aF1`，正是 `BOPAlgo_PaveFiller::MakeBlocks` 绑定的 `myDS->Shape(nF1)`（`BOPAlgo_PaveFiller_6.cxx:747-748`、`:914`）。端口只做了排序、没做换回 ⇒ 平面拿到柱面的 pcurve。
- **落地**：`int_face_face.rs::perform` 在 `intersect_surfaces` 之后按 `reverse` 交换 `pcurve1/pcurve2`（端口 `FaceFaceCurve` 只带 pcurves，无交点集）；文档注释改为记录 OCCT 的真实约定。同批删除 `pave_ff_make.rs::valid_block_point_for_face` 里端口自创的「pcurve 值必须能重建 3D 中点」检查——OCCT `IntTools_Context.cxx:738-752` 只要 pcurve 非空就走 2D 分支，空才回落到 `IsValidPointForFace`。
- **探针实测（临时，已删；`git grep zzprobe`/`ZZ_DBG` = 0）**：`box[-1,1]³` vs `cyl(r=0.4,z∈[0,2])`
  - 改前：`pc1 on f1=2.086 / on f2=2.8e-17`、`pc2 on f1=6.2e-17 / on f2=0.709` ⇒ 平面面 `reconstructed=false` → 3D 分支 false ⇒ 节线块被丢；`FUSE=6 面全平面 vol 8.0000`、`CUT` 不变、`COMMON` 空。
  - 改后：`pc1 on f1=6.2e-17`、`pc2 on f2=2.8e-17`，两面 2D 分支都 true ⇒ 节线 edge 建出、盒体顶面 `Plane[4 外 + 1 内]`（圆孔 ✅）、`FUSE=Solid 7 面（含 Cylinder）closed chi=2 vol 8.2132`；`CUT=8 面`；`COMMON` 仍空。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16**（顶点/面数逐位一致）。

**批 47（T-80 根因第二段：`BRep_Tool::Parameter(V,E,S,L)` 移植；登记 T-81）—— 2026-09-20 第 69 轮**

- **缺口**：`BOPAlgo_WireSplitter::SplitBlock` 里 `Coord2d`/`Angle2D` 都调 `BRep_Tool::Parameter(V,E,S,L)`（`BOPAlgo_WireSplitter_1.cxx:315`、`:412`、`:417`、`BRep_Tool.cxx:301-352`），端口用的是自写的 `vertex_parameter`：按**形状同一性**取 first/last（`stored_vertices(e).first().same_tshape(v)`）⇒ 闭合边（环）两端是**同一个顶点形状**，两端都取到起始 UV。
- **实测后果（探针轨迹）**：柱面 walk 一步就"回到起点" ⇒ `Path` 在 i=0 处闭环成**单边 loop**（4 条闭合边各成一个 loop），柱面 UV 矩形从未走通；柱面最终只剩「以节线圆为唯一 wire」的畸形面。
- **落地**：`wire_splitter_block.rs` 新增 `vertex_parameter_on_face`（逐行对 `BRep_Tool.cxx:301-352`：在**存储的**顶点方向上找 `VF`；第二次同形出现时记 `rev = (E.Orientation()==REVERSED)`，仅当朝向与 `v` 相同才替换；再按 `FORWARD→first`/`REVERSED→last`（各受 `rev` 翻转）；`INTERNAL`/未命中支回落到 3D `BRep_Tool::Parameter(V,E)` 并注明 `BRep_TVertex::Points` 未建模），`coord2d` 与 `angle_2d` 改用之；把 `(FORWARD, REVERSED)` 归一化抽成 `stored_vertices_oriented`，`oriented_vertices` 在其上再合成边朝向。
- **实测（同探针）**：柱面 walk 变为 `(2π,0)→(0,0)→(0,1)→(2π,1)→(0,1)→(0,2)→(2π,2)`（正确遍历 UV 矩形），柱面侧面 + 两端盖保留；`FUSE` 变成 2 实体 / 10 面（**仍不正确**）。**新的断点已定位**：环带在缝处闭不上——`SplitBlock` 的 `bIsClosed` 候选过滤（`cxx:470-500`）比较 `Coord2dVf(候选)` 与当前 `aPb`，走到 u=2π 时缝边候选仍是 (0,·)，2π > tol 被拒 ⇒ 只剩"同形反向副本"可选（`cxx` 里角为 `TwoPI` 的那支）⇒ 记 `priz=false` 不闭环。⇒ 立项 **T-81**：`BOPTools_AlgoTools2D::AdjustPCurveOnFace`（`:209-243`→`:247-400`）只做**周期平移**（`GeomInt::AdjustPeriodic` + 圆柱面 `dFi=tol/R` 支路 + `Translate`），端口在 `pave_blocks/make_blocks.rs:641` 换成了 `trim_pcurve_to_face`（裁剪，OCCT 无此动作）⇒ 周期面上缝边缺 u=2π 的第二表示。
- **门禁（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3（红 = T-05）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`export_data_obj` **16/16** 计数逐位一致。

**批 45（T-80 决定性测量：pcurves 与 DS 面索引错位 + 解析臂 pcurve 的无穷声明区间）—— 2026-09-20 第 68 轮**

- **探针（临时 `examples/zz_probe_pc2.rs`，已删，`git grep zz_probe`=0）**：跑 `PaveFiller` over `[box, cyl]`，对每个 `interf_ff` 的每条 `BopdsCurve` dump `range`、3D 曲线区间、`pcurve1/pcurve2` 的**声明区间**与在 `r0/mid/r1` 处的 `d0` 值。
- **box∪cyl（`ff=0 pair=(11,36) curve=0`）实测**：
  - `range=(0.000000, 6.283185)`、`c3d=(0.0, 6.283185)` ✅ 同参；
  - **`pc1`：声明区间 `(-inf, inf)`**，`d0(π)=(3.141593, 1.000000)`、`d0(r0)=(0,1)`、`d0(r1)=(6.283185,1)` ⇒ 这组 UV 是**柱面**的 (angle, z)（v≡1）——取值正确，但**声明区间是无穷**；
  - **`pc2`**：声明 `(0, 2π)`，`d0(π)=(0.600000, 1.000000)`、两端点同为 `(1.400000, 1.000000)`（固定点/退化）。
- **两个结论**：① 解析臂产出的 pcurve 带**无穷声明区间**（`make_pcurve_full` 的缺陷）⇒ 上一轮"按声明区间判定是否重参数化"的修法② 在无穷区间上失效（这正是它无效的直接原因）；② **pcurves 与 DS 面索引错位**：`pc1` 的取值属于柱面，而被先判定、报 `uv=(2.715211,1.000000)` 的失败面是**平面**。`int_face_face.rs:248-249` 写明"结果曲线一律以**排序后**的 `face1_idx=0`/`face2_idx=1` 携带 pcurves"，而 `pave_ff_perform.rs:464-466` 却按**未排序**的 DS 索引 `(n_f1, n_f2)` 存入 `BopdsInterfFf` ⇒ 平面拿到柱面的 pcurve ⇒ `IsValidBlockForFaces` 的重建检查失败 ⇒ 节线块被丢弃 ⇒ 无节线 edge ⇒ 面不分割 ⇒ GF 结果 = 原参数 ⇒ 只剩盒体。
- **下一轮入口（已写入 T-80 行）**：按 OCCT `IntTools_FaceFace::Face1()/Face2()`（排序后的面）与 `BOPAlgo_PaveFiller` 的用法，在 `pave_ff_perform` 里依据 `ff.face1()/ff.face2()` 与 `(n_f1,n_f2)` 的对应关系（必要时交换 `pcurve1/pcurve2`）修正；顺带修 `make_pcurve_full` 产出的**无穷声明区间**。判据：`box∪cyl` 的 `reconstructed` 变 true → 建出节线 edge → `b.result()` 出现被切面、体积 ≈8.5。
- **本批无行为代码改动**（探针已删，工作树干净）；门禁与基线一致（`occt-topo --lib` 1291/1、四道 STEP 门禁、`phase3/4/6`、`export_data_obj` 16/16）。

**批 44（T-80：确认分支 = 解析臂；问题收敛到 pcurve 的"声明区间 vs 实际取值"）—— 2026-09-20 第 67 轮**

- **读码求证（本批只读 + 记录）**：① `int_face_face.rs::intersect_surfaces` 的 `(Plane, Cylinder) | (Cylinder, Plane) => self.plane_cylinder(...)`（`:349-352`）⇒ 该对走**解析臂**；② `plane_cylinder → curve_from_ic → finish_curve`（`pcurve1/2 = pcurve_of_curve(&curve, range, fa/fb)`）；③ `pave_ff_perform.rs:464-470` **原样**把 `c.pcurve1/c.pcurve2` 存进 `BopdsCurve`（`nc.set_pcurves(...)`）⇒ 上一轮的修法②（在 `pcurve_of_curve` 里重参数化）**理应在该对上生效**却无效果。
- **唯一自洽解释**：`make_pcurve_full` 返回的 pcurve **声明的 `first/last` 就等于 3D 的 `[range.first, range.last]`**（故包裹条件 `|f-range.first|<1e-12 && |l-range.last|<1e-12` 命中、直接返回未包裹），但其**取值**是按索引/弧长参数化的。再结合第 64 轮插桩的 `uv=(2.715211, 1.000000)`：`u=2.715` 对**平面面**（盒体顶面，UV 即 x,y ∈ [−1,1]）**越界** ⇒ 失败面就是盒体顶面，其 pcurve 的 u 不是该点的 x ⇒ `surf.d0(u,v)` 距 3D 中点 1.7 左右 ⇒ `reconstructed=false`。
- **下一轮入口（已写入 T-80 行）**：进 `make_pcurve_full`（`ShapeConstruct_ProjectCurveOnSurface` 移植体）看它为"平面上的圆"构造的 2D 曲线——若用"采样点 + 索引结点"建 2D B 样条，需按 OCCT 把 2D 曲线与 3D 曲线**同参**；或在 `pcurve_of_curve` 里改按**取值**判定（例如比较 `pc.d0(range.first)`/`pc.d0(range.last)` 与 `curve.d0(...)` 重建误差）而不是按声明区间决定是否重参数化。
- **本批无行为代码改动**：门禁与基线一致（`occt-topo --lib` 1291/1、四道 STEP 门禁、`phase3/4/6`、`export_data_obj` 16/16）；工作树干净。

**批 43（T-80：修法② 试做并回退；下一步转向"该对的 BopdsCurve 由哪条分支产出"）—— 2026-09-20 第 66 轮**

- **试做**：在 `int_face_face_helpers.rs::pcurve_of_curve` 里加一个 2D 线性重参数化适配器（`ReparamCurve2d`，镜像 3D 的 `occt_geom::curve_reparam::ReparamCurve`：d0/d1/d2 的链式缩放、first/last/period/continuity/is_line/gp_circ2d 转发、`transform`/`reverse` clone-on-write），把通用投影 `make_pcurve_full` 的结果映到 3D 的 `[range.first, range.last]`，意图让 `pc.D0(t)` 落在 `C(t)`（`IntTools_Curve::SetCurves` 语义）。
- **结果：无效**（判据探针，临时 `examples/zz_probe_pc.rs`，已删）：`box∪cyl` 仍 `FUSE=6 平面面/vol 8.0000`、`CUT=6 平面面/vol 8.0000`、`COMMON=空(0 面)`；平面参考 `FUSE box∪box(重叠)=14 面/vol 10.0000` 不变。⇒ **该对的 pcurve 并非（或不只是）由 `finish_curve`/`pcurve_of_curve` 提供**。
- **处置**：按"失配即停"`git checkout -- crates/occt-topo/src/int_face_face_helpers.rs` **回退**，复跑 `occt-topo --lib` **1291/1** 确认恢复基线；探针已删（`git grep zz_probe` = 0），工作树干净。
- **下一步（已写入 T-80 行）**：先确认该对（`ff=0 pair=(11,36)`）的 `BopdsCurve` 由哪条分支产出——解析臂 `int_face_face_analytic.rs::plane_cylinder → curve_from_ic → finish_curve`，还是通用 tracer（`int_face_face_make_curve.rs::make_bspline`/`make_bspline2d`/`wline_to_curve`）；再 dump 其 `pcurve1/pcurve2` 的参数区间与 3D 曲线 `[first,last]` 的差异。若来自 tracer，要修的是 **WLine 的 2D 点参数化**（`make_bspline2d` 与 `make_bspline` 的结点来源是否同参）。
- **本批无行为代码改动**（改动静止于试做并回退）：门禁与基线一致（`occt-topo --lib` 1291/1、四道 STEP 门禁、`phase3/4/6`、`export_data_obj` 16/16）。

**批 42（T-80 追到 pcurve 生成处：解析臂的 pcurves 走通用投影，与 3D 曲线不同参）—— 2026-09-20 第 65 轮**

- **对读链路（本批只读 + 记录，无行为改动）**：`int_face_face_analytic.rs::plane_cylinder`（`IntAna_QuadQuadGeo` 的精确解）→ `int_face_face.rs::curve_from_ic`（`crate::int_face_face_analytic` 的解析臂）→ **`finish_curve`**：`pcurve1 = pcurve_of_curve(&curve, range, fa)`、`pcurve2 = pcurve_of_curve(&curve, range, fb)`；而 `int_face_face_helpers.rs::pcurve_of_curve`（`:277`）把 3D 曲线**临时做成 edge**（参数 `[range.first, range.last]`）后调 **通用投影 `make_pcurve_full`**（`ShapeConstruct_ProjectCurveOnSurface` 式近似）⇒ 得到的 2D 曲线参数系是"采样/裁剪"参数，**与 3D 曲线（解析圆，0..2π）不同参**；这正好解释第 64 轮插桩看到的 `pc.d0(t_mid=1.0) → uv=(2.715211, 1.000000) reconstructed=false`。
- **OCCT 侧对照**：`IntTools_FaceFace` 的解析分支由 `IntPatch` 直接产出**与 3D 曲线同参**的 pcurves（`IntTools_Curve::SetCurves` 只做装载），因此 `IsValidBlockForFaces` 的 `pc.D0(aMidPar)` 必然落在 3D 中点，节线块有效、节线 edge 建出、面被切开。
- ****第 66 轮：修法② 试做 → 无效，已回退**（`int_face_face_helpers.rs::pcurve_of_curve` 里加 2D 线性重参数化适配器，把 `make_pcurve_full` 的结果映到 `[range.first, range.last]`）：`box∪cyl` 三操作仍为 `FUSE=6 平面面/vol 8.0`、`CUT=6 平面面/vol 8.0`、`COMMON=空`（平面参考 `FUSE=14 面/vol 10.0` 不变）⇒ **该对的 pcurve 并非由 `finish_curve`/`pcurve_of_curve` 提供**（或不只是它），改动已 `git checkout --` 回退（无残留）。**第 67 轮：分支已定，问题收敛到 pcurve 的"声明区间 vs 实际取值"**：读码确认 `(Plane, Cylinder)` 走**解析臂**（`int_face_face.rs::intersect_surfaces` 的 `(Plane,Cylinder) => self.plane_cylinder`，`:349-352`），其 pcurves 由 `finish_curve → pcurve_of_curve` 产出，并在 `pave_ff_perform.rs:464-466` **原样**存进 `BopdsCurve`（`nc.set_pcurves(c.pcurve1.clone(), c.pcurve2.clone())`）⇒ 修法② 理应生效却无效果，唯一自洽的解释是：**该 pcurve 声明的 `first/last` 就等于 3D 的 `[range.first, range.last]`（故包裹条件未触发），但其取值是按索引/弧长参数化的**。再结合上一轮插桩的 `uv=(2.715211, 1.000000)`：u=2.715 对**平面面**（盒体顶面，UV 即 x,y ∈ [-1,1]）来说**越界** ⇒ 失败面就是盒体顶面，其 pcurve 的 u 不是点的 x。**第 68 轮决定性测量（探针，已删）**：对 `ff=0 pair=(11,36) curve=0` 直接 dump `BopdsCurve` 的几何：`range=(0.000000,6.283185)`、`c3d=(0.0,6.283185)`；**`pc1` 声明区间 `(-inf, inf)`**、`d0(π)=(3.141593,1.000000)`、`uv_at_r0=(0,1)`、`uv_at_r1=(6.283185,1)`（这组 UV 正是**柱面**的 (angle, z)，v≡1）；`pc2` 声明 `(0,2π)`、`d0(π)=(0.600000,1.000000)`、两端点都是 `(1.400000,1.000000)`（退化/固定点）。⇒ 两个结论：① **`pc1` 的声明区间是无穷**（`make_pcurve_full` 对解析臂给出的 pcurve 声明有缺陷；我上一轮按声明区间判定的重参数化因此在无穷区间上失效——这正是"修法②无效"的原因）；② 更关键：**`pc1` 的取值属于柱面**，而失败面（先被判定、`uv=(2.715211,1.000000)`）是**平面** ⇒ **两片 pcurve 与 DS 的面索引对不上**：`int_face_face.rs:248-249` 明确写"结果曲线一律以**排序后**的 `face1_idx=0/face2_idx=1` 携带 pcurves"，而 `pave_ff_perform.rs:464-466` 把它们按**未排序**的 DS 索引 `(n_f1,n_f2)` 存进 `BopdsInterfFf` ⇒ 平面拿到柱面的 pcurve ⇒ `IsValidBlockForFaces` 重建失败 ⇒ 节线块被丢 ⇒ 无节线 edge ⇒ 面不分割。**下一轮入口**：按 OCCT `IntTools_FaceFace::Face1()/Face2()`（排序后的面）与 `BOPAlgo_PaveFiller` 的用法，在 `pave_ff_perform` 里按 `ff.face1()/ff.face2()` 与 `(n_f1,n_f2)` 的对应关系（必要时交换 `pcurve1/pcurve2`）修正；顺带修 `make_pcurve_full` 产出 pcurve 的**无穷声明区间**。

两条修法（已写入 T-80 行，下一步择一）**：① 解析臂按交点解析式**直接构造 2D 圆/直线**（与 3D 同参）——最忠实；② 在 `pcurve_of_curve` 里把 `make_pcurve_full` 的结果**重参数化到 `[range.first, range.last]`**（端口目前只有 3D 的 `occt_geom::curve_reparam::reparameterize_curve`，需补 2D 版）。
- **判据（修好后必须复跑）**：`box∪cyl` 探针应从 `reconstructed=false / total section edges 0 / b.result=两个未分割实体` 变为 `reconstructed=true / 建出节线 edge / b.result 出现被切面且体积 ≈8.5`，随后复跑四道 STEP 门禁 + `export_data_obj` 与 `bop_curved` 直连用例。
- **本批无行为代码改动**：门禁与基线不变（`occt-topo --lib` 1291/1、`phase3/4/6`、四道 STEP 门禁、`export_data_obj` 16/16）。

**批 41（T-80 定位到具体谓词 + `IntTools_Context::ProjPS` 改走忠实投影）—— 2026-09-20 第 64 轮**

- **插桩（临时，已 `git checkout --` 复原；`git grep ZZ_DBG`/`zz_probe` = 0）**：`pave_ff_make::make_blocks_ff` 逐对 dump（`ff/pair/nb_c/nb_p`）＋在每个 `continue;` 前打行号＋`valid_block_point_for_face` dump 每个面（`has_pc/uv/reconstructed/2D-or-3D/result`）；再配一个 driver example 跑 `box∪cyl` 与平面参考。
- **实测（决定性）**：
  - `box∪cyl`：`ff=0 pair=(11,36) nb_c=1` → `face=<key> has_pc=true uv=(2.715211,1.000000) reconstructed=false` → 3D 分支 `result=false` → `continue at L394` → **`total section edges before micro removal: 0`**（节线 edge 一条都没建）。
  - 平面参考：每对 `reconstructed=true` → 2D 分支 `result=true` → 正常建边、正常分割。
  ⇒ FF 找到了交线，但 **`IsValidBlockForFaces` 把它判为无效而丢弃**：`pc.d0(t)` 得到的 UV 无法重建 3D 中点（`uv.y=1.0` 恰是 3D 曲线中点参数，而柱面 UV 的 v 也应为 1.0、u 却对不上）⇒ **根因是 FF 曲线 pcurve 与 3D 曲线参数不同系**（OCCT 的 `IntTools_Curve` 保证同参）。下一轮入口：在 `int_face_face`/`make_pcurve_full` 生成"平面×柱面"pcurve 处对照 OCCT 修参数化。
- **顺带落地一处忠实化（本批唯一行为改动）**：`int_tools_full/context.rs::IntToolsContext::project_point_on_face` 的非平面分支由自创的 `surface_closest_params(surf, p, 32, 32)`（A1 替代件）改为忠实 `occt_geom::geom_api::project_point_on_surface`（`GeomAPI_ProjectPointOnSurf` = `Extrema_ExtPS`，初等面精确；未 done 时仍退回网格兜底并在注释写明出处）。该谓词在布尔管线上被大量调用（`IntTools_Context::ProjPS`，`IntTools_Context.cxx:617-630`）。
- **验证（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（T-05）、`export_data_obj` **16/16**。（注：本批的投影改动**未**修复 T-80 —— 探针复跑后仍 `continue at L394`，说明该对的失败不是投影精度所致，而是 pcurve 参数系。）

**批 40（T-80 再订正：`b.result()` 就是两个未分割的实体 —— 断点在"FF 之后到面分割"这一段，不在 `build_solid`）—— 2026-09-20 第 63 轮**

- **插桩（临时，已 `git checkout --` 复原，`git grep ZZ_DBG` = 0）**：在 `bop_bop::build_rc` 的 Fuse 分支 dump `b.result()` 的类型/子形状数与 `explore(Solid)` 的每个实体的面数与键；`build_solid` 里 dump `mfs`/`sfs`（含面类型）/`dmsts`/`BuilderSolid::areas`。
- **实测（FUSE，`box[-1,1]³` vs `cyl(r=0.4,z∈[0,2])`）**：
  - `build_rc(Fuse): b.result type=Compound children=2 explore(Solid)=2 kept=["Solid/6f", "Solid/3f"]` ⇒ **GF 结果就是参数本身的 6 面盒体 + 3 面柱体，面从未被节线切开**（3 面柱体 = 柱面 + 两片平面端盖，这也解释了批 39 看到的"8 Plane + 1 Cylinder = 9 面"：那是 6+2 个平面面加 1 个柱面，不是"被切过的面"）；
  - 于是 `mfs=3`、`sfs=3`（Plane,Plane,Cylinder）、`BuilderSolid areas=**0**`，`dmsts=1` ⇒ 最终输出只剩那个**未改动的盒体**（6 平面面）；
  - 对照平面参考：`kept=["Solid/6f","Solid/6f","Solid/6f"]`、`mfs=16`、`sfs=14`、`BuilderSolid areas=1 ["Solid/14f"]` ✅。
- **两次订正的意义（方法论留痕）**：批 38 误读 `pave_blocks`（边分块）为"FF 没找到交线"；批 39 又据此把责任推给 `build_solid`；实际链条是 **FF 找到节线（`nb_curves=1`）→ 但没有产出节线 edge / 没有发生面分割 → GF 结果 = 原参数 → `build_solid` 无从装配**。凡涉及管线，必须**逐段 dump 形状组成**（面数/类型/键）而不是只看计数式统计。
- **下一轮入口（已写入 T-80 行）**：在「FF 之后到面分割」这一段插桩——① `pave_ff_post`/`pave_ff_make`/`pave_blocks` 里该对的 section edge 是否被创建、其 pcurve 是否落在盒体顶面与柱面上；② `FillImages`/`BuilderFace`（`bop_build_common`/`builder_face_occt`）之后 `b.result()` 的面数（现在恒等于原参数面数 ⇒ 分割从未发生），对照 OCCT `BOPAlgo_PaveFiller::PostTreatFF`/`MakeBlocks` 与 `BOPAlgo_Builder::FillImagesFaces`。
- **本批不改行为代码**，门禁与基线不变（`occt-topo --lib` 1291/1、四道 STEP 门禁、`phase3/4/6`、`export_data_obj` 16/16）。

**批 39（T-80 精确定位并**订正**批 38 结论：断点在 `build_shape` 的实体装配，不在 FF/GF）—— 2026-09-20 第 62 轮**

- **插桩（临时，均已删；`git checkout --` 复原，`git grep zz_probe` = 0）**：① `pave_ff_perform::perform_ff` 里按 `ZZ_DBG_FF` dump 每一对的 `kinds/tol/starts/failed/nb_curves`；② `bop_builder2::arguments::perform_internal` 里按 `ZZ_DBG_GF` dump GF 结果与 `build_shape` 之后的组成（面数 + 每面 `classify_surface`）。
- **实测数据（`box[-1,1]³` vs `cyl(r=0.4,z∈[0,2])`，FUSE）**：
  - FF 阶段：**只有一对需要求交** —— `pair (11,36) kinds=(Plane,Cylinder) failed=false nb_curves=1` ⇒ 交线**找到了**；
  - GF 装配：`GF result: faces=9 kinds=[Plane×8, Cylinder]`（顶面被 r=0.4 的圆切开、柱面被 z=1 切开）⇒ 分割**正常**；
  - `build_shape` 之后：`type=Compound faces=6 kinds=[Plane×6]` ⇒ **柱面 + 另两块平面被实体装配丢掉**。
  - 对照（平面参考 `box ∪ box(重叠)`）：FF 16 对全部 `failed=false`、GF `faces=16`、`build_shape` 后 `faces=14` ✅ 正常。
- **订正批 38 的错误结论（留痕）**：批 38 据"`box vs cyl` 只有 1 个带 pave block 的形状 / 2 条 pave block（对照 16/24）"判定"面–面阶段没产出交线"——**这是误读**：`pave_blocks` 是**边**的分块计数，一条节线只落在新生成的 section edge 上，与"FF 是否找到交线"无关。正确结论见上：FF 与 GF 都正常，**断点在 `bop_bop::build_solid`（`BOPAlgo_BOP::BuildSolid`）/`builder_solid`（`BOPAlgo_BuilderSolid`）**。
- **下一轮入口（已写入 T-80 行）**：在 `build_solid` 里 dump `mfs`/`mu_sols`/`dmsts` 与壳分组，看 9 面 RC 为何只装配出 6 面盒体（怀疑点：只返回了第一个/最大的壳、或柱面片未进入任何壳而被静默丢弃），对照 `BOPAlgo_BOP.cxx` 的 `BuildSolid` 与 `BOPAlgo_BuilderSolid.cxx`。
- **本批不改行为代码**，门禁与基线不变（`occt-topo --lib` 1291/1、四道 STEP 门禁、`phase3/4/6`、`export_data_obj` 16/16）。

**批 38（T-80 定位：断点在 pave filler 的面–面阶段，不在曲面求交本身）—— 2026-09-20 第 61 轮**

- **探针 1（曲面求交，临时，已删）**：`intpatch::intersection_curve_points(plane z=1, cylinder r=0.4)` → **4 条折线**，逐点半径 0.399999~0.400001、z 全为 1.000000 ⇒ 采样型求交（`intpatch_trace`）能给出这个圆（注意：同一圆被返回 4 次，属另一处待查的重复，未追）。
- **探针 2（pave filler，临时，已删）**：构造 `PaveFiller::new()` + `set_arguments([box, cyl])` + `perform()`，再统计 `ds()` 里带 pave block 的形状数：
  - `box vs cyl` → `shapes=41, with_pave_blocks=1, total_pave_blocks=2`（唯一带块的是一个 Edge，2 条）；
  - 对照 `box vs box(重叠)` → `shapes=68, with_pave_blocks=16, total_pave_blocks=24`（16 条边各带块）。
  ⇒ **面–面（平面∩柱面）阶段没有产出可用交线**：GF 无法分割面、`build_result`/`build_solid` 随后把工具整体丢掉，于是 `FUSE` 只剩 box（6 面、vol 8.0000）、`CUT` 的 box 不变、`COMMON` 为空。三条 BOP 通道（`boolean` / 裸 `builder_bop_with_fuzzy` / `curved_boolean_full`）表现一致，印证问题在共享的 filler+GF 上，而不在各自的分派或 `BuildRC`。
- **结论与下一轮入口**：T-80 的修复点在 **pave filler 的 FF 阶段**（`pave_ff_perform::perform_ff` → `pave_intersect`），需 dump "平面面 × 柱面面" 这一对的相交器选择与返回曲线数（候选：`IntAna_QuadQuadGeo`（平面∩二次曲面解析）与 `IntPatch_ImpImpIntersection`），对照 OCCT `BOPAlgo_PaveFiller::PerformFF`（`BOPAlgo_PaveFiller_*.cxx`）与 `IntTools_FaceFace`。已写入 T-80 行。
- **本批不改行为代码**（两次探针均已删，`git grep zz_probe` = 0；`git status` 干净），门禁与基线不变（`occt-topo --lib` 1291/1，四道 STEP 门禁、phase3/4/6、`export_data_obj` 16/16）。

**批 37（T-79 判定性实验：曲面分支整体迁移到 `BOPAlgo_BOP` → 探针证否 + 回退；登记 T-80）—— 2026-09-20 第 60 轮**

- **做了什么（两步，均已回退）**：① `boolean_dispatch` 的最后一行由 `bop_curved::curved_boolean_full` 改为 `bop_builder::boolean`（= 已移植的 `BOPAlgo_BOP`）；② `curved_boolean_full` 自身也改为直接调 `bop_builder::boolean`（让它对外的 4 处调用点一并忠实化）。
- **门禁结果（危险地"全绿"）**：两步之下 `occt-topo --lib` **1291/1**（唯一红仍 T-01）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、四道 STEP 门禁 14/14、13/13、11/11、2/3、`export_data_obj` **16/16 计数逐位一致**。若就此收工，会得出"迁移完成"的错误结论。
- **判定性探针（临时 `examples/zz_probe_bop79.rs`，已删，`git grep zz_probe` = 0）**：对 `box[-1,1]³` 与 `cyl(r=0.4, z∈[0,2])` 分别走三条路径，统计结果的面类型与体积：
  - `FUSE box ∪ cyl`：`bop_builder::boolean` / 裸 `builder_bop_with_fuzzy` / `curved_boolean_full` **三者都给 6 面全平面、vol 8.0000**（= 只剩 box，柱面段与柱面全丢；正确应 ≈8.5 且保留柱面）；
  - `CUT box − cyl` → box 不变（6 面、vol 8.0000，应有内孔）；
  - `COMMON box ∩ cyl` → **空**（0 面，应非空）；
  - 同引擎**平面**用例 `FUSE box ∪ box(重叠)` → **14 面、vol 10.0000**（正确）⇒ **缺口在曲面分支**（`BuildRC`/面状态分类/分割），与 T-01 同族。
- **为什么门禁没抓到**：会抓到这个缺口的用例是 `bop_curved::tests::{trimmed_face_box_cylinder_fuse, curved_boolean_full_quadric_unchanged}`，它们**直接调 `curved_boolean_full`**；第①步只改 dispatch，故它们仍走旧体而全绿；第②步改到它们头上后两例立刻红（`has_cylinder` 断言失败 / "quadrics route to curved_boolean" 的分派断言失败）。⇒ 已在 §9 记一条坑：**"门禁绿 ≠ 可迁移"，必须先看会被改动的路径有没有测试直接绑在旧体上。**
- **处置**：按"失配即停"**回退两步**（`git checkout -- bop_builder_dispatch.rs bop_curved/region_trim.rs`），复跑 `occt-topo --lib` **1291/1** 确认恢复；把探针证据登记为 **T-80**（曲面 operand 下 `BOPAlgo_BOP` 移植体失效），并标注 **T-79 与 A1 的最终迁移被 T-80 阻塞**。第 59 轮的 `voxel_fallback` 摘除**保留**（其改动只涉及空/平面输入，探针已证 BOP 在这些用例上正确）。
- **本批不改行为代码**（净效果 = 新增 T-80 登记 + T-79/A5 依赖说明 + §9 坑记录），门禁与基线一致。

**批 36（A5/T-41：删掉 dispatch 级体素布尔 `voxel_fallback`，空/非平面输入改走忠实 `BOPAlgo_BOP`）—— 2026-09-20 第 59 轮**

- **对读的 OCCT 事实**：OCCT 没有体素/网格布尔——`BRepAlgoAPI_*` 一律走 `BOPAlgo_BOP`（`BOPAlgo_Builder` + `IntPatch`）。空实参由 `BOPAlgo_BOP::CheckData` 记 `BOPAlgo_AlertEmptyShape` 后**跳过**（`BOPAlgo_BOP.cxx:162-167`，空组取另一组维数 `:203-209`），平面/曲面输入走**同一套**引擎。
- **改动**：① `bop_builder_core::voxel_fallback`（`32³` `boolean_ops::voxel_boolean` + `mesh_to_brep` 重建）**删除**；② `bop_builder::boolean` 的空面输入（原 `:48-50`）与 `bop_builder_planar` 的空/非平面输入（原 `:44-52`）改为调 `crate::bop_builder2::builder_bop_with_fuzzy`（= 已移植的 `BOPAlgo_BOP::Perform`+`BuildShape`），并把 `bop_builder_dispatch::to_bool_op2` 提升为 `pub(crate)` 复用；两处就地写明 OCCT 出处；③ `bop_builder_planar_geom/trace` 的死导入清理；④ `feature.rs::boolean_feature` 的 `voxel_boolean` 就地标 **UNPORTED**（OCCT 的特征走 `BRepFeat_MakePrism/MakeRevol` + `BOPAlgo_BOP`，属 A12/T-48）。
- **验证（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01；`bop_builder_tests*` 中所有以 `boolean()` 为入口的用例仍全绿）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（T-05）、`export_data_obj` **16/16** 计数逐位一致 ⇒ dispatch 级摘除安全。
- **余项（如实登记）**：`bop_curved` 的**网格布尔**（`boolean_mesh`/`boolean_mesh_curves`/`region_inside_other`/`general_boolean_trimmed`）仍在 `boolean_dispatch` 的曲面分支上、**门禁当前正走此路**，迁移到忠实 BOP 需逐项回归分诊 ⇒ 立项 **T-79**。

**批 35（登记，无代码改动）：IGES 写侧替身立项 T-78 —— 2026-09-20 第 58 轮**

- **核查结论**：`crates/occt-topo/src/iges.rs` 的**写侧**是部分替身——发射器只有 `emit_line`(110)、`emit_circular_arc`(100)、`emit_face_surface` 里为球面**自造**的 120 回转面（`iges.rs:409-434`，OCCT 是 **196**）、`emit_solid`(186 MSBO)；曲线族靠 6 点二阶导采样 `(max_d2-min_d2)/max_d2 < 0.02` 猜圆（`:158-172`），`CurveKind::Other` 直接退化成弦线（`:389-395`）。OCCT 对应件是 `GeomToIGES_GeomCurve`（按 `GetType()` 发 100/104/106/108/110/126）与 `GeomToIGES_GeomSurface`（192/194/196/198/120）+ `IGESData_*` 参数布局 + 144 裁剪面。
- **处置**：按纪律登记为 **T-78**（后续补发射器时逐类对照上述文件），本轮**不改代码**；A26 的其余项是 PLY 每面 node 块（需 `RWPly_PlyWriterContext` 的逐面重复 node+偏移语义）。⇒ A26 仍未完全结项，但每一项都有了可跟踪出处。

**批 34（未登记自创：`int_tools_curve_box` 的 48 点采样曲线包围盒改为忠实 `GeomBndLib_Curve`；删掉无调用者的采样距离函数）—— 2026-09-20 第 58 轮**

- **对读的 OCCT 事实**：`BndBuildBox`（`IntTools_EdgeEdge.cxx:1410-1419`）就是 `BndLib_Add3dCurve::Add(theBAC, aT1, aT2, theTol, theBox)`；8.0.0 里该重载经 `GeomBndLib_Curve` 按 `GetType()` 分派（`GeomBndLib_Curve.cxx:101-172`：Line/Circle/Ellipse/Hyperbola/Parabola/Bezier/BSpline/Offset，其余 `OtherCurve`），容差即盒的 gap —— **不是采样**。
- **改动（`crates/occt-topo/src/int_tools_curve_box.rs`）**：① `add_curve_to_box` 由"48 点均匀采样 + `enlarge(tol)`"改为调用仓内已有的忠实件 `geom_bnd_lib_curve3d::box_curve(curve, a, b, tol)`（= `GeomBndLib_Curve::Box`，含解析臂与 gap 语义；该忠实件此前只有 `brep_bnd_lib`/`geom_bnd_lib_surface3d`/`pcurve_full` 在用），保留下方 `finite_or_unit` 的无界曲线保护（端口曲线可能是 `±inf`，OCCT 的 adaptor 是 `±Precision::Infinite()`）；② 删除 `CURVE_BOX_SAMPLES` 常量与**无任何调用者**的 `sampled_distance_to_curve`（48 点采样距离，属同类未登记自创）。⇒ 该文件的"采样替代件"清零，同时为 **T-77** 备好 `BndBuildBox` 前置。
- **验证（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（T-05）。`add_curve_to_box` 的生产调用方（`bean_face_localize.rs:68,442`、`pave_ee_aux.rs:29`、本文件的 `curve_box_unexpanded`/`CheckCurve`）全部走新路径且无需改签名。

**批 33（A11/T-47 子项 3：删掉 `edge_edge` 的"采样解当补集"，A11 结项；忠实 bbox 递归立项 T-77）—— 2026-09-20 第 57 轮**

- **对读的 OCCT 事实**：`IntTools_EdgeEdge::FindSolutions` **不枚举 `Extrema_ExtCC` 的解**——它按参数盒递归细分（`IntTools_EdgeEdge.cxx:290-549`），用 `BndBuildBox`（`:1410-1419`）、`FindParameters`（`:553-671`，沿曲线按 `Resolution` 自适应步长找"落在对方盒内的参数"）、`IsIntersection`（`:1060-1146`，含角度准则与 `FindDistPC` 复核）、`CheckCoincidence`（`:1150-1206`）与 `SplitRangeOnSegments`（`:1366-1406`）逐步收敛，最后 `MergeSolutions`（`:675-779`）把区间归并成 `CommonPrt`/点。**OCCT 里没有任何"把另一个采样器的解并进来"的分支**。
- **改动（`crates/occt-topo/src/edge_edge/edge_edge.rs`）**：删掉 `find_solutions` 里"用 `inttools::edge_edge_intersections` 采样结果当补集并入"的那 7 行，只保留替代判据"距离在容差内的极值即交点"（其来源已是忠实 `Extrema_ExtCC`/`GGenExtCC`，见 T-66）；函数文档改写为**就地 UNPORTED**，写明 OCCT 的递归及其所需各函数与行号。`self.edge1/edge2` 的 clone 随之不再需要，一并去掉。
- **验证（删除后逐项仍等于基线 ⇒ 该补集是冗余的自创）**：`occt-topo --lib` **1291/1**（唯一红 = T-01；`edge_edge` 全部单测绿，含 `line_bspline_crossing_two_hits` 两个交点）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（T-05）、`export_data_obj` **16/16** 计数逐位一致。
- **结项与遗留**：A11 的三项自创（`boolean_degenerate` 质心规则、丢"顶点<3"面、采样解当补集）全部处理完毕 ⇒ **A11 ✅**。OCCT 真正的 bbox 递归未移植，按纪律登记为 **T-77**（属"替代件升级为忠实件"，不是未登记自创），其完全忠实化还需 `BndLib_Add3dCurve` 的忠实件（端口目前只有 `int_tools_curve_box` 的 48 点采样近似）。

**批 32（A7/T-66 步 b-2：`Extrema_GGenExtCC::Perform` 逐段移植并接线，自创种子集与采样分类器全部删除；A7 结项）—— 2026-09-20 第 56 轮**

- **对读的 OCCT 事实**：`Extrema_ECC` **就是** `Extrema_GGenExtCC<…>`（`Extrema_ECC.hxx:23-34`），`Extrema_ExtCC` 对"直线+初等曲线"或"两圆"用解析 `Extrema_ExtElC`，其余一律走它（`Extrema_ExtCC.cxx:247-316`，并按 `GetType()` 分派 `cxx:251-305`）。`Perform`（`Extrema_GGenExtCC.hxx:453-802`）分五段：① 区间集（`GeomAbs_C2` 起步，乘积 >100 退 `C1`；长度比 `mult=20` 细分；闭合单区间曲线 3 分；`hxx:461-541`）；② Lipschitz 估计（`aMaxDer = max(1/C1.Resolution(1), 1/C2.Resolution(1))*√2`、21×21 梯度扫描、`aLC` 与 `isConstLockedFlag`、直线特例；`hxx:543-616`）；③ 每个区间对 `SetLocalParams` + `Perform(GetSingleSolutionFlag())`（`hxx:639-649`），用 `aSameTol*aValueTol` 维护全局最优、`PointsInspector` 按 `aCellSize` 去重（`hxx:650-688`）；④ `aNbSol==0` 未 done、`==1` 直接给解、否则按 `comp`（先 X 后 Y）排序（`hxx:690-707`）；⑤ 平行/同向分析与 `ProjPOnC` 复核，最后组装 `myPoints1/2`（`hxx:707-802`）。
- **改动**：① 新增 `occt-geom/src/extrema_cc/general_extrema.rs`（`GGenExtCC`）：上述五段逐段移植（含 `change_intervals`=`Extrema_GGenExtCC_ChangeIntervals`、`proj_p_on_c`=`Extrema_GGenExtCC_ProjPOnC`（用忠实 `extrema_ext_pc_range` 取所有解的最小平方距离）、`is_closed`（`Geom_Circle.cxx:26-29` 恒真 + `Geom_TrimmedCurve.cxx` 整周期支））；② `occt-math/globoptmin.rs` 拆出 `set_global_params`（=ctor 的 `SetGlobalParams`，含 `ComputeInitSol`）与 `perform_local`（=OCCT `Perform`，**不重置**最优值与解表 ⇒ 跨区间累积），`perform` 保留为二者的便捷组合（原测试不动）；为此 `occt-geom` 增加 `occt-math` 依赖（无环：`occt-math → occt-core`）；③ **接线**：`curve_curve_extrema_all` 拆出 `curve_curve_extrema_all_range(c1,c2,u1,u2,v1,v2)`（= `Extrema_ExtCC` 的 `Initialize`/`SetParams` 范围形式），解析臂改由 **类型查询** `gp_line`/`gp_circ` 分派；`edge_edge::find_solutions` 改用**边范围**调用（OCCT 那边是 `BRepAdaptor_Curve`，`IntTools_EdgeEdge.cxx:94-95`），因为端口曲线的自身范围可能是无穷；④ **删除全部自创**：`newton_curve_curve_all`（网格+变号+局部极值+边界种子的种子集）、`build_samples`（含扩窗采样）、`is_line`/`classify_circle`/`line_of_curve`/`circumcenter`（采样分类器，A7/A16 同族）。
- **过程中修掉的两处真实缺陷（留痕）**：① **无穷边界**：端口无界曲线报 IEEE ±inf，而 OCCT 的 `Adaptor3d_Curve` 报 `±Precision::Infinite()`（有限大数，`Geom_Line::FirstParameter()`），于是 `math_GlobOptMin` 的 `(A+B)/2` 会得 NaN 而整个引擎报错；现按 OCCT 数值约定把边界与区间中的无穷替换为 `precision::INFINITE`（已就地注明）。② **区间必须取"受限 adaptor"**：OCCT 的 `Intervals` 来自范围受限的 adaptor，端口直接用曲线自身区间（直线是 ±inf）会让局部盒退化成 2e100 宽、根本分辨不出小参数处的极值；新增 `restricted_intervals`（自身区间与边界盒求交）后，`edge_edge::tests::line_bspline_crossing_two_hits` 由 **0 命中 → 2 命中**，且与解析值逐位吻合（`u1 = 2∓√2 = 0.585786/3.414214`，`u2 = (1∓1/√2)/2 = 0.146447/0.853553`，`d = 0.000000`）。
- **探针实测（临时，均已删，`git grep zz_probe` = 0）**：直线/二次 Bezier 且直线为无界曲线时报错 → 修无穷边界后得 `u=2.0, v=0.5, d=2.0`（解析最小距离 2，Bezier 参数 0.5）；两圆（r=1，(0,0) 与 (3,0)）经解析臂给 min/max；边范围内的直线/Bezier 交叉给两个 d=0 的解（见上）。
- **验证（逐项等于基线）**：`occt-math --lib` **215/215**、`occt-geom --lib` **151/151**、`occt-topo --lib` **1291/1**（唯一红 = T-01）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（T-05）、`export_data_obj` **16/16** 且计数逐位一致。

**批 31（A7/T-66 步 b-1：`Extrema_GlobOptFuncCCC2` 忠实移植——曲线对距离函数的值/梯度/Hessian）—— 2026-09-20 第 55 轮**

- **对读的 OCCT 事实**：`Extrema_GlobOptFuncCC.hxx:25-99` 声明三个类 `CCC0`/`CCC1`/`CCC2`（"C0/C1/C2 continuity"），三者**数学相同**，只是分别实现 `math_MultipleVarFunction` / `WithGradient` / `WithHessian` 三个接口——存在的唯一理由是 `math_GlobOptMin::computeLocalExtremum` 用 `dynamic_cast` 选局部引擎（`math_GlobOptMin.cxx:266-339`）。2 变量（C1 参数 u、C2 参数 v）的静态式在 `Extrema_GlobOptFuncCC.cxx:24-185`：`_Value` 越界返回 false（`cxx:38-42`）、否则 `F = |C2(v)-C1(u)|`（`cxx:44`）；`_Gradient` 为**平方距离**的梯度 `(2[(C1-C2)·C1'], 2[(C2-C1)·C2'])`（`cxx:87-91`）；`_Hessian` 为平方距离的 Hessian（`cxx:140-151`）。⇒ **OCCT 自身的不一致：值返回距离、导数返回平方距离的导数**（极小点相同，只影响步长），端口按原样保留并注明。
- **改动**：新增 `crates/occt-geom/src/extrema_cc/glob_opt_func.rs`（`GlobOptFuncCCC2`）：`new`/`nb_variables`/`value`（距离，越界 `None` = OCCT 的 `false`）/`gradient`/`values`/`hessian`，公式逐行对应上述静态式（含 `GpVec::from_pnts` 的 Δ、`d1/d2` 的一二阶导），并在模块头写明"三类的差别只是接口、非数学"与"端口 `GlobOptMin` 只吃值闭包 ⇒ 引擎选择仍为 UNPORTED（见 `globoptmin.rs::compute_local_extremum`）"。`mod.rs` 注册 `glob_opt_func`。
- **探针实测（临时 `examples/zz_probe_globfunc.rs`，已删，`git grep zz_probe` = 0）**：C1 = 原点沿 X 的直线、C2 = `(0,1,0)` 沿 Y 的直线（解析：`d² = u² + (1+v)²`）⇒ `value` 与解析 `sqrt(u²+(1+v)²)` **逐位一致**（1.000000000000 / 1.581138830084 / 4.472135955000）；`gradient` = `(2u, 2(1+v))` 与解析完全一致，并与 `d²` 的中心差分吻合到 `4.0e-11`/`1.5e-11`；`hessian` = `[[2,0],[0,2]]` 与解析一致；越界（圆 `u=7 > 2π`）返回 `None`（= OCCT 的 `false`）、`nb_variables()=2`。
- **验证（逐项等于基线）**：`occt-geom --lib` 151/151、`occt-topo --lib` **1291/1**（唯一红 = T-01）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（T-05）、`export_data_obj` **16/16**。（本件是 T-66 的前置步，尚无生产调用方 ⇒ 门禁为等价性兜底。）

**批 30（A7/T-66 步 a：`math_GlobOptMin` 补齐四个 setter，并对齐 OCCT 的"全局盒 → 局部盒 → Perform"次序）—— 2026-09-20 第 54 轮**

- **对读的 OCCT 事实**：`math_GlobOptMin` 的生命周期是 `ctor`→`SetGlobalParams`（`cxx:112-148`：写 `myGlobA/myGlobB/myA/myB`、`myMaxV=(b-a)/3`、`myTol/mySameTol`、`initCellSize()`、**`ComputeInitSol()`**、`myDone=false`）⇒ `SetLocalParams`（`cxx:154-171`：写 `myA/myB`、`myMaxV=(b-a)/3`、`myZ=-1`、`myDone=false`，**不重跑 `ComputeInitSol`**）⇒ `SetTol`/`SetFunctionalMinimalValue`/`SetLipConstState`/`SetContinuity` ⇒ `Perform`（`cxx:192-262`：由**当前** `myA/myB` 算 `minLength/maxLength`、`myV=0`、退化返回、未锁则 `computeInitialValues()` 估 Lipschitz、`myE1/E2/E3`、`CheckFunctionalStopCriteria`（`cxx:600-604`：`myIsFindSingleSolution && |myF - myFunctionalMinimalValue| < mySameTol*0.01`）、`computeGlobalExtremum(myN)`）。`Extrema_GGenExtCC::Perform` 正好按这个次序用：**整个边框建一次 finder**（`hxx:617`），然后每个区间对 `SetLocalParams`+`Perform`（`hxx:639-649`）。
- **改动（`crates/occt-math/src/globoptmin.rs`）**：① 新增 `set_lip_const_state`/`lip_const_state`（`hxx:116-118`）、`set_continuity`/`continuity`（`hxx:102-104`）、`set_functional_minimal_value`/`functional_minimal_value`（`hxx:107-112`）、`set_local_params`/`local_params`（`cxx:154-171`）与 `cont` 字段（默认 2，`cxx:69`）；② `perform` 的次序改为：全局盒建立 → `done=false`/`z=-1` → **`compute_init_sol`（全局盒，对应 ctor 里的 `ComputeInitSol`）** → 应用待定局部盒（`a/b/max_v`，`z=-1`，`done=false`）→ `minLength/maxLength`+`v=0`+退化检查 → Lipschitz → `e1/e2/e3` → 停止判据 → `compute_global_extremum`；此前端口把 `compute_init_sol` 排在 `e1/e2/e3` 之后（无局部盒时二者等价，有局部盒时才分道）⇒ 现在与 OCCT 逐点对应。
- **如实标注的 UNPORTED**：`set_continuity` 在 OCCT 里用于选局部精化引擎（`computeLocalExtremum`，`cxx:266-339`：`myCont>=2`+Hessian → `math_NewtonMinimum`；`myCont>=1`+Gradient → `math_BFGS`；否则 `math_Powell`），而端口的 `GlobOptMin` 只接收**纯值闭包**、无函数类层次 ⇒ 只保留 BFGS 臂且梯度为数值差分，`cont` 仅存储不切换引擎（已在该函数文档与 setter 文档写明）。
- **探针实测（临时 `examples/zz_probe_globopt.rs`，已删，`git grep zz_probe` = 0）**：双井函数 `f = -1/(1e-4+d²(P,(0.25,0.25))) - 1/(1e-4+d²(P,(0.75,0.75)))` 在单位方盒上：**无局部盒** → 解 `(0.2546, 0.1900)`、`f=-270.52`；`set_local_params([0,0],[0.5,0.5])` → 解 `(0.2799, 0.1950)`（**落在给定子盒内**）、`f=-250.77`；`set_local_params([0.5,0.5],[1,1])` → 解 `(0.7814, 0.6950)`（子盒内）、`f=-245.18` ⇒ 局部盒确实限制搜索范围；setter 往返（`locked=true cont=1 fmin=-7.5 local=Some(...)`）与"锁定 Lipschitz 常数"路径均按预期（锁定后同一双井函数得 `f=-4041.16`，与 `SetLipConstState` 语义一致）。
- **验证**：`occt-math --lib` **215/215**（含既有 globoptmin 单测，次序调整后仍全绿）、`occt-topo --lib` **1291/1**（唯一红 = T-01）、`phase3` 4/4、`phase4` 9/9、`phase6` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（T-05）、`export_data_obj` **16/16** 计数逐位一致。（`GlobOptMin` 目前**无生产调用方**，本批是 T-66 的前置步，行为等价性由上述门禁兜底。）

**批 29（A4/T-40：`AlgoTools::compute_state` 的 F 分支改走忠实 `BRepClass_FaceClassifier`，删掉 OCCT 没有的"距离 > tol → Out"门）—— 2026-09-20 第 53 轮**

- **对读的 OCCT 事实**：`BRepClass_FaceClassifier::Perform(F, P, Tol)`（`BRepClass_FaceClassifier.cxx:76-125`）= ① `BRepAdaptor_Surface aSurf(theF, false)`（**不受限**于面的 UV 边界，`cxx:90`）；② `BRepTools::UVBounds(theF, U1,U2,V1,V2)`；③ `aExtrema.Initialize(aSurf, U1,U2,V1,V2, theTol, theTol)` + `Perform(P)`（`cxx:91-97`）；④ `!IsDone()` 或 `NbExt()==0` ⇒ **直接返回**（状态留 `UNKNOWN`、`rejected = true`，`cxx:98-107`）；⑤ 取**最小平方距离**的解（`cxx:109-117`）；⑥ 用该解 `(u, v)` 走 2D 分类器（`cxx:121-123`）。**该函数完全不测 3D 投影距离** —— 面外 0.1 的点只要投影落在面的 UV 域内就是 `In`。
- **改动**：① `algo_tools/construct.rs::AlgoTools::compute_state` 的 `ShapeType::Face` 分支按上表逐行重写（UV 窗口来自忠实件 `crate::brep_uv_bounds::uv_box_of_face` = `BRepTools::AddUVBounds`/`UVBounds`，投影走 T-67 步 2 的 `ExtPs` 分派层）；**删除自创的两处**：`surface_closest_params(surf, p, 32, 32)` 网格投影（A1 替代件）与 `surf.d0(u,v).distance(p) > tol → FaceState::Out` 距离门；② 同族 `algo_tools/queries.rs` 的 `get_edge_off`（33 点采样里逐点 32×32 投影）与 `face_normal_at_point` 改用忠实 `geom_api::project_point_on_surface`（T-52），未 done 的采样点跳过；③ `algo_tools/mod.rs` 去掉不再用的 `surface_closest_params` 转出。**E 分支无需改**：`brep_extrema::closest_point_on_edge` 早已基于忠实 `point_curve_extrema_all`（`Extrema_ExtPC`），其 `samples` 形参是遗留物（已注明）。
- **一处既有单测的期望按 OCCT 订正（留痕）**：`algo_tools/tests.rs::compute_state_face_in_out_on` 的第 2 条断言原本期望"面上方 0.1 的点 → `Out`"，那正是被删掉的自创门的行为；按 `BRepClass_FaceClassifier`（投影 + UV 分类、不测距离）改为 `In`，并在测试里写明出处与订正原因。这不是"为对齐新写测试"，而是修正编码了自创语义的既有断言（同类先例：T-06 的期望订正、T-42/T-47 删除自创单测）。
- **验证（逐项等于基线）**：`occt-topo --lib` **1291/1**（唯一红 = T-01；F 分支的两处自创删除后除上述断言外无其他用例受影响）、`phase3_integration` 4/4、`phase4_integration` 9/9、`phase6_integration` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05）、`export_data_obj` **16/16** 且计数逐位一致。`compute_state` 的生产调用方（`algo_tools_face.rs:546,562`、`bopalgo_tools_class.rs:156-175`）全部传**实体** ⇒ 走未改动的 Solid 分支。
- **未做**：`compute_state` 的容器分支仍是端口自创的"递归取第一个已知状态"（OCCT 没有这个签名，其 `ComputeState` 是"子形状 vs 参考实体"；本函数的忠实对应件是 `BRepClass_FaceClassifier`（面）与 `BRepClass3d_SolidClassifier`（体））——已在函数文档写明；`algo_tools/queries.rs::get_edge_off` 的 33 点采样结构本身仍是自创（OCCT `BOPTools_AlgoTools::GetEdgeOff` 走 pcurve），下一步可另立 T-xx。

**批 28（A16 投影 half：`geom_api` 的点–曲线/点–面投影改走忠实 Extrema；删掉 124 行自创网格+黄金分割）—— 2026-09-20 第 52 轮**

- **对读的 OCCT 事实**：`GeomAPI_ProjectPointOnCurve::Perform`（`GeomAPI_ProjectPointOnCurve.cxx:135-155`）= `Extrema_ExtPC` 在曲线**自身区间**上（`myExtPC.Initialize(myC, FirstParameter, LastParameter)`，`cxx:58`），`IsDone = IsDone() && NbExt() > 0`（`cxx:139`），答案取**最小**平方距离（`cxx:143-152`）；**没有**端点修正（那是 `ShapeAnalysis_Curve::Project` 的事，`ShapeAnalysis_Curve.cxx:161-182`）。`GeomAPI_ProjectPointOnSurf::Perform`（`GeomAPI_ProjectPointOnSurf.cxx:214-246`）= `Extrema_ExtPS`（`Initialize(adaptor, Umin, Usup, Vmin, Vsup, Tol, Tol)`，`cxx:128`），`IsDone` 规则同上（`cxx:83`），答案取最小 `SquareDistance`（`cxx:88-100`）。
- **改动（`crates/occt-geom/src/geom_api.rs`）**：① `project_point_on_curve` 改调 `crate::extrema_pc::point_curve_extrema`（= 已忠实化的 `Extrema_ExtPC`，T-43；`None` 即 OCCT 的 `IsDone()==false`）；② `project_point_on_surface` 改调 `crate::extrema_surf::ExtPs::with_surface`（T-67 步 2 的分派层）+ 取最小，未 done/无解返回 `None`；③ 删掉自创件 `refine_closest_curve`（64 点扫描 + 黄金分割）、`closest_params_in_window`（16×16 栅格 + 6 轮爬山 + 逐轴黄金分割）、`surface_closest_params`（-1/1 起、×8 扩窗直到"解在内部"）共约 124 行，并让 `dist_curve_surface`/`pcurve_of_curve_on_surface` 走同一忠实投影；④ 模块头改写为"投影已忠实 / 求交 UNPORTED"，两个求交采样器就地标 `UNPORTED` 并写明真出处（`GeomAPI_IntCS`→`IntCurveSurface_Intersection`+`IntPatch_Intersection`；曲线–曲线走 `IntTools_EdgeEdge`，**8.0.0 无 `GeomAPI_IntCC`**）与"生产仍在使用、只能替换不能删"。
- **影响面（实测）**：`project_point_on_curve` 有 10 处生产调用（全在 `edge_edge/edge_edge.rs`+`edge_edge/solvers.rs`，对应 OCCT `IntTools_EdgeEdge::FindSolutions` 里的 `GeomAPI_ProjectPointOnCurve`，`IntTools_EdgeEdge.cxx:486-506`）；`project_point_on_surface` 有 4 处（`geom_int_intss_pcurve.rs`、`geom_int_quadric.rs`、`int_face_face_bounds.rs`）。
- **验证（逐项等于基线）**：`occt-geom --lib` 151/151、`occt-topo --lib` **1291/1**（唯一红 = T-01；`edge_edge` 自身单测全绿）、`phase3_integration` 4/4、`phase4_integration` 9/9、`phase6_integration` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05）、`export_data_obj` **16/16** 且 16 个模型顶点/面数与基线**逐位一致** ⇒ 在现有语料上"自创近似 → 忠实精确"的结果一致，改动仅去掉了自创算法本身（初等面的投影现在是 `ExtPElS` 解析解，一般曲面仍是 `GenExtPS` 替代件，已在 `extrema_surf/point_surface_extrema.rs` 标 UNPORTED）。
- **未做（T-52 求交 half）**：`curve_surface_intersections`（曲线采样 + 距离凹陷检测 + 黄金分割）与 `curve_curve_intersections`（256×256 双网格 + 交替一维极小 + 去重），两者被 `inttools/intersections.rs`（`:278/:285/:292/:295/:450`）生产调用；忠实化需要 `IntCurveSurface_Intersection`/`IntPatch_Intersection` 或 `IntTools_EdgeEdge::FindSolutions`（A11 子项 3 同一前置）。同批未动 `occt-geom2d/curve_ops.rs:68-69`（2D 侧，属 A16 的另一半）。

**批 27（A15/T-51 余项：`curve_reparam` 的基函数改走忠实 `BSplCLib::EvalBsplineBasis`；订正"缺 `BSplCLib::BasisFuns`"这一假前提）—— 2026-09-20 第 51 轮**

- **订正审查结论（假前提）**：A15 余项原写"`curve_reparam.rs:189` 的 `1e-15`（缺忠实 `BSplCLib::BasisFuns`）"。**全树 grep 实测：OCCT 8.0.0 里不存在 `BasisFuns`**（`BSplCLib.cxx`/`.hxx`/`BSplCLib_2.cxx` 均无此名，全仓 `*.cxx/*.hxx/*.lxx` 0 命中）⇒ 老名字已被 8.0.0 移除。8.0.0 的真出处是 **`BSplCLib::EvalBsplineBasis`（`BSplCLib_2.cxx:429-563`）**，由 `BSplCLib::Eval`（`BSplCLib.cxx:3526`）与 `BSplSLib` 调用；其退化判据是 `std::abs(aScale) < gp::Resolution()`（`cxx:509-512`/`:538-541`），**不是**端口自创的 `1e-15`；错误码 `2` 时 OCCT 调用方**直接放弃**（`BSplCLib.cxx:3532-3535`）。
- **落地 1（新增忠实件 `occt-core/src/bspl/eval_basis.rs`）**：逐行移植 `EvalBsplineBasis`——`LocateParameter`（`BSplCLib.cxx:189-214` 的 `Degree/FromK1/ToK2` 重载，非周期走 `(0.,1.)`、周期走 `[k_degree, k_{n}]`）→ `FirstNonZeroBsplineIndex` → 两段增量循环（含 `math_Matrix` 的**行主序**平坦索引，`NCollection_Array2.hxx:317` 已核对）；矩阵尺寸按 `BSplCLib_LocalMatrix(LocalRequest, Order)` 分配，故 OCCT 的"矩阵过小"错误码 1 不可达（已注明）；返回 `(FirstNonZeroBsplineIndex, basis)` 或 OCCT 错误码 `2`。
- **落地 2（`occt-geom/src/curve_reparam.rs`）**：删掉 The NURBS Book A2.2 的 `basis_values`（含 `1e-15`），改调忠实件并把窗口起点由 `FirstNonZeroBsplineIndex` 散射到逐极点向量；`resample_bspline` 在 `EvalBsplineBasis` 报错时返回显式错误（对应 OCCT 调用方的"放弃"），不再伪造基值。
- **探针实测（临时 `examples/zz_probe_basis.rs`，已删，`git grep zz_probe` = 0）**：三次夹持曲线（6 极点、结点 `0⁴ 1 2 3⁴`）⇒ 六个采样点 `ΣB_i(u) = 1.000000000000000`（分区性）、`FirstNonZeroBsplineIndex` 随区间正确推进 `1→2→3`、散射后与端口 de Boor 求值 **≤5.6e-16**；一阶导行与中心差分 **≤8.3e-10**（差分尺度）且 `ΣB'_i(u) = 0`；全同结点 → `Err(2)`（= OCCT 的退化跨度码，判据为 `gp::Resolution()`）；周期结点（order 3，周期窗 `[k₂,k₄]`）在 `u=0.5/2.5/4.5/5.5` 归一后给同一组值 `[0.125, 0.75, 0.125]`（二次 B 样条跨度中点解析值）。
- **验证（逐项等于基线）**：`occt-core --lib` 290/290、`occt-geom --lib` 151/151、`occt-topo --lib` **1291/1**（唯一红 = T-01）、`phase3_integration` 4/4、`phase4_integration` 9/9、`phase6_integration` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05）、`export_data_obj` **16/16** 且计数与基线逐位一致。
- **旁支观察（未改）**：探针一开始喂了**不自洽**的输入（6 极点 + 长度 8 的平坦结点），端口的两处 port-internal 求值器 `bspl/eval.rs::eval_curve`（`:26`）与 `eval_curve_d1`（`:84`）随即**越界 panic**（`len is 8, index is 8`），而 OCCT 的 `BSplCLib::Eval` 在不自洽输入下同样只是读越界（无契约检查）⇒ 属**输入契约**问题（端口所有生产调用点都传自洽结点，实测 u=0/0.25/1/1.5/2.5/3 均正常），仅登记不再追。

**批 26（A11/T-47 子项 1：删掉自创的 `boolean_degenerate` 整套，非实体输入改由 `BOPAlgo_BOP` 判定）—— 2026-09-20 第 50 轮**

- **OCCT 事实（逐行核对 `BOPAlgo_BOP.cxx:140-210`）**：OCCT **没有**"非实体/退化输入的尽力而为层"。`CheckData` 是唯一判据：① 空实参（`BOPTools_AlgoTools3D::IsEmptyShape`）记**警告** `BOPAlgo_AlertEmptyShape` 并 `continue` 跳过（`cxx:162-167`），该组随后**取另一组的维数**（`cxx:203-209`）；② `FUSE` 要求两组 `iDimMax` 相等、且组内不自相混维（`cxx:181-186`/`:193-195`）⇒ 否则错误 `BOPAlgo_AlertBOPNotAllowed`；③ `CUT` 要求 `iDimMax(objects) <= iDimMin(tools)`（`cxx:194`）；④ `COMMON` 任意维数都允许。⇒ 端口里"face∩solid 按质心、solid−face 原样返回"的规则在 OCCT 中**不存在**，是自创几何。
- **改动（`crates/occt-topo/src/bop_builder_dispatch.rs`）**：① 删掉 `boolean_degenerate` + `degenerate_empty` + `degenerate_non_solid` + 只剩它用的 `shape_centroid`（共约 140 行）与 `use crate::brep_extrema::is_inside`；② `boolean_dispatch` 的非实体分支改调新的 `boolean_non_solid`：按 OCCT 把两个操作数**原样**交给 `bop_builder2::builder_bop_with_fuzzy`（= `BOPAlgo_BOP::Perform` + `BuildShape`），合法性由已移植的 `bop_bop::check_data`（`BOPAlgo_BOP::CheckData`）裁决，空实参补上 `BOPAlgo_AlertEmptyShape` 警告（`cxx:165`），错误按 OCCT 的 alert 名向上抛；③ 删掉那两个**为自创行为写的**单测（`boolean_degenerate_face`/`boolean_degenerate_empty`）；④ 更新 `boolean_multi` 折叠分支的过时注释（空实参由 `CheckData` 的维数回退自然吸收，不是特例）。
- **探针实测（临时 `examples/zz_probe_degen2.rs`，已删，`git grep zz_probe` = 0）**：`empty ∪ box → Compound[Solid] vol=1.000000`（空实参被跳过，结果就是 box）、`box − empty → box`、`empty − box → 空`、`empty ∩ box → 空`、`empty ∪ empty → Err BOPAlgo_Builder: too few arguments`（两组都空）；**非法组合按 OCCT 报错**：`face ∪ box → Err BOPAlgo_AlertBOPNotAllowed`（维数不等）、`box − face → Err BOPAlgo_AlertBOPNotAllowed`（3 > 2）；**合法且几何正确**：`face ∩ box → Compound[Face]`，`z=0.5` 的平面面 ∩ box 给出**真正裁剪出的面**（旧自创层只会按质心返回整张面或整个丢弃）。⇒ 这正是"删掉自创规则"的收益：以前静默产出错误几何的用例，现在要么得到 OCCT 的正确答案，要么得到 OCCT 的显式错误。
- **验证**：`occt-topo --lib` **1291/1**（唯一红仍是 T-01 `groove_cuts_cylinder`；总数 1293→1291 = 删掉的两个自创单测，属**减少自创测试**而非回归）、`phase3_integration` 4/4、`phase4_integration` 9/9、`phase6_integration` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05）、`export_data_obj` **16/16** 且 16 个模型顶点/面数与基线逐位一致。`occt-core`/`occt-geom` 未改动。
- **旁支观察（未改）**：新路径下非实体的 `FUSE`/`CUT` 会**返回 Err**（OCCT 语义），因此任何仍以"面 ∪ 实体"形式调用布尔的调用方必须像 OCCT 一样处理该错误；当前全部门禁语料都不走这条支路（否则会红）。另：`shape_volume` 作用在**单张面**上会给出无意义的数值（探针里 `mid-face ∩ box` 的面报 `vol=0.166667`）——开放壳上的体积公式，与本次改动无关，已登记不再追。

**批 25（A1/T-67 分步 2：`Extrema_ExtPS` 分派层忠实落地 —— 窗口/IsoIsDeg/`TreatSolution`/周期归一）—— 2026-09-20 第 49 轮**

- **对读的 OCCT 事实**（`src/ModelingData/TKGeomBase/Extrema/`，8.0.0 实测行数）：`Extrema_ExtPS.cxx` **428 行**（不是先前记的 376），`Extrema_GenExtPS.cxx` **1195 行**（不是 1056）。`Perform`（`cxx:269-367`）是**三引擎 switch**：五个初等面 → `myExtPElS`（`cxx:276-290`）；`SurfaceOfExtrusion` → `Extrema_ExtPExtS`（`cxx:292-317`）；`SurfaceOfRevolution` → `Extrema_ExtPRevS`（`cxx:319-343`）；`default:` → `Extrema_GenExtPS`（`cxx:345-356`），三者结果都过 `TreatSolution`。**初等面臂不接收窗口**（只传 `Precision::Confusion()`），窗口只由 `TreatSolution`（`cxx:97-135`：周期 `ElCLib::InPeriod` 归一 + 允许一个周期越界的裁剪面例外 + `[uinf∓tolu, usup±tolv]` 测试）施加；`Initialize`（`cxx:202-265`）做 ±1e10 无穷钳制、`nbU/nbV = 44`（B 样条/Bezier）否则 `32`、`IsoIsDeg`（`cxx:32-93`，11 个采样点、步长 < `PConfusion` 即返回 false）命中则 `=300`。另外 `IsDone` 来自**引擎**（初等臂里点落在轴/顶点上时 `myDone=false`，`ExtPElS.cxx:83`/`:228`/`:279`/`:299`），而 `NbExt` 是 `TreatSolution` 过滤**之后**的长度 ⇒ `done=true, NbExt=0` 是合法状态；`d11..d22/P11..P22`（`cxx:401-418`）在 8.0.0 里**从未被写入**（构造即 0，`cxx:148-152`）。
- **新增 `occt-geom/src/extrema_surf/point_surface_extrema.rs`（`ExtPs`）**：逐条落地上述内容——`ext_ps_surface_type`（`GetType()` 等价，来自 `gp_pln`/`gp_cylinder`/`gp_cone`/`gp_sphere`/`gp_torus`/`is_surface_of_linear_extrusion`/`is_surface_of_revolution`/`is_bspline_surface`/`is_bezier_surface`）、`iso_is_deg`、`initialize`（含 ±1e10 钳制与 `nbU/nbV`/300 规则，暴露 `sample_counts()`/`iso_degenerate()` 供探针取证）、`perform`（三引擎 switch + `TreatSolution`）、`is_done`/`nb_ext`/`square_distance`/`point`（引擎序，不排序）、`trimmed_square_distances`、`set_flag`/`set_algo`。**UNPORTED 就地标注**：`Extrema_ExtPExtS`/`Extrema_ExtPRevS`（`cxx:292-343`，两种面型现落通用臂）与 `Extrema_GenExtPS`（`cxx:346`，含 `math_FunctionSetRoot`/`GeomGridEval_Surface`/`Bnd_Sphere` UBTree，端口仍用 `extrema_surf/numeric_extrema.rs` 的 24×24 网格 + 数值 Jacobian 替代；`set_flag`/`set_algo` 因而只存值不影响结果）。
- **配套**：`occt-core/src/elib/clib.rs` 补 3D `ElCLib::InPeriod`（`ElCLib.cxx:95-111`，此前只有 2D 版）——`TreatSolution` 需要它；**接线**：`point_surface_extrema_all`（自然区间，`TolU=TolV=Precision::PConfusion()`，与 `ShapeAnalysis_Surface.cxx:1349` 一致）与 `point_surface_extrema_box`（显式窗口，含 `TreatSolution` 裁剪）都改走 `ExtPs`；窗口内无解时保留两级：先在同一窗口上跑通用替代件（OCCT 的 `ValueOfUV` 此时走 `SurfaceNewton`+`UVFromIso`，`ShapeAnalysis_Surface.cxx:1449-1459`），再退到自然区间的 `fallback_point_surface`（仍标 UNPORTED）。
- **探针实测（临时 `examples/zz_probe_extps.rs`，已删，`git grep zz_probe` = 0）**：柱面 `r=2`、点 `(3,0,5)` → `done=true, nb=2`，解 `u=0 → (2,0,5) d=1`、`u=π → (-2,0,5) d=5`（解析精确）；同柱面窗口 `u∈[1,2]` → **`done=true, nb=0`**（两个驻点被 `TreatSolution` 滤掉，正是 OCCT 语义）；窗口内两解都被滤掉时 `_box` 落回窗口内的替代件（`_box(plane,[5,6]×[0,10]) → (5.0,2.0) d=5.0`，未回到自然区间）；球面极点 `iso=(false,true) → counts=(32,300)`；平面 `(1,2,3) → u=1,v=2,d=3`。旁支观察（未改）：端口的**未裁剪**柱面 `v_range=(-inf,inf)` 会让 `IsoIsDeg` 的有限性守卫整体跳过而返回 true（OCCT 同一输入亦然，OCCT 侧之所以不触发是因为喂进来的是裁剪面的 `Adaptor3d_Surface`）⇒ `counts=(300,32)` 只是该差异的可见结果，对端口行为无影响（`nbU/nbV` 只被未移植的 `GenExtPS` 使用）。
- **验证（逐项等于基线）**：`occt-core --lib` 290/290、`occt-geom --lib` 151/151、`occt-topo --lib` **1293/1**（唯一红 = T-01）、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05）、`phase3_integration` 4/4、`phase4_integration` 9/9、`phase6_integration` 5/5、`export_data_obj` **16/16** 且十项顶点/面数与基线逐位一致（17745/22119、24/12、146/140、3343/4336、3105/4792、642/1244、1369/2592、3494/5078、104/92、600/790）。
- **未做（T-67 余项 → 分步 3）**：`Extrema_GenExtPS` 引擎（1195 行）+ `Extrema_ExtPExtS`/`Extrema_ExtPRevS`（+`math_FunctionSetRoot`，仓内未移植）。T-37 的 39–40 处调用点迁移**仍不做**——迁到 `extrema_surf` 只会让初等面变精确，一般曲面仍会落到同一个替代件，故按"先有忠实引擎再一次性迁"的次序继续挂起（`brep_surface::surface_closest_params` 仍保留）。

**批 24（A3/T-71：STEP 写侧曲线分派忠实化 —— 未识别曲线不伪造实体，B 样条/Bezier 直写）—— 2026-09-20 第 48 轮**

- **落地 1（`occt-topo/src/step/write_context.rs::emit_curve_entity`）**：删除 `self.splines` 的 B 样条**采样拟合**兜底与"起点切线 LINE"兜底，未识别曲线返回 `None`；`emit_edge` 相应写 `EDGE_CURVE('',#v1,#v2,$,.T.)`——对应 OCCT `GeomToStep_MakeCurve.cxx:100-103` 置 `done=false` 后 `TopoDSToStep_MakeStepEdge.cxx:345` 拿到的**空曲线句柄**（STEP 里即几何属性未设）。签名收窄为 `emit_curve_entity(&mut self, c: &dyn Curve) -> Option<usize>`（OCCT 的 `GeomToStep_MakeCurve(C)` 本就不带参数区间）。
- **落地 2（补齐 `Geom_BoundedCurve` 臂，否则会大面积写出 `$`）**：按 `MakeCurve.cxx:94-99` → `GeomToStep_MakeBoundedCurve.cxx:37-80`：B 样条直接 `write_bspline_curve`（有 weights 走 `rational`），**Bezier 先按 `GeomConvert::CurveToBSplineCurve` 转成夹持 B 样条**（degree = 极点数−1、结点 `0`/`1` 各 `degree+1` 重）再写；为此给 `Curve` trait 加 `bspline_weights()`（默认 `None`；`GeomBSplineCurve` 返回 `weights.as_deref()`，`ReparamCurve`/`IsoCurve`/两侧 `GeomTrimmedCurve`/`shhealing::transfer_params` 的包装类型转发）。**UNPORTED**：`Geom_BSplineCurve::SetNotPeriodic`（`cxx:46-51`）未移植 ⇒ 周期 B 样条保留周期性；Bezier 在本端口只有非有理表示（`GeomBezierCurve{poles}`），有理 Bezier 仍缺。
- **落地 3（退化边路由）**：端口为球极点这类**退化边**存了一个"点曲线"占位，而 OCCT 的 `BRep_Tool::Curve` 在那里返回**空**（`MakeStepEdge.cxx:194-196` 的分支）；现在 `emit_edge_curve` 先看退化标志（`BRepTool::is_degenerated` / 注册表 `is_degenerated_edge`），命中就走 `MakeStepEdge.cxx:263-330` 的 "edge without 3d curve; creating" 支。端口在该支仍以"两端点连线的 LINE"代表 OCCT 的"平面+直线 ⇒ `Geom_Line`；否则采样拟合 B 样条"，后者已就地标 `UNPORTED`。
- **顺带删除**死代码 `fit_bspline_curve`（8 点采样插值），即 T-62 余项里的"STEP 写侧 B 样条采样重拟"。
- **探针（临时，已删）**：`write_step` 后统计实体——`Shape-2.step` 82 条 `B_SPLINE_CURVE_WITH_KNOTS`、0 条 `$`；`Shape.step` 14 条 B 样条、0 条 `$`；`ATU01038.step` 163 条 B 样条 / 391 `LINE` / 337 `CIRCLE` / **2 条 `$`**（真正的 `done=false` 用例）。三份文件写回后 `read_step` 均无语法错误（解析 ok；**旁支观察**：`write_step(&model)` 的输出再读回得到 0 个 shape，属写/读表示层不对称的既有问题，与本批无关，已登记待另案）。
- **验证（全等于基线）**：`occt-core --lib` 290/290、`occt-geom` 151/151、`occt-geom2d` 72/72、`occt-topo --lib` **1293/1**（唯一红 = T-01）、`phase3_integration` **4/4**、`phase4_integration`（写侧）**9/9**、`phase6_integration` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05）、`export_data_obj` **16/16**。
- **过程中的一次回归与订正（留痕）**：只删兜底不加 `Geom_BoundedCurve` 臂时 `occt-topo --lib` 变 **1287/7**（6 个 STEP 写/读回环与 xcaf 回环红）；探针定位到两类受害者——① 球极点的"点曲线"占位（已按退化标志路由到 `MakeStepEdge` 的 creating 支解决）；② B 样条/Bezier 曲线（已按 `MakeBoundedCurve` 直写解决）。二者修好后回到 1293/1。

**批 23（A10/T-73 + 新 T-76：`choose_left_way` 3D/2D 最小夹角 + `GpDir2d::Angle` 有符号化）—— 2026-09-20 第 47 轮**

- **落地 1（3D，`occt-core/src/poly/make_loops.rs`）**：`MakeLoopsHelper` 补 `get_normal`/`get_first_tangent`/`get_last_tangent`（默认返回 `None`，即 OCCT 基类 `Helper` 返回 `false` 的语义：调用方回落到 `theLstIndS.First()`）；`choose_left_way` 逐行移植 `Poly_MakeLoops3D::chooseLeftWay`（`cxx:614-678`）——`Normal.CrossCrossed(TgtRef, Normal)` 投影出参考方向、逐候选投影后取 `gp_Dir::AngleWithRef` 的**最小有符号角**、`angle < 1e-4 - PI` 归为 `PI`、访问器失败/参考方向退化/全候选被拒三处都回落 `First()`。
- **落地 2（2D）**：新增 `PolyMakeLoops::new_2d(helper, left_way)`（`myRightWay = !theLeftWay`，`Poly_MakeLoops.cxx:682-687`）与 `two_d` 判别（端口把 `Poly_MakeLoops3D`/`Poly_MakeLoops2D` 合成一个结构体，已在文档注释说明），`choose_left_way` 按 flavour 分派到新的 `choose_left_way_2d`（逐行移植 `cxx:692-738`，含 `myRightWay` 取反）。
- **落地 3（新登记 T-76，同批修复）**：`GpDir2d::Angle` 原为 `acos(dot)`（**无符号** [0,π]），OCCT `gp_Dir2d::Angle`（`gp_Dir2d.cxx:26-63`）是**有符号** ]−π,π] 的 `acos`/`asin` 分段；现按 OCCT 重写，并把 `is_normal`/`is_parallel`/`is_opposite` 改成 `gp_Dir2d.hxx:393-430` 的 `abs(Angle)` 形式（对原语义等价，故门禁不变）。
- **探针验证（临时探针，已删）**：3D —— normal=+Z、候选 −60°/−30° ⇒ 选 −30° 候选（索引 3）；normal=−Z ⇒ 符号翻转 ⇒ 选索引 2；候选序 `[3,2]` 结果不变（证明不再是"取首候选"）；无 normal ⇒ 回落首元素。2D —— `theLeftWay=true` ⇒ 选 −60° 候选（索引 2）；`myRightWay=true` ⇒ 取反后选 +30° 候选（索引 3）。`GpDir2d::Angle` 输出 +30°/−30°/180°/−135°，`is_parallel(+170°)/is_opposite(−179°)/is_normal(90°)` 全 true。
- **验证（全等于基线）**：`occt-core --lib` **290/290**、`occt-math` 215/215、`occt-geom` 151/151、`occt-geom2d` 72/72、`occt-topo --lib` **1293/1**（唯一红 = T-01）、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05）。
- **影响面（如实）**：`PolyMakeLoops` 只被 `make_loops.rs` 自身测试使用（`git grep` 全仓仅本文件），`GpDir2d::angle` 的调用者只有三个谓词与 `GpLin2d::angle`（OCCT 的 `gp_Lin2d::Angle` 本身也是有符号），故本批为**潜在缺陷**忠实化：门禁数字不变，但一旦 2D/3D 的 `Poly` 拓扑（`poly` 分层、`BOPAlgo` 的 wire 重建）接线，行为将与 OCCT 一致。

**批 22（A11/T-47 子项 2：`BuildDraftSolid` 不再丢"顶点<3"的面）—— 2026-09-20 第 46 轮**

- **问题**：`bop_draft_solid_occt.rs::add_draft_face_once`（= `BOPAlgo_Builder::BuildDraftSolid` 的落地点）先用自创的 `face_is_degenerate(face) = vertices_of(face).len() < 3` 把面丢掉再 Add，而 OCCT 在 `BOPAlgo_Builder_3.cxx:329`（same-domain 镜像分支）、`:342`（非 SD 镜像分支）、`:357`（无镜像分支）**一律** `iFlag = 1; aBB.Add(aShD, aFx);`——没有任何顶点数/面积谓词。文件里那条注释也自述"OCCT `BRep_Builder::Add` would still add them"，属已登记的自创规则（audit A11）。
- **落地**：删除该谓词与 `face_is_degenerate` 辅助函数，`add_draft_face_once` 只保留 `as_face` 判定与 `same_tshape` 去重（去重不是几何规则：OCCT 的 `myImages` 列表本身不含重复 TShape，这里作为不变量校验保留并已在注释中说明）。
- **验证（全等于基线）**：`occt-topo --lib` **1293/1**（唯一红 = T-01）、`phase3_integration` **4/4**、`phase6_integration` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05）、`export_data_obj` **16/16** 且 16 个模型 v/f 与 §14.1 记录**逐位一致**（⇒ 删谓词在当前语料上不改变任何输出，是"潜在缺陷"的忠实化）。
- **旁支观察（未动手）**：`bop_builder_splitapi.rs:351-391` 的 `heal_tolerance` 还带一条**端口自创**的修复链（`weld_coincident_vertices` → `remove_small_edges(tol*8)` → `remove_degenerate_faces`，其中 `face_is_degenerate_tol` 也是"边数<3 或近零面积"），OCCT 对应的是 `ShapeFix_Shape`/`BRepBuilderAPI_Sewing` 的容差修形——与 A5/A11 同族，属另一处独立缺口，仅登记。

**批 21（A0/T-75：offset 曲线自身的 `EvalD3`/`EvalDN` 补齐，A0 收口）—— 2026-09-20 第 45 轮**

- **落地 1（3D，`crates/occt-geom/src/offset.rs`）**：新增 `calculate_d3` = `Geom_OffsetCurveUtils::CalculateD3`（`pxx:215-306`）的逐行转写：`Ndir/DNdir/D2Ndir/D3Ndir = d1..d4 ^ Direction`，`R2..R7`、`Dr/D2r/D3r = Ndir·D3Ndir + 3·DNdir·D2Ndir`，`R7 <= gp::Resolution()` 的**不稳定分支**（内层 `R6` 判据、`R4 = R2²`）与 `else` 的 IICURV **稳定分支**，尾部 `isDirChange ⇒ D3.Reverse()`（`pxx:301-305`）。补 `GeomOffsetCurve::d3`（`EvalD3` `cxx:342-380` → `EvaluateD3` `pxx:495-529`：基曲线 `EvalDN(u,4)` 作 `theD4`，基曲线 `D1` 奇异时走 `AdjustDerivative(..., 4, ...)`）与 `eval_dn`（`cxx:386-410`：1/2/3 阶取自 `EvalD1/2/3`，**更高阶转发基曲线** `cxx:409`）。
- **落地 2（2D，`crates/occt-geom2d/src/offset.rs`）**：同构的三件——`calculate_d3`（`pxx:189-279`，法向 `(D.y, -D.x)`，`GP_RESOLUTION` 判据）、`Geom2dOffsetCurve::d3`（`cxx:289-327` → `EvaluateD3` `pxx:460-484`）、`eval_dn`（`cxx:332-356`）。
- **解析验证（临时探针，已删）**：基曲线 = 半径 2 的 XY 圆、偏移 0.5（方向 +Z）⇒ 偏移曲线应是半径 2.5 的同心圆；实测 `off.d3(u)` 与解析 r=2.5 圆的 `D0/D1/D2/D3` 在 `u = 0, 0.7, 2.5, 5.0` 上逐分量差 **最大 1.11e-16**，`eval_dn(4)` 与基曲线 `eval_dn(4)` 完全一致（差 0）。
- **验证（全等于基线）**：`occt-geom --lib` **151/151**、`occt-geom2d --lib` **72/72**、`occt-topo --lib` **1293/1**（唯一红 = T-01）、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05）；五 crate `--all-targets` 编译 exit 0。
- **A0 结项口径（如实）**：offset 曲线的 D0–D3 与 `EvalDN` 已全部忠实；**遗留 UNPORTED** 是非初等/非 B 样条基曲线（如 `Geom_BezierCurve::EvalDN`，`Geom_BezierCurve.cxx:601-617`）在 4/5 阶上仍是 trait 默认零值 —— 两处 offset 文件头已如此登记。

**批 20（A2/T-38：点–实体分类改走忠实 SolidClassifier，并补上 `SolidExplorer::myMapEV`）—— 2026-09-20 第 44 轮**

- **落地 1（`occt-topo/src/brep_extrema.rs::is_inside`）**：Solid/CompSolid 改调 `crate::brep_class3d::SolidClassifier::classify(shape, p, CONFUSION)`（= `BRepClass3d_SClassifier::Perform` 的移植），`In` 即"在体内"；原来的 **7×7 UV 网格 + 抖动 +X 射线奇偶**实现改名为 `is_inside_mesh`，只留给**非实体**形状（`BRepClass3d_SolidClassifier` 只对实体/闭壳有定义，OCCT 对 shell 走 `BRepClass3d_SolidExplorer`，端口未接线）并就地标 `UNPORTED`。严格内部语义不变（`On` 不算 In）。
- **落地 2（同批发现的真实缺陷，`occt-topo/src/brep_class3d.rs`）**：端口 `SClassifier::perform` 的 ON 前置判定用 `on_vertex_or_edge(expl.shape(), …)`，即**整个形状的全部顶点/边**；OCCT 用的是 `myMapEV`——`BRepClass3d_SolidExplorer::Init`（`cxx:930-982`）只把**非 INTERNAL/EXTERNAL 的面**上的**非 INTERNAL/EXTERNAL、非退化边及其顶点**收进 BVH 树（`cxx:217-227` 用它做 ON 判定）。差别实测：对"单位盒 + 一个位于 (0.5,0.5,0.5) 的内部孤立顶点"，旧实现把查询点 (0.5,0.5,0.5) 判成 **`On`**（⇒ `bop_build_common::tests::inside_vertex_settles_into_original_solid_as_copy` 回归），OCCT 语义应为 **`In`**。现补 `SolidExplorer::edge_vertex_map`（逐面→逐 wire→逐边，按上述过滤）并让 ON 判定只用该列表。
- **验证（全部等于或优于基线）**：`occt-topo --lib` **1293/1**（唯一红 = T-01 `brepfeat::tests::groove_cuts_cylinder`，与基线同）、`phase3_integration` **3/4 → 4/4**（`primitives_measure_correctly` 转绿）、`phase6_integration` 5/5、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05，报错文本仍逐字 `Offset divergence 2208.0 vs shape_volume 1612.9`）、`export_data_obj` **16/16** 且 16 个模型的 v/f 与 §14.1 记录**逐位一致**；探针（`examples/zz_probe_a2.rs`）已删除，工作树只含上述两个文件的行为改动。
- **未完成（如实登记）**：非实体（shell/compound）的 `is_inside` 仍是端口自创的网格射线奇偶（`is_inside_mesh`，已标 UNPORTED + OCCT `BRepClass3d_SolidExplorer` 出处）；A4（`algo_tools::compute_state` 的 32×32 投影）与 A16（256×256 采样求交）仍需 T-67 **分步 3** 的 `Extrema_GenExtPS` 引擎（分派层与初等面解析臂已在第 49 轮落地）。

**T-69 第 3 轮（2026-09-20 第 43 轮 goal round）：判据性实验否证"跨 wire 周期对齐"假设；face 20 的几何/分片参数与成功的 face 8 对比记录在案**

- **实验（临时探针，已删）**：在 `DelaunayNodeInsertionMeshAlgo::perform` 开头对周期面（`gp_torus()` 且 `is_u_periodic()`）做"把每条 wire 的 pcurve 整链按 `k·2π` 平移到第一条 wire 的 u 窗口"的归一，然后数 `triangulate_model_faces` 的 `mapping_failed`（= 忠实路径 0 三角的面）：**基线 162 ⇒ 归一后 161**，且 **face 20 仍是 `tris=0 nodes=74`**。⇒ 第 2 轮提出的"两 wire 的 u 窗口相差一个周期导致 face 20 失败"**被否证**（对齐后仍失败），u 窗口差异只是伴生现象，不是成因。
- **face 20 与 face 8 的实测参数（探针，已删）**：
  - face 20：Torus **major=1.699445、minor=1.5（非 spindle，major>minor）、u_periodic=true**；splitter `range_u=(−4.7124, 1.5708)`（诚实周期钳制后的窗口）、`range_v=(0.05088, 1.5708)`、`delta=(0.39482, 0.66683)`、`tol_uv=(6.283e−7, 1.520e−7)`；两条 wire 各 2 条边（`FixLacking` 复制的闭合圆），v 分别 0.05088 与 π/2。
  - face 8（忠实成功，74v/72t）：Torus major=15.35、minor=5.5；`range_u=(1.5708, 7.8540)`、`range_v=(2.65050, 3.63269)`、`delta=(0.09765, 0.18184)`。
  ⇒ 两面的链形状/边数/v 结构同类，差别只落在**具体参数值**（fat torus 的 minor/major 比 0.88、v 带 `[0.05, π/2]`、u 窗口偏移）与 splitter 派生量（`delta`、`tol_uv`、`cells_count`）上；第 2 轮的"u 窗口"已排除，故下一步须**直接对这张面对拍 Delaunay 内部**。
- **本轮的另一条负结果（重要，避免重走）**：`build_shape_mesh` 开始时逐面扫"边界边的离散点是否经过 `(−0,−11.823,−424.742)`"，**0 命中**（`edge.discretization()` 在网格化前为空，且即便在网格化后也 0 命中）⇒ 该参考顶点**不是**端口任何面的边界顶点，而是 OCCT 在某张面内部放下的网格顶点；"共享该顶点的两个面"这一推断作废（参考里 v962/v965 重合可能来自两张**共面/贴面**的面各自网格化）。
- **验证**：本轮**未改行为代码**（探针全部 `git checkout` 还原 + `git grep dbg-t69/DSH_T69_ALIGN` 于 `crates/` 0 残留）；`occt-topo --lib` **1293/1**。
- **下一步（第 44 轮，聚焦单面对拍）**：把第 9 轮那套 Delaunay 分支计数/链 dump 仪器**只对 face 8（成功）与 face 20（失败）开启**：`init_data_structure` 的链节点 UV、`mesh_polygon` 的入参数多边形（长度/面积/首尾节点）、`decompose_simple_polygon` 的 `NO_EAR` 原因（`skip_prec`/`skip_neg`/`filter`/`isect`）与 `Glued/Same/PointOnSegment/Cross` 事件、以及 `cells_count`/`set_cell_size`/`set_tolerance` 的取值，逐项对照 OCCT `BRepMesh_Delaun.cxx`，找出**成功面与失败面在这条控制流上的第一处分叉**（第 9 轮只证明了"T0M 全局无自创分支"，未针对这对同类面做过差分）。

**T-68 第九轮（2026-09-20）：Delaunay 层逐行核对全部忠实；`NO_EAR` 清空来自"拼接多边形 + 真实重叠共线前沿链"，根因在 Delaunay 之外 ⇒ 立项 T-69；A25 的"删除补丁致清空"结论被实测推翻**

- **分支级计数（临时探针，已删）**：T0M 上 `decompose_simple_polygon` 共 **975 次 `NO_EAR`**（`used_link_id == 0` → `thePolygon.Clear()`），而这些失败里 **`skip_filter = 0`、`isect = 0`**——没有任何候选是被"角度过滤"或"相交测试"拒掉的；**全部**落在 `anAbsDist < Precision`（`skip_prec`）或 `aDist < 0.`（`skip_neg`）。两类：
  - **负面积（CW）**：`poly_area = -1032.496` / `-1.7226`，`skip_neg=3`（三个 pivot 的叉积全负）；
  - **零面积退化**：`poly_area = 0.0`，`skip_prec=6..7`，所有 pivot 的 `dist = -0.0`、`angle = π`（全部落在参考边反向射线上）。
- **链一致性检查（决定性）**：`NO_EAR#2..#5` 的多边形起于节点 84、止于 42（首≠尾）；**`NO_EAR#200` 的入参是闭合的 26 链**（`incoming_chain = (20,19),(19,18),(18,17)…`，`self_closed=true`），修正后变成 `[19, 26, 25, 24, …, 20]`，链在索引 1 处断成 `(20,19) | (1,26)`（`breaks=[(1,19,1)]`）⇒ 修正循环把一个**不成环的拼接多边形**交给了分解器，分解器"判无耳并清空"是**正确的**。
- **修正循环事件跟踪**：所有命中都是 **`Glued`**，且**几何上真实**：命中两端 4 点共线且线段重叠——例 `(20.3780056797, 25.6094339642)→(19.4517326942, 25.6094339642)` 对 `(19.5024277029, 25.6094339642)→(20.4350375644, 25.6094339642)`，两对端点各差 ~0.05、**v 完全相同**。
- **重复点普查**：这些面 `exact_pairs = 0`（无同位置节点，`add_node` 的 `index_of_node` 合并正常）⇒ 不是"重复插点"，而是**前沿链本身沿同一条直线重叠**。
- **逐行核对（本轮，全部一致，行号取 OCCT 原文）**：`meshPolygon`(1818-2076)↔`mesh_polygon`、`processLoop`(1756-1776)↔`process_loop`（`Prepend` 倒序遍历后仍是升序，与端口 `polygon[link_from + i]` 的 1-based 等价性逐项核对）、`createAndReplacePolygonLink`(1783-1814)、`decomposeSimplePolygon` 头/耳循环/尾(2120-2314)、`getOrientedNodes`(1735-1749)、`checkIntersection`(1324-1372)、`findNextPolygonLink`(1206-1316)、`meshLeftPolygonOf`(1066-1196)、`classifyPoint`(516-552)、`IntLinLin`/`IntSegSeg`(302-461)、`AddLink`/`SubstituteLink`/`RemoveLink`(71-145，含 `myDelLinks` FIFO 复用)、`cleanupPolygon`(1404-1529)、`AngDeviation90Deg = π/2`(cxx:39-40)、`Precision = PConfusion`(cxx:43)。**Delaunay 层没有找到自创/偏差。**
- **推翻上一轮结论（重要）**：第八轮把 T0M 缺口归因于 **A25 的自创"删邻三角形"补丁**（`meshing/delaun/polygon_meshing.rs:346-375`）。本轮实测 `DELETE_NEIGHBOUR` 事件 **0 次**（该分支在 T0M 上从不触发），而 975 次清空全部由 `skip_prec`/`skip_neg` 触发 ⇒ **删除补丁不是 T0M 缺口的成因**。A25 仍是应修的自创项（按原卡处理），改它**不会**修好 T0M（这点直接决定 T-61 不再是 A13/A18/T-55/T-59 的前置）。
  - **订正（批 12，2026-09-20）**：该"0 次"是**采样口径**的结论（只打印前若干次与每 200 次 `NO_EAR`）。批 12 删掉补丁后 `step_obj_parity` 的 T0M 网格由 46514/46945 变 **46516/46962**（v+2/f+17，其余模型逐位不变）⇒ 补丁在 T0M 上**确实触发过**；但"它不是 166 面缺口的成因"这一结论不变（缺口由 `NO_EAR` 耳选择失败造成）。
- **T0M 缺口的直接机制（已闭环到可修点）**：带孔面（`wires=2`）的前沿链里有**重叠共线链接** → 修正循环 `Glued` 大批删段 → 交出拼接/开放多边形 → 分解器正确清空 → 该面 0 三角（166 面）。**来源在 Delaunay 之外**，故立项 **T-69**（`node_insertion.rs`/`model_builder`/`shape_tool` 的重复/偏移插入）。
- **未决**：无 OCCT 运行时可对拍（`DRAWEXE` 因缺 DLL 无法启动，见 §9），故"OCCT 在同输入下是否也会 `Glued` 清空"无法实测；但本轮证明端口在这条路径上与 `.cxx` 逐行一致，且输入（重叠前沿链）本身可疑。
- **验证**：本轮**未改行为代码**；`cargo check` exit 0；插桩与探针（`zz_probe_t0m6.rs`、`meshing/delaun/polygon_meshing.rs` 全部 `[dbg-*]`）已 `git checkout` 还原 + `git grep` 复核 0 残留，工作树干净。**T-55 的逐面 UV 栅格回退仍然承重**，继续阻塞到 T-69 落地。

**批 2（T-56 / A20 + A22 + 新 A29）已完成 —— 2026-09-20**

- **落地**：
  1. `crates/occt-core/src/elib/clib.rs`：新增 `ellipse_parameter` / `hyperbola_parameter` / `parabola_parameter`（逐行对 `ElCLib.cxx:1226-1272`）与 `parameter_elips` / `parameter_hypr` / `parameter_parab`（`ElCLib.lxx:335-351`）。
  2. `crates/occt-geom/src/curve.rs` + `line.rs` + `hyperbola.rs` + `parabola.rs`：补 `Curve::{gp_line,gp_hyperbola,gp_parabola}`（`Adaptor3d_Curve::Line/Hyperbola/Parabola`）。
  3. `crates/occt-topo/src/step/read_geometry.rs`：`shape_analysis_project_act` 的 `!ok` 分支改为 `ShapeAnalysis_Curve.cxx:355-477` 的忠实 switch（Circle→Hyperbola→Parabola→Line→Ellipse，`default:` 才分段搜索）；`Project`/`ProjectAct` 去掉 `Option` 假失败通道。
  4. `crates/occt-topo/src/step/read_geometry.rs::edge_params_for_curve`：**删除全部三个自创回退**（端点距捷径 `p1.distance(p2) < 1e-3`、脱靶 `(a.0-b.0).abs() <= PConfusion` 退回整条结点域、`classify_curve` 采样族——含三点外心 `circle_edge_params`/`circle_params_from_circ` 共 −127 行），改为 `TranslateEdge.cxx:442-444` 的两次 `Project`。
  5. `crates/occt-topo/src/step/read_topology.rs`：补 `TranslateEdge.cxx:437-441` 的 `V1.IsSame(V2)` → `GetCartesianPoints`。**注**：端口里 `IsSame ⟺ 同一实体引用`，两点必然相同 ⇒ 该分支当前不可达，属潜伏正确性（同 A0 的性质），已在代码里写明。
  6. `crates/occt-geom/src/extrema_pc/point_curve.rs`（**根因修复，A29**）：直线解析臂原先把 `Extrema_ExtPElC` 的解析 `gp_Lin` 建在 `d0(uinf)` 上，参数因此整体平移 `-uinf`；改为取曲线自身的 `gp_line()`（OCCT `Extrema_GGExtPC` 用 `theCurve.Line()`），无 `gp_line` 时回退 `d0(0)` 原点重建。
- **根因链（本轮最重的一段，可复现）**：删掉自创回退后 `step_obj_area` 由 11/11 → 10/11，`data/occ/OffsetPlaneHoleEdge.step` 面积 **280.00 → 204.00**；`step_geometry_parity` 1/3、`step_to_obj` 12/13。临时插桩（已删、grep 复核 0 残留）逐步定位：
  1. 该模型 6 条边全是 `Geom_Line`，忠实臂给出 `(2, 11)`/`(2, 3)`/`(2, 5)`，而 `ElCLib::LineParameter` 给出真值 `(0, 10)`/`(0, 2)`/`(0, 4)` ⇒ `d(w1)-p1 = 2.0`；
  2. `shape_analysis_project_act` 的 extrema 臂直接返回 `Some((2.0, 0.0))`——**参数与距离不自洽**（`d0(2)` 距 p1 为 2）；
  3. `interval_perform`/`perform_general` 从未被调用 ⇒ 走的是解析臂 `ext_pelc_all`；
  4. 偏移量恰为**窗口下界**：上游对无界曲线用 `t_est ± (|t_est|+1)` 窗口（−2 或 −1），而 `ext_pelc_all` 用 `loc = d0(uinf)` 重建 `gp_Lin`，于是 `line_all` 的参数 = 真参数 − `uinf` ✓ 与观测的 +2/+1 完全吻合。
- **验证（全绿，与基线逐项一致）**：`occt-topo --lib` **1293/1**、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**（`OffsetPlaneHoleEdge` 回到 280.00）、`step_geometry_parity` **2/3**；`occt-geom` 151/151、`occt-core` 290/290（1 ignored）、`occt-geom2d` 72/72；`cargo check` 四个 crate exit 0；插桩与探针 0 残留。
- **旁支**：`classify_curve`/`CurveKind` 仍被 STEP **写**侧使用（`step/format.rs:670`、`step/write_context.rs:157`），其删除归 **A3/T-39**；`step/read_geometry.rs` 内已无 `classify_curve` 调用。

**批 3（A3/T-39：`classify_curve` → `IsKind` 分派）已完成 —— 2026-09-20**

- **落地**：`crates/occt-topo/src/step/format.rs`、`crates/occt-topo/src/step/write_context.rs`。
  1. **删除** `CurveKind` 与 `classify_curve`（6 点 `|d²|` 采样 + `(max−min)/max < 0.02` 自创阈值）及只为其服务的 `midpoint`。
  2. `step/write_context.rs::emit_curve_entity` 改为 `GeomToStep_MakeCurve.cxx:50-104` 的 `IsKind` 顺序：`gp_line` → `gp_circ` → `gp_ellipse` → `gp_hyperbola` → `gp_parabola` → trimmed 基曲线递归（`cxx:66-92`）→（登记为 UNPORTED 的兜底臂）。
  3. `write_conic_params` 改为 `GeomToStep_MakeConic.cxx` 的顺序：Circle → Ellipse → Hyperbola → Parabola（`Geom_Line`/有界曲线由调用方处理，故直线返回 `None`，对应 `MakeCurve.cxx:54-59/94-99`）。
  4. 四个 `emit_*_entity` 由**采样重建**改为 OCCT 精确值：`GeomToStep_MakeCircle/Ellipse/Parabola/Hyperbola` 用 `Position()` 建 `AXIS2_PLACEMENT_3D`，半径/半轴/焦距直接取 `Radius/MajorRadius/MinorRadius/Focal`（此前的三点外心、`|d²|` 反推焦距、`d0(±1)` 中点全部删除）。
- **保留的两个 UNPORTED（新立项 T-71）**：`MakeCurve.cxx:100-103` 的 `done = false`（不写实体）未移植，端口仍写 B-spline 拟合或起点切线；trimmed 的 BSpline/Bezier 基曲线走端口自己的重映射结点（等价 OCCT `Segment`），已在代码处写明。
- **验证（全绿，与基线逐项一致）**：`step_to_obj` **13/13**、`step_obj_parity` **14/14**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**、`occt-topo --lib` **1293/1**；`cargo check` exit 0；现有 `step::tests::circle_ellipse_params`（断言 `CIRCLE … ,2.0)` 与 `ELLIPSE … 3.0,1.5)`）在精确实现下仍通过。
- **旁支**：`step` 模块内仅剩 `classify_surface`（`brep_surface`）这一处采样分类器 = **A24/T-60**，属下一批（面型分派）。

**批 4（A24/T-60：`classify_surface` → `GetType()` 分派）已完成 —— 2026-09-20**

- **落地**：`crates/occt-topo/src/meshing/range_splitter/param_set.rs`、`crates/occt-geom/src/{surface.rs,bezier_surface.rs}`。
  1. `classify_surface` 改为 `GeomAdaptor_Surface::Load`（`GeomAdaptor_Surface.cxx:422-513`）的逐行转写：**`DynamicType` 精确类判定**，顺序 RTS（递归基面）→ Plane → Cylinder → Cone → Sphere → Torus → SurfaceOfRevolution → SurfaceOfLinearExtrusion → **Bezier**（`cxx:480`，先于 BSpline）→ BSpline → Offset → `OtherSurface`。
  2. **删除**周期标志 `match (up, vp)`（`(:110-123)`，其中 `(true,false) if v_fin => Sphere` 会把任何 V 有限的柱/锥判成 Sphere）与自创判别式 `is_cylinder_like`（`(:128-133)`，两点半径采样 + `1e-7` 相对阈值）。
  3. 新增 `Surface::is_bezier_surface()`（`Adaptor3d_Surface::GetType() == GeomAbs_BezierSurface`），`GeomBezierSurface` 置 `true`——此前 `SurfaceType::BezierSurface` 不可达。
- **实测效应（网格密度，`step_obj_parity --nocapture` 对拍）**：Shape-2、ATU01038、T0M、`occ/bottom`、`occ/top` 均有变化（例 T0M **46647/46754 → 46514/46945**、ATU01038 17767/22160 → 17745/22119），其余 11 个模型逐位不变 ⇒ 证实 A24 描述的"错判会换掉整张内部节点网格"确实在 live 路径上发生。
- **验证（全绿，与基线逐项一致）**：`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**、`occt-topo --lib` **1293/1**、`occt-geom --lib` **151/151**；`cargo check` exit 0。
- **旁支（已立项 T-72/A30）**：`brep_surface::classify_surface`（`brep_surface.rs:78-101`）是**另一处**采样分类器（8×8 采样 + `1e-6`/`1e-4` 阈值判 Plane/Sphere/Other），服务于 STEP **读入**面型；本轮只报告未动手。

**批 5（A30/T-72：`brep_surface::classify_surface` → 精确 `GetType()`）已完成 —— 2026-09-20**

- **落地**：`crates/occt-topo/src/brep_surface.rs`。
  - `classify_surface` 改为 `GeomAdaptor_Surface::Load`（`GeomAdaptor_Surface.cxx:422-513`）的精确类判定：`rectangular_trimmed_basis` 递归 → `gp_pln` → `gp_cylinder` → `gp_cone` → `gp_sphere` → `gp_torus` → `Other`。
  - **删除** 8×8 采样分类（`is_planar(1e-6)`、等距球心判定、`1e-4*r` 阈值）；该实现**只能**报 Plane/Sphere/Other，柱/锥/环一律 `Other`。
- **可观测修复**：`step/read_topology.rs::unit_conversion_2d`（对应 `GeomConvert_Units.cxx:191-227`，OCCT 用 `IsKind` 分派）此前对柱/锥面拿到 `Other` ⇒ 走 `return c2d`，**整段单位换算被跳过**；现在柱面取 `(AngleFact, LengthFact)`、锥面取 `(AngleFact, LengthFact/cos α)`（锥面本就由 `cone_ref()` 先命中）。本轮门禁数字不变是因为测试模型用单位长度因子。
- **副作用（如实登记，含纪律说明）**：精确化后 `--lib` 由 1293/1 变 **1285/9**（8 个新增红）。逐个追因后确认**全部**来自"断言分类器看不见非 plane/sphere 面"这一**缺陷行为**，不是真实回归：
  1. 3 处**查找面的谓词**把 `SurfaceKind::Other` 当作"环面/锥面"——`fillet_curved/tests/helpers.rs::find_blend_face`（`:279-288`）、`tests/cone_blends.rs::cone_extraction_geometric`（`:28`）、`tests/cone_blends.rs::plane_cone_blend`（`:54`）⇒ 改为 `Torus`（前两处为锥面的改为 `Cone`），语义不变（"找到那张环面/锥面"）且更精确。
  2. 1 处**断言缺陷本身**——`fillet_edge/tests.rs:100-102` 注释原文 `// The vanilla classifier cannot see cylinders; the extended one can.` 并 `assert_eq!(vanilla, SurfaceKind::Other)` ⇒ 改为断言正确事实 `classify_surface(...) == SurfaceKind::Cylinder`。
  - **声明**：本轮**未新增任何测试**、**未放宽任何门禁断言**，只订正了 4 处"把缺陷写进测试"的位置；这与"为对齐新写单元测试"是两回事，故在此与 `_index.md` A30 行双重登记，便于后人复核。
- **验证（全绿，与基线逐项一致）**：`occt-topo --lib` **1293/1**（唯一红 = T-01 `brepfeat::groove_cuts_cylinder`）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3、`occt-geom --lib` 151/151；`cargo check` exit 0。
- **旁支（未动手）**：`fillet_edge::classify_surface_full`（`fillet_edge/chain.rs:245-255`）与 `fillet_curved::classify_surface_analytic`（`fillet_curved/surface_info.rs:138-148`）仍各自带一层 6×6 采样柱面/锥面检测，现在对真实柱/锥已成冗余（`classify_surface` 直接给出 `Cylinder`/`Cone`）；`gprop_analytic/analytic_props.rs:47` 还有一份同名实现。⇒ 建议立项（下一轮可选）：把这三处收敛为直接使用精确判定，删掉采样臂。

**批 6（A7/T-65：`Extrema_ExtPElC` 椭圆/双曲/抛物臂接线 + `myIsMin`）已完成 —— 2026-09-20**

- **落地**：`crates/occt-geom/src/extrema_pc/point_curve.rs`、`crates/occt-geom/src/trimmed.rs`。
  1. `ext_pelc_all` 补上三条解析臂（此前只接 Line/Circle，椭圆/双曲/抛物落 `Extrema_GGExtPC` 的数值 `default:` 臂）：
     - 椭圆 → `ellipse_all` + `myIsMin = sqDist(Us) < |P − C(Us + 0.1)|²`（`Extrema_ExtPElC.cxx:270-279`）；
     - 双曲 → `hyperbola_all` + `myIsMin` 用步长 `+1`（`cxx:377-384`）；
     - 抛物 → `parabola_all` + `myIsMin` 用步长 `+1`（`cxx:469-476`）。
  2. `GeomTrimmedCurve` 补 `gp_line`/`gp_ellipse`/`gp_hyperbola`/`gp_parabola` 转发（`GeomAdaptor_Curve::load` `cxx:252-254` 会解包 trimmed 并保留基曲线，故适配器的 `GetType()`/`Line()`/… 都是基曲线的）——此前只转发了 `gp_circ`，导致 trimmed 椭圆/双曲/抛物拿不到类型标签。
- **验证（全绿，与基线逐项一致）**：`occt-geom --lib` **151/151**（其中 `extrema_pc/tests.rs` 的椭圆/双曲/抛物极值断言在**解析臂**下仍通过 ⇒ 解析臂与数值路径在这些用例上一致）、`occt-topo --lib` 1293/1（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3；网格密度与批 4 后逐位相同（ATU01038 17745/22119、Shape-2 3105/4792）⇒ 本批是**分派忠实化**，在这些模型上数值不敏感。
- **A7 剩余**：仅 `extrema_cc` 侧（曲线–曲线种子集）⇒ **T-66**。

**批 7（A10/T-46：`occt-core` 静默默认值 4 处）已完成 3/4 —— 2026-09-20**

- **落地**：`crates/occt-core/src/elib/surface_eval.rs`、`bnd/bsphere.rs`、`elib/intersect.rs`、`poly/make_loops.rs`。
  1. `surface_d2`：球/环面的 `duu/dvv/duv` 由**恒为 0** 改为 `ElSLib::SphereD2`（`cxx:975-1037`）/`TorusD2`（`cxx:1039-1100`）的公式（以 `Vxy = cosU·X + sinU·Y`、`DVxy = -sinU·X + cosU·Y` 表示）：球 `Vuu=-R·cosV·Vxy`、`Vvv=-R·cosV·Vxy-R·sinV·Z`、`Vuv=-R·sinV·DVxy`；环面 `Vuu=-R·Vxy`、`Vvv=-r·cosV·Vxy-r·sinV·Z`、`Vuv=-r·sinV·DVxy`（与端口既有 `sphere_d1`/`torus_d1` 的参数约定一致）。
  2. `BndSphere::add_sphere`：按 `Bnd_Sphere.cxx:73-101` 补两个分支——**被包含 ⇒ 整体替换（圆心也换）**、包含对方 ⇒ 忽略；原实现只放大半径、圆心不动（偏心/偏小会漏检）。
  3. `BndSphere::distance`：按 `Bnd_Sphere.cxx:63-66` 改为**到球心**的距离（原名下实现的是到球面，同名不同义）；面距单独提供 `distances`（`cxx:45-50`）与 `square_distance`（`cxx:68-71`）。
  4. `circle_plane_intersection`：**删除**共面时凭空造的 `[center, C(π/2)]` 两点（OCCT 从不产生），改为返回空（共面 ⇒ 整圆、无孤立交点），整个函数标 `UNPORTED` 并写明忠实件是 `IntAna_Quadric`/`GeomAPI_IntCS`（`ElCLib` 无求交函数）。
- **未完成 ⇒ T-73（已就地标 UNPORTED）**：`poly/make_loops.rs::choose_left_way` 仍取首候选。OCCT 的最小夹角选择在 `Poly_MakeLoops.cxx:611-676`（3D，用 `myHelper->GetNormal`/`GetLastTangent`）与 `:688-700`（2D，加 `myRightWay`），仅在取不到法向/切线时才 `return theLstIndS.First()`；端口既无 helper 的法向/切线访问器也无 `myRightWay`，故永远走该兜底分支。该路径目前仅被自身测试调用。
- **验证（全绿，与基线逐项一致）**：`occt-core --lib` **290/290**（1 ignored）、`occt-geom --lib` 151/151、`occt-topo --lib` **1293/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3；`cargo check` exit 0。本批三处均为**潜伏**修复（`surface_d2`/`bsphere`/`intersect` 在仓内均无调用者），故门禁数字不变。

**批 8（A14/T-50 的高风险项：`validate.rs` 的"恒 true"假门禁）已完成声明 —— 2026-09-20**

- **落地**：`crates/occt-topo/src/validate.rs`。
  - 模块头新增 `UNPORTED` 段：说明本文件**不是** `BRepCheck_Analyzer` 的翻译——OCCT 的 `IsValid(S)`（`BRepCheck_Analyzer.cxx:458-478`）遍历 `myMap` 里每个子形状的 `BRepCheck_Status`，任一 `!= BRepCheck_NoError` 即 `false` 并递归；本文件的检查（Euler 特征、计数合理性、朝向统计）是本地不变量，**不等价**。
  - `Analyzer::is_valid()` 与 `check_geometry()` 加显式告警：前者**恒返回 `true`**、后者**恒返回 `Report::ok()`**，**不得作为校验门禁**使用。
  - **无行为改动**（纯声明）：`occt-topo --lib` 仍 **1293/1**（唯一红 = T-01）。
- **未完成（仍在 T-50 内）**：A14 其余包（`geom/`（csg/delaunay/triangulate/fit*/polygon_*）、`elib/measure.rs`、`hlr.rs`、`viz_scene/`、`draw/`、`xcaf/`、`render_svg.rs`）的模块头声明。

**批 9（A26/T-62 的 3 个可收口小项）已完成 —— 2026-09-20**

- **落地**（三处逐行对 OCCT 原文）：
  1. `crates/occt-core/src/io/ply.rs:99`：写侧属性类型 `property list uchar int vertex_indices` → **`uchar uint`**（`RWPly_PlyWriterContext.cxx:214`）。读侧解析器未动（其测试夹具 `:178` 仍用旧拼写，解析器两者都能读）。
  2. `crates/occt-core/src/io/stl.rs::compute_normal`：退化判据由自创的 `|cross| < 1e-12` 改为 OCCT 的**平方量**比较 `sq > gp::Resolution()`（`RWStl.cxx:325`、`:407`；`gp::Resolution() == RealSmall() == DBL_MIN`，端口用 `crate::precision::REAL_SMALL`）。
  3. `crates/occt-topo/src/vrml.rs:92`：`solid TRUE` → **`solid FALSE`**（`VrmlData_ShapeConvert.cxx:356-361` 构造 `VrmlData_IndexedFaceSet(..., IsCCW=true, IsSolid=false, IsConvex=false)`）。
- **验证（等于基线）**：`occt-core --lib` **290/290**（PLY/STL 的读写测试全过，说明改属性类型与阈值未破坏解析/往返）、`occt-topo --lib` **1293/1**（唯一红 = T-01，VRML 用例通过）、四道 STEP 门禁 14/14、13/13、11/11、2/3；`cargo check` exit 0。
- **未修（仍在 T-62）**：PLY `weld_vertices(1e-9)` 与"每面重复自身 node+偏移"（OCCT 每面重复 ⇒ 盒体 24 vs 8 顶点，改它要动 `PlyMesh` 结构）、STL 80 字节头与格式嗅探、IGES 曲线族 `<2%` 采样与整球自造 120 回转面、STEP 写侧 B-spline `n=8`/6×6 采样重拟、OBJ 写侧恒空 `vn`。

**批 9 旁的侦察结论（T-44/A8：`comp_curve_to_bspline`）——本轮未动手，已记录前置**

- `occt-geom/src/convert_bspl.rs::comp_curve_to_bspline` 现为**采样 0.1 折线 + `resample_bspline(pts, 1)`**（`:(326)` 的 `((b-a)/0.1).ceil().clamp(2,64)` 与 `:349`），而 OCCT `GeomConvert_CompCurveToBSplineCurve::Add`（`cxx:135-215`）是 **`IncreaseDegree` 齐次化 + 结点/极点/权重拼接**（不采样）。
- **前置缺口（本轮实测）**：端口**没有** `GeomBSplineCurve::increase_degree`、`knots()`、`multiplicities()`、`set_knots`，也没有 `GeomConvert_CurveToBSpline`（把 line/conic 转成 B-spline，`GeomConvert_CurveToBSpline.cxx`）。⇒ 忠实移植需先补这两件，故 A8 不能在本批收口（已把结论写进任务卡；`comp_curve_to_bspline` 当前**无生产调用者**，仅自身测试用）。

**批 10（T-63：`EvalD3` 与两处假 `d2`）已完成步 1 —— 2026-09-20**

- **落地**：`occt-core/src/elib/clib.rs`、`occt-geom/src/{circle,ellipse,hyperbola,parabola,bspline_curve}.rs`、`occt-geom2d/src/{circle,ellipse,hyperbola}.rs`。
  1. 忠实件：`clib::{circle_d3, ellipse_d3, hyperbola_d3}`（3D，逐行对 `ElCLib::CircleD3` `cxx:435-460`、`EllipseD3` `cxx:464-491`、`HyperbolaD3` `cxx:494-516`——注意双曲的 OCCT 实现把 `V3` 设成与 `V1` 相同的线性式，圆/椭圆则是 `-V1`）；`clib::{circle2d_d3, ellipse2d_d3, hyperbola2d_d3}`（2D，`cxx:809/843/878`，同样验证过 `V3` 的符号约定）。
  2. 3D 具体曲线补 `d3`：`GeomCircle`/`GeomEllipse`/`GeomHyperbola` → clib；`GeomBSplineCurve` → `eval_dn(u,3)`（`Geom_BSplineCurve::D3` = `BSplCLib::DN(...,3)`）；`GeomLine`/`GeomParabola` 保持零（`Geom_Parabola::EvalD3` `Geom_Parabola.cxx:195-200` 明确把 `V3` 置零）。
  3. 2D 具体曲线补 `d3`：`Geom2dCircle`/`Geom2dEllipse`/`Geom2dHyperbola`。
  4. **顺带修掉两处"假 `d2`"**（同一 A0 家族，也让 offset 的 `CalculateD1` 拿到正确二阶导）：3D `GeomHyperbola::d2`、`GeomParabola::d2` 原先返回 `GpVec::zero()`；现新增 `clib::{hyperbola_d2, parabola_d2}`（`ElCLib::HyperbolaD2` `cxx:378-401`、`ParabolaD2` `cxx:404-431`，含 `|Focal| <= gp::Resolution()` 时 `P=Loc, V1=YDir, V2=0` 的退化支路）并接线。
- **验证（全绿，与基线逐项一致）**：`occt-core --lib` **290/290**、`occt-geom --lib` **151/151**、`occt-geom2d --lib` **72/72**、`occt-topo --lib` **1293/1**（唯一红 = T-01）、四道 STEP 门禁 11/11、2/3、14/14、13/13；`cargo check` 四个 crate exit 0。
- **如实说明（避免过度声明）**：本批是**结构性忠实化**（把 OCCT 里存在的 `EvalD3`/`D2` 分支补齐）。A0 当时的临时探针（圆基 R=2 + 偏移 0.5 对解析圆）显示 `D0/D1/D2` 误差已为 0，说明该用例上 `d3` 项对结果不敏感；本批**未**新增数值验证（offset 类在仓内仍无生产调用者，属 A0 的潜在缺陷面），故不主张"数值已改善"，只主张"分派与公式已与 OCCT 一致"。
- **未完成（仍在 T-63）**：`Geom_OffsetCurve::EvalD2`（`cxx:311-330`）/`Geom2d_OffsetCurve::EvalD2`（`cxx:265-280`）的 `AdjustDerivative` 奇异支路——`isDirectionChange` 目前恒为 `false`（`offset.rs` 已就地标 UNPORTED）。

**批 11（A15/T-51 的可收口子集：`Surface::d2` 的有限差分默认值）已完成 —— 2026-09-20**

- **问题（审查 A15 第 5 条"失败/缺实现被兜底吞掉"）**：`Surface` trait 的 `d2` 默认实现用 `h = 1e-6` 的中心差分近似（`occt-geom/src/surface.rs:12-22`），而 OCCT 的 `Geom_Surface::D2` 是**纯虚**——每个具体类都解析计算（初等面 `ElSLib::*D2`、B 样条 `BSplSLib::D2`、偏移面 `Geom_OffsetSurfaceUtils`），OCCT 里根本没有有限差分兜底。实测 14 个 `impl Surface` 中只有 3 个（`bspline_surface`/`rectangular_trimmed`/`surface_of_linear_extrusion`）定义了自己的 `d2`，**五个初等面（plane/cylinder/cone/sphere/torus）都在吃这个近似**，而 `d2` 在 live 网格路径上被 `meshing/edge_discret.rs:1218` 与 `meshing/geomlib_norm.rs:30` 使用。
- **落地**：
  1. `occt-core/src/elib/surface_eval.rs`：把 D2 拆成 per-type 忠实件 `plane_d2` / `cylinder_d2` / **新增 `cone_d2`**（`ElSLib::ConeD2` `ElSLib.cxx:867-…`：`Vuu = -R·Vxy`、`Vvv = 0`、`Vuv = sinA·DVxy`）/ `sphere_d2`（`cxx:975-1037`）/ `torus_d2`（`cxx:1039-1100`），`surface_d2(&SurfaceRef, …)` 改为委托（`SurfaceRef` 无 Cone 变体，锥面直接调 `cone_d2`）。
  2. `occt-geom/src/{plane,cylinder,cone,sphere,torus}.rs`：各自实现 `d2` 调用上述忠实件（`Geom_Plane::D2` → `ElSLib::PlaneD2` 等），不再走有限差分。
  3. `occt-geom/src/surface.rs`：默认实现就地标 `UNPORTED`，并把"仍继承该近似的 6 个实现者"（`GeomBezierSurface`、`GeomOffsetSurface`、`GeomSurfaceOfRevolution`、`GeomSurfaceOfLinearExtrusion`、`surface_fit`、`surface_to_grid`）逐一列名，避免被当成"已对齐"。
- **验证（全绿，与基线逐项一致）**：`occt-core --lib` **290/290**、`occt-geom --lib` **151/151**、`occt-geom2d --lib` **72/72**、`occt-topo --lib` **1293/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3；四个 crate `cargo check` exit 0。⇒ 在现有门禁模型上，用解析 D2 替换 1e-6 差分**未引起可见变化**（差分误差 1e-6 量级，低于网格偏转与 bbox 判据的敏感度）。
- **同批追加（A15 的自创阈值两处 + 一处出处订正）**：
  4. `occt-geom/src/curve_reparam.rs::gauss_solve` 的奇异判据 `best < 1e-14` → **`1.0e-20`**（OCCT `math_Gauss.hxx:45-49` 的 `MinPivot` 默认值，文档原文"If the largest pivot found is less than MinPivot the matrix A is considered singular"）。
  5. `occt-math/src/matrix.rs::MathMatrix::solve`（自述 `Source: math_Gauss.cxx`）的 `pivot.abs() < 1e-30` → 同样改 **`1.0e-20`**（同上依据）。
  6. `occt-geom/src/curve_reparam.rs::basis_values`（The NURBS Book A2.2 教科书算法）加**非 OCCT 出处**注记：OCCT 等价件是 `BSplCLib::BasisFuns`/`bspl`；文件头的 `Geom_Curve`/`Geom_BSplineCurve`/`GCPnts_AbscissaPoint` 只是调用侧来源。
  - 验证追加：`occt-math --lib` **215/215**（1 ignored）+ 上述全套门禁仍等于基线。

**批 12（A25/T-61：删除"先删邻三角形"自创补丁，改为 OCCT 的失败语义）已完成 —— 2026-09-20**

- **问题**：`decompose_simple_polygon` 在加耳三角形前会把"已挂两个三角形的耳链接"上的邻居三角形**删掉**再 `AddElement`（原 `delaun/polygon_meshing.rs:343-375`，注释自述 "free the slot before `AddElement`"）。OCCT `BRepMesh_Delaun.cxx:2259-2274` 只做 `AddLink` × 2 + `addTriangle`，**不碰既有三角形**；链接满时 `BRepMesh_PairOfIndex::Append`（`BRepMesh_PairOfIndex.hxx:41`）抛 `Standard_OutOfRange`，被 `BRepMesh_BaseMeshAlgo.cxx:52-62` 空 catch 吞掉 ⇒ 该面**没有网格**。删除补丁是仓库明令禁止的"OCCT 里不存在的规则"。
- **落地**：
  1. `Delaun` 新增 `failed: bool`（`meshing/delaun/constants.rs` 结构体 + 三处构造点）。
  2. `add_triangle`（`meshing/delaun/triangulation.rs`）由"直接 `add_element`（链接满则由 `DelaunPairOfIndex::append` panic）"改为**先判满**：任一链接 `elements_connected_to(..).extent() >= 2` ⇒ 置 `failed = true` 并返回 `false`（= OCCT 的抛出条件，不再 panic、也不动既有三角形）；否则照常建三角形并返回 `true`。`add_triangle_by_info`（`meshing/delaun/frontier.rs`）同样返回 `bool`。
  3. `decompose_simple_polygon`（`meshing/delaun/polygon_meshing.rs`）：**删除整段删除补丁**；`add_triangle_by_info` 返回 `false` 时清空多边形并返回（等价于 OCCT 异常把该面留成无网格）。
  4. 失败即止的传播：`mesh_polygon` 入口与分解循环、`create_triangles_on_new_vertices` 的尾部（`process_constraints` 前）都检查 `self.failed`。
- **验证（全绿，与基线逐项一致）**：`occt-topo --lib` **1293/1**（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3；`cargo check` exit 0。
- **实测网格变化（重要，订正 T-68 第九轮的"0 次触发"结论）**：`step_obj_parity --nocapture` 对比——Cube 24/12、Cylinder 146/140、Shape-2 3105/4792、ATU01038 17745/22119 **逐位不变**；**T0M 由 46514/46945 变为 46516/46962（v +2 / f +17）**。⇒ 该自创删除在 T0M 上**确实触发过**（第九轮只对前若干个及每 200 个 `NO_EAR` 采样打印，漏掉了这些事件），删掉后网格略有增加、门禁仍全绿；这是"移除自创行为"的**正向证据**，不是回归。
- **保留的差异（如实标注）**：OCCT 抛出后是**整个面**没有网格；端口现在只中止当前多边形的分解（此前已加入的三角形保留）。要让整面清零需要把失败标志上抛到面级管线（`IMeshData_Failure` 语义），属 T-61 的收尾项，已写进任务卡。

**批 13（A14/T-50 收尾：整包非移植件的 `UNPORTED` 声明）已完成 —— 2026-09-20**

- **落地**（10 处模块头，逐处写"不是谁的移植 + 真实出处"）：
  1. `occt-core/src/geom/mod.rs`：整包声明——OCCT 对应为 `GCPnts_*`/`BRepMesh_*`（离散化与面网格）、`GeomAPI_PointsToBSpline`/`math_*`（拟合）、`IntTools`/`BOPTools_AlgoTools`（平面谓词）、`BRepAlgoAPI_*`+`BOPAlgo_Builder`（布尔）。
  2. `geom/csg.rs`：OCCT **没有**体素布尔；真实件是 `BRepAlgoAPI_*`+`BOPAlgo_Builder`（消费方 `solid_union.rs` 已自述为近似）。
  3. `geom/delaunay.rs`：**出处订正**——原头部写 `Source: OCCT Poly_Triangulation and math_Recipes`，但 Bowyer–Watson 非 OCCT（`Poly_Triangulation` 只是三角形容器，`math_Recipes` 是书目）；OCCT 的 2D 约束面网格是 `BRepMesh_Delaun`+`BRepMesh_BaseMeshAlgo`（已移植在 `meshing/delaun`）。
  4. `geom/triangulate.rs`：耳切非 OCCT（OCCT 走 `BRepMesh` + `Poly_MakeLoops`/`Poly_Connect`）。
  5. `elib/measure.rs`：自造 Simpson/弦长/中点矩形面积 ⇒ OCCT 等价件 `GCPnts_AbscissaPoint` / `GProp_*`。
  6. `hlr.rs`：非 `HLRBRep_*`（TKHlr），只有轻量投影/线框。
  7. `render_svg.rs`：OCCT 不产 SVG（其可视化为 TKV3d）。
  8. `viz_scene/mod.rs`：非 `AIS_*`/`V3d_View`/OpenGl。
  9. `draw/mod.rs`：非 TKDraw（`Draw_Interpretor`）。
  10. `xcaf/mod.rs`：非 `XCAFDoc_ShapeTool`/`TDocStd_Document`。
  （`validate.rs` 的"`is_valid()` 恒 true 不得作门禁"已在批 8 声明。）
- **验证（纯声明，无行为改动）**：`occt-core --lib` 290/290、`occt-geom` 151/151、`occt-geom2d` 72/72、`occt-topo --lib` 1293/1（唯一红 = T-01）、四道 STEP 门禁 11/11、2/3、14/14、13/13；四个 crate `cargo check` exit 0。
- **旁支**：`geom/` 内其余文件（`fit.rs`/`fit2.rs`/`polygon_*`/`polyline*.rs`/`curve_*`/`mesh_analysis.rs`）不再逐一加头，由 `geom/mod.rs` 的整包声明覆盖（已在该声明里给出 OCCT 对应件）。

**批 14（A15/T-51 续：`IntAna_IntLinTorus` 的根校验 + `Extrema_ExtElSS` 本地扩展声明）已完成 —— 2026-09-20**

- **落地**：
  1. `occt-geom/src/intana/line_torus.rs::line_torus_intersect`：根校验由自创的**隐式方程残差** `|(ρ−R)² + z² − r²| < 1e-7` + `1e-7` 去重，改为 OCCT `IntAna_IntLinTorus.cxx:98-106` 的做法——用 `ElSLib::Parameters(gp_Torus, P, u, v)`（端口 `slib::torus_parameters`）反求参数、`ElSLib::Value(u, v, T)` 回代（`slib::torus_value`），比较**平方距离**与 **`1.0e-10`**（OCCT 的 `a0 > 0.0000000001`）；**去掉去重**（OCCT 把每个有效根原样存入）。
  2. `occt-geom/src/extrema_ss.rs::{sphere_sphere_extrema, plane_sphere_extrema}`：加 `UNPORTED / 本地扩展` 注记——OCCT `Extrema_ExtElSS::Perform(gp_Pln, gp_Sphere)`（`Extrema_ExtElSS.cxx:62-69`）与 `(gp_Sphere, gp_Sphere)`（`:77-83`）在设好 `myDone/myNbExt` 后 **`throw Standard_NotImplemented`**，OCCT 根本没有闭式解 ⇒ 端口这两个函数是**本地扩展**，不能与 OCCT 行为对比（模块头原本已提到，这次补到函数级）。
- **验证**：`occt-core --lib` 290/290、`occt-geom --lib` **151/151**（其中 `intcurvesurface/tests.rs` 的线–环交用例 `x = ±2, ±4` 仍通过 ⇒ 新根校验与原实现给出同样的根）、`occt-topo --lib` 1293/1（唯一红 = T-01）、四道 STEP 门禁 11/11、2/3、14/14、13/13；`cargo check` 全绿。

**批 15（A10 收尾子集 + 新发现 A31：删除伪造的回转面类）已完成 —— 2026-09-20**

- **A31（新发现，已修）**：`occt-geom/src/revolved.rs::GeomRevolvedSurface` 是一份**伪造的**回转面实现——`d1` 的 `dv` 直接取 basis 点坐标 `GpVec::new(pt.x(), pt.y(), pt.z())`（注释自承 "simplified derivative"，既不是导数也不在轴坐标系），`transform` 是空实现；而忠实件 `surface_of_revolution.rs::GeomSurfaceOfRevolution`（含解析 `d2`）早已存在。**实测 `git grep GeomRevolvedSurface` 只有定义与 `lib.rs` 的 re-export ⇒ 零消费者**（`revolution_basis_curve` 的两处消费方 `approx_same_parameter.rs:460`、`cos_locate.rs:463` 走的是 trait 方法，忠实件已实现）。⇒ 按"自创代码摘除"删除该文件与 `lib.rs` 的两行接线；四 crate `--lib` 与门禁全部等于基线，证明它确是死代码。
- **A10 收尾（`Surface::d2` 的最后一个可解析实现者）**：`occt-geom/src/bezier_surface.rs` 补 `d2`——通过 `osculating_bspline()` 的**精确 B 样条等价表示**（Bezier = `[0,1]` 上满重结点的 B 样条，各阶导数相同）委托求值，取代 trait 的 `h = 1e-6` 中心差分；仅退化面（`n_u < 2 || n_v < 2`）走解析 `d1` + 零二阶导。
- **A10 现状**：`Surface::d2` 的解析实现已覆盖 **plane / cylinder / cone / sphere / torus / Bezier / BSpline / RectangularTrimmed / SurfaceOfRevolution / SurfaceOfLinearExtrusion / OffsetSurface**（偏移面的 D2 走 `offset_surface_utils::compute_derivatives` 路径，trait 默认只余 `surface_fit`/`surface_to_grid` 两个**非 OCCT** 的网格工具类，已在 trait 文档中列名标 `UNPORTED`）。
- **验证（全绿，与基线逐项一致）**：`occt-core --lib` 290/290、`occt-geom` **151/151**、`occt-geom2d` 72/72、`occt-topo --lib` 1293/1（唯一红 = T-01）、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 2/3；`cargo check` 全绿。

## 14. A0–A31 整改归档（查号表：A# → 任务 ID → 当前状态）

> 本表只用于 **A 编号 ↔ 任务 ID** 查号。**状态一律以 §3.1–3.5 为准**；本节旧的状态矩阵、旧"基线"行、以及 §14.1–14.4 的四份中途交接已作废（见 §14.1–14.4 说明）。逐项证据在 `specs/_audit/_index.md` 与 §3.6 批次明细。

| A# | 任务 ID | 一句话 | 当前状态（指向 §3） |
|---|---|---|---|
| A0 | T-35 / T-63 / T-75 | offset 曲线 D0–D3、基曲线 `EvalD3`、`AdjustDerivative`、初等曲线 `EvalDN` 至 5 阶 | ✅ 归档（§3.5）；**遗留 UNPORTED**：非初等、非 B 样条基曲线的 `EvalDN` 4/5 阶仍为 trait 默认零值 |
| A1 | T-37 / T-67 | `Extrema_ExtPS/GenExtPS/ExtPElS` 移植与 39 处调用点迁移 | ◐ 分步 1+2 done；**步 3 = §3.2 序 8** |
| A2 | T-38 | `is_inside` 改 `SolidClassifier` + `SolidExplorer::myMapEV` | ✅ 归档 |
| A3 | T-39 / T-71 | 写侧 `IsKind` 分派；未识别曲线写 `$`、`Geom_BoundedCurve` 臂 | ✅ 归档 |
| A4 | T-40 | F 分支改 `BRepClass_FaceClassifier`（删距离门） | ✅ 归档 |
| A5 | T-41 / T-79 / T-80 | 网格/体素布尔摘除；曲面 operand 的 `BOPAlgo_BOP` | ◐ dispatch 级已摘（第 59 轮）；**余项 = §3.2 序 7 + §3.3 T-41** |
| A6 | T-42 | `gcpnts` 二分族改忠实 `gcpnts_perform` | ✅ 归档 |
| A7 | T-43 / T-65 / T-66 | 点–曲线、曲线–曲线极值忠实化（`GGenExtCC::Perform`） | ✅ 归档（就地 UNPORTED：局部引擎仅 BFGS、`CellFilter` 线性表、长度比用 `gcpnts::curve_length`） |
| A8 | T-44 | `GeomConvert_CompCurveToBSplineCurve` | pending（触发式）→ §3.3 T-44 |
| A9 | T-45 | 假出处改正（7 处） | ✅ 归档 |
| A10 | T-46 / T-73 / T-76 | `surface_d2`/`bsphere`/`choose_left_way`/`GpDir2d::Angle` | ✅ 归档 |
| A11 | T-47 / T-77 | 三项自创删除 + `FindSolutions` bbox 递归 | ✅ 归档（T-77 = R2-1，批 80） |
| A12 | T-48 | `brepfeat` 解析体积覆盖、`features.rs:501` 分辨率 | **blocked（T-80）** → §3.3 |
| A13 | T-49 / T-54 | 未裁剪 UV 窗口栅格 → `BRepMesh_FaceDiscret`；平面耳切 → 约束 Delaunay | **blocked（T-68/T-69）** → §3.2 序 2/3、§3.3 T-54 |
| A14 | T-50 | 模块头 `UNPORTED` + 真实出处 | ✅ 归档 |
| A15 | T-51 | 失败通道/积分忠实化；**反解 UNPORTED** | ◐ → §3.3 T-51 余项 |
| A16 | T-52 / R2-18 | 投影 + 求交采样器删除 | ◐（非解析 2D 组合）→ §3.3 R2-18 余项 |
| A17 | T-53 | `brep_exchange` 静默回退删除 | ✅ 归档 |
| A18 | T-54 | 平面耳切/桥洞 → 约束 Delaunay | pending → §3.3 |
| A19 | T-55 | `WIREFRAME_FALLBACK_RATIO_MAX` 摘除 | 阈值已摘；**逐面回退 blocked（T-69）** → §3.3 |
| A20 | T-56 | 边域两次 `Project`（含 A29 根因） | ✅ 归档 |
| A21 | T-57 | `model_healer` 退化支路 + 平方距离判据 | ✅ 归档 |
| A22 | T-58 | `ProjectAct` 解析臂；删三点外心族 | ✅ 归档 |
| A23 | T-59 | `face_uv_bounds` → `BRepTools::AddUVBounds` | **blocked（T-49/T-69）** → §3.3 |
| A24 | T-60 | `classify_surface` 精确类判定 | ✅ 归档 |
| A25 | T-61 | 删除补丁 + 失败标志上抛到面级 | ✅ 归档 |
| A26 | T-62 / T-78 / T-84 / T-85 | IGES/PLY/STL/VRML 写侧与 STEP 写侧 | ◐：STEP 采样重拟已删、IGES 写侧 1–11 步 done；**余项 = §3.3（R2-19、T-85 余项、T-62 余项）** |
| A27 | T-64 | 15 处 `continuity()` 改 `CN(6)` | ✅ 归档 |
| A28 | （随 T-68） | `FaceGauss` 边界弧 pcurve 优先 | ✅ 归档 |
| A29 | T-70 | `ExtPElC` 直线臂参数系平移 | ✅ 归档 |
| A30 | T-72 | `brep_surface::classify_surface` 精确类型 | ✅ 归档 |
| A31 | T-74 | 伪造 `GeomRevolvedSurface` 删除 | ✅ 归档 |


**批 16（A26/T-62 续：STL 头/嗅探对齐 + A12 与 OBJ `vn` 的实地结论）—— 2026-09-20**

- **落地（`occt-core/src/io/stl.rs`，逐行对 `RWStl.cxx`）**：
  1. 二进制写头：自定义的 `"solid generated by occt-core"` → OCCT 的 80 字节字面量 **`"STL Exported by Open CASCADE Technology [dev.opencascade.org]"`**（零填充；`RWStl.cxx:374-375`）。
  2. ASCII 写头：`"solid\n"` → **`"solid \n"`**（OCCT 原文注释：`// note that space after 'solid' is necessary for many systems`，`RWStl.cxx:304-305`；且 OCCT 不写 solid 名称）。
  3. 格式嗅探：`starts_with("solid") && contains("facet")` → **只看前 5 字节 `solid`**（`RWStl.cxx:523-531`：`isAscii = (strncmp(aHeader, "solid", 5) == 0)`）——`contains("facet")` 是自创的额外条件。
- **验证**：`occt-core --lib` **290/290**（STL 读写与 ASCII/binary 往返测试全过）、`occt-topo --lib` 1293/1（唯一红 = T-01）、四道 STEP 门禁 11/11、2/3、14/14、13/13。
- **T-48/A12 的实地结论（本轮未动手，已写入任务卡）**：`brepfeat/features.rs:107-109` 的 `result.volume = (v0 + πr²h − overlap).max(0)` 是**解析体积覆盖**，使 `brepfeat/tests.rs:172` 的 `after.volume > before` 恒真；但删掉它会让该断言依赖真实网格体积，而代码注释已自述"the boolean does not close through-hole topology"（体积不可靠）⇒ **删覆盖会立刻把 `boss_thru_all_pierces` 打红（1293→1292/2，劣于基线）**。⇒ A12 的体积部分**被 BOP 缺口阻塞**（与 A5/A11 同源），本轮只确认了因果关系，未改代码。
- **OBJ `vn` 的实地结论（订正审查措辞）**：写侧**并非**"恒空 `vn`"——`occt-core/src/io/obj.rs:130` 在 `mesh.normals` 非空时会写 `vn`；OCCT 也是 `if (theFace.HasNormals())` 才写（`RWObj_CafWriter.cxx:245,281,317`）。真正的缺口是**网格侧不产出法向**（`brep_to_obj` 的 `ObjMesh.normals` 恒空），属网格管线（A13/A18）而非 OBJ 写侧。⇒ A26 该子项应改写为"网格法向未产出"，不再是写侧偏离。

**批 17（A15/T-51 的失败通道：曲线族改为可失败语义）—— 2026-09-20**

- **背景（消费者实测，决定改动边界）**：`point_curve_extrema`/`curve_curve_extrema`（返回 `ExtremaPair` 的**带兜底**版本）在 `occt-topo` **无生产调用者**——生产路径用的是 `point_curve_extrema_all`（`brep_extrema.rs:18`、`edge_edge/mod.rs:29`）与 `extrema_surf::{curve_surface_extrema_all, point_surface_extrema}`；曲面族**有**生产消费方（`bean_face_exact.rs:278` 调 `extrema::curve_surface_extrema`），其忠实化要等 T-67 的 `Extrema_ExtPS/GenExtPS` ⇒ 本轮只动曲线族、曲面族仅就地声明。
- **落地**：
  1. `occt-geom/src/extrema_pc/point_curve.rs::{point_curve_extrema, point_curve_max_extrema}` → 返回 **`Option<ExtremaPair>`**，删除"退化时黄金分割兜底"（`refine_curve_point` 那条自创路径）；`None` 即 OCCT 的 `IsDone() == false`（`Extrema_GGExtPC.hxx:531`、`:545-550`）。
  2. `occt-geom/src/extrema_cc/curve_curve.rs::curve_curve_extrema` → **`Option<ExtremaPair>`**，删除 16×16 参数网格兜底与 `pair_cc(zero, 0, zero, 0)` 伪造对；`None` 即 `Extrema_GGenExtCC.hxx:691-695` 的 `if (aNbSol == 0) { myDone = false; return; }`。
  3. `occt-geom/src/extrema.rs` 的两个公开包装同步改为 `Option`。
  4. **曲面族兜底就地标 `UNPORTED`**：`extrema_ss.rs::fallback_ss`、`extrema_surf/numeric_extrema.rs::{fallback_point_surface, fallback_curve_surface}`——注明 OCCT 无此兜底（`Extrema_ExtElSS.cxx:62-83` 抛 `Standard_NotImplemented`、`GGExtPC.hxx:531` 置 `done=false`），因 `bean_face_exact.rs:278` 依赖而暂留，移除需先做 T-67。
- **测试适配（仅机械加 `.expect(...)`，不放宽断言）**：`extrema.rs` 3 处、`extrema_pc/tests.rs` 3 处、`extrema_cc/tests.rs` 1 处、`occt-topo/tests/phase6_integration.rs` 1 处。
- **验证（全绿，与基线逐项一致）**：`occt-core --lib` 290/290、`occt-geom` **151/151**、`occt-geom2d` 72/72、`occt-topo --lib` 1293/1（唯一红 = T-01）、`phase6_integration` **5/5**、四道 STEP 门禁 11/11、2/3、14/14、13/13；`cargo check --tests` 两 crate exit 0。

**批 18（A26/T-62：PLY 去掉全局焊接；新登记一处预先存在的红）—— 2026-09-20**

- **落地（`occt-topo/src/brep_exchange.rs::brep_to_ply`）**：删除 `weld_vertices(&mut mesh, 1e-9)`。OCCT 的 PLY 写侧 `RWPly_PlyWriterContext` 按**每个 CAD 面自己的 node 块**写出，用递增的 `myVertOffset` 给面内索引加偏移（`RWPly_PlyWriterContext.cxx:276`（三角）、`:300`（四边）），**从不跨面合并节点** ⇒ 盒体应是 24 条 vertex 记录（6 面 × 4），与 OBJ 路径（`brep_to_obj`，注释已述 "24 nodes / 12 triangles"）一致；改前焊接成 8 节点共享网格，OCCT 不会这么写（audit A26）。
- **未动（有依据）**：`brep_to_stl_binary` 的焊接保留——STL 每个 facet 都自带 3 个顶点，索引是否共享**不影响输出字节**；OBJ 写侧是共享顶点表 + 面索引（`RWObj_ObjWriterContext` 每节点写一次 `v`），保持现状正确。
- **新登记（预先存在的红，非本批引入）**：`occt-topo/tests/phase3_integration.rs::primitives_measure_correctly` 在 `:79` 的 `is_inside(&sphere.solid.0, &Gpnt::new(0.5,0.5,0.5))` 断言失败——这是**球面点分类**（A2/T-38 的 7×7 采样 + 射线奇偶）的已知缺口，与本批（仅改 PLY 导出）无关；该集成测试不在既有基线清单内，故计入"旁支红"而不影响门禁。⇒ 建议并入 A2/T-38 修复范围。
- **验证（门禁与 lib 逐项等于基线）**：`occt-core --lib` 290/290、`occt-topo --lib` 1293/1（唯一红 = T-01）、`phase3_integration` 3/4（红 = 上述预先存在项）、四道 STEP 门禁 11/11、2/3、14/14、13/13；`cargo check --tests` 两 crate exit 0。

**批 19（T-63 步 2：offset 的 `AdjustDerivative` 奇异支路 + 初等曲线 `EvalDN`）—— 2026-09-20 第 41 轮**

- **背景**：T-63 步 1 只补齐了基曲线的 `EvalD3`；`EvalD2` 的奇异支路仍把 `isDirectionChange` 硬编码为 `false`（两处 `offset.rs` 已标 UNPORTED）。本批按 OCCT 原文逐行移植该支路。
- **落地（3D，`Geom_OffsetCurveUtils.pxx:322-390` / `Geom_OffsetCurve.cxx:311-330`）**：
  1. 新增自由函数 `occt-geom/src/offset.rs::adjust_derivative`：`aTol = gp::Resolution()`（`pxx:331`）、`aMinStep = 1e-7`、`aMaxDerivOrder = 3`、`DivisionFactor = 1e-3`；无界参数域（`u_sup >= RealLast()` 或 `u_inf <= RealFirst()`，端口用 `f64::MAX`/`f64::MIN`）取 `du = 0` ⇒ `aDelta = 1e-7`；`do { V = EvalDN(u, ++anIndex) } while (|V|² <= aTol && anIndex < 3)`；`u = (u − u_inf < aDelta) ? u + aDelta : u − aDelta`；`V1 = P(min) → P(max)`；`isDirectionChange = V·V1 < 0`；`D1 = V·sign`，`for i in 1..maxDerivative { D_{i+1} = EvalDN(u, anIndex+i)·sign }`；恒返回 `true`（OCCT 唯一 `false` 路径是 `EvalDN` 抛异常）。
  2. `GeomOffsetCurve::d2` 接线：当**基曲线** `D1` 的平方模 `<= gp::Resolution()` 时以 `maxDerivative = 3` 调用（`aDummyD4` 丢弃），把 `is_direction_change` 传给 `calculate_d2`（`pxx:301-305` 的 `theD3.Reverse()` 语义不变，因为 `d3` 由本函数覆写）。
- **落地（2D，`Geom2d_OffsetCurveUtils.pxx:296-364` / `Geom2d_OffsetCurve.cxx:258-280`）**：`occt-geom2d/src/offset.rs::adjust_derivative` 为同一算法的逐行转写（`pxx:305-313` 的常数、`:326-333` 的 do-while、`:337-344` 的 `u` 偏移、`:351-360` 的符号与 `D2..D4`），`Geom2dOffsetCurve::d2` 同样接线。
- **依赖补齐（`AdjustDerivative` 要 `EvalDN` 到 5 阶；此前 3 阶以上一律返回零）**：
  1. 3D 新增忠实件 `occt-core/src/elib/clib.rs::{line_dn, circle_dn, ellipse_dn, hyperbola_dn, parabola_dn}`（`ElCLib::LineDN` `ElCLib.cxx:911-918`、`CircleDN` `:922-953`、`EllipseDN` `:957-992`、`HyperbolaDN` `:996-1016`、`ParabolaDN` `:1020-1045`，含圆/椭圆按 `N mod 4` 的四类残数循环与抛物 `|Focal| <= gp::Resolution()` 退化臂）；5 个 3D 初等曲线（line/circle/ellipse/hyperbola/parabola）新增 `eval_dn` 覆写（`Geom_Line.cxx:207-219`、`Geom_Circle.cxx:184-191`、`Geom_Ellipse.cxx:230-237`、`Geom_Hyperbola.cxx:269-276`、`Geom_Parabola.cxx:205-212`）。
  2. 2D `Curve2d` 新增 `eval_dn`（`Geom2d_Curve::EvalDN`，默认 1..3 阶走 `d1/d2/d3`、`N<1`/`N>3` 返回零，与 3D trait 同一约定），5 个 2D 初等曲线覆写并接已有的 `clib2d::{line,circle,ellipse,hyperbola,parabola}_dn_*`（`ElCLib.cxx:1049-1188`）。
  3. 两侧 `GeomTrimmedCurve` 补 `d3`/`eval_dn` 转发基曲线（`Geom_TrimmedCurve.cxx:231-243`、`Geom2d_TrimmedCurve.cxx:273-283`）；3D 那个视图在端口里重参数化到 `[0,1]`，故按既有 `d1/d2` 的约定乘 `s^N`（`s = last − first`）。
- **验证（门禁与基线逐项一致）**：`occt-core --lib` **290/290**、`occt-geom --lib` **151/151**、`occt-geom2d --lib` **72/72**、`occt-math --lib` **215/215**、`occt-topo --lib` **1293/1**（唯一红 = T-01 `brepfeat::tests::groove_cuts_cylinder`）、`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **2/3**（红 = T-05，报错文本仍逐字为 `Offset divergence 2208.0 vs shape_volume 1612.9`，与本批前记录一致）、`phase6_integration` **5/5**、`phase3_integration` 3/4（红 = 既有 A2 `primitives_measure_correctly`）；`cargo check` 四 crate exit 0。
- **支路的实证（临时探针，已删除）**：用一条在 `u=0` 处一阶导为零的 B 样条基曲线（`P0 = P1`、degree 2、夹持结点）驱动 `GeomOffsetCurve::d2(0)`，插桩输出 `[dbg-adjust] u=0 index=2 sign=1 v=(2,2,0)`——即支路被命中一次（`V = D2`、弦方向判定为"无换向"），且偏移点/切线为解析值（法向 `D1 × dir` 归一后偏移 0.5 ⇒ `p = (0.353553, −0.353553, 0)`、`D1 = (2,2,0)`）；`u=0.5`（基曲线 `D1 ≠ 0`）不命中该支路。探针文件与 `eprintln` 已删除，`git grep` 复核 0 残留。
- **未完成 ⇒ T-75（新登记，零消费者）**：offset 类自身的 `EvalD3`（`EvaluateD3`+`CalculateD3`）与 `EvalDN` 未移植 ⇒ `GeomOffsetCurve`/`Geom2dOffsetCurve` 的 `d3` 仍是 trait 默认零值；两处 offset 文件头已把 UNPORTED 清单刷新为当前真实缺口（去掉"无曲线覆写 `d3`"这条已过期项，新增该两条并注明 `Bezier` 的 `EvalDN` 仍未移植）。
- **旁支观察（未动手）**：`occt-core/src/elib/clib.rs::parabola_d2` 的退化判据用 `precision::RESOLUTION`（`1e-12`），而 OCCT `ElCLib::ParabolaD2` 用的是 `gp::Resolution() == RealSmall() == DBL_MIN`（本批新增的 `parabola_dn` 用的是后者）；2D `Geom2dParabola::d1/d2` 亦未带 `parabola_d1_ax22d`/`parabola_d2_ax22d` 里的退化臂。两处都只在 `|Focal| <= 1e-12` 的退化抛物线上可见，属旁支，仅登记。

### 14.1–14.4 历史交接（已作废，2026-09-21 删除）

> 原 §14.1（第 40 轮复核）、§14.2（第 63 轮交接）、§14.3（第 82 轮交接）、§14.4（第 89 轮交接）是会话中期的交接文档，其状态与"基线数字"已被 **§2（第 111 轮 A/B 实测）** 与 **§3.1–3.5** 取代，故整体删除。批次细节看 §3.6 的批 1–99 明细与 §7 的轮次日志（含第 111/112 轮）。

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

### 6.1 门禁波次的第一步：R2-18-gate 差分（批 85–90 前/后）

收到"导出 obj 测试"指令后，**先做本差分再看任何红项**（批 85/86 改 2D 投影与求交、批 87 改圆锥包围盒、批 88 改 STEP 边域投影回退、批 89 改 healing 的 pcurve 结点吸附、批 90 改周期 B 样条包围盒；四者都在 P0/P1 门禁的下游，见 R2-18-gate 卡）。

```powershell
cd D:\source\repos\dogs
$M = 'crates/occt-topo/Cargo.toml'
# ① 现工作区（batch 89 之后）
cargo test --manifest-path $M --no-fail-fast 2>&1 | Select-String 'test result:|Running tests|panicked at' | Tee-Object $env:TEMP\gate_after.txt
cargo test --manifest-path crates/occt-geom2d/Cargo.toml --lib 2>&1 | Select-String 'test result:' | Tee-Object -Append $env:TEMP\gate_after.txt
# ② 批 84 结束点（= 批 85 之前）作为对照基线
git worktree add --detach .target-pre85 9d8e596
cd .target-pre85
cargo test --manifest-path $M --no-fail-fast 2>&1 | Select-String 'test result:|Running tests|panicked at' | Tee-Object $env:TEMP\gate_before.txt
cargo test --manifest-path crates/occt-geom2d/Cargo.toml --lib 2>&1 | Select-String 'test result:' | Tee-Object -Append $env:TEMP\gate_before.txt
cd ..; git worktree remove --force .target-pre85
# ③ 逐项比对两文件；差异项按 §3 R2-18-gate 的"第一嫌疑"顺序归因（85→86→87→88→89）
```

> 纪律复核：差分前不要动任何源码；`git worktree` 路径命中 `.gitignore` 的 `/.target-*/`；**禁止**在共享树 `git stash`。

### 6.2 门禁波次 playbook（收到指令后按序执行；勿跳步）

> **执行状态（2026-09-21）**：用户指令"导出 obj / 跑门禁" ⇒ **步骤 1–5 已执行完毕**（步骤 2 差分 = 零回归；步骤 3 全门禁见 §2；步骤 5 报告见 §7 第 111 轮）。**步骤 6 已执行 = 批 99（R2-17 结项，门禁逐项相同，见 §7 批 99）**。**步骤 7（R2-7…R2-10 分诊）见 §7 第 112 轮报告。**

| 步 | 动作 | 命令/判据 | 产出 |
|---|---|---|---|
| 1 | 工作区确认 | `git status --porcelain` 只应有未跟踪的 `data/iges/`；`git log -1` = 最新批 | — |
| 2 | **R2-18-gate 差分**（§6.1） | 现工作区 vs `9d8e596`（批 84 末）逐项比对 `test result:` | `gate_before/after.txt` + 差异清单 |
| 3 | 全门禁 + IGES 校验 + 导出门禁 | §6 的 5 条 + `cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example iges_check -- <18 模型>`（批 84/4 改过 IGES 写侧） | 新 §2 快照（R2-15 收口） |
| 4 | 差异/回归归因 | 若某用例由绿转红：按 **85 `e6c12c2` → 86 `2eccb5b` → 87 `ae67e92` → 88 `0ff1ac2` → 89 `d017d9f` → 90** 顺序用 `git worktree add --detach .target-bN <commit>` 二分，定位首个 offender | 首个 offender + 差异证据 |
| 5 | 停止并报告（门禁 3/4） | 失配即停：写清"哪条门禁、哪条用例、期望 vs 实际、首个 offender"，**不**顺手修、不新写测试 | §7 报告 |
| 6 | 无回归时：R2-17 单独一批 | 删 `perform` 的 `(Circle, Circle)` 分派与 `compute_circle_circle_full`（OCCT 无此分支）→ 复跑第 3 步 | R2-17 结项或回退 |
| 7 | 之后才动 R2-7…R2-10 | 按第 3 步的红项顺序（T-01/T-03/T-04/T-05）分诊 | 更新 §3/§14 |

**预期影响表（第 2/3 步解读用，避免把"预期变化"误当回归）**：

| 批次 | 预期会动的门禁 | 理由 |
|---|---|---|
| 85 `e6c12c2` | 2D 相关用例（`geom2d_api`/`pave_de` 消费） | 投影改 `Extrema_ExtPC2d` |
| 86 `2eccb5b` | 同上 + 重合元素语义变化（0 个孤立点） | 求交改 `IntAna2d` |
| 87 `ae67e92` | `occt-topo --lib` 中依赖 bbox 的用例、`bop_builder2_boss`/`phase19`、四道 STEP 门禁（间接） | 圆锥盒由采样变解析（更紧）⇒ 盒递归剪枝变化 |
| 88 `0ff1ac2` | 四道 STEP 门禁（边域投影回退臂） | `ProjectAct` 的 `!OK` 五臂 |
| 89 `d017d9f` | `occt-topo --lib` 的 healing/notch 用例、网格门禁（T0M 等） | `CorrectParameter` 结点吸附经 `FixNotchedEdges` |
| 90（§7 批 90） | **订正（批 91）：实际无影响** ——`box_bspline` 的周期臂判 `curve.is_periodic()`，而端口读入的 B 样条恒 `periodic=false`（见 R2-21）；待该标志能进来时才进本表 | 周期 B 样条盒由 33 点采样改为 OCCT 的 knot-span 采样 + 极点归约（实现忠实，当前不可达） |
| 91（§7 批 91） | **无**（只新增 `occt-core/src/convert/{conic_to_bspline,conic_curves}.rs` 与三个 `gp` 访问器，仓内零调用点） | 不加进差分；批 92 把 `GeomConvert::CurveToBSplineCurve` 接到写侧后再并入本表 |
| 92（§7 批 92） | **读侧无**；仅 STEP 写侧 Bezier 产物改走 `GeomConvert::CurveToBSplineCurve`（逐字段同结果，写侧不进 OBJ 门禁） | 新增入口 + `gp_Trsf::SetTransformation(From,To)` + `knots::nb_poles`；本批不动读路径，故不进差分 |
| 93（§7 批 93） | **无**（周期表示/求值只可能经 `CurveToBSplineCurve` 进入，而仓内唯一调用点只在 Bezier 分支 ⇒ 圆锥周期臂不可达；读侧未动） | 新增 `knot_sequence_periodic`/`set_periodic`/周期 `d0,d1,d2`；批 94 读侧接线后本行需改为"读侧周期曲线求值与包围盒"并进差分 |
| 94（§7 批 94） | **读侧首次变为可达的周期曲线**：`occt-topo --lib` 中依赖 bbox/求值的用例、`bop_builder2_boss`/`phase19`、四道 STEP 门禁、`export_data_obj`（语料中 linkrods/screw/Shape/ATU01038/bottom/motoc/top 含 `closed_curve=.T.` 的 B 样条） | ① 闭曲线（`closed=.T.` ∧ degree>1 ∧ `IsClosed()`）按 `StepToGeom.cxx:920-926` 强制周期化 ⇒ 极点数 −1、改走周期求值；② 描述子"像周期"时按周期表示构造（`:873-894`）；③ 重数 > degree+1 由"报错"改为"钳制+裁剪首尾极点"；④ 重复结点按 `Epsilon` 归并；⇒ **批 90 的 `box_bspline` 周期 arm 自此可达** |
| 95（§7 批 95） | **无（读侧不动）**：仅 IGES **写**侧整周椭圆由 104 改 126 ⇒ 只影响 `iges_check` 的实体统计，OBJ 门禁不读 IGES | 新增 `reparametrize`/`set_knots`/`set_not_periodic`/`distinct_knots_and_mults`/`GpAx2::Rotated` + `emit_whole_period_ellipse` |
| 96（§7 批 96） | **读侧（2D/pcurve）**：`occt-topo --lib` 中 pcurve 相关用例、`phase*`、四道 STEP 门禁、`export_data_obj`/网格门禁（pcurve 变周期会改 pcurve 求值与后续网格化） | ① 2D 臂接入 `MakeBSplineCurveCommon`（重数归并/钳制/极点裁剪/周期判定/闭曲线强制周期化）；② 周期 pcurve 走 `curve_dn` 周期求值（③④ 同 §6.2 第 94 行的重数钳制/重复结点归并语义） |
| 97（§7 批 97） | **无（写侧）**：IGES 周期 B 样条边由"弦线"改 126 ⇒ 仅 `iges_check` 的实体统计变化 | `emit_bspline_curve` 前置 `SetNotPeriodic` 副本（`GeomToIGES_GeomCurve.cxx:294-307`） |
| 98（§7 批 98） | **写侧为主**：① IGES——窄区间 B 样条边由"弦线"改 126 子曲线（`iges_check` 实体统计），且所有 126 的 `UMin/UMax` 由请求值改为钳制值（全区间边下二者相同 ⇒ 多数无变化）；② STEP 写侧——trimmed 的 B 样条基曲线改由 `Segment` 产出（此前写"重映射结点"替身 ⇒ `B_SPLINE_CURVE_WITH_KNOTS` 的极点/结点数会变）；③ **读侧只有一个变化**：无（`InsertKnots/Segment/SetOrigin` 只被上面两处写侧调用，仓内无读侧调用点）⇒ 四道 STEP/OBJ 门禁与 `export_data_obj` 预期不变 | `Geom_BSplineCurve::{InsertKnots,Segment,SetOrigin}` + `BSplCLib::{PrepareInsertKnots,InsertKnots,BoorScheme}` |

> **注意**：批 80–84 在本差分对照点之前，故其影响不在第 2 步里；若要分开验证它们（R2-1 的盒递归、R2-2 的求交接线、R2-4 的 IGES 重编号、R2-5 的 SRR 组合、批 84 的 128 反周期化），用同一二分法把对照点前移到 `f4d4d0b^`（= 批 79 末）即可，共 5 个对照点。

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
- **探针 2 轮**（均已撤除；`bop_bop.rs`、`bop_builder2/arguments.rs`、`brepfeat/tests.rs` 三文件复核为未修改状态）：
  1. `build_shape`/`build_rc`：`op=Cut dim0=3 dim1=3 open=false gf_solids=2 gf_faces=256`；`build_rc obj_src=1 tool_src=1 it_shapes=1 check_keys=1 it_exp=1 tool_sets=1`；
  2. `fill_images_solids`：`fill_in3d in_parts=1`，`build_split_solids_occt` 后 **两个 argument 都有镜像，但每个只有 1 个**。
- **结论（链路已闭合）**：T-01 的表现是"Cut 只切 0.414/1.6、Common/Fuse 空"，根因是**每个源 solid 只产出一个镜像（=未被切割）**，即 `BuildSplitSolids` 拿不到柱侧片被工具面切开的碎片 ⇒ 上溯到**面级切割**（`FillImagesFaces` / PaveFiller 的 planar-planar FF 干涉）。
- **门禁复核**：topo lib 1,293/1,294（仍仅 T-01），探针零残留。
- **下一轮入口（明确）**：`bop_build_faces::fill_images_faces` → 统计柱侧片 `history.image(face)` 条数；查 PaveFiller 是否给出该面对的 FF 干涉；对 OCCT `BOPAlgo_Builder_2.cxx::BuildSplitFaces`。此链同时是 T-03（boss 合并）、T-04（phase19 DS 计数）的共同上游，建议按"BOP 簇"一次收敛。

### 2026-09-20 03:00–04:00 · DSH 会话（goal round 2：面级切割定位）
- **探针 3 轮**（全部撤除；`bop_builder2/arguments.rs`、`bop_split_faces_occt.rs` 复核为未修改、零 `eprintln` 残留）：
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
- **下一轮入口**：`MakeSplitEdges`/pave 顶点身份 —— OCCT `BOPAlgo_PaveFiller_7.cxx:371`（新截面边必须用 DS 的 pave 顶点建立，才能与边界边分割点共享 vertex）、`:589 MakePCurves`。**已核对**：Rust 的边界碎片侧是对的（`pave_blocks/make_blocks.rs:314 make_split_edges` 用 `f.ds().shape(t.v1/v2)` 建边），所以嫌疑收窄到**截面边（FF 新建边）一侧的顶点/端点**：需验证 FF 曲线端点是否复用了边界 pave 的同一 DS vertex（对应 OCCT 的 `BOPDS_DS::Index` 复用 / `IsNewShape`），而不是在同一点再建一个 vertex。探针配方：取一个 `avoid_bnd=1` 的坏面，打印被弃的 1 条边界碎片与 2 条截面边的端点 `shape_key` 与 3D 坐标，看是否"同点不同 vertex"。
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

### 2026-09-21 第 111 轮 · DSH 会话（门禁波次：`9d8e596` ↔ `2c9a66b` 全量 A/B + `iges_check` + 导出）

> 用户指令："导出 obj / 跑门禁" ⇒ 按 §6.2 playbook 步骤 1–5 执行（步骤 6=R2-17、步骤 7=R2-7…R2-10 分诊，按纪律 4 等下一步指令）。

- **步骤 1（工作区确认）**：`git status --porcelain` = 仅未跟踪 `data/iges/`；`git log -1` = `2c9a66b`（批 98）✓。
- **步骤 2（R2-18-gate 差分，§6.1）**：`git worktree add --detach .target-pre85 9d8e596`（已 `worktree remove --force` 撤除；共享树未 `stash`）。两端各跑：`occt-topo --no-fail-fast`（全 16 个 test target + doc-tests）、`occt-geom/geom2d/core/math --lib`。**逐项比对结果：零差异**（唯一变化是 `step_to_obj` **12/13 → 13/13**，属既有红转绿）⇒ 批 85–98 对 live 路径（2D 求交/投影、圆锥盒、`ProjectAct`、`CorrectParameter`、周期曲线族、IGES 写侧）**没有引入任何回归**，也**没有**改变 T-01/T-03/T-04/T-05/phase10 任一红项的数字。日志：`.target-gate/{gate_before,gate_after,after_topo_raw,before_topo_raw}.txt`（`.target-gate/` 命中 `/.target-*/` 忽略，未入库）。
- **步骤 3（全门禁 + IGES + 导出）**：见 §2 新实测表。要点：五 crate `--all-targets` exit 0/0 error；`step_obj_parity` 14/14、`step_obj_area` 11/11、`step_to_obj` **13/13**、`step_geometry_parity` 2/3（T-05）；topo `--lib` **1286/2**；phase20/3/4/5/6/7/8/9 全绿、phase10 7/8、phase19 3/5、boss 1/2；core 290(1i)/geom 143/geom2d 72/math 215(1i)；`iges_check` **18/18 ok**；`export_data_obj` **16/16**。
- **步骤 3 附带（IGES 预期影响，§6.2 第 95/97/98 行的验收）**：`iges_check` 逐模型 A/B —— **DE 数**与**unreferenced 统计**（含类型分布）在 18 个模型上**完全相同**；只有 6 个模型（ATU01038、screw、Shape、Shape-2、occ/bottom、occ/top）的 **P 段长度变化**（= 126 参数块变大：周期/窄区间 B 样条由弦线改精确 126）。与 §6.2 的"仅实体统计变化"预期一致。
- **步骤 4（归因）**：无失配项 ⇒ 无需二分（`.target-pre85` 只用于步骤 2 的 A/B）。
- **步骤 5 报告（两处板内口径订正 + 一条新红）**：
  1. **§2 旧表数字错误**（R2-15 口径问题收口）：旧记"第 90 轮 lib 1287/1、geom 146/146"；实测两端 lib 都是 **1286/2**、geom 都是 **143**。`iges.rs` 在批 84(`9d8e596`)→批 95(`83ed3f8`) 之间无提交改动 ⇒ 旧表的第二条红（`sphere_iges_has_arc_and_solid`）当时也在失败，只是没被记账。
  2. **新登记 T-86**（板内此前未登记的红）：`iges::tests::sphere_iges_has_arc_and_solid` 在批 84 与批 98 **都失败**（单独跑也失败），断言 `iges.contains("100") || iges.contains("128")`；临时探针（已删）dump 出球面 IGES 的 DE 类型 = `[144, 196, 116, 123, 123]`（无 100/110/128/142）⇒ **归因 = T-68 的下游症状**（端口球面 face 没有边界 wire，T-68 探针 `wires/face=[0]`），不是本波次引入，也不应通过改断言"修"。
  3. **旁支（未动手）**：`iges_check` 的 `unreferenced` 在**两端逐模型相同**（Cube 等仅 `{402:1}` = 根自身；但 `Sphere`/`Torus`/`Shape`/`linkrods`/`ATU01038`/`occ/*` 仍列出 `144/196/123/128/198/120/124` 等）⇒ 与 R2-4 卡"孤儿应归零"的预期不符，但**非本波次引入**；下一轮需判定是"统计口径（根 402 未计入引用）"还是"写侧确实漏引用"。
- **未做（按纪律 4 停在报告）**：未执行 §6.2 步骤 6（R2-17 删 `(Circle, Circle)` 快路径）与步骤 7（R2-7…R2-10 分诊）；未改任何源码（本轮只有板内文档改动）。- **导出计数的一处待确认（旁支，已登记未动手）**：本轮 `export_data_obj` 的 16 个计数与 §13 记录**逐项一致**（Cube 24/12、Sphere 642/1244、Torus 1369/2592、Shape-1 3343/4336、linkrods 3494/5078），只有 **ATU01038 = 17,745 v / 22,119 f** 与 §3-P3 表的 09-19 记录（17,767 / 22,160）差 **−22 v / −41 f**；其 bbox parity 仍绿（14/14）。该模型含 `closed_curve=.T.` 的 B 样条，正是 §6.2 第 94 行预告的"读侧周期化会改 `export_data_obj`"影响面 ⇒ 判为**批 94（+96）的预期变化**；本轮**未对 export 做 A/B**（对照 worktree 按 §6.1 已撤除），下一轮如要钉死可补一次 `9d8e596` 对照。**批 98 不可能影响它**（该批只改 IGES/STEP 写侧，OBJ 导出不经过）。

### 2026-09-21 第 112 轮 · DSH 会话（§6.2 步骤 6 = 批 99 ✅（见上）；步骤 7 = T-01/T-03/T-04/T-05 分诊）

> 用户指令：执行步骤 6 与步骤 7。步骤 6 已按 playbook 复跑全门禁（结果逐项相同）并结项 R2-17。本节是**步骤 7 的分诊报告**：每项给出"现状实测（本轮新测）+ 根因/阻塞点 + 可执行下一步"，**未改任何源码**（纪律 4：报告即停）。全部数据来自临时探针 `examples/triage_probe.rs`（**已删**，`git status` 复核无残留）。

**T-01（R2-7）＝ 夹具症状，根因在上游 T-32（mesh→BRep），不在 BOP 链**
- **现状**：`brepfeat::tests::groove_cuts_cylinder` → `grooved volume 8.903438339491371 vs expected 7.708990185052769`（阈值 0.5；removed 仅 0.414 / 应 1.6085）。
- **本轮决定性实验**：把同一 `groove()` 调用喂给**忠实 BRep 圆柱**（`BRepPrimCylinder::make_cylinder(1.0,3.0)`，即 `BRepPrimAPI_MakeCylinder` 的端口对应件）而不是夹具的 `faceted_cyl`（`Solid::wrap(mesh_cylinder(r,h,24))`，mesh→BRep）：`before 9.3754 → after 7.8049`，**removed 1.5705**（解析 1.6085），`|after-(before-groove_vol)| = 0.0380 < 0.5` ⇒ **该断言在忠实输入下通过**。
- **结论**：T-01 的红**不是** BOP/groove 算法的缺陷（13 轮探针结论一致：PaveFiller/BuildSplitFaces/BuilderFace 都忠实），而是夹具经 **mesh→BRep 造出的面 wire 不成链**（T-32：`mesh_to_brep.rs` 未做 `BRepBuilderAPI_MakeWire::Add` 的 `CorrectEdgeOrientation`；round 8 实测 112 面里边全 Forward、source/sink）。**旁支订正**：夹具本身是有效输入（OCCT 布尔能吃），所以**不建议改夹具**，而应按 round 8–9 的结论把 T-32 做**成套**修复（链向 + `BRepBuilderAPI_MakeFace` 面朝向 + shell/solid 装配一致性），修复后本用例应自然转绿。
- **下一步（可执行）**：T-32 一批：`triangulation_to_brep` 补 `CorrectEdgeOrientation` 与面朝向规则，对照回放 `groove_cuts_cylinder` + round 8 的两个回归用例（`groove_negative_volume_delta`、`brepfeat::tests::groove_*` 减材断言）+ 四道 STEP 门禁。

**T-03（R2-8）＝ 输入有效，BOP 装配真的失败（与 T-01 不同族）**
- **现状**：`boss_single_disc_base_merges_one_solid` → `left: 0, right: 1`（Fuse 结果里**一个 solid 都没有**）。
- **关键区分（本轮读测例）**：该测例的工具**不是** mesh→BRep，而是 `single_disc_cylinder(0.25,0.8,24)`（用 `TopoBuilder` 手工装配：24 个侧面 + **整块多边形底盘** + 顶盘）⇒ 输入 wire 天然成链；而**同一断言的 tri-fan 变体 `boss_tri_fan_base_merges_one_solid` 通过**（差异只在底盘/顶盘是"整块多边形"还是"三角扇"）。
- **结论**：T-03 是**引擎侧**缺陷，集中在"底面盘与该 box 顶面**共面**、且盘是大顶点数多边形"的 Fuse 装配：`build_split_solids_full` / `BuilderSolid` 这一路把 solid 全丢了（0 solid），与 T-01 的 wire 不成链**不是**同一根因。→ 与 R2-8 卡原判（"未跨 solid 合并共面片"）一致，但要把范围钉在**共面盘 + 大顶点数环**这一配置。
- **下一步（可执行）**：探测该 Fuse 的 `BuildRC`（`obj_src/tool_src/it_shapes/check_keys/it_exp/tool_sets`）与 `fill_images_solids`/`build_split_solids_full` 的每 solid 镜像数、以及 box 顶面（应变成 2-wire 环）的 `rebuild_split_areas` 块数；对照 `BOPAlgo_BOP.cxx:583-711` 与 `BOPAlgo_Builder_3.cxx::BuildSplitSolids`。若镜像数=1（如 T-01 的早期症状），则瓶颈同样在**面级切割**（`perform_loops`/WireSplitter 对"大顶点数多边形 + 共面重叠"的处理）。

**T-04（R2-9）＝ 断言阈值是端口自造常数，且与端口自身 DS 规模不自洽**
- **现状**：`overlapping_boxes_pipeline_runs_without_errors` → `expected more than the two raw boxes, got 84`；`disjoint_boxes_pipeline_runs` → `expected at least the two boxes, got 88`（断言分别是 `> 2*56` 与 `>= 2*56`）。
- **本轮实测（同一 `PaveFiller`，`set_fuzzy_value(1e-7)`）**：**单盒 `ds.nb_shapes() = 28`**（= 1 solid + 6 faces + 12 edges + 8 vertices + 1）；2 盒相离 = **88**；2 盒重叠 = **84**。
- **两个发现**：① 断言里的 **56/盒** 与端口实际的 **28/盒** 差 2 倍 ⇒ `2*56 = 112` 是**端口自造且不可达**的阈值（两条用例因此恒红，与实现是否忠实无关）；② **重叠(84) < 相离(88)** 属反直觉——重叠本应多出截面边/顶点；说明重叠路径上**有形状没被登记或少登记**（与 T-80/T-81 的"截面边未建/未绑定"同族）。
- **下一步（可执行）**：先按 OCCT `BOPDS_DS::Init`（`BOPDS_DS.cxx`，`Init`/`Append`：按 `TopExp_Explorer` 注册哪些子形状、以及 `Interf` 之后新增的 DS 形状）核对"每盒应登记多少 DS 形状"，据此判定 28 是缺注册还是正确；再对 84/88 做**按 `ShapeType` 的组成分解**（本轮已确认 `count_kind` 可用）定位重叠路径少登记的类别。**禁止**把 56 直接改成 28 了事（那正是规则 3 禁止的"改断言对齐"）。

**T-05（R2-10）＝ 不是 offset 几何/公式问题（T-35 假设被证伪），是"网格朝向 + 体积估计器"问题**
- **现状**：`offset_geometry_is_consistent` → `Offset divergence 2208.0 vs shape_volume 1612.9`，**与 §3 记录逐位相同 ⇒ 批 T-35（A0 offset 曲线公式修复）未改变任何数字** ⇒ 卡内"很可能与 A0 同源"的假设**证伪**。
- **本轮实测（Offset.step，同一个形状，五种独立度量）**：
  | 度量 | 值 | 说明 |
  |---|---|---|
  | `obj_volume(brep_to_obj(sh, 0.01))`（断言的左操作数） | **2208.0** | 散度公式，但**网格有 6 个三角反向** ⇒ 该值无效 |
  | `shape_volume(sh, 0.01)`（右操作数） | **1612.9** | 镶嵌型估计器 |
  | `brep_gprop::volume(sh, 0.01)` | **1665.7** | 另一个镶嵌型估计器（与上面差 3%） |
  | `brep_gprop_full::volume_properties` / `volume_properties_gk` | **1400.0 / 1400.0** | 两个独立精确积分器给出同一个值 |
  | 网格"射线奇偶"体积（与绕向无关，24³ 网格） | **≈2556.6** | 网格 bbox 14³ |
  | 解析（10³ 盒子向外偏移 r=2 的 Minkowski 和，Steiner） | **≈2610.5** | 6 平面 1000 + 1200 + 377 + 33.5 |
- **三条判据**：① 网格按**位置焊接**后 `verts 712→448、edges 1338、boundary(used once)=0、non-manifold=0` ⇒ **网格闭合流形**；② 以内部点 `(5,5,5)` 为顶点的逐三角形有符号体积：**6 个为负**（‑116.67，其余 +116.67）⇒ **OBJ 有 6 个三角绕向反了**（形状是凸的，正确绕向下应全为正）⇒ 断言左操作数 2208 无效；③ 精确积分器的 **1400.0 = (1/3)·Σ(6 个平面面 A·d) = (1/3)·600·7** ⇒ 精确积分**只贡献了 6 个平面面**，20 个曲面面（12 柱 + 8 球）贡献 0 ⇒ 与 **T-68/T-49 同族**（这些面没有边界 wire/pcurve，裁剪面积分路径返回 0）。
- **结论**：T-05 的红**不属于 offset 几何族**，属**网格/属性族**：① 新登记 **T-87**（OBJ 逐面绕向不一致，6 个三角反向）；② 镶嵌型估计器（`shape_volume`/`brep_gprop::volume`）比网格自身奇偶体积低 36%，需单独核对公式；③ 精确积分器对无 wire 的曲面面返回 0（T-68 下游）。测试当前的"两个坏数互相比较"在任一侧修好前不可能收敛；正确 oracle 是解析 2610.5（要不要写进断言 = 需按规则 3 决策，不由本批自定）。
- **下一步（可执行）**：先修 **T-87**（`brep_to_obj`/面朝向：为什么 6 个三角反了——与 round 9 的"面法向·wire 绕向 112/128 为负"同源），再核对 `shape_volume` 的公式，并复核 T-68 落地后 `gprop_full` 的曲面面贡献是否恢复 1200+377+33.5。

### 2026-09-21 第 113 轮 · DSH 会话（**任务状态清理 + 画板监控清单重写**；无代码改动）

> 用户指令："重新清理任务状态，整理未完成或未执行任务，重写画板监控任务列表"。

- **做法（不丢证据的搬迁）**：① **§3 整体重写**为唯一监控清单 —— `§3.1 红门禁`（7 条，逐条给实测数字与归属根因批）、`§3.2 根因批`（8 条，按"修完能让哪条红转绿"排序：T-32 → T-68 → T-69 → T-87 → R2-8 探测 → R2-9 对照 → T-80 链 → T-67 步 3 + T-37）、`§3.3 移植缺口池`、`§3.4 收尾/卫生/需拍板`、`§3.5 归档`；② 旧分块（`R2-xx` 表、`P0–P5` 表、批 1–99 与 T-6x/T-8x 诊断）**原样下沉**为 `§3.6 历史明细（留档）`，并在其开头写明"状态列已过期、监控以 §3.1–3.5 为准"；③ **§14 由 190 行压成 34 行查号表**（A# → 任务 ID → 当前状态指向 §3），§14.1–14.4 四份中途交接（第 40/63/82/89 轮）整体删除（其内容已被 §2 实测与 §3.6 批次明细覆盖）；④ 头部"基线/读法"、§8 维护协议、§10 五问、§11 范围与门禁全部改写为指向新清单；⑤ §3.2 序 3 内嵌的旧 `§14.1 复核表` 引用改为 `§2 实测`（避免把过期数字当基线）。
- **本次清理纠出的状态错误（已在新清单订正）**：① `T-35` 在旧 P0 表仍写 pending，实际 A0 早已 done（见 §3.5 归档）；② `T-01` 旧状态 "in_progress" 与 R2-7 的 "待做" 矛盾 ⇒ 统一为 `blocked（T-32）`；③ `T-77` 旧行仍写 pending，实际已由 R2-1（批 80）结项；④ 旧 P1 的 T-07/T-08/T-09/T-10/T-24（未提交 wave/数据入库/output 归档）在**工作树已干净、`data/occ*`+`data/occ-ref/` 已入库、`output/` 按 §11 不入库**之后已无事项 ⇒ 归档；⑤ **phase10 `curved_face_fillet_sphere_plane`（自 `7178861` 起同红）板内从未登记** ⇒ 新登记 **T-88**；⑥ `T-85`（IGES 可达性）步 1+2 已落地但旧行没写结论 ⇒ 定为 ◐ + 余项（`unreferenced` 未归零，需判口径/漏引用）。
- **统计（第 113 轮末精确复核，取代本行初稿的口径）**：在办**行数 42** = `§3.1` 7 + `§3.2` 8 + `§3.3` 18 + `§3.4` 9；**去重后可追踪 ID 38**（`T-87`/`T-03`/`T-04` 在 3.1 与 3.2 各出现一次、`T-67` 在 3.2 与 3.3 各一次）；`§3.4` 的 9 行里 **4 行是架构级/同族项**（`T-25`/`T-26`/`T-27`/`T-28`，不是额外条目）。归档 **55 个 ID** = 门禁收敛 2（T-02/T-06）+ 覆盖扩展 4（T-21/T-22/R2-14/R2-15）+ A0–A26 已完成 30 + 模块自报 PARK 4（T-12/T-15/T-17/T-31）+ R2 已完成 10（R2-1/2/4/5/17/18/18b/20/21/22）+ 卫生 5（T-07/08/09/10/24）。
- **验证**：仅文档改动（`git status` 只有 `specs/_board.md` 修改 + 未跟踪 `data/iges/`）；未跑代码门禁（无代码改动 ⇒ 门禁数字不会变，最后一次实测仍是第 112 轮的 §2 表）；文件 UTF-8 无 BOM、2008 行、`§3` 各表列数校验通过。
- **下一步**：按 §3.2 序 1（**T-32**）等指令开工。


---

### 2026-09-26 · DSH 会话（round 81→：§3.2 余项 + §3.3 缺口池；提交 `a647836`/`2a28b8c`/`e8deb3d`/`da61a30`/`63ae2b9`）

- **① OCCT ground-truth 能力**：确认工具链可用（VS2022 Pro `vcvars64` + `cl.exe`），既有 GT 探针 `.target-gate/occt_probe/occt_probe.cpp` 重新构建并**跑通**（`--fuse`/`--ds`/`--cylface`/`--cyl`/`--sph`/`--facek`/`--perface`）。运行配方：`set THIRDPARTY_DIR=D:\source\occt-8.0.0\3rdparty-vc14-64 && call env.bat vc14 64`（`env.bat` 的 `THIRDPARTY_DIR` 默认相对路径在本机不成立，必须显式给绝对路径；缺少 `tbb12.dll`/`jemalloc.dll` 时 exe 报 `0xC0000135`）。
  实测：`data/Offset.step` → `BRepCheck valid=1`、`BRepGProp volume=2610.501440`；`box[-1,1]³ ∪ cyl(r=0.4,z∈[0,2])` → **8.50265 / 8 faces**；`--ds`：相离 68=68、重叠 80。
- **② §3.2 T-80 根因三处**（`da61a30`，见 §3.0）：wire 朝向、顶点容差、FClass2d 取 pcurve。GF 结果由「9P+1C」到「7P+2C+2P」，柱面切成 2 rings/2 areas；剩余在 `BOPAlgo_BuilderSolid` 面集合的 avoid/loops（T-82），已记录两处候选控制流语义与「实验版→Fuse 空」的回退证据。
- **③ §3.3 缺口池**：T-44 / R2-6 的 `RationalC1 ∧ U2-U1>=6` 臂（`IncreaseDegree`+`CompCurveToBSplineCurve`，GT 逐字段一致）；R2-18 的 2D 弧长（`GCPnts_AbscissaPoint` 全分派）；T-51/T-67 共同前置 `math_FunctionSetRoot`+`math_FunctionRoot`。
- **④ parity 覆盖**：`data/occ` 的 a3n00/acs10/TDB 断言并入 `step_obj_parity`（14/14）。
- **⑤ 门禁**：见 §3.0 表；**逐项无回归**（含 `export_data_obj` 16/16 逐位计数）。
- **仍红**：`phase10` 7/8（T-88）、`bop_builder2_boss` 1/2（R2-8 定案 (b)）、以及 `--no-fail-fast` 并发整跑时的 `output/` 竞争假红（与引擎无关）。
- **下一步**：T-67 步 3（`Extrema_GenExtPS`）与 T-51 反解（子代理进行中）；R2-19；T-80/T-82 续做。

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
   (a) **最高危 A0 = offset 曲线公式写错**（3D `occt-geom/src/offset.rs:14-19` 用 `p + Offset·Dir` 平移，OCCT `Geom_OffsetCurveUtils.pxx:53-61` 是沿法向；2D `occt-geom2d/src/offset.rs:19-41` 法向 `(-dy,dx)` ↔ OCCT `(dy,-dx)` **反号**；两侧 `d1/d2` 均缺 `DNdir`；3D 版被 STEP 读入 `step/read_topology.rs:1049` 直接使用）⇒ **公式错误而非近似**，先修（并很可能解释 T-05 的 offset 体积偏差）。
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
| 开工前 | 读 **§3.1–3.4**（在办清单）+ §2（门禁基线）；需要背景时再查 §3.6/§7 |
| 每完成一条 | 更新 **§3** 该行（状态/证据/下一步）；若已 done ⇒ 从在办表移入 **§3.5 归档**并写提交号；§7 追加轮次/批次日志 |
| 发现新事实 | 写入 **§5**；在 **§3** 新开一行并续 `T-xx` 编号（**不要**再新建 P0–P5/R2 分块）；若属既有项则并入该项 |
| 遇到错误 | 写入 §9（错误、尝试次数、结论），**禁止重复同一失败动作** |
| 收尾 | 更新头部"最后实测/最后重梳" + §10 五问 + §11 范围；结论是"门禁满足即停，旁支只报告" |
| 状态语义 | `pending`/`in_progress`/`blocked`/`waived`；`done` 必须附门禁输出（新口径以 §2 为准），`blocked` 必须点名阻塞物 ID |

## 9. 错误 / 坑记录

| 现象 | 尝试 | 结论 |
|---|---|---|
| 根目录 `cargo check --workspace` 失败：`could not find Cargo.toml in D:\source\repos\dogs` | 1 | 本仓库**没有根 workspace**，必须 `--manifest-path crates/<crate>/Cargo.toml` |
| `git status` 看不到某文件的改动，但测试数字变了（`builder_solid.rs`） | 1 | 用 `git hash-object` vs `git rev-parse HEAD:<path>` 确认；行为来自共享依赖 |
| 用 `Get-ChildItem output -Recurse` 判"新模型无 Rust 产物" | 1（误判） | 大文件其实是 `output/occ/*.obj`（OCCT 参考）；`output/` 根只有 16+1 个 Rust 产物。判产物要用**非递归**列表 |
| 依赖 mtime 判断"最后动作在 9/18 00:55" | 1 | 全树 mtime 被一次改写刷平，不可作证据 |
| 共享树 `git stash` 造成 44 处冲突标记（历史） | — | 一律改用 `git worktree` A/B |
| 插探针时 `old_string` 跨到循环头，`for` 头被吃掉 → 编译失败 | 1 | 长 old_string 只覆盖"锚点行"；改完先 `cargo test` 编译确认，再继续 |
| **`read` 返回的行号 ≠ 物理行号**（例：`specs/_audit/mesh-exchange.md` 物理 88 行，`read` 报 225 行） | 1（2026-09-20 发现） | `read` 会把超过 ~88 列的长行**折行后逐段编号**（物理第 4 行 260 字符 → 显示 5–7 行），故其行号是**显示行号**；`grep` 输出才是**物理行号**。⇒ ① 报 `文件:行`、定位代码、写探针一律以 `grep` 为准；② 引用他人报告的行号时，若目标文件有长行（CJK 文档、长链式调用），须先 `grep` 复核；③ 已抽查 7 处关键引用（`occt-geom/src/offset.rs:14-19`、`wireframe.rs:392/407-408/462/502`、`step/format.rs:114/121/132`、`model_healer.rs:279-290`、`Geom_OffsetCurveUtils.pxx:53-61`、`BRepMesh_Delaun.cxx:2263-2274`、`BRepMesh_ModelHealer.cxx:491-504`）**均与物理行号一致**，本轮审计的行号未被污染 |
| **`data/occ-ref/` 整目录在会话中途消失**（4 个已提交基线 `T0M/bottom/motoc/top.obj` 报 ` D`） | 1（2026-09-20 发现，未定位） | 出现时机在 T-56 门禁通过之后、仅跑 `git grep`/`git status`/编辑 specs 的指令之间；`crates/` 内**无** `remove_dir_all`/对 `data/` 的删除代码，`git grep occ-ref` 只有 `step_obj_parity.rs` 的**读**引用。已 `git checkout -- data/occ-ref/` 恢复（4 文件回来，`git status` 干净）。⇒ 纪律：**提交用显式路径**（不要 `git add -A`），提交前必看 `git status --porcelain` 里有无 `data/` 的删除；下次若复现，先记录"上一条指令 + 时间"再动手 |
| **T-79 实验：把曲面分支整体迁到 `BOPAlgo_BOP` 后门禁/lib 用例全绿，但几何其实变差** | 1（2026-09-20 第 60 轮） | 门禁（4 道 STEP + phase + export）与 `occt-topo --lib` 都不覆盖该路径的真实几何：会抓到的 `bop_curved::tests::{trimmed_face_box_cylinder_fuse, curved_boolean_full_quadric_unchanged}` **直接调 `curved_boolean_full`**，第一步只改 dispatch 时它们仍走旧体 ⇒ 全绿是假象。⇒ 纪律：**改某函数的分派/实现前，先 `git grep` 该函数的直接调用者，确认哪些现存断言绑在旧体上**；并用探针（box∪cyl 的面类型 + 体积）做判定性对照，别只看测试数字 |

| **批量重命名脚本把新文件名（含 `.rs`）当替换值**（2026-09-20 第 93 轮，批 79） | 1（编译当场抓到） | 映射表存的是 `新名.rs`，但替换模式的匹配段只是旧名 `pNN`（`.rs` 在**前瞻**里、不属于匹配段）⇒ 产出 `mod mesh_queries.rs;`、`poly_roots.rs::f_and_jac`、`foo.rs.rs` 三类畸形，`occt-core` 报 3 个语法错。⇒ 纪律：**批量替换的"匹配段"与"替换值"必须同粒度**（要么都含扩展名，要么都不含）；大批量改名后第一件事仍是门禁 1 编译，且改名后要 `git grep` 三类畸形（`\.rs\.rs`、`(mod|use) \w+\.rs`、`\w+\.rs::`） |

## 10. 五问重启检查（2026-09-21 第 113 轮刷新）

| 问题 | 答案 |
|---|---|
| 我在哪？ | **§3.1 有 7 条红/缺陷**：`--lib` 2 条（T-01 = mesh→BRep 夹具症状、T-86 = 球面无边界 wire）、集成 3 条（T-03 boss Fuse 无 solid、T-04 phase19 DS 断言自造、T-88 phase10 fillet 未分诊）、属性/网格 2 条（T-05 体积度量、T-87 OBJ 绕向）。**§3.2 有 8 条根因批**（最高优先 T-32） |
| 我要去哪？ | §3.2 序 1 **T-32**（mesh→BRep 链向 + 面朝向成套）→ 序 2/3 **T-68/T-69**（球面/带孔面前沿链）→ 序 4 **T-87**（OBJ 绕向 + 体积估计器）→ 序 5/6 **R2-8/R2-9 探测批** → 序 7/8 **T-80 链 / T-67 步 3 + T-37**；缺口池 §3.3 随批夹带；卫生 §3.4 |
| 目标是什么？ | STEP→OBJ 几何/网格对齐 OCCT 8.0.0（DRAWEXE 8.0.0p1 参考），门禁不劣于 §2 基线，无特例补丁、无 OCCT 之外规则 |
| 我学到了什么？ | §5 + 第 111/112 轮：① 批 85–98 对 live 路径**零回归**（唯一变化 `step_to_obj` 12/13→13/13）；② **T-01 与 T-03 必须分族**（前者夹具 mesh→BRep 症状、后者输入有效=引擎缺陷）；③ **T-05 与 A0/T-35 无关**（数字逐位不变），是网格绕向+估计器问题；④ 板内两条断言阈值是自造常数（T-04 的 `2*56`、T-05 的 `2208 vs 1612.9`）；⑤ 旧 §2/§14 的数字口径不可信，一律重测 |
| 我做了什么？ | §7：第 111 轮（门禁波次 A/B + iges_check + 导出）、第 112 轮（批 99=R2-17 + 步骤 7 分诊）、第 113 轮（§3 重写为唯一监控清单、§14 压缩为查号表） |

## 11. 范围与门禁（2026-09-21 第 113 轮重定；不再回问）

**范围（= §3.2 顺序 + §3.3 池）**：① 按 §3.2 序 1→8 修"根因批"（T-32 → T-68/T-69 → T-87 → R2-8/R2-9 → T-80 链 → T-67 步 3 + T-37），每批以"能让 §3.1 的哪条红转绿"为准；② §3.3 缺口池（R2-19 / R2-23 / R2-6 余项 / R2-18 余项 / T-44 / T-51 / T-41 / T-54 / T-62 余项 / T-13/T-14/T-16/T-18/T-19/T-20/T-30）随批夹带；③ §3.4 卫生与需拍板项；④ 架构级项（T-25/T-26/T-27/T-28）先出设计再动手。

**门禁口径**：只引用现有基线 —— `step_obj_parity`(14/14) / `step_to_obj`(13/13) / `step_obj_area`(11/11) / `step_geometry_parity`(2/3) + 五个 crate `--lib`（topo **1283/2**、core 290(1i)、geom 143、geom2d 72、math 215(1i)）+ `export_data_obj`(16/16) + `iges_check`(18/18)；新增 parity 用例只能引用**既有**参考 OBJ；不得劣于 **§2 实测基线**；不为对齐新写算法单测；不加 OCCT 里不存在的规则/阈值；**禁止**为让用例转绿而改断言（须先给出 OCCT 侧依据，如 T-04 的 `BOPDS_DS::Init` 对照）。

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
