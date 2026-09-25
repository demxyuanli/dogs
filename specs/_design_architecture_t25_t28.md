# 架构级四项设计（T-25 / T-26 / T-27 / T-28）

> 依据 `specs/_board.md` §3.4 的要求：「**架构级项（T-25/T-26/T-27/T-28）先出设计再动手**」。
> 本文件只做设计，**不改行为**。每项给出：现状（实测）、OCCT 权威控制流（文件:行号）、影响面、分步方案、验收门禁、风险。
>
> 基线（本文件撰写时实测，HEAD `4db3743`）：topo `--lib` **1287/0**；`bop_builder2_boss` 1/2（R2-8 已定案 (b)）；
> `phase10` 7/8（T-88）；`phase19` **5/5**；`phase20/3/4/5/6/7/8/9` 全绿；`step_geometry_parity` **3/3**；
> `step_obj_area` 11/11；`step_obj_parity` 14/14；`step_to_obj` 13/13；`export_data_obj` 15 模型；`iges_check` Cube/Sphere ok。
> 总纪律：每步先编译；对照 `.cxx` 审控制流；缺分支补忠实移植并注 OCCT 文件:行号，缺失处标 `UNPORTED`；
> 不自创判据/阈值；只跑现有基线；有回归立即回退；提交只用显式路径。

---

## T-25 `GeometryRegistry` 侧表 → 几何进 `TShape`（架构级重构）

### 现状（实测）

- 唯一定义 `crates/occt-topo/src/tgeometry.rs`（612 行）：
  `GeometryRegistry { vertices/edges/faces: RwLock<HashMap<usize, _>>, ids: RwLock<HashMap<usize, u64>> }`（`:162-170`），
  键 = `TShape` 在 `RwLock` 内的堆地址 `std::ptr::addr_of!(*lock) as usize`（`:176-179`）；`global()` 为 `OnceLock` 单例（`:189-197`）。
- 条目：`VertexGeom{point, tolerance}`（`:26-30`）；`EdgeGeom{curve: Arc<dyn Curve>, first, last, tolerance,
  same_parameter, same_range, degenerated, pcurves: HashMap<usize, Vec<Arc<dyn Curve2d>>>, pcurve_ranges}`（`:64-83`）；
  `FaceGeom{surface: Arc<dyn Surface>, tolerance, natural_restriction}`（`:143-147`）。读侧封口在 `brep_tool.rs:24-95`。
- `shape_key`：`tgeometry.rs:473-475`；生命周期 = `TShape::drop → remove_by_ptr`（`tshape.rs:73-82` + `tgeometry.rs:534-547`），
  `ids` 用于防地址复用（`tshape.rs:18-24`）。
- 侧表的自述理由：`tgeometry.rs:1-14`（让 `TShape` 保持轻量、不含 trait object），并已在
  `specs/_coverage.md:19/31/48/81`、`specs/_tasks.md:42`、`brep_builder_full/mod.rs:56-61` 登记为**语义偏差**。
- 半成品：`tshape.rs:84-131` 已有 `VertexShape/EdgeShape/FaceShape` 壳，**生产零消费**（仅 `lib.rs:24` re-export 与 `tshape.rs:142` 测试）。

### OCCT 权威

几何**就在 TShape 子类内**：`BRep_TEdge.hxx:39`（类）、`:66` `Curves()/ChangeCurves()`、`:79-81` `myTolerance/myFlags/myCurves`；
`BRep_TVertex.hxx:33/:47/:64-66`（`Pnt()`、`myPnt/myTolerance/myPoints`）；`BRep_TFace.hxx:51/:59/:141-144`
（`Surface()`、`mySurface/myLocation/myTolerance/myNaturalRestriction`）。曲线参数与 pcurve 在 **representation** 上：
`BRep_GCurve.hxx:62-63`（`myFirst/myLast`）、`BRep_CurveOnSurface.hxx:78-79`（`myPCurve/mySurface`）。

### 影响面（量化）

- 引用 `GeometryRegistry` 的文件 **127 个 / 501 处**；`global()` **279 处 / 88 文件**；`shape_key` **101 处 / 38 文件**。
  `global()` 最多：`shhealing/wire_fix.rs` 22、`sweep.rs` 14、`boptools_2d.rs` 12、`brep_tool.rs` 11、`brep_surface.rs` 10
  （其后 `builder.rs` 9、`step/write_context.rs` 8、`pcurve_full/make_pcurve.rs` 8、`bop_split_seam.rs` 8）。
- 清/拷贝依赖：`shape_ops.rs:24-73`（`transform_shape`）、`:100-218`（`transformed_copy`，`:180-217` 按 `key_map` 重键 pcurve）、
  `:228-259`（`has_geometry`）、`:262-275`（`clear_shape_geometry`）；`clear_shape` **34 处 / 29 文件**，其中多数是
  测试夹具在**共享全局表**上互相清理 —— 迁移后这些夹具应删除。
- pcurve 的**面表回退**是硬依赖：`tgeometry.rs:304-337`、`:364-387` 用 `Arc::ptr_eq(face.surface)` 在同 surface 的另一 face key 上回退；
  删面表后裸 `usize` 无法解析出 surface ⇒ 须按 `BRep_CurveOnSurface.hxx:78-79` 改为**以 surface 句柄为键**；`remove_pcurves_on_surface :439-469` 同理。
- 序列化读表：`bincaf/format.rs:140-243`、`xmlcaf/document.rs:212/293/365`。
- 可行性：`TShape` 已在 `Arc<RwLock<..>>` 内（`shape.rs:12`），`Curve`/`Surface` 已 `Send + Sync`（`occt-geom/src/curve.rs:4`、`surface.rs:8`），
  故 `Arc<dyn Curve>` 可入 `TShape`；`TShape` 已手写 `Debug`（`tshape.rs:39-45`）。

### 分步方案（每步后跑门禁 1–4 + 导出门禁）

1. `TShape` 加 3 个 `Option` 几何槽与类型访问器（复用 `tshape.rs:84-131` 的壳），**不动**现有 API。
2. 把 `GeometryRegistry` 改为 **TShape 槽的转发层**，保留 `set_*`/`xxx_geom`/`edge_*`/`face_*` 全部签名，`shape_key` 仍返回 TShape 地址
   ⇒ 101 处调用**零改动**；跑门禁 1/2/4。
3. 按 `global()` Top5 文件分批把直读点改走 `BRepTool`/访问器（每批跑门禁 1–4）。
4. pcurve/`pcurve_ranges` 键由**面指针**改为 **surface 句柄**（对齐 `BRep_CurveOnSurface.hxx:78-79`），同步 `remove_pcurves_on_surface`。
5. 删 `ids`/`remove_by_ptr`/`TShape::drop` 钩子与 `global()`；`clear_shape` 收敛为 TShape 自身；删除 29 个测试文件里的 `clear_shape` 夹具调用。
6. 全量门禁 + `export_data_obj` 逐位比对（网格数须与基线一致）。

**可停在任一中间步**：若第 5 步无法一次完成，保留 `global()` 作兼容外壳（第 2 步之后它已是纯转发，语义已对齐 OCCT）。

### 风险

- 第 2 步是**语义等价**的转发，但会引入一次全局锁竞争面（`RwLock` 从注册表移到每个 TShape）⇒ 需盯 `--lib` 的耗时（当前 ~13 s）。
- 第 4 步改变 pcurve 的键空间 ⇒ 必须同时改 `transformed_copy` 的重键逻辑（`shape_ops.rs:180-217`，本会话 `2fcec12` 刚修好）。
- `clear_shape` 的 34 处里有测试夹具之间**互清全局表**的隐含顺序依赖（§2 已登记 `groove_cuts_cylinder` 的偶发差异）⇒ 删除它们本身是**净收益**（消掉一类并行抖动）。

---

## T-27 `GetFaceOff` / 角法向：`is_covering_face` 死路径的去留

### 现状（实测）

- `is_covering_face`（`bop_build_solids_leftover.rs:394-403`）自创谓词：「面在实体边界上、至少一条边压在自身分割面上、且边集不与已有面完全相同」；
  唯一调用点 `:285`（`classify_faces_in_solid` 内），无任何测试直接断言它。
- 值链在模块内**确实影响出体**（`:285 → :286-288 → :289-296 in_faces → :107 → :121 → :124 → :135/:151 → :161 → :165 → :188 add_image`），
  但**整条链对生产是死路**：`build_split_solids_full` 的调用者只有合成测试（`bop_build_solids_leftover_tests.rs:263`、`:418`）；
  该模块头注 `:1-11` 自述 *Not called by `BopBuilder`*，声明处 `bop_build_solids.rs:503-504`（`#[path=…] mod leftover;`，**非 pub**、无 re-export）⇒
  grep `leftover::` 生产 **0 命中**。
- 绿路径另有其道：`bop_builder2/arguments.rs:470 fill_in_3d_parts_builder → :472 build_split_solids_occt（bop_split_solids_occt.rs:79）
  → :126/:131 SplitSolid::perform → builder_solid.rs:329-333`。
- 内容来源：`FillIn3DParts`（`BOPAlgo_Builder_3.cxx:97-263`）+ `BuildSplitSolids`（同文件 `:413-618`），并用 `ShellSplitter` 顶替 `BOPAlgo_SplitSolid`（`:17-34`、`:53-54`）。

### OCCT 权威

- `BOPTools_AlgoTools::GetFaceOff` = **角法向**选择：`BOPTools_AlgoTools.cxx:994-1095`（声明 `.hxx:231`）。
  取边中点与切矢（`:1014-1019`）→ 过点垂直切矢的平面（`:1022`）→ 对 `theF1` 与 `theLCSOff` 每对用 `GetFaceDir` 求法向/bi-normal（`:1028-1038`、`:1052-1054`）→
  `aAngle = AngleWithRef(aDBF, aDBF2, aDTF)`（`:1063`；同面退化 π/2π `:1065-1075`）→ 取绕边角序**最小正角**之面为 `theFOff`（`:1088-1092`）；
  角退化（`< Confusion` 或最小角相等）⇒ 返回 `false`（`:1042`、`:1077-1081`）。
- 消费方三条：① `IsInternalFace`（`:939-990`，在 `:977` 调用；`theFace.IsEqual(theFOff)` ⇒ internal，`:983-987`）；
  ② `BOPAlgo_Tools::ClassifyFaces`（`BOPAlgo_Tools.cxx:1622`）→ `BOPAlgo_FillIn3DParts::Perform`（`:1334`）→ `IsInternalFace`（`:1505`），
  再被 `BOPAlgo_BuilderSolid::PerformInternalShapes`（`BOPAlgo_BuilderSolid.cxx:673`）与 `BOPAlgo_Builder::FillIn3DParts`（`BOPAlgo_Builder_3.cxx:201`）调用；
  ③ 直接调用：`BOPAlgo_ShellSplitter::SplitBlock`（`BOPAlgo_ShellSplitter.cxx:359`，函数 `153-620`）、`GetEdgeOnFace`（`:1906`）。`BOPAlgo_BuilderArea` 不直接调用。

### 端口对应件（**角法向已移植且活着**）

- `BuilderSolid` = `crates/occt-topo/src/builder_solid.rs`；`perform :328-336` ↔ OCCT `Perform`（`cxx:76`），
  四段 `perform_shapes_to_avoid(:349↔cxx:129)`、`perform_loops(:430↔cxx:223)`、`perform_areas(:495↔cxx:397)`、`perform_internal_shapes(:591↔cxx:602)`。
- 两条活链：a) `perform_loops :435/:451 → ShellSplitter → shell_splitter_block.rs:359-363 get_face_off`（↔ `ShellSplitter.cxx:359`）；
  b) `perform_internal_shapes :629 classify_faces_occt → bop_classify_occt.rs:289 is_internal_face → algo_tools_face.rs:433 get_face_off`
  （线性版 `bopalgo_tools_class.rs:292`）↔ `ClassifyFaces`（`cxx:673`）。growth/hole（`perform_areas:495`、`is_growth_shell:165`）在 OCCT 侧走
  `IsGrowthShell`/`BRepClass3d`（`cxx:864`），与面偏移无关。

### 方案：**摘除**（不接线）

**理由**：`is_covering_face` 在 OCCT 里**没有对应物**，把它接进绿路径等于引入 OCCT 不存在的规则（违反总纪律）；
而 OCCT 的 `GetFaceOff → IsInternalFace → ClassifyFaces → BuilderSolid` 与 `ShellSplitter` 两条链在端口**都已存在且活着**。

1. 冻基线：`cargo check --all-targets`、topo `--lib`、`fuse_box_cylinder_is_closed_solid`、`step_obj_parity`、`export_data_obj`，记数。
2. **判定性探针**（不改行为代码）：box∪cyl 下打印 `is_covering_face` 收下的 On-边界面集合，与 live `is_internal_face`（`algo_tools_face.rs:479`）对拍。
3. 若一致 ⇒ 删 `:394-403` 与 `:285-288` 的 `covering` 门；`build_split_solids_full` 随之失去与 live 的唯一语义区别 ⇒ 整体摘除 leftover 模块
   （`bop_build_solids.rs:503-504` 的 `#[path]` 声明、`bop_build_solids_leftover_tests.rs`）及仅它使用的私有助手
   （`classify_faces_in_solid`、`connexity_blocks`、`face_state_in_solid`、`collect_all_candidate_faces`、`piece_bbox*`、`BopOp`）。
4. 同批登记 `build_split_solids`（`bop_build_solids.rs:155`）与 `close_open_shells`（`:240`）的去留（调用者同为测试：`:666`/`:698`、`builder_solid.rs:452`）；
   保留则标 `UNPORTED` + 理由。
5. 文档化：在 `algo_tools_face.rs`/`builder_solid.rs` 注释写清上述两条活链的 OCCT 行号；更新 `specs/_coverage.md:33` 与 `specs/_board.md` 对应行。
6. 验证：`cargo check` → 复跑步 1 的全部门禁（逐位一致）→ `graphify update .`；任一计数变动即回退并报告。

### 风险

- 删除会移除 `bop_build_solids_leftover_tests.rs` 里的合成用例（**既有测试**）⇒ 属"删死代码带走的测试"，须在提交信息里明确列出被删用例名与理由（不是"删测试逃避红"）。
- `build_split_solids`/`close_open_shells` 仍被 `builder_solid.rs:452` 等调用 ⇒ 摘除范围必须**限定**在 leftover 模块，不能顺手删。

---

## T-28 ImpPrm 的 HVertex 合并 + `intana`/`intpatch` 重叠合并

> 真实路径：OCCT 在 `src/ModelingAlgorithms/TKGeomAlgo/IntPatch/`（板内写的 `src\IntPatch\` 不存在）。

### 现状与缺口（实测）

- 端口 `intpatch_impprm.rs:57 ImpPrmIntersection::perform` → `:110` 调 `:243 compute_tangency`；该函数**只实现 `IsNew` 分支**，
  `:240-242` 自注 *"Vertex/TopTrans_CurveTransition merge (cxx 329–465) is not-ported: UV-box TopolTool has no HVertex"*。
- 具体缺失：① `PathPoint`（`intpatch_impimp_sonb.rs:14-20`）只有 `p/u/v/param_on_arc/arc`，**无 `is_new`/`vertex`**；
  ② **无同顶点去重** —— `intpatch_impprm.rs:250-287` 每点各产一个 `WalkStart`，`dest[i]` 恒为 `i+1`；
  ③ `TopolTool`（`geom_int_topol.rs:19-25`）只有 UV 盒 + `FClass2d`，**无 `Identical`/顶点表**；
  ④ 全仓无任何 `HVertex` 类型（"HVertex" 只出现在 `intpatch_impprm.rs:242` 的注释里），也无 `TopTrans_CurveTransition` 对应件；
  ⑤ `cxx:1149-1151`/`1255-1257` 的 `SetVertex` 端点绑定同样缺失（`intpatch_impprm_ends.rs:210` 自注 *"UV-box PathPoints are IsNew (one UV)"*）。

### OCCT 权威：`ComputeTangency` 的 `IsNew == false` 分支

- `IntPatch_ImpPrmIntersection.cxx:221-469`（`ComputeTangency`，被 `Perform` 在 `:764`/`:768` 按 D2/D1 域各调一次）；
  本段从 `:329` 起（*"traiter la transition complexe"*，解点落在**已存在的边顶点**上）：
  `:365-418` 扫 `k > i`，`:373 Domain->Identical(vtx, vtxbis)` 判同顶点 ⇒ `:414 Destination(k) = Destination(i)`、
  `:379`/`:382 PPoint.AddUV(...)` 叠加多组 UV；`TopTrans_CurveTransition`（`:333` 构造，`:362`/`:412 Compare`）取 StateBefore/After；
  `:419-454` 定 `ispassing` 与 `vectg/dirtg` 是否反向，`:455-464` 不保留时把 `Destination` 取负回滚 —— 即「**逐顶点去重 + 过渡态 pass/fail**」。
  helper：`TopTrans_CurveTransition`、`IntPatch_TheSOnBounds`、`IntPatch_ThePathPointOfTheSOnBounds`（`IsNew()`/`Vertex()`）、
  `Adaptor3d_TopolTool::Identical/Orientation`、`IntPatch_TheSurfFunction`、`IntSurf_PathPoint::AddUV`。
  等价语义亦见 `cxx:1149-1151`、`1255-1257`（`Perform` 内的 `SetVertex` 绑定）。

### `intana` 与 `intpatch` 的重叠（同一几何两套实现）

- **closed form 双份**：`intana` 侧 `crates/occt-geom/src/intana/analytic_intersections.rs:379`（quadric×plane sphere）、`:404`（sphere×sphere）、
  `:477`（plane×cylinder）、`:540`（plane×cone）、`intana_torus.rs:79`；`intpatch` 侧重写：`intpatch_analytic.rs:69`（plane×sphere）、
  `:92`（sphere×sphere）、`:123`（plane×cylinder）、`:217`（plane×cone）、`:244`（plane×torus）。
- **同一分派器内混用**：`int_face_face_analytic.rs` 的 `plane_sphere:120`/`sphere_sphere:152`/`plane_cylinder:174` 走 intpatch（`:141`/`:163`/`:198`），
  而 `plane_cone:212`（`:232`）、`cylinder_cylinder:281`、`sphere_cone:324`、`cylinder_cone:346`、`cone_cone:368` 走 intana；
  `brep_face_intersect.rs:163-195` 又走一遍 intana（`:177`/`:189`）。
- **plane∩plane 四路**：`intpatch.rs:161`、`int_face_face_analytic.rs:30`、`brep_face_intersect.rs:139`、`intana::quadric_quadric_planes`（`:352`）；
  **plane∩torus 三路**：`int_face_face_helpers.rs:161`、`intana_torus.rs:79`、`intpatch_analytic.rs:244`。
- **两套 surface×surface 分派器**：`intpatch.rs:138 surface_surface_intersection`（解析 arm + 网格 tracer，**不含 ImpPrm**）与
  `intpatch_intersection.rs:41 PatchIntersection::perform`（ImpImp/ImpPrm/PrmPrm）；`int_face_face.rs:414` 先走前者、`:448` 才退后者；
  BOP 侧 `bop_builder_dispatch.rs:412`、`bop_builder_repair.rs:120`/`:511`、`bop_curved/general_mesh.rs:182` **只走前者**。
- **OCCT 权威**：`IntPatch_ImpImpIntersection.cxx` 的 `:9489`/`:7902`/`:8137+8266`/`:8697+9024`/`:9208+9351`/`:3119+3173`/`:3800+9574-9628`
  全部调 `IntAna_QuadQuadGeo`/`IntAna_IntQuadQuad` ⇒ **intpatch 自己没有解析 closed form** ⇒ 合并方向 = **删 intpatch 版、统一到 intana**。

### 分步方案（每步跑：`cargo check -p occt-topo` + 相关套件 + 四道 STEP 门禁）

1. **数据层**：`PathPoint`（`intpatch_impimp_sonb.rs:14`）加 `is_new`/`vertex_id`，`search_on_bounds` 按
   `IntPatch_TheSOnBounds_0.cxx` + `IntPatch_ThePathPointOfTheSOnBounds_0.cxx:34` 填值；`TopolTool`（`geom_int_topol.rs:19`）加顶点表 + `identical`。
   验证：`--lib` 不低于基线（现 **1287/0**）。
2. **平移 `cxx:329-465`**（同顶点合并 + `AddUV` + `TopTrans_CurveTransition` 忠实新增）。验证：`int_face_face_tests`、phase 套件不回归。
3. **补 `attach_wline_ends` 的 `SetVertex`**（`cxx:1149-1151`/`1255-1257`），删 `intpatch_impprm_ends.rs:210` 的 `IsNew` 假设。
   验证：`fuse_box_cylinder_is_closed_solid`、phase5/20。
4. **合并 closed form**：删 `intpatch_analytic.rs:69-260` 的五个重复件，改由 intana 提供，保留 `IntersectionCurve` 适配；
   调用点 `intpatch.rs:77` 的 re-export 与 `int_face_face_analytic.rs:141/163/198` 改接 intana。验证：`intpatch.rs:374-507`、`int_face_face_tests.rs:219-244`。
5. **合并分派器**：让 `intpatch.rs:138` 内部改调 `PatchIntersection::perform`，单一入口覆盖 analytic + ImpPrm + tracer。
   验证：BOP 三处调用点（`bop_builder_dispatch.rs:412`、`bop_builder_repair.rs:120`/`:511`）的 phase 套件。
6. **全门禁**：四道 STEP + phase3–10/19/20 + `export_data_obj` 逐位对齐基线；失配即停并报告。

### 风险

- ImpPrm 在**生产路径可达**（`int_face_face.rs:414 → geom_int_intss.rs:180 → intpatch_intersection.rs:72 → intpatch_impprm.rs:57`；
  `FaceFace` 被 `pave_ff_perform.rs:33`、`int_tools_lines.rs:32` 使用）⇒ 第 1–3 步会直接影响 BOP 结果，
  必须逐模型对拍 `export_data_obj` 与 `step_obj_parity`（本会话 `export_data_obj` = 15 模型，网格数须与基线一致）。
- 第 4/5 步的"参考基线换源"要谨慎：`intpatch.rs:374-507` 的既有断言若与被删实现耦合，应按 goal ④ 以 OCCT 为准订正，而不是改回。

---

## T-26 无面 / 不可用 operand 的判据（`voxel_fallback` 订正）

### 重要订正：标题已过期

`voxel_fallback` **已不存在** —— 第 59 轮 T-41 由提交 **`393ebb8`** 摘除。历史形态（`393ebb8^`）：定义 `bop_builder_core.rs:79`，
调用 `bop_builder.rs:49`、`bop_builder_planar.rs:46`/`:51`；触发条件 `if fa.is_empty() || fb.is_empty() { return voxel_fallback(a, b, op); }`
（任一 operand `faces_of()` 为空）。**今天同一条件已改走忠实引擎**：`bop_builder.rs:50-61`、`bop_builder_planar.rs:50-58`/`:61-69`、
`bop_builder_dispatch.rs:195-197` → `boolean_non_solid:214-231` → `bop_builder2::builder_bop_with_fuzzy`（`bop_builder2/data_access.rs:150-169`）。

残留的体素实现：`boolean_ops::voxel_boolean`（`boolean_ops.rs:173`），唯一 live 调用 = `brepfeat/features.rs:335`（**曲面输入**，属 A12/T-48 另批）；
`solid_union::voxel_union`（`solid_union.rs:237`）只被 `bop_builder_tests.rs:174` 用。旧 warning 串全仓已无。

### OCCT 权威

- 「参数不可用」的判据是 **`BOPTools_AlgoTools3D::IsEmptyShape`**（`BOPTools_AlgoTools3D.cxx:732-741`；`Add:745-786`、`HasGeometry:790-854`）——
  语义是**无几何**（EDGE 无曲线/pcurve/polygon、FACE 无 surface/三角化），**不是"无面"**。
- `BOPAlgo_BOP::CheckData`（`BOPAlgo_BOP.cxx:106-210`）：`:162-167` 空参数 ⇒ `AddWarning(BOPAlgo_AlertEmptyShape)` + `continue`；
  `:181-186`/`:191-201` 维数不合法 ⇒ `AddError(BOPAlgo_AlertBOPNotAllowed)`；`:203-209` 空组继承维数。
  **`TreatEmptyShape`**（`:214-256`）：全空 ⇒ 空结果；一组全空 ⇒ 直接取另一组。
- `BOPAlgo_ArgumentAnalyzer`：`Prepare:115-126` 置 `myEmpty1/2`；`TestTypes:275-352`（空形 `:305-328`）追加
  `BOPAlgo_CheckResult{BOPAlgo_BadType}`；`HasFaulty():261-264`。8.0.0 **没有** `ArgumentAnalyzer::CheckFace`、**没有** `myFaultyShapes`
  （后者在 `BRepAlgoAPI_Check.cxx:89`/`.hxx:159`）；自交在 `TestSelfInterferences:356` ← `BOPAlgo_CheckerSI::CheckFaceSelfIntersection`（`BOPAlgo_CheckerSI.cxx:413`）。
  `BOPAlgo_Builder::CheckData`（`BOPAlgo_Builder.cxx:129-139`）与 `BOPAlgo_PaveFiller.cxx:179` 只报 `AlertTooFewArguments`。

### 端口今天的行为（实测）

- 空 operand **不 Err**，只在 `bop_builder_dispatch.rs:216-221` 加 warning `"BOPAlgo_AlertEmptyShape (objects/tools)"`；
  `bop_bop::check_data`（`bop_bop.rs:71-98`）返回 `"BOPAlgo_AlertTooFewArguments"`（`:73`）/`"BOPAlgo_AlertBOPNotAllowed"`（`:79`/`:82`/`:92`）。
- 真正 `Err("无面")` 的三处（`bop_curved/face_meshing.rs:503`、`general_mesh.rs:209`、`region_trim.rs:337`）经 `dispatch:195` 分流后**对无面 operand 不可达**。
- **关键偏差**：端口 `is_empty_shape`（`bop_bop.rs:42-44`）＝「无顶点且无面」，与 OCCT 的「无几何」**不等价** ⇒ 本设计要钉的就是这一点。

### 分步方案

1. 按 `IsEmptyShape` 语义（`BOPTools_AlgoTools3D.cxx:732-854`）重写 `bop_bop.rs:42-44`（EDGE 无曲线/pcurve/polygon、FACE 无 surface/三角化都算空）。
2. 在 `bop_builder_dispatch.rs:195` / `bop_builder.rs:50` 之前加 `Prepare + TestTypes` 等价检查（空 ⇒ 记 `BadType` 警告，**不是** `Err`）。
3. 把 `TreatEmptyShape`（`BOPAlgo_BOP.cxx:214-256`）语义落进 `boolean_non_solid:214-231`（全空 ⇒ 空结果；一组全空 ⇒ 取另一组）。
4. 摘 `brepfeat/features.rs:329-348` 的体素臂（A12/T-48 另批，与本项解耦）。
5. 改前/改后逐项对拍：`--lib` + phase5/9/10/20 + boss/phase19 + 四道 STEP 门禁 + `export_data_obj`；零失配才继续。
6. T-79/T-80 未解前**只验收到 dispatch 层**，不碰 `bop_curved` 的曲面分支。

### 风险

- 门禁里**没有任何用例传无面 operand**，也没有断言旧 fallback ⇒ 本项属"语义订正 + 清零偏差"，**不会**带来红转绿，但会消掉一处已登记的语义偏差（`specs/_coverage.md`）。
- 唯一体素耦合断言是 `bop_builder_tests.rs:170-177`（与 `solid_union::union_volume` 对拍 < 15%），第 4 步前不要动它。

---

## 汇总：四项的执行顺序建议

| 顺序 | 项 | 类型 | 预计轮次 | 先决条件 | 验收 |
|---|---|---|---|---|---|
| 1 | **T-27** 摘除 leftover 死模块 | 删代码 + 文档化 | 1 | 判定性探针（`is_covering_face` vs live `is_internal_face`）一致 | 门禁逐位不变；被删用例在提交信息中列名 |
| 2 | **T-26** `IsEmptyShape` + `TreatEmptyShape` 语义 | 语义订正 | 1–2 | 无 | 门禁不变；`_coverage.md` 偏差条目消账 |
| 3 | **T-25** 几何进 `TShape` | 架构重构（可分步停在中间态） | 5–8 | 第 2 步的转发层先落地 | 每步门禁 1–4 + `export_data_obj` 逐位 |
| 4 | **T-28** ImpPrm HVertex 合并 + 分派器合并 | 算法移植 | 6–10 | 第 1–3 步各自独立 | 四道 STEP + phase3–10/19/20 + `export_data_obj` |

> 四项都**不需要新测试**（总纪律）；被删/被改的既有断言一律按 goal ④ 以 OCCT 为依据订正。



