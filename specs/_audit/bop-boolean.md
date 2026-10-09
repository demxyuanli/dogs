# BOP/布尔 移植忠实度审查（2026-09-20）

范围：`crates/occt-topo/src/` 的 `bop_*`/`pave_*`/`bopds*`/`int*`/`fclass2d*`/`builder_solid.rs`/`builder_area.rs`/
`shell_splitter*`/`wire_splitter*`/`brep_extrema.rs`/`bop_classify_occt.rs`/`algo_tools*`/`boptools_*`/
`intcurvesurface*`/`bean_face*`/`edge_face*`。只读审查：未改任何源码（`specs/_audit/` 为新建）。OCCT 侧对照 `D:\source\OCCT-src`（V8_0_0）。

## 结论摘要（哪条 live 路径上仍是非忠实实现）

忠实内核**存在**：`bop_builder2::BopBuilder` 走 OCCT `BOPAlgo_Builder`/`BOPAlgo_BOP` 的
`FillImages(V/E/W/F/Sh/Sol) → BuildSplitSolids → FillInternalShapes → BuildShape`
（`bop_builder2/arguments.rs:373-475`、`bop_bop.rs:548`、`bop_split_solids_occt.rs:79`）；点分类用
`brep_class3d::SolidClassifier`（`algo_tools/construct.rs:191`）；同域用 `BopToolsSet`（`bop_tools_set.rs:107`）；
growth/hole 用 `FClass2d::is_hole`（`builder_face_occt.rs:318`）。**但三条 live 入口仍是非忠实的**：

1. **`bop_builder::boolean`/`boolean_dispatch` 家族**（`bop_builder.rs:45`；被 `brep_pattern.rs:81`、
   `boolean_multi/with_check/compound/...` 使用）：BOPAlgo 之前加了**包围盒不相交短路**（现 `bop_builder.rs:88-89`，**仍在**）。
   **✅ 已修**：空面输入的**体素布尔**与**非实体输入的 `boolean_degenerate` 尽力而为**（face∩solid 用**质心**判内外、
   `solid−face` **原样返回 solid**）**均已删除** —— `bop_builder_core::voxel_fallback` 删于 T-41、
   `boolean_degenerate` 整套（含 `degenerate_empty`/`degenerate_non_solid`/`shape_centroid`）删于 T-47；
   空/非实体输入改走 `bop_builder2::builder_bop_with_fuzzy`（= `BOPAlgo_BOP::Perform`+`BuildShape`），
   合法性只由已移植的 `bop_bop::check_data`（`BOPAlgo_BOP::CheckData`，`cxx:140-210`）判定。仅余 `boolean_with_check` 的 tol×10/×100 重算（见发现 16）。
2. **`curved_boolean_full` 曲面分派**（`bop_curved/region_trim.rs:409` ← `boolean_dispatch:215`、`curved_boolean_ext`
   ← `draw/interpreter.rs:422`、`brepfeat/features.rs:320`）：非平面非球面一律走 `general_boolean_trimmed`（样点**多数表决**
   分区 + UV 网格重建 + 空间哈希焊接，`region_trim.rs:42-70,325-401`），球面走解析球冠网格组装（`face_meshing.rs:487`）。
3. **`brepfeat::boolean_feature`**（`brepfeat/features.rs:318`，公开特征 API）：任何曲面输入**直接体素布尔**
   （`features.rs:335`），产出阶梯状 faceted BRep。

忠实内核内部的近似项（**2026-10-06 复核**）：① 面/实体分类 —— `brep_extrema::is_inside`（`brep_extrema.rs:273`）对
**Solid/CompSolid 已改调忠实 `brep_class3d::SolidClassifier`**（T-38），只有**非实体**形状仍走 7×7 网格射线奇偶
（就地标 **UNPORTED**，`is_inside_mesh:290`）；② 点投影 —— 32×32 网格 + 坐标下降**已删除**
（`int_tools_full/context.rs:194-195`，T-52）；③ draft shell 丢弃「顶点<3」的面 —— `bop_draft_solid_occt.rs::face_is_degenerate`
**已删除**（T-47；`face_is_degenerate` 现仅存在于 `bop_builder_repair.rs`，语义不同）。**仍未修**：`bop_curved` 的网格曲面布尔
与 `bop_curved/region_mesh.rs::classify_face` 的 5×5/8×8 采样分类（发现 3/4/7），以及 `boolean_with_check` 的容差重算（发现 16）。

统计：**自创 13 条 + 等价替换未登记 3 条 + 已登记 2 条 = 18 条**（其中 3 条仅在测试/死代码路径）。
> **2026-10-06 复核**：发现 **1/2/5/6/8/9/10** 已修（见各行「处置」行）；发现 **3/4/7/12/16** 仍为自创。

## 发现

### 1. `boolean_degenerate`：非实体输入"尽力而为"（质心判内外、原样返回输入）
- 判定：**自创**
- 证据：`bop_builder_dispatch.rs:572` `if is_inside(b, &shape_centroid(a)) { Ok(empty_result(op)) } else { single_shape_result(a) }`；`:561` `let mut r = single_shape_result(a); r.warnings.push("degenerate cut: 'b' is not a solid; 'a' returned unchanged")`；`:590` `if is_inside(solid_s, &shape_centroid(other)) { single_shape_result(other) }`
- 是否在 live 路径：**是**。`boolean_dispatch:212-213`（任一操作数非 Solid/闭合 Shell）→ `boolean_multi`/`boolean_with_check`/`boolean_compound`/`boolean_cut_many`/`boolean_fold`（`bop_builder_dispatch.rs:232/439/275/651/627`）。
- OCCT 对应：**未找到**。`BOPAlgo_BOP::CheckData`（`BOPAlgo_BOP.cxx:106-137`）只校验参数个数报 `BOPAlgo_AlertTooFewArguments`；OCCT 对 face∩solid 是**真求交**得到被裁剪的面片，face−solid 得到裁剪后的面。
- 影响：**静默产生错误几何**（Cut 面/开壳 → 返回整个未裁剪输入），只写 warning。门禁：`bop_builder_tests_api.rs:158-183`。
- 建议：标为**未移植**并删除该分支——非实体输入交给 `bop_builder2`（`bop_builder2/arguments.rs:251` 只要求 ≥2 个 argument）。

- **处置（✅ 已修：T-47，2026-09-20 第 50 轮）**：`boolean_degenerate`（含 `degenerate_empty`/`degenerate_non_solid`/`shape_centroid`）**整套删除**（2026-10-06 复核：`crates/occt-topo/src` 内已无该符号）；`boolean_dispatch` 的非实体分支改调 `boolean_non_solid` → `bop_builder2::builder_bop_with_fuzzy`（= `BOPAlgo_BOP::Perform`+`BuildShape`），合法性只由 `bop_bop::check_data`（`BOPAlgo_BOP::CheckData`）判定，空实参按 OCCT 记 `BOPAlgo_AlertEmptyShape` 警告。探针实测：`empty∪box → box(vol=1)`、`box−empty → box`、`empty−box`/`empty∩box → 空`、`empty∪empty → Err too few arguments`、`face∪box` 与 `box−face → Err BOPAlgo_AlertBOPNotAllowed`（OCCT 维数规则）、`z=0.5 平面面 ∩ box → 真正裁剪出的面`。

### 2. `voxel_fallback`：空面输入 → 32³ 体素布尔 → faceted BRep
- 判定：**自创**
- 证据：`bop_builder_core.rs:85` `let mesh = crate::boolean_ops::voxel_boolean(a, b, 32, vop)?;`；`:86` `crate::mesh_to_brep::shape_mesh_to_brep(&mesh)`；`:87` `"non-planar input: fell back to voxel boolean"`
- 是否在 live 路径：**是**。`bop_builder.rs:48-50`（任一操作数无 face）；亦作 `bop_builder_planar.rs:46,51` 的平面性门。
- OCCT 对应：**未找到**。OCCT 无体素/网格布尔；空 shape 由 `CheckData` 报错或按 GF 正常处理。
- 影响：**静默产生错误几何**：阶梯状、分辨率相关的实体替换真解（测试自承 "voxel staircase inflates"，`bop_builder_tests.rs:170-177`）。
- 建议：删该分支，标为未移植（空输入返回 OCCT alert/空 compound）。

- **处置（✅ 已修：T-41，2026-09-20 第 59 轮）**：`bop_builder_core::voxel_fallback` 已删除（`bop_builder.rs:3-8` 文件头登记）；`bop_builder::boolean` 与 `bop_builder_planar` 的空/非平面输入改走忠实 `bop_builder2::builder_bop_with_fuzzy`（`bop_builder_planar.rs:47-58`），门禁逐项不变。**余项**：`boolean_ops.rs:173` 的 `voxel_boolean` 仍在，但只作为 port-internal mesh 布尔工具（`bop_builder_planar.rs:46-48` 已写明"OCCT 无对应物"），且不在 `boolean_dispatch` 上；`brepfeat::boolean_feature` 的曲面体素路径仍在（见发现 11）。

### 3. `curved_boolean`/`boolean_mesh`：自创网格曲面布尔（球冠 + 网格分类 + 焊接）
- 判定：**自创**
- 证据：`bop_curved/face_meshing.rs:549` `let mesh = boolean_mesh(a, b, op, tol)?;`；`face_meshing.rs:456` `parts.push(mesh_crossing_face(f, a, keep_in, rim_shared, tol));`；`face_meshing.rs:115` `let inside = point_in_solid(centroid, other, tol);`
- 是否在 live 路径：**是**。`curved_boolean_full`（`region_trim.rs:418`）← `boolean_dispatch:215`/`boolean_degenerate:487`；`curved_boolean_ext`（`general_mesh.rs:298`）← `draw/interpreter.rs:422`（draw `bool`）；`brepfeat/features.rs:320`；`curved_boolean_volume`。
- OCCT 对应：**未找到**。OCCT 曲面布尔 = `IntTools_FaceFace` 求交曲线 → `BOPAlgo_Builder` 分片 → `BOPAlgo_SplitSolid`，全程解析曲面、无球-球特例。
- 影响：**静默产生错误几何**（结果被 faceted 化、丢解析面；面数随 24×24/32×32 网格而变）。门禁：`tests/phase9_integration.rs:39-43`、draw `bool`。
- 建议：标为未移植；曲面分派改走 `bop_builder2::builder_bop_with_fuzzy`，删网格装配。

- **处置（2026-10-06 复核：仍未修）**：曲面网格布尔仍在 `boolean_dispatch` 的曲面分支上（T-79，被 T-80 阻塞——第 60 轮把曲面分支指向 `bop_builder::boolean` 后门禁全绿但探针显示曲面结果错误：`FUSE box∪cyl` 只剩 box/vol 8.0、`CUT` box 不变、`COMMON` 空，故已回退）。**行号已漂移**：现文件为 `crates/occt-topo/src/bop_curved/{face_meshing,general_mesh,region_mesh,region_trim}.rs`。

### 4. `general_boolean_trimmed` + `region_inside_other` 多数表决
- 判定：**自创**
- 证据：`bop_curved/region_trim.rs:65` `} else if inside * 2 >= total { FaceRegion::Inside } else { FaceRegion::Outside }`；`:379` `let reg_class = region_inside_other(f, surf.as_ref(), region, other, tol);`
- 是否在 live 路径：**是**。`curved_boolean_full:413-414`（任何含 `Other`/圆/锥/环面的实体）。
- OCCT 对应：**未找到**。"过半点即 Inside" 不属任何 OCCT 分支；OCCT 由 `BOPAlgo_BuilderFace::PerformAreas` + `FClass2d::IsHole` + `IsInside`（`BOPAlgo_BuilderFace.cxx:441-447,518-527`）做拓扑判定。
- 影响：**静默产生错误几何**：薄壁/临界区域整块误判，切掉或保留整片面。
- 建议：标为未移植。

- **处置（2026-10-06 复核：仍未修）**：`inside * 2 >= total` 计票仍在（现 `bop_curved/region_trim.rs:74`，`:388` 调用 `region_inside_other`），文件内 `:49` 自承 "port-invented"。

### 5. `builder_solid::classify_point`/`shell_inside_solid`/`solid_inside_solid`（live 空腔归属）
- 判定：**等价替换未登记**（实为自创分类）
- 证据：`builder_solid.rs:121` `if is_inside(solid, p) { FaceState::In } else { FaceState::Out }`；`:171` `let Some(p) = face_sample_point(&f) else { return false };`（只取首面中心点）；`:521` `if !shell_inside_solid(&self.loops[i], &solids[j], tol) { continue; }`
- 是否在 live 路径：**是**。`BuilderSolid::perform_areas` ← `bop_split_solids_occt.rs:126-131` ← `bop_builder2/arguments.rs:472`（所有 Fuse/Cut/Common）。
- OCCT 对应：`BOPAlgo_BuilderSolid.cxx:835-860`（`IsInside`）→ `BOPTools_AlgoTools::ComputeState(Face,Solid)`（`BOPTools_AlgoTools.cxx:889`）→ 点分类 `BOPTools_AlgoTools.cxx:790-802`（`BRepClass3d_SolidClassifier::Perform`）。Rust 换算法（网格奇偶）+ 换输入（单点 vs 整 face 边/角方法）。
- 影响：**静默产生错误几何**：空腔归错 growth → 多算/少算内部空腔。门禁：phase20 / boss / 任意带孔 Cut。
- 建议：改用 `AlgoTools::compute_state`（`algo_tools/construct.rs:191`），并按 `algo_tools_face.rs:533-567` 移植 `ComputeState(Face,Solid)`。

- **处置（✅ 主体已修：T-38，2026-09-20 第 44 轮）**：`builder_solid::classify_point` 现为 `crates/occt-topo/src/builder_solid.rs:158`，其 `:167` 改调 `brep_extrema::is_inside`（实体 ⇒ 忠实 `SolidClassifier::classify`，见发现 6），"单点 vs 整 face 边/角方法"的输入差异仅剩在 `face_sample_point` 的首面中心点选择上（`BOPAlgo_BuilderSolid.cxx:835-860` 的 `IsInside` 语义已由分类器对齐）；`algo_tools::compute_state` 已是 F 分支的忠实入口（T-40）。

### 6. `brep_extrema::is_inside`：7×7 网格三角化射线奇偶（#1/#3/#5 的共同根）
- 判定：**自创**
- 证据：`brep_extrema.rs:267` `let meshes = mesh_faces(shape, 7, 7);`；`:283` `ray_plus_x_hits(&verts[*i], &verts[*j], &verts[*k], &jittered)`；`:269-273` 网格化失败时**用包围盒包含兜底**
- 是否在 live 路径：**是**（经 #1/#3/#5：`bop_builder_dispatch.rs:572,590`、`builder_solid.rs:121`、`bop_curved/region_mesh.rs:90,110`、`bop_builder_splitapi.rs:255`）。
- OCCT 对应：`BRepClass3d_SolidClassifier`（`BOPTools_AlgoTools.cxx:797`）；OCCT 无"网格失败→包围盒"兜底。
- 影响：**静默产生错误几何**：边界附近/薄壳/曲面体大量误判；兜底会把整 box 判为内部。
- 建议：删掉该函数在 BOP 路径的用法，改调 `SolidClassifier`（`brep_class3d.rs`）。

- **处置（✅ 已修：T-38，2026-09-20 第 44 轮）**：`is_inside`（`brep_extrema.rs:273`）对 **Solid/CompSolid** 改调 `brep_class3d::SolidClassifier::classify`（`BRepClass3d_SClassifier::Perform`）；原 7×7 网格 + 抖动射线奇偶实现更名 `is_inside_mesh`（`:290`），**仅留给非实体形状**并就地标 **UNPORTED**（OCCT 对 shell 走 `BRepClass3d_SolidExplorer`，端口未接线）。同批修掉同族真实缺陷：`SolidExplorer` 原先缺 `myMapEV`，现按 `BRepClass3d_SolidExplorer::Init`（`cxx:930-982`）只收非 INTERNAL/EXTERNAL 面的、非退化的边及其顶点。证据：`phase3_integration` 3/4 → 4/4、`export_data_obj` 16/16 逐位一致。

### 7. `bop_curved::classify_face`/`point_in_solid`：5×5 UV 样点 + 网格射线构成的面分类
- 判定：**自创**
- 证据：`bop_curved/region_mesh.rs:246` `let (nu, nv) = (5usize, 5usize);`；`:254` `let p = surf.d0(u, v);`；`:258` `if point_in_solid(p, solid, tol) { inside += 1; }`（`classify_face_general` 用 8×8，`:274`）
- 是否在 live 路径：**是**。`classify_face` ← `face_meshing.rs:412-413,506-507`；`classify_face_general` ← `region_trim.rs:356`。
- OCCT 对应：面-实体状态见 `BOPTools_AlgoTools::ComputeState(Face,Solid)`（`BOPTools_AlgoTools.cxx:660-717`，边/角方法）；同一分类在 `algo_tools_face.rs:533` 已有忠实实现，本文件是**重复且不等价的第二套**。
- 影响：**静默产生错误几何**：跨边界曲面片可被整体当 Inside/Outside，网格密度决定结果。
- 建议：标为未移植，改用 `algo_tools_face::is_internal_face`/`compute_state_face_in_solid`。

- **处置（2026-10-06 复核：仍未修）**：`bop_curved/region_mesh.rs:246` 的 `(5usize, 5usize)` 与 `:254` `point_in_solid` 仍在（`classify_face:240` / `classify_face_general:274` 用 8×8）；忠实件 `algo_tools_face::compute_state_face_in_solid` 已存在但仍未接线。

### 8. `IntToolsContext::project_point_on_face`：32×32 网格 + 坐标下降代替 `GeomAPI_ProjectPointOnSurf`
- 判定：**等价替换未登记**（注释只承认平面的问题）
- 证据：`int_tools_full/context.rs:197` `Ok(surface_closest_params(surf.as_ref(), p, 32, 32))`；`brep_surface.rs:238` 双重循环网格扫描；`:252-267` 6 轮四方向坐标下降（非 Newton）
- 是否在 live 路径：**是**。`pave_ef::compute_vf_uv`（`pave_ef.rs:53`）← `pave_vf.rs:73,149,206`（VF 干涉）；亦用于 `bop_build_solids.rs:350,439`。
- OCCT 对应：`IntTools_Context::ComputeVF`（`IntTools_Context.cxx:545-560`）→ `ProjPS` = `GeomAPI_ProjectPointOnSurf`（`IntTools_Context.cxx:247-260`，`Extrema_ExtFlag_MIN`），极值/Newton 精确投影。
- 影响：**静默产生错误几何**：落到错误局部极小 → UV/容差错 → 顶点误并入面（`pave_ef.rs:70` 用 `dist` 放大 vertex tolerance）。
- 建议：改用曲面极值投影（`brep_extrema.rs:203-208` 已用 `point_surface_extrema`）。

- **处置（✅ 已修：T-52，2026-09-20 第 52 轮）**：`IntToolsContext::project_point_on_face`（`int_tools_full/context.rs:195`）已改走忠实 `Extrema_ExtPS`（`ExtPs`），`:194` 就地写明"`surface_closest_params` 兜底已删：`ProjPS` 没有这样的回退"；`brep_surface::surface_closest_params` 的自创 64 点扫描/16×16 栅格/扩窗爬山已删除（求交 half 见 T-52 的行记录）。

### 9. `bop_draft_solid_occt::face_is_degenerate`（顶点<3 丢面）+ 同 TShape 去重
- 判定：**自创**（注释自述 OCCT 会加该面）
- 证据：`bop_draft_solid_occt.rs:41` `crate::topo_tools_full::vertices_of(&face.0).len() < 3`；`:49` `if face_is_degenerate(&fc) { return false; }`；`:54` `if iter_children(&shell.0).iter().any(|c| c.same_tshape(face)) { return false; }`
- 是否在 live 路径：**是**。`build_draft_solid_occt` ← `bop_fill_in3d.rs:75` ← `bop_builder2/arguments.rs:470`（每个源实体）。
- OCCT 对应：`BOPAlgo_Builder_3.cxx:329`（`aBB.Add(aShD, aFx);`）与 `:341-343` **无条件 Add**，无退化门、无 TShape 去重。
- 影响：**静默产生错误几何**：合法 split 面被丢弃 → draft shell 开口/少面 → `SplitSolid` 出错误 pieces；去重可吞掉合法重复面。
- 建议：改回无条件 `builder_add`（对齐 `_3.cxx:329`）；须判退化时用 `BRep_Tool::Degenerated`/面积并登记。

- **处置（✅ 已修：T-47，2026-09-20 第 46 轮）**：`bop_draft_solid_occt.rs::add_draft_face_once` 已删掉自创的"顶点 < 3 即丢面"（含 `face_is_degenerate` 辅助），与 `BOPAlgo_Builder::BuildDraftSolid`（`_3.cxx:329/342/357` 无条件 `aBB.Add`）一致（2026-10-06 复核：`face_is_degenerate` 现仅存在于 `bop_builder_repair.rs`，语义为退化面计数，与此条无关）。门禁 topo lib 1293/1、phase3 4/4、phase6 5/5、四道 STEP 门禁、`export_data_obj` 16/16 逐位一致。

### 10. `edge_edge::find_solutions`：采样求解器作为"补集"并入精确解
- 判定：**自创**
- 证据：`edge_edge/edge_edge.rs:368-371` `for h in edge_edge_intersections(&e1, &e2, tol.max(1e-9)) { if r1.contains(h.u1) && r2.contains(h.u2) { solutions.push((h.u1, h.u2, h.point)); } }`；`:292-303` `is_coincident` 用端点投影判重合
- 是否在 live 路径：**是**。`inttools::edge_edge_intersections` 是采样/近似实现（`inttools/mod.rs:8` "fall back to a sampling+refine solver"），而 `EdgeEdge` 是 PaveFiller 的 EE 求解器：`pave_ee_perform.rs:23`、`pave_force_ee.rs:17`、`pave_intersect/mod.rs:46`。
- OCCT 对应：`IntTools_EdgeEdge::FindSolutions`（`IntTools_EdgeEdge.cxx:290,353`）= 包围盒递归 + 区间细分；`IsCoincident`（`:247`）用 `GeomAPI_ProjectPointOnCurve` 成段比较。
- 影响：**可能静默产生错误几何**：采样补集引入伪交（伪 pave/section 边）或漏窄区间真交（实体不闭合）；重合段缺失使 `CommonPrt`/重合边处理失效。
- 建议：按 `IntTools_EdgeEdge.cxx:290-353` 移植（含 `ComputeLineLine:902`、`FindBestSolution:826`），删采样补集与投影式 `is_coincident`。

- **处置（✅ 主体已修：T-47 删补集 → T-77 移植递归）**：T-47（第 57 轮）删除"把采样器 `edge_edge_intersections` 的解当补集并入"那段，只保留"距离在容差内的极值即交点"；**T-77（批 80 / R2-1）**把 `IntTools_EdgeEdge::FindSolutions` 的 bbox 递归与全部 helper（`DistPC`/`FindDistPC`/`AddSolution`/`FindBestSolution`/`TypeToInteger`/`ResolutionCoeff`/`Resolution`/`CurveDeflection`/`IsClosed`）逐段移植进新文件 `crates/occt-topo/src/edge_edge/find_solutions.rs`（`find_solutions:566`、`find_solutions_range:629`），替代判据与其 Newton polish/`merge_solutions` 去重整体删除。**余项 R2-17**：`perform` 的 `(Circle, Circle)` 快路径 OCCT 不存在。

### 11. `brepfeat::boolean_feature`：曲面输入直接用体素布尔
- 判定：**自创**
- 证据：`brepfeat/features.rs:334` `let resolution = resolution_for(solid, tool, tol);`；`:335` `let mesh = crate::boolean_ops::voxel_boolean(&solid.0, tool, resolution, vop)?;`
- 是否在 live 路径：**是**（公开特征 API：凸台/凹槽/加强筋/拔模）。
- OCCT 对应：**未找到**（`BRepFeat_MakePrism`/`BRepFeat_MakeDPrism` 内部走 `BRepAlgoAPI` 真布尔）。
- 影响：**静默产生错误几何**：曲面特征的体积/拓扑由体素分辨率与自定 `resolution_for` 决定。
- 建议：标为未移植；接 `bop_builder2::builder_bop_with_fuzzy`。

- **处置（2026-10-06 复核：仍未修）**：`brepfeat/features.rs` 的 `resolution_for` 与 `crate::boolean_ops::voxel_boolean` 曲面路径仍在（见 `topo-rest.md` 发现 4/7，行号已漂移：`r.clamp(16, 64)` 现为 `features.rs:483`）。

### 12. `bop_builder::boolean` 包围盒不相交短路
- 判定：**等价替换未登记**
- 证据：`bop_builder.rs:54-56` `if !bbox_a.is_void() && !bbox_b.is_void() && bbox_a.is_out_box(&bbox_b) { return Ok(disjoint_result(a, b, op)); }`；`bop_builder_core.rs:62` 该结果 `faces: vec![], solid: None`
- 是否在 live 路径：**是**（`brep_pattern.rs:81`、`brepfeat`、`curved_boolean*` 的平面分支）。
- OCCT 对应：**未找到**。`BOPAlgo_BOP.cxx` 无 `IsOut`/`Bnd_Box`；`Perform`（`:431` 起）直接 `CheckData` → PaveFiller。
- 影响：几何通常等价（Fuse=compound、Cut=原对象、Common=空），但 `faces` 为空 → 调用方按 face 继续处理会静默丢面；`shape_bbox` 本身是几何采样实现的包围盒，短路条件比 OCCT 乐观。
- 建议：登记为等价替换或删除，让 PaveFiller 自判无干涉。

- **处置（2026-10-06 复核：仍未修；行号已漂移）**：短路仍在 `crates/occt-topo/src/bop_builder.rs:88-89`（`bbox_a.is_out_box(&bbox_b)` → `disjoint_result`），`:3-8` 文件头已登记"Disjoint bounding boxes short-circuit"与 `boolean_planar_legacy` 的 leftover 性质。

### 13. `close_open_shells` + `geometrically_open`：OCCT 没有的 section 面重挂/几何开口丢壳
- 判定：**等价替换已登记**（`builder_solid.rs:36-40` 已写明 translation boundary）
- 证据：`builder_solid.rs:426` `let shells = crate::bop_build_solids::close_open_shells(splitter.shells(), &self.shapes);`；`:431` `if crate::bop_build_solids::geometrically_open(&sh) { continue; }`；`bop_build_solids.rs:284-290` 按"开边集被包含"重挂 loose face
- 是否在 live 路径：**是**（`BuilderSolid::perform_loops` ← `SplitSolid::perform` ← 所有布尔）。
- OCCT 对应：`BOPAlgo_BuilderSolid.cxx:223` 取闭合 shell 后不再重挂；`BOPAlgo_Builder_3.cxx:502-511` 的 IN 面本就成对 FWD/REV 交给 `BOPAlgo_SplitSolid`。Rust 以"共享 section 面 TShape"代替 OCCT "每片各自持一份面"。
- 影响：`geometrically_open` 会**丢弃**几何开口的 shell（OCCT 不会）→ 少算区域；重挂判据为几何包含，可能挂错面。
- 建议：保持登记，补"每源实体 piece 数/总体积恒定"门禁；丢壳改为报 alert 而非 `continue`。

### 14. `bop_geomlib_closed::geomlib_is_closed`：23 点 iso 采样代替 BSpline 极点比较
- 判定：**等价替换已登记**
- 证据：`bop_geomlib_closed.rs:14-15` `* Other / offset / revolution → sampled iso-curves (23 samples, matching OCCT nbp = 23)`；`:17-20` `BSpline/Bezier pole comparison (IsBSplUClosed / IsBzUClosed) is replaced by the same sampled iso test`
- 是否在 live 路径：**是**（`bop_build_faces` 判断 closed edge 的 split 是否为 seam）。
- OCCT 对应：`GeomLib.cxx:2693` `GeomLib::IsClosed`，BSpline 分支为极点比较（精确）。
- 影响：可能漏"两端 iso 只在极点相接"的闭合（漏 seam 分裂）或误判（伪 seam 多切面）→**改变拓扑面数**。
- 建议：登记；待 `dyn Surface` 暴露极点/`IsUClosed` 后改回精确判据。

### 15. `bop_split_seam`：`Geom2dAPI_ProjectPointOnCurve` 的 AABB/采样替代
- 判定：**等价替换已登记**
- 证据：`bop_split_seam.rs:7-8` `//! AABB / sampling stand-in for Geom2dAPI_ProjectPointOnCurve. Dual p-curves are stored with GeometryRegistry::set_edge_pcurves.`
- 是否在 live 路径：**是**（`BOPTools_AlgoTools3D::DoSplitSEAMOnFace`，由 `bop_build_faces` 调用）。
- OCCT 对应：`BOPTools_AlgoTools3D.cxx:275` `Geom2dAPI_ProjectPointOnCurve aProjPC1, aProjPC2;`
- 影响：seam pcurve 的 2D 参数偏移（AABB 近似）→ seam 边分裂参数偏；曲面本身不变，拓扑点参数变。
- 建议：登记；后续接 `occt_geom2d` 投影 API。

### 16. `boolean_with_check` 的 tol×10/×100 重算 + `heal_tolerance` 的 8·tol 滤边
- 判定：**自创**
- 证据：`bop_builder_dispatch.rs:444-446` `for &scale in &[10.0f64, 100.0] { let nt = tol * scale; let candidate = boolean_dispatch(a, b, op, nt)?; ... }`；`bop_builder_splitapi.rs:340` `crate::shhealing::remove_small_edges(&s1, tol * 8.0)`
- 是否在 live 路径：**否（公开 API，不在门禁/导出路径）**。`boolean_with_check` 仅被 `bop_builder_report.rs:369` 与测试调用。
- OCCT 对应：**未找到**。"跑三遍取问题最少者"不是 OCCT 控制流（OCCT 只有单一 `myFuzzyValue`）；`8·tol` 不是 `ShapeFix` 参数。
- 影响：若被外部调用会**静默给出与请求容差不一致的几何**（以 100×tol 算出的结果）。
- 建议：标为未移植，或限定为容差修复工具并禁止默认使用。

- **处置（2026-10-06 复核：仍未修）**：`bop_builder_dispatch.rs:505` 的 `for &scale in &[10.0f64, 100.0]` 重算仍在；`bop_builder_splitapi.rs` 的 `heal_tolerance` 自创修形链仍在（`weld_coincident_vertices`/`remove_small_edges`/`remove_degenerate_faces`）。该函数不在门禁/导出路径（只被 `bop_builder_report.rs:369` 与测试调用）。

### 17. `bop_build_solids_leftover`：包围盒量化去重 + 已删的 `merge_sharing_faces`
- 判定：**自创**（模块自述 UNPORTED）
- 证据：`bop_build_solids_leftover.rs:176-186` `UNPORTED: a merge_sharing_faces pass used to run here. ... OCCT has no such stage`；`:452` `let q = |x: f64| (x / 1e-6).round() as i64;`；`:135` `if piece_bbox_volume(&p) < 1e-9 { continue; }`
- 是否在 live 路径：**否（死代码，仅测试）**。`bop_build_solids_leftover_tests.rs:263,418`；模块头 `:3-6` 自述 "Not called by `BopBuilder`"。
- OCCT 对应：`BOPAlgo_Builder_3.cxx:579-616`（只有 `aMST` face-set intern，无合并、无包围盒去重）。
- 影响：无（死代码）；但 `piece_center`（`:465`）按包围盒中心分类、包围盒体积丢件是第二批自创规则，若误接线会截断合法区域。
- 建议：保留 UNPORTED 标注，或删除以免误用。

### 18. `bop_build_faces::group_wires_as_areas`：用"严格包含"替代 `FClass2d::IsHole`
- 判定：**自创**
- 证据：`bop_build_faces/builder_like.rs:197-203` `/// FClass2d::IsHole is deliberately *not* used for the growth/hole decision`；`:227` `if j != i && loop_contains(oes, &wire_edges[i]) { hole = true; break; }`
- 是否在 live 路径：**否（仅测试）**。调用点仅 `bop_build_faces/mod.rs:336-388` 两个 `#[cfg(test)]`；live 的是 `builder_face_occt.rs:292`（`:318-319` 正确用 `FClass2d::is_hole`）。
- OCCT 对应：`BOPAlgo_BuilderFace.cxx:441-447`（`IsGrowthWire` + `aClsf.IsHole()`）。
- 影响：无（live 已用忠实实现）；风险是这套"故意不用 FClass2d"的实现被接线回主路径，会把 growth/hole 判成"嵌套即孔"，对分离片/相切环出错。
- 建议：删除该实现，或让测试改调 `builder_face_occt`。

## 抽样核实为忠实的区域（无自创判据）

- `bop_tools_set.rs:107-123`（退化边跳过、INTERNAL 存两次）与 `bop_same_domain_faces.rs:279-304`
  （planar-bounded 快捷、`aFaceToParent` 传播、代表选取）对齐 `BOPAlgo_Builder_2.cxx:562-925`。
- `bop_classify_occt.rs`（AABB 取代 `BOPTools_BoxTree`）、`bop_connexity_faces.rs`（`MakeConnexityBlock`）、
  `bop_pair_sd.rs`（`PairOfShapeBoolean`/`VFI` 串行化）属**已登记**等价替换（模块头均写明）。
- `bop_cells.rs`（`BOPAlgo_CellsBuilder` 索引/`FindParts`）、`bop_pair_selector.rs`、`bop_box2d_tree.rs`、
  `bop_aabb_faces.rs`、`builder_face_occt.rs` 未见自创判据。
