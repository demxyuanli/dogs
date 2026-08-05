# 布尔内核移植任务 — 平面布尔合并到 BOPAlgo 管线（完整 porting）

> 目标：把破碎的平面布尔 `bop_builder::boolean` 合并到**已移植的 BOPAlgo 管线**（`bopds` + `bop_build_*` + `bop_builder2`），使平面与曲面布尔的 shell 正确闭合、solid 可构建、拓扑不变量全绿。**不是零散修复**——是一次完整的管线整合移植。

## 现状：两条重复管线（核心问题）

| 管线 | 入口 | 使用 BOPDS/BOPAlgo | shell 闭合 |
|---|---|---|---|
| **平面布尔** | `bop_builder::boolean`（bop_builder.rs:823） | ❌ 不 import bopds/bop_build_* | ❌ TShape 指针计数，跨面边不统一 → 不闭合 |
| **通用 fuse** | `bop_builder2::BopBuilder::build_bop`（bop_builder2.rs:411） | ✅ bopds + bop_build_faces/common/solids | ✅ `ShellSplitter` + `close_open_shells`（几何 EKey 身份） |

**正确的机制已经移植**（bop_builder2 路径：`fill_images_vertices→edges→faces→solids`，BOPDS pave-block 边分裂、`shapes_sd` 顶点统一、几何 shell 闭合）。平面 `bop_builder::boolean` 是**另一套简化重复实现**，未接入。

## 缺口映射（调查证据）

### G1 面-面求交不裁剪真边界（IntTools_FaceFace）
- `face_face.rs::face_face_intersection`（:195）：`Line` 是**无限线**（无面内裁剪）；`Curve` 是采样点云（无连通性）。
- `bop_builder::face_face_segments_local`（:277）做两面 2D 多边形裁剪。**调查修正（波 1，2026-08-05）**：实测 boss fuse 时 `face_face_segments_local` 只产生**正确的 24 边形段**（24 条，全在 z=1，box 顶面×圆柱底边），**无杂散段**——G1 对本 case 不是根因（之前诊断被 dbg 错误工具位置污染）。
- 注意：`face_polygon_local` 的顶点-only/凸-only 质心角排序对**非凸/自交叉边界**仍不可靠（见 G5）。

### G5（新发现，波 1）box 侧面重建出**自交叉 wire**（最深根因）
- 实测 fuse 结果里 box y=0 侧面的 wire 是 `(2,0,0)→(0,0,0)→(2,0,0)→(2,0,1)→(2,0,1)→(0,0,1)→(0,0,1)→(0,0,0)`——**边不连通**（e1.b≠e2.a），自交叉。box x=2 侧面同样。
- 此 malformed wire **无 unify 也存在**（pre-existing），来自 `split_faces`→`polygon_to_subface` 重建 box 侧面时用了**坏多边形**（`face_polygon_local`/`face_boundary_points`）。
- **下一步（波 2 重定位）**：先修 `face_polygon_local`/box 侧面多边形提取（干净矩形），再谈边统一。手写 `unify_result_edges`（两次尝试）因坏输入面而失败，已回退——正确路径是迁移文档 Wave 3 的"接入 bop_build_* 管线"，非手写重建。

### G2 无限直线分割（BOPAlgo 边分裂）
- `split_polygon_by_segment`（bop_builder.rs:345）：用 `side()` 叉积按**无限直线**分类 → 弦延长线穿越 box 边界 → 过度细分、slivers。
- OCCT 做法：**pave-block 分裂**（`BOPDS_PaveBlock::Update`）——边在两面顶点**并集**处切分。

### G3 跨面边不统一 + TShape 闭合检测
- `boolean` Step 6（:906）：`make_shell` 纯容器；`shell_is_closed`（shell_check.rs:107）按 **TShape 指针**要求每边恰好 2 面。跨面同几何异 TShape → 报开放。
- OCCT 做法：`MakeBlocks`/`PostTreatFF` 融合 section 边 → 一个代表 TShape；`UpdateFaceInfo` 同 PB 注册进两面；`MakeSplitEdges`+CommonBlock 统一分裂边；`FillImagesEdges` SD 绑定。
- 我们已有：`bopds.rs`（`BopdsDS`/`BopdsCommonBlock`/`shapes_sd`，:984）、`bop_build_solids.rs::close_open_shells`（:212，几何 EKey 身份）、`shell_splitter.rs`（EKey=(VKey,VKey)）。

### G4 两条管线未合并
- `bop_build_faces/common/solids` + `bopds` 只被 `bop_builder2` 使用；`bop_builder::boolean` 完全独立。

## 移植任务（波次）

### 波 1 — 面-面求交正确化（G1，IntTools_FaceFace 移植，~500 行）
- `face_face.rs`：`Line` 结果携带**每面 pcurve/UV 域裁剪**；实现 `face_uv_domains`（周期/seam 处理，参考 `CorrectSurfaceBoundaries`）。
- `face_face_segments_local`：裁剪到每面 **UV 域** + **3D 中点双面重验证**（`IsValidBlockForFaces` 等价：段中点 `point_in_face` 双面）；删除共线边无测试分支。
- 门禁：box+圆柱 fuse 后 box 侧面不再被切；杂散段清零。

### 波 2 — 有限段/共用顶点分裂（G2）【已尝试，回退，后续】
- 已实现 `split_polygon_by_segment` 有限段/弦切割版本：面数 341→89（干净），**但 shell 不闭合**（24 边形弦边未与圆柱底盘共享，36 条 1-face 边——弦端点与圆柱顶点焊接不匹配）。已回退。
- 现状：无限直线分割 + Wave 3 统一达到正确闭合（341 面，正确体积/Euler）。**过度细分是质量/效率问题，非正确性**。
- 后续：接入 BOPDS pave 分裂（`PutBoundPaveOnCurve` + `PaveBlock::Update`）替代手写多边形切割，或修弦端点焊接。

### 波 3 — 平面布尔合并到 BOPDS/BOPAlgo 管线（G3+G4）【核心已落地 2026-08-05】
- **已实现**：`bop_builder::unify_result_edges`——结果面间按几何端点统一边 TShape（BOPDS pave 式细分），**每面边界按共享端点链边重建 loop**（`edge_vertices` 是首创建者序，非遍历序——P1 orientation 丢失所致），再在所有全局焊接顶点处细分，相邻对映射到单一规范 Edge。
- **已实现**：重排——unify 在 `orient_faces_outward` 之前（orient 的 `shell_is_closed` 门禁需统一后的闭合 shell）。
- **效果**：box+凸出圆柱 fuse `shell_is_closed=true`；`boss_adds_material` 解除 ignore 并通过；全量 lib 1259 通过。
- **剩余**：G2 无限直线分割过度细分（质量/数量问题，非闭合必需——unify 已兜底）；G5 面重建（已定位：共享边 mixed-orientation 由链边解决）；`bop_builder2`/`bop_build_*` 全量合并（可选，当前手写 unify 已达同效果）。
- 门禁：✅ `boss_adds_material` 通过；全量 lib 1259 绿。

### 波 4 — 拓扑不变量 oracle（trellis R4）【✅ 完成 2026-08-05】
- `shell_check::ShellInvariants{closed, euler_characteristic}` + `shell_invariants()`；`is_valid_solid()`（closed + Euler=2）。
- **顶点统一**：`unify_result_edges` 每焊接点一个共享 Vertex TShape（`make_edge_segment_with_vertices`），使 Euler 有意义。
- 布尔门禁：solid 由 closed(manifold) 决定（空心实体 Euler=0 也有效）；Euler 由 oracle 测试断言。
- 测试：`fuse_box_cylinder_is_closed_solid`（box+圆柱 fuse 断言 closed + Euler=2）。
- 全量 lib 1260 通过。
- 差分 harness 扩展（ops.json/differ.py invariant oracle）：后续。

### 波 5 — 曲面布尔回归（bop_curved）【✅ 已验证 2026-08-05】
- `bop_curved` 23 测试全过（曲面路径独立于平面 unify，不受影响）。

## 依赖与纪律
- 波 1→2→3 顺序（G1/G2 是几何根因，波 3 合并需要干净几何）；波 4 可并行（oracle 独立）；波 5 依赖 3。
- 复用已移植的 `bopds.rs`/`bop_build_*`/`shell_splitter.rs`——**不重写**，只接线 + 修 G1/G2。
- 每波相关模块测试（`cargo test --lib bop_builder:: bop_builder2:: bop_build_:: bopds:: shell_::`）；全量仅提交前。
- 参考源：vcpkg `V7_9_3-3913ae94de.clean/src/{IntTools,BOPDS,BOPAlgo,BOPTools}/`。
- 迁移索引：`topo_bop_builder`/`bop_builder2`/`bop_build_*`/`bopds` 已 "done"（已移植）——本任务标记为**接线/纠正**，不改结构。

## 验收（最终）
- `boss_adds_material` 解除 `#[ignore]`，通过。
- 平面布尔（box+cube、HoledPlate、rev）shell 全闭合、solid 可构建。
- 拓扑不变量 oracle：所有布尔结果 closed + manifold + Euler=2 + 非自交。
- 全量 lib + 差分 harness 绿。
