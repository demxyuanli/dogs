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
- `bop_builder::face_face_segments_local`（:277）虽做两面 2D 多边形裁剪，但：`face_polygon_local` 是**顶点-only、凸-only、质心角排序**近似；无 3D 中点重验证；共线边分支（:234）无 inside 测试 → **杂散段**（box 侧面 × 圆柱侧片）。
- OCCT 做法：`IntTools_FaceFace` 裁剪到 **UV 域**（`CorrectSurfaceBoundaries`）+ 稀疏点测试（`dom.Classify`），靠 `PutBoundPaveOnCurve`（面边界穿越点）+ `IsValidBlockForFaces`（块中点 in-face 过滤）兜底。

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

### 波 2 — 有限段/共用顶点分裂（G2，BOPAlgo pave 分裂，~400 行）
- `split_polygon_by_segment` 改**有限段**切割（交点在段内才切），或接入 BOPDS pave 分裂：边在**两面顶点并集**处切分（`PutBoundPaveOnCurve` + `PaveBlock::Update`）。
- 门禁：box 顶面只切环+盘，外边界不被弦延长线细分。

### 波 3 — 平面布尔合并到 BOPDS/BOPAlgo 管线（G3+G4，~800 行）
- `bop_builder::boolean` 的 Step 6 替换为 `bop_builder2`/`bop_build_*` 装配（`fill_images_faces` + `build_split_solids` + `close_open_shells`），获得：跨面边统一（`shapes_sd`/CommonBlock）、几何 shell 闭合（`ShellSplitter` EKey）。
- 分类/选择逻辑（`classify_face`/`select`）保留，输出改走 BOPAlgo 重建。
- `shell_is_closed` 增加**几何身份**闭包检查（`shell_splitter::EKey`）作为兜底，不再只看 TShape 指针。
- 门禁：`boss_adds_material` 通过（shell 闭合、体积 ≈ 4.157±0.6）；现有布尔测试全绿。

### 波 4 — 拓扑不变量 oracle（trellis R4，~300 行）
- 实现/接入验证门禁：`shell_is_closed`（几何）、`shell_manifold_check`、`shell_euler_characteristic`、`bop_checker` 非自交。
- 扩展差分 harness（engine/occt/ops.json）加布尔不变量 op；`differ.py` 加 `invariant` oracle 策略。
- 门禁：每波提交前全量 `cargo test` + 不变量绿。

### 波 5 — 曲面布尔回归（bop_curved，~400 行）
- `bop_curved` 路径走同一装配（FaceRegion 分类 + BOPAlgo 重建）；rev/曲面 fuse 的 shell 闭合。
- 门禁：rev、曲面布尔测试全绿。

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
