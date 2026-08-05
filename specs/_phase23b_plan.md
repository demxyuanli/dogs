# Phase 23b Plan — 布尔结果 shell 闭合：跨结果面的边统一（BOPAlgo 边合并）

> 审查结论（已复现并定位根因）：
> - **现象**：`bop_builder::boolean` 对 box + 突出 24 片刻面圆柱的 Fuse，结果 shell 不闭合（`warnings: ["result shell is not closed"]`，exact 体积 4.674 vs 实际 3.88）。boss_adds_material 因比被 ignore。
> - **开放位置**：`shell_manifold_check` 报告 16 条边被 1 个面引用，**全部在 z=1**（box 顶面 = 圆柱底面平面）——混合 box 顶面外边界段（`(2,0,1)->(0,0,1)` 应与 y=0 侧面共享）与 24 边形细分段。
> - **根因**：`split_faces`（bop_builder.rs:503）给分割碎片经 `polygon_to_subface` → `edge_map.edge`（按焊接点索引 (min,max) 去重，L483-497）**新建边**；但**未分割面（box 侧面、圆柱侧片）保留原始 TShape 边**。两组边几何位置相同但 TShape 不同 → `shell_is_closed`（按 TShape 计数）不共享 → shell 开放。
> - OCCT 参考：`BOPAlgo_Builder::Build` 经 `BOPDS_DS` 在**整个结果**上做边/顶点统一（intersection edges + original edges 合并为一个共享 edge），非局部去重。

## 根因证据链（Wave A 调查深化，2026-08-05）
1. 复现：`boolean(box, faceted_cyl, Fuse, 1e-4)` → `shell_is_closed: false`，16 条 1-face 边全在 z=1。
2. **三层根因**（逐一实测）：
   - **#1 杂散分割段**：`face_face_segments_local` 对 box 侧面 × 圆柱侧片算出的交段，来自**面片平面无限延伸穿越**——圆柱在 (1,1) 半径 0.25 完全不触 box 侧面，但侧面被切成碎片（y=0 面被切出三角形缺 (0,0,0)）。→ **先修 #1**。
   - **#2 无限直线分割**：`split_polygon_by_segment` 用 `side()` 叉积按**无限直线**分割，24 边形每条弦的延长线都穿越 box 边界 → box 顶面被过度细分，环外边界与侧面细分不一致。
   - **#3 无跨结果面边统一**：`edge_map` 只在 split_faces 内部去重；跨面（环碎片 vs 侧面、环 vs 圆柱底盘）不共享边。
3. 已尝试 #3（`unify_result_edges`，BOPDS 式边细分+统一）——**单独不足且自身脆弱**（重建 loop 丢顶点、退化环），已回退。**必须从 #1/#2 修起**：#1/#2 消除后几何干净，#3 才能兜底。

## 移植任务（波浪式，修订）

### 波 A0：修杂散分割段（#1，~400 行）【先做】
- `face_face_segments_local`（bop_builder.rs:277）：交段必须**同时落在两面的实际多边形内**（段端点投影到两面的 2D 多边形做 point-in-polygon），面片平面无限延伸产生的段要剔除。
- 门禁：fuse 后 box 侧面不再被切（y=0 面保持矩形）；offending 边数从 16 下降。

### 波 A：修无限直线分割（#2，~600 行）
- `split_polygon_by_segment` 改为**按有限段**分割（或对切割多边形做正确布尔），box 顶面只切成环+盘，外边界不再被弦延长线细分。
- 门禁：环外边界与侧面全等（可共享边）；offending 边数显著下降。

### 波 B：跨结果面边统一（#3，~800 行）
- 在 `boolean`（`orient_faces_outward` 后、`make_shell` 前）做 BOPDS 式边统一（几何端点对→共享边 + 边细分）。
- 门禁：上述 fuse 后 `shell_is_closed == true`；`shell_manifold_check` 无违规边；体积 ≈ 4.157±0.6。

### 波 C：回归（~400 行）
- 全量 lib 回归；现有 boolean 测试全绿；boss_adds_material 解除 `#[ignore]`；新增 shell 闭合+体积断言。

## 依赖与纪律
- **A0 → A → B → C 顺序**（#1/#2 是几何根因，B 只在干净几何上兜底）。
- 每波相关模块测试（`cargo test --lib bop_builder:: bop_build_faces::`）；全量仅提交前。
- 参考源：OCCT `BOPAlgo_Builder.cxx` + `BOPDS_DS.cxx` + `IntTools_FaceFace`。
