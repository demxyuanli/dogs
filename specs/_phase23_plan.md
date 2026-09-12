# Phase 23 Plan — TShape.children 存完整 TopoShape（orientation+location 语义对齐 OCCT）

> 审查结论（已核实 OCCT 8.0.0 源码，rustocc/target/occt-8_0_0_rev3.../inc + vcpkg v7.9 src）：
> - **存储**：`TopoDS_TShape.hxx` → `private: NCollection_List<TopoDS_Shape> myShapes`。children 是**完整 TopoDS_Shape**（TShape+location+orientation）列表，不是 TShape 句柄。
> - **Add**：`TopoDS_Builder::Add`（TopoDS_Builder.cxx）`L.Append(aComponent)` 追加完整组件；**父 REVERSED → 存子反向**（`S.Reverse()`）；**父位置非恒等 → 存子位置 = 组件位置 × 父位置⁻¹**（`S.Move(aLoc.Inverted())`）；并 `aComponent.TShape()->Free(Standard_False)` 冻结组件。
> - **读**：`TopoDS_Iterator::Initialize` 遍历时 `myShape.Orientation(TopAbs::Compose(父方向, 子存储方向))` 复合父方向。
> - 我们的 `builder::add`（builder.rs:25）只 `t.add_child(sub.tshape.clone())` → orientation **和 location 全丢**。

## 现状与根因

`TShape.children: Vec<HandleTShape>`（tshape.rs:19）；所有读取端 `TopoShape::from_handle`（shape.rs:19，恒 Forward）重建 → 同 TShape 方向信息全灭。

**实测后果**（Cylinder.step 侧壁 wire）：4 边 3 唯一 TShape、全 Forward（#77 Reversed 塌缩）。`map_shapes` 去重键 `(TShape, orientation)`（topo_tools_full.rs:31）**设计已对**，但数据恒 Forward → 去重失效。STEP 写出硬编码 `.T.`（step.rs:945）。`pc_curve_orientation`（pcurve_full.rs:553）已实现但未接入 `make_pcurve_full`。

## 实现范围（已量化）

| 面 | 规模 | 位置 |
|---|---|---|
| 写侧 | 6 文件 / 11 `add_child` | builder.rs、brep_assembly.rs、shape_ops.rs、bbox_from_geometry.rs、xcaf.rs、tshape.rs |
| 中央读取器 | ~7 文件 | topo_tools_full、topexp、iterator（已是 Vec\<TopoShape\> 设计）、shell_check、shape_checks、wireframe |
| 长尾读取 | 44 文件 / 168 `.children` / 39 文件 94 `from_handle` | 机械改，编译器逐个暴露 |
| 新行为消费 | ~4 处 | step writer、make_pcurve_full+pc_curve_orientation、bop_build_faces loop_signed_area、bopds |
| 已就绪基建 | — | map_shapes 去重键、pc_curve_orientation、iterator.rs API |

**选型：方案 A `children: Vec<TopoShape>`**（OCCT 证实，见上）。A/B 成本相同（都改 168 读点），A 多修 location + 匹配 iterator.rs + 读取端更短。

## 波次

### 波 A：结构 + 写侧（tshape.rs + 6 写文件）
- `children: Vec<HandleTShape>` → `Vec<TopoShape>`；`add_child(c: TopoShape)`。
- `builder::add` 按 OCCT `Add` 语义：存完整 sub；**父 Reversed → 子反向**；**父位置非恒等 → 子位置 = 子位置 × 父位置⁻¹**；**冻结子（Free(false)）**。
- 验证：核心编译过；`add` 单测覆盖 F/R 父、恒等/非恒等父位置。

### 波 B：中央读取器（topo_tools_full/topexp/iterator/shell_check/shape_checks/wireframe）
- 已落地（2026-08-26）：`cumulated_children` 为 `TopoDS_Iterator` 默认 Compose/Move；`map_shapes`/`all_subshapes`/`Explorer`/`edges_of_wire`/`wires_of_face`/`edge_vertices`/`structure_is_valid`/`shell_check`/`wireframe` 面遍历共用。
- `map_shapes` 去重仍按 TShape 指针（对齐 `IsSame`，忽略 orientation）；seam 双边靠 `edges_of_wire`（不去重）。
- 验证：Cylinder 侧壁 `edges_of_wire` 返回 **4 边**（seam 双出现）；wire 内 #77 为 Reversed。

### 波 C：长尾 168 读取点（并行 agent 按文件分包，机械改）
- 验证：全 lib 编译 + 全量测试绿。

### 波 D：新行为消费
- `step.rs:945` 写真实 orientation（读 children 的 TopoShape.orientation）→ 往返保真。
- `make_pcurve_full` 接入 `pc_curve_orientation`（Reversed 边翻转 pcurve 方向）。
- 布尔 `loop_signed_area` / seam 处理验证（bop_build_faces）。
- 验证：Cylinder STEP 读→写→读 orientation 保真；布尔 seam 不塌缩。

### 波 E：回归门禁
- 新增测试：seam 双出现、orientation 往返、Reversed 边 pcurve 方向。
- 全量 `cargo test`（提交前一次，测试纪律）。

## 依赖与纪律
- 波 A→B→C 顺序依赖；D 依赖 C；E 最后。
- 每波：相关模块测试（`cargo test --lib <mod>::`）；全量仅提交前。
- 编译错误是迁移向导：改类型后 168 处逐一报错，机械修复。
- 匹配现有风格；不做无关重构。
