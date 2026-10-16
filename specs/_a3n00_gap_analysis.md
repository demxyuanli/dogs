# a3n00 网格/几何缺口分析（摘要）

> 完整逐轮记录（约 26000 行）已精简移除，原文保存在 git 提交 `c8060dd1`：
> `rtk git show c8060dd1:specs/_a3n00_gap_analysis.md`。
> 相关画板：`specs/occt-scope-drift.canvas.tsx`（旧 `specs/board.canvas.tsx` 已于 2026-10-16 退役删除），OCCT 参考数据：`specs/_occt_mesh_gt.md`。

---

## 1. 问题与起点（2026-09-26）

- 现象：`data/occ/a3n00.stp` 在 OCCT 中渲染正确；端口导出的 `a3n00.obj` 法兰与大曲率面破碎。
- 起点对拍（226 面，面数一致）：

| 量 | OCCT 8.0.0p1 | 端口（起点） |
|---|---|---|
| 顶点 / 三角 | 11052 / 12324 | 9091 / 8922 |
| 形状 bbox maxComp | 269.03 | 418.267 |
| Prs3d 偏转 lin | 1.07612 | 1.673068 |
| 无三角的面 | 0 | 4 |

- 表面类型分类（Plane 95 / Cylinder 62 / Cone 30 / Torus 7 / BSpline 32）两侧完全一致，问题不在曲面识别。

## 2. 根因链与结论

| 编号 | 结论 | 状态 |
|---|---|---|
| G1 | STEP 读入缺少 `StepToTopoDS_TranslateFace` 的 `ConvertToPeriodic`（`StepToTopoDS_TranslateFace.cxx:557-567`）：B 样条闭合面未设周期，pcurve 越域，`check_same_parameter` 报出 8.956，顶点容差被抬到 115.68，经 `BRepBndLib` 放大 bbox，lin 变成 1.673 | **已修**（§3） |
| G3 | 4 个空面 | 随 G1 消失 |
| G4 | 平面/柱面过密（1.16× / 1.27×） | 是 G1 造成的假象，修复后 1.01× / 0.98× |
| G2 | Cone/Torus/BSpline 三角数只有 OCCT 一半 | Cone、Torus 随 G1 恢复到约 1.07×；剩余的 7 个周期 B 样条面只有 1–11 个三角 |
| G2 余项 | 这 7 个面缺少 seam 边，`FixMissingSeam` 要用 `ShapeFix_ComposeShell` 插入缝边 | 由 T-92 完成（§4） |

- 已排除：STEP 曲面类型解析、`lin` 公式（`brep_bnd_lib`/`prs3d` 与 OCCT 控制流一致）、`surface_newton`、`value_of_uv`、`next_value_of_uv`、采样循环。

## 3. G1 修复（round 8）

- `crates/occt-geom/src/bspline_surface.rs`：新增 `set_u_periodic` / `set_v_periodic`（`Geom_BSplineSurface_1.cxx:940-1022`），`d0/d1/d2` 在周期时走 `prepare_eval::dn`。
- 新增 `crates/occt-topo/src/shape_custom_surface.rs::convert_to_periodic`（`ShapeCustom_Surface.cxx:475-616`）。
- `step/read_topology/translate_face.rs::resolve_face` 对 `B_SPLINE_SURFACE[_WITH_KNOTS]` 调用。

| 量 | 修前 | 修后 | OCCT |
|---|---|---|---|
| bbox gap | 74.633 | 0.000815 | — |
| maxComp | 418.267 | 269.00 | 269.03 |
| lin | 1.673068 | 1.076007 | 1.07612 |
| 有三角的面 | 222/226 | 226/226 | 226/226 |

## 4. T-92：ShapeFix_ComposeShell 移植（round 1–194，已收束）

- 新模块 `crates/occt-topo/src/shape_fix_compose_shell/`，对应 `ShapeFix_ComposeShell.cxx`（3603 行）。
- 覆盖：`SplitWire`、`BreakWires`、`CollectWires`、`SplitByLine`/`SplitByGrid`、`MakeFacesOnPatch`、`DispatchWires`、`FixMissingSeam` 尾段、`FixDegenerated` 回写等。
- `FixMissingSeam` 接入 STEP 读入（round 96），之后逐段修正拓扑：`WireSegment` 端点朝向（§9.157）、`MapReShape` 归一化（§9.103）、`CollectWires` 候选搜索（§9.493–9.532）等。
- §9.194–9.195 宣布 T-92 目标完成。

## 5. 读入链的后续移植（round 195 之后）

| 轮次 | 内容 | 结果 |
|---|---|---|
| T-101（F113 螺母斜切面） | `MapReShape::apply_impl`（§9.579–9.580） | a3n00 面积比 0.9963 → **1.0000**，WIREHIST 与 OCCT 逐项相同 |
| T-102 / T-103 | 环面碎片、共享 ReShape、分类器 XOR 反向（§9.585–9.599） | 用于 T0M、motoc 等模型 |
| T-104 | `Geom2d` 曲线 Segment、圆锥转 BSpline、SameRange 各臂（§9.604–9.612） | motoc 圆柱横向环带消失 |
| T-106 | `gp_Trsf`、`gp_Mat`、`gp_GTrsf` 等变换类正本清源（§9.614–9.623） | 门禁不劣化 |
| §9.624 | `ShapeFix_Shape` 共享 `ReShape` 改为唯一路径，去掉 `OCCT_SHARED_HEAL` 开关 | T0M 开放边 3643 → 250，a3n00 571 → 0 |
| §9.625 | `ShapeFix_Face::FixOrientation` 单环分支 | top.step 20042 → 22794（OCCT 22800） |
| §9.626 | `ShapeFix_ComposeShell::MakeFacesOnPatch` 全段 | 多环路径真实触发，门禁不变 |
| §9.627–9.631 | `IntCurve_IntConicConic` 专用重载、`IntPolyPolyGen`、`Extrema_ExtElC2d`、自交检查 | 门禁不变，部分为潜伏分支 |
| §9.632 | `wire_fix.rs` 拆为 `wire_fix/` 子模块（3379 行 → 多个文件） | 逻辑不变 |
| §9.633 | `FixIntersectingWires` 与 `FixWiresTwoCoincEdges` 接入 `ShapeFix_Face::Perform` | 结构对齐 |
| §9.634 | 提交前门禁与提交内容整理 | `step_obj_gates` 4/4 |

## 6. 当前基线（面积比，门禁口径）

| 模型 | 面积比 | 备注 |
|---|---|---|
| top | 1.00000 | |
| motoc | 1.00000 | |
| a3n00 | 1.0000 | §9.580 |
| ATU01038 | 0.99999 | |
| acs10 | 1.00004 | |
| T0M | 0.99998 | §9.625 实测 |
| TDB | 0.99978 | |
| **bottom** | **≈0.976** | 未闭环，见 §7 |

- 门禁：`occt-topo --lib` 1281/1281，`step_obj_gates` 4/4（23 个模型表），`fuse_box_cylinder_is_closed_solid` 通过。

## 7. 遗留项

1. **`data/occ/bottom.step`（面积比约 0.976，缺约 2.4%）**：
   - 包围盒与 OCCT 完全一致；323 个面可按包围盒一一配对，三角数合计只差 53。
   - 偏差主要集中在 `face 316`（竖直圆柱带孔，三通管区域）：端口 184 个顶点只用到 110 个，输出 104 个三角，OCCT 为 184 个。
   - 下一步方向：对照 `BRepMesh_Delaun` 与 `FixMissingSeam` 的孔洞处理，逐段核对 face 316 的环与边界离散。
2. 周期面 `ShapeFix_ComposeShell` 的个别分支仍有未对拍项（见 §4 的 round 194 之后记录）。
3. 历史曾关闭的 `FixMissingSeam` 接线已在 round 96 打开；如要回退，需同时保留 `bottom.step` 门禁。

## 8. 复现入口

```powershell
cd D:\source\repos\dogs
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj -- data/occ
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_probe_a3n00 -- data/occ/a3n00.stp --fstats
specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --facestats 0.451774 0.349066
```
