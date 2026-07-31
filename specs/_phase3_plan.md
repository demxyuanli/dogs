# Phase 3 Plan — 从"几何引擎"到"可操作 BRep 模型"

## 现状 (2026-07-31)
- 5 crates · 179 模块 · 18,959 行 · 215 tests
- ✅ 几何求值 / 网格交换 / 数学内核 / 基本体创建 / 测量
- ⚠️ **TShape 无几何数据** → 任意形状不可网格化/不可查询/不可导入 STEP
- ⚠️ 无曲线曲面逼近 → B样条无法落到网格

## 核心洞察
下一阶段瓶颈 = **TShape 几何数据存储**。所有"真正做 CAD"的能力
(可视化/导入导出/分析) 都依赖它。优先打通这条主干。

---

## Phase 3A — BRep 几何数据内核 (最高优先, ~3000 行)

| # | 模块 | 内容 | 预估 | 依赖 |
|---|---|---|---|---|
| A1 | `tgeometry` | TShape 挂接几何：VertexShape 存 GpPnt+tol、EdgeShape 存曲线句柄 Arc<dyn Curve>+参数区间、FaceShape 存 Arc<dyn Surface>+UV 边界 | 400 | geom |
| A2 | `brep_builder` | 真实顶点/边/面构建：MakeVertex(point)、MakeEdge(curve,a,b)、MakeFace(surface, wires)，把几何写入 TShape | 500 | A1 |
| A3 | `brep_tool_full` | 完成 BRep_Tool：Pnt(v)、Curve(e)、Surface(f)、Tolerance、ParameterRange、UVBounds | 300 | A1 |
| A4 | `topexp_full` | 真实子形状存储：TShape 持有 children Vec<HandleTShape>，Explorer 走真树 | 400 | A1 |
| A5 | `wireframe` | 边→折线逼近（B样条→Bezier→采样），面→三角化网格 | 500 | A4 + gcpnts |
| A6 | `shape_mesh` | 任意 BRep 形状网格化：遍历边/面 → 组装 Triangulation + 法线 + 缝合 | 600 | A5 |
| A7 | `bbox_from_geometry` | 从几何计算真实包围盒（曲线采样 + 曲面采样） | 300 | A5 |

**里程碑 A**: `make_box` 产生**真几何** box，`brep_tool` 能取回顶点坐标和面曲面，`shape_mesh` 能把任意 BRep 三角化 → **可视化/导出任意形状**。

## Phase 3B — 曲线曲面转换 (支撑 A5, ~1500 行)

| # | 模块 | 内容 | 预估 |
|---|---|---|---|
| B1 | `bspline_to_bezier` | 完整 B样条→Bezier 段提取（含有理、周期性） | 400 |
| B2 | `bezier_to_polyline` | Bezier → 折线 (de Casteljau 细分 + 弦偏差) | 300 |
| B3 | `curve_approx` | 参数曲线 → 折线（GCPnts 适配 Curve trait） | 300 |
| B4 | `surface_to_grid` | 曲面 → UV 网格 (等参线采样) | 300 |
| B5 | `convert_geom` | gp↔Geom、Geom↔Geom2d 投影转换 | 200 |

## Phase 3C — 基础布尔与拓扑分析 (~1500 行)

| # | 模块 | 内容 | 预估 |
|---|---|---|---|
| C1 | `bbox_intersect` | 包围盒级剔除（加速 A6/A7 已有） | — |
| C2 | `face_face` | 面-面求交（平面-平面 / 平面-曲面 采样近似） | 500 |
| C3 | `edge_split` | 边在交点分裂 → 拓扑重建 | 400 |
| C4 | `shell_check` | 闭合壳检查（每条边被两面共享） | 300 |
| C5 | `solid_union` | 基于 A6 网格的体素/网格布尔融合 | 300 |

## Phase 3D — 模型导入通道 (~1000 行)

| # | 模块 | 内容 | 预估 |
|---|---|---|---|
| D1 | `mesh_to_brep` | Triangulation → BRep (顶点/边/面从网格重建) | 500 |
| D2 | `brep_exchange` | BRep → OBJ/STL/PLY 导出（走 A6 网格化） | 300 |
| D3 | `brep_scene` | BRepModel ↔ 场景图，多形状组合导出 | 200 |

## 里程碑
- **M1 (A1-A4)**: BRep 可存几何、可查询 → 为 STEP 导入打基础
- **M2 (A5-A7 + B)**: 任意形状可网格化/可视化/导出 → **首个"真模型"**
- **M3 (C+D)**: 网格导入 + 基础布尔 → 可编辑模型闭环

## 顺序建议
```
A1→A2→A3→A4  (几何内核, 必须先)
   ↓
B1→B3→B5     (曲线逼近, A5 依赖)
   ↓
A5→A6→A7     (网格化, M2 完成)
   ↓
C4→C2→C3     (拓扑分析)
   ↓
D1→D2→D3     (导入通道)
```

## 并行建议
- A1-A4 串行（强依赖）
- B 系与 A2-A4 可并行（B 只依赖 geom+bspl）
- C5/D 系等 A6 完成后并行

## 风险
- A1 需要把 `Arc<dyn Curve>` 塞进 TShape — 打破当前 TShape 的 Clone 要求（trait object 不可 Clone）→ 用 Arc 持有，TShape 本身保留 Clone
- 曲线/曲面 trait 对象在 TShape 内需要 `Debug` — 当前 Curve 无 Debug → 用自定义 Display 或去掉 TShape 的 Debug
