# Phase 21 Plan — STEP→OBJ 管线对拍移植（对齐 OCCT）

> 背景：本仓库翻译自 OCCT，但 STEP 读入→网格化→OBJ 导出走的路径与 OCCT 不一致，对拍门禁（`tests/step_obj_parity.rs`，以 `data/occ-*.obj` 为 oracle）量化出以下缺口。目标：按 OCCT 语义逐项移植，收紧门禁容差至对拍一致。
> 门禁现状：7 形状 bbox 通过（Cube/Cone 精确、Sphere/Torus 采样差、Cylinder/rev 宽容差、Offset 无参考）。密度差距记录未断言。
> 参考源：`D:\source\occt-src\src\DataExchange\TKDESTEP\StepToTopoDS\` + `ModelingAlgorithms\TKMesh\BRepMesh\` + `DataExchange\TKDEOBJ\RWObj\`。
> 基线：5 crate · 158,239 行 · 2,045 测试（topo 1333）。

## 缺口清单（对拍门禁量化）

| # | 缺口 | 现象 | OCC 机制 | 波及 | 影响 |
|---|---|---|---|---|---|
| 1 | face UV 域推导 | Cylinder 侧壁 v[-1,0] 应为 [0,10]；rev base 偏 | `BRep_Tool::UVBounds`/`BRepAdaptor_Surface` 用**边界边 pcurve** 定界，非 surface 无限 range | Cylinder/rev/Torus/Offset | 收紧 Cylinder 容差 1.05→0.05 |
| 2 | SURFACE_OF_REVOLUTION 未解析 | Shape/Shape-1/Shape-2 读入 0 形状；rev base 偏移 | `StepToTopoDS` 的旋转面构造（生成线绕轴） | Shape/Shape-1/Shape-2/rev | 解锁 4 文件完整读入 |
| 3 | ORIENTED_EDGE 方向忽略 | `resolve_oriented_edge` 丢 `.T./.F.` | ORIENTED_EDGE 的 orientation 影响边在 wire 的方向 | seam 面/所有面 | UV 域正确性的前置 |
| 4 | 网格化密度 | Cube 30k 面 vs OCCT 12 | `BRepMesh_IncrementalMesh` 平面精确+曲面 deflection | 全部 | 密度对齐（记录→断言） |
| 5 | pcurve 存储缺失 | 边无 pcurve → UV 域无法从 pcurve 推导 | `BRep_Tool::CurveOnSurface` 每边一条 pcurve | 缺口 1 的前置 | EdgeGeom.pcurves 已加字段 |

## 波次方案（每波独立可验证，用门禁收紧验证）

### 波 A：pcurve 基建 + face UV 域（前置，~2,500 行）
- **A1 `pcurve 域推导`**：`BRepTool::uv_bounds` 对无限 surface range 时，用每条边在面上的 pcurve 包围盒定界（OCCT `BRepAdaptor_Surface::UVBounds` 语义）。需边有 pcurve——复用 Phase 16 的 `pcurve_full::make_pcurve_full` + `EdgeGeom.pcurves`。
- **A2 门禁收紧**：Cylinder 容差 1.05→0.05，rev 4.3→0.05。
- 模块：改 `brep_tool.rs`/`brep_surface.rs`，新 `pcurve_bounds.rs`。

### 波 B：SURFACE_OF_REVOLUTION + 曲面类型补全（~3,000 行）
- **B1 `SURFACE_OF_REVOLUTION`**：生成线（B样条）绕轴旋转 → `GeomSurfaceOfRevolution`（StepToTopoDS 语义）。读入 Shape/Shape-1/Shape-2/rev。
- **B2 其它缺失曲面**：对拍揭示的 SURFACE_* 变体逐一补 `resolve_surface` 分支。
- 模块：`occt-geom` 加 `surface_of_revolution.rs`，`step.rs` 加解析分支。
- 门禁：Shape/Shape-1/Shape-2 从"0 形状"→读入成功 + bbox 对拍。

### 波 C：ORIENTED_EDGE 方向 + 网格化对齐（~2,500 行）
- **C1 `ORIENTED_EDGE` 方向**：`resolve_oriented_edge` 保留 `.T./.F.`，影响 face wire 边方向（seam 面）。
- **C2 网格化密度**：`mesh_shape` 接入 `brepmesh::incremental_mesh`（平面精确+曲面自适应），修 `point_triangle_distance` 交互（盒内点距离 0.5 回归）。
- 模块：改 `step.rs`/`shape_mesh.rs`/`geometry_query.rs`。
- 门禁：Cube f-ratio 从 ~4000 降到 ~1；密度记录→断言。

## 验证（每波）

- **波 A**：Cylinder 侧壁 v[0,10]、bbox z[0,10] 与 occ 一致（容差 0.05）；rev base 修正。
- **波 B**：Shape.step 读入 ≥1 形状、bbox 对拍 occ-shape.obj；rev 读入完整。
- **波 C**：Cube f-ratio ≈1（12 面 vs 12 面）；`point_distance_inside_and_outside` 恢复 d_in=0.5。
- 全量回归：提交前跑一次（测试纪律）。

## 依赖与纪律

- 每波独立模块/改动，可并行 agent（A/B/C 波内拆 agent）。
- 匹配现有风格（`Result<T,String>`、`mod tests`）；pcurve 复用 Phase 16。
- 对拍门禁是每波的验收标准（容差收紧 = 移植完成）。
- 不碰 bop_builder/fillet 等布尔模块（无涉）。
