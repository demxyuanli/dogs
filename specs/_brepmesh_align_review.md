# BRepMesh 对齐审查文档 — seam pcurve 双面 + B 样条边参数域

> ⚠️ **仅历史（T-29 登记）**：本文件记录的是提交 `66615d7` 当时的审查结论；此后网格/BOP 管线已有多次变更，
> 现行基线与门禁以 `specs/_board.md` §2 / §3 为准，本文件不参与任何门禁。
>
> 本次提交 `66615d7` 的详细审查。目标：把 STEP→OBJ 网格管线的三处 bbox 漂移
> （Shape-2 / Cone / Torus）逐条追到 OCCT 源码对应，修复并验证。
> 参考源：`D:\source\occt-src\src\`（TKDESTEP / TKShHealing / TKMesh）。

---

## 1. 结论摘要

| 症状 | 根因（一句话） | OCCT 对应 | 修复 |
|---|---|---|---|
| Shape-2 bbox ±135→±60 | 边参数域启发式把有理 B 样条圆弧误判 parabola | `StepToTopoDS_TranslateEdge::MakeFromCurve3D` 投影顶点取参数 | 有界非周期曲线用自然 knot 域 |
| Cone bbox min z 9.28→0 | seam 边两侧 pcurve 被单槽覆盖 + 去重丢掉一侧朝向 | `B.UpdateEdge(E,C2d1,C2d2,Face)` 存两条 + `SelectForwardSeam` | pcurve 存 Vec + 按朝向选边 |
| Torus bbox min x -8→-12 | 闭合圆 minor seam 顶点重合，pcurve 参数域塌缩成单点 | SameParameter 边 pcurve 与 3D 同参数域 | 有界 3D 曲线用 `[a,b]` 作 pcurve 域 |

验证：parity 门禁 **12/12**，occt-topo lib **1291/1291**，无回归。

---

## 2. 逐条根因链

### 2.1 Shape-2（±135 → ±60）

**症状**：Shape-2.step 是 B 样条曲面模型（有理曲面 + knots）。顶缘边界圆是
`B_SPLINE_CURVE(6, 7 poles, knots(0., 0.7633))` 的有理圆弧。bbox 漂移到 ±128.15
（parity 注释记为 ±135），OCCT 是 ±60。

**根因链**：
1. `edge_params_for_curve`（`step.rs`）用 `classify_curve` 几何启发式判曲线类型：
   采样 6 个 `d2`（曲率），`|d²|` 近恒定（`(max-min)/max < 0.02`）且 `is_periodic()==false`
   → 判成 **Parabola**（`step.rs:205-206`）。
2. Parabola 分支做点积，把参数域算成 `[0, 28.82]`，而曲线真实 knot 域是 `[0, 0.7633]`。
3. 网格在 t=28.82（超域 ~37.8 倍）处求值曲线 → 多项式外推鼓包 → bbox ±128。

**OCCT 对应**：`StepToTopoDS_TranslateEdge::MakeFromCurve3D`
（`StepToTopoDS_TranslateEdge.cxx:443-485`）把每个顶点投影到曲线上
（`ShapeAnalysis_Curve::Project`）取参数 U1/U2，再 `BRepLib_MakeEdge(C1,V1,V2,U1,U2)`。
**不分类、不启发式**——有界曲线的自然参数域就是边域。

**修复**（`step.rs::edge_params_for_curve` 加 6 行 early-return）：
```rust
// 有界非周期曲线（B 样条/裁剪曲线）：自然域即边域
if f.is_finite() && l.is_finite() && !curve.is_periodic() {
    return if l > f { (f, l) } else { (0.0, 1.0) };
}
```
无界（line/parabola）、周期（circle/ellipse）仍走旧启发式，不碰 analytic 曲线。

**验证**：bbox `[-128.15,-128.15,0]~[128.15,128.15,150]` → `[-60,-60,0]~[60,60,150]`，
与 OCCT 逐字节一致。parity 门禁 shape2 从 `100.0` 收紧到 `0.05`。

---

### 2.2 Cone（min z 9.28 → 0）

**症状**：`CONICAL_SURFACE`（顶点 z=0、底面 z=10 半径 4）。网格 bbox min z=9.28，
顶点 z=0 整段丢失，只剩顶上一圈 v∈[10,10.77]。

**根因链**（三层）：
1. **pcurve 存储覆盖**：seam 边 #55 有两条 pcurve——#63（u=2π）和 #70（u=0），
   是锥切缝的两侧。`EdgeGeom.pcurves` 是 `HashMap<face, 单条>`，`set_pcurve`
   用 `insert` 覆盖，第二条冲掉第一条，只剩 u=0。
2. **wire 去重丢朝向**：`orient_edges_on_wire`（`algo_tools.rs:398`）按 **TShape**
   去重，把 seam 边的两个朝向（Reversed + Forward）合并成一个。实测
   `stored=2 orient=["Reversed","Reversed"]`——Forward 侧被丢弃。
3. 结果：网格模型里 seam 只有**一个** pcurve slot（Reversed），UV 多边形缺了
   Forward 侧（u=2π），Delaunay 只铺出顶上一条。

**OCCT 对应**：
- `StepToTopoDS_GeometricTool::IsSeamCurve`（`StepToTopoDS_GeometricTool.cxx:77-121`）：
  同一 wire 里两个 ORIENTED_EDGE 共享同一边 → seam。
- `StepToTopoDS_TranslateEdgeLoop`（`StepToTopoDS_TranslateEdgeLoop.cxx:689-766`）：
  `B.UpdateEdge(E, C2d1, C2d2, Face, 0.)` 存**两条** pcurve（forward 在前）。
- `ShapeAnalysis_Curve::SelectForwardSeam`（`ShapeAnalysis_Curve.cxx:852+`）：
  对 +v 方向的 seam，取 u 坐标较大者为 forward。

**修复**（四步联动）：
1. `tgeometry.rs`：`EdgeGeom.pcurves` 改 `HashMap<face, Vec<Curve2d>>`，加
   `set_pcurves` / `edge_pcurves`。
2. `pcurve_full.rs`：移植 `select_forward_seam`（返回 1/2），`make_pcurve_full`
   按 `edge.orientation()` 选 forward/reversed 侧。
3. `step.rs::associate_edge_pcurve`：收集两条 pcurve，`select_forward_seam` 定序，
   `set_edge_pcurves` 整组替换。
4. `algo_tools.rs::orient_edges_on_wire`：去重条件改「同 TShape 且同朝向」。
5. `model_builder.rs::add_wire` + `incremental_mesh.rs::discretize_pcurves`：
   传边**实际朝向**而非强制 Forward。

**验证**：bbox min z 9.28 → 0，Cone 恢复完整锥体。

---

### 2.3 Torus（min x -8 → -12）

**症状**：`TOROIDAL_SURFACE`（大半径 10、小半径 2）。网格 bbox min x=-8，外圈
（半径 8→12）整圈丢失。

**根因链**：
1. minor seam（#52 圆半径 2，圆心 (10,0,0)）是**闭合圆**，起点=终点=内点 (8,0,0)。
2. `discretize_pcurves` 的 pcurve 参数域 `(pf,pl)` 对无界 pcurve（LINE）用**顶点投影**：
   两个顶点都投影到 v=π（内点）→ `(pf,pl)=(π,π)` 塌缩成单点。
3. 所有 pcurve 点都得到 v=π，外圈（v=0/2π）整段丢失。

**OCCT 对应**：SameParameter 边（闭合圆）pcurve 与 3D 曲线同参数域——3D 圆
`[0,2π]` 就是 pcurve 的 `[0,2π]`。

**修复**（`incremental_mesh.rs::discretize_pcurves`）：
```rust
let (pf, pl) = if pc.first_parameter().is_finite() && ... {
    (pc.first_parameter(), pc.last_parameter())
} else if a.is_finite() && b.is_finite() && b > a {
    (a, b)   // 有界 3D 曲线：pcurve 同参数域，顶点投影会塌缩
} else {
    // 无界 3D 曲线（line）：顶点投影取弧长
    ...
};
```

**验证**：bbox min x -8 → -12，Torus 恢复完整外圈。

---

## 3. OCCT 源码对应速查

| 功能 | OCCT 位置 | Rust 目标 |
|---|---|---|
| B 样条曲线权重/knots/极点 | `StepToGeom.cxx:751-928` `MakeBSplineCurveCommon` | `step.rs::resolve_curve` |
| 边参数域（顶点投影） | `StepToTopoDS_TranslateEdge.cxx:443-485` | `step.rs::edge_params_for_curve` |
| seam 检测 | `StepToTopoDS_GeometricTool.cxx:77-121` `IsSeamCurve` | `step.rs` |
| seam 双 pcurve 赋值 | `StepToTopoDS_TranslateEdgeLoop.cxx:689-766` | `step.rs::associate_edge_pcurve` |
| forward pcurve 选择 | `ShapeAnalysis_Curve.cxx:852+` `SelectForwardSeam` | `pcurve_full.rs::select_forward_seam` |
| pcurve 存储（两条） | `BRep_Builder::UpdateEdge(E,c1,c2,face)` | `tgeometry.rs::set_edge_pcurves` |

---

## 4. 架构改动一览

| 文件 | 改动 | 类型 |
|---|---|---|
| `tgeometry.rs` | pcurve 存储 `Vec` + `set_pcurves`/`edge_pcurves`/`get_pcurves` | 存储 |
| `pcurve_full.rs` | `select_forward_seam` 移植 + `make_pcurve_full` 按朝向选边 | 解析 |
| `step.rs` | `edge_params_for_curve` early-return + `associate_edge_pcurve` 排序 | 解析 |
| `algo_tools.rs` | `orient_edges_on_wire` 去重改「TShape+朝向」 | 拓扑 |
| `meshing/model_builder.rs` | `add_wire` 传实际朝向 | 网格 |
| `meshing/incremental_mesh.rs` | `discretize_pcurves` 按朝向 + 闭合圆 `[a,b]` 域 | 网格 |

---

## 5. 验证

- **parity 门禁** `cargo test --test step_obj_parity`：**12/12 通过**（此前 cone/torus 失败）。
- **occt-topo lib** `cargo test --lib`：**1291/1291 通过**。
- **Shape-2 门禁收紧**：`100.0 → 0.05`，锁入修复。

---

## 6. 残留差距（来自 `_brepmesh_align_tasks.md`，非本次修复）

1. **四叉树兜底**（`export_mesh` 的 `brepmesh` 路径）仍用旧的点向投影，未对齐 Delaunay。
2. **NURBS RangeSplitter 的移除/thinning pass 未做**（纯密度优化，不影响 bbox）。
3. **ModelHealer 的 `amplifyEdges`/`FaceChecker`** — 已接到 `IncrementalMesh::heal_self_intersecting_wires`（`BRepMesh_ModelHealer.cxx:234-211`）。
4. **FaceDiscret 的 uniform 内部 UV 网格未加**（`Triangulator` 的 deflection 细化替代）。

---

## 7. 岔路记录（供追溯）

排查中验证过、但**不是**根因的三条岔路，均已还原：

- **曲面 d0 域外 clamp**（`occt-geom/src/bspline_surface.rs::basis_funs`）：
  clamp 不在网格路径上，加/去 clamp bbox 逐字节不变。
- **有理曲线权重丢失**：顶缘曲线权重 ≈1.007 近等权，丢权重只差 ~1%，不会鼓 2.3×。
- **曲面周期 wrap**：曲面 #40 是 `closed_u=FALSE` 非周期，OCCT 也不 wrap。
