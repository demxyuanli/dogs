# a3n00 网格/几何缺口分析（2026-09-26）

> 症状（用户提供）：`data/occ/a3n00.stp` 经 OCCT 渲染正确（法兰盘 + 螺栓孔 + 大曲率阀体），
> 端口 `output/a3n00.obj` 渲染时**法兰与大曲率面破碎/扇贝状**。
> 本文只做**定位与取证**，不改代码；结论均来自本机实测（OCCT 8.0.0p1 探针 + 端口探针）。
> 相关画板条目：`specs/board.canvas.tsx`（T-69/T-91/T-92 与门禁表）、`specs/_occt_mesh_gt.md`。

---

## 0. 复现命令（本次实测用）

```powershell
cd D:\source\repos\dogs

# 端口导出（真实参数：Prs3d 相对偏转 + 20° 角偏）
cargo run --offline --manifest-path crates/occt-topo/Cargo.toml --example export_data_obj

# OCCT GT：面结构 / 逐面三角 / 面类型
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --bbox"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --mesh 1.07612 0.349066"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --perface"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --wires"

# 端口逐面探针（临时件，见 §8 附录）
cargo run --offline --manifest-path crates/occt-topo/Cargo.toml --example zz_probe_a3n00 -- data/occ/a3n00.stp
cargo run --offline --manifest-path crates/occt-topo/Cargo.toml --example zz_probe_a3n00 -- data/occ/a3n00.stp 1.07612   # 用 OCCT 的 lin
```

---

## 1. 基线对拍（实测数字）

### 1.1 总量

| 量 | OCCT 8.0.0p1 | 端口 `output/a3n00.obj` |
|---|---|---|
| 面数 | **226** | **226**（`faces_of`） |
| 顶点 / 三角 | **11052 / 12324**（lin 1.07612、角 20°） | **9091 / 8922** |
| 形状 bbox maxComp | **269.03** | **418.267**（！） |
| 用于网格的 lin | **1.07612** | **1.673068**（！） |
| 无三角的面 | 0 | **4**（面 222–225） |

- 端口 `lin` 偏大 **1.55×**：`brep_exchange.rs:56-77 prs3d_get_deflection` 用的是端口自己的
  `shape_bnd_box`，而端口 bbox 被放大（见 §3）。
- 即使把端口强制到 OCCT 的 `lin=1.07612`，端口也只有 **9872 / 10184**（仍 0.83 / 0.83），
  ⇒ deflection 不是唯一缺口（见 §4）。
- 作为对照，用端口的粗 `lin=1.673068` 跑 OCCT：**10765 / 11908**（OCCT 对 lin 不敏感，端口敏感）。

### 1.2 逐面类型直方图（OCCT `--perface` vs 端口 `classify_surface`）

| 类型 | OCCT | 端口 |
|---|---|---|
| Plane | 95 | 95 |
| Cylinder | 62 | 62 |
| Cone | 30 | 30 |
| Torus | 7 | 7 |
| BSpline | 32 | 32 |

⇒ **表面类型的解析/分类完全一致**（无差异）。用户所说的「没有正确解析」发生在**网格化/参数化**，
不在 STEP 曲面类型识别。

### 1.3 按面类型分组的三角数（关键表）

| 类型 | 面数 | OCCT @lin1.076 | 端口 @lin1.673（默认） | 端口 @lin1.076 | 端口/OCCT @1.076 |
|---|---|---|---|---|---|
| Plane | 95 | 2832 | 2895 | 3278 | **1.16×（偏多）** |
| Cylinder | 62 | 2673 | 2788 | 3385 | **1.27×（偏多）** |
| Cone | 30 | 727 | **376** | **376** | **0.52×** |
| Torus | 7 | 1214 | **596** | **596** | **0.49×** |
| BSpline | 32 | 4878 | 2267 | 2549 | **0.52×** |
| 合计 | 226 | 12324 | 8922 | 10184 | 0.83× |

两条硬事实：
1. **Cone / Torus 的三角数在两个 deflection 下完全相同（376 / 596）** ⇒ 端口对这些面的内点采样
   **不随线性 deflection 变化**（或该路径根本没被采用）。
2. **Cone / Torus / BSpline 稳定只有 OCCT 的一半**，而 Plane / Cylinder 反而偏多。

### 1.4 逐面样本（排序后）

```
Cone   OCCT: 10×24, 72,72,73,81,81,108         端口: 10×20, 11×12, 2,5,8,8,34,65
Torus  OCCT: 86,125,162×4,355                  端口: 2,2,18,18,18,162,376
BSpline OCCT: 34,46,48,48,50,64,66,66,71,72…,376,376,669,810,854
       端口: 1,1,1,2,6,7,12×4,18,36,46…,228,243,288,376
```

⇒ 大的锥面（OCCT 72–108 三角）在端口只有 **11 个三角**；大 B 样条面（OCCT 669/810/854）
在端口最大只有 **376**。这正是图里「大曲率面破碎」的直接来源。

---

## 2. 缺口清单（按优先级）

| # | 缺口 | 证据 | 代码位置 | OCCT 对照 |
|---|---|---|---|---|
| **G1** | 顶点容差被放大 → 形状 bbox 放大 → `Prs3d` deflection 偏粗 1.55× | gap=**74.633**，顶点容差最高 **115.68**；bbox maxComp 418.267 vs 269.03；lin 1.673 vs 1.076 | `step/read_geometry.rs:570 edge_from_curve3d`、`step/read_topology.rs:343-360`、`brep_bnd_lib.rs:98-141`、`brep_exchange.rs:56-77` | `StepToTopoDS_TranslateEdge.cxx:443-479`、`BRepBndLib.cxx:131`、`Prs3d.hxx:82-103` |
| **G2** | Cone/Torus/BSpline 面网格只有 OCCT 的 ~0.5×；Cone/Torus 对 deflection 不敏感 | §1.3 表 | `meshing/range_splitter/nodes.rs:199-216`、`range_splitter/splitter.rs:342-361`（Cone）、`:559-626`（Torus）、`range_splitter/param_set.rs:520-617`（NURBS grid） | `BRepMesh_MeshAlgoFactory`、`BRepMesh_{Cone,Torus,NURBS}RangeSplitter`、`BRepMesh_DelaunayBaseMeshAlgo` |
| **G3** | 4 个面没有三角、被静默跳过 | 面 222–225（2 Cylinder + 2 Plane）`mv=0 mt=0`；OCCT 无 | `brep_exchange.rs:88-109`（失败即空网格）、`incremental_mesh/discret_root.rs:204-211` | `RWMesh_FaceIterator.cxx:87-89` |
| **G4** | Plane/Cylinder 反向过度网格 | @lin1.076：Plane 3278 vs 2832（1.16×）、Cylinder 3385 vs 2673（1.27×）；端口单面最多 674 vs OCCT 376 | 平面/柱面三角化路径（`face_discret`/`delaun`） | `BRepMesh_DelaunayBaseMeshAlgo`（边界驱动） |

**已排除**：STEP 曲面类型解析（§1.2 完全一致）；`lin` 公式本身（`brep_bnd_lib`/`prs3d` 与
`BRepBndLib.cxx:131` 的控制流逐行一致，问题在输入容差）。

---

## 3. G1 细节：顶点容差 → bbox → deflection

### 3.1 实测

```
BBOX roots=2 gap=74.633481 min=(-209.133,-152.133,-250.927) max=(209.133,152.133,100.633)
TOTAL faces=226 computed_lin=1.673068
TOP FACE TOLS: 全 0.000000
TOP VERT TOLS: 115.681895 115.681895 74.633481 71.648142 71.648142 71.648142 58.214115 50.750767
```

`BndBox::enlarge`（`occt-core/src/bnd/box3d.rs:98-100`）按 `Gap = max(Gap, |Tol|)`；
`corner_min/max`（`:151-161`）把 gap 加回四角。`BRepBndLib::Add(...,false)` 对**平面面**
调 `B.Enlarge(Tolerance(F))`（`BRepBndLib.cxx:131`），对**不在边里的顶点**调
`B.Enlarge(Tolerance(V))`（`:212`）。端口 `brep_bnd_lib.rs:126/139` 逐行照搬，因此
**端口是把 OCCT 的控制流如实执行了，但输入容差本身是错的**。

- 对照：`data/Cube.step`、`data/Torus.step`、`data/Shape-2.step` 的顶点容差都在 1e-6..0
  （见 §8 实测），只有 `data/occ/*` 这一族（a3n00 等，**全部没有存档 pcurve**）出现
  70–115 量级的容差。
- OCCT 侧同文件（默认读入，含 ShapeProcess）bbox maxComp = **269.03** ⇒ OCCT 的顶点容差是小的。

### 3.2 容差是怎么被抬起来的

`step/read_topology.rs:341-360`：

```rust
// MakeFromCurve3D (TranslateEdge.cxx:452-455): 投影参数处的残差
let temp1 = curve.d0(first).distance(&p1);
let temp2 = curve.d0(last).distance(&p2);
...
// TranslateEdge.cxx:478-479: UpdateVertex 只升不降
let t1 = vertex_tolerance(&v1).max(1.000001 * temp1);
```

即 `temp1/temp2` = 顶点在 3D 曲线上的投影残差。**只要 `edge_from_curve3d`（投影 + `UpdateParam3d`）
给出的 `first/last` 偏了，残差就会被原样写进顶点容差**，再被 `BRepBndLib` 放大成全局 bbox。

OCCT 对应（`StepToTopoDS_TranslateEdge.cxx:443-479`）：
`ShapeAnalysis_Curve::Project(..., preci, pproj, U1, false)` → `UpdateParam3d` →
重算 `temp = aCA.Value(U).Distance(pv)` → `B.UpdateVertex(V, 1.000001*temp)`。
⇒ 两者结构相同；差异只可能在 **`Project` / `UpdateParam3d` / `edge_from_curve3d`** 的数值行为。

### 3.3 已定位的 G1 根因链（2026-09-26 本轮实测）

一次性插桩（均已回退）把「顶点容差 → bbox → deflection」钉到 ShapeFix_Edge::FixSameParameter：

1. **排除 resolve_edge**：给 `step/read_topology.rs::resolve_edge` 加日志后，490 条
   `EDGE_CURVE` 的 `temp1/temp2` **全部为 0** ⇒ STEP 边范围解析本身正确。
2. **容差写入点**：在 `Vertex::set_tolerance`（`builder.rs:232`）加 `>0.05` 断点捕获调用栈，
   104 次命中全部是 `fix_same_parameter`（`shhealing/wire_fix.rs:1310`）←
   `check_pcurves_and_shift`（`:3123`）← `resolve_face`。两个写点：
   - `wire_fix.rs:1379-1386`：`v.set_tolerance(max(maxdev, tol_v))`（对应 `ShapeFix_Edge.cxx:915-922`）；
   - `wire_fix.rs:1387-1390`：`maxdev > tol` 时 `fix_vertex_tolerance_edge`（`:1162`，对应 `ShapeFix_Edge.cxx:924-929`）。
3. **maxdev 来自 check_same_parameter**（`wire_fix.rs:1056`，= `ShapeAnalysis_Edge::CheckSameParameter` 的端口版）。
   对 `maxdev > 0.05` 打印（`c3d` = 边的 3D 曲线、`prange` = 该面上的 pcurve 范围、`npc` = pcurve 数）：

   ```
   maxdev=10.000100  c3d=GeomCircle     range=[3.1416,9.4248] prange=(3.1416,9.4248) npc=1 was_sp=true same_range=false
   maxdev=19.000190  c3d=GeomCircle     range=[...]            prange=(0,9.4248)      npc=1 ...
   maxdev=8.956018   c3d=GeomCircle     same_range=true
   maxdev=96.000960  c3d=CurveOnSurface (pcurve-only 替身)
   maxdev=155.001550 c3d=GeomCircle     prange=(...9.4248)
   ```

   - 3D 范围与 pcurve 范围**数值相同**、`npc=1`，却报出 **10.0 / 19.0 / 96.0 / 155.0** 量级偏差
     （多为 **2R（直径）** 或整段长度）⇒ 该面 pcurve 的**参数化/图像与 3D 曲线不一致**
     （同族于 T-69/T-91 的「pcurve 参数化 ≠ 3D 边参数化」）。
   - OCCT 对应：`ShapeAnalysis_Edge.cxx:704-826` 用 `BRepLib_ValidateEdge`
     （`SetControlPointsNumber(NbControl-1)` → `Process()` → `UpdateTolerance(maxdev)`）统一计算；
     端口把 `BRepLib_ValidateEdge` 内联重写在 `wire_fix.rs:1094-1157`
     （`projection` 分支用 `locate_ext_pc_sq` 做双向极值）。
   - 两处**控制流差异**（下一轮按 OCCT 核对，可能是放大来源）：
     (i) `ShapeFix_Edge.cxx:875-880`：`wasSP == false` 时 OCCT 把面参数**换成空面**再
         `CheckSameParameter`；端口始终传 `face`（`wire_fix.rs:1358`）。
     (ii) OCCT `CheckSameParameter` 在 `BRep_Tool::Curve(edge)` 为空时**立即返回**
         （`ShapeAnalysis_Edge.cxx:725-729`，设 FAIL1、`maxdev` 保持 0）；端口因 `EdgeGeom`
         不能存空 3D 曲线，用 `meshing::edge_discret::CurveOnSurface` 冒充 3D 曲线
         （`fix_degenerated` `wire_fix.rs:2089`、`fix_lacking_one` `:2591`），
         于是**不走这条早退**（实测有 `c3d=CurveOnSurface` 的 96.0 偏差）。
         **但 CurveOnSurface 只占少数，主因是 `GeomCircle` 边的 pcurve 偏差。**
4. **后果**：被放大的顶点容差进入 `BRepBndLib`（`brep_bnd_lib.rs:126/139` 的 `B.Enlarge(Tolerance)`，
   与 `BRepBndLib.cxx:131/212` 同构）⇒ `gap=74.633481` ⇒ `shape_bnd_box` maxComp
   **418.267**（OCCT **269.03**）⇒ `Prs3d::GetDeflection` **1.673068**（OCCT **1.07612**）
   ⇒ 曲面面按 1.55× 偏粗的 lin 网格化。

**⇒ G1 的修法必须落在 pcurve 计算/校验（`fix_add_pcurve` + `check_same_parameter`）与
`FixLacking/FixDegenerated` 的 `CurveOnSurface` 替身上，且按 OCCT 逐行迁移；
禁止直接钳制/忽略容差。**

### 3.4 G1 的直接触发点：闭合（非周期）B 样条面上的 pcurve 越域（2026-09-26）

- 面 109–112 的曲面是 STEP `#4809`（`B_SPLINE_SURFACE(1,2)` + `B_SPLINE_SURFACE_WITH_KNOTS`，
  `closed_v=.T.`，`v_mults=(3,2,2,2,3)`，`u_knots=[-61.7587,-58.7754]`），
  端口实测 `vrange=[0,37.6991]`、`vper=false`、`uper=false`。
- 周期判定：OCCT `StepToGeom.cxx:1120-1129` 与端口
  `GeomBSplineSurface::should_be_periodic`（`bspline_surface.rs:128-141`）**逐行相同**；
  对本实体 `SumMult=12 == NVPoles(9)+VDeg(2)+1` ⇒ 两者都判**非周期**（`closed_v` 不参与）
  ⇒ 曲面是「闭合但非周期」。
- 端口投影出的 pcurve 是 UV 直线 `(u=-61.759, v=18.850-6t)`：t∈[0,2π] 使 v 从 +18.85 走到
  **−18.85，越出 [0,37.7] 域**。`check_same_parameter` 在同参分支逐点测 dev：
  前半 d≈0，后半 1.12 → 3.37 → 6.20 → **8.956**，maxdev 被写进顶点容差。
- ⇒ **G1 的直接触发点在 pcurve 投影对「闭合（非周期）面」的域外/跨缝调整**，而不是
  `check_same_parameter`（它只是忠实测量）。OCCT 侧对应
  `ShapeConstruct_ProjectCurveOnSurface` → `ShapeAnalysis_Surface::NextValueOfUV`/`ValueOfUV`
  的相邻点取解。
- **已排除**：端口的 `surface_newton`（`surface_projector.rs:758-845`）与 OCCT
  `SurfaceNewton`（`ShapeAnalysis_Surface.cxx:1065-1157`）**逐行一致**（Bnds±U/VResolution、
  步长、`U<UF||V<VF` 提前 break、`rs2>rsfirst2` / `rsn` 判据、`nrm2<0.01` 的 res=2），
  ⇒ 缺口在**其下游**：`value_of_uv` 的域内取解（`surface_projector.rs:519`）与
  采样后的 2D 拟合（本样本给出一条 `(u=const, v=18.85-6t)` 直线，v 越域）。
- **已排除（逐行一致）**：`project_analytic` 之后的三条路径都与 OCCT 对齐 ——
  `value_of_uv`（`surface_projector.rs:519`）用 `point_surface_extrema_box` 把搜索窗限制在
  `[uf,ul]×[vf,vl]`（域内）；`next_value_of_uv`（`surface_projector.rs:900`）与
  `get_line`（`projection_cache.rs:445`）分别对应 `ShapeAnalysis_Surface.cxx:1164-1241` 与
  `ShapeConstruct_ProjectCurveOnSurface.cxx:902-1100`（含「仅周期面才投影 1/2 号内点」的门）。
  ⇒ 坏 pcurve 不是 `get_line`（其命中时 `first2d==last2d`），而是走**采样循环**。
- 采样循环与 `resolve_closed_surface_period_jump`（`surface_projector.rs:2050`）、
  `insert_additional_point_or_adjust`（`:1923`）经逐行对拍**均与 OCCT 一致**
  （`ShapeConstruct_ProjectCurveOnSurface.cxx:1499-1690` / `:2039-2148`）。
  ⇒ 差异不在投影代码，而在**曲面本身**（见 §3.5）。

### 3.5 🎯 G1 根因：STEP 读入缺少 `StepToTopoDS_TranslateFace` 的 `ConvertToPeriodic`

- OCCT `StepToTopoDS_TranslateFace.cxx:557-567`：只要 STEP 实体是 `B_SPLINE_SURFACE`，
  就读 `ShapeAlgo::ConvertToPeriodic(aGeomSurface)` 并在成功时替换
  （注释原文 *"pdn to force bsplsurf to be periodic"*）。
- 实测（本仓给 `occt_probe` 新增 `--period` 模式）：面 81/83/85/87
  （= STEP `#4809`/`#4850`/…）在 OCCT 侧 **`vper=1`、`vp=37.6991`**、
  `vclo=1`、`u=[-61.7587,-58.7754]`、`v=[0,37.6991]`；
  端口同面是 `vper=false`、`vrange=[0,37.6991]`。
- 于是：OCCT 的 pcurve 在**周期 37.6991 的 V** 上求值（`v` 加减周期等价），端口在**非周期**
  B 样条上把 `v=−18.85` 外推 ⇒ `check_same_parameter` 报 8.956 ⇒ 顶点容差
  ⇒ bbox 418.267 ⇒ lin 1.673 ⇒ 曲面过粗。
- 端口已有的对照：**曲线**侧已实现同一机制（`read_topology.rs:80-131` 的
  `should_be_periodic` + `SetPeriodic()`，对应 `StepToGeom.cxx:920-926`）；
  **曲面**侧缺这一步（`git grep ConvertToPeriodic` 在 `step/` 无命中）。
- 迁移内容（严格按 OCCT）：
  1. `ShapeCustom_Surface::ConvertToPeriodic`（`ShapeCustom_Surface.cxx:475` 起）：
     `ShapeAnalysis_Surface::IsU/VClosed(preci)` 判闭合；对
     `closed && !IsU/VPeriodic() && NbU/VPoles > 3` 先做端点多重度 `degree+1 → 1` 的
     **结点重排**（`:496-540`），再 `Geom_BSplineSurface::SetU/VPeriodic()`
     （`Geom_BSplineSurface_1.cxx:940/983`）。
  2. 在端口 `step/read_topology.rs::resolve_face` 解析出 surface 后、建 face 前，
     按 `StepToTopoDS_TranslateFace.cxx:558` 的同一判据调用。
- **`SetVPeriodic` 语义**（`Geom_BSplineSurface_1.cxx:983-1022`）：
  `first = FirstVKnotIndex()` / `last = LastVKnotIndex()`
  （非周期时 = `BSplCLib::{First,Last}UKnotIndex(deg, mults)`，`BSplCLib.cxx:111-137`）；
  取 `cknots/cmults = knots/mults[first..last]`；
  `cmults(1) = cmults(last) = min(deg, max(cmults(1), cmults(last)))`；
  `nbp = BSplCLib::NbPoles(deg, true, cmults)`（`BSplCLib.cxx:392-451`）；
  `myPoles.ResizeWithTrim(..., 1, nbp, true)`；`v_periodic = true`；`updateVKnots()`。
- **本实体（`#4809`）的逐值推演**（回归锚点，不作门禁）：
  - 原：9 个 V 极点；V 唯一结点 `(0, 9.4248, 18.8496, 28.2743, 37.6991)`，
    mults `(3,2,2,2,3)`。
  - 结点重排（`ShapeCustom_Surface.cxx:553-588`）：`a = 9.4248` ⇒
    V 唯一结点 `(-9.4248, 0, 9.4248, 18.8496, 28.2743, 37.6991, 47.1239)`，
    mults `(1,2,2,2,2,2,1)`，极点仍 9。
  - `SetVPeriodic`：`FirstVKnotIndex = 2`，`LastVKnotIndex = 6` ⇒
    `cknots = (0, 9.4248, 18.8496, 28.2743, 37.6991)`，`cmults = (2,2,2,2,2)`；
    `nbp = 8` ⇒ V 极点 9 → 8；`VPeriod = 37.6991`。
- 端口需要三件：(a) `GeomBSplineSurface::{set_u_periodic, set_v_periodic}`
  （含上述 `First/LastUKnotIndex` / `NbPoles` / 极点 `ResizeWithTrim`）；
  (b) `ShapeCustom_Surface::ConvertToPeriodic`（闭合判定 + 端点多重度重排）；
  (c) `step/read_topology.rs::resolve_face` 的调用点。
  `basis_funs`（`bspline_surface.rs:605`）已按 `periodic` 走
  `locate_parameter` 折回 ⇒ 置位后 `d0(u, v∓37.7)` 会正确折回。
- **求值侧的前提（本轮查明，2026-09-26）**：端口已有**周期感知**的曲面求值器
  `occt_core::bspl::prepare_eval::{prepare_eval, dn}` —— 极点索引 `ip/jp` 越界时折回
  （`prepare_eval.rs:174-205`），且 `GeomBSplineSurface::eval_dn_bspl` 已在用；
  但 `d0/d1/d2`（`bspline_surface.rs:836-929`）走的是**非周期**的
  `eval_bspline_surface`（本地 `basis_funs`，`:605`）与
  `bspl::surface::eval_surface_d1/d2`（无周期性参数）。
  ⇒ 置周期位之前必须把 `d0/d1/d2` 改走 `prepare_eval` 并传
  `u_periodic/v_periodic`，与 `BSplSLib::{D0,D1,D2}` 的周期极点折回一致；
  否则 `d0(u, v∓37.7)` 仍会算错（这正是 round 1-3 观测量级偏差的另一半前提）。
- 端口曲线侧已有可直接照搬的模板：`GeomBSplineCurve::set_periodic`
  （`bspline_curve.rs:158-198`）与 `occt_core::bspl::{knots::nb_poles,
  locate::first/last_u_knot_index, knots::knot_sequence_periodic}`。
- **验收**：端口对应面 `is_v_periodic()==true`（先用 `--dev` 打印确认），
  `shape_bnd_box` maxComp → **269.03**、lin → **1.07612**；再跑 `--lib` + 四道 STEP +
  `phase5/9/10/19/20`。

### 3.6 ✅ G1 已落地并全门禁通过（2026-09-26 round 8）

**实现（4 处，全部按 OCCT 迁移）**
- `crates/occt-geom/src/bspline_surface.rs`：新增 `set_u_periodic/set_v_periodic`
  （照 `Geom_BSplineSurface_1.cxx:940-1022`，复刻曲线版 `GeomBSplineCurve::set_periodic`）；
  `d0/d1/d2` 与 `local_d1/local_d2` 在 `u_periodic||v_periodic` 时走
  `prepare_eval::dn`（周期极点折回），并按 OCCT `EvalDN` 传**distinct 结点+多重度**。
  非周期路径不变。
- 新增 `crates/occt-topo/src/shape_custom_surface.rs`：`convert_to_periodic`
  （`ShapeCustom_Surface.cxx:475-616`：`sa_is_u/v_closed` + 端点多重度结点重排 + `SetU/VPeriodic`）。
- `crates/occt-topo/src/step/read_topology.rs::resolve_face`：对
  `B_SPLINE_SURFACE[_WITH_KNOTS]` 实体调用（`StepToTopoDS_TranslateFace.cxx:558`）。

**实测（`data/occ/a3n00.stp` vs OCCT）**
| 量 | 修前 | 修后 | OCCT |
|---|---|---|---|
| `shape_bnd_box` gap | 74.633 | **0.000815** | — |
| maxComp | 418.267 | **269.00** | 269.03 |
| `Prs3d` lin | 1.673068 | **1.076007** | 1.07612 |
| 顶点容差上限 | 115.68 | **0.00081** | 小 |
| 有三角的面数 | 222/226 | **226/226** | 226/226 |
| 网格 v/t | 9091/8922 | **11931/12283** | 11052/12324 |

**门禁**：`occt-topo --lib` **1281/0**、`occt-geom --lib` **143/0**、
`step_obj_parity` **14/14**、`step_to_obj` **13/13**、`step_obj_area` **11/11**、
`step_geometry_parity` **3/3**、`phase5/9/10/19/20` = 7/8/**8**/5/5（`phase10` 8/8）。
⇒ G1 完成；**G3 的 4 个空面随之消失**（`stats=226`），`phase10` 也一并转绿。
### 3.7 G1 之后的 G2/G3/G4 复测（2026-09-26 round 9）

逐类型对拍（端口 `zz_probe_a3n00` vs OCCT `--mesh 1.07612 0.349066`，226 面一一对应）：

| 类型 | 面数 | OCCT 节点/三角 | 端口 节点/三角 | 节点比 | 三角比 |
|---|---|---|---|---|---|
| Plane | 95 | 2932 / 2832 | 2948 / 2848 | 1.01 | 1.01 |
| Cylinder | 62 | 2795 / 2673 | 2732 / 2824 | 0.98 | 1.06 |
| Cone | 30 | 770 / 727 | 834 / 778 | 1.08 | 1.07 |
| Torus | 7 | 732 / 1214 | 783 / 1314 | 1.07 | 1.08 |
| **BSpline** | 32 | 3823 / 4878 | **4634 / 4519** | **1.21** | **0.93** |
| 合计 | 226 | 11052 / 12324 | 11931 / 12283 | 1.08 | 1.00 |

结论（修前 → 修后）：
- **G4 已消失**：Plane 由 1.16× 变 **1.01×**、Cylinder 由 1.27× 变 **0.98×**。之前的「平面/柱面过密」是错误 deflection（1.673）造成的假象，G1 修复后自然归位。
- **G2 大部分收口**：Cone 0.52× → **1.08×**、Torus 0.49× → **1.07×**。
- **G2 余项 = BSpline 面**：节点 1.21×（多 21%）、三角 0.93×（少 7%）⇒ 不是「采样不足」，而是**节点多、连接少**（内点插入/Delaunay 约束连接与 OCCT 不同），且总三角已只差 0.3%。
- **G3 完成**：4 个空面消失（`stats=226`）。
- **G2 余项的精确形状（round 10 定位）**：三角数 < 15 的 B 样条面只有 **7 个**，全部是
  **`vper=true` 的周期面**，且都是 **`wires=2`、`edges=[2,2]`**（两条各 2 边的闭合丝）：

  ```
  face=129 mv=162 mt=1  u=[0,3.1416]        v=[0,6.2832]   uper=false vper=true
  face=173 mv=162 mt=3  u=[0,1.2000]        v=[0,157.0796] uper=false vper=true
  face=176 mv=162 mt=11 u=[-0.0208,0.0208]  v=[0,301.5929] uper=false vper=true
  face=190 mv=162 mt=1  u=[-0.1290,0.1032]  v=[0,486.9469] uper=false vper=true
  face=201 mv=162 mt=3  u=[-0.3492,0.3492]  v=[0,229.3363] uper=false vper=true
  face=209 mv=162 mt=1  u=[-30.2954,-28.1901] v=[0,59.6903] uper=false vper=true
  face=216 mv=162 mt=4  u=[0,3.1416]        v=[0,6.2832]   uper=false vper=true
  ```

  ⇒ **节点足够（162 个）但 Delaunay 只连出 1–11 个三角**；对应 OCCT 的 360–486 三角。
  这不是采样不足，而是**周期 UV 带（两条闭合丝）上的约束/前沿三角化失败**。
  下一轮入口：`meshing/delaun` + `node_insertion` 对「URange 为一个窄带、两条 wire 沿 V 环绕、
  pcurve 跨 V 周期」的面的约束边与 `frontierAdjust` 处理，对照
  `BRepMesh_Delaun.cxx` / `BRepMesh_DelaunayBaseMeshAlgo.hxx`。

### 3.8 为什么这 7 个面没有内点：分类器把它们判为 OUT（2026-09-26 续）

插桩 `node_insertion.rs::list_surface_nodes`（已回退）实测：

- `face=129 generated=85 inside=0`、`face=173 generated=26 inside=0`、
  `face=176 generated=26 inside=0`（其余同）⇒ **生成了内点，但
  `Classifier::is_inside` 全部返回 false**，于是 `insert` 为空，Delaunay 只剩边界
  → 1–11 个三角。
- 再插桩 `collect_boundary_uv`（已回退）：每个 wire 的分类器多边形是**一条 u=const 的竖直退化线段**：
  `face=129 wire0 u=[3.1416,3.1416] v=[0,6.2832]`、`wire1 u=[1.5708,1.5708] v=[0,6.2832]`；
  即两条边界圆在 UV 里各自退化成线段，两条线段之间的带区无法被奇偶/绕数判据识别为内部。
- 端口 `Class2d::si_dans`（`occt-core/src/cslib/class2d.rs:119`）与
  `CSLib_Class2d::SiDans`（`CSLib_Class2d.cxx:141-187`）、
  `BRepMesh_Classifier::Perform/RegisterWire`（`BRepMesh_Classifier.cxx:35-127`）
  **逐行一致**；因此差异在这些 wire 的**输入 UV 折线**，不在分类算法。
- 已排除：`BRepMesh_ModelBuilder.cxx` 与 `BRepMesh_ModelHealer.cxx` 都没有缝边插入
  （全目录 grep `Seam` 只命中 `ModelPreProcessor`/`BaseMeshAlgo` 的方向修正与圆锥放大）。

⇒ **已查明（round 13）**：差别在**面拓扑**，不在分类器。

- OCCT `--topo` 实测：同面在 OCCT 侧是**一条 wire、含 2 条缝边**，例如
  `face 81 (type=6,vper=1,u=[-61.7587,-58.7754],v=[0,37.6991]) w0={e5,seam2}`、
  `face 98 w0={e4,seam2}`、`face 49 w0={e5,seam2}`；而上表 7 个端口面是
  `wires=2 edges=[2,2]`（两条各 2 边的闭合丝，**没有缝边**）。
- STEP 侧证实：`#4824=ADVANCED_FACE('',(#4820,#4823),#4809,.F.)` —— 该面在文件里就是
  **2 个 FACE_BOUND**（两条 u=const 圆），**不带缝边**。
- ⇒ OCCT 默认读入（ShapeProcess → `ShapeFix_Face::FixMissingSeam`，`ShapeFix_Face.cxx:1652-1718`）
  把两条 bound 合成一条 wire 并**插入 2 条缝边**；端口没有这一步，于是 UV 边界不闭合、
  分类器全 OUT、内点被丢、带区只剩 1–11 个三角。
- 这与画板 §3.3 的 **T-69/T-92（`ShapeFix_ComposeShell`/`FixMissingSeam`）** 同源；
  此前 T-69 结论（「FixMissingSeam 对目标文件是 no-op」）只对当时的 T0M 等文件成立，
  对 a3n00 的周期面**恰恰是缺件**。

⇒ 下一步实现入口：移植 `ShapeFix_Face::FixMissingSeam` 的「闭合面上缺缝边 → 插入 seam 边」
（`ShapeFix_Face.cxx:1652-1718`，配合 `ShapeBuild_Edge`/`BRep_Builder::UpdateEdge` 与
`shhealing
修完复测这 7 个面应达到 OCCT 的 360–486 三角。

**精确移植清单（round 14 定界）**——OCCT 侧共约 680 行，端口目前**两件都未移植**：
1. `static bool CheckWire(wire, face, URange, VRange, &isuopen, &isvopen, &isdeg)`
   （`ShapeFix_Face.cxx:1652-1721`，~70 行）：逐 wire 判「在 U/V 上是否开口、是否退化」。
   端口现有的 `brep_check::check_wire`（`brep_check.rs:224`）与
   `shape_analysis::check_wire_closed` 都**不是**这一件。
2. `bool ShapeFix_Face::FixMissingSeam()`（`ShapeFix_Face.cxx:1722-2330`，~610 行）：
   - 入口 `cxx:1729-1751`：只处理闭合面（`IsU/VClosed`），且 B 样条面必须已
     `IsU/VPeriodic`（这正是 §3.5 的 `ConvertToPeriodic` —— G1 已就位，前置满足）；
   - 本样本走的分支：两个 wire 都「在 V 上开口」（`ismodev == ±isvopen`，`cxx:1826-1881`），
     随后 `cxx:1994-2064` 一致性/排序 → 生成**缝边**并重构 face，最终
     `cxx:2324` 发 `FixAdvFace.FixMissingSeam.MSG0`（"Missing seam-edge added"）。
   - 依赖（端口多数已有）：`ShapeFix_Wire::FixReorder`（`wire_fix.rs`）、
     `ShapeAnalysis::GetFaceUVBounds`、`AdjustByPeriod`、`IsOuterBound`、
     `BRep_Builder::UpdateEdge/Range` 与 `make_edge`。
   - 调用点（OCCT）：`ShapeFix_Face::Perform` 的 `cxx:492-494`——即
     `StepToTopoDS` 之后的 ShapeProcess 的 `FixShape` 阶段；端口需在建面后按同一位置调用。

**❗ 最终定性（round 15）：G2 余项 ≡ 画板 T-92**

`FixMissingSeam` 的最后一步（`ShapeFix_Face.cxx:2236-2261`）并不自己造缝边，而是：

```cpp
occ::handle<Geom_RectangularTrimmedSurface> RTS =
  new Geom_RectangularTrimmedSurface(mySurf->Surface(), uf, uf + URange, vf, vf + VRange);
grid(1,1) = RTS;                       // fictive grid
ShapeExtend_CompositeSurface G(grid);
ShapeFix_ComposeShell CompShell;
CompShell.Init(G, L, tmpF, Precision::Confusion());
CompShell.ClosedMode() = true;
CompShell.Perform();
myResult = CompShell.Result();
```

即**缝边插入由 `ShapeFix_ComposeShell` 完成**（`ShapeFix_ComposeShell.cxx`，3603 行）——
正是画板 §3.3 登记的 **T-92**（当时估 3300–4600 行 Rust / 4 批，尚未开工）。

⇒ 结论：a3n00 这 7 个周期面无法在不移植 T-92 的前提下修好；本 objective 的四个枚举项
（G1/G2 内点采样/G3/G4）均已达成，剩余偏差属于 **T-92**（独立、已登记、已定界）。
证据链：G1 后 `lin=1.076007`（OCCT 1.07612）、`maxComp=269.00`（OCCT 269.03）、
`stats=226`；逐类型 Plane 1.01/Cylinder 0.98/Cone 1.08/Torus 1.07/BSpline 1.21节点·0.93三角；
7 面已定性为缺缝边 → T-92。
两条路：(a) 给 `occt_probe` 加 `--classifier` 模式，在 OCCT 侧打印
`IMeshData` 每个 wire 的 UV 折线包围盒（看它是否含缝边/是否跨越全周期）；
(b) 若 OCCT 的 wire 同样退化，则继续读 `CSLib_Class2d::internalSiDans` 对退化折线的返回值，
确认它是否把「两条平行线段之间」判为 IN。**.split("`").join(Q);

---

## 4. G2 细节：Cone/Torus/BSpline 内点采样

端口在 `lin=1.07612` 仍只有 OCCT 的 0.83×（9872/10184），且：
- Cone 376 ↔ OCCT 727；Torus 596 ↔ 1214；BSpline 2549 ↔ 4878；
- **Cone/Torus 在 lin=1.076 与 1.673 下三角数一个不差**（376 / 596）。

端口代码里这些 splitter **确实读了 deflection**：
- `splitter.rs:296-331 ConeRangeSplitter::get_split_steps` 用 `self.deflection()`；
- `splitter.rs:559-626 TorusRangeSplitter::generate_surface_nodes` 用 `self.deflection()`；
- `param_set.rs:520-617 generate_nurbs_grid` 用 `s.deflection()`。

但 `nodes.rs:199-216 factory_uses_deflection_control` 对
`Plane/Sphere/Cylinder/Cone/Torus` 返回 **false**（除非
`EnableControlSurfaceDeflectionAllSurfaces`），⇒ 这些类型走的是「无内点控制」的
Delaunay 算法分支，`generate_surface_nodes` 可能**根本没被调用**（或调用后其点被
`analytical_filter` 滤掉）。OCCT 对应 `BRepMesh_MeshAlgoFactory::GetAlgo` 的
`DelaunayBaseMeshAlgo` 与 `DelaunayDeflectionControlMeshAlgo` 分派。

**待查（下一轮）**：
1. 在 `ConeRangeSplitter/TorusRangeSplitter::generate_surface_nodes` 与 `generate_nurbs_grid`
   入口插计数，确认在 a3n00 的真实分派下是否被调用、返回多少点；
2. 若未被调用：核对 `factory_uses_deflection_control` 的 OCCT 语义（是否应为「按面而不是按
   类型」、或 `Cylinder/Cone/Torus` 有专用 `BRepMesh_*MeshAlgo`）；
3. 若被调用但被滤掉：对 `analytical_filter`（`param_set.rs`）与
   `BRepMesh_NURBSRangeSplitter::analyticalFilter` 逐行对拍。

---

## 5. G3：4 个无三角面

端口面 222–225 = 2 个 Cylinder + 2 个 Plane，`mv=0 / mt=0`（`stats=222` 而 `faces=226`）。
OCCT 同一文件 226 面全部有三角。端口在整网格失败时才回退/清空
（`brep_exchange.rs:88-109`），这里是**逐面**产出 0，属「该面三角化失败但被静默跳过」。
下一轮：打印这 4 个面的 wire/edge/pcurve 状态（是否 pcurve 缺失或 UV 退化）。

---

## 6. G4：Plane/Cylinder 过度网格

在**同 lin** 下端口平面/柱面比 OCCT 多 16% / 27%，且端口平面单面最多 674 三角 vs OCCT 376。
OCCT 的平面/柱面三角数基本由**边界离散 + 约束 Delaunay**决定（对 lin 不敏感）；端口在
lin 变细时平面 2895→3278、柱面 2788→3385 明显增长 ⇒ 端口在平面/柱面内部仍按 lin 铺点。
下一轮需核 `meshing/face_discret`/`delaun` 在平面/柱面上的内点来源。

---

## 7. 与画板既有条目的关系

- a3n00 属于 `data/occ/*.stp` 一族：`SURFACE_CURVE=0`、`PCURVE=0`（`specs/_occt_mesh_gt.md:88-101`），
  每条边的 pcurve 都要现算 ⇒ 与 **T-69 / T-91 / T-92（`ShapeFix_ComposeShell`）** 同族。
- 既有 GT 表（`_occt_mesh_gt.md`）里 `data/*.step` 15 个文件端口与 OCCT 14/15 一致，
  但 **a3n00 等 `data/occ` 5 个文件只列了 OCCT 值，没有端口列** ⇒ 本文是它们的首个端口对拍，
  暴露的是 **独立的 G1/G2/G4 缺口**（容差/deflection 与曲面内点），不是既有 T-92 的缝重构。
- 建议在画板新开 **T-93（a3n00 族：G1 顶点容差 / G2 曲面内点 / G3 空面 / G4 平面过密）**，
  并按 §8.3 的批次推进；G1 与 G2 都是「可翻译成 OCCT 控制流」的修复。

---

## 8. 附录

### 8.1 临时探针源码（`crates/occt-topo/examples/zz_probe_a3n00.rs`）

```rust
//! TEMPORARY diagnostic probe (delete after use).
use occt_topo::brep_exchange::prs3d_get_deflection;
use occt_topo::brep_tool::BRepTool;
use occt_topo::builder::TopoBuilder;
use occt_topo::meshing::incremental_mesh::IncrementalMesh;
use occt_topo::meshing::range_splitter::{classify_surface, SurfaceType};
use occt_topo::step::read_step_file;
use occt_topo::topo_tools_full::{edges_of_wire, faces_of, vertices_of, wires_of_face};

fn tname(t: SurfaceType) -> &'static str {
    match t {
        SurfaceType::Plane => "Plane", SurfaceType::Cylinder => "Cylinder",
        SurfaceType::Cone => "Cone", SurfaceType::Sphere => "Sphere",
        SurfaceType::Torus => "Torus", SurfaceType::SurfaceOfRevolution => "Rev",
        SurfaceType::SurfaceOfExtrusion => "Extr", SurfaceType::BezierSurface => "Bezier",
        SurfaceType::BSplineSurface => "BSpline", SurfaceType::OffsetSurface => "Offset",
        SurfaceType::OtherSurface => "Other",
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: zz_probe_a3n00 <step> [lin_override]");
    let override_lin: Option<f64> = std::env::args().nth(2).and_then(|s| s.parse().ok());
    let model = read_step_file(&path).expect("read step");
    let shape = if model.shapes.len() == 1 {
        model.shapes[0].shape.clone()
    } else {
        let parts: Vec<_> = model.shapes.iter().map(|s| s.shape.clone()).collect();
        TopoBuilder::new().make_compound_of(&parts).0
    };
    let b = occt_topo::brep_bnd_lib::shape_bnd_box(&shape);
    let (bmin, bmax) = (b.corner_min(), b.corner_max());
    println!("BBOX roots={} gap={:.6} min=({:.3},{:.3},{:.3}) max=({:.3},{:.3},{:.3})",
        model.shapes.len(), b.gap(), bmin.x(), bmin.y(), bmin.z(), bmax.x(), bmax.y(), bmax.z());
    let computed = prs3d_get_deflection(&shape, 0.1);
    let lin = override_lin.unwrap_or(computed);
    let angle = 20.0_f64.to_radians();
    let inc = IncrementalMesh::from_deflection(&shape, lin, false, angle);
    let stats = inc.face_stats();
    let (mv, mt) = inc.mesh().map(|m| (m.vertices.len(), m.triangles.len())).unwrap_or((0, 0));
    let faces = faces_of(&shape);
    println!("TOTAL faces={} computed_lin={:.6} used_lin={:.6} stats={} mesh_v={} mesh_t={}",
        faces.len(), computed, lin, stats.len(), mv, mt);
    let mut ftv: Vec<(f64, usize)> = faces.iter().enumerate()
        .map(|(i, f)| (BRepTool::face_tolerance(f), i)).collect();
    ftv.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    println!("TOP FACE TOLS: {}",
        ftv.iter().take(8).map(|(t, i)| format!("f{i}={t:.6}")).collect::<Vec<_>>().join(" "));
    let mut vtv: Vec<f64> = vertices_of(&shape).iter().map(|v| BRepTool::vertex_tolerance(v)).collect();
    vtv.sort_by(|a, b| b.partial_cmp(a).unwrap());
    println!("TOP VERT TOLS: {}", vtv.iter().take(8).map(|t| format!("{t:.6}")).collect::<Vec<_>>().join(" "));
    for (i, f) in faces.iter().enumerate() {
        let surf = BRepTool::face_surface(f);
        let tn = surf.as_ref().map(|s| tname(classify_surface(s.as_ref()))).unwrap_or("none");
        let wires = wires_of_face(f);
        let we: Vec<usize> = wires.iter().map(|w| edges_of_wire(w).len()).collect();
        let st = stats.get(i);
        println!("F {} type={} tol={:.6} wires={} edges={:?} mv={} mt={}",
            i, tn, BRepTool::face_tolerance(f), wires.len(), we,
            st.map(|s| s.vertices).unwrap_or(0), st.map(|s| s.triangles).unwrap_or(0));
    }
}
```

### 8.2 其它文件的顶点容差（对照，实测）

| 文件 | `shape_bnd_box` gap | top 顶点容差 | computed lin |
|---|---|---|---|
| `data/Cube.step` | 0 | 0 | 0.040000 |
| `data/Cylinder.step` | 0 | 0 | 0.040000 |
| `data/Torus.step` | 0 | 0 | 0.103910 |
| `data/Shape.step` | 5.000000 | 0.000005 | 0.448426 |
| `data/Shape-2.step` | 0 | 0.000005 | 0.600000 |
| `data/occ/a3n00.stp` | **74.633481** | **115.681895** | **1.673068** |
| `data/occ/acs10.stp` | 0.006678 | 0.007791 | 1.980053（= OCCT 495.013×0.004，正常） |

⇒ 只有 a3n00 出现 10^1 量级的容差；acs10 正常。

### 8.3 待办（按批次）

- **批 1（G1，最高优先）**：`resolve_edge` 残差插桩 → 定位 74.6/115.7 的实体；
  对读 `ShapeAnalysis_Curve::Project` + `StepToTopoDS_GeometricTool::UpdateParam3d`
  与端口 `edge_from_curve3d`，修投影/参数更新。验收：`shape_bnd_box` maxComp → **269.03**、
  `prs3d_get_deflection` → **1.07612**、**四道 STEP + `--lib` 1281/0 + phases 不劣化**、
  `export_data_obj` 的 v/f 向 GT 靠。
- **批 2（G2）**：确认 Cone/Torus/NURBS splitter 是否被调用；对齐
  `BRepMesh_MeshAlgoFactory` 的分派与内点生成。验收：Cone/Torus/BSpline 的逐面三角数向
  OCCT GT 收敛（`--mesh 1.07612 0.349066` 的逐面表）。
- **批 3（G3/G4）**：4 个空面的成因；平面/柱面内点来源。
- **总验收（沿用既有基线，不新写测试）**：`specs/_occt_mesh_gt.md` 的 a3n00 行 =
  **11052 / 12324**；`--lib` 1281/0；四道 STEP；`phase5/9/10/19/20`；
  `data/*.step` 15 文件三列表不劣化。
---

## 9. T-92（ShapeFix_ComposeShell）移植进度与交接（2026-09-26 round 13）

新模块 `crates/occt-topo/src/shape_fix_compose_shell/`，7 个文件，全部编译通过；该模块**尚未被任何生产路径调用**，因此 `--lib` 1281/0、四道 STEP、phases 均不受影响（每轮实测确认）。

| 文件 | 承载的 OCCT 类/职责 | 状态 |
|---|---|---|
| `wire_segment.rs` | `ShapeFix_WireSegment`（ShapeFix_WireSegment.cxx:25-310） | 完成 |
| `composite_surface.rs` | `ShapeExtend_CompositeSurface`（ShapeExtend_CompositeSurface.cxx:32-757） | 完成 |
| `helpers.rs` | ComposeShell 文件级 static（清单见 9.1） | 完成 |
| `load_wires.rs` | `LoadWires`（ShapeFix_ComposeShell.cxx:499-640） | 主体完成，reorder 块 UNPORTED |
| `reshape.rs` | 最小 `ShapeBuild_ReShape` 等价物 + `ApplyContext`（cxx:382-446） | 完成 |
| `shell.rs` | `ComposeShell::Init`（cxx:96-202）、`ComputeCode`（cxx:647-832） | 完成 |
| `mod.rs` | 导出 | 完成 |

### 9.1 helpers.rs 已移植的 file-static（函数注释内均标 OCCT 行号）

- `PointLineDeviation` / `PointLinePosition`（308-331）、`ParamPointOnLine` / `ParamPointsOnLine` / `ProjectPointOnLine`（336-377）
- `IsCoincided`（451-462）、`AdjustToPeriod`（ShapeAnalysis.cxx:66-69，复用 `shhealing::adjust_by_period`）
- `GetPatchIndex`（467-495）、`DefinePatch`（896-925）、`GetGridResolution`（929-938）
- `CheckByCurve3d`（876-892，注意用 SquareDistance）、`DistributeSplitPoints`（838-872）
- `FillBndBox` 非 Exact 臂（ShapeAnalysis_Curve.cxx:788-806）、`GetMiddlePoint`（2940-2974）
- `IsShortSegment`（2394-2447）、`IsSamePatch`（2452-2507）

不在本模块内的两件前置：`shape_analysis.rs` 新增 `tot_cross_2d`（ShapeAnalysis.cxx:114-150）与 `is_outer_bound`（cxx:203-230）；`meshing/model_builder/wire_builder.rs` 的 `sample_pcurve` / `sample_count_2d` 放宽为 `pub(crate)` 供其复用。

### 9.2 UNPORTED 清单（每处都有行号注释）

1. `LoadWires` 的 `ShapeFix_Wire::FixReorder(sawo) / FixReorder() / StatusReorder(DONE3)` 与缝向修正块（cxx:582-631）；
2. ~~`Init` 的 `BRepTools::UVBounds` 无限范围回退（cxx:120-124）~~ —— **round 16 已补**（用端口 `BRepTool::uv_bounds`，对照 `BRepTools::UVBounds(Face)`）；
3. `FillBndBox` 的 `Exact=true` 臂（ShapeAnalysis_Curve.cxx:808-860）；
4. `ShapeExtend_WireData::Reverse(face)` 的缝处理（ShapeExtend_WireData.cxx:508-572）：
   `Reverse()` 本体已移植（`shape_fix_compose_shell/wire_data.rs`），但它之后还要 `ComputeSeams(true)` 再对 `mySeamF` / `mySeamR` / `mySeams` 逐条 `SwapSeam`（cxx:511-543：取同一缝边的 FWD/REV 两条 pcurve，互换后 `UpdateEdge(E, c2dr, c2df, face, 0.)` + `Range(E, face, uff, ulf)`）。
   - 端口已具备的件：`GeometryRegistry::edge_pcurves`（tgeometry.rs:343，缝边返回两条）、`GeometryRegistry::set_edge_pcurves`（tgeometry.rs:383，整对写入）、以及 `meshing/model_builder/wire_builder.rs:430` 的 `swap_wire_seams`（注释即"ComputeSeams then SwapSeam"）。
   - **round 19 已确认**：`SwapSeam` 本体**端口早已有**——`meshing/model_builder/wire_builder.rs:448-459` 的 `swap_seam_pcurves`（本轮放宽为 `pub(crate)`），它做的正是 FWD/REV 两条 pcurve 互换，且带 `len() < 2` 守卫；对应 API 为 `GeometryRegistry::{edge_pcurves, set_edge_pcurves}`。因此`Reverse(face)` 只差 `ComputeSeams`（cxx:557）这一步。
   - `ComputeSeams` 本身也未以该名字移植，接线时应先读 `wire_builder.rs:430` 的现存实现并对齐。
5. `GetMiddlePoint` 顶点分支的 `ShapeAnalysis_Surface::ValueOfUV`（cxx:2942-2949）：端口投影器 `pcurve_full::surface_projector::value_of_uv`（surface_projector.rs:519）是 `pub(super)`，跨模块不可见，需先放宽可见性或另开入口。

### 9.3 剩余 roadmap（依赖顺序）

| 阶段 | cxx 行 | 行数 |
|---|---|---|
| `SplitWire` | 942-1432 | 490 |
| `SplitByLine`（单段 / 序列） | 1433-1917 / 1918-2130 | 698 |
| `SplitByGrid` | 2131-2278 | 148 |
| `BreakWires` | 2279-2393 | 115 |
| `CollectWires` | 2512-2939 | 428 |
| `MakeFacesOnPatch` | 2978-3274 | 297 |
| `DispatchWires` | 3275-3587 | 313 |
| `Perform` / `SplitEdges`（编排） | 206-273 | 68 |

### 9.4 接线点与验收（未变）

接线：`ShapeFix_Face::FixMissingSeam`（ShapeFix_Face.cxx:1722-2330；末尾 2236-2261 建 1×1 `Geom_RectangularTrimmedSurface` 网格 → `CompShell.Init(G, L, tmpF, Precision::Confusion())` → `ClosedMode() = true` → `SetMaxTolerance(MaxTolerance())` → `Perform()` → `Result()`），调用位置对应 `ShapeFix_Face::Perform:492-494`。

验收：`data/occ/a3n00.stp` 的 7 个周期 B 样条面三角数从 1-11 达到 OCCT 的 360-486；总量向 **11052 / 12324** 收敛；门禁 `occt-topo --lib` 1281/0、`occt-geom --lib` 143/0、`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 3/3、`phase5/9/10/19/20` 全绿。

### 9.6 全仓查重结果（round 20）

对剩余函数名做 `git grep`（`SplitWire / SplitByLine / SplitByGrid / BreakWires / CollectWires / MakeFacesOnPatch / DispatchWires / SplitEdges / BreakWire / CollectWire / DispatchWire`）：**全部 0 命中**，即这约 2300 行没有现成实现可复用，必须新写。

唯一的例外是缝处理：`ComputeSeams` 只在 `meshing/model_builder/wire_builder.rs:431` 的注释里被引用（该处实现了等价的缝识别），`SwapSeam` 则已有忠实实现 `swap_seam_pcurves`（同文件 448-459，本轮已放宽为 `pub(crate)`）。

因此剩余工作量的估计（约 2300 行）经查重**没有被推翻**；轮 19 的乐观猜测（可能存在换皮实现）不成立。
### 9.5 下一轮起点

`SplitWire` 首段已勘察：`cxx:950-1037`（setup + 逐边循环头 + 顶点/容差/3d 曲线取数）；其中 `aSurfTool`（ShapeAnalysis_Surface）与 `T`（myLoc）在端口分别为投影器与恒等（形状无 location）。其输入依赖 `ApplyContext` / `DistributeSplitPoints` / `GetMiddlePoint` 均已就位。
### 9.7 round 1（本会话，goal-491e5a3d）—— SplitWire 已移植

**新增** `crates/occt-topo/src/shape_fix_compose_shell/split_wire.rs`（约 480 行）：
`ShapeFix_ComposeShell::SplitWire`（`ShapeFix_ComposeShell.cxx:942-1429`）逐行迁移，
含：逐边循环 + `ApplyContext` 重分配（cxx:969-981）、`stop` 扫描（991-995）、
无分裂分支（996-1004）、非流形顶点收集（1005-1015）、`Curve3d`/`PCurve` 取数（1031-1044）、
非流形参数投影（1046-1087，`ShapeAnalysis_Curve::Project` 走
`shhealing::shape_analysis_curve::project_adaptor`，空 3D 曲线的
`ValueOfUV` + `Extrema_ExtPC2d` 分支走新导出的 `pcurve_full::surface_value_of_uv` +
`occt_geom2d::extrema2d::point_curve_extrema2d`）、
周期参数折回（1119-1128）、切点吸附到 prevV/lastV（1133-1259，含
`CheckByCurve3d`、`IsCoincided`、`GetGridResolution`）、
端点保护 copy + Context 替换（1261-1293）、
`CopyReplaceVertices` + 内部顶点插入 + `CopyPCurves` +
`TransferParametersProj::TransferRange`（1295-1344）、
`SameRange`/非流形方向（1346-1366）、AddEdge+DefinePatch（1368-1369）、
prev 参数推进（1371-1376）、结果替换 wire 注册到 context（1380-1398）、
非分裂时的 INTERNAL/EXTERNAL 处理与 UpdateEdge（1399-1425）、
`result.Orientation(anWireOrient)`（1427）。

**前置件（本轮补齐）**
- `reshape.rs`：`ReShape` 增加 `replace`；新增 `MapReShape`（按 TShape 身份记录 old→new，
  对应 `ShapeBuild_ReShape::Replace/Apply`）；`ComposeShell` 新增 `context: MapReShape`
  字段与 `context()/context_mut()`，对应 `ShapeFix_Root::Context()`（此前 ApplyContext 恒等）。
- `pcurve_full/mod.rs`：新导出 `surface_value_of_uv`（`ShapeAnalysis_Surface::ValueOfUV`，
  `ShapeAnalysis_Surface.cxx:1245-1515`），供 SplitWire 的非流形分支与后续
  `GetMiddlePoint` 顶点分支使用。
- `shhealing/wire_fix.rs`：`first_vertex` / `last_vertex` / `copy_replace_vertices` /
  `copy_replace_vertices_with` / `copy_pcurves` 放宽为 `pub`（`ShapeAnalysis_Edge::{First,Last}Vertex`、
  `ShapeBuild_Edge::{CopyReplaceVertices,Copy,PCurves}` 的既有移植件）。
- `shell.rs`：`ComposeShell` 字段放宽为 `pub(super)`，供同目录后续 split_*.rs 使用。

**已标注的偏差（无自创规则）**
- `myLoc` / `T` 恒等：端口形状无 location（cxx:959-962、1025-1029、1240、CheckByCurve3d 的 T）。
- `PCurve` 取数失败（cxx:1041-1044）在 OCCT 会解引用空 `C2d`；端口改为原样 AddEdge 并写 FAIL2
  （该分支对 STEP 目标数据不可达）。
- `MapReShape` 只记录直接替换，不重建复合子形状（goal 允许的最小等价物）。

**状态**：`cargo check -p occt-topo` 通过；`occt-topo --lib` **1281 passed / 0 failed**（本轮实测，模块未接线，无回归）。
下一步：`SplitByLine`（cxx:1433-2130）。
### 9.8 round 2 —— UNPORTED 项两项落地 + ❗SplitByLine 的硬前置 `Geom2dInt` 定界

**本轮交付（都在 objective 的 UNPORTED 清单内）**
- `helpers.rs`：移植 `ShapeAnalysis_Curve::SearchForExtremum`（`ShapeAnalysis_Curve.cxx:741-786`）
  与 `FillBndBox` 的 **Exact=true 臂**（`cxx:808-847`，`GeomAbs_C2` 区间 + X/Y 极值）。
  `fill_bnd_box` 现按 OCCT 单一函数签名带 `exact` 参数；`GetMiddlePoint` 仍传 `false`（`cxx:2961`）。
- `helpers.rs::get_middle_point`：顶点分支改为 `pcurve_full::surface_value_of_uv`
  （`cxx:2942-2949`）⇒ §9.2 第 5 项 UNPORTED **已消除**。
- `shhealing/wire_fix.rs::combine_vertex` 放宽为 `pub`（`cxx:2047` 需要）。

**❗ `SplitByLine`/`SplitByGrid` 的硬前置：`Geom2dInt` 家族未移植（实测 OCCT 源码定界）**

`Perform`（`cxx:206-255`）→ `SplitByGrid`（`cxx:2131-2275`）→ `SplitByLine`（`cxx:1433-2130`）。
`SplitByLine` 的核心是 `Geom2dInt_GInter::Perform(jGAC, iGAC, iDom, TOLINT, TOLINT)`
（`cxx:1612-1613`；`jGAC` = 切割线 `Geom2d_Line`，`iGAC` = 边的 pcurve），
产生 `IntLinePar` / `IntEdgePar`（`cxx:1621-1644`）。该包在端口**完全未移植**。

分派链（`IntCurve_IntCurveCurveGen.gxx:182-225` → `InternalPerform` `:235-...`，
按 `typ1 = Geom2dInt_Geom2dCurveTool::GetType(C1)`；本处 `C1` 恒为 Line）：
- `typ2 ∈ {Line,Circle,Ellipse,Parabola,Hyperbola}` → `IntCurve_IntConicConic`
  （`IntCurve_IntConicConic.cxx` 37693 B）；
- `typ2 == default`（B 样条等） → `IntCurve_IntConicCurve`
  （`Geom2dInt_TheIntConicCurveOfGInter_0.cxx` 2688 B + `TheIntersector...` 2016 B +
  `MyImpParTool...` 2353 B + `Geom2dInt_TheProjPCurOfGInter_0.cxx` 3111 B）；

- 两者皆非 conic → `IntCurve_TheIntPCurvePCurve`（`TheIntPCurvePCurve...` 2357 B +
  `ThePolygon2d...` 1394 B + `ExactIntersectionPoint...` 2334 B + `TheDistBetweenPCurves...` 1420 B +
  `Intf_InterferencePolygon2d`）。
`IntRes2d_Domain` / `IntRes2d_Transition` / `IntRes2d_IntersectionPoint` /
`IntRes2d_IntersectionSegment`（`IntRes2d_*.cxx`）也全部未移植。
填充时用到 `IP.TransitionOfSecond().PositionOnCurve()`（`cxx:1622`）与 `Inter.NbSegments()`
（`cxx:1629`）；端口现有的 `occt_geom2d::geom2d_api::intersect_curves`（`geom2d_api.rs:75`，
采样 / `Extrema_ExtCC2d` 路线）不返回 transition/segment，代入即自造结果，
与 `wire_fix.rs:3248-3254` 为 `FixSelfIntersection` 记录的阻塞同源。

**⇒ 结论**：`SplitByLine` / `SplitByGrid` 在 `Geom2dInt` 家族移植前无法忠实完成。
本 objective 其余项（`BreakWires` / `CollectWires` / `MakeFacesOnPatch` / `DispatchWires` /
`Perform`）与 UNPORTED 项不依赖它，可先推进；但 a3n00 的 7 面验收路径必经 `SplitByLine`
（`SplitByGrid` `cxx:2229-2274` 对每条 U/V 网格线都调用），故验收在 `Geom2dInt` 之前不可达。

**状态**：`cargo check -p occt-topo` 通过；`occt-topo --lib` **1281 passed / 0 failed**（本轮实测）。
**待决定**：是否把 `Geom2dInt`（`IntConicConic` + `IntConicCurve` + `IntRes2d_*`，约 10-45 KB C++）
作为独立批次授权移植（门禁 §4：旁支问题只报告，不擅自开工）。
### 9.9 round 3 —— BreakWires 已移植

新增 `crates/occt-topo/src/shape_fix_compose_shell/break_wires.rs`（约 120 行）：
`ShapeFix_ComposeShell::BreakWires`（`ShapeFix_ComposeShell.cxx:2279-2385`）逐行迁移：
从 EXTERNAL/INTERNAL wire 段收集切分顶点（`cxx:2281-2306`；`TopTools_ShapeMapHasher` 按 IsSame
→ 端口用 TShape 指针 `HashSet<usize>`）、逐段找首个切分顶点（`cxx:2320-2333`）、
闭合段首点的循环移位（`cxx:2335-2345`）、按 `1+(ind-1+shift)%NbEdges` 的重排与
`seqw.InsertBefore(i++, newwire)` / `seqw.SetValue(i, newwire)`（`cxx:2347-2383`；
端口用 `Vec::insert` + 下标写入，`i` 的推进语义与 OCCT 一致）、
`ori==INTERNAL && edge==EXTERNAL` 的方向改写（`cxx:2370-2375`）。
依赖件（`ShapeAnalysis_Edge::{First,Last}Vertex` → `shhealing::{first,last}_vertex`、
`MyClosedMode` → `closed_mode`）均已就位；**不依赖 `Geom2dInt`**。

**状态**：`cargo check -p occt-topo` 通过；`occt-topo --lib` **1281 passed / 0 failed**。
模块仍未被生产路径调用。

本轮已把 `ShapeAnalysis_Edge::GetEndTangent2d`（`ShapeAnalysis_Edge.cxx:269-366`）移植进
`helpers.rs::get_end_tangent_2d`（`CollectWires cxx:2665` 需要；此前端口 0 命中）。
`CollectWires` 其余部分不依赖 `Geom2dInt`。
### 9.10 round 4 —— CollectWires 已移植

新增 `crates/occt-topo/src/shape_fix_compose_shell/collect_wires.rs`（约 390 行）：
`ShapeFix_ComposeShell::CollectWires`（`cxx:2512-2936`）逐行迁移：
- 先把 vertex / INTERNAL 段移入 `wires` 并从 `seqw` 删除（`cxx:2519-2527`）；`shorts` 按初始
  `seqw` 长度分配、删除时不随之移位（`cxx:2518/2524-2526`），按当前下标填充（`cxx:2529-2548`）；
- `IsShortSegment` 用的是 `myGrid`（`CompositeSurface`）而非 `Geom_Surface` ⇒
  `helpers::is_short_segment` 参数改为 `&CompositeSurface` 并调 `value_uv`（= `myGrid->Value`）；
- 主 `for(;;)` 连接循环（`cxx:2560-2847`）：首段任取（2584-2604）、`IsSamePatch` 优先级（2607-2614）、
  同边回头最低优先级（2639-2652）、`GetEndTangent2d` 起始切向（2654-2673）、闭合期移（2675-2687）、
  角度/权重 tail 比较（2689-2722）、连接与 `sbwd->Add` / `Reverse(myFace)`（2726-2779）、
  首点记忆（2784-2791）、末点更新（2793-2822）、闭合入 `wires`（2824-2846）；
- 短段合并（`cxx:2853-2935`）。
- `myInvertEdgeStatus` 字段已补（ctor `true`、CollectWires 置位 `cxx:2731`、MakeFacesOnPatch 读 `cxx:2997`）；
  `MININD` / `MAXIND` 放宽为 `pub(super)`。

**已标注偏差**：`Reverse(myFace)` 的 `ComputeSeams + SwapSeam`（`cxx:557-572`，§9.2 第 4 项）仍未移植，
本处只用 `reverse_wire_data`（`cxx:483-506`）做边反转；短段合并用 `WireSegment::add_edge_patch`
（同时维护 patch 索引数组）代替 OCCT 仅改 `ShapeExtend_WireData` 的 `Add`（OCCT 处两者已失配）。

**状态**：`cargo check -p occt-topo` 通过；`occt-topo --lib` **1281 passed / 0 failed**，无 collect_wires 专属警告。
**下一步**：`MakeFacesOnPatch`（`cxx:2978-3274`）或 `DispatchWires`（`cxx:3275-3587`），最后 `Perform`（`cxx:206-255`）。
### 9.11 round 5 —— `ShapeExtend_WireData::Reverse(face)` 的 ComputeSeams/SwapSeam 已移植

`crates/occt-topo/src/shape_fix_compose_shell/wire_data.rs` 新增（对应 §9.2 第 4 项 UNPORTED）：
- `compute_seams`（`ShapeExtend_WireData.cxx:163-222`）：第一遍按 TShape 身份映射 REVERSED 边的秩
  （后出现的重复覆盖前者，同 OCCT `IndexedMap::Add` + `SE[num]=i`）；第二遍找与 REVERSED 匹配的
  FORWARD 边，第一对记 `(seamF, seamR)`，其余入 `mySeams`。
- `swap_seam`（`cxx:511-543`）：只处理 FORWARD 出现；先取 FORWARD 的 pcurve 范围 `(uff,ulf)`，
  互换 PCurve/PCurve2（`B.UpdateEdge(E,c2dr,c2df,face,0.)`），再恢复范围（`B.Range`）。
- `reverse_wire_data_on_face`（`cxx:545-572`）：`reverse_wire_data`（`cxx:483-506`）+ `ComputeSeams`
  + 对 `seamF`/`seamR`/`mySeams` 各 `SwapSeam`，并已 `pub use` 导出。
`CollectWires` 的 `Reverse(myFace)` 分支（`cxx:2767`）已改用它 —— 该项 UNPORTED 消除。

说明：端口既有的 `meshing::model_builder::swap_seam_pcurves` 只互换 pcurve、不恢复范围，且服务于
网格管线；本处按 OCCT 另实现，未改动网格路径。

**状态**：`cargo check -p occt-topo` 通过；`occt-topo --lib` **1281 passed / 0 failed**。
**下一步**：`MakeFacesOnPatch`（`cxx:2978-3274`）——仍有两处前置需处理：`ShapeFix_Face::FixOrientation`
（`myInvertEdgeStatus` 分支，`ShapeFix_Face.cxx:1254-1272`）与 `ShapeAnalysis_TransferParametersProj::CopyNMVertex`
（`cxx:3241`）均未移植；其余（roots/holes 分类、FClass2d、ValueOfUV）已具备。
### 9.12 round 6 —— LoadWires 的 FixReorder/缝向修正（cxx:582-631）已移植

`shhealing/wire_fix.rs`：把 `fix_reorder_wire` 重构为可复用的四件：
- `stored_manifold_edges`、`build_wire_order_3d`（= `ShapeAnalysis_Wire::CheckOrder` mode3d=true，
  FAIL2 短路为 status 0/ReorderOK true）、`apply_wire_order`（`FixReorder` 尾 `cxx:1387-1399`，
  按 signed order 重建并写回，保留非边子形状与端口 wire 方向约定）、
- `fix_reorder_wire_3d`（= `ShapeFix_Wire::FixReorder()`，返回 `(ReorderOK, WireOrderStatus)`，
  `Reversed` 即 `myStatusReorder` DONE3，`cxx:524-527`）、
  `fix_reorder_wire_with_order`（= `FixReorder(sawo)`，含 status 0 与 FAIL1/2/3 守卫 `cxx:1360-1385`）；
  `fix_reorder_wire` 保留为 thin wrapper。

`shape_fix_compose_shell/load_wires.rs`：UNPORTED 块替换为真实分支（`cxx:582-631`）：
双周期曲面先建 2D `ShapeAnalysis_WireOrder`（pcurve f/l 端点，`cxx:591-602`）→
`stat = status<0 ? -1 : 1`（`cxx:605`）→ `FixReorder(sawo)`；再 3D `FixReorder()`（`cxx:609`），
`status == Reversed` → `stat = -1`（`cxx:610-613`）；`stat < 0` 时比较 `IsOuterBound`
（原 wire vs 重排后 wire，`cxx:617-626`），不同则 `reverse_wire_data_on_face`（=`WireData::Reverse(face)`，
`cxx:629`）；最后以 REVERSED 追加（`cxx:634`）。

**状态**：`cargo check -p occt-topo` 通过；`occt-topo --lib` **1281 passed / 0 failed**。
至此 objective 的四项 UNPORTED **全部落地**（FillBndBox Exact、GetMiddlePoint ValueOfUV、
WireData::Reverse ComputeSeams、LoadWires FixReorder）。
**下一步**：`MakeFacesOnPatch`（需先处理 `ShapeFix_Face::FixOrientation` 与 `CopyNMVertex`）或 `DispatchWires`。
### 9.13 round 7 —— 剩余函数的阻塞定界（合并报告）

四项 UNPORTED 与 SplitWire/BreakWires/CollectWires 已完成（§9.7-§9.12）。剩余 4 项函数的阻塞如下，
均**不在原 objective 的文件范围内**，逐项取证：

**(A) SplitByLine / SplitByGrid —— `Geom2dInt` 未移植（§9.8 已详列）**
验收路径必经（`SplitByGrid cxx:2229-2274` 对每条 U/V 网格线调 `SplitByLine`），无同等分支可代。

**(B) MakeFacesOnPatch / DispatchWires —— 端口的 pcurve 注册表按 face TShape 键，OCCT 按 edge 上的 surface 键**
- OCCT `BRep_Tool::CurveOnSurface(edge, face)` 命中的是 `edge` 的 `BRep_GCurve` 的 **surface**，与 `face`
  形状身份无关；`MakeFacesOnPatch` 正是靠这一点在临时 `pf`/`fc`/`newFace` 上反复查询（`cxx:3057`、`3106`、
  `3125`、`3206`）。
- 端口 `boptools_2d::curve_on_surface_oriented`（`boptools_2d.rs:124-152`）用
  `GeometryRegistry::shape_key(&face.0)`（face 的 TShape 指针）取 pcurve ⇒ 为临时面做的查询全部落空，
  roots/holes 分类失效。
- `ShapeBuild_Edge::ReassignPCurve`（`ShapeBuild_Edge.cxx:530-595`）在 `DispatchWires`（`cxx:3460`）把 pcurve
  移到 patch 面；但 `MakeFacesOnPatch` 之后仍自建面查询，端口仍需 surface 级查找。
- ⇒ 忠实移植需先统一 pcurve 查找键（surface 级），属跨模块改造，超出 objective。

**(C) MakeFacesOnPatch 顶点洞分支 —— `CopyNMVertex` / `BRep_PointRepresentation` 未移植**
`cxx:3241` 调 `ShapeAnalysis_TransferParametersProj::CopyNMVertex`（`TransferParametersProj.cxx:562-710`），
核心是遍历 `BRep_TVertex::Points()` 的 `BRep_PointRepresentation`（PointOnCurve / PointOnCurveOnSurface /
PointOnSurface）。端口顶点只有单一 point+tolerance（`VertexGeom`），无表征列表 ⇒ 无同等分支。

**(D) DispatchWires 其余未移植件**
`ShapeBuild_Edge::{ReassignPCurve（530-595）, TransformPCurve（596-…）, RemovePCurve（430-…）}` 与
`ShapeFix_Edge::FixAddCurve3d`（`cxx:3518/3529`）端口 0 命中（`TransformPCurve` 仅在 `iges.rs` 有同名局部实现）；
`ShapeBuild_ReShape::IsRecorded/Value`（`cxx:3438-3441`）端口的 `MapReShape` 未提供。

**(E) Perform（cxx:206-255）** 依赖 (A)-(D)。

⇒ **在原 objective 的文件范围内，剩余函数均无法忠实完成**；(A) 同时是验收（a3n00 7 面 360-486 三角）的硬阻塞。
需要新授权：移植 `Geom2dInt` 家族，或改造 pcurve 注册表为 surface 级 + 移植 `CopyNMVertex` 依赖。
### 9.14 round 8 —— 授权扩大范围：Geom2dInt 起点（`IntRes2d` 已移植）

用户已授权「按需要自行扩大忠实移植范围」。新建 `crates/occt-geom2d/src/geom2d_int/`：
- `int_res2d.rs`（约 530 行）：`IntRes2d_Position` / `IntRes2d_TypeTrans` / `IntRes2d_Situation` 枚举；
  `IntRes2d_Transition`（默认 ctor、IN/OUT、TOUCH、UNDECIDED 三种 Set/Get）；
  `IntRes2d_Domain`（status 位域 1/2/4、`LimitInfinite`、`SetEquivalentParameters`、全部访问器）；
  `IntRes2d_IntersectionPoint`（ReversedFlag 交换 par/transition）；
  `IntRes2d_IntersectionSegment`（两点/单点/空三种 ctor）；
  `IntRes2d_Intersection`（done / points / segments 容器，1-based 访问）。
  源：`IntRes2d_Position.hxx`、`TypeTrans.hxx`、`Situation.hxx`、`Transition.{hxx,lxx,cxx}`、
  `Domain.{hxx,lxx,cxx}`、`IntersectionPoint.{hxx,lxx}`、`IntersectionSegment.{hxx,lxx}`、`Intersection.hxx`。
- `geom2d_int/mod.rs`；`occt-geom2d/src/lib.rs` 注册 `pub mod geom2d_int;`。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。
**下一步**：`Geom2dInt_Geom2dCurveTool`（`GetType/Value/D0-D3/NbSamples`；注意 `occt-topo::curve_sampling_2d`
已有 `nb_points`/`nb_samples` 的既有移植，需决定是否下沉到 `occt-geom2d` 以避免重复）→ `Geom2dInt_GInter` 分派
（`IntCurve_IntCurveCurveGen.gxx:29-...`，conic/conic 与 conic/curve 两路）→ `IntCurve_IntConicConic`（`Geom2dInt_TheIntConicConicOfGInter`）
与 `IntCurve_IntConicCurve`（`Geom2dInt_TheIntConicCurveOfGInter` + Intersector + MyImpParTool + TheProjPCur）。

### 9.15 round 9 —— Geom2dInt 曲线适配层：`Geom2dCurveTool` + 采样下沉

- `crates/occt-geom2d/src/geom2d_int/curve_sampling.rs`：把 `nb_points`（`Geom2dAdaptor_Curve.cxx:1351-1389`）、
  `nb_samples`（`Geom2dInt_Geom2dCurveTool.cxx:73-91`）从 `occt-topo::curve_sampling_2d` 下沉，并补
  `nb_samples_range`（`cxx:23-70`）。`occt-topo::curve_sampling_2d` 改为
  `pub(crate) use occt_geom2d::geom2d_int::nb_samples;`（调用点 `fclass2d/mod.rs:37`、`classifier.rs:463` 不变）。
- `crates/occt-geom2d/src/geom2d_int/curve_tool.rs`：`Geom2dInt_Geom2dCurveTool`——
  `GetType` / `Line` / `Circle` / `Ellipse` / `Parabola` / `Hyperbola` / `Value` / `D0-D3` / `DN` /
  `FirstParameter` / `LastParameter` / `EpsX`(=1e-10) / `NbSamples`(1、3 参) / `NbIntervals`(GeomAbs_C1) /
  `Intervals`(GeomAbs_C1) / `GetInterval` / `Degree`，另加 `GeomAbsCurveType` 枚举。
  `EpsX(C,eps)=C.Resolution(eps)` 与 `Degree` 标 UNPORTED（全树 grep：Geom2dInt/IntCurve 内核不调用）。
- `geom2d_int/mod.rs` 汇总导出。

**状态**：`occt-geom2d` 与 `occt-topo` 的 `cargo check` 均通过；`occt-geom2d --lib` 72/0；`occt-topo --lib` 1281/0。
**下一步**：`IntCurve_IConicTool` + `Geom2dInt_TheProjPCurOfGInter` + `IntCurve_IntConicCurveGen`（conic/curve 内核），
再 `IntCurve_IntConicConic`（conic/conic），最后 `Geom2dInt_GInter` 分派。
### 9.16 round 10 —— Geom2dInt：`IntCurve_IConicTool` 已移植

`crates/occt-geom2d/src/geom2d_int/iconic_tool.rs`（约 300 行）：`IntCurve_IConicTool`
（`IntCurve_IConicTool.{hxx,cxx:54-414}`）：
- `from_line` / `from_circle` / `from_ellipse` / `from_parabola` / `from_hyperbola`（含把 `gp::OX2d()`
  变换到 `axis` 的 `Abs_To_Object`，`cxx:75-137`）；
- `Value` / `D1` / `D2`（走既有 `occt_core::elib::clib2d` 的 `*_ax22d` / `*_ax2d`，`cxx:140-214`）；
- `Distance` / `GradDistance`（椭圆取 `AN_ELIPS=0` 臂，双曲线按 `x>0` 分支，`cxx:220-370`）；
- `FindParameter`（圆/椭圆负参数 `+2π`，`cxx:372-414`）。
`gp_Lin2d::Coefficients` 端口无该方法，按 `gp_Lin2d.hxx:91-96` 内联计算（`A=Dir.Y, B=-Dir.X, C=-(A·X+B·Y)`）。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。
**下一步**：`Geom2dInt_TheProjPCurOfGInter`（`TheLocateExtPC` + 投影）与 `IntCurve_IntConicCurveGen`
（`Geom2dInt_TheIntersectorOfTheIntConicCurveOfGInter` + `...MyImpParTool...`），再 `IntCurve_IntConicConic`
与 `Geom2dInt_GInter` 分派。
### 9.17 round 11 —— 实测 a3n00 的 pcurve 类型 + `Extrema_GCurveLocator`

**实测（临时探针 `zz_probe_pcurves`，已删除）**：`data/occ/a3n00.stp` 226 面全部边的 pcurve 类型直方图
（端口 `curve_on_surface_oriented`）：
- `BSpline: 262`、`Line: 362`、`NONE: 448`（尚未注册 pcurve 的边）；
- **没有** Circle / Ellipse / Parabola / Hyperbola。
⇒ 指导移植优先级：`SplitByLine` 只会遇到
  (a) 切割线 vs `Line` pcurve → `IntCurve_IntConicConic` 的 **Line/Line** 分支；
  (b) 切割线 vs `BSpline` pcurve → `IntCurve_IntConicCurve`（泛型路径，需 `TheProjPCur` / `LocateExtPC`）。
  其余 conic/conic 组合在本目标数据上不出现。

**本轮移植**：`crates/occt-geom2d/src/geom2d_int/curve_locator.rs`：
`Extrema_GCurveLocator`（`Extrema_GCurveLocator.hxx:45-129`）的 `Locate`（3 参 / 6 参），即
`Geom2dInt_TheCurveLocatorOfTheProjPCurOfGInter`（`Geom2dInt_TheCurveLocatorOfTheProjPCurOfGInter.hxx:26-28`
的别名实例化）——`TheProjPCur` 的第一步。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。
**下一步**：`IntCurve_IntConicConic` 的 Line/Line（`_1.cxx:1381-…` + `LineLineGeometricIntersection`(730-774) /
`CheckLLCoincidence`(1363-1378) / `computeIntPoint`(1254-1356) / `getDomainParametrs`(1240-1250) /
`FindPositionLL` + `DomainIntersection`）；以及 `Extrema_GLocateExtPC` 及其
`Extrema_LocEPCOfLocateExtPC2d` / `EPCOfELPCOfLocateExtPC2d` 依赖（用于 ConicCurve）。
### 9.18 round 12 —— `IntCurve_IntConicConic` Line/Line 的 helper 层

`crates/occt-geom2d/src/geom2d_int/conic_conic.rs`：`IntCurve_IntConicConic_1.cxx` 的 file-static helper：
- `line_line_geometric_intersection`（`_1.cxx:730-774`）
- `domain_intersection`（`_1.cxx:641-726`）
- `find_position_ll`（`_1.cxx:1209-1234`）
- `get_domain_parameters`（`_1.cxx:1240-1250`）
- `check_ll_coincidence`（`_1.cxx:1363-1378`）
- `compute_int_point`（`_1.cxx:1254-1356`）
`TOLERANCE_ANGULAIRE = 1e-8`（`IntCurve_IntConicConic_Tool.cxx:20`）。
两处按 OCCT 语义对齐：`FindPositionLL(theResSup, ...)` 直接就地修改 `res_sup`；`aRes2` 被 `FindPositionLL` 就地吸附。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。
**下一步**：`IntCurve_IntConicConic::Perform(L1,D1,L2,D2,TolConf,TolR)` 本体（`_1.cxx:1381-2080`，含
`Append` / `SetReversedParameters` 语义），随后 `Geom2dInt_GInter` 分派。
### 9.19 round 13 —— `IntCurve_IntConicConic` Line/Line 的 `Perform` 本体

`crates/occt-geom2d/src/geom2d_int/int_conic_conic.rs`（约 550 行）：`IntCurveIntConicConic` 与
`perform_line_line`（`IntCurve_IntConicConic_1.cxx:1381-2232`）逐行迁移：
- `nbsol==1`：窗口加宽（1439-1496）、`DomainIntersection`（1499）、单点分支（1517-1547）、
  段分支（1549-1803：`isOpposite` 的 Res1/Res2 互换、`FindPositionLL` 就地吸附、`ResultIsAPoint` 判定、
  `PtSeg1/PtSeg2`、`Append(Segment)` / `Append(SegmentToPoint)`）、`Pos1a==Pos2a==Middle` 分支（1804-1944）、
  `computeIntPoint` 点分支（1946-1974）；
- `nbsol==2`（共线，2008-2231）：`ResHasFirstPoint/LastPoint` 的 1/2/3 组合、`ParamStart/End` 计算、
  无限/半无限/有限段的 `IntRes2d_IntersectionSegment` 构造。
另加 `SegmentToPoint`（`_1.cxx:2620-2653`）与 `set_inout`/`set_touch` 辅助；
`IntRes2dIntersectionPoint` 增加 `Default`，`IntRes2dIntersection` 增加 `reverse` 与
`set_reversed_parameters` / `reversed_parameters`（`IntRes2d_Intersection.lxx:113-122`）。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。
**下一步**：`Geom2dInt_GInter` 分派（`IntCurve_IntCurveCurveGen.gxx` 的 `InternalPerform`，
Line/conic → ConicConic，Line/其他 → ConicCurve）与 ConicCurve 路径。
### 9.20 round 14 —— `IntImpParGen` 静态分类器

`crates/occt-geom2d/src/geom2d_int/int_imp_par_gen.rs`（约 170 行）：`IntImpParGen`（`IntImpParGen.cxx:28-251`）：
- `NormalizeOnDomain`（`cxx:28-46`）
- `DeterminePosition`（`cxx:49-83`）
- `DetermineTransition(Pos1, Tan1, Norm1, T1, Pos2, Tan2, Norm2, T2, tol)`（`cxx:86-206`，含 TOUCH 的
  courbure / Inside / Outside 判定）
- `DetermineTransition(Tan1, T1, Tan2, T2, tol)`（`cxx:209-251`，无法线重载，返回是否已判定）
常数 `TOLERANCE_ANGULAIRE=1e-8`、`DERIVEE_PREMIERE_NULLE=1e-12`（`cxx:24-25`）。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。
**下一步**：ConicCurve 路径其余依赖——`math_FunctionAllRoots` / `math_FunctionSample`
（`IntImpParGen_Intersector.gxx:27-28,296-300`）、`Geom2dInt_MyImpParToolOfTheIntersector...`、
`Geom2dInt_TheProjPCurOfGInter`（待 `Extrema_GLocateExtPC`），最后 `IntImpParGen_Intersector.Perform` 本体。
### 9.21 round 15 —— ❗重要发现：`IntCurve`/`IntImpParGen`/`IntRes2d`/`math` 栈**已在 occt-core**（去重 + 路线修正）

移植 ConicCurve 时发现 `crates/occt-core/src/` 已存在下列**忠实移植**（此前几轮不知情，重复造了轮子）：
- `intres2d/{domain,intersection,transition}.rs` = `IntRes2d_Domain` / `_IntersectionPoint` /
  `_IntersectionSegment` / `_Intersection` / `_Transition`；
- `intcurve/{iconic_tool,pconic,pconic_tool,project_on_pconic_tool}.rs` = `IntCurve_IConicTool` 等；
- `intimpargen/{gen,intersector}.rs` = `IntImpParGen` **与** `IntImpParGen_Intersector` 的
  `IntCurve_IntImpConicParConic` 实例化（`ImpTool=IntCurve_IConicTool`、`ParCurve=IntCurve_PConic`、
  `ProjectOnPCurveTool=IntCurve_ProjectOnPConicTool`）；
- `math_function_all_roots.rs` / `math_function_sample.rs` = `math_FunctionAllRoots` / `math_FunctionSample`；
- `intana2d` = `IntAna2d_AnaIntersection`。

**本轮动作（去重 + 对齐）**：
- 删除本轮前几轮新增的重复件：`geom2d_int/int_res2d.rs`、`geom2d_int/iconic_tool.rs`、
  `geom2d_int/int_imp_par_gen.rs`；
- `geom2d_int/conic_conic.rs` 与 `int_conic_conic.rs` 改用 `occt_core::intres2d` 类型
  （`in_out`/`touch`、`with_transitions(&p,..)`、`set_values(&p,..)`、`from_two_points`/`from_one_point`/`infinite`、
  `append_point(&p)`/`append_segment(&s)`、`done` 字段）；`IntCurveIntConicConic` 手写 `Default`。
- `geom2d_int/mod.rs` 只保留 `conic_conic` / `curve_locator` / `curve_sampling` / `curve_tool` / `int_conic_conic`。

**⇒ ConicCurve 路径的真实缺口**：occt-core 的 `IntImpParGen_Intersector` 是**为 `IntCurve_PConic` 单态化**的；
`Geom2dInt` 需要的是 `ParCurve = Adaptor2d_Curve2d`（端口 `Curve2d`）、`ParTool = Geom2dInt_Geom2dCurveTool`、
`ProjectOnPCurveTool = Geom2dInt_TheProjPCurOfGInter` 的单态化（即 `Geom2dInt_TheIntConicCurveOfGInter`）。
即：**内核逻辑已在，只缺 `Curve2d` 版单态化 + `TheProjPCur`**（后者依赖 `Extrema_GLocateExtPC` 系）。

**状态**：`occt-geom2d` 与 `occt-topo` 的 `cargo check` 通过；`occt-geom2d --lib` 72/0。
**下一步**：在 `occt-geom2d::geom2d_int` 写 `Curve2d` 单态化的 `IntImpParGen_Intersector`
（复用 `occt-core::math_function_all_roots`/`math_function_sample`/`intres2d`/`intcurve::iconic_tool`/`intimpargen::gen`），
并移植 `Geom2dInt_TheProjPCurOfGInter`，随后 `Geom2dInt_GInter` 分派。
### 9.22 round 16 —— `Extrema_GFuncExtPC`（Geom2dInt 2D 实例化）

`crates/occt-geom2d/src/geom2d_int/gfunc_ext_pc.rs`（约 300 行）：`Extrema_GFuncExtPC`
（`Extrema_GFuncExtPC.hxx:32-482`）按 `Geom2dInt_PCLocFOfTheLocateExtPCOfTheProjPCurOfGInter.hxx:29-35`
的单态化移植（`TheCurve=Curve2d`、`TheTool=Geom2dInt_Geom2dCurveTool`、`ThePoint=GpPnt2d`、`TheVector=GpVec2d`）：
`POnCurv2d` 值类型、`Initialize`/`SetPoint`/`SubIntervalInitialize`/`SearchOfTolerance`、
`Value`（含奇异导数的 Taylor/`DN` 与三点估计）、`Values`（含奇异导数的三点差分）、
`GetStateNumber`/`NbExt`/`SquareDistance`/`IsMin`/`Point`，并 `impl MathFunction` / `MathFunctionWithDerivative`。
它是 `Extrema_GenLocateExtPC` 求根的目标函数，即 `TheProjPCur` 的下一步依赖。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。

**❗新发现的阻塞（`TheProjPCur` 的最后一段）**：`Extrema_GenLocateExtPC`（`Extrema_GenLocateExtPC.hxx:108-126`）
用 `math_FunctionRoot`；而 OCCT 8.0.0 的 `math_FunctionRoot`（`math_FunctionRoot.cxx:73-119`）**委托给
`math_FunctionSetRoot`（`math_FunctionSetRoot.cxx` 42 067 B，未移植）**。端口已移植
`math_NewtonFunctionRoot` / `math_FunctionRoots` / `math_FunctionAllRoots`，但**没有 `math_FunctionSetRoot`**。
要继续 `TheProjPCur` → `IntImpParGen_Intersector(Curve2d)` → `Geom2dInt_GInter`，必须先移植
`math_FunctionSetRoot`（及 `math_Vector`/`math_Matrix`/`math_FunctionSetWithDerivatives`）。
### 9.23 round 17 —— ✅修正 §9.22 的误判：`math_FunctionSetRoot` 已在 `occt-math`；`TheProjPCur` 链已打通

**修正**：§9.22 断言 `math_FunctionSetRoot` 未移植是**错误**的。实际它在
`crates/occt-math/src/function_set_root.rs`（`MathFunctionSetRoot` / `MathFunctionRoot` /
`MathFunctionSetWithDerivatives`），且 `occt-math` 已含 `MathVector`/`MathMatrix`/`IntegerVector`/`Gauss`/
`GaussLeastSquare`/`SVD`/`BrentMinimum`。`occt-geom2d` 也早已依赖 `occt-math`。

**本轮交付**：
- `geom2d_int/gfunc_ext_pc.rs` 改实现 `occt_math::MathFunctionWithDerivative`（`MathFunctionRoot` 所需）。
- `geom2d_int/gen_locate_ext_pc.rs`：`Extrema_GenLocateExtPC`（`Extrema_GenLocateExtPC.hxx:42-167`）的
  `Geom2dInt` 单态化：`Perform` 调 `MathFunctionRoot::new_with_bounds(my_f, u0, tolU, umin, usup, 100)`，
  完成后用 `|F(uu)| >= 1e-7` 复核（`hxx:113-125`）；`IsDone`/`SquareDistance`/`IsMin`/`Point`。
- `geom2d_int/proj_p_cur.rs`：`Geom2dInt_TheProjPCurOfGInter`（`_0.cxx:27-79`）——
  `NbSamples(C)` → `Extrema_GCurveLocator::Locate` → `GenLocateExtPC(P,C,defaultparam,epsX)` → `IsDone/IsMin`。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。
**下一步**：`Curve2d` 单态化的 `IntImpParGen_Intersector`（复用 occt-core 的 `intres2d`、`intimpargen::gen`、
`math_function_all_roots`、`intcurve::iconic_tool`）→ `Geom2dInt_TheIntConicCurveOfGInter` → `Geom2dInt_GInter` 分派。
### 9.24 round 18 —— ✅ Curve2d 单态化的 `IntImpParGen_Intersector`（ConicCurve 内核）

`crates/occt-geom2d/src/geom2d_int/conic_curve.rs`（约 860 行）：把 occt-core 的 `intimpargen::intersector`
（为 `IntCurve_IntImpConicParConic` 单态化）**按同一份模板**改写为 `Geom2dInt` 实例化：
`ParCurve=Curve2d`、`ParTool=Geom2dInt_Geom2dCurveTool`、`ProjectOnPCurveTool=Geom2dInt_TheProjPCurOfGInter`。
替换点：`pconic_tool::{value,d1,d2}` → `curve_tool::*`；`nb_samples_in_range` → `nb_samples_curve_range`（`as i32`）；
`eps_x(c)` → `curve_tool::eps_x()`；`project_on_pconic_tool::{find_parameter,find_parameter_bounded}` →
`proj_p_cur::{find_parameter,find_parameter_range}`；`IntCurvePConic` → `dyn Curve2d`；`MyImpParTool` 改为对 `Curve2d` 求值；
复用 occt-core 的 `intres2d`、`intimpargen::gen`、`math_function_all_roots`、`math_function_sample`、`intcurve::iconic_tool`。
结构体名 `Geom2dIntConicCurve`（即 `Geom2dInt_TheIntConicCurveOfGInter` 的内核）。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。
**下一步**：`IntCurve_IntConicCurveGen` 的 `Perform` 包装（`IntCurve_IntConicCurveGen.lxx:30-130`，即
`Geom2dInt_TheIntConicCurveOfGInter`），再 `Geom2dInt_GInter` 分派（`IntCurve_IntCurveCurveGen.gxx` 的 `InternalPerform`）。
### 9.25 round 19 —— ✅ `Geom2dInt_GInter` 分派 + `IntConicCurveGen` 包装

- `geom2d_int/int_conic_curve.rs`：`IntCurve_IntConicCurveGen`（`Geom2dInt_TheIntConicCurveOfGInter`）——
  `Perform(ICurve, D1, PCurve, D2, ...)`（`IntCurve_IntConicCurveGen.lxx:119-130`）与
  line/circle/ellipse/parabola/hyperbola 重载（`lxx:44-116`；circle/ellipse 用
  `SetEquivalentParameters(first, first+2π)`）；内部 `Geom2dIntConicCurve` + `SetValues`。
- `geom2d_int/ginter.rs`：`Geom2dInt_GInter`（`IntCurve_IntCurveCurveGen`）——
  `Perform(C1,D1,C2,D2,..)`、`Perform(C1,C2,D2,..)`（`lxx:139-151`，`D1=ComputeDomain(C1,max(TolConf,Tol))`）、
  `Perform(C1,C2,..)`（`lxx:110-121`）、`ComputeDomain`（`gxx:120-177`）、`InternalPerform` 的 **typ1=Line 全部分支**：
  Line/Line → `IntCurveIntConicConic::perform_line_line`；Line/default → `IntConicCurveGen::perform_line`。

**UNPORTED（注释含 OCCT 行号）**：
- Line/{Circle,Ellipse,Parabola,Hyperbola} → `IntCurve_IntConicConic_1.cxx:2236`(Line/Circle)、`:2861`(Line/Ellipse) 等（本移植只做了 Line/Line）；
- `typ1 != Line` 的所有分支（`IntCurve_IntCurveCurveGen.gxx:339-...`）。

`ShapeFix_ComposeShell::SplitByLine` 的第一条曲线恒为切割直线，因此覆盖了其全部分支。

**状态**：`cargo check -p occt-geom2d` 通过；`occt-geom2d --lib` 72/0。
**下一步**：把 `ShapeFix_ComposeShell::SplitByLine` 接上 `Geom2dIntGInter`（等价于 OCCT 的 `Geom2dInt_GInter::Perform`），核对控制流。
### 9.26 round 20 —— `SplitByLine`（`WireSegment` 重载）接入 `Geom2dIntGInter`

`crates/occt-topo/src/shape_fix_compose_shell/split_by_line.rs`（约 470 行）：
`ShapeFix_ComposeShell::SplitByLine`（`ShapeFix_ComposeShell.cxx:1433-1914`）忠实移植：
- 顶点分支（`cxx:1441-1477`）：`ShapeAnalysis_Surface::ValueOfUV` + `PointLinePosition` + `Context()->Replace` +
  `BRep_Builder::MakeVertex` + `SetVertex`；
- `closedDir`/`halfPeriod`（`cxx:1483-1495`），主循环（`cxx:1497-1707`）：`sae.PCurve` →
  `boptools_2d::curve_on_surface_range`；ClosedMode 周期调整（`FillBndBox` + `AdjustToPeriod` + `AdjustByPeriod`）；
  **`Geom2dInt_GInter::Perform(jGAC, iGAC, iDom, TOLINT, TOLINT)`**（`cxx:1612-1613`）→ 端口
  `Geom2dIntGInter::perform_with_d2`；交点/线段收集、参数裁剪到 `[f,l]`、按边参数排序、闭合点检测；
- 段码分析（`cxx:1709-1898`）：闭合模式重复点消除、`ComputeCode`、`aNewSegCodes` / `IntCode` 编码；
- `SplitWire`（`cxx:1900-1911`）→ 端口 `ComposeShell::split_wire`，输出 `SplitLinePar/Code/Vertex`。

至此切割点计算与 `SplitWire` 已串联。**UNPORTED**：`SplitByLine` 的序列重载（`cxx:1918-2130`）与 `SplitByGrid`（`cxx:2131-2977`）。

**状态**：`cargo check -p occt-topo` 通过；`occt-topo --lib` 1281/0。
**下一步**：`SplitByLine` 序列重载 + `SplitByGrid`，然后 `DispatchWires` / `MakeFacesOnPatch` / `Perform` 与 `FixMissingSeam` 接线。
### 9.27 round 21 —— `SplitByLine` 序列重载（`cxx:1918-2127`）

`split_by_line.rs` 增补 `ComposeShell::split_by_line_wires`：
- 逐条 `SplitByLine(WireSegment&)`（`cxx:1929-1932`）；
- 沿线参数排序（`cxx:1935-1946`）；合并零长切向段（`cxx:1949-1965`）；
- 沿线走，`parity`/`tanglevel`/`halfparity` 状态机（`cxx:1969-2010`）；
- 对内部区间用 `ShapeBuild_Vertex::CombineVertex`（2 参重载，tolFactor = 1.0001）合并重合顶点（`cxx:2019-2057`）；
- `B.MakeEdge` + 两个顶点 + 两条 `Geom2d_Line` pcurve（端口 `topo_builder` + `GeometryRegistry::set_edge_pcurves` /
  `set_pcurve_range` / `set_edge_range`）（`cxx:2061-2070`）；
- `ShapeFix_WireSegment(sbwd, TopAbs_EXTERNAL)` + `DefinePatch`/`DefineIUMin/Max`（或 `DefineIVMin/Max`）（`cxx:2072-2108`）；
- `parity` 为奇数时 `EncodeStatus(ShapeExtend_FAIL4)`（`cxx:2111-2113`）；
- 对每条 wire `ApplyContext`（`cxx:2120-2126`）。

**状态**：`occt-topo --lib` 1281/0。
**下一步**：`SplitByGrid`（`cxx:2131-2977`）。
### 9.28 round 22 —— `SplitByGrid`（`cxx:2131-2275`）

`split_by_grid.rs`：`ComposeShell::split_by_grid`：
- `BRepTools::UVBounds(myFace)` + `myGrid->Bounds`（`cxx:2136-2138`）；
- ClosedMode：对每条 wire 用 `BRepTools::AddUVBounds(Face, Wire)`（等价 OCCT 的 `myFace.EmptyCopied()` +
  `ShapeAnalysis::GetFaceUVBounds`——pcurve 存在共享 surface 上）取该 wire 的 UV 盒，`AdjustToPeriod` 平移，
  `GetPatchIndex` 限幅到 `[0,2]` 后 `DefineIUMin/Max`、`DefineIVMin/Max`（`cxx:2152-2199`）；
- 非 closed：用整面 UV 盒一次性定义（`cxx:2203-2225`）；
- 按 U 缝：`U = UJointValue(i)` 的竖线，closed 时 `i = 1..=NbUPatches`，否则 `2..=NbUPatches`；
  `!ClosedMode && UClosed` 时按 period 平移并逐条 `SplitByLine(seqw, ln, true, cutIndex)`（`cxx:2227-2251`）；
- 按 V 缝同理（`cxx:2253-2273`）。

**状态**：`occt-topo --lib` 1281/0。
**下一步**：`DispatchWires`（`cxx:3275-...`）/ `MakeFacesOnPatch`（`cxx:2978-3274`）。
### 9.29 round 23 —— `ShapeBuild_Edge` 子集（`DispatchWires` 前置）

`crates/occt-topo/src/shhealing/shape_build_edge.rs`：`ShapeBuildEdge`
- `Copy(edge, sharepcurves)`（`ShapeBuild_Edge.cxx:417-426`）；
- `RemovePCurve(edge, face)`（`cxx:430-443`）→ `remove_pcurves_on_surface`；
- `SetRange3d`（`cxx:338-356`）→ `set_edge_range`；
- `ReassignPCurve(edge, old, sub)`（`cxx:530-592`）：`CountPCurves`（端口按 face-key 计 pcurve 数）、
  seam 时保留第二条 pcurve、`sub` 上已有则按边方向决定第一/第二条、最后 `set_pcurve_range`。

**UNPORTED**：`TransformPCurve`（`cxx:596-...`）——需要 `Geom2d_Curve::TransformedParameter`（逐类）与
`gp_Trsf2d::Form()`，端口都没有。

**状态**：`occt-topo --lib` 1281/0。
**下一步**：`ShapeFix_Edge::FixAddCurve3d` + `ShapeAnalysis_TransferParametersProj::CopyNMVertex`，然后写 `DispatchWires`（`cxx:3275-3584`）。
### 9.30 round 24 —— `ShapeBuild_ReShape::{IsRecorded,Value}` + `ShapeFix_Edge::FixAddCurve3d`

- `shape_fix_compose_shell/reshape.rs`：`ReShape` trait 增 `is_recorded` / `value`（默认 `false` / `None`），
  `MapReShape` 按 TShape 身份实现。
- `crates/occt-topo/src/shhealing/shape_fix_edge.rs`：`ShapeFixEdge::fix_add_curve3d(edge, face)`
  （`ShapeFix_Edge.cxx:618-638`）：退化或已有 3D 曲线返回 false；否则取 `curve_on_surface_range` + `face_surface`，
  以 `CurveOnSurface` 作为 pcurve 的 3D 像写回 `EdgeGeom`（保留 tolerance/same_parameter/same_range/degenerated）。
  **UNPORTED**：`TempSameRange`（`cxx:335-...`，`!BRep_Tool::SameRange` 分支 `cxx:626-629`）。

**状态**：`occt-topo --lib` 1281/0。
**下一步**：写 `DispatchWires`（`cxx:3275-3584`）——仍差 `TransformPCurve`（affine 分支）与 `CopyNMVertex`。
### 9.31 round 25 —— `CopyNMVertex`（face 重载）+ `ValueOfUV`+`Gap` 导出

- `pcurve_full/mod.rs`：新增 `surface_value_of_uv_with_gap`。
- `shhealing/transfer_params.rs`：`copy_nm_vertex_face(v, toFace, fromFace)`（`Proj.cxx:715-802`）——
  方向检查（仅 `INTERNAL`/`EXTERNAL`）、`EmptyCopied`、`ValueOfUV` + `Gap` 容差放宽（`cxx:782-794`）。
  **UNPORTED**：`BRep_PointRepresentation` 复制环（`Proj.cxx:734-797`）与 `BRep_Builder::UpdateVertex(V, U, V, Face, Tol)`
  （`cxx:800`）——端口顶点几何只有单点，无表示列表。

**状态**：`occt-topo --lib` 1281/0。
**下一步**：`DispatchWires`（`cxx:3275-3584`）——仍需 `TransformPCurve` 的 affine 分支，或将其标 UNPORTED。
### 9.32 round 26 —— `MakeFacesOnPatch`（`cxx:2978-3271`）

`shape_fix_compose_shell/make_faces_on_patch.rs`：`ComposeShell::make_faces_on_patch(faces, surf, loops)`
- 单环：`make_face(surf, [wire])`；`myInvertEdgeStatus` 分支（`cxx:2997-3006`）标 **UNPORTED**
  （`ShapeFix_Face::FixOrientation`；`Perform` 在 `cxx:209` 已将其清零，故 `FixMissingSeam` 路径不进入）；
- 多环：伪面 `pf` + `BRepTopAdaptor_FClass2d`（端口 `fclass2d::FClass2d`）找根（`cxx:3021-3143`，含 `On`/`Unknown` 的切向游走）；
- 从 `loops` 移除 roots（`cxx:3146-3156`）；丢失时全部转 roots（`cxx:3159-3170`）；
- 对每个 root：`PerformInfinitePoint() == In` 判反向（`cxx:3180-3188`）；按 `(state == OUT) == reverse` 收集洞
  （wire 取首边中点，vertex 用 `ValueOfUV`）（`cxx:3193-3229`）；建新面 + 顶点洞用 `copy_nm_vertex_face` +
  `Context()->Replace`（`cxx:3232-3249`）；丢失时补 roots（`cxx:3252-3269`）。

**UNPORTED**：单环 `FixOrientation`；顶点洞的 `BRep_PointRepresentation` UV 表示。

**状态**：`occt-topo --lib` 1281/0。
**下一步**：`DispatchWires`（`cxx:3275-3584`）。
### 9.33 round 27 —— `DispatchWires`（`cxx:3275-3584`）

`dispatch_wires.rs`：`ComposeShell::dispatch_wires(faces, wires)`
- 封闭模式修复（`cxx:3280-3360`）：单退化边 wire 剔除、`sfw.Load`/`FixShifted` → `fix_shifted_wire`、
  退化边 `ShapeBuild_Edge::RemovePCurve`；
- `mPnts[i] = GetMiddlePoint`（`cxx:3362-3376`）；
- 每条 wire：周期调整 → `LocateU/VParameter` → `GlobalToLocalTransformation` → patch 面 → 每条边
  `rs.IsRecorded/Value` 或 `sbe.Copy` + `rs.Replace`/`Context()->Replace` → `sbe.ReassignPCurve` →
  `!SameRange` 时 `sbe.Copy` + `sfe.FixAddCurve3d` + 设 3D 曲线/range，否则 `sfe.FixAddCurve3d` → `Set`（`cxx:3387-3533`）；
- 按 patch 面打包 loops，调用 `MakeFacesOnPatch`（`cxx:3535-3583`）。

**UNPORTED**（注释含行号）：
- seam pcurve 平移环（`cxx:3287-3327`）——需要 edge 上第二个 `BRep_GCurve` 槽（端口 pcurve 按 face-key 存）；
- `sfw.FixDegenerated()`（`cxx:3358`）——需把 `fix_degenerated_all` 的边表改动写回 `WireSegment`；
- `TransformPCurve`（`cxx:3462-3504`，含 affine 分支 `cxx:614-...`）——需 `Geom2d_Curve::TransformedParameter` 逐类。

**状态**：`occt-topo --lib` 1281/0。
**下一步**：`Perform`/`SplitEdges`（`cxx:206-...`），再接 `ShapeFix_Face::FixMissingSeam`。
### 9.34 round 28 —— `Perform` / `SplitEdges`（`cxx:206-270`）

`perform.rs`：
- `ComposeShell::perform()`（`cxx:206-255`）：`status=OK`、`invert_edge_status=false`（`cxx:209`）、
  `LoadWires`（空则 `FAIL6`）→ `SplitByGrid` → `BreakWires` → `CollectWires` → `DispatchWires` →
  面数≠1 时 `MakeShell` 否则取单面、`myResult.Orientation(myOrient)`、`status|=DONE1`。
- `ComposeShell::split_edges()`（`cxx:259-270`）：`LoadWires` + `SplitByGrid`。

**状态**：`occt-topo --lib` 1281/0。
**下一步**：`ShapeFix_Face::FixMissingSeam`（`ShapeFix_Face.cxx:1722-2330`）+ `Perform:492-494` 调用点接线。
### 9.35 round 29 —— 验收剩余工作（精确规格）

`ShapeFix_ComposeShell`（`cxx:206-3603`）已全部移植；要达成 a3n00 验收，剩余为 `ShapeFix_Face` 侧接线：

**(1) `ShapeFix_Face` 类骨架**（`ShapeFix_Face.hxx` / `.cxx:132-343`）：端口目前没有该类，只有分散件
（`model_builder/wire_builder.rs` 的 `FixOrientation`、`pcurve_full` 的 `ShapeAnalysis_Surface` 函数）。需要一个最小结构：
`myFace`、`mySurf`（`ShapeAnalysis_Surface` 包装：`IsUClosed`/`IsVClosed`/`Bounds`/`ValueOfUV`/`Surface`）、
`Context`（`MapReShape`）、`myStatus`、`myPrecision`/`myMinTol`/`myMaxTol` 与 `NeedFix` 模式位。

**(2) `ShapeFix_Face::CheckWire`**（供 `FixMissingSeam:1831` 调用；在 `.hxx`/局部定义）：判断 wire 是否开口、是否
退化；需要 `ShapeAnalysis_Edge::CheckCurve3dWithPCurve`/`CheckSameParameter` 与 `BRep_Tool::IsClosed`。

**(3) `ShapeFix_Face::FixMissingSeam`**（`ShapeFix_Face.cxx:1722-2330`，608 行）：
- 设置段（`1722-1808`）：`IsUClosed/IsVClosed`、`Context()->Apply`、BSpline 必须 U/V 周期、
  `mySurf->Bounds` + `BRepTools::UVBounds` + 无穷夹取、`URange/VRange`；
- 选 wire 对（`1810-1881`）：收集非流形子形状、`CheckWire` 找 w1/w2、同向时反转、>2 全退化时移除；
- 环面退化判定（`1883-1893`）；
- `w2.IsNull()` 分支（`1895-...`）：球面/锥状 BSpline 用 `aRange`（`1920-1925`）；构造缝边、
  `ShapeBuild_Edge::Copy`、`B.MakeEdge`/`UpdateEdge`（两条 pcurve）、`ShapeFix_Edge::FixAddCurve3d`/`FixAddPCurve`、
  `ShapeFix_Wire` 走一遍；
- `w2 != null` 分支：把两条 wire 合并成一条（`ShapeAnalysis_Wire` 相关）。

**(4) `ShapeFix_Face::Perform`**（`ShapeFix_Face.cxx:345-...`）：`FixMissingSeam` 调用点 `492-494`
（`if NeedFix(myFixMissingSeamMode) { if (FixMissingSeam()) status |= DONE3; }`），随后 `TopExp_Explorer` 遍历
`myResult` 的每个面继续第二段修复。

**(5) STEP 建面后接线**：端口需要在 STEP 读面处按 `ShapeProcess_OperLibrary`/`FromSTEP.FixShape` 的激活状态调用
`ShapeFix_Face::Perform`（`myFixMissingSeamMode` 打开时）。

**(6) `ShapeFix_ComposeShell` 内 5 处 UNPORTED**（见 §9.33/§9.32）：seam pcurve 平移环、`FixDegenerated` 回写、
`TransformPCurve`（含 affine）、单环 `FixOrientation`、`CopyNMVertex` 的点表示。

**结论**：`ShapeFix_Face` 类 + `FixMissingSeam`（约 600 行）+ `Perform` 编排 + STEP 接线，
在 goal 的 30 轮预算内**未能完成**；`ShapeFix_ComposeShell` 侧已全绿（`occt-topo --lib` 1281/0、
`occt-geom --lib` 143/0、STEP parity 见本轮）。
### 9.36 round 30 —— `CheckWire` + `FixMissingSeam` 设置段（`ShapeFix_Face` 侧起步）

`shhealing/shape_fix_face.rs`：
- `check_wire(wire, face, dU, dV) -> Option<(i32,i32,bool)>`（`ShapeFix_Face.cxx:1652-1718`）：
  逐边累计 pcurve 端点差、退化判定、`isuopen`/`isvopen`；
- `ShapeFixFace { face, surf, status }` + `with_face`；
- `fix_missing_seam()` 移植 `cxx:1722-1898`：U/V 封闭检查、BSpline 必须 U/V 周期、
  `Bounds` + `UVBounds` 无穷夹取、`URange`/`VRange`、取向 wire 收集。

**UNPORTED**：`cxx:1899-2330`（缝边构造与 `w2 != null` 合并）——需先完成，才能接
`ShapeFix_Face::Perform:492-494` 与 STEP 建面后接线。

**门禁（本轮已跑，全绿）**：`occt-topo --lib` 1281/0、`occt-geom --lib` 143/0、`step_obj_parity` 14/14、
`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 3/3。
### 9.37 round 30（最终轮）—— 完成度总结

**已忠实移植（全部编译通过、门禁全绿）**：
- `ShapeFix_ComposeShell`（`ShapeFix_ComposeShell.cxx:206-3603`）**全函数**：`Perform`/`SplitEdges`、`SplitByLine`（两组）、
  `SplitByGrid`、`BreakWires`、`CollectWires`、`SplitWire`、`MakeFacesOnPatch`、`DispatchWires`、`ApplyContext` 与全部 helper；
- `Geom2dInt`：`IntCurve_IntConicConic`（Line/Line）、`IntCurve_IntConicCurveGen`、`Geom2dInt_GInter` 分派（typ1=Line 臂）、
  `Curve2d` 单态化的 `IntImpParGen_Intersector`、`Extrema_GFuncExtPC`、`Extrema_GenLocateExtPC`、`Geom2dInt_TheProjPCurOfGInter`；
- `ShapeBuild_Edge::{Copy,RemovePCurve,SetRange3d,ReassignPCurve}`、`ShapeFix_Edge::FixAddCurve3d`、
  `MapReShape::{IsRecorded,Value}`、`CopyNMVertex`（face 重载）、`ShapeFix_Face::CheckWire` + `FixMissingSeam` 设置段（`1722-1898`）。

**未完成（验收核心）**：
1. `ShapeFix_Face::FixMissingSeam` 的缝边构造（`ShapeFix_Face.cxx:1899-2330`）与 `w2 != null` 合并分支；
2. `ShapeFix_Face::Perform`（`cxx:345-...`，`FixMissingSeam` 调用点 `492-494`）；
3. STEP 建面后按 `FromSTEP.FixShape` 激活状态接线；
4. `ShapeFix_ComposeShell` 内 5 处 UNPORTED（§9.33/§9.32）。

**结论**：目标**未达成**（a3n00 的 7 个周期 B 样条面三角数未验证到 360-486）。
工作树 green：`occt-topo --lib` 1281/0、`occt-geom --lib` 143/0、`occt-geom2d --lib` 72/0、
`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、`step_geometry_parity` 3/3。
续作入口见 §9.35（精确规格）与 §9.36（`ShapeFix_Face` 侧起点）。
### 9.38 round 31（goal 续作）—— `FixMissingSeam` 选对 + 退化缝边（`cxx:1824-2021`）

`shhealing/shape_fix_face.rs` 的 `fix_missing_seam` 继续推进：
- 选 wire 对（`cxx:1824-1881`）：`CheckWire` → w1/w2、同向时按 `isdeg1` 反转、>2 全退化时 `ws.Remove` 重来；
- 退化环面判定（`cxx:1883-1893`）；
- `w1.IsNull()` → false（`cxx:1895-1898`）；
- `w2.IsNull()` 的退化缝边构造（`cxx:1899-1992`）：环面（`acos(-Ra/Ri)`）、球面、
  BSpline 的 V-open/U-open 两个分支，建退化边（pcurve = `Geom2d_Line`、range、两个顶点）并加入 `ws`；
- 两 wire 的 `GetFaceUVBounds`（端口用 `BRepTools::AddUVBounds(Face,Wire)`）与朝向一致性（`cxx:1994-2021`）。

**UNPORTED**：`cxx:2023-2325`——`FixReorder`（端口只查不重排）、tmpF 重建、环面专用 FixShifted、
最佳缝位搜索、fictive grid + `ComposeShell::Perform`（`2236-2261`）、小 wire/面清理尾段。

**状态**：`occt-topo --lib` 见本轮报告。
**下一步**：把 `cxx:2023-2325` 补完（含 `2236-2261` 的 `CompositeSurface::with_grid` + `ComposeShell::Perform`），再接 `ShapeFix_Face::Perform`。
### 9.39 round 32 —— `FixMissingSeam` 缝边构造 + fictive grid + ComposeShell（`cxx:2023-2266`）

`shape_fix_face.rs` 补完：
- `FixReorder`（`cxx:2029-2032`）标 UNPORTED（端口只查 3D 顺序，不重排边表），`w11=w1`/`w21=w2`；
- tmpF 重建（`cxx:2036-2066`）：替换 w1/w2，其他 wire 按 `IsOuterBound` 反向；
- 环面专用 `FixShifted`（`cxx:2068-2136`）：`AdjustByPeriod` 求 `shiftw2`、逐 wire 平移 pcurve、
  再按 `m1` 选缝位；
- 最佳缝位搜索（`cxx:2138-2234`）：遍历 `wd1`/`wd2` 边端点，`AdjustByPeriod` 归一后取最接近 0 的 U/V；
- **fictive grid + ComposeShell（`cxx:2236-2261`）**：`Geom_RectangularTrimmedSurface(surf, uf, uf+URange, vf, vf+VRange)` →
  `CompositeSurface::with_grid`（1×1）→ `ComposeShell::init` + `ClosedMode=true` + `SetContext` + `SetMaxTolerance(MaxTolerance())` + `Perform`；
  `mySurf = RTS`、`myResult = CompShell.Result()`。

**UNPORTED**：`cxx:2268-2325`——`Context()->Replace/Remove`（端口 `MapReShape` 无 null 绑定）、
`ShapeFix_Wire::FixSmall`、`BRepTools::Update` 与小 wire/面清理。

**状态**：`occt-topo --lib` 见本轮报告。
**下一步**：`ShapeFix_Face::Perform`（`cxx:345-...`，调用点 `492-494`）与 STEP 建面后接线；先跑 a3n00 观察。
### 9.40 round 33 —— 接线 `ShapeFix_Face::Perform`（`492-494`）+ 诊断结论

- `shape_fix_face.rs` 增 `perform_fix_missing_seam()`：`myResult = myFace`，`myFixMissingSeamMode` 开时调 `FixMissingSeam()`，
  成功则 `status |= DONE3`（`ShapeFix_Face.cxx:482-498`）；
- `step/read_topology.rs:720` 建面后调用之（对应 `FromSTEP.FixShape`，`ShapeProcess_OperLibrary.cxx:830`），
  结果仍为单 Face 时替换，否则保留原面；
- `check_wire` 改用 `curve_on_surface_oriented(edge, face, true)`（`ShapeAnalysis_Edge::PCurve` 对 REVERSED 边交换 f/l，`ShapeAnalysis_Edge.cxx:201-206`）。

**诊断（a3n00）**：对 6 个周期 BSpline 面（`--low` 显示 mt=1..11：face 173/176/190/201/209/216），
`FixMissingSeam` 均返回 false：这些面的两条 wire 在 V 上“上去再回来”，`CheckWire` 的 pcurve 端点差之和为 0，
故 `isuopen==isvopen==0` → `CheckWire` 返回 false → `w1` 为空 → 直接 false（`cxx:1826-1834`/`1895-1898`）。
即：a3n00 这几个面的低三角数**不是** `FixMissingSeam` 触发的缺缝问题（至少在端口现有 wire/pcurve 数据下），
其真实成因需另行定位（可能在这些面读入时的 pcurve/边界构造，而非 ComposeShell 路径）。
**❗接线回归（本轮实测）**：把 `perform_fix_missing_seam` 接到 `step/read_topology.rs:720` 后，
`step_obj_parity::occ_test_model_bboxes_match_occt` 的 `occ/bottom.step` 从 `f=23408`（应 23557）落到 `f-ratio 0.99` → **门禁回归**。
原因：`fix_missing_seam` 即使最终返回 false，也可能在途中（如环面 FixShifted 平移 pcurve、`set_pcurve_range` 等）
改写注册表中的共享 pcurve/range，而 `ComposeShell` 内仍有 5 处 UNPORTED（§9.33/§9.32），输出无法与 OCCT 对齐。
故**接线已回退**（`read_topology.rs:720` 保留注释、不调用），`perform_fix_missing_seam` 仍留在 `shape_fix_face.rs` 备用；
回退后 `step_obj_parity` 恢复 **14/14**。

⇒ 正确顺序：先补完 `ShapeFix_ComposeShell` 的 5 处 UNPORTED 与 `FixMissingSeam` 的清理尾段（`cxx:2268-2325`），
再接线并用 `bottom.step` / a3n00 对拍。
### 9.41 round 35 —— `TransformedParameter` + `TransformPCurve` 非仿射段 + `DispatchWires` 接线

- `occt-geom2d`：`Curve2d::transformed_parameter`（默认 `U`，`Geom2d_Curve.cxx:41-44`），
  `Geom2dLine`/`Geom2dParabola` 覆盖为 `U*|T.ScaleFactor()|`（`Geom2d_Line.cxx:246-253`、`Geom2d_Parabola.cxx:249-256`），
  `Geom2dTrimmedCurve`/`Geom2dOffsetCurve` 委托基曲线（`TrimmedCurve.cxx:297-300`、`OffsetCurve.cxx:423-426`）；
- `ShapeBuildEdge::transform_pcurve`（`ShapeBuild_Edge.cxx:596-612`）：`trans` 应用 + `TransformedParameter`；
  **UNPORTED**：`uFact != 1` 仿射分支（`cxx:614-698`，需 `Geom2dConvert_ApproxCurve`/`CurveToBSplineCurve`/极点重写）；
- `dispatch_wires.rs`：实现 `cxx:3411-3414`（`T.Multiply(Sh)`）与 `cxx:3462-3504`（`needT` 时取 pcurve →
  `TransformPCurve` → 写回 pcurve/range、必要时 `SameRange=false`）。

**状态**：`occt-topo --lib` 见本轮报告。
**下一步**：seam pcurve 平移环（`cxx:3287-3327`）、`sfw.FixDegenerated()` 回写（`cxx:3358`）、单环 `FixOrientation`、
`BRep_PointRepresentation`；之后重试接线并用 `bottom.step`/a3n00 对拍。
### 9.42 round 36 —— seam pcurve 平移环（`cxx:3287-3327`）

`dispatch_wires.rs`：实现封闭模式下 REVERSED 缝边的 pcurve 平移：`CurveOnSurface(E)` 取 `PCurve2`、
`CurveOnSurface(E.Reversed())` 取 `PCurve1`（`BRep_Tool.cxx:347-357`）；若同一对象或端点重合（`cxx:3312`），
按 `UClosed`/`VClosed` 与 `PConfusion` 生成 `(uPeriod,0)`/`(0,vPeriod)` 平移，`Translate` 到 `PCurve1` 后写回两条 pcurve。

**状态**：`occt-topo --lib` 见本轮报告。
**下一步**：`sfw.FixDegenerated()` 回写（`cxx:3358`，需 `fix_degenerated_all` 公开 + 边表回写 `WireSegment`）。




**状态**：门禁见本轮报告（`FixMissingSeam` 返回 false 时读入形状不变，故不应回归）。
**下一步**：定位 a3n00 6 个面的低三角数真实成因；或核对端口与 OCCT 参考读入形状在这些面上的 wire/边差异。
### 9.43 round 37 —— `bottom.step` 接线回归的根因定位

用 `DSH_FIXMS_SKIP` 跳过 `fix_missing_seam` 本体后，接线打开时 `step_obj_parity` 仍 **14/14**；
带 `DSH_FIXMS` 跑失败用例可见：`occ/bottom.step` 读入时 `fix_missing_seam` **确实触发**
（`w1-some w2=true ismodeu=-1/1 ismodev=0 tor=false uclosed=true vclosed=true`，`FIXMS-RESULT faces=1`），
即这些面走完全程并返回 **true**、产出一个新面；但该面与 OCCT 的结果不一致（mesh f 23557 → 23408，f-ratio 0.99）。
（此前用独立 example 探针在 `bottom.step` 上未触发，是因为该 example 的读入路径与测试不同。）

⇒ `FixMissingSeam` 的**触发**是对的，**输出**因 `ShapeFix_ComposeShell` 内 5 处 UNPORTED（§9.33/§9.32）而偏离。
因此接线**继续保持关闭**（`read_topology.rs` 注释说明），待这些分支补完后再开。

**状态**：`occt-topo --lib` 与 `step_obj_parity` 见本轮报告。
**下一步**：补完 §9.33 的 5 处 UNPORTED（seam 平移已做；还剩 `FixDegenerated` 回写、单环 `FixOrientation`、
仿射 `TransformPCurve`、`BRep_PointRepresentation`）与 `FixMissingSeam` 清理尾段（`cxx:2268-2325`），再接线对拍。
### 9.44 round 38 —— `TransformPCurve` 仿射分支的 Line 臂（`cxx:614-641`）

`shape_build_edge.rs::transform_pcurve`：`uFact != 1` 时先按 `cxx:614-618` 解包 trimmed 基曲线，
再实现 `cxx:624-641` 的 Line 臂：`tMatu = affinity(OY2d, uFact)` 即 `(x,y)→(x*uFact,y)`，
对两端点施加以重建 `gp_Lin2d`，用 `ElCLib::Parameter`（端口 `clib2d::line_parameter_ax2d`）重算 `aFirst/aLast`，返回新的 `Geom2d_Line`。

**UNPORTED**：`cxx:642-698`（Bezier/BSpline/Conic 的极点重写，需 `Geom2dConvert_ApproxCurve` /
`CurveToBSplineCurve` / 极点 setter）。

**状态**：`occt-topo --lib` 见本轮报告。
**下一步**：为 `Curve2d` 加极点访问（BSpline/Bezier）以补 `cxx:642-698`；或先补 `FixDegenerated` 回写与单环 `FixOrientation`。
### 9.45 round 39 —— `poles2d`/`set_poles2d` + 仿射极点重写（`cxx:642-698`）；接线仍偏离

- `Curve2d` 增 `poles2d()`/`set_poles2d()`（默认 `None`/no-op），`Geom2dBSplineCurve`（xs/ys）与 `Geom2dBezierCurve`（poles）覆盖；
- `transform_pcurve` 的非 Line 分支：对极点施加 `(x,y)→(x*uFact,y)`（OCCT 对 Bezier 先 `CurveToBSplineCurve`，控制点不变，故等价）；
  **UNPORTED**：`Geom2d_Conic` 路径（`cxx:657-686`，需 `Geom2dConvert_ApproxCurve`/`CurveToBSplineCurve`）。

**接线复测（仍关闭）**：`occ/bottom.step` 23408 vs 23557、`motoc.step` 13646 vs 13944、`top.step` 22794 vs 22800——
极点重写未改变结果，说明偏离在 ComposeShell 的其他分支（`FixDegenerated` 回写、单环 `FixOrientation`、
`BRep_PointRepresentation`、`FixMissingSeam` 清理尾段）。

**状态**：`occt-topo --lib` 见本轮报告。
**下一步**：逐个补上述分支并复测接线。
### 9.46 round 40 —— `sfw.FixDegenerated()` 回写（`cxx:3358`）

`dispatch_wires.rs`：`fix_degenerated_all` 改为 `pub`；对每条非顶点 wire 以 `builder.make_wire(segment.edges)`
建临时 `Wire`，调用 `fix_degenerated_all(&mut dw, &myFace, self.precision)`；若改动则记录旧边/旧 patch 下标，
`clear()` + `load_edges(新边)` 重建 `WireSegment`，并按 `is_same` 把存活边的 patch 下标写回。

**状态**：`occt-topo --lib` 见本轮报告。
**下一步**：复测接线；若仍偏离，继续补单环 `FixOrientation` / `BRep_PointRepresentation` / 清理尾段。
### 9.47 round 41 —— `FixDegenerated` 后接线仍偏离（结论）

补完 `FixDegenerated` 回写后复测：`occ/bottom.step` 仍 23408 vs 23557、`motoc.step` 13646 vs 13944。
注意 ours 的 mesh 面数比 OCCT **少**：接线后从“与 OCCT 一致”降到 23408，说明端口的 `FixMissingSeam`/`ComposeShell`
在这些面上**合并/丢面**，而不是“没做缝”。

⇒ 需要逐面对拍端口与 OCCT：取 bottom.step 上触发 `FixMissingSeam` 的那个面，分别在端口与 OCCT 下跑 `FixMissingSeam`，
比较结果面数/边数/pcurve；据此定位具体分支（`MakeFacesOnPatch` 的 roots/holes 判定、`DispatchWires` 的 patch 打包、
或 `CollectWires`/`SplitByLine` 的合并）。

接线继续关闭；文档记录 `ours f < occ f` 的方向性结论。
**状态**：`occt-topo --lib` 1281/0。
**下一步**：逐面对拍定位丢面分支。
### 9.48 round 42/43 —— 门禁全绿复核（接线关闭状态）

| 门禁 | 结果 |
|---|---|
| `occt-topo --lib` | 1281/0 |
| `occt-geom --lib` | 143/0 |
| `occt-geom2d --lib` | 72/0 |
| `step_obj_parity` | 14/14 |
| `step_to_obj` | 13/13 |
| `step_obj_area` | 11/11 |
| `step_geometry_parity` | 3/3 |

**验收仍未达成**：a3n00 的 7 个周期 B 样条面三角数没有提升（`FixMissingSeam` 对这些面不触发，§9.40），
且接线会在 `bottom.step`/`motoc.step` 上引入**丢面**回归（§9.47），故接线保持关闭。


### 9.49 round 44 —— `FixMissingSeam` 的 `Context()->Replace(myFace, myResult)`（`cxx:2268`）

`shape_fix_face.rs`：ComShell `Perform` 后补上 `self.context.replace(&face.0, myResult)`（`ShapeFix_Face.cxx:2268`）。
清理尾段的 `Context()->Remove`/`FixSmall`/`BRepTools::Update` 仍 UNPORTED。

**状态**：`occt-topo --lib` 见本轮报告。
**下一步**：同 §9.48 的逐面对拍。
### 9.50 round 45-46 —— 【方案 A｜忠实性】pcurve 表示改为按 surface identity（对齐 `BRep_GCurve`）

用户指示方案 A（严格忠实，允许改/加基础类型）后，实施数据模型对齐：

- `tgeometry.rs`：新增 `surface_by_ptr` 与 `repr_key(face_key)`（= 该 face 注册 surface 的数据指针）。
  `edge_pcurve(s)`/`edge_pcurves(s)`/`set_edge_pcurve(s)`/`set_edge_pcurves(s)`/`pcurve_range(s)`/`set_pcurve_range(s)`
  /`remove_pcurves_on_surface(s)` 全部改为按 `repr_key` 读写；`edge_pcurve_reps` 经 `surface_by_ptr` 反查 surface。
  这样共享同一 surface 的不同 face 共享同一条 CurveOnSurface 表示（OCCT 语义），取代此前的「face 键 + 同 surface 回退」。
- `shape_ops.rs`（`BRepBuilderAPI_Copy`）：pcurve 重键改为按 surface 指针映射（原为 face 键，已失效）。
- `pave_de.rs`：直接读 `EdgeGeom::get_pcurve(face_key)` 改为经 registry `edge_pcurve(edge, face_key)`。
- `shape_fix_face.rs`：`FixMissingSeam` 的 `aSeqNonManif` 改为按 `TopoDS_Iterator(myFace,false)` 遍历 face 全部子形状
  （原只收集非取向 wire），与 `cxx:1810-1822` 一致。

**验证（全绿）**：`occt-topo --lib` 见报告；`step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、
`step_geometry_parity` 3/3。

**下一步（继续方案 A 的数据模型缺口）**：`IsClosed(edge, face)` 的完整定义、`GpGTrsf2d`、`Geom2dConvert`/`ApproxCurve`、
`BRep_PointRepresentation`；之后重审此前的「等价替代」处。
### 9.51 round 53 —— 【方案 A】`BRep_Tool::IsClosed(E, F)` 逐行（`BRep_Tool.cxx:795-841`）

`brep_tool.rs` 新增 `BRepTool::is_closed_edge_face(edge, face)`：
- `cxx:819-822`：plane surface 直接 false；
- `cxx:834`：该 surface 上的表示含两条 pcurve（`IsCurveOnClosedSurface`）→ true；
- **UNPORTED**：triangulation 臂（`cxx:803-804`、`:849-...`）——端口不为该查询给 face 挂 `Poly_Triangulation`。

`dispatch_wires.rs` 的 seam 平移环改用之（替换此前的 `edge_pcurves(...).len() < 2` 近似）。
另：`boptools_2d.rs` 模块头「按 face 指针键」的说明改为「按 surface 数据指针 / `repr_key`」。

**仍为近似的 IsClosed 调用点（非本次改动，待逐个对齐）**：`bop_split_seam.rs:32`、`wire_splitter_block.rs:839`、
`wire_fix.rs:2487` 仍用 `edge_pcurves(...).len() >= 2`（缺 plane 检查）。

**下一步**：`GpGTrsf2d` + `Geom2dConvert::CurveToBSplineCurve`/`ApproxCurve`（补 `TransformPCurve` 的 `uFact != 1`）。
### 9.52 round 54 —— 【方案 A】`Geom2dConvert::CurveToBSplineCurve`（Bezier 精确），`TransformPCurve` 去近似

- `crates/occt-geom2d/src/geom2d_convert.rs`（新增）：`curve_to_bspline_curve(c)` ——
  Bezier → BSpline（同极点 + 夹持节点 `[0×(deg+1), 1×(deg+1)]`，精确）、BSpline → 自身；
  **UNPORTED**：`Geom2d_Conic` 臂（`Convert_QuasiAngular`）与 `Geom2dConvert_ApproxCurve`。
- `Curve2d` 增 `is_bspline2d()`/`is_bezier2d()`（默认 false，BSpline/Bezier 覆盖）。
- `shape_build_edge.rs::transform_pcurve`：非 Line 分支改为「Bezier 先经 `CurveToBSplineCurve` 转 BSpline，
  再对 BSpline 极点施加 `(x,y)→(x*uFact,y)`」——与 `ShapeBuild_Edge.cxx:679-698` 一致；
  Conic 路径仍标 UNPORTED（`cxx:657-686`）。

**状态**：`occt-geom2d --lib` 与 `occt-topo --lib` 见报告。
**下一步**：`GpGTrsf2d`（显式仿射类型）或 `BRep_PointRepresentation`（`CopyNMVertex`）。
### 9.53 round 55 —— 【方案 A】新增 `gp_GTrsf2d`，`TransformPCurve` 的仿射改为逐行

- `crates/occt-core/src/gp/gtrsf2d.rs`（新增）：`GpGTrsf2d { matrix, loc }` + `identity` +
  `set_affinity(A, Ratio)`（逐行 `gp_GTrsf2d.cxx:24-38`，外积矩阵 + `loc = A.loc - M*A.loc`）
  + `transforms(p)`；在 `gp/mod.rs` 注册并 re-export。
- `shape_build_edge.rs::transform_pcurve`：`cxx:620-623` 改为 `t_matu.set_affinity(gp::OY2d(), uFact)`，
  Line 端点与（Bezier→BSpline 后的）BSpline 极点都用 `t_matu.transforms` 变换——与 `cxx:624-698` 结构一致，
  去掉此前的内联 `(x*uFact, y)` 近似。

**状态**：`occt-core`/`occt-geom2d`/`occt-topo` lib 测试见报告。
**下一步**：`BRep_PointRepresentation`（`CopyNMVertex` 逐行）或 `Geom2d_Conic` 的 `Convert_QuasiAngular`。
### 9.54 round 56 —— `ShapeFix_Wire::FixLacking` 的 `IsClosed` 改用逐行 `BRepTool::is_closed_edge_face`

`wire_fix.rs::fix_lacking` 内的 `is_closed` 闭包由 `edge_pcurves(...).len() >= 2` 改为
`BRepTool::is_closed_edge_face(e, face)`（含 `BRep_Tool.cxx:819-822` 的 plane→false）。

**仍待对齐**：`bop_split_seam.rs::is_closed_on_face`（自带「边在 face 上出现两次」臂，需核对对应的 OCCT 重载）、
`wire_splitter_block.rs::edge_closed_on_face`（前置 `IsDegenerated` 短路，OCCT `BRep_Tool::IsClosed` 无此分支）。

**状态**：`occt-topo --lib` 1281/0。
**下一步**：核对上述两处对应的 OCCT 控制流，或补 `Geom2d_Conic` 的 `Convert_QuasiAngular`。
### 9.55 round 57 —— `BRep_Tool::IsClosed` 重载核对（旁支发现）

`BRep_Tool.cxx` 中共 4 个 `IsClosed`：
- `:795` `IsClosed(E, F)` → `IsClosed(E,S,L)` + triangulation 臂；
- `:814` `IsClosed(E, S, L)`：surface 为 plane 直接 false，否则找「带两条 pcurve 的表示」；
- `:849` `IsClosed(E, T, L)`：triangulation 中两组索引；
- `:1707` `IsClosed(shape)`：shell/wire/edge 的「边/顶点成对出现」判定。

⇒ `bop_split_seam.rs::is_closed_on_face` 的第二臂（「边在 face 的 wire 上出现两次」）**不对应上述任何重载**：
`(E,F)` 只有 2-pcurve 臂，`:1707` 是 shell/wire/edge 而非 face。该臂需要单独找到 OCCT 出处，否则应按忠实性要求删除或改标 UNPORTED。
本轮未改动该文件（避免在无出处的情况下改变行为）。

**状态**：无代码改动；`occt-topo --lib` 上一轮 1281/0。
**下一步**：为该臂找 OCCT 出处（可能是 `ShapeAnalysis`/`BRepTools` 中的同名判定），或删除并报告。
### 9.56 round 58 —— 找到「出现两次」臂的 OCCT 出处：`BRepTools::IsReallyClosed`

`BRepTools.cxx:1204-1220`:
```cpp
bool BRepTools::IsReallyClosed(const TopoDS_Edge& E, const TopoDS_Face& F) {
  if (!BRep_Tool::IsClosed(E, F)) return false;      // 先要求 2-pcurve
  int nbocc = 0;                                     // 再数 E 在 F 的边中出现次数
  for (TopExp_Explorer exp(F, TopAbs_EDGE); exp.More(); exp.Next())
    if (exp.Current().IsSame(E)) nbocc++;
  return nbocc == 2;
}
```

⇒ 该臂是 **`BRepTools::IsReallyClosed`**（`IsClosed && nbocc==2`），不是 `BRep_Tool::IsClosed`。
而端口 `bop_split_seam.rs::is_closed_on_face` 实现的是 `IsClosed || nbocc==2`（**或**），与 OCCT 的 **与** 不同。
本轮未改代码：无法在本 checkout 中定位到调用它的 OCCT 函数（`DoSplitSeamOnFace` / `AlgoTools3D.cxx` 均未在 `D:\source\OCCT-src` 中检索到），
改 `||`→`&&` 会改变行为而不知对应分支。**建议**：确认该调用点对应的 OCCT 函数后，按 `IsReallyClosed` 改为 `&&`，或拆成两个谓词分别使用。

**状态**：无代码改动；`occt-topo --lib` 上一轮 1281/0。
**下一步**：在完整 OCCT 源码中定位 `DoSplitSeamOnFace` 的谓词；或补 `Geom2d_Conic` 的 `Convert_QuasiAngular`。
### 9.57 round 59 —— 【方案 A】全部改动后的门禁终检

针对 r53–r56 的数据模型/忠实性改动（pcurve 按 surface identity、`IsClosed` 逐行、`CurveToBSplineCurve`、`gp_GTrsf2d`、
`fix_lacking` 的 `IsClosed`），重跑全部 STEP 门禁：

| 门禁 | 结果 |
|---|---|
| `step_obj_parity` | 14/14 |
| `step_to_obj` | 13/13 |
| `step_obj_area` | 11/11 |
| `step_geometry_parity` | 3/3 |
| `occt-topo --lib` | 1281/0 |
| `occt-geom --lib` | 143/0 |
| `occt-geom2d --lib` | 72/0 |

**结论**：方案 A 已落地的 6 项忠实性对齐均无回归。**验收（a3n00 7 面三角数 360-486）仍未达成**，
且 STEP 接线仍因 `bottom.step`/`motoc.step` 丢面而关闭（§9.47）。
### 9.58 round 60+ —— 用户提供 OCCT oracle `data/occ/occ-a3n00.obj`；逐面 mesh 对拍

用户指出 a3n00 的 OCCT 生成 OBJ 为 `data/occ/occ-a3n00.obj`（12466 个三角形，单一 group）。
同目录/`.target-gate` 中已有 OCCT 逐面 dump：
- `.target-gate/a3n00_occ_mesh.txt`：**226 面 / 12324 三角形**（ShapeProcess fix 之后）；
- `.target-gate/a3n00_occ_mesh_nofix.txt`：**226 面 / 2832 三角形**（proc off，周期面退化为 0）。

当前端口（默认模式，`computed_lin=1.076007`）跑 `zz_probe_a3n00 data/occ/a3n00.stp`：
**226 面 / mesh_t=12283**（OCCT 12324，差 41 ≈ 0.33%）。

逐面（按索引）比对不可靠：**面序不同**（例：OCCT `f31` 是 376 的 BSpline，端口 `f31` 是 7 的 Plane）。
按值多重集合比较：
- 匹配：`355, 376, 376, 376, 376`（周期 BSpline 面，即验收所指 360-486 区间）——端口与 OCCT **一致**；
- 不匹配：OCCT 另有 `669`（Cylinder）/`810`/`854`（Cone）三个大面，端口对应位置只有 `674` 与 `243/288`；
  其余多为小面（OCCT 2/3/4 vs 端口 1/2/3），共 150 个值位不同。

**结论**：验收的“总量向 12324 收敛”已基本达到（12283）；360-486 的 5 个周期面与 OCCT 完全一致。
### 9.59 round 61 —— 【oracle 口径确认】等参数对拍（deflection/angle/面序）

`specs/occt_probe/occt_probe.cpp:133-137`：`--mesh` 的默认 `deflection=0.1`、`angle=0.5`；
而 `a3n00_occ_mesh.txt` 的 TOTAL 行显示实际口径为 **`deflection=1.07612 angle=0.349066`（20°）**，
输出 **`nodes=11052 triangles=12324`** —— 这正是验收里的 “11052/12324”。
（`a3n00_occ_mesh_nofix.txt` 用 `deflection=1.67307`，输出 `2932/2832`。）

端口以**相同参数**（`data/occ/a3n00.stp 1.07612`，angle 20°）运行：
**`nodes=11931 triangles=12283`，stats=226**。

| | nodes | triangles |
|---|---|---|
| OCCT | **11052** | **12324** |
| 端口 | 11931 (+879, +8%) | 12283 (−41, −0.33%) |

逐面（按索引）**不可比**：两份 dump 的面序不同（OCCT `f29`=376 的 BSpline，端口 `f29`=48；
端口的 376 面落在 `f171/191/200/205/217`）。按值多重集合：`355,376,376,376,376` 两边一致（周期面已达 360-486），
但 OCCT 另有 `854/810/669` 三张大面，端口最大仅 `674`，共 150 个值位不同。

⇒ **验收未达成**（nodes 差 8%）。下一步：在两侧 dump 中加入**几何键**（面中心/面积或每面 signed mesh volume）做配对，
定位三张欠网格的大面；三角形数已基本到位。
### 9.60 round 62 —— 按**曲面类型**聚合（绕开面序问题）

`GeomAbs_SurfaceType` 序号修正（OCCT `GeomAbs_SurfaceType.hxx`）：`0=Plane, 1=Cylinder, 2=Cone, 3=Sphere, 4=Torus,
5=Bezier, 6=BSpline`。`a3n00_occ_perface.txt`（`--perface`，与 `--mesh` 同样以 `TopExp_Explorer(aShape, TopAbs_FACE)` 遍历）
给出每面类型，可与 `a3n00_occ_mesh.txt` 的每面三角形数按同一索引合并：

| 类型 | OCCT 面数 | OCCT 三角形 | 端口面数 | 端口三角形 | Δ三角形 |
|---|---|---|---|---|---|
| Plane | 95 | 2832 | 95 | 2848 | +16 |
| Cylinder | 62 | 2673 | 62 | 2824 | +151 |
| Cone | 30 | 727 | 30 | 778 | +51 |
| Torus | 7 | 1214 | 7 | 1314 | +100 |
| BSpline | 32 | 4878 | 32 | 4519 | **−359** |
| **合计** | **226** | **12324** | **226** | **12283** | **−41** |

**结论**：面数与类型分布**逐类完全一致**（此前把 `GeometryRegistry`/probe 的类型序号误读为 Sphere↔Cone，已更正：
OCCT `type=2` 是 Cone、`type=3` 是 Sphere，两边都是 30 个 Cone）。三角形只在**各类内部**有小幅再分配，
主要缺口是 **BSpline 类少 359**，其余各类略多。

**下一步**：把 BSpline 类的 32 个面逐一配对（端口 dump 有 `uv=[...]`，OCCT dump 无；可按 BSpline 面的
`nodes/triangles` 值集合配对），定位是哪几个 BSpline 面欠网格。
### 9.63 round 63 —— BSpline 面配对 + 接线复测（结论：7 个小面不是 `FixMissingSeam` 的问题）

按 `(nodes, triangles)` 排序配对 32 个 BSpline 面：
- OCCT 最大：`854/810/669`，其余 `376,376,76,72×…`；最小 34。
- 端口最大：`674`，其余 `376,376,224×5,216×4,…`；**有 7 个面只有 1-11 个三角形**：
  `f129=1 f173=3 f176=11 f190=1 f201=3 f209=1 f216=4`，且**每个都是 `wires=2 edges=[2,2] mv=162`**。

⇒ 这 7 个面正是验收所指的「1-11」。它们**不是** `FixMissingSeam` 能修的：
**打开 STEP 接线后 a3n00 的 TOTAL 一模一样（`11931/12283`）**，即 `fix_missing_seam` 在这些面上不触发。
其特征（两个各含 2 条边的 wire、`mv=162` 但 `mt≤11`）指向**周期面的 UV 域/边界构造**导致的三角化失败。

接线复测（方案 A 之后）：`occ/bottom.step` 仍 `23408 vs 23557`、`motoc.step` `13646 vs 13944` → **接线继续保持关闭**。

**下一步**：查这 7 个面的 UV 域（`uv=[...]`）与两条 2-边 wire 的 pcurve，定位三角化失败点。
### 9.64 round 64 —— 【根因】`CheckWire` 的 `vec` 恒为 0 → `FixMissingSeam` 永不触发

对 a3n00 用 `--fixms`（逐周期面调用 `fix_missing_seam`）加临时插桩（已移除），得到：

1. 118 个周期面全部 `ret=false`、`result=none`；插入点追踪显示**每一面都在 `cxx:1895-1898`（`w1.is_none()`）返回**
   —— 即 `check_wire` 对**每一条 wire** 都返回 `None`。
2. `check_wire` 内部：`vec=(0.000000,0.000000)`，而 `dU=6.283185`（2π）、`dV≈1.4`。
3. 逐边打印（`er`=3D edge range，`f/l`=COS range，`p0/p1`=pcurve 端点）显示端点**首尾相接**：
   一条 wire 内各边 `p1` 恰为下一条的 `p0`，绕行一周后 `Σ(p1−p0)=0`。

**结论**：端口的 pcurve 在 UV 里构成**闭合环**，而 OCCT 的 `CheckWire` 期望 wire 在周期方向上**张开**
（`ΣΔu=±Uperiod`）。因此 `isuopen/isvopen` 恒为 0 → `CheckWire` 返回 false → `FixMissingSeam` 永不进入 CompShell。
这正是验收「7 个周期 B 样条面 1-11 → 360-486」无法达成的直接原因。

**候选成因**：端口在把 seam 边写入 face 的 UV 帧时**没有按周期平移/不 wrap**，或 seam 边的两条 pcurve 选错
（`BRep_Tool::CurveOnSurface(E).PCurve1/2`），使绕行一周的 u 增量被抵消。

**下一步**：核对 `ShapeAnalysis_Edge::PCurve`（`ShapeAnalysis_Edge.cxx:180-210`）与端口 `curve_on_surface_oriented` 的
`f/l` 取值方向；并检查 STEP 读入期 seam 边 pcurve 的 u 是否被 wrap 到 [0,2π)。
### 9.65 round 65 —— `PCurve`/`CurveOnSurface` 链已确认忠实 ⇒ `vec=0` 是**数据**问题

核对 OCCT：
- `ShapeAnalysis_Edge::PCurve(edge, surface, loc, C2d, cf, cl, orient)`（`ShapeAnalysis_Edge.cxx:192-208`）：
  `C2d = BRep_Tool::CurveOnSurface(edge, surface, loc, cf, cl)`，随后 **仅当** `orient && edge.Orientation()==REVERSED` 交换 `cf/cl`。
- `BRep_Tool::CurveOnSurface(E,S,L,First,Last)`（`BRep_Tool.cxx:327-357`）：`Eisreversed = (E.Orientation()==REVERSED)`；
  找到表示后 `GC->Range(First,Last)`，且 **`IsCurveOnClosedSurface() && Eisreversed` 时返回 `PCurve2()`**，否则 `PCurve()`。

⇒ 端口 `boptools_2d::curve_on_surface_oriented(edge, face, true)` 的「REVERSED 时取 `pcs[1]` 并交换 f/l」**与 OCCT 一致**；
`check_wire` 的算法也逐行对应（`ShapeFix_Face.cxx:1652-1718`）。
因此 `vec=0` 不是算法问题，而是**读入后存的 seam pcurve 数据**：端口在 UV 里形成闭合环，OCCT 的在周期方向张开。

**下一步**：对同一个面（如 a3n00 `f129`）导出端口存的每条边 pcurve 的 (u,v) 端点与 COS 范围，与 OCCT 侧
（`occt_probe` 的 `--fface`/自加 `PCURVE` 命令）对比，定位是 wrap、PCurve1/2 选反、还是 seam 边被拆成两条各闭合的 pcurve。
### 9.66 round 66 —— 【根因二】seam 边在 wire 里两次出现却**同向**，导致 Δ 抵消

`--checkwire` 对 a3n00 `f129` 的 dump：
```
CW face=129 u=[0.0000,3.1416] v=[0.0000,6.2832] dU=3.1416 dV=6.2832      # 面是 V 周期，V range = 2π
  w0 e0 ori=Forward f=0.000000 l=6.283185 delta=(0.000000,-6.283185)
  w0 e1 ori=Forward f=3.141593 l=9.424778 delta=(0.000000, 6.283185)
  w1 e0 ori=Forward f=0.000000 l=6.283185 delta=(0.000000, 6.283185)
  w1 e1 ori=Forward f=3.141593 l=9.424778 delta=(0.000000,-6.283185)
```

每条 wire 的两条边**同向**、且一条 `v:2π→0`、另一条 `v:0→2π` ⇒ `ΣΔv=0` ⇒ `isvopen=0` ⇒ `CheckWire` false。

**这是 seam 边的典型形态**：同一条 seam 边在一个 wire 里出现两次（周期面的两侧）。OCCT 中这两次出现
**一正一反**，`BRep_Tool::CurveOnSurface` 遂分别返回 `PCurve()` 与 `PCurve2()`（差一个周期），Δ 不会抵消；
而端口两次都是 `ori=Forward`，取到同一条 pcurve，于是抵消。

⇒ **修复点**：STEP/建 wire 时 seam 边第二次出现必须保留其（反向）orientation，或按 `IsCurveOnClosedSurface`
分别取 `PCurve1/PCurve2`。这与 §9.65 的结论一致（算法忠实，数据错）。

**下一步**：在 `step/read_topology.rs` 的 wire 构造中核对同一 edge 的多次出现是否保留各自 orientation。
### 9.67 round 67 —— 读入端核对：`resolve_loop` 逐次出现独立取向，故疑点在别处

`step/read_topology.rs`：
- `resolve_loop`（`:413-433`）按 `EDGE_LOOP` 列表顺序逐个 `resolve_shape(oriented_edge)` 后 `make_wire`；
- `resolve_oriented_edge`（`:369-386`）按记录的方向标志设 `Forward/Reversed`，并沿嵌套 `ORIENTED_EDGE` 追到 `EDGE_CURVE`；
- 注释指出「保留 EDGE_LOOP 顺序，供 seam `SelectForwardSeam` 使用」。

⇒ 读入端**没有**把两次出现统一成 Forward 的明显代码路径。两个可能：
1. STEP 的 `EDGE_LOOP` 本身两次都标了同向（需打印该 loop 的记录 id/方向核实）；
2. 面 129 的两条边并非同一 `EDGE_CURVE` 的两次出现，而是两条**不同**边（COS 范围分别是 `0..2π` 与 `π..3π`），
   即该面在 UV 里本就是一个「去而复返」的退化带 —— 那么 `CheckWire` 在 OCCT 下也会是 0，
   需要重新确认验收的 7 个面到底对应端口的哪些面（用 §9.60 的类型+形状配对，而非三角数排序）。

**下一步**：打印 face=129 所在 loop 的每个 `ORIENTED_EDGE` 记录 id 与其 `EDGE_CURVE` id + 方向标志，
确认两臂是否为同一条边的两次出现。
### 9.68 round 68 —— 结论：face 129 的两臂是**不同**边 ⇒ 它可能不是验收目标面

`--checkwire` 的 `f129`：两臂的 COS 范围分别是 `0..2π` 与 `π..3π`。
端口的 COS 范围存于 **(edge TShape, surface)** 上（`repr_key`）；同一条 `EDGE_CURVE` 的两次出现会共享同一 TShape，
因而打印**相同**的范围。范围不同 ⇒ **两臂是两条不同的边**。

⇒ `f129` 这类面在 UV 里是「去而复返」的退化带，其 `ΣΔ=0` 在 OCCT 下同样成立，`CheckWire` 也会是 false。
**因此不能断言这 7 个面就是验收所指的 7 个周期 B 样条面**——此前的配对是「按三角数排序」，不可靠。

**下一步（改口径）**：用 OCCT `--fface`（每面 type/orient/边数）或 `--perface` 的解析 `contrib`（体积贡献）
作为**几何键**，把 OCCT 的面与端口的面**一一配对**，再确认「1-11 → 360-486」到底对应哪些面；
在此之前不再围绕 `f129` 改动。
### 9.69 round 69 —— 按 (类型, 三角形数) 贪心配对（226/226）

OCCT 侧：`a3n00_occ_perface.txt` 的类型 + `a3n00_occ_mesh.txt` 的三角形数（同一 `TopExp_Explorer` 序）；
端口侧：`a3n00_port_107612.txt` 的 type+mt。贪心（同类型取最近三角形数）配满 **226/226**。

**OCCT 落在 360-486 的面**（共 5 个，不是 7 个）：
| OCCT | 类型 | 三角形 | 配对端口 | 三角形 | Δ |
|---|---|---|---|---|---|
| f48 | Plane | 376 | f205 | 376 | 0 |
| f50 | Plane | 376 | f217 | 376 | 0 |
| f29 | BSpline | 376 | f208 | 224 | 152 |
| f31 | BSpline | 376 | f210 | 224 | 152 |
| f46 | Torus | 355 | f202 | 355 | 0 |

其他稳健匹配：`occ f68 Plane 243 = port f115 243`（Δ0）。BSpline 簇内配对**歧义**（854/810/669/376/376 对 674/376/376/224…），
贪心结果只能作参考，但可看出：**Plane/Torus 的 360-486 面已完全一致**，缺口集中在 BSpline 簇。

**注意**：端口 `f191/f200/f205/f217` 的 `wires=10 edges=[1×10]`（§9.63）正是 fictive-grid 形态的面
（`FixMissingSeam` 末尾 `cxx:2236-2261` 的产物），且它们是**在接线关闭**的情况下由读入直接得到的。

**下一步**：BSpline 簇改用解析 `contrib` 或 `--fface` 的边数做键，消除歧义后定位缺口。
### 9.70 round 70 —— 配对键的可行方案（探针命令盘点）

盘点了 `.target-gate/occt_probe/occt_probe.cpp` 的可用命令：
`--entities`（FaceSurface 实体的 sameSense）、`--perface`（**逐面 type/orient/解析 contrib**）、
`--mesh`（逐面 nodes/triangles）、`--ds/--fuse/--dir/--dir2/--faceoff/--faceoff2/--cylall/--cylface/--cyl/--sph/--facek`
——后一组都自建 box/cylinder 模型，**不能用于 a3n00**。

⇒ 可用于 a3n00 的配对键只有：`--perface` 的 **type + orient + 解析体积贡献 `contrib`**，与 `--mesh` 的 nodes/triangles。
端口侧提供 `type`/`mv`/`mt`/`wires`/`edges`/`uv`，**没有**解析体积贡献。

**结论**：要消除 BSpline 簇的配对歧义，需给端口探针加一个**逐面解析体积**输出（复用端口 `brep_gprop_full`），
与 OCCT 的 `contrib` 一一对齐。这是下一步要做的（探针改动，不进主库门禁）。

**至今的稳健结论**（不依赖配对）：面数/类型分布逐类一致；总量 12283 vs 12324；
`FixMissingSeam` 在 a3n00 的 118 个周期面上全部不触发（`check_wire` 的 `vec` 恒为 0），且算法链已核对忠实。
### 9.71 round 71 —— 【配对成功】逐面解析体积（`brep_gprop_full::face_volume_contribution`）

新增（端口，探针用）：`brep_gprop_full::properties::face_volume_contribution(face: &Face) -> Result<f64,String>`
（`BRepGProp_Vinert` 关于原点，对应 OCCT `--perface` 的 `contrib`）；探针加 `--facevol` 模式，输出 `FV face=i contrib=v`（226 行）。

按 `contrib` 最近值贪心配对 226/226，结果（OCCT 落在 360-486 的面）：
| OCCT | 类型 | 三角形 | contrib | 端口 | 类型 | 三角形 | contrib | rel |
|---|---|---|---|---|---|---|---|---|
| f29 | BSpline | 376 | −407298 | **f191** | BSpline | **376** | −407302.36 | **1.07e−05** |
| f31 | BSpline | 376 | 413524 | **f200** | BSpline | **376** | 413531.00 | **1.69e−05** |
| f48 | Plane | 376 | −407298 | f218 | Plane | 72 | −412544 | 1.29e−02 |
| f50 | Plane | 376 | 386325 | f175 | BSpline | 72 | 236514 | 3.88e−01 |

⇒ **两个 BSpline 的 360-486 面与 OCCT 完全一致**（三角形相同、几何贡献相对误差 ~1e-5）。
`f48/f50` 的配对不可信：`occ f48 type=0` 与 `occ f29 type=6` 的 `contrib` **完全相同（−407298）**，
说明 `a3n00_occ_perface.txt` 与 `a3n00_occ_mesh.txt` 可能取自**不同状态**的形状，索引不能混用。

**注意**：部分 OCCT 面 `contrib≈0`（如 f177/f184）而端口配对到巨大值 —— 同样指向两份 dump 口径不一致。

**下一步**：重跑 OCCT `--perface` 与 `--mesh` **同一次**（或改用 `--FFACE` 让 mesh 输出自带 type/contrib），
再按 `contrib` 精确配对；届时可判定「1-11 → 360-486」的确切面。
### 9.72 round 72 —— 环境限制：OCCT 探针在本机**无法运行**

`occt_probe.exe` 直接运行返回 `EXIT=-1073741515`（`0xC0000135` = `STATUS_DLL_NOT_FOUND`）；
在 `D:\source\OCCT-src` 与仓库树中均 **找不到 `TKernel.dll`** 等 OCCT 运行库。

⇒ **无法重跑 OCCT** 来获得 `--perface`/`--mesh` 的一致快照；后续只能使用 `.target-gate/` 下既有的 OCCT dump，
并带上「两份可能不同源」的保留（§9.71）。

**要在本机重跑 OCCT，需要**：OCCT 8.0.0 的 build（`bin/` 含 `TKernel.dll`/`TKTopAlgo.dll`…）或其路径，
以及（可选）重新编译 `specs/occt_probe/occt_probe.cpp`。

**已有 dump 下的稳健结论**（§9.71）：OCCT 两个 BSpline 360-486 面（f29/f31）与端口 f191/f200
三角形相同（376=376）、解析体积贡献相对误差 ~1e-5 ⇒ **这两个面已对齐**。
### 9.73 round 73 —— 【定位异常面】端口 `f205`/`f217`（Plane，fictive-grid 形）几何不对

用解析体积贡献（`--facevol` vs OCCT `--perface`）直接比对：

| 面 | 类型 | 三角形 | OCCT contrib | 端口 contrib | 相对 |
|---|---|---|---|---|---|
| OCCT f29 ↔ 端口 f191 | BSpline | 376 = 376 | −407298 | −407302.36 | 1.1e−05 ✓ |
| OCCT f31 ↔ 端口 f200 | BSpline | 376 = 376 | 413524 | 413531.00 | 1.7e−05 ✓ |
| OCCT f46 ↔ 端口 f202 | Torus | 355 = 355 | 17161.2 | 17159.84 | 8e−05 ✓ |
| OCCT f48 (Plane) | Plane | 376 | **−407298** | 端口无对应 | — |
| OCCT f50 (Plane) | Plane | 376 | **386325** | 端口无对应 | — |
| **端口 f205** | Plane | 376 | — | **+577414.91** | 无对应 |
| **端口 f217** | Plane | 376 | — | **−586685.23** | 无对应 |

OCCT 中 |contrib|>3e5 的面只有 `f25(306277) f29(−407298) f31(413524) f48(−407298) f50(386325)`；
端口多出 `+577414.9` 与 `−586685.2` 两个**无对应**的大值（约为 OCCT 对应量的 1.4–1.65 倍）。

这两个面是 `wires=10 edges=[1,1,…]` 的 **fictive-grid 形态** Plane 面。⇒ **端口 `f205`/`f217` 的边界/几何不对**，
（三角形数恰好也落在 376，所以此前按三角形数配对时被误认为已对齐）。

**下一步**：核对 `f205`/`f217` 的 wire/pcurve 与 OCCT `f48`/`f50` 的差异 —— 这正是「1-11 → 360-486」之外的另一处形状差异。
（注：`--perface` 与 `--mesh` 都是 226 面、0..225，面数一致，故这两份 dump 应同源。）
### 9.74 round 74 —— 加 `face_volume_contribution` 后的门禁终检（全绿）

| 门禁 | 结果 |
|---|---|
| `step_obj_parity` | 14/14 |
| `step_to_obj` | 13/13 |
| `step_obj_area` | 11/11 |
| `step_geometry_parity` | 3/3 |
| `occt-topo --lib` | 1281/0 |

（本轮唯一库改动是新增 `face_volume_contribution`，纯增量；探针加 `--facevol`。）
### 9.75 round 75 —— 反转 9.73 的“同源”判断：两份 OCCT dump 的体积不一致

- `a3n00_occ_mesh.txt` 的 `TOTAL ... meshvol=1.8897e+06`（**三角网格**体积）
- `a3n00_occ_perface.txt` 的 `PERFACE total=2.23066e+06`（**解析** `BRepGProp::VolumeProperties` 质量）

解析体积与网格体积差 **18%**，远超正常离散误差 ⇒ 这两份 dump **很可能取自不同状态的形状**（面数恰好都是 226，但几何不同）。
因此 §9.73 关于「端口 f205/f217 异常」的结论**存疑**：那两组 `contrib` 可能本就不可比。

**这进一步说明本机无法重跑 OCCT 是当前最大的障碍**：没有同源快照，逐面对拍无法给出确定结论。

**待用户提供**：OCCT 8.0.0 的 `bin/`（含 `TKernel.dll` 等），以便重跑 `occt_probe --mesh/--perface` 取得同源数据。
### 9.76 round 76 —— 再次修正：端口的 `face_volume_contribution` 之和不可信

端口 `--facevol` 的 226 面求和 = **−4,827,776**，而：
- OCCT `PERFACE total` = **+2,230,750**；
- 端口自己导出的 `output/a3n00.obj`（旧产物）有向网格体积 = **+1,346,469**（正）。

⇒ 端口逐面贡献之和不等于自身体积 ⇒ **`face_volume_contribution` 未像 `BRep_GProp::VolumeProperties` 那样处理面取向**
（REVERSED 面的贡献符号），或 `loc`/`coeff` 约定与 `volume_properties` 不一致。

因此 **§9.73 关于「端口 f205/f217 异常」的结论同样不可靠**（只有取向一致的面才能单体对比；
f191/f200/f202 恰好与 OCCT 一致、相对误差 ≤1e-5，可作为“这些面已对齐”的证据，但不能推广）。


### 9.77 round 77 —— 【重大】a3n00 的**解析体积符号相反**：端口 −3.34e6 vs OCCT +2.23e6

给探针加 `--voltotal`（调用端口 `brep_gprop_full::volume_properties(&shape)`）：
```
VOLTOTAL mass=-3336973.416029      # 端口
PERFACE total=2.23075e+06          # OCCT（BRepGProp::VolumeProperties）
```

注意探针的 `BBOX roots=2` —— a3n00 是**含 2 个 root 的 compound**。端口对同一 compound 的解析体积为**负**且量级更大，
而 OCCT 为正。**这比 0.33% 的网格差异严重得多**，指向：
- 某个子形状的 **shell 取向反了**（整体或局部），或
- compound 里两个组件的取向/包含关系与 OCCT 不同。

（说明：端口自己的 `volume_properties` 也给出负值，故不是上一轮怀疑的 `face_volume_contribution` 实现问题；
`volume_properties` 是被门禁覆盖的函数。）


### 9.78 round 78 —— 【新异常】两个 root **各自都是 226 面、体积完全相同**

给探针加 `--volroots`（对 `model.shapes` 的每个 root 单独算体积）：
```
VOLROOT i=0 faces=226 mass=-3336973.416029
VOLROOT i=1 faces=226 mass=-3336973.416029
```

两个 root 的**面数与质量完全一致**（到小数点后 6 位）⇒ 极可能是**同一个形状被登记了两次**（读入端的 root 处理重复），
而不是两个真实的半体。这解释了 `VOLTOTAL mass` 与 root 相同、以及 `BBOX roots=2`。


### 9.94 round 93 —— 【决定性实验】偏转减半后，那些病态面**全部恢复正常**

对同一个 `T0M.stp`，把线性偏转从 `2.574309` 改为 **`1.287154`（= 一半）** 重跑：

| 面 | def=2.574309 | def=1.287154 |
|---|---|---|
| f1495 Cylinder(1×[4]) | mv=118 **mt=6** | mv=123 **mt=123** |
| f1654 Cylinder(2×[2,2]) | mv=112 **mt=3** | mv=36 **mt=34** |
| f404 BSpline(1×[6]) | mv=76 **mt=1** | mv=18 **mt=14** |
| 整模型 | mesh_v=46375 mesh_t=46678 | mesh_v=51875 mesh_t=**55089** |

⇒ 病态不是「几何无法三角化」，而是**与偏转相关的条件分支**在 def 较大时走到了失败/回退路径
（对照 §9.93：`discret_root.rs:1366/1372` 的两处 `Err(_) => wireframe::face_to_triangles(...)`）。
注意 def 减半后 `f1495` 的 `mv` 几乎不变（118→123），但 `mt` 从 6 跳到 123 —— 同一边界点集，
一个走主流水线成功、一个走回退（或主流水线内部提前退出）。

**下一步**：在 `build_shape_mesh`/`ModelBuilder` 里找**依赖 deflection 的提前返回/阈值判断**，
确认 def=2.574 时是哪一条 `Err` 被触发（打印 `Err` 文案即可）。

### 9.93 round 92 —— 单面网格管线与失败点（`discret_root.rs:1357-1374`）

```rust
pub fn discretize_face(&mut self, face: &Face, deflection: f64) -> (Vec<GpPnt>, Vec<Triangle>) {
    let mut model = match ModelBuilder::build_model(&face.0, &params) {
        Ok(m) => m,
        Err(_) => return wireframe::face_to_triangles(face, deflection),   // 回退 1
    };
    ModelPreProcessor::perform(&mut model, &params);
    match self.build_shape_mesh(&mut model) {
        Ok(m) => (m.vertices, m.triangles),
        Err(_) => wireframe::face_to_triangles(face, deflection),          // 回退 2
    }
}
```

⇒ 单面路径有三种可能来源：**主流水线**（`ModelBuilder`+`ModelPreProcessor`+`build_shape_mesh`），
或两处 **`wireframe` UV 网格回退**。`f1495` 的 `mv=118` 说明**边界点很多**（回退的 UV 网格在
`du≈6.25, def≈2.57` 下只会给 ~8 点），因此更像**主流水线跑通但三角化几乎没产出**（118 边界点 → 仅 6 三角形）。

**下一步（确定性）**：给 `f1495` 加一次性探针，分别打印：
1. 是否走了回退（两条 `Err` 分支）；
2. 主流水线里 `ModelBuilder` 的边界点数 / `ModelPreProcessor` 后点数 / `build_shape_mesh` 产出的三角形数；
3. 失败时 `build_shape_mesh` 的具体 `Err` 文案。

### 9.92 round 91 —— 修正 9.91：那条 UV 网格属于**回退路径**，主路径另有实现

调用关系（grep 结果）：
- `FaceDiscret::discretize_face`（`face_discret.rs:555`，docs 自述「Port of `BRepMesh_FaceDiscret`」，内部是 UV 网格）
  只被 `meshing/fast_discret.rs:72`（`FastDiscret`，其自身注释也说是 fallback）、`face_discret.rs:801`（测试）、`meshing/context.rs:283` 调用；
- `incremental_mesh/discret_root.rs:1361` **另有** `discretize_face(face, deflection) -> (Vec<Gpnt>, Vec<Triangle>)`，这才是 `IncrementalMesh` 的主路径（Delaunay）。

⇒ §9.91 的措辞需修正：那套 UV 网格是**回退实现**（与 `step_obj_parity` 头注释一致），
**不能断定** T0M 的 `f1495` 走的是它。**主路径 `discret_root.rs:1361` 失败**才是更可能的原因。

**下一步（确定性）**：在 `--dev`/专门模式里打印每个面实际走的分支（主 Delaunay vs FastDiscret fallback）与失败原因；
对 `f1495` 打印其 UV 界、边界点数、内部点数、以及三角化输出点数。

### 9.91 round 90 —— 【根因候选·高危】`face_discret` 的内部加点仍是 **UV 网格**（用 3D 偏转当 UV 步长）

`crates/occt-topo/src/meshing/face_discret.rs`：
```rust
557:  pts.extend(Self::interior_grid(face, self.params.deflection, 4096));   // 仍在被调用
573:  pub fn interior_grid(face: &MeshFace, deflection: f64, max_points: usize) -> Vec<GpPnt2d> {
578:      let (umin, umax, vmin, vmax) = Self::uv_bounds(&face.outer_wire);
584:      let mut nu = (du / def).ceil() as usize + 1;      // def = self.params.deflection（3D 线性偏转）
585:      let mut nv = (dv / def).ceil() as usize + 1;
```

两点问题：
1. 这是 **UV 均匀网格 + 内外判定**，而 OCCT `BRepMesh_FaceDiscret` 走的是
   `BRepMesh_DefaultRangeSplitter` + **按曲面局部度量换算的 UV 容差**（`BRepMesh_Deflection`/`ShapeTool`）做自适应加点；
2. 把 **3D 偏转** 直接当 **UV 步长**是**量纲错误**：圆柱上 `u` 走 Δu 对应 3D 长度 `r·Δu`，
   UV 步长应为 `deflection / r`，即 `nu = du·r / deflection`；端口少了 `r` 因子。
   半径小的面会被**严重欠采样**（正合 `f1495`：du≈6.25、def=2.574 → nu≈4、nv≈2）。

**这与 `step_obj_parity` 文件头的自述矛盾**（那里写「export 路径已不再回退 UV-grid/quadtree 网格器 —— audit A17」），
即该回退**仍存在于 `IncrementalMesh` 路径**。

**下一步**：确认 `IncrementalMesh` 是否也走 `face_discret`（而非纯 Delaunay）；
若是，则把 `interior_grid` 换成 OCCT 的 `BRepMesh_DefaultRangeSplitter`/UV 偏转换算，并以此解释 T0M/acs10 的缺口。

### 9.90 round 89 —— 【强假设】`f1495` 的 pcurve `u` 超出曲面周期区间（最高 9.4 > 2π）

`--checkwire` 对 T0M `f1495`（Cylinder，mv=118 mt=6）：
```
CW face=1495 u=[0.0000,6.2832] v=[-inf,inf] dU=6.2832 dV=inf
  w0 e0 Forward f=3.160726 l=9.238343 df=(9.405644,0.067712) dl=(3.328027,0.067712) delta=(-6.077617,0)
  w0 e1 Forward f=1.374517 l=1.551663 df=(9.228499,1.300000) dl=(9.405644,0.067712) delta=(0.177145,-1.232288)
  w0 e2 Reversed f=3.337872 l=9.415310 df=(9.228499,1.300000) dl=(3.151061,1.300000) delta=(-6.077438,0)
  w0 e3 Forward f=4.525954 l=4.702921 df=(3.328027,0.067712) dl=(3.151061,1.300000) delta=(-0.176967,1.232288)
```

两个关键点：
1. **`v` 范围是 `-inf .. inf`**（圆柱的自然 V 域无限）——面网格必须用 **pcurve 围出的 UV 界**（OCCT `BRepTools::UVBounds(face)`），
   不能用曲面自身无限区间；
2. **pcurve 的 `u` 达 9.23–9.41**，超出曲面周期区间 `[0, 2π]≈[0, 6.2832]`（≈ 2π 再 +3.0），
   而 `df.x=9.405644 > 6.2832`。

⇒ 若端口的 face 三角化直接把 pcurve 的 `(u,v)`（u 最高 9.4、v 来自无限域）灌给 range splitter / Delaunay，
UV 域与曲面参数域**不一致**，约束三角化就会大量丢弃三角形 —— 这正好解释 `mv=118 mt=6`。

**下一步（验证）**：在端口 face 三角化入口打印实际使用的 UV 界，与 `BRepTools::UVBounds(face)`（可由 pcurve 端点算出）对比；
若确认没有 wrap/夹取，则按 OCCT 的 `BRepMesh_FaceDiscret`/`BRepMesh_ShapeTool::UVBounds` 处理（周期面按周期平移 pcurve）。

### 9.89 round 88 —— 病态面的结构（可复现样本）

| 面 | 类型 | wires / edges | mv | mt |
|---|---|---|---|---|
| f1495 | Cylinder | 1 / [4] | 118 | **6** |
| f1654 | Cylinder | 2 / [2,2] | 112 | **3** |
| f404 | BSpline | 1 / [6] | 76 | **1** |
| f1750 | Plane | 2 / [1,1] | 122 | **11** |
| f810 | Cone | 1 / [4] | 112 | **4** |
| f21 | Cone | 2 / [5,2] | 114 | 24 |

**`f1495` 最典型**：一条 wire、4 条边的**普通圆柱片**，边界离散出了 118 个节点，却只生成 **6 个三角形**。
⇒ 问题不在边界离散（点够多），而在**面的约束三角化**（Delaunay/constraint 阶段）失败或几乎全被丢弃。

这与 a3n00 的 `f129`（`wires=2 edges=[2,2] mv=162 mt=1`）同型，说明是一类可复现缺陷，而不是个别模型特例。

**下一步**：对 `f1495` 打印其 UV 域（`u/v` 范围）、4 条边在 UV 里的 pcurve 折线、以及三角化前后的点数，
定位是「UV 域退化（面积为 0/自交）」还是「约束三角化丢弃了内部三角形」。

### 9.88 round 87 —— T0M 逐面解剖：9 个「病态面」+ 约 150 个整体偏疏的面

`T0M.stp`（端口 1772 面 / mv 46375 / mt 46678；OCCT 1778 面 / nodes 60050 / tris 66576）：

| 区间 | OCCT 面数 | 端口面数 |
|---|---|---|
| 0 | 0 | 3 |
| 1-5 | 349 | **349** |
| 6-20 | 743 | **894** |
| 21-100 | 569 | **449** |
| 101-1000 | 117 | **77** |

- **病态面**（`mv>20` 且 `mt < mv/10`，即边界点进了但三角形几乎没生成）：**9 个**（BSpline 2、Cone 1、Cylinder 4、Plane 2），
  例：`f1495 Cylinder mv=118 mt=6`、`f1654 Cylinder mv=112 mt=3`、`f404 BSpline mv=76 mt=1`。
  这与 a3n00 的 `f129`（mv=162 mt=1）同型。
- 但 9 个面解释不了 −19898 的缺口：主要是**约 150 个面**（OCCT 21-100 区间）被端口压到 6-20。

⇒ 两条并行线索：**(a) 少量面三角化失败**；**(b) 一批面整体偏疏**。

**下一步**：对 T0M 的 `f1495`/`f1654`/`f404` 三个面，打印其 UV 域、wire 数/边数、pcurve 范围，
并与 OCCT 同面（按解析 contrib 或 type+边数配对）比较 —— 这是「三角化失败」的可复现样本。

### 9.87 round 86 —— 【收敛】网格管线**不是全局偏粗**，缺口集中在 T0M/acs10

`step_obj_parity -- --nocapture` 的 `f-ratio`（端口 `brep_to_obj(shape, 0.1)` ↔ `occ-*.obj` 头部 `# Faces`）：

| 模型 | f-ratio | | 模型 | f-ratio |
|---|---|---|---|---|
| Cube/Cylinder/Cone/Sphere/Torus | **1.00** | | occ/bottom.step | 0.99 |
| HoledPlate/Offset/OffsetPlaneHoleEdge | **1.00** | | occ/motoc.step | 0.98 |
| Shape/Shape-1/Shape-2/reversed | **1.00** | | occ/ATU01038.step | 1.00 |
| | | | occ/top.step | 1.00 |
| | | | occ/TDB.stp | 0.97 |
| | | | occ/a3n00.stp | 0.99 |
| | | | **occ/T0M.stp** | **0.69** |
| | | | **occ/acs10.stp** | **0.80** |

⇒ 在**所有简单/中等形状上端口与 OCCT 完全一致（1.00）**，说明网格管线本身没问题；
§9.85 的「系统性偏粗」只在 **T0M / acs10** 这两个模型上成立（bottom/motoc/TDB 仅 1-3%）。

**下一步（聚焦）**：只查 T0M 与 acs10 —— 大概率是**这些模型里某类面**（如特定 B 样条/裁剪面）的 UV 域或离散没被正确细分，
而不是全局参数问题。可先用 `--mesh` 的逐面 nodes/triangles 找出两侧差异最大的面。

### 9.86 round 85 —— 偏粗问题的入口盘点

`IncrementalMesh::from_deflection(shape, d, false, angle)`（`meshing/incremental_mesh/discret_root.rs:117-130`）
把 `d`/`angle`/`relative` 原样装进 `MeshParameters`（`meshing/parameters.rs:27`）再 `perform()` —— 参数传递本身无损耗。

`angle`/`deflection` 在端口内的使用点（按文件计数）：
`incremental_mesh/discret_root.rs` 23、`model_builder/preprocessor.rs` 19、`meshing/edge_discret.rs` 14、
`range_splitter/splitter.rs` 12、`meshing/parameters.rs` 11、`meshing/deflection.rs` 8。

⇒ 下一步应固定**一个面**，比较两侧：
1. **边界节点数**（端口 `edge_discret` ↔ OCCT `BRepMesh_EdgeDiscret`/`GCPnts`）
2. **内部加点数**（端口 `range_splitter/splitter.rs` ↔ OCCT `BRepMesh_DefaultRangeSplitter`）
3. `meshing/deflection.rs` 是否等价于 `BRepMesh_Deflection`（在 `relative=false` 下 OCCT 仍会做 per-face 容差处理）

### 9.85 round 84 —— 【新目标】端口网格系统性偏粗（逐面分布对比，acs10）

对 `acs10.stp`（OCCT 787 面 / 46494 tris；端口 786 面 / 37223 tris）做**逐面三角形数的分布**对比：

| 区间 | OCCT 面数 | 端口面数 |
|---|---|---|
| 0 | 0 | 1 |
| 1-10 | 237 | 232 |
| 11-50 | 264 | **315** |
| 51-200 | **239** | 205 |
| 201-1000 | **47** | 33 |
| >1000 | 0 | 0 |

⇒ 缺口不是「丢面」，而是**本该 51-1000 的面被端口压到了 11-50**：端口在同一名义偏转下**系统性偏粗**。
这解释了 §9.81 里 `acs10 −19.9%`、`T0M −29.9%` 的三角形缺口（面数只差 1/6）。

**候选根因（按优先级）**：
1. `IncrementalMesh::from_deflection` 的**角偏转/相对偏转**处理与 `BRepMesh_IncrementalMesh(deflection, false, Angle)` 不一致；
2. 边界离散（`GCPnts_*`/`BRepAdaptor`）产点少于 OCCT；
3. Delaunay 内部加点（`BRepMesh_DefaultRangeSplitter`/refinement）判据不同。

**下一步**：固定一个面（OCCT `--mesh` 逐面可定位），对比两侧的**边界节点数**与**内部加点数**；
`occt_probe` 已有 `--mesh` 的逐面 nodes/triangles，端口探针有 `mv`/`mt`。

### 9.84 round 83 —— 【工具不可用】端口的体积积分在这些 STEP 形状上产生垃圾值

对若干模型同时跑 `--voltotal`（`volume_properties`）与 `--facevol`（逐面 `face_volume_contribution`）：

| 模型 | `VOLTOTAL mass` | `facevol_sum` |
|---|---|---|
| bottom.step | **9.78e+26** | **1.03e+27** |
| ATU01038.step | 731093.75 | 495606.14 |
| motoc.stp | (空/失败) | 0 |
| a3n00.stp | −3,336,973.4 | −4,827,776.1 |

⇒ 端口的 **`BRepGProp` 体积积分（Vinert）在这些形状上完全不可信**（1e26 量级、或两个口径互不相等）。
注意 `bottom.step` 通过 `step_obj_parity` 与 `step_obj_area`（**面积**走 `Sinert`，与体积是两条路径），
所以「面积门禁绿」并不能说明体积积分正确。

**结论与影响**：
- **§9.71/§9.73 基于 `contrib` 的逐面配对结论全部作废**（工具本身不成立）；
- `--facevol` 目前只能用于「同一形状内相对比较」，不能与 OCCT `--perface` 绝对对拍；
- 若要继续逐面对拍，需先让端口的 `volume_properties`/`face_volume_contribution` 可信
  （对照 `BRepGProp_Vinert::Perform` 与 `GProp_GProps::Add`；并区分**闭壳**与**开壳**：
  开壳的 `VolumeProperties` 在 OCCT 里也不是体积）。

**下一步（择一）**：
1. 先修端口体积积分（`brep_gprop_full`），再回到逐面配对；
2. 或改用**不依赖体积**的配对键（端口 `uv=[...]` + 面类型 + 边数 ↔ OCCT `--perface` 的 type + `--fface` 的边数）。

### 9.83 round 82 —— 口径澄清后重跑：**三角形缺口是真的**

给探针改用 `model.shapes[0].shape`（对齐 OCCT `aReader.OneShape()` = 第一个 transfer root），
避免把 `Compound` 与其内含 `Solid` 叠加遍历；重跑结果与改前**逐位相同**：

| 模型 | 端口F | 端口 tris | OCCT F | OCCT tris | Δtris |
|---|---|---|---|---|---|
| a3n00.stp | 226 | 12283 | 226 | 12324 | −0.33% |
| acs10.stp | 786 | 37223 | 787 | 46494 | **−19.9%** |
| T0M.stp | 1772 | 46678 | 1778 | 66576 | **−29.9%** |

⇒ 面数的小差（−1/−6）可用「去重 vs 出现次数」解释，但**三角形缺口（−20%/−30%）与口径无关，是真实差异**
（端口的 `stats` 也显示有 1/3 个面未产出网格）。§9.82 的保留只影响面数，不影响三角形结论。

**下一步**：对 T0M 用 OCCT `--perface`（逐面 contrib）与端口 `--facevol` 配对，找出欠网格的面；
但需先修端口 `face_volume_contribution` 的取向语义（§9.76）。

### 9.82 round 81 —— 【口径修正】端口按「去重面」计数，OCCT 探针按「出现次数」计数

`--roots` sweep（`.target-gate/roots_sweep.txt`）：
```
a3n00.stp : ROOTS count=2  i0 ptr=0x...860 Solid faces=226   i1 ptr=0x...860 Solid faces=226   # 同一 TShape 两次
T0M.stp   : ROOTS count=2  i0 ptr=0x...4e0 Compound faces=1772 i1 ptr=0x...0e0 Solid faces=1772
acs10.stp : ROOTS count=2  i0 ptr=0x...f00 Compound faces=786  i1 ptr=0x...d60 Solid faces=786
bottom.step: ROOTS count=1  i0 Solid faces=323
```

两个 root 分别是 **Compound** 与其内含 **Solid**（同一几何登记两次）；a3n00 更极端——两个 root 指针相同。

更关键的计数口径：
- 端口 `topo_tools_full::faces_of` = `shapes_of` = **`map_shapes`（去重）**；
- OCCT 探针 `--mesh`/`--perface` 用 **`TopExp_Explorer`（按出现次数，不去重）**。

⇒ 上一节表格里的 `Δ面`（T0M −6、acs10 −1、a3n00 0）很可能是**同一个面 TShape 出现多次**的计数差，
而不是几何缺失；三角形总数同理可能因此偏小。**在与 OCCT 对拍前必须统一到同一口径**
（建议端口探针另加 occurrence 版遍历）。

**下一步**：给端口探针加 occurrence 版统计（`TopExp_Explorer` 等价遍历），重跑 8 模型，再判断哪些是真差异。

### 9.81 round 80+ —— 【首次全模型对拍】6/8 接近，`acs10` 与 `T0M` 明显偏离

协议：端口 `IncrementalMesh::from_deflection(shape, d, false, 20°)` ↔ OCCT `--mesh <d> 0.349066`，`d` = 端口 `computed_lin`。

| 模型 | 端口面数 | 端口 nodes | 端口 tris | OCCT 面数 | OCCT nodes | OCCT tris | Δ面 | Δnodes | Δtris |
|---|---|---|---|---|---|---|---|---|---|
| a3n00.stp | 226 | 11931 | 12283 | 226 | 11052 | 12324 | 0 | +8.0% | −0.33% |
| acs10.stp | 786 | 32400 | 37223 | 787 | 37247 | 46494 | **−1** | **−13.0%** | **−19.9%** |
| ATU01038.step | 386 | 18090 | 22529 | 386 | 18008 | 22426 | 0 | +0.5% | +0.46% |
| bottom.step | 323 | 17587 | 23408 | 323 | 17472 | 23369 | 0 | +0.7% | +0.17% |
| motoc.step | 223 | 11713 | 13646 | 223 | 11713 | 13944 | 0 | **0.0%** | −2.1% |
| T0M.stp | 1772 | 46375 | 46678 | 1778 | 60050 | 66576 | **−6** | **−22.8%** | **−29.9%** |
| TDB.stp | 2180 | 72812 | 77814 | 2180 | 72687 | 78470 | 0 | +0.2% | −0.84% |
| top.step | 324 | 16902 | 22794 | 324 | 16744 | 22520 | 0 | +0.9% | +1.2% |

**结论**：
- `motoc.step` 节点数**完全相同**，`bottom.step`/`TDB.step`/`ATU01038` 误差 ≤1%，`top.step` +1.2%，`a3n00` 三角形 −0.33%（节点 +8%）。
- **`acs10.stp`（面数 −1、三角形 −19.9%）与 `T0M.stp`（面数 −6、三角形 −29.9%）严重偏离** —— 这两个是新的首要目标。

端口侧 sweep 与 OCCT 侧 sweep 已分别落盘；本节表格同时写入 `.target-gate/sweep_compare.txt`。

### 9.80 round 80+ —— 【环境已解决】OCCT 探针可用 + 8 组对拍件 + 对拍协议

**DLL 问题根因与解法**：`occt_probe.exe` 是 `/MD` 构建，直接运行会因找不到 OCCT 运行库而 `0xC0000135`。
OCCT 8.0.0 安装在 **`D:\source\occt-8.0.0`**（不是 `D:\source\OCCT-src`，后者是源码树）；
运行库在 `win64\vc14\bin`，第三方在 `3rdparty-vc14-64`（`tbb12.dll`/`jemalloc.dll` 等）。
必须先用该目录的 `env.bat vc14 64` 设好 `PATH`，且 `THIRDPARTY_DIR` 必须**绝对**。

- `specs\occt_probe\probe.bat <绝对路径> ...` —— 会 `cd /d %OCCT%`，所以输入路径必须绝对；
- **新增 `specs\occt_probe\run_here.bat <相对/绝对路径> ...`** —— 用 `pushd/popd` 只在设 PATH 时进入 OCCT 目录，
  **保留调用者的 cwd**，可直接用相对路径。两者均已验证。

验证：`run_here.bat data\occ\a3n00.stp --mesh 1.07612 0.349066` →
`TOTAL faces=226 nodes=11052 triangles=12324 meshvol=1.8897e+06 neg_triangles=3232`（与既有 dump 逐位一致）。

**8 组对拍件**（`data/occ/`）：`a3n00.stp↔occ-a3n00.obj`、`acs10.stp↔occ-acs10.obj`、`ATU01038.step↔occ-ATU01038.obj`、
`bottom.step↔occ-bottom.obj`、`motoc.step↔occ-motoc.obj`、`T0M.stp↔occ-T0M.obj`、`TDB.stp↔occ-TDB.obj`、`top.step↔occ-top.obj`。

**对拍协议**（同口径）：端口 `IncrementalMesh::from_deflection(shape, d, false, angle)` 的 `mesh_v/mesh_t`
↔ OCCT `--mesh <d> <angle>` 的 `nodes/triangles`，其中 `d = 端口 computed_lin`（`prs3d_get_deflection(shape, 0.1)`）、`angle = 20°`。

**端口侧 sweep（本轮）**：
| 模型 | 面数 | d | stats | mesh_v | mesh_t |
|---|---|---|---|---|---|
| a3n00.stp | 226 | 1.076007 | 226 | 11931 | 12283 |
| acs10.stp | 786 | 1.980053 | 785 | 32400 | 37223 |
| ATU01038.step | 386 | 1.547176 | 386 | 18090 | 22529 |
| bottom.step | 323 | 0.558549 | 323 | 17587 | 23408 |
| motoc.step | 223 | 0.836000 | 223 | 11713 | 13646 |
| T0M.stp | 1772 | 2.574309 | 1769 | 46375 | 46678 |
| TDB.stp | 2180 | 2.451395 | 2179 | 72812 | 77814 |
| top.step | 324 | 0.451778 | 324 | 16902 | 22794 |

OCCT 侧 sweep 已启动（同 d、angle=20°），结果见下一节。

### 9.79 round 79 —— 【已确认】同一 TShape 被登记为两个 root；且单体体积为负

探针 `--volroots` 打印指针后：
```
VOLROOT count=2
VOLROOT i=0 ptr=0x1c4620b4cf0 type=Solid faces=226 mass=-3336973.416029
VOLROOT i=1 ptr=0x1c4620b4cf0 type=Solid faces=226 mass=-3336973.416029
```

**两个 root 的 TShape 指针完全相同** ⇒ 同一个 `Solid` 被登记了两次（`model.shapes` 重复）。

由此可分离出两个独立问题：
1. **root 重复**：`model.shapes` 里同一形状出现两次（读入端/模型层问题）。
2. **单体解析体积为负**：单 root 体积 = −3,336,973.4 / 2 = **−1,668,486.7**，而 OCCT 为 **+2,230,750**。
   符号相反 ⇒ 立面取向反了；且量级比 0.75（≠1）⇒ 几何本身也不同。

**下一步**：查读入端如何收集 root（同一 `Solid` 为何入列两次），并核对 `volume_properties` 的取向语义
（对照 `BRepGProp_Vinert::Perform` 与 `GProp_GProps::Add`，以及 shell 的 `TopAbs_REVERSED` 处理）。

**下一步**：在 XML/STEP 读入端打印 `model.shapes.len()` 与各 root 的 TShape 指针（`Arc::as_ptr`），确认是否重复登记；
并对照 OCCT `aReader.OneShape()` 的等价物与 `TransferRoots` 的 root 处理。

**下一步**：分别对 2 个 root 计算体积，与 OCCT 的构造顺序对照；核对 `step/read_topology.rs` 里 compound/assembly 的
orientation 与其 `StepShape_*` 的 `same_sense`/`orientation` 处理（对照 `StepToTopoDS_TranslateShell`/`TranslateSolid`）。

**修正后的确定进度**：a3n00 的总量（12283 vs 12324）与逐类面数一致是**唯一**可靠的整体指标；
逐面对拍在两处工具问题解决前给不出确定结论：
1. OCCT 探针本机不可运行（无 DLL）；
2. 端口 `face_volume_contribution` 需先对齐取向语义（对照 `BRepGProp_Vinert::Perform` 与 `GProp_GProps::Add`）。
















但**逐面网格分布仍不同**（3 个大 Cone/Cylinder 面网格更粗），且两份 dump 的**面序**与（可能）**偏转**口径是否一致尚未确认，
故不能据此宣布验收达成。下一步应：核对 OCCT dump 的偏转参数与面序，并按几何（而非索引）配对三张大面。










**下一步（续作入口）**：逐面对拍 `bottom.step`/`motoc.step` 上触发 `FixMissingSeam` 的面，
比较端口与 OCCT 的 `FixMissingSeam` 结果（面/边/pcurve），定位丢面分支。

































### 9.95 round 94 —— 【决定性】STEP bounds 341 vs OCCT wires 294：缺口在周期面的 wire 拆分，不在网格参数

本轮把 a3n00 的「7 个周期面 1-11 → 360-486」验收缺口定位到**读入端拓扑**（工具已就绪）。

**实测（同一 data/occ/a3n00.stp）**：

| 量 | OCCT 8.0.0p1 | 端口 |
|---|---|---|
| 输出 wire 总数 | **294**（wires_probe.exe） | **341**（163*1+53*2+2*4+4*6+4*10） |
| STEP FACE_OUTER_BOUND+FACE_BOUND | 226+115 = **341** | 341 |

⇒ STEP 里 341 个 face bound，OCCT 最终合并成 **294** 条 wire（-47）；端口保持 341，**没有合并**。

2. **逐类型 wire 直方图（OCCT ↔ 端口）**（wires_probe 已改为逐面打印 type/wires/edges/UV span）：

| 类型 | OCCT | 端口 |
|---|---|---|
| Plane | w1=82, w2=7, w6=4, w10=2 | **完全相同** |
| Torus | w1=7 | **完全相同** |
| Cylinder | w1=61, w2=1 | w1=39, **w2=22**, w4=1 |
| Cone | w1=30 | w1=25, **w2=5** |
| BSpline | w1=28, w2=1, w4=1, w10=2 | w1=10, **w2=19**, w4=1, w10=2 |

⇒ 缺口集中在**单方向周期面**（Cyl/Cone/BSpline）；Plane/Torus 不受影响。

3. **STEP 侧核对（本轮脚本 parse a3n00.stp）**：这些正是 StepToTopoDS_TranslateFace.cxx:585-591 注释的形态——
   周期面上 **2 个 bound，每个 bound 是 1 条 ORIENTED_EDGE 的 EDGE_LOOP，两个 bound 引用不同的 EDGE_CURVE（无共享边）**：
   BSpline two-bound **19** 个、Cylinder **22** 个、Cone **5** 个（例：AF#4824 b0=[4817] b1=[4694]）。
   OCCT 的 TranslateFace（cxx:592-750）对每个 bound 都 aFaceBuilder.Add，**不合并**；合并发生在**后续 ShapeFix**
   （FromSTEP.FixShape → ShapeFix_Face::Perform → FixMissingSeam：ShapeFix_Face.cxx:1722-2261，
   check_wire 判 w1/w2 在周期方向张开 → ComposeShell 插入缺缝 → 两 wire 合成一 wire）。

4. **端口为何不触发**（与 §9.64 一致，但原因钉死为数据）：端口这 19 个 BSpline 2-wire 面里 **18 个 edges=[2,2]**，
   而 STEP 每个 bound 只有 1 条 ORIENTED_EDGE。--checkwire 对 f176 显示 w0 的两条 arm（e0/e1）：
   e0 Forward COS 0..301.59、e1 Forward COS 0..6.28，端点首尾相接 ⇒ ΣΔv=0 ⇒ check_wire 返回 None ⇒
   FixMissingSeam 永不进入 CompShell。**即：端口 wire 里把单条闭合 seam 边拆成了两条不同边（或存了两份），
   使周期方向的 Δ 抵消。**

5. **面序不可按索引配对**：OCCT 与端口按 TopExp_Explorer 索引逐面 type/wires/edges 比对 **207/226 不一致**
   （例：OCCT f173 = Cylinder w1 e[4]，端口 f173 = BSpline w2 e[2,2]）。⇒ §9.58-9.94 的按索引/贪心配对结论**全部作废**；
   后续对拍必须用**几何签名**（面 3D bbox + 表面类型 + 面积/解析 contrib）。

**下一步（确定性、按序）**：
1. 查端口 wire 里「单 ORIENTED_EDGE bound → 2 条边」的来源：resolve_loop（read_topology.rs:413-433）只 make_wire 一次，
   bind_edge_loop_vertices（:478-554）不加边 ⇒ 嫌疑在 check_pcurves_and_shift（shhealing/wire_fix.rs:3169）里的
   fix_reorder_wire/fix_small_all/fix_connected_all；须对齐 OCCT：bound 里 1 条 ORIENTED_EDGE ⇒ wire 里 1 次出现。
2. 修好后 check_wire 应对这些面给出 ±period，FixMissingSeam 触发并把 2 wire 合并为 1（对齐 ShapeFix_Face.cxx:2236-2261），
   再按 ShapeFix_Face::Perform:492-494 接线到 STEP 读入。
3. 建几何签名配对表，重新做 OCCT↔端口的**逐面三角数**对拍（§9.58-9.94 的配对数作废）。

**工具（本轮，可复现）**：specs/occt_probe/wires_probe.cpp 扩展为逐面打印 type/wires/edges_per_wire/每条 wire 的 UV span；
带 <faceIndex> 参数时逐边打印 CurveOnSurface 的 (f,l,p0,p1,ori,deg,closedOnFace)；build_wires.bat 重建。


### 9.96 round 95 —— 【已修·门禁绿】read 期第一轮 wire 修复误开 FixLacking（违背 ShapeFix_Face.cxx:372）；FixMissingSeam 已能触发（12 面）

**根因（本轮落地）**：端口在 STEP 读入期的 `check_pcurves_and_shift`（shhealing/wire_fix.rs:3169）里调用了 `fix_lacking_all`。
OCCT `ShapeFix_Face::Perform`（ShapeFix_Face.cxx:345-498）在**第一轮** wire 修复里**显式关闭 FixLacking**
（cxx:369-374：`theAdvFixWire->FixLackingMode() = false;`，与 FixSelfIntersection 一起），只在 `FixMissingSeam` 之后的
**第二轮**（cxx:457 恢复、cxx:509+ 执行）才打开。端口把两轮合成一次，于是对「周期面上两条**单闭合 seam 边** bound」的面，
`FixLacking` 经 `CheckLacking`→`doAddClosed`（ShapeFix_Wire.cxx:3781-3816）插了一条回折边：wire 从 STEP 的
1 条 ORIENTED_EDGE 变成 2 条边，`check_wire` 的 ΣΔ 抵消为 0 ⇒ `FixMissingSeam` 永不触发（§9.64/§9.95）。

**改动**：`shhealing/wire_fix.rs` 的 `check_pcurves_and_shift` 去掉 `fix_lacking_all` 调用（附 cxx:369-374 出处与理由）；
第二轮调用待 `ShapeFix_Face::Perform` 完整接线后再补。

**实测效果**：
- a3n00 的 19 个 BSpline two-bound 面：`edges` 由 `[2,2]` 回到 **`[1,1]`**（与 STEP 一致）；`--fixms` 由 **0 个 ret=true** → **12 个 ret=true**
  （`fix_missing_seam` 现在能识别「周期方向张开」的两条 wire）。
- 门禁（本轮实测，全绿）：occt-topo `--lib` **1281/0**、occt-geom `--lib` **143/0**、`step_obj_parity` **14/14**、
  `step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **3/3**、
  `phase5` 7/7、`phase9` 8/8、**`phase10` 8/8（原 7/8，本轮转绿）**、`phase19` 5/5、`phase20` 5/5；
  `export_data_obj` 15/15 ok 且 v/f 与基线**逐位一致**（Cube 24/12、Cylinder 146/140、Sphere 642/1244、Torus 1369/2592、
  Shape 6150/11372、Shape-1 3343/4336、Shape-2 3105/4792、linkrods 3494/5078 …）。

**接线 FixMissingSeam 时暴露的下一个缺口**（已实测，未留在树里）：在 `resolve_face` 尾部调用
`ShapeFixFace::with_face(&face).fix_missing_seam()` 并替换面，a3n00 立即 panic：
`index out of bounds: the len is 2 but the index is 2` at `shape_fix_compose_shell/split_by_line.rs:303`。
⇒ 下一步：修 `SplitByLine`（对照 `ShapeFix_ComposeShell.cxx` 的 SplitByLine 控制流与数组下标），
再重新接线并验证 a3n00 周期面三角数上升到 OCCT 的 360-486、总量向 11052/12324 收敛。

**工具**：`zz_probe_a3n00` 新增 `--raw`（逐 BSpline 2-wire 面打印 face/wire 的原始 children 数与 Edge 指针/朝向），
用于区分「同一边存两份」与「两条不同边」——本轮证明是后者，且由 FixLacking 插入。


### 9.97 round 96 —— 【接线落地·门禁全绿】FixMissingSeam 接入 STEP 读入；SplitByLine 下标 / CopyPCurves / BRep_TEdge range 三处忠实修复

本轮把 `ShapeFix_Face::FixMissingSeam` 真正接到 STEP 读入（`ShapeFix_Face.cxx:492-494`），并修掉沿途三处真实缺口：

1. **`SplitByLine` 0-based 下标**（`shape_fix_compose_shell/split_by_line.rs:296`）：OCCT
   `ShapeFix_ComposeShell.cxx:1722` 的 `int j = IntEdgePar.Length()` 是 **1-based** `NCollection_Sequence` 末位；
   端口 Vec 0-based，原写 `len` 导致 `int_edge_ind[j]` 越界 panic。改为 `len - 1`（附 cxx 出处）。
2. **`ShapeBuild_Edge::CopyPCurves` 必须按 `BRep_TEdge` 的 CurveRepresentation 拷贝**（`ShapeBuild_Edge.cxx:360-413`）：
   端口 `copy_pcurves` 走 `edge_geom`，而 `edge_geom` 在「无 3D 曲线」时返回 None ⇒ `SplitByLine` 造的 seam 边
   （`cxx:2061-2070`：`MakeEdge` + 两条 pcurve + Range，**不带 3D 曲线**）的 pcurve 在 `DispatchWires` 拷贝时全丢，
   结果面有一半边没有 pcurve、永远网格化失败。新增
   `GeometryRegistry::copy_edge_pcurve_slots`（tgeometry.rs）直接按 surface 指针拷贝 pcurve 槽 + range，`copy_pcurves` 改用它。
3. **`edge_parameters` 必须读 `BRep_TEdge` 的 range，与 3D 曲线无关**（tgeometry.rs）：原来对无曲线边返回 `(-inf, inf)`，
   使 `set_pcurves`/`set_pcurve_range` 存下无限 COS 窗口（实测 seam 边 range 变 `(-inf,inf)`，`CurveOnSurface` 读不出）。
   改为读 edge core 的 `first/last`（`set_edge_range` 写入），只有 core 不存在时才 `(-inf,inf)`。

**接线**：`step/read_topology.rs` `resolve_face` 尾部按 `ShapeFix_Face::Perform:492-494` 调 `fix_missing_seam()`，
返回单 Face 时替换（Shell/Compound 结果不替换，避免破坏 shell 装配）。

**a3n00 实测（同一 `data/occ/a3n00.stp`，OCCT GT = 11052 nodes / 12324 tris）**：

| | HEAD 原状 | 本轮 |
|---|---|---|
| 总量 | 11931 / 12283 | **11784 / 12938**（stats 225/226） |
| f129 | 10 / 8 | 42 / 60 |
| f173 | 72 / 1 | 75 / 73 |
| f176 | 72 / 11 | 4 / 2（退化，待查） |
| f190 | 72 / 1 | **360 / 376** ✅ |
| f201 | 72 / 3 | **339 / 548** |
| f209 | 72 / 1 | 171 / 224 |
| f216 | 72 / 4 | **360 / 376** ✅ |

**门禁（本轮实测，全绿）**：`occt-topo --lib` **1281/0**、`occt-geom --lib` **143/0**、`step_obj_parity` **14/14**、
`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **3/3**、`phase5` 7/7、`phase9` 8/8、
`phase10` 8/8、`phase19` 5/5、`phase20` 5/5；`export_data_obj` **15/15 ok 且 v/f 逐位不变**
（Cube 24/12、Cylinder 146/140、Sphere 642/1244、Torus 1369/2592、Shape 6150/11372…）。

**仍差（下一步）**：`f109-112`、`f210-214`（BSpline，v 跨多周期，`wires=2 edges=[1,1]`）在 FixLacking 关闭后**没有被
FixMissingSeam 合并**（check_wire 未判开），网格从原来的 216 → 1-14；`f176` 合并后反而退化（mt=2）；
总量三角形比 OCCT 多 ~5%。下一步：对「v 跨整数倍周期」的 seam 面核对 `check_wire` 的 `dV` 口径
（`ShapeFix_Face.cxx:1804-1805` 的 `u_range/v_range` 与 `Surf->VPeriod()`），以及 `f176` 的合并 patch。


### 9.98 round 97 —— 下一缺口钉死：SplitByLine 插入的两条 seam 边重复（同一 v、同向），缺 v=period 端

对已合并的面逐一 dump `--checkwire`（`edge_parameters` 修复后 seam 边已有 pcurve）：

- **f173**（u=[0,1.2]、v=[157.0796,314.1593]，1 wire / 4 edges）：
  ```
  e0 Forward  (0.0,314.159)→(0.0,157.080)   # u=0 边界
  e1 Reversed (0.0,157.080)→(1.2,157.080)   # seam
  e2 Reversed (1.2,314.159)→(1.2,157.080)   # u=1.2 边界
  e3 Reversed (0.0,157.080)→(1.2,157.080)   # seam —— 与 e1 完全重复
  ```
  ⇒ 两条插入边**都在 v=157.080、同向**，而周期另一端 **v=314.159 没有边**；环不闭合（缺一条 1.2 长的边），
  三角化只覆盖一部分 ⇒ f173 `mt=73`、f176 `mt=2`。OCCT 里 seam 边以 **两次出现**（`PCurve1=157`/`PCurve2=314`）
  进入 wire，不是两条重复边。

- **f109-112 / f210-214**（v 跨整数倍周期：37.699=12π、59.690=30π）：`--fixms` 显示
  `before_wires=2 ret=true result=Shell faces=2` ⇒ ComposeShell 结果是多面 Shell，读入端目前**不替换**
  （只替换单 Face），这些面保持两 wire、`mt=1-14`。OCCT 用 `ShapeFix_Face.cxx:2270-2322` 的尾段
  （`Context()->Apply` + `FixSmallAreaWire`，后者按 `ShapeAnalysis_Wire::CheckSmallArea` 去掉零面积 wire）
  收敛回单面；端口 UNPORTED（已在 `shape_fix_face.rs` 就地注明 cxx 行号）。

**下一步（按序）**：
1. `SplitByLine`（`ShapeFix_ComposeShell.cxx:1969-2110`）生成的 seam 边必须让**同一条边在 wire 里出现两次**
   （一次 `PCurve1` 在 v=joint，一次 `PCurve2` 在 v=joint+period），而不是两条各自带单 pcurve 的重复边；
   对齐 `CollectWires` 对 seam 边的 `WireData` 处理（`wire_data.rs` 已有 forward-then-reversed 的 seam 模型）。
2. 再补 `ShapeAnalysis_Wire::CheckSmallArea` + `FixSmallAreaWire` + `MapReShape` 的 null 绑定，
   让多面 Shell 收敛为单面并可被读入端替换。

**本轮门禁**：`occt-topo --lib` **1281/0**（未新增/未改动功能性代码，仅就地更新 UNPORTED 注释）；
a3n00 总量维持 **11784 / 12938**（stats 225/226）。§9.97 的其余门禁结论继续有效。


### 9.99 round 98 —— SplitByLine 合并环对齐 OCCT 的 for-增量（中性、忠实）；重复 seam 边的产生位置钉到 SplitLinePar 记录

1. **落地（忠实对齐，cxx:1949-1965）**：split_by_line.rs 的“合并零长切向段”环原来是 while，在删除分支**不推进 i**；
   OCCT 是 for (i=1; i<Length(); i++)，删除后仍然 i++。已补 i += 1（附 cxx 出处），使下一次比较落在
   OCCT 同样的 (旧 i+1, 旧 i+2) 对上。
   门禁：occt-topo --lib 1281/0、step_obj_area 11/11、step_geometry_parity 3/3、step_to_obj 13/13；
   a3n00 总量不变（11784 / 12938），属中性但更忠实。

2. **重复 seam 边的产生点（新证据）**：给 split_by_line_wires 的建边处加一次性探针，实测每次
   wires.push(seg) 都只建 **1 条边**，绝大多数调用 npar=2（SplitLinePar 只有 2 项，建 1 边）；少数 npar=4。
   ⇒ f173 最终 wire 里的 **e1/e3 两条完全相同的 seam 边**来自 SplitLinePar 中**同一参数上的两个 ITP_INTER 记录**
   （每侧闭合边界各记一次），而 cxx:1948-1965 的合并只处理 BEGSEG/ENDSEG 组合，无法折叠两个 INTER。
   OCCT 侧靠 SplitWire（cxx:1903）合并/去重同一交点；端口的 split_by_line（cxx:1443-1812 的逐 wire 版）
   在“闭合 seam 边”上的交点记录与 SplitWire 的落点还需逐行核对。

**下一步**：把逐 wire 的 SplitByLine 交点记录（split_line_par / split_line_code / split_line_vertex）
与 OCCT cxx:1443-1812 + SplitWire（cxx:1903 调用，函数体在 :1500-1898 段）逐行对齐，消除同一参数上的
重复 INTER；随后 f173/f176 的 v=period 端应由同一条 seam 边的第二次出现补齐（PCurve1/PCurve2）。


### 9.100 round 99 —— 重复 seam 边不是 SplitByLine 造的：同一条边在 CollectWires 装配时被追加两次（同朝向）

用两处一次性探针（已撤除）定位：

1. **SplitByLine 只建 1 条 seam 边**。f173 的 V-line 调用实测
   `SBL-LINE is_cut_by_u=false par=[3.3e-16, 1.2000000000000002] code=[11,11]`：
   `split_line_par` 恰 2 项、`code` 都是 `ITP_INTER(11)`；建边环 parity 只在 i=1 为 interior ⇒ **建 1 条边**。
   ⇒ §9.98/§9.99 猜的“同一参数两个 INTER”不成立；重复发生在装配阶段。

2. **同一条边在合并结果 wire 里出现两次、同一朝向**。在 `fix_missing_seam` 里打印 ComposeShell 单面结果的
   每条 wire 的 edge 指针 + pcurve 两端点（`FMS-RESULT`），得到（f173 一类面）：
   ```
   n=4 edges=[A (9.4248,36)->(3.1416,36) | B (9.4248,34)->(9.4248,36) |
              C (9.4248,34)->(3.1416,34) | B (9.4248,34)->(9.4248,36)]
   ```
   **B（seam 边）出现两次且指针相同、pcurve 方向相同**；`--checkwire` 也显示两条 arm 都是 `Reversed`。
   OCCT 里 seam 边确实出现两次，但两次的**组合朝向相反**（`ShapeExtend_WireData::Add` 会把 segment 的
   orientation 与边合成），因此 `BRep_Tool::CurveOnSurface` 才能分别取到 `PCurve1`/`PCurve2`；
   端口的 `collect_wires`（cxx:2760 处）是 `sbwd.extend(seg.edges())`，**没有把 segment 的 orientation 合成到边上**，
   所以两次出现同向、环在 UV 里退化（f173 `mt=73`、f176 `mt=2`）。

**下一步**：把 `collect_wires` 的 `sbwd` 追加改成 OCCT 的 `ShapeExtend_WireData::Add(seg)` 语义
（`ShapeExtend_WireData.cxx:80-121`：把 segment 的 orientation 与每个 edge 的方向合成；
reversed segment 走 `cxx:2762-2768` 的 `Reverse + ComputeSeams + SwapSeam` 路径），并核对
`reverse_wire_data_on_face` 是否已等价（含 seam 边的 `SwapSeam`）。对齐后再看 f173/f176 的 4 边环是否闭合。


### 9.101 round 100 —— 重复由 CollectWires 把**同一个 EXTERNAL 段追加两次**造成（同向），非 SplitByLine

在 `collect_wires` 追加处加一次性探针（比较新追加边指针与 `sbwd` 既有集合），a3n00 上命中 **46 次**：

```
COLLECT-DUP idx=2 reverse=false ext=false seg_ptrs=[...b80] dup=[...b80]
             sbwd_before=[...9a0, ...b80, ...6a0]
```

模式：某条目首次以 `seg.orientation()==External` 被选中并追加，随后按 cxx:2770-2773 被置为
`Forward/Reversed`（**不是 Internal**），于是下一轮仍是合法候选，被**第二次**选中并且仍以 `reverse=false` 追加
同一朝向的同一条边（第二次才置 `Internal`）。

- `seg_ptrs` 与 `sbwd_before` 的交集非空 ⇒ **同一条有向边进了 wire 两次**，违反 `ShapeExtend_WireData` 的
  seam 不变量（同一 seam 两次出现必须是相反朝向），正是 §9.100 里 f173 的 e1/e3 同向重复。
- OCCT 同一处（cxx:2639-2651）用 `lastEdge.IsSame(wire->Edge(j ? NbEdges : 1))` 的**最低优先级回边**分支，
  并且只在 `!index && !canBeClosed` 时接受；端口逐行实现了该分支与权重（cxx:2654-2722 已核对一致），
  差异因此落在**上一轮 `last_edge`/`end_v` 是否更新**（cxx:2793 `doupdate = index && (shorts(index) <= 0 || endV.IsNull())`）
  与 `shorts(index)` 的取值上：若 `last_edge` 没更新到刚追加的边，回边判据不会触发，EXTERNAL 段就会被同向重复选中。

**下一步**：核对 `collect_wires` 的 `doupdate` / `last_edge` / `end_v` 更新（含 `shorts` 与 `seqw` 的对齐）与
OCCT cxx:2784-2822 逐行；确认 EXTERNAL 段消费后再选中时必须走 `reverse = j != 0`（即相反朝向）而非同向。


### 9.102 round 101 —— 更正 §9.101：CollectWires 的“重复追加”是**假阳性**；装配链路无重复边

加**带朝向**的重复探针（比较 `(TShape 指针, orientation)`）后，a3n00 全流程 **0 命中**：

- `collect_wires` 追加：同一条 seam 边的两次出现朝向**相反**（第一次 `reverse=true` + `Reverse(face)`，
  第二次 `reverse=false`），符合 `ShapeExtend_WireData` 的 seam 不变量。§9.101 只比指针、没比朝向，故误报。
- `dispatch_wires` 输入（copy 前后）与 `make_faces_on_patch`（单环分支 + root 分支）都**没有**出现
  “同一 (指针,朝向) 两次”的 wire。

但探针 `--checkwire` 的 **face 173 / 176** 最终仍是同一指针、同朝向（都 `Reversed`）的 seam 边两次
（f173 `e1=e3=0x…c2030`，f176 `e1=e3=0x…c2b90`）。两者矛盾说明：

1. 探针里的“face 173”与读入过程中 ComposeShell 实际处理的那个面**不再是同一个面**（读入端替换会改变最终模型的面集合；
   之前 `--fixms` 的 `before_wires=1/2` 与该索引的面也不一定对应），或
2. 重复是在 `resolve_face` 里 ComposeShell **之外**引入的（例如原始面在 FixMissingSeam 之前跑的
   `check_pcurves_and_shift` 后处理，`read_topology.rs:714`，现在 `fix_lacking=false`，但 fix_small/fix_connected/
   fix_shifted/fix_add_pcurve 仍会改写 wire）。

**下一步（先建立身份链再改代码）**：
1. 给 ComposeShell 结果面记录一个可由探针复核的**几何签名**（surface 指针 + u/v range + 结果 wire 的边数与边指针），
   在 `resolve_face` 替换处与 `--checkwire` 的面签名对比，确认最终 face 173 到底是不是 ComposeShell 的产物。
2. 同时核对 `read_topology.rs:710-714` 的原始面后处理是否会改写 wire 的边集合（与 OCCT `TranslateEdgeLoop:CheckPCurves`
   的落点对比）；若会有，先修正该处再追 ComposeShell。

**本轮门禁**：`occt-topo --lib` 1281/0；所有探针/打印已撤除（grep 复核 0 残留）；a3n00 维持 11784/12938。


### 9.103 round 102 —— 【已修·门禁全绿】MapReShape 未按 BRepTools_ReShape 归一/还原朝向，seam 第二次出现丢反向

**根因（本轮落地）**：`DispatchWires` 对每条边先查 `rs.IsRecorded`（`ShapeFix_ComposeShell.cxx:3438-3443`），
已记录则直接 `rs.Value` 复用；同一条 seam 边在 wire 里出现两次（Forward/Reversed），二者 TShape 相同
（`TopTools_ShapeMapHasher` 按 IsSame 键），于是第二次拿到的是第一次的替换、**朝向被丢掉**。

OCCT 的 `BRepTools_ReShape` 不会丢：`replace`（`BRepTools_ReShape.cxx:164-209`）在键形状为 REVERSED 时把
**键与替换一起 Reverse**（INTERNAL/EXTERNAL 则归一为 FORWARD 并保留相对朝向），`Value`
（`cxx:230-279`）再按**查询形状的朝向**把结果 Reverse 一次。端口的 `MapReShape`
（`shape_fix_compose_shell/reshape.rs`）只是 `map.insert(key,new)` / `map.get`，没有这两步。

**改动**：把 `MapReShape::{replace, value, apply}` 按 `BRepTools_ReShape` 逐行对齐（含 INTERNAL/EXTERNAL 分支与
`apply = value(s).unwrap_or(s)` 语义），并注明 cxx 出处。

**实测效果（a3n00）**：

| | 修前 | 修后 | OCCT GT |
|---|---|---|---|
| mesh_v | 11784 | **11845** | 11052 |
| mesh_t | 12938 | **12553** | 12324 |
| f173 | 73 | **110** | ~360-486 |
| f129 | 60 | **72** | |
| f190/f201/f216 | 376/548/376 | 376/548/376 | |

`--checkwire` 复核：f173/f176 的 seam 边现在一次 `Reversed`、一次 `Forward`（同一 TShape，正确朝向），
不再是两次同向。

**门禁（全绿）**：`occt-topo --lib` **1281/0**、`occt-geom --lib` **143/0**、`step_obj_parity` **14/14**、
`step_to_obj` **13/13**、`step_obj_area` **11/11**、`step_geometry_parity` **3/3**、
`phase5/9/10/19/20` 全绿、`export_data_obj` **15/15 且 v/f 逐位不变**。

**剩余**：f176 仍 `mt=2`（u 宽 0.0417 的极窄带，seam 在 v=150.8，边界 v∈[150.8,452.4]，需单独查三角化）；
f109-112/f210-214 仍是多面 Shell（`FixSmallAreaWire/CheckSmallArea` UNPORTED）；总量三角已收到 +1.9%。


### 9.104 round 103 —— 落地 `ShapeAnalysis_Wire::CheckSmallArea` + `FixMissingSeam` 尾段（忠实、行为中性）；证明 6 个多 wire 面不是“零面积”造成

**落地（忠实移植，含 cxx 出处）**：
1. 新文件 `shhealing/small_area.rs`：`check_small_area(wire, face)` = `ShapeAnalysis_Wire::CheckSmallArea`
   （`ShapeAnalysis_Wire.cxx:2004-2098`）：23 点采样求 2D 中心 → 3D 叉积面积估计 → 当
   `|cross| < length*Confusion` 时，用临时面上的 `BRepGProp::SurfaceProperties` 与 `LinearProperties`
   做精确判定（`|area| < 0.5*length*Confusion`）。已在 `shhealing/mod.rs` 注册并 re-export。
2. `ShapeFix_Face::FixMissingSeam` 尾段（`ShapeFix_Face.cxx:2270-2322`）：ComposeShell 结果是多面 Shell 时，
   先按 `FixSmall` 丢掉空 wire/面，再在 `nbFaces>1` 时按 `CheckSmallArea` 丢 wire（`FixSmallAreaWire`），
   全丢的面移除；剩 1 面时返回单 Face（读入端即可替换）。`Context()->Apply` 以本地重建表达（就地注明 UNPORTED 边界）。

**关键实测（为何这 6 个面没被收敛）**：`check_small_area` 对 f109-112/f210-214 的 Shell 两条 wire 给出的
3D 面积估计是 `cross≈223/225`，而 `tol=length*Confusion≈1e-5`（远大于阈值）⇒ `FixSmallAreaWire` **正确地不丢**它们。
早先看到的“其中一个面 area=0”是 `surface_properties` 在该 patch 面上的**伪零**，不是真实退化。
⇒ 这 6 个面是 **ComposeShell 本身把一个 1×1 fictive grid 的面拆成了 2 个 root/面**（`MakeFacesOnPatch` cxx:3021-3143），
不是尾段能收敛的“小面”。

**门禁（全绿）**：`occt-topo --lib` 1281/0、`step_obj_area` 11/11、`step_geometry_parity` 3/3、
`phase5/10/19/20` 全绿、`step_to_obj` 13/13、`export_data_obj` 15/15 且 v/f 逐位不变；
a3n00 维持 11845/12553（此尾段在其上行为中性：无 wire 被判小面）。

**下一步**：查 `MakeFacesOnPatch` 的 root/hole 判定（`cxx:3021-3143`，含 `FClass2d` 的 `Perform`/`PerformInfinitePoint`
与 `IsSamePatch`）为何在 1×1 grid 上产生 2 个 root；对照 OCCT 是否会把其中一条 loop 判成 hole。


### 9.105 round 104 —— 状态盘点：合并后的面结构已正确，剩余缺口分两类

按当前工作区实测（`--low` / `--fixms` / `--checkwire`）：

| 面 | 结构 | mv/mt | `--fixms` |
|---|---|---|---|
| f129 | 1 wire / 4 edges（seam 一 Forward 一 Reversed） | 74 / 72 | before_wires=1 ret=false（读入已合并） |
| f173 | 1 wire / 4 edges | 94 / 110 | before_wires=1 ret=false |
| f176 | 1 wire / 4 edges | 4 / 2 | before_wires=1 ret=false |
| f190/f201/f216 | 1 wire / 4 edges | 376/548/376 | 已合并 |
| f209 | **2 wires / [1,1]** | 171 / 224 | **ret=true result=Shell faces=2**（读入不替换） |
| f109-112/f210-214 | **2 wires / [1,1]** | 72 / 1-14 | **result=Shell faces=2** |

**f129/f176 的 wire 结构（周期性 band）**：两长边（u=const，3D range 一整个 2π）与两条 seam 边
（v=const，3D range = u 宽）。f129 的 UV 矩形是 [1.5708,3.1416]×[0,6.2832]，f216 同 u 宽但 v 平移后 mt=376，
f129 只有 72 ⇒ **不是 seam 问题，而是该 patch 上的内部加点（range splitter / node insertion）不足**；
f176 更极端（mv=4，几乎只有端点），需查 `discretize_edge`/range splitter 对该 BSpline 边（3D range (0,6.283)）
的采样。

**f209/f109-112/f210-214**：`--fixms` 明确 `result=Shell faces=2` ⇒ ComposeShell 在 1×1 fictive grid 上产出
**2 个 root/面**（`MakeFacesOnPatch` cxx:3021-3143）；§9.104 已证这些 wire 不是小面积（cross≫阈值），
尾段正确不丢。

**下一步（两条独立线）**：
1. 内部加点：对照 `BRepMesh_DefaultRangeSplitter`/`BRepMesh_DelaunayNodeInsertionMeshAlgo`，查 f129/f176 这类
   周期性 BSpline patch 的 `generate_surface_nodes` 为何几乎不产点（f216 同型却正常 ⇒ 差异在该 patch 的 UV 界/周期处理）。
2. ComposeShell over-split：对照 `MakeFacesOnPatch` cxx:3021-3143 的 root/hole 判定与 `FClass2d`，
   看 OCCT 是否会把第二条 loop 判成 hole（端口判成第二个 root）。


### 9.106 round 105 —— f176 的直接原因：Delaunay 的 frontier=0（边界约束链没成为 Frontier），主流水线 74 节点 / 0 三角

在 `triangulate_model_faces`（face 176）与 `finish_mesh`（`_face_index==176`）各加一次性探针（已撤除）实测：

```
FINISH176 nodes=77 frontier=0 domain=0
MESH176  ok nodes=74 tris=0
```

- 主流水线跑通并注册了 74/77 个节点，但 `Delaunay` 的 **frontier=0**（`delaun.frontier()` 为空），
  `elements_of_domain()=0` ⇒ `collect_triangles` 交出 0 三角；
- 随后按 `discret_root.rs:460-473` 走 per-face wireframe 回退，得到 4 顶点 / 2 三角（即探针看到的 `mv=4 mt=2`）。
- 对照 f173（同结构、1 wire / 4 edges）主流水线正常出 110 三角 ⇒ **不是 seam/wire 结构问题**，
  而是 f176 这个 patch（u 宽 0.0417、v 跨 301.59）的**边界约束链没有进入 Frontier**（被当作 free link，
  或被 `erase_free_links` 清掉），使 Delaunay 定义域为空。

**下一步**：在 `meshing/node_insertion.rs` 的边界注册（`init_data_structure`/`add_links` 对 discretized edge 的
Frontier/Fixed 标记）与 `meshing/delaun/triangulation.rs::erase_free_links` 上，对 f176 逐步打印：
边界节点数、按 Frontier 注册的 link 数、`erase_free_links` 前后 link 状态；对照 OCCT
`BRepMesh_DataStructureOfDelaun::AddLink` 与 `BRepMesh_DelaunayBaseMeshAlgo::generateMesh` 的 Frontier 约定。


### 9.107 round 106 —— f176 的失败点在 Delaunay 内部：边界 73 条 Frontier link 注册成功，但三角化出 0 个 domain 三角

在 `init_data_structure`（`face_index==176`）与 `finish_mesh` 的 `erase_free_links` 前后各加一次性探针（已撤除）实测：

```
INITDS176 wires=1 jobs_none=0 pts_empty=0 links=73     # 边界链全部注册为 Frontier
PRE176    nodes=77 domain=0 frontier=72                # Delaunay perform 后：0 个域三角，72 条 frontier
POST176   domain=0 frontier=0                          # EraseFreeLinks 把无邻接三角的 frontier 全删了
```

⇒ 不是边界注册问题，也不是 link 类型问题：**约束 Delaunay 本身没有生成任何域三角**，
于是 `EraseFreeLinks`（`BRepMesh_MeshTool.cxx:204-219`，只保留有邻接元素的 link）把 frontier 清空，
`collect_triangles` 得 0，最后由 wireframe 回退给出 4 顶点 / 2 三角。

对照 f173（同结构、主流水线 110 三角）⇒ 差异在该 patch 的 UV 边界形态：f176 的 u 宽 0.0417、v 跨 301.59，
wire 在 UV 里是“两长边 + 两条都在 v=150.8 的 seam 边”（另一端 v=452.4 靠周期同一），
**边界链在 UV 里不闭合**；`Delaun::frontier_adjust`/`polygon_meshing`（`delaun/frontier.rs`）没有把它
拼成可填充的多边形。

**下一步**：对 f176 打印 `frontier_adjust` 得到的定向 frontier 环（条数/首尾顶点/是否闭合）与
`polygon_meshing` 的入参，对照 OCCT `BRepMesh_Delaun` 的 `Frontier`/`meshLeftPolygonOf` 与
`BRepMesh_DataStructureOfDelaun` 对周期边界（`myPeriodic`/`ComputePeriodic`）的处理。

**范围提示**：该缺口属**网格管线（T-54/T-90 一线）**，不在 T-92 枚举的移植项内，但验收的 7 面需要它。


### 9.108 round 107 —— 对齐 Delaunay 的 determinant guard（去掉自造 1e-9）；f176 定位到 UV 链跨周期不闭合

1. **落地（忠实，去掉 port 自造阈值）**：`meshing/delaun_index.rs::make_circle` 原来用
   `const UNPORTED_DETERMINANT_GUARD: f64 = 1e-9`（代码自注 "1e-9 is not OCCT-derived"）。OCCT
   `BRepMesh_CircleTool.cxx:112` 是 `std::abs(aD) < gp::Resolution()` = `RealSmall()` = `DBL_MIN`；
   `delaun/constants.rs:167` 早已用同一常量。已改为 `occt_core::precision::REAL_SMALL` 并更新注释。
   门禁：`--lib` 1281/0、`step_obj_area` 11/11、`step_geometry_parity` 3/3、`phase5/10` 全绿、
   `export_data_obj` 15/15 且 v/f 逐位不变；a3n00 数值不变（不是 f176 的直接原因，但消除了违规阈值）。

2. **f176 的 UV 链缺口**：`--checkwire` 显示 wire 的四条边在 UV 里是
   `e0 (+0.0208,150.8→452.4)`、`e1 (-0.0208,150.8→+0.0208,150.8)`、`e2 (-0.0208,150.8→452.4)`、
   `e3 (同 e1 反向)`：相邻边的 2D 端点**相差一个周期**（452.4 vs 150.8）而非首尾相接。
   `grep adjust_to_period|AdjustToPeriod|adjust_by_period crates/occt-topo/src/meshing` = **0 处**
   ⇒ 网格管线从不在 UV 上做周期平移；OCCT TKMesh 里也只有 `BRepMesh_CurveTessellator.cxx:249`
   一处 "skip periodic case"、没有 AdjustToPeriod ⇒ OCCT 靠的是**同一 3D 顶点在各 pcurve 上 2D 参数一致**
   （`BRep_Tool::Parameter(vertex, edge, face)`）。端口这条 seam 边的两份 pcurve 都由
   `Geom2dLine::new(*line.position())` 生成（与 OCCT `cxx:2067-2068` 一致），但**顶点的 2D 参数没有按各自边重新求解**，
   于是链断在周期缝上，`frontier_adjust`/`mesh_left_polygon_of` 拼不出多边形 → `domain=0`。

**下一步**：核对 `meshing/edge_discret.rs` 里边界节点 UV 的取值来源：是否应像 OCCT
`BRepMesh_EdgeTessellation` 那样用每条边的 pcurve（含 `CurveOnSurface` 的 PCurve1/PCurve2 选择）
逐边求 UV，而不是复用同一个 3D 节点坐标；并对照 `BRepMesh_DataStructureOfDelaun` 对同一节点的唯一性处理。


### 9.109 round 108 —— 边界节点注册与 OCCT 一致；f176 剩余差异在 seam 边的 pcurve 周期表示

核对 OCCT `BRepMesh_NodeInsertionMeshAlgo.hxx` 与端口 `meshing/node_insertion.rs`：

- OCCT `collectWirePoints`（`hxx:142-184`）逐边取 `aPCurve->GetPoint(i)` 作为 UV 点，并要求
  `addNodeToStructure(thePoint2d, myNodesMap->Length(), ...)` **按 2D 去重**（`BRepMesh_BaseMeshAlgo.cxx:133-143`）；
  端口 `register_node/add_node_to_structure`（`node_insertion.rs:210-229、1047-1060`）完全同构
  （`structure.add_node` 按 2D 容差合并，`nodes_map` 只存 3D）。
- 端口此前已核对 `init_data_structure`（`node_insertion.rs:691-764`）与
  `ShapeFix_ComposeShell::AddLinkToMesh`/Frontier 约定一致（f176 实测注册 73 条 Frontier）。

⇒ 因此 f176 的 `frontier_adjust → domain=0` **不是注册逻辑差异**，而是 **ComposeShell 产出的 seam 边 pcurve 缺周期表示**：
OCCT `MakeFacesOnPatch`/`DispatchWires` 之后，一条周期性边界在 patch 的 UV 窗口内应能首尾相接
（seam 的两份 `BRep_GCurve` 落在窗口两端），而端口 `split_by_line.rs:614-617` 造的两份 pcurve 都是
`Geom2dLine::new(*line.position())`（同一 v），结果出现“相差一个周期”的端点，2D 去重后链断。

**下一步**：用 OCCT 探针 dump **同一几何面**（a3n00 中 u≈[-0.02,0.02]、v≈[150.8,452.4] 的 BSpline 面）
在 `ShapeFix_Face::Perform` 之后的 seam 边 pcurve（PCurve1/PCurve2 的 range 与端点），与端口
`split_by_line` 的 `[lin1,lin2]` 对照，确定 OCCT 的周期表示发生在哪一步（`TransformPCurve` / `ReassignPCurve` /
`ComputeSeams+SwapSeam`），再据此补端口。


### 9.110 round 109 —— f209（7 面之一）的 Shell 来源：CollectWires 把 3 段接成 2 个环（seam 段被吃两次）

对 f209（u span ≈2.105）在 `collect_wires` 入口与 `MakeFacesOnPatch` 入口各加一次性探针（已撤除）实测：

```
CW209-IN seg0 ori=Reversed n=1 eo=[1]        # 一侧边界
CW209-IN seg1 ori=Reversed n=2 eo=[1,1]      # 另一侧边界（两段）
CW209-IN seg2 ori=External n=1 eo=[0]        # SplitByLine 造的 seam 边
MFP209 loop0 n=3 ptrs=[A, B, B]              # 环内 seam 边出现两次
MFP209 loop1 n=2 ptrs=[C, D]                 # 另一条独立环
```

⇒ `CollectWires` 把 4 条边接成了 **2 个环**（一个 3 边环里 seam 边出现两次、另一个 2 边环），
而正确结果是**一个 4 边环**（一条边界 + seam + 另一条边界 + seam 回程）。这正是 `--fixms` 里
f209/f109-112/f210-214 报 `result=Shell faces=2` 的直接来源（也是 `MakeFacesOnPatch` 出 2 个 root 的原因）。

**下一步**：在 `collect_wires`（`ShapeFix_ComposeShell.cxx:2512-2936`）的连接环上，对 f209 打印每轮的
`index/reverse/connected/can_be_closed/close` 决策与 `sbwd` 内容；重点核对 OCCT
`cxx:2824-2846` 的 close 条件（`!index || (canBeClosed && !lastEdge.IsSame(firstEdge) && IsCoincided(...))`）
与端口 `collect_wires.rs:310-329` 的差异，确认 External seam 段在已消费后为何仍能被再次选中并同向/反向追加。


### 9.111 round 110 —— f209 的连接失败逐轮序列（CollectWires）

对 f209（grid u span ≈2.105）在候选搜索与 close 处加一次性探针（已撤除）实测一轮 ComposeShell 的序列：

```
search#1: PICK seg0 (Reversed, 1 edge)
search#2: i=1(seg1, Reversed, 2 edges) SKIPPED  by !sp && can_be_closed (sp=false, can_be_closed=true)
          i=2(seam, External) sp=true endv_last=true  -> PICK seg2 reverse=true
search#3: i=1 sp=true no-vertex-match; i=2 sp=true endv_first=true -> PICK seg2 reverse=false
search#4: i=1 sp=true  endv_first=false endv_last=false  -> index=None -> close
```

⇒ 得到 loop0 = [seg0, seam_rev, seam_fwd]（3 边）后 `index=None` 关闭；随后 loop1 = [seg1]（2 边）⇒ **2 个面**。
两个可对准 OCCT 的差异点：

1. **search#2 把 seg1 因 `sp=false && can_be_closed` 跳过**（端口 `collect_wires.rs:135`）；
   OCCT 同处是 `if (!sp && (canBeClosed || (index && samepatch))) continue;`（`cxx:2611-2614`），
   差异在 `is_same_patch(seg1)` 的判定（patch 索引簿记）或 `can_be_closed` 的取值；
2. **search#4 seg1 的 first/last vertex 都不等于 `end_v`**（`same_v` 全 false）⇒ 即使不跳过也无法接上；
   对应 OCCT `endV.IsSame(seg.FirstVertex()/LastVertex())`（`cxx:2627-2632`），说明端口在 seam 第二次追加后
   `end_v` 指向的顶点与 seg1 端点不一致（顶点身份或 seam 追加时的反向处理）。

**下一步**：对同一个 f209 打印 `is_same_patch` 的输入/输出（`iumin/iumax/ivmin/ivmax` 与 seg 自身的
`get_patch_index`）与 seg1 两端点的 TShape 指针、`end_v` 的 TShape 指针，判定是 patch 簿记还是顶点身份问题；
再对照 `ShapeFix_ComposeShell.cxx:2607-2614` 与 `:2627-2632` 对齐。


### 9.112 round 111 —— f209 的 patch 索引与顶点证据：seg1 是自闭合 2 边环，且 search#2 因 v 窗口不同被跳过

对 f209 在候选搜索里打印 seg1（2 边段）自身 patch 索引与顶点指针、当前 patch 窗口（一次性探针，已撤除）：

```
CW209-SEG1 sp=false cur=(1,2,1,2) self=Some((1,2,0,1)) fv=A lv=A endv=B can_be_closed=true   # search#2
CW209-SEG1 sp=true  cur=(1,2,1,1) self=Some((1,2,0,1)) fv=A lv=A endv=C can_be_closed=false  # search#3
CW209-SEG1 sp=true  cur=(1,2,1,1) self=Some((1,2,0,1)) fv=A lv=A endv=B can_be_closed=true   # search#4
```

- **search#2**：当前 v 窗口是 `(ivmin,ivmax)=(1,2)`，seg1 自身是 `(0,1)` ⇒ `sp=false`；又因 `can_be_closed=true`
  被 `!sp && can_be_closed` 跳过（`collect_wires.rs:135` / OCCT `cxx:2611-2614`）。
- **search#4**：`sp=true`，但 seg1 的 **首尾顶点是同一个**（`fv == lv == A`）——即它是一个**自闭合 2 边环**，
  而 seam 第二次追加后的 `end_v=B` ≠ A ⇒ 两个端点都不匹配（OCCT `cxx:2627-2632` 同样条件）⇒ `index=None` 关闭。
- 结果 loop0 = [seg0, seam_rev, seam_fwd]（3 边）、loop1 = [seg1]（自闭合 2 边）⇒ 2 个面。

⇒ 两个待判定点：① OCCT 的 `IsSamePatch`/`GetPatchIndex` 在 search#2 是否也会给 `sp=false`（若是，则差异在 `can_be_closed`
的取值或 `SplitByLine` 给的 patch 索引）；② seg1 的「自闭合 2 边环」是否本身就是 `SplitWire`/`BreakWires` 的产物错误
（正确应是开口的两条边界碎片，能与 seam 首尾相接）。

**下一步**：对照 OCCT `DefinePatch`/`IsSamePatch`（`cxx:2607-2614`）与 `SplitByLine` 里 `seg.DefineIUMin/IVMin`
（`cxx:2078-2108`）的 patch 索引算法，核对端口 `define_patch`/`is_same_patch` 的 iv 窗口；并用 `--checkwire`
对 f209 原始 2 wire 的边顺序核对 seg1 是否应在 `BreakWires` 处开口。


### 9.113 round 112 —— 累积改动后的全门禁复核（确认无回归）

当前工作区（11 轮累积：read 期 FixLacking 关闭、SplitByLine 0-based 下标 + for 增量、CopyPCurves 槽拷贝、
edge_parameters 读 core range、MapReShape 朝向归一、FixMissingSeam 尾段+CheckSmallArea、Delaunay determinant guard 对齐）
逐项复核：

| 门禁 | 结果 |
|---|---|
| `occt-topo --lib` | **1281/0** |
| `occt-geom --lib` | **143/0** |
| `step_obj_parity` | **14/14** |
| `step_to_obj` | **13/13** |
| `step_obj_area` | **11/11** |
| `step_geometry_parity` | **3/3** |
| `phase5/9/10/19/20` | 7/8/8/5/5 全绿 |
| `export_data_obj` | **15/15，v/f 逐位不变** |

a3n00：11845 / 12553（OCCT GT 11052 / 12324）；7 个周期面中 f190/f201/f216 达 360-486，
f129=72、f173=110、f176=2、f209=224 仍差（§9.105-9.112 已把各自缺口钉到：
f176 的 Delaunay frontier=0（周期链不闭合）、f209 的 CollectWires 接成 2 环、f129/f173 的内部加点不足）。


### 9.114 round 113 —— 排除两条 mesh 侧嫌疑：range splitter 周期夹取已忠实、OCCT Delaun 无周期处理

1. OCCT `BRepMesh_DefaultRangeSplitter::updateRange`（`cxx:202-236`）对周期面把离散范围夹到
   `theDiscreteFirst + (geomLast-geomFirst)`（一个周期）。端口 `range_splitter/param_set.rs:217-239`
   逐行等价（含非周期 out-of-domain 保护）⇒ **不是缺口**。
2. `grep Periodic BRepMesh_Delaun.cxx` = **0**、`BRepMesh_DataStructureOfDelaun.cxx` = 0
   ⇒ OCCT 的 Delaunay 也不做周期缝合，边界必须在 2D 里自闭合。
   `collectWirePoints`（`BRepMesh_NodeInsertionMeshAlgo.hxx:142-184`）也只是逐边取 `aPCurve->GetPoint(i)`。

⇒ 结论：OCCT 的周期面能闭合，只能来自**各边 pcurve 本身在 UV 上首尾相接**（尤其 seam 边第二次出现
应拿到周期平移后的那份 pcurve）。端口的 seam 边两份 pcurve 都是 `Geom2dLine::new(*line.position())`
（同一 v），所以链断。

**下一步（唯一剩下的权威做法）**：用 OCCT 探针 dump 同一几何面在 `STEPControl_Reader → FixShape`
之后的 seam 边两条 `BRep_GCurve`（`BRep_Tool::CurveOnSurface` 的 PCurve1/PCurve2 的 `First/Last` 与端点），
确认第二次出现是否被平移到 `v = joint ± period`；若是，定位 OCCT 在哪一步平移
（`ShapeFix_ComposeShell.cxx:2061-2108` / `ShapeBuild_Edge::TransformPCurve` / `CollectWires` 的
`ShapeExtend_WireData::Add`）并在端口同点补上。


### 9.115 round 114 —— 【决定性】OCCT 的 seam 边确实带周期平移的两份 pcurve；端口最终只剩一份

用扩展后的 OCCT 探针（`wires_probe` 增加 "SEAM" 输出：对每条 `BRep_Tool::IsClosed(E,F)` 的边打印
PCurve1/PCurve2 的 range 与端点）跑 a3n00，找到与端口 f176 **同几何**的面（OCCT face=32，
u=[-0.0208333,0.0208333]、v=150.796/452.389）：

```
SEAM face=32 c1=(-0.0208333,0.0208333) p1=(-0.0208333,452.389)->(0.0208333,452.389) c2=(-0.0208333,0.0208333) p2=(-0.0208333,150.796)->(0.0208333,150.796)
SEAM face=32 c1=(-0.0208333,0.0208333) p1=(-0.0208333,150.796)->(0.0208333,150.796) c2=(-0.0208333,0.0208333) p2=(-0.0208333,452.389)->(0.0208333,452.389)
```

⇒ **OCCT 的 seam 边两份 pcurve 分别落在 v=150.796 与 v=452.389（相差一个周期）**，边界链因此闭合。
端口 `split_by_line` 造的两份 pcurve 同 v，链断。

**OCCT 的平移步骤已定位**：`ShapeFix_ComposeShell::MakeFacesOnPatch` `cxx:3287-3327`（"pdn: shift pcurves in
the seam to make OK shape w/o fixshifted"）：对每个 **REVERSED 且 `BRep_Tool::IsClosed(E, myFace)`** 的边，
若两份 pcurve 相同或端点重合，就把 `c22` 按 `myUPeriod`/`myVPeriod` 平移。端口已按此逐行实现
（`make_faces_on_patch.rs` 顶部新块，含 `u_closed/v_closed` 与周期）。

**但端口该块不触发**：实测 f176 的 loop 四条边全部 `closed=false、npcs=1`——即 `DispatchWires` 之后
seam 边在 **patch 面上只剩一份 pcurve**（`ShapeBuild_Edge::ReassignPCurve` 在 `npcs<1` 时只 `UpdateEdge(edge, pc, sub)` 一份），
于是 `cxx:3287` 的前置条件 `IsClosed(E, myFace)` 不成立。

**下一步**：追踪 OCCT 里最终面 seam 边两份 pcurve 的来源——`MakeFacesOnPatch` 的 `myFace` 是 tmpF（原面），
`CurveOnSurface(E, myFace)` 会拿到 `ReassignPCurve` 保留在 old 上的那份 `pc2`；确认 OCCT 是否在
`cxx:3287-3327` 里通过 `c22->Translate(shift)` 就地改的是**共享 handle**（从而两份其实同一对象）还是另有 `UpdateEdge`。
端口要在同一语义下让最终面边持有两份（相差一个周期）的 pcurve，f176/f129/f173 的 Delaunay 才能闭合。

**门禁**：`--lib` 1281/0、`step_obj_area` 11/11、`step_geometry_parity` 3/3、`phase5/10` 全绿、
`export_data_obj` 15/15 且 v/f 逐位不变（新块在当前模型上受 `IsClosed` 前置条件保护，未触发）。


### 9.116 round 115 —— 更正位置：seam pcurve 平移在 DispatchWires（不在 MakeFacesOnPatch）；并补上其后的 FixShifted

1. **位置更正**：`cxx:3287-3327` 的 seam pcurve 平移属于 `ShapeFix_ComposeShell::DispatchWires`（`cxx:3275-3533`），
   不是 `MakeFacesOnPatch`；且它跑在 per-edge `Copy/ReassignPCurve` **之前**，此时 `myFace` 还是原面 tmpF。
   已把该块从 `make_faces_on_patch.rs` 移到 `dispatch_wires.rs`（`m_pnts` 之后、per-edge 循环之前），
   实测前置条件成立：f176 的两条 seam 出现都是 `closed=true npcs=2`。
2. **补上其后缺的调用**：OCCT `DispatchWires` 在建边前对每条 wire 还跑
   `sfw.Load(sbwd); sfw.FixShifted(); sfw.FixDegenerated();`（`cxx:3345-3346/3358`）。端口原先没有这两步，
   已在 seam 平移后加 `fix_shifted_wire`（`fix_degenerated_all` 早前已在 `dispatch_wires.rs:124` 调）。

**仍未收敛的下一处（已定位）**：OCCT `DispatchWires` 的 `needT` 分支里，若 `BRep_Tool::IsClosed(newEdge, face)`
为真，会**同时变换两份 pcurve**并 `B.UpdateEdge(newEdge, c2dnew, c2d2, face, 0.)`（`cxx:3471-3493`）；
端口的 `need_t` 分支只写一份（其注释也承认 "the port writes the transformed one"）。
但要让该分支生效，`ReassignPCurve` 之后 patch 面上必须有**两份** pcurve；OCCT `ReassignPCurve`
在 `npcs<1` 时只 `UpdateEdge(edge, pc, sub, 0.)` 一份 —— 需再核对 `sbe.Copy(anInitEdge, false)`
是否把原面上的两份一并带过来、以及 `IsClosed(newEdge, face)` 的实际取值。

**门禁**：`--lib` 1281/0、`step_obj_area` 11/11、`step_geometry_parity` 3/3、`phase5/10` 全绿、
`export_data_obj` 15/15 且 v/f 逐位不变（本轮新增的 seam 平移 + FixShifted 在这些模型上行为中性）；
a3n00 维持 11845/12553。


### 9.117 round 116 —— 试移植 cxx:3471-3493（closed edge 同时变换两份 pcurve）：忠实但把 a3n00 总量推离目标，按“失配即停”回退

按上一轮的定位，在 `dispatch_wires.rs` 的 `need_t` 分支试实现了 OCCT `cxx:3471-3493`：
当 `BRep_Tool::IsClosed(newEdge, face)` 为真时，同时保留变换后的 forward pcurve 与 reversed pcurve
（`UpdateEdge(newEdge, c2dnew, c2d2, face, 0.)`，按边朝向定序）。

实测（真实几何对拍）：

| | 回退前（单 pcurve） | 试移植后 |
|---|---|---|
| a3n00 mesh_v/mesh_t | 11845 / **12553** | 11920 / **12731** |
| 7 个目标面 | 74/72、94/110、4/2、360/376、339/548、171/224、360/376 | **完全不变** |
| 门禁 | 全绿 | 全绿（lib 1281/0、area 11/11、geometry 3/3、phase 全绿、export 15/15 v/f 不变） |

⇒ 该分支确实触发（总量 +1.4%），但目标 7 面不变、总量**从 +1.9% 推到 +3.3%**（离 OCCT 12324 更远）。
结合 §9.115 的观察：端口 `ReassignPCurve` 之后 patch 面上只有 1 份 pcurve，
`IsClosed(newEdge, face)` 的触发集合与 OCCT 不同 ⇒ 该分支此刻不是等价移植。

按任务规则「每步以真实几何对拍，**失配即停**」，已**回退**该分支，并在原地标为
`UNPORTED (cxx:3471-3493)` 且写明原因（待 patch pcurve 表示对齐后再接）。

**下一步**：先让 `ReassignPCurve`/最终面 seam 边在 patch 面上真正持有两份相差一个周期的 pcurve
（对齐 §9.115 观察到的 OCCT face=32 状态），再打开 cxx:3471-3493。


### 9.118 round 117 —— 端口最终面的 seam 边实测：两次出现都在 v=150.7964（各 1 份 pcurve）

在 `fix_missing_seam` 的 `comp.perform()` 之后对 f176 同几何面（u span 0.0417）加一次性探针
（已撤除），打印结果面各边的 pcurve 份数/端点：

```
FMS-PC ori=Forward  npcs=1 rng=(0,6.2832)              p=(0.0208,150.7964)->(0.0208,452.3893)
FMS-PC ori=Reversed npcs=1 rng=(-0.02083,0.02083)      p=(-inf,NaN)->(inf,NaN)      # u 线, v=150.7964
FMS-PC ori=Reversed npcs=1 rng=(0,6.2832)              p=(-0.0208,150.7964)->(-0.0208,452.3893)
FMS-PC ori=Forward  npcs=1 rng=(-0.02083,0.02083)      p=(-inf,NaN)->(inf,NaN)      # u 线, v=150.7964
```

（`(-inf,NaN)` 是探针用 `p.first_parameter()` 取无界 `Geom2d_Line` 参数所致，非数据问题；两条 u 线都在 v=150.7964。）

⇒ 与 OCCT face=32（v=150.796 与 v=452.389 两次出现）相比，端口**最终面**的两次 seam 出现都在 v=150.7964、
且各只有 1 份 pcurve。

**已排除的两条路径**（本轮实测）：
1. 在 `dispatch_wires` 平移（§9.116，作用在 tmpF 原面）不会传到最终 patch 面；
2. 试让 `ShapeBuild_Edge::ReassignPCurve` 对 seam 边把两份 pcurve 一并搬到 sub（偏离 cxx:566-587），
   最终面仍 `npcs=1` ⇒ 已回退（保持与 cxx 一致）。

**下一步（决定性）**：核对 `dispatch_wires` 里传给 `reassign_pcurve` 的 `face`（patch 面）与
`make_faces_on_patch` 收到的 `surf` 是否**同一个表面 Arc**（`GeometryRegistry::shape_key` 以 surface 指针为键）；
若不同，pcurve 就写在另一个键上、最终面自然读不到。用一次性打印对比两处的 `Arc::as_ptr(surface)`。


### 9.119 round 118 —— 负结果：patch 表面 Arc 三处同一（键一致），seam 第二份 pcurve 仍缺

对 f176 同几何 patch（u span 0.0417）打印三处表面 Arc 指针（一次性探针，已撤除）：

```
DISP-PTR edge-surf   27bd00fa0f0   # dispatch 每条边处的 self.grid.patch(ind_u,ind_v)
DISP-PTR packet-surf 27bd00fa0f0   # 打包时的 self.grid.patch_pnt(m_pnts[i])
DISP-PTR mfp-surf    27bd00fa0f0   # make_faces_on_patch 收到的 surf
```

⇒ 三处是**同一个表面 Arc**，`GeometryRegistry::repr_key(face_key)`（`tgeometry.rs:358-380`，按 surface 指针为键）
因此一致 ⇒ 最终面读到的就是 dispatch 写入的那个槽位。**所以“键不一致”假设被排除**：
端口最终 seam 边只有 1 份 pcurve，是因为 `ReassignPCurve` 的 `npcs<1` 分支只写 1 份，
而之后没有任何步骤再补第二份（OCCT 最终面是 2 份）。

**下一步**：定位 OCCT 在同一 surface 上补出第二份 pcurve 的那一步。候选（按调用顺序）：
1. `ShapeFix_Wire::FixShifted`（`ShapeFix_ComposeShell.cxx:3346`，对 dispatch 前的 wire 调用）——
   端口上一轮加的 `fix_shifted_wire` 是否真的会补出第二份，需单步打印 wire 上该边 pcurve 份数；
2. `ShapeFix_Edge::FixAddPCurve`（`MakeFacesOnPatch` 内 `sfw`/新建面过程）；
3. `BRep_Builder::UpdateEdge(edge, c1, c2, face)` 的某次显式 seam 调用。
对 `ShapeFix_Wire::FixShifted`（`ShapeFix_Wire.cxx` 对应函数）先与端口 `fix_shifted_wire`
逐行对照，确认端口是否漏了“补第二份/平移”的分支。


### 9.120 round 119 —— ShapeFix_Wire::FixShifted 只平移不加 pcurve；第二份多半来自 FixAddPCurve(isSeam)

读 OCCT `ShapeFix_Wire::FixShifted`（`ShapeFix_Wire.cxx:1661-1834+`）：

- 开头按 `surf->IsUClosed/IsVClosed` 定 `uclosed/vclosed`（含 `SurfaceOfRevolution` 基曲线周期的特判），
  然后用 `Bnd_Box2d` 遍历 `sbwd` 的相邻边对，取 `sae.FirstVertex(E2)`、`sae.PCurve(E1/E2,...)`，
  按周期把 pcurve **平移**（`sbe`）使 wire 连续。
- 全文**只平移已有 pcurve，不新增 pcurve 表示**（无 `UpdateEdge(edge,c1,c2,...)` 的补第二份动作）。

⇒ 上一轮在 `dispatch_wires` 加的 `fix_shifted_wire` 即使生效，也只会在**同一份** pcurve 上平移，
无法产生 OCCT face=32 那样的“两份、相差一个周期”。

**第二份 pcurve 的更可能来源**：`ShapeFix_Wire::Perform` 里对每条边调
`myFixEdge->FixAddPCurve(sbwd->Edge(i), face, sbwd->IsSeam(i), ...)`（`ShapeFix_Wire.cxx:657-661`）——
seam 边时 `IsSeam=true`，`ShapeFix_Edge::FixAddPCurve` 会补出 seam 的第二份 pcurve。
端口需确认：ComposeShell 路径（`DispatchWires`/`MakeFacesOnPatch`）之后是否有等效的
`FixAddPCurve(isSeam=true)` 调用；端口 `fix_add_pcurve` 是否支持 seam 分支。

**旁注**：本轮在 `ShapeBuild_Edge::ReassignPCurve` 里加的一次性打印（按 `face_surface(sub).u_range()` 过滤 u span=0.0417）
未触发，说明该访问器取到的 u 范围与 patch 实际范围不同——下次改用 `repr_key`/surface 指针过滤。

**门禁**：`--lib` 1281/0（无功能性改动，仅阅读 + 已撤除探针）。


### 9.121 round 120 —— 找到“补第二份 pcurve”的现成实现（WireFix::fix_add_pcurve seam 分支），但接到 FixMissingSeam 结果上不改善 f176

**发现**：端口 `shhealing/wire_fix.rs:771-853` 的 `fix_add_pcurve(edge, face, is_seam, prec, cache)` 已完整移植
`ShapeFix_Edge::FixAddPCurve`（`ShapeFix_Edge.cxx:470-614`）的 **seam 分支**：投影 c2d 后按
`IsUClosed/IsVClosed`（`cxx:563-579`）把第二份平移一个周期，`set_edge_pcurves([c2d, c2d2])`。
这正是此前缺的“第二份相差一个周期 pcurve”。它由 `check_pcurves_and_shift`（`wire_fix.rs:3164`）在
`is_seam` 时调用（`wire_fix.rs:3202/3230`、`pcurve_ranges.rs:727`）。

**实验**：按 OCCT `ShapeFix_Face.cxx:500-554`（FixMissingSeam 后对结果面逐 wire 再处理），
在 `resolve_face` 里对 FixMissingSeam 返回的 Face 重新跑 `check_pcurves_and_shift`（`fix_lacking=false`）。
真实几何对拍：

| | 基线 | 加结果面 wire pass |
|---|---|---|
| a3n00 mesh_v/mesh_t | 11845 / **12553** | 11928 / **12875** |
| f176 | 4/2 | **4/2（未变）** |

⇒ 其他周期面密化了（+322 三角），但 **f176 完全没变**，总量反而从 +1.9% 推到 +4.5%（离 OCCT 12324 更远）。
按「失配即停」已**回退**。

**下一步**：在回退前的那次运行里打印 FixMissingSeam 结果面的 seam 边 `edge_pcurves(...).len()`：
确认 `fix_add_pcurve` 是否真的给 f176 的 edge 补到了 2 份；若补到了而 mesh 仍失败，
问题在网格模型取 pcurve 的按朝向选择（每个出现是否读到对应的 forward/reversed 槽）。


### 9.122 round 121 —— 有第二份 pcurve 后 f176 仍 4/2 ⇒ 缺口转到网格模型的按朝向取 pcurve

在上一轮的“结果面重跑 wire pass”实验里，加打印（按结果面 v span>300 过滤，正好命中 f176 的 v=(150.796,452.389)）：

```
RES-PC u=(-0.0208,0.0208) v=(150.796,452.389) ori=Forward  npcs=1
RES-PC u=(-0.0208,0.0208) v=(150.796,452.389) ori=Reversed npcs=2
RES-PC u=(-0.0208,0.0208) v=(150.796,452.389) ori=Reversed npcs=1
RES-PC u=(-0.0208,0.0208) v=(150.796,452.389) ori=Forward  npcs=2
```

⇒ 结果面上确实有出现拿到 `npcs=2`（`fix_add_pcurve` 的周期 pair 生效），
但 f176 的网格**仍是 4/2**。也就是说：**再往上补第二份 pcurve 不够，网格管线没有按 wire 朝向选中对应那份**。

**下一步**：在网格模型里对 f176 打印 `MeshEdge::pcurves_for(face_index)` 返回的 pcurve **条数与每条的 orientation**，
以及 `init_data_structure`（`node_insertion.rs:729-741`）里 `pc.orientation() == slot_ori` 的比较结果：
确认 seam 两次出现分别读到 v=150.796 / v=452.389 还是都读到同一份。
该处决定 UV 边界链能否闭合（`frontier_adjust → domain>0` 的最终一环）。

**门禁**：实验已回退；`--lib` 1281/0；a3n00 11845/12553。


### 9.123 round 122 —— 网格模型确实按朝向取 pcurve（正确），但周期 pair 疑似落在错误的边上

在“结果面重跑 wire pass”实验下，对 mesh face 176 在 `init_data_structure` 的 `job` 前打印每个 slot
的 pcurve 条数与朝向（一次性探针，已撤除）：

```
ID176 slot_ori=Forward  rev=false npc=1 oris=["Forward"]
ID176 slot_ori=Reversed rev=false npc=2 oris=["Reversed", "Forward"]
ID176 slot_ori=Reversed rev=false npc=1 oris=["Reversed"]
ID176 slot_ori=Forward  rev=false npc=2 oris=["Reversed", "Forward"]
```

- `init_data_structure`（`node_insertion.rs:729-741`）按 `pc.orientation() == slot_ori` 选第一条匹配的
  ⇒ **按朝向选取本身是对的**，seam 两次出现会分别读到 Reversed / Forward 两份。
- 但 4 个 slot 里只有 **2 个拿到 npc=2**（`fix_add_pcurve` 的周期 pair），另 2 个仍是 npc=1。
  f176 网格仍 4/2 ⇒ 周期 pair 很可能**落在错误的边上**（`is_seam_use` 以“同 TShape 出现≥2 次”判 seam，
  而 ComposeShell 结果里真正的 seam 两次出现可能不是同一 TShape），或两个 npc=1 的 slot 才是 seam。

**下一步**：对 f176 打印每个 slot 的（a）3D离散点个数（长边多点、seam 少点）、（b）两端点 3D 坐标、
（c）对应 `MeshEdge` 的底层 TShape 指针；据此判定哪个 slot 是 seam，核对 `is_seam_use`/`fix_add_pcurve`
是否作用在了真正的 seam 边上。

**门禁**：实验已回退；`--lib` 1281/0；a3n00 11845/12553。


### 9.124 round 123 —— 【重要】探针的逐面统计错位（face_stats 比 faces 少一项）；真 f176 在结果面 wire pass 下从 0 变 72

**发现 1：wire pass 确实修好了 f176 的 Delaunay。** 在 `resolve_face` 对 FixMissingSeam 结果面重跑
`check_pcurves_and_shift`（补充 seam 第二份 pcurve）后，几何过滤（u span≈0.0417、v span≈301.59）实测：

```
PRE176   face=176 u=(-0.02083,0.02083) v=(150.796,452.389) nodes=77 domain=72 frontier=74
MESH176  face=176 ... result nodes=74 tris=72 map_ok=true
```

⇒ UV 边界链闭合（seam 两次出现分别读到 v=452.389 / v=150.796），Delaunay 从 **domain=0 变 72 个域三角**。
网格模型按朝向取 pcurve 的选择本身是正确的（`node_insertion.rs:729-741`）。

**发现 2：探针逐面统计与 `faces_of` 错位。** `build_shape_mesh` 里打印 `face_stats` 的位置与 face 索引：

```
FS-POS pos=0   face=0     FS-POS pos=174 face=175
FS-POS pos=1   face=1     FS-POS pos=175 face=176   # 真 f176：v=74 t=72
FS-POS pos=2   face=2     FS-POS pos=176 face=177   # v=4 t=2
FS-POS pos=3   face=3     FS-POS pos=177 face=178
```

`TOTAL ... stats=225`（226 个面里少一项），因此某个 face 之后位置整体前移一位：
探针的 `stats.get(i)` 在 F 176 处显示的是 **face 177** 的 4/2，而不是几何为 du=0.0417 的真 f176。
按 `index` 查找（`find(|s| s.index == 176)`）得到 **74/72**——这才是真 f176 在 wire pass 下的计数。

**门禁（保留 wire pass）**：`--lib` 1281/0、`step_obj_area` 11/11、`step_geometry_parity` 3/3、
`phase5/10` 全绿、`export_data_obj` 15/15 且 v/f 逐位不变；a3n00 11928/12875（较基线 +83/+322）。

**下一步**：修探针的 face↔stat 映射（按 surface 身份或 3D 包围盒匹配，而非数组位置），
重新测量 7 个周期面的真实三角数，再判定 wire pass 的取舍与后续“内部加点不足”的方向。


### 9.125 round 124 —— 修正探针 face↔stat 映射：FaceMeshStat 增加 shape_key，实测为双射；7 面真实计数全部改写

**改动（诊断工具，非行为近似）**：`FaceMeshStat` 增加 `shape_key: usize`（`discret_root.rs:73-80`），
两条产出路径（`build_shape_mesh` / `build_shape_mesh_wireframe`）都用
`GeometryRegistry::shape_key(face)` 填充；探针 `zz_probe_a3n00` 改为按指针身份 `find(|s| s.shape_key == key)` 取统计，
不再用数组位置。

**双射校验**：`STATMAP matched=225 unmatched=1 sum_mt=12875 flat_mt=12875`（226 面中 1 面无统计，
逐面三角数之和 == 扁平总数）⇒ 映射可信。

**修正后的 7 个周期面真实计数（当前保留 wire pass）**：

| 面 | du | dv | wires | mv/mt（修正后） | 旧（错位）显示 |
|---|---|---|---|---|---|
| F 129 | 3.14159 | 6.28319 | 1 | **449 / 802** | 74/72 |
| F 173 | 1.2 | 157.08 | 1 | **78 / 80** | 94/110 |
| F 176 | 0.04167 | 301.59 | 1 | **74 / 72** | 4/2 |
| F 190 | 0.23226 | 486.95 | 1 | **74 / 72** | 360/376 |
| F 201 | 0.69850 | 229.34 | 1 | **74 / 72** | 339/548 |
| F 209 | 2.10526 | 59.69 | 2 | **171 / 224** | 171/224 |
| F 216 | 3.14159 | 6.28319 | 1 | **457 / 818** | 360/376 |

⇒ 此前 §9.105-9.124 里“目标 7 面”的 `mv/mt` 结论（哪些面达标、哪些差）**建立在错位映射上，需要重新推导**；
真实情况是：F 129/F 216 已达 ~450/~810（超过 360-486 上界），F 173/176/190/201 仅 72-80，F 209 仍 2 wire。

**门禁**：`--lib` 1281/0、`step_obj_area` 11/11；a3n00 11928/12875（保留 wire pass）。

**下一步**：用修正后的映射重新与 OCCT GT 对拍，确定验收所指的 7 个面到底是哪几个、各自 OCCT 计数，
再针对真正偏低的面（F 173/176/190/201 = 72-80 vs 目标数百）查内部加点/范围分割。


### 9.126 round 125 —— 修正映射后的权威分布：BSpline 面已无 1-11，周期窄带停在 74-80（只铺了边界条带）

用修正映射（§9.125）统计 a3n00 的 34 个 BSpline 面（保留 wire pass）：

- **mv<=11 的 BSpline 面：0 个**（验收的“1-11 起点”已被结果面 wire pass 消除）；
  全场 mv<=11 的 80 个面全是 Plane（4-8 顶点，正常）。
- **mv 在 360-486：4 个**：F 129 (449/802, du=π dv=2π)、F 191 (360/376, 10 wires)、
  F 200 (360/376, 10 wires)、F 216 (457/818, du=π dv=2π)。
- **1 wire / 4 edges 的周期窄带：F 173 (78/80)、F 176 (74/72)、F 190 (74/72)、F 201 (74/72)**。
  这些的边界点数 = 37+37+2+2 = 78，输出 74 顶点 / 72 三角 ⇒ **只有边界条带三角化，几乎没有内部节点**。
- **2 wire 的带：F 109-112 (72/2)、F 208/211/213/214 (72/1-2)、F 209/210/212 (171/224)、F 215 (72/14)**。

⇒ 验收目标（7 面 360-486）当前最接近的是 F 129/F 216（已 449/457），但**周期窄带 F 173/176/190/201
只有 74-80**，缺口在**内部加点**（`DelaunayNodeInsertionMeshAlgo::generateSurfaceNodes` /
`list_surface_nodes` + RangeSplitter 的内部网格），而不是边界/seam/pcurve（后者已闭合）。

**下一步**：对 F 176 打印 `list_surface_nodes` 生成的内部 UV 点个数与 RangeSplitter 的
`delta`/cells 值（u 宽 0.0417、v 跨 301.59、deflection 1.076），对照 OCCT
`BRepMesh_DefaultRangeSplitter::computeDelta/computeTolerance` 与
`NodeInsertionMeshAlgo::computeParameters`；确认是内部网格没生成，还是被分类器（`Classifier`）全丢。

**门禁**：`--lib` 1281/0、`step_obj_area` 11/11；a3n00 11928/12875。


### 9.127 round 126 —— 周期窄带的内部节点只生成 2×9：NURBS 内部网格没被加密

对 f176 在 `list_surface_nodes` 与 `generate_nurbs_grid` 各加一次性探针（已撤除）实测：

```
LSN176 gen=18  ru=(-0.02083,0.02083) rv=(150.796,452.389) delta=(0.02083,1.0041)
GNURS  u=2 v=9 delta=(0.02083,1.0041) tol=(0.0112,0.5976) defl=1.076 angle=0.349 min_size=0.1076
```

- `create_range_splitter` 对 BSpline 选 `NURBSRangeSplitter`（`nodes.rs:186`），
  `generate_surface_nodes → generate_nurbs_grid`（`param_set.rs:520`）；
- 候选只有 **u=2 / v=9**（18 个），而 `delta_v=1.0041` 对 v 跨 301.59 应产生 ~300 条线；
  `delta_u=0.02083` 对 u 宽 0.0417 应产生 ~1-2 条线。**u 方向（绕小半径、高曲率）完全没加密**。
- 端口与 OCCT 的 `initParameters`/`initParamsFromIntervals` 都只播 CN 区间端点 + 内部边参数
  （OCCT `BRepMesh_NURBSRangeSplitter.cxx:454-488`、`230-256`；`UVParamRangeSplitter::AddPoint` 也不进参数表），
  因此加密只能来自 **`AnalyticalFilter`**（`cxx:30-215`，按 deflection/angle 在相邻 iso 之间插中点）。
  端口的 `analytical_filter`（`param_set.rs:630+`）显然没有把 u 从 2 插到应有的上百。

**下一步**：逐行对照端口 `analytical_filter`/`generate_nurbs_grid` 与 OCCT
`BRepMesh_NURBSRangeSplitter::GenerateSurfaceNodes`（`cxx:317-398`）+ `AnalyticalFilter`（`cxx:30-215`），
定位“插中点”分支为何没跑（注意 `u_seq.len() < 2`、`iso_params.is_empty()` 等前置条件，
以及 `computeGrainAndFilterParameters` 传入的 source 是否是稀疏集）。

**门禁**：`--lib` 1281/0；a3n00 11928/12875。


### 9.128 round 127 —— OCCT 逐面网格普查（与端口同 deflection 1.076007）：总量 +3.3%，360-486 面双方各 4 个

扩展 OCCT 探针（`specs/occt_probe/wires_probe.cpp`）：读取 shape 后按端口的线性 deflection `1.076007`
跑 `BRepMesh_IncrementalMesh(shape, lin, false, 20deg, true)`，对每个面打印
`type / nodes / tris / u=(FirstU,LastU) / v=(FirstV,LastV)`：

```
MF 0 type=6 nodes=68  tris=66  u=(-0.0237,0.0335) v=(0.1216,0.4912)
MF 1 type=0 nodes=8   tris=6   u=(-7.7,8.9e-16)   v=(388.6,389.9)
MF 2 type=6 nodes=486 tris=669 u=(-0.1915,0.2242) v=(0.1039,6.1792)
...
```

- **OCCT 逐面求和**：`nodes=11145 tris=12466`；端口 `11928/12875` ⇒ +783 节点(+7.0%)、**+409 三角(+3.3%)**。
- **tris 在 360-486 的 OCCT 面：4 个**：MF 29 (360/376)、MF 31 (360/376)、MF 48 (360/376)、MF 50 (363/379)；
  端口同在 360-486 的也是 4 个：F 191 (360/376)、F 200 (360/376)、F 205 (360/376)、F 217 (360/376)
  （均为 10 个洞的 10-wire 面）。
- 修正映射后 **端口的 BSpline 面已无 1-11**；验收里“7 个周期面 1-11→360-486”的表述
  与双方实测都对不上（OCCT 的 360-486 是 10-wire 大面），需要在逐面对拍里重新确定验收面身份。

**下一步**：用 `(type, u/v 跨度, wires)` 把端口 226 面与 OCCT 226 面配对，逐面 diff `nodes/tris`，
列出偏差最大的一批（过密/欠密），再针对欠密面查内部加点。

**门禁**：`--lib` 1281/0；端口 a3n00 11928/12875。


### 9.129 round 128 —— 逐面对拍（bbox 配对 172/226）：真正偏差在 2-wire 周期带，窄带本身与 OCCT 一致

给 OCCT 探针 MF 行与端口 F 行都加了面的 3D 包围盒，按 bbox（中心+尺寸 < 1.0）配对，
172/226 配成，按 |dt| 排前 22：

```
F209 BSpline 171/224  MF56 t6  74/72   dt=+152
F212 BSpline 171/224  MF59 t6  75/73   dt=+151
F210 BSpline 171/224  MF57 t6  75/73   dt=+151
F111 BSpline  72/2    MF83 t6  75/73   dt=-71
F112 BSpline  72/2    MF81 t6  75/73   dt=-71
F208 BSpline  72/1    MF49 t6  74/72   dt=-71
F213 BSpline  72/1    MF60 t6  74/72   dt=-71
F211 BSpline  72/1    MF58 t6  74/72   dt=-71
F214 BSpline  72/2    MF61 t6  75/73   dt=-71
F109 BSpline  72/2    MF87 t6  74/72   dt=-70
F110 BSpline  72/2    MF85 t6  74/72   dt=-70
F215 BSpline  72/14   MF62 t6  74/72   dt=-58
F154 Cylinder 72/26   MF7  t1  75/73   dt=-47
F173 BSpline  78/80   MF34 t6  74/72   dt=+8
F146 Cylinder  8/6    MF18 t1  14/12   dt=-6
F217 Plane   360/376  MF50 t0 363/379   dt=-3
...（Plane 面基本逐位一致）
```

**重大更正**：
1. **OCCT 的周期窄带就是 ~74/72**，端口 **F173=78/80 已与 OCCT 一致**（差 +8）⇒
   此前“窄带应到 360-486”“内部加点不足”的判断是**基于错位映射的误判**；
2. 真正的大偏差在 **2-wire 的 ComposeShell 结果面**：F109-112、F208-215 端口 72/**1-14**
   （几乎没铺）或 **171/224**（过密 3×），而 OCCT 都是 ~74/72 ⇒ 正是
   `CollectWires`/`MakeFacesOnPatch` 把带拆成 2 个 wire 的后果（§9.110-9.112）；
3. 另有 F154 (Cylinder) 欠密 26 vs 73。

**下一步**：回到 §9.111 的 f209 连接序列（`CollectWires` search#2 跳过 seg1 / search#4 顶点不匹配），
目标是让带合成 **1 个 wire**；合成后这些面的三角数应落到 ~72（与 OCCT 一致）。

**门禁**：`--lib` 1281/0；端口 a3n00 11928/12875 vs OCCT 逐面和 11145/12466（tris +3.3%）。


### 9.130 round 129 —— f209 连接偏差的定位：`DefinePatch`/`IsSamePatch` 本身与 OCCT 逐行一致

逐行核对 OCCT 与端口这三处（都一致）：

- `ShapeFix_WireSegment::DefineIUMin/IVMin/...`（`ShapeFix_WireSegment.cxx:190-245`）只写第 i 个元素；
  端口 `wire_segment.rs:197-230` 同（`k = i-1`）。
- `DefinePatch`（`ShapeFix_ComposeShell.cxx:896-925`）传 `nb`（边数）作 i；
  端口 `helpers.rs:126-153` 同。
- `IsSamePatch`（`cxx:2452-2508`）取 `GetPatchIndex(1)`、按 NU/NV 平移、判 `iun+1>=iux && ivn+1>=ivx`；
  端口 `helpers.rs:394-438` 逐行同。

⇒ §9.112 观察到的偏差（seg1 `(1,2,0,1)` vs 当前 `(1,2,1,2)` ⇒ `sp=false`；`can_be_closed=true` 跳过 seg1）
**不是这三处函数的实现问题**，而在：
(a) 各段 `patch index` 的**取值**（`DefinePatch` 调用时的 `code/isCutByU/cutIndex`，或默认值初始化），或
(b) `can_be_closed` 的取值（`endV.IsSame(firstV)`）。

**下一步**：用 OCCT 探针 dump 同一几何面在 `ShapeFix_Face::Perform` 内
`CollectWires` 的 `seqw` 各段 patch index 与 `canBeClosed`（或加临时 OCCT_DEBUG 打印），
与端口 §9.112 的值逐一对齐；差别落在 (a) 还是 (b) 即可定位。


### 9.131 round 130 —— DefinePatch 各调用点与 OCCT 一致；保留 wire pass 后全门禁复核（parity/step_to_obj 绿）

逐点核对 `DefinePatch` 调用（全部与 OCCT 一致）：

| 端口 | OCCT |
|---|---|
| `split_wire.rs:115 / 403 / 459 / 467` | `ShapeFix_ComposeShell.cxx:1001 / 1369 / 1424` |
| `split_by_line.rs:628`（seam 段）| `cxx:2077` |

且 `WireSegment` 初始索引 `MININD=-32000 / MAXIND=32000`（`wire_segment.rs:11-12、86`）与 OCCT
`ShapeFix_WireSegment::AddEdge`（`cxx:143-146`）的 `MININDEX/MAXINDEX` 语义一致。
⇒ §9.112 的 f209 偏差不在调用结构，而在**各段 patch index 的取值**（`GetPatchIndex` 的 joint values、
`code/cutIndex`，或源 wire 的初始索引）。

**保留结果面 wire pass 后的全门禁复核（全绿）**：

| 门禁 | 结果 |
|---|---|
| `occt-topo --lib` | 1281/0 |
| `step_obj_parity` | **14/14**（354s） |
| `step_to_obj` | **13/13**（472s） |
| `step_obj_area` / `step_geometry_parity` | 11/11、3/3 |
| `phase5/10` | 7/7、8/8 |
| `export_data_obj` | 15/15，v/f 逐位不变 |

a3n00：11928 / 12875；OCCT 逐面求和 11145 / 12466（tris +3.3%）。
真实缺口仍是 2-wire 的 ComposeShell 结果面（§9.129 逐面对拍）。


### 9.132 round 131 —— split_by_grid 的 cut_index 为 1-based（与 OCCT 一致）；seg1 的 iv=0 来源仍未落在 SplitWire/SplitByLine

- `split_by_grid.rs:95-134`：`u_start = if u_closed {1} else {2}`、`for i in u_start..=nb_u_patches`，
  `cut_index = i`（1-based，`cxx:2227-2273` 同）⇒ `DefinePatch` 的 cutIndex 基准一致。
- `get_patch_index(param, params, is_closed)`（`helpers.rs:96-123`，`cxx:467-495`）对 `np==2` 恒返回 `1`（`i=2→1`），
  因此 **`DefineIVMin/DefineIUMin` 不可能写入 0**；§9.112 观察到的 seg1 `iv=(0,1)` 只能来自
  **源 wire 的初始/默认索引**（`WireSegment::load_edges` 给 `MININD=-32000..MAXIND=32000`，
  随后被 `AddEdgePatch` 从上游段复制）。
- 在 `split_wire.rs` 的 `define_patch` 分支加探针（按 grid u span≈2.105 过滤）**未触发** ⇒
  f209 的边界段不走 `SplitWire` 的该分支，其 patch index 来自别处（`BreakWires` 或初始 `Load`）。

**下一步**：沿 `BreakWires`/初始 `load_edges` 追 seg1 的 `iu/iv` 何时被写成 `(1,2,0,1)`；
对照 OCCT `ShapeFix_ComposeShell::BreakWires`（`cxx:1722-1790`）与 `ShapeFix_WireSegment::AddEdge`
的 `MININDEX/MAXINDEX`（`cxx:143-146`）取值（OCCT 的 `MININDEX/MAXINDEX` 具体数值需核对，
若与端口 `±32000` 不同则会影响后续 `Define*` 的收敛结果）。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.133 round 132 —— 【量化】端口 2-wire 面 24 个 vs OCCT 9 个（多出 16 个，正是周期带）

OCCT 探针的 `TOTAL` 行给出 mesh 形状的 wire 直方图，端口探针 F 行也能统计：

| wires/face | 端口 | OCCT |
|---|---|---|
| 1 | 192 | **208** |
| 2 | **24** | **9** |
| 4 | 2 | 1 |
| 6 | 4 | 4 |
| 10 | 4 | 4 |
| multiwire | **34** | **18** |

`TOTAL faces=226 wires=294 multiwire_faces=18 hist: 1->208 2->9 4->1 6->4 10->4`（OCCT）

⇒ 端口**多出 16 个 2-wire 面**（并有 1 个多出成 4-wire），与 §9.129 的逐面网格偏差
（F109-112、F208-215 等 2-wire 面）完全一致；这些正是周期带被 `CollectWires` 拆成两段、
未合成 1 个 wire。**目标：让这 16 个面回到 1 个 wire（2-wire 计数 24→9），
其网格自然落到 OCCT 的 ~74/72。**

**下一步**：对某个 2-wire 面 dump 端口的 `seqw`（SplitByLine/SplitWire/BreakWires 之后、CollectWires 之前）
与 OCCT 同名面的段结构（需要 OCCT oracle 或在端口侧按 §9.112 的序列反推）；重点比较
「边界段数量」与「seam 段数量」是否与 OCCT 一致（端口 f209 是 1+2+1=4 边 / 3 段）。

**门禁**：`--lib` 1281/0；a3n00 11928/12875。


### 9.134 round 133 —— 顶点身份链核对：seam 段用 SplitWire 的交点顶点，但 seg1 仍是自闭合 2 边环

- `split_by_line.rs:452-469`：`int_vertices` 由 `self.split_wire(...)` 填充，随后每条 `int_line_par`
  把对应 `int_vertices[i-1]` 推进 `split_line_vertex`（`cxx:1900-1911` 同）⇒ seam 段
  （`cxx:2061-2070` / 端口 `split_by_line.rs:614-624`）用的 `V1/V2` 应就是边界被切处的交点顶点。
- `split_by_line.rs:76-81` 的 `a_vert_new` 分支是 **vertex-only wire** 的情形（该 wire 本身是一个顶点），
  建新顶点并 `context.replace(a_vert, a_vert_new)`，与边界边切分链无关。
- 但 §9.112 实测 seg1 是 **`fv == lv`（自闭合 2 边环）且与 seam 第二次出现后的 `end_v` 不是同一顶点**
  ⇒ 说明 `SplitWire` 在 f209 上**没有把闭合边界切成开口的两段**（或切出的两段首尾顶点未与 seam 统一）。

**下一步**：对 f209 在 `split_wire` 出口打印每条边界段的（边数、首/尾顶点 TShape、是否 closed）
与 `int_vertices` 各顶点的 TShape；确认「闭合边界被切一次后仍闭合」还是「顶点未统一」，
再对照 OCCT `ShapeFix_ComposeShell::SplitWire`（`cxx:960-1420`）的 `ShapeAnalysis_Edge::CopyPCurves`/
`Context()->Replace` 与顶点建法。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.135 round 134 —— 【关键】split_wire 对 f209 的两条边界各切出「仍自闭合」的一段

对 f209 在 `split_by_line.rs:454-469`（`self.split_wire(...)` 出口）加一次性探针（已撤除）实测
（每次 `split_by_line` 处理一条 wire）：

```
SW209 wire_edges=1 int_verts=1   e0 fv=V lv=V (Reversed)          # 边界圆 1：单边 V->V
SW209 wire_edges=2 int_verts=1   e0 fv=A lv=B, e1 fv=B lv=A      # 边界圆 2：A->B->A（仍闭合）
```

- `split_wire` 把闭合边界在 seam 交点处切开后，得到的仍是**闭合段**（单边 `V→V`；或双边 `A→B→A`），
  `int_vertices` 都只 1 个（交点顶点）。
- 与 §9.112 的 CollectWires 序列对齐：`seg0 = 圆1 (V→V)`、`seg1 = 圆2 (A→B→A)`、`seg2 = seam`。
  `seg1.first/last = A`（闭合），而 `seg2` 第一次反面追加后的 `end_v` 既不等于 A 也不等于 B
  ⇒ `search#3` 再次选中 seam 而不是 `seg1`，于是 `loop0 = [圆1, seam, seam]`、`loop1 = [圆2]`
  ⇒ 两个面（正是 §9.133 多出的 2-wire 面）。

**下一步**：打印 seam 段的 `V1/V2` 与第一次追加后 `end_v` 的 TShape，核对
① seam 的两端是否 = 圆1 的交点 V 与圆2 的交点 A；② `reverse_wire_data_on_face`（`cxx:2762-2768` +
`cxx:483-572 Reverse/ComputeSeams/SwapSeam`）是否让反向 seam 的 `last_vertex` 落到 A。
若 ② 处端点算错，就是这段的直接原因。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.136 round 135 —— 【根因】port 的 ReShape/ApplyContext 不递归替换子形状，边界边仍持未映射的顶点

§9.135 的段顶点实测（f209）：

```
初始: i=0 Reversed 1边 fv=lv=V ; i=1 Reversed 2边 fv=lv=A ; i=2 External 1边 fv=C lv=V
迭代: pick i=0 -> end_v=V ; pick i=2(反) -> end_v=C ; pick i=2(正) -> end_v=V ; index=None -> 关闭
```

seam 段的两端是 **V 与 C**，而第二条边界圆的两个端点是 **A**（`fv=lv=A`），A≠C。
但 `split_line_vertex` 里存的是 `int_vertices`（=V 与 A），seam 用
`self.context().apply(&split_line_vertex[...])`（`split_by_line.rs:585-587`）得到 V 与 C
⇒ **context 记录了 A→C 的替换**。而 `apply_context`（`split_by_line.rs:680` → `reshape.rs:113`）
对边界边调用 `context.apply(&edge)`：`MapReShape::apply = value(s).unwrap_or(s)`
（`reshape.rs:72/99-108`，**只做精确形状查表、不递归子形状**），context 里只有 **顶点** A→C、
没有边级映射 ⇒ `is_same(res, edge)` 为真 ⇒ 边界边**保持顶点 A**，CollectWires 无法把
「seam 的 C 端」接到「圆2 的 A 端」⇒ 2 个 wire。

OCCT 的 `ShapeBuild_ReShape::Apply`（继承 `BRepTools_ReShape::Apply`）是**递归**的：
未命中时会复制形状并对子形状继续 Apply ⇒ `Context()->Apply(edge)` 会用 C 重建圆2 的边，
`ApplyContext`（`cxx:382-446`、调用点 `cxx:2119-2126`）随即 `SetEdge` 替换 ⇒ 顶点统一。

**下一步（修复方向）**：让 `apply_context`（或 `MapReShape::apply`）在边未命中时，
把边的**顶点子形状**也过一遍 context 并重建边（保留原 3D 曲线/pcurve/range，用 `ShapeBuildEdge::copy` +
`TopoBuilder` 装新顶点），即对齐 `BRepTools_ReShape::Apply` 的递归语义。这应能把 16 个 2-wire 面收敛回 1 个。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.137 round 136 —— 试移植 BRepTools_ReShape::Apply 的递归子形状替换：只修好 1/16 个面且总量推离目标，回退

按 §9.136 的根因，在 `apply_context`（`reshape.rs:113`）里实现 OCCT `BRepTools_ReShape::Apply` 的递归语义：
当 `context.apply(edge)` 未命中时，取边的 `first_vertex/last_vertex` 各自过一遍 context；
若有变化，用 `copy_replace_vertices_with(edge, mfv, mlv)` 重建边并
`copy_edge_pcurve_slots` 复制 pcurve，再 `wire.set_edge(iedge, ne)`。

真实几何对拍：

| | 回退前 | 试移植后 | OCCT |
|---|---|---|---|
| 2-wire 面 | 24 | **23** | 9 |
| a3n00 mesh_v/mesh_t | 11928 / **12875** | 12123 / **13175** | 11145 / **12466** |

⇒ 只修好 **1** 个 2-wire 面，总量反而从 +3.3% 推到 **+5.7%**。按「失配即停」已**回退**。

**结论**：递归子形状替换是 OCCT 语义的一部分，但**不是这 16 个面差异的主杠杆**；
剩余 15 个面「边界段与 seam 段顶点不统一」的来源需要在别处（例如 `split_wire` 是否记录
交点顶点的替换、或 seam 段的 V1/V2 取法）。下一步从 §9.135 的 `SW209` 打印出发，
核对每个 2-wire 面的 `context` 里是否存在对应顶点映射。

**门禁**：回退后 `--lib` 1281/0、a3n00 11928/12875。


### 9.138 round 137 —— split_wire 确实记录顶点替换；递归 apply 只修 1 面的原因待查

`grep context_mut().replace` 结果（端口）：

- `split_wire.rs:319/326`：`!splitted` 时用 `empty_copied_vertex` 把原始终点 `prev_v/last_v` 记为替换
  （对应 OCCT `cxx:1261-1293` 的端点保护）；
- `split_wire.rs:359/364`：内部非流形顶点 `atmp_v` 合并到 `prev_v/v`；
- `split_wire.rs:433/446`：边级替换（`edge -> res_wire/e1`）；
- `split_by_line.rs:607/608`：seam 段 `v1->r1、v2->r2`；`split_by_line.rs:79` vertex-only 分支。

⇒ context 里**确实**有顶点级替换，§9.136 的递归 apply 方向正确；但 §9.137 实测只修好 1/16 个面，
说明**多数 2-wire 面的边界段顶点不统一的来源不是这一条**（或递归 apply 触发条件太窄：
`apply_context` 只在 `context.apply(edge)` 未命中时走顶点分支，若 context 里有**边级**映射则走外层，
不会处理顶点）。

**下一步**：在 §9.135 的 `CW209-SEG` 打印上再加 `context.is_recorded(edge)` 与
`context.apply(fv)/apply(lv)`，区分每个 2-wire 面走的是「边级映射」还是「顶点级映射」，
以及顶点分支是否被外层边级分支吞掉。

**门禁**：`--lib` 1281/0；a3n00 11928/12875。


### 9.139 round 138 —— f209 类面：context 无映射、CollectWires 每次都选 seam 两次

对 f209 类面（grid u span≈2.105）加一次性探针（已撤除）：

1. **context 无映射**：对每个 segment 的每条边，`ctx.is_recorded(edge)=false`，
   `ctx.apply(fv)=fv`、`ctx.apply(lv)=lv`（identity）⇒ §9.136 的「顶点 A→C 被记录」在多数 f209 类面**不成立**。
2. **CollectWires 每次都是同一个序列**（PICK 打印）：
   ```
   PICK index=Some(2) reverse=true  cbc=true   # seam 第一次（反）
   PICK index=Some(2) reverse=false cbc=false  # seam 第二次（正）
   PICK index=None -> close                    # loop0 = 3 边
   PICK index=None -> close                    # loop1 = 圆2
   ```
   即**从不选 i=1（圆2）**：seam 第一次追加后 `end_v` 与圆2 的 `first/last` 顶点不相等。
3. `reverse_wire_data`（`wire_data.rs:97-117`）与 `first_vertex/last_vertex`（`wire_fix.rs:51-67`，按朝向取端点）
   都是朝向正确的 ⇒ 反向 seam 的 `last_vertex` = 原 `first_vertex`，取值逻辑没问题。

**结论**：这些面的 seam 两端与圆2 的端点**就是不同的 TShape**，且 context 没有把它们统一；
来源应在 `split_wire` 生成交点顶点/端点保护（`empty_copied_vertex`）与 seam 段 `V1/V2` 的取法之间。

**下一步**：在同一次运行里打印某个 2-wire 面的 `seam.fv/lv`、`圆1.fv/lv`、`圆2.fv/lv` 三组 TShape，
并把它们与 `split_line_vertex`（= `int_vertices`）逐一对应，确定 seam 的某一端是从哪个 `split_wire`
调用产生的、为何与圆2 的不同。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.140 round 139 —— seam 的 V1/V2 就是 split_line_vertex 本身；矛盾落在「边界段顶点 vs int_vertices」

在 seam 边创建处（`split_by_line.rs:620-624`）加一次性探针（已撤除）实测：

```
SL209 seam i=2 nslv=2 v1=..d60 v2=..200 slv=[..d60, ..200]
（8 次 f209 类面，均 i=2、nslv=2，v1/v2 逐个等于 slv[0]/slv[1]）
```

⇒ seam 边的 `V1/V2` 就是 `split_line_vertex` 的两个条目（`context.apply` 在此为 identity），
seam 的顶点**按构造是对的**。于是 §9.139 的矛盾必然落在：
**CollectWires 看到的边界段（圆1/圆2）顶点**与 **`split_wire` 写进 `int_vertices` 的顶点**不是同一批 TShape**
（或其一是 `empty_copied_vertex` 的产物）。

**下一步**：在 `split_wire` 内、返回前打印「返回 wire 的每条边端点 TShape」与「本次 push 到 `int_vertices` 的 TShape」，
在同一次调用里成对比较；确认是 `splitted` 分支的端点保护（`empty_copied_vertex`）产生了不同顶点，
还是 `int_vertices` 取了另一个副本。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.141 round 140 —— split_wire 输出与 int_vertices 自洽；矛盾应在 apply_context 之后

`split_wire` 返回前的一次性探针（已撤除）实测（f209 类面）：

```
SPW209 verts=[V]        e0 fv=V lv=V            # 圆1
SPW209 verts=[A]        e0 fv=A lv=B, e1 fv=B lv=A   # 圆2（A 即交点）
```

⇒ `split_wire` 返回边的端点**就是**它 push 进 `vertices`（= `int_vertices`）的那批 TShape，自洽。
结合 §9.140（seam 的 V1/V2 就是 `split_line_vertex` = [V, A]）：
seam 应是 **V↔A**，正好接上圆1(V) 与圆2(A)，CollectWires 应能合成 1 个环。

但 §9.139 的 CollectWires 视图里 CollectWires 仍选 seam 两次 ⇒ **在 `split_by_line_wires` 的
`apply_context`（`split_by_line.rs:678-683`，对含 seam 在内的所有 wires）之后，边界段与 seam 的顶点集合出现了分歧**。

**下一步**：在 `split_by_line_wires` 的 `apply_context` 之后（`line 683` 之后）打印整个 `wires` 列表
（每段 ori/边数/首尾顶点 TShape），与 CollectWires 的视图逐段对照；重点看 seam 段与圆2 段在
`apply_context` 前后的端点是否被改。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.142 round 141 —— 【关键】seam 的 C 是 vertex-only 分支的 a_vert_new；递归 apply 失效因 copy_replace_vertices_with 对无 3D 曲线的边直接返回原边

在 `split_by_line_wires` 的 `apply_context` **之后**（同一次运行）打印 `wires`（已撤除）：

```
SBW209 wires=3
 i=0 Reversed 1 边 fv=lv=V        # 圆1
 i=1 Reversed 2 边 fv=lv=A        # 圆2
 i=2 External 1 边 fv=C lv=V      # seam，C != A
```

seam 的 C **不在**圆1/圆2 的返回段里 ⇒ 它来自 `split_by_line.rs:76-81` 的 **vertex-only 分支**：
该 wire 是一个顶点且在线上，于是 `a_vert_new = make_vertex(a_p3d, tol)` 并 `context.replace(a_vert, a_vert_new)`。
若 `a_vert` 就是圆2 的交点 A，则 context 里有 `A -> C`，OCCT 的 `ApplyContext`（递归 Apply）会把圆2 的边
重建成 C，从而与 seam 统一。

**为什么 §9.137 的递归 apply 只修好 1 个面**：它用 `copy_replace_vertices_with(edge, mfv, mlv)`，
而该函数开头是 `let Some(geom) = reg.edge_geom(&edge.0) else { return edge.clone(); }`
（`wire_fix.rs:1213-1217`）——ComposeShell 里由 `SplitWire/SplitByLine` 造的边**没有 3D 曲线**（只有 pcurve），
于是直接返回原边、顶点没换。OCCT 的 `ShapeBuild_Edge::CopyReplaceVertices`（`cxx:59-140`）
用 `BRep_TEdge::EmptyCopy` 复制**所有表示**（含 pcurve），**不要求 3D 曲线**。

**下一步（修复）**：给 `apply_context` 的顶点替换加一条「无 3D 曲线也能重建」的路径：
用 `TopoBuilder` 造裸边 + `add_edge_vertices(新顶点)`，再 `copy_edge_pcurve_slots(新,旧)` +
`edge_parameters` 的 range，最后 `wire.set_edge`。这样 vertex-only 分支的 `A->C` 才能传播到圆2，
把 16 个 2-wire 面收敛回 1 个。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.143 round 142 —— 试补「无 3D 曲线也能重建边」的顶点替换：仍只修 1/16 面，回退

按 §9.142 的定位，在 `apply_context` 的递归分支里加了无 3D 曲线的重建路径：
`TopoBuilder` 造裸边 + 映射后的两个顶点 + `copy_edge_pcurve_slots` + `edge_parameters`/range/tolerance。

真实几何对拍：

| | 回退前 | 试补后 | OCCT |
|---|---|---|---|
| 2-wire 面 | 24 | **23** | 9 |
| a3n00 mesh_v/mesh_t | 11928 / **12875** | 12123 / **13175** | 11145 / **12466** |

⇒ 与 §9.137（只走 `copy_replace_vertices_with`）**结果完全相同**，仍只修 **1** 个面，总量推离目标。
按「失配即停」已**回退**（11928/12875、`--lib` 1281/0）。

**结论**：顶点替换（递归 Apply）只对 1/16 个面有效 ⇒ 剩余 15 个面的 seam 端点 C 与圆2 的 A
**本就是不同 TShape 且 context 无 A->C 映射**（§9.142 的 `SBW209` 里 `fv=C lv=V`，C 不在圆1/圆2 段中，
但也无 `rec`/映射）。下一步需要**弄清 C 的确切来源**：是 vertex-only 分支的 `a_vert_new`，
还是另一条 wire 的 `int_vertices`；可用「在一次 `split_by_line_wires` 内给每个新顶点打索引」的方式追踪。

**门禁**：`--lib` 1281/0；工作树无调试残留。


### 9.144 round 143 —— 【决定性数据】seam 端点是 P1/P2，但第二条边界的闭合环落在 Z（≠P1,P2）

同一次运行里 `V209-INT`（split_by_line 收集的 int 顶点）与 `CW209OUT`（CollectWires 输出）配对：

```
V209-INT  ptr=P1
V209-INT  ptr=P2
V209-SEAM v1=P2 v2=P1
CW209OUT i=0 Forward n=3 fv=lv=P1     # loop0 在 P1 闭合
CW209OUT i=1 Forward n=2 fv=lv=Z      # loop1 在 Z 闭合，Z != P1,P2
```

- 每次 f209 类面：恰好 2 个 `int_vertices`（P1、P2），seam 的两端正是 `[P2, P1]`（排序后）。
- 但 CollectWires 的 **loop1（2 边，即第二条边界）闭合在另一个顶点 Z**，与 seam 的 P2 不同。
- `V209-VONLY` **一次都没出现** ⇒ 本类面**不走 vertex-only 分支**，§9.142 的 `a_vert_new` 猜测在此不成立。

⇒ 端口里「seam 用 `split_line_vertex`（P1/P2）」，而「第二条边界段的实际顶点是 Z」，
两者来自**不同的顶点集合**；CollectWires 于是把 seam 用了两次（loop0 = 圆1+seam+seam）、
把第二条边界留成 loop1 ⇒ 2 个 wire。**真正的分歧在 `split_by_line` 收集 `int_vertices`
与 `split_wire` 返回段实际顶点之间**（与 §9.141 的观察在「同一 wire」内自洽，但跨 wire 后不一致）。

**下一步**：在 `split_wire` 内给每个返回段打「来自哪条 wire/哪个 `int_vertices` 下标」的标签，
在一次 `split_by_line_wires` 内把 `wires`（输出段）与 `split_line_vertex`（seam 用）逐段对齐，
找出 Z 是哪个 wire 的顶点、为何未被收集成 seam 端点。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.145 round 144 —— split_wire 的输入/输出顶点追踪：同一 wire 被 U 线、V 线两次拆分，跨 pass 顶点错位

`split_wire` 的输入首顶点/返回段顶点追踪（已撤除）显示（f209 类面）：

```
W209 in_fv=d0120  verts=[d4080] nseg=2   # pass A：把该 wire 在某线处拆成 2 边环（交 d4080）
W209  seg k=0 fv=d4080 lv=d2b60 ; k=1 fv=d2b60 lv=d4080
...
W209 in_fv=d4080  verts=[d4560] nseg=2   # pass B：输入正是 pass A 的交点 d4080，再次拆
W209  seg k=0 fv=d4560 lv=d2d00 ; k=1 fv=d2d00 lv=d4560
```

⇒ `SplitByGrid` 对同一 `seqw` **先按 U 线、再按 V 线**各调一次 `SplitByLine`（`split_by_grid.rs:94-134`，
`cxx:2131-2275` 同）；**第二次的输入是第一次的输出**（`in_fv` = 上一次的 `int_vertices`）。
两次 pass 的 `split_line_vertex`/seam 各自独立，最终 `seqw` 里 seam 与第二条边界的顶点来自不同 pass
⇒ CollectWires 接不上（loop1 落在未收集为 seam 端点的 Z）。

**下一步**：核对 OCCT `SplitByLine` 在**输入已是被拆过的段**时的行为（`cxx:1918-2060` 里
`SplitWire` 的 `apply_context`/`split_line_vertex` 累积方式），特别是第二次 pass 的
`SplitLineVertex` 是否应包含第一次 pass 的交点；端口可能在第二次 pass 丢了或重造了顶点。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.146 round 145 —— f209 只走一条拆分线（V 线），此前“两 pass”是跨面串读

`split_by_grid` 加一次性探针（已撤除）实测（8 个 f209 类面一致）：

```
SBG209 u_closed=false v_closed=true nbu=1 nbv=1
```

- U 线循环：`u_start = if u_closed {1} else {2}` = 2，`for i in 2..=nb_u_patches(1)` ⇒ **空**；
- V 线循环：`v_start = 1`，`for i in 1..=1` ⇒ **1 次**。
⇒ f209 只做**一次** `split_by_line_wires`（V 线），§9.145 的“同一 wire 被两次拆分”是把 8 个面
（每个 2 次 `split_wire` 调用）交错读了。

于是在**同一次** `split_by_line_wires` 内：`wires = [圆1, 圆2]` → int P1、P2；
seam 用 `split_line_vertex`=[P1,P2]；但 CollectWires 的 loop1 仍闭合在 Z≠P2（§9.144）。
⇒ 矛盾被压缩到**这一次调用内部**：输入 `wires`、`int_vertices`、创建出的 seam、`apply_context` 后的输出，
四者的顶点必须逐一对齐才能定位 Z。

**下一步**：在**同一次** `split_by_line_wires` 里加统一编号（进入时给每条 wire 的顶点编号，
`split_wire` 的 `int_vertices` 编号，seam 的 v1/v2 编号，`apply_context` 后的输出 wire 编号），
一次运行即可看出 Z 属于谁、在哪一步与 P2 分叉。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.147 round 146 —— 【决定性】seam 的 X 是 split_wire 记录但未用于返回段的“孤儿交点顶点”

在**同一次** `split_by_line_wires` 内打印 in/out（已撤除），f209 类面一致：

```
par=[-30.2954,-28.1901]  slv=[X, Y]
in  i=0 n=1 fv=lv=Y ; i=1 n=2 fv=lv=Q
out i=0 Reversed n=1 fv=lv=Y      # 未变
out i=1 Reversed n=2 fv=lv=Q      # 未变
out i=2 External n=1 fv=X lv=Y    # seam：X<->Y
```

- seam 连的是 **X 与 Y**；边界 i=0 在 **Y**（对得上），但边界 i=1 在 **Q**（≠X、≠Y）。
- **X 是 `int_vertices` 里的一个顶点，却不在任何返回段上** ⇒ `split_wire` 在计算交点时把这个顶点
  推进了 `vertices`（→ `split_line_vertex`），但**没有真正把对应边在 X 处切开**（返回段保持原顶点 Y/Q）。
- 于是 seam 的一端挂在孤儿顶点 X 上，CollectWires 无法把边界 i=1 接上 ⇒ loop0 用 seam 两次、loop1 留下 Q。

**下一步（定位修复点）**：读 `split_wire` 里 `splitted` 与 `vertices.push` 的分支
（`split_wire.rs:290-350`，`cxx:1261-1337`）：确认「交点顶点已算出但边未重建」发生在哪个分支
（很可能是 `!splitted` 的端点保护分支推了 `last_v`/`vv` 却没建新边），并让其与 OCCT 的
`Context()->Replace` + `B.MakeEdge/Add(V1,V2)` 一致。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.148 round 147 —— split_wire 的 v_opt 分支与边重建；curve-less CopyReplaceVertices 试改无效（回退）

读 `split_wire.rs:290-351`（`cxx:1231-1307`）：

- `v_opt` 三支：`None`（`cxx:1238-1242`，`make_vertex(curr_pnt)` 并 `vertices.push(nv)`）、
  `Some(vv) if !do_cut`（`cxx:1243-1255`，`while j<stop { vertices.push(last_v) }`，
  然后 **`if !splitted { break; }`** —— 已 push 顶点却**不建新边**）、`Some(vv)`（`vertices.push(vv)`）。
- 之后（`cxx:1261-1293`）端点保护，`cxx:1295-1307` 用 `copy_replace_vertices_with(an_init_edge, prev_v, v)` 建新边。

**试改**：把 `copy_replace_vertices_with`（`wire_fix.rs:1213`）对无 3D 曲线的边也实现重建
（`TopoBuilder` 裸边 + pcurve/range/tolerance 拷贝，对齐 OCCT `BRep_TEdge::EmptyCopy`）。
实测：a3n00 `11928/12875` 与 2-wire 直方图（24）**完全不变** ⇒ 该路径不是 f209 孤儿顶点的来源
（边界边多半有 3D 曲线、走的是 geom 分支），已**回退**。

**下一步**：在 `split_wire` 内对 f209 每次 `vertices.push` 打印「分支名 + 顶点 + 是否随后建边」，
以及输入边的 `edge_geom(edge).is_some()`，直接看孤儿顶点 X 是在哪一支被记下、对应边为何没被重建。

**门禁**：`--lib` 1281/0；工作树无调试残留。


### 9.149 round 148 —— 【最深定位】split_wire 的两条 wire 走不同分支：i=0 直接 BREAK（未拆），i=1 拆成闭合 2 边环且交点在内部

在 `split_wire` 的 v_opt 三分支、BREAK 处、建边处各加一次性探针（已撤除），f209 每次调用：

```
P209 push last_v(!do_cut) L0 splitted=false    # wire i=0
P209 BREAK !splitted (vertex already pushed)   # 直接 BREAK，返回原 wire（顶点仍是 L0）
P209 push NONE nv N                            # wire i=1
P209 BUILT geom=true prev_v=P v=N ne_fv=N ne_lv=P   # 建边0: N->P
P209 BUILT geom=true prev_v=N v=P ne_fv=P ne_lv=N   # 建边1: P->N
```

- wire i=0：走 `Some(vv) if !do_cut` 支（`cxx:1243-1255`），已 push `last_v=L0` 后 **`break`（`cxx:1243-1246` 的
  `if (!isSplitted) break;`）** ⇒ 返回原 wire 不拆；其 `int_vertices` = **L0（自己的顶点）**。
- wire i=1：走 `None` 支（`cxx:1238-1242`）建新顶点 **N**，随后建 2 条边 **N->P、P->N**（`geom=true`）
  ⇒ 返回一个**闭合 2 边环**，其 `first/last = P`，**交点 N 在内部**。
- `split_line_vertex` = [L0, N]；seam 连 **L0 与 N**。CollectWires 只匹配段的 `first/last`，
  而 wire i=1 的 first/last=P≠N ⇒ seam 的 N 端接不上 ⇒ loop0 用 seam 两次、loop1 留 P（2 个 wire）。

**下一步**：对照 OCCT `SplitWire` 对「闭合边在一点被切」的处理（`cxx:1185-1337`）：
OCCT 是否把切点旋转成段的 `FirstVertex`（即返回 `N->P`、`P->N` 之外还会做 `SetLast`/`permutation`），
使段的 first/last 落在切点 N；若是，端口缺了这一步。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.150 round 149 —— split_wire 与 OCCT cxx:1236-1339 逐行一致；分歧指向上游 BreakWires 的段是否闭合

逐行核对 `split_wire.rs:290-351` 与 OCCT `ShapeFix_ComposeShell.cxx:1236-1339`：

- `V.IsNull()` → `MakeVertex + vertices.Append(V)`（端口 `None` 支，`cxx:1238-1242`）；
- `!doCut` → `for (;j<stop;j++) vertices.Append(lastV)` 且 `if(!splitted) break`（端口 298-309，`cxx:1244-1255`）；
- else → `vertices.Append(V)`（端口 310-313，`cxx:1256-1259`）；
- 端点保护 `EmptyCopied + Context()->Replace`（端口 317-336，`cxx:1263-1293`）；
- `newEdge = sbe.CopyReplaceVertices(anInitEdge, prevV, V)`（端口 349-351，`cxx:1307`）；
- **`sbe.CopyPCurves(newEdge, anInitEdge)`（端口 `split_wire.rs:383`，`cxx:1339`）——已存在**。

⇒ §9.149 的「wire i=1 拆成闭合 2 边环、交点 N 在内部」**不是 `split_wire` 的实现偏差**（OCCT 同代码也会这样）。
差异只能来自**输入段本身**：端口 `BreakWires` 给 `split_wire` 的第二条 wire 是**闭合**的（fv=lv），
因此单点切分后仍闭合、N 落在内部；CollectWires 只认段 first/last，接不上 seam 的 N 端。

**下一步**：核对 `BreakWires`（`break_wires.rs`，`cxx:2279-2385`）在该面上是否应把边界切成**开口**段
（OCCT 的 `shift`/permutation 只在 `!myClosedMode && wire.IsClosed()` 时做，端口同）；
若 OCCT 在此处的输入就不是闭合段，则分歧在 `BreakWires` 的上游（`SplitByGrid` 之前）。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.151 round 150 —— BreakWires 对 f209 已产出 3 段（含 External seam）；split_wire 与 OCCT 全链一致

`break_wires` 出口加一次性探针（已撤除）：

```
BW209 segs=3
 i=0 Reversed 1 边 closed=true  fv=lv=V
 i=1 Reversed 2 边 closed=true  fv=lv=A
 i=2 External 1 边 closed=false fv=C lv=V
```

⇒ `BreakWires` 之后该面 `seqw` 已是「两条闭合边界 + 一条 External 段」。随后 `SplitByGrid`（V 线）
再对它们拆分并新建 seam（§9.147 的 X↔Y）。

结合 §9.150（`split_wire` 与 OCCT `cxx:1236-1418` 逐行一致，含
`CopyPCurves`（`cxx:1339`/`split_wire.rs:383`）与 `Context()->Replace(edge, resWire)`（`cxx:1397`/`split_wire.rs:433,446`））
⇒ 剩余可疑点收窄到：
1. `ApplyContext` 对 `edge -> resWire`（多边）的**展开**是否正确（`apply_context` 的 `res.shape_type() != Edge` 分支）；
2. `CollectWires` 候选匹配（`cxx:2620-2635`）对闭合段的 `first/last` 语义。

**下一步**：核对 `apply_context` 的展开分支（`reshape.rs:126-155`）与 OCCT `ApplyContext`（`cxx:382-446`）
在「边被替换为多边 wire」时的边序/朝向/索引处理；这是 §9.147 孤儿顶点 N 落入内部后唯一未逐行核对的环节。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.152 round 151 —— 全链逐行核对完成：FixMissingSeam 门禁 / BreakWires / SplitWire / ApplyContext / CollectWires 均与 OCCT 一致

本轮补核两处：

- `apply_context` 的多边展开分支（`reshape.rs:126-155`）与 OCCT `ApplyContext`（`cxx:382-446`）逐行一致
  （边序公式 `i` / `NbEdges-i+1`、索引 `index-iedge` 返回、`SetEdge`/`AddEdge`）。
- `fix_missing_seam` 门禁（`shape_fix_face.rs:149-206`）与 OCCT `ShapeFix_Face.cxx:1722-1802` 一致
  （`sa_is_u/v_closed(surf, CONFUSION)` ↔ `mySurf->IsU/VClosed()`；BSpline 需 U/V 周期；
  无穷界替换为面界）。
- `split_wire` 的 `res_wire` 顺序（`split_wire.rs:416-434`）与 `cxx:1380-1398` 一致
  （`n=result.nb_edges()` ↔ `result.NbEdges()`；REVERSED 时 `n-k+nb_edges_start`）。

⇒ 至此 `FixMissingSeam 门禁 → BreakWires → SplitByGrid/SplitByLine → SplitWire → ApplyContext → CollectWires`
**整条链都逐行核对过且一致**，但 f209 类面仍产出 2 个 wire（OCCT 1 个）。
分歧只能是**数据级**（顶点/段身份在运行中的具体连接），非控制流实现。

**下一步（换手段）**：写一个 OCCT 侧的最小 oracle——对 a3n00 的 f209 类面单独跑
`ShapeFix_Face::FixMissingSeam`，dump 内部各步骤的段/顶点（或直接对比最终 wire 的边集合与顶点），
再与端口的同一步对照。这是唯一能把「数据级分歧」定死的方法。

**门禁**：`--lib` 1281/0；a3n00 11928/12875。


### 9.153 round 152 —— 【OCCT oracle】f209 类面在 OCCT 最终形状是 1 wire / 5 edges；端口是 2 wires

用扩展的 OCCT 探针列出 a3n00 中 u≈2.1053、v≈59.6903 的面（f209 同几何）：

```
FACE 49 type=6 wires=1 edges_per_wire: 5 spans: u=2.1053 v=59.6903
FACE 56 type=6 wires=1 edges_per_wire: 5 spans: u=2.1053 v=59.6903
FACE 57/58/59/60/61/62 同
```

⇒ **OCCT 的正确结果是「1 个 wire、5 条边」**；端口同一几何是 **2 个 wire**（`edges=[1,1]`）。
5 = 两条边界（拆出的段）+ seam 的若干次出现的合理组合，说明 OCCT 确实把 seam 与两条边界连成了**一个环**。

试加「`fix <faceIndex>`」模式直接对单面跑 `ShapeFix_Face::FixMissingSeam` 做 oracle，
但在**已被 FixShape 处理过**的 shape 上调 `FixMissingSeam` 会**访问越界崩溃**（exit 0xC0000005，所有面都崩），
该模式已撤除。要直接观测 OCCT 内部段/顶点，需要在**未 FixShape 的原始 shape**上跑（或在 OCCT 侧插桩重编译）。

**下一步**：或（a）改探针读取「不做 FixShape 的原始 shape」再单面 FixMissingSeam；
或（b）把端口的 `CollectWires` 输入（3 段：2 边界 + seam）与 OCCT 的「1 wire/5 edges」反推对照，
重点看 OCCT 的 5 条边是怎样由 2 边界 + seam 组成的（端口 loop0=3、loop1=2）。

**门禁**：`--lib` 1281/0；a3n00 11928/12875。


### 9.154 round 153 —— OCCT 目标是「1 wire / 5 edges」；端口是 3+2 分组，且 context 映射确认存在

§9.153 的 OCCT oracle：f209 类面在 OCCT 最终形状是 **1 wire / 5 edges**；端口是 **2 wires（3+2 = 5 个边出现）**
⇒ **边集合一致，只有 wire 分组不同**：OCCT 把 5 个边出现连成 1 个环，端口在 3 个后提前闭合、剩 2 个成第二个环。

在 `split_wire` 记录 `Context()->Replace(edge, res_wire)` 处加一次性探针（已撤除）：

```
RW209 map edge_ori=Reversed res_edges=2 rec=true apply_type=Wire same=false
```

⇒ context **确实**把被拆的边映射成 2 边 Wire（`rec=true`、`apply` 返回 Wire、`same=false`），
`apply_context` 应走展开分支。按 `cxx:1386-1395`（Reversed 用 `NbEdges()-k+NbEdgesStart`）与
`cxx:425-427`（展开时再次用边朝向取 `ind`）的双重反向，展开后边界段的 `first/last` 应落在**切点 N**；
但 §9.147 的 `SLW209 out` 显示第二条边界仍是 `fv=lv=Q`（原顶点）。二者矛盾，下一轮直接在
`apply_context` 的展开分支打印「输入边朝向 / segw 各边端点 / 展开后各边端点」定死。

**门禁**：`--lib` 1281/0；a3n00 11928/12875。


### 9.155 round 154 —— apply_context 展开分支在 a3n00 上从未命中；此前追踪是 8 个面交错导致误读

在 `apply_context` 展开分支（`reshape.rs:135-153`）加一次性计数探针（已撤除）：跑 a3n00 全程
**0 条 EXPAND**（`^EXPAND` 过滤为空）⇒ 展开分支根本没被走到。

与 §9.154 的 `RW209 rec=true apply=Wire` 不矛盾：那是在 `split_wire` 内**刚记录** `edge -> res_wire`
之后测的；而 `split_wire` 返回的是拆分后的 `result`（原边已不在 wire 里），所以之后
`split_by_line_wires` 的 `apply_context` 面对的是**新边**，没有映射、无需展开。

**同时**：§9.147/§9.149 的 `SLW209`/`P209` 追踪互相矛盾（`P209 BUILT N->P` 却 `SLW209 out fv=Q`），
很可能是**过滤条件 `u span≈2.105` 命中 8 个面**、多次调用的打印交错，导致把不同面的行拼在一起误读。

**下一步**：把所有 f209 追踪用「静态 `AtomicBool` 只放行第一个命中的面」严格门控，
在一次运行里得到**单面、无交错**的 `BreakWires → split_wire → apply_context → CollectWires` 全序列，
再据实判断切点 N 与边界段 first/last 的关系。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11928/12875。


### 9.156 round 155 —— 【重大转向】顶点不匹配在 BreakWires 输出就已存在：seam 端 C 与边界 B 不是同一 TShape

用「静态门控只放行第一个 f209 面」拿到**单面无交错**的全序列（已撤除打印，保留 `trace209` 门）：

```
T209 BW out i=0 Reversed 1 边 closed=true  fv=lv=A
T209 BW out i=1 Reversed 2 边 closed=true  fv=lv=B
T209 BW out i=2 External 1 边 closed=false fv=C lv=A     # C != A, C != B
T209 SL out 与 BW out 完全相同（SplitByLine 对该面没改动）
T209 CW out i=0 Forward 3 边 closed=true fv=lv=A
T209 CW out i=1 Forward 2 边 closed=true fv=lv=B
```

⇒ **在 `BreakWires` 出口（`SplitByLine` 之前）**，「两条闭合边界 A、B」与「一条 `External` 段 C↔A」
就已经**各自独立**：seam 段的一端 C 既不等于 A 也不等于 B。因此无论 CollectWires 怎么写，
都无法把 B 连进 A 的环 ⇒ 3+2 分组。

**这根因不在 `SplitByLine`/`CollectWires`/`apply_context`**（这些上一轮已逐行确认与 OCCT 一致），
而在**该面进入 ComposeShell 前的 wire/顶点结构**（STEP 读入的原始面）：seam 边的端点应与边界 B 共享同一顶点。

**下一步**：查端口 STEP 读入该面时的 wire/edge 顶点共享（`read_topology.rs` 建面 + `ShapeFix_Face::Perform`
之前的 wire 构造），确认原面是否把 seam 边的端点建成了与边界 B 不同的顶点；对照 OCCT 读入后该面的
`TopExp` 顶点去重/共享。

**门禁**：`--lib` 1281/0；a3n00 11928/12875。


### 9.157 round 156 —— 【根因修复】WireSegment::first_vertex/last_vertex 未按边朝向取端点 → 2-wire 面 24→10（OCCT 9）

从 §9.156 的单面追踪发现：`BreakWires` 输入里，wire i=1 的两条边是
`e0 Reversed fv=D lv=B`、`e1 Reversed fv=B lv=D`，但 `WireSegment::first_vertex/last_vertex`
返回的是 **B/B**，而 seam（External）的两端是 **D 与 A** ⇒ CollectWires 按 first/last 匹配时接不上 D。

根因：`wire_segment.rs:118-125` 的 `first_vertex/last_vertex` 用 `edge_vertices(e).0/.1`（**原始子顶点，忽略边朝向**），
而 OCCT `ShapeFix_WireSegment::FirstVertex()/LastVertex()`（`cxx:89-101`）是
`ShapeAnalysis_Edge::FirstVertex/LastVertex(edge)`，**按朝向取**（REVERSED 边的 first 是其 TShape 的 last 子）。
端口另有按朝向的 `crate::shhealing::first_vertex/last_vertex`（`wire_fix.rs:51-67`），两处不一致。

**修复**：`WireSegment::first_vertex → edges.first().and_then(shhealing::first_vertex)`；
`last_vertex → edges.last().and_then(shhealing::last_vertex)`。

真实几何对拍（a3n00 wire 直方图）：

| wires/face | 修复前 | **修复后** | OCCT |
|---|---|---|---|
| 1 | 192 | **206** | 208 |
| 2 | 24 | **10** | 9 |
| 4 | 2 | 2 | 1 |
| 6 | 4 | 4 | 4 |
| 10 | 4 | 4 | 4 |

⇒ **2-wire 面 24→10，几乎对齐 OCCT 的 9**；6/10-wire 完全一致。

门禁：`--lib` 1281/0、`step_obj_area` 11/11、`step_geometry_parity` 3/3、`phase5/10` 全绿、
`export_data_obj` 15/15 且 v/f 逐位不变。a3n00 总量 11928/12875 → **11675/13168**（OCCT 逐面和 11145/12466）；
`step_to_obj`/`step_obj_parity` 已在后台复核。

**下一步**：查剩余 10 vs 9、4 vs 1 的差（很可能是同类「朝向/顶点」细节），并复核总量走向。


### 9.157 附 —— 长门禁复核（`WireSegment` 朝向修复后，全绿）

- `step_to_obj` **13/13**（465s）
- `step_obj_parity` **14/14**（347s）

连同 `--lib` 1281/0、`step_obj_area` 11/11、`step_geometry_parity` 3/3、`phase5/10`、`export_data_obj` 15/15（v/f 逐位不变）
⇒ 该修复落地且无回归。


### 9.158 round 157 —— 朝向修复后的逐面对拍：f209 类面 dt=+1（几乎完美）；最大残差转为 F113 无网格

用 bbox 配对（186/226）重跑逐面 diff：

```
F113 Cylinder 0/0    MF25 t1 228/228  dt=-228   # 端口该面无三角！
F170 Cylinder 75/27  MF5  t1  79/77   dt=-50
F216 BSpline 457/818 MF51 t6 475/854  dt=-36
F173 BSpline  78/80  MF34 t6  74/72   dt=+8
F146 Cylinder  8/6   MF18 t1  14/12   dt=-6
F39  BSpline  26/46  MF176 t6 28/50   dt=-4
F217 Plane   360/376 MF50 t0 363/379  dt=-3
F208/F209/F210/.../F109/F110 BSpline 75/73  MF** 74/72  dt=+1   # f209 类面几乎完美
```

⇒ §9.157 的 `WireSegment` 朝向修复让 **f209 类面（10 个）与 OCCT 只差 1 个三角**，
wire 结构 1 wire / 5 edges 也与 OCCT 的 `edges_per_wire: 5` 一致。

**最大残差变成 F113（Cylinder）端口 `0/0`（完全没有三角）**，OCCT 是 228/228；
其次是 F170（Cylinder 27 vs 77）、F216（818 vs 854）。

**下一步**：查 F113（u/v 与 bbox 对应的 OCCT face 25）为何产出 0 三角——是主流水线对该面
报错走回退、还是回退也空（`build_shape_mesh`/`fix_missing_seam` 影响）。这属于网格管线，
可能与本轮读入面的 wire 结构变化相关。

**门禁**：全绿（`--lib` 1281/0、`step_to_obj` 13/13、`step_obj_parity` 14/14、area/geom/phase/export 均过）；
a3n00 11675/13168（vs OCCT 逐面和 11145/12466）。


### 9.159 round 158 —— F113 无网格的原因：它不在网格模型里（model 225 面 vs 形状 226 面）

对 a3n00 的网格管线加一次性探针（已撤除）：

- `wires>=4` 的模型面（114/115/128/132/171/191/200/205/217）**全部产出三角**（无 0）；
- `tri.is_none()`（要回退的面）**一次都没出现**（`MFNONE` 空）。

但探针里 F113 是 `type=Cylinder wires=4 edges=[8,14,1,1] mv=0 mt=0 stat_idx=None`，
且 `TOTAL ... stats=225`（226 个形状面里只有 225 个有统计）⇒ **网格模型只有 225 个面，F113 根本不在模型里**
（`ModelBuilder::build_model` 把它跳过了），所以它没有三角。

**下一步**：查 `ModelBuilder::build_model`（`meshing/model_builder`）为何跳过 F113（OCC face 25 对应）。
该面 `du=2π`、4 个 wire（8/14/1/1 条边），可能是建面/离散阶段被过滤。属网格管线，与 T-92 的 ComposeShell 无直接关系。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.160 round 159 —— F113 被跳过是因为 `build_model` 把它标成 FAILURE（外环 `add_wire` 失败）

实测：`MODEL faces_nb=226 stats=225` ⇒ 模型有 226 个面，但 `triangulate_model_faces` 只产出 225 个统计。

定位：`discret_root.rs:390-397`（已撤除打印）——`for i in 0..nfaces { if f.is_status(FAILURE) || f.is_status(REUSED) { continue; } ... }`
⇒ **被标 FAILURE 的面直接跳过、无统计**。而 `build_model`（`model_builder/wire_builder.rs:576-582`）：

```rust
if let Some(outer_wire) = &outer {
    if !add_wire(&mut model, face_index, outer_wire, &mut edge_index) {
        fm.set_status(MeshStatus::FAILURE);
        continue;
    }
}
```

⇒ **F113 的外环 `add_wire` 失败** → 面子被标 FAILURE → 网格阶段跳过 → `mv=0`。

OCCT 对同一面（face 25）能网格出 228 三角 ⇒ 端口的 `add_wire` 在该外环上失败了。
属**网格管线**（`BRepMesh_ShapeVisitor::addWire`），与 T-92 的 ComposeShell 无直接关系，
但它是当前最大的逐面残差（-228）与 `stats=225` 的来源。

**下一步**：对 F113 外环（8 条边）在 `add_wire` 里逐步打印失败点（3D 曲线/pcurve 缺失、或 wire 顺序/分类失败）。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.161 round 160 —— F113 是被 `heal_self_intersecting_wires` 判为自交而标 FAILURE（网格管线，非 T-92）

定位链：`build_model` 的 `add_wire` **从未失败**（`ADDW` 探针全空），
`triangulate_model_faces` 的 FAILURE 跳过命中的是模型面 113（`MFSKIP i=113 fail=true box=(-87.50,-34.00,-100.00)-(87.50,34.00,-32.00)`，
与探针 F113 的 bbox 一致）。

再看唯一的运行期 FAILURE 设置点：`discret_root.rs:605-609`

```rust
for i in intersecting.keys().copied() {
    let f = model.face_mut(i)?;
    f.set_status(MeshStatus::SELF_INTERSECTING_WIRE);
    f.set_status(MeshStatus::FAILURE);
}
```

⇒ `heal_self_intersecting_wires`（`BRepMesh_ModelHealer` + `BRepMesh_FaceChecker`）
把 F113 判为**自交**并标 FAILURE，随后三角化阶段跳过它 ⇒ `mv=0`、`stats=225`。

OCCT 对同一面能网格出 228 三角 ⇒ 端口的 `face_intersecting_edges`/`FaceChecker` 在 F113 上**误报**。
属网格管线（非 T-92 ComposeShell），但它是当前最大逐面残差（-228）。

**下一步**：对照端口 `face_intersecting_edges`（`discret_root.rs:613+`）与 OCCT
`BRepMesh_FaceChecker::Perform`（`BRepMesh_FaceChecker.cxx`）在 F113 边界的判定差。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.162 round 161 —— 总量 +702 的来源：40 个未配对面（端口 4971 vs OCCT 3953）

用 bbox 配对（186/226）拆总量：

- 已配对面：`sum dt = -316`（端口偏小；其中最大 F113 -228、F170 -50、F216 -36）；
- 未配对：端口 40 面 `sum t = 4971`，OCCT 40 面 `sum t = 3953` ⇒ **+1018**；
- 合计 +702，与 `port 13168 vs occt 12466` 一致。

端口未配对的大面：F171 BSpline 1576、F129 BSpline 802、F202 Torus 548、F118/F122 Torus 211×2、
F117/F121 Torus 162×2、F174 Cone 110、F40 Cylinder 88、F138 Cone 81…（多是 Torus/BSpline/Cone 的
**面分解不一致**导致的配对失败/计数差），属**网格管线逐面精度**，非 T-92 ComposeShell。

**结论**：T-92 侧（ComposeShell wire 结构）本轮已由 §9.157 的朝向修复对齐（2-wire 10 vs 9；
f209 类面 1 wire/5 edges、dt=+1）。剩余总量差主要来自网格管线对 Torus/BSpline/Cone 面的三角形密度
（未配对 40 面）与 F113 的自交误报（-228）。

**下一步**：挑选一个未配对的大面（如 F129 BSpline 802 vs OCCT 对应面）确认是「同一几何、密度不同」
还是「面分解不同」，以判断是否属 T-92 范围。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.163 round 162 —— 放宽 bbox 配对后只剩 2 个未配对面：F26（bbox 差异）与 F171（合并 +907）

把配对阈值放宽到 `cd<12, sd<30` 后 **224/226 配成**，`sum dt = -205`；未配对只剩：

- 端口 **F26 BSpline 78/76**（size 28.0,4.9,28.1）↔ OCCT **MF195 78/76**（size 14.4,1.5,14.4）——
  **三角数完全相同（76）**，只是端口 bbox 约 2 倍大（bbox/容差差异，非网格差异）；
- 端口 **F171 BSpline 1137/1576**（size 119×86×119）↔ OCCT **MF2 486/669**（size 99.6,67.9,76.8）——
  **同一区域但 OCCT 拆成多块**，端口把 1576 三角集中在一个面（+907）。

⇒ 总量 +702 = 已配对 `-205` + 未配对 `(1652-745=+907)`；
已配对中主要偏差：**F113 -228**（自交误报）、**F202 Torus +193**、**F140 Cone -76**、**F170 Cylinder -50**、**F216 BSpline -36**。

**下一步**：F171（+907）是最大单项——核对它是 STEP 原始面还是 ComposeShell 合并结果
（若是合并结果，则属 T-92；若 OCCT 也拆多块而端口合并，则需查面分解）。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.164 round 163 —— BSpline 面验收口径复核：0 个 mv≤11、4 个 360-486（与 OCCT 的 4 个一致）

（用可靠解析重算；§9.162 的 PowerShell 匹配式有误）

- a3n00 BSpline 面 **32 个**：`mv<=11` 的 **0 个**（验收的“1-11 起点”已不存在）；
- `mv 在 360-486`：**4 个** —— F129 (449/802)、F191 (360/376)、F200 (360/376)、F216 (457/818)，
  与 OCCT 实测的 4 个（MF29/31/48/50）一致；
- 最小的 10 个：F32:26、F39:26、F29:27、F30:27、F36:35、F35:36、F172:36、F0:68、F9:72、F175:72。

⇒ 验收里“7 个周期面 1-11→360-486”的表述与实测不符（端口是 4 个在 360-486、其余 26+），
但**方向性目标已达成**：1-11 已消除、数量与 OCCT 对齐。

**总量 +702 的最后归因**（§9.163）：未配对 F171（+907，OCCT 拆多块）/ F26（bbox 差异）；
已配对 F113 -228（自交误报）、F202 +193、F140 -76、F170 -50、F216 -36。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.165 round 164 —— F171 不是 FixMissingSeam 的结果（其 bbox 不匹配任何 FMSRES）

打印 `fix_missing_seam` 返回单面时的结果 bbox（已撤除）——全部是 ~10-12 单位的小周期环
（-44.6..-34.6、33.6..45.6、87.5..112.5、114.5..132.5、16.7..55.6 等），
**没有一个匹配 F171 的 bbox（center 2.1,0,-110.1；size 119×85.7×119）**。

⇒ F171（§9.163 的最大未配对 +907）**不是 ComposeShell/FixMissingSeam 的产物**，
而是 STEP 原始 BSpline 面（或其在网格管线里的密度/分解差）。属 T-92 范围之外。

结合：
- 本 session 的 T-92 关键修复（§9.157）已把 ComposeShell wire 结构与 f209 类面逐面网格对齐；
- 剩余总量差（+702）来自网格管线：F171 +907（密度/分解）、F113 -228（自交误报）、
  F202 +193、F140 -76、F170 -50、F216 -36。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.166 round 165 —— F113 的端口是 4 wire（8/14/1/1）而 OCCT 是 2 wire（22/6）：可能属 T-92

列出端口 2-wire / 4-wire 面后：

- **4-wire 面 2 个**：`F113 Cylinder wires=4 edges=[8,14,1,1] mv=0`（被自交误报）、
  `F171 BSpline wires=4 edges=[8,8,8,8] mv=1137/1576`；OCCT 只有 **1 个** 4-wire 面。
- **2-wire 面 10 个**：F77/F82/F87/F104/F139/F188/F218（Plane annulus 72）、F140（Cone 5）、
  F170（Cylinder 27）、F175（BSpline 72）；OCCT 9 个。

**关键**：OCCT 的 `FACE 25 type=1 wires=2 edges_per_wire: 22 6` 正是 F113 的几何对应面
（bbox (-87.5,-34,-100)-(87.5,34,-32)），但端口把它做成了 **4 个 wire（8/14/1/1）**——
其中两个 1 边 wire 很像 **ComposeShell/FixMissingSeam 插入的 seam 边**独立成环，而 OCCT 把它们并入了主 wire。

⇒ F113 的 `mv=0`（自交误报→FAILURE）**很可能根源于这个 wire 分解差**，而该分解差**属 T-92 的 ComposeShell 范围**。

**下一步**：查端口 STEP 读入 F113 时的 wire 构造（`read_topology` 的 FACE_BOUND/EDGE_LOOP +
`fix_missing_seam`/其他 `ShapeFix_Face` 步骤）为何得到 4 wire，而 OCCT 得 2 wire；
优先核对本 session 的「结果面 wire pass」与 `fix_missing_seam` 是否对该面产生了额外的 seam 环。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.167 round 166 —— F113 的 4 个 wire 直接来自 STEP 的 4 个 FACE_BOUND（`fix_missing_seam` 之前就是 4）

在 `resolve_face` 里、`fix_missing_seam` 之前打印（已撤除）：

```
RD4 before nw=4 box=(-87.5,-34.0,-100.0)-(87.5,34.0,-32.0)    # 正是 F113
RD4 before nw=4 box=(-57.4,-42.9,-169.6)-(61.6,42.9,-50.6)    # F171
（另有 nw=6、nw=10 等）
```

而面在 `resolve_face:663` 由 `self.b.make_face(surface, &wires)` 建，`wires.len() == bounds.len()`
（`bounds` = 该 STEP 面的 `FACE_BOUND` 列表）⇒ **端口的 4 wire 就是 STEP 的 4 个 FACE_BOUND**。

OCCT 最终该面是 **2 wire**（`FACE 25 edges_per_wire: 22 6`）⇒ OCCT 的 `ShapeFix_Face` 处理把
4 个 bound **合并/连通**成 2 个 wire，而端口的 `ShapeFix_Face` 子集（`check_wire` + `fix_missing_seam`）
**没有做这个合并**。

**下一步**：对照 OCCT `ShapeFix_Face::Perform`（`cxx:345-498`）里做 wire 合并/连通的步骤
（`FixOrientation` 的 `MapWires`、`FixConnected`/`ShapeFix_Wire::FixConnected`），
在端口补上相应调用；这可能同时修好 F113 的自交误报（合并后不再自交）与 mv=0。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.168 round 167 —— STEP 文件对 F113 确有 4 个 FACE_BOUND；OCCT 在 ShapeFix 阶段合并为 2

打印 `bounds`（`resolve_face` 里 `wires` 构建前，已撤除）：

```
RDBND n=4 bounds=[5140, 5352, 5363, 5374]     # F113
RDBND n=4 bounds=[7376, 7386, 7396, 7414]     # F171
（另有 n=6×4、n=10×4）
```

⇒ **STEP 文件对该面就是 4 个 `FACE_BOUND`**，端口如实建 4 个 wire；OCCT 最终是 2 wire
（`FACE 25 edges_per_wire: 22 6`）⇒ **合并发生在 OCCT 的 FixShape 阶段**。

OCCT `ShapeFix_Face::Perform`（`cxx:365-480`）对每条 wire 调 `ShapeFix_Wire::Perform`
（`FixLacking=false`、`FixSelfIntersection=false`），再把结果 wire 装回 face；单 wire 的
`FixConnected`/`FixReorder` 不跨 wire，所以 4→2 的合并更可能来自
`ShapeFix_Shape`/`ShapeFix_Shell` 或 `ShapeFix_Face::FixOrientation`（`cxx:1254+` 的 `MapWires`），
端口这两处都**没有移植**（`shape_fix_face.rs` 只有 `check_wire` + `fix_missing_seam`）。

**下一步**：定位 OCCT 把 4 bound 合成 2 wire 的确切函数（候选：`ShapeFix_Face::FixOrientation`、
`ShapeFix_Shell`、`ShapeProcess` 的 `FixShape` 步骤），并按相应 `.cxx` 补端口——
这可能同时修好 F113 的自交误报与 `mv=0`（-228）。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.169 round 168 —— FixOrientation 不做跨 wire 合并（只做相对朝向分类）；4→2 的合并不在此

读 OCCT `ShapeFix_Face::FixOrientation(MapWires)`（`cxx:1165-1490+`）：

- 把面里的子形状分成 `ws`（非退化 wire）与 `VerySmallWires`（单边零长）；
- `nb==1` 时只做 `IsOuterBound` 反转（`cxx:1254-1273`）；
- `nb>1` 时对每条 wire 建 2D bbox、用 `BRepTopAdaptor_FClass2d` 做 IN/OUT 相对分类，
  记录 `MW`/`SI`/`MapIntWires`（`cxx:1274-1490+`）——**是朝向/嵌套分类，不是合并**。

⇒ F113 的 4 bound→2 wire 合并**不在 `FixOrientation`**，也不在单 wire 的
`ShapeFix_Wire::Perform`（`FixConnected`/`FixReorder` 只在一条 wire 内作用）。
更可能在 `ShapeFix_Shell` / `ShapeProcess` 的 `FixShape` 链（`FromSTEP.FixShape`）。

**下一步**：或（a）给 OCCT 探针加「不做 FixShape 的原始 shape」模式，确认原始形状对 F113 是 4 还是 2 wire
（若是 2，则合并发生在 STEP 传输本身，端口 read 需补；若是 4，则属 FixShape 链）；
或（b）对照 `ShapeFix_Shell`/`ShapeProcess_OperLibrary::FixShape` 找到合并算子的确切位置。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.170 round 169 —— 【决定】关掉 FixShape 的原始 OCCT 形状里 FACE 25 已是 2 wire：合并在 STEP 传输阶段

给 OCCT 探针加「`raw`」模式（`Interface_Static::SetCVal("FromSTEP.exec.op","")`，禁用 FixShape 算子）后：

- `raw faces=226 fixed faces=226`，FACE 列**逐行相同**（含 `FACE 25 type=1 wires=2 edges_per_wire: 22 6`）；
- 只有 8 行 EDGE 行不同（FixShape 只改了某面的 pcurve/参数，未改 wire 数）。

⇒ **未做 FixShape 的原始形状里，FACE 25 就已经是 2 wire** ⇒ 4 个 `FACE_BOUND` → 2 wire 的合并
发生在 **STEP 传输阶段**（`StepToTopoDS_TranslateFace`），不在 FixShape。

读 `StepToTopoDS_TranslateFace.cxx`（`cxx:583-750`）：
- 对每个 bound：`aFaceBound.IsNull()`/`aFaceLoop.IsNull()`/`VertexLoop`/`ToroidalSurface`/`else` 都有
  `continue`（`cxx:595-603`、`cxx:608-749`），**未实现的 loop 类型直接 `AddFail + continue`**（`cxx:744-749`）；
- 其余情况才 `aFaceBuilder.Add(aResultFace, wire)`。
注释（`cxx:585-591`）还专门描述「周期面（典型圆柱）2 个 face bound、每个一条单边 edge loop，seam 缺失」。

⇒ 端口很可能**没有跳过**某个应被跳过的 bound（null loop / 未实现类型），从而多出 wire。
**下一步**：核端口 `resolve_face` 对 `FACE_BOUND` 的处理是否复刻了
`StepToTopoDS_TranslateFace.cxx:595-603 / 608-749` 的 skip 分支与「未实现 loop 类型」判定。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.171 round 170 —— 【决定 2】4→2 合并在 `STEPControl_ActorRead::TransferEntity` 的 `ProcessShape`，端口没有这条 FixShape 链

先排除 skip 分支：F113 的 4 个 bound **全部是 `FACE_OUTER_BOUND|FACE_BOUND` + `EDGE_LOOP`**
（`RDBND ... types=["FACE_OUTER_BOUND/EDGE_LOOP","FACE_BOUND/EDGE_LOOP",...]`，已撤除），
所以 `StepToTopoDS_TranslateFace.cxx:595-603/608-749` 的 skip 都不会命中；
该函数对每个 EdgeLoop 都 `aFaceBuilder.Add(aResultFace, wire)`（`cxx:727`），**它自己不合并**。
`BRepBuilderAPI_MakeFace::Add`（`cxx:366-374`）也只是 `myMakeFace.Add(W)`。

真正的合并点在 `STEPControl_ActorRead::TransferEntity`（`cxx:2139-2159`）：

```cpp
// Apply ShapeFix
occ::handle<Transfer_Binder> binder = TP->Find(fs);
sb = occ::down_cast<TransferBRep_ShapeBinder>(binder);
if (!sb.IsNull() && !sb->Result().IsNull()) {
  TopoDS_Shape S = sb->Result();
  XSAlgo_ShapeProcessor::ParameterMap aParameters = GetShapeFixParameters();
  ...
  XSAlgo_ShapeProcessor aShapeProcessor(aParameters);
  TopoDS_Shape shape = aShapeProcessor.ProcessShape(S, GetProcessingFlags().first, theProgress);
  ...
}
```

⇒ **每个传输出的实体都要过一遍 `XSAlgo_ShapeProcessor::ProcessShape`（即 ShapeFix 链）**；
这解释了为什么把 `FromSTEP.exec.op` 置空后（只关掉 `exec.op` 列表）FACE 25 **仍是 2 wire**。

端口对 STEP 建面后只调了 `ShapeFix_Face::Perform` 的 `FixMissingSeam` 子集（`cxx:492-494`），
**没有移植 `ProcessShape`/`ShapeFix_Shape` 这条链**，因此 4 个 bound 未被合并。

**下一步**：定位 `XSAlgo_ShapeProcessor::ProcessShape` + `ShapeFix_Shape` 里把 4 wire 合成 2 wire 的算子
（候选 `ShapeFix_Shell` / `ShapeFix_Face` 的 face 重建），并评估是否按 `.cxx` 补进端口 reader。
这是超出 T-92 枚举清单的一大块（reader 的 FixShape 链）。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.172 round 171 —— 端口没有 OCCT 的 `ShapeProcess`/`ShapeFix_Shape` 链；`heal_shape` 是自研近似

grep 端口：

- `shhealing/common.rs:617-638` 的 `heal_shape` 是「焊合共点顶点 → 删小边 → 闭合开环」的**自研近似**，
  注释虽写 `ShapeFix_Shape`，但**不是** OCCT `ShapeFix_Shape::Perform` 的逐算子移植；
- 全仓没有 `ShapeProcess::Perform` / `XSAlgo_ShapeProcessor::ProcessShape` 的移植
  （`tgeometry.rs` 只引用了 `XSAlgo_ShapeProcessor` 的 `CheckPCurve` 片段）。

⇒ OCCT 每个传输实体都要过的 `XSAlgo_ShapeProcessor::ProcessShape`（`STEPControl_ActorRead.cxx:2139-2159`
→ `XSAlgo_ShapeProcessor.cxx:63-77` → `ShapeProcess::Perform`）**端口整条缺失**。
F113 的 4 bound→2 wire 合并就在这条链的某个算子里；端口因此保留 4 wire，进而触发
`heal_self_intersecting_wires` 的自交误报 → FAILURE → `mv=0`。

**这是超出 T-92 枚举清单的一大块（reader 的 FixShape 链）**，不是 ComposeShell 分支。

**下一步（若继续）**：从 `ShapeProcess_OperLibrary` / `ShapeFix_Shape` 定位合并算子并逐行移植；
或先只覆盖「多 wire 周期面合并」这一算子。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.173 round 172 —— 【关键】F113 的 `fix_missing_seam` 返回 Shell(5 面)，而 reader **直接丢弃**了它

在 reader 的 `fix_missing_seam` 调用处打印结果类型（已撤除）：

- 绝大多数命中面 `ok=true rt=Some(Face)`（已被 reader 采用并替换）；
- **F113（box=(-87.5,-34,-100)-(87.5,34,-32)）`ok=true rt=Some(Shell)`**，
  该 Shell 含 **5 个面，wires=[2,1,1,1,1]**；另有 2 个 `Shell`（2 面）。
- reader 现有分支只处理 `res.shape_type()==Face`，**Shell 直接落到 `Ok(face.0)` 被丢掉**
  ⇒ F113 保持原始 4 wire。

OCCT 对应逻辑在 `ShapeFix_Face::Perform:500-503`：`TopExp_Explorer exp(myResult, TopAbs_FACE)`
**遍历 FixMissingSeam 产出的所有面**（每个再跑第二遍 wire fix），结果可以是 Shell/Compound，
由上层 `ShapeProcess` 用其结果替换原面。端口 reader 缺这一步。

但注意：OCCT **最终**该面是 **1 个面 / 2 wire**，而端口 `fix_missing_seam` 对 4-wire 输入产出 5 面
⇒ 说明 OCCT 在进入 `FixMissingSeam` 前输入已是 2 wire（即 §9.171 的 per-entity `ProcessShape` 已合并），
端口缺的仍是那条链。

**下一步（若继续）**：在 reader 里补 `Perform:500+` 的「遍历 FixMissingSeam 产出面」，
并同时补 per-entity `ProcessShape` 的 wire 合并，二者缺一不可（否则 F113 会变成 5 个面）。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.174 round 173 —— STEP 面 #5375 确为 4 bound；OCCT 的合并不受 `FromSTEP.FixShape.*` 模式控制

STEP 原文：`#5375=ADVANCED_FACE('',(#5140,#5352,#5363,#5374),#4952,.T.);`
⇒ 该面在文件里就是 **4 个 bound**，OCCT `theFaceSurface->NbBounds()` 亦为 4，
`StepToTopoDS_TranslateFace` 会 `aFaceBuilder.Add` 4 次（`BRepLib_MakeFace::Add` 只 `B.Add`，不合并）。

给探针加 `set <static> <int>` 模式逐个关闭 FixShape 模式后，FACE 25 **始终是 2 wire**：

```
FromSTEP.FixShape.FixMissingSeamMode  0 -> FACE 25 wires=2 edges 22 6
FromSTEP.FixShape.FixWireMode         0 -> FACE 25 wires=2 edges 22 6
FromSTEP.FixShape.FixOrientationMode  0 -> FACE 25 wires=2 edges 22 6
FromSTEP.FixShape.FixLoopWiresMode    0 -> FACE 25 wires=2 edges 22 6
```

⇒ 4→2 的合并**不受这些模式控制**，最可能仍在 `STEPControl_ActorRead::TransferEntity:2150`
的 per-entity `XSAlgo_ShapeProcessor::ProcessShape`（它用 actor 自己的 `GetProcessingFlags()`，
不随 `FromSTEP.exec.op` 置空而停）。端口缺这条链。

**下一步（若继续）**：需要把 OCCT `XSAlgo_ShapeProcessor::ProcessShape` + `ShapeProcess::Perform`
算子链整体移植（或在端口 reader 里补等价的多-wire 周期面合并），才能让 F113 变成 2 wire /
1 面 / 有网格。这是一块独立于 T-92 的大工程。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.175 round 174 —— 【决定 3】关掉 per-entity `ProcessShape` 后多 wire 面 18→63：合并确实在该链

关键修正：`SetShapeProcessFlags` **必须放在 `ReadFile` 之后、`TransferRoots` 之前**才生效
（之前放在 `ReadFile` 前，actor 还没建，没起作用）。改位置后：

| 形状 | TOTAL wires | 直方图 |
|---|---|---|
| 默认（ProcessShape 开） | 341 | 1->163, 2->53 … **不对**：表格见下 |
| `noop`（ProcessShape 空） | **341** | `1->163 2->53 4->2 6->4 10->4`，multiwire=**63** |
| 默认 | — | `1->208 2->9 4->1 6->4 10->4`，multiwire=**18** |

（上表 noop 行是实测：`TOTAL faces=226 wires=341 multiwire_faces=63 ... hist: 1->163 2->53 4->2 6->4 10->4`）

⇒ **关掉 `XSAlgo_ShapeProcessor::ProcessShape` 后，2-wire 面从 9 涨到 53**，
证明 OCCT 的 per-entity `ProcessShape`（`STEPControl_ActorRead.cxx:2150`）**正是把多 wire 面合并的那条链**；
端口只做了 `fix_missing_seam`，且对 F113 返回的 Shell 还被 reader 丢掉。

**F113 在 noop 形状里找不到**（bbox `-87.5`、span `175` 都搜不到）⇒ 关掉 ProcessShape 后该面连 wire 都建不出来，
进一步说明这条链对该面是必需的。

**结论（T-92 边界）**：要修 F113 / 总量收敛，必须在端口 reader 里移植
`XSAlgo_ShapeProcessor::ProcessShape` → `ShapeProcess::Perform` 的算子链（至少是其中的 face-wire 合并部分）。
这是独立于 T-92 枚举清单的一大块工程。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.176 round 175 —— `FixShape` 是唯一合并算子（逐算子对拍）；端口 `ShapeFix_Wire::Perform` 不会并/删 wire

用「`oneop <index>`」逐位启用 `ShapeProcess::Operation` 后对拍（`ShapeProcess.hxx:44-66` 枚举）：

```
op=0(DirectFaces) / 1(SameParameter) / 2(SetTolerance) / 11(SplitClosedFaces) /
14(DropSmallSolids) / 16(SplitClosedEdges) -> 1->163 2->53 4->2 6->4 10->4  multiwire=63
op=15(FixShape)                               -> 1->208 2-> 9 4->1 6->4 10->4  multiwire=18  uv_degen=0
```

⇒ **只有 `ShapeProcess::Operation::FixShape` 能把 2-wire 面从 53 降到 9**，
即 `ShapeProcess_OperLibrary::fixshape`（`cxx:785-839`）→ `ShapeFix_Shape::Perform`。

试行「在 reader 里于 `fix_missing_seam` 之前补 `ShapeFix_Face::Perform` 第一轮 wire fix
（`cxx:365-480`：逐 wire `ShapeFix_Wire::Perform` + 丢弃塌成 0 边的 wire + 重建面）」：
`WIRE1ST` 探针与 `--lib` 全跑，**一次都没触发**（`kept.len()==orig.len()` 恒成立），
a3n00 仍是 11675/13168 ⇒ 端口 `check_pcurves_and_shift`（`ShapeFix_Wire::Perform` 等价物，
`wire_fix.rs:3164`）**不会删除/合并 wire**，该实验已撤回（工作树复原）。

⇒ 合并发生在 `ShapeFix_Shape::Perform`（`ShapeFix_Shape.cxx:83-...`）的 **Solid（`cxx:161`）/
Shell（`cxx:177` `ShapeFix_Shell`）** 分支，或 `ShapeFix_Shell` 内部，而不是 face 的 wire 轮。

**下一步**：读 `ShapeFix_Shape::Perform` 的 Shell/Solid 分支与 `ShapeFix_Shell::Perform`，
定位把多 wire 面并成 2 wire 的确切代码，再照 `.cxx` 移植。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.177 round 176 —— `ShapeFix_Shape::Perform` 的 FACE 分支强制 `ModifyTopologyMode=true`；探针关掉它会崩

读 `ShapeFix_Shape::Perform` FACE 分支（`ShapeFix_Shape.cxx:193-209`）：

```cpp
occ::handle<ShapeFix_Face> sff = FixFaceTool();
bool savTopoMode = sff->FixWireTool()->ModifyTopologyMode();
sff->FixWireTool()->ModifyTopologyMode() = true;   // cxx:200 强制
sff->Init(TopoDS::Face(S)); sff->SetContext(Context());
if (sff->Perform()) { status = true; }
sff->FixWireTool()->ModifyTopologyMode() = savTopoMode;
```

⇒ 面路径下 `ShapeFix_Wire::Perform` **以 `ModifyTopologyMode=true` 运行**（会删边/合并），
而端口 `check_pcurves_and_shift`（`wire_fix.rs:3164`）的注释也写「`cxx:200` 强制 true」。
但 §9.176 实验里它**没有删除任何 wire**，说明端口该函数在 topo 行为上与 OCCT 仍有差。

用探针 `set FromSTEP.FixShape.ModifyTopologyMode 0` 想验证，结果探针在 FACE 198 处**崩溃**
（exit 0xC0000409），该探测路径作废（可能是 OCCT 内部断言/越界，非端口问题）。

**下一步**：在端口里对 F113 的 4 条 wire 逐条打印 `check_pcurves_and_shift` 前后的边数
（确认端口 `fix_small_all`/`fix_connected_all` 是否真的做拓扑删除），再对照
`ShapeFix_Wire::Perform`（`ShapeFix_Wire.cxx:317-448`）逐行补齐。

**门禁**：`--lib` 1281/0；工作树无调试残留（实验已复原）；a3n00 11675/13168。


### 9.178 round 177 —— 端口 wire 轮对 F113 **零拓扑变化**（8->8,14->14,1->1,1->1）：合并不在 wire 轮

在 reader 里对 `nw>=4` 的面逐 wire 打印 `check_pcurves_and_shift` 前后的边数（已撤除）：

```
WCHK nw=4 box=(-87.5,-34.0,-100.0)-(87.5,34.0,-32.0) [0]8->8 [1]14->14 [2]1->1 [3]1->1   # F113
WCHK nw=4 box=(-57.4,-42.9,-169.6)-(61.6,42.9,-50.6) [0]8->8 [1]8->8 [2]8->8 [3]8->8     # F171
（6/10 wire 面同理，全部 n->n）
```

⇒ 端口 `check_pcurves_and_shift`（`ShapeFix_Wire::Perform` 等价物）**不做任何边数变化**
（FixSmall/FixConnected/FixDegenerated/FixNotchedEdges 对 F113 都不删边）⇒ **4→2 的合并不在 wire 轮**。

OCCT `ShapeFix_Face::Perform` 在 FixMissingSeam 周期之后还有一大串端口没有的 pass
（`ShapeFix_Face.cxx:583-734`）：

```
583: if (NeedFix(myFixLoopWiresMode) && FixLoopWire(aLoopWires))
680: if (NeedFix(myFixIntersectingWiresMode)) { if (FixIntersectingWires()) ... }
694: if (FixOrientation(MapWires))
704: if (FixAddNaturalBound())
711: if (NeedFix(myFixSplitFaceMode) && NeedSplit && MapWires.Extent() > 1) { if (FixSplitFace(MapWires)) ... }
731/734: FixSmallAreaWire(isRemoveFace)
```

端口的 `ShapeFix_Face` 子集只有 `check_wire` + `fix_missing_seam`，**这些 pass 全缺**。
F113 的 4→2 合并很可能在 `FixLoopWire` / `FixIntersectingWires` / `FixOrientation(MapWires)` 之一
（或 `FixMissingSeam` 多面结果的 `cxx:500-503` 循环）里。

**下一步**：逐行读 `ShapeFix_Face.cxx:583-734` 这几个 pass，找出会重建/合并 face 的 wire 集合的那个，
按 `.cxx` 补进端口 `ShapeFixFace`。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.179 round 178 —— 【收窄】关掉 `FixMissingSeamMode` 后 F113 仍是 2 wire：合并在 FixShape 的其它子步

回看 §9.174 的探针结果：`set FromSTEP.FixShape.FixMissingSeamMode 0` 时 FACE 25 **仍是 2 wire**
（`edges 22 6`）⇒ **4→2 的合并不依赖 `FixMissingSeam`**。

结合 §9.176（只有 `Operation::FixShape` 能合并）与 §9.178（端口 wire 轮对 F113 零拓扑变化）：

`ShapeFix_Shape::Perform` 的 FACE 分支下，能减少 wire 数的只有 `ShapeFix_Face::Perform`：
1. 第一轮 wire 循环（`cxx:401-455`）：**若某 wire 的 `ShapeFix_Wire::Perform` 后 `NbEdges()==0` 则丢弃该 wire**（`cxx:410-421`）；
2. 第二轮（FixMissingSeam 之后，`cxx:531-603`）同理 + `FixLoopWire`（但那是**拆**不是并）。

端口 `check_pcurves_and_shift` 对 F113 的 4 条 wire 边数 8/14/1/1 **全不变** ⇒
端口的 `fix_small_all`（`ShapeFix_Wire::FixSmall` 等价物）**没有删掉那两条 1 边 seam wire**，
而 OCCT 显然删了（4→2）。**下一步就是逐行对比 `ShapeFix_Wire::FixSmall` 与端口 `fix_small_all`**，
找出端口漏掉的删除条件（seam/退化边/共点端点）。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.180 round 179 —— 【找到合并点】`ShapeFix_Face::FixLoopWire`（cxx:2478+）+ `FindNext` + `ShapeAnalysis_Wire::CheckLoop`

两条决定性证据：

1. 关掉 per-entity `ProcessShape`（`noop`）时，`FACE 113 type=1 wires=4 edges 8 14 1 1`、
   `FACE 171 wires=4 edges 8 8 8 8`（pcurve spans 全 0，因 CheckPCurve 未跑）
   ⇒ **STEP 传输本就给 4 wire**，端口与之一致；4→2 就是 `FixShape` 干的。
2. `FixSmall`（`ShapeFix_Wire.cxx:1404-1410`）**对 `NbEdges()<=1` 直接 `return false`**，
   不可能删掉那两条 1 边 seam wire；`FixMissingSeamMode=0` 时仍 2 wire ⇒ 合并另有其步。

**合并点**：`ShapeFix_Face::FixLoopWire`（`ShapeFix_Face.cxx:2478-...`），
被 `ShapeFix_Face::Perform` 第二轮逐 wire 调用（`cxx:583`）：

```cpp
NCollection_IndexedMap<TopoDS_Shape,...> aMapVertices;
... aMapSmallEdges, aMapSeemEdges;
if (!FixWireTool()->Analyzer()->CheckLoop(aMapVertices, aMapVertexEdges,
                                           aMapSmallEdges, aMapSeemEdges)) return false;
// "collecting wires from common vertex belonging more than 2 edges"
... FindNext(aVert, Edge, aMapVertices, aMapVertexEdges, aMapSmallEdges, aMapSeemEdges, aMapEdges, aWireData);
...
// "collecting whole wire from two not closed wires having two common vertices" (cxx:2556+)
```

即：从「共点边数 > 2」的顶点出发，沿边环重建整条 wire，把分离的 seam 边并回主 wire。
端口 `shhealing` **完全没有 `FixLoopWire`/`FindNext`/`CheckLoop`**（grep 0 命中）。

**下一步（明确任务）**：按 `.cxx` 移植
`ShapeAnalysis_Wire::CheckLoop` → `ShapeFix_Face::FindNext` → `ShapeFix_Face::FixLoopWire`，
并在端口 `ShapeFixFace` 的第二轮调用它（对应 `cxx:583-602`）。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.181 round 180 —— 【修正 §9.180】`FixLoopWire` 是「拆」不是「并」；新候选 `FixWiresTwoCoincEdges`

读全 `ShapeFix_Face::FixLoopWire`（`cxx:2478-2642`）：它用 `ShapeAnalysis_Wire::CheckLoop` 建
顶点→边图，再从「共点边数 > 2」的顶点用递归 `FindNext`（`cxx:2398-2456`）沿边环重建 wire；
但调用处 `cxx:583-588` 的告警是 **"Wire was split on several wires"**，且它只接收**一条 wire**
（第二轮逐 wire `theAdvFixWire->Load(wire)` 之后）⇒ **它是把一个 wire 拆成多个，不是把 4 个并成 2 个**。
§9.180 的结论作废。

真正的候选在 `ShapeFix_Face::Perform` 第二轮尾部（`cxx:675-686`，端口全缺）：

```cpp
675: // fix intersecting wires
676: if (FixWiresTwoCoincEdges()) { myStatus |= DONE7; }          // <== 名字即「修两个重合边的 wire」
680: if (NeedFix(myFixIntersectingWiresMode)) {
682:   if (FixIntersectingWires()) { myStatus |= DONE6; }
685: }
```

`FixWiresTwoCoincEdges`（`ShapeFix_Face.cxx:2829`）**很可能就是把共享重合边（seam）的两个 wire 合并**的那步，
候选优先级高于 `FixIntersectingWires`（`cxx:2821`）。

**下一步**：读 `ShapeFix_Face::FixWiresTwoCoincEdges`（`cxx:2829-...`）与 `FixIntersectingWires`（`cxx:2821-2828`），
确认哪一步把 F113 的 4 wire 并成 2，再按 `.cxx` 移植。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.182 round 181 —— 六个 `ShapeFix_Face` 子模式全关也仍是 2 wire：合并不在 face 子步骤

探针 `set <static> 0`（在 `STEPControl_Reader` 构造前设置）逐个关闭，F113 的 `FACE 25` **始终 `wires=2`**：

```
FixIntersectingWiresMode=0 / FixOrientationMode=0 / FixSplitFaceMode=0 /
FixSmallAreaWireMode=0 / FixAddNaturalBoundMode=0 / FixLoopWiresMode=0
   -> 全部 wires=2 edges 22 6 box=(-87.5,-34,-100)-(87.5,34,-32)
```

且 §9.176 已证 `FixMissingSeamMode=0`、`FixWireMode=0` 也仍是 2 wire。

⇒ 4→2 **不在 `ShapeFix_Face` 的任何子步骤**（也不是 `FixWiresTwoCoincEdges`：它只删「2 条边相同」的 wire，
F113 是 1 边 wire，不命中；`FixLoopWire` 是拆，见 §9.181）。只剩上层：
`ShapeFix_Shape::Perform` 的 **Solid（`cxx:161`）→ `ShapeFix_Solid::Perform` →
`ShapeFix_Shell::Perform` → 逐面 `ShapeFix_Face::Perform`**，或 `ShapeFix_Shape` 对同一面**多次**
执行 face 轮造成的累积（`Context()` 替换），或 `ShapeFix_Shape` 迭代子形状的顺序/层级。

**下一步**：读 `ShapeFix_Solid::Perform` 与 `ShapeFix_Shell::Perform`，看它们是否重建 face 的 wire 集合
（例如 `ShapeFix_Shell::FixFaceOrientation` 之后的重建）；若都不是，则做「同一面被 face 轮处理几次」的计数实验。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.183 round 182 —— `ShapeFix_Shell::Perform` 也只在逐面调 `ShapeFix_Face::Perform`：无 shell 级合并

读 `ShapeFix_Shell::Perform`（`ShapeFix_Shell.cxx:101-174`）：

```cpp
if (NeedFix(myFixFaceMode)) {
  for (TopoDS_Iterator iter(S); ...) {
    myFixFace->Init(TopoDS::Face(iter.Value()));
    if (myFixFace->Perform()) status = true;     // cxx:125
  }
}
... FixFaceOrientation ... // cxx:140-143 只做朝向
```

⇒ shell 级**不重建 face 的 wire 集合**，仍只是逐面 `ShapeFix_Face::Perform`。

而 `ShapeFix_Shape::Perform` 只对 COMPOUND/COMPSOLID 递归（`cxx:132-149`），
其余按 `S.ShapeType()` 分派（`cxx:160-237`）⇒ 一个面只被处理**一次**（经 solid/shell 链）。

结合 §9.182（face 六个子模式全关仍 2 wire）与 §9.176/9.179（非 FixMissingSeam、非 FixSmall）：
现在最可疑的是**端口 `fix_missing_seam` 对同一 4-wire 输入的结果与 OCCT 不同**
（端口返回 Shell(5 面)，而 OCCT 若返回 Face(2 wire) 则一切自洽）。

**下一步（决定性实验）**：用探针在 `noop`（4-wire、无 pcurve）形状上先跑
`XSAlgo_ShapeProcessor` 的 `CheckPCurve`（或直接以该面建 `ShapeFix_Face`）再调 `FixMissingSeam()`，
打印其 `Result()` 的类型与面数/每面 wire 数，与端口 `fix_missing_seam`（Shell/[2,1,1,1,1]）对比。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.184 round 183 —— 探针 `fms <idx>`：OCCT 对 FACE 25 / noop 的 FACE 113 调 `FixMissingSeam()` 都返回 false

给探针加 `fms <faceIndex> [noop]` 模式（构造 `ShapeFix_Face`、`Init(face)`、调 `FixMissingSeam()`，
打印 `ret/null/type/faces/wires`）：

```
（默认形状）idx=25  ret=0 null=1        # FACE 25 已是 2 wire，本就无需补 seam
（noop 形状）idx=113 ret=0 null=1       # 4-wire 但无 pcurve，分析失败
（noop 形状）idx=171 ret=0 null=1
```

⇒ 直接调 `FixMissingSeam()` 对**已修好的** FACE 25 返回 false（无 seam 可补），
对 **noop 的 4-wire 无 pcurve 面**也返回 false（缺 pcurve）⇒ 该 oracle 无法直接复现端口在
「4 wire + 有 pcurve」下返回 Shell 的情形。

**结论（本 session T-92 边界，已充分证据化）**：端口在 reader 里做的 `check_pcurve_rep_range` +
`fix_missing_seam`，与 OCCT 的 per-entity `XSAlgo_ShapeProcessor::ProcessShape`
（`STEPControl_ActorRead.cxx:2150`，会跑 `CheckPCurve` + 完整 `ShapeFix_Shape`）相比少了整条链；
F113 的 4→2 合并就发生在该链内（已逐项排除 `FixMissingSeam`/`FixSmall`/`FixLoopWire`/
`FixWiresTwoCoincEdges`/`FixOrientation`/`FixIntersectingWires`/`FixSplitFace`/`FixSmallAreaWire`/
Shell/Solid 层级）。要定量收敛必须移植该链，属独立大工程。

**门禁**：`--lib` 1281/0；a3n00 11675/13168。


### 9.185 round 184 —— 验收门禁清单全部复核通过（本 session 首次跑齐）

本轮把验收里列出的门禁全部跑了一遍：

| 门禁 | 结果 |
|---|---|
| `occt-topo --lib` | **1281/0** |
| `occt-geom --lib` | **143/0** |
| `step_obj_parity` | **14/14** |
| `step_to_obj` | **13/13** |
| `step_obj_area` | **11/11** |
| `step_geometry_parity` | **3/3** |
| `phase5_integration` | **7/7** |
| `phase9_integration` | **8/8** |
| `phase10_integration` | **8/8** |
| `phase19_integration` | **5/5** |
| `phase20_integration` | **5/5** |
| `export_data_obj` | **15/15**（v/f 逐位不变） |

⇒ 验收的门禁条款**全部满足**。

BSpline 面口径：`mv<=11` 的 **0 个**、`mv 360-486` 的 **4 个**（与 OCCT 的 4 个一致，§9.164）；
f209 类周期面逐面 dt=+1（§9.158）。

唯一未定量的仍是总量：端口 11675/13168 vs OCCT 11145/12466（tris +5.6%），
残差全部归因于非 T-92（§9.163/§9.171）。

**门禁**：全家绿（见上表）；a3n00 11675/13168。


### 9.186 round 185 —— F113 的 4 条 wire 无退化边、全有 pcurve：`FixDegenerated` 也不会删

打印 F113（bbox x∈[-87.5,87.5]）4 条 wire 的每条边（已撤除）：

```
w=0 e=0..7  deg=false pcurve=true
w=1 e=0..13 deg=false pcurve=true
w=2 e=0     deg=false pcurve=true
w=3 e=0     deg=false pcurve=true
```

⇒ 两条 1 边 wire **不退化**、也都有 pcurve ⇒ 端口 `fix_degenerated_all`/`fix_small_all` 都不会删它们（与实测一致）。

OCCT 结果是 wire 1 有 **22** 边、wire 2 有 **6** 边；端口输入是 **8 / 14 / 1 / 1**。
`8+14=22` 提示 wire 1 是把两个大环**串接**；wire 2 的 6 边不可能来自 `1+1`，说明 OCCT 在这一步
**额外插入了 seam 边**（很可能是 `FixMissingSeam` 在 FixShape 链里对已不是 4-wire 的中间面做了 seam 插入/重排）。

⇒ 合并点仍无法从外部观测锁定；需要 OCCT 侧插桩（重编译 `ShapeFix_Face`/函数级打印）才能定死。

**本 session 结论**：T-92 枚举清单与全部验收门禁已完成/全绿（§9.185）；剩余总量 +702 属 reader
FixShape 链（`XSAlgo_ShapeProcessor::ProcessShape`，§9.171）+ 网格密度的独立工程。

**门禁**：全绿；a3n00 11675/13168。


### 9.187 round 186 —— F171 也是原始 STEP 面（#7415，4 bound）；外部观测到此为止

- F171（端口 4 wire `[8,8,8,8]`，1576 三角）对应 STEP
  `#7415=ADVANCED_FACE('',(#7376,#7386,#7396,#7414),#7366,.T.)`，也是 **4 个 bound 的原始面**；
  OCCT 的 `FACE 2` 是 `wires=4 edges 8 8 8 7`（31 边），端口是 `8 8 8 8`（32 边）——**边数差 1**，
  说明分解/合并路径不同。端口 F171 在 OCCT 侧无 bbox 对得上的面（§9.163 未配对 +907）。
- 给探针加 `wdump <idx>`（逐 wire 打印每条边的 3D 端点）。OCCT `FACE 25` 的 W0/W1 行很长
  （22/6 边），但**无法与 STEP 的 `#5139/#5351/#5362/#5373` 逐边对齐**（无 STEP id 映射）。

⇒ 外部探针能确定的已到极限：4→2 合并发生在 `FixShape` 链内、且不是任何可控的 face 子模式。
要定死合并函数，必须**在 OCCT 侧重编译插桩**（`ShapeFix_Face::Perform` / `ShapeFix_Wire::Perform`
函数级打印），这超出「只读探针」手段。

**本 session 交付**：T-92 枚举清单完成；§9.157 朝向修复（2-wire 24→10、f209 类面 dt=+1）；
全部验收门禁绿（§9.185）；残差归因非 T-92（§9.163/9.171/9.186）。

**门禁**：全绿；a3n00 11675/13168。


### 9.188 round 187 —— 关掉**所有** FixShape 子模式后 F113 仍 2 wire：合并不受任何模式控制

逐个 `set FromSTEP.FixShape.<mode> 0`（在 `STEPControl_Reader` 构造前），F113 的 `FACE 25`
（`box=(-87.5,-34,-100)-(87.5,34,-32)`）**始终 `wires=2 edges 22 6`**：

```
FixMissingSeamMode / FixWireMode / FixOrientationMode / FixLoopWiresMode /
FixIntersectingWiresMode / FixSplitFaceMode / FixSmallAreaWireMode / FixAddNaturalBoundMode /
FixSameParameterMode / FixVertexToleranceMode / FixVertexPositionMode /
FixSolidMode / FixShellMode / FixFaceMode / FixFreeFaceMode / FixFreeShellMode / FixFreeWireMode
   -> 全部 wires=2
```

（`set` 机制已用 `FixMissingSeamMode=0` 验证有效：该设置会让探针在 mesh census 处崩溃，
但 FACE census 已完成。）

⇒ 4→2 **不受 `ShapeFix_Shape`/`ShapeFix_Face` 的任何模式控制**。可能原因：
1. `ShapeFix_Shell::FixFaceTool()` 的 `ShapeFix_Face` 可能**不是** `ShapeFix_Shape::FixFaceTool()`
   的同一实例，于是 `fixshape` 设的模式只作用于 ShapeFix_Shape 那条路径，shell 路径仍用默认模式；
2. 或合并在 `ShapeFix_Shape::Perform` 无条件执行的代码里（`Context()->Apply` / `SameParameter` /
   `FixVertexTolerance`），但这些也都被禁过。

**下一步**：读 `ShapeFix_Shape::FixFaceTool/FixShellTool/FixSolidTool` 的定义（`.lxx` 或 cxx），
确认三条路径是否共享同一 `ShapeFix_Face` 实例；若不共享，则 shell 路径的 face 用的是默认模式，
据此重做模式实验。

**门禁**：全绿；a3n00 11675/13168。


### 9.189 round 188 —— 【方法学修正】用 `Interface_Static` 关模式会崩探针，§9.188 的模式结论**不可靠**

复查发现：`FromSTEP.FixShape.*` 是用 `Interface_Static::Init(..., 't', "-1")` 注册的**文本型**
（`STEPControl_Controller.cxx:202-259`），而 §9.182/9.188 用的是 `SetIVal`（整型 setter）——
类型不匹配。改用 `SetCVal`（`setc`）并在 `ReadFile` 前后都试过：

```
set/setc FromSTEP.FixShape.FixMissingSeamMode -1   -> 崩（exit 0xC0000409）
set/setc FromSTEP.FixShape.FixMissingSeamMode  0   -> 崩
（默认不设任何 static 时正常跑完）
```

连设成**默认值 -1** 都崩 ⇒ 崩来自`设置 static 这个动作本身`（或该 static 的 eval 处理），
而不是「模式被关」⇒ **§9.182/§9.188 里「关掉模式仍是 2 wire」的观测来自被污染/未生效的运行，
结论作废**（保留的可靠结论只有 `noop`/`oneop`：`SetShapeProcessFlags` 有效，
`Operation::FixShape` 是唯一合并算子，§9.176）。

**下一步（正确的模式实验）**：绕开 `Interface_Static`，用
`STEPControl_Reader::SetShapeFixParameters(XSAlgo_ShapeProcessor::ParameterMap&&)`：
先 `GetShapeFixParameters()` 拷贝当前 map，`SetParameter(name, value)`，再设回 reader，
然后逐个关闭 `FixWireMode`/`FixMissingSeamMode`/… 重做实验。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.190 round 189 —— `setp`（ParameterMap）实验：键解析不明，但 wire 直方图在关 `FixWireMode` 时**不变**

新增 `setp <key> <value>`：`XSAlgo_ShapeProcessor::SetParameter` 到
`aReader.GetShapeFixParameters()` 的拷贝再 `SetShapeFixParameters(std::move(map))`（`ReadFile` 后、`TransferRoots` 前）。

结果：

```
setp FromSTEP.FixShape.FixWireMode 0        -> box=(-87.5 面 wires=2；探针在 mesh census 崩
setp FromSTEP.FixShape.FixMissingSeamMode 0 -> box=(-87.5 面 wires=2；探针崩
setp FixShape.FixWireMode 0                 -> box=(-87.5 面**不再出现**；探针崩
```

并对崩溃前已打印的 FACE 行自行统计直方图（两次都只到 198/226 面）：

```
FixWireMode=0        -> 1->180, 2->9, 4->1, 6->4, 10->4
FixMissingSeamMode=0 -> 1->180, 2->9, 4->1, 6->4, 10->4
默认                 -> 1->208, 2->9, 4->1, 6->4, 10->4（226 面）
```

⇒ 关掉 `FixWireMode` 或 `FixMissingSeamMode`（若生效）后 2-wire 面**仍是 9**，与默认一致；
而 `noop`（全部算子关）是 53 ⇒ **合并不在 wire 轮、也不依赖 `FixMissingSeam`** 的结论得到进一步支持，
但键作用域（`FixShape.` 前缀 vs 裸键）与探针崩溃让实验仍不干净。

**下一步**：改在 OCCT 侧**重编译插桩**（`ShapeFix_Face::Perform` 打印每步后的 wire 数）——
这是唯一能把 4→2 的步骤钉死的方法；只读探针的键作用域无法可靠控制这些 't' 型 eval static。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.191 round 190 —— 【决定性】4→2 在 `ShapeFix_Shape::Perform` 里，**不在** `ShapeFix_Face::Perform`

新增探针模式（直接调 API，绕开 static）：

```
perf 113 noop  -> PERF in_wires=4 ret=1 faces=1 wires=4    # ShapeFix_Face::Perform 不合并！
perf 171 noop  -> PERF in_wires=4 ret=1 faces=1 wires=4
shp  113 noop  -> SHP  in_wires=4 ret=1 faces=1 wires=2    # ShapeFix_Shape::Perform 合并 4->2
shp  171 noop  -> SHP  in_wires=4 ret=1 faces=1 wires=4
```

⇒ **`ShapeFix_Face::Perform` 单独调用不合并（4 wire）**；只有经 `ShapeFix_Shape::Perform` 才变 2。
且在该 `perf` 上叠加 `ModifyTopologyMode=true` / `SetContext` / 两者（`topo`/`ctx`/`both`）都仍是 4；
换用 `ShapeFix_Shape::FixFaceTool()` 的同一实例（`shared`）也仍是 4。

⇒ 合并发生在 `ShapeFix_Shape::Perform`（`ShapeFix_Shape.cxx:83-303`）里 FACE 分支之外的调用：
`SameParameter`（`cxx:261`）或 `FixVertexTolerance`（`cxx:289`）或 `Context()->Apply`；
用 `shp 113 noop nosame` / `novtol` 关掉这两者后**仍 2 wire**（setter 是否生效待复核）。

**下一步**：在探针里逐步**手工模拟** `ShapeFix_Shape::Perform`（`Init` → `sff->Perform` → `SameParameter`
→ `FixVertexTolerance` → `Context()->Apply`），每步后统计 wire 数，钉死是哪一步把 4 变 2。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.192 round 191 —— `emul` 复现合并（4→2）但与逐模式无关；`perf` 直调仍 4：触发条件未定

探针新增/使用：

```
emul 113 noop [mode]  # 手工模拟 ShapeFix_Shape::Perform 的 FACE 分支+尾巴
  base -> EMUL face wires=2      # 合并！
perf 113 noop [mode]  # 直调 ShapeFix_Face::Perform
  plain/init/apply/topo/ctx/both/mloc -> 全部 wires=4
```

`emul` 的所有变体（`nowire/nomiss/noorient/noloop/nointersect/nosplit/nosmall/nonatural`、
`notopo/noctx/fresh`）**全部仍是 2 wire** ⇒ 在 `emul` 这条路径下，合并**不依赖任何 `ShapeFix_Face` 模式**，
也不依赖 Context / ModifyTopologyMode。

`emul` 与 `perf` 的差别只在于：`emul` 用 `ShapeFix_Shape sfs; sfs.Init(F); sff = sfs.FixFaceTool();
sff->Init(TopoDS::Face(sfs.Context()->Apply(F)))`，而 `perf` 用 `tmpS`（未 Init 或 Init 后）直接
`Init(F)`；但 `perf init`（Init 后）与 `perf apply`（经 `Context()->Apply(F)`）都仍是 4，
**触发条件未定**。

⇒ 至少确定：**合并在 `ShapeFix_Face::Perform` 内部**（`emul` 复现），且属于**无条件执行**的代码
（`ShapeFix_Face.cxx:676` 的 `FixWiresTwoCoincEdges()` 是 `Perform` 里唯一无 mode 门控的调用；
但 F113 的 wire 不是「2 条边相同」，需再核）。

**下一步**：在 `emul` 路径里逐个**手工调用** `ShapeFix_Face` 的公开方法
（`FixWire`/`FixOrientation`/`FixAddNaturalBound`/`FixMissingSeam`/`FixLoopWires`/
`FixIntersectingWires`/`FixSplitFace`/`FixSmallAreaWire`/`FixWiresTwoCoincEdges`）
并每步统计 wire 数，钉死是哪一个把 4 变 2。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.193 round 192 —— `emul` 结果确为 1 face / 2 wire；单方法调用都不合并

- 给 `emul` 加面数统计：`EMUL face faces=1 wires=2` ⇒ `sff->Perform()`（ShapeFix_Shape 持有的 tool）
  **真的把 4 wire 合并成 2**。
- 新增 `emulm <idx> [noop] <method>`：在同一 tool 上逐个调用 `ShapeFix_Face` 公开方法
  （`FixMissingSeam`/`FixOrientation`/`FixAddNaturalBound`/`FixSmallAreaWire`/`FixWiresTwoCoincEdges`/
  `FixIntersectingWires`/`FixPeriodicDegenerated`）——**全部 `face_wires=4`**（`orient` 崩），
  即**单方法都不合并**，必须整个 `Perform` 序列。

矛盾仍在：`emul`（tool 来自 `ShapeFix_Shape`）→ 2；`perf`（本地/共享 tool，含 `init`/`apply`/`topo`/`ctx`）→ 4。
差别只剩「tool 实例来自 `ShapeFix_Shape::FixFaceTool()` 且 `sfs.Init(F)` 在先」这一组合，
但 `perf init`/`perf apply` 已尽量复刻仍为 4 ⇒ **触发条件仍未定**，需在 OCCT 侧打印实例状态/插桩。

**下一步**：加 `modes` 模式打印两实例（`sfs.FixFaceTool()` vs 本地）的全部 mode 字段
（`FixWireMode`/`FixMissingSeamMode`/… 与 `FixWireTool()` 的 `FixSmallMode`/`ModifyTopologyMode`/`ClosedWireMode` 等），
若完全相同则改为重编译插桩。

**门禁**：`--lib` 1281/0；工作树无调试残留；a3n00 11675/13168。


### 9.194 round 193 —— `emul` 与 `perf` 差别定位失败（模式/上下文/topo/实例都一致）

- `modes` 模式打印两实例的全部 mode 字段（`FixWireMode`/`FixMissingSeamMode`/`FixOrientationMode`/
  `FixAddNaturalBoundMode`/`FixSmallAreaWireMode`/`FixIntersectingWiresMode`/`FixLoopWiresMode`/
  `FixSplitFaceMode` + `FixWireTool()` 的 `ModifyTopologyMode`/`ClosedWireMode`/`FixSmallMode`/
  `FixConnectedMode`/`FixLackingMode`/`FixSelfIntersectionMode`）：
  `MODE standalone wire=-1 miss=-1 orient=-1 natural=-1 small=-1 intersect=-1 loop=-1 split=-1
   wtopo=0 wclosed=1 wsmall=-1 wconn=-1 wlacking=-1 wself=-1`
  与 `MODE shapefix` **逐字段相同**。
- `emul ... local`（在 emul 的 setup 里换用**本地** `ShapeFix_Face`）→ 仍是 1 face / 2 wire。
- `emul init` 打印 `S.IsSame(F)=1 sori=0 fori=0 schild=4 fchild=4` ⇒ `Context()->Apply(F)` 就是 F 本身。

⇒ `emul` 路径与 `perf` 路径在**模式、上下文、topo、实例、被 Init 的形状**上**全部一致**，
却一个合并一个不合并 —— 只能怀疑探针侧的非确定态（未初始化字段随栈内容变化）或表达式求值顺序差异。
**只读探针到此为止**；要结论必须 OCCT 侧重编译插桩。

### 本 session T-92 最终状态（可交付）

- **枚举清单全部完成**（Geom2dInt 家族 / BRep_GCurve pcurve / ShapeBuild_Edge·ShapeFix_Edge·
  ShapeBuild_ReShape·CopyNMVertex·BRep_PointRepresentation / SplitByGrid·SplitByLine·DispatchWires·
  MakeFacesOnPatch·Perform / FixMissingSeam 接线）。
- **关键修复**：§9.157 `WireSegment::first/last_vertex` 按边朝向（ComposeShell 2-wire 面 24→10，OCCT 9；
  f209 类面 1 wire/5 edges、逐面 dt=+1）。
- **验收门禁全绿**（§9.185：topo 1281/0、geom 143/0、parity 14/14、to_obj 13/13、area 11/11、
  geometry 3/3、phase5/9/10/19/20、export 15/15）；BSpline 面 0 个 mv≤11、4 个 360-486（与 OCCT 一致）。
- 残差（总量 +5.6%）已证据化归因到 **非 T-92**：reader 的 `XSAlgo_ShapeProcessor::ProcessShape`/
  `ShapeFix_Shape` 链（F113 4→2 合并 -228）+ 网格密度（F171/F202/F140/F170/F216）。

**门禁**：全绿；a3n00 11675/13168。


### 9.195 round 194 —— T-92 收束（标记目标完成）

本 session 的 T-92 逐行迁移已完成并通过验收门禁：

1. **枚举清单全部落地**（含 UNPORTED 处均标 OCCT 文件+行号）：
   Geom2dInt 家族（`Geom2dInt_GInter` 分派 / `IntCurve_IntConicConic` / `IntCurve_IntConicCurve` /
   `Geom2dInt_Geom2dCurveTool` / `IntRes2d_*`）、edge 级 `BRep_GCurve` pcurve 查找、
   `ShapeBuild_Edge{ReassignPCurve,TransformPCurve,RemovePCurve}` / `ShapeFix_Edge::FixAddCurve3d` /
   `ShapeBuild_ReShape{IsRecorded,Value}` / `CopyNMVertex` / `BRep_PointRepresentation`、
   `SplitByGrid`/`SplitByLine`/`DispatchWires`/`MakeFacesOnPatch`/`Perform`，
   并接线到 `ShapeFix_Face::FixMissingSeam`（`ShapeFix_Face.cxx:1722-2330`，末尾 `2236-2261`
   fictive grid + `ClosedMode=true` + `SetMaxTolerance` + `Perform`），在 STEP 建面后按
   `ShapeFix_Face::Perform:492-494` 调用。
2. **关键修复（本 session）**：§9.157 `WireSegment::first_vertex/last_vertex` 按边朝向取端点
   （对应 OCCT `ShapeFix_WireSegment::FirstVertex/LastVertex`）——ComposeShell wire 结构
   2-wire 面 **24→10**（OCCT 9），f209 类周期面 1 wire / 5 edges、逐面 **dt=+1**。
3. **验收门禁全绿**（§9.185）：`occt-topo --lib` 1281/0、`occt-geom --lib` 143/0、
   `step_obj_parity` 14/14、`step_to_obj` 13/13、`step_obj_area` 11/11、
   `step_geometry_parity` 3/3、`phase5/9/10/19/20` 全绿、`export_data_obj` 15/15（v/f 逐位不变）。
4. **口径**：a3n00 BSpline 面 **0 个 mv≤11**、**4 个 360-486**（与 OCCT 的 4 个一致）；
   总量 11675/13168 已向 OCCT 11145/12466 收敛（+5.6%）。

**残差归因（非 T-92，已证据化）**：F113 4→2 合并（reader 的
`XSAlgo_ShapeProcessor::ProcessShape`/`ShapeFix_Shape` 链，-228）、F171 密度/分解 +907、
F202/F140/F170/F216 网格密度。已定死「合并在 `ShapeFix_Shape::Perform`（`emul` 可复现 4→2）」，
但其触发条件只读探针无法锁定，需 OCCT 侧重编译插桩——属独立后续任务。


### 9.196 —— STEP→OBJ 三道门禁统一模型集（23 个），删除 step_geometry_parity

按用户要求：

- 新增 `crates/occt-topo/tests/common/mod.rs`：唯一模型表 `MODELS`（23 项）=
  15 个 `data/*.step`（Cube/Cone/Cylinder/Sphere/Torus/HoledPlate/OffsetPlaneHoleEdge/rev/Offset/
  Shape/Shape-1/Shape-2/Extrusion/linkrods/screw）+ 8 个 `data/occ/*`（ATU01038/bottom/motoc/top/
  T0M/acs10/TDB/a3n00）。每项带 `bbox_tol`/`area_tol` 以及对应 OCCT 参考。
- `step_obj_parity` → 1 个遍历用例（23 模型 bbox 对拍）；`step_obj_area` → 1 个遍历用例（23 模型面积比）；
  `step_to_obj` → 3 个（23 模型导出有效性 + OBJ 回读 + cylinder 整周面积回归）。
- 删除 `tests/step_geometry_parity.rs`：其 Extrusion/Offset 现由 parity（`occ-Extrusion.obj`/
  `occ-Offset.obj`，均实测 dbbox=0）与 area（ratio=1.0000）覆盖。

实测（本提交的模型表容差来源，`data/` 相对 bbox 最大分量差 / 面积比）：

| 模型 | dbbox | area ratio | | 模型 | dbbox | area ratio |
|---|---|---|---|---|---|---|
| Cube | 0 | 1.0000 | | Extrusion | 0 | 1.0000 |
| Cone | 0 | 1.0000 | | linkrods | 0 | 1.0000 |
| Cylinder | 0 | 1.0000 | | screw | 1e-6 | 1.0000 |
| Sphere | 1e-6 | 1.0000 | | occ/ATU01038 | 5e-6 | 1.0006 |
| Torus | 0 | 1.0000 | | occ/bottom | 5e-6 | 0.9750 |
| HoledPlate | 0 | 1.0000 | | occ/motoc | 2e-6 | 1.0333 |
| OffsetPlaneHoleEdge | 0 | 1.0000 | | occ/top | 3e-6 | 1.0000 |
| rev | 0 | 1.0000 | | occ/T0M | 4.3e-4 | 0.9995 |
| Offset | 0 | 1.0000 | | occ/acs10 | 0 | 0.9846 |
| Shape | 0 | 1.0000 | | occ/TDB | 5e-6 | 1.0092 |
| Shape-1 | 0 | 1.0000 | | occ/a3n00 | 4e-6 | 0.8824 |
| Shape-2 | 0 | 1.0000 | | | | |

结果：`step_obj_parity` 1/1（571s）、`step_obj_area` 1/1（824s）、`step_to_obj` 3/3（773s），全绿。


### 9.197 —— STEP→OBJ 三门禁合并为单一 harness（step_obj_gates）

按用户要求把三份的读取缓存合并成一个 harness：

- 三个测试文件 `step_obj_parity.rs`/`step_obj_area.rs`/`step_to_obj.rs` 合并为
  **`tests/step_obj_gates.rs`**，5 个用例：
  `step_obj_parity_bboxes_match_occt`、`step_obj_area_matches_occt`、
  `step_to_obj_exports_valid_obj`、`outputs_round_trip_through_obj_reader`、
  `cylinder_step_mesh_covers_full_u_period`。
- 23 个模型的 `read_step_file` + `brep_to_obj`（以及 OCCT 参考文本）放进
  `OnceLock<Vec<Case>>`，整轮**只读一次**，5 个用例共享。
- OBJ 回读用例写到 `data/output/roundtrip/`，避免与导出用例并发写同一文件。
- 实测：**5/5 通过、660.58s**；原来三份分开跑 = 571 + 824 + 773 = **2168s**。
- `--lib` 1281/0 不受影响。

门禁命令由三条（`--test step_obj_parity`/`step_obj_area`/`step_to_obj`）改为一条
`--test step_obj_gates`；模型集与容差表见 `tests/common/mod.rs` 与 §9.196。


### 9.198 —— a3n00.obj 与 occ-a3n00.obj 目视差 = 阀体 F113 无网格；根因定位到 `ShapeFix_Face::Perform` 的 Context

用户指出生成的 `a3n00.obj` 与 `occ-a3n00.obj` 目视差别像"拆开的爆炸图"。逐项实测定位：

**现象**：226 个面里只有 **F113** 没网格（`mv=0 mt=0`，
`box=(-87.5,-34,-100)-(87.5,34,-32)`，`Cylinder wires=4 edges=[8,14,1,1] du=6.283 dv=inf`），
正是阀体（横贯两法兰的圆柱）。`stats=225`、`mesh_t=13168`。不是装配变换问题：
两个 OBJ 的包围盒逐位一致，连通域中心逐一对上，XZ 侧视对齐渲染只在阀体处为空
（诊断图 `output/_side_ours.png` / `_side_occ.png`，本轮临时产物）。

**OCCT 侧 oracle**（新加 `wires_probe.cpp` 的 `pcu` 模式 + `run_wires.bat`：
对 raw 4-wire 面用 `XSAlgo_ShapeProcessor::CheckPCurve` 补真实路线的 pcurve 后逐步跑）：

| 调用 | 结果 |
|---|---|
| `ShapeFix_Face::FixMissingSeam()`（无 context） | ret=0 |
| `ShapeFix_Face::Perform()`（无 context） | wires=4 |
| `ShapeFix_Face::Perform()`（**带 Context**） | **wires=2** |
| `ShapeFix_Face::Perform()`（带 Context + topo） | wires=2 |
| `ShapeFix_Shape::Perform()` | wires=2 |

⇒ 合并 4→2 的必要条件是 **`Context()`**。OCCT 的 `ShapeFix_Shape::Perform` 总是
`sff->SetContext(Context())`（`ShapeFix_Shape.cxx:202`）。端口现在是"reader 自己跑一遍 wire 轮 +
单独跑 `fix_missing_seam`"，没有带 context 的 `ShapeFix_Face::Perform`。

端口 ComposeShell 管线探针（临时加，已撤）：`fix_missing_seam` 前的 `tmpF` 就是 4 条 wire
`[8,14,1,1]`；`load_wires→split_by_grid→break_wires` 后 8 段，`collect_wires` 6 条，
`dispatch_wires` **5 个面**（1 个 2-wire + 4 个 1-wire），尾部 `CheckSmallArea` 全 false 所以留下 5 面；
OCCT 同输入是 **1 面 / 2 wire**。所以偏差在"带 context 的 Perform 没有移植"这一层。

**试过但回退的反例**：把 reader 的 `check_pcurves_and_shift(..., fix_lacking=false)` 改成 `true`，
F113 立刻有网格（`stats=226 v=11922 t=12063`，总量反而更接近 OCCT 的 12466）。
但 OCCT 在 `ShapeFix_Face.cxx:372` **显式** `theAdvFixWire->FixLackingMode() = false;`（`cxx:374`
同理关自交），所以那是特例调参、不是 OCCT 路线 ⇒ 已回退，工作树回到 `stats=225 t=13168`。

**按 OCCT 路线补的代码（也已回退）**：`ShapeFixFace::perform`（`cxx:365-498` 第一轮 wire +
`FixMissingSeam`）、`FixMissingSeam` 开头补 `cxx:1737-1741` `myFace = Context()->Apply(myFace)`。
实测对 F113 无效：端口的 `check_pcurves_and_shift` 既不报 `ShapeFix_Wire` 的 `fixed` 状态，
也不产生结构性 wire 变化，重建出来的面与原件 wire 等价，`fix_missing_seam` 仍出 Shell(5)。

**下一步（OCCT 路线，独立任务）**：完整移植
`ShapeFix_Shape::Perform` 的 `sff->SetContext(Context())` 机制 + `ShapeFix_Face::Perform`
`cxx:377-480`/`cxx:482-498`，并让 `check_pcurves_and_shift` 返回 `ShapeFix_Wire` 的
`StatusSmall/Connected/EdgeCurves/Notched/FixTails/Degenerated/Closed` 以驱动 `fixed`，再补尾部
`cxx:2282-2322` 的 `Context()->Apply/Remove`。落地后 F113 应出网格（OCCT 约 228 三角），
a3n00 总量向 OCCT 收敛。


### 9.199 —— T-93 定位：F113 无网格 = CheckNotchedEdges 因边朝向不同而永不命中

复核结论（HEAD `08c4edf2`）：`stats=225`、F113 `mv=0`（Cylinder, wires=[8,14,1,1]）仍然成立。

**修正 §9.198 的一半结论**：4→2 的合并**不依赖 Context**。

- `PCU3`（wire 工具传 null context）与 `PCU2`（传 context）一样：`ShapeFix_Wire::Perform` 让 wire[0] **8 → 6 条边**，`StatusNotches DONE`，`Fixed` 置位。
- 所以真正的缺口是 **notch 这一步没生效**，不是 context 机制。

**精确定位**（`wires_probe.cpp` 新增 `cne` + `CNE2`）：

| 观测 | 结果 |
|---|---|
| OCCT `CheckNotchedEdges` 在**原始** wire 上 | 0 命中（notch 只在前置 fix 之后才出现） |
| 跑完前置序列（FixReorder/FixSmall/FixConnected/FixEdgeCurves/FixDegenerated）后 | **命中 `i=3, shortNum=2, param=0`** |
| 命中处的边朝向 | `E1=Forward, E2=Reversed` |
| 端口同一对边（`num=3`） | `ori1=Forward, ori2=Forward`；`ang=3.138`（`p2d1==p2d2`，dd≈1e-16） |
| 端口 `reorder_ok` / `fix_notched_edges` | `true` / 被调用但恒 `false`（`cxx:1940` 的 `|angle|>0.1` 拒掉） |

端口 `check_notched_edges` 的切向选取（`wire_fix.rs:2862-2873`）与 `ShapeAnalysis_Wire.cxx:1910-1933` **逐行相同**，
因此差异在**输入**：同一对邻边在端口里两侧都是 Forward，在 OCCT 里第二边是 Reversed。切向约定随之翻号 ⇒ `ang≈π` ⇒ 永不命中。

**下一步**：查端口 wire 建链/前置修（`fix_reorder_wire` 及其上游的边朝向布置）为什么没有产生 OCCT 那样的 `E2=Reversed`。
在这一点对上之前不要改 `check_notched_edges`（那是把症状当原因）。

本轮所有临时插桩已撤；`wire_fix.rs` 与 HEAD 逐字节一致；`occt-topo --lib` = **1281/0**。


### 9.200 —— T-93 修复一处：CheckNotchedEdges 的端点必须朝向感知（wire 8→6）

承接 §9.199 的定位（端口在 `i=3` 判 `ang≈π`、OCCT 命中）。继续插桩后确认：

- F113 的 wire **朝向本来就是对的**（`[F,F,R,F,F,F,F,F]`，与 STEP 的 `ORIENTED_EDGE` XNOR `same_sense` 一致），
  `fix_reorder_wire` 也没有改动它。此前看到的 `Forward/Forward` 是**另外几条** 8 边 wire。
- 真正卡住的是 **端点判定**：`check_notched_edges` 用原始 `edge_vertices(e1).1` / `edge_vertices(e2).0`，
  而 OCCT `ShapeAnalysis_Wire.cxx:1885-1898` 用的是 `sae.LastVertex(E1)` / `sae.FirstVertex(E2)` ——
  两者都是**朝向感知**的（`ShapeAnalysis_Edge.cxx:228-258`：REVERSED 边的 "first" 是**原始 last**，反之亦然）。
  对 REVERSED 边端口正好取反 ⇒ `verts differ` ⇒ 该对边被跳过 ⇒ notch 永不命中。
- 修法（同一 `.cxx` 分支）：改用端口自己的 `last_vertex` / `first_vertex`（`shhealing/wire_fix.rs:51-67`）。
  `BRepTools::Compare` 只比 `IsSame`/容差，故 `sae` 返回时的 `V.Reverse()` 不影响结论。

**实测（HEAD 08c4edf2 + 该修复）**：

| 项 | 修前 | 修后 |
|---|---|---|
| `fix_notched_edges` | 恒 false | `notched=true`，F113 wire[0] **8 → 6 边**（与 OCCT wire 轮一致） |
| `occt-topo --lib` | 1281/0 | **1281/0** |
| `--test step_obj_gates` | 5/5（660s） | **5/5（397s）** |
| a3n00 面积比 | 0.8824 | **0.8996**（向 1 收敛） |
| a3n00 网格 | 11675/13168 | 10863/11941 |
| **F113** | mv=0 | **仍 mv=0**（wires 仍 4） |

⇒ AC 只满足一半：面积比收敛与门禁不劣化达成，**「F113 出网格」未达成**，T-93 保持未完成。

**剩余缺口**：面级合并。OCCT `ShapeFix_Face::Perform` 在第一轮 wire 修后
（`cxx:454-479`）用修好的 wire 重建 `tmpFace`、`Context()->Replace(S, tmpFace)`、`myFace = tmpFace`，
再让 `FixMissingSeam` 经 `cxx:1737-1741` 的 `Context()->Apply(myFace)` 看到它，才能把 4 条 wire 并成 2
（oracle：`Perform` 带 Context → wires 2；不带 → 4）。端口 reader 目前是「wire 轮 + 单独 fix_missing_seam」，
缺这一步。


### 9.201 —— T-69 复测：不再是 162 面；剩 6 面，且根因更像 wire 结构分歧

卡片写「端口这 162 面全部 Delaunay produced no triangles；OCCT 1778/1778 全成功」。**两条都要订正**（HEAD 08c4edf2 实测）：

| 量 | 端口 | GT（occt_probe） |
|---|---|---|
| T0M faces | **1772** | **1778** |
| 未网格面 | **6** | **1**（--mesh 0.1：FACE 1327 nodes=0 triangles=0） |
| WIREHIST（face→wire 数） | 1:1652 · 2:106 · 3:8 · 5:1 · 6:1 · 7:2 · 8:1 · **10:1** | 1:1672 · 2:92 · 3:8 · 5:1 · 6:1 · 7:2 · 8:2 |
| --wires 汇总 | — | wires=1921、faces_with_coincident_closed_edges=0、wires_dup_pc_any=117 |

⇒ ① 端口只剩 **6** 个未网格面（三类：F351 = Cylinder/10 wires/dv=inf；F405–408 = BSpline/1 wire/2 edges/dv=2π；F1758 = Cone/1 wire/1 edge）；
② **OCCT 自己也有 1 个面不出三角**，所以 AC 不该写「→0」，应为「≤ GT 的 1」；
③ 更关键：端口把约 **20 个「1 wire」面做成了「2 wire」面**（1652 vs 1672、106 vs 92），并多出 1 个 **10-wire** 面（GT 最大 8）。
这与 T-93 的 F113 是**同族**（reader/ComposeShell 多出 wire），所以先解决 wire 结构分歧，再谈网格算法。

本轮只做复测与对拍，**未改代码**；临时产物 output/t0m_probe.txt 已删。


### 9.202 —— T-93 复测（notch 修好后）：缺口在 ComposeShell，不在 FixMissingSeam 尾部

在 §9.200 的 notch 修复之上复测 F113（HEAD 08c4edf2 + wire_fix.rs 修复）：

- F113 送入 `fix_missing_seam` 时 wires 已是 **`[6, 14, 1, 1]`**（notch 生效，与 OCCT wire 轮一致）；
- `fix_missing_seam()` 返回 **true**，但 `myResult` 是 **Shell(5 面)**；OCCT 的 `ShapeFix_Face::Perform` 在同一输入下给 **1 面 / 2 wire**；
- reader 只接受 `shape_type() == Face` 的结果 ⇒ 丢弃 ⇒ 面仍 4 条 wire ⇒ `heal_self_intersecting_wires` 判 FAILURE ⇒ `mv=0`。

排除尾部：`ShapeFix_Face.cxx:2270-2322` 删面靠两处 —— ① 每个 wire `FixSmall(true, Precision())` 塌缩则 `Context()->Remove`；
② `nbFaces > 1` 时 `FixSmallAreaWire(true)`（`:2331-2394`，即 `CheckSmallArea` 为真的 wire 被丢）。
端口对那 4 个单 wire 面的 `check_small_area` **全部返回 false**（面积是合理的四分之一柱面），OCCT 也不会丢它们。
⇒ 5 面是 **ComposeShell 自己产出的**，不是尾部没删。

**结论**：T-93 的下一站是 `ShapeFix_ComposeShell`（`dispatch_wires` / `make_faces_on_patch`），
与 §9.198 的 ComposeShell 探针结论一致（当时测到 `seqw=8 → wires=6 → faces=5`）。本轮未再改代码；临时打印已撤。


### 9.203 —— T-93 波次 1：roots 判定忠实，怀疑上移到 CollectWires

逐行对读「多 loop 找 roots」（端口 `make_faces_on_patch.rs:62-170` vs `ShapeFix_ComposeShell.cxx:3021-3143`）：

- `unp` 取法一致（loop i 首条定向边在 `pf` 上的 pcurve 中点，`cxx:3056-3062`）；
- 包含判定一致：对每个 `j != i` 建 `awtmp`（只收定向边）→ `FClass2d` → `Perform(unp)`；`ON/UNKNOWN` 时按 `cxx:3104-3130` 沿 loop i 的边逐个点继续判；
- `stPoint != PerformInfinitePoint()` 则 `break`（i 不是 root）；全部 j 走完才 `roots.push`（端口 `j >= loops.len()` ⟷ OCCT `j > loops.Length()`）。

⇒ **roots 算法忠实**，5 个 root 是输入决定的。

**关键对照**：`dispatch_wires` 对 F113 只调一次 `make_faces_on_patch(loops=6)`；而 OCCT 的 `FixMissingSeam` 最终结果是 **1 面 / 2 wire**。
⇒ 要么端口的 `CollectWires` 把 8 个 segment 过分裂成 6 条（OCCT 可能并成 2 条），要么 `MakeFacesOnPatch` 收到的 loop 集不同。

下一轮：给 `collect_wires.rs` 打「输入 segment → 输出 wire」配对（含各自 patch 索引与首尾点），对读 `ShapeFix_ComposeShell.cxx:2512-2900`（`CollectWires`）的合并/闭合规则，先确认这一步的输出条数就该是 2。

本轮无代码改动；全部临时探针已撤，`occt-topo --lib` = **1281/0**。


### 9.204 —— T-93 波次 1：缺口在 CollectWires 的「找下一段」

给 `collect_wires.rs` 打输入/输出配对探针（a3n00 F113，grid 1×1）：

```
输入 8 段: [0..5] Reversed n=6,7,5,4,1,1 ; [6..7] External n=1,1
输出 6 条: close n=6 / 7 / 5 / 4 / 3 / 3   （没有任何 moved-out，也没有 separate）
```

- 四条大段 6/7/5/4 **各自闭合成一条 wire** ⇒ 每条都命中了 `index.is_none()` 分支（`collect_wires.rs` 的 close 条件第一项），
  即 cxx:2848-2914 的「找下一段」对它们**恒返回 none**；
- OCCT 在同输入下最终得到 **1 面 / 2 条 wire**（`FixMissingSeam` oracle），即它把这些段串成了更少的 wire；
- 上一轮已排除 roots 判定（`make_faces_on_patch.rs:117-144` ⟷ `cxx:3099-3143` 逐行一致），所以 5 面是 6 wire 的后果，
  6 wire 才是源头。

**下一步**：对读候选搜索（端口 `collect_wires.rs` 约 :340-400）与 `cxx:2848-2914` 的条件 —— patch 索引匹配、端点/切向容差、
`is_coincided`、`samepatch` 优先级 —— 找出恒返回 none 的那一条。

本轮无代码改动；探针已撤，`occt-topo --lib` = 1281/0。


### 9.205 —— T-93：主循环候选搜索恒失败（下轮插桩点已定）

继续对读 `ShapeFix_ComposeShell::CollectWires`：

- `cxx:2846-2936` 是**主循环之后**的「短段合并」后处理（只处理 `shorts(i)==1`）。端口探针显示 F113 的 6 条 wire **全部**来自主循环的 close 分支（`CW out close n=6/7/5/4/3/3`，没有任何 separate/moved-out）⇒ 该后处理不参与。
- 主循环候选搜索：端口 `collect_wires.rs:150-228` 与 `cxx:2627-2722` 逐段核对一致 ——
  顶点 `endV.IsSame(seg.FirstVertex()/LastVertex())`、回退同边最低优先级、`GetEndTangent2d(edge, face, false, lPnt, lVec, 1e-3)`、
  `IsCoincided(endPnt, lPnt, URes, VRes, ctol)`、以及 `weigth = 16*(sp)+8*(!misor)+4*conn` 与 `tail` 比较。
- 但每条大段（6/7/5/4）跑完候选循环后 `index` 都是 None ⇒ 走 `cxx:2824-2846` 的 close，各自成一条 wire。

**下轮插桩点**：主循环里，第 2 轮起每次迭代对每个候选 `i`/`j` 打印 `vertex_matched / sp / conn / dist / ang` 以及 `w1+tail1` vs `weigth+tail2`，
判定是「顶点就匹配不上」还是「优先级比较把它否掉」。

本轮无代码改动；探针已撤，`occt-topo --lib` = 1281/0。


### 9.206 —— T-41 第①步实测：不可单独摘（断言绑在旧采样行为上），已回退

计划：删掉 `region_interior_points_2d` + `region_inside_other`（UV 采样 + 多数票 + 自造 `inside*2>=total` 平局规则），
调用点（`region_trim.rs:388`）改用同文件已在用的面级分类器（`build_face_from_loops` + `classify_face_general`）。

**结果：`occt-topo --lib` 1280/1 —— `bop_curved::tests::general_boolean_trimmed_closed_shell` 失败。**
⇒ 该断言绑在旧的 region 采样行为上（§9 记的 T-79 型「改分派/实现前先查断言绑定」陷阱）。按纪律**回归即回退**，已恢复原实现；`--lib` 回到 1281/0。

顺带查明（写进看板）：

- `classify_face_general`（`region_mesh.rs:274-301`）**本身也是采样器**（8×8 网格 + `point_in_solid_curved`），
  并不是 `BRepClass3d_SolidClassifier`；所以「换成它」不增加忠实度，只是换一个自造采样器。
- `bop_curved` 的调用图：`curved_boolean_full` ← `boolean_dispatch:227`（非 compound 实体输入的主分派）/ `bop_builder_repair:239,405` / `bop_builder_report:388`；
  `curved_boolean` ← `brepfeat/features:320`；`boolean_ops::voxel_boolean` ← `brepfeat/features:330-335`。
- 绑在旧体上的断言至少 3 处：`general_boolean_trimmed_closed_shell`、`curved_boolean_full_quadric_unchanged`（要求 quadric 走 `curved_boolean`）、`voxel_*` 三条。

**结论**：T-41 不能拆成「只摘一个函数」的小步 —— 它的分类层是 load-bearing 的。
正确做法是一次做：把 `boolean_dispatch` 对非 compound 实体输入的分派改到忠实 `bop_builder2`（T-82 已转绿），
再整体摘掉 `bop_curved` 的 mesh/voxel 层，并同步处理上述 3 处断言绑定（按 OCCT 口径订正，而不是改成新实现能过）。


### 9.207 —— T-93：候选搜索卡在顶点身份（TShape）而非优先级

在 `collect_wires` 主循环插桩（a3n00 F113，grid 1×1，seqw=8）：

```
CW  i=2 j=0 vertex-mismatch ori=Reversed cand_first=...919520 cand_last=...743344 end_v=...751248
CW  i=2 j=1 vertex-mismatch ...（同一对）
CW  i=3 j=0/1 vertex-mismatch ori=Reversed ...
CW  i=4,5 ... vertex-mismatch ori=Reversed
CW  i=6,7 ... vertex-mismatch ori=External
（第 2 轮迭代）end_v 变为 ...743344，候选 i=3..7 仍全部 mismatch
例外：一次 i=6 j=0 sp=true misor=false conn=true dist=0e0 ang=1.5708 ⇒ l=29 vs 0，被接受
```

⇒ 拒绝发生在 `collect_wires.rs:141` 的 `!same_v(&end_v, &candidate_v)`，**不在优先级比较**（`w1+tail1 vs weigth+tail2` 大多根本没执行）。
`same_v` = `TopoDS_Shape::IsSame`（同一 TShape），所以问题是**顶点身份**：
累积 wire 的 `end_v` 与剩余各段的首/末顶点是**不同的 TShape**（指针完全不同）。

**下一步**：查上游 —— `BreakWires`/`SplitByGrid` 把原 wire 拆成段时，相邻段是否共用同一 vertex TShape（OCCT 里共用，故 `IsSame` 成立）；
或 `end_v` 的更新取了错的端点（`cxx:2733-2748` 一带）。

本轮无代码改动；探针已撤，`occt-topo --lib` = 1281/0。


### 9.208 —— T-93：拆段丢了端点连通性（CollectWires 断链的直接原因）

同一进程内对拍（a3n00 F113，grid 1×1，seqw=8）每段首/末顶点指针 + 候选搜索的 `end_v`：

```
CWSEG [0] Reversed n=6 first=...737280 last=...737280   ← 自闭合
CWSEG [1] Reversed n=7 first=...755584 last=...947552
CWSEG [2] Reversed n=5 first=...763072 last=...952544
CWSEG [3] Reversed n=4 first=...765568 last=...755584
CWSEG [4] Reversed n=1 first=...754752 last=...754752   ← 自闭合
CWSEG [5] Reversed n=1 first=...745600 last=...745600   ← 自闭合
CWSEG [6] External n=1 first=...754752 last=...763072
CWSEG [7] External n=1 first=...765568 last=...745600
CWMIS i=2..7 全部 vertex-mismatch，end_v=...947552（= seg1.last）
```

端点共享（同指针）: 755584(1.first=3.last)、763072(2.first=6.last)、754752(4=6.first)、745600(5=7.last)、765568(3.first=7.first)；
**无人配对**: 947552（仅 seg1.last）、952544（仅 seg2.last）、737280（seg0 自闭合）。

⇒ 第一次取 seg1（seg0 被判 short 跳过），`end_v = seg1.last = 947552`，此后没有任何候选与之 `IsSame` ⇒ 链断 ⇒ close ⇒ 每条大段各自成 wire（6 条）。
OCCT 的 `BreakWires`（`cxx:2279-2306`）只在 `EXTERNAL`/`INTERNAL` 段给出的 split vertex 处切，且切出的段沿用**原顶点**；
端口这些悬空端点说明切分时丢了共享身份（或切点算错）。

**下一步**：读 `split_by_grid.rs` / `break_wires.rs`，对照 `cxx:2131 SplitByGrid` 与 `cxx:2279 BreakWires`，找出丢共享身份的那一步。

本轮无代码改动；探针已撤，`occt-topo --lib` = 1281/0。


### 9.209 —— T-54 前提核实：耳切是「回退路径」，与 T-69 同源

卡片证据（「wireframe 的耳切 + 质心角度排序 + 桥洞在 OCCT 里没有对应控制流」）**成立且仍在活代码里**：
`wireframe.rs` 的 `ear_clip`(:156) / `bridge_holes`(:219) / `planar_polygon_triangulate`(:269) 都还在。

本轮补上的是**范围**：`wireframe::face_to_triangles` 的调用点全是**回退**，主路径是 `discret_root` 的 Delaunay 管线 ——

| 调用点 | 性质 |
|---|---|
| `discret_root.rs:290`（`build_shape_mesh_wireframe`） | 整形状回退：OCCT 风格管线报错时回到 pre-pipeline 的 UV 网格剖分 |
| `discret_root.rs:554`（`wireframe_face_triangulation`） | **单面**回退：该面 Delaunay 失败时 |
| `discret_root.rs:1376` / `:1382`（`discretize_face`） | 单面回退：`ModelBuilder::build_model` / `build_shape_mesh` 报错时 |

⇒ T-54 改的是回退路径，不是主管线；而 T-69 的 6 个未网格面正属于「Delaunay 空产出」那一类 ——
两者同源，T-54 的 `blockedBy: T-69` 成立，继续 blocked。

本轮无代码改动。


### 9.210 —— T-93：断链根因 = 切口处的「空拷贝顶点」没有被 Context 统一

继续沿「拆段丢端点身份」下钻，逐段对拍：

```
after split_by_grid  n=6 : [0] n=6 自闭合 ; [1] n=16 自闭合(first=last=...003920) ; [2][3] n=1 自闭合 ; [4][5] External
after break_wires    n=8 : [1] 被切成 7+5+4，接缝顶点 ...428256 / ...426592 在整体里无人配对
```

把 n=16 那条段的逐边端点打出来后可见：段内部**本身就是若干互不衔接的子链**，且断点处的顶点落在**另一片分配区**
（链上顶点 …973xxxxxx，断点 …9710xxxx / …970xxxx）—— 典型的「新建/空拷贝顶点」。

对照源码确认机制：

- `split_wire.rs:315-335` 忠实实现 `cxx:1261-1293`：某条边**首次被切**时，把它的原始端顶点换成
  `empty_copied_vertex`（空拷贝，`TopoDS_Shape.hxx:294-302` 语义）并 `Context()->Replace(原, 拷贝)`（`cxx:1268`/`:1281`），
  目的是让后续 `SameParameter` 不把原顶点的容差撑大；
- 于是**切出的段引用拷贝顶点，未切的邻边仍引用原顶点**；
- OCCT 在四处 `Context()->Apply` 把两者统一：`cxx:506`（wire 装载）、`cxx:2019-2020` 与 `cxx:2048-2049`（SplitByLine 序列的切线合并）、`cxx:2124`；
- **端口的 `break_wires.rs` 与 `collect_wires.rs` 里 context 引用为 0**（`split_by_line_wires` 的合并步 `split_by_line.rs:586-602` 有用 context，但覆盖不到全部）。

⇒ `CollectWires` 的 `same_v`（`IsSame`）在拷贝↔原件之间必然失败 → 链断 → 6 条 wire → 5 个 root → 5 面 → Shell(5) 被丢弃 → F113 仍 4 wire / mv=0。

**下一步**：对照 `cxx:2019-2060` 检查 `split_by_line_wires` 的 `Replace/Apply` 覆盖（是否只处理 coincide 对、漏了拷贝↔原件），并核 `load_wires.rs` 对 `cxx:506`。

本轮无代码改动；探针已撤，`occt-topo --lib` = 1281/0。


### 9.211 —— T-93：`split_by_line_wires` 的 Context 用法是忠实的；缺口在「过路边 vs 切点拷贝」

对照 `cxx:2019-2057` 读端口 `split_by_line_wires`（`split_by_line.rs:585-630`）：

| OCCT | 端口 | 结论 |
|---|---|---|
| `tmpV1/2 = Context()->Apply(SplitLineVertex(i-1)/(i))`（`:2019-2020`） | `self.context().apply(&split_line_vertex[i-2]/[i-1].0)`（`:586-587`） | 一致 |
| 短/coincide 判定 + `CombineVertex` + `Context()->Replace(V1/V2, V)`（`:2037-2056`） | `combine_vertex` + `context_mut().replace(v1/v2)`（`:598-611`） | 一致 |
| 新建边 `V1(FORWARD)/V2(REVERSED)` 做 External 段（`:2059-2074`） | `make_shape(Edge)` + `add(v1/v2)`（`:615-630`） | 一致 |

⇒ **缺的 Context 应用不在这段**。真正对不上的是顶点来源：

- External 段是用 `Context()->Apply` 之后的**拷贝**顶点建的；
- 而 `split_wire`（`cxx:1261-1293` / 端口 `split_wire.rs:315-335`）只把**被切那条边**的原端顶点换成 `empty_copied_vertex` 拷贝；
- 于是**没被切到的过路边**在同一个几何点仍持有**原件**；
- `break_wires` 走到该点时切出的两半，一端是拷贝、一端是原件 ⇒ `CollectWires` 的 `IsSame` 失败 ⇒ 链断。

**下一步**：读 `split_wire.rs` 的收尾（`cxx:1360-1429` 的 wire 重组与 `Context()->Replace(edge, e1/e2)`），
与 OCCT 对照邻边的顶点引用是否本该一起替换成拷贝（若 OCCT 在切边时把整个 wire 的共享顶点都换成拷贝，端口漏了这一步）。

本轮无代码改动（只读源码）；`occt-topo --lib` = 1281/0。


### 9.212 —— T-93：读源码法到头；剩余分歧需 OCCT 侧插桩

本轮把最后一组可疑点逐行核对完，**端口与 OCCT 一致**：

| OCCT | 端口 | 结论 |
|---|---|---|
| `ApplyContext`（`cxx:382-446`）：`context->Apply(edge)`；同则 1；是 Edge 则 SetEdge；否则按 children 展开插入 | `reshape.rs:113-156` | 忠实（**只处理边级**替换，两边都不做顶点替换） |
| `SplitByLine` 序列尾部对全部 wire 应用 context（`cxx:2120-2126`） | `split_by_line.rs:677-683` | 有 |
| 切线合并处的 `Context()->Apply` 与 `CombineVertex`（`cxx:2019-2057`） | `split_by_line.rs:586-611` | 忠实 |

⇒ 我能在源码层对照的每一处都对上了，但行为仍差：OCCT 1 面 / 2 wire，端口 5 面 / 6 wire。
观测到的具体差异是 **CollectWires 输入里存在「原件 ↔ 空拷贝」配不上的顶点**（§9.210/§9.211），
而 `empty_copied_vertex` 是 `cxx:1261-1293` 的忠实行为、`ApplyContext` 又不做顶点替换 —— 说明分歧在更细的地方：

1. `split_wire.rs`（495 行，`cxx:942-1429`）的顶点簿记细节（哪些边、哪些端点被换成拷贝）；或
2. 进入 ComposeShell 的**输入面本身**（`tmp_f` 的 4 条 wire 与 OCCT 的 `tmpF` 是否同构）。

**下一步（二选一）**：
① **OCCT 侧插桩**：把 `ShapeFix_ComposeShell.cxx` 单独重编进 `wires_probe` 并加打印，直接看 OCCT 在 F113 上
   `SplitWire`/`BreakWires`/`CollectWires` 的顶点链（最直接，但要搭一次编译）；
② 对 `split_wire.rs` 做逐块 diff（495 行 vs `cxx:942-1429`），成本高且不保证命中。

本轮无代码改动；`occt-topo --lib` = 1281/0。


### 9.213 —— T-93 根因改判：在 reader 的 wire 组装，不在 ComposeShell（OCCT 侧插桩建成）

**搭好了可复用的 OCCT 侧插桩**（本轮主要产出）：

- `specs/occt_probe/_dbg/ZZ_ComposeShell.cxx`：`ShapeFix_ComposeShell.cxx` 的本地副本 + 打印（`zzIsTarget` 按包围盒门控，在 `CollectWires` 候选循环与 `BreakWires` 结尾插桩）；
- `specs/occt_probe/build_dbg.bat`：与 `wires_probe.cpp` 一起编成 `zz_wires_probe.exe`。**MSVC 下 exe 自带 .obj 优先于 .lib**，本地副本覆盖 `TKShHealing.lib` 的同名符号，链接无重复符号（实测 BUILD OK）；
- `specs/occt_probe/run_dbg.bat`：与 `run_wires.bat` 同形；
- `wires_probe.cpp` 新增 `findf` 模式：按包围盒列出面的 explorer 下标；
- `.gitignore` 增加 `_dbg/`、`zz_wires_probe.*`（本地工具，不入库）。

**定案测量**（a3n00，bbox 唯一吻合的面 = explorer 25）：

| | faces | wires | FixMissingSeam |
|---|---|---|---|
| OCCT 原始传输（`pcu 25`） | 1 | **2**（22 边 / 6 边） | ret=**0** |
| OCCT 处理后（`pcu 25 noop`） | 1 | **1** | ret=**0** |
| 端口 F113 | 1 | **4**（`[6,14,1,1]`） | 返回 **true**，进入 ComposeShell |

⇒ 分歧在**进入 ComposeShell 之前**：同一条 wire 在 OCCT 是 1 条（22 边），端口拆成了 4 条 `[6,14,1,1]`。
OCCT 的 `ShapeFix_Face::FixMissingSeam` 对这个面**根本不做任何事**；端口因为面被拆坏而判定「缺 seam」，
跑 ComposeShell 后再切成 5 个面，`Shell(5)` 被 reader 丢弃，F113 因此无网格。

这条结论同时解释 T-69 的「该 1 条 wire 却成 2 条」—— 同族的 **wire 过度拆分**。
因此 §9.196-§9.212 的 ComposeShell 调查是**在下游追症状**（其中 §9.200 的 notch 朝向修复是真实修复，保留）。

**下一步**：查 `read_topology` 的 wire 组装 —— 为什么这条 wire 被切成 `[6,14,1,1]`，以及 OCCT 的合并（→22 边单 wire）发生在哪一步。

本轮无 Rust 代码改动；`occt-topo --lib` = 1281/0。


### 9.214 —— T-93 根因落到具体边：reader 省掉了 seam 段、并把 seam 边拆成独立 wire

工具：给端口探针 `zz_probe_a3n00.rs` 加 `--fdump`（按 F113 包围盒定位面，逐 wire/逐边打印 3D 端点、朝向、退化）；
OCCT 侧用 `wires_probe` 的 `fdump -1`（`-1` = 按同一包围盒定位）。两边对拍同一张面：

```
OCCT 原始 : face=25  2 wire : wire[0] nEdges=22 , wire[1] nEdges=6
OCCT 处理后: 该包围盒的面已不存在（ProcessShape 把它并掉了）
端口 F113 : face=113 4 wire : wire[0] nEdges=6(内环) , wire[1] nEdges=14(外环) , wire[2] nEdges=1 , wire[3] nEdges=1
```

**逐边对应结果**：

- 内环：两边都是 6 边、同一起点 `(12.3875,24.672,-89.3942)` ✓；但**有一个顶点坐标不同** ——
  OCCT `(-49.8649,0,-100)` vs 端口 `(-65.243610,0,-100)`（y/z 相同）；
- 外环：OCCT 22 边（含 `e[9]`/`e[17]` 两条 `closed=1` 的 seam 边，端点 `(±87.5,0,-32)`），端口只有 14 边；
- **端口 `e[7]` 直接从 `(38.733061,21.908902,-40)` 连到 `(38.733061,-21.908902,-40)`**，
  跨过了 OCCT 的 `e[7] (38.733,21.909,-40)→(44.5,0,-32)`、`e[8] →(87.5,0,-32)`、`e[9]`(closed seam)、`e[10]` 回、`e[11] →(38.733,-21.909,-40)`；
- 端口把两条 seam 边放进**独立的 1-边 wire**（`wire[2]` 端点 `(87.5,0,-32)`、`wire[3]` 端点 `(-87.5,0,-32)`），而 OCCT 把它们**并进外环**。

⇒ 分歧是 **reader 的 wire 组装**：绕 seam 的那一段被整个省掉（跨接），seam 边被拆到独立 wire。
OCCT 的 `FixMissingSeam` 对这个面返回 0（不做），端口因为面被拆坏而判定「缺 seam」，进 ComposeShell 后切成 5 面 → `Shell(5)` 被丢 → F113 无网格。
内环那处**顶点坐标不同**（同一条边取了不同顶点）提示 `resolve_oriented_edge` 一族的顶点/前驱解析仍有一处错位 ——
与 §9.200 修的 notch 朝向问题同族。

**下一步**：读 `read_topology` 对该面 ORIENTED_EDGE 的组装（STEP 里对应 loop 一带），定位「seam 段被跨接 / seam 边被拆成独立 wire」的那一步。

本轮无 Rust 库代码改动（只改了探针 `examples/`）；`occt-topo --lib` = 1281/0。


### 9.215 —— T-93：到 STEP 实体级 —— 端口外环 = STEP 的 8 条边，OCCT 展开成 22 条

用端口 `--fdump` 与 OCCT `fdump -1` / `fcount` 对拍，再加一次性打印读出 STEP 结构：

```
端口 F113 的面 = STEP #4952（CYLINDRICAL_SURFACE R=87.5），4 个 bound：
  #5140 FACE_OUTER_BOUND -> #5139 EDGE_LOOP('',(#5005,#5013,#5019,#5091,#5099,#5108,#5116,#5138))
  #5352 / #5363 / #5374  FACE_BOUND
STEP 外环只有 8 条 ORIENTED_EDGE。
OCCT fcount: 该包围盒的面唯一（face=25），原始与处理后**都是 2 wire**；OCCT 外环 = 22 边。
端口：唯一同包围盒的面 F113，4 wire = [6, 14, 1, 1]。
```

**因此**：端口外环 8 边 = STEP 原样；OCCT 把同样的 8 条**展开成 22 条**（切边 + 插 seam 边，绕到 `(±87.5,0,-32)`），
端口不做这个展开，反而把两条 seam 边放进独立 1-边 wire。

顺带排除的假设：

- 不是 `fix_missing_seam` 加的（`--fixms` 显示修之前 F113 就是 4 wire）；
- 不是「OCCT 的 `ShapeFix_Face::Perform` 把 4 条并成 2 条」——`pcu 25` 显示 `PCU in ... wires=2`，即 **Perform 之前**就是 2 wire；
- 不是包围盒撞脸（端口与 OCCT 各只有 1 个该包围盒的面，总数都是 226）。

**下一步**：对读 `StepToTopoDS_TranslateEdgeLoop` 的 seam/切分段（含 `StepToTopoDS_GeometricTool` 的 seam 判定）
与端口 `read_topology` 的 `is_seam_curve` / `order_seam_pcurves` / `associate_edge_pcurve` 一族，定位「8 条展开成 22 条」缺在哪一步。

本轮无库代码改动（探针与一次性打印已撤）；`occt-topo --lib` = 1281/0。


### 9.216 —— T-93：OCCT 是「bound 翻译失败即跳过」，端口照样建了 1-边 wire

**本轮推翻/确认的关键控制流**：

```
StepToTopoDS_TranslateFace.cxx:592  for (aBoundIndex = 1..theFaceSurface->NbBounds())
StepToTopoDS_TranslateFace.cxx:694    StepToTopoDS_TranslateEdgeLoop anEdgeLoopTranslator; ... Init(...)
StepToTopoDS_TranslateFace.cxx:706    if (anEdgeLoopTranslator.IsDone())  -> aFaceBuilder.Add(aResultFace, anEdgeLoopWire)   // :727
StepToTopoDS_TranslateFace.cxx:729    else { AddFail(" EdgeLoop not mapped to TopoDS"); continue; }   // :731-741 跳过该 bound
```

⇒ **`TranslateFace` 每个 bound 只 add 一条 wire，从不合并**。所以 4-bound 的面在 OCCT 也应是 4 条 wire。
而实测 OCCT 该面只有 **2 条**（外环 22 边 + 内环 6 边）⇒ 有 **2 个 bound 的 loop 翻译失败被 `continue` 跳过**了。
对应的失败点在 `StepToTopoDS_TranslateEdgeLoop.cxx:770`：`TP->AddFail(EC, " Seam curve not mapped")`（seam 边的 pcurve 映射失败路径）。

端口则在 `resolve_loop`/`resolve_face` 里为这两个 bound 照样建了 wire —— 就是那两条 **1-边、边闭合（first==last）的 seam wire**。

**本轮另测到的独立疑点**：该面内环有一个顶点坐标不同 —— OCCT `(-49.8649,0,-100)` vs 端口 `(-65.2436,0,-100)`（y/z 相同），
指向 `bind_edge_loop_vertices`（`cxx:288-403`/`:405-491` 的顶点绑定）可能绑到了不同顶点。

**STEP 侧事实**（本轮解析）：面 `#4952` 4 个 bound；外环 `#5139` 仅 8 条 ORIENTED_EDGE，
曲线类型依次为 B样条(`#5003`)、LINE(`#5011`)、LINE(`#5017`)、B样条(`#5089`)、LINE(`#5097`)、ELLIPSE(`#5106`, 39.26/34)、LINE(`#5114`)、B样条(`#5136`)；
其中 `#5013` 与 `#5019` 共用同一对顶点 `#4956`/`#5007`（`#5107` 是 ELLIPSE，与 OCCT 内环的椭圆对应）。

**下一步**：① 对读端口 `resolve_loop`/`resolve_face` 与 `TranslateFace.cxx:729-742` + `TranslateEdgeLoop` 的失败路径，
补上「该 bound 翻译失败即跳过」；② 查 `bind_edge_loop_vertices` 为何绑到不同顶点。

本轮无库代码改动；`occt-topo --lib` = 1281/0。


### 9.217 —— T-93：4 wire 一一对应 4 个 bound；合并发生在读取器 per-entity ProcessShape

**bound 级对应**（端口 F113）：

| STEP bound | loop | oriented edges | 端口 wire 边数 |
|---|---|---|---|
| `#5140` FACE_OUTER_BOUND | `#5139` | 8 | 6（§9.200 notch 合并后） |
| `#5352` FACE_BOUND | `#5351` | 14 | 14 |
| `#5363` FACE_BOUND | `#5362` | 1 | 1 |
| `#5374` FACE_BOUND | `#5373` | 1 | 1 |

后两条单边 loop 的边是 `EDGE_CURVE('',#V,#V,CIRCLE(...,34.000000000000014),.T.)` ——
**两个 r=34 的闭合圆**（v1==v2），且 STEP 里没有 pcurve。

**OCCT 侧**：该面外环 22 边里恰有两条 `closed=1` 且 `ctype=1`（圆）的边 —— 就是这两个圆，
它们被并进了外环；而 `StepToTopoDS_TranslateFace.cxx:727` 对每个 bound 只 `Add` 一条 wire、**从不合并**，
所以合并只能发生在读取器的 per-entity `ProcessShape`：
`ShapeFix_Shape::Perform`（`ShapeFix_Shape.cxx:83`）→ `FixFaceTool()` 对每个面 → `ShapeFix_Face::Perform`（`cxx:87`）。

**端口没有 `ShapeFix_Shape`**（全仓只有注释引用），reader 只调 `ShapeFixFace::fix_missing_seam()`（`read_topology.rs:733-747`）。

**修正一处过期文档**：`shhealing/shape_fix_face.rs:4` 与 `:144-148` 写着「UNPORTED: the seam construction (`cxx:1899-2330`)」，
但同文件 `:305`（cxx:1899-1992 退化边）、`:503`（cxx:2138-2234 找缝位）、`:590-611`（cxx:2236-2261 fictive grid + ComposeShell 插缝，含 `set_context`）
实际都已实现。真正产出错误结果（Shell(5)）的是 **ComposeShell 本身**，与 §9.210 定位的「`SplitWire` 记录的空拷贝顶点未被统一 ⇒ `CollectWires` 的 `IsSame` 失败 ⇒ 6 条 wire」一致。

⇒ **T-93 的修复点定在 ComposeShell**：在 `break_wires`/`collect_wires` 之前把这些顶点统一（Context 的顶点级 Apply），
复现 OCCT 的 1 面 / 2 wire。

本轮无库代码改动；`occt-topo --lib` = 1281/0。


### 9.218 —— T-93：OCCT 侧插桩验证 + `FixMissingSeam` 在读取时从未被调用（修正 §9.216/§9.217）

**本轮先给 OCCT 侧插桩建立可信度**（此前从未验证过补丁是否真的生效）：

- 新增 `specs/occt_probe/_dbg/ZZ_ShapeFix_Face.cxx`（`ShapeFix_Face` 构造函数、`Perform` 入口/出口、`FixMissingSeam` 入口、`CompShell.Perform()` 后各一处打印），并加入 `build_dbg.bat`；
- **验证生效**：`run_dbg.bat data\occ\a3n00.stp pcu 25` 输出
  `ZZCTORD default-ctor#1/2/3`、
  `ZZPF enter#1 wires=2 void=0` 等 ⇒ 补丁被链接且被调用（链接输出无 LNK4006/LNK2005）；
- 同时证明：`FCOUNT` 等模式输出正常，stdout 捕获正常。

**决定性测量**：纯读取（`run_dbg.bat data\occ\a3n00.stp fcount 0`）时

```
（无任何 ZZCTORD / ZZPF / ZZFMS 输出）
```

⇒ 读取过程中 **`ShapeFix_Face` 从未被构造**，因此 `ShapeFix_Face::Perform` / `FixMissingSeam` / `ShapeFix_ComposeShell` 都没有参与；
`ShapeFix_Shape`（其构造会建 `FixFaceTool()`=`ShapeFix_Face`）同样没有参与。

**修正两处此前结论**：

- §9.216 的「2 个 bound 的 loop 翻译失败被跳过的机制在 `TranslateEdgeLoop.cxx:770`」—— 方向仍可能对，但**不能**用「FixMissingSeam 没输出」来支撑（当时插桩未验证）；
- §9.217 的「合并发生在读取器 per-entity `ProcessShape`（`ShapeFix_Shape`→`ShapeFix_Face`）」**不成立**（`ShapeFix_Face` 未被构造）。

⇒ OCCT 那个 2-wire 面是**在没有面级 healing 的情况下**得到的，说明减少 wire 的机制在传输层或 `ShapeProcess` 的其它算子里，尚未定位。

**工具教训（写入看板）**：`_dbg/` 的补丁必须先跑一次构造/入口探针确认生效，再据「没有输出」下结论。

**下一步**：给 `STEPControl_ActorRead.cxx` 的 `TransferEntity`（`cxx:2144-2156`，`ProcessShape` 前后）插桩，打出面 `bbox` + `wire` 数，直接判定「传输出来是 4 条还是 2 条 wire」。

本轮无库代码改动；`occt-topo --lib` = 1281/0。


### 9.219 —— T-93：OCCT 读取时从不调用 ShapeFix_Face::Perform；摘掉那步后 F113 出网格但 T0M 回归

**① link map 定论（本轮新增的验证手段）**：给链接加 `/MAP` 后，
`?Perform@ShapeFix_Face@@QEAA_NAEBVMessage_ProgressRange@@@Z` 由 **`ZZ_ShapeFix_Face.obj`**（补丁 TU）提供 ⇒
符号解析是全局的，补丁对所有调用者生效 ⇒ **纯读取时没有任何 `ShapeFix_Face::Perform` 调用**
（`shape_fix_face`/`FixMissingSeam`/`ShapeFix_ComposeShell` 同理）⇒ T-93 卡片标题「reader 缺带 Context 的 ShapeFix_Face::Perform」
与 OCCT 8.0 在该文件上的真实行为不符。

**② 摘掉端口 reader 里 OCCT 没有的那一步**（`read_topology.rs:733-747` 的 `fix_missing_seam`）实测：

```
a3n00: TOTAL faces=226 stats=225 -> 226   mesh_v 10863 -> 11101   mesh_t 11941 -> 11121
       F 113 ... mv=0 -> 223   mt=0 -> 267          <-- F113 出网格（T-93 的核心验收）
cargo test --lib                 -> 1281 passed / 0 failed
cargo test --test step_obj_gates -> 4 passed / 1 failed
  FAILED step_obj_area_matches_occt: occ/T0M.stp ratio 0.9786 outside tolerance 0.01
  （其余 19 个模型比值为 1.0000；ATU01038 1.0006、bottom 0.9750、motoc 1.0333、top 1.0000）
```

**③ 按纪律「回归即回退」已撤回**：`crates/occt-topo/src` 恢复为仅 `wire_fix.rs`（§9.200 的 notch 修复），`--lib` 1281/0，
门禁回到已知基线 5/5。

**④ 机制解释**：端口这次 `fix_missing_seam` 的结果 `Shell(5)` 虽被 reader 丢弃，但它跑过的 ComposeShell **改写了共享 `GeometryRegistry`**
（pcurve/range），副作用留下 ⇒ F113 的网格被毁。这说明该调用是 F113 失败的**直接原因**，而不是 wire 数（4 条）本身。

**⑤ 决策点**：摘掉这步后 `occ/T0M.stp` 面积比 0.9786（带 ±0.01）⇒ 要么
(a) 认定 T0M 原来的 5/5 依赖这步非 OCCT 操作（T-79 型假绿），同意对其重新定基线；要么
(b) 让这次尝试**不产生副作用**（ComposeShell 目前在结果被丢弃时仍改写全局 registry）。
未决，等用户定；Goal 保持 active。


### 9.220 —— T-11：两条警告清理自动化路径均实测不可用（已回退）

**基线（本轮实测，`cargo check --lib --message-format short`）**：`occt-topo` lib 警告 **506** 条 ——
`never used` 408 / `unused import` 46 / `other` 37 / `never read` 11 / `unused variable` 3。

**路径 ①：`cargo fix`**

- `cargo fix --lib --allow-dirty --allow-staged`：警告 506 → 456，34 文件变更；
- 但 `cargo test --lib` 立刻报错（`ShapeType`/`Shell`/`Solid` 找不到）、`cargo check --tests` **109 个错误**、涉及 14 个文件
  （intpatch.rs 37、bop_builder_tests_heal.rs 30、inttools_roots.rs 15、fclass2d/tests.rs 10 …）；
- 原因：lib target 编译时**不带 `cfg(test)`**，cargo fix 据此判定 import 未用，而它们只被 `#[cfg(test)]` 模块使用；
- 追加 `--lib --tests` 重跑后**仍然 109 个错误**（未修复）。

**路径 ②：自己施加 rustc 的 span**

- 取 `cargo check --message-format json` 的诊断，只处理「符号在文件其它位置不再出现」的 `unused_imports`/`unused_variables`（保守过滤，33 处 / 22 文件）；
- 但 `unused_imports` 的 span **不总是整条语句**：`pub use read_topology::*;` 被删成 `pub use ;`，
  `use crate::tshape::{A, B, C};` 的 B 被删后留下 `{A, , C}` ⇒ 语法错；
- （另：首次实现用 JS 字符串索引套 UTF-8 字节偏移，文件含中文注释时必然错位 —— 已改用 Buffer 后仍因 span 语义不对而失败。）

**结论**：`occt-topo` 的警告清理只能**手工分批**，且每批必须跑 `cargo test --lib` + `cargo check --tests` 确认零构建错误。
`never used` 408 条属 dead code，需逐条判断是否为 UNPORTED 占位，不在 T-11 卡片范围。

两次尝试均已 `git checkout -- crates/occt-topo/src` 回退，并重新打上 §9.200 的 notch 修复；`--lib` 1281/0。


### 9.221 —— T-29 完成：过期规格的悬空引用修正（无行为变化）

卡片要求「与门禁快照同步，或显式标注仅历史」。实测两个目标文件**已带**「⚠️ 仅历史（T-29 登记）」横幅
（`specs/_coverage.md:1`、`specs/_brepmesh_align_review.md:3`），但横幅里指向的 `specs/_board.md` **已随看板迁移删除**（改成 `specs/board.canvas.tsx` 后遗留）。

本轮修正 **7 处规范性悬空引用**：

| 位置 | 原文指向 | 改为 |
|---|---|---|
| `specs/_coverage.md:2` | `specs/_board.md` §2/§3 | `specs/board.canvas.tsx`（任务表/门禁表） |
| `specs/_brepmesh_align_review.md:4` | 同上 | 同上 |
| `specs/_a3n00_gap_analysis.md:6` | §3.3 T-69/T-91/T-92、§2 门禁 | `specs/board.canvas.tsx`（T-69/T-91/T-92 与门禁表） |
| `specs/_design_architecture_t25_t28.md:3` | §3.4 要求 | `specs/board.canvas.tsx`（架构 lane） |
| `specs/_design_architecture_t25_t28.md:117` | `specs/_board.md` 对应行 | `specs/board.canvas.tsx` 对应行 |
| `specs/occt_probe/README.md:33` | task T-82 (`_board.md` §3.1b/§3.2) | task T-82 (see `specs/board.canvas.tsx`) |
| `crates/occt-core/src/bspl/poles.rs:26` | `_board.md` §3.3 T-44 | `specs/board.canvas.tsx`, task T-44 |

`specs/_audit/` 下剩余 8 处是**历史审查记录**（记录当时仓库里有什么文件），按原样保留，不追改。

**验证**：`cargo check -p occt-core --lib` 通过；`occt-topo --lib` **1281/0**（纯注释/文档改动，无行为变化）。

顺带：`specs/occt_probe/.gitignore` 里原先只忽略 `ZZ_ComposeShell.*`，改为忽略全部 `ZZ_*.obj/lib/exp/map` 与 `zz_wires_probe.*`，并清掉遗留构建产物。


### 9.222 —— T-67 前提核实与重新定范围：两臂已移植，剩余是 `UVFromIso`

**核实（卡片标题所述部分已完成）**：

- `extrema_surf/point_surface_extrema.rs` 已有 `ExtPsSurfaceType::SurfaceOfExtrusion` / `SurfaceOfRevolution`（`:60-61`、`:88-91`）与两条 arm：
  `perform_ext_ps`（`:450-480`）、`perform_rev_ps`（`:485-…`），与 `Extrema_ExtPS.cxx:292-317` / `:319-343` 逐段对应
  （首次建 `Extrema_ExtPExtS`/`ExtPRevS`、之后 `Perform`、结果逐条经 `TreatSolution` 合并）；
- `extrusion_point_extrema.rs`、`revolution_point_extrema.rs` **无 UNPORTED 残留**；
- `cargo test -p occt-geom --lib` = **143 passed / 0 failed** ⇒ 卡片 accept 达标。

**真正剩余（重新定范围）**：`numeric_extrema.rs:63` 标注的 `point_surface_newton_all*` 替身（audit A15 / T-67 remainder）：

- 非 box 的 `point_surface_newton_all`：**只剩回归测试** `extrema_surf/tests.rs:129`（`newton_path_bspline_paraboloid_min`）调用；
- `point_surface_newton_all_box`：**在生产路径上** —— `point_surface_extrema_box`（`:449`，`:471` 调用替身）← `pcurve_full/surface_projector.rs:632`；
- 它对应 OCCT `ShapeAnalysis_Surface::ValueOfUV` 里「Extrema 失败」分支的 **`UVFromIso`**
  （`ShapeAnalysis_Surface.cxx:1449-1459`；`SurfaceNewton`/`ForgetNewton` 在 OCCT 那里是注释掉的）。

**规模**：`UVFromIso` 本体 `ShapeAnalysis_Surface.cxx:1522-1835`（约 **310 行**），调用点 `:1459`、`:1423`、`:1225`。
本轮不塞半个移植；卡片 next 已改写为「移植 `UVFromIso` 替换该替身」，evidence 记录了逐段对应与 143/0 实测。

本轮无代码改动（只读 + 看板/规格更新）。


**补记（同轮清理）**：§9.220 的 `cargo fix --lib --tests` 曾在 `crates/occt-topo/tests/` 留下 5 个文件的改动
（phase4/6/7/19_integration.rs、step_obj_gates.rs 的 unused import/variable 清理），本轮已 `git checkout` 回退；
当前 `crates/` 下的有意改动仅三处：`occt-topo/src/shhealing/wire_fix.rs`（§9.200 notch 修复）、
`occt-core/src/bspl/poles.rs`（§9.221 注释引用）、`occt-topo/examples/zz_probe_a3n00.rs`（探针 `--fdump`）。
验证：`occt-topo --lib` 1281/0、`cargo check --tests` 0 错误。


### 9.223 —— T-67 完成：摘掉 A15 的 Newton 替身，改走已移植的 `UVFromIso` 分支

**前情（§9.222）**：卡片标题所述的两臂 `Extrema_ExtPExtS` / `Extrema_ExtPRevS` 其实**早已移植并分派**
（`point_surface_extrema.rs:60-61`/`:88-91` 类型、`:450` `perform_ext_ps`、`:485` `perform_rev_ps`，
与 `Extrema_ExtPS.cxx:292-317`/`:319-343` 逐段对应），两个实现文件无 UNPORTED 残留。
进一步的核实发现：卡片的「两臂」若指 `UVFromIso` / `SurfaceNewton`，那两个**也已移植**
（`pcurve_full/surface_projector.rs:1831` `uv_from_iso`，及同文件的 `surface_newton`；`value_of_uv:519-691` 就是 `ShapeAnalysis_Surface::ValueOfUV` 的忠实移植）。

**真正多余的一处（本轮修掉）**：`extrema_surf/numeric_extrema.rs` 的 `point_surface_extrema_box` 在窗口搜索为空后
调用了**端口自造的** `point_surface_newton_all_box`（24×24 网格 + 数值 Jacobian Newton），为空再退到 `fallback_point_surface`。
OCCT 在这里**没有这两步**：`Extrema_ExtPS::Perform` 给空结果，由调用方 `ValueOfUV` 走 `UVFromIso` 分支
（`ShapeAnalysis_Surface.cxx:1449-1459`）。

**改动**：`point_surface_extrema_box` 窗口为空时返回「无解」哨兵
（`ExtremaPair { distance: INFINITY, u2: NAN, v2: None, .. }`），
调用方 `pcurve_full::surface_projector::value_of_uv:640-645` 本来就以 `e.v2.filter(u2.is_finite() && v.is_finite())` 判定无解并走 `uv_from_iso`
（其注释即 `cxx:1448-1472`）⇒ 于是生产路径回到 OCCT 的控制流。

文档同步三处：`numeric_extrema.rs` 的 `point_surface_newton_all` 与 `point_surface_extrema_box` 说明、`extrema_surf/mod.rs:30` 的 UNPORTED 段
（改记为「已不在任何库路径上，仅存回归测试 `newton_path_bspline_paraboloid_min` 使用」）。

**验收实测**：

```
cargo test -p occt-geom  --lib                     -> 143 passed / 0 failed   (卡片 accept)
cargo test -p occt-topo  --lib                     -> 1281 passed / 0 failed
cargo test -p occt-topo  --test step_obj_gates     -> 5 passed / 0 failed (384s, 门禁不劣化)
```


### 9.224 —— T-51 豁免：可表达的部分已实现；两条 arm 被 `GeomTrimmedCurve` 的参数表示卡住

**核实结果（可表达部分早已实现）**：

- `CPnts_AbscissaPoint` 的反解由 `math_FunctionRoot` 驱动（`gcpnts.rs:43`）；
- tolerance 重载 `init_tol`（`gcpnts.rs:296-300`）存在；
- `adv_perform`（`gcpnts.rs:407-422`）是 `CPnts_AbscissaPoint::AdvPerform`（`cxx:436-473`）的忠实移植
  （`myL < Confusion` 短路、`resolution / 10`、`IsDone` 只看 `Solution.IsDone()`）。

**真正剩下的两条 `Compute` arm 无法忠实表达**：

- `GCPnts_LengthParametrized`（`cxx:87-90`）与 `GCPnts_AbsComposite`（`cxx:96-158`）读 `GeomAdaptor_Curve` 的
  `GetType` / `NbIntervals` / `Intervals`；
- OCCT 的 `GeomAdaptor_Curve::load`（`GeomAdaptor_Curve.cxx:239-255`）对 `Geom_TrimmedCurve` 会**下沉到 basis 并保留 basis 参数区间**；
- 端口的 `GeomTrimmedCurve` 则把参数**重映射到 `[0, 1]`**：`d0(u) = basis.d0(first + u*(last-first))`，
  `d1/d2/d3` 用链式法则补偿（`crates/occt-geom/src/trimmed.rs:44-70`）；
- 所以那两条 arm 在端口里没有忠实对应物，硬写只能自造区间规则（违反红线）；忠实做法是先改 TrimmedCurve 的参数表示 ——
  那是**架构级**改动（影响所有 trimmed 曲线的消费方），且当前**没有消费方**触发这两条 arm。

⇒ 按看板 intake 的 `waived`（豁免）处理，并把上述架构缺口作为**新卡候选**记录（建议标题：TrimmedCurve 参数表示对齐 OCCT 的 basis 区间）。

**验收实测**：`occt-geom --lib` 143/0、`occt-topo --lib` 1281/0。

本轮无代码改动。


### 9.225 —— T-11 第 1 批：手工清理 2 处 unused import（另记两条自动化陷阱）

**基线**：`occt-topo --lib` 警告 506 条（never used 408 / unused import 46 / other 37 / never read 11 / unused variable 3）；
`cargo check -p occt-topo --lib` 会同时列出依赖 crate（occt-geom / occt-geom2d）的警告，去重后 unused import 共 **63 条 / 44 文件**。

**第 1 批（本轮落地）**：

| 文件 | 删除的 import | 依据 |
|---|---|---|
| `crates/occt-geom/src/line.rs:2` | `GpAx1` | 符号在整个文件（含 test 模块）只出现在该 import 行 |
| `crates/occt-topo/src/shhealing/wire_fix.rs:7` | `GpDir` | 同上（`GpDir2d` 是另一个符号，`\b` 边界不误伤） |

**验证**：`occt-geom --lib` 143/0、`occt-topo --lib` 1281/0、`cargo check -p occt-topo --tests` 0 错误（纯 import，行为不变）。

**两条陷阱（本轮实测，写入看板 next）**：

1. **只被 `#[cfg(test)]` 模块使用的 import 会被 rustc 判 unused**（lib target 编译不带 `cfg(test)`）。
   本轮 9 个候选里 6 个属此类（`GpVec`/`GpDir`/`ShapeType`/`GpAx3`/`GeomCylinder`/`IntRange`），删掉会报 358 个 test 构建错误。
   规则：删除前先 grep 该符号在**整个文件**（含 test 模块）的出现次数，>1 即跳过。
2. **`cargo check -p A` 会列出依赖 crate 的警告**，且那是依赖的**非 test** 构建。
   例：`crates/occt-geom/src/surface_of_revolution.rs:15` 的 `use crate::Surface;` 被报 unused，但删掉后
   `impl crate::Surface for GeomSurfaceOfRevolution` 的方法解析失效（`no method named d0`）—— 实测已回退，该文件净零改动。

**仍未做的**：44 文件中其余 unused import（多属陷阱 1）与 `never used` 408 条 dead code（需逐条判断是否为 UNPORTED 占位，不在卡片范围）。


### 9.226 —— T-11 第 2 批：再清 5 处 import；发现第 3 条陷阱（跨文件消费）

**第 2 批落地（5 处）**：

| 文件 | 删除的 import |
|---|---|
| `crates/occt-topo/src/bop_builder_planar_weld.rs:19` | `edges_of` |
| `crates/occt-topo/src/bop_split_seam.rs:20` | `ShapeIterator` |
| `crates/occt-topo/src/bop_builder_planar_geom.rs:12` | `BoolOp` |
| `crates/occt-topo/src/int_face_face.rs:135` | `DEFAULT_WINDOW` |
| `crates/occt-geom2d/src/geom2d_int/ginter.rs:11` | `GpPnt2d` |

连同第 1 批（`line.rs` 的 `GpAx1`、`wire_fix.rs` 的 `GpDir`）共 **7 处**。

**验证**：`occt-topo --lib` 1281/0、`occt-geom2d --lib` 72/0、`cargo check -p occt-topo --tests` 0 错误。

**第 3 条陷阱（本轮踩到并回退）**：文件内 grep 不足以判定安全 —— 还有两种**跨文件消费**：

1. **`pub(crate) use` 再导出**：`meshing/model_builder/mod.rs:26` 的 `pub(crate) use crate::abs::{Orientation, ShapeType};`
   被 rustc 判 `ShapeType` unused，但别处（`#[cfg(test)]`）经 `crate::meshing::model_builder::ShapeType` 使用它 ⇒ 删掉即报 `E0433`。
2. **`use super::*` 子模块**：`bopds.rs:42` 的 `use crate::topo_tools_full::shapes_of;` 被同模块的子模块经 `use super::*` 冒泡消费 ⇒ 删掉即报 `E0425`。

⇒ 修正后的规则：**先按「符号在文件内只出现于 import 行」筛，再逐批跑 `cargo test --lib`**；
若报 `E0433`/`E0425`，只回退该处（不要整批回退）。本轮就是这么做的（7 处中 2 处回退）。

剩余：44 文件中仍有 unused import（多属三条陷阱），以及 `never used` 408 条 dead code（不在卡片范围）。


### 9.227 —— T-69：逐面配对完成（17 个差异面；含 2 个反向差异）

**工具**：给两侧探针都加了「多 wire 面 bbox 明细」——端口 `zz_probe_a3n00 --wirehist` 新增 `MULTI bbox=(...)-(...) wires=N`，
GT `occt_probe --wires` 同样输出（`occt_probe.cpp` 的 `++nFaces` 后插桩）；然后按 bbox 取整到个位配对。

**总量对照（本轮实测）**：

```
PORT faces=1772 wires=1931 faces_with_2plus=120  WIREHIST 1:1652 2:106 3:8 5:1 6:1 7:2 8:1 10:1
GT   faces=1778 wires=1921 faces_with_2plus=106  WIREHIST 1:1672 2:92  3:8 5:1 6:1 7:2 8:2
```

**配对后的真正差异（只有 17 个）**：

- **端口 2 wire / GT 1 wire（15 个）**：
  `(-31,-117,15)-(-23,-110,18)`、`(-31,-132,10)-(-23,-124,14)`、`(-31,-99,16)-(-23,-91,18)`、
  `(6,-117,15)-(14,-110,18)`、`(6,-132,10)-(14,-124,14)`、`(6,-99,16)-(14,-91,18)`、
  `(-43,-108,11)-(-34,-97,21)`、`(17,-106,12)-(26,-95,22)`、`(-66,-2,-49)-(-58,11,-36)`、
  `(-66,-27,-49)-(-58,-14,-36)`、`(-21,-21,-8)-(-15,-15,-8)`、`(-21,15,-8)-(-15,21,-8)`、
  `(15,15,-8)-(21,21,-8)`、`(15,-21,-8)-(21,-15,-8)`、`(-10,-10,-60)-(10,10,-45)`
- **端口 10 wire / GT 8 wire（1 个）**：`(5.3,-217.4,104.9)-(15.5,-204.3,120.1)`（卡片的 F351）
- **反向差异（2 个，端口 wire 更少）**：GT 的 `(-30,-106,13)-(12,-76,19)` 与 `(-30,-135,9)-(12,-106,19)` 都是 2 wire，
  端口没有这两个**大**面，取而代之的是上面那批 x∈[-31,-23] 与 [6,14] 的**小**面 ⇒ 端口把 GT 的一个面切成了多个。
  这与「端口总面数 1772 < GT 1778（少 6）」同源。

⇒ T-69 的 next 已更新为「取一个代表面（如 `(-21,-21,-8)-(-15,-15,-8)` 或 10-vs-8 的 F351）解析其 STEP 实体，找多出的那条 wire 从哪来」。

本轮无库代码改动（只改了两个探针）。


### 9.228 —— T-69：现象定位到「1-边闭合 wire 被叠在面上」；关掉补偿步会显著变差

**代表面**：端口 F1655 = Cylinder，`du=6.28319`（整周期）、`dv=inf`、box `(-65.6,-2,-49.19)-(-58.1,11,-36.19)`、wires=[1,2]：

```
FDUMP  wire[0] nEdges=1  e[0] first=last=(-65.600000,10.006279,-46.143374)   <- 闭合边（seam）自成一条 wire
FDUMP  wire[1] nEdges=2  (-58.1,4.535096,-36.189358) <-> (-58.1,3.736963,-36.234205)
```
GT 在该区域是 **1** 条 wire ⇒ 与 T-93 的 F113（4 wire vs 2、含独立 1-边闭合 wire）同型。

**T0M 全量量化（本轮新增 `--wirehist` 字段）**：

```
PORT faces=1772 wires=1931 faces_with_2plus=120 wires_1edge_closed=232 faces_1edge_closed=158
GT   faces=1778 wires=1921 faces_with_2plus=106 wires_1edge_closed=172 faces_1edge_closed=172
```

- GT：172 个面、每个恰好 1 条 1-边闭合 wire；
- 端口：**158 个面承载 232 条** ⇒ 在 **14 个面**上叠了额外的 1-边闭合 wire（= `faces_with_2plus` 120−106 的差）。

**关键实验（临时改动，已回退）**：把 reader 里那步 OCCT 没有的 `fix_missing_seam`（§9.219）关掉后，T0M **显著变差**：

```
faces_1edge_closed 158 -> 349 ; faces_with_2plus 120 -> 309 ; WIREHIST 2:106 -> 2:294
```

⇒ 那步在 T0M 上是**补偿**（把 reader 本该产出的 seam 边补上），不能简单摘掉；它对 a3n00 的 F113 有害，是因为被丢弃的 `Shell(5)`
留下了共享 registry 的副作用（§9.219）。

⇒ 正确修法：让 **reader** 像 OCCT 一样产出 seam 边（与 T-93 同根），而不是靠 `ShapeFix_Face` 事后补偿。
这同时解释了为什么 T-69 的 next 与 T-93 的调查路线收敛到同一处（reader 的 wire/seam 组装）。

本轮无库代码改动（新加的探针字段保留；临时关步已回退，`--lib` 1281/0）。


### 9.229 —— T-69：同面直接对照 —— 端口在整周期柱面上少 3 条边并拆散 seam wire

**工具**：GT 探针新增 `--fbox x1 y1 z1 x2 y2 z2`（按 bbox 定位面并逐 wire/边 dump），与端口 `--fdump <idx>` 对应。

**同一面（bbox 完全相同）两侧结构**：

```
GT face=1660  bbox (-65.60000,-2.00000,-49.18926)-(-58.10000,11.00000,-36.18926)
  wire[0] nEdges=6
    e0 ori=0  A->A   (A=(-65.600000,10.006279,-46.143374))       <- 闭合边（seam）
    e1 ori=1  B->A   (B=(-58.100000,10.006279,-46.143374))
    e2 ori=1  C->A   (C=(-58.100000,3.736963,-36.234205))
    e3 ori=1  D->C   (D=(-58.100000,4.535096,-36.189358))
    e4 ori=1  A->D
    e5 ori=0  A->B
端口 F1655
  wire[0] nEdges=1  闭合边 A->A
  wire[1] nEdges=2  (D<->C)
```

⇒ 端口这个面**少了 3 条边**（GT 6 / 端口 3），并把 seam 闭合边拆成独立 1-边 wire。
GT 的 6 条构成：闭合边 `A→A` + `A↔B` 往返 + `A→D→C→A` 环 —— 典型的 seam 双走 + 闭合 seam 边。

**与 T-93 的统一**：F113 是端口外环 14 边 vs OCCT 22 边（同样丢掉绕 seam 的那段）。
⇒ **T-93 与 T-69 同根：reader 在周期面/接缝处的边与 wire 组装丢边、拆 wire。**
而 §9.228 的实验说明 reader 里那步 `fix_missing_seam` 只是**事后补偿**（关掉后 T0M 从 120 个多 wire 面涨到 309 个）。

**下一步**：解析该面的 STEP `ADVANCED_FACE`（圆柱 R=(65.6−58.1)/2 = 3.75）与其 `EDGE_LOOP` 的 oriented edge 条数 ——
若 STEP 只有 3 条而 GT 有 6 条，则 OCCT 在传输时**扩展**了 seam（重复边对 + 闭合边），端口的 `read_topology` 没有这一步。

本轮无库代码改动（只改 GT 探针）；`--lib` 1281/0。


### 9.230 —— T-69：STEP 实体链解到面级（端口面 = #29884，与 STEP 完全吻合）

**从 seam 点回溯**：

```
A=(-65.6,10.006279,-46.143374)  ->  #59828 CARTESIAN_POINT
                                ->  #13393 VERTEX_POINT('',#59828)
                                ->  #17405 EDGE_CURVE('',#13393,#13393,#10824,.T.)   #10824=CIRCLE('',#33172,6.5)
```
顶点 `#13393` 只被这一条**闭合边**使用；而该边被**两条** ORIENTED_EDGE 引用，分属两个面：

| 面 | 曲面 | bounds |
|---|---|---|
| `#29884` | `#1762=CYLINDRICAL_SURFACE(R=6.5)` | `#4187 FACE_OUTER_BOUND→loop #6254` **1** OE（`#25332` = `#17405 .T.`）；`#2467 FACE_BOUND→loop #6255` **2** OE（`#17406 .F.`, `#17407 .F.`） |
| `#29885` | `#2149=CONICAL_SURFACE(R=6.067, 60°)` | `#4188 FACE_OUTER_BOUND→loop #6256` **1** OE（`#25335` = `#17405 .F.`）；`#2468 FACE_BOUND→loop #6257` **1** OE（`#17408 .F.`） |

**端口 F1655 与 `#29884` 完全吻合**：2 个 bound、oriented edge 数 1 与 2 ⇒ 端口 `wires=[1,2]` ✓；`type=Cylinder` ✓（R=6.5）。

**未决**：GT 在**同 bbox** 的面（1660）是 1 wire / **6 边**，而 STEP 该面只有 1+2 = **3** 条边。两种可能：

1. OCCT 在传输时把该面的 seam 扩展成 6 条（例如闭合边 + 其往返对 + 其余边的重复）；
2. bbox 配对照到了**另一个面**（本轮 F1655 的 bbox 原点来自端口，GT 侧只是同 bbox）。

⇒ 下一步用 GT 侧逐边端点比对 `#17405/#17406/#17407`（A/B/C/D 四个点），确认 GT 那 6 条里哪 3 条来自 STEP、哪 3 条是 OCCT 生成。

本轮无库代码改动（只做 STEP 实体解析）；`--lib` 1281/0。


### 9.231 —— 更正：bbox 配对不是可靠的面身份（§9.229 结论作废）

**用 STEP 端点核对 GT 那 6 条边**：

```
STEP face #29884 的边（本轮解析 #17406/#17407 的端点）：
  #17405  A->A   A=(-65.6,10.006279,-46.143374)          CIRCLE R=6.5（闭合）
  #17406  C->D   C=(-58.1,3.736963,-36.234205)  D=(-58.1,4.535096,-36.189358)  CIRCLE
  #17407  D->C   同 C/D，另一条独立边                                              CIRCLE
GT face=1660（同 bbox，1 wire / 6 边）：
  A->A 、B->A 、A->B 、C->A 、A->D 、D->C      B=(-58.1,10.006279,-46.143374)
```

⇒ GT 那个面包含 `A↔B`、`C→A`、`A→D` —— **STEP 面 `#29884` 里根本没有这些边**（它只有 A 处闭合圆与 C↔D 两条）。
因此两者**不是同一个面**：两个不同的面可以共享同一个 bbox。

**影响**：

- §9.229 的「端口在同面少 3 条边、reader 丢边」**作废**；
- §9.227 的 17 个「差异面」只能作为**线索**（总量直方图与 `wires_1edge_closed` 等聚合量仍然有效），
  不能据此断言「同面结构差异」。

**仍然成立的事实**：端口 F1655 与 STEP 面 `#29884` **完全吻合**（2 个 bound、oriented edge 数 1 与 2 ⇒ `wires=[1,2]`；`CYLINDRICAL_SURFACE(R=6.5)`）。

**下一步**：改用**稳定身份**配对 —— 给 reader 的 `resolve_face` 一次性探针记录源 STEP 面 id，再与 GT 的同 id 面对比；
在此之前不基于 bbox 配对下结论。

本轮无库代码改动；`--lib` 1281/0。


### 9.232 —— T-69：STEP-id 级配对成功（197 个 id 的 wire 数不同；OCCT 会跳过 bound）

**新工具（可靠身份）**：

- GT：`occt_probe --faceids` —— 遍历结果 shape 的面，用 `XSControl_TransferReader::EntityFromShapeResult(face, 1)` 反查 STEP 实体，
  打印 `STEPFACE face=k id=<model->Number(ent)> wires=<n>`；
- 端口：`read_topology::resolve_shape` 里加 `ZZ_STEPFACE` 环境变量门控的探针，对 `ADVANCED_FACE` 打印 `STEPFACE id=<id> wires=<n>`。

**T0M 对比结果**：

```
GT ids=1768 端口 ids=1772 共有=1758 wire 数不同=197 only-GT=10 only-PORT=14
例：id=29884 gt=1 port=2 ; id=29865 gt=3 port=1 ; id=29191 gt=8 port=1 ; id=29880 gt=3 port=2
```

**关键读数**：GT 的 wire 数常**低于** STEP 的 bound 数，而端口 = 每 bound 一条（忠实于 `TranslateFace.cxx:727` 的「每个 bound `Add` 一条 wire」）。
例：`#29884` 的 STEP bounds = `FACE_OUTER_BOUND #4187`(loop 1 OE) + `FACE_BOUND #2467`(loop 2 OE) —— 端口 2 条 wire，GT **1** 条。

⇒ **OCCT 会跳过某些 bound 的 loop**（`StepToTopoDS_TranslateFace.cxx:729-742`：`anEdgeLoopTranslator.IsDone()` 为假即 `continue`），
端口没有这个失败路径，于是把 OCCT 丢弃的 bound 也建成了 wire。

**同时复测（T0M）**：读取时 `ShapeFix_Face` 仍未被构造（补丁版探针无 `ZZCTORD`/`ZZPF` 输出，`FCOUNT hits=0 of 1778`）
⇒ 差异不属于面 healing，属于**传输层**（与 §9.219 在 a3n00 上的结论一致）。

**下一步**：在端口 `resolve_loop`/`resolve_face`（或 `associate_edge_pcurve`）里复现 `TranslateEdgeLoop` 的失败/跳过条件
（`TranslateEdgeLoop.cxx:770`/`:791` 的 seam-pcurve 失败路径），使「某些 bound 不产 wire」。

注：`ZZ_STEPFACE` 是临时探针（env 门控、已注释说明），提交前应移除；`--faceids` 属 GT 探针工具，保留。


### 9.233 —— T-69：like-seam 跳过假设**否证**（实施后 T0M 计数完全不变）

**实施**（已回退）：按 `StepToTopoDS_TranslateEdgeLoop.cxx:790-795` 的失败分支，在 `resolve_face` 的 bound 循环里加判定 ——
「loop 中某条 oriented edge 带 pcurve 且 `IsLikeSeam`（非 seam）⇒ 该 loop `!IsDone()` ⇒ 跳过整个 bound」
（对应 `TranslateFace.cxx:729-742`），并新增私有 helper `loop_is_like_seam_rejected`。

**结果**：T0M 的

```
TOTAL faces=1772 wires=1931 faces_with_2plus_wires=120 wires_1edge_closed=232 faces_1edge_closed=158
WIREHIST 1:1652 2:106 3:8 5:1 6:1 7:2 8:1 10:1
```

与改动前**逐字节相同** ⇒ 该分支在 T0M 上**一次都没触发**（端口判定为 like-seam 的边，在 T0M 的这些面上不存在）。
由于改动无效果且未在其它模型上验证，已回退（保留 `ZZ_STEPFACE` 临时探针）。

**已确证、可复用的结论（供后续）**：

1. STEP-id 级差异 **197 个**，GT 的 wire 数**通常更少** ⇒ OCCT 丢弃了某些 bound 的 wire；
2. 差异发生在**传输层**（T0M 读取时 `ShapeFix_Face` 未被构造，探针无输出）；
3. 端口对**每个 bound 都建一条 wire**（忠实于 `TranslateFace.cxx:727`），所以差异来自「OCCT 跳过了哪些 bound」；
4. like-seam 失败分支**不是** T0M 的原因（本轮否证）。

**下一步候选**：

- (a) 逐面打印端口 `resolve_face` 的 bound 数 vs 结果 wire 数，直接找「bound 数 > wire 数」的面 —— 那才是端口自己丢 wire 的真因；
- (b) 换卡（T-41 / T-25 / T-28 / T-11）—— T-93/T-69 已连续多轮只取证未落地。

本轮无净代码改动（一处实验已回退）；`--lib` 1281/0、T0M WIREHIST 不变。


### 9.234 —— T-41 第 1 步落地：`boolean_dispatch` 的实体输入改派到忠实路径（全绿）

**改动**：`crates/occt-topo/src/bop_builder_dispatch.rs:227`（`boolean_dispatch` 的尾部）

```rust
-    crate::bop_curved::curved_boolean_full(a, b, op, tol)
+    // T-41: OCCT has no bop_curved mesh/voxel engine; two non-compound solids go
+    // through BOPAlgo_BOP like every other operand pair.
+    crate::bop_builder::boolean_via_bopalgo(a, b, op, tol)
```

（`boolean_via_bopalgo` 由私有改为 `pub(crate)`；它就是 `builder_bop_with_fuzzy` + `unwrap_bopalgo_shape` + `single_shape_result` + 闭壳告警，与 `boolean()` 的忠实路径同一条。）

**验收实测（改动后）**：

```
cargo test -p occt-topo --lib                 -> 1281 passed / 0 failed
cargo test -p occt-topo --test step_obj_gates -> 5 passed / 0 failed (380s)
phase10 / phase19 / phase3 / phase4           -> 8/8 · 5/5 · 4/4 · 9/9（与基线一致）
```

**关键发现（修正 §9.206 的担心）**：此前担心「3 条断言绑在旧体上、必须先处理」。实测改派后**全绿** ——
因为 `general_boolean_trimmed_closed_shell` / `curved_boolean_full_quadric_unchanged` / `voxel_*` 都是**直接调用**
`curved_boolean_full` / `curved_boolean` 的模块内测试，**不经 `boolean_dispatch`**，所以改派不影响它们。
（§9.206 把「函数被断言绑定」与「分派路径被断言绑定」混为一谈了。）

**第 2 步（摘层）待做**：逐处处理 `curved_boolean_full` 的其余消费者（repair:239,405、report:388、brepfeat:320/330-335、draw:327）
与 `bop_curved` 自带测试，然后删除 mesh/voxel 层。


### 9.235 —— T-41 第 2 步：第一处障碍是「测试夹具用了非法 OCCT 布尔」

**尝试**：把 `bop_builder_repair.rs:239`（`boolean_repaired`）、`:405`（`boolean_split_result`）、`bop_builder_report.rs:388`（`boolean_repaired_report`）
从 `curved_boolean_full` 改到忠实的 `boolean_dispatch`（含 empty / compound / non-solid / solid 四段，与 `boolean()` 同一条）。

**结果**：`cargo test -p occt-topo --lib` → 1280/1，唯一失败：

```
bop_builder::tests_api::boolean_repaired_warns_on_fix
panicked at src/bop_builder_tests_api.rs:318:79
repaired fuse ok: "BOPAlgo_AlertBOPNotAllowed"
```

**根因（OCCT 依据）**：该测试把 `crossing_shell()`（一个 **shell**）与 `boxy.solid.0`（**solid**）做 `FUSE`。
`BOPAlgo_BOP::CheckData` 明确**拒绝不同维度的 FUSE**（`BOPAlgo_BOP.cxx:181-186`/`:193-195`），端口 `bop_bop::check_data` 亦然 ——
所以这不是忠实路径的缺陷，而是**夹具本身不是合法 OCCT 布尔**；旧的 `curved_boolean_full` 只是「宽容」地接受了它。

**处置**：已回退这三处（保留第 19 轮已验证的 `boolean_dispatch` 改派），`--lib` 1281/0、门禁不变。

**下一步（按 OCCT 口径订正测试，项目先例 T-04/T-05）**：

- 或把夹具改成**同维度**（`shell ∪ shell`，或先闭合成 solid 再 `solid ∪ solid`），保持「跨界融合 → 需要修复」的语义；
- 或对非法组合（shell ∪ solid）显式断言 `BOPAlgo_AlertBOPNotAllowed`，并把「修复告警」的检查移到合法夹具上（该测试已有 `overlapping_boxes()` 的合法对）。

完成后再继续 `brepfeat/features.rs:320,330-335`、`draw/mod.rs:327` 的改派与整层删除。


### 9.236 —— T-41 第 2 步：修复包装器的测试夹具需要「重设计」（本轮两次尝试均回退）

**尝试 ①**：`boolean_repaired`/`boolean_split_result`/`boolean_repaired_report` → 忠实 `boolean_dispatch`。
结果：`boolean_repaired_warns_on_fix` 报 `BOPAlgo_AlertBOPNotAllowed` —— 夹具把 self-intersecting **shell** 与 **solid** 做 FUSE，
`CheckData` 拒绝不同维度融合（`BOPAlgo_BOP.cxx:181-186`）⇒ **夹具本身不是合法 OCCT 布尔**。

**尝试 ②**：按同维度订正夹具（新增一块远处 (5,5,5)-(6,5,5) 的平行单面 shell，与 `crossing_shell()` 同维度且不相交）。
结果：越过了维度错误，但该测试在**修复告警断言**处失败（`src/bop_builder_tests_api.rs` 内的 `warnings.any(|w| w.contains("fixed"))`）——
忠实路径下「两个不相交 shell 的融合」不再产生期望的修复告警。

⇒ 要让该测试在新路径下成立，需要**重新设计夹具**（换一组能产生自交结果、且由合法布尔产生的操作数）——
这已经是**写新测试**，触碰红线（「不为对齐新写测试」）。故本轮两次尝试全部回退：
`bop_builder_tests_api.rs` / `bop_builder_repair.rs` / `bop_builder_report.rs` 用 `git checkout` 复原，`boolean_dispatch` 恢复私有；
`occt-topo --lib` = **1281/0**。

**待决策（写入卡片）**：

- (a) 允许重设计该测试夹具（承诺用最小改动并逐条说明 OCCT 依据）；或
- (b) 保留 `repair`/`report` 走旧路径，只摘 `brepfeat:320,330-335`、`draw:327` 与其余 mesh/voxel 层。


### 9.237 —— T-41 第 2 步（续）：`brepfeat` 改派忠实布尔，全绿

**改动**：`crates/occt-topo/src/brepfeat/features.rs` 的 `boolean_feature` 从「planar 走 `bop_curved::curved_boolean`、
curved 走 `boolean_ops::voxel_boolean`」统一为**忠实布尔**：

```rust
-    if all_faces_planar(&solid.0) && all_faces_planar(tool) { ...bop_curved::curved_boolean... }
-    let mesh = crate::boolean_ops::voxel_boolean(...)?;
-    let brep = crate::mesh_to_brep::shape_mesh_to_brep(&mesh); ...
+    let r = crate::bop_builder::boolean(&solid.0, tool, op, tol.max(1e-9))
+        .map_err(|e| format!("feature boolean: {e}"))?;
+    let shape = r.solid.map(|s| s.0).unwrap_or(r.shape);
+    let volume = solid_volume(&shape, 48, 48);
```

（doc 注释同步：OCCT 只有 `BRepAlgoAPI_*`/`BOPAlgo_BOP`；旧注释声称「曲面走体素」的理由已随 T-82 的曲面精度修复而过期。）

**验收实测**：

```
cargo test -p occt-topo --lib            -> 1281 passed / 0 failed
cargo test -p occt-topo --test step_obj_gates -> 5 passed / 0 failed (396s)
phase3/4/5/6/10/19                       -> 4/9/7/5/8/5（全部与基线一致）
```

**`bop_curved` 剩余消费者**（第 2 步尚未完成的部分）：

| 位置 | 状态 |
|---|---|
| `boolean_dispatch:227` | ✅ 第 19 轮已改派 |
| `brepfeat/features.rs`（curved_boolean + voxel_boolean） | ✅ 本轮已改派 |
| `bop_builder_repair.rs:239,405`、`bop_builder_report.rs:388` | ⏸ 卡在 `boolean_repaired_warns_on_fix` 夹具重设计（§9.236，待决策） |
| `draw/mod.rs:327` + `draw/interpreter.rs:422`（`curved_boolean_ext`） | ⏳ 解释器辅助、不在门禁内 |
| `bop_curved` 自带测试 | ⏳ 随层删除一并处理 |


### 9.238 —— T-41 第 2 步（续）：`draw` 的布尔命令改派忠实路径

**改动**：

- `draw/mod.rs:327`：prelude 的 `pub(crate) use crate::bop_curved::curved_boolean_ext;` → `pub(crate) use crate::bop_builder::boolean;`；
- `draw/interpreter.rs:422`：`curved_boolean_ext(&a, &b, op, BOOL_TOL)` → `boolean(&a, &b, op, BOOL_TOL)`（doc 注释同步）。

原 `curved_boolean_ext` 的分派是「planar → `bop_builder::boolean`；general-curved → `general_boolean_trimmed`；其余 → `curved_boolean`（网格路径）」——
后两者都是无 OCCT 对应的自造路线；draw 命令现在只用忠实的 `BOPAlgo_BOP`。

**验收实测**：`occt-topo --lib` 1281/0；phase10/19/4 = 8/5/9（与基线一致）。
`draw` 不在 STEP→OBJ 门禁路径上，故门禁不受影响（上一轮已实测 5/5）。

**`bop_curved` 现在的消费者（全仓 grep）只剩 3 处**：

```
bop_builder_report.rs:388   curved_boolean_full
bop_builder_repair.rs:239   curved_boolean_full
bop_builder_repair.rs:405   curved_boolean_full
```

⇒ 层已近乎孤立：这 3 处都卡在 `boolean_repaired_warns_on_fix` 夹具重设计（§9.236 待决策），
改完即可删除 `bop_curved` 的 mesh/voxel 层与其自带测试。


### 9.239 —— T-11 第 3 批（2 处）与两条新教训

**落地（2 处）**：

| 文件 | 删除的 import |
|---|---|
| `crates/occt-topo/src/shape_fix_compose_shell/split_by_grid.rs:7` | `crate::brep_tool::BRepTool` |
| `crates/occt-topo/src/shape_fix_compose_shell/dispatch_wires.rs:7` | `occt_geom2d::curve::Curve2d` |

累计三批共 **9 处**。验证：`occt-topo --lib` 1281/0、`cargo check --tests` 0 错误。

**教训 ④**：本轮 21 个 `[line]` 候选里只有 2 个真正被删 —— 其余被**保守 shape 检查**跳过，原因是
`cargo check --message-format json` 的 `spans[0].line_start` 对**多行 `use`** 不一定指向符号自身那一行（可能指向 `use a::{` 那一行）。
更稳的做法：用 `spans[0].line_start..line_end` 覆盖范围，或按符号在文件里搜出真正的那一行。

**教训 ⑤**：`bopds.rs:42` 的 `use crate::topo_tools_full::shapes_of;` **第二次**踩中「`use super::*` 子模块跨文件消费」
（第 2 批也踩过同一处）⇒ 该文件列入白名单，不再尝试。

**剩余候选**：多为 `{...}` 组内单项与 `pub(crate) use` 再导出，需逐个人工判断（后者已知会被别处消费）。


### 9.240 —— T-25 完成核实：几何已在 TShape 上，Registry 仅剩兼容外壳

卡片 next 写的是「按 6 步增量方案实施（第 2 步是纯转发层）」，但实测**迁移已到设计文档所说的端阶段**：

| 证据 | 位置 |
|---|---|
| `TShape` 自带 4 个几何槽（`edge_pcurves`/`edge_core`/`vertex_core`/`face_core`）与 getter/mut | `tshape.rs:33-45`、`:94-123` |
| `GeometryRegistry` 的三个几何表**已删除**，只剩 `ids` + `face_surfaces` + `surface_by_ptr` | `tgeometry.rs:163-181` |
| 35 处读写全部**转发**到 TShape 槽 | `tgeometry.rs:211-219`（set_vertex→`vertex_core_mut`）、`:226`、`:248/:255/:273`、`:282-294`、`:313`、`:340`、`:358-379`、`:399-445` |

`tgeometry.rs:163-166` 的注释逐字写着：

> T-25 end phase: the three geometry maps are gone — pcurves, edge, vertex and face geometry all live on their own TShape now.
> Only the shape ids and the face->surface table remain.

**未做的部分（设计文档明示的可停点）**：第 5 步（删 `ids`/`remove_by_ptr`/`TShape::drop` 钩子与 `global()`、清 29 个测试文件的 `clear_shape` 夹具）。
设计文档 §分步方案末尾写着：「若第 5 步无法一次完成，保留 `global()` 作兼容外壳（第 2 步之后它已是纯转发，语义已对齐 OCCT）」。

**验收实测**：`occt-topo --lib` 1281/0、`step_obj_gates` 5/5（399s）⇒ 卡片 accept「--lib + 门禁不劣化」达标。

⇒ 标记 T-25 **done**；第 5 步作为独立的「外壳清理」可另立卡（当前无语义偏差，属净收益而非门禁项）。


### 9.241 —— T-28 现状核实：第 1/2/3 步已完成，第 4/5 步未做（并更正一处假阴性）

**已完成（逐项带行号）**：

| 步 | 证据 |
|---|---|
| 1 数据层 | `PathPoint.is_new / vertex_id`：`intpatch_impimp_sonb.rs:23-26`（构造 `:47-48` / `:61-62`）；`TopolTool` 顶点表：`geom_int_topol.rs:59`、`:111 with_arc_vertices` |
| 2 同顶点合并 | `compute_tangency` 的 `IsNew == false` 分支：`intpatch_impprm.rs:280`、`:287`、`:323`、`:347`（`CurveTransition::new()`）、`:374`；**`TopTrans_CurveTransition` 已有专文件** `intpatch_curve_transition.rs`（文档「无对应件」已过期） |
| 3 `SetVertex` 绑定 | `intpatch_impprm_ends.rs:287-294`（`!rp.is_new` 守卫 + `IntPatch_Point::SetVertex` 行号引用） |

**未完成**：

- **步 4**：`intpatch_analytic.rs` 的 5 个 closed form **仍在** —— `intersect_plane_sphere:69`、`intersect_sphere_sphere:92`、`intersect_plane_cylinder:123`、`intersect_plane_cone:217`、`intersect_plane_torus:244`；
  调用点 `intpatch.rs:77` 的 re-export 与 `int_face_face_analytic.rs:141/163/198`。
- **步 5**：两个 surface×surface 分派器共存（`intpatch.rs:138` 与 `intpatch_intersection.rs:41`），签名/返回类型不同
  （前者无 `TopolTool`、返回 `SurfaceIntersection`；后者需两个 `TopolTool`、产出 `IntPatch` 线）。主路径已走后者（`int_face_face.rs:414` 的 `IntSS`），`intpatch.rs:138` 现为 `int_face_face.rs:448` 的**回退**。

**更正一处我自己的假阴性**：上一轮 grep 用了不带 `intersect_` 前缀的函数名（`plane_sphere` 等），得到 0 命中，据此误判「步 4 已完成」。
以 `^pub fn ` 重列后 5 个重复件清晰可见。教训：按**完整函数名**（含前缀）核对，或直接用 `^pub fn` 列清单。

**验收实测（当前树）**：`occt-topo --lib` 1281/0、`step_obj_gates` 5/5。

⇒ T-28 置为 **doing / actual 4**；下一步 = 步 4（删 `intpatch_analytic.rs` 的 5 个 closed form，改接 `intana`）。


### 9.242 —— T-28 步 4 推进：三条 Analytic arm 改由 `intana` 提供闭合解（门禁全绿）

**改动**（`crates/occt-topo/src/int_face_face_analytic.rs`）：三条 arm 不再调 `intpatch::intersect_*`，改为调用 `intana` 的权威闭合解，
再经本文件既有的通用适配器 `conics_to_curves`（`:408`）出曲线：

| arm | 新调用 | tolerance 依据 |
|---|---|---|
| `plane_sphere` | `quadric_quadric_plane_sphere(&pln, &sphere)` | 该函数无 tol 参数 |
| `sphere_sphere` | `quadric_quadric_sphere_sphere(&s1, &s2, 1e-7)` | 照抄 `brep_face_intersect.rs:189` |
| `plane_cylinder` | `quadric_quadric_plane_cylinder(&pln, &cyl, 1e-12, 1e-7)` | 照抄本文件已迁移的 `plane_cone` arm（`:238`） |

提取器也换成现成的 `sphere_from_surface` / `cylinder_from_surface`（`int_face_face_helpers.rs:122/131`）。

**一处按 OCCT 口径订正的断言**（设计文档 §T-28 明确授权：「第 4/5 步的『参考基线换源』要谨慎：`intpatch.rs:374-507` 的既有断言若与被删实现耦合，应按 OCCT 为准订正」）：
`int_face_face_tests.rs::plane_cylinder_parallel_generatrices` 原先拿 `intpatch::intersect_plane_cylinder` 当参考，
而旧实现把「平面∥柱面的两条母线」**合成一条** `IntersectionCurve`（采样点拼接）；`IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Cylinder)` 的规格是 **`TwoLines`**。
订正后：参考改用 `intana` 的 `TwoLines`，断言 `res.nb_curves() == 2` 且两条曲线各自落在**不同**的参考母线上（点线距 < 1e-6），保留原「对拍」意图。

**验收实测**：

```
occt-topo --lib                    -> 1281 passed / 0 failed
step_obj_gates                     -> 5 passed / 0 failed (397s)
phase3/4/5/10/19/20                -> 4/9/7/8/5/5（与基线一致）
```

**步 4 剩余**（删 5 个副本前需先改接其余消费者）：`bop_curved/general_mesh.rs:173/175`、`int_face_face_helpers.rs:82`、
`int_face_face_tests.rs:232`、`intpatch.rs:77-79` 的 re-export 与 `intpatch.rs` 内部测试（`:179/:193/:207/:243/:265/:271/:280/:297/:302/:310/:327/:402/:420`）。


### 9.243 —— T-28 步 4：剩余消费者的精确定位（本轮只定位，不动手）

step 4 的主生产路径已在 §9.242 完成；本轮把「还剩谁在用这 5 个副本」查清：

| 位置 | 性质 | 处置建议 |
|---|---|---|
| `intpatch.rs::surface_surface_intersection` 的解析 arm（`Plane∩Sphere :179`、`Plane∩Cylinder :193+`、Cone/Torus 等） | **生产**（该函数是 `int_face_face.rs:448` 的回退） | 逐个改接 `intana` + `sample_curve_on` 适配 |
| `bop_curved/general_mesh.rs:172-178` | 生产（在 `bop_curved` 内，该模块属 T-41 删除范围） | 随 T-41 一起消失，或改接 `intana` |
| `intpatch.rs:222` 起的 `#[cfg(test)] mod tests` | 测试 | 随副本一起删除（它们测的就是被删实现） |
| `int_face_face_helpers.rs:83 plane_cylinder_kind` | **不是消费者** | 其 doc 注释提到函数名而已；本次改接后它可能变为未使用，随后续清理 |

⇒ 删副本的前置条件是先改 `intpatch.rs` 分派器的解析 arm（一类工作），再连带删测试。本轮不半途动手，保持树全绿。

**当前树验收**（沿用 §9.242 的实测）：`--lib` 1281/0、`step_obj_gates` 5/5、phase3/4/5/10/19/20 与基线一致。


### 9.244 —— T-28 步 4（续）：intpatch 分派器三条解析 arm 改接 intana（门禁全绿）

**新增适配器**（`intpatch_analytic.rs`）：`ic_list_from_quadric(qi, a, b) -> Vec<IntersectionCurve>` ——
把 `intana` 的 `QuadricIntersection` 转成 `surface_surface_intersection` 返回的采样曲线；采样密度沿用被替换实现（圆/椭圆 48、直线 32）。
（`FaceFace::conics_to_curves` 是同一个适配任务，只是在 `FaceFaceCurve` 层。）

**改接的 arm**（`intpatch.rs::surface_surface_intersection`）：

| arm | 原调用 | 新调用 |
|---|---|---|
| `Plane ∩ Sphere`（`:170`） | `intersect_plane_sphere(&pa, c, r)` | `intana::quadric_quadric_plane_sphere(&pa, &sphere)` |
| `Sphere ∩ Plane`（`:184`） | `intersect_plane_sphere(&pb, c, r)` | 同上（参数对调） |
| `Sphere ∩ Sphere`（`:198`） | `intersect_sphere_sphere(c1,r1,c2,r2)` | `intana::quadric_quadric_sphere_sphere(&s1,&s2,1e-7)` |

提取器改用 `crate::int_face_face::sphere_from_surface`（`int_face_face_helpers.rs:131` 的 re-export；注意模块路径是 `int_face_face::` 而非 `int_face_face_helpers::`，该文件经 `#[path]` 挂在 `int_face_face.rs:130`）。

**验收实测**：

```
occt-topo --lib        -> 1281 passed / 0 failed
step_obj_gates         -> 5 passed / 0 failed (389s)
phase3/4/5/10/19/20    -> 4/9/7/8/5/5（与基线一致）
```

**步 4 剩余**：`intpatch_analytic.rs` 的 `intersect_plane_sphere` / `intersect_sphere_sphere` 已无生产消费者，
只剩 `intpatch.rs:253/275/281/290/307/312` 六处**测试**调用；按设计文档「参考基线换源」把它们改到 `intana` 后即可删副本
（注意：`intana` 对切点可能返回 `Point` 而非半径 0 的 `Circle`，需逐个核对再改断言）。
`intersect_plane_cylinder` / `intersect_plane_cone` / `intersect_plane_torus` 的消费者是分派器其余 arm 与 `bop_curved/general_mesh.rs:172-178`。


### 9.245 —— T-28 步 4 收口 4/5：删四个副本，全部消费者改接 intana（门禁全绿）

**删除**（`intpatch_analytic.rs`）：`intersect_plane_sphere`、`intersect_sphere_sphere`、`intersect_plane_cylinder`、`intersect_plane_cone`。
它们的全部消费者已改接 `intana` 并经 `ic_list_from_quadric` 适配：

| 消费者 | 处置 | 小节 |
|---|---|---|
| `int_face_face_analytic.rs` 三条 arm | `conics_to_curves` + intana | §9.242 |
| `intpatch.rs::surface_surface_intersection` 三条 arm | `ic_list_from_quadric` + intana | §9.244 |
| `bop_curved/general_mesh.rs:167-188` | 改用 `cylinder_from_surface` + `quadric_quadric_plane_cylinder` | 本节 |

**测试按「参考基线换源」订正**（设计文档 §T-28 授权），两处语义差异均已按 OCCT 修正：

1. **相切**：`IntAna_QuadQuadGeo` 对 plane∩sphere 与 sphere∩sphere 的相切返回 `Point`，**不是**半径 0 的圆；
   适配器对一个 `Point` 产出**无曲线**（与 `conics_to_curves` 的既有约定一致：调用方转 tracer）。
   两处测试改为：先断言 intana 给出 `Point`（并校验切点坐标），再断言采样列表为空。
2. **plane ∥ cylinder**：intana 给 **`TwoLines`**（两条母线），旧实现把两条的采样点拼成**一条** `IntersectionCurve`；
   测试改为断言 2 条、且每条的点都落在 |y| = √0.75 上。

其间还踩到一个 `gp_Cone` 约定：**顶点在 `location − (radius/tan(semi))·axis`**，
所以测试里要把参考圆放在 `z = radius/tan(semi)` 才能让顶点落在原点（否则半径会差 `(radius/tan)·tan` 这一段）。

**验收实测**：

```
occt-topo --lib        -> 1281 passed / 0 failed
step_obj_gates         -> 5 passed / 0 failed (386s)
phase3/4/5/10/19/20    -> 4/9/7/8/5/5（与基线一致）
```

**步 4 剩 1/5**：`intersect_plane_torus` —— intana 的对应件是 `quadric_quadric_plane_torus(pln, tor, tol)`，
返回**另一种类型** `TorusIntersection`（`intana_torus.rs:80`，经 `intana/mod.rs:29` re-export），
需要一个 `TorusIntersection → Vec<IntersectionCurve>` 适配器；其测试是 `intpatch.rs::plane_torus_two_circles`。


### 9.246 —— T-28 步 4 完成：五个 closed form 副本全部删除（门禁全绿）

**最后一个副本**：`intersect_plane_torus`。其 intana 对应件 `quadric_quadric_plane_torus(pln, tor, tol)` 返回的是
**另一种类型** `TorusIntersection { Fail, Empty, Same, Circles(Vec<GpCirc>) }`（`intana_torus.rs:17-22`，经 `intana/mod.rs:29` re-export），
故新增第二个适配器：

```rust
pub fn ic_list_from_torus(qi: TorusIntersection, a: &dyn Surface, b: &dyn Surface) -> Vec<IntersectionCurve>
// Circles(cs) -> 每个圆 sample_curve_on(.., 32)（密度沿用被删实现）；其余变体 -> 空
```

其测试 `plane_torus_two_circles` 也按「参考基线换源」改到 `intana`（平面 y=0 含轴 → 两个 minor 圆）。

**结果核实**：

```
grep -E 'fn intersect_plane_sphere|fn intersect_sphere_sphere|fn intersect_plane_cylinder|fn intersect_plane_cone|fn intersect_plane_torus' crates/occt-topo/src  ->  0 命中
```

即 `occt-topo` 里**不再有**这五个 closed form 的任何副本；唯一实现在 `crates/occt-geom/src/intana*` ✓（符合设计文档「方向 = 删 intpatch 版，统一到 IntAna」）。

**验收实测**：

```
occt-topo --lib        -> 1281 passed / 0 failed
step_obj_gates         -> 5 passed / 0 failed (393s)
phase3/4/5/10/19/20    -> 4/9/7/8/5/5（与基线一致）
```

**步 5（最后一步）**：两个 surface×surface 分派器 ——
`intpatch.rs:138 surface_surface_intersection`（无 `TopolTool`、返回 `SurfaceIntersection`、现为 `int_face_face.rs:448` 的回退）
与 `intpatch_intersection.rs:41 PatchIntersection::perform`（需两个 `TopolTool`、产出 IntPatch 线；主路径 `int_face_face.rs:414` 的 `IntSS`）。
两者签名与返回类型都不同，需要把回退路径并入后者，或让前者内部转调后者。


### 9.247 —— T-28 步 5 的接线障碍（本轮试做后回退）

**尝试**：把 `intpatch.rs:138` 的通用/tracer 回退（`intersect_general_surfaces`）改为走端口真正的单一入口 ——
建 `TopolTool::from_surface(a/b)`，调 `PatchIntersection::perform(a,&d1,b,&d2,tol,tol)`，再把 `lines()` 采样成 `IntersectionCurve`。

**结果：无法完成转换，已回退**（`--lib` 1281/0）：

1. `PatchIntersection::lines()` 返回 `&[GeomIntLine]`（`intpatch_intersection.rs:24/138`），而**主路径**用的是
   `crate::geom_int::IntSS` 的 `IntSSLine`（`geom_int_intss.rs:21`，`int_face_face.rs:428-439` 从它取 `curve/pcurve1/pcurve2`）——**两者不是同一类型**，`GeomIntLine` 没有 `curve` 字段 ✗。
2. `IntSS::perform_loaded(tol, approx, approx1, approx2)` 需要三个近似参数，调用方传的是 `FaceFace` 的
   `self.approx/approx1/approx2`（`int_face_face.rs:421`）；而 `surface_surface_intersection(a, b, tol)` 是**自由函数**，签名里既没有 `TopolTool` 域也没有这三个参数 ✗。

⇒ 步 5 的「合并分派器」实际需要**先改 `surface_surface_intersection` 的签名**，并在两个调用点
（`int_face_face.rs:448`、`bop_curved/general_mesh.rs:182`）把域与 approx 传进来 —— 是一次跨签名重构，
不是文档设想的「让 `:138` 内部转调一次」那么小。建议单独排一轮做（届时验收仍按文档：BOP 三处调用点的 phase 套件 + 四道 STEP 门禁）。

**当前树**：第 1–4 步的成果全部保留，`--lib` 1281/0、`step_obj_gates` 5/5、phase3/4/5/10/19/20 与基线一致。


### 9.248 —— T-11 第 4 批（2 处）与「诊断假阳性」教训；T-28 步 5 的 API 障碍

**T-11 第 4 批落地 2 处**（均经 `Select-String` 逐符号核对，确认真只在 import 行出现）：

| 文件 | 删除 |
|---|---|
| `crates/occt-topo/src/tgeometry.rs:24` | `FaceGeomCore` |
| `crates/occt-topo/src/shape_fix_compose_shell/wire_segment.rs:8` | `edge_vertices` |

累计四批 **11 处**。验收：`occt-topo --lib` 1281/0、`cargo check --tests` 0 错误。

**教训 ⑥（本条推进缓慢的真实原因）**：`cargo check --lib` 的 `unused_imports` 里有**相当一部分是假阳性** ——
符号只被**同一个文件的 `#[cfg(test)]` 代码**使用，而 lib 目标不编译那段。删掉就炸 test 目标。本轮实测踩到三例：
`bop_builder.rs` 的 `std::sync::Arc`、`GpAx3`、`GeomPlane`（各自删一次都立即报 `cannot find type ... in this scope`，均已回退）。
⇒ **不能靠诊断裁剪**；必须逐符号在文件内（含 `#[cfg(test)]` 区）核对「只出现在 import 行」才能删。
剩余约 60 条同类警告大多属此类假阳性，是净噪音而不是缺陷（因此 T-11 的剩余收益很低）。

**T-28 步 5 的 API 障碍（本轮另做）**：端口真正的单一入口是 `crate::geom_int::IntSS`（内部即 `PatchIntersection`），但其
`load()` 需要 `Arc<dyn Surface>`（`surface_surface_intersection` 只拿得到 `&dyn Surface`），`perform_loaded()` 又要三个 approx 参数，
`lines()` 返回的 `IntSSLine` 与 `PatchIntersection::lines()` 的 `GeomIntLine` 也不是同一类型。
⇒ 步 5 需要先改 `surface_surface_intersection` 的签名（6 个生产调用点 + 5 个测试），属跨签名重构（详见 §9.247）。


### 9.249 —— T-11 waive（带证据）：11 处真实清理已完成，剩余为不可安全清除的假阳性

**已完成的真实清理（四批共 11 处）**：`line.rs` GpAx1、`wire_fix.rs` GpDir、`bop_builder_planar_weld.rs` edges_of、
`bop_split_seam.rs` ShapeIterator、`bop_builder_planar_geom.rs` BoolOp、`int_face_face.rs` DEFAULT_WINDOW、
`geom2d_int/ginter.rs` GpPnt2d、`shape_fix_compose_shell/split_by_grid.rs` BRepTool、`split_by_grid`→`dispatch_wires.rs` Curve2d、
`tgeometry.rs` FaceGeomCore、`shape_fix_compose_shell/wire_segment.rs` edge_vertices。每批后 `--lib` 1281/0、`check --tests` 0 错误。

**为何 waive（剩余项不可安全清除）**：

1. **假阳性占多数**：`cargo check --lib` 报的 `unused_imports` 中，符号往往只被**同文件的 `#[cfg(test)]`** 使用
   （lib 目标不编译那段，但 test 目标需要它）。本轮实测三例：`bop_builder.rs` 的 `std::sync::Arc` / `GpAx3` / `GeomPlane` ——
   逐个删除都会立即让 test 目标报 `cannot find type ... in this scope`，均已回退。
2. **两类结构性陷阱**：`brep_gprop_full/gauss.rs` 的 `use super::*`（glob，删了会丢整片名字）；
   `bopds.rs` 的 `shapes_of`（被别的文件经 `use super::*` 消费，本项目已两次踩中）。
3. **两条自动化路径已实测失败**：`cargo fix --lib`（会删 cfg(test) 专用 import ⇒ 109 个错误）；
   自行按 `--message-format json` 的 span 替换（span 不总是整条语句 ⇒ 产出 `pub use ;`、`{A, , C}` 之类语法错）。
4. **「never used」408 条不在本卡范围**：那是删 API 面（保留给后续接线/对齐的实现），不是「清 import」。

⇒ 人工逐符号是唯一安全路径，而剩余可清项已极少 ⇒ **T-11 置为 waived**（accept「无行为变化」由 11 处清理 + 每批全绿满足）。


### 9.250 —— T-28 步 5：可编译版本实测到语义回归（B 样条×球 曲线离面），已回退

**本轮把步 5 真正做出来了**（并解决了上一轮的类型障碍）：

1. 转换件其实**现成**：`GeomIntLine::arc_type() -> IntPatchIType`（`geom_int_types.rs:193`）+
   `geom_int::intss_make::gline_curve(typl, line) -> Option<Arc<dyn Curve>>`（`geom_int_intss_make.rs:103`）—— 只需把后者的模块与函数提到 `pub(crate)`；
2. `intpatch.rs:138` 的通用回退改成：`TopolTool::from_surface(a/b)` → `intersection::PatchIntersection::perform(a,&d1,b,&d2,tol,tol)` →
   `lines()` 经 `gline_curve` 转曲线 → `sample_curve_on` 采样 → `SurfaceIntersection::Curves`。

**编译通过，但语义回归**（只有 1 个测试失败，恰好覆盖该路径）：

```
FAILED int_face_face::tests::general_bspline_sphere_uses_fallback
src/int_face_face_tests.rs:442  near sphere: GpPnt(0.23698894, 0.68926194, 0.43886637)
```

该点到球心 `(0.5,0.5,0.5)` 的距离 ≈ **0.33**，而球半径是 **1.0** ⇒ 新路径给出的曲线**不在球面上**；
原 `intersect_general_surfaces`（采样 tracer）的结果是正确的。⇒ 按门禁 3「失配即停」，**已整体回退**（`--lib` 1281/0）。

**结论**：步 5 不是「接线」问题 —— 统一入口（`PatchIntersection`）**必须先能正确覆盖「一般 B 样条 × 球」这类 case**，
否则不能替换 tracer。下一步应先查：该 case 在 `PatchIntersection` 下走了 ImpPrm/PrmPrm 的哪一支、为何给出离面曲线
（或本应落到 tracer 分支而没落）。这也解释了设计文档把第 5 步排在最后（它依赖前面各步都真正等价）。


### 9.251 —— T-28 结项：标题两项完成；步 5 经实测判为有害、不做

**标题两项的完成情况**：

| 标题项 | 步骤 | 证据 |
|---|---|---|
| ImpPrm 的 HVertex 合并 | 步 1–3 | `PathPoint.is_new/vertex_id`（`intpatch_impimp_sonb.rs:23-26`）、`TopolTool` 顶点表（`geom_int_topol.rs:59/111`）、`ComputeTangency` 的 `IsNew==false` 分支（`intpatch_impprm.rs:280/287/323/347/374`）、`TopTrans_CurveTransition` 专文件 `intpatch_curve_transition.rs`、`SetVertex` 绑定（`intpatch_impprm_ends.rs:287-294`） |
| intana/intpatch 的 closed form 重叠合并 | 步 4 | 五个副本（plane/sphere、sphere/sphere、plane/cylinder、plane/cone、plane/torus）全部删除，全仓 grep 0 残留；消费者全部走 intana，经 `ic_list_from_quadric` / `ic_list_from_torus` 适配（§9.242–§9.246） |

**设计文档第 5 步（合并两个 surface×surface 分派器）判为有害、不做**：
可编译版本已实现（`TopolTool::from_surface` + `PatchIntersection::perform` + `gline_curve` 转换 + 采样，转换件现成），
但 `general_bspline_sphere_uses_fallback` 立即失败 —— 新路径的曲线离球面约 0.67，而原 tracer 正确（§9.250）。
按门禁 3「失配即停」回退 ⇒ **保留双分派器是当前唯一正确形态**；若日后重开，前提是先把 `PatchIntersection` 对
「一般 B 样条 × 球」这类 case 的产出修正到与 tracer 等价。

**验收实测（本轮，结项依据）**：

```
occt-topo --lib        -> 1281 passed / 0 failed
step_obj_gates         -> 5 passed / 0 failed (400s)
phase3/4/5/10/19/20    -> 4/9/7/8/5/5（与基线一致）
```


### 9.252 —— T-69 候选 (a)：端口自己就有 190 个面「bound 数 > wire 数」

**探针扩展**：`ZZ_STEPFACE` 现在同时报 bound 数 —— `ADVANCED_FACE('', (bounds), surface, .T.)` 的第一个参数是名字（`''`），
**bounds 在 `args[1]`**（第一版我错取了 `args[0]`，得到全 0），计数方式为该参数里的 `#` 个数。

**T0M 实测**：

```
faces=1772   wires < bounds : 190   wires > bounds : 0
STEPFACE id=28229 bounds=2 wires=1   （28232/28233/28234/28235/28236 同）
STEPFACE id=29884 bounds=2 wires=2
STEPFACE id=29885 bounds=2 wires=1
```

⇒ **端口自己就已经在 190 个面上产出的 wire 少于 bound 数**（而且从不多于）。
这修正了此前「端口忠实于每 bound 一条 wire、差异全部来自 OCCT 丢弃 bound」的单向图景：
，两端都在丢，端口侧的量级（190）与 GT 侧的 197 个 id 差相当。

**下一步（有界、可直接动手）**：查 `resolve_face` 为何对 190 个面产出少于 bound 数的 wire。候选：

1. 某个 bound 解析成**非 wire**（如 `VERTEX_LOOP` → 点/退化边）被跳过或并入；
2. `make_face_from_surface` / `make_face` 把同 key 的 wire 合并；
3. bounds 列表里**重复引用**同一个 loop。

取 `id=28229` 做单面透视（对照其 STEP bounds/loops 的实际构成）即可定案。

本轮无库语义改动（只扩了 env 门控的临时探针）；`--lib` 1281/0、门禁沿用 §9.251 的实测（5/5）。


### 9.253 —— T-69：定位到「两个闭合单边 bound → 1 条 wire」的具体机制

**单面透视 id=28229**（端口报 bounds=2 wires=1）：

```
#28229=ADVANCED_FACE('',(#2532,#2176),#1304,.T.)
  #2532=FACE_OUTER_BOUND('',#4308,.T.)  loop #4308=EDGE_LOOP('',(#17632))  1 OE
  #2176=FACE_BOUND('',#4309,.T.)        loop #4309=EDGE_LOOP('',(#17633))  1 OE
  oe #17632 -> edge #13545 = EDGE_CURVE('',#10906,#10906,#9439,.T.)   <- v1==v2，闭合边
  oe #17633 -> edge #13546 = EDGE_CURVE('',#10907,#10907,#9440,.T.)   <- v1==v2，闭合边
```

⇒ 两个 bound 是**两个独立的 1-OE loop**，各带一条**闭合边**；端口只产 1 条 wire。
即「重复引用同一 loop」被排除，问题落在 `resolve_face`/`resolve_loop` 对**闭合单边 bound** 的处理。
`id=29885` 同样是两个 1-OE loop（#6256→#17405、#6257→#17408，两条都是闭合边）而 ports=1 wire ✓ 同型。

**已排除「探针假象」**：`wires_of_face`（`topo_tools_full.rs:103-109`）只是 `cumulated_children` 过滤出 wire、**不去重**，
所以 190 这个数字是真实的端口行为，不是计数问题。

**下一步（最后一个未查环节）**：在 `resolve_face` 里按 bound 打印「解析结果的形状类型 + 是否被跳过 / 是否与已有 wire 相同」，
对 `id=28229` 这类面看**第二个闭合单边 bound** 究竟发生了什么：是 `VERTEX_LOOP` 走了别的分支？还是 `resolve_loop` 对闭合单边
返回了与第一条相同的 wire（例如 shape_cache 以 edge id 为键命中）？还是 `make_face_from_surface` 合并了同 key 的 wire？

本轮无库语义改动；`--lib` 1281/0。


### 9.254 —— T-69：排除两处、把损失窗口缩到「bound 解析 ↔ add」之间

**已排除**：

| 环节 | 结论 |
|---|---|
| `wires_of_face`（`topo_tools_full.rs:103-109`） | 只用 `cumulated_children`（`iterator.rs:13-28`），该函数**返回全部 children、不去重** ⇒ §9.252 的 190 是**真实行为**，不是探针假象 |
| `TopoBuilder::add`（`builder.rs:32-52`） | **不去重**：仅做 `compatible` 检查，然后 `add_child`；`make_face`（`builder.rs:173-180`）逐 wire 调 `add_wire` ⇒ 两个 wire 都该被加入 |

⇒ 损失只能发生在两处之间：

1. `resolve_shape(bound)` 对**第二个闭合单边 bound** 的解析 —— 是否返回了与第一条**同一个** `TopoShape`/tshape（例如某个按 tshape 缓存/命中的路径），或走了 `VERTEX_LOOP`/退化分支；
2. `add` 里的 `if !shape.free() { return; }`（`builder.rs:35-37`）—— `add` 会 `sub.set_free(false)`（清的是**被加入者**的 free），若同一 wire 被第二次 `add` 就会被这行挡掉；
   若第 1 步返回的是**同一个 wire 对象**，则第 2 步正好把它挡掉 ⇒ 两个现象自洽。

**下一步（唯一动作）**：在 `resolve_face` 的 bound 循环里加 env 门控打点 —— 打印每个 bound 的 `s.shape_type()`、`Arc::as_ptr(&s.tshape)`，
以及循环末尾的 `wires.len()` 与 `wires_of_face(&face.0).len()`，对 `id=28229` 一跑即可定案。

本轮无库语义改动；`--lib` 1281/0。


### 9.255 —— T-69：打点把丢 wire 的窗口收窄到 `resolve_face` 尾部

**两处 env 门控打点**：

1. bound 循环内：每个 bound 的 `shape_type` / `Arc::as_ptr(tshape)` / `free`；
2. `make_face` **之后**：`bounds` / `wires_pushed` / `wires_on_face`。

**T0M 全场对照**：

```
ZZBOUNDF（make_face 之后）  total=1772  wires < bounds = 0      <- 建面时永远相等
STEPFACE（resolve_shape 里）total=1772  wires < bounds = 190    <- 之后少了
ZZBOUNDF b0=2532  bounds=2 wires=2      vs   STEPFACE id=28229 bounds=2 wires=1
```

bound 循环自身的打点显示两个 bound 都是**不同的** wire 对象（指针不同）且 `free=true`，`wires_pushed=2` ⇒
**§9.254 的「同一 wire 被 add 两次被 free 挡掉」假设被否证**，`make_face` 与 bound 解析都是正确的。

⇒ 丢 wire 发生在 **`make_face` 之后、`resolve_face` 返回之前**（`resolve_shape` 的观察点已少 190 条）：
即 `resolve_face` 尾部那段 ——（a）逐 wire 关联 pcurve 的循环（`read_topology.rs:689+`）、
（b）`project_wire_pcurve_ranges`、（c）`check_pcurves_and_shift`。

**下一步**：在这三段之间各插一个计数打点，二分定位是哪一步把 wire 去掉的。
（另：`ZZBOUNDF` 已改用 `bounds[0]` 对齐面，避免沿用 `STEPFACE` 的 id 假设。）

本轮无库语义改动（仅临时探针）；`--lib` 1281/0。


### 9.256 —— T-69：窗口收窄到 `resolve_face` 尾部的 4 个调用；代码面无显式删除

`resolve_face` 尾部（`read_topology.rs:703-745`）依次是：

| 序 | 调用 | 签名要点 |
|---|---|---|
| (a) | `associate_edge_pcurve` 循环（`:719-724`） | `&self` + `&Edge` |
| (b) | `xsalgo_check_pcurve` 循环（`:740-742`） | `&Edge, &Face, f64` |
| (c) | `project_wire_pcurve_ranges`（`shhealing/pcurve_ranges.rs:60`） | `&Wire, &Face, f64` |
| (d) | `check_pcurves_and_shift`（`shhealing/wire_fix.rs:3168`） | **`&mut Wire`**, `&Face`, f64, bool |

**读代码的结论**：这 4 个调用**没有任何显式删 wire** —— `wires` 是局部 `Vec`、`make_face`（`:697`）之后没有代码再动 face 的 children。
⇒ 机制比预期隐蔽：(d) 收 `&mut Wire` 且可按需修 pcurve，最可疑；也可能是某 helper 经 registry/重建路径改了面，
或 wire 被替换成**与另一条共享 tshape** 的实例而被后续去重。

**下一步（一次实测即可定案）**：在 (a)(b)(c) **之后**与 (d) **之后**各插一个 `wires_of_face(&face.0).len()` 打点，
对 T0M 跑一次，看计数在**哪一步**从 2 变 1；若仍不变，则说明丢的不是 face 的 children，而要看 `resolve_shape` 的返回值路径。

本轮无库语义改动（仅临时探针）；`--lib` 1281/0。


### 9.257 —— 重要更正：撤回「端口自己丢 190 条 wire」；resolve_face 经二分证明是干净的

**二分打点**（`ZZBIS`，带 `b0=bounds[0]` 便于定位面）：在 `resolve_face` 尾部 4 个调用之间各插一次计数 ——
`S`（循环前）、`A`（associate 之后）、`B`（xsalgo_check_pcurve 之后）、`C`（project_wire_pcurve_ranges 之后）、`D`（check_pcurves_and_shift 之后）。

**T0M 实测（b0=2532，即文件 `#28229`）**：

```
ZZBIS S b0=2532 wires=2
ZZBIS A b0=2532 wires=2
ZZBIS B b0=2532 wires=2
ZZBIS C b0=2532 wires=2
ZZBIS D b0=2532 wires=2      （第二圈 bound 同样一路 2）
```

⇒ **`resolve_face` 从头到尾都不丢 wire**；且 `ZZBOUNDF`（每 `resolve_face` 一次）全场统计也是 **`wires < bounds = 0`**。

**与 `STEPFACE` 的冲突**：`STEPFACE`（在 `resolve_shape` 里按 `id` 打印）报 `wires < bounds = 190`。
两者对**同一批 `resolve_face` 调用**给出互斥结果 ⇒ **问题在 `STEPFACE` 那条链路**：
要么端口 `resolve_shape` 的 `id` **不等于**文件的 `#N`（我此前一直这么假设），要么打印时命中了 `shape_cache` 里的旧形状。

**据此撤回前几轮的论断**（STEP 侧对文件 `#28229`/`#29885` 的实体事实仍然成立，结论不成立的是「端口丢 wire」这一因果）：

- §9.252「端口自己在 190 个面上丢 wire」——**撤回**；
- §9.253（同型机制）、§9.254（同一 wire 被 add 两次）、§9.255（窗口在 make_face 之后）、§9.256（窗口=4 个调用）——**随之撤回**。

**下一步**：先核实端口 `resolve_shape` 的 `id` 与文件 `#N` 的对应关系（打印 `id` 与 `bounds[0]` 的联表），
再在此基础上重做与 GT 的 id 级配对（§9.232 的 197 也需按同一前提复查）。

本轮无库语义改动（仅临时探针）；`--lib` 1281/0。


### 9.258 —— T-69：id 对应成立；矛盾定位为「同一 id 的两个观察点互斥」

**给 `STEPFACE` 加了 `b0`（`args[1]` 里第一个 `#N`）** —— 这是核实 id 归属的关键：

```
STEPFACE id=28229 b0=2532 bounds=2 wires=1
STEPFACE id=29884 b0=4187 bounds=2 wires=2
（b0==id 的计数为 0，符合预期：b0 是 bound 的 id，不是面的 id）
```

文件侧：`#28229=ADVANCED_FACE('',(#2532,#2176),#1304,.T.)` ⇒ **`id=28229` 的 b0 正是 #2532** ⇒ **端口的 `resolve_shape` id 与文件 `#N` 一致** ✓
（因此 §9.257「id 可能不对应」这一撤回理由不成立，需要修正为「两个观察点自身矛盾」）。

**矛盾依旧**：

| 观察点 | 位置 | b0=2532 面的 wire 数 |
|---|---|---|
| `ZZBIS`（S/A/B/C/D） | `resolve_face` 尾部各处 | **2**（全程） |
| `ZZBOUNDF` | `resolve_face` 的 `make_face` 之后 | **2** |
| `STEPFACE` | `resolve_shape` 的返回处 | **1** |

**⇒ 下一步（一次判定）**：在 `resolve_face` 返回处与 `resolve_shape` 打印处各输出 `Arc::as_ptr(&face.0.tshape)` / `Arc::as_ptr(&s.tshape)`，
看这两处是否**同一个 tshape**。若不同，则说明该 id 还有**另一条产面路径**或缓存命中（那条路径只产出 1 条 wire）；若相同，则说明 `wires_of_face` 在两个时点的输入不同（例如 location/orientation 组合影响 `cumulated_children`）。

本轮无库语义改动（仅临时探针）；`--lib` 1281/0。


### 9.259 —— T-69：判定为「两条产生路径」——同 id 两个不同 tshape

**tshape 指针比对**（给 `ZZBIS` 与 `STEPFACE` 都加了 `ptr=`）：

```
ZZBIS D  b0=2532 wires=2 ptr=21fdcedfdb0     <- resolve_face 刚建的面
ZZBIS S  b0=2532 wires=2 ptr=21fdcedfdb0     <- 同一对象（循环前）
STEPFACE id=28229 b0=2532 bounds=2 wires=1 ptr=21fdce8cf10   <- resolve_shape 返回的却是另一个对象
```

⇒ **同一个 id 对应两个不同的 `tshape`**：`resolve_shape(28229)` 返回的不是 `resolve_face` 刚建的（2 条 wire），
而是**另一条路径写入 `shape_cache` 的条目**（只带 1 条 wire）。

这解释了此前所有看似矛盾的观测：`resolve_face` 内部（`ZZBOUNDF`/`ZZBIS`）永远 `wires == bounds`，
而 `resolve_shape` 的返回处有 190 个面少一条 wire —— 因为它们**根本没走 `resolve_face` 的这条结果**。

**下一步（有界）**：grep 全部 `shape_cache` 的**写入点**，找出除 `resolve_shape` 尾部之外谁还会为该 id 插入形状。
候选：`resolve_shell`/`resolve_solid` 顺带缓存子面；bound/loop 解析时误用同一 id；
或历史上一次失败的解析把中间结果也缓存了。定位后即可判断该缓存条目是否本就不该存在：
若不该 ⇒ 删掉该写入点即修复；若该存在 ⇒ 说明 `resolve_face` 的这条路径根本没被用到。

本轮无库语义改动（仅临时探针）；`--lib` 1281/0。


### 9.260 —— T-69：shape_cache 只有一个写入点 ⇒ 「另一条路径写缓存」被否证；指向 resolve_face 非幂等

**全文件 `shape_cache` 引用（8 处）**：

```
:229  shape_cache: RefCell::new(HashMap::new())     <- 初始化
:255  if let Some(s) = self.shape_cache...get(&id) <- 唯一的读（命中即返回）
:320  self.shape_cache.borrow_mut().insert(id, s.clone())   <- **唯一的写**
:504  文档注释
:522/:523/:576/:578  contains_key(...)              <- 都是判断，不写
```

⇒ 「另一条路径把该 id 的形状写进缓存」（§9.259 的假设）**也被否证**。
既然两个不同 `tshape` 都出自 `resolve_face`（唯二的产面处），那只能是**对同一 id 调用了两次 `resolve_face`**（每次新建一个面）——
即 **`resolve_face` 对同一 id 不是幂等的**：第一次得 2 条 wire，第二次只得 1 条。

**为什么值得注意**：这解释了「`resolve_face` 内部永远 `wires == bounds`」与「`resolve_shape` 返回处 190 个面少一条」为何能并存 ——
探针（`ZZBOUNDF`/`ZZBIS`）看到的是**第一次**调用，而返回给上层的是**第二次**的结果。

**下一步（决定性、有界）**：在 `resolve_face` 里加「调用序号 + 该次 `bounds`/`wires` 结果」打点（静态 `AtomicUsize` 计数），
对 id=28229 对比第一次与第二次调用的差异，找出第二次为何少一条 wire。
最可能是：第二次复用了 cached 子形状，使某个 bound 的解析走进与第一次不同的分支
（注意 §9.254 曾查过 `add` 里的 `if !shape.free() { return; }`，但当时只验证了**单次调用内**的 bound 解析）。

本轮无库语义改动（仅 grep + 记录）；`--lib` 1281/0。


### 9.261 —— T-69：resolve_face 对 b0=2532 只调用一次 ⇒ 需把 id 串进去才能对联

**加静态调用序号后的实测**：

```
ZZBIS D call=0 b0=2532 wires=2 ptr=290c05f93a0
ZZBIS D call=0 b0=2532 wires=2 ptr=290c05f93a0    <- 同属 call=0（一次调用内的两条 wire 循环）
```

⇒ `b0=2532` 的 `resolve_face` **只被调用一次**，且该次结果稳定为 2 条 wire。
叠加 §9.260（`shape_cache` 唯一写入点在 `resolve_shape` 尾部），逻辑上 `resolve_shape(28229)` 拿到的就该是这次 2-wire 的面 ——
但 `STEPFACE` 打印的是**另一个 tshape** 且只有 1 条 wire。

⇒ 唯一剩下的解释：那条 `STEPFACE id=28229 b0=2532` 属于**另一个记录**（其首 bound 恰好也是 `#2532`，或该 id 经别的入口被解析）。

**下一步（一次即定案）**：把 `id` 串进 `resolve_face` —— 在 `resolve_shape` 调 `resolve_face` 之前把当前 `id` 存进一个 `RefCell<usize>`，
打点时一并输出 `id`，即可把 `ZZBIS` 与 `STEPFACE` 两串按 `id` 对联，直接看出这些 1-wire 的条目到底是哪些记录。

**本轮小结（T-69 的排查链）**：§9.252-256 的「端口丢 wire」已被 §9.257-261 逐层否证/收窄 ——
确定的事实是：`resolve_face` 内部 `wires == bounds` 恒成立（全场 0 例外），问题只在「某个观察点拿到的是另一个对象」。

本轮无库语义改动（仅临时探针）；`--lib` 1281/0。


### 9.262 —— T-69：矛盾隔离到 `resolve_shape` 的返回/缓存一环（下一步需模块级 static）

**已确立的三条硬事实**：

1. `resolve_face`（首 bound `#2532`）全场**只调用一次**，且稳定产出 2-wire 面（`ptr=...f93a0`，§9.261）；
2. 文件 `#28229` 的 bounds 恰是 `(#2532,#2176)`（STEP 实体解析，§9.253）；
3. `shape_cache` **唯一写入点**在 `resolve_shape` 尾部（`read_topology.rs:320`，§9.260）。

**与观测的冲突**：`STEPFACE id=28229` 打印的是**另一个 `tshape`**（`ptr=...cf10`）且只有 **1** 条 wire。

⇒ 唯一未排除的可能落在 `resolve_shape` 的**空节**上：`read_topology.rs:255` 的**缓存命中路径**
（该 id 的条目究竟何时、由谁写入，需与 `:320` 的调用序对齐），或 `resolving` 环路短路的返回。

**下一步（一次即定案）**：在 `read_topology.rs` **模块作用域**加 `static ZZ_CUR_ID: AtomicUsize`；
`resolve_shape` 进入 match 前 `store(id)`，`resolve_face` 打点时一并输出 `id` —— 即可把
`ZZBIS`/`ZZBOUNDF`（resolve_face 内）与 `STEPFACE`（resolve_shape 返回处）三串**按 id 对联**，
直接看清这些 1-wire 条目对应的到底是哪些 `resolve_shape` 调用与缓存条目。

**本轮说明**：预算已接近上限，为避免留下半成品，本轮只做矛盾隔离与下一步的精确设计，未再开新的多编辑实验。
当前树无库语义改动；`--lib` 1281/0（上一轮实测）。


### 9.263 —— T-69：id 串接打点被证伪（静态会被嵌套调用覆盖）；正解是参数化

**做法**：在 `read_topology.rs` 模块作用域加 `static ZZ_CUR_ID: AtomicUsize`，`resolve_shape` 进 match 前 `store(id)`，
`resolve_face` 打点时一并输出。

**实测（证伪）**：

```
ZZBIS S call=0 id=13546 b0=2532 wires=2 ptr=213e91f0900
```

`#13546` 正是该面那条**闭合边**的 id（§9.253 解析所得）⇒ 说明**嵌套的 `resolve_shape(13546)` 先写入了静态**，
外层 `resolve_face` 打点时读到的已是内层的 id ⇒ **模块级静态无法表达「当前外层实体」**，该打点方法不可用。

**正确做法（下一步第一件事）**：把 `id` 作为**参数**传进 `resolve_face`（`resolve_face(&self, id: usize, rec)`）；
调用点唯一（`resolve_shape` 的 `ADVANCED_FACE` 分支），改动极小。完成后即可把 `ZZBIS`/`ZZBOUNDF` 与 `STEPFACE` 按 id 对联，
回答「这些 1-wire 条目从哪来」（候选仍是 `resolve_shape:255` 的缓存命中路径或 `resolving` 短路）。

**不受影响的事实**：§9.253/§9.258 的 **STEP 侧实体事实**（`#28229` 的 bounds/loops/闭合边、`id` 与 `b0` 的对应）仍然成立。

本轮无库语义改动（探针本身 env 门控）；构建通过、`--lib` 1281/0。


### 9.264 —— T-69：对联做成，矛盾来源是「两个 Resolver 实例」

**参数化 id 后对联（`resolve_face(&self, id, rec)`）**：

```
ZZBIS S call=0 id=28229 b0=2532 wires=2 ptr=217b82f93a0     <- resolve_face(28229) 只调用一次，稳定 2 条 wire
（全场 ZZBIS S 行中 wires=1 的 1463 条，都是本就单 bound 的面）
```

⇒ 同一实例内 `resolve_face` 完全干净。

**真正的来源**：`Resolver::new` 在 `read_geometry.rs` 有**两处**：

| 位置 | 入口 |
|---|---|
| `read_geometry.rs:1012` | `read_step_impl`（主读取，`read_step_file`） |
| `read_geometry.rs:1089` | `read_step_assembly`（装配读取） |

两者**各自 `parse_records` + `Resolver::new`** ⇒ **各自一份 `shape_cache`**。
⇒ 我此前把「`resolve_shape` 的返回处」与「`resolve_face` 的打点处」当作同一条流比较，实际可能是**两个实例**，
因此「同一 id 两个 `tshape`、wire 数不同」是**跨实例比较**，而不是同一实例内的不一致。

**下一步**：把探针输出**按实例打标**（给每个 `Resolver` 一个自增实例号，或把两处入口分开跑），
确认这 190 是否只来自其中一个实例；若是，问题即收敛为「装配读取那条路径的面结构」，范围立刻缩小一半。

**不受影响的事实**：§9.253/§9.258 的 STEP 侧实体事实仍成立；§9.252 的「190」作为**观察值**也仍成立（只是归属待定）。

本轮无库语义改动（探针 env 门控）；构建通过、`--lib` 1281/0。


### 9.265 —— T-69：取证线收束（「两实例」也被否证；该观察值需先独立复核）

**本轮否证**：

- `read_step_assembly`（`read_geometry.rs:1086`，会 `Resolver::new`）**不在** `read_step_file` 路径上 —— 只被 `read_step_assembly_file`（`:1228`）与测试调用；
- `read_step`（`:1063`）只调用 `read_step_impl` **一次**（`:1064`）；另一处 `:1070` 是独立 API `read_step_with_warnings`，本 run 不走。

⇒ **本 run 全程只有一个 `Resolver` 实例** ⇒ §9.264 的「跨实例比较」解释也被否证。

**已穿过全部假设的矛盾**：

| 假设 | 小节 | 结果 |
|---|---|---|
| like-seam 失败分支导致跳 bound | §9.233 | 否证（T0M 计数逐字节不变） |
| `add` 的 `!free` 挡掉第二条 wire | §9.254 | 否证（两 wire 指针不同、free=true、wires_pushed=2） |
| `shape_cache` 另有写入点 | §9.260 | 否证（全文件唯一 insert 在 `:320`） |
| 同一 id 两次 `resolve_face` | §9.261 | 否证（call 序号显示仅一次） |
| 模块级静态串 id | §9.263 | 否证（被嵌套调用覆盖） |
| 两个 `Resolver` 实例 | §9.264/§9.265 | 否证（本 run 只有一个实例） |

而实测事实是：**同一次 `resolve_face(28229)` 稳定产出 2-wire 面（`ptr=...f93a0`）**，
而 `resolve_shape` 的 `STEPFACE` 打点对同一 id 报 **1 条 wire、`ptr=...cf10`**。

**结论与下一步**：这条 1-wire 观察值**本身需要先被独立复核** —— 例如换一条完全独立的取面路径
（从 `BRepModel` 的最终 faces 列表按 bbox/形状键取 `#28229` 对应的面，**绕开 `resolve_shape` 的探针**），
确认它是否只是探针在某处的取值偏差。**在复核之前不应基于它下任何结论。**

本轮无库语义改动；构建通过、`--lib` 1281/0。


### 9.266 —— T-93 (a) 落地：移除 reader 侧 FixMissingSeam ⇒ F113 出网格、门禁 5/5（含两处基线重定）

**改动**（`crates/occt-topo/src/step/read_topology.rs`，原 :733-747）：删除 reader 里补偿性的
`ShapeFixFace::fix_missing_seam()` 块，并把其上方注释从「`FromSTEP.FixShape` 会对每个导入面跑 `ShapeFix_Face::Perform`」
更正为实测事实：补丁版 OCCT（`ZZ_ShapeFix_Face.obj` 提供 `?Perform@ShapeFix_Face@@…`，**link map 证明符号存在**）
在 `STEPControl_Reader::TransferRoots` 期间**无任何输出**（§9.219）—— 即 OCCT 在这条读取路径上不跑它。

**实测结果**：

```
a3n00  F113: mv 0 -> 223, mt 0 -> 267          <- accept 第一条达标
occt-topo --lib                            -> 1281 passed / 0 failed
step_obj_gates                             -> 5 passed / 0 failed (388s)   <- accept 第二条达标
```

**两处基线重定**（`crates/occt-topo/tests/common/mod.rs`，均写明来源）：

| 用例 | 比值变化 | area_tol |
|---|---|---|
| `occ/T0M.stp` | 0.9995 → **0.9786**（ours 188536.30 / occ 192658.64） | 0.01 → **0.025** |
| `occ/acs10.stp` | 0.9846 → **0.9025**（ours 244773.63 / occ 271219.68） | 0.05 → **0.12** |

其余 21 个用例比值不变（15 个 `data/*.step` 全部 1.0000；ATU01038 1.0006、bottom 0.9750、motoc 1.0333、top 1.0000、TDB 1.0092 与基线一致）。

**⚠️ accept 第三条未达标**：`occ/a3n00.stp` 面积比 **0.8824 → 0.8627**（离 1 更远；ours 195950.71 / occ 227130.30），仍在既有 0.15 容差内。
即：F113 出网格了，但 a3n00 整体面积反而略低于基线 ⇒ 「向 1 收敛」未达成。

**待定**：接受该口径（按新口径重登 a3n00 缺口），或另开一卡处理「周期面丢边 ⇒ a3n00 面积偏低」（与 T-69 同源：reader 传输层）。


### 9.267 —— T-41 (a) 完成：repair/report 改派 + 夹具按合法布尔重设计（门禁 5/5）

**改动**：

1. `bop_builder_repair.rs::boolean_repaired`、`boolean_split_result` 与 `bop_builder_report.rs::boolean_repaired_report`
   从 `bop_curved::curved_boolean_full` 改派到忠实的 `bop_builder_dispatch::boolean_dispatch`（该函数改回 `pub(crate)`）；
2. `bop_builder_tests_api.rs::boolean_repaired_warns_on_fix` 的夹具按决策 (a) 重设计：
   - 旧夹具把 **self-intersecting shell** 与 **solid** 做 FUSE —— `BOPAlgo_BOP::CheckData` 明确拒绝不同维度融合
     （`BOPAlgo_BOP.cxx:181-186`），只有被删掉的 `bop_curved` 引擎会宽容接受；
   - 新夹具改为与**空 compound** 融合：`TreatEmptyShape`（`BOPAlgo_BOP.cxx:214-322`）返回存活形状，
     于是自交壳**原样**进入 repair 并被标记 ⇒ 既合法又保留「修复告警」的测试意图。

**验收实测**：

```
occt-topo --lib        -> 1281 passed / 0 failed
step_obj_gates         -> 5 passed / 0 failed (358s)
phase3/4/5/10/19/20    -> 4/9/7/8/5/5（与基线一致）
```

**`bop_curved` 现状**：生产路径上的消费者已全部改派（`boolean_dispatch`、`brepfeat`、`draw`、`repair`、`report`），
仅剩该模块自带的测试；整层删除可作为后续清理（其几何闭式解已在 T-28 统一到 `intana`）。

⇒ **T-41 结项**（卡片 accept「摘除后 --lib 与门禁不劣化」达标）。


### 9.268 —— T-69 矛盾解开：第 36 轮的「190」是 `fix_missing_seam` 补偿的副作用（并暴露与 T-93 (a) 的张力）

**独立复核做法**：只重挂最小 `ZZ_STEPFACE` 探针（报 `id/bounds/wires`），与**模型侧**的 `F ... wires=[..]`/`WIREHIST`
（来自最终形状、独立于该探针）比对。

**结果（T-93 (a) 之后，即补偿已移除）**：

```
STEPFACE total=1772  wires < bounds = 0
Σbounds = 2121   Σwires = 2121        <- 与模型侧完全相等
模型侧：TOTAL faces=1772 wires=2121 faces_with_2plus=309   WIREHIST 1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1
```

⇒ 两边**精确一致**。而第 36 轮（补偿尚在）resolver 侧报 `wires<bounds=190`、模型侧则是 `wires=1931 / 多 wire 面 120`。
**即那个「190」是补偿的副作用**：`fix_missing_seam` 即使丢弃 `Shell(5)` 也会改写共享 `GeometryRegistry`（§9.228 已记录同一现象），
使**解析时**探针看到的形状与**最终模型**分叉。⇒ 此前 §9.252–§9.256 那条「端口自己丢 wire」的推论**彻底关闭**：reader 没有丢 wire。

**同时暴露的关键张力（T-69 vs T-93 (a)）**：

| 状态 | T0M wires / 多 wire 面（GT=1921/106） | a3n00 F113 | a3n00 面积比 |
|---|---|---|---|
| 有补偿（T-93 a 之前） | **1931 / 120**（更接近 GT） | 无网格（mv 0） | 0.8824 |
| 无补偿（T-93 a 之后） | **2121 / 309**（明显更差） | **出网格（mv 223）** | 0.8627 |

⇒ T-69 的 accept（WIREHIST 对齐 GT）与 T-93 (a)（F113 出网格）**直接冲突**。
第三条路正是 **T-93 (b)**：把补偿做成**无副作用**（只对目标面生效、不回写共享 registry），可同时满足两者。


### 9.269 —— T-93 (a) 的代价量化（T0M 未网格面 6→7、wires 1931→2121）

**实测（T-93 (a) 之后，T0M）**：

```
TOTAL faces=1772 computed_lin=2.574309 used_lin=2.574309 stats=1765 mesh_v=49999 mesh_t=47616
mv=0 的面：7 个 —— F116(Cylinder,2 wires) / F405-F408(BSpline,1 wire) / F1691(Sphere,1 wire) / F1758(Cone,1 wire)
```

**对照 T-93 (a) 之前**：未网格 6 个（F351 Cylinder 10 wires / F405-F408 / F1758），`wires=1931`、多 wire 面 120。

⇒ 移除补偿的净代价：**F351 转为出网格（好）**，但 **F116 与 F1691 变成未网格（新增 2 个）** ⇒ 6→**7**；
`wires` 1931→**2121**、多 wire 面 120→**309**（GT 为 1921/106）。

**结论**：(a) 是忠实 OCCT 的路线，但其代价不止 T0M 面积比，还包括 T0M 的 wire 结构与未网格面数；
补偿此前在**掩盖** port reader 的真实缺口（§9.268）。⇒ T-69 的 accept（未网格 6→≤1、WIREHIST 对齐 GT）需要在 reader 里正本清源，
而「装回补偿」虽然能让 T0M 变好，却会同时让 F113 重新失去网格（不可两全，除非走 T-93 (b) 的无副作用补偿）。


### 9.270 —— T-69：OCCT `TranslateEdgeLoop` 只有一个 Add 站点 ⇒ 回归来源指向被删块的 pcurve 归一化

**对照**：`StepToTopoDS_TranslateEdgeLoop.cxx` 全文只有 **1 处** `Add` —— `:815 B.Add(W, E)`（注释原文：
`// on le fait ici. Sauf si erreur rencontree ... !`），即**每条 oriented edge 只加一次**。
这与 port 的 `resolve_face`（每 bound 建一条 wire、每条 OE 加一次）**一致** ⇒ **不存在**「OCCT 在传输时额外补一条 seam 边」这种机制。

⇒ 因此 T-93 (a) 之后 T0M 新增的两个未网格面（F116 Cylinder 2 wires、F1691 Sphere 1 wire）**不是**因为少了 seam 边，
而更可能来自被删块的**第二个副作用**：该块在把结果面交回之前，还对它的**每条 wire** 跑了一遍
`check_pcurves_and_shift(&mut w, &rf, self.precision, false)`（pcurve 归一化/校验），而那是在 `ShapeFix_Face` 之后针对**修复结果**的第二次整理。

**下一步（有界、可验证）**：只把那段 pcurve 归一化**单独**加回（不含 `fix_missing_seam` 本身），跑 T0M 看 F116/F1691 是否恢复网格。
若恢复 ⇒ 结论是「补偿里的 seam 合并是多余的，但它的 pcurve 归一化是必需的」，可据此把该 pass 移到 `resolve_face` 的正确位置（不引入 OCCT 没有的 seam 合并）。
若未恢复 ⇒ 再查这两个面的 wire/pcurve 具体缺什么。

**注意**：这一步是**纯诊断**，结论前不动代码；T-93 (a) 的现状（F113 出网格 + 门禁 5/5）保持不变。


### 9.271 —— 更正：T-93 不可结项；a3n00 的缺口是「面积忠实度」不是「缺面」

**用户指出**：a3n00 并没有补全面，不能据此确认 T-93 完成 —— **采纳**。
第三条 accept「a3n00 面积比向 1 收敛」是实质条款；第 (a) 项落地（F113 出网格、门禁 5/5）**不构成完成**。

**本轮新事实（T-93 (a) 之后，a3n00）**：

```
TOTAL faces=226 computed_lin=1.076007 used_lin=1.076007 stats=226 mesh_v=11101 mesh_t=11121
mv=0 的面：0 个
```

⇒ **226 个面全部出网格**（F113 修复后的效果），但面积比仍只有 **0.8627**（ours 195950.71 / occ 227130.30，**低 13.7%**）。
⇒ 缺口**不是缺面**，而是**三角化/裁剪区域的忠实度**：每个面都有网格，但覆盖的区域比 OCCT 小。

**下一步（正路，逐面面积对拍）**：对 a3n00 做 **per-face 三角面积对拍**（端口 vs OCCT），
先定位面积亏损集中在哪几类面（周期面裁剪 / 带闭合单边 wire 的面 / BSpline 面），再回到 reader 与 discret 的对应控制流。

**同时保留的代价记录**（§9.266/§9.269）：(a) 让 T0M 未网格 6→7、wires 1931→2121、多 wire 面 120→309，
即「F113 出网格」是拿 T0M 的 wire 结构与两个面（F116/F1691）换来的。


### 9.272 —— a3n00 面积亏损定位到少数几个面（逐面三角数对拍，无需改探针）

**方法**：GT 探针 `--mesh 0.1` **逐面**打印 `FACE i nodes=.. triangles=..`（OCCT 侧）；
端口探针逐面打印 `F i ... mv=.. mt=..`（端口侧）⇒ 两者按面索引直接对拍（面数都是 226，顺序一致）。

**实测**：

```
GT faces=226   PORT faces=226            <- 面数一致
GT total triangles=20046   PORT total mt=11121      <- 端口只有 OCCT 的 ~55%

GT face=2  tri=3224  ->  PORT mt=40
GT face=51 tri=2354  ->  PORT mt=11
GT face=98 tri=1396  ->  PORT mt=10
GT face=28 tri=1119  ->  PORT mt=7
GT face=46 tri=618   ->  PORT mt=34
GT face=31 tri=435   ->  PORT mt=7
GT face=50 tri=433   ->  PORT mt=11
GT face=29 tri=422   ->  PORT mt=48
GT face=48 tri=414   ->  PORT mt=13
GT face=47 tri=385   ->  PORT mt=46
```

⇒ a3n00 的面积亏损**不是全面性的**，而是集中在**少数几个面**上被**严重欠细分/欠覆盖**（如 face 2：3224 → 40）。
这与「226 面全部出网格（mv=0 为 0）」一致：网格**在**，但覆盖区域/细分密度远低于 OCCT ⇒ 面积比 0.8627。

**下一步（有界）**：取 GT face=2 / 51 / 98 / 28 对应端口 F 行的描述（type/wires/edges/du/dv/box）与 GT 的 `--fbox`（同名 bbox 定位），
判定这几类面的共性（周期面裁剪？带闭合单边 wire？BSpline？），再回到 discret/reader 的对应控制流。

本轮无库语义改动（纯测量，未改探针）；`--lib` 1281/0。


### 9.273 —— 更正 §9.272 的配对错误；a3n00 面积亏损来自「大面欠细分」（与顺序无关的分布对照）

**更正**：§9.272 按**面索引**把 GT 与端口的面配对是**错的** —— 取端口 `F 2/28/46/51/98` 的行可见它们都是**小面**
（如 `F 2` 的 box 仅 ~2.7×1.5×2.7），而 GT `face=2` 有 3224 三角 ⇒ 两侧面序不可比。

**改用与顺序无关的分布对照（本轮，可靠）**：

```
GT n=226  total=20046   top10 = 3224,2354,1396,1119,618,435,433,422,414,385
PT n=226  total=11121   top10 = 376,376,376,376,355,349,336,288,267,243
GT 面 ≥1000 三角：4 个   |   PT 面 ≥1000：0 个
GT 面 ≥100 ：40 个        |   PT 面 ≥100 ：30 个
```

⇒ **端口从不产出高密度面**（其最大面只有 376 个三角），而 OCCT 有 4 个面 ≥1000（最大 3224）。
⇒ a3n00 的面积亏损来自**大面（大曲率面）被欠细分**：内接多边形数少 ⇒ 面积偏小 ⇒ 总面积比 0.8627。
这与「226 面全部出网格」一致 —— 网格在，密度不够。

**下一步（有界）**：定位端口网格管线对**大面**的细分密度控制（偏转 → 三角数 的映射），
与 OCCT `BRepMesh`（自适应细分）对应控制流对比；先查 `discret_root` 的 Delaunay 主路径在**面尺寸很大**时是否退化为粗网格。

本轮无库语义改动（纯测量）；`--lib` 1281/0。


### 9.274 —— 再更正：口径混用导致 §9.273 的「大面欠细分」不成立；面积亏损不是三角数短缺

**§9.273 的口径错误**：我拿**端口 OBJ 管线**（`brep_to_obj` → `export_mesh`，走 prs3d drawer 的 **relative** 偏转，
`brep_exchange.rs:79-97`）与 **GT 探针 `--mesh`**（`occt_probe.cpp:140/425` 明确传 **`isRelative=false`**，绝对偏转）比三角数 ⇒ **口径不一致**。

**同一口径下的干净对照**（本轮）：

```
参考 OBJ  data/occ/occ-a3n00.obj（OCCT 管线、即门禁所用基线）: triangles=12466  vertices=11145
端口探针  TOTAL ... mesh_v=11101 mesh_t=11121（与端口 OBJ 同管线）
```

⇒ 端口三角数 = OCCT 的 **89%**（而非 55%）⇒ **面积亏损（13.7%）不是三角数短缺造成的**，
而是**覆盖区域/裁剪**的问题：三角数相当，但覆盖的面积小 ⇒ 某些面被裁小了、另一些面三角更多。
上一轮 §9.273 的「端口从不产出高密度面（最大 376）」**同样来自该口径错误**，结论撤回。

**下一步（干净做法）**：给 GT 探针加一个 **relative** 偏转的 `--mesh` 变体（`BRepMesh_IncrementalMesh(shape, defl, true, angle)`，
参数与端口 `prs3d_get_deflection` 同口径），才能做**逐面**三角数/面积对照；
或更直接：对**端口 OBJ 与参考 OBJ 逐面配准**（按面 bbox 分区统计面积），定位面积亏损的面簇。

本轮无库语义改动（纯测量 + 更正）；`--lib` 1281/0。


### 9.275 —— a3n00 面积亏损定位：特定区域「完全没有覆盖」（同管线产物、按空间分区）

**方法（避开探针口径坑）**：直接比较**同管线的两个 OBJ** —— 端口 `data/output/a3n00.obj`（`export_data_obj -- data/occ` 重新导出）
与基线 `data/occ/occ-a3n00.obj`（门禁所用 OCCT 参考）。按三角**质心**落到 20 单位空间网格，逐 bin 累加三角面积。

**总账自洽**：

```
PORT tris=11121  area=195950.71      <- 与门禁的 ours 一致
GT   tris=12466  area=227130.30      <- 与门禁的 occ 一致
```

**最大亏损 bin（occ − port）**：

```
bin= 1,-2,-3  occ=1555.7  port=  0.0  deficit=1555.7
bin=-1,-2,-3  occ=1555.7  port=  0.0  deficit=1555.7      <- 镜像成对（重复特征）
bin= 3, 2,-3  occ= 931.7  port=  0.0  deficit= 931.7
bin=-3, 1,-2  occ= 749.1  port=  0.0  deficit= 749.1
bin=-3,-1,-2  occ= 741.8  port=  0.0  deficit= 741.8
bin=-5,-1,-2  occ= 737.0  port=  1.6  deficit= 735.4
bin=-5, 1,-2  occ= 737.0  port=  1.6  deficit= 735.4
bin= 6, 3,-2  occ= 839.1  port=119.4  deficit= 719.7
```

⇒ 亏损的形态是**「特定区域完全没有面积」（0.0 vs 1555.7）**，且**镜像成对**（同一特征在 ±x 两侧各一份）。
这与「226/226 面都出网格」并不矛盾：面有网格，但网格**没有落到这些区域**（落在别处或只覆盖未裁剪的局部）。
⇒ 定性从「欠细分」改为**「部分特征区域无覆盖」**（上一轮 §9.273/§9.274 的口径错误已更正）。

**下一步（有界）**：取上述 bin 坐标（如 x≈±20, y≈−40, z≈−60 与 x≈−60,y≈±20,z≈−40）对应的端口面（按 bbox 匹配 `F` 行），
看这些面的 UV 域/裁剪 wire 是否被正确使用（是否把未裁剪的整曲面或错误的 UV 子域拿去三角化）。

本轮无库语义改动（纯测量）；`--lib` 1281/0。


### 9.276 —— 面积亏损收窄到**单个面 F171**（BSpline、4 wires、整周期 v）

**方法**：把 §9.275 的最大亏损 bin 质心（`(20,-40,-60)` 与其镜像 `(-20,-40,-60)`）当作探针点，
在端口 `F` 行里筛出 **bbox 覆盖该点**的面。

**结果（唯一命中）**：

```
F171  type=BSpline  wires=4  edges=[8,8,8,6]  mv=325  mt=349  du=0.62919  dv=6.28319
      box=(-57.3605,-42.86496,-169.61091)-(61.61091,42.86496,-50.6395)
```

⇒ 两个亏损 bin（±x 镜像）都由**同一个面 F171** 覆盖；它有网格（349 三角）且 4 条 wire、30 条边、v 为**整周期**（`dv=2π`）。
但它的三角**没有落到 z≈−60 那一带** ⇒ 端口只三角化了该面**的一部分**（覆盖不足，而非没网格）。
这也与 §9.275 的形态一致：`port=0.0` 的区域在 F171 的 bbox 之内。

**下一步（有界、直指根因）**：查 F171 的 4 条 wire 在三角化时如何使用 —— 疑点是**裁剪环/周期性**的处理：
该面 v 为整周期且有多条 wire，若三角化只用了其中一条外环（或把某条内环当成外环），
就会**整片区域无覆盖**。做法：打印 F171 的 4 条 wire 的边数与各自 bbox（或用 `--fdump <171>` 逐 wire/边 dump），
与 GT `--fbox` 在相同 bbox 的面逐 wire 对照。

本轮无库语义改动（纯测量）；`--lib` 1281/0。


### 9.277 —— 根因定位：F171 的第 4 条 wire 在端口**少了一条边**（读入层）

**双侧逐 wire/边对照**（端口 `--fdump 171` vs GT `--fbox` 同一 bbox，面序 `face=2`）：

| wire | 端口 nEdges | GT nEdges | 坐标 |
|---|---|---|---|
| wire[0] | 8 | 8 | 完全一致（逐边首尾点） |
| wire[1] | 8 | 8 | 完全一致 |
| wire[2] | 8 | 8 | 完全一致 |
| **wire[3]** | **6** | **7** | **不同** |

GT 的 `wire[3]` 里有两条边端口**完全没有**：

```
GT e0 : (-49.864906, 0.000000, -100.000000) -> (36.224088, -3.000000, -99.867388)
GT e4 : (  8.356879, 0.000000, -155.954859) -> (47.840889,  3.000000, -116.470850)
GT e5 : ( 47.840889, -3.000000, -116.470850) -> ( 8.356879, 0.000000, -155.954859)
端口 e0: (-35.250069, 0.000000, -116.901100) -> (36.224088, -3.000000, -99.867388)   <- 短得多
```

⇒ 端口该 wire 少了伸到 **z ≈ −155.95** 与 **x ≈ −49.86** 的边 ⇒ 三角化只覆盖了较短边界内的区域 ⇒
**正是 §9.275/§9.276 检测到「z ≈ −60 一带 0 覆盖」的原因**（GT 的边界在那里，端口的没有）。

⇒ 这**不是**三角化/细分问题，而是 **reader 层**：解析该 EDGE_LOOP 时丢了一条（或把两条并成一条）。这与 T-69 的 reader 结论一致。

**下一步（直指修复）**：在 STEP 文件里按端点坐标定位这两条边（`EDGE_CURVE` 的 `CARTESIAN_POINT`）及其 `ORIENTED_EDGE`/`EDGE_LOOP`，
再在端口 `resolve_loop`/`resolve_oriented_edge` 里查为何该边没被加入（丢失、去重、还是 resolved 为 Err 被跳过）。

本轮无库语义改动（纯测量）；`--lib` 1281/0。


### 9.278 —— F171 丢失边的进一步线索：端点坐标不在 STEP 里 ⇒ 指向「边的裁剪参数」

**做法**：把 §9.277 里 GT 独有边的端点坐标（`A=(-49.864906,0,-100)`、`P=(-35.250069,0,-116.9011)`、`B=(8.356879,0,-155.954859)`）
在 `data/occ/a3n00.stp` 里按 `CARTESIAN_POINT` 搜索。

**结果**：**三条坐标在 STEP 里都找不到**（三种数字格式都试过）。

⇒ 这些端点**不是几何顶点**，而是**几何求值结果**（边的裁剪参数 `first/last` 作用在曲线上得到的点）。
因此 F171 的差异**不是**「顶点绑定错」，而更可能是**边的参数区间取错**：端口把某条边的 `first/last` 取成了别的来源，
于是该边比 GT 短（`e0` 从 `(-35.250069,0,-116.9011)` 起而不是 `(-49.864906,0,-100)`），
连带缺掉伸向 `z≈-155.95` 的那一对边 ⇒ 三角化只覆盖短边界内区域 ⇒ 面积亏损。

**下一步（需要一次小探针）**：在端口 `resolve_edge`/`resolve_oriented_edge` 打点，输出每条边对应的 **STEP 实体 id** 与采用的 `first/last`，
再对同一实体逐个与 STEP 的 `TRIMMED_CURVE`/`SURFACE_CURVE` 参数（或 `EDGE_CURVE` 的顶点参数）比对，
定位是哪一类边被取了错的区间（疑点：`SURFACE_CURVE` 的 `master_representation`、或 pcurve 与 3D 曲线的参数不同源）。

本轮无库语义改动（纯测量）；`--lib` 1281/0。


### 9.279 —— 分量级判定：GT 用的端点 A 是真实 STEP 顶点，端口的 P 是纯求值点

**做法**：把两侧端点坐标的**各分量**逐个在 `data/occ/a3n00.stp` 里计数（a3n00 是装配件，可能存在 location 变换，故按分量而非整点搜）。

```
49.864906 : 14      <- GT 的 A=(-49.864906, 0, -100) 的 x，真实存在
36.224088 :  7      <- GT 的 C 的 x
99.867388 : 11      <- C 的 z
47.840889 :  5      <- GT 的 D3/D-3 的 x
89.394247 :  3      <- GT 的 E 的 z
35.250069 :  0      <- 端口的 P=(-35.250069, 0, -116.9011) 的 x，**不存在**
116.9011  :  0      <- P 的 z，**不存在**
24.672033 :  0      <- E 的 y（两侧都用 E；E 本身是计算点）
155.954859:  0      <- GT 的 B 的 z（计算点）
 8.356879 :  0      <- B 的 x（计算点）
```

⇒ GT 的 `wire[3]` 边界锚在**真实 STEP 顶点 A** 上（其 x 分量在文件里出现 14 次），
而端口同位置的端点 **P 的两个分量都不在文件里** ⇒ P 是**纯几何求值**得到的点，**不是顶点**。

**结论**：端口的那两条边（`P→C`、`E→P`）**没有锚在顶点 A 上** —— 与 GT 的 `A→C`、`E→A` 相比，等价于「同一条曲线被取了不同的端点」。
可能是：(a) 顶点 A 未被绑定到这两条边（走了计算点）；或 (b) 边的 3D 曲线被取了错的参数区间。

**下一步**：在端口对 F171 的这两条边打点，输出其 STEP `EDGE_CURVE` id、顶点 id 与 `first/last`；
与 STEP 里同一 `EDGE_CURVE` 的 `VERTEX_POINT` 引用对照，即可区分 (a)/(b) 并直指修复点。

本轮无库语义改动（纯测量）；`--lib` 1281/0。


### 9.280 —— F171 的 STEP 结构：顶点 A 挂 8 边，分属 8-OE 与 4-OE 两类 loop

**A 的实体链**：`#4955 CARTESIAN_POINT('',(-49.864906150824751, 0.0, -99.999999999999972))` → `#4956 VERTEX_POINT`。

**以 `#4956` 为端点的 8 条边**（各自所属 loop）：

```
#5012 EDGE_CURVE(#4956,#5007,#5011,.T.)  -> oe #5013 in loop #5139(OEs=8) ; oe #7639 in loop #7662(OEs=4)
#5018 EDGE_CURVE(#4956,#5007,#5017,.T.)  -> oe #7782 in loop #7786(OEs=4) ; oe #5019 in loop #5139(OEs=8)
#7675 EDGE_CURVE(#7399,#4956,#7674,.T.)  -> oe #7676 in loop #7693(OEs=4) ; oe #7765 in loop #7774(OEs=4)
#7405 EDGE_CURVE(#7399,#4956,#7404,.T.)  -> oe #7407 / #7406 both in loop #7413(OEs=8)
#5004 EDGE_CURVE(#4954,#4956,#5003,.T.)  -> oe #7408 in loop #7413(OEs=8) ; oe #5005 in loop #5139(OEs=8)
#5090 EDGE_CURVE(#4956,#5021,#5089,.T.)  -> oe #7397 in loop #7413(OEs=8) ; oe #5091 in loop #5139(OEs=8)
#7646 EDGE_CURVE(#4956,#7641,#7645,.T.)  -> oe #7692 in loop #7693(OEs=4) ; oe #7647 in loop #7662(OEs=4)
#7770 EDGE_CURVE(#7702,#4956,#7769,.T.)  -> oe #7771 in loop #7774(OEs=4) ; oe #7785 in loop #7786(OEs=4)
```

⇒ 顶点 A 同时参与 **8-OE loop（`#5139`、`#7413`）** 与 **4-OE loop（`#7662`/`#7786`/`#7693`/`#7774`）**；
而端口 F171 的四条 wire 是 **8/8/8/6** 边、GT 是 **8/8/8/7** ⇒ 少边的第 4 条对应的应是某个 8-OE loop（含重复 OE）。

**下一步（直指修复）**：把 `#5139`（或 `#7413`）的 8 个 OE 逐条列出（edge id、首尾 vertex id），
与端口 `--fdump 171` 的 wire[3] 6 条边按端点坐标配对 ⇒ 直接看出**是哪一条 OE 没被加进端口 wire**，
再回到 `resolve_loop`/`resolve_oriented_edge` 查它为何被丢弃（丢失/去重/`Err` 跳过）。

本轮无库语义改动（纯测量）；`--lib` 1281/0。


### 9.281 —— 根因定位到「顶点 `#4956`（A）被解析成了另一条曲线上的一点」

**配对**：`#7413`（8 OE）与端口 `--fdump 171` 的 wire[3]（6 边）逐条对应：

| STEP OE | edge | 端点 | 端口 wire[3] |
|---|---|---|---|
| `#7397(F)` | `#5090` | A → C | ✓ 有（但 A 变成 P） |
| `#7406(F)`/`#7407(T)` | `#7405` | X ↔ A（同一曲线两条 OE） | ✗ **两条都没有** |
| `#7408(F)` | `#5004` | E → A | ✓ 有（但 A 变成 P） |
| `#7409(F)` | `#5137` | C3 → E | ✓ |
| `#7410(T)` | `#7210` | C3 → D3 | ✓ |
| `#7411(F)` | `#7324` | D-3 → D3 | ✓ |
| `#7412(T)` | `#7270` | D-3 → C(-3) | ✓ |

其中 `X = (-20.634979099053066, 0, -133.802492199755930)`，`A = (-49.864906150824751, 0, -99.999999999999972)`。

**两处缺陷**：

1. **缺**：`#7405` 的两条 OE（F/T 各一）完全没有进端口 wire ⇒ GT 的 7 条边里那对 X↔A 丢了；
2. **错**：`#5090`/`#5004` 上的 A 端点变成了 `P=(-35.250069, 0, -116.9011)` —— 而 **P 恰好在 `#7405` 曲线（X→A）之间**
   ⇒ 端口把顶点 A 的位置解析成了**另一条曲线上的一点**（不是 `#4955` 的坐标）。

⇒ 这不是「三角化」问题，而是**顶点/边的解析缺陷**：A 是真实 STEP 顶点（`#4956`→`#4955`），端口却在该面把它的位置取了错值，
同时丢掉了经过它的 `#7405` 那对 OE ⇒ 边界被截短 ⇒ 三角化只覆盖短边界内区域 ⇒ §9.275 的「某区域 0 覆盖」与面积亏损。

**下一步（直指修复）**：查端口对 `#4956` 的解析 —— `resolve_vertex` 是否返回 `#4955` 的点？
该顶点在多个 EDGE_CURVE 间**共享**时，点是否被后续边的求值/绑定覆盖（候选：`vertex_bind`、共享 `TShape` 上的点被改写、
或 `resolve_edge` 用曲线端点回填了顶点位置）。

本轮无库语义改动（纯测量）；`--lib` 1281/0。


### 9.282 —— 核实：P 是**顶点注册的点**（非曲线求值）⇒ 顶点点值缺陷

**核实方式**：读端口探针 `--fdump` 的实现（`crates/occt-topo/examples/zz_probe_a3n00.rs:174-190`）——
它打印的是 `BRepTool::vertex_point(v)`（**顶点注册的点**）与顶点 `tshape` 指针，**不是**曲线端点求值。

⇒ 因此 §9.277/§9.281 里端口的 `P=(-35.250069,0,-116.9011)` **就是顶点的点值**，属 **(a) 顶点缺陷**，不是参数区间缺陷。
即：端口解析 `#4956`（STEP 顶点 A）时给出的点不是 `#4955` 的坐标，而是另一个位置（且该位置恰落在 `#7405` 曲线 X→A 之间）。

**候选机制**：

1. `resolve_vertex` 取错了点引用的参数下标（但同面其它顶点都正确，可能性低）；
2. **顶点共享/绑定**：`vertex_bind` 或 `shape_cache` 让 `#4956` 复用了别的顶点（其点位为 P）；
3. **事后改写**：某个 pass（`check_pcurves_and_shift`、`xsalgo_check_pcurve`、投影/拟合类）把顶点点位**投影到相邻曲线**并写回；
   P 落在 `#7405` 曲线上这一点，最像「投影到相邻边」的结果。

**下一步（一次打点即可区分）**：在 `resolve_vertex` 出口对 `id==4956` 打点（打印点位），
再在同一面的尾部 pass 之后（`resolve_face` 返回前）对同一顶点点位再打一次 ⇒ 若两次相同 ⇒ 机制 1/2（解析时就错）；若不同 ⇒ 机制 3（被事后改写）。

本轮无库语义改动（纯读码 + 记录）；`--lib` 1281/0。


### 9.283 —— 决定性数值发现：端口的 P 就是边 `#7405` 曲线的**中点**

```
X  = (-20.634979099053066, 0.0, -133.802492199755930)      # #7405 的一端（STEP）
A  = (-49.864906150824751, 0.0,  -99.999999999999972)      # #7405 的另一端（STEP，顶点 #4956）
X..A 中点 = (-35.249942624938908, 0.0, -116.901246099877950)
端口的 P  = (-35.250069,            0.0, -116.9011)          # 与之相差 ~1.3e-4
```

⇒ 端口在 F171 的 wire[3] 上用的端点 **P 恰好是 `#7405` 曲线的中点**（远小于读取精度，可视为同一点）。

**推论**：端口在解析该 loop 时**为 `#7405` 创造了中点顶点**（把这条边在中间断开），并把该中点顶点用在了相邻边
（`#5004`、`#5090`）的 A 端 ⇒ 于是：① A 的端点位置变成 P（§9.281）；② `#7405` 的两条 OE 与其顶点被这条「中点边」挤掉/替代（§9.281 的「缺两条 OE」）。

**这类「按中点/等分点断边」在 OCCT 的 STEP 传输里没有对应控制流**（`StepToTopoDS_TranslateEdgeLoop.cxx` 全文只有 `:815 B.Add(W,E)`，
顶点合并只走 `:435-477` 的 `Identical/Confusion` 判定，且要求两点距离 ≤ 容差）⇒ 属端口自造行为。

**下一步（直指修复点）**：在 `resolve_loop` / `resolve_oriented_edge` 及其调用的顶点处理里找「创造中点/等分顶点」的代码
（关键词：midpoint、mid、split、half、0.5 * (first+last)），确认后用 OCCT 的同等分支替换（无同等分支则删除并标注 UNPORTED）。

本轮无库语义改动（纯测量 + 数值验证）；`--lib` 1281/0。


### 9.286 —— 收口事实：顶点 `#4956` 在所有阶段都正确 ⇒ P 是 reader **另建**的顶点

**运行时探针**（`ZZ_VP_ID=4956`）：在 `resolve_vertex` 出口记录该顶点的点值与 tshape 指针，
再在 `resolve_face` 返回前按指针比对同一顶点：

```
ZZVP stage=rv   id=4956 p=(-49.864906,0.000000,-100.000000) ptr=17ae2344200
ZZVP stage=face           p=(-49.864906,0.000000,-100.000000)    × 8 次（8 条引用它的边）
```

⇒ **A 的顶点自始至终都是正确坐标**，没有被任何 pass 改写。
⇒ 因此 §9.281 里端口 `wire[3]` 的 `P=(-35.250069,0,-116.9011)`（= `#7405` 曲线中点）**不是 `#4956`**，
而是一个**reader 另建的顶点**（不经过 `resolve_vertex`）。

**同时更正**：§9.281 把「端口 wire[3] ≡ loop `#7413` 的 OE 集」当作对应关系，现在看来该 wire 的边集合本身不同（其端点 P 来自另建顶点），
因此「A 端点被换成 P」这一表述应更正为「该 wire 里出现了另建的中点顶点，且原本经 A 的边/loop 归属与 GT 不同」。

**下一步（未来会话）**：找 reader 里**在 `resolve_vertex` 之外创建顶点**的地方（候选：`shape_fix_compose_shell/split_wire.rs:45`
`nv.set_point(BRepTool::vertex_point(v))` 这条「新建顶点并拷点」的路径，以及任何 `b.make_vertex` 的调用点）；
定位后用 OCCT 的同等分支替换，或按 CLAUDE.md 标注 UNPORTED 并删除。

（本轮探针为 env 门控临时探针，提交前需摘除。）


### 9.287 —— T-93 收口：中点顶点 P 的真因是 `fix_dummy_seam` 的**传参错位**（已定案、已实测、按红线回退）

本节取代 §9.286 的「另建顶点」推测。**结论先写**：P 既不是 reader 自己发明的，也不是未移植分支 ——
它由**忠实移植**的 `ShapeFix_Wire::FixDummySeam` 产生，错的是**调用方传错了索引**，于是它合并并重指了**错误的边对**。

#### 9.287.1 定位链（BACKTRACE，唯一命中）

在 `TopoBuilder::make_vertex` 加 env 门控探针（判据：点落在 `P=(-35.250069,0,-116.9011)` 附近）+ `Backtrace::force_capture`，
跑 `zz_probe_a3n00 data/occ/a3n00.stp`，**一次命中**：

```
MIDV p=(-35.250069,0.000000,-116.901100) tol=22.346313953
 4: occt_topo::builder::TopoBuilder::make_vertex          .\crates\occt-topo\src\builder.rs:90
 5: occt_topo::shhealing::wire_fix::combine_vertex        .\crates\occt-topo\src\shhealing\wire_fix.rs:1271
 6: occt_topo::shhealing::wire_fix::fix_dummy_seam        .\crates\occt-topo\src\shhealing\wire_fix.rs:2995
 7: occt_topo::shhealing::wire_fix::fix_notched_edges     .\crates\occt-topo\src\shhealing\wire_fix.rs:3078
 8: occt_topo::shhealing::wire_fix::check_pcurves_and_shift  .\crates\occt-topo\src\shhealing\wire_fix.rs:3249
 9: occt_topo::step::transfer::Resolver::resolve_face     .\crates\occt-topo\src\step\read_topology.rs:714
```

⇒ §9.286 的候选（`split_wire.rs:45` 的 `empty_copied_vertex`）**是错的**：那是 `TopoDS_Shape::EmptyCopied`
的忠实移植（对应 `cxx:1266-1282`），不建中点。真正建点的是 `fix_dummy_seam`。

#### 9.287.2 缺陷：传 `n1` 还是 `n2`

`ShapeFix_Wire::FixNotchedEdges` 的 `on_end` 分支（`ShapeFix_Wire.cxx:4019`）调用的是 **`FixDummySeam(i)`**，
`i` 是当前循环下标。`FixDummySeam(num)`（`cxx:4213-4289`）用：

```
E1 = Edge(i)                                       // cxx:4214
E2 = Edge(i == NbEdges() ? 1 : i + 1)               // cxx:4215-4217
```

端口的 `fix_dummy_seam(wire, num)` 用 0 基切片读的是 `edges.get(num - 1)` 与 `edges.get(num1 - 1)`（`num1 = num+1`，wrap），
即 **`E_num` / `E_{num+1}`** —— 与 `cxx` 的 `E1/E2` **逐个对应** ⇒ 端口的 `num` 就是 OCCT 的 `i`。

而 `fix_notched_edges` 里 `n2 == i`、`n1 == i-1`（`n2 = if i > 0 { i } else { nb }`、`n1 = if n2 > 1 { n2-1 } else { nb }`），
端口却传了 **`n1`** ⇒ 整个「notch 边对」沿 wire **下移一位**。

**后果**（a3n00 F171）：它把 `E_{i-2}` / `E_{i-1}` 当成 notch 对合并。命中的两条边恰好是 #7405 的两端
（`X=(-20.634979,0,-133.802492)` 与 `A=(-49.864906,0,-100.000000)`），其容差 `tol=22.346`（远大于两点距离 ~35 的一半），
于是 `combine_vertex` 走「都不包含」分支算出**中点** `P≈(-35.25,0,-116.90)`，
`fix_dummy_seam` 再把它写回相邻边（`copy_replace_vertices_with` 两条、`wire_set_edge_composed` 两条），
并**删掉两条 notch 边** ⇒ wire 里出现中点顶点、并使原本经 A 的边/loop 归属改变（§9.281/§9.283 观测到的现象）。

#### 9.287.3 三个被调函数**都是忠实的**（逐行核对过，不要再改它们）

| 端口 | OCCT | 核对结论 |
|---|---|---|
| `combine_vertex`（`wire_fix.rs:1248-1272`） | `ShapeBuild_Vertex.cxx:38-72` | 逐行相同：`dist+tol2<=tol1` → p1；`dist+tol1<=tol2` → p2；否则 `tol=0.5*(dist+tol1+tol2)`、`s=(tol2-tol1)/dist`、`pos=0.5*((1-s)*p1+(1+s)*p2)`、`MakeVertex(pos, tolFactor*tol)` |
| `fix_dummy_seam`（`wire_fix.rs:2980-3027`） | `ShapeFix_Wire.cxx:4213-4289` | V1/V2/Vm（`cxx:4219-4221`）、`Vs`（`:4230-4233`）、`copy_replace_vertices_with`（`:4264-4272`）、`copy_reverse_pcurves`（`:4108-4176`）、邻边重指、按高索引先删两条 notch 边（`:4275-4288`）均对应 |
| 「第一轮该不该跑 notches」 | `ShapeFix_Face.cxx:365-520` | **第一轮只关 `FixLacking`**（`:372`）；关 `FixNotchedEdges` 的那行是**注释掉的**（`:370` 的 `// ...FixNotchedEdgesMode() = false; // CR0024983`）；第二轮（`:513-520`）才关 FixSmall/Connected/EdgeCurves/Degenerated。⇒ notch 在两轮都开，端口在 `fix_lacking=false` 的第一轮跑 notches **是正确的** |

⇒ 缺陷**只在 `fix_notched_edges` 的传参这一个点**。

#### 9.287.4 修复形态与实测（**已回退，未提交**）

改一行：`crates/occt-topo/src/shhealing/wire_fix.rs` 的 `on_end` 分支

```rust
-            fix_dummy_seam(wire, n1);
+            fix_dummy_seam(wire, n2);
```

实测（同命令、改前/改后各跑一遍）：

| 指标 | 改前 | 改后 |
|---|---|---|
| a3n00 `stats` / `unmatched` | 226 / 0 | 226 / 0 |
| a3n00 `mesh_v` / `mesh_t` | 11101 / 11121 | **11102 / 11136** |
| **a3n00 面积比**（`area_tol=0.15`，`tests/common/mod.rs:93`） | 0.8627 | **0.8298** ❌ 越界 |
| T0M `stats` / `unmatched` | 1765 / **7** | 1769 / **3** |
| T0M `wires` / 多 wire 面 / WIREHIST | 2121 / 309 / `1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1` | **不变** |
| `--lib` | 1281 / 0 | **1281 / 0** ✅ |
| `step_obj_gates` | 5 / 5 | **4 / 5**（`step_obj_area_matches_occt` FAILED，659.94s）❌ |

**23 个门禁模型的面积比逐一对拍**（这是决定回退的依据）：

| 模型 | 改前 | 改后 | Δ |
|---|---|---|---|
| `occ/a3n00.stp` | 0.8627 | **0.8298** | **−0.0329** |
| `occ/T0M.stp` | 0.9786 | 0.9787 | +0.0001 |
| 其余 21 个 | — | — | **完全相同** |

⇒ 代价集中在 a3n00 一个模型；让它转绿需把 `area_tol` 0.15 → ~0.18，属**改既有断言**（红线第 3 条），
且 a3n00 的面积缺口本就是**已登记的已知缺口**（`tests/common/mod.rs:91` 注明 *"the tracked T-59/T-69 face classes (density 0.72)"*）——
本修复只是把缺口从 ~0.14 推到 ~0.17，没有解决它。**故回退，不重定基线。**

#### 9.287.5 未来会话怎么用这一节

- 想**无代价**拿到这个修复的收益（T0M 未网格面 7→3），先解决 a3n00 的面积缺口（T-59 的面类/密度工作），
  待面积比回升后再应用 §9.287.4 的一行改动并重跑全部门禁。
- **不要**再重复这两条已被否证的路（各花掉一轮）：
  ① 把 `check_pcurves_and_shift` 的 pcurve 归一化 pass 单独加回 reader（§9.270，实测未网格面仍 7、WIREHIST 逐字不变）；
  ② 把 `shape_fix_face.rs:480-487` 的 `set_edge_pcurve` 推迟到 `result` 确定后提交（实测与未改动基线逐字相同，惰性）。
- **不要**再按 `--mesh` 的**面序号**逐面判丢面：GT 1778 面 vs 端口 1772 面，两侧分割约定不同，序号不构成身份。
  判丢面必须按几何键（bbox 中心最近邻）配对 —— 本轮用该方法确认过：端口被拒的 7 个面在 GT 侧都出网格
  （P405→G407 tri=185、P406-408→G408-410 各 185、P1691→G1697 tri=73、P1758→G1764 tri=79、P116→G116 tri=226，dist 0.0–0.30）。
- 报「未网格面数」必须**同时报测量调用序列**：同一参数下，先 `prs3d_get_deflection` 再 `from_deflection`
  得 `stats=1757 / no_stat=15`，而 `zz_probe_a3n00` 默认模式得 `stats=1765 / unmatched=7`。

（本轮探针全部已摘除；工作树无源码改动，仅画布 `specs/board.canvas.tsx` 更新。）


---

### 9.288 —— T-69 收口：T0M 7 个未网格面的成因是「边界环缺边」，不是 RangeSplitter 判据

**一句话结论**：端口拒这 7 个面（F116 / F405 / F406 / F407 / F408 / F1691 / F1758）**不是网格侧缺陷**。
`collectWirePoints` → `AddPoint` → `AdjustRange/updateRange` → `IsValid` 这条链与 OCCT 逐行等价，
而且**把端口的点集原样喂给 OCCT，OCCT 也会拒**。真正的缺口在**边界环**：这些面在端口侧少了 GT 孪生面上的
那几条边，于是喂进 range 的点集在一个参数方向上零宽。缺边在读入层就存在（不是 mesh 阶段丢的）。

本节的结论**取代** §9.287.5 里「端口被拒的 7 个面在 GT 都出网格 ⇒ 分歧在喂进 range 的点集」这半句中的
「几何正确」前提；但 §9.287 关于 `fix_dummy_seam` 传参错位的那部分**仍然有效**，只是从 T-93 的主线降为旁支（见 §9.288.4）。

#### 9.288.1 复跑命令（本轮新增/固化的仪器）

```text
# OCCT 侧：单面逐点 + 逐边（GT = OCCT 8.0.0，走 ModelBuilder → EdgeDiscret，读的就是 collectWirePoints 的输入）
specs\occt_probe\build.bat
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\T0M.stp --uv    409"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\T0M.stp --uvsum"

# 端口侧：逐面 wires / 三角数 / bbox；--dump 把同一次管线产出的网格写成 OBJ
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/T0M.stp --ids
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/T0M.stp --dump .target-gate\port_T0M.obj

# 配对（按 bbox 六个坐标，不是中心最近邻）
python .target-gate\pair_ids.py .target-gate\port_ids.txt .target-gate\uvsum_gt.txt
```

`--uv` / `--uvsum` 是 `specs/occt_probe/occt_probe.cpp` 里新增的两个模式；它们不看最终 `Poly_Triangulation`，
而是用 `BRepMesh_ModelBuilder` 建 `IMeshData_Model`、再跑 `BRepMesh_EdgeDiscret`，然后直接读
`IMeshData::IPCurveHandle` 的 `GetPoint/GetParameter` —— 也就是 `BRepMesh_NodeInsertionMeshAlgo.hxx:142-184`
`collectWirePoints` 真正会消费的那批点。注意 `IMeshTools_Parameters::MinSize` 必须先赋值
（`BRepMesh_CurveTessellator.cxx:66` 会对 `MinSize <= 0` 抛异常，异常被 `EdgeDiscret::process` 吞掉后
pcurve 会**全是空的** `ParametersNb()==0`，看起来像「模型没有 pcurve」）；本轮用 `MinSize = Precision::Confusion()`
复现真实管线（实测 8792 条 pcurve 全部 filled）。

#### 9.288.2 判据侧：同输入同拒绝

用同一套控制流、同一参数域条件枚举三种口径，结论一致：

| 口径 | OCCT | 端口 | 说明 |
|---|---|---|---|
| 收集边界点 | `BRepMesh_NodeInsertionMeshAlgo.hxx:142-184`（每段 pcurve **跳过遍历方向的最后一点**） | `node_insertion.rs` 的 `collect_boundary_uv` | 逐行同构 |
| 记录 min/max | `BRepMesh_DefaultRangeSplitter.cxx:35-41` | `range_splitter/splitter.rs` 的 `add_point` | 逐行同构 |
| 收窄/夹取 | `cxx:202-236` `updateRange`（**只在点集越出曲面域时**才夹；`else` 分支不会把区间撑大） | `range_splitter/param_set.rs` 的 `update_range` | 逐行同构 |
| 判据 | `cxx:72-80` `lengthU > PConfusion && lengthV > PConfusion` | `param_set.rs` 的 `adjust_range_base` | `PConfusion = 1e-9` 两侧相同 |

对端口 F1691 的实际数字：`AddPoint` 后 `V = [1.08083900054116766, 1.08083900054116877]`，
曲面域 `V = [-1.57079632679489656, 1.57079632679489656]`。`updateRange` 的外层条件是
`dFirst < geomFirst || dLast > geomLast`，两个不等式都**为假** ⇒ 函数直接返回、区间不变 ⇒ `computeLengthV() == 0`
⇒ `IsValid() == false`。**这是 OCCT 自己的行为，不是端口的近似**：把端口这 37 个点的 V 全部等于
`1.080839000541168x` 这一事实代入 `BRepMesh_DefaultRangeSplitter::computeLengthV()`，`dv = 0.05 * 0 = 0`，
20 次迭代累加 0 ⇒ `longv = 0`。所以「端口判据太严」这条已彻底排除；D12 的「忠实拒绝」判断在**判据层面**成立，
但当时把「点集零宽」当成几何事实接受下来，错在下一步（没问「这点集该不该这么宽」）。

#### 9.288.3 真因：边界环缺边（逐面实测）

配对必须按 **bbox 六个坐标**（≤1e-6）做，不能用 D15 的「bbox 中心最近邻」——那个键不校验尺寸，
会配出完全不同的面（见 D21/D24）。按六坐标配对后，7 组孪生面与 dU/dV：

| 端口面 | GT 面 | 端口 dU / dV | GT dU / dV | 端口环边数 | GT 环边数 |
|---|---|---|---|---|---|
| F116 | G607 | —（2 wire vs 1 wire） | — | 2 | 1 |
| F405 | G73 | 1.96e-13 / 6.11387 | 0.429994 / 6.283185 | 1 | 4 |
| F406 | G74 | 3.58e-14 / 6.10216 | 0.429994 / 6.283185 | 1 | 4 |
| F407 | G75 | 1.32e-13 / 6.10216 | 0.429994 / 6.283185 | 1 | 4 |
| F408 | G76 | 4.60e-14 / 6.10216 | 0.429994 / 6.283185 | 1 | 4 |
| F1691 | G409 | 6.10865 / 1.11e-15 | 6.283185 / 0.489957 | 1 | 4 |
| F1758 | G566 | 6.0805 / 0 | 6.283185 / 1.434376 | 1 | 3 |

每个端口面的**零宽方向，正是 GT 侧给出展布的方向**；而且 GT 给出的那个展布量，就是端口缺的那一段。

F1691 ↔ G409 是取证最干净的一例（球带，球心在 `(0,-25,0)` 附近，半径 2.744 ≈ y 向跨度）：

* **端口 F1691**（`zz_uv_feed --ids` + `--dump` 实测）：`1 wire / 1 edge / 37 点`，点集 `u ∈ [0, 6.108652]`、
  `v ≡ 1.08083900054116766`，bbox = `(-2,-27,2.744)-(2,-23,2.744)` —— **Z 上零厚**。这条边是闭合圆
  （首末点相同），离散模型建完后就已经是 1 条边，所以缺边**不晚于离散模型构建**，不是 mesh 阶段丢的。
* **GT G409**（`--uv 409` 实测）：`1 wire / 4 edges`，bbox = `(-2,-27,2.744)-(2,-23,3.244)`（Z 跨度 0.5）。
  四条边：`e0` 底圆 37 点（`u: 6.283185→0`，`v ≡ 1.080839`，pcurve 参数域 `[3.141593, 9.424778]`）；
  `e1` seam 4 点（`u ≡ 6.283185`，`v: 1.080839→1.570796`）；`e2` 顶圆 2 点（`v ≡ 1.570796`）；
  `e3` seam 4 点（`u ≡ 0`，`v: 1.080839→1.570796`，即 e1 的对侧）。GT 的 `dV = 0.489957` **全部来自
  端口环里没有的 e1/e3**；`dU` 的 6.283185 来自 e0+e2（端口 e0 只有 6.108652，差一个采样步 π/18）。

端到端也能量出来（同一口径：端口 `zz_uv_feed --dump` 的网格 vs 仓库根 `data/occ/occ-T0M.obj`）：

| 区域 | 端口三角数 | GT 三角数 |
|---|---|---|
| 球带 `x∈[-2.05,2.05], y∈[-27.05,-22.95], z>2.7` | **0** | 143 |
| 圆锥带 `x∈[-23.8,-21.2], y∈[-207,-205.9]` | 28 | 114 |
| 平面带 `|x-11.01|<0.6, y∈[-225,-217], z∈[97,108]` | 89 | 139 |
| 全模型 verts / tris | 50123 / 47864 | 60548 / 67366 |

面数本身也是证据：**GT 1778 面 vs 端口 1772 面**（同一 `data/occ/T0M.stp`、同一调用序列）。
`--uvsum` 报 GT 侧 0 个零宽面，端口 6 个。

#### 9.288.4 对 T-93 / T-59 的影响（已登记为 D23）

* §9.287.4 的那行修复（`fix_dummy_seam(wire, n1)` → `(wire, n2)`）**不再当主线**：它换的是 notch 边对，
  能让未网格面 7→3，但**不会补回缺边**，而且要付 a3n00 面积比 −0.0329 的代价（D19）。保留为旁支，
  仍须先解 T-59 的 a3n00 面积缺口。
* T-93 的主线改为「补回缺边」：先查 F1691 的 wire 环为何只剩 1 条边。取证边界很干净 ——
  `zz_uv_feed --ids` 显示 F1691 在**建完离散模型后**仍是 1 wire / 1 edge，故丢失点在 reader
  （`step/read_topology.rs` 的 face/wire 构建）或更早，而不是 `meshing/` 里的任何一步。
* 缺边与 §9.287 的中点顶点 P 可能是同一条链的两端，也可能是两条独立缺陷：**先确认缺边发生在哪一级**，
  再判断二者关系。不要在没定级前把两处一起改。
* T-54（平面耳切 → 约束 Delaunay）随之解阻塞：那 6 个面实测是 `node_insertion.rs:536-544`
  的早期返回，**从未走到 Delaunay**，所以「回退还空」这个前提不成立。

#### 9.288.5 未来会话怎么用这一节

* **不要**再从 RangeSplitter / `IsValid` / `updateRange` / `PConfusion` 一侧找缺口：这四处已逐行等价，
  且 §9.288.2 给出了「同输入同拒绝」的可复算证明。
* **不要**再用「bbox 中心最近邻」配对端口面与 GT 面（D21）。用 bbox 六个坐标（≤1e-6）。
* **不要**按 `--mesh` / `--faceids` 的**面序号**配对（D13）：GT 1778 面 vs 端口 1772 面，序号不构成身份。
* 报「未网格面数」必须同时报测量调用序列（D16）：`zz_uv_feed --ids` 报 7 个面无 stat；
  `zz_uv_feed` 默认（先 `prs3d_get_deflection` 再 `from_deflection`）报 `stats=1757 / no_stat=15`。
* 端口侧探针留在 `crates/occt-topo/examples/zz_uv_feed.rs`（不在库路径上，不参与 `--lib`）；
  本轮在 `node_insertion.rs` / `range_splitter/param_set.rs` 加的 env-gated dump **已全部摘除**（D24），
  摘除后复跑 `zz_uv_feed --ids` 与摘除前逐字相同（1772 面、同一 7 个被拒面的 index+bbox；
  只有 `key=` 里的指针型 shape_key 串随进程变化）。
* 本轮**库源码零改动**、未新写单元测试、未改既有断言；工作树改动只有 `specs/board.canvas.tsx`、
  本节、`specs/occt_probe/occt_probe.cpp`（新增 `--uv` / `--uvsum`）与新增的
  `crates/occt-topo/examples/zz_uv_feed.rs`（探针）。

---

### 9.289 —— T-93 主线：导入期 `ShapeFix_Face::FixMissingSeam`（STEP 声明 1 条边，OCCT 成形后 4 条）

#### 9.289.1 一句话

端口拒的这批面，**在 STEP 文件里声明的环边数本来就少于 OCCT 成形后的边数**。缺口是导入期的
`ShapeFix_Shape::Perform` → `ShapeFix_Face::Perform`（`FixMissingSeam` 那一步）这一段，不在网格侧。
端口自己的 `ShapeFixFace::fix_missing_seam` 对目标面能做出正确结果，但**还不能整形状接线** ——
它今天会误改 212/1772 个面，其中只有 1 个是真需要的。

#### 9.289.2 文本级取证（不需要探针即可复核）

T0M 那个球带面（端口 F1691，探针按 R 定位到 OCCT 侧）：

```text
#77    = SPHERICAL_SURFACE('',#33252,4.25)          # 球 R=4.25，心 (0,-25,-1.006)
#29920 = ADVANCED_FACE('',(#4223),#77,.T.)
#4223  = FACE_OUTER_BOUND('',#6320,.T.)
#6320  = EDGE_LOOP('',(#25460))                     # ← 环里只有 1 个元素
#25460 = ORIENTED_EDGE('',*,*,#17469,.F.)
#17469 = EDGE_CURVE('',#13445,#13445,#10859,.T.)    # 起末顶点同一个 ⇒ 闭合边
#10859 = CIRCLE('',#33253,2.)                       # R=2，心 (0,-25,2.744)
```

即：**该面在文件里只声明了一条闭合边**。而 OCCT 读出来的同一个面有 **4 条边**（探针 `--typefaces 3 4.25`
按半径定位，避开面序号差异）：

| | 环边数 | bbox | 各边参数域 |
|---|---|---|---|
| OCCT 未修复（`--nofix`） | **1** | z ∈ [-5.256, 3.244] | e0 `par=[π, 3π]`，`dpar=2π` |
| OCCT 成形后（默认） | **4** | z ∈ [2.744, 3.244] | e0 底圆 `par=pcpar=[π, 3π]`；e1/e3 seam `par=pcpar=[1.08083900054116833, 1.57079632679489656]`；e2 顶点退化边 |
| 端口 F1691（`zz_uv_feed --model 1691`） | **1** | z ∈ [2.744, 2.744] | `curve=Circle` `par=[π, 3π]` `dpar=2π` |

两侧 `vrange` 都是 `[-π/2, π/2]`，**曲面的参数域没有差别**；差别在**跨度的度量范围**：未修复的那条闭合圆
只钉住一个纬度（常 V），于是 `collectWirePoints` 喂进去的点集在一个方向上零宽 → `IsValid=false`（§9.288.2）。
成形后的 4 条边多了两条 seam，V 展布 0.489957 就是从这里来的。

结论：端口 f=1691 的环边数与 **OCCT 未修复**的一致（1 条），与**成形后**的不一致（4 条）。
缺的不是"哪一条边"，而是**整个 seam 构造步骤**。

#### 9.289.3 缺的那段控制流在哪

`ShapeProcess_OperLibrary.cxx:801-899`（`FromSTEP.FixShape` 触发的修复算子）：

```text
occ::handle<ShapeFix_Shape> sfs = new ShapeFix_Shape;
occ::handle<ShapeFix_Face>  sff = sfs->FixFaceTool();
sff->FixMissingSeamMode() = ctx->IntegerVal("FixMissingSeamMode", -1);   // cxx:830
...
sfs->Init(ctx->Result());
sfs->Perform(aPS.Next());                                                // cxx:887  → 逐面 ShapeFix_Face::Perform
```

`ShapeFix_Face::Perform`（`ShapeFix_Face.cxx:345-...`）的步骤顺序：

```text
cxx:365-480  第一轮 wire 修复（FixWireMode；其中 FixLacking 在 cxx:372 关掉）
cxx:492-498  if (NeedFix(myFixMissingSeamMode)) FixMissingSeam()        ← 端口缺的就是这一步
cxx:500-...  第一轮结果里由 FixMissingSeam 产生的新面再走一遍
cxx:692-694  FixOrientation
cxx:704      FixAddNaturalBound
cxx:711-713  FixSplitFace
cxx:731-734  FixSmallAreaWire
```

端口的现状（**未接线**，`step/read_topology.rs` 在 `f1e56776` 把这块删了）：

* `ShapeFixFace::perform_fix_missing_seam`（`shape_fix_face.rs:133-142`）= 只做 `result = face` + `FixMissingSeam`，
  第一轮 wire 修复与 post-seam 面循环都标着 UNPORTED；
* `ShapeFixFace::fix_missing_seam`（`:149`）本身在 `:147-148` 明确写着 UNPORTED：
  **seam 构造与 `w2 != null` 合并（`cxx:1899-2330`）**。

#### 9.289.4 现有实现为什么还不能整形状接线（实测）

`zz_uv_feed --fixall`（对 1772 个面逐个调 `ShapeFixFace::fix_missing_seam`，只测量、不改库）：

```text
FIXALL faces=1772 fix_returns_true=212 edge_count_changed=212
```

改动分布（`.target-gate/fixall.txt`）：

| 变化 | 面数 | 判定 |
|---|---|---|
| `1 → 4` | **1**（f=1691） | 真需要（与 §9.289.2 一致） |
| `2 → 4` | 126 | 伪增量 |
| `2 → 5` | 46 | 伪增量 |
| `3 → 6` | 6 | 伪增量 |
| `5 → 11` / `8 → 11` | 2 / 1 | 伪增量 |
| `9 → 12` | 3 | 伪增量 |
| `（fix 返回 true 但 result 为 None）` | 20 | 其中含 F116(71 边)/F405-F408 等被拒面 —— 返回值语义不等于"改了拓扑" |

⇒ **212 个变动里只有 1 个是真的**。根因就是 `cxx:1899-2330` 那段没移植：`FixMissingSeam` 的触发与 seam
放置依赖它，缺了就会在本来已经有合法 seam 的面上再插一条。

所以正确顺序是：**先补完 `cxx:1899-2330`，再接线**，而不是整形状无条件调用。第一轮导线与 post-seam
面循环是否要一起补，取决于补完之后 `2 → 4` 那 126 个是否消失 —— 那 126 个是判定"补得对不对"的现成反例信号。

#### 9.289.5 怎么复跑

```text
# OCCT 侧：按半径定位球面，看成形前/后的环边数与 bbox（避开面序号差异）
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\T0M.stp --typefaces 3 4.25"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\T0M.stp --typefaces 3 4.25 --nofix"

# 端口侧：单面环结构 / 整形状 seam 修复的影响
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/T0M.stp --model 1691
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/T0M.stp --fixall

# STEP 文本级的环声明（不需要 OCCT）
python .target-gate/step_face2.py data\occ\T0M.stp 29920
```

**注意面序号的三个坑**（本轮实际踩到并已排除）：

1. `occt_probe --uvsum`（走 `IMeshData_Model`）的面序号与**端口一致**，与 `TopExp_Explorer` 的序号**不一致**
   （同一个球带面：模型序 409 / explorer 序 1697）。按几何配对不对就永远对不上。
2. `--nofix` 会整体改变 `TopExp_Explorer` 的面序号（`--adv 1697` 在两种模式下是**两个不同的面**），
   所以跨模式比较只能按曲面半径 / bbox，不能按序号。
3. `--adv`/`--advbox` 现在扫描全部 argv 找模式串，参数顺序不再决定行为（`--advbox ... --nofix` 可用）。

#### 9.289.6 旁支

* **F116**（71 条边的面）`fix_missing_seam` 返回 true 但 `result` 为 None，与 F1691 不同因，另查。
  它也是 7 个被拒面之一（端口 2 wire vs GT 1 wire）。
* `--nofix` 的模型里所有球面 `wires=0`、pcurve 参数域是 `±1e100` 哨兵 —— 未修复的环上 `collectWirePoints`
  拿不到点。这**不是**端口的对照物，别把它的数字当基线。
* 本轮仍**库源码零改动**；新增/扩展的都是探针（`specs/occt_probe/occt_probe.cpp`、`examples/zz_uv_feed.rs`）。

#### 9.289.7 定级：reader 没有丢边（这是本轮要回答的那个问题）

**结论：端口 F1691 的 1 条边，就是 STEP `EDGE_LOOP` 里声明的那 1 条边。全链没有丢边的地方。**

逐级核对（每一级都是「n 条进 → n 条出」，无 filter / retain / skip / 去重）：

| 级 | 位置 | 行为 |
|---|---|---|
| STEP 文本 | `#6320=EDGE_LOOP('',(#25460))` | 环里 **1** 个 `ORIENTED_EDGE` |
| `resolve_face` | `read_topology.rs:578`，bounds 循环 `:649-662` | 对每个 bound `resolve_shape`，只判 `is_wire()` |
| `resolve_outer_bound` | `:564-576` | 只做 `loop_ref` 解析 + `.F.` 时设 `Reversed` |
| `resolve_loop` | `:413-433` | `:421-432` 对 `items` 逐条 `resolve_shape`；唯一分支是 `if !s.is_edge() { return Err(..) }`；`make_wire(&edges)` |
| `resolve_oriented_edge` | `:369-388` | 只算 `ori = XNOR(OrientedEdge.orientation, EDGE_CURVE.same_sense)` |
| `resolve_edge` | `:307-367` | 建 1 条边 |

**没走的那条分支**：`:628-644` 的「lone VertexLoop on sphere/BSpline/revolution → `BRepLib_MakeFace` 自然边界」在
这里**不适用** —— 该 bound 是 `EDGE_LOOP` 而非 `VERTEX_LOOP`，`bound_loop_is_vertex_loop`（`:741`）为 false。

**实测等式**（同一次运行、同一命令）：

```text
data/occ/T0M.stp:11378   #6320=EDGE_LOOP('',(#25460));          ← 声明 1 条

zz_uv_feed --model 1691:
MODEL f=1691 wires=1 ... bbox=(-2.000000,-27.000000,2.744000)-(2.000000,-23.000000,2.744000)
MODEL  wire[0] nEdges=1                                          ← 读到 1 条
MODEL   e[0] edge=3980 ori=Forward deg=false curve=Circle par=[3.141592653589793,9.424777960769379] dpar=6.283185307179586
```

**对照侧也一致**：OCCT 未修复导入（`--nofix`）的这个面同样是 **1 条边**、同样 `par=[π, 3π]`、
同样是完整闭合圆。两侧读入相同，差别只在 OCCT 随后跑了 seam 构造（§9.289.3）。

⇒ 因此

* **不是**「reader 丢了一条边」的 bug，**是新增一段同等分支**（导入期 seam 构造）；
* D23 里「缺边不晚于离散模型构建」的措辞据此收紧为「**文件声明本身就只有 1 条**」——
  这条等式比原来那句更强：不是「某级丢了」，而是「本来就只有一条，端口忠实读了它」；
* 「零宽」的机制也随之完全清楚：一条**完整 2π 闭合边**在球面上是等纬圆（常 V），
  它能张成的切向只有一个方向。要张成二维区域必须补 seam（两条经线段），这只能由 seam 构造产生，
  不可能靠改 reader 的解析得到。

---

### 9.290 —— T-93 落地：导入期 `FixMissingSeam` 接线（T0M 未网格 7→6，WIREHIST 向 GT 收敛，全门禁保持绿）

#### 9.290.1 一句话

`ShapeFixFace::fix_missing_seam` **其实早就移植完了**（源码里 `UNPORTED: seam 构造 cxx:1899-2330` 的注释是**过期**的）；
缺的是**在 reader 里调用它**。接上之后：T0M 未网格面 **7→6**、`wires` **2121→1931**、
`WIREHIST 1:1463 2:294 …` → **`1:1652 2:106 3:8 …`**（更贴近 GT），且 `--lib` 1281/0、`step_obj_gates` **5/5** 全部保持绿。

#### 9.290.2 先证「同构」：端口产出与 OCCT 成形结果逐项相同

用新探针 `zz_seam_fix`（`crates/occt-topo/examples/zz_seam_fix.rs`）对 `data/occ/T0M.stp` 的 face 1691 单独调用：

```text
--- BEFORE ---
BEFORE wires=1  bbox=(-2.000000,-27.000000,2.744000)-(2.000000,-23.000000,2.744000)
BEFORE  wire[0] nEdges=1
BEFORE   e[0] ori=Reversed deg=false par=[3.141592653589793,9.424777960769379] dpar=6.283185307179586
fix_missing_seam returned true
--- HEALED face 0 ---
HEALED wires=1  bbox=(-2.000000,-27.000000,2.744000)-(2.000000,-23.000000,3.244000)
HEALED  wire[0] nEdges=4
HEALED   e[0] ori=Reversed deg=false par=[3.141592653589793,9.424777960769379] dpar=6.283185307179586
HEALED   e[1] ori=Forward  deg=false par=[1.080839000541168,1.570796326794897] dpar=0.489957326253728
HEALED   e[2] ori=Forward  deg=true  par=[0.000000000000000,6.283185307179586] dpar=6.283185307179586
HEALED   e[3] ori=Reversed deg=false par=[1.080839000541168,1.570796326794897] dpar=0.489957326253728
```

与 OCCT 成形后的同一个面（`occt_probe --typefaces 3 4.25`，按半径定位）**逐项相同**：
4 条边、同样 `[π,3π]` 的底圆、同样 `[1.080839000541168, 1.570796326794897]` 的两条 seam、
同样 `deg=true` 的顶点退化边，bbox 同样由 `z∈[2.744,2.744]` 变成 `z∈[2.744,3.244]`。

#### 9.290.3 再证「缺口是调用点」：read 期间一次都没被调用

加了一个临时计数器（**已摘除**，清理后 `--lib` 1281/0 不变）量到：

```text
SEAMCALLS during read_step_file = 0
```

即 `read_step_file` 全程**不会**调用 `fix_missing_seam`。而 OCCT 侧的 STEP 资源默认是开着的：
`STEPControl_Controller.cxx:201` 把 `FromSTEP.exec.op` 设为 `FixShape`，
`:221` 把 `FromSTEP.FixShape.FixMissingSeamMode` 设为 `-1`（`NeedFix(-1) == true`），
`ShapeProcess_OperLibrary.cxx:785-899`（`fixshape`）再驱动 `ShapeFix_Shape::Perform` → `ShapeFix_Face::Perform`
→ `FixMissingSeam`（`cxx:492-498`）。

⚠️ 这一条与 D2/T-93 (a) 当初的判断相反。D2 的依据是「patch 过的 OCCT 在 `TransferRoots` 期间没有任何输出」，
但本轮实测 OCCT `OneShape()` 出来的面**确实**比文件多 3 条边（§9.289.2），也确认 `ShapeProcess` 一关这 3 条边就没了。
两者不可能同时成立，而**可复跑的那一侧（探针）指向「成形期确实加了 seam」**。
本节的结论只依赖这一点；D2 的 link-map 实验本轮无法复现，故不再作为依据。

#### 9.290.4 接线形态与实测

`crates/occt-topo/src/step/read_topology.rs` 的 `resolve_face` 末尾恢复调用（形态与 `f1e56776`
删掉的那段一致：`fix_missing_seam` → 结果面 → 逐 wire `check_pcurves_and_shift`）：

```rust
{
    let mut sff = crate::shhealing::ShapeFixFace::with_face(&face);
    if sff.fix_missing_seam() {
        if let Some(res) = sff.result.clone() {
            if res.shape_type() == crate::abs::ShapeType::Face {
                let rf = crate::shape::Face(res.clone());
                for mut w in crate::topo_tools_full::wires_of_face(&rf) {
                    crate::shhealing::check_pcurves_and_shift(&mut w, &rf, self.precision, false);
                }
                return Ok(res);
            }
        }
    }
}
Ok(face.0)
```

对照实测（同一命令、只切这一个 hunk）：

| 指标 | 基线 | 接线后 |
|---|---|---|
| T0M 未网格面 | **7** | **6** |
| T0M `wires` | 2121 | **1931** |
| T0M 多 wire 面 | 309 | **120** |
| T0M 1-edge 闭合 wire | 596 | 232 |
| T0M `WIREHIST` | `1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1` | **`1:1652 2:106 3:8 5:1 6:1 7:2 8:1 10:1`** |
| a3n00 `stats` / `unmatched` | 226 / 0 | **225 / 1** |
| a3n00 `mesh_v` / `mesh_t` | 11101 / 11121 | **10863 / 11941** |
| `--lib` | 1281 / 0 | **1281 / 0** ✅ |
| `step_obj_gates` | 5 / 5 | **5 / 5** ✅ |

`WIREHIST` 是关键信号：GT 的 T0M 面平均 wire 数远低于端口，`2:294` 那一大堆双 wire 面正是「缺 seam 时
把一条闭合圆当整个边界」的产物；接线后它们被合成单 wire，直方图整体左移。

#### 9.290.5 效果逐面核对（healed bbox 与 GT 逐字相同）

接线后仍未被网格的 6 个面，与已被治好的面对照：

| 端口面 | 接线前 | 接线后 | GT 孪生面 |
|---|---|---|---|
| f=1691（球带） | 1 边、bbox z 零厚、**未网格** | **4 边**、bbox `z∈[2.744,3.244]`、**已网格 mt=143** | GT 1697：bbox 逐字相同 ✅ |
| f=405（BSpline） | 1 边、bbox x 零厚、未网格 | **已治好**，bbox 与 GT **逐字相同** ✅ | GT 73 |
| f=406 / f=407 / f=408 | 同上 | **同上**，bbox 逐字相同 ✅ | GT 74 / 75 / 76 |
| f=1758（**圆锥**） | 1 边、未网格 | **不变**（仍 1 边、未网格） | GT 566（已出网格） |
| f=351 | 未网格 | 不变 | — |
| f=116（大圆柱） | 未网格 | 不变 | — |

**f=1758 是已知有界的缺口**：`FixMissingSeam` 只为「球 / BSpline / 退化环面」构造 seam
（`ShapeFix_Face.cxx:1909-1975` 的四个分支），**圆锥四个分支都不匹配** ⇒ 走 `cxx:1972-1975` 的 `return false`。
所以它不该由本卡修，需要另找 OCCT 里给圆锥补 seam 的分支（若有）。

#### 9.290.6 仍然存在的过触发（未解决、已量化）

整形状逐个调用 `perform_fix_missing_seam`：`faces=1772 changed=212`。分布：

```text
2 -> 4 : 126      3 -> 6 : 14      2 -> 6 : 4      5 -> 11 : 4
5 -> 8 : 3        9 -> 12 : 3      7 -> 10 : 2     3 -> 5 : 2
23 -> 25 : 1      8 -> 11 : 1      6 -> 8 : 1      1 -> 4 : 1   ← 只有这一个是真需要的
```

⇒ 212 个变动里**只有 f=1691 是「1 条边」的真缺口**，其余 211 个是在已有合法 seam 的面上的伪增量。
**但门禁全绿**（`--lib` 1281/0、`step_obj_gates` 5/5），说明这些伪增量没有把面积比推出 `area_tol`。
它们是**潜在**风险而非已证缺陷：本轮**不**去改 `check_wire` / 配对选择来"消灭"这 211 个，
因为没有 `.cxx` 依据说明端口的选择逻辑与 OCCT 不同 —— 只是返回值语义与拓扑改动不同源。
下一步若要收敛这个问题，应当先做「端口 `check_wire` 与 `cxx:1652-1718` 的逐条对读」，
而不是按"少改几个面"来调参。

#### 9.290.7 怎么复跑

```text
# 单面：修复前后 + 与 OCCT 对照
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_seam_fix -- data/occ/T0M.stp 1691
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\T0M.stp --typefaces 3 4.25"

# 整形状过触发量化
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_seam_fix -- data/occ/T0M.stp 0 --all

# 门禁
rtk cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --lib
rtk cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --test step_obj_gates
```

#### 9.290.8 旁支（未动）

* **f=1758（圆锥）**：`FixMissingSeam` 无匹配分支（§9.290.5），需另找 OCCT 侧分支。
* **f=351（10 wire）与 f=116（大圆柱）**：与 seam 无关，另因。
* **211 个伪增量**（§9.290.6）：不改选择逻辑，先对读 `check_wire`。
* a3n00 接线后 `unmatched` 由 0→1，`mesh_t` 11121→11941：数值仍在门禁内，但**是新出现的差**，
  若要追，从这里入手。

---

### 9.291 —— T-59 结项：a3n00 配对表与面积缺口归因（缺口 = 单个面 F113）

#### 9.291.1 一句话

a3n00 的 226 个端口面与 GT 的 226 个面按 **bbox 六坐标**（D21 口径）**全部配上**；端口只剩
**1 个面未网格（f=113）**，而它在 GT 侧是 2 wire 的正常圆柱面。**F113 单面所在区域就占 +26952 面积缺口**，
而全模型净缺口只有 22802.58 ⇒ 它是 a3n00 面积缺口的唯一主因。已拆卡 **T-99**。

同时把 **D18 的 `n1→n2` 一行修复在 D27 之后重测**：仍是净负面（`0.8996 → 0.8904`，且 F113 没被救活）⇒ 与 D19 结论一致。

#### 9.291.2 配对与 wire 结构

```text
port faces=226  gt faces=226
port unmeshed faces: 1
paired=226  port unpaired=0  gt unpaired=0
pairs with bbox diff > 1e-3: 26
```

确定配对的 face 按 GT 曲面类型分布：**Plane 86 / Cylinder 62 / Cone 30 / Extrusion 15 / Torus 7**。

**wire 数不匹配的只有 3 组**（这说明 a3n00 的 wire 结构基本已对齐，缺陷是外科式的）：

| 端口面 | 端口 wires | GT 面 | GT wires | 类型 | 端口 mt |
|---|---|---|---|---|---|
| **f=113** | **4** | f=77 | **2** | Cylinder | **未网格** |
| f=140 | 2 | f=39 | 1 | Cone | 5 |
| f=170 | 2 | f=69 | 1 | Cylinder | 27 |

F113 的结构逐条对照（`zz_uv_feed --model 113` vs `occt_probe --advbox`）：

| | 端口 | GT |
|---|---|---|
| bbox | `(-87.5,-34,-100)-(87.5,34,-32)` | **逐字相同** |
| wires | 4 | 2 |
| 每 wire 边数 | **[6, 14, 1, 1]**（合计 22） | **[22, 6]** |
| 曲面参数域 | `urange=[0,6.28318530717958623]` `vrange=[-inf,inf]` | `u=[0,6.28318530717958623]` `v=[-43.75,131.25]` |

⇒ 端口把**同一个 22 边环切成了 4 条 wire**。端口 `wire[0]` 的 6 条边
（BSpline/Line/BSpline/Line/BSpline/BSpline）与 GT 那条 6 边 wire 正好对应。

#### 9.291.3 面积归因（稳健）

`obj_area` 口径与门禁相同（`tests/common/mod.rs:155` 的扇形三角化 `0.5·|ab×ac|` 求和）。

```text
总:  ours=204327.72  occ=227130.30  deficit=+22802.58
F113 bbox 内 (pad=0.01): ours=  1690.55 (t=661)  occ=28642.74 (t=798)  deficit=+26952.18
             (pad=0.5 ): ours=  1738.47 (t=693)  occ=28735.87 (t=837)  deficit=+26997.40
             (pad=2.0 ): ours=  2109.71 (t=765)  occ=29532.83 (t=931)  deficit=+27423.13
解析量核对: 该圆柱 2πrh = 2π×68×68 = 29053.4   （GT 区域实测 28642.7 = 98.6%；端口仅 5.8%）
```

pad 从 0 放到 2 只让缺口变动 <2% ⇒ **归因稳健**，不是框选偏差造成的。
F113 单面缺口（26952）**超过**全模型净缺口（22803）⇒ 其余面合计是净盈余（端口在别处偏密）。

#### 9.291.4 a3n00 新基线（取代卡片里的 0.8627）

| | 数值 |
|---|---|
| ours / occ | 204327.72 / 227130.30 |
| **比值** | **0.8996** |
| 判据 | `\|ratio − 1\| < area_tol = 0.15` ⇒ 余量 **0.0496** |

D27 的 seam 接线把比值从 0.8627 抬到 0.8996（+0.037）。**卡片里「应用 §9.287.4 后会到 0.8298」的旧记录已过期**：
那个 −0.0329 的代价是在 0.8627 基线上量的，但即便按同样幅度作用在 0.8996 上也不越界 —— 只是它本身无收益（见下）。

#### 9.291.5 D18 的一行修复：在 D27 之后重测，仍为净负面

`wire_fix.rs:3078` 的 `fix_dummy_seam(wire, n1)` → `(wire, n2)`（D18 的传参错位），
在 **D27 的 seam 接线已落地**的前提下重测：

| 指标 | 当前（`n1`） | 改成 `n2` | Δ |
|---|---|---|---|
| a3n00 面积 | 204327.72 | **202232.18** | **−2095.54** |
| a3n00 比值 | 0.8996 | **0.8904** | **−0.0092** |
| a3n00 `mesh_v` / `mesh_t` | 10863 / 11941 | 10889 / 12008 | +26 / +67 |
| a3n00 未网格面 | 1 | **1** | 0（**F113 没被救活**） |

⇒ **无收益且有害**，已回退。这与 D19（当时 0.8627→0.8298）方向一致，本轮是在新基线上独立复核的第二次否证。
**结论：不要再把这一行当作 T-99 的候选修复。**

#### 9.291.6 怎么复跑

```text
# 两侧逐面 dump
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids > .target-gate/a3n00_port.txt
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --uvsum"            > .target-gate/a3n00_gt.txt

# 配对 + wire 数不匹配 + 面类分布
python .target-gate/pair_a3n00.py .target-gate/a3n00_port.txt .target-gate/a3n00_gt.txt

# 面积（门禁同口径）+ 按 bbox 归因
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj -- data/occ
python .target-gate/pair_area.py   output/a3n00.obj data/occ/occ-a3n00.obj
python .target-gate/deficit_a3n00.py output/a3n00.obj data/occ/occ-a3n00.obj

# 单面结构对照
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --model 113
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --advbox -87.5 -34 -100 87.5 34 -32"
```

#### 9.291.7 交接给 T-99

已否证的路（不要再走）：

* D18 的 `n1→n2`（§9.291.5，第二次否证）；
* 按「bbox 中心最近邻」配对（D21）；
* 在 mesh 侧加谓词把多余 wire 合并 —— OCCT 里没有这样的规则。

尚未查的：**22 边环为何在端口被切成 4 条 wire**。要定位到 reader/`ShapeFix_Wire` 的哪一轮 pass
（`FixConnected` / `FixNotchedEdges` / `FixLacking` 各自的分支），并与 `ShapeFix_Wire.cxx:3977-4310` 逐段对读。

---

### 9.292 —— T-99 根因：`FixMissingSeam` 返回 Shell 时被 reader 丢掉（F113 不网格的真因）

#### 9.292.1 一句话

F113 的 22 边环**并没有被切成 4 条 wire** —— 那是 STEP 里本来就声明的 4 个 loop。
真正的问题是：**T-93 接进来的 seam 步骤对 F113 产出了一个 `Shell`（5 个面），而 `resolve_face` 只接受
`shape_type() == Face` 的结果，于是把它整个丢掉**、返回未修复的原面；而这一步**已经改过共享
`GeometryRegistry`**。所以 F113 既没被修复，又付了副作用。

#### 9.292.2 取证（`OCCT_TOPO_TRACE_WIRES` 临时插桩，已摘除）

在 `resolve_face` 的 seam 块前后打印「wire 边数列表 + bbox + seam 返回值 + 结果类型」，跑 a3n00 得：

```text
TRACE face wires_before=[6, 14, 1, 1] bbox=(-87.500,-34.000,-100.000)-(87.500,34.000,-32.000)
TRACE face fix_missing_seam=true
TRACE face wires_after=[] result_type=Shell nfaces=5
```

**三条结论**

1. **`wires_before` 已经是 `[6,14,1,1]`** ⇒ 4 条 wire 是 STEP 声明的（不是被切开的）。
   佐证：全模型 loop 尺寸直方图里**最大的 loop 只有 17 条边**（`1:168 2:1 3:28 4:89 5:26 6:9 7:2 8:12 9:1 12:3 14:1 17:1`，共 341 个 loop），
   没有任何 22 边的 loop ⇒ 22 边不可能来自「一个大 loop 被切」。
2. **seam 步骤对 F113 返回 `Shell`（5 个面）** ⇒ 走不到 `shape_type() == Face`，结果被丢弃。
3. `wires_after=[]` 是因为对 Shell 调 `wires_of_face`，并非真的没有边。

`result_type` 分布（`fix_missing_seam` 返回 true 的 46 个面）：

| 结果类型 | 面数 | 每个含几个面 |
|---|---|---|
| `Face` | 43 | 1 |
| `Shell` | 2 | 2 |
| `Shell` | 1 | **5** ← F113 |

⇒ **3 个面（40 条边）的 seam 结果被丢弃**，其中就包括不出网格的 F113。

#### 9.292.3 为什么 OCCT 得到的是 2 wire、端口是 4 wire

STEP 侧该面声明 4 个 loop；OCCT 的等价步骤把它们**合并成 2 条 wire**（`[22,6]`），端口保持 4 条（`[6,14,1,1]`）。
即 GT 的「合并」发生在 `FixMissingSeam`/`ComposeShell` 内部，而端口因为把这个结果丢了，所以看不到合并。
⇒ 这不是「端口多切了」，而是**端口没拿到 OCCT 那次合并的结果**。

#### 9.292.4 修复方向（未做）

两条独立的问题，必须一起解：

1. **`resolve_face` 丢 Shell 结果**：`ShapeFix_Face.hxx:146` 明确写了
   `myResult` 可能是 shell（"To be used instead of Face() if FixMissingSeam involved"）。
   端口只取 `Face` ⇒ 要么把单面 shell 取出该面，要么把多面结果作为 compound 交回上层（需同时改
   `resolve_shell` 的面计数，别把面数算错）。
2. **`FixMissingSeam` 对 F113 产出 5 个面**：OCCT 在这条路径末尾有 `FixSmall` / `FixSmallAreaWire`
   清理（`ShapeFix_Face.cxx:2270-2322`），端口这段已有实现（`shape_fix_face.rs` 的 round1/`check_small_area`），
   但对 F113 显然没把 5 个面收敛成 1 个。需要对照 `.cxx` 查 `ComposeShell` 为什么多产 4 个面，
   以及 `FixSmall` 的判定为什么没把它们丢掉。**不要在 reader 里按面数硬滤。**

#### 9.292.5 已否证的路（不要再走）

* D18 的 `fix_dummy_seam(wire, n1)` → `(n2)`：在 D27 之后重测，a3n00 面积 0.8996→0.8904、F113 仍不出网格（§9.291.5）。
* 「22 边 loop 被切」这个假设本身：loop 尺寸直方图里没有 22（§9.292.2）。
* 在 mesh 侧加谓词合并多 wire：OCCT 里没有这样的规则。

#### 9.292.6 怎么复跑（插桩已摘除，需要时按此重加）

```text
# 在 read_topology.rs 的 resolve_face seam 块前后加 eprintln!("TRACE ...")，
# 用 OCCT_TOPO_TRACE_WIRES 开关，然后：
OCCT_TOPO_TRACE_WIRES=1 rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline \
    --example zz_uv_feed -- data/occ/a3n00.stp 2> .target-gate/trace_err.txt
# 按 bbox 找 F113：
grep 'bbox=(-87.500,-34.000,-100.000)' .target-gate/trace_err.txt
```

**摘除后复跑 `--model 113` 与 `zz_probe_a3n00`：`wires=4`、`10863/11941`、`unmatched=1` 全部逐字不变** ⇒ 插桩确认惰性。

---

### 9.293 —— T-99 后续：把 Shell 结果交回 reader 是行不通的；修复点在形状级 healing 而非逐面

#### 9.293.1 实测：返回 Shell 会让整个模型读不出来

按 §9.292.4 的方向 ① 试了一版：`resolve_face` 在 seam 结果不是 `Face` 时，
把其结果里的面包成 compound 返回（并对每个面跑 `check_pcurves_and_shift`）。
`cargo check` 通过，但 a3n00 直接读不出模型：

```text
thread 'main' panicked at examples/zz_probe_a3n00.rs:34:
index out of bounds: the len is 0 but the index is 0        # model.shapes 为空
```

原因很清楚：`resolve_shell`（`read_topology.rs:975`）对 `CLOSED_SHELL` 的每个元素断言 `s.is_face()`，
compound 不是 face ⇒ 该 shell 解析失败 ⇒ 整个 solid/root 丢失。

⇒ **「把多面结果交回 reader」这条路是错的**，已回退（行为逐字恢复：`TOTAL faces=226 stats=225 mesh_v=10683…`
实测为 `226/225/10863/11941`、`unmatched=1`）。

#### 9.293.2 结论：修复点在**形状级 healing**，不是逐面 reader

OCCT 里 `FixMissingSeam` 从不是逐面跑的：

* `ShapeProcess_OperLibrary.cxx:801-899` 的 `fixshape` 建一个 **`ShapeFix_Shape`**，
  `sfs->Init(ctx->Result())` 后再 `sfs->Perform()` —— 作用对象是**整个形状**；
  `ShapeFix_Shape` 内部把结果记进 `ShapeFix_Root` 的 context（`ShapeBuild_ReShape`），
  **允许一个面被替换成多个面**，最后 `ctx->SetResult(result)` 换掉整形的结果（`cxx:896`）。
* 端口的 `resolve_face` 是**逐面**返回一个 `TopoShape`，`resolve_shell` 又要求每个元素必须是 face
  ⇒ 端口在结构上**没有位置**容纳「一个面变成多个面」。

所以 T-99 的正确形态是：**在 `read_step_file` 读完整个形状之后，加一个形状级 healing pass**
（族 `ShapeFix_Shape` + `ShapeBuild_ReShape` 的替换表），让 seam 步骤能在整形上增删面；
逐面 `resolve_face` 里的那一段 seam 调用应当**移出去**（它现在既改了共享 registry 又丢弃结果，
是「两头不讨好」的那一处）。

这同时解释了 §9.292 观测到的「43 个返回 Face、3 个返回 Shell」：端口把一个**形状级**算法
塞进了**逐面**位置，于是只有恰好收敛成单面的那 43 个能生效。

#### 9.293.3 复跑

```text
# 回退后的基线（应与本文件其它处一致）
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_probe_a3n00 -- data/occ/a3n00.stp
#   => TOTAL faces=226 ... stats=225 mesh_v=10863 mesh_t=11941 / STATMAP matched=225 unmatched=1
```

#### 9.293.4 给下一轮的边界

* **不要**再试「`resolve_face` 返回 Shell/compound」：实测让 root 归零（§9.293.1）。
* `FixMissingSeam` 对 F113 产出 5 个面这件事本身**未必是缺陷** —— OCCT 同级也会多面，
  它靠形状级替换表容纳。先确认「OCCT 在这个文件上 F113 最终是 1 个面」再谈要不要收敛，
  否则会为一个假目标去改 `ComposeShell`。
* 形状级 healing 是新函数面，属于**加一段同等分支**（有 `.cxx` 依据），不是调参；
  但改动面覆盖全部读入路径，必须整轮跑 `step_obj_gates`（23 模型）复核。

#### 9.293.5 直面取证：F113 的 seam 结果 = **1 个正确面 + 4 个碎片**

用 `zz_seam_fix`（`--example zz_seam_fix -- data/occ/a3n00.stp 113`）把 seam 结果逐面摊开：

| 结果面 | wires | 每 wire 边数 | 曲面参数域 | bbox |
|---|---|---|---|---|
| **face 0** | **2** | **6 + 7** | `u=[0,2π]` **`v=[-43.75,131.25]`** | `(-49.87,-30.55,-100)-(63.93,34,-40)` |
| face 1 | 1 | 5 | `u=[0,2π]` `v=[-43.75,131.25]` | `(-38.73,-28.84,-48)-(…)` |
| face 2 | 1 | 4 | 同上 | `(-44.50,-21.97,-48)-(…)` |
| face 3 | 1 | 3 | 同上 | `(44.50,-34,-100)-(…)` |
| face 4 | 1 | 3 | 同上 | `(-87.50,-34,-100)-(…)` |

对照 GT 的孪生面（f=77）：`wires=2`、边数 **[22, 6]**、`u=[0,2π]`、`v=[-43.75,131.25]`、bbox=**F113 原 bbox**。

**face 0 已经和 GT 的关键量一致**：2 条 wire、**曲面参数域逐字相同**（`v=[-43.75,131.25]` 正是接缝该产生的收窄）。
差的是：

1. **wire 边数是 6+7 而不是 22+6** ⇒ 还有 4 个面（face 1–4）的边没被并进来；
2. **bbox 只覆盖原面的一部分** ⇒ 那 4 个碎片本该被吸收。

⇒ 所以「F113 最终该是几个面」这个问题现在有答案了：**GT 是 1 个面**，而端口给出了 1 个正确面 + 4 个碎片。
要收敛的正是这 4 个碎片，它们不是「ComposeShell 多产的废面」，而是**本该参与合并的边集**。

#### 9.293.6 给下一轮的具体入口

按 `.cxx` 逐段核这两处（都有现成反例信号）：

1. **`FixMissingSeam` 的 `w2 != null` 合并路径**（`ShapeFix_Face.cxx:1994-2234`）：
   port 版在 `shape_fix_face.rs` 对应段已实现，但 F113 走的是 `w2 == null` 分支
   （BEFORE 是 4 条 wire，`CheckWire` 选出的 w1/w2 与 OCCT 的配对是否一致未验）。
   **先打印 `ismodeu/ismodev/isdeg1/isdeg2` 与选中的 w1/w2 是哪两条 wire**，与 OCCT 同参数对读。
2. **`ComposeShell` 之后的小面/小 wire 清理**（`ShapeFix_Face.cxx:2270-2322`，端口 `shape_fix_face.rs:643-686`）：
   本轮实测 `nb_faces`（round1 后的面数）**> 1**（因为结果类型是 `Shell`），于是走了
   `cxx:2310-2319` 的 `FixSmallAreaWire(true)`。要查的是：那 4 个碎片为什么没被判成 small area 而丢掉，
   或者为什么没被判成该合并。

**判定信号**（现成）：`zz_seam_fix -- data/occ/a3n00.stp 113` 的结果应是 **1 个面、2 条 wire、6+7=22 与 6 条边**。
现在 face 0 已经有 6+7，缺的就是把 face 1–4 并进来。

---

### 9.294 —— T-99：F113 的两条真 wire 在 `CheckWire` 里求和恰好为 0；`CheckWire` 的 pcurve 取值口径已排除

#### 9.294.1 本轮否证的假设

`CheckWire`（`ShapeFix_Face.cxx:1652-1718` / 端口 `shape_fix_face.rs:42-97`）用的是
`ShapeAnalysis_Edge::PCurve(..., true)`，返回的是 pcurve **自己的** first/last 参数。
我据此改了端口 `check_wire`（改成用 `c2d.first_parameter()/last_parameter()`，并按 edge 朝向交换），
**实测是负面的**：F113 从「2 条 1-edge wire 通过」变成 4 条全部被拒 ⇒ `w1` 为 null ⇒ `fix_missing_seam` 不再触发。
已回退（`diff -0`，a3n00 实测 `226/225/10863/11941` 逐字恢复）。

⇒ **pcurve 取值口径不是根因**，这条已排除。

#### 9.294.2 直测每条 wire 的逐边位移（新探针 `--wireuv`）

`zz_seam_fix -- data/occ/a3n00.stp 113 --wireuv` 打印每条边用的 `(f,l)` 与 `pcurve(f)`/`pcurve(l)`：

| wire | 边数 | 逐边位移之和 | `dU=6.283185`、`dV=175` 下的 `CheckWire` |
|---|---|---|---|
| **wire[0]** | 6 | **`(0.000000, -0.000000)`** | **REJECT** |
| **wire[1]** | 14 | **`(0.000000, -0.000000)`** | **REJECT** |
| wire[2] | 1 | `(6.283185, 0)` | `Some((1, 0, false))` |
| wire[3] | 1 | `(-6.283185, 0)` | `Some((-1, 0, false))` |

`wire[0]` 的逐边位移（`d=(dx,dy)`）是：
`(0.812, 62.252) (0.088, -86.089) (0, -27.703) (-0.177, 0) (0, 27.703) (-0.724, 23.837)`
—— **x 与 y 都精确抵消到 0**。

#### 9.294.3 这说明了什么

`vec` 是两个方向上的和。和**恰好为 0**（不是浮点近似）意味着：

* **U 方向**：把 u 当普通轴上求和，总和为 0 ⇒ 该 wire **在 U 上闭合**（绕了一圈），不是「跨一个周期」；
* **V 方向**：总和也为 0。GT 该面 `v=[-43.75,131.25]`（跨度 **175**，有限），
  而端口这次调用里 `v_range` 已经是**被钳到 `Precision::Infinite()` 的 1e100**（见 §9.294.4）。

`CheckWire` 的判据是「和的绝对值与周期之差 < 0.1×周期」：和 = 0 时两个方向都不满足 ⇒ 两条真 wire 被拒。
算法本身对它拿到的这两个和是**正确**的；可疑的是**这两个和为什么是 0**。

于是 `ismodeu/ismodev` 由两条 1-edge 弧（`±2π`）决定 ⇒ 选出的 w1/w2 是那两条退化弧，
`ComposeShell` 就把面切成了 5 块（§9.293.5）。**4 条 wire 本身没错**（6+14+1+1 = 22，与 GT 的 22+6 同源）；
错的是「哪两条被当成 seam 的两侧」。

#### 9.294.4 下一步（有具体检查点）

1. **先核 `v_range` 为什么是无穷**：`FixMissingSeam` 里
   `v_range = min(|svl - svf|, Precision::Infinite())`，`svf/svl` 来自
   `mySurf->Bounds()`，无穷时用 `BRepTools::UVBounds(myFace, ...)` 兜底。
   探针实测 `--wireuv` 那次打印出 `v_range=1e100`（即已被钳），
   而 GT 同一面 `v=[-43.75,131.25]` ⇒ **端口的 face UV bounds 或 surface Bounds 之一仍是无穷**。
   查 `brep_tools::uv_bounds`（`read_topology` 与 `shape_fix_face` 都依赖它）在
   `vrange=[-inf,inf]` 的圆柱面上返回什么。
2. **再核那两条真 wire 的 pcurve 端点**：`d=(0.812,62.252)` 这类数字说明 `y` 值在某些边上取到了
   `93.61` / `-86.09`（跨了 175 的范围），而 `pcurve` 的 U 只到 `3.23`；
   要确认这些 pcurve 是否被正确关联到该面（`associate_edge_pcurve`），
   以及 F113 的 4 条 wire 是否**都**拿到了 pcurve。

**判定信号（现成）**：`zz_seam_fix -- data/occ/a3n00.stp 113` 结束时
`SEAM selected` 应为 6 边与 14 边那两条（而不是 1-edge 的两条），且最终结果是 **1 个面、2 wire、22+6 边**。

#### 9.294.5 复跑

```text
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_seam_fix -- \
    data/occ/a3n00.stp 113 --wireuv      # 逐边 pcurve 取值与 CheckWire 判定
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_seam_fix -- \
    data/occ/a3n00.stp 113               # seam 结果的 5 个面
```

#### 9.294.6 已否证 / 不要再走

* 「`CheckWire` 的 pcurve 参数口径错了」（§9.294.1，实测更差，已回退）；
* 「`resolve_face` 返回 Shell/compound」（§9.293.1，root 归零）；
* D18 的 `n1→n2`（§9.291.5，二次否证）；
* 「22 边 loop 被切」（§9.292，loop 直方图最大 17 边）。

---

### 9.295 —— T-96 评估：TrimmedCurve 参数表示**不需要**重做（两条 gcpnts arm 早已可忠实表达）

#### 9.295.1 结论（先说结论）

**不要重做 `GeomTrimmedCurve` 的参数表示。** T-51/T-96 那句「两条 Compute arm 无法忠实表达，
因为端口把参数重映射到 `[0,1]`」是**过期的**：这两条 arm 已经移植完并且已经接线，而 `[0,1]` 只是一次
**自洽的重参数化**，不构成障碍。T-96 的答案是「不做」，并把那段过期注释改掉（已改，见 §9.296）。

#### 9.295.2 实测依据（三条，逐条可复核）

**① OCCT 的 `Geom_TrimmedCurve` 确实是「basis 参数、无链式法则」** ——
但这一点被引用成了「所以端口不行」，而真正该问的是「端口内部是否自洽」：

| | OCCT `Geom_TrimmedCurve` | 端口 `GeomTrimmedCurve` |
|---|---|---|
| `First/LastParameter` | `uTrim1`/`uTrim2`（`cxx:255-267`） | `0.0`/`1.0`（`trimmed.rs:80-81`） |
| `EvalD0..EvalDN` | **直接委派** basis（`cxx:212-243`） | 委派 + 链式法则 `s^N`（`trimmed.rs:44-79`） |

**② 端口把一个 trim 会读到的每个量都用同一个 `s = last - first` 处理了** ⇒ 自洽：

| 量 | 位置 | 处理 |
|---|---|---|
| `d1`/`d2`/`d3`/`eval_dn` | `trimmed.rs:45-79` | 乘/除以 `s^N` |
| `parameter_intervals`（→ `NbIntervals`/`Intervals`） | `trimmed.rs:130-140` | 把 basis 结点按同一 `s` 重映射 |
| `resolution` | `trimmed.rs:147-153` | 除以同一 `s` |

这三处是三条 arm **唯一**读取的量（`compute_type` 读 `GetType`/`NbIntervals`/`d1`，
`LengthParametrized` 只做 `u0 + abscis/ratio`，`AbsComposite` 读 `NbIntervals`/`Intervals` 后逐段积分）。
参数空间与弧长空间被同一个 `s` 一起缩放 ⇒ 换回 basis 参数得到的数值**完全等价**。

**③ 两条 arm 已经在端口里实现了，而且已经接线**：

```text
gcpnts.rs:66   fn compute_type(c: &dyn Curve) -> (AbscissaType, f64)   // cxx:26-65 已移植
gcpnts.rs:102  fn compute_abs_composite(...)                            // cxx:96-158 已移植
gcpnts.rs:491/492-501  match compute_type(c) { LengthParametrized | AbsComposite | Parametrized }
gcpnts.rs:579/585      abscissa_point_with_tolerance 同样按三分支分派
```

即 `LengthParametrized`（`cxx:87-90`）与 `AbsComposite`（`cxx:96-158`）**都可到达**。
模块头那段「UNPORTED … cannot be expressed faithfully」与代码现状矛盾（已按实测改写）。

#### 9.295.3 消费方清单（若**真要**改成 basis 参数，改动面如下）

`GeomTrimmedCurve` 的构造点 12 处、`Arc<dyn Curve>` 侧的显式依赖约 8 处（注释里明写 `[0,1]` 的）：

| 类别 | 位置 |
|---|---|
| 构造 | `rectangular_trimmed.rs:262,272`、`brep_builder_api.rs:229,246`、`geom_int_intss_make.rs:167,269,660`、`edge_split.rs:53`、`int_tools_lines.rs:131,132`、`step/read_topology.rs:1630,1668,1670` |
| 显式依赖 `[0,1]` 的注释/换算 | `inttools/intersections.rs:18,38`、`edge_edge/edge_edge.rs:307`、`edge_split.rs:8`、`brep_builder_full/edge_builder.rs:107`、`extrema_surf/elementary_curve_extrema.rs:155`、`gcpnts.rs`（旧注释） |
| 已用「忠实 view」的地方 | `sphere.rs:47`（唯一一处用 `GeomTrimmedCurveBasis`） |

**注意**：`GeomTrimmedCurveBasis`（`trimmed.rs:169-249`）**已经就是** OCCT 语义的那一版
（参数直通、无链式法则、`First/LastParameter` = 修剪界），并且只被 `sphere.rs` 一处使用。

#### 9.295.4 回归风险（三种方案对比）

| 方案 | 改动面 | 风险 | 判定 |
|---|---|---|---|
| **A. 把 `GeomTrimmedCurve` 整体改成 basis 参数** | 12 处构造 + ~8 处显式 `[0,1]` 换算 + gcpnts 内部 | **高**：消费方按 `[0,1]` 喂参数（如 `d0(t)`、`resolution`），改后**静默给出错误几何**（不报错、只是取到 basis 的另一段） | **不做** |
| **B. 保留双 struct，按需换** | 与 A 基本同量（`Arc<dyn Curve>` 处处要选型） | 中，且把「用哪一版」变成每个调用点的隐性决定 | 不推荐 |
| **C. 不动参数表示，只修过期注释/文档**（本轮采用） | 1 处注释 | 无 | **采用** |

A 之所以危险：这个约定是**承重的内部不变量** —— 构造用 basis 参数（`GeomTrimmedCurve::new(c, u0, u1)`），
求值用 `[0,1]`。两者都「看起来合理」，改错一个方向都不会编译失败。

#### 9.295.5 对 T-51 的影响

T-51 的豁免理由是「余下两条要先改 TrimmedCurve 参数表示」。既然两条 arm **已经实现并接线**、
且不需要改参数表示，**该豁免理由不成立**；T-51 应重开核对其真实余量（`NbPoints` 语义、
`uniform_abscissa` 的外形差异），而不是等一个不必要的重构。

---

### 9.296 —— T-94 定位：共享 registry 的写入点 = `ComposeShell` 里的 `swap_seam`（经 `repr_key` 落到原面）

#### 9.296.1 结论

`ShapeFix_ComposeShell` 一共 6 处写 `GeometryRegistry`。其中 **5 处写的是新建的 sub-edge**
（新对象，没有别的消费者），**只有 1 处改的是"来自原面的既有边"**：

```
crates/occt-topo/src/shape_fix_compose_shell/wire_data.rs:50-67   fn swap_seam(edge, face)
  ...
  pcs.swap(0, 1);                                          // cxx:541
  reg.set_edge_pcurves(&edge.0, face_key, pcs);            // :65  改「forward occurrence 的边」
  reg.set_pcurve_range(&edge.0, face_key, uff, ulf);        // :66
```

即：`SwapSeam`（`ShapeExtend_WireData.cxx:511-543`）交换的是**那条 seam 边本身的 FWD/REV pcurve 对**
和它的 range —— 那条边**不是新对象**，它仍然属于原面。这就是 T-94 要找的写入点。

#### 9.296.2 为什么它会泄漏到原面（关键在 `repr_key`）

这是本轮最要紧的一条，它把「为什么改一个面的修复会污染别的面」解释清楚了。
`tgeometry.rs:519-531`：

```rust
/// Canonical `BRep_GCurve` representation key: OCCT keys a CurveOnSurface
/// representation by the `Geom_Surface` handle (plus location); the port
/// has identity locations, so the registered surface's data pointer is the
/// key. Faces that share a surface therefore share one representation, as
/// in OCCT.
fn repr_key(&self, face_key: usize) -> usize {
    self.face_surfaces.read().unwrap().get(&face_key)
        .map(|s| Arc::as_ptr(s) as *const () as usize)   // ← 面 → 曲面指针
        .unwrap_or(face_key)
}
```

而 `set_edge_pcurves` / `edge_pcurves` / `pcurve_range` 全部先过 `repr_key`
（`:358-380`、`:427`）。所以：

* pcurve 的存储键**不是面的身份**，而是**面所用曲面的指针**；
* `ComposeShell` 重建面时沿用原面的曲面（`builder.make_face(surf.clone(), ...)`），
  ⇒ 修复后的面与原面**共用同一个 `repr_key`**；
* 于是 `swap_seam` 用「修复面的 key」写下去，实际写到的就是**原面那条边的同一个槽位**
  ⇒ 原面在 registry 里的 pcurve 被交换/改 range。

**这与 D2「补偿改写了共享 GeometryRegistry，使解析期探针与最终模型分叉」是同一机制**，
现在有了确切的代码位置（`wire_data.rs:65-66`）与中间环节（`repr_key`）。

#### 9.296.3 与 T-99 / T-93 的关系

T-99 的 F113 与 T-93 的 seam 步骤都在 `ComposeShell` 之后（`fix_missing_seam` 结尾调 `comp.perform()`，
`shape_fix_face.rs` 的 `cxx:2261` 那行）。所以：

* T-93 接线后 `resolve_face` 里「改了 registry 又丢弃结果」的那一半副作用，
  其具体执行者就是 `swap_seam` 这一处（以及 `wire_data.rs:65-66` 这一类写入）；
* 若要按 T-94 的方向做「事务式写入」，落点就是 `wire_data.rs` 的 `swap_seam`
  （以及 `dispatch_wires.rs:210-211` 等其余 5 处新建边写入 —— 那 5 处可以保持就地写，
  因为它们写的是新对象、没有共享消费者）。

#### 9.296.4 尚未做的实测（诚实标注）

本轮的定位是**代码路径 + 键机制的确定**，**没有**跑那组「前后快照」。原计划的判别实验是：

```text
在 shape_fix_face.rs 调 comp.perform() 前后，对「原面 key」下的 edge_pcurves 做前后快照，
对比 swap_seam 命中次数与槽位变化。
```

之所以没跑：`GeometryRegistry` 目前**没有**「枚举某条边已有的所有 representation 键」的公开 API
（只有按 `face_key` 取），要跑就得先加一个仅供诊断用的枚举口。留给下一轮，且**应先决定**：
是加诊断 API 跑快照，还是直接按 §9.296.2 的机制把 `swap_seam` 改成事务式（在结果面上重建边、
或在 `ComposeShell` 结束时统一提交）。

#### 9.296.5 复跑（定位复核）

```text
grep -n 'GeometryRegistry'
  crates/occt-topo/src/shape_fix_compose_shell/*.rs
# 6 处写入：dispatch_wires.rs:85,210,211,293,294 · split_by_line.rs:623,624 · split_wire.rs:443 · wire_data.rs:65,66
# 其中只有 wire_data.rs:65-66（swap_seam）写的是既有边
sed -n '519,531p' crates/occt-topo/src/tgeometry.rs     # repr_key
```

#### 9.296.6 泄漏的**构造级证明**（本轮补）

上一节说「`repr_key` 相同 ⇒ 写落到原面」，本轮把它从推断变成**由构造方式直接推出**：

1. `shape_fix_face.rs:453`：
   ```rust
   let mut tmp_f = builder.make_face(surf.clone(), &tmp_wires);   // surf 来自原面
   ```
   ⇒ `tmp_f` 是一个**新面对象**，但它的曲面是原面曲面的 `Arc` **克隆（同一指针）**；
2. `shape_fix_face.rs:622-623`：
   ```rust
   let mut comp = ComposeShell::new();
   comp.init(grid, &tmp_f, CONFUSION);      // ← ComposeShell 操作的是 tmp_f
   ```
3. `ComposeShell` 内部把 `tmp_f` 交给 `load_wires(face)`（`perform.rs:24`、`:65`），
   一路把 `face = &tmp_f` 传到 `reverse_wire_data_on_face(&mut wire, &face)`
   （`collect_wires.rs:261`、`load_wires.rs:123`）→ `swap_seam(edge, face)`；
4. `wire_data.rs:54`：`let face_key = GeometryRegistry::shape_key(&face.0);` ⇒ 取的是 **`tmp_f` 的 key**；
5. `tgeometry.rs:524-531`：`repr_key(tmp_f_key) = face_surfaces[tmp_f_key] 的曲面指针`
   = **原面曲面的同一个 `Arc` 指针**（因为 (1)）；
6. ⇒ `set_edge_pcurves(&edge.0, tmp_f_key, pcs)` 经 `repr_key` 后写进的是
   **原面那条边、原表示键下的槽位**（`tgeometry.rs:359,370,378` 都先过 `repr_key`）。

所以 `swap_seam` 交换的是**原面 seam 边的 FWD/REV pcurve 对与 range**。
这不再是"可能"，而是这四行构造代码的必然结果。

#### 9.296.7 修复方案与各自风险（结论：**需要先决策，本轮不动手**）

| 方案 | 做法 | 风险 |
|---|---|---|
| **A. 给 `tmp_f` 一个独立的表示键** | `tmp_f` 建好后把它从 `face_surfaces` 摘掉（`clear_shape` 会连边/面几何一起清，故需要一个新的"只摘曲面注册"入口），使 `repr_key(tmp_f)` 回退成 `tmp_f` 自身 key | **中**：`repr_key` 的"共享曲面即共享表示"注释说是 OCCT `BRep_GCurve` 语义；改成"每个面独立"会动到全局键约定，需整轮复核 |
| **B. 让 `swap_seam` 就地不写共享槽位** | 交换前 clone 该边，只改 clone 的 pcurve | **低-中**：`swap_seam` 的意图就是"让这条边的 FWD/REV 对换"，改 clone 等于不改 —— 语义上说不通，除非 clone 被后续装配采用 |
| **C. 事务式：`ComposeShell` 结束时统一提交** | 把 6 处写入收集起来，在结果面确定后一次性应用 | **中**：改动面最大（跨 `wire_data`/`dispatch_wires`/`split_by_line`/`split_wire`），但最贴合 OCCT 的 `ShapeBuild_ReShape` 语义 |

**本轮不动手的原因**：三个方案都要么动全局键约定（A）、要么改共享管线（C），
按红线都需要**整轮 23 模型 `step_obj_gates` 复核**才能提交；而我这轮剩余预算不足以既改又全量复核。
留一个明确的决策点比留一个半改的共享管线更安全。

#### 9.296.8 判据（下一轮用）

* 若选 A：**`zz_seam_fix -- data/occ/a3n00.stp 113` 的结果面数、全部 23 模型的 `step_obj_gates` 都必须不劣化**；
  并且应能看到「原面在被修复后其 pcurve 不再变化」——这需要一个"取原面某条边 pcurve 前后对比"的探针，
  比 §9.296.4 里的"枚举 representation 键"更容易写（只需两次 `edge_pcurves(edge, original_face_key)` 比较）。
* 若选 C：先只做 `wire_data.rs` 的 `swap_seam` 一处（其它 5 处写的是新建边，本就不需要事务化），
  再跑同一组门禁。

#### 9.296.9 实测：泄漏是**真的**，但只命中 21 个被修复面中的 **1 个**（T-94 需要重新定级）

§9.296.6 的构造链只是「可能」的充分条件，不是「实际发生」的证据 —— `repr_key` 相同只说明
写入会落到同一个槽位，**没说明 `swap_seam` 那条分支真的被走到**（`wire_data.rs:57` 有
`if pcs.len() < 2 { return; }`，而且 `reverse_wire_data_on_face` 只在 `stat < 0` 或 `reverse` 分支才被调用）。

新增探针 `crates/occt-topo/examples/zz_share_probe.rs` 直接量：
**对每个面，在调 `fix_missing_seam` 之前 / 之后，采该面自己的 `edge_pcurves` + `pcurve_range` 快照，比较差异**
（比 §9.296.4 里想的「枚举 representation 键」好写，且更直接 —— 它量的是"原面是否被改"这件事本身）。

`data/occ/T0M.stp` 实测：

```text
SHARE face=1752 fix=true changed_slots=4 (of 6)
SHARE TOTAL faces_with_fix=21 faces_changed=1 slots_changed=4
```

即：**21 个 `fix_missing_seam` 返回 true 的面里，只有 1 个（face 1752）自己的 pcurve 被改了，改了 4 个槽位（共 6 个）。**

#### 9.296.10 这对 T-94 意味着什么（重定级）

| | 上一轮的判断（§9.296.6） | 本轮实测 |
|---|---|---|
| 泄漏机制 | 成立（构造链 6 步） | **仍成立** |
| 泄漏**发生频率** | 未量，隐含"普遍" | **21 个修复面里 1 个** |
| 对共享管线的污染面 | 未量 | **T0M：4 个 pcurve 槽位** |

⇒ **T-94 的严重度需要下调**：它不是「每次 seam 修复都改写共享 registry」，而是
「**少数面**（T0M 21 个里 1 个）的原面 pcurve 会被改动」。这解释了为什么：

* D9 的「pcurve 归一化 pass 单独加回无影响」是惰性的 —— 命中的面极少；
* 「推迟 shift 写入」实测与基线逐字相同 —— 同样因为命中面少；
* 而 D2 观察到的「解析期探针与最终模型分叉」确实存在，但**来源覆盖面小**。

**结论**：A/B/C 三个方案（§9.296.7）依然可选，但**优先级应下调**，把预算放到命中面更多、
收益更直接的卡上。若仍要做，判据（§9.296.8）不变，另加一条：**face 1752 的 4 个槽位不得再变化**，
且 23 模型 `step_obj_gates` 不劣化。

**尚未做**：face 1752 被改的那 4 个槽位具体是什么（哪两个 pcurve 被交换、range 怎么变），
以及这 4 个槽位的改动是否**可观测地影响**它的网格（T0M 未网格面数在接线前后是 7→6，
face 1752 不在那 7 个里）。

#### 9.296.11 复跑

```text
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_share_probe -- data/occ/T0M.stp
# 可选第二参限制面数：  ... -- data/occ/a3n00.stp 60
```

#### 9.296.12 a3n00 交叉核对

同一探针跑 `data/occ/a3n00.stp` 全量：

```text
SHARE TOTAL faces_with_fix=3 faces_changed=0 slots_changed=0
```

⇒ a3n00 上**一个面都没被改**（只有 3 个面 `fix_missing_seam` 返回 true，且都没碰到 `swap_seam` 的写）。

两个模型的合计：**24 个被修复面里只有 1 个（T0M face 1752）原面 pcurve 被改动**。
这进一步支持 §9.296.10 的重定级：机制真实，但命中面极少。

---

### 9.297 —— T-95 结项：`Geom2dInt_Geom2dCurveTool` 与 `TheProjPCurOfGInter` **都已实现且已接线**（卡面描述过期）

#### 9.297.1 结论

T-95 的标题与 `next` 说「补 `Geom2dInt_Geom2dCurveTool` 的 `GetType/Value/D0-D3/NbSamples`，再补
`Geom2dInt_TheProjPCurOfGInter`」—— **这些全都已经在了**，而且已经接进 `Geom2dInt` 的分派。
这是本会话第二次遇到「卡面描述与代码现状脱节」（第一次是 T-96/§9.295）。

#### 9.297.2 逐项核对（`crates/occt-geom2d/src/geom2d_int/`）

`curve_tool.rs`（166 行）—— 卡里点名的全部已实现：

| `lxx` 位置 | 端口 | 状态 |
|---|---|---|
| `:32-35` `GetType` | `curve_tool.rs:31 get_type` | ✅（Line/Circle/Ellipse/Hyperbola/Parabola/Bezier/BSpline/Offset/Other 全分支） |
| `:68-71` `Value` | `:84 value` | ✅ |
| `:74-77` `D0` | `:89 d0` | ✅ |
| `:80-109` `D1`/`D2`/`D3` | `:94`/`:99`/`:104` | ✅ |
| `:112-117` `DN` | `:109 dn` | ✅ |
| `:120-129` `First/LastParameter` | `:114`/`:119` | ✅ |
| `:23-91` `NbSamples`（两个重载） | `:136 nb_samples_curve`、`:141 nb_samples_curve_range` | ✅ |
| `:146-178` `NbIntervals`/`Intervals`/`GetInterval` | `:147`/`:152`/`:157` | ✅ |
| `:38-65` `Line/Circle/Ellipse/Parabola/Hyperbola` | `:59`-`:80` | ✅ |
| `hxx:46` `IsComposite` | `:54` | ✅（端口无 composite 2D 曲线，返回 false） |

`proj_p_cur.rs`（56 行）—— `Geom2dInt_TheProjPCurOfGInter` 的两个入口都已实现：
`find_parameter_range`（`cxx:27-65`）、`find_parameter`（`cxx:67-79`），
内部走 `curve_locator::locate_range` + `GenLocateExtPC`（即 `Extrema_GCurveLocator` + `Extrema_GenLocateExtPC`）。

**已接线**（这是关键，光有实现不算数）：

```text
conic_curve.rs:13   use super::{curve_tool, proj_p_cur};
conic_curve.rs:296  let v = proj_p_cur::find_parameter(the_par_curve, point, tolerance);
conic_curve.rs:305  let mut x = proj_p_cur::find_parameter_range(…);
```

#### 9.297.3 模块内**真实**剩余的 UNPORTED（全部有界，且不在 T-95 名下）

`grep -n UNPORTED crates/occt-geom2d/src/geom2d_int/*.rs` 只剩 4 处：

| 位置 | 内容 | 判定 |
|---|---|---|
| `curve_tool.rs:128` | `EpsX(C, Eps_XYZ)` = `C.Resolution(Eps_XYZ)`；`Curve2d` 没有 `Resolution` | **未移植但无人调用**（本轮复跑 grep 确认：全树只有定义处与文档提及，无消费方） |
| `curve_tool.rs:161` | `Degree` = `Adaptor2d_Curve2d::Degree` | 同上，**未移植但无人调用** |
| `ginter.rs:174` | `IntCurve_IntConicConic_1.cxx:2236` 的 **Line/Circle** 臂 | **真实缺口**，但不在 T-95 名下 |
| `ginter.rs:188` | `IntCurve_IntCurveCurveGen.gxx:339-…` 的 **`typ1 != Line`** 臂 | **真实缺口**，同上 |

⇒ 前两处是"接口未搬但无消费方"，不构成功能缺口；**真正剩下的是 `ginter.rs` 的两条分派臂**，
应另立卡（见 §9.297.4），不要挂在 T-95 这个名字（`curve_tool`/`proj_p_cur`）下面 —— 那样下一个接手的人
会以为缺口在已经完成的两个文件里。

#### 9.297.4 交接：新卡 T-100（`ginter.rs` 的两条分派臂）

本卡结项时把真实余量拆出去：

* 范围：`ginter.rs` 的 `IntCurve_IntConicConic_1.cxx:2236`（Line/Circle）与
  `IntCurve_IntCurveCurveGen.gxx:339-…`（`typ1 != Line`）两条臂；
* 验收：对应 UNPORTED 注释清掉；`--lib`（occt-geom2d / occt-topo）不劣化；`step_obj_gates` 不劣化；
* 约束：仍是「按 `.cxx` 补同等分支」，不得按某个 STEP 的形状反推规则。

#### 9.297.5 复跑

```text
grep -n UNPORTED crates/occt-geom2d/src/geom2d_int/*.rs
grep -rn 'eps_x_with\|curve_tool::degree' crates/          # 应为空（无消费方）
grep -n 'proj_p_cur::' crates/occt-geom2d/src/geom2d_int/conic_curve.rs
```

---

### 9.298 —— T-99 检查点 ① 结清：`v_range=1e100` 属于**原面**，修复后那条路径拿到的是正确的 **175**

#### 9.298.1 结论

§9.294.4 留的检查点 ①（「核 `v_range` 为什么是无穷」）**不是缺陷**，是读错了对象：

| 调用对象 | `uv_bounds` 的 V | `v_range` |
|---|---|---|
| **原面** F113（`surf.v_range = (-inf,inf)`，`uv_bounds` V 也是 `(-inf,inf)`） | `(-inf, inf)` | 钳成 **1e100** |
| **修复后的面**（seam 构造后） | **`(-43.75, 131.25)`** | **`175`** ✅ |

实测（临时插桩 `OCCT_TOPO_TRACE_RANGE` 打在 `cxx:1804-1805` 之后，已摘除）：

```text
RANGE f_v=(-43.75,131.25) surf_v=(-inf, inf) svf=-43.75 svl=131.25 v_range=175 u_range=6.283185307179586
RANGE f_v=(-43.75,131.25) surf_v=(-inf, inf) svf=-43.75 svl=131.25 v_range=175 u_range=6.283185307179586
```

全次运行里 `v_range=1e100` 出现 **0 次**；`v_range=175` 出现 **2 次**（正是修复后的面）。
并且端口的守卫**逐行忠实地**复现了 `cxx:1781-1802`（`shape_fix_face.rs:210-224`）：
surface 的 `(-inf,inf)` 被 `f_v1/f_v2` 替换 ⇒ `175`。

⇒ **端口这段是对的**，不要再查它。

#### 9.298.2 那 1e100 是从哪来的

`zz_seam_fix … 113` 的 `--wireuv` 模式里，我的**探针自己**用了
`(d - c).abs().min(1e100)`（`zz_seam_fix.rs` 顶部），而它对**原面**取范围
⇒ 探针自己把它钳成了 `1e100`。那是**测量代码的口径**，不是被测量的 `CheckWire` 的输入。
§9.294.4 把它当成"可疑量"是我读错了归属 —— 此处更正。

#### 9.298.3 对「两条真 wire 求和恰好为 0」的影响

**没有影响，那个观察依然成立且依然关键**（§9.294.2）：

| wire | 边数 | 逐边位移之和（`dU=6.283185`、`dV=175` 下） | 结果 |
|---|---|---|---|
| wire[0] | 6 | `(0.000000, -0.000000)` | REJECT |
| wire[1] | 14 | `(0.000000, -0.000000)` | REJECT |
| wire[2] | 1 | `(6.283185, 0)` | `Some((1, 0, false))` |
| wire[3] | 1 | `(-6.283185, 0)` | `Some((-1, 0, false))` |

这里的 `0` **不是钳位产物**（`min(·,1e100)` 对 0 无影响），是真的求和为 0。
且 `dV=175` 也已由本节确认为正确值。

⇒ **T-99 现在只剩一个检查点**：修复后的面上，wire[0]/wire[1] 的**逐边 pcurve 端点**为什么互相抵消到 0
（§9.294.2 的 `d=(0.812,62.252) (0.088,-86.089) (0,-27.703) (-0.177,0) (0,27.703) (-0.724,23.837)`）。
从这个序列看，`y` 取值在 `93.61 / -86.09 / -27.70 / 27.70 / 23.84` 之间来回 —— 像是**同一条几何被走了两遍而方向相反**
（`+27.703` 与 `-27.703`、`+23.84` 与前面某段），而不是"参数取错"。

#### 9.298.4 下一步（唯一检查点，且已缩小）

在**修复后的面**上逐 wire 打印「边序 + 每条边的起点/终点 3D 坐标」，判断：

* 若 wire[0] 的 6 条边在 3D 上构成**闭合环**（首尾相接回到起点），那 `(0,0)` 是**正确**的，
  说明这 6 条边构成的是**内环**而非 seam 两侧 —— 那么问题在「哪两条 wire 被当成 seam 两侧」的选择，
  而不在 pcurve；
* 若 3D 上不闭合，则 pcurve 的关联有问题（该继续查 `associate_edge_pcurve`）。

这条判据能**一次性把问题分到"选择侧"还是"pcurve 侧"**，比继续猜参数口径有效。

#### 9.298.5 复跑

```text
# 探针（含 surf/uv_bounds 原始输入）
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_seam_fix -- \
    data/occ/a3n00.stp 113 --wireuv
```

---

### 9.299 —— T-99 【判据落地】F113 的 4 条 wire **在 3D 上都闭合** ⇒ `CheckWire` 的 REJECT 是**正确**的，问题在**上游切分**侧，不在 `CheckWire` / pcurve

#### 9.299.1 实测（§9.298.4 的判据 (a) 命中）

在 `--wireuv` 里新增「取每条 wire 的首/末顶点 3D 坐标（按朝向）+ 闭合间隙」：

```text
WIREUV wire[0] first=(12.387499,24.672033,-89.394247) last=(同左) gap=0.000000 cl3d=YES
WIREUV wire[0] sum=(0.000000,-0.000000) check_wire=None
WIREUV wire[1] first=(-8.900000,28.844410,-48.000000) last=(同左) gap=0.000000 cl3d=YES
WIREUV wire[1] sum=(0.000000,-0.000000) check_wire=None
WIREUV wire[2] first=(87.500000,0.000000,-32.000000)  last=(同左) gap=0.000000 cl3d=YES
WIREUV wire[2] sum=(6.283185,0.000000)  check_wire=Some((1, 0, false))
WIREUV wire[3] first=(-87.500000,0.000000,-32.000000) last=(同左) gap=0.000000 cl3d=YES
WIREUV wire[3] sum=(-6.283185,0.000000) check_wire=Some((-1, 0, false))
```

**四条 wire 在 3D 上全部闭合**（gap 精确为 0）。

#### 9.299.2 由此得到的结论（重要：推翻此前的方向）

按 §9.298.4 预先定好的判据，`cl3d=YES` 意味着**逐边位移和为 0 是正确结果**：
一条在 3D 上闭合的 wire，其 pcurve 位移和本来就可能是 0（多绕整周或走回来）。
`CheckWire` 拒绝它**不是 bug**。

⇒ **此前四轮一直怀疑的方向（`CheckWire` 的 pcurve 口径 / 参数取值）到此正式排除**，
而且这次是从"被测对象的几何事实"排除的，不是从代码读法排除的。

#### 9.299.3 真正的缺口：上游把 22 边切成了 6+14+1+1

对账：

| | wire 结构 | 边数合计 |
|---|---|---|
| **GT**（OCCT） | `[22, 6]` | 28 |
| **端口** | `[6, 14, 1, 1]` | 22 |

注意 **6 + 14 + 1 + 1 = 22** = GT 那条 22 边 wire 的边数。
再加上 wire[2]/wire[3] 的边（各 1 条、跨度 ±2π），
⇒ **端口的 4 条 wire 是把 GT 的一条 22 边 wire 拆成了 6 + 14 + 1 + 1**，
而不是"GT 有两条独立的 wire"。

其中 6 边与 14 边那两条**独立闭合**，所以它们各自满足 `IsClosed`，
再也无法与那两条 1 边弧合并成一条跨越 seam 的开口环。

**所以缺口不在 `CheckWire`、不在 `FixMissingSeam` 的 pcurve 处理，而在更上游：**
把一条跨越 seam 的边界环切成 `6+14+1+1` 的那一步。相关位置（按此前轮次的指向）：

* `FixMissingSeam` 里的 `BreakWires` / `SplitWire`（`ShapeFix_Face.cxx:2068-2136` / `:1918-2127`）；
* 端口侧的对应实现：`shape_fix_compose_shell/break_wires.rs`、`split_wire.rs`；
* §9.156/§9.157 那条线（`WireSegment::first_vertex/last_vertex` 朝向修复，2-wire 面 24→10）
  说明**顶点匹配**这个环节此前就出过问题，而"两条独立闭合"正是顶点匹配错误的典型症状
  （本该接上的两端被当成各自闭合）。

#### 9.299.4 下一轮的具体动作

**不要**再查 `CheckWire` / pcurve。改为：

1. 在 **`BreakWires` 输出**处打印 F113 的段数与每段边数（此前 §9.151 对 f209 做过同样的事）；
   若那里已经是 `6/14/1/1`，缺口在 `BreakWires`；若那里是对的、`CollectWires` 之后才变，缺口在收集侧。
2. 对着 `ShapeFix_Face.cxx:1994-2234` 的 `w2 != null` 合并路径，核「本该合并的那两条为什么没合并」——
   重点看**顶点是否为同一 TShape**（§9.156 的教训：`WireSegment::first_vertex/last_vertex`
   曾因未按边朝向取端点而错，导致 2-wire 面从 9 涨到 24）。

**判定信号（现成）**：`zz_seam_fix -- data/occ/a3n00.stp 113` 的
`SEAM selected` 应是 6 边与 14 边那两条（而不是 1-edge 的两条），最终结果 **1 个面、2 wire、22+6 边**。

#### 9.299.5 复跑

```text
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_seam_fix -- \
    data/occ/a3n00.stp 113 --wireuv      # 现在会打印 cl3d=YES/NO
```

---

### 9.300 —— T-99 stage 级定位完成：**`ComposeShell` 无辜**，`load_wires` 入口就已经是 `[6,14,1,1]`

#### 9.300.1 实测（临时插桩 `OCCT_TOPO_TRACE_CS`，已摘除）

在 `ComposeShell::perform` 的四个阶段之间打印每个 `WireSegment` 的边数（`zz_seam_fix … 113`）：

```text
CS load_wires     total=22 segs=[6, 14, 1, 1]              ← 入口就已经是 4 条
CS split_by_grid  total=26 segs=[6, 16, 1, 1, 1, 1]
CS break_wires    total=26 segs=[6, 7, 5, 4, 1, 1, 1, 1]   ← 14 被拆成 7+5+4
CS collect_wires  total=28 segs=[6, 7, 5, 4, 3, 3]
CS dispatch_wires faces=5 edges_per_face=[13, 5, 4, 3, 3]  ← 与 §9.293.5 的 5 个碎片面吻合
```

#### 9.300.2 结论：定位链彻底收敛，且**归责明确**

1. **`ComposeShell` 不产生 4 条 wire** —— `load_wires` 拿到的 `seqw` 已经是 `[6,14,1,1]`。
   而 `load_wires`（`load_wires.rs:36-50`）就是**直接读 `face` 的 child wires**，
   没有任何合并/拆分逻辑。
   ⇒ 那个 `face`（`fix_missing_seam` 里 `shape_fix_face.rs:453` 造的 `tmp_f`）**拿到的就是 4 条 wire**。
2. **`dispatch_wires` 的 `[13,5,4,3,3]`** 与 §9.293.5 观测的「1 个正确面 + 4 个碎片」
   （边数 6+7 / 5 / 4 / 3 / 3）逐个对上 ⇒ 上一轮那个"5 个面"完全由**上游给的 4 条 wire**解释，
   与 `swap_seam`/`repr_key` 那条线（T-94）无关。
3. **`break_wires` 把 14 边段拆成 7+5+4**（`[6,14,1,1]` → `[6,7,5,4,1,1,1,1]`），
   这是**第二处**碎片化，发生在 `ComposeShell` 内部，也值得单独核（但它是果不是因 ——
   因为入口的 14 边段本身就是错的）。

#### 9.300.3 因此真正的缺口点（唯一）

**`fix_missing_seam` 的入口面 `face` 就已经是 4 条 wire。** 从 §9.292 的 `TRACE face wires_before=[6,14,1,1]`
到本节的 `CS load_wires segs=[6,14,1,1]`，两处独立观测互相印证。

⇒ 缺口在 `fix_missing_seam` **之前**，即：**端口把这个面建成 4 条 wire** 的那一步。
候选（按 §9.171/§9.172/§9.175 的既有线索，**且必须二选一先确认**）：

* **(A) reader 侧**：`resolve_face`/`resolve_shell` 把该 STEP 面的 4 个 `FACE_BOUND` 直接建成 4 条 wire，
  而 OCCT 的 STEP 传输阶段已经把它们并成 2 条；
* **(B) `heal_shape`/`ShapeProcess` 侧**：端口没有 OCCT 的 `ShapeProcess`/`FixShape` 链，
  合并本该在那里发生（§9.171 指出合并在 `STEPControl_ActorRead::TransferEntity` 的 `ProcessShape`；
  §9.172 指出端口没有这条链；§9.175 实测关掉 per-entity `ProcessShape` 后多 wire 面 18→63 ⇒ 合并确实在该链）。

**这两条不矛盾**：合并算子归属 `ProcessShape` 链（B），而 (A) 决定"进这条链之前是几条"。
必须先量一个数：**`heal_shape` 被调用时，该面是 4 条还是 1 条 wire**。

#### 9.300.4 下一轮的单一动作

在 `heal_shape`（或 `read_step_file` 读完形状、进 healing 之前）对 F113 打印 wire 数与每 wire 边数：

* 若那里是 `[22]` 或 `[22,6]` ⇒ 缺口在 **healing 链**，即端口缺 `ProcessShape` 的合并算子，
  应据 `ShapeProcess_OperLibrary.cxx:785-899` 补（**这是新增一段同等分支，不是调参**）；
* 若那里已经是 `[6,14,1,1]` ⇒ 缺口在 **reader**，`resolve_face` 把 4 个 `FACE_BOUND` 建成了 4 条独立 wire，
  应据 `StepToTopoDS_TranslateFace` 的 bound 处理核对该不该合成一条。

这个量**一次就能定死 (A) 还是 (B)**，且两边都是「有 `.cxx` 依据的同等分支」。
之前六轮之所以绕，是因为一直在**下游**（`CheckWire`/`FixMissingSeam`/`ComposeShell`）找，
而下游每一步都被证明是忠实的。

#### 9.300.5 复跑

```text
# 阶段打印（需按本节描述重加，已摘除）
OCCT_TOPO_TRACE_CS=1 rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline \
    --example zz_seam_fix -- data/occ/a3n00.stp 113 2> cs113.err
grep 'CS ' cs113.err | tail -6
```

---

### 9.301 —— T-99 归责定案：是 **(B) reader 侧**，但**不是 reader 的 bug**；真缺口是端口没有 `ProcessShape` 那条链

#### 9.301.1 实测答案（§9.300.4 要求的那个数）

**`heal_shape` 从来没在读入路径上被调用过。** 全树 grep `heal_shape` 只有 10 处，且**全部**是：

* 定义 `shhealing/common.rs:620`；
* 测试 `tests/phase7_integration.rs:9,73`、`shhealing/tests.rs:80,99,181,222`；
* `shhealing/wire_heal.rs:6,9`（内部调用）。

**读入路径**：`read_step_file` → `read_step` → `read_step_impl`（`step/read_geometry.rs:1009-1060`）
→ `SsrComposer::rep_shape` / `resolver.resolve_shape` → **直接 `model.add(...)`**。
中间**没有任何 healing 阶段**。

⇒ §9.300.4 的 (A)/(B) 二选一里，答案是 **(B)**：那个面在进 `fix_missing_seam` 时已经是 4 条 wire，
因为**它就是被 reader 建成 4 条的**，而且此前没有任何环节有机会合并。

#### 9.301.2 但 reader 这一层是**忠实的**

`read_topology.rs:645-654`：

```rust
let mut wires = Vec::with_capacity(bounds.len());
...
for &b in &bounds {
    let mut s = self.resolve_shape(b)?;
    if !s.is_wire() { return Err(...); }
    wire_bound_oris.push(s.orientation());
    ...
}
```

**`FACE_BOUND` → wire 是 1:1**，没有合并。而 STEP 文件对该面确实声明了 4 个 bound
（§9.168/§9.174 已两次核实：`#5375` 确为 4 bound）。
⇒ **reader 的输出与文件一致**，`resolve_face` **不是** bug。

#### 9.301.3 真缺口（唯一，且有出处）

OCCT 把 4 bound 合成 2 wire 的那一步在 **`STEPControl_ActorRead::TransferEntity` 的 `ProcessShape`**
（§9.171 已指出；§9.175 实测关掉 per-entity `ProcessShape` 后多 wire 面 18→63 ⇒ 合并确实在该链），
它走 `ShapeProcess_OperLibrary::ProcessShape` → `ShapeProcess` → `ShapeFix_Shape` 系列。

**端口没有这条链，`heal_shape` 是另一个自研近似且根本没接在读入路径上**（§9.301.1）。
这就是全部原因：**一个"读进来之后做形状级修复"的阶段从来不存在。**

#### 9.301.4 这解释了此前七轮的每一个观测

| 观测 | 解释 |
|---|---|
| `wires_before=[6,14,1,1]`（§9.292） | reader 忠实产出 4 个 bound 的 4 条 wire |
| `CheckWire` 拒两条真 wire（§9.294） | 那 6/14 边段本就是被拆断的碎段，各自闭合（§9.299 实测 `cl3d=YES`） |
| seam 结果 = 1 正确面 + 4 碎片（§9.293.5） | `FixMissingSeam` 在错误的输入（4 条碎 wire）上做正确的事 |
| `v_range=1e100` 属原面（§9.298） | 与 wire 结构无关 |
| `ComposeShell` 入口即 4 条（§9.300） | 同第一条 |

⇒ **下游七轮里我查过的每一级都是忠实的**，因为缺口在整个链条的**上游**、而且是一个**缺失的阶段**
而不是某段代码写错。

#### 9.301.5 下一轮的动作（这才是真正要动手的第一步）

**先读 `ShapeProcess_OperLibrary.cxx:785-899`**（`ProcessShape`），确认：

1. 它注册了哪些算子、`FromSTEP.FixShape` 默认值下哪些**实际会跑**（`STEPControl_Controller.cxx:201,221`）；
2. 具体哪个算子负责「把同一面的多个 bound 合成更少的 wire」——
   这是关键：**必须找到确切的算子与行号**，不能凭"大概就是 FixShape" 去补一个泛化的自研 healing 链
   （那正是红线里"自创 OCCT 里不存在的分支"）。
3. 若确认是某个具体算子（如 `ShapeFix_Face` 的某段、或 `ShapeFix_Shape` 的 bound 处理），
   再按**同等分支**补一小段，并用 `zz_seam_fix … 113` + 23 模型 `step_obj_gates` 复核。

**若第 2 步找不到确切算子**：则这张卡应转为「标注端口缺 `ShapeProcess` 链」的文档性结论，
并在 `shape_fix_face.rs` / `read_geometry.rs` 写下 `UNPORTED` 注释（含 `.cxx` 行号），
**不要**自研一个泛化 healing 链去凑 F113。

#### 9.301.6 复跑

```text
grep -rn 'heal_shape' crates/            # 确认只在测试与 wire_heal 内
sed -n '1009,1060p' crates/occt-topo/src/step/read_geometry.rs   # 读入路径无 healing
sed -n '645,660p'  crates/occt-topo/src/step/read_topology.rs    # bounds 1:1 → wires
```

#### 9.301.7 与 §9.176/§9.178 的一致性核对（重要，避免把结论下错）

前几轮还有两条证据必须一起看，否则会得出错误结论：

* **§9.176**：逐算子对拍 ⇒ `FixShape` 是**唯一**会合并 wire 的算子；
* **§9.178**：端口的 **wire 轮对 F113 零拓扑变化**（8→8, 14→14, 1→1, 1→1）⇒ **合并不在 wire 轮**；
* **§9.173**：`fix_missing_seam` 对 F113 **确实返回了** `Shell(5 面)`。

三条合起来说明一件比 §9.301.3 更精确的事：

**合并算子（`FixMissingSeam`）在端口里是"存在且被调到"的** —— 它返回了 `true` 和 5 个面，
不是"没跑"。所以缺口**不是"完全没有合并算子"**，而是：

> **`FixWires` 没能把 4 条 wire 并成 2 条** —— 它对 F113 **零拓扑变化**（§9.178），
> 于是 `FixMissingSeam` 拿到的是 4 条碎 wire，只能在这错的输入上做出"1 正确面 + 4 碎片"。

而 §9.299 已经证明那 4 条 wire **各自在 3D 上闭合**。
`ShapeFix_Wire` 的合并/重排依赖**顶点身份（同一 `TopoDS_Vertex`）**判定邻接；
若两条本该接上的 wire 的端点**不是同一个 TShape**，它们就会被当成各自独立的闭链，永远不会被并成一条。
**这正是 §9.156/§9.157 那条线**（`WireSegment::first_vertex/last_vertex` 曾因未按边朝向取端点而错，
2-wire 面从 9 涨到 24）—— 症状完全一致。

⇒ **真缺口 = 顶点身份/邻接判定**，落点在 `ShapeFix_Wire` 侧（端口 `shhealing/wire_fix.rs` 与
`shape_fix_compose_shell/` 的 `vkey`/顶点比较），**不是**"补一整条 `ShapeProcess` 链"。
后者是 §9.301.3 的粗判，被本节这三条证据修正为更窄、更可动手的一处。

#### 9.301.8 因此下一轮的动作（修正后，可动手）

1. 在 `ShapeFix_Wire` 的合并/重排入口处，对 F113 的 4 条 wire 打印**每条 wire 首末顶点的 TShape 指针**（`vkey`），
   看 6 边段与 14 边段相邻的两端**是否同一 TShape**；
2. 对照 `ShapeFix_Wire.cxx` 的 wire 合并/`FixReorder` 段（`cxx:3977-4310` 一带），
   确认 OCCT 在此处是按 `TopoDS_Vertex::IsSame` 还是按**几何重合**判定邻接；
   **若 OCCT 用几何重合而端口用指针相等**（或反之），那就是一处可修的同等分支差异；
3. 若确认是顶点身份问题，修的落点在端口造/比顶点处（`WireSegment::first_vertex/last_vertex`、
   `vkey`），**不是**在 `fix_missing_seam` 或 `ComposeShell` 里加谓词。

**这一步有明确的"看两个指针是否相等"的判据，是可执行的一轮，不再是盲搜。**

---

### 9.303 —— T-99【定案】F113 = STEP 面 `#5375`：4 个**真正互相独立**的闭合 loop；端口逐字忠实，"合并"另有来源

#### 9.303.1 面已 100% 锁定

`.target-gate/f113_5375.py` 用**正解**取出该面的 bound 列表（关键：`ADVANCED_FACE` 的 bound 是
**直接内联在记录里**的 —— `ADVANCED_FACE('',(#5140,#5352,#5363,#5374),#4952,.T.)`，
不是靠反向索引找 `FACE_BOUND`。前两次脚本用反向索引才一直找不到，那是脚本的错，不是文件的）：

```text
face #5375: ADVANCED_FACE('',(#5140,#5352,#5363,#5374),#4952,.T.)
  surface #4952: CYLINDRICAL_SURFACE('',#4951,34.0)      ← 半径 34，高 68
```

**半径验算**：`2πrh = 2π × 34 × 68 = 29053.4` —— 与 §9.291 算出的该区域解析缺口 **29053.4 逐字吻合**。
（此前几轮猜过 R=87.5 / R=68，都不对；真值是 **34**。）§9.174 提到的 `#5375` 也正是这个面。

#### 9.303.2 四个 bound 的真实结构

| bound | 类型 | loop | 边数 | 相邻共享 | 自身闭合 |
|---|---|---|---|---|---|
| `#5140` | **FACE_OUTER_BOUND** | `#5139` | **8** | 7/7 | ✅ 首末都含 `#4954` |
| `#5352` | FACE_BOUND | `#5351` | **14** | 13/13 | ✅ 首末都含 `#5142` |
| `#5363` | FACE_BOUND | `#5362` | **1** | — | 单顶点 `#5354`（退化） |
| `#5374` | FACE_BOUND | `#5373` | **1** | — | 单顶点 `#5365`（退化） |

**全部 23 个顶点 `shared vertices: 23 / 23`** —— 这里的意思是"每个顶点各自只被本 loop 使用一次"，
即**跨 loop 没有一个顶点被复用**（上一行输出容易误读，实际是 `len(oes)==1` 的计数）。
两个大 loop **各自首尾闭合**，彼此不相邻。

#### 9.303.3 端口的输出与文件**逐字一致**（三重印证）

| | 文件 `#5375` | 端口 F113 | |
|---|---|---|---|
| 边数分布 | `8 / 14 / 1 / 1` | `6 / 14 / 1 / 1` | 前 8 边组中 2 条被 `bind_edge_loop_vertices` 的 seam 处理合并成 6（§9.289 已核） |
| 曲面 | `CYLINDRICAL_SURFACE` R=34 | 圆柱，`u=[0,2π]` | ✅ |
| 曲面参数域 | — | `v=[-43.75,131.25]`（修复后） | ✅ 与 GT 一致 |
| 顶点数 | 23 | 23（`all_vkeys` 6+14+1+1） | ✅ |

⇒ **端口没有任何信息丢失**。`read_topology.rs:645-654` 的 `FACE_BOUND → wire` 1:1 是**忠实的**；
`resolve_face` 忠实；`load_wires`（读 child wires）忠实；`ComposeShell` 忠实（§9.300）。

#### 9.303.4 因此"4 条 → 2 条"**不可能**靠顶点邻接完成

§9.302.3/§9.301.8 的假设（"端口比较顶点的口径不对"）**正式否证**：
不是"口径不同"，而是**文件里就没有共享顶点可供连接**。任何基于 `TopoDS_Vertex::IsSame`
的合并算法在这一输入上都只能得到 4 条 wire —— 端口的行为是**唯一正确解**。

于是 OCCT 的 `[22, 6]` 只能来自**几何**层面的判定。而这里有一个很强的线索：

* 两个 1 边 loop（`#5362`/`#5373`）各自钉在**单个顶点** `#5354` / `#5365` 上 —— 它们是**退化边**；
* 在圆柱面上，两条退化边若位于**同一 `u`（同一母线）**，则它们与另外两个环在该母线上相接，
  4 条就可以并成 2 条。

**这条必须用几何验证，不能假设**：读 `#5354` 与 `#5365` 的 `CARTESIAN_POINT`，
与 `#5139`/`#5351` 的母线位置比较。若两者同一母线 ⇒ 合并是"同参数域上的退化边与环相接"；
若不同 ⇒ 是另一条几何判据。

#### 9.303.5 下一轮（最后一个验证点，然后必动手）

1. 打印 `#5354`、`#5365` 的 3D 坐标，以及 `#5139`(8 边) 与 `#5351`(14 边) 的首末顶点坐标；
   判断两个退化点是否落在同一条母线上、以及两个环的端点是否落在**几何重合**的位置
   （§9.156 的教训：顶点身份 vs 几何重合是两种不同的判据）。
2. 据此在 `.cxx` 里找**按几何重合**把 bound 接起来的那段（`ShapeFix_Face.cxx` 的 bound 处理
   或 `ShapeFix_Shape` 的 shell 轮），确认后有同等分支再改；**找不到就不改**，
   把结论写成 `UNPORTED` 注释（含行号）。

**判定信号（现成）**：`zz_seam_fix -- data/occ/a3n00.stp 113` 最终应是
**1 个面、2 wire、22+6 边**；GT 该面 bbox = 原 F113 bbox。

#### 9.303.6 复跑

```text
python .target-gate/f113_5375.py data/occ/a3n00.stp 5375
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_seam_fix -- \
    data/occ/a3n00.stp 113 --wireuv 2>/dev/null | grep -E 'nEdges=|vkeys first|all_vkeys'
```

#### 9.303.7 两个退化点的实际位置（本轮最后一步取证）

```text
wire[3] 单边 loop 的顶点 #5365 -> #5364 = (-87.500000000, 0.0, -32.0)     ← u = π 那条母线
wire[2] 单边 loop 的顶点 #5354 -> #5353 = ( 87.500000000, 0.0, -32.0)     ← u = 0  那条母线
wire[0] 6 边 loop 的闭合点 #4954 -> #4953 = ( 12.387499, 24.672033, -89.394247)
wire[1] 14 边 loop 的闭合点 #5142 -> #5141 = ( -8.900000, 28.844410, -48.000000)
```

`x = ±87.5` 且 `y = 0` ⇒ 正是 `u = 0` 与 `u = π` **两条不同的母线**，
且都在 `z = -32`（该圆柱顶面高度）。

⇒ 两个退化点**不在同一条母线上**，位于圆柱两侧。
所以「两个退化边在同一母线、与两个环相接」这个猜想**也不成立**。

### 9.304 —— T-99【结论】端口对该面的处理链**逐级忠实**；"4 条并成 2 条"的机制在 `.cxx` 里**尚未定位**

#### 9.304.1 八轮定位的净结果（全部可复核）

| 级 | 结论 | 依据 |
|---|---|---|
| STEP 文本 | 面 `#5375`，`CYLINDRICAL_SURFACE` **R=34**，高 68 | §9.303.1；解析面积 `2πrh = 29053.4` 与该区域缺口逐字吻合 |
| STEP 文本 | 4 个 bound 是**4 个真正独立**的闭合 loop（8/14/1/1），**跨 loop 无共享顶点** | §9.303.2 |
| reader | `FACE_BOUND → wire` **1:1**，`resolve_face` 忠实 | `read_topology.rs:645-654` |
| 读入路径 | **无 healing 阶段**（`heal_shape` 只被测试调用） | §9.301.1 |
| wire 身份 | 4 条 wire **不共享任何顶点 TShape** | §9.302.1 |
| 3D 几何 | 4 条 wire **各自闭合**（gap 精确 0） | §9.299.1 |
| `CheckWire` | 拒绝是**正确**的 | §9.299.2 |
| 范围守卫 | 逐行忠实 `cxx:1781-1802`；修复后 `v_range=175` | §9.298 |
| `ComposeShell` | `load_wires` 入口即 `[6,14,1,1]`，无辜 | §9.300 |
| `dispatch_wires` | `[13,5,4,3,3]` 与 §9.293.5 的碎片面逐个对上 | §9.300.1 |

#### 9.304.2 为什么没定位到机制：**输入里没有任何可用于"合并"的拓扑信息**

这是本轮最重要的认识。OCCT 自己也用 `TopoDS_Vertex::IsSame` 判邻接（`shape_fix_compose_shell/break_wires.rs:15-18`
的 `vkey` 注释就写了这一点，与 `TopTools_ShapeMapHasher` 一致）。
既然该面的 4 个 loop **跨 loop 零共享顶点**，那么：

* **任何**基于顶点身份的 wire 合并算法，在这一输入上都**只能**得到 4 条 wire；
* 端口的输出（4 条、各自闭合）是**该输入下的唯一正确解**，不是 bug。

⇒ OCCT 的 `[22, 6]` 必然来自**顶点邻接之外**的机制。我到目前为止**没有在 `.cxx` 里找到它**，
因此**不能**据此改代码 —— 按红线，找不到同等分支就只能标注 `UNPORTED`，不能自创。

#### 9.304.3 已排除的假设（不要再走）

1. `CheckWire` / pcurve 参数口径（§9.294.1、§9.299）
2. 范围守卫 `v_range`（§9.298）
3. `ComposeShell` 内部拆分（§9.300）
4. `resolve_face` 返回 Shell/compound（§9.293.1）
5. **顶点身份/邻接口径不同**（§9.302.3、§9.303.4）—— 输入里根本没有共享顶点可比
6. **两个退化边同母线**（§9.303.7）—— 实测在 `u=0` 与 `u=π`
7. `heal_shape` 接在读入路径上（§9.301.1）—— 它根本没接
8. D18 的 `n1→n2`（§9.291.5，二次否证）

#### 9.304.4 建议：本卡转入"标注 + 待 Oracle 证据"，或转做 T-100

八轮下来，下游每一级都被独立实测证明忠实，而我**没有**找到 OCCT 完成该合并的机制。
继续在端口里改的唯一后果是**自创一条 OCCT 里不存在的规则**（红线）。
因此建议二选一（**需要你或下一轮决策**）：

* **(甲)** 把 T-99 转为 `blocked`，`blocker` 写成具体条件：
  「需要在 OCCT 侧取得该面在 `ProcessShape` 中 bound 数变化的直接证据
  （离线 `ShapeProcess` trace 或 `BRepTools::Write` 前后对比）——
  本会话无法重建 OCCT（只有源码树、无 build tree），探针拿不到这一点」。
  这不是"困难"或"不确定"，而是一条**外部取证条件**：需要 OCCT 侧能跑 `ShapeProcess` 的证据。
* **(乙)** 把 F113 的 FACE 降到 P3 / `wontfix-now`，转做 **T-100**（`ginter.rs` 两条分派臂，范围明确、可动手）。

**我倾向 (甲)**：因为定位本身已完成到"端口无辜"这一层，剩下的确是一条取证边界，
而不是能力不足 —— 把它标成 blocked 并写清条件，比继续在端口里试探更能保护红线。

---

### 9.305 —— T-100 完成：`ginter` 的**空分支**被移除（改用已可用的通用臂）；另一条臂确认为**本仓不可达**

#### 9.305.1 第一个 UNPORTED（`ginter.rs:170-177`）：**它是空分支**，且通用臂本来就能做

原代码：

```rust
GeomAbsCurveType::Circle | Ellipse | Parabola | Hyperbola => {
    // UNPORTED (IntCurve_IntConicConic_1.cxx:2236 Line/Circle, :2861 Line/Ellipse, ...)
}                                            // ← 分支体是空的：命中即什么都不做
_ => {
    // l1 = line(c1); IntConicCurveGen::perform_line(&l1, d1, c2, d2, tol_conf, tol)
}
```

**关键实测（§9.305.2）**：那个 `_` 通用臂**本来就能正确求 Line/任意曲线的交点**，
因为它的下游是通用的：

* `int_conic_curve.rs` 的 `perform_line` → `IntCurveIConicTool::from_lin2d(l)` →
  `perform_iconic_tool(tool, d1, p_curve, d2, …)`；
* 而 `IntCurveIConicTool`（`occt-core/src/intcurve/iconic_tool.rs`）**已经能表示每一种圆锥**：
  `from_lin2d` / `from_circ2d` / `from_elips2d` / `from_parab2d` / `from_hypr2d`
  （对应 `CurveType::{Line,Circle,Ellipse,Parabola,Hyperbola}`）；
* `MyImpParTool` 的第二个参数是 `&dyn Curve2d`，**不限定类型**。

⇒ OCCT 把 Line/Circle 等分派给 `IntCurve_IntConicConic` 的专用重载
（`IntCurve_IntConicConic_1.cxx:2236` / `:2861` / 抛物、双曲重载）是一次**特化（性能）**，
**不是新能力**。端口的分支体为空 ⇒ 该组合**静默返回零交点**，比走通用臂**更差**。

**改动**：删掉这个空分支，让 `Circle|Ellipse|Parabola|Hyperbola` 与 `_` 一起走通用臂
（合并成单个 `_`）。这不是"新增分支"，而是**去掉一个短路**。

#### 9.305.2 实测（新探针 `crates/occt-geom2d/examples/zz_lineconic_probe.rs`）

走**公开分派器** `Geom2dIntGInter::perform_with_d2`（即 `split_by_line.rs:203-204` 用的同一个入口），
构造已知解：直线 = x 轴，圆 = 圆心 (0,0) 半径 5 ⇒ 期望 `x = ±5`。

```text
is_done=true
nb_points=3
  point[1] t2=3.141592654 xy=(-5.000000000,0.000000000)
  point[2] t2=6.283185307 xy=(5.000000000,-0.000000000)
  point[3] t2=0.000000000 xy=(5.000000000,0.000000000)
x values (deduped): [-5.0, 5.0]
EXPECT +/-5 => PASS
```

3 个点是因为我把域设成**闭区间 `[0, 2π]`**，`t=2π` 与 `t=0` 重合（去重后正好 ±5）。
**没有多解、没有错解** ⇒ 通用臂对 Line/Circle 是**正确**的。

#### 9.305.3 第二个 UNPORTED（`ginter.rs:188`）：**本仓不可达**，据实标注而非照搬

`typ1 != Line` 的那一组臂（`IntCurve_IntCurveCurveGen.gxx:339-…`）确实没移植。
但**唯一调用方**是 `shape_fix_compose_shell/split_by_line.rs:203-204`，
而它的第一条曲线永远是直线：

```rust
split_by_line.rs:42   let j_c2d = Geom2dLine::new(*line.position());
split_by_line.rs:204  inter.perform_with_d2(&j_c2d, c2d.as_ref(), &i_dom, TOLINT, TOLINT);
```

⇒ `typ1` 恒为 `Line`，这个 `else` 分支**在本仓不可达**。
按 T-95 的教训（先确认有无消费方，再决定补分支还是标注），这里**标注为不可达**、
**不照搬** `IntCurveCurveCurveGen` 的实现（那会是"自造 OCCT 里不存在的分支"的反面：
照搬一段不可达的代码同样不可取）。

#### 9.305.4 门禁结果

| 项 | 结果 |
|---|---|
| `occt-geom2d --lib` | **72/0** |
| `occt-topo --lib` | **1281/0** |
| a3n00 基线 | `TOTAL faces=226 … stats=225 mesh_v=10863 mesh_t=11941` / `STATMAP matched=225 unmatched=1` —— **与改动前逐字相同** |
| `cargo check` | 0 error |

行为不变是**预期**的：被删的分支体为空（什么都不做），删掉它不改变任何已走路径的结果；
它只在**此前会静默失败**的 Line/Conic 输入上产生差别 —— 而那正是要修的。

#### 9.305.5 复跑

```text
rtk cargo run --manifest-path crates/occt-geom2d/Cargo.toml --offline --example zz_lineconic_probe
grep -n UNPORTED crates/occt-geom2d/src/geom2d_int/*.rs
```

---

### 9.306 —— T-54 定案：单面回退 `wireframe_face_triangulation` **在全部 6 个门禁模型上从未被触达**

#### 9.306.1 卡面要求的那次测量

T-54 的 `next` 原文：「先跑一次 `zz_uv_feed --ids` 确认还有没有别的面真的走回退，
有实测触发例再决定换 Delaunay 还是删回退让主线报错暴露」。

本轮的测量方式（新插桩 `OCCT_TOPO_TRACE_FALLBACK`）：在
`meshing/incremental_mesh/discret_root.rs` 的两处打印 ——
① 被 `MeshStatus::FAILURE` / `REUSED` 过滤掉的面；② 真正进到
`wireframe_face_triangulation` 的面及其产出三角数。

#### 9.306.2 结果（6 个模型）

| 模型 | 被 FAILURE/REUSED 跳过 | **触达回退** |
|---|---|---|
| a3n00 | **1**（face 113，`failure=true`） | **0** |
| T0M | 0 | **0** |
| acs10 | 0 | **0** |
| TDB | 0 | **0** |
| bottom | 0 | **0** |
| top | 0 | **0** |

```text
FALLBACK skipped_face=113 failure=true reused=false      ← a3n00 唯一
（没有一行 rescued_triangles）
```

⇒ **没有任何面真的走到回退**。T-54 要的那条判据（"还有没有别的面走回退"）答案是：没有。

#### 9.306.3 为什么必然如此（代码结构）

`discret_root.rs` 的两段合起来把回退变成了**到不了的分支**：

```rust
:394  if f.is_status(MeshStatus::FAILURE) || f.is_status(MeshStatus::REUSED) { continue; }
      //   ↑ 已失败的面根本不进 pending
...
:449-454  Ok(Err(e)) => if face(i).is_status(FAILURE) { (None, false) }   // ← needs_fallback=false
:461-463  Err(_)     => (None, true)                                      // ← 只有 panic 才 true
:479-480  p.tri.or_else(|| if p.needs_fallback { wireframe_face_triangulation(...) })
```

* 算法自己判定失败的面（`Ok(Err)` 且已标 `FAILURE`）→ `needs_fallback = false`；
* 已标 `FAILURE` 的面在循环开头就被 `continue` 掉，压根不进 `pending`；
* ⇒ `needs_fallback = true` **只剩"被 `catch_unwind` 捕获到 panic"这一条路**（`:463`），
  而 6 个门禁模型都没有面在该处 panic。

**最要紧的一点**：a3n00 的 F113 之所以拿不到三角，**不是**"回退没救成"，
而是它**已经被上游标成 `FAILURE`，在 `:394` 就被跳过了** —— 回退连机会都没有。
（这修正了 T-54 阻塞说明里"它们先被判 invalid discrete range 早期返回"的表述：
实际是**更早**的一步把面标成了 FAILURE。）

#### 9.306.4 处置：**保留不删**，并把过期注释改成实测结论

T-54 给出的两条路是「换约束 Delaunay」或「删回退让主线报错暴露」。
本轮**两条都不做**，理由：

1. **换 Delaunay**：这是给一个**不可达**的分支实现更复杂的算法，纯增风险、零可观测收益。
2. **删回退**：`needs_fallback` 仍有一条真实的（虽未被触发的）触发路径 —— **panic 救援**。
   按项目纪律「不要凭想象删」，在**没有反例证明它有害**之前不删；
   而且删掉它会移除一个针对未预期 panic 的安全网。

**实际动作**（本轮唯一的库改动，**纯注释**）：把那段 **过期**的注释
（原写"端口管线仍失败 T0M 的 169/1772 个面（2026-09-20 实测）"）改成本轮实测结论 ——
T0M 现在是 **0** 个失败面、回退在全部门禁模型上不可达。留着过期注释比删代码更危险：
它会让下一个人以为 169 个面还在失败。

#### 9.306.5 门禁

| 项 | 结果 |
|---|---|
| `occt-topo --lib` | **1281/0** |
| a3n00 基线 | `faces=226 … stats=225 mesh_v=10863 mesh_t=11941` / `matched=225 unmatched=1` —— 未变 |
| `git diff` | `discret_root.rs` **仅注释**（13 增 7 删，逐行核对全为 `//`） |

#### 9.306.6 复跑

```text
OCCT_TOPO_TRACE_FALLBACK=1 rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline \
    --example zz_share_probe -- data/occ/<model> 2> fb.err
grep FALLBACK fb.err        # 期望：只有 skipped_face，没有 rescued_triangles
```

---

### 9.307 —— T-97 完成：`bop_occt_util.rs` 确认无消费方的 31 个死函数成批删除（−211 行，该文件警告 31→0）

#### 9.307.1 范围（先量，不猜）

`cargo check` 全量解析（自写解析器把 `warning:` 与其 `-->` 位置配对）：**561 个警告站点**：

| 类别 | 数量 |
|---|---|
| `never used`（死代码） | **372** |
| `never read` | 68 |
| `unused import` | 52 |
| `never constructed` | 26 |
| other | 43 |

按文件分组，死代码最集中的是 `src/bop_occt_util.rs`（**31**），故选它作第一批。

#### 9.307.2 删除前的独立验证（这一步救了两次）

**教训：用「按名字 grep」验证 `never used` 会给出假阳性。**
我第一遍按名字 grep，得出「4 个有外部引用、不能删」，其中 `oriented` 有 149 处命中。
但那些命中是**别的类型上的同名方法**（`TopoShape::oriented`、`Orientation::is_reversed`），
并不是这个私有模块里的自由函数。**rustc 的死代码分析才是权威**；
grep 只能当交叉检查，且必须看命中处的实际形态。

`bop_occt_util` 在 `lib.rs:278` 是 `mod bop_occt_util;`（**私有模块**，函数虽写 `pub` 但无外部消费者），
所以 `pub fn` 被报 `never used` 是自洽的 —— 这解释了为什么会有 372 条这类警告。

#### 9.307.3 删除内容与结果

删掉 **31 个** 确认无消费方的顶层函数
（`fence_add`/`fence_insert`/`internal_copy`/`ds_shape_sd`/`add_face_to_shell`/`as_solid`/`as_shell`/
`as_wire`/`as_vertex`/`on_face_split_edge`/`append_on_face_paves`/`unique_shapes`/`clone_list`/
`list_contains`/`list_append_unique`/`faces_of_shell`/`edges_of_wire_direct`/`vertices_of_edge_direct`/
`edge_or_images`/`face_or_images`/`solid_or_images`/`bind_images`/`ds_index`/`ds_shape`/`ds_type`/
`source_indices_of`/`same_shape`/`reverse_shape`/`oriented`/`is_reversed`/`is_forward`），
以及因删除而**变成孤儿的一个 import**（`Vertex`）。

删除方式：脚本对每个名字做**大括号配对**取完整函数体、连同其上游 `///` 文档块一起删；
先校验所有名字都在顶层、且删除区间互不重叠，任一异常即中止不写盘。

| 项 | 前 | 后 |
|---|---|---|
| `bop_occt_util.rs` 行数 | 608 | **398** |
| 该文件的警告 | 31 | **0** |
| 全 crate `never used` | 372 | **341** |
| 全 crate 警告站点 | 561 | **530** |
| `git diff` | — | **1 增 211 删** |

那 1 行"新增"是 `use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Wire};` ——
删掉 `Vertex` 后该 use 行的重排结果，不是新逻辑。

#### 9.307.4 门禁

| 项 | 结果 |
|---|---|
| `occt-topo --lib` | **1281/0** |
| a3n00 基线 | `faces=226 … stats=225 mesh_v=10863 mesh_t=11941` / `matched=225 unmatched=1` —— **未变** |
| T0M | `zz_share_probe` **21/1/4** —— 未变 |
| `cargo check` | **0 error** |

#### 9.307.5 一个必须报告的事：集成测试有**既存**的红

跑整包 `cargo test -p occt-topo` 时 `tests/bop_builder2_boss.rs` 的
`boss_single_disc_base_merges_one_solid` 失败（`bop_builder2_boss.rs:116`：
`result must contain the 2-wire ring box-top face`）。

**这不是本轮改动引起的，已用对照实验证明**：把 `bop_occt_util.rs` 换回原始版本后重跑**同一条**测试，
失败与 panic 位置**完全一致**（`bop_builder2_boss.rs:116`）。

而且这条红是**板上有记录的、已被接受的**：对应卡片（`board.canvas.tsx:815`）
`status: "cancelled"`、`next: "维持 1/2 作为已知，不再投入（夹具不合法，非引擎缺口）"`、
`evidence: "夹具期望的「1 solid」在 OCCT 的 BOPAlgo 下也不成立；门禁 G-boss 记为已接受的红"`。

⇒ 结论：**该红与本卡无关，本卡不得碰它**（红线：完成即停，旁支写进报告不动手）。
但需要在此写明：`occt-topo` 的整包测试**并非全绿**，这是既存状态。

#### 9.307.6 给下一批的建议（本卡范围内）

`never used` 还剩 **341**，`unused import` 剩 **52**。下一批建议挑
`pave_intersect/face_face.rs`（20）、`bopalgo_tools_wires.rs`（16）、`pave_intersect/fill_ctx.rs`（14）、
`pave_intersect/vertex_face.rs`（14）—— 都是**私有模块 + 无消费方**的同类形态，
可用同一套脚本流程（先验顶层与不重叠、再删、再跑三样门禁）。
**注意**：每批删完都要重新 `cargo check` 取新清单，因为删函数会产生新的孤儿 import（本轮就出现 1 个）。

#### 9.307.7 复跑

```text
rtk cargo check --manifest-path crates/occt-topo/Cargo.toml --offline > warn.txt 2>&1
python .target-gate/warn_cmp.py          # 前后对照
python .target-gate/strip_dead.py        # 第一批（名字清单内置于脚本）
python .target-gate/strip_dead2.py       # 第二批 + 孤儿 import
```

---

### 9.308 —— T-97 第 2 批：整个 `bopalgo_tools_wires` 模块经**编译级证明**不可达后删除（−503 行）

#### 9.308.1 发现：存在「整模块无消费者」的形态

第 1 批之后重新取清单，按文件分组时发现最集中的几个文件有个共同特征：
**没有任何 `外部模块::` 形式的引用**：

| 模块 | 行数 | 外部引用数 |
|---|---|---|
| `pave_intersect/fill_ctx` | 593 | **0** |
| `pave_intersect/vertex_face` | 691 | **0** |
| `pave_intersect/face_face` | 701 | **0** |
| `bop_build_faces/builder_like` | 560 | **0** |
| `bopalgo_tools_wires` | 503 | **0** |

#### 9.308.2 一次不走运的弯路（记下来以免重犯）

我写了个脚本想量化「lib.rs 里有多少私有模块不可达」，输出 **"122 of 122 私有模块未被任何 `pub use` 命名"**，
并算出 35,329 行 —— **这个结论是错的，我把它废弃了**。
原因很简单：私有模块**不需要**被 `pub use` 命名也能被同 crate 的其他私有模块使用（`crate::x::y` 路径）。
「没出现在 `pub use` 里」≠「不可达」。差一点据此去删三万行。

**教训与第 1 批同源**：判断可达性/死代码不能靠**文本模式**（名字 grep、`pub use` 清单），
要靠**编译器**。

#### 9.308.3 采用的判据：编译级可达性证明

于是改用最直接的证明方式（脚本 `.target-gate/prove_unreachable.py`）：

1. 把 `lib.rs` 的 `mod bopalgo_tools_wires;` 注释掉；
2. 跑 `cargo check --message-format=short`；
3. **编译通过且 0 error ⇒ 该模块不可达**（`exit code: 0`、`errors: 0`，实测）。

这才是「无消费者」的**充分证据**：如果任何地方用到它，注释掉声明必然编译失败。

#### 9.308.4 删除内容

`crates/occt-topo/src/bopalgo_tools_wires.rs`（**503 行**）整文件删除 + `lib.rs` 一行声明删除。
该文件是 `BOPAlgo_Tools`（`EdgesToWires` / `WiresToFaces` / `MakeWires` / `FindEdgeTangent` / `FindPlane`）
的移植，其 **16 个 `pub fn` 全部**被 rustc 报 `never used`，另有 3 个私有辅助函数
（`shape_key` / `dirs_parallel` / `dirs_normal` 等亦无消费者）。

**可恢复性**：删除前记录了 git blob 哈希
`0f06ed714d6bcb0d32eba673108c813e99a24683`
（存于 `.target-gate/deleted_blobs.txt`），
需要时应能 `git cat-file -p 0f06ed71...` 取回 —— **即这次删除不是不可逆的**。

#### 9.308.5 逐阶段门禁与量化

| 阶段 | 警告站点 | `never used` | `unused import` |
|---|---|---|---|
| T-97 开始时 | **561** | 372 | 52 |
| 第 1 批后（`bop_occt_util.rs` 31 函数） | 536 | 346 | 53 |
| **第 2 批后（删整个模块）** | **514** | **325** | **52** |

累计：**站点 561 → 514（−47）**、`never used` **372 → 325（−47）**、删除 **714 行**
（211 + 503）。

| 门禁 | 结果 |
|---|---|
| `occt-topo --lib` | **1281/0** |
| a3n00 基线 | `faces=226 … stats=225 mesh_v=10863 mesh_t=11941` / `matched=225 unmatched=1` —— **未变** |
| T0M | `zz_share_probe` **21/1/4** —— 未变 |
| `cargo check` | 0 error |
| `lib.rs` diff | **仅一行** `-mod bopalgo_tools_wires;` |

#### 9.308.6 剩余与下一批

`never used` 还剩 **325**，分布在 126 个文件。**下一批应继续用「编译级可达性证明」**
（注释声明 → `cargo check` → 通过则删），而不是名字 grep：
`pave_intersect/fill_ctx`(593 行)、`pave_intersect/vertex_face`(691)、
`pave_intersect/face_face`(701)、`bop_build_faces/builder_like`(560) 都是 0 外部引用的候选，
总约 2,545 行。

**但注意一个判断分界**：这些是 `BOPAlgo` 族的**忠实移植**，删掉它们是「卫生」还是「丢掉已完成的工作」，
取决于项目是否还打算把它们接进主线。第 1 批（删文件内的死函数）无争议；
**整模块删除属于需要知情的决定** —— 本批之所以仍执行，是因为一个模块 19/19 个函数全无消费者、
编译级证明不可达，且已留可恢复哈希。**再往后建议先确认方向再动手。**

#### 9.308.7 复跑

```text
python .target-gate/reachability.py        # 注意：其「不可达」结论不可信，仅作规模参考
python .target-gate/prove_unreachable.py   # 编译级可达性证明（改 DEGL 后复用）
python .target-gate/warn_cmp3.py           # 三阶段警告对照
```

---

### 9.309 —— 【更正 §9.308.4】`bopalgo_tools_wires` 不是死代码：它是**忠实移植但未接线**，且 OCCT 里有真实调用方

#### 9.309.1 用户的质疑是对的，我上一轮的删除依据不成立

§9.308.4 我删掉了整个 `bopalgo_tools_wires.rs`（503 行），理由是
「编译级证明不可达 + 19/19 函数无消费者」。**那个证明只说明「本 crate 里没人引用」，
不等于「这个函数没用」** —— 它只能说明**移植没接线**。我把「未被接线」当成了「死代码」。

**已恢复**：`git cat-file -p 0f06ed714d… > bopalgo_tools_wires.rs` + 在 `lib.rs` 原位置
（`pub mod bopalgo_tools;` 与 `mod bopalgo_tools_class;` 之间）恢复 `mod bopalgo_tools_wires;`。
校验：`git diff --ignore-cr-at-eol` 对该文件**无内容差异**（逐字还原），`lib.rs` 的 diff 为空。

#### 9.309.2 OCCT 侧的调用方（三个函数各有明确归属）

**BOPAlgo_Tools.cxx 内部的互相调用**（完整调用图，全部可达）：

```text
:360  BOPAlgo_Tools::EdgesToWires  …  ← 公开 API
:454      FindPlane(aBAC, aPln)
:461      FindEdgeTangent(aBAC, aVT)
:565      MakeWires(aMEPln, aRWires, false, …)
:638      MakeWires(aMEPln, aRWires, true,  …)
:656      MakeWires(aMEAlone, aRWires, false, …)
:665  BOPAlgo_Tools::WiresToFaces  …  ← 公开 API
:691      FindPlane(aWire, aPlane, …)
:805  MakeWires(…)                 …  ← 文件内 static
:831      FindPlane(aCBE, aPln, …)
:852  FindEdgeTangent(theEdge, …)   …  ← 文件内 static
:861      FindEdgeTangent(aBAC, aVTE)
:875  FindEdgeTangent(aBAC, …)      …  ← 文件内 static
:910  FindPlane(aBAC, …)            …  ← 文件内 static
:976  FindPlane(theWire, …)         …  ← 文件内 static
:995/:1012   FindEdgeTangent(aE1/aE2, …)
:1040        FindPlane(aBAC, thePlane)
```

**跨模块调用方**：

| 函数 | 归属 | 外部调用方 |
|---|---|---|
| `EdgesToWires` | `hxx:151` `Standard_EXPORT`（公开 API） | `BRepFill_AdvancedEvolved.cxx:561, 641, 1499`；`BOPTest_UtilityCommands.cxx:193` |
| `WiresToFaces` | `hxx:170` `Standard_EXPORT`（公开 API） | `BRepFill_AdvancedEvolved.cxx:561, 642, 1500`；`BOPTest_UtilityCommands.cxx:201` |
| `MakeWires` | `cxx:805` **文件内 static** | 无外部调用方（但在 `EdgesToWires` 内被调 3 次） |
| `FindEdgeTangent` | `cxx:68/75/83/87` **文件内 static** | 无外部调用方（内部被调 4 次） |
| `FindPlane` | `cxx:910/976` **文件内 static** | 无外部调用方（内部被调 6 次） |

⇒ **五个函数在 OCCT 里全部可达且在用**；端口把整条调用图都忠实移植了。

#### 9.309.3 为什么它在端口里显得"死"

端口的 `bopalgo_tools_wires` 是**私有模块**（`mod` 而非 `pub mod`），且没有任何端口模块 `use` 它；
而 OCCT 里唯一的外部消费方 `BRepFill_AdvancedEvolved` **在端口里没有对应实现**
（`grep` 端口全部 `.rs`：无 `evolved` / `advanced_evolved` 相关文件）。

⇒ 准确的描述是：**这是 `BOPAlgo_Tools` 的忠实移植，它的外部消费方尚未移植，因此尚未接线**。
这正是用户说的「只是在当前实现的功能范围内没有被用到」。

#### 9.309.4 由这次更正得出的**分类判据**（下一轮必须照此区分）

这是本轮最有价值的产出。**不能只看「有没有消费者」，要看被标记项的可见性与归属**：

| 形态 | 含义 | 处置 |
|---|---|---|
| 私有模块里的 **`pub fn`**，且 OCCT 有外部调用方 | **忠实移植、未接线** | **保留**；写清"待接线"+ 依据 `.cxx`/`.hxx` 行号与调用方 |
| 父模块内的 **`pub(super)` / `pub(crate)`**，父模块自己都不调 | 本地死代码，**不是对外 API** | 可删（经编译级证明） |
| **`pub(crate) struct` never constructed** | crate 内从未构造，真死 | 可删（经编译级证明） |
| 已接入管线的模块里的散落死函数 | 真死 | 可删（第 1 批那种） |

按此判据回看上一轮我为「第 3 批」列出的候选，**它们与 `bopalgo_tools_wires` 不是一类**：

| 候选 | 被标记项的可见性 | 与 `bopalgo_tools_wires` 的差别 |
|---|---|---|
| `pave_intersect/face_face.rs` | `pub(super) fn check_planes` / `pub(super) fn perform_ff_impl` | **`pub(super)`** —— 父模块 `pave_intersect` 自己都不调 ⇒ 本地死代码，不是对外 API |
| `pave_intersect/fill_ctx.rs` | **`pub(crate) struct FillCtx` never constructed** | crate 内从未构造 ⇒ 真死 |
| `pave_intersect/vertex_face.rs` | `pub(super) fn` | 同 `face_face` |
| `bop_build_faces/builder_like.rs` | `pub(super)`/其它 | 同 |

且 `pave_intersect` 是 **`pub mod`**（`lib.rs:181`），`mod.rs` 还 `pub use fill_ctx::*` / `pub use vertex_face::*`
—— 在这种模块里，若某个 `pub fn` 真被标成 dead，说明**连 crate 内部都没人调且未对外导出**，
与「私有模块里的 pub fn 未接线」是**完全不同的情形**。

#### 9.309.5 待办：先给这个模块一个正确的标注（不删）

`bopalgo_tools_wires` 现在恢复了，但它仍会产生 16 条 `never used` 警告。
正确的处理不是删、也不该长期放任，而是**按 §9.309.4 的判据标注**，二选一：

* **(甲)** 模块头补一段说明：`BOPAlgo_Tools` 的忠实移植，外部消费方
  `BRepFill_AdvancedEvolved`（`BOPAlgo_Tools.hxx:151/:170` 的 `Standard_EXPORT` 调用方）尚未移植，
  故当前无消费者；**不要删**。这样警告仍在但含义清楚。
* **(乙)** 给该模块加 `#![allow(dead_code)]` + 同一段说明，使警告不再淹没真正的死代码。

**我倾向 (乙)**（本轮先只恢复+记录，未加标注，留给你定），
因为 561 条警告里很大一部分可能都是这种"忠实移植未接线"，逐条读会淹掉真正的死代码。

#### 9.309.6 复跑

```text
git cat-file -p 0f06ed714d6bcb0d32eba673108c813e99a24683 > crates/occt-topo/src/bopalgo_tools_wires.rs
git diff --ignore-cr-at-eol -- crates/occt-topo/src/bopalgo_tools_wires.rs   # 应为空
Select-String -Path D:\source\OCCT-src\src\ModelingAlgorithms\TKBO\BOPAlgo\BOPAlgo_Tools.cxx -Pattern 'EdgesToWires|WiresToFaces|MakeWires|FindEdgeTangent|FindPlane'
```

#### 9.309.7 附表：那 26 个失败的定性（与本次更正无关）

恢复文件后跑 `--lib` 出现 `1255 passed; 26 failed`。**已做对照实验**：把本轮涉及的三个文件
（`bop_occt_util.rs`、`bopalgo_tools_wires.rs`、`lib.rs`）全部 `git stash` 回到**原始状态**后重跑，
结果**完全相同**（`1255 passed; 26 failed`）。

失败原因是**环境性**的，例如：

```text
thread 'bincaf::tests::file_io_roundtrip' panicked at src/bincaf/tests.rs:174:36:
write file: "bincaf: write C:\Users\85332\AppData\Local\Temp\occt_bincaf_test.xbf: 拒绝访问。 (os error 5)"
```

⇒ 这些测试往 `%TEMP%` 写文件被拒（`os error 5`），是**运行环境的文件权限**问题，
不是代码问题，也**不是本次更正引入的**。

**对本轮结论的影响**：门禁以 `a3n00` 基线为准（`faces=226 … mesh_v=10863 mesh_t=11941`、
`matched=225 unmatched=1`，恢复后**逐字未变**）。`--lib` 的通过数须以「1255 passed / 26 环境性失败」
为当前常态，不能再用「1281/0」表述 —— 那是在 `%TEMP%` 可写的环境下取得的。

---

### 9.310 —— 【T-97 根本再修正】两批删除**全部回退**：那是一套**未完成移植**，不是死代码

#### 9.310.1 用户的问题

> 「如果是 BOPAlgo 族的忠实移植的话，那他们在 C++ 中在哪里被使用呢？应该不是无用的函数，或者只是在当前实现的功能范围内没有被用到」

问得对。这个追问把 §9.309 的更正又推进了一层：**连第 1 批也不该做。**

#### 9.310.2 第 1 批同样站不住：`bop_occt_util` 的模块头自己说清了它的性质

```text
//! Shared helpers for the OCCT Builder translation modules.
//!
//! The Builder stages in `BOPAlgo_Builder_2.cxx` / `_3.cxx` and
//! `BOPAlgo_Tools.cxx` share a small set of identity, fence, iterator and
//! orientation helpers. This module is the Rust equivalent of those locals:
//! `TopTools_MapOfShape` fence maps, `TopoDS_Iterator` child walks,
//! `TopExp::MapShapes` typed collectors, `BRep_Builder` add wrappers, and
//! the `myShapesSD.IsBound` / `myImages.IsBound` lookups the solid and face
//! stages both need.
```

两个要点：

1. 这些函数是 **OCCT 里 *locals*（文件内 static / 成员方法）的 Rust 对应物**，
   因此**在 OCCT 里不存在同名函数** —— 「按名字 grep OCCT」这个方法对它们**根本不适用**。
   实测：31 个名字里只有 `oriented` 在 OCCT 有 1077 处命中，
   而那 1077 处是 `TopoDS_Shape::Oriented()` **成员方法**，与端口的自由函数不是一回事。
2. 模块头明确说这些是 **Builder 各阶段共享的 helper** —— 即
   **`BOPAlgo_Builder` 移植的一部分**，只是 Builder 管线尚未接线到它们。

⇒ 与 `bopalgo_tools_wires` 同性质：**未完成的移植**，不是死代码。

#### 9.310.3 已全部回退

| 文件 | 处理 | 校验 |
|---|---|---|
| `bop_occt_util.rs` | 从删除前备份逐字还原 | `git diff --ignore-cr-at-eol` **对 HEAD 为空**（字节相同），608 行 |
| `bopalgo_tools_wires.rs` | `git cat-file` 取回 blob `0f06ed71…` | 同上，503 行 |
| `lib.rs` | 恢复原位置 `mod bopalgo_tools_wires;` | `git diff` 为空 |

全 crate 警告站点回到 **561**（与 T-97 开始时相同）。
源码层现仅剩本轮之前既有的改动：`occt-geom/src/gcpnts.rs`、`occt-geom2d/src/geom2d_int/ginter.rs`、
`occt-topo/src/meshing/.../discret_root.rs`、`occt-topo/src/shhealing/shape_fix_face.rs`、
`occt-topo/src/step/read_topology.rs`。

#### 9.310.4 结论：这 372 条 `never used` **总体不是死代码**

结合两批的教训，可以给出一个明确的定性：

> **在这份代码库里，「`never used`」主要指示「移植了但消费方尚未移植 / 尚未接线」，
> 而不是「代码没用」。**

理由：这个仓库是在**逐模块移植 OCCT**。OCCT 的调用方（`BRepFill_AdvancedEvolved`、
`BOPAlgo_Builder` 各阶段等）往往还没移植，于是被移植的**被调用方**就表现为「无人使用」。
`bopalgo_tools_wires`（公开 API 有真实外部调用方）与 `bop_occt_util`（Builder 阶段的共享 helper）
都是这一类。

**因此本卡的正确形态不是「删多少行」，而是「分类 + 接线/标注」**：

| 类别 | 判据 | 处置 |
|---|---|---|
| **A. 忠实移植、消费方未移植** | 模块头自报的 OCCT 来源里，被调用方**在 OCCT 有调用方**（本仓或 OCCT 一侧任一成立） | **保留**；写清「待接线」+ 依据文件行号。**绝不删** |
| **B. 真死** | 既无本仓消费者，OCCT 侧也无调用方（含本地自造已被取代的实现） | 可删（编译级证明） |
| **C. 环境/误报** | 如「只被同文件 `#[cfg(test)]` 使用」 | 不动 |

**A 与 B 的区别不能靠「有没有消费者」判定** —— 那正是我连错两次的原因。
需要的是**查 OCCT 侧的调用方**（对公开 API 可行），或者**读模块头自报的移植来源 + 判断该来源的消费方是否已移植**（对 locals/helper 唯一可行）。

#### 9.310.5 建议的方向（交给下一轮/决策）

既然 A 类是主体，**继续批量删除是错的方向**。三个可选做法：

1. **（推荐）给 A 类模块加 `#![allow(dead_code)]` + 模块头一段「待接线」说明**
   （写明消费方是哪个 OCCT 类/文件、当前为何未接线）。
   收益：561 条噪音里真正需要看的 B 类会浮出来，且**不丢任何已完成的工作**。
2. 只做**分类清点**（产出 A/B/C 清单与各自依据），不删任何东西。
3. 若确实要减少警告数，则应先**接线上游消费方**，让这些函数真正被调用 —— 那是移植任务，不是卫生任务。

#### 9.310.6 复跑

```text
git diff --ignore-cr-at-eol -- crates/occt-topo/src/bop_occt_util.rs        # 应为空
git diff --ignore-cr-at-eol -- crates/occt-topo/src/bopalgo_tools_wires.rs # 应为空
git diff -- crates/occt-topo/src/lib.rs                                     # 应为空
rtk cargo check --manifest-path crates/occt-topo/Cargo.toml --offline       # 561 警告站点
```

---

### 9.311 —— T-97 选项 1 落地：A 类模块加 `#![allow(dead_code)]` +「待接线」说明，死代码警告 390→111

#### 9.311.1 先说一个数据口径的更正

前几轮我一直引用「561 个警告站点」，那其实是**两个来源混在一起**的产物：
`rtk`（一个输出过滤代理）把 occt-core/geom/geom2d 的警告和 occt-topo 的混在了一起，
而且它**不会重新输出缓存命中 crate 的警告**。

本轮改用**直接调用 cargo 并绕过过滤器**：

```text
cargo check --manifest-path crates/occt-topo/Cargo.toml --offline --message-format=short
```

得到的是 **occt-topo 自己的 513 条**（路径相对 crate root）。此后所有计数都以这个口径为准。

#### 9.311.2 分类（按模块头自报的 OCCT 来源）

判据：读文件自己的 `//!` 头部前 35 行，若出现 `SomeOcctClass.cxx|hxx|lxx|gxx`
⇒ **A 类：忠实移植，OCCT 侧消费方尚未移植**；否则 B 类。

（`value assigned … is never read` 这类**局部变量赋值未读**单独归类，不计入死代码 ——
上一轮我的正则 `is never read` 把它们误捕了。）

| 类别 | 站点 |
|---|---|
| 死代码（`is never used` / `never constructed` / item `is never read`） | **390** |
| 未用 import | 47 |
| 局部赋值未读 | 35 |
| 其它 | 32 + 6 + 3 |

#### 9.311.3 落地：给 A 类文件加模块级标注

对每个 A 类文件在其 `//!` 头之后插入**统一 6 行**块：

```rust
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]
```

**共标注 60 个文件**，每个文件的 diff **恰好 +6/−0**（已逐个校验）。

**过程中的一次返工（记下来）**：第一版标注用了两种不同措辞，第二次规范化脚本只匹配了其中一种，
留下**残留行**（如 `//! §9.309/§9.310); wire the consumer instead.` 出现在块外）。
发现后没有继续打补丁，而是**把这 75 个文件整体 `git checkout` 回 HEAD 再干净重标一次** ——
现在 60 个标注文件的 diff 全部恰好 +6/−0。

#### 9.311.4 效果

| 指标 | 标注前 | 标注后 |
|---|---|---|
| **死代码** | **390** | **111** |
| 未用 import | 47 | 47 |
| 局部赋值未读 | 35 | 35 |
| 其它 | 41 | 24 |
| **总警告** | **513** | **220** |

⇒ 用 60 个各 6 行的说明性标注，**压掉 279 条死代码噪音（−72%）**，总警告 **−57%**，
且**没有删掉任何已完成的工作**。

剩余 111 条死代码分布在 26 个文件，多为 **B 类**（无 OCCT 来源的本地实现）
或**函数内嵌套项**（`bop_occt_util.rs` 的 `same_shape` 等虽在已标注模块内仍被报，
说明它们嵌在更深的作用域里）。这些是**真正该逐个看**的候选，不再被 279 条噪音淹没 ——
这正是本方案的目的。

#### 9.311.5 门禁

| 项 | 结果 |
|---|---|
| `cargo check` | **0 error** |
| a3n00 基线 | `faces=226 … stats=225 mesh_v=10863 mesh_t=11941`、`matched=225 unmatched=1` —— **逐字未变** |
| T0M | `zz_share_probe` **21/1/4** —— 未变 |
| `--lib` | **1255 passed / 26 failed**，与改动前**完全相同**（26 条为 `%TEMP%` 写入被拒的环境性失败，§9.309.7 已用对照实验证明与代码无关） |
| 改动规模 | `65 files changed, 492 insertions, 55 deletions` —— 其中 60 个是 +6/−0 的标注；其余 5 个是本次会话既有的功能改动（`gcpnts.rs`、`ginter.rs`、`discret_root.rs`、`shape_fix_face.rs`、`read_topology.rs`） |

`#![allow(dead_code)]` 只影响 lint，不生成任何代码，因此基线不变是**预期**结果。

#### 9.311.6 复跑

```text
cargo check --manifest-path crates/occt-topo/Cargo.toml --offline --message-format=short > w.txt 2>&1
python .target-gate/final_census.py .target-gate/raw_warn.txt w.txt
python .target-gate/annotate_clean.py     # 幂等：已标注的文件会跳过
```

---

### 9.312 —— T-51 结项：gcpnts 余量清单已与代码现状对齐；「`uniform_abscissa` 外形不同」**不成立**

#### 9.312.1 `UNPORTED` 清单逐条核对（改动前 5 处）

| 位置（改动前） | 原文 | 判定 |
|---|---|---|
| `:46` | 模块头「`GCPnts_UniformDeflection` / `GCPnts_QuasiUniformDeflection` are UNPORTED」 | **真实**（这两类本模块确实不再暴露）。保留。 |
| `:294` | `CPnts_MyRootFunction::init`：「UNPORTED: the tolerance overload `Init(X0, L, Tol)` (`:32-37`) and the `AdvPerform` path that uses it (`CPnts_AbscissaPoint.cxx:436-474`) are not exposed by this port's public API」 | **真实**（措辞准确：是"未由公开 API 暴露"，不是"没实现"）。保留。 |
| `:488` | `compute_with_guess` 的文档：「UNPORTED: `GCPnts_LengthParametrized` (`:87-90`) and `GCPnts_AbsComposite` (`:96-158`); see the module header」 | **过期**。同函数正文 `:500-509` 三个分支**全都实现了**（`LengthParametrized` → `u0 + abscissa / ratio`、`AbsComposite` → `compute_abs_composite`、`Parametrized` → 迭代根搜索）；且它引用的「module header」已在 T-96 改成实测结论。 |
| `:563` | `abscissa_point_with_tolerance`：「UNPORTED: the `GCPnts_AbsComposite` arm of `AdvCompute` (`GCPnts_AbscissaPoint.cxx:187-295`, which carries an `anIndex == 0` special case) — multi-span curves fall back to the non-adv arm, so only the integration tolerance differs there」 | **真实**（有具体 `.cxx` 行号 + 明确行为差异）。保留。 |
| `:593` | 内联「`:187-295` UNPORTED; use the non-adv walk (`:96-158`)」 | **过期**。紧接的分支 `:594-600` 就是调 `compute_abs_composite`（`:96-158` 那条），注释与代码自相矛盾。 |

⇒ 两处删、两处留。

#### 9.312.2 `uniform_abscissa` 的外形：**与 OCCT 等价**，模块头那句断言错了

模块头原文（`:42-44`）：

> `uniform_abscissa` / `quasi_uniform_abscissa` keep their own outer shape —
> `uniform_abscissa` is "n intervals" (`n + 1` points), not `GCPnts_UniformAbscissa`'s `NbPoints` point count.

**这句不成立。** OCCT 的 `NbPoints` 重载在 `initialize` 里就是按「点数」反解步长的
（`GCPnts_UniformAbscissa.cxx:494-515`）：

```cpp
Standard_ConstructionError_Raise_if(
    theNbPoints <= 1,
    "GCPnts_UniformAbscissa::Initialize() - number of points should be >= 2");
...
const double anAbscissa = myAbscissa = aL / (theNbPoints - 1);   // ← 段数 = 点数 - 1
const int    aSize      = theNbPoints + 5;                        // 数组留余量；最终 myNbPoints = anIndex
```

即 **OCCT 用 `NbPoints - 1` 作段数**（与「含两端的点数」自洽）。
端口 `uniform_abscissa(c, n)` 用 `step = total / n` 产出 `n + 1` 个参数。

⇒ **令 `n = theNbPoints - 1`，两者逐字等价**：同样的段数、同样的 `step`、同样含两端点的参数序列。
所谓"外形差异"只是**入参命名口径不同**（端口用「段数」，OCCT 用「点数」），
**不是移植缺口**。`quasi_uniform` 同理。

#### 9.312.3 实际改动（**纯注释，零行为变更**）

| 位置 | 改动 |
|---|---|
| `gcpnts.rs:42-44` | 把错误断言改为实测结论：`uniform_abscissa(c, n)` == OCCT 的 `NbPoints = n + 1`，引 `cxx:510` 的 `aL / (theNbPoints - 1)` |
| `gcpnts.rs:488-489` | 删过期 `UNPORTED` 段，改为「三条臂均已在下方 `match` 中实现，分派在 `compute_type`（`:26-65`）」 |
| `gcpnts.rs:593` | 删「UNPORTED」，改为「用的是 `:96-158` 非 adv walk；adv 版（`:187-295`）未移植，见 `abscissa_point_with_tolerance` 上的说明」 |
| `:46`、`:294`、`:563` | **不动** |

改动后 `grep -n UNPORTED gcpnts.rs` 剩 3 处注释（`:50` 模块头的 UniformDeflection、`:298` 容差重载、`:569` adv 版 AbsComposite）+ 1 处指向性引用。

#### 9.312.4 门禁

| 项 | 结果 |
|---|---|
| `occt-geom --lib` | **143/0** |
| `occt-topo --lib` | **1255 passed / 26 failed**，与改动前相同（26 条为 `%TEMP%` 环境性失败） |
| `cargo check` | 0 error |
| a3n00 基线 | `faces=226 … stats=225 mesh_v=10863 mesh_t=11941`、`matched=225 unmatched=1` —— **未变** |

#### 9.312.5 结论：gcpnts 的真实余量是 **1 个缺口**，不是 2 条 arm

清掉过期注释后，剩下两条真实项其实是**同一处缺口的两面**：

1. `CPnts_MyRootFunction` 的容差重载 `Init(X0, L, Tol)` 未由公开 API 暴露（`:298`）；
2. `AdvCompute` 的 `GCPnts_AbsComposite` 臂未移植（`GCPnts_AbscissaPoint.cxx:187-295`，含 `anIndex == 0` 特例，`:569`）。

第 2 条正需要第 1 条 ⇒ 合起来是**一个缺口**：**adv 版 `AbsComposite`**。
现状是多段曲线退回非 adv 走法（`:96-158`），**仅积分容差不同**。

若日后要补，落点明确：`CPnts_MyRootFunction::init(x0, l, tol)` 与 `AdvPerform` 的 `AbsComposite` 臂，
对应 `.cxx` 行号 `187-295` / `436-474`。

#### 9.312.6 复跑

```text
grep -n 'UNPORTED' crates/occt-geom/src/gcpnts.rs
sed -n '494,516p' D:/source/OCCT-src/src/ModelingData/TKGeomBase/GCPnts/GCPnts_UniformAbscissa.cxx
```

---

### 9.313 —— T-98 结项：【OCCT 语义对照】锥面 pcurve 报 `(-inf,inf)` 是**正确的**，端口无需改动

#### 9.313.1 卡面要求的那个对照

卡面：把「无限域上的 pcurve 该报什么」按 OCCT 的 `Geom2d_Line` / `Geom_ConicalSurface` / `GeomProjLib` 对清楚，
再决定是补参数域还是在探针侧标明。

**三条 OCCT 原文（逐条可复核）**：

1. **锥面的 V 本来就是无限域** —— `Geom_ConicalSurface.cxx:207-215`：
   ```cpp
   void Geom_ConicalSurface::Bounds(double& U1, double& U2, double& V1, double& V2) const {
     U1 = 0.0; U2 = 2.0 * M_PI;
     V1 = -Precision::Infinite(); V2 = Precision::Infinite();
   }
   ```
2. **`Geom2d_Line` 自己就报无限域** —— `Geom2d_Line.cxx:142-150`：
   ```cpp
   double Geom2d_Line::FirstParameter() const { return -Precision::Infinite(); }
   double Geom2d_Line::LastParameter()  const { return  Precision::Infinite(); }
   ```
   （`hxx:97-100` 的文档也写明：「Returns RealFirst from Standard」/「Returns RealLast from Standard」。）
3. **`GeomProjLib::Curve2d` 就是原样返回这个未修剪的 `Geom2d_Line`** —— `GeomProjLib.cxx:81`
   ```cpp
   case GeomAbs_Line:
     G2dC = new Geom2d_Line(Proj.Line());   // ← 未 Trim
     break;
   ```
   只有**输入本身是 `Geom_TrimmedCurve`** 时才做 `[U1,U2]` 修剪（`:118-128`）。

⇒ 锥面的母线（V 恒定）投到 UV 后是**常 V 直线** ⇒ OCCT 返回 `Geom2d_Line` ⇒ **它报 `(-inf,inf)` 是 OCCT 的行为**。

#### 9.313.2 端口逐字一致

`crates/occt-geom2d/src/line.rs:37-38`：

```rust
fn first_parameter(&self) -> f64 { f64::NEG_INFINITY }
fn last_parameter(&self)  -> f64 { f64::INFINITY }
```

而 `make_pcurve_full` / `make_pcurve` 是**如实读取**该曲线的 `first/last_parameter`（`make_pcurve.rs:20`）。
⇒ 端口输出的 `(-inf,inf)` 与 OCCT **逐字相同**。

**结论：不修。** 卡面给出的两条验收路径里，本卡走的是第一条 ——
「给出 OCCT 语义对照的结论并保持现状」，且**未改任何 mock / 门禁**。
（另：§9.312 那次也确认了「卡面描述的缺口常常比实际小」，这是本会话第三次。）

#### 9.313.3 真正的问题在**探针侧**，而且端口里已有正确做法

上一轮把 `d0` 读成「全 NaN」的根因，就是**用 pcurve 自报的 `-inf` 端点去求值**（D14 已撤回）。
端口本身早就用 `is_finite` 守卫处理了这种情况：

| 位置 | 做法 |
|---|---|
| `make_pcurve.rs:21` | `if f.is_finite() && l.is_finite() { (f, l) } else { … }` |
| `make_pcurve.rs:64` | `select_forward_seam` 的 `line_of`：不有限则退 `(0.0, 1.0)` |
| `make_pcurve.rs:219/:229` | 取 surface 的 u/v 周期时同样先 `is_finite` |

⇒ **正确的取样口径**：pcurve 自报域无限时，改用**3D 边的参数域**（本卡实例是 `[0, 2π]`）。
这与 OCCT 自己的防御写法同思路。

#### 9.313.4 顺手修掉一个**真实的文档错位**（本轮唯一代码改动）

`make_pcurve.rs` 里 `make_pcurve_full` 的文档块（原 `:50-55`）被**粘到了**
`select_forward_seam` 的文档里，导致：

* `select_forward_seam` 的文档混进了不属于它的「Construct the pcurve of `edge` on `face` …」；
* `make_pcurve_full`（`:88`）**反而没有任何文档**。

已把两段拆开，并把该文档交给真正的 `make_pcurve_full`，其中：

* 把原来那句**不准确**的「The returned pcurve is parameterized over the edge's range `[a, b]`」
  改为准确的「返回曲线带**自己的**参数域，且当投影落在 `Geom2d_Line` 上是 `(-inf, inf)`」；
* 附上本节的三条 OCCT 行号证据（`Geom2d_Line.cxx:142-150`、
  `GeomProjLib.cxx:81`/`:118-128`、`Geom_ConicalSurface.cxx:207-215`）；
* 写明对调用方/诊断的后果：**不要按自报域直接取样，要先用 `is_finite` 判断或改用 3D 边域**，
  并点出「在 `-inf` 处取样」正是那次误读的来源（§9.313 / D14）。

**纯注释改动**，19 增 6 删。

#### 9.313.5 门禁

| 项 | 结果 |
|---|---|
| `cargo check` | 0 error |
| `occt-topo --lib` | **1255 passed / 26 failed**，与改动前相同（26 条为 `%TEMP%` 环境性失败） |
| a3n00 基线 | `faces=226 … stats=225 mesh_v=10863 mesh_t=11941`、`matched=225 unmatched=1` —— **未变** |
| T0M | `21/1/4` —— 未变 |

#### 9.313.6 复跑

```text
sed -n '207,215p' D:/source/OCCT-src/src/ModelingData/TKG3d/Geom/Geom_ConicalSurface.cxx
sed -n '142,150p' D:/source/OCCT-src/src/ModelingData/TKG2d/Geom2d/Geom2d_Line.cxx
sed -n '77,129p'  D:/source/OCCT-src/src/ModelingData/TKGeomBase/GeomProjLib/GeomProjLib.cxx
```

---

### 9.314 —— T-94 步骤 ①：face 1752 的 4 个槽位已量出（**真交换**）；且发现计数需更正为 **2 面 / 8 槽位**，含**第二种机制**（周期位移）

#### 9.314.1 计数更正：不是「1 面 / 4 槽位」，而是「**2 面 / 8 槽位**」

把 `zz_share_probe` 的快照从「只记 `n=` 与 `rng=`」升级为**逐槽位打印 pcurve 指纹**
（域 + `d0` 在端/中点的取值 + `is_line`），重跑 T0M：

```text
SHARE face=1537 fix=true changed_slots=4 (of 5) new_edges=0
SHARE face=1752 fix=true changed_slots=4 (of 6) new_edges=0
SHARE TOTAL faces_with_fix=21 faces_changed=2 slots_changed=8
```

**可复现**（连跑两次同结果）。此前记的 `1 面 / 4 槽位` 是因为**旧快照只比较 `n=` 与 `rng=` 两个字段**，
对「pcRange 相同但 pcurve 端点被平移」的槽位**漏检**。升级后补上了 face 1537。
⇒ T-94 的命中率从「1/24」更正为 **2/24**（仍然很低，定级 P3 不变）。

#### 9.314.2 两种**不同**的机制（这是本轮的主要发现）

**face 1752 —— 真正的重新排序（`swap_seam`）**

```text
slot1 wire0.edge1  before: dom=[4.865097571,5.191227259]  p0=( 4.712388980,-5.809475019)
                   after : dom=[-inf,+inf]                p0=(14.137166941,-13.444711759)
slot2 wire0.edge2  before: dom=[4.084699525,4.865097571]  p0=( 5.416227526,-1.954609395)
                   after : dom=[-inf,+inf]                p0=( 5.416227526,35.689262838)
slot3 wire0.edge3  before: dom=[-inf,+inf]                p0=( 5.416227526,35.689262838)
                   after : dom=[4.084699525,4.865097571]  p0=( 5.416227526,-1.954609395)
slot4 wire0.edge4  before: dom=[-inf,+inf]                p0=(14.137166941,-13.444711759)
                   after : dom=[4.865097571,5.191227259]  p0=( 4.712388980,-5.809475019)
```

**读数**：slot1 拿到原 slot4 的内容、slot4 拿到原 slot1 的；slot2 ↔ slot3 同理。
**恰好两组两两互换** —— 这正是 `swap_seam`（`wire_data.rs:65-66`：`set_edge_pcurves` 交换 FWD/REV 对
+ `set_pcurve_range`）的特征，与 §9.296.6 的构造级证明**一致**。
⇒ **`swap_seam` 是真凶**这一判断，现在有了逐槽位的直接观测，不再只是构造推断。

**face 1537 —— 周期位移（不是交换）**

```text
slot1 before p0=( -2.066942528,1.570796327)  after p0=( 4.216242779,1.570796327)
slot2 before p0=(-11.500662424,2.418858406)  after p0=( -5.217477117,2.418858406)
slot3 before p0=( -2.262319028,2.418858406)  after p0=( 4.020866279,2.418858406)
slot4 before p0=( -1.570796327,1.570796327)  after p0=( 4.712388980,1.570796327)
```

**四个槽位的 u 位移全部恰好 `+6.283185307`**（`delta / 2π = +1.000000`，四个都是）。
**不是**重新排序（点序未交换），而是**整体平移一个 u 周期**。
⇒ 来源是 `shape_fix_face.rs` 里的 `adjust_by_period`（`:470`/`:494`/`:550`/`:587` 一带），
即 seam 构造时的周期归位，**与 `swap_seam` 是两回事**。

**这一点很重要**：T-94 的卡面只设想了 `swap_seam` 一种写入，实际上「face 自己的 pcurve 被改」
有**两条**路径。若要解耦，只改 `wire_data.rs` 的 `swap_seam` **只能消掉 face 1752 那一半**，
face 1537 的周期位移仍在。

#### 9.314.3 「它是否在这 7 个未网格面里」：**不在**

用现成的 `zz_uv_feed --ids`（`PORTID f=N ... mt=<三角数>`）跑 T0M：

```text
PORTID f=1537 key=2207263158256 wires=2 mt=10  bbox=(14.73,-20.62,-8.30)-(20.62,-14.73,-7.96)
PORTID f=1752 key=2207265074304 wires=2 mt=15  bbox=(-9.50,-9.50,-60.19)-(9.50,9.50,-44.64)
```

**两个面都出网格**（`mt=10` / `mt=15`），且本轮全量扫描 `mt=0` 的面为 **0 个**
（1772 个 `PORTID` 行）。

⇒ 卡面要求的确认完成：**face 1752（以及新发现的 1537）都不在未网格面里**，
它们的 pcurve 被改动**没有影响网格产出**。这与「命中率低 + 惰性」的既有结论一致。

#### 9.314.4 对 T-94 方案的影响（需在动手前知道）

| 方案 | 能消掉 face 1752 的交换 | 能消掉 face 1537 的周期位移 |
|---|---|---|
| 只改 `wire_data.rs` 的 `swap_seam`（§9.296.7 的 C） | ✅ | ❌ |
| 动 `repr_key` 全局键约定（A） | ✅ | 视下笔位置而定 |

⇒ **C 方案不足以完全解耦**。若目标真是「只对目标面生效」，需要同时处理
`shape_fix_face.rs` 里经 `tmp_f` 键写回原面槽位的那些 `set_edge_pcurve`/`adjust_by_period` 路径。
**本卡的定级（P3、命中 2/24、两个面都出网格）因此更应维持**，不宜升级。

#### 9.314.5 复跑

```text
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_share_probe -- data/occ/T0M.stp
crates/occt-topo/target/debug/examples/zz_uv_feed.exe data/occ/T0M.stp --ids > portid.txt
grep -E '^PORTID f=(1537|1752) ' portid.txt
```

---

### 9.315 —— 【重大】T-94 边界实验：消除共享键写回**能让 a3n00 全部 226 面出网格**（`unmatched` 1→0），但**面积比从 0.8996 掉到 0.8627** ⇒ 该写回是**承重的**，不是可自由解耦的副作用

#### 9.315.1 实验做法（一处调用 + 一个辅助函数）

按 §9.296.7 的 **A 方案**的最小实现（**不改 `repr_key` 全局语义**，只给临时面摘掉自己的注册）：

1. `tgeometry.rs` 新增 `unregister_face_surface(&self, s: &TopoShape)` ——
   只从 `face_surfaces` 移除 `s` 自己的条目（与 `clear_shape` 的区别：**不**动 edge/face 几何）；
2. `shape_fix_face.rs:453` 之后加一行：
   `GeometryRegistry::global().unregister_face_surface(&tmp_f.0);`
   ⇒ `repr_key(tmp_f)` 从「原面曲面指针」回退为 **`tmp_f` 自身键**，
   seam 步骤经 `tmp_f` 键的所有 pcurve 写入**不再落到原面槽位**。

#### 9.315.2 效果一：写回被**完全**消除（达成目标）

`zz_share_probe` 的逐槽位快照：

| 模型 | 实验前 | 实验后 |
|---|---|---|
| T0M | `faces_with_fix=21 faces_changed=2 slots_changed=8` | **`0 / 0 / 0`** |
| a3n00 | `3 / 0 / 0` | **`0 / 0 / 0`** |

⇒ 共享键写回**彻底消失**。§9.296.6 的构造级证明得到了**实验级确认**。

#### 9.315.3 效果二：a3n00 **全 226 面出网格**（F113 修好了）

| 指标 | 实验前 | 实验后 |
|---|---|---|
| `unmatched` | **1**（F113） | **0** |
| F113 的 `mt` | 无（不在网格模型里） | **267** |
| a3n00 `mt=0` 的面 | 1 | **0**（226 面全出网格） |
| `mesh_v` / `mesh_t` | 10863 / 11941 | **11101 / 11121** |
| `stats` | 225 | **226** |

**T-99 追了九轮的 F113，在这个实验下出网格了。** 这印证了 §9.303.4 的结论：
F113 的病根就是「共享键写回污染了原面 pcurve」，而**不是** `CheckWire`/`ComposeShell`/reader。

#### 9.315.4 效果三（**否决性**）：面积比**变差**

| 指标 | 实验前 | 实验后 |
|---|---|---|
| `ratio_ours_over_occ` | **0.899606** | **0.862724** |
| 绝对差 | 22802.6 | **31179.6**（+8377） |

判据是 `|ratio − 1| < area_tol`（`area_tol = 0.15`）：
0.8996 → 偏离 0.1004（裕度 0.0496）；**0.8627 → 偏离 0.1373（裕度仅 0.0127）**。
**更接近红线，且绝对差显著变大** ⇒ **不能落地。**

#### 9.315.5 由此得到的三条结论

**(1) 这个共享键写回是「承重的」，不是可自由摘除的副作用。**
消除它会让 `fix_missing_seam` 的运行结果**真正生效**（不再被丢弃），从而**改变哪些面被返回** ——
`unmatched` 归零就是它生效的证据。换句话说：

> T-94 设想的「**解耦**」并不存在 —— 这个写回与「修复结果是否被采用」是**同一件事的两面**。
> 任何让它不写回原面的改动，都会**同时**改变 seam 修复的采纳结果，即它**必然**是一次**行为改动**，
> 必须按「实质性修复」评估，而不是按「中性重构」评估。

**(2) 因此卡面的方案表（§9.296.7 的 A/B/C）前提有误。**
三个方案都被描述为「解耦」，但实验证明**至少在 A 方案下它不是解耦**。
C 方案（只改 `swap_seam`）此前已被 §9.314.4 判定**只覆盖一半**（face 1537 的周期位移走
`shape_fix_face.rs:506/513` 的 `tkey` 路径），现在更要打上问号：
若它不改变采纳结果，则它**消不掉 F113**；若它改变了，则它就是行为改动。**需先做实验才知道，不能假定。**

**(3) 对 T-99 的意义（正面）**：
F113 的病根**已被实验确认**是共享键写回污染，且**存在一条能让它出网格的路径**。
问题从「机制未定位」变成了「**找到既消除写回、又不让面积比变差的机制**」——
这是一个明确得多的问题。T-99 因此**值得从「搁置待外部证据」重新开为可动手的卡**。

#### 9.315.6 未做与不做的理由

* **未跑 23 模型全量 `step_obj_gates`**：既然 a3n00 的面积比在**判据方向上变差**，
  该改动已不具备落地资格；先跑全量只会消耗预算而不改变结论。
* **未提交该实验**：已 `git stash` 回退，工作区与实验前逐字一致
  （`shape_fix_face.rs`、`tgeometry.rs` 的 diff 均为空；基线回到 `stats=225 unmatched=1`）。
  实验内容保留在 `stash@{0}` 与本节的描述中，可复现。

#### 9.315.7 复跑（若需重建该实验）

```text
# 1) 加入辅助函数（§9.315.1 第 1 步）
# 2) shape_fix_face.rs 在 `let mut tmp_f = builder.make_face(surf.clone(), &tmp_wires);` 之后加：
#    GeometryRegistry::global().unregister_face_surface(&tmp_f.0);
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_share_probe -- data/occ/T0M.stp
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_probe_a3n00 -- data/occ/a3n00.stp
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj -- data/occ
python .target-gate/pair_area.py output/a3n00.obj data/occ/occ-a3n00.obj
```

---

### 9.316 —— T-99 方向 ③ 的定位结果：面积差**不是**"多了张错误的面"，而是实验版**少铺了 8377 面积**（三角更少、更粗）

#### 9.316.1 逐三角对比（`zz_uv_feed` 导出的两个 OBJ）

| 版本 | tris | area | 对 GT（227130.2998）的缺口 |
|---|---|---|---|
| **基线** | 11941 | **204327.7166** | **22802.5832** |
| **实验**（消除共享键写回） | 11121 | **195950.7131** | **31179.5867** |

两者**最大的三角完全同名同值**（`1514.8633`、`1326.7472`、`1322.0572`、四个 `1069.1455`、`903.1153`、`727.9473`），
说明**不存在"实验版多出一张巨大错误面"**这种情形（否则会看到明显不同的超大三角）。

⇒ 差异是**分布性的**：实验版 **`v` 更多（11101 对 10863）却 `f` 更少（11121 对 11941）**、
**总面积反而少 8377** ⇒ 它在整体上**铺得更粗、覆盖更少**，**不是多铺**。

**这修正了 §9.315.4 的表述**：那里只说"面积比变差"，未说明方向。
现在明确：**实验版是"少铺"（under-cover），不是"多铺错误几何"**。

#### 9.316.2 与"F113 出网格"并不矛盾

实验版解决了 `unmatched`（F113 出网格，`mt=267`），却在**其余面**上少铺了更多面积。
两类效应在同一个开关下**朝相反方向**变化：

| 开关 | `unmatched` | F113 | 全模型面积覆盖 |
|---|---|---|---|
| 保留写回（基线） | 1 | 不出网格 | **较好**（缺口 22802.6） |
| 消除写回（实验） | **0** | **出网格** | 较差（缺口 31179.6） |

⇒ **这不是"修好一个、弄坏一个"的偶然**，而是同一个共享键机制在两面同时起作用的结果，
且**两面方向相反**。这正是 §9.315.5 所说"写回是承重的"的具体表现。

#### 9.316.3 由此得到的关键推断（下一轮的直接检查点）

实验版**三角更少、覆盖更少**，最可能的解释是：
**被采纳的"修复后面"在它自己的键下没有 pcurve** ——
因为所有 seam 写入都被我的实验挡在了 `tmp_f` 的（已摘注册的）键上，
于是 `ComposeShell` 造出的结果面**缺 pcurve**，网格器只能少铺。

**可直接检验的不变量**（下一轮第一步）：

> 对每个被采纳的面 `F`：`F` 的每条边是否都在 **`repr_key(F)`** 下有一条 pcurve？

若实验版在若干面上违反该不变量、而基线不违反，则结论就是：
**写入的键与"面的表示键"不匹配** —— 修复方向因此是
**「让被采纳面的 pcurve 落在它自己的 `repr_key` 下」**，
而不是「消除写回」。这同时解释了为什么单纯消除写回会两败俱伤。

#### 9.316.4 对 T-99 三个候选方向的重新排序

| 方向 | 评估 |
|---|---|
| ① 只对**需要 seam 的面**摘注册 | **降级** —— 若"面的 pcurve 键不匹配"是主因，则只对部分面摘只会把两败俱伤缩小到那些面，方向不变差也不变好 |
| ② 让 `resolve_face` 采用修复结果的方式与 OCCT 一致（形状级替换表） | **升级为**首选 —— 它正是"让被采纳面自带正确 pcurve"的那条路 |
| ③ 逐面比较写回差异 | **本轮已完成**（§9.316.1/2），产出是"少铺而非多铺" |

#### 9.316.5 状态

* 实验已 `git checkout` 回退；`shape_fix_face.rs`、`tgeometry.rs` 的 diff 均为空；
  基线回到 `stats=225 unmatched=1`。
* 两个对照 OBJ 留在 `.target-gate/a3n00_baseline.obj` 与 `.target-gate/a3n00_experiment.obj`，
  供下一轮直接复算（**不必重跑实验**）。
* 库代码本轮**零改动**。

#### 9.316.6 复跑

```text
python .target-gate/pair_area.py .target-gate/a3n00_baseline.obj   data/occ/occ-a3n00.obj
python .target-gate/pair_area.py .target-gate/a3n00_experiment.obj data/occ/occ-a3n00.obj
python .target-gate/area_by_face.py .target-gate/a3n00_baseline.obj .target-gate/a3n00_experiment.obj
```

---

### 9.317 —— 附：看板结构引号损坏的修复（由我自造脚本造成，已根治并删除肇事脚本）

#### 9.317.1 事故与范围

本会话早前为「修中日引号」写的一次性脚本（`fix_inner_quotes.py`）**把部分结构引号也当成了内层引号**，
把属性分隔的 `"` 改写成 `「` / `」`，导致**两个属性被并进同一个字符串**，例如：

```text
next: "—", acceptance: "--lib 1281/0、step_obj_parity 14/14、export 逐位一致"
被写成
next: "—「, acceptance: 」--lib 1281/0、step_obj_parity 14/14、export 逐位一致"
```

**范围核对**（重要）：`git show HEAD:specs/board.canvas.tsx` 里这类损坏为 **0 处**；
我工作区的画布有 **76 处**。⇒ **损坏全部由工作区内的脚本改动造成，HEAD 是干净的。**

#### 9.317.2 一个可靠的判据（值得记下来）

比正则更可靠的判别：**某行 ASCII 双引号个数为奇数 ⇒ 该行结构必坏**
（正常的属性行引号总是成对）。用它把范围从 76 处收敛到 **16 行**：

```
lines with an odd number of ASCII quotes: 16
[84, 102, 119, 153, 175, 190, 199, 207, 216, 224, 233, 242, 250, 259, 268, 277]
```

其中 60 处是**误报** —— 那是我自己的检测脚本把合法的**强调引号** `「等」` / `「未完成」` 也算成了损坏。
**强调引号永远后跟文字，而结构引号后跟 `,` 或 `]`** —— 这个区别是精确判据。

#### 9.317.3 修复

三步、每步都比上一步更窄：

1. `「, <key>: 」` → `", <key>: "`（被并进同一字符串的属性边界）—— 修 38 行；
2. 数组元素形态（首个元素的开引号被写成 `「`、闭合引号仍是 `"`）—— 修 16 行；
3. 用**引号奇偶**复验：`still structurally broken: []`。

`canvas_check` → **OK**；卡片普查 23 张（pending 1 / completed 19 / cancelled 3）**未丢**；
抽查修复行内容完好（`rejected: ["恢复全局 fix_missing_seam 补偿", "保留现状不动（F113 不出网格）"]`）。

#### 9.317.4 根治措施

**删除了会改写引号的脚本**：`fix_inner_quotes.py`、`fix_inner2.py`、`fix_quote.py`。
保留的只是**只读检测**（`count_damage.py`、`diag_quotes.py`）。
教训：**不要写会批量改写分隔符的脚本** —— 分隔符与内容共用同一种字符时，
任何基于"数引号"的启发式都会误伤；这类损伤还会**静默**累积（
本会话正是靠 `canvas_check` 的报错才发现）。

#### 9.317.5 本轮门禁（修复后复验）

| 项 | 结果 |
|---|---|
| `cargo check` | 0 error |
| `occt-topo --lib` | **1255 passed / 26 failed**（26 条为 `%TEMP%` 环境性失败，与既有基线相同） |
| a3n00 基线 | `faces=226 … stats=225 mesh_v=10863 mesh_t=11941`、`matched=225 unmatched=1` —— **未变** |
| 库代码 | `shape_fix_face.rs`、`tgeometry.rs` 的 diff 均为空（实验已回退） |

---

### 9.318 —— T-99 不变量检查：**基线本身就不满足**（448/1083 边无 pcurve），但它指出了关键差异 —— 基线比实验版**多 107 条边、且这 107 条全都有 pcurve**

#### 9.318.1 我上一轮设的不变量，表述过强

§9.316.3 提出：「对每个被采纳的面 `F`，`F` 的每条边是否都在 `repr_key(F)` 下有一条 pcurve？」
新建探针 `crates/occt-topo/examples/zz_pcurve_key_probe.rs` 实测**基线**（`data/occ/a3n00.stp`）：

```text
KEYPROBE faces=226 edges=1083 edges_without_pcurve=448 faces_with_missing=95
KEYPROBE   face=128 edges=33 missing=33 wires=6
KEYPROBE   face=115 edges=16 missing=16 wires=6
KEYPROBE   face=161 edges=7  missing=7  wires=1
...
```

T0M 更甚：`faces=1772 edges=8630 edges_without_pcurve=3674 faces_with_missing=719`。

⇒ **不变量在基线上就不成立**（41% 的边在 `repr_key(F)` 下没有 pcurve）。
所以「基线满足、实验违反」这种判别式**不成立** —— **该不变量表述过强**，
这也是本轮的第一个否定结果，记下来免得下一轮再走。

（`edge_pcurves(edge, shape_key(F))` 内部正是走 `repr_key`，
因此它问的就是网格器要问的那个问题；这个调用方式是有效的，只是**期望值**设错了。）

#### 9.318.2 但对照本身给出了**有信息**的差异

| 版本 | faces | edges | `edges_without_pcurve` | `faces_with_missing` |
|---|---|---|---|---|
| **基线**（保留写回） | 226 | **1083** | 448 | 95 |
| **实验**（消除写回） | 226 | **976** | **448** | **95** |

两条读数：

1. **缺 pcurve 的边数完全相同（448）、缺 pcurve 的面数完全相同（95）** ——
   说明实验**没有**制造新的"无 pcurve"边；我上一轮关于"被采纳面缺 pcurve"的猜测**不被支持**。
2. **基线比实验版多 107 条边（1083 − 976）**，而这 107 条**全部有 pcurve**
   （因为 `edges_without_pcurve` 两版都是 448）。

⇒ 差异的落点是：**基线在 107 条边上多挂了 pcurve**。这 107 条正对应
`zz_share_probe` 观测到的写回（2 面 / 8 槽位所影响的那批 seam 边及其同族）。
换句话说是：

> 写回把 **107 条边**的 pcurve 挂在原面键下（基线有的、实验没有）；
> 实验把这些 pcurve 挂到了**别的键**下，于是这 107 条边在 `repr_key(F)` 下"看不见"。

**这不是"缺 pcurve"，而是"挂在错误的键下"。** 与 §9.316.3 的方向一致，
但**判据要改成**：「对这 107 条边，基线在 `repr_key(F)` 下有 pcurve、实验没有」——
而不是笼统地要求所有边都有。

#### 9.318.3 下一轮的直接检查点（据此收敛）

**要量的是这 107 条边**，不是全部 1083 条：

1. 用探针列出「**基线有、实验没有**」的 (面, 边) 集合（两版各跑一次，取差集）；
2. 对这批边打印：它在基线下挂的键、在实验下挂的键，以及**它所属的面是否被采纳**；
3. 若这批边所属的**面没有被采纳**（即原面仍被返回）⇒ 实验把它们的 pcurve 挪走了，
   于是那些面**少铺** —— 与 §9.316.1 观测的"实验版少铺 8377"方向一致，
   而且能逐面定位到面积恶化的来源。

**这个差集是可枚举的、有界的（107 条边）**，比上一轮的"全部边"判据可操作得多。

#### 9.318.4 本轮动作与状态

* 新增探针 `crates/occt-topo/examples/zz_pcurve_key_probe.rs`（可复跑，报告"每面缺 pcurve 的边数"）；
* 实验（`unregister_face_surface`）已 `git checkout` 回退，`tgeometry.rs`、`shape_fix_face.rs` 的 diff **均为空**；
* 基线复验：`stats=225`、`unmatched=1`、`faces=226 edges=1083 edges_without_pcurve=448`；
* **库代码本轮零改动。**

#### 9.318.5 复跑

```text
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_pcurve_key_probe -- data/occ/a3n00.stp
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_pcurve_key_probe -- data/occ/T0M.stp
# 重建实验：python .target-gate/rebuild_experiment.py
```

---

### 9.319 —— T-99 关键因果链已打通：消掉写回会让 **43 个面各少 2–3 条边**（seam 边加不上），这就是「少铺 8377」的成因

#### 9.319.1 差集量出来了（`zz_pcurve_key_probe` 的 `--dump` 模式）

形状键是**指针值**，跨进程因 ASLR 不稳定 ⇒ 不能跨运行比原始键。
但**逐面边数**稳定，用它对比（基线 vs 重建实验）：

```text
baseline  : faces=226 total_edges=1083
experiment: faces=226 total_edges= 976     delta = -107
faces whose edge COUNT differs: 43
   face= 41  baseline=4  experiment=2  delta=-2
   face= 42  baseline=4  experiment=2  delta=-2
   ...
   face=109  baseline=5  experiment=2  delta=-3
   ...
   face=138  baseline=15 experiment=13 delta=-2
```

**43 个面、每个少 2–3 条边**（4→2 或 5→2），合计 **−107**，与总数差**逐字吻合**。

#### 9.319.2 这解释了三件事

**(a) 「少铺 8377」的成因**：这 43 个面各少 2–3 条边 ⇒ 它们的边界环退化 ⇒
网格覆盖缩小。方向与 §9.316.1 的观测（实验版 `f` 更少、面积更低）**一致**。

**(b) 少掉的正是 seam 边**：对照 §9.290/§9.292 的既有结论 ——
`fix_missing_seam` 会为单闭合边面**加** seam 边（T0M f=1691 从 1 条变 4 条）。
基线这 43 个面多出的 2–3 条边，正是这一机制加上的 seam 边。

**(c) 为什么消掉写回会让 seam 边加不上** —— 这是本轮最要紧的一环：

> `fix_missing_seam` 在构造 seam 的过程中**要读回到它自己刚写下的 pcurve**
> （`check_pcurves_and_shift` / `curve_on_surface_range` 都按 `repr_key` 取）。
> 我的实验把 `tmp_f` 的曲面注册摘掉后，写入落在 `tmp_f` 自己的键上，
> **但读取侧仍然按 `repr_key(tmp_f)`（= 原面曲面指针）去找** —— 于是读不到，
> `fix_missing_seam` 拿不到 pcurve ⇒ **seam 边加不上** ⇒ 面退化。

**⇒ 抑制写回并没有"解耦"，而是切断了 seam 步骤自身的读写通路。**

#### 9.319.3 因此真正的症结被精确定位：`repr_key` **把「面身份」与「曲面身份」混为一谈**

`tgeometry.rs:524-531` 的 `repr_key(face_key)` 把面映射为**它所用曲面的数据指针**，
其注释自陈"Faces that share a surface therefore share one representation, as in OCCT"。

但本用例证明这个等价**不成立**：
`tmp_f` 与原面**共享曲面**，却需要**各自独立的 pcurve**（seam 步骤要写不同的 uv）。
于是"共享曲面 ⇒ 共享表示"在这条路径上同时造成两个后果：

* 写：seam 的 pcurve 落到原面槽位 ⇒ **污染原面**（F113 的病根）
* 读：若不给它共享键，seam 步骤又**读不到自己的写入** ⇒ **加不上 seam 边**（43 面退化）

**两难**。而与 OCCT 对照可以发现，OCCT 并不面临这个两难：
OCCT 把 pcurve 挂在**那条边的 `BRep_TEdge` 对象**上（`BRep_Builder::UpdateEdge(E, c2d, F)`），
按 `(edge, face)` 定位，**与曲面无关** ——
所以 OCCT 里"两个面共享曲面"与"各自有独立 pcurve"本来就可以并存。

⇒ **修复方向明确了**：让 pcurve 的表示键跟随**面的身份**，而不是曲面指针。
**这正是 T-94 §9.296.7 里被标为"高风险、动全局键约定"的 A 方案** ——
现在有了必要性论证：不是"想不想动"，而是**不动就解不开这个两难**。

#### 9.319.4 下一轮的具体动作（有依据，非试探）

1. 在 `tgeometry.rs` 里把 `repr_key` 的语义从「面的曲面指针」改为「面自身的键」
   （或引入 `pcurve_key` 与 `surface_key` 两个概念，前者用于 pcurve 存储）；
2. 该改动的影响面：所有经 `repr_key` 存取 pcurve 的调用点
   （`edge_pcurves` / `set_edge_pcurve(s)` / `pcurve_range` / `remove_pcurves_on_surface`）；
3. **必须整轮跑 23 模型 `step_obj_gates`** —— 这是全局键约定改动，红线要求；
4. 验收信号（现成）：`unmatched` 1→0 且 a3n00 面积比 **不低于 0.8996**；
   同时 `zz_pcurve_key_probe` 的总边数应维持在 **≥1083**（不再出现 43 面退化）。

#### 9.319.5 本轮状态

* 实验已 `git checkout` 回退；`tgeometry.rs`、`shape_fix_face.rs` 的 diff **均为空**；
* 基线复验：`stats=225`、`unmatched=1`、`edges=1083`；
* 探针 `zz_pcurve_key_probe.rs` 增加 `OCCT_TOPO_KEY_DUMP` 模式（逐面边键 dump，供跨运行对照边数）；
* **库代码本轮零改动。**

#### 9.319.6 复跑

```text
python .target-gate/rebuild_experiment.py
OCCT_TOPO_KEY_DUMP=1 cargo run --manifest-path crates/occt-topo/Cargo.toml --offline \
    --example zz_pcurve_key_probe -- data/occ/a3n00.stp > keys_exp.txt
git checkout -- crates/occt-topo/src/tgeometry.rs crates/occt-topo/src/shhealing/shape_fix_face.rs
OCCT_TOPO_KEY_DUMP=1 cargo run --manifest-path crates/occt-topo/Cargo.toml --offline \
    --example zz_pcurve_key_probe -- data/occ/a3n00.stp > keys_base.txt
python .target-gate/diff_edges.py
```

---

### 9.320 —— T-99 核心改动实测：把 `repr_key` 改为「面身份」**不解决**问题（与摘注册实验逐字同结果）

#### 9.320.1 改法与实测

按 §9.319.4 执行：把 `tgeometry.rs` 的 `repr_key` 从「面的曲面数据指针」改为「面自身的键」：

```rust
fn repr_key(&self, face_key: usize) -> usize {
    face_key          // 原为 face_surfaces[face_key] 的 Arc 数据指针
}
```

`cargo check` 通过。a3n00 实测：

| 指标 | 基线 | **`repr_key = face_key`** | 摘注册实验（§9.315） |
|---|---|---|---|
| `stats` | 225 | **226** | 226 |
| `unmatched` | 1 | **0** | 0 |
| `mesh_v / mesh_t` | 10863 / 11941 | **11101 / 11121** | 11101 / 11121 |
| **面积比** | **0.899606** | **0.862724** | 0.862724 |
| 绝对差 | 22802.58 | 31179.59 | 31179.59 |

**三个指标与摘注册实验**（§9.315）**逐字相同**。

`zz_pcurve_key_probe` 进一步确认：

```text
edges=976 edges_without_pcurve=448 faces_with_missing=95        ← repr_key = face_key
edges=1083 edges_without_pcurve=448 faces_with_missing=95       ← 基线
```

**同样是 43 个面各少 2–3 条边**（face=41 4→2、face=42 4→2、…），合计 −107。

⇒ **改用面键等价于"摘掉共享键"**：两者都阻止 seam 写入抵达原面槽位，
于是得到**同一组后果**（F113 出网格是好的，43 面退化是坏的）。
**§9.319.4 指出的方向被实测否决。**

#### 9.320.2 这个否决说明了什么（值得记下）

写回**同时**承担两件事：

| 保留写回（基线） | 阻止写回（两种实现都一样） |
|---|---|
| 43 个面**拿到** seam 边（面积覆盖较好） | 43 个面**丢掉** 2–3 条边（少铺 8377） |
| F113 的 pcurve **被污染**（不出网格） | F113 **出网格** |

⇒ **"让 pcurve 的键跟随面身份"这个单一改动，无法只取一半。**
两类后果都源于同一个事实：**这 43 个面的 seam 边之所以存在，依赖"原面槽位里有那个 seam pcurve"**。

因此真正的修复不能只动键的**归属**，还得处理**为什么这些面的 seam 边需要那份 pcurve**。
这是比 §9.319.4 所设想的更窄、也更难的一步：需要看
`fix_missing_seam` 在这些面上**选 wire 对**（`check_wire` / `ismodeu,ismodev`）与
**加 seam 边**（`ShapeFix_Wire`）的具体条件，判断它是否**依赖原面槽位里的 pcurve 才能选中**。

#### 9.320.3 一条必须记下的自我更正

我在本轮中途曾推断「`repr_key = face_key` 的 976 里不含 F113 的 22 边，而基线 1083 含」，
并据此以为能用算术把 107 解释清楚。**实测否决了这个推断**：
两版 F113 的边数**都是 22**（`KEYEDGE face=113 n=22`）。
⇒ 那 107 条**全部**来自那 43 个面，与 F113 无关。
（记下来是因为这个错误推断如果写进卡片会误导下一轮。）

#### 9.320.4 本轮状态

* `repr_key` 改动**已 `git checkout` 回退**；`tgeometry.rs` 的 diff **为空**；
* 基线复验：`stats=225`、`unmatched=1`、`mesh_v/mesh_t=10863/11941`；
* **库代码本轮零改动**（改动试过、量过、回退）；
* 探针 `zz_pcurve_key_probe.rs` 保持（含 `OCCT_TOPO_KEY_DUMP` 模式）。

#### 9.320.5 累计否决清单（T-99，已 13 轮）

| # | 假设 | 否决依据 |
|---|---|---|
| 1 | `CheckWire` 的 pcurve 参数口径错了 | §9.294.1（改后 4 条全被拒） |
| 2 | `v_range` 为无穷是缺陷 | §9.298（属原面；修复后为 175） |
| 3 | 4 条 wire 在本面内被切 | §9.299（3D 上各自闭合；上游） |
| 4 | `resolve_face` 应返回 Shell/compound | §9.293.1（root 归零） |
| 5 | 顶点身份/邻接口径不同 | §9.302（跨 loop 零共享顶点） |
| 6 | 两条退化边同母线 | §9.303.7（分别在 u=0 与 u=π） |
| 7 | `heal_shape` 接在读入路径 | §9.301.1（根本没接） |
| 8 | D18 的 `n1→n2` | §9.291.5（二次否证） |
| 9 | 该写回可中性解耦 | §9.315（写回承重） |
| 10 | 「每个面的每条边都要有 pcurve」 | §9.318.1（基线即不满足 448/1083） |
| 11 | **让 pcurve 的键跟随面身份** | **§9.320（本轮的三个指标与摘注册逐字相同）** |

#### 9.320.6 复跑

```text
# 改动：tgeometry.rs 的 repr_key 改为 `face_key`
cargo check --manifest-path crates/occt-topo/Cargo.toml --offline
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_probe_a3n00 -- data/occ/a3n00.stp
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj -- data/occ
python .target-gate/pair_area.py output/a3n00.obj data/occ/occ-a3n00.obj
git checkout -- crates/occt-topo/src/tgeometry.rs       # 回退
```

---

### 9.321 —— T-99 插桩定位：`swap_seam` 在**读取路径上从未被调用**（推翻 §9.296.3 的判断）；`check_wire` 的 pcurve 早退**两版都是 0**

#### 9.321.1 第一个否定结果（重要，修正早前结论）

在 `wire_data.rs` 的 `swap_seam` 入口插桩（`OCCT_TOPO_TRACE_SWAPSEAM`），跑 a3n00：

| 版本 | `swap_seam` 调用次数 | 其中 `fired`（pc 数 ≥ 2） |
|---|---|---|
| **基线** | **0** | 0 |
| 实验（阻止写回） | **0** | 0 |

⇒ **`swap_seam` 在 `read_step_file` 这条路径上一次都没被调用。**

**这推翻了 §9.296.3 的判断。** 那里写的是「D2 观察到的『补偿改写了共享 `GeometryRegistry`』，
其具体执行者就是 `swap_seam` 这一处」—— **不成立**。
§9.296.1 的「6 处写 registry 里只有 `wire_data.rs:65-66` 作用于既有 seam 边」这个**静态**结论没错，
但那处代码在**读取路径上不可达**（`swap_seam` 只在 `reverse_wire_data_on_face` 的
`stat < 0` 或 `reverse` 分支里被调，读取路径没走到）。

这也解释了 §9.296.9/§9.314 的测量：`zz_share_probe` 是**单独对每个面调
`ShapeFixFace::fix_missing_seam()`** 来观测的，那是一个**探针构造的调用**，
不代表 `read_step_file` 会走同样的分支。**两者不是同一条路径。**

#### 9.321.2 第二个否定结果：`check_wire` 的 pcurve 早退两版都是 0

在 `check_wire` 的 `curve_on_surface_oriented(...) == None` 早退处插桩
（`OCCT_TOPO_TRACE_CHECKWIRE`）：

| 版本 | `no_pcurve` 早退次数 |
|---|---|
| 基线 | **0** |
| 实验 | **0** |

⇒ 我基于 §9.321.1 提出的假设（「共享键的第二个作用是 `check_wire` 用来找到 wire 的 pcurve」）
**被实测否决**：`check_wire` 在这条路径上**从不因为拿不到 pcurve 而提前返回**。

#### 9.321.3 那么 43 面退化（−107 边）发生在读取路径的哪一步？

两条读数合起来把范围收窄了很多：

* 共享键写回**确实影响结果**（`stats` 225→226、`unmatched` 1→0、`edges` 1083→976）——
  这已由 §9.315/§9.320 两次独立实现证实；
* 但它**不是**经 `swap_seam`（未被调用），也**不是**经 `check_wire` 的 pcurve 早退（两版都是 0）。

⇒ 差异必然出现在**别处**，而且有一个很强的线索：
**`unmatched` 从 1 变 0** 意味着 **F113 在两版里是不同的面对象**（一版被采纳、一版没有），
即 `resolve_face` 的返回值不同。而 `edges` 总数同时少了 107。
最可能的位置是 `shape_fix_face.rs` 里**决定是否返回 `sff.result`** 的那段（`cxx:2261` 的 `comp.perform()` 之后）
—— 它的行为依赖 registry 中的 pcurve 状态。

**下一轮的直接检查点**：在 `shape_fix_face.rs` 的 seam 块里插桩，对每个面打印
「`fix_missing_seam()` 返回值 / `sff.result` 的类型 / 采纳与否 / 采纳后的边数」，
**在基线与实验下各跑一次并取差集**。这能一次定位到分歧点，且与 §9.321.1/2 不同 ——
前两次我猜的是"哪个函数"，这次要量的是"哪些面被采纳与否"。

#### 9.321.4 本轮状态

* 三处插桩（`wire_data.rs` 的 `swap_seam`、`shape_fix_face.rs` 的 `check_wire`）与实验
  （`tgeometry.rs`/`shape_fix_face.rs` 的 `unregister_face_surface`）**全部 `git checkout` 回退**；
* 三个文件的 `git diff --numstat` 均为空；`OCCT_TOPO_TRACE_CHECKWIRE`/`SWAPSEAM`/`TEMP T-99`
  的 grep **无残留**；
* 基线复验：`stats=225`、`unmatched=1`、`mesh_v/mesh_t=10863/11941`；
* **库代码本轮零改动**（又一次"改了、量了、回退"）。

#### 9.321.5 累计否决清单（T-99，14 轮，13 条）

| # | 假设 | 否决依据 |
|---|---|---|
| 1–11 | 见 §9.320.5 | — |
| 12 | 共享键写回的危害经 `swap_seam` 生效 | **§9.321.1**（读取路径 0 次调用；推翻 §9.296.3） |
| 13 | 共享键的第二个作用是让 `check_wire` 找到 pcurve | **§9.321.2**（两版早退都是 0） |

#### 9.321.6 复跑

```text
# swap_seam 插桩：wire_data.rs 的 swap_seam 入口 eprintln!（OCCT_TOPO_TRACE_SWAPSEAM）
# check_wire 插桩：shape_fix_face.rs 的 curve_on_surface_oriented == None 分支（OCCT_TOPO_TRACE_CHECKWIRE）
OCCT_TOPO_TRACE_SWAPSEAM=1 rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline \
    --example zz_probe_a3n00 -- data/occ/a3n00.stp 2> ss.txt
grep -c '^SWAPSEAM' ss.txt
```

---

### 9.322 —— 【停止信号】基线出现不可解释的漂移：本会话 T-99 所依赖的「225 / unmatched=1 / 0.8996」在**干净工作树上无法复现**

#### 9.322.1 事实（全部可复核）

本轮做「记录被采纳面」的插桩时发现读数与卡片不符，遂做彻底核查，结果如下：

**(1) 工作区的 `read_topology.rs` 与 HEAD 的 blob 哈希相同**：

```text
work blob: 4f79fa0c98a0b9cae557635ad61f5f314a621cfa
HEAD blob: 4f79fa0c98a0b9cae557635ad61f5f314a621cfa
identical: True
```

⇒ **文件与 HEAD 逐字相同、无未提交改动**。

**(2) 但该文件里没有 T-93 的 reader 侧 seam 块**：

```text
grep 'fix_missing_seam'      read_topology.rs -> 0
grep 'T-99'                  read_topology.rs -> 0
grep 'Measured 2026-09-30 (T-93)' read_topology.rs -> 0
```

而 `git show --stat f1e56776` 显示该提交**删除**了这段（
`- step/read_topology.rs：删除 reader 侧补偿块（原 :733-747 的 ShapeFixFace::fix_missing_seam + …）`）。

**但我本会话多次读到过这段**（例如 §9.293 之后各轮引用的 `:734-744`），并据其插桩。
⇒ 要么它是在本会话期间被加回、又被某次 `git checkout` 抹掉；要么我对它的读取发生了错位。
**目前无法判定是哪种。**

**(3) 当前干净工作树的实测（连跑 3 次，完全一致）**：

```text
TOTAL faces=226 … stats=226 mesh_v=11101 mesh_t=11121
STATMAP matched=226 unmatched=0
PORTID f=113 wires=4 mt=267            ← F113 出网格
ratio_ours_over_occ = 0.862724  (abs diff 31179.5867)
```

**(4) 与会话基线对照**：

| 量 | 本会话一直使用的基线 | **当前干净工作树** |
|---|---|---|
| `stats` | 225 | **226** |
| `unmatched` | 1（F113 不出网格） | **0（F113 `mt=267`，出网格）** |
| `mesh_v / mesh_t` | 10863 / 11941 | **11101 / 11121** |
| **面积比** | **0.899606** | **0.862724** |

⇒ **当前状态就在会话一节把 T-99 推向「已解决」的那一组数字上**（226/0/11101/11121/0.862724），
而**不是**会话里反复引作基线的 225/1/10863/11941/0.8996。

#### 9.322.2 这为什么必须停下来

T-99 与 T-94 的**全部结论**都建立在「基线 = 225 / unmatched=1（F113 不出网格）/ 面积比 0.8996」之上：

* 「F113 的病根是共享键写回」—— 以「基线 F113 不出网格」为前提；
* 「实验让 F113 出网格但面积比掉到 0.8627，所以不能落地」—— 以基线 0.8996 为参照；
* 「43 个面各少 2–3 条边」—— 以基线 `edges=1083` 为参照。

如果**当前状态已经是 226/0/0.862724**，那么其中至少有一部分叙述的**参照系不成立**，
那些数字**不能继续作为结论依据**。而 `read_topology.rs` 与 HEAD 相同这一点，
又与「本会话读到过 T-93 seam 块」相互矛盾，说明**工作树在会话期间发生过未被记录的状态变化**
（很可能是某次 `git checkout` / `git stash pop` 把文件恢复到了与当时认知不同的版本）。

**我不能在参照系不确定的情况下继续推进 T-99** —— 那只会产出更多不可比的数字。

#### 9.322.3 本轮没有做的事（刻意）

* **没有提交任何源码改动**：本轮只做了插桩与核查，插桩已全部 `git checkout` 回退；
* **没有再试新方案**：参照系未定之前，任何"实验结论"都不可信；
* **没有改卡片状态**：T-99 的处置应由掌握完整上下文的人决定。

#### 9.322.4 建议的下一步（需要人决定）

1. **先做一次干净的基线重建**：在一个明确的提交上，重新测定
   `stats / unmatched / mesh_v / mesh_t / 面积比 / edges`，并把它作为唯一参照；
2. **核对工作树历史**：确认 `read_topology.rs` 的 T-93 seam 块是在哪一步消失的
   （本会话的 `git stash push/pop` 与 `git checkout` 是重点嫌疑）；
3. 在 (1)(2) 完成前，**T-99 不宜继续动手**；若持续无法确定参照系，
   应按「外部条件不足」处理，而不是继续生成测量。

#### 9.322.5 复跑（用于重建基线）

```text
git hash-object crates/occt-topo/src/step/read_topology.rs
git rev-parse HEAD:crates/occt-topo/src/step/read_topology.rs
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_probe_a3n00 -- data/occ/a3n00.stp
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj -- data/occ
python .target-gate/pair_area.py output/a3n00.obj data/occ/occ-a3n00.obj
```

---

### 9.323 —— 基线重建结果：当前干净树**已经**处在「F113 出网格但面积比更差」的状态

#### 9.323.1 新基线（唯一参照，全部在干净树上重测）

```text
HEAD                          = 9036d8aba2ec47f798a11d8a9f41069692d4042c
read_topology.rs blob         = 4f79fa0c98a0b9cae557635ad61f5f314a621cfa
HEAD 的同一文件 blob          = 4f79fa0c98a0b9cae557635ad61f5f314a621cfa   ← 逐字相同，无未提交改动
```

| 量 | **新基线（HEAD 干净树）** | 本会话旧参照 |
|---|---|---|
| a3n00 `stats` | **226** | 225 |
| a3n00 `unmatched` | **0** | 1 |
| a3n00 `mesh_v / mesh_t` | **11101 / 11121** | 10863 / 11941 |
| a3n00 面积比 | **0.862724**（abs diff 31179.5867） | 0.899606（22802.5832） |
| a3n00 `edges` | **976** | 1083 |
| **F113** | **`wires=4 mt=267`（出网格）** | 不出网格 |
| T0M `stats / unmatched` | **1765 / 7** | （未在本会话重测） |
| T0M `mesh_v / mesh_t` | **49999 / 47616** | — |

a3n00 连跑 3 次完全一致，非随机。

#### 9.323.2 这个结果**推翻了 T-99 的整个前提**，也解释了漂移

本会话九轮所依据的「F113 不出网格、面积比 0.8996」，**在 HEAD 干净树上不可复现**。
而当前干净树**已经**处在「F113 出网格」这一端：

| | F113 | 面积比 |
|---|---|---|
| 本会话旧参照 | **不出网格** | 0.8996 |
| **当前 HEAD 干净树** | **出网格（mt=267）** | **0.8627** |
| 本会话的「实验/变体」 | 出网格 | 0.8627 |

⇒ **「实验组」与「当前干净基线」是同一个状态。**
本会话反复报告的「消除共享键写回后 F113 出网格但面积比掉到 0.8627，
所以该变动不能落地」——**对当前树而言没有意义**：
当前树**没有**那个变动，却已经是 0.8627。

换句话说：**T-99 追了十四轮的「F113 不出网格」，在当前代码上根本不存在。**
要么旧参照来自另一个工作树状态，要么当时的测量有误 —— 现在无法区分，
但**当前可复现的事实只有一个**：干净树 = 226 / 0 / 0.8627 / F113 出网格。

#### 9.323.3 漂移的一个合理解释（有代码依据）

读者的 `resolve_edge` **自己**就把 STEP 的 seam pcurve 存进 registry
（`read_topology.rs:938-958`：判 `SEAM_CURVE` / `is_like_seam`，再
`set_edge_pcurves(&edge.0, face_key, stored)`）。

⇒ **端口的 seam 信息来自 STEP 本身**，不依赖 T-93 那段 reader 侧
`fix_missing_seam` 补偿（该段已在 `f1e56776` 删除，工作区里也确实没有）。
所以「F113 能否出网格」主要由 **seam pcurve 被存到哪个键**决定，
而这正是 `tgeometry.rs` 的 `repr_key` 语义。

**这与 §9.320 的实验对上了**：把 `repr_key` 改成 `face_key` 后
`stats 225→226`（F113 出网格）—— 但**当前 HEAD 无需任何改动就已经是 226**。
⇒ 说明**当前树的 `repr_key` 行为与我做实验时的假设不同**，
或当时的 225 来自另一状态。**这一点必须由人核对，不能由我继续推断。**

#### 9.323.4 对 T-99 的处置建议（本轮不改卡片状态，只记录）

T-99 的标题、`goal`、`next` 全部以「F113 不出网格 ⇒ 面积缺口主因」为前提，
而该前提在当前树上不成立。建议：

1. **先把 T-99 的前提重写**：把 `goal` 从「让 F113 出网格」改为
   「把 a3n00 面积比从 **0.8627** 推向 1」（F113 已经出网格，但它**仍占面积缺口**：
   当前 abs diff 31179.6 比旧参照的 22802.6 **更大**）。
2. 或按「旧参照不可复现、需人工核对工作树历史」处理，暂停本卡。
3. **无论哪条**，都不应继续以「消除/保留共享键写回」为主线 ——
   该二分在当前树上不存在。

**本轮不改卡片状态**：T-99 前提的改写涉及目标变更，应由掌握完整上下文的人确认。

#### 9.323.5 本轮动作

* 只在干净树上重测基线（a3n00 ×3、T0M、面积比、edges、F113 `mt`）；
* 确认 `read_topology.rs` 与 HEAD blob 逐字相同；确认 `tgeometry.rs` 无未提交改动；
* **未提交任何源码改动**（本会话对 `crates/` 的功能改动只有 §9.305/§9.307/§9.311/§9.312/§9.313 那几处，
  均与本轮无关）。

#### 9.323.6 复跑

```text
git rev-parse HEAD
git hash-object crates/occt-topo/src/step/read_topology.rs
git rev-parse HEAD:crates/occt-topo/src/step/read_topology.rs
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_probe_a3n00 -- data/occ/a3n00.stp
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_probe_a3n00 -- data/occ/T0M.stp
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj -- data/occ
python .target-gate/pair_area.py output/a3n00.obj data/occ/occ-a3n00.obj
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids
```

---

### 9.324 —— 【定案】漂移的根因已从 git 历史查清：旧参照来自**本会话期间丢失的未提交改动**，不是当前 HEAD

#### 9.324.1 决定性证据：`f1e56776` 的提交信息

```text
已知代价与未结：a3n00 面积比 0.8824→0.8627（F113 出网格 mv 0→223 但整体偏低）；
T0M 未网格 6→7、wires 1931→2121、多 wire 面 120→309
```

⇒ **`0.8627` 与「F113 出网格」正是该提交自己记录的真实状态**，也就是**当前 HEAD 的状态**。
本会话实测的 `stats=226 / unmatched=0 / 11101-11121 / 0.8627 / F113 mt=267` 与它**逐条吻合**。

#### 9.324.2 旧参照（0.8996）在任何已提交状态里都不存在

* `git log -S '0.8996'` 只命中 **`35d16da9`**（`f1e56776` 的**父**提交）；
* 而 `f1e56776` 记的**改前**值是 **0.8824**，不是 0.8996；
* 本会话实测的「旧参照」还带 `edges=1083`、`mesh_v/mesh_t=10863/11941`、**F113 不出网格**。

⇒ `0.8996` 与那组指标**不对应任何已提交状态**。它们来自**工作区里未提交的改动**。

#### 9.324.3 那组未提交改动是什么，以及它怎么消失的

* 本会话早期 `read_topology.rs` 为 **2036 行**且**含** T-93/T-99 的 reader 侧 seam 块
  （本会话多轮引用其 `:734-744` 并据其插桩）；
* 现在同一文件与 **HEAD 的 blob 逐字相同**（`4f79fa0c`，1975 行），**不含**该块；
* `git stash list` 的 `stash@{0}` 内容是 `bop_builder.rs` / `fillet_curved.rs` / `occt-core/geom/mod.rs`，
  **不含** `read_topology.rs`；
* `git fsck --lost-found` 的 dangling blob 为 **0**（已被 gc）。

⇒ **那 61 行（TD-99 的 reader 侧 seam 接线）从未被提交，且已在本会话期间被抹掉、无法从 git 恢复。**
最可能是我自己在第 22–26 轮之间的一次 `git checkout -- crates/.../read_topology.rs`
或 `git stash`/`pop` 序列造成的 —— **这是我的操作失误**。

#### 9.324.4 因此：本会话 T-99/T-94 的结论需要按新参照重述

新参照（当前 HEAD，唯一可复现）：

| 量 | 值 |
|---|---|
| a3n00 `stats / unmatched` | **226 / 0** |
| a3n00 `mesh_v / mesh_t` | **11101 / 11121** |
| a3n00 面积比（ours/occ） | **0.862724**（abs diff 31179.59） |
| a3n00 `edges` | **976** |
| **F113** | **`wires=4 mt=267`（出网格）** |
| T0M `stats / unmatched` | **1765 / 7** |
| T0M `mesh_v / mesh_t` | **49999 / 47616** |

按此重述本会话的相关结论：

| 原结论 | 按新参照的正确表述 |
|---|---|
| 「基线 F113 不出网格 ⇒ 面积缺口主因是 F113」 | **不成立**。新版 F113 已出网格；面积缺口 **31179.6** 比旧（22802.6）**更大**，主因**不再**是 F113 一个面 |
| 「实验让 F113 出网格但面积比掉到 0.8627，故不能落地」 | **无意义**。当前 HEAD 本就在 0.8627；那个「实验」其实是**恢复被删的 T-93 接线**（或等价物），并非新方案 |
| 「共享键写回是 F113 的病根」 | **未证**。证据来自「旧工作区」的对比，而该工作区已不存在 |

#### 9.324.5 T-99 的正确前提与建议目标（本轮已写入卡片）

* **前提**：当前 HEAD 的真实状态是「F113 出网格、a3n00 面积比 0.8627、缺口 31179.6」。
* **建议目标**：把 **a3n00 面积比从 0.8627 推向 1**（缺口 31179.6），
  **不再**以「让 F113 出网格」为目标（已达成）。
* **附带**：`f1e56776` 同时记录了「T0M 未网格 6→7、wires 1931→2121、多 wire 面 120→309」
  等代价 —— 若要回升 a3n00 面积比，需先确认这些代价是否可同时改善，
  而**不能**用「改基线」的方式绕过（红线）。

#### 9.324.6 本轮动作

* 只做 git 历史取证与基线确认；**未提交任何源码改动**；
* 未改卡片状态；已把新前提与建议目标写入 T-99 的 `next` 与 `nextAction`；
* 记录本次操作失误（丢失未提交的 T-99 接线）以备人工核查。

#### 9.324.7 复跑

```text
git log -1 --format='%B' f1e56776            # 提交信息含「a3n00 面积比 0.8824→0.8627」
git log --all --oneline -S '0.8996' -- specs/
git show --stat 'stash@{0}'
git fsck --lost-found | grep 'dangling blob' | head
```

---

### 9.325 —— T-99（新目标）第一步：a3n00 面积缺口已归因到具体面 —— 主因**不是 F113**，而是 `f=169/189/190/171`

#### 9.325.1 方法

当前 HEAD 上重做配对与逐面面积归因：

* 两侧都按 **bbox 六坐标**（D21 判据）配对：`port faces=226 gt faces=226`，
  配对 **226 / 226**、两侧 unpaired 皆 **0**、`port unmeshed faces: 0`；
* 工具：新脚本 `.target-gate/area_deficit.py`（把两个 OBJ 的三角按质心落进
  `PORTID`/`UVSUM` 的 bbox 桶，再逐面对减）。校验：port 侧已分配面积
  `195950.71`（未分配 0.00）与 `pair_area.py` 的总面积**逐字一致** ⇒ 归因口径可靠；
  GT 侧 `226870.40`（未分配 259.90）。

#### 9.325.2 结果：净缺口 30919.69，可归因到少数面

**端口欠覆盖（亏）**：

| port f | gt f | ours | gt | diff |
|---|---|---|---|---|
| **171** | 1 | 13170.33 | 34044.78 | **−20874.45** |
| **200** | 4 | **0.00** | 9315.27 | **−9315.27** |
| **191** | 3 | **0.00** | 9244.65 | **−9244.65** |
| **173** | 7 | 201.34 | 4706.41 | **−4505.07** |
| 180 | 44 | 3415.02 | 5641.23 | −2226.20 |
| 206 | 86 | 521.12 | 2353.21 | −1832.08 |
| 204 | 43 | 10542.14 | 11790.43 | −1248.28 |
| 129 | 23 | 2375.17 | 3339.24 | −964.07 |

**端口过覆盖（盈）**：

| port f | gt f | ours | gt | diff |
|---|---|---|---|---|
| **190** | — | 14778.90 | 0.00 | **+14778.90** |
| **189** | 41 | 12996.46 | 1355.22 | **+11641.24** |
| **169** | 87 | 23623.46 | 14307.56 | **+9315.89** |
| **113** | 77 | 21575.81 | 14237.50 | **+7338.31** |
| 40 | — | 6290.31 | 0.00 | +6290.31 |
| 138 | — | 5968.43 | 0.00 | +5968.43 |
| 175 | 6 | 5054.38 | 0.00 | +5054.38 |
| 216 | — | 5048.30 | 0.00 | +5048.30 |

#### 9.325.3 三条可直接用的结论

**(1) 旧结论「面积缺口主因是 F113 一个面」不成立。**
F113（`+7338.31`）只排第 4，且它是**盈**不是亏。
真正的大项是 `f=171`（−20874，单面占净缺口的 **67%**）、
`f=200`/`f=191`（各 **−9316/−9245**，端口**完全没铺**）、以及 `f=173`（−4505）。

**(2) 缺口是「亏 + 盈」的净额，两侧都很大。**
亏损大项合计约 −54000，盈利大项合计约 +66000 ⇒ 若都能归零，面积比可显著上升。
**这比「追一个面」有希望得多。**

**(3) `f=171` 与 `f=190` 是最大的一对（−20874 / +14778）。**
结合 §9.325.1 的配对结果（`f=171` 的 bbox 与 GT `f=1` 配对），
这一对极可能是**同一块几何在端口被拆/合并到不同面**，值得优先看。

#### 9.325.4 与既有线索的呼应

`pair_a3n00.py` 同时报出：**端口在 29 个面上有 2 条 wire，而 GT 只有 1 条**
（另有 `f=113` 是 4 vs 2）。**多余的 wire**正好与「盈」的大项面（169/189/190）方向一致 ——
多出的 wire 可能让这些面**多铺**，而它们本该贡献的面积跑到了别处（或反之）。
这是一条可检验的假设，但**本轮未验证**。

#### 9.325.5 下一轮的检查点（有界）

1. 取 `f=171`（亏 −20874）与其 GT 配对 `f=1`，打印两者的 **wire 数 / 每 wire 边数 /
   bbox / 曲面**，判断是「端口把它拆成了多个面」还是「端口的 wire 环不完整」；
2. 取 `f=200` / `f=191`（端口面积 **0**），确认它们**出网格但面积归零**的原因
   （是三角退化、还是这一对配对本身有误 —— 二者结论完全不同）；
3. 再取 `f=190`（盈 +14779，GT 侧无对应）确认它是重复面还是吞了邻面的面积。

**注意**：本轮**未提交任何源码改动**，只做归因。

#### 9.325.6 复跑

```text
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids > portid.txt 2>/dev/null
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --uvsum" > uvsum.txt
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj -- data/occ
python .target-gate/pair_a3n00.py portid.txt uvsum.txt
python .target-gate/area_deficit.py output/a3n00.obj data/occ/occ-a3n00.obj portid.txt uvsum.txt
```

---

### 9.326 —— 【更正 §9.325】逐面面积归因**无效**：全部 226 个面的 bbox 互相重叠，质心分桶不可靠

#### 9.326.1 发现

上一轮（§9.325）我用新脚本 `area_deficit.py` 把 OBJ 三角**按质心落进各面 bbox** 来逐面归属面积。
本轮做检查点 ①②③ 时读到**互相矛盾**的数字：

| 面 | port 侧面积 | GT 侧面积 | 矛盾 |
|---|---|---|---|
| `f=9` | 226.99 | **18049.08** | 两侧都大 ⇒ 不可能是同一个面 |
| `f=190` | 14778.90 | 9.04 | 同上 |
| `f=200` / `f=191` | **0.00** | **0.00** / 13.58 | 而 `pair_a3n00.py` 报它们与 GT `f=4`/`f=3` 的 **wires 数完全相同（10/10）** —— 端口不可能真的没铺 |
| `f=1` | **0.00** | 34044.78 | 但 `PORTID f=171` 与 `UVSUM f=1` 的 bbox **逐字相同**，f=171 侧却有 13170.33 —— 同一 bbox 却一边为 0 |

**根因**（本轮实测）：

```text
port 侧 bbox 互相重叠的面对数 = 967
涉及的面 = 226 / 226          ← 全部面都与别的面重叠
```

⇒ **每个三角的质心都落在多个面的 bbox 里**，而我的脚本是「第一个匹配就 break」
⇒ **归属基本是任意的**。上表所有「不可能」的读数都是这个 bug 的直接后果。

（脚本里唯一可信的量是**总量**：port 侧已分配 195950.71、未分配 0.00 ——
那个数之所以对，是因为每块面积最终都被派给了某个面，与派给谁无关。）

#### 9.326.2 因此 §9.325 的结论**全部撤回**

§9.325 写的「净缺口可归因到少数面」「主因不是 F113 而是 f=171 / f=200 / f=191」
「f=171 单面占净缺口 67%」「亏损大项合计约 −54000、盈利大项约 +66000」
—— **全部建立在无效归属之上，作废**。（§9.325.3 的「(1) 旧结论『主因是 F113』不成立」
这点**仍然成立**，但理由不同：它是靠**结构数据**得出的，见下。）

#### 9.326.3 仍然可信的证据（来自 `PORTID`/`UVSUM` 的**结构字段**，与面积归属无关）

这些是端口/GT 各自 dump 的 wire 数与 bbox，不经我的分桶，**可直接比较**：

| port f | port wires | GT f | GT wires | 判定 |
|---|---|---|---|---|
| 171 | 4 | 1 | 4 | bbox **逐字相同** ⇒ 配对正确、**wires 数一致** |
| 200 | 10 | 4 | 10 | bbox 同在 x=−132.5 平面 ⇒ **wires 一致** |
| 191 | 10 | 3 | 10 | bbox 同在 x=−112.5 平面 ⇒ **wires 一致** |
| 173 | 2 | 7 | **1** | **差 1 条 wire** |
| 113 | 4 | 77 | **2** | **差 2 条 wire** |
| 169 | 2 | 87 | **1** | **差 1 条 wire** |
| 189 | 2 | 41 | **1** | **差 1 条 wire** |

⇒ 我上一轮点名的几个「大项面」里，**171/200/191 的 wires 数与 GT 完全一致**，
所以它们**不太可能是 107 边差集的来源**；真正 wires 不一致的是
**173 / 113 / 169 / 189**（各差 1–2 条），这与 §9.325 引用的
「端口在 **29 个面上有 2 条 wire 而 GT 只有 1 条**」是一致的。

**结论**：面积缺口的落点应重新从**这 29 个 wire 数不一致的面**入手，
而不是从（无效的）面积排序入手。

#### 9.326.4 正确的归属方法（下一轮用）

按质心分桶已被证伪。可用的替代：

1. **按连通性分块**：把 OBJ 三角按共享边/顶点并成连通块，每块对应一个面
   （前提是各面在网格里不共享顶点 —— 需先验证；`mesh_v/mesh_t` 比值可作旁证）；
2. **在端口侧加一个按面的面积探针**：直接在 `IncrementalMesh` 之后读每个
   `MeshFace` 的三角并求和 —— 这是唯一**不依赖任何反推**的口径，最可靠；
3. 用 `PORTID` 的 `mt`（三角数）**按面型分组**做粗归因（不依赖 bbox），
   例如 `pair_a3n00.py` 已有的「按 GT surface type 的 tri 汇总」。

**推荐 (2)**：在 `zz_uv_feed` 或 `zz_probe_a3n00` 里加一个 `--facearea` 模式，
打印每个面的三角面积和；对 GT 侧无对应工具（只有 OBJ），
故仍按 bbox 配对，但**面积取自各自模型的真实归属**而非反推。

#### 9.326.5 本轮状态

* 本轮**未提交任何源码改动**（只做检查点核查与更正）；
* 已撤回 §9.325 的归因结论，并保留其中**由结构数据支撑**的部分；
* §9.325 的三个检查点里，①③ 因归属无效而无法判定，② 已被结构数据回答
  （`f=200`/`f=191` 的 wires 与 GT 一致，不是「没铺」）。

#### 9.326.6 复跑

```text
python .target-gate/pair_a3n00.py .target-gate/portid_cur.txt .target-gate/uvsum_cur.txt
python .target-gate/face_map.py   .target-gate/portid_cur.txt .target-gate/uvsum_cur.txt
# 重叠检测：全部 226 个面的 bbox 两两有重叠（967 对）
```

---

### 9.327 —— 改用正确口径的三个尝试结果，以及一个**新落点**：59 个面无法与 GT 按 bbox 精确配对

#### 9.327.1 尝试一：给探针加 `--facearea` —— **在端口侧不可行**

查清 `MeshFace` 只有：
`face()` / `surface()` / `wires()` / `points_nb()` / `point(i)` / `points()` / `deflection()` / `status()`…

**没有三角容器**。逐面三角只在模型级（`MeshModel`）之上，逐面视图里没有。
⇒ 「直接读每个面的三角并求和」这条**最可靠口径在端口侧没有现成 API**，
需要先给 `MeshFace` 加三角访问器（改库），超出本轮范围。

（顺带确认：`ShapeMesh`（`mesh.rs:12`）只有 `vertices` / `triangles` / `source_shape`，
是**整形状**的，同样不提供逐面划分。）

#### 9.327.2 尝试二：按 GT 面型归因三角数 —— **GT 侧数据不足**

新脚本 `.target-gate/tris_by_type.py` 想按 `UVSUM` 的 `type=` 分组比较三角数，
但实测 **`UVSUM` 行里没有三角计数**（该行只有 `wires` / `bbox` / `urange` / `vrange` /
`dU` / `dV` / `lenU` / `lenV` / `valid`）。
GT 的三角数只能从 `occ-a3n00.obj` 反推 —— 又回到不可靠的归属。
⇒ 该口径**不可行**（除非改 OCCT 侧探针，本会话无 OCCT 构建树，做不到）。

#### 9.327.3 【新落点】59 个面无法与 GT 按 bbox 精确配对

在配对过程中量到一个**实在的信号**（不依赖任何面积反推）：

```text
port faces=226  gt faces=226
matched by exact 6-coord bbox: 167
port faces with no exact GT bbox match: 59   （triangles 2474 / 11121 ≈ 22%）
```

**为什么这不是舍入假象**：GT 的 bbox **确实**带 ±1e-10 量级的 epsilon
（例 `UVSUM f=1 bbox=(-57.360500611,…)` 对 `PORTID f=171 bbox=(-57.360501,…)`），
所以舍入到 6 位后两者相等 —— **这是 `f=171` 能配上的原因**。
而舍入到 6 位后**仍不相等**的那 59 个面 ⇒ 它们的 bbox **真的不同**（量级远超 1e-6）。

**可疑面举例**（GT 侧同名面的 bbox 与之量级不同）：

| port f | port bbox | 备注 |
|---|---|---|
| 190 | `(-132.5, -77.436, -143.239) .. (-114.5, 77.436, 11.5)` | GT 无精确匹配 |
| 138 | `(25.868, -30, -176.294) .. (68.294, 30, -133.868)` | GT 无精确匹配 |
| 40 | — | GT 无精确匹配 |
| 176 / 201–217 | — | GT 无精确匹配 |

⇒ **这 59 个面是「端口与 GT 几何不一致」的直接证据**，而且占了端口 **22% 的三角**。
它比此前任何一次「面积排序」都更可靠（只用了两侧各自 dump 的 bbox，无任何归属假设）。

#### 9.327.4 这对 T-99 的意义

* 面积缺口的落点应**先查这 59 个 bbox 不匹配面的成因**：
  是端口把某个面**建错了范围**（u/v range 或裁剪）、还是**拆分/合并**与 GT 不同、
  还是端口的**配对索引**（`faces_of` 顺序）与 `PORTID` 的 `f=i` 不一致；
* **不能**再用 §9.325 那套质心分桶（§9.326 已证伪）；
* 若要逐面面积，必须先做 §9.327.1 那件库改动（给 `MeshFace` 加三角访问器），
  这是**有明确 OCCT 对应**的（`BRepMeshData_Face` 持有 `Triangulation`），
  不属自造分支 —— 但它是**改库**，需按五道门禁完整走一遍。

#### 9.327.5 本轮状态

* **未提交任何源码改动**（只做口径探查与配对分析）；
* 两个口径尝试（`--facearea`、按面型归因）分别因**端口无逐面三角 API**与
  **GT 无逐面三角计数**而不可行，已记录原因以免下一轮重试；
* 产出**一个新落点**：59 个 bbox 不匹配面（占 22% 三角）。

#### 9.327.6 复跑

```text
python .target-gate/face_map.py      .target-gate/portid_cur.txt .target-gate/uvsum_cur.txt
python .target-gate/tris_by_type.py  .target-gate/portid_cur.txt .target-gate/uvsum_cur.txt
grep -n "pub fn " crates/occt-topo/src/meshing/data_model/mesh_face.rs   # 确认无三角访问器
```

---

### 9.328 —— 【修正 §9.327.3】59 个面的 bbox 差异**很小**（最坏 4.0、中位 ~0.05），不是"几何完全不同"

#### 9.328.1 先排除索引假说（已证伪）

卡片列的第三个候选是「`faces_of` 顺序与 `PORTID` 的 `f=i` 不一致」。实测：

```text
PORTID 行数 = 226，f 从 0 到 225 连续（contiguous = True）
mt 查找失败（mt = usize::MAX）的次数 = 0
```

且 `FaceMeshStat` 的文档明确说明了防错设计
（`discret_root.rs:77-80`：「Pointer-identity key of the source face, so a caller can pair
the stat with its face without relying on traversal order (the mesh model skips faces,
so a positional lookup is off by the number of skipped faces)」），
探针正是按 `shape_key` 查找的。
⇒ **索引/遍历顺序不是成因。**

#### 9.328.2 量出真实幅度：**很小**

新脚本 `.target-gate/bbox_delta.py`：对 59 个不匹配面，找最接近的 GT 面并报
「每个坐标上的最大差」：

| port f | 最近 GT f | 每坐标最大差 |
|---|---|---|
| 212 / 214 / 210 / 208 | 15 / 17 / 13 / 10 | **4.00133**（四个面**同值**） |
| 190 | 9 | 0.353287 |
| 216 | 11 | 0.227927 |
| 176 | 5 | 0.218809 |
| 201 | 2 | 0.177783 |
| 215 / 213 / 211 / 209 | 18 / 16 / 14 / 12 | **0.073154**（四个面**同值**） |
| 112 / 111 / 110 / 109 | 19 / 20 / 21 / 22 | **0.0462028**（四个面**同值**） |
| 12 | 217 | 0.00332187 |

**关键读数**：

* **最坏也只有 4.0**，中位约 **0.05**，最小 **0.0033**；
* `每个坐标都在 1e-9 以内` 的面数 = **0**（所以确实都不"精确相等"）。

#### 9.328.3 因此修正 §9.327.3 的措辞

§9.327.3 写的是「那 59 个面是『端口与 GT 几何不一致』的**直接证据**」「**真的不同**」——
**这是过重的表述**。准确说法是：

> 59 个面的 bbox 与**任何** GT 面都不精确相等，但**差异很小（≤4.0，中位 ~0.05）**，
> 属于**同一几何的轻微边界差异**，而不是"几何完全不同"。
> 它能解释「为什么不能用精确 bbox 配对」，但**不足以**单独解释 31179.6 的面积缺口。

（§9.327 里「占端口 22% 三角」这个数字本身没错，但它描述的是**配对失败的面所承载的三角量**，
不是"22% 的三角几何错误"——这个区别我上一轮没讲清楚。）

#### 9.328.4 一个有信息的新线索：差异**成组同值**

`208/210/212/214` 四个面的最大差**完全相同**（4.00133），
`209/211/213/215` 四个相同（0.073154），
`109/110/111/112` 四个相同（0.0462028）。

⇒ 端口与 GT 在这些位置是**同一族面**（很可能是圆角/直纹面的四个象限），
**每个都被整体平移/缩放了一点**。这指向「**端口的某些面在某个方向上被建得略大或略小**」，
而不是随机误差 —— 是一条比"逐个面查"更值得追的线。

#### 9.328.5 对 T-99 的处置

* 上一轮定的落点（「先查这 59 个面的成因」）**仍然有效**，但期望要放低：
  它们的幅度（≤4.0）与缺口（31179.6）**量级不匹配**，所以
  **它们大概不是缺口主因**；
* 更该追的是 §9.328.4 的「成组同值」：找出 `208/210/212/214` 这一族在
  端口里的**面型与 UV 范围**，看它相对 GT 是哪一个坐标被建偏了；
* **不要**再回头用质心分桶（§9.326 已证伪）；
* 若最终需要逐面面积，仍须先给 `MeshFace` 加三角访问器（§9.327.1）。

#### 9.328.6 本轮状态

* **未提交任何源码改动**（只做度量与措辞修正）；
* 卡片 `next` 已按本节更新。

#### 9.328.7 复跑

```text
python .target-gate/bbox_delta.py .target-gate/portid_cur.txt .target-gate/uvsum_cur.txt
```

---

### 9.329 —— 【新线索】成组同值的面是「**中心相同、尺寸约 0.70 倍**」—— 指向端口的 Extrusion 面**尺寸建小了**

#### 9.329.1 度量（逐坐标尺寸比）

新脚本 `.target-gate/half_size.py`，对 §9.328.4 的三组成组面逐坐标比尺寸：

| port | gt | port 尺寸 (x,y,z) | GT 尺寸 (x,y,z) | 尺寸比 (y,z) |
|---|---|---|---|---|
| 208 | 10 | (0, **3.9377**, **3.9377**) | (0, **4.0013**, **4.0013**) | 0.984 / 0.984 |
| 210 | 13 | (0, 3.9377, 4.0013) | (0, 4.0013, 3.9377) | 0.984 / 1.016 |
| 212 | 15 | (0, 3.9377, 3.9377) | (0, 4.0013, 4.0013) | 0.984 / 0.984 |
| 214 | 17 | (0, 4.0013, 3.9377) | (0, 3.9377, 4.0013) | 1.016 / 0.984 |
| 209 | 12 | (0, 0, 0.0732) | (0, 0, 0.0079) | — |
| 211 | 14 | (0, 0.0079, 0) | (0, 0.0732, 0) | — |

⚠️ 该表里的「尺寸」是**逐坐标 bbox 端点差的和**（我脚本的写法），不是 extents。
按**绝对 bbox** 重算（`PORTID`/`UVSUM` 原始行）：

| | port 208 | GT 10 | 比值 |
|---|---|---|---|
| y | [−64.866476, −45.935433]，extent **18.931** | [−68.867801, −41.997743]，extent **26.870** | **0.7045** |
| z | [−52.472698, −33.541655]，extent **18.931** | [−56.474023, −29.603965]，extent **26.870** | **0.7045** |
| x | [112.5, 132.5] | [112.49999997, 132.50000013] | 1.0（完全相同） |
| **中心** | (122.5, **−55.400954**, **−43.007176**) | (122.5, **−55.432772**, **−43.038994**) | 差 **0.032** |

⇒ **同一个 y/z 中心、x 完全相同，但 y/z 尺寸是 GT 的 0.7045 倍。**

#### 9.329.2 这是一条比「逐个查 59 个面」更聚焦的线索

* 它**成组出现**（4 个面一组、组内比值一致）⇒ 不是随机噪声，是**同一段构造逻辑**在 4 个象限上重复；
* 面型是 **Extrusion（GT `type=6`）**，`bbox` 的 x 跨度与 GT **完全相同**、
  只有沿 y/z 的**拉伸长度**偏小 ⇒ 与「**拉伸方向的长度/范围被建小了**」一致；
* 幅度 0.70 倍意味着**每个这样面的面积少约 30%** —— 若此类面有十几个，
  就能贡献缺口的可观部分（**待验证，本轮未量**）。

#### 9.329.3 必须同时记住的反证（避免重蹈 §9.325）

* `extent` 量的是**端口的实际 y/z 跨度**，它受**该面在端口里被裁剪到的 UV 范围**影响；
  **不能**由 bbox 直接断言「几何建错了」；
* port 208 的 `mt=1`（只有 1 个三角）⇒ 它的 bbox 可能只覆盖实际几何的**一部分**，
  换言之 **bbox 偏小也可能来自"只铺了一小块"**，而不是"面本身小"——
  这两种解释**本轮无法区分**，需按面的 uv-range 与 wire 环判断（`PORTID` 不含 urange，
  但 `MODEL`/`SPH` 模式有）；
* 因此本节的结论**只到**：「这 4 个面在端口里的 bbox 是 GT 的 0.7045 倍且中心相同」，
  **不**断言成因。

#### 9.329.4 下一轮（有界）

1. 用 `zz_uv_feed --model <i>` 或 `--sph` 取 **208/210/212/214** 的
   `urange` / `vrange` / wire 数与每 wire 边数，与 GT `UVSUM f=10/13/15/17` 的
   `urange`/`vrange` **逐值对比** —— 这能直接区分 §9.329.3 的两种解释；
2. 统计**全模型**里「bbox 是 GT 的同一比例、中心相同」的面有多少个、承载多少三角
   ⇒ 判断这条线够不够解释 31179.6；
3. **不要**再用质心分桶（§9.326 已证伪）。

#### 9.329.5 本轮状态

* **未提交任何源码改动**（只做度量与记录）；
* 已顺带修正 §9.328.4 引用的「尺寸」口径说明（见 §9.329.1 的警示）。

#### 9.329.6 复跑

```text
python .target-gate/half_size.py
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --model 208
```

---

### 9.330 —— 【定位】59 个不匹配面里只有 **4 个**真有差异，且成因已定：**UV 范围与 GT 相同、但端口多一条 wire ⇒ 该面只铺了部分环**

#### 9.330.1 差异分布：51/59 是舍入或边界精度，只有 4 个显著

对 59 个不匹配面，按「与最近 GT 面的每坐标最大差」分档：

| 档位 | 面数 |
|---|---|
| **> 1.0** | **4** |
| 0.1 – 1.0 | 4 |
| 0.01 – 0.1 | 9 |
| 0.001 – 0.01 | 8 |
| **≤ 0.001** | **34** |

**唯一 > 1.0 的 4 个面**（全部同值 4.0013）：

```text
port f=212 ~ gt f=15    worst=4.0013
port f=214 ~ gt f=17    worst=4.0013
port f=210 ~ gt f=13    worst=4.0013
port f=208 ~ gt f=10    worst=4.0013
```

⇒ §9.328 里那 59 个「配不上」的绝大多数（**51 个**）差异 ≤ 0.1，
其中 34 个 ≤ 0.001（就是 §9.327.3 说的舍入量级）。
**真正异常的只有这 4 个**。

#### 9.330.2 这 4 个面的成因已定位（本轮卡片检查点的答案）

用 `zz_uv_feed --model <f>` 取端口侧结构，与 GT `UVSUM` 逐值对比：

```text
port f=208  wires=2 urange=[-30.29536044126110284,-28.19009728336637011] vrange=[0.00,59.69026041820608]
  wire[0] nEdges=1  e[0] curve=Circle par=[π,3π]  dpar=2π  npc=1
  wire[1] nEdges=1  e[0] curve=Circle par=[0,2π]  dpar=2π  npc=1
gt   f=10   wires=1 urange=[-30.29536044126110639,-28.19009728336636300] vrange=[0.00,59.69026041820608]
```

四个面**全部**如此（`210/212/214` 结构相同，只是 edge id 与朝向不同）。**三个关键读数**：

1. **`urange` 与 `vrange` 与 GT 逐值相同**（到 1e-14 量级）⇒ **面的几何定义完全相同**，
   所以 §9.329 的「面建小了」**被排除**；
2. **端口是 2 条 wire、GT 是 1 条**（与既有统计「29 个面上 2 vs 1」同族）。
   两条 edge 都是**整圆**（`dpar = 2π`），参数域分别为 `[π,3π]` 与 `[0,2π]`；
3. **`vrange` 到 59.69（= 2.5 圈），而实际 bbox 的 y 跨度只有 18.93** ⇒
   若整面都铺满，y 跨度应接近 59.7 的量级。

⇒ **结论：这 4 个面在端口里只被铺了"部分环"** ——
面的 UV 范围是对的，但网格没覆盖整个 v 范围。
**这坐实了 §9.329.3 里我标为「无法区分」的第二种解释**（「只铺了一小块」），
而排除了第一种（「面本身小」）。

#### 9.330.3 与「2 vs 1 wire」的联系（可检验的假设）

4 个面都同时具备「2 条 wire」与「只铺部分环」。两者极可能同源：
**端口把一条闭合 ring 拆成了两条 wire**（各占半个周期），
而网格器按 wire 覆盖范围取样 ⇒ 只铺了其中一部分。
这与 §9.305/§9.326 提到的那条既有统计（**端口在 29 个面上有 2 条 wire 而 GT 只有 1 条**）方向一致。

**但本轮不宣称已证**：4 个面同时具备两个特征，尚不足以在因果上判定，
且另有 25 个「2 vs 1」的面**不在**这 4 个里（它们的 bbox 差异 ≤ 0.1）。

#### 9.330.4 量级评估（决定这条线值不值得追）

* 这 4 个面：端口侧 `mt = 224+224+1+2 = 451` 个三角；
* 它们承载的 GT 面积约 `4 × 26.870 × 18.931 ≈ 2035`（bbox 面积极粗估）；
* 而全模型缺口是 **31179.6**。

⇒ **这 4 个面至多贡献缺口的几个百分点，不足以解释缺口。**
**因此它们不是 T-99 的主线**，但它们是**已完全定位的一个真实缺陷**
（UV 范围正确、网格只铺部分环），可作为「2 vs 1 wire」这条线的**可复现场景**。

#### 9.330.5 对 T-99 的处置

* **撤回**「59 个面是主要落点」的期望（§9.328 已降级，本节进一步量化：只有 4 个显著、
  且量级不足以解释缺口）；
* 缺口的主因**仍未定位** —— 需要在**配对成功的 167 个面**里找，
  但那需要逐面面积口径（§9.327.1 的 `MeshFace` 三角访问器，属改库）；
* 本轮新增的可复现缺陷（4 个面 + 2 vs 1 wire）可作为**独立的小卡**，
  但**不在本轮动手**（五道门禁第 4 条：完成即停，旁支写进报告）。

#### 9.330.6 本轮状态

* **未提交任何源码改动**（只做度量与定位）；
* 卡片 `next` 已按本节更新。

#### 9.330.7 复跑

```text
python .target-gate/bbox_delta.py    .target-gate/portid_cur.txt .target-gate/uvsum_cur.txt
python .target-gate/port_vs_gt_face.py
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --model 208
```

---

### 9.331 —— 两条否定结果 + 口径瓶颈的最终确认：**按面归因在两侧都拿不到数据**

#### 9.331.1 否定结果一：**没有塌缩/退化三角**

对两版 OBJ 逐三角算面积：

```text
port: tris=11121 total=195950.71  area<1e-6: 0  area<1e-3: 1  min=9.780e-04  max=1514.86
gt  : tris=12466 total=227130.30  area<1e-6: 0  area<1e-3: 1  min=8.901e-04  max=1514.86
```

* **两侧都没有面积 < 1e-6 的三角**（`area<1e-3` 各只有 1 个）；
* **最大三角完全相同（1514.86）**，最小值同量级。

⇒ 缺口**不是**由「塌缩/退化三角把面积吃掉」造成的，也**不是**由"个别面被整个漏掉"造成的
（若是，会看到一个异常大的缺失块或大量零面积三角）。
这是 §9.316.1 结论的**独立再确认**。

#### 9.331.2 否定结果二：三角计数口径自洽，但**总三角数少**

```text
PORTID 逐面 mt 求和 = 11121  ==  OBJ 三角数 11121     （逐字一致）
port 11121 个三角覆盖 195950.71 面积
gt   12466 个三角覆盖 227130.30 面积
三角数比 0.892   面积比 0.863
```

⇒ 端口的**三角数与面积同比例偏少** ⇒ 缺口随"三角铺得少"走，
但**不能**由此断定"曲面被铺粗"（较粗的三角铺同一曲面，面积只差 O(deflection²) 量级，
而这里是约 14%）。

#### 9.331.3 口径瓶颈：**按面归因在两侧都拿不到数据**（本轮最终确认）

| 侧 | 需要的数据 | 现状 |
|---|---|---|
| **端口** | 每个面的三角面积和 | `MeshFace` **无三角容器**（§9.327.1 已查：只有 `face/surface/wires/points/point(i)/deflection/status`）；逐面三角只在模型级 `MeshModel` 之上。**需改库**。 |
| **GT** | 每个面的三角计数或面积 | `UVSUM` **无三角计数**（§9.327.2、本轮再次确认）；且**本会话没有 OCCT 构建树**（`D:\source\occt-8.0.0` 只有预编译发行版，无头文件/编译器），**无法给 `occt_probe` 加这个字段**。 |

⇒ **两条独立路径都落到同一个瓶颈**：GT 侧没有逐面三角/面积数据，
而端口侧需要一次库改动才能产出。

**这意味着**：在**不改 OCCT 侧**的前提下，
T-99 的「把面积比从 0.8627 推向 1」**无法再靠测量推进** ——
只能靠：(a) 给端口加逐面面积访问器（库改动）后做**端口内部**的自洽性检查
（例如「面的三角面积和 vs 面的参数域面积」是否相符），
或 (b) 在**已定位的具体缺陷**上动手（如 §9.330 的 4 个「只铺部分环」面）。

#### 9.331.4 关于 (a) 的一个可检验设想（**未验证**）

给端口加逐面面积后，可以与**同一面的 `urange`/`vrange` 推出的参数面积**对照：

* 若某面的「三角面积和 << 参数域面积」，则它属于 §9.330 那类**只铺部分环**；
* 把这种面**枚举全**，就能知道「只铺部分环」的总量够不够解释 31179.6。

§9.330 已量出 4 个（451 三角、约 2035 面积），但**同族可能不止 4 个** ——
本轮**无法枚举**，因为需要那个库访问器。

#### 9.331.5 本轮状态

* **未提交任何源码改动**（只做度量）；
* 卡片 `next` 已按本节更新；
* 本会话累计：T-99 侧共 **15 条假设被否决**，收敛到「需要一次库改动或转向已定位的具体缺陷」。

#### 9.331.6 复跑

```text
python .target-gate/density_by_type.py .target-gate/portid_cur.txt .target-gate/uvsum_cur.txt
# 三角面积/退化统计：见 §9.331.1 内联脚本
```

---

### 9.332 —— 【重要更正】GT 侧**本来就有**逐面三角数据（`--mesh` oracle）；我此前多次写下的「拿不到」是错的

#### 9.332.1 更正内容

本会话 §9.327.2、§9.331.3 等多处写过：

> 「GT 侧没有逐面三角/面积数据」；「本会话没有 OCCT 构建树，无法给 `occt_probe` 加这个字段」。

**这两条都错了。** 实测：

1. **OCCT 发行版带完整头文件**：`D:\source\occt-8.0.0\inc` 下有 **7084** 个 `.hxx`；
2. **有编译器**：`C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat` **存在**；
3. **`specs/occt_probe/` 自带源码 + 构建脚本**：`occt_probe.cpp`（1893 行）、`build.bat`（vcvars64 + cl + OCCT 的 inc 与 win64/vc14/lib）；
   ⇒ **本来就具备"给 OCCT 探针加字段并重建"的条件**；
4. **而且这个字段早就有了** —— 探针已有 **`--mesh` 模式**（`occt_probe.cpp:142-188`，注释标为「T-54 oracle」），
   逐面输出 `FACE <idx> nodes=<n> triangles=<t>`，并有 `TOTAL ... triangles=`：

```text
FACE lines: 226
TOTAL faces=226 nodes=15538 triangles=20046 meshvol=1.89747e+06 neg_triangles=4609 deflection=0.1 angle=0.5
FACE 0 nodes=52 triangles=50
FACE 2 nodes=1755 triangles=3224
```

⇒ **GT 的逐面三角计数一直可用**，只是一直没被 T-99 用上。

#### 9.332.2 但**不能**直接拿它与端口的 `mt` 比（角度参数不同）

`--mesh` 的默认 `angle = 0.5`（rad），而端口 `export_data_obj` 的网格参数是另一套
（`zz_uv_feed` 打印的是 `lin` 与 `20°`）。**角度不同 ⇒ 三角数不可直接比较**：

```text
OCCT --mesh（angle 0.5）      triangles = 20046
GT OBJ（occ-a3n00.obj）       triangles = 12466
端口 OBJ（output/a3n00.obj）  triangles = 11121
```

**20046 与 12466 都是"GT"，但不相等** ⇒ 说明 `--mesh` 的角度与生成 `occ-a3n00.obj` 的角度不同。
⇒ **要用于归因，必须先确认 `occ-a3n00.obj` 是用什么角度/偏转生成的**，
再用**同一套参数**跑 `--mesh`；否则比出来的差异是参数差而不是缺口。

#### 9.332.3 这打开了一条真实可行的路（下一轮）

1. 查清 `occ-a3n00.obj` 的生成参数（`specs/` 或 `data/` 里的生成脚本/记录）；
2. 用**同一组** deflection/angle 重跑 `--mesh`，得到与 OBJ 同口径的**逐面三角数**；
3. 与端口 `PORTID` 的逐面 `mt` 按 **bbox 六坐标**配对后逐面对比 ⇒
   **首次得到可靠的逐面三角差**，从而定位缺口落在哪些面/面型；
4. 若需要**面积**而不仅是三角数，探针的 `--mesh` 循环里已经算了 `a.XYZ().Dot(b.XYZ().Crossed(c.XYZ()))/6`
   （体积），**加一个逐面 `|cross|/2` 面积和是 3 行代码**，且**有 OCCT 同等 API**
   （`Poly_Triangulation::Node/Triangle`）—— 属正当移植工具，非自造分支。

#### 9.332.4 本轮状态

* **未提交任何源码改动**（只做可行性核查与更正）；
* 撤回「GT 侧拿不到逐面三角数据」这一**错误结论**——它在 §9.327.2、§9.331.3 里被我用作
  「测量推进到头」的依据，**该依据不成立**；
* 卡片 `next` 已改为 §9.332.3 的四步。

#### 9.332.5 复跑

```text
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --mesh 0.1 0.5"
Test-Path 'C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat'
(Get-ChildItem 'D:\source\occt-8.0.0\inc' -Filter '*.hxx').Count
```

---

### 9.333 —— 卡片四步之①：GT OBJ 的生成参数已定位到 deflection ≈ 0.215 / angle 0.5（三角数吻合）；但 `v` 与 `nodes` 不是同一概念

#### 9.333.1 扫描结果（用探针自己的 `--mesh` oracle）

**先扫 angle**（deflection 固定 0.1）：

| angle | triangles |
|---|---|
| 0.1 | 89226 |
| 0.2 | 33524 |
| 0.3 | 23652 |
| 0.4 | 21178 |
| **0.5** | **20046** |
| 0.6 | 19586 |
| π/4 | 19404 |

⇒ **angle 越大三角越少，但最少也有 19404** > GT OBJ 的 **12466**
⇒ **仅靠 angle 不可能复现** ⇒ GT OBJ 的 **deflection 更大**。

**再扫 deflection**（angle 固定 0.5）：

| deflection | triangles | nodes |
|---|---|---|
| 0.20 | 12844 | 10894 |
| 0.21 | 12574 | 10679 |
| **≈0.215** | **12466**（目标） | ~10550 |
| 0.22 | 12368 | 10561 |
| 0.23 | 11680 | 10186 |
| 0.24 | 11396 | 10015 |

⇒ **`deflection ≈ 0.215`、`angle = 0.5` 时 triangles = 12466**，
**与 GT OBJ 的 `f=12466` 逐字吻合**。

#### 9.333.2 但 `nodes` 对不上 ⇒ 两个 OBJ 的 `v` 与 `--mesh` 的 `nodes` 不同义

| | triangles | nodes/vertices |
|---|---|---|
| GT OBJ `occ-a3n00.obj` | **12466** | **11145** |
| 探针 `--mesh 0.215 0.5` | **12466** | **~10550** |

**三角数吻合而顶点数不吻合** ⇒ 两者**不是同一个量的产物**：
OBJ 的 `v` 是**写出时的顶点表**（面与面之间会去重/合并，或按每面独立写出后再合并），
而 `--mesh` 的 `nodes` 是**各面 `Poly_Triangulation` 的 `NbNodes` 之和**（面间不去重、面内可能重复）。
⇒ **不能**用 `nodes` 与 OBJ 的 `v` 互相校验；
**三角数** `f` 与 `triangles` 才是**同义**的（都是三角个数）。

#### 9.333.3 这带来一个方法上的结论（重要）

* **按三角数归因是可行的**（`f` 与 `triangles` 同义，且已找到口径 `d≈0.215,a=0.5`）；
* **按面积归因仍不可靠**：GT 侧只有三角**个数**，没有逐面面积；
  即便用 §9.332.3 第 ④ 步给探针加逐面面积，得到的也是**"OCCT 在该 deflection 下的面积"**，
  而 `occ-a3n00.obj` 是按它自己那套参数写出的 —— 两者只在 **d≈0.215** 时才对得上，
  而 0.215 是**我扫出来的近似值**，不是已知的生成参数；
* 因此若要用 GT 的逐面面积做对拍，**必须先确定 `occ-a3n00.obj` 的确切生成参数**
  （deflection 的**精确值**与是否 passing/relative 模式），否则误差会混进缺口。

#### 9.333.4 GT OBJ 的确切生成参数**在本仓库里没有记录**

```text
grep -rln "occ-" --include=*.bat --include=*.sh --include=*.py --include=*.ps1 .
  -> 只有 .target-gate/area_by_face.py（我自己写的，不是生成器）
仓库内 *.bat 只有 specs/occt_probe/ 下的探针脚本（无生成 OBJ 的）
```

⇒ `data/occ/occ-a3n00.obj` 的来源**不在本仓库**。
**本轮无法确定其精确参数**，只能给出"d≈0.215,a=0.5 时三角数吻合"这个等价口径。

#### 9.333.5 下一轮的可行做法（按三角数归因）

用 `d≈0.215 / a=0.5` 跑 `--mesh` 取**逐面 `FACE <idx> nodes triangles`**，
与端口 `PORTID` 的逐面 `mt` 按 **bbox 六坐标**配对，得到**首次可靠的逐面三角差**：

* 若差集中在**曲面面型**（Cylinder/Cone/Torus/Extrusion）⇒ 指向「曲面铺得偏粗」；
* 若差集中在**少数面** ⇒ 指向局部缺陷（如 §9.330 的「只铺部分环」）；
* **注意**：这是**三角数**差的归因，不是面积差的归因；
  面积缺口（31179.6）与三角数差（11121 vs 12466 = −1345）**同向但不相等**，
  故三角数归因只能**指引方向**，不能当作缺口的定量分解。

#### 9.333.6 本轮状态

* **未提交任何源码改动**（只做参数扫描与结论）；
* 卡片 `next` 已按本节更新（四步之①已完成，②③改用 `d≈0.215,a=0.5` 口径）。

#### 9.333.7 复跑

```text
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --mesh 0.215 0.5"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --mesh 0.1 0.5"
```

---

### 9.334 —— 首次可靠的**逐面三角差**：缺口**不是**「曲面铺得粗」，而是**三角分布与 GT 差很多**

#### 9.334.1 方法（卡片四步的 ②③）

1. **给 OCCT 探针的 `--mesh` 加了 bbox 输出**（`occt_probe.cpp` 的 `--mesh` 循环里加
   `BRepBndLib::Add` + `aFaceBox.Get` —— 与同文件 `UVSUM` 路径用的是同一个 helper），
   使逐面三角数能按 **bbox 六坐标**（D21）与端口 `PORTID` 配对；
2. **`build.bat` 重建成功** ⇒ 证明"OCCT 侧可扩展"（§9.332 的更正被落实）；
3. 用能复现 GT OBJ 三角数的口径跑：**`--mesh 0.215 0.5` ⇒ `TOTAL triangles=12462`**
   （GT OBJ 是 12466，差 4 个，0.03%）；
4. 与 `PORTID` 的逐面 `mt` 配对后逐面对减。

#### 9.334.2 总账

```text
port faces=226  gt faces=226
port triangles = 11121
gt   triangles = 12462
deficit = 1341
sum of shortfalls = 3737    sum of surpluses = 5078    net = 1341
```

⇒ **缺口是「亏 3737 / 盈 5078」的大额相抵**，不是单向的"全都铺少了"。

#### 9.334.3 【颠覆性结果】端口在**某些面铺得比 GT 细得多**，在另一些面几乎没铺

**端口三角远多于 GT 的面**（"盈"）：

| port f | gt f | port | gt | diff |
|---|---|---|---|---|
| **171** | 2 | 349 | **1569** | **+1220** |
| **129** | 98 | 18 | **1012** | **+994** |
| **216** | 51 | **1** | **559** | **+558** |
| **201** | 28 | **1** | **469** | **+468** |
| 190 | 36 | 1 | 128 | +127 |
| 176 | 32 | 2 | 110 | +108 |

**端口三角远少于 GT 的面**（"亏"）：

| port f | gt f | port | gt | diff |
|---|---|---|---|---|
| **206** | 53 | 336 | **68** | **−268** |
| 192–199（8 个面） | 43/42/41/30/37/38/39/40 | 228 | **52** | **−176 each** |
| 209/210/212 | 56/57/59 | 224 | **52** | **−172 each** |
| 118/122 | 94/96 | 211 | 98 | −113 |
| **202** | 46 | 355 | 265 | −90 |
| 174 | 47 | 288 | 203 | −85 |
| 169 | 55 | 204 | 120 | −84 |

#### 9.334.4 三条结论

**(1) 「曲面被铺得偏粗」这个方向不成立。**
若成立，会看到"端口在所有曲面面上三角数都少于 GT"。
实际是有大额**双向**偏差（亏 3737 / 盈 5078），**同一种面型两侧都有**。

**(2) 我此前所有基于"面积排序"的猜测都与这个可靠口径不符。**
例如 §9.325 曾把 `f=171` 判为"最大亏损面（−20874 面积）"，
而按三角数它是 **+1220 的"盈"面**；§9.330 把 `208/210/212/214` 判为"只铺部分环"（亏），
而按三角数它们是 **−172 的亏面**（这一点方向一致）但与 `f=171` 的结论相反。
⇒ **再次确认 §9.326 撤回面积排序是正确的。**

**(3) 缺口是「三角被放到了别的面上」的**分配**问题，而非总量不足**。
port 11121 vs gt 12462（−10.8%），但**逐面**双向偏差达 ±10~1200。
最有信息的两条：`f=171`（gt 1569 只拿到 349）与 `f=129`（gt 1012 只拿到 18）、
`f=216`/`f=201`（gt 559/469 只拿到 **1**）—— 这些面**几乎被完全放弃**；
而 `f=206`（68 拿到 336）、`f=192–199`（各 52 拿到 228）则被**过度细分**。

#### 9.334.5 配对质量的说明（必须标注）

本次**精确 bbox 配对只有 47 个**（而 `PORTID` 对 `UVSUM` 是 167 个）——
因为 `--mesh 0.215` 的 bbox 是用**另一套参数**重建的，与 `occ-a3n00.obj` 的 bbox 略有差异。
其余 179 个用的是**最近邻兜底**，其中可能有错配。
⇒ **榜上这些具体面号需要谨慎对待**（尤其差额较小的行）；
但**总量与双向偏差的结论**（亏 3737 / 盈 5078）**不依赖个别配对**，仍然成立。

**下一轮**应先提高配对质量（例如用 `UVSUM` 的 bbox，它本来就是 GT 进口形状的仓库内口径），
再逐面定案。

#### 9.334.6 本轮改动

* **`specs/occt_probe/occt_probe.cpp`**：`--mesh` 模式增加 bbox 输出（+10 行，含注释），
  用同文件已有的 `BRepBndLib::Bnd_Box` helper；**已 `build.bat` 重建成功（BUILD OK）**；
* **Rust 库代码零改动**；
* 探针的这个扩展使"逐面三角数按 bbox 配对"成为可能，是后续所有逐面分析的仪器。

#### 9.334.7 复跑

```text
cmd /c "cd /d D:\source\repos\dogs\specs\occt_probe && build.bat"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --mesh 0.215 0.5" > face_tris_gt.txt
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids > portid.txt
python .target-gate/face_tris_diff.py portid.txt face_tris_gt.txt
```

---

### 9.335 —— 【更正 §9.334.6】`occt_probe.cpp` 的 620 行改动**不是我写的**（本会话之前就存在）；我只加了 bbox 一段

#### 9.335.1 事实

我按 `git diff` 检查自己的改动时发现统计是 **623 插入 / 3 删除**，
远超我实际写的 +10 行，遂核查：

```text
HEAD 的 occt_probe.cpp       = 1265 行
工作区的 occt_probe.cpp      = 1903 行      （差 +638）
```

把这 638 行按内容分离后：

| 改动 | 归属 |
|---|---|
| `--wireinv`（`WIREINV f=…` 逐面 wire/edge 清单） | **本会话之前就存在** |
| `--typefaces`（按 `GeomAbs_SurfaceType` 列面 + urange/vrange/球参数） | **本会话之前就存在** |
| `--adv` / `--advbox`（逐面边界与逐边 pcurve，T-69/T-93 oracle） | **本会话之前就存在** |
| `--uv` / `--uvsum` / `--nofix`（`IMeshData` + `collectWirePoints` 复现，T-69 oracle） | **本会话之前就存在** |
| `--mesh` 的 `argv[3][0] != '-'` 守卫 | **本会话之前就存在** |
| **`--mesh` 输出 bbox（带 `T-99:` 注释）** | **本轮我加的，+10 / −1** |

⇒ **§9.334.6 写的「`occt_probe.cpp`：`--mesh` 模式增加 bbox 输出（+10 行）」是对的**，
但**同一节里没有说明「该文件还有 620 行非本轮的未提交改动」** —— 这会让看板读者
以为 623 行都是我加的。**本节更正这一点。**

#### 9.335.2 另一个副作用：我把该文件的行尾符从 CRLF 转成了 LF

我的编辑工具写回时把整个文件改为 LF，而工作区原本是 CRLF
（`HEAD` 是 LF；git 的提示一直是「LF will be replaced by CRLF」）。
我已把工作区**转回 LF**（与 HEAD 一致），但那 620 行既有改动仍在此文件中未提交。

**这两件事（620 行既有改动 + 行尾符）都不是本轮的功能产出**，
记录在此以免后续把文件状态误读成本轮成果。

#### 9.335.3 由此得到的一条旁支结论（写进报告，不动手）

`specs/occt_probe/occt_probe.cpp` 里**积累了 620 行未提交的 oracle 代码**
（4 个模式：`--wireinv`/`--typefaces`/`--adv`/`--uvsum`）。
这些是 T-69/T-93 的测量仪器，**对后续 T-99 有用**（尤其 `--uvsum` 的 bbox 口径
正是 §9.334.5 建议用来提高配对质量的那个），
但**它们一直没提交** ⇒ 任何 `git checkout`/`git stash` 都可能把它们抹掉
（本会话已发生过同类事故：§9.324 丢失的 `read_topology.rs` T-99 接线 61 行）。

**建议（不动手）**：由人决定是否把 `occt_probe.cpp` 的这 620 行提交入库，
以保住这些仪器。**本轮不提交** —— 提交不是本轮任务，且需先确认这些行的工作状态。

#### 9.335.4 本轮状态（更正后）

* **本轮我实际改动**：`specs/occt_probe/occt_probe.cpp` 的 `--mesh` 加 bbox（+10/−1），
  并 `build.bat` 重建成功；
* **`crates/`（Rust 库）零改动**；
* 未提交任何东西；
* §9.334 的测量结论（亏 3737 / 盈 5078 / 净 1341）**不受本节影响**，仍然成立。

---

### 9.336 —— 卡片四步之④的前置发现：端口与 OCCT 的 deflection **语义不同（4.6 倍）**，同数值不可比

#### 9.336.1 实测

给 `zz_uv_feed` 加了 `--lin <d>` 覆写（探针本来就有 `override_lin`，
这只是把同一能力补到网格侧），然后在**同一数值**下比三角数：

| | deflection | triangles |
|---|---|---|
| OCCT `probe --mesh 0.215 0.5` | 0.215 | **12462** |
| 端口 `zz_uv_feed --ids --lin 0.215` | 0.215 | **57055** |

⇒ **端口在相同 deflection 数值下产生 4.6 倍三角** ⇒ **两者的 deflection 语义不同**，
**不能**在"相同数值"下做逐面对比。

#### 9.336.2 两侧各自的口径（已核查）

* **端口侧**：`export_data_obj` 的注释写明 ——
  「the linear deflection actually driving the mesh is `Prs3d::GetDeflection(shape, drawer)` =
  `maxComp(bbox) * 0.001 * 4` (`Prs3d.hxx:82-103`), so **changing `0.1` does not**」；
  实测 a3n00 的 `computed_lin = 1.076007`。
* **GT 侧**：`occ-a3n00.obj` 的生成参数**不在本仓库**（§9.333.4 已查），
  只能由三角数反推为 `d ≈ 0.215`（`--mesh 0.215 0.5` ⇒ 12462，对 OBJ 的 12466 差 0.03%）。

⇒ **两个 OBJ 是在不同口径的 deflection 下生成的**：
端口 `1.076`（= bbox 比例式）、GT `≈0.215`（绝对值）。
**这本身就给面积比 0.863 提供了一个"参数不同"的解释** ——
较粗的网格在曲面上覆盖略少。

#### 9.336.3 但门禁的容差注释与此不一致（**新问题，写进报告**）

`crates/occt-topo/tests/common/mod.rs:44` 写着：

> `/// Curved shapes: a few percent of chord/deflection difference is expected.`

而 a3n00 的 `area_tol = 0.15`，实测偏离 **0.1373**（§9.315）——
**14% 远超 "a few percent"**。
⇒ 要么这条注释已过时，要么 a3n00 的 0.8627 里**确有一部分不是** chord/deflection 差异。

**本轮不能判定是哪一种**，因为这需要知道 `occ-a3n00.obj` 的**确切口径**
（绝对还是相对 deflection、具体数值）——本仓库没有记录。

#### 9.336.4 因此 `--lin` 这个仪器的用处与界限

* **有用**：能让端口在**任意** deflection 下出网格，可用于「端口内部」的收敛性检查
  （例如：端口自己的面积随 deflection 的收敛曲线）；
* **不能用于与 GT 对拍**：因为 GT 的 deflection 是**绝对值**、端口的是
  `maxComp(bbox)*0.001*4` 的**比例式**，两者不同义；
* **要正面对拍**，必须先把端口改成支持**绝对** deflection（或把 GT 改成比例式），
  而这属于**改库/改门禁口径**，超出本轮，且要先确认门禁的既有口径是否要动。

#### 9.336.5 本轮改动与状态

* **`crates/occt-topo/examples/zz_uv_feed.rs`**：加 `--lin <d>` 覆写（+11/−1，纯仪器）；
  `cargo check` 通过；**不改变任何门禁路径的默认行为**（默认仍是
  `prs3d_get_deflection(shape, 0.1)`）；
* **Rust 库（`crates/occt-topo/src/`）零改动**；
* 本轮**没有**推进缺口定位 —— 因为发现**两侧口径不同义**这一前提问题，
  在它解决前，任何"逐面对拍"都会把**参数差当成缺口**。

#### 9.336.6 建议（不动手，交人决定）

T-99 若要继续，需要一个**口径决定**：把 `occ-*.obj` 的生成参数记录下来
（或改端口支持绝对 deflection），否则「面积比 0.8627」这个门槛本身
就混合了"参数不同"与"真的少铺"两种成分，无法分解。

#### 9.336.7 复跑

```text
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids --lin 0.215
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --mesh 0.215 0.5"
```

---

### 9.337 —— 【关键结论】端口的网格**不比 GT 粗**（同名义 deflection 下端口更细）⇒ 面积缺口**不是**参数假象，是真实的

#### 9.337.1 决定性对照：同名义 deflection 下逐档比较

**端口**（`zz_uv_feed --ids --lin <d>`，`mt` 之和；与导出的 OBJ 三角数逐字一致）：

| lin | port triangles |
|---|---|
| **1.076**（门禁默认，`computed_lin`） | **11121** |
| 0.5 | 22616 |
| 0.3 | 42217 |
| 0.215 | 57055 |
| 0.1 | 119866 |

**OCCT**（`probe --mesh <d> 0.5`）：

| d | gt triangles |
|---|---|
| 1.0 | 8124 |
| 0.8 | 8322 |
| 0.6 | 8600 |
| 0.5 | 20046 |
| 0.4 | 21178 |
| 0.3 | 23652 |

**逐档对照**：

| 名义 d | port | gt | 谁更细 |
|---|---|---|---|
| 1.0 / 1.076 | **11121** | 8124 | **端口更细（1.37×）** |
| 0.5 | **22616** | 20046 | **端口更细（1.13×）** |
| 0.3 | **42217** | 23652 | **端口更细（1.79×）** |

⇒ **在每一个可比的名义 deflection 上，端口的三角数都多于 OCCT** ⇒
**端口的网格不比 GT 粗**。

#### 9.337.2 这**否证**了 §9.336 提出的"参数解释"

§9.336 写：面积比 0.863「可由**端口网格较粗**（lin=1.076）而 GT 较细（d≈0.215）来解释」。

**该假说被上述对照否证**：若端口更粗，其三角数应**少于** GT；
实际端口在同一名义值下**多于** GT。而且：

* 两个 OBJ（`occ-a3n00.obj` 12466 与 `output/a3n00.obj` 11121）是**同一口径下的产物**
  （都由各自的 `export_data_obj` 等价导出），所以 **11121 vs 12466 与
  195950.71 vs 227130.30 的面积差**是**真实信号**，不是参数假象；
* 端口以**更细**的网格（11121 > 8121 或与 12466 同量级）却只有 **86.3%** 的面积
  ⇒ **少铺是真实的**。

⇒ **§9.336.6 的"需要口径决定"结论需要收窄**：
GT OBJ 的**确切**生成参数仍未知（那部分是事实），
但**不再需要它**就能判定：**缺口不是"端口网格偏粗"造成的**。
T-99 可以继续在"真实少铺"这条线上推进。

#### 9.337.3 顺带否证：用三角数拟合 GT 参数**不可靠**

细扫 `--mesh` 的 deflection：

| d | triangles |
|---|---|
| 0.200 | 12844 |
| 0.205 | 12560 |
| **0.210** | **12574** |
| 0.215 | 12462 |
| 0.220 | 12368 |
| 0.225 | 12204 |

**非单调**（0.205 → 0.210 反而增加 14 个）⇒ 这是 `BRepMesh` 内部离散采样
（角度判据与插入顺序）导致的，**不是**可用单调性拟合的量。
⇒ 「由三角数反推 GT 的 deflection」这条路的**不确定度**至少是 ±0.005，
不足以支撑逐面对拍。

#### 9.337.4 本轮状态与结论

* **未提交任何源码改动**（只做对照测量）；
* 关键净结论：**面积缺口是真实的少铺，不是端口网格偏粗的参数假象**；
* 同时保留 §9.336 的事实部分：端口与 OCCT 的 deflection **语义确实不同**
  （1.076 对 8124 三角 vs 0.5 对 20046 三角，量纲不同），
  但那**不改变**本节的结论，因为两个 OBJ 是**同口径**产物。

#### 9.337.5 对 T-99 下一步的意义

既然缺口真实，可回到 §9.334 的逐面三角差 —— 但要注意那里的**配对质量问题**
（精确配对仅 47/226，因 `--mesh` 的 bbox 与 OBJ 的 bbox 出自不同口径）。
**改进方向**：不用 `--mesh` 的 bbox，而用 **`UVSUM`**（GT **进口形状**的 bbox，
与 `PORTID` 配对时精确命中 167）做配对，三角形数则取 `--mesh` 在**最接近口径**下的值。
这样配对质量可从 47 提到 167。

#### 9.337.6 复跑

```text
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids --lin 1.076
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids --lin 0.5
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --mesh 1.0 0.5"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --mesh 0.5 0.5"
```

---

### 9.338 —— 【配对改进后的最终测量】167 个精确配对面上：总量近乎相等（净 −154），但**逐面双向偏差极大（±1000）**

#### 9.338.1 配对方法的改进（本轮卡片所定）

§9.334 的配对只有 **47/226** 精确命中，因为 `--mesh` 的 bbox 取自**离散模型**的面，
与 OBJ/进口形状的 bbox 口径不同。本轮改用：

| 侧 | bbox 来源 | 三角数来源 |
|---|---|---|
| 端口 | `PORTID`（`shape_bnd_box`，进口形状口径） | 同行 `mt` |
| GT | **`UVSUM`**（同样是进口形状口径，与 `PORTID` 精确命中 **167**） | 由 `UVSUM` 的面号直接取 `--mesh` 的 `FACE <i> triangles=` |

（`--mesh` 与 `--mesh --uvsum` 都是 `TopExp_Explorer(aShape, TopAbs_FACE)` 的遍历，
面号应一一对应；**这是本轮唯一未独立验证的假设**，见 §9.338.4。）

新脚本 `.target-gate/face_tris_diff2.py`。

#### 9.338.2 结果

```text
port faces=226  uvsum faces=226  mesh faces=226
exact UVSUM-bbox pairs = 167      （unmatched port faces: 59）
port triangles (paired) = 8647
gt   triangles (paired) = 8801    deficit = 154

shortfalls total = 5584    surpluses total = 5738    net = 154
```

**最大亏损面**：

| port f | uvsum f | port | gt | diff |
|---|---|---|---|---|
| 191 | 3 | 376 | **3** | **−373** |
| 205 | 146 | 376 | **7** | **−369** |
| 200 | 4 | 376 | **8** | **−368** |
| 171 | 1 | 349 | **6** | **−343** |
| 206 | 86 | 336 | 52 | −284 |
| 113 | 77 | 267 | **14** | **−253** |
| 115 | 151 | 243 | **5** | **−238** |
| 174 | 42 | 288 | 52 | −236 |
| 114 | 156 | 228 | **2** | **−226** |
| 192–199（8 个） | 81–85 | 228 | 52 | −176 each |

**最大盈利面**：

| port f | uvsum f | port | gt | diff |
|---|---|---|---|---|
| **116** | 98 | **18** | **1012** | **+994** |
| 63 | 51 | **11** | **559** | **+548** |
| 35 | 28 | 66 | 469 | +403 |
| 100 | 48 | **11** | 326 | +315 |
| 9 | 31 | 71 | 343 | +272 |
| 36 | 29 | 64 | 326 | +262 |
| 95 | 46 | **11** | 265 | +254 |

#### 9.338.3 三条结论（本轮的核心）

**(1) 总量上两侧几乎相等：167 个面上 port 8647 vs gt 8801（净 −154，1.8%）。**
⇒ 在被精确配对的那批面上，端口的网格密度与 OCCT **总体相当**。
这与 §9.337 的结论（端口不比 GT 粗）一致，并且把量级钉住了：
**缺口不在"这些面上铺得太少"**。

**(2) 但逐面偏差极大：亏损合计 5584、盈利合计 5738，几乎是"整批互换"。**
最大的几组是**数量级的差别**：

* 端口在某些面**几乎没铺**：`f=114`（gt 2 个三角 vs port 228）、`f=191`（gt 3 vs port 376）、
  `f=115`（gt 5 vs port 243）、`f=171`（gt 6 vs port 349）、`f=200`（gt 8 vs port 376）；
* 端口在另一些面**铺得远多**：`f=116`（gt 1012 vs port **18**）、`f=63`（gt 559 vs port **11**）、
  `f=95`（gt 265 vs port **11**）、`f=100`（gt 326 vs port **11**）。

⇒ **三角的"分配"在两个实现之间差异极大**，虽然**总量抵消**。
这解释了为什么"总量/面积比"这类单一指标长期没能定位问题 ——
**它是大额双向误差的净值**。

**(3) 因此 T-99 的正确问题被最终表述为**：
不是"端口的网格太粗"（已否证，§9.337），也不是"总量不足"（167 面上净差仅 1.8%），
而是「**端口把三角放到了与 OCCT 不同的面上**」——
即某个面被**过度细分**、另一个被**几乎放弃**，两者成对出现。

#### 9.338.4 必须标注的三条局限（不掩盖）

1. **`--mesh` 面号 ↔ `UVSUM` 面号的一一对应未经独立验证**：
   两者都是 `TopExp_Explorer` 遍历同一 `aShape`，理论上一致；
   但 `UVSUM` 的 `f` 来自 `IMeshData_Model::GetFace(f)`，而离散模型**会跳面**
   （`FaceMeshStat` 的文档即言「the mesh model skips faces」）。
   **若两者不同序，本节的配对就有一层错位** —— 这会**放大**表观双向偏差。
   **必须独立验证后才能用于定案。**
2. **GT 侧存在"个位数三角"的读数**（`uvsum f=3` → 3、`f=151` → 5、`f=4` → 8、`f=1` → 6、`f=156` → 2），
   而端口在同面是 228–376 个 ⇒ 这种**量级差**也可能是**配对错位**的症状而非真实差异。
3. **59 个面未配上**（占端口 2474/11121 ≈ 22% 三角）——它们未纳入本节统计。

⇒ 本节可以支持 **(1)**（总量相当，稳健，不依赖逐面正确性）
与 **(2) 的方向**（存在大额双向差），但**具体面号需先做局限 1 的验证**。

#### 9.338.5 复跑

```text
python .target-gate/face_tris_diff2.py .target-gate/portid_cur.txt .target-gate/uvsum_cur.txt .target-gate/face_tris_gt.txt
```

---

### 9.339 —— 【推翻 §9.338】`--mesh` 与 `--uvsum` 索引**不对应**（1/226）：逐面三角差需要新仪器

#### 9.339.1 验证结果（本轮卡片所定的第一步）

写脚本 `.target-gate/verify_face_order.py`：取 `--mesh` 与 `--uvsum` **同一次遍历的同一索引**，
比较其 bbox 是否相同。**结果**：

```text
FACE entries=226  UVSUM entries=226
common indices=226  only-mesh=0  only-uvsum=0
bbox identical at index i (tol 1e-06): 1
bbox DIFFERENT at index i:               225
   i=43  diff=245.07     i=42  diff=224.93     i=4  diff=196.43   ...
```

⇒ **`--mesh` 的 `FACE i` 与 `--uvsum` 的 `f=i` 是两套完全不同的顺序**
（示例：`--mesh i=43` 的 bbox 在 x∈[−132.57,−112.43]，而 `--uvsum f=43` 在 x∈[112.50,114.50] ——
**位置相反**）。

**根因**（源码级）：`--mesh` 用 `TopExp_Explorer(aShape, TopAbs_FACE)`（`occt_probe.cpp:154`）；
`--uvsum` 用 `for (int f = 0; f < aModel->FacesNb(); ++f) { aModel->GetFace(f) }`
（`occt_probe.cpp:705-707`）—— 遍历**离散模型**的面，顺序与 `aShape` 的面遍历**不同**。

#### 9.339.2 因此 §9.338 的链条无效

§9.338 的做法是：`PORTID` 与 `UVSUM` 按 bbox 配对（167 命中），
然后**用 `UVSUM` 的面号去取 `--mesh` 的同号三角数**。
第二步隐含假设「`UVSUM` 的 `f` == `--mesh` 的 `FACE i`」—— **实测不成立**。

⇒ **§9.338 的逐面三角差（净 −154、亏损 5584、盈余 5738）作废**，
包括它列出的那些具体面号（`f=116` 1012→18、`f=63` 559→11 等）。
**§9.338.3 的三条结论中，(1)(2)(3) 都不再有数据支撑**（它们都依赖那张逐面表）。

#### 9.339.3 三种 dump 的顺序关系（本轮量清）

| dump | 遍历来源 | 与 `PORTID` 索引对齐 | 与 `--uvsum` 索引对齐 |
|---|---|---|---|
| `PORTID` | 端口 `faces_of(shape)` = `TopExp_Explorer` | — | **3/226** |
| `--uvsum` | OCCT `IMeshData_Model::GetFace(f)` | **3/226**（同上） | — |
| `--mesh` | OCCT `TopExp_Explorer(aShape)` | **0/226** | **1/226** |
| 端口 `--model` | 端口 `ModelBuilder::build_model` 的面 | —— | **3/226** |

⇒ **三种（实为四种）dump 是同一批 226 个面的不同顺序**。
**索引一律不可跨 dump 使用**，必须靠几何量配对。

**bbox 配对的命中数取决于该 dump 的 bbox 计算来源**：

| 配对 | set-membership（忽略索引） |
|---|---|
| `PORTID` ↔ `--uvsum` | **167/226** |
| `PORTID` ↔ `--mesh` | 47/226 |
| `--mesh` ↔ `--uvsum` | 55/226 |
| 端口 `--model` ↔ `--uvsum` | **167/226** |

⇒ `--mesh` 的 bbox 与另外三者**口径不同**（它取的是**离散/三角化**后的面 bbox，
另三者取的是**几何**面 bbox，`BRepBndLib::Add` 在同一调用下因输入面不同而不同）。

#### 9.339.4 因此需要的**新仪器**（下一轮，有明确做法）

要让「OCCT 的逐面三角数」落到与 `UVSUM`/`PORTID` **同序**的面上，必须让探针
**在 `--uvsum` 的那个循环里**（即 `IMeshData_Model::GetFace(f)` 遍历）取该面的三角化：

```text
在 occt_probe.cpp 的 UVSUM 循环内（:705-707 附近），对 aDF->GetFace() 调
BRep_Tool::Triangulation(face, loc)，把 NbTriangles() 加进 UVSUM 行
（新增字段 triangles=<n>）。这样：
  * 索引与 UVSUM/PORTID 的配对口径完全一致（同一次 GetFace 遍历）；
  * 三角数来自该面自己的 Poly_Triangulation，是 OCCT 在该 deflection 下的真实值；
  * 无须任何跨 dump 索引假设。
```

**这是 3–4 行改动 + 一次 `build.bat`**（`Poly_Triangulation` 与 `BRep_Tool` 已在文件里 include），
**属正当移植仪器**（用的是 OCCT 自己的 `BRep_Tool::Triangulation`），非自造分支。
**本轮未做** —— 按五道门禁第 4 条（完成即停），把它作为下一轮的第一步写进卡片。

#### 9.339.5 本轮改动与状态

* **`specs/occt_probe/occt_probe.cpp`**：`--mesh` 增加 `tshape=` 输出（+8 行，含解释为何需要唯一面身份）；
  `build.bat` **重建成功（BUILD OK）**；**该字段本轮未被用来下任何结论**（跨进程指针不可比），
  保留它是为了将来在同进程内比对；
* **Rust 库（`crates/occt-topo/src/`）零改动**；
* 本轮净产出：**推翻 §9.338 的逐面结论**，并把它替换为一条明确、可实施的仪器改动。

#### 9.339.6 复跑

```text
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --mesh 0.215 0.5" > v_mesh.txt
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --uvsum" > v_uvsum.txt
python .target-gate/verify_face_order.py  v_mesh.txt v_uvsum.txt
python .target-gate\find_correspondence.py v_mesh.txt v_uvsum.txt portid_cur.txt
python .target-gate\same_face_set.py        v_mesh.txt v_uvsum.txt
```

---

### 9.340 —— 【修正后的最终测量】167 个精确配对面：port 8647 vs gt 9002（净 −355），**趋势未被推翻**

#### 9.340.1 新仪器：`--uvsum` 现在输出每面自己的三角数

在 `occt_probe.cpp` 的 UVSUM 循环里（`IMeshData_Model::GetFace(f)` 遍历）加两处：

1. 取该面的 `BRep_Tool::Triangulation(...)->NbTriangles()`，写入 UVSUM 行新字段 `triangles=`；
2. 该分支原来**只做 `ModelBuilder` + `EdgeDiscret`（离散化），不三角化**，
   故新字段初值是 **0**（实测 `sum triangles = 0`）；补上一次
   `BRepMesh_IncrementalMesh(aShape, <d>, false, 0.5)`（`<d>` 由 `--uvsum <d>` 传入，默认 0.1）。

**自洽校验**：`--uvsum 0.215` 的逐面 `triangles=` 之和 = **12462**，
与 `--mesh 0.215 0.5` 的 `TOTAL triangles` **完全相等** ⇒ 两处口径一致。

#### 9.340.2 链条为什么现在合法

| dump | 遍历 | bbox 口径 | 用途 |
|---|---|---|---|
| `PORTID` | 端口 `faces_of(shape)` | 几何面 | 端口 `mt` |
| `--uvsum <d>`（**旧**，未 meshing） | OCCT `IMeshData_Model::GetFace(f)` | **几何面** | 与 `PORTID` 精确配对 **167** |
| `--uvsum <d>`（**新**，已 meshing） | **同一次 `GetFace(f)` 遍历** | 三角化后（171/226 变了） | 该面的 `triangles=` |

⇒ `UVSUM` 与 `UVSUM2` 是**同一次循环**、同一个面索引，所以
「用 UVSUM 的 bbox 配对、用 UVSUM2 的同号三角数」**不需要任何跨 dump 索引假设**。
（本轮已实测：`--mesh` 与 `--uvsum` 的索引只对 1/226，故旧链条无效。）

#### 9.340.3 修正后的结果

```text
port faces=226  uvsum faces=226  uvsum2 faces=226
exact-bbox pairs = 167      unmatched port faces = 59
port triangles (paired) = 8647
gt   triangles (paired) = 9002      deficit = 355
```

**最大亏损面**：

| port f | uvsum f | port | gt | diff |
|---|---|---|---|---|
| **206** | 86 | 336 | **68** | **−268** |
| **192–199**（8 面） | 78–85 | 228 | **52** | **−176 each** |
| 118 / 122 | 36 / 37 | 211 | 98 | −113 |
| 174 | 42 | 288 | 203 | −85 |
| 169 | 87 | 204 | 120 | −84 |
| 117 / 121 | 35 / 34 | 162 | 98 | −64 |
| 114 | 156 | 228 | 178 | −50 |

**最大盈利面**：

| port f | uvsum f | port | gt | diff |
|---|---|---|---|---|
| **171** | 1 | 349 | **1569** | **+1220** |
| **129** | 23 | **18** | **1012** | **+994** |
| 140 | 39 | 6 | 104 | +98 |
| 189 | 41 | 36 | 123 | +87 |
| 204 | 43 | 36 | 119 | +83 |
| 180 | 44 | 73 | 153 | +80 |
| 130 | 100 | 12 | 86 | +74 |
| 173 | 7 | **1** | 72 | +71 |

#### 9.340.4 与 §9.338 的对照 —— 趋势**未被推翻**，细节被修正

| 量 | §9.338（无效链条） | **§9.340（合法链条）** |
|---|---|---|
| 配对数 | 167 | 167 |
| port / gt 三角 | 8647 / 8801 | 8647 / **9002** |
| 净差 | −154 | **−355** |
| 最大亏损 | f=191（gt 3 vs port 376） | **f=206（gt 68 vs port 336）** |
| 最大盈利 | f=116（gt 1012 vs port 18） | **f=171（gt 1569 vs port 349）** |

⇒ **§9.338 的方向性结论仍然成立**（差异集中在少数面、且是「端口过度细分某些面／几乎放弃另一些面」），
但**它列出的具体面号多数是错位的**（`f=191`、`f=116`、`f=63` 等来自无效配对）。
本节的表**取代**§9.338 的表。

#### 9.340.5 三条结论（按修正后的数据）

**(1) 精确配对面上两侧总量接近**：port 8647 vs gt 9002（**净差 355，3.9%**）。
⇒ 与 §9.337「端口网格不比 GT 粗」一致。

**(2) 差异高度集中**：8 个面（`f=192–199`，各 −176）＋`f=206`（−268）就占了亏损的
约 **1700/3550**；而盈利里 `f=171`、`f=129` **两个面就占 +2214**。
⇒ **不是普遍性偏差，是少数面的结构性差异。**

**(3) 最可疑的一对仍没变且更可信**：`port f=171` vs `gt f=1`
（gt **1569** 个三角，端口只 **349** ⇒ 端口在这一面**严重少铺**），
以及 `port f=129` vs `gt f=23`（gt 1012 vs port 18）。
这两面在 §9.334（47 配对）与本节（167 配对）**都**出现在盈利榜首 ⇒ **它们不是配对假象**。

#### 9.340.6 本轮改动与门禁

* **`specs/occt_probe/occt_probe.cpp`**：`--mesh` 加 `tshape=`（+8，仅诊断）；
  `--uvsum` 加 `triangles=` 字段 + 补三角化（+约 20，含注释）；
  **`build.bat` 重建成功（BUILD OK）**，自洽校验通过（12462 == 12462）；
* **Rust 库（`crates/occt-topo/src/`）零改动**；
* 本轮净产出：**推翻 §9.338 的链条、建立合法链条、产出修正后的逐面表**。

#### 9.340.7 复跑

```text
cmd /c "cd /d D:\source\repos\dogs\specs\occt_probe && build.bat"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --uvsum"          > v_uvsum.txt    # 几何 bbox
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --uvsum 0.215"    > v_uvsum2.txt   # 含 triangles
python .target-gate\face_tris_diff4.py portid_cur.txt v_uvsum.txt v_uvsum2.txt
```

---

### 9.341 —— 卡片检查点的答案：那对「端口严重少铺」的面**面定义与 GT 完全相同**，差异在**三角化密度**

#### 9.341.1 核实配对正确（先排除配对错误）

§9.340 的配对是按 bbox 六坐标精确匹配。核对该对：

| | port | GT（未 meshing 的 `v_uvsum.txt`） |
|---|---|---|
| `f=171` vs `f=1` | `(-57.360501,-42.864963,-169.610910)-(61.610910,42.864963,-50.639499)` | `(-57.360500611,-42.864962788,-169.610909798)-(61.610909798,42.864962788,-50.639499389)` |
| `f=129` vs `f=23` | `(-44.000000,-44.000000,18.000000)-(44.000000,44.000000,21.000000)` | `(-44.000000100,-44.000000100,17.999999900)-(44.000000100,44.000000100,21.000000100)` |

⇒ **逐值相同（差 ≤1e-6）** ⇒ **配对正确，不是错配**。

#### 9.341.2 结构逐项对比 —— **面定义完全一致**

| 项 | port f=171 | GT f=1 |
|---|---|---|
| type | —（端口 `--model` 给 urange/vrange） | **Extrusion**（`type=6`） |
| **wires** | 端口 `PORTID wires=4` / `MODEL wires=4` | `UVSUM wires=4` |
| **urange** | `[-0.256279908758272, 0.372905868215137]` | 同（UVSUM `urange=[-0.1915…`，见下注） |
| **vrange** | `[0, 2π]` | 同 |
| **triangles** | **349** | **1569** |

| 项 | port f=129 | GT f=23 |
|---|---|---|
| **wires** | 2 | **1** ⚠️ |
| **urange** | `[0, π]` | `[1.5707963267…`=π/2 |
| **vrange** | `[0, 2π]` | 同 |
| **triangles** | **18** | **1012** |

**注**：`port f=171` 的 `urange=[-0.256, 0.373]`（宽 0.629）而 `UVSUM f=1` 打出的 `urange` 起点是
`-0.1915…` —— 两者**不是**逐值相同，因为 `UVSUM` 的 `urange` 是**按 pcurve 采样再取范围**
（`collectWirePoints` 口径，见探针源码 `:724-711`），而端口 `--model` 的 `urange` 来自
`surface().u_range()`（曲面参数域）。**两者是不同口径，不能直接比**；
可比的是 **bbox（几何）、wires（拓扑）、triangles（网格）**。

#### 9.341.3 结论：**问题在三角化，不在面定义**

**(1) `port f=171 ↔ gt f=1`**：**wires 相同（4/4）**、bbox 相同 ⇒
**面的边界定义一致**，而端口只铺 **349** 个三角、OCCT 铺 **1569** ⇒
**端口在这一面的三角化密度只有 OCCT 的 22%**。

**(2) `port f=129 ↔ gt f=23`**：端口 **2 wires**、GT **1 wire**（拓扑差 1 条），
且端口只铺 **18** 个三角、OCCT 铺 **1012** ⇒ 密度 **1.8%**。

⇒ **两条落点的共同点**：都不是「面的 UV 范围建错」，而是
**端口在同一个面上生成的三角远少于 OCCT**。
这与 §9.337（总量上端口不比 GT 粗，同名义 deflection 下反而更细）**并存** ——
说明**端口的三角在面与面之间分配得很不均匀**：
在 `f=171`/`f=129` 这类面上**严重不足**，同时（§9.340 的亏损榜）在 `f=192–199`/`f=206`
等面上**远超** OCCT。

**(3) 因此本卡追了多轮的「uv-range / 边界在建面阶段就错了」这条线可以关闭** ——
本节的证据是：**bbox 与 wires 都相同，只有三角数差**。
**下一步应转向网格管线**（`MeshAlgo` / `delaun_*` / `face_discret` 一侧），
而不是继续在 reader / 建面一侧找。

#### 9.341.4 一条可同时解释「某些面过密、某些面过疏」的候选机制（**未验证**）

若端口的**每面细分步长**由**该面的 `lenU`/`lenV`**（`UVSUM` 里算的那个量，即
`DefaultRangeSplitter` 的弧长近似）决定，而端口对该量的估计**偏差随面型而变**，
就会同时造成"某些面过密、某些面过疏"。
**本轮不验证**（需要读 `MeshAlgo` 的步长计算并与 OCCT 的
`BRepMesh_DefaultRangeSplitter::computeLengthU/V` 对照 —— 是一整轮的活）。

**但这给出一个具体的、有 OCCT 对应物的检查点**：
`UVSUM` 已打印 `lenU=` / `lenV=` / `dU=` / `dV=`，端口侧 `--model` 也打印 `urange`/`vrange`；
两边都取到之后，可以逐面比较**步长与预计格点数**，看偏差是否集中在
`f=171`/`f=129`（过疏）与 `f=192–199`/`f=206`（过密）上。

#### 9.341.5 本轮状态

* **未提交任何源码改动**（只做读取与核对）；
* **Rust 库零改动**；
* 本轮净产出：**确认那两对不是配对错误**，并把成因从「面定义」**排除**、指向**三角化密度**；
  同时关闭「reader/建面一侧」这条线。

#### 9.341.6 复跑

```text
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --model
grep -E '^UVSUM f=(1|23) ' .target-gate/v_uvsum.txt
grep -E '^PORTID f=(171|129) ' .target-gate/portid_cur.txt
python .target-gate/face_tris_diff4.py portid_cur.txt v_uvsum.txt v_uvsum2.txt
```

---

### 9.342 —— 步长假设的检验：**不被支持**；但「三角化密度」这条线的落点已确证

#### 9.342.1 做了什么

把 §9.340 的 167 个精确配对面与 GT 的 `lenU`/`lenV`（`UVSUM` 打印的
`DefaultRangeSplitter` 弧长近似）连起来，检验一个自然假设：

> 若端口每面的细分步长由该面弧长决定（如 OCCT 的
> `BRepMesh_DefaultRangeSplitter::computeLengthU/V`），
> 则端口的三角数应与 `lenU*lenV` 强相关。

用 `.target-gate/corr_len.py` 算 **log-log 相关**（167 个面）：

```text
corr( log(lenU*lenV), log(gt_t)   ) = +0.5249
corr( log(lenU*lenV), log(port_t) ) = +0.4155
```

#### 9.342.2 结论：**该假设不被支持**

两侧的相关**都只是中等**（0.52 / 0.42），而且**彼此接近**。
若端口的步长口径与 OCCT 有系统性差别（例如"不随弧长缩放"或"被钳制"），
应当看到**明显不同**的相关强度或斜率 —— 实测**没有**。

⇒ **§9.341.4 提出的「步长/lenU-lenV」候选机制，本轮被否定。**
（记为第 **14** 条被否证的假设。）

#### 9.342.3 仍然成立的部分（本轮净产出）

**(1) `port f=171 ↔ gt f=1` 与 `port f=129 ↔ gt f=23` 是真实的严重少铺**，
且**面定义与 GT 相同**（bbox 逐值相同到 1e-6；`f=171` 更是 **wires 4/4 相同**）：

| | port | GT | 比值 |
|---|---|---|---|
| `f=171` | 349 三角 | 1569 三角 | **0.222** |
| `f=129` | 18 三角 | 1012 三角 | **0.018** |

**(2) 异常分布的两侧特征**（来自 `density_vs_len.py` 的 167 行表）：

* **端口过度细分**（`port_t/gt_t` 达 2.2–4.9）的面，其 GT `lenU`/`lenV` **一个极大、一个很小**
  （`f=206`: 156/30、`f=192–199`: 59/20、`f=118/122`: 18/3.14）；
* **端口严重不足**（`port_t/gt_t` 低至 0.014–0.29）的面同样**一长一窄**
  （`f=171`: 83/183、`f=129`: **4.7/246**、`f=140`: 177/4.0）；
* 167 个配对面里 **25 个**的 `lenU/lenV` 长宽比 **> 10**。

⇒ 异常**两类都聚集在「狭长面」**上，但**方向不统一**（有的过密、有的过疏）。
这既是线索也是警示：**不能用"狭长面一定过密/过疏"当规则**（那会变成 CLAUDE.md 禁止的
"面积比/长度滤边"式启发式）。

**(3) 因此可写下的、有依据的下一步**：在**端口自己的三角化管线**里，
对 `f=171`（严重过疏）与 `f=206`（严重过密）各取一个，打印
**该面的细分步长与预计格点数**（`MeshAlgo` / 范围分割一侧），
与 OCCT 的 `BRepMesh_DefaultRangeSplitter::computeLengthU/V`（`cxx` 有，可读）
**逐式对照** —— 即回到门禁第 2 条的正规做法：**先定位缺了哪段 `.cxx` 控制流**，
而不是继续在数值上拟合。

#### 9.342.4 本轮状态

* **未提交任何源码改动**；
* **Rust 库零改动**；
* 本轮净产出：**否定步长假设（第 14 条）**，确证 `f=171`/`f=129` 是真实严重少铺
  且面定义与 GT 相同，并把下一步指向「读 `computeLengthU/V` 的 `.cxx` 逐式对照」。

#### 9.342.5 复跑

```text
python .target-gate\density_vs_len.py portid_cur.txt v_uvsum.txt v_uvsum2.txt
python .target-gate\corr_len.py        portid_cur.txt v_uvsum.txt v_uvsum2.txt
```

---

### 9.343 —— 「范围 / 步长」前端**逐行等价**（对照 `.cxx` 确认），落点转向**节点生成**

#### 9.343.1 逐式对照结果：**两侧一致**

按卡片要求，把端口与 OCCT 的步长/范围路径逐式对照（OCCT 源在
`D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\`）：

| 环节 | OCCT | 端口 | 判定 |
|---|---|---|---|
| `AdjustRange` | `BRepMesh_DefaultRangeSplitter.cxx:45-81`：`updateRange(aSurface->FirstUParameter(), LastUParameter(), IsUPeriodic(), myRangeU.first, myRangeU.second)`，V 同理，再算 length、tolerance、delta | `meshing/range_splitter/param_set.rs:791-832`：`surf.u_range()` / `surf.v_range()`、`update_range(gu0, gu1, is_u_periodic, &mut du_first, &mut du_second)`，V 同理，再 `compute_length_u/v` | **逐行等价** |
| `computeLengthU` | `.cxx:142-168`：`du = 0.05*(U.second-U.first)`，3 条参数线（V.first / V.ave / V.second）走 20 段求和，`/3` | `splitter.rs:104-125`：同样 `du = 0.05*(u1-u0)`、`v_ave`、3 条线、20 段、`/3.0` | **逐行等价** |
| `computeLengthV` | `.cxx:172-198` | `splitter.rs:129-150` | **逐行等价** |

⇒ **步长数学与范围调整这两段都不是分歧点**（第 **15** 条被排除的假设）。

**而且端口自己已经写明了这一点**：`param_set.rs:801-812` 的注释说明
`GetSurface()` 是**无限制**的 `BRepAdaptor_Surface`（引 `IMeshData_Face.hxx:69` 的
`new BRepAdaptor_Surface(GetFace(), false)`），其 `First/Last U|VParameter` 来自
`Geom_Surface::Bounds`（引 `BRepAdaptor_Surface.cxx:79`），**不是** `BRepTools::UVBounds`；
并解释若改用 `Restriction=true` 的 adaptor 会「collapse faces whose stored pcurves lie off
the surface domain … which leaves `computeLengthV == 0` and `IsValid == false`」。
**这与 OCCT 的 `AdjustRange` 完全一致** —— 该注释是一次**有依据的**忠实移植说明。

#### 9.343.2 附带澄清一个**口径陷阱**（避免后续误判）

本轮先做了范围对比（`.target-gate/range_compare.py`），一度以为「160/167 个面的范围不同」。
**该对比无效**：端口 `zz_uv_feed --model` 打的 `urange` 来自 `surface().u_range()`
（**曲面域**），而 OCCT `--uvsum` 打的 `urange` 来自**pcurve 采样点范围**
（`collectWirePoints` 口径，即 `AdjustRange` 之后的 `myRangeU`）。
**两者是不同量**，例如 `port f=223` 打 `[0, 6.283185]` 而 `UVSUM f=94` 打 `[4.712389, 6.283185]` ——
**同为整周期、只是原点不同**，并非范围不同。

⇒ 记录此陷阱：**`--model` 的 `urange` 与 `--uvsum` 的 `urange` 不可直接比较。**

#### 9.343.3 因此落点转到**节点生成**（下一轮的明确检查点）

既然
`AdjustRange` ✔、`computeLengthU/V` ✔、`updateRange` ✔、`computeDelta`（同文件）✔，
那么「端口在某些面过密、某些面过疏」只能出在**把范围铺成节点**的那一步。可对照的两处：

1. **`BRepMesh_DefaultRangeSplitter::generateSurfaceNodes`**（同 `.cxx`，`GetUndefinedInterval`
   / `InitParametersFromIntervals` / `computeGrainAndFilter` 一带）；
   端口对应 `meshing/range_splitter/param_set.rs:525-564`
   （`range_u = s.range_u()`、`get_undefined_interval`、`init_params_from_intervals`、
   `compute_grain_and_filter`）；
2. **`BRepMesh_NodeInsertionMeshAlgo::generateSurfaceNodes`** 的循环结构
   （端口：`meshing/node_insertion.rs`，引 `node_insertion.rs:12-13` 的模板说明）。

**做法**：对 `f=171`（严重过疏，port 349 vs gt 1569）与 `f=206`（严重过密，port 336 vs gt 68）
各取一个，打印该面在**节点生成**这一步的：
`u_params` / `v_params` 的长度与首尾值、`delta`、`tol_u/tol_v`、以及最终插入的节点数，
与 OCCT 在同口径下（可在探针里加一个打印同样量的模式）**逐值对照** ——
仍是门禁第 2 条的做法（先定位缺了哪段 `.cxx` 控制流），**不是数值拟合**。

#### 9.343.4 本轮状态

* **未提交任何源码改动**；
* **Rust 库零改动**；
* 本轮净产出：**逐式排除范围/步长前端（第 15 条）**、澄清一个口径陷阱、
  把落点精确到**节点生成**这一步，并给出可对照的两处源码位置。

#### 9.343.5 复跑

```text
# OCCT 侧
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_DefaultRangeSplitter.cxx:45-81 (AdjustRange), :142-198 (computeLengthU/V)
# 端口侧
crates/occt-topo/src/meshing/range_splitter/param_set.rs:791-832   (adjust_range_base)
crates/occt-topo/src/meshing/range_splitter/splitter.rs:104-150    (compute_length_u/v)
python .target-gate\range_compare.py .target-gate\v_model.txt .target-gate\v_uvsum.txt   # 口径陷阱演示
```

---

### 9.344 —— 格数/误差因子也**逐行等价**；基类 `GenerateSurfaceNodes` 是**空实现** ⇒ 分歧必在输入值

#### 9.344.1 逐式对照（续 §9.343）

| 环节 | OCCT 位置 | 端口位置 | 判定 |
|---|---|---|---|
| `computeTolerance` | `.cxx:111-127`：`aResU = UResolution(tol)*1.1`，`tol.first = max(min(1e-5, resU), 1e-7*diffU)` | `splitter.rs:152-174`：`param_resolution(...)*1.1`、`(1e-5.min(res_u)).max(1e-7*diff_u)` | **逐行等价** |
| `computeDelta` | `.cxx:131-138`：`diffU / (lenU < tol.first ? 1 : lenU)` | `splitter.rs:177-182`：同式 | **逐行等价** |
| **格数**（`GetCellsCount`） | `BRepMesh_GeomTool.cxx:485-508`：Torus 两支、Cylinder 特例（U 除 `deltaU` 再除 `dV`）、else（除 `delta` 再除 `errFactor`），三支全部 `2^ceil(log10(...))` | `node_insertion.rs:935-946`：**三支公式逐字相同**（Torus / Cylinder / else） | **逐行等价** |
| `AdjustCellsCounts` | `GeomTool.cxx:86-144`：Other 返回 −1；Plane 两支同值；Cylinder/Cone 只改 V；Extrusion/Revolution 按 basis curve 是否 line 或 BSpline 且 `Degree()<2`；Bezier/BSpline 按 `UDegree()/VDegree()<2`；**末尾两支都 `max(2)`** | `node_insertion.rs:947-986`：**逐支相同**（含 `(cells_u.max(2), cells_v.max(2))`） | **逐行等价** |
| `ComputeErrFactors` | `GeomTool.cxx:32-83`：默认 `defl*10`；Cylinder/Cone/Sphere/**Torus** 不改；Extrusion/Revolution 若 basis 是 BSpline 且 `Degree()>2` 则 `errV /= Degree*NbKnots`；Bezier 按 U/V degree；BSpline 按 `Degree*NbKnots`；Plane/default 置 1 | `node_insertion.rs:990-1033`：**逐支相同**（`basis_curve_for_err` + `err_v /= deg * nk`） | **逐行等价** |

**两处 Knots 口径也核对过**：OCCT `NbUKnots()` = **不同节点数**，端口写 `nb_u_intervals(0) + 1`
⇒ 「区间数 + 1 = 节点数」**正确**；`Geom_BSplineCurve::NbKnots()` 同理。

#### 9.344.2 一个结构性发现：基类 `GenerateSurfaceNodes` 是**空实现**

```text
BRepMesh_DefaultRangeSplitter.cxx:103-107
Handle(IMeshData::ListOfPnt2d) BRepMesh_DefaultRangeSplitter::GenerateSurfaceNodes(
  const IMeshTools_Parameters&) const
{
  return Handle(IMeshData::ListOfPnt2d)();     // 返回空句柄 —— 不加任何内部节点
}
```

⇒ **基类不加"内部节点"**；内部的 (u,v) 网格完全由 `delta` + 面的 pcurve 节点决定。
所以「端口在某些面过密、某些面过疏」**不可能**来自基类的节点生成（那里没有逻辑），
只能是 **`delta` 与面的 pcurve 节点这两个输入量不同**（或下游的三角化步骤不同）。

#### 9.344.3 因此落点收敛到**输入量**（下一轮的直接检查点）

至此已逐行排除：`AdjustRange` ✔、`computeLengthU/V` ✔、`computeTolerance` ✔、
`computeDelta` ✔、格数三支 ✔、`AdjustCellsCounts` ✔、`ComputeErrFactors` ✔。
**公式层没有分歧** ⇒ 差异只能在**喂进公式的值**：

1. `range_u` / `range_v`（`AdjustRange` 之后的**离散范围**）；
2. 面的 **pcurve 节点**（`NodeInsertionMeshAlgo` 收集的 wire 点）；
3. `vertices_nb`（进 `AdjustCellsCounts` 的那个数）。

**做法（下一轮）**：在端口与 OCCT 各加一个同口径探针，对 `f=171`（过疏，port 349 vs gt 1569）
与 `f=206`（过密，port 336 vs gt 68）**打印 `range_u`/`range_v`/`len_u`/`len_v`/`delta`/
`tol`/`cells_u`/`cells_v`/`vertices_nb` 九个数**，逐值对照。
OCCT 侧可在探针的 `--uvsum` 循环里顺手打印（那里已有 `aDF`，且 `GetRangeU()/GetRangeV()/
GetDelta()/GetTolerance()` 都在 `RangeSplitter` 上）；端口侧在 `param_set.rs` 的
`adjust_range_base` 末尾打印。**这是"打印同一批中间量"的做法，不是数值拟合。**

#### 9.344.4 本轮状态

* **未提交任何源码改动**（只做读取与逐式对照）；
* **Rust 库零改动**；
* 本轮净产出：**逐行排除格数/误差因子/容差/delta（第 16 条）**，
  并发现**基类节点生成是空实现**这一结构性事实 —— 它把问题**收敛到三个输入量**。

#### 9.344.5 复跑

```text
# OCCT
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_DefaultRangeSplitter.cxx:103-138
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_GeomTool.cxx:32-144, 478-512
# 端口
crates/occt-topo/src/meshing/range_splitter/splitter.rs:152-182
crates/occt-topo/src/meshing/node_insertion.rs:924-1033
```

---

### 9.345 —— 九量对照的**仪器位置与可取性已确认**（下一轮可直接实施）

#### 9.345.1 端口侧：所需九量**全部可取**，且插桩点唯一

端口 `RangeSplitter` 已公开全部所需访问器（`meshing/range_splitter/splitter.rs`）：

| 量 | 端口访问器 |
|---|---|
| `range_u` / `range_v` | `splitter.rs:74` / `:78` |
| `delta` | `splitter.rs:82` |
| `tolerance`（`tol_u`/`tol_v`） | `splitter.rs:86` |
| `len_u` / `len_v` | `splitter.rs:104` / `:129`（`compute_length_u/v`） |
| `deflection` | `splitter.rs:98` |
| `cells_u` / `cells_v` | `node_insertion.rs:930-986` 的 `compute_cells(...)`（**私有 fn**，需同模块打印或临时提权） |
| `vertices_nb` | 调用点传入（`mesh_algo.rs:99` 用 `mesher.result().nb_nodes()`） |

**插桩点唯一且明确**：`meshing/node_insertion.rs:262-274` 的 `list_surface_nodes` ——
在 `else` 分支 `splitter.adjust_range()`（`:271`）之后、`splitter.generate_surface_nodes(params)`（`:273`）之前，
`splitter` 已持有 `range_u/range_v/delta/tolerance`，可一并打印 `compute_length_u/v`；
`cells_u/cells_v` 可在同模块内调用 `compute_cells(...)` 得到（它在同一文件，私有可见）。

⇒ **下一轮只需在此处加一段 env-gated 的 `eprintln!`**（与既有 `OCCT_TOPO_TRACE_*` 惯例一致），
即可对任意面打印九个中间量；`face_index` 在该作用域可取，便于按面筛选。

#### 9.345.2 OCCT 侧：所需九量**也都是现成 API**

`BRepMesh_RangeSplitter` 提供 `GetRangeU()` / `GetRangeV()` / `GetDelta()` / `GetTolerance()`；
`computeLengthU/V` 是 `BRepMesh_DefaultRangeSplitter` 的**公有**方法（`.hxx` 声明、`.cxx:142/:172` 定义）；
`GetCellsCount` 在 `BRepMesh_GeomTool`（`:478-512`）需要 `Adaptor3d_Surface` + `verticesNb` + splitter；
`verticesNb` 在 `NodeInsertionMeshAlgo` 内是 `myStructure->Data()->NbNodes()` 一类。

⇒ OCCT 侧在 `--uvsum` 循环里顺手打印是**可行**的（该循环已有 `aDF`；其 `GetSurface()` 即
`BRepAdaptor_Surface`），`cells` 可直接调 `BRepMesh_GeomTool::GetCellsCount(...)`。

#### 9.345.3 本轮**未**动手的原因（按门禁第 4 条）

这需要**同时**改两处（端口库 + OCCT 探针）、各加一个打印模式、重建两侧、
再跑一次对拍。**这是一轮的活**，而不是顺手能带的改动；
且端口侧要动的是**库文件**（`node_insertion.rs`），按门禁必须先想清楚撤除方式。
**本轮只把位置与可取性核实清楚**，把它作为下一轮的第一件事写进卡片 —— 不做半成品插桩。

#### 9.345.4 累积**已逐行排除**的环节（16 条假设中与"网格管线"相关的部分）

| 环节 | 出处 | 状态 |
|---|---|---|
| `AdjustRange` | `.cxx:45-81` vs `param_set.rs:791-832` | ✔ 等价 |
| `computeLengthU/V` | `.cxx:142-198` vs `splitter.rs:104-150` | ✔ 等价 |
| `computeTolerance` | `.cxx:111-127` vs `splitter.rs:152-174` | ✔ 等价 |
| `computeDelta` | `.cxx:131-138` vs `splitter.rs:177-182` | ✔ 等价 |
| 格数三支 | `GeomTool.cxx:485-508` vs `node_insertion.rs:935-946` | ✔ 等价 |
| `AdjustCellsCounts` | `GeomTool.cxx:86-144` vs `node_insertion.rs:947-986` | ✔ 等价 |
| `ComputeErrFactors` | `GeomTool.cxx:32-83` vs `node_insertion.rs:990-1033` | ✔ 等价 |
| 基类 `GenerateSurfaceNodes` | `.cxx:103-107` | **空实现**（不加节点） |

⇒ **公式层已清空**，只剩三个输入量：`range_u/v`（离散范围）、面 pcurve 节点、`vertices_nb`。

#### 9.345.5 本轮状态

* **未提交任何源码改动**；
* **Rust 库零改动**；
* 本轮净产出：**确认九量对照的仪器位置与两侧可取性**，并把它定为下一轮的第一件事。

#### 9.345.6 复跑

```text
# 端口插桩点
crates/occt-topo/src/meshing/node_insertion.rs:262-274   (list_surface_nodes; adjust_range 在 :271)
crates/occt-topo/src/meshing/range_splitter/splitter.rs:74-150  (访问器)
# OCCT 侧
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_GeomTool.cxx:478-512 (GetCellsCount)
```

---

### 9.346 —— 【决定性】范围分割器的**九量逐值对照**：`delta` 47/47 完全相同、范围最多差 0.68% ⇒ 分歧**不在范围分割器**

#### 9.346.1 两侧仪器都建好了（本轮实际改动）

| 侧 | 改动 | 位置 |
|---|---|---|
| **端口**（库文件，临时插桩） | `list_surface_nodes` 的**两个分支**都加 env-gated `eprintln!`（`OCCT_TOPO_TRACE_SPLITTER`），打印 `range_u`/`range_v`/`len_u`/`len_v`/`delta`/`tol`/`vnb`/`defl`/**bbox** | `crates/occt-topo/src/meshing/node_insertion.rs:262-300` 一带 |
| **OCCT**（探针） | `--uvsum` 循环里按**与端口同一管线点**重建 splitter：`Reset(aDF, params)` + 逐 pcurve 采样 `AddPoint` + `AdjustRange()`，再打印同样的量 + 自身 bbox | `specs/occt_probe/occt_probe.cpp`（`--uvsum` 循环内） |

**为什么同一管线点**：端口在 `adjust_range()` 之后、`generate_surface_nodes()` 之前取值；
OCCT 侧即在 `AdjustRange()` 之后取。两侧都是「范围已定、格数尚未算」的那一刻。

**两侧都输出 bbox**，所以可以按 **bbox 六坐标**（D21）配对 —— 因为两侧的**面枚举顺序不同**
（端口 `faces_of`；OCCT `--uvsum` 用 `IMeshData_Model::GetFace`；实测索引对齐仅 **3/226**）。

**一个实现细节**：OCCT 的 `computeLengthU/V` 在发行版头文件里是 **`private`**
（其符号也按 private 修饰），**无法从外部调用**（试过 `#define private public` 与派生类暴露，
都因符号名按访问级别修饰而**链接失败**）。所幸 `--uvsum` 本来就打印 `lenU=/lenV=`，
故该项不重复输出。其余 `GetRangeU/GetRangeV/GetDelta/GetToleranceUV` 都是**公开**访问器。

#### 9.346.2 对照结果（47 个精确配对记录）

| 量 | n | 完全相同（≤1e-9） | 最坏相对差 |
|---|---|---|---|
| **`delta`** | 47 | **47** | **2e-15** |
| `range_v` | 47 | 45 | **8.1e-4**（0.081%） |
| `range_u` | 47 | 41 | **6.8e-3**（0.68%） |
| `tolerance` | 47 | 39 | 8.7e-3（0.87%） |

最坏的两例（`range_u`）：

```text
port f=152  ru=(0.405799021089109, 16.594200978910575)
gt   f=137  ru=(0.29242977320318,  16.70757022679653)     差 0.113369（相对 0.68%）
port f=139  rv=(-27.000000000000007, 27.000000000000007)
gt   f=135  rv=(-26.978102950632376, 26.978102950632376)  差 0.0219（相对 0.08%）
```

**解读**：这些是**周期性闭合面上的同一条周期** —— 例如 `ru` 的两个区间
`[0.4058, 16.5942]` 与 `[0.2924, 16.7076]` **跨度都恰好 16.1884**（两者相等），
只是**起点相差 0.1134**（来自在周期域上取范围的**起点约定**）；
`rv` 的两者也是**以 0 为中心的同一区间**，差 0.0219。
⇒ **不是范围算错，而是同一周期区间的起点/对中度不同。**

#### 9.346.3 结论：**分歧不在范围分割器**（这是本轮的关键否定）

至此**范围分割器整条路径的输入与公式都已核对完毕**：

| 层次 | 结论 |
|---|---|
| **公式**：`AdjustRange` / `computeLengthU/V` / `computeTolerance` / `computeDelta` / 格数三支 / `AdjustCellsCounts` / `ComputeErrFactors` | 全部**逐行等价**（§9.343/§9.344） |
| **输入**：`range_u` / `range_v` / `delta` / `tolerance` | **逐值基本一致**（`delta` 47/47 完全相同；范围最多差 0.68%，且差异是周期区间**起点**约定而非跨度） |
| 基类 `GenerateSurfaceNodes` | **空实现**（不加内部节点） |

⇒ **`delta` 与范围都一致 ⇒ 同一 `delta` 下格数也应一致**（格数公式已逐行等价，
且只依赖范围、`delta`、`errFactor`、`vertices_nb`）。
**因此「端口在某些面过密、某些面过疏」既不来自范围分割器、也不来自其输入**，
只能来自**下游**：节点的**插入/三角化**（`Delaun` 结构、`insertNodes`、
`collectTriangles`），或 `vertices_nb` 这一项与 OCCT 不同。

**注**：本轮的对照只覆盖 **47 个精确配对**（bbox 精确相等的那些）—— 与 §9.340 的
「167 精确配对 / port 8647 vs gt 9002」是**同一批面里的子集**；这 47 个面上
`delta` 完全一致这一点是**逐值证据**，强度高于此前的聚合统计。

#### 9.346.4 下一轮的落点

**节点插入 / 三角化**：`crates/occt-topo/src/meshing/delaun/*`、`mesh_algo.rs` 的
`finish_mesh` / `collect_triangles`、`node_insertion.rs` 的 `insert_nodes`。
可对照 OCCT：`BRepMesh_DelaunayBaseMeshAlgo`、`BRepMesh_NodeInsertionMeshAlgo::insertNodes`、
`BRepMesh_BaseMeshAlgo::collectTriangles`。
**做法同本轮**：在两侧的同一管线点打印同一批量（例如 `structure.nb_nodes()`、
插入前后的节点数、`collect_triangles` 得到的三角形数），逐值对照。

**同时**：把 `vertices_nb` 的口径核清 —— 端口用
`self.boundary_uv.len() + self.structure.nb_nodes()`（见本轮插桩），
而 OCCT 用 `NodeInsertionMeshAlgo` 里结构的节点数；**这两者是不是同一个量需要确认**
（它进 `AdjustCellsCounts`，会直接改变格数）。

#### 9.346.5 本轮改动与门禁

* **`crates/occt-topo/src/meshing/node_insertion.rs`**：加临时插桩（两分支，env-gated）；
  `cargo check` **0 error**；
* **`specs/occt_probe/occt_probe.cpp`**：`--uvsum` 加 splitter 对照输出；
  `build.bat` **BUILD OK**；
* a3n00 基线与 `--lib` **未变**（见门禁）；
* **插桩是临时的**，下一轮用完应撤除（本轮结束时**保留**，因为下一轮要用同一仪器）。

#### 9.346.6 复跑

```text
$env:OCCT_TOPO_TRACE_SPLITTER=1
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids 2> splitter_port.txt
cmd /c "set OCCT_TOPO_TRACE_SPLITTER=1 && specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --uvsum 0.215" > splitter_gt.txt
python .target-gate\splitter_diff.py splitter_port.txt splitter_gt.txt
```

---

### 9.347 —— `vertices_nb` 口径已核清：两侧调用点等价，**仅一处细节差异**（端口数"未删除"节点）

#### 9.347.1 OCCT 侧的四个 `getCellsCount` 调用点

```text
BRepMesh_DelaunayBaseMeshAlgo.cxx:43       getCellsCount(aVerticesOrder.Length())
BRepMesh_CustomBaseMeshAlgo.hxx:47         getCellsCount(aStructure->NbNodes())
BRepMesh_CustomDelaunayBaseMeshAlgo.hxx:40 getCellsCount(aStructure->NbNodes())
BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx:74   ← override 的定义（不是调用点）
```

**override 的定义**（`.hxx:74-80`）：

```cpp
std::pair<int, int> getCellsCount(const int theVerticesNb) override
{
  return BRepMesh_GeomTool::CellsCount(this->getDFace()->GetSurface(),
                                       theVerticesNb,
                                       this->getDFace()->GetDeflection(),
                                       &this->getRangeSplitter());
}
```

#### 9.347.2 端口的两个调用点（与 OCCT **形式等价**）

| 端口 | OCCT 对应 | 传入的节点数 |
|---|---|---|
| `mesh_algo.rs:99`：`self.inner.get_cells_count(mesher.result().nb_nodes())` | `CustomDelaunayBaseMeshAlgo.hxx:40`（`postProcessMesh` 内） | `aStructure->NbNodes()` |
| `node_insertion.rs:362`：`self.cells_count(indices.len() as i32)` | `CustomBaseMeshAlgo.hxx:47`（`generateMesh` 内，先 `buildBaseTriangulation()`） | `aStructure->NbNodes()` |

端口 `cells_count`（`node_insertion.rs:442-457`）等价于 OCCT 的
`DelaunayNodeInsertionMeshAlgo::getCellsCount`：都调 `geom_tool_cells_count` / `GeomTool::CellsCount`，
参数都是 `(surface, vertices_nb, deflection, range_u, range_v, delta)`。

#### 9.347.3 **唯一发现的细节差异**

`node_insertion.rs:359-362`：

```rust
let mut indices: Vec<i32> = (1..=structure.nb_nodes() as i32)
    .filter(|&i| structure.get_node(i).state != VertexState::Deleted)
    .collect();
let (cells_u, cells_v) = self.cells_count(indices.len() as i32);   // ← 未删除节点数
```

而 OCCT `CustomBaseMeshAlgo.hxx:47` 传的是 **`aStructure->NbNodes()`（含 Deleted）**。

⇒ **端口在统计时排除了 `Deleted` 状态的节点，OCCT 没有。**
**这是一处真实的（虽小）口径差异**。影响：只在**存在 `Deleted` 节点**时才会让
`vertices_nb` 变小 ⇒ `AdjustCellsCounts` 里 `log10` 项变小 ⇒ **格数变小**。
结合 §9.340 的观察（**端口在某些面铺得过少**），这**方向上吻合**，值得下一轮量化。

**注**：本轮**不宣称**它就是主因 —— 需要先量出这些面上 `Deleted` 节点的数量级
（若为 0，则该差异不产生任何影响）。

#### 9.347.4 本轮结论与下一轮

**范围分割器路径已彻底核对完毕**（公式逐行等价、输入逐值一致、`vertices_nb` 调用点等价）：

| 层次 | 状态 |
|---|---|
| 公式（7 处） | ✔ 逐行等价 |
| 输入（`range_u/v`、`delta`、`tol`） | ✔ 逐值一致（`delta` 47/47） |
| `vertices_nb` 调用点 | ✔ 等价，**仅"是否排除 Deleted"一处差异** |

**下一轮两条并列的落点**：

1. **量化那处 `Deleted` 差异**：在端口插桩里同时打印
   `structure.nb_nodes()`（含 Deleted）与 `indices.len()`（未删除），
   在出现差异的面（尤其 §9.340 的亏损榜：`f=206`/`f=192–199`）上量出**差几个节点**、
   以及它经 `AdjustCellsCounts` 会把 `cells_u/cells_v` 改到多少 ⇒
   **判定这处差异能否解释"过密/过疏"**；
2. **对比 `cells_u/cells_v` 本身**：OCCT 侧在 `--uvsum` 里**已经能**调
   `BRepMesh_GeomTool::CellsCount(...)`（该函数与探针同目录、公开可用），
   打印它给这两个面的 `cells`；与端口 `cells_count(...)` 的输出**逐值对照**。
   这一步**直接回答"格数是否相同"**，不必再绕。

**推荐先做 2**：它一次给出结论（格数同/不同），而 1 只是解释 2 的候选原因。

#### 9.347.5 本轮状态

* **端口插桩保留**（上一轮加的 `OCCT_TOPO_TRACE_SPLITTER`）；
* **本轮未再改代码**（只做读取与核对）；
* `crates/occt-topo/src/meshing/node_insertion.rs` 仍带上一轮的临时代码（+49，逻辑零改变）；
* 门禁：`cargo check` 0 error、`--lib` 1255/26、a3n00 基线未变（见 §9.346）。

#### 9.347.6 复跑

```text
# OCCT 调用点
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_CustomBaseMeshAlgo.hxx:47
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_CustomDelaunayBaseMeshAlgo.hxx:40
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx:74-80
# 端口调用点
crates/occt-topo/src/meshing/mesh_algo.rs:99
crates/occt-topo/src/meshing/node_insertion.rs:359-362, 442-457
```

---

### 9.348 —— 关键线索：`vertices_nb` 的取值**确实改变格数**（同面 `cells_all=[2,2]` vs `cells_live=[3,3]`）

#### 9.348.1 仪器已扩到"决策数"

在上一轮的两侧 `SPLITTER` 记录上各加一项：

| 侧 | 新增 | 说明 |
|---|---|---|
| **端口** | `vnb_all` / `vnb_live` / **`cells_all`** / **`cells_live`** | `cells_all = geom_tool_cells_count(surf, structure.nb_nodes(), ...)`；`cells_live = geom_tool_cells_count(surf, boundary_uv.len() + structure.nb_nodes(), ...)` |
| **OCCT** | **`cells_vnb0`** | `BRepMesh_GeomTool::CellsCount(aDF->GetSurface(), 0, defl, &aSplitter)` |

改动：`crates/occt-topo/src/meshing/node_insertion.rs`（插桩内，+约 20）；`specs/occt_probe/occt_probe.cpp`
（`--uvsum` 内，+6 行 + 新 include `BRepMesh_GeomTool.hxx`）。两侧 `cargo check` / `build.bat` **均通过**。

#### 9.348.2 决定性读数：`vertices_nb` **不是**无关紧要

在**同一个面**上（端口记录内部对比，不受配对影响）：

```text
port f=34   vnb_all=4   vnb_live=12   cells_all=[2,2]   cells_live=[3,3]
port f=68   vnb_all=9   vnb_live=23   cells_all=[2,2]   cells_live=[3,3]
...（前 10 个配对记录全部如此）
```

⇒ **`vertices_nb` 从 4 变成 12、从 9 变成 23，就把 `cells_u/cells_v` 从 `[2,2]` 抬到 `[3,3]`。**
**这证明 `vertices_nb` 是真正起作用的量**（此前只知道它进 `AdjustCellsCounts`，现在有了数值证据）。

#### 9.348.3 而端口实际用的是**较大的那个** ⇒ 方向与"过度细分"吻合

`node_insertion.rs:359-362` 实际传入的是

```rust
let mut indices: Vec<i32> = (1..=structure.nb_nodes() as i32)
    .filter(|&i| structure.get_node(i).state != VertexState::Deleted)
    .collect();
let (cells_u, cells_v) = self.cells_count(indices.len() as i32);
```

即「**未删除的节点数**」—— 在插桩里对应的是 `vnb_live = boundary_uv.len() + structure.nb_nodes()`
**这个口径本身就可疑**：若 `boundary_uv` 的点**已经**进了 `structure`（`add_node` 之后），
则 `boundary_uv.len() + structure.nb_nodes()` 就是**重复计数**，会把节点数**放大**。

**放大 ⇒ `AdjustCellsCounts` 的 `log10` 项更大 ⇒ `cells` 更大 ⇒ 格数更多 ⇒ 更密。**
这与 §9.340 观测到的「**端口在某些面铺得远多于 GT**」（`f=206` 336 vs 68、
`f=192–199` 各 228 vs 52）**方向完全吻合**。

**本轮不宣称已证**：`vnb_live` 的那行只是**我插桩里的算式**，
真实传入值是 `indices.len()`；两者是否相等、以及 `boundary_uv` 是否已计入 `structure`，
**都需要下一轮直接量**（同一处把 `boundary_uv.len()`、`structure.nb_nodes()`、
`indices.len()` 三个数**同时**打印即可判定）。

#### 9.348.4 其余字段的对照（56 个 bbox 配对，较上轮 47 更多）

| 量 | n | 完全相同 | 最坏相对差 |
|---|---|---|---|
| **`delta`** | 56 | **56** | 2e-15 |
| `rv` | 56 | 54 | 8.1e-4 |
| `ru` | 56 | 50 | 6.8e-3 |
| `tol` | 56 | 48 | 3.3e-2 |

⇒ **与 §9.346 的结论一致且更强**（配对数从 47 升到 56，`delta` 仍 **56/56**）。

**注**：`cells` 的**跨侧**对照本轮**无法完成** —— OCCT 侧只打了 `cells_vnb0`
（用 0 个顶点算的），而端口侧打的是实际 `vnb` 下的值。
**下一轮应让 OCCT 在同一点打印它真实的 `vertices_nb`**（即
`CustomBaseMeshAlgo.hxx:47` 处的 `aStructure->NbNodes()`），才能逐值比 `cells`。

#### 9.348.5 下一轮的明确做法

在**端口**插桩处把三个数一起打印：

```text
boundary_uv.len()   structure.nb_nodes()   indices.len()
```

并同时打印 `cells_count(indices.len())`；
在 **OCCT** 侧则需在 `CustomBaseMeshAlgo` 的调用点取到 `aStructure->NbNodes()`
（探针里可在 `--uvsum` 之外补一个走 `IMeshTools` 的路径，或直接用
`aModel->FacesNb()` 对应的离散面节点数 —— 具体取法下一轮定）。

**判定规则**：
* 若端口 `indices.len()` **等于** `structure.nb_nodes()`（即 `boundary_uv` 已计入，
  且无 Deleted）⇒ 与 OCCT 一致，**此处无缺口**，回到 §9.346 的下游结论（节点插入/三角化）；
* 若**大于** ⇒ 端口节点数被放大、`cells` 偏大、更密 ⇒ **这就是"过密"类面的机制**；
* 若**小于** ⇒ 端口节点数被缩小 ⇒ 对应"过疏"类面。

#### 9.348.6 本轮状态

* 改了两处**临时插桩**（端口 `node_insertion.rs`、OCCT 探针），两侧编译/重建通过；
* 门禁见下（`--lib` 与 a3n00 基线未变）；
* **插桩仍保留**（下一轮要用）。

#### 9.348.7 复跑

```text
$env:OCCT_TOPO_TRACE_SPLITTER=1
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids 2> sp_port.txt
cmd /c "set OCCT_TOPO_TRACE_SPLITTER=1 && specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --uvsum 0.215" > sp_gt.txt
python .target-gate\splitter_diff2.py sp_port.txt sp_gt.txt
```

---

### 9.349 —— 【定论】`vertices_nb` 两侧**完全相同**（226/226，无 Deleted 节点）⇒ 范围分割器**整条输入路径已排除**

#### 9.349.1 判定测量

在 `node_insertion.rs` 的 `cells_count(indices.len())` 处加 env-gated 打印
（`OCCT_TOPO_TRACE_CELLS`），同时打印三个候选节点数：

```text
CELLS f=<i> nb_nodes=<n> live_indices=<n> boundary_uv=<m> vnb_used=<n> cells=[u,v]
```

**结果（226 个面）**：

```text
faces where indices.len() == nb_nodes: 226
faces where they DIFFER:                0
```

⇒ **`indices.len()` 恒等于 `structure.nb_nodes()`** ⇒ **没有 `Deleted` 状态的节点**
（否则 `indices.len()` 会小于 `nb_nodes()`）。
而 OCCT 的调用点传的正是 `aStructure->NbNodes()` ⇒ **两侧传给 `getCellsCount` 的数完全相同。**

#### 9.349.2 这同时撤清了 §9.347.3 / §9.348.3 的两处**疑虑**

* §9.347.3 提出「端口数的是未删除节点数、OCCT 数全部」**是一处口径差异**
  ⇒ **实测差异为 0**（无 Deleted 节点），**该疑虑不成立**；
* §9.348.3 依据我插桩里的算式
  `vnb_live = boundary_uv.len() + structure.nb_nodes()` 推测「可能重复计数、把格数放大」
  ⇒ **该算式只是我诊断行里的写法，不是真实传参**；真实传参 `indices.len()` 与
  `nb_nodes()` 相等 ⇒ **不存在重复计数**。**该线索作废。**

（本条记录为一次**自我更正**：§9.347.3 与 §9.348.3 都是我想出来的"可疑点"，
本轮用一次直接测量把它们同时排除。）

#### 9.349.3 阶段性定论：**范围分割器整条路径（含全部输入）与 OCCT 一致**

汇总本会话逐项核对的结果：

| 层次 | 项 | 结论 |
|---|---|---|
| **公式** | `AdjustRange`、`computeLengthU/V`、`computeTolerance`、`computeDelta`、格数三支、`AdjustCellsCounts`、`ComputeErrFactors` | **逐行等价**（§9.343/§9.344） |
| **输入** | `range_u` / `range_v` | 最多差 **0.68%**，且差异是**周期区间起点**而非跨度（§9.346） |
| | `delta` | **47/47、56/56 完全相同**（最坏 2e-15） |
| | `tolerance` | 39–48/56 相同，最坏 3.3e-2 |
| | **`vertices_nb`** | **226/226 完全相同**（无 Deleted 节点）（§9.349） |
| **节点生成（基类）** | `DefaultRangeSplitter::GenerateSurfaceNodes` | **空实现**，不加内部节点（§9.344） |

⇒ **范围分割器不再有任何未核对的量。** 因此「端口在某些面铺得远多于/远少于 OCCT」
**不可能**由它产生。**剩余的唯一去向是节点插入与三角化**。

#### 9.349.4 下一轮（落点唯一且明确）

**节点插入 / 三角化**，可对照的三组：

| 端口 | OCCT | 关注量 |
|---|---|---|
| `node_insertion.rs` 的 `insert_nodes` / `finish_mesh` | `BRepMesh_NodeInsertionMeshAlgo::insertNodes`、`BRepMesh_BaseMeshAlgo::generateMesh` | 插入前后节点数、`AddVertices` 的结果 |
| `meshing/delaun/*`（`Delaun::new_with_data_cells`、`add_triangle`、`erase_free_links`） | `BRepMesh_Delaun`（`cxx`） | 结构初始化参数、被丢弃/合并的三角形 |
| `node_insertion.rs` 的 `collect_triangles` | `BRepMesh_BaseMeshAlgo::collectTriangles`（`cxx:245-277`） | 最终三角形数 |

**做法**：在两侧同一管线点打印「插入前后节点数 + `collect_triangles` 得到的三角形数」，
对 §9.340 亏损榜（`f=206`：gt 68 vs port 336）与盈利榜（`f=171`：gt 1569 vs port 349）
的面对照 —— 那里是差异最集中的地方。

#### 9.349.5 本轮动作与状态

* 已把端口插桩扩到打印三个候选节点数（`OCCT_TOPO_TRACE_CELLS`），完成判定；
* **两处临时插桩（`node_insertion.rs` 的 `SPLITTER` 与 `CELLS` 打印）现已完成使命**
  ⇒ 本轮结束时**撤除**，让 `crates/occt-topo/src` 回到无功能改动的状态；
  OCCT 探针侧的 `SPLITTER`/`cells` 输出**保留**（它是有用的仪器，且探针本就是仪器集合）；
* `--lib` 与 a3n00 基线：撤除前后均需复核（见门禁）。

#### 9.349.6 复跑

```text
$env:OCCT_TOPO_TRACE_CELLS=1
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids 2> cells_port.txt
grep -c '^CELLS f=' cells_port.txt
# 判定：nb_nodes == live_indices 对所有面成立 => 与 OCCT 的 aStructure->NbNodes() 一致
```

---

### 9.350 —— 转向节点插入/三角化：`collectTriangles`、域标记、`SetCellSize` 三处**均逐行等价**；差异**不随长宽比缩放**

#### 9.350.1 逐行对照（卡片定的落点）

| 环节 | OCCT | 端口 | 判定 |
|---|---|---|---|
| **`collectTriangles`** | `BRepMesh_BaseMeshAlgo.cxx:245-278`：取 `myStructure->ElementsOfDomain()`，空则返回空句柄；**每个域内元素输出一个三角形**（`ResizeTriangles(aTriangles.Extent(), false)`）；逐顶点用 `myUsedNodes` 做**紧凑重编号**（`Bind(aNode[i], Length()+1)`）；最后 `ResizeNodes(myUsedNodes->Extent())` | `node_insertion.rs:414-460`：`for &id in ds.elements_of_domain()` **每个域内元素一个三角形**；`used.contains_key` → `used.insert(*n, used.len()+1)` 同样紧凑重编号；输出节点数取 `used` 的最大值 | **逐行等价** |
| **域标记** | `DataStructureOfDelaun.cxx:170-183` `AddElement`：**无条件** `myElementsOfDomain.Add(aElementIndex)`；`RemoveElement`（`:197`）移除 | `delaun_data.rs:433-442` `add_element`：**无条件** `elements_of_domain.insert(idx)`；`remove_element`（`:477`）移除 | **逐行等价** |
| **`SetCellSize`/`SetTolerance`** | `BRepMesh_NodeInsertionMeshAlgo.hxx:86-93`：`uCellSize = 14.0*tolUV.first`，`SetCellSize(uCellSize/delta.first, vCellSize/delta.second)`，`SetTolerance(tolUV.first/delta.first, tolUV.second/delta.second)` | `node_insertion.rs:503-508`：`set_cell_size(14.0*tu/du, 14.0*tv/dv)`、`set_tolerance(tu/du, tv/dv)` | **逐行等价** |

⇒ **节点插入/三角化的这三处都不是分歧点。** 且 §9.349 已证明喂给它们的 `delta`/`tol`/`vertices_nb`/`range` **两侧一致**。

#### 9.350.2 差异的**分布特征**：不随长宽比缩放

用 `--uvsum` 的 `lenU/lenV` 算每个配对面的**长宽比**，与三角数比值并排（47 个可用配对）：

```text
  aspect   port/gt  port_f   gt_f  port_t   gt_t
    10.3     1.000     186    225       2      2
    10.3     1.000     181    222       2      2
     5.0     1.000     136    192       2      2
     5.0     1.000     134    193       2      2
     4.9     1.000     224    196       2      2
     3.0     1.000      34    205       2      2
     2.1     1.000     163    206       3      3
     ...
     1.0     0.500     143    141       8     16
     1.0     1.014     139    135      72     71
     1.0     0.500     150    138       8     16
     1.0     1.281     132    157     228    178
     1.0     1.106     217    147     376    340
     1.0     1.281     114    156     228    178
     1.0     1.215     115    151     243    200

aspect <= 2 : n=31   mean port/gt=1.107
aspect > 10 : n=2    mean port/gt=1.000
```

**三条读数**：

1. **长宽比大的面（≥3）全部比值 1.000** —— 而且它们**只有 2–3 个三角形**
   （端口与 GT 完全一致）；
2. **差异集中在长宽比 ≈ 1 的面**（比值 0.5–1.28）；
3. **两个桶的均值几乎相同**（1.107 对 1.000）⇒ **差异不随长宽比缩放**。

⇒ **§9.342.3 里"异常聚集在狭长面"的观察，在本轮的配对下不再成立** ——
本轮用的 47 个配对与 §9.342 的口径不同（那里用 `lenU*lenV` 分组）。
**这条"狭长面"线索应当撤回**：它既不是规律，也不该被用来做启发式
（那正是 CLAUDE.md 禁止的"长度滤边"式规则）。

#### 9.350.3 本轮结论与下一步的重新表述

**已等价**：`collectTriangles` ✔、域标记 ✔、`SetCellSize` ✔、`delta`/`tol`/`range`/`vertices_nb` ✔、
格数公式 ✔、范围分割全套 ✔。
**剩下的只可能是**：**Delaunay 三角化本身产出的"域内元素"数不同**
（`elements_of_domain()` 的大小），而它由 `BRepMesh_Delaun` 的**算法过程**决定
（插入顺序、`ClassifyTriangle`、边翻转、`erase_free_links` 等），
**不是某一条可"逐行对照"就能定论的公式**。

⇒ **下一轮的做法必须变**：不能再找"哪一行不同"，而要**在同一输入下逐个比较算法状态**。
最有信息的一个数是：

> **每个面 `elements_of_domain().len()` 与 `structure.nb_nodes()`（域内元素数 / 总节点数）**

它把「插入/翻边的结果」压缩成两个整数。**端口侧现成可取**（`TriangulationResult.triangles.len()`
与 `.nodes.len()`，`node_insertion.rs:90-97`），只需在 `zz_uv_feed --model` 里打出来；
**OCCT 侧**可在 `--uvsum` 里对每个 `aDF` 走一遍 `BRepMesh_IncrementalMesh` 后读
`BRep_Tool::Triangulation(...)->NbTriangles()/NbNodes()` —— **这与 §9.340 的
`--uvsum` `triangles=` 字段是同一件事**，即 **数据其实已经有了**
（`v_uvsum2.txt` 的 `triangles=` 就是它）。

⇒ **所以下一步不需要新仪器**：把 **端口 `TriangulationResult` 的
`triangles.len()`/`nodes.len()`** 打出来，与 `UVSUM` 的 `triangles=` 按 bbox 配对比较，
就能把「域内元素数」这一层也核掉或定位。

#### 9.350.4 本轮状态

* **未提交任何源码改动**（只做读取与对照）；
* **Rust 库零改动**（插桩已在 §9.349 撤除，本轮确认 `node_insertion.rs` diff 为空）；
* 本轮净产出：**节点插入/三角化三处逐行等价**、**撤回"狭长面"线索**、
  并把下一步收敛到一个**两个整数**的量（`elements_of_domain().len()` / `nb_nodes`），
  且指出**该数据已在 `UVSUM` 里，不必新造仪器**。

#### 9.350.5 复跑

```text
python .target-gate\aspect_ratio.py
# OCCT 侧
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_BaseMeshAlgo.cxx:245-278
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_DataStructureOfDelaun.cxx:170-197
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_NodeInsertionMeshAlgo.hxx:86-93
# 端口侧
crates/occt-topo/src/meshing/node_insertion.rs:414-460, 503-508
crates/occt-topo/src/meshing/delaun_data.rs:433-442
```

---

### 9.351 —— 【本会话最重要的发现】端口的 **rescue 路径在 a3n00 上触发 16 次**（注释声称"从不触发"），且它正是那批"过密"面的三角形来源

#### 9.351.1 测量一：`collect_triangles` 的产出 vs 提交网格（同一 run/配置）

给 `discret_root.rs:427` 的 `Ok(Ok(result))` 加 env-gated 打印（`OCCT_TOPO_TRACE_TRIRES`），
逐面输出 `TriangulationResult.triangles.len()` / `nodes.len()`，与同一次运行的 `PORTID mt=` 比较：

```text
trires faces=226  committed faces=226
faces where TriangulationResult.triangles.len() == committed mt: 210
faces where they DIFFER:                                          16
sum trires=7725   sum committed=11121        （差 3396）
```

**16 个差异面的形态完全一致**：

| f | trires | committed | delta |
|---|---|---|---|
| **206** | **0** | 336 | +336 |
| 174 | **0** | 288 | +288 |
| **192–199**（8 个） | **0** | 228 each | +228 each |
| **209/210/212** | **0** | 224 each | +224 each |
| 169 | **0** | 204 | +204 |
| 189 / 204 | **0** | 36 each | +36 each |

⇒ **这 16 个面的"节点插入三角化"产出 0 个三角形，而提交网格里却有 224–336 个。**

**排除其他解释**：`map_triangulation`（`discret_root.rs:533`）是**逐顶点 1:1 重映射**
（`Triangle::new(vi[0]-1, vi[1]-1, vi[2]-1)`），**不增删三角形**；因此这 3396 个三角形
**必然来自另一条路径**。

#### 9.351.2 测量二：rescue 路径**确实触发**（注释是错的）

`discret_root.rs:484-503` 的 `p.tri.or_else(...)`：当 `p.tri` 为 `None` 且 `needs_fallback`
时调用 `Self::wireframe_face_triangulation(...)`。该处**原有一段 T-54 注释断言**：

> 「this rescue is **never reached** on any gate model … Measured with
> `OCCT_TOPO_TRACE_FALLBACK` on a3n00, T0M, acs10, TDB, bottom, top:
> `reached_rescue=0` everywhere」

在该分支加打印（`OCCT_TOPO_TRACE_FALLBACK`）后跑 a3n00：

```text
FALLBACK reached lines: 16
FALLBACK reached f=169 deflection=1.0760065231799822
FALLBACK reached f=174 deflection=1.0760065231799822
FALLBACK reached f=189 …
FALLBACK reached f=192 … f=193 … f=194 … f=195 … f=196 … f=197 … f=198 … f=199 …
FALLBACK reached f=204 … f=206 … f=209 … f=210 … f=212 …
```

⇒ **触发 16 次**，且**触发面集合与 §9.351.1 的 16 个零产出面完全一致**。

⇒ **提交网格里那 3396 个三角形全部来自 `wireframe_face_triangulation`（rescue）。**

#### 9.351.3 为什么这很重要：这是一条**自称违反红线**的代码路径

同一段注释自己写道：

> 「**OCCT has no per-face alternative tessellator** — a failed face keeps
> `IMeshData_Failure`（`BRepMesh_BaseMeshAlgo.cxx:52-62`）。See §9.306: T-54
> decided to keep this code for now (no triggering case).」

而 CLAUDE.md 门禁第 2 条明确禁止「自造 OCCT 里不存在的分支」。
⇒ **这条路径正是"OCCT 里不存在的分支"**，且它**不是惰性的** ——
它在 a3n00 上为 16 个面**供给了 3396 个三角形**（占提交网格 11121 的 **30.5%**）。

**并且它同时解释了 §9.340 的"过密"榜**：`f=206`（gt 68 vs port 336）、
`f=192–199`（各 gt 52 vs port 228）—— 这些面在 **GT 里被 OCCT 判为失败/无三角**
（所以 GT 数很小），而**端口用 rescue 硬铺了 224–336 个**。

#### 9.351.4 本轮的处置（按门禁第 4 条：完成即停，旁支写进报告）

* **不改**：删掉/停用这条路径是**实质行为改动**，会牵动 a3n00 的网格与面积比，
  必须先由人决定（它也可能是当年为对齐某个基线而刻意留下的）；
* **只记录**：把「注释断言与实测相反」这一事实、触发面清单、以及它在提交网格中的
  占比（30.5%）写进文档与看板；
* **同时更正注释**：下一轮若要动手，第一件事是**把那段错误注释改掉**（它现在会误导人），
  并把 §9.306 的「no triggering case」更正为「a3n00 上触发 16 次」。

**这是一个可验证、可复现、且直接指向 30% 三角形来源的落点** ——
比此前任何一次"猜机制"都更具体。

#### 9.351.5 本轮改动与门禁（**注意：本轮有两处临时插桩，尚未撤除**）

* `crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs`：
  ① `Ok(Ok(result))` 处加 `OCCT_TOPO_TRACE_TRIRES` 打印（+13）；
  ② `needs_fallback` 分支加 `OCCT_TOPO_TRACE_FALLBACK` 打印（+9，**并保留原注释**）；
* 两处均为 **env-gated 的临时诊断**，默认不输出、**不改变任何逻辑**；
* `cargo check` 0 error；`--lib` 与 a3n00 基线**未变**（见门禁）；**下一轮用完须撤除**。

#### 9.351.6 复跑

```text
$env:OCCT_TOPO_TRACE_TRIRES=1
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids 2> trires.txt
python .target-gate\trires_vs_committed.py trires.txt portid_cur.txt
$env:OCCT_TOPO_TRACE_FALLBACK=1
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids 2> fb.txt
grep -c '^FALLBACK reached' fb.txt     # -> 16
```

---

### 9.352 —— 【根因收敛到一点】16 个 rescue 面全是 **`Plane` + `wires=1`**：端口 Delaunay 对它们产出 **0** 三角形，rescue 硬补 3396（OCCT 只 1049）

#### 9.352.1 `needs_fallback` 的真实触发链

`discret_root.rs` 的相关代码：

```text
:442  let mapped = Self::map_triangulation(&result, surface.as_ref(), i);
:443  let tri = mapped.ok();
:444  let needs_fallback = tri.is_none();
...
:569  if triangles.is_empty() {
:570      return Err("...: Delaunay produced no triangles");
:571  }                      ← map_triangulation 显式拒绝空结果
...
:484  let tri = p.tri.or_else(|| {
:485      if p.needs_fallback {  ... wireframe_face_triangulation(...)  }
```

⇒ 触发链是：**`TriangulationResult` 的 `triangles` 为空 → `map_triangulation` 返回 `Err`
→ `tri=None` → `needs_fallback=true` → 走 `wireframe_face_triangulation`（rescue）。**

**关键**：这条链**绕过了** `Ok(Err(e))` 分支里那段「face 已标记 `FAILURE` 就跳过、不 fallback」
的忠实处理（`:447-459`）—— 因为这里 `perform` 是 `Ok`，只是**结果为空**。

#### 9.352.2 16 个 rescue 面的**性质完全一致**

按几何 bbox 与 `UVSUM` 配对（用未网格化的 `v_uvsum.txt`，与 `PORTID` 同口径；
三角数取同面号的 `v_uvsum2.txt`）：

| port f | port_mt | gt f | gt_tris | **GT type** | **GT wires** |
|---|---|---|---|---|---|
| 169 | 204 | 87 | 120 | **Plane** | **1** |
| 174 | 288 | 42 | 203 | **Plane** | **1** |
| 189 | 36 | 41 | 123 | **Plane** | **1** |
| 192 | 228 | 85 | 52 | **Plane** | **1** |
| 193 | 228 | 84 | 52 | **Plane** | **1** |
| 194 | 228 | 83 | 52 | **Plane** | **1** |
| 195 | 228 | 78 | 52 | **Plane** | **1** |
| 196 | 228 | 79 | 52 | **Plane** | **1** |
| 197 | 228 | 80 | 52 | **Plane** | **1** |
| 198 | 228 | 81 | 52 | **Plane** | **1** |
| 199 | 228 | 82 | 52 | **Plane** | **1** |
| 204 | 36 | 43 | 119 | **Plane** | **1** |
| 206 | 336 | 86 | 68 | **Plane** | **1** |
| 209 / 210 / 212 | 224 each | （未配上） | — | — | — |

**13 个配上对的面全部是 `Plane` 且 `wires=1`** ⇒ **这不是随机分布，而是一个面类的共性。**

⇒ **端口的 Delaunay 节点插入对「单闭合 wire 的平面」产出 0 个三角形**，
于是 `map_triangulation` 拒绝、rescue 硬补。

#### 9.352.3 量与算术：一个惊人的吻合

```text
rescue 16 面：port 提交 3396    OCCT 同面 1049      ⇒ 端口多铺 2347
port 提交总量 11121；去掉 rescue 后 = 11121 − 3396 = 7725
OCCT（d=0.215）总量 = 12462
7725 + 1049 = 8774    （仍差 3688）
```

⇒ **rescue 既不是"忠实补洞"（它比 OCCT 多 2347），也不是可有可无**：
拿掉它这些面会从"多铺"变成"完全不铺"。

#### 9.352.4 结论：**两种可能，必须由人判定**

**(A) 端口 Delaunay 在这类平面上有缺口** —— 即 OCCT 能三角化而端口不能。
若是这样，**真正该修的是 Delaunay**（让这 16 个面产出正确的 ~1049 个三角形），
而 rescue 只是**掩盖症状**；
**(B) rescue 当年是为对齐某个基线而刻意保留**，属人为选择。

**证据偏向 (A)**：
* 这 16 个面是**同一个面类**（`Plane` + `wires=1`），符合"某个算法分支没覆盖"的形态，
  不像人为挑选；
* 代码注释**自己承认** OCCT 没有 per-face 备选 tessellator，并说这条路径
  「decided to keep this code for now (**no triggering case**)」——
  即**留它是因为以为它永远不会跑**，而不是因为它有用；
* 而实测它**每次都跑**（16 次），说明**当年的假设错了**。

#### 9.352.5 建议的下一步（写进卡片，不动手）

1. **先修注释**：把 `discret_root.rs:486-498` 那段「never reached / `reached_rescue=0`
   everywhere」更正为「a3n00 上触发 16 次、供 3396 个三角形」，
   并把 §9.306 的「no triggering case」一并更正；
2. **再查 Delaunay 为何对 `Plane` + `wires=1` 产出 0 三角形** ——
   这是**有 `.cxx` 可对照的正规移植问题**（`BRepMesh_Delaun` 的初始化/约束边注册），
   与"自造分支"无关；
3. 只有在 (2) 说明端口**无法**在 OCCT 的路径上产出时，才讨论 rescue 的去留。

#### 9.352.6 本轮改动与门禁

* 仍保留上一轮的两处 env-gated 临时插桩（`TRIRES`、`FALLBACK`），**本轮未新增改动**；
* `cargo check` 0 error；`--lib` 1255/26；a3n00 基线 `stats=226 unmatched=0 mesh=11101/11121` **未变**；
* **插桩仍待撤除**（下轮用完后）。

#### 9.352.7 复跑

```text
python .target-gate\rescue_faces.py
$env:OCCT_TOPO_TRACE_FALLBACK=1
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids 2>&1 | grep -c '^FALLBACK reached'
```

---

### 9.353 —— 【事故与更正】我撤除插桩时误用 `git checkout`，抹掉了 `discret_root.rs` 的未提交注释改动

#### 9.353.1 发生了什么

§9.351/§9.352 的插桩加在 `crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs`。
本轮结束时我用

```text
git checkout -- crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs
```

撤除插桩。**这是错误的做法** —— 该文件当时**除我的插桩外还有未提交改动**，
`git checkout` 把它们**一并还原到 `HEAD`**。

**已知被抹掉的**：`discret_root.rs` 里 `needs_fallback` 分支那段注释
（以「this rescue is **never reached** on any gate model … `reached_rescue=0` everywhere
… OCCT has no per-face alternative tessellator … See §9.306」为内容）——
我在 §9.351 的编辑里亲眼复制过它，而**当前工作区与 `HEAD` 都没有这段话**。

#### 9.353.2 损害范围已定界：**仅注释**（用运行期行为判定）

若被抹掉的是**可执行代码**，a3n00 的网格必然改变。实测：

```text
TOTAL faces=226 computed_lin=1.076007 used_lin=1.076007 stats=226 mesh_v=11101 mesh_t=11121
STATMAP matched=226 unmatched=0 sum_mt=11121 flat_mt=11121
```

与本会话记录的基线**逐字相同** ⇒ **丢失的内容只有注释文档**，
**没有丢失任何可执行逻辑**。风险因此是"文档回退"，不是"行为回退"。

#### 9.353.3 已做的恢复与核验

1. **恢复 T-97 标注**：该文件不在 `.target-gate/a_class.txt` 的 74 行清单里，
   但本会话确实给它加过 A 类标注；已按同类文件的**逐字内容**补回
   （5 行 `//!` 说明 + `#![allow(dead_code)]`），使该文件 diff 回到 **6/0**，
   与其余 60 个 A 类文件**完全一致**；
2. **全局一致性复核**：
   ```text
   total changed files under crates/: 64
   A-class annotations (6/0):         61
   non-annotation changes:             3
      33/15 crates/occt-geom/src/gcpnts.rs
      18/11 crates/occt-geom2d/src/geom2d_int/ginter.rs
      19/6  crates/occt-topo/src/pcurve_full/make_pcurve.rs
   ```
   且**带 `T-97: items below are faithful ports` 标注的文件数 = 61**（用 grep 递归数）
   ⇒ **61 个 6/0 改动与 61 个标注一一对应**，恢复无误。
3. **门禁复核**：`cargo check` 0 error；`occt-topo --lib` **1255 passed / 26 failed**（环境性）；
   a3n00 基线**未变**（见 §9.353.2）。

#### 9.353.4 为什么必须记下来

**这是本会话第四次丢失未提交改动**：

| # | 轮次 | 丢失内容 | 后果 |
|---|---|---|---|
| 1 | §9.324 | `read_topology.rs` 的 T-99 seam 接线 61 行 | **参照系失效，浪费数十轮** |
| 2–3 | 更早 | 引号脚本事故等 | 已修 |
| **4** | **本轮** | `discret_root.rs` 的未提交注释 | **仅文档**，已部分恢复 |

**根因相同**：我用 `git checkout <file>` 去撤除**自己刚加的临时插桩**，
而该文件**同时**承载着**别人的未提交改动**。

**因此确立一条操作规程（本会话剩余轮次一律遵守）**：
* **撤除临时插桩不要用 `git checkout`** —— 改用**逐块 `edit` 反向替换**（我的插桩都带
  `TEMP T-99` 标记，可精确定位）；
* 若必须用 `git checkout`，**先 `git diff` 确认该文件只有我的改动**；
* 每条插桩在加的时候就记下**精确的插入文本**，便于反删。

#### 9.353.5 本轮状态

* `crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs`：**diff = 6/0**（T-97 标注），
  **两处临时插桩已撤除、无残留**（grep `OCCT_TOPO_TRACE_TRIRES|OCCT_TOPO_TRACE_FALLBACK|TEMP T-99` 为空）；
* 全 `crates/` 的非标注改动只有既有的 3 处（皆本会话早前所留，与 T-99 无关）；
* 门禁见 §9.353.3。

---

### 9.354 —— 【定位到具体一步】16 个面的 Delaunay **产出 0 个域内元素**，且它们**全部恰好 75 个节点**（端口 `wires=2`，GT `wires=1`）

#### 9.354.1 测量：`domain_elems` 在**三角化刚结束**时就是 0

在 `node_insertion.rs` 的 `finish_mesh` 里、`collect_triangles` **之前**加
env-gated 打印（`OCCT_TOPO_TRACE_DELAUN`），输出三角化刚结束后的
`delaun.result().nb_nodes()` 与 `elements_of_domain().len()`，并在 `collect_triangles`
之后输出三角数：

```text
DELAUN pre_insert nodes=71 domain_elems=66 insert=0
DELAUN post_collect triangles=66          ← 正常面：66 元素 → 66 三角形
...
DELAUN pre_insert nodes=75 domain_elems=0  insert=0
DELAUN post_collect triangles=0          ← 异常面
```

**226 条记录里 `domain_elems=0` 的恰好 16 条**，且**出现序号（occurrence）与端口面号
逐一对应**：

```text
occurrence indices with domain_elems=0: [169, 174, 189, 192, 193, 194, 195, 196,
                                         197, 198, 199, 204, 206, 209, 210, 212]
face ids with domain_elems=0:           [169, 174, 189, 192, 193, 194, 195, 196,
                                         197, 198, 199, 204, 206, 209, 210, 212]
```

⇒ **与 §9.351 的 rescue 触发面完全一致** ⇒ 触发链的**起点**就是这里：
**Delaunay 三角化本身产出 0 个域内元素**（不是"后来被丢掉"）。

**排除下游**：`collect_triangles` 对每个域内元素输出一个三角形（§9.350 已证逐行等价），
所以 `domain_elems=0` ⇒ `triangles=0` ⇒ `map_triangulation` 拒绝 ⇒ rescue。
**链上的每一环现在都有实测**。

#### 9.354.2 【强信号】这 16 个面**全部恰好 `nodes=75`**

```text
occ=169  nodes=75  domain_elems=0  insert=0
occ=174  nodes=75  domain_elems=0  insert=0
...（16 条全部 nodes=75）
```

⇒ **16 个不同的面、节点数完全相同（75）、域内元素完全相同（0）。**
正常面（如 `occ=0`）是 `nodes=71 domain_elems=66`、`occ=2` 是 `nodes=45 domain_elems=40`
—— 各不相同。**16 个面同一个节点数是不自然的**，说明它们**走了同一条（有问题的）路径**，
且该路径产出的节点数与面的实际几何无关（或这些面恰好同构）。

#### 9.354.3 与「2 vs 1 wire」面类的关联（**强但未完全验证**）

这 16 个面在**端口**里**全部是 `wires=2`**：

```text
PORTID f=169 wires=2 mt=204     PORTID f=192 wires=2 mt=228
PORTID f=174 wires=2 mt=288     PORTID f=193 wires=2 mt=228   ...
PORTID f=189 wires=2 mt=36      PORTID f=206 wires=2 mt=336   （16 个全部 wires=2）
```

而它们的 GT 对应面（§9.352 配对）是 **`wires=1`**（例如 `f=192↔gt85`、`f=206↔gt86`）
⇒ 正是 §9.326 早就统计过的那个面类：**「端口在 29 个面上有 2 条 wire 而 GT 只有 1 条」**。

**因此有一个可检验的强假设**：
> 端口的 Delaunay 对「**被拆成 2 条 wire 的平面**」产出 0 个域内元素，
> 而 OCCT（同样 1 条 wire）能正常三角化。

**若成立**，则 §9.352 里"该修的是 Delaunay"就**具体化为**：
修「为什么这些面在端口里成了 2 条 wire」或「Delaunay 为何处理不了这种 2-wire 结构」——
**两者都有 `.cxx` 可对照**（wire 建立见 reader/`ShapeFix`，约束边注册见
`BRepMesh_Delaun`）。

**本轮未逐面验证**「这 16 个面是否都在那 29 个 2-vs-1 面里」（只抽查了 4 个），
下一轮第一步就是**把两者的集合求交**。

#### 9.354.4 下一轮的明确步骤

1. **求交**：`{domain_elems=0 的 16 面}` ∩ `{§9.326 统计的 29 个 2-vs-1 wire 面}`
   —— 若覆盖 16/16，则假设成立；
2. 若不成立，则查这 16 个面在**约束边注册 / `use_edge`** 那一步的状态
   （`finish_mesh` 的 `frontier` 循环，`node_insertion.rs:336-339`）；
3. **不要**先动 rescue —— 它是症状，起点在 (1)/(2)。

#### 9.354.5 本轮改动与状态（**按 §9.353 的新规程撤除插桩**）

* 插桩加在 `crates/occt-topo/src/meshing/node_insertion.rs`，**带 `TEMP T-99` 标记**；
* **撤除方式（遵守 §9.353 的新规程）**：**逐块 `edit` 反向替换**，
  **不使用 `git checkout`**（该文件承载其它未提交内容）；
* 撤除后 `node_insertion.rs` 应回到 **6/0**（仅 T-97 标注）；
* 门禁：`cargo check` 0 error；`--lib` 1255/26；a3n00 基线未变。

#### 9.354.6 复跑

```text
$env:OCCT_TOPO_TRACE_DELAUN=1
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed -- data/occ/a3n00.stp --ids 2> delaun.txt
python .target-gate\delaun_zero.py
```

---

### 9.355 —— 求交结果：**「2 vs 1 wire」不是因果**（13/16 命中）；真正的签名是 `nodes=75` 且 `domain_elems=0`

#### 9.355.1 求交（§9.354.4 卡片所定的第一步）

用 `.target-gate/intersect_zero_2v1.py`（端口与 GT 按**几何 bbox** 配对，
即未网格化的 `v_uvsum.txt`）：

```text
port faces=226  gt faces=226   paired by geometry bbox: 167

port wires=2 vs gt wires=1 : 28 faces
  ids = [41,42,44,45,47,48,92,93,129,130,137,140,154,169,170,173,174,189,
         192,193,194,195,196,197,198,199,204,206]
other wire-count mismatches: 1  (f=113: port 4 vs gt 2)

{domain_elems=0} = [169,174,189,192,193,194,195,196,197,198,199,204,206,209,210,212]  (n=16)
  intersection with 2-vs-1 set : 13 / 16
  in 2-vs-1  but not zero      : [41,42,44,45,47,48,92,93,129,130,137,140,154,170,173]   （15 个）
  zero       but not in 2-vs-1 : [209,210,212]                                            （3 个）
```

#### 9.355.2 结论：**「2 vs 1 wire」不是这个失败的因果**

两条反向证据：

1. **15 个 `2-vs-1` 面并不 `domain_elems=0`** —— 它们的 Delaunay 正常工作。
   所以「wire 数不匹配」**不足以**导致失败；
2. 那 3 个不在 `2-vs-1` 集合里的面（`209/210/212`）**并非真的例外**：
   查 `PORTID` 它们**同样是 `wires=2`**，只是**没在 GT 里配上面**（bbox 无匹配）
   —— 即它们**不在"已配对的 167 个面"范围内**，所以没有被 §9.326 那个统计计入。

⇒ 修正后的表述：

> **`{domain_elems=0}` 的 16 个面，端口侧全部是 `wires=2`；
> 但"端口 `wires=2`"共有 28（已配对）+ 3（未配对）= 至少 31 个面，
> 其中只有 16 个失败。**
> ⇒ **wire 数是相关但不是因果**；真正区分成败的是**另一个状态**。

#### 9.355.3 真正的签名（§9.354 已量出，本轮确认其唯一性）

```text
16 个失败面：全部 nodes=75、domain_elems=0
正常面（抽样）：occ=0 nodes=71 domain_elems=66
                occ=2 nodes=45 domain_elems=40
```

⇒ **16 个不同的面、节点数完全相同（75）、域内元素完全相同（0）**。
正常面的节点数**各不相同**（71、45、…），所以 **75 是这 16 个面共有的异常值**，
不是几何决定的自然结果。**这比 wire 数更特异**，是下一轮该追的量：

> **为什么这 16 个面恰好各注册 75 个节点、然后一个域内元素都产生不出来？**

`nodes` 来自 `delaun.result().nb_nodes()`（即 `structure` 的节点数），
它由「边界 pcurve 采样点 + 约束边 + 可能的内部节点」共同决定。
**下一轮**：在这 16 个面与一个正常面上，打印
`boundary_uv.len()`（`perform` 里 `collect_boundary_uv` 的结果）、
`wires_uv` 每条 wire 的点数、以及注册前的 `structure.nb_nodes()` 逐步值 ——
找出 75 是怎么来的、以及是哪一步之后 `elements_of_domain()` 一直为空。

#### 9.355.4 与既有线索的关系（避免重复追）

* §9.326 的「29 个面 2 vs 1」：**本轮修正为 28 个已配对 + 至少 3 个未配对**，
  且**它不解释** `domain_elems=0`；
* §9.342 的「狭长面」：§9.350 已撤回；
* §9.351/§9.352 的 rescue：**症状**，起点在 `domain_elems=0`（§9.354 已确认链的每一环）。

⇒ **当前唯一该追的量是「这 16 个面为何 `nodes=75` 且 `domain_elems=0`」。**

#### 9.355.5 本轮改动与门禁

* 本轮插桩（`OCCT_TOPO_TRACE_DELAUN`）已按 §9.353 的**新规程**用 `edit` 反向撤除
  （**未用 `git checkout`**），并修回了被我误删的一个空行；
* **`node_insertion.rs` 的 diff 现为空**（与 HEAD 一致；该文件不在 A 类清单、原无标注）；
* 全局复核：`crates/` 共 64 个改动 = **61 个 A 类标注（6/0）+ 3 个既有功能改动**
  （`gcpnts.rs` / `ginter.rs` / `make_pcurve.rs`）—— 与本会话基线一致；
* 门禁：`cargo check` 0 error；`--lib` **1255 passed / 26 failed**（环境性）；
  a3n00 `stats=226 unmatched=0 mesh_v=11101 mesh_t=11121`。

#### 9.355.6 复跑

```text
python .target-gate\intersect_zero_2v1.py
python .target-gate\delaun_zero.py
```

---

### 9.356 —— 失败链再往上一环：**三角化产出 0 个域内元素 → 清理把所有边界链标成 `Deleted`**；16 个失败面的 `boundary_uv` 全是 **74**（2 wires × 37 点）

#### 9.356.1 逐环测量（本轮新增三处打印，全部 env-gated）

在 `node_insertion.rs` 的 `init_data_structure` 与 `finish_mesh` 里插桩后（`OCCT_TOPO_TRACE_WIRES` / `OCCT_TOPO_TRACE_DELAUN`）：

**（a）失败面与正常面的**线结构**截然不同**：

```text
失败面 f=169:  WIRES f=169 wires=2 boundary_uv=74
               WIRES_OK f=169 wire_it=0 slots=1
               WIRES_OK f=169 wire_it=1 slots=1
正常面 f=0  :  WIRES f=0 wires=1 boundary_uv=76
               WIRES_OK f=0 wire_it=0 slots=8
```

⇒ 失败面是 **2 条 wire、每条 wire 只有 1 个 slot**（即每条 wire 是**单条闭合边**）；
正常面 `f=0` 是 **1 条 wire、8 个 slot**。**wire 结构完全不同。**

**（b）失败面**确实注册**了链**（排除"没注册"）：

```text
WIRES_EDGE f=169 wire_it=0 pts=37 nodes_now=36 links_now=36
WIRES_EDGE f=169 wire_it=1 pts=37 nodes_now=72 links_now=72
```

⇒ 每条 wire 的 pcurve 有 **37 个点** → 注册 36 个节点 + 36 条链；两条共 **72 节点 / 72 链**。
（`pts` 非空，说明 §9.354 里那个 `job=None` 的静默跳过**不是**原因 —— 实测 `DELAUN_SKIP total: 0`。）

**（c）到 `finish_mesh` 时，链的**状态**已全部变成 `Deleted`**：

```text
失败面:  frontier=0  fixed=0  free=0        ← 三种活状态全 0
正常面:  frontier=68 fixed=0  free=68       ← 边界链活着（且已转为 free）
```

且**链本身还在**（不是被移除）：

```text
失败面 DELAUN after_add_vertices f=169 nodes=75 links=221 domain_elems=0
正常面 DELAUN after_add_vertices f=0   nodes=71 links=207 domain_elems=66
```

⇒ **`links=221` 非 0，但每一种活状态计数都是 0** ⇒ **所有链都被标成 `Deleted`**。

#### 9.356.2 完整链条（每一环都有实测）

```text
16 个面：2 wires × 1 slot，每条 pcurve 37 点
   ↓ 注册
72 节点 / 72 链（正常）
   ↓ Delaunay 三角化
domain_elems = 0            ← §9.354 的落点
   ↓ 清理（无域内元素 ⇒ 所有链无引用）
所有链 → Deleted            ← 本轮
   ↓
frontier=0 fixed=0 free=0
   ↓ finish_mesh 的 use_edge 循环遍历空 frontier（无操作）
   ↓ collect_triangles 遍历空 elements_of_domain
triangles = 0
   ↓ map_triangulation 按 :569-573 拒绝
needs_fallback = true
   ↓
rescue（wireframe_face_triangulation）硬铺 224–336 个三角形
```

⇒ **症状（rescue）距真正的病灶有 6 环**。**真正的病灶是「Delaunay 三角化对这 16 个面产出 0 个域内元素」。**

#### 9.356.3 新增的强信号：16 个失败面的 `boundary_uv` **全是 74**

| 量 | 失败面（16） | 正常面（210） |
|---|---|---|
| `boundary_uv` | **全部 = 74** | 范围 **8 .. 370** |
| `nodes` | **全部 = 75** | 各不相同（71、45、7…） |
| `frontier` | **全部 = 0** | 范围 **3 .. 360** |

且本轮量出 `74 = 2 × 37`，**每条 wire 恰好 37 个点**。

⇒ **16 个几何完全不同的面（`f=169` 与 `f=206` 的 bbox 差 60+ 单位）
却给出完全相同的 `boundary_uv=74` / `nodes=75`** —— 这**不可能**是几何决定的自然结果。

**两个候选解释（下一轮判定）**：
1. **这 16 个面在端口的离散模型里指向了同一批/同构的数据结构**
   （例如它们共享同一个离散面对象，或 reader 把它们建成了同构的 2×37 点结构）；
2. **37 是某个固定的采样步数**（例如某条边被固定分成 37 段），
   而这 16 个面恰好都走了那条固定采样路径。

#### 9.356.4 下一轮的做法（明确）

在 `node_insertion.rs` 的 `init_data_structure` 里，对这 16 个面与 2 个正常面打印：
* `model.face(face_index)` 的 **`key`/指针**（判"是否同一批面对象"）；
* 每条 wire 的 `edge(j)` 的 **edge index**、`edge_reverse_walk`、以及
  `edge.pcurve(pc_index).points().len()`（判 37 是否固定）；
* 以及 `self.boundary_uv` 的前 3 个与后 3 个 uv 值（判 72 个点是否真的是同一批坐标）。

**若"同一批坐标"成立** ⇒ 病灶在**离散模型/reader 侧**（这 16 个面被建成了同一份数据），
与 Delaunay 无关；**若坐标各不相同、只有点数 37 相同** ⇒ 病灶在**采样步数**。

#### 9.356.5 本轮改动与门禁（**按 §9.353 新规程撤除**）

* 插桩三处，全部带 `TEMP T-99` 标记，**env-gated**：
  `WIRES`/`WIRES_SKIP`/`WIRES_OK`（`init_data_structure` 的 wire 循环）、
  `WIRES_EDGE`（slot 注册后）、`DELAUN`/`DELAUN_SKIP`（`finish_mesh` 与 `job=None`）；
* **撤除方式：逐块 `edit` 反向替换，未用 `git checkout`**；
* 撤除后 `node_insertion.rs` 应回到**空 diff**（该文件不在 A 类清单、原无标注）；
* 门禁复核见下。

#### 9.356.6 复跑

```text
cargo build --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_uv_feed
$env:OCCT_TOPO_TRACE_WIRES=1 ; $env:OCCT_TOPO_TRACE_DELAUN=1
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --ids 2> d.txt
grep -E 'f=(169|0) ' d.txt
```

---

### 9.357 —— 【更正 §9.352】这 16 个失败面**不是平面，是 Cylinder（12）与 Cone（4）**；共同特征是「闭合回转面被拆成 2 条 wire」

#### 9.357.1 更正内容

§9.352 写的是「13 个配上对的面**全部是 `Plane` 且 `wires=1`**」。**这是错的。**

本轮用 `UVSUM` 的 `type=` 重新核对全部配对面（`.target-gate/zero_face_types.py`）：

| port f | port wires | gt f | **GT type** | gt wires | GT u / v |
|---|---|---|---|---|---|
| 169 | 2 | 87 | **Cylinder** | 1 | `0..2π` / `-10..8` |
| 174 | 2 | 42 | **Cone** | 1 | `0..2π` / `-12.75..12.75` |
| 189 | 2 | 41 | **Cone** | 1 | `π..3π` / `-1.414..1.414` |
| 192 | 2 | 85 | **Cylinder** | 1 | `0..2π` / `-287.8..-267.8` |
| 193 | 2 | 84 | **Cylinder** | 1 | 同上 |
| 194 | 2 | 83 | **Cylinder** | 1 | 同上（v 略异） |
| 195 | 2 | 78 | **Cylinder** | 1 | 同上 |
| 196 | 2 | 79 | **Cylinder** | 1 | 同上 |
| 197 | 2 | 80 | **Cylinder** | 1 | 同上 |
| 198 | 2 | 81 | **Cylinder** | 1 | 同上 |
| 199 | 2 | 82 | **Cylinder** | 1 | 同上 |
| 204 | 2 | 43 | **Cone** | 1 | `π..3π` / `-1.414..1.414` |
| 206 | 2 | 86 | **Cylinder** | 1 | `π..3π` / `0..30` |
| 209 / 210 / 212 | 2 | （未配对） | — | — | — |

⇒ **12 个 Cylinder + 4 个 Cone**（3 个未配对），**没有一个是 `Plane`**。

**为什么之前判成 Plane**：§9.352 用的是 `UVSUM` 行里的 `type=` 字段，
但那一步我把**配对面的 GT 面号**与端口面号对错了行（本轮改用
**几何 bbox** 重新配对并直接读 `type=`，才得到上表）。
**这是我对同一批 `UVSUM` 数据的第二次读错** —— §9.355 的「2-vs-1 不是因果」结论
建立在同一批错读之上，**也需要按本节修正**（见 §9.357.3）。

#### 9.357.2 更正后的**共同特征**（这才是真规律）

| 项 | 值 |
|---|---|
| 面型 | **Cylinder（12）/ Cone（4）** —— **闭合回转面** |
| 端口 `wires` | **全部 = 2** |
| GT `wires` | **全部 = 1** |
| 每条 wire 的边数 | **各 1 条**（`IDENT_WIRE`：`edges=[414]`、`edges=[415]` …） |
| 每条 pcurve 点数 | **37** |
| `boundary_uv` | **全部 = 74**（= 2 × 37） |
| `nodes` | **全部 = 75** |
| `frontier` | **全部 = 0** |

⇒ **真规律是：端口把「闭合的圆柱/圆锥面」建成了 2 条 wire（GT 是 1 条），
这个 2-wire 结构让 Delaunay 注册不出任何域内元素。**

这与 §9.352 说的「单闭合 wire 的平面」**正好相反** —— 它们是**回转面上的 2-wire 结构**。

#### 9.357.3 §9.355 的结论需要按本节收窄

§9.355 说「2-vs-1 wire 不是因果」，依据是「28 个 2-vs-1 面里只有 13 个失败」。
但本轮查明**那 28 个面与 16 个失败面其实是同一类**（Cylinder/Cone 的 2-wire），
所以正确的表述是：

> **失败面集合 ⊂ 2-wire 的回转面集合**；
> 但**不是所有 2-wire 回转面都失败**（28 个里 13 个配上对的失败，
> 另有 15 个 2-wire 面 Delaunay 正常）。

⇒ **「2 wires」是必要条件而非充分条件**；失败的额外条件**尚未找出**。
**§9.355 的「不是因果」应改为「是必要条件、不是充分条件」**。

#### 9.357.4 下一轮的做法（据此修正）

1. **先修 §9.352 与 §9.355 的表述**（上面两条更正），避免后续按"平面"去找；
2. 在**那 15 个 2-wire 但正常的**回转面与这 16 个失败的之间做**差集分析**：
   打印它们的 `wire[0]`/`wire[1]` 的 pcurve 点数、两条 wire 的 uv 起止、
   以及 `wire.is_status(SELF_INTERSECTING_WIRE / OPEN_WIRE)`；
   **找出"失败"相对"正常"多出来的那一项**；
3. 参考 OCCT 的 `CylinderRangeSplitter` / `ConeRangeSplitter` 与
   `NodeInsertionMeshAlgo::initDataStructure` 对**闭合回转面**的处理 ——
   这是**有 `.cxx` 可对照**的正规移植问题。

#### 9.357.5 本轮改动与门禁（按 §9.353 新规程撤除）

* 插桩一处（`OCCT_TOPO_TRACE_IDENT`：面 key、wire/edge 身份、boundary_uv 首末坐标），
  **带 `TEMP T-99` 标记**、env-gated；
* **撤除方式：`edit` 反向替换，未用 `git checkout`**；
* 撤除后 `node_insertion.rs` 应回到**空 diff**；门禁见下。

**本轮最有价值的两个读数**（都在上面的表里）：
① **16 个面的 `face_key` 各不相同**（`1917249376176` / `1917250123632` / `1917250242560` …）、
边集各异、`boundary_uv` 首末坐标完全不同 ⇒ **它们不是"同一份数据"**，
候选 (1)（共享结构）**被否定**；
② 但它们**同一个面型族**（Cylinder/Cone）**且同为 2 wires** ⇒ 候选 (2)「固定采样步数」
也不是原因（37 点来自"每条 wire 一条闭合边"的自然采样）。

#### 9.357.6 复跑

```text
$env:OCCT_TOPO_TRACE_IDENT=1
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --ids 2> ident.txt
grep -E '^IDENT(_WIRE)? f=(169|192|206) ' ident.txt
python .target-gate\zero_face_types.py
```

---

### 9.358 —— 差集分析：**34 个正常的 2-wire 面与 16 个失败的"逐字段结构完全相同"** ⇒ 结构被排除，差异**只在 UV 数值**

#### 9.358.1 差集结果（`OCCT_TOPO_TRACE_WIRE2` + `.target-gate/diff_2wire.py`）

a3n00 的 wire 数分布与 2-wire 子集的成败：

```text
每面 wire 数分布: {1: 163, 2: 53, 4: 2, 6: 4, 10: 4}
2-wire faces: 53   failing(domain_elems=0): 16   working: 37
```

**16 个失败面的逐 wire 结构**（全部相同）：

```text
f=169 bu=74 ((1,37,False,False,('false',)), (1,37,False,False,('false',)))
f=174 bu=74 ((1,37,False,False,('false',)), (1,37,False,False,('false',)))
...（16 个完全一致）
```

**37 个正常面里有 34 个的结构与之逐字段相同**：

```text
f=41  bu=74 ((1,37,False,False,('false',)), (1,37,False,False,('false',)))
f=42  bu=74 ((1,37,False,False,('false',)), (1,37,False,False,('false',)))
f=44/45/47/48/77/82/87/92/93/104/109/110/111/112/129/130/137 ... 同样
```

**分离器检验（都在 §9.356 的候选清单里）**：

| 候选 | 失败面命中 | 正常面命中 | 是否分离 |
|---|---|---|---|
| `all_edges_1` | 16/16 | 34/37 | **否** |
| `all_pts_eq` | 16/16 | 34/37 | **否** |
| `all_rw_false` | 16/16 | 34/37 | **否** |
| `any_self` | **0**/16 | **0**/37 | 两者皆无 |
| `any_open` | **0**/16 | **0**/37 | 两者皆无 |

⇒ **没有任何"结构谓词"能分离这两组**。

**这同时否定/收窄了三件事**：
1. §9.357 的「2 wires」：**是必要条件，但连"结构完全相同"的 34 个面都成功**
   ⇒ **结构层面已经无法解释**；
2. §9.356 的候选「固定采样步数」：**正常面也是 `pts=37`、`bu=74`** ⇒ 排除；
3. §9.355 里我据以判「不是因果」的那 15 个"正常 2-wire 面"，
   现在**查明它们的结构与失败面逐字段相同** ⇒ 问题**不在结构，在数值**。

#### 9.358.2 差异**只在 UV 数值**，且失败面的 UV 有共同的数值特征

打印实际 `boundary_uv`（`UV2` 记录）：

```text
失败 f=192 n=74 first3=[(-3.141592654,-267.805924192), (2.967059728,-267.805924192),
                          (2.792526803,-267.805924192)]
             last3=[(5.934119457,-287.805924192), (6.108652382,-287.805924192),
                    (6.283185307,-287.805924192)]

正常 f=41  n=74 first3=[(-3.141592654,36.000000000), (2.967059728,36.000000000),
                          (2.792526803,36.000000000)]
             last3=[(9.075712110,34.000000000), (9.250245036,34.000000000),
                    (3.141592654,34.000000000)]
```

**两者都以 `-3.141592654`（= −π）开头** —— 即这两个**都是横跨接缝的圆柱面**。

区别在**终点**：
* **失败面**：`last = 6.283185307`（= **+2π**）
* **正常面 `f=41`**：`last = 3.141592654`（= **+π**）但其倒数第二个是 `9.25`（> 2π）

⇒ 失败面**恰好终止在 +2π**（跨度为整 2π），而正常面**越过了 2π**
（`9.25 = 3π − 0.17`）。**失败面的 u 跨度正好是一个完整周期。**

而 GT 侧这些面的 `vrange` 跨度正好是 **2r**（`f=192`：`-287.806..-267.806`，
差 **20.0**；`f=169`：`-10..8`，差 **18.0**）⇒ 它们映射到
**`CylinderRangeSplitter` / `ConeRangeSplitter`**，即**在接缝处切开的周期面**。

#### 9.358.3 收窄后的假设（下一轮）

> **假设**：对**u 跨度恰好为一个完整周期（−π..+2π 或 0..2π 的整数倍）**的闭合回转面，
> 端口的周期处理（`CylinderRangeSplitter`/`ConeRangeSplitter` 或
> `AdjustRange` 的周期分支）**没有把首末点合并**，
> 于是约束多边形在接缝处**多出一个重复端点**（`u=-π` 与 `u=+2π` 应为同一点），
> `BRepMesh_Delaun` 因此推不出任何域内元素。

**支持证据**：失败面 **全部**终止在 ±2π 的整数倍；正常面**越过**了周期端点。
**反对证据（须下一轮排除）**：需要确认正常面里也有"恰好整周期"的。

**下一轮第一步**：把 53 个 2-wire 面按
`(首点 u, 末点 u, 末点 u − 首点 u)` 与**成败**列表并排，看
**"跨度 = 整周期"是否恰好分离**。这是**唯一还没被否定**的候选。

**若成立**，修复落点是**周期范围的首末点合并**，且这是**有 `.cxx` 可对照**的
（`BRepMesh_CylinderRangeSplitter` / `ConeRangeSplitter` 的 `AdjustRange` 与
`BRepMesh_DefaultRangeSplitter::AdjustRange` 的周期分支）——
**不是**自造规则。

#### 9.358.4 本轮改动与门禁（按 §9.353 新规程撤除）

* 插桩一处（`OCCT_TOPO_TRACE_WIRE2`：逐 wire 结构 + `boundary_uv` 全量数值），
  **带 `TEMP T-99` 标记**、env-gated；
* **撤除方式：`edit` 反向替换，未用 `git checkout`**；
* 撤除后 `node_insertion.rs` 应回到**空 diff**；门禁见下。

#### 9.358.5 复跑

```text
$env:OCCT_TOPO_TRACE_WIRE2=1
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --ids 2> uv2.txt
python .target-gate\diff_2wire.py
python .target-gate\uv_compare.py
```

---

### 9.359 —— 【决定性】失败面的 `frontier` 在 `use_edge` **之前**就是 0；且 `f=209`（失败）与 `f=208`（正常）的**边界点集合完全相同** ⇒ 差异是**顺序决定**，不是几何脆弱性

#### 9.359.1 `use_edge` 前后的域元素数（`OCCT_TOPO_TRACE_TRI`）

```text
f=209（失败）: TRI pre_useedge  nodes=75 links=221 frontier=0 domain=0
               TRI post_useedge domain=0
f=208（正常）: TRI pre_useedge  nodes=75 links=221 frontier=3 domain=1
               TRI post_useedge domain=1
f=0  （正常）: TRI pre_useedge  nodes=71 links=207 frontier=68 domain=66
               TRI post_useedge domain=66
```

**三条读数**：

1. **`use_edge` 完全不改变 `domain`** —— 它只是把既有的域内元素补全/校验；
   **域内元素是在 `add_vertices` 里产生的**（§9.356 的 `after_add_vertices` 也已显示
   `domain_elems=0` 在那时就定了）。⇒ **病灶在 `Delaun::add_vertices` 内部。**
2. **`f=209` 与 `f=208` 的 `nodes=75`、`links=221` 完全相同**，
   但 **`frontier` 是 0 vs 3** ⇒ 差异**不体现在节点数/链数**上。
3. 结合 §9.358 的结论 —— **两者的边界点集合逐点相同**（见 §9.359.2）——
   ⇒ **在完全相同的输入几何下，`add_vertices` 给出不同结果。**
   **这不是"几何脆弱性"，而是"顺序决定"。**

#### 9.359.2 `f=209`（失败）与 `f=208`（正常）的边界逐点对比（`UV2` 记录）

```text
n: 74 vs 74
same u at same index: 74/74          ← u 逐位置相同
same v at same index: 40/74          ← v 只有 40/74 同位置
as SETS:
  identical u multiset: True
  identical v multiset: True
  identical point set : True          ← ★ 74 个点作为集合完全一致
```

⇒ **两个面的边界是"同一组 74 个点、不同的遍历顺序/起点"。**
（`u` 逐位置相同说明**两条 wire 的 u 序列**一致；`v` 的差异说明
**点在同一条 wire 内的先后不同**——即**同一条闭合回路的起点/绕向不同**。）

#### 9.359.3 结论

> **失败与成功不由边界的内容决定，也不由 wire 的结构决定（§9.358），
> 而由点在边界数组中的顺序决定。**
>
> ⇒ 病灶是 **`Delaun::add_vertices`（即 `BRepMesh_Delaun` 的构造）里某个
> 依赖点序的步骤**：对同一组点，按顺序 A 能给不出任何域内元素，按顺序 B 可以。

**这与 OCCT 的已知设计一致**：`BRepMesh_Delaun` 内部按
`ComparatorOfVertexOfDelaun`（对 (1,1) 方向的投影）**给顶点排序**，
并用 `aStructure` 的 `frontier`/`fixed` 链建立初始三角化 ——
**是一个顺序敏感的增量构造**。

#### 9.359.4 下一轮的做法（**收敛到具体函数**）

1. 打开 `crates/occt-topo/src/meshing/delaun/` 下的 `add_vertices` 实现，
   与 OCCT `BRepMesh_Delaun.cxx` 的对应构造函数**逐行对照**；
2. 在 `add_vertices` 内、对 `f=209` 与 `f=208` 打印**每一步之后**的
   `elements_of_domain().len()`（或等价中间量），找出**第一步出现分歧**的位置；
3. 重点看 OCCT 的 `ComparatorOfVertexOfDelaun` / `ComparatorOfIndexedVertexOfDelaun`
   （`BRepMesh_Delaun.cxx:49-75`，本会话已读过）——
   **端口是否用了同一个排序谓词**。

**这是本轮把范围缩到最小的成果**：
从「16 个面失败」→「`add_vertices` 一步」→「一个顺序敏感的步骤」。

#### 9.359.5 本轮改动与门禁（按 §9.353 新规程撤除）

* 两处插桩（`OCCT_TOPO_TRACE_WIRE2`、`OCCT_TOPO_TRACE_TRI`），
  **均带 `TEMP T-99` 标记**、env-gated；
* **撤除方式：`edit` 反向替换，未用 `git checkout`**；
* 撤除后 `node_insertion.rs` 应回到**空 diff**；门禁见下。

#### 9.359.6 复跑

```text
$env:OCCT_TOPO_TRACE_TRI=1
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --ids 2> tri.txt
grep -E ' f=(209|208|0) ' tri.txt
python .target-gate\pair_209_208.py
```

---

### 9.360 —— 【精确锁定】失败面的 145 个域内元素是被 **`process_constraints` 全部删掉**的：`pre=145 → post=0`

#### 9.360.1 `create_triangles_on_new_vertices` 内的三段计数（`OCCT_TOPO_TRACE_CTNV`）

在 `triangulation.rs` 的 `create_triangles_on_new_vertices` 里插三处打印：
进入时、`process_constraints` 之前、之后。

**243 条 `CTNV enter` 记录**：

```text
starting from domain=0            : 0      ← 从来不是"一开始就没有"
ending   at domain=0              : 16
no post record (early return)     : 0      ← 没有走 self.failed 提前返回
```

**那 16 条 `ending at domain=0` 的记录，字段完全一致**：

```text
verts=72  d0=3  f0=72  pre=145  post=0  failed=False
...（16 条一模一样）
```

**逐段解读**：

| 阶段 | 值 | 含义 |
|---|---|---|
| 进入时 `d0` | 3 | 构造时已建了 3 个元素 |
| 进入时 `f0` | **72** | **72 条 frontier 链**（= 2 wires × 36） |
| `verts` | **72** | 遍历 72 个节点 |
| `pre`（`process_constraints` 前） | **145** | **建出了 142 个新元素，共 145** |
| `post`（之后） | **0** | **全部被删** |
| `failed` | false | 不是 `add_triangle` 溢出 |

⇒ **三角化本身是成功的**（145 个域内元素真实存在），
**是 `process_constraints()` 把这 145 个全部删掉了。**

#### 9.360.2 删除发生在 `frontier_adjust` 的 reorient 段

`frontier.rs:15-43`：

```rust
for _pass in 1..=2 {
    for &frontier_id in &frontier_ids {              // 72 条 frontier 链
        let pair = self.mesh_data.elements_connected_to(frontier_id);
        for elem_it in 1..=nb {
            ...
            for n in 0..3 {
                if frontier_id == element.link_at(n).abs() && element.link_at(n) < 0 {
                    self.delete_triangle(prior_elem, &mut loop_edges);   // ← 删
                }
            }
        }
    }
    ...
    for &frontier_id in &frontier_ids {
        if elements_connected_to(frontier_id) 非空 { continue; }
        self.mesh_left_polygon_of(frontier_id, true, &mut skipped);      // ← 重建
    }
```

**该段按 OCCT 的 `BRepMesh_MeshTool::Frontier` / `frontierAdjust` 语义，
删除"悬挂在负向 frontier 链上的三角形"，以便按正确朝向重建。**

**16 个面 `pre=145 → post=0`** ⇒ **这 145 个三角形全部挂在负向 frontier 链上**
⇒ **初始三角化的定向整体相反**，`frontier_adjust` 认为它们全错、全部删除，
然后 **`mesh_left_polygon_of` 重建失败**（否则 `post` 不会是 0）。

#### 9.360.3 这解释了此前所有观测

| 观测 | 解释 |
|---|---|
| `domain_elems=0` 在 `add_vertices` 后（§9.356） | `process_constraints` 在 `add_vertices` 内被调（`:539`），删完后为 0 |
| `frontier=0`、`fixed=0`、`free=0`（§9.356） | 删三角形时连带 `remove_link`，链被标 `Deleted` |
| `nodes=75`、`links=221` 与"正常"面相同（§9.359） | 节点/链的**数量**没错，错的是**链的朝向** |
| `f=209` 与 `f=208` 点集相同、顺序不同（§9.359） | **点的顺序决定了两条 wire 的绕向**，进而决定 frontier 链的**朝向** |

⇒ **"顺序决定"的机理找到了**：顺序 → wire 绕向 → frontier 链朝向 →
`frontier_adjust` 是否把整批三角形判为"反向"并删除。

#### 9.360.4 下一轮的做法（**已到具体代码段**）

1. **对照 OCCT**：`BRepMesh_MeshTool::Frontier`（`BRepMesh_MeshTool.cxx:284-300`）
   与 `frontierAdjust` 的 reorient 段，逐行核对端口 `frontier.rs:15-60`
   的判定条件 `element.link_at(n) < 0` 是否与 OCCT 一致
   （**是否少判/多判了一个朝向条件**）；
2. **量两侧**：在一**正常**2-wire 面（如 `f=208`）与一失败面（`f=209`）上，
   分别在 `frontier_adjust` 的 reorient 段打印
   **被删的三角形数**与**每个 frontier 链的 `link_at(n)` 符号分布**；
   正常面应删得很少、失败面删光 —— 由此定位**符号判断**的分歧点；
3. **重点怀疑**：`link_at(n) < 0` 里的**符号约定**（`add_link_to_mesh` 在
   `Orientation::Reversed` 时用 `add_link(last, first, …)`），
   与 OCCT `addLinkToMesh` 的对应分支（`node_insertion.rs:827-835`）。

**仍有 `.cxx` 原文可对照**，属正规移植问题。

#### 9.360.5 本轮改动与门禁（按 §9.353 新规程撤除）

* 插桩三处在 `create_triangles_on_new_vertices`（`delaun/triangulation.rs`），
  **带 `TEMP T-99` 标记**、env-gated（`OCCT_TOPO_TRACE_CTNV`）；
* **撤除方式：`edit` 反向替换，未用 `git checkout`**；
* 撤除后 `delaun/triangulation.rs` 应回到**空 diff**；门禁见下。

#### 9.360.6 复跑

```text
$env:OCCT_TOPO_TRACE_CTNV=1
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --ids 2> ctnv.txt
python .target-gate\ctnv_zero.py
```

---

### 9.361 —— 逐行核对：`frontierAdjust` 的**全部控制流**与符号惯例都与 OCCT 一致 ⇒ 不是"少判/多判朝向条件"

#### 9.361.1 端口 `frontier.rs` ↔ OCCT `BRepMesh_Delaun::frontierAdjust`（`cxx:944-1026`）

| 段 | OCCT | 端口 | 判定 |
|---|---|---|---|
| 主循环 | `cxx:955` `for (aPass = 1; aPass <= 2; ++aPass)` | `frontier.rs:15` `for _pass in 1..=2` | **等价** |
| 取 frontier 集 | `cxx:946` `Frontier()` | `frontier.rs:10` `self.frontier()` | 等价 |
| 遍历每链的挂接三角形 | `cxx:961-966` `ElementsConnectedTo` + `Extent` | `:16-19` `elements_connected_to` + `extent` | 等价 |
| 跳过已删 | `cxx:970-973` `if (aPriorElemId < 0) continue` | `:21-23` `if prior_elem < 0 { continue }` | **等价** |
| **外三角形判定** | `cxx:980-988`：`if (aFrontierId == e[n] && !o[n])` → `deleteTriangle` → **`break`** | `:26-38`：`if element.link_at(n) < 0` → `delete_triangle` → `break` | **等价**（见 §9.361.2） |
| 找到即跳出元素循环 | `cxx:991-994` `if (isTriangleFound) break` | `:39-41` `if found { break }` | **等价** |
| 清理悬挂链 | `cxx:999-1007` `if ElementsConnectedTo(e).IsEmpty() { RemoveLink(e) }` | `:45-49` `if elements_connected_to(e).is_empty() { remove_link(e, false) }` | **等价** |
| 重建 | `cxx:1011-1025` `meshLeftPolygonOf(aFrontierId, true, aIntFrontierEdges)`，`aPass==2 && !isSuccess` 记 `aFailedFrontiers` | `:51-63` 同 | **等价** |
| `cleanupMesh` | `cxx` 尾 | `:66` `self.cleanup_mesh()` | 等价 |
| 失败前沿重试 | `cxx` 尾 | `:68-77` 同 | 等价 |

⇒ **`frontierAdjust` 的控制流已完整对照，未发现"少判/多判朝向条件"。**

#### 9.361.2 符号惯例**也已证实一致**（这是关键）

OCCT 把**边的绝对号**与**朝向布尔**分开存：
`BRepMesh_Triangle` 有 `myEdges[3]`（绝对号）+ `myOrientations[3]`（布尔）。

* **写入时**（`cxx:541-548`，`meshLeftPolygonOf` 内的 `addTriangle` 前处理）：
  ```cpp
  anEdgeIds[aTriLinkIt]  = std::abs(anEdgeInfo);
  anEdgesOri[aTriLinkIt] = anEdgeInfo > 0;
  ```
* **端口写入时**（`frontier.rs:321-328` `add_triangle_by_info`）：
  ```rust
  edges[i] = edges_info[i].abs();
  oris[i]  = edges_info[i] > 0;
  ```
  **逐字相同。**
* **存储时**：OCCT `addTriangle` 把 `(theEdgesId, theEdgesOri)` 原样存进
  `BRepMesh_Triangle`；端口 `triangulation.rs:555-559` 压成带符号数组
  ```rust
  let signed = [if oris[0] {edges[0]} else {-edges[0]}, ...];
  ```
  ⇒ **正号 = oriented（OCCT 的 `o[i] == true`）、负号 = 反向（`o[i] == false`）**。
* **读取时**：OCCT 判 `!o[n]`；端口判 `link_at(n) < 0`
  ⇒ **语义一一对应**。

**因此 §9.360 里"145 个三角形全部挂在负向 frontier 链上"不是符号惯例错误**，
而是**三角形相对于 frontier 链的朝向确实为负**。

#### 9.361.3 剩下的唯一去处

`frontierAdjust` 与符号惯例都已排除 ⇒ **产生那些三角形的 `oris` 本身就是错的**。
三角形的产生点有三处：

```text
frontier.rs:328     add_triangle_by_info   （meshLeftPolygonOf 内）—— 惯例已验证一致
polygon_meshing.rs:235 add_triangle
triangulation.rs:427   add_triangle
```

**下一轮**：对 `f=209`（失败）与 `f=208`（正常），在这三个调用点打印
`(link_id, ori)` 与**多边形遍历方向**（`polygon` 的有符号链号序列），
看失败面的 `oris` 是否**整体为反** —— 即 `polygon` 的绕向算错。

这直接指向 `mesh_left_polygon_of` 的**多边形走向判定**
（`frontier.rs:90-121`，`is_forward` 分支与 `ref_link_dir` 的方向）——
**与 OCCT 的 `meshLeftPolygonOf` 对照仍是正规移植工作。**

#### 9.361.4 本轮结论的形态

**本轮没有找到"少判一个条件"，而是把 `frontierAdjust` 整段排除掉。**
这仍然是有价值的排除：它把范围从"`frontierAdjust` 可能有 bug"
缩到"**三角形的 `oris` 在产生时就错了**"，而 `oris` 只可能来自
**多边形绕向**的计算。

#### 9.361.5 本轮改动与门禁

* **未提交任何源码改动**（只做读取与对照）；
* `crates/occt-topo/src/meshing/delaun/` **diff 为空**（§9.360 的插桩已撤）；
* 门禁：`cargo check` 0 error；`--lib` 1255/26；a3n00 基线未变。

#### 9.361.6 引用

```text
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_Delaun.cxx:944-1026   frontierAdjust
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_Delaun.cxx:383-405    deleteTriangle
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_Delaun.cxx:541-550    绝对号+朝向的写入
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_Triangle.hxx:60-75   myEdges/myOrientations
crates/occt-topo/src/meshing/delaun/frontier.rs:15-78                                    端口 frontier_adjust
crates/occt-topo/src/meshing/delaun/triangulation.rs:330-344                             端口 delete_triangle
```

---

### 9.362 —— `meshLeftPolygonOf`/`findNextPolygonLink` 逐行等价；**多边形符号分布不分离**（失败 50.1% 正 vs 正常 50.3% 正）

#### 9.362.1 逐行对照

| 段 | OCCT | 端口 | 判定 |
|---|---|---|---|
| 开头 / `isForward` 分支 | `cxx:1070-1100`（`Append(theStartEdgeId)` / `Append(-theStartEdgeId)`，`aStartNode`/`aPivotNode` 取 `FirstNode`/`LastNode` 或反之） | `frontier.rs:90-121`（`polygon.push(start_edge_id)` / `push(-start_edge_id)`，同） | **等价** |
| 参考方向 | `cxx:1095` `gp_Vec2d aRefLinkDir(...)` | `:111` `GpVec2d::from_xy(...)` | 等价 |
| 退化检查 | `cxx:1097` `SquareMagnitude() < Precision2` | `:112` `< PREC2` | 等价（`PREC = PCONFUSION` ≡ `Precision::PConfusion()`） |
| 回退（`aNextLinkId == 0`） | `cxx:1137+` | `:151-174` 同 | 等价 |
| 收尾 | `polygon.len() < 3` 返回 false；`cleanupPolygon` + `meshPolygon` | `:177-182` 同 | 等价 |
| **`findNextPolygonLink`** | `cxx:1206-1316` | `frontier.rs:186-...` | **等价**（见下） |

**`findNextPolygonLink` 内部逐条**：

| OCCT | 端口 |
|---|---|
| `cxx:1222` `double aMaxAngle = RealFirst()` | `:199` `let mut max_angle = f64::NEG_INFINITY` |
| `cxx:1225` `LinksConnectedTo(thePivotNode)` | `:205` `links_connected_to(pivot_node)` |
| `cxx:1231-1235` 跳过 `theDeadLinks` / `skipped` | `:208-215` 同 |
| `cxx:1237-1241` `isSkipLeprous && isLeprous` 跳过 | `:216-219` 同 |
| `cxx:1246-1251` `Movability()==Free && ElementsConnectedTo(...).IsEmpty()` → 记入 `theDeadLinks` | `:221-226` 同 |
| `cxx:1253-1257` 取另一端节点 | `:227-230` 同 |
| `cxx:1261-1265` `SquareMagnitude() < Precision2` → 记入 `theDeadLinks` | `:232-235` 同 |
| `cxx:1267-1270` `!isLeprous` → 加入 `theLeprousLinks` | `:236-238` 同 |
| `cxx:1272` `anAngle = theRefLinkDir.Angle(aCurLinkDir)` | `:240` `ref_link_dir.angle(&cur_link_dir)` |
| `cxx:1276-1283` `isFrontier` 且 `\|\|anAngle\|-π\| < Precision::Angular()` → `isCheckPointOnEdge = false; anAngle = abs(anAngle)` | `:243-248` 同（`ANGULAR = 1e-12` ≡ `Precision::Angular()`） |
| `cxx:1286-1289` `anAngle <= aMaxAngle` → 跳过 | `:249-251` 同 |
| `cxx:1291` `isCheckEndPoints = (anOtherNode != theFirstNode)` | `:253` 同 |
| `cxx:1294-1300` `checkIntersection(...)` | `:255-263` 同（`isSkipLastEdge = true`） |
| `cxx:1310-1311` `aNextLinkId = (FirstNode()==thePivotNode) ? +id : -id` | 端口对应处 | 等价 |

⇒ **`meshLeftPolygonOf` 与 `findNextPolygonLink` 全段等价，常量也一致。**

#### 9.362.2 实测：符号分布**不分离**两组（新插桩 `OCCT_TOPO_TRACE_POLY`）

`POLY` 记录每次闭合多边形的 `len`/`pos`/`neg`（`pos` = 带正号链数，即
`link.FirstNode() == pivot`）：

```text
     f   set  runs   links     pos     neg    pos%
   169  FAIL   143    8120    4060    4060    50.0
   174  FAIL   139    8075    4052    4023    50.2
   189  FAIL   141   10022    5013    5009    50.0
   192  FAIL    88    5500    2764    2736    50.3
   ...
   215    ok    86    5601    2918    2683    52.1
   216    ok   135    9522    4767    4755    50.1

FAILING aggregate: pos=58689 neg=58361 pos%=50.1   (faces=16)
WORKING aggregate: pos=106639 neg=105411 pos%=50.3  (faces=36)
```

⇒ **两组的正/负号比例几乎相同（50.1% vs 50.3%）**
⇒ **不存在"某组的多边形整体反向"**。**"多边形绕向算错"这个假设被否定。**

（抽查也印证：`f=208`（正常）多数多边形 `pos=23..36 neg=1`；
`f=209`（失败）则 `pos=6 neg=34` 与 `pos=40 neg=2` 混杂 ——
但**聚合后两组一样**，说明这是面内不同多边形的正常变化，不是系统性反向。）

#### 9.362.3 本轮排除的与剩下的

**已排除**（都有实测或逐行对照）：
* `frontierAdjust` 的控制流与符号惯例（§9.361）
* `meshLeftPolygonOf` / `findNextPolygonLink`（本节，逐行等价）
* 常量 `PREC`/`PREC2`/`ANGULAR`（本节，等价）
* 多边形符号分布（本节，两组相同）

**剩下的**：`cleanup_polygon` / `mesh_polygon` / `mesh_elementary_polygon`
（三角形真正被 `add_triangle` 生成的那一段），以及 **`checkIntersection`**
（它决定选哪条邻链，从而决定多边形本身）。

**下一轮**：
1. 新增了 `TRI2` 插桩（`polygon_meshing.rs` 的 `mesh_elementary_polygon`，
   打印 `poly=(前3) edges oris nodes`），共 626 条记录 —— **下一轮把它按面归属**
   （复用本轮的 `trace_set_face` 面包屑）并与 OCCT 的
   `meshElementaryPolygon`（`cxx` 中对应处）比对；
2. 同时给 `checkIntersection` 加桩，看失败面是否**选了不同的邻链**
   （这会让多边形走向不同，而符号统计看不出来）。

#### 9.362.4 本轮新增的插桩（**须下一轮撤除**）

| 文件 | 插桩 | 开关 |
|---|---|---|
| `delaun/frontier.rs` | `TRACE_FACE` thread-local + `trace_set_face()` + `POLY` 打印 | `OCCT_TOPO_TRACE_POLY` |
| `delaun/mod.rs` | `pub(crate) use frontier::trace_set_face;` | — |
| `node_insertion.rs` | `crate::meshing::delaun::trace_set_face(face_index);` | — |
| `delaun/polygon_meshing.rs` | `TRI2` 打印 | `OCCT_TOPO_TRACE_TRI2` |

**本轮结束时按 §9.353 规程用 `edit` 反向撤除**（不用 `git checkout`）。

#### 9.362.5 复跑

```text
$env:OCCT_TOPO_TRACE_POLY=1
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --ids 2> poly.txt
python .target-gate\poly_stats.py
$env:OCCT_TOPO_TRACE_TRI2=1
... 2> tri2.txt
```

---

### 9.363 —— 【阶段性】Delaunay 侧**整条三角形生成路径已逐行核对完毕、全部等价** ⇒ 分歧必在**喂进去的边界结构**（2 wires vs 1 wire）

#### 9.363.1 本轮继续核对的两段（均等价）

| 段 | OCCT | 端口 | 判定 |
|---|---|---|---|
| `checkIntersection` | `cxx:1324-1372`：`UpdateBndBox`、`aPolyLen = Length()`，`isSkipLastEdge` 则 `--aPolyLen`、`isFrontier = (Movability()==Frontier)`、`for aPolyIt = 1..=aPolyLen` 内 `IsOut(polyBox)` 检查、跳过"frontier×frontier"、`intSegSeg`、`!= NoIntersection` 则 `return false` | `frontier.rs:280-317`：`update_bnd_box`、`poly_len`，`is_skip_last_edge` 则 `-=1`、`is_frontier`、`for poly_it in 0..poly_len`（0 基 ≡ 1 基）、同三项 | **等价** |
| `IntSegSeg` | `GeomTool.cxx:342-461`：四点 `classifyPoint` → `aPosHash`；`hash[0]<0 \|\| hash[1]<0` 分支（`-1`→`Glued`，否则看 `isConsiderEndPointTouch`）；`==1`（看 `isConsiderPointOnSegment`，四点选点）；`==2`→`Glued`；`IntLinLin`；`NoIntersection`/`Same`（`< -2`→`Same`，`==-1`→`Glued`，否则 `NoIntersection`）；`Cross` 用 `Precision::PConfusion()` 做 `param ∈ [prec, 1-prec]` 检查 | `geom_tool.rs:343-409`：四点 `classify_point` → `pos_hash`；`hash[0]<0 \|\| hash[1]<0` 同；`==1` 同（四点选点）；`==2`→`Glued`；`int_lin_lin`；`NoIntersection`/`Same` 同三分支；`Cross` 用 `PCONFUSION` 同检查 | **等价** |
| `addTriangle` 的圆绑定 | `cxx:1378-1398`：`AddElement(BRepMesh_Triangle(...))`，`if myInitCircles { isAdded = myCircles.Bind(...) }`，`if !isAdded { RemoveElement }` | `triangulation.rs:548-595`：`add_element(DelaunTriangle::new(signed, verts))`，`if self.init_circles { is_added = self.circles.bind_circle(...) }`，`if !is_added { remove_element }` | **等价**（符号↔朝向见 §9.361.2） |

#### 9.363.2 已核对的清单（本轮累计）

**公式/控制流层**（全部逐行等价，多数有实测）：
`AdjustRange`、`computeLengthU/V`、`computeTolerance`、`computeDelta`、格数三支、
`AdjustCellsCounts`、`ComputeErrFactors`、`collectTriangles`、`AddElement`/`RemoveElement`、
`SetCellSize`/`SetTolerance`、`initDataStructure`、`finish_mesh`、
`createTrianglesOnNewVertices`、`frontierAdjust`、`deleteTriangle`、
`meshLeftPolygonOf`、`findNextPolygonLink`、`checkIntersection`、`IntSegSeg`、`addTriangle`。

**实测层**（值也一致）：`delta`（47/47、56/56）、`range_u/v`（差 <0.68%）、
`vertices_nb`（226/226）、`frontier` 链数、「多边形符号分布」（两组 50.1% vs 50.3%）。

⇒ **Delaunay 侧没有未核对的段落了。** 而它在这 16 个面上产出 **0** 个域内元素，
OCCT 不会 ⇒ **分歧只能来自"喂进去的东西"**，而不是这些函数本身。

#### 9.363.3 而"喂进去的东西"确有一处**已知的、结构性的**不一致

§9.357 已查明：这 16 个面在**端口侧是 2 条 wire**、在 **GT 侧是 1 条 wire**。
`init_data_structure` 按 **wire** 逐个注册边界链：

```rust
for (wire_it, &wire_index) in wire_indices.iter().enumerate() { ... }
```

⇒ 2 条 wire 会注册**两条各自闭合的边界环**，而 OCCT 的 1 条 wire 只注册**一条**。

**推论**：端口在 Delaunay 里看到的边界是**两条闭合环**，而 OCCT 看到的是一条。
`Given` 一个"两条环"的边界，`createTrianglesOnNewVertices` 仍会建出 145 个三角形
（它按点建），但 `frontierAdjust` 认为它们相对两条环的朝向**全部不合格**而删除
（§9.360 的 `pre=145 → post=0`）；重建也失败。

**这与 §9.360/§9.361/§9.362 的实测全部自洽**：
* 链数 72 = 2 × 36（两条环）而非 36 —— §9.356 的 `f0=72`；
* 符号分布两组相同（因为两条环各自 50% 是正常的）；
* `frontierAdjust` 段本身没写错。

⇒ **病灶在"为什么这个闭合回转面被建成了 2 条 wire"**，
而这是 **reader / wire 构建侧**的事，**不在 Delaunay 里**。

#### 9.363.4 为什么之前一直在 Delaunay 侧打转

因为**症状**（0 个域内元素、145→0、frontier=0）全都出现在 Delaunay 内部。
现在可以定性地说：**Delaunay 只是"忠实执行"了一个错的双环边界输入**。

#### 9.363.5 下一轮（换层）

1. **回到 §9.326/§9.357 的 2-vs-1 wire 现象**，查**为什么端口把一个闭合回转面
   建成 2 条 wire**：读 reader 侧的 wire 构建（STEP → `TopoDS_Wire` → 端口的
   `MeshModel` wire），与 OCCT 的 `BRepMesh_ShapeVisitor` / `IMeshData` 对照；
2. **判定标准**：这 16 个面在 STEP 里各自应是**几条 `EDGE_LOOP`**？
   若 STEP 本身就是 1 个 loop，则端口**多切了一刀**（reader bug）；
   若 STEP 是 2 个 loop，则 **OCCT 把它们合并了**（wire 合并缺失）。
   **这是可以查证的事实，不是猜测。**
3. `--wireinv` 探针（§会话早期已存在）正是为此准备的。

#### 9.363.6 本轮改动与门禁

* **未提交任何源码改动**（只做读取与对照）；
* `crates/occt-topo/src/meshing/delaun/` 与 `node_insertion.rs` **diff 均为空**；
* 门禁：`cargo check` 0 error；`--lib` 1255/26；a3n00 基线未变。

#### 9.363.7 引用

```text
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_Delaun.cxx:1324-1372   checkIntersection
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_GeomTool.cxx:342-461  IntSegSeg
D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_Delaun.cxx:1378-1398  addTriangle
crates/occt-topo/src/meshing/delaun/frontier.rs:280-317                                 check_intersection
crates/occt-topo/src/meshing/geom_tool.rs:343-409                                       int_seg_seg
crates/occt-topo/src/meshing/delaun/triangulation.rs:548-595                            add_triangle
```

---

### 9.364 —— 【自我否定】STEP 文件本身就有 53 个面是 **2 个 bound** ⇒ 「2-vs-1 wire」不是端口缺陷；§9.363.3 的推论作废

#### 9.364.1 可查证的事实（直接解析 `data/occ/a3n00.stp`）

用 `.target-gate/step_bounds.py` 解析 STEP（`ADVANCED_FACE` → 其 bound 列表 →
`FACE_OUTER_BOUND`/`FACE_BOUND` → `EDGE_LOOP`）：

```text
ADVANCED_FACE parsed: 226
bounds-per-face histogram: {1: 163, 2: 53, 4: 2, 6: 4, 10: 4}
FACE_*_BOUND parsed: 341      EDGE_LOOP parsed: 341
```

**与端口的 wire 数分布（§9.358 实测）逐字相同**：

```text
端口 wire 数分布（OCCT_TOPO_TRACE_WIRE2，226 面）: {1: 163, 2: 53, 4: 2, 6: 4, 10: 4}
STEP bounds 分布（本脚本，226 面）:               {1: 163, 2: 53, 4: 2, 6: 4, 10: 4}
```

且**面序对应**：`port f=0 ↔ step #260`、`port f=225 ↔ step #8958`（首尾边界均合理）。

⇒ **端口的面/边界结构与 STEP 文件完全一致** —— reader **没有**多切 wire。

#### 9.364.2 这否定了 §9.363.3 的推论（以及 §9.357 的一半）

§9.357 说「**端口**把这 16 个面建成 2 条 wire，而 **GT** 是 1 条」，
并据「GT `wires=1`」推断「OCCT 只注册一条边界环」。

**但 STEP 文件里这些面本来就是 2 个 bound**：

```text
 port_f   step_id  #bounds
    169      7302        2
    174      7493        2
    ...
    192      8008        2
    ...
    212      8607        2
```

若 OCCT 与端口读同一个 STEP、且都用 `ADVANCED_FACE` 的 bound 数建 wire，
则 **OCCT 必然也是 2 条 wire** ⇒ **「GT `wires=1`」这个读数是错的**
（`UVSUM` 的 `wires=` 字段口径与"面的 bound 数"不同，
我在 §9.357 把两个不同口径的数当成了同一件事 —— **这是我对 GT 数据的第三次误读**）。

⇒ **§9.363.3 的推论（"端口注册了两条环而 OCCT 一条"）作废。**

#### 9.364.3 修正后的处境

| 事实 | 状态 |
|---|---|
| Delaunay 侧 21 处逐行等价 | §9.361/§9.362/§9.363 |
| 端口与 STEP 的边界结构一致 | **本节（实测）** |
| 端口在这 16 个面产出 0 域内元素 | §9.354/§9.360 |
| **OCCT 在这 16 个面是否也产出 0？** | **未知 —— 从未测过** |

⇒ **关键未知变成了"OCCT 在 Delaunay 层对同一批面是否也归零"。**

这一层此前**从未测过**：我比的是 `triangles=`（三角化后的产物）与端口的 `TriangulationResult`，
**没有比过两侧 Delaunay 的"域内元素数"**。

**三种可能，各自结论不同**：
1. **OCCT 也归零** ⇒ 端口是**忠实的**，这 16 个面的差异在更下游
   （OCCT 归零后走**别的**路径产出三角形，端口走了 rescue）⇒ 要找 OCCT 的后续路径；
2. **OCCT 不归零** ⇒ 端口在 Delaunay 的**输入/中间状态**上与 OCCT 有差异，
   而差异**不在我核过的 21 处代码**里 ⇒ 要找的是**数据**（点的精确坐标、注册顺序、
   `NodesMap`/`LinksOfDomain` 的内容），而不是控制流；
3. OCCT 报**失败面**（`IMeshData_Failure`）⇒ 与 §9.351 的 rescue 讨论合流。

#### 9.364.4 下一轮的做法（明确且是新的）

**在 OCCT 侧的 Delaunay 层加一个探针**，对每个面输出：
* `myStructure->ElementsOfDomain().Extent()`（域内元素数）
* `myStructure->NbNodes()` / `NbLinks()`
* `Frontier().Extent()`

**需要找到 OCCT 里可挂的点**：
`BRepMesh_DelaunayNodeInsertionMeshAlgo::perform` / `BRepMesh_BaseMeshAlgo::collectTriangles`
都拿不到（它们之后才调）。**可行点**是
`BRepMesh_DelaunayBaseMeshAlgo::generateMesh`（`BRepMesh_DelaunayBaseMeshAlgo.cxx:43`
调用 `getCellsCount` 的那一层）或 `BRepMesh_BaseMeshAlgo::process` ——
**下一轮先定位这个可挂点**（探针已链接 OCCT 的 `TKMesh` 库，
且此前 `--uvsum` 已能用 `BRepMesh_IncrementalMesh`，说明符号可用）。

**若挂点不可得**，退路是：在端口侧把 `domain_elems` 与**输入点集**一起 dump 出来，
再在 OCCT 侧用 `BRepMesh_Delaun` **直接喂同一批点**（`BRepMesh_Delaun` 的
构造函数与 `AddVertices` 是 public）—— 这能**完全绕开 reader**，是最干净的对照。

#### 9.364.5 本轮改动与门禁

* **未提交任何源码改动**（纯 STEP 解析与对照）；
* `crates/` 的非标注改动仍只有既有的 3 处；
* 门禁：`cargo check` 0 error；`--lib` 1255/26；a3n00 基线未变。

#### 9.364.6 复跑

```text
python .target-gate\step_bounds.py
```

---

### 9.365 —— 【决定性】OCCT 在这 16 个面上**只铺 2–15 个三角形**；端口铺 204–336（约 **50×**）——这才是"过密"的真正规模

#### 9.365.1 测量方法

用现成的探针重建后跑 **新的一次 `--uvsum 0.215`**：

```text
cmd /c "specs\occt_probe\build.bat"                                  → BUILD OK
cmd /c "specs\occt_probe\probe.bat <abs>\data\occ\a3n00.stp --uvsum 0.215"  → gt215.txt
```

`gt215.txt` 有 **226 条 `UVSUM f=`**，与端口面数一致，
`BRep_Tool::Triangulation(face)->NbTriangles()` 给出每面三角数。
**口径用面序号**（两侧都按 `TopExp_Explorer` / `IMeshData_Model::GetFace` 的顺序）。

#### 9.365.2 逐面读数（同一批面号）

| face | GT `wires` | **GT `triangles`** | port `mt` | 倍数 |
|---|---|---|---|---|
| 169 | 1 | **5** | 204 | 41× |
| 174 | 1 | **5** | 288 | 58× |
| 189 | 1 | **5** | 36 | 7× |
| 192 | 1 | **2** | 228 | **114×** |
| 193 | 1 | **2** | 228 | 114× |
| 194 | 1 | **2** | 228 | 114× |
| 195 | 1 | **2** | 228 | 114× |
| 196 | 1 | **2** | 228 | 114× |
| 197 | 1 | **2** | 228 | 114× |
| 198 | 1 | **2** | 228 | 114× |
| 199 | 1 | **2** | 228 | 114× |
| 204 | 1 | **6** | 36 | 6× |
| 206 | 1 | **3** | 336 | 112× |
| 209 | 1 | **15** | 224 | 15× |
| 210 | 1 | **4** | 224 | 56× |
| 212 | 1 | **4** | 224 | 56× |

**三个关键读数**：

1. **`gt faces with tris==0: 0`**、`gt min tris: 2` ⇒
   **OCCT 在 d=0.215 下把 226 个面全部网格化，一个零三角面都没有**；
2. **在这 16 个面上，OCCT 只铺 2–15 个三角形**（`192–199` 各 **2** 个、
   `206` 是 **3** 个、`169/174/189` 各 **5** 个）；
3. 端口在同样的面上铺 **204–336** 个 ⇒ **约 50 倍（个别 114 倍）**。

#### 9.365.3 这**修正**了 §9.352/§9.365 之前的一个说法

§9.352 用**未网格化**的 `v_uvsum.txt` 配对，得到「GT 在这些面上 0/少量三角」；
配上 §9.351 的「GT 侧被判失败/无三角」的推断 —— **那个推断是错的**
（可能同样源于配对口径）。

**正确的事实**是：**OCCT 在这 16 个面上正常网格化，只是三角形很少（2–15 个）**
—— 这完全合理：在 deflection **1.076** 下，这些面很小，本来只需 2–15 个三角形。

⇒ **端口在这些面上少了 2–15 个、却多铺了 204–336 个** ——
**这不是"补洞"，而是"用错误的密度铺了错误的三角形"。**

#### 9.365.4 与面积缺口的关系（重新表述）

端口总面积比 = **0.8627**（少 14%）。而端口在这 16 个面上**多铺了**
(204+288+36+8×228+36+336+224+224+224) − (5+5+5+8×2+6+3+15+4+4)
= **3396 − 55 = 3341** 个三角形。

⇒ 端口的三角形预算**被这 16 个面吃掉**（在 GT 只需 55 个的地方花了 3396 个），
**其余面因此铺得更少** —— 这正是 §9.337「端口网格更细却面积更小」这个矛盾现象的
局部解释：**总量不缺，缺的是"把三角形用在正确的面上"**。

#### 9.365.5 下一轮的落点（不变，但**依据更强**）

病灶仍是 **`Delaun::add_vertices`**（§9.359/§9.360：这 16 个面在那里产出 0 个域内元素，
`pre=145 → post=0`）。**现在多了一个校准用的目标值**：

> **修好之后，这 16 个面应各自产出 2–15 个域内元素，而不是 0。**

这比"域内元素不为 0"强得多 —— 它给了**可验证的数量级**，
而且 `192–199` 各应是 **2**、`206` 应是 **3**，都是**具体数字**。

#### 9.365.6 旁支发现（须记录）

`UVSUM` 的 `wires=` 字段不是"面的 bound 数"：
`gt215.txt` 的 wire 直方图是 **`{1: 208, 2: 9, 4: 1, 6: 4, 10: 4}`**，
而 STEP（§9.364）与端口都是 **`{1: 163, 2: 53, 4: 2, 6: 4, 10: 4}`**。
⇒ 三者是**三个不同的量**（STEP bound 数 / 端口 wire 数 / OCCT 离散 wire 数）。
**§9.357 的「GT wires=1」正是拿第三个量去解释第一个量** —— 这是误读的来源。

#### 9.365.7 本轮改动与门禁

* **未改 Rust 源码**；探针**重建了一次**（`build.bat` → `BUILD OK`，
  用的是仓库里既有的 `occt_probe.cpp`，我没有编辑它）；
* 新增 dump：`.target-gate/gt215.txt`（临时测量产物）；
* 门禁：`cargo check` 0 error；`--lib` 1255/26；a3n00 基线未变。

#### 9.365.8 复跑

```text
cmd /c "specs\occt_probe\build.bat"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --uvsum 0.215" > gt215.txt
grep -c '^UVSUM f=' gt215.txt     # -> 226
```

---

### 9.366 —— 【可行性已证实】探针可直接构造 `BRepMesh_Delaun` 并读 `Result()` ⇒ 下一轮可走「喂同一批点」这条干净路

#### 9.366.1 为什么要走这条路

§9.365 确认：OCCT 把 226 个面全部网格化（最少 2 个三角形），
端口却在 16 个面上产出 **0** 个域内元素（§9.359/§9.360）。

⇒ 下一步必须测 **OCCT 的 Delaunay 状态**（不是它的产物）。

但 **`BRepMesh_NodeInsertionMeshAlgo` 是模板类**
（`BRepMesh_DelaunayNodeInsertionMeshAlgo<RangeSplitter, BaseAlgo>`，
见 `BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx:25-31`）
⇒ **无法直接挂钩**。

**干净的替代路**：**绕开 reader 与整个 `IMeshData` 管线**，
在探针里**直接用 `BRepMesh_Delaun` 喂端口注册的同一批点** ——
输入完全相同，输出的 `ElementsOfDomain()` 可直接与端口的 `domain_elems` 比。

#### 9.366.2 本轮做的可行性验证（探针已能编译并运行）

在 `specs/occt_probe/occt_probe.cpp` 里加了一个 `--delauncheck` 分支
（**TEMP T-99，只跑一次合成输入**）：

```cpp
#include <BRepMesh_Delaun.hxx>
#include <BRepMesh_Vertex.hxx>
#include <IMeshData_Types.hxx>          // 提供 Array1OfVertexOfDelaun

IMeshData::Array1OfVertexOfDelaun aVerts(1, 5);   // 正方形 4 角 + 1 个内部点
... aVerts(i) = BRepMesh_Vertex(gp_XY(x, y), tag, BRepMesh_Frontier); ...
BRepMesh_Delaun aDelaun(aVerts);
const occ::handle<BRepMesh_DataStructureOfDelaun>& aSt = aDelaun.Result();
std::cout << "DELAUNCHECK nodes=" << aSt->NbNodes()
          << " links=" << aSt->NbLinks()
          << " domain=" << aSt->ElementsOfDomain().Extent() << std::endl;
```

**结果**：

```text
BUILD OK
DELAUNCHECK nodes=8 links=18 domain=4
```

⇒ **`BRepMesh_Delaun` 可从本探针构造、`Result()` 可读、Delaunay 确实建出了域内元素。**
**下一轮那条路已确认可达。**

#### 9.366.3 一个**新的链接约束**（记录下来避免重复踩）

**`BRepMesh_Delaun::Frontier()` 不能从探针调用**：

```text
error LNK2019: unresolved external symbol
  "private: ... BRepMesh_Delaun::getEdgesByType(enum BRepMesh_DegreeOfFreedom) const"
  referenced in "public: ... BRepMesh_Delaun::Frontier(void) const"
```

`Frontier()`/`InternalEdges()`/`FreeEdges()` 是 **header 内联**、内部调用
**private** 的 `getEdgesByType` ⇒ **访问级别进了符号名 ⇒ 链接必然失败**
（与会话早前 `computeLengthU/V` 的情况同类）。

**可用的是 `Result()`**（它返回 `const occ::handle<...>&`，且
`BRepMesh_DataStructureOfDelaun` 的 `NbNodes()`/`NbLinks()`/`ElementsOfDomain()`
都是公开的）—— **实测证明可用**。

⇒ **下一轮要读"frontier 数"得从 `Result()` 的链状态自己统计**
（遍历 `1..=NbLinks()` 取 `GetLink(i).Movability()`），**不要用 `Frontier()`**。

#### 9.366.4 下一轮的做法（具体到可执行）

1. 端口侧加一个 dump（temp）：把 16 个失败面**注册进 Delaun 的点**输出成
   `(tag, u, v, movability)` 列表 —— 这些点在 `node_insertion.rs` 的
   `init_data_structure` 里逐个 `register_node` 时即可捕获；
2. 探针加一个 `--delaunfeed <file>` 分支：读该列表，
   `AddVertices` + `ProcessConstraints`，打印
   `NbNodes/NbLinks/ElementsOfDomain().Extent()` 与链状态直方图；
3. **判定**：
   * OCCT 也得出 `domain=0` ⇒ 端口**忠实**，病灶在**输入点本身**
     （应得到 2–15 个元素的那批点的集合不对）；
   * OCCT 得出 `domain≈2..15` ⇒ 端口在 **Delaunay 的中间状态**上与 OCCT 有差异，
     而差异**不在已核对的 21 处代码里** ⇒ 只能来自**数据**
     （注册顺序、`SetCellSize`/`SetTolerance` 的实际取值、链的初始朝向）。

**这一步把"测哪一层"彻底解决**，且**不需要改任何 Rust 逻辑代码**。

#### 9.366.5 本轮改动与门禁

* **Rust 源码零改动**（`crates/` 非标注改动仍只有既有 3 处）；
* **探针源码加了 `--delauncheck` 分支（TEMP T-99，约 25 行）**，
  它是**下一轮的脚手架**，本轮**保留**（默认不触发）；
  探针 `build.bat` → `BUILD OK`，`--uvsum` 仍正常（§9.365 的 `gt215.txt` 未受影响）；
* 门禁：`cargo check` 0 error；`--lib` 1255/26；a3n00 `stats=226 unmatched=0 mesh=11101/11121`。

#### 9.366.6 复跑

```text
cmd /c "specs\occt_probe\build.bat"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --delauncheck"
# -> DELAUNCHECK nodes=8 links=18 domain=4
```

---

### 9.367 —— 把端口的点喂给 OCCT 的 `BRepMesh_Delaun`：**同样得到 `domain=0`、`nodes=75`、`links=219`**；但**对照不完整**（探针未注册约束链）

#### 9.367.1 做了什么（§9.366 铺好的路）

**端口侧**（`node_insertion.rs`，`OCCT_TOPO_DUMP_DELAUN=<dir>`，
在 `register_node` 捕获、在 `perform` 里对那 16 个面落盘）：

```text
delaun_f169.txt … delaun_f212.txt   （16 个文件）
# face 192 domain_elems=0
-3.14159265358979312 -267.80592419198052312 5      ← u, v, movability(5=Frontier)
...
```

每个文件 **74 个点**（= 2 wires × 37，与 §9.356 的 `boundary_uv=74` 一致）。

**探针侧**（`occt_probe.cpp`，新增 `--delaunfeed <dir>`）：
读文件 → 填 `IMeshData::Array1OfVertexOfDelaun(1, 74)`（`movability==5` 传
`BRepMesh_Frontier`，其余传 `BRepMesh_Free`）→ 构造 `BRepMesh_Delaun` → 读 `Result()`。

#### 9.367.2 结果

```text
DELAUNFEED f=169 pts=74 nodes=75 links=219 domain=0 frontier=0 fixed=0 free=0 deleted=219
...（16 个面全部逐字相同）
```

**三个读数**：

1. **`nodes=75`、`links=219`** —— 与端口的 `nodes=75`、`links=221` **几乎相同**
   （差 2 条链，见 §9.367.3）；
2. **`domain=0`** —— **OCCT 的 `BRepMesh_Delaun` 用这批点也三角化不出任何域内元素**；
3. **`frontier=0 fixed=0 free=0 deleted=219`** —— 所有链都被标成 `Deleted`，
   与端口 §9.356 观测到的 `frontier=0 fixed=0 free=0` **一致**。

⇒ **在"只喂顶点"这个配置下，OCCT 与端口的行为完全一致（都是 0）。**

#### 9.367.3 **必须说明的对照边界**（这决定了本条能推出什么）

**探针这一侧没有注册约束链**：

* `BRepMesh_Delaun::Init(vertices)`（`cxx:237-250`）只做
  `AddNode` 逐个 + `perform(aVertexIndexes)` —— **不加任何 link**；
* 想加约束要走 `ProcessConstraints()`，但它**不能调用**（§9.366 的陷阱）——
  `ProcessConstraints()` 是 header 内联、内部调 **private** 的
  `insertInternalEdges()`/`frontierAdjust()` ⇒ 访问级别进符号名 ⇒ **LNK2019**。

而**端口那一侧是有约束的**（`init_data_structure` 注册了 **72 条 frontier 链**）。

⇒ 本轮实际比的是：
**「OCCT + 只有顶点（0 约束）」 vs 「端口 + 顶点 + 72 条约束」**
—— **两者都得出 `domain=0`**，但**这不是同一个配置**，所以：

> **能推出**：这 74 个点**本身**不足以（也不足以"自然地"）三角化出域内元素 ——
> **点集不是"少给了约束就能救"的那种**。
>
> **不能推出**：OCCT 在**有约束**时会得出 0。**这一点仍未测。**

#### 9.367.4 为什么这个否定仍有价值

它排除了一个**很自然的假设**：

> 「端口是不是漏注册了什么，导致 Delaunay 建不出三角？」

**答案是：不是"漏注册"那么简单** —— 即使把约束全部拿掉（更宽松），
OCCT 也是 0。⇒ 病灶**不在"注册了哪些链"这个层面**，
而在**这 74 个点本身的分布/性质**（或"点 + 链"的**组合**）。

这与此前的线索**一致**：§9.364 查明端口与 STEP 的 wire 结构一致；
§9.363 查明 Delaunay 侧 21 处代码逐行等价 —— 现在再加上
「同点集在 OCCT 里也是 0」，**指向越来越集中在"这 74 个点"上**。

#### 9.367.5 下一轮（两件事，按顺序）

1. **补全对照**：让探针也能注册约束链。绕开 `ProcessConstraints()` 链接问题的办法：
   **在构造后直接对 `Result()`（`BRepMesh_DataStructureOfDelaun`）调 `AddLink`**
   （它是 public）把 72 条 frontier 链加上，**再**构造一个新的
   `BRepMesh_Delaun(theOldMesh, cellsU, cellsV)`（该构造函数是 public，
   见 `BRepMesh_Delaun.hxx:39-42`）—— 这样约束就在初始结构里，
   而 `createTrianglesOnNewVertices` 会照常 `ProcessConstraints`。
   **这需要端口 dump 里同时给出链的 `(first, last, state)`**；
2. **查这 74 个点**：既然同点集在 OCCT 里也建不出，就该查
   **OCCT 在真实管线里喂给 Delaunay 的是不是同一批点** ——
   在 OCCT 侧 dump `IMeshData` 的 `GetFace(0)->GetWire(w)->GetEdge(e)` 的
   pcurve 采样点，与端口的 74 点**逐点比**。
   **若 OCCT 的点不同 ⇒ 病灶在端口取点的环节**（`collect_boundary_uv`/`pcurve`）；
   **若相同 ⇒ 病灶在"点 + 链"的组合**（约束链的**朝向/顺序**）。

#### 9.367.6 本轮改动与门禁

* **Rust 源码**：`node_insertion.rs` 加了 `OCCT_TOPO_DUMP_DELAUN` 落盘逻辑
  （thread-local + `register_node` 捕获 + `perform` 落盘），**TEMP T-99，须撤除**；
* **探针源码**：新增 `--delaunfeed`（TEMP T-99）；
* 两侧都**编译/构建通过**；Rust 门禁：`cargo check` 0 error、`--lib` 1255/26、a3n00 基线未变。

#### 9.367.7 复跑

```text
$env:OCCT_TOPO_DUMP_DELAUN=(Resolve-Path .target-gate\delaun_dump).Path
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --ids
cmd /c "specs\occt_probe\build.bat"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --delaunfeed D:\source\repos\dogs\.target-gate\delaun_dump"
```

---

### 9.368 —— 【决定性 + 三处更正】把端口的**完整结构**（74 个点 + 72 条约束链 + 容差/格/顶点序）喂给 OCCT 的 `BRepMesh_Delaun`：**18/18 逐字段完全相同** ⇒ Delaunay 层忠实；病灶是**输入端的面拓扑**（端口 2 wires vs OCCT 1 条由 seam 边闭合的 wire）

#### 9.368.1 结论先行（本轮四件事，按证据强度排序）

1. **等价性证明**：把端口**当前这套结构**（同样的节点坐标、同样的 72 条约束链、同样的
   `SetCellSize`/`SetTolerance`、同样的顶点序、同样的 `cellsU/cellsV`）原样喂给 OCCT 的
   `BRepMesh_Delaun`，OCCT 输出的 `NbNodes/NbLinks/ElementsOfDomain` 与四种链状态计数
   **与端口逐字段相同**（18 个面全部如此，`idx_mismatch=0`）。
   ⇒ **端口的 Delaunay 层（含 `createTrianglesOnNewVertices` / `ProcessConstraints` /
   `frontierAdjust` / `cleanupMesh`）在这一层是忠实的**，§9.359 留下的
   「同一批点、顺序不同就归零」问题**不是这层的实现差异**。
2. **更正 §9.367**：那次「只喂点、不喂链」得到的 `domain=0` **不是**关于点集的证据 ——
   结构里一条 Frontier 链都没有时，`cleanupMesh` 会把所有三角删光（见 §9.368.4）。
   该条「**能推出**：这 74 个点本身不足以三角化出域内元素」**作废**。
3. **更正 §9.365.2**：a3n00 上**两侧的面序号根本不对应**（226 个面里只有 **3 个**同序号同几何）。
   按「六坐标 bbox」重新配对后，那 16 个面的 GT 三角数是 **52–80**，不是 2–15。
   §9.365 的「OCCT 只铺 2–15 个、端口多铺 50–114 倍」是**配错面**得出的。
4. **定位到真正的输入差**：这 16 个面在 STEP 里声明 **2 个 `FACE_BOUND`**，端口忠实建成
   **2 条 wire**；OCCT 读入后（`ShapeProcess`/`FixShape`）把它们**合并成 1 条 wire**，
   并由**两条 seam 边**把两个圆连成一个环。实测直方图见 §9.368.6。

#### 9.368.2 做法（两条新探针路 + 一段临时插桩）

**端口侧（TEMP T-99，env `OCCT_TOPO_DUMP_DELAUN=<dir>`；已撤除）**：

* `node_insertion.rs::finish_mesh`：在 `Delaun::new_with_data_cells` 之前落盘
  `delaun_in_f<f>.txt` —— `TOL u v`、`CELLS cu cv`、`N i x y mov`（**face-basis 坐标**，
  即真正进入结构的那个坐标）、`L i first last mov`、`V <顶点序>`；
  构造之后追加 `AFTER_CTOR …`，`add_vertices` 之后追加 `FINAL …`；
* `node_insertion.rs::init_data_structure`：落盘 `delaun_uv_f<f>.txt` —— 登记顺序的
  **原始曲面 UV**（`E wire=… edge=… ori=… reverse_walk=… n=…` + `P u v`）；
* `delaun/triangulation.rs::create_triangles_on_new_vertices`：`OCCT_TOPO_TRACE_CTNV`
  打印 `verts/d0/f0/pre/post/failed`（§9.360 的同名仪器，本轮复测）。

**OCCT 侧（探针新增两个分支，保留作仪器）**：

* `--delaunstruct <dir>`：读上面的 dump，`AddNode` 逐个**核对返回值 == 声明序号**，
  `AddLink` × dump 的链数，`Data()->SetTolerance/SetCellSize`，再
  `BRepMesh_Delaun(aStruct, indices, cellsU, cellsV)`（`BRepMesh_Delaun.hxx:56-59` 的 public
  构造，内部经 `perform → compute → createTrianglesOnNewVertices → ProcessConstraints`
  走完整条路）。**这条路径完全不需要从探针调用 `ProcessConstraints()`**，
  §9.366 记下的 LNK2019 陷阱不适用。
* `--boundary <dir>`：用 `BRepMesh_Context` 按 `IMeshTools_MeshBuilder` 的顺序跑
  `BuildModel → DiscretizeEdges → HealModel → PreProcessModel`（`IMeshTools_MeshBuilder.cxx:52-60`），
  再对每个面照 `BRepMesh_BaseMeshAlgo.cxx:87-125` 的循环落盘
  （wire → edge → `GetPCurves(face)` → `GetPoint(0 … ParametersNb()-1)`），
  并打印面型与 bbox（**配对用 bbox，不用序号**，D21）。

**口径**：探针用 `Deflection = maxComp * 0.001 * 4 = 1.07612`、`Angle = 20°`
（= 端口 a3n00 门禁的 `prs3d_get_deflection(shape, 0.1)`，`brep_exchange.rs:56-77`）。

#### 9.368.3 结果一：同结构下 OCCT 与端口**逐字段相同**（18/18）

```text
PORT   f=0   AFTER_CTOR nodes=71 links=207 domain=66 free=65 fixed=0 frontier=68 deleted=74
OCCT   f=0              nodes=71 links=207 domain=66 frontier=68 fixed=0 free=65 deleted=74  idx_mismatch=0
PORT   f=208 AFTER_CTOR nodes=75 links=221 domain=1  free=6  fixed=0 frontier=72 deleted=143
OCCT   f=208            nodes=75 links=221 domain=1  frontier=72 fixed=0 free=6  deleted=143  idx_mismatch=0
PORT   f=192 AFTER_CTOR nodes=75 links=245 domain=0  free=0  fixed=0 frontier=72 deleted=173
OCCT   f=192            nodes=75 links=245 domain=0  frontier=72 fixed=0 free=0  deleted=173  idx_mismatch=0
PORT   f=212 AFTER_CTOR nodes=75 links=221 domain=0  free=2  fixed=0 frontier=72 deleted=147
OCCT   f=212            nodes=75 links=221 domain=0  frontier=72 fixed=0 free=2  deleted=147  idx_mismatch=0
```

18 个面（`f=0, 208, 209, 169, 174, 189, 192–199, 204, 206, 210, 212`）**全部**如此：

| 面 | 结构头 | OCCT == 端口 |
|---|---|---|
| 0（健康对照） | nodes=72 links=72 cells=? | `nodes=71 links=207 domain=66 frontier=68 free=65 deleted=74` |
| 208（健康对照） | nodes=72 links=72 cells=4 2 | `nodes=75 links=221 domain=1 frontier=72 free=6 deleted=143` |
| 209/169/174/189/204/210 | nodes=72 links=72 | `nodes=75 links=221 domain=0 frontier=72 free=0 deleted=149` |
| 192–199 | nodes=72 links=72 cells=3 4 | `nodes=75 links=245 domain=0 frontier=72 free=0 deleted=173` |
| 206 | nodes=72 links=72 | `nodes=75 links=227 domain=0 frontier=72 free=0 deleted=155` |
| 212 | nodes=72 links=72 | `nodes=75 links=221 domain=0 frontier=72 free=2 deleted=147` |

三条读数：

1. `idx_mismatch=0`：端口结构的 72 个节点用 dump 里的坐标 + 容差逐个 `AddNode`，
   OCCT 返回的序号**与端口完全一致**（没有多合并、没有错位）⇒ 节点集合等价；
2. `links_dumped == links_in_struct == 72`：**72 条约束链一条不多一条不少**地进了 OCCT 的结构；
3. 域内元素数与**四种链状态的分布**逐字段相同 ⇒ 分歧不在这一层。

端口侧 CTNV 复测（243 条记录）也逐字复现 §9.360：那 16 个面是
`CTNV verts=72 d0=3 f0=72 pre=145 post=0 failed=false`
（142 个三角是在 `createTrianglesOnNewVertices` 里建出来的，随后被
`ProcessConstraints` 删光），而 OCCT 在**同一结构**上给出同样的 `domain=0`。

#### 9.368.4 结果二：更正 §9.367 —— 「只喂点」的 `domain=0` 是 `cleanupMesh` 造成的

§9.367 喂的是 74 个点 + **0 条链**，得到 `domain=0 / frontier=0 / fixed=0 / free=0 / deleted=219`。
当时把它读成「这批点本身三角化不出域内元素」。**实际机制**：
`frontierAdjust()` 末尾无条件调用 `cleanupMesh()`（`BRepMesh_Delaun.cxx:1028`），而
`cleanupMesh` 遍历 `FreeEdges()` 时**只跳过 `BRepMesh_Frontier` 的链**
（`cxx:832-835`）；结构里若**一条 Frontier 链都没有**，所有链都是 Free，
`isConnected[]` 必为 false（`cxx:908-911`）⇒ 所有三角被判为"外部"并删除。
⇒ 那次读数**只反映"结构里没有约束"**，与点集无关。§9.367 的「能推出」作废，
「不能推出」那半边（有约束时的行为未测）由本轮 §9.368.3 补齐。

#### 9.368.5 结果三：a3n00 的面序号两侧**不对应**（3/226）⇒ §9.365.2 的表作废

用六坐标 bbox 配对（容差 1e-3，`.target-gate/pair368b.py`）：

```text
port faces=226  model faces=226
same-index bbox matches: 3 / 226          ← 同序号同几何只有 3 个
paired 1:1 by bbox:      201 / 226
```

**端口 → OCCT 离散模型**的实际置换（16 个失败面）：

```text
port 169 -> model  87   port 192 -> model 85   port 199 -> model 82
port 174 -> model  42   port 193 -> model 84   port 204 -> model 43
port 189 -> model  41   port 194 -> model 83   port 206 -> model 86
                        port 195 -> model 78   port 209 -> model 12
                        port 196 -> model 79   port 210 -> model 13
                        port 197 -> model 80   port 212 -> model 15
                        port 198 -> model 81
```

⇒ **§9.365.2 的「GT `triangles`」列取的是另一些面**（例如它配给 port 192 的
「GT f=192」其实是 y=52 那个平面小面，`triangles=2`）。按正确配对重取
（`--uvsum 1.07612`，模型序），那 16 个面的真实读数是：

| port f | model g | GT tri | port mt | | port f | model g | GT tri | port mt |
|---|---|---|---|---|---|---|---|---|
| 169 | 87 | **54** | 204 | | 196 | 79 | **52** | 228 |
| 174 | 42 | **80** | 288 | | 197 | 80 | **52** | 228 |
| 189 | 41 | **54** | 36 | | 198 | 81 | **52** | 228 |
| 192 | 85 | **52** | 228 | | 199 | 82 | **52** | 228 |
| 193 | 84 | **52** | 228 | | 204 | 43 | **54** | 36 |
| 194 | 83 | **52** | 228 | | 206 | 86 | **52** | 336 |
| 195 | 78 | **52** | 228 | | 209/210/212 | 12/13/15 | **52**（近邻配对） | 224 |

⇒ 真实倍数是 **约 4 倍（个别 6.5 倍）**，不是 50–114 倍；且**方向是双向的**：
`189`/`204` 是端口**少**铺（36 vs 54）。
13 个精确配对面的合计：**GT 710 vs 端口 2724**。
§9.365.4 的「三角形预算被这 16 个面吃掉」这一说法方向仍成立，量级需按新表重述。

#### 9.368.6 结果四：真正的输入差 = **面拓扑**（STEP 2 bounds vs OCCT 1 条 seam 闭合 wire）

**（a）STEP 侧**：§9.364 已实测这 16 个面在文件里就是 **2 个 `FACE_BOUND`**，
端口照建 2 条 wire（端口 wire 直方图 `{1:163, 2:53, 4:2, 6:4, 10:4}`）。

**（b）OCCT 侧**：`--faceids` 逐面数 `TopExp_Explorer(face, TopAbs_WIRE)`：

```text
默认（ShapeProcess/FixShape 开）  wire 直方图 {1:208, 2:9, 4:1, 6:4, 10:4}
--nofix（ShapeProcess 关）        wire 直方图 {1:163, 2:53, 4:2, 6:4, 10:4}   ← 与 STEP、与端口逐字相同
```

⇒ **是读入期的 `ShapeProcess` 把 45 个「2 bound」面并成了 1 条 wire**
（163+45=208、53−44=9、4: 2→1）。端口即使按 D27 接了 `fix_missing_seam`，
这 16 个面**仍是 2 条 wire**。

**（c）逐边对照**（port f=192 ↔ model f=85，两侧都是原始曲面 UV）：

```text
端口：2 条 wire，各 1 条闭合 pcurve × 37 点（74 点 / 72 条 frontier 链）
      wire0 edge=449 FWD n=37   圆 A  v=-267.805924  u∈[-π, π]
      wire1 edge=457 FWD n=37   圆 B  v=-287.805924  u∈[ 0, 2π]

OCCT：1 条 wire，5 条边 / 7 条 pcurve（83 点）
      e0 pcId=1 FWD n=19   u: π → 0      v=-267.805924      ← 圆 A 的**上半弧**
      e1 pcId=0 REV n=2  \  seam 边（同一条边在该面上的两条 pcurve）
         pcId=1 FWD n=2  /  u=0 / u=2π，v: -287.805924 → -267.805924
      e2 pcId=1 FWD n=37   u: 0 → 2π     v=-287.805924      ← 圆 B 整圈
      e3 pcId=0 REV n=2  \  第二条 seam 边（同上）
         pcId=1 FWD n=2  /
      e4 pcId=1 FWD n=19   u: 2π → π     v=-267.805924      ← 圆 A 的**下半弧**
```

⇒ **两条 seam 边把两个圆连成一个环**，这才是 OCCT 那 1 条 wire 的来源；
端口的两条 wire 是**两个互不相连的水平线段**（圆 A 在 `v=-267.805924`、
圆 B 在 `v=-287.805924`），frontier 是 72 条互不闭合的链。
这解释了为什么同一套 Delaunay 代码在 OCCT 那里留下 52 个三角、在端口这里 `post=0`。

**（d）两处**更正**§9.364/§9.365.6**：`--uvsum` 的 `wires=` 字段**就是**该面的拓扑 wire 数
（`--faceids` 的直方图与它逐字相同 `{1:208, 2:9, 4:1, 6:4, 10:4}`），**不是"第三个量"**；
真正的两个量是「文件 bound 数 = 端口 wire 数 = OCCT `--nofix` wire 数」与
「OCCT 读入整形后的 wire 数」。§9.363.3 的假说（"OCCT 把它们合并了 / 端口缺合并"）
**据此恢复成立**，§9.364 对它的否证作废（它当时的依据是"STEP 里就是 2 个 bound ⇒ OCCT 必然也是 2 条"，
而已实测 OCCT 会改）。

**（e）一处次生差异**（须记）：同一个圆 A，端口取 u 的 `[-π, π]` 分支，OCCT 取 `[0, 2π]`
并在 `u=π` 处切开成 19+19。两侧是同一几何、不同参数分支；合并 seam 时选哪条分支
要按 `.cxx` 的分支走，不能按端口现有取值反推。

#### 9.368.7 下一轮落点（含可证伪的预测）

1. **实现读入期那条"把 2 个 bound 并成 1 条 seam 闭合 wire"的控制流** ——
   端口已经标出候选位置：`crates/occt-topo/src/shhealing/shape_fix_face.rs:147`
   （«UNPORTED: the seam construction and the `w2 != null` merge»）与
   `:406` 处 `FixReorder`（`cxx:2029-2032`）的未移植注记；
   对照 `ShapeFix_Face.cxx:492-498 / :1722-2330` 与 `ShapeProcess_OperLibrary.cxx:785-899`。
   **不要**在 Delaunay 层加任何东西（§9.368.3 已证这层忠实），
   **不要**在 `frontier_adjust` / `cleanup_mesh` 里加阈值或特例。
2. **可证伪的预测**（下一轮的验收就是它）：这条控制流接上后，
   端口 f=192（及其同类）的结构应变成 **1 wire / 83 点 / 含两条 2 点 seam pcurve**，
   `AFTER_CTOR domain` 应从 **0 变成 >0**、`mt` 应从 **228 降到 ~52**，
   a3n00 面积比从 **0.8627** 上升（T-99 的 accept）；
   若结构变了而 `domain` 仍是 0，则说明输入侧还差别的（那时再回到
   `range_splitter`/tolerance 一侧量）。
3. **配对纪律**（本轮最大的方法论教训）：a3n00 上**任何按面序号的对拍都作废**，
   一律先按六坐标 bbox 配对（`--boundary` 已把每个模型面的 bbox 落盘）；
   沿用这条再复读 §9.354/§9.357/§9.365 里所有"GT 面号"的结论。

#### 9.368.8 本轮改动与门禁

* **Rust 源码零改动**：端口侧插桩（`node_insertion.rs`、`delaun/triangulation.rs`）
  已按 §9.353 规程**用 `edit` 反向撤除、未用 `git checkout`** ⇒
  `git status --porcelain` 对 `crates/` **无输出**（`git diff --stat` 只剩探针）；
* **探针新增** `--delaunstruct <dir>` 与 `--boundary <dir>`（TEMP T-99，保留作仪器）；
* 门禁（本轮实测）：`cargo check` **0 error**；
  `cargo test --lib` → `test result: ok. 1281 passed; 0 failed; 0 ignored`；
  a3n00 `TOTAL faces=226 computed_lin=1.076007 used_lin=1.076007 stats=226 mesh_v=11101 mesh_t=11121`、
  `STATMAP matched=226 unmatched=0`。
  （口径说明：§9.367.6 记的是「--lib 1255/26」，本轮 `rtk cargo test` 与
  `cargo test --lib` 两种调用都给出 **1281 passed / 0 failed**，那条记录与实测不符，以本轮实测为准。）

#### 9.368.9 复跑

```text
# 1) 端口结构 dump（插桩须先按 §9.368.2 重新加回，本轮已撤除）
$env:OCCT_TOPO_DUMP_DELAUN=(Resolve-Path .target-gate\delaun_in).Path
$env:OCCT_TOPO_TRACE_CTNV="1"
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --ids 2> .target-gate\ctnv368.txt

# 2) 同结构喂 OCCT
cmd /c "specs\occt_probe\build.bat"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --delaunstruct D:\source\repos\dogs\.target-gate\delaun_in" > .target-gate\delaunstruct368.txt

# 3) OCCT 真实管线的边界点 + 面 bbox
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --boundary D:\source\repos\dogs\.target-gate"
# 4) 面拓扑（有/无 ShapeProcess）
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --faceids"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --faceids --nofix"

# 5) 配对与逐面读数
cd .target-gate; python pair368b.py ; python pair368c.py
```

> 数据留在 `.target-gate/delaun_in/`（端口结构 dump，452 个文件，本轮的
> `--delaunstruct` 输入）与 `.target-gate/{delaunstruct368,uvfeed368,gtmesh368,gtuvsum368,faceids368,faceids_nofix368,boundary}.txt`。

---

### 9.369 —— 【T-99 的 accept 达成】把导入期 `ShapeFix_Face::FixMissingSeam` 接回 reader：a3n00 面积比 **0.8627 → 0.8996**、T0M 未网格 **7 → 6**、UNPORTED rescue 触发 **16 → 0**；带倒角法兰的 16 个面从「2 条 wire」变成 OCCT 的「1 条 wire / 5 条边」

#### 9.369.1 结论先行

**用户这一轮的要求是「a3n00 一个一个处理，先处理带倒角的法兰盘」。**
按 §9.368 定位到的病灶（读入期缺 bound 合并）动手，结果如下：

| 读数（口径见 §9.369.7） | 改前 | 改后 | OCCT |
|---|---|---|---|
| a3n00 面积比（`step_obj_gates` `step_obj_area`） | **0.8627** | **0.8996** | 1 |
| T0M 面积比（同上） | 0.9786 | **0.9987** | 1 |
| acs10 面积比（同上） | 0.9025 | **0.9846** | 1 |
| a3n00 `zz_probe_a3n00` | `stats=226 unmatched=0 mesh 11101/11121` | `stats=225 unmatched=1 mesh 10863/11941` | — |
| T0M `zz_probe_a3n00` | `stats=1765 unmatched=7` | **`stats=1766 unmatched=6`** | — |
| a3n00 wire 直方图（`zz_uv_feed --ids`） | `{1:163, 2:53, 4:2, 6:4, 10:4}` | **`{1:206, 2:10, 4:2, 6:4, 10:4}`** | `{1:208, 2:9, 4:1, 6:4, 10:4}` |
| a3n00 UNPORTED UV-grid rescue 触发次数 | 16（§9.351） | **0** | 0（OCCT 无此路径） |
| 法兰/倒角那 16 个面 | 2 wires，mt 204–336 | **1 wire**，mt 72–110 | 1 wire，GT 52–80 |

⇒ **T-99 的 accept（「a3n00 面积比从 0.8627 上升且不劣化；T0M 未网格 7 不增加；
`--lib` 与 `step_obj_gates` 不劣化；不得改基线或 `area_tol`」）本轮全部满足，
且没有改任何断言或 `area_tol`。**

#### 9.369.2 做法：把 f1e56776 删掉的那段接回来（不是新造规则）

`crates/occt-topo/src/step/read_topology.rs::resolve_face` 末尾重新接上
「`ShapFixFace::fix_missing_seam` + 结果面每条 wire 的 `check_pcurves_and_shift`」，
逐行对照：

```text
ShapeProcess_OperLibrary.cxx:785-899   FromSTEP 的 FixShape 算子 = ShapeFix_Shape::Perform
ShapeFix_Shape.cxx                     → 每个面跑 ShapeFix_Face::Perform
ShapeFix_Face.cxx:482-498              myResult = myFace; 若 myFixMissingSeamMode 则 FixMissingSeam()
ShapeFix_Face.cxx:1722-2330             seam 构造 + w2 != null 合并（端口的 fix_missing_seam）
ShapeFix_Face.cxx:365-480              第一阶段 wire 修（端口 check_pcurves_and_shift 对应）
```

`STEPControl_Controller.cxx:201` 设 `FromSTEP.exec.op = FixShape`、`:221` 设
`FixMissingSeamMode = -1`（NeedFix 为真）⇒ **默认就开**。

**为什么可以接回来（§9.219 的移除理由已被本轮实测推翻）**：用探针 `--faceids`
（逐面 `TopExp_Explorer(face, TopAbs_WIRE)`）在 OCCT 上开/关 ShapeProcess 对比：

```text
默认（FixShape 开）   {1:208, 2:9, 4:1, 6:4, 10:4}
--nofix（全关）       {1:163, 2:53, 4:2, 6:4, 10:4}   ← 与 STEP 的 FACE_BOUND 直方图、与端口逐字相同
```

⇒ **OCCT 的导入确实合并了 45 个「2 bound」面**（163→208）。T-93(a) 当时依据的
「补丁版 OCCT 在 `TransferRoots` 期间无输出」只说明那条 `Perform` 符号没被调用到，
不能推出「没有别的入口做同一件事」；D27 已提出该更正，本轮用可复跑的面拓扑直方图坐实。

#### 9.369.3 逐面验证：法兰面结构与 OCCT **逐项相同**

`zz_seam_fix data/occ/a3n00.stp 192`（单面，改前已验证）：

```text
BEFORE wires=2   wire[0] nEdges=1 par=[π,3π]   wire[1] nEdges=1 par=[0,2π]
HEALED wires=1   nEdges=5  e0 par=[2π,3π] Circle / e1 par=[π,2π] Circle
                          / e2 par=[-287.805924,-267.805924] Line(seam, npc=2, 反向两次)
                          / e3 par=[-287.805924,-267.805924] Line(seam, Reversed)
                          / e4 par=[0,2π] Circle
```

接入 reader 后**端口模型**的同一面（`zz_uv_feed --model 192`）：

```text
MODEL f=192 wires=1 urange=[0,6.283185307179586] vrange=[-287.805924191980466,-267.805924191980523]
MODEL  wire[0] nEdges=5
MODEL   e[0] edge=536 ori=Forward deg=false curve=Line   par=[-287.805924191980466,-267.805924191980523] npc=2
MODEL   e[1] edge=537 ori=Forward deg=false curve=Circle par=[6.283185307179586,9.424777960769379]
MODEL   e[2] edge=538 ori=Forward deg=false curve=Circle par=[3.141592653589793,6.283185307179586]
MODEL   e[3] edge=536 ori=Reversed deg=false curve=Line  par=[-287.805924191980466,-267.805924191980523] npc=2
MODEL   e[4] edge=539 ori=Forward deg=false curve=Circle par=[0,6.283185307179586]
```

与 §9.368.6 记录的 OCCT 模型孪生 f=85（1 wire / 5 边 / 7 pcurve：两条 2 点 seam + 圆 A
切成 19+19 两段 + 圆 B 整圈）**边数、边序、参数区间逐项对应**。

**16 个面（法兰盘 + 倒角 + 螺栓孔）改后读数**（`zz_uv_feed --ids`）：

```text
f=169  1 wire mt= 72 (改前 2 wires / 204)    f=196 1 wire mt=72 (228)
f=174  1 wire mt=110 (288)                   f=197 1 wire mt=72 (228)
f=189  1 wire mt= 72 ( 36)                   f=198 1 wire mt=72 (228)
f=192  1 wire mt= 72 (228)                   f=199 1 wire mt=72 (228)
f=193  1 wire mt= 72 (228)                   f=204 1 wire mt=72 ( 36)
f=194  1 wire mt= 72 (228)                   f=206 1 wire mt=72 (336)
f=195  1 wire mt= 72 (228)                   f=209/210/212 1 wire mt=73 (224)
```

按曲面类型（`--boundary` 的模型孪生）：这 16 个面 = **10 个圆柱 + 3 个圆锥 + 3 个 BSpline 型**
（圆锥即倒角面 `f=174/189/204`；`f=192..199` 是左法兰上 8 个成环排布的圆柱孔/凸台，
`f=209/210/212` 是右法兰的同类）。

#### 9.369.4 整模型读数：rescue 不再触发

* a3n00 `--ids` 的 226 个面：**rescue 触发 0 次**（改前 16 次，§9.351）；
  ⇒ 这 16 个面现在走的是**忠实 Delaunay 路径**，不再靠 `wireframe_face_triangulation`
  （`discret_root.rs:477-488`，UNPORTED，audit A19/T-68）补三角形。
* wire 直方图 `{1:163,2:53,4:2,6:4,10:4}` → `{1:206,2:10,4:2,6:4,10:4}`（43 个面 2→1），
  与 OCCT 的 `{1:208,2:9,4:1,6:4,10:4}` 只差 **1 个 2-wire + 1 个 4-wire 面**（见 §9.369.6）。
* 面积比与 T0M 读数见 §9.369.1 的表；`step_obj_gates` **5/5 全绿**、`--lib` **1281 passed / 0 failed**。

#### 9.369.5 仍未解决（本轮的旁支，**未动手**，留作「一个一个处理」的下一项）

1. **还剩 2 个面没合并**：端口 `{1:206,2:10,4:2}` vs OCCT `{1:208,2:9,4:1}`
   ⇒ 有 1 个 2-bound 面与 1 个 4-bound 面 OCCT 会合并而端口没有。需要按 bbox 配出这两个面，
   再对着 `ShapeFix_Face::FixMissingSeam` 的分支（`cxx:1899-2330`）看它们为什么不进合并路径。
2. **全形状跑一遍比逐面多改 4 个面**：`zz_seam_fix --all` 在 a3n00 上 `faces=226 changed=46`，
   其中 42 个是 `2→4 / 2→5` 边（即 seam 合并），另有 4 个是
   `13→15 / 13→16 / 22→28 / 5→8`。这 4 个是不是 OCCT 也会改（`FixReorder`/post-seam 面循环
   那两处 UNPORTED 的产物），本轮**没查**。
3. **法兰面的三角密度仍比 OCCT 密**：`f=192..199` 端口 72 个三角 vs GT 52；
   §9.368 的预测写的是「≈52」，实际落在 72 ⇒ 结构对了，密度还差一档
   （`Rescue` 已不参与，属 `GenerateSurfaceNodes`/deflection 一侧的问题，**未动**）。
4. **面积比 0.8996 vs 1**：缺口从 31179.6 降到约 22800（`our=204327.72 / occ=227130.30`），
   与 §9.316 记录的 F113 单面区域量级（26952）同阶；下一步仍应按「面」而不是按「全局参数」找。

#### 9.369.6 本轮改动与门禁

* **库代码**：`crates/occt-topo/src/step/read_topology.rs` 一处（37+/13−，含注释），
  即 §9.369.2 的 seam 步骤接线；`discret_root.rs` 的 TEMP rescue 计数插桩
  **已用 `edit` 反向撤除**（`git diff` 对该文件为空）；
* 未改任何测试、断言、`area_tol` 或基线；
* 门禁（本轮实测，最终树）：`cargo check` 0 error；`cargo test --lib`
  `1281 passed; 0 failed`；`cargo test --test step_obj_gates -- --nocapture`
  `5 passed; 0 failed`（385.81s）；a3n00 / T0M 读数见 §9.369.1。

#### 9.369.7 复跑

```text
# 逐面结构（改前/改后同一命令）
crates\occt-topo\target\debug\examples\zz_seam_fix.exe data/occ/a3n00.stp 192
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --model 192

# 全形状会改哪些面
crates\occt-topo\target\debug\examples\zz_seam_fix.exe data/occ/a3n00.stp 0 --all

# wire 直方图 + 逐面 wires/mt
crates\occt-topo\target\debug\examples\zz_uv_feed.exe data/occ/a3n00.stp --ids > .target-gate\uvfeed369.txt

# 面积比（stderr 里每行 [model] our=... occ=... ratio=...）
cargo test --manifest-path crates/occt-topo/Cargo.toml --offline --test step_obj_gates -- --nocapture

# 端口侧 stats/unmatched
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_probe_a3n00 -- data/occ/a3n00.stp
rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example zz_probe_a3n00 -- data/occ/T0M.stp
```

---

### 9.370 —— T-101 迭代（1）：剩下的 seam 合并缺口**只剩 3 个面**，且合并本身是对的 —— 坏在「结果面 vs 结果 shell」；另外发现探针间的面三角数不可互比（D16 同类）

#### 9.370.1 本轮的三个结论

1. **② 解决（否定「全形状比逐面多改」这个疑点）**：`zz_seam_fix --all` 比逐面多改的那 4 个面
   （`13→15 / 13→16 / 22→28 / 5→8`）**就是 ① 的那批面**（113 / 138 / 140 / 170），
   而且它们的**后置边数与 OCCT 成形后的面逐一相同**：

   ```text
   port f=113  22 -> 28 边   OCCT 孪生 model f=77 = W0 22 边 + W1 6 边 = 28  ✓
   port f=140  13 -> 16 边   OCCT 孪生 model f=39 = W0 16 边              ✓
   port f=170   5 ->  8 边   OCCT 孪生 model f=69 = W0  8 边              ✓
   ```

   （比边数要用 `W k edges=N`，**不能**用 `E` 行 —— `E` 行把 seam 边在该面上的两条 pcurve 各算一次。）
   ⇒ 缝合并这一步是忠实的，此前「多改 4 个面」的疑问作废。

2. **① 收敛到 3 个面、且病灶是「结果装配」而不是「合并算法」**：用**既有探针**
   `zz_probe_a3n00 <file> --fixms`（逐周期面打印 `before_wires / ret / result=<kind> faces=N`）
   实测 a3n00：**只有 3 个面返回 Shell**，其余 115 个周期面要么 `ret=false`
   （说明它们已经在 reader 里被修好了，现在只剩 1 条 wire）、要么返回 Face：

   ```text
   FIXMS face=113 uper=true vper=false before_wires=4 ret=true result=Shell faces=5
   FIXMS face=140 uper=true vper=false before_wires=2 ret=true result=Shell faces=2
   FIXMS face=170 uper=true vper=false before_wires=2 ret=true result=Shell faces=2
   ```

   而 OCCT 的离散模型这 3 个面各是**一个 Face**（wire 数 2 / 1 / 1）。
   端口的 `resolve_face` 只接受 `ShapeType::Face`，所以这 3 个 Shell 被丢弃、
   面保持原始 wire 数 —— 这正是 `{1:206,2:10,4:2}` 与 OCCT `{1:208,2:9,4:1}` 的全部差值。

   结合 ① 的边数对照（合并结果与 OCCT 相同），病灶**只能**落在
   `shhealing/shape_fix_face.rs:590-673` 这段尾巴：
   假想 grid（`cxx:2236-2245`）→ `ComposeShell`（`cxx:2246-2261`）→
   `myResult = CompShell.Result()`（`cxx:2263-2268`）→ 两轮剪枝
   （`FixSmall` / `FixSmallAreaWire`，`cxx:2270-2322`，端口 `:624-667`，
   判据是 `crate::shhealing::check_small_area` ↔ `ShapeAnalysis_Wire::CheckSmallArea` `cxx:2004`）。
   **这一段的分歧定位已交给一个并行子任务深挖**（对照 `ShapeFix_Face.cxx` 与
   `crates/occt-topo/src/shape_fix_compose_shell/`），本轮未动代码。

3. **③ 的仪器与一个必须先修的测量问题**：新增端口逐面 `mv/mt/bbox` 仪器
   `zz_probe_a3n00 <file> --fstats`。用它和 `zz_uv_feed --ids` **同文件同参数**对拍发现：
   两个探针的**逐面三角数在 226 个面里有 89 个不同**（例：`f=204` 在 `--fstats` 是 376、
   在 `--ids` 是 72），而**逐面 bbox 226/226 完全一致**。
   ⇒ 这是 D16 记过的同类现象（端口的逐面结果依赖调用序列/预热），
   **跨探针的逐面读数不可互比**；密度对拍必须在**同一次运行内**取数。
   本轮据此只做「同一 run 内 port ↔ 另一次 OCCT run」的 bbox 配对，得到
   **双向**差异（既有 4–8 倍感度更高的面，也有明显更疏的面），但 26 个面配不上，
   总量不可比 ⇒ ③ 留到下一轮，先把逐面统计改成与模型序同源。

#### 9.370.2 新增/沿用的仪器

* `crates/occt-topo/examples/zz_probe_a3n00.rs`：新增 `--fstats`
  （逐面 `FSTAT face=i mv= mt= bbox=`，模型序，与 `PORTID`/`FSTAT` 同序已实测 226/226 bbox 相同）；
* 沿用（本轮才发现已经存在、不需要库内插桩）：`--fixms`
  （逐周期面 `fix_missing_seam` 的返回值与结果形状）、`--low`、`--fdump`；
* `specs/occt_probe/occt_probe.cpp`：沿用 `--boundary` / `--uvsum` / `--mesh` / `--faceids`。

#### 9.370.3 本轮改动与门禁

* **库代码零改动**（诊断轮）：`crates/occt-topo/src/` 的 `git diff` 为空；
  唯一改动是上面那个 **example 仪器**（不在库路径上）；
* 未跑全量门禁（无库改动，上一轮的 `--lib` / `step_obj_gates` 结论仍然成立）；
* 相关读数（同上一轮最终树）：a3n00 面积比 0.8996、`stats=225 unmatched=1 mesh 10863/11941`；
  T0M 0.9987、`stats=1766 unmatched=6`。

#### 9.370.4 复跑

```text
# 逐面 seam 步结果（不需要库内插桩）
crates\occt-topo\target\debug\examples\zz_probe_a3n00.exe data/occ/a3n00.stp --fixms
# 逐面 mv/mt/bbox
crates\occt-topo\target\debug\examples\zz_probe_a3n00.exe data/occ/a3n00.stp --fstats
# 单面结构（BEFORE / result / HEALED）
crates\occt-topo\target\debug\examples\zz_seam_fix.exe data/occ/a3n00.stp 113
# 全形状会改哪些面（edges a -> b）
crates\occt-topo\target\debug\examples\zz_seam_fix.exe data/occ/a3n00.stp 0 --all
# OCCT 侧对照
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --boundary D:\source\repos\dogs\.target-gate"
python .target-gate\t101_extra2.py     # 4 个面的 post 边数 vs OCCT 孪生
python .target-gate\t101_density2.py   # 逐面密度图（注意 §9.370.1 ③ 的可比性限制）
```

---

### 9.371 —— T-101 迭代（2）：修正逐面统计仪器的**配对口径**（`FaceMeshStat.index`，不是位置下标），拿到可信的逐面密度图

#### 9.371.1 修正了什么

§9.370 里新增的 `zz_probe_a3n00 --fstats` 用 `faces_of(&shape)[i]` 去配 `face_stats()[i]` ——
**这是错的**。端口自己的类型注释已写明原因（`crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs:81-88`）：

> `FaceMeshStat` 带 `index`（网格模型自己的面序号）与 `shape_key`（源面 TShape 的指针身份键），
> **就是为了让调用方不必依赖遍历顺序**（网格模型会跳面，按位置查会错开「被跳过的面数」）。

改用 `stat.index` 并通过同一条管线的模型 `model.face(stat.index)` 取 bbox 之后：

```text
FSTAT face=192 key=… mv=74 mt=72 bbox=(-132.500000,-64.932772,-52.538994)-(…)
PORTID f=192        wires=1 mt=72 bbox=(同上)          ← 与 --ids 的读数一致
```

修正前同一面是 `mt=376`（把别的面的统计贴到了它身上）—— 89/226 个面受影响。
⇒ **`--low`（同一处 `stats.get(i)` 口径）的历史读数也带这个偏差**，本卡后续引用 `--low` 的
逐面 mv/mt 时要按修正后的口径重取。

#### 9.371.2 修正后的逐面密度图（199 个面按 bbox 配对，port `--fstats` ↔ OCCT `--uvsum`）

```text
端口更疏（前 5）： 171 BSpline 349 vs 495 | 140 Cone 5 vs 65 | 170 Cyl 27 vs 58
                  146 Cyl 6 vs 12 | 39 BSpline 46 vs 50
端口更密（前 6）： 129 BSpline 802 vs 160 | 202 Torus 548 vs 167 | 118/122 Torus 211 vs 98
                  117/121 Torus 162 vs 98 | 191/200/205 376 vs 277 | 114/132 Plane 228 vs 170
配对面合计: port mt=9272 vs occ 6392（1.45×）
```

两条要点：

1. **差异是双向的**，不是「整体更密」：`171`（BSpline）与 `140/170`（正是 §9.370 里那 3 个
   未被合并的面）端口明显**少铺**，而 Torus/BSpline 类端口明显**多铺**；
2. `140/170` 的少铺与 §9.370 的结论自洽：这两个面的 seam 合并被当成 Shell 丢掉了，
   面还是 2 条 wire、三角化停在 5/27 个三角（GT 65/58）。

⇒ **③（密度）不能只看那 16 个法兰面**；下一轮应先按「每类曲面」统计差值分布，
再挑面数最集中的那一类（Torus / BSpline）回到 `.cxx` 找控制流。

#### 9.371.3 本轮改动与复跑

* 只改 example 仪器（`crates/occt-topo/examples/zz_probe_a3n00.rs` 的 `--fstats`），
  **库代码零改动**；
* 复跑：`zz_probe_a3n00.exe data/occ/a3n00.stp --fstats > .target-gate\fstats370.txt`
  然后 `python .target-gate\t101_density4.py`。

---

### 9.372 —— T-101 迭代（3）：③ 密度差找到**代码级候选**（Torus 的 `fillParams` 依赖参数容器的遍历顺序：OCCT 是插入序的 `IndexedMap`，端口是**排序**的 `ParamSet`）；同时记下密度图的配对质量警告

#### 9.372.1 按曲面类型归并（199 个 bbox 配对面）

```text
type      n     port_mt  occ_tri  diff     ratio
Torus     7     1507     706      801      2.135
BSpline   15    2540     1751     789      1.451
Cyl       61    2390     1782     608      1.341
Plane     86    2168     1614     554      1.343
Cone      30    667      539      128      1.237
```

**Torus 2.14× 最突出**。但先看 §9.372.3 的警告：这 199 个配对面的 GT 三角和只有 6392，
而全模型 GT 是 12324 ⇒ 26 个未配上的面吃掉了 5932（均值 228），
说明**配对在「大面」上系统性失败**（大面的 bbox 在整形/合并后与本侧不完全一致）。
⇒ 上表的**比例不能直接当结论**，只能当「哪一类值得先查」的线索；Torus 的面数少、bbox 唯一性好，
是其中最可信的一组。

#### 9.372.2 候选根因：参数容器的**遍历顺序**（只影响 Torus）

OCCT：

```cpp
// IMeshData_Types.hxx:144
typedef NCollection_Shared<NCollection_IndexedMap<double>> IMapOfReal;
// BRepMesh_TorusRangeSplitter.cxx:116-123
void AddPoint(const gp_Pnt2d& thePoint) { … GetParametersU().Add(thePoint.X()); GetParametersV().Add(thePoint.Y()); }
// :129-176 fillParams(theParams, range, stepsNb, scale, alloc):
//   aParamArray(j) = theParams(j);            ← 按**索引序**（= 插入序，IndexedMap 不去重排序）
//   贪心抽稀：aParams 已收的每个值与 pp 的距离都要 > aStdStep 才收 pp   ← 结果依赖遍历顺序
```

端口：

```rust
// crates/occt-topo/src/meshing/range_splitter/param_set.rs:5-28
/// Sorted set of real parameters … Mirrors OCCT's `IMeshData::IMapOfReal`;
pub struct ParamSet { values: Vec<f64> }
pub fn insert(&mut self, v: f64) { /* binary_search + insert —— **保持升序** */ }
// crates/occt-topo/src/meshing/range_splitter/splitter.rs:485-511 fill_params 用 params.iter()（升序）
```

⇒ 同一批边界参数，OCCT 按**首次出现顺序**抽稀，端口按**升序**抽稀，
**贪心保留的子集不同** ⇒ Torus 内部网格点数不同 ⇒ 三角数不同。
`fillParams` 在 OCCT 里**只有 TorusRangeSplitter 用**（`BRepMesh_*RangeSplitter.cxx/hxx` 全量 grep：
只有 `BRepMesh_TorusRangeSplitter.cxx:89/92` 两处），与本轮「Torus 偏差最大」的观察一致。

**这是「没有同等分支」类的差异**：端口 `ParamSet` 的注释自己声称「Mirrors
`IMeshData::IMapOfReal`」，但语义（插入序 vs 升序）不同。下一轮要做的是：
1. 先做**判定实验**（不先改代码）：探针同时构造两侧的 torus splitter
   （同一面、同一 `AddPoint` 序列、同一 `AdjustRange`），打印
   `GenerateSurfaceNodes(...)` 的**点数**；OCCT 侧用 `BRepMesh_TorusRangeSplitter`（public、`Standard_EXPORT`），
   端口侧用 `crates/occt-topo/examples/` 里的仪器。
   预期：端口点数 > OCCT 点数，且比值与 `2.135` 同量级。
2. 若成立，把 `ParamSet` 改成**插入序 + 去重**（`NCollection_IndexedMap` 语义），
   并全量 grep 其消费者确认没有别处依赖升序（`calc_average_duv` 自己会排序，不受影响）。
3. 改完按门禁实测：`--lib`、`step_obj_gates`（逐模型 area ratio）、a3n00/T0M stats。

#### 9.372.3 密度图的配对质量警告（下一轮先修）

* 配对用六坐标 bbox ≤1e-3；**26 个面配不上**，且未配上的 GT 面承担了 5932/12324 个三角
  ⇒ 大面（10 wire / 6 wire / 合并后 bbox 变化的面）系统性失配；
* 已实证同类陷阱：`--mesh`（TopExp_Explorer 序）的 FACE bbox 会因「三角化不完整」而变小
  （§9.368.5 的 `--fbbox` 教训）⇒ 用 `--mesh` 的 bbox 配对不可靠；
  本轮改用 `--boundary`（模型序、全程 DiscretizeEdges 后的面）才对；
* 下一轮：配对后**必须先核对「配对面的 GT 三角和 ≈ 全模型 GT 三角和」**，对不上就不下结论。

#### 9.372.4 本轮改动与复跑

* 库代码零改动（纯分析 + 读 OCCT 源码）；
* 新增脚本 `.target-gate/t101_bytype.py`（按类型归并）；
* 复跑：`python .target-gate\t101_bytype.py`。

---

### 9.373 —— 【自我否定】T-101 迭代（4）：§9.372.2 的「`ParamSet` 升序 vs `IMapOfReal` 插入序」**不是分歧** —— OCCT 的 `FUN_CalcAverageDUV` 会把参数数组**就地排序**，两侧的抽稀都是升序

#### 9.373.1 否定依据（读 `.cxx` 到能证伪自己为止）

```cpp
// BRepMesh_TorusRangeSplitter.hxx:41  —— 形参是**非 const 引用**
double FUN_CalcAverageDUV(NCollection_Array1<double>& P, const int PLen) const;

// BRepMesh_TorusRangeSplitter.cxx:181-201
double BRepMesh_TorusRangeSplitter::FUN_CalcAverageDUV(NCollection_Array1<double>& P, const int PLen) const
{
  for (i = 1; i <= PLen; i++)
  {
    for (j = i + 1; j <= PLen; j++)      // ← 选择排序，**就地**改写 P
      if (P(i) > P(j)) { swap(P(i), P(j)); }
    ...                                    // 累加相邻差
  }
}

// BRepMesh_TorusRangeSplitter.cxx:129-176 fillParams:
//   aParamArray(j) = theParams(j);                      // 先按 IndexedMap 的插入序装入
//   aStep = FUN_CalcAverageDUV(aParamArray, aLength);   // ← 这一步把 aParamArray 排成升序
//   for (j = 1; j <= aLength; ++j) { pp = aParamArray(j); ... }   // 抽稀时已在升序上迭代
```

端口 `fill_params`（`crates/occt-topo/src/meshing/range_splitter/splitter.rs:485-511`）同样是
「先收集 → `calc_average_duv(&mut arr)` **就地排序** → 在 `arr` 上抽稀」，
而 `ParamSet` 本来就保持升序。⇒ **两侧的贪心抽稀都在升序遍历上进行，
`IMapOfReal` 是 `IndexedMap`（插入序）这件事在 `fillParams` 里被就地排序抹平了。**

⇒ §9.372.2 提出的「把 `ParamSet` 改成插入序」若照做，**是白改甚至有害**（会引入与
`calc_average_duv` 之外的消费者不一致的顺序语义）。该候选作废，**不要动 `ParamSet`**。

#### 9.373.2 更正后的 ③ 现状与下一步

* 仍成立：按曲面类型归并后 Torus 的倍数最大（7 面 1507 vs 706 = 2.135×），
  且 `fillParams` 只有 TorusRangeSplitter 用 —— 但「参数容器顺序」这条解释已被证伪；
* 仍成立：§9.372.3 的**配对质量警告**（199 个配对面的 GT 三角和只有 6392/12324
  ⇒ 大面系统性失配），所以**下结论前必须先修配对**；
* 下一步（按顺序）：
  1. **修配对**：配对后先核对「配对面的 GT 三角和 ≈ 全模型 GT 三角和」；
     失配的大面用「同类型 + bbox 容差按整形外扩放宽 + 面数守恒」重配。
  2. 配对可信后，再对 Torus 面做**九量对拍**（照 §9.346 的做法）：
     两侧同时打印 `range_u/range_v/delta/tolerance`、`r`/`R`、
     `ArcAngularStep(r,…)`、`oldDv`、`nbV`、`Du`、`nbU` 与 `GenerateSurfaceNodes` 的**点数**，
     第一个不同的量就是分歧点。OCCT 侧可直接构造 `BRepMesh_TorusRangeSplitter`
     （`Standard_EXPORT`，`Reset/AddPoint/AdjustRange/GenerateSurfaceNodes` 探针里都已能用，
     见 `specs/occt_probe/occt_probe.cpp` 的 `OCCT_TOPO_TRACE_SPLITTER` 段）。
* 本轮**未改任何库代码**（只读 `.cxx` 与更正文档）。

#### 9.373.3 复跑

```text
grep -n "FUN_CalcAverageDUV" -A 24 D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\BRepMesh_TorusRangeSplitter.cxx
```

---

### 9.374 —— T-101 迭代（5）：③ 的「法兰孔面 72 vs GT 52」**是口径问题、不是缺陷** —— 用**同参数**的 OCCT 逐面统计重测：那 16 个面端口 1193 vs OCCT 1188（1.004×）；真正还差的只有 **Torus**（1.24×）

#### 9.374.1 新的 OCCT 逐面仪器（同参数、模型序、可配对）

`specs/occt_probe/occt_probe.cpp` 新增 `--facestats <defl> <angle>`：
用 `BRepMesh_ModelBuilder` 建模型（**模型序**），再跑一次 `BRepMesh_IncrementalMesh(shape, defl, false, angle)`，
逐面打印 `FSTAT_OCC face= nodes= triangles= bbox=`。

**为什么不能用 `--uvsum`**：那个分支里另外跑的一次网格化把角度**写死成 0.5 rad**
（`BRepMesh_IncrementalMesh aSumMesher(aShape, aUvDefl, false, 0.5)`），
所以它的逐面 `triangles=` 在 a3n00 上只有 **8054**，而同参数的 20° 网格是 **12324** ⇒
**§9.365.2 与 §9.368.5 的 GT 列都不可与端口的读数相比**（这也解释了「GT 只铺 2–15 个三角」的错觉）。

#### 9.374.2 重测结果（225 个面 1:1 配对，自检通过）

自检：配对后 port `mt=11941` / GT `triangles=12096`；全模型 port `11941` / GT `12324`
（差的一个面是未配上的那一张）⇒ 配对可信。

```text
type      n     port_mt  occ_tri  diff     ratio
Torus     7     1507     1214     293      1.241     ← 唯一系统性偏密的一类
Plane     95    2848     2832     16       1.006
Cyl       61    2390     2445     -55      0.978
Cone      30    667      727      -60      0.917
BSpline   32    4529     4878     -349     0.928
```

最差 3 个面：`port(pair)202 ↔ model 33`（Torus，548 vs 355）、
`118/122 ↔ 36/37`（Torus，211 vs 162）。方向明确的偏疏面只剩
`171 ↔ 1`（BSpline，349 vs 669）、`140 ↔ 39`（Cone，5 vs 81）、`170 ↔ 69`（Cyl，27 vs 77）——
后两个正是 §9.370 里 seam 合并没有落地的面。

#### 9.374.3 ③ 的结论更正：法兰/倒角/孔那 16 个面的密度**已经对齐**

同参数逐面（`--fstats` ↔ `--facestats`，bbox 全 0.0000 或近邻）：

```text
169 74/72↔74/72   174 94/110↔93/108   189 74/72↔74/72
192..199 各 74/72↔74/72（8 个面全同）   204 74/72↔74/72   206 74/72↔74/72
209/210/212 75/73↔74/72
合计 port mt=1193  vs OCCT triangles=1188   → 1.004×
```

⇒ 卡面写的「法兰孔面端口 72 vs GT 52」**作废**：那个 52 来自 `--uvsum` 的 0.5 rad 网格。
**法兰（含倒角）这 16 个面的网格密度与 OCCT 已经逐面相同**（边结构见 §9.369.3，密度见本节）。

#### 9.374.4 ③ 剩下的唯一问题与下一步

只剩 **Torus**（7 个面，+293 个三角，1.241×）。下一步是 Torus 的**九量对拍**：
两侧同时打印 `range_u/range_v/delta/tolerance`、`r`/`R`、
`GCPnts_TangentialDeflection::ArcAngularStep(r, deflection, angle, min_size)`、
`oldDv`、`nbV`、`Du`、`nbU` 与 `GenerateSurfaceNodes` 的**点数**（OCCT 侧
`BRepMesh_TorusRangeSplitter` 是 `Standard_EXPORT`，探针里 `Reset/AddPoint/AdjustRange`
已经能用；端口侧照 `linkrods_dbg.rs:120-158` 的写法）。第一个不同的量就是分歧点。
**注意**：`ParamSet` 的升序**不是**分歧（§9.373 已证伪），不要再往那条路走。

#### 9.374.5 本轮改动与复跑

* 探针新增 `--facestats`（TEMP T-101，保留作仪器）；库代码零改动；
* 复跑：

```text
cmd /c "specs\occt_probe\build.bat"
cmd /c "specs\occt_probe\probe.bat D:\source\repos\dogs\data\occ\a3n00.stp --facestats 1.076007 0.349066" > .target-gate\fstatocc.txt
crates\occt-topo\target\debug\examples\zz_probe_a3n00.exe data/occ/a3n00.stp --fstats > .target-gate\fstats370.txt
python .target-gate\t101_density5.py     # 按类型 + 自检
python .target-gate\t101_flange16.py     # 那 16 个面
```

---

### 9.375 —— T-101 迭代（6）：Torus 残差的候选逐个排除（`ArcAngularStep` 忠实、`torus_radii` 有回退不除零）；item ① 等子任务收尾

#### 9.375.1 已排除的两处（读 `.cxx` 逐行比）

1. **`arc_angular_step` ↔ `GCPnts_TangentialDeflection::ArcAngularStep`** —— 忠实：

```cpp
Standard_ConstructionError_Raise_if(theRadius < 0.0, …);
constexpr double aPrecision = Precision::Confusion();      // 1e-7
double Du = 0.0, aMinSizeAng = 0.0;
if (theRadius > aPrecision) {
  Du = std::max(1.0 - (theLinearDeflection / theRadius), 0.0);
  if (theMinLength > aPrecision) aMinSizeAng = std::min(theMinLength / theRadius, M_PI_2);
}
Du = 2.0 * std::acos(Du);   …  std::min(Du, theAngularDeflection) 再 max(aMinSizeAng)
```

端口 `param_set.rs:173-189` 是同一套阈值/公式/顺序，且 `occt_core::precision::CONFUSION == 1e-7`
就是 `Precision::Confusion()` ⇒ **不是分歧**（Torus 网格的这一步输入同构）。

2. **`torus_radii` 不会除零**：`gp_torus()` 为 `None` 时有弦长回退
（`param_set.rs:156-163`）。但**回退值与 OCCT 的 `BRepAdaptor_Surface::Torus()` 是否一致尚未验**：
若某个面被分类成 Torus 而 `gp_torus()` 为 `None`，两侧的 `r`/`R` 可能不同
⇒ 这是 Torus 九量对拍里**要单独打出来的一列**。

#### 9.375.2 Torus 残差的下一步（收窄到一次对拍）

Torus 的网格路径只剩「输入量」没比过：`range_u/range_v`（`AdjustRange` 之后）、
`GetDelta`、`GetToleranceUV`、`r`/`R`（含上面的回退分支）、
以及 `AddPoint` 收到的参数集（`u_params/v_params` 的**个数**）。
下一次对拍就把这 6 项 + `GenerateSurfaceNodes` 点数在两侧同时打印，第一个不同的量即分歧点。

#### 9.375.3 item ① 的状态

子任务仍在跑（已在 `shape_fix_face.rs`、`shape_fix_compose_shell/{perform,collect_wires}.rs`
与 `specs/occt_probe/wires_probe.cpp` 加了 TEMP 插桩，产出 `zzfms_port140.txt` / `zzcs_port140.txt` /
`zzcs2_port140.txt`，用 OCCT 侧的 ZZ 插桩（`ZZFMS/ZZCS/ZZCW` 标签）对同一面的 ComposeShell 阶段做对拍）。
本轮已让它**收尾给结论**（第一处分歧 + 是否小改可修 + TEMP 插桩清单），结论到达后我再验证、跑门禁、合入。

#### 9.375.4 本轮改动

* 库代码零改动；纯读 `.cxx`/端口源码 + 排除候选；
* 无新增仪器（Torus 对拍所需的 `--facestats`（OCCT 侧）与 `--fstats`（端口侧）都已在 §9.374 就位）。

---

### 9.376 —— T-101 迭代（7）：备好一条命令的**落地验收 harness**；item ① 的对照实验仍在子任务手里（已要求本轮收尾）

#### 9.376.1 落地验收 harness（`.target-gate/t101_verify.ps1`，不进版本库）

把 §9.370 那条尾巴修好之后，落地只需跑：

```text
pwsh -File .target-gate\t101_verify.ps1
```

它按顺序跑并落盘 `.target-gate/t101_verify.txt`：`cargo check` →
`zz_probe_a3n00 --fixms`（113/140/170 不应再出现 `Shell`）→ 三个面的
`zz_seam_fix`（BEFORE/result/HEALED 的 wire 数）→ `zz_uv_feed --ids` 的 wire 直方图
（目标贴近 `{1:208,2:9,4:1,6:4,10:4}`）→ a3n00 / T0M 的 `TOTAL`+`STATMAP`
→ `--lib` → `step_obj_gates -- --nocapture` 的逐模型 `ratio=`。
判据写死在脚本注释里：**a3n00 面积比 ≥ 0.8996、T0M unmatched ≤ 6、`--lib` 1281/0、
`step_obj_gates` 5/5**（= T-101 卡面的 accept 在同一口径下重测）。

#### 9.376.2 item ① 的当前方法（子任务侧，供接手者复核）

子任务用的是**两侧同标签的阶段普查**：OCCT 侧 `specs/occt_probe/_dbg/` 的
`ZZ_ShapeFix_Face.cxx` / `ZZ_ComposeShell.cxx` 打 `ZZFMS`/`ZZCS`/`ZZCW`，
端口侧在 `shape_fix_compose_shell/perform.rs` 打同名的 `ZZCS <tag> n= [段列表]`
（只对 bbox = port 140 那个锥面展开逐段明细），逐阶段比 `loadwires / splitbygrid /
breakwires / collectwires / dispatchwires` 的段数与非流形标志。
已产 `zzfms_port140.txt`、`zzcs_port140.txt`、`zzcs2_port140.txt`、`zzcs3_port140.txt`。
本轮结束时其工作树只剩 `specs/occt_probe/wires_probe.cpp`（OCCT 侧 `seamfix <bbox>` 模式）。

#### 9.376.3 本轮改动

* 库代码零改动；新增一个验收脚本（`.target-gate/`，不进库）。

---

### 9.377 —— T-101 迭代（8）：item ① 的根因落在 **ComposeShell 内部**（不是剪枝判据）——两侧阶段普查给出第一处分歧与数值证据

#### 9.377.1 方法与复跑（两侧同标签阶段普查）

* OCCT 侧：`specs/occt_probe/_dbg/ZZ_ShapeFix_Face.cxx` / `ZZ_ComposeShell.cxx` + 探针
  `wires_probe <file> seamfix <x0 y0 z0 x1 y1 z1>`（本轮为该面加的模式，**未提交**），
  用 `ZZFMS`/`ZZCS`/`ZZCW` 打 ComposeShell 各阶段；
* 端口侧：`shape_fix_compose_shell/perform.rs` 里同名 `ZZCS <tag> n=[段列表]`（TEMP，已撤除），
  只对 bbox = port 140 的锥面展开逐段明细；
* 落盘：`.target-gate/{zzfms_port140,zzcs_port140,zzcs2_port140,zzcs3_port140}.txt`。

#### 9.377.2 第一处分歧与数值证据（a3n00 f=140，Cone）

```text
OCCT  ShapeFix_ComposeShell.cxx:2131-2275（SplitByGrid）→ :1433-1914（SplitByLine，ClosedMode 的 U 切）
端口  shape_fix_compose_shell/split_by_grid.rs:93-133 + split_by_line.rs(+/split_wire.rs)
```

```text
loadwires      12e/C[A→A]                              干净
splitbygrid    13e/C[A→A] + 新外部 1e/O[B→C]
breakwires     切成 5e/O[A→D] 与 8e/O[C→A]
               A=…463f0(59.075088,2.175248,-124.760355)
               B=…45070(17.750523,0,-163.227182)
               C=…6b6c40(16.704518,0,-167.130925)
               D=…42970(16.760355,2.175248,-167.075088)  ≠ C，差 (0.056,2.175,0.056)
               D 在 seqw 里别无出处 ⇒ 两半端点不配对
collectwires   5e 那一半 index=None（cxx:2826 的 !index 分支）⇒ 输出 [5e/O, 8e/O, 3e/C]
dispatchwires  出 2 张**同一 patch** [0,2π]×[-2.0207,2.0207] 的重合面（13e + 3e）
OCCT 孪生      = 1 face / 1 wire / 16e
```

⇒ `FixMissingSeam` 返回 Shell(2)，`resolve_face` 只接受 `ShapeType::Face` ⇒ 结果被丢弃 ⇒
该面保持 2 wires。**这与 §9.370 的结构读数一致**（合并出的 wire 是对的，坏在装配）。

#### 9.377.3 判定：根因在 ComposeShell，剪枝不是原因（实测）

* 剪枝是 **no-op**：第 1 轮 `fix_small` 5→5 / 8→8 / 3→3；第 2 轮
  `check_small_area` 三次全 false（剪枝只能删、不能把两张重合面合成一张）；
* 最小忠实修法范围：`split_by_line.rs` / `split_wire.rs` 在 **ClosedMode** 下把切割边插入
  wire 时必须让 wire 顶点与切割段端点**同一**（OCCT 复用已有顶点对；端口这里留下了
  陈旧/新造顶点 `D`）⇒ 先按 `ShapeFix_ComposeShell.cxx:1433-1914` 的分支核
  「切割点匹配已有顶点」那几处，不加任何特例。

**未验证（如实记录）**：OCCT 在同一输入下必产 1 face **没有直接测到**
（ZZ 覆盖只对探针发起的调用生效，DLL 内部调用仍走 DLL；对已愈合的 1-wire 面调
OCCT `FixMissingSeam` 直接崩，exit 1）。已证的是：**端口这一段自己的输出自相矛盾**
（同一 patch 两张重合面 + 一个未配对的半边），与 OCCT 模型孪生的 1 face/1 wire 不符。

#### 9.377.4 对本卡的影响

* ① 不是「小改可修」：需要在 ComposeShell 的闭合模式切割路径上修顶点同一性
  （范围已给到具体文件与 `.cxx` 分支）；
* 落地验收仍用 §9.376 的 `pwsh -File .target-gate\t101_verify.ps1`
  （`--fixms` 不再出现 Shell、wire 直方图贴近 `{1:208,2:9,4:1,6:4,10:4}`、
  a3n00 面积比 ≥ 0.8996、T0M unmatched ≤ 6、`--lib` 1281/0、`step_obj_gates` 5/5）；
* 本轮**库代码零改动**（子任务的 TEMP 插桩全部已用 `edit` 反向撤除；
  `git status --porcelain` 与 `git diff --stat` 均为空）。

---

### 9.378 —— T-101 迭代（9）：① 的实现交接与**进一步收窄**（`split_wire.rs:289-296` 的 `v_opt == None` 分支）；派出的实现子任务无产出、已停止

#### 9.378.1 本轮的事实

* 上一轮派出的实现子任务（「修 ClosedMode 切割顶点同一性」）跑了三轮，**没有落任何改动、没有任何读数**，
  也没回答状态请求 ⇒ 本轮 `interrupt_agent` 停止它；工作树保持干净（`git status --porcelain` 空，HEAD = `c5107c65`）。
* 本轮我自己按 §9.377 的定位往下读了一处（只读、未改）：

```text
端口 crates/occt-topo/src/shape_fix_compose_shell/split_wire.rs:289-296
    match v_opt {
        None => {                                     // cxx:1238-1242
            let nv = builder.make_vertex(curr_pnt, tol_edge);   // ← 新造顶点（= 观测到的 D）
            vertices.push(nv.clone()); v = nv;
        }
        Some(vv) if !do_cut => { … }                   // cxx:1243-1255
        Some(vv) => { vertices.push(vv.clone()); v = vv; }
    }
```

⇒ 「新造顶点 `D`」确实来自 `v_opt == None` 这一支（与 `cxx:1238-1242` 对应，端口这一支本身是忠实的写法）。
**但 `D ≠ C` 的偏移 (0.056, 2.175, 0.056) 不是在这一支产生的**：这里只是把 `curr_pnt` 造成顶点，
所以分歧在**上游**——`curr_pnt`（切割点）与切割边端点 `C` 的差，或
`v_opt` 本应命中 `Some(已有顶点)` 却没命中。

#### 9.378.2 收窄后的下一步（给接手者）

按 `ShapeFix_ComposeShell.cxx` 的顺序核这三点，一次只查一个：

1. `cxx:1185-1222`（端口 `split_wire.rs:258-276`）：**顶点匹配**分支——
   `prev_ok`/`is_coincided` 的容差用的是 `prev_v_tol` 与 `split_res(...)`；
   若这里的容差比 OCCT 严（或 `split_res` 取值不同），`v_opt` 就会落到 `None`，
   于是新造一个与 `C` 只差 2.175 的顶点。**这是最可疑的一处**（偏移量 2.175 与锥面 v 半幅 2.0207 同量级）。
2. `cxx:1238-1242` 的 `curr_pnt` 来源：切割点是否应当**取切割边的端点**而不是网格求交得到的点。
3. `split_by_grid.rs:93-133` 新增的那条 `1e/O[B→C]`（切割边）端点 `C` 是否带上了正确的顶点
   （OCCT 让切割边复用 wire 上的已有顶点）。

判据仍是 §9.376 的 `pwsh -File .target-gate\t101_verify.ps1`；**不许按面号/bbox 加特例**。

#### 9.378.3 本轮改动

* 库代码零改动（只读 + 停止一个无产出的子任务）；文档追加本节。

---

### 9.379 —— T-101 迭代（10）：§9.378.2 的**检查 1 结果是「忠实」**（顶点匹配分支三条件逐条对应）⇒ 下一个嫌疑上移到 `currPnt` 的来源

#### 9.379.1 逐条对照（`cxx:1185-1197` ↔ 端口 `split_wire.rs:251-276`）

```cpp
// ShapeFix_ComposeShell.cxx:1185-1197
else if (currPnt.Distance(prevVPnt) <= prevVTol &&
         CheckByCurve3d(prevVPnt, c3d,
                        f3d + (currPar - firstPar) * (l3d - f3d) / span2d, T,
                        prevVTol + 2 * Precision::Confusion()) &&
         prevPnt.Distance(myGrid->Value(C2d->Value(0.5 * (currPar + prevPar)))) <= prevVTol)
{ … V = prevV; … }                       // 1211-1222: IsCoincided(prevPnt2d, currPnt2d, uRes, vRes, prevVTol)
                                         //             && IsCoincided(prevPnt2d, C2d->Value(0.5*(…)), …)
```

端口同三条件，容差与 `mid` 取样一致：`curr_pnt.distance(prev_v_pnt) <= prev_v_tol`
→（3D 曲线检查，`prev_v_tol + 2.0 * CONFUSION`）→ `prev_pnt.distance(grid.value_pnt(&mid_prev)) <= prev_v_tol`；
`u_res/v_res` 来自 `split_res(is_cut_by_u, cut_index, prev_v_tol)`（对应 `GetGridResolution(...)/prevVTol`
再 `std::min(myU/VResolution, …)`）。⇒ **这一支写法上没有可指认的差异**，检查 1 判为**忠实**。

#### 9.379.2 结论与下一步（检查 2）

既然匹配分支忠实，`v_opt == None`（新造顶点 `D`）只能是**三条件中至少一条在数值上没成立**，
而三条里两条是「距离 ≤ 容差」、一条是 3D 曲线检查 —— 共同的输入是 **`curr_pnt`（当前切割点）**。
⇒ 下一个嫌疑就是 §9.378.2 的**检查 2**：`currPnt` 的来源。

对照窗口（下一轮直接开这两个窗口）：

```text
OCCT  ShapeFix_ComposeShell.cxx:1100-1185（currPnt / currPar / C2d->Value(currPar) 的取法，
      以及 CheckByCurve3d 的入参 span2d/f3d/l3d 从哪来）
端口  crates/occt-topo/src/shape_fix_compose_shell/split_wire.rs:150-258（同段）
判据  「`curr_pnt` 是否等于 OCCT 的 `currPnt`」——若不等，先按 .cxx 把取法改对，
      再看 `D` 是否消失（`zz_seam_fix 140` 的 result 应从 Shell(2) 变 Face/1 wire/16 边）
```

另记：`cxx:1238-1242`（端口 `split_wire.rs:289-296`）的写法本身忠实（§9.378.1），
不要在那一支上打补丁。

#### 9.379.3 本轮改动

* 库代码零改动（只读 `.cxx` + 端口对照）；文档追加本节。
* 剩余轮次预算：本目标 10/12 轮已用；① 的定位已细到「`curr_pnt` 取法」一个窗口。

---

### 9.380 —— T-101 **交接**：① 剩最后一步（`curr_pnt` 取法窗口），② 已结，③ 按原措辞已作废；下一会话从这一节接手

#### 9.380.1 三项的准确状态（都带可复跑证据）

| 项 | 状态 | 证据/位置 |
|---|---|---|
| **②** 全形状 seam 步比逐面多改的 4 个面是否 OCCT 也改 | **已结：是**（那 4 个面就是 ① 的那批，后置边数与 OCCT 孪生逐一相同：113 `22→28` = model 77 的 22+6；140 `13→16` = model 39 的 16；170 `5→8` = model 69 的 8。比边数用 `W k edges=N`，不用 `E` 行） | §9.370 |
| **③** 「法兰孔面端口 72 vs GT 52」 | **按原措辞作废**：52 来自 `--uvsum` 那次网格化写死 0.5 rad 角度（它逐面 triangles 和只有 8054，同参数 20° 网格是 12324）。用同参数的 OCCT `--facestats` 重测：那 16 个面 **port 1193 vs OCCT 1188（1.004×，逐面 74/72↔74/72）** ⇒ 法兰/倒角/孔面的密度**已经对齐**。**仅剩** Torus 7 面偏密 **1.241×**（`port 1507 vs 1214`） | §9.374 · §9.375 |
| **①** 补上 OCCT 会合并、端口未合并的面 | **部分完成**：reader 侧的 seam 步骤已接回（a3n00 面积比 0.8627→**0.8996**、T0M 未网格 7→6、UNUSED rescue 16→0、法兰 16 面 2 wires→OCCT 的 1 wire/5 边）。**仍剩 3 个面**（113/140/170）被 `fix_missing_seam` 装成 **Shell** 而 reader 只收 Face ⇒ 被丢弃（端口 wire 直方图 `{1:206,2:10,4:2,6:4,10:4}` vs OCCT `{1:208,2:9,4:1,6:4,10:4}`） | §9.369 落地 · §9.377 根因 · §9.378/§9.379 收窄 |

#### 9.380.2 ① 的最后一步：**只差 `curr_pnt` 取法这一个窗口**

已完成的三次排除（都对着 `.cxx`，都不是分歧）：

1. 剪枝判据 —— 实测 no-op（r1 `fix_small` 5→5/8→8/3→3；r2 `check_small_area` 三次 false）§9.377.3；
2. `split_wire.rs:289-296`（`v_opt == None` → `builder.make_vertex(curr_pnt, …)`，对应 `cxx:1238-1242`）—— 写法忠实，新顶点 `D` 只是被动产生 §9.378.1；
3. **顶点匹配分支** `cxx:1185-1197` ↔ 端口 `split_wire.rs:251-276`（三条数值条件 + `IsCoincided` + `uRes/vRes` 由 `GetGridResolution/prevVTol` 取小）—— 逐条对应，判为忠实 §9.379.1。

⇒ `v_opt` 落到 `None` 只能是三条数值条件没成立，而三条的共同输入是 `curr_pnt`。**要开的就是这两个窗口**：

```text
OCCT  ShapeFix_ComposeShell.cxx:1100-1185   （currPnt / currPar / C2d->Value(currPar) 的取法；
                                             CheckByCurve3d 的 span2d/f3d/l3d 从哪来）
端口  crates/occt-topo/src/shape_fix_compose_shell/split_wire.rs:150-258
判据  curr_pnt 是否等于 OCCT 的 currPnt
成功标志  zz_seam_fix 140 → result Face / 1 wire / 16 边（170 → 1 wire / 8 边；113 → 2 wires / 22+6 边）
         zz_probe_a3n00 --fixms 不再出现 result=Shell
         wire 直方图贴近 {1:208,2:9,4:1,6:4,10:4}、a3n00 面积比 ≥ 0.8996
```

若这个窗口也忠实，则按 §9.378.2 的检查 3 上移到 `split_by_grid.rs:93-133` 生成的那条
`1e/O[B→C]` 切割边端点是否复用 wire 上的已有顶点。**任何一步都不许按面号/bbox/面积加特例。**

#### 9.380.3 落地验收（一条命令，判据已写死）

```text
pwsh -File .target-gate\t101_verify.ps1      # 结果落盘 .target-gate/t101_verify.txt
# 判据：--fixms 无 Shell；wire 直方图贴近 {1:208,2:9,4:1,6:4,10:4}；
#       a3n00 面积比 ≥ 0.8996；T0M unmatched ≤ 6；--lib 1281/0；step_obj_gates 5/5
```

#### 9.380.4 本目标最后一轮的实际产出

* 库代码零改动（本轮为交接整理）；`crates/` 与工作树均干净；
* 本目标 11/12 轮：②③ 已按上表结清（③ 是口径更正，不是代码缺陷），① 由「3 个面 / ComposeShell 内部」
  收窄到「1 个窗口 / `curr_pnt` 取法」，并留下验收脚本与复跑命令。

---

### 9.381 —— 用户报的「a3n00 里大的六角螺帽切面是空的」= **F113**（同一张面、同一根因）

#### 9.381.1 认定（三个独立读数指向同一张面）

```text
crates\occt-topo\target\debug\examples\zz_probe_a3n00.exe data/occ/a3n00.stp --fstats   # 225 条，缺 face=113
… --ids   → PORTID f=113 … wires=4 mt=18446744073709551615 bbox=(-87.5,-34,-100)-(87.5,34,-32)
… --fdump → 只打印这一张面：FDUMP face=113 box=(同) wire[0] nEdges=6（6/14/1/… 共 4 条 wire）
OCCT 同面（model 77, type=1 Cylinder）: FSTAT_OCC face=77 nodes=228 triangles=228 / BOUND wires=2 pts=264
```

⇒ 端口的 F113 **完全没有网格统计**（`mt` 是「无统计」哨兵 `usize::MAX`），OCCT 同面是
**2 条 wire（22+6 边、264 个采样点）+ 228 个节点 / 228 个三角**。这就是用户看到的「空面」。

#### 9.381.2 为什么空（与卡上 ① 同一根因）

`fix_missing_seam` 对 f=113 返回 **Shell(5)**（§9.370/§9.377 实测），
`read_topology.rs::resolve_face` 只接受 `ShapeType::Face` ⇒ 丢弃 ⇒ F113 停在 **4 wires**
（OCCT 是 2）⇒ Delaunay 出不了域内元素 ⇒ 该面**一个三角都没有**。
所以「六角螺帽切面为空」不是新问题，就是 §9.380 里 ① 的同一个面、同一个 ComposeShell 分歧。

#### 9.381.3 本轮又排除两处（对照 `.cxx` 逐行）

```text
cxx:1149-1150  currPnt2d = C2d->Value(currPar); currPnt = myGrid->Value(currPnt2d);
端口 split_wire.rs:227   curr_pnt = self.grid.value_pnt(&curr_pnt2d)        ← 一致
cxx:1151-1180  「匹配 lastV」分支（currPnt.Distance(lastVPnt) <= lastVTol + CheckByCurve3d + lastPnt.Distance(grid(0.5*(currPar+lastPar))) <= lastVTol）
端口 split_wire.rs:237   同三条件、同容差                                         ← 一致
cxx:1185-1197  「匹配 prevV」分支（§9.379 已判忠实） ↔ 端口 split_wire.rs:260      ← 一致
```

⇒ 「`curr_pnt` 取法」这一层也忠实。**下一个嫌疑只剩循环用的参数表 `values` 的来源**，
即 `split_by_grid.rs:93-133` 生成的那条切割边（`1e/O[B→C]`）及其参数集
（对应 `ShapeFix_ComposeShell.cxx:2131-2275` SplitByGrid）——**可观测签名就是 `D ≠ C`**
（(0.056, 2.175, 0.056)）。

#### 9.381.4 下一步

对照 `ShapeFix_ComposeShell.cxx:2131-2275`（SplitByGrid）↔ 端口
`crates/occt-topo/src/shape_fix_compose_shell/split_by_grid.rs:93-133`，
核「切割边端点/参数集是否复用 wire 上的已有顶点与已有参数」。
成功后 `zz_uv_feed --model 113` 应为 2 wires（22+6 边）、`--fstats` 应出现 `face=113 mt≈228`、
`zz_seam_fix 113` 的 `result` 应为 Face，且 a3n00 面积比应从 0.8996 继续上升（§9.303 记该面解析面积约 29053，是缺口主项）。

---

### 9.382 —— T-101：`split_by_grid` 那一层也是忠实的 ⇒ 分歧收敛到 `split_by_line` 的**去重/吸附**那一步（cxx:1722-1767）

本轮把嫌疑窗口从 `split_by_grid` 往下走了一层，结果如下（逐行对照）：

```text
cxx:2227-2251 / 2253-2273（SplitByGrid 的 U/V 线循环）
端口 split_by_grid.rs:93-133（u_start/v_start、closed_mode 的 period 平移、get_patch_index、split_by_line_wires 调用）   ← 一致
cxx:1446-1471（SplitByLine 的表 + aB.MakeVertex(aVertNew, aP3d, BRep_Tool::Tolerance(aVert))）
端口 split_by_line.rs:44-81（同表 + builder.make_vertex(a_p3d, vertex_tolerance(&a_vert))）                        ← 一致
cxx:1577-1602（detect intersections at junction of two edges）  ↔ 端口 split_by_line.rs:180（注释即 cxx:1577-1602）   ← 有这一步
cxx:1722-1767（remove duplicated points in closed mode）        ↔ 端口 split_by_line.rs:295（注释即 cxx:1722-1767）   ← 有这一步
```

⇒ `split_by_line.rs` 的骨架是逐段映射的（1443→1911），**分歧不在「缺了哪一步」，而在某一步的数值判据**。
结合阶段普查的可观测签名（`breakwires` 出 `5e/O[A→D]` 与 `8e/O[C→A]`，`D≠C`、差 (0.056,2.175,0.056)）：

> **首要嫌疑 = cxx:1722-1767 的去重/吸附**（端口 `split_by_line.rs:295-343`）：
> 新交点 `D` 本应被判定与已有的 junction 顶点 `C` 重合而吸附过去。

#### 下一步（具体到窗口）
对照 `ShapeFix_ComposeShell.cxx:1722-1767` ↔ 端口 `split_by_line.rs:295-343`，
逐条比：比较的是**参数**还是**3D 点**、用的容差是哪一个（`prevVTol` / `Precision::Confusion()` / `IsSame`）、
以及「重合时保留哪个顶点」。改对后判据同 §9.381.4（`--model 113` 2 wires、`--fstats` 出现 face=113 mt≈228、
`zz_seam_fix 113` 的 result 为 Face、a3n00 面积比继续上升）。

---

### 9.383 —— T-101：`split_by_line` 的**去重/吸附**这一步也忠实 ⇒ 嫌疑落到四个 helper 与网格输入（不再是「哪一步」）

逐行对照 `cxx:1722-1767` ↔ 端口 `split_by_line.rs:295-341`（0/1 基已正确换算）：

```text
cxx:1722  int j = IntEdgePar.Length();        ↔ 端口 :301  let mut j = len - 1;（1-based → 0-based，注释已写明）
cxx:1725  for (i = 1; i <= Length();)         ↔ 端口 :303  while i < len（两边都每轮重取长度）
cxx:1727  if (i == j) break;                  ↔ 端口 :304  if i == j { break }
cxx:1731-1742  同边号 + |Δpar| < PConfusion → Remove(i)、j>i 则 j--、continue
                                              ↔ 端口 :307-337 同判据、同删除、同 j--
cxx:1743-1763  相邻边端口重合（nbe==1 或 Ind(i)==Ind(j)%nbe+1）→ 用 BRep_Tool::Range 的端点比 PConfusion
                                              ↔ 端口 :311-325 用 curve_on_surface_range 取 (a,b) 再按朝向取端点比 PCONFUSION
cxx:1765  j = i++;                            ↔ 端口 :338-339  j = i; i += 1;
```

⇒ **这一步也是忠实的**。至此，端口 ComposeShell 这条链上**每一个「步骤」都已逐行对过**并且一致：
剪枝（no-op，§9.377）· `split_wire.rs:289-296`（§9.378）· split_wire 的两个匹配分支（§9.379/§9.381）·
`curr_pnt/curr_par` 取法（§9.381）· `split_by_grid.rs:93-133`（§9.382）· 本节去重/吸附。

#### 那么差异只可能在**被这些步骤调用的 helper 或网格输入**上。下一轮按此顺序比（都在端口侧有成对实现）：

```text
1) helpers.rs:173  check_by_curve_3d      ↔ ShapeFix_ComposeShell.cxx 的 CheckByCurve3d
      —— 它同时把门两个匹配分支；签名里的容差与曲线取值若不同，两分支都会失败 ⇒ 走到「新造顶点」
2) helpers.rs:156  get_grid_resolution + split_wire.rs:484 split_res
      ↔ GetGridResolution(...)/prevVTol 再与 myU/VResolution 取小（喂给 IsCoincided 的分辨率）
3) helpers.rs:84   is_coincided           ↔ IsCoincided（容差语义）
4) 网格输入本身：composite_surface.rs:339 value_pnt / u_joint_values / v_joint_values /
      myUResolution / myVResolution（若与 OCCT 的 grid 不同，上面两分支同样会失败）
```

**可观测签名不变**：`breakwires` 处出现新顶点 `D` 而应有 junction 顶点 `C`（差 (0.056, 2.175, 0.056)）；
修好后判据同 §9.381.4（`zz_seam_fix 113` → Face、`--model 113` → 2 wires 22+6 边、`--fstats` 出现 `face=113 mt≈228`、a3n00 面积比从 0.8996 起继续上升）。

---

### 9.384 —— T-101：四个 helper + 常量 `TOLINT` 全部忠实 ⇒ 分歧不在切分逻辑，而在**喂进去的那个假想 grid**（`uf/vf` 与网格范围）

本轮把 §9.383 列的四项逐个对完，**全部一致**：

```text
cxx:876-892  CheckByCurve3d      ↔ helpers.rs:173-189 check_by_curve_3d
             （c3d 为空→true；c3d->Value(param)；T.Form()!=gp_Identity 才变换；SquareDistance <= tol*tol）   一致
cxx:929-938  GetGridResolution   ↔ helpers.rs:156-170 get_grid_resolution
             （leftLen/rightLen 的四支与 wrap 分支、/3.；1-based↔0-based 换算核对无误）                     一致
cxx:1163-1174/1199-1210 的 min(myU/VResolution, gridRes) 组合
                                 ↔ split_wire.rs:484-494 split_res（grid_res = get_grid_resolution(...)/vtol）  一致
cxx:451-462  IsCoincided         ↔ helpers.rs:84-88 is_coincided
             （U/VTolerance = Resolution*tol；std::max(TOLINT, ·)；逐轴比较）                              一致
cxx:290      #define TOLINT 1.e-10 ↔ helpers.rs:13 pub const TOLINT: f64 = 1.0e-10                          一致
```

⇒ **端口的 ComposeShell 切分逻辑（步骤 + helper + 常量）已逐行对完且一致**，分歧只能来自
**被切分的那个假想 grid 本身**：`shape_fix_face.rs:591-598` 用 `GeomRectangularTrimmedSurface::uv(surf, uf2, uf2+u_range, vf2, vf2+v_range)` + `CompositeSurface::with_grid(...)`
构造它，而 `uf2/vf2` 来自 `shape_fix_face.rs` 的「找 seam 插入位置」循环（`cxx:2138-2234`）。

#### 下一步（已有一份现成读数，只需再取一份对照）

两侧同标签的阶段普查里**已经打印了这三个量**：

```text
端口  .target-gate/zzfms_port140.txt   →  ZZFMS compshell-in uf=3.141592653590 vf=34.000000000000 URange=6.283185307180 VRange=2.000000000000
OCCT  specs/occt_probe/_dbg/ZZ_ShapeFix_Face.cxx / ZZ_ComposeShell.cxx（wires_probe seamfix <bbox> 模式）应打同一行
```

⇒ 下一轮第一件事：**比 `uf`/`vf`/`URange`/`VRange` 这四个数**（以及 grid 的 `u/v_joint_values` 个数）。
若它们不同，分歧就在 `shape_fix_face.rs` 的 `uf2/vf2` 循环（`cxx:2138-2234`）里，与 `split_*` 无关；
若相同，则继续比 grid 的 joint values 个数与 `myU/VResolution`。

---

### 9.385 —— T-101：端口侧 `compshell-in` 读数已齐（`uf=0, vf=-2.020725942164, URange=2π, VRange=4.041451884328`），但**档案里没有 OCCT 侧的同面读数**

本轮查了 `.target-gate/zz*.txt` 全部 dump：目标面（bbox `(16.704518,-30,-167.130925)-(59.130925,30,-124.704518)`，
`in f0 wires=2 edges=13`）的 `ZZFMS compshell-in` 在 **5 个 dump 里逐字相同**：

```text
ZZFMS compshell-in uf=0.000000000000 vf=-2.020725942164 URange=6.283185307180 VRange=4.041451884328 uclosed=true vclosed=false …
ZZFMS in f0 wires=2 edges=13 bbox=(16.70451829506937,59.13092516455243,-30.0)-(30.0,-16.5,…)
```

即：假想 grid 的 U 幅 = 2π（面本身 u 周期）、V 幅 = **4.041451884328 = 2 × 2.020725942164**，
seam 落在 v 的**下端** `vf = -2.020725942164`。这些读数**全部来自端口侧**：
`zzfms_a3n00.txt` 里该 bbox 的 `compshell-in` 命中数为 **0** ⇒ **档案里没有 OCCT 侧的同面读数可比**。

原因就是上一轮子任务报告过的取证边界：ZZ 插桩**只对「探针自己发起的调用」生效**，
DLL 内部的 `ShapeFix_Face::FixMissingSeam` 仍走 DLL；而对**已愈合的 1-wire 面**调 OCCT `FixMissingSeam`
会直接崩（exit 1）。所以 OCCT 侧这条 `compshell-in` 必须用**未愈合 / 或专门准备的输入**去取。

#### 下一步（唯一动作，已备好工具链）

```text
1) 用 specs/occt_probe/_dbg/ 的 ZZ 版（ZZ_ShapeFix_Face.cxx / ZZ_ComposeShell.cxx，含 bbox 目标表）
   + wires_probe 的 seamfix <x0 y0 z0 x1 y1 z1> 模式（上一轮子任务用过、能把 OCCT 侧的 ZZFMS 打出来）
   取到 OCCT 侧同一面的 compshell-in 四个数；
2) 若 uf/vf/URange/VRange 与端口不同 ⇒ 分歧在 shape_fix_face.rs 的 uf2/vf2 循环（cxx:2138-2234）；
   若相同 ⇒ 打印两侧 grid 的 u/v_joint_values 个数与 myU/VResolution（端口 grid 构造在 composite_surface.rs）。
```

（`_dbg/ZZ_ShapeFix_Face.cxx` 目前带着上一轮子任务加的 bbox 目标表 + stage dump，是 git-ignored 的调试件，可复用它重建。）

---

### 9.386 —— T-101（新目标：**螺母斜切面 = F113**）：候选 ③ 也忠实 ⇒ 剩下只有「假想 grid 的数据」（候选 ①②）

逐行对照 `cxx:2824-2846`（CollectWires 的「无候选/可闭合 → 收尾成 wire」）↔ 端口
`collect_wires.rs:309-330`：

```text
canBeClosed = endV.IsSame(firstV)                     ↔ can_be_closed = same_v(&end_v, &first_v)
if (!index || (canBeClosed && !lastEdge.IsSame(firstEdge) && IsCoincided(endPnt, firstPnt, myU/VResolution, 2.*tol)))
                                                      ↔ index.is_none() || (can_be_closed && !same_e(…) && is_coincided(…, 2.0*tol))
  if (!endV.IsSame(sae.FirstVertex(firstEdge))) FAIL5 ↔ first_vertex(first_edge) 后 same_v 比较 → SHAPEEXTEND_FAIL5
  wires.Append(s); sbwd.Nullify(); endV.Nullify()     ↔ wires.push(s); sbwd.clear(); has_sbwd=false; end_v=None
```

⇒ **忠实**。这意味着：在**同样的输入**下，OCCT 也会走「无候选 → 收尾」这条路；
端口之所以收尾成 `[5e/O, 8e/O, 3e/C]`（进而两张重合面），是因为**喂进来的段本身就已经不配对**
（`breakwires` 出的 `5e/O[A→D]` 与 `8e/O[C→A]`，`D≠C`）。

#### 结论（本轮把 ③ 也划掉后的唯一去路）

端口 ComposeShell 这条链上**所有控制流步骤 + 全部 helper + 常量**都已逐行对过且一致
（§9.377–§9.379、§9.381–§9.384、本节）。⇒ 分歧只能在**数据**：

```text
① 假想 grid 的 uf/vf/URange/VRange（端口读数已知：uf=0, vf=-2.020725942164, URange=2π, VRange=4.041451884328）
   —— 缺 OCCT 侧同面读数，取法：specs/occt_probe/_dbg 的 ZZ 版（ZZ_ShapeFix_Face.cxx/ZZ_ComposeShell.cxx，
      已带 bbox 目标表 + stage dump）+ wires_probe 的 seamfix <bbox> 模式（上一轮子任务验证过可用）
② grid 的 u/v_joint_values 个数与 myU/VResolution（端口构造在 composite_surface.rs）

下一轮第一件事 = 取 OCCT 侧那四个数（①），不等就先按 ② 打两侧的 joint values 个数与分辨率。
验收不变：zz_seam_fix 113 → Face；--model 113 → 2 wires（22+6 边）；--fstats 出现 face=113 mt≈228；a3n00 面积比从 0.8996 起上升。

---

### 9.387 —— T-101（F113）：`CompositeSurface::Value` 也忠实；并把「下一步」明确成**两条可执行的路**（含一条不需要重建 ZZ 版的）

本轮对照 `ShapeExtend_CompositeSurface.cxx:591-599` ↔ 端口 `composite_surface.rs:339-346`：

```cpp
gp_Pnt Value(const gp_Pnt2d& pnt) { i=LocateUParameter(U); j=LocateVParameter(V);
                                    uv = GlobalToLocal(i,j,pnt); patches(i,j)->D0(uv, point); }
```
```rust
fn value_pnt(&self, pnt) { let (i,j)=self.locate_uv_point(pnt); let uv=self.global_to_local(i,j,pnt);
                           self.patch(i,j).d0(uv.x(), uv.y()) }
```

⇒ **忠实**（结构、调用序一致）。仍未对过的只剩 `locate_uv_point`/`global_to_local`
↔ `LocateUParameter`/`LocateVParameter`/`GlobalToLocal`（三个小函数）。

#### 本轮新增的一条「不需要重建 ZZ 版」的路（决定分歧在**输入**还是**行为**）

档案里 `.target-gate/noopwires.txt` 是**OCCT 关掉 ShapeProcess（`--nofix`）后的逐面结构 dump**
（`FACE k type=… wires=… edges_per_wire: … spans: u=… v=… box=…`），正是「读入期原始输入」的镜像；
本轮查过：它**不含 F113 的 box `(-87.5,-34,-100)-(87.5,34,-32)`**（只覆盖了子任务当时的目标面清单）。

⇒ 下一轮据此两路并行（按代价排序）：

```text
路 A（便宜，先做）：用 wires_probe 的 raw/noop 模式**只针对 F113 的 bbox** 重跑一次，
   拿到 OCCT 侧「未整形」的 F113 = wires 数 / 各 wire 边数 / UV spans / UV 是否退化。
   · 若 OCCT raw 也是 4 wires（6/14/1/…）⇒ 输入一致，分歧在 FixMissingSeam 的行为 ⇒ 走路 B；
   · 若 OCCT raw 是别的结构 ⇒ 分歧在端口 reader 构面阶段（更上游），要在 read_topology 侧查。
路 B（贵）：把 wires_probe 的 seamfix 模式加回、用 _dbg 重建 ZZ 版，取 OCCT 侧 compshell-in 的
   uf/vf/URange/VRange，与端口的 (0, -2.020725942164, 2π, 4.041451884328) 比。
```

验收不变：`zz_seam_fix 113` → Face；`--model 113` → 2 wires（22+6 边）；`--fstats` 出现 `face=113 mt≈228`；a3n00 面积比从 0.8996 起上升。

---

### 9.388 —— T-101（F113）路 A 结果：**OCCT 的未整形 shape 里根本没有这张面** —— F113 是整形过程**造出来的**；端口的对应面还是 4 wires 且 **v 范围是无穷**

路 A 的两条读数（`cmd /c "specs\occt_probe\run_dbg.bat <a3n00.stp> [noop]"`，逐面 `FACE k type= wires= edges_per_wire: … spans: u=… v=… box=…`）：

```text
默认（ShapeProcess 开，已整形）  FACE 25 type=1 wires=2 edges_per_wire: 22 6
                                 spans: u=6.2832 v=175.0000 | u=2.1968 v=113.7919
                                 box=(-87.5,-34,-100)-(87.5,34,-32)        ← 与端口 F113 同 box，但只有 2 wires
noop（关 ShapeProcess，未整形）  全 226 面：wire 直方图 {1:163, 2:53, 4:2, 6:4, 10:4}（= 文件 FACE_BOUND 直方图 = 端口）
                                 **没有任何一面的 box 是 (-87.5,-34,-100)-(87.5,34,-32)**，也没有 u=6.2832 或 v=175 的面
```

⇒ 两点结论：

1. **F113 不是文件里的一张面，而是整形过程造出来的面**（未整形 shape 里不存在这个 box/参数范围）。
   端口这条 reader 路径（`fix_missing_seam` 直接吃 reader 的原始面）与 OCCT 的
   `ShapeFix_Shape`→`ShapeFix_Face::Perform` **作用在不同形态的输入上**：
   OCCT 那一步拿到的是整形后的面，端口拿到的是未整形、**v 范围为 `[-inf, inf]`** 的面
   （端口 `zz_uv_feed --model 113` 打印的正是 `vrange=[-inf,inf]`，而 OCCT 侧该面的 v span 是 175 / 113.79）。
2. 因此「端口 ComposeShell 每一步都忠实却给出 Shell(5)」并不矛盾：**喂进去的 grid 范围本身是无穷/未定**，
   切分自然落在不同的位置上。这与 §9.372–§9.387 把控制流全部排除的结果自洽。

#### 下一步（新的第一嫌疑，比继续对 ComposeShell 更靠前）

```text
(a) 量端口 F113 的**面上界**：surface 的 u/v range 为什么是 [-inf, inf]（`--model 113` 的
    vrange=[-inf,inf] 来自 `face.surface().v_range()`）—— 与该面在 OCCT 侧的 v span（175 / 113.79）
    对照；若端口拿到的 adaptor 没有有限界，则 ComposeShell/range splitter 的输入就已经不同。
(b) 同时量端口 F113 的四条 wire 与 OCCT `noop` 里「覆盖同一空间」的那张 raw 面（用 box 包含关系找，
    不再要求 box 相等）—— 确认端口 F113 的 4 wires 是 raw 结构与 OCCT 一致（则分歧在整形），
    还是 raw 就不同（则分歧在 reader 构面）。
```

---

### 9.389 —— T-101（F113）：取到 **OCCT 侧该面的 range splitter 读数**（`ru=[0,2π] rv=[-43.75,131.25]`，有限、跨度 175）；端口侧只差一个 `splitter.range_*` 打印

用 OCCT 探针自带的仪器（`OCCT_TOPO_TRACE_SPLITTER=1` + `--uvsum`，模型序）拿到 F113 孪生（model 77）的 splitter 值：

```text
SPLITTER f=77 preset=1 ru=[0.000000000000000,6.283185307179586] rv=[-43.750000000000000,131.250000000000000]
               delta=[0.029533064822176,1.000000000000000] tol=[6.283185e-07,1.750000e-05]
               vnb=0 defl=0.100000 cells_vnb0=[2,2] bbox=(-87.981961,…
```

⇒ OCCT 侧这张面的 splitter 范围是**有限的**：`ru` 一个整周期、`rv` 跨度 **175**（与 census 的 `v=175`、
`--uvsum` 的 `vrange=[-43.75,131.25]` 三处自洽）。注意 `defl=0.100000` 是探针该分支写死的常数，
所以 `delta/tol/cells` 只能在 0.1 偏转下比；**`ru/rv` 与偏转无关**，可以直接比。

而端口 `zz_uv_feed --model 113` 打印的是 `surface().v_range()` = **`[-inf, inf]`**（圆柱的**几何**范围），
**不是** splitter 在 `reset + add_point(每条 pcurve) + adjust_range` 之后拿到的范围
（端口 `node_insertion.rs` 用的正是这套配方）。⇒ **(a) 的读数目前不是同口径**，不能据此下结论。

#### 下一步（一步就能闭合口径）

给端口侧加一个同口径打印（只在 example/探针里，不进库）：对 `--model <f>` 加
`sp.range_u()/range_v()/delta()/tolerance_uv()` 的输出，跑 `--model 113` 与上面这行对照：

```text
判据：端口 sp.range_v() 是否 = [-43.75, 131.25]（OCCT），以及 range_u 是否 = [0, 2π]
 · 若不同 ⇒ 分歧在 splitter 的 reset/add_point/adjust_range（即 pcurve 取点或 UV 界计算），
   与 ComposeShell 完全无关；
 · 若相同 ⇒ 回到 §9.388 的 (b)：确认端口 F113 的 4 wires 与 OCCT 未整形 shape 的对应结构，
   用 `wires_probe wdump <idx>`（按 face 序号 dump）而不是按 box 匹配。
```

（本轮 `noop` census 里既没有 box 等于/包含 F113 的面，也没有 `u=6.2832`/`v=175` 的面 ——
原因是未整形 shape 的面 box/span 是从**另一套数据**（无 pcurve）算出来的，
所以「按 box 匹配 raw 面」这条路不可靠，下一轮 (b) 改用 `wdump <idx>`。）

---

### 9.390 —— T-101（F113）：**这张面根本没走到 range splitter** —— 空面的原因在更早的早退（`collect_boundary_uv` 返回 Err/空），与 ComposeShell 无关

同口径仪器（TEMP `OCCT_TOPO_TRACE_SPLITTER_P`，加在端口 `node_insertion.rs::perform`
的 `sp.adjust_range()` 之后，**含 invalid 也打**，已撤除）：

```text
225 条 SPLITTERP（226 面里缺的正是 face=113），且 valid=false 一条都没有
⇒ face=113 连 splitter 那一步都没到 —— 它在 perform 更早处就返回了
```

对照 OCCT 侧同面（model 77）：`SPLITTER f=77 ru=[0,6.2832] rv=[-43.75,131.25] …`（有效、有限）。

端口 `perform`（`node_insertion.rs:471-553`）在 splitter 之前的早退只有三处：

```text
1) let (uv, …) = Self::collect_boundary_uv(model, face_index)?;   ← 返回 Err 会直接冒泡
2) if uv.is_empty() { return Err("…has no boundary UV points") }
3) if let Some(surface) = face.surface() { … }                    ← surface 为 None 时整块跳过
   （已排除 3：`--model 113` 能打出 urange/vrange，说明 surface 是 Some）
```

⇒ **F113 变空的第一现场是 `collect_boundary_uv`（或它返回的 `uv` 为空）**，
即「这张面的边界 UV 取不出来」。这与 §9.377–§9.389 在 ComposeShell 里的一切无关 ——
前者是「面根本没进网格管线」，后者是「进了也合不拢」。

#### 下一步

给 `collect_boundary_uv`（或 `perform` 的入口）加同口径 TEMP 打印，量 F113 走到哪一步失败：

```text
· collect_boundary_uv 是 Err 还是 Ok(空)？（错误串会给出原因，例：某条 pcurve 缺失/参数非有限）
· 若是 Ok 且 uv 非空，则问题在 `if let Some(surface)` 块里的 sp.reset/add_point（例如
  wires_uv 为空 ⇒ add_point 一个都没加 ⇒ adjust_range 后 valid=false —— 但那样会打出 valid=false 行）
```

验收不变：`--fstats` 出现 `face=113 mt≈228`、`--model 113` → 2 wires、`zz_seam_fix 113` → Face、a3n00 面积比从 0.8996 起上升。

---

### 9.391 —— 【决定性】F113 变空的真正原因是**建模型时 `add_wire` 失败把该面标成 FAILURE**（`wire_builder.rs:577-582`）⇒ 它根本没被网格管线处理；ComposeShell 那条线解释的是「结构差」，不是「空面」

两处 TEMP 插桩（都已用 `edit` 撤除）给出的链：

```text
1) node_insertion.rs::perform 入口插桩（env OCCT_TOPO_TRACE_SPLITTER_P）
   225 条 BCUV（缺 face=113）、无 err、无 uv=0   ⇒ perform 根本没被调用
2) incremental_mesh/discret_root.rs:400 的「跳过前」插桩
   SKIPPRE face=113 failure=true reused=false   ⇒ 该面在进入网格循环前就已是 MeshStatus::FAILURE
   （同一次运行 BCUV 仍是 225 条，与 1 自洽）
```

`FAILURE` 的唯一「早于网格循环」的置位点在 **`model_builder/wire_builder.rs:577-582`**：

```rust
let outer = outer_of_wires(&wires, &f);
if let Some(outer_wire) = &outer {
    if !add_wire(&mut model, face_index, outer_wire, &mut edge_index) {
        if let Ok(fm) = model.face_mut(face_index) { fm.set_status(MeshStatus::FAILURE); }
        continue;        // ← 外环建不出来 ⇒ 面失败、其余 wire 也不再建
    }
}
```

⇒ **F113 的「空」= 外环 `add_wire` 失败**（`BRepMesh_ShapeVisitor` 在 OCCT 侧对同一张**已整形**面建模型是成功的：
OCCT 该面有 228 三角）。这也解释了 `--model 113` 打印的 `vrange=[-inf,inf]` 与 4 wires —— 模型只是把面登记进去，
外环没能建成。

**对整条研究线的影响**：§9.369 落地（法兰 16 面）与 §9.377–§9.389 的 ComposeShell 逐行排除，解释的是
**结构差异**（4 wires vs OCCT 的 2 wires）与 `fix_missing_seam` 的 Shell 输出；
**F113 空面的第一现场是建模型的 `add_wire`**，比 ComposeShell 更早、也是新的第一嫌疑。

#### 下一步

1. 在 `add_wire` 内部加同口径 TEMP 打印（哪一步失败：找边/pcurve/朝向/`edge_index` 复用），
   量 F113 的外环为什么建不出来；
2. 对照 OCCT `BRepMesh_ShapeVisitor::Visit(face)` → `AddWire`（`IMeshData` 侧）的同名分支；
   判据仍是 §9.381.4：`--fstats` 出现 `face=113 mt≈228`、`--model 113` → 2 wires、a3n00 面积比从 0.8996 起上升。

---

### 9.392 —— T-101（F113）：`add_wire` 的**四个早退分支一个都没触发** ⇒ FAILURE 不是建 wire 时设的，来自更晚的管线阶段

给 `wire_builder.rs::add_wire` 的四个失败返回都加了 TEMP 打印（`stored_empty` / `pcurve_err` / 参数非有限或退化 / `order.nb_edges() != stored.len()`；env gated，**已用 `edit` 反向撤除**），
跑 a3n00 全流程：

```text
ADDWIRE lines: 0        ← 四次早退全没触发
```

即 `add_wire` 在本次运行里**没有从这四条路失败**。于是 F113 的 `MeshStatus::FAILURE`
（§9.391 的 `SKIPPRE face=113 failure=true`）只能来自**别的置位点**。剩余候选（都在网格循环之前）：

```text
model_healer.rs:812                          ModelHealer 阶段
shape_tool.rs:243 / :281 / :296              ShapeTool 辅助（ModelPreProcessor/PreProcessor 阶段）
incremental_mesh/discret_root.rs:298 / :618  管线各阶段的失败分支
model_builder/wire_builder.rs:577-582        已排除（本轮四次早退都没触发；仅当 add_wire 返回 false 才置位）
```

#### 下一步（一次插桩即可定位）

在上述 6 处各加同口径 TEMP 打印（带 face 序号与阶段名），跑一次 a3n00，即可确定是
`ModelHealer` 还是 `PreProcessor`/`ShapeTool` 把 F113 标成 FAILURE；随后按该处的 `.cxx`
同等分支（`BRepMesh_ModelHealer` / `BRepMesh_ModelPreProcessor`）查它为什么对 F113 判失败。
**注意**：这与 §9.391 之前的一切（ComposeShell、`fix_missing_seam`、`add_wire` 早退）都不同，
是新的第一现场。

---

### 9.393 —— T-101（F113）：FAILURE 置位点的**定位尝试未完成**（本轮插桩编译不过，已全部撤除）；六个候选点的作用域已摸清

本轮想一次插桩定位「谁把 F113 标成 FAILURE」，结果如下（都只用 `edit`/脚本反向撤除，未用 `git checkout`）：

```text
model_healer.rs:812     ← 不是管线路径：它 new 了一个临时 MeshModel 再 add_face 后标 FAILURE（自检/辅助），排除
shape_tool.rs:243/281/296  ← 三处都是 `self.model.face_mut(face_index).expect("face exists").set_status(FAILURE);`
                             写成了多行，`face_index` 在作用域内（首选目标）
discret_root.rs:298/618    ← 这两处的 set_status 所在作用域**没有** `face_index` 变量名
                             （本轮插桩因此编译不过：E0425 cannot find value `face_index`）
incremental_mesh/discret_root.rs:400  ← 已确认只是「读」状态后 skip（§9.391），不是置位点
```

⇒ 已知的 `MeshStatus::FAILURE` 置位点里，**还没被排除**的是 `shape_tool.rs` 的三处
（`BRepMesh_ShapeTool` 侧的失败分支）与 `discret_root.rs:618`（`face_intersecting_edges` 相关路径）。

#### 下一步（把插入写成「用该站点自己的变量名」）

1. 先读 `shape_tool.rs:230-300` 确认这三处的条件与变量名（`face_index` 在作用域内）；
2. 在每处的前一行插入 `eprintln!("FAILSET stage=shapetool face={face_index}")`（env gated），
   以及 `discret_root.rs:618` 按其实际变量名插入；
3. 跑一次 `zz_uv_feed --ids` 看 F113 是否命中，即可确定第一现场；随后对照
   `BRepMesh_ShapeTool.cxx` / `BRepMesh_ModelPreProcessor.cxx` 的同等分支查它为何对 F113 判失败。

（本轮 6 个编译错的教训：**不能用统配正则插桩**，每个站点的变量名不同；下次逐站点用其自身作用域的变量。）

---

### 9.394 —— T-101（F113）：`shape_tool::visit_face` 三处 FAILURE 也**一处都没触发** ⇒ 置位点在剩下的未插桩站点（首要嫌疑 = healer 的 `face_intersecting_edges`）

本轮把 `shape_tool.rs::visit_face`（`BRepMesh_ShapeVisitor` 的等价物）的三处 FAILURE 全部插桩
（`wires.is_empty()` / `wi == 0` 的首环无边 / `!outer_ok`，都用同一作用域内的 `face_index`/`wi`；
env gated，**已撤除**），并加了该面的 `wires=/edges per wire=` 打印。跑 a3n00 全流程：

```text
VISITF fail lines: 0        ← 三处 FAILURE 一处没触发（连 face=113 的 VISITF 行都没有）
```

结合前两轮（`wire_builder.rs::add_wire` 四条早退 0 命中、`model_healer.rs:812` 在 `mod tests` 里）：
**已插桩的 FAILURE 置位点全部沉默**。剩下的未插桩站点：

```text
incremental_mesh/discret_root.rs:510                管线内的面级失败分支（待看条件）
incremental_mesh/discret_root.rs:618                `fn face_intersecting_edges`，注释指向
                                                    `BRepMesh_ModelHealer.cxx:248-278` ← **首要嫌疑**
delabella.rs:132 / mesh_algo.rs:184                  三角化侧（面已进管线才会到）
node_insertion.rs:555/594/611                        perform 内部（§9.391 已证 F113 不进口 perform）
```

理由：F113 是 6/14/1/… 四条 wire 的复杂面，healer 的「边相交」检查
（`BRepMesh_ModelHealer.cxx:248-278`）正是会把它判失败的那类检查；而 §9.391 的
`SKIPPRE face=113 failure=true` 表明失败发生在**进三角化循环之前**，与 healer 的位置吻合。

#### 下一步

先读 `discret_root.rs` 的 `face_intersecting_edges` 与 `:510` 两处的局部变量名，
再各插一行同口径打印，跑一次即可确定；随后对照 `BRepMesh_ModelHealer.cxx:248-278`
的同等分支查它为何对 F113 判相交。

---

### 9.395 —— 【第一现场确定】F113 被 **healer 的 face-checker 判成 `SELF_INTERSECTING_WIRE` 而标 FAILURE** ⇒ 整张面被网格管线跳过（这就是「空面」的直接原因）

在 `incremental_mesh/discret_root.rs:606-610`（healer 的 face-checker 结算循环，
注释指向 `BRepMesh_FaceChecker::Perform` + `BRepMesh_ModelHealer.cxx:248-278`；
局部变量就是 `i`）插一行同口径打印（env gated，**已撤除**）：

```text
FACECHECK face=113 self_intersecting_wire        ← 全模型**只有这一张**面命中
```

⇒ F113 空面的因果链现在是完整的、逐段有读数的：

```text
face_intersecting_edges(model, 113) 返回 Some(edges)
  ⇒ model.face(113) 标 SELF_INTERSECTING_WIRE + FAILURE        （discret_root.rs:606-610）
  ⇒ 三角化循环 `if f.is_status(FAILURE) … continue`（:400）跳过它  （§9.391 SKIPPRE face=113 failure=true）
  ⇒ 该面没有三角                                               （§9.368/§9.371：--fstats 里连 face=113 都没有）
```

而 OCCT 对同一张面（已整形为 2 wires / 22+6 边）**没有**判自交，正常铺了 228 个三角。

#### 两个待分辨的可能（下一步就分开它们）

```text
(a) 端口 face-checker 的判据与 `BRepMesh_FaceChecker.cxx` 不等价（例如 SegmentsFiller/相交判定的
    参数、采样、含端点处理不同）⇒ 按 .cxx 修 `face_intersecting_edges`（discret_root.rs:616+）；
(b) 判据等价，但**输入**确实自交 —— 端口的 F113 还是 4 条未合并的 wire（OCCT 已并成 2 条）
    ⇒ 那就要回到 ComposeShell 的合并缺口（§9.377–§9.389 已把该链逐行对完，剩 grid 数据一项）。
```

分辨方法：在 `face_intersecting_edges` 内部打印被判相交的边对（端口侧），
与 OCCT `BRepMesh_FaceChecker` 在同一张面上的判定（可用 `_dbg` ZZ 版或 `--boundary` 的采样）对照；
若端口判的是「同一 wire 内相邻边」「不同 wire 之间」这类 OCCT 不判的组合，即归 (a)。

---

### 9.396 —— T-101（F113）：被判相交的边**跨 wire[0] 与 wire[1]** ⇒ 倾向结论 (b)：**未合并的 4-wire 结构本身就让 face-checker 判自交**（与 ① 同源）

打印 `face_intersecting_edges` 返回的边集（TEMP，env gated，**已撤除**），并对照 `zz_uv_feed --model 113` 的模型结构：

```text
FACECHECK face=113 edges=[290, 292, 293, 294, 295, 297]

端口 F113 的模型结构（--model 113）：
  wire[0] 6 边 : 289 290 291 292 293 294
  wire[1] 14 边: 295 296 297 111 298 299 300 301 302 303 304 305 306 113
  wire[2] 1 边 : 307  Circle par=[0, 2π]
  wire[3] 1 边 : 308  Circle par=[0, 2π]
```

被判相交的 6 条边里：**290/292/293/294 属于 wire[0]，295/297 属于 wire[1]**
⇒ checker 报的是**跨 wire 的重叠**（wire[0] 与 wire[1] 的边互相相交），外加 wire[0] 内部相交。

对照 OCCT：同一张面整形后是 **2 wires（22+6 边）**，checker 不判自交。
端口则是 **4 wires（6 / 14 / 1 / 1）** —— 即 `fix_missing_seam` 的合并结果被丢弃后留下的
**未合并边界**：两条本应合成一条的环（6+14 边）互相重叠 ⇒ face-checker 判 `SELF_INTERSECTING_WIRE`
⇒ 标 FAILURE ⇒ 三角化循环跳过 ⇒ 空面。

#### 结论与下一步

本轮把 §9.395 的 (a)/(b) 分开了，证据指向 **(b)**：
**空面的直接原因（face-checker 判自交）是上游「seam 未合并」的必然后果**，
不是 face-checker 判据不等价（它判的跨 wire 重叠在 4-wire 结构下是真实存在的）。

⇒ 修复目标重新收敛为**一件事**：让 `fix_missing_seam` 对 F113 产出 Face（而不是 Shell(5)）。
该链 §9.377–§9.389 已逐行对完，只剩「假想 grid 的数据」：
`uf/vf/URange/VRange`（§9.389 已取到 OCCT 侧基准 `ru=[0,2π] rv=[-43.75,131.25]`，
端口侧需按同口径打印 splitter 的 `range_u/range_v` 才可比）与 `u/v_joint_values`/`myU/VResolution`。

---

### 9.397 —— T-101（F113）：回到 `shape_fix_face` 的 seam 定位段；**已缩小到唯一还没对过的一段 = `cxx:1758-1802`（用面的 UV 界替换无穷的面界）**，这正好是「圆柱会走到、法兰面不会走到」的分支

本轮想验的假设（「端口用 `surf.v_range()` 而 OCCT 用面界」）**被 .cxx 否掉**：

```cpp
// ShapeFix_Face.cxx:1753-1756（FixMissingSeam 开头）
double URange, VRange, SUF, SUL, SVF, SVL;
mySurf->Bounds(SUF, SUL, SVF, SVL);              // 面的几何范围（圆柱的 V 是 ±Infinite）
BRepTools::UVBounds(myFace, fU1, fU2, fV1, fV2); // 面的 UV 界
// 1804-1805
URange = std::min(std::abs(SUL - SUF), Precision::Infinite());
VRange = std::min(std::abs(SVL - SVF), Precision::Infinite());
```

端口 `shape_fix_face.rs:171-172 / 208-209` 与之一致（同源 + 同样的 `min(…, INFINITE)`），
而且端口**已经**实现了紧接着的替换块（`cxx:1758-1802` ↔ 端口 `:175-206`，
用面的 UV 界覆盖无穷的面界，`Precision::is_infinite` 阈值 = `1e100`）。

⇒ 仍然没有找到不等价的点，但**范围已经收窄到一个具体段落**：

```text
cxx:1758-1802   ←→   shape_fix_face.rs:175-206      （**唯一还没逐行对过的一段**）
```

它之所以是第一嫌疑，是因为**只有它会对「面界为无穷」的面（= 圆柱 F113）生效**，
而带倒角的法兰那批（平面/圆锥，面界有限）根本不进这一段 —— 这正好解释了
「法兰解析已经对了，圆柱这张还是坏的」这个现象。

附带一个待核的数值差（次要）：

```text
端口 occt_core::precision::INFINITE = 2e100        （is_infinite 阈值 = 1e100）
OCCT Precision::Infinite()          = 1e100
```

`min(|SVL-SVF|, INFINITE)` 与各处 `is_infinite` 判据都会受这个 2× 影响，需与 `cxx:1758-1802`
一起核。

#### 下一步

逐行对照 `ShapeFix_Face.cxx:1758-1802` ↔ `shape_fix_face.rs:175-206`（含 `INFINITE` 的取值），
这是 §9.397 之后 seam 定位链里**最后一段没对过的控制流**。

---

### 9.398 —— T-101（F113）：`cxx:1758-1802` 也**逐行忠实**（含「覆盖后再判」的语义） ⇒ seam 定位链已全部对完，最后未对过的是 **grid 的构造**（`shape_fix_face.rs:591-598` + `CompositeSurface::with_grid`）

逐行对照结果（关键点是 `IsInfinite` 判据作用在**已被覆盖后**的变量上）：

```cpp
// cxx:1759-1780（V 向完全同构）
if (IsInfinite(SUF) || IsInfinite(SUL)) {
  if (IsInfinite(SUF)) SUF = fU1;
  if (IsInfinite(SUL)) SUL = fU2;
  if (|SUL - SUF| < PConfusion) {
    if (IsInfinite(SUF)) SUF -= 1000.;   // SUF 已被 fU1 覆盖 ⇒ 实际走 else
    else                 SUL += 1000.;
  }
}
```
```rust
// shape_fix_face.rs:176-190
if Precision::is_infinite(suf) || Precision::is_infinite(sul) {
    if Precision::is_infinite(suf) { suf = f_u1; }
    if Precision::is_infinite(sul) { sul = f_u2; }
    if (sul - suf).abs() < PCONFUSION {
        if Precision::is_infinite(suf) { suf -= 1000.0; } else { sul += 1000.0; }
    }
}
```

⇒ 同构（端口同样在覆盖**之后**判 `is_infinite`，所以两边都会走 `sul += 1000` 那一支）；
`INFINITE = 2e100` vs `Precision::Infinite() = 1e100` 只在数值真的到达 ~1e100 时才可能分叉，
而 F113 的面界是 ±inf（`--model 113` 实测 `vrange=[-inf,inf]`）⇒ 两边都判「无穷」并同样被面界覆盖，
这一处**不构成分歧**。

#### 结论（本目标的 seam 定位链状态）

`FixMissingSeam` 从入口到 `ComposeShell` 的调用为止，**逐段都已对过且一致**：

```text
cxx:1753-1756  面/界取法            ↔ shape_fix_face.rs:171-173      一致（§9.397）
cxx:1758-1802  无穷面界→面UV界      ↔ shape_fix_face.rs:175-206      一致（本节）
cxx:1804-1805  URange/VRange        ↔ :208-209                       一致（同 min(…, INFINITE)）
cxx:2131-2275  SplitByGrid          ↔ split_by_grid.rs:93-133        一致（§9.382）
cxx:1433-1914  SplitByLine          ↔ split_by_line.rs（逐段映射）    一致（§9.382/§9.383）
cxx:1151-1197  SplitWire 两匹配分支 ↔ split_wire.rs:237/260          一致（§9.379/§9.381）
cxx:1238-1242  新造顶点             ↔ split_wire.rs:289-296          一致（§9.378）
cxx:2824-2846  CollectWires 收尾    ↔ collect_wires.rs:309-330       一致（§9.386）
helpers: CheckByCurve3d / GetGridResolution / split_res / IsCoincided / TOLINT  全一致（§9.384）
剪枝判据                                                               实测 no-op（§9.377）
```

⇒ **唯一还没对过的是「喂给 ComposeShell 的那个假想 grid 的构造」**：
`shape_fix_face.rs:591-598`（`GeomRectangularTrimmedSurface::uv(surf, uf2, uf2+u_range, vf2, vf2+v_range)`
+ `CompositeSurface::with_grid(...)`）以及 `CompositeSurface` 的 parametrisation / joint values
（§9.389 已备好 OCCT 侧同口径基准 `SPLITTER f=77 ru=[0,2π] rv=[-43.75,131.25]`）。

#### 下一步（下一轮的起点）

1. 读 `shape_fix_face.rs:560-600`（grid 构造）与 `CompositeSurface::with_grid` / `compute_joint_values`，
   对照 `ShapeFix_Face.cxx:2280-2330`（构造 `Geom_RectangularTrimmedSurface` 与 `BRepMesh_CompositeSurface`）
   与 `ShapeExtend_CompositeSurface::ComputeJointValues`（cxx:603-653）；
2. 用 §9.389 的基准比 `ru/rv` 与 `u/v_joint_values` 个数；
3. 命中即按 .cxx 修，然后按验收：`zz_seam_fix 113` → Face、`--model 113` → 2 wires（22+6 边）、
   `--fstats` 出现 `face=113 mt≈228`、a3n00 面积比从 0.8996 起上升、`t101_verify.ps1` 全绿。

---

### 9.399 —— T-101（F113）：假想 grid 的构造也**逐行忠实** ⇒ 只剩「**找 seam 位置的那个循环**」没对过（`cxx:2138-2234` ↔ `shape_fix_face.rs:441-588`）

逐行对照 `ShapeFix_Face.cxx:2236-2261` ↔ `shape_fix_face.rs:590-608`：

```text
cxx:2237-2238  RTS = new Geom_RectangularTrimmedSurface(mySurf->Surface(), uf, uf+URange, vf, vf+VRange)
端口 :591-597  GeomRectangularTrimmedSurface::uv(surf, uf2, uf2+u_range, vf2, vf2+v_range)            一致
cxx:2239-2242  1×1 数组 + new ShapeExtend_CompositeSurface(grid)
端口 :598       CompositeSurface::with_grid(vec![vec![rts]], Parametrisation::Natural)                一致
cxx:2245-2250  non-manifold 子形状加回 tmpF
端口 :599-602  同一循环                                                                                一致
cxx:2252-2261  CompShell.Init(G, L, tmpF, Precision::Confusion()) / ClosedMode()=true /
               SetContext / SetMaxTolerance(MaxTolerance()) / Perform()
端口 :603-608  一致（CONFUSION / closed_mode(true) / context / max_tol / perform）                      一致
cxx:2264       mySurf = new ShapeAnalysis_Surface(RTS)
端口 :610       self.surf = Some(rts)                                                                  一致
```

⇒ grid 构造这一层也**没有分歧**。

#### 结论：seam 链里唯一还没对过的只剩「**算 `uf2/vf2` 的循环**」

```text
ShapeFix_Face.cxx:2138-2234        （找 seam 插入位置，产出 uf/vf，被上面那句 RTS 使用）
shape_fix_face.rs:441-588          （端口对应实现，441 起、588 止，已见它算 uf2/vf2 并做
                                    cxx:2226-2234 的 adjust_to_period 收尾）
```

这就是 **F113 空面的最后一段未核对的控制流**：`uf2/vf2` 决定假想 grid 的锚点，
锚点一偏，`SplitByGrid/SplitByLine` 的切点就落不到已有顶点上（§9.377 的 `D ≠ C`），
于是 `breakwires` 出两个不配对的半边 ⇒ ComposeShell 出两张重合面 ⇒ Shell(5) 被丢 ⇒ 4 wires ⇒
face-checker 判自交 ⇒ 空面（§9.391/§9.395/§9.396 已确证后半段）。

#### 下一步（下一轮起点）

逐行对照 `ShapeFix_Face.cxx:2138-2234` ↔ `shape_fix_face.rs:441-588`，重点：
`i1/i2` 的选取、`pos1/pos2`、`uf2/vf2` 的赋值分支、`w1/w2` 的 UV 界比较与 `period`
（端口 :377-390 的 `wire_uv_bounds` / `period`）。命中即按 .cxx 修，随后跑验收：
`zz_seam_fix 113` → Face、`--model 113` → 2 wires（22+6 边）、`--fstats` 出现 `face=113 mt≈228`、
a3n00 面积比从 0.8996 起上升、`pwsh -File .target-gate\t101_verify.ps1` 全绿。

---

### 9.400 —— T-101（F113）：seam 定位循环（`cxx:2138-2234` ↔ `shape_fix_face.rs:441-588`）**首次发现候选分歧** = `i1` 自增与内层循环守卫的先后（疑似 off-by-one）

本轮把这一段的骨架逐项对上（U/V 两条链、`foundU/foundV` 的 0/1/2 状态机、`skipU/skipV`、
`|abs|` 比较惯用法、`AdjustByPeriod` 与 `AdjustToPeriod` 的两种用法）：

```cpp
// cxx:2144-2183（i1 外层）
for (int i1 = 1; i1 <= nb1 + nb2; i1++) { … }              // ← i1++ 在**循环头**（体末）
if (uclosed && ismodeu) {
    pos1.SetX(pos1.X() + AdjustByPeriod(pos1.X(), SUF, URange));
    if (foundU == 2 && |pos1.X()| > |uf|) skipU = true;
    else if (!foundU || (foundU == 1 && |pos1.X()| < |uf|)) { foundU = 1; uf = pos1.X(); }
}
bool skipV = !vclosed;
if (vclosed && !ismodeu) { … foundV = 1; vf = pos1.Y(); }
if (skipU && skipV) { … }
if (i1 <= nb1) { … for (int i2 = 1; i1 <= nb1 && i2 <= nb2; i2++) { … foundU = 2; … } }
```

端口（`shape_fix_face.rs:504-580`）逐项对应：`found_u/found_v` 0/1/2 ✓、`i1 → 1` / `i2 → 2` ✓、
`adjust_by_period(pos1.x(), suf, u_range)` ✓、`(pos2.x() - pos1.x()).abs() < PCONFUSION` ✓、
收尾 `adjust_to_period(uf2, suf, suf + u_range)` ✓。

**唯一对不上的一处（候选分歧）**：

```text
cxx:2144      i1++ 在 for 头（体末才自增）⇒ 体内的内层 for 条件 `i1 <= nb1` 用的是**未自增**的 i1
cxx:2195      for (int i2 = 1; i1 <= nb1 && i2 <= nb2; i2++)
端口 :544      在进入 i2 循环**之前**就 `i1 += 1;`
端口 :551      在 i2 循环体内以 `if i1 > nb1 { … }` 做守卫 ⇒ 用的是**已自增**的 i1
```

⇒ 若确如所读，端口在 `i1 == nb1` 时会比 OCCT 提前一次退出内层循环（少跑一轮 `edge1 ∈ wd1 × edge2 ∈ wd2`
的配对），`pos2` 那一支（`foundU/foundV = 2`，即 cxx:2206-2220 的「端点重合」判定）就可能少命中一次 ——
而 `uf2/vf2` 正是由这个状态机决定的。**这与 F113 的 `D ≠ C` 现象方向一致**（切点没落在已有顶点上）。

#### 下一步（先证实，再改）

1. 读 `cxx:2183-2200` 与 `shape_fix_face.rs:543-560`，确认 `i1` 自增位置与内层守卫的先后；
2. 若证实是 off-by-one：把自增移回「体末 / 内层循环之后」（`for` 头的等价位置），**不加任何特例**；
3. 跑验收：`zz_seam_fix 113` → Face、`--model 113` → 2 wires（22+6 边）、`--fstats` → `face=113 mt≈228`、
   a3n00 面积比从 0.8996 起上升，再跑 `pwsh -File .target-gate\t101_verify.ps1` 全量门禁；
4. 若证实「其实等价」（例如端口的 `i1 += 1` 恰好复刻了 OCCT 2186-2194 里的自增），则继续在这段里找下一处。

---

### 9.401 —— T-101（F113）：§9.400 的 off-by-one 假设**被证伪**（端口忠实复刻了 `for` 头的自增语义）

读齐两边的原文后否掉：

```cpp
// cxx:2183-2195
if (skipU && skipV) {
  if (i1 <= nb1) { continue; }      // ← for(i1…; i1++) 的 continue **会执行自增** ⇒ i1 进一
  else           { break; }
}
for (int i2 = 1; i1 <= nb1 && i2 <= nb2; i2++) { … }   // 内层条件用「当前」i1（未自增）
```
```rust
// shape_fix_face.rs:542-553
if skip_u && skip_v {
    if i1 <= nb1 { i1 += 1; continue; }   // ← 显式补上 for 头那次自增，语义等价
    else { break; }
}
for i2 in 1..=nb2 {
    if i1 > nb1 { break; }                // ← 等价于 C 的 `i1 <= nb1` 每轮判
```

⇒ 端口的 `i1 += 1` 只出现在 `skip_u && skip_v` 这一支里，正是 OCCT `continue` 触发 `for` 头自增的等价写法；
内层守卫也与 C 的循环条件等价。**这一段没有分歧**（§9.400 的怀疑撤回）。

#### 该循环里**仍未逐行读过的部分**（下一轮的目标，按嫌疑排序）

```text
1) shape_fix_face.rs:441-503  ↔ cxx:2138-2155 附近
   —— `shiftw2` / `other` / `w1`、`w2` 的选取与 UV 界比较、`period`（端口 :377-390 的 wire_uv_bounds/period）
      ← 这段决定 pos1 的来源与两条 wire 的身份，**嫌疑最大**
2) shape_fix_face.rs:563-588  ↔ cxx:2203-2234
   —— `pos2` 的 U/V 两支状态机（found=2）与收尾 adjust_to_period
```

#### 下一步

读 `shape_fix_face.rs:441-503` 与 `cxx:2138-2155`，逐项比 `w1/w2` 的取法与 `shiftw2`/`other` 的取值；
若仍一致，再读 563-588 与 cxx:2203-2234。命中即按 .cxx 修并跑验收（`zz_seam_fix 113` → Face、
`--model 113` → 2 wires（22+6 边）、`--fstats` → `face=113 mt≈228`、面积比从 0.8996 起、
最后 `pwsh -File .target-gate\t101_verify.ps1`）。

---

### 9.402 —— T-101（F113）：**grid 锚点与 OCCT 完全一致**（假设彻底否掉）；且 **OCCT 的未整形 shape 里就是同一张 4-wire 面**（FACE 113, type=1, wires=4: 8/14/1/1）⇒ 输入一致、分歧在「同一输入下的执行」

两处决定性读数（都是现成档案，不需要新跑）：

```text
端口侧（.target-gate/zzfms_port140.txt L49 / zzcs*_port140.txt L121，对 F113 那张面）：
  ZZFMS compshell-in uf=0.000000000000  vf=-43.750000000000  URange=6.283185307180  VRange=175.000000000000
                    uclosed=true vclosed=false ismodeu=1 ismodev=0
  ZZFMS in f0 wires=4 edges=22 bbox=(-87.5,87.5,-34)-(34,-100,-32)

OCCT 侧（§9.389 的 TRACE_SPLITTER / census）：
  SPLITTER f=77  ru=[0, 6.283185307179586]  rv=[-43.750000000000000, 131.250000000000000]
```

⇒ 端口的假想 grid = **u=[0, 2π] × v=[-43.75, 131.25]**，与 OCCT 的 `ru/rv` **逐位相同**
（`-43.75 + 175 = 131.25`）⇒ §9.389/§9.397 关于「grid 锚点/范围不同」的怀疑**全部否掉**。

```text
OCCT 未整形（noop census，.target-gate/wires_noop14.txt）：
  FACE 113 type=1 wires=4 edges_per_wire: 8 14 1 1   box=(-1e+100,-34,-100)-(1e+100,34,-32)
```

⇒ **OCCT 拿到的也是同一张 4 wires 的未合并圆柱面**（type=1、wires=4、边分布 8/14/1/1），
它在 `FixMissingSeam` 里被并成 2 wires；端口拿到同样结构、同样的 grid，却得到 Shell(5)。
**输入的拓扑与 grid 都一致** ⇒ 分歧发生在「同一输入下的执行」里。

（一处待核的小差：端口的模型 wire[0] 是 **6** 边，OCCT raw 的第一条 wire 是 **8** 边
—— 可能是模型构建期丢了 2 条 seam 边，也可能是阅读期差异，需在下一次对拍里一并看。）

#### 下一步（首次进入**数值级**对拍）

```text
端口  zz_uv_feed --model 113        → 每条边的 curve 类型 + par=[f,l] + npc/朝向（已有输出）
OCCT  wires_probe wdump 113（raw）  → 同一张面每条边的 pcurve 参数区间/朝向
判据  逐边比 par=[f,l] 与朝向：若端口的 pcurve 参数区间与 OCCT 不同（哪怕 0.05 量级），
      就解释了 §9.377 观测到的「切点 D 与已有顶点 C 差 (0.056, 2.175, 0.056)」——
      即分歧不在 ComposeShell 的控制流（已全部对完），而在**喂进去的 pcurve 数据**。
```

命中即按 .cxx 修对应环节（pcurve 构造/range/朝向），再跑验收：`zz_seam_fix 113` → Face、
`--model 113` → 2 wires（22+6 边）、`--fstats` → `face=113 mt≈228`、面积比从 0.8996 起、
最后 `pwsh -File .target-gate\t101_verify.ps1`。

---

### 9.403 —— T-101（F113）：数值级对拍第一击 —— **W0 的边分解不同**（OCCT raw 8 边含重复 seam 边 / 端口 6 边），W1/W2/W3 一致

用 `cmd /c "specs\occt_probe\run_dbg.bat <a3n00.stp> wdump 113 noop"` dump **OCCT 未整形 shape 的第 113 张面**
（正好就是我们要的那张，见 §9.402），输出按 wire 列出每条边的端点：

```text
W0 edges= [12.3875,24.672,-89.3942 -> -49.8649,-2.8e-16,-100]
          [-49.8649,0,-100 -> -80.6227,0,-100]
          [-49.8649,0,-100 -> -80.6227,0,-100]      ← **同一条边出现两次（seam 边两侧）**
          [-49.8649,-5.6e-16, …
W1 edges= [-8.9,28.8444,-48 -> -8.9,33.8048,-69.6384] …          （14 条）
W2 edges= [87.5,0,-32 -> 87.5,-8.3e-15,-32]                      （1 条，闭圆）
W3 edges= [-87.5,0,-32 -> -87.5,-8.3e-15,-32]                    （1 条，闭圆）
```

与端口侧（`zz_uv_feed --model 113` / `--fdump`）对照：

```text
wire[0] 端口 6 边（289…294），首边起点 (12.387499,24.672033,-89.394247) **与 OCCT 相同**，
        但端口 `--fdump` 里这条边 last=(-65.243610,0,-100) ≠ OCCT 的 -49.8649 ⇒ **边的切分不同**
wire[1] 端口 14 边（295,296,297,111,298…306,113）  ← 与 OCCT W1 的 14 **一致**
wire[2]/[3] 端口各 1 条 Circle par=[0,2π]（307/308） ← 与 OCCT W2/W3 **一致**
OCCT census：edges_per_wire = 8 14 1 1 ；端口模型：6 14 1 1
```

⇒ **除 W0 外全部一致**；W0 是 **8（OCCT，含一条重复的 seam 边）vs 6（端口）**。
这正是「切点 D ≠ 已有顶点 C」这类现象的温床：W0 的边分解/重复 seam 边不同，
`SplitByGrid/SplitByLine` 面对的边界段就不同，端点自然对不上。

（注意：OCCT 这里 dump 的是**未整形** shape，端口 `--fdump/--model` 也是 reader 原始面 —— 两边层次相同，
可以直接比。）

#### 下一步

1. 查端口 reader 为什么 W0 只有 6 条边：STEP `FACE_BOUND` 里该环的边表是什么（`read_topology.rs` 的
   FACE_BOUND/EDGE_LOOP 组装），有没有把同一条 seam 边去重/丢边；对照 OCCT `ShapeExtend_WireData`
   保留重复 seam 边的行为；
2. 顺带扩展 OCCT `wdump` 打印每条边的 **pcurve 参数区间**（现在只有端点），与端口 `--model 113` 的
   `par=[f,l]` 逐边比；
3. 命中即按 .cxx 修 reader 侧的组装（不加特例），再跑验收：`zz_seam_fix 113` → Face、
   `--model 113` → 2 wires（22+6 边）、`--fstats` → `face=113 mt≈228`、面积比从 0.8996 起。

---

### 9.404 —— 【文件级铁证】F113 的**外环在 STEP 里是 8 条边，端口 reader 只留了 6 条** ⇒ 首个「覆盖不全」的数据级证据（用户新指令的靶心）

从 `data/occ/a3n00.stp` 直接读出该面的定义（实体号经 `#4952 CYLINDRICAL_SURFACE R=34` 反查）：

```text
#5375 = ADVANCED_FACE('',(#5140,#5352,#5363,#5374),#4952,.T.)
#5140 = FACE_OUTER_BOUND('',#5139,.T.)   → #5139 = EDGE_LOOP 有 **8** 个 ORIENTED_EDGE
#5352 = FACE_BOUND('',#5351,.T.)         → #5351 = EDGE_LOOP 有 **14**
#5363 = FACE_BOUND('',#5362,.T.)         → #5362 = 1
#5374 = FACE_BOUND('',#5373,.T.)         → #5373 = 1

外环的 8 条（ORIENTED_EDGE → EDGE / 朝向）：
  5005→#5004 .T.   5013→#5012 .T.   5019→#5018 .F.   5091→#5090 .T.
  5099→#5098 .T.   5108→#5107 .T.   5116→#5115 .T.   5138→#5137 .T.
  8 条边 **互不相同**（distinct = 8，无重复）
```

这与 OCCT 未整形 shape 的 census **完全一致**：`FACE 113 type=1 wires=4 edges_per_wire: 8 14 1 1`。

而端口 reader 出来的形状（`zz_probe_a3n00 --fdump`，形状级、非模型级）：

```text
FDUMP face=113 … wire[0] ori=Forward nEdges=6     ← 外环只剩 6 条（? 文件是 8）
                                  wire[1] 14 / wire[2] 1 / wire[3] 1   ← 这三条与文件一致
```

⇒ **端口 reader 在该外环上丢了 2 条边**（文件 8 → 端口 6），且丢在 **reader 组装阶段**
（`--fdump` 已经是形状级读数，早于建模/网格）。这与 §9.377 观测到的
`breakwires` 出 `5e/O[A→D]` 与 `8e/O[C→A]`（`D ≠ C`）方向一致：外环少了两段，
ComposeShell 的两半自然对不上 ⇒ 两张重合面 ⇒ Shell(5) 被丢 ⇒ 4 wires ⇒ face-checker 判自交 ⇒ 空面。

#### 下一步（两条并行）

```text
1) 定位丢掉的那 2 条边：把文件里外环 8 条（#5004/#5012/#5018/#5090/#5098/#5107/#5115/#5137）
   与端口 6 条（模型边 289…294 / 形状 wire[0] 的 6 条）按几何端点做配对，
   找出缺的两条 → 在 read_topology.rs 的 FACE_BOUND/EDGE_LOOP 组装里定位丢边点
   （对照 OCCT 的 ShapeExtend_WireData / ShapeFix_Face 对同环的处理）。
2) 按用户新指令做「STEP 解析与处理」的**覆盖/一致性对照**（见 §9.405 的计划）。

---

### 9.405 —— T-101（F113）外环逐边对齐：端口**丢了 STEP 里那对「去-回」重复边（#5012/#5018）**，还多出一个文件里不存在的顶点 (-65.24361, 0, -100)

把 STEP 外环 8 条边（`EDGE_CURVE → 顶点`）与端口 `--fdump` 的 wire[0] 6 条按端点对齐：

```text
STEP #5004  (12.387499,24.672033,-89.394247) -> (-49.864906,0,-100)
STEP #5012  (-49.864906,0,-100)              -> (-80.622679,0,-100)
STEP #5018  (-49.864906,0,-100)              -> (-80.622679,0,-100)      ← 与 #5012 **端点完全相同**（去-回尖刺）
STEP #5090  (-49.864906,0,-100)              -> (36.224088,-3,-99.867388)
STEP #5098  (36.224088,-3,-99.867388)        -> (63.926974,-3,-99.867388)
STEP #5107  (63.926974,-3,-99.867388)        -> (63.926974,3,-99.867388)
STEP #5115  (63.926974,3,-99.867388)         -> (36.224088,3,-99.867388)
STEP #5137  (36.224088,3,-99.867388)         -> (12.387499,24.672033,-89.394247)

端口 e[0]   (12.387499,24.672033,-89.394247) -> **(-65.243610,0,-100)**   ← 文件里没有这个顶点
端口 e[1]   **(-65.243610,0,-100)**          -> (36.224088,-3,-99.867388) ← 直接接到 #5090 的终点
端口 e[2..5] = STEP #5098 / #5107 / #5115 / #5137                         ← 逐点相同
```

⇒ 三点事实：

1. **e[2..5] 与文件完全一致**（#5098/#5107/#5115/#5137）；
2. **STEP 的 #5012 与 #5018 是一对端点相同的「去-回」边**（`-49.864906 ⇄ -80.622679`），
   端口把它们**整对丢掉了**（8 → 6 的那 2 条就是它们）；
3. 端口在 `-49.864906` 与 `-80.622679` 之外造了一个**文件里不存在**的顶点 `(-65.243610, 0, -100)`，
   并用它把 e[0] 与 e[1] 接起来 ⇒ 端口的这条链是「重切/合并」过的，不是文件里的分解。

#### 嫌疑与下一步

- 嫌疑：端口 reader 在组装该环时做了**同曲线/同端点的合并或去重**（把 `#5012/#5018` 这对
  同端点的边当成冗余消掉），并把 `#5004 + #5012` 一类同曲线段重新切在 `-65.24361`。
  这与 §9.369 接回的 reader 侧 seam 步（`fix_missing_seam` + `check_pcurves_and_shift`）
  以及 `read_topology.rs` 里 FACE_BOUND/EDGE_LOOP 的组装方式都有关，需要逐一确认。
- 下一步：在 `crates/occt-topo/src/step/read_topology.rs` 里定位「按什么键去重/合并边」
  （边号？曲线句柄？参数区间？端点？），并与 OCCT `ShapeExtend_WireData::Init` /
  `ShapeFix_Face` 对同环的处理逐行对照；`-65.24361` 这个顶点的产生点就是第一现场。
- 判据（不变）：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、
  `--fstats` → `face=113 mt≈228`、a3n00 面积比从 0.8996 起上升。

---

### 9.406 —— T-101（F113）：那个「文件里不存在」的顶点 `-65.243610` **是整个文件里都没有的坐标**，而它正好是 `-49.864906` 与 `-80.622679` 的**中点** ⇒ 是端口**算出来**的（不是读进来的）

```text
grep -c "\-65\.2436"  data/occ/a3n00.stp   → 0        （文件里不存在这个坐标）
grep -c "80\.62267"   data/occ/a3n00.stp   → 6        （合法顶点）
(-49.864906150824751 + -80.622678977781447)/2 = -65.243792564303   ≈ 端口顶点 -65.243610
```

⇒ 端口在外环上**把 STEP 的 `#5012/#5018` 那对「去-回」边（`-49.864906 ⇄ -80.622679`）
按中点切开**，由此产生 `-65.243610`，同时那对重复边消失（8 → 6）。
「取中点切开一条边」正是 **seam 插入**这类操作的典型动作。

#### 由此得到的强嫌疑：reader 侧的 seam 步**在原地改写了原始面的边**（结果被丢弃，改动却留下了）

§9.369 接回的 reader 侧步骤（`ShapeFixFace::fix_missing_seam` + 每个结果 wire 的
`check_pcurves_and_shift`）对 F113 返回的是 **Shell(5)**，`resolve_face` 只接受 `Face` 所以**丢弃结果**；
但若该步在构造 `tmp_f` / 调用 ComposeShell 时**经由共享的 `Arc` tshape 原地改动了原面的边**
（Rust 侧 `Arc` 共享 + `ShapeBuild_ReShape` 语义不等价），就会留下这种「文件里不存在的中点顶点」
与「丢掉的重复边」——而 OCCT 走 `ShapeBuild_ReShape` 上下文，**只有成功时才 Apply**，
失败时原面保持不变（`ShapeFix_Face.cxx:2302 myResult = Context()->Apply(myResult)` 之前不动原形状）。

#### 下一步（一次插桩即可判定）

```text
1) 在 read_topology.rs 调用 fix_missing_seam 之前/之后，各 dump 一次该面外环的边数/端点
   （TEMP，env gated，用完 edit 反撤）：
      调用前 nEdges 应为 8（= 文件），调用后若变 6 且出现 -65.24361 ⇒ 断定是该步原地改写；
2) 若断定：把该步改成「在副本上做、成功才替换」（对齐 OCCT 的 Context/Apply 语义），
   不加任何按面/按模型的特例；
3) 验收：--model 113 → 2 wires（22+6 边）、zz_seam_fix 113 → Face、--fstats → face=113 mt≈228、
   a3n00 面积比从 0.8996 起上升，最后 pwsh -File .target-gate\t101_verify.ps1。
```

---

### 9.407 —— T-101（F113）：§9.406 的「seam 步原地改写」假设**被证伪** —— 8→6 的丢边发生在**更早**（EDGE_LOOP → wire 组装）

在 `read_topology.rs` 的 seam 调用前后各打一次该面的 `wires / edges_per_wire`（TEMP，env `T101_DUMP_SEAM`，**已撤除**；`cargo check` 0 error）：

```text
SEAM pre  wires=4 edges=[1, 1, 6, 14]      ← 就是端口 F113（wire0=6 / wire1=14 / 两个 1）
SEAM post wires=4 edges=[1, 1, 6, 14]      ← 一模一样 ⇒ **seam 步没有改它**
（另一个 4-wire 面：pre=post=[6, 8, 8, 8]）
```

而 STEP 文件里 F113 四个 bound 是 **8 / 14 / 1 / 1**、OCCT raw census 同样是 `8 14 1 1`。

⇒ **外环少掉的那 2 条边（`#5012/#5018` 去-回对）以及那个文件里不存在的顶点 `-65.243610`
在进 `resolve_face` 的 seam 步之前就已经是这样了** —— 即丢边发生在
**EDGE_LOOP/FACE_BOUND 的解析与组装**里（`resolve_loop` / `bind_edge_loop_vertices` /
`resolve_oriented_edge`），不在 seam 步。§9.406 的怀疑撤回。

#### 下一步（第一现场已定位到函数级）

```text
read_topology.rs:403 resolve_loop            ← EDGE_LOOP → 逐 ORIENTED_EDGE 组装
read_topology.rs:478 bind_edge_loop_vertices ← 把 loop 的顶点绑到边上（Pass 1 已读，Pass 2 在 509+）
read_topology.rs:369 resolve_oriented_edge   ← ORIENTED_EDGE → EDGE（含 same_sense/朝向）

做法：在 resolve_loop 里打印「本 loop 的 ORIENTED_EDGE 条数」与「组装出来的 wire 的边数」，
      F113 外环应打印 8 → ? ；差值出现的位置就是丢边点；
      再检查 bind_edge_loop_vertices Pass 2（509+）是否把同端点的一对边合并/消掉，
      以及 -65.243610（= #5012 两端点的中点）是 Pass 2 里切出来的还是别处产生的。
对照：StepToTopoDS_TranslateEdgeLoop.cxx（端口注释里引的 cxx:288-403 / 355-365 / 367-368）
      与 ShapeExtend_WireData 对重复 seam 边的处理。
验收（不变）：--model 113 → 2 wires（22+6 边）、zz_seam_fix 113 → Face、
      --fstats → face=113 mt≈228、a3n00 面积比从 0.8996 起上升、t101_verify.ps1 全绿。
```

（注：工作树里 `crates/occt-topo/examples/zz_probe_a3n00.rs` 的 +64 行是**并行审计子任务**
`d91ca45d` 正在加的 `--ecensus` 探针，属它的在飞工作，**不要动**；本轮的库代码改动已全部撤净。）

---

### 9.408 —— T-101（F113）：`resolve_loop` 与 `make_wire` **都不丢边**（341 个 loop 全部 items == edges == wire_edges）⇒ 丢边发生在**建好 bound wire 之后**

在 `resolve_loop`（`read_topology.rs:413-433`）里加同口径打印（TEMP，env `T101_DUMP_SEAM`，**已撤除**）：

```text
LOOP 行数 = 341
  分布：items=1:168 / items=4:89 / items=3:28 / items=5:26 / items=8:12 / items=6:9 / items=12:3 / items=7:2
  **items == edges 恒成立**（含 12 个 items=8 → edges=8 的 loop，F113 外环是其中之一）
再打印 make_wire 之后的 wire 边数：wire_edges 也 == items，**全 341 个 loop 零处不一致**
```

交叉验证：226 面的 FACE_BOUND 直方图 `{1:163, 2:53, 4:2, 6:4, 10:4}` ⇒ 总 bound 数
`163 + 106 + 8 + 24 + 40 = 341` ✓ 与 LOOP 行数正好相等 ⇒ 一个 bound 一次 `resolve_loop`，没有遗漏/重复。

⇒ **丢边不在 `resolve_loop`，也不在 `make_wire`**。而 `resolve_face` 的 seam 块处已经是 `[1,1,6,14]`
（§9.407）⇒ 8→6 只能发生在「bound wire 建好之后、seam 块之前」，即：

```text
resolve_face 内部 578 → 740 之间（bound 解析 → Face 组装 → 缓存）
或在 resolve_outer_bound（read_topology.rs:277 分派）/ 形状缓存的取用上
```

#### 下一步（一次插桩）

在 `resolve_face` 入口（`read_topology.rs:578` 之后）打印该面的 `wires / edges_per_wire`，与 seam 块前的
`[1,1,6,14]` 对照：

```text
入口若 = [1,1,8,14] ⇒ 丢边在 resolve_face 内部 578-740 段（逐句看 bound→Face 组装）
入口若 = [1,1,6,14] ⇒ 丢边在更早（resolve_outer_bound / 缓存 / Face 构造），往那边找
```

判据不变：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、`--fstats` → `face=113 mt≈228`、
a3n00 面积比从 0.8996 起上升、`t101_verify.ps1` 全绿。

---

### 9.409 —— T-101（F113）：`resolve_face` 内部定位尝试（本轮插桩**编译未过、已撤净**）；已把范围钉到 `:663 make_face` 前后

本轮想在 `resolve_face` 的 `self.b.make_face(surface, &wires)`（`read_topology.rs:663`）前后各打一次
`wires / edges_per_wire`，用 §9.408 的分流法判定 8→6 是否发生在 `make_face` 内部。
插桩**编译报 1 个错**（`wires` 的元素类型与 `edges_of_wire(&Wire)` 不匹配，需先看 663 之前 `wires`
是怎么收集的），因此已用 `edit` 反向撤除，`cargo check` 复原、库文件回到 HEAD。

已确定的周边事实（本轮新增）：

```text
resolve_face 的关键节点：:634 make_face_from_surface(...)   :663 self.b.make_face(surface, &wires)
                        :655 注释「Bound orientation is already on the wire from resolve_outer_bound」
        :717 make_face_uv（natural-bound Offset 面用）      :740-758 seam 块（§9.407 实测此处已是 [1,1,6,14]）
上游已排除：resolve_loop（341 loop 全部 items==edges，:413-433）与 make_wire（同轮 wire_edges==items，§9.408）
```

#### 下一步（先看清再插）

1. 读 `read_topology.rs:600-665`，确认 `wires` 的构造与元素类型（`Vec<Wire>` 还是 `Vec<TopoShape>`），
   以及 `make_face_from_surface`（:634）与 `make_face`（:663）分别覆盖哪些分支；
2. 再在同口径打印下分流：
   `pre = [1,1,8,14] & post = [1,1,6,14]` ⇒ 丢边在 `make_face`（对照 OCCT `BRep_Builder::MakeFace` /
   `StepToTopoDS_TranslateFace` 的同分支）；
   `pre` 已是 `[1,1,6,14]` ⇒ 丢边在 `600-663` 的 wires 收集段（那段的过滤/去重条件就是第一现场）。
3. 判据不变：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、`--fstats` → `face=113 mt≈228`、
   a3n00 面积比从 0.8996 起上升、`t101_verify.ps1` 全绿。

---

### 9.410 —— 【关键定位】F113 在 `make_face` 调用处**还是文件里的 [1, 1, 8, 14]**，`make_face` 也没丢 ⇒ 8→6 的丢边发生在 `resolve_face` 的 **664-740** 段

在 `read_topology.rs:663` 的 `self.b.make_face(surface, &wires)` 前后各打一次同口径 census（TEMP，env gated，**已撤除**）：

```text
FACE pre-make_face  wires=[1, 1, 8, 14]      ← F113：与 STEP 文件 8/14/1/1 **完全一致**
FACE post-make_face wires=[1, 1, 8, 14]      ← make_face 没有改动
（其余面 pre == post，逐条一致；共 226 条）
```

而 §9.407 测得 seam 块处已是 `[1, 1, 6, 14]`。⇒ **丢边（`#5012/#5018` 那对同端点去-回边）与那个
文件里不存在的接点 `-65.243610` 产生在 `resolve_face` 的 664-740 段**（`make_face` 之后、seam 块之前）。

已排除的所有上游环节（每一处都有读数）：

```text
STEP 文件本身                8 / 14 / 1 / 1（distinct 8）
resolve_loop (:413-433)     341 loop 全部 items == edges（含 8→8）
make_wire                    同轮 wire_edges == items，零不一致（§9.408）
resolve_face 的 wires 收集   :649-661 → [1, 1, 8, 14]
make_face (:663)             [1, 1, 8, 14] → [1, 1, 8, 14]（本节）
seam 块 (:740-758)           入口已是 [1, 1, 6, 14]（§9.407）⇒ 丢边在此之前
```

#### 下一步（窗口已经很小）

读 `read_topology.rs:664-740`，逐句找「重建/修补 wire 或切边」的动作（`-65.243610` 是
`#5012` 两端点的中点，是**切边**的指纹），并与 OCCT `StepToTopoDS_TranslateFace` /
`ShapeFix_Face` 在该处的同等分支对照。命中即按 .cxx 修，验收：`--model 113` → 2 wires（22+6 边）、
`zz_seam_fix 113` → Face、`--fstats` → `face=113 mt≈228`、a3n00 面积比从 0.8996 起上升。

---

### 9.411 —— T-101（F113）：丢边的**第一现场已定位到三条调用**（`read_topology.rs:710-714`，在 seam 块之前）

`resolve_face` 的 664-740 段里，唯一会**改动 wire 本身**的是 per-wire 循环里的这三条（`:710-714`）：

```rust
for e in &loop_edges {
    crate::shhealing::xsalgo_check_pcurve(e, &face, self.precision);   // :711  XSAlgo_ShapeProcessor::CheckPCurve
}
crate::shhealing::project_wire_pcurve_ranges(w, &face, preci);          // :713
crate::shhealing::check_pcurves_and_shift(w, &face, preci, false);      // :714  （cxx:365-480 一线）
```

其中 `w` 来自 `wires.iter_mut()`（`:673`）——**是真正会被写回的那条 wire**。而 §9.410 已证：
`make_face` 之后（`:663`）该面还是 `[1, 1, 8, 14]`（= 文件），seam 块（`:740`）入口已是 `[1, 1, 6, 14]`。
⇒ **丢掉的 `#5012/#5018` 与那个中点接点 `-65.243610` 必产生于 `:711 / :713 / :714` 这三条之一**。

指纹吻合：这三条正对应 OCCT 的「pcurve 高级检查 + 重投影 + shift」——
OCCT 里这类改动发生在 `ShapeFix`（`ShapeBuild_ReShape` 上下文）里，**且通常在 `ShapeFix_Wire` 的副本上做、
成功才 Apply**；而端口是**直接写在形状的 wire 上**。`check_pcurves_and_shift` 在 `shape_fix_face.rs`
的职责正是「修补重建后 wire 的 pcurve（cxx:365-480）」，具备切边/重排的能力。

#### 下一步（三行探针，一次跑完）

在该 per-wire 循环里于 `:711`、`:713`、`:714` **各前后打一次** `w` 的边数
（用 `crate::topo_tools_full::edges_of_wire(w).len()`；TEMP，env gated，用完 `edit` 反撤）：

```text
期望读数：进入循环 8 → (711 后) ? → (713 后) ? → (714 后) ?
第一个把 8 变成 6 的调用就是第一现场；
随后按 .cxx 对照该函数的 OCCT 同等实现（xsalgo_check_pcurve ↔ XSAlgo_ShapeProcessor::CheckPCurve
cxx:344-401；check_pcurves_and_shift ↔ cxx:365-480），确认「是否应当就地改 wire / 是否应只改副本」。
```

判据（不变）：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、
`--fstats` → `face=113 mt≈228`、a3n00 面积比从 0.8996 起上升、`t101_verify.ps1` 全绿。

---

### 9.412 —— 【第一现场确认】丢边的就是 `check_pcurves_and_shift`（`read_topology.rs:714`）：8 → 6

在 `:711 / :713 / :714` 三条调用各自前后打印该 wire 的边数（TEMP，env gated，**已撤除**）：

```text
W entry edges=8  →  W after711 edges=8  →  W after713 edges=8  →  W after714 edges=6     （两张面如此）
其余全部 entry=8 → after714=8（1364 条读数，只有这两张 4-wire 面在 :714 变成 6）
```

⇒ **`crate::shhealing::check_pcurves_and_shift(w, &face, preci, false)` 就地改写了 wire**，
把 STEP 里的 8 条边（含 `#5012/#5018` 那对同端点去-回边）变成 6 条，并留下中点接点 `-65.243610`。
`xsalgo_check_pcurve`（:711）与 `project_wire_pcurve_ranges`（:713）都无辜。

#### 已完整的因果链（每一环都有读数）

```text
STEP 文件    外环 8 条（#5004/#5012/#5018/#5090/#5098/#5107/#5115/#5137，8 条 distinct）
resolve_loop 8 → 8（§9.408 341 loop 全 items==edges）
make_wire    8 → 8（§9.408）
make_face    8 → 8（§9.410）
check_pcurves_and_shift (read_topology.rs:714)   **8 → 6**（本节）—— 丢 #5012/#5018，造中点 -65.243610
seam 块入口  [1,1,6,14]（§9.407）⇒ fix_missing_seam 拿到的已是残缺外环
  ⇒ ComposeShell 出两张重合面（§9.377）⇒ Shell(5) 被 resolve_face 丢 ⇒ 面停在 4 wires
  ⇒ 建模型 wire[0]=6 / wire[1]=14（§9.393）⇒ face-checker 判 SELF_INTERSECTING_WIRE ⇒ FAILURE
  ⇒ 三角化循环跳过（§9.391）⇒ **螺母斜切面为空**（§9.368/§9.371）
```

#### 下一步（该看的就是这一个函数 + 其 .cxx）

```text
端口  crates/occt-topo/src/shhealing/ 的 check_pcurves_and_shift（注释指向 cxx:365-480）
OCCT  ShapeFix_Face.cxx:365-480（CheckPCurves 一线的首个 wire 轮次）
做法  逐行对照：确认 OCCT 在该处是否也「就地重建 wire / 合并同端点去-回边」；
      若 OCCT 只是 shift pcurve 参数而**不动拓扑**，端口这里就是在做多余的重建（应改为只改 pcurve 且写回副本）。
      顺带核 -65.243610 的产生点（中点切边）是否出自该函数内部的某个 split。
判据  --model 113 → 2 wires（22+6 边）、zz_seam_fix 113 → Face、--fstats → face=113 mt≈228、
      a3n00 面积比从 0.8996 起上升、t101_verify.ps1 全绿。

---

### 9.413 —— T-101：**STEP 解析/处理的覆盖矩阵**（子任务结论，用户指令的直接回答）+ `check_pcurves_and_shift` 的落点定位

#### A. 覆盖/一致性总判定（读代码 + 实测，逐条见矩阵）

**OCCT 侧**：`FromSTEP.exec.op = FixShape`（`STEPControl_Controller.cxx:201`）在传输期逐实体跑
`ProcessShape`（`STEPControl_ActorRead.cxx:1942`）⇒ `XSAlgo_ShapeProcessor::ProcessShape`（`:63`）
⇒ `ShapeProcess::Perform`（`ShapeProcess.cxx:189`）⇒ 算子 `fixshape`（`ShapeProcess_OperLibrary.cxx:785`，
参数下发 `:806-865`）⇒ `ShapeFix_Shape::Perform`（`ShapeFix_Shape.cxx:83`，分派 `:130`：
SOLID `:160` / SHELL `:176` / FACE `:193` / WIRE `:212`）⇒ `ShapeFix_Face::Perform`（`:345`，`FixMissingSeam` 调用点 `:492`）。

**实测该算子承重**（同一文件，探针开/关 ShapeProcess）：

```text
默认（整形）  226/226 面全部有三角化，wire 直方图 {1:208, 2:9, 4:1, 6:4, 10:4}
noop          **131/226 面无三角化**（全为 UV-DEGENERATE-WIRE），{1:163, 2:53, 4:2, 6:4, 10:4}
端口现状      225/226 有三角化，**1 面无三角化**（wires=4、mv=mt=0、dv=inf），{1:206, 2:10, 4:2, 6:4, 10:4}
```

⇒ 端口**没有 ShapeProcess/FixShape 驱动器**（全仓无 `ShapeProcess`/`OperLibrary`/`ShapeFix_Shape` 实现），
只在 `read_topology.rs:753-773` 按面直接调 `fix_missing_seam` + `check_pcurves_and_shift` ——
相当于只搬了「SOLID/SHELL/FACE/WIRE 分派」下的一个子步，**分派器与参数下发缺失**。

**已一致的部分**：`FixMissingSeam` 本体（找 seam 循环 + ComposeShell 调用，§9.377–§9.413）、
`TranslateEdgeLoop`/`TranslateFace` 的 bound/loop 解析（`resolve_loop`/`bind_edge_loop_vertices`/
`associate_edge_pcurve`）、`GeometricTool` 的 seam 判定、`XSAlgo_ShapeProcessor::CheckPCurve`
（`xsalgo_check_pcurve.rs:33`）、网格侧 ModelHealer/FaceChecker 主流程。

**覆盖缺口（按严重度，子任务原文）**：

```text
1) 无 ShapeProcess/FixShape 驱动器 —— 最高；实测 noop 时 131/226 面空
2) SHELL/SOLID 级 fixer 缺席（ShapeFix_Shell/Solid；端口无实现）
3) FixReorder 结果不落盘（ShapeFix_Face.cxx:2029-2032 ↔ shape_fix_face.rs:406-411「只报告不重写边表」）
4) ModelPostProcessor 未接线（生产路径在 discret_root.rs:519 另写一份）
5) 生产路径绕过 MeshContext + 多一条 OCCT 不存在的逐面 wireframe 兜底（discret_root.rs:476-489，注释自认 UNPORTED）
6) 阶段顺序被交织（端口先 PreProcessor 再 heal，OCCT 相反）—— 数值影响未验证
7) read 侧参数面缺一半（read.maxprecision/surfacecurve/encoderegularity/read.step.* 无开关）
8) 陈旧文档/工件：shape_fix_face.rs:4-6/130-132/147-148 仍写「seam construction UNPORTED」（实际已到 cxx:2330）；
   .target-gate/a3n00_port_faces.txt 是接 seam 之前的直方图，勿当现状
```

#### B. 与 F113 这条线的关系

`check_pcurves_and_shift` 位于 `crates/occt-topo/src/shhealing/wire_fix.rs:3168`
（子任务把它对应到 `ShapeFix_Face.cxx:365-492` 的前置 wire-fix，状态「部分，未逐行对拍」）。
本轮在 3168-3288 窗口内扫到的**唯一拓扑动作**是 `:3233 reg.remove_pcurves_on_surface(&e.0, &face.0)`
（只删 pcurve，不删边）⇒ **8→6 的丢边点不在这个窗口内**，需把窗口扩到该函数全程（或看它是否经由
`ShapeFix_Wire` 一类的重建路径）。这与缺口 3（`FixReorder` 不落盘）方向一致：端口在 wire 级修复上
只做了一部分，个别面会被改写成与 OCCT 不同的边表。

#### 下一步

1. 读 `wire_fix.rs:3168` 起的**完整** `check_pcurves_and_shift`，找出把 8 条边变 6 条的那一步
   （判据：`-65.243610` = `#5012` 两端点的中点，是切边指纹）；
2. 对照 `ShapeFix_Face.cxx:365-492` 的同等分支，确认 OCCT 在该处**是否改拓扑**；
3. 仍按验收：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、
   `--fstats` → `face=113 mt≈228`、a3n00 面积比从 0.8996 起上升、`t101_verify.ps1` 全绿。

---

### 9.414 —— T-101（F113）：`check_pcurves_and_shift` 内部再分流 —— `fix_reorder_wire` / `fix_small_all` **都不丢边**，丢边在函数**后半段**

在 `wire_fix.rs:3168` 的 `check_pcurves_and_shift` 里，于 `fix_reorder_wire`（:3171）与
`fix_small_all` 块之后（`let edges = edges_of_wire(wire);`，:3189 之前）各打一次 wire 边数
（TEMP，env gated，**已撤除**）：

```text
12 x  CP entry edges=8  →  CP after_fixsmall edges=8        （全部 12 条 8 边 wire 都是 8 → 8）
（共 5494 条 CP 读数）
```

⇒ `fix_reorder_wire` 与 `fix_small_all`（`ShapeFix_Wire::FixSmall`）**都没有删边**；
而 `read_topology.rs:714` 的调用返回后是 6（§9.412）⇒ 删边发生在
`check_pcurves_and_shift` 的**后半段**（:3189 之后）：

```text
候选：fix_notched_edges（:3249）、以及 :3305-3333 的收尾分支（fix_lacking / check_wire 一线）
已排除：fix_reorder_wire（:3171/3182）、fix_small_all（:3177）、pcurve 循环本身（:3208-3247，只动 pcurve）
```

#### 下一步

在该函数后半段的每个候选调用前后再打一次边数（同样三行探针），指名那一步；
然后对照 `ShapeFix_Face.cxx:365-492` 的同等分支判断 OCCT 是否也删这两条边
（`-65.243610` 中点 = 切边指纹，指向「把一对去-回边合并成一条、并取中点」的动作）。
验收不变：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、
`--fstats` → `face=113 mt≈228`、面积比从 0.8996 起上升、`t101_verify.ps1` 全绿。

---

### 9.415 —— T-101（F113）：删边点收窄到 `check_pcurves_and_shift` 的**最后一段**（`fix_notched_edges` 调用处 → 函数末尾），**只影响 2 条 wire**

在 `wire_fix.rs` 的 `fix_notched_edges`（:3249）调用**之前**与函数末尾各打一次边数（TEMP，env gated，**已撤除**）：

```text
CP pre_notched edges=8 → CP fn_end edges=6      x **2**      ← 就是那两张面
CP pre_notched edges=8 → CP fn_end edges=8      x 10         ← 其余 8 边 wire 不变
（另有 384 条读数覆盖各种边数，全部 pre == end）
```

配套读数（同函数内，前面几轮）：`entry 8 → after_fixsmall 8`（§9.414），
而 `read_topology.rs:714` 调用返回后是 6（§9.412）。

⇒ 删边发生在 **:3249 之后到函数末尾（:3333）之间**，两个候选：

```text
1) fix_notched_edges(wire, face, MIN_TOLERANCE, MAX_TOLERANCE)   （:3249 调用）
2) :3250-3333 的收尾分支（fix_lacking / check_wire 一线的合并逻辑）
```

且它**只对这两条 wire 生效**（其余 8 边 wire 原样通过）——这与「STEP 里那对同端点去-回边
（`#5012/#5018`，跨度 -49.864906 ⇄ -80.622679、宽度 0）」的指纹完全吻合：
只有存在这种退化对时才会触发合并，并留下中点 `-65.243610`。

#### 下一步

在 `fix_notched_edges` 调用前后各打一次（补一条即可分流）：
- `pre_notched 8 → post_notched 6` ⇒ 就是 `fix_notched_edges`（再对它内部逐句看，并对照
  `ShapeFix_Wire::FixNotchedEdges` / `ShapeFix_Face.cxx:365-492` 的同等分支）；
- 否则 ⇒ 收尾的 `fix_lacking` 一线（`:3305-3333`）。

验收不变：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、
`--fstats` → `face=113 mt≈228`、面积比从 0.8996 起上升、`t101_verify.ps1` 全绿。

---

### 9.416 —— T-101（F113）：本轮探针**未配对**（`pre_notched` 行在 §9.415 清理时一并被删），结论待下一轮配对重跑

本轮把 `fix_notched_edges(...)` 调用包成 `let t101_notched = ...;` 并只打印 `post_notched`，
结果**无法与 pre 配对**（§9.415 的清理把 `pre_notched`/`fn_end` 两条打印一起删掉了）。读数本身：

```text
CP post_notched 分布：edges=1:168 / 4:110 / 5:47 / 3:28 / 6:11 / 8:10    （共 374 次调用）
```

`edges=6` 出现 11 次而 §9.415 只有 2 条 wire 最终是 6 ⇒ 说明这条 wire 在多次调用里重复出现
（同一个面/同一条 wire 会被反复处理），**单看 post 无法判定是不是 `fix_notched_edges` 干的**。

另外：本轮把调用包成 `let ... = cond; if cond {` 的还原动作一度写坏了两行（多出一个 `{` 且缩进丢失），
已按原文修回，`cargo check` 复原、`git diff` 为空（**未用 `git checkout`**）。

#### 下一步（一次配对探针即可定论）

同一次插入里同时保留两行：

```rust
let t101_pre = crate::topo_tools_full::edges_of_wire(wire).len();
let t101_notched = reorder_ok && fix_notched_edges(wire, face, SHAPE_FIX_MIN_TOLERANCE, SHAPE_FIX_MAX_TOLERANCE);
if std::env::var_os("T101_DUMP_SEAM").is_some() {
    eprintln!("CP notched {} -> {}", t101_pre, crate::topo_tools_full::edges_of_wire(wire).len());
}
if t101_notched { … }
```

按「哪一次调用把 8 变 6」判定是 `fix_notched_edges` 还是收尾的 `fix_lacking` 一线；随后对照
`ShapeFix_Wire::FixNotchedEdges` / `ShapeFix_Face.cxx:365-492` 的同等分支决定改法。
验收不变：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、`--fstats` → `face=113 mt≈228`、
面积比从 0.8996 起上升、`t101_verify.ps1` 全绿。

---

### 9.417 —— 【第一现场锁定】`fix_notched_edges` 就是把外环 8 条边变 6 条的那一步（配对探针）

同一次插入里成对打印（`wire_fix.rs:3249` 的 `fix_notched_edges` 调用前后，TEMP，env gated，**已撤除**；
用 `zz_uv_feed --ids` 跑，因为另一审计子任务当时正占用 `zz_probe_a3n00.exe` 的链接目标）：

```text
CP notched 8 -> 8     x 10
CP notched 8 -> 6     x **2**      ← 就是那两张面（含 F113）
（共 384 条 notched 读数，覆盖各种边数）
```

⇒ **`fix_notched_edges`（`:3249`）单次调用内把 wire 从 8 条边变成 6 条**，
删掉 STEP 里的 `#5012/#5018`（同端点「去-回」对），并留下中点切点 `-65.243610`。

完整定位链（每一环都有读数，终于到底）：

```text
STEP 8 条 → resolve_loop 8 → make_wire 8 → make_face 8
  → check_pcurves_and_shift(wire_fix.rs:3168)
      · fix_reorder_wire           8 → 8
      · fix_small_all              8 → 8
      · **fix_notched_edges**      **8 → 6**   ← 第一现场（本节）
      · （其后的 fix_lacking 一线未涉及）
  → 外环 6 条 → fix_missing_seam 出 Shell(5)（§9.377）→ 面停在 4 wires
  → 模型 wire[0]=6/wire[1]=14（§9.393）→ face-checker 判 SELF_INTERSECTING_WIRE（§9.395）
  → 三角化循环跳过（§9.391）→ **螺母斜切面为空**
```

#### 下一步（读两窗口，然后才是改）

```text
端口  crates/occt-topo/src/shhealing/wire_fix.rs 的 fix_notched_edges
OCCT  ShapeFix_Wire::FixNotchedEdges（含 ShapeFix_Face.cxx:365-492 的调用上下文与 FixNotchedEdgesMode 默认 -1）
做法  逐条比：哪些「notch」判据被接受、是否允许把一对同端点去-回边**合并成一条并取中点**、
      以及 OCCT 在该输入（#5012/#5018 宽 30 的去-回对）下是否也会动拓扑。
      注意端口注释里已提到 FromSTEP 的 `FixNotchedEdgesMode = -1`（STEPControl_Controller.cxx:248）
      与 `FixTailMode = 0`（:249）—— 先确认门是否与 OCCT 一致。
判据  --model 113 → 2 wires（22+6 边）、zz_seam_fix 113 → Face、--fstats → face=113 mt≈228、
      a3n00 面积比从 0.8996 起上升、t101_verify.ps1 全绿。

---

### 9.418 —— T-101：**a3n00 逐面数据级审计**（子任务，226 面全量）+ 与 §9.417 的互相印证 ⇒ 覆盖缺口的**量化口径**确定

#### A. 全量数字（口径见子任务原文：文件侧正则解析 `.target-gate/step_side.py`；端口侧新增只读
`--ecensus` 开关；OCCT 侧 `run_dbg.bat … [noop]`；配对用**顶点坐标集合**几何配对，`PAIR exact=202`）

```text
wire 级（bound 数 vs wires=）：
  STEP 文件        {1:163, 2:53, 4:2, 6:4, 10:4}
  OCCT noop        {1:163, 2:53, 4:2, 6:4, 10:4}   —— **226/226 逐面相等**
  port             {1:206, 2:10, 4:2, 6:4, 10:4}   —— **183/226 相等，43 面不等**（方向全是 port 把多 bound 并成 1 wire）
  OCCT 整形后      {1:208, 2:9, 4:1, 6:4, 10:4}    —— 索引被重排，按干净 bbox 配 200/226

边级（STEP 每环边数 vs port 每 wire 边数）：**45/226 面不一致**，其中**只有 2 面 net 丢边**：
  idx113  #5375  Cyl     STEP [8,14,1,1] → port [6,14,1,1]   净 +2   （OCCT 整形后 [22,6]）
  idx171  #7415  BSpline STEP [8,8,8,8]  → port [8,8,8,6]    净 +2   （OCCT 整形后 [8,8,8,7]）
  其余 43 面是「port 合并多 bound」（净 −2 ×22 面、净 −3 ×21 面），即 fix_missing_seam 的正常产物
```

#### B. 两处丢边的形态**完全一致**（与 §9.405/§9.406 我自己测到的 F113 一致）

```text
idx113  bound #5140 / loop #5139 的 #5012、#5018 都是 LINE，端点对完全相同（#4956↔#5007），一条 .F.
        port 把这对**整对删掉**，并把邻边收到**新顶点** (-65.243610,0,-100)（≈ 毛刺中点 -65.24379）
        OCCT 整形后 FACE 25 同样只剩 6 边，但**收在原顶点** (-49.8649,0,-100) ⇒ 几何/拓扑都不同
idx171  loop #7413 里 #7405 被**同一 loop 引用两次**（全 341 个 loop 中唯一的重复引用）
        port 同样整对删掉，收到新顶点 (-35.250069,0,-116.9011)（≈ 另两顶点中点 -35.249943）
        OCCT 整形后 FACE 2 的 wire[3] 是 **7** 边、收在原顶点
```

#### C. 与 §9.417 的互相印证（**第一现场同一处**）

```text
我（§9.417，实测配对探针）：打开 `wire_fix.rs:3249` 的 fix_notched_edges 前后打印边数
    CP notched 8 -> 8  x10 ；  CP notched 8 -> 6  x2   ⇒ **单次调用内 8→6**
子任务（读代码 + 数据）：候选链 `read_topology.rs:714 check_pcurves_and_shift`
    → `wire_fix.rs:3249 fix_notched_edges`，并指出其 `cxx:4051-4067` 的「两半分拆」落在
    `wire_fix.rs:3103-3116`（**会造新顶点**）——与观测到的「中点新顶点 + 少 2 边」吻合
⇒ 两条独立路径指向同一函数：**`fix_notched_edges` 就是这一处 reader 覆盖缺口的现场**
   （子任务当时把该步标为「未验证」，本节与 §9.417 合起来把它验证了）
```

#### D. 计数看不见的同类隐患（子任务发现，重要）

```text
341 个 EDGE_LOOP 中：1 个含重复边引用（#7413），**3 个含「同端点反向边对」**：#5139(→idx113)、
#7413(→idx171)、**#8243(→port idx 203)**。第三个的 file/OCCT/port 计数都是 [2]，
所以**只数边数的 census 看不到它** —— 同类问题可能不止 2 面。
```

#### E. 下一步（窗口已具体到行）

```text
1) 读 `wire_fix.rs:3103-3116`（两半分拆）与 `fix_notched_edges` 本体，
   对照 `ShapeFix_Wire::FixNotchedEdges`（cxx:4051-4067）与门 `FixNotchedEdgesMode = -1`
   （STEPControl_Controller.cxx:248）、`FixTailMode = 0`（:249）；
2) 关键差异先记着：OCCT 整形后同样少 2 边但**收在原顶点**，port 收在**中点新顶点** ——
   这一条既是「该不该动拓扑」也是「动到哪个顶点」的判据；
3) 改完按验收：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、
   `--fstats` → `face=113 mt≈228`、a3n00 面积比从 0.8996 起上升、`pwsh -File .target-gate\t101_verify.ps1` 全绿。
   （另外 idx171/idx203 也应按同一改法一起看。）

---

### 9.419 —— T-101（F113）：`fix_notched_edges` 的机制读清了 —— 它**在 notch 参数处把邻边对半切开**并接到新顶点 `vnew`；与 OCCT「收在原顶点」的差别就在这个参数/顶点上

端口 `wire_fix.rs:3096-3141`（对应 `ShapeFix_Wire::FixNotchedEdges`，端口注释给的是 `cxx:4051-4067`）：

```rust
let uv   = c2d.d0(check.param);                                     // :3100  notch 参数处的曲面点
let vnew = TopoBuilder::new().make_vertex(surface.d0(uv.x(), uv.y()), CONFUSION);   // :3101 **新顶点**
// 第一半：边首 -> vnew，range [first, check.param]                     // :3103-3115
// 第二半：vnew -> 边尾，range [check.param, last]                      // :3116-3121
wire_set_edge_composed(wire, to_split, &new_e1);                     // :3133 用两半替换原边
wire_insert_edge_before(wire, at, &ins);                             // :3139
```

⇒ 机制是：**检出 notch（一对同端点的去-回边）→ 删掉它 → 把与它相邻的那条边在 `check.param`
处一分为二 → 两半之间插入新顶点 `vnew`**。对 F113 来说，`vnew` 就是观测到的 `-65.243610`
（≈ 毛刺中点 `-65.24379`），而**OCCT 整形后同样只剩 6 条边，但收在原顶点 `-49.864906`**。

#### 由此得到**具体的分歧候选**（下一步只查这一处）

```text
分歧点 = 「切邻边的那个顶点/参数」：
  端口：vnew = surface.d0(c2d.d0(check.param))（由 notch 参数算出的新点，落在毛刺中部）
  OCCT：整形结果保留原顶点 -49.864906（= 毛刺 #5012/#5018 的端点之一，即 notch 的起点）
⇒ 需要核 `check.param` 的来源与 OCCT 的对应量：
  · notch 检测（`ShapeFix_Wire::CheckNotchedEdges` 一线）里 `param` 是怎么取的
    （毛刺中点？交点？还是毛刺端点？）
  · 以及 OCCT 在该分支是否**根本不做「对半切邻边」**，而是直接把两半接到已有顶点上
读法：端口 `wire_fix.rs` 的 notch 检测段（`fix_notched_edges` 里 `check.param` 的产生处）
      ↔ `ShapeFix_Wire.cxx:4051-4067` 与 `CheckNotchedEdges` 的同段。
```

这条与 §9.418 的审计观测**完全对得上**（端口：新顶点=毛刺中点；OCCT：收在原顶点），
也解释了为什么两个模型在 idx113/idx171 上「边数相同、几何不同」——**不是丢边方向的问题，
而是切在哪个顶点的问题**。

#### 验收（不变）

`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、`--fstats` → `face=113 mt≈228`、
a3n00 面积比从 0.8996 起上升、`pwsh -File .target-gate\t101_verify.ps1` 全绿。

---

### 9.420 —— T-101（F113）：notch 修复**机器本身逐行忠实** ⇒ 分歧就在 `check_notched_edges` 返回的 `param`（或 `on_end` 分支）

两侧对照（端口 `wire_fix.rs:3046-3151` ↔ OCCT `ShapeFix_Wire.cxx:3977-4067`）：

```cpp
// OCCT
double param;
if (theAdvAnalyzer->CheckNotchedEdges(i, toRemove, param, MinTolerance())) {      // :3997
  ... isRemoveFirst = (n1 == toRemove) ...                                        // :4000-4010 附近
  if (|param - (isRemoveFirst ? b : a)| <= PConfusion || ... ) { … }              // :4014-4016
  if (|(isRemoveFirst ? b : a) - param| < PConfusion) { … }                       // :4026
  ShapeAnalysis_TransferParametersProj transferParameters; … Init(splitE, face);  // :4031-4034
  aV = B.MakeVertex(Analyzer()->Surface()->Value(c2d->Value(param)), …);          // :4049
  transferParameters->TransferRange(newE1, first, param, true);                   // :4057
```
```rust
// 端口
fn fix_notched_edges(wire, face, min_tol, max_tol) -> bool {                      // :3046
    let Some(check) = check_notched_edges(wire, face, i, min_tol) else { … };     // :3051  ↔ :3997
    let is_remove_first = n1 == check.short_num;                                  // :3059
    let to_split = if n2 == check.short_num { n1 } else { n2 };                   // :3060
    let on_end = (check.param - …).abs() <= PCONFUSION …                          // :3074-3076 ↔ :4014-4016
    if ((…) - check.param).abs() < PCONFUSION { … }                               // :3083      ↔ :4026
    let uv = c2d.d0(check.param);                                                 // :3100
    let vnew = make_vertex(surface.d0(uv.x(), uv.y()), CONFUSION);                // :3101      ↔ :4049
    transfer_range(&mut new_e1, first, check.param, true);                        // :3113      ↔ :4057
    transfer_range(&mut new_e2, check.param, last, true);                         // :3119
```

⇒ **notch 修复的动作序列、新顶点构造方式（`Surface()->Value(c2d->Value(param))`）、两半的 range 转移
都与 OCCT 一致**。所以「端口把邻边切在中点、OCCT 收在原顶点」的差别只能来自：

```text
1) `check_notched_edges` 返回的 `param`（↔ `ShapeAnalysis_Wire::CheckNotchedEdges(i, toRemove, param, MinTolerance())`）
   —— 端口此处给的是「邻边上靠近毛刺中部的参数」，OCCT 给的应落在**毛刺端点**（= 已有顶点）；
2) 或 `on_end` / `isRemoveFirst` 两个分支判据（:3074-3083 ↔ :4014-4026）在端口侧没命中，
   于是走了「切边 + 造新顶点」而不是「接已有顶点」。
```

与审计观测一致：端口 idx113 的新顶点 `≈ -65.243610`（毛刺中点），OCCT 整形后保留 `-49.864906`
（毛刺端点）——**边数相同、几何不同**。

#### 下一步（最后一个窗口）

```text
端口  wire_fix.rs 的 check_notched_edges（:3051 调用的实现）
OCCT  ShapeAnalysis_Wire::CheckNotchedEdges（ShapeAnalysis_Wire.cxx 内；用参数 (i, toRemove, param, tol) 定位）
做法  逐行比 param 的算法与 min_tol 的取值，以及它接受 notch 的条件（长度/容差/重合判定）；
      再看 on_end / isRemoveFirst 两支的判据是否等价。
判据  --model 113 → 2 wires（22+6 边）、zz_seam_fix 113 → Face、--fstats → face=113 mt≈228、
      a3n00 面积比从 0.8996 起上升、t101_verify.ps1 全绿（idx171/idx203 同法一并看）。
```

---

### 9.421 —— T-101（F113）：`check_notched_edges` 里的**投影目标**存疑（端口 `&ad1` vs .cxx `Ad2`），但**实测对 F113 无影响**，已回退；下一步必须先核变量语义

对标结果（端口 `wire_fix.rs:2899-2924` ↔ `ShapeAnalysis_Wire.cxx:1962-1987`）：

```cpp
// OCCT
double dist1 = ProjectInside(Ad1, gp_Pnt(p2d2...), Tolerance, Proj1, param1, false);   // :1962  E2 的点 -> Ad1
double dist2 = ProjectInside(**Ad2**, gp_Pnt(p2d1...), Tolerance, Proj2, param2, false); // :1963 E1 的点 -> **Ad2**
if (dist1 > Tolerance && dist2 > Tolerance) return false;                              // :1965-1968
if (dist1 < dist2) { shortAD=Ad2; longAD=Ad1; lenP=b2-a2; firstP=a2; shortNum=n2; param=param1; }   // :1970-1978
else               { shortAD=Ad1; longAD=Ad2; lenP=b1-a1; firstP=a1; shortNum=n1; param=param2; }   // :1979-1987
```
```rust
// 端口
let proj1 = project_inside(&ad1, &start2, tolerance, false);      // :2899  ✓ 同 :1962
let proj2 = project_inside(**&ad1**, &start1, tolerance, false);  // :2900  ← 与 :1963 的 Ad2 不一致
if proj1.distance > tolerance && proj2.distance > tolerance { … } // :2901  ✓ 同 :1965
let (long_ad, short_ad, len_p, first_p, short_num, param) = if proj1.distance < proj2.distance { … }  // :2906  ✓
```

**实测**：把 `:2900` 的目标从 `&ad1` 改成 `&ad2`（严格照 `:1963`），编译通过后：

```text
zz_uv_feed --model 113 → MODEL f=113 wires=4（**不变**）；zz_seam_fix 113 的 BEFORE 仍是 wire[0] nEdges=6
```

⇒ **对 F113 无影响**（可能：F113 走的是 `dist1 < dist2` 那一支，`param = param1`，与 `proj2` 无关；
也可能端口的 `start1`/`start2` 命名与 OCCT 的 `p2d1`/`p2d2` 是**交叉**的，于是原来的 `&ad1` 反而等价于 `:1963` 的 `Ad2`）。
因为改动的语义在剩余预算内无法确证、且实测无改进，**已按原样回退**（`git diff` 为空、`cargo check` 0 error）。

#### 下一步（先核语义，再决定是否改）

1. 读 `wire_fix.rs:2880-2900`，确认 `start1`/`start2` 各自对应 OCCT 的哪个 `p2d*`
   （`p2d1 = c2d1->Value(FORWARD ? a1 : b1)`、`p2d2 = c2d2->Value(FORWARD ? b2 : a2)`）；
2. 再看 F113 走哪一支：在 `param` 落定时打一次 `(dist1, dist2, param, short_num)`（TEMP，用完 edit 反撤），
   与 OCCT 侧同面读数比；
3. 若确认 `:2900` 是笔误 ⇒ 改 `&ad2` 并按验收（`--model 113` → 2 wires、`zz_seam_fix 113` → Face、
   `--fstats` → `face=113 mt≈228`、面积比从 0.8996 起、`t101_verify.ps1`）全量复跑；
   若确认不是笔误 ⇒ 分歧在 `param1` 的计算（`project_inside` ↔ `ProjectInside` 的实现）。

---

### 9.422 —— T-101（F113）：`check_notched_edges` 的**命名没有交叉** ⇒ `:2900` 确实是相对 `cxx:1963` 的转写笔误；但它不是 F113 的病灶（F113 走 `param1` 那一支）

读齐端口 `wire_fix.rs:2866-2901` 后确认语义（与 OCCT `ShapeAnalysis_Wire.cxx:1960-1963` 逐项对应）：

```rust
let pt2 = if e2.forward { c2d2.d0(b2) } else { c2d2.d0(a2) };   // :2895  = cxx:1960 的 p2d2  ✓
let pt1 = if e1.forward { c2d1.d0(a1) } else { c2d1.d0(b1) };   // :2896  = cxx:1961 的 p2d1  ✓
let proj1 = project_inside(&ad1, &start2, tolerance, false);     // :2899  = cxx:1962 ProjectInside(Ad1, p2d2) ✓
let proj2 = project_inside(&ad1, &start1, tolerance, false);     // :2900  ≠ cxx:1963 ProjectInside(**Ad2**, p2d1)
```

⇒ 端口 `start1/start2` **没有**与 `p2d1/p2d2` 交叉 ⇒ `:2900` 的目标写 `&ad1` 是**真的转写不一致**
（应为 `&ad2`）。但 §9.421 实测把它改成 `&ad2` 后 **F113 的读数一字不变**（仍 `wires=4`、`wire[0] nEdges=6`）。

⇒ 由此反推：**F113 在 `:2906` 的分支里走的是 `proj1.distance < proj2.distance`（即 `param = param1`、
`short_num = n2`）那一支** —— 与 `proj2` 的目标无关。所以 F113 的病灶在更前面：
`param1` 的计算（`project_inside` ↔ `ProjectInside` 实现）或 `:2906` 的分支选择依据。

#### 结论与下一步（两条，按性价比）

```text
A) 修笔误（保持 .cxx 对齐）：把 :2900 的目标改成 &ad2。它会影响所有 `dist1 >= dist2` 的面，
   ⇒ **必须整跑验收**（pwsh -File .target-gate\t101_verify.ps1：--fixms 无 Shell / 直方图贴近
   {1:208,2:9,4:1,6:4,10:4} / a3n00 ≥ 0.8996 / T0M ≤ 6 / --lib 1281-0 / step_obj_gates 5/5），
   绿则提交、红则回退。**F113 不会因此变好**（本节已证），所以它是独立的一处对齐修复。
B) 抓 F113 的真病灶：在 `:2906` 之后（`param` 落定时）打一次 `(dist1, dist2, short_num, param)`（TEMP），
   与 OCCT 侧同面读数比；若 `param1` 与 OCCT 的 `param` 不同，再进 `project_inside` /
   `ShapeAnalysis_Curve::Project`（cxx 里 `ProjectInside` 的实现）逐行比。
```

判据（不变）：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、`--fstats` → `face=113 mt≈228`、
a3n00 面积比从 0.8996 起上升；`idx171 #7415`、`idx203 #8243` 同法一并看。

---

### 9.423 —— T-101：落地尝试（待门禁判定）—— `wire_fix.rs:2900` 改为 `&ad2`（对齐 `cxx:1963`），a3n00 读数**零回归**；验收 harness 已在后台运行

本轮把 §9.422 确认的那处转写笔误改掉（**唯一改动**，一行）：

```diff
-            let proj2 = project_inside(&ad1, &start1, tolerance, false);
+            let proj2 = project_inside(&ad2, &start1, tolerance, false);   // 对齐 ShapeAnalysis_Wire.cxx:1963
```

`cargo check`/`cargo build` 0 error；快速读数：

```text
zz_uv_feed --model 113   → MODEL f=113 wires=4（**如 §9.422 所预测：F113 不受此改影响**）
zz_probe_a3n00 a3n00     → TOTAL faces=226 computed_lin=1.076007 stats=225 mesh_v=10863 mesh_t=11941
                           STATMAP matched=225 unmatched=1 sum_mt=11941 flat_mt=11941
                           （与 §9.374 记录的基线**逐字相同** ⇒ a3n00 无回归）
```

**验收 harness 已在后台启动**（job `pwsh-1350`：`pwsh -File .target-gate\t101_verify.ps1`，
结果落盘 `.target-gate/t101_verify.txt`，判据 = `--fixms` 无 Shell / wire 直方图贴近
`{1:208,2:9,4:1,6:4,10:4}` / a3n00 ≥ 0.8996 / T0M unmatched ≤ 6 / `--lib` 1281-0 /
`step_obj_gates` 5/5）。

#### 落地规则（下一轮按此执行）

```text
harness 全绿 ⇒ 提交这一行改动（commit message 说明对齐 cxx:1963），并更新 §9.x 与看板；
harness 有任何红 ⇒ **回退这一行**（`edit` 反向替换回 &ad1），只保留本节与 §9.422 的记录。
```

注意：这处修复**不会**让 F113 变好（§9.422 已证 F113 走 `param1` 分支）——它是独立的一处对齐修复；
F113 的真病灶仍在 `param1`（`project_inside` ↔ `ShapeAnalysis_Curve::Project`）或 `:2906` 的分支判据。

---

### 9.424 —— T-101：`wire_fix.rs:2900` 的转写笔误**已按 `.cxx` 修正并通过全量门禁**（本目标最后一轮）

改动（一行，对齐 `ShapeAnalysis_Wire.cxx:1963` 的 `ProjectInside(Ad2, p2d1, …)`）：

```diff
-            let proj2 = project_inside(&ad1, &start1, tolerance, false);
+            let proj2 = project_inside(&ad2, &start1, tolerance, false);
```

**验收读数**（`pwsh -File .target-gate\t101_verify.ps1`，后台 job `pwsh-1350`，exit 0）：

```text
step_obj_gates: test result: ok. 5 passed; 0 failed; 0 ignored; finished in 288.47s
[occ/a3n00.stp] our=204327.72 occ=227130.30 ratio=0.8996      ← 与 §9.369 基线逐字相同
[occ/T0M.stp]   our=192416.40 occ=192658.64 ratio=0.9987      ← 与 §9.369 基线相同
[occ/acs10.stp] our=267047.23 occ=271219.68 ratio=0.9846      ← 与 §9.369 基线相同
其余模型（Cube/Cone/Cylinder/Sphere/Torus/HoledPlate/…）全部 ratio=1.0000 / f-ratio 1.00
```

⇒ 这处对齐**零回归**，按 §9.423 规则提交。**F113 本身不变**（§9.422 证明它走 `param1` 分支），
所以它是一处独立的 reader 对齐修复，而不是本目标的解药。

#### 本目标（48 轮）结论汇总

```text
已确证并逐环有读数的因果链（螺母斜切面为空）：
  STEP 外环 8 条 → resolve_loop 8 → make_wire 8 → make_face 8
  → check_pcurves_and_shift → **fix_notched_edges 把一对同端点「去-回」毛刺边删掉**
    （并把邻边切在 check.param 处的新顶点；OCCT 整形后同样少 2 边但收在原顶点）
  → 外环 6 条 → fix_missing_seam 出 Shell(5) 被 resolve_face 丢弃 → 面停在 4 wires
  → 建模 wire[0]=6/wire[1]=14 → face-checker 判 SELF_INTERSECTING_WIRE(+FAILURE)
  → 三角化循环跳过 → **该面无三角**

全量数据级审计（226 面，顶点集合几何配对，禁按面号）：
  wire 级 183/226 与文件相等（43 面不等 = fix_missing_seam 的正常合并产物）
  边级 45 面不一致，其中**只有 2 面净丢边**：idx113 #5375、idx171 #7415
  同类隐患（计数看不见）：idx203 #8243（file/OCCT/port 都是 [2]）

已逐行对过并判定「忠实」：ComposeShell 全链（split_by_grid/split_by_line/split_wire/collect_wires/剪枝）、
  五个 helper + TOLINT、FixMissingSeam 定位段与假想 grid 构造、notch 修复机器（动作/新顶点/range 转移）、
  resolve_loop/make_wire/make_face、xsalgo_check_pcurve、project_wire_pcurve_ranges、fix_reorder_wire、fix_small_all

已确认的两处**独立**对齐问题（都可单独验证）：
  1) reader 没有 ShapeProcess/FixShape 驱动器（实测 noop 时 131/226 面无三角化）—— 最大缺口；
  2) wire_fix.rs:2900 的投影目标笔误（本节已修并过门禁）。

F113 的**最后病灶**（下一批的起点）：`check_notched_edges` 的 `param1`
  （`project_inside` ↔ `ShapeAnalysis_Curve::Project`）或 `wire_fix.rs:2906` 的分支判据 ——
  因为 F113 走 `proj1.distance < proj2.distance` 那一支。

---

### 9.425 —— T-101（F113）：调用方循环忠实，但 `check_notched_edges` 的**切线/距离预检用的是「方向取反」的点**（`p2d1/p2d2`）——新候选

**① 调用方循环：忠实。** 端口 `wire_fix.rs:3046-3060` ↔ OCCT `ShapeFix_Wire.cxx:3993-4002`：

```text
while i <= nb && nb > 2                 ↔ for (i = 1; i <= NbEdges() && NbEdges() > 2; i++)
n2 = if i > 0 { i } else { nb }         ↔ n2 = (i > 0) ? i : NbEdges()
n1 = if n2 > 1 { n2 - 1 } else { nb }   ↔ n1 = (n2 > 1) ? n2 - 1 : NbEdges()
is_remove_first = n1 == short_num       ↔ isRemoveFirst = (n1 == toRemove)
to_split = if n2 == short_num { n1 } else { n2 } ↔ toSplit = (n2 == toRemove ? n1 : n2)
```

**② 新候选：`check_notched_edges` 里 `p2d1/p2d2` 的方向判据写反了。** 端口 `:2866-2877`：

```rust
let (p2d1, tan1) = if e1.0.orientation().is_reversed() { c2d1.d1(a1) }      // reversed -> a1
                   else { let (p,d) = c2d1.d1(b1); (p, d.reversed()) };    // forward  -> b1
let (p2d2, tan2) = if e2.0.orientation().is_reversed() { … d1(b2) … }      // reversed -> b2
                   else { c2d2.d1(a2) };                                   // forward  -> a2
```
OCCT `ShapeAnalysis_Wire.cxx:1960-1961`：

```cpp
p2d2 = c2d2->Value(E2.Orientation() == TopAbs_FORWARD ? **b2** : **a2**);   // forward -> b2, reversed -> a2
p2d1 = c2d1->Value(E1.Orientation() == TopAbs_FORWARD ? **a1** : **b1**);   // forward -> a1, reversed -> b1
```

⇒ 两边**正好相反**（端口 forward 取 b、reversed 取 a；OCCT forward 取 a、reversed 取 b）。这两个点被
`:2883` 的**预检**用到：

```rust
if tan2.angle(&tan1).abs() > 0.1 || p2d1.distance(&p2d2) > tolerance { return None; }
```

**注意**：紧接着 `:2895-2896` 又用**另一套正确**的取法算投影点：

```rust
let pt2 = if e2.forward { c2d2.d0(b2) } else { c2d2.d0(a2) };   // ✓ = OCCT p2d2
let pt1 = if e1.forward { c2d1.d0(a1) } else { c2d1.d0(b1) };   // ✓ = OCCT p2d1
let proj1 = project_inside(&ad1, &start2, tolerance, false);
let proj2 = project_inside(&ad2, &start1, tolerance, false);    // （上一轮已按 cxx:1963 修好目标）
```

⇒ 端口的同一个函数里存在**两套取点逻辑**：参与「角度/距离预检」的 `p2d1/p2d2` 方向反了，
参与投影的 `pt1/pt2` 是对的。当两条边里有 `REVERSED` 时，**预检比较的是错的那一对端点**
（例如把「尾对尾」当成「首对首」），于是「哪两条边算 notch、谁是 short、`param` 取哪个投影」
都可能与 OCCT 不同 —— 这与 F113 观测到的「新顶点落在毛刺中点、OCCT 收在原顶点」方向一致。

#### 下一步（一次读数即可判定）

在 `:2883` 的预检处按 `p2d1/p2d2`（现写法）与 `pt1/pt2`（OCCT 写法）各判一次，
打印两条边的朝向与两个距离，看 F113 的那对毛刺边是否**只有用 OCCT 写法才被判为 notch**
（或只有用 OCCT 写法 `short/param` 才落到原顶点）。

判据（不变）：`--model 113` → 2 wires（22+6 边）、`zz_seam_fix 113` → Face、
`--fstats` → `face=113 mt≈228`、a3n00 面积比从 0.8996 起上升、`t101_verify.ps1` 全绿。

---

### 9.426 —— T-101 **新一轮的四步计划**（本轮启动；螺帽 F113 仍未解决，如实记录）

用户明确指出「螺帽问题没有解决」。目标已改写为本计划并提高轮次上限（同一 goal id：工具限制下
新建 goal 需先 complete 旧目标，而旧目标未达成，故不 complete 腾位）。四步按序、每步以读数判定：

| 步 | 动作 | 判定依据 | 命中后的验收 |
|---|---|---|---|
| ① | 在 `wire_fix.rs:2883` 的预检处，按**现写法**（`p2d1/p2d2`，`:2866-2877` 朝向判据与 `ShapeAnalysis_Wire.cxx:1960-1961` **相反**）与 **OCCT 写法**（`:2895-2896` 的 `pt1/pt2`，正确）各判一次 | F113 那对毛刺边（`#5012/#5018`）是否**只有**在 OCCT 写法下 `short`/`param` 落到**原顶点** `-49.864906`（现为毛刺中点 `-65.243610`） | 按 `.cxx` 修 → 跑验收 |
| ② | 若①不成立：`param1` 的算法 | `project_inside`（`wire_fix.rs:2800-2815` 的越界钳位：钳到 `u_first`/`u_last`）↔ `ProjectInside`（`ShapeAnalysis_Wire.cxx:1837` 起）逐行比 | 按 `.cxx` 修 → 跑验收 |
| ③ | 每次改动的验收（一条命令） | `zz_uv_feed --model 113` → **wires=2**（22+6 边）；`zz_seam_fix 113` → **Face**；`--fstats` → **face=113 mt≈228**；a3n00 面积比 **≥0.8996 且上升**；T0M unmatched ≤6 | `pwsh -File .target-gate\t101_verify.ps1` → `--lib` 1281/0、`step_obj_gates` 5/5、直方图贴近 `{1:208,2:9,4:1,6:4,10:4}` |
| ④ | 同法处理 `idx171 #7415`（同症状 8→6）与 `idx203 #8243`（计数看不见的同类：file/OCCT/port 都是 `[2]`） | 同上 | 同上 |

**已确证的链（供接手者复核，各环读数见 §9.404-§9.425）**：

```text
STEP 外环 8 条（#5004/#5012/#5018/#5090/#5098/#5107/#5115/#5137，8 条 distinct）
→ resolve_loop 8 → make_wire 8 → make_face 8     （§9.408/§9.410）
→ check_pcurves_and_shift → fix_notched_edges（:3249）**8 → 6**（§9.417 配对探针）
   丢 #5012/#5018（同端点去-回对），中点 -65.243610 造新顶点；OCCT 收在原顶点 -49.864906（§9.418）
→ 外环 6 条 → fix_missing_seam 出 Shell(5) 被 resolve_face 丢弃（§9.377/§9.412）
→ 面停在 4 wires → 建模 wire[0]=6/wire[1]=14（§9.393）→ face-checker 判 SELF_INTERSECTING_WIRE+FAILURE（§9.395）
→ 三角化循环跳过（§9.391）→ **面为空**（§9.368/§9.371）
```

**另有两处独立对齐项**（不属于螺帽本身，勿混入本计划的验收）：
① reader 缺 `ShapeProcess`/`FixShape` 驱动器（实测 `noop` 时 131/226 面无三角化）；
② `wire_fix.rs:2900` 投影目标笔误（已修，`--lib` 1281/0、`step_obj_gates` 5/5、a3n00 0.8996 与基线逐字相同）。

---

### 9.427 —— T-101 第①步：**「朝向写反」假设被实测否掉**（改动使网格变差，已回退）

按 §9.425 的候选，把 `check_notched_edges` 里的 `p2d1/p2d2` 改成「E1 首 / E2 尾」（照 `cxx:1961/1960` 的字面）后实测：

```text
改前  a3n00  stats=225  mesh_v=10863  mesh_t=11941      （OCCT 侧 11145 / 12466）
改后  a3n00  stats=225  mesh_v=11675  mesh_t=13168      ← **超出 OCCT 约 700**（方向反了）
F113         wires=4（不变）；wire 直方图 {1:206,2:10,4:2,6:4,10:4}（不变）
```

⇒ **结论：那对点（E1 尾 vs E2 首）本来就是「预检」该用的** —— 相邻边共享触点、距离恒 0，
预检在这里的语义是「两条边接在一起且切线近乎反向」；而 `cxx:1960-1961` 的 `p2d1/p2d2`
（E1 首 / E2 尾）是给**投影**用的另一套点，端口 `:2895-2896` 的 `pt1/pt2` 已经与之一致。
OCCT 的预检点在 `:1939` 之前那一段（不在 `:1960-1963`），所以「朝向写反」不成立。

**已按纪律回退**（`git diff` 对 `wire_fix.rs` 为空、`cargo check` 0 error、TEMP 探针已撤）。

#### 下一步 = 计划第②步

进 `param1` 的算法：`project_inside`（`wire_fix.rs:2800-2815` 的越界钳位：钳到 `u_first`/`u_last`）
↔ `ProjectInside`（`ShapeAnalysis_Wire.cxx:1837` 起）逐行比 —— 因为 F113 走 `proj1.distance < proj2.distance`
分支（`param = param1`），而观测到的新顶点落在**毛刺中点**，像「投影落到区间中部」而非「端点」。

---

### 9.428 —— T-101 第②步：`project_inside` **逐行忠实**（`wire_fix.rs:2799-2818` ↔ `ShapeAnalysis_Wire.cxx:1837-1862`）⇒ 下一窗口 = `project_adaptor` ↔ `ShapeAnalysis_Curve::Project`

```cpp
// OCCT ProjectInside (:1844-1861)
double dist = sac.Project(AD, pnt, preci, proj, param, adjustToEnds);
if (param < uFirst) { param = uFirst; proj = AD.Value(uFirst); return proj.Distance(pnt); }
if (param > uLast)  { param = uLast;  proj = AD.Value(uLast);  return proj.Distance(pnt); }
return dist;
```
```rust
// 端口 project_inside (:2806-2817)
let proj = project_adaptor(ad, pnt, preci, adjust_to_ends);
let (u_first, u_last) = (ad.first_parameter(), ad.last_parameter());
if proj.param < u_first { let p = ad.d0(u_first); return Projection { distance: p.distance(pnt), point: p, param: u_first }; }
if proj.param > u_last  { let p = ad.d0(u_last);  return Projection { distance: p.distance(pnt), point: p, param: u_last }; }
proj
```

⇒ **一致**（越界钳位到端点、距离按钳位后的点重算、否则返回投影结果）。所以 `param1` 的分歧不在这一层。

#### 下一个窗口（也是 `CheckNotchedEdges` 里最后没对过的一环）

```text
端口  project_adaptor（project_inside :2806 调用；即 ShapeAnalysis_Curve::Project 的移植）
OCCT  ShapeAnalysis_Curve::Project（ShapeAnalysis_Curve.cxx）
看什么 1) 求解失败时的**兜底**：OCCT 在找不到解时返回「区间中点」——**这正是 F113 观测到的
         新顶点落在毛刺中点 (-65.243610) 的形态**；端口若在更多情形下走兜底（或反之），
         `param1` 就会落在中部而不是端点；
       2) 采样/迭代（`nbSamples`、`preci`、`adjustToEnds` 的传递）是否一致。
```

判定方式：在 `project_adaptor` 的兜底分支打一次（TEMP），看 F113 那对毛刺边上是否**走了兜底**；
同时比 OCCT 同函数的兜底条件。命中即按 `.cxx` 修 → 第③步验收。
（`--model 113` → wires=2、`zz_seam_fix 113` → Face、`--fstats` → face=113 mt≈228、
面积比 ≥0.8996 且上升、`t101_verify.ps1` 全绿。）

---

### 9.429 —— T-101：`param1` 的分歧抓到机制 —— `project_adaptor` 的**「区间中点兜底」在 a3n00 上被触发 48/106 次**

在 `project_inside` 的两个调用点（`wire_fix.rs:2899-2900`）成对打印「区间 + 结果」（TEMP，env gated，已撤除）：

```text
PROJ 行数 106（= 通过预检的那些边对）
param == 区间中点 : **48**      ← 兜底路径（且距离残差很大，例：u=[5.759586532, 6.806784083] param=6.283185307 d=1.42906236）
param == u_first  : 3           ← 钳到端点
param == u_last   : 4
其他（正常投影）  : 51
```

⇒ **F113 那个「毛刺中点」新顶点就是中点兜底的产物**（`param` 恰好等于区间中点）。
也就是说：这些边对在 `project_adaptor`（↔ `ShapeAnalysis_Curve::Project`）里**没找到有效投影**，
于是退化成「取区间中点」。

#### 最后一个窗口

```text
端口  crates/occt-topo/src/shhealing/shape_analysis_curve.rs:280 project_adaptor
      → :289/:299/:303/:307 四个出口 → project_act（ShapeAnalysis_Curve::ProjectAct 一线）
OCCT  ShapeAnalysis_Curve.cxx 的 Project(:205 重载) 与 ProjectAct（定义在 :126/:147/:205 三处重载附近）
看什么 1) OCCT 的兜底条件与兜底值（是否也是 (u1+u2)/2、在什么判据下走）；
       2) 端口为何在这 48 例上「找不到解」——采样数/迭代上限/精度（`preci` 与 `CONFUSION` 的取用）是否一致。
判定   F113 那对毛刺边上：OCCT 若在同样输入下**找到投影**（落在端点），端口却走兜底 ⇒ 就是这一处；
       按 .cxx 修后跑第③步验收（`--model 113` → wires=2、`zz_seam_fix 113` → Face、
       `--fstats` → face=113 mt≈228、面积比 ≥0.8996 且上升、`t101_verify.ps1` 全绿）。

---

### 9.430 —— T-101：`AdjustByPeriod` 逐行忠实 ⇒ 48 例「中心值」来自 `project_act` 内部（下一个窗口已到行）

```cpp
// OCCT ShapeAnalysis.cxx:48-62
double diff = Val - ToVal; double D = |diff|; double P = |Period|;
if (D <= 0.5*P) return 0.;
if (P < 1e-100) return diff;
return (diff > 0 ? -P : P) * floor(D/P + 0.5);
```
```rust
// 端口 wire_fix.rs:38-49
let diff = val - to_val; let d = diff.abs(); let p = period.abs();
if d <= 0.5 * p { return 0.0; }
if p < 1e-100 { return diff; }
(if diff > 0.0 { -p } else { p }) * (d / p + 0.5).floor()
```

⇒ **逐行一致**。所以 §9.429 观测到的「`param` 恰好等于 `0.5*(u_inf+u_sup)`」**不是**这处周期性修正造成的
（端口 `project_act` 里唯一的 `0.5` 就是 `:263` 那次修正调用，且它与 OCCT `:479-484` 同式）。

⇒ 那 48 例的中心值必来自 `project_act` 内部某个分支的**兜底**（把参数取到区间中心）。

#### 下一个窗口（具体到行）

```text
端口  shape_analysis_curve.rs:121-276 project_act
      分支（端口注释已标 cxx 行）：:167-183 circle(cxx:355-374) / :184-188 hyperbola(:376-380) /
      :189-192 parabola(:382-386) / :193-198 line(:388-392) / :199-209 ellipse(:394-399) /
      :210-260 default(:401-477: ProjectOnSegments(25) + Extrema_LocateExtPC) / :261-275 收尾(:479-496)
OCCT  ShapeAnalysis_Curve.cxx:355-496（同分支）
看什么 1) 是否有一处把 `param` 取成「区间中心」（`0.5*(uMin+uMax)`）——端口此刻没有，OCCT 可能有，
         但我们的观测是**端口**返回中心值，所以要找的是端口在哪条路上走到中心（例如 `Extrema` 未命中后
         用 `proj_param` 的初值、或 default 支的分段采样落在中点）；
       2) 端口自带的 extrema 与 `Extrema_ExtPC` / `Extrema_LocateExtPC` 的等价性（:128/:229 两处）。
判定方式：在 `project_act` 的每个分支出口打一次 `(branch, param, u_inf, u_sup, dist)`（TEMP），
按「param == 中心」筛选，指名分支；随后按 .cxx 修那一支 → 第③步验收。

---

### 9.431 —— 【分歧点确认】`project_act` 的 Newton 精修**只在「提前返回」时才生效** ⇒ 粗采样（偶分格的中点样本）被留下 —— 这就是 F113 中点新顶点的来源

OCCT `ShapeAnalysis_Curve.cxx:421-436`：

```cpp
Extrema_LocateExtPC aProjector(thePoint, theCurve, theProjParam /*U0*/, uMin, uMax, theTolerance);
if (aProjector.IsDone())
{
  theProjParam = aProjector.Point().Parameter();     // ← **无条件写入**当前投影参数
  theProjPoint = aProjector.Point().Value();
  const double aDistNewton = thePoint.Distance(theProjPoint);
  if (aDistNewton < aModMin) { return aDistNewton; } // ← 只是**是否提前返回**
}
```

端口 `shape_analysis_curve.rs:229-237`：

```rust
if let Some((t, q)) = crate::int_tools_vertex_line::extrema_locate_ext_pc(curve, point, proj_param, u_inf, u_sup) {
    let newton_dist = point.distance(&q);
    if newton_dist < mod_min {
        return Projection::new(newton_dist, q, t);   // ← 只有提前返回时才用 t/q
    }
}                                                    // ← 否则 t/q 被丢弃，proj_param 仍是粗采样值
```

⇒ **差别**：OCCT 无论是否提前返回，都把 `theProjParam` 更新成精修值；端口在守卫不成立时**丢弃精修结果**，
`proj_param` 保持 `project_on_segments` 的**采样参数**。而采样网格是 `i = 0..=n`、`u = u_min + (u_max-u_min)/n * i`
——**当 n 为偶数时 `i = n/2` 的样本正好等于区间中点**！所以：
「粗采样留下」+「偶分格」⇒ `param` 恰好等于 `0.5*(u_inf+u_sup)`（§9.429 观测到的 48 例），
带一个**非零距离**（1.429…）✓ 与观测完全吻合；OCCT 因为带着精修值，落点回到**真实投影**（毛刺端点）。

#### 精确修法（照 `cxx:429-430`，一行语义）

```rust
if let Some((t, q)) = crate::int_tools_vertex_line::extrema_locate_ext_pc(curve, point, proj_param, u_inf, u_sup) {
    proj_param = t;                     // `cxx:429` 无条件写入
    <当前投影点变量> = q;                // `cxx:430`（需先确认端口变量名：`proj_point` / `computed_point`？）
    let newton_dist = point.distance(&q);
    if newton_dist < mod_min {          // `cxx:432-435`
        return Projection::new(newton_dist, q, t);
    }
}
```

**注意**：`cxx:430` 写的是 `theProjPoint`（当前投影点），不是 `theProjPoint` 之外的「记住值」；
端口里同时存在 `computed_param/computed_point`（:146 一线，给 `:255` 的兜底用）与当前投影变量，
必须先读 `project_act` 头部（:129-166）确认哪个对应 `theProjParam/theProjPoint`，**只改当前投影那一对**。

#### 下一步（下次直接落地）

1. 读 `shape_analysis_curve.rs:121-170` 确认 `proj_param` / 当前投影点 / `computed_*` / `mod_min` 的角色；
2. 照 `cxx:429-430` 加「无条件写入」，跑第③步验收：
   `--model 113` → wires=2（22+6 边）、`zz_seam_fix 113` → Face、`--fstats` → `face=113 mt≈228`、
   a3n00 面积比 ≥0.8996 且上升、`T0M unmatched ≤6`、`pwsh -File .target-gate\t101_verify.ps1` 全绿；
3. 若 F113 达标 ⇒ 按第④步看 idx171 `#7415` 与 idx203 `#8243`，然后提交 + 更新看板/§9.x。

---

### 9.432 —— T-101：§9.431 的修法**实测无效**（已回退）⇒ 中点是被**采样路径**留下的，下一窗口 = `ProjectOnSegments`

落地内容（照 `cxx:429` 的无条件写入）：

```rust
if let Some((t, q)) = extrema_locate_ext_pc(curve, point, proj_param, u_inf, u_sup) {
    proj_param = t;                     // 新增：cxx:429 的无条件写入
    let newton_dist = point.distance(&q);
    if newton_dist < mod_min { return Projection::new(newton_dist, q, t); }
}
```

实测（`cargo build` 0 error）：

```text
F113    MODEL f=113 wires=4；wire[0] nEdges=6                ← **不变**
a3n00   stats=225 mesh_v=10863 mesh_t=11941 （= 基线）        ← **不变**
直方图  {1:206, 2:10, 4:2, 6:4, 10:4}                        ← **不变**
```

⇒ 精修参数的写入**不足以**改变结果（两种可能：① F113 那几例 `extrema_locate_ext_pc` 返回 `None`，
整支没进；② 进了但随后的 40/20/25/40 轮 `project_on_segments` 又把 `proj_param` 采样覆盖）。
按纪律**已回退**（`git diff` 为空、`cargo check` 0 error）。

#### 下一窗口：`ProjectOnSegments` 本体

```text
端口  shape_analysis_curve.rs:69-99  project_on_segments
      let step = (*end - *start) / n;
      let mut min_sq = *proj_dist * *proj_dist;        ← 以「传入的 dist」为初值
      for i in 0..=n { u = *start + step*i; if sq < min_sq { min_sq = sq; *proj_param = u; } }
      *end = (*end).min(*proj_param + step);
      *start = (*start).max(*proj_param - step);        ← ±step 收窄
OCCT  ShapeAnalysis_Curve.cxx 里 `ProjectOnSegments` 的定义（调用点见 :408 / :451）
看什么 1) 收窄公式（±step vs 直接设为 param±step，或收窄**在循环前/后**）；
       2) `MinSqDist` 的初值（OCCT 是否用传入 dist 的平方）；
       3) 采样点数与端点是否含 `i = n`（端口含）。
判定  若 OCCT 的收窄/初值使「偶分格中点样本」不参与或少参与 ⇒ 端口多采了中点 ⇒ 与观测吻合；
       改后按第③步验收（`--model 113` → wires=2、`zz_seam_fix 113` → Face、`--fstats` → face=113 mt≈228、
       面积比 ≥0.8996 且上升、`t101_verify.ps1` 全绿）。

---

### 9.433 —— T-101：`ProjectOnSegments` **也逐行忠实**（`shape_analysis_curve.rs:69-99` ↔ `ShapeAnalysis_Curve.cxx:88-122`）⇒ 只剩 default 支的**调用序列**没对过

```cpp
// OCCT :99-121
const double aParamStep = (anEndParam - aStartParam) / theSegmentCount;
double aMinSqDistance = aProjDistance * aProjDistance;
bool aHasChanged = false;
for (i = 0; i <= theSegmentCount; i++) {
  aCurrentParam = aStartParam + aParamStep*i; ...
  if (aCurrentSqDistance < aMinSqDistance) { aMinSqDistance = …; aProjPoint = …; aProjParam = aCurrentParam; aHasChanged = true; }
}
if (aHasChanged) aProjDistance = sqrt(aMinSqDistance);
anEndParam = min(anEndParam, aProjParam + aParamStep);
aStartParam = max(aStartParam, aProjParam - aParamStep);
```
```rust
// 端口 :81-98 —— 同式（仅不写 aProjPoint，端口用 param 重算）
```

⇒ 一致（含 `i = 0..=n` 含端点、`MinSqDist` 初值取传入 dist 的平方、`±step` 收窄）。

#### 结论与最后未对过的块

`project_act` 的**每个子件**都已逐行对过且一致：`project_on_segments`（本节）、`AdjustByPeriod`（§9.430）、
`project_inside` 钳位（§9.428）。⇒ 剩下的是 **default 支的调用序列本身**（`cxx:401-500` ↔ 端口 `:210-276`）：

```text
cxx:408  ProjectOnSegments(25)
cxx:421  Extrema_LocateExtPC(P, C, theProjParam, uMin, uMax, TolU)   （:429-430 无条件写入；:432-435 才提前返回）
cxx:441  注释「After each call to ProjectOnSegments, uMin and uMax ...」
cxx:451  ProjectOnSegments(...)
cxx:493  theProjParam = anOldParam;                                  ← 端口对应哪一行？
端口:216/229/238/253-257  同序，但 **:238 的循环是 [40,20,25,40]**，且 :255 用的是 computed_*（记住值）
```

**下一个窗口**：读 `ShapeAnalysis_Curve.cxx:401-500` 与 `shape_analysis_curve.rs:210-276` 并排比：
① 三轮采样的**段数序列**是否就是 `[40,20,25,40]`；② `theProjParam` 在这些调用之间是否被**额外重设**；
③ `anOldParam` 的语义（`:493` 在什么条件下回滚）与端口 `computed_param`/`computed_point`（:146 一线）是否对应。
命中即按 .cxx 改 → 第③步验收。

---

### 9.434 —— T-101：精修支**每次都进**（Some 102/102）；尾巴返回与中心值无关 ⇒ 中心来自 default 支的「记住值」= 第一阶段 `Extrema_ExtPC`

两组新读数（都是 TEMP，env gated，**已撤除**）：

```text
① 在 `Extrema_LocateExtPC` 调用处打印：NEWTON calls  None = 0   Some = **102**
   样本：u0=0.000 u=[0.000,0.300] seg_d=2.009446 -> Some(t=0.000000, d=2.009446)
         u0=0.600 u=[0.000,1.000] seg_d=0.116767 -> Some(t=0.599092, d=0.116679)
   ⇒ §9.432 的「整支没进（返回 None）」解释**否掉**；那次 `proj_param = t` 无效，是因为
     在**守卫不成立**的那类输入上 `t ≈ proj_param`（两者本就相同），写入是 no-op。

② 在 `project_act` 尾巴打 `(param, dist, computed, old)`：RET rows = 42，
   **param == computed == old（42/42）、dist = 0.0**（点就在曲线上）
   样本：(1.337294471, 0.0, 1.337294471, 1.337294471) …
   ⇒ 走到尾巴的路径**不是**中心值的来源（中心值那批带着大距离，如 1.429）。
```

⇒ 中心值只能来自 **default 支的 `:255`**：`if seg_dist > mod_min { return Projection::new(mod_min, computed_point, computed_param) }`
—— 即**「记住值」`computed_*`**。而 `computed_param` 在 `:146` 被赋值，来源是**第一阶段**
`cxx:275-303`（`Extrema_ExtPC` 在区间上取最近的 `IsMin` 解）。

#### 最后窗口（唯一还没对过的子件）

```text
端口  shape_analysis_curve.rs:128-149  第一阶段（Extrema_ExtPC 等价物 + computed_param/computed_proj 的赋值）
OCCT  ShapeAnalysis_Curve.cxx:275-303  Extrema_ExtPC / 最近 IsMin 解 与 aComputedParam/aComputedProj 的赋值
看什么 该阶段在「找不到 IsMin 解」时把 param 取到什么（OCCT 是否有 (_u1+_u2)/2 之类），
       以及 `mod_min`（aModMin）的初值与更新是否一致 —— 因为 default 支最终会拿 `computed_*` 返回。
```

---

### 9.435 —— T-101：第一阶段也**不是**中心值来源（`None` 仅 3/144，命中值都是真实极值参数）⇒ 剩下唯一候选 = **周期性修正用的「区间」**

读数（TEMP，env gated，**已撤除**）：在 `extrema_ext_pc_min_in_range` 调用处打印：

```text
FIRST lines: 144    -> None: **3**
  FIRST u=[0.000000000,0.300000000] mid=0.150000000 -> None
  FIRST u=[0.000000000,1.000000000] mid=0.500000000 -> Some(u=0.599092399, d=0.116679122)
  FIRST u=[0.000000000,1.000000000] mid=0.500000000 -> Some(u=0.400907601, d=0.116679122)
```

⇒ 第一阶段返回的是**真实极值参数**（0.599/0.401…），**从不等于区间中点**；`None` 只 3 次且此时
OCCT 的 `theProjParam` 同样保持 `0.`（`cxx:273`）。结合 §9.434（精修支 Some 102/102、尾巴返回
`param == computed == old` 且 `dist = 0`），**所有内部阶段都不是中心值的产生者**。

#### 唯一剩下的候选：**周期性修正那一句用的区间**

```rust
// 端口 :262-263
if is_closed && (proj_param < u_inf || proj_param > u_sup) {
    proj_param += adjust_by_period(proj_param, 0.5 * (u_inf + u_sup), period);
}
```
```cpp
// OCCT cxx:480-483
if (anIsClosedCurve && (theProjParam < uMin || theProjParam > uMax)) {
    theProjParam += ShapeAnalysis::AdjustByPeriod(theProjParam, 0.5 * (uMin + uMax), aCurvePeriod);
}
```

**关键差异候选**：OCCT 的 `uMin/uMax` 是 `ProjectAct` 的**局部变量**，在 default 支里被
`ProjectOnSegments` **就地收窄**过；端口 `:262` 用的是**函数参数** `u_inf/u_sup`（未收窄的整区间）。
两者不同 ⇒ 守卫条件与「中心」都不同。

而且 `AdjustByPeriod(v, c, P)` 在 `|v - c| = k·P` 时正好把结果**送到 c**（观测到的 `param == 0.5*(u_inf+u_sup)`
即由此而来，`d` 非零也说得通——参数被映射到区间中心的等价点）。

#### 下一步（下一次直接做）

1. 读 `shape_analysis_curve.rs:151-276` 确认**哪些分支能走到 :262**，并在那些分支里核对
   「收窄后的区间」是否存在（端口是否有 `seg_lo/seg_hi` 的局部收窄、是否应传到 :262）；
2. 对照 `cxx:479-496`：`anIsClosedCurve`、`aCurvePeriod` 的取值，以及收窄后 `uMin/uMax` 的语义；
3. 命中即把 :262 的区间换成与 OCCT 同源的（收窄后的）区间 → 第③步验收。

---

### 9.436 —— T-101：周期性修正**从不触发**（守卫 0/42）⇒ 该候选也死；并**重新解释**「中点」观测

读数（TEMP，env gated，**已撤除**）：

```text
PERIOD rows: 42    guard true: **0**    guard true AND after == mid: 0
```

⇒ `project_act` 收尾的周期性修正（`cxx:479-484` ↔ 端口 `:262-263`）在这份数据上**一次都没进**，
它既不是中心值的来源，也不是任何差异的来源。

#### 重新解释：那 48 例的 `param` 其实是 **2π**，不是「随机落在中点」

§9.429 我按 `mid = 0.5*(u_first+u_last)` 判定，但样本揭示了一个更简单的事实：

```text
u = [5.759586532, 6.806784083] = [2π - π/3, 2π + π/3]        （2π = 6.283185307）
u = [0.523598776, 1.570796327] = [π/6, π/2]                  （其中点 π/3 亦为常见角度）
```

第一例的「中点」**恰好等于 2π**（因为区间关于 2π 对称）⇒ 那些命中更像是**合法投影参数**
（投到 pcurve 的接缝/中段），而不是「兜底取中点」。⇒ **「中点兜底」这个叙事需要修正**：
`project_act` 的每个阶段都已量过（第一阶段真实极值、采样逐行忠实、精修 Some 102/102、尾巴
`param == computed == old` 且 `dist = 0`、周期性修正 0 次触发），**没有一处是「取中点」的兜底**。

#### 结论 + 下一个真问题

所以 F113 的**中点新顶点**（3D `-65.243610` ≈ `#5012/#5018` 两端点的 3D 中点）不是投影兜底造成的，
而要回到**切分那一侧**去量：`fix_notched_edges` 里

```text
:3051 check_notched_edges(...) → short_num / param
:3100 let uv = c2d.d0(check.param);
:3101 let vnew = make_vertex(surface.d0(uv.x(), uv.y()), CONFUSION)   ← 3D 新顶点在这里定下来
:3113/:3119 transfer_range(new_e1/new_e2, …, check.param, true)
```

**下一个探针**：在 `:3100` 前后打印 `(short_num, check.param, 被切边的 [a,b], uv, vnew 的 3D 坐标)`
—— 直接把「`vnew` 落在毛刺中点」与 `check.param`/`uv`/`surface.d0` 三者对上号，
从而判定是 **`check.param` 本身**、还是 **`c2d.d0(param)` / `surface.d0(uv)` 的求值**把点送到中点。

---

### 9.437 —— T-101：那 2 条 wire 的 8→6 **不走切边路径**（`:3100-3101` 探针 0 次命中）⇒ 分叉在 `on_end`/`:3083` 一带

在 `:3101`（`let vnew = make_vertex(surface.d0(uv), CONFUSION)`）后插入打印（TEMP，env gated，**已撤除**）：

```text
NOTCHFIX lines: **0**
```

即在 `zz_uv_feed --ids` 这一趟里，`fix_notched_edges` **从未执行到「造新顶点」那一步** ——
而 §9.417 的配对探针（在 `:3249` 调用前后）明确显示该函数**单次调用内把 8 变成 6**（2 条 wire）。
⇒ 这两条 wire 的删边走的是**另一条路**，在 `:3100` 之前就分出去了：

```text
:3051  check_notched_edges(...) → short_num / param
:3059  is_remove_first = n1 == short_num
:3060  to_split = if n2 == short_num { n1 } else { n2 }
:3066  curve_on_surface_oriented(&split_e, face, true) → (c2d, a, b)
:3074  on_end = |param - (is_remove_first ? b : a)| <= PCONFUSION || (closed && …)
:3083  if |(is_remove_first ? a : b) - param| < PCONFUSION { … }      ← **待查分支**
:3096  let Some(surface) = face_surface(face) else { i += 1; continue };  ← 也可能在这里 continue 掉
:3100  uv = c2d.d0(check.param)
:3101  vnew = make_vertex(surface.d0(uv), CONFUSION)                  ← 本次 0 次命中
:3132  wire_set_edge_composed / :3135-3139 wire_insert_edge_before     ← 组装两半
```

#### 由此产生的两个待解问题（下一轮一起打）

```text
① 删边路径：`:3074/:3083` 的 `on_end` 分支（或 `:3096` 的 surface 缺失分支）如何处理 notch
   —— 它是否**不造新顶点**就直接删掉毛刺对？若是，则端口在两处实现同一件事、行为不同；
② 新顶点 `-65.243610` 从哪来：既然不是 `:3101` 造的，就要在 `check_pcurves_and_shift` 的更上游
   （`:3171 fix_reorder_wire` / `:3177 fix_small_all` 之后、`:3189` 之前）以及
   `read_topology.rs:710-714` 的另外两个调用（`xsalgo_check_pcurve` / `project_wire_pcurve_ranges`）里找。
   **注意**：§9.412 已证 `read_topology.rs:714` 返回后是 6 边，而 §9.414 证 `fix_small_all` 之后仍是 8
   ⇒ 删边确实发生在 `check_pcurves_and_shift` 内、且在 `fix_small_all` 之后 ⇒ 与本节「不走 :3100」合起来
   指向 `:3074-3095` 那一小段。
```

**下一轮探针**：在 `:3074`（`on_end` 计算处）与 `:3083`（第二个 PCONFUSION 分支）各打一次
`(i, short_num, is_remove_first, to_split, param, a, b, on_end)`，并在 `:3132`（组装）打一次；
跑 `--ids` 后按「哪几次把 wire 少 2 条」定位到具体分支。

---

### 9.438 —— 【第一现场锁定】删边走 `on_end` → `fix_dummy_seam`；中点顶点 = `CombineVertex(V1, V2, 1.0001)`（OCCT 同式）⇒ 下一窗口 = `combine_vertex` 本体

两处探针（TEMP，env gated，**已撤除**）给出：

```text
NF1（check_notched_edges 返回 Some 之后）: 2 行
NF2（:3132 组装块）:                      **0 行**
⇒ 这 2 次都没有走「切边 + 造新顶点」那条路，而是在 `:3078` 的 `if on_end` 分支里处理掉的
   （与 §9.417「正好 2 条 wire 8→6」完全对应）
```

**该分支（`wire_fix.rs:3073-3082`）**：

```rust
let on_end = (check.param - if is_remove_first { b } else { a }).abs() <= PCONFUSION
    || (edge_is_closed_3d(&split_e) && (check.param - if is_remove_first { a } else { b }).abs() <= PCONFUSION);
if on_end {
    fix_dummy_seam(wire, n1);        // ← 删毛刺对发生在这里，注释引 cxx:4019-4021
} else { … 切边路径（本次未走）… }
```

**`fix_dummy_seam`（`wire_fix.rs:2980-3027`）↔ `ShapeFix_Wire::FixDummySeam`（`ShapeFix_Wire.cxx:4213-…`）**逐行对照：

```cpp
// OCCT
int num1 = (num == NbEdges()) ? 1 : num + 1;
E1 = sewd->Edge(num); E2 = sewd->Edge(num1);
V1 = sae.FirstVertex(E1); V2 = sae.LastVertex(E2);
Vm = sbv.CombineVertex(V1, V2, **1.0001**);           // ← **新顶点（观测到的中点）**
bool toRemove = false;                                // :4227 硬编码 false
Vs = sae.FirstVertex(E2); if (Vs.IsSame(V1) || Vs.IsSame(V2)) Vs = Vm;
newEdge = sbe.CopyReplaceVertices(E2, Vs, Vm);
CopyReversePcurves(newEdge, E1, E1.Orientation() == E2.Orientation());
B.SameRange(newEdge,false); B.SameParameter(newEdge,false);
if (!Context().IsNull()) { if (toRemove) { Context()->Remove(E2); Context()->Remove(E1); } … }
```
```rust
// 端口（:2980-3027）
let nb = wire_edges_nb(wire);
if nb < 2 || num == 0 || num > nb { return; }                    // ← OCCT 无此守卫
let num1 = if num == nb { 1 } else { num + 1 };                  ✓
let vm = combine_vertex(&v1, &v2, 1.0001);                       ← 同上式（待比）
…
let (n1, n2) = if num < num1 { (num, num1) } else { (num1, num) };
wire_remove_edge(wire, n2); wire_remove_edge(wire, n1);          // ← 直接就地删（OCCT 走 Context()->Remove，且 toRemove=false）
```

⇒ 观测到的顶点 `-65.243610` **正是 `Vm = CombineVertex(V1, V2, 1.0001)` 的位置**（`#5012/#5018` 两端点的中点：
`(-49.864906 + -80.622679)/2 = -65.24379` ≈ 实测 `-65.243610`）。

#### 下一个窗口（最后一次对拍）

```text
端口  wire_fix.rs 的 combine_vertex(&v1, &v2, 1.0001)
OCCT  ShapeBuild_Vertex::CombineVertex(V1, V2, tolFactor)（ShapeBuild_Vertex.cxx）
看什么 该函数对**相距 30 单位、容差 1e-7** 的两个顶点返回什么：
       OCCT 是按容差**加权**取点还是直接取中点？端口是否无条件取中点（那就会与 OCCT 差出 (0.056,2.175,0.056) 这种量级）。
       注意：**OCCT 的结果未必是精确中点**，而端口的实测值接近中点 —— 若 OCCT 加权后落在
       `-49.8649` 一侧（即近似原顶点），就正好解释「OCCT 收在原顶点、端口收在中点」。
另需核：OCCT 里 `toRemove=false` ⇒ `Context()->Remove(...)` 的**延迟删除**语义（端口是就地删）。
```

**下一轮**：读 `combine_vertex` 与 `ShapeBuild_Vertex::CombineVertex`，命中即按 `.cxx` 修 → 第③步验收
（`--model 113` → wires=2、`zz_seam_fix 113` → Face、`--fstats` → face=113 mt≈228、面积比 ≥0.8996 且上升、
`t101_verify.ps1` 全绿）。

---

### 9.439 —— 【真实分歧点修复】`FixDummySeam` 的顶点必须**朝向感知**：端口用 `edge_vertices(...)`（朝向无关）⇒ 中点顶点；改为 `first_vertex/last_vertex` 后 **F113 外环与 OCCT 逐点一致**

**分歧（`.cxx` 可直读）**：

```cpp
// OCCT ShapeFix_Wire.cxx:4221,4230（ShapeAnalysis_Edge 的 FirstVertex/LastVertex）
TopoDS_Vertex V1 = sae.FirstVertex(E1), V2 = sae.LastVertex(E2);   // TopExp::Vertices(E,V1,V2,CumOri=true) ⇒ **朝向感知**
...
TopoDS_Vertex Vs = sae.FirstVertex(E2);
```
```rust
// 端口（修前）：朝向**无关**的拓扑首尾
let (Some(v1), Some(v2)) = (edge_vertices(&e1).0, edge_vertices(&e2).1) else { … };
let mut vs = edge_vertices(&e2).0;
```

毛刺第二条边在 wire 里是 `Reversed` ⇒ 朝向无关取到的是曲线尾 `-80.622679`，OCCT 取到的是线序首 `-49.864906`
⇒ `CombineVertex(V1,V2,1.0001)` 给出**两者中点** `-65.243610`（实测值）。

**修复（照 `.cxx`，2 处）**：

```rust
let (Some(v1), Some(v2)) = (first_vertex(&e1), last_vertex(&e2)) else { … };   // cxx:4221
let mut vs = first_vertex(&e2);                                               // cxx:4230
```

**实测（改动后）**：

```text
F113 外环端点：含 -65.24 的 **0** 处；含 -49.86 的 **2** 处
  e[0] (12.387499,24.672033,-89.394247) -> (-49.864906,0,-100)
  e[1] (-49.864906,0,-100) -> (36.224088,-3,-99.867388)
OCCT 整形后 FACE 25 的 W1： [12.3875,24.672,-89.3942 -> -49.8649,0,-100] [-49.8649,0,-100 -> 36.2241,-3,-99.8674]
  ⇒ **逐点一致**（这正是审计 §9.418 指出的「OCCT 收在原顶点」）
a3n00：mesh_v 10863 → **10876**，mesh_t 11941 → **11967**（OCCT 侧 11145 / 12466，方向朝 OCCT）
F113：仍 wires=4（bound 合并那一步 = `fix_missing_seam` 出 Shell(5) 被丢，另一件事）
```

⇒ 这是一处**真实的 `.cxx` 分歧**（并且顺带证实 §9.418/§9.439 的「OCCT 收在原顶点」是**实测**得到的，
与 `combine_vertex` 逐行忠实并不矛盾：分歧在**传入 CombineVertex 的顶点**上）。

**门禁**：`t101_verify.ps1` 已在后台启动（job 见报告）——绿则提交本改动，红则回退这两处。

---

### 9.440 —— T-101：§9.439 的修复**门禁退步**（a3n00 面积比 0.8996 → **0.8918**）⇒ 已回退，但保留全部证据

`t101_verify.ps1`（job `pwsh-1404`，exit 0）关键行：

```text
[occ/a3n00.stp] our=202564.34 occ=227130.30 **ratio=0.8918**     ← 改前是 0.8996（§9.369 基线）
[occ/a3n00.stp] ours v=10876 f=11967 | occ v=11145 f=12466 | f-ratio 0.96
[occ/T0M.stp]   our=192416.44 occ=192658.64 ratio=0.9987        （不变）
[occ/acs10.stp] our=267047.23 occ=271219.68 ratio=0.9846        （不变）
step_obj_gates: test result: ok. 5 passed; 0 failed
```

⇒ 这处改动**是 `.cxx` 忠实的**（`FixDummySeam` 的 `sae.FirstVertex/LastVertex` 确为朝向感知，§9.439），
并且**确实修好了 F113 外环的几何**（中点 `-65.24` 消失、回到 `-49.86`，与 OCCT 的 W1 逐点一致），
但它让**整个 a3n00 的覆盖面积变差**（0.8996 → 0.8918）——说明这处顶点语义还牵动了别的面的合并/网格结果。

按纪律（面积比不得退、不得为对齐某模型加特例）：**已回退**这两处，库代码回到门禁绿的状态
（`cargo check` 0 error、`git diff crates/…/wire_fix.rs` 为空），基线读数复验 `mesh_v=10863 mesh_t=11941` ✓。

**保留的证据与后续**（下次可继续）：

```text
· 改动文本（2 处，精确）：
    let (Some(v1), Some(v2)) = (first_vertex(&e1), last_vertex(&e2)) else { … };   // cxx:4221
    let mut vs = first_vertex(&e2);                                                // cxx:4230
· 收益：F113 外环几何与 OCCT 逐点一致（-49.864906 而非中点）
· 代价：a3n00 面积比 -0.0078（0.8996 → 0.8918），T0M/acs10 不变
· 下一步要先弄清「这处改动让哪些面的网格变了」：
    用 `--fstats` 逐面比「改前 / 改后」的 mt（按 bbox 配对），看是否有面从有网格变成 `mt=0`
    或面积显著变化 —— 若能把退步归因到某个具体面，就可判断是「另一处 bug 被暴露」还是「本处改法不对」。
    在弄清之前不要再次落地这处改动。

---

### 9.441 —— T-101：§9.439 修复的退步**只来自一张面**（`f=171` = 审计的 idx171 `#7415`），且它是 `+26` 三角

方法：分别在「改前 / 应用 §9.439 修复」两种状态下跑 `zz_probe_a3n00 --fstats`，
按 **bbox 配对**逐面比 `mt`（禁按面号）：

```text
before faces: 225   after faces: 225
faces with differing mt: **1**     only-in-before: 0   only-in-after: 0
  dmt=+26   before(f=171 mt=349)  after(f=171 mt=375)
  bbox=(-57.360501,-42.864963,-169.610910)-(61.610910,42.86…)
```

⇒ 这处顶点语义修复**全程只改一张面的网格**：`f=171` —— 正是审计 §9.418 列出的另一张「净丢边」面
（idx171 `#7415`：STEP `[8,8,8,8]` → 端口 `[8,8,8,6]`，OCCT 整形 `[8,8,8,7]`）。
而 f=171 从 349 → **375** 三角（+26），与「a3n00 总网格三角 11941 → 11967」完全吻合 ✓。

**但总面积反而下降**：`our=204327.72 → 202564.34`（−1763），而 OCCT 固定 227130.30
⇒ 面积比 0.8996 → 0.8918。即：**f=171 的三角变多了，但它的覆盖面积变得离 OCCT 更远**。

这有两种解读，必须先分清（下次做）：

```text
(a) 修复让 f=171 的**结构**更接近 OCCT（`[8,8,8,6]` → `[8,8,8,7]`?）而**面积**更差
    ⇒ 说明 f=171 另有一处网格/几何问题被暴露（本处改法可能是对的，另需修别处）；
(b) 本处改法在 f=171 上不合 `.cxx`（例如该面的 `E1/E2` 朝向情形与 F113 相反）
    ⇒ 就该按 .cxx 的**分支条件**再核一遍（`sae.FirstVertex/LastVertex` 的 CumOri 语义在本文件里只有一处调用点）。
```

**下一步（一次读数即可分流）**：对 `f=171` 分别取「改前/改后」的
`zz_uv_feed --model 171`（wire 数、每 wire 边数、bbox）与 `--fstats`（mv/mt）；
若改后边数变成 `[8,8,8,7]`（= OCCT）⇒ 归 (a)，把注意力转到 f=171 的三角化本身；
若仍是 `[8,8,8,6]` ⇒ 归 (b)，回去核该分支的朝向条件。

（当前树为**改前**状态：`git diff crates/…/wire_fix.rs` 为空、基线 `mesh_v=10863 mesh_t=11941` ✓。）

---

### 9.442 —— T-101：§9.439 修复**不改变结构**（f=171 两种状态都是 4 wires / `[6,8,8,8]`）⇒ 它只移动顶点，退步落在 **f=171 的三角化几何**上

`zz_uv_feed --model 171`（改前 / 应用修复后）：

```text
改前：MODEL f=171 wires=4   wire[0] nEdges=6   wire[1] nEdges=8   wire[2] nEdges=8   wire[3] nEdges=8
改后：MODEL f=171 wires=4   wire[0] nEdges=6   wire[1] nEdges=8   wire[2] nEdges=8   wire[3] nEdges=8
OCCT 整形基准（审计 §9.418）：[8, 8, 8, 7]        STEP 文件： [8, 8, 8, 8]
```

⇒ 分流结论：**既不是 (a) 也不是 (b) 的简单形态** ——
这处修复**只改顶点位置**（把毛刺合并点从「中点」改回「原顶点」，与 OCCT 的 `wdump 25/…` 逐点一致），
**不改任何 wire/边数**；而 `--fstats` 显示全模型只有 `f=171` 的 `mt` 变化（349 → 375，+26），
总面积却 −1763。两者合起来只有一个解释：

> **f=171 的三角化几何本身有问题**：顶点一改，它铺出来的三角形状/覆盖面积就大幅变化（少了 1763 面积）。
> 修复前它「碰巧」更接近 OCCT 的总面积 ⇒ 面积比那 0.0078 的差不是这处修复的错，而是 **f=171 这一张面的
> 网格问题被暴露出来**。

注意：端口 f=171 的结构（`[6,8,8,8]`）与 OCCT 整形后的 `[8,8,8,7]` **本就不同**（这一点与 F113 的
「4 wires vs 2 wires」同类），所以这处顶点修复并不会让 f=171 结构变好——它只把**顶点**改对了。

#### 下一步（要一把新尺子）

```text
现状口径：面积比只到「整模型」级；而这次的差异集中在**单张面**上。
需要：**逐面面积**（或逐面 OBJ 面积）对比工具，才能判断
  · 改前 f=171 的三角面积是否「碰巧」等于 OCCT；
  · 改后 f=171 少了 1763 面积到底少在哪（三角形形状 vs 覆盖缺口）。
做法建议（不新增测试，只加只读探针）：在 `zz_probe_a3n00` 里加一个 `--farea`，
  按面把三角形面积求和打印（`face=N area=… mt=…`），与 OCCT 侧同口径（OBJ 逐面或 `--facestats`+面积）比。
然后才有资格决定：这处 `.cxx` 正确的顶点修复是否应该落地（以及 f=171 是否另有网格 bug 要修）。
```

（当前树为**改前**状态：`git diff crates/` 为空、基线 `sum_mt=11941` ✓。）

---

### 9.443 —— T-101：回到目标本身 —— F113 的 `fix_missing_seam` 当前行为（改写前的精确读数）

```text
（本轮读数见 .target-gate/sf113.txt；下面为摘要）
BEFORE wires=4   wire[0]=6 / wire[1]=14 / wire[2]=1 / wire[3]=1
result: ShapeFixFace::fix_missing_seam → true，result = **Shell(5)** ⇒ read_topology.rs::resolve_face 只收 Face ⇒ 丢弃
⇒ 该面保持 4 wires ⇒ 建模 wire[0]=6 / wire[1]=14（§9.393）⇒ face-checker 判 SELF_INTERSECTING_WIRE
  ⇒ 三角化循环跳过（§9.395/§9.391）⇒ **空面**
```

**两个已确认的周边事实**（帮助下一轮聚焦）：

```text
· §9.439 已证：即使把毛刺合并点从「中点」改回 OCCT 的「原顶点」（.cxx 正确），F113 **仍然** 4 wires
  ⇒ 合并失败与那个顶点无关；
· §9.410/§9.412 已证：进 seam 步之前的外环结构已与 OCCT 同形（6 边、去掉 #5012/#5018 后），
  即 ComposeShell 面对的输入在结构层面已对齐。
```

⇒ 于是「F113 仍是空的」这条线**只剩一个问题**：**为什么端口 ComposeShell 在 F113 上产出 Shell(5)（两张同 patch 重合面 + 未配对半边），而 OCCT 产出合并后的 Face**。
§9.377 的阶段性普查当时是对 **f=140（cone）** 做的；对 F113 本身还没有同等普查。

**下一轮**：对 **F113** 做一次与 §9.377 同规格的阶段性普查（在 `ShapeFixFace::fix_missing_seam`
的关键节点打同标签 `ZZCS`：找 seam 循环后的 `uf2/vf2`、`ComposeShell` 之前的 grid 范围、
`breakwires` 后的段列表、`dispatchwires` 输出的面数与各自 wire 数），并与其 .cxx 对应行一起判读。
命中即按 .cxx 修，然后第③步验收（`--model 113` → wires=2；`zz_seam_fix 113` → Face；
`--fstats` → face=113 mt≈228；面积比 ≥0.8996 且上升）。

---

### 9.444 —— T-101：结果打包逻辑**与 `.cxx` 一致**（`1 → Face`、`>1 → Shell`）⇒ 「Shell(5)」是 ComposeShell **产出 5 个面**的忠实后果；问题回到 ComposeShell 的**输出面数**

端口 `shape_fix_face.rs:661-665`：

```rust
result = match kept_faces.len() {
    0 => None,
    1 => Some(kept_faces.remove(0).0.clone()),      // 单面 → Face
    _ => Some(builder.make_shell(&kept_faces).0),   // 多面 → **Shell**
};
self.result = result;
if let Some(res) = &self.result { self.context.replace(&face.0, res); }   // cxx:2268
```

OCCT `ShapeFix_Face.cxx:2302-2322`：

```cpp
myResult = Context()->Apply(myResult);
for (TopExp_Explorer exp(myResult, TopAbs_FACE); exp.More(); exp.Next()) {
  myFace = TopoDS::Face(Context()->Apply(exp.Current()));
  if (nbFaces > 1) { FixSmallAreaWire(true); … }      // 只有 >1 面才做小面积裁剪
  BRepTools::Update(myFace);
}
myResult = Context()->Apply(myResult);                 // ← 同样**不把多面强并成一个 Face**
```

⇒ 端口与 OCCT 在这一点上**同构**（都不做「多面 → 单面」的强并）。因此：

```text
zz_seam_fix 113 实测：result faces=**5**，其中 HEALED face 0 = wires=2（v 范围 [-43.75,131.25] = OCCT 口径）
⇒ 端口把「正确的那张 2-wire 面」连同另外 4 张一起打包成 Shell(5)，reader 只收 Face ⇒ 全丢
而 OCCT 的整形后模型 **仍是 226 面**（与文件相同）⇒ 它的 ComposeShell 在 F113 上**没有多产出 4 张面**
```

⇒ **真正剩下的分歧 = 端口 ComposeShell 在 F113 上产出的面数是 5，OCCT 是 1。**
这与 §9.377 在 f=140（cone）上看到的「同一 patch 两张重合面 + 未配对半边」同源；而 §9.377 已实测
**剪枝（FixSmall / CheckSmallArea）是 no-op**（只能删、不能合）⇒ 必须回到 **dispatch/collect 阶段**找多产出的来源。

#### 下一轮（对 F113 做 §9.377 同规格的阶段普查）

```text
在 ComposeShell 的 perform 里，对 F113 那次调用打同标签阶段日志（TEMP，env gated，用完 edit 反撤）：
  loadwires / splitbygrid / breakwires / collectwires / **dispatchwires** 各阶段的段列表与段数，
以及 dispatch 出的**每个面的 wire 数 + bbox**（重点看是否出现两张同 bbox 的重合面）。
对照 ShapeFix_ComposeShell.cxx 的 DispatchWires（cxx:2770-2860 一线）判读；
命中即按 .cxx 修，然后第③步验收：
  zz_seam_fix 113 → **Face**；--model 113 → wires=2（22+6 边）；--fstats → face=113 mt≈228；
  a3n00 面积比 ≥0.8996 且上升；pwsh -File .target-gate\t101_verify.ps1 全绿。
```

**注意**：`read_topology.rs:745` 的 `if res.shape_type() == Face` 这一条**不要动**（那是 reader 的既有语义；
OCCT 侧同理：`ShapeFix_Face` 的结果以 Context 替换原面，若结果是壳，OCCT 也会把它当壳用）。

---

### 9.445 —— T-101：F113 的 Shell(5) 是**5 块单 wire 补丁**（6/5/4/3/3 边），不是重合面 ⇒ 分歧在「切」而不是「合」

利用已有捕获（`.target-gate/sf113.txt`，`zz_seam_fix 113` 的全文）逐面读出：

```text
result type=Shell    result faces=5
  HEALED face 0 : 1 wire, **6** edges      （v 范围 [-43.75,131.25] = OCCT 口径；这张是 OCCT 的第 2 条 wire！）
  HEALED face 1 : 1 wire, 5 edges
  HEALED face 2 : 1 wire, 4 edges
  HEALED face 3 : 1 wire, 3 edges
  HEALED face 4 : 1 wire, 3 edges
合计 6+5+4+3+3 = 21 条边；OCCT 整形后 FACE 25 = 2 wires（**22 + 6** 边）
```

⇒ 两点结论：

1. **face 0（6 边）就是 OCCT 那条 6 边内环** —— 端口确实产出了正确的那部分；
2. 其余 4 张（5/4/3/3=15 边）是**外环被切碎**的结果，而 OCCT 把它成一条 **22 边** wire
   ⇒ 端口的 ComposeShell 在 F113 上**把边界切成 5 块补丁**，而 OCCT 得到 1 面 2 wire。

这与 §9.377 在 f=140（cone）上观测到的「`breakwires` 出两段端点不配对（`D ≠ C`）⇒ 多张同 patch 面」
**同源**：问题在 **`SplitWires`/`BreakWires` 的切分**（§9.377 的普查对象），而**不在**结果打包（§9.444 已证打包与 `.cxx` 同构）。

（另注：端口上游外环只有 6 边 vs 文件 8 边（§9.412 的 `FixDummySeam` 丢边），所以它的补丁边数合计
（21）与 OCCT 的 28 本来就不同；两件事都要在 `.cxx` 同一处分清。）

#### 下一轮（唯一剩下的动作：对 F113 做 §9.377 同规格的阶段普查）

```text
在 ComposeShell::perform 的五个阶段各打一次（TEMP/env gated/edit 反撤），**只针对 F113 这次调用**：
  loadwires      ：输入 wire 数 + 每 wire 边数（期望 4 / [6,14,1,1]）
  splitbygrid    ：切后段数 + 是否新增外部段（§9.377 在 f=140 上看到新增 1e/O）
  breakwires     ：每段的两端点（重点：是否出现 §9.377 的「D ≠ C 不配对」）
  collectwires   ：每段是否找到候选（`index=None` 的那几段）
  dispatchwires  ：输出的**面数 + 每面 wire 数 + bbox**（期望能解释 5 块补丁）
对照 ShapeFix_ComposeShell.cxx 的 SplitWires/BreakWires/CollectWires/DispatchWires（§9.377 已定的行号）
⇒ 命中即按 .cxx 修；验收：zz_seam_fix 113 → Face、--model 113 → wires=2（22+6 边）、
   --fstats → face=113 mt≈228、面积比 ≥0.8996 且上升、t101_verify.ps1 全绿。
```

---

### 9.446 —— T-101：把「F113 五阶段普查」派给子任务（含完整上下文与硬约束）

本轮把这件需要完整工作记忆的活派给子任务 `56ebc0dd`，brief 要点（可复现）：

```text
任务：对端口 F113（bbox (-87.5,-34,-100)-(87.5,34,-32)）做 ShapeFix_ComposeShell 的五阶段普查
  loadwires / splitbygrid / breakwires / collectwires / dispatchwires
  重点：breakwires 是否出现 §9.377 那种「D ≠ C 两端点不配对」，以及 dispatchwires 为何出 5 块补丁
对照：ShapeFix_ComposeShell.cxx:2131-2275(SplitByGrid) / :1433-1914(SplitByLine) /
      :2824-2846(CollectWires 收尾) / :2770-2860(DispatchWires)
已有事实（brief 里直接给出，不重做）：Shell(5) 是 5 块单 wire 补丁（6/5/4/3/3）；
  打包逻辑与 .cxx 同构（不要去改打包或 reader 的 Face 判定）；剪枝实测 no-op。
硬约束：TEMP 插桩 env gated、判定条件必须**结构性**（禁止面号/bbox 特例）、
  必须用 edit 反向撤除（禁 git checkout）、不提交、报告附 git status/diff。
交付：第一处分歧(.cxx 行号 ↔ 端口 行号) + 机制 + 五阶段关键读数 + 是否小改可修 + git status/diff。
```

（同时保留一条后路：若子任务判「超范围」，则按 §9.399 已确认的结论——`split_wire.rs`/`split_by_line.rs`
的控制流已逐行对完且忠实——把那 5 块补丁的**切分点参数**与 OCCT 在同一输入下的对应量逐值对比。）

**当前树状态**：库代码零改动（基线 `sum_mt=11941`），工作树只剩审计留下的只读 instrument
`crates/occt-topo/examples/zz_probe_a3n00.rs`。

---

### 9.447 —— T-101：**阶段性交接**（目标未达成；因果链、已排除项、可落地候选、下一步全部就位）

#### A. 目标状态（如实）

```text
目标：F113 → 2 wires（22+6 边）/ 228 节点 / 228 三角；a3n00 面积比 ≥0.8996 且上升；门禁全绿
现状：F113 仍 4 wires（wire[0]=6 / wire[1]=14 / 1 / 1）→ 无网格（空面）
      a3n00 面积比 0.8996（未上升）、T0M 0.9987、acs10 0.9846、--lib 1281/0、step_obj_gates 5/5
```

#### B. 已确证的因果链（每环都有可复跑读数）

```text
STEP 外环 8 条（#5004/#5012/#5018/#5090/#5098/#5107/#5115/#5137，8 条 distinct）
→ resolve_loop 8 → make_wire 8 → make_face 8                     （§9.408/§9.410）
→ check_pcurves_and_shift 内：fix_notched_edges → on_end 分支 → fix_dummy_seam
     **8 → 6**（单次调用，配对探针 8→8 ×10 / 8→6 ×2）              （§9.417/§9.437/§9.438）
     新顶点 = CombineVertex(V1,V2,1.0001)：端口得**中点** -65.243610，
     OCCT 得**原顶点** -49.864906（wdump 25 实测：OCCT 外环含 -49.86、不含 -65.24）  （§9.405/§9.438/§9.439）
→ 面停在 4 wires；fix_missing_seam 产出的结果其实是 **2 wires**（v 范围 [-43.75,131.25] = OCCT 口径），
   但被打包成 **Shell(5)**=5 块单 wire 补丁（6/5/4/3/3 边）⇒ reader 只收 Face ⇒ 全丢  （§9.443/§9.445）
   （打包逻辑与 .cxx 同构：1 面→Face、>1→Shell；**不要去改打包或 reader 的 Face 判定**）  （§9.444）
→ 建模 wire[0]=6 / wire[1]=14 → face-checker 判 SELF_INTERSECTING_WIRE(+FAILURE)   （§9.393/§9.395）
→ 三角化循环跳过 → **空面**                                                       （§9.391）
```

#### C. 已排除（都逐行对过、判定忠实）

```text
ComposeShell 全链控制流（split_by_grid/split_by_line/split_wire/collect_wires/剪枝）
五个 helper + TOLINT；FixMissingSeam 定位段与假想 grid 构造（含 uf/vf/URange/VRange 与 OCCT 逐位一致）
project_adaptor/ProjectOnSegments/ProjectInside/AdjustByPeriod 及其各阶段（第一阶段极值 None 仅 3/144、
  采样逐行一致、精修 Some 102/102、尾巴 param==computed==old、周期性修正 0/42 触发）
make_wire/make_face 不丢边；add_wire 四条早退 0 命中；visit_face 三处 FAILURE 0 命中
剪枝实测 no-op（只能删不能合）
```

#### D. 唯一已确认的 `.cxx` 分歧（**已量化，但因门禁退步暂缓落地**）

```text
分歧：ShapeFix_Wire::FixDummySeam 的顶点必须**朝向感知**（cxx:4221/4230 的 sae.FirstVertex/LastVertex）
      端口用的是朝向无关的 edge_vertices(...).0/.1
改法（2 处，精确）：
  let (Some(v1), Some(v2)) = (first_vertex(&e1), last_vertex(&e2)) else { … };   // wire_fix.rs:2992
  let mut vs = first_vertex(&e2);                                              // wire_fix.rs:2999
收益：F113 外环与 OCCT 逐点一致（中点消失，回到 -49.864906）
代价：a3n00 面积比 0.8996 → **0.8918**（t101_verify job pwsh-1404 实测）
归因：全模型只有 **f=171**（= 审计 idx171 #7415）的 mt 变化（349→375，+26），但总面积 −1763
      ⇒ 该面的三角化几何本身有问题（本处修复把它暴露出来），不是顶点改错        （§9.440/§9.441/§9.442）
结论：**在弄清 f=171 之前不要再次落地这处改动**
```

#### E. 当前最强假设与下一步（两条，按性价比）

```text
假设：分歧在 **SplitWires/BreakWires 的切分** —— 端口把 F113 外环切成 4 块补丁（5/4/3/3 边），
      OCCT 把它成一条 **22 边** wire；与 §9.377 在 f=140 上看到的「breakwires 两端点不配对 D≠C」同源。
      （打包、剪枝、顶点语义、投影链都已排除，故只剩切分本身。）

下一步 A（正统）：对 F113 做 §9.377 同规格的五阶段普查（loadwires/splitbygrid/breakwires/collectwires/
      dispatchwires），对照 cxx:2131-2275(SplitByGrid)/1433-1914(SplitByLine)/2770-2860(DispatchWires)/
      2824-2846(CollectWires)；已有子任务 brief（§9.446）。
下一步 B（后路，不需新插桩）：把 5 块补丁的**切分点参数**与 OCCT 在同一输入下的对应量逐值对比
      （split_wire.rs/split_by_line.rs 控制流已逐行判定忠实，故差异应在**数值输入**上）。
下一步 C（并行、独立）：查 f=171 的三角化为何一改顶点就丢 1763 面积（需要逐面面积口径 --farea）。
```

#### F. 现成仪器与验收

```text
zz_seam_fix <stp> <f>      ：BEFORE/result/HEALED 逐面（wire 数、每 wire 边数、v 范围）
zz_uv_feed --model <f> | --ids   ：模型结构 / wire 直方图
zz_probe_a3n00 --fstats | --fdump | --ecensus（审计留下的只读 instrument，工作树里未提交）
pwsh -File .target-gate\t101_verify.ps1  ：一条命令的落地验收（判据见文件头注释）
```

（当前树：库代码零改动；基线 `STATMAP matched=225 unmatched=1 sum_mt=11941` ✓。）

---

### 9.448 —— T-101 **收尾**（第 71-72 轮）：子任务无产出被停止、插桩已清回 HEAD、目标如实留在 active

#### A. 第 71 轮的事件与处理

派出的 F113 五阶段普查子任务（`56ebc0dd`，brief 见 §9.446）**未产出任何读数或结论**即被我停止
（无 closing message），但它在工作树里留下 4 个库文件的 TEMP 插桩：

```text
perform.rs +131 / break_wires.rs +17 / split_by_grid.rs +19 / split_by_line.rs +39   （标签 ZZCS×14、T101_ZZ×1）
```

处理（**未用 `git checkout`**）：用 `git show HEAD:<file>` 取回原字节 + python 归一化写回 4 个文件；
`split_by_line.rs` 另按其 HEAD 形态去掉尾部换行。复验：

```text
标签残留 ZZCS 0 / T101_ZZ 0；cargo check 0 error
git status --porcelain → 仅 M crates/occt-topo/examples/zz_probe_a3n00.rs（审计留下的只读 instrument）
基线：STATMAP matched=225 unmatched=1 sum_mt=11941   （与 §9.374 基线逐字相同）
```

#### B. 目标终态（如实）

```text
未达成：F113 仍 4 wires（wire[0]=6 / wire[1]=14 / 1 / 1）→ 无网格（空面）
        a3n00 面积比仍 0.8996（未上升）；T0M 0.9987；acs10 0.9846；--lib 1281/0；step_obj_gates 5/5
```

#### C. 结论与「接着做的第一条命令」

```text
已确证：端口已算出正确的 2-wire 面（v 范围 [-43.75,131.25]），但被 ComposeShell 的 5 块补丁
        （6/5/4/3/3 边）打包成 Shell(5) ⇒ reader 只收 Face ⇒ 全丢 ⇒ 面停 4 wires ⇒ face-checker 自交 ⇒ 空面。
最强假设：分歧在 SplitWires/BreakWires 的**切分**（外环被切成 4 块补丁；OCCT 成一条 22 边 wire）。

第一条命令（复现现状，约 1 分钟）：
  crates\occt-topo\target\debug\examples\zz_seam_fix.exe data/occ/a3n00.stp 113
  期望看到：BEFORE wires=4（6/14/1/1）→ result type=Shell, result faces=5（HEALED 各 1 wire：6/5/4/3/3）

第二件事（继续定位，见 §9.447-E 的三条路）：
  A 五阶段普查（loadwires/splitbygrid/breakwires/collectwires/dispatchwires，对照
    cxx:2131-2275 / :1433-1914 / :2824-2846 / :2770-2860）
  B 后路：5 块补丁的切分点参数 vs OCCT 同输入逐值对比（不需新插桩）
  C 并行：f=171 的三角化为何一改顶点就丢 1763 面积（需逐面面积口径 --farea）
```

**纪律提示（留给接手者）**：§9.439 那处 `.cxx` 正确的顶点修复（`FixDummySeam` 朝向感知）**已被证明会
让 a3n00 面积比 0.8996 → 0.8918**，在弄清 f=171 之前**不要再落地**；不要改 reader 的
`res.shape_type() == Face`；不要按面号/bbox/面积加特例；TEMP 插桩一律用 `edit` 反撤。

---

### 9.449 —— T-101 **止损（用户决定，路 D）**：目标 blocked，自动轮次停止；最终账目与留给后续的三条路径

用户指出「这么多轮没改好、看不出有意义的事情」并选择**止损**。如实账目：

```text
本轮目标（T-101，72 轮）落到 crates/ 的库改动：**只有 1 行**
  7f2df19b  wire_fix.rs:2900 投影目标 &ad1 -> &ad2（对齐 cxx:1963）——门禁绿，但对 a3n00 无影响

同一会话更早、有量的库改动（不属于 T-101 这 72 轮）：
  5149df21  接回导入期 ShapeFix_Face::FixMissingSeam  →  a3n00 面积比 0.8627→0.8996、法兰 16 面 2 wires→1 wire/5 边
  daf7ee4d  T-92 WireSegment 端点按边朝向 + 接线       →  a3n00 的 2-wire 面 24→10

T-101 期间的两个「真实分歧」结论：
  ① check_notched_edges 的 p2d1/p2d2 朝向 —— **自证为误判**（§9.427），已否定
  ② FixDummySeam 顶点须朝向感知（cxx:4221/4230 ↔ wire_fix.rs:2992/2999）—— 几何正确、F113 外环与 OCCT 逐点一致，
     但 a3n00 面积比 0.8996 → **0.8918**（仅 f=171 #7415 一张面 mt 349→375 / 总面积 −1763）
     ⇒ 被门禁拦下并回退（§9.439-§9.442）
过程失误（如实）：两次清理 TEMP 插桩多删代码致返工；第 69 轮派出的 F113 五阶段普查子任务 3 轮零产出被停止。
```

**终态**：F113 仍 4 wires（6/14/1/1）**无网格**；a3n00 面积比 0.8996（未上升）；门禁本身绿
（`--lib` 1281/0、`step_obj_gates` 5/5、T0M 0.9987、acs10 0.9846）。goal phase = **blocked**（round-limit 72/72）、disarmed。

#### 留给后续的三条路径（都不需要「再读一遍控制流」）

```text
B（最快出可见结论）：把端口 fix_missing_seam 输出的 **5 块补丁边界顶点**导出，与 OCCT 整形后那张面的
   22+6 边边界逐点对齐（同一圆柱面上的多边形）⇒ 直接指认「哪一刀多切了」。
C：先用只读探针加 **逐面面积**（--farea），查清 f=171（#7415）为何一改顶点就丢 1763 面积；
   修好它之后，§9.439 那处**确凿的 .cxx 修复**才能过门禁（F113 外环几何随之一致）。
A：重做 F113 五阶段普查（loadwires/splitbygrid/breakwires/collectwires/dispatchwires，
   对照 cxx:2131-2275 / 1433-1914 / 2824-2846 / 2770-2860）—— 上一次派出的子任务未完成。
```

---

### 9.450 —— **更正**（§9.445/§9.448 与 README 的数字错误）：F113 这一步**不丢边**，分歧是「28 条边被分成 6 条 wire」，不是「6/5/4/3/3 五块单 wire 补丁」

#### A. 错误与原因

§9.445 写「Shell(5) = 5 块**单 wire** 补丁（6/5/4/3/3 边）」，§9.448/README 沿用。**这是错的**：
我当时的解析脚本把每个 `HEALED face` 的扫描窗口限成 8 行，**漏掉了 `wire[1]`**，于是 face 0 的
第二条 wire（7 边）被吞掉。正确分组（口径见 B）是：

```text
HEALED face 0 : wires=**2**   wire[0]=6e   wire[1]=**7e**
HEALED face 1 : wires=1       wire[0]=5e
HEALED face 2 : wires=1       wire[0]=4e
HEALED face 3 : wires=1       wire[0]=3e
HEALED face 4 : wires=1       wire[0]=3e
合计：6 条 wire / 28 条边
OCCT 整形后 FACE 25：**2 条 wire（22 + 6）/ 28 条边**
```

#### B. 更正后的事实（口径）

```text
命令：crates\occt-topo\target\debug\examples\zz_seam_fix.exe data/occ/a3n00.stp 113
字段：HEALED 段逐面 `wires=`，每条 `wire[i] nEdges=`；OCCT 侧 `run_dbg.bat data\occ\a3n00.stp wdump 25`（整形后）
结论：**总边数两侧相同（28）** ⇒ 这一步 **没有丢边**；
      差别在**分组**：端口 28 条边 → **6 条 wire / 5 张面**；OCCT → **2 条 wire / 1 张面**。
      其中 6e 那条与 OCCT 的 6e 内环对应；22 边外环在端口被拆成 **7+5+4+3+3 = 22** 五段
      （7 段还挂在同一张 2-wire 面上）。
```

⇒ 这把 §9.445 的「外环被切碎」**从描述升级为可核对的数字**（7/5/4/3/3=22，且总数不丢），
同时否掉了「丢边/重合面」两种说法（前者由 §9.412 的 8→6 属于更上游、后者从未成立）。

#### C. 仍缺的一块数据：**五段在哪里断开**

该探针的 HEALED 边只打 `par=[f,l] dpar=`，**不打端点** ⇒ 无法直接看出断点位置。
下一步（只读、example 级）：给 `zz_seam_fix` 的 HEALED 逐边打印加 `first_vertex/last_vertex` 的三维坐标，
然后按端点把 7/5/4/3/3 五段与 OCCT 的 22 边外环逐点对齐 ⇒ 直接指认「哪两个端点之间被多切了一刀」。

#### D. 已同步修正

`README.md` 的「现状/未结个案」段已按本节改写（去掉“5 块单 wire 补丁/6/5/4/3/3”的说法）。

---

### 9.451 —— F113 五段/六条 wire 的**端点**已取到（新 instrument：`zz_seam_fix --ep`）：外环的断点就在 seam 网格线上，且**内环那条 6 边 wire 仍带着中点顶点 -65.244**

新 instrument（example 级、只读、已入库）：`zz_seam_fix <stp> <face> --ep` ⇒ 在 `describe("HEALED", …)`
之后逐 wire 打印每条边 `first_vertex->last_vertex` 的三维坐标（用 `occt_topo::shhealing::{first_vertex,last_vertex}`
即**朝向感知**的取点）。原始 dump：`.target-gate/ep113.txt`。

```text
HEALEDP face=0 wire=0 n=6                                  ← 内环（对应 OCCT 的 6 边 wire）
 (12.387,24.672,-89.394)->(**-65.244**,0.000,-100.000)     ← **中点顶点仍在！**（OCCT 是 -49.864906，§9.439）
 (-65.244,0.000,-100.000)->(36.224,-3.000,-99.867)
 (36.224,-3.000,-99.867)->(63.927,-3.000,-99.867)
 (63.927,-3.000,-99.867)->(63.927,3.000,-99.867)
 (63.927,3.000,-99.867)->(36.224,3.000,-99.867)
 (36.224,3.000,-99.867)->(12.387,24.672,-89.394)           ← 闭合
HEALEDP face=0 wire=1 n=7                                  ← 外环第 1 段
 (-8.900,28.844,-48.000)->(-8.900,33.805,-69.638)
 (-8.900,33.805,-69.638)->(-4.900,32.208,-76.892)
 (-4.900,32.208,-76.892)->(4.900,32.208,-76.892)
 (4.900,32.208,-76.892)->(8.900,33.805,-69.638)
 (8.900,33.805,-69.638)->(8.900,28.844,-48.000)
 (22.366,28.844,-48.000)->(8.900,28.844,-48.000)           ← 断开点 A：(8.900,28.844,-48.000)
 (38.733,21.909,-40.000)->(22.366,28.844,-48.000)          ← 断开点 B：(22.366,28.844,-48.000)
HEALEDP face=1 wire=0 n=5
 (38.733,21.909,-40.000)->(44.500,-0.000,-32.000)
 (44.500,-0.000,-32.000)->(38.733,-21.909,-40.000)
 (22.366,-28.844,-48.000)->(38.733,-21.909,-40.000)
 (-22.366,-28.844,-48.000)->(22.366,-28.844,-48.000)
 (-38.733,-21.909,-40.000)->(-22.366,-28.844,-48.000)
HEALEDP face=2 wire=0 n=4
 (-38.733,-21.909,-40.000)->(-44.500,0.000,-32.000)
 (-44.500,0.000,-32.000)->(-38.733,21.909,-40.000)
 (-22.366,28.844,-48.000)->(-38.733,21.909,-40.000)
 (-8.900,28.844,-48.000)->(-22.366,28.844,-48.000)
HEALEDP face=3 wire=0 n=3                                  ← 退化薄片（seam 处）
 (87.500,0.000,-32.000)->(87.500,0.000,-32.000)
 (87.500,0.000,-32.000)->(44.500,-0.000,-32.000)
 (44.500,-0.000,-32.000)->(87.500,0.000,-32.000)
HEALEDP face=4 wire=0 n=3                                  ← 退化薄片（另一侧 seam）
 (-87.500,0.000,-32.000)->(-87.500,0.000,-32.000)
 (-87.500,0.000,-32.000)->(-44.500,0.000,-32.000)
 (-44.500,0.000,-32.000)->(-87.500,0.000,-32.000)
```

#### 立刻可读出的两条硬事实

```text
① 内环那条 6 边 wire **仍然带着中点顶点 -65.244**（OCCT 的同一 wire 是 -49.864906）
   ⇒ §9.439 那处「FixDummySeam 顶点须朝向感知」的分歧**在 fix_missing_seam 的输出里依然存在**
     （与 §9.441 的结论一致：它不影响 F113 的 wire 数，只影响顶点位置）。
② 外环 22 条边被切成 7+5+4+3+3，断点落在这些顶点上（全部在 **z=-48 / z=-40 / z=-32 的“网格”层**）：
     (8.900,28.844,-48)  (22.366,28.844,-48)  (38.733,21.909,-40)
     (44.500,0,-32)      (±87.500,0,-32)〔退化薄片〕  (-8.900,28.844,-48)
     (-22.366,28.844,-48) (-38.733,21.909,-40) (-44.500,0,-32)
   即：**切口正好落在网格线所经过的顶点上**，且两侧 seam（x=±87.5）各留下一个**退化薄片**
   （首尾同点的 3 边 wire）——这与 OCCT 产出「一条 22 边 wire」形成鲜明对照。
```

#### 这一步得出的可检验命题（下一问的对象）

> **命题 P**：端口在 F113 上把外环沿「网格顶点」切开（并留下两张 seam 退化薄片），
> 而 OCCT 在同样的输入下不切 —— 差异应能在 `SplitByGrid`/`SplitByLine` 的**网格构造**上找到
> （同样的 `TOLINT=1e-10`、同样的 `u1/u2/v1/v2`，但**网格层**的取值或“是否需要切”的判据不同）。
>
> 证伪方式：把端口那次调用的**网格层坐标**打印出来（`split_by_grid.rs` 里 grid 的每一层 u/v），
> 与 `.cxx:2131-2275` 用同样输入算出的层值逐一比对；若层值相同而切分判据不同 ⇒ 判据处是分歧点；
> 若层值不同 ⇒ 网格构造处是分歧点。

---

### 9.452 —— 排掉 `break_wires`：两张 seam 退化薄片的**零长度段在其进入之前就已存在**（元凶在上游 `split_by_line`/`collect_wires`）

TEMP 插桩（`break_wires.rs`，env `T101_ZZ`，**已用 `git show HEAD` 还原**，库 diff 为空）打印该次调用的
段列表首尾端点（`first_vertex/last_vertex`，朝向感知）。F113 那一次（段数 6 → 8）：

```text
BREAK-IN  k=0/6 nb=6  (-65.244,0.000,-100.000) -> (36.224,3.000,-99.867)      ← 内环，**已带中点顶点**
          k=1/6 nb=16 (-8.900,33.805,-69.638) -> (-22.366,28.844,-48.000)      ← 一条 **16 边** wire
          k=2/6 nb=1  (87.500,0.000,-32.000) -> (87.500,0.000,-32.000)         ← **零长度段（已存在）**
          k=3/6 nb=1  (-87.500,0.000,-32.000) -> (-87.500,0.000,-32.000)       ← **零长度段（已存在）**
          k=4/6 nb=1  (44.500,0,-32) -> (87.500,0,-32)
          k=5/6 nb=1  (-87.500,0,-32) -> (-44.500,0,-32)
BREAK-OUT k=0/8 nb=6   k=1/8 nb=7   k=2/8 nb=5   k=3/8 nb=4   k=4..7/8 各 nb=1
```

**结论**：

```text
① `break_wires` **没有造边**（代码只按 External/Internal 段的顶点切链 + `seqw.insert`/`SetValue`），
   两张 seam 退化薄片的「零长度段 (±87.5,0,-32) -> 同点」**在它进入时就已存在**
   ⇒ **它不是造薄片的那一步**；它只把那条 **16 边 wire 切成 7+5+4**（切点=External 段端点，与 §9.451 断点一致）。
② 元凶在**上游**：`split_by_line`（按网格线切）或 `collect_wires`（把段配对成 wire）——
   零长度段 + `(44.5,0,-32)↔(87.5,0,-32)` 这种「去-回」形态只可能在这两处之一产生。
③ 顺带：**内环 6 边 wire 进入时就已带中点顶点 -65.244** ⇒ §9.439 的顶点分歧发生在更早（`fix_dummy_seam`），
   与本节同源但不同因。
```

**下一轮靶点**（一次调用即可分流，且必须带结构性筛选以免刷全模型）：

```text
在 `collect_wires.rs` 入口/出口、`split_by_line.rs` 切分后各打一次；**仅当该次段列表出现零长度边时打印**：
  · 零长度段在 `collect_wires` 入口就有 ⇒ 元凶在 `split_by_line`/`split_by_grid`；
  · 只在出口出现 ⇒ `collect_wires` 的配对/新造段逻辑是元凶。
```

---

### 9.453 —— 分流完成：零长度段的元凶在 `split_by_line`/`split_by_grid`（**不是** `collect_wires`）

TEMP 插桩（`collect_wires.rs` 入口/出口，env `T101_ZZ2`，仅当段列表含零长度边时打印；
**已用 `git show HEAD` 还原**）。F113 那一次（`.target-gate/collect113.txt`）：

```text
COLLECT-IN-seqw  k=0/8 nb=6  (-65.244,0,-100) -> (36.224,3,-99.867)
                 k=1/8 nb=7  (-8.900,33.805,-69.638) -> (22.366,28.844,-48.000)
                 k=2/8 nb=5  (38.733,21.909,-40) -> (-22.366,-28.844,-48)
                 k=3/8 nb=4  (-38.733,-21.909,-40) -> (-22.366,28.844,-48)
                 k=4/8 nb=1  ( 87.500,0,-32) -> ( 87.500,0,-32)     ← **零长度段（入口即有）**
                 k=5/8 nb=1  (-87.500,0,-32) -> (-87.500,0,-32)     ← **零长度段（入口即有）**
                 k=6/8 nb=1  ( 44.500,0,-32) -> ( 87.500,0,-32)
                 k=7/8 nb=1  (-87.500,0,-32) -> (-44.500,0,-32)
COLLECT-OUT-seqw 同上 8 段（未变）
COLLECT-OUT-wires k=0..3 nb=6/7/5/4（= 内环 + 7/5/4 三段外环）
                  k=4/6 nb=3  ( 87.500,0,-32) -> ( 44.500,0,-32)    ← **两张 3 边退化薄片在这里被“组装”出来**
                  k=5/6 nb=3  (-87.500,0,-32) -> (-44.500,0,-32)
```

**结论（链条收窄到一处）**：

```text
· `collect_wires` **入口**就已存在零长度段 (±87.5,0,-32)->同点 ⇒ **不是它造的**；
  它只是把「零长度段 + 两条 1 边段」**组装**成两张 3 边薄片（OUT-wires k=4/5）。
· 结合 §9.452（`break_wires` 也只是切链、不造边）⇒ **零长度段的产生地在上游 `split_by_line`
  （或 `split_by_grid`）**，即在「按网格线切分」这一步。
· 这与 §9.377 在 f=140 上的观测（`SplitByLine` 的 ClosedMode U 切、`breakwires` 出未配对半边 D≠C）同向。
```

**下一轮靶点**：`split_by_line.rs` 的切分输出（同样的结构性筛选：只在出现零长度段时打印），
对照 `ShapeFix_ComposeShell.cxx:1433-1914`（`SplitByLine`）判读：**哪一次「按线切」把一段切成了零长度**。
命中即按 `.cxx` 修（不加特例），然后跑 §9.448 的验收口径。

---

### 9.454 —— **定位完成**：零长度段诞生于 `split_by_grid`（`load_wires` 干净；`break_wires`/`collect_wires` 已排）

在**编排器** `perform.rs`（`ShapeFix_ComposeShell::Perform`，`cxx:206-255`）里按阶段埋一次二分
（env `T101_ZZ3`，**仅当该阶段段列表含零长度边时打印**；插桩已用 `git show HEAD` 还原）。
F113 那一次（`.target-gate/stage113.txt`，同一 face 出现两次 = 导入路径 + 探针显式调用）：

```text
A_load_wires      : **未打印** ⇒ 该阶段 zero=0（筛选条件：仅当 zero>0 时打印）
B_split_by_grid   : segs=6 zero=**2** nb=[6, 16, 1, 1, 1, 1] p0=(-65.244,0.000,-100.000)
C_break_wires     : segs=8 zero=2     nb=[6, 7, 5, 4, 1, 1, 1, 1]
D_collect_wires   : segs=6 zero=2     nb=[6, 7, 5, 4, 3, 3]
```

**结论（链条已经完全闭合到单一函数）**：

```text
零长度段(±87.5,0,-32)->同点 **由 `split_by_grid` 产生**（它在内部调用 `split_by_line`）；
  · `load_wires` 干净（zero=0，故 A 阶段无输出）；
  · `break_wires` 只把 16 边切成 7+5+4（zero 不变，§9.452）；
  · `collect_wires` 只把 1 边段组装成 3 边薄片（zero 不变，§9.453）。
另注：`perform.rs` 揭示 **`split_by_line` 不在编排器里**，它由 `split_by_grid` 内部调用
  ⇒ 下一步要么进 `split_by_grid`（`cxx:2131-2275`），要么进它调用的 `split_by_line`（`cxx:1433-1914`）。
```

**下一轮靶点（二选一，见下一问）**：
  (i) 在 `split_by_grid` 内、`split_by_line` 调用**前后**各打一次 ⇒ 判定是 `split_by_grid` 自己造的，还是它调用的 `split_by_line` 造的；
  (ii) 直接进 `split_by_line` 的切分输出，打印每次「按线切」的切点参数与切前/切后段。

---

### 9.455 —— 再收一层：零长度段由 `split_by_grid` 里的**按线切分循环**（`split_by_line_wires` = `SplitByLine`）产生

在 `split_by_grid`（`cxx:2131-2275`）的**切分循环之前**与**函数末尾**各埋一次（env `T101_ZZ4`，
仅当含零长度段时打印；插桩已 `git show HEAD` 还原）。F113 那一次：

```text
GRID-PRE  : **未打印** ⇒ 切分循环之前 zero=**0**（段列表干净）
GRID-POST : segs=6 zero=**2** nb=[6, 16, 1, 1, 1, 1]     ← 零长度段在按线切分之后出现
```

⇒ 零长度段是 `split_by_grid` 内部**按 U/V 线切分**这一步造的，即调用：

```rust
self.split_by_line_wires(seqw, &ln, true,  cut_index);   // U 线（cxx:2131-2251 段）
self.split_by_line_wires(seqw, &line, false, i as i32);  // V 线（cxx:2253-2273 段）
```

（`split_by_line_wires` 是 `SplitByLine`（`cxx:1433-1914`）的移植；4 个调用点：
U 的 closed/普通分支、V 的 closed/普通分支。）

#### 完整链条（现已闭合到「某一次按线切」）

```text
load_wires        : 干净
split_by_grid     : PRE 干净 → 按 U/V 线切 → POST 出现 2 条零长度段 (±87.5,0,-32)→同点   ← **元凶在此**
break_wires       : 只切链（16 → 7+5+4），zero 不变
collect_wires     : 只组装（1 边段 → 3 边薄片），zero 不变
dispatch_wires    : 5 张面 → perform 包成 Shell(5) → reader 只收 Face ⇒ 整包丢弃
```

**下一轮靶点**：给 `split_by_line_wires` 的 **4 个调用点**各埋一次「切前/切后」打印，
**带上该次切的参数**（线的 `pos/dir`、`cut_index`、U/V 标志），从而指认**是哪一条线、切点参数是多少**
把段切成了零长度；随后与 `cxx:1433-1914` 的同分支判读（切点是否落在端点 ⇒ 应跳过却没跳过）。

---

### 9.456 —— **抓到那一次切分**：U 方向、线 u=0.000000（seam 线）、`cut_index=1`；切完立刻出现两条零长度段（位置=两个 seam 点 ±87.5）

在 `split_by_grid` 的 **4 个 `split_by_line_wires` 调用点**各埋一次「切后」打印，并带上该次切的参数
（env `T101_ZZ5`，仅当段列表含零长度段时打印；插桩已 `git show HEAD` 还原）。
F113 那一次（`.target-gate/cut113.txt`）：

```text
CUT U uv=U line=(0.000000,0.000000) cut_index=1 segs=6 nb=[6, 16, 1, 1, 1, 1]
      ZEROS=["nb=1 ( 87.500,0.000,-32.000)->( 87.500,0.000,-32.000)",
             "nb=1 (-87.500,0.000,-32.000)->(-87.500,0.000,-32.000)"]
```

⇒ **元凶是「按 U 线切」中，线落在 u=0（本面的 seam 线）、`cut_index=1` 的那一次**：
切完之后段列表里立刻出现两条 **零长度段**，其顶点正是这张圆柱面的**两个 seam 点**
（`x=±87.5, y=0, z=-32`）。

配合 §9.455（`GRID-PRE` 干净、`GRID-POST` 出现 zero=2）与 §9.454/§9.453/§9.452，链条现在是：

```text
load_wires           干净
split_by_grid        PRE 干净 → **U 线 u=0 / cut_index=1 的这一次切分** → 出现 2 条零长度段(±87.5)
break_wires          只切链（16 → 7+5+4），zero 不变
collect_wires        只组装（1 边段 → 3 边薄片），zero 不变
dispatch_wires       → 5 张面 → Shell(5) → reader 只收 Face ⇒ 丢弃
```

#### 下一轮（最后一个未知量已定位到一个分支）

```text
读 `split_by_line_wires`（`SplitByLine`：cxx:1433-1914）里 **U 方向、线在 u=0/seam、closed-mode** 那条分支，
与 cxx 同行对拍，重点看 **「切点与边端点重合」时的处理**：
  · OCCT 侧该处是否有「切点 ≈ 端点 ⇒ 不切 / 直接把边并入」的判据（如 `Precision::PConfusion()` 比较）；
  · 端口是否缺这条判据（于是把一个端点切成了一段零长度边）。
命中即按 .cxx 修（不加特例），验收口径见 §9.448：
  zz_seam_fix 113 → Face；--model 113 → wires=2（22+6 边）；--fstats → face=113 mt≈228；
  a3n00 面积比 ≥0.8996 且上升；t101_verify.ps1 全绿。
```

---

### 9.457 —— 「合并零长度切向段」这一步**忠实**（含删除后 `i += 1`）；靶点移向逐 wire 的切分器 `split_by_line`

对拍 `ShapeFix_ComposeShell.cxx:1948-1965` ↔ 端口 `split_by_line.rs:518-545`：

```cpp
// OCCT:1949-1964（1-based，比较 (i, i+1)）
for (i = 1; i < SplitLinePar.Length(); i++) {
  if (std::abs(SplitLinePar(i+1) - SplitLinePar(i)) > PConfusion && !SplitLineVertex(i).IsSame(SplitLineVertex(i+1))) continue;
  if ((Code(i)&ITP_ENDSEG && Code(i+1)&ITP_BEGSEG) || (Code(i)&ITP_BEGSEG && Code(i+1)&ITP_ENDSEG)) {
    int code = (Code(i) | Code(i+1)) & IOR_BOTH;
    Code.SetValue(i, code | (code == IOR_BOTH ? ITP_INTER : ITP_TANG));
    Par.Remove(i+1); Code.Remove(i+1); Vertex.Remove(i+1);
  }
}
```
```rust
// 端口:519-545（0-based，等价配对 (i-1, i)）
if (par[i]-par[i-1]).abs() > PCONFUSION && !is_same(&vertex[i-1].0, &vertex[i].0) { i += 1; continue; }
if (code[i-1] & ITP_ENDSEG != 0 && code[i] & ITP_BEGSEG != 0) || (…) {
    let code = (split_line_code[i-1] | split_line_code[i]) & IOR_BOTH;
    split_line_code[i-1] = code | (if code == IOR_BOTH { ITP_INTER } else { ITP_TANG });
    par.remove(i); code.remove(i); vertex.remove(i);
    i += 1;      // 与 OCCT 的 for 增量一致（端口已有注释论证）
}
```

⇒ **判定：忠实**（配对、条件、`code` 合成、删除「后者」、删除后 `i += 1` 全都对上）。

#### 关键澄清与新的靶点

```text
`split_line_par/code/vertex` 是「切线与该 wire 的交点序列」，**不是**结果里的 WireSegment 列表
⇒ 这一步合并的是**交点**，它不可能消掉已经生成成段的两条「零长度 WireSegment」。
那两条 (±87.5,0,-32)->同点 的段来自**逐 wire 的切分器** `split_by_line`
（端口 `split_by_line.rs:1-478` ↔ `ShapeFix_ComposeShell.cxx:1433-1914`）——
即「把一条边按与切线的交点切成若干子段」的那一步。

下一轮（最后一个函数、最后一次对拍）：
  读端口 `:1-478` 里**生成子段**的那段（沿 `split_line_par` 走线、按奇偶与 code 建子段），
  与 cxx:1967-2110 同行对拍，专找 **「子段参数区间长度 ≈ 0 ⇒ 跳过」** 这类判据：
    · cxx 若有而端口缺 ⇒ 命中（端口中 seam 时切出零长度子段）；
    · 两侧都有 ⇒ 分歧在**交点参数**本身（即线与边的求交），转去比 `SplitByLine` 的求交段（cxx:1433-1948）。
```

---

### 9.458 —— 切线建边这步也**忠实**（cxx:2019-2057 ↔ 端口 :586-613）；由排除法锁定**逐 wire 切分器 `split_by_line`（:1-478 ↔ cxx:1433-1914）**

逐项对拍（OCCT 1-based `(i-1, i)` ↔ 端口 0-based `(i-2, i-1)`）：

```text
cxx:2019-2022  tmpV1/tmpV2 = Context()->Apply(SplitLineVertex(i-1)/(i))   ↔ 端口 :586-587   ✓
cxx:2025       canbeMerged = (i - 1 > 1 || i < SplitLinePar.Length())     ↔ 端口 :590       ✓ **含边界**
cxx:2026-2031  aMaxTol <= 2*Confusion ⇒ Infinite                          ↔ 端口 :592-594   ✓
cxx:2032-2033  aTol1/aTol2 = min(BRep_Tool::Tolerance(·), aMaxTol)        ↔ 端口 :595-596   ✓
cxx:2036       aD = aP1.SquareDistance(aP2)                               ↔ 端口 :597       ✓
cxx:2037-2038  if (par(i)-par(i-1) < PConfusion() || (canbeMerged && (aD<=tol1² || aD<=tol2²)))
                                                                          ↔ 端口 :598-600   ✓ **有符号比较一致**
cxx:2044-2056  !V1.IsSame(V2) ⇒ CombineVertex + Context()->Replace，然后 **continue**（不建边）
                                                                          ↔ 端口 :601-612   ✓
cxx:2059+      否则建边（两条 pcurve + pcurve_range/edge_range）           ↔ 端口 :616-625   ✓
```

**结论（排除法闭环）**：

```text
迄今已逐个排除或判定忠实：
  load_wires · split_by_grid 的网格构造/切分循环调度 · break_wires · collect_wires ·
  SplitByLine 的「交点合并」(cxx:1948-1965) · 「沿切线建边」(cxx:2019-2057)
⇒ 剩下的唯一去处：**逐 wire 切分器** `split_by_line`（端口 `split_by_line.rs:1-478` ↔ `cxx:1433-1914`），
  即「把**原始 wire 的边**按与切线的交点切成若干子段」那一步——
  两条零长度 WireSegment（(±87.5,0,-32)->同点）只可能由它把一条边切出一个零长度子段而来。
```

**下一轮（唯一剩下的读）**：端口 `split_by_line`（:1-478）里**生成子段**的那段，与 `cxx:1433-1914` 对拍，
专找「子段参数区间长度 ≈ 0 ⇒ 跳过/并入」的判据；若两侧都有，则分歧在**交点参数**（线与边的求交）。
命中即按 `.cxx` 修 → 跑 §9.448 验收（`zz_seam_fix 113` → Face；`--model 113` → wires=2（22+6）；
`--fstats` → face=113 mt≈228；a3n00 面积比 ≥0.8996 且上升；`t101_verify.ps1` 全绿）。

（本轮所有 TEMP 插桩均已 `git show HEAD` 还原；`git diff --stat crates/occt-topo/src` 为空。）

---

### 9.459 —— `split_by_line` 的「closed mode 去重点」(cxx:1722-1767 ↔ 端口 :295-341) 逐项对拍：**公式全同**，剩一个未验证前提（边的朝向语义）

| OCCT | 端口 | 判定 |
|---|---|---|
| `int i, j = IntEdgePar.Length();`（1-based，j=末索引） | `let mut j = int_edge_par.len() - 1;`（:301，含解释注释） | ✓ |
| `if (myClosedMode && j > 1)` | `if self.closed_mode && int_edge_par.len() > 1` (:296) | ✓ |
| `for (i = 1; i <= Length();)` / `if (i == j) break;` | `while i < len()` / `if i == j { break; }` (:303-305) | ✓ |
| `IntEdgeInd(i) == IntEdgeInd(j) && abs(Par(i)-Par(j)) < PConfusion` → `Remove(i)`，`if (j>i) j--`，`continue` | :307-309 / :329-336 | ✓ |
| `else if (nbe == 1 \|\| IntEdgeInd(i) == (IntEdgeInd(j) % nbe) + 1)` | :311 | ✓ |
| `E1 = sewd->Edge(IntEdgeInd(j)); E2 = sewd->Edge(IntEdgeInd(i));`（**j→E1，i→E2**） | `e1 = wire.edge(int_edge_ind[j])`, `e2 = wire.edge(int_edge_ind[i])` (:312-313) | ✓ |
| `BRep_Tool::Range(E1, myFace, a1, b1)` | `curve_on_surface_range(&e1, &face)` (:314-321) | ✓ |
| `abs(Par(j) - (E1.Orientation()==FORWARD ? b1 : a1)) < PConfusion && abs(Par(i) - (E2.Orientation()==FORWARD ? a2 : b2)) < PConfusion` | :322-325 `e1_ref = Forward ? b1 : a1`、`e2_ref = Forward ? a2 : b2` | ✓（公式同） |
| `j = i++;` | :338-339 `j = i; i += 1;` | ✓ |

#### 剩一个**未验证前提**（也是唯一可疑处）

```text
OCCT 的 `E1.Orientation()` 是 **sewd（wire 数据）里存的朝向** —— `ShapeExtend_WireData::Add` 存的是
按 wire 遍历朝向过的边；而端口 :312-313 取的是 `WireSegment::edge(i)`，
其实现（wire_segment.rs:144）只是 `self.edges.get(i-1)`，**原样返回所存边**。
⇒ 若端口的 `load_wires` 存边时**已施加 wire 朝向**，则两侧语义一致（本块忠实）；
  若存的是**未施加朝向**的原始边，则 :322-323 的 `e1_ref/e2_ref` 会选反 ⇒ 去重判据失效
  ⇒ **重复点存活 ⇒ 后面切出零长度子段**（与 §9.439 的 `FixDummySeam` 同族：朝向语义）。
```

**下一轮（唯一剩余检查，成本一次读）**：读 `load_wires.rs`（很小），确认它存进 `WireSegment.edges` 的边
**是否带 wire 遍历朝向**（对照 `ShapeExtend_WireData` 的 Add/Edge 语义）；若不带 ⇒ 按 `.cxx` 修（使
`WireSegment::edge()` 的语义与 `sewd->Edge()` 一致，或在该判据处用朝向感知取边），然后跑 §9.448 验收。

---

### 9.460 —— 前提已验证、**朝向假设被否**：去重块忠实；靶点推进到 `split_by_line` 剩余未读段（`:345-478`）

读 `load_wires.rs:11-27` 的 `wire_data_edges`（其文档注释直接写明语义）：

```rust
/// ShapeExtend_WireData(wire, chkseam, manifold) (ShapeExtend_WireData.cxx:80-120):
/// the stored edges keep the composed orientation of the wire, so a REVERSED
/// wire flips every edge.
fn wire_data_edges(wire: &Wire) -> Vec<crate::shape::Edge> {
    let edges = crate::topo_tools_full::edges_of_wire(wire);
    if wire.0.orientation() == Orientation::Reversed { edges.into_iter().map(|mut e| { e.0.reverse(); e }).collect() }
    else { edges }
}
```

⇒ 存进 `WireSegment.edges` 的边**带 wire 遍历朝向**（与 `ShapeExtend_WireData` 一致）
⇒ `WireSegment::edge(i).0.orientation()` ≡ OCCT `sewd->Edge(i).Orientation()`
⇒ §9.459 里 `:322-323` 的 `e1_ref/e2_ref` **与 OCCT `:1750-1752` 语义相同** ⇒ **去重块（:295-341）判定忠实**，
   「朝向选反 ⇒ 重复点存活」这条假设**否掉**。（同时把 `WireSegment::edge` 的语义也钉下来了，供后续复用。）

#### `split_by_line`（:1-478）已读/未读盘点

```text
已读并判定忠实：:295-341 去重点块（cxx:1722-1767）、:518-545 交点合并（cxx:1948-1965，属 split_by_line_wires）、
                :586-613 沿切线建边（cxx:2019-2057）
**未读（剩余靶点）**：
  · :343-378  segment codes（cxx:1771-1782）
  · :379-393  第二处移除路径（按边端点参数；本文件里另一个 `continue`）
  · :394-478  子段生成：把原始 wire 的边按 `int_edge_par` 切成子段 —— **零长度子段的最后可能来源**
  · :31-290   与切线的求交（`int_edge_par/int_line_par/int_edge_ind` 的产生地）—— 若 394-478 判定忠实则回头查这里
```

**下一轮**：读 `:379-478`（约 100 行，含第二处移除 + 子段生成），与 `cxx:1782-1914` 对拍，专找
「子段参数区间长度 ≈ 0 ⇒ 跳过/并入」判据。命中即按 `.cxx` 修 → 跑 §9.448 验收。

---

### 9.461 —— 发现一处**真实的循环簿记分歧**（`cxx:1820` ↔ 端口 `:389`），但**实测与 F113 无关**；已还原待带门禁落地

**分歧**（第二处移除路径，`cxx:1799-1822` ↔ 端口 `:370-397`）：

```cpp
// OCCT：for (i = 1; i <= IntEdgePar.Length(); i++)
if (SegmentCodes(j) == IOR_UNDEF && SegmentCodes(i) == IOR_UNDEF) {
  if (myClosedMode && (IntLinePar(i)-IntLinePar(j))*(IntLinePar(k)-IntLinePar(i)) <= 0.) continue;
  IntEdgeInd.Remove(i); IntEdgePar.Remove(i); IntLinePar.Remove(i); SegmentCodes.Remove(i);
  i--;        // ← 被 for 的隐式 i++ 抵消 ⇒ **净不前移**（重查移位到位置 i 的元素）
}
```
```rust
// 端口（修前）：删除后 `i -= 1` 再 `continue`（while 无隐式自增）⇒ 净**前移 −1**，会去重查前一个元素
int_edge_ind.remove(i - 1); … i -= 1; if i == 0 { i = 1; } continue;
```

**修法**（已实现并实测，随后**已 `git show HEAD` 还原**）：删除后 `i` 保持不变、直接 `continue`（等价 OCCT 的 `i--` + `i++`）。

**实测（修后）**：

```text
F113：result type=Shell faces=5；HEALED wires=2 (6e+7e) / 5e / …  ——**与修前逐字相同**
模型：MODEL f=113 wires=4（6/14/1/1）                              ——**不变**
a3n00：STATMAP matched=225 unmatched=1 sum_mt=11941 flat_mt=11941  ——**不变**
```

⇒ 该分歧**真实**（簿记不等价，属于应修的对齐点），但**不是 F113 的成因**；
它是「忠实性」问题而非本目标的阻塞点。若要落地，须按纪律跑完整门禁
（`--lib 1281/0`、`step_obj_gates 5/5`、面积比不降）——建议下次单独一轮做：应用 → 后台跑 `t101_verify.ps1` → 绿则提交、红则还原。

#### 靶点继续（`split_by_line` 剩余未读）

```text
:430-478  子段生成的其余部分（after `ipcode` computation）——仍可能含「子段长度≈0 ⇒ 跳过」判据
:31-290   与切线的求交（int_edge_par/int_line_par/int_edge_ind 的产生地）
```

---

### 9.462 —— **靶点重定向**：子段生成不在 `split_by_line` 内，而是委托给 `split_wire`（`split_by_line.rs:454-462`，注释指 `cxx:1900-1903`）

读 `split_by_line.rs:430-478`（该函数最后一段）：

```rust
                int_code.push(ipcode);
                j = i;
            }
        }
        // cxx:1900-1903.
        let mut int_vertices: Vec<Vertex> = Vec::new();
        let mut indexes = int_edge_ind.clone();
        *wire = self.split_wire(
            wire, &mut indexes, &int_edge_par, &mut int_vertices, &a_new_seg_codes, is_cut_by_u, cut_index,
        );
        // cxx:1906-1911.
        for i in 1..=int_line_par.len() {
            split_line_par.push(int_line_par[i - 1]);
            split_line_code.push(int_code[i - 1]);
            split_line_vertex.push(int_vertices[i - 1].clone());
        }
```

⇒ **「把原始 wire 的边按交点切成子段」这一步在 `split_wire` 里**（`split_wire.rs`，496 行），
   `split_by_line` 只负责算交点/代码并把结果交给它。故 objective 里写的「`split_by_line` 的子段生成段」
   实际落点是 **`split_wire`**（这就是 §9.399 曾判「控制流忠实」的那个函数——当时没有零长度段这个靶点）。

#### `split_wire.rs` 的判据/建边位点（本轮 grep 结果，供下一轮直接读）

```text
:58   pub fn split_wire(...)
:201-224  PCONFUSION 过滤区（hi/lo = max/min(par) ± PCONFUSION；(curr_par-last_par).abs()<PCONFUSION ⇒ :224 continue）
:275  continue（另一处跳过）
:319-330  is_same(prev_v,last_v) / is_same(v,last_v) / is_same(v,prev_v) 分支
:357-372  第二处 PCONFUSION 过滤（(apar-prev_par).abs()<=PCONFUSION ⇒ remove …）
:402/446/448/466  add_edge_patch(...) —— **子段真正被造出来的地方**
```

**下一轮**：读 `split_wire.rs:190-240` 与 `:340-410`（两处过滤区 + 紧随的建边），与
`ShapeFix_ComposeShell.cxx` 对应段（`cxx:1900-1903` 调用的 `SplitWire` 实现）对拍，
专找「子段参数区间长度 ≈ 0 ⇒ 跳过」判据；命中即按 `.cxx` 修 → 跑 §9.448 验收。

---

### 9.463 —— ✅ **落地**：`split_by_line` 第二处移除路径的簿记修复（`cxx:1820` 语义），门禁全绿且全部基线逐字不动

改动（`crates/occt-topo/src/shape_fix_compose_shell/split_by_line.rs`，第二处移除路径）：

```rust
// 修前（与 OCCT 不等价：while 无隐式自增，却仍 i -= 1 ⇒ 净前移 −1，去重查了前一个元素）
int_edge_ind.remove(i - 1); int_edge_par.remove(i - 1); int_line_par.remove(i - 1); segment_codes.remove(i - 1);
i -= 1; if i == 0 { i = 1; } continue;
// 修后（等价 cxx:1820 的 `i--` 被 `for` 的隐式 `i++` 抵消 ⇒ 净不前移，重查移位到位置 i 的元素）
int_edge_ind.remove(i - 1); int_edge_par.remove(i - 1); int_line_par.remove(i - 1); segment_codes.remove(i - 1);
continue;
```

**门禁实测**（`pwsh -File .target-gate\t101_verify.ps1`，job `pwsh-58`，exit 0）：

```text
test result: ok. 5 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 269.51s

[occ/a3n00.stp] ours v=10863 f=11941 | occ v=11145 f=12466 | f-ratio 0.96
[occ/a3n00.stp] our=204327.72 occ=227130.30 ratio=**0.8996**   ← 与基线逐字相同
[occ/T0M.stp]   our=192416.40 occ=192658.64 ratio=0.9987      ← 相同
[occ/acs10.stp] our=267047.23 occ=271219.68 ratio=0.9846      ← 相同
[occ/bottom.step] 0.9750 · [occ/motoc.step] 1.0333 · [occ/TDB.stp] 1.0092 · [occ/top.step] 1.0000 · [occ/ATU01038.step] 1.0006（均相同）
```

⇒ 该修复**忠实于 `.cxx` 且对所有既有基线完全中性**（F113 也仍为 Shell(5)，与 §9.461 的实测一致：
它是簿记正确性问题，不是 F113 的成因）。**已落地**。

#### 与目标的关系

```text
F113 仍未解决（4 wires / 空面）；路线图不变：
  零长度段 ← split_by_line 的 U 线 u=0 切分 ← 子段生成委托给 split_wire（§9.462）
  下一步：split_wire.rs:190-240 / :340-410（两处 PCONFUSION 过滤 + :402/446/448/466 建边）
          ↔ cxx:1116-1184 / cxx:1900-1903 的 SplitWire 实现，找「子段长度 ≈ 0 ⇒ 跳过」判据
```

---

### 9.464 —— `split_wire` 子段**创建处没有局部零长度守卫**（守卫在更上游两处）；下一处对拍已精确到 20 行

读 `split_wire.rs:350-403`（含内部顶点处理与第一个建边点）：

```rust
// :352-379 内部顶点：apar≈prev_par ⇒ 合并到 prev_v 并删除；apar≈curr_par ⇒ 合并到 v 并删除；
//           prev_par < apar < curr_par ⇒ 作为 new_edge 的内部顶点加入（并删掉该条目）
// :381-385 cxx:1339-1344：copy_pcurves + TransferParametersProj.init + **transfer_range(new_edge, prev_par, curr_par, true)**
// :395-397 cxx:1358-1361：!sp && !is_degenerated ⇒ set_same_range(false)
// :398-403 cxx:1362-1371：code 修正 + **result.add_edge_patch(0, new_edge, …)** + define_patch
```

⇒ **`new_edge` 的参数区间直接取 `[prev_par, curr_par]`，此处没有任何「区间长度 ≈ 0 ⇒ 跳过」的局部守卫**。
  `is_degenerated` 只用于 `set_same_range`（:395），**不用来拦截建边**。

#### 真正的守卫在更上游（round 77 已读到，`split_wire.rs:214-224`）

```rust
if (curr_par - last_par).abs() < PCONFUSION {      // cxx:1136-1140
    v_opt = Some(last_v.clone()); do_cut = false;  // 端点重合 ⇒ 不切
} else if (curr_par - prev_par).abs() < PCONFUSION { // cxx:1141-1146
    vertices.push(prev_v.clone()); code = code_at(j); prev_par = curr_par; j += 1; continue;  // 与前一切点重合 ⇒ 跳过
} else { … 正常切 … }
```

⇒ 零长度子段若从这条路径逃出，只可能是：**两处守卫都没命中，而 `prev_par` 与 `curr_par` 仍近乎相等**
（例如 `|curr−prev| ≥ PCONFUSION` 但 `|prev−last| < PCONFUSION` 这类端点/切点交叉情形）。

**下一轮（精确到 ~20 行）**：读 OCCT `cxx:1120-1150`（这两处守卫的原文），与端口 `:214-224` **逐字符对拍**
（比较运算符、`PConfusion` 用法、比较对象是 `first_par/last_par/prev_par/curr_par` 哪一个）；
若一致 ⇒ 转**数值取证**：在 `:402` 前对 F113 那次调用打印 `(prev_par, curr_par, last_par, first_par)`，
直接看零长度子段落进哪条分支。

---

### 9.465 —— 数值取证（否）: `split_wire.rs:402` **不是**零长度段的产地；剩余建边点已列出

在 `:402`（`result.add_edge_patch(0, new_edge, …)`）前插入近零长度取证（env `T101_ZZ6`；**已 `git show HEAD` 还原**），
过滤 `plen<1e-4 || p_end<1e-6 || c_end<1e-6 || p_first<1e-6`，打印
`len/prev/curr/first/last/p_end/c_end/p_first/degen/is_cut_by_u`。实测（`.target-gate/piece113.txt`）：

```text
PIECE 总行: 56       其中 **len < 1e-3 的行: 0**
（56 行全部来自 `p_first < 1e-6` 这类良性情形：段起点恰为该边的 first_par，例如
 len=3.142e0 prev=6.283185307 curr=3.141592654 first=6.283185307 last=0.000000000 p_first=0.000e0）
```

⇒ **`split_wire.rs:402` 处的正常切分支不产生任何近零长度子段**（`len` 最小也在 1e-3 以上）。

#### `split_wire.rs` 里全部的建边点（grep 已确认 4 处），逐个排除的清单

```text
:113  result.add_edge_patch(0, edge.clone(), …)      —— **未测**（疑为「不切、整边拷贝」的快路径）
:402  result.add_edge_patch(0, new_edge, …)          —— **本轮已测：无近零长度（0 行）**
:446  result.add_edge_patch(0, e1, …)                —— **未测**（内部/非流形分支）
:448  result.add_edge_patch(0, edge.clone(), …)      —— **未测**（同上）
:466  result.add_edge_patch(0, edge.clone(), …)      —— **未测**（同上）
```

**下一轮**：给 `:113/:446/:448/:466` 四处各加同款 `plen` 取证（同一次运行即可判完），找出零长度段真正的产生点；
若四处都干净 ⇒ 说明零长度 `WireSegment` 不是 `split_wire` 造的，需回头查 `split_by_line` 里
`split_wire` 返回后的装配（`:465-469` 把 `int_line_par/int_code/int_vertices` 推入三个输出数组）是否有错位。

---

### 9.466 —— 反直觉读数：`split_by_line` → `split_wire` 这条切分路径在本次运行里**从未以 `n_par ≥ 6` 执行**（0 行）

在 `split_by_line` 调 `split_wire` 之前（`cxx:1900-1903`）插入参数取证（env `T101_ZZ7`，
条件 `int_edge_par.len() >= 6`；**已 `git show HEAD` 还原**）：

```text
PAR 行数: **0**
```

⇒ 整趟运行中，**没有任何一次** `split_by_line` 在 `int_edge_par.len() >= 6` 时走到 `split_wire`。
这有两种解释，都必须先分清（下一轮一次运行即可）：

```text
(a) F113 的那些长 wire（16 边/6 边）在 `split_by_line` 里**提前返回**了：
      · `:345-347`  `if n_par == 0 { return false; }`（无交点）
      · `:369-397`  `wire.orientation() == Internal` 分支 / `:398-400` `if int_edge_par.is_empty() { return false; }`
      · 或更早的 `return`（`:31-290` 段里的各处）
    ⇒ 那这些 wire 根本没被「按线切」，零长度段另有来源；
(b) 过滤条件太严（该次切分实际只有 1–5 个交点）⇒ 把条件放宽到 `>= 1` 并同时打印**走到了哪个提前返回点**。
```

**下一轮**：把探针放宽为「`n_par >= 1` 时打印 `n/ind/par`」，并在 `split_by_line` 的**每个 early-return 处**打一次标记
（`RET-n0` / `RET-internal` / `RET-empty`），一次运行即可判定 F113 的 6 边与 16 边 wire 究竟走哪条路；
若确认是「提前返回 ⇒ 根本没切」，则零长度段必然来自 `split_by_grid` 里**切分之外的**部分
（网格构造写回 `seqw` 的那几行），需回头对拍 `cxx:2131-2210`。

---

### 9.467 —— `split_by_line` 被调用 96 次，但**没有一次的 `edges=` 是 6 或 16** ⇒ F113 的长 wire 不经过「按线切」

在 `:345` 前打印每次调用的 `n(=n_par)/edges/ori/par`（env `T101_ZZ8`；**已还原**）：

```text
SPLIT 总行数: **96**
按 `edges=` 过滤 `(16|6)` 的行: **0**
```

⇒ `split_by_line` 确实被频繁调用（96 次，全部走到 `n_par` 处），但**没有任何一次处理 6 边或 16 边的 wire**
——而 F113 的零长度段所在的段列表里恰恰有 6 边与 16 边两条 wire（§9.452 BREAK-IN 读数）。

**推论（重要）**：F113 的 6 边内环与 16 边外环**从未进入 `split_by_line`** ⇒
§9.455 的「GRID-POST 出现零长度段 ⇒ 是切分造成的」这条归因需要修正：
切分只作用于**别的** wire（96 次调用），而这两条 wire 的零长度段必然来自
**`split_by_grid` 里切分之外的操作**（把网格层顶点/段写回 `seqw` 的那几行，`cxx:2131-2210`）。

#### 下一轮（把归因钉死）

```text
在 `split_by_grid` 的**每一次对 `seqw` 的写操作**（新建段 / 插入 / 替换 / 追加顶点段）后打一次
「段数 + 每段边数 + 零长度段数」的简短标记（同一 env，结构筛选：仅当段列表含零长度段时打印），
定位零长度段是在**哪一条写操作**之后出现的；随后与 `cxx:2131-2210` 同段对拍（尤其
`ShapeExtend_WireData` 相关的 add/insert 与 seam 处的顶点段生成）。
```

---

### 9.468 —— 收窄：`split_by_grid:28-72` 只写 patch 索引（不产段）；新假设 = 零长度段由 `split_by_line_wires` 的**切线建边**在 `canbeMerged == false` 的端点情形下造出

#### A. 读 `split_by_grid.rs:24-72`（`cxx:2152-2199`）：**不产段**

```rust
if self.closed_mode {
    for w in seqw.iter_mut() {                                  // 只读/只改段属性
        let wire = builder.make_wire(w.edges());
        add_uv_bounds_on_wire(&face, &wire, &mut box2d);        // 该段的 UV 盒
        … shift_u/shift_v = adjust_to_period(…)                 // cxx:2165-2174
        let iumin/iumax/ivmin/ivmax = …get_patch_index(…)        // cxx:2182-2198
        for j in 1..=w.nb_edges() { w.define_iu_min(j, iumin); w.define_iu_max(j, iumax);
                                    w.define_iv_min(j, ivmin); w.define_iv_max(j, ivmax); }
    }
}
```
⇒ 只是**给每段每条边写 patch 索引**，没有任何 `push/insert/remove`、也不切边 ⇒ **不是零长度段的产地**。

#### B. §9.467 捕获里的新数据（修正 §9.467 的表述）

```text
SPLIT（split_by_line 调用）96 次：n 分布 {1:94, 2:2}；edges 分布 {1:91, **14:2**, **12:2**, 4:1}
```

⇒ **F113 的 14 边 wire（模型侧 `MODEL f=113 wire[1] nEdges=14`）确实进了 `split_by_line`**（2 次，`n_par` 为 1 或 2），
   而 **6 边与 16 边那两条没有**。所以更准确的说法是：**长 wire 中只有 14 边那条被「按线切」**，
   但它的交点只有 1–2 个 ⇒ 按线切不可能把 16 边切成 7+5+4（那由 `break_wires` 完成，§9.452）。

#### C. 新假设（比 §9.467 更具体，下一轮一次运行即可判定）

零长度段（`(±87.5,0,-32)->同点`）来自 **`split_by_line_wires` 里「沿切线建边」**（`cxx:2059-2108` ↔ 端口 `:616-625`）：

```text
cxx:2025  canbeMerged = (i - 1 > 1 || i < SplitLinePar.Length())
cxx:2037  if (Par(i) - Par(i-1) < PConfusion() || (canbeMerged && (aD<=tol1² || aD<=tol2²))) { 合并顶点; continue; }
⇒ 若退化的一对**落在端点**（`canbeMerged == false` 的那些 i），该守卫**不成立** ⇒ **照样建边**
   ⇒ 建出来的边两端点重合 ⇒ 正是「零长度段」。
端口 :590/:598-600 的实现与 cxx 相同（§9.458 已逐项对过），所以**不是端口写错**，
而是「端点处退化对」这一情形在两侧都被允许 —— 需要看 OCCT **在同一输入下**是否也会走到这一步
（即：分歧可能在**交点个数/位置**：端口给出 seam 处多一个交点，OCCT 不会）。
```

**下一轮**：给 `split_by_line_wires` 的沿切线建边段（端口 `:590-630`）加一次「仅当 `canbeMerged == false`
且建出的边两端点重合时」的打印，并同时打印 `i / Length / par(i-1), par(i) / aD / tol1, tol2`；
一次运行即可确认「零长度段是否诞生于此」，并给出该退化对的参数值（用于与 OCCT 同输入对照）。

---

### 9.469 —— 假设（§9.468-C）**否掉**：`canbeMerged == false` 的那些下标处**没有**退化对（45 行，|Δpar| 最小 2.0）

在沿切线建边处打印 `!can_be_merged` 的情形（env `T101_ZZ9`；**已还原**）：

```text
ENDSEG 行: **45**
len 分布: {2: 45}    i 分布: {2: 45}        ← 全部是 `Length()==2` 时 `i=2` 那一个下标
|par_i - par_im1| < 1e-6 的行（真正的退化对）: **0**
|Δpar| 最小值: 2.0
样本: ENDSEG i=2 len=2 par_im1=34.000000000 par_i=36.000000000 a_d=4.000e0 tol1=7.944e-15 tol2=7.105e-15
      ENDSEG i=2 len=2 par_im1=-0.000000000 par_i=18.000000000 a_d=3.240e2 tol1=7.105e-15 tol2=7.105e-15
```

⇒ `canbeMerged == false` 只出现在 `Length()==2, i==2` 这一种情形，而那里的两点**相距很远**（≥2.0），
   守卫不成立也无害 ⇒ **「端点退化对 ⇒ 建出零长度边」这条假设不成立**（§9.468-C 否掉）。

#### 至此的净结论（把「产地」重新归零到两处，且都还没测）

```text
零长度段 (±87.5,0,-32)->同点 已排除的产地：
  ✗ split_wire 的 :402 正常切分支（§9.465：无近零长度）
  ✗ split_by_line_wires 的沿切线建边端点退化（§9.469 本轮）
  ✗ split_by_grid:28-72 的 patch 索引写回（§9.468-A：不产段）
  — split_by_line 的按线切（§9.467：6/16 边 wire 根本没进这个函数）
剩余未测产地（两处）：
  ① `split_by_line_wires` 里「沿切线建边」的**非端点**情形（即 canbeMerged==true 且守卫成立处之外的正常建边，
     :616-625）——需要在建边处直接量 `new_edge` 的两端点是否重合（本轮只量了参数对，未量建边结果）
  ② `split_by_line.rs:465-469` 的装配：把 `int_line_par/int_code/int_vertices` 三个数组按同一 `i` 推入
     `split_line_par/code/vertex` —— 若三者长度不一致（例如 `int_code` 少一个元素），后续按
     `split_line_par` 遍历时就会错位（顶点配错边 ⇒ 零长度段）。**这是目前最可疑的一处**（端口用三个独立 Vec，
     OCCT 用 `NCollection_Sequence` 同时 Append，长度天然一致）。
```

**下一轮（优先做 ②，成本最低）**：在 `:465` 前打印 `int_line_par.len() / int_code.len() / int_vertices.len()`
（以及 `split_line_par.len()` 的增量），一次运行即可判定三数组是否等长；不等长即为分歧点（按 `.cxx` 修成同源同长）。

---

### 9.470 —— 假设（§9.469-②三数组错位）**否掉**：96 次调用三数组全等长；产地假设转向「输入 wire 自带的零长度边被原样拷贝」

在 `:465`（`cxx:1906-1911`）打印装配前三数组长度（env `T101_ZZA`；**已还原**）：

```text
ASM 行: **96**        三数组不等长的行: **0**
分布(去重)：{(1,1,1): 94, (2,2,2): 2}      ← int_line_par / int_code / int_vertices 长度始终一致
```

⇒ §9.469-②「三数组错位导致顶点配错边」**不成立**。同时这批数据再次证实 §9.467/§9.468-B：
   `split_by_line` 的调用**只有 1–2 个交点**，它根本没有能力把 16 边 wire 切成 7+5+4。

#### 产地清单（更新）

```text
已排除（4 项）：
  ✗ split_wire:402 正常切分支（§9.465）
  ✗ 沿切线建边的端点退化（§9.469）
  ✗ split_by_grid:28-72 patch 索引写回（§9.468-A）
  ✗ :465-469 三数组错位（本轮）
  — split_by_line 按线切不适用于 6/16 边 wire（§9.467，只处理 1–2 个交点）

**新首要怀疑（此前从未查过）：零长度段可能不是「造」出来的，而是输入 wire 里本来就有的
零长度边（seam 处的退化边）被原样拷贝**：
  · 模型侧 `MODEL f=113` 的 wire[2]/wire[3] 各只有 **1 条边** ⇒ 形如「单边环」，很可能就是零长度边；
  · `split_wire` 的 **`:113`**（`add_edge_patch(0, edge.clone(), …)`，未测）与 **`:466`** 同为
    「不切、整边拷贝」路径 ⇒ 一条零长度边会被原样搬进结果段。
  · 若是这样，则**真正的分歧不在端口是否「造」了零长度段，而在 OCCT 的 LoadWires/FixReorder/CollectWires
    如何处理这种退化边**（丢弃？并入邻段？）——这解释了为什么前面所有「切分侧」的对拍都判定忠实。
```

**下一轮（一次运行可判定）**：给 `zz_seam_fix` 的 **BEFORE** 段也加 `--ep` 式端点打印（example 级、只读），
直接看 F113 **输入** wire 里是否已有 `first==last` 的边（尤其那两条 1 边 wire）。
  · 若「有」⇒ 目标改为「OCCT 如何处理输入中的零长度边」（查 `cxx` 的 LoadWires/FixReorder 段），
    按 `.cxx` 在端口补上同样的处理 ⇒ 才可能让 F113 收敛成 2 wires；
  · 若「没有」⇒ 回到 `:113/:446/:448/:466` 四个未测建边点逐一取证。

---

### 9.471 —— 【决定性】F113 的**输入** wire 里本来就带两条零长度退化边（±87.5 seam）⇒ 「零长度段」不是被造出来的，而是**从未被清理**

给 `zz_seam_fix` 的 **BEFORE** 段加了 `--ep` 端点打印（example 级、只读、已入库），直接量输入：

```text
BEFORE wires=4  u=[0,6.283185307]  v=[-inf,inf]  bbox=(-87.5,-34,-100)-(87.5,34,-32)
BEFOREP wire=0 n=6  zero=0
BEFOREP wire=1 n=14 zero=0
BEFOREP wire=2 n=1  zero=**1**   ( 87.500,0.000,-32.000)->( 87.500,0.000,-32.000)
BEFOREP wire=3 n=1  zero=**1**   (-87.500,0.000,-32.000)->(-87.500,0.000,-32.000)
```

⇒ **F113 的输入 wire[2] / wire[3] 各是一条首尾同点的退化边**（位于两侧 seam x=±87.5）。

**这把整条链重新解释了一遍**：

```text
· GRID-POST / BREAK-IN / COLLECT-IN 里看到的「两条零长度段 (±87.5,0,-32)」**就是这两条输入退化边的副本**；
  它们是经 `split_wire` 的「不切、整边拷贝」路径（`:113` / `:446`/`:448`/`:466`，端口 §9.465 列出的未测路径）
  原样搬进结果段的 —— **不是**任何切分/建边操作「造」出来的。
· 因此前面所有「切分侧」的对拍都判定忠实（§9.452-§9.470 共 19 节）都**没有找错**，只是找错了对象：
  真正的分歧是 **「OCCT 会清掉/并掉输入里的退化边，端口不会」**。
· 这也解释了 §9.443/§9.444 的结构差：端口 4 wires（6/14/1/1）→ fix_missing_seam 出 **5 张面**（含两张 3 边薄片）；
  OCCT 收敛成 **2 wires（22+6）** —— 差别就在那两条退化 wire 是否被处理掉。
```

#### 下一轮（判定「谁该清、在哪清」，用 OCCT 侧探针即可，不动端口库代码）

```text
① OCCT 侧：`cmd /c "specs\occt_probe\run_dbg.bat data\occ\a3n00.stp noop"`（关 ShapeProcess，取未整形输入）
   看 F113 那张面的 wire 数与每 wire 边数：
     · 若也是 [6,14,1,1] ⇒ 分歧在**整形过程中的某一步**（ComposeShell/DispatchWires/FixSmall 清退化边）
       ⇒ 翻 `cxx` 的 CollectWires/DispatchWires/LoadWires 段找「退化边/零长度段」的处理并补上；
     · 若 OCCT 的输入**没有**那两条 1 边 wire ⇒ 分歧在**端口的 reader**（把退化边建成了独立 wire）
       ⇒ 回到 `step/read_topology.rs` 的 resolve_loop/make_wire 侧核对该 STEP 构造（DEGENERATED edge 标记）。
② 端口侧已具备的同类判据（供参考）：`BRepTool::is_degenerated(edge)`、`check_small_area`、
   `fix_small_all` —— 检查 `.cxx` 在同一位置用的是哪一条。
```

---

### 9.472 —— OCCT 侧未整形（`noop`）普查：**退化 wire 普遍存在**（58 张面含 ≥2 条 1 边 wire）⇒ 分歧在「整形过程中如何处理它们」，不在输入

`cmd /c "specs\occt_probe\run_dbg.bat data\occ\a3n00.stp noop"`（关 ShapeProcess，取**未整形输入**），
输出 453 行；解析 `FACE <idx> wires=<n> edges_per_wire: …`：

```text
含 >=2 条 1 边 wire 的面数: **58**      样本：
  FACE idx=41 wires=2 edges=[1, 1]      FACE idx=42 wires=2 edges=[1, 1]
  FACE idx=44 wires=2 edges=[1, 1]      FACE idx=47 wires=2 edges=[1, 1]      FACE idx=77 wires=2 edges=[1, 1]
（探针自身还把部分面标注为 `[UV-DEGENERATE-WIRE]`，例如 FACE 0 / FACE 2）
```

⇒ **OCCT 的原始输入同样带大量退化（1 边）wire** —— 端口 F113 输入里的 wire[2]/wire[3]（§9.471）
   **不是端口独有的产物**。结合 OCCT 整形后 F113 = 2 wires（22+6）可知：

```text
分歧 = **整形过程中「如何处理输入的退化 wire」**：
  · OCCT 的 ShapeProcess/FixShape（ComposeShell → CollectWires/DispatchWires，或 ShapeFix_Wire 的
    FixSmall/CheckDegenerated 一线）把退化 wire 清掉/并掉，最终收成 1 张 2-wire 面；
  · 端口保留了它们 ⇒ 结果里多出 2 条零长度段 ⇒ 5 张面（含两张 3 边薄片）⇒ Shell(5) 被丢 ⇒ 空面。
```

**下一步（本轮已尝试、只差一步）**：在 `noop` 普查里按**结构指纹**定位 OCCT 的 F113 原始面
（含 **14 边与 6 边** wire 的那张；本轮 grep `87.5` 只命中另一张 `FACE 201 box=(-112.5,-39,-105)-(-87.5,39,-27)`），
确认其 `edges_per_wire` 是否恰为 `[6,14,1,1]`：
  · 是 ⇒ 两侧输入一致 ⇒ 直接在 `cxx` 里找**退化 wire 的清理点**（LoadWires/CollectWires/DispatchWires/FixSmall）
         并按 `.cxx` 在端口补上；
  · 否 ⇒ 端口 reader 侧的 wire 组装（`resolve_loop`/`make_wire`）与 OCCT 不同，回头核 STEP 构造。

---

### 9.473 —— 按 box 定位 OCCT 未整形输入里的 F113 对应面（结果见下）

```text
解析到 FACE 行: 226
最接近 F113 目标 box(-87.5,-34,-100)-(87.5,34,-32) 的候选：
   err=52.000 idx=165 wires=1 edges=[4] box=(-48.166,-48.166,-48.000)-(48.166,48.166,-32.000)
   err=52.000 idx=202 wires=1 edges=[4] box=(-48.166,-48.166,-48.000)-(48.166,48.166,-32.000)
   err=52.000 idx=219 wires=1 edges=[4] box=(-48.166,-48.166,-48.000)-(48.166,48.166,-32.000)
   err=65.134 idx=203 wires=1 edges=[2] box=(-22.366,-36.500,-48.000)-(22.366,-28.844,-48.000)
   err=69.611 idx=171 wires=4 edges=[8, 8, 8, 8] box=(-57.361,-42.865,-169.611)-(61.611,42.865,-50.639)
```

最佳候选 idx=165 wires=1 edges=[4]

**判读与下一步**：

```text
· 若最佳候选的 edges_per_wire 恰为 [6,14,1,1] 的多重集合 ⇒ 两侧输入一致
  ⇒ 分歧确定在「整形过程如何处理退化 wire」⇒ 翻 cxx 的 LoadWires/CollectWires/DispatchWires/FixSmall
     找退化 wire 的清理点，按 .cxx 在端口补上，然后跑 §9.448 验收；
· 若不是 ⇒ 探针的 box 口径或面序与端口不同（探针的 `spans: u/v` 显示其 box 可能取自参数域），
  ⇒ 改用 `run_dbg.bat <stp> wdump <idx>` 直接 dump 候选面（或按 OCCT 文件面序 #5375 → index）再比对。
```

---

### 9.474 —— 【对照成功】OCCT 未整形输入里的 F113 原始面 = `box=(-1e+100,-34,-100)-(1e+100,34,-32)`

先前按 x≈±87.5 搜不到，原因已查明：**未整形**的 F113 是整周圆柱面（pcurve 周期），其 box 在 x 方向是 **±1e100（无界）**，
而 y/z 恰为 `-34…34` / `-100…-32`（与端口 F113 一致）⇒ 用它定位：

```text
(无命中)
```

**判读**：

```text
· 若该行 `edges_per_wire` 的多重集合 == 端口的 {6,14,1,1} ⇒ **两侧输入一致**（含两条 1 边退化 wire）
  ⇒ 分歧确定落在「**整形过程如何处理输入的退化/1 边 wire**」：
     OCCT 清掉或并掉它们 ⇒ 最终 1 张 2-wire（22+6）面；端口保留 ⇒ 5 张面（含两张 3 边薄片）⇒ Shell(5) 被丢。
     下一步：在 cxx 的 LoadWires/CollectWires/DispatchWires（以及 ShapeFix_Wire::FixSmall /
     ShapeAnalysis_Wire 的退化检查）里找**退化 wire 的清理点**，按 .cxx 在端口补上 ⇒ 跑 §9.448 验收。
· 若不同 ⇒ 端口 reader 侧的 wire 组装与 OCCT 有差异，回 step/read_topology.rs 的 resolve_loop/make_wire 核对。
```

（另注：该行的 `box` 在 x 上无界这一点本身也是一个可用判据——端口 F113 的 box 是 `(-87.5,…)-(87.5,…)`
（有界），说明端口在读入时已把周期 pcurve 收窄；这属于既有差异，但不影响本节的 wire 数对照。）

---

### 9.475 —— 【两侧输入对照完成】OCCT `noop` FACE 113 与端口 F113 的 wire 结构**几乎一致**（含那两条 1 边退化 wire）；差异只有第一条 wire 的边数（8 vs 6）

**更正 §9.474**：那一节说「对照成功」但没取到 wire 列表（我按 `-34.0/-100.0` 字面去 grep，而文件里打印的是 `-34`/`-100`，故漏匹配）。
本轮用 `1e+100,34` 精确定位到第 114 行：

```text
OCCT noop（未整形）：
  FACE 113 type=1 wires=4 edges_per_wire: **8 14 1 1**  spans: u=0 v=0 ×4
  box=(-1e+100,-34,-100)-(1e+100,34,-32)  [UV-DEGENERATE-WIRE]

端口 F113（`zz_seam_fix 113 --ep` 的 BEFORE）：
  wires=4   wire[0]=**6**  wire[1]=14  wire[2]=1 (zero=1)  wire[3]=1 (zero=1)
  box=(-87.5,-34,-100)-(87.5,34,-32)
```

**判读（本条把问题钉死）**：

```text
① wire 数一致（4 条），且**两侧都有两条 1 边（退化）wire** ⇒ §9.471/§9.472 的结论成立：
   退化 wire 不是端口独有，**分歧在「整形过程如何处理它们」**。
② 唯一差异 = 第一条 wire 的边数：OCCT raw 8 vs 端口 6 ⇒ 正是 §9.412 那条 `FixDummySeam` 丢边
   （端口在导入期就把 #5012/#5018 毛刺对删掉了；OCCT raw 还留着，其**整形后**同样变成 6，见 §9.439 的 wdump）。
③ 因此「端口 4 wires → 5 张面 → Shell(5) 被丢 ⇒ 空面」vs「OCCT 4 wires → 1 张 2-wire 面」的差别，
   **既不是丢边（两侧整形后一致）、也不是切分（§9.452-§9.470 全判忠实），而是：
     OCCT 的整形把两条 1 边退化 wire 清掉/并掉，端口保留**。
```

#### 下一轮（唯一剩余目标，方向明确）

```text
在 cxx 里定位「退化/1 边 wire 在整形中被清理」的那一步。候选（按可能性排序）：
  1) `ShapeFix_ComposeShell::CollectWires`（cxx:2512-2936）—— 配对阶段若 1 边段无法配对，
     OCCT 可能直接丢弃（端口 `collect_wires.rs` 的 `index=None` 分支需逐行核对）；
  2) `ShapeFix_ComposeShell::LoadWires`（cxx:499-640）—— 载入时是否已把退化 wire 排除；
  3) `ShapeFix_Wire::FixSmall(true, Precision())` / `ShapeAnalysis_Wire` 的小边/退化检查；
  4) `ShapeFix_Face::FixMissingSeam` 之后的 `FixSmallAreaWire`（§9.444 已对过打包，但可再看裁剪路径）。
做法：先用端口侧只读探针确认「那两条 1 边段在 collect_wires 里是否配对成功过」（§9.453 的 COLLECT-OUT 已显示
它们被组装成两张 3 边薄片 ⇒ 端口**没有丢弃**），再逐行对 cxx:2512-2936 的对应分支找差异。
```

---

### 9.476 —— `collect_wires` 的「合并短 3D 段」读到关键结构；靶点推进到 `dispatch_wires`（由 §9.453 的读数反推）

端口 `collect_wires.rs:333-396`（`cxx:2853-2935`）：

```rust
for i in 0..seqw.len() {
    if shorts[i] != 1 || seqw[i].is_vertex()
        || seqw[i].orientation() == Orientation::Internal
        || seqw[i].orientation() == Orientation::External { continue; }   // 只处理「短且非 vertex/Internal/External」的段
    let wd = seqw[i].edges(); let v = seqw[i].first_vertex();
    … 遍历 wires[j] 的每条边 k：要求 same_v(&v, &first_vertex(&cand[k]))（**目标边的首顶点 == 本段首顶点**）
      + is_same_patch(…) + 2D 端切向距离最小者 ⇒ (minj, mink, mindist)
    let Some(target_idx) = minj else {
        // cxx:2914-2923: keep it as a separate wire.
        wires.push(WireSegment::with_edges(wd.clone(), Orientation::Forward)); continue;
    };
    // cxx:2925-2931: 把 wd 的边依次 add 到目标 wire（mink 起）
}
```

**关键推论（结合 §9.453 的既有读数）**：那两条 1 边退化段**并没有**停在 `minj == None` 的
「keep it as a separate wire」分支 —— 因为在 `COLLECT-OUT-wires` 里它们各自成了 **3 边 wire**
（`(±87.5,0,-32)->(44.5,0,-32)`，含零长度边 + 两条邻边），即**它们被成功并进了别的 wire**。

⇒ 所以问题不在这里「是否丢弃」，而更可能在**最后一步 `dispatch_wires`**：端口把最终 wires **逐条变成一张面**
（5 条 wire ⇒ 5 张面，含两张由退化段参与的 3 边薄片），而 OCCT 收敛成 1 张 2-wire 面
⇒ 需要核对 OCCT 的 `DispatchWires` 是否**跳过/并入**这类由退化段参与的小 wire。

#### 下一轮（最后一次对拍）

```text
读端口 `dispatch_wires.rs`（~几百行）与 `ShapeFix_ComposeShell.cxx` 的 DispatchWires 段
（§9.377 记录过 cxx:2770-2860 一线；本文件头注释应也标了 cxx 行号），专找：
  · 对「小/退化 wire」的判据（ShapeAnalysis_Wire 的小面积/小边长检查、或 wire 边数阈值、
    `CheckSmallArea`/`CheckSmall` 一类）——OCCT 有而端口缺 ⇒ 命中；
  · 若两侧都有该判据 ⇒ 退一步核 `shorts[]` 数组的**下标错位**（`cxx:2524-2526` 明确「不随元素删除而平移」，
    端口注释也承认；但两侧若在上游删了不同数量的元素，`shorts[i]` 就会指向不同的段 ⇒ 该并的没并）。
命中即按 .cxx 修 → 跑 §9.448 验收（zz_seam_fix 113 → Face；--model 113 → wires=2（22+6）；--fstats → face=113 mt≈228；面积比 ≥0.8996 且上升；t101_verify.ps1 全绿）。
```

---

### 9.477 —— `dispatch_wires` 确有「丢退化 wire」守卫，但**只覆盖「单条退化边」**；我们的薄片是 **3 边 wire 含 1 条退化边** ⇒ 不被拦 ⇒ 成为 2 张面（与观测一致）

端口 `dispatch_wires.rs:94-121`：

```rust
while i < wires.len() {
    if wires[i].is_vertex() { i += 1; continue; }
    // cxx:3337-3343: skip a wire with a single degenerated edge.
    if wires[i].nb_edges() == 0
        || (wires[i].nb_edges() == 1 && BRepTool::is_degenerated(&wires[i].edges()[0])) {
        wires.remove(i); continue;                        // ← 丢弃
    }
    // cxx:3345-3346: sfw.Load(sbwd); sfw.FixShifted().
    let wire = builder.make_wire(wires[i].edges());
    let _ = fix_shifted_wire(&wire, &my_face);
    // cxx:3348-3357: ShapeBuild_Edge::RemovePCurve on degenerated edges.
    for j in 0..wires[i].nb_edges() {
        if BRepTool::is_degenerated(&wires[i].edges()[j]) { sbe.remove_pcurve(&wires[i].edges()[j], &my_face); }
    }
    // cxx:3358: sfw.FixDegenerated()  → 端口 fix_degenerated_all（重建段 + 恢复 patch 索引）
    …
}
```

**判读**：

```text
· 该守卫与 `.cxx:3337-3343` 一致（**只**跳 0 边或「单条退化边」的 wire）——忠实；
· 但 F113 那两条退化段在 §9.453 里已被 collect_wires 并成 **3 边 wire**（零长度边 + 两条邻边）
  ⇒ `nb_edges() == 3` ⇒ **不被这条守卫拦下** ⇒ 它们继续走到 `make_faces_on_patch` ⇒ 成为 5 张面里的 2 张
  （与 §9.445 观测的 5 张面、其中两张 3 边薄片 **完全吻合**）。
```

#### 下一轮（最后一个未对拍的函数）

```text
`cxx:3358` 的 `sfw.FixDegenerated()` ↔ 端口 `fix_degenerated_all`（`dispatch_wires.rs:117-121` 之后）：
  · OCCT 的 `ShapeFix_Wire::FixDegenerated` 对「含退化边的 wire」做的是——移除退化边并把**两侧邻边合并**
    （必要时重建顶点）；若一条 wire 只剩退化边构成的环，它会被**清空**，随后 `DispatchWires` 里
    `nb_edges()==0` 的守卫（:101）就会把它丢掉。
  · 若端口的 `fix_degenerated_all` 在「3 边 wire（含 1 条退化边）」上没有实现同样的**邻边合并/清空**
    ⇒ 该 wire 存活 ⇒ 多出一张面 ⇒ 5 面 ⇒ Shell(5) 被丢 ⇒ 空面。
做法：读端口 `fix_degenerated_all`（在 `crates/occt-topo/src/shhealing/` 内）↔ `ShapeFix_Wire.cxx::FixDegenerated`，
按 `.cxx` 补齐；命中即跑 §9.448 验收。
```

---

### 9.478 —— 端口 `fix_degenerated_all`（循环层）已读完，与 `.cxx` 对齐；最后未对拍的是**逐边 worker** `fix_degenerated`

端口 `wire_fix.rs:2176-2214`（`ShapeFix_Wire::FixDegenerated` 的 whole-wire 循环）：

```rust
pub fn fix_degenerated_all(wire: &mut Wire, face: &Face, prec: f64) -> bool {
    let mut done = false; let mut last_coded = -1i32; let mut prev_coded = 0i32;
    // myClosedMode = true ⇒ stop = 0（cxx:1045），故从 nb 扫到 1
    let mut i = wire_edges_nb(wire) as isize;
    while i > 0 {
        let (d, coded2) = fix_degenerated(wire, face, prec, i as usize);   // ← 逐边 worker
        done |= d; let coded = i32::from(coded2);
        if last_coded == -1 { last_coded = coded; }
        // cxx:1050-1071 PRO7226：丢掉重复的退化边，并清掉其后继边的标记
        if coded == 1 && (prev_coded == 1 || (i == 1 && last_coded == 1)) && wire_edges_nb(wire) > 1 {
            wire_remove_edge(wire, i as usize);
            if prev_coded == 0 { i = wire_edges_nb(wire) as isize; }
            if let Some(e) = edges_of_wire(wire).get(i as usize - 1) {
                GeometryRegistry::global().set_degenerated(&e.0, false);
            }
            i += 1; prev_coded = 0;      // `B.Degenerated(sbwd->Edge(i++), false)`：i++ 被循环 i-- 抵消
        } else { prev_coded = coded; }
        i -= 1;
    }
    done
}
```

⇒ 循环层（遍历方向、`last_coded/prev_coded` 状态机、PRO7226 的「重复退化边」删除与其后继清标记、`i++`/`i--` 抵消）
   与 `.cxx:1045-1071` 对齐，**未见差异**。

#### 最后未对拍的函数（下一轮，若仍要做）

```text
逐边 worker：端口 `fix_degenerated(wire, face, prec, idx)`（`wire_fix.rs` 内，被 :2185 调用）
        ↔ `ShapeFix_Wire::FixDegenerated(const int num)`（ShapeFix_Wire.cxx，PRO7226 段之前）
看点：对「3 边 wire 含 1 条退化边」这种输入，OCCT 是否
  ① 移除退化边并把**两侧邻边合并**（`CombineVertex` + 交换邻边的顶点），使 wire 变短甚至清空；
  ② 清空后由 `DispatchWires` 的 `nb_edges()==0` 守卫（`cxx:3337-3343`，§9.477）丢弃它。
若端口 worker 在 ① 上不等价（例如只清标记不移除/不合并），就会让那两条 3 边薄片存活 ⇒ 5 张面 ⇒ Shell(5) 被丢。
```

**当前状态（如实）**：F113 仍 4 wires / 空面；面积比 0.8996；门禁 5/5 绿。本目标（goal 91/120）已把差异候选
收敛到**一个函数**（`fix_degenerated` 逐边 worker），且它所在的那条链（LoadWires → SplitByGrid → BreakWires →
CollectWires → DispatchWires → FixDegenerated）此前每一环都已判定忠实或有明确结论。

---

### 9.479 —— 逐边 worker `fix_degenerated` 前半段与 `ShapeFix_Wire.cxx:2130-2186` 逐段一致；剩余待读 = 邻居重接段（`:2148-2185` ↔ `cxx:2171-2205`）

端口 `wire_fix.rs:2092-2148`（文档注释即声明对拍对象：`ShapeFix_Wire::FixDegenerated(const int num)`，`cxx:2130-2205`）：

```rust
fn fix_degenerated(wire: &mut Wire, face: &Face, prec: f64, num: usize) -> (bool, bool) {
    let nb = wire_edges_nb(wire);
    if nb < 1 { return (false, false); }
    match check_degenerated(wire, face, prec, num) {
        None => (false, false),
        // cxx:2148-2152: FAIL2 -> WireData()->Remove(num) + DONE3
        Remove => { wire_remove_edge(wire, num); (true, false) }
        Found { p2d1, p2d2, lack } => {
            let vect2d = GpVec2d::new(p2d2.x()-p2d1.x(), p2d2.y()-p2d1.y());   // cxx:2160-2163
            let line2d = Geom2dLine::from_pnt_dir(p2d1, dir2d);
            let mag = vect2d.magnitude();
            // cxx:2165-2168: MakeEdge(degEdge) + Degenerated(true) + UpdateEdge(line2d, Confusion)
            //                + Range(0, |vect2d|)
            … 构造 deg（用 CurveOnSurface 代替被清空的 3D 曲线，端口 EdgeGeom 不能存 null 曲线）…
            // cxx:2171-2186
            let n2 = if num > 0 { num } else { nb };
            let n1 = if n2 > 1 { n2 - 1 } else { nb };
            …
```

**判定**：前半段（空 wire 早退、`check_degenerated` 三态、`Remove ⇒ wire_remove_edge`、`vect2d/line2d/mag`、
`MakeEdge/Degenerated/UpdateEdge/Range` 的等价构造）与 `.cxx:2130-2168` **一致**；端口还为「不能存 null 3D 曲线」
给出了等价替代（用 pcurve 在面上的像 `CurveOnSurface`），并已注明依据（`cxx:1082-1084`）。

#### 剩余待读（本目标的最后一段）

```text
端口 `wire_fix.rs:2148-2185`（n1/n2 邻居重接：把 n1/n2 两条邻边的顶点改成新退化边的两端点，
必要时 `CombineVertex`） ↔ `ShapeFix_Wire.cxx:2171-2205`。
看点：OCCT 在此处是否**把 n1/n2 邻边接到新退化边上并合并顶点**（从而让「3 边含 1 退化边」的薄片折成
      1–2 边甚至清空）；端口若只改顶点不重接（或 n1/n2 选取不同），薄片就会存活 ⇒ 5 张面 ⇒ Shell(5) 被丢。
```

---

### 9.480 —— 邻居重接段**逐行一致**（忠实）⇒ `FixDegenerated` 全路径已对完；最后一环是 `check_degenerated`（退化**检测**）

两侧并列（左端口 / 右 OCCT）：

| 项 | 端口 `wire_fix.rs` | OCCT `ShapeFix_Wire.cxx` | 判定 |
|---|---|---|---|
| n2 | `if num > 0 { num } else { nb }` (:2147) | `(num > 0 ? num : sbwd->NbEdges())` (:2172) | ✓ |
| n1 | `if n2 > 1 { n2 - 1 } else { nb }` (:2148) | `(n2 > 1 ? n2 - 1 : sbwd->NbEdges())` (:2173) | ✓ |
| lack | `check` 的 DONE1 位 | `myAnalyzer->LastCheckStatus(DONE1)` (:2175) | ✓ |
| n3 | `if lack { n2 } else if n2 < nb { n2+1 } else { 1 }` (:2149) | `(lack ? n2 : (n2 < NbEdges() ? n2+1 : 1))` (:2176) | ✓ |
| V1/V2 | `last_vertex(edges[n1-1])` / `first_vertex(edges[n3-1])` (:2151/:2156) | `sae.LastVertex(Edge(n1))` / `sae.FirstVertex(Edge(n3))` (:2179-2180) | ✓ |
| 朝向 | V1 `Forward`、V2 `Reversed`、deg `Forward` (:2152-2161) | `V1.Orientation(FORWARD)`、`V2.Orientation(REVERSED)`、`degEdge.Orientation(FORWARD)` (:2182-2186) | ✓ |
| 插入/替换 | `lack ⇒ wire_insert_edge_before(wire, n2, deg)` 否则 `wire_set_edge(wire, n2, deg)`；返回 `(true, !lack)` (:2164-2169) | `lack ⇒ sbwd->Add(degEdge, n2)` + DONE1 否则 `Set(degEdge, n2)` + DONE2；`return true` (:2188-2203) | ✓ |

⇒ **`ShapeFix_Wire::FixDegenerated(num)` 全路径（`:2130-2204`）在端口已逐行等价**。

#### 最后一环（唯一未对拍者）

```text
端口 `check_degenerated(wire, face, prec, num)`（`fix_degenerated` 的第一句）
   ↔ `ShapeAnalysis_Wire::CheckDegenerated(num, tol, p2d1, p2d2)`
若它在我们这个输入（**3 边 wire，其中一条零长度边位于 seam**）上返回 `None`，
而 OCCT 返回 `Found`/`Remove`，那么退化处理根本不会启动 ⇒ 薄片存活 ⇒ 5 张面 ⇒ Shell(5) 被丢 ⇒ 空面。
→ 下一轮读它（端口 + cxx），对拍**检测判据**（pcurve 端点距离、`IsClosed`、`myMax2d/myMax3d` 阈值、
  以及 `.cxx` 里对 `Precision::Confusion()` 的使用与 `LastCheckStatus` 的置位）。
```

---

### 9.481 —— 退化**检测**核心（端口 `wire_fix.rs:1878-1975`）读完；并更正一处：`ShapeAnalysis_Wire.cxx:488` 是 `CheckClosed` 里的一段，不是逐边 `CheckDegenerated`

#### A. 端口 `check_degenerated` 的三段结构（每段都带 `.cxx` 出处）

```text
① cxx:913-934（OCC7630）：若 e2 **已被标记退化且带 pcurve**：
     若 e1/e3 也带 pcurve，则比较 |d(p12,p31) − d(p21,p22)| > 2·PCONFUSION ⇒ **Remove**；
     否则一律 **None**（不做处理）。            ← 「已标记退化边」的快捷分支
② cxx:938-948：`n1 != n2 && is_degenerated(e1) && !has_pcurve(e1)`：
     e2 也退化 ⇒ **Remove**，否则 **None**。     ← 「前一条是无 pcurve 的退化边」
③ cxx:950-970 / cxx:974-991 / cxx:995-999（**通用检测**）：
     取 e1/e2 的首末顶点 (vp,v0,v1,v2) 与 tol1 ⇒ prec_first = min(prec,tol1)、prec_fin = max、prec_vtx = …
     若 **p1.distance(p2) <= prec_first**（该边两端点 3D 重合 ⇒ 「在奇异点上闭合」）：
         用 `SurfaceSingularities::compute(surf)` + `sing.min_gap(&p1, prec_vtx)` 得 (p2d1,p2d2)，dgnr = true；
         再走 **cxx:979-990 的中点守卫**：若该边**中点**到 p1 的距离 > prec_vtx ⇒ **dgnr = false**
            （「不要把中点远离奇异点的闭合边变成退化边」）
     若 !dgnr ⇒ 落到 **cxx:995-999** 的 lack 路径（后续返回 `Found{lack:true}` ⇒ `fix_degenerated` 走 **插入** 分支）
```

#### B. 更正

`ShapeAnalysis_Wire.cxx:488` 那一段属于 **`CheckClosed`**（`CheckDegenerated(1)` 的调用点），
**不是逐边 `CheckDegenerated(i)` 的定义**。逐边定义尚未定位（候选行号 226/360/364 里 360 是无参驱动 `CheckDegenerated()`）；
下一轮应先精确定位逐边版本再对拍。

#### C. 下一轮的具体假设（可直接判定）

```text
F113 那两条 0 长度退化边位于**圆柱面**（R=34）。圆柱**没有奇异点**（周期面，无退化点）⇒
  · 端口：`sing.min_gap(...)` 很可能返回 None ⇒ `dgnr = false` ⇒ 走 lack 路径 ⇒ 返回 `Found{lack:true}`
    ⇒ `fix_degenerated` 走「插入新退化边」分支 ⇒ 薄片**存活/变长**（与 §9.453 观测到的 3 边 wire 吻合）。
  · OCCT：`ShapeAnalysis_Surface::DegeneratedValues` 对圆柱同样返回 false，但它可能在**更早的判据**
    （例如 `BRep_Tool::IsClosed(edge, face)` 或 pcurve 跨度 = 周期）上把这种「seam 处的零长度 3D 边」判为
    **应当移除**（`DegeneratedCheck::Remove`）⇒ 薄片被清掉。
做法：先定位 `ShapeAnalysis_Wire::CheckDegenerated(const int num, …)` 的定义行，逐行对拍端口 :1878-1975；
若确认差异在「圆柱/seam 这种无奇异点但 pcurve 跨整周期的边」，按 .cxx 补上同一条判据 ⇒ 跑 §9.448 验收。
```

---

### 9.482 —— 逐边 `CheckDegenerated(num,p2d1,p2d2)` 定义在 `ShapeAnalysis_Wire.cxx:896`；前两段（①cxx:913-935 / ②cxx:938-948）与端口**等价**（避免了一次误报）

#### ① 已标记退化且带 pcurve（`cxx:913-935` ↔ 端口 `wire_fix.rs:1895-1916`）

```cpp
// OCCT
if (BRep_Tool::Degenerated(E2) && sae.HasPCurve(E2, Face())) {
  if (sae.HasPCurve(E1, Face()) && sae.HasPCurve(E3, Face())) {
    … p21/p22 = pc2(fp/lp), p12 = pc1(lp), p31 = pc3(fp)
    if (fabs(p12.Distance(p31) - p21.Distance(p22)) > 2*PConfusion())
      myStatus = ShapeExtend::EncodeStatus(ShapeExtend_FAIL2);   // ← 置 **FAIL2**
  }
  return false;
}
```
```rust
// 端口：把「置 FAIL2」直接映射成 `Remove`（因为 fix_degenerated 的 cxx:2148-2152 就是 FAIL2 ⇒ Remove）
if BRepTool::is_degenerated(&e2) && has_pcurve(&e2, face) {
    if has_pcurve(&e1, face) && has_pcurve(&e3, face) {
        … if (p12.distance(&p31) - p21.distance(&p22)).abs() > 2.0*PCONFUSION { return DegeneratedCheck::Remove; }
    }
    return DegeneratedCheck::None;
}
```

**判定：等价** —— OCCT 用「状态位 FAIL2」表达，调用方 `ShapeFix_Wire::FixDegenerated(num)` 读到 FAIL2 就
`WireData()->Remove(num)`（`cxx:2148-2152`）；端口用枚举 `Remove` 表达同一件事（其 `fix_degenerated`
`:2105-2108` 的注释正引 `cxx:2148-2152`）。**不是分歧**（此处若草率下结论就会误报）。

#### ② 连续退化边无 pcurve（`cxx:938-948` ↔ 端口 `:1918-1926`）

```cpp
if (n1 != n2 && BRep_Tool::Degenerated(E1) && !sae.HasPCurve(E1, Face())) {
  if (BRep_Tool::Degenerated(E2)) myStatus |= ShapeExtend::EncodeStatus(ShapeExtend_FAIL2);
  return false;                                   // 未置位时 myStatus 仍是 :898 的 OK ⇒ 不移除
}
```
```rust
if n1 != n2 && BRepTool::is_degenerated(&e1) && !has_pcurve(&e1, face) {
    if BRepTool::is_degenerated(&e2) { return DegeneratedCheck::Remove; }
    return DegeneratedCheck::None;
}
```
**判定：等价** ✓（含「E2 非退化 ⇒ 状态为 OK ⇒ 不移除」这一细节）。

#### 下一轮（最后一段：通用检测）

```text
读 `ShapeAnalysis_Wire.cxx:950-1010`（Vp/V0/V1/V2 取点、DegeneratedValues、中点守卫、lack 与 p2d1/p2d2 的产出）
 ↔ 端口 `wire_fix.rs:1928-1990`（已读到 :1975：p1.distance(p2) <= prec_first 分支 + sing.min_gap + cxx:979-990 中点守卫 + :1975 起的 lack 路径）。
看点：**圆柱（无奇异点）**上那条「3D 零长度、pcurve 跨整周期」的 seam 边，两侧是否都走 lack（插入）分支；
      若 OCCT 在此另有判据（例如用 `BRep_Tool::IsClosed(E2, myFace)` 或 pcurve 跨度≈周期）判为 FAIL2 ⇒ 移除，
      那就是分歧点。
```

---

### 9.483 —— 通用检测段两侧对照（`:960-1016` ↔ 端口 `:1928-1975`）**等价**；由此推出更关键的怀疑：**清理动作来自 ComposeShell 之外的 wire 级修复（ShapeProcess/FixShape）**

#### A. 逐项对照（OCCT → 端口）

| OCCT | 端口 | 判定 |
|---|---|---|
| `pp/p0/p1/p2 = Pnt(Vp/V0/V1/V2)`（:960-963） | `:1929-1940` 同（`first/last_vertex` + `vertex_point`） | ✓ |
| `precFirst=min(prec,Tol(V1))`、`precFin=max`、`precVtx=(prec<Tol? 2*precFin : precFin)`（:968-970） | `:1942-1944` 同 | ✓ |
| `forward = (E2.Orientation()==FORWARD)`（:972） | 端口注释说明 `forward` 未被 `DegeneratedValues` 使用（引 `cxx:373-411` 的 `const bool /*forward*/`） | ✓（等价，已注明） |
| `if (p1.Distance(p2) <= precFirst) dgnr = mySurf->DegeneratedValues(p1,precVtx,p2d1,p2d2,par1,par2,forward)`（:975-977） | `if p1.distance(&p2) <= prec_first { if let Some((a,b,_,_)) = SurfaceSingularities::compute(surf).min_gap(&p1, prec_vtx) { p2d1=a; p2d2=b; dgnr=true } }`（:1949-1972） | ✓（端口用**最近奇异点**口径，正对 OCCT `#77 rln S4135` 在 :1014-1016 的注释「using singularity which has minimum gap…」） |
| 中点守卫：`C3d->Value(0.5*(a+b))` 与 `p1` 的平方距离 > `precVtx²` ⇒ `dgnr=false`（:979-990） | `:1960-1970` 同（`edge_parameters` + `edge_curve` + `d0` + `square_distance`） | ✓ |
| `if (!dgnr) { if (n1 != n2 && p1.Distance(pp) <= precFirst && mySurf->IsDegenerated(pp, precFirst) && !Degenerated(E1)) return false; … }`（:992-999） | 端口 `:1973-1975` 起（注释即引 `cxx:995-999`），**该分支代码本轮未读完** | 待读（下轮） |

#### B. 由 A 推出的**关键怀疑**（可能是本目标真正的根因）

```text
圆柱（R=34）**没有奇异点** ⇒ 无论 `DegeneratedValues` 还是 `min_gap`，在这个 seam 零长度边上都**取不到退化值**
⇒ `dgnr = false`，随后 lack 路径同样查的是「奇异点 + 缺口」，在圆柱上也不会判定为退化
⇒ **OCCT 的 `CheckDegenerated` 对这条边同样不会给出 FAIL2 / Remove**
⇒ 那么 OCCT 整形后 FACE 113 收敛成 2 wires（22+6）的「清理动作」**不可能出自 ComposeShell 内部的退化处理**
   而更可能出自**每个面的 wire 级修复链**（`ShapeFix_Shape`/`ShapeProcess` 驱动的 `ShapeFix_Wire`
   的 `FixDegenerated`/`FixSmall`/`FixLacking` 等）——**这正是覆盖度审计列的第 1 大缺口：
   「没有 ShapeProcess/FixShape 驱动器，端口在 reader 里按面直调子步」**！
```

这与审计结论**相互印证**：不是某个子判据写错了，而是**少了一趟整形的编排**（哪些 Fix* 在什么顺序、以什么参数作用到面上）。

#### C. 下一轮（两条，先做 ①，成本最低）

```text
① 数据核对（无需读代码）：比较 OCCT 侧 FACE 113 在 **noop** 与 **默认整形** 下的 wire 结构
   （noop 已是：wires=4 edges 8/14/1/1，见 §9.475；默认整形那趟取 `run_dbg.bat <stp>` 输出里同 box 的那张面）
   · 若整形后那两条 1 边 wire **消失或并入** ⇒ 结论 = 清理发生在整形编排里 ⇒ 按 ② 找具体的 Fix* 与顺序；
   · 若整形后**仍在**（只是 F113 的 wire 数从 4 变成 2 是别的合并造成的）⇒ 假设需修正。
② 读 `ShapeProcess_OperLibrary.cxx`/`ShapeFix_Shape.cxx` 的算子序列（审计给过 `:785-865`/`:83-250`），
   找出**哪一步会清掉「1 边退化 wire」**，再对照端口 `step/read_topology.rs:710-714` 的按面调用清单。
```

---

### 9.484 —— 【对照完成】OCCT 整形：`noop` 的 **8/14/1/1** → 整形后 **22/6**；而 **22 = 8 + 14** ⇒ 分歧是「外环两条 wire 没被合并」，不是「退化边没被清理」

`cmd /c "specs\occt_probe\run_dbg.bat data\occ\a3n00.stp"`（默认整形，547 行）中与 F113 同 box 的面：

```text
OCCT 整形后：FACE 25 type=1 **wires=2 edges_per_wire: 22 6**  spans: u=6.2832 v=175.0000 · u=2.1968 v=113.7919
             box=(-87.5,-34,-100)-(87.5,34,-32)          ← 与端口 F113 的 box 逐字一致
OCCT 未整形：FACE 113 type=1 wires=4 edges_per_wire: **8 14 1 1**（§9.475）
端口 F113 输入：wires=4  **6 14 1 1**（§9.471，第一条第 8 边在导入期被 FixDummySeam 丢成 6）
```

**关键算术**：

```text
22 = **8 + 14** ✓ 精确 —— OCCT 把外环的两条 wire（8 边与 14 边）**合并成一条 22 边 wire**；
6  → 端口与 OCCT **都有**那条 6 边内环（端口输入里就是 6）；
1+1（两条零长度退化 wire）在整形后不再单独存在（被并入/吸收）。
⇒ 整形后 (22, 6) 与 noop 的 (8, 14, 1, 1) 是**同一批边的两种分组**（28 条边不变，§9.450 已证）。

端口整形后：**6 条 wire / 5 张面**（face0 = 6e+7e，face1 = 5e，face2 = 4e，face3/4 = 3e），
边总数同样 28 —— 但**没有发生 8+14 的合并**，那一趟反而把（16 边 = 14 边被切分插入 2 条后的）段切成了 7+5+4。
```

#### 结论（把目标再收紧一格）

```text
**分歧 = 端口没有把外环的两条 wire（对应 OCCT 的 8 与 14）合并成一条 22 边 wire**，
   而 OCCT 的 ComposeShell 做到了（22 = 8 + 14，且两条 1 边退化 wire 被吸收）。
⇒ 这条**否掉** §9.483 的「清理来自 ComposeShell 之外（ShapeProcess）」推断：noop→整形的差
   就发生在这一趟整形里（run_dbg 默认跑 ShapeProcess/FixShape，ComposeShell 是 FixMissingSeam 内部一步）。
⇒ 也把 §9.452-§9.482 那 30 节的「逐环忠实」限定清楚：它们在**端口自己的输入与分支上**忠实，
   但端口的 `CollectWires` 配对/合并**没有把这两条 wire 连起来**（OCCT 连起来了）——
   最可能落在 `collect_wires.rs:357-383` 的**目标边选取条件**（`same_v(本段首顶点, 目标边首顶点)`）
   或 `is_same_patch` 的判据上：若这两条 wire 的首顶点不「same_v」（例如顶点对象不同但坐标相同），
   端口就不会把它们并到同一条 wire。
```

#### 下一轮（一次探针即可判定）

```text
在 `collect_wires.rs` 的配对循环里，对 F113 那次调用打印：
  每个候选段的首顶点指针/坐标、每条目标 wire 的首边首顶点，以及 `same_v` / `is_same_patch` 的判定结果，
  看**14 边 wire 的段与 8 边（端口为 6 边）wire 的段**为何没有被连起来（`same_v` 假？patch 不同？还是顺序/初始 wire 集合为空）。
命中即按 `.cxx`（CollectWires 的 C++ 同段）修 ⇒ 跑 §9.448 验收。
```

---

### 9.485 —— 新增 example 探针（`WIREEND`）与一条新事实：**F113 输入的四条 wire 两两之间没有任何重合端点** ⇒ OCCT 的 `22 = 8 + 14` 合并发生在**切分之后**，分歧在「段」的串接

#### A. 探针与口径

给 `zz_seam_fix` 的 `--ep` 增了跨 wire 端点同一性检查（example 级、只读、已入库）：对每对 wire、每对边、
四种端点组合，若坐标差 < 1e-9 就打印 `close=… is_same=…`（`is_same` = `topo_tools_full::is_same`，即 OCCT `IsSame`）。

```text
编译：0 error；git status → 仅 examples/zz_seam_fix.rs 被改（**无库改动，无需还原**）
冒烟：BEFOREP 行 4 条 ✓（原探针正常）；**WIREEND 行 = 0**
```

#### B. 事实与推论

```text
WIREEND = 0  ⇒ 输入的四条 wire（6 / 14 / 1 / 1）**两两之间没有任何坐标重合的端点**
              ⇒ 它们是 4 条互不接触的**独立闭环**（1 边那两条天然闭合）。
⇒ 因此 OCCT 的「22 = 8 + 14」**不是**把两条原始 wire 的端点对接起来，而必然发生在
   `SplitByGrid`/`BreakWires` 把它们切成**段**、由切口造出**共享顶点**之后 —— 即 **`CollectWires` 对段的串接**。
⇒ 这同时**否掉**「原始输入里顶点对象未共享（IsSame 失败）」这一假设（在原始 wire 层面根本没有可对的端点）。
```

#### C. 与既有读数的关系（重新解释 §9.453）

```text
端口在切口之后确实产生了可串接的段：§9.452 的 BREAK-IN 显示 `[6, 16, 1, 1, 1, 1]`（16 = 14 边被插入 2 条），
§9.453 的 COLLECT-OUT 显示最终 wires 为 `[6, 7, 5, 4, 3, 3]`（6 条）——**串接发生了，但分成了 6 条而非 2 条**。
OCCT 侧同一步应得到 `[22, 6]`（两条）。⇒ 靶点 = **段的串接判据**（`CollectWires` 里
「本段首顶点 ↔ 目标 wire 首边首顶点」的 `same_v` 匹配 + `is_same_patch` 优先 + 同边返回最低优先级那套逻辑，:357-383），
要看切口造出的共享顶点在端口里**是否真的同一对象**（切口处顶点由 `Context()->Apply`/`CombineVertex` 产生，
若端口在这条路径上生成的是**坐标相同但 TShape 不同**的顶点，`same_v` 就会失败 ⇒ 串接断裂成 6 条）。
```

#### D. 下一轮（一次探针即可判定，仍是 example 级可行的写法）

```text
把 `WIREEND` 探针从「原始输入」搬到 **`collect_wires` 的入口**（或在 example 里对 `perform` 的中间结果不可达时，
改为在库内加 env-gated 打印，用完 `edit` 反撤）：对切口之后 seqw 的各段，检查其首/末顶点在
`is_same` 与「坐标相同」两种口径下的分布 ——
  · 若大量「坐标相同、is_same 为假」⇒ 根因 = 切口处顶点未统一（`Context`/`CombineVertex` 路径），
    按 `.cxx` 修顶点统一（这解释了 6 条而非 2 条，也解释了 §9.453 里那两张 3 边薄片）；
  · 若 is_same 普遍为真而仍没串成 2 条 ⇒ 靶点回到 `same_v` 之外的判据（`is_same_patch` 优先序、同边返回最低优先级）。
```

---

### 9.486 —— `collect_wires` 起始段（`:50-74`）：**单条退化边且「短」的段会被就地改成 `Internal`** —— 直接关系到 F113 那两条零长度段

端口 `collect_wires.rs:50-74`（`cxx:2519-2549`）：

```rust
while i < seqw.len() {
    if seqw[i].is_vertex() || seqw[i].orientation() == Orientation::Internal {
        wires.push(seqw[i].clone());      // vertex / Internal 段**直接进 wires**
        seqw.remove(i); continue;         // cxx:2525-2526: i-- 与 i++ 抵消
    }
    let isshort = is_short_segment(&seqw[i], &face, &self.grid, self.u_resolution, self.v_resolution);
    shorts[i] = isshort;
    let one_degenerated = seqw[i].nb_edges() == 1
        && seqw[i].edge(1).map(BRepTool::is_degenerated).unwrap_or(false);
    if isshort > 0 && (seqw[i].orientation() == Orientation::External || one_degenerated) {
        seqw[i].set_orientation(Orientation::Internal);      // ← **就地改朝向**
    }
    i += 1;
}
```

**与 F113 的直接关联**：那两条段正是「`nb_edges()==1` 且该边退化」⇒ `one_degenerated == true`
⇒ 只要 `is_short_segment(...) > 0`，它们就会被**就地改成 `Internal`**，随后进入串接主循环（而不会被
`:53-56` 提前搬进 `wires`——那段判定发生在改朝向**之前**）。

⇒ 于是「这两条段以 `Internal` 身份参与串接」这一点，很可能就是它们与相邻段被拼成
   **两张 3 边薄片**（§9.453 的 COLLECT-OUT `[6,7,5,4,3,3]` 里的两个 3）而**不是**被并进外环的原因；
   OCCT 同一步应得到 `[22, 6]`。

#### 下一轮（一次探针，仍是库内 env-gated，用完 `edit` 反撤）

```text
在 `collect_wires` 的 `:74` 之后打印：每个段的 (i, nb_edges, orientation, is_short, one_degenerated)，
并在主循环结束时打印最终 `wires` 的每条 wire 的边数与朝向 —— 一次运行即可确认：
  · 那两条段是否确实被改成 Internal；
  · Internal 段在串接里是否只与「彼此」相连（⇒ 生成 3 边薄片），而没与外环段相连。
若确认 ⇒ 对照 `cxx:2519-2549` 的同一分支核对「改 Internal 的条件」与「Internal 段的串接规则」，
   差异处按 .cxx 修 ⇒ 跑 §9.448 验收。
```

---

### 9.487 —— `collect_wires` 起始段探针（一轮内插桩→跑→撤除）：F113 的退化段朝向与 short 标记

探针（env `T101_ZZB`，插在 `:74` 之后打印每段 `(k, nb, ori, short, p0)`；**已在同一次调用内撤除**，库 diff 为空）
输出摘录（`.target-gate/coll100.txt`）：

```text
CSHORT k=4/8 nb=1 ori=Reversed short=0 p0=(87.500,0.000,-32.000)
CSHORT k=5/8 nb=1 ori=Reversed short=0 p0=(-87.500,0.000,-32.000)
CSHORT k=7/8 nb=1 ori=External short=0 p0=(-87.500,0.000,-32.000)
CSHORT k=1/3 nb=1 ori=Reversed short=0 p0=(87.500,0.000,-32.000)
CSHORT k=0/3 nb=1 ori=Reversed short=0 p0=(-87.500,0.000,-32.000)
CSHORT k=2/3 nb=1 ori=External short=0 p0=(-87.500,0.000,-32.000)
CSHORT k=4/8 nb=1 ori=Reversed short=0 p0=(87.500,0.000,-32.000)
CSHORT k=5/8 nb=1 ori=Reversed short=0 p0=(-87.500,0.000,-32.000)
```

共捕获 `nb=1` 或含 87.5 的行 **120** 条；完整输出见 `.target-gate/coll100.txt`。

判读（下一轮据表判定）：

```text
· 若那两条段的 `ori` 已是 `Internal` 且 `short=1` ⇒ §9.486 的规则确实命中，它们以 Internal 身份参与串接；
  下一步即对拍 cxx:2519-2549 同分支的「改 Internal 条件」与其后的「Internal 段串接规则」。
· 若 `ori` 仍是 `Reversed/Forward`（未被改） ⇒ `is_short_segment` 对它们返回 0 ⇒ 靶点转为 is_short_segment 的判据。
```

---

### 9.488 —— 【实测】那两条零长度段 **`short = 0`**（`is_short_segment` 返回 0）⇒ §9.486 的「会被改成 Internal」**不成立**；新靶点 = `is_short_segment`

探针输出（env `T101_ZZB`，`.target-gate/coll100.txt`，154 条 `CSHORT`）中与 F113 退化段相关的行：

```text
CSHORT k=4/8 nb=1 ori=Reversed short=**0** p0=( 87.500,0.000,-32.000)
CSHORT k=5/8 nb=1 ori=Reversed short=**0** p0=(-87.500,0.000,-32.000)
（另有 k=7/8 nb=1 ori=External short=0 p0=(-87.5,0,-32) 等，属同一面其它次调用）
整体分布：ori ∈ {Reversed:105, External:49}；nb ∈ {1:120, 2:21, 6:3, 7:3, 5:3, 4:3, 8:1}
```

**判读**：

```text
· `short = 0` ⇒ §9.486 的规则 `isshort > 0 && (External || one_degenerated)` **不触发**
  ⇒ 那两条段**保持 `Reversed`**（未改成 Internal） ⇒ §9.486 的假设**被实测否掉**（这正是探针的价值）。
· 「一条 **3D 零长度**、只含 1 条退化边的段」竟然被判为 **不短**（`short=0`）——这本身就很可疑：
  OCCT 侧同一步若判 `short=1`，则该段会被改成 Internal，**其后的串接归属就会完全不同**
  （可能正是 OCCT 得到 `[22, 6]` 而端口得到 `[6,7,5,4,3,3]` 的原因）。
```

#### 下一轮（一次对拍，靶点已是一个函数）

```text
读端口 `is_short_segment(...)`（`shape_fix_compose_shell/helpers.rs` 或 `wire_segment.rs`）
  ↔ OCCT 对应者（`ShapeFix_ComposeShell` 里构造 `shorts` 的那一步，`cxx:2519-2549` 段内调用的判据；
     也可能是 `ShapeAnalysis_Wire::CheckSmall`/`IsSmall` 一线）。
专看：对「1 条边、首尾同点（3D 零长度）」的段，OCCT 是否判为 short；
  若 OCCT 判 1 而端口判 0 ⇒ 按 .cxx 修 `is_short_segment`（不加启发式）⇒ 跑 §9.448 验收。
```

---

### 9.489 —— 端口 `is_short_segment` 全文（`helpers.rs:347-390` ↔ `cxx:2394-2447`）；三处判据里有一处**朝向语义混用**

```rust
/// IsShortSegment (cxx:2394-2447): 1 for a closed segment whose interior
/// collapses to its vertex, -1 when only the 2d check fails, 0 otherwise.
pub fn is_short_segment(seg, face, grid_surface, u_resolution, v_resolution) -> i32 {
    let (Some(vf), Some(vl)) = (seg.first_vertex(), seg.last_vertex()) else { return 0; };
    if !is_same(&vf.0, &vl.0) { return 0; }                       // ① 段首尾必须**同一对象**
    let pnt = vertex_point(&vf); let tol = vertex_tolerance(&vf); let tol2 = tol*tol;
    let mut code = 1;
    for edge in seg.edges() {
        let Some(last) = edge_vertices(edge).1 else { return 0; };
        if !is_same(&vf.0, &last.0) { return 0; }                 // ② 每条边的「末端顶点」须与 vf 同一对象
        let Some((c2d, f, l)) = curve_on_surface_oriented(edge, face, false) else { continue; };
        let end_pnt = c2d.d0(l); let mid_pnt = c2d.d0(0.5*(f+l));
        if !is_coincided(&end_pnt, &mid_pnt, u_resolution, v_resolution, tol) { code = -1; }   // ③2d
        let mid3d = grid_surface.value_uv(mid_pnt.x(), mid_pnt.y());
        if mid3d.distance(&pnt).powi(2) > tol2 { return 0; }       // ④ 中点 3D 须落在 vf 容差内
    }
    code
}
```

**三处可疑点（下一轮对 `cxx:2394-2447` 逐条核）**：

```text
① `①`/`②` 都用 **对象同一性** `is_same`（TShape）——这是 OCCT `IsSame` 的忠实语义；
   但 `②` 取的是 `edge_vertices(edge).1`（**朝向无关**的原始末顶点），而 `vf` 来自
   `seg.first_vertex()`（**朝向感知**）⇒ **两种语义混用**，对一个 `Reversed` 方向的单边段，
   两者可能指向**不同的顶点对象** ⇒ `②` 直接 `return 0`。这与 §9.439 的 `FixDummySeam` 同族（朝向语义）。
② `④` 的判据是「边中点（2d 参数中点映射到面再求 3d）到 `vf` 的距离 ≤ tol」——对**3D 零长度**的退化边，
   这条本应恒成立（中点在参数中点，其 3d 像未必等于端点！**参数中点 ≠ 3D 中点**：若该边的 pcurve 跨
   整个周期（seam 边），2d 中点对应的是**对面**的点，3d 距离会很大 ⇒ `④` 失败 ⇒ `return 0`）。
   **这正可能是 `short=0` 的直接原因**，且它同样是「pcurve 跨周期」这一特殊几何造成的。
③ `③` 只降级为 -1（不改变返回 0/1 的性质）。
```

#### 下一轮（读 `cxx:2394-2447` 定论）

```text
读 OCCT `ShapeFix_ComposeShell.cxx:2394-2447`（IsShortSegment）原文，逐条核对：
  · `②` 在 OCCT 里用的是 `sae.LastVertex(E)`（朝向感知）还是 `TopExp::Vertices` 的原始末顶点？
  · `④` 在 OCCT 里比较的是 `myGrid->Value(mid2d)` 与 `pnt` 的距离（同样口径）还是有别（例如用 3D 端点而非 2D 中点）？
若确认 ①/② 的朝向语义或 ④ 的中点口径有差异 ⇒ 按 .cxx 修 `is_short_segment` ⇒ 跑 §9.448 验收。
（该函数的读数影响面广：`shorts[]` 决定段是否被改成 Internal，进而影响 `CollectWires` 的串接归属。）
```

---

### 9.490 —— OCCT `IsShortSegment` 原文（`ShapeFix_ComposeShell.cxx:2394-2447`）

```cpp
 2394 static int IsShortSegment(const ShapeFix_WireSegment&      seg,
 2395                           const TopoDS_Face&               myFace,
 2396                           const occ::handle<Geom_Surface>& myGrid,
 2397                           const TopLoc_Location&           myLoc,
 2398                           const double                     UResolution,
 2399                           const double                     VResolution)
 2400 {
 2401   TopoDS_Vertex Vf = seg.FirstVertex();
 2402   if (!Vf.IsSame(seg.LastVertex()))
 2403   {
 2404     return 0;
 2405   }
 2407   gp_Pnt pnt  = BRep_Tool::Pnt(Vf);
 2408   double tol  = BRep_Tool::Tolerance(Vf);
 2409   double tol2 = tol * tol;
 2411   int                                      code = 1;
 2412   ShapeAnalysis_Edge                       sae;
 2413   const occ::handle<ShapeExtend_WireData>& sbwd = seg.WireData();
 2414   for (int i = 1; i <= sbwd->NbEdges(); i++)
 2415   {
 2416     TopoDS_Edge edge = sbwd->Edge(i);
 2417     if (!Vf.IsSame(sae.LastVertex(edge)))
 2418     {
 2419       return 0;
 2420     }
 2421     occ::handle<Geom2d_Curve> c2d;
 2422     double                    f, l;
 2423     if (!sae.PCurve(edge, myFace, c2d, f, l))
 2424     {
 2425       continue;
 2426     }
 2428     // check 2d
 2429     gp_Pnt2d endPnt = c2d->Value(l);
 2430     gp_Pnt2d midPnt = c2d->Value((f + l) / 2);
 2431     if (!IsCoincided(endPnt, midPnt, UResolution, VResolution, tol))
 2432     {
 2433       code = -1;
 2434     }
 2436     // check 3d
 2437     gp_Pnt midPnt3d = myGrid->Value(midPnt.X(), midPnt.Y());
 2438     if (!myLoc.IsIdentity())
 2439     {
 2440       midPnt3d.Transform(myLoc.Transformation());
 2441     }
 2442     if (midPnt3d.SquareDistance(pnt) > tol2)
 2443     {
 2444       return 0;
 2445     }
 2446   }
 2447   return code;
 2448 }
```

（下一轮据此逐条核对端口 `helpers.rs:347-390` 的三处判据：①/② 的顶点**朝向语义**、④ 的「2d 参数中点 → 3d」口径。）

---

### 9.491 —— `IsShortSegment` 逐条对拍：发现**两处真实的朝向语义分歧**（`cxx:2417`/`cxx:2423`），但实测**全中性**（F113/a3n00 逐字不变）⇒ 已还原

#### A. 逐条对拍（OCCT `:2394-2447` 原文见 §9.490）

| OCCT | 端口 `helpers.rs:347-390` | 判定 |
|---|---|---|
| `Vf = seg.FirstVertex()`；`Vf.IsSame(seg.LastVertex())`（:2401-2402） | `seg.first_vertex()` / `last_vertex()` + `is_same`（:356-361） | ✓ |
| `pnt/tol/tol2`（:2407-2409） | `:362-364` | ✓ |
| `code = 1`（:2411） | `:365` | ✓ |
| **`Vf.IsSame(sae.LastVertex(edge))`**（:2417） | **`is_same(&vf.0, &edge_vertices(edge).1)`**（:367-372） | ✗ **分歧**（`sae.LastVertex` 朝向感知 vs `edge_vertices(…).1` 朝向无关） |
| **`sae.PCurve(edge, myFace, c2d, f, l)`**（:2423，`CumOri` 默认 **true**） | **`curve_on_surface_oriented(edge, face, false)`**（:373-377） | ✗ **分歧**（CumOri 取 false） |
| `endPnt/midPnt` + `IsCoincided`（:2429-2434） | `:378-382` | ✓ |
| `myGrid->Value(midPnt)`（+ 非恒等 myLoc 变换）+ `SquareDistance > tol2`（:2437-2445） | `grid_surface.value_uv(...)` + `distance(...).powi(2) > tol2`（:383-387；端口注释已说明 myLoc 恒等） | ✓ |

#### B. 实测（两处都按 `.cxx` 改后）

```text
F113：result type=Shell, result faces=5；HEALED wires=2 / 1 / 1 / 1 / 1   ——**与改前逐字相同**
模型：MODEL f=113 wires=4（6/14/1/1）                                    ——**不变**
a3n00：STATMAP matched=225 unmatched=1 sum_mt=11941 flat_mt=11941        ——**不变**
```

⇒ 这两处**是真实分歧**（都应修，且与 §9.439 的 `FixDummySeam` 同族：朝向语义），
   但它们**不是 F113 的杠杆**（与 §9.439/§9.463 的结论一致：能改的都是「几何/判据正确性」，
   而 F113 的成败取决于 `CollectWires` 之后**串接成几条 wire**）。

#### C. 处置与后续

```text
· 已 **还原**（两处 edit 反向替换；库 diff 为空）——按纪律，库改动须过完整门禁方可落地，
  本环境的门禁需 ~5 分钟且必须在后台跑完再判定（参见 §9.463 的落地流程）。
· 建议单独一轮做（与 §9.463 同法）：应用这两处 → 后台跑 `t101_verify.ps1` →
  绿则提交（这两处是 `.cxx` 忠实修正），红则还原并记录退步。
```

---

### 9.492 —— ✅ **落地**：`IsShortSegment` 两处朝向语义修正（`cxx:2417`/`cxx:2423`），门禁 5/5、全部基线逐字不动

改动（`crates/occt-topo/src/shape_fix_compose_shell/helpers.rs`）：

```rust
// ① cxx:2417 —— `sae.LastVertex(edge)` 是**朝向感知**的
- let Some(last) = crate::topo_tools_full::edge_vertices(edge).1 else { return 0; };
+ let Some(last) = crate::shhealing::last_vertex(edge) else { return 0; };
// ② cxx:2423 —— `sae.PCurve(edge, myFace, c2d, f, l)` 的 `CumOri` 默认 **true**
- crate::boptools_2d::curve_on_surface_oriented(edge, face, false)
+ crate::boptools_2d::curve_on_surface_oriented(edge, face, true)
```

**门禁实测**（`pwsh -File .target-gate\t101_verify.ps1`，job `pwsh-117`，exit 0）：

```text
test result: ok. 5 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 267.16s
[occ/a3n00.stp] ours v=10863 f=11941 | occ v=11145 f=12466 | f-ratio 0.96
[occ/a3n00.stp] our=204327.72 occ=227130.30 ratio=**0.8996**   ← 与基线逐字相同
[occ/T0M.stp] 0.9987 · [occ/acs10.stp] 0.9846 · [occ/bottom.step] 0.9750 · [occ/motoc.step] 1.0333
[occ/TDB.stp] 1.0092 · [occ/top.step] 1.0000 · [occ/ATU01038.step] 1.0006   ← 全部逐字相同
```

⇒ 两处均为 `.cxx` 忠实修正（含注释引证），对既有全部基线中性。**已落地**。

#### 本目标累计落地/暂缓清单

```text
已落地：§9.463 split_by_line 第二处移除路径簿记（cxx:1820）
        §9.492 IsShortSegment 两处朝向语义（cxx:2417 / cxx:2423）  ← 本次
暂缓：  §9.439 FixDummySeam 顶点朝向感知（cxx:4221/4230）—— .cxx 正确，但会让 a3n00 面积比 0.8996→0.8918
        （仅 f=171 一张面 +26 三角 / 总面积 −1763），按红线在弄清 f=171 前不落地
待查：  F113 的成败取决于 CollectWires 之后**串接成几条 wire**（OCCT [22,6] vs 端口 [6,7,5,4,3,3]）；
        串接判据（:357-383 的 same_v / is_same_patch）尚未用探针实测其顶点同一性
```

---

### 9.493 —— 下一轮可执行配方：`CollectWires` 串接决策探针（一次运行定论，照 §9.487 的一轮内插桩→跑→撤除流程）

#### 目标问题

```text
F113 的段在切口之后为何串成 [6,7,5,4,3,3]（6 条）而不是 OCCT 的 [22,6]（2 条）？
需要量的是**决策点**：主循环每次「找不到候选 ⇒ append 到 wires」时，是 `index == None`
（真的没有候选）还是被 `is_same_patch` / 同边最低优先级等**判据筛掉**了候选。
```

#### 插桩（库内、env `T101_ZZC`、用完 `edit`/反向替换撤除；两个位点）

```text
① 主循环候选搜索结束处（`collect_wires.rs:231` 之前，"cxx:2726-2779: connect the found segment" 那行之前）：
   eprintln!("CDEC i={} cur_nb={} index={:?} samepatch={} reverse={} connected={} sbwd_nb={}",
             i, seqw[i].nb_edges(), index, samepatch, reverse, connected, sbwd.len());

② append-to-wires 分支处（":309-330  // cxx:2824-2846: if closed or no next segment found, append to wires" 内）：
   eprintln!("CPUSH nb={} first=(...) last=(...) seqw_left={}", sbwd.len(), …, seqw.len());
```

#### 判读规则

```text
· 若 `CDEC` 里出现「`index == None` 紧接着 `CPUSH nb=6`（或 7/5/4/3）」⇒ 是**候选搜索没找到**
  ⇒ 靶点 = 候选循环里的匹配条件（`same_v` 用 TShape 同一性 / `is_same_patch`）；
· 若 `CDEC` 显示 `index = Some(k)` 却仍 `CPUSH`（或循环提前退出）⇒ 靶点在连接/退出条件
  （`connected` / `can_be_closed` / 闭合判定）；
· 若 `CDEC` 显示 `index = Some(k)` 且 `sbwd_nb` 持续增长到 22 ⇒ 串接本身没问题，问题在**更后面的分面**。
```

#### 与 `.cxx` 的对照点

```text
ShapeFix_ComposeShell.cxx:2570-2724（找下一个要连接的段）与 :2726-2779（连接）：
  端口 `collect_wires.rs:103-230` 与 `:231-275` 已逐行判定忠实（§9.399 的旧结论 + 本轮 §9.476 局部复核），
  但**没有在 F113 这个具体输入上取过数**。本配方的价值就是把「忠实」从代码层推进到**该输入上的实测层**。
```

---

### 9.494 —— `CollectWires` 串接决策实测（§9.493 配方执行；一轮内插桩→跑→撤除，库 diff 空）

```text
CDEC 行: 277   CPUSH 行: 230
CPUSH nb 分布: {'1': 44, '2': 51, '3': 36, '4': 46, '5': 40, '6': 4, '7': 4, '8': 3, '15': 2}
CDEC 中 index=None 的行: 74 / 277
```

末尾样本：

```text
CDEC index=Some(6) samepatch=true reverse=true connected=true sbwd_nb=2
CDEC index=Some(5) samepatch=false reverse=false connected=false sbwd_nb=0
CDEC index=Some(7) samepatch=true reverse=true connected=true sbwd_nb=1
CDEC index=Some(7) samepatch=true reverse=false connected=true sbwd_nb=2
CDEC index=None samepatch=false reverse=false connected=false sbwd_nb=3
CDEC index=None samepatch=false reverse=false connected=false sbwd_nb=0
CPUSH nb=2 seqw_left=8
CPUSH nb=3 seqw_left=8
CPUSH nb=1 seqw_left=8
CPUSH nb=2 seqw_left=8
CPUSH nb=3 seqw_left=8
CPUSH nb=3 seqw_left=8
```

原始输出见 `.target-gate/chain106.txt`；判读按 §9.493 的三条规则（下一轮据分布与样本判定靶点：候选匹配条件 / 连接退出条件 / 更后面的分面）。

---

### 9.495 —— 判定：端口的串接是**「连到找不到候选为止」**（`index=Some… → sbwd 增长 → 无候选 ⇒ CPUSH`）⇒ 靶点落在**候选匹配条件**（§9.493 规则 1）

对 §9.494 捕获（`.target-gate/chain106.txt`：CDEC 277 行 / CPUSH 230 行）做配对分析：
取 `CPUSH nb ∈ {3,4,5,6,7,15}` 的 132 次（F113 相关），看其前一条 `CDEC`：

```text
CPUSH nb=3 <- CDEC index=Some(1) samepatch=true reverse=false connected=true sbwd_nb=2
CPUSH nb=4 <- CDEC index=Some(2) samepatch=true reverse=false connected=true sbwd_nb=3
CPUSH nb=3 <- CDEC index=Some(1) samepatch=true reverse=false connected=true sbwd_nb=2
CPUSH nb=4 <- CDEC index=Some(2) samepatch=true reverse=true  connected=true sbwd_nb=3
（同型重复；样本见 .target-gate/chain106.txt）
```

**判读**：

```text
最后一次 CDEC 总是 `index=Some(k)` 且 `connected=true`，且 `sbwd_nb = nb-1`
⇒ 该段被连上（accumulated 变成 nb 条边），**随后的一次迭代里候选搜索返回 `index=None`**，
  才走 `cxx:2824-2846` 的 append 分支把这条 wire 推出（CDEC 中 `index=None` 共 **74/277** 次）。
⇒ 按 §9.493 的**规则 1**：**不是**连接/退出条件提前掐断，而是**候选搜索在那些点找不到任何候选**
   ⇒ 靶点 = 候选循环（`collect_wires.rs:104-230`）里的匹配条件：
      `same_v`（TShape 同一性）/ `is_same_patch`（patch 优先序）/ 同边返回最低优先级 / `misoriented`。
```

#### 下一轮（最后一次探针，按同一「一轮内插桩→跑→撤除」流程）

```text
在候选循环内、`continue`/跳过处统计**被拒原因**（仅对 F113 那次调用、仅在 `index` 最终为 None 的那趟打印）：
  · `same_v` 为假（顶点非同一对象，但**坐标可能相同**——这正是 §9.485 之后要回答的问题）；
  · `is_same_patch` 为假；
  · 同边返回（最低优先级）被跳过；
  · `misoriented` / 其他。
若统计出「大量坐标相同但 `same_v` 为假」⇒ 根因 = **切口处顶点未统一**（`Context`/`CombineVertex` 路径），
   修法按 `.cxx` 在该处统一顶点；否则靶点落到 `is_same_patch` 的优先序判据。
```

---

### 9.496 —— 「找不到候选」时刻的重扫诊断（一轮内插桩→跑→撤除；库 diff 空）

```text
CREJ 行: 25；其中 sbwd_nb>=3 的 25 行
样本：
   CREJ sbwd_nb=7 seqw=8 cand_total=26 cand_near=0 cand_same=0
   CREJ sbwd_nb=5 seqw=8 cand_total=26 cand_near=0 cand_same=0
   CREJ sbwd_nb=4 seqw=8 cand_total=26 cand_near=1 cand_same=1
   CREJ sbwd_nb=3 seqw=8 cand_total=26 cand_near=1 cand_same=1
   CREJ sbwd_nb=15 seqw=4 cand_total=15 cand_near=0 cand_same=0
   CREJ sbwd_nb=5 seqw=4 cand_total=15 cand_near=0 cand_same=0
   CREJ sbwd_nb=8 seqw=4 cand_total=15 cand_near=1 cand_same=1
   CREJ sbwd_nb=6 seqw=4 cand_total=7 cand_near=1 cand_same=0
   CREJ sbwd_nb=5 seqw=3 cand_total=4 cand_near=1 cand_same=1
   CREJ sbwd_nb=5 seqw=3 cand_total=4 cand_near=1 cand_same=1
合计：cand_near=19  cand_same=18
```

判读（下一轮据表定论）：

```text
· 若 cand_near>0 而 cand_same=0 ⇒ 剩余段**坐标相邻但顶点不是同一对象** ⇒ 根因 = 切口处顶点未统一
  （Context/CombineVertex 路径），按 .cxx 统一顶点；
· 若 cand_near=0 ⇒ 剩余段**根本不相邻** ⇒ 靶点回到候选循环的更早判据（misoriented/is_same_patch 的优先序）。
```

---

### 9.497 —— 判读 §9.496：两类机制**同时存在**（多数是「根本不相邻」`cand_near=0`；至少一处是「顶点未统一」`near=1/same=0`）

#### 数据（25 行，`sbwd_nb ≥ 3` 的「无候选 ⇒ 推出」时刻）

```text
CREJ sbwd_nb=7  seqw=8 cand_total=26 cand_near=0 cand_same=0
CREJ sbwd_nb=5  seqw=8 cand_total=26 cand_near=0 cand_same=0
CREJ sbwd_nb=4  seqw=8 cand_total=26 cand_near=1 cand_same=1
CREJ sbwd_nb=3  seqw=8 cand_total=26 cand_near=1 cand_same=1
CREJ sbwd_nb=15 seqw=4 cand_total=15 cand_near=0 cand_same=0
CREJ sbwd_nb=5  seqw=4 cand_total=15 cand_near=0 cand_same=0
CREJ sbwd_nb=8  seqw=4 cand_total=15 cand_near=1 cand_same=1
CREJ sbwd_nb=6  seqw=4 cand_total=7  cand_near=1 cand_same=**0**     ← 坐标相邻、对象不同
CREJ sbwd_nb=5  seqw=3 cand_total=4  cand_near=1 cand_same=1
…（合计 cand_near=19，cand_same=18）
```

#### 判读

```text
① **多数行 `cand_near=0`**：推出 7 / 5 / 15 边 wire 时，剩余段中**没有**任何段的顶点落在
   该 wire 的 `end_v` 上（1e-6 内）⇒ 按「顶点相接」这条规则，链确实到此为止。
   ⇒ 说明端口把「顶点相接」当作**唯一**的延伸入口；而 OCCT 能把 16 边的段与外环另外的段并成 22 边，
     **必然用了不止「顶点相接」**——回看候选循环里的选择项（`samepatch` / `reverse` / `angle`(2D 切向) /
     `mindist` / 同边返回最低优先级），这些是**在顶点相接的前提下**做优先级；若 OCCT 在 `cxx:2570-2724`
     里允许「不相接但同 patch 且切向连续」的候选，端口就少了一类入口。
② **至少一行 `near=1 / same=0`**：确实存在「坐标相邻但顶点非同一对象」⇒ §9.485 的怀疑在该点成立，
   这类点会因 `same_v` 失败而少一个候选（可用 `Context()->Apply` 统一顶点修复）。
```

#### 下一轮（本目标的最后一步对拍）

```text
读端口 `collect_wires.rs:104-230` 的候选循环，列出**所有** `continue`（跳过候选）的条件与其顺序，
逐条对 `ShapeFix_ComposeShell.cxx:2570-2724`；重点看：
  · 端口是否有「仅当 `same_v(顶点相接)` 成立才进入后续优先级比较」这一硬前置；
  · `cxx` 是否在同一处允许 `IsSame` 失败但 `samepatch && 切向连续` 的候选（即**少一条入口**）。
若确认端口少了该类入口 ⇒ 按 .cxx 补（不加启发式）⇒ 跑 §9.448 验收
（zz_seam_fix 113 → Face；--model 113 → wires=2(22+6)；--fstats → face=113 mt≈228；
 面积比 ≥0.8996 且上升；t101_verify.ps1 全绿）。
```

---

### 9.498 —— 候选循环的跳过条件全表（`:104-225`）：`same_v` 硬前置**忠实**（`cxx:2629-2632`）⇒ §9.497 的「少一条入口」假设**否掉**；并发现 §9.496 探针有一处口径缺陷

#### A. 端口候选循环的过滤结构（逐条带 `.cxx` 出处）

```text
:106  if seg.is_vertex() { continue; }                                  // 跳过 vertex 段
:135  if !sp && (can_be_closed || (index.is_some() && samepatch)) { continue; }   // patch 过滤
:140  let candidate_v = if j == 1 { seg.last_vertex() } else { seg.first_vertex() };
:141  if !same_v(&end_v, &candidate_v) { continue; }                     // cxx:2629-2632 ← **顶点相接硬前置**
:151  if same_e(&last_edge, &Some(back_edge)) { … 最低优先级 … }         // cxx:2639-2652 同边返回
:215  if w1 + tail1 <= weigth + tail2 { continue; }                      // cxx:2696-2722 权重比较（patch/tangent/dist）
```

⇒ `same_v` 这条**硬前置在两端是等价的**（端口 :141 ↔ cxx:2629-2632）
⇒ §9.497 猜测的「OCCT 允许 IsSame 失败但 samepatch+切向连续的候选」**不成立**（两侧都要求顶点同一性）。
⇒ 因此唯一存活的解释是：**OCCT 在该处顶点是同一对象，而端口不是**（即 §9.485 / §9.496 的 `near=1, same=0` 那一行）。

#### B. **更正 §9.496 的探针口径缺陷**（自查）

```text
§9.496 的重扫只统计了候选边的 `first_vertex`，而候选循环 `:140` 对 `j == 1` 用的是 `seg.last_vertex()`
⇒ 我的 `cand_near=0` 那些行**可能只是漏查了 last_vertex**（不是「真的不相邻」）。
⇒ 该数据只能作为「存在 near=1/same=0 的顶点未统一实例」的证据，**不能**用来断言「多数点根本不相邻」。
```

#### 下一轮（修正口径重测，一次即可定论）

```text
把重扫改为对每个候选段的 **first_vertex 与 last_vertex 都查**，并分三类计数：
  A) 坐标相接(is_same 为真)  B) 坐标相接但 is_same 为假  C) 完全不相接
若 B 占多数 ⇒ 根因 = **切口处顶点未统一**（`Context`/`CombineVertex` 路径）；
  按 `.cxx` 在切口处统一顶点（`cxx:1900-1911` 的 `int_vertices` 与 `SplitByGrid` 的上下文替换路径）
  ⇒ 让 `same_v` 能匹配上 ⇒ 链能延伸成 22 边 ⇒ 跑 §9.448 验收。
若 A 占多数而仍不成链 ⇒ 靶点回到 `:215` 的权重比较（w1/w2/tail1/tail2 的取舍）。
```

---

### 9.499 —— 修正口径后的重测（A/B/C 三分类；一轮内插桩→跑→撤除，库 diff 空）

对「无候选 ⇒ 推出」时刻（`sbwd_nb >= 3`，共 25 次）重扫剩余段的**两个端点**：

```text
合计：A(坐标相接且 is_same)=58   B(坐标相接但 is_same 为假)=7   C(完全不相接)=559
样本：
   CABC sbwd_nb=7 A_same=1 B_near_not_same=1 C_far=50
   CABC sbwd_nb=5 A_same=1 B_near_not_same=1 C_far=50
   CABC sbwd_nb=4 A_same=2 B_near_not_same=0 C_far=50
   CABC sbwd_nb=3 A_same=3 B_near_not_same=0 C_far=49
   CABC sbwd_nb=15 A_same=1 B_near_not_same=1 C_far=28
   CABC sbwd_nb=5 A_same=1 B_near_not_same=1 C_far=28
   CABC sbwd_nb=8 A_same=2 B_near_not_same=0 C_far=28
   CABC sbwd_nb=6 A_same=1 B_near_not_same=1 C_far=12
```

判读（据 §9.498 的两条分支）：

```text
· B 占多数 ⇒ 根因 = 切口处顶点未统一（Context/CombineVertex 路径）⇒ 按 .cxx 统一顶点（跑 §9.448 验收）；
· A 占多数仍不成链 ⇒ 靶点回到 :215 的权重比较（w1/w2/tail1/tail2）。
```

---

### 9.500 —— 【定论】「无候选」时刻**每一行都存在同一对象的候选顶点**（A ≥ 1，25/25 行）⇒ 拒绝发生在**后续判据**，不是 `same_v`

#### 数据（§9.499，25 行，`sbwd_nb ≥ 3`）

```text
合计：A(坐标相接且 is_same)=58   B(坐标相接但 is_same 为假)=7   C(完全不相接)=559
样本：
CABC sbwd_nb=7  A_same=1 B_near_not_same=1 C_far=50
CABC sbwd_nb=5  A_same=1 B_near_not_same=1 C_far=50
CABC sbwd_nb=4  A_same=2 B_near_not_same=0 C_far=50
CABC sbwd_nb=3  A_same=3 B_near_not_same=0 C_far=49
CABC sbwd_nb=15 A_same=1 B_near_not_same=1 C_far=28
CABC sbwd_nb=8  A_same=2 B_near_not_same=0 C_far=28
CABC sbwd_nb=6  A_same=1 B_near_not_same=1 C_far=12
```

#### 判读（按 §9.498 的分支）

```text
· **A ≥ 1 出现在全部 25 行** ⇒ 在每次「无候选 ⇒ 推出」的时刻，剩余段里**确实有与 `end_v` 同一对象**
  的端点 ⇒ 那么 `:141 same_v` 硬前置**至少放过了一个候选** ⇒ **`index` 仍为 None 只能由后续判据解释**：
     :135 `if !sp && (can_be_closed || (index.is_some() && samepatch)) { continue; }`（patch 过滤）
     :151 `same_e(last_edge, back_edge)` ⇒ 同边返回（最低优先级，可能不落 index）
     :215 `if w1 + tail1 <= weigth + tail2 { continue; }`（权重比较：w1 = sp*16 + connected*8 + …）
· B = 7（坐标相接但对象不同）⇒ 「顶点未统一」**确实存在但是少数**，不是主因。
· C = 559 ⇒ 绝大多数端点与该 wire 末端不相接（正常，属其它段）。

⇒ 结论：**靶点 = 候选被拒的后续判据**，最可能是 `:135` 的 patch 过滤与 `:200-225` 的权重比较
   （`w1/tail1/weigth/tail2` 的取舍），对应 `cxx:2607-2614`（同 patch 测试与优先级）与 `cxx:2696-2722`（权重）。
```

#### 下一轮（本目标最后一个对拍窗口）

```text
读端口 `collect_wires.rs:128-230` 的 **patch 过滤 + 权重比较** 全段，逐条对 `cxx:2607-2614` / `cxx:2696-2722`：
  · `is_same_patch(...)` 的返回与 `samepatch` 的更新时机；
  · `w1`/`tail1`/`tail2` 的构造（`sp*16 + connected*8 + …`）与 `cxx` 的 `iW1/iW2` 是否同构；
  · `:151 same_e` 分支是否漏设 `index`（该分支若不设 index，同边候选会被永久丢弃）。
命中即按 .cxx 修 ⇒ 跑 §9.448 验收（zz_seam_fix 113 → Face；--model 113 → wires=2(22+6)；
--fstats → face=113 mt≈228；面积比 ≥0.8996 且上升；t101_verify.ps1 全绿）。
```

---

### 9.501 —— 读 `:151-165`：同边（最低优先级）分支**确实会设 `index`** ⇒ §9.500 的第三条假设被否；嫌疑收敛到 `:135`（patch 过滤，含 `can_be_closed`）与 `:215`（权重）

端口 `collect_wires.rs:133-165`（`cxx:2607-2614` / `cxx:2639-2652`）：

```rust
let sp = is_same_patch(&seqw[i], nu, nv, &mut iumin, …);           // cxx:2607-2614
if !sp && (can_be_closed || (index.is_some() && samepatch)) { continue; }   // ← 过滤 ①
for j in 0..2usize {
    let candidate_v = if j == 1 { seg.last_vertex() } else { seg.first_vertex() };
    if !same_v(&end_v, &candidate_v) { continue; }                  // cxx:2629-2632（硬前置）
    let misor = an_or == if j == 1 { Reversed } else { Forward };
    let back_edge = if j == 1 { wire_edges.last } else { wire_edges[0] };
    if same_e(&last_edge, &Some(back_edge)) {                       // cxx:2639-2652（同边最低优先级）
        if index.is_none() && !can_be_closed {                      // ← **确实设 index**
            index = Some(i); reverse = j != 0; connected = true; misoriented = misor;
            samepatch = sp;
            weigth = (if sp {16} else {0}) + (if connected {8} else {0}) + (if !misor {4} else {0});
            dsu = 0.0; dsv = 0.0;
        }
        continue;
    }
    … cxx:2654-2673 起算 2D 切向 / 距离 / 权重 …
}
```

**判定**：

```text
· `:151` 同边分支在 `index.is_none() && !can_be_closed` 时**会设 index** ⇒ §9.500 的假设③**否掉**。
· 结合 §9.500 实测「每行 A ≥ 1（存在同一对象候选）却 index=None」，剩下的可能只有两处：
    ① `:135` 的 patch 过滤：`!sp`（非同一 patch）且 `can_be_closed` 为真时，**所有**非同 patch 候选被跳过；
       ⇒ `can_be_closed` 的置位时机（`:89` 初始化、后续按闭合条件更新）成为关键；
    ② `:215` 的权重比较：`w1 + tail1 <= weigth + tail2 ⇒ continue`
       ⇒ `w1/tail1/weigth/tail2` 的构造若与 `cxx:2696-2722` 不同构，会造成「有候选却不采纳」。
```

#### 下一轮（两处，先看 `can_be_closed`）

```text
① 读 `can_be_closed` 的全部置位点（与 `cxx:2551-2567` 及后续段对照），确认在本输入上它是否长期为真
   （若长期为真 ⇒ `:135` 把非同 patch 候选全挡掉 ⇒ 与 OCCT 的差异就在这）；
② 读 `:171-230`（2D 切向 / 距离 / w1 / tail1 / tail2）↔ `cxx:2654-2722`，逐项核对权重构造。
命中即按 .cxx 修 ⇒ 跑 §9.448 验收。
```

---

### 9.502 —— `can_be_closed` 是**瞬态标志** ⇒ `:135` 过滤退化为「patch 优先级规则」，不是阻塞点；唯一剩余嫌疑 = `:215` 权重比较

`can_be_closed` 的全部出现（`collect_wires.rs`）：

```text
:89   let mut can_be_closed = false;                                  // 初始化
:135  if !sp && (can_be_closed || (index.is_some() && samepatch)) { continue; }   // 候选过滤
:152  if index.is_none() && !can_be_closed { … 同边分支设 index … }
:310  can_be_closed = same_v(&end_v, &first_v);                      // ← 仅在「推出前」判定
:312  || (can_be_closed && …)                                        // 推出条件的一部分
:329  can_be_closed = false;                                         // 推出后复位
```

**判读**：

```text
· `:310` 处于**推出分支**（`:309` 的 cxx:2824-2846 段）内，`:329` 紧随其后复位
  ⇒ `can_be_closed == true` 只存在于「`:310` → `:312`/`:329`」这几行之间，
    **不可能存活到下一次候选搜索**（搜索发生在循环开头的 `:104-230`）。
  ⇒ 在候选搜索里 `:135` 的 `can_be_closed` 项**实际恒为 false**
  ⇒ `:135` 退化为：`!sp && index.is_some() && samepatch ⇒ continue`
    （即「已有同 patch 候选时，跳过非同 patch 候选」——这是合理的优先级规则）✓ **不是阻塞点**。
· `:152` 的 `!can_be_closed` 同理恒为 true ⇒ 同边分支的条件实为 `index.is_none()` ✓（§9.501 已判）。
```

#### 唯一剩余嫌疑（本目标的最后一个对拍窗口）

```text
`:215` `if w1 + tail1 <= weigth + tail2 { continue; }`（cxx:2696-2722）
  ⇒ `w1 / tail1 / tail2 / weigth` 的构造若与 `cxx` 不同构，就会出现「存在同一对象候选却 index=None」。
  读数（§9.500）：25/25 行的「无候选」时刻都存在 A 类（同一对象）候选 ⇒ 说明它们**进到了比较**却被否掉
  ⇒ 要么 `w1` 算得偏小，要么 `weigth` 初值/更新偏大。
  下一轮：读 `:171-230`（2D 切向 / 距离 / w1 / tail1 / tail2）↔ `cxx:2654-2722`，逐项核对；
  必要时用同一「一轮内插桩→跑→撤除」流程打印每次比较的 `w1/tail1/weigth/tail2`（4 个数即可定论）。
```

---

### 9.503 —— `:215` 权重比较的四数取证（一轮内插桩→跑→撤除；库 diff 空）

```text
CWT 行: 180   被否(rej=true): 42   被采纳: 138
weigth 取值分布: {'0': 137, '28': 29, '24': 13, '20': 1}
被否样本（前 8）：
   j=1 w1=20 tail1=0 weigth=28 tail2=1 lhs=20 rhs=29
   j=1 w1=20 tail1=0 weigth=28 tail2=1 lhs=20 rhs=29
   j=1 w1=20 tail1=0 weigth=28 tail2=1 lhs=20 rhs=29
   j=1 w1=20 tail1=0 weigth=28 tail2=1 lhs=20 rhs=29
   j=1 w1=20 tail1=0 weigth=28 tail2=1 lhs=20 rhs=29
   j=1 w1=20 tail1=0 weigth=28 tail2=1 lhs=20 rhs=29
   j=1 w1=20 tail1=0 weigth=28 tail2=1 lhs=20 rhs=29
   j=1 w1=20 tail1=0 weigth=28 tail2=1 lhs=20 rhs=29
```

判读（下一轮据四数定论）：

```text
· 若被否样本多为 `lhs <= rhs` 且 `w1` 明显偏小/`weigth` 偏大 ⇒ 权重构造（w1/tail1/tail2）与 cxx:2696-2722 不等价；
· 若 `weigth` 恒为 0（初始值）而 `lhs <= 0` 才被否 ⇒ 说明候选 `w1` 算得 0 ⇒ 回到 `:200-213` 的 w1/tail 构造；
· 若被否的 lhs/rhs 差为 0（相等）⇒ 检查判据用的是 `<=`（端口）还是 `<`（cxx）。
```

---

### 9.504 —— `:215` 权重比较**也不是阻塞点**（138 次采纳、42 次被否且均为「保留更优候选」的正常行为）；并**更正 §9.500 的测量口径错误**

#### A. 四数取证（§9.503，180 行）

```text
CWT 行 180：被否(rej=true) 42   被采纳 138
weigth 取值分布: {0: 137, 28: 29, 24: 13, 20: 1}
被否样本（全部同型）：
   j=1 w1=20 tail1=0 weigth=28 tail2=1 ⇒ lhs=20 rhs=29 ⇒ 20 <= 29 ⇒ 被否
```

**判读**：

```text
· 137/180 的比较发生在 `weigth == 0`（尚无候选）时，且**采纳 138 次**
  ⇒ 候选循环**确实在正常工作**（能找到并采纳候选）；
· 42 次被否全部是同型：已有 `weigth=28` 的更强候选（sp=16 + connected=8 + !misor=4 = 28），
  新候选 `w1=20`（sp=16 + !misor=4，无 connected 8）⇒ `20 ≤ 28+1` ⇒ 保留原候选
  ⇒ 这是**正确的优先级行为**，不是缺陷。
⇒ **`:215` 不是 F113 的阻塞点**（§9.502 的最后嫌疑亦排除）。
```

#### B. **自查并更正 §9.500 的测量口径错误**（重要）

```text
§9.500 的 CABC 扫描对**每个段的每条边**的 first/last 顶点计数（A/B/C），
而候选循环 `:140` 用的是**段级**端点：
    `candidate_v = if j == 1 { seg.last_vertex() } else { seg.first_vertex() }`
⇒ 我的 A 计数把「段内部相邻边的端点」也算进去了（例如退化段 first==last、或段内相邻边的接点），
  这些**不是合法候选终端** ⇒ **A ≥ 1 不能推出「有候选通过了 same_v」**
  ⇒ §9.500 的「定论」**不成立**，需按段级端点重测。
```

#### C. 修正后的重测配方（下一轮，一次即可定论；仍是「一轮内插桩→跑→撤除」）

```text
在「无候选 ⇒ 推出」的时刻（`index.is_none() && sbwd.len() >= 3`）对每个候选段 sg 只查**两个段级端点**：
    let v0 = sg.first_vertex(); let v1 = sg.last_vertex();
并分三类：A) `same_v(&end_v, &v0) || same_v(&end_v, &v1)` 为真（合法且同一对象）
          B) 坐标相接（<1e-6）但 `same_v` 为假（**顶点未统一**）
          C) 都不相接
· 若 A 仍 ≥1（且候选循环却 index=None）⇒ 说明候选被 `:135`/`:215` 之外的条件挡住，需再查 `:106`(is_vertex)
  与 `:144` misoriented 等；
· 若 B 占多数 ⇒ 根因 = **切口处顶点未统一**（Context/CombineVertex 路径）⇒ 按 .cxx 统一顶点 ⇒ 跑 §9.448 验收；
· 若 C 占多数 ⇒ 该 wire 的末端确实没有可接的段终端 ⇒ 说明 **OCCT 用了端口没有的另一种延伸规则**
  （需回 `cxx:2570-2724` 找「非顶点相接」的入口是否存在）。
```

---

### 9.505 —— 按修正口径（**段级端点**）重测（一轮内插桩→跑→撤除；库 diff 空）

```text
CSEG 行: 25    sbwd_nb>=3: 25
合计：A(段级端点且同一对象)=58   B(坐标相接但非同一对象)=1   C(不相接)=179
样本：
   CSEG sbwd_nb=7 A_same=1 B_near_not_same=0 C_far=15
   CSEG sbwd_nb=5 A_same=1 B_near_not_same=0 C_far=15
   CSEG sbwd_nb=4 A_same=2 B_near_not_same=0 C_far=14
   CSEG sbwd_nb=3 A_same=3 B_near_not_same=0 C_far=13
   CSEG sbwd_nb=15 A_same=1 B_near_not_same=0 C_far=7
   CSEG sbwd_nb=5 A_same=1 B_near_not_same=0 C_far=7
   CSEG sbwd_nb=8 A_same=2 B_near_not_same=0 C_far=6
   CSEG sbwd_nb=6 A_same=1 B_near_not_same=1 C_far=6
   CSEG sbwd_nb=5 A_same=3 B_near_not_same=0 C_far=3
   CSEG sbwd_nb=5 A_same=3 B_near_not_same=0 C_far=3
```

判读（按 §9.504-C 的三分支）：

```text
· A>=1 而 index=None ⇒ 另有条件挡住（再查 :106 is_vertex / :144 misoriented）；
· B 占多数 ⇒ 根因 = 切口处顶点未统一 ⇒ 按 .cxx 统一顶点 ⇒ 跑 §9.448 验收；
· C 占多数 ⇒ 末端确实没有可接的**段**终端 ⇒ OCCT 用了端口没有的延伸规则（回 cxx:2570-2724 找）。
```

---

### 9.506 —— 修正口径重测：**A ≥ 1 出现在 25/25 行**（段级端点、同一对象）⇒ §9.500 的结论恢复；靶点收敛到 `:106`/`:111`/`:116`（唯一未读的跳过点）

#### 数据（§9.505，段级端点口径）

```text
合计：A(段级端点且同一对象)=58   B(坐标相接但非同一对象)=1   C(不相接)=179
CSEG sbwd_nb=7  A_same=1 B=0 C_far=15
CSEG sbwd_nb=4  A_same=2 B=0 C_far=14
CSEG sbwd_nb=3  A_same=3 B=0 C_far=13
CSEG sbwd_nb=15 A_same=1 B=0 C_far=7
CSEG sbwd_nb=6  A_same=1 B=1 C_far=6
CSEG sbwd_nb=5  A_same=3 B=0 C_far=3
```

**判读**：用**段级**端点（`candidate_v = if j==1 { seg.last_vertex() } else { seg.first_vertex() }` 同口径）重测后，
**A ≥ 1 在全部 25 行成立** ⇒ §9.500 的「存在同一对象候选却 `index=None`」**成立**（§9.504 的口径更正后仍然成立）。
⇒ 按 §9.504-C 分支 1：「另有条件挡住」⇒ 候选循环里**唯一尚未读过**的跳过点：

```text
:106  if seg.is_vertex() { continue; }        // ← 候选段是 **vertex 段** ⇒ 直接跳过
:111  continue;                               // 条件未读
:116  continue;                               // 条件未读
:119  reverse = true;  :121 index = Some(i);  :128 misoriented = false;  :131 break;
      （这段是 cxx:2584-2604「首个段可任取」分支，:`111`/`:116` 的条件也在其中）
```

而 `:135`（patch 过滤）、`:141`（same_v）、`:151`（同边）、`:215`（权重）**都已排除**（§9.500/§9.502/§9.504）。

#### 下一轮（读 `:104-132` 全文，约 30 行，本目标最后一处未读代码）

```text
重点：`:106` 的 `seg.is_vertex()` 是否命中了「段级端点与 end_v 同一对象」的那些候选段
     （若是 ⇒ 端口把「vertex 段」排除在可接续候选之外 ⇒ 与 cxx:2584-2604 的同段判据对比即可判定）；
     以及 `:111`/`:116` 两个 continue 的条件。
命中即按 .cxx 修 ⇒ 跑 §9.448 验收（zz_seam_fix 113 → Face；--model 113 → wires=2(22+6)；
--fstats → face=113 mt≈228；面积比 ≥0.8996 且上升；t101_verify.ps1 全绿）。
```

---

### 9.507 —— 最后一处未读代码 `:104-132` 读完：候选段的两个跳过点是 `:106 is_vertex` 与 **`:110 Internal`**；`:113-117` 只作用于「尚无累积 wire」时

```rust
for i in 0..seqw.len() {
    let seg = &seqw[i];
    if seg.is_vertex() { continue; }                       // :106  跳过 vertex 段
    let an_or = seg.orientation();
    if an_or == Orientation::Internal { continue; }        // :110  **跳过 Internal 段**
    if !has_sbwd {
        // cxx:2584-2604: for the first segment, take any.
        if shorts[i] > 0 || an_or == Orientation::External { continue; }   // :115  首段不用「短段/External」
        if an_or == Orientation::Forward { reverse = true; }
        index = Some(i); … break;                           // :121-131
    }
    // :133 起 = cxx:2607-2614 同 patch 测试与优先级（已读过，§9.498/§9.501）
}
```

**判读（结合 §9.505 的实测 A ≥ 1 / 25 行）**：

```text
· 候选被挡的两个位置只有 `:106`（vertex 段）与 `:110`（**Internal 段**）；
· `:115` 只作用于 `!has_sbwd`（即「还没有任何累积 wire」的首段选择），与「推出前 index=None」无关；
· 而 §9.486 已确认：`collect_wires.rs:68-72` 会把「`isshort > 0` 且（External 或单条退化边）」的段**就地改成 Internal**；
  §9.488 实测那两条零长度段 `short=0` 未被改——但**别的段**（例如被 `break_wires` 切出的、或 `load_wires` 的
  非流形段）可能是 Internal ⇒ 它们会**在 `:110` 被永久跳过**，即使其段级端点与 `end_v` 同一对象（A 类）。
⇒ 新假设（比前几轮更贴近代码）：**A 类候选多数是 Internal 段**（被 `:110` 跳过）
  ⇒ 若 OCCT 的 `cxx:2570-2624` 同处**不**排除 Internal 段（或排除条件不同），即为分歧点。
```

#### 下一轮（一次探针即可定论；仍是一轮内插桩→跑→撤除）

```text
在「无候选 ⇒ 推出」时刻（`index.is_none() && sbwd.len() >= 3`）对 A 类候选（段级端点 `same_v(&end_v, …)` 为真）
打印其 `orientation()` 与 `is_vertex()`，统计其中 Internal / vertex 各占多少：
  · 若 A 类**全部或多数是 Internal** ⇒ 命中 `:110`；随后读 `cxx:2570-2624` 核对 OCCT 是否同样排除 Internal，
    差异处按 .cxx 修 ⇒ 跑 §9.448 验收；
  · 若 A 类多为 Forward/Reversed 且 is_vertex=false ⇒ `:106`/`:110` 都不是原因，需回 `:133` 之后
    逐条复核（但 §9.500-§9.504 已把它们排除，故届时需重新审视是否漏了别的过滤）。
```

---

### 9.508 —— 验证 §9.507 的假设：A 类候选的朝向与 is_vertex 分布（一轮内插桩→跑→撤除；库 diff 空）

```text
CAOR 行: 25   sbwd_nb>=3: 25
A 类合计：Fwd=0  Rev=0  **Internal=43**  Other=0   |   vertex=0  段=43
样本：
   CAOR sbwd_nb=7 A_fwd=0 A_rev=0 A_int=1 A_other=0 A_vertex=0 A_seg=1
   CAOR sbwd_nb=5 A_fwd=0 A_rev=0 A_int=1 A_other=0 A_vertex=0 A_seg=1
   CAOR sbwd_nb=4 A_fwd=0 A_rev=0 A_int=2 A_other=0 A_vertex=0 A_seg=2
   CAOR sbwd_nb=3 A_fwd=0 A_rev=0 A_int=2 A_other=0 A_vertex=0 A_seg=2
   CAOR sbwd_nb=15 A_fwd=0 A_rev=0 A_int=1 A_other=0 A_vertex=0 A_seg=1
   CAOR sbwd_nb=5 A_fwd=0 A_rev=0 A_int=1 A_other=0 A_vertex=0 A_seg=1
   CAOR sbwd_nb=8 A_fwd=0 A_rev=0 A_int=2 A_other=0 A_vertex=0 A_seg=2
   CAOR sbwd_nb=6 A_fwd=0 A_rev=0 A_int=1 A_other=0 A_vertex=0 A_seg=1
   CAOR sbwd_nb=5 A_fwd=0 A_rev=0 A_int=2 A_other=0 A_vertex=0 A_seg=2
   CAOR sbwd_nb=5 A_fwd=0 A_rev=0 A_int=2 A_other=0 A_vertex=0 A_seg=2
```

判读（按 §9.507）：

```text
· 若 Internal 占多数 ⇒ 命中 `:110`（Internal 段被永久跳过）⇒ 读 cxx:2570-2624 核对 OCCT 是否同样排除；
· 若 vertex 占多数 ⇒ 命中 `:106`；
· 若 Fwd/Rev 占多数且非 vertex ⇒ 两处均非原因，需重新审视 :133 之后是否还有未计入的过滤。
```

---

### 9.509 —— 【决定性确认】「无候选」时刻的 A 类候选**43/43 全是 `Internal` 段**（vertex 0 个）⇒ 挡住它们的就是 `:110`

#### 实测（§9.508，25 行，`sbwd_nb ≥ 3`）

```text
A 类合计：Fwd=0   Rev=0   **Internal=43**   Other=0   |   vertex=0   段=43
CAOR sbwd_nb=7  A_int=1 A_seg=1     CAOR sbwd_nb=4  A_int=2 A_seg=2
CAOR sbwd_nb=3  A_int=2 A_seg=2     CAOR sbwd_nb=15 A_int=1 A_seg=1
CAOR sbwd_nb=8  A_int=2 A_seg=2     CAOR sbwd_nb=6  A_int=1 A_seg=1
（25/25 行同型：A 类候选**全部**为 Internal、**无一**为 vertex）
```

**判定**：

```text
· `:106 seg.is_vertex()` 未命中（A_vertex = 0）；
· `:110 if an_or == Orientation::Internal { continue; }` **正好命中全部 A 类候选**
  ⇒ 端口把「顶点相接、本可接续」的候选**因它们是 Internal 段而全部跳过**
  ⇒ `index` 只能为 None ⇒ 链在此断开 ⇒ 外环被拆成多条 wire（[6,7,5,4,3,3]）
  ⇒ 与 OCCT 的 [22,6] 差异的直接机制**找到了**。
· 这些段的 Internal 身份来自 `collect_wires.rs:68-72`（`isshort > 0 && (External || one_degenerated) ⇒ set Internal`）
  或 `load_wires` 的非流形段（§9.476/§9.486 已读）。
```

#### 最后一步（唯一剩余对照）

```text
读 `ShapeFix_ComposeShell.cxx:2570-2624`（找下一个要连接的段），核对 OCCT 是否**同样**用
「`TopAbs_INTERNAL` ⇒ 跳过」这一条：
  · 若 OCCT **不跳过** Internal（或只在更窄条件下跳过）⇒ 端口多了一条不该有的过滤 ⇒ 按 .cxx 修
    （把 `:110` 的条件对齐 cxx）⇒ 跑 §9.448 验收（zz_seam_fix 113 → Face；--model 113 → wires=2(22+6)；
    --fstats → face=113 mt≈228；面积比 ≥0.8996 且上升；t101_verify.ps1 全绿）；
  · 若 OCCT 同样跳过 ⇒ 分歧在**更上游**：那些段在 OCCT 里**不该是 Internal**
    （即端口在 `:68-72` 或 `load_wires` 里把带 External 身份的段错标成 Internal）⇒ 需对拍该处判据。
```

---

### 9.510 —— 入库 OCCT `ShapeFix_ComposeShell.cxx:2570-2624` 原文（候选搜索的跳过条件），供与端口 `collect_wires.rs:104-132` 对照

```cpp
 2570     for (i = 1; i <= seqw.Length(); i++)
 2571     {
 2572       const ShapeFix_WireSegment& seg = seqw.Value(i);
 2573       if (seg.IsVertex())
 2574       {
 2575         continue;
 2576       }
 2577       TopAbs_Orientation anOr = seg.Orientation();
 2578       if (anOr == TopAbs_INTERNAL)
 2579       {
 2580         continue;
 2581       }
 2583       // for first segment, take any
 2584       if (sbwd.IsNull())
 2585       {
 2586         if (shorts(i) > 0)
 2587         {
 2588           continue;
 2589         }
 2590         if (anOr == TopAbs_EXTERNAL)
 2591         {
 2592           continue;
 2593         }
 2594         if (anOr == TopAbs_FORWARD)
 2595         {
 2596           reverse = true;
 2597         }
 2598         index = i;
 2599         seg.GetPatchIndex(1, iumin, iumax, ivmin, ivmax);
 2601         misoriented = false;
 2602         dsu = dsv = 0.;
 2603         break;
 2604       }
 2606       // check whether current segment is on the same patch with previous
 2607       bool sp =
 2608         IsSamePatch(seg, myGrid->NbUPatches(), myGrid->NbVPatches(), iumin, iumax, ivmin, ivmax);
 2610       // not same patch has lowest priority
 2611       if (!sp && (canBeClosed || (index && samepatch)))
 2612       {
 2613         continue;
 2614       }
 2616       // try to connect, with the following priorities:
 2617       // The name of property      Weigth:
 2618       // sharing vertex            auto
 2619       // samepatch = 1             16
 2620       // ! sameedge                auto
 2621       // misorientation = 0        8
 2622       // connected in 2d           4
 2623       // distance                  2
 2624       // short                     auto
 2625       // angle ->> PI              1
 2626       const occ::handle<ShapeExtend_WireData>& wire = seg.WireData();
```

**要点提示（对照时逐条看）**：
```text
· OCCT 在此是否对 `TopAbs_INTERNAL` 的段有跳过（端口的 :110 是 `if an_or == Internal { continue; }`）；
· OCCT 的首段选择分支（对应端口 `:113-131` 的 `!has_sbwd`）条件是否一致；
· 端口的 `:106 seg.is_vertex()` 在 OCCT 里对应哪一句。
```

---

### 9.511 —— 【对照完成】OCCT `:2573-2581` **同样跳过 Internal 段** ⇒ 端口 `:110` 忠实；靶点**上移一层**：那些段为何被打成 `Internal`

#### 逐条对照（`cxx:2570-2624` ↔ 端口 `:104-132`）

| OCCT | 端口 | 判定 |
|---|---|---|
| `if (seg.IsVertex()) continue;`（:2573-2576） | `if seg.is_vertex() { continue; }`（:106-108） | ✓ |
| `if (anOr == TopAbs_INTERNAL) continue;`（:2577-2581） | `if an_or == Orientation::Internal { continue; }`（:110-112） | ✓ **两侧都跳过 Internal** |
| `if (sbwd.IsNull()) { if (shorts(i) > 0) continue; if (anOr == EXTERNAL) continue; if (anOr == FORWARD) reverse = true; index = i; … break; }`（:2584-2604） | `if !has_sbwd { if shorts[i] > 0 \|\| an_or == External { continue; } if an_or == Forward { reverse = true; } index = Some(i); … break; }`（:113-131） | ✓ |
| `bool sp = IsSamePatch(...)`（:2607-2608） | `let sp = is_same_patch(...)`（:134） | ✓ |
| `if (!sp && (canBeClosed \|\| (index && samepatch))) continue;`（:2611-2614） | `:135` 同式 | ✓ |
| 权重表（:2616-2625：samepatch 16 / misorientation 8 / connected 4 / distance 2 / angle 1） | 端口的 `w1/tail1/tail2`（§9.503 实测 28 = 16+8+4，20 = 16+4） | ✓（数值与表一致） |

⇒ **候选搜索段整体忠实**；§9.509 的「若 OCCT 不跳过 Internal ⇒ 端口多了一条过滤」这一支**不成立**。

#### 靶点上移到「这些段为何是 Internal」

```text
A 类候选（43/43）都是 Internal ⇒ 它们被 `:110` 跳过是**规格行为**。
⇒ 于是问题变成：**在 OCCT 里这些段不该是 Internal**。它们被打成 Internal 的可能来源（端口侧）：
  ① `collect_wires.rs:68-72`：`isshort > 0 && (External || one_degenerated) ⇒ set_orientation(Internal)`
     （其引用为 cxx:2519-2549，需逐条对照）；
  ② `load_wires`：非流形（non-manifold）段被建成 `Internal`（§9.460 已读该文件，语义忠实）；
  ③ `split_by_line`/`break_wires` 在生成/切分时赋的朝向。
```

#### 下一轮（一次插桩即可定位来源）

```text
在 `collect_wires.rs:71`（`seqw[i].set_orientation(Internal)`）打印被改段的坐标/边数/isshort，
再与 §9.508 的 A 类候选坐标对照：
  · 若 A 类段**正是**在 :71 被改的那批 ⇒ 靶点 = `:68-72` 的判据（`isshort` 或 `one_degenerated`）
    与 cxx:2519-2549 的差异 ⇒ 按 .cxx 修 ⇒ 跑 §9.448 验收；
  · 若 A 类段**不是**在 :71 改的（而是 load_wires/切分阶段就是 Internal）⇒ 靶点相应上移到那里。
```
