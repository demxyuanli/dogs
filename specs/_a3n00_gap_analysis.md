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

