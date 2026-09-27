# a3n00 网格/几何缺口分析（2026-09-26）

> 症状（用户提供）：`data/occ/a3n00.stp` 经 OCCT 渲染正确（法兰盘 + 螺栓孔 + 大曲率阀体），
> 端口 `output/a3n00.obj` 渲染时**法兰与大曲率面破碎/扇贝状**。
> 本文只做**定位与取证**，不改代码；结论均来自本机实测（OCCT 8.0.0p1 探针 + 端口探针）。
> 相关画板条目：`specs/_board.md` §3.3 T-69/T-91/T-92、§2 门禁、`specs/_occt_mesh_gt.md`。

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
































