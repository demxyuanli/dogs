# occt-geom / occt-geom2d 移植忠实度审查（2026-09-20）

## 结论摘要

只读审查 `crates/occt-geom/src`（~14.4k 行）与 `crates/occt-geom2d/src`（~4.6k 行），逐函数对照 `D:\source\OCCT-src`（tag V8_0_0）。
共列 15 条：**自创 14 条**、**等价替换未登记 1 条**、**已登记 0 条**（仓库内不存在等价替换登记表，`specs/` 下只有 `_board.md`/`_tasks.md` 等）。

最高危：**`Geom_OffsetCurve` 的 D0/D1 公式写错**（`offset.rs:16`）——OCCT 是 `P = p + Offset·(p′×Dir)/‖p′×Dir‖`（沿法向偏移），
Rust 写成 `p + Offset·Dir`（沿参考方向平移）。该类型被 STEP 读入路径直接使用（`occt-topo/src/step/read_topology.rs:1049`），
且 2D 版本法向还反号（`occt-geom2d/src/offset.rs:24` vs `Geom2d_OffsetCurveUtils.pxx:50`），偏移落在另一侧。

次高危：**求根/分类整体被采样网格 + 阈值 + 牛顿/黄金分割替代**。讽刺的是同一分支**已有忠实移植**（`extrema_pc/general_extrema_pc.rs` 逐行对上
`Extrema_GGExtPC.hxx`、`GGenExtPC.hxx`、`GFuncExtPC.hxx`），但被消费的是并存的网格版：`occt-topo/src/brep_extrema.rs:175`、
`edge_edge/edge_edge.rs:359`、`shhealing/shape_analysis_curve.rs:136`、`occt-geom/src/approx_same_parameter.rs:227`。

方法：`rg` 扫候选（`not ported|UNPORTED|simplif|approx|fallback|heuristic|sampled|iterat|clamp|Newton`）后逐条打开 OCCT
`.cxx/.hxx/.pxx` 比对分支、常量、收敛判据。注意 V8_0_0 已把 `Extrema_ExtPC`/`Extrema_ECC` 改为模板别名
（`Extrema_GGExtPC.hxx:51`、`Extrema_ECC.hxx:29` → `Extrema_GGenExtCC` → `math_GlobOptMin`），Rust 注释引用的
`Extrema_ExtPC.cxx`/`Extrema_GExtPC.hxx` 在 V8_0_0 已不存在。

> **2026-10-06 补记（T-104 / §9.588、T-106 / §9.601、§9.604–§9.613）**：本报告是 2026-09-20 基线；此后的批 91–93 与 T-104/T-106 又落地一批忠实件 —— 3D `GeomConvert::CurveToBSplineCurve` + `GeomConvert_CompCurveToBSplineCurve`（`occt-geom/src/convert_bspl.rs::curve_to_bspline_curve`，批 92），**2D 侧** `Geom2dConvert::CurveToBSplineCurve`（`occt-geom2d/src/geom2d_convert.rs`，`Geom2dConvert.cxx:181-449`）、`Geom2dConvert_CompCurveToBSplineCurve`（新文件 `occt-geom2d/src/comp_curve_to_bspline.rs`）、`Geom2dConvert_ApproxCurve`（新文件 `occt-geom2d/src/convert_approx_curve.rs`），`Geom2d_BSplineCurve` 的操作族（`bspline_curve.rs` 的 `insert_knots:247`/`insert_knot:306`/`set_origin:319`/`segment:489`/`set_knots:116`/`increase_degree:384`/`remove_knot:434`/`u_knot_index_range:210`/`update_flat_knots:230`/`bspline_copy2d:918`，`Geom2d_BSplineCurve.cxx:236-292/343/410-478/707/989`），`Geom2d_BezierCurve::segment`（`bezier_curve.rs:57`），以及 `AdvApprox` 机械件整体迁到 `crates/occt-core/src/adv_approx/`（新增 `ApproxAFunction2d`；`occt-geom/src/lib.rs:76` 改 `pub use occt_core::adv_approx;`）。发现 6 据此订正、发现 15 前半据此订正。

## 发现

### 1. `Geom_OffsetCurve` D0/D1/D2 公式错误（沿 Dir 平移，非沿法向偏移）
- 判定：自创
- 证据：`crates/occt-geom/src/offset.rs:16,18-19` — `p.coord.added(&self.direction.xyz().multiplied(self.offset))` / `fn d1(..) { let (p,d) = self.basis.d1(u); (self.d0(u), d) }`（D2 同）
- OCCT 对应：`Geom_OffsetCurveUtils.pxx:53,60-61`（`Ndir = D1.XYZ().Crossed(theDirXYZ)`，`P = p + Offset·Ndir/R`）、`:86-114`（D1 必须追加 `DNdir` 法向旋转项）
- 影响：改几何结果。`occt-topo/src/step/read_topology.rs:1049` 直接构造该类型，STEP 读入的 OFFSET_CURVE 几何全错；D1/D2 恒等于基曲线导数，切线方向错
- 建议：按 `Geom_OffsetCurveUtils.pxx:47-201` 重写 D0/D1/D2，或标未移植
- **处置（✅ 已完成：T-35 重写 D0/D1/D2；T-63 步 1 补基曲线 `EvalD3`；T-63 步 2 / 2026-09-20 第 41 轮批 19 移植 `AdjustDerivative` `pxx:322-390` 并接进 `EvalD2` 奇异支路，`isDirectionChange` 不再恒 `false`）**。余项 **T-75**：offset 类自身的 `EvalD3`（`EvaluateD3` `pxx:494-529` / `CalculateD3` `pxx:203-307`）与 `EvalDN`（`cxx:386-410`）未移植，`d3` 仍走 trait 默认零值；仓内无消费者，已在 `offset.rs` 文件头登记。

### 2. `Geom2dOffsetCurve` 法向反号 + D1/D2 缺 `DNdir`
- 判定：自创
- 证据：`crates/occt-geom2d/src/offset.rs:24,33,41` — `let nx = -d.y() / len; let ny = d.x() / len;` / `(pt, d) // first derivative same as basis` / `(pt, d1, d2)`
- OCCT 对应：`Geom2d_OffsetCurveUtils.pxx:50`（`gp_Dir2d aNormal(theD1.Y(), -theD1.X())`）、`:99`（`theD1.Add(gp_Vec2d(DNdir))`）
- 影响：改几何结果。正 offset 落到相反一侧；D1/D2 不含法向旋转项
- 建议：按 `Geom2d_OffsetCurveUtils.pxx:43-101` 直译；被 `occt-topo/src/geom_bnd_lib_offset2d.rs:101` 消费，需同步对拍

### 3. 点–曲线通用路径被换成「span/0.1 网格 + 变号 + 牛顿」
- 判定：自创
- 证据：`crates/occt-geom/src/extrema_pc/poly_roots.rs:621,644,563` — `let n = ((span / 0.1).ceil() as usize).clamp(24, 256);` / `newton_point_curve_all` / `for _ in 0..60`
- OCCT 对应：`Extrema_GGExtPC.hxx:391`（`aMaxSample = 17`）、`:424`（`mysample = max(RealToInt(aMaxSample*(sup-inf)/maxint), 3)`）、`:390-470`（按 `NbIntervals(C2)`/`DeflCurvIntervals` 分区）；`Extrema_GGenExtPC.hxx:166`（`math_FunctionRoots`）
- 影响：改几何结果。OCCT 无 `span/0.1`、无 `clamp(24,256)`、无 60 次上限
- 建议：删 p01 网格路径，调用方改走 p03（已对齐 `GGExtPC.hxx:390-470`）

- **处置（✅ 已完成：T-43，2026-09-20）**：`poly_roots` 的网格 + 60 步牛顿路径、采样分类器（`is_line`/`classify_circle`）与 `point_curve::param_for_point` 全部删除；点–曲线改走 `extrema_pc::point_curve_extrema_all` = `Extrema_ExtPC` 的真实类型分派（line/circle 走 `Extrema_ExtPElC` 解析臂，其余走 `general_extrema_pc` 的 `default:` 臂 = `Extrema_GGExtPC`）。同批按 OCCT 语义去掉了"无条件把区间端点当极值"。

### 4. 曲线–曲线被「2D 网格 + 边最佳点种子 + 16×16 兜底」替代 `math_GlobOptMin`
- 判定：自创
- 证据：`crates/occt-geom/src/extrema_cc/curve_curve.rs:83,214-239` — `let mut edges: [Option<(f64,f64)>; 4] = [None,None,None,None];` / `// Degenerate fallback: coarse grid best pair` / `let n = 16;`
- OCCT 对应：`Extrema_GGenExtCC.hxx:617`（`math_GlobOptMin aFinder(...)`）、`:639-688`（区间对 + `NCollection_CellFilter` 去重）；`Extrema_ExtCC.cxx:312-316`
- 影响：改几何结果。"边最佳点"种子与 16×16 兜底在 OCCT 中不存在
- 建议：按 `Extrema_GGenExtCC.hxx:617` 直译 GlobOptMin，或标未移植

- **处置（✅ 已完成：T-66 步 b-2，2026-09-20 第 56 轮）**：`GGenExtCC::Perform`（`Extrema_GGenExtCC.hxx:453-802`）逐段移植（区间集 / `math_GlobOptMin` / Lipschitz / 逐区间 `SetLocalParams`+`Perform` 累积最优 / `aCellSize` 去重 / 平行检测）；自创种子集 `newton_curve_curve_all`、采样分类器与 16×16 兜底**全删**，解析臂改由 `gp_line`/`gp_circ` 类型查询分派。

### 5. `IsDone`/`myDone` 被全部吞掉（失败返回兜底数值）
- 判定：自创（静默丢弃失败）
- 证据：`crates/occt-geom/src/extrema_cc/curve_curve.rs:210-241`、`extrema_ss.rs:437,494`、`extrema_surf/numeric_extrema.rs:142`、`extrema_pc/point_curve.rs:120-126` — `None => { … best.unwrap_or_else(|| pair_cc(GpPnt::zero(), 0.0, GpPnt::zero(), 0.0)) }`
- OCCT 对应：`Extrema_GGExtPC.hxx:531,545-550`（未完成 `throw StdFail_NotDone()`）；`Extrema_GGenExtCC.hxx:691-695`（`if (aNbSol == 0) { myDone = false; return; }`）
- 影响：改几何结果。曲线自交时返回正距离而非"未完成"，调用方无法区分
- 建议：删全部 `fallback_*`/黄金分割/坐标下降，按 OCCT 传 done

- **处置（◐ 主体已修）**：批 17 起 `extrema_pc::{point_curve_extrema,point_curve_max_extrema}` 与 `extrema_cc::curve_curve_extrema` 改为 `Option<ExtremaPair>`（`None` = OCCT 的 `IsDone()==false`），黄金分割兜底与 16×16 伪造对删除；**曲面族**（`extrema_surf`/`extrema_ss` 的 `fallback_*`）就地标 **UNPORTED**（OCCT 抛 `Standard_NotImplemented`/置 done=false），因有生产消费方而暂留。

### 6. `convert_bspl.rs::comp_curve_to_bspline` 把真 B 样条拼接降级为 1 次折线插值
- 判定：自创
- 证据：`crates/occt-geom/src/convert_bspl.rs:326,349` — `let n = ((b - a) / 0.1).ceil().clamp(2.0, 64.0) as usize;` / `crate::curve_reparam::resample_bspline(&pts, 1)`
- OCCT 对应：`GeomConvert_CompCurveToBSplineCurve.cxx:135-215`（`IncreaseDegree` 统一次数 + 结点/重数拼接，全程不采样折线）
- 影响：改几何结果。圆弧链 → 折线，次数/结点全变
- 建议：按 `:135-215` 直译，或标未移植

- **处置（◐ 主体已修）**：**3D** 侧 `GeomConvert_CompCurveToBSplineCurve`（`GeomConvert_CompCurveToBSplineCurve.cxx:135-215`，`IncreaseDegree` + 结点/重数拼接、不采样折线）已在批 91–92 落地并接线，`convert_bspl.rs::comp_curve_to_bspline` 的 `(b-a)/0.1` 采样与 `resample_bspline(&pts,1)` 已删除。**2D** 侧 T-104 补齐：`Geom2dConvert_CompCurveToBSplineCurve`（`occt-geom2d/src/comp_curve_to_bspline.rs`，`Geom2dConvert_CompCurveToBSplineCurve.cxx:28-244`）、`Geom2dConvert::CurveToBSplineCurve`（`occt-geom2d/src/geom2d_convert.rs`，`Geom2dConvert.cxx:181-449` 的 trimmed/非 trimmed/圆锥臂）、`Geom2dConvert_ApproxCurve`（`occt-geom2d/src/convert_approx_curve.rs`，`Geom2dConvert_ApproxCurve.cxx:31-203`）；`GeomLib.cxx:842-969` 的 `SameRange` 链由 `occt-topo` 侧接线。**仍 UNPORTED**：`Convert_Polynomial` 臂（缺 2D `BSplCLib::Trimming`）。

### 7. 3D/2D 求交与投影整体是采样器（已弃用求根）
- 判定：自创
- 证据：`crates/occt-geom/src/geom_api.rs:5-7,100-103`（`//! ponytail: sampler-based approximations throughout`、`na=256/nb=256`、`coarse = tol*100.0 + 1e-9`）；`crates/occt-geom2d/src/curve_ops.rs:68-69,82,104,117`；`geom2d_api.rs:19-44,163-199`（8 次牛顿越界 break；64×64+黄金分割，无界直接 `None`）
- OCCT 对应：`GeomAPI_ProjectPointOnCurve.cxx:65`（`Extrema_ExtPC`）；`Geom2dAPI_InterCurveCurve.cxx:66,108`（`Geom2dInt_GInter` + `NbSegments`）；`IntAna2d_AnaIntersection_1.cxx:22,77-90`（lin/lin 的 `iden`/`para`）；`GeomLib_Tool.cxx:157-184`
- 影响：改几何结果。相切/重合漏解；重合直线在 OCCT 是 `iden`（一段），Rust 既非点也非段；投影的 `tol` 被忽略（`_tol`）
- 建议：改调本仓库 extrema 端口或按 `IntAna2d_AnaIntersection_1.cxx:22-129` 直译

- **处置（✅ 已完成：批 81 / R2-2，2026-09-21）**：见文末「其他已核实条目」最后一条 —— 求交 half 的自创采样器已整体摘除，生产调用点改走忠实件；投影 half 在 T-52 已改走 `Extrema_ExtPC`/`ExtPS`。

### 8. `intana` 求交族：四次方程族自创阈值 + 锥/锥阈值 + 未移植返回 None
- 判定：自创（多处）+ 未移植
- 证据：`crates/occt-geom/src/intana/analytic_intersections.rs:120,145,152`（`two_m_p < 1e-12` / `poly4(..).abs() < 1e-6` / 去重 `1e-8`；`:116` 失败返空根且 done 仍真）；`intana/line_torus.rs:336,340,370`（`1e-14` 绝对阈值替代自适应 `aTolAng`）；`intana/line_torus.rs:396,490`、`intana/analytic_intersections.rs:519`（直接 `None`）
- OCCT 对应：`math_DirectPolynomialRoots.cxx:344-393`（`ShouldReduceDegreeQuartic`）、`:583-587`（`if (!aSuccess) { myDone = false; return; }`）、`:678-705`、`:92-120`（`RefineRoot`）；`IntAna_QuadQuadGeo.cxx:1440,1458-1475,1524`（`TOL_APEX_CONF=1e-10`、`EstimDist` 求 `aTolAng`）；`:1324,1603`、`:1060-1200`（未移植的完整分类）
- 影响：改几何结果（近轴近平行锥漏解；本应 not-done 被当成"无解曲线"）
- 建议：按 `math_DirectPolynomialRoots.cxx` 重译并传 done；补 `EstimDist`/`aTolAng`

### 9. 直线–环面：物理残差 `1e-7` 取代参数回代，且无法表达 `done=false`
- 判定：自创（静默失败）
- 证据：`crates/occt-geom/src/intana/line_torus.rs:574-576` — `let err = ((rho - r)*(rho - r) + dz*dz - rr2).abs(); if err < 1e-7 { if seen.iter().all(|s| s.distance(&p) > 1e-7) {`
- OCCT 对应：`IntAna_IntLinTorus.cxx:99-103`（回代 `ElSLib::Parameters/Value`，`if (a0 > 0.0000000001) aNbBadSol++;`）、`:116-120`（全解被拒 → `nbpt = 0; done = false;`）
- 影响：改几何结果 + 全解被拒时返回空 `Vec`，调用方误判"无交点"
- 建议：改 `Result`/done，按 `:99-125` 直译

### 10. 类型判定：采样几何不变量取代 `GetType()`
- 判定：自创
- 证据：`crates/occt-geom/src/extrema_pc/poly_roots.rs:434,452`（6 点切向 + `1e-7*m0*m`）、`:505,515`（共面 `1e-4*nmag*ni`、等距 `1e-4*r`）；复制到 `extrema_cc/poly_roots.rs:476,523`、`extrema_surf/analytic_solvers.rs:537,567`、`extrema2d/analytic_solvers.rs:459,497`
- OCCT 对应：`Extrema_GGExtPC.hxx:123,160-181`（`type = TheCurveTool::GetType(theC)` + switch）；`Extrema_ExtCC.cxx:188-189`、`Extrema_ExtCC2d.cxx:105`、`Extrema_ExtSS.cxx:122-127`
- 影响：改几何结果（近似圆/浅弧可能进错分支）
- 建议：`Curve` trait 已有 `is_line()/gp_circ()/gp_ellipse()/circle_radius()`（`curve.rs:34-50`），直接用于分派并删采样分类器

- **处置（✅ 已完成）**：`Curve` trait 的 `is_line()`/`gp_circ()`/`gp_ellipse()` 类型查询已作为唯一分派依据；`classify_curve`（STEP 曲线族采样器）随 T-39 整体删除，`brep_surface::classify_surface` 随 T-72 改为 `GeomAdaptor_Surface::Load` 的精确类判定（详见 `_index.md` A3 / A24 / A30）。

### 11. 解析圆锥曲线的 D2 未移植：双曲线/抛物线静默返回零向量
- 判定：等价替换未登记（未移植臂静默给默认值）
- 证据：`crates/occt-geom/src/hyperbola.rs:16`、`parabola.rs:16` — `(self.d0(u), self.d1(u).1, GpVec::zero())`
- OCCT 对应：`Geom_Hyperbola.cxx:244`（`ElCLib::HyperbolaD2`，`:253` 还有 D3）；`Geom_Parabola.cxx:186`（`ElCLib::ParabolaD2`）
- 影响：改几何结果。抛物线 D2 恒为 `Yd/(2f)` 而非 0；消费方 `approx_same_parameter.rs:85`、`occt-core/src/gcpnts_perform.rs:384`、`extrema_cc/poly_roots.rs:590` 的曲率/逼近/求交全被污染
- 建议：在 `occt-core/src/elib/clib.rs` 补 `hyperbola_d2/parabola_d2`（对照 `ElCLib.cxx`）后替换

- **处置（✅ 已完成：T-63，2026-09-20）**：`occt-core/src/elib/clib.rs` 已补 `hyperbola_d2`/`parabola_d2`（`ElCLib.cxx:378-401` / `:404-431`）并接线，`hyperbola.rs`/`parabola.rs` 的零 `d2` 不复存在。

### 12. `Surface::d2` 默认用一阶导前向差分，解析曲面全走这条自创路径
- 判定：自创
- 证据：`crates/occt-geom/src/surface.rs:12-20`（`let h = 1e-6; … let d2u = du_u.subtracted(&du).divided(h);`）；`sphere.rs:18-29`、`cylinder/cone/torus/plane.rs` 均未重写 d2
- OCCT 对应：`Geom_SphericalSurface.cxx:230`、`Geom_Plane.cxx:214`、`Geom_CylindricalSurface.cxx:232`、`Geom_ToroidalSurface.cxx:241`、`Geom_ConicalSurface.cxx:277`（均 `EvalD2` 解析）
- 影响：改几何结果。`occt-topo/src/meshing/geomlib_norm.rs:30` 正是 `GeomLib::NormEstim` 的移植（OCCT `GeomLib.cxx:2592` 用解析 D2 交 `CSLib::Normal`），差分 D2 在极点/奇点附近改变法向判定
- 建议：给 5 个解析曲面补 ElSLib 版 d2，或把默认 d2 改为显式 panic 防静默近似

- **处置（✅ 已完成：T-51，2026-09-20）**：五个初等面（plane/cylinder/cone/sphere/torus）已改调解析 `surface_eval::{plane,cylinder,cone,sphere,torus}_d2`（`ElSLib::*D2`）；trait 默认的 `h=1e-6` 差分就地标 **UNPORTED**，仍继承它的 6 个实现者（Bezier/Offset/SurfaceOfRevolution/LinearExtrusion/`surface_fit` 等）已在注释列出。

### 13. 弧长用自适应 Simpson + 自创容差，OCCT 是 Gauss 积分
- 判定：自创
- 证据：`crates/occt-geom/src/gcpnts.rs:16,25-70,100`（`MAX_DEPTH=20`、`tol.max(span * 1e-10)`、`let tol = CONFUSION * 0.1;`）、`:218-230`（`quasi_uniform_abscissa` 对所有曲线用 `2*n` 折线表）
- OCCT 对应：`CPnts_AbscissaPoint.cxx:178`（`math_GaussSingleIntegration(FG, U1, U2, order(C), Tol)`）、`GCPnts_AbscissaPoint.cxx:441`（容差 `theC.Resolution(Precision::Confusion())`）；`GCPnts_QuasiUniformAbscissa.cxx:117-128`（非 Bezier/BSpline 委托 `GCPnts_UniformAbscissa`）、`GCPnts_AbscissaPoint.cxx:77-90`（Line/Circle 闭式）
- 影响：改几何结果（所有弧长/等弧长采样值）
- 建议：按 `CPnts_AbscissaPoint.cxx` 直译 Gauss 积分，补类型分派与闭式臂
- **处置（2026-10-06 复核：仍未修）**：`occt-geom/src/gcpnts.rs` 的自适应 Simpson 与 `CONFUSION*0.1` 容差仍在（`gcpnts.rs:16,25-70,100`），忠实件 `CPnts_AbscissaPoint`/`GCPnts_AbscissaPoint`（`CPnts_AbscissaPoint.cxx:178` 的 `math_GaussSingleIntegration` + Line/Circle 闭式 + 类型分派）**尚未移植**；该项与 A15 的「`gcpnts.rs` 积分」同条，仍为 live 缺口。

### 14. `curve_reparam.rs` 用 1e-15/1e-14/0.0 判据做 de Boor
- 判定：自创
- 证据：`crates/occt-geom/src/curve_reparam.rs:37`（`if den.abs() <= 0.0 { return self.new_a; }`）、`:189`（`if denom.abs() > 1e-15`）、`:228`（`if best < 1e-14 { return None; }`）
- OCCT 对应：`BSplCLib::Deboor`/`Geom_BSplineCurve::Reparametrize` 无这些阈值；OCCT 分母由结点序列保证非零
- 影响：改几何结果（阈值内被静默置零 → 曲率/参数映射偏差）
- 建议：按 `BSplCLib` 直译，删阈值
- **处置（✅ 已完成：第 51 轮）**：`curve_reparam.rs` 的 `1e-15`/`1e-14`/`0.0` 判据已删除，改走忠实 `BSplCLib::EvalBsplineBasis`（`BSplCLib_2.cxx:429-563`），判据为 `gp::Resolution()`；**订正**：OCCT 8.0.0 全树**无** `BSplCLib::BasisFuns`，老名字已移除。

### 15. `Geom2dBSplineCurve` 权重/周期**已补齐**（◐ 订正）；`Geom_RectangularTrimmedSurface` 仍缺 `SetTrim` 控制流
- 判定：自创（缺失控制流）；**前半已于 T-104/T-106 订正**
- 证据（**订正后**）：~~`crates/occt-geom2d/src/bspline_curve.rs:9-14,42-64`（仅 `xs/ys/knots/degree`）~~ —— 现该结构已带 `weights: Option<Vec<f64>>`（`bspline_curve.rs:16`）、`periodic: bool`（`:32`）、`set_periodic`（`:139`）与 `weights` 变更传播（`:165,172`），`Geom2d_BSplineCurve.cxx:136-178` 的 `myWeights`/`myPeriodic` 已可表达。**仍缺**：`crates/occt-geom/src/rectangular_trimmed.rs:54-61`（`uv`）/`:63-73`（`uv_raw`）直接存参，无 `ElCLib::AdjustPeriodic`、无 U1>U2/Sense 反转、无越界检查；`rectangular_trimmed.rs` 内无 `set_trim`/`SetUTrim`/`SetVTrim` 符号。
- OCCT 对应：`Geom2d_BSplineCurve.cxx:136-178`（`myWeights`/`myPeriodic` + `CheckCurveData`，**已对齐**）；`Geom_RectangularTrimmedSurface.cxx:210-336`（`:237` U1==U2 抛异常、`:248` `ElCLib::AdjustPeriodic`、`:269-273` 越界 `Standard_ConstructionError`、`:328-335` `UReverse/VReverse`）
- 影响：U1>U2 或 Sense=false 时 OCCT 反转基面参数方向（法向/du 符号翻转），Rust 原样存参
- 建议：按 `cxx:210-336` 补 `SetTrim`/`SetUTrim`/`SetVTrim` 全部分支与 `AdjustPeriodic`/越界检查；`rectangular_trimmed.rs:119-160` 的 D0/D1/D2 转发已与 `cxx:391-407` 一致，不用动
- **处置（2026-10-06 订正）**：前半（2D 样条权重/周期）**不再是缺口**；后半（`SetTrim` 族）**仍为 UNPORTED**，需按 `Geom_RectangularTrimmedSurface.cxx:210-336` 移植。
## 已核对为忠实移植（反证，供"未发现"依据）

- `extrema_pc/general_extrema_pc.rs` 整体：BSpline 结点臂 `:476-611` ↔ `Extrema_GGExtPC.hxx:190-388`；default 臂 `:699-774` ↔ `:390-470`；`AddSol :347` ↔ `:617-632`；`SearchOfTolerance`/`MaxOrder`/奇异 DN/三点 DF ↔ `Extrema_GFuncExtPC.hxx`；`FunctionRoots :366-375` ↔ `Extrema_GGenExtPC.hxx:166`；`defl_curv_intervals :632-695` ↔ `Extrema_CurveTool.cxx:41-91`
- `extrema_pc/poly_roots.rs:263-295` 点–圆 ↔ `Extrema_ExtPElC.cxx:125-190`；`extrema_cc/poly_roots.rs:257-277` 线–线 ↔ `Extrema_ExtElC.cxx:327-338`；`:285-338` 线–圆系数 `A1..A5` ↔ `:560-564`（逐字一致）；`circle_circle_extrema ↔ :982-1110`
- `approx_same_parameter.rs:703-737` 后处理兜底 ↔ `Approx_SameParameter.cxx:482-539`（含 `11/40`、`U/VResolution`、`anApproxTol < myTolReached` 比较）
- `osculating_surface.rs:54-78` 十等分采样 ↔ `Geom_OsculatingSurface.cxx:784-832`；`build_osculating :308-467` ↔ `cxx:532-780`（无曲率阈值类启发式）
- `intana_curve.rs` 系数/容差 ↔ `IntAna_Curve.cxx:95-217,279-374,440-568`；`intana_intquadquad.rs:139-193` ↔ `IntAna_IntQuadQuad.cxx:61-118`；`intana_torus.rs` 常量 ↔ `IntAna_QuadQuadGeo.cxx:351-359,2539-2560,2606`；`intana_trig.rs` 主体 ↔ `math_TrigonometricFunctionRoots.cxx:75-504`
- `occt-geom2d/trimmed.rs:8-96` ↔ `Geom2d_TrimmedCurve.cxx:98-154`；`occt-geom2d/offset.rs:48` ↔ `Geom2d_OffsetCurve.cxx:414-418`；2D 椭圆/双曲/抛物线 D1/D2 ↔ ElCLib；`plane/cylinder/cone/torus.rs` D0/D1 与 UIso/VIso ↔ `ElSLib.cxx:1705-1810`
- `bspline2d_to_bezier.rs:126-132` 的 `clamp(0,1)` 仅为守卫，与 `BSplCLib::InsertKnots` 等价、无数值后果

## 未覆盖

`bspline_surface.rs`（de Boor/有理求值内部）、`surface_of_revolution.rs`、`surface_to_grid.rs`、`transform.rs`、`bezier_curve.rs`、
`offset_surface_utils.rs` 的其余分支未逐行穷举，建议下一轮复查。**订正（2026-10-06）**：`bezier_curve.rs` 已不再是空白 —— `Geom2d_BezierCurve::segment` 已移植（`bezier_curve.rs:57`，`Geom2d_BezierCurve.cxx:356-391`）并被 `Geom2dConvert` 链使用。**订正（A31，已修）**：`revolved.rs` 是伪造的回转面类（`d1` 把 basis 点坐标当 `dv`），**已整文件删除**（`crates/occt-geom/src/revolved.rs` 不存在），正确实现为 `surface_of_revolution.rs`（`Geom_SurfaceOfRevolution`）。

## 其他已核实条目（超出 15 条上限，按需处理）

- `osculating_surface.rs:99-107,136-143`：失败时 `clear_flags(); return;`，OCCT `Geom_OsculatingSurface.cxx:206-213` 只清标志不 return（`:245` 继续建 `myOsculSurf2`）→ 自创（提前退出）
- `interp_curve.rs:211,96,252-254`：均匀结点 + 有限差分导数基（`h=1e-7`）+ 逐点切矢，OCCT `GeomAPI_Interpolate.cxx:622-646,664-669,795` 用弦长参数与 `BSplCLib::Interpolate` → 自创（零生产调用点）
- `crates/occt-core/src/adv_approx/approx.rs:257-258`（原 `occt-geom/src/adv_approx/`，T-104 已整体迁到 `occt-core`）：`DichoCutting`（`tmil=0.5*(a+b)` + `20.0*PCONFUSION`）**已不再是唯一选项** —— `adv_approx/cutting.rs:40` 的 `PrefCutting`（`AdvApprox_PrefCutting.cxx:20-44`）与 `:75` 的 `PrefAndRec`（`AdvApprox_PrefAndRec.cxx:21-75`，`Weight<=1` 抛 `Standard_DomainError`）**均已移植**，且 `Geom2dConvert_ApproxCurve` 已按 `Geom2dConvert_ApproxCurve.cxx:145-151` 传入推荐/优先切分点（`occt-geom2d/src/convert_approx_curve.rs`）。**订正**：本条原称"缺 `AdvApprox_PrefAndRec`"已不成立；`CutPnts_C2/C3` 在 `cutting.rs` 中无同名符号。生产路径 `offset_surface.rs:130` 是否已切到 `PrefAndRec` **未核**。
- `surface_fit.rs:13-40`：`fit_plane` 无 `GeomLib_IsPlanarSurface.cxx:44-45` 的 `gz < Tol` 门，调用点 `occt-topo/src/brep_offset/curve_face_offset.rs:491-500` 每边仅采 4 点；`fit_sphere` 无 OCCT 对应 → 等价替换未登记
- `curve_approx.rs:23-25,49-55`：无界曲线回退 `curve_to_polyline_uniform(c,64)`/代理区间 `(0,1)`，OCCT `GCPnts_UniformDeflection` 无此分支；`bezier_to_polyline.rs:7,48-68` 的 `MAX_DEPTH=24` 与内部极点平坦度判据 OCCT 中不存在 → 自创（仅 `brepmesh.rs:385` 调用）
- `offset_surface.rs:160-164,185,209-212` 与 `projlib.rs:102`：前者 `transform` 空实现、`EvaluateD0/D1` 失败静默返回基面点/向量（`:63` 兜底写死 `(0,0,1)`），OCCT `Geom_OffsetSurface.cxx:338-343,839-842` 抛 `Geom_UndefinedValue`/执行 `Transform`；后者固定 32 采样而 OCCT `ProjLib_ProjectedCurve.cxx:639-659` 走自适应 `Approx_CurveOnSurface` → 均自创
- `approx_same_parameter.rs:436,442` 的 `UNPORTED: Geom_BezierSurface::Resolution` 是真实缺口：OCCT `GeomAdaptor_Surface.cxx:1872-1875,1934-1937` 有该分支，Rust 走 `default: Precision::Parametric`，U/V 容差不同 → 改几何结果，需登记
- `geom_api.rs` 采样器**已整体摘除（订正，2026-10-06）**：`curve_surface_intersections`/`curve_curve_intersections`/`curve_curve_distance`/`surface_bbox`/`minimize_1d`/`dist_curve_surface`/`curve_window_for_bbox`/`sample_curve` 已删除（A16 求交 half，批 81 / R2-2），生产调用点改走忠实件（曲线–曲线用 `edge_edge::EdgeEdge`、边–面用 `intcurvesurface::perform_curve_surface`）；`geom_api.rs` 现存函数为 `project_point_on_curve/surface`（忠实 `Extrema_ExtPC`/`ExtPS`）、`distance_point_*`、`tangent_at`、`curve_bbox`、`pcurve_of_curve_on_surface`（`:114`/`:133` 的 `samples` 参数仅为诊断采样，非求交）。**2D 余项 R2-18 仍在**：`occt-geom2d/curve_ops.rs` 的 `curve2d_intersections`/`curve2d_closest_point` 采样族仍被 `geom2d_api` 生产使用。
