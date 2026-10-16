/** @canvas
 * title: OCCT→Rust 移植范围与偏移审查
 * description: 全库范围量测（OCCT vs Rust）、语义/文档偏移残余、本轮偏移收敛读数、门禁现状与下一步
 * icon: search
 */
import {
  BarChart,
  Callout,
  Divider,
  Grid,
  H1,
  H2,
  H3,
  Pill,
  Row,
  Stack,
  Stat,
  Table,
  Text,
} from "dsh/canvas";

// ---------------------------------------------------------------------------
// 量测口径（2026-10-16 实测；OCCT 侧与 10-13 / 10-11 逐位相同）
//   OCCT 侧: D:\source\OCCT-src\src  .hxx/.cxx/.lxx/.gxx
//   Rust 侧: crates/{occt-core,occt-math,occt-geom,occt-geom2d,occt-topo}/src  .rs
//   基线  : docs/goals/drift-and-gap-autofix/state.yaml (board truth) · specs/_audit/_index.md (A0-A27 各行「处置」,
//           逐包同步; 最新的 P033 记录在 specs/_audit/occt-geom.md §9)
//           （旧任务看板 specs/board.canvas.tsx 已于 2026-10-16 退役删除，本卡与其门禁表是唯一维护中的板）
//   对拍  : crates/occt-topo --test step_obj_gates（`--nocapture` 逐模型 v/f 与面积比）
// ---------------------------------------------------------------------------
const MEASURE = {
  asOf: "2026-10-16",
  occt: { files: 13803, lines: 2439487 },
  // rust/crates: `src/**/*.rs` only (tests/ and examples/ excluded), matching the
  // `git ls-files "crates/*/src/*.rs"` file set (1111 files at P034).
  rust: { files: 1111, lines: 338784 },
  cadKernel: 1397840, // Foundation 252505 + ModelingData 347970 + ModelingAlgorithms 797365
  // case-sensitive `UNPORTED` occurrences over `crates/**/*.rs` (176 case-insensitive).
  unported: 155,
};

const GROUPS = [
  { group: "FoundationClasses", files: 1057, lines: 252505, rust: "core + math + geom/geom2d", depth: "深" },
  { group: "ModelingData", files: 1431, lines: 347970, rust: "geom/2d + topo(BRep/Poly)", depth: "深—中" },
  { group: "ModelingAlgorithms", files: 3247, lines: 797365, rust: "topo(BOPAlgo/BRepMesh/IntTools/ShHealing)", depth: "中" },
  { group: "DataExchange", files: 4699, lines: 524997, rust: "topo(step/iges/rwmesh/gltf/bincaf/xcaf)", depth: "浅—中" },
  { group: "Visualization", files: 1103, lines: 209817, rust: "topo(viz_scene/render_svg/hlr)", depth: "lite 替代" },
  { group: "Draw", files: 429, lines: 172420, rust: "topo(draw)", depth: "lite 替代" },
  { group: "ApplicationFramework", files: 865, lines: 105786, rust: "topo(xcaf/xmlcaf)", depth: "lite 替代" },
  { group: "Deprecated", files: 972, lines: 28627, rust: "—", depth: "未移植" },
];

const CRATES = [
  { crate: "occt-topo", files: 638, lines: 213964 },
  { crate: "occt-core", files: 220, lines: 45107 },
  { crate: "occt-geom", files: 133, lines: 45074 },
  { crate: "occt-geom2d", files: 63, lines: 19747 },
  { crate: "occt-math", files: 57, lines: 14892 },
];

// 本轮已收敛的偏移（均已跑现基线核实；行号为当前工作树实测）
const CLOSED = [
  {
    item: "缝边（seam）只注册单条 pcurve，且锥面缝未按 `FirstParameter` 交换",
    id: "T-124",
    where: "topo/meshing/node_insertion.rs:init_data_structure · topo/meshing/model_builder/preprocessor.rs:split_cone_seam",
    evidence:
      "OCCT 建链用边的**全部** pcurve 表示（`BRep_Tool::CurveOnSurface` 的各 slot），端口原先只按 `CheckOrder` 槽位取一条；锥面缝另需按 `FirstParameter()` 交换（`BRepMesh_ModelHealer` 的 `splitConeSeam` 语义，`FirstParameter` 在 `-Precision::INFINITE` 处求值）。接线后 TDB 三角 81,721 → **79,926**（OCCT 79,977，Δ −51），面积比 **1.0042 → 1.0000**；逐面配对 2,050 面与 OCCT 逐位相同，OCCT 未成网格的 12 面在端口同样为 0（提交 `bc7249fa`）。",
  },
  {
    item: "STEP 读入缺 `ShapeFix_Shape::Perform` 尾段（`SameParameter` + `FixVertexTolerance`）",
    id: "T-93",
    where: "topo/step/read_topology/resolver_context.rs:96 close_heal_context · shhealing/shape_fix_same_parameter.rs",
    evidence:
      "尾段按 `ShapeFix_Shape.cxx:259-261`（`SameParameter(sh,false,0.0)`，`NeedFix(-1)` 为真）与 `:267-291`（`FixVertexTolerance`）接线，再走 `:294` 的收口 `Apply`。根因：尾段缺失时环面 u-iso 边留在 `same_param=false`，`BRepMesh_EdgeDiscret::CreateEdgeTessellator`（`cxx:57-59`）改取 `Adaptor3d_CurveOnSurface` 重载，`EvalKPart` 判出 `GeomAbs_Circle` 后走 `PerformCircular` 的定角步长（17 个 V 参数）而 OCCT 自适应 `PerformCurve` 只给 11（face 4: nvp 28→11、nodes 525→315、mesh 629/1152→407/720，与 OCCT 探针逐位一致）。全模型 T0M 三角 73,150 → 67,424（OCCT 67,366，缺口 +5,784 → +58），面积比 1.0000。",
  },
  {
    item: "面网格走 UV 栅格 + 耳切/桥洞（非 OCCT 网格器）",
    id: "A13 / A18 / A23",
    where: "topo/wireframe.rs:110（`face_to_triangles` 已改委托 `IncrementalMesh::discretize_face`）· `face_uv_bounds` 委托 `brep_tools::uv_bounds`",
    evidence:
      "UV 栅格采样、`ear_clip`/`bridge_holes`/`planar_polygon_triangulate`、`WIREFRAME_FALLBACK_RATIO_MAX` 失败率阈值、`legacy` 四叉树两级兜底全部删除（2026-10-12，提交 `93b076b4`）；平面面片改走同一套 `BRepMesh_Delaun`。同批删除 `topo/brepmesh.rs`（−602 行）。",
  },
  {
    item: "共享 `ShapeBuild_ReShape` 被开关默认关闭（重合边 ⇒ T 型接点）",
    id: "A33",
    where: "topo/step/read_topology/resolver_context.rs",
    evidence:
      "删除端口私有的 `OCCT_SHARED_HEAL` 开关，共享 `ReShape` 成为唯一路径。T0M 开放边 3,643 → 250、a3n00 571 → 0、acs10 646 → 4。",
  },
  {
    item: "网格/几何自创件摘除（offset 公式、`classify_*` 采样族、`gcpnts` 二分族、`extrema` 采样族、IGES/STEP 写侧采样重拟）",
    id: "A0 / A3 / A6 / A7 / A20-A24 / A26",
    where: "见 specs/_audit/_index.md 各行「处置」",
    evidence:
      "均已改为忠实件并跑现基线（`step_obj_gates` 4/4、`export_data_obj` 23 ok / 0 err、逐面回退计数 0）。",
  },
  {
    item: "IGES 写侧最后两格：周期面反周期化（`SetUOrigin`/`SetVOrigin` 重定原点 + `SetUNotPeriodic`）与全周椭圆支（104→126）",
    id: "A26 收口",
    where: "topo/iges（`u_reorigin`/`v_reorigin` 块、`emit_whole_period_ellipse`、`emit_circle_transfer`）· geom/bspline_surface.rs::set_u_origin/set_v_origin · core/bspl/unperiodize.rs",
    evidence:
      "对拍现有基线：`iges_check` 14/14 `ok`；参数段记录多重集 vs `_igesdump/oracle/*.iges` 10/14 逐字节一致（Shape-2 的 28 条 128 记录逐字节一致），余 3 个模型为浮点噪声/读侧极点密度，`OffsetPlaneHoleEdge` 仅缺 314/406（颜色/名字）。`iges/emit_surface.rs` 内「`SetUOrigin` 未移植」旧注释已清除（2026-10-13）。",
  },
  {
    item: "`Geom_RectangularTrimmedSurface::SetTrim` 族（周期折叠 / Sense 反转 / 越界抛错 / `UReverse`·`VReverse`）",
    id: "A15 姊妹项",
    where: "geom/rectangular_trimmed.rs · surface.rs:rectangular_trimmed_spec · plane/cylinder/cone/sphere/torus/bezier/bspline/offset/revolution/extrusion 的 `UReversedParameter`·`VReversedParameter`",
    evidence:
      "按 `Geom_RectangularTrimmedSurface.cxx:210-378` 全量移植：`SetTrim` 8 参版（`U1==U2`/`V1==V2` 抛错、周期向 `ElCLib::AdjustPeriodic(udeb,ufin,min(|u2-u1|/2,PConfusion),u1,u2)`、非周期向 `U1<U2` 取/反 Sense 与 `PConfusion` 越界检查）、6 参版转发、4 参单方向版转发、`UReverse`/`VReverse`（基面 `UReversedParameter(utrim2)`/`(utrim1)` 映射后重走 `SetTrim`）与 `Geom_Surface::UReversed`/`VReversed`（`Copy()` + `U/VReverse`）。各具体面按 `Geom_*.cxx` 补 `U/VReversedParameter`（Plane `-u`；Cylinder/Cone/Sphere `2π-u`；Torus `2π-u`、`2π-v`；Bezier `1-u`；BSpline `Flat(1)+Flat(Upper)-u`；Revolution `2π-u`、基曲线 `ReversedParameter`；Extrusion 基曲线 `ReversedParameter`、`-v`；Offset 转发基面）。STEP 读侧 `RECTANGULAR_TRIMMED_SURFACE` 取 `Usense`/`Vsense` 走 `uv_sense`；`Adaptor3d_CurveOnSurface::EvalFirstLastSurf` 的 `UTrim`/`VTrim` 是 adaptor 窗口（`GeomAdaptor_Surface.cxx:411-427` 原样存窗口、不折叠），故另设 `window()` 供 `cos_locate` 使用。`step_obj_gates` 4/4、`data/occ` 八模型面积比 1.0000（2026-10-13 复核）。",
  },
  {
    item: "偏移面分支用采样判面型（`is_planar` 8×8 + 外接圆心估 `sphere_center`）代替 downcast",
    id: "发现 14",
    where: "topo/brep_surface.rs:107,112,136 · topo/brep_offset/curve_face_offset.rs:517 · topo/edge_face_kind.rs:34,59,207,217,271,284 · topo/inttools_range.rs · topo/algo_tools_range.rs · topo/bean_face.rs · topo/wire_splitter_block.rs · topo/int_tools_vertex_line.rs · geom2d/curve.rs",
    evidence:
      "（已全部收口，2026-10-15）(1) `sphere_center` 改 `GeomAdaptor_Surface::Sphere()` 等价（`s.gp_sphere()`，`GeomAdaptor_Surface.cxx:1997-2006`；`gp_sphere` 已按 `cxx:423-425` 递归裁剪面基面，`Geom_OffsetSurface` 给 `GeomAbs_OffsetSurface` 不作答，`cxx:496-523`），`det3`/`solve3`/`circumcenter` 采样族整段删除；`offset_face` 平面臂改 `gp_pln()`（法向取 `gp_Pln::Position().Direction()`，与 `BRepOffset.cxx:54-61` 的 `XDirection ^ YDirection * Offset` 同源）、球臂改 `gp_sphere()` 并直接取 `Radius()`。原实现把**几何上平的 B 样条面**判成平面臂（OCCT 给 `GeomAbs_BSplineSurface`，落一般臂），且 `sphere_center` 的 `Some` 兜底把**锥/柱面**拖进球臂做一次以伪圆心为中心的 `SetScale`，面边界被错误缩放。(2) 新增 `is_plane_type`/`face_is_plane_type`（`brep_surface.rs:107,112`，即 `GetType() == GeomAbs_Plane`）后分两批替换采样调用点：批次 2（2026-10-15）换 BOP 交叉与 EF 家族 9 处（`pave_ff_perform` / `pave_intersect/face_face` / `pave_ef_perform` / `pave_force_ef` / `pave_pcurves` / `edge_face`），并把 `is_plane_ff` 按 `_6.cxx:84-102` 重写为真正的 `IsPlaneFF`；批次 3（T-129）换完余 19 处（`face_plane`、`pcurve.rs:249`、`fillet_curved`、`gprop_analytic`、`brep_offset/shell_offset.rs`、`brep_faces` 及测试夹具），并把 `brep_surface::is_planar`/`face_is_planar` 的 8×8 采样体整段删除，全仓已无该符号。(3) 批次 3 同族把 `edge_face_kind.rs` 的曲线/曲面族也从采样换成 `GetType()` 同等 downcast：`curve_kind`＝`gp_line()`/`gp_circ()`/`nurbs_degree()`（与 `GeomAdaptor_Curve::load` `GeomAdaptor_Curve.cxx:252-311` 同源）、`circle_geometry`＝`CurveTool::Circle`（`IntTools_EdgeFace.cxx:800,829`）、`surface_kind`→`classify_surface`、`is_eq_distance`＝`cxx:240-299`（frame 取 `gp_cylinder()`/`gp_cone()`/`gp_torus()`，锥顶按 `gp_Cone::Apex()`）、`is_coplanar`/`is_radius`＝`cxx:787-811`/`cxx:816-840`；采样件 `is_line_like`/`is_circle_like`/`det3`/`solve3`/`circumcenter`/`cylinder_matches`/`cone_matches`/`torus_matches`/`quadric_u0`/`midpoint` 删除。门禁：`--lib` 1251/0、`--tests` 全 `ok`（`step_obj_gates` 4/4、`brep_offset` 12/12）、`export_data_obj` 15/15 且 v/f 逐位一致（`Sphere` 642/1244、`Torus` 1369/2592、`Offset` 712/892、`Shape` 6150/11372）。(4) 批次 4（T-130，2026-10-15）把 `Resolution` 采样族整体换成 `Adaptor::Resolution`：`occt-geom2d` 新增 `Curve2d::resolution`＝`Geom2dAdaptor_Curve::Resolution`（`Geom2dAdaptor_Curve.cxx:1186-1219`，Line/Circle/Ellipse/Bezier/BSpline/Trimmed 逐型，default `Precision::Parametric`）；`edge_face_kind::curve_resolution`(`:73`)、`inttools_range::curve_resolution`（3 处）、`algo_tools_range` 的 `curve_is_line`+`curve_is_freeform`+`curve_resolution`、`bean_face::resolution`、`int_tools_vertex_line::adaptor_resolution`、`wire_splitter_block::angle_2d` 全部改走解析闭式（`IntTools_EdgeFace.cxx:326-327,373`、`IntTools_ShrunkRange.cxx:162`、`BRepLib_1.cxx:61,202`、`BOPTools_AlgoTools_2.cxx:299-341/385-415`、`IntTools_BeanFaceIntersector.cxx:178,200,223,583`、`BOPAlgo_WireSplitter_1.cxx:792,795`、`BOPAlgo_PaveFiller_6.cxx:3001-3002`），Circle/Ellipse/Parabola/Hyperbola 从 D1 臂放回 `Resolution` 臂；`ReparamCurve2d`/`PlaneKeepParam2d` 两个装饰器补 `resolution` 以免退回 trait default。门禁：`occt-geom2d` 72/0、`--lib` 1251/0、`--tests` 1317/0（`step_obj_gates` 4/4）、`export_data_obj` 15/15 且 v/f 逐位一致。**同批残留（未改）**：`inttools_range::find_nearest_valid_point` 仍缺 `BRepLib_1.cxx:69-83` 的 `aD1Mag` 局部奇点加速分支（原注释已标 `ponytail:`）。**该残留已于 T-133（2026-10-15）收口**（见 B5 / T-133 行）：阈值与 `D1` 加速倍增循环落地，判型复用 `edge_edge::curve_kind`＝`GeomAdaptor_Curve::GetType()`（`GeomAdaptor_Curve.cxx:252-311`），`ponytail:` 注记删除。",
  },
  {
    item: "`findNearestValidPoint` 缺 `aD1Mag` 局部奇点加速臂（T001 的 B5）",
    id: "B5 / T-133",
    where: "topo/inttools_range.rs:295 `find_nearest_valid_point`（调用方 `ShrunkRange`/`find_valid_range` <- `pave_ff_make` / `pave_ff_pave_bound` / `pave_shrunk` / `pave_ff_se` / `pave_common`）· topo/edge_edge/mod.rs（`pub(crate)` 再导出）",
    state: "**已修（2026-10-15）**：原实现只按 `aStep = Resolution(theTol)*1.01` 逐步推进，并在注释里自注 `ponytail:` 有意省略 OCCT 的局部奇点逃逸臂。现按 `BRepLib_1.cxx:69-81` 补 `aD1Mag`：`GeomAbs_OffsetCurve` 先取 `OffsetCurve()->BasisCurve()` 再重取类型，仅 Bezier/BSpline 取 `aD1Mag = (1./Resolution(1.)*0.01)^2`（本移植 `curve.resolution(1.0)`），其余类型保持 `0.` 不进该臂；判型用单一忠实件 `edge_edge::curve_kind`＝`GeomAdaptor_Curve::GetType()`（`GeomAdaptor_Curve.cxx:252-311`），故在 `edge_edge/mod.rs` 加 `pub(crate) use find_solutions::{curve_kind, GeomCurveKind};`（纯可见性再导出，避免第二份 `GetType()` 实现）。再按 `cxx:108-141` 补加速循环：`aD1Mag > 0.` 时用 `D1` 同时取点与切矢，若点仍在球内且 `|D1|^2 < aD1Mag` 则 `aStepLocal *= 2.` 继续推进，仍在范围内就继续加倍，越界则回落到端点并以端点判 `isOut`（端点也在球内则整段无有效点，返回 `None`）；`aD1Mag == 0.` 保持原 `Value(u)` 步进臂。门禁：`--all-targets` 0 error（警告 373）、`--lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 且 v/f 与包围盒逐位一致。**范围偏差（已披露）**：`edge_edge/mod.rs` 不在 T005 预登记的 `allowed_files` 内，改动为纯可见性再导出（逻辑 0 处）。",
  },
  {
    item: "`GCPnts_AbscissaPoint::AdvCompute` 缺 `GCPnts_AbsComposite` 臂；两处 adv 入口的 `theEPSILON` 误传 `Resolution(Confusion)`",
    id: "A15",
    where: "geom/gcpnts.rs:compute_abs_composite_adv · abscissa_point_from_guess · abscissa_point_with_tolerance",
    evidence:
      "按 `GCPnts_AbscissaPoint.cxx:187-295` 补 `compute_abs_composite_adv`：逐步 `CPnts_AbscissaPoint::Length(..., theEPSILON)`、`anIndex == 0 && aDirection > 0` 前置块（`:202-229`，仅 adv 有）、`Init(C, U1, U2, theEPSILON)` + `AdvPerform`、循环内零长判据用 `Precision::PConfusion()`（`:232`）而前置块用 `theEPSILON`（`:206`）、尾段按 `aSign` 偏移并按 `FirstParameter`/`LastParameter` 钳制（`:273-293`）。同时订正 `.cxx:486`/`:540` 的 `theEPSILON` 语义：原先传 `c.resolution(Confusion)`，现传调用方 `tol`（`Approx_CurvlinFunc::GetUParameter` 传 `myTolLen`，`Approx_CurvlinFunc.cxx:560`）。`occt-geom` 143 测试通过、`step_obj_gates` 4/4。",
  },
  {
    item: "`occt-geom` 曲面族 Extrema 失败伪造（`fallback_point_surface` / `fallback_curve_surface` / `fallback_ss`，T001 的 B6）",
    id: "B6 / T-136",
    where: "geom/extrema_surf/numeric_extrema.rs（40×40 栅格 + 坐标下降、20×16×16 栅格）· geom/extrema_ss.rs（12⁴ 栅格）· geom/extrema.rs 三个薄包装 · geom/extrema_surf/analytic_solvers.rs · geom2d/extrema2d/analytic_solvers.rs · geom2d/extrema2d/curve_curve.rs · geom2d/extrema2d/ext_pc2d.rs · geom2d/extrema2d/gfunc_ext_pc2d.rs · topo/bean_face_exact.rs `ComputeUsingExtremum` · topo/brep_extrema.rs:209 `closest_point_on_face` · topo/gprop_analytic/edge_length.rs:492 `extrema_point_surface_exact` · topo/shape_fix_compose_shell/split_wire.rs",
    state: "**B6a 已修（2026-10-15）**：三处兜底整体删除；`point_surface_extrema`/`curve_surface_extrema`/`surface_surface_extrema` 改返回 `Option<ExtremaPair>`（`None` = `Extrema_ExtPS.cxx:307-371` / `Extrema_ExtCS.cxx:122,532` / `Extrema_ExtSS.cxx:226-260` 的 `!IsDone()` 或 `NbExt() == 0`），`numeric_extrema::no_solution_pair` 收拢原有的「无解」标记对。消费方 `bean_face_exact::ComputeUsingExtremum` 回到 `IntTools_BeanFaceIntersector.cxx:956` 的守卫（`Extrema_ExtCS` 未 done 即整块跳过、不置 range flag），`closest_point_on_face` 无解时跳过面内投影，`extrema_point_surface_exact` 非初等面臂改 `Err`。旧登记「有生产消费方故暂留、前置 T-67」不成立。门禁：`--all-targets` 0 error、`occt-geom --lib` 143/0、`--lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 且 v/f 逐位一致。**B6b 已修（2026-10-15）**：两侧采样式分类器全部换成 `GetType()`/`gp_*` 下转。`geom/extrema_surf/analytic_solvers.rs`：`classify_plane`＝`s.gp_pln()`、`classify_sphere`＝`s.gp_sphere()`（`Extrema_ExtSS.cxx:117-129` 的 `GetType()` 分派 + `Adaptor3d_Surface::Plane()/Sphere()`）、`reconstruct_line`＝`c.gp_line()`（与 `Extrema_ExtCS.cxx:129-140` 的 `C.Line()` / `myS->Sphere()` / `myS->Plane()` 同源），采样件 `is_planar`（8×8 法向并行判）、`sphere_center`（中纬度三点外心）、`circumcenter`、`surface_normal`、`is_line`（6 点切向）、`sample_bounds` 整段删除。`geom2d/extrema2d`：`point_curve_extrema2d_all` 线臂改 `c.gp_lin2d()` 并按 `Extrema_ExtPElC2d.cxx:63-79` 用 `Mydist = L.Direction().Dot(L.Location()→P)` 取**线上参数**（原实现取的是相对 `c.d0(first)` 的参数增量，裁剪线/带偏置起点的线上参数直接失真）、圆臂改 `c.gp_circ2d()`；`curve_curve_extrema2d_all` 的配对分派改 `gp_lin2d()`/`gp_circ2d()`（`Extrema_ExtCC2d.cxx:104-145` 的 `GetType(C1)`×`GetType(C2)` 与 `Extrema_Curve2dTool::Line/Circle`）；采样件 `is_line2d`/`circumcenter2d`/`classify_circle2d`/`line_of_curve2d` 删除。**可观测修复**：几何上直的 B 样条/近直线曲线原先被 6 点切向判据拖进线臂（OCCT 给 `GeomAbs_BSplineCurve`，落一般臂），裁剪线的点-线极值参数亦由增量纠正为线上参数。门禁：`--all-targets` 0 error、`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（117s）、`export_data_obj` 15/15 且 15 个模型 v/f 与包围盒逐位一致。**B6c 已修（2026-10-15）**：`extrema2d` 一般臂的种子栅格与其求解器整体删除 —— `analytic_solvers.rs::build_samples2d`（`span/0.1` clamp[24,256]）、`newton_point_curve_all2d`、`curve_curve::newton_curve_curve_all2d` 及 `fval2d`/`fprime2d`/`dist22d`/`solve_f_zero2d`/`refine_seed2d`/`pair2d_c`/`dedupe_sort2d` 全部下线；替身对应的 cxx 是 `Extrema_ExtPC2d` = `Extrema_GGExtPC`（`Extrema_ExtPC2d.hxx:31-38`；`Initialize` 取 `mytolu = Resolution(Confusion)`、`type = GetType(theC)`、`mysample = 17`，`Extrema_GGExtPC.hxx:119-128`；`Perform` 先把两端点距离入解，再按 `type` 分派：Line/Circle/Ellipse/Parabola/Hyperbola → `Extrema_ExtPElC2d`，Bezier → `NbPoles()*2` 个样本 + `IntervalPerform`，BSpline → 按度数与节点区间，`default:` → `NbIntervals(GeomAbs_C2)` / `DeflCurvIntervals` 逐区间 `IntervalPerform`）加 `IntervalPerform`（`hxx:590-620`，`myExtPC.Initialize(mysample, uinf, usup, tolu, tolf)` + `Perform`，`AddSol` 以 `PConfusion()` 归并）与 `Extrema_GFuncExtPC`（`Value`/`Values`，含 `D1`/`DF` 奇点处的 Taylor 与三点差分）；本轮落成单一忠实件 `extrema2d/ext_pc2d.rs`（`ExtremaExtPC2d`）与 `extrema2d/gfunc_ext_pc2d.rs`（`FuncExtPc2d`），`point_curve_extrema2d_all` / `point_curve_extrema2d`（改返回 `Option`，对齐 `Geom2dAPI_ProjectPointOnCurve::Perform` 的 `NbExt() > 0` 且取最小 `SquareDistance`）与 `curve_curve_extrema2d_all`（`Extrema_ExtCC2d`）接线到位；`general_extrema::proj_p_on_c` 改用 `ExtremaExtPC2d::new_range`（`Extrema_GGenExtCC.hxx:249-263`），其 **UNPORTED (audit A16)** 标记随之清除；`topo/shape_fix_compose_shell/split_wire.rs` 按 `ShapeFix_ComposeShell.cxx:1046-1087` 恢复 `apar` 默认 `firstPar`、仅有解时覆盖。四条既有 lib 断言按 cxx 次序订正（`Extrema_ExtPElC2d.cxx:113-150` 圆臂远点在前、`Extrema_ExtElC2d.cxx:334-364` 圆/圆 `P1[NoSol]×P2[kk]` 次序、平行直线 `IsParallel()` 且 0 极值），未新增测试。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（23 个面积比与基线逐位一致）、`export_data_obj` 15/15 且 v/f 与包围盒逐位一致。**B6e 已修（2026-10-15，包 P017 / T-147，commit `e3a4335e`）**：`extrema2d::locate_extcc2d`（`curve_curve.rs:255-262`）的自创 Newton 整体删除 —— 原 `refine_curve_curve2d_newton` 是夹取 32 次迭代 + `1e-12` 收敛判据，且 `f_and_jac2d` 的 Jacobian 四个带符号项相对 `Extrema_GFuncExtCC::Values` 全反；现按 OCCT 链落地：`Extrema_LocateExtCC2d.cxx:26-82` → `Extrema_LocECC2d` = `Extrema_GenLocateExtCC`（`Extrema_GenLocateExtCC.hxx:47-110`：逐轴量程规范化到 `[FirstParameter, LastParameter]`、`TolU`/`TolV` 取 `Adaptor::Resolution(Confusion)`、`MathFunctionSetRoot::perform_with_bounds`，随后 `F.NbExt() > 0` 与 `F.Value < 1e-10` 复核）→ `Extrema_CCLocFOfLocECC2d` = `Extrema_GFuncExtCC` 2D 实例。新忠实件 `extrema2d/gfunc_ext_cc2d.rs`（`GFuncExtCC2d`，`MathFunctionSetWithDerivatives` 两方程两变量：`SearchOfTolerance` `hxx:150-197`、`GetType` 分派 `hxx:263-291`、`Value` `hxx:359-550` 含 `D1` 奇点处的三点 / `DN` 重建、`Values` `hxx:568-827` 含奇点导数的扰动三点 Jacobian、`GetStateNumber` `hxx:832-870` 存 `(POnCurv2d, POnCurv2d)`、`Points`/`SubIntervalInitialize`）与 `extrema2d/gen_locate_ext_cc2d.rs`（`GenLocateExtCC2d` + `ExtremaLocateExtCC2d`）。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（248.5s）、`export_data_obj` 15/15 逐位一致。**B6g 已修（2026-10-15，包 P018 / T-148，commit `8cad6a01`）**：`extrema2d/curve_curve.rs` 里 `line_line_extrema2d`/`line_circle_extrema2d`/`circle_circle_extrema2d` 三个模块级副本（连带 `line2d_value`/`circle_param2d`/`dedupe_sort_cc2d` 三个只服务它们的辅助件）整体删除，测试与调用点改走唯一忠实件 `extrema2d/ext_el_c2d.rs` 的 `Extrema_ExtElC2d`（`new_line_line`/`new_line_circle`/`new_circle_circle`）。副本已与 cxx 分叉：平行直线返回了带点的「距离对」而 cxx 只置 `IsParallel()` 与 `SquareDistance(1)`（`Extrema_ExtElC2d.cxx:63-70`）、线/圆额外追加了 cxx 解析臂从不产出的交点（`cxx:105-167` 只给 ±Y 两对）、圆/圆多一路自创求交。`extrema2d/tests.rs` 五条断言按 cxx 次序与计数重写（平行线 `NbExt()==1` 且 `SquareDist 9`、正交线 `SquareDist 0` 且在原点、线 y=3 对单位圆 `sq 4/16`、线 y=0 过圆心给 ±Y 两对与 `sq 1`、圆/圆四对 `9,1,25,9`），未新增测试。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（150.7s）、`export_data_obj` 15/15 逐位一致。**B6e 的 3D 孪生已修（2026-10-15，包 P019 / T-149，commit `616ea0c1`）**：`occt-geom/src/extrema_cc/curve_curve.rs::locate_extcc`（`:152`）改调 `ExtremaLocateExtCC::new`，自创 `refine_curve_curve`（夹取 32 次迭代 + `1e-12` 判据）与 `f_and_jac`（Jacobian 四个带符号项相对 `Extrema_GFuncExtCC::Values` 全反）从 `poly_roots.rs` 删除；新增 `extrema_cc/p_on_curv.rs`（`POnCurv`＝`Extrema_POnCurv.hxx:25-60`）、`extrema_cc/gfunc_ext_cc.rs`（`GFuncExtCC`＝`Extrema_CCLocFOfLocECC`＝`Extrema_GFuncExtCC` 3D 实例，`MathFunctionSetWithDerivatives` 两方程两变量：`SearchOfTolerance`、`GetType` 分派、`Value` 的 `D1` 奇点三点 / `DN` 重建、`Values` 的奇点扰动三点 Jacobian、`GetStateNumber`、`Points`、`SubIntervalInitialize`）与 `extrema_cc/gen_locate_ext_cc.rs`（`GenLocateExtCC`＝`Extrema_LocECC`＝`Extrema_GenLocateExtCC`（`Extrema_GenLocateExtCC.hxx:49-110`：逐轴量程规范化、`TolU`/`TolV`＝`Curve::resolution(Confusion)`、`perform_with_bounds`、`NbExt() > 0` 与 `Value < 1e-10` 复核）+ `ExtremaLocateExtCC`＝`Extrema_LocateExtCC.cxx:24-74`）。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（181.8s）、`export_data_obj` 15/15 逐位一致。**B6f 已修（2026-10-15，包 P020 / T-150，commit `2452a51c`）**：`extrema2d::curve_curve_intersections2d` 的 33×33 采样栅格整体删除（原实现按 `bound()` 在两条曲线上取均匀参数格点、把距离小于 `tol` 的点对全部推入并加本地 `1e-6` 去重），改走 `Geom2dAPI_InterCurveCurve.cxx:54-106` 的忠实链：`Geom2dInt_GInter(AC1, AC2, Tol, Tol)`（已移植的 `Geom2dIntGInter::perform_curves`）→ 未 `IsDone()` 则 `NbPoints() == 0` → `i ∈ 1..=NbPoints()` 的 `Point(i).Value()`；只服务该栅格的 `analytic_solvers.rs::bound()` 删除。**可观测修复**：交点只要不落在 33 等分格点上就会被栅格静默丢弃或替换为邻近采样点，`Geom2dInt_GInter` 走曲线真实域（`ComputeDomain`，`Circle`/`Ellipse` 取整周期），计数与位置即求交器自身结果。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（205.1s）、`export_data_obj` 15/15 逐位一致。**B6d 已修（2026-10-15，包 P021 / T-151，commit `11b12e0e`）**：`occt-geom/src/extrema_ss.rs` 的单文件自创求解器换成 OCCT 忠实件。删除球/球、球/面、面/球三条解析臂（`sphere_sphere_extrema`/`plane_sphere_pairs`/`plane_sphere_extrema`/`sphere_plane_extrema`）与只服务它们的 `perpendicular`/`project_on_plane`，以及一般臂的自创 Newton（`ss_f`/`solve4`/`solve_surface_surface`）；新增 `extrema_ss/mod.rs`（`ExtremaExtSS`＝`Extrema_ExtSS`，`cxx:100-255` 按 `S1.GetType()` 分派：`Plane`×`Plane` 只走 `myExtElSS.Perform(S1.Plane(), myS2->Plane())`（`cxx:122-129`），其余类型对一律 `Extrema_GenExtSS(S1, *myS2, 20, 20, mytolS1, mytolS2)`，随后 `ElCLib::InPeriod` 归约周期参数 + 八条落窗判据过滤 `mySqDist`/`myPOnS1`/`myPOnS2`）、`extrema_ss/ext_el_ss.rs`（`ExtremaExtElSS`＝`Extrema_ExtElSS.cxx:40-54`）、`extrema_ss/gen_ext_ss.rs`（`ExtremaGenExtSS`＝`Extrema_GenExtSS.cxx:119-470`：S1/S2 各 20×20 抽样轴、全点对取 min/max 种子、`math_BFGS(4)` 改进 min、max 走 `math_FunctionSetRoot`）、`extrema_ss/func_ext_ss.rs`（`ExtremaFuncDistSS`/`ExtremaFuncExtSS`＝`Extrema_FuncExtSS.cxx:52-241`）、`extrema_ss/tests.rs`，`gp_Pln` 补 `signed_distance_pnt`/`signed_distance_pln`/`square_distance_pln`（`gp_Pln.hxx:345-379`、`:173-177`）。球/面这对在 OCCT 无解析臂（`Extrema_ExtElSS::Perform(gp_Pln, gp_Sphere)` 抛 `Standard_NotImplemented`，`cxx:62-69`），无界平面量程 ±2e100 由 `GeomAdaptor_Surface` 原样传入 `GenExtSS`，故旧断言的 4.0 是被删解析臂编出来的数；`occt_probe --ssx` 重编译实跑给出真值（`SSX sph-pln nbExt=1 sqDist=7.992002e+200`、`SSX gen-sph-pln nbExt=2 [1] 2.2138509695290867e+198 [2] 7.992002e+200`、`SSX api-sph-pln dist=2.827012911183817e+100`、`SSX sph-sph nbExt=2 sqDist=4.0000000000003126/36`），三条既有断言随之订正。门禁：`occt-core --lib` 283/0、`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（114.1s）、`export_data_obj` 15/15 逐位一致。**同批登记（未改）**：**B6h** —— `extrema_ss/mod.rs::surface_surface_extrema` 取 `mySqDist` 首项，而 OCCT 对应件 `GeomAPI_ExtremaSurfaceSurface` 取最小 `SquareDistance` 的下标 `myIndex`（`GeomAPI_ExtremaSurfaceSurface.cxx:71-91`，`NearestPoints`/`LowerDistance` 都用它），`NbExt() > 1` 时两者不同；同处 `point_on_s1`/`point_on_s2` 多暴露了 OCCT 的私有序列（`Extrema_ExtSS.hxx:99-102` 私有，对外只有 `Points`），两者均无生产调用方。B6 家族（B6a/B6b/B6c/B6d/B6e 2D+3D/B6f/B6g）至此全部收口。",
  },
];

// 语义偏移残余（对照 specs/_audit/_index.md A 编号；位置/行号为 2026-10-13 工作树实测）
const SEMANTIC = [
  {
    item: "Delaunay 插点密度：`BRepMesh_DelaunayDeflectionControlMeshAlgo` 的插点与最小边判据",
    id: "A19 余项",
    where: "topo/meshing/deflection_control/ · topo/meshing/node_insertion.rs",
    state: "**已移植（2026-10-15 行内复核，原「未移植」条目标题本身有误）**：OCCT 8.0 的该类没有 `checkAndAdd`/`addInternalNodes` 这两个符号（全树 grep 无命中），行内把 `optimizeMesh`/`splitTriangleGeometry`/`splitLinks`/`usePoint`/`checkDeflectionOfPointAndUpdateCache`/`rejectByMinSize`/`rejectSplitLinksForMinSize`/`checkLinkEndsForAngularDeviation`（`BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx:75-454`）与 `insertNodes`/`registerSurfaceNodes`/`getCellsCount`（`BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx:74-155`）全部落地，并接线（`node_insertion.rs:519-679` 的 `post_process` + `insert_parametric`）；`A_ITERATIONS_NB = 11` 与 `BRepMesh_GeomTool::CellsCount` 亦已对齐。原 TDB +1,744 面证据经逐面配对复核归属缝边 pcurve 注册（2.1 T-124）而非本项；2026-10-16 复测残余为双向小量：T0M −469、top −30、ATU01038 −32、a3n00 −24、acs10 −20、bottom −7、motoc 0、TDB +17（面数；读数见 `.target-gate/t157/gates-raw.txt`），面积比 20 个 1.0000 + `Shape-2` 0.9999 + `T0M` 1.0017 + `TDB` 1.0015 —— 无 4~9% 的系统性偏密（2026-10-13 记的 T0M +58 / TDB −51 与「面积比全 1.0000」是当轮读数，已被 10-15～10-16 各包的共享路径忠实化改写）",
  },
  {
    item: "bop_curved 网格级布尔（曲面 operand）",
    id: "A5 / T-79",
    where: "topo/bop_curved/region_trim.rs:51 · face_meshing.rs:418",
    state: "已闭环（2026-10-14）：网格引擎已从 `boolean_dispatch` 的曲面分支摘除（改调 `bop_builder::boolean_via_bopalgo`），忠实路径的曲面布尔对齐 —— `zz_bop_curved_probe`（box[-1,1]^3 与 `cylinder(r=0.4,h=2)`）`FUSE` faces=8/solids=1/closed/vol_gauss=8.50265、`CUT` faces=8/vol_gauss=7.49735、`COMMON` faces=3/vol_gauss=0.50265，三项全中 OCCT ground truth（此前 `FUSE` faces=11/vol=10.720）；门禁 `--lib` 1276/0 + 12 个集成套件全绿（1345 passed / 0 failed）",
  },
  {
    item: "`brepfeat` 自创件：网格夹具冒充解析件 + `solid_volume` 线框扇形和（OCCT 无此规则）",
    id: "A12 余项 / T-125",
    where: "topo/brepfeat/features.rs translated_cylinder（`BRepPrimAPI_MakeCylinder`）· :revolve_profile_about 车削网格（余项）· solid_volume 已改 `brep_gprop_full::volume_properties`",
    state: "**2026-10-15 收口第 (2) 项的一半**：(1) 体积积分已换忠实件（2026-10-13，`solid_volume` 调 `BRepGProp::VolumeProperties`，自创扇形和/Newell 法向与 `polygon_fan_volume`/`polygon_normal`/`wire_polygon_3d` 已删）。(2) **解析圆柱工具已落地** —— `boss`/`boss_thru_all` 的工具由 `translated_mesh_cylinder`+`mesh_cylinder(..,24)` 改为 `BRepPrimAPI_MakeCylinder`（`translated_cylinder`），死代码 `mesh_cylinder`/`mesh_volume_signed`/`resolution_for`/`all_faces_planar` 一并删除；实测 `box ∪ boss` = 4.157080（= 4 + π·0.0625·0.8，精确；原 24 边形 4.155291）、`box ∪ 贯通 boss` = 2.565487（裁剪值，不再是 2.706858 的未裁剪并集）—— 原「被 `BOPAlgo_BOP` 阻塞」的结论随 T-79 闭环而失效。(3) **解析旋转体工具仍不可用**，`revolve_profile_about` 的车削网格夹具保留并就地标 UNPORTED：直接换解析体后 `box ∪ 解析 neck` 丢工具（4.000000 应 4.489835）、`faceted cylinder − 解析 groove 环` 反增体（9.764392 > 9.424778）；再按 `BRepPrim_OneAxis::BottomFace` 补面向外法向反转后变为 neck 零体积、groove 6.911504（应 7.816283）⇒ 缺口不在 `BOPAlgo_BOP`，在 `sweep_revolve` 未按 `BRepSweep_Revol`/`BRepPrim_OneAxis` 的 wire/seam/pcurve 约定建面（即 `BRepPrimCylinder` 为让 `BOPAlgo_WireSplitter::Tolerance2D` 闭合缝环所必须满足的那组约定）。门禁：`brepfeat --lib` 17/17、`occt-topo --lib` 1251/0",
  },
  {
    item: "`heal_tolerance` 自创修形链 + `boolean_with_check` 容差重算",
    id: "A11 / 发现 16",
    where: "topo/bop_builder_splitapi.rs:301 · topo/bop_builder_dispatch.rs:504",
    state: "已修：`remove_small_edges` 改走 `ShapeFix_Wire::FixSmall`+`ShapeAnalysis_Wire::CheckSmall`（无弧长阈值），`remove_degenerate_faces` 换 `ShapeFix_Face::FixSmallAreaWire`（`CheckSmallArea`），`weld_coincident_vertices` 端点替换改 `ShapeBuild_Edge::CopyReplaceVertices`，`boolean_with_check` 的 tol×10/×100 重算删除（OCCT 单一 `myFuzzyValue`）",
  },
  {
    item: "自创自适应 Simpson 余项：`occt-math/integrate_adaptive.rs`（无 live 调用者）· `occt-core/elib/measure.rs::simpson`（唯一 live 消费者 `topo/inttools_range.rs:174`）",
    id: "A15 余项",
    where: "core/elib/measure.rs:58 · math/integrate_adaptive.rs:5",
    state: "已闭环（2026-10-13）：`GCPnts_AbscissaPoint::Length` 类型分派已补（`geom/gcpnts.rs:length`，`GCPnts_AbscissaPoint.cxx:371-423` 的 `GCPnts_LengthParametrized`/`GCPnts_Parametrized`/`GCPnts_AbsComposite` 三臂），`inttools_range.rs:174` 改调 `gcpnts::length(curve, s1, s2, Some(a_ptol_e))`（`IntTools_ShrunkRange.cxx:170`）；两个自造积分件已删（`occt-math/integrate_adaptive.rs`、`occt-core/elib/measure.rs`，零残余消费者）",
  },
  {
    item: "`BRepCheck_Edge` 的 SameParameter/SameRange 检查",
    id: "—",
    where: "topo/brep_check.rs",
    state: "已移植 `Minimum`（2026-10-13，`BRepCheck_Edge.cxx:74-258`）：`!SameRange && SameParameter` → `InvalidSameParameterFlag`、无 3D 表示 → `No3DCurve`、空曲线表示 + `Degenerated` → `InvalidDegeneratedFlag`、`Last <= First` 与基曲线周期/范围比对 → `InvalidRange`（`Geom_TrimmedCurve`/`Geom2d_TrimmedCurve` 按 BasisCurve 取 f/l/IsPeriodic/Period）；`Status` 已按 `BRepCheck_Status.hxx:20-59` 重排补齐。UNPORTED：`InContext`（依赖未移植的 `BRepLib_ValidateEdge`，`BRepCheck_Edge.cxx:263-556`）、`Tolerance`（`:598-707`）、`CheckPolygonOnTriangulation`（`:711-839`）、表示的 `TopLoc_Location`（`:163-171`）",
  },
  {
    item: "`surface_closest_params` 采样替代件",
    id: "A1",
    where: "topo/brep_surface.rs:261",
    state: "已删（2026-10-13）：256 采样 + 6 次对分的网格搜索整段移除，活调用点 0",
  },
  {
    item: "`brep_surface::is_planar` / `face_is_planar` 的 8×8 法向采样（发现 14 余项）",
    id: "发现 14 余项",
    where: "topo/brep_surface.rs（原 :61,:188 采样体，已删）",
    state: "已闭环（2026-10-15，批次 3/T-129 + 批次 4/T-130）：余 19 处调用点全部换成 `is_plane_type` / `face_is_plane_type`（`GetType() == GeomAbs_Plane` 等价）——`brep_surface.rs:153`（`face_plane`）、`pcurve.rs:249`、`fillet_curved/surface_info.rs:139,220`、`fillet_curved/torus_blend.rs:603`、`gprop_analytic/analytic_props.rs:378,587`、`gprop_analytic/mass_properties.rs:17,61,325,443`、`brep_offset/shell_offset.rs:122`、`brep_faces.rs:203` 及测试夹具（`brep_builder_full/mod.rs:247`、`brep_sketch.rs:234`、`brep_builder_api.rs:372,405-406`、`fclass2d/tests.rs:160`、`brep_offset/tests.rs:132`、`fillet_curved/tests/helpers.rs:693`、`brep_surface.rs:288,304`）；随后 `is_planar`/`face_is_planar` 采样体整段删除，全仓已无该符号。同批把 `edge_face_kind.rs` 的曲线/曲面族（`curve_kind`/`circle_geometry`/`surface_kind`/`is_eq_distance`/`is_coplanar`/`is_radius`）换成 `GetType()`/`gp_*` 同等 downcast。批次 4（T-130）再把 `Resolution` 采样族换成 `Adaptor::Resolution` 解析闭式：`edge_face_kind::curve_resolution`、`inttools_range::curve_resolution`（3 处）、`algo_tools_range` 的 `curve_is_line`+`curve_is_freeform`+`curve_resolution`、`bean_face::resolution`、`int_tools_vertex_line::adaptor_resolution`、`wire_splitter_block::angle_2d`，并在 `occt-geom2d` 新增 `Curve2d::resolution`（`Geom2dAdaptor_Curve::Resolution` 等价）与两个装饰器（`ReparamCurve2d`/`PlaneKeepParam2d`）的 `resolution`。门禁：`occt-geom2d` 72/0、`--lib` 1251/0、`--tests` 1317/0、`step_obj_gates` 4/4（23 模型面积比 1.0000）、`export_data_obj` 15/15 且 v/f 逐位一致。**新登记的同族残留**：`inttools_range::find_nearest_valid_point` 缺 `BRepLib_1.cxx:69-83` 的 `aD1Mag` 局部奇点加速分支（原注释已标 `ponytail:`）。**该残留已于 T-133（2026-10-15）收口**（见 B5 / T-133 行）：阈值与 `D1` 加速倍增循环整体落地，判型复用 `edge_edge::curve_kind`＝`GeomAdaptor_Curve::GetType()`，`ponytail:` 注记删除。",
  },
  {
    item: "`BOPAlgo_WireSplitter` 的 `Angle2D` / `Tolerance2D` / `RefineAngle2D` 三臂（`_audit/bop-boolean.md` 发现 19）",
    id: "发现 19 / T-131",
    where: "topo/wire_splitter_block.rs（原 `uv_tolerance_2d:580` · `curve2d_radius:543` · `refine_angle_2d:771`）",
    state: "**已修（2026-10-15）**：(1) `uv_tolerance_2d` 原先手抄 `GeomAdaptor_Surface::U/VResolution`，回转(V)/拉伸(U) 两臂退回有限差分 `fd_uv_tolerance`（`t3/|dS/du|`），且缺 `RectangularTrimmedSurface` 下基（`cxx:423-430`）、缺 `OffsetSurface` 臂（`cxx:1882-1885`/`:1944-1947`）、缺锥面无界守卫（`cxx:1855-1858`）、锥半径用采样代 `S->VIso(V)->Radius()`（`cxx:1861-1866`）——现改转发唯一忠实件 `approx_same_parameter::{u,v}_resolution`，手抄 switch 与差分件删除。(2) `angle_2d` 的内联 `|D1 x D2|/|D1|^3`（守卫 `n2 <= PConfusion` 量纲不符）改 `GeomLProp_CLProps2d` + `is_tangent_defined()`/`curvature()`（`BOPAlgo_WireSplitter_1.cxx:797-807`）。(3) `refine_angle_2d` 的 `aTV` 由 3D 边参数改 pcurve 参数 `BRep_Tool::Parameter(aV,aE,aFace)`（`cxx:1060`）；UNPORTED 注记清除。门禁：`--lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 且 v/f 逐位一致。",
  },
  {
    item: "`pcurve::pc_curve_kind` 几何不变量启发式分类器（`_index` 未登记，T001 的 B2）",
    id: "B2 / T-132",
    where: "topo/pcurve.rs（原 `CurveKind` / `pc_curve_kind`，已删）· topo/geom_bnd_lib_curve2d.rs · topo/geom_bnd_lib_offset2d.rs · topo/geom_bnd_lib_bezier2d.rs · topo/int_face_face_bounds.rs · geom2d/adaptor_curve.rs（新增）· geom2d/extrema2d/curve2d_tool.rs",
    state: "**已修（2026-10-15）**：原先 4 类分类（`is_periodic() && |period-2*PI|<1e-9` ⇒ Circle、无界 ⇒ Line、否则 BSpline）被 `Geom2dAdaptor_Curve::GetType()` 同等件取代 —— 新增 `occt_geom2d::adaptor_curve::{basis2d,get_type2d}`（`Geom2dAdaptor_Curve.cxx:96-120` + `load` 类型表 `:279-346`，含 `Geom2d_TrimmedCurve` 下基 `:285-288`），`extrema2d/curve2d_tool` 的 `adaptor_basis2d`/`adaptor_type2d` 改为转发该唯一实现。(1) `geom_bnd_lib_curve2d::curve2d_abs_type` 改读 `get_type2d`，`box_curve2d` 按 `GeomBndLib_Curve2d.cxx:97-171` 的 downcast 链补全九臂：Line/Circle 优先取 `gp_lin2d`/`gp_circ2d`（采样重建件 `lin2d_from_curve`/`circ2d_from_curve` 保留为 `PlaneKeepParam2d` 这类未转发 `gp_*` 的 wrapper 替身）、Ellipse/Hyperbola/Parabola 取 `gp_elips2d`/`gp_hypr2d`/`gp_parab2d` 走解析盒（原先一律落采样）、Bezier 走 `GeomBndLib_BezierCurve2d::Box` 极点包络（`cxx:30-52`，子区间无 2D `Segment` 时退回采样）、BSpline 走 `bspline_copy2d`、Offset 走 `offset_basis`+`offset_value`。(2) `geom_bnd_lib_offset2d` 的 `try_offset_line2d`/`try_offset_circle2d` 由 `pc_curve_kind` 改用 `GeomBndLib_OffsetCurve2d.cxx:36,63` 的精确 `occ::down_cast<Geom2d_Line>/<Geom2d_Circle>` 等价（`gp_lin2d()`/`gp_circ2d()`），基圆圆心半径不再靠 `d0(0)`/`d0(pi)` 采样反推；`box_offset_as_curve` 由恒等 `box_other` 改为按 `cxx:123-143` 先试线/圆快路径。(3) `int_face_face_bounds.rs` 的 `Geom2d_Line` 判定改精确 `gp_lin2d()` 并直接用 `Position().Direction()`（`IntTools_FaceFace.cxx:2149,2161-2162`），原 `d1(0)` 切向改判方向。(4) `pcurve_full/tests.rs`、`pcurve.rs`、`boptools_2d.rs` 中引用 `pc_curve_kind` 的既有断言改断言 `curve2d_abs_type`。门禁：`occt-geom2d` `cargo check` 0 error；`--lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15（15 个模型 v/f 与包围盒与基线逐位一致）。**同批发现、未改 → 已由 P026 / T-156 收口（2026-10-16）**：`ReparamCurve2d`（`GeomLib::SameRange` 的线性替身）不再转发内层 kind 给 `bspline_copy2d`（该 override 已删除 —— wrapper 只在 `Geom2dConvert::CurveToBSplineCurve` 拒绝时存在，`Geom2dConvert.cxx:442-444`，故 trait default `None` 即正确），且 `pcurve_full/surface_projector/reparam.rs` 与 `shhealing/wire_fix/pcurve_range.rs` 的两份 `SameRange` 转换尾副本合并为唯一 `trimmed_reparam_bspline`（`GeomLib.cxx:912-921` / `:962-968`）与 `same_range_window`（`cxx:934-959`）；`wire_fix` 副本原先在 cxx 不判型的不等跨臂（`cxx:924-969`）按曲线类型分派（`line_trim_to_reparam_bspline` / `unequal_trim_window` → `same_range_unequal_line` / `same_range_unequal_other`），该分派已删除。",
  },
  {
    item: "`face_sample_point` 采样式 BOP 面代表点（`_audit/bop-boolean.md` 发现 5 的残留点选择，T001 的 B4）",
    id: "B4 / T-135",
    where: "topo/builder_solid.rs · topo/bop_build_solids.rs · topo/bop_build_common/{build_ops,draft_solid}.rs · topo/bopalgo_tools_class.rs · topo/algo_tools_face.rs",
    state: "**已修（2026-10-15）**：4 处重复的 `face_sample_point`（UV 中心 / 首边中点）与本地 `is_split_to_reverse` 替身全部改走唯一忠实件 `algo_tools_face::{compute_state_face_point,compute_state_face}` = `BOPTools_AlgoTools::ComputeState(face,solid,tol,bounds,ctx)`（`BOPTools_AlgoTools.cxx:666-715`：先取不属于 solid 的非退化边的中点（边重载 `:734-788`，含无界参数 `aT2-10`/`aT1+10`/`0` 臂与空曲线取首顶点），否则 `PointInFace`（`:692`），其失败再对非退化边逐个 `PointNearEdge`（`:694-707`））与 `algo_tools_face::is_split_to_reverse` = `IsSplitToReverse(Sp,Sr,ctx,err)`（`:1255-1291` / 面 `:1316-1427` / 边 `:1432-1500`）。(1) `builder_solid::{solid_inside_solid,shell_inside_solid}` 合并为 `shape_inside_solid` = `BOPAlgo_BuilderSolid.cxx:835-861` 的静态 `IsInside`：取 S1 首面用固定 `Precision::Confusion()` 分类（`cxx:854-858`，不再用 `self.fuzzy`），无面走无限点（`cxx:846-849`）；删 `face_sample_point`/`classify_point`。(2) `bop_build_common::is_split_to_reverse` 由本地 `face_split_to_reverse`/`edge_split_to_reverse`（单中点 + 32 点 `closest_point_on_edge`）与有限差分法向 `surface_closest_normal`（`1e-6` 步长）改转发忠实件，四个采样件删除。(3) `bop_build_solids::is_split_to_reverse` 删除、调用点改走忠实件；`classify_solid_state` 由 16×16 采样改 `AlgoTools::compute_state`（= `BRepClass3d_SolidClassifier`，`cxx:790-822`）。(4) `draft_solid::compute_state_by_one_point` 改为转发 `bopalgo_tools_class::compute_state_by_one_point`（单一实现），采样件 `point_solid_state` 删除。(5) `algo_tools_face::is_split_to_reverse_face` 的同面快路径由 `shape_key` 相等改为 `Arc::ptr_eq`（`Geom_Surface` 句柄相等，`cxx:1330-1333`，对被分割面与源面成立）。`bop_build_solids::face_sample_point` 仅剩已登记的 port-only 启发式 `orient_face_outward` 一个用户，函数头已加注。门禁：`--lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 且 v/f 逐位一致。**同批发现、未改**：OCCT `ComputeStateByOnePoint` 的 `default` 臂递归首子形状后丢弃返回值（`cxx:648-655`）⇒ 容器恒 `TopAbs_UNKNOWN`，本仓两处实现对子形状循环取首个非 UNKNOWN 状态，属已登记偏差（牵动 `FillInternalShapes` 容器分支）。",
  },
  {
    item: "`Analyzer::is_valid` 恒 true",
    id: "A14",
    where: "topo/validate.rs:150",
    state: "已声明非校验、不得作门禁",
  },
  {
    item: "`Geom_OsculatingSurface` 失败臂提前 `return`（T001 的 B8 osculating 条）",
    id: "B8 / T-139",
    where: "geom/osculating_surface.rs（`init`，旧 `:141,202,238,299`）· 可达 `offset_surface.rs:56` → `brep_offset`",
    state: "**已修（2026-10-15）**：`init` 原先在 `build_osculating` 失败时四处 `clear_flags(); return;`，把余下的 `for` 结点迭代连同后面整段 V 向扫描一起丢掉，并把五个臂的成功与否挤进一个外层 `oscul_ok`。OCCT `Geom_OsculatingSurface.cxx:179-378` 每臂各持局部 `OsculSurf`（下轮迭代重新赋值），失败只 `clearOsculFlags()`（`:206-213`、`:273-276`、`:308-312`、`:373-376`）后继续；只有 `if (myAlong[1] && OsculSurf)` 这类守卫因标志被清而跳过。现新增 `build_chain`＝`cxx:187-206` 的 `while (IsQPunc)` 链（返回 `Option<(GeomBSplineSurface, k)>`，`None` ＝ 空 `OsculSurf`，`k` 为成功次数），两个 `for` 循环改为 `if self.along[0] { A; if self.along[1] && ok { B } } else if let Some(..) { C } else { clear_flags() }` 的逐臂判定；四处 `return` 与共享 `oscul_ok` 删除。可观测修复：一次失败后余下的结点索引不再被整段跳过。门禁：`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 逐位一致。**同批订正**：B1（`wire_splitter_block.rs` 三臂）实已在 T-131（`67d7f01f`）收口，清单原标 LIVE 属陈旧；`extrema_cc/mod.rs:9`、`occt-math/globoptmin.rs:590` 的 UNPORTED 注释亦为陈旧（代码已用 `math_GlobOptMin`／已用 `newton_minimum.rs`）。**旁支已收口**：本项点名的 `surface_fit.rs::fit_plane` 已在包 P024（T-154，2026-10-16）处理：真实对应件不是 `GeomLib_IsPlanarSurface` 而是 `BRepLib_FindSurface`（`BRepLib_FindSurface.cxx:52-628`）经 `BRepLib_MakeFace(W, OnlyPlane)`（`BRepLib_MakeFace.cxx:188-256`），两者已忠实移植（`crates/occt-topo/src/brep_lib_find_surface.rs`、`brep_lib_make_face.rs::make_face_from_wire`）并接线 `brep_offset/curve_face_offset.rs::offset_wire_3d`；`GeomLib::IsPlanarSurface`/`Controle`/`GeomLib::Inertia` 在仓内无对应件也无可达分支。`occt-geom/src/surface_fit.rs` 的 `fit_plane`/`fit_surface_kind`/`fit_surface` 为 port-only 且已无生产调用方，是否删除留待后续裁决。",
  },
  {
    item: "STEP 写侧周期 B-spline 未收周期（T001 的 B7-writer 第 1 条）",
    id: "B7-writer / T-140",
    where: "topo/step/write_context.rs（`emit_curve_entity` 的 BoundedCurve 臂 + 裁剪 B-spline 臂）",
    state: "**已修（2026-10-15）**：`GeomToStep_MakeBoundedCurve.cxx:44-51` 写周期 B-spline 时先走非周期副本（`BS = B->Copy(); BS->SetNotPeriodic()`），本仓直接按曲线自身数据建实体，周期结点向量原样写出。现新增 `bounded_bspline_of_curve`（＝`cxx:42-62` 整臂：曲线自身 poles/knots/weights/degree 后接收周期步，键为 `C->IsPeriodic()`，`cxx:45`）与 `open_periodic_bspline`（＝`Copy()` + `SetNotPeriodic()` 对；副本必须显式恢复 `periodic`，因为移植的 `GeomBSplineCurve::new`/`rational` 构造器一律建非周期形，否则 `set_not_periodic()` 立即返回）。同文件裁剪 B-spline 臂（`MakeCurve.cxx:71-76`）在 `segment(...)` 前补 `bs.periodic = basis.is_periodic()`，使 `Geom_BSplineCurve::Segment` 进入周期分支（`Geom_BSplineCurve.cxx:606-618`：`SetOrigin` + `SetNotPeriodic`）。该处整段 UNPORTED 注释删除。门禁：`occt-topo` `--all-targets` 0 error、`--lib` 1251/0、`phase4_integration` 9/9、`step_obj_gates` 4/4、`export_data_obj` 15/15（v/f+bbox 与 T002 冻结表逐位一致；`screw` 601/792、`Shape-1` 3335/4328、`Shape-2` 3081/4744 即冻结值）。**验证边界**：四道门禁都不直接覆盖 STEP 写侧（写侧只由 `step/tests.rs` 与 `phase4_integration` 的回读用例覆盖），故本轮靠 cxx 控制流对齐 + 既有回读基线不变作证，未新增测试。**同族未修（登记 B7-writer-3）**：`step/format.rs:472` 的曲线 `closed` 恒 `.F.`（应 `BS->IsClosed()`，`GeomToStep_MakeBSplineCurveWithKnots_gen.pxx:42-50`）、`:526` 与曲线侧的 `KnotSpec` 恒 `.UNSPECIFIED.`（`gen.pxx:71-87` / `MakeBSplineSurfaceWithKnots.cxx:127-150`）；依赖 `Geom_BSplineCurve::IsClosed`、`BSplCLib::KnotDistribution`（枚举 `BSplKnotDistribution` 已在 `occt-core/src/kernel/geomabs.rs:37`）。**该同族余项已收口（2026-10-15，T-141 / 包 P011）**：新增 `occt-core/src/bspl/knot_analysis.rs`（`BSplCLib::KnotForm`/`MultForm`/`KnotAnalysis`，`BSplCLib.cxx:602-753`）与 `GeomBSplineCurve::{is_closed,knot_distribution}`（`Geom_BSplineCurve_1.cxx:146-149,363-366`）、`GeomBSplineSurface::{u,v}_knot_distribution`（`Geom_BSplineSurface_1.cxx:1851-1862`），`step/format.rs` 两个写器改由几何推导（曲面仅 U/V 分布相同时具名，`MakeBSplineSurfaceWithKnots.cxx:127-150`）；偏离有 OCCT 自身输出可对照（`data/occ/ATU01038.step` 245 `PIECEWISE_BEZIER_KNOTS` + 106 `QUASI_UNIFORM_KNOTS`、`bottom.step` 934）。门禁：`occt-core --lib` 285/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`phase4_integration` 9/9、`step_obj_gates` 4/4、`export_data_obj` 15/15 逐位一致；`step/tests.rs` 两条既有记录布局断言按 cxx 推导值重定（未新增测试）。",
  },
  {
    item: "`intana` / `occt-math` 的第二份自创多项式解算器（T001 的 B8-intana 第 1 半）",
    id: "B8-intana / T-142",
    where: "geom/intana/analytic_intersections.rs（本地二次/三次/Ferrari 四次 + `poly4`/`polish4`）· occt-math/polynomial.rs（闭式 + 自创 Durand-Kerner）",
    state: "**已修（2026-10-15，包 P012，commit `642a0c88`）**：`math_DirectPolynomialRoots` 其实早已移植为 `occt_core::math_direct_poly_roots::DirectPolynomialRoots`（`math_DirectPolynomialRoots.cxx:66-506`、`:558-625`、`:627-708`、`:712-764`、`:768-786`，已被 `intana_intconicquad`、`intana2d/outils`、`math_trig_roots` 使用），但 `intana/analytic_intersections.rs` 与 `occt-math/polynomial.rs` 各带一份第二自创解算器（阈值 `1e-14`/`1e-9`/`1e-12`/`1e-6`/`1e-8`/`1e-300` + 本地牛顿抛光），且头注声称 \"occt-geom 不依赖 occt-math\"（实际两者都依赖）。现已为 `occt-core/src/math_direct_poly_roots.rs` 增 `quadratic_roots`/`cubic_roots`/`quartic_roots` 三个 `Value(1..NbSolutions)` 适配器（即 OCCT 调用点内联写的循环，自身无数值判定），`intana` 侧整段本地解算器与 `poly4`/`polish4` 删除改为 re-export，`occt-math` 侧闭式改委托、自创 Durand-Kerner `polynomial_roots`（OCCT 无对应）与其复数辅助件、单测删除。**唯一边界**：本半只覆盖 `intana` 与 `occt-math` 两份拷贝；同一解算器仍在 `extrema_pc/poly_roots.rs`、`extrema_cc/poly_roots.rs`、`extrema2d/analytic_solvers.rs`、`occt-core/numeric.rs` 各有一份，`line_torus.rs` 的绝对 `1e-14` 与锥/锥 `EstimDist`/`aTolAng` 亦未动（已写入下一步）。门禁：四 crate `--all-targets` 0 error、`occt-core --lib` 285/0、`occt-math --lib` 210/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 逐位一致；未改任何断言。**该边界已收口（2026-10-15，T-143 / 包 P013，commit `8db13d5c`）**：`extrema_pc`/`extrema_cc`/`extrema2d` 三份改 `pub(super) use occt_core::math_direct_poly_roots::{cubic_roots, quadratic_roots, quartic_roots}`，`occt-core/numeric.rs` 的第四份自创件与其两条测试删除（`occt-core --lib` 285→283）；`extrema_pc/tests.rs:146-152` 的升序断言改为按集合比较（未新增测试）。",
  },
  {
    item: "锥/锥角度容差由调用方常数决定（T001 的 B8-intana 第 2 半）",
    id: "B8-intana / T-144",
    where: "geom/intana/line_torus.rs（`quadric_quadric_cone_cone` 多收 `tol_ang`）· 调用方 intpatch_impimp_pairs.rs / int_face_face_analytic.rs",
    state: "**已修（2026-10-15，包 P014，commit `dd31babc`）**：`IntAna_QuadQuadGeo::Perform(gp_Cone, gp_Cone, Tol)`（`IntAna_QuadQuadGeo.cxx:1433-1476`）只有一个入参容差并自行推出角度容差，本仓却让每个调用方各带常数（BOP 对照表 `TOL_ANG = 1e-8`、别处 `Precision::Angular()`）决定同轴判定。现入参 `tol_ang` 删除：`aTolAng` 起点 `myEPSILON_ANGLE_CONE`（`= Precision::Angular()`，`InitTolerances` `cxx:351-359`），仅当 `|tg1 - tg2| < Tol`、`A1A2.Parallel()`、`A1A2.Distance() > 100 * Tol` 时放宽为 `Tol / EstimDist(Con1,Con2)` 并钳到 `[Precision::Angular(), Tol]`（`cxx:1440`）；`EstimDist`（`cxx:260-326`）全量移植（过两锥顶与锥 1 方向的平面 `gce_MakePln(P1,P2,P3)` `gce_MakePln.cxx:76-92`、两轴投影 `ProjLib_Plane::Project` `ProjLib_Plane.cxx:102-107`、±`SemiAngle`、四对直线用 `IntAna2d_AnaIntersection::Perform(gp_Lin2d,gp_Lin2d)` 取 `max(min|ParamOnFirst|, min|ParamOnSecond|)`），新增唯一忠实件 `occt-geom/src/projlib.rs::project_pln_lin2d`；`cxx:1555` 嵌套平面/锥 reduce 订正为 `(Tol, Tol)`。**同批订正审计**：`line_torus.rs:336,340` 的两个 `1e-14` 是 `AxeOperator` 构造默认值（`cxx:69-72`）与 `A1A2.Same()`（`cxx:83`），忠实而非自创 —— 原先记作「绝对阈值替代自适应 `aTolAng`」的是 370 行，已在本包收口。门禁：`occt-geom`/`occt-topo` `--all-targets` 0 error、`occt-geom --lib` 143/0、`occt-geom2d --lib` 72/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（190.0s，23 面积比一致）、`export_data_obj` 15/15 逐位一致；未新增测试。**同文件未修（登记）**：`IntAna_QuadQuadGeo.cxx` 圆柱/锥（`cxx:1324-1342`）与圆柱/圆柱（`cxx:1050-1293`）分类臂**经 P015 逐臂复核实为已移植**（分别由 `quadric_quadric_cylinder_cone` / `quadric_quadric_cylinder_cylinder` 覆盖，audit 原记「未移植」属陈旧）；`IntAna_IntLinTorus.cxx:116-120` 的 `done=false` 与 `line_torus.rs:669` 的输出重排（OCCT `cxx:104-108` 按根序号写入，不重排）两项未改。",
  },
  {
    item: "锥/锥平行臂的嵌套结果原样返回（T001 的 B8-intana 尾项）",
    id: "B8-intana / T-145",
    where: "geom/intana/line_torus.rs（`quadric_quadric_cone_cone` 分支 2）· 消费方 topo/int_face_face_analytic.rs::cone_cone",
    state: "**已修（2026-10-15，包 P015，commit 见 `_tasks.md` T-145）**：分支 2 把嵌套 `quadric_quadric_plane_cone` 的结果原样返回，而 `IntAna_QuadQuadGeo.cxx:1557-1603` 对它做 `switch` 重映射 —— `Ellipse`/`Circle`/`Hyperbola` 保留，`IntAna_Line`（cxx 自身 nbint 1 或 2 皆然）报成同一条线两次（`pt1 = pt2 = H.Location()`、`dir1 = dir2 = H.Position().Direction()`、`nbint = 2`）⇒ 本仓 `TwoLines(l, l)`，其余（`Parabola`/`Point`/`Empty`）⇒ `IntAna_NoGeometricSolution` ⇒ `None`。原先 `Parabola`/`Point` 会被当成精确交集曲线用掉，线数也与 cxx 不符。同轮把圆柱/圆柱的 `(r1+r2-dist) <= f64::EPSILON` 订正为 `cxx:1104` 的 `RealSmall()`（`RealEpsilon()` 只用于 `cxx:1394`），并给 `rmr_relative <= 1e-13` 补上 `myEPSILON_CYLINDER_DELTA_RADIUS`（`cxx:356`）与 `A1A2.Intersect()`（`cxx:1224`）的出处注释。门禁：`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（159.4s）、`export_data_obj` 15/15 逐位一致。**仍未移植（下一片）**：锥/锥分支 4「两锥共母线」（`cxx:1759-1886`，`myCommonGen`/`myPChar` + 最大平面内再 reduce 并按 `cxx:1836-1881` 重映射），现走 `None` → 消费方回落一般数值追踪。**该尾项已收口（2026-10-15，T-146 / 包 P016，commit `a01b0bd1`）**：分支 4 整臂落地（顶点互含校验 `cxx:1762-1778`、`aGen`/`myPChar` `cxx:1782-1794`、最大平面法向三支 `cxx:1795-1830`、嵌套 `plane∩cone` 重映射 `cxx:1833-1881`），`common_gen` 由 `cone_cone_analytic` 返回；消费方 `intpatch_impimp_glines::with_cone_common_gen` 按 `IntPatch_ImpImpIntersection.cxx:9149-9171` 追加两锥顶母线，并按 `cxx:8869-9017` 给每条 `IntPatch_GLine` 在参数 `0.` 处挂 `aPChar` 顶点。门禁：`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（143.4s）、`export_data_obj` 15/15 逐位一致。",
  },
];

// `data/*` + `data/occ/*` 逐模型对拍（2026-10-16 `step_obj_gates --nocapture` 实测，23 个模型全过；
// 读数存 `.target-gate/t157/gates-raw.txt`，与 `.target-gate/t156/gates-raw.txt` 逐位相同）
const PARITY = [
  { model: "motoc", occV: 11713, occF: 13944, ourV: 11713, ourF: 13944, area: "1.0000" },
  { model: "a3n00", occV: 11145, occF: 12466, ourV: 11138, ourF: 12442, area: "1.0000" },
  { model: "acs10", occV: 37296, occF: 46558, ourV: 37280, ourF: 46538, area: "1.0000" },
  { model: "ATU01038", occV: 18102, occF: 22586, ourV: 18085, ourF: 22554, area: "1.0000" },
  { model: "bottom", occV: 17587, occF: 23557, ourV: 17581, ourF: 23550, area: "1.0000" },
  { model: "top", occV: 16905, occF: 22800, ourV: 16888, ourF: 22770, area: "1.0000" },
  { model: "T0M", occV: 60548, occF: 67366, ourV: 60199, ourF: 66897, area: "1.0017" },
  { model: "TDB", occV: 73544, occF: 79977, ourV: 73586, ourF: 79994, area: "1.0015" },
];

// 状态/文档偏移：2026-10-16 复核定案（第十至第二十七包逐条复核）
const DOCDRIFT = [
  {
    what: "看板基准日",
    detail: "board asOf 2026-10-11 · r41 → 2026-10-13 · r42 → 2026-10-16 · r43",
    gap: "✅ 已同步：`a91` 已收录 10-12 网格批次与 10-13 缝边注册修复。**本卡（偏移审查）已推进到 2026-10-16 并接管唯一维护中的板**：旧任务看板 `specs/board.canvas.tsx`（此前同步到 `asOf 2026-10-16` / `revision r43`，新增活动记录 `a92` = T-139..T-156 / P009..P026 收口 + T999 `full_outcome_complete: true` + 本轮的 T0M/TDB 面积比订正，门禁数组订正 `G-lib` 1276/0 → **1251/0**、`G-core` 290·215·143·72 → **283·210·143·72**）已于 2026-10-16 退役删除，其门禁表与 `nextAction` 职责并入本画布的 `GATES` / `CYCLES` 两段（板真相另见 `docs/goals/drift-and-gap-autofix/state.yaml`）；未复跑的 `G-fuse` / `G-phase` / `G-boss` / `G-iges` 仍逐条标注「本窗口未复跑」；T-121..T-159 的逐包证据按画布 0.2.0 窗口口径不内联，指向 `specs/_tasks.md`（已在「已归档」区说明）",
    tone: "success" as const,
  },
  {
    what: "G-lib 门禁基线过期",
    detail: "board baseline 1281 / 0 → 1276 / 0 → 1251 / 0",
    gap: "✅ 已同步：2026-10-15 实测 1251 / 0（`--lib`）。`1276` 是 2026-10-13 的读数，之后删除了自创单测并合并了测试文件，读数随之下降；`notes/T002-gate-baseline.md` 已留档",
    tone: "success" as const,
  },
  {
    what: "G-export 门禁基线过期",
    detail: "board baseline 15 / 15 → 23 / 23",
    gap: "✅ 已同步：`-- data/occ data` 口径实测 23 ok / 0 err（15 `data/*.step` + 8 `data/occ/*`）；无参跑是 15 ok / 0 err（2026-10-15 实测，15 个模型 v/f 与包围盒逐位一致）",
    tone: "success" as const,
  },
  {
    what: "G-core 门禁读数过期",
    detail: "canvas 290·215·143·72 → 实测 285·210·143·72",
    gap: "✅ 已同步：2026-10-15 实测 `--lib`（occt-core 285/0、occt-math 210/0、occt-geom 143/0、occt-geom2d 72/0）。旧的 290/215 与 `G-lib` 的 1276 同源——自创单测删除后读数下降，`notes/T002-gate-baseline.md` 已留档",
    tone: "success" as const,
  },
  {
    what: "SEMANTIC 表 LIVE 方向过期",
    detail: "T001 清单登记的自创/采样行其代码已闭环（topo-rest 1/4/6/11 + 18 的 classifier 半、bop 1-4/7-16-17、mesh 1-14/19）",
    gap: "✅ 已订正：2.2 的 SEMANTIC 各行随 P008..P026 逐条换成收口文案（2026-10-13 判词「本表尚未按方向订正」已不成立）。**仍未同步的只剩 `specs/_audit/topo-rest.md` 的 9 / 10 / 15 / 17 四条**：标题仍写 `live：**是**`，代码也确实在（`shape_metrics.rs:146` 曲率对角采样、`validate.rs:150` 恒真 `is_valid`、`hlr.rs` 画家算法、`pipe_along_polyline` 只平移剖面）—— 四者都是无 OCCT 类的 port-only / lite 替代件，不属本 tranche 的 drift 行，故只留档、不改代码",
    tone: "warning" as const,
  },
  {
    what: "未上板的自创件：`bbox_from_geometry.rs` 采样包围盒",
    detail: "64 点边包围盒 + 16×16 UV 面包围盒（`bbox_from_geometry.rs:46,56-67`）",
    gap: "✅ 已闭环（2026-10-15，T-134）：`shape_bbox` 改委托既有忠实件 `brep_bnd_lib::shape_bnd_box`（= `BRepBndLib::Add(S,B,false)`，`BRepBndLib.cxx:81-215`：面走 `BndLib_AddSurface::Add` 并在平面面补边界边 + `Enlarge(Tolerance(F))`，边/顶点走 `BndLib_Add3dCurve::Add`，容差按 `BRep_Tool.cxx:881-895` 抬到 `Precision::Confusion()`）；采样版 `add_to_bbox` 与 `face_uv_bounds` 依赖整体删除，5 处既有 lib 断言按 cxx 盒 `[-t, 1+t]` 订正；已登记 `_audit/topo-rest.md` §19。门禁 1251/0、4/4、15/15 逐位不变",
    tone: "success" as const,
  },
  {
    what: "未提交工作树",
    detail: "`git diff HEAD` 空（0 个已跟踪文件修改）",
    gap: "✅ 已提交：`bc7249fa`（缝边按 OCCT 全 pcurve 注册 + 锥面缝 `FirstParameter` 交换）、`4a20e078`（撤除 3 处 `IGES_UVDBG` 临时插桩，−59 行；并入 goalbuddy 板产物忽略与 IGES 目标文档）",
    tone: "success" as const,
  },
  {
    what: "10-10 未上板批次",
    detail: "brep_class3d：SClassifier::Perform 段循环 / OtherSegment / FindAPointInTheFace",
    gap: "✅ 已上板：`a90b`（2026-10-10，2026-10-13 补记）收录该批次（10 次提交、25 文件 +2,556 / −381、结构 1,825 行），并记明它清掉了 §9.634 E.1/E.2 的 6 条 lib 失败",
    tone: "success" as const,
  },
  {
    what: "_audit §9.3 过期断言",
    detail: "把 geom2d/curve_ops.rs 采样族列为未修",
    gap: "已订正：该采样 2026-09-26（ok267451e8）已删除，改走 Extrema_ExtCC2d / IntAna2d",
    tone: "info" as const,
  },
];

// 看板门禁（读数日期 2026-10-16；G-check / G-lib / G-gates / G-export 覆盖到 P034 已提交工作树 `c272c08a`
// —— G-check / G-lib 走 P034 门禁日志（`occt-geom --all-targets` 0 error = 105 warnings，`occt-topo` 0 error = 413 warnings，
//    `occt-geom --lib` 143/0、`occt-geom2d --lib` 72/0、`occt-core --lib` 283/0、`occt-topo --lib` 1251/0；五 crate 计数沿用 P029 的 12/30/105/82/413），G-gates 走 `.target-gate/t164/gates-raw.txt`（4/4，直跑测试二进制 275.75s：23 面积比与 23 条 ours v=/f= 与 P033 的 t163 逐字一致，仅并发 println 交错顺序与耗时不同）、G-export 15/15 逐位一致（`.target-gate/t164/export-obj.txt`，15 行模型读数与 `t163` 逐字一致，整档仅 `Finished ... in` 耗时行不同），判词与 P026..P033 不变 ——
// G-fuse / G-phase / G-boss / G-iges 为 T-130..T-133 轮读数，其后改动不触及其路径）
const GATES = [
  { id: "G-check", name: "五 crate 编译", base: "exit 0", now: "exit 0", tone: "success" as const },
  { id: "G-lib", name: "topo 单测 --lib", base: "1276 / 0", now: "1251 / 0", tone: "success" as const },
  { id: "G-gates", name: "STEP->OBJ 门禁", base: "4 / 4", now: "4 / 4", tone: "success" as const },
  { id: "G-fuse", name: "fuse 闭环实体", base: "1 / 1", now: "1 / 1", tone: "success" as const },
  { id: "G-phase", name: "phase3-10 / 19 / 20", base: "61 / 0", now: "61 / 0", tone: "success" as const },
  { id: "G-boss", name: "boss 合并", base: "2 / 2", now: "2 / 2", tone: "success" as const },
  { id: "G-export", name: "导出总检", base: "23 / 23", now: "15 ok / 0 err", tone: "success" as const },
  { id: "G-iges", name: "IGES 结构自洽", base: "23 / 23", now: "23 / 23", tone: "success" as const },
  { id: "G-core", name: "core/math/geom/geom2d", base: "290·215·143·72", now: "283·210·143·72", tone: "success" as const },
];

// 目标看板跟踪（`docs/goals/drift-and-gap-autofix/state.yaml`；收口前 active_task 恒为 T005，2026-10-16 随 T999 收口后曾为 T999 / `goal.status: done`，同日后由 P027 复核发现漏审臂而重开为 T005 / `goal.status: active`；P027..P034 逐包在其上滚动，每个包 = 一个已核实并提交的切片）
const CYCLES = [
  {
    pkg: "P012",
    task: "T-142",
    commit: "642a0c88",
    scope: "`intana` / `occt-math` 的第二份自创多项式解算器",
    state: "✅ 闭合：`intana/analytic_intersections.rs` 与 `occt-math/polynomial.rs` 的本地二次/三次/Ferrari 四次 + `poly4`/`polish4` + 自创 Durand-Kerner 删除，改经 `occt-core/src/math_direct_poly_roots.rs` 新增的 `Value(1..NbSolutions)` 适配器委托 `DirectPolynomialRoots`",
  },
  {
    pkg: "P013",
    task: "T-143",
    commit: "8db13d5c",
    scope: "同族解算器的最后四份拷贝（`extrema_pc` / `extrema_cc` / `extrema2d` / `occt-core/numeric.rs`）",
    state: "✅ 闭合：三处改 `pub(super) use occt_core::math_direct_poly_roots::{cubic_roots, quadratic_roots, quartic_roots}`，`numeric.rs` 第四份自创件与其两条测试删除（`occt-core --lib` 285→283）；全仓只剩一条根求解路径",
  },
  {
    pkg: "P014",
    task: "T-144",
    commit: "dd31babc",
    scope: "锥/锥角度容差由调用方常数决定",
    state: "✅ 闭合：`quadric_quadric_cone_cone` 的 `tol_ang` 形参删除，改按 `IntAna_QuadQuadGeo.cxx:1458-1476` 内部推导 `aTolAng`（起点 `Precision::Angular()`，必要时放宽为 `Tol / EstimDist` 并钳制）；`EstimDist`（`cxx:260-326`）全量移植 + 新忠实件 `projlib::project_pln_lin2d`；`cxx:1555` 容差订正",
  },
  {
    pkg: "P015",
    task: "T-145",
    commit: "2da4eaee",
    scope: "锥/锥平行臂的嵌套结果重映射 + 圆柱/圆柱阈值出处",
    state: "✅ 闭合：分支 2 按 `cxx:1557-1603` 重映射嵌套 `plane_cone` 结果（`IntAna_Line` → 同一条线两次，`Parabola`/`Point` → `NoGeometricSolution`）；`(r1+r2-dist)` 改 `RealSmall()`（`cxx:1104`）；同轮订正 audit「圆柱/锥、圆柱/圆柱分类臂未移植」的陈旧登记（逐臂复核实为已移植）",
  },
  {
    pkg: "P016",
    task: "T-146",
    commit: "a01b0bd1",
    scope: "锥/锥分支 4「两锥共母线」（`IntAna_QuadQuadGeo.cxx:1759-1886`）+ 消费方 `IntCoCo` 的 `HasCommonGen` 接线",
    state: "✅ 闭合：`cone_cone_analytic` 新增返回 `common_gen`（`cxx:1782-1794` 的 `aGen` + `myPChar`），分支 4 按 `cxx:1762-1778` 的 `ElSLib::Parameters`/`Value` 顶点互含校验、`cxx:1795-1830` 的最大平面法向三条支路、`cxx:1833-1881` 的嵌套 `plane∩cone` 重映射（无 `IntAna_Line` 臂 ⇒ `None`）落地；消费方 `intpatch_impimp_glines::with_cone_common_gen` 按 `IntPatch_ImpImpIntersection.cxx:9149-9171` 把 `aMkLin(aPApex1, aPApex2)` 生成的母线加入解，并按 `cxx:8869-9017` 给每条二次曲线在参数 `0.` 处挂 `aPChar` 顶点（容差 `Tol`、非多点）。`cxx:1805-1811` 的 `Something wrong ???` 出口 `myCommonGen` 已为真、但其前置的 `IntAna_NoGeometricSolution` 臂（`cxx:9022-9140`）未移植，已在源码注明",
  },
  {
    pkg: "P017",
    task: "T-147",
    commit: "e3a4335e",
    scope: "B6e：`extrema2d::locate_extcc2d` 的自创 Newton 换 `Extrema_LocateExtCC2d` / `Extrema_GFuncExtCC` 2D + `math_FunctionSetRoot`",
    state: "✅ 闭合：`locate_extcc2d`（`curve_curve.rs:255-262`）改调 `ExtremaLocateExtCC2d::new`，`refine_curve_curve2d_newton`（32 次夹取迭代 + `1e-12` 判据）与 `f_and_jac2d`（Jacobian 四个带符号项相对 `Extrema_GFuncExtCC::Values` 全反）删除；新增 `extrema2d/gfunc_ext_cc2d.rs`（`GFuncExtCC2d`：`SearchOfTolerance`、`GetType` 分派、`Value` 的 `D1` 奇点三点 / `DN` 重建、`Values` 的奇点扰动三点 Jacobian、`GetStateNumber`）与 `extrema2d/gen_locate_ext_cc2d.rs`（`GenLocateExtCC2d` = `Extrema_GenLocateExtCC` 的逐轴量程规范化 + `TolU`/`TolV` = `Adaptor::Resolution(Confusion)` + `perform_with_bounds` + `NbExt() > 0` / `Value < 1e-10` 复核；`ExtremaLocateExtCC2d` = `Extrema_LocateExtCC2d.cxx:26-82`）。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（248.5s）、`export_data_obj` 15/15 逐位一致",
  },
  {
    pkg: "P018",
    task: "T-148",
    commit: "8cad6a01",
    scope: "B6g：`extrema2d` 的 `line_line` / `line_circle` / `circle_circle_extrema2d` 三个 `Extrema_ExtElC2d` 臂副本删除",
    state: "✅ 闭合：`crates/occt-geom2d/src/extrema2d/curve_curve.rs` 的三个模块级副本与 `line2d_value`/`circle_param2d`/`dedupe_sort_cc2d` 一并删除（副本已分叉：平行直线多返点、线/圆多返 cxx 不产出的交点、圆/圆多一路自创求交），测试改走唯一忠实件 `extrema2d/ext_el_c2d.rs` 的 `Extrema_ExtElC2d`（`Extrema_ExtElC2d.cxx:63-70` 平行直线只置 `IsParallel()`+`SquareDistance(1)`、`:105-167` 线/圆 ±Y 两对、`:274-368` 圆/圆四对 `P1[NoSol]×P2[kk]` 次序）；五条既有断言按 cxx 计数与次序重写，未新增测试。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（150.7s）、`export_data_obj` 15/15 逐位一致",
  },
  {
    pkg: "P019",
    task: "T-149",
    commit: "616ea0c1",
    scope: "B6e 的 3D 孪生项：`occt-geom/src/extrema_cc/curve_curve.rs::locate_extcc` 的自创 `refine_curve_curve` 换 `Extrema_LocateExtCC` / `Extrema_GFuncExtCC` 3D + `math_FunctionSetRoot`",
    state: "✅ 闭合：`locate_extcc`（`curve_curve.rs:152`）改调 `ExtremaLocateExtCC::new`，`refine_curve_curve`（夹取 32 次迭代 + `1e-12` 判据）与 `f_and_jac`（Jacobian 四个带符号项相对 `Extrema_GFuncExtCC::Values` 全反）从 `poly_roots.rs` 删除；新增 `extrema_cc/p_on_curv.rs`（`POnCurv`＝`Extrema_POnCurv.hxx:25-60`）、`extrema_cc/gfunc_ext_cc.rs`（`GFuncExtCC`＝`Extrema_CCLocFOfLocECC`＝`Extrema_GFuncExtCC` 3D 实例：`SearchOfTolerance`、`GetType` 分派、`Value` 的 `D1` 奇点三点 / `DN` 重建、`Values` 的奇点扰动三点 Jacobian、`GetStateNumber`）与 `extrema_cc/gen_locate_ext_cc.rs`（`GenLocateExtCC`＝`Extrema_LocECC`＝`Extrema_GenLocateExtCC`（`Extrema_GenLocateExtCC.hxx:49-110`：逐轴量程规范化 + `TolU`/`TolV`＝`Curve::resolution(Confusion)` + `perform_with_bounds` + `NbExt() > 0` / `Value < 1e-10` 复核）；`ExtremaLocateExtCC`＝`Extrema_LocateExtCC.cxx:24-74`）。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（181.8s）、`export_data_obj` 15/15 逐位一致",
  },
  {
    pkg: "P020",
    task: "T-150",
    commit: "2452a51c",
    scope: "B6f：`occt-geom2d/src/extrema2d/curve_curve.rs::curve_curve_intersections2d` 的 33×33 采样栅格换 `Geom2dAPI_InterCurveCurve` / `Geom2dInt_GInter`",
    state: "✅ 闭合：`curve_curve_intersections2d` 原实现对两条曲线在 `bound()` 上取 33×33 均匀参数栅格、把距离小于 `tol` 的点对全部推入（另加本地 `1e-6` 去重），OCCT 无此控制流，且交点只要不落在 33 等分格点上就会被静默丢弃或被替换为邻近采样点。现按 `Geom2dAPI_InterCurveCurve.cxx:54-106` 落地：`Geom2dInt_GInter(AC1, AC2, Tol, Tol)`（走已移植的 `Geom2dIntGInter::perform_curves`），未 `IsDone()` 时 `NbPoints() == 0`，随后 `i` 取 `1..=NbPoints()` 的 `Point(i).Value()`。只服务该栅格的 `extrema2d/analytic_solvers.rs::bound()` 一并删除。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（205.1s）、`export_data_obj` 15/15 逐位一致",
  },
  {
    pkg: "P021",
    task: "T-151",
    commit: "11b12e0e",
    scope: "B6d：`occt-geom/src/extrema_ss.rs` 的球/球、球/面解析臂与自创 Newton 换 `Extrema_ExtSS` / `Extrema_ExtElSS` / `Extrema_GenExtSS` / `Extrema_FuncExtSS` 忠实件",
    state: "✅ 闭合：`Extrema_ExtSS::Perform`（`Extrema_ExtSS.cxx:100-255`）只有 `Plane`×`Plane` 一条解析臂（`cxx:122-129` → `Extrema_ExtElSS`），其余类型对一律 `Extrema_GenExtSS(S1, *myS2, 20, 20, tol1, tol2)`（`cxx:131-175` / `cxx:180-221`）并按 `ElCLib::InPeriod` + 落窗判据过滤。删除三条解析臂、只服务它们的辅助件与一般臂的自创 Newton `solve_surface_surface`；新增 `extrema_ss/` 五文件（`ExtremaExtSS` / `ExtremaExtElSS` / `ExtremaGenExtSS` / `ExtremaFuncDistSS`·`ExtremaFuncExtSS` / `tests.rs`）与 `gp_Pln` 的 `signed_distance_*`·`square_distance_pln`（`gp_Pln.hxx:345-379`、`:173-177`）。球/面这对在 OCCT 无解析臂（`Extrema_ExtElSS::Perform(gp_Pln, gp_Sphere)` 抛 `Standard_NotImplemented`），无界平面量程 ±2e100 由 `GeomAdaptor_Surface` 原样透传，故旧断言的 4.0 是被删解析臂编出来的数；`occt_probe --ssx` 重编译实跑钉住真值（`SSX sph-pln nbExt=1 sqDist=7.992002e+200`、`SSX gen-sph-pln nbExt=2 [1] 2.2138509695290867e+198 [2] 7.992002e+200`、`SSX api-sph-pln dist=2.827012911183817e+100`、`SSX sph-sph nbExt=2 sqDist=4.0000000000003126/36`），三条既有断言随之订正（未新增测试）。门禁：`occt-core --lib` 283/0、`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（114.1s）、`export_data_obj` 15/15 逐位一致。同批登记 **B6h**（`surface_surface_extrema` 取首解 vs `GeomAPI_ExtremaSurfaceSurface` 的最小下标 `myIndex`；`point_on_s1`/`point_on_s2` 多暴露私有序列）。B6 家族至此全部收口",
  },
  {
    pkg: "P022",
    task: "T-152",
    commit: "3b5728b3",
    scope: "B7-writer-4：`occt-geom2d/src/geom2d_convert.rs:78` 的第二份残缺 `curve_to_bspline_curve` 副本删除，只留 `Geom2dConvert.cxx:181-449` 单一入口的忠实件",
    state: "✅ 闭合：副本只处理 Bezier（同极点 + `[0,1]` 夹持结点，等价 `cxx:399-419`）与 BSpline（`clone_dyn()`），对全部圆锥、全部裁剪输入与 Offset 一律 `None`，而 cxx 的单一入口 `CurveToBSplineCurve(C, Parameterisation)`（`Geom2dConvert.cxx:181-449`，`Parameterisation` 默认 `Convert_RationalC1`）覆盖裁剪基曲线的圆锥/Bezier/BSpline/Offset（`cxx:189-374`）、非裁剪周期椭圆/圆（`cxx:379-397`）、Bezier（`cxx:399-419`）、BSpline `Copy()`（`cxx:420-423`）、Offset 近似（`cxx:425-440`），其余抛 `Standard_DomainError`（`cxx:442-444`）；同文件 `curve_to_bspline_curve_bspl` 已是忠实件（7 个生产调用点）。副本两条臂与忠实件逐字同构、且缺 `trimmed_basis` 与全部圆锥臂，故删除不改任何可达行为；全仓 grep 该名字只出现在自身定义与 doc（无调用方）。模块头去掉「无参数化入口」的叙述。门禁：`--all-targets` 0 error / 82 warnings（与改前同数）、`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（114.7s）、`export_data_obj` 15/15 逐位一致。旁支登记：2D 侧名 `..._bspl` 与 3D 侧 `convert_bspl::curve_to_bspline_curve`（带 `parameterisation`）不对称，改名需动 `occt-topo` 7 个调用点，超出本包 allowed_files",
  },
  {
    pkg: "P023",
    task: "T-153",
    commit: "38a13eae",
    scope: "B6h：`occt-geom` 的 `GeomAPI_ExtremaSurfaceSurface` 落地。`extrema_ss/mod.rs::surface_surface_extrema` 原取 `mySqDist` 首项，而 OCCT 对应件在 `Init` 里记录最小 `SquareDistance` 的下标 `myIndex`（`GeomAPI_ExtremaSurfaceSurface.cxx:75-89`），`NearestPoints`（`cxx:202-208`）、`LowerDistanceParameters`（`cxx:210-221`）与 `LowerDistance`（`cxx:223-231`）都读它；同处 `ExtremaExtSS::point_on_s1`/`point_on_s2` 按序号暴露 OCCT 的私有序列（`Extrema_ExtSS.hxx:99-102`，对外只有 `Points(N, P1, P2)`，`hxx:92`）",
    state: "✅ 已修（2026-10-15，包 P023，commit `38a13eae`）：`crate::geom_api::ExtremaSurfaceSurface` 按 `cxx:57-231` 落地（`Init` 用 `PConfusion` 建 `Extrema_ExtSS`、`myIsDone = IsDone && NbExt() > 0`、`myIndex` 全表最小下标、`NbExtrema`/`Points`/`Parameters`/`Distance`/`IsParallel`/`NearestPoints`/`LowerDistanceParameters`/`LowerDistance`/`Extrema` 逐条对应；10 参 `Init` 与构造无调用方，标注未移植）。消费方 `extrema::surface_surface_extrema` 改由该类组对（即最小解），首解门面 `extrema_ss::surface_surface_extrema` 删除，`ExtremaExtSS::point_on_s1`/`point_on_s2` 删除（grep 确认无调用方）。可达性：`Extrema_GenExtSS::Perform` 先入 BFGS 最小解（`Extrema_GenExtSS.cxx:377-386`）、`SR2` 再追加最大解（`cxx:410-414`）、`Plane`×`Plane` 只有 1 条（`Extrema_ExtSS.cxx:229-233`），probe 读数亦以首项为最小，故四门禁数值按构造不变，闭合的是结构缺口。门禁：`--all-targets` 0 error / 105 warnings（同改前）、`occt-geom --lib` 143/0、`occt-geom2d --lib` 72/0、`occt-topo --lib` 1251/0、`phase6_integration` 5/0、`step_obj_gates` 4/4（115.6s）、`export_data_obj` 15/15 逐位一致。",
  },
  {
    pkg: "P024",
    task: "T-154",
    commit: "549e5761",
    scope: "`BRepLib_FindSurface` / `BRepLib_MakeFace(W, OnlyPlane)` 落地：`brep_offset` 的 `offset_wire_3d` 原用「每边采 4 点」的自创最小二乘拟合求平面，而其 OCCT 对应件是 `BRepOffsetAPI_MakeOffset::BuildDomains`（`BRepOffsetAPI_MakeOffset.cxx:184-192`）→ `BRepBuilderAPI_MakeFace(WorkWires.First(), true)`（`BRepBuilderAPI_MakeFace.cxx:215-224`）→ `BRepLib_MakeFace(W, OnlyPlane)`（`BRepLib_MakeFace.cxx:188-256`）→ `BRepLib_FindSurface(W, -1, OnlyPlane, true)`（`BRepLib_FindSurface.cxx:52-628`）",
    state: "✅ 已修（2026-10-16，包 P024）：(1) **Scout 复核订正两处误标** —— ① `GeomLib::IsPlanarSurface` / `GeomLib::Inertia` / `GeomLib::AxeOfInertia` / `Controle` 在仓内既无对应件也无可达分支（`grep -rn \"IsPlanarSurface\\|AxeOfInertia\" crates/**/*.rs` 零命中；命中的 `Inertia`/`Controle` 是 `BRepGProp_Gauss` 累加器、`fillet_var` 注释与 `math_NewtonMinimum.cxx:230`），故这一切片不是「补 `Controle`/`GeomLib::Inertia`」，而是「按消费方真实调用链移植 `BRepLib_FindSurface`」；② `occt-geom/src/surface_fit.rs` 的 `fit_plane` 在改写后已无生产调用方（只剩自带测试）。(2) **新增忠实件** `crates/occt-topo/src/brep_lib_find_surface.rs`（588 行）＝`BRepLib_FindSurface.cxx:52-628`：`Controle`（`cxx:52-68`，点到面最大距离）、`Is2DConnected`（`cxx:75`）、`Is2DClosed`（`cxx:107`）、`fillParams`（`cxx:178`）、`fillPoints`（`cxx:217`）、`Init`（`cxx:246`）与 `Found`/`Surface`/`Tolerance`/`ToleranceReached`/`Existed`/`Location`（`cxx:583-621`）；先按 `CurveOnSurface` 复用既有曲面（`OnlyClosed` 时过 `Is2DClosed`），否则按边曲线类型取点（Bezier/B-spline 结点、直线 2 点、圆锥 4 点、其余 `15 + NbIntervals(GeomAbs_C3)`），按 `W_i = |P_{i-1}P_i| + |P_iP_{i+1}|` 加权求重心与协方差，`math_Jacobi` 取最小特征值对应特征向量（奇异时取另两轴叉积）拟合 `gp_Pln`，`Controle` 判 `tol_reached`，再用临时面的 `FClass2d::perform_infinite_point() == In` 定向法向。(3) **新增 `brep_lib_make_face.rs::make_face_from_wire`**＝`BRepLib_MakeFace.cxx:188-256`：跑 `FindSurface(W, -1, OnlyPlane, true)`，`tol = max(1.2 * ToleranceReached, Tolerance)`（`cxx:200`），`OnlyPlane` 时按 `cxx:211-242` 剔除退化边（保持线框朝向与闭合），加线建面；未移植项（`UpdateTolerances`/`SameParameter` `cxx:250-252`、`CheckInside` `cxx:255`）在源码按行号标注。(4) **新增 `occt-core/src/gp/ax3.rs::GpAx3::from_point_dir`**＝`gp_Ax3.cxx:29-80`（按规范坐标分支选正交向量）。(5) **接线**：`brep_offset/curve_face_offset.rs::offset_wire_3d` 改调 `make_face_from_wire(wire, true)` + `face_plane(&face)`，自创拟合不再在 `brep_offset` 路径上。门禁：`cargo check --manifest-path crates/occt-topo/Cargo.toml` 0 error、`occt-core --lib` 283/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 逐位一致。旁支登记：`occt-geom/src/surface_fit.rs` 的 `fit_plane`/`fit_surface_kind`/`fit_surface` 为 port-only 且已无生产调用方，是否删除留待后续裁决",
  },
  {
    pkg: "P025",
    task: "T-155",
    commit: "8d77b13a",
    scope: "B7-mesh_to_brep：`occt-topo/src/mesh_to_brep.rs` 原先按自创规则建 B-Rep（顶点按 `1e-9` 位置哈希归并、`Vertex` 容差写 `0.0`、每三角形手工按 `(pb-pa) x (pc-pa)` 定面法向并逐边定向），换成 `BRepBuilderAPI_MakeShapeOnMesh`（`BRepBuilderAPI_MakeShapeOnMesh.cxx:74-263`）忠实件",
    state: "✅ 已修（2026-10-16，包 P025）：新增 `MakeShapeOnMesh`＝`BRepBuilderAPI_MakeShapeOnMesh.cxx:74-263` 全段 —— 空 mesh 的 not-done 出口（`cxx:80-85`）、节点 1:1 建点于 `BRepLib::Precision()`（`cxx:96-111`，`BRepLib_MakeVertex.cxx:26-31`、`BRepLib.cxx:84`，即 `Precision::Confusion()` = 1e-7）、按规范节点下标对登记共享边（`cxx:31-55`、`cxx:150-179`，`NCollection_IndexedDataMap::Add` 首次值语义）、下标退化与 `gp::Resolution()` 退化三角形跳过（`cxx:126-141`）、每三角形加平面面并按 `cxx:213-228` 逐边定向、`BRepAdaptor_Curve::Line()` 的 `d1 x d2` 法向按 `Precision::SquareConfusion()` 判非退化（`cxx:238-247`）、过 `Node(anIdx[0])` 建面（`cxx:255`）、收尾 `Done()` + `TopoDS_Compound`（`cxx:182-184`、`cxx:257-263`）。`cxx:248-253` 的 `aN.Reverse()` 折入朝向（存储边恒为规范 `lo -> hi`，面朝向承载三角形走向）。`triangulation_to_brep` 只保留 `LOCAL:` 标记的扩展（`weld_nodes`、闭合 2-manifold 判据、`OrientClosedSolid`、散度定理绕向）。`Vertex` 容差改 `CONFUSION` 是承重的：`pave_ff.rs:158`、`pave_split_blocks.rs:74`、`brep_fill.rs:213` 以 `max(...)` 把顶点容差折进工作精度。门禁：`cargo check --manifest-path crates/occt-topo/Cargo.toml --offline` 0 error、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（115.7s）、`export_data_obj` 15/15 逐位一致；未新增测试、未改既有断言。旁支登记：`mesh_to_brep` 仍被 `brepfeat/features.rs:367`（`shape_mesh_to_brep`）与 `rwmesh/mod.rs:24` / `rwmesh/scene.rs:307`（`triangulation_to_brep`）消费，这些入口及 `LOCAL:` 焊接/闭合判据/定向无 OCCT 类；另订正登记条目的指针 —— `_audit/bop-boolean.md` 第 2 项记的 `bop_builder_core.rs:85-87` `32^3 voxel_boolean` 兜底在当前树已不存在（该处现为仅告警的网格体积后检，`voxel_boolean` 只剩 `boolean_ops.rs` 自身与 `bop_builder_planar.rs:47` 注释），该审计条目待复核",
  },
  {
    pkg: "P026",
    task: "T-156",
    commit: "fbd59532",
    scope: "`ReparamCurve2d` / `GeomLib::SameRange` 转换尾去重：`crates/occt-topo/src/pcurve_full/surface_projector/reparam.rs` 与 `shhealing/wire_fix/pcurve_range.rs` 各带一份「`Geom2d_TrimmedCurve` → `Geom2dConvert::CurveToBSplineCurve` → `BSplCLib::Reparametrize` → `SetKnots`」转换尾副本（`GeomLib.cxx:912-921` / `:962-968`），且 `pcurve_range.rs` 的副本在 cxx 不判型的不等跨臂（`cxx:924-969`）按曲线类型分派（`line_trim_to_reparam_bspline` / `unequal_trim_window` → `same_range_unequal_line` / `same_range_unequal_other`）",
    state: "✅ 已收口（2026-10-16，包 P026，commit `fbd59532`）。判词先行：该文件此前已委托曲线 kind 且 `reparam_curve2d` 已返回具体 `Geom2dBSplineCurve`，故**不是** stale-closed —— 真实缺口是两份转换尾副本的重复与漂移。<br>(1) 两份合并为唯一忠实件 `trimmed_reparam_bspline`（`GeomLib.cxx:912-921` + `:962-968`：`Geom2dTrimmedCurve(CurvePtr, U1, U2)` → `Geom2dConvert::CurveToBSplineCurve` → 对 `BS->Knots()` 做 `BSplCLib::Reparametrize(RequestedFirst, RequestedLast, Knots)` → `BS->SetKnots(Knots)`）与 `same_range_window`（`cxx:934-959`：`aCCheck` 周期判据 `:929-932`、周期窗 `:936-938` / 回落 `:940-943`、`Udeb/Ufin` 夹取 `:948-949`、退化回落 `:950-957`）。<br>(2) 删除的重复件：`same_range_unequal_line` / `same_range_unequal_other` / `trim_to_reparam_bspline` / `line_trim_to_reparam_bspline` / `unequal_trim_window`（cxx 的不等跨臂对直线、圆锥、B 样条走同一条路，没有判型分支），以及 `ReparamCurve2d::bspline_copy2d`（wrapper 只在 `Geom2dConvert::CurveToBSplineCurve` 拒绝时存在，`Geom2dConvert.cxx:442-444`，故 trait 默认 `None` 才是正确答案）。`pcurve_range.rs` 只保留 cxx 的分支选择（早退 `cxx:854-859`、等跨判据 `:862`、直线 `:864-871`、圆 `:872-889`、裁剪递归 `:890-901`、一般臂 `:908-922`、不等跨 `:924-969`），`wire_fix/mod.rs` 的未用 `reparameterize` 导入与 `lacking.rs` 指向已删件的注释同步清理，全部引用行号按 grep 物理行号订正。<br>(3) 门禁：`cargo check --manifest-path crates/occt-topo/Cargo.toml --offline --all-targets` **0 error / 413 warnings**、lib **0 error / 375 warnings**（与基线同数）；`occt-topo --lib` **1251/0**（等于基线，单调不降）；`step_obj_gates` **4/4**（118.1s，23 个模型面积比 1.0000，Shape-2 0.9999 / T0M 1.0017 / TDB 1.0015 均在既有容差内）；`export_data_obj` **15/15 ok** 且逐模型 v/f 与包围盒与 T002 冻结表逐位一致。未新增测试、未改既有断言。",
  },
  {
    pkg: "P027",
    task: "T-157",
    commit: "76893533",
    scope: "B8-intana 的**平面/平面臂**（P026 的 T999 只核了锥/锥与柱族臂，平面臂与 `IntAna_Int3Pln` 臂从未与 cxx 对拍）：`IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Pln, TolAng, Tol)`（`IntAna_QuadQuadGeo.cxx:389-511`）与 `IntAna_Int3Pln::Perform`（`IntAna_Int3Pln.cxx:40-64`）",
    state: "✅ 已收口（2026-10-16，包 P027，commit `76893533`）。三处自创 + 一段未移植：<br>(1) `plane_plane_line`（`analytic_intersections.rs:97`）的交点原点用 cxx 中不存在的最小范数公式 `p = (d1·(n2×v) + d2·(v×n1))/|v|²`（`v = n1×n2`）配 `v2 < 1e-20` 叉积长度门 —— cxx 是 `aEps = 1.e-16` 的分母钳制 `denom = (|1 - (A1A2+B1B2+C1C2)²| <= aEps) ? aEps : ddenom` 加两次法向投影 `inter1 = aVN1×vd` / `inter2 = aVN2×vd` 的中点（`cxx:423-445`），落地为常量 `A_EPS`（`:156`）。<br>(2) `cxx:449-508` 的**小角原点重解整段未移植**：`aMVD < aTreshAng = 2.e-6` 且 `|aDist1| > aTreshDist = 1.e-12 || |aDist2| > aTreshDist`（常量 `A_TRESH_ANG` `:158` / `A_TRESH_DIST` `:160`）时用两次 `IntAna_IntConicQuad::Perform(gp_Lin, gp_Pln, TolAng, Tol)` 重算 `pt1`（`aL1 = gp_Lin(pt1, gp_Dir(aVN1))` 对 P1 `cxx:477-487`；`aDL2 = gp_Dir(dir1.Crossed(aDN1))`、`aL2 = gp_Lin(aPnt1, aDL2)` 对 P2 `cxx:489-506`）；`!IsDone()` 或 `IsParallel()` 时 cxx 直接 `return`（`done = false` ⇒ 本仓 `None`）。<br>(3) `solve3`（`:71`）是手写部分选主元消元加 `1e-14` 主元门，cxx 走 `math_Gauss Resol(M, gp::Resolution())`（`IntAna_Int3Pln.cxx:51`）—— 现调已移植的 `occt_math::Gauss::new(&M, REAL_SMALL)` 并按 `cxx:53-61` 判 done / `V = -V` / `Solve`。<br>(4) 自创的 `Intersect3Pln::Line` 三平面共线分支（`v2 < 1e-20` + `|a·q + d| < 1e-9`）删除 —— `IntAna_Int3Pln.hxx:29-70` 只有点/空两态，`three_planes_intersect` 只剩 `Point`/`NoResult`。<br>调用面：`plane_plane_line` 多收 `(tol_ang, tol)`，`cone_cone_common_apex` 内嵌套平面/平面求解按 `cxx:1709` 传 `(Tol, Tol)`。<br>门禁：五 crate `--all-targets` 0 error；`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0；`step_obj_gates` 4/4（251.49s，23 面积比与 P026 逐位一致）；`export_data_obj` 15/15 逐位一致；未新增测试。<br>旁支登记（未改）：`quadric_quadric_plane_sphere`（`:266`）的相切判据 `(dist.abs() - radius).abs() < radius.abs() * f64::EPSILON` 与 cxx 的 `std::abs(std::abs(dist) - radius) < Epsilon(radius)`（`IntAna_QuadQuadGeo.cxx:1000`）不符；球/球臂（`cxx:2034-2163`）的 `0.01 * CONFUSION` 圆判定待核。",
  },
  {
    pkg: "P028",
    task: "T-158",
    commit: "1d595f8d",
    scope: "`IntAna_QuadQuadGeo` 的两条 P027 旁支：`quadric_quadric_plane_sphere`（`IntAna_QuadQuadGeo.cxx:980-1042`）的相切判据，与 `quadric_quadric_cylinder_sphere`（`cxx:1373-1405`）/ `quadric_quadric_sphere_cone`（`cxx:1917-2005`）的球心-轴判定；后者在 `occt-topo/src/intpatch_impimp_pairs.rs::int_cy_sp` / `int_co_sp` 另有两份副本",
    state: "✅ 已收口（2026-10-16，包 P028，commit `1d595f8d`）。<br>(1) **`Epsilon(radius)`**：`quadric_quadric_plane_sphere` 的相切判据由 `(dist.abs() - radius).abs() < radius.abs() * f64::EPSILON` 改为 `(dist.abs() - radius).abs() < epsilon(radius)`，即 cxx:997 的 `std::abs(std::abs(dist) - radius) < Epsilon(radius)`；`Epsilon(x) = nextafter(x, +inf) - x`（`Standard_Real.hxx:242-246`，本仓 `occt_core::precision::epsilon`，经 `intana::prelude` 重导出）与 `|x| * DBL_EPSILON` 只在 `|x| in [1,2)` 相等，球半径 5/10 上相差 4~8 倍。该臂其余控制流（`|dist| < radius` 圆交、`== radius` 切点、`|dist| > radius` 空，`cxx:984-1042`）逐行一致，未动。<br>(2) **新增 `axes_meet_at_point(a1, a2, pt)`**（`intana/line_torus.rs`，`pub`）：逐行实现 cxx:1378-1379 / :1924-1925 的 `(A1A2.Intersect() && Pt.Distance(A1A2.PtIntersect()) == 0.0) || A1A2.Same()` —— `Intersect() = thecoplanar && !theparallel`（`cxx:87`）、`Same() = theparallel && (thedistance < myEPSILON_DISTANCE)`（`cxx:83`）、默认 `theEpsDistance = 1.e-14` / `theEpsAxesPara = Precision::Angular()`（`cxx:72-73`、`:124-131`），复用已移植的 `axe_operator`；`Pt.Distance(...) == 0.0` 是 OCCT 的精确比较，本实现同样不加容差。两臂改调它，`dist_point_axis` 因无调用方而删除。<br>(3) **occt-topo 侧两份副本**：`intpatch_impimp_pairs.rs::int_cy_sp` / `int_co_sp` 的自创 `dist_point_axis(...) > 1e-9` 门改为 `!axes_meet_at_point(...)`（`IntPatch_ImpImpIntersection.cxx:8137` / `:9208` 是无条件调 `IntAna_QuadQuadGeo`，失败才落 `cxx:8263-8363` 的 `IntAna_NoGeometricSolution` 臂），回落臂保持原样；`intpatch_impimp_quad.rs` 的 `dist_point_axis` 同步删除。<br>(4) **核对结论（无改动）**：球/球臂的 `0.01 * CONFUSION` **不是**漂移 —— `myEPSILON_MINI_CIRCLE_RADIUS` 在 `InitTolerances()` 里即 `0.01 * Precision::Confusion()`（`cxx:355`），该臂逐行核对一致（含 `d <= Tol` 时 `done` 保持 true、`typeres` 保持 `IntAna_Empty` 的空返回语义）。<br>门禁：五 crate `--all-targets` **0 error**（occt-core 12 / occt-math 30 / occt-geom 105 / occt-geom2d 82 / occt-topo 413 warnings）；`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0；`step_obj_gates` 4/4（23 面积比与 P027 逐位一致）；`export_data_obj` 15/15 逐位一致；未新增测试、未改既有断言。<br>**逐臂普查（本节给出的事实）**：`IntAna_QuadQuadGeo` 共 15 个 `Perform` 重载，已移植 10 个（`cxx:389` Pln/Pln、`:543` Pln/Cylinder、`:752` Pln/Cone、`:980` Pln/Sphere、`:1050` Cylinder/Cylinder、`:1324` Cylinder/Cone、`:1373` Cylinder/Sphere、`:1433` Cone/Cone、`:1917` Sphere/Cone、`:2034` Sphere/Sphere），其中仅 5 个做过逐行 cxx 对拍；**未移植 5 个，全部是 Torus 家族**（`cxx:2163` Pln/Torus、`:2278` Cylinder/Torus、`:2357` Cone/Torus、`:2496` Sphere/Torus、`:2588` Torus/Torus 在 `occt-geom/src/intana` 里无对应函数），须由下一片从消费方判定是可达缺口还是无消费者的未移植标记，不得用替代件冒充。",
  },
  {
    pkg: "P029",
    task: "T-159",
    commit: "6330c6d7",
    scope: "`IntAna_QuadQuadGeo` 环面五臂逐行核对（`cxx:2163-2276` Pln/Torus、`:2278-2355` Cylinder/Torus、`:2357-2494` Cone/Torus、`:2496-2586` Sphere/Torus、`:2588-2666` Torus/Torus），并订正 P028 普查「Torus 族未移植」的错误结论",
    state: "✅ 已收口（2026-10-16，包 P029，commit `6330c6d7`）。<br>(1) **普查订正**：P028 的「已移植 10 / 未移植 5，全部是 Torus」不成立 —— 那次只列了 `crates/occt-geom/src/intana/` 目录，而 `crates/occt-geom/src/intana/mod.rs:30-35` 用 `#[path = \"../intana_torus.rs\"] mod torus;` + `pub use torus::*;` 把 `crates/occt-geom/src/intana_torus.rs` 并入同一命名空间。五臂齐备且各有生产调用方：`quadric_quadric_plane_torus`（`occt-topo/src/intpatch_impimp_pairs.rs:117`）、`quadric_quadric_cylinder_torus`（`:155`）、`quadric_quadric_cone_torus`（`:162`）、`quadric_quadric_sphere_torus`（`:169`，另 `occt-topo/src/intpatch.rs:576`）、`quadric_quadric_torus_torus`（`:175`）。**正确普查：15 个 `Perform` 重载全部在仓**。方法教训写进 `notes/T001-inventory.md`：点名某个 OCCT 类在仓内的存在性要按 **`mod` 声明**枚举（可能带 `#[path]`、指向目录外），不能按目录列举。<br>(2) **真自创（出现两次，已删）**：`quadric_quadric_sphere_torus`（`cxx:2496-2562`）与 `quadric_quadric_torus_torus`（`cxx:2588-2666`）各带 `if a_dist <= f64::EPSILON { if |r1-r2| > tol { Empty } else { Circles(<单圆>) } }`。cxx 无此段：`IntAna_Empty` 早退（`cxx:2529` / `cxx:2636`）后即 `anAlpha = 0.5*(r1*r1 - r2*r2 + aDist*aDist)/aDist`（`cxx:2534` / `cxx:2640`），随后**无条件**构造 `gp_Dir aDir12(aVec12)`（`cxx:2542` / `cxx:2647`），零向量抛 `Standard_ConstructionError`。`f64::EPSILON` 即 `RealEpsilon()`，这两臂实际用的是 `myEPSILON_DISTANCE`（`cxx:353`）与 `myEPSILON_CYLINDER_DELTA_RADIUS`（`cxx:356`）；自创的同心圆答案与阈值一并删除，保留已有 `GpDir::from_vec(..)` 的 `Err` 分支作为 OCCT 异常的替身并注释 `cxx:2542` / `cxx:2647`。<br>(3) **同轮复核（无改动）**：`EPS_DISTANCE = 1e-14`（`cxx:353`）、`EPS_CYL_DELTA_RADIUS = 1e-13`（`cxx:356`）、`dirs_parallel`/`dirs_normal` = `gp_Dir::IsParallel`/`IsNormal`；平面臂 `aTolNum` 判据（`cxx:2163-2276`）与锥臂双解循环（`cxx:2357-2494`）逐行一致。<br>门禁：五 crate `--all-targets` **0 error**（12/30/105/82/413 warnings）；`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0；`step_obj_gates` 4/4（151.50s，23 面积比与 P027/P028 逐位一致）；`export_data_obj` 15/15 逐位一致；未新增测试、未改既有断言。<br>**剩余未知（下一片 P030）**：仍未逐行对拍的重载 `cxx:543-745` Pln/Cylinder、`:752-972` Pln/Cone、`:1050-1293` Cylinder/Cylinder、`:1324-1342` Cylinder/Cone、`:1433-1886` Cone/Cone；P015 关于柱/柱与柱/锥的笔记只是引文核对。",
  },
  {
    pkg: "P030",
    task: "T-160",
    commit: "0a185eb5",
    scope: "`IntAna_QuadQuadGeo` 余下五臂逐行对拍（`cxx:543-745` Pln/Cylinder、`:752-972` Pln/Cone、`:1050-1293` Cylinder/Cylinder、`:1324-1342` Cylinder/Cone、`:1433-1886` Cone/Cone），对照 `quadric_quadric_plane_cylinder` / `quadric_quadric_plane_cone`（`intana/analytic_intersections.rs`）与 `quadric_quadric_cylinder_cylinder` / `quadric_quadric_cylinder_cone` / `quadric_quadric_cone_cone`（`intana/line_torus.rs`）",
    state: "✅ 已收口（2026-10-16，包 P030，commit `0a185eb5`）。两处缺口一修一标，另三臂核对为忠实。<br>(1) **Pln/Cylinder 缺两段，已修**：`quadric_quadric_plane_cylinder` 的签名没有 `H`，旧注释自认 `(H = 0, so the near-parallel direction-refinement branch is skipped)`，即把 `cxx:568-598` 整段声明为省略；并行与否也只用内联的 `|axis . normal| < Tolang`，而 cxx 用 `IntAna_IntConicQuad(axec, P, tolang, toltang, H).IsParallel()`（`cxx:601`）—— 后者在 `H != 0` 时还要过 `Len` 端点距离检查，两者并不同义。现新增 `h: f64`（= `cxx:546` 的 `H`，`hxx:94` 默认 0），移植 `cxx:568-598`（`dA > pi/4`、`dangle = abs(dA - pi/2) > Tolang`、`sinda = abs(sin(dangle))`、`dif = abs(sinda - Tol) < Tol` 或 `H > 0 && sinda * H < 2 * Tol` 时放宽 `tolang = sinda * 2`、`toltang = max(Tol, sinda * H * 1.01)` 并置 `newparams`）；并行判定改由忠实件 `IntAnaIntConicQuad::line_plane(&axec_lin, p, tolang, toltang, h)` 回答（`cxx:559` 的 `gp_Lin axec(Cl.Axis())` 对 `cxx:601`）；`newparams` 的两条重算方向臂按 `cxx:617-631`（切线单直线）与 `cxx:647-668`（双直线，含 `cxx:654-655` 的 `ht = max(arg, 0).sqrt()` 钳位）落地；非并行臂圆心改取 `inter.point(1)`（`cxx:687`）。`gp_Dir` 的 `Standard_ConstructionError`（`cxx:629` / `:659` / `:665`）由既有 `GpDir::from_vec(..)` 的 `Err` 分支承担。<br>(2) **消费方 `IntPCy` 同步接面高**：`occt-topo/src/intpatch_impimp_pairs.rs::int_p_cy` 原先不传 `H`（等于恒按 `H = 0`），现按 `IntPatch_ImpImpIntersection.cxx:2579-2587` 取 `aSCyl`（由 `Reversed` 决定是 `S1` 还是 `S2`）的 `VMax - VMin`，任一端无限时为 `0`；其余 6 处调用点显式传 `0.0`（`hxx:94` 默认实参）并注明理由。<br>(3) **Pln/Cone 尾部一段未移植，按门禁第 2 条标 `UNPORTED`（未改行为）**：`cxx:922-951` 的 `EllipseLimit = 1.0E+9`（`cxx:926,929`）与 `HyperbolaLimit = 2.0E+6`（`cxx:935-950`）回绝段置 `done = false` 并返回，消费方翻成 `return false`（`IntPatch_ImpImpIntersection.cxx:3480-3483`）= 本仓 `PairOutcome::Fail`；`cone_cone_analytic` 分支 2 亦用 `INTER_QUAD_PLN.IsDone()`（`cxx:1556`）判该嵌套调用。本仓 `QuadricIntersection` 无「未完成」态，映射成 `None` 等于 `IntAna_Empty`（语义相反），故只标行号，理由写在 `quadric_quadric_plane_cone` 的函数文档里。<br>(4) **三臂逐行核对为忠实**：`cxx:1050-1293`（柱/柱：`AxeOperator` 容差、五条并行臂含 `RealSmall()` `cxx:1104`、`RmR_Relative <= myEPSILON_CYLINDER_DELTA_RADIUS`、等半径双椭圆参数与交换、外切点）、`cxx:1324-1342`（柱/锥：默认 `AxeOperator` 实参 `1.e-14`/`Precision::Angular()`、`Same()` → `apex +/- r/tan(semiAngle)` 两圆、否则 `NoGeometricSolution`）、`cxx:1433-1886`（锥/锥：容差块 `cxx:1436-1476` 与四分支入口 `cxx:1478`/`:1524`/`:1608`/`:1760` 复核，内部由 P014/P015/P016 收口）。<br>门禁：五 crate `--all-targets` **0 error**（105/413 warnings，与 T002 同数）；`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0；`step_obj_gates` 4/4（131.26s，`.target-gate/t160/gates-raw.txt`，23 面积比与 P029 逐位一致）；`export_data_obj` 15/15 逐位一致；未新增测试、未改既有断言。<br>**诚实边界**：新增的 `H` 臂在 15 个导出模型与 23 个对拍模型里都没有覆盖（需要一个轴略偏离平面法向、面高很小的柱面），门禁只证明无回归，不证明该臂；`int_face_face_analytic.rs` 的两条自设解析快速路径仍按 `H = 0` 并已注释。<br>**下一片 P031**：给 `QuadricIntersection` 补 declined 态并翻译 `cxx:922-951`。",
  },
  {
    pkg: "P031",
    task: "T-161",
    commit: "9bef99ad",
    scope: "`IntAna_QuadQuadGeo` 平面/锥臂的极端圆锥 tail 落地（`IntAna_QuadQuadGeo.cxx:922-951`）：给 `occt_geom::intana::QuadricIntersection` 加 declined 态 `NotDone`，按 cxx 翻译 `EllipseLimit = 1.0E+9`（`cxx:926-934`）与 `HyperbolaLimit = 2.0E+6`（`cxx:935-951`）的回绝，并接上三条消费路径（`occt-topo/src/intpatch_impimp_glines.rs::from_quadric_intersection`、`intana/line_torus.rs` 的两处嵌套平面/锥约化、`int_face_face_analytic.rs` / `intpatch_analytic.rs` 的解析快速路径）",
    state: "✅ 已收口（2026-10-16，包 P031，commit `9bef99ad`）。P030 留下的 `UNPORTED cxx:922-951` 已翻译，靠新增的枚举态。<br>(1) **declined 态**：`QuadricIntersection`（`crates/occt-geom/src/intana/analytic_intersections.rs:188-210`）增 `NotDone` 变体，文档注释写明它即 `done = false`（调用方必须改用双参数方法，而不是读一条曲线），`PartialEq` 与自身相等；cxx 的两条判据按原值落成具名常量 `ELLIPSE_LIMIT = 1.0e9`（`cxx:927` 的 `EllipseLimit`）与 `HYPERBOLA_LIMIT = 2.0e6`（`cxx:928` 的 `HyperbolaLimit`）。<br>(2) **tail 口径**：`quadric_quadric_plane_cone` 由 `if/else if/else` 改成先建 `res` 再判 —— 椭圆臂 `|param1| > EllipseLimit || |param1bis| > EllipseLimit`（`cxx:926-934`）、双曲线臂 `|param2| / |param2bis| > HyperbolaLimit`（`cxx:935-942`）与 `|param1| / |param1bis| > HyperbolaLimit`（`cxx:943-951`）命中即 `NotDone`。本仓没有 `param1/param2` 这些局部量：椭圆臂存的是 `GpElips::new(pos, param1 /*major*/, param1bis /*minor*/)`，两条双曲线臂存的是 `GpHypr::new(pos, param1 = param2, param1bis = param2bis)`（`cxx:850-851` 的轴含平面臂与 `cxx:901-902` 的 `cost < sina` 臂都置 `nbint = 2` 且两对参数相等），故 `cxx:935-942` 与 `:943-951` 在本题上等价，合并成一次 `major/minor` 判据，理由写进函数文档与行内注释；回绝只可能在椭圆/双曲线臂触发（顶点平面臂只产 `Line`/`TwoLines`/`Point`），所以「在返回的曲线上判」与 cxx「在函数尾判」等价。<br>(3) **消费路径**：`from_quadric_intersection` 增 `NotDone => PairOutcome::Fail`，即 `IntPatch_ImpImpIntersection.cxx:3480-3483` 的 `if (!inter.IsDone()) return false`（`PairOutcome::Fail` 在 `intpatch_impimp.rs:115` 不置 `status`、不置 `empty`，整对曲面交给双参数方法）；`line_torus.rs::map_plane_cone_common_axis` 的分支 2 显式增 `NotDone => None`、分支 4 沿用 `_ => None`，因为 cxx 两处都用 `if (INTER_QUAD_PLN.IsDone())` 把重映射 switch 整体括住（`cxx:1556`、`:1835`），嵌套未完成时不进 switch，`typeres` 保持对象初值 `IntAna_Empty`（`cxx:960`）= 本仓 `None`；本仓自设的两条解析快速路径 `int_face_face_analytic.rs::conics_to_curves` 与 `intpatch_analytic.rs::ic_list_from_quadric` 把 `NotDone` 归入「无曲线」臂，`plane_cone` 随之回落 `self.general`（数值追踪 = cxx 的双参数方法）。<br>门禁：五 crate `--all-targets` **0 error**（12/30/105/82/413 warnings，与 T002 同数，无新增告警）；`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0；`step_obj_gates` 4/4（直接跑测试二进制 251.50s，`.target-gate/t161/gates-raw.txt`，23 面积比与 23 条 `ours v=/f=` 行与 P030 的 `.target-gate/t160/gates-raw.txt` 逐字一致）；`export_data_obj` 15/15 逐位一致；未新增测试、未改既有断言。<br>**诚实边界**：新 tail 在 15 个导出模型与 23 个对拍模型上都没有触发（触发需要平面几乎与母线平行、且截面参数超过 1e9/2e6 的锥对），四道门禁只证明无回归，不证明该 tail 本身 —— 与 P030 的 `H` 臂同性质。<br>**下一片 P032**：`IntAna_IntLinTorus::Perform`（`IntAna_IntLinTorus.cxx:88-129`）的两处 `done = false` —— `math_DirectPolynomialRoots` 未完成（`cxx:124-129`）与「有根但全部被回代判为坏解」（`cxx:116-120` 的 `n > 0 && nbsolvalid == 0 && aNbBadSol == n`）都表达不了：本仓 `line_torus_intersect` 返回 `Vec<GpPnt>`，而便利包装 `occt_core::math_direct_poly_roots::quartic_roots`（`math_direct_poly_roots.rs:724-726`）丢掉了 `DirectPolynomialRoots::is_done()`；消费方 `crates/occt-topo/src/intcurvesurface/solvers.rs:456` 的 cxx 对照点是 `IntCurveSurface_InterUtils.pxx:1290-1294` 的 `if (!intlintorus.IsDone()) return false;`。见 audit `occt-geom.md` §9。",
  },
  {
    pkg: "P032",
    task: "T-162",
    commit: "55b0a6d6",
    scope: "`IntAna_IntLinTorus` 的 `done` 语义（`IntAna_IntLinTorus.cxx:88-133`）：`line_torus_intersect`（`crates/occt-geom/src/intana/line_torus.rs`）改返回携带 `done` 的 `LineTorusIntersection`，函数体改持 `DirectPolynomialRoots::new5(..)` 以分清 `mdpr.IsDone()`（`cxx:89`）与 `NbSolutions() == 0`，消费方 `crates/occt-topo/src/intcurvesurface/solvers.rs` 按 `IntCurveSurface_InterUtils.pxx:1290-1294` 回落",
    state: "✅ 已收口（2026-10-16，包 P032，commit `55b0a6d6`）。<br>(1) **`done` 载荷**：新增 `LineTorusIntersection { done, points }`（`line_torus.rs:826-843`），`is_done()` = `IntAna_IntLinTorus::IsDone()`（`IntAna_IntLinTorus.lxx:18-21`），`points()` = `NbPoints()`/`Value(i)` 的载荷（`lxx:23-44`，`!done` 时为空表）；`line_torus_intersect`（`line_torus.rs:851-924`）改返回该结构。<br>(2) **不再经 `quartic_roots` 适配器**：`occt_core::math_direct_poly_roots::quartic_roots`（`math_direct_poly_roots.rs:724-726`）把 `!IsDone()` 与 `NbSolutions() == 0` 都折成空 `Vec<f64>`，而 cxx 对两者给出相反的 `done` —— 于是函数体直接持有 `DirectPolynomialRoots::new5(a4, a3, a2, a1, a0)`，先按 `cxx:89` 的 `if (mdpr.IsDone())` 分流：未完成即 `done = false` 早退（`cxx:126-131`，原文 `nbpt = 0; done = false;`），完成则走 `NbSolutions()`/`Value(i)` 循环（`cxx:90-91`）。<br>(3) **两支 `done = false` 都落地**：`cxx:116-120` 的 `if (n > 0 && nbsolvalid == 0 && aNbBadSol == n) { nbpt = 0; done = false; }` 现按 `n > 0 && out.is_empty() && nb_bad_sol == n` 实现（坏解计数在回代判据 `PSolT.SquareDistance(PSolL) > 1e-10` 处分增），即「有根但全被回代判为坏解」不再与「真无根」混为同一个空表；`cxx:121-125` 的完成态（含 `n == 0`）仍回 `done = true`。<br>(4) **消费方回落**：`solvers.rs:456-465` 收 `inter` 后先判 `is_done()`，未完成即 `return None`，即 `IntCurveSurface_InterUtils.pxx:1290-1294` 的 `if (!intlintorus.IsDone()) { return false; }`，其后果为 `IntCurveSurface_Inter.pxx:578-581` 清 `isAnaProcessed` 并把该对送多面体干涉（本仓 `perform_curve_surface` 的 `polyhedron_path` 臂）。<br>门禁：五 crate `--all-targets` **0 error**（12/30/105/82/413 warnings，与 T002 同数，无新增告警）；`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0；`step_obj_gates` 4/4（直跑测试二进制，`.target-gate/t162/gates-raw.txt` 与 P031 的 `.target-gate/t161/gates-raw.txt` 逐字节相同）；`export_data_obj` 15/15 逐位一致（`.target-gate/t162/export-obj.txt` 与 `t156` 逐字节相同）；仅更新 `intana/tests.rs` 两条既有断言以适配新返回类型，未新增测试。<br>**诚实边界**：四道门禁逐位一致说明两处 `done = false` 在 15 个导出模型与 23 个对拍模型上都未触发 —— 只证明无回归，不证明这两支本身。<br>**下一片候选**：audit `occt-geom.md` §9 的 (b) 项 —— `out.sort_by(|a, b| a.x()...)`（`line_torus.rs:922`）在 cxx 无对应（`cxx:104-108` 按根序号写入 `thePoint[nbsolvalid]`，不重排）。",
  },
  {
    pkg: "P033",
    task: "T-163",
    commit: "a4eeac65",
    scope: "删掉 `line_torus_intersect` 里 cxx 不存在的点重排（audit `occt-geom.md` §9 的 (b) 项）：`out.sort_by(|a, b| a.x()...)`（`crates/occt-geom/src/intana/line_torus.rs:922`）在 `IntAna_IntLinTorus.cxx:104-108` 无对应（cxx 按根序号写入 `thePoint[nbsolvalid]`，函数体内无任何排序），并随删除撤掉 `crates/occt-geom/src/intana/mod.rs:12` 已无使用点的 `Ordering` 再导出",
    state: "✅ 已收口（2026-10-16，包 P033，commit `a4eeac65`）。<br>(1) **无对应当处**：`IntAna_IntLinTorus.cxx:104-108` 把每个通过回代判据的根按序号写进 `thePoint[nbsolvalid]`（`theFi`/`theTheta`/`theParam` 同步），函数体内没有任何排序；本仓 `line_torus.rs:922` 的 `out.sort_by(|a, b| a.x().partial_cmp(&b.x()).unwrap_or(Ordering::Equal))` 随 `09e227c2` 引入，属门禁第 2 条禁止的「OCCT 里不存在的规则」，已删除；原处改为注记 cxx 的写入口径。<br>(2) **消费方无依赖**（audit §9 的判断经复核）：两处调用点都不读该顺序 —— `crates/occt-topo/src/intcurvesurface/solvers.rs:210` 的 `perform_curve_surface` 统一按曲线参数重排（`points.sort_by(|a, b| a.param...)`），`crates/occt-geom/src/intana/tests.rs:232` 与 `crates/occt-topo/src/intcurvesurface/tests.rs:232` 各自再排序，故删除不改任何既有断言。<br>(3) **顺带撤掉无用再导出**：删除该排序后 `crates/occt-geom/src/intana/mod.rs:12` 的 `pub(crate) use std::cmp::Ordering;` 在 `occt-geom` 全 crate 内再无使用点（`intana` 其余模块与 `intana_torus.rs` 都用全限定 `std::cmp::Ordering`），一并删除以免留下 `unused_imports`；五 crate 告警数 12/30/105/82/413 与 T002 同数，无新增。<br>门禁：五 crate `--all-targets` **0 error**（12/30/105/82/413 warnings）；`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0；`step_obj_gates` 4/4（直跑测试二进制 124.16s，`.target-gate/t163/gates-raw.txt`，23 面积比与 23 条 `ours v=/f=` 行与 P032 的 `.target-gate/t162/gates-raw.txt` 逐字一致，仅并发 println 交错顺序与耗时不同）；`export_data_obj` 15/15 且 15 行模型读数与 `t162` 逐字一致（`.target-gate/t163/export-obj.txt`，整档仅 `Finished ... in` 耗时行不同）；未新增测试、未改既有断言。<br>**诚实边界**：四道门禁一致只证明该重排在当前语料上无消费方依赖（与 audit §9 的事前判断一致），不构成「排序会算错」的证据 —— 删除的依据是 cxx 无对应，不是结果差异。<br>**下一片候选**：`occt-geom.md` §9 的 (a)/(b) 两项至此全部闭环，代码侧登记只剩那条无基线的 STEP 写侧 harness（T999 裁定本 tranche 范围外），T999 复跑的前置条件回到 clean。",
  },
  {
    pkg: "P034",
    task: "T-164",
    commit: "c272c08a",
    scope: "3D 点-曲线**局部**极值 `Extrema_LocateExtPC` / `Extrema_GLocateExtPC`（`Extrema_GLocateExtPC.hxx:40-399`）在 `crates/occt-geom/src/extrema_pc/` 整体缺件：补 `curve_tool.rs`（`Extrema_CurveTool`）、`gfunc_ext_pc.rs`（`Extrema_GFuncExtPC`）、`gen_locate_ext_pc.rs`（`Extrema_GenLocateExtPC`）、`locate_ext_pc.rs`（`Extrema_GLocateExtPC` 的 3D 单态化），并把三处自创夹取 Newton 调用点改接：`crates/occt-topo/src/int_tools_vertex_line.rs`（`IntTools_Context.cxx:842-871`）、`crates/occt-topo/src/shhealing/wire_fix/same_parameter.rs`（`BRepLib_ValidateEdge.cxx:168-227`）、`crates/occt-topo/src/brep_lib_continuity.rs`（`BRepLib.cxx:2341-2355`）",
    state: "✅ 已收口（2026-10-16，包 P034，commit `c272c08a`）。<br>(1) **四个忠实件落地**：`curve_tool.rs` = `Extrema_CurveTool`（`hxx:38-142`：`GetType`/`FirstParameter`/`LastParameter`/`Resolution`/`NbIntervals`/`Intervals`/`D0`/`D1`/`D2`/`DN`）；`get_type` 按 `GeomAdaptor_Curve::load`（`GeomAdaptor_Curve.cxx:239-312`）的 `DynamicType` 顺序判类，**不用几何判据**（几何上笔直的 B 样条在 OCCT 里仍是 `GeomAbs_BSplineCurve`）。`gfunc_ext_pc.rs` = `Extrema_GFuncExtPC`（`hxx:146-415`）：`theF = PPc.Dot(D1c) / Ndu`（`hxx:251`）、`theDF = Ndu + PPc.Dot(D2c)/Ndu - theF*(D1c.Dot(D2c))/(Ndu*Ndu)`（`hxx:346`）、`IsMin` 即 `myD1f > 0.0`（`hxx:371-376`）、`|D1c|` 过小时的 `PPc` 三点/`DN` 兜底，`myMaxDerivOrder` 的 `GetType` 开关（`hxx:87-99`、`hxx:115-128`：Bezier/BSpline/Offset/Other -> `MaxOrder` + `SearchOfTolerance`，其余 -> `0` + `MinTol`）。`gen_locate_ext_pc.rs` = `Extrema_GenLocateExtPC`（`hxx:94-159`）：`Initialize` 存 `myumin/myusup/mytolU`（`hxx:94-105`），`Perform` 走已移植的 `math_FunctionRoot`（`hxx:111`），`myDone = S.IsDone()`（`hxx:112`）后按 `hxx:116-126` 用 `|ff| >= 1.e-07` 复核并清 `myDone`。<br>(2) **`Extrema_GLocateExtPC` 三臂**（`hxx:125-307`）：`Initialize`（`hxx:102-122`）按 `GetType` 分派（`tolu = TheCurveTool::Resolution(theC, Precision::Confusion())`）；`OtherCurve`/`OffsetCurve`/`BSplineCurve` 臂（`hxx:132-272`）先夹 `U0` 进窗口（`hxx:142-150`）、在 `NbIntervals(GeomAbs_C2)` 里找含 `U0` 的区间（`hxx:152-162`，按 `hxx:164-165` 的 `IFV 16.06.00` 注释把 `inter` 回退一格）、在该区间 `Initialize`+`Perform`（`hxx:167-176`），未命中则按 `hxx:183-270` 以 `k` 递增向两侧邻区间走 —— 入/出口切矢点积 `s1sup*s2sup <= RealEpsilon()`（`hxx:206`，不是 `<= 0`）或 `s1inf*s2inf <= RealEpsilon()`（`hxx:243`）即记为该区间端点极值（`myismin = (s1sup <= 0.0)` / `(s1inf <= 0.0)`，`hxx:211`/`:248`），否则对该区间再做一次 `Perform(theP, (myintuinf+myintusup)*0.5)`（`hxx:217-219`/`:252-254`），两向都要求 `Precision::IsInfinite` 时 `break`；`BezierCurve` 单列一臂（`hxx:275-281`）只做 `Perform(theP, theU0)`（既不夹 `U0` 也不走区间）；初等 `default` 臂（`hxx:282-303`）先判 `myExtremPC.IsDone()`（`hxx:285`）再按 `|u - U0|` 取最近解、`numberext == 0` 清 `myDone`（`hxx:299-300`）。<br>(3) **全局臂补 `IsDone` 载荷（承重）**：`general_extrema_pc.rs` 新增 `extrema_ext_pc_range_done(..) -> (Vec<ExtPcSolution>, bool)`，`extrema_ext_pc_range` 改为它的 `.0` —— `Extrema_GGExtPC` 自由臂的 `mydone = IntExtIsDone`（`Extrema_GGExtPC.hxx:468`），解表可在任何 `IntervalPerform` 成功前就非空（`hxx:453-455` 追加区间端点），而初等臂必须先判 `IsDone()`（`hxx:285`）。<br>(4) **三处调用点改接、自创件全删**：`int_tools_vertex_line.rs` 的 `LOCATE_ITERS`/`newton_step`/`clamp_param`/本地 `locate_ext_pc` 删除，改 `LocateExtPC::new_range(p, curve, t_seed, first, last, EXTREMA_EPS)`；`wire_fix/same_parameter.rs` 的 12 步 Newton `locate_ext_pc_sq` 改为「持久 `LocateExtPC` + `perform` + `!is_done()` 返 `None` + `square_distance()`」，每曲线在循环外 `Initialize` 一次（`BRepLib_ValidateEdge.cxx:168-177` 口径）；`brep_lib_continuity.rs` 的 `locate_on_cos`/`cos_d1`/`cos_d2`/`extrema_f`/`extrema_df` 与 `has_analytic_d2` 前置条件删除，改为 `CurveOnSurface::new(c2, s2_world, f, l)` + `LocateExtPC::new_range(..)`，只在 `is_done() && is_min()` 时按 `BRepLib.cxx:2347-2355` 重估第二面导数。`extrema_cc/mod.rs` 增 `pub use p_on_curv::POnCurv;`（共享类，两个模板都实例化它）。<br>门禁：五 crate `--all-targets` **0 error**（12/30/105/82/413 warnings，与 T002 同数，无新增告警）；`occt-geom --lib` 143/0、`occt-geom2d --lib` 72/0、`occt-core --lib` 283/0、`occt-topo --lib` 1251/0；`step_obj_gates` 4/4（直跑测试二进制 275.75s，`.target-gate/t164/gates-raw.txt`，23 面积比与 23 条 `ours v=/f=` 行与 P033 逐字一致）；`export_data_obj` 15/15 且 15 行模型读数与 `t163` 逐字一致；未新增测试、未改既有断言。<br>**诚实边界**：四道门禁与 P033 逐字一致只证明无回归，不证明自由曲线臂新增的邻区间行走与 `|F| >= 1.e-07` 清 `IsDone` 两支被覆盖（旧 Newton 与它们在这些语料上给出相同读数）。<br>**登记（未改）**：`CurveOnSurface` 缺 `Adaptor3d_CurveOnSurface::Resolution` 覆写（`Adaptor3d_CurveOnSurface.cxx:1364-1370`，覆写落在共享离散适配器上）与 `CurveOnSurface` 不暴露 `gp_line`（`EvalKPart` 的 `GeomAbs_Line` case，`cxx:1565`/`:1687`/`:1701`），两处均在源码就地标 UNPORTED。<br>**下一片候选**：P035 —— 从 `specs/_audit/occt-geom.md` / `specs/_audit/topo-rest.md` / canvas SEMANTIC 里找一条未复核的登记项，或判定 `CurveOnSurface::Resolution` 覆写是否属本轮可改。",
  },
  {
    pkg: "T999",
    task: "Judge",
    commit: "7b164e34",
    scope: "收口审计：T001 清单逐行对账（每个登记行要么有包可指、要么有记录在案的结论），四条门禁在冻结 HEAD 上重跑（五 crate `cargo check --all-targets`、`occt-topo --lib`、`step_obj_gates`、`export_data_obj`）",
    state: "✅ 结论 `full_outcome_complete: true`（`docs/goals/drift-and-gap-autofix/state.yaml` T999 收据）：B1..B8、B6a-h、B7-writer 两族、B8-intana、`surface_fit`、`B7-mesh_to_brep`、`ReparamCurve2d` 全部有包可指；冻结 HEAD `fbd59532` 上 `--all-targets` 0 error、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（118.02s）、`export_data_obj` 15/15 逐位一致；工作树在收口提交后为 clean。收据同时记明两条诚实边界：本 tranche 无一条门禁读数来自新写测试（四道基线都是 T002 冻结读数逐包重跑），唯一改动过的三条基线断言是以 `occt_probe --ssx` 实测值重定（非新测试）。<br>**该收据随后被 P027 与 P028 两次推翻（2026-10-16）**：其「B8-intana 全部有包可指」只覆盖锥/锥与柱族臂 —— `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Pln)` 与 `IntAna_Int3Pln::Perform` 两臂从未与 cxx 对拍，P027 在其中查出三处自创 + 一段未移植（最小范数原点解 + `1e-20` 叉积门、`cxx:449-508` 小角重解缺失、`solve3` 的 `1e-14` 手写消元）并收口；P028 又在该类的球臂查出两处（`plane_sphere` 的 `|radius| * DBL_EPSILON` 冒充 `Epsilon(radius)`，`cylinder_sphere` / `sphere_cone` 的 `1e-9` 轴距启发式，含 `occt-topo` 两份副本）并收口。故本 tranche 以 `active_task: T005` / `goal.status: active` 重开；教训：`IntAna_QuadQuadGeo` 这类多臂类的清单行必须**逐臂**核对，不能按文件记账 —— P028 的普查因此把 15 个 `Perform` 重载逐条点名（10 已移植 / 5 Torus 未移植），**而该普查本身被 P029 第三次推翻**：`intana/mod.rs` 用 `#[path = \"../intana_torus.rs\"] mod torus;` 把环面五臂放在目录外，按目录列举必然漏掉它们，正确读数是 15 个重载全部在仓、五臂各有生产调用方；P029 还在 `quadric_quadric_sphere_torus` / `quadric_quadric_torus_torus` 删掉 cxx 不存在的 `a_dist <= f64::EPSILON` 自创同心分支（`cxx:2496-2562` / `2588-2666`）。三轮下来该类的记账已被推翻三次（按文件 → 按目录 → 按 `mod`），T999 必须等剩余 5 个从未逐行对拍的重载走完再重跑。<br>**P030 已于 2026-10-16 走完这 5 个重载（commit `0a185eb5`，T-160）**：`cxx:543-745` Pln/Cylinder 查出并修复两处缺口（缺 `H` 参数使 `cxx:568-598` 近平行细化臂整段省略；并行判定用内联 `|axis . normal| < Tolang` 冒充 `IntAna_IntConicQuad(axec, P, tolang, toltang, H).IsParallel()`，消费方 `int_p_cy` 同步按 `IntPatch_ImpImpIntersection.cxx:2579-2587` 传面高），`cxx:752-972` Pln/Cone 的尾部 `cxx:922-951` 极端圆锥曲线回绝段按门禁第 2 条标 `UNPORTED`（本仓 `QuadricIntersection` 无「未完成」态），`cxx:1050-1293` 柱/柱、`cxx:1324-1342` 柱/锥、`cxx:1433-1886` 锥/锥 逐行核对为忠实。**至此该类 15 个重载全部逐行对拍且全部在仓，剩余 1 条代码侧登记项（declined 态）列为 P031**；T999 复跑的前置条件只剩 P031，另有一条无基线的 STEP 写侧 harness 保持登记。<br>**P031 已于 2026-10-16 清掉那条前置条件（commit `9bef99ad`，T-161）**：`QuadricIntersection` 新增 declined 态 `NotDone`（= cxx 的 `done = false`），`quadric_quadric_plane_cone` 按 `cxx:922-951` 翻译 `EllipseLimit = 1.0E+9` / `HyperbolaLimit = 2.0E+6` 的回绝（`ELLIPSE_LIMIT` / `HYPERBOLA_LIMIT` 两个具名常量），三条消费路径同步接上（`from_quadric_intersection` 的 `PairOutcome::Fail` = `IntPatch_ImpImpIntersection.cxx:3480-3483`；`line_torus.rs` 两处嵌套约化按 `cxx:1556` / `:1835` 保持 `None` = `IntAna_Empty`；两条自设解析快速路径归入「无曲线」并回落 `self.general`）。门禁：五 crate 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（`.target-gate/t161/gates-raw.txt`，23 面积比与 v/f 行与 P030 逐字一致）、`export_data_obj` 15/15 逐位一致。**`IntAna_QuadQuadGeo` 至此除那条无基线的 STEP 写侧 harness 外已无代码侧登记项，T999 复跑的前置条件已清空**；随后由 P032 接手同族的下一条（`IntAna_IntLinTorus` 的 done 语义，见 `occt-geom.md` §9）。<br>**P032 已于 2026-10-16 收口同族那条（commit `55b0a6d6`，T-162）**：`line_torus_intersect` 改返回 `LineTorusIntersection { done, points }`（`line_torus.rs:826-843`），函数体改持 `DirectPolynomialRoots::new5(..)` 以分清 `mdpr.IsDone()`（`cxx:89`）与 `NbSolutions() == 0`，`cxx:116-120`（有根但全为坏解）与 `cxx:126-131`（求解器未完成）两支都置 `done = false`；消费方 `solvers.rs:456-465` 在 `!is_done()` 时 `return None` = `IntCurveSurface_InterUtils.pxx:1290-1294` 的 `return false`。门禁：五 crate 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15，两者均与 P031 逐字节相同。**代码侧登记项现只剩本条无基线的 STEP 写侧 harness 与 `occt-geom.md` §9 的排序余项（(b) 项）两项**。<br>**P033 已于 2026-10-16 清掉其中的排序余项（commit `a4eeac65`，T-163）**：`line_torus_intersect` 里 `out.sort_by(|a, b| a.x()...)`（`line_torus.rs:922`）在 `IntAna_IntLinTorus.cxx:104-108` 无对应（cxx 按根序号写入 `thePoint[nbsolvalid]`，函数体内不排序），随 `09e227c2` 引入，已删除；`intana/mod.rs:12` 的 `Ordering` 再导出失去唯一使用点、一并撤掉；消费方 `solvers.rs:210` 按曲线参数重排、两处测试各自再排序，故无依赖语义。门禁：五 crate 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15，两者均与 P032 逐字一致（仅耗时行不同）。**`occt-geom.md` §9 的 (a)/(b) 两项至此全部闭环，代码侧登记只剩那条无基线的 STEP 写侧 harness（T999 裁定范围外），T999 复跑的前置条件回到 clean**。<br>**P034 已于 2026-10-16 第八次推翻本收据（commit `c272c08a`，T-164），而这次落在本清单从未提过的类上**：3D 点-曲线**局部**极值族在仓内整体缺件 —— `crates/occt-geom/src/extrema_pc/` 只有全局臂 `Extrema_GGExtPC`，没有 `Extrema_CurveTool` / `Extrema_GFuncExtPC` / `Extrema_GenLocateExtPC` / `Extrema_GLocateExtPC`，于是三个需要「已知参数附近求点-曲线最近极值」的调用点各自写了夹取 Newton：`int_tools_vertex_line.rs`（12 次迭代 + `\\|C'\\|^2 <= EXTREMA_EPS^2` 门 + 位移判据）对 `IntTools_Context.cxx:842-871`、`wire_fix/same_parameter.rs`（12 次迭代 + `1e-20` 速度门 + `1e-10` 位移判据）对 `BRepLib_ValidateEdge.cxx:168-227`、`brep_lib_continuity.rs`（100 次迭代 + `has_analytic_d2` 前置，使不满足的边整条静默跳过）对 `BRepLib.cxx:2341-2355`。四类已全部移植、三处自创循环已删除、`extrema_ext_pc_range_done` 导出初等臂所需的 `IsDone`。门禁：五 crate 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、`occt-geom2d --lib` 72/0、`occt-core --lib` 283/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15，两者均与 P033 逐字一致。**T999 复跑再次被推迟**：本收据的清单行必须按「类族」逐族重核（`Extrema_*` 局部极值族此前从未出现），而不是只核 `IntAna_*`。",
  },
  {
    pkg: "登记项",
    task: "—",
    commit: "—",
    scope: "缺失的 `data/occ/*.st*` 驱动式 STEP 写侧 harness：`data/occ/` 只有读侧 `step_obj_gates` 口径，写侧（`StepToGeom_Make*` / `GeomToStep_Make*`）没有可对拍的既有基线",
    state: "⬜ 保持登记（T999 裁定「本 tranche 范围外」）：无既有基线可对拍且本 tranche 不允许新写测试，故不改代码；它与 P010/P011/P012 的写侧缺陷不同 —— 那三条落在读侧共享管线上，能由既有门禁判定，本项不能。推进需先裁定「无基线不可对齐」还是另找既有写侧口径。代码侧登记项在 P026 收口时曾清零，P027 复核后登记两条、**P028 已全部收口**；P028 同时给出本类的新登记（见次行）。",
  },
  {
    pkg: "登记项",
    task: "—",
    commit: "—",
    scope: "audit `occt-geom.md` §9 的登记项（`IntAna_IntLinTorus`）：(a) 已由 P032 闭环、(b) 点重排已由 P033 闭环",
    state: "✅ 已闭环（2026-10-16，包 P033，commit `a4eeac65`）。(a) 项**已由 P032 闭环** —— `line_torus_intersect`（`line_torus.rs:851-924`）改返回 `LineTorusIntersection { done, points }`（`line_torus.rs:826-843`），直接持有 `DirectPolynomialRoots::new5(..)` 并按 `IntAna_IntLinTorus.cxx:89` 的 `mdpr.IsDone()` 分流，`cxx:116-120`（`n > 0 && nbsolvalid == 0 && aNbBadSol == n`）与 `cxx:126-131`（求解器未完成）两支都置 `done = false`，消费方 `intcurvesurface/solvers.rs:456-465` 在 `!is_done()` 时 `return None`（= `IntCurveSurface_InterUtils.pxx:1290-1294` 的 `return false`）。(b) 项**已由 P033 闭环** —— `line_torus.rs:922` 的 `out.sort_by(|a, b| a.x()...)` 在 `IntAna_IntLinTorus.cxx:104-108` 无对应（cxx 按根序号写入 `thePoint[nbsolvalid]`，函数体内不排序），随 `09e227c2` 引入，属门禁第 2 条禁止的「OCCT 里不存在的规则」，已删除（该处的原注释换成 cxx 写入口径的注记）；删除后 `intana/mod.rs:12` 的 `Ordering` 再导出失去唯一使用点，一并撤掉；两处调用点都不读该顺序（`intcurvesurface/solvers.rs:210` 的 `perform_curve_surface` 按曲线参数重排、`intana/tests.rs:232` 与 `intcurvesurface/tests.rs:232` 各自再排序），故删除不改任何既有断言。门禁：五 crate 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15，两者均与 P032 逐字一致。代码侧登记项沿革：P026 收口时清零 → P027 后两条 → P028 全部收口 → P029 后五条 → P030 后一条 → P031 清掉 `IntAna_QuadQuadGeo` 的 declined 态 → P032 清掉 `IntAna_IntLinTorus` 的 done 语义 → **P033 清掉 §9 的排序余项**，现只剩本节上条（无基线的 STEP 写侧 harness）。",
  },
];

const pct = (a: number, b: number) => `${(100 * (a / b)).toFixed(1)}%`;

export default function OcctScopeDrift() {
  return (
    <Stack gap={16}>
      <Stack gap={6}>
        <H1>OCCT→Rust 移植范围与偏移审查</H1>
        <Text tone="secondary">
          对照 OCC 8.0.0 源码树 <Text weight="semibold">D:\source\OCCT-src\src</Text> 与工作树 5 个 crate 实测。
          所有读数取 <Text weight="semibold">{MEASURE.asOf}</Text>；门禁基线引用{" "}
          <Text weight="semibold">docs/goals/drift-and-gap-autofix/state.yaml</Text>（板真相）与{" "}
          <Text weight="semibold">specs/_audit/_index.md</Text>（A0–A27 各行「处置」，逐包同步；最新的 P031 记录在{" "}
          <Text weight="semibold">specs/_audit/occt-geom.md</Text> §8）。
        </Text>
      </Stack>

      <Grid columns={4} gap={12}>
        <Stat value={pct(MEASURE.rust.lines, MEASURE.occt.lines)} label="全库行数比（CAD 全量）" />
        <Stat value={pct(MEASURE.rust.lines, MEASURE.cadKernel)} label="CAD 内核行数比" tone="info" />
        <Stat value={MEASURE.rust.lines.toLocaleString("en-US")} label={`Rust 行数 · ${MEASURE.rust.files.toLocaleString("en-US")} 文件`} />
        <Stat value={MEASURE.occt.lines.toLocaleString("en-US")} label="OCCT 行数 · 13,803 文件" tone="info" />
      </Grid>

      <Divider />

      <H2>一、范围</H2>
      <Text tone="secondary">
        模块组按 OCCT 源码分组统计；「对应 crate」是当前承载该层职责的 Rust 模块，「深度」是已核实的移植饱和度。
      </Text>

      <Table
        headers={["OCCT 模块组", "文件", "行数", "Rust 对应", "深度"]}
        rows={GROUPS.map((g) => [
          g.group,
          g.files.toLocaleString("en-US"),
          g.lines.toLocaleString("en-US"),
          g.rust,
          g.depth,
        ])}
        columnAlign={["left", "right", "right", "left", "left"]}
        striped
      />

      <H3>OCCT 各组行数（单位：源码行）</H3>
      <BarChart
        categories={GROUPS.map((g) => g.group.replace("Application", "App"))}
        series={[{ name: "OCCT 行数", data: GROUPS.map((g) => g.lines) }]}
        height={220}
        showValues
      />
      <Text size="small" tone="tertiary">
        来源：D:\source\OCCT-src\src（8 个模块组，.hxx/.cxx/.lxx/.gxx）· 2026-10-16 复测（与 10-13 / 10-11 逐位相同）· 单位：源码行。
      </Text>

      <H3>Rust 侧分布</H3>
      <Table
        headers={["crate", "文件", "行数"]}
        rows={CRATES.map((c) => [c.crate, String(c.files), c.lines.toLocaleString("en-US")])}
        columnAlign={["left", "right", "right"]}
      />
      <Text size="small" tone="tertiary">
        较 2026-10-13（1,060 文件 / 335,743 行）：<Text weight="semibold">+47 文件 / +2,297 行</Text>；较 2026-10-11
        （1,080 文件 / 337,891 行）为 <Text weight="semibold">+27 文件 / +149 行</Text>。增量来自 P024..P031
        等包的忠实件新增（`brep_lib_find_surface.rs`、`brep_lib_make_face.rs`、`extrema_ss/` 与 `extrema_cc/`
        的 `*_func_*` / `gen_locate_*` / `p_on_curv`、`extrema2d/` 的 `gfunc_ext_*` / `gen_locate_*`、
        `bspl/knot_analysis.rs`、P030 的 `H` 臂与其 `IntAna_IntConicQuad` 并行判定、P031 的 `NotDone` 态与平面/锥 tail），
        减量来自删除面（`mesh_to_brep` 的 `1e-9` 归并格、`SameRange` 转换尾副本、
        `surface_fit` 调用点、四份自创多项式解算器、`brepmesh.rs` −602 行与 `wireframe.rs` 平面机器 −288 行）。
      </Text>

      <Callout tone="info" title="范围口径的要点">
        实际在办范围 ≠ 全库移植：看板目标已收敛为「STEP→OBJ 几何/网格管线对齐 OCCT」，
        即在 CAD 内核里只对 STEP 读/写 → Geom/BRep → BOPAlgo/BRepMesh/ShapeFix 这条共享链做逐行忠实化；
        Visualization / Draw / OCAF / Deprecated 是有意 lite 替代，不计入门禁。
      </Callout>

      <Divider />

      <H2>二、偏移</H2>

      <H3>2.1 本轮偏移收敛（已跑现基线核实）</H3>
      <Table
        headers={["已修正项", "编号", "位置（当前工作树）", "核实读数"]}
        rows={CLOSED.map((c) => [c.item, c.id, c.where, c.evidence])}
        columnAlign={["left", "left", "left", "left"]}
        rowTone={CLOSED.map(() => "success" as const)}
        striped
      />

      <H3>2.2 语义偏移残余（相对 OCCT 的忠实度残余）</H3>
      <Row gap={8} align="center" wrap>
        <Pill size="md">自创已收敛</Pill>
        <Text tone="secondary">
          全仓共 <Text weight="semibold">{MEASURE.unported} 处 UNPORTED 标记</Text>；下表为仍可核实到的残余项。
        </Text>
      </Row>
      <Table
        headers={["残余项", "审计编号", "位置（当前工作树）", "状态"]}
        rows={SEMANTIC.map((s) => [s.item, s.id, s.where, s.state])}
        columnAlign={["left", "left", "left", "left"]}
        rowTone={SEMANTIC.map((s) =>
          s.state.startsWith("UNPORTED") || s.state.startsWith("未移植")
            ? ("danger" as const)
            : s.state.startsWith("活调用点")
              ? ("info" as const)
              : ("warning" as const),
        )}
        striped
      />

      <H3>2.3 逐模型对拍（`data/occ/*`，2026-10-16 `step_obj_gates --nocapture` 实测）</H3>
      <Table
        headers={["模型", "OCCT v / f", "端口 v / f", "Δf", "面积比"]}
        rows={PARITY.map((p) => [
          p.model,
          `${p.occV} / ${p.occF}`,
          `${p.ourV} / ${p.ourF}`,
          String(p.ourF - p.occF > 0 ? `+${p.ourF - p.occF}` : `${p.ourF - p.occF}`),
          p.area,
        ])}
        columnAlign={["left", "right", "right", "right", "right"]}
        rowTone={PARITY.map((p) => (p.ourF === p.occF ? ("success" as const) : Math.abs(p.ourF - p.occF) > 100 ? ("danger" as const) : ("info" as const)))}
        striped
      />
      <Text size="small" tone="tertiary">
        `data/*.step`（15 个）中 12 个逐位一致；`Shape-1` 差 −8 面、`Shape-2` 差 −34 面、`screw` 差 +2 面。八个 `data/occ/*` 模型里 `motoc` 逐位一致，其余为双向小量：T0M −349 v / −469 f、top −17 v / −30 f、ATU01038 −17 v / −32 f、a3n00 −7 v / −24 f、acs10 −16 v / −20 f、bottom −6 v / −7 f、TDB +42 v / +17 f。**2026-10-16 复测订正（`.target-gate/t157/gates-raw.txt`）**：23 个面积比中 20 个为 `1.0000`，`Shape-2` `0.9999`、`T0M` `1.0017`、`TDB` `1.0015`；各模型的 `area_tol` 是登记值、本卡未改过，故无门禁放宽。T0M/TDB 的网格数相对 2026-10-13 读数（T0M f 67,424 / 面积比 1.0000；TDB f 79,926 / 面积比 1.0000）发生了位移，属 10-15～10-16 各包在共享求交/极值路径上的忠实化后果，已登记为旁支（不调参对齐）；剩下面级差集中在 2.2 首行登记的 Delaunay 插点密度分支。
      </Text>

      <H3>2.4 状态 / 文档偏移（2026-10-16 复核：canvas 侧全部收口；`topo-rest.md` 四条 lite 标注留档）</H3>
      <Table
        headers={["偏差", "文档现状", "实际"]}
        rows={DOCDRIFT.map((d) => [d.what, d.detail, d.gap])}
        columnAlign={["left", "left", "left"]}
        rowTone={DOCDRIFT.map((d) => d.tone)}
        striped
      />

      <H3>2.5 门禁现状</H3>
      <Table
        headers={["门禁", "说明", "基线", "现测", "状态"]}
        rows={GATES.map((g) => [
          g.id,
          g.name,
          g.base,
          g.now,
          g.tone === "success" ? "green" : "red",
        ])}
        columnAlign={["left", "left", "right", "right", "left"]}
        rowTone={GATES.map((g) => g.tone)}
      />
      <Text size="small" tone="tertiary">
        `G-check` / `G-lib` / `G-gates` / `G-export` 读数首次测量于 2026-10-15（覆盖到 T-147 已提交工作树 `e3a4335e`），末次为其在 2026-10-16 收口审计所对的冻结 HEAD `fbd59532`：五 crate `--all-targets` 0 error、`occt-geom2d --lib` 72 / 0、`occt-geom --lib` 143 / 0、`occt-topo --lib` 1251 / 0、`step_obj_gates` 4 / 4（118.02s；23 个面积比全部落在登记的 `area_tol` 内：20 个 `1.0000`、`Shape-2` `0.9999`、`T0M` `1.0017`、`TDB` `1.0015`）、`export_data_obj` 15 / 15 且逐模型 v/f 与包围盒与 T002 冻结表逐位一致。中间每个包都在自己的提交上重跑同样四条门禁，四项**判词**一直不变：P018 `8cad6a01`（gates 4/4，150.7s）、P019 `616ea0c1`（181.8s）、P020 `2452a51c`（205.1s）、P021 `11b12e0e`（114.1s）、P022 `3b5728b3`（114.7s）、P023 `38a13eae`（115.6s）、P024 `549e5761`、P025 `8d77b13a`（115.7s）、P026 `fbd59532`（118.1s）、P027 `76893533`（251.49s）、P028 `1d595f8d`、P029 `6330c6d7`（151.50s）、P030 `0a185eb5`（131.26s）、P031 `9bef99ad`（251.50s）、P032 `55b0a6d6`（直跑测试二进制，日志与 P031 逐字节相同）；`--lib` 恒为 `occt-topo` 1251 / 0、`export_data_obj` 恒为 15 / 15 逐位一致。**P030・P031・P032 的 `step_obj_gates` 读数除耗时外逐字相同**（`.target-gate/t160`、`.target-gate/t161` 与 `.target-gate/t162` 的 23 条面积比与 23 条 `ours v=/f=` 行逐字一致，t162 的 diff 判词为「Files are identical」）。**2026-10-16 独立复测**（`.target-gate/t157/`，货真基于 `fbd59532` 之后的工作树）：`step_obj_gates` 4 / 4（cargo 177.06s / 直跑测试二进制 163.58s）、`export_data_obj -- data/occ data` 23 ok 且 15 个 `data/*.step` 逐位一致、`--lib` 1251 / 0；逐模型 v/f 与 `.target-gate/t156/gates-raw.txt` 逐位相同，仅面积比读数按上表订正。`G-fuse` / `G-phase` / `G-boss` / `G-iges` 为 T-130..T-133 轮读数，其后的包未触及这些路径。`G-core` 的 `occt-core` 读数由 285 变 283 是 T-143 删除两条自创单测所致（同 `G-lib` 的 1276 → 1251，`notes/T002-gate-baseline.md` 已留档）；`occt-math` 210、`occt-geom` 143、`occt-geom2d` 72 未变。
      </Text>

      <H3>2.6 目标看板跟踪（`docs/goals/drift-and-gap-autofix/state.yaml`）</H3>
      <Text size="small" tone="tertiary">
        `docs/goals/drift-and-gap-autofix/state.yaml` 于 2026-10-16 随 T999 收口：`goal.status: done`、`active_task: T999`、T999 收据 `full_outcome_complete: true`（唯一登记项为缺失的 STEP 写侧 harness，裁定为本 tranche 范围外）。**同日的 P027 复核推翻了该收据对 `B8-intana` 的记账**：`IntAna_QuadQuadGeo` 的平面/平面臂与 `IntAna_Int3Pln` 臂从未与 cxx 对拍，P027 在其间查出三处自创 + 一段未移植并收口，故 board 以 `active_task: T005` / `goal.status: active` 重开（T005 为 PM 滚动推进；每个包 = 一个已跑现基线核实并单独提交的切片，提交信息点名 `.cxx` 文件与行号范围）。
      </Text>
      <Table
        headers={["包", "任务", "提交", "切片范围", "状态"]}
        rows={CYCLES.map((c) => [c.pkg, c.task, c.commit, c.scope, c.state])}
        columnAlign={["left", "left", "left", "left", "left"]}
        rowTone={CYCLES.map((c) => (c.pkg === "登记项" ? ("warning" as const) : ("success" as const)))}
        striped
      />

      <Divider />

      <H2>三、结论与下一步</H2>
      <Row gap={16} wrap>
        <Stat value="13.9%" label="全库行数覆盖" />
        <Stat value="24.2%" label="CAD 内核行数覆盖" tone="info" />
        <Stat value="2 处" label="登记未收代码项（P034 新增两条 `CurveOnSurface` 覆写缺失：`Resolution`（`Adaptor3d_CurveOnSurface.cxx:1364-1370`）与 `gp_line`（`EvalKPart` 的 `GeomAbs_Line` case，`cxx:1565`/`:1687`/`:1701`），均已在源码就地标 UNPORTED；另 1 条无基线的 STEP 写侧 harness 保持留档）" tone="warning" />
        <Stat value="1 处" label="留档待办（`topo-rest.md` 9/10/15/17 的 live 标注）" tone="warning" />
      </Row>
      <Stack gap={6}>
        <Text>
          <Text weight="semibold">范围</Text>：全库 13.9%、CAD 内核约 24.2% 的行数覆盖；有效范围是 STEP→OBJ 共享管线，
          在此链上已做逐行对 `.cxx` 的深度移植（occt-topo 占 21.4 万行）。
        </Text>
        <Text>
          <Text weight="semibold">跟踪</Text>：2.6 的滚动看板把最近九个包（P025 `mesh_to_brep` 的 `1e-9` 归并格 / `Vertex` 容差 0.0 / 自创三角形法向换 `BRepBuilderAPI_MakeShapeOnMesh` 忠实件、P026 `GeomLib::SameRange` 转换尾两份副本合并为唯一 `trimmed_reparam_bspline`、P027 `IntAna_QuadQuadGeo` 平面/平面臂与 `IntAna_Int3Pln` 的忠实化、P028 球/平面 `Epsilon(radius)` 与球心-轴精确交点判定、P029 环面五臂逐行核对并订正 Torus 普查、P030 余下五臂逐行对拍：Pln/Cylinder 补 `H` 近平行细化臂与 `IntPCy` 面高、Pln/Cone 尾部 `cxx:922-951` 标 `UNPORTED`、柱族三臂核对为忠实、P031 把该 `UNPORTED` 落地成 `QuadricIntersection::NotDone` 的 declined 态并按 `cxx:922-951` 翻译回绝、P032 把 `IntAna_IntLinTorus` 的 `done` 语义补齐：`line_torus_intersect` 返回携带 `done` 的 `LineTorusIntersection`、函数体改持 `DirectPolynomialRoots::new5(..)` 以分清 `cxx:89` 的 `mdpr.IsDone()` 与零根、消费方按 `IntCurveSurface_InterUtils.pxx:1290-1294` 回落、P033 删掉 `line_torus_intersect` 里 cxx 不存在的点重排（`line_torus.rs:922` 的 `out.sort_by`，`cxx:104-108` 只按根序号写入 `thePoint[nbsolvalid]`）并撤掉随之失去使用点的 `Ordering` 再导出）各自落到任务号与提交号；B6 家族（B6a/B6b/B6c/B6d/B6e 2D+3D/B6f/B6g/B6h）已随 P008/P017/P018/P019/P020/P021/P023 全部收口，B7 写侧两族亦随 P010/P011/P022 收口，`surface_fit` 一项随 P024 收口（真实对应件为 `BRepLib_FindSurface`，非 `GeomLib::IsPlanarSurface`），B7-mesh_to_brep 一项随 P025 收口、`ReparamCurve2d` 的 `SameRange` 一项随 P026 收口；P026 的 T999 收据（`full_outcome_complete: true`）与那一条保持登记的 STEP 写侧 harness 曾收口本 tranche（`state.yaml` 的 `goal.status` 置 `done`），但同日 **P027 的复核推翻了对 `B8-intana` 的记账**：`IntAna_QuadQuadGeo` 的多臂性质导致平面/平面臂与 `IntAna_Int3Pln` 臂漏审，P027 在其间收口三处自创 + 一段未移植（最小范数原点解 + `1e-20` 叉积门、`cxx:449-508` 小角重解缺失、`solve3` 的 `1e-14` 手写消元、自创 `Intersect3Pln::Line`），board 随之重开为 `active_task: T005` / `goal.status: active`（登记项：STEP 写侧 harness + P027 登记的 `plane_sphere` `Epsilon(radius)` 与球/球 `0.01 * CONFUSION` 两臂）；P028 随后收口了 P027 登记的两臂（`plane_sphere` 改 `epsilon(radius)`、`cylinder_sphere` / `sphere_cone` 的自创 `1e-9` 轴距门换逐行 `axes_meet_at_point`，含 `occt-topo` 两份副本；球/球 `0.01 * CONFUSION` 经核对为忠实，无改动），并给出本类唯一剩下的结构未知：15 个 `Perform` 重载中 10 已移植 / 5 Torus 未移植，且已移植臂里仅 5 个做过逐行对拍 —— 下一片按逐臂口径把这两件事一次判定完。**P029 紧接着推翻了 P028 那份普查本身**：`intana/mod.rs:30-35` 的 `#[path = "../intana_torus.rs"] mod torus;` 让环面五臂位于 `crates/occt-geom/src/intana_torus.rs`（目录之外、同命名空间），五臂全部在仓且各有生产调用方（`intpatch_impimp_pairs.rs:117/155/162/169/175`、`intpatch.rs:576`），正确读数是 **15 个重载全部在仓**；P029 在同轮删掉 `quadric_quadric_sphere_torus` / `quadric_quadric_torus_torus` 里 cxx 不存在的 `a_dist &lt;= f64::EPSILON` 自创同心分支（cxx 在 `IntAna_Empty` 后无条件构造 `gp_Dir aVec12`，`cxx:2542` / `cxx:2647`），并把「按 `mod` 声明而非目录列举」写进 `notes/T001-inventory.md`。于是本类的记账与 `B8-intana` 这一行已被推翻三次：P027（按文件记账漏了平面臂）、P028（按目录记账漏了 `#[path]` 模块）、P029（订正普查并清掉环面自创）；`IntAna_QuadQuadGeo` 现已知 15/15 在仓，但其中仅 10 个做过逐行对拍，剩余 5 个（`cxx:543-745`、`:752-972`、`:1050-1293`、`:1324-1342`、`:1433-1886`）是下一片。**P030 随即把这一片走完**：`cxx:543-745`（Pln/Cylinder）查出两处缺口并修复 —— 签名缺 `H` 使 `cxx:568-598` 的近并行方向细化臂整段省略（旧注释自认 `H = 0` 跳过），并行判定用内联 `|axis . normal| &lt; Tolang` 冒充 `IntAna_IntConicQuad(axec, P, tolang, toltang, H).IsParallel()`（`cxx:601`，`H != 0` 时还有 `Len` 端点距离检查，两者不同义），消费方 `intpatch_impimp_pairs.rs::int_p_cy` 原先不传面高（等于恒按 `H = 0`），现按 `IntPatch_ImpImpIntersection.cxx:2579-2587` 取柱面 `VMax - VMin`；`cxx:752-972`（Pln/Cone）尾部 `cxx:922-951` 的 `EllipseLimit = 1.0E+9` / `HyperbolaLimit = 2.0E+6` 回绝段（命中置 `done = false`，消费方翻成 `return false` = 本仓 `PairOutcome::Fail`）因本仓 `QuadricIntersection` 无「未完成」态而按门禁第 2 条标 `UNPORTED`，理由写进 `quadric_quadric_plane_cone` 的函数文档；`cxx:1050-1293` 柱/柱、`:1324-1342` 柱/锥、`:1433-1886` 锥/锥 三臂逐行核对为忠实（P015 的旧笔记只是引文核对，本轮补齐线级 diff）。于是 `B8-intana` 这一行的记账已被推翻并修正四次：P027（按文件记账漏平面臂）、P028（按目录记账漏 `#[path]` 模块）、P029（订正普查并清掉环面自创）、P030（清掉最后五臂，查出 `H` 臂缺失与一处 `UNPORTED`）；该类 15 个重载至此全部逐行对拍且全部在仓。**P031 随即将那条唯一的 `UNPORTED` 落地**：`QuadricIntersection` 新增 declined 态 `NotDone`（= cxx 的 `done = false`），`quadric_quadric_plane_cone` 按 `cxx:922-951` 翻译 `EllipseLimit = 1.0E+9`（`cxx:926-934`）与 `HyperbolaLimit = 2.0E+6`（`cxx:935-951`）的回绝，三条消费路径同步接上（`from_quadric_intersection` 的 `PairOutcome::Fail` = `IntPatch_ImpImpIntersection.cxx:3480-3483`；`line_torus.rs` 两处嵌套约化的 `None` = `IntAna_Empty`，`cxx:1556` / `:1835`；两条自设解析快速路径归入「无曲线」并回落 `self.general`），门禁 4/4 且与 P030 逐字一致。**至此 `IntAna_QuadQuadGeo` 无代码侧登记项，代码侧登记转由同族的 `IntAna_IntLinTorus` done 语义接手（列为 P032）**，另有 1 条无基线的 STEP 写侧 harness 保持留档。**P032 随即把那条接手项落地**：`line_torus_intersect`（`line_torus.rs:851-924`）改返回携带 `done` / `points` 载荷的 `LineTorusIntersection`（`line_torus.rs:826-843`），即 `IntAna_IntLinTorus::IsDone()` 与 `NbPoints()/Value(i)` 载荷；函数体不再走把 `!IsDone()` 折成空根表的 `quartic_roots` 适配器，而直接持有 `DirectPolynomialRoots::new5(a4, a3, a2, a1, a0)`，先按 `cxx:89` 的 `mdpr.IsDone()` 分流（未完成即 `done = false`，`cxx:126-131`），再走 `NbSolutions()`/`Value(i)` 循环（`cxx:90-91`），并按 `cxx:116-120` 的「根数非零却全被回代判为坏解」（`nbsolvalid == 0 && aNbBadSol == n`）把该情形也置 `done = false`；消费方 `intcurvesurface/solvers.rs:456-465` 在 `!is_done()` 时 `return None`，即 `IntCurveSurface_InterUtils.pxx:1290-1294` 的 `return false` → `IntCurveSurface_Inter.pxx:578-581` 清 `isAnaProcessed` 转多面体干涉。门禁：五 crate 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4 与 `export_data_obj` 15/15 均与 P031 逐字节相同。代码侧登记项沿革再进一格：P031 清掉 `IntAna_QuadQuadGeo` 的 declined 态 → **P032 清掉 `IntAna_IntLinTorus` 的 done 语义**，现剩 `occt-geom.md` §9 的 (b) 排序余项（`line_torus.rs:922` 的 `out.sort_by` 在 cxx 无对应）与那条无基线的 STEP 写侧 harness。**P033 随即清掉其中的排序余项**：`line_torus_intersect` 里 `out.sort_by(|a, b| a.x()...)`（`line_torus.rs:922`）在 `IntAna_IntLinTorus.cxx:104-108` 无对应（cxx 按根序号写入 `thePoint[nbsolvalid]`，函数体内不排序），随 `09e227c2` 引入，已删除并换成 cxx 写入口径的注记；`intana/mod.rs:12` 的 `Ordering` 再导出失去唯一使用点、一并撤掉；消费方 `solvers.rs:210` 按曲线参数重排、`intana/tests.rs:232` 与 `intcurvesurface/tests.rs:232` 各自再排序，故无依赖语义。门禁：五 crate 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4 与 `export_data_obj` 15/15 均与 P032 逐字一致（`.target-gate/t163/`，仅耗时行不同）。代码侧登记项沿革再进一格：**P033 清掉 §9 的排序余项**，现只剩那条无基线的 STEP 写侧 harness。**P034 随即发现在本清单从未点名的一个类族上整体缺件**：3D 点-曲线**局部**极值族 —— `crates/occt-geom/src/extrema_pc/` 只有全局臂 `Extrema_GGExtPC`，没有 `Extrema_CurveTool` / `Extrema_GFuncExtPC` / `Extrema_GenLocateExtPC` / `Extrema_GLocateExtPC`，于是三处需要「已知参数附近求点-曲线最近极值」的调用点各自写了夹取 Newton：`int_tools_vertex_line.rs`（12 次迭代 + `|C'|^2 &lt;= EXTREMA_EPS^2` 门 + 位移判据）对 `IntTools_Context.cxx:842-871`、`wire_fix/same_parameter.rs`（12 次迭代 + `1e-20` 速度门 + `1e-10` 位移判据）对 `BRepLib_ValidateEdge.cxx:168-227`、`brep_lib_continuity.rs`（100 次迭代 + `has_analytic_d2` 前置使不满足的边整条静默跳过）对 `BRepLib.cxx:2341-2355`；cxx 的控制流是 `math_FunctionRoot` 再用 `|F| &gt;= 1.e-07` 复核（`Extrema_GenLocateExtPC.hxx:108-127`），自由曲线按 `GeomAbs_C2` 区间行走、初等类型跑 `Extrema_ExtPC` 取离 `U0` 最近解（`Extrema_GLocateExtPC.hxx:125-307`）。四类已全部移植、三处自创循环已删除，`extrema_ext_pc_range_done` 补上初等臂所需的 `IsDone`（`Extrema_GGExtPC.hxx:462-469`：解表可非空而 `mydone == false`）；门禁五 crate 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、`occt-geom2d --lib` 72/0、`occt-core --lib` 283/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4 与 `export_data_obj` 15/15 均与 P033 逐字一致（`.target-gate/t164/`）。本类另登记两条 `CurveOnSurface` 覆写缺失（`Resolution`、`gp_line`），均就地标 UNPORTED、未改。至此 `B8-intana` 一行的记账被推翻并修正五次（P027 按文件、P028 按目录、P029 订正普查、P030 清掉最后五臂、P034 才发现同类「按类族盘点」的口径对 `Extrema_*` 局部极值族也不成立）。
        </Text>
        <Text>
          <Text weight="semibold">偏移</Text>：本轮把缝边 pcurve 注册按 OCCT 全量重写（锥面缝再按 `FirstParameter()` 交换），
          TDB 三角 81,721 → <Text weight="semibold">79,926</Text>（OCCT 79,977）、面积比 1.0042 → <Text weight="semibold">1.0000</Text>，
          逐面配对 2,050 面逐位相同（此行为 2026-10-13 当轮读数）；`data/occ` 八模型 2026-10-16 复测为 `motoc` 逐位一致 + 双向小量
          （T0M −469、top −30、ATU01038 −32、a3n00 −24、acs10 −20、bottom −7、TDB +17），面积比 20×1.0000 / `Shape-2` 0.9999 / `T0M` 1.0017 / `TDB` 1.0015
          （10-13 记的「全 1.0000、T0M +58、TDB −51」已被后续共享路径忠实化改写，`.target-gate/t157/gates-raw.txt` 为准）。2.4 的七条状态/文档偏移（看板基准日、`G-lib`、`G-export`、`G-core`、
          未上板自创件、未提交工作树、10-10 未上板批次）已全部收口；**本卡于 2026-10-16 按新口径复测刷新**（Rust 侧 1,111 文件 /
          338,784 行、`UNPORTED` 155 条（大小写敏感计数，`crates/**/*.rs`），见 2.1 与 2.4）。
        </Text>
      </Stack>
      <Callout tone="success" title="T999 Judge 审计（本 tranche 结论：full_outcome_complete）">
        <Text>
          <Text weight="semibold">剩余登记项</Text>：`ReparamCurve2d` 的 `SameRange` 线已作为 P026 / T-156（commit `fbd59532`）收口
          —— 判词是「非 stale-closed，真实缺口是转换尾重复」：`pcurve_full/surface_projector/reparam.rs` 与 `shhealing/wire_fix/pcurve_range.rs`
          各带一份「`Geom2d_TrimmedCurve` → `Geom2dConvert::CurveToBSplineCurve` → `BSplCLib::Reparametrize` → `SetKnots`」副本
          （`GeomLib.cxx:912-921` / `:962-968`），且 `pcurve_range.rs` 的副本在 cxx 不判型的不等跨臂（`cxx:924-969`）按曲线类型分派
          （`line_trim_to_reparam_bspline` / `unequal_trim_window` → `same_range_unequal_line` / `same_range_unequal_other`）。
          现两份合并为唯一 `trimmed_reparam_bspline` 与 `same_range_window`（`cxx:934-959`），`ReparamCurve2d::bspline_copy2d` 删除
          （wrapper 只在 `Geom2dConvert::CurveToBSplineCurve` 拒绝时存在，`Geom2dConvert.cxx:442-444`），`pcurve_range.rs` 只留 cxx 的分支选择
          （早退 `cxx:854-859`、等跨判据 `:862`、直线 `:864-871`、圆 `:872-889`、裁剪递归 `:890-901`、一般臂 `:908-922`、不等跨 `:924-969`）。
          门禁：`--all-targets` 0 error / 413 warnings（lib 375 与基线同数）、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（118.1s，23 面积比在登记容差内：20×1.0000、`Shape-2` 0.9999、`T0M` 1.0017、`TDB` 1.0015）、
          `export_data_obj` 15/15 逐位一致；未新增测试。**T999 Judge 审计已在冻结 HEAD `fbd59532` 上跑完**（四条门禁读数：`--all-targets` 0 error、
          `occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 逐位一致），逐条对照 T001 清单后判定
          `full_outcome_complete: true`：B1..B8、B6a-h、B7-writer、B8-intana、`surface_fit`、`ReparamCurve2d` 全部有包可指。唯一仍登记的是缺失的
          `data/occ/*.st*` 驱动式 STEP 写侧 harness —— 写侧没有既有基线、本 tranche 不允许新写测试，故其结论是「本 tranche 范围外」（已登记，未静默丢弃）。
          **该结论随后被 P027 部分推翻（2026-10-16，commit `76893533`）**：T999 对 `B8-intana` 的「有包可指」只覆盖锥/锥与柱族臂，
          `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Pln)` 与 `IntAna_Int3Pln::Perform` 两臂从未与 cxx 对拍 —— P027 在其中查出三处自创 + 一段未移植
          （最小范数原点解 `p = (d1·(n2×v) + d2·(v×n1))/|v|²` 配 `v2 &lt; 1e-20` 门、`cxx:449-508` 小角原点重解缺失、`solve3` 的 `1e-14` 手写消元、
          自创 `Intersect3Pln::Line` 共线分支）并按 cxx 全部收口（门禁：五 crate 0 error、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、
          `step_obj_gates` 4/4（251.49s，23 面积比与 P026 逐位一致）、`export_data_obj` 15/15 逐位一致）。board 随之重开为
          `active_task: T005` / `goal.status: active`；新登记两条旁支（`plane_sphere` 的 `Epsilon(radius)` 相切判据、球/球臂的 `0.01 * CONFUSION` 圆判定），
          并由 T005 滚动推进。教训已记入 2.6：`IntAna_QuadQuadGeo` 这类多臂类必须**逐臂**核对，不能按文件记账。
          **该类的四条后续包随后把本行彻底重建**（2026-10-16）：P028 `1d595f8d` 收口 P027 登记的两条旁支；P029 `6330c6d7`
          订正 P028 的「5 臂未移植」普查（`intana/mod.rs:30-35` 的 `#[path = "../intana_torus.rs"] mod torus;` 让环面五臂在目录之外，
          正确读数是 15/15 在仓）并删掉环面两臂的 `a_dist &lt;= f64::EPSILON` 自创同心分支；P030 `0a185eb5` 逐行走完最后五臂，
          查出 Pln/Cylinder 缺 `H` 与 `cxx:568-598` 近平行细化臂、消费方 `int_p_cy` 恒按 `H = 0`，并把 Pln/Cone 的 `cxx:922-951`
          标 `UNPORTED`；**P031 `9bef99ad` 把那条 `UNPORTED` 落地** —— `QuadricIntersection::NotDone` 承载 cxx 的 `done = false`，
          `ELLIPSE_LIMIT = 1.0e9` / `HYPERBOLA_LIMIT = 2.0e6` 按原值，三条消费路径同步接上（`PairOutcome::Fail`、两处嵌套约化的 `None`、
          两条解析快速路径的无曲线回落）。四包门禁逐条：五 crate `--all-targets` 0 error（12/30/105/82/413 warnings）、
          `occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 逐位一致（P031 的 gates 日志
          与 P030 逐字一致，仅耗时不同）。**至此该类 15 个重载全部逐行对拍、全部在仓、无代码侧登记项**，代码侧登记转由同族的
          `IntAna_IntLinTorus` done 语义接手（列为 P032，audit `occt-geom.md` §9）；**P032 `55b0a6d6` 已于同日收口该项** ——
          `line_torus_intersect` 返回携带 `done` 的 `LineTorusIntersection`、函数体改持 `DirectPolynomialRoots::new5(..)`
          以分清 `cxx:89` 的 `mdpr.IsDone()` 与零根、`cxx:116-120` 与 `cxx:126-131` 两支都置 `done = false`，消费方
          `solvers.rs:456-465` 在 `!is_done()` 时 `return None`（= `IntCurveSurface_InterUtils.pxx:1290-1294`），门禁与 P031
          逐字节相同（`.target-gate/t162/`）。**P033 `a4eeac65` 随即清掉 §9 的 (b) 排序余项** —— `line_torus.rs:922` 的
          `out.sort_by(|a, b| a.x()...)` 在 `cxx:104-108` 无对应（OCCT 按根序号写入 `thePoint[nbsolvalid]`，函数体内不排序），
          删除后 `intana/mod.rs:12` 的 `Ordering` 再导出失去使用点、一并撤掉；消费方 `solvers.rs:210` 按曲线参数重排、
          两处测试各自再排序，故无依赖语义。门禁 4/4 与 15/15 均与 P032 逐字一致（`.target-gate/t163/`，仅耗时行不同）。
          **`occt-geom.md` §9 的 (a)/(b) 两项至此全部闭环**，代码侧登记项只剩那条无基线的 STEP 写侧 harness。
          **P034 `c272c08a` 随即在另一个类族上整体补件（2026-10-16）**：3D 点-曲线**局部**极值族此前完全缺席 ——
          `crates/occt-geom/src/extrema_pc/` 只有全局臂 `Extrema_GGExtPC`，没有 `Extrema_CurveTool` / `Extrema_GFuncExtPC` /
          `Extrema_GenLocateExtPC` / `Extrema_GLocateExtPC`，于是三处需要「已知参数附近求点-曲线最近极值」的调用点各自写了夹取
          Newton（`int_tools_vertex_line.rs` 12 次迭代 + `|C'|^2 &lt;= EXTREMA_EPS^2` 门 + 位移判据，对 `IntTools_Context.cxx:842-871`；
          `wire_fix/same_parameter.rs` 12 次迭代 + `1e-20` 速度门 + `1e-10` 位移判据，对 `BRepLib_ValidateEdge.cxx:168-227`；
          `brep_lib_continuity.rs` 100 次迭代 + `has_analytic_d2` 前置使不满足的边整条静默跳过，对 `BRepLib.cxx:2341-2355`）。
          四类已按 cxx 全部移植（自由曲线走 `math_FunctionRoot` 再以 `|F| &gt;= 1.e-07` 复核、初等类型跑 `Extrema_ExtPC` 取离 `U0`
          最近解）、三处自创循环与它们的私有助手全删，`extrema_ext_pc_range_done` 补上初等臂所需的 `IsDone`（`Extrema_GGExtPC.hxx:462-469`：
          解表可非空而 `mydone == false`）。门禁：五 crate `--all-targets` 0 error（12/30/105/82/413 warnings）、`occt-geom --lib` 143/0、
          `occt-geom2d --lib` 72/0、`occt-core --lib` 283/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15
          （后两者与 P033 逐字一致，`.target-gate/t164/`）。本类另登记两条 `CurveOnSurface` 覆写缺失（`Resolution`、
          `Adaptor3d_CurveOnSurface.cxx:1364-1370`；`gp_line`、`EvalKPart` 的 `GeomAbs_Line` case `cxx:1565`/`:1687`/`:1701`），
          均就地标 UNPORTED、未改。**至此 `B8-intana` 一行的记账被推翻并修正五次**，教训升级为：清单必须按**类族**盘点，
          不按文件、不按目录、也不按单一类（P034 找到了一个此前从未被点名的族）。
        </Text>
        面级残余的定位受「可信的 OCCT 面级预言机」未钉住所限：实测同一条 `--facestats 2.45132 0.349066`
        在探针上是 `72683 / 78464`，而冻结参考是 `73544 / 79977`、端口是 `73516 / 79926` —— 探针用的是
        `STEPControl_Reader` 直读，与冻结参考的产出者（未记录）以及端口的不同源，故探针总量不可直接当日面级基线；
        要么先钉住与参考同源的读取+参数组合，要么只以冻结 OBJ 的聚合口径为规格。
        在此之前可推进的独立项是 2.2 里的自创残余：`heal_tolerance` 自创修形链（本轮已收口：`FixSmall`/`CheckSmall`/`FixSmallAreaWire`/`CopyReplaceVertices` 落地，`boolean_with_check` 容差重算删除）、`Geom_RectangularTrimmedSurface::SetTrim` 族与 `gcpnts` 的 `GCPnts_AbsComposite` adv 臂（两轮均本轮收口，见 2.1）均已闭环；`bop_curved` 网格级布尔（A5/T-79）已随引擎删除闭环，`brepfeat` 的圆柱夹具（A12 余项 / T-125）本轮已换 `BRepPrimAPI_MakeCylinder` 且体积精确（`box ∪ boss` 4.157080、`box ∪ 贯通 boss` 2.565487），现存自创项只剩 `revolve_profile_about` 的解析旋转体工具 —— 缺口经本轮实测不在 `BOPAlgo_BOP`，而在 `sweep_revolve` 未按 `BRepSweep_Revol`/`BRepPrim_OneAxis` 的 wire/seam/pcurve 约定建面（见 2.2 与 `_audit/_index.md` A12）。本轮第二个独立项是「用采样识别二次曲面代替 downcast」（`_audit/topo-rest.md` 发现 14）：偏移面分派已换精确 downcast —— `offset_face` 平面/球臂改 `gp_pln()`/`gp_sphere()`（见 2.1 首条），`brep_surface::sphere_center` 改 `GeomAdaptor_Surface::Sphere()` 等价并删除采样外接圆心族；残留的 `brep_surface::is_planar`/`face_is_planar` 调用点第二批 9 处也已换 `is_plane_type`（BOP 交叉 + EF 家族，逐条 `.cxx` 行号见 2.2 发现 14 余项），第三批（T-129）换完余 19 处并删除采样体，`edge_face_kind` 的曲线/曲面族亦换成 `GetType()` 同等 downcast ⇒ 发现 14 已整体收口；批次 4（T-130）把同族的 `Resolution` 采样估计（`edge_face_kind::curve_resolution`、`inttools_range` 3 处、`algo_tools_range` 的共线/自由形判据、`bean_face::resolution`、`int_tools_vertex_line::adaptor_resolution`、`wire_splitter_block::angle_2d`）全部换成 `Adaptor::Resolution`／`Curve2d::resolution` 解析闭式，并在 `occt-geom2d` 补上 `Curve2d::resolution`（`Geom2dAdaptor_Curve::Resolution` 等价，逐型 Line/Circle/Ellipse/Bezier/BSpline/Trimmed），`ReparamCurve2d`/`PlaneKeepParam2d` 两个装饰器补 `resolution` 以免退回 trait default。剩余同族残留只剩 `inttools_range::find_nearest_valid_point` 缺 `BRepLib_1.cxx:69-83` 的 `aD1Mag` 局部奇点加速分支（原注释已标 `ponytail:`）—— 该项亦已于 T-133（2026-10-15）收口（见上 B5 / T-133 行）。
      </Callout>
    </Stack>
  );
}
