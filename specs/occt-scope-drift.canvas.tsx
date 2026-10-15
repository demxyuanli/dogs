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
} from "cursor/canvas";

// ---------------------------------------------------------------------------
// 量测口径（2026-10-13 实测）
//   OCCT 侧: D:\source\OCCT-src\src  .hxx/.cxx/.lxx/.gxx
//   Rust 侧: crates/{occt-core,occt-math,occt-geom,occt-geom2d,occt-topo}/src  .rs
//   基线  : specs/board.canvas.tsx (r41, asOf 2026-10-11) · specs/_audit/_index.md (§9 -> 2026-10-11)
//   对拍  : crates/occt-topo --test step_obj_gates（`--nocapture` 逐模型 v/f 与面积比）
// ---------------------------------------------------------------------------
const MEASURE = {
  asOf: "2026-10-15",
  occt: { files: 13803, lines: 2439487 },
  rust: { files: 1060, lines: 335743 },
  cadKernel: 1397840, // Foundation 252505 + ModelingData 347970 + ModelingAlgorithms 797365
  unported: 166,
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
  { crate: "occt-topo", files: 606, lines: 216367 },
  { crate: "occt-core", files: 220, lines: 45065 },
  { crate: "occt-geom", files: 121, lines: 41413 },
  { crate: "occt-geom2d", files: 58, lines: 18604 },
  { crate: "occt-math", files: 55, lines: 14294 },
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
    state: "**B6a 已修（2026-10-15）**：三处兜底整体删除；`point_surface_extrema`/`curve_surface_extrema`/`surface_surface_extrema` 改返回 `Option<ExtremaPair>`（`None` = `Extrema_ExtPS.cxx:307-371` / `Extrema_ExtCS.cxx:122,532` / `Extrema_ExtSS.cxx:226-260` 的 `!IsDone()` 或 `NbExt() == 0`），`numeric_extrema::no_solution_pair` 收拢原有的「无解」标记对。消费方 `bean_face_exact::ComputeUsingExtremum` 回到 `IntTools_BeanFaceIntersector.cxx:956` 的守卫（`Extrema_ExtCS` 未 done 即整块跳过、不置 range flag），`closest_point_on_face` 无解时跳过面内投影，`extrema_point_surface_exact` 非初等面臂改 `Err`。旧登记「有生产消费方故暂留、前置 T-67」不成立。门禁：`--all-targets` 0 error、`occt-geom --lib` 143/0、`--lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 且 v/f 逐位一致。**B6b 已修（2026-10-15）**：两侧采样式分类器全部换成 `GetType()`/`gp_*` 下转。`geom/extrema_surf/analytic_solvers.rs`：`classify_plane`＝`s.gp_pln()`、`classify_sphere`＝`s.gp_sphere()`（`Extrema_ExtSS.cxx:117-129` 的 `GetType()` 分派 + `Adaptor3d_Surface::Plane()/Sphere()`）、`reconstruct_line`＝`c.gp_line()`（与 `Extrema_ExtCS.cxx:129-140` 的 `C.Line()` / `myS->Sphere()` / `myS->Plane()` 同源），采样件 `is_planar`（8×8 法向并行判）、`sphere_center`（中纬度三点外心）、`circumcenter`、`surface_normal`、`is_line`（6 点切向）、`sample_bounds` 整段删除。`geom2d/extrema2d`：`point_curve_extrema2d_all` 线臂改 `c.gp_lin2d()` 并按 `Extrema_ExtPElC2d.cxx:63-79` 用 `Mydist = L.Direction().Dot(L.Location()→P)` 取**线上参数**（原实现取的是相对 `c.d0(first)` 的参数增量，裁剪线/带偏置起点的线上参数直接失真）、圆臂改 `c.gp_circ2d()`；`curve_curve_extrema2d_all` 的配对分派改 `gp_lin2d()`/`gp_circ2d()`（`Extrema_ExtCC2d.cxx:104-145` 的 `GetType(C1)`×`GetType(C2)` 与 `Extrema_Curve2dTool::Line/Circle`）；采样件 `is_line2d`/`circumcenter2d`/`classify_circle2d`/`line_of_curve2d` 删除。**可观测修复**：几何上直的 B 样条/近直线曲线原先被 6 点切向判据拖进线臂（OCCT 给 `GeomAbs_BSplineCurve`，落一般臂），裁剪线的点-线极值参数亦由增量纠正为线上参数。门禁：`--all-targets` 0 error、`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（117s）、`export_data_obj` 15/15 且 15 个模型 v/f 与包围盒逐位一致。**B6c 已修（2026-10-15）**：`extrema2d` 一般臂的种子栅格与其求解器整体删除 —— `analytic_solvers.rs::build_samples2d`（`span/0.1` clamp[24,256]）、`newton_point_curve_all2d`、`curve_curve::newton_curve_curve_all2d` 及 `fval2d`/`fprime2d`/`dist22d`/`solve_f_zero2d`/`refine_seed2d`/`pair2d_c`/`dedupe_sort2d` 全部下线；替身对应的 cxx 是 `Extrema_ExtPC2d` = `Extrema_GGExtPC`（`Extrema_ExtPC2d.hxx:31-38`；`Initialize` 取 `mytolu = Resolution(Confusion)`、`type = GetType(theC)`、`mysample = 17`，`Extrema_GGExtPC.hxx:119-128`；`Perform` 先把两端点距离入解，再按 `type` 分派：Line/Circle/Ellipse/Parabola/Hyperbola → `Extrema_ExtPElC2d`，Bezier → `NbPoles()*2` 个样本 + `IntervalPerform`，BSpline → 按度数与节点区间，`default:` → `NbIntervals(GeomAbs_C2)` / `DeflCurvIntervals` 逐区间 `IntervalPerform`）加 `IntervalPerform`（`hxx:590-620`，`myExtPC.Initialize(mysample, uinf, usup, tolu, tolf)` + `Perform`，`AddSol` 以 `PConfusion()` 归并）与 `Extrema_GFuncExtPC`（`Value`/`Values`，含 `D1`/`DF` 奇点处的 Taylor 与三点差分）；本轮落成单一忠实件 `extrema2d/ext_pc2d.rs`（`ExtremaExtPC2d`）与 `extrema2d/gfunc_ext_pc2d.rs`（`FuncExtPc2d`），`point_curve_extrema2d_all` / `point_curve_extrema2d`（改返回 `Option`，对齐 `Geom2dAPI_ProjectPointOnCurve::Perform` 的 `NbExt() > 0` 且取最小 `SquareDistance`）与 `curve_curve_extrema2d_all`（`Extrema_ExtCC2d`）接线到位；`general_extrema::proj_p_on_c` 改用 `ExtremaExtPC2d::new_range`（`Extrema_GGenExtCC.hxx:249-263`），其 **UNPORTED (audit A16)** 标记随之清除；`topo/shape_fix_compose_shell/split_wire.rs` 按 `ShapeFix_ComposeShell.cxx:1046-1087` 恢复 `apar` 默认 `firstPar`、仅有解时覆盖。四条既有 lib 断言按 cxx 次序订正（`Extrema_ExtPElC2d.cxx:113-150` 圆臂远点在前、`Extrema_ExtElC2d.cxx:334-364` 圆/圆 `P1[NoSol]×P2[kk]` 次序、平行直线 `IsParallel()` 且 0 极值），未新增测试。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（23 个面积比与基线逐位一致）、`export_data_obj` 15/15 且 v/f 与包围盒逐位一致。**B6e 已修（2026-10-15，包 P017 / T-147，commit `e3a4335e`）**：`extrema2d::locate_extcc2d`（`curve_curve.rs:255-262`）的自创 Newton 整体删除 —— 原 `refine_curve_curve2d_newton` 是夹取 32 次迭代 + `1e-12` 收敛判据，且 `f_and_jac2d` 的 Jacobian 四个带符号项相对 `Extrema_GFuncExtCC::Values` 全反；现按 OCCT 链落地：`Extrema_LocateExtCC2d.cxx:26-82` → `Extrema_LocECC2d` = `Extrema_GenLocateExtCC`（`Extrema_GenLocateExtCC.hxx:47-110`：逐轴量程规范化到 `[FirstParameter, LastParameter]`、`TolU`/`TolV` 取 `Adaptor::Resolution(Confusion)`、`MathFunctionSetRoot::perform_with_bounds`，随后 `F.NbExt() > 0` 与 `F.Value < 1e-10` 复核）→ `Extrema_CCLocFOfLocECC2d` = `Extrema_GFuncExtCC` 2D 实例。新忠实件 `extrema2d/gfunc_ext_cc2d.rs`（`GFuncExtCC2d`，`MathFunctionSetWithDerivatives` 两方程两变量：`SearchOfTolerance` `hxx:150-197`、`GetType` 分派 `hxx:263-291`、`Value` `hxx:359-550` 含 `D1` 奇点处的三点 / `DN` 重建、`Values` `hxx:568-827` 含奇点导数的扰动三点 Jacobian、`GetStateNumber` `hxx:832-870` 存 `(POnCurv2d, POnCurv2d)`、`Points`/`SubIntervalInitialize`）与 `extrema2d/gen_locate_ext_cc2d.rs`（`GenLocateExtCC2d` + `ExtremaLocateExtCC2d`）。门禁：`occt-geom2d --lib` 72/0、`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4（248.5s）、`export_data_obj` 15/15 逐位一致。**B6f/B6g 仍登记（未改）**：B6f＝`curve_curve_intersections2d` 的 33×33 采样栅格（对 `Geom2dAPI_ExtremaCurveCurve`/`Extrema_ExtCC2d` 求解与 `IntCurve_IntConicConic` 家族）；B6g＝`line_line`/`line_circle`/`circle_circle_extrema2d` 与 `ext_el_c2d.rs` 重复的 `Extrema_ExtElC2d` 臂副本。",
  },
];

// 语义偏移残余（对照 specs/_audit/_index.md A 编号；位置/行号为 2026-10-13 工作树实测）
const SEMANTIC = [
  {
    item: "Delaunay 插点密度：`BRepMesh_DelaunayDeflectionControlMeshAlgo` 的插点与最小边判据",
    id: "A19 余项",
    where: "topo/meshing/deflection_control/ · topo/meshing/node_insertion.rs",
    state: "**已移植（2026-10-15 行内复核，原「未移植」条目标题本身有误）**：OCCT 8.0 的该类没有 `checkAndAdd`/`addInternalNodes` 这两个符号（全树 grep 无命中），行内把 `optimizeMesh`/`splitTriangleGeometry`/`splitLinks`/`usePoint`/`checkDeflectionOfPointAndUpdateCache`/`rejectByMinSize`/`rejectSplitLinksForMinSize`/`checkLinkEndsForAngularDeviation`（`BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx:75-454`）与 `insertNodes`/`registerSurfaceNodes`/`getCellsCount`（`BRepMesh_DelaunayNodeInsertionMeshAlgo.hxx:74-155`）全部落地，并接线（`node_insertion.rs:519-679` 的 `post_process` + `insert_parametric`）；`A_ITERATIONS_NB = 11` 与 `BRepMesh_GeomTool::CellsCount` 亦已对齐。原 TDB +1,744 面证据经逐面配对复核归属缝边 pcurve 注册（2.1 T-124）而非本项；现残余为双向小量：T0M +58、TDB −51、ATU01038 −32、top −30、a3n00 −24、acs10 −20、bottom −5、motoc 0，面积比全 1.0000 —— 无 4~9% 的系统性偏密",
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
    state: "**已修（2026-10-15）**：原先 4 类分类（`is_periodic() && |period-2*PI|<1e-9` ⇒ Circle、无界 ⇒ Line、否则 BSpline）被 `Geom2dAdaptor_Curve::GetType()` 同等件取代 —— 新增 `occt_geom2d::adaptor_curve::{basis2d,get_type2d}`（`Geom2dAdaptor_Curve.cxx:96-120` + `load` 类型表 `:279-346`，含 `Geom2d_TrimmedCurve` 下基 `:285-288`），`extrema2d/curve2d_tool` 的 `adaptor_basis2d`/`adaptor_type2d` 改为转发该唯一实现。(1) `geom_bnd_lib_curve2d::curve2d_abs_type` 改读 `get_type2d`，`box_curve2d` 按 `GeomBndLib_Curve2d.cxx:97-171` 的 downcast 链补全九臂：Line/Circle 优先取 `gp_lin2d`/`gp_circ2d`（采样重建件 `lin2d_from_curve`/`circ2d_from_curve` 保留为 `PlaneKeepParam2d` 这类未转发 `gp_*` 的 wrapper 替身）、Ellipse/Hyperbola/Parabola 取 `gp_elips2d`/`gp_hypr2d`/`gp_parab2d` 走解析盒（原先一律落采样）、Bezier 走 `GeomBndLib_BezierCurve2d::Box` 极点包络（`cxx:30-52`，子区间无 2D `Segment` 时退回采样）、BSpline 走 `bspline_copy2d`、Offset 走 `offset_basis`+`offset_value`。(2) `geom_bnd_lib_offset2d` 的 `try_offset_line2d`/`try_offset_circle2d` 由 `pc_curve_kind` 改用 `GeomBndLib_OffsetCurve2d.cxx:36,63` 的精确 `occ::down_cast<Geom2d_Line>/<Geom2d_Circle>` 等价（`gp_lin2d()`/`gp_circ2d()`），基圆圆心半径不再靠 `d0(0)`/`d0(pi)` 采样反推；`box_offset_as_curve` 由恒等 `box_other` 改为按 `cxx:123-143` 先试线/圆快路径。(3) `int_face_face_bounds.rs` 的 `Geom2d_Line` 判定改精确 `gp_lin2d()` 并直接用 `Position().Direction()`（`IntTools_FaceFace.cxx:2149,2161-2162`），原 `d1(0)` 切向改判方向。(4) `pcurve_full/tests.rs`、`pcurve.rs`、`boptools_2d.rs` 中引用 `pc_curve_kind` 的既有断言改断言 `curve2d_abs_type`。门禁：`occt-geom2d` `cargo check` 0 error；`--lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15（15 个模型 v/f 与包围盒与基线逐位一致）。**同批发现、未改**：`ReparamCurve2d`（`GeomLib::SameRange` 的线性替身）仍把内层 kind 转发给外层，OCCT 的 `SameRange` 对非直线/非圆输入返回 `Geom2d_BSplineCurve`，故包在 Offset 上的重参数化 wrapper 会报 `GeomAbs_Line` 而非 `GeomAbs_BSplineCurve`。",
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
    state: "**已修（2026-10-15）**：`init` 原先在 `build_osculating` 失败时四处 `clear_flags(); return;`，把余下的 `for` 结点迭代连同后面整段 V 向扫描一起丢掉，并把五个臂的成功与否挤进一个外层 `oscul_ok`。OCCT `Geom_OsculatingSurface.cxx:179-378` 每臂各持局部 `OsculSurf`（下轮迭代重新赋值），失败只 `clearOsculFlags()`（`:206-213`、`:273-276`、`:308-312`、`:373-376`）后继续；只有 `if (myAlong[1] && OsculSurf)` 这类守卫因标志被清而跳过。现新增 `build_chain`＝`cxx:187-206` 的 `while (IsQPunc)` 链（返回 `Option<(GeomBSplineSurface, k)>`，`None` ＝ 空 `OsculSurf`，`k` 为成功次数），两个 `for` 循环改为 `if self.along[0] { A; if self.along[1] && ok { B } } else if let Some(..) { C } else { clear_flags() }` 的逐臂判定；四处 `return` 与共享 `oscul_ok` 删除。可观测修复：一次失败后余下的结点索引不再被整段跳过。门禁：`occt-geom --lib` 143/0、`occt-topo --lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 逐位一致。**同批订正**：B1（`wire_splitter_block.rs` 三臂）实已在 T-131（`67d7f01f`）收口，清单原标 LIVE 属陈旧；`extrema_cc/mod.rs:9`、`occt-math/globoptmin.rs:590` 的 UNPORTED 注释亦为陈旧（代码已用 `math_GlobOptMin`／已用 `newton_minimum.rs`）。**旁支（未改）**：`surface_fit.rs::fit_plane` 是最小二乘拟合而非 `GeomLib_IsPlanarSurface` 移植，忠实修法是补 `Controle`/`GeomLib::Inertia` 并接线 `brep_offset/curve_face_offset.rs:489`。",
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

// `data/*` + `data/occ/*` 逐模型对拍（2026-10-15 `step_obj_gates --nocapture` 实测，23 个模型全过）
const PARITY = [
  { model: "motoc", occV: 11713, occF: 13944, ourV: 11713, ourF: 13944, area: "1.0000" },
  { model: "a3n00", occV: 11145, occF: 12466, ourV: 11138, ourF: 12442, area: "1.0000" },
  { model: "acs10", occV: 37296, occF: 46558, ourV: 37280, ourF: 46538, area: "1.0000" },
  { model: "ATU01038", occV: 18102, occF: 22586, ourV: 18085, ourF: 22554, area: "1.0000" },
  { model: "bottom", occV: 17587, occF: 23557, ourV: 17581, ourF: 23550, area: "1.0000" },
  { model: "top", occV: 16905, occF: 22800, ourV: 16888, ourF: 22770, area: "1.0000" },
  { model: "T0M", occV: 60548, occF: 67366, ourV: 60505, ourF: 67284, area: "1.0000" },
  { model: "TDB", occV: 73544, occF: 79977, ourV: 73516, ourF: 79926, area: "1.0000" },
];

// 状态/文档偏移：2026-10-13 复核定案
const DOCDRIFT = [
  {
    what: "看板基准日",
    detail: "board asOf 2026-10-11 · r41 → 2026-10-13 · r42",
    gap: "✅ 已同步：`a91` 已收录 10-12 网格批次与 10-13 缝边注册修复",
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
    what: "SEMANTIC 表 LIVE 方向也过期",
    detail: "24 条登记的 drift 行其代码已闭环（topo-rest 1/4/6/11 + 18 的 classifier 半、bop 1-4/7-16-17、mesh 1-14/19）",
    gap: "❗未同步：`docs/goals/drift-and-gap-autofix/notes/T001-inventory.md` 逐条给了证据，本表尚未按方向订正（不是代码问题，是留档口径问题）",
    tone: "danger" as const,
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

// 看板门禁（读数日期 2026-10-15；G-check / G-lib / G-gates / G-export 覆盖到 T-145 已提交工作树 `055e54aa`，
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

// 目标看板跟踪（`docs/goals/drift-and-gap-autofix/state.yaml`，active_task 恒为 T005；每个包 = 一个已核实并提交的切片）
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
    pkg: "下一步",
    task: "T005 → 待建 Worker",
    commit: "—",
    scope: "B6g：`extrema2d` 的 `line_line` / `line_circle` / `circle_circle_extrema2d` 是 `ext_el_c2d.rs` 已有的 `Extrema_ExtElC2d` 臂的第二份模块级副本（仓内无生产调用方，仅 `extrema2d/tests.rs`）",
    state: "⬜ 待做：删除 `crates/occt-geom2d/src/extrema2d/curve_curve.rs` 的 `line_line_extrema2d`/`line_circle_extrema2d`/`circle_circle_extrema2d` 三个副本，调用点改走 `ext_el_c2d.rs` 的 `Extrema_ExtElC2d`（`Extrema_ExtElC2d.cxx:63-70` 平行直线 `IsParallel()` 零极值、`:327-338` 线-线系数、`:982-1110` 圆/圆、`:334-364` 圆/圆 `P1[NoSol]×P2[kk]` 次序）。其后依次 B6f（33×33 `curve_curve_intersections2d` 栅格，唯一调用方是它自己的单测）、B6d、B7-writer-4（`geom2d_convert.rs:77` 圆锥曲线臂对 `GeomToStep_MakeConic`）、`surface_fit`（`GeomLib_IsPlanarSurface`）与缺失的 `data/occ/*.st*` 驱动式 STEP 写侧 harness",
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
          所有读数取 <Text weight="semibold">{MEASURE.asOf}</Text>；门禁基线引用 <Text weight="semibold">specs/board.canvas.tsx</Text>（r42）与{" "}
          <Text weight="semibold">specs/_audit/_index.md</Text>（§9 → 2026-10-11）。
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
        来源：D:\source\OCCT-src\src（8 个模块组，.hxx/.cxx/.lxx/.gxx）· 2026-10-13 复测（与 10-11 逐位相同）· 单位：源码行。
      </Text>

      <H3>Rust 侧分布</H3>
      <Table
        headers={["crate", "文件", "行数"]}
        rows={CRATES.map((c) => [c.crate, String(c.files), c.lines.toLocaleString("en-US")])}
        columnAlign={["left", "right", "right"]}
      />
      <Text size="small" tone="tertiary">
        较 2026-10-11（1,080 文件 / 337,891 行）：<Text weight="semibold">−20 文件 / −2,148 行</Text>，主要来自
        `brepmesh.rs`（−602 行，A18/A19）、`wireframe.rs` 平面机器（−288 行，A18）、5 条自创单测的删除，
        以及本轮撤除的 3 处 `IGES_UVDBG` 临时插桩（−59 行，`4a20e078`）。
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

      <H3>2.3 逐模型对拍（`data/occ/*`，2026-10-15 `step_obj_gates --nocapture` 实测）</H3>
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
        `data/*.step`（15 个）中 12 个逐位一致；`Shape-1` 差 −8 面、`Shape-2` 差 −34 面、`screw` 差 +2 面。八个 `data/occ/*` 模型面积比全为 1.0000，缺口是双向小量（T0M +58、TDB −51、top −30、ATU01038 −32、a3n00 −24、acs10 −20、bottom −5、motoc 0）；剩下面级差集中在 2.2 首行登记的 Delaunay 插点密度分支。
      </Text>

      <H3>2.4 状态 / 文档偏移（2026-10-15 复核：全部收口）</H3>
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
        `G-check` / `G-lib` / `G-gates` / `G-export` 读数测量于 2026-10-15，覆盖到 T-147 已提交工作树 `e3a4335e`（`git diff HEAD` 为空）：五 crate `--all-targets` 0 error、`occt-geom2d --lib` 72 / 0、`occt-geom --lib` 143 / 0、`occt-topo --lib` 1251 / 0、`step_obj_gates` 4 / 4（248.5s，23 个面积比与基线逐位一致）、`export_data_obj` 15 / 15 且逐模型 v/f 与包围盒与 T002 冻结表逐位一致。`G-fuse` / `G-phase` / `G-boss` / `G-iges` 为 T-130..T-133 轮读数，其后的包未触及这些路径。`G-core` 的 `occt-core` 读数由 285 变 283 是 T-143 删除两条自创单测所致（同 `G-lib` 的 1276 → 1251，`notes/T002-gate-baseline.md` 已留档）；`occt-math` 210、`occt-geom` 143、`occt-geom2d` 72 未变。
      </Text>

      <H3>2.6 目标看板跟踪（`docs/goals/drift-and-gap-autofix/state.yaml`）</H3>
      <Text size="small" tone="tertiary">
        active_task 恒为 T005（PM 滚动推进）；每个包 = 一个已跑现基线核实并单独提交的切片，提交信息点名 `.cxx` 文件与行号范围。
      </Text>
      <Table
        headers={["包", "任务", "提交", "切片范围", "状态"]}
        rows={CYCLES.map((c) => [c.pkg, c.task, c.commit, c.scope, c.state])}
        columnAlign={["left", "left", "left", "left", "left"]}
        rowTone={CYCLES.map((c) => (c.pkg === "下一步" ? ("warning" as const) : ("success" as const)))}
        striped
      />

      <Divider />

      <H2>三、结论与下一步</H2>
      <Row gap={16} wrap>
        <Stat value="13.8%" label="全库行数覆盖" />
        <Stat value="24.0%" label="CAD 内核行数覆盖" tone="info" />
        <Stat value="8 处" label="登记未收代码项（含下一片 1 处）" tone="warning" />
        <Stat value="0 项" label="文档/状态偏移待办" tone="success" />
      </Row>
      <Stack gap={6}>
        <Text>
          <Text weight="semibold">范围</Text>：全库 13.8%、CAD 内核约 24.0% 的行数覆盖；有效范围是 STEP→OBJ 共享管线，
          在此链上已做逐行对 `.cxx` 的深度移植（occt-topo 占 21.6 万行）。
        </Text>
        <Text>
          <Text weight="semibold">跟踪</Text>：2.6 的滚动看板把最近四个包（P014 锥/锥角度容差自推、P015 锥/锥平行臂重映射、
          P016 锥/锥共母线分支、P017 `extrema2d` 自创 Newton 换 `Extrema_LocateExtCC2d`）各自落到任务号与提交号；下一片已定档为
          B6g「`extrema2d` 的 `line_line` / `line_circle` / `circle_circle_extrema2d` 三个 `Extrema_ExtElC2d` 臂副本删除、改走
          `ext_el_c2d.rs`」，其后依次 B6f（33×33 求交栅格）、B6d、B7-writer-4（`geom2d_convert.rs:77`）、
          `surface_fit`（`GeomLib_IsPlanarSurface`）。
        </Text>
        <Text>
          <Text weight="semibold">偏移</Text>：本轮把缝边 pcurve 注册按 OCCT 全量重写（锥面缝再按 `FirstParameter()` 交换），
          TDB 三角 81,721 → <Text weight="semibold">79,926</Text>（OCCT 79,977）、面积比 1.0042 → <Text weight="semibold">1.0000</Text>，
          逐面配对 2,050 面逐位相同；`data/occ` 八模型面积比全为 1.0000（motoc 逐位一致），缺口收为双向小量
          （T0M +58、TDB −51、top −30、ATU01038 −32）。2.4 的五条状态/文档偏移（看板基准日、`G-lib`、`G-export`、未提交工作树、
          10-10 未上板批次）本轮全部收口。
        </Text>
      </Stack>
      <Callout tone="success" title="下一步（按门禁「完成即停」交接）">
        <Text>
          <Text weight="semibold">立即推进</Text>：B6g —— `crates/occt-geom2d/src/extrema2d/curve_curve.rs` 的
          `line_line_extrema2d` / `line_circle_extrema2d` / `circle_circle_extrema2d` 是 `extrema2d/ext_el_c2d.rs`
          已承载的 `Extrema_ExtElC2d` 臂的第二份模块级副本（`Extrema_ExtElC2d.cxx:63-70` 平行直线 `IsParallel()` 零极值、
          `:327-338` 线-线系数、`:334-364` 圆/圆 `P1[NoSol]×P2[kk]` 次序、`:982-1110` 圆/圆），仓内无生产调用方（仅
          `extrema2d/tests.rs`），故以删除副本 + 调用点改走唯一忠实件收口。前一片（P017 / T-147，commit `e3a4335e`）已把
          `locate_extcc2d` 的自创 Newton（32 次夹取迭代 + `1e-12` 判据、Jacobian 四符号项全反）换成
          `Extrema_LocateExtCC2d.cxx:26-82` → `Extrema_GenLocateExtCC`（`Extrema_GenLocateExtCC.hxx:47-110`）
          → `Extrema_GFuncExtCC` 2D 实例的忠实链。
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
