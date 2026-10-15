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
      "（已全部收口，2026-10-15）(1) `sphere_center` 改 `GeomAdaptor_Surface::Sphere()` 等价（`s.gp_sphere()`，`GeomAdaptor_Surface.cxx:1997-2006`；`gp_sphere` 已按 `cxx:423-425` 递归裁剪面基面，`Geom_OffsetSurface` 给 `GeomAbs_OffsetSurface` 不作答，`cxx:496-523`），`det3`/`solve3`/`circumcenter` 采样族整段删除；`offset_face` 平面臂改 `gp_pln()`（法向取 `gp_Pln::Position().Direction()`，与 `BRepOffset.cxx:54-61` 的 `XDirection ^ YDirection * Offset` 同源）、球臂改 `gp_sphere()` 并直接取 `Radius()`。原实现把**几何上平的 B 样条面**判成平面臂（OCCT 给 `GeomAbs_BSplineSurface`，落一般臂），且 `sphere_center` 的 `Some` 兜底把**锥/柱面**拖进球臂做一次以伪圆心为中心的 `SetScale`，面边界被错误缩放。(2) 新增 `is_plane_type`/`face_is_plane_type`（`brep_surface.rs:107,112`，即 `GetType() == GeomAbs_Plane`）后分两批替换采样调用点：批次 2（2026-10-15）换 BOP 交叉与 EF 家族 9 处（`pave_ff_perform` / `pave_intersect/face_face` / `pave_ef_perform` / `pave_force_ef` / `pave_pcurves` / `edge_face`），并把 `is_plane_ff` 按 `_6.cxx:84-102` 重写为真正的 `IsPlaneFF`；批次 3（T-129）换完余 19 处（`face_plane`、`pcurve.rs:249`、`fillet_curved`、`gprop_analytic`、`brep_offset/shell_offset.rs`、`brep_faces` 及测试夹具），并把 `brep_surface::is_planar`/`face_is_planar` 的 8×8 采样体整段删除，全仓已无该符号。(3) 批次 3 同族把 `edge_face_kind.rs` 的曲线/曲面族也从采样换成 `GetType()` 同等 downcast：`curve_kind`＝`gp_line()`/`gp_circ()`/`nurbs_degree()`（与 `GeomAdaptor_Curve::load` `GeomAdaptor_Curve.cxx:252-311` 同源）、`circle_geometry`＝`CurveTool::Circle`（`IntTools_EdgeFace.cxx:800,829`）、`surface_kind`→`classify_surface`、`is_eq_distance`＝`cxx:240-299`（frame 取 `gp_cylinder()`/`gp_cone()`/`gp_torus()`，锥顶按 `gp_Cone::Apex()`）、`is_coplanar`/`is_radius`＝`cxx:787-811`/`cxx:816-840`；采样件 `is_line_like`/`is_circle_like`/`det3`/`solve3`/`circumcenter`/`cylinder_matches`/`cone_matches`/`torus_matches`/`quadric_u0`/`midpoint` 删除。门禁：`--lib` 1251/0、`--tests` 全 `ok`（`step_obj_gates` 4/4、`brep_offset` 12/12）、`export_data_obj` 15/15 且 v/f 逐位一致（`Sphere` 642/1244、`Torus` 1369/2592、`Offset` 712/892、`Shape` 6150/11372）。(4) 批次 4（T-130，2026-10-15）把 `Resolution` 采样族整体换成 `Adaptor::Resolution`：`occt-geom2d` 新增 `Curve2d::resolution`＝`Geom2dAdaptor_Curve::Resolution`（`Geom2dAdaptor_Curve.cxx:1186-1219`，Line/Circle/Ellipse/Bezier/BSpline/Trimmed 逐型，default `Precision::Parametric`）；`edge_face_kind::curve_resolution`(`:73`)、`inttools_range::curve_resolution`（3 处）、`algo_tools_range` 的 `curve_is_line`+`curve_is_freeform`+`curve_resolution`、`bean_face::resolution`、`int_tools_vertex_line::adaptor_resolution`、`wire_splitter_block::angle_2d` 全部改走解析闭式（`IntTools_EdgeFace.cxx:326-327,373`、`IntTools_ShrunkRange.cxx:162`、`BRepLib_1.cxx:61,202`、`BOPTools_AlgoTools_2.cxx:299-341/385-415`、`IntTools_BeanFaceIntersector.cxx:178,200,223,583`、`BOPAlgo_WireSplitter_1.cxx:792,795`、`BOPAlgo_PaveFiller_6.cxx:3001-3002`），Circle/Ellipse/Parabola/Hyperbola 从 D1 臂放回 `Resolution` 臂；`ReparamCurve2d`/`PlaneKeepParam2d` 两个装饰器补 `resolution` 以免退回 trait default。门禁：`occt-geom2d` 72/0、`--lib` 1251/0、`--tests` 1317/0（`step_obj_gates` 4/4）、`export_data_obj` 15/15 且 v/f 逐位一致。**同批残留（未改）**：`inttools_range::find_nearest_valid_point` 仍缺 `BRepLib_1.cxx:69-83` 的 `aD1Mag` 局部奇点加速分支（原注释已标 `ponytail:`）。",
  },
  {
    item: "`GCPnts_AbscissaPoint::AdvCompute` 缺 `GCPnts_AbsComposite` 臂；两处 adv 入口的 `theEPSILON` 误传 `Resolution(Confusion)`",
    id: "A15",
    where: "geom/gcpnts.rs:compute_abs_composite_adv · abscissa_point_from_guess · abscissa_point_with_tolerance",
    evidence:
      "按 `GCPnts_AbscissaPoint.cxx:187-295` 补 `compute_abs_composite_adv`：逐步 `CPnts_AbscissaPoint::Length(..., theEPSILON)`、`anIndex == 0 && aDirection > 0` 前置块（`:202-229`，仅 adv 有）、`Init(C, U1, U2, theEPSILON)` + `AdvPerform`、循环内零长判据用 `Precision::PConfusion()`（`:232`）而前置块用 `theEPSILON`（`:206`）、尾段按 `aSign` 偏移并按 `FirstParameter`/`LastParameter` 钳制（`:273-293`）。同时订正 `.cxx:486`/`:540` 的 `theEPSILON` 语义：原先传 `c.resolution(Confusion)`，现传调用方 `tol`（`Approx_CurvlinFunc::GetUParameter` 传 `myTolLen`，`Approx_CurvlinFunc.cxx:560`）。`occt-geom` 143 测试通过、`step_obj_gates` 4/4。",
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
    state: "已闭环（2026-10-15，批次 3/T-129 + 批次 4/T-130）：余 19 处调用点全部换成 `is_plane_type` / `face_is_plane_type`（`GetType() == GeomAbs_Plane` 等价）——`brep_surface.rs:153`（`face_plane`）、`pcurve.rs:249`、`fillet_curved/surface_info.rs:139,220`、`fillet_curved/torus_blend.rs:603`、`gprop_analytic/analytic_props.rs:378,587`、`gprop_analytic/mass_properties.rs:17,61,325,443`、`brep_offset/shell_offset.rs:122`、`brep_faces.rs:203` 及测试夹具（`brep_builder_full/mod.rs:247`、`brep_sketch.rs:234`、`brep_builder_api.rs:372,405-406`、`fclass2d/tests.rs:160`、`brep_offset/tests.rs:132`、`fillet_curved/tests/helpers.rs:693`、`brep_surface.rs:288,304`）；随后 `is_planar`/`face_is_planar` 采样体整段删除，全仓已无该符号。同批把 `edge_face_kind.rs` 的曲线/曲面族（`curve_kind`/`circle_geometry`/`surface_kind`/`is_eq_distance`/`is_coplanar`/`is_radius`）换成 `GetType()`/`gp_*` 同等 downcast。批次 4（T-130）再把 `Resolution` 采样族换成 `Adaptor::Resolution` 解析闭式：`edge_face_kind::curve_resolution`、`inttools_range::curve_resolution`（3 处）、`algo_tools_range` 的 `curve_is_line`+`curve_is_freeform`+`curve_resolution`、`bean_face::resolution`、`int_tools_vertex_line::adaptor_resolution`、`wire_splitter_block::angle_2d`，并在 `occt-geom2d` 新增 `Curve2d::resolution`（`Geom2dAdaptor_Curve::Resolution` 等价）与两个装饰器（`ReparamCurve2d`/`PlaneKeepParam2d`）的 `resolution`。门禁：`occt-geom2d` 72/0、`--lib` 1251/0、`--tests` 1317/0、`step_obj_gates` 4/4（23 模型面积比 1.0000）、`export_data_obj` 15/15 且 v/f 逐位一致。**新登记的同族残留**：`inttools_range::find_nearest_valid_point` 缺 `BRepLib_1.cxx:69-83` 的 `aD1Mag` 局部奇点加速分支（原注释已标 `ponytail:`）。",
  },
  {
    item: "`BOPAlgo_WireSplitter` 的 `Angle2D` / `Tolerance2D` / `RefineAngle2D` 三臂（`_audit/bop-boolean.md` 发现 19）",
    id: "发现 19 / T-131",
    where: "topo/wire_splitter_block.rs（原 `uv_tolerance_2d:580` · `curve2d_radius:543` · `refine_angle_2d:771`）",
    state: "**已修（2026-10-15）**：(1) `uv_tolerance_2d` 原先手抄 `GeomAdaptor_Surface::U/VResolution`，回转(V)/拉伸(U) 两臂退回有限差分 `fd_uv_tolerance`（`t3/|dS/du|`），且缺 `RectangularTrimmedSurface` 下基（`cxx:423-430`）、缺 `OffsetSurface` 臂（`cxx:1882-1885`/`:1944-1947`）、缺锥面无界守卫（`cxx:1855-1858`）、锥半径用采样代 `S->VIso(V)->Radius()`（`cxx:1861-1866`）——现改转发唯一忠实件 `approx_same_parameter::{u,v}_resolution`，手抄 switch 与差分件删除。(2) `angle_2d` 的内联 `|D1 x D2|/|D1|^3`（守卫 `n2 <= PConfusion` 量纲不符）改 `GeomLProp_CLProps2d` + `is_tangent_defined()`/`curvature()`（`BOPAlgo_WireSplitter_1.cxx:797-807`）。(3) `refine_angle_2d` 的 `aTV` 由 3D 边参数改 pcurve 参数 `BRep_Tool::Parameter(aV,aE,aFace)`（`cxx:1060`）；UNPORTED 注记清除。门禁：`--lib` 1251/0、`step_obj_gates` 4/4、`export_data_obj` 15/15 且 v/f 逐位一致。",
  },
  {
    item: "`Analyzer::is_valid` 恒 true",
    id: "A14",
    where: "topo/validate.rs:150",
    state: "已声明非校验、不得作门禁",
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
    what: "SEMANTIC 表 LIVE 方向也过期",
    detail: "24 条登记的 drift 行其代码已闭环（topo-rest 1/4/6/11 + 18 的 classifier 半、bop 1-4/7-16-17、mesh 1-14/19）",
    gap: "❗未同步：`docs/goals/drift-and-gap-autofix/notes/T001-inventory.md` 逐条给了证据，本表尚未按方向订正（不是代码问题，是留档口径问题）",
    tone: "danger" as const,
  },
  {
    what: "未上板的自创件：`bbox_from_geometry.rs` 采样包围盒",
    detail: "64 点边包围盒 + 16×16 UV 面包围盒（`bbox_from_geometry.rs:46,56-67`）",
    gap: "❗新发现（T001）：无 `BRepBndLib`/`BndLib` 解析极值；被 `bop_builder.rs:86`、`bop_builder_dispatch.rs:109`、`bop_classify_occt.rs:78`、`builder_solid.rs:128` 使用，尚无任何审计文档登记",
    tone: "danger" as const,
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

// 看板门禁（读数日期 2026-10-15，覆盖 T-130 / T-131 已提交工作树）
const GATES = [
  { id: "G-check", name: "五 crate 编译", base: "exit 0", now: "exit 0", tone: "success" as const },
  { id: "G-lib", name: "topo 单测 --lib", base: "1276 / 0", now: "1251 / 0", tone: "success" as const },
  { id: "G-gates", name: "STEP->OBJ 门禁", base: "4 / 4", now: "4 / 4", tone: "success" as const },
  { id: "G-fuse", name: "fuse 闭环实体", base: "1 / 1", now: "1 / 1", tone: "success" as const },
  { id: "G-phase", name: "phase3-10 / 19 / 20", base: "61 / 0", now: "61 / 0", tone: "success" as const },
  { id: "G-boss", name: "boss 合并", base: "2 / 2", now: "2 / 2", tone: "success" as const },
  { id: "G-export", name: "导出总检", base: "23 / 23", now: "15 ok / 0 err", tone: "success" as const },
  { id: "G-iges", name: "IGES 结构自洽", base: "23 / 23", now: "23 / 23", tone: "success" as const },
  { id: "G-core", name: "core/math/geom/geom2d", base: "290·215·143·72", now: "290·215·143·72", tone: "success" as const },
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

      <H3>2.3 逐模型对拍（`data/occ/*`，2026-10-13 实测）</H3>
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

      <H3>2.4 状态 / 文档偏移（2026-10-13 复核：全部收口）</H3>
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
        九道门禁全绿，读数测量于 2026-10-13，覆盖当前工作树（`git diff HEAD` 为空，提交 `4a20e078`）。`G-lib` / `G-export` 的看板基线已在 r42 同步为 `1276 / 0` 与 `23 / 23`（见 2.4）；`G-fuse` / `G-phase` / `G-boss` / `G-core` 本轮未复跑（改动不触及其路径）。
      </Text>

      <Divider />

      <H2>三、结论与下一步</H2>
      <Row gap={16} wrap>
        <Stat value="13.8%" label="全库行数覆盖" />
        <Stat value="24.0%" label="CAD 内核行数覆盖" tone="info" />
        <Stat value="2 处" label="阻塞型语义残余" tone="warning" />
        <Stat value="0 项" label="文档/状态偏移待办" tone="success" />
      </Row>
      <Stack gap={6}>
        <Text>
          <Text weight="semibold">范围</Text>：全库 13.8%、CAD 内核约 24.0% 的行数覆盖；有效范围是 STEP→OBJ 共享管线，
          在此链上已做逐行对 `.cxx` 的深度移植（occt-topo 占 21.6 万行）。
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
        面级残余的定位受「可信的 OCCT 面级预言机」未钉住所限：实测同一条 `--facestats 2.45132 0.349066`
        在探针上是 `72683 / 78464`，而冻结参考是 `73544 / 79977`、端口是 `73516 / 79926` —— 探针用的是
        `STEPControl_Reader` 直读，与冻结参考的产出者（未记录）以及端口的不同源，故探针总量不可直接当日面级基线；
        要么先钉住与参考同源的读取+参数组合，要么只以冻结 OBJ 的聚合口径为规格。
        在此之前可推进的独立项是 2.2 里的自创残余：`heal_tolerance` 自创修形链（本轮已收口：`FixSmall`/`CheckSmall`/`FixSmallAreaWire`/`CopyReplaceVertices` 落地，`boolean_with_check` 容差重算删除）、`Geom_RectangularTrimmedSurface::SetTrim` 族与 `gcpnts` 的 `GCPnts_AbsComposite` adv 臂（两轮均本轮收口，见 2.1）均已闭环；`bop_curved` 网格级布尔（A5/T-79）已随引擎删除闭环，`brepfeat` 的圆柱夹具（A12 余项 / T-125）本轮已换 `BRepPrimAPI_MakeCylinder` 且体积精确（`box ∪ boss` 4.157080、`box ∪ 贯通 boss` 2.565487），现存自创项只剩 `revolve_profile_about` 的解析旋转体工具 —— 缺口经本轮实测不在 `BOPAlgo_BOP`，而在 `sweep_revolve` 未按 `BRepSweep_Revol`/`BRepPrim_OneAxis` 的 wire/seam/pcurve 约定建面（见 2.2 与 `_audit/_index.md` A12）。本轮第二个独立项是「用采样识别二次曲面代替 downcast」（`_audit/topo-rest.md` 发现 14）：偏移面分派已换精确 downcast —— `offset_face` 平面/球臂改 `gp_pln()`/`gp_sphere()`（见 2.1 首条），`brep_surface::sphere_center` 改 `GeomAdaptor_Surface::Sphere()` 等价并删除采样外接圆心族；残留的 `brep_surface::is_planar`/`face_is_planar` 调用点第二批 9 处也已换 `is_plane_type`（BOP 交叉 + EF 家族，逐条 `.cxx` 行号见 2.2 发现 14 余项），第三批（T-129）换完余 19 处并删除采样体，`edge_face_kind` 的曲线/曲面族亦换成 `GetType()` 同等 downcast ⇒ 发现 14 已整体收口；批次 4（T-130）把同族的 `Resolution` 采样估计（`edge_face_kind::curve_resolution`、`inttools_range` 3 处、`algo_tools_range` 的共线/自由形判据、`bean_face::resolution`、`int_tools_vertex_line::adaptor_resolution`、`wire_splitter_block::angle_2d`）全部换成 `Adaptor::Resolution`／`Curve2d::resolution` 解析闭式，并在 `occt-geom2d` 补上 `Curve2d::resolution`（`Geom2dAdaptor_Curve::Resolution` 等价，逐型 Line/Circle/Ellipse/Bezier/BSpline/Trimmed），`ReparamCurve2d`/`PlaneKeepParam2d` 两个装饰器补 `resolution` 以免退回 trait default。剩余同族残留只剩 `inttools_range::find_nearest_valid_point` 缺 `BRepLib_1.cxx:69-83` 的 `aD1Mag` 局部奇点加速分支（原注释已标 `ponytail:`）。
      </Callout>
    </Stack>
  );
}
