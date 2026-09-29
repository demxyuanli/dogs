/** @canvas
 * title: 项目看板 · OCCT→Rust 移植
 * description: 支撑「下一批做哪条 / 门禁能否提交 / 哪些红是已知非缺陷」；任务行为 T-xx，门禁单独一段，指标全部现算
 * icon: board
 */
import {
  H1,
  H2,
  Text,
  Code,
  Stack,
  Row,
  Grid,
  Divider,
  Card,
  CardBody,
  CardHeader,
  Callout,
  Stat,
  Table,
  BarChart,
  Pill,
  Button,
  CollapsibleSection,
  useMemo,
  useCanvasState,
  useCanvasOverlay,
  useCanvasAction,
} from "dsh/canvas";

// ---- 5 值状态模型（intake 决议）------------------------------------------
type Status = "todo" | "doing" | "blocked" | "done" | "waived";
type Tone = "neutral" | "info" | "success" | "warning" | "danger";

// 快照日与阈值都在 DATA 里，派生量一律现算，禁止手抄。
// status 五种取值都至少出现一次（overlay.set 的类型前提），owner 三种取值也都出现一次。
export const DATA = {
  meta: {
    project: "OCCT 8.0.0 → Rust 移植（STEP→OBJ 几何/网格对齐）",
    northStar: "把 STEP→OBJ 几何/网格管线对齐 OCCT 8.0.0（源码树 D:\\source\\OCCT-src @ V8_0_0）",
    repo: "D:\\source\\repos\\dogs",
    head: "08c4edf2",
    snapshot: "2026-09-28",
    staleAfterDays: 7,
    wipLimit: 2,
    gate: "红线：无同等 OCCT 分支就标 UNPORTED（写 文件+行号）；不加自造谓词/阈值；不为对齐新写测试或改断言。",
    archive: "done/waived 移出当前窗口，归档到 git 与 specs/_a3n00_gap_analysis.md。",
    measured: "门禁数字必须当轮实测（本轮 2026-09-28 @ 08c4edf2）；从旧文档手抄会立刻过期——T-88 就是这样被误记为红。",
    goal: {
      id: "goal-cb46d692-98fc-4805-845e-8e205504595c",
      rounds: 50,
      objective: "把本板 10 个未完成任务逐条推进到 done（或带实测证据 waive），验收以各条 accept 为准；遵守红线与五道门禁。",
      waves: [
        { w: 1, ids: "T-93", why: "P0：面 #4952 的 8→22 展开缺失（specs/_a3n00_gap_analysis.md §9.215）" },
        { w: 2, ids: "T-69", why: "同族 wire 过度拆分；T-93 修好大概率一起推进" },
        { w: 3, ids: "T-54", why: "等 T-69 解锁：回退路径耳切 → 约束 Delaunay（§9.209）" },
        { w: 4, ids: "T-41", why: "boolean_dispatch 改派 bop_builder2 + 摘 bop_curved（§9.206）" },
        { w: 5, ids: "T-67, T-51", why: "occt-geom 两支补移植（Extrema / gcpnts）" },
        { w: 6, ids: "T-25, T-28", why: "架构项：侧表进 TShape / HVertex 合并" },
        { w: 7, ids: "T-11, T-29", why: "清理项：编译警告 / 过期规格" },
      ],
    },
  },
  lanes: ["bop", "mesh", "port-gap", "arch", "hygiene"],
  tasks: [
    {
      id: "T-93",
      lane: "port-gap",
      title: "reader 缺带 Context 的 ShapeFix_Face::Perform（a3n00 阀体 F113 无网格）",
      status: "doing",
      priority: "P0",
      owner: "agent",
      est: 4,
      actual: 4,
      start: "2026-09-28",
      updated: "2026-09-28",
      done: "",
      blockedBy: "",
      next: "**根因已定位到「reader 另建顶点」，最后一步待做**（§9.281–§9.286）。链条：a3n00 面积比 0.8627 → 亏损集中在特定区域（0 覆盖，镜像成对）→ 该区域由**单个面 F171**（BSpline/4 wires/整周期 v）覆盖 → 其第 4 条 wire 在端口**少 1 条边**（6 vs GT 7）→ 端口该 wire 用了一个**不在 STEP 里的端点 P**，而 P 恰是边 `#7405` 曲线的**中点** → 运行时探针证明顶点 `#4956`(A) 在解析时与面处理后**都是正确坐标**（未被改写）⇒ **P 是 reader 在 `resolve_vertex` 之外另建的顶点**。**最后一步**：搜 reader 里所有 `make_vertex` / 新建顶点并拷点的路径（首候选 `shape_fix_compose_shell/split_wire.rs:45`），找出在 F171 该 loop 上「按中点断边/建点」的那一处；有 OCCT 同等分支就替换，无则标 UNPORTED 并删除",
      accept: "F113 出网格；--test step_obj_gates 5/5 不劣化；a3n00 面积比向 1 收敛",
      evidence: "已修一处：check_notched_edges 取端点用原始 edge_vertices，而 OCCT 用朝向感知的 sae.LastVertex(E1)/FirstVertex(E2)（ShapeAnalysis_Wire.cxx:1885-1898）⇒ 对 REVERSED 边取反，notch 永不命中。改后 F113 wire[0] 8→6、notched=true；--lib 1281/0、step_obj_gates 5/5。根因定案（本轮，OCCT 侧插桩）：bbox 唯一吻合的面 = explorer 25，OCCT 原始 2 wire（22/6 边）/ 处理后 1 wire，FixMissingSeam 返回 0；端口 4 wire [6,14,1,1] 且 fix_missing_seam 返回 true 才进 ComposeShell ⇒ 分歧在 reader 的 wire 组装（同一条 wire 被拆成 4 条）。工具：specs/occt_probe/_dbg/ + build_dbg.bat + run_dbg.bat + wires_probe 的 findf 模式",
      write: "crates/occt-topo/src/shhealing/shape_fix_face.rs",
      ref: "ShapeFix_Face.cxx:365-498 · ShapeFix_Shape.cxx:202",
      commit: "",
    },
    {
      id: "T-82",
      lane: "bop",
      title: "BOPAlgo_BOP::BuildSolid 曲面布尔精确化",
      status: "done",
      priority: "P0",
      owner: "agent",
      est: 3,
      actual: 3,
      start: "2026-09-26",
      updated: "2026-09-28",
      done: "2026-09-28",
      blockedBy: "",
      next: "—",
      accept: "--lib 不劣化；box∪cyl 单实体、体积≈8.5；phase10 / boss 不劣化",
      evidence: "HEAD 08c4edf2 实测 box∪cyl FUSE = 1 solid / 8 面 / volume 8.502655，逐面面积 0.502655·2.513274·3.497345·4×5 与 GT（occt_probe --fuse：8.50265 / 8 面 / 同一组面积）逐一相同；卡片的 8.737 与 loops=0 是更早快照，未再复现。AC 未覆盖的其他曲面布尔组合未复核",
      write: "crates/occt-topo/src/bop_bop.rs",
      ref: "BOPAlgo_BOP.cxx:583-711 · BOPAlgo_Builder_3.cxx::BuildSplitSolids",
      commit: "",
    },
    {
      id: "T-88",
      lane: "bop",
      title: "phase10 curved_face_fillet_sphere_plane",
      status: "done",
      priority: "P0",
      owner: "agent",
      est: 1,
      actual: 1,
      start: "2026-09-26",
      updated: "2026-09-28",
      done: "2026-09-28",
      blockedBy: "",
      next: "—",
      accept: "--test phase10_integration 8/8",
      evidence: "HEAD 08c4edf2 实测 phase10 8/8（curved_face_fillet_sphere_plane ok）；旧板记的 7/8 是 2026-09-21 快照，已过期；原诊断（内部盘未丢 + 圆边重复）未再复现",
      write: "crates/occt-topo/src/bop_bop.rs",
      ref: "BOPAlgo_BuilderSolid · BOPTools_AlgoTools",
      commit: "",
    },
    {
      id: "T-41",
      lane: "bop",
      title: "摘除 bop_curved 的体素/网格布尔（OCCT 无对应）",
      status: "done",
      priority: "P1",
      owner: "agent",
      est: 2,
      actual: 3,
      start: "2026-09-20",
      updated: "2026-09-28",
      done: "2026-09-28",
      blockedBy: "",
      next: "已结项。三处改派 + repair/report 改派 + 夹具重设计（按决策 (a)，§9.267）",
      accept: "摘除后 --lib 与门禁不劣化",
      evidence: "A5：采样计票在 OCCT 里没有对应物（OCCT 用 BRepClass3d_SolidClassifier / IntTools_Context::ComputeState；面区间来自 BOPAlgo_BuilderFace 的 pcurve）。调用图：curved_boolean_full ← boolean_dispatch:227（本轮已改派）/ bop_builder_repair:239,405 / bop_builder_report:388；curved_boolean ← brepfeat/features:320；voxel_boolean ← brepfeat/features:330-335。本轮实测（改派后）：--lib 1281/0、step_obj_gates 5/5（380s）、phase10/19/3/4 = 8/5/4/9 ⇒ 无劣化；classify_face_general 本身也是 8×8 采样 + point_in_solid_curved，不是 BRepClass3d_SolidClassifier（§9.206 已记）",
      write: "crates/occt-topo/src/bop_curved/region_trim.rs",
      ref: "—",
      commit: "",
    },
    {
      id: "T-69",
      lane: "mesh",
      title: "T0M 未网格面（6）+ wire 结构与 GT 分歧",
      status: "doing",
      priority: "P1",
      owner: "agent",
      est: 3,
      actual: 3,
      start: "2026-09-24",
      updated: "2026-09-28",
      done: "",
      blockedBy: "",
      next: "**矛盾已解开（§9.268）**：T-93 (a) 移除补偿后重测，resolver 探针与模型侧完全一致 —— STEPFACE 全场 wires<bounds=0、Σbounds=Σwires=2121 且等于模型 TOTAL wires=2121。⇒ 此前第 36 轮那个「190 个面 wires<bounds」（即「端口自己丢 wire」）是**当时那个 `fix_missing_seam` 补偿的副作用**：它即使丢弃 Shell(5) 也会改写共享 GeometryRegistry，使解析时探针看到的形状与最终模型分叉（§9.228 已记同一现象）。**同时暴露关键张力**：移除补偿后 T0M 的 wire 结构明显变差（wires 1931→2121、多 wire 面 120→309），而 GT 是 1921/106 ⇒ T-69 的 accept「WIREHIST 对齐 GT」反而更远。下一步（需决策）：T-69 与 T-93 (a) 处于直接张力中 —— 要么恢复补偿（T0M 好、F113 不出网格），要么保留现状（F113 出网格、T0M 变差），要么把补偿做成**无副作用**（只对目标面生效、不回写共享 registry）——第三条路正是 T-93 (b)，可同时满足两者",
      accept: "T0M 未网格面 6→≤1（GT 自己也有 1 个）；WIREHIST 对齐 GT；T0M 面类计数 + 门禁",
      evidence: "本轮实测（HEAD 08c4edf2）：端口 T0M faces=1772 / stats=1766 ⇒ 未网格 6 个（不是卡片写的 162），三类：F351（Cylinder, 10 wires, dv=inf）、F405-408（BSpline, 1 wire/2 edges, du=0.43 dv=2π）、F1758（Cone, 1 wire/1 edge）。GT 对照：occt_probe --mesh 0.1 ⇒ faces=1778 / triangles=75743 / triangles=0 的面 1 个（FACE 1327）；occt_probe --wires ⇒ faces=1778 wires=1921 faces_with_coincident_closed_edges=0 wires_dup_pc_any=117 WIREHIST 1:1672 2:92 3:8 5:1 6:1 7:2 8:2",
      write: "crates/occt-topo/src/meshing/model_healer.rs",
      ref: "BRepMesh_ModelHealer.cxx:491-504 · BRepMesh_Delaun.cxx:2263-2274",
      commit: "",
    },
    {
      id: "T-54",
      lane: "mesh",
      title: "平面耳切 → 约束 Delaunay",
      status: "blocked",
      priority: "P2",
      owner: "agent",
      est: 4,
      actual: 1,
      start: "2026-09-20",
      updated: "2026-09-21",
      done: "",
      blockedBy: "T-69",
      next: "仍 blocked → T-69（本轮已核实前提与范围）：要改的是回退路径。等 T-69 把那 6 个「Delaunay 空产出」的面定案后，先判断它们是否本应走 wireframe_face_triangulation 回退（:554）、为何没走或走了也空 —— 再决定是把回退换成约束 Delaunay（BRepMesh_DelaunayBaseMeshAlgo/BRepMesh_Delaun），还是直接删回退让主线报错暴露",
      accept: "网格门禁",
      evidence: "wireframe 的耳切（ear_clip:156）+ 桥洞（bridge_holes:219）+ planar_polygon_triangulate:269 在 OCCT 里没有对应控制流。本轮实测调用点全是回退：discret_root.rs:290（build_shape_mesh_wireframe 整形状回退）/ :554（wireframe_face_triangulation 单面回退，Delaunay 失败时）/ :1376+:1382（discretize_face 单面回退）；主路径是 discret_root 的 Delaunay 管线",
      write: "crates/occt-topo/src/wireframe.rs",
      ref: "BRepMesh_DelaunayBaseMeshAlgo.cxx · BRepMesh_Delaun.cxx",
      commit: "",
    },
    {
      id: "T-67",
      lane: "port-gap",
      title: "Extrema_ExtPExtS / ExtPRevS 两臂",
      status: "done",
      priority: "P2",
      owner: "agent",
      est: 3,
      actual: 3,
      start: "2026-09-26",
      updated: "2026-09-28",
      done: "2026-09-28",
      blockedBy: "",
      next: "—",
      accept: "occt-geom --lib 143/0 + 门禁不劣化",
      evidence: "两臂已移植并分派（point_surface_extrema.rs:60-61/:88-91 类型、:450 perform_ext_ps、:485 perform_rev_ps，与 Extrema_ExtPS.cxx:292-343 逐段对应；两个实现文件无 UNPORTED）。本轮摘掉 A15 替身：point_surface_extrema_box 不再调用非 OCCT 的 point_surface_newton_all_box / fallback_point_surface，窗口为空时返回「无解」（v2=None），由调用方 pcurve_full::surface_projector::value_of_uv 走已移植的 uv_from_iso（OCCT ShapeAnalysis_Surface.cxx:1449-1459）。文档（numeric_extrema.rs 两处、extrema_surf/mod.rs）已同步。验收实测：occt-geom --lib 143/0、occt-topo --lib 1281/0、cargo test --test step_obj_gates 5/5（384s）",
      write: "crates/occt-geom/src/extrema_surf",
      ref: "Extrema_ExtPS.cxx:292-343",
      commit: "",
    },
    {
      id: "T-51",
      lane: "port-gap",
      title: "gcpnts 两处 UNPORTED",
      status: "waived",
      priority: "P3",
      owner: "agent",
      est: 2,
      actual: 1,
      start: "2026-09-26",
      updated: "2026-09-28",
      done: "2026-09-28",
      blockedBy: "",
      next: "—",
      accept: "occt-geom / occt-topo --lib",
      evidence: "本轮核实：可表达的那部分**已实现** —— CPnts_AbscissaPoint 的反解由 math_FunctionRoot 驱动，tolerance 重载 init_tol（gcpnts.rs:296-300）与 AdvPerform（adv_perform :407-422，对应 CPnts_AbscissaPoint.cxx:436-473，含 Resolution/10）都在。余下两条 Compute arm（GCPnts_LengthParametrized :87-90、GCPnts_AbsComposite :96-158）**无法忠实表达**：它们读 GeomAdaptor_Curve 的 GetType/NbIntervals/Intervals，而 OCCT 的 GeomAdaptor_Curve::load 对 Geom_TrimmedCurve 会下沉到 basis 并保留 basis 参数区间（GeomAdaptor_Curve.cxx:239-255），端口的 GeomTrimmedCurve 却把参数重映射到 [0,1]（trimmed.rs:44-70 的 d0/d1/d2/d3 链式法则）。要忠实移植必须先改 TrimmedCurve 的参数表示（架构级、影响所有 trimmed 曲线消费方），且无消费方触发这两条 arm ⇒ 豁免。实测：occt-geom --lib 143/0、occt-topo --lib 1281/0",
      write: "crates/occt-geom/src/gcpnts.rs",
      ref: "CPnts_AbscissaPoint.cxx:32-37 :436-474 :87-90 :96-158",
      commit: "",
    },
    {
      id: "T-25",
      lane: "arch",
      title: "GeometryRegistry 侧表 → 几何进 TShape",
      status: "done",
      priority: "P2",
      owner: "agent",
      est: 6,
      actual: 6,
      start: "2026-09-20",
      updated: "2026-09-28",
      done: "2026-09-28",
      blockedBy: "",
      next: "—",
      accept: "--lib + 门禁不劣化",
      evidence: "本轮核实：迁移已到设计文档所说的**端阶段**。① TShape 自带几何槽与访问器（tshape.rs:33-45 的 edge_pcurves/edge_core/vertex_core/face_core，:94-123 的 getter/mut）；② GeometryRegistry 只剩 ids + face_surfaces + surface_by_ptr，三个几何表已删除（tgeometry.rs:163-181 注释逐字写着重明：「the three geometry maps are gone … all live on their own TShape now」），35 处读写全部转发到 TShape 槽（如 :211-219 set_vertex、:226 读 vertex_core、:248/:255/:273 edge_core、:358-379 edge_pcurves、:399-445 face_core）；③ 设计文档允许「停在任一中间步，保留 global() 作兼容外壳」（第 5 步未做，属文档明示的可停点）。验收实测：occt-topo --lib 1281/0、step_obj_gates 5/5（399s）",
      write: "specs/_design_architecture_t25_t28.md",
      ref: "specs/_design_architecture_t25_t28.md",
      commit: "",
    },
    {
      id: "T-28",
      lane: "arch",
      title: "ImpPrm HVertex 合并 + intana/intpatch 重叠合并",
      status: "done",
      priority: "P2",
      owner: "agent",
      est: 8,
      actual: 8,
      start: "2026-09-20",
      updated: "2026-09-28",
      done: "2026-09-28",
      blockedBy: "",
      next: "已结项。标题两项均完成：① ImpPrm 的 HVertex 合并（步 1–3）；② intana/intpatch 的 closed form 重叠合并（步 4，五个副本清零）。设计文档的第 5 步「合并两个 surface×surface 分派器」经实测**判为有害、不做**（见 §9.250）：把 intpatch.rs:138 的通用回退换成统一入口 PatchIntersection 后，general_bspline_sphere_uses_fallback 立即失败（曲线离球面 ≈0.67），而原 tracer 正确 ⇒ 保留双分派器是当前唯一正确形态。若日后要重开步 5，前提是先把 PatchIntersection 对「一般 B 样条 × 球」这类 case 的产出修正到与 tracer 等价",
      accept: "--lib + 门禁不劣化",
      evidence: "设计已出（specs/_design_architecture_t25_t28.md §T-28，6 步）。本轮逐项核实（校正了文档的两处过期描述）：步 1/2/3 已完成（见 next 的行号）；步 4 未做 —— intpatch_analytic.rs 的 5 个 closed form 仍在（我先前按不带 intersect_ 前缀的名字 grep 得到假阴性，已更正）；步 5 未做，但主路径已走 PatchIntersection（IntSS），intpatch.rs:138 只剩回退角色。验收实测（当前树）：occt-topo --lib 1281/0、step_obj_gates 5/5",
      write: "specs/_design_architecture_t25_t28.md",
      ref: "IntPatch_ImpPrmIntersection.cxx:221-469",
      commit: "",
    },
    {
      id: "T-11",
      lane: "hygiene",
      title: "occt-topo 编译警告清理",
      status: "waived",
      priority: "P3",
      owner: "agent",
      est: 1,
      actual: 4,
      start: "2026-09-26",
      updated: "2026-09-28",
      done: "",
      blockedBy: "",
      next: "已 waived（§9.249）：四批共清 11 处真实 unused import，剩余项经实验判定为不可安全清除（见 evidence）。若日后要动「never used」那批死代码，应另立卡（那属于删 API 面，不是清 import）",
      accept: "无行为变化",
      evidence: "【waive 依据（2026-09-28）】四批清掉 **11 处真实 unused import**（line.rs GpAx1、wire_fix.rs GpDir、bop_builder_planar_weld.rs edges_of、bop_split_seam.rs ShapeIterator、bop_builder_planar_geom.rs BoolOp、int_face_face.rs DEFAULT_WINDOW、geom2d_int/ginter.rs GpPnt2d、shape_fix_compose_shell/split_by_grid.rs BRepTool、dispatch_wires.rs Curve2d、tgeometry.rs FaceGeomCore、wire_segment.rs edge_vertices），每批后 occt-topo --lib 1281/0、check --tests 0 错。**剩余 ~60 条同类警告经实验判定不可安全清除**：它们多是「符号只被同文件 #[cfg(test)] 代码使用」的假阳性（lib 目标不编译 cfg(test)，诊断却报 unused），本轮实测三例 —— bop_builder.rs 的 std::sync::Arc / GpAx3 / GeomPlane，各自删除后 test 目标立即 cannot find type ⇒ 已回退；另有 brep_gprop_full/gauss.rs 的 use super::*（glob）与 bopds.rs 的 shapes_of（跨文件 use super::* 消费）两类陷阱。两条自动化路径亦已实测失败（见下）。⇒ 人工逐符号是唯一安全路径而剩余收益见底，故 waive；「never used」408 条属删 API 面，需另立卡。\n【原始记录】实测基线（本轮）：occt-topo --lib 警告 506 条 = never used 408 / unused import 46 / other 37 / never read 11 / unused variable 3。两条自动化路径都已实测失败并回退：① cargo fix --lib（及 --lib --tests）会删掉只被 #[cfg(test)] 模块使用的 import（lib target 不带 cfg(test) 编译）⇒ 报 109 个构建错误；② 用 cargo check --message-format json 的 span 自行施加替换也不可靠（unused_imports 的 span 不总是「整条语句」，出现 pub use ; 与 {A, , C} 之类语法错）",
      write: "crates/occt-topo/src",
      ref: "—",
      commit: "",
    },
    {
      id: "T-29",
      lane: "hygiene",
      title: "过期规格刷新",
      status: "done",
      priority: "P3",
      owner: "unassigned",
      est: 1,
      actual: 1,
      start: "2026-09-21",
      updated: "2026-09-28",
      done: "2026-09-28",
      blockedBy: "",
      next: "—",
      accept: "—",
      evidence: "两个目标文件已带「仅历史（T-29 登记）」横幅（specs/_coverage.md:1、specs/_brepmesh_align_review.md:3）；本轮修正其指向已删除的 specs/_board.md ⇒ specs/board.canvas.tsx，并同步 4 处规范性悬空引用（_a3n00_gap_analysis.md:6、_design_architecture_t25_t28.md:3/:117、specs/occt_probe/README.md:33、crates/occt-core/src/bspl/poles.rs:26）；specs/_audit/ 下 8 处为历史审查记录，保留原样。occt-core check 通过、occt-topo --lib 1281/0，无行为变化",
      write: "specs/_coverage.md",
      ref: "—",
      commit: "",
    },
    {
      id: "T-23",
      lane: "hygiene",
      title: "ATU01038 顶点/面密度 −1.8%",
      status: "waived",
      priority: "P3",
      owner: "human",
      est: 1,
      actual: 0,
      start: "2026-09-21",
      updated: "2026-09-21",
      done: "",
      blockedBy: "",
      next: "不改，作为已知记录",
      accept: "—",
      evidence: "UV-grid 与 deflection-adaptive 的采样差异；parity 已明确不断言密度",
      write: "crates/occt-topo/tests/common/mod.rs",
      ref: "—",
      commit: "",
    },
    {
      id: "T-03",
      lane: "bop",
      title: "boss_single_disc_base_merges_one_solid（门禁 1/2）",
      status: "waived",
      priority: "P3",
      owner: "human",
      est: 1,
      actual: 1,
      start: "2026-09-19",
      updated: "2026-09-21",
      done: "",
      blockedBy: "",
      next: "维持 1/2 作为已知，不再投入",
      accept: "—",
      evidence: "夹具期望的「1 solid」在 OCCT 的 BOPAlgo 下也不成立 ⇒ 夹具不合法，非引擎缺口",
      write: "crates/occt-topo/tests/bop_builder2_boss.rs",
      ref: "BOPAlgo_BOP.cxx:583-711",
      commit: "",
    },
    {
      id: "T-92",
      lane: "bop",
      title: "ComposeShell / FixMissingSeam 逐行迁移",
      status: "done",
      priority: "P1",
      owner: "agent",
      est: 4,
      actual: 4,
      start: "2026-09-25",
      updated: "2026-09-28",
      done: "2026-09-28",
      blockedBy: "",
      next: "—",
      accept: "--lib 1281/0、step_obj_parity 14/14、export 逐位一致",
      evidence: "a3n00 wire 直方图 2-wire 面 24→10；mesh 11675/13168",
      write: "crates/occt-topo/src/shape_fix_compose_shell",
      ref: "ShapeFix_ComposeShell.cxx · ShapeFix_Face.cxx:1722-2330",
      commit: "daf7ee4d",
    },
    {
      id: "T-80",
      lane: "bop",
      title: "BOP 根因三处（AddWireEdge 朝向 / MakeVertex 容差 / REVERSED 边取 PCurve2）",
      status: "done",
      priority: "P1",
      owner: "agent",
      est: 3,
      actual: 3,
      start: "2026-09-25",
      updated: "2026-09-26",
      done: "2026-09-26",
      blockedBy: "",
      next: "—",
      accept: "box∪cyl FUSE = 8 面（GT 8.50265）",
      evidence: "GT 探针 --fuse：volume 8.50265 / 8 faces",
      write: "crates/occt-topo/src/bop_bop.rs",
      ref: "BRepPrim_Builder.cxx:184-192 · BRep_Tool.cxx:301-315",
      commit: "da61a30",
    },
    {
      id: "T-05",
      lane: "hygiene",
      title: "Offset 解析体积对齐 GT",
      status: "done",
      priority: "P1",
      owner: "agent",
      est: 2,
      actual: 2,
      start: "2026-09-20",
      updated: "2026-09-21",
      done: "2026-09-21",
      blockedBy: "",
      next: "—",
      accept: "step_geometry_parity 3/3（断言按 OCCT 规格改对 GT 实测值）",
      evidence: "2610.501436 vs GT 2610.501440（相对 1.5e-9）",
      write: "crates/occt-geom/src/bspline_curve.rs",
      ref: "BSplCLib::Reverse · BRepGProp_Face::Bounds",
      commit: "59fefad",
    },
    {
      id: "T-04",
      lane: "bop",
      title: "phase19 断言按 OCCT 口径订正",
      status: "done",
      priority: "P1",
      owner: "agent",
      est: 1,
      actual: 1,
      start: "2026-09-20",
      updated: "2026-09-21",
      done: "2026-09-21",
      blockedBy: "",
      next: "—",
      accept: "--test phase19_integration 5/5",
      evidence: "GT 探针 --ds：不相交 68 / 相交 80；去掉自造 2*56 阈值",
      write: "crates/occt-topo/tests/phase19_integration.rs",
      ref: "BOPDS_DS.cxx::Init/Append",
      commit: "—",
    },
    {
      id: "T-01",
      lane: "mesh",
      title: "groove_cuts_cylinder 转绿",
      status: "done",
      priority: "P2",
      owner: "agent",
      est: 2,
      actual: 2,
      start: "2026-09-19",
      updated: "2026-09-20",
      done: "2026-09-20",
      blockedBy: "",
      next: "—",
      accept: "occt-topo --lib 该用例转绿",
      evidence: "BRepLib_MakeWire 移植 + OrientClosedSolid；removed 0.414 → 0.909",
      write: "crates/occt-topo/src/brep_lib_make_wire.rs",
      ref: "BRepLib_MakeWire.cxx:123-453",
      commit: "113b4d1",
    },
  ],
  gates: [
    { id: "G-check", name: "五 crate 编译", cmd: "cargo check --all-targets", result: "exit 0", status: "green", known: false, owner: "", note: "core · math · geom · geom2d · topo" },
    { id: "G-lib", name: "topo 单测", cmd: "--lib", result: "1281 / 0", status: "green", known: false, owner: "", note: "T-01 / T-86 已闭" },
    { id: "G-gates", name: "STEP→OBJ 门禁", cmd: "--test step_obj_gates", result: "5 / 5", status: "green", known: false, owner: "", note: "23 模型一次读入 660.58s（原三份 2168s）" },
    { id: "G-phase", name: "phase3-10 / 19 / 20", cmd: "--test phaseN_integration", result: "全绿", status: "green", known: false, owner: "", note: "4·9·7·5·5·5·8·8·5·5（phase10 8/8，T-88 已闭）" },
    { id: "G-boss", name: "boss 合并", cmd: "--test bop_builder2_boss", result: "1 / 2", status: "red", known: true, owner: "T-03", note: "夹具不合法（定案 b），非引擎缺口" },
    { id: "G-export", name: "导出总检", cmd: "--example export_data_obj", result: "15 / 15", status: "green", known: false, owner: "", note: "内部回归哨兵（旧板记 16/16 为过期）；保真度见 specs/_occt_mesh_gt.md" },
    { id: "G-iges", name: "IGES 结构自洽", cmd: "--example iges_check -- <18 模型>", result: "18 / 18", status: "green", known: false, owner: "", note: "15 个 data/*.step + data/occ/{ATU01038,bottom,top}.step；全部 unreferenced=1（只剩根）" },
    { id: "G-core", name: "core / math / geom / geom2d", cmd: "--lib", result: "290·215·143·72", status: "green", known: false, owner: "", note: "逐项与基线相同" },
  ],
} as const;

const STATUS_TONE: Record<Status, Tone> = {
  todo: "neutral",
  doing: "info",
  blocked: "danger",
  done: "success",
  waived: "warning",
};

const STATUS_LABEL: Record<Status, string> = {
  todo: "todo",
  doing: "doing",
  blocked: "blocked",
  done: "done",
  waived: "waived",
};

const OPEN_STATUS: readonly Status[] = ["todo", "doing", "blocked"];

// 纯整数民用日历 → 日内序数（Howard Hinnant 算法），不引 Date，避免非纯调用。
const toDays = (s: string) => {
  const parts = s.split("-").map(Number);
  const y = parts[0];
  const m = parts[1];
  const d = parts[2];
  const yy = m <= 2 ? y - 1 : y;
  const era = Math.floor(yy / 400);
  const yoe = yy - era * 400;
  const mp = m + (m > 2 ? -3 : 9);
  const doy = Math.floor((153 * mp + 2) / 5) + d - 1;
  const doe = yoe * 365 + Math.floor(yoe / 4) - Math.floor(yoe / 100) + doy;
  return era * 146097 + doe;
};

export default function ProjectBoard() {
  const dispatch = useCanvasAction();
  const [filter, setFilter] = useCanvasState<string>("filter", "open");
  const [activeId, setActiveId] = useCanvasState<string>("active", "T-93");
  const overlay = useCanvasOverlay("tasks", DATA.tasks);
  const tasks = overlay.items;

  const snapshot = DATA.meta.snapshot;
  const snapDay = toDays(snapshot);

  const open = useMemo(
    () => tasks.filter((t) => OPEN_STATUS.includes(t.status as Status)),
    [tasks],
  );
  const closed = useMemo(
    () => tasks.filter((t) => !OPEN_STATUS.includes(t.status as Status)),
    [tasks],
  );
  const visible = useMemo(() => {
    if (filter === "all") return tasks;
    if (filter === "open") return open;
    if (filter === "closed") return closed;
    if (filter === "blocked") return tasks.filter((t) => t.status === "blocked");
    return tasks.filter((t) => t.lane === filter);
  }, [tasks, open, closed, filter]);

  const active = tasks.find((t) => t.id === activeId) ?? visible[0] ?? tasks[0];

  // 派生 1：加权进度（权重 = est；completion：done=1，其余 min(1, actual/est)）
  const completion = (t: (typeof DATA.tasks)[number]) => {
    if (t.status === "done") return 1;
    if (t.est <= 0) return 0;
    return Math.min(1, t.actual / t.est);
  };
  const weightSum = open.reduce((s, t) => s + Math.max(1, t.est), 0);
  const weighted = weightSum === 0 ? 0 : open.reduce((s, t) => s + Math.max(1, t.est) * completion(t), 0) / weightSum;

  // 派生 2：WIP
  const doing = tasks.filter((t) => t.status === "doing");
  const wipOver = doing.length > DATA.meta.wipLimit;

  // 派生 3：陈旧
  const staleOf = (t: (typeof DATA.tasks)[number]) => snapDay - toDays(t.updated);
  const stale = open.filter((t) => staleOf(t) > DATA.meta.staleAfterDays);

  // 派生 4：超估算
  const overrun = open.filter((t) => t.actual > t.est);

  // 派生 5：分组成度
  const byLane = DATA.lanes.map((lane) => {
    const l = tasks.filter((t) => t.lane === lane);
    const lo = l.filter((t) => OPEN_STATUS.includes(t.status as Status));
    const w = lo.reduce((s, t) => s + Math.max(1, t.est), 0);
    const p = w === 0 ? 0 : lo.reduce((s, t) => s + Math.max(1, t.est) * completion(t), 0) / w;
    return { lane, open: lo.length, blocked: lo.filter((t) => t.status === "blocked").length, done: l.length - lo.length, pct: Math.round(p * 100) };
  });

  // 派生 6：下一步 3 条（优先级 → 陈旧优先）
  const rank: Record<string, number> = { P0: 0, P1: 1, P2: 2, P3: 3 };
  const next3 = open
    .slice()
    .sort((a, b) => rank[a.priority] - rank[b.priority] || staleOf(b) - staleOf(a))
    .slice(0, 3);

  // 派生 7：风险
  const risks: { id: string; sev: Tone; what: string }[] = [];
  tasks.forEach((t) => {
    if (t.status === "blocked" && !t.blockedBy) risks.push({ id: t.id, sev: "danger", what: "blocked 但没写等待对象" });
    if (t.status === "doing" && t.owner === "unassigned") risks.push({ id: t.id, sev: "warning", what: "doing 但无人认领" });
    if (OPEN_STATUS.includes(t.status as Status) && staleOf(t) > DATA.meta.staleAfterDays) risks.push({ id: t.id, sev: "warning", what: "陈旧 " + staleOf(t) + " 天（> " + DATA.meta.staleAfterDays + "）" });
    if (t.actual > t.est && t.status !== "done") risks.push({ id: t.id, sev: "warning", what: "超估算 " + t.actual + "/" + t.est });
  });
  if (wipOver) risks.push({ id: "WIP", sev: "warning", what: "doing " + doing.length + " > 上限 " + DATA.meta.wipLimit });
  DATA.gates.forEach((g) => {
    if (g.status === "red" && !g.known) risks.push({ id: g.id, sev: "danger", what: "门禁红且未标已知非缺陷" });
  });

  const blockers = DATA.gates.filter((g) => g.status === "red" && !g.known);

  const startTurn = (t: (typeof DATA.tasks)[number]) =>
    dispatch({
      type: "startTurn",
      prompt:
        "看板决定 —— 处理 " + t.id + "：" + t.title +
        "\n状态：" + t.status + " · 优先级：" + t.priority +
        "\n下一步：" + t.next +
        "\n验收：" + t.accept +
        "\n证据：" + t.evidence +
        "\nOCCT 参考：" + t.ref +
        "\n遵守看板红线（无同等分支标 UNPORTED，不加自造阈值，不为对齐新写测试）。",
    });

  const handoff = () =>
    dispatch({
      type: "startTurn",
      prompt:
        "看板决定 —— 按优先级推进：" +
        next3.map((t) => "\n- " + t.id + "（" + t.priority + "）：" + t.title + " → " + t.next).join("") +
        "\n未澄清红门禁：" + (blockers.length === 0 ? "无" : blockers.map((g) => g.id).join(", ")) +
        "\n请先给出执行计划，再动手。",
    });

  return (
    <Stack gap={24}>
      <Stack gap={6}>
        <H1>项目看板 · OCCT→Rust 移植</H1>
        <Text tone="secondary">{DATA.meta.northStar}</Text>
        <Text size="small" tone="tertiary">
          快照 {DATA.meta.snapshot} · HEAD <Code>{DATA.meta.head}</Code> · 仓库 <Code>{DATA.meta.repo}</Code>
        </Text>
      </Stack>

      <Callout tone="info" title={"Goal " + DATA.meta.goal.id + " · " + DATA.meta.goal.rounds + " 轮上限 · armed"}>
        <Stack gap={3}>
          <Text>{DATA.meta.goal.objective}</Text>
          {DATA.meta.goal.waves.map((wv) => (
            <Text key={wv.w} size="small" tone="secondary">
              {"波次 " + wv.w + "：" + wv.ids + " —— " + wv.why}
            </Text>
          ))}
        </Stack>
      </Callout>

      <Grid columns={4} gap={16}>
        <Stat value={Math.round(weighted * 100) + "%"} label="加权进度（在办）" tone="info" hint="Σ(est×completion)/Σest，只在办条目" />
        <Stat value={doing.length + " / " + DATA.meta.wipLimit} label="WIP" tone={wipOver ? "danger" : "success"} hint="doing 数 / 上限" />
        <Stat value={open.filter((t) => t.status === "blocked").length} label="blocked" tone="danger" hint="等待前置的条目" />
        <Stat value={stale.length} label={"陈旧 > " + DATA.meta.staleAfterDays + " 天"} tone={stale.length ? "warning" : "success"} hint="快照日 − 更新日" />
      </Grid>

      <Callout tone={blockers.length === 0 ? "success" : "danger"} title={blockers.length === 0 ? "可以提交：红项都是已知非缺陷" : "不可提交：" + blockers.length + " 条未澄清红门禁"}>
        {blockers.length === 0
          ? "门禁红项均已登记为已知非缺陷（" + DATA.gates.filter((g) => g.status === "red").map((g) => g.id).join(", ") + "），不阻塞提交。"
          : "未澄清红门禁：" + blockers.map((g) => g.id + "（" + g.name + "）").join("、") + "。"}
      </Callout>

      <Grid columns="minmax(0, 1fr) minmax(0, 1fr)" gap={20} align="start">
        <Stack gap={10}>
          <H2>下一步（派生 Top 3）</H2>
          <Text size="small" tone="tertiary">按 优先级 → 陈旧天数 排序，全部由 DATA 现算。</Text>
          {next3.map((t) => (
            <Card key={t.id}>
              <CardHeader trailing={<Pill size="sm" tone={STATUS_TONE[t.status]}>{STATUS_LABEL[t.status]}</Pill>}>
                {t.id} · {t.priority}
              </CardHeader>
              <CardBody>
                <Stack gap={6}>
                  <Text weight="semibold">{t.title}</Text>
                  <Text size="small">{t.next}</Text>
                  <Row gap={8} wrap>
                    <Button size="sm" onClick={() => setActiveId(t.id)}>看详情</Button>
                    <Button size="sm" variant="primary" onClick={() => startTurn(t)}>交给 agent</Button>
                  </Row>
                </Stack>
              </CardBody>
            </Card>
          ))}
          <Button variant="primary" onClick={handoff}>把 Top 3 一起交给 agent</Button>
        </Stack>

        <Stack gap={10}>
          <H2>风险清单（派生）</H2>
          <Table
            headers={["对象", "级别", "风险"]}
            columnAlign={["left", "left", "left"]}
            striped
            emptyText="无派生风险"
            rows={risks.map((r) => [
              <Code>{r.id}</Code>,
              <Pill size="sm" tone={r.sev}>{r.sev === "danger" ? "high" : "warn"}</Pill>,
              r.what,
            ])}
            rowTone={risks.map((r) => r.sev)}
          />
          <H2>分组成度（派生）</H2>
          <BarChart
            categories={byLane.map((b) => b.lane)}
            series={[{ name: "加权进度 %", data: byLane.map((b) => b.pct), tone: "info" }]}
            beginAtZero
            height={160}
          />
          <Table
            headers={["lane", "在办", "blocked", "已完成", "进度%"]}
            columnAlign={["left", "right", "right", "right", "right"]}
            striped
            rows={byLane.map((b) => [b.lane, b.open, b.blocked, b.done, b.pct + "%"])}
          />
        </Stack>
      </Grid>

      <H2>门禁基线</H2>
      <Table
        headers={["门禁", "命令", "结果", "红项归属", "说明"]}
        columnAlign={["left", "left", "right", "left", "left"]}
        striped
        stickyHeader
        rows={DATA.gates.map((g) => [
          g.name,
          <Code>{g.cmd}</Code>,
          <Pill size="sm" tone={g.status === "green" ? "success" : g.known ? "warning" : "danger"}>{g.result}</Pill>,
          g.owner ? <Code>{g.owner}</Code> : <Text size="small" tone="tertiary">—</Text>,
          <Text size="small" tone="tertiary">{g.note}</Text>,
        ])}
        rowTone={DATA.gates.map((g) => (g.status === "green" ? "success" : g.known ? "warning" : "danger"))}
      />

      <H2>任务</H2>
      <Row gap={8} wrap>
        <Pill active={filter === "open"} onClick={() => setFilter("open")}>Open {open.length}</Pill>
        <Pill active={filter === "blocked"} onClick={() => setFilter("blocked")}>Blocked {open.filter((t) => t.status === "blocked").length}</Pill>
        <Pill active={filter === "closed"} onClick={() => setFilter("closed")}>Closed {closed.length}</Pill>
        <Pill active={filter === "all"} onClick={() => setFilter("all")}>All {tasks.length}</Pill>
        <Divider />
        {DATA.lanes.map((lane) => (
          <Pill key={lane} active={filter === lane} onClick={() => setFilter(lane)}>{lane}</Pill>
        ))}
      </Row>

      <Grid columns="minmax(0, 1.25fr) minmax(0, 0.75fr)" gap={20} align="start">
        <Table
          headers={["ID", "lane", "标题", "状态", "进度", "更新", "阻塞在等"]}
          columnAlign={["left", "left", "left", "left", "right", "right", "left"]}
          striped
          stickyHeader
          emptyText="没有条目"
          onRowClick={(index) => {
            const row = visible[index];
            if (row) setActiveId(row.id);
          }}
          rows={visible.map((t) => [
            <Code>{t.id}</Code>,
            t.lane,
            t.title,
            <Pill size="sm" tone={STATUS_TONE[t.status]}>{STATUS_LABEL[t.status]}</Pill>,
            t.actual + "/" + t.est,
            t.updated,
            t.blockedBy ? <Code>{t.blockedBy}</Code> : <Text size="small" tone="tertiary">—</Text>,
          ])}
          rowTone={visible.map((t) => STATUS_TONE[t.status])}
        />

        {active ? (
          <Card>
            <CardHeader trailing={<Pill size="sm" tone={STATUS_TONE[active.status]}>{STATUS_LABEL[active.status]}</Pill>}>
              {active.id} · {active.priority}
            </CardHeader>
            <CardBody>
              <Stack gap={10}>
                <Text weight="semibold">{active.title}</Text>
                <Text size="small">lane <Code>{active.lane}</Code> · owner <Code>{active.owner}</Code> · 估算 {active.est} / 实际 {active.actual}</Text>
                <Text size="small">开始 {active.start} · 更新 {active.updated}{active.done ? " · 完成 " + active.done : ""} · 陈旧 {staleOf(active)} 天</Text>
                {active.blockedBy ? <Text size="small" tone="danger">阻塞：在等 <Code>{active.blockedBy}</Code></Text> : null}
                <Divider />
                <Text size="small"><Text weight="semibold">下一步 </Text>{active.next}</Text>
                <Text size="small"><Text weight="semibold">验收 </Text>{active.accept}</Text>
                <Text size="small"><Text weight="semibold">证据 </Text>{active.evidence}</Text>
                <Text size="small">改动 <Code>{active.write}</Code></Text>
                <Text size="small">OCCT <Code>{active.ref}</Code></Text>
                {active.commit ? <Text size="small">提交 <Code>{active.commit}</Code></Text> : null}
                <Divider />
                <Row gap={8} wrap>
                  <Button variant="primary" onClick={() => startTurn(active)}>交给 agent</Button>
                  <Button onClick={() => overlay.set(active.id, { status: "doing" })}>标记 doing</Button>
                  <Button onClick={() => overlay.set(active.id, { status: "done" })}>标记 done</Button>
                  <Button onClick={() => overlay.set(active.id, { status: "blocked" })}>标记 blocked</Button>
                  <Button onClick={() => overlay.set(active.id, { status: "waived" })}>豁免</Button>
                  <Button variant="ghost" onClick={() => overlay.set(active.id, { owner: "human" })}>认领</Button>
                  <Button variant="ghost" onClick={() => overlay.clear(active.id)}>还原</Button>
                  <Button variant="ghost" onClick={() => dispatch({ type: "openFile", path: active.write })}>打开文件</Button>
                </Row>
              </Stack>
            </CardBody>
          </Card>
        ) : null}
      </Grid>

      <CollapsibleSection title="已关闭（最近完成 / 已豁免）" count={closed.length}>
        <Table
          headers={["ID", "lane", "标题", "状态", "提交", "证据"]}
          columnAlign={["left", "left", "left", "left", "left", "left"]}
          striped
          rows={closed.map((t) => [
            <Code>{t.id}</Code>,
            t.lane,
            t.title,
            <Pill size="sm" tone={STATUS_TONE[t.status]}>{STATUS_LABEL[t.status]}</Pill>,
            <Code>{t.commit || "—"}</Code>,
            <Text size="small" tone="tertiary">{t.evidence}</Text>,
          ])}
          rowTone={closed.map((t) => STATUS_TONE[t.status])}
        />
      </CollapsibleSection>

      <CollapsibleSection title="口径与节奏（intake 决议）" count={7}>
        <Stack gap={6}>
          <Text size="small">· 决定：下一批做哪条 + 门禁能否提交 + 哪些红是已知非缺陷。</Text>
          <Text size="small">· 行单位：任务 T-xx（lane = 模块/族）；门禁在「门禁基线」段。</Text>
          <Text size="small">· 状态：todo / doing / blocked / done / waived；blocked 必填 blockedBy。done 或 waived 后移出当前窗口。</Text>
          <Text size="small">· 字段：状态 · 优先级 · owner · est/actual · start/updated/done · blockedBy · next · accept · evidence · write · ref · commit。</Text>
          <Text size="small">· 派生（全部现算）：加权进度 · WIP 超限 · 陈旧 · 超估算 · 分组成度 · 下一步 3 条 · 风险清单。</Text>
          <Text size="small">· 人的动作：标记 doing/done/blocked、豁免、认领、还原、打开文件、交给 agent —— 全部走 overlay，落 .canvas/board.state.json。</Text>
          <Text size="small">· 节奏：agent 每批改 DATA；人只动状态/认领/决定；归档到 git 与 specs/_a3n00_gap_analysis.md。</Text>
          <Text size="small" tone="tertiary">{DATA.meta.gate}</Text>
          <Text size="small" tone="tertiary">{DATA.meta.measured}</Text>
          <Text size="small" tone="tertiary">{DATA.meta.archive}</Text>
        </Stack>
      </CollapsibleSection>

      <Text size="small" tone="tertiary">
        每个数字都从 DATA 现算；明细表的「证据」列指向 文件:行 / 测试名 / 提交号。结论性描述用引用，不粘贴段落。
      </Text>
    </Stack>
  );
}
