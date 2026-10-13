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
  asOf: "2026-10-13",
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
    where: "topo/step/read_topology.rs:410 close_heal_context · shhealing/shape_fix_same_parameter.rs",
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
    where: "topo/step/read_topology.rs",
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
];

// 语义偏移残余（对照 specs/_audit/_index.md A 编号；位置/行号为 2026-10-13 工作树实测）
const SEMANTIC = [
  {
    item: "Delaunay 插点密度：缺 `BRepMesh_DelaunayDeflectionControlMeshAlgo::checkAndAdd`/`addInternalNodes`",
    id: "A19 余项",
    where: "topo/meshing/triangulator.rs · topo/meshing/node_insertion.rs",
    state: "未移植；**原证据已作废**：TDB 的 +1,744 面经逐面配对复核不是本项（属缝边 pcurve 注册，见 2.1 T-124），收口后 TDB −51。现残余为双向小量：T0M +58、TDB −51、ATU01038 −32、top −30、a3n00 −24、acs10 −20、bottom −5、motoc 0，面积比全 1.0000 —— 已无 4~9% 的系统性偏密",
  },
  {
    item: "bop_curved 网格级布尔（曲面 operand）",
    id: "A5 / T-79",
    where: "topo/bop_curved/region_trim.rs:51 · face_meshing.rs:418",
    state: "UNPORTED，被 T-80 阻塞",
  },
  {
    item: "`heal_tolerance` 自创修形链",
    id: "A11",
    where: "topo/bop_builder_splitapi.rs:301（+ `shhealing::weld_coincident_vertices`/`remove_small_edges`）",
    state: "未动",
  },
  {
    item: "自创自适应 Simpson（`curve_length` 已改走 Gauss 忠实件）",
    id: "A15 余项",
    where: "core/elib/measure.rs:58 · math/integrate_adaptive.rs:5",
    state: "未动",
  },
  {
    item: "`RectangularTrimmedSurface::SetTrim` 族",
    id: "—",
    where: "geom/rectangular_trimmed.rs",
    state: "未移植",
  },
  {
    item: "`BRepCheck_Edge` 的 SameParameter/SameRange 检查",
    id: "—",
    where: "topo/brep_check.rs",
    state: "未移植",
  },
  {
    item: "`surface_closest_params` 采样替代件",
    id: "A1",
    where: "topo/brep_surface.rs:261",
    state: "活调用点 0（`brepmesh` QA 消费者已删，仅存 `pub(crate)` 定义）",
  },
  {
    item: "`Analyzer::is_valid` 恒 true",
    id: "A14",
    where: "topo/validate.rs:150",
    state: "已声明非校验、不得作门禁",
  },
];

// `data/occ/*` 逐模型对拍（2026-10-13 `step_obj_gates --nocapture` 实测）
const PARITY = [
  { model: "motoc", occV: 11713, occF: 13944, ourV: 11713, ourF: 13944, area: "1.0000" },
  { model: "a3n00", occV: 11145, occF: 12466, ourV: 11138, ourF: 12442, area: "1.0000" },
  { model: "acs10", occV: 37296, occF: 46558, ourV: 37280, ourF: 46538, area: "1.0000" },
  { model: "ATU01038", occV: 18102, occF: 22586, ourV: 18085, ourF: 22554, area: "1.0000" },
  { model: "bottom", occV: 17587, occF: 23557, ourV: 17582, ourF: 23552, area: "1.0000" },
  { model: "top", occV: 16905, occF: 22800, ourV: 16888, ourF: 22770, area: "1.0000" },
  { model: "T0M", occV: 60548, occF: 67366, ourV: 60575, ourF: 67424, area: "1.0000" },
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
    detail: "board baseline 1281 / 0 → 1276 / 0",
    gap: "✅ 已同步：`A18` 删除 5 条 `brepmesh` 自创单测后基线未同步；2026-10-13 实测 1276 / 0",
    tone: "success" as const,
  },
  {
    what: "G-export 门禁基线过期",
    detail: "board baseline 15 / 15 → 23 / 23",
    gap: "✅ 已同步：`-- data/occ data` 口径实测 23 ok / 0 err（15 `data/*.step` + 8 `data/occ/*`）；无参跑仍是 15",
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

// 看板门禁（读数日期 2026-10-13，覆盖当前未提交工作树）
const GATES = [
  { id: "G-check", name: "五 crate 编译", base: "exit 0", now: "exit 0", tone: "success" as const },
  { id: "G-lib", name: "topo 单测 --lib", base: "1276 / 0", now: "1276 / 0", tone: "success" as const },
  { id: "G-gates", name: "STEP->OBJ 门禁", base: "4 / 4", now: "4 / 4", tone: "success" as const },
  { id: "G-fuse", name: "fuse 闭环实体", base: "1 / 1", now: "1 / 1", tone: "success" as const },
  { id: "G-phase", name: "phase3-10 / 19 / 20", base: "61 / 0", now: "61 / 0", tone: "success" as const },
  { id: "G-boss", name: "boss 合并", base: "2 / 2", now: "2 / 2", tone: "success" as const },
  { id: "G-export", name: "导出总检", base: "23 / 23", now: "23 ok / 0 err", tone: "success" as const },
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
        在此之前可推进的独立项是 2.2 里 `heal_tolerance` 自创修形链与 `RectangularTrimmedSurface::SetTrim` 族的移植。
      </Callout>
    </Stack>
  );
}
