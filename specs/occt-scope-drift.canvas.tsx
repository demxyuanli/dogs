/** @canvas
 * title: OCCT→Rust 移植范围与偏移审查
 * description: 全库范围量测（OCCT vs Rust）、语义/文档偏移残余清单、门禁现状与下一步
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
// 量测口径（2026-10-11 实测）
//   OCCT 侧: D:\source\OCCT-src\src  .hxx/.cxx/.lxx/.gxx
//   Rust 侧: crates/{occt-core,occt-math,occt-geom,occt-geom2d,occt-topo}/src  .rs
//   基线  : specs/board.canvas.tsx (r41, asOf 2026-10-11) · specs/_audit/_index.md (§9 -> 2026-10-11)
// ---------------------------------------------------------------------------
const MEASURE = {
  asOf: "2026-10-11",
  occt: { files: 13803, lines: 2439487 },
  rust: { files: 1080, lines: 337891 },
  cadKernel: 1397840, // Foundation 252505 + ModelingData 347970 + ModelingAlgorithms 797365
  unported: 169,
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
  { crate: "occt-topo", files: 626, lines: 218852 },
  { crate: "occt-core", files: 219, lines: 44720 },
  { crate: "occt-geom", files: 121, lines: 41381 },
  { crate: "occt-geom2d", files: 59, lines: 18644 },
  { crate: "occt-math", files: 55, lines: 14294 },
];

// 语义偏移残余（对照 specs/_audit/_index.md A 编号；行号为当前工作树实测）
const SEMANTIC = [
  { item: "面网格走 UV 栅格 + 耳切/桥洞", id: "A13/A18/A23", where: "topo/wireframe.rs:156,219,403", state: "UNPORTED，被 T-69 阻塞" },
  { item: "bop_curved 网格级布尔", id: "A5 / T-79", where: "bop_curved/region_trim.rs:44 · face_meshing.rs:410", state: "UNPORTED，被 T-80 阻塞" },
  { item: "heal_tolerance 自创修形链", id: "A11", where: "topo/bop_builder_splitapi.rs", state: "未动" },
  { item: "gcpnts 自适应 Simpson", id: "A15", where: "geom/gcpnts.rs", state: "未动" },
  { item: "RectangularTrimmedSurface::SetTrim 族", id: "—", where: "geom/rectangular_trimmed.rs", state: "未移植" },
  { item: "BRepCheck_Edge SameParameter/SameRange", id: "—", where: "topo/brep_check.rs", state: "未移植" },
  { item: "surface_closest_params 采样替代", id: "A1", where: "topo/brep_surface.rs:244", state: "活调用点降至 1（brepmesh QA）" },
  { item: "Analyzer::is_valid 恒 true", id: "A14", where: "topo/validate.rs:150", state: "已声明非校验、不得作门禁" },
];

// 状态/文档偏移：2026-10-11 复核定案 —— 前三条为已订正的历史偏移，后两条仍开
const DOCDRIFT = [
  { what: "看板基准日", detail: "board asOf = 2026-10-07 · r40", gap: "已订正：r41 · asOf 2026-10-11（c8d58992）", tone: "info" as const },
  { what: "G-boss 基线过期", detail: "board G-boss = 1 / 2（red，已接受）", gap: "已订正：2026-10-11 实测 2 / 2；4244ab19 修好夹具后红消解", tone: "info" as const },
  { what: "_audit §9.3 过期断言", detail: "把 geom2d/curve_ops.rs 采样族列为未修", gap: "已订正：该采样 2026-09-26（ok267451e8）已删除，改走 Extrema_ExtCC2d / IntAna2d", tone: "info" as const },
  { what: "10-10 未上板批次", detail: "brep_class3d：SClassifier::Perform 段循环 / OtherSegment / FindAPointInTheFace", gap: "a90 现覆盖基线重建；该批次的实现内容仍未上板", tone: "warning" as const },
  { what: "未提交工作树", detail: "106 个已跟踪文件修改", gap: "+6,567 / −27,578（topo 50 · geom 27 · core 22 · geom2d 5 · specs 2）；删量主要来自 _a3n00_gap_analysis.md 精简", tone: "danger" as const },
];

// 看板门禁（读数日期 2026-10-11，已覆盖当前未提交工作树；G-boss 基线由 1 / 2 订正为 2 / 2）
const GATES = [
  { id: "G-check", name: "五 crate 编译", base: "exit 0", now: "exit 0", tone: "success" as const },
  { id: "G-lib", name: "topo 单测 --lib", base: "1281 / 0", now: "1281 / 0", tone: "success" as const },
  { id: "G-gates", name: "STEP->OBJ 门禁", base: "4 / 4", now: "4 / 4", tone: "success" as const },
  { id: "G-fuse", name: "fuse 闭环实体", base: "1 / 1", now: "1 / 1", tone: "success" as const },
  { id: "G-phase", name: "phase3-10 / 19 / 20", base: "全绿 61", now: "全绿 61", tone: "success" as const },
  { id: "G-boss", name: "boss 合并", base: "2 / 2", now: "2 / 2", tone: "success" as const },
  { id: "G-export", name: "导出总检", base: "15 / 15", now: "15 / 15", tone: "success" as const },
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
          所有读数取 <Text weight="semibold">{MEASURE.asOf}</Text>；门禁基线引用 <Text weight="semibold">specs/board.canvas.tsx</Text>（r41）与 <Text weight="semibold">specs/_audit/_index.md</Text>（§9 → 2026-10-11）。
        </Text>
      </Stack>

      <Grid columns={4} gap={12}>
        <Stat value={pct(MEASURE.rust.lines, MEASURE.occt.lines)} label="全库行数比（CAD 全量）" />
        <Stat value={pct(MEASURE.rust.lines, MEASURE.cadKernel)} label="CAD 内核行数比" tone="info" />
        <Stat value={MEASURE.rust.lines.toLocaleString("en-US")} label="Rust 行数 · 1,080 文件" />
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
        来源：D:\source\OCCT-src\src（8 个模块组，.hxx/.cxx/.lxx/.gxx）· 2026-10-11 实测 · 单位：源码行。
      </Text>

      <H3>Rust 侧分布</H3>
      <Table
        headers={["crate", "文件", "行数"]}
        rows={CRATES.map((c) => [c.crate, String(c.files), c.lines.toLocaleString("en-US")])}
        columnAlign={["left", "right", "right"]}
      />

      <Callout tone="info" title="范围口径的要点">
        实际在办范围 ≠ 全库移植：看板目标已收敛为「STEP→OBJ 几何/网格管线对齐 OCCT」，
        即在 CAD 内核里只对 STEP 读/写 → Geom/BRep → BOPAlgo/BRepMesh/ShapeFix 这条共享链做逐行忠实化；
        Visualization / Draw / OCAF / Deprecated 是有意 lite 替代，不计入门禁。
      </Callout>

      <Divider />

      <H2>二、偏移</H2>

      <H3>2.1 语义偏移（相对 OCCT 的忠实度残余）</H3>
      <Row gap={8} align="center" wrap>
        <Pill size="md">自创已收敛</Pill>
        <Text tone="secondary">全仓共 <Text weight="semibold">{MEASURE.unported} 处 UNPORTED 标记</Text>；下表为仍可核实到的残余项。</Text>
      </Row>
      <Table
        headers={["残余项", "审计编号", "位置（当前工作树）", "状态"]}
        rows={SEMANTIC.map((s) => [s.item, s.id, s.where, s.state])}
        columnAlign={["left", "left", "left", "left"]}
        rowTone={SEMANTIC.map((s) =>
          s.state.startsWith("UNPORTED") ? "danger" : s.state.startsWith("未") ? "warning" : "info",
        )}
        striped
      />

      <H3>2.2 状态 / 文档偏移（2026-10-11 复核：3 条已订正，2 条仍开）</H3>
      <Table
        headers={["偏差", "文档现状", "实际"]}
        rows={DOCDRIFT.map((d) => [d.what, d.detail, d.gap])}
        columnAlign={["left", "left", "left"]}
        rowTone={DOCDRIFT.map((d) => d.tone)}
        striped
      />
      <Callout tone="info" title="文档偏移已订正 3 条（2026-10-11）">
        `board.canvas.tsx` 与 `_audit/_index.md §9.3` 的过期断言已按实测改定（r40 → r41；curve_ops 采样族订正；
        G-boss 由过期读数 1 / 2 订正为 2 / 2，看板唯一 red 清零）。剩余偏移只在未提交工作树与 10-10 未上板批次。
      </Callout>

      <H3>2.3 门禁现状</H3>
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
        原唯一 red（G-boss 1/2）已消解：`4244ab19` 修好夹具后 2026-10-11 实测 2/2，板上基线同步订正。所有读数测量于 2026-10-11，已覆盖当前未提交工作树。
      </Text>

      <Divider />

      <H2>三、结论与下一步</H2>
      <Row gap={16} wrap>
        <Stat value="13.9%" label="全库行数覆盖" />
        <Stat value="24.2%" label="CAD 内核行数覆盖" tone="info" />
        <Stat value="2 处" label="阻塞型语义残余" tone="warning" />
        <Stat value="5 项（3 已订正）" label="文档/状态偏移" tone="warning" />
      </Row>
      <Stack gap={6}>
        <Text>
          <Text weight="semibold">范围</Text>：全库 13.9%、CAD 内核约 24% 的行数覆盖；有效范围是 STEP→OBJ 共享管线，在此链上已做逐行对 `.cxx` 的深度移植（occt-topo 占 21.9 万行）。
        </Text>
        <Text>
          <Text weight="semibold">偏移</Text>：语义自创已大幅收敛，残余集中在网格（A13/A18 UV 栅格）与曲面布尔（T-79/T-80）两处，且二者互相阻塞；文档 / 看板偏移本轮已订正 3 条（看板基准日、G-boss 基线、§9.3 断言），仅剩未提交工作树与 10-10 未上板批次。
        </Text>
      </Stack>
      <Callout tone="success" title="本轮已闭环（2026-10-11）">
        基线重建完成：`cargo check --all-targets`（修 1 处 E0063 后五 crate 0 error）+ `occt-topo --lib` 1281 / 0
        + 九道门禁全绿，并据实同步 `specs/board.canvas.tsx`（r41）与 `specs/_audit/_index.md` §9.3（`curve_ops` 过期断言订正）。
      </Callout>
    </Stack>
  );
}
