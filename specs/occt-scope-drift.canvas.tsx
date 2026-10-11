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
//   基线  : specs/board.canvas.tsx (r40, asOf 2026-10-07) · specs/_audit/_index.md (§9 -> 2026-10-06)
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

// 状态/文档偏移：文档落后于代码与 HEAD
const DOCDRIFT = [
  { what: "看板基准日", detail: "board asOf = 2026-10-07 · r40", gap: "HEAD 已到 2026-10-10 c8060dd1（约 +39 提交）；tasks 段 31 完成 / 3 取消 / 0 在办" },
  { what: "10-10 未上板批次", detail: "brep_class3d：SClassifier::Perform 段循环 / OtherSegment / FindAPointInTheFace", gap: "活动流最新仅到 a89（2026-10-08）" },
  { what: "_audit §9.3 过期断言", detail: "把 geom2d/curve_ops.rs 采样族列为未修", gap: "该采样已于 2026-09-26（ok267451e8）删除，改走 Extrema_ExtCC2d / IntAna2d" },
  { what: "未提交工作树", detail: "93 个已跟踪文件修改", gap: "+6,082 / −27,470；删量主要来自 _a3n00_gap_analysis.md 精简（约 2.6 万行 → 摘要）" },
];

// 看板门禁（读数日期 2026-10-06 … 10-08，未覆盖当前未提交工作树）
const GATES = [
  { id: "G-check", name: "五 crate 编译", base: "exit 0", now: "exit 0", tone: "success" as const },
  { id: "G-lib", name: "topo 单测 --lib", base: "1281 / 0", now: "1281 / 0", tone: "success" as const },
  { id: "G-gates", name: "STEP->OBJ 门禁", base: "4 / 4", now: "4 / 4", tone: "success" as const },
  { id: "G-fuse", name: "fuse 闭环实体", base: "1 / 1", now: "1 / 1", tone: "success" as const },
  { id: "G-phase", name: "phase3-10 / 19 / 20", base: "全绿 61", now: "全绿 61", tone: "success" as const },
  { id: "G-boss", name: "boss 合并", base: "1 / 2", now: "1 / 2", tone: "danger" as const },
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
          所有读数取 <Text weight="semibold">{MEASURE.asOf}</Text>；门禁基线引用 <Text weight="semibold">specs/board.canvas.tsx</Text>（r40）与 <Text weight="semibold">specs/_audit/_index.md</Text>（§9 → 2026-10-06）。
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

      <H3>2.2 状态 / 文档偏移（文档落后于代码与 HEAD）</H3>
      <Table
        headers={["偏差", "文档现状", "实际"]}
        rows={DOCDRIFT.map((d) => [d.what, d.detail, d.gap])}
        columnAlign={["left", "left", "left"]}
        rowTone={["warning", "warning", "danger", "info"]}
        striped
      />
      <Callout tone="warning" title="本轮最紧的隐患">
        看板与审查文档落后于 HEAD（2026-10-10）与当前未提交工作树：`_audit/_index.md` §9.3 仍把
        `geom2d/curve_ops.rs` 采样族列为未修，而该采样早在 2026-09-26 已删除并改走忠实 `Extrema_ExtCC2d` / `IntAna2d`。
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
        唯一 red（G-boss 1/2）已在板上定案为「夹具不合法，非引擎缺口」的已接受红。所有读数测量于 2026-10-06…10-08，未覆盖当前未提交工作树。
      </Text>

      <Divider />

      <H2>三、结论与下一步</H2>
      <Row gap={16} wrap>
        <Stat value="13.9%" label="全库行数覆盖" />
        <Stat value="24.2%" label="CAD 内核行数覆盖" tone="info" />
        <Stat value="2 处" label="阻塞型语义残余" tone="warning" />
        <Stat value="4 项" label="文档/状态偏移" tone="danger" />
      </Row>
      <Stack gap={6}>
        <Text>
          <Text weight="semibold">范围</Text>：全库 13.9%、CAD 内核约 24% 的行数覆盖；有效范围是 STEP→OBJ 共享管线，在此链上已做逐行对 `.cxx` 的深度移植（occt-topo 占 21.9 万行）。
        </Text>
        <Text>
          <Text weight="semibold">偏移</Text>：语义自创已大幅收敛，残余集中在网格（A13/A18 UV 栅格）与曲面布尔（T-79/T-80）两处，且二者互相阻塞；文档 / 看板偏移是当前更紧的隐患。
        </Text>
      </Stack>
      <Callout tone="info" title="下一步（建议入口）">
        先对当前未提交工作树重建基线 —— `cargo check --all-targets` + `occt-topo --lib` + `--test step_obj_gates`，
        再据实同步 `specs/board.canvas.tsx` 与 `specs/_audit/_index.md` §9.3（含 `curve_ops` 过期断言订正）。
      </Callout>
    </Stack>
  );
}
