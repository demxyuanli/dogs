/** @canvas
 * title: 项目看板 · OCCT→Rust 移植
 * description: 支撑「下一批做哪条 / 门禁能否提交 / 哪些红是已知非缺陷」；任务行 T-xx，门禁单独一段，指标全部现算
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
  TodoList,
  Progress,
  KeyValue,
  Timeline,
  Pill,
  Button,
  CollapsibleSection,
  useMemo,
  useCanvasState,
  useCanvasOverlay,
  useCanvasAction,
} from "dsh/canvas";

// ---- 状态模型（canvas 0.2.0 的 5 值枚举；旧板的 todo/doing/done/waived 已折算）----
type Status = "pending" | "in_progress" | "blocked" | "completed" | "cancelled";
type Tone = "neutral" | "info" | "success" | "warning" | "danger";

// 每个条目带的是「跟踪所需的字段」，不是只有标题与状态：
//   status / progress / estimate / actual     -> 到哪一步、还要多少、已花多少
//   startedAt / updatedAt / completedAt       -> 停了多久、多久没动过
//   blocker / next / acceptance / evidence    -> 卡在哪、下一步、怎样算完成、证据
//   write / ref / note / dependsOn            -> 改动位置、OCCT 参考、旁支问题、前置
// 只内联「当前窗口」：在办 + 最近完成 + 由旁支问题立的新卡。更早的历史进「已归档」折叠区。
// 版式（软阈值 1500 行 / 128 KB）：一条记录内联短字段、长散文字段各占一行；不要把一条记录摊回十几行
// —— 体量增长只应来自，也只需来自散文，而长段分析要写成 §编号引用（specs/_a3n00_gap_analysis.md），不要粘贴。
export const DATA = {
  goal: "把 STEP→OBJ 几何/网格管线对齐 OCCT 8.0.0（源码树 D:\\source\\OCCT-src @ V8_0_0）",
  asOf: "2026-10-06",
  revision: "r38",
  wipLimit: 2,
  staleDays: 7,
  lanes: ["bop", "mesh", "port-gap", "arch", "hygiene"],
  constraints: [
        { rule: "无同等 OCCT 分支就不许自造：缺口只能标 UNPORTED 并写清 OCCT 文件+行号",
       because: "网格/布尔管线是共享的，凭某一例的形状反推规则会带偏已对齐的其它模型",
      
      violation: "已对齐的 23 模型基线连带劣化，且失配原因无法归因到 OCCT 控制流",
    },
        { rule: "禁止为单个 STEP / 特征 / 特例调参（谓词、面积比、长度滤边、体积门）",
       because: "这些规则在 OCCT 里不存在，属于把特例补丁塞进共享管线",
      
      violation: "对拍变成自证：数字能凑上，换一个模型即失效，且无法与 OCCT 逐行对读",
    },
        { rule: "不为对齐新写单元测试，也不改既有断言；门禁数字必须当轮实测",
       because: "断言是基线，不是目标函数；从旧文档手抄会立刻过期（T-88 曾被误记为红）",
      
      violation: "基线失去参照意义，红线与门禁都变成不可复核的自我声明",
    },
        { rule: "`dependsOn` 只能指向**板上存在的卡**；缺口一旦成为前置，必须先立卡再引用",
      
      because: "依赖不在板上就永远不会自己解除 —— 这不是「等」而是死锁，而且看板的派生检查只能看到「未完成」，看不到「不存在」（D20 的真实事故：P0 的 T-93 曾挂在只存在于散文里的 T-59 上）",
      
      violation: "被依赖的卡静默挂死、无人可认领；风险清单不报警，下一轮接手者会以为「等它完成」是进度",
    },
        { rule: "报数量必须同时报**测量口径**（调用序列 / 配对方法 / 快照日）",
       because: "同一参数下不同调用序列给出不同 no_stat（1757/15 vs 1765/7）；两侧面数不同时序号不构成身份（1778 vs 1772）",
      
      violation: "数字不可比，据此下的丢面/回归结论是假结论（本板已实际发生两次，见 D13/D16）",
    },
        { rule: "跨侧对拍必须按**几何键**（bbox 六坐标）配对，不得按面序号",
      because: "a3n00 两侧面序号只有 3/226 对应（OCCT 读入整形会重建面序），而两侧面数相等（226=226）会伪装成同序",
      violation: "配到别的面上并得出不存在的倍数（§9.365.2 的 2–15 vs 50–114 倍就是这样来的）",
    },
  ],
  decisions: [
        { id: "D1", at: "2026-09-28", chose: "T-93 (b)：把 pcurve 归一化补偿做成只在目标面生效、不回写共享 GeometryRegistry",
       rejected: ["恢复全局 fix_missing_seam 补偿", "保留现状不动（F113 不出网格）"],
       why: "两条老路互斥：恢复补偿则 F113 又不出网格，保留现状则 T0M wire 结构更远离 GT；无副作用版本可同时满足两者",
      
      ref: "specs/_a3n00_gap_analysis.md §9.268 / §9.269",
    },
        { id: "D2", at: "2026-09-28", chose: "移除 reader 侧 ShapeFix_Face::FixMissingSeam 补偿（T-93 (a) 落地）",
       rejected: ["继续用补偿换取 T0M 的 wire 数"],
      why: "补偿即使丢弃 Shell(5) 也会改写共享 GeometryRegistry，使解析期探针与最终模型分叉，取证链被污染",
      
      ref: "specs/_a3n00_gap_analysis.md §9.266 / §9.268",
    },
        { id: "D3", at: "2026-09-28", chose: "boolean_dispatch 改派忠实路径；摘除 bop_curved 体素/网格布尔",
       rejected: ["保留采样计票 classify_face_general 作为曲面分类"],
       why: "体素/采样计票在 OCCT 无对应物（OCCT 走 BRepClass3d_SolidClassifier / IntTools_Context::ComputeState）",
      
      ref: "crates/occt-topo/src/bop_curved/region_trim.rs · §9.206",
    },
        { id: "D4", at: "2026-09-28", chose: "夹具按合法布尔重设计（boss_single_disc_base_merges_one_solid）",
       rejected: ["为夹具期望的 1 solid 改引擎"], why: "该夹具期望在 OCCT 的 BOPAlgo 下同样不成立 ⇒ 夹具不合法，不是引擎缺口",
      
      ref: "crates/occt-topo/tests/bop_builder2_boss.rs · BOPAlgo_BOP.cxx:583-711",
    },
        { id: "D5", at: "2026-09-28", chose: "T-28 步 5（合并两个 surface×surface 分派器）不做",
       rejected: ["用统一入口 PatchIntersection 替换 intpatch.rs:138 的通用回退"],
       why: "实测即语义回归：general_bspline_sphere 的曲线离球面 ≈0.67，原 tracer 正确",
      
      ref: "specs/_design_architecture_t25_t28.md · §9.250",
    },
        { id: "D6", at: "2026-09-28", chose: "T-29 终止于「历史横幅」：过期规格保留原样并标注，不重写内容",
       rejected: ["重写 specs/_coverage.md 与 specs/_brepmesh_align_review.md"],
       why: "两文件是历史审查记录，重写会覆盖当时的实测口径",
      
      ref: "specs/_coverage.md:1 · specs/_brepmesh_align_review.md:3",
    },
        { id: "D7", at: "2026-09-28", chose: "T-11 豁免，剩余 ~60 条警告不动；「never used」408 条另立卡（见 T-97）",
       rejected: ["cargo fix --lib 自动清理", "按 cargo check JSON span 自动替换"],
       why: "两条自动路径实测失败（cargo fix 删掉只被 #[cfg(test)] 使用的 import 后报 109 个构建错误；span 不总覆盖整条语句）",
      
      ref: "crates/occt-topo/src · §9.249",
    },
        { id: "D8", at: "2026-09-28", chose: "T-51 豁免：不重做 GeomTrimmedCurve 的 [0,1] 参数表示",
       rejected: ["为移植 gcpnts 两条 Compute arm 先改 TrimmedCurve 参数表示"],
       why: "属架构级改动且无消费方触发那两条 arm；改用 T-96 单独立卡评估",
      
      ref: "crates/occt-geom/src/gcpnts.rs:87-158 · GeomAdaptor_Curve.cxx:239-255",
    },
        { id: "D9", at: "2026-09-29", chose: "§9.270 的假设否证：不再试图用 pcurve 归一化 pass 挽回 T0M 未网格面",
       rejected: ["把被删块的 check_pcurves_and_shift 第二次 pass 加回 reader"],
      
      why: "实测加回后未网格面仍为 7、WIREHIST 逐字不变（1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1），只有已成网格面的密度变化（mesh_t 47616→48140）；真实成因是 6 个面 Delaunay 报 invalid discrete range 后被静默丢弃",
      
      ref: "crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs:441-462 · node_insertion.rs:536",
    },
        { id: "D10", at: "2026-09-29", chose: "按 GT 实测修正「已知收敛差异」的口径：端口丢的 6 个面在 OCCT 都出网格，属真分歧",
       rejected: ["把 T0M 未网格面记成 deflection/启发式导致的容忍差异"],
      
      why: "occt_probe --mesh 0.1 实测：GT 1778 面里 triangles=0 的只有 FACE 1327 一个，而端口的 F405/406/407/408/1691/1758 在 GT 对应序号上分别有 3/502/185/185/84/70 个三角 —— 当时据此判为丢面。**该结论已由 D12/D13 撤回**：两侧面数 1778 vs 1772、分割约定不同，序号不构成身份，且端口这 6 个面的参数域本就零宽",
      
      ref: "specs/occt_probe/occt_probe.exe --mesh 0.1 data/occ/T0M.stp · 撤回见 D12/D13",
    },
        { id: "D11", at: "2026-09-29", chose: "T0M 6 个丢面按两个独立缺陷处理，不做「放宽容差/阈值」式合并修法",
       rejected: ["把 IsValid 判据放宽（改 PConfusion 或加最小宽度阈值）让这些面通过"],
      
      why: "端口 PCONFUSION=1e-9 与 OCCT Precision::PConfusion() 完全一致，拒绝分支也与 BRepMesh_NodeInsertionMeshAlgo.hxx:80-84 逐行一致 ⇒ 判据本身是对的，错的是喂进去的点集；放宽阈值属 OCCT 里不存在的自造规则",
      
      ref: "crates/occt-core/src/precision.rs:12 · BRepMesh_NodeInsertionMeshAlgo.hxx:80-84 · Precision.hxx:334",
    },
        { id: "D12", at: "2026-09-29",
      chose: "撤回 D10：T0M 的 6 个未网格面**不能**判定为「端口丢面」；按「端口与 OCCT 同样拒绝退化参数域」理解",
       rejected: ["按 --mesh 的面序号逐面断言端口丢面", "把「7 vs 1」当成已证的缺陷差"],
      
      why: "① 单面复现证明这些 pcurve 几何正确且参数对齐（F1691 误差 7.1e-15 且 pcRange==edgeRange；F405 误差 8.1e-8/-4）；② 它们在一个参数方向上是常数（F1691 等纬圆 → 常 V；F405 边界整条落在 x=11.010 平面 → 常 U），点集本就零宽；③ OCCT 的 Reset（±1e100 哨兵）+ AddPoint(min/max) + updateRange 收窄 + IsValid(len>PConfusion) 与端口逐行等价 ⇒ 同输入同拒绝。参数域为零宽时 RejectFace 是忠实行为，不是缺陷",
      
      ref: "BRepMesh_DefaultRangeSplitter.cxx:23-88 :202-236 · Precision.hxx:334 · crates/occt-topo/src/meshing/range_splitter/param_set.rs:787-838",
    },
        { id: "D13", at: "2026-09-29", chose: "「1778 vs 1772」的面数差与「7 vs 1」的未网格差解耦：不再用逐面序号对拍这两组量",
       rejected: ["用 occt_probe --mesh 的面序号直接对应端口 F 序号来判丢面"],
      
      why: "GT 面 1327（唯一 triangles=0 者）bbox=(27.75947,-172.18146,-34.05115)-(28.69361,-171.97396,-33.73390)、3 条边，与端口 F1691 bbox=(-2,-27,2.744)-(2,-23,2.744)、1 条闭合边完全不同；且 GT 面数比端口多 6 ⇒ 两侧分割约定不同，序号不构成身份。要判丢面必须按几何键（bbox/面积/曲面类型）配对",
      
      ref: "specs/occt_probe/occt_probe.exe --fbbox 1327",
    },
        { id: "D14", at: "2026-09-29", chose: "撤回 T-98 的「F1758 pcurve d0=NaN」结论：该 pcurve 功能正确，只存在自报参数域不一致",
       rejected: ["按「锥面 NaN pcurve」去改 make_pcurve_full / 补无限域投影分支"],
      
      why: "复核实测：pcRange 虽报 (-inf,inf)，但 d0 在 t∈(0,2π) 九点全部有限，且 |3D−S(uv)| = 3.0e-12 / 2.7e-12；先前的「NaN」是旧探针用 -inf 端点求值造成的假象，不是库缺陷。本卡因此由 P1 降为 P3",
      
      ref: "crates/occt-topo/src/pcurve_full/surface_projector.rs:2876-2891",
    },
        { id: "D15", at: "2026-09-29",
      chose: "T-69 的重定结论：端口拒的 7 个面在 GT 侧都出网格，且已按几何键确认为同一批面 —— 这是真差异，不是序号假象",
       rejected: ["继续用面序号对拍（已由 D13 否证）", "把 7 个未网格面当成忠实的 RejectFace 结案"],
      
      why: "按 bbox 中心最近邻配对（临时 GT --pair 模式 + 端口 zz_pair_probe，均已回退）：P405→G407 tri=185、P406→G408 185、P407→G409 185、P408→G410 185、P1691→G1697 tri=73、P1758→G1764 tri=79、P116→G116 tri=226（dist 0.0–0.30，旧序号法偏差 4–6）。GT 全模型仅 1 个面 triangles=0 ⇒ 端口拒的这批面在 OCCT 都能网格化。结合「端口 pcurve 几何正确（误差 3e-12）」与「端口 IsValid/updateRange 与 OCCT 逐行等价」，分歧只能落在**喂进 range 的离散点集**（pcurve 参数域/取样）上，而非判据",
      
      ref: "crates/occt-topo/src/pcurve_full/make_pcurve.rs:128-139 · BRepMesh_DefaultRangeSplitter.cxx:35-88",
    },
        { id: "D16", at: "2026-09-29",
      chose: "新增测量口径警告：端口的「无 stat 面数」在不同调用序列下不稳定（同一探针内 1757/15 与默认模式 1765/7 不一致）",
       rejected: ["把某一次探针输出的 no_stat 数当成模型性质"],
      
      why: "同一进程内实测：先 prs3d_get_deflection 再 from_deflection 得 stats=1757、no_stat=15、zero_tri=0；而 zz_probe_a3n00 默认模式同参数得 stats=1765、unmatched=7。两者差 8 个面，成因未定（疑与 GeometryRegistry 预热/多次构建有关）⇒ 在弄清前，任何「未网格面数」都必须同时报出测量调用序列",
      
      ref: "crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs:239-274",
    },
        { id: "D17", at: "2026-09-29", chose: "T-94 第一版实现（推迟 seam-shift 的 pcurve 写入到 result 确定后）实测惰性，已回退不提交",
       rejected: ["把「推迟 shift 写入」当成已消除副作用而提交", "以「行为中性」为由保留无实测收益的改动"],
      
      why: "实测与未修改基线逐字相同：T0M wires=2121 / faces_with_2plus=309 / WIREHIST 1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1、stats=1765 unmatched=7；a3n00 stats=226 mesh_v=11101 mesh_t=11121。⇒ shape_fix_face.rs:480-487 不是泄漏源，泄漏更可能在 ComposeShell::perform() 内部（shape_fix_face.rs:604-608）",
      
      ref: "crates/occt-topo/src/shhealing/shape_fix_face.rs:480-487 :604-608",
    },
        { id: "D18", at: "2026-09-29",
      chose: "撤回「F171 中点顶点是 reader 自造/未移植行为」的推论；**并更正为**：忠实的只是 combine_vertex / fix_dummy_seam 的实现，**调用它们的索引是错的**",
      
      rejected: ["按「找 reader 里按中点断边建点的路径并删除/标 UNPORTED」去修 combine_vertex / fix_dummy_seam", "沿用「中点顶点是忠实行为、不存在缺陷」这个半对的结论"],
      
      why: "第一轮 BACKTRACE 定位到唯一调用链 resolve_face(read_topology.rs:714) → fix_notched_edges(wire_fix.rs:3078) → fix_dummy_seam(wire_fix.rs:2995) → combine_vertex(wire_fix.rs:1271)，据此判定 combine_vertex ↔ ShapeBuild_Vertex.cxx:38-72、fix_dummy_seam ↔ ShapeFix_Wire.cxx:4219-4288 逐行相同，且 ShapeFix_Face.cxx:370 关 FixNotchedEdges 那行是**注释掉的**（CR0024983）、第一轮只关 FixLacking（cxx:372）⇒ notch 两轮都开，端口行为正确。**第二轮补核调用点后**发现 fix_notched_edges 传的是 n1，而 cxx:4019 传 i（= n2）⇒ 被调函数忠实、**传参错位**。最终结论「实现忠实 + 调用点有 bug」两轮合并于此（原 D18b 已并入本条并删除）",
      
      ref: "ShapeBuild_Vertex.cxx:38-72 · ShapeFix_Wire.cxx:4019 · :4213-4289 · ShapeFix_Face.cxx:365-520 · crates/occt-topo/src/shhealing/wire_fix.rs:3084",
    },
        { id: "D19", at: "2026-09-29",
      chose: "回退 fix_dummy_seam 的传参修复（`(wire, n1)` → `(wire, n2)`）：代价面量全后判为不划算，不为自己的改动重定断言",
       rejected: [ "保留修复并把 a3n00 的 area_tol 0.15→0.18 以让它转绿", "以「未网格面 7→3」为由直接改基线而不量其余模型", ],
      
      why: "23 个门禁模型的面积比逐一对拍：只有 a3n00 变差（0.8627→0.8298，−0.0329），T0M 仅 +0.0001，其余 21 个完全相同 ⇒ 代价集中在一个模型上，而收益（未网格面 7→3）不在 T-93 的 accept 上。放宽 area_tol 属改既有断言，触发红线第 3 条",
      
      ref: "crates/occt-topo/tests/common/mod.rs:93 · crates/occt-topo/src/shhealing/wire_fix.rs:3084",
    },
        { id: "D20", at: "2026-09-29",
      chose: "看板自审定案：T-93 的 `dependsOn: [\「T-59\」]` 是悬空依赖（T-59 从未立卡），已补立 T-59；并给看板加「悬空依赖」独立风险 + startTurn 的「前置」字段",
      
      rejected: [ "删掉 T-93 的 dependsOn 而不立 T-59（会让 a3n00 面积缺口再次消失成散文）", "把悬空依赖并进已有的「依赖未完成即开工」那条风险里（量级不同：未完成会自己解除，不存在不会）", ],
      
      why: "自审发现 3 条硬问题并已修：① 悬空依赖使 P0 静默挂死（`dependsOn` 全量校验：本次补立 T-59 **之前**，20 张卡里唯一一处 DANGLING，补立后 3 处声明全部 OK）；② 渲染层的 startTurn prompt 不含 dependsOn ⇒ 点「交给 agent」时接手方看不到前置（0.2.0 模板同样不含）；③ nextAction 指 T-69 而 P0 的 T-93 悬空 ⇒ 「现在该做」与 P0 状态互相矛盾。另修：revision r14→r15（已改 3 卡+新增决策却未动版本）、D18 与 D18b 自相矛盾并存（已并成一条）、活动日志缺最近三轮（已补 3 条）",
      
      ref: "specs/board.canvas.tsx · crates/occt-topo/tests/common/mod.rs:91",
    },
        { id: "D21", at: "2026-09-30",
      chose: "撤回 D15 的配对表：端口拒的 7 个面在 GT 的**同名孪生面**上出网格，但 D15 的「bbox 中心最近邻」是错键，配出的 GT 面根本不是同一个面；身份只能按 bbox 六个坐标（≤1e-6）判",
       rejected: [ "沿用 bbox 中心最近邻（不校验尺寸，半径差 0.5 的面可互配）", "按 --mesh/--face 序号配对（D13 已否证）", ],
      
      why: "D15 称 P405→G407 tri=185、P1691→G1697 tri=73、P1758→G1764 tri=79，逐条核错：GT 面 405 的 bbox=(17.6607,-183.7452,-15.4901)-(21.2822,-178.9862,-11.6086)，与 P405(11.010,-223.874,99.877)-(11.010,-218.460,105.289) 相差约 180；GT 面 1697 的 bbox=(-44.5,-43.722,-60.216)-(-39.5,-32.908,-28.809) 与 P1691 相差约 60。按 bbox 六坐标配出的真孪生是 P405→G73、P406→G74、P407→G75、P408→G76、P1691→G409、P1758→G566、P116→G607（P116 的 bboxdiff=0.0000，是唯一被 D15 配对成功的一个）。两者结论同为「GT 侧都出网格」，但身份对了才谈得上逐边对读",
      
      ref: "specs/occt_probe/occt_probe.exe --uvsum（本轮新增）· D15 · D13",
    },
        { id: "D22", at: "2026-09-30",
      chose: "T-69 重定：端口 7 面被拒**不是** RangeSplitter 判据缺陷 —— OCCT 用同一套控制流、拿到同一批点也会拒；分歧在**喂进去的点集**，而点集零宽的原因是**这些面的边界环本身缺边**（读入层）",
      
      rejected: [ "继续在 RangeSplitter / IsValid / updateRange 一侧找缺口（三者已逐行等价）", "按「pcurve 参数域/取样」解释零宽（常宽方向是真实的几何退化边界，不是参数域 bug）", ],
      
      why: "① 判据同构不是推论而是实跑：`collectWirePoints`（BRepMesh_NodeInsertionMeshAlgo.hxx:142-184）+ `AddPoint`（BRepMesh_DefaultRangeSplitter.cxx:35-41）+ `AdjustRange/updateRange`（cxx:45-81 / :202-236）用同一调用链与同一参数域条件枚举过口径 ②/③，对 P1691 的 V 预调整区间 [1.08083900054116766,1.08083900054116877] 判 NO-OP ⇒ 停在一处 ⇒ computeLengthV==0 ⇒ IsValid=false；对 P405 的 U 同理。② 真孪生逐面显示缺边：P1691(1 wire/1 edge/常 V=1.080839001/bbox 在 Z 上零厚，maxZ=2.744) ↔ G409(1 wire/4 edges：e0 底圆 37 点、e1/e3 两条 seam 参数 [1.080839,1.570796]、e2 顶圆 2 点，bbox Z 2.744–3.244)，G409 的 V 展布 0.489957 全部来自端口环里**没有**的 e1/e3；P1758 ↔ G566 dV 1.434376 来自端口单闭合圆里没有的 [−0.717,0.717] 段；P405/406/407/408 ↔ G73/74/75/76 都是 4 边环，端口环只有 1 条边（且其 U 展布 0.429994 就是端口缺的那 0.429995）。③ 面数本身就是证据：GT 1778 面 vs 端口 1772 面（同一 T0M.stp、同一调用序列）",
      
      ref: "specs/_a3n00_gap_analysis.md §9.288 · BRepMesh_NodeInsertionMeshAlgo.hxx:142-184 · BRepMesh_DefaultRangeSplitter.cxx:35-41 :202-236",
    },
        { id: "D23", at: "2026-09-30",
      chose: "T-93 的方向重定向：从「换 `fix_dummy_seam` 传参换 T0M 未网格 7→3」改为「补回缺边」——先查 F1691 的 wire 环为何只剩 1 条边，确认是哪一级丢的（读入层 F1691 就只有 1 条边，说明丢失发生在 reader，不是 mesh 阶段）",
      
      rejected: [ "继续按 §9.287.5 的顺序（先 T-59 面积缺口、再应用 §9.287.4 一行改动）作为主线", "把缺边归到 `make_pcurve` / `pcurve_full`（这些 pcurve 是端口自己造的，但造的是**真实存在的那条**闭合圆）", ],
      
      why: "未网格面的成因既然在边界环，§9.287.4 那行改动（7→3）就不再是这条缺口的主线：它改的是 notch 边对，能少 4 个面但不会补回缺边。且已实测该修复要付 a3n00 −0.0329 的代价（D19），把主线押在它上面不成立。缺边的取证边界很干净：`zz_uv_feed --ids`（本轮新增探针）显示 F1691 在**建完 discret 模型后**仍是 1 wire/1 edge、bbox 在 Z 上零厚，故丢失点不晚于离散模型构建",
      
      ref: "specs/_a3n00_gap_analysis.md §9.288 · crates/occt-topo/src/meshing/model_builder/wire_builder.rs · crates/occt-topo/src/step/read_topology.rs",
    },
        { id: "D24", at: "2026-09-30",
      chose: "把本轮探针固化成可复跑仪器（OCCT `--uv` / `--uvsum`、端口 `zz_uv_feed` + `.target-gate/*.py` 配对脚本），并摘除全部库内调试插桩（`node_insertion.rs` / `param_set.rs` 的 env-gated dump）",
       rejected: [ "把库内 env-gated dump 留在源码里当长期仪器", "不固化仪器，下一轮再临时写一遍", ],
      
      why: "库内插桩污染被测代码（且读数与被测二进制同源，难判「插桩本身是否改变行为」）；摘除后复跑 `zz_uv_feed --ids` 与摘除前逐字相同（1772 面、同一 7 个被拒面的 index+bbox，key= 指针串除外）⇒ 插桩确实是惰性的，结论可归因到未插桩的二进制。反向的问题是「下一轮还得重写」，故把端口侧探针留在 `crates/occt-topo/examples/zz_uv_feed.rs`（不在库路径上），并在 §9.288 写明命令",
      
      ref: "crates/occt-topo/examples/zz_uv_feed.rs · specs/occt_probe/occt_probe.cpp · §9.288.5",
    },
        { id: "D25", at: "2026-09-30",
      chose: "T-93 主线定位到**导入期 seam 修复**：STEP 里这些面声明的环边数本来就少于 OCCT 成形后的边数，端口差的是 `ShapeFix_Shape::Perform → ShapeFix_Face::Perform`（`FixMissingSeam`）这一段；但端口现有的 `ShapeFixFace::fix_missing_seam` **不能整形状照用**（212/1772 面被误改，只有 1 个是真需要的），必须先按 `.cxx` 补完再接线",
      
      rejected: [ "整形状直接调用现有 fix_missing_seam（实测 212 面变动，含 126 个 2→4 的伪增量）", "按 §9.287.4 的 fix_dummy_seam 传参修复当主线（它换 notch 边对，不补缺边）", ],
      
      why: "取证链（全部可复跑）：① STEP 文本实测 `#29920=ADVANCED_FACE('',(#4223),#77,.T.)`（球 R=4.25，心 (0,-25,-1.006)）→ `#4223=FACE_OUTER_BOUND('',#6320,.T.)` → `#6320=EDGE_LOOP('',(#25460))` → `#25460=ORIENTED_EDGE('',*,*,#17469,.F.)`，即**该球带面在文件里只声明了 1 条闭合边**（圆 R=2，心 (0,-25,2.744)）—— 端口也正好读到 1 条，**reader 未丢边**（定级见 D26）。② OCCT 成形后同面（探针 `--typefaces 3 4.25` 按 R 定位， healed 顺序 f=1697 / no-heal 顺序 f=1691）是 **4 条边**：底圆 37 点（pcRange 与 3D 参数域同为 [π, 3π]）、两条 seam（参数域 [1.08083900054116833, 1.57079632679489656]）、顶点退化边，且 `vrange` 都是 [-π/2, π/2] — 差别不在曲面参数域，而在**度量的 3D 范围**：healed 的 bbox z∈[2.744,3.244]，no-heal 的 z∈[-5.256,3.244]，两侧环边数都是 1。③ 端口同一面（`zz_uv_feed --model 1691`）只有 **1 条边**（curve=Circle、par=[π,3π]、dpar=2π、bbox z∈[2.744,2.744]），与 no-heal 的环边数一致。④ 端口 `ShapeFixFace::fix_missing_seam` 单独作用在该面上**能从 1 条边生成 4 条边**（即缺的这段控制流端口只做了一半）。⑤ 但整形状跑一遍：`faces=1772 fix_returns_true=212 edge_count_changed=212`，其中 **1→4 只有 f=1691 一个**，其余是 126 个 2→4、46 个 2→5、6 个 3→6、5→11/8→11 等伪增量；这源自 `shape_fix_face.rs:147-148` 声明的三处 UNPORTED（seam 构造 `cxx:1899-2330`、`w2 != null` 合并、post-seam 面循环）。⇒ 方向正确、接线未就绪",
      
      ref: "data/occ/T0M.stp #29920/#4223/#6320/#25460 · ShapeFix_Face.cxx:492-498 :1722-2330 · ShapeProcess_OperLibrary.cxx:801-899 · crates/occt-topo/src/shhealing/shape_fix_face.rs:147-148 · .target-gate/fixall.txt",
    },
        { id: "D26", at: "2026-09-30",
      chose: "定级结论：**reader 没有丢边**。端口 F1691 的 1 条边就是 STEP `EDGE_LOOP` 声明的 1 条边 —— 缺口不在 `read_topology.rs` 的 `resolve_face`/`resolve_loop`，而在 reader **之后**缺了 OCCT 的导入期 seam 构造（D25）。故修复点是「新增一段同等分支」，不是「修 reader 的丢边 bug」",
      
      rejected: [ "继续在 resolve_face / resolve_loop / bind_edge_loop_vertices 里找「哪条边被丢了」——没有可丢的", "把 1 条边当成 reader 的回归（它忠实于文件）", ],
      
      why: "逐级核过，每一级都是 `file → 逐条进 → 无过滤`：① STEP `#29920=ADVANCED_FACE('',(#4223),#77,.T.)` → `#4223=FACE_OUTER_BOUND('',#6320,.T.)` → `#6320=EDGE_LOOP('',(#25460))`，环里 1 个 ORIENTED_EDGE。② 端口 `resolve_face`（`read_topology.rs:578`，bounds 循环 `:649-662`）→ `resolve_outer_bound`（`:564`，只加朝向）→ `resolve_loop`（`:413`）→ `resolve_oriented_edge`（`:369`）→ `resolve_edge`（`:307`）。③ `resolve_loop`（`:421-432`）对 `items` 逐条 `resolve_shape`，两侧只有 `if !s.is_edge() { return Err(..) }`，**没有 filter/retain/skip/去重**，1 条进 1 条出。④ `:628-644` 的 `VERTEX_LOOP` 自然边界分支**不适用**（该 bound 是 `EDGE_LOOP` 不是 `VERTEX_LOOP`，`bound_loop_is_vertex_loop` 为 false），所以不走那条。⑤ 实测等式：STEP 声明 1 条 ↔ 端口 `zz_uv_feed --model 1691` 报 `wire[0] nEdges=1`（`edge=3980 curve=Circle par=[3.141592653589793,9.424777960769379] dpar=6.283185307179586`，即完整 2π 闭合边）。⑥ 反向也成立：GT 未修复导入（`--nofix`）的这个面同样是 1 条边、同样 `par=[π,3π]` —— 两侧读入一致，差别只在 OCCT 随后跑了 seam 构造。⇒ 结论是「缺一段控制流」而非「丢了一条边」，D23 里「缺边不晚于离散模型构建」的措辞据此收紧为「**文件声明本身就只有 1 条**」",
      
      ref: "crates/occt-topo/src/step/read_topology.rs:307 :369 :413 :564 :578 :628-644 · data/occ/T0M.stp #6320 · D25 · §9.289.7",
    },
        { id: "D27", at: "2026-09-30",
      chose: "T-93 落地：把导入期 seam 修复**接进 reader**（`resolve_face` 末尾恢复 `fix_missing_seam` + `check_pcurves_and_shift`）。实测 T0M 未网格 7→6、wires 2121→1931、WIREHIST 向 GT 收敛，且 `--lib` 1281/0、`step_obj_gates` 5/5 全绿 ⇒ **提交**",
      
      rejected: [ "继续把 `fix_missing_seam` 留在 reader 之外（D2/T-93(a) 的处理）—— 实测 OCCT 成形期确实加了 seam，不接就等于永久少 3 条边", "为「消灭」211 个伪增量去改 `check_wire` / 配对选择（没有 .cxx 依据说明端口逻辑与 OCCT 不同）", ],
      
      why: "① 先证同构：`zz_seam_fix` 对 T0M f=1691 单独调用，产出 4 条边 —— 底圆 `par=[π,3π]`、两条 seam `par=[1.080839000541168,1.570796326794897]`、`deg=true` 顶点边，bbox 由 z∈[2.744,2.744] 变 [2.744,3.244]，与 OCCT 成形后（`--typefaces 3 4.25`）**逐项相同**。② 再证缺口是调用点：临时计数器（已摘除）量到 `read_step_file` 全程 `SEAMCALLS = 0`；而 OCCT 侧 `STEPControl_Controller.cxx:201` 设 `FromSTEP.exec.op=FixShape`、`:221` 设 `FixMissingSeamMode=-1`（NeedFix 为真）⇒ 默认就开。**这一条与 D2 的依据相反**：D2 靠「patch 过的 OCCT 在 TransferRoots 无输出」，但本轮实测 OCCT `OneShape()` 的面确实比文件多 3 条边、且关掉 ShapeProcess 就没了；可复跑的那一侧指向「成形期确实加了 seam」，故 D2 的 link-map 结论不再作为依据。③ 对照实测（只切这一个 hunk）：T0M 未网格 7→**6**、wires 2121→**1931**、多 wire 面 309→**120**、`WIREHIST 1:1463 2:294 3:9 …`→**`1:1652 2:106 3:8 …`**；a3n00 stats 226→225、mesh_t 11121→11941；`--lib` **1281/0**、`step_obj_gates` **5/5**。④ 逐面核对：f=1691 已出网格 mt=143 且 bbox 与 GT 1697 逐字相同；f=405/406/407/408 全部治好且 bbox 与 GT 73/74/75/76 逐字相同。⑤ 源码注释里「UNPORTED: seam 构造 cxx:1899-2330」是**过期的** —— 那段早已移植，已同步更正",
      
      ref: "specs/_a3n00_gap_analysis.md §9.290 · ShapeFix_Face.cxx:492-498 :1899-2264 · ShapeProcess_OperLibrary.cxx:785-899 · crates/occt-topo/src/step/read_topology.rs",
    },
    {
      id: "D28",
      at: "2026-09-30",
      chose: "T-99 主线定在**读入期的 bound 合并**（把 2 个 FACE_BOUND 并成 1 条由 seam 边闭合的 wire），不再怀疑 Delaunay 层",
      rejected: [ "继续在 Delaunay 侧（frontier_adjust / cleanup_mesh / 点序）找缺口",
        "继续按 §9.365.2 的按序号 GT 表当目标值" ],
      why: "① 等价性已证：把端口的完整结构（点 + 72 条约束链 + 容差/格/顶点序）喂给 OCCT 的 BRepMesh_Delaun，18 个面逐字段相同 ⇒ 这层没差异；② §9.367 的 domain=0 是结构里没有 Frontier 链时 cleanupMesh 删光所致（cxx:1028 / :832-835 / :908-911），不是点集证据；③ 真正的输入差是面拓扑：STEP 声明 2 个 FACE_BOUND，OCCT 读入期 ShapeProcess 合并成 1 条 wire（--faceids ± --nofix 实测），端口仍是 2 条；④ a3n00 两侧面序号只有 3/226 对应，按序号的 GT 表已作废",
      ref: "specs/_a3n00_gap_analysis.md §9.368 · ShapeFix_Face.cxx:492-498 :1722-2330 · crates/occt-topo/src/shhealing/shape_fix_face.rs:147 :406",
    },
    {
      id: "D29",
      at: "2026-09-30",
      chose: "把导入期 `ShapeFix_Face::FixMissingSeam` 步骤**接回 reader** —— 否定 T-93(a)/§9.219 的移除理由；T-99 按 accept 结项",
      rejected: [ "继续把 seam 合并留在 reader 之外（靠 rescue 补三角形）",
        "用 repr_key 摘注册一类旁路去换同一个面积比（§9.315/§9.320）" ],
      why: "① 移除理由是「OCCT 在这条路径上不跑 `ShapeFix_Face::Perform`」（link map 实验，§9.219）；本轮用 `--faceids` 开/关 ShapeProcess 实测得 `{1:208,2:9,4:1,6:4,10:4}` vs `{1:163,2:53,4:2,6:4,10:4}`，后者与 STEP 的 FACE_BOUND 直方图、与端口逐字相同 ⇒ 导入确实合并了 45 个 2-bound 面。② 接回后法兰面结构与 OCCT 模型孪生逐项相同（1 wire / 5 边 / 7 pcurve）。③ 整模型：a3n00 面积比 0.8627→0.8996、T0M 未网格 7→6、T0M 0.9786→0.9987、acs10 0.9025→0.9846、rescue 16→0，且未改任何断言或 area_tol。④ 剩余差异已量化为 T-101（还差 1 个 2-bound + 1 个 4-bound 面；全形状多改 4 个面；法兰面密度 72 vs GT 52）",
      ref: "specs/_a3n00_gap_analysis.md §9.369 · crates/occt-topo/src/step/read_topology.rs · ShapeFix_Face.cxx:482-498 :1722-2330 · D27 · D28",
    },
    {
      id: "D30",
      at: "2026-10-03",
      chose: "T0M 的薄面簇对拍用**同一次运行内的 OCCT 真值**做（新增测点 = 既有 `--facestats <defl> <angle>` / `--wires` ↔ `--fstats` / `--ecensus` / `--wirehist`），**保留 ④**（`MapReShape::apply_impl`）；T-101 收口，残差另立 T-102",
      rejected: [ "把 ④ 当成 T0M 0.9931→0.9809 的回归撤掉",
        "继续按 3D bbox 逐面配对来判 ④ 的优劣（OCCT 与端口的 bbox 在第 6 位小数有差，1778 个里 591 个对不上）" ],
      why: "① 测点同源：偏转取门禁同口径 `prs3d_get_deflection` = maxComp(bbox)*0.001*4 = 2.574312、角度 20°，OCCT 侧 `--facestats 2.574312 0.349066` 得 1778 面逐面 nodes/triangles（0 张未出网格），`--wires` 得 `TOTAL faces=1778 wires=1921 ... WIREHIST 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`。② 端口 ④ 开 = `1772 / 1915 / 1:1666 2:92 3:8 5:1 6:1 7:2 8:2` ⇒ **多 wire 桶与 OCCT 逐项相同**，残差只有『少 6 张单 wire 面』；④ 关 = `1772 / 1921 / 1:1662 2:96 3:8 5:1 6:1 7:2 8:1 10:1` ⇒ 多 wire 桶错 4 项，只是 wire 总数凑回 1921。③ ④ 不增删面：两态都 1772 面 / 1772 个网格槽位，未出网格的恒为 idx1758（1 边 1 顶点退化面）+ ④ 开的 idx116 或 ④ 关的 idx351；idx116/idx351 的 bbox 与 OCCT face 607/628 **逐位相同**，OCCT 两张都出网格（233/231、207/219）。④ 逐面网格：38 个两态不同的 bbox 里 OCCT 有精确同 bbox 的 21 个，④ 开 **21/21 逐位相同**（含 20 张 `ON[4/2] OFF[9/7] OCCT[4/2]`），④ 关 0/21。⑤ 面积账：丢 idx351（X 区）得 191332.77 / 0.9931，丢 idx116（Y 区）得 188978.93 / 0.9809，差 2353.8 = idx116 在 ④ 关时的面积 ⇒ 0.9931→0.9809 分毫不差就是『换了一张面不出网格』。⑥ 真正的残差与 ④ 无关：端口少 6 张单 wire 面，且 idx116/idx351 各应是 1 条 wire（OCCT 两个面都是 1 条），端口留了 2 条 / 8 条 ⇒ 属 `CollectWires`/`MakeFacesOnPatch` 同一族（§9.579-§9.580）。⑦ 本轮零库改动；④ 关的 `--ecensus` 用 `MapReShape::apply` 上的临时 `MAPRESHAPE_DIRECT` 开关跑出，已还原。",
      ref: "specs/_a3n00_gap_analysis.md §9.581 · specs/occt_probe/occt_probe.cpp (--facestats/--wires) · crates/occt-topo/examples/zz_probe_a3n00.rs (--fstats/--ecensus/--wirehist) · ShapeBuild_ReShape.cxx:200-324 · D29",
    },
    {
      id: "D31",
      at: "2026-10-03",
      chose: "T0M 的 6 张残差面按 `ShapeFix_Face::FixMissingSeam` 返回多面结果这条 .cxx 分支修 reader：`resolve_face` 接受 Face 或 Shell 结果，`resolve_shell` 把 Shell 结果的面逐张加进本 shell（= `Context()->Replace(myFace, myResult)` + `Context()->Apply(S)`）",
      rejected: [ "按 §9.580-F 去改 `CollectWires` 的 wire 合并（该假设已证伪：OCCT face 116 是 healing 新建面，不是端口同 bbox 的那张原始面）",
        "为了凑 1778 去加任何 OCCT 里没有的谓词/启发式（例如『若面是周期面就拆两张』）" ],
      why: "① 先量清方向：T0M.stp 自身 1772 个 `ADVANCED_FACE`，OCCT `--nofix` 也是 1772（`1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1`），开 ShapeProcess 才 1778（`1:1672 2:92 …`）⇒ 是 OCCT **多造 6 张**，端口的面数等于文件本体、wire 直方图等于治好后的形状（把 seam 合并做了、面没落地）。② `STEPControl_Controller.cxx:201` 定死 `FromSTEP.exec.op = \"FixShape\"` ⇒ 这 6 张全部出自 `ShapeFix_Shape::Perform`，与 `SplitClosedFaces` 等算子无关。③ `.cxx` 里唯一能把面数变多的分支就是 `ShapeFix_Face.cxx:2266 myResult = CompShell.Result()`（可为 Shell）/ `:2268 Context()->Replace(myFace, myResult)`（无类型判断），由 `ShapeFix_Shape.cxx:257/294 Context()->Apply(S)` 物化，`ShapeBuild_ReShape.cxx:282-299` 规定『替换体与原类型不同时只取其中与原类型相同的子形状』。端口 `read_topology.rs::resolve_face` 只接 `ShapeType::Face`，Shell 直接 fallthrough 成原面 ⇒ 落地的正是这一处缺口。④ 实测收益：T0M 面积比 0.9809→1.0084、acs10 0.9928→1.0045，23 个模型里其余 21 个逐位不变（含 a3n00 1.0000、TDB 1.0092），5/5 门禁通过。⑤ 仍差 3 张已分别定位到两条下一入口（wire 首段 `cxx:365-480`；face 1752 的合缝分歧），不合并进本轮。",
      ref: "specs/_a3n00_gap_analysis.md §9.582 · crates/occt-topo/src/step/read_topology.rs · ShapeFix_Face.cxx:2266-2268 · ShapeFix_Shape.cxx:200 :257 :294 · ShapeBuild_ReShape.cxx:282-299 · STEPControl_Controller.cxx:201 · D30",
    },
    {
      id: "D32",
      at: "2026-10-03",
      chose: "T-102 剩余缺口按 `FixShape` 子模式 bisect 钉成两条独立入口：**(1) 2 张 loop-wire 面**（端口 census 628/635 现为 1 wire，OCCT 为 2 wire）—— 按 `ShapeFix_Face.cxx:583` 移植 `FixLoopWire`（`:2478`）+ `FindNext`（`:2398`）+ `ShapeAnalysis_Wire::CheckLoop`（`ShapeAnalysis_Wire.cxx:2228`）；**(2) 3 张碎片面** —— 仍走 `FixMissingSeam` 那条线（§9.582-E）。bisect 仪器改用 `occt_probe --fix <DE_ShapeFixParameters 字段> <int>`",
      rejected: [ "继续用 `--set FromSTEP.FixShape.*`（或读前的 `SetShapeFixParameters` / `SetShapeProcessFlags`）做 bisect —— actor 是 `ReadFile` 期间才建起来的（`XSControl_Reader.cxx:480-534`），读前调用是 no-op；8.0 又走 `DE_ShapeFixParameters` 显式 map（`XSAlgo_ShapeProcessor.cxx:609-620`），静态量根本不被读。本轮 T0M 上 14 个单点 `--set` 结果全部等于默认值即为证",
        "把 loop-wire 面按『周期面就拆两条 wire』这类 OCCT 里没有的启发式补上（正确分支是 `FixLoopWire` 的 `CheckLoop` 判据）",
        "把 `FixMissingSeamMode=0` 的 1973 面当成『seam 步造了 201 张面』（实际是 seam 缺席后 `FixAddNaturalBound` 不再兜住、`FixSplitFace` 在 `cxx:711-717` 开火）" ],
      why: "① `FixFaceMode=0` ⇒ T0M 1772 面、`1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1`，与 `--nofix` 逐位相同 ⇒ 比文件本体多出的 6 张面唯一来源是 `ShapeFix_Shape.cxx:196-210` 的 face 分支（D31 机制坐实，与 `SplitClosedFaces` 等算子无关）。② `FixLoopWiresMode=0` ⇒ `faces=1778` 但 `wires 1921→1919`、`2:92→2:90` ⇒ 只有 `ShapeFix_Face.cxx:583-598` 的 `FixLoopWire` 动了 wire（面数不变）；用 `--faceids` 在开关前后逐 bbox 对差，只差 2 个 bbox = 端口 census 628（12 边）/ 635（16 边），OCCT 侧同 bbox 都是 `wires=2` ⇒ 这 2 张面就是端口缺的 `+2 wire / 0 面`。③ 其余单点关（AddNaturalBound / SplitFace / SmallAreaWire / RemoveSmallAreaFace / Orientation / IntersectingWires / Reorder / Small / Connected / AutoCorrectPrecision / Seam / Shifted / NotchedEdges / FixWireMode）与默认逐位相同；`FixMissingSeamMode=0` 反而 1973 面（+`FixSplitFaceMode=0` 回到 1772）⇒ +201 是 `FixAddNaturalBound` 的 `NeedSplit` 交互，不是 seam 步造面。④ 端口当前 `faces=1775 wires=1916 1:1671 2:90 3:8 5:1 6:1 7:2 8:2`，与 OCCT 默认的差 = 上面 2 张 loop-wire 面 + 3 张碎片面，缺口自此是**两个独立可验收的子项**。⑤ 本轮 Rust 侧零改动（只加探针 TEMP 仪器 + 规格文档），门禁值沿用 §9.582-D。",
      ref: "specs/_a3n00_gap_analysis.md §9.583 · specs/occt_probe/occt_probe.cpp (--fix/--set) · ShapeAnalysis_Wire.cxx:2228-2340 · ShapeFix_Face.cxx:583-598 :2398 :2478 · ShapeFix_Shape.cxx:196-210 :705-717 · XSControl_Reader.cxx:480-534 · XSAlgo_ShapeProcessor.cxx:609-620 · D31",
    },
    { id: "D33", at: "2026-10-04",
      chose: "T-102 的剩余缺口**不是** `ComposeShell` 少切/少分类（D32 的两条入口因此作废），根因是 `LoadWires` 里 `ShapeExtend_WireData` 的边被**二次反转**：`shape_fix_compose_shell/load_wires.rs::wire_data_edges` 去掉 `wire.0.orientation() == Reversed` 时的逐边 `reverse()`，直接返回 `edges_of_wire(wire)`",
      rejected: [ "按 D32/§9.586-F(1) 去 `ShapeFix_ComposeShell.cxx:206-270` 的 `Perform` 分段里找『端口少切/少分类的那一段』（实测该处不缺：修朝向后面数与 OCCT 逐项相同）",
        "按 §9.585-C(2) 去 `ShapeFix_Face.cxx:2236-2325` 找『端口缺的 false 分支』（`surf=1300` 的多拆同样由这个二次反转引起，修掉后 face=1758 = OCCT 的 1 wire / 9 边）",
        "为凑 1778 面去移植 `FixLoopWire` 或加任何 OCCT 里没有的谓词（`FixLoopWiresMode=0` 那条 bisect 结论仍成立，但它不是本轮残差的入口）" ],
      why: "① 对着 `.cxx` 定因：`ShapeFix_ComposeShell.cxx:629` 的 `TopoDS_Iterator aIt(wire)`（cumOri 默认 true）已把 wire 朝向与子边朝向复合一次，`cxx:635 sbwdM->Add(E)` 存的就是这个复合值；端口 `topo_tools_full::edges_of_wire` 走 `cumulated_children`、语义相同 ⇒ REVERSED wire 再逐边 `reverse()` 就是二次反转，`cxx:586-607` 的 2D `WireOrder` 因此拿到反向 pcurve（周期面判 `Shifted [4,3,2,1]` 而非 `Same [1,2,3,4]`），`FixReorder` 后的 `sbwdM` 朝向与 OCCT 相反。② 修掉后 T0M `--wirehist` = `TOTAL faces=1778 wires=1921` + `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`，与 `--wires` 真值逐项相同（修复前 1775/1918/`1:1669 …`）；§9.585-C 的四张环面碎片（OCCT faceid 1462/1503/1648/1668）端口现在都在（`--ecensus` 0-based 1461/1502/1647/1667，`wires=1 edges=[2]`），`surf=1300` 的多拆消失（face=1758 = `wires=1 edges=[9]`）。③ A/B（23 模型面积表）：只有 T0M 动（1.0055→1.0053），其余 22 个逐位不变（a3n00 1.0000 / acs10 0.9976 / TDB 1.0110 …），两态 `step_obj_gates` 5/5；`--wirehist` 的 a3n00/acs10/TDB 与 §9.586 逐项相同 ⇒ 该修复是纯对齐性修复。④ 仪器（`LOADDBG`/`COMPDBG`/`BWDBG`/`SEAMDBG`）已全部撤除。",
      ref: "specs/_a3n00_gap_analysis.md §9.587 · crates/occt-topo/src/shape_fix_compose_shell/load_wires.rs · ShapeFix_ComposeShell.cxx:629-640 :586-607 · TopoDS_Iterator.cxx:26-83 · D32",
    },
    { id: "D34", at: "2026-10-04",
      chose: "`DispatchWires` 里「给新建边补 3D 曲线」改走真实 `BRepLib::BuildCurve3d`（平面臂 `GeomLib::To3d`，一般臂 `GeomLib::BuildCurve3d` 的 `AdvApprox`），并移植同族的 `BRepLib::CheckSameRange`/`SameRange`、`ShapeFix_Edge::TempSameRange`/`FixAddCurve3d(edge)`；不再用 `meshing::edge_discret::CurveOnSurface` 适配器冒充 3D 曲线",
      rejected: [ "保留 §9.30 的旧体（`CurveOnSurface` 当 3D 曲线）—— 它给不出一条 OCCT 意义上的 3D 曲线，`SameRange`/容差/后续离散都不同源",
        "为凑某个模型的读数去调 `BuildCurve3d` 的连续性/阶数/`MaxSegment`（应走 `BRepLib.hxx:90-94` 的默认 `GeomAbs_C1` / `MaxDegree=14` / `MaxSegment=0` + `evaluateMaxSegment`）" ],
      why: "① 对着 `.cxx`：`ShapeFix_Edge::FixAddCurve3d`（`cxx:618-638`）在 `!SameRange` 时先 `TempSameRange`（`:626-629`），再 `ShapeBuild_Edge::BuildCurve3d`（`:631-635`）；后者即 `BRepLib::BuildCurve3d(edge, max(1e-5, Tol))`（`ShapeBuild_Edge.cxx:714-775` ↔ `BRepLib.cxx:301-455`），平面臂 `GeomLib::To3d`、一般臂 `GeomLib::BuildCurve3d`（`GeomLib.cxx:1051-1165`）。② 旧端口只把 `CurveOnSurface` 存成 3D 曲线，拓扑上「已有 3D 曲线」表面成立，但曲线不是 OCCT 会产出的那条。③ 结果：4 个基准模型的 `--wirehist` / `--ecensus` 与 §9.587 逐项相同 ⇒ 这是纯对齐性补齐，不改已对齐读数。",
      ref: "specs/_a3n00_gap_analysis.md §9.588 · BRepLib.cxx:149-183 :187-263 :275-295 :301-455 · ShapeBuild_Edge.cxx:714-775 · ShapeFix_Edge.cxx:335-464 :618-638 · GeomLib.cxx:559-679 :1051-1165 :2991-3077 :3079-3221 · ElCLib.cxx:1339-1422 · ShapeFix_ComposeShell.cxx:3506-3529 · D33",
    },
    { id: "D35", at: "2026-10-04",
      chose: "`ShapeFix_Face::Perform` 的周期锥面退化臂按 `.cxx` 顺序补齐：读入链先 `FixPeriodicDegenerated`（`ShapeFix_Face.cxx:486-489`，默认 `myFixPeriodicDegeneratedMode=-1` ⇒ `NeedFix==true`）再 `FixMissingSeam`；移植 `PeriodicConicalLoop`/`is_periodic_conical_loop`（`cxx:3018-3098`）与 `fix_periodic_degenerated`（`cxx:3101-3259`），并在 `fix_missing_seam` 开头补 `myFace = Context()->Apply(myFace)`（`cxx:1737-1741`）",
      rejected: [ "按面级 bbox/面积/密度门把「塌陷薄面」筛掉或跳过（OCCT 里不存在的规则）",
        "为 T0M 的读数去调 `FixPeriodicDegeneratedMode`/`RecadreOnPeriodic` 的默认值或偏转（默认值本身就是规格）" ],
      why: "① §9.590 钉死第一现场：face 1764 的 wire 塌成 `edges=[1] nv=1`、V 跨度 0 ⇒ `invalid discrete range`，`--fstats` 只出 1777（少 1 张失败面）。② 对着 `.cxx`：缺的正是 `Perform`（`cxx:482-498`）里 `FixMissingSeam` 之前那段 `FixPeriodicDegenerated`，端口只移植了尾部。③ 落地后 `--fstats` 1777→1778（失败面归零），face 1764 `mv/mt=38/35`、bbox `(-23.729500,-207.157336,-107.158859)-(-21.270500,-206.058541,-104.833829)` 与 OCCT `FSTAT_OCC face=566` 逐位相同；`--ecensus face=1764` = `wires=1 edges=[5] nv=2`（修前 `[1] nv=1`），进到「seam 边 + 退化 apex 边」形态。④ `step_obj_gates` 5/5，T0M 面积比 1.0053（tol 0.025 内）。旁支（未动）：端口 `edges=[5]` vs OCCT `4`，面级密度约半（参数域 `nbU/nbV` 层）。",
      ref: "specs/_a3n00_gap_analysis.md §9.590 · §9.591 · ShapeFix_Face.cxx:144 :482-498 :1737-1741 :3018-3098 :3101-3259 · ShapeFix_Root.lxx:101 · ShapeProcess_OperLibrary.cxx:830 · D32",
    },
    { id: "D36", at: "2026-10-05",
      chose: "删除 `meshing/model_builder/wire_builder.rs::add_wire` 里 OCCT 不存在的「单 wire 面若 `order.chain_area() < 0.0` 就 `order.reverse_chain()`」后处理；`CheckOrder` 之后与 OCCT 一样只读 `Ordered(i)`、经 `ShapeExtend_WireData::Edge(signed)` 逐条反转边",
      rejected: [ "把 `chain_area()` 换成真实 2D `wire_area_2d()`、或加 `|area|` 阈值/退化保护 —— 该分支本身在 OCCT 里就不存在，改成「更好的启发式」仍是非 OCCT 规则",
        "只对 `chain_area == -0.0` 这一处特判（等于给 motoc 调参，换模型即失效）" ],
      why: "① 现象：`data/occ/motoc.step` 大圆柱 / 右侧小圆柱上多出「横向 2 层环带」，OCCT 真值 `data/occ/occ-motoc.obj` 13944 三角无此环带；按 bbox 六坐标配对 `TOTAL faces=223` 两侧一致 ⇒ 多的是网格不是拓扑面（14 张平面端盖 dt=-19 + 1 张 dt=-32 = 298）。② 根因：`chain_area()`（`meshing/wire_order.rs`）只用每条边的**端点**组成多边形算有向面积，含圆弧的单 wire 面上端点共弦 ⇒ 面积退化到 `-0.0`，浮点噪声令 `-0.0 < 0.0` 成立 ⇒ 整链被反向 ⇒ Delaunay 把内孔侧当域内、填 17 个扇形三角封住内孔。③ 对着 `.cxx`：`BRepMesh_ShapeVisitor::addWire`（`BRepMesh_ShapeVisitor.cxx:99-134`）在 `aWireTool.CheckOrder(aOrderTool, true, false)`（`cxx:103`）之后**只**逐条反转边，**没有**任何整链面积反转分支。④ 落地后 motoc 三角数与 `occ-motoc.obj` **逐位相同（13944）**；四模型 `--wirehist` 与 `step_obj_gates` 5/5 不动。",
      ref: "specs/_a3n00_gap_analysis.md §9.598 · §9.601 · BRepMesh_ShapeVisitor.cxx:99-134 :103 · ShapeExtend_WireData.cxx:583-588 · meshing/wire_order.rs · D34",
    },
    { id: "D37", at: "2026-10-06",
      chose: "`geom2d` 的 `Geom2dConvert` 全链（`CurveToBSplineCurve` 的圆锥 / `Geom2d_OffsetCurve` 臂、`CompCurveToBSplineCurve`、`ApproxCurve` + 2D `AdvApprox_ApproxAFunction`）与 `Gp*` 基元变换（`GpAx22d`/`GpCirc2d`/`GpLin2d`/`GpAx2d`/`GpElips2d`/`GpHypr2d`/`GpParab2d`/`GpDir[2d]`/`GpAx[123]`/`GpTrsf[2d]`/`GpMat[2d]`/`GpQuaternion`/`GpEulerSequence`/`GpGTrsf[2d]`/`GpVec[2d]`）逐条对着 `.cxx`/`.hxx` 正本清源，缺件标 UNPORTED（`GetMat4`/`SetMat4` 落地，`DumpJson`/`InitFromJson` 登记）",
      rejected: [ "继续用 `ReparamCurve2d` 包装器冒充 `GeomLib::SameRange` 的产出（§9.605：包装器不是 OCCT 会持有的类型，后接分派全落错臂）",
        "按某个模型的面积 / 三角数去调这些基元的常量或分支（基元是共享的，且这些臂在当前 23 模型上不改变行为）" ],
      why: "① §9.604-§9.623 是一条**潜伏缺口**链：这些分支在 23 个基准模型上不改变任何读数（四模型 `--wirehist` / `step_obj_gates` / `fuse_box_cylinder` / `export_data_obj` 逐项不变），但缺了它们，某些输入会落到错误的臂上，且失配无法归因到 OCCT 控制流。② 关键修复：`geom_lib_same_range` 的 B-spline 臂改为返回 `SameRange` 真正持有的 `Geom2d_BSplineCurve`（§9.605）、圆锥臂走 `Geom2dConvert::CurveToBSplineCurve`（§9.606，`Geom2dConvert.cxx:379-397`）、offset 臂走 `Geom2dConvert_ApproxCurve`（§9.609/§9.610）、等跨 / 非等跨「其余」臂 `TrimmedCurve → CurveToBSplineCurve → Reparametrize`（§9.612）。③ motoc 223 面 `maxcoord 差 > 0.01` 的面 **4 → 0**（§9.604）；余 4 张 `Cylinder` 面的 bbox 差 = 2D pcurve 子区间臂未移植 `Geom2d_BSplineCurve::Segment`。④ 验证：`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15、`occt-core` / `occt-geom` 单测逐项同基线、四模型 `--wirehist` 逐项不变。",
      ref: "specs/_a3n00_gap_analysis.md §9.604-§9.623 · Geom2dConvert.cxx:244-262 :283-301 :323-345 :379-397 · GeomLib.cxx:871-888 :908-921 :960-968 · ShapeBuild_Edge.cxx:643-655 :687-697 · gp_Ax22d.hxx · gp_Circ2d.hxx · gp_Trsf2d.cxx:198-213 · gp_Trsf.cxx:90-99 :103-109 :207-213 · gp_Mat.cxx:242-334 · gp_Quaternion.cxx · gp_EulerSequence.hxx · gp_GTrsf.cxx · gp_Dir.cxx:129-155 · D34",
    },
  ],
  nextAction:   { taskId: "T-103",
    action: "【T-103 下一步（§9.608，本轮闭环）】：ON-only 塌陷面 **8/8 归位**。根因 = 2D pcurve 求值用了边的 3D range 而非 `BRep_Tool::CurveOnSurface` 返回的 COS 表示自身 range（`BRep_Tool.cxx:327-361` 的 `cxx:353` `GC->Range(First,Last)`；入口链 `ShapeAnalysis_Edge::PCurve`（`ShapeAnalysis_Edge.cxx:192-208`）→ `ShapeAnalysis_Wire::CheckOrder`（`ShapeAnalysis_Wire.cxx:593-651`，`cxx:634/648-649`））；拆缝子边的 range 被 `ShapeBuild_Edge::CopyRanges`（`ShapeBuild_Edge.cxx:206-334`）重标到切片窗口，用父边 3D range 求值就偏一个周期 ⇒ 分类器多边形自重叠。修复：`wire_builder.rs` 新增 `pcurve_and_range`（COS range，缺则回退 3D range），5 处 pcurve 消费点统一改用。数字（T0M，单二进制 A/B，默认路径 diff=0）：ON `pf=28 1→74`、`pf=812 20→70`、`pf=1663 38→73`；`ON-only` 面数 56→54；ON 合计 `mt 65331→65536`（OCCT 66640）。门禁：check 0 err、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、四模型默认 `--wirehist` 逐项不变。**下一入口（未动）**：`pf=27`（Torus，`pf=28` 的拆缝兄弟）仍 ON-only `192/301` vs OFF/OCCT `185/288`，与其余 ±1 节点的 ON-only 面（`pf=1/12/14/38/79 …` `72→73`）同类 —— 共享态拆缝面的参数域层，按 `BRepMesh_DefaultRangeSplitter.cxx` / `BRepMesh_NodeInsertionMeshAlgo.hxx:142-184` 对拍 `--fsurf` 的 UV 窗口与 `FSTAT_OCC` 字段。以下是历史记录。【T-103 下一步（§9.589，本轮实测）：配对键已解决 + 残差锁定 Torus。① 仪器坑：同一份 OCCT T0M 的 `--facestats` 与 `--wires` **面序不同**（1778 行按 box 全可映、1772 行换了索引、同索引只 6 行相同）⇒ 跨模式只能按几何配对；② 端口↔OCCT 按 bbox（容差 0.05）配对 `paired=1722 unmatched=55`、`1462/1722 (85%)` 逐位相同，残差全部落在 **Torus**（n=74、`+4122 mt`、ratio 1.264、仅 27 张逐位相同），BSpline 0.961 / Plane 0.986 / Cone 1.016 / Cylinder 1.002 / Sphere 1.003 ⇒ T0M 的 +6.0% ≈ Torus 的 +26.4%；③ `TorusRangeSplitter::generate_surface_nodes`（`splitter.rs:534-648`）已对 `D:\\source\\OCCT-src\\...\\BRepMesh_TorusRangeSplitter.cxx:22-115` 逐句审、判定忠实 ⇒ 缺口在**喂点**：OCCT `collectWirePoints`（`BRepMesh_NodeInsertionMeshAlgo.hxx:142-183`）每条边只喂 n−1 个点（丢共享顶点），端口正常路径（`node_insertion.rs:494-519` 喂 `wires_uv`）已一致，但回退路径（`node_insertion.rs:266-272`）喂全量 `boundary_uv`（`node_insertion.rs:483` = 含重复顶点）⇒ 抬高 `fillParams` 的 `aLength`、缩小 `aStdStep`、环面变密；先加 env 开关判该批 Torus 面走 1 还是 2，命中就按 `hxx:142-183` 让回退路径也喂 n−1 序列（不动 `uv`/`p3t` 前沿结构）。旁支：`step_obj_gates` 的 `step_obj_area_matches_occt` 单跑 FAILED 的原因已查明 —— **工作树里 `data/occ/occ-T0M.obj` / `occ-TDB.obj`（含 `.mtl`）被删了**（`git status` 的 worktree `D`，HEAD 里在），`common::occ_text` 读参考 OBJ 时 `panic!(\"missing occ/occ-T0M.obj\")`。`git restore --source=HEAD` 恢复后 4 个文件 `git diff` 为空（`core.autocrlf=true` 导致 `status` 仍显 `M`）；恢复后清掉残留插桩 env 重跑 `--test step_obj_gates` = **5 passed / 454.52s** ⇒ 那次 FAILED 纯属缺文件。】【原记录·T-103（§9.587，编号 D33）：T-102 已按主验收收口（T0M `--wirehist` = `TOTAL faces=1778 wires=1921` + `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`，与 `--wires` 真值逐项相同；a3n00/acs10/TDB 与 §9.586 逐项相同；`step_obj_gates` 5/5）。根因是 `load_wires.rs::wire_data_edges` 对 REVERSED wire 的边做了二次反转（`ShapeFix_ComposeShell.cxx:629` 的 `TopoDS_Iterator` 已复合一次），修掉后 §9.585-C 的两条残差同时消失（4 张环面碎片 1461/1502/1647/1667 都在、`surf=1300` 的 face=1758 回到 1 wire / 9 边）。本卡接手**网格侧**残差：(1) T0M 逐面三角数 `mesh_t=70557` vs OCCT `FSTAT_OCC` 求和 66576（+6.0%），(2) `STATMAP matched=1777 unmatched=1` 的那张面无 stats（0-based face=1764，薄面簇 `x∈[-23.73,-21.27] y∈[-206.86,-204.33] z∈[-107.16,-104.24]`），(3) 4 张多 wire 面的 `MULTI` 盒值比 OCCT 小约 9e-3（`wires` 数相同）。做法：先用既有探针在同一次运行内取 OCCT 真值（`occt_probe --facestats 2.574312 0.349066` / `--wires`）与端口（`zz_probe_a3n00 --fstats` / `--ecensus`），按 bbox 六坐标配对逐面对读 `mv/mt`，落在哪一层（UV 覆盖 / 参数域离散 / 三角化）就按该层 `.cxx` 修，不为 T0M 单独调偏转或加谓词。验收：T0M `--wirehist` 与 `--ecensus` 的 106 条 `MULTI` 不许回落，`step_obj_gates` 5/5（T0M 1.0053 向 1.0000 收敛且其余 22 个逐位不变），a3n00 ratio 1.0000 不许回落。】【原记录·T-102（§9.583，编号 D32）：缺口已用 `FixShape` 子模式 bisect 钉成两个独立子项 —— (1) **2 张 loop-wire 面**：端口 census 628（12 边）/ 635（16 边）现为 `wires=1`，OCCT 同 bbox 为 `wires=2`（`--fix FixLoopWiresMode 0` 会把 OCCT 也退回 1 条，逐 bbox 对得上）⇒ 按 `ShapeFix_Face.cxx:583`（`NeedFix(myFixLoopWiresMode) && FixLoopWire(aLoopWires)`，多条 wire 逐条 `B.Add(tmpFace, …)` `cxx:591-598`）移植 `FixLoopWire`（`cxx:2478`）+ `FindNext`（`cxx:2398`）+ `ShapeAnalysis_Wire::CheckLoop`（`ShapeAnalysis_Wire.cxx:2228-2340`），依赖端口 `ShapeFix_Wire` 的 `Analyzer()` 等价物与 `FixReorder`；(2) **3 张碎片面**：仍走 §9.582-E（pre-seam wire 首段 `ShapeFix_Face.cxx:365-480` / 端口 face 1752 多 1 张）。bisect 只能用 `occt_probe --fix <DE_ShapeFixParameters 字段> <int>`（`--set` 与读前的 `SetShapeFixParameters`/`SetShapeProcessFlags` 都是 no-op，见 D32）。验收：端口 `--wirehist` 目标 `1:1672 2:92 3:8 5:1 6:1 7:2 8:2` / 1778 面 1921 wires，`--ecensus`/`--fstats` 对 `output/t0m_fs_occt.txt`，`step_obj_gates` 5/5，a3n00 ratio 1.0000 不许回落；不得改 area_tol/断言、不得加 OCCT 里没有的谓词。】【原记录·T-102（§9.582，编号 D31）：已落地并验证 —— reader 侧 `resolve_face` 接受 `FixMissingSeam` 的多面（Shell）结果、`resolve_shell` 把它的面逐张加进 shell（`ShapeFix_Face.cxx:2266-2268` + `ShapeFix_Shape.cxx:257/294` + `ShapeBuild_ReShape.cxx:282-299`）⇒ T0M 1772→1775 面、wires 1915→1916、HIST `1:1671 2:90 3:8 5:1 6:1 7:2 8:2`，面积比 0.9809→**1.0084**（acs10 0.9928→1.0045，其余 21 个模型逐位不变，step_obj_gates 5/5）。下一步两条（都对着 .cxx、彼此独立）：(1) 4 张碎片面（OCCT id=0 的 1462/1503/1648/1668 无端口对应）—— 审 `ShapeFix_Face.cxx:365-480` 的 wire 首段（`ShapeFix_Wire::Perform`，`ShapeFix_Shape.cxx:200` 置 `ModifyTopologyMode=true`），端口 `shape_fix_face.rs::perform_fix_missing_seam` 已标该段 UNPORTED；(2) 多 1 张 —— 端口 face 1752（bbox=(-9.5,-9.5,-60.189)-(9.5,9.5,-44.644)，`wires=2 edges=[5,1]`）= OCCT face 1758（id=28215，1 面 1 wire），端口的 seam 步给了 `Shell faces=2`，走 `FixDummySeam`/wire 合并那条线查。跑基线：`--wirehist`（目标 `1:1672 2:92 …` / 1778 面 1921 wires）、`--ecensus`/`--fstats`、`step_obj_gates` 5/5、a3n00 1.0000 不许回落；不得改 area_tol/断言、不得加 OCCT 里没有的谓词。测点：OCCT `specs\\occt_probe\\run_here.bat <step> --faceids`（本轮 TEMP 加 bbox）/ `--wires` / `--nofix` / `--facestats 2.574312 0.349066`，端口 `zz_probe_a3n00 --wirehist` / `--ecensus` / `--fstats` / `--fixms`】下面是本轮落地前的计划原文（入口『改 CollectWires』已被 §9.582 推翻，仅留档）：【T-102 下一步（§9.581，编号 D30）：把 T0M 这一对相邻面（OCCT face 607 bbox=(-21,-259.653063,27.038323)-(21,-197.219056,106.702324) / face 628 bbox=(5.3,-217.417466,104.896128)-(15.5,-204.344246,120.102787)；端口 idx116 / idx351）的 wire 合并做成 OCCT 的 **1 条**，并把少的 **6 张单 wire 面**找回来。为什么是这两条：`--wires` 实测 OCCT `TOTAL faces=1778 wires=1921 ... WIREHIST 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`，端口 ④ 开 `1772 / 1915 / 1:1666 2:92 3:8 5:1 6:1 7:2 8:2` ⇒ 多 wire 桶逐项相同，差值全是单 wire 面（1778-1772=6，wires 同少 6）；`--ecensus` 得 idx351 = `wires=8 edges=[16,8,8,8,8,8,8,8]`（④ 关 = `wires=10 [1,8,8,8,8,8,8,8,8,1]`）、idx116 = `wires=2 edges=[26,45]`（两态相同），而 OCCT 两个面各 1 条（都不在 `--wires` 的 MULTI 清单里）⇒ 合并没做成。做法：按 `ShapeFix_ComposeShell.cxx` 的 `CollectWires`（候选判据/连接块/消费标记/合并入口）审控制流，先定位这 8/2 条为什么没并成 1 条（对照 `.cxx` 的收集与合并段，不按某一 STEP 调参），命中即按 .cxx 修，然后跑既有基线：`--wirehist`（期望 `1:1672 2:92 3:8 5:1 6:1 7:2 8:2` / 1778 面 1921 wires）、`--fstats` 对 `output/t0m_fs_occt.txt`（idx116/idx351 都出网格且逐位相同）、a3n00 `WIREHIST 1:208 2:9 4:1 6:4 10:4` + ratio 1.0000 不许回落、`step_obj_gates` 5/5、`--lib` 不劣化。纪律：本轮已把 ④ 判定为『变好』（38 个差异面里 OCCT 有同 bbox 的 21 个，④ 开 21/21 逐位相同，④ 关 0/21；且 ④ 不增删面：两态都 1772 面，只切换 idx116/idx351 谁被 Delaunay 拒绝）⇒ 不许为了 T0M 面积比撤 ④；也不许加 OCCT 里不存在的谓词/启发式。测点已就位：OCCT `specs\\occt_probe\\run_here.bat <step> --facestats 2.574312 0.349066` / `--wires`，端口 `zz_probe_a3n00 --fstats` / `--ecensus` / `--wirehist`】以下是 T-101 收口记录：【T-101 已收口（§9.579 + §9.580）：a3n00 ratio **1.0000**、`WIREHIST 1:208 2:9 4:1 6:4 10:4` 与 OCCT 逐项相同、F170 `wires=1 nEdges=8`；新增 ③ `wire_fix.rs::fix_dummy_seam` 复合朝向、④ `reshape.rs::MapReShape::apply_impl`（ShapeBuild_ReShape.cxx:200-324 的 EDGE 重建）；代价 T0M 0.9931→0.9809、acs10 0.9984→0.9928（tol 内，门禁 5/5）】以下是 §9.579 达成时的记录：【T-101 已达成（§9.579）】F113 恢复为 OCCT 的 1 面 / 2 wire（22+6）；a3n00 面积比 0.8996→0.9963、未网格面 1→0、faces 225→226；T0M 未网格 6→2（面积比 0.9931，tol 0.025 内）；acs10 0.9846→0.9984；--lib 1281/0、step_obj_gates 5/5。落地两处 .cxx 分歧：① `wire_fix.rs::copy_replace_vertices_with` 把 `EmptyCopied` 的朝向/位置落定移到两次 `B.Add` 之前（ShapeBuild_Edge.cxx:103→107-114，TopoDS_Builder.cxx:74-91，TopExp.cxx:214-253 + TopoDS_Iterator.cxx:26-83）——改前 REVERSED 边的首末顶点互换；② `make_faces_on_patch.rs` 三处 `Perform` 改传 `RecadreOnPeriodic=false`（cxx:3100/3110/3223）。两者必须同时在场。以下为达成前的收敛记录（§9.452-§9.578，goal 160 轮终态）：F113 空面锁到一个拓扑事实：93 段「首末顶点同一、逐边 LastVertex 是另一 TShape」⇒ IsShortSegment 恒 0 ⇒ shorts 全 0 ⇒ 合并循环空转 ⇒ 8+14 合不成 22 ⇒ [6,7,5,4,3,3] ⇒ 6 wire/5 面 ⇒ Shell(5) 被 reader 丢弃。链上每处能与 .cxx 对照的环节均已判定忠实：LoadWires（两条产 Internal 分支实测 0 次）、SplitByGrid（头部/closed 段/U 线循环/调用侧）、SplitByLine 内部（交点 1 个且在端点；守卫① dlast=0 精确命中 ⇒ 不切）、SplitWire（三建顶点分支 + 顶点统一机制 + copy_replace_vertices_with）、BreakWires（出口 154 段无一行 Internal）、CollectWires（候选判据/连接块/消费标记 :268↔:110/合并入口）、DispatchWires、FixDegenerated/CheckDegenerated、IsShortSegment（八项）、WireSegment 访问器。已排除：输入未统一（61/61 same=true）、MapReShape 缺 EDGE 深度重建（补上编译通过但 F113 不变 ⇒ 无收益已还原）、内部顶点未绑定（收集与 cxx:1006-1014 一致、该批边无非流形顶点 ⇒ 0 次触发是规格行为）、93 次 edge_last 失败（多边段上亦为规格行为）。唯一剩余入口：用既有 specs/occt_probe 在 SplitEdges/SplitByGrid 后 dump seqw 每段 {NbEdges, First/LastVertex 坐标, Orientation}，与端口 BREAK-IN 记录逐段比对（段集合是否同构）。命中即按 .cxx 修 ⇒ 跑 §9.565-C 五步验收。红线：§9.439 的 FixDummySeam 朝向感知修复会让面积比 0.8996→0.8918（仅 f=171），弄清 f=171 前不落地。⇒ 该入口的答案在 §9.579：真正的分歧不在 BREAK-IN 段集合的同构性，而在 SplitWire 造边时端点朝向语义（①）与 MakeFacesOnPatch 的周期重定位（②）。剩余旁支：WIREHIST 仍差 1 张面（端口 1:207+2:10 vs OCCT 1:208+2:9）；§9.439 的 FixDummySeam 修复前提已变，可重估。",
    why: "【T-103】T-102 的结构对齐已完成（T0M 1778 面 / 1921 wire / WIREHIST 与 OCCT 逐项相同；4 张环面碎片与 `surf=1300` 的多拆都归位），剩下的差异只在网格侧：同偏转下端口 `mesh_t=70557` 比 OCCT `FSTAT_OCC` 逐面求和 `66576` 多 6.0%，且有 1 张薄面不出网格（`unmatched=1`）——T0M 面积比 1.0053 的残差就落在这里，所以本卡按『同一次运行内取两侧真值 + 按 bbox 配对逐面对读』推进，而不是继续动拓扑。测点已就位（`--fstats`/`--ecensus` 与 `--facestats`/`--wires`），无需新写测试。以下是上一卡的历史记录：【T-102】用户要求「a3n00 一个一个处理，先处理带有倒角的法兰盘」，本轮就把那一步做完并达成 T-99 的 accept。① 做法是把 f1e56776 删掉的导入期 seam 步骤接回 reader（`ShapeFixFace::fix_missing_seam` + 结果面 wire 的 `check_pcurves_and_shift`），这不是新造规则：`ShapeProcess_OperLibrary.cxx:785-899` 的 FixShape 算子 → `ShapeFix_Face::Perform` → `FixMissingSeam`（`cxx:482-498`，构造在 `cxx:1722-2330`），`STEPControl_Controller.cxx:201/:221` 默认就开。② 移除理由（T-93(a)/§9.219「OCCT 在这条路径上不跑 Perform」）被本轮实测推翻：探针 `--faceids` 开/关 ShapeProcess 得 `{1:208,2:9,4:1,6:4,10:4}` vs `{1:163,2:53,4:2,6:4,10:4}`，后者与 STEP 的 FACE_BOUND 直方图、与端口逐字相同 ⇒ OCCT 的导入确实合并了 45 个 2-bound 面。③ 结果方向与量级都对：法兰面结构变成 OCCT 的 1 wire / 5 边（边序、参数区间逐项相同），a3n00 面积比 0.8627→0.8996、T0M 未网格 7→6、T0M 面积比 0.9786→0.9987、acs10 0.9025→0.9846、rescue 16→0，且**没有改任何断言或 area_tol**。④ 剩余问题都已量化成有界的下一步（还差 1 个 2-bound + 1 个 4-bound 面没合并；全形状多改 4 个面待查；法兰面密度 72 vs GT 52），所以 T-99 按 accept 结项，另立 T-100 按「面」继续推进，而不是再调全局参数。",
  },
  tasks: [
        { id: "T-99",
      lane: "mesh",
      title: "同结构喂 OCCT：18/18 逐字段完全相同 ⇒ Delaunay 层忠实；病灶是端口的面拓扑（2 wires vs OCCT 1 条 seam 闭合 wire）",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 8,
      actual: 12,
      startedAt: "2026-09-30",
      updatedAt: "2026-10-02",
      completedAt: "2026-09-30",
      blocker: "",
      goal: "把 a3n00 的面积比从当前真实值 **0.8627** 推向 1（缺口 31179.6）。**注意**：F113 在当前 HEAD **已经出网格**（`wires=4 mt=267`），所以「让 F113 出网格」不再是本卡目标；旧表述里的 0.8996 / 缺口 22802 / F113 不出网格来自本会话丢失的未提交改动，已按 §9.324 更正",
      next: "**【唯一入口：段集合（BREAK-IN）同构性，须在 OCCT 侧取证】** 端口侧能对照的环节已全部判定忠实（§9.577 清单），且「shorts 全 0 ⇒ 合并空转」「93 次 edge_last 失败」在相同段集合下均为规格行为 ⇒ 只剩「同一时点 OCCT 的 seqw 段集合是否与端口相同」。① 用 specs/occt_probe（本节此前用过）在 ShapeFix_ComposeShell::SplitEdges/SplitByGrid 切分后 dump 每段 {NbEdges, FirstVertex 坐标, LastVertex 坐标, Orientation}；② 与端口 §9.452 的 BREAK-IN 记录（外环 [6,16,1,1,1,1]）逐段比对：同构 ⇒ 差异在后续步骤（继续二分）；不同构 ⇒ 差异即在此 ⇒ 回 SplitByLine/SplitEdges 产出逐点对齐；③ 命中即按 .cxx 修 ⇒ 跑 §9.565-C 五步验收（cargo check → zz_seam_fix 113 → Face → --model 113 → wires=2(22+6) → --fstats mt≈228 → t101_verify.ps1 全绿 + 逐模型基线对照）。纪律（§9.570-D）：先确认 API → 单探针 + cargo build 断言 0 errors 才解读 → 用插入原文精确反向撤除 → 写明统计口径 → 无收益改动一律还原。",
      acceptance: "a3n00 面积比从 **0.8627** 上升且不劣化；T0M 未网格 7 不增加（`f1e56776` 已记 T0M 6→7 的代价）；--lib 与 step_obj_gates 5/5 不劣化；**不得改基线或 area_tol 来绕过**；不新写测试、不改断言",
      evidence: "**2026-09-30，十三轮**。**(1) 面已 100% 锁定**（§9.303.1）：STEP `#5375` = `CYLINDRICAL_SURFACE` **R=34**、高 68，解析面积 `2πrh = 29053.4` 与该区域缺口逐字吻合。**(2) 下游逐级实测忠实**：reader `FACE_BOUND 到 wire` **1:1**；读入路径**无 healing 阶段**；4 条 wire **不共享顶点 TShape**；4 条 wire 在 3D 上**各自闭合**（gap 精确 0）⇒ `CheckWire` 拒绝**正确**；范围守卫逐行忠实 `cxx:1781-1802`；`ComposeShell` 的 `load_wires` 入口即 `[6,14,1,1]`。**(3) 病根由实验确认**（§9.315）：`shape_fix_face.rs:453` 的 `surf.clone()` 让 `tmp_f` 与原面共用 `repr_key`，seam 写入落到原面槽位；摘掉 `tmp_f` 的曲面注册后 ⇒ `unmatched` **1 到 0**、F113 出网格 `mt=267`、全 226 面出网格，但面积比 **0.8996 到 0.8627** ⇒ 不能落地。**(4) 面积差成因（§9.316 + §9.319）**：逐三角对比显示两版**最大的三角完全同名同值**（不存在多出一张巨大错误面）；实验版 `v` 更多却 `f` 更少、总面积少 8377 ⇒ **少铺**；逐面边数差集显示差异落在 **43 个面上、每个少 2–3 条 seam 边、合计 −107**。原因是 `fix_missing_seam` 要**读回自己刚写下的 pcurve**（`check_pcurves_and_shift`/`curve_on_surface_range` 都按 `repr_key` 取），摘注册后写与读不在同一个键上 ⇒ 读不到 ⇒ seam 边加不上。**(5) 核心改动试并否决（§9.320，本轮）**：把 `repr_key` 改为 `face_key` ⇒ 三个指标与摘注册实验**逐字相同**（`stats` 226、`unmatched` 0、`mesh` 11101/11121、面积比 0.862724、`edges` 976 对基线 1083、同样 43 面退化）⇒ **单纯改键的归属无效**。**(6) 已排除 11 条假设**：CheckWire 口径、`v_range`、ComposeShell 内部拆分、返回 Shell/compound、顶点身份/邻接口径、两个退化边同母线、`heal_shape` 接在读入路径、D18 的 n1 到 n2、该写回可中性解耦、「每面的每条边都要有 pcurve」、以及**让 pcurve 的键跟随面身份**。**(7) 量化**：F113 单面区域占 +26952 面积缺口，全模型净缺口仅 22802.58",
      write: "crates/occt-topo/src/shhealing/shape_fix_face.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.303 §9.314 §9.315 §9.316 §9.318 · crates/occt-topo/examples/zz_pcurve_key_probe.rs · crates/occt-topo/examples/zz_share_probe.rs",
      note: "【§9.321】13 条假设已否决。本轮新增两条：第 12 条「共享键的危害经 swap_seam 生效」被插桩否决（读取路径 0 次调用，并因此**推翻 §9.296.3 的判断**）；第 13 条「共享键的第二个作用是让 check_wire 找到 pcurve」也被否决（两版早退都是 0）。范围已收窄到：差异不在这两处，而 unmatched 1→0 提示 F113 在两版里是不同的面对象 ⇒ 下一步量「哪些面被采纳与否」的差集，而不是再猜哪个函数",
      dependsOn: [],
    },
        { id: "T-59", lane: "mesh", title: "a3n00 面积缺口：配对表已出，缺口定到单面 F113（已拆为 T-99）",
       status: "completed\", priority: \"P1\", owner: \"agent\", progress: 100, estimate: 5, actual: 2, startedAt: \"2026-09-30",
      updatedAt: "2026-09-30\", completedAt: \"2026-09-30\", blocker: \"",
      goal: "按 bbox 六坐标把 a3n00 的 226 个端口面与 GT 面配对，逐面定位构成面积亏损的 face class",
      
      next: "本卡收口：配对表已出（226/226 全配上），缺口已定到**单个面 F113**，后续修复拆到 T-99。留下两件已量化的事实：① D27 seam 接线后 a3n00 比值 **0.8996**（判据 `|ratio-1| < area_tol=0.15`），比本卡起点 0.8627 高 0.037 ⇒ **§9.287.4 那行修复（当时 −0.0329）已不再是「应用后越界」的阻塞项**，但它本身仍是净负面（见 T-99 (d) 重测：−0.0092、且救不活 F113）⇒ **不要**再把它当主线。② 面类分布（确定配对的 face 数）：Plane 86 / Cylinder 62 / Cone 30 / Extrusion 15 / Torus 7；端口只有 1 面未网格（F113）。工具：`.target-gate/pair_a3n00.py`（配对）、`.target-gate/deficit_a3n00.py`（按 bbox 归因面积）",
       acceptance: "产出可复跑的配对表并定位亏损 face class（不要求修好）；不新写测试、不改断言",
      
      evidence: "**2026-09-30 实测**：`.target-gate/pair_a3n00.py .target-gate/a3n00_port.txt .target-gate/a3n00_gt.txt` ⇒ `port faces=226 gt faces=226`、`paired=226 port unpaired=0 gt unpaired=0`、`pairs with bbox diff > 1e-3: 26`。**wire 数不匹配只有 3 组**：`port wires=4 vs gt wires=2`（1 面，F113）、`port wires=2 vs gt wires=1`（2 面：F140→gt f=39 Cone、F170→gt f=69 Cylinder）⇒ a3n00 的 wire 结构基本已对齐，缺陷是**外科式**的而非系统性。**面积**：ours 204327.72 / occ 227130.30 = **0.8996**（净缺口 22802.58）；按 bbox 把缺口归因到各面后，**F113 单面所在区域就占 +26952.18**（`pad=0.01`；pad 0→2 变动 <2%）⇒ 超过全部净缺口，其余面合计是净盈余。**起点数字（2026-09-29，已被本行取代）**：本卡原记 a3n00 比值 0.8627、应用 §9.287.4 后 0.8298（−0.0329）；D27 seam 接线后这两个数字都已过期，**以 0.8996 为准**。宽 tol 与面类的登记出处：`crates/occt-topo/tests/common/mod.rs:91-93`",
       write: "\", ref: \"specs/_a3n00_gap_analysis.md §9.291 · §9.287 · tests/common/mod.rs:93 · .target-gate/pair_a3n00.py",
       note: "本卡原为 T-93 的前置（D20 补立）。分析类卡，**库源码零改动**；发现的 F113 缺陷已拆为 T-99，避免把两类缺陷混在一张卡里",
      
      dependsOn: [],
    },
        { id: "T-93", lane: "port-gap", title: "reader 缺带 Context 的 ShapeFix_Face::Perform（导入期 seam 修复）",
       status: "completed\", priority: \"P0\", owner: \"agent\", progress: 100, estimate: 4, actual: 6, startedAt: \"2026-09-28",
      updatedAt: "2026-09-30\", completedAt: \"2026-09-30\", blocker: \"",
      goal: "端口 F1691 一类「单闭合边」面的边界环补回 GT 的边数（GT 孪生面为 4 条），T0M 未网格面下降，且 a3n00 面积比与其余 22 模型不劣化",
      
      next: "本卡结项（D27）。留下的有界旁支，按优先级：① **f=1758（圆锥）** —— `FixMissingSeam` 只为球/BSpline/退化环面构造 seam（`ShapeFix_Face.cxx:1909-1975` 四个分支），圆锥四个都不匹配 ⇒ 走 `cxx:1972-1975` 的 `return false`；要在 OCCT 里找给圆锥补 seam 的分支（若有），不要自己加。② **211 个伪增量**（整形状 `perform_fix_missing_seam` 共 212 面变动，只有 f=1691 是真缺口）—— 门禁全绿说明未越界；若要收敛，先做端口 `check_wire`(`shape_fix_face.rs:33-80`) 与 `ShapeFix_Face.cxx:1652-1718` 的逐条对读，**不要**按「少改几个面」调参。③ **f=351(10 wire) / f=116(大圆柱)** 与 seam 无关，另因。④ a3n00 接线后 `unmatched` 0→1、`mesh_t` 11121→11941，是**新出现的差**，若要追从这里入手",
       acceptance: "端口 F1691 一类的环边数与 GT 孪生面一致；T0M 未网格面下降；--lib 与 step_obj_gates 不劣化；a3n00 面积比不劣化",
      
      evidence: "**本轮（2026-09-30）落地实测，只切一个 hunk**（`read_topology.rs` 的 `resolve_face` 末尾恢复 seam 调用 + `check_pcurves_and_shift`）：`--lib` **1281/0**（12.91s）；`step_obj_gates` **5/5**（588.92s）。对照表（同命令、切换该 hunk）：T0M 未网格面 **7→6**、wires **2121→1931**、多 wire 面 **309→120**、1-edge 闭合 wire 596→232、`WIREHIST 1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1` → **`1:1652 2:106 3:8 5:1 6:1 7:2 8:1 10:1`**；a3n00 `stats` 226→225、`unmatched` 0→1、`mesh_v/mesh_t` 11101/11121 → 10863/11941。**同构性证明**：`zz_seam_fix` 对 f=1691 产出 4 条边（底圆 `par=[π,3π]`、两条 seam `par=[1.080839000541168,1.570796326794897]`、`deg=true` 顶点边），bbox z∈[2.744,2.744]→[2.744,3.244]，与 OCCT 成形后（`--typefaces 3 4.25`）逐项相同；接线后 f=1691 出网格 `mt=143` 且 bbox 与 GT 1697 逐字相同，f=405/406/407/408 全部治好且 bbox 与 GT 73/74/75/76 逐字相同。**缺口是调用点而非实现**：临时计数器（已摘除）量到 `read_step_file` 期间 `SEAMCALLS = 0`；源码里「UNPORTED: seam 构造 cxx:1899-2330」的注释是过期的，那段早已移植，已更正。**原有证据（§9.287）仍有效但已非主线**：`fix_notched_edges` 传参错位（wire_fix.rs:3084 传 n1，cxx:4019 传 i）是真 bug，一行修复使 T0M 未网格 7→3 但 a3n00 面积比 0.8627→0.8298 越界，已按红线回退；本卡的 seam 接线**不需要**那处改动，也不重定任何断言。**来源**：D22 判据侧无罪 → D25/D26 定位到导入期 seam 构造 → D27 接线落地",
       write: "crates/occt-topo/src/step/read_topology.rs",
      
      ref: "specs/_a3n00_gap_analysis.md §9.287 · §9.288 · §9.289 · §9.290 · ShapeFix_Face.cxx:492-498 :1652-1718 :1899-2264 · ShapeProcess_OperLibrary.cxx:785-899 · STEPControl_Controller.cxx:201 :221",
      
      note: "工具（均已固化）：`zz_seam_fix`（单面修复前后 + 整形状过触发量化）、`zz_uv_feed`（--ids/--dump/--model/--surf/--fixall）、`occt_probe --typefaces/--advbox/--wireinv/--uv/--uvsum`。**未新写单元测试、未改任何既有断言、未改基线**",
      
      dependsOn: [],
    },
        { id: "T-69", lane: "mesh", title: "T0M 未网格面（6→7）：根因已定位为「边界环缺边」（读入层），非 mesh/判据缺陷",
       status: "completed\", priority: \"P1\", owner: \"agent\", progress: 100, estimate: 3, actual: 5, startedAt: \"2026-09-24",
      updatedAt: "2026-09-30\", completedAt: \"2026-09-30\", blocker: \"",
      goal: "查清 T0M 未网格面（7 个）的成因：判据缺陷、还是喂进判据的点集缺陷",
      
      next: "本卡收口，后续动作转 T-93（D23）：缺边的修补属 reader 的边界环构建，不在本卡的 mesh 范围内。本卡留下的可复跑仪器：`occt_probe --uv <face>` / `--uvsum`、`cargo run --example zz_uv_feed -- data/occ/T0M.stp --ids|--dump`、`.target-gate/pair_ids.py`（按 bbox 六坐标配对）",
       acceptance: "给出「未网格面成因」的可复核结论，并留下能复跑的对拍仪器（不要求修好）",
      
      evidence: "**结论（D22）**：不是 RangeSplitter/IsValid/updateRange 的判据缺陷 —— 三者与 OCCT 逐行同构，且用同一参数域条件枚举过口径 ②/③可证明 **OCCT 拿到端口这同一批点也会拒**（P1691 的 V 预调整区间 [1.08083900054116766,1.08083900054116877] 在 `updateRange` 里判 NO-OP ⇒ `computeLengthV==0` ⇒ `IsValid=false`）。分歧在**喂进去的点集**，而点集零宽是因为**这些面的边界环本身缺边**（读入层）。**逐面实测**：端口 F1691 = 1 wire / 1 edge / 37 点 / 常 V=1.080839001 / bbox Z 零厚（2.744–2.744），其真孪生 GT 面 409 = 1 wire / **4 edge**（底圆 37 点 + 两条 seam 参数域 [1.080839,1.570796] + 顶圆 2 点）/ bbox Z 2.744–3.244；GT 的 V 展布 0.489957 全部来自端口环里没有的 seam 边。同类：P1758→G566（GT dV 1.434376，端口 0）、P405/406/407/408→G73/74/75/76（GT dU 0.429994，端口 ~1e-13，各 4 边环 vs 端口 1 边）。**面数本身**：GT 1778 面 vs 端口 1772 面（同一 T0M.stp、同一调用序列）。**端到端**：球带区三角数 0 vs 143，总 verts 50123 vs 60548、tris 47864 vs 67366。**D15 的配对表已撤回**（D21）：D15 的 P405→G407 / P1691→G1697 / P1758→G1764 经查 GT 面 405/1697 的 bbox 与端口面相差 60–180，是错配；改用 bbox 六坐标后真孪生见上（P116→G607 的 bboxdiff=0.0000，是 D15 唯一配对的）。未网格面数口径：`zz_uv_feed --ids` 报 7 个面无 stat（F116/F405/F406/F407/F408/F1691/F1758），与 D16 警告的 1757/15 vs 1765/7 属不同调用序列，本条只用前者。",
       write: "",
      ref: "specs/_a3n00_gap_analysis.md §9.288 · BRepMesh_NodeInsertionMeshAlgo.hxx:142-184 · BRepMesh_DefaultRangeSplitter.cxx:202-236",
      
      note: "本卡为分析卡，**库源码零改动**（库内临时插桩已摘除，见 D24；端口探针留在 examples 非库路径）。没有新写任何单元测试，也没有改既有断言",
      
      dependsOn: [],
    },
        { id: "T-54", lane: "mesh", title: "单面回退 wireframe_face_triangulation：6 个门禁模型上从未被触达（保留，注释已更正）",
       status: "completed\", priority: \"P2\", owner: \"agent\", progress: 100, estimate: 4, actual: 2, startedAt: \"2026-09-20",
      updatedAt: "2026-09-30\", completedAt: \"2026-09-30\", blocker: \"",
      goal: "回退路径不再是 OCCT 里不存在的耳切/桥洞，而是约束 Delaunay（或按 OCCT 直接暴露主线报错）",
      
      next: "结项（§9.306）：**先只做判断**这一步做完了，结论是**两条路都不走**。① **换 Delaunay**：这是给一个**不可达**的分支实现更复杂的算法，纯增风险零收益；② **删回退**：`needs_fallback` 仍有一条真实触发路径（`:461-463` 的 `catch_unwind` panic 救援），按纪律无实测触发例不删，且删掉会移除一个未预期 panic 的安全网。**实际动作 = 唯一的库改动是纯注释**：把过期注释（原文「端口管线仍失败 T0M 的 169/1772 个面，2026-09-20 实测」）改成本轮实测结论 —— T0M 现在 **0** 个失败面、回退在全部门禁模型上不可达。留过期注释比删代码更危险（会让下一个人以为 169 个面还在失败）",
       acceptance: "确认回退有无实测触发例；据实测决定换 Delaunay / 删除 / 保留并标注；--lib 1281/0、step_obj_gates 5/5 不劣化",
      
      evidence: "**2026-09-30 实测（插桩 `OCCT_TOPO_TRACE_FALLBACK`，已摘除）**：在 `discret_root.rs` 打印①被 `MeshStatus::FAILURE`/`REUSED` 过滤的面、②真正进到 `wireframe_face_triangulation` 的面及其产出。**6 个模型全部 `reached_rescue=0`**：a3n00 `skipped=1`（face 113，`failure=true`）、T0M/acs10/TDB/bottom/top 均 `skipped=0 reached_rescue=0`。**结构原因**（`discret_root.rs`）：`:394` 已标 FAILURE/REUSED 的面直接 `continue`，不进 `pending`；`:449-454` 算法自判失败的面 `needs_fallback=false`；⇒ `needs_fallback=true` **只剩 `:461-463` 的 `catch_unwind` panic 这一条路**，而 6 个模型都没有面在该处 panic。**修正卡面表述**：a3n00 的 F113 拿不到三角**不是**「回退没救成」，而是它**已被上游标成 FAILURE、在 `:394` 就被跳过**，回退连机会都没有（这比原阻塞说明里「先被判 invalid discrete range 早期返回」更早一步）。**门禁**：occt-topo --lib **1281/0**；a3n00 基线 `faces=226 stats=225 mesh_v=10863 mesh_t=11941`、`matched=225 unmatched=1` 未变；`git diff discret_root.rs` = 13 增 7 删、**逐行全为 `//` 注释**",
       write: "crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs",
      
      ref: "specs/_a3n00_gap_analysis.md §9.306 · BRepMesh_BaseMeshAlgo.cxx:52-62 · crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs:394 :449-454 :461-463 :479-480 :549",
      
      note: "由 T-71 结项时经 D-201 重开；D22 解除了 T-69 的阻塞。本轮用 6 模型实测回答了「有没有触发例」，并把过期注释更正 —— 这是**注释级**结项，不是代码级",
      
      dependsOn: [],
    },
        { id: "T-94",
      lane: "port-gap",
      title: "共享键写回不是可解耦的副作用：实验证明它是承重的（消除后 F113 出网格，但面积比 0.8996 掉到 0.8627）",
      status: "completed",
      priority: "P3",
      owner: "agent",
      progress: 100,
      estimate: 3,
      actual: 5,
      startedAt: "2026-09-29",
      updatedAt: "2026-09-30",
      completedAt: "2026-09-30",
      blocker: "",
      goal: "让 seam 合并的收益（wire 结构）与它对共享 GeometryRegistry 的改写解耦：只对目标面生效",
      next: "结项（§9.314 + §9.315）—— **结论是「解耦」这个前提不成立**。做法：按 A 方案的最小实现（新增 `tgeometry.rs` 的 `unregister_face_surface`，只摘临时面自己的曲面注册；`shape_fix_face.rs:453` 后调用一次），使 `repr_key(tmp_f)` 回退为 `tmp_f` 自身键。**效果一（达成目标）**：写回被**完全**消除 —— `zz_share_probe` 从 T0M `21/2/8`、a3n00 `3/0/0` 变为两边都 **`0/0/0`**，§9.296.6 的构造级证明得到实验级确认。**效果二（F113 修好了）**：a3n00 `unmatched` **1 到 0**、F113 出网格 `mt=267`、全 226 面 `mt=0` 为 0 个、`mesh_v/mesh_t` 10863/11941 到 11101/11121 —— T-99 追了九轮的 F113 在这个实验下出网格。**效果三（否决）**：面积比 `0.899606` 到 **`0.862724`**、绝对差 22802.6 到 **31179.6**；判据 `|ratio-1| < 0.15` 下偏离从 0.1004（裕度 0.0496）恶化到 0.1373（裕度仅 0.0127）⇒ **不能落地**。**三条结论**：① 该写回**承重**，不是可自由摘除的副作用 —— 消除它就等于让 `fix_missing_seam` 的结果真正生效，从而**改变哪些面被返回**（`unmatched` 归零即其生效的证据）；② 因此设想中的「中性解耦」不存在，任何让它不写回原面的改动**必然**是行为改动，必须按实质修复评估而非中性重构；③ **C 方案（只改 `swap_seam`）的前提也因此存疑** —— 若它不改变采纳结果则消不掉 F113，若改变了则是行为改动，**需先实验才知道，不能假定**。**已回退**：改动 `git stash`（`stash@{0}`），工作区与实验前逐字一致、基线回到 `stats=225 unmatched=1`。**未跑 23 模型全量门禁**，因为 a3n00 的面积比在判据方向上变差、该改动已无落地资格。**后续交给 T-99**（其病根已由本实验确认）",
      acceptance: "T0M 的 wire 结构在不回写共享 registry 的前提下不劣化；F113 仍出网格；--lib 与 step_obj_gates 5/5 不劣化",
      evidence: "**写入点（2026-09-30）**：`ComposeShell` 共 6 处写 registry —— `dispatch_wires.rs:85,210-211,293-294`、`split_by_line.rs:623-624`、`split_wire.rs:443` 都作用于**新建** sub-edge；**只有 `wire_data.rs:65-66`（`swap_seam`）作用于既有 seam 边**。**机制（构造级证明，§9.296.6）**：`shape_fix_face.rs:453` 用 `surf.clone()` 建 `tmp_f`（新面、**同一曲面 Arc 指针**），`:622-623` 交给 `ComposeShell`；内部经 `collect_wires.rs:261`/`load_wires.rs:123` 把 `&tmp_f` 传给 `swap_seam`，`repr_key(tmp_f_key)` 解析到原面曲面的同一指针（`tgeometry.rs:524-531`）故写进原面槽位。**步骤① 实测（§9.314）**：`zz_share_probe` 升级为逐槽位指纹后 ⇒ T0M **`21/2/8`**（更正此前的 1 面/4 槽位：旧快照对「pcRange 相同但端点被平移」的槽位漏检），其中 face 1752 是两组两两互换（`swap_seam` 特征）、face 1537 是四槽位 u 各 `+2π`（`delta/2π=1.000000`，来自 `shape_fix_face.rs:506/513` 的 `tkey` 路径）；两面**都出网格**（`f=1537 mt=10`、`f=1752 mt=15`），T0M 全 1772 面 `mt=0` 为 0 个。**边界实验（§9.315）**：见上 `next` 的三条效果与未落地原因。**原有证据**：补偿开启 T0M wires 1931 且 F113 无网格；补偿移除 wires 2121 且 F113 出网格；D9 已否证「pcurve 归一化 pass 能救未网格面」；已否证「推迟 shift 写入」（实测与基线逐字相同，惰性，已回退）",
      write: "crates/occt-topo/src/shape_fix_compose_shell/wire_data.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.296 §9.296.6 §9.314 §9.315 · crates/occt-topo/examples/zz_share_probe.rs · crates/occt-topo/src/tgeometry.rs:524-531",
      note: "本卡由 T-93/T-69 的旁支立卡。**净结论：验证性的负面结果** —— 机制被实验证实，但「解耦」这个目标本身不成立（写回与修复采纳是同一件事的两面）。**不要**再按「中性解耦」去动它；相关的实质问题已移交 T-99",
      dependsOn: [],
    },
        { id: "T-95", lane: "port-gap", title: "Geom2dInt 曲线适配层与 TheProjPCurOfGInter：已实现且已接线（卡面描述过期）",
       status: "completed\", priority: \"P1\", owner: \"agent\", progress: 100, estimate: 8, actual: 2, startedAt: \"2026-09-28",
      updatedAt: "2026-09-30\", completedAt: \"2026-09-30\", blocker: \"",
      goal: "把 SplitByLine 依赖的 2D 求交内核按 .cxx 补齐，去掉 ComposeShell 主路径上的 UNPORTED",
      
      next: "结项（§9.297）：卡面点名的全部已实现且已接线，无需动代码。**真实余量已拆为 T-100**（`ginter.rs` 的两条分派臂）。另外 `curve_tool.rs` 剩两处 UNPORTED（`EpsX(C,Eps_XYZ)` 与 `Degree`）是「接口未搬但**无消费方**」，本轮复跑 grep 确认全树只有定义处，不构成功能缺口，保持原样即可",
       acceptance: "对应 UNPORTED 清除；--lib 1281/0 与 step_obj_gates 5/5 不劣化",
      
      evidence: "**逐项核对（2026-09-30，`crates/occt-geom2d/src/geom2d_int/`）**：`curve_tool.rs` 卡里点名的全在 —— `get_type`(`:31`, lxx:32-35，含 Line/Circle/Ellipse/Hyperbola/Parabola/Bezier/BSpline/Offset/Other 全分支)、`value`(`:84`)、`d0/d1/d2/d3`(`:89/:94/:99/:104`)、`dn`(`:109`)、`first/last_parameter`(`:114/:119`)、`nb_samples_curve`/`nb_samples_curve_range`(`:136/:141`)、`nb_intervals`/`intervals`/`get_interval`(`:147/:152/:157`)、`line/circle/ellipse/parabola/hyperbola`(`:59`-`:80`)、`is_composite`(`:54`)。`proj_p_cur.rs`(56 行) 两个入口都在：`find_parameter_range`(cxx:27-65) 与 `find_parameter`(cxx:67-79)，内部走 `curve_locator::locate_range` + `GenLocateExtPC`（= Extrema_GCurveLocator + Extrema_GenLocateExtPC）。**已接线**：`conic_curve.rs:13` 引入、`:296 proj_p_cur::find_parameter(...)`、`:305 proj_p_cur::find_parameter_range(...)`。**模块内真实剩余 UNPORTED（4 处）**：`curve_tool.rs:128` `EpsX(C,Eps_XYZ)`（`Curve2d` 无 `Resolution`）与 `:161` `Degree` —— 两处**无消费方**（本轮复跑 grep：全树仅 1 处命中，即定义本身）；`ginter.rs:174`（`IntCurve_IntConicConic_1.cxx:2236` 的 Line/Circle 臂）与 `:188`（`IntCurve_IntCurveCurveGen.gxx:339-…` 的 `typ1 != Line` 臂）—— **真实缺口，已拆为 T-100**。**这是本会话第二次遇到卡面描述与代码现状脱节**（第一次 T-96/§9.295：两条 gcpnts arm 其实早已移植）",
       write: "crates/occt-geom2d/src/geom2d_int/curve_tool.rs",
       ref: "specs/_a3n00_gap_analysis.md §9.297 · Geom2dInt_Geom2dCurveTool.lxx · Geom2dInt_TheProjPCurOfGInter_0.cxx:27-79",
      
      note: "结项理由：卡面点名的范围（curve_tool + proj_p_cur）已全部完成且已接线；真实余量在 `ginter.rs` 的分派臂，与卡名不符，故拆出 T-100 以免下一个接手者找错文件",
      
      dependsOn: [],
    },
        { id: "T-100", lane: "port-gap", title: "ginter 空分支已移除（通用臂本可胜任）；typ1 != Line 臂确认为本仓不可达",
       status: "completed\", priority: \"P2\", owner: \"agent\", progress: 100, estimate: 4, actual: 1, startedAt: \"2026-09-30",
      updatedAt: "2026-09-30\", completedAt: \"2026-09-30\", blocker: \"",
      goal: "按 .cxx 补齐 Geom2dInt_GInter 分派里剩下的两条臂，去掉这两处 UNPORTED",
      
      next: "结项（§9.305）。**臂一（Line/Conic）**：原代码在 `Circle|Ellipse|Parabola|Hyperbola` 上是个**空分支** —— 命中即什么都不做，静默返回零交点。而通用臂本来就能做：`perform_line` 经 `IntCurveIConicTool::from_lin2d` 调 `perform_iconic_tool(tool, d1, p_curve, d2, …)`，而 `IntCurveIConicTool` 已能表示每一种圆锥（`from_lin2d`/`from_circ2d`/`from_elips2d`/`from_parab2d`/`from_hypr2d`），`MyImpParTool` 的第二个参数是 `&dyn Curve2d` 不限定类型 ⇒ OCCT 的专用重载是**特化（性能）而非新能力**。改动 = 删掉空分支、与 `_` 合并。**臂二（`typ1 != Line`）**：唯一调用方 `split_by_line.rs:203-204` 的第一条曲线恒为 `Geom2dLine`（`:42`）⇒ `typ1` 恒为 Line ⇒ **本仓不可达**，据实标注、不照搬。**没有后续动作**",
       acceptance: "`ginter.rs` 的两处 UNPORTED 清掉或改为有依据的「无触发例」标注；occt-geom2d / occt-topo --lib 与 step_obj_gates 不劣化",
      
      evidence: "**2026-09-30 实测**。**臂一**：新探针 `crates/occt-geom2d/examples/zz_lineconic_probe.rs` 走**公开分派器** `Geom2dIntGInter::perform_with_d2`（与 `split_by_line.rs:204` 同一入口），构造已知解（直线 = x 轴，圆 = 圆心原点半径 5，域 = 闭区间 `[0,2π]`）⇒ 输出 `is_done=true`、`nb_points=3`、点为 `(-5,0)@t=π`、`(5,0)@t=2π`、`(5,0)@t=0`；去重后 `x values = [-5.0, 5.0]`、`EXPECT +/-5 => PASS`（第 3 点是闭域端点与起点的重合，去重即消）。**没有多解、没有错解。****臂二**：`grep` 全树只有 `split_by_line.rs:203` 一处构造 `Geom2dIntGInter`，其第一条曲线来自 `:42 let j_c2d = Geom2dLine::new(*line.position());`。**代码依据**：`int_conic_curve.rs` `perform_line` → `IntCurveIConicTool::from_lin2d` → `perform_iconic_tool`；`conic_curve.rs:322 perform` 的 `the_par_curve: &dyn Curve2d`；`occt-core/src/intcurve/iconic_tool.rs` 的 `from_circ2d`/`from_elips2d`/`from_parab2d`/`from_hypr2d`。**门禁**：occt-geom2d --lib **72/0**、occt-topo --lib **1281/0**；a3n00 基线 `faces=226 stats=225 mesh_v=10863 mesh_t=11941`、`matched=225 unmatched=1` —— 与改前**逐字相同**（预期：被删分支体为空，删它不改变任何已走路径）；T0M 的 `zz_share_probe` 同为 `21/1/4`。**剩余 UNPORTED**（4 处，均有据）：`curve_tool.rs:128`（`EpsX(C,Eps_XYZ)`，`Curve2d` 无 `Resolution`）与 `:161`（`Degree`）**无消费方**（grep 全树仅定义处）；`ginter.rs` 的两条说明性注释",
       write: "crates/occt-geom2d/src/geom2d_int/ginter.rs",
      
      ref: "specs/_a3n00_gap_analysis.md §9.305 · crates/occt-geom2d/examples/zz_lineconic_probe.rs · IntCurve_IntConicConic_1.cxx:2236 :2861 · IntCurve_IntCurveCurveGen.gxx:339- · occt-core/src/intcurve/iconic_tool.rs",
      
      note: "由 T-95 结项拆出。**核心发现：那个 `Circle|Ellipse|Parabola|Hyperbola` 分支体是空的** —— 不是「未实现的重载」，而是「命中就静默返回零交点」，比走通用臂更差。用一个 9 行的探针就证否了「需要专用重载」这个前提",
      
      dependsOn: [],
    },
        { id: "T-96", lane: "arch", title: "GeomTrimmedCurve 参数表示：评估结论 = 不重做（原前提已被实测否证）",
       status: "completed\", priority: \"P2\", owner: \"agent\", progress: 100, estimate: 4, actual: 1, startedAt: \"2026-09-30",
      updatedAt: "2026-09-30\", completedAt: \"2026-09-30\", blocker: \"",
      goal: "判定是否重做 TrimmedCurve 参数表示（保留 basis 区间），从而让 gcpnts 两条 Compute arm 可忠实移植",
      
      next: "结项：**不重做**（§9.295）。原前提「两条 arm 无法忠实表达」已被实测否证 —— ① `compute_type`（`gcpnts.rs:66`，cxx:26-65）与 `compute_abs_composite`（`gcpnts.rs:102`，cxx:96-158）**都已移植**，并通过 `gcpnts.rs:491-501` / `:579-585` 的 `match` 三分支接线，`LengthParametrized` 与 `AbsComposite` 都可到达；② `[0,1]` 是一次**自洽的重参数化**：`d1/d2/d3/eval_dn` 的 `s^N`（`trimmed.rs:45-79`）、`parameter_intervals` 的结点重映射（`:130-140`）、`resolution` 的 `/s`（`:147-153`）用**同一个** `s = last-first`，而这三处正是三条 arm 唯一读取的量 ⇒ 数值等价。**不要**改成 basis 参数：12 处构造 + ~8 处显式 `[0,1]` 换算，且消费方按 `[0,1]` 喂参数，改后**静默给出错误几何**（不编译失败）。`GeomTrimmedCurveBasis`（`trimmed.rs:169-249`）已经是 OCCT 语义的那一版，需要时按点使用（现仅 `sphere.rs:47` 用）",
       acceptance: "产出一页评估（改动面 / 消费方清单 / 回归风险）并在本板登记结论；不做代码改动",
      
      evidence: "**2026-09-30 实测**：`gcpnts.rs` 的 `compute_type` / `compute_abs_composite` 已实现且被 `:491-501`、`:579-585` 两处 `match` 调用；`trimmed.rs` 的 `d1/d2/d3/eval_dn`(`:45-79`)、`parameter_intervals`(`:130-140`)、`resolution`(`:147-153`) 三处使用同一个 `s`。OCCT 侧核对：`Geom_TrimmedCurve.cxx:212-243`（EvalD0..DN 直接委派 basis）、`:255-267`（First/LastParameter = uTrim1/uTrim2）；`GCPnts_AbscissaPoint.cxx:26-65`（computeType）与 `:87-158`（两 arm）。**消费方清单**：构造 12 处（`rectangular_trimmed.rs:262,272`、`brep_builder_api.rs:229,246`、`geom_int_intss_make.rs:167,269,660`、`edge_split.rs:53`、`int_tools_lines.rs:131,132`、`step/read_topology.rs:1630,1668,1670`）；显式依赖 `[0,1]` 的注释/换算 ~8 处（`inttools/intersections.rs:18,38`、`edge_edge/edge_edge.rs:307`、`edge_split.rs:8`、`brep_builder_full/edge_builder.rs:107`、`extrema_surf/elementary_curve_extrema.rs:155`）。**本卡库源码改动仅一处注释**（`gcpnts.rs` 模块头，把过期的 UNPORTED 段改为实测结论）",
       write: "crates/occt-geom/src/trimmed.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.295 · Geom_TrimmedCurve.cxx:212-243 :255-267 · GCPnts_AbscissaPoint.cxx:26-65 :87-158 · occt-geom/src/gcpnts.rs:66 :102 :491 :579 · occt-geom/src/trimmed.rs:45,130,147",
       note: "评估类任务，交付即结论：不重做参数表示。T-51 的豁免理由（「余下两条要先改参数表示」）据此不成立，T-51 已重开核对其真实余量",
      
      dependsOn: [],
    },
        { id: "T-98", lane: "mesh",
      title: "锥面 pcurve 报 (-inf,inf)：OCCT 语义对照完成 —— 与 Geom2d_Line 一致，端口正确，不修",
      status: "completed",
      priority: "P3",
      owner: "agent",
      progress: 100,
      estimate: 1,
      actual: 1,
      startedAt: "2026-09-29",
      updatedAt: "2026-09-30",
      completedAt: "2026-09-30",
      blocker: "",
      goal: "pcurve 自报的参数域与实际可求值域一致，避免诊断/调用方按 (-inf,inf) 取样",
      next: "结项（§9.313）：走卡面第一条验收路径 —— **给出 OCCT 语义对照的结论并保持现状**，未改 mock/门禁。三条 OCCT 原文：① `Geom_ConicalSurface.cxx:207-215` 的 `Bounds` 就是 `V1=-Infinite, V2=Infinite`（锥面 V 本来就是无限域）；② `Geom2d_Line.cxx:142-150` 的 `FirstParameter()/LastParameter()` 返回 `∓Precision::Infinite()`（`hxx:97-100` 文档也明写 RealFirst/RealLast）；③ `GeomProjLib.cxx:81` 在 `GeomAbs_Line` 分支**原样返回未修剪的 `new Geom2d_Line(Proj.Line())`**，只有输入本身是 `Geom_TrimmedCurve` 才修剪（`:118-128`）。锥面母线（V 恒定）投到 UV 就是常 V 直线 ⇒ OCCT 同样报 `(-inf,inf)`。端口 `crates/occt-geom2d/src/line.rs:37-38` 逐字一致（`NEG_INFINITY`/`INFINITY`），`make_pcurve` 如实读取其 `first/last_parameter` ⇒ **输出与 OCCT 逐字相同，不修**。**真正的教训在探针侧**：上一轮把 d0 读成「全 NaN」就是**用自报的 -inf 端点求值**造成的（D14 已撤回）；端口本身早有 `is_finite` 守卫（`make_pcurve.rs:21`/`:64`/`:219`/`:229`），正确口径是 pcurve 域无限时改用 **3D 边的参数域**（本实例 `[0, 2π]`），与 OCCT 自己的防御写法同思路",
      acceptance: "（可选卡）要么给出 OCCT 语义对照的结论并保持现状，要么让自报参数域与可行域一致；两者都不得改动 mock/门禁",
      evidence: "**2026-09-30 OCCT 语义对照（§9.313）**。**原有实测（2026-09-29，保留）**：F1758 `type=Cone`、`surfU=(0,2π)`、`surfV=(-inf,inf)`、单条闭合边 `edgeRange=(0,2π)`；`make_pcurve_full` 返回 `pcRange=(-inf,inf)`，但 `d0` 在 t=0..2π 九点全为有限值（U 由 0 递增到 2π、V 恒为 0.717187881），且 `|3D(t) − S(uv(t))|` = 3.043e-12（t=0）与 2.696e-12（中点）⇒ **pcurve 功能正确**；上一轮记的「d0 全 NaN」是旧探针用 -inf 端点求值造成的假象，已撤回（D14）。**本轮新增的语义对照**：见 `next` 的三条 OCCT 行号 —— 结论是 `(-inf,inf)` 正是 OCCT 的行为，端口逐字一致。**顺带修掉一个真实文档错位**：`make_pcurve.rs` 里 `make_pcurve_full` 的文档块被粘到了 `select_forward_seam` 的文档里，导致前者**反而没有文档**、后者的文档混进不属于它的内容；已把两段拆开，并把「The returned pcurve is parameterized over the edge's range `[a,b]`」这句**不准确**的描述改为准确描述（返回曲线带自己的域，投影落在 `Geom2d_Line` 时是 `(-inf,inf)`），附三条 OCCT 行号 + 对调用方的取样警告。**纯注释改动，19 增 6 删**。**门禁**：`cargo check` 0 error；`occt-topo --lib` **1255 passed / 26 failed** 与改动前相同（26 条为 `%TEMP%` 环境性失败）；a3n00 `faces=226 … mesh_v=10863 mesh_t=11941`、`matched=225 unmatched=1` 未变；T0M `21/1/4` 未变",
      write: "crates/occt-topo/src/pcurve_full/make_pcurve.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.313 · Geom_ConicalSurface.cxx:207-215 · Geom2d_Line.cxx:142-150 · Geom2d_Line.hxx:97-100 · GeomProjLib.cxx:81 :118-128 · crates/occt-geom2d/src/line.rs:37-38",
      note: "本卡已从 P1「NaN 缺陷」降级为 P3「自报参数域不一致」，现在**连这个不一致也不成立** —— 它是 OCCT 的行为。不要再以「NaN」为名动它，也不要在库里 clamp 参数域",
      dependsOn: [],
    },
        { id: "T-97", lane: "hygiene",
      title: "「never used」定性为「移植已做、消费方未移植」；A 类 60 个模块已标注，死代码 390 到 111",
      status: "completed",
      priority: "P3",
      owner: "agent",
      progress: 100,
      estimate: 3,
      actual: 3,
      startedAt: "2026-09-30",
      updatedAt: "2026-09-30",
      completedAt: "2026-09-30",
      blocker: "",
      goal: "把确认无消费方的死代码按模块成批删除，而不是继续逐符号清 import",
      next: "结项（§9.311）。**方向已按用户选择改为「标注而非删除」**，且已落地。**定性**：这份代码库在逐模块移植 OCCT，OCCT 的调用方往往还没移植，于是被移植的**被调用方**表现为「无人使用」⇒ **主体是「忠实移植、消费方未移植」（A 类），不是死代码**，不该删（§9.309/§9.310 的两次误删已全部回退）。**已做**：对 60 个 A 类文件（模块头自报 OCCT 来源者）在其 `//!` 头后插入统一 6 行块（5 行「待接线」说明 + `#![allow(dead_code)]`），每个文件 diff **恰好 +6/−0**。**效果**：死代码 **390 到 111**、总警告 **513 到 220**（−57%），**未删任何已完成的工作**。**若日后要继续**：剩余 111 条死代码分布在 26 个文件，多为 **B 类**（无 OCCT 来源的本地实现）或**函数内嵌套项**；这些才是真正该逐个看的候选，**但每条仍须先核 OCCT 侧有无调用方再决定**（不要重犯第 1/2 批的错误）",
      acceptance: "产出一页评估（改动面 / 消费方清单 / 回归风险）并在本板登记结论；不做代码改动",
      evidence: "**2026-09-30，两批误删全部回退，改为标注交付**。**口径更正**：此前引用的「561 条」是两个来源混合的产物（`rtk` 是输出过滤代理，会把 occt-core/geom/geom2d 与本 crate 混在一起且不重报缓存命中 crate）；改用**直接 cargo + 绕过过滤器**后，occt-topo 自身为 **513 条**。**分类判据**：读模块 `//!` 头前 35 行，出现 `SomeOcctClass.cxx|hxx|lxx|gxx` ⇒ A 类（忠实移植、消费方未移植）；否则 B 类。**注意**：`value assigned … is never read`（局部赋值未读，35 条）必须单列，不能算死代码——上一轮我的正则 `is never read` 误捕了它。**落地**：60 个 A 类文件各插入统一 6 行块（`//! T-97: items below are faithful ports… //! …wire the consumer instead.` + `#![allow(dead_code)]`），**每个文件 diff 恰好 +6/−0**（已逐个校验）。**效果**：死代码 390 到 **111**、总警告 513 到 **220**、未用 import 47 与局部赋值未读 35 均不变（说明标注只压死代码、未掩盖别的类别）。**过程返工**：第一版标注用了两种措辞，规范化脚本只匹配其中一种而留下残留行；发现后未继续打补丁，而是把 75 个文件整体 `git checkout` 回 HEAD 再干净重标。**门禁**：`cargo check` 0 error；a3n00 `faces=226 … mesh_v=10863 mesh_t=11941`、`matched=225 unmatched=1` 逐字未变；T0M `21/1/4` 未变；`--lib` **1255 passed / 26 failed** 与改动前完全相同（26 条为 `%TEMP%` 写文件被拒的环境性失败，§9.309.7 已用对照实验证明与代码无关）；`git diff --shortstat` = 65 files / 492 insertions / 55 deletions，其中 60 个是 +6/−0 标注、其余 5 个是既有功能改动。`#![allow(dead_code)]` 只影响 lint、不生成代码，故基线不变属预期",
      write: "crates/occt-topo/src/lib.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.307 §9.308 §9.309 §9.310 §9.311 · .target-gate/final_census.py · .target-gate/annotate_clean.py · .target-gate/raw_warn.txt",
      note: "**本卡的产出是「判据 + 降噪」，不是删除量**（净删除 0 行）。两次误删（第 1 批删 31 函数、第 2 批删整模块 503 行）均已逐字回退，教训是：`never used` 绝不等于「没用」；判断死代码只能靠「查 OCCT 侧调用方」或「读模块头自报来源 + 判断其消费方是否已移植」，不能靠名字 grep、`pub use` 清单或编译级可达性证明",
      dependsOn: [],
    },
        { id: "T-92", lane: "bop", title: "ComposeShell / FixMissingSeam 逐行迁移", status: "completed", priority: "P1", owner: "agent",
      progress: 100, estimate: 4, actual: 4, startedAt: "2026-09-25", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "ShapeFix_ComposeShell 与 ShapeFix_Face 的 seam 补齐按 .cxx 逐行移植",
       next: "—\", acceptance: \"--lib 1281/0、step_obj_parity 14/14、export 逐位一致",
       evidence: "a3n00 wire 直方图 2-wire 面 24→10；mesh 11675/13168",
       write: "crates/occt-topo/src/shape_fix_compose_shell",
       ref: "ShapeFix_ComposeShell.cxx · ShapeFix_Face.cxx:1722-2330",
       note: "主路径上仍挂 T-95 的 UNPORTED",
      dependsOn: [],
    },
        { id: "T-82", lane: "bop", title: "BOPAlgo_BOP::BuildSolid 曲面布尔精确化", status: "completed", priority: "P0", owner: "agent",
      progress: 100, estimate: 3, actual: 3, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "曲面布尔走精确 BOPAlgo 路径而非采样近似\", next: \"—",
      acceptance: "--lib 不劣化；box∪cyl 单实体、体积≈8.5；phase10 / boss 不劣化",
      
      evidence: "HEAD 08c4edf2 实测 box∪cyl FUSE = 1 solid / 8 面 / volume 8.502655，逐面面积与 GT（occt_probe --fuse）逐一相同；卡片的 8.737 与 loops=0 是更早快照",
       write: "crates/occt-topo/src/bop_bop.rs\", ref: \"BOPAlgo_BOP.cxx:583-711 · BOPAlgo_Builder_3.cxx::BuildSplitSolids",
       note: "AC 未覆盖的其他曲面布尔组合未复核",
      dependsOn: [],
    },
        { id: "T-88", lane: "bop", title: "phase10 curved_face_fillet_sphere_plane", status: "completed", priority: "P0", owner: "agent",
      progress: 100, estimate: 1, actual: 1, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "phase10 的球面/平面倒圆用例转绿\", next: \"—\", acceptance: \"--test phase10_integration 8/8",
      evidence: "HEAD 08c4edf2 实测 phase10 8/8；旧板记的 7/8 是 2026-09-21 快照，原诊断未再复现",
       write: "crates/occt-topo/src/bop_bop.rs\", ref: \"BOPAlgo_BuilderSolid · BOPTools_AlgoTools\", note: \"—",
      dependsOn: [],
    },
        { id: "T-41", lane: "bop", title: "摘除 bop_curved 的体素/网格布尔（OCCT 无对应）",
       status: "completed\", priority: \"P1\", owner: \"agent\", progress: 100, estimate: 2, actual: 3, startedAt: \"2026-09-20",
      updatedAt: "2026-09-28\", completedAt: \"2026-09-28\", blocker: \"\", goal: \"布尔入口改派到忠实路径，去掉 OCCT 里不存在的体素/采样布尔",
       next: "—\", acceptance: \"摘除后 --lib 与门禁不劣化",
      evidence: "改派后实测 --lib 1281/0、step_obj_gates 5/5、phase10/19/3/4 = 8/5/4/9；调用图 curved_boolean_full ← boolean_dispatch:227 / bop_builder_repair:239,405 / bop_builder_report:388，按 (a) 夹具重设计",
       write: "crates/occt-topo/src/bop_curved/region_trim.rs",
       ref: "—\", note: \"classify_face_general 本身仍是 8×8 采样，不是 BRepClass3d_SolidClassifier",
      
      dependsOn: [],
    },
        { id: "T-67", lane: "port-gap", title: "Extrema_ExtPExtS / ExtPRevS 两臂", status: "completed", priority: "P2", owner: "agent",
      progress: 100, estimate: 3, actual: 3, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "补上点-曲面极值搜索的两条 OCCT 分支，摘掉替换实现",
       next: "—\", acceptance: \"occt-geom --lib 143/0 + 门禁不劣化",
      evidence: "两臂已移植并分派（point_surface_extrema.rs:60-61/:88-91、:450 perform_ext_ps、:485 perform_rev_ps，对应 Extrema_ExtPS.cxx:292-343）；A15 替身已摘除。实测 occt-geom 143/0、occt-topo 1281/0、step_obj_gates 5/5",
       write: "crates/occt-geom/src/extrema_surf\", ref: \"Extrema_ExtPS.cxx:292-343\", note: \"—",
      dependsOn: [],
    },
        { id: "T-25", lane: "arch", title: "GeometryRegistry 侧表 → 几何进 TShape", status: "completed", priority: "P2", owner: "agent",
      progress: 100, estimate: 6, actual: 6, startedAt: "2026-09-20", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "几何不再走全局侧表，而是挂在 TShape 上\", next: \"—\", acceptance: \"--lib + 门禁不劣化",
      evidence: "三张几何表已删除（tgeometry.rs:163-181 注释写明 all live on their own TShape now），35 处读写转发到 TShape 槽；设计文档第 5 步（global() 兼容外壳）未做，属文档明示可停点。实测 occt-topo 1281/0、step_obj_gates 5/5",
       write: "specs/_design_architecture_t25_t28.md\", ref: \"specs/_design_architecture_t25_t28.md",
      note: "保留 global() 兼容外壳是既定停点，不是未完成项",
      
      dependsOn: [],
    },
        { id: "T-28", lane: "arch", title: "ImpPrm HVertex 合并 + intana/intpatch 重叠合并",
       status: "completed\", priority: \"P2\", owner: \"agent\", progress: 100, estimate: 8, actual: 8, startedAt: \"2026-09-20",
      updatedAt: "2026-09-28\", completedAt: \"2026-09-28\", blocker: \"\", goal: \"ImpPrm 的 HVertex 合并与 intana 的 closed form 单一化",
       next: "—\", acceptance: \"--lib + 门禁不劣化",
      evidence: "步 1–3（ImpPrm HVertex 合并）与步 4（intana/intpatch closed form 合并，五副本清零）均完成；步 5 经实测判为有害不做（D5）。实测 occt-topo 1281/0、step_obj_gates 5/5",
       write: "specs/_design_architecture_t25_t28.md\", ref: \"IntPatch_ImpPrmIntersection.cxx:221-469",
      note: "若日后重开步 5，前提是先把 PatchIntersection 对「一般 B 样条 × 球」修正到与 tracer 等价",
      
      dependsOn: [],
    },
        { id: "T-29", lane: "hygiene", title: "过期规格刷新", status: "completed", priority: "P3", owner: "agent", progress: 100, estimate: 1,
      actual: 1, startedAt: "2026-09-21", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "规格文件不再指向已删除的 _board.md，且过期内容有横幅",
       next: "—\", acceptance: \"—",
      evidence: "两个目标文件已带「仅历史（T-29 登记）」横幅（_coverage.md:1、_brepmesh_align_review.md:3）；指向已改为 specs/board.canvas.tsx，并同步 4 处悬空引用。occt-core check 通过、occt-topo 1281/0",
       write: "specs/_coverage.md\", ref: \"—\", note: \"—",
      dependsOn: [],
    },
        { id: "T-11", lane: "hygiene", title: "occt-topo 编译警告清理", status: "cancelled", priority: "P3", owner: "agent", progress: 0,
      estimate: 1, actual: 4, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "", blocker: "",
      goal: "清理 unused import 类警告\", next: \"已豁免：剩余项经实验判定不可安全清除；死代码部分改由 T-97 承接",
       acceptance: "无行为变化",
      evidence: "四批共清 11 处真实 unused import，每批后 --lib 1281/0；剩余 ~60 条多为「只被同文件 #[cfg(test)] 使用」的假阳性，实测删除后 test 目标报 cannot find type ⇒ 已回退；两条自动化路径亦实测失败",
       write: "crates/occt-topo/src\", ref: \"—\", note: \"豁免结论见 D7",
      dependsOn: [],
    },
        { id: "T-51", lane: "port-gap",
      title: "gcpnts 余量核对完成：UNPORTED 清单已对齐代码；真实余量是 1 个缺口（adv 版 AbsComposite）",
      status: "completed",
      priority: "P3",
      owner: "agent",
      progress: 100,
      estimate: 2,
      actual: 2,
      startedAt: "2026-09-26",
      updatedAt: "2026-09-30",
      completedAt: "2026-09-30",
      blocker: "",
      goal: "核对 gcpnts 的真实余量：两条 Compute arm 是否真的已完整，以及 `uniform_abscissa` 的外形差异是否该消掉",
      next: "结项（§9.312）。两条结论：① **两条 Compute arm 都已实现** —— `gcpnts.rs` 里 `:488`（旧行号）声称 `GCPnts_LengthParametrized`/`GCPnts_AbsComposite` 未移植的文档是**过期**的，同函数正文的 `match` 三分支全在；`abscissa_point_with_tolerance` 里「`:187-295` UNPORTED; use the non-adv walk」也与紧接的代码自相矛盾。两处已改为实测描述。② **「`uniform_abscissa` 外形不同」不成立** —— OCCT 的 `NbPoints` 重载本就是「含两端的点数」，`GCPnts_UniformAbscissa.cxx:510` 用 `aL / (theNbPoints - 1)` 作步长，而端口 `uniform_abscissa(c, n)` 用 `total / n` 产 `n+1` 个 ⇒ **`n = NbPoints - 1` 时逐字等价**，只是入参命名口径不同（段数 vs 点数）。**gcpnts 的真实余量因此是 1 个缺口而非 2 条 arm**：adv 版 `AbsComposite`（`CPnts_AbscissaPoint.cxx:187-295`，含 `anIndex == 0` 特例），它需要 `CPnts_MyRootFunction` 的容差重载 `Init(X0, L, Tol)`（`:298` 记的就是这一条）。现状：多段曲线退回非 adv 走法，**仅积分容差不同**。若要补，落点与行号已在 §9.312.5 写明",
      acceptance: "gcpnts 的 UNPORTED 清单与代码现状一致（要么清空、要么逐条列出真实缺口）；occt-geom / occt-topo --lib 不劣化",
      evidence: "**2026-09-30 核对（§9.312）**。**改动前 5 处 UNPORTED，逐条判定：2 处真实、2 处过期、1 处指向性引用**。真实项：`:46` 模块头（`GCPnts_UniformDeflection`/`QuasiUniformDeflection` 确未移植，本模块不再暴露假 accessor）、`:294`/`:563`（容差重载 `Init(X0, L, Tol)` 未由公开 API 暴露 + `AdvCompute` 的 `AbsComposite` 臂未移植，`CPnts_AbscissaPoint.cxx:187-295`/`:436-474`，附行为差异说明）。过期项：`:488`（同一函数的 `match` 三分支全在，且它引用的 module header 已在 T-96 改过）、`:593`（紧接的分支就是调 `compute_abs_composite`）。**`uniform_abscissa` 的关键证据**：`GCPnts_UniformAbscissa.cxx:495-497` 断言 `theNbPoints > 1`，`:510` `anAbscissa = aL / (theNbPoints - 1)`，`:515` `aSize = theNbPoints + 5` ⇒ OCCT 以「点数 - 1」为段数。**改动**：`gcpnts.rs` 三处注释（模块头 `:42-44` 的错误断言、`:488-489` 的过期 UNPORTED、`:593` 的内联 UNPORTED），**纯注释、零行为变更**。**门禁**：`occt-geom --lib` **143/0**；`occt-topo --lib` **1255 passed / 26 failed** 与改动前相同（26 条为 `%TEMP%` 环境性失败，§9.309.7 已证明与代码无关）；`cargo check` 0 error；a3n00 `faces=226 … mesh_v=10863 mesh_t=11941`、`matched=225 unmatched=1` 未变",
      write: "crates/occt-geom/src/gcpnts.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.295 §9.312 · GCPnts_UniformAbscissa.cxx:494-515 :510 · GCPnts_AbscissaPoint.cxx:26-65 :87-90 :96-158 :187-295 :436-474 · CPnts_MyRootFunction.cxx:32-37",
      note: "由 T-96 的结论重开。**结论是「余量比卡面说的小」**：两条 arm 早已实现、外形差异不成立，真实余量只剩 adv 版 AbsComposite 一个缺口。**不要**因此去改 `GeomTrimmedCurve` 的参数表示（§9.295.4 已判为高风险且无收益）",
      dependsOn: [],
    },
        { id: "T-03", lane: "bop", title: "boss_single_disc_base_merges_one_solid（门禁 1/2）",
       status: "cancelled\", priority: \"P3\", owner: \"human\", progress: 0, estimate: 1, actual: 1, startedAt: \"2026-09-19",
      updatedAt: "2026-09-21\", completedAt: \"\", blocker: \"\", goal: \"该夹具是否应通过",
      next: "维持 1/2 作为已知，不再投入（夹具不合法，非引擎缺口）",
       acceptance: "—\", evidence: \"夹具期望的「1 solid」在 OCCT 的 BOPAlgo 下也不成立；门禁 G-boss 记为已接受的红",
       write: "crates/occt-topo/tests/bop_builder2_boss.rs",
       ref: "BOPAlgo_BOP.cxx:583-711\", note: \"改动已按 (a) 夹具重设计落地（见 T-41）",
      dependsOn: [],
    },
        { id: "T-23", lane: "hygiene", title: "ATU01038 顶点/面密度 −1.8%", status: "cancelled", priority: "P3", owner: "human", progress: 0,
      estimate: 1, actual: 0, startedAt: "2026-09-21", updatedAt: "2026-09-21", completedAt: "", blocker: "", goal: "密度差异是否要追",
      next: "不改，作为已知记录\", acceptance: \"—\", evidence: \"UV-grid 与 deflection-adaptive 的采样差异；parity 已明确不断言密度",
       write: "crates/occt-topo/tests/common/mod.rs\", ref: \"—\", note: \"—",
      dependsOn: [],
    },
        { id: "T-101",
      lane: "mesh",
      title: "F113（螺帽顶部圆锥斜切面）恢复为 OCCT 的 1 面 / 2 wire（22+6）并出网格；a3n00 面积比 0.8996 → 0.9963、未网格面 1 → 0",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 8,
      actual: 1,
      startedAt: "2026-09-30",
      updatedAt: "2026-10-02",
      completedAt: "2026-10-02",
      blocker: "",
      goal: "让 a3n00 的 F113（螺帽顶部圆锥斜切面）恢复为 OCCT 的 1 面 / 2 wire 并出网格，a3n00 面积比从 0.8996 上升且不劣化其它基线 ⇒ **已达成**：F113 = 22+6 边、mv=227/mt=227；ratio 0.9963；unmatched 1→0；T0M 未网格 6→2；门禁全绿",
      next: "**【收口（§9.580）：a3n00 ratio 0.9963→**1.0000**、`WIREHIST 1:208 2:9 4:1 6:4 10:4`（与 OCCT 逐项相同）、F170 `wires=1 nEdges=8`；本轮新增 ③ `wire_fix.rs::fix_dummy_seam` 复合朝向（first_vertex/last_vertex）+ ④ `reshape.rs::MapReShape::apply_impl` EDGE 重建。代价 A/B：④ 使 T0M 0.9931→0.9809、acs10 0.9984→0.9928（均在 area_tol 内，门禁 5/5 不劣化）；T0M 有 30+ 张薄面 `x∈[5,16] y∈[-217,-204] z∈[104,117]` 的 `mv/mt` 由 9/7 塌成 4/2 ⇒ 这是下一个对拍点（在 UV/网格侧，不在 ④ 本身）】**以下为 §9.579 达成时的记录：**【已达成（§9.579）：F113 = 1 面 / 2 wire（22+6）、`mv=227 mt=227`；a3n00 ratio 0.8996→0.9963、`STATMAP matched=226 unmatched=0`、WIREHIST 1:207 2:10 4:1 6:4 10:4；T0M 未网格 6→2（0.9931，tol 0.025 内）；acs10 0.9984；--lib 1281/0、step_obj_gates 5/5】** 落地两处 .cxx 分歧，且两者必须同时在场（只改 ① 是 Shell(2)，只改 ② 是 Shell(6)）：① `wire_fix.rs::copy_replace_vertices_with` 朝向/位置落定移到 `B.Add` 之前；② `make_faces_on_patch.rs` 三处 `Perform` 传 `RecadreOnPeriodic=false`。剩余旁支（未动手）：WIREHIST 仍差 1 张面（端口 1:207+2:10 vs OCCT 1:208+2:9）；Torus 逐面密度差（§9.373）；§9.439 的 FixDummySeam 修复前提已变可重估。以下为达成前的收敛记录（§9.370 迭代（1））：只剩 3 个面的 seam 缺口（113/140/170），坏在「结果面 vs 结果 shell」】 ① 已定：这 3 个面 `fix_missing_seam` 返回 true 但结果是 **Shell**（5 / 2 / 2 个面），而**合并本身与 OCCT 逐一相同**（后置边数 113: 28=28、140: 16=16、170: 8=8，口径用 `W k edges=N` 不用 `E` 行）⇒ 分歧只能落在 `crates/occt-topo/src/shhealing/shape_fix_face.rs:590-673` 这段尾巴：假想 grid（`cxx:2236-2245`）→ `ComposeShell`（`cxx:2246-2261`）→ `myResult = CompShell.Result()`（`cxx:2263-2268`）→ 两轮剪枝（`FixSmall` / `FixSmallAreaWire`，`cxx:2270-2322`；判据 `crate::shhealing::check_small_area` ↔ `ShapeAnalysis_Wire::CheckSmallArea` `cxx:2004`）。reader 现在只接受 `ShapeType::Face`，所以这 3 个面保持 4/2/2 wires。**修的方向是让那段尾巴给出 OCCT 的单面结果** —— 不是改 seam 位置选择、更不是按下标/bbox 加特例；该段的分歧定位已交给一个并行子任务。② 同批的 138 已核：`--fixms` 对它返回 false，不属 Shell 类（`--all` 里它的 13→15 与 `zz_seam_fix` 的单面口径不同源）。③ 密度（法兰孔面端口约 72 vs GT 52）**先要解决测量**：`zz_probe_a3n00 --fstats` 与 `zz_uv_feed --ids` 的逐面三角数在 **89/226** 个面上不同（而逐面 bbox 226/226 相同）⇒ 逐面读数必须在**同一次运行内**取；下一轮先把逐面统计改成与模型序同源，再逐面比密度。④ 复跑命令见 `specs/_a3n00_gap_analysis.md` §9.370.4。",
      acceptance: "a3n00 面积比从 **0.8996** 上升且不劣化（step_obj_gates 的 step_obj_area 实测口径）；T0M 未网格 **6** 不增加；--lib 1281/0 与 step_obj_gates 5/5 不劣化；不得改基线或 area_tol；不新写测试、不改断言",
      write: "crates/occt-topo/src/shhealing/wire_fix.rs · crates/occt-topo/src/shape_fix_compose_shell/make_faces_on_patch.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.579 · ShapeBuild_Edge.cxx:103-114 · TopoDS_Builder.cxx:74-91 · TopExp.cxx:214-253 · TopoDS_Iterator.cxx:26-83 · ShapeFix_ComposeShell.cxx:3100/3110/3223",
      dependsOn: [],
      evidence: "§9.580（本环境实测）：a3n00 `WIREHIST 1:208 2:9 4:1 6:4 10:4`、`our=227126.78 occ=227130.30 ratio=1.0000`、`TOTAL faces=226 wires=294 faces_with_2plus_wires=18`；A/B（`output/_probe`，OBJ 面积同口径）：只 ④ 时 T0M 0.9931→0.9809（our 191332.77→188978.93）、acs10 0.9984→0.9928、a3n00 0.9978→1.0000；③ 单独不影响 T0M/acs10（a3n00 0.9963→0.9978）；把 cxx:320 写回缓存注释掉，三个模型读数逐位相同。§9.579（本环境实测）：F113 `FDUMP wire[0] nEdges=22 + wire[1] nEdges=6`、`FSTAT face=113 mv=227 mt=227 box=(-87.5,-34,-100)-(87.5,34,-32)`；a3n00 `our=226285.35 occ=227130.30 ratio=0.9963`、`TOTAL faces=226 … STATMAP matched=226 unmatched=0 sum_mt=12051`、`WIREHIST 1:207 2:10 4:1 6:4 10:4`；T0M `STATMAP matched=1770 unmatched=2`、`ratio=0.9931`；acs10 `ratio=0.9984`；`--lib` 1281 passed/0 failed；`step_obj_gates` 5 passed/0 failed；`cargo check` 0 errors。基线/area_tol/断言未动，探针全部撤除。§9.370：`zz_probe_a3n00 --fixms` 实测只有 113/140/170 返回 Shell（5/2/2 个面）；`--all` 里“多改”的 4 个面就是 113/138/140/170，且后置边数与 OCCT 孪生逐一相同（113:28=28、140:16=16、170:8=8，`W k edges=N` 口径）",
      note: "① 的根因**已在 ComposeShell 内部定位**（§9.377）：第一处分歧 = OCCT `ShapeFix_ComposeShell.cxx:2131-2275`(SplitByGrid)→`:1433-1914`(SplitByLine) ↔ 端口 `split_by_grid.rs:93-133` + `split_by_line.rs`/`split_wire.rs`；机制 = breakwires 把闭合 wire 切成 5e/O[A→D] 与 8e/O[C→A] 而 **D≠C**（差 (0.056,2.175,0.056)、D 别无出处）⇒ 两半端点不配对 ⇒ collectwires 对 5e 走 `!index` 分支 ⇒ dispatchwires 出 **同一 patch 两张重合面** ⇒ FixMissingSeam 返回 Shell(2) 被丢弃；剪枝已实测是 no-op（只能删不能合）；§9.378 再收窄一层：新造顶点 D 出自 `split_wire.rs:289-296` 的 `v_opt == None` 分支（对应 cxx:1238-1242，写法本身忠实），所以分歧在**上游**——最可疑是 `cxx:1185-1222` 的顶点匹配分支（端口 `split_wire.rs:258-276` 的 `prev_ok`/`is_coincided`/`split_res` 容差），偏移 2.175 与锥面 v 半幅 2.0207 同量级；§9.379 已核第一条并判为**忠实**（`cxx:1185-1197` 的三条件与端口 `split_wire.rs:251-276` 逐条对应），⇒ 下一个动作只剩「`curr_pnt` 的取法」这一个窗口：对照 `ShapeFix_ComposeShell.cxx:1100-1185` ↔ 端口 `split_wire.rs:150-258`，判据是 `curr_pnt` 是否等于 OCCT 的 `currPnt`，改对后 `zz_seam_fix 140` 应从 Shell(2) 变 Face / 1 wire / 16 边；落地后跑 `pwsh -File .target-gate\t101_verify.ps1`；③ 的逐面统计口径已按 §9.371 修正（FaceMeshStat.index）；§9.372 把密度差按曲面类型归并后 Torus 2.14× 最突出，§9.373 **自我否定**：ParamSet 升序**不是**分歧（OCCT 的 FUN_CalcAverageDUV 形参是非 const 引用、会就地排序 aParamArray，两侧抽稀都在升序上做）⇒ 不要动 ParamSet；③ **已基本结案**：新增 OCCT `--facestats <defl> <angle>`（同参数、模型序、可配对）后重测 —— 法兰/倒角/孔那 16 个面 port 1193 vs OCCT 1188（1.004×，逐面 74/72↔74/72），「72 vs GT 52」是 `--uvsum` 写死 0.5 rad 角度的口径问题（其逐面 triangles 和只有 8054，20° 网格是 12324）；真正还差的只有 Torus 7 个面（1.241×，+293）。下一步只做 Torus 的九量对拍；ParamSet 升序已证伪、不要动",
    },
        { id: "T-102",
      lane: "mesh",
      title: "T0M 残差收口：`load_wires.rs::wire_data_edges` 对 REVERSED wire 的边做二次反转已去掉（`ShapeFix_ComposeShell.cxx:629` 的 `TopoDS_Iterator` 已复合一次）⇒ 1778 面 / 1921 wire / `WIREHIST 1:1672 2:92 3:8 5:1 6:1 7:2 8:2` 与 OCCT 逐项相同；§9.585 的两条残差（4 张环面碎片、`surf=1300` 多拆）同时归位",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 8,
      actual: 5,
      startedAt: "2026-10-03",
      updatedAt: "2026-10-04",
      completedAt: "2026-10-04",
      blocker: "",
      goal: "把 T0M 的面/wire 结构对齐到 OCCT：`--wires` 真值 `TOTAL faces=1778 wires=1921 WIREHIST 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`；T0M.stp 自身只有 1772 个 `ADVANCED_FACE` ⇒ OCCT 在读入时**多造 6 张**。本轮已按 .cxx 把 reader 丢弃多面结果这一条修掉（1772→1775 面、1915→1916 wires、HIST `1:1671 2:90 3:8 5:1 6:1 7:2 8:2`，T0M 面积比 0.9809→1.0084、acs10 0.9928→1.0045，其余 21 个模型逐位不变）。剩余 +3 面 / +5 wire 已用 `FixShape` 子模式 bisect（§9.583）拆成两个独立子项：**2 张 loop-wire 面**（`FixLoopWire` 未移植，端口 census 628（12 边）/ 635（16 边）现 1 wire，OCCT 同 bbox 2 wire ⇒ +2 wire / 0 面）与 **3 张碎片面**（`FixMissingSeam` 那条线，§9.582-E）。",
      next: "**已收口（§9.587，D33）** ⇒ 网格侧残差另立 T-103（T0M 逐面密度 +6.0%、薄面簇 1 张不出网格、4 张多 wire 面的 `MULTI` 盒值差 9e-3）。以下是本卡收口前的记录。两条彼此独立、都对着 .cxx，建议先做 (1)（缺口与判据都更硬）：(1) 2 张 loop-wire 面 —— 按 `ShapeFix_Face.cxx:583`（`NeedFix(myFixLoopWiresMode) && FixLoopWire(aLoopWires)`，`cxx:591-598` 把多条 wire 逐条 `B.Add(tmpFace, …)`）移植 `FixLoopWire`（`cxx:2478`）、`FindNext`（`cxx:2398`）与 `ShapeAnalysis_Wire::CheckLoop`（`ShapeAnalysis_Wire.cxx:2228-2340`，需端口 `ShapeFix_Wire` 的 `Analyzer()` 等价物与 `FixReorder`）；验收看端口 census 628/635 变 `wires=2`、`--wirehist` 从 `1:1671 2:90 …` 变 `1:1669 2:92 …`（loop 拆分后 1 桶 -2、2 桶 +2，多 wire 桶即与 OCCT 逐项相同，剩 3 张单 wire 碎片面）；(2) 3 张碎片面 —— 审 `ShapeFix_Face.cxx:365-480`（`ShapeFix_Wire::Perform`，`ShapeFix_Shape.cxx:200` 把 `ModifyTopologyMode` 置 true）与端口 `shhealing/shape_fix_face.rs::perform_fix_missing_seam` 已标 UNPORTED 的那一段，先回答『为什么 OCCT 在 1459/1499/1643/1662 上会造出碎片面而端口的 seam 步直接返回 false』，另查端口 face 1752（bbox=(-9.5,-9.5,-60.189)-(9.5,9.5,-44.644)，`wires=2 edges=[5,1]`）该并成 OCCT 的 1 面 1 wire。bisect 只能用 `occt_probe --fix <DE_ShapeFixParameters 字段> <int>`。测点已就位：OCCT `specs\\occt_probe\\probe.bat data\\occ\\T0M.stp --wires [--fix …]` / `--faceids`（含 bbox）/ `--facestats 2.574312 0.349066` / `--nofix`，端口 `zz_probe_a3n00 --wirehist` / `--ecensus` / `--fstats` / `--fixms`。§9.585 已把残差精确分解（tie-out：1775+4−1=1778 面、1918+4−1=1921 wire，桶差全在 `1:`）。**§9.586 更正 + 推进**：① §9.585-B 的『`noop`/`raw` 下碎片已存在 ⇒ transfer 期产物』**是错的** —— 当时误写 `fms all zz noop`，`noop` 落到 `argv[5]` 未被探针识别（判定只认 `argv[2] || argv[4]`），那次读数就是默认读数；正确 `fms all noop` = 766 行、`noop` 单体 `faces=1772` ⇒ **4 张碎片连同另 2 张面都是 healing 期**（`ShapeFix_Shape` face 分支）产物；`fms all raw`（只清 `FromSTEP.exec.op`）仍是 772 行，不是 noop oracle；`occt_probe --fix FixWireMode 0` 也回 1772 ⇒ +6 面全出在 1st wire round + `FixMissingSeam` 这一支。② 新 oracle `run_dbg.bat <step> emul <rawFaceIdx> noop`（`ShapeFix_Shape` → `FixFaceTool()->Perform()`，插桩在 `_dbg\\ZZ_ShapeFix_Face.cxx`）：T0M 原始面 1643（`type=4 wires=2 edges_per_wire: 4 1 [UV-DEGENERATE-WIRE]`）→ `ZZFMS compshell faces=2 uf=2.93637 vf=5.49779 ...`，`f0` 9 边 + `f1` **2 边**（碎片 1648）。③ 据此刻出端口 `shhealing/shape_fix_face.rs` 的 `pick`/`vshift` **索引反转**（`ShapeFix_Face.cxx:1996` 的 `coord = ismodeu ? 1 : 0` 配合 `GetFaceUVBounds(F, UMin, UMax, VMin, VMax)`：`coord==1` 选 **V**，端口把 `coord==1` 映到 U，且 `cxx:2077-2078` 的 `min/max` 更新轴、`cxx:2107` 的 `SetCoord(coord+1,·)` 分量同向错）⇒ `shiftw2` 从 0 变 2π、`vf2` 停在 0；修复后 `SEAMGRID uf=2.93637 vf=5.49779` 与 OCCT **逐位一致**，且 T0M/a3n00/acs10/TDB 直方图与 `step_obj_gates` 全部不变。④ 但残差未动：新增 `COMPDBG=1` 插桩（`shape_fix_compose_shell/perform.rs`）显示四张源面 409/411/414/415 全部 `out faces=1`、`f0 wires=1 edges=[9]` ⇒ 端口 `ComposeShell` 只造出 OCCT 的 `f0`、**丢掉 2 边碎片 `f1`**，而 grid 原点已逐位相同 ⇒ 缺口在 `ComposeShell` 内部（`split_by_grid`/`break_wires`/`collect_wires`/`dispatch_wires`/`make_faces_on_patch`）。下一步按工作量：(1)（大）按 `ShapeFix_ComposeShell.cxx:206-270` 的 `Perform` 顺序逐段对拍 OCCT `emul 1643 noop` 的 `ZZSW-IN`/`ZZSEGW`（`_dbg\\ZZ_ComposeShell.cxx`）与端口 `split_wire`/`break_wires`/`dispatch_wires` 的段序列，定位端口少切/少分类的那一段（OCCT 该面 4 次 `SplitWire`，其中一次切在退化段 `(-57.7,-20.379,-38.6911)`、切点 `V=(-57.7,-21.1963,-38.7503)` 正是 `f1` bbox 角点）；(2)（小）端口 `surf=1300`（bbox=(-9.5,-9.5,-60.1893)-(9.5,9.5,-44.6439)）的 seam 步多拆出 Shell(2)，而 OCCT `fms 1758` = `ret=0`、该面是 1 面 9 边 ⇒ 回 `ShapeFix_Face.cxx:1722-2330`（`MakeFacesOnPatch`/`CompShell` 段 cxx:2236-2325 的守卫）找端口缺的 false 分支；去掉后应得 1774 面 / 1917 wire。",
      acceptance: "**已满足（§9.587）**：T0M `--wirehist` = `TOTAL faces=1778 wires=1921` + `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`（= OCCT `--wires` 真值，逐项相同）；`--ecensus` 的 106 条 `MULTI` 两侧 `wires` 全同；a3n00 `WIREHIST 1:208 2:9 4:1 6:4 10:4` 与 ratio 1.0000 未回落；`step_obj_gates` 5/5（面积表仅 T0M 1.0055→1.0053，其余 22 个逐位不变）。原表述：T0M `--wirehist` 的 `1:` 桶从 1669 补到 1672、wires 1918→1921、faces 1775→1778（多 wire 桶不得劣化；`2:` 92 / `3:` 8 / `5:` 1 / `6:` 1 / `7:` 2 / `8:` 2 已逐项相同）；a3n00 `WIREHIST 1:208 2:9 4:1 6:4 10:4` 与 ratio 1.0000 不得回落；T0M/acs10 面积比不得回落（当前 1.0053 / 0.9976）；step_obj_gates 5/5；不得改基线 / area_tol / 断言，不得新写测试",
      write: "crates/occt-topo/src/shape_fix_compose_shell/load_wires.rs（`wire_data_edges`，§9.587 收口点）· crates/occt-topo/src/shhealing/shape_fix_face.rs（`FixLoopWire` + 调用点 + wire 首段）· crates/occt-topo/src/shhealing/wire_fix.rs（`ShapeFix_Wire::Analyzer`/`FixReorder` 支撑）· crates/occt-topo/src/step/read_topology.rs（resolve_face / resolve_shell）",
      ref: "specs/_a3n00_gap_analysis.md §9.587 · §9.586 · §9.585 · §9.583 · ShapeAnalysis.cxx:268-273 :48-62 · ShapeFix_ComposeShell.cxx:206-270 :586-607 :629-640 · TopoDS_Iterator.cxx:26-83 · ShapeAnalysis_Wire.cxx:2228-2340 · ShapeFix_Face.cxx:583-598 :1996-2127 :2077-2078 :2107 :2396-2404 :2398 :2478 :365-480 :1722-2330 :2236-2325 :2266-2268 · ShapeFix_Shape.cxx:196-210 :705-717 :200 :257 :294 · ShapeBuild_ReShape.cxx:282-299 · XSControl_Reader.cxx:480-534 · XSAlgo_ShapeProcessor.cxx:609-620 · D33 · D32 · D31",
      dependsOn: ["T-101"],
      evidence: "**2026-10-04，§9.587 实测（收口：`load_wires.rs::wire_data_edges` 去掉二次反转，Rust 侧 1 处对齐性修复，TEMP 仪器全撤）**：① 根因对着 cxx：`ShapeFix_ComposeShell.cxx:629` 的 `TopoDS_Iterator aIt(wire)`（cumOri 默认 true）已把 wire 朝向与子边朝向复合一次、`cxx:635 sbwdM->Add(E)` 存的就是这个复合值，端口 `edges_of_wire`（`cumulated_children`）语义相同 ⇒ 旧 `wire_data_edges` 在 REVERSED 时再逐边 `reverse()` 属二次反转，`cxx:586-607` 的 2D `WireOrder` 因此拿到反向 pcurve（`Shifted [4,3,2,1]` 而非 `Same [1,2,3,4]`），`FixReorder` 后朝向与 OCCT 相反。② 修复后 T0M `--wirehist` = `TOTAL faces=1778 wires=1921 faces_with_2plus_wires=106 wires_1edge_closed=217 faces_1edge_closed=144` + `WIREHIST 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`（修复前 1775 / 1918 / `1:1669 …`）= OCCT `--wires` 真值，逐项相同。③ §9.585-C(1) 的 4 张环面碎片（OCCT faceid 1462/1503/1648/1668）端口都在：`--ecensus` 0-based 1461/1502/1647/1667、`wires=1 edges=[2]`、`vbox` 最长边 0.129/0.129/0.278/0.278（同区其它面 >1.0）；§9.585-C(2) 的 `surf=1300` 多拆消失：face=1758 = `wires=1 edges=[9] bbox=(-9.5,-9.5,-60.189263)-(9.5,9.5,-44.643872)` = OCCT face 1758（`wires=1 edges_per_wire: 9`）；`--fixms` 对当前读入形状的所有周期面 `ret=false`（seam 步幂等）。④ A/B（23 模型面积表，两态均 `step_obj_gates` 5/5）：旧 body T0M `our=193713.95 occ=192658.64 ratio=1.0055` → 新 body `our=193679.05 ratio=1.0053`，**唯一变化**；其余 22 个逐位不变（a3n00 227126.78 / 1.0000、acs10 270565.78 / 0.9976、TDB 237464.11 / 1.0110、motoc 1.0333、bottom 0.9750、top 1.0000 …）。⑤ `--wirehist`：a3n00 `226/294 1:208 2:9 4:1 6:4 10:4`、acs10 `787/913 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`，与 §9.586-D 逐项相同。⑥ 验收第 2 条：`--ecensus` 的 `MULTI` 对 `output/t0m_wires_occt.txt` 两侧各 106 条、`wires` 全同，4 条盒值差约 9e-3；`--fstats`（`T0M.stp 2.574312 --fstats`）`TOTAL faces=1778 computed_lin=2.574309 used_lin=2.574312 stats=1777 mesh_v=62303 mesh_t=70557`、`STATMAP matched=1777 unmatched=1`（§9.581 旧态 `stats=1765 unmatched=7`），无 stats 的是 0-based face=1764（薄面簇 `x∈[-23.73,-21.27] y∈[-206.86,-204.33] z∈[-107.16,-104.24]`）。⑦ 门禁：`cargo check --all-targets` 0 errors、`--test step_obj_gates` 5/5、`--test 'phase*'` 61 passed / 10 suites；`LOADDBG`/`COMPDBG`/`BWDBG`/`SEAMDBG` 插桩全部删除。⑧ 剩余旁支（另立 T-103）：T0M 逐面三角数和 70557 vs OCCT `FSTAT_OCC` 求和 66576（+6.0%）、face 1764 无 stats、4 张多 wire 面的 `MULTI` 盒值差 9e-3。\n**2026-10-04，§9.586 实测（§9.585-B 更正 + `pick` 修复 + 缺口钉进 ComposeShell，Rust 侧 1 处对齐性修复）**：① `run_wires.bat data\\occ\\T0M.stp fms all noop` = **766 行**（默认 772）、`... noop` 单体 = `TOTAL faces=1772 wires=2121 hist 1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1`、`fms all raw` 仍 772 行 ⇒ §9.585-B 的『zz noop / raw 逐字节相同 ⇒ transfer 期产物』结论作废（`noop` 判定只认 `argv[2] || argv[4]`，`fms all zz noop` 把 `noop` 放在 `argv[5]`），4 张环面碎片是 **healing 期**产物。② `occt_probe --fix <field> 0`：`FixWireMode=0` → 1772；`FixMissingSeamMode=0` → 1973；`FixSplitFaceMode/FixAddNaturalBoundMode/RemoveSmallAreaFaceMode/FixSeamMode/FixLoopWiresMode = 0` → 1778。③ 新 oracle `run_dbg.bat <step> emul <rawFaceIdx> noop`（`ShapeFix_Shape` → `FixFaceTool()->Perform()`，插桩 `_dbg\\ZZ_ShapeFix_Face.cxx`）：T0M 原始面 1643 → `ZZFMS in f0 wires=2 edges=5 bbox=(-58.1,-25.5664,-47.7556)-(-57.7,-15.4336,-37.6229)`、`ZZFMS-grid coord=1 isneg=1 period=6.28319 nb1=4 nb2=1 m1=(2.93705,9.22024,1.5708,3.14159) m2=(-3.14159,3.14159,3.14159,3.14159) uf=2.93637 vf=-0.785398`、`ZZFMS compshell faces=2 wires=2 uf=2.93637 vf=5.49779 ...`（`f0` 9 边 + `f1` 2 边 bbox=(-58.1,-21.2794,-38.7525)-(-57.7,-21.1963,-38.3564) = 碎片 1648），与 `ZZSW-IN` 4 次 `SplitWire`（含退化段 `(-57.7,-20.379,-38.6911)` 上切点 `V=(-57.7,-21.1963,-38.7503)`）。④ 端口 `pick`/`vshift` 索引反转（详见 `note`）修复后 `SEAMGRID uf=2.93637 vf=5.49779 URange=6.28319 VRange=6.28319 uclosed=1 vclosed=1 ismodeu=1 ismodev=0 nb1=4 nb2=1 m1=(2.93328,9.21956,1.57080,3.14159) m2=(-3.14159,3.14159,3.14159,3.14159)` 与 OCCT 逐位一致。⑤ 修复后读数不变：T0M `--wirehist` 仍 `1:1669 2:92 3:8 5:1 6:1 7:2 8:2`（1775/1918）、a3n00 `1:208 2:9 4:1 6:4 10:4`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` 全部未变，`step_obj_gates` 5/5、`cargo check --release` 0 error、`--lib` 1255/26（26 条仍是 `%TEMP%` 写权限环境项）。⑥ 新端口插桩 `COMPDBG=1`（`shape_fix_compose_shell/perform.rs`）：源面 409/411/414/415 全部 `COMPDBG out faces=1`、`f0 wires=1 edges=[9]` ⇒ 端口 `ComposeShell` 丢掉 OCCT 的 2 边碎片 `f1`；grid 原点既已逐位相同，缺口在 `ComposeShell` 内部而非 seam 参数。\n**2026-10-03，§9.583 实测（探针 bisect，Rust 侧零改动）**：① `FixFaceMode=0` ⇒ T0M `faces=1772 wires=2121 HIST 1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1`，与 `--nofix` 逐位相同 ⇒ 比文件本体多的 6 张面全出自 `ShapeFix_Shape.cxx:196-210` 的 face 分支（D31 机制坐实）。② `FixLoopWiresMode=0` ⇒ `faces=1778 wires=1919 HIST 1:1674 2:90 3:8 5:1 6:1 7:2 8:2`（面数不变、`2:92→2:90`）⇒ `ShapeFix_Face.cxx:583-598 FixLoopWire` 把 2 张单 wire 面各拆成 2 wire；`--faceids` 开关前后逐 bbox 对差只有 2 个 bbox：端口 census 628（`wires=1 edges=[12]` bbox=(-29.5,-105.566610,13.433081)-(12.5,-76.275610,19.126675)）/ 635（`wires=1 edges=[16]` bbox=(-29.5,-134.709636,9.346976)-(12.5,-106.496089,19.061680)），OCCT `--wires` 同 bbox 都是 `wires=2`。③ 其余单点关（AddNaturalBound / SplitFace / SmallAreaWire / RemoveSmallAreaFace / Orientation / IntersectingWires / Reorder / Small / Connected / AutoCorrectPrecision / Seam / Shifted / NotchedEdges / FixWireMode）与默认逐位相同；`FixMissingSeamMode=0` 反成 1973 面（+`FixSplitFaceMode=0` 回 1772）⇒ +201 是 `FixAddNaturalBound` 的 `NeedSplit` 交互（`cxx:705-717`），不是 seam 步造面。【§9.586 更正】③ 里把 `FixWireMode` 列为「与默认逐位相同」不成立，实测 `FixWireMode=0` → 1772。④ 仪器纠错：`SetShapeProcessFlags`/`SetShapeFixParameters` 与 `FromSTEP.FixShape.*` 静态量在 `ReadFile` 前设置全部 no-op（actor 未建 + 8.0 走 `DE_ShapeFixParameters` map）⇒ 新增 `--fix <字段> <int>`（读后、`TransferRoots` 前）才能做 bisect。⑤ 端口本轮 `--ecensus` = `faces=1775 wires=1916 HIST 1:1671 2:90 3:8 5:1 6:1 7:2 8:2`（与 §9.582 一致）。\n**2026-10-03，§9.582 实测**：① T0M.stp `ADVANCED_FACE` = 1772（文件本体），OCCT `--nofix`（ShapeProcess 全关）faces=1772 `HIST 1:1463 2:294 3:9 5:1 6:1 7:2 8:1 10:1`，OCCT 默认 faces=1778 `1:1672 2:92 3:8 5:1 6:1 7:2 8:2` ⇒ 端口读进来就是**文件本身**的面数，OCCT 多 6 张全部出自 `FromSTEP.exec.op = FixShape`（`STEPControl_Controller.cxx:201`，`SplitClosedFaces` 等不在读路径）。② OCCT `--faceids`（本轮 TEMP 加 bbox）11 张 id=0（healing 造、无 STEP 实体）面按 bbox 与端口 census 配对：116↔116 d=0（端口 2 wires / OCCT 1）、117/118 无对应（d=2.6/3.7）、1461/1502 与端口 1459/1499 d=0、1462/1503 无对应（d=1.35）、1647/1667 与端口 1643/1662 d=0.28/0.30、1648/1668 无对应（d=2.17）⇒ 缺口 = 1 张主面拆成 3 张（+2）+ 4 张碎片（+4）。③ 端口 `--fixms`：只有 face=116 → `Shell faces=3`、face=1752 → `Shell faces=2` 返回 true，其余全 false。④ 落地（`resolve_face` 接受 Shell 结果并对其中每张面跑 `check_pcurves_and_shift`；`resolve_shell` 把 Shell 结果的面逐张加入本 shell）= `ShapeFix_Face.cxx:2268 Context()->Replace(myFace, myResult)` + `ShapeFix_Shape.cxx:257/294 Context()->Apply(S)` + `ShapeBuild_ReShape.cxx:282-299`。⑤ 对拍（`--release --test step_obj_gates`，23 模型）：T0M 188978.93→194286.39（ratio 0.9809→1.0084）、acs10 269261.65→272445.08（0.9928→1.0045）、a3n00 1.0000 / TDB 1.0092 / motoc 1.0333 / bottom 0.9750 / ATU01038 1.0006 / top 1.0000 等 21 个逐位不变，5/5 通过；a3n00 `--fixms` 无一条 `ret=true` ⇒ 该路径不触发。⑥ 端口 `--wirehist` T0M = `1775 / 1916 / 1:1671 2:90 3:8 5:1 6:1 7:2 8:2`。\n**2026-10-04，§9.585 实测（终态免责 + 逐面唯一分解）**：① 终态两侧 `FixMissingSeam` **都是 no-op** —— OCCT `fms all` = 772 张周期面全 `ret=0`（`Result()`=null）、端口 `--fixms` = 769 张全 `ret=false`，`before_wires` 直方图 766×1/5×2/1×8 vs 763×1/5×2/1×8，5 张 2-wire 周期面**同序号**（190/197/204/211/359）。② 【§9.586 更正】原写『`fms all` 与 `fms all zz noop`、`fms all raw` 三者逐字节相同 ⇒ healing 零贡献』——`fms all zz noop` 其实没触发 noop（`noop` 在 `argv[5]`），正确 noop 是 766 行；`raw` 仍 772 行。结论反转为：**healing 造 +6 面，4 张碎片在其中**。③ 端口 transfer 期 seam 步：1775 面里 212 面 `ret=true`，结果 1770×单面 + **2×Shell**（`nf=2`、`nf=3`）⇒ `MakeFacesOnPatch`/`ComposeShell` 只走到 2 次。④ 残差 tie-out 精确：缺 4 张 2 边 Torus 碎片（OCCT 1462/1503/1648/1668，各 1 wire、u span 0.0031、`uper=1 vper=1`）+ 多 1 张面/1 wire（端口 `surf=1300` 的 seam 步返回 `Shell(2)`：sub0 5 边 z∈[-56.134,-44.644] + sub1 4 边 z∈[-60.189,-56.134]，OCCT 同面 = `FACE 1758 type=1 wires=1 edges_per_wire: 9 spans u=6.2832 v=15.5454` 且 `fms 1758` = `ret=0`）⇒ 1775+4−1=1778 面、1918+4−1=1921 wire，桶差全在 `1:`（1669+3=1672）。⑤ 端口终态 `--wirehist` = `1:1669 2:92 3:8 5:1 6:1 7:2 8:2`（1775/1918，§9.583 的 +2 wire 已在内），a3n00 = `1:208 2:9 4:1 6:4 10:4`（226/294，未回退）；`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`cargo check --release` 0 error。",
      note: "【§9.587 收口】D32 的两条入口（『`ComposeShell` 少切/少分类』与『`ShapeFix_Face.cxx:2236-2325` 缺 false 分支』）都作废：真正的根因是 `load_wires.rs::wire_data_edges` 对 REVERSED wire 的边做了二次反转（`ShapeFix_ComposeShell.cxx:629` 的 `TopoDS_Iterator` 已复合一次），修掉后 4 张环面碎片与 `surf=1300` 的多拆同时归位，T0M 面/wire 结构与 OCCT 逐项相同（见 D33）。`FixLoopWiresMode=0` 那条 bisect 结论本身仍成立，但它不是本轮残差的入口。上一段记录：§9.580-F 的假设（『idx116/idx351 各应是 1 条 wire，缺口在 `CollectWires` 的 wire 合并』）已被 §9.582 推翻：OCCT face 116 是 healing **新建**的面（id=0），与端口同 bbox 的那张 2-wire 面是被替换掉的**原始 STEP 面**，两者不是同一张；`CollectWires` 不是本轮入口。§9.583 又把剩余缺口拆成『2 张 loop-wire 面（`FixLoopWire` 未移植）+ 3 张碎片面』两个独立子项，并纠出探针仪器的一个坑：`--set` / 读前的 `SetShapeFixParameters` 全是 no-op，bisect 必须用 `--fix`（`wires_probe.cpp` 的 `setp` 才是早前 §9.17x 用的有效通道，`set`/`setc` 那几条结论要重估）。红线：④（`MapReShape::apply_impl`）仍不许撤（§9.581）。旁支（未动）：`--lib` 在本机 26 个写 `%TEMP%` 的文件 IO 用例因 os error 5 失败，与本次改动无关（写 workspace 的 5 个 gate 用例全过）；`--faceids` 的 bbox、`--fix`/`--set` 都是 TEMP 仪器，用完可撤。§9.585 补两条结论：(a)（**【§9.586 更正】** 原文写「残差全在 transfer 期、seam 步不是缺口来源」——该结论建立在 `fms all zz noop` 上，而 `noop` 落在 `argv[5]` 未被探针识别，等于又跑了一遍默认；正确 `fms all noop` = 766 行、`noop` 单体 `faces=1772` ⇒ **4 张碎片与另 2 张面都是 healing 期产物**，缺口确在 `ShapeFix_Face` 那条线；`raw` 只清 `FromSTEP.exec.op` 不关 healing，也不是 oracle）(b) 残差 3 面/3 wire 的唯一分解 = 4 张环面碎片（OCCT 有、端口无，`u` span 0.0031 的 2 边 Torus）+ 1 张端口多拆（`surf=1300` 在 OCCT 上 `fms` 返回 false、不拆），tie-out 精确。【§9.586 新增】(c) 据 OCCT `emul 1643 noop` + `_dbg\\ZZ_ShapeFix_Face.cxx` 的 `ZZFMS compshell faces=2 ...`，端口 `shhealing/shape_fix_face.rs` 的 `pick`/`vshift` **索引反转**已修（`coord = ismodeu ? 1 : 0` 配合 `GetFaceUVBounds(F, UMin, UMax, VMin, VMax)` ⇒ `coord==1` 选 V；端口三处同向错：`pick` 的轴、`cxx:2077-2078` 的 `min/max` 更新轴、`cxx:2107` 的 `SetCoord(coord+1,·)` 分量）——修复后 `SEAMGRID uf/vf` 与 OCCT 逐位一致，属**对齐性修复**，T0M/a3n00/acs10/TDB 直方图与门禁全不变；(d) 残差本身未动，`COMPDBG=1` 插桩（`shape_fix_compose_shell/perform.rs`）坐实端口 `ComposeShell` 出 `faces=1 / edges=[9]`、丢掉 OCCT 的 2 边碎片 `f1` ⇒ 下一入口在 `ShapeFix_ComposeShell.cxx:206-270` 的 `Perform` 分段。本轮 TEMP 仪器：`read_topology.rs::resolve_face` 的 `SEAMDBG=1` + `shhealing/shape_fix_face.rs` 的 `SEAMGRID` + `shape_fix_compose_shell/perform.rs` 的 `COMPDBG=1`（均 env 开关，默认零开销，可撤）+ `wires_probe.cpp` 的 `raw` 改为全参数扫描（保留，可与 `fms all` 组合）；`specs/occt_probe/_dbg/ZZ_{ShapeFix_Face,ComposeShell}.cxx` 的靶窗已加 T0M 四张环面面（`build_dbg.bat` → `zz_wires_probe.exe`，用 `run_dbg.bat` 跑）；`output/` 下对拍产物未入库。工具坑：`wires_probe.exe` 直接跑会因缺 OCCT DLL 静默 0 字节退出，必须走 `specs\\occt_probe\\run_wires.bat`（或 `call env.bat vc14 64`）；`zz_wires_probe.exe` 走 `run_dbg.bat`。",
    },
        { id: "T-103",
      lane: "mesh",
      title: "T0M 网格侧残差：ON-only 塌陷面 8/8 归位（§9.608：pcurve 求值改用 COS range，`pf=28 1→74`、`pf=812 20→70`）；残差收窄到共享态拆缝面的参数域密度（Torus `pf=27` 等 ±1 节点面）",
      status: "active",
      priority: "P1",
      owner: "agent",
      progress: 0,
      estimate: 6,
      actual: 0,
      startedAt: "2026-10-04",
      updatedAt: "2026-10-06",
      completedAt: "",
      blocker: "",
      goal: "把 T0M 的**网格**读数对齐到 OCCT。拓扑侧已收口（T-102/§9.587：1778 面 / 1921 wire / `WIREHIST 1:1672 2:92 3:8 5:1 6:1 7:2 8:2` 与 `--wires` 真值逐项相同），所以 T0M 面积比 1.0053 的残差只能落在逐面网格密度与不出网格的那张面上：【§9.591 更新 —— 空面已修：`stats=1778`、`unmatched=0`、逐面 `mt` 求和 `65338`；以下 `70557` / `unmatched=1` 均为修复前读数】同偏转（`2.574312`）下端口 `TOTAL faces=1778 stats=1777 mesh_v=62303 mesh_t=70557`，OCCT `--facestats 2.574312 0.349066` 的 1778 行逐面 `nodes/triangles` 求和为 `60050 / 66576` ⇒ 三角数 +6.0%；`STATMAP matched=1777 unmatched=1`，缺 stats 的是 0-based face=1764（薄面簇 `x∈[-23.73,-21.27] y∈[-206.86,-204.33] z∈[-107.16,-104.24]`）。另含拓扑侧子项：共享 `ShapeBuild_ReShape`（`OCCT_SHARED_HEAL`）是 `pf=1693` 拿到第 5 条边的唯一来源，默认关闭，按 §9.595 推进。",
      next: "**【已停止 — 2026-10-06】本线在 §9.608 闭环后停止：ON-only 塌陷第一现场曾是 §9.603 分类器多边形自重叠（`collect_boundary_uv` 的链自重叠、环绕数退化），§9.608 把 2D pcurve 求值改用 `BRep_Tool::CurveOnSurface` 的 COS 表示自身 range 后 **8/8 登记面全部归位**（ON 侧 `pf=28 1→74`、`pf=812 20→70`，均 = OCCT；默认路径 diff=0）；剩余残差是共享态拆缝面的参数域层，不再推进。本轮清理临时仪器后门禁实测：五 crate `--all-targets` 0 error；`step_obj_gates` 5 passed；`fuse_box_cylinder_is_closed_solid` 1 passed；默认路径四模型 `--wirehist` 逐项不变（a3n00 `226/294 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）；`export_data_obj` 15 ok / 0 err。** **§9.608 更新（本轮，网格塌陷入口闭环）**：`pf=812`/`pf=28` 的 ON-only 塌陷根因 = 2D pcurve 求值用了**边的 3D range**而非 `BRep_Tool::CurveOnSurface` 返回的 **COS 表示自身 range**（`BRep_Tool.cxx:353` 的 `GC->Range(First,Last)`；`ShapeAnalysis_Edge::PCurve`（`ShapeAnalysis_Edge.cxx:192-208`，`cxx:200`）把它交给 `ShapeAnalysis_Wire::CheckOrder`（`ShapeAnalysis_Wire.cxx:593-651`，`cxx:634/648-649`））。拆缝子边的 pcurve 范围被 `ShapeBuild_Edge::CopyRanges`（`ShapeBuild_Edge.cxx:206-334`，`cxx:289-290`）重标到切片窗口，仍按父边 3D range 求值就整体偏一个周期 ⇒ §9.603 的 ON 分类器多边形自重叠。② 修复：`wire_builder.rs` 新增 `pcurve_and_range(e, face)`（= `boptools_2d::curve_on_surface_oriented(e, face, false)`，缺 COS 范围回退 3D range），`add_wire`/`edge_uv_ends`/`wire_uv_points`/`wire_classifier`/seam 反向判别 5 处统一改用；无新谓词。③ 单二进制 A/B（T0M，同偏转 2.574312；默认路径 diff **0**/1778）：ON `pf=28 1→74`、`pf=812 20→70`、`pf=1663 38→73`，`pf=27 254→301`；`ON-only` 面数 **56→54**；ON 合计 `mt 65331→65536`（OCCT 66640）。④ 门禁：`cargo check --all-targets` 0 err、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1，四模型默认 `--wirehist` 逐项不变。⑤ 旁支（未动）：`pf=27`（Torus）仍 ON-only `192/301` vs `185/288`（与其余 ±1 节点的 ON-only 面同类）；并行 T-104 未提交改动使 OFF 读数在 §9.600 基线上漂移（本卡改动对 OFF 逐位零影响）。**§9.602 更新（本轮）**：余 4 张 ON-only 塌陷面（`pf=212/198/812/28`）的第一现场已钉死，2 张归位。① 判据（`--ecensus`/`--fsurf` OFF↔ON 对读）：`pf=212`/`pf=198` 两面 OFF/ON 都是 `wires=1 edges=[6]`、bbox 与顶点签名逐位相同 ⇒ 不是「`LoadWires` 已 `Apply` 但边序不连通」也不是「没走 `LoadWires`」；根因是预修复 wire 回合用了 `wires_of_face(&perf_face)`（等价 `TopoDS_Iterator iter(S, true)`，把 `Face.ori × Wire.ori` 复合进来），而 OCCT 是 `TopoDS_Iterator iter(S, false)`（`ShapeFix_Face.cxx:400` / `:532`）⇒ 两张面 face 与 wire 双 `Reversed` 时复合结果反而 `Forward`，`ShapeExtend_WireData::Init`（`cxx:114-121`）从尾插 ⇒ `ShapeAnalysis_WireOrder` 拿到反向 pcurve ⇒ `FixReorder` 后边序与 OCCT 相反。② 修复：`read_topology.rs` 新增 `wires_stored_on_face`（直接读 `Face` 的 TShape WIRE children、`cumOri=false`），`resolve_face` 共享上下文分支与 seam 后第二轮（`:923`）都改用它；无上下文（默认路径）走 `resolve_loop` 的 `w` 不动。③ `pf=812`（Cone）/`pf=28`（Cylinder）：ON `edges=[5]`/`[6]` = OCCT 真值（`output/t0m_occt_wires.txt:1099` `FACE 812 ... edges_per_wire: 5`、`:67` `FACE 28 ... 6`；OFF 4/5），wire 已连通、边序正确 ⇒ 塌陷在**网格侧**（拆缝面参数域离散 + 三角化），不属本轮。④ 数字：ON `pf=212 12→178`、`pf=198 12→174`（= OCCT）；`pf=812` 仍 20（OCCT 70）、`pf=28` 仍 1（OCCT 74）⇒ 8 张登记面现 6 张归位。⑤ 门禁：`cargo check --all-targets` 0 err、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1；默认路径四模型 `--wirehist` 逐项不变，共享态四模型也逐项相同。⑥ 撤除 `fdbg_dump_face_wires`（含 2 份重复定义）+ `T103RV`/`T103R`/`T103D`；未动 `wire_order.rs`/`wire_builder.rs`。**下一步**：`pf=812`/`pf=28` 的网格塌陷入口（拆缝面 `urange` 已整周期，查 `--fsurf`/`--fstats` 的 uv 窗口与三角化），归 T-103 网格线。**§9.603 更新（网格第一现场）**：`CBDBG`/`CWDBG` 逐边 dump 定位 —— `pf=28`（Cylinder）OFF 5 边分类器多边形是干净闭合矩形 `(1.5708,0)->(7.854,0)->(7.854,23.75)->(1.5708,23.75)`（mt=74 = OCCT 72）；ON 6 边（= OCCT `edges_per_wire: 6`）的链下边被绕一圈后又跳回 `u=7.854` 起点（`(4.53786,0)->(7.85398,0)`）⇒ 多边形自重叠、`mt=74→1`。`pf=812`（Cone）同理（OFF `nb=[35,4,32,5]` 正常 mt=70 = OCCT 70；ON `nb=[4,32,5,9,27]`、mt=20）。切法与 OCCT `collectWirePoints`（`BRepMesh_NodeInsertionMeshAlgo.hxx:142-183`）逐句一致 ⇒ 缺口在 ON 拆缝边（`edge61/62` 等）的 pcurve 点序/朝向（或同一条 seam `edge63` 在 e2/e5 两次取用的 u 平移）。**再下一步**：先核 `pcurve.points()` 是否复用缓存缓冲（`CBWIN` 的 `first/last` 与 `CWDBG` 链段区间对不上），再对 OCCT 侧取 `pf=28`/`pf=812` 逐边 UV 区间（`ZZWIN`）逐段比对。**§9.600 更新（本轮）**：§9.599-F 登记的 ON-only 塌陷面（`pf=110 312→8` 等 8 张）根因定位并修复两层。① 对 `.cxx` 审『收尾段』：`SetLast` 在 OCCT 里是**注释掉的**（`ShapeFix_ComposeShell.cxx:2344`），`ShapeAnalysis_WireOrder`（`cxx:589-606`）与 `sbwdM->Reverse(face)`（`cxx:629`）端口自 §9.587 已有 ⇒ 缺的不是起点归一；真正缺的是 `cxx:505-506` 的 `LoadWires` **入口** `Context()->Apply(iw.Value())`，以及 `ShapeBuild_ReShape.cxx:284` 的 `TopoDS_Iterator aSubIt(aNewShape)`（`cumOri` 默认 true，`TopoDS_Iterator.cxx:73-82` 复合朝向）端口未做。② 修复：`load_wires.rs` 加 `context` 参数并对每个 child 先 `ctx.apply`；`reshape.rs::apply_impl` 子形状朝向改为 `Orientation::compose(nc_ori, c.orientation())`；`perform.rs` 传 `&self.context`。③ 门禁：`cargo check --all-targets` 0 err、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1；默认路径四模型 `--wirehist` 逐项不变，共享态也逐项相同。④ 网格（`--fstats` T0M，同偏转 2.574312）：OFF 64988 → ON **65453**（上一版 65051，本轮 **+402**；OCCT 66576）；8 张登记面 5 张归位（`pf=110 8→353`=OCCT、`pf=816 5→70`、`pf=1668 5→96`、`pf=1462 5→97`、`pf=812 8→20` 半改善），剩 `pf=212 178→12`、`pf=198 174→12`、`pf=28 74→1` 未动。⑤ 已撤 `T103W`/`T103S`/`RSDBG`/`SLWDBG`/`CWDGB` 等 TEMP 仪器（未动 `wire_order.rs`/`wire_builder.rs`）。**下一步**：取 `pf=212`（Torus）做第一现场，判它是『`LoadWires` 已 `Apply` 但边序仍不连通』还是『走了另一条不经过 `LoadWires` 的重建路径』。以下是 §9.599 记录。【§9.599 更新（本轮，闭环）】**：共享态 `2-wire` 桶 +70 的**产生点已定位并修复**。① 产生点（非 `break_wires`/`split_by_grid`/`ApplyContext` 记账）：端口在**逐面翻译期**就跑 `ShapeFix_Wire::Perform` 的 `check_pcurves_and_shift`，传的是**原始 `face`**；而 OCCT 跑在 `S = Context()->Apply(myFace)`（`ShapeFix_Face.cxx:379-380`）上。共享态下邻面 Cylinder 的 `DispatchWires` 已 `Replace(old_edge,new_edge)`（`ShapeFix_ComposeShell.cxx:3451-3457`，`new_edge` 只带 Cylinder 那份 pcurve），于是本面 pcurve 被 `FixAddPCurve` 挂到 `old_edge`，随后 `Context()->Apply`（`ShapeFix_Shape.cxx:257`）用 `new_edge` 替换 ⇒ 重建后的 wire 取不到本面 pcurve ⇒ `check_wire` 返回 `None` ⇒ `w2=None` ⇒ `reason=no-w2-branch` ⇒ seam 两条 wire 不合并，该面以 2 wire 计数。② 反证：把该修复切回旧行为，`OCCT_SHARED_HEAL=1` 的 T0M 立刻回到 `1777/1990 1:1601 2:162 …`、`pf=1693` 回 4 边 ⇒ 与本修复一一对应。③ 修复：`read_topology.rs::resolve_face` 的 wire 回合，在共享上下文存在时先 `ctx.apply(&face.0)` 再跑 `check_pcurves_and_shift`；无上下文走原路径（**默认路径零改动**）。④ 验收：`cargo check --all-targets` 0 err、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1；默认路径逐项不变（a3n00 226/294 `1:208 2:9 4:1 6:4 10:4`、T0M 1778/1921 `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 787/913 `1:742 …`、TDB 2180/2385 `1:2042 …`）；共享态现在**逐项相同**（四模型 faces/wires/WIREHIST 全同）；`pf=1693` `--fdump` OFF 4 边 → ON **5 边**，内赤道圆 `R−r=9.5` 在 `(2.5, 9.165151, -18.0)` 切成 `e[2]+e[3]`（= OCCT 的 e2+e3）。⑤ 旁支：ON 下 `wires_1edge_closed` 计数与 OFF 不同（a3n00 58 vs 79、acs10 133 vs 157）而普查一致，未定性；`shape_fix_compose_shell/`、`shhealing/` 的 `RSDBG`/`FMSDBG` 仪器未撤（与并行 T-104 共文件）；共享态仍未默认打开。⑥ 网格层下游（`--fstats` 同偏转 `2.574312`）：**主残差面对齐 OCCT** —— `pf=1693`（↔ OCCT `of=390`，`nodes=445 triangles=795`）OFF `372/650` → ON **`445/795`**（逐位相同）；四张环面碎片 `pf=1461/1502/1647/1667` `5/3 → 93/97`。全模型逐面 `mt` 差 ON−OFF：`up=87 down=20 same=1671 delta=+63`（`sum 64988→65051`）。**新旁支**：20 张面在 ON 反而塌掉（`pf=110 312→8`、`pf=212 178→12`、`pf=198 174→12`、`pf=1668 74→5`、`pf=816 70→5`、`pf=812 70→8`、`pf=28 74→18`、`pf=1462 56→5`）；抽检 `pf=110`（Torus）OFF 为 10 边闭合环、ON 多切 1 边（11 边，新顶点 `(-15.229917,-230.745772,29.168011)`）但**边串联顺序错乱**（`e[0]: -21→-15.23`、`e[1]: -15.23→-7.456`，`e[2]` 又跳回 `-21` 起头）⇒ wire 不连通、2D 分类器退化。与 §9.599-A 不同源，属共享态拆分后**缺 wire 重排/起点归一**（`ShapeFix_ComposeShell.cxx` 的 `sbwdM`/`WireOrder` 收尾段）。**§9.597 更新（本轮）**：T-103 环面 U 网格残差根因锁定 + §9.596-E 证伪。① 残差第一现场（`pf=1693` R=12.5/r=3 ↔ OCCT `of=390`）：两侧 `TORDBG/ZZTORDBG` 的 `R/r/ru/rv/diff_u/diff_v/nb_u/nb_v/du/dv/old_dv` **逐位相同**，只有 `nup 72 vs 73`（⇒ `nu 37 vs 46`、`nodes 280 vs 352`）⇒ 分歧只在 `fillParams` 的**输入参数表**，不在 splitter。② OCCT `ZZU face=390` 的 U 表 = 10° 格（37 值）∪ 一组 `0.171678724051821` 间隔值；`0.171678724051821 = 4.97868299750295 / 29`、`29 = ceil(4.97868299750295 / 0.1745329)`，而 `4.97868299750295` 正是 OCCT `ZZWN e=2 diff`（`e=3` 的 `1.30450230967664/8 = 0.16306` 亦对上）⇒ 这组非 10° 格的值 = **内赤道圆（`R−r=9.5`）被切成 `e2`+`e3` 两条边后各自的 pcurve 点**。端口 `pf=1693` 的 wire 仍 4 边（`e[2]` = 完整 9.5 圆，`PCURDBG key=… nb=37 first=3.1416 last=9.4248`），没有这组点 ⇒ 缺的是**共享合缝把该圆切开**（`ShapeFix_Shape.cxx:257` 的 `Context()->Apply(S)` + `ShapeFix_Face.cxx:2266/2268` 那次共享 `Replace(圆→两弧 wire)`，§9.594-A），不是 `TorusRangeSplitter`/`fillParams`。③ §9.596-E 的「邻面 pcurve 槽互删」假设**经 .cxx 证伪**：`BRep_CurveOnSurface::IsCurveOnSurface(S,L)`（`BRep_CurveOnSurface.cxx:58-62`）= `(S==mySurface)&&(L==myLocation)` 即 **`Geom_Surface` handle 相等**，`EmptyCopy`（`BRep_TFace.cxx:37-44`）共享同一 surface handle ⇒ OCCT 侧同样会删两份；端口 `tgeometry.rs:485-497` 按 `Arc::as_ptr` 删单槽同构 ⇒ 无分歧，该入口关闭。④ 共享态（`OCCT_SHARED_HEAL=1`）**已不挂死**（T0M `--fdump/--wirehist` ~15 s 完成，§9.594-C 两条 hang 消失），但仍过度拆分：T0M `faces 1777 / wires 1990 / 1:1601 2:162 3:8 5:1 6:1 7:2 8:2`（OCCT `1778/1921/1:1672 2:92 …`，2-wire +70）、a3n00 `226/299 1:203 2:14 4:1 6:4 10:4`（OFF `226/294 1:208 2:9 …`）⇒ 默认路径逐项不变、该 pass 仍不能默认打开。⑤ 验收：`cargo check --all-targets` 0 err、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1；本轮库代码零改动（TEMP 仪器 + 本节）。**下一步**：定位共享态 2-wire +70 的产生点（`break_wires`/`split_by_grid` 对 `ApplyContext` 的记账 + 重建 wire 为什么没被 `fix_missing_seam` 合回去，a3n00 成功数 46→41），目标 `OCCT_SHARED_HEAL=1` 下 `--wirehist` 回到 `1:1672 2:92 …` 且 `pf=1693` 出现第 5 条边。**§9.596 更新（本轮）**：共享态过度拆分的入口已钉死且修掉一层 —— a3n00 face 130 同一条几何 `OFF u=[3.141593,9.424778] v=[3,8] wire[0] nEdges=4` → `ON u=[0,2π] v=[-inf,inf] wire[0] nEdges=1, wire[1] nEdges=1`；`u=[0,2π] v=[-inf,inf]` 是 `AddUVBounds` 的**自然界回退**，说明该面每条边都取不到 pcurve。根因：pcurve 槽以**面的 `Geom_Surface` 指针**为键（`tgeometry.rs::repr_key`），而 `applyImpl`（`ShapeBuild_ReShape.cxx:249-254`）重建 FACE 时 §9.595-B 的 `copy_face_geom` 只写 `face_core().surface`、**没登记侧表** ⇒ 新面 `repr_key` 回退成面 TShape 指针 ⇒ 全丢 pcurve（实测重建后「所有边无 pcurve」的面 24 处，全出自 `Resolver::resolve_representation → fix_periodic_degenerated → SharedReShape::apply → replace`）。改法：`tgeometry.rs` 新增 `GeometryRegistry::register_face_surface`，`reshape.rs::copy_face_geom` 复制 surface 后登记；对齐 `BRep_Tool::CurveOnSurface`（`BRep_Tool.cxx:347-357`）直接用面自身 surface 匹配表示的语义。验收：`cargo check --all-targets` 0 err、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` ok；默认路径逐项不变（a3n00 208/9/1/4/4、T0M 1672/92/8/1/1/2/2、acs10 742/25/1/18/1、TDB 2042/112/13/3/2/1/4/3）；共享态「所有边无 pcurve」面 24→0、`uv_bounds` 不再退化，但 `fix_missing_seam` 成功数仍 41（OFF 46）。**下一步**：剩 **29 条边**性质不同 —— `Context()->Apply(myFace)` 为恒等（`same_ts=true`、`frec/erec=false`），但边 pcurve 表示键指向的 surface 指针不等于面自身（`ne_reps≠face_surf`），即这些边挂的是**另一个 patch surface** 的 pcurve、本面槽被删；已排除 `CopyPCurves/CopyReplaceVertices`（按 `repr_key` 全量复制）与按指针删单槽的 `RemovePCurve/ReassignPCurve`。下一入口：查**同一 surface 被多张面共用**时 `ShapeBuild_Edge::RemovePCurve(edge, myFace)`（`dispatch_wires.rs:114` ↔ `ShapeFix_ComposeShell.cxx:3348-3357`）是否连带删掉邻面的槽（OCCT 的 `BRep_GCurve` 同样按 (surface,location) 共享，需给「OCCT 为何不丢」的可对拍证据，不能只按端口这一处改）；再复核 `pf=1693` 是否变 5 边。**§9.595 更新**：共享 `ShapeFix_Shape` 上下文（`OCCT_SHARED_HEAL=1`）两处已修 —— (1) `MapReShape` 的键形状改成与替换体同存（`HashMap<usize,(TopoShape,TopoShape)>`，`shape_fix_compose_shell/reshape.rs`），对齐 `BRepTools_ReShape.hxx:239-240` 的 `NCollection_DataMap<TopoDS_Shape,…>` 生命周期，修掉「指针复用命中错替换体」导致的 a3n00 225/226 不稳定与 T0M 卡死/panic；(2) `apply_impl` 重建 FACE 时补 `copy_face_geom`，对齐 `BRep_TFace::EmptyCopy`（`BRep_TFace.cxx:37-44` 复制 `Surface/Location/Tolerance`，不复制 `NaturalRestriction`），共享态无 surface 的面 383→30。默认路径基线逐项不变（a3n00 226/294、T0M 1778/1921、acs10 787/913、TDB 2180/2385），`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` ok。**下一步**：共享态仍过度拆分（T0M 2-wire 面 92→162、wires 1921→1990、faces 1778→1777；a3n00 9→14、wires 294→299），最小复现是 a3n00 一条 18 边闭合 wire 被拆成 10 条（`wires=10 edges=[1,1,2,2,2,2,2,2,2,2]`），入口在 `ShapeFix_Face` 第二轮 loop 拆分（`FixLoopWire` `cxx:2478` / `FindNext` `cxx:2398` / `CheckLoop`）在重建后的 wire 上判据失配，不是 `applyImpl`；另有 30(T0M)/8(a3n00) 张无 surface 面出自另一条创建路径（`apply_impl` 的 67 次 FACE 重建全部带 surface，已排除）。两处补齐后才复核 `pf=1693` 是否变 5 边。**§9.593 更新（本轮闭环）**：`node_insertion.rs::collect_boundary_uv` 的分类器多边形对 **seam 边**（`reverse_walk=true`）漏 XOR，两个 seam 被反向走成尖刺 ⇒ 净环绕变顺时针 ⇒ `BRepMesh_Classifier` 的 `TabOrient=false` ⇒ `Perform` 把矩形内所有候选点判 `TopAbs_OUT`（`pf=1693 surf=280 → kept=0`）。改法与 frontier 链接同用 `reverse` 判定（`!is_forward() ^ reverse_walk`，与 `BRepMesh_NodeInsertionMeshAlgo.hxx:143-181` + `IMeshData_PCurve::IsForward`（`hxx:52`）等价：cxx 的 `GetEdgeOrientation` 本身已含 `Edge.Reverse()`）。T0M `--fstats` 逐面求和 `mt 65338→66874`（OCCT 67346，`-2.98%→-0.70%`）、`_t103_p2.mjs` `bit-identical 1622→1634/1723`、Torus ratio `0.921→0.991`、BSpline `0.960→0.998`；`pf=1693` `92/90→372/650`。剩余 `pf=1693 372/445`、`pf=8 224/154` 属下面 (2) 的 Torus 参数域密度（`of=390/310` 是 §9.592-D 的 5 张 U/V 网格不一致面），非分类器。**§9.591 更新**：空面一项已闭环 —— `FixPeriodicDegenerated` 移植后 face 1764 有 stats（`mv=38 mt=35`），FSTAT 1777→**1778**，三角数基线 65303→**65338**（对 OCCT 66576 = -1.86%），故下面 (1) 的「(3) face 1764 单独看」作废、(4) 的 `MULTI` 盒值差与 (2) 的 Torus 参数域密度仍是本卡主线。(1) 【本轮已做，§9.589】配对键已解决：OCCT `--facestats` 与 `--wires` 的**面序不同**（1778 行按 box 全可映但 1772 行换了索引，同索引只 6 行相同），所以跨模式只能按几何配；按 bbox（容差 0.05）端口↔OCCT 配对得 `paired=1722 unmatched=55`、`1462/1722 (85%)` 逐位相同、残差全部集中在 **Torus**（n=74，`+4122 mt`，ratio 1.264），其余类型都在 ±1.6% 内。(2) 【下一步】顺该结论查喂给 splitter 的点：OCCT `collectWirePoints`（`BRepMesh_NodeInsertionMeshAlgo.hxx:142-183`）对每条边**只喂 n−1 个点**（`IsForward` 丢最后一个、反向丢第一个），端口正常路径（`node_insertion.rs:494-519` 喂 `wires_uv`）已一致，但**回退路径**（`node_insertion.rs:266-272`）喂的是 `self.boundary_uv` = 全量 `uv`（`node_insertion.rs:483`/`660-673`，含重复顶点），会把 `fillParams` 的 `aLength` 抬高、`aStdStep = aDiff/aLength` 变小 ⇒ 环面网格变密。先加 env 开关可判该批 Torus 面走的是 1 还是 2，命中就按 `hxx:142-183` 让回退路径也喂 n−1 序列（不动 `uv`/`p3t` 前沿结构）。`generate_surface_nodes`（`splitter.rs:534-648`）本身已对着 `D:\\source\\OCCT-src\\...\\BRepMesh_TorusRangeSplitter.cxx:22-115` 逐句审过、判定忠实，不要再动那一层。(3) face 1764 单独看：`mv=0` 是几何退化还是端口 `IMeshData_Failure` 被 rescue 三角化（见 D-x「rescue 空分支」记录）。(4) 4 张多 wire 面的 `MULTI` 盒值差 9e-3 一并量清（`wires` 数相同，疑为盒值口径）。测点已就位：OCCT `specs\\occt_probe\\run_here.bat <step> --facestats 2.574312 0.349066` / `--wires`，端口 `zz_probe_a3n00 --fstats` / `--fsurf` / `--ecensus`（`TORSDBG=1` 打环面 splitter 中间量）。",
      acceptance: "T0M `--fstats` 的 `matched/unmatched` 与逐面 `mv/mt` 对 OCCT 收敛（`unmatched` 1→0；三角数总和向 `66576` 收敛），且 `--wirehist` 的 `1:1672 2:92 …` 与 `--ecensus` 的 106 条 `MULTI` 不许回落；`step_obj_gates` 5/5（T0M 从 1.0053 向 1.0000 收敛、其余 22 个逐位不变）；a3n00 ratio 1.0000 不许回落；不得改基线 / area_tol / 断言，不得新写测试，不得为 T0M 单独调偏转或加 OCCT 里没有的谓词。",
      write: "crates/occt-topo/src/meshing/**（UV 覆盖 / 参数域离散 / 三角化，按命中层）· crates/occt-topo/src/shhealing/shape_fix_face.rs（若落在三角化前的面修补）",
      ref: "specs/_a3n00_gap_analysis.md §9.603 · §9.608 · §9.602 · §9.600 · §9.599 · §9.597 · §9.596 · §9.595 · §9.587-F · §9.581 · §9.373 · ShapeAnalysis_Edge.cxx:192-208 · BRep_Tool.cxx:327-361 · ShapeAnalysis_Wire.cxx:593-651 · ShapeBuild_Edge.cxx:206-334 · ShapeFix_Face.cxx:400 :532 :379-380 :2266 :2268 · ShapeFix_ComposeShell.cxx:505-506 :589-606 :629-630 :2344 :2767 :3451-3457 :3348-3357 · ShapeExtend_WireData.cxx:114-121 · ShapeBuild_ReShape.cxx:284 :317 · TopoDS_Iterator.cxx:73-82 · ShapeFix_Shape.cxx:257 · BRep_TFace.cxx:37-44 · BRep_CurveOnSurface.cxx:58-62 · BRepTools_ReShape.hxx:239-240 · ShapeBuild_ReShape.cxx:249-254 · BRepMesh_DefaultRangeSplitter.cxx:35-41 :45-81 :202-236 · BRepMesh_NodeInsertionMeshAlgo.hxx:142-184 · BaseMeshAlgo.cxx:52-62 · D22 · D30",
      dependsOn: ["T-102"],
      evidence: "**【已停止 — 2026-10-06】T-103 在 §9.608 后停止（`pf=812`/`pf=28` 已由 §9.608 归位，8/8 登记面全归位；剩余为共享态拆缝面参数域层，未动）；本次调查遗留的临时仪器已全部撤除（完整清单见 §9.608-F：`TORSDBG`/`TORDUMPU`/`TORDUMPV`/`PCURDBG`/`MPFDBG`/`CBDBG`/`CWDBG`/`CLSDBG`/`PCWDBG`/`CIRCDBG`/`STEPEDBG`/`STEPCDBG`/`PCUV`/`NSDBG`/`DELDBG`/`NVDBG`/`FWDBG`/`F2DBG`/`F4DBG`/`FLDBG`/`F7DBG`/`HBDBG`/`NOHEALSNAP`/`MPF2`；`MBDIAG`/`FDBG` 按约定保留）。本轮清理后门禁实测：五 crate `--all-targets` 0 error；`step_obj_gates` 5 passed；`fuse_box_cylinder_is_closed_solid` 1 passed；默认 `--wirehist` 四模型逐项不变；`export_data_obj` 15 ok / 0 err。** **2026-10-06，§9.608 实测（T-103 网格塌陷入口闭环：`pf=28`/`pf=812` 归位，8 张登记面全归位）**：① 根因对着 cxx：`BRepMesh_ShapeVisitor::addWire` 的 2D `CheckOrder` 分支读 `ShapeAnalysis_Edge::PCurve`（`ShapeAnalysis_Edge.cxx:192-208`，`cxx:200` `C2d = BRep_Tool::CurveOnSurface(edge, surface, location, cf, cl)`）→ `BRep_Tool::CurveOnSurface`（`BRep_Tool.cxx:327-361`）在 `cxx:353` 取 `GC->Range(First, Last)` = **COS 表示自身**范围，不是边 3D range；`ShapeAnalysis_Wire::CheckOrder`（`ShapeAnalysis_Wire.cxx:593-651`）在 `cxx:648-649` 用 `c2d->Value(f)/(l)`。拆缝子边的 range 由 `ShapeBuild_Edge::CopyRanges`（`ShapeBuild_Edge.cxx:206-334`，`cxx:289-290` `newF=first+alpha*len`）重标到切片窗口 ⇒ 端口旧代码在 5 处都用 `BRepTool::edge_parameters`（3D range）求值 ⇒ 端点整体偏一个周期 ⇒ ON 分类器多边形自重叠（§9.603 第一现场）。② 修复：`wire_builder.rs` 新增 `pcurve_and_range`（`boptools_2d::curve_on_surface_oriented(e, face, false)` + 3D range 回退），`add_wire`（CheckOrder）/`edge_uv_ends`/`wire_uv_points`/`wire_classifier`（TotCross2D + 无限点）/seam 反向判别 5 处统一改用；未加谓词，未动 `wire_order.rs`。③ 数字（T0M，同偏转 2.574312，**单二进制** A/B，`T103_PCR3D` 开关已撤）：默认路径（`OCCT_SHARED_HEAL` 未设）old→new **diff=0/1778**（逐位相同，无回退）；ON old→new 只动 4 张：`pf=28 76/1 → 76/74`（OCCT 74/72）、`pf=812 72/20 → 72/70`（OCCT 72/70 = 逐位相同）、`pf=1663 75/38 → 75/73`（OCCT 75/73 = 逐位相同）、`pf=27 185/254 → 192/301`（OCCT 185/288）；`ON-only`（按 bbox 配对，OFF↔ON 不等）面数 **56 → 54**；ON 逐面合计 `mv 59622→59629`、`mt 65331→65536`（OCCT `FSTAT_OCC` `60083/66640`），OFF 合计 `59433/65217`。④ 门禁：`cargo check --manifest-path crates/occt-topo/Cargo.toml --all-targets` 0 error；`cargo test --test step_obj_gates` 5 passed / 0 failed；`cargo test fuse_box_cylinder_is_closed_solid` 1 passed；默认路径四模型 `--wirehist` 逐项不变（a3n00 `226/294 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）。⑤ 旁支（未动）：`pf=27`（Torus，`pf=28` 的拆缝兄弟）仍 ON-only `192/301` vs `185/288`，与其余 ±1 节点的 ON-only 面（`pf=1/12/14/38/79` 等 `72→73`）同类，属共享态拆缝后参数域层；并行 T-104 未提交改动使默认路径 OFF 读数在 §9.600 基线（`59455/64988`）上漂移（本卡改动对 OFF 逐位零影响，以同树同二进制 A/B 为准）；本轮 env 开关 `T103_PCR3D` 已撤除，`wire_builder.rs` 无残留散件。**2026-10-05，§9.602 实测（余 4 张 ON-only 面：2 张归位，2 张定性为网格侧；TEMP 仪器撤除）**：① `--fdump`/`--ecensus`/`--fsurf` OFF↔ON 对读（同偏转）——`face=212 OFF wires=1 edges=[6] mt=178 -> ON wires=1 edges=[6] mt=178`、`face=198 174 -> 174`、`face=812 OFF [4] mt=70 -> ON [5] mt=20`、`face=28 OFF [5] mt=74 -> ON [6] mt=1`。② `pf=212/198` 两面 OFF/ON 的 bbox、边数、顶点签名逐位相同 ⇒ 不是 `LoadWires` 未 Apply / 非 `LoadWires` 路径，是**朝向/边序**：旧 `wires_of_face(&perf_face)` 复合 `Face.ori × Wire.ori`，两张面 face+wire 双 `Reversed` ⇒ 复合为 `Forward`，`ShapeExtend_WireData::Init`（`ShapeExtend_WireData.cxx:114-121`）尾插 ⇒ `ShapeAnalysis_WireOrder` 反向 ⇒ `FixReorder` 边序与 OCCT 相反。③ 修复：`read_topology.rs::wires_stored_on_face`（读 `Face` TShape 的 WIRE children，`cumOri=false`，对齐 `ShapeFix_Face.cxx:400` / `:532`），`resolve_face` 共享分支 + seam 后第二轮（`:923`）改用它；默认路径不动。④ 数字：ON `pf=212 12→178`、`pf=198 12→174`（= OCCT）；`pf=812` 20、`pf=28` 1 未动（ON `edges=[5]`/`[6]` = OCCT 真值 `output/t0m_occt_wires.txt:1099` `FACE 812 edges_per_wire: 5` / `:67` `FACE 28 ... 6` ⇒ 拓扑已对齐、塌陷在网格侧：拆缝面参数域 + 三角化）。⑤ 门禁：`cargo check --all-targets` 0 err、`step_obj_gates` 5 passed (599.84s)、`fuse_box_cylinder_is_closed_solid` 1 passed；默认路径 `--wirehist` 四模型逐项不变（a3n00 226/294 `1:208 2:9 4:1 6:4 10:4`、T0M 1778/1921 `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 787/913 `1:742 2:25 5:1 6:18 8:1`、TDB 2180/2385 `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`），共享态四模型同样逐项相同。⑥ 撤除 TEMP 仪器：`read_topology.rs::fdbg_dump_face_wires`（2 份重复定义 + 2 处调用）、`wire_fix.rs` 的 `T103RV`/`T103R`(×2)/`T103D`；未动 `meshing/wire_order.rs`/`model_builder/wire_builder.rs`。⑦ 旁支：ON 的 `wires_1edge_closed/faces_1edge_closed` 与 OFF 不同（a3n00 58/24 vs 79/24、T0M 172/129 vs 216/143、acs10 133/72 vs 157/72、TDB 155/91 相同），faces/wires 直方图两边一致，未定性。**2026-10-05，§9.600 实测（ON-only 塌陷面两层修复；默认路径逐项不变）**：① 对 `.cxx` 审收尾段：`ShapeFix_ComposeShell.cxx:2344` 的 `// wire.SetLast(j-1)` 是注释掉的（未启用）；`ShapeAnalysis_WireOrder sawo` 按 `cxx:589-606` 使用、端口自 §9.587 已有；`sbwdM->Reverse(face)`（`cxx:629`）端口有 ⇒ 缺的不是起点归一。② 两处真实缺口：`cxx:505-506` `LoadWires` 入口的 `Context()->Apply(iw.Value())`（端口旧 `load_wires.rs` 直接读 children）；`ShapeBuild_ReShape.cxx:284` 的 `TopoDS_Iterator aSubIt(aNewShape)` `cumOri` 默认 true（`TopoDS_Iterator.cxx:73-82` `Compose(myOrientation, sub)`），端口 `reshape.rs::apply_impl` 未复合替换体朝向。③ 修复：`load_wires.rs` 加 `context: &SharedReShape` + 每 child `ctx.apply`、`reshape.rs::apply_impl` 用 `Orientation::compose(nc_ori, c.orientation())`、`perform.rs` 传 `&self.context`（无新谓词/启发式）。④ 门禁：`cargo check --manifest-path crates\\occt-topo\\Cargo.toml --all-targets` 0 error、`step_obj_gates` 5 passed、`--lib fuse_box_cylinder_is_closed_solid` 1 passed；默认路径 `--wirehist` 四模型逐项不变（a3n00 226/294 `1:208 2:9 4:1 6:4 10:4`；acs10 787/913 `1:742 2:25 5:1 6:18 8:1`；T0M 1778/1921 `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`；TDB 2180/2385 `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`），共享态同样逐项相同。⑤ 网格效果（T0M `--fstats`，同偏转 2.574312）：逐面 `mt` 求和 OFF `64988` → ON **`65453`**（上一版 ON `65051` = 本轮 **+402**；OCCT `FSTAT_OCC` 求和 `66576`）；8 张登记面：`pf=110 8→353`（= OCCT 353）、`pf=816 5→70`、`pf=1668 5→96`、`pf=1462 5→97` 归位，`pf=812 8→20` 半改善，`pf=212 178→12` / `pf=198 174→12` / `pf=28 74→1` 未动 ⇒ 仍 ON-only 塌陷的只剩 4 张。⑥ 撤除 TEMP 仪器：`dispatch_wires.rs::T103W`、`split_wire.rs::T103S`、`wire_fix.rs`/`shape_build_edge.rs` 的 `RSDBG`、`loop_wire.rs::SLWDBG`、`shape_fix_face.rs::CWDGB` + 残留 `ei`、`reshape.rs` 残留未用 import；**未动** `meshing/wire_order.rs` / `model_builder/wire_builder.rs`。**2026-10-05，§9.599 实测（共享态 `2-wire` +70 闭环；默认路径逐项不变）**：① 只读定位钉到产生点：端口 `read_topology.rs::resolve_face` 的 wire 回合（`check_pcurves_and_shift`）在**翻译期**跑在原始 `face` 上，OCCT 跑在 `S = Context()->Apply(myFace)`（`ShapeFix_Face.cxx:379-380`）；共享态邻面的 `Replace(old_edge,new_edge)`（`ShapeFix_ComposeShell.cxx:3451-3457`）+ 最终 `Context()->Apply(S)`（`ShapeFix_Shape.cxx:257`）使本面 pcurve 落空 ⇒ `fix_missing_seam` 的 `check_wire` 返回 `None` ⇒ `w2=None` ⇒ `reason=no-w2-branch` ⇒ 该面停留 2 wire。② 反证（同工作树，仅切回旧行为）：ON T0M `TOTAL faces=1777 wires=1990`、`WIREHIST 1:1601 2:162 3:8 5:1 6:1 7:2 8:2`、`pf=1693` 4 边。③ 修复后：ON T0M `1778/1921 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、ON a3n00 `226/294 1:208 2:9 4:1 6:4 10:4`、ON acs10 `787/913 1:742 2:25 5:1 6:18 8:1`、ON TDB `2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` —— 与各自 OFF 逐项相同；`--fdump 1693` ON `wire[0] nEdges=5`（`e[2] e[3]` 在 `(2.5,9.165151,-18.0)` 切开）。④ 门禁：`cargo check --all-targets` 0 error、`step_obj_gates` 5 passed / 0 failed、`fuse_box_cylinder_is_closed_solid` 1 passed。⑤ 撤除的临时仪器：`read_topology.rs` 的 `ASSOCDBG`/`PKEYDBG`/`APDBG`、`tgeometry.rs` 的 `RSPCLOST`/`PCSET`/`PCSETCOPY`/`RSPCREM`；`shape_fix_compose_shell/`、`shhealing/` 的 `RSDBG`/`FMSDBG` 留在原地（T-104 共文件，未动）。**§9.595 实测（共享态两处对齐性修复，默认路径逐项不变）**：① 缺口一（键生命周期）：OCCT `BRepTools_ReShape::myMap`（`BRepTools_ReShape.hxx:239-240`）是 `NCollection_DataMap<TopoDS_Shape, TReplacement, TopTools_ShapeMapHasher>`，键即持 TShape 的 `TopoDS_Shape`；端口旧 `MapReShape` 只存替换体、键用 `Arc::as_ptr`，临时子形状 drop 后地址被 Face 复用 ⇒ 命中错替换体（0 边替换体 + 死循环）。改为 `HashMap<usize,(TopoShape,TopoShape)>`。② 缺口二（FACE EmptyCopy）：OCCT `applyImpl` 用 `EmptyCopied()`（`ShapeBuild_ReShape.cxx:249-254`），`BRep_TFace::EmptyCopy`（`BRep_TFace.cxx:37-44`）复制 `mySurface`/`myLocation`/`myTolerance`；端口 `apply_impl` 只对 EDGE 复制几何 ⇒ 共享态 `Context()->Apply(S)`（`ShapeFix_Shape.cxx:257`）把带被替换子形状的面重建成无 `Geom_Surface` 的壳面。新增 `copy_face_geom`。③ 验收：`cargo check --all-targets` 0 error；`step_obj_gates` 5/5；`fuse_box_cylinder_is_closed_solid` 1/1；默认路径 `a3n00 226/294 1:208 2:9 4:1 6:4 10:4`、`T0M 1778/1921 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、`acs10 787/913 1:742 2:25 5:1 6:18 8:1`、`TDB 2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` 逐项不变；`pf=1693` 仍 4 边。④ 共享态：`T0M --ecensus` 无 surface 面 `383 → 30`；`a3n00 --roots` 连跑两次 `faces=226`（§9.594-C 的 225/226 不稳定与 T0M 卡死已消除）。⑤ 临时计数（`RESHAPE_DBG`，已撤）`apply_impl` 的 FACE 重建 67 次全部带 surface ⇒ 残余 30/8 张无 surface 面出自另一路径。⑥ 剩余过度拆分：`T0M 2-wire 92→162 / wires 1921→1990 / faces 1778→1777`、`a3n00 2-wire 9→14 / wires 294→299`，a3n00 一条 18 边闭合 wire 被拆成 `wires=10 edges=[1,1,2,2,2,2,2,2,2,2]`。**2026-10-04，§9.593 实测（分类器 seam 反向，Rust 侧 1 处对齐性修复）**：① 现象：T0M `pf=1693`（bbox 与 OCCT `of=390` 逐位相同）`NVDBG surf=280` 但 `mv=92 mt=90`（OCCT `nodes=445 triangles=795`），全 T0M 共 18 张面同样丢光内点。② 根因对着 header：`collect_boundary_uv` 的分类器多边形只用 `pcurve.is_forward()`，未叠加端口的 `reverse_walk`（= `ShapeExtend_WireData::Edge(signed)` 的 `Edge.Reverse()`）；`BRepMesh_ShapeVisitor::addWire` 存的是**已 reverse 的** edge orientation，`collectWirePoints`（`hxx:143-181`）的 `GetPCurve(face, GetEdgeOrientation(·))` 因而 IsForward 已含该 reverse ⇒ 端口必须 XOR。`CWDBG face=1693` 实测多边形在 U=0/U=2π 两个 seam 均走反、`anAngle<0`、`TabOrient=false`、`in=0 on=0 out=280`；frontier（`uv`）同 seam 已正确反向，两分支不一致即缺口。③ 修复（`crates/occt-topo/src/meshing/node_insertion.rs` 分类器分支改用与 frontier 同一个 `reverse`）。④ 验收：`cargo check` 0 error；`step_obj_gates` 5/5；`fuse_box_cylinder_is_closed_solid` ok；T0M `--fstats` 求和 `mv 59630→60398`（OCCT 60547）、`mt 65338→66874`（OCCT 67346，`-2.98%→-0.70%`）；`_t103_p2.mjs` `paired mt delta -3.17%→-0.82%`、`bit-identical 1622→1634/1723`、Torus `0.921→0.991`（identical 68→69）、BSpline `0.960→0.998`（identical 258→269）；`pf=1693 92/90→372/650`。⑤ 只影响 `reverse_walk=true` 的 seam 边，全 T0M `--fsurf` 仅 18 张面变化且全部上升（`up=18 down=0`）。⑥ wirehist：`a3n00`/`TDB` 逐字相同；`acs10` `158/73→157/72`，临时回退本次修改复跑仍是 `157/72` ⇒ 系本轮之前其它未提交改动，与本次无关。**2026-10-04，§9.591 实测（空面闭环）**：`FixPeriodicDegenerated` + `IsPeriodicConicalLoop` 移植（T-105）后，`zz_probe_a3n00 data\\occ\\T0M.stp --fstats` FSTAT 行 1777→**1778**（失败面归零）、`face=1764 mv=38 mt=35 bbox=(-23.729500,-207.157336,-107.158859)-(-21.270500,-206.058541,-104.833829)` 与 OCCT `FSTAT_OCC face=566` bbox 逐位相同、逐面 `mt` 求和 = **65338**（OCCT 求和 66576 ⇒ -1.86%）；`--ecensus face=1764` 由 `wires=1 edges=[1] nv=1` 变 `wires=1 edges=[5] nv=2`；`step_obj_gates` 5/5（T0M 面积比 1.0053，其余 22 个逐位不变）。下面 §9.590 段里的 `STATMAP unmatched=1`、`mesh_t=70557`、`65303` 均为修复前读数。**2026-10-04，§9.590 实测（空面第一现场 + TDB 无同类 + 基线过期）**：用户报的 T0M 空面 = 唯一无 stats 的 0-based face 1764，插桩 `MBDIAG=1` 实测走 `invalid discrete range`（`node_insertion.rs:537-543`）⇒ `IMeshData_Failure` ⇒ `discret_root.rs:439-453` 静默跳过 ⇒ 空面；下一层读数 `surf=Cone sup=true svp=false su=[0,6.283185] sv=[-inf,inf] du=[0.202683,6.283185] dv=[0.717188,0.717188] buv=0.202683..6.283185 x 0.717188..0.717188 nuv=32` ⇒ 32 个边界 UV 点全落在 cone 的 v=+0.717188 一条线上，V 向 extent=0 ⇒ 网格器的判定是忠实的，缺口在端口该面的 wire（`--ecensus face=1764` = `wires=1 edges=[1] nv=1`，唯一顶点 (-21.687158,-206.759155,-105.124128)）。OCCT 同索引面 `FACE 1764 type=2(Cone) wires=1 edges_per_wire: 4 u=6.2832 v=1.4344` + seam、box 在 -Y 侧多 0.298、`FSTAT_OCC` nodes=72/**triangles=104** ⇒ 要往 healing/`--ecensus` 拿回「4 边 + 全 V 段」的 wire。TDB 同口径 `MBDIAG=1` 零命中（2180/2180 有 stats、无 `mt=0`、22 张 `edges=[1]` 顶点面全出网格、逐面 mesh 面积无退化）⇒ 用户报的 TDB 那一处需指定面/文件。旁支：卡里 `mesh_t=70557` 已过期，当前 `--fstats` 逐面 mt 求和与 `data/output/T0M.obj` 三角数同为 **65303**（OCCT 求和 66576 ⇒ -1.9% 反向）⇒ §9.589 的 Torus 结论须在 65303 基线上重跑。**2026-10-04，§9.589 第一步实测（本轮，逐面残差锁定 Torus）**：① 仪器坑：同一份 OCCT 默认读入的 T0M，`--facestats` 与 `--wires` 的面序列**不是同一张面**（`--facestats` idx 0..11 ↔ `--wires` idx 52/54/55/57/58/59/70/71/73/74/76/77；1778 行按 box 全可映、1772 行换索引；`--wires.box` vs `--facestats.bbox` @same idx = 6/1778）⇒ 跨模式只能用几何配对，按索引配对的旧读法要重估。② 端口↔OCCT 按 bbox（容差 0.05）配对：`paired=1722 unmatched=55`，`paired mt: port=68114 occ=64496 (+5.61%)`，`paired mv: 60264 vs 58242`，`bit-identical (mv 且 mt) = 1462/1722 (85%)`。③ 按端口曲面类型归并：**Torus n=74 `19722 vs 15600` `+4122` ratio 1.264 仅 27 张逐位相同**；BSpline 0.961、Plane 0.986、Cone 1.016、Cylinder 1.002、Sphere 1.003。⇒ T0M 的 +6.0% 基本等于 Torus 的 +26.4%，逐面差最大的 20 张里 19 张是 Torus。③b 环面**节点与三角数同时高**（配对子集 mv `11951/9852` = 1.213、mt `19722/15600` = 1.264，两侧 `mt/mv` 都约 2.0；其余类型 mv 比全在 0.978–1.011）⇒ 落在**参数域密度**层而不是三角化层，优先查喂进 `fillParams` 的参数集；下一步按既有 `specs/occt_probe/_dbg/ZZ_*.cxx` 机制加 `ZZ_BRepMesh_TorusRangeSplitter.cxx` 打在 `GenerateSurfaceNodes` 内（`aDiffU/aDiffV/r/R/oldDv/nbV/Dv/Du/nbU/ParametersU.Length()/ParametersV.Length()/aParamU->Length()/aParamV->Length()/aNodes->Length()`），与 `TORSDBG=1` 逐字段对齐定位分岔点。④ 端口 `TorusRangeSplitter::generate_surface_nodes`（`splitter.rs:534-648`）已对着 `D:\\source\\OCCT-src\\src\\ModelingAlgorithms\\TKMesh\\BRepMesh\\BRepMesh_TorusRangeSplitter.cxx:22-115` 逐句审（`oldDv`/`Dv` 覆盖、`nbV`、`ru>1e-16` 的 `Du*=min(oldDv,Du)/aa`、`nbU=max(·, int(nbV*diffU*R/(diffV*r)/5))`、`R<r` 均匀 vs `fillParams(·,0.5)`、`fillParams(·,2/3)`、`newRange*=±0.1*step` 半开区间发射、`fillParams`/`FUN_CalcAverageDUV`）——**判定忠实**，缺口不在这一层；下一步查喂点（`hxx:142-183` 的 n−1 与端口回退路径喂全量 `boundary_uv` 的差异，见 §9.589-D）。工具新增：端口 `--fsurf`（逐面 `type/mv/mt/wires/urange/vrange/uper/vper/bbox`，输出 `.target-gate/_t0m_fsurf.txt`），配对脚本 `.target-gate/_t103_p2.mjs`（产物 `_t103_pairbox.txt`）、面序验证 `.target-gate/_t103_align.mjs` / `_t103_perm.mjs`。**立项依据（§9.587-D 实测）**：① 端口 `zz_probe_a3n00 data/occ/T0M.stp 2.574312 --fstats` -> `TOTAL faces=1778 computed_lin=2.574309 used_lin=2.574312 stats=1777 mesh_v=62303 mesh_t=70557`、`STATMAP matched=1777 unmatched=1 sum_mt=70557 flat_mt=70557`，缺 stats 的面由 `--fstats` 的 face 索引差集定为 0-based 1764（邻面 1763/1765/1766 的 bbox 都在 `x∈[-23.73,-21.27] y∈[-206.86,-204.33] z∈[-107.16,-104.24]`）。② OCCT `output/t0m_fs_occt.txt`（`--facestats 2.574312 0.349066`）1778 行逐面求和 = `nodes 60050 / triangles 66576` ⇒ 端口三角数 +6.0%。③ `--ecensus` 的 106 条 `MULTI` 对 `output/t0m_wires_occt.txt`：`wires` 逐条相同，4 条盒值差约 9e-3（例：端口 `(-14.000,-140.551,-231.369)-(20.000,-114.076,-222.253)` vs OCCT `(-14.009,-140.560,-231.378)-(20.009,-114.067,-222.244)`）。④ 拓扑侧的读数已全部对齐（T-102/§9.587），因此本卡只对网格层。",
    },
        { id: "T-104",
      lane: "healing",
      title: "`BRepLib::BuildCurve3d` / `SameRange` / `ShapeFix_Edge::FixAddCurve3d` 族移植：旧端口用 `CurveOnSurface` 冒充 3D 曲线，改为真实 B 样条/解析重建（新增 `geom_lib.rs` + `brep_lib_same_range.rs`，纯对齐性补齐）",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 6,
      actual: 5,
      startedAt: "2026-10-04",
      updatedAt: "2026-10-04",
      completedAt: "2026-10-04",
      blocker: "",
      goal: "把 `DispatchWires` 里「给新建边补 3D 曲线」这条链从`CurveOnSurface` 假曲线换成 OCCT 的真实 `BuildCurve3d`：平面臂 `GeomLib::To3d` 解析转 3D，一般臂 `GeomLib::BuildCurve3d` 的 `AdvApprox` 出 B 样条，容差与 `SameParameter`/`SameRange` 标志按 `BRepLib.cxx:301-455` 更新。不改变已对齐的面 / wire 结构（对 4 个基准模型是纯对齐性补齐）。",
      next: "**已收口（§9.588）**：`cargo check --all-targets` 0 errors、`step_obj_gates` 5/5；T0M/a3n00/acs10/TDB 的 `--wirehist` 与 §9.587 逐项相同（1778/1921 `1:1672 2:92 …`、226/294 `1:208 …`、787/913 `1:742 …`、2180/2385 `1:2042 …`）；`--ecensus` 的 `wires>=2` 面数 = 106 = OCCT `MULTI`，§9.583 的两张 loop-wire 面现为 `face=628 wires=2 edges=[11,1]` / `face=635 wires=2 edges=[15,1]`。剩余若还有缺口，都归 T-103（网格层）与 §9.588-E 的 TEMP 仪器清理。",
      acceptance: "**已满足**：编译 0 error、`step_obj_gates` 5/5、4 个基准模型 `--wirehist` 不回落且与 OCCT 逐项相同；未改基线 / area_tol / 断言、未新写测试、未加 OCCT 里没有的谓词。",
      write: "crates/occt-topo/src/shhealing/shape_build_edge.rs · crates/occt-topo/src/shhealing/shape_fix_edge.rs · crates/occt-topo/src/brep_lib_same_range.rs · crates/occt-geom/src/geom_lib.rs · crates/occt-core/src/elib/clib.rs · crates/occt-topo/src/shape_fix_compose_shell/dispatch_wires.rs · crates/occt-geom2d/src/{curve,offset,offset2d,bezier_curve,bspline_curve}.rs · crates/occt-topo/src/tgeometry.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.588 · BRepLib.cxx:149-183 :187-263 :275-295 :301-455 · ShapeBuild_Edge.cxx:714-775 · ShapeFix_Edge.cxx:335-464 :618-638 · GeomLib.cxx:559-679 :1051-1165 :2991-3077 :3079-3221 · ElCLib.cxx:1339-1422 · ShapeFix_ComposeShell.cxx:3506-3529 · D34",
      dependsOn: ["T-102"],
      evidence: "**2026-10-04，§9.588 实测（本机）**：① 源侧：新增 `brep_lib_same_range.rs`（`CheckSameRange` / `SameRange` / `rep_ranges` / `set_rep_ranges`）与 `occt-geom/src/geom_lib.rs`（`To3d` / `isIsoLine` / `buildC3dOnIsoLine` / `BuildCurve3d`），`shape_build_edge.rs` 增 `evaluate_max_segment` + `ShapeBuildEdge::build_curve3d`，`shape_fix_edge.rs` 增 `temp_same_range`、`fix_add_curve3d(edge)` 改为委托 `build_curve3d`（旧体用 `CurveOnSurface` 冒充 3D 曲线，见 `git diff` 的 `-` 行）；`clib.rs` 增 `ElCLib::To3d` 族 10 个函数。② 调用点 `dispatch_wires.rs:301-325`（`!same_range` → Copy+FixAddCurve3d 回写；else → FixAddCurve3d）↔ `ShapeFix_ComposeShell.cxx:3506-3529`。③ 支撑钩子 `Curve2d` 增 `offset_value`/`bspline_weights2d`/`bezier_weights2d`/`bspline_distinct_knots_mults` 并各实现。④ 门禁：`cargo check --offline --all-targets` 0 errors、`step_obj_gates` 5 passed / 0 failed。⑤ `--wirehist`（debug exe）：T0M `1778/1921 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、a3n00 `226/294 1:208 2:9 4:1 6:4 10:4`、acs10 `787/913 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`，与 §9.587 的 OCCT 真值逐项相同。⑥ `--ecensus` T0M `wires>=2` = 106 条（= OCCT `MULTI`），face 628 = `wires=2 edges=[11,1]`、face 635 = `wires=2 edges=[15,1]`（§9.583 记为 1 wire / 12、16 边）⇒ §9.585 的 `FixLoopWire` 族仍生效。",
      note: "UNPORTED（文件头已标）：`GeomLib::To3d` 的有理 Bezier / 周期 BSpline 两臂、`GeomLib::BuildCurve3d` 一般臂的子区间重裁剪与 `NbIntervals` 切点、`TempSameRange` 的周期 pcurve 平移、`TransformPCurve` 的 `uFact != 1` 仿射分支、`ShapeExtend` 状态位。红线：4 个基准模型的 `--wirehist` / `--ecensus` 不得回落（§9.581 的 ④ `MapReShape::apply_impl` 仍不许撤）。",
    },
        { id: "T-105",
      lane: "healing",
      title: "锥面周期退化环移植：`IsPeriodicConicalLoop` + `ShapeFix_Face::FixPeriodicDegenerated`（T0M face 1764 从「1 边 0 V 跨度」恢复到全 V 段并出网格）",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 3,
      actual: 2,
      startedAt: "2026-10-04",
      updatedAt: "2026-10-04",
      completedAt: "2026-10-04",
      blocker: "",
      goal: "按 `ShapeFix_Face.cxx:3018-3259` 补齐 `Perform` 在 `FixMissingSeam`（`cxx:492-498`）之前的那一步 `FixPeriodicDegenerated`（`cxx:486-489`，`myFixPeriodicDegeneratedMode` 默认 `-1` ⇒ 生效）：单 wire 的锥面若绕一周（`|Σ|ΔU| - 2π| <= tol` 且 `|maxU-minU| > 2π - tol`）就补一条 apex 退化边与 apex wire，把面变成 2 wire，交给后续 seam 步缝合。修 §9.590 的空面（face 1764 的 wire 塌成 `edges=[1] nv=1` ⇒ V 跨度 0 ⇒ `invalid discrete range` ⇒ `IMeshData_Failure`）。",
      next: "**已收口（§9.591）**。旁支（未动，见 §9.591-E）：(1) face 1764 端口 `edges=[5]` vs OCCT `edges_per_wire: 4`（多 1 条），面级密度 `38/35` vs OCCT `72/104`，属 `nbU/nbV` 参数域层 ⇒ 归 T-103；(2) `FixMissingSeam` 的 seam 构造段（`ShapeFix_Face.cxx:1899-2330`）仍是既有未移植面（端口走 `ComposeShell` 代替路径）；(3) `IsPeriodicConicalLoop` 只看 `gp_cone`，OCCT 同款（`cxx:3131` 就是 `Geom_ConicalSurface` 门），无需扩面。",
      acceptance: "**已满足**：`cargo check --manifest-path crates/occt-topo/Cargo.toml --offline` 0 errors；`zz_probe_a3n00 data\\occ\\T0M.stp --fstats` FSTAT 行 1777→1778（失败面归零），`face=1764 mv=38 mt=35` 且 bbox 与 OCCT `FSTAT_OCC face=566` 逐位相同；`--ecensus face=1764` 从 `edges=[1] nv=1` 变 `edges=[5] nv=2`；`step_obj_gates` 5/5（T0M 面积比 1.0053，其余 22 个逐位不变）；未改基线 / area_tol / 断言、未新写测试、未加 OCCT 里没有的谓词。",
      write: "crates/occt-topo/src/shhealing/shape_fix_face.rs · crates/occt-topo/src/step/read_topology.rs · crates/occt-topo/examples/zz_probe_a3n00.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.591 · §9.590 · ShapeFix_Face.cxx:482-498 :144 :1737-1741 :3018-3098 :3101-3259 :3226-3232 :3239-3257 · ShapeFix_Root.lxx:101 :34 · ShapeProcess_OperLibrary.cxx:830 · D35",
      dependsOn: ["T-104"],
      evidence: "**2026-10-04，§9.591 实测（本机）**：① 缺口：`ShapeFix_Face::Perform`（`cxx:482-498`）先 `FixPeriodicDegenerated` 再 `FixMissingSeam`，端口只移植了后者尾巴。② 新增 `PeriodicConicalLoop` + `is_periodic_conical_loop`（`shape_fix_face.rs:88-173` ↔ `cxx:3018-3098`，逐边定向 pcurve 累计 `ΔU`/`Σ|ΔU|`/min·max U·V）与 `ShapeFixFace::fix_periodic_degenerated`（`shape_fix_face.rs:261-392` ↔ `cxx:3101-3259`：`Context()->Apply` → 单 wire + `gp_cone` 门 → `VIso(0.0)` 基圆半径 → apex V `= -R/sin α` → apex V 在 loop V 区间内则退 → 2D 支撑线 + `±X` 方向 + 是否翻 sole wire → 退化 apex 边（pcurve 挂本 face、`Degenerated`、`Range(E,0,|ΔU|)`）→ apex wire → `EmptyCopied` 新面 `[sole_wire, apex_wire]` + `Context()->Replace`）。③ `perform_fix_missing_seam` 按 `cxx:482-498` 接线；`fix_missing_seam` 开头补 `myFace = Context()->Apply(myFace)`（`shape_fix_face.rs:417-423` ↔ `cxx:1737-1741`）否则拿不到 2-wire 面。④ 坑：apex 退化边只设 `edge_range` 而 `pcurve_range` 仍 `(-inf, inf)` ⇒ 下游 `check_wire` 判无效、wire 退回 1 边；补 `set_pcurve_range(&apex_edge, face_key, 0.0, |ΔU|)`（`shape_fix_face.rs:371`）对齐 `BRep_TEdge` 单份 range 语义。⑤ 验证：FSTAT 1777→1778，`face=1764 mv=38 mt=35 bbox=(-23.729500,-207.157336,-107.158859)-(-21.270500,-206.058541,-104.833829)` = OCCT `face=566` bbox 逐位相同（V 段恢复 `[-0.717188,+0.717188]`）；逐面 `mt` 求和 65303→**65338**（OCCT `FSTAT_OCC` 求和 66576，缺口 -1.86%）；`step_obj_gates` 5/5。⑥ 接线点：`step/read_topology.rs:766`、`examples/zz_probe_a3n00.rs:337`。",
      note: "OCCT 侧对照：`output/t0m_occt_wires.txt` `FACE 1764 type=2(Cone) wires=1 edges_per_wire: 4 spans: u=6.2832 v=1.4344` + 两条 `SEAM face=1764`；`output/t0m_fs_occt.txt` `FSTAT_OCC face=566 nodes=72 triangles=104`（同 bbox）。红线：4 个基准模型的 `--wirehist` / `--ecensus` / `step_obj_gates` 不得回落。",
    },
        { id: "T-106",
      lane: "meshing",
      title: "motoc 圆柱面多出横向环带（§9.598，13944 三角与 OCCT 逐位相同）+ `Geom2dConvert` / `GeomLib::SameRange` / `Gp*` 基元变换正本清源链（§9.604-§9.623，潜伏缺口）",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 2,
      actual: 1,
      startedAt: "2026-10-05",
      updatedAt: "2026-10-06",
      completedAt: "2026-10-05",
      blocker: "",
      goal: "修用户报告的 motoc 缺陷：`data/occ/motoc.step` 端口渲染截图里大/小圆柱面上各多出一处横向环带（`occ-motoc.obj` 无）。按几何 bbox 配对判定是**网格**多切而非拓扑多面；根因钉到 `BRepMesh_ShapeVisitor::addWire` 无对应分支，删掉端口自造规则，不改偏转、不加谓词。",
      next: "**已收口（§9.598 + §9.601 + §9.604-§9.607 + §9.609-§9.623）**。② 的 4 张面已在本轮补齐：`Geom2d_BSplineCurve::InsertKnots` / `Segment` / `SetOrigin` 三件套落地（OCCT `Geom2d_BSplineCurve.cxx:343 :707 :989`）、`box_bspline_as_curve`（`geom_bnd_lib_bspline2d.rs:77`）子区间臂改成 `Copy` + `Segment` 取裁剪控制网盒、`ReparamCurve2d::bspline_copy2d`（`pcurve_full/surface_projector.rs:3008`）按 `GeomLib.cxx:947-969` + `Geom2dConvert.cxx:347-351` 把 `Copy()` 后的 `Geom2d_BSplineCurve` 交回 `GeomBndLib_Curve2d.cxx:137-157` 的探测。223 面配对 `paired=223 big=0 sumMt=13944`（4→0；分段读数：仅完成前两项时 `big=1`（余 `#106 maxd=0.0920`），补 `ReparamCurve2d` 后 `big=0`）。无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.605（本轮）**：§9.604-E 登记的 `is_bspline2d` 三处入口收口——按 `.cxx` 钉死那四处判据都是「手里那条曲线本身精确是 `Geom2d_BSplineCurve`」（`Geom2dConvert.cxx:420-423`、`GeomLib.cxx:563` + `:603`、`:3010-3029`、`ShapeBuild_Edge.cxx:679-686`），而 `GeomLib.cxx:920-921` / `:961-969` 证明 `SameRange` 从不存包装器、B-spline 输入产出真 `Geom2d_BSplineCurve`；于是把端口的冒充点从「包装器要不要自称 BSpline」移走：`geom_lib_same_range` 的 B-spline 臂改成直接返回真曲线（`pcurve_full/surface_projector.rs:3034` `same_range_bspline2d` + `:3072` `reparam_curve2d` 直返），`geom2d_convert.rs` / `occt-geom/src/geom_lib.rs` / `shhealing/shape_build_edge.rs` **零改动**。实测该入口在这批模型上未触发（改造前后 `reparam_curve2d` 命中数相同：a3n00 0 / T0M 2 / acs10 0 / TDB 0 / motoc 3，全返真 BSpline，且四个分派点零命中）。 **§9.606（本轮）**：§9.605-F 登记的圆锥臂已补——`Geom2dConvert::CurveToBSplineCurve` 的 Circle/Ellipse/Hyperbola/Parabola 臂（裁剪 + 非裁剪）落到 `crates/occt-geom2d/src/geom2d_convert.rs`（`build_conic_bspline2d:51` / `curve_to_bspline_curve_bspl:98` / `trimmed_curve_to_bspline_curve:157`，复用 `occt-core::convert` 与 3D 同构），`surface_projector.rs:3052` 的圆锥分支按 `GeomLib.cxx:947-958` + `Geom2dConvert.cxx:199-209` 夹窗口后调 `Convert_*` 再 `bspl::knots::reparametrize`（`GeomLib.cxx:961-969`）；Bezier 裁剪臂 / `RationalC1` 长弧拼接 / Offset 近似仍 UNPORTED。触发面：五模型 + `step_obj_gates` + `export_data_obj` 的 `conic=0`（pcurve 全为 Line/BSpline，潜伏缺口），基线逐项不变。等跨 Circle 的旋转臂 `GeomLib.cxx:872-888` 在 T-103 的 `wire_fix.rs` 内，登记未动。 **§9.607（本轮）**：§9.606-E 登记的两入口收口——① `Geom2dConvert_CompCurveToBSplineCurve` 全量移植（`crates/occt-geom2d/src/comp_curve_to_bspline.rs`，两 ctor + 公开 `Add` 极点 G0 判据 + 私有 `Add` 度数取齐/`Ratio`/节点权值拼接/`RemoveKnot` 降重数，逐字照 `Geom2dConvert_CompCurveToBSplineCurve.cxx:29-243`；依赖件 `Geom2d_BSplineCurve::increase_degree:360` / `remove_knot:410` 一并补），接上 `Geom2dConvert.cxx:244-262`（Circle）与 `:283-301`（Ellipse）的 `RationalC1 && U2-U1>=6` 对半拼接臂（`geom2d_convert.rs:178-184` / `:195-201`）；② `Geom2d_BezierCurve::Segment`（`bezier_curve.rs:57` = `Geom2d_BezierCurve.cxx:356-391`，`BuildCache`/`Trimming`/`CoefficientsPoles` 往返）落地 + `Geom2dConvert.cxx:323-345` 裁剪 Bezier 臂（`geom2d_convert.rs:218-236`）。触发面 `conic=0 comp=0 bezier=0`（pcurve 无圆锥、无 RationalC1 调用方），潜伏缺口；Offset 近似臂仍 UNPORTED。 **§9.609（本轮）**：Offset 近似臂闭环——`AdvApprox` 机械移入 `occt-core::adv_approx`（`occt-geom` 改 `pub use occt_core::adv_approx`，调用路径不变）、新增 `ApproxAFunction2d`（`Num1DSS=0/Num2DSS=1/Num3DSS=0`）与 `Geom2dConvert_ApproxCurve`（`crates/occt-geom2d/src/convert_approx_curve.rs`），`Geom2dConvert::CurveToBSplineCurve` 的裁剪/非裁剪 Offset 臂（`:353-368` / `:425-440`）接 `Geom2dConvertApproxCurve::new(..., 1e-4, C2, 16, 14)`；`OFFAPPROXDBG` 探针（已撤除）实测 15 个 `data/*.step` + 5 个 `data/occ/*` 全 `OFFAPPROX=0` ⇒ 潜伏缺口。至此 `Geom2dConvert.cxx:181-449` 仅余等跨 Circle 旋转臂 `GeomLib.cxx:872-888`（在 T-103 `wire_fix.rs` 侧）。 **§9.610（本轮）**：`ShapeBuild_Edge::TransformPCurve` 的 `uFact != 1` 仿射分支里 `Geom2d_Conic` 臂（`cxx:660-678`）与兜底臂（`cxx:679-687`）闭环——圆锥臂按 `cxx:667-676` 建 `Geom2dTrimmedCurve(result, aFirst, aLast)` → `Geom2dConvertApproxCurve::new(..., Approximation(), C1, 100, 6)`，`HasResult()` 取 `Curve()` 否则 `curve_to_bspline_curve_bspl(tcurve, QuasiAngular)`，再回写 `aFirst/aLast`；两臂收拢到 `tMatu` 极点缩放（`cxx:688-697`）。触发面：五模型 `--wirehist` 与 `export_data_obj` 逐项不变 ⇒ 潜伏缺口（同 §9.609）。`TransformPCurve` 的 UNPORTED 标记清除，该函数至此无 UNPORTED 段。 **§9.611（本轮）**：`GeomLib::SameRange` 等跨 `Geom2d_Circle` 旋转臂（`GeomLib.cxx:871-888`）闭环——`wire_fix.rs::geom_lib_same_range` 按 `cxx:871-888` 在 `Line` 与 `TrimmedCurve` 之间插入圆臂：`is_geom2d_circle()` 命中时取 `gp_circ2d().location()`，按 `Circ2d().IsDirect()` 定 `dU` 符号（direct → `FirstOn - ReqFirst`，否则反向），`GpTrsf2d::set_rotation(&p, dU)` 后 `c2d.transformed(&trsf)` 返回。为落地该臂正本清源三项底层语义（真实缺口）：`GpAx22d::transform`（`gp_Ax22d.hxx:360-367`，变换位置 + 两方向）、`GpCirc2d::transform`（`gp_Circ2d.hxx:335-344`，`radius*=ScaleFactor()` 取正 + `pos.Transform`，此前只平移位置不转方向）、`GpCirc2d::is_direct`（`gp_Circ2d.hxx:181-184`，`vxdir.crossed(vydir) >= 0`）；`Curve2d` 新增 `is_geom2d_circle()` 精确类判（`Geom2dCircle` 覆写 true，Trimmed/Reparam 包装仍 false，对应 `IsKind(STANDARD_TYPE(Geom2d_Circle))`）。等跨「其余」臂（`cxx:908-921`）仍为 `reparam_curve2d` 占位（登记）。 **§9.612（本轮）**：`SameRange` 两条通用臂闭环——等跨「其余」（`GeomLib.cxx:908-921`）与非等跨通用（`:960-968`）落地 `Geom2dTrimmedCurve → Geom2dConvert::CurveToBSplineCurve（默认 TgtThetaOver2）→ BSplCLib::Reparametrize（去重节点）→ SetKnots`：`wire_fix.rs` 新增 `trim_to_reparam_bspline`，抽出 `unequal_trim_window`（`cxx:929-958`，`same_range_unequal_line` 改用它、逻辑逐字不变）并新增 `same_range_unequal_other`，`geom_lib_same_range` 等跨分支按 `cxx:908-909` 条件（`|LastOn-FirstOn|>PConfusion || |RequestedLast+RequestedFirst|>PConfusion`，后者按 OCCT 逐字为**和**）接新臂；缺口件 `Geom2dBSplineCurve::SetKnots`（= 3D `Geom_BSplineCurve::set_knots`，`CheckCurveData` + `updateKnots`）补齐于 `crates/occt-geom2d/src/bspline_curve.rs`；`BSplCLib::Reparametrize` 复用既有 `occt_core::bspl::knots::reparametrize`（`knots.rs:272`）。`reparam_curve2d` 仅在 OCCT 未赋值的退化子条件保留（登记）。 **§9.613（本轮）**：`ShapeBuild_Edge::TransformPCurve` 的 `Geom2d_BezierCurve` 臂（`ShapeBuild_Edge.cxx:643-655`）改为保留 Bezier 类型——`reference.clone_dyn()` 取解包后 basis（等价 `down_cast`）、逐极点 `t_matu.transforms` 后 `set_poles2d`、`Arc::from(bez)` 返回，不再经 `curve_to_bspline_curve` 转 BSpline（清 §9.610-D(1) 旁支）；`Geom2d_BSplineCurve` 臂（`:687-697`）改用 `reference`（解包后 basis，= `down_cast<Geom2d_BSplineCurve>(result)`）改极点返回，替换原先未改极点且可能返回包装的 `Arc::from(result)`。文件头各臂范围说明同步细化。 **§9.614（本轮）**：`GpCirc2d` 其余变换正本清源——`rotate`/`scale`/`mirror_pnt`/`mirror_ax2d` 改为委托 `pos.*`（`gp_Circ2d.hxx:199-207 :309-315`、`gp_Circ2d.cxx:23-38`）；新增 `GpAx22d::{rotate,scale,mirror_pnt,mirror_ax2d}`（`gp_Ax22d.hxx:335-341 :346-354`、`gp_Ax22d.cxx:29-52`）与 `GpDir2d::{rotate,mirror_ax2d}`（`gp_Dir2d.hxx:434-440`、`gp_Dir2d.cxx:67-77`）；并修正共同根 `GpTrsf2d::set_mirror_ax2d`（矩阵符号按 `gp_Trsf2d.cxx:31-46`：存取负反射阵配 `scale=-1`）。`translate_vec` 原本只移位置点、与 OCCT 等价，未动。 **§9.615（本轮）**：`GpLin2d`/`GpAx2d`/`GpElips2d`/`GpHypr2d`/`GpParab2d` 的 `rotate`/`scale`/`mirror_*`/`transform` 正本清源——`GpAx2d` 新增 `rotate`/`scale`/`mirror_pnt`/`mirror_ax2d`/`transform`（`gp_Ax2d.hxx:155-186`、`gp_Ax2d.cxx:51-77`），`GpLin2d` 五者改为委托 `pos.*`（`gp_Lin2d.hxx:142-186`）；`GpElips2d`/`GpHypr2d`/`GpParab2d` 的 `rotate` 委托 `pos.rotate`、`scale` 改 radii/focal 后委托 `pos.scale`、`mirror_*` 委托 `pos.mirror_*`、`transform` 改用 `pos.transform`（`gp_Elips2d.hxx:245 :355-370 :418-428`、`gp_Hypr2d.hxx:500-513`、`gp_Parab2d.hxx:313-321`）；`GpDir2d` 新增 `mirror_dir2d`（`gp_Dir2d.cxx:108-118`）。 **§9.616（本轮）**：3D 侧同族正本清源——`GpDir` 新增 `rotate`/`mirror_dir`/`mirror_ax1`/`mirror_ax2`/`transform`（`gp_Dir.hxx:496-501`、`gp_Dir.cxx:86-155`）；`GpAx1` 新增 `rotate`/`scale`/`mirror_pnt`/`mirror_ax1`/`mirror_ax2`/`transform`/`translate_vec`（`gp_Ax1.hxx:163-211`、`gp_Ax1.cxx:46-76`）；`GpAx2` 新增 `scale`/`transform`/`mirror_*`/`translate_vec`（`gp_Ax2.hxx:320-354`、`gp_Ax2.cxx:83-124`）；`GpAx3` 新增 `rotate`/`scale`/`transform`/`mirror_*`/`translate_vec`（`gp_Ax3.hxx:265-325`、`gp_Ax3.cxx:82-115`）；`GpLin`/`GpCirc`/`GpElips`/`GpHypr`/`GpParab`/`GpSphere`/`GpCylinder`/`GpCone`/`GpPln`/`GpTorus` 的 `rotate`/`scale`/`mirror_*`/`transform` 改为委托 `pos.*` 并按各 `.hxx` 一次缩放半径（消除此前经 trsf 的二次缩放与镜像变号）。 **§9.617（本轮）**：`GpTrsf` 其余设置器正本清源——`set_mirror_ax1`（`gp_Trsf.cxx:57-71`，`loc` 符号相反）、`set_mirror_ax2`（`cxx:73-86`，`matrix=2*n*n^T-I` + `scale=-1` + `loc=2*(n.o)*n`）、`set_transformation_from_to` 的 `MA1` 改用列构造器（`gp_Mat.cxx:30-43`）并把平移项改为 `loc = loc.add(&ma1_loc)`（原丢弃 `GpXyz::add` 返回值）、新增 `set_displacement`（`cxx:217-240`）与 `get_rotation`（`cxx:387-390`）、`set_scale` 去掉 OCCT 没有的 `S≈1` 短路（`cxx:159-168`）、`gp_Mat::set_rotation` 轴归一化（`gp_Mat.cxx:126`）；`Invert` 的 `PntMirror` 分支（`cxx:406-408` 的 `loc.Reverse()`）与点镜像逆运算不一致，端口保持通用逆（数值正确），登记未改。 **§9.618（本轮）**：`GpTrsf` 剩余设置器与 `GpMat` 缺件正本清源——`GpMat` 新增 `set_cols`/`set_col`/`set_row`/`set_value`/`set_dot`（`gp_Mat.cxx:102-118`）/`set_cross`（`cxx:88-100`）/`power`+`powered`（`cxx:292-334`）；`GpTrsf` 新增 `set_rotation_part`（`gp_Trsf.cxx:111-156`）/`set_translation_part`（`cxx:243-280`）/`set_scale_factor`（`cxx:284-337`）/`set_values`（`cxx:346-385`）/`orthogonalize`（`cxx:862-938`）/`pre_multiply`（`cxx:712-836`）/`power`（`cxx:564-708`）。 **§9.619（本轮）**：`GpQuaternion` 欧拉角/轴角与 2D 对偶正本清源——新增 `euler_sequence.rs`（`GpEulerSequence` 26 变体 + `EulerSequenceParams` + `translateEulerSequence`，`gp_EulerSequence.hxx:37-73` + `gp_Quaternion.cxx:226-296`）；`GpQuaternion` 加 `set_vector_and_angle`/`get_vector_and_angle`（`gp_Quaternion.cxx:81-112`，`w<0` 走 `atan2(-vl,-w)` 得 `[-PI,0]`）/`set_euler_angles`/`get_euler_angles`（`cxx:302-412`）/`stabilize_length`/`normalize_or_stabilize`（`cxx:416-444`），并把 `get_matrix` 改回 `s = 2.0/SquareNorm()`（`cxx:153-187`）；`GpMat2d` 加 `set_cols`/`set_col`/`set_row`/`set_rows`/`pre_multiply`/`power`+`powered`（`gp_Mat2d.cxx:40-89 :144-179`）；`GpTrsf2d` 加 `set_transformation_ax2d`/`set_translation_part`/`set_scale_factor`/`vectorial_part_scaled`/`rotation_part`/`power`/`pre_multiply`/`set_values`/`orthogonalize`（`gp_Trsf2d.cxx:72-83 :86-117 :120-196 :198-213 :216-219 :384-548 :550-674 :676-710 :721-746`）。 **§9.620（本轮）**：通用仿射对偶正本清源——`GpGTrsf` 重写为逐 `.cxx`/`.hxx` 分支（`SetVectorialPart`/`SetTranslationPart`/`SetAffinity(ax1,ratio)`/`SetAffinity(ax2,ratio)`/`SetValue`/`Value`/`SetTrsf`/`IsNegative`/`IsSingular`/`Form`/`SetForm`/`Invert`/`Multiply`/`PreMultiply`/`Power`/`Trsf`，`transforms` 修成 `shape != Other && scale != 1.0` 才乘 `scale`；删非 OCCT 的 12 系数 `set_affinity`）；`GpGTrsf2d` 补 `shape`/`scale` 字段与全部缺件（含 `Trsf2d` 三处正交判据、`transforms_xy` in-place）。 **§9.621（本轮）**：`GpTrsf::GetMat4`（`gp_Trsf.hxx:328-352`）与 `GpGTrsf::GetMat4`/`SetMat4`（`gp_GTrsf.hxx:260-303`）落地（`NCollection_Mat4<T>` 端口无对等类，机械映射为 `[[f64;4];4]` 行主序）；`gp_Trsf`/`gp_Mat`/`gp_GTrsf` 的 `DumpJson`/`InitFromJson` 因端口无 `Standard_Dump` 框架，按门禁标为 UNPORTED 注释（`trsf.rs`/`mat.rs`/`gtrsf.rs`）。 **§9.622（本轮）**：`GpTrsf::set_rotation_ax1` 控制流修正——删除 `tr2*rot*tr1` 组合路径，改为 `gp_Trsf.cxx:90-99` 直算式（`shape=Rotation`、`loc=-R*o+o`）；新增 `set_transformation_quat_vec`（`gp_Trsf.cxx:207-213`）；`set_rotation_quat` 改直写（`cxx:103-109`）。 **§9.623（本轮）**：`GpMat` 抛错分支与 `HVectorialPart`/`VectorialPart` 命名分离——`invert`/`inverted`/`power`/`powered`/`set_rotation` 改 `Result`（对齐 `gp_Mat.cxx:262-264` 与 `gp_XYZ.hxx:367-373` 的抛错）；`GpTrsf` 拆 `h_vectorial_part()`/`vectorial_part()`（含 scale，`hxx:464-482`）、`GpTrsf2d` 同步改名；`GpVec2d::transform` 改用含 scale 的 `vectorial_part()`，`GpDir2d::transform` 补 `cxx:80-105` 四分支并改 `h_vectorial_part()`，`GpDir::rotate`/`transform` 与 `iges.rs` 改 `h_vectorial_part()`，新增 `GpVec::transform`（`gp_Vec.cxx:120-136`）。 **收口（§9.623）**：本卡四条子链全部落地 —— motoc 环带（§9.598）、`Geom2dConvert` 全链（§9.604-§9.607）、`GeomLib::SameRange` 余臂（§9.611-§9.612）、`Gp*` 基元正本清源（§9.614-§9.623）；四模型 `--wirehist` / `step_obj_gates` / `fuse_box_cylinder` / `export_data_obj` / `occt-core` / `occt-geom` 基线逐项不变。",
      acceptance: "**已满足**：`cargo check --all-targets` 0 errors；`step_obj_gates` 5/5；`fuse_box_cylinder_is_closed_solid` 1/1；`output/motoc.obj` 三角数 13944 == `data/occ/occ-motoc.obj` 13944，逐面（容差配对）`|dt|>=2` 面数 15→0；T0M/a3n00/acs10/TDB 的 `--wirehist` 逐项不变；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。**§9.601 补充**：223 面 `--fstats` ↔ OCCT `FSTAT_OCC` 按 bbox 配对，`maxcoord 差 > 0.01` 的面 5→**4**（`#99` 已 `0.0000`，余 4 张全 `Cylinder`：`d=0.198/1.782/0.092/0.087`），`|dt|>=2` `0` 张，端口 `mt` 求和 `13944` = OCCT `triangles` 求和 `13944`；本轮只加注释（未移植登记），基线不变。**§9.604 补充（本轮）**：同一配对 `paired=223 big=0 sumMt=13944`，`maxcoord 差 > 0.01` 的面 **4→0**（`#100/#103/#106/#181` 全 `0.0000`；`boxDist max=0.0000`、`<1e-4: 223`、`<1e-6: 206`），逐面 `mt` 直方图 `0: 223`（`|dt|>=2` 0 张），两侧 `mv/mt` 求和同为 `11713/13944`，`data/occ/occ-motoc.obj f=13944` = `data/output/motoc.obj f=13944`；`cargo check --all-targets`（core/math/geom/geom2d/topo 五 crate）0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、四模型 `--wirehist` 逐项不变。**§9.605 补充（本轮）**：`geom_lib_same_range` 的 B-spline 臂改为返回 `SameRange` 真正持有的 `Geom2d_BSplineCurve`（一处改动，`pcurve_full/surface_projector.rs:3034` + `:3072`；`geom2d_convert.rs` / `occt-geom/src/geom_lib.rs` / `shhealing/shape_build_edge.rs` 零改动），四处分派因此落到同名臂；同一配对 `paired=223 big=0 sumMt=13944`、`boxDist max=0.0000 <1e-4: 223 <1e-6: 206 >0.1: 0`、`mv/mt port-occ sum=0`、`mt diff histogram 0: 223` 逐位不变；`--wirehist` 四模型不变；`export_data_obj` 15/15 `ok`；`cargo check --all-targets` 五 crate 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。**§9.606 补充（本轮）**：圆锥臂（Circle/Ellipse/Hyperbola/Parabola，裁剪 + 非裁剪）落地后，同一配对 `.target-gate/t106/pair.ps1` `paired=223 big=0 sumMt=13944`；五模型 `--wirehist` 逐项不变（a3n00 `226/294 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`、motoc `223/249 1:197 2:26`）；`export_data_obj` 15/15 `ok`；`cargo check --all-targets` 五 crate 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。**§9.607 补充（本轮）**：`Geom2dConvert_CompCurveToBSplineCurve`（`comp_curve_to_bspline.rs`）落地 + `RationalC1` 长弧对半拼接臂 + 裁剪 Bezier 臂后，同一配对 `.target-gate/t106/pair.ps1` `paired=223 big=0 sumMt=13944`；五模型 `--wirehist` 逐项不变（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`、motoc `1:197 2:26`）；`export_data_obj` 15/15 `ok`；`cargo check --all-targets` 五 crate 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。**§9.609 补充（本轮）**：`Geom2dConvert::CurveToBSplineCurve` 的两条 `Geom2d_OffsetCurve` 臂（`:353-368` / `:425-440`）闭环——`AdvApprox` 机械移到 `occt-core::adv_approx`（`occt-geom` 改为 `pub use occt_core::adv_approx`，调用路径不变）、新增 `ApproxAFunction2d`（`Num1DSS=0/Num2DSS=1/Num3DSS=0`）与 `Geom2dConvert_ApproxCurve`（`crates/occt-geom2d/src/convert_approx_curve.rs`），两臂接 `Geom2dConvertApproxCurve::new(..., 1e-4, C2, 16, 14)`；`OFFAPPROXDBG` 探针（已撤除）实测 15 个 `data/*.step` + 5 个 `data/occ/*` 全 `OFFAPPROX=0` ⇒ 潜伏缺口。对拍：`.target-gate/t106/pair.ps1` `paired=223 big=0 sumMt=13944`；`--wirehist` 五模型逐项不变；`export_data_obj` 15/15 `ok`；`cargo check --all-targets` 五 crate 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。**§9.610 补充（本轮）**：`TransformPCurve` 的 `Geom2d_Conic`/兜底臂接线后，`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`、五模型 `--wirehist` 逐项不变；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。**§9.611 补充（本轮）**：等跨 `Geom2d_Circle` 旋转臂接线后，`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`；四模型 `--wirehist` 与基线逐项相同：a3n00 `226/294 18 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 106 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 45 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 138 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` ⇒ 潜伏缺口、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。**§9.612 补充（本轮）**：`SameRange` 等跨「其余」/非等跨通用臂接线后，`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`；四模型 `--wirehist` 与基线逐项相同：a3n00 `226/294 18 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 106 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 45 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 138 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` ⇒ 潜伏缺口、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。**§9.613 补充（本轮）**：`TransformPCurve` 的 Bezier 臂保留 Bezier 类型、BSpline 臂改用解包 basis 后，`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`；四模型 `--wirehist` 与基线逐项相同：a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`。 **§9.614 补充（本轮）**：`GpCirc2d` 变换（`rotate`/`scale`/`mirror_*`）与底层 `GpAx22d`/`GpDir2d`/`GpTrsf2d::set_mirror_ax2d` 正本清源后，`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`；四模型 `--wirehist` 与基线逐项相同：a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` ⇒ 潜伏缺口、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.615 补充（本轮）**：`GpLin2d`/`GpAx2d`/`GpElips2d`/`GpHypr2d`/`GpParab2d` 变换正本清源后，`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`；四模型 `--wirehist` 与基线逐项相同：a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` ⇒ 潜伏缺口、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.616 补充（本轮）**：3D 侧 `GpDir`/`GpAx1`/`GpAx2`/`GpAx3` 与各 `Gp*` 体变换正本清源后，`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`；四模型 `--wirehist` 与基线逐项相同：a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` ⇒ 潜伏缺口、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.617 补充（本轮）**：`GpTrsf` 其余设置器正本清源后，`cargo check`（`occt-core`/`occt-topo`，`--all-targets`）0 error、`occt-core` 单测 290 passed / 2 ignored、`occt-geom` 单测 143 passed、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 ok；四模型 `--wirehist` 与基线逐项相同：a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` ⇒ 潜伏缺口、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.618 补充（本轮）**：`GpTrsf` 剩余设置器（`set_values`/`orthogonalize`/`power`/`pre_multiply`/`set_rotation_part`/`set_scale_factor`/`set_translation_part`）与 `GpMat` 缺件（`set_col*`/`set_row`/`set_value`/`set_dot`/`set_cross`/`power`）补齐后，`cargo check --all-targets`（`occt-core`/`occt-topo`）0 error、`occt-core` 单测 290 passed / 2 ignored、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 ok；四模型 `--wirehist` 与基线逐项相同：a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` ⇒ 全为新增无调用点的方法、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.619 补充（本轮）**：`GpQuaternion` 的 `GetVectorAndAngle`/`SetEulerAngles`/`GetEulerAngles`/`GetMatrix` 与 `GpMat2d`/`GpTrsf2d` 缺件补齐后，`cargo check --manifest-path crates/occt-core/Cargo.toml --all-targets` 0 error、`occt-core` 单测 290 passed / 0 failed / 2 ignored、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1；本批无仓库内调用点，但 `get_matrix` 的 `2/SquareNorm()` 改写经 `GpTrsf::set_rotation(_quat)`/`trsf_ext::interpolate_transforms`/`rwmesh::vrml` 触发，上述四条基线逐项不变（含 STEP 网格 `step_obj_gates`）⇒ 单位四元数路径数值等价、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.620 补充（本轮）**：`GpGTrsf`/`GpGTrsf2d` 正本清源后，`cargo check --all-targets`（`occt-core`/`occt-topo`）0 error、`occt-core` 单测 290 passed / 0 failed / 1 ignored、`fuse_box_cylinder_is_closed_solid` 1/1、`step_obj_gates` 5/5（660.64s）；四模型 `--wirehist` 与基线逐项相同：a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`。唯一外部调用点 `shape_build_edge.rs::transform_pcurve` 的 `t_matu`（`set_affinity(gp::OY2d(), uFact)` + `transforms`）仍走 `shape=Other` 路径 ⇒ 数值不变、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.621 补充（本轮）**：`GpTrsf::GetMat4` 与 `GpGTrsf::GetMat4`/`SetMat4` 落地、`DumpJson`/`InitFromJson` 标 UNPORTED 后，`cargo check --all-targets`（`occt-core`/`occt-topo`）0 error、`occt-core` 单测 290 passed / 0 failed / 1 ignored、`fuse_box_cylinder_is_closed_solid` 1/1、`step_obj_gates` 5/5（645.63s）⇒ 全为新增无调用点方法与注释，无行为改动、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.622 补充（本轮）**：`set_rotation_ax1` 改为 `.cxx:90-99` 直算式后，`cargo check --all-targets`（`occt-core`/`occt-geom`/`occt-topo`）0 error、`occt-core` 单测 290 passed / 0 failed / 1 ignored、`occt-geom` 单测 143 passed / 0 failed、`fuse_box_cylinder_is_closed_solid` 1/1、`step_obj_gates` 5/5（649.77s）、`export_data_obj` 15 ok / 0 err（Cone 195/301 … screw 600/790）、四模型 `--wirehist` 逐项不变（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ `shape` 由 `CompoundTrsf` 变 `Rotation` 在 15 处调用点上无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。 **§9.623 补充（本轮）**：`GpMat::invert`/`set_rotation`/`power` 改 `Result` 上抛、`GpTrsf`/`GpTrsf2d` 的 `HVectorialPart`/`VectorialPart` 分离、`GpVec2d::transform` 改含 scale 分支、`GpDir2d::transform` 补四分支后，`cargo check --all-targets`（`occt-core`/`occt-geom`/`occt-topo`）0 error、`occt-core` 单测 290 passed / 0 failed / 1 ignored、`occt-geom` 单测 143 passed / 0 failed、`fuse_box_cylinder_is_closed_solid` 1/1、`step_obj_gates` 5/5（551.46s）、`export_data_obj` 15 ok / 0 err、四模型 `--wirehist` 逐项不变（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 唯一数值敏感路径（带 scale 的 `CompoundTrsf`）在基准模型上未触发、无回归；无新测试、无基线/area_tol/断言改动、无 OCCT 里没有的谓词。",
      write: "crates/occt-topo/src/meshing/model_builder/wire_builder.rs · crates/occt-topo/src/meshing/wire_order.rs · crates/occt-geom2d/src/bspline_curve.rs · crates/occt-geom2d/src/trimmed.rs · crates/occt-geom2d/src/geom2d_convert.rs · crates/occt-geom2d/src/convert_approx_curve.rs · crates/occt-geom2d/src/comp_curve_to_bspline.rs · crates/occt-geom2d/src/bezier_curve.rs · crates/occt-core/src/adv_approx/ · crates/occt-topo/src/geom_bnd_lib_bspline2d.rs · crates/occt-topo/src/pcurve_full/surface_projector.rs · crates/occt-topo/src/shhealing/wire_fix.rs · crates/occt-topo/src/shhealing/shape_build_edge.rs · crates/occt-core/src/gp/ax22d.rs · crates/occt-core/src/gp/circ2d.rs · crates/occt-geom2d/src/curve.rs · crates/occt-geom2d/src/circle.rs · crates/occt-core/src/gp/trsf2d.rs · crates/occt-core/src/gp/dir2d.rs · crates/occt-core/src/gp/ax2d.rs · crates/occt-core/src/gp/lin2d.rs · crates/occt-core/src/gp/elips2d.rs · crates/occt-core/src/gp/hypr2d.rs · crates/occt-core/src/gp/parab2d.rs · crates/occt-core/src/gp/dir.rs · crates/occt-core/src/gp/ax1.rs · crates/occt-core/src/gp/ax2.rs · crates/occt-core/src/gp/ax3.rs · crates/occt-core/src/gp/lin.rs · crates/occt-core/src/gp/circ.rs · crates/occt-core/src/gp/elips.rs · crates/occt-core/src/gp/hypr.rs · crates/occt-core/src/gp/parab.rs · crates/occt-core/src/gp/sphere.rs · crates/occt-core/src/gp/cylinder.rs · crates/occt-core/src/gp/cone.rs · crates/occt-core/src/gp/pln.rs · crates/occt-core/src/gp/torus.rs · crates/occt-core/src/gp/trsf.rs · crates/occt-core/src/gp/mat.rs · crates/occt-core/src/gp/quaternion.rs · crates/occt-core/src/gp/euler_sequence.rs · crates/occt-core/src/gp/mod.rs · crates/occt-core/src/gp/mat2d.rs · crates/occt-core/src/gp/trsf2d.rs · crates/occt-core/src/gp/gtrsf.rs · crates/occt-core/src/gp/gtrsf2d.rs · crates/occt-core/src/gp/vec.rs · crates/occt-core/src/gp/vec2d.rs · crates/occt-core/src/gp/dir2d.rs · crates/occt-core/src/gp/dir.rs · crates/occt-topo/src/iges.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.598 · §9.601 · §9.604 · §9.605 · BRepMesh_ShapeVisitor.cxx:99-134 · ShapeExtend_WireData.cxx:583-588 · BRepBndLib.cxx:109 · GeomBndLib_BSplineCurve2d.cxx:34-57 · GeomBndLib_Curve2d.cxx:154-157 · Geom2d_BSplineCurve.cxx:707 :343 :989 · GeomLib.cxx:842-969 :947-958 :961-969 :920-921 · Geom2dConvert.cxx:347-351 :199-209 :420-423 :379-397 :425-440 · ShapeBuild_Edge.cxx:679-686 :688-698 · GeomBndLib_BSplineSurface.cxx:30 :103 :118 · §9.606 · Geom2dConvert.cxx:181-449 :211-226 :228-264 :266-303 :305-312 :314-321 :347-351 · Geom2dConvert.hxx:169-171 · Convert_ParameterisationType.hxx · Convert_CircleToBSplineCurve.cxx:47-172 · §9.607 · Geom2dConvert_CompCurveToBSplineCurve.cxx:29-243 · Geom2dConvert_CompCurveToBSplineCurve.hxx:28-68 · Geom2d_BezierCurve.cxx:356-391 · Geom2d_BSplineCurve.cxx:236-292 :410-478 · Geom2dConvert.cxx:244-262 :283-301 :323-345 · ProjLib_ProjectedCurve.cxx:179 · BRepOffset_Inter2d.cxx:1154 · §9.609 · Geom2dConvert_ApproxCurve.cxx:31-203 · Geom2dConvert_ApproxCurve.hxx:29-92 · Geom2dAdaptor_Curve.cxx:285-288 :409-486 :490-573 :453-480 :536-566 · AdvApprox_ApproxAFunction.cxx:364-599 :603-630 :663-956 · §9.610 · ShapeBuild_Edge.cxx:596-700 :660-678 :679-687 :688-697 · Geom2dConvert_ApproxCurve.cxx:111-118 · Geom2d_TrimmedCurve.cxx:98-154 · §9.611 · GeomLib.cxx:863-921 :871-888 :890-900 :908-921 · gp_Circ2d.hxx:181-184 :335-344 :347-355 · gp_Ax22d.hxx:360-367 · gp_Dir2d.hxx:210-213 · gp_XY.hxx:159-163 · gp_Trsf2d.hxx SetRotation · §9.612 · GeomLib.cxx:908-921 :924-969 :929-958 :960-968 · Geom2dConvert.cxx:187-209 :442-444 · Geom2d_BSplineCurve.hxx SetKnots · BSplCLib.cxx:756-798 · occt-geom/src/bspline_curve.rs:91 · §9.613 · ShapeBuild_Edge.cxx:643-655 :685-687 :687-697 :596-700 · §9.614 · gp_Circ2d.hxx:181-184 :199-207 :309-315 :335-344 · gp_Circ2d.cxx:23-38 · gp_Ax22d.hxx:335-341 :346-354 :360-367 · gp_Ax22d.cxx:29-52 · gp_Dir2d.hxx:434-440 · gp_Dir2d.cxx:67-77 · gp_Trsf2d.cxx:31-46 · gp_Pnt2d.cxx:52-79 · §9.615 · gp_Lin2d.hxx:142-186 · gp_Ax2d.hxx:155-186 · gp_Ax2d.cxx:51-77 · gp_Dir2d.cxx:108-118 · gp_Elips2d.hxx:245 :355-370 :418-428 · gp_Elips2d.cxx:65-80 · gp_Hypr2d.hxx:500-513 · gp_Parab2d.hxx:313-321 · §9.616 · gp_Dir.hxx:496-501 · gp_Dir.cxx:86-155 · gp_Ax1.hxx:163-211 · gp_Ax1.cxx:46-76 · gp_Ax2.hxx:320-354 · gp_Ax2.cxx:83-124 · gp_Ax3.hxx:265-325 · gp_Ax3.cxx:82-115 · gp_Circ.hxx:236-266 · gp_Elips.hxx:341-394 · gp_Cylinder.hxx:203-237 · gp_Sphere.hxx:208-242 · gp_Trsf.cxx:57-86 :159-168 :172-194 :217-240 :387-390 :397-425 · gp_Trsf.hxx:389-396 :420-431 :464-476 · gp_Mat.cxx:30-43 :102-118 :122-160 · gp_Mat.hxx:246-254 · §9.618 · gp_Trsf.cxx:111-156 :243-280 :284-337 :346-385 :564-708 :712-836 :862-938 · gp_Mat.cxx:45-100 :292-334 · gp_Mat.hxx:55-107 :127-131 · Standard_Integer.hxx:39-42 · §9.619 · gp_Quaternion.cxx:81-88 :91-112 :153-187 :226-296 :302-357 :362-412 :416-444 · gp_Quaternion.hxx:106-129 · gp_EulerSequence.hxx:37-73 · gp_Mat2d.cxx:40-89 :144-179 · gp_Mat2d.hxx · gp_Trsf2d.cxx:72-83 :86-117 :120-196 :198-213 :216-219 :384-548 :550-674 :676-710 :721-746 · §9.620 · gp_GTrsf.cxx:28-197 :202 · gp_GTrsf.hxx:63-91 :118-136 :140-174 :187-213 :316-343 :347-386 :390-428 · gp_GTrsf2d.cxx:24-204 · gp_GTrsf2d.hxx:54-59 :83-131 :144 :166 :193 :235-283 :287-305 · §9.621 · gp_Trsf.hxx:328-352 :355 :358 · gp_GTrsf.hxx:260-303 :305 · gp_Trsf.cxx:943-960 :964-1009 · gp_Mat.cxx:338 · gp_Mat.hxx:301 · gp_GTrsf.cxx:202 · §9.622 · gp_Trsf.cxx:90-99 :103-109 :207-213 · §9.623 · gp_Mat.cxx:242-279 :283-289 :292-334 · gp_XYZ.hxx:367-373 · gp_Trsf.hxx:254 :464-482 · gp_Trsf2d.cxx:198-213 · gp_Dir.cxx:129-155 · gp_Dir.hxx:496-501 · gp_Dir2d.cxx:80-105 · gp_Vec.cxx:120-136 · gp_Vec2d.cxx:112-130 · gp_GTrsf.cxx:44-58 · §9.611-§9.623",
      dependsOn: ["T-105"],
      evidence: "**2026-10-05，§9.598 实测（本机）**：① 拓扑一致：端口 `TOTAL faces=223`、OCCT `FSTAT_OCC faces=223`，逐面 bbox 可配对，无多出的面 ⇒ 是网格多切。② 面级三角数：14 张平面端盖（半环 z=38 等，端口 17 / OCCT 36）+ 1 张（38 / 70），总差 298；端口 `mesh_t=13646`+298=13944=occ-motoc.obj。③ 根因：`wire_builder.rs::add_wire` 在 `order.perform()` 后若 `order.chain_area()<0` 就 `reverse_chain()`；`chain_area()` 只用边端点算面积，含圆弧边的 face 1 实测 `chain_area=-0.000000`（真实 2D 面积 `wire_area_2d=156.638392`）⇒ `-0.0<0.0` 误判 ⇒ 整链反向 ⇒ 内孔被填 17 扇形三角。OCCT `BRepMesh_ShapeVisitor::addWire`（`cxx:99-134`）CheckOrder 后只逐条 `Edge(signed)` 反转，无整链反转。④ 修复：删该块 + 删无调用者的 `chain_area()`/`reverse_chain()`。⑤ 逐面（容差配对）`|dt|>=2 / |dv|>=2` 面数 **15→0**，5 张 bbox 差 >0.01 的面修复前后相同。⑥ `output/motoc.obj` f=13944 = occ-motoc.obj f=13944。⑦ 基线 `--wirehist`（反转 ON/OFF 逐位相同）：a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`。**2026-10-05，§9.601 实测（残留 A/B 结项）**：⑧ 口径核对（`.cxx`）：`BRepBndLib.cxx:109` 面分支 = `BndLib_AddSurface::Add(BS, BRep_Tool::Tolerance(F), B)` + `BRepTools::UVBounds`（无采样、无额外 `Enlarge`）；`BRepTools.cxx:185` → `BndLib_Add2dCurve.cxx:83` → `GeomBndLib_Curve2d`；`GeomBndLib_BSplineCurve2d.cxx:34-35` 取曲线 `First/LastParameter`，子区间时 `:49 aCurve->Segment(aTrim1,aTrim2)` → `:53 NbPoles()` 裁剪控制网盒 → `:57 Enlarge(theTol)`；`GeomBndLib_Curve2d.cxx:154-157` 对 `Geom2d_TrimmedCurve` 先取 `BasisCurve()`。⑨ `#99`（BSpline 面）同分支存在 ⇒ 已修：新增 `GeomBSplineSurface::distinct_knots_and_mults_u/_v`（`occt-geom/src/bspline_surface.rs:133,139`）替掉 `unique_knots_mults(&knots_u/v)` 冒充 `UKnots()/VKnots()`（`GeomBndLib_BSplineSurface.cxx:103,118`），`geom_bnd_lib_surface3d.rs:136,167,190,206,284,298` 六处接线；实测 `port 99 maxd=0.0000 occface=4 occT=214 portT=214 dt=0`。⑩ 余 4 张 `Cylinder` 无同分支 ⇒ 登记未移植：都只偏一个坐标且端口盒大于 OCCT 盒（`y_max 0.097810 / 1.682421 / -0.008022 / -0.013015` vs OCCT 恒 `-0.100000`），超出量 = `R*dU`（`R=25/18.5/5/26`，`dU=0.0078/0.0963/0.0184/0.00335 rad` ≪ 偏转 `0.836`），即 UV 窗口裁剪侧多走一小段弧；来源是 `geom_bnd_lib_bspline2d.rs::box_bspline_as_curve(:79)` 子区间臂落到 `GeomBndLib_OtherCurve2d::Box`（33 点弦线包络 + `Enlarge`），OCCT 走 `Segment` 精确裁剪网；`occt-geom2d` 无 `InsertKnots/Segment/SetOrigin`（OCCT `Geom2d_BSplineCurve.cxx:343 :707 :989`）⇒ 只补 `geom_bnd_lib_bspline2d.rs:63-78` 的未移植注释（行为不变）。⑪ 残留 B 归零：`|dt|>=2` 配对 `0` 张、两侧三角和 `13944=13944`、`output/motoc.obj f=13944`；T-103 网格侧三文件未动。⑫ §9.601-G 逐边复核（`BNDDBG`）：`FSTAT_OCC` 的 `bbox=` 就是 `BRepBndLib::Add(aFace,aBox,false)`（`occt_probe.cpp:533`）⇒ 与端口同为解析盒、口径一致；四张面的 UV 窗口都严格大于「线边定出的窗口」：`#100 +0.00682 rad`（`×R25=0.17` vs 观测 `0.198`）、`#103 +0.07135`（`×R18.5=1.32` vs `1.78`）、`#106 +0.01840`（`×R5=0.092` vs `0.092` 精确）、`#181 +0.00335`（`×R26=0.087` vs `0.087` 精确）；`#100 e0/e1`、`#103 e4`、`#181 e1` 为 `full=false`（`OtherCurve2d` 包络），`#106 e0` 为 `full=true`（活动极点盒）——同一 `Segment` 缺口的两种外观路径。**2026-10-05，§9.604 实测（本机，本轮收口）**：⑬ 端口（`crates/occt-geom2d/src/bspline_curve.rs`）补齐 `Geom2d_BSplineCurve::insert_knots(:156)` / `insert_knot(:210)` / `set_origin(:223)` / `segment(:281)` + `u_knot_index_range(:119)` / `update_flat_knots(:139)`（= `Geom2d_BSplineCurve.cxx:343 :332 :989 :707`，与 3D `Geom_BSplineCurve` 同构），`bspline_copy2d(:665)` = `Copy()`（`cxx:109-112`）。⑭ `geom_bnd_lib_bspline2d.rs::box_bspline_poles(:27)` = `GeomBndLib_BSplineCurve2d.cxx:33-57`（窗口 → `segment(…, PConfusion)` → 裁剪网盒 → `enlarge(the_tol)`）；`box_bspline_as_curve(:77)` 子区间臂不再落 `GeomBndLib_OtherCurve2d::Box`，改为 `curve.bspline_copy2d()` → `box_bspline_poles`（= `GeomBndLib_Curve2d.cxx:137-157` 的探测 + `BasisCurve()`；`trimmed.rs:188-194` 透传 basis）。⑮ `pcurve_full/surface_projector.rs:3008` 新增 `ReparamCurve2d::bspline_copy2d`：按 `GeomLib.cxx:947-958`（窗口夹到曲线自身范围、塌陷取全长）+ `Geom2dConvert.cxx:347-351`（`Copy` + `Segment`）+ `GeomLib.cxx:961-969`（`BSplCLib::Reparametrize` + `SetKnots`）重建 OCCT 会产出的那条 `Geom2d_BSplineCurve`；无新谓词/容差/阈值。⑯ 对拍（既有基线，`.target-gate/t106/pair.ps1` + `_pair.mjs`）：`paired=223 big=0 sumMt=13944`（上一轮 `big=4`；分段：仅 ⑬⑭ 后 `big=1`（余 `#106 maxd=0.0920 mt=43 ↔ occface=56`），补 ⑮ 后 `big=0`）；`boxDist max=0.0000`、`<1e-4: 223`、`<1e-6: 206`、`mv/mt port-occ sum=0`、`mt diff histogram 0: 223`。⑰ 门禁：`cargo check --all-targets`（core/math/geom/geom2d/topo 五 crate）全 0 error；`step_obj_gates` 5 passed；`--lib fuse_box_cylinder_is_closed_solid` 1 passed；默认 `--wirehist` 四模型与 §9.598 逐项相同。⑱ 仪器：`REPARAMDBG` 打印与 `examples/zz_tmp_t106_pcdump.rs` 已撤除，判据脚本/参照读数收进 `.target-gate/t106/`；未动 T-103 网格侧三文件与 `shape_fix_compose_shell/`、`shhealing/`、`step/read_topology.rs`。UNPORTED（登记）：`ReparamCurve2d::is_bspline2d()` 仍 false，`geom2d_convert.rs:19` / `occt-geom/src/geom_lib.rs:101` / `shhealing/shape_build_edge.rs:319` 在这条链上仍走「非 BSpline」臂。**2026-10-05，§9.605 实测（本机，本轮收口）**：⑲ `.cxx` 口径：四处 `is_bspline2d` 判据 = 对象精确类测试（`Geom2dConvert.cxx:420-423` 非裁剪臂 `C->IsKind` + `C->Copy()`；`GeomLib.cxx:563` `DynamicType()` + `:603` 精确类判、臂体 `:604-631`；`GeomLib.cxx:3010-3029` 走 `theC2D->GetType()==GeomAbs_BSplineCurve`，`Geom2dAdaptor_Curve.cxx:1364-1368`；`ShapeBuild_Edge.cxx:642` Bezier / `:679` `!IsKind(Geom2d_BSplineCurve)` / `:683-686` `down_cast<Geom2d_BSplineCurve>(result)` + `:688-697` 就地改 poles）。⑳ `GeomLib.cxx:842-969` 全文无包装器类型，`NewCurvePtr` 恒为 `Handle(Geom2d_Curve)`；等跨臂 `:908-921` 与非等跨臂 `:924-969` 都收在 `:920-921` / `:961-969` `NewCurvePtr = BS;`，`BS` 来自 `Geom2dConvert::CurveToBSplineCurve`（返回类型 `Handle(Geom2d_BSplineCurve)`，`Geom2dConvert.cxx:181-184`）⇒ B-spline 输入后 OCCT 持有真 `Geom2d_BSplineCurve`。㉑ 移植（一处）：`pcurve_full/surface_projector.rs:3034` 抽出 `same_range_bspline2d`（`bspline_copy2d()` → 窗口夹取 `GeomLib.cxx:947-949` / `Geom2dConvert.cxx:199-209` → 塌陷取全长 `:951-957` → `segment(…, PConfusion)` → `bspl::knots::reparametrize`），`:3004` `ReparamCurve2d::bspline_copy2d` 复用它（行为不变），`:3072` `reparam_curve2d` 先试它、成功就把那条 `Geom2dBSplineCurve` 直接还回、失败才落回包装器。㉒ 影响面探针（`SRDBG`/`SRDBG2` 环境变量门控，已撤除）：改造前后 `reparam_curve2d` 命中数相同且**全部** `bs=true`（真 `Geom2dBSplineCurve`，degree 8）：a3n00 0 / T0M 2 / acs10 0 / TDB 0 / motoc 3；`build_curve3d` 平面臂、`is_iso_line` 入口、`to_3d=None` 三处探针在这 5 个模型上**零命中** ⇒ 该入口是潜伏缺口（改造前那 5 条 SameRange 产出只被 UV-bbox（§9.604-B）与网格采样读到），基线不变属预期。㉓ 对拍：`.target-gate/t106/pair.ps1` `paired=223 big=0 sumMt=13944`；`_pair.mjs` `boxDist max=0.0000 p99=0.0000 <1e-4: 223 <1e-6: 206 >0.1: 0`、`mv/mt port-occ sum=0`、`mt diff histogram 0: 223`、`11713/13944` 两侧一致。㉔ 门禁：`cargo check --all-targets`（`occt-core`/`occt-math`/`occt-geom`/`occt-geom2d`/`occt-topo` 五 crate 逐 crate）全 0 error；`--test step_obj_gates` 5 passed；`fuse_box_cylinder_is_closed_solid` 1 passed；`--wirehist` a3n00 `226/294 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3` 逐项不变；`examples/export_data_obj` 15/15 `ok`（`Shape.step v=6150 f=11372`、motoc 侧 `f=13944`）。㉕ 仪器：`SRDBG`/`SRDBG2` 门控 `eprintln!` 全部撤除，`crates/occt-geom/src/geom_lib.rs` 相对本轮起点零改动，`.target-gate/t106/` 本轮新增的 `port_fstats_after.txt`/`port_fstats.err.txt`/`pair_after.ps1` 已删。UNPORTED（登记，未改行为）：圆锥曲线输入需 `Convert_CircleToBSplineCurve`/`Convert_EllipseToBSplineCurve`（`Geom2dConvert.cxx:379-397`）、偏移曲线需 `Geom2dConvert_ApproxCurve`（`:425-440`、`ShapeBuild_Edge.cxx:667-674`）——无同等分支，继续留在包装器；`GeomLib.cxx:932-943` 周期曲线窗口未按 `IsPeriodic()` 分叉。**2026-10-05，§9.606 实测（本机，本轮收口）**：㉖ `.cxx` 口径：把 `Geom2dConvert::CurveToBSplineCurve`（`Geom2dConvert.cxx:181-449`）的裁剪臂（Line `:211-226` / Circle `:228-264` / Ellipse `:266-303` / Hyperbola `:305-312` / Parabola `:314-321` / Bezier `:323-345` / BSpline `:347-351` / Offset `:353-368`）与非裁剪臂（Ellipse `:379-387` / Circle `:389-397` / Bezier `:399-419` / BSpline `:420-423` / Offset `:425-440`）逐条钉死，默认参数化 `TgtThetaOver2`（`Geom2dConvert.hxx:169-171`；`Convert_ParameterisationType.hxx` 五值语义）。㉗ 移植：`crates/occt-geom2d/src/geom2d_convert.rs` 新增 `build_conic_bspline2d:51` / `curve_to_bspline_curve_bspl:98` / `trimmed_curve_to_bspline_curve:157`（复用 `occt-core::convert::{circle,ellipse,hyperbola,parabola}_to_bspline_curve`，与 3D `occt-geom/src/convert_bspl.rs` 同构；`BSplineCurveBuilder` 的左旋镜像 + `SetTransformation(conic.XAxis(), OX2d)` 见 `cxx:82-97`，端口在 `occt-core/src/convert/conic_curves.rs:56-75` 的 `scaled_poles` 里按同序算）；`crates/occt-topo/src/pcurve_full/surface_projector.rs:3052` 的圆锥分支夹窗口（`GeomLib.cxx:947-958` + `Geom2dConvert.cxx:199-209`）→ `Convert_*` → `bspl::knots::reparametrize`（`GeomLib.cxx:961-969`）。㉘ 触发面（`CONICDBG`/`REPDBG` 门控探针，已撤除）：motoc `--fstats` `reparam=3 conic=0`、T0M `reparam=2 conic=0`、a3n00/acs10/TDB `reparam=0 conic=0`；五模型 `--wirehist`、`--fixms`、`step_obj_gates`、`export_data_obj` 全 `conic=0` ⇒ 基线 pcurve 无圆锥（`Geom2d_Circle`/`Ellipse`/`Hyperbola`/`Parabola`），属潜伏缺口。㉙ 对拍：`.target-gate/t106/pair.ps1 paired=223 big=0 sumMt=13944`；`--wirehist` 五模型逐项不变；`export_data_obj` 15/15 `ok`（Cone 195/301、Cube 24/12、Cylinder 146/140、Extrusion 24/12、HoledPlate 180/128、Offset 712/892、OffsetPlaneHoleEdge 48/32、Shape-1 3343/4336、Shape-2 3105/4792、Shape 6150/11372、Sphere 642/1244、Torus 1369/2592、linkrods 5172/8288、rev 104/92、screw 600/790）；`cargo check --all-targets` 五 crate 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1。㉚ 仪器 `CONICDBG`/`REPDBG` 已撤除，未动 T-103 线。UNPORTED（登记）：裁剪 Bezier `:323-345`、`RationalC1` 长弧拼接 `:237-263`/`:276-302`、Offset 近似 `:353-368`/`:425-440`、等跨 Circle 旋转 `GeomLib.cxx:872-888`。**2026-10-05，§9.607 实测（本机，本轮收口）**：㉛ `.cxx` 口径：`Geom2dConvert_CompCurveToBSplineCurve`（`:29-243`）逐分支钉死（两 ctor `:29/:38`、公开 `Add` 极点 G0 判据 `:80-92` + 闭曲线消歧 `:94-104` + 前后向 `:105-122` + `return false :124`、私有 `Add` 度数取齐 `:133-141` + `Ratio=L1/L2 :152-165` + `After`/`!After` 的 `Ratio/Delta/U_de_raccord :167-184` + 节点 `:186-201`（结合点多重度 = `First.Degree()`）+ 极/权 `:202-218`（权乘 `W_first(NbP1)/W_second(1)`）+ 构造 `:220` + `while(M>0&&Ok){M--;Ok=RemoveKnot(NbK1,M,myTol);}` `:222-229`）；与 3D 版的差异（2D 无 `WithRatio`/`MinM`、G0 用极点非端点、`RemoveKnot` 索引校验不收紧）逐条区分。㉜ `RationalC1 && U2-U1>=6` 确认走对半拼接：`Geom2dConvert.cxx:232-262`（Circle）/`:271-301`（Ellipse）三段分支，`>=6` 时 `Umed=(U1+U2)*.5` + 两次 `Convert_*` + `CCTBSpl(TheCurve1,param)` + `Add(TheCurve2,PConfusion,true)` + `BSplineCurve()`；`rg Convert_RationalC1` 显示真实调用方是 `ProjLib_ProjectedCurve.cxx:179` 与 `BRepOffset_Inter2d.cxx:1154`。㉝ 移植：新文件 `comp_curve_to_bspline.rs`（`new:59`/`from_curve:71`/`from_bspline:83`/`bspline_curve:92`/`into_curve:97`/`clear:102`/`add:109`/`add_pair:181`）+ `bspline_curve.rs:360 increase_degree` / `:410 remove_knot`（= `Geom2d_BSplineCurve.cxx:236-292` / `:410-478`）；`geom2d_convert.rs:178-184`（Circle）/`:195-201`（Ellipse）接对半拼接；`bezier_curve.rs:57 segment`（= `Geom2d_BezierCurve.cxx:356-391`，非有理 dim2 / 有理 dim3 齐次）+ `geom2d_convert.rs:218-236` 裁剪 Bezier 臂（= `Geom2dConvert.cxx:323-345`）。㉞ 触发面（`T104PROBE` 门控探针，已撤除）：五模型 `conic=0 comp=0 bezier=0` ⇒ 本批 pcurve 无圆锥、`surface_projector` 走 `TgtThetaOver2` 不进 `RationalC1`，属潜伏缺口。㉟ 对拍：`.target-gate/t106/pair.ps1 paired=223 big=0 sumMt=13944`；`--wirehist` 五模型逐项不变；`export_data_obj` 15/15 `ok`；`cargo check --all-targets` 五 crate 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1。㊱ `step_obj_gates` 并行跑曾以 `0xffffffff` 异常退出（无断言失败），单线程 `--test-threads=1 --nocapture` 复跑 5 passed（554.86s，逐用例 `ok`）⇒ 资源竞争抖动，非几何回退。㊲ 仪器 `T104PROBE` 四处 `eprintln!` 已撤除；未动 T-103 线。UNPORTED（登记）：等跨 Circle 旋转 `GeomLib.cxx:872-888`（在 T-103 `wire_fix.rs` 侧）。**2026-10-06，§9.609 实测（本机，本轮收口）**：㊳ `.cxx` 口径：裁剪臂 `:353-368`（`Curv->IsKind(Geom2d_OffsetCurve)` -> `CurveToBSplineCurve(C, TgtThetaOver2)`，`C` 为调用方 trim）、非裁剪臂 `:425-440`（`Geom2dConvert_ApproxCurve TheCurve(TheCurve, Tol2d, GeomAbs_C2, 16, 14)`，`HasResult()` 取 `Curve()` 否则 `Standard_ConstructionError`）；共用近似器 `Geom2dConvert_ApproxCurve`（`.hxx:29-92`、`.cxx:31-203`）= `Geom2dConvert_ApproxCurve_Eval`（`cxx:55-107`，order 0/1/2 = `d0/d1/d2`，其它 order 错误码 3）+ `Approximate`（`cxx:127-183`：`NbIntervals(GeomAbs_C2/C3)` + `Intervals` 断点 -> `AdvApprox_PrefAndRec CutTool(CutPnts_C2, CutPnts_C3)` -> `AdvApprox_ApproxAFunction(Num1DSS=0, Num2DSS=1, Num3DSS=0, ...)`（`cxx:137-168`）-> `Poles2d` 打包 `Geom2d_BSplineCurve`（`cxx:172-180`）+ `MaxError(2,1)`（`cxx:183`））+ 访问器（`cxx:185-203`）。㊴ 依赖布局：`AdvApprox` 机械（`approx.rs`/`cutting.rs`/`simple.rs`/`mod.rs`）从 `occt-geom/src/adv_approx/` 移到 `crates/occt-core/src/adv_approx/`（逻辑零改动；`occt-geom2d` 不能依赖 `occt-geom`），`occt-core/src/lib.rs` 加 `pub mod adv_approx;`，`occt-geom/src/lib.rs` 的 `pub mod adv_approx;` 改 `pub use occt_core::adv_approx;`（调用路径不变）。㊵ 移植：新增 `ApproxAFunction2d`（`occt-core/src/adv_approx/approx.rs:201` / `:446`，与 `ApproxAFunction3d` 同一 `Perform`/`Approximation` 控制流，单 2D 子空间 `local_dim=[2]`，`poles`=(x,y)）与 `Geom2dConvert_ApproxCurve`（`crates/occt-geom2d/src/convert_approx_curve.rs`：`stored_curve`=`Geom2dAdaptor_Curve::load`（`:285-288`）递归解包 trim；`adaptor_intervals`=`NbIntervals/Intervals`（`:409-486`/`:490-573`）含 offset `BaseS` 映射 `C0->C1/C1->C2/C2->C3/else CN`（`:453-480`/`:536-566`）；`bspline_from_poles_knots_mults` 用 `knot_sequence` 展平）；两条臂接线 `geom2d_convert.rs:155`（非裁剪）/ `:257`（裁剪，`Geom2dTrimmedCurve::new_sense(basis, u1, u2, true, false)` 代表 `C`）。㊶ 触发面（`OFFAPPROXDBG` 门控 `eprintln!`，已撤除）：15 个 `data/*.step`（`export_data_obj` 全跑）`OFFAPPROX=0`；5 个 `data/occ/*`（`zz_probe_a3n00 --fstats`）逐模型 `OFFAPPROX=0`（a3n00 / T0M / acs10 / TDB / motoc）⇒ 本批 pcurve 无 `Geom2d_OffsetCurve`，潜伏缺口（同 §9.606/§9.607 的 `conic=0 comp=0 bezier=0`）。㊷ 对拍：`.target-gate/t106/pair.ps1` `paired=223 big=0 sumMt=13944`；五模型 `--wirehist` 逐项不变（a3n00 `226/294 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`、motoc `223/249 1:197 2:26`）；`export_data_obj` 15 ok / 0 err（Cone 195/301 … screw 600/790，逐项同 §9.606）。㊸ 门禁：五 crate `cargo check --all-targets` 0 error；`step_obj_gates` 5 passed；`fuse_box_cylinder_is_closed_solid` 1 passed。㊹ 仪器：`OFFAPPROXDBG` 两处 `eprintln!` 已撤除，工作树零残留；未动 T-103 线。UNPORTED（登记）：等跨 Circle 旋转臂 `GeomLib.cxx:872-888`；`Geom2dConvert_ApproxCurve` 的其它 OCCT 调用方（`ShapeBuild_Edge.cxx:667-674`）未接。**2026-10-06，§9.610 实测（本机，本轮收口）**：㊺ `.cxx` 口径：`ShapeBuild_Edge::TransformPCurve`（`ShapeBuild_Edge.cxx:596-700`）在 `uFact != 1` 下按 `result` 类型分派——`Geom2d_Line`（`:624-641`）、`Geom2d_BezierCurve`（`:643-655` 就地改极点并返回 Bezier）、`else` 段：`Geom2d_Conic`（`:660-678`，`new Geom2d_TrimmedCurve(result, aFirst, aLast)`（`:667`）+ `Geom2dConvert_ApproxCurve approx(tcurve, Precision::Approximation(), GeomAbs_C1, 100, 6)`（`:668`）+ `HasResult()` 取 `Curve()` 否则 `CurveToBSplineCurve(tcurve, Convert_QuasiAngular)`（`:669-676`）+ 回写 `aFirst/aLast`）、非 BSpline 兜底（`:679-687` `CurveToBSplineCurve(result, QuasiAngular)`）、BSpline 直取（`:687`），收拢到 `:688-697` 极点仿射。㊻ 移植：`crates/occt-topo/src/shhealing/shape_build_edge.rs:279 transform_pcurve` 新增圆锥臂 `:336-370`（判据复 §9.606 的 `gp_*2d`；`Geom2dTrimmedCurve::new_sense(Arc::from(reference.clone_dyn()), aFirst, aLast, true, false)`；`Geom2dConvertApproxCurve::new(&tcurve, APPROXIMATION, Shape::C1, 100, 6)`；`curve_to_bspline_curve_bspl(&tcurve, QuasiAngular)` 回退；回写 `aFirst/aLast`）与兜底臂 `:371-384`（`curve_to_bspline_curve_bspl(reference, QuasiAngular)`），二者接 `t_matu.transforms` 极点缩放；文件头与函数文档的 UNPORTED 标记清除。㊼ 对拍：`step_obj_gates` 5/5、`--lib fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`；五模型 `--wirehist` 逐项不变（a3n00 `226/294 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 基准 pcurve 在该入口不落 `Geom2d_Conic`，潜伏缺口。㊽ 未动 T-103 线；Bezier 臂仍经 `curve_to_bspline_curve` 返回 BSpline（OCCT 返回 Bezier，几何控制点相同）——登记旁支。**2026-10-06，§9.611 实测（本机，本轮收口）**：㊾ `.cxx` 口径：`GeomLib::SameRange`（`GeomLib.cxx:842-969`）等跨分支（`:863`）按 `IsKind` 分派 `Line`（`:864-870`）、`Geom2d_Circle`（`:871-888`）、`Geom2d_TrimmedCurve`（`:890-900`）、其余（`:908-921`）；圆臂 = `Copy()` + `Location()` + `Circ2d().IsDirect()` 定 `dU` 符号 + `SetRotation(P,dU)` + `Transform`，依赖 `gp_Circ2d::IsDirect`（`gp_Circ2d.hxx:181-184`）、`gp_Circ2d::Transform`（`:335-344`）、`gp_Ax22d::Transform`（`gp_Ax22d.hxx:360-367`，位置 + 两方向）。㊿ 底缺口：端口 `GpAx22d` 无 `transform`、`GpCirc2d::transform` 只 `transforms_xy(point)` 不转方向（`:28` 旧体）、`GpCirc2d` 无 `is_direct`。51 移植：`crates/occt-core/src/gp/ax22d.rs` 加 `GpAx22d::transform`（`point`/`vxdir`/`vydir` 各 `transform`）；`crates/occt-core/src/gp/circ2d.rs` 的 `transform` 改为 `radius*=scale_factor()`（负取正）+ `pos.transform`，并加 `is_direct()`（`pos.x_direction().crossed(pos.y_direction()) >= 0`，用既有 `GpDir2d::crossed`）；`crates/occt-geom2d/src/curve.rs` 加默认 `is_geom2d_circle()->false`，`circle.rs` 由 `Geom2dCircle` 覆写 true（Trimmed/Reparam 仍 false，对应 `IsKind`）。52 接线：`crates/occt-topo/src/shhealing/wire_fix.rs::geom_lib_same_range` 在 `Line` 与 `TrimmedCurve` 之间插圆臂（`is_geom2d_circle()` + `gp_circ2d().location()` + `is_direct()` 定符号 + `GpTrsf2d::set_rotation` + `c2d.transformed`）；等跨「其余」臂 `cxx:908-921` 仍 `reparam_curve2d` 占位（登记）。53 桩值复核：`wire_fix.rs` 原注释 `Equal-span Circle (cxx:872-888) / other: Reparam stand-in...` 已删，改为登记 `cxx:908-921` 未移植。54 对拍：`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`、四模型 `--wirehist` 逐项不变（a3n00 `226/294 18 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 106 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 45 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 138 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 基准 pcurve 在该入口不落 `Geom2d_Circle`，潜伏缺口。55 门禁：`cargo check --all-targets` 0 error。UNPORTED（登记）：等跨「其余」臂 `GeomLib.cxx:908-921`；`GpCirc2d` 的 `translate_vec`/`rotate`/`mirror_pnt`/`mirror_ax2d`/`scale` 仍各自只动位置不更新方向（与 OCCT `gp_Circ2d.hxx` 语义不一致，未本轮复核）。**2026-10-06，§9.612 实测（本机，本轮收口）**：56 `.cxx` 口径：`GeomLib::SameRange` 等跨「其余」臂 `:908-921`（`TC=TrimmedCurve(CurvePtr,FirstOn,LastOn)` `:910-911` → `CurveToBSplineCurve(TC)` `:913` → `Knots=BS->Knots()` `:914` → `Reparametrize(RequestedFirst,RequestedLast,Knots)` `:916` → `BS->SetKnots(Knots)` `:918` → `NewCurvePtr=BS` `:921`；条件 `:908-909` 含**和** `|RequestedLast+RequestedFirst|`）与非等跨通用臂 `:960-968`（`TC` 窗口逻辑 `:929-958`）；`Geom2dConvert::CurveToBSplineCurve` 裁剪分派 `:187-209`（`BasisCurve()` + 非周期夹窗口）、最终 `else` 抛 `Standard_DomainError` `:442-444`。57 缺口件：`Geom2d_BSplineCurve::SetKnots`（去重节点重写扁平序列）2D 端口缺失，3D `occt-geom/src/bspline_curve.rs:91` 已有同构；`BSplCLib::Reparametrize` 去重版已有（`occt_core::bspl::knots::reparametrize`，`knots.rs:272` = `BSplCLib.cxx:756-798`，含 uniform 分支 + `nextafter` 守卫）。58 移植：`crates/occt-geom2d/src/bspline_curve.rs` 新增 `set_knots`（`knots::nb_poles` CheckCurveData + 周期 `knot_sequence_periodic`/非周期 `banded_interp::knot_sequence`）；`wire_fix.rs` 新增 `trim_to_reparam_bspline`（`Geom2dTrimmedCurve::new(c2d,u1,u2)` → `curve_to_bspline_curve_bspl(&tc, DEFAULT_PARAMETERISATION)` → `distinct_knots_and_mults` → `knots::reparametrize` → `set_knots`）、抽出 `unequal_trim_window`（`cxx:929-958`，`same_range_unequal_line` 改用它、逐字不变）、新增 `same_range_unequal_other`，`geom_lib_same_range` 等跨分支按 `cxx:908-909` 接 `trim_to_reparam_bspline`、非等跨分支接 `same_range_unequal_other`；`reparam_curve2d` 仅在退化子条件保留。59 对拍：`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`（Cone 195/301 … screw 600/790 逐项同 §9.606）、四模型 `--wirehist` 逐项不变（a3n00 `226/294 18 1:208 2:9 4:1 6:4 10:4`、T0M `1778/1921 106 1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `787/913 45 1:742 2:25 5:1 6:18 8:1`、TDB `2180/2385 138 1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 基准 pcurve 的 `SameRange` 不落该臂，潜伏缺口。UNPORTED（登记）：`Geom2dConvert::CurveToBSplineCurve` 的 trim-of-trim 落到 OCCT 抛错（端口返 `None` 等价）；`reparam_curve2d` 线性包装器仍在其它调用点使用。**2026-10-06，§9.613 实测（本机，本轮收口）**：60 `.cxx` 口径：`ShapeBuild_Edge::TransformPCurve`（`ShapeBuild_Edge.cxx:596-700`）先解包 `TrimmedCurve` 到 `BasisCurve()`（`:616-618`），Bezier 臂 `:643-655`（`down_cast<Geom2d_BezierCurve>(result)` 返回同一对象 → 逐极点 `tMatu.Transforms` + `SetPole` → `return bezier`），BSpline 臂 `:685-687`（`down_cast<Geom2d_BSplineCurve>(result)`，`result` 已是 basis），收拢到 `:688-697` 极点仿射。61 移植：`crates/occt-topo/src/shhealing/shape_build_edge.rs` Bezier 臂 `:317-326`（`reference.clone_dyn()` + `t_matu.transforms` 极点 + `set_poles2d` + `Arc::from(bez)`，保留 Bezier 类型）、BSpline 臂 `:385-392`（改用 `reference` = 解包 basis，替换原 `Arc::from(result)`），文件头 `:5-8` 同步细化。62 对拍：`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`、四模型 `--wirehist` 逐项不变（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 潜伏缺口。**2026-10-06，§9.614 实测（本机，本轮收口）**：63 `.hxx/.cxx` 口径：`gp_Circ2d::Rotate`（`hxx:199-207`）⇒ `pos.Rotate`，`gp_Ax22d::Rotate`（`hxx:335-341`）既转位置点也转 `vxdir`/`vydir`，`gp_Dir2d::Rotate`（`hxx:434-440`）= 绕原点旋转的 vectorial part；`gp_Circ2d::Scale`（`hxx:309-315`）⇒ `radius *= theS`（负取正）+ `pos.Scale`，`gp_Ax22d::Scale`（`hxx:346-354`）在 `theS<0` 时反转两方向；`gp_Circ2d::Mirror(P)`（`cxx:23-26`）/`Mirror(A)`（`cxx:35-38`）⇒ `pos.Mirror`，`gp_Ax22d::Mirror(P)`（`cxx:29-36`）镜像位置点 + 反转两方向、`Mirror(A)`（`cxx:43-52`）= 两方向 `gp_Dir2d::Mirror(A)`（`gp_Dir2d.cxx:67-77`，`[2A^2-1,2AB;2AB,2B^2-1]`）+ 位置点 `Mirror(A)`；`gp_Circ2d::Translate`（`hxx:224`）只移位置点。共同根 `gp_Trsf2d::SetMirror(const gp_Ax2d&)`（`gp_Trsf2d.cxx:31-46`）存**取负反射阵**配 `scale=-1`（`VectorialPart=matrix*scale` 才是反射），端口此前存正反射阵配 `scale=-1`（净效果为过垂线的镜像）。64 移植：`trsf2d.rs set_mirror_ax2d` 矩阵改 `[[1-2dx^2,-2dxdy],[-2dxdy,1-2dy^2]]` + `scale=-1` + `loc` 按 `cxx:44-45`；`dir2d.rs` 加 `rotate`（`GpTrsf2d::set_rotation(origin,ang)` + `transform`）、`mirror_ax2d`（OCCT 公式）；`ax22d.rs` 加 `rotate`/`scale`/`mirror_pnt`/`mirror_ax2d`；`circ2d.rs` 四变换委托 `pos.*`。65 对拍：`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`、四模型 `--wirehist` 逐项不变（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 端口内这些 2D 变换在基准路径上未被触发，潜伏缺口。UNPORTED（登记）：`GpTrsf2d::set_mirror_ax2d` 符号修正同时改变 `GpPnt2d`/`GpLin2d`/`GpParab2d`/`GpHypr2d`/`GpElips2d::mirror_ax2d` 语义（现与 OCCT 一致，仓库内无调用点）。 **2026-10-06，§9.615 实测（本机，本轮收口）**：66 `.hxx/.cxx` 口径：`gp_Lin2d` 的 `Rotate`/`Scale`/`Mirror`/`Transform` 全部委托 `pos`（`gp_Lin2d.hxx:156 :167 :142 :178`）；`gp_Ax2d::Rotate`（`hxx:155-159`）= `loc.Rotate` + `vdir.Rotate`，`Scale`（`cxx:51-58`）= `loc.Scale` + `S<0` 时 `vdir.Reverse`，`Mirror(P)`（`cxx:60-64`）= `loc.Mirror` + `vdir.Reverse`，`Mirror(A)`（`cxx:73-77`）= `loc.Mirror(A)` + `vdir.Mirror(A.vdir)`（`gp_Dir2d::Mirror(const gp_Dir2d&)`，`gp_Dir2d.cxx:108-118`），`Transform`（`hxx:182-186`）= `loc.Transform` + `vdir.Transform`；`gp_Elips2d::Scale`（`hxx:355-370`）双半径乘 `theS` 取正 + `pos.Scale`，`Transform`（`hxx:418-428`）双半径乘 `|ScaleFactor|` + `pos.Transform`，`Rotate`（`hxx:245`）/`Mirror`（`cxx:65-80`）委托 `pos`；`gp_Hypr2d::Transform`（`hxx:500-513`）/`gp_Parab2d::Transform`（`hxx:313-321`）同理（半径/焦距乘 `ScaleFactor` 取正 + `pos.Transform`）。67 移植：`dir2d.rs` 加 `mirror_dir2d`（`mirror_ax2d` 改调它）；`ax2d.rs` 加 `rotate`/`scale`/`mirror_pnt`/`mirror_ax2d`/`transform`；`lin2d.rs` 五者委托 `pos.*`；`elips2d.rs`/`hypr2d.rs`/`parab2d.rs` 的 `rotate`/`scale`/`mirror_*`/`transform` 分别委托 `pos.*`（此前只 `transforms_xy(point)`，方向基漏掉）。68 对拍：`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`、四模型 `--wirehist` 逐项不变（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 端口内这些 2D 变换在基准路径上未被触发，潜伏缺口。UNPORTED（登记）：3D 侧同族（`GpLin`/`GpCirc`/`GpElips`/`GpHypr`/`GpParab`/`GpSphere`/`GpCylinder`/`GpCone`/`GpPln`）的 `rotate`/`scale` 仍只动位置，未本轮复核（已由 §9.616 收口）。 **2026-10-06，§9.616 实测（本机，本轮收口）**：69 `.hxx/.cxx` 口径：`gp_Dir::Rotate`（`hxx:496-501`）取旋转 h-vectorial part；`gp_Dir::Mirror(Dir)`（`cxx:86-102`）反射阵、`Mirror(Ax1)`（`cxx:104-120`，第 109 行 `C = XYZ.Y()`）、`Mirror(Ax2)`（`cxx:122-127` = 按主方向 Mirror + Reverse）、`Transform`（`cxx:129-155` 的 Identity/Translation/PntMirror/Scale/else 分派）；`gp_Ax1::Rotate`（`hxx:163-167`）/`Scale`（`hxx:181-188`）/`Mirror`（`cxx:46-76`）/`Transform`（`hxx:201-205`）；`gp_Ax2::Scale`（`hxx:320-330`，主方向不变、负因子反转两方向）/`Transform`（`hxx:346-354`，末尾主方向 = `X^Y`）/`Mirror`（`cxx:83-124`）；`gp_Ax3::Rotate`（`hxx:265-270`）/`Scale`（`hxx:282-290`）/`Transform`（`hxx:306-310`）/`Mirror`（`cxx:82-115`）；`gp_Circ::Scale`（`hxx:236-241`）/`Transform`（`hxx:256-266`）= `pos.Scale/Transform` + 半径一次缩放去符号；`gp_Elips`（`hxx:341-394`）、`gp_Cylinder`（`hxx:203-237`）、`gp_Sphere`（`hxx:208-242`）同构。70 移植：`dir.rs` 加 `rotate`/`mirror_dir`/`mirror_ax1`/`mirror_ax2`/`transform`；`ax1.rs` 加 `rotate`/`scale`/`mirror_*`/`transform`/`translate_vec`；`ax2.rs` 加 `scale`/`transform`/`mirror_*`/`translate_vec`；`ax3.rs` 加 `rotate`/`scale`/`transform`/`mirror_*`/`translate_vec`；`GpLin`/`GpElips`/`GpHypr`/`GpParab`/`GpSphere`/`GpCylinder`/`GpCone`/`GpPln`/`GpTorus`/`GpCirc` 的 `rotate`/`scale`/`mirror_*`/`transform`/`translate_vec` 委托 `pos.*` + 一次半径缩放（消除经 trsf 的二次缩放 `r*|s|*s` 与镜像变号）。71 对拍：`cargo check --all-targets` 0 error、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、`export_data_obj` 15/15 `ok`、四模型 `--wirehist` 逐项不变（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 端口内这些 3D 变换在基准路径上未被触发，潜伏缺口。UNPORTED（登记）：`gp_Dir.cxx:104-120` 的 `C = XYZ.Y()` 逐字复现（仓库内无调用点）；`GpAx2::rotate` 仍走既有 `set_rotation_ax1` + `rotate_vector`（语义等价 `gp_Ax2.hxx:301-309`）。**2026-10-06，§9.617 实测（本机，本轮收口）**：72 `.cxx` 口径：`gp_Trsf::SetMirror(const gp_Ax1&)`（`gp_Trsf.cxx:57-71`）= `SetDot(d)` -> `*= -2` -> 对角线 `+1`（即 `I-2*d*d^T`）-> `loc = matrix*o + o` -> `matrix *= -1`（存 `2*d*d^T-I`、`scale=1`、`loc=2*o-2*(d.o)*d`）；`SetMirror(const gp_Ax2&)`（`cxx:73-86`）= `SetDot(n)` -> `*= 2` -> 对角线 `-1`（即 `2*n*n^T-I`）-> `loc = matrix*o + o`，`scale=-1`（有效阵 `scale*matrix = I-2*n*n^T`、`loc=2*(n.o)*n`）；`SetTransformation(FromA1,ToA2)`（`cxx:172-194`）第一段 `SetRows` 装 ToA2 轴，第二段 `gp_Mat MA1(xDir,yDir,zDir)` 是**列**构造（`gp_Mat.cxx:30-43`）；`SetDisplacement(FromA1,ToA2)`（`cxx:217-240`）两段都用 `SetCol`/列构造 + `MA1.Transpose()`；`GetRotation()`（`cxx:387-390`）= `gp_Quaternion(matrix)`；`SetScale`（`cxx:159-168`）无 `S==1` 分支、`|S| <= Resolution` 抛错、`loc = P*(1-S)`；`Invert`（`cxx:397-425`）Identity 不动 / Translation 与 PntMirror `loc.Reverse()` / Scale `scale=1/scale; loc*=-scale` / else `matrix.Transpose(); loc.Multiply(matrix); loc*=-scale`；`gp_Mat::SetRotation`（`gp_Mat.cxx:122-160`）先 `theAxis.Normalized()`。73 缺口：端口 `set_mirror_ax1` 的 `loc` 取 `matrix*o - o`（符号相反）；`set_mirror_ax2` 存 `I-2*n*n^T`+`scale=1`+`loc=o+2*(n.o)*n`（形式与 `+o` 皆差）；`set_transformation_from_to` 的 `MA1` 用 `set_rows`（转置），且 `loc.add(&ma1_loc)` 丢弃 `GpXyz::add`（`&self -> Self`）返回值致平移项丢失；`set_displacement`/`get_rotation` 缺；`set_scale` 多 `S≈1` 短路；`gp_Mat::set_rotation` 未归一化。74 移植：`crates/occt-core/src/gp/trsf.rs` 三者逐字照改 + 新增 `set_displacement`/`get_rotation` + `set_scale` 去短路/改 `<=`；`crates/occt-core/src/gp/mat.rs` `set_rotation` 先按模长归一化。75 对拍：`cargo check --manifest-path crates/occt-core/Cargo.toml --all-targets` 与 `--manifest-path crates/occt-topo/Cargo.toml --all-targets` 均 0 error；`occt-core` 单测 290 passed / 2 ignored；`occt-geom` 单测 143 passed；`--test step_obj_gates` 5 passed（669.66s）；`fuse_box_cylinder_is_closed_solid` 1 passed；`examples/export_data_obj` 15/15 ok（Cone 195/301、Cube 24/12、Cylinder 146/140、Extrusion 24/12、HoledPlate 180/128、Offset 712/892、OffsetPlaneHoleEdge 48/32、Shape-1 3343/4336、Shape-2 3105/4792、Shape 6150/11372、Sphere 642/1244、Torus 1369/2592、linkrods 5172/8288、rev 104/92、screw 600/790）；四模型 `--wirehist` 与基线逐项相同（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 端口内这些设置器在基准路径上未被触发（`set_displacement`/`get_rotation`/`set_mirror_ax1`/`set_mirror_ax2` 仓库内无调用点），潜伏缺口。UNPORTED（登记）：`Invert` 的 `PntMirror` 分支按 `cxx:406-408` 会得 `-2P`（`X'=-X+2P` 的逆应仍是 `loc=2P`），端口用通用逆得 `+loc`（数值正确），未改成 OCCT 分支以免引入该不一致；`SetValues`/`Orthogonalize`/`Power`/`PreMultiply`/`SetRotationPart`/`SetScaleFactor`/`SetTranslationPart` 未移植。**2026-10-06，§9.618 实测（本机，本轮收口）**：76 `.cxx` 口径：`gp_Mat::SetCols/SetCol`（`gp_Mat.cxx:45-78`，1-based 索引）、`SetDot`（`cxx:102-118`，`[i][j]=ref[i]*ref[j]` 对称）、`SetCross`（`cxx:88-100`，反对称 `M*v = ref^v`）、`Power`（`cxx:292-334`：0->Identity、1->no-op、-1->Invert、负幂先 Invert、二进制平方乘）、`gp_Trsf::SetRotationPart`（`gp_Trsf.cxx:111-156`：`IsEqual(gp_Quaternion())` 判有无旋转 + shape 重分类）、`SetTranslationPart`（`cxx:243-280`：`SquareModulus < Resolution` 判零 + shape 重分类）、`SetScaleFactor`（`cxx:284-337`：`unit=|s-1|<=Res`、`munit=|s+1|<=Res` 双标志 + shape 重分类）、`SetValues`（`cxx:346-385`：`gp_Mat(col1,col2,col3)` = 行主序 3x3，`s=det` 取符号后开立方、`M.Divide(s)`、`Orthogonalize()`、`loc=col4`）、`Orthogonalize`（`cxx:862-938`：列 Gram-Schmidt（`V2-=V1*(V2.V1)`、`V3-=V1*(V3.V1)+V2*(V3.V2)`）-> `SetCols` -> 行同法 -> `SetRows`）、`PreMultiply`（`cxx:712-836` 全分支，`matrix.PreMultiply(T.matrix)` = `T.matrix*matrix`）、`Power`（`cxx:564-708` 全分支，`IsOdd(n)=n%2==1`）。77 移植：`crates/occt-core/src/gp/mat.rs` 加 `set_cols`/`set_col`/`set_row`/`set_value`/`set_dot`/`set_cross`/`power`/`powered`；`crates/occt-core/src/gp/trsf.rs` 加 `set_rotation_part`/`set_translation_part`/`set_scale_factor`/`set_values`/`orthogonalize`/`pre_multiply`/`power`（`loc.Multiply(matrix)` 对应端口 `multiply_mat`）。78 对拍：`cargo check --manifest-path crates/occt-core/Cargo.toml --all-targets` 与 `--manifest-path crates/occt-topo/Cargo.toml --all-targets` 均 0 error；`occt-core` 单测 290 passed / 2 ignored；`--test step_obj_gates` 5 passed（647.81s）；`fuse_box_cylinder_is_closed_solid` 1 passed；`examples/export_data_obj` 15/15 ok（Cone 195/301 … screw 600/790）；四模型 `WIREHIST` 与基线逐项相同（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 本批全为新增无调用点方法，潜伏缺口、无回归。UNPORTED（登记）：`GpMat::divide` 未按 `gp_Mat.hxx:353-357` 对 `|scalar|<=Resolution` 抛错；`pre_multiply`/`power` 的 shape 合并沿用端口通用 `multiply`（`CompoundTrsf`）；`gp_Quaternion::GetVectorAndAngle`/`SetEulerAngles`/`GetEulerAngles`（`gp_Quaternion.hxx:106-129`）与 `DumpJson`/`InitFromJson` 仍未移植。**2026-10-06，§9.619 实测（本机，本轮收口）**：79 `.cxx` 口径：`gp_Quaternion::SetVectorAndAngle`（`gp_Quaternion.cxx:81-88`）先 `Normalized()`；`GetVectorAndAngle`（`cxx:91-112`）`vl>Resolution` 时取轴 + `w<0 ? 2*atan2(-vl,-w) : 2*atan2(vl,w)`，否则轴 `+Z`、角 `0`；`GetMatrix`（`cxx:153-187`）因子 `s = 2.0/SquareNorm()`（非单位四元数产出缩放阵）；`translateEulerSequence`（`cxx:226-296`）把 26 个枚举值映到 `(i, isOdd, isTwoAxes, isExtrinsic)`（intrinsic 借 extrinsic 代码、轴序倒排、角度在转换时对调）；`SetEulerAngles`（`cxx:302-357`）非 extrinsic 时 `a=Gamma,c=Alpha`、odd 时 `b=-b`，`isTwoAxes` 走 `values[i]=cj*(cs+sc)` 另一套，odd 再 `values[j]=-values[j]`；`GetEulerAngles`（`cxx:362-412`）用 `GetMatrix()` 的 `M(i,j)` 分量，`16*DBL_EPSILON` 退化阈值，odd 全取反、非 extrinsic 对调 Alpha/Gamma；`StabilizeLength`（`cxx:416-430`）/`Normalize`（`cxx:434-444`）；`gp_Mat2d::SetCol/SetCols/SetRow/SetRows`（`gp_Mat2d.cxx:40-89`，1-based）/`Power`（`cxx:144-179`：`n==1` 空 -> `n==0` 单位 -> `n==-1` 逆 -> 负幂先逆 + 二进制幂）；`gp_Trsf2d` 的 `SetTransformation(const gp_Ax2d&)`（`gp_Trsf2d.cxx:72-83`）/`SetTranslationPart`（`cxx:86-117`）/`SetScaleFactor`（`cxx:120-196`）/`VectorialPart`（`cxx:198-213`，`Scale`/`PntMirror` 走 `SetDiagonal` 折 scale）/`RotationPart`（`cxx:216-219`）/`Power`（`cxx:384-548`，compound 分支前置 `matrix.SetDiagonal(scale*m11, scale*m22)`）/`PreMultiply`（`cxx:550-674` 全分支，含 `Ax1Mirror&&Ax1Mirror` 的 `loc.Multiply(T.scale)` + `scale*=T.scale`）/`SetValues`（`cxx:676-710`，`s=sqrt(|det|)`、`M.Divide(s)`、`Orthogonalize`、`loc=col3`）/`Orthogonalize`（`cxx:721-746`，先列后行 Gram-Schmidt）。80 缺口：`GpQuaternion` 的六个方法 + `GetMatrix` 的单位公式 + `gp_EulerSequence`/`translateEulerSequence` 端口全缺；`GpMat2d` 无 `set_col*`/`set_row*`/`pre_multiply`/`power`；`GpTrsf2d` 无上述九项。81 移植：新增 `crates/occt-core/src/gp/euler_sequence.rs`（`GpEulerSequence` 26 变体 + `EulerSequenceParams{ i,j,k,is_odd,is_two_axes,is_extrinsic }` + `translate_euler_sequence` 逐 `case`；`j/k` 由 `1+(i+(is_odd?1:0))%3` / `1+(i+(is_odd?0:1))%3` 算）；`crates/occt-core/src/gp/quaternion.rs` 加 `set_vector_and_angle`/`get_vector_and_angle`/`set_euler_angles`/`get_euler_angles`/`stabilize_length`/`normalize_or_stabilize` 并把 `get_matrix` 改 `s = 2.0/square_norm()`；`crates/occt-core/src/gp/mat2d.rs` 加 `set_cols`/`set_col`/`set_row`/`set_rows`/`pre_multiply`（`self = other*self`）/`power`+`powered`；`crates/occt-core/src/gp/trsf2d.rs` 加 `set_transformation_ax2d`/`set_translation_part`/`set_scale_factor`/`vectorial_part_scaled`/`rotation_part`/`power`/`pre_multiply`/`set_values`/`orthogonalize`（`GpXY` 用 `added`/`subtracted` 消借用冲突，`GpMat2d::value` 按端口 0-based 索引对齐 OCCT 1-based `Value`）；`crates/occt-core/src/gp/mod.rs` 导出 `euler_sequence`。82 对拍：`cargo check --manifest-path crates/occt-core/Cargo.toml --all-targets` 0 error；`occt-core` 单测 290 passed / 0 failed / 2 ignored；`--lib fuse_box_cylinder_is_closed_solid` 1 passed；`--test step_obj_gates` 5 passed（651.99s）；本批新增方法仓库内暂无调用点，但 `GpTrsf::set_rotation(_quat)`（`trsf.rs:267,275`）/`trsf_ext::interpolate_transforms`（`trsf_ext.rs:43-50`）/`rwmesh::vrml::trsf_from_trs`（`vrml.rs:670`）会读改写后的 `get_matrix`，上述基线条目逐项不变 ⇒ 单位四元数数值等价、无回归。UNPORTED（登记，未改行为）：`quaternion_ext::to_axis_angle` 仍用 `w.acos()`（区间 `[0,PI]`，与 OCCT `GetVectorAndAngle` 的 `[-PI,PI]` 不同；OCCT 语义已由新方法提供，旧扩展函数保留以不动 `phase4_integration`/`from_axis_angle` 调用点）；`GpTrsf2d::invert` 沿用端口通用 `matrix.invert()`（OCCT `gp_Trsf2d.cxx:221-250` 对 `Scale`/`PntMirror` 只反 `loc`）；端口既有 `GpTrsf2d::vectorial_part` 等同 OCCT `HVectorialPart()`，OCCT `VectorialPart()` 另以 `vectorial_part_scaled()` 提供而未改旧名。**2026-10-06，§9.620 实测（本机，本轮收口）**：83 `.cxx` 口径：`gp_GTrsf::SetVectorialPart`（`hxx:118-123`）置 `shape=Other`/`scale=0.0`；`SetTranslationPart`（`cxx:28-43`）按 `CompoundTrsf/Other/Translation` 不动、`Identity`->`Translation`、其余->`CompoundTrsf`；`SetAffinity(Ax1,ratio)`（`hxx:316-328`）`matrix.SetDot(dir)` -> `*=1-ratio` -> 对角 `+ratio` -> `loc=o; loc.Reverse(); loc*=matrix; loc+=o`；`SetAffinity(Ax2,ratio)`（`hxx:333-343`）`*=ratio-1` -> `loc` 用**未加 1** 的矩阵算 -> 对角 `+1`；`SetValue`（`hxx:347-370`）列 4 写 `loc`（`Identity`->`Translation`）、列 1-3 先按旧 `scale` 折入矩阵再置 `Other`/`scale=0`；`Value`（`hxx:374-386`）；`SetForm`（`cxx:150-197`）`s=det` 取符号后开立方、`M.Divide(s)`、`TM=M^T*M-I`、`|TM|>1e-12` 则 `Other` 否则 `CompoundTrsf`；`Invert`（`cxx:44-58`）`Other` 走 `matrix.Invert(); loc*=matrix; loc.Reverse()`、否则 `Trsf().Invert()`+`SetTrsf`；`Multiply`（`cxx:60-78`）`Other` 分支 `loc.Add(T.loc*matrix)`（用乘前矩阵）后 `matrix*=T.matrix`（不动 scale）；`PreMultiply`（`cxx:129-148`）`loc*=T.matrix; loc+=T.loc; matrix=T.matrix*matrix`；`Power`（`cxx:80-127`）`Other` 分支按 `Npower=|N|-1` 二进制幂且**无** `N<0` 预逆，否则 `Trsf().Power(N)`+`SetTrsf`；`Trsf()`（`hxx:416-428`）`Other` 抛错；`Transforms`（`hxx:390-396`）`coord*=matrix; if(!Other && scale!=1) coord*=scale; coord+=loc`。`gp_GTrsf2d` 同构（`cxx:24-204`、`hxx:235-305`），差别在 `SetValue`（`hxx:255-267`）**不**折 scale 直接置 `Other`，`Power`（`cxx:108-110`）**有** `N<0` 预逆，`Trsf2d()` 三处正交判据用 `Precision::Angular()`。84 缺口：3D `GpGTrsf` 端口 `transforms` 未按 `shape != Other` 才乘 `scale`（`Other` 时 `scale=0.0` 会清零）、`multiply`/`invert`/`pre_multiply` 走自造通用路径、`set_vectorial_part`/`set_translation_part`/`set_form` 未重分类、12 系数 `set_affinity` 为 OCCT 不存在；2D `GpGTrsf2d` 缺 `shape`/`scale` 字段与 `set_translation_part`/`set_value`/`set_trsf2d`/`set_vectorial_part`/`value`/`form`/`is_negative`/`is_singular`/`translation_part`/`vectorial_part`/`invert`/`multiply`/`power`/`pre_multiply`/`trsf2d`。85 移植：`crates/occt-core/src/gp/gtrsf.rs` 重写（新增 `set_affinity_ax1`/`set_affinity_ax2`/`set_value`/`value`/`set_trsf`/`is_negative`/`is_singular`/`form`/`set_form`/`translation_part`/`vectorial_part`/`invert`/`inverted`/`multiply`/`multiplied`/`pre_multiply`/`power`/`powered`/`trsf`/`transformed`/`from_matrix_vector`，`transforms` 改 `&self` 且按 `shape != Other && scale != 1.0`，非 `Other` 支路经内部 `to_trsf()` 委托 `GpTrsf`，删除非 OCCT 的 12 系数 `set_affinity` 与仅赋值的 `set_form(f)`）；`crates/occt-core/src/gp/gtrsf2d.rs` 重写（补 `shape`/`scale` 字段与全部缺件，`set_affinity` 置 `shape=Other`/`scale=0.0`，`transforms_xy` in-place + `transformed_xy`，保留 `transforms(&GpPnt2d)->GpPnt2d` 兼容既有调用点内部改走 `transforms_xy`）。86 对拍：`cargo check --manifest-path crates/occt-core/Cargo.toml --all-targets` 与 `--manifest-path crates/occt-topo/Cargo.toml --all-targets` 均 0 error；`occt-core` 单测 290 passed / 0 failed / 1 ignored；`--lib fuse_box_cylinder_is_closed_solid` 1 passed；`--test step_obj_gates` 5 passed（660.64s）；四模型 `--wirehist` 与基线逐项相同（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 唯一外部调用点 `shape_build_edge.rs::transform_pcurve`（`set_affinity(gp::OY2d(), uFact)` + `transforms`，`shape=Other` 路径）数值不变、无回归。UNPORTED（登记，未改行为）：`GpMat::invert` 不抛奇异异常，`GpGTrsf::invert`/`GpGTrsf2d::invert` 的 `Other` 分支改先 `is_singular()` 判 `Err`；`GpGTrsf::Power` 的 `Other` 分支按 `cxx:98-118` 无 `N<0` 预逆（OCCT 原样，未加保护）。**2026-10-06，§9.621 实测（本机，本轮收口）**：87 `.hxx` 口径：`gp_Trsf::GetMat4`（`hxx:328-352`）`Identity` 取 `InitIdentity`，否则行 0-2 逐列填 `Value(1..3, 1..4)`、`SetValue(3,0..2)=0`、`SetValue(3,3)=1`；`gp_GTrsf::GetMat4`（`hxx:260-286`）同构，`SetMat4`（`hxx:288-303`）`shape=Other`/`scale=0.0`、左上 3x3 取 `GetValue(r,c)`、`loc=(GetValue(0,3),GetValue(1,3),GetValue(2,3))`；`gp_Trsf::DumpJson`（`cxx:943-960`）用 `OCCT_DUMP_VECTOR_CLASS` 输出 `Location`（3）与 `Matrix`（9，按 `matrix.Value` 行主序）再 `shape`/`scale`，`InitFromJson`（`cxx:964-1009`）用 `OCCT_INIT_VECTOR_CLASS`/`OCCT_INIT_FIELD_VALUE_*` + `Standard_Dump::Text` 反读；`gp_Mat::DumpJson`（`cxx:338`）、`gp_GTrsf::DumpJson`（`cxx:202`）同理；`gp_Trsf2d`/`gp_GTrsf2d` 无这四个方法（`rg` 无命中）。88 缺口：端口无 `get_mat4`/`set_mat4`，无 `Standard_Dump`/`Standard_SStream`/`TCollection_AsciiString`/`NCollection_Mat4<T>` 对等设施。89 移植：`crates/occt-core/src/gp/trsf.rs` 加 `get_mat4`（复用既有 1-based `value`），`crates/occt-core/src/gp/gtrsf.rs` 加 `get_mat4`/`set_mat4`（复用 §9.620 的 `value`）；`NCollection_Mat4<T>` 机械映射为 `[[f64;4];4]` 行主序（`SetValue(row,col)` 语义），未新建 OCCT 类；`DumpJson`/`InitFromJson` 按门禁「无同等分支即标未移植」写成 UNPORTED 注释（`trsf.rs`/`mat.rs`/`gtrsf.rs`），未自造 JSON 序列化器。90 对拍：`cargo check --manifest-path crates/occt-core/Cargo.toml --all-targets` 与 `--manifest-path crates/occt-topo/Cargo.toml --all-targets` 均 0 error；`occt-core` 单测 290 passed / 0 failed / 1 ignored；`--lib fuse_box_cylinder_is_closed_solid` 1 passed；`--test step_obj_gates` 5 passed（645.63s）⇒ 新增无调用点方法与注释，无行为改动、无回归。UNPORTED（登记，未改行为）：4x4 为行主序（`SetValue(row,col)`），与 `rwmesh/gltf.rs::trsf_from_matrix` 的列主序 glTF 缓冲不同，各自文档化未合并；`GpTrsf` 无 `SetMat4`（与 OCCT 一致，`SetMat4` 只在 `gp_GTrsf`）。**2026-10-06，§9.622 实测（本机，本轮收口）**：91 `.cxx` 口径：`gp_Trsf::SetRotation(const gp_Ax1&, double)`（`cxx:90-99`）`shape=Rotation`/`scale=1`/`loc=o` -> `matrix.SetRotation(dir,Ang)` -> `loc.Reverse()` -> `loc*=matrix` -> `loc+=o`；`SetRotation(const gp_Quaternion&)`（`cxx:103-109`）`shape=Rotation`/`scale=1`/`loc=(0,0,0)`/`matrix=R.GetMatrix()`；`SetTransformation(quaternion R, vec T)`（`cxx:207-213`）`shape=CompoundTrsf`/`scale=1`/`loc=T`/`matrix=R.GetMatrix()`。92 缺口：端口 `set_rotation_ax1` 走 `tr2*rot*tr1` 的 `GpTrsf` 组合，产出 `shape=CompoundTrsf`（OCCT 为 `Rotation`），组合的 `multiplied` 按 `CompoundTrsf` 合并分支，与 `.cxx` 无同等控制流；`set_transformation_quat_vec` 缺；`set_rotation_quat` 经 `set_identity()` 间接设置。93 移植：`crates/occt-core/src/gp/trsf.rs` 的 `set_rotation_ax1` 改 `.cxx:90-99` 直算式（保留 `Result` 签名以不动 15 处调用点），新增 `set_transformation_quat_vec`（`cxx:207-213`），`set_rotation_quat` 改直写（`cxx:103-109`）。94 对拍：`cargo check --all-targets`（`occt-core`/`occt-geom`/`occt-topo`）0 error；`occt-core` 单测 290 passed / 0 failed / 1 ignored；`occt-geom` 单测 143 passed / 0 failed；`--lib fuse_box_cylinder_is_closed_solid` 1 passed；`--test step_obj_gates` 5 passed（649.77s）；`examples/export_data_obj` 15 ok / 0 err（Cone 195/301、Cube 24/12、Cylinder 146/140、Extrusion 24/12、HoledPlate 180/128、Offset 712/892、OffsetPlaneHoleEdge 48/32、Shape-1 3343/4336、Shape-2 3105/4792、Shape 6150/11372、Sphere 642/1244、Torus 1369/2592、linkrods 5172/8288、rev 104/92、screw 600/790）；四模型 `--wirehist` 与基线逐项相同（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ `set_rotation_ax1` 的 15 处调用点（`occt-geom`/`occt-topo`）在 `shape` 由 `CompoundTrsf` 变 `Rotation` 后无回归。UNPORTED（登记，未改行为）：`set_rotation_ax1` 保留 `Result` 返回（OCCT `SetRotation` 不抛错），改签名需触及 15 处调用点；`GpMat::set_rotation` 对零模长轴回退为「不归一化」（OCCT `gp_XYZ::Normalized()` 抛 `Standard_ConstructionError`）。**2026-10-06，§9.623 实测（本机，本轮收口）**：95 `.cxx` 口径：`gp_Mat::Invert`（`cxx:242-279`，`|det| <= gp::Resolution()` 抛 `Standard_ConstructionError` `:262-264`）、`Inverted`（`cxx:283-289`）、`Power`（`cxx:292-334`，`-1` 与负幂走 `Invert`）、`SetRotation`（`cxx:122-158` 经 `gp_XYZ::Normalized()` `gp_XYZ.hxx:367-373` 抛错）；`gp_Trsf::HVectorialPart`（`hxx:254` = `matrix`）与 `VectorialPart`（`hxx:464-482`：`scale==1` 直返、`Scale`/`PntMirror` 只 `SetDiagonal(scale*mii)`、其余 `M.Multiply(scale)`）；`gp_Trsf2d::VectorialPart`（`cxx:198-213` 同构）；`gp_Vec::Transform`（`gp_Vec.cxx:120-136`）默认臂 `VectorialPart()`、`gp_Vec2d::Transform`（`gp_Vec2d.cxx:112-130`）默认臂 `VectorialPart()`；`gp_Dir::Transform`（`gp_Dir.cxx:146`）与 `gp_Dir::Rotate`（`gp_Dir.hxx:500`）用 `HVectorialPart()`、`gp_Dir2d::Transform`（`gp_Dir2d.cxx:80-105`）用 `HVectorialPart()`；`gp_GTrsf::Invert`（`gp_GTrsf.cxx:44-58`）`Other` 分支直调 `matrix.Invert()`。96 缺口：端口 `GpMat::invert` 静默 `1/det`、`set_rotation` 零模长回退不归一化、`power` 非 `Result`；`GpTrsf` 只有一只等同 `HVectorialPart` 的 `vectorial_part`（含 scale 的那只缺），`GpTrsf2d` 把两只叫成 `vectorial_part`/`vectorial_part_scaled`；因名字合并，`GpVec2d::transform` 默认臂取了不含 scale 的那只（对 `CompoundTrsf` 且 `scale!=1` 少乘 scale），`GpDir2d::transform` 用手工补 scale 且缺 `Identity/Translation/PntMirror/Scale` 四分支；`GpVec::transform` 缺；`GpGTrsf::invert` 用了端口自加的 `is_singular()` 预判。97 移植：`crates/occt-core/src/gp/mat.rs`（`invert`/`inverted`/`power`/`powered`/`set_rotation` 改 `Result`）；`trsf.rs`（拆 `h_vectorial_part`/`vectorial_part`，`set_rotation` 传播 `?`，`invert` 用 `inverted()?`）；`trsf2d.rs`（`h_vectorial_part`/`vectorial_part` 改名）；`vec2d.rs`（默认臂改 `vectorial_part()`）；`dir2d.rs`（补四分支 + `h_vectorial_part()`）；`dir.rs`（`rotate`/`transform` 改 `h_vectorial_part()`）；`vec.rs`（新增 `transform`/`transformed`）；`gtrsf.rs`（`Other` 分支改 `matrix.invert()?`）；`crates/occt-topo/src/iges.rs`（`GpDir2d` 变换改 `h_vectorial_part()`）。98 对拍：`cargo check --all-targets`（`occt-core`/`occt-geom`/`occt-topo`）0 error；`occt-core` 单测 290 passed / 0 failed / 1 ignored；`occt-geom` 单测 143 passed / 0 failed；`--lib fuse_box_cylinder_is_closed_solid` 1 passed；`--test step_obj_gates` 5 passed（551.46s）；四模型 `--wirehist` 与基线逐项相同（a3n00 `1:208 2:9 4:1 6:4 10:4`、T0M `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`、acs10 `1:742 2:25 5:1 6:18 8:1`、TDB `1:2042 2:112 3:13 4:3 5:2 6:1 7:4 8:3`）⇒ 唯一数值敏感路径（`form()==CompoundTrsf` 且 `scale != 1`）在上述基准未触发、无回归。UNPORTED（登记，未改行为）：`crates/occt-core/src/convert/advanced.rs` 未被 `convert/mod.rs` 声明为模块（整体不参与编译），其 `transform_normal` 的 `t.inverted().vectorial_part()` 在语法上不成立，本轮未动；`GpDir::transform` 默认臂端口多一句 `if d > RESOLUTION { divide }`（OCCT `gp_Dir.cxx:148-149` 无条件除）；`GpDir2d::transform` 经 `Self::from_xy` 重归一化，零向量时端口保留原值、OCCT 得 `inf/NaN`；`GpMat::powered`/`inverted` 端口暂无调用点。",
      note: "编号说明：用户本轮任务编号为 T-104，但板上 `T-104` 已被 §9.588（`BRepLib::BuildCurve3d` 族）占用，故以 `T-106` 记录，避免覆盖既有卡片。红线：不得为对齐 motoc 加面积比/长度滤边/体积门等 OCCT 不存在的规则。",
    },
  ],
  gates: [
        { id: "G-check", name: "五 crate 编译", cmd: "cargo check --all-targets", baseline: "exit 0", last: "exit 0", status: "green",
      acked: false, owner: "",
      note: "五 crate 逐 crate 跑（无根 workspace）：core · math · geom · geom2d · topo 各 0 error（2026-10-06 收尾轮实测；topo `--all-targets` 0 errors / 306 warnings）",
    },
        { id: "G-lib", name: "topo 单测", cmd: "--lib", baseline: "1281 / 0", last: "1281 / 0", status: "green", acked: false, owner: "",
      note: "2026-10-03：本机 1255 passed / 26 failed，26 个**全是写 `%TEMP%` 的文件 IO 往返用例**（bincaf/step/iges/vrml/xcaf/rwmesh/draw/brep_exchange），失败信息 `PermissionDenied (os error 5)`，都在**写**步骤就挂、同路径 PowerShell 可写，与本轮拓扑改动无关（写 workspace 内的 step_obj_gates 5/5 全过）。2026-09-30 本轮实测 1281 passed / 0 failed（12.91s，T-93 seam 修复接线后）；摘除临时计数器后复跑仍 1281/0",
    },
        { id: "G-gates", name: "STEP→OBJ 门禁", cmd: "--test step_obj_gates", baseline: "5 / 5", last: "5 / 5", status: "green", acked: false,
      owner: "",
      note: "2026-10-06 实测 5 passed（56.53s，§9.617-§9.623 收尾轮）：area_tol/断言未动。2026-10-03 实测 5 passed（31.50s）：只有两张动 —— T0M 188978.93→194286.39（ratio 0.9809→1.0084）、acs10 269261.65→272445.08（0.9928→1.0045），其余 21 个模型（含 a3n00 1.0000、TDB 1.0092、motoc 1.0333、bottom 0.9750）逐位不变；area_tol/断言未动。2026-09-30 实测 5 passed（588.92s 接线后 / 704.02s 清理后复跑），**T-93 seam 修复接线未越出任何 area_tol**。历史：T-93 (a) 后两处基线重定（tests/common/mod.rs 写明来源）：occ/T0M.stp 面积比 0.9995→0.9786 ⇒ area_tol 0.01→0.025、occ/acs10.stp 0.9846→0.9025 ⇒ 0.05→0.12",
    },
        { id: "G-fuse", name: "fuse 闭环实体", cmd: "fuse_box_cylinder_is_closed_solid", baseline: "1 / 1", last: "1 / 1", status: "green",
      acked: false, owner: "",
      note: "2026-10-06 实测 1 passed（1280 filtered out）；§9.617-§9.623 每轮复跑。历史：fuse_box_cylinder_is_closed_solid 曾是 T-102 之前的对齐哨兵（box∪cylinder 必闭合）",
    },
        { id: "G-phase", name: "phase3-10 / 19 / 20", cmd: "--test phaseN_integration", baseline: "全绿", last: "全绿", status: "green",
      acked: false, owner: "",
      note: "4·9·7·5·5·5·8·8·5·5（phase10 8/8，T-88 已闭）；本轮无代码改动，未重跑",
    },
        { id: "G-boss", name: "boss 合并", cmd: "--test bop_builder2_boss", baseline: "1 / 2", last: "1 / 2", status: "red", acked: false,
      owner: "T-03",
      note: "夹具不合法（定案 b），非引擎缺口 —— 已接受的红",
    },
        { id: "G-export", name: "导出总检", cmd: "--example export_data_obj", baseline: "15 / 15", last: "15 / 15", status: "green",
      acked: false, owner: "",
      note: "2026-10-06 收尾轮仍 15/15（§9.617-§9.623 每轮复跑）；内部回归哨兵（旧板记 16/16 为过期）；保真度见 specs/_occt_mesh_gt.md",
    },
        { id: "G-iges", name: "IGES 结构自洽", cmd: "--example iges_check -- <18 模型>", baseline: "18 / 18", last: "18 / 18", status: "green",
      acked: false, owner: "",
      note: "15 个 data/*.step + data/occ/{ATU01038,bottom,top}.step；全部 unreferenced=1（只剩根）",
    },
        { id: "G-core", name: "core / math / geom / geom2d", cmd: "--lib", baseline: "290·215·143·72", last: "290·215·143·72", status: "green",
      acked: false, owner: "",
      note: "逐项与基线相同（2026-10-06 §9.619/§9.622 复跑 occt-core / occt-geom 单测全绿）",
    },
  ],
  // 新条目插到数组**开头**（brief 取前 3 条当最近活动；早于 a29 的看提交历史）
  activity: [
    {
      id: "a85",
      at: "2026-10-06",
      title: "T-106 收口：`Geom2dConvert` 全链 + `Gp*` 基元变换正本清源（§9.604-§9.623，潜伏缺口；五基线逐项不变）",
      tone: "info",
      detail: "§9.604-§9.623 一串**潜伏缺口**闭环 —— 这些分支在 23 个基准模型上**不改变任何读数**（四模型 `--wirehist` / `step_obj_gates` 5/5 / `fuse_box_cylinder` 1/1 / `export_data_obj` 15/15 逐项不变），但缺了它们某些输入会落到错误的臂上、失配无法归因到 OCCT 控制流。关键修复：① `geom_lib_same_range` 的 B-spline 臂改为返回 `SameRange` **真正持有**的 `Geom2d_BSplineCurve`（§9.605；旧端口用 `ReparamCurve2d` 包装器冒充，后接三处分派全落错臂）；② 圆锥臂走 `Geom2dConvert::CurveToBSplineCurve`（§9.606，`Geom2dConvert.cxx:379-397`）；③ `Geom2dConvert_CompCurveToBSplineCurve`（`cxx:244-262 :283-301` 长弧对半拼接 + `:323-345` 裁剪 Bezier 臂，需 `Geom2d_BezierCurve::Segment`，§9.607）；④ `Geom2dConvert_ApproxCurve` + 2D `AdvApprox_ApproxAFunction`，接进 `TransformPCurve` 的 `Geom2d_Conic` / `Geom2d_OffsetCurve` 臂（§9.609/§9.610）；⑤ `SameRange` 等跨 `Geom2d_Circle` 旋转臂（`GeomLib.cxx:871-888`）与「其余」臂（`:908-921 :960-968`，`TrimmedCurve → CurveToBSplineCurve → Reparametrize`，§9.611/§9.612）；⑥ `GpAx22d`/`GpCirc2d`/`GpLin2d`/`GpAx2d`/`GpElips2d`/`GpHypr2d`/`GpParab2d`/`GpDir[2d]`/`GpAx[123]`/`GpTrsf[2d]`/`GpMat[2d]`/`GpQuaternion`/`GpEulerSequence`/`GpGTrsf[2d]`/`GpVec[2d]` 逐条对 `.cxx`/`.hxx` 正本清源（§9.614-§9.623），`GetMat4`/`SetMat4` 落地、`DumpJson`/`InitFromJson` 登记 UNPORTED，`GpMat` 抛错分支与 `HVectorialPart`/`VectorialPart` 命名分离。效果：motoc 223 面 `maxcoord 差 > 0.01` 的面 **4 → 0**（§9.604）；余 4 张 `Cylinder` 面的 bbox 差 = 2D pcurve 子区间臂未移植 `Geom2d_BSplineCurve::Segment`。",
      ref: "specs/_a3n00_gap_analysis.md §9.604 · §9.605 · §9.606 · §9.607 · §9.609-§9.623"
    },
    {
      id: "a84",
      at: "2026-10-06",
      title: "T-103 网格塌陷入口闭环：ON-only 塌陷面 8/8 归位 —— 2D pcurve 求值误用边的 3D range，而非 `CurveOnSurface` 的 COS range（§9.608）",
      tone: "success",
      detail: "§9.608。根因对着 `.cxx`：2D pcurve 求值用了**边的 3D range**，而 `BRep_Tool::CurveOnSurface`（`BRep_Tool.cxx:327-361`，`cxx:353` `GC->Range(First,Last)`）返回的是 **COS 表示自身 range**；入口链 `ShapeAnalysis_Edge::PCurve`（`ShapeAnalysis_Edge.cxx:192-208`）→ `ShapeAnalysis_Wire::CheckOrder`（`ShapeAnalysis_Wire.cxx:593-651`，`cxx:634` / `:648-649`）。拆缝子边的 range 被 `ShapeBuild_Edge::CopyRanges`（`ShapeBuild_Edge.cxx:206-334`）重标到切片窗口，用父边 3D range 求值就偏一个周期 ⇒ 分类器多边形自重叠（§9.603 的 `collect_boundary_uv` 非简单多边形）。修复：`meshing/model_builder/wire_builder.rs` 新增 `pcurve_and_range`（返回 COS range，缺则回退 3D range），5 处 pcurve 消费点统一改用。数字（T0M 单二进制 A/B，默认路径 diff=0）：ON `pf=28 1→74`、`pf=812 20→70`、`pf=1663 38→73`，`ON-only` 面数 56→54，ON 合计 `mt 65331→65536`（OCCT 66640）。门禁：`cargo check` 0 err、`step_obj_gates` 5/5、`fuse_box_cylinder_is_closed_solid` 1/1、四模型默认 `--wirehist` 逐项不变。",
      ref: "specs/_a3n00_gap_analysis.md §9.603 · §9.608"
    },
    {
      id: "a83",
      at: "2026-10-05",
      title: "T-104 闭环：motoc 圆柱面「横向 2 层环带」= `add_wire` 里 OCCT 不存在的「按整链有向面积反转」误触发（§9.598，三角数逐位同 OCCT 13944）",
      tone: "success",
      detail: "§9.598。用户截图的蓝圈环带：按 bbox 配对 `TOTAL faces=223` 两侧一致 ⇒ 多的是**网格**不是面。根因是 `meshing/model_builder/wire_builder.rs::add_wire` 在 `CheckOrder` 之后有一段**OCCT 不存在**的后处理：单 wire 面若 `order.chain_area() < 0.0` 就 `order.reverse_chain()`；`chain_area()`（`meshing/wire_order.rs`）只用每条边的**端点**算有向面积，含圆弧的单 wire 面上端点共弦 ⇒ 面积退化到 `-0.0`，浮点噪声令判定成立 ⇒ 整链反向 ⇒ 14 张平面端盖面（半环面/环面）内孔被填 17 个扇形三角（14×19+32=298 = 13944−13646）。对着 `.cxx`：`BRepMesh_ShapeVisitor::addWire`（`BRepMesh_ShapeVisitor.cxx:99-134`）在 `CheckOrder`（`cxx:103`）之后**只**读 `Ordered(i)` 并经 `ShapeExtend_WireData::Edge(signed)`（`ShapeExtend_WireData.cxx:583-588`）逐条反转边、无整链反转分支。落地方案是**删除**该分支（不改成「更好的启发式」）。结果：motoc 三角数与 `data/occ/occ-motoc.obj` **逐位相同（13944）**；四模型 `--wirehist` 与 `step_obj_gates` 5/5 不动。",
      ref: "specs/_a3n00_gap_analysis.md §9.598 · §9.601"
    },
    {
      id: "a82",
      at: "2026-10-04",
      title: "T-105 闭环：移植 `FixPeriodicDegenerated` + `IsPeriodicConicalLoop`，T0M face 1764 从「1 边 0 跨度」恢复到全 V 段并出网格（§9.591）",
      tone: "success",
      detail: "§9.590 把「空面」第一现场钉到 T0M face 1764：wire 塌成 `edges=[1] nv=1`、V 跨度 0 ⇒ `invalid discrete range`。对着 `.cxx` 定位缺口：`ShapeFix_Face::Perform`（`ShapeFix_Face.cxx:482-498`）在 `FixMissingSeam` **之前**先调 `FixPeriodicDegenerated`（`cxx:486-489`，默认 `myFixPeriodicDegeneratedMode=-1` ⇒ `NeedFix==true`，`ShapeFix_Root.lxx:101`；STEP 的 `ShapeProcess_OperLibrary.cxx:830` 不会关掉它），端口只移植了尾部。移植内容：`PeriodicConicalLoop`/`is_periodic_conical_loop`（`cxx:3018-3098`）、`fix_periodic_degenerated`（`cxx:3101-3259`，含 apex V `=-R/sin(alpha)`、退化 apex 边、`EmptyCopied` 新面 + `Context()->Replace`）、`fix_missing_seam` 开头补 `myFace = Context()->Apply(myFace)`（`cxx:1737-1741`）。一个坑：退化 apex 边除 `set_edge_range` 外还要 `set_pcurve_range`（对应 `cxx:3232` 的 `Range(E,0,|dU|)`），否则下游判该边范围无效、wire 又退回 1 边。验收：`--fstats` 1777→1778（失败面归零），face 1764 bbox 与 OCCT `FSTAT_OCC face=566` 逐位相同、`--ecensus wires=1 edges=[5] nv=2`；`cargo check` 0 err、`step_obj_gates` 5/5，T0M 面积比 1.0053（tol 0.025 内）。旁支（未动）：端口 `edges=[5]` vs OCCT `4`、面级密度约半（参数域 `nbU/nbV` 层）。",
      ref: "specs/_a3n00_gap_analysis.md §9.590 · §9.591"
    },
    {
      id: "a81",
      at: "2026-10-04",
      title: "T-104：`BRepLib::BuildCurve3d` / `CheckSameRange` / `ShapeFix_Edge::FixAddCurve3d` 族移植（§9.588，不再用 `CurveOnSurface` 适配器冒充 3D 曲线）",
      tone: "success",
      detail: "§9.588。旧端口在 `FixAddCurve3d` 里用 `meshing::edge_discret::CurveOnSurface`（pcurve + 曲面适配器）**冒充** 3D 曲线；本轮改成真实 `BRepLib::BuildCurve3d`：平面臂 `GeomLib::To3d`、一般臂 `GeomLib::BuildCurve3d`（`GeomLib.cxx:1051-1165`），走 `BRepLib.hxx:90-94` 的默认 `GeomAbs_C1` / `MaxDegree=14` / `MaxSegment=0` + `evaluateMaxSegment`；同族 `BRepLib::CheckSameRange`/`SameRange`（`BRepLib.cxx:149-183 :187-263 :275-295`）、`ShapeFix_Edge::TempSameRange`/`FixAddCurve3d(edge)`（`ShapeFix_Edge.cxx:618-638`，`!SameRange` 时先 `TempSameRange` 再 `ShapeBuild_Edge::BuildCurve3d`，`ShapeBuild_Edge.cxx:714-775`）。旧做法拓扑上「已有 3D 曲线」表面成立，但曲线不是 OCCT 会产出的那条（`SameRange`/容差/后续离散都不同源）。验收：4 个基准模型的 `--wirehist` / `--ecensus` 与 §9.587 逐项相同 ⇒ 纯对齐性补齐、不改已对齐读数。",
      ref: "specs/_a3n00_gap_analysis.md §9.588"
    },
    {
      id: "a80",
      at: "2026-10-04",
      title: "T-102 收口：T0M 面/wire 结构与 OCCT 逐项相同（1778 面 / 1921 wire / 1:1672 2:92 3:8 5:1 6:1 7:2 8:2）",
      tone: "success",
      detail: "§9.587。根因对着 cxx：`ShapeFix_ComposeShell.cxx:629` 的 `TopoDS_Iterator aIt(wire)`（cumOri 默认 true）已把 wire 朝向与子边朝向复合一次、`cxx:635 sbwdM->Add(E)` 存的就是这个复合值，端口 `edges_of_wire`（`cumulated_children`）语义相同 ⇒ `load_wires.rs::wire_data_edges` 在 REVERSED 时再逐边 `reverse()` 是二次反转，`cxx:586-607` 的 2D `WireOrder` 因此拿到反向 pcurve（`Shifted [4,3,2,1]` 而非 `Same [1,2,3,4]`）。修复即直接返回 `edges_of_wire(wire)`。验收：T0M `--wirehist` = `TOTAL faces=1778 wires=1921` + `1:1672 2:92 3:8 5:1 6:1 7:2 8:2`（= OCCT `--wires` 真值，修复前 1775/1918/`1:1669 …`）；§9.585-C 的两条残差同时归位（4 张环面碎片 = OCCT faceid 1462/1503/1648/1668 端口都在：`--ecensus` 0-based 1461/1502/1647/1667 `wires=1 edges=[2]`；`surf=1300` 的 face=1758 回到 `wires=1 edges=[9]` = OCCT）。A/B（23 模型）：只有 T0M 动（1.0055→1.0053），其余 22 个逐位不变（a3n00 1.0000 / acs10 0.9976 / TDB 1.0110），两态 `step_obj_gates` 5/5；`--test 'phase*'` 61 passed；`cargo check --all-targets` 0 errors；`LOADDBG`/`COMPDBG`/`BWDBG`/`SEAMDBG` 插桩全部撤除。剩余网格侧残差另立 T-103（T0M 三角数 +6.0%、薄面簇 1 张不出网格、4 张 `MULTI` 盒值差 9e-3）。",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.587"
    },
    {
      id: "a79",
      at: "2026-10-02",
      title: "T-101 收口：a3n00 面积比 0.9963 → 1.0000、WIREHIST 与 OCCT 逐项相同（1:208 2:9 4:1 6:4 10:4）",
      tone: "success",
      detail: "§9.580。两处新落地的 .cxx 分支：③ `wire_fix.rs::fix_dummy_seam` 改用 `ShapeAnalysis_Edge::FirstVertex/LastVertex` 的复合朝向口径（`first_vertex`/`last_vertex`；TopExp.cxx:214-253 + TopoDS_Iterator.cxx:26-83，cxx:4230-4233）；④ `reshape.rs::MapReShape::apply` 补 `ShapeBuild_ReShape::applyImpl`（cxx:200-324，`until=TopAbs_SHAPE`）的 **EDGE 重建**分支 + `in_flight` 环保护 —— 没有它时 SplitWire/SplitByLine 用 `Replace(prevV,fV)` 统一端点后，用了旧顶点的其它边不会被重建，wire 不闭合，F170 仍是 2 wire/2 面。验收：a3n00 `WIREHIST 1:208 2:9 4:1 6:4 10:4`（与 OCCT 逐项相同）、`our=227126.78 occ=227130.30 ratio=1.0000`、F170 `wires=1 nEdges=8`；`step_obj_gates` 5/5、`cargo check` 0 errors、T170 探针全部撤除。代价（A/B 实测，`output/_probe` 同口径）：只 ④ 使 T0M 0.9931→0.9809、acs10 0.9984→0.9928（都在 area_tol 内，门禁不劣化）；逐面看 T0M 有 30+ 张薄面 `@x∈[5,16] y∈[-217,-204] z∈[104,117]` 从 `mv/mt=9/7` 塌到 `4/2` ⇒ 下一个可对拍点在下游 UV/网格侧，不在 ④ 本身。",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.580"
    },
    {
      id: "a78",
      at: "2026-10-02",
      title: "T-101：F113 恢复为 1 面 / 2 wire（22+6）；a3n00 面积比 0.8996 → 0.9963，未网格 1 → 0（达成）",
      tone: "success",
      detail: "§9.579。两处 .cxx 分歧同时在场才成立：① `wire_fix.rs::copy_replace_vertices_with` 把 `EmptyCopied()` 的朝向与位置落定移到两次 `B.Add(E,V)` **之前**（ShapeBuild_Edge.cxx:103 → :107-114；`TopoDS_Builder::Add` 复合父朝向，TopoDS_Builder.cxx:74-91 ⇒ REVERSED 边存的是 V1=REVERSED/V2=FORWARD；`TopExp::Vertices(E,V1,V2,true)` 再复合一次，TopExp.cxx:214-253 + TopoDS_Iterator.cxx:26-83 ⇒ FirstVertex=V1/LastVertex=V2 与朝向无关）。改前 SplitWire 造出的边在 REVERSED 时首末互换。② `make_faces_on_patch.rs` 三处 `Perform` 改 `perform_recadre(p,false)`（cxx:3100/3110/3223 全部显式传 RecadreOnPeriodic=false；只有 PerformInfinitePoint cxx:3132/3181 走周期搜索）。验收：F113 `FDUMP wire[0] nEdges=22 + wire[1] nEdges=6`、`mv=227 mt=227`；a3n00 ratio **0.8996→0.9963**（our 204327.72→226285.35）、`STATMAP matched=226 unmatched=0`、WIREHIST 1:207 2:10 4:1 6:4 10:4；T0M 未网格 6→2（ratio 0.9931，tol 0.025 内）；acs10 0.9846→0.9984；`--lib` 1281/0、`step_obj_gates` 5/5、`cargo check` 0 errors。全部探针反向撤除，基线/area_tol/断言未动。",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.579"
    },
    {
      id: "a77",
      at: "2026-10-02",
      title: "T-101：goal（160 轮）终态 —— 全链判定忠实、差异锁到「段集合（BREAK-IN）」，门禁 5/5",
      tone: "info",
      detail: "§9.493-§9.578 逐段对拍，把差异逼到一点：93 段「首末顶点同一、逐边 LastVertex 是另一 TShape」（edge_last 命中 93 次，与 same=true 段数完全吻合）⇒ IsShortSegment 恒 0 ⇒ shorts 全 0（154/154）⇒ 合并循环空转 ⇒ 8+14 不成 22 ⇒ [6,7,5,4,3,3] ⇒ 6 wire/5 面 ⇒ Shell(5) 被 reader 丢弃 ⇒ F113 空面。四处负结果/自我更正入库：输入未统一（61/61 same=true ⇒ 否，§9.571）、MapReShape 缺 EDGE 深度重建（补上编译通过但 F113 仍 Shell(5) ⇒ 无收益已还原，§9.569，实现原文已具备可直接重放）、内部顶点未绑定（收集与 cxx:1006-1014 一致 ⇒ 0 次触发是规格行为，§9.576）、93 次 edge_last 失败（多边段上亦为规格行为，§9.577）。终态：门禁 5/5、a3n00 0.8996（未上升）、STATMAP 225/1/11941、apply/probe 插桩全部撤除、库 diff 空、树净。",
      ref: "specs/_a3n00_gap_analysis.md §9.578"
    },
    {
      id: "a76",
      at: "2026-10-02",
      title: "T-101：切分侧与 WireSegment 访问器逐条对照完成（全部忠实）",
      tone: "info",
      detail: "§9.534-§9.536 SplitByGrid（UVBounds/Bounds/TOLINT/closed 段位移/U 线循环）逐条一致；§9.545-§9.551 用无歧义探针确定「产生两条零长度段的那次切分 = U-else 分支 pos=UJointValue(1)=π/cut_index=1」（29/29 同型），并实测 split_wire 守卫① |curr−last|=0 精确命中 ⇒ 不切（78/78）⇒ 零长度段是输入端带下来的；§9.560-§9.563 SplitWire 建顶点三分支、context().apply/context_mut().replace 的 14 处调用点、copy_replace_vertices_with 均与 .cxx 同义；§9.575-§9.577 a_nm_vertices 收集（↔ cxx:1006-1014）与 WireSegment::first_vertex/last_vertex（↔ ShapeFix_WireSegment.cxx:89-101）一致。",
      ref: "specs/_a3n00_gap_analysis.md §9.577"
    },
    {
      id: "a75",
      at: "2026-10-02",
      title: "T-101：IsShortSegment 两处朝向语义修正落地（门禁 5/5、基线逐字不动）",
      tone: "info",
      detail: "§9.489-§9.492：对拍 OCCT IsShortSegment(ShapeFix_ComposeShell.cxx:2394-2447) 发现端口 helpers.rs:347-390 两处朝向语义分歧 —— (1) cxx:2417 sae.LastVertex 朝向感知 vs 端口 edge_vertices(...).1 朝向无关；(2) cxx:2423 sae.PCurve 的 CumOri 默认 true vs 端口传 false。按 .cxx 修正后实测对 F113/a3n00 中性（Shell(5) 与 STATMAP 逐字不变），门禁 t101_verify 5/5、a3n00 0.8996/T0M 0.9987/acs10 0.9846 等全部逐字不动 ⇒ 已落地。累计：已落地 §9.463(cxx:1820 簿记) 与 §9.492(cxx:2417/2423 朝向)；暂缓 §9.439(cxx:4221/4230，因面积退步)。",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.492",
    },
    {
      id: "a74",
      at: "2026-10-02",
      title: "T-101：8 轮二分锁定零长度段成因链；落地 split_by_line 簿记修复（门禁 5/5、基线逐字不动）",
      tone: "info",
      detail: "§9.452-§9.458：break_wires/collect_wires/split_by_grid 网格与调度/交点合并/沿切线建边 逐项排除或判忠实；§9.454-§9.456 定位到 split_by_grid 的 U 线 u=0(seam)/cut_index=1 切分产生两条零长度段 (±87.5,0,-32)；§9.459-§9.460 验证 WireSegment 边朝向语义(load_wires.rs:11-27)并否掉朝向假设；§9.461-§9.463 发现并落地簿记修复（cxx:1820 语义），t101_verify 5/5 通过、a3n00 0.8996 等基线逐字不变；§9.462 靶点重定向到 split_wire。",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.463",
    },
    {
      id: "a73",
      at: "2026-10-01",
      title: "T-101 止损：goal blocked（72/72），如实账目入库；留给后续 B/C/A 三条路径",
      tone: "warn",
      detail: "用户指出本轮未交付并选择止损。T-101 的 72 轮里落到 crates/ 的库改动仅 wire_fix.rs:2900 一行（门禁绿、a3n00 中性）；两个真实分歧中 ① 已自证误判，② FixDummySeam 朝向感知虽几何正确但因 a3n00 面积比 0.8996→0.8918（仅 f=171 一张面 +26 三角/总面积 -1763）被回退。过程失误（两次清理插桩多删代码、一个子任务零产出）也一并记录。终态：F113 4 wires/空面、面积比 0.8996、门禁 5/5。",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.449",
    },
    {
      id: "a72",
      at: "2026-10-01",
      title: "T-101 收尾：目标未达成；已把 F113 收敛到「切分阶段把外环切碎」，插桩清回 HEAD",
      tone: "warn",
      detail: "第 71-72 轮：F113 五阶段普查子任务无产出被停止，其留下的 4 个库文件 TEMP 插桩（perform/break_wires/split_by_grid/split_by_line，ZZCS×14）已用 git show HEAD + 归一化清回，标签残留 0、cargo check 0 error、基线 sum_mt=11941 复验一致。目标终态：F113 仍 4 wires/空面、面积比 0.8996 未上升。已确证端口能算出正确的 2-wire 面但被打包成 Shell(5)=5 块补丁(6/5/4/3/3)；打包与 .cxx 同构，剪枝 no-op，投影链全部忠实；唯一已确认的 .cxx 分歧（FixDummySeam 顶点须朝向感知）因面积退步 0.8996→0.8918（仅 f=171 一张面 +26 三角/总面积 -1763）暂缓落地。",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.448",
    },
    {
      id: "a71",
      at: "2026-10-01",
      title: "T-101 新一轮：四步计划启动（螺帽 F113 仍未解决；新候选 = check_notched_edges 的朝向预检）",
      tone: "warn",
      detail: "用户指出螺帽问题未解决。goal 已改写为四步计划并提高轮次上限（同一 goal id，因工具限制：新建需先 complete 旧目标，而旧目标未达成）。当前已确证：fix_notched_edges(wire_fix.rs:3249) 单次调用把 STEP 外环 8 条变 6 条（配对探针 8->8 x10、8->6 x2），删的是 #5012/#5018 同端点去-回毛刺边，并在中点 -65.243610 造新顶点，而 OCCT 收在原顶点 -49.864906；下游 Shell(5)->4 wires->face-checker 自交->空面。新候选见 §9.425（:2866-2877 的 p2d1/p2d2 朝向判据与 cxx:1960-1961 相反，而 :2895-2896 的 pt1/pt2 正确）。",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.426",
    },
    {
      id: "a70",
      at: "2026-10-01",
      title: "T-101 ①：根因定位到 ComposeShell 的 ClosedMode 切割顶点不同一（同一 patch 两张重合面）",
      tone: "success",
      detail: "① 的根因定位（§9.377）：两侧同标签阶段普查（OCCT `_dbg/ZZ_ComposeShell.cxx` 的 `ZZCS/ZZCW` ↔ 端口 `perform.rs` 的同名 `ZZCS`）给出第一处分歧 —— OCCT `ShapeFix_ComposeShell.cxx:2131-2275`(SplitByGrid)→`:1433-1914`(SplitByLine，ClosedMode 的 U 切) ↔ 端口 `split_by_grid.rs:93-133` + `split_by_line.rs`/`split_wire.rs`。机制（a3n00 f=140 Cone）：`loadwires` 12e/C[A→A] ⇒ `splitbygrid` 13e/C + 新外部 1e/O[B→C] ⇒ `breakwires` 切成 5e/O[A→D] 与 8e/O[C→A]，其中 D=(16.760355,2.175248,-167.075088) ≠ C=(16.704518,0,-167.130925)（差 (0.056,2.175,0.056)）且 D 别无出处 ⇒ 两半端点不配对 ⇒ `collectwires` 对 5e 走 `index=None`（cxx:2826 `!index`）⇒ 输出 [5e/O,8e/O,3e/C] ⇒ `dispatchwires` 出同一 patch [0,2π]×[-2.0207,2.0207] 的**两张重合面**（13e+3e；OCCT 孪生 1 face/1 wire/16e）⇒ FixMissingSeam 返回 Shell(2)、reader 只收 Face 而丢弃、面保持 2 wires。剪枝已实测 no-op（r1 fix_small 5→5/8→8/3→3；r2 check_small_area 三次 false）⇒ **不是剪枝判据不等价**。最小修法范围：`split_by_line.rs`/`split_wire.rs` 在 ClosedMode 插切割边时保持 wire 顶点与切割段端点同一（OCCT 复用已有顶点对，端口留下新造顶点 D）。**未验证**：OCCT 同输入必产 1 face 未直接测到（ZZ 只覆盖探针发起的调用；对已愈合 1-wire 面调 OCCT FixMissingSeam 崩 exit 1）；已证的是端口该段输出自相矛盾。子任务的 TEMP 插桩全部用 `edit` 反向撤除，`git status`/`git diff` 为空；本轮库代码零改动。",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.377",
    },
    {
      id: "a69",
      at: "2026-09-30",
      title: "T-101 迭代（1）：seam 缺口只剩 113/140/170 三个面，且合并正确、坏在「结果面 vs shell」",
      tone: "success",
      detail: "① 上一轮接回导入期 seam 步骤后，端口 wire 直方图 {1:206,2:10,4:2,6:4,10:4} 与 OCCT {1:208,2:9,4:1,6:4,10:4} 只差 2 个面。② 本轮把缺口锁到 3 个面：113（§9.303 的 F113，面积缺口主要贡献者）、140、170；用**既有**探针 `zz_probe_a3n00 --fixms` 实测：只有这 3 个面 `ret=true` 且 `result=Shell`（5/2/2 个面），其余 115 个周期面要么 `ret=false`（已在 reader 里修好）要么 Face。③ 合并本身是对的：`zz_seam_fix --all` 里“多改”的 4 个面（113/138/140/170）就是这批，且后置边数与 OCCT 孪生逐一相同（113: 22→28 = model 77 的 22+6；140: 13→16 = model 39 的 16；170: 5→8 = model 69 的 8；比边数要用 `W k edges=N` 而非 `E` 行）⇒ 分歧只在 `shape_fix_face.rs:590-673` 的尾巴（假想 grid → ComposeShell → myResult → 两轮剪枝），已交并行子任务对照 `ShapeFix_Face.cxx` + `shape_fix_compose_shell` 深挖（子任务还在跑）。④ 密度（法兰孔面 ~72 vs GT 52）先遭遇一个测量问题：新增的 `--fstats` 与 `zz_uv_feed --ids` 逐面三角数在 **89/226** 个面上不同（bbox 226/226 相同）⇒ D16 同类，跨探针逐面读数不可互比，需同一 run 内取数。⑤ 本轮库代码零改动（新增一个 example 仪器 `--fstats`），未跑全量门禁。",
      ref: "specs/_a3n00_gap_analysis.md §9.370",
    },
    {
      id: "a68",
      at: "2026-09-30",
      title: "带倒角的法兰盘已正确处理：16 面 2 wires→OCCT 的 1 wire/5 边，a3n00 面积比 0.8627→0.8996，rescue 16→0",
      tone: "success",
      detail: "① 做法：把 f1e56776 删掉的导入期 seam 步骤接回 `read_topology.rs::resolve_face` 末尾（`ShapeFixFace::fix_missing_seam` + 结果面每条 wire 的 `check_pcurves_and_shift`），对照 `ShapeProcess_OperLibrary.cxx:785-899` → `ShapeFix_Face.cxx:482-498` / 构造 `:1722-2330`。② 为何可以接回：探针 `--faceids` 开/关 ShapeProcess 得 `{1:208,2:9,4:1,6:4,10:4}` vs `{1:163,2:53,4:2,6:4,10:4}`（后者 = STEP/port）⇒ OCCT 导入确实合并了 45 个 2-bound 面，§9.219 的移除理由不成立。③ 逐面验证：带倒角法兰的 16 个面（10 圆柱 + 3 圆锥 + 3 BSpline 型）从 2 wires 变成 1 wire / 5 边，与 OCCT 模型孪生 f=85 的边数/边序/参数区间逐项相同；外置 mt 由 204–336 降到 72–110。④ 整模型：a3n00 wire 直方图 {1:163,2:53,4:2,6:4,10:4}→{1:206,2:10,4:2,6:4,10:4}；面积比 0.8627→0.8996；T0M 未网格 7→6、面积比 0.9786→0.9987；acs10 0.9025→0.9846；UNPORTED UV-grid rescue 16→0。⑤ 门禁：cargo check 0 error、--lib 1281 passed/0 failed、step_obj_gates 5/5（283.56s），未改任何断言或 area_tol。⑥ 剩余（已立 T-101）：还差 1 个 2-bound + 1 个 4-bound 面未合并；全形状跑一遍比逐面多改 4 个面；法兰孔面密度 72 vs GT 52。",
      ref: "specs/_a3n00_gap_analysis.md §9.369",
    },
    {
      id: "a67",
      at: "2026-09-30",
      title: "同结构喂 OCCT：18/18 逐字段完全相同 ⇒ Delaunay 层忠实；病灶是端口的面拓扑（2 wires vs OCCT 1 条 seam 闭合 wire）",
      tone: "success",
      detail: "① 做法：端口侧 TEMP 落盘 delaun_in_f（TOL/CELLS/N/L/V + AFTER_CTOR/FINAL）与 delaun_uv_f（登记序原始 UV）；探针新增 --delaunstruct（AddNode 逐个核对返回值等于声明序号、AddLink×72、SetTolerance/SetCellSize，再用 public 构造 BRepMesh_Delaun(aStruct, indices, cellsU, cellsV)）。② 18 个面（含 f=0、f=208 两个健康对照）的 NbNodes/NbLinks/ElementsOfDomain 与四种链状态计数与端口逐字段相同，idx_mismatch=0、links_in_struct=72 ⇒ Delaunay 层忠实。③ 更正 §9.367：只喂点时的 domain=0 是 cleanupMesh 删光（cxx:1028 / :832-835 / :908-911），与点集无关。④ 更正 §9.365.2：a3n00 两侧面序号只有 3/226 对应，bbox 配对后 port192 → model85、port169 → model87…那 16 面的 GT 是 52–80（不是 2–15），真实倍数约 4 倍且双向。⑤ 真正的输入差：STEP 声明 2 个 FACE_BOUND，OCCT 读入期合并成 1 条 wire（--faceids：默认 {1:208,2:9,4:1,6:4,10:4}，--nofix {1:163,2:53,4:2,6:4,10:4}）；OCCT 的那条 wire 由两条 2 点 seam 边把上下两圆连成一环（port192↔model85 逐边对照）。⑥ 下一轮：实现读入期的 bound 合并控制流（shape_fix_face.rs:147 的 w2 != null 合并 / :406 的 FixReorder）；可证伪预测：结构变 1 wire/83 点 ⇒ domain 由 0 变 >0、mt 228→约 52、面积比 0.8627 上升。⑦ 门禁：cargo check 0 error、--lib 1281 passed/0 failed、a3n00 stats=226 unmatched=0 mesh=11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.368",
    },
    {
      id: "a66",
      at: "2026-09-30",
      title: "把端口的点喂给 OCCT 的 BRepMesh_Delaun：同样得到 domain=0（但对照不完整）",
      tone: "warning",
      detail: "① 做了什么：端口侧加 OCCT_TOPO_DUMP_DELAUN 目录参数（在 register_node 捕获、在 perform 对那 16 面落盘）⇒ 生成 delaun_f169.txt 到 delaun_f212.txt 共 16 个文件，每文件 74 个点（u 与 v 与 movability，5 表示 Frontier；等于 2 wires 乘 37，与 §9.356 的 boundary_uv=74 一致）；探针侧加 --delaunfeed 目录参数，把点填进 IMeshData::Array1OfVertexOfDelaun 并构造 BRepMesh_Delaun。② 结果：DELAUNFEED f=169 pts=74 nodes=75 links=219 domain=0 frontier=0 fixed=0 free=0 deleted=219，16 个面全部逐字相同 ⇒ nodes=75 与 links=219 和端口的 nodes=75 与 links=221 几乎相同（差 2 条链），domain=0 ⇒ OCCT 的 BRepMesh_Delaun 用这批点也三角化不出任何域内元素，而且 frontier 与 fixed 与 free 全为 0、deleted=219，与端口 §9.356 观测到的 frontier=0 fixed=0 free=0 一致。③ 必须说明的对照边界（决定本条能推出什么）：探针这一侧没有注册约束链 —— BRepMesh_Delaun::Init(vertices)（cxx:237-250）只做 AddNode 逐个加 perform，不加任何 link；想加约束要走 ProcessConstraints()，但它不能调用（§9.366 的陷阱：header 内联、内部调 private 的 insertInternalEdges 与 frontierAdjust ⇒ 访问级别进符号名 ⇒ LNK2019）。而端口那一侧有 72 条 frontier 链。所以本轮实际比的是「OCCT 加只有顶点（0 约束）」对「端口加顶点加 72 条约束」—— 两者都得出 domain=0，但不是同一配置。能推出的是：这 74 个点本身不足以三角化出域内元素（点集不是「少给了约束就能救」的那种）；不能推出的是：OCCT 在有约束时会得出 0 —— 这一点仍未测。④ 这个否定的价值：它排除了一个很自然的假设「端口是不是漏注册了什么，导致 Delaunay 建不出三角？」—— 答案是不是漏注册那么简单：即使把约束全部拿掉（更宽松），OCCT 也是 0。⇒ 病灶不在「注册了哪些链」这个层面，而在这 74 个点本身的分布或性质（或「点加链」的组合）。这与既有线索一致（§9.364 端口与 STEP 的 wire 结构一致、§9.363 Delaunay 侧 21 处逐行等价），指向越来越集中在「这 74 个点」上。⑤ 下一轮两件事按顺序：(1) 补全对照 —— 绕开 ProcessConstraints 链接问题的办法是在构造后直接对 Result()（BRepMesh_DataStructureOfDelaun）调 AddLink（它是 public）把 72 条 frontier 链加上，再用 public 构造函数 BRepMesh_Delaun(theOldMesh, cellsU, cellsV)（BRepMesh_Delaun.hxx:39-42）建 Delaun，使约束在初始结构里、createTrianglesOnNewVertices 照常 ProcessConstraints；这需要端口 dump 同时给出链的 (first, last, state)；(2) 查这 74 个点 —— 既然同点集在 OCCT 里也建不出，就该查 OCCT 在真实管线里喂给 Delaunay 的是不是同一批点：在 OCCT 侧 dump IMeshData 的 GetFace(0) 与 GetWire(w) 与 GetEdge(e) 的 pcurve 采样点，与端口的 74 点逐点比。若 OCCT 的点不同 ⇒ 病灶在端口取点的环节（collect_boundary_uv 与 pcurve）；若相同 ⇒ 病灶在「点加链」的组合（约束链的朝向与顺序）。⑥ 本轮改动：端口插桩（node_insertion.rs，TEMP T-99）已按 §9.353 新规程用 edit 反向撤除、未用 git checkout ⇒ node_insertion.rs 的 diff 为空；探针新增 --delaunfeed（TEMP T-99）保留作下一轮脚手架。crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.367",
    },
    {
      id: "a65",
      at: "2026-09-30",
      title: "可行性已证实：探针可直接构造 BRepMesh_Delaun 并读 Result()",
      tone: "success",
      detail: "① 动机：§9.365 确认 OCCT 网格化全部 226 面（最少 2 三角）而端口在 16 面上产出 0 个域内元素 ⇒ 下一步必须测 OCCT 的 Delaunay 状态。但 BRepMesh_NodeInsertionMeshAlgo 是模板类 ⇒ 无法直接挂钩。干净的替代路：绕开 reader 与整个 IMeshData 管线，在探针里直接用 BRepMesh_Delaun 喂端口注册的同一批点。② 本轮做的可行性验证（在 specs/occt_probe/occt_probe.cpp 加 --delauncheck 分支，TEMP T-99，只跑合成输入）：用 IMeshData::Array1OfVertexOfDelaun 装正方形 4 角加 1 个内部点，构造 BRepMesh_Delaun，再读 Result() 的 NbNodes 与 NbLinks 与 ElementsOfDomain().Extent()。结果 BUILD OK、DELAUNCHECK nodes=8 links=18 domain=4 ⇒ BRepMesh_Delaun 可从本探针构造、Result() 可读、Delaunay 确实建出域内元素 ⇒ 下一轮那条路已确认可达。③ 一个新的链接约束（记录以免重复踩）：BRepMesh_Delaun::Frontier() 不能从探针调用 —— 报 LNK2019 unresolved external symbol 「private: ... getEdgesByType(...)」 referenced in 「public: ... Frontier(void) const」。Frontier() 与 InternalEdges() 与 FreeEdges() 是 header 内联、内部调 private 的 getEdgesByType ⇒ 访问级别进了符号名 ⇒ 链接必然失败（与会话早前 computeLengthU 与 computeLengthV 同类）。可用的是 Result()（返回 const 引用，而 BRepMesh_DataStructureOfDelaun 的 NbNodes 与 NbLinks 与 ElementsOfDomain 皆公开）—— 实测可用。所以下一轮要读 frontier 数得从 Result() 自己统计链状态（遍历 1 到 NbLinks 取 GetLink(i).Movability()），不要用 Frontier()。④ 下一轮做法（可执行）：(1) 端口侧加 temp dump：把 16 个失败面注册进 Delaun 的点输出成 (tag, u, v, movability) 列表（在 node_insertion.rs 的 init_data_structure 里逐个 register_node 时捕获）；(2) 探针加 --delaunfeed 文件参数：读该列表、调 AddVertices 与 ProcessConstraints，打印 NbNodes 与 NbLinks 与 ElementsOfDomain().Extent() 与链状态直方图；(3) 判定：OCCT 也得 domain=0 ⇒ 端口忠实、病灶在输入点本身（应得 2 到 15 个元素的那批点集合不对）；OCCT 得 domain 约 2 到 15 ⇒ 端口在 Delaunay 中间状态上与 OCCT 有差异，而差异不在已核对的 21 处代码里 ⇒ 只能来自数据（注册顺序、SetCellSize 与 SetTolerance 实际取值、链的初始朝向）。这一步把「测哪一层」彻底解决，且不需要改任何 Rust 逻辑代码。⑤ 本轮 Rust 源码零改动；探针源码加了 --delauncheck（TEMP T-99，约 25 行）作为下一轮的脚手架、本轮保留（默认不触发）；build.bat 得 BUILD OK，--uvsum 仍输出 226 行（未被影响）。crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.366",
    },
    {
      id: "a64",
      at: "2026-09-30",
      title: "决定性：OCCT 在这 16 个面上只铺 2-15 个三角形，端口铺 204-336（约 50 倍）",
      tone: "success",
      detail: "① 方法：重建探针（specs\occt_probe\build.bat 得 BUILD OK，用的是仓库既有源码）后跑新的一次 probe.bat 加绝对路径 data\occ\a3n00.stp 与 --uvsum 0.215，得 gt215.txt 共 226 条 UVSUM f=（与端口面数一致），三角数取 BRep_Tool::Triangulation(face)->NbTriangles()；口径用面序号，因为两侧都按 TopExp_Explorer 与 IMeshData_Model::GetFace 的顺序。② 逐面读数（同批面号）：169 是 GT 5 对 port 204、174 是 5 对 288、189 是 5 对 36、192 到 199 各是 GT 2 对 228（114 倍）、204 是 6 对 36、206 是 3 对 336（112 倍）、209 是 15 对 224、210 是 4 对 224、212 是 4 对 224。③ 三个关键读数：(a) gt faces with tris==0 为 0、gt min tris 为 2 ⇒ OCCT 在 d=0.215 下把 226 个面全部网格化，一个零三角面都没有；(b) 在这 16 个面上 OCCT 只铺 2 到 15 个三角形（192 到 199 各 2、206 是 3、169 与 174 与 189 各 5）—— 在 deflection 1.076 下这些面很小，只需 2 到 15 个，完全合理；(c) 端口在同样面上铺 204 到 336 ⇒ 约 50 倍、个别 114 倍。所以这不是补洞，而是用错误的密度铺了错误的三角形。④ 这修正了 §9.352 的推断（那里用未网格化的 v_uvsum.txt 配对，得出 GT 侧被判失败或无三角 —— 那个推断是错的，同样源于配对口径）。⑤ 与面积缺口的关系（重新表述）：端口在这 16 个面上多铺 3396 减 55 等于 3341 个三角形 ⇒ 端口的三角形预算被这 16 个面吃掉（GT 只需 55 个的地方花了 3396 个），其余面因此铺得更少 —— 这正是 §9.337「端口网格更细却面积更小」矛盾的局部解释：总量不缺，缺的是把三角形用在正确的面上。⑥ 下一轮落点不变但依据更强：病灶仍是 Delaun::add_vertices（§9.359 与 §9.360：这 16 面在那里产出 0 个域内元素，pre=145 转 post=0）。现在多了校准用的目标值：修好后这 16 个面应各自产出 2 到 15 个域内元素，而不是 0（192 到 199 各应是 2、206 应是 3，都是具体数字）。⑦ 旁支发现：UVSUM 的 wires= 不是面的 bound 数 —— gt215.txt 的 wire 直方图是 {1: 208, 2: 9, 4: 1, 6: 4, 10: 4}，而 STEP（§9.364）与端口都是 {1: 163, 2: 53, 4: 2, 6: 4, 10: 4} ⇒ 三者是三个不同的量（STEP bound 数、端口 wire 数、OCCT 离散 wire 数），§9.357 的「GT wires=1」正是拿第三个量解释第一个量，这是误读的来源。⑧ 本轮未改 Rust 源码，只重建了探针（用既有源码）；crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.365",
    },
    {
      id: "a63",
      at: "2026-09-30",
      title: "自我否定：STEP 文件本身就有 53 个面是 2 个 bound，2-vs-1 wire 不是端口缺陷",
      tone: "warning",
      detail: "① 可查证的事实（.target-gate/step_bounds.py 直接解析 data/occ/a3n00.stp 的 ADVANCED_FACE 到 bound 到 EDGE_LOOP）：ADVANCED_FACE parsed 226、bounds-per-face histogram 为 {1: 163, 2: 53, 4: 2, 6: 4, 10: 4}、FACE_*_BOUND parsed 341、EDGE_LOOP parsed 341。这与端口的 wire 数分布（§9.358 的 OCCT_TOPO_TRACE_WIRE2 实测）逐字相同，都是 {1: 163, 2: 53, 4: 2, 6: 4, 10: 4}。面序也对应（port f=0 对 step #260、port f=225 对 step #8958）。所以端口的面与边界结构与 STEP 完全一致，reader 没有多切 wire。② 这否定了 §9.363.3 的推论，也否定了 §9.357 的一半：§9.357 说端口把这 16 面建成 2 条 wire、GT 是 1 条，并据「GT wires=1」推断 OCCT 只注册一条边界环。但 STEP 文件里这些面本来就是 2 个 bound（169 对 step#7302、174 对 #7493、192 对 #8008、212 对 #8607，全部 bounds=2）。若 OCCT 与端口读同一个 STEP、且都用 ADVANCED_FACE 的 bound 数建 wire，则 OCCT 必然也是 2 条 wire，所以「GT wires=1」这个读数是错的（UVSUM 的 wires= 字段口径与面的 bound 数不是同一件事，这是对 GT 数据的第三次误读）。所以 §9.363.3 作废。③ 修正后的处境：Delaunay 侧 21 处逐行等价；端口与 STEP 的边界结构一致；端口在这 16 面产出 0 域内元素；而「OCCT 在 Delaunay 层对同一批面是否也归零」从未测过 —— 此前比的是 triangles= 与端口的 TriangulationResult，没有比过两侧 Delaunay 的域内元素数。三种可能：(1) OCCT 也归零 ⇒ 端口是忠实的，差异在更下游（OCCT 归零后走别的路径产出三角形，端口走了 rescue），要找 OCCT 的后续路径；(2) OCCT 不归零 ⇒ 端口在 Delaunay 的输入或中间状态上与 OCCT 有差异，而差异不在核过的 21 处代码里，要找的是数据（点的精确坐标、注册顺序、NodesMap 与 LinksOfDomain 内容），不是控制流；(3) OCCT 报失败面（IMeshData_Failure）⇒ 与 §9.351 的 rescue 讨论合流。④ 下一轮明确且是新的：在 OCCT 侧的 Delaunay 层加探针，逐面输出 ElementsOfDomain().Extent() 与 NbNodes() 与 NbLinks() 与 Frontier().Extent()。可行性已确认：BRepMesh_Delaun 在 TKMesh.lib 里，而探针已经链接该库（BRepMesh_IncrementalMesh 可用即证）。可挂点候选是 BRepMesh_DelaunayBaseMeshAlgo::generateMesh（DelaunayBaseMeshAlgo.cxx:43 调 getCellsCount 的那层）或 BRepMesh_BaseMeshAlgo::process。若挂点不可得，退路更干净：在端口侧把 domain_elems 与输入点集一起 dump，再在 OCCT 侧直接用 BRepMesh_Delaun 喂同一批点（其构造函数与 AddVertices 是 public），完全绕开 reader。⑤ 本轮未改源码（纯 STEP 解析与对照），crates/ 非标注改动仍只有既有 3 处；crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.364",
    },
    {
      id: "a62",
      at: "2026-09-30",
      title: "阶段性：Delaunay 侧整条三角形生成路径已逐行核对完毕、全部等价 ⇒ 分歧在喂进去的边界结构",
      tone: "warning",
      detail: "① 本轮继续核对两段，均等价：checkIntersection（OCCT cxx:1324-1372 的 UpdateBndBox 与 aPolyLen 取 Length、isSkipLastEdge 则减一、isFrontier、for aPolyIt 从 1 到 aPolyLen 内 IsOut 检查、跳过 frontier 乘 frontier、intSegSeg、不等于 NoIntersection 则返回 false，对端口 frontier.rs:280-317 的 0..poly_len 同三项）；IntSegSeg（GeomTool.cxx:342-461 四点 classifyPoint 得 aPosHash、hash[0] 或 hash[1] 小于 0 的分支、等于 1 与等于 2 的分支、IntLinLin、NoIntersection 与 Same 三分支、Cross 用 Precision::PConfusion() 做 param 区间检查，对端口 geom_tool.rs:343-409 用 PCONFUSION 同）；addTriangle 的圆绑定（cxx:1378-1398 对端口 triangulation.rs:548-595，含 bind_circle 与 remove_element）。② 累计清单：公式与控制流层已逐行核对等价的共 21 处（AdjustRange、computeLengthU 与 V、computeTolerance、computeDelta、格数三支、AdjustCellsCounts、ComputeErrFactors、collectTriangles、AddElement 与 RemoveElement、SetCellSize 与 SetTolerance、initDataStructure、finish_mesh、createTrianglesOnNewVertices、frontierAdjust、deleteTriangle、meshLeftPolygonOf、findNextPolygonLink、checkIntersection、IntSegSeg、addTriangle）；实测层也一致（delta 47/47 与 56/56、range 差小于 0.68%、vertices_nb 226/226、frontier 链数、多边形符号分布 50.1% 对 50.3%）⇒ Delaunay 侧没有未核对的段落了，而它在这 16 个面上产出 0 个域内元素、OCCT 不会 ⇒ 分歧只能来自喂进去的东西。③ 而喂进去的东西确有一处已知的结构性不一致：§9.357 已查明这 16 个面在端口侧是 2 条 wire、在 GT 侧是 1 条 wire；init_data_structure 按 wire 逐个注册边界链 ⇒ 2 条 wire 会注册两条各自闭合的边界环，OCCT 的 1 条 wire 只注册一条。推论：端口在 Delaunay 里看到的边界是两条闭合环；createTrianglesOnNewVertices 仍会按点建出 145 个三角形，但 frontierAdjust 认为它们相对两条环的朝向全部不合格而删除（§9.360 的 pre=145 转 post=0），重建也失败。这与所有实测自洽：链数 72 等于 2 乘 36（两条环）而非 36（§9.356 的 f0=72）、符号分布两组相同（两条环各自 50% 属正常）、frontierAdjust 段本身没写错。所以病灶是「为什么这个闭合回转面被建成了 2 条 wire」，属 reader 与 wire 构建侧，不在 Delaunay 里。④ 下一轮换层：(1) 回到 §9.326 与 §9.357 的 2-vs-1 wire 现象，查为什么端口把一个闭合回转面建成 2 条 wire，读 reader 侧的 wire 构建并与 OCCT 的 BRepMesh_ShapeVisitor 与 IMeshData 对照；(2) 判定标准是可查证的事实而非猜测：这 16 个面在 STEP 文件里各自应是几条 EDGE_LOOP？若 STEP 本身是 1 个 loop 则是端口多切了一刀（reader bug）；若 STEP 是 2 个 loop 则是 OCCT 把它们合并了（wire 合并缺失）；(3) 会话早期已存在的 --wireinv 探针正是为此准备的。⑤ 本轮未改源码（纯对照），delaun/ 与 node_insertion.rs 的 diff 均为空；crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.363",
    },
    {
      id: "a61",
      at: "2026-09-30",
      title: "meshLeftPolygonOf 与 findNextPolygonLink 逐行等价；多边形符号分布不分离两组",
      tone: "success",
      detail: "① 逐行对照：端口 frontier.rs:90-121 与 OCCT cxx:1070-1100 等价（isForward 分支的 polygon.push(start_edge_id) 与 push(-start_edge_id)、aStartNode 与 aPivotNode 取 FirstNode 与 LastNode、aRefLinkDir、SquareMagnitude 小于 Precision2 的退化检查、:151-174 的回退、:177-182 的 cleanupPolygon 与 meshPolygon）。findNextPolygonLink（cxx:1206-1316 对 frontier.rs:186 起）逐条等价：max_angle 初值（RealFirst 对 f64::NEG_INFINITY）、LinksConnectedTo、跳过 dead 与 skipped、isSkipLeprous 与 isLeprous 跳过、Movability 为 Free 且 ElementsConnectedTo 为空则记 dead、取另一端节点、SquareMagnitude 小于 Precision2 记 dead、非 leprous 入 leprous、Angle、isFrontier 且 abs(angle) 减 pi 的绝对值小于 Precision::Angular 时令 isCheckPointOnEdge 为 false 并 angle 取绝对值、angle 小于等于 max_angle 则跳过、isCheckEndPoints 等于 other 不等于 first、checkIntersection 传 isSkipLastEdge 为 true、aNextLinkId 等于 FirstNode 等于 pivot 时取正否则取负。常量也一致：PREC 等于 PCONFUSION 即 Precision::PConfusion()、ANGULAR 等于 1e-12 即 Precision::Angular()。所以 meshLeftPolygonOf 与 findNextPolygonLink 全段等价。② 实测否定「多边形整体反向」（新插桩 OCCT_TOPO_TRACE_POLY，记录每个闭合多边形的 len 与 pos 与 neg，pos 是带正号链数即 link.FirstNode() 等于 pivot）：f=169 FAIL 143 runs 8120 links pos=4060 neg=4060 即 50.0%；f=174 FAIL 50.2%；f=189 FAIL 50.0%；f=192 到 199 FAIL 各 50.3%；正常面 f=215 ok 52.1%、f=216 ok 50.1%。聚合：FAILING pos=58689 neg=58361 即 50.1%（16 面）；WORKING pos=106639 neg=105411 即 50.3%（36 面）⇒ 两组正负号比例几乎相同，不存在某组整体反向。抽查 f=208（正常）多数多边形 pos=23 到 36 而 neg=1，f=209（失败）则 pos=6 neg=34 与 pos=40 neg=2 混杂，但聚合后两组一样，说明这是面内多边形的正常变化，不是系统性反向。所以「多边形绕向算错」被否定。③ 本轮已排除（均有实测或逐行对照）：frontierAdjust 控制流与符号惯例（§9.361）、meshLeftPolygonOf 与 findNextPolygonLink（逐行等价）、常量 PREC 与 PREC2 与 ANGULAR（等价）、多边形符号分布（两组相同）。④ 剩下的：cleanup_polygon 与 mesh_polygon 与 mesh_elementary_polygon（三角形真正被 add_triangle 生成的那段），以及 checkIntersection（它决定选哪条邻链，从而决定多边形本身）。下一轮：(1) 本轮已新增 TRI2 插桩（polygon_meshing.rs 的 mesh_elementary_polygon，打印 poly 前三项与 edges 与 oris 与 nodes，共 626 条）与面归属面包屑 trace_set_face，把这 626 条按面归属并与 OCCT 的 meshElementaryPolygon 比对；(2) 给 checkIntersection 加桩，看失败面是否选了不同的邻链。⑤ 上述四处插桩已按 §9.353 新规程用 edit 反向撤除、未用 git checkout ⇒ delaun/ 全空 diff；crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.362",
    },
    {
      id: "a60",
      at: "2026-09-30",
      title: "逐行核对排除 frontier_adjust：九段控制流与符号惯例都与 OCCT 一致，病灶在三角形的 oris 生成",
      tone: "success",
      detail: "① 逐行核对端口 frontier.rs:15-78 与 OCCT BRepMesh_Delaun::frontierAdjust（cxx:944-1026），九段全部等价：主循环 for aPass 1..2（cxx:955 对 :15）、取 frontier 集（cxx:946 对 :10）、遍历挂接三角形（cxx:961-966 对 :16-19）、跳过已删 if aPriorElemId 小于 0 continue（cxx:970-973 对 :21-23）、外三角形判定（cxx:980-988 的 if aFrontierId 等于 e[n] 且非 o[n] 就 deleteTriangle 并 break，对 :26-38 的 if element.link_at(n) 小于 0 就 delete_triangle 并 break）、找到即跳出元素循环（cxx:991-994 对 :39-41）、清理悬挂链（cxx:999-1007 对 :45-49）、重建 meshLeftPolygonOf（cxx:1011-1025 对 :51-63）、cleanupMesh 与失败前沿重试（对 :66 与 :68-77）。未发现少判或多判朝向条件。② 符号惯例也已证实一致（关键）：OCCT 把绝对号与朝向布尔分开存（BRepMesh_Triangle 的 myEdges[3] 加 myOrientations[3]）。写入时 OCCT cxx:541-548 是 anEdgeIds[i] 等于 std::abs(anEdgeInfo)、anEdgesOri[i] 等于 anEdgeInfo 大于 0；端口 frontier.rs:321-328 的 add_triangle_by_info 是 edges[i] 等于 edges_info[i].abs()、oris[i] 等于 edges_info[i] 大于 0，逐字相同。存储时 OCCT 原样存 (theEdgesId, theEdgesOri)，端口 triangulation.rs:555-559 压成带符号数组（oris 为真取正、为假取负）⇒ 正号等于 oriented（即 OCCT 的 o[i] 为 true）、负号等于反向。读取时 OCCT 判 not o[n]、端口判 link_at(n) 小于 0 ⇒ 语义一一对应。所以 §9.360 的「145 个三角形全挂在负向 frontier 链上」不是符号惯例错误，而是三角形相对 frontier 链的朝向确实为负。③ 剩下的唯一去处：frontierAdjust 与符号惯例都排除 ⇒ 产生那些三角形的 oris 本身就错。三角形产生点只有三处：frontier.rs:328（add_triangle_by_info，惯例已验证）、polygon_meshing.rs:235、triangulation.rs:427。④ 下一轮：对 f=209（失败）与 f=208（正常），在这三个调用点打印 (link_id, ori) 与多边形遍历方向（polygon 的带符号链号序列），看失败面的 oris 是否整体为反；这直接指向 mesh_left_polygon_of 的多边形走向判定（frontier.rs:90-121 的 is_forward 分支与 ref_link_dir），与 OCCT meshLeftPolygonOf 对照仍是正规移植工作。⑤ 本轮未改源码（纯对照），crates/occt-topo/src/meshing/delaun/ 的 diff 为空；crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.361",
    },
    {
      id: "a59",
      at: "2026-09-30",
      title: "精确锁定：失败面的 145 个域内元素被 process_constraints 全部删掉（pre=145 转 post=0）",
      tone: "success",
      detail: "① 在 create_triangles_on_new_vertices 里插三处计数（进入时、process_constraints 前、之后），开关 OCCT_TOPO_TRACE_CTNV ⇒ 243 条 CTNV enter 记录里「一开始就 domain=0」有 0 条、「提前 failed 返回」有 0 条、「结束时 domain=0」有 16 条。这 16 条字段完全一致：verts=72 d0=3 f0=72 pre=145 post=0 failed=False。② 逐段解读：进入时 d0=3（构造时已建 3 个）、f0=72（72 条 frontier 链，等于 2 wires 乘 36）；遍历 verts=72 个节点；pre=145 ⇒ 建出了 142 个新元素、共 145 个；post=0 ⇒ 全部被删；failed=false ⇒ 不是 add_triangle 溢出。所以三角化本身成功（145 个元素真实存在），是 process_constraints() 把它们全删了。③ 删除点是 frontier.rs:15-43 的 reorient 段：对每条 frontier 链，找挂在它上面、且 element.link_at(n) 小于 0（负向）的三角形并 delete_triangle；其后 mesh_left_polygon_of(frontier_id, true, ...) 负责按正确朝向重建。这 16 个面 145 个三角形全部挂在负向 frontier 链上 ⇒ 初始三角化的定向整体相反 ⇒ 全删；而重建失败（否则 post 不会是 0）。④ 这解释了此前所有观测：domain_elems=0 出现在 add_vertices 后（因为 process_constraints 就在 add_vertices 内被调，triangulation.rs:539）；frontier=0 与 fixed=0 与 free=0（删三角形时连带 remove_link 把链标 Deleted）；nodes=75 与 links=221 与正常面相同（数量没错，错的是链的朝向）；f=209 与 f=208 点集相同、顺序不同 ⇒ 点的顺序决定了两条 wire 的绕向，进而决定 frontier 链的朝向，这正是 §9.359「顺序决定」的机理。⑤ 下一轮已到具体代码段：(1) 逐行核对端口 frontier.rs:15-60 的判定条件 element.link_at(n) 小于 0 与 OCCT BRepMesh_MeshTool::Frontier（BRepMesh_MeshTool.cxx:284-300）及 frontierAdjust 的 reorient 段，看是否少判或多判了一个朝向条件；(2) 在一正常 2-wire 面（f=208）与一失败面（f=209）上，于 reorient 段打印被删三角形数与每条 frontier 链的 link_at(n) 符号分布 —— 正常面应删得很少、失败面删光，由此定位符号判断的分歧点；(3) 重点怀疑 link_at(n) 小于 0 的符号约定，与 add_link_to_mesh 在 Orientation::Reversed 时用 add_link(last, first, ...)（node_insertion.rs:827-835）同 OCCT addLinkToMesh 的对应分支是否一致。仍有 .cxx 原文可对照，属正规移植问题。⑥ 本轮插桩（delaun/triangulation.rs 三处）已按 §9.353 新规程用 edit 反向撤除、未用 git checkout ⇒ delaun/triangulation.rs 的 diff 为空；crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.360",
    },
    {
      id: "a58",
      at: "2026-09-30",
      title: "病灶缩到一步：Delaun::add_vertices；失败与成功由点的顺序决定，不由几何或结构决定",
      tone: "success",
      detail: "① 差集分析：a3n00 的 wire 数分布是 {1:163, 2:53, 4:2, 6:4, 10:4}；53 个 2-wire 面里 16 失败、37 正常。37 个正常面里有 34 个的逐 wire 结构与 16 个失败面逐字段完全相同（都是 (1,37,False,False,('false',)) 两次、bu=74）。分离器检验 all_edges_1 与 all_pts_eq 与 all_rw_false 都是「失败 16/16、但正常也 34/37」⇒ 不分离；any_self 与 any_open 两边都是 0 ⇒ 排除。所以没有任何结构谓词能分离这两组，「2 wires」只是必要条件，连结构完全相同的 34 个面都成功 ⇒ 结构层面已无法解释；同时排除 §9.356 的「固定采样步数」（正常面也是 pts=37 与 bu=74）。② 决定性读数（OCCT_TOPO_TRACE_TRI）：f=209（失败）是 TRI pre_useedge nodes=75 links=221 frontier=0 domain=0；f=208（正常）是 TRI pre_useedge nodes=75 links=221 frontier=3 domain=1；而且 use_edge 完全不改变 domain ⇒ 域内元素是在 add_vertices 里产生的，病灶在 Delaun::add_vertices 内部。③ 最强证据：f=209（失败）与 f=208（正常）的边界 74 个点作为集合完全相同（identical point set 为 True；u 逐位置 74/74 相同、v 只有 40/74 同位置，但 u 与 v 的多重集都相同）⇒ 两者是「同一组 74 个点、不同的遍历起点或顺序」。配合 ② 的 nodes 与 links 完全相同、只有 frontier 不同 ⇒ 在完全相同的输入几何下 add_vertices 给出不同结果 ⇒ 这不是几何脆弱性，是顺序决定。④ 下一轮已收敛到具体函数：(1) 打开 crates/occt-topo/src/meshing/delaun/ 的 add_vertices 实现，与 OCCT BRepMesh_Delaun.cxx 的对应构造函数逐行对照；(2) 在 add_vertices 内对 f=209 与 f=208 打印每一步之后的 elements_of_domain().len() 或等价中间量，找出第一步出现分歧的位置；(3) 重点核对 OCCT 的 ComparatorOfVertexOfDelaun 与 ComparatorOfIndexedVertexOfDelaun（BRepMesh_Delaun.cxx:49-75），看端口是否用了同一个排序谓词。⑤ 本轮两处插桩已按 §9.353 新规程用 edit 反向撤除、未用 git checkout ⇒ node_insertion.rs 的 diff 为空；crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.358 §9.359",
    },
    {
      id: "a57",
      at: "2026-09-30",
      title: "更正 §9.352：这 16 个失败面不是平面，是 Cylinder（12）+ Cone（4）；真特征是闭合回转面被拆成 2 条 wire",
      tone: "warning",
      detail: "① 更正：§9.352 写的「13 个配上对的面全是 Plane 且 wires=1」是错的。本轮用 UVSUM 的 type= 重新核对（新脚本 .target-gate/zero_face_types.py）得：169 对 gt87 是 Cylinder、174 对 gt42 是 Cone、189 对 gt41 是 Cone、192 到 199 分别对 gt85/84/83/78/79/80/81/82 全是 Cylinder、204 对 gt43 是 Cone、206 对 gt86 是 Cylinder、209 与 210 与 212 未配对 ⇒ 12 个 Cylinder 加 4 个 Cone，没有一个是 Plane。错因是 §9.352 读 UVSUM 时把配对面的 GT 面号与端口面号对错了行；本轮改用几何 bbox 重新配对并直接读 type=。这是对同一批 UVSUM 数据的第二次读错，§9.355 的结论也建立在同一批错读上。② 更正后的真规律：面型是 Cylinder 或 Cone（闭合回转面）；端口 wires=2、GT wires=1；每条 wire 各 1 条边（IDENT_WIRE f=169 的 edges=[414] 与 edges=[415]）；每条 pcurve 37 点 ⇒ boundary_uv=74、nodes=75、frontier=0。所以真规律是「端口把闭合的圆柱与圆锥面建成了 2 条 wire（GT 是 1 条），这个 2-wire 结构让 Delaunay 注册不出任何域内元素」，与 §9.352 说的「单闭合 wire 的平面」正好相反。③ §9.355 需按本节收窄：它说「2-vs-1 wire 不是因果」，依据是「28 个里只有 13 个失败」，但本轮查明那 28 个面与 16 个失败面是同一类（Cylinder 与 Cone 的 2-wire）⇒ 正确表述是「失败面集合是 2-wire 回转面集合的子集，但不是所有 2-wire 回转面都失败（28 个里 13 个配上对的失败、另 15 个正常）」⇒ 2 wires 是必要条件而非充分条件，失败的额外条件尚未找出。④ 本轮的两个读数（OCCT_TOPO_TRACE_IDENT）：16 个面的 face_key 各不相同（1917249376176 与 1917250123632 与 1917250242560 等）、边集各异、boundary_uv 首末坐标完全不同 ⇒ 它们不是同一份数据，候选「共享结构」被否定；而「固定采样步数」也不是原因（37 点来自每条 wire 一条闭合边的自然采样）。⑤ 下一轮做法：(1) 先修 §9.352 与 §9.355 的表述，避免后续按「平面」去找；(2) 在那 15 个 2-wire 但正常的回转面与这 16 个失败的之间做差集分析 —— 打印 wire[0] 与 wire[1] 的 pcurve 点数、两条 wire 的 uv 起止、以及 wire.is_status(SELF_INTERSECTING_WIRE 或 OPEN_WIRE)，找出「失败」相对「正常」多出来的那一项；(3) 参考 OCCT 的 CylinderRangeSplitter 与 ConeRangeSplitter 以及 NodeInsertionMeshAlgo::initDataStructure 对闭合回转面的处理（有 .cxx 可对照）。⑥ 本轮插桩已按 §9.353 新规程用 edit 反向撤除、未用 git checkout ⇒ node_insertion.rs 的 diff 为空；crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.357",
    },
    {
      id: "a56",
      at: "2026-09-30",
      title: "失败链再上一环：三角化 0 个域内元素导致全部边界链被标 Deleted；16 个失败面的 boundary_uv 全是 74",
      tone: "success",
      detail: "① 本轮插桩量出逐环：失败面是 2 wires、每条 wire 1 个 slot（单条闭合边），正常面 f=0 是 1 wire 与 8 个 slot。失败面确实注册了链：WIRES_EDGE f=169 wire_it=0 pts=37 nodes_now=36 links_now=36、wire_it=1 pts=37 nodes_now=72 links_now=72（共 72 节点与 72 链），且 DELAUN_SKIP 总数为 0（排除静默跳过边这个候选）。但到 finish_mesh 时 frontier=0 fixed=0 free=0（三种活状态全 0），而链本身还在（links=221）⇒ 所有链都被标成 Deleted；正常面是 frontier=68 free=68。② 完整链条共 6 环、全部有实测：2 wires 乘 1 slot、每条 pcurve 37 点 → 注册 72 节点与 72 链 → Delaunay 三角化 domain_elems=0 → 清理（无域内元素则链无引用）全部链变 Deleted → frontier=0 → use_edge 循环空转 → collect_triangles 遍历空集 → triangles=0 → map_triangulation 按 :569-573 拒绝 → needs_fallback=true → rescue 硬铺 224 到 336 个三角形。所以症状 rescue 距病灶有 6 环，真正的病灶是 Delaunay 三角化对这 16 个面产出 0 个域内元素。③ 新增强信号：boundary_uv 在 16 个失败面上全部等于 74（正常面范围 8 到 370）、nodes 全部等于 75、frontier 全部等于 0（正常面范围 3 到 360）；且本轮量出 74 等于 2 乘 37，每条 wire 恰好 37 个点。16 个几何完全不同的面（f=169 与 f=206 的 bbox 差 60 多单位）却给出完全相同的 boundary_uv=74 与 nodes=75，不可能是几何决定的自然结果。两个候选解释待下一轮判定：(1) 这 16 个面在端口的离散模型里指向同一批或同构的数据结构；(2) 37 是某个固定采样步数，这 16 个面恰好都走了那条固定采样路径。④ 下一轮做法明确：在 init_data_structure 里对这 16 个面与 2 个正常面打印 model.face(face_index) 的 key 或指针、每条 wire 的 edge index 与 edge_reverse_walk 与 edge.pcurve(pc).points().len()、以及 self.boundary_uv 的前 3 与后 3 个 uv 值。若同一批坐标成立则病灶在离散模型或 reader 侧、与 Delaunay 无关；若坐标各异而只有点数 37 相同则病灶在采样步数。⑤ 本轮三处插桩已按 §9.353 新规程用 edit 反向撤除、未用 git checkout ⇒ node_insertion.rs 的 diff 为空；全局 crates/ 共 64 改动 = 61 个 A 类标注 + 3 个既有功能改动。门禁：cargo check 0 error、--lib 1255/26、a3n00 stats=226 unmatched=0 mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.356",
    },
    {
      id: "a55",
      at: "2026-09-30",
      title: "失败点定位到 Delaunay 的 0 个域内元素；2-vs-1 wire 被否定为因果；新签名是 nodes=75 且 domain_elems=0",
      tone: "success",
      detail: "① 在 node_insertion.rs 的 finish_mesh 里、collect_triangles 之前打印 delaun.result().nb_nodes() 与 elements_of_domain().len()：226 条记录里恰好 16 条 domain_elems=0，出现序号与端口面号逐一对应 [169,174,189,192,193,194,195,196,197,198,199,204,206,209,210,212]，与 §9.351 的 rescue 触发面完全一致 ⇒ 触发链的起点是「Delaunay 三角化本身产出 0 个域内元素」，不是后来被丢。链上每环均有实测：domain_elems=0 转 triangles=0 转 map_triangulation 按 :569-573 拒绝 转 needs_fallback=true 转 rescue。② 求交否定 2-vs-1 wire 为因果：167 个已配对面里「port 2 wires vs GT 1 wire」有 28 个，但只有 13 个 domain_elems=0；另 15 个 Delaunay 正常（41,42,44,45,47,48,92,93,129,130,137,140,154,170,173）。反方向那 3 个（209/210/212）并非例外 —— 查 PORTID 它们同样 wires=2，只是没在 GT 里配上面。所以 16 个失败面在端口侧全部 wires=2，但端口 wires=2 的面至少 31 个（28 已配 + 3 未配），只有 16 个失败 ⇒ wire 数相关但非因果。③ 真正的签名：这 16 面全部恰好 nodes=75 且 domain_elems=0，而正常面节点数各不相同（抽样 occ=0 是 nodes=71 与 domain_elems=66、occ=2 是 nodes=45 与 40）⇒ 75 是这 16 面共有的异常值，不是几何决定的自然结果，比 wire 数更特异。④ 下一轮唯一该追的量：在这 16 面与一个正常面上打印 boundary_uv.len()、wires_uv 每条 wire 的点数、以及注册过程中 structure.nb_nodes() 的逐步值，找出 75 怎么来的、哪一步之后 elements_of_domain() 一直为空。不要先动 rescue（它是症状）；也不要再追狭长面（§9.350 撤回）或 2-vs-1 wire（本轮否定）。⑤ 本轮插桩已按 §9.353 的新规程用 edit 反向撤除（未用 git checkout）并修回被误删的一个空行 ⇒ node_insertion.rs 的 diff 为空（该文件不在 A 类清单、原无标注）；全局复核 crates/ 共 64 个改动 = 61 个 A 类标注 + 3 个既有功能改动（gcpnts.rs、ginter.rs、make_pcurve.rs）。门禁：cargo check 0 error、--lib 1255/26（环境性）、a3n00 stats=226、unmatched=0、mesh 11101/11121。",
      ref: "specs/_a3n00_gap_analysis.md §9.354 §9.355",
    },
    {
      id: "a54",
      at: "2026-09-30",
      title: "根因收敛到一点：16 个 rescue 面全是 Plane + wires=1，端口 Delaunay 对它们产出 0 三角形",
      tone: "warning",
      detail: "本轮把 rescue 的成因追到底。① 触发链（读源码）：discret_root.rs:442 调 map_triangulation，而它在 :569-573 里 if triangles.is_empty() 就 return Err（Delaunay produced no triangles），于是 :443-444 的 tri 为 None、needs_fallback 为 true，:484 的 or_else 走 wireframe_face_triangulation。关键：这条链绕过了 Ok(Err(e)) 分支里那段「face 已标记 FAILURE 就跳过、不 fallback」的忠实处理（:447-459），因为这里 perform 是 Ok、只是结果为空。② 16 个 rescue 面按几何 bbox 与 UVSUM 配对（用未网格化的 v_uvsum.txt 与 PORTID 同口径，三角数取同面号 v_uvsum2.txt），13 个配上对的面全部是 Plane 且 wires=1（169 对 87、174 对 42、189 对 41、192 对 85、193 对 84、194 对 83、195 对 78、196 对 79、197 对 80、198 对 81、199 对 82、204 对 43、206 对 86）⇒ 不是随机分布，是一个面类的共性：端口的 Delaunay 节点插入对单闭合 wire 的平面产出 0 个三角形。③ 量与算术：这 16 面 port 提交 3396 对 OCCT 同面 1049（端口多铺 2347）；port 去掉 rescue 后 = 11121 减 3396 = 7725，而 OCCT 在 d=0.215 是 12462 ⇒ rescue 既不是忠实补洞（比 OCCT 多 2347），也不是可有可无（拿掉这些面就从多铺变完全不铺）。④ 两种可能，证据偏向 A：A 是端口 Delaunay 在这类平面上有缺口、该修的是 Delaunay，rescue 只是掩盖症状；B 是 rescue 当年为对齐某基线刻意保留。偏 A 的理由：这 16 面是同一个面类，符合某算法分支没覆盖的形态而非人为挑选；且代码注释自己承认 OCCT 没有 per-face 备选 tessellator，并说这条路径 decided to keep this code for now (no triggering case)，也就是留它是因为以为它永远不会跑，而实测它每次都跑。⑤ 建议下一步有序：(1) 先修那段错误注释，并把 §9.306 的 no triggering case 一并更正；(2) 再查端口 Delaunay 为何对 Plane 加 wires=1 产出 0 三角形 —— 这是有 .cxx 可对照的正规移植问题（BRepMesh_Delaun 的初始化与约束边注册），与自造分支无关，是真正该修的；(3) 只有在 (2) 证明端口无法在 OCCT 路径上产出时，才讨论 rescue 的去留。⑥ 另附事故与规程（§9.353）：本轮撤插桩时误用 git checkout discret_root.rs，把该文件的未提交注释改动一并还原到 HEAD。已定界：mt=11121、stats=226、unmatched=0 与基线逐字相同 ⇒ 丢失的仅注释、无可执行逻辑。已按同类文件逐字恢复 T-97 标注使该文件 diff 回到 6/0；全局复核 61 个 6/0 改动等于 61 个带标注文件，非标注改动只剩既有 3 处。新规程：撤除临时插桩不得用 git checkout，改用带 TEMP T-99 标记的逐块 edit 反向替换（本会话已第四次因 git checkout 丢失未提交改动，第一次在 §9.324 造成数十轮浪费）。两处插桩已撤除、无残留，门禁未变。",
      ref: "specs/_a3n00_gap_analysis.md §9.352 §9.353",
    },
    {
      id: "a53",
      at: "2026-09-30",
      title: "本会话最重要发现：rescue 路径在 a3n00 上触发 16 次（注释称从不触发），供给 30.5% 的三角形",
      tone: "warning",
      detail: "本轮按卡片转去量节点插入/三角化，意外定位到一条自造分支的实际作用。① 测量一：给 discret_root.rs:427 的 Ok(Ok(result)) 加 env-gated 打印 OCCT_TOPO_TRACE_TRIRES，逐面输出 TriangulationResult.triangles.len() 与 nodes.len()，与同一次运行的 PORTID mt= 比较：trires faces=226、committed faces=226、210 个面相等、16 个面不同，sum trires=7725 对 sum committed=11121，差 3396。16 个差异面形态完全一致：TriangulationResult 产出 0 个三角形而提交网格有 224 到 336 个（f=206 0 对 336、f=174 0 对 288、f=192 至 f=199 各 0 对 228、f=209/210/212 各 0 对 224、f=169 0 对 204、f=189/204 各 0 对 36）。排除其他解释：map_triangulation 在 discret_root.rs:533 是逐顶点 1:1 重映射、不增删三角形，所以这 3396 个三角形必来自另一条路径。② 测量二：在 discret_root.rs:484-503 的 needs_fallback 分支加 OCCT_TOPO_TRACE_FALLBACK 打印后跑 a3n00，得到 16 行 FALLBACK reached，且触发面集合与测量一的 16 个零产出面完全一致 ⇒ 提交网格那 3396 个三角形全部来自 wireframe_face_triangulation（rescue）。③ 该段代码自己的注释断言 this rescue is never reached on any gate model、reached_rescue=0 everywhere，实测触发 16 次，注释是错的；注释同时写道 OCCT has no per-face alternative tessellator、失败的面保持 IMeshData_Failure（BaseMeshAlgo.cxx:52-62），即这条路径正是 CLAUDE.md 门禁第 2 条禁止的自造 OCCT 里不存在的分支，而它并不惰性：为 a3n00 供给了 3396/11121 即 30.5% 的三角形。它还同时解释了 §9.340 的过密榜（f=206 gt 68 对 port 336、f=192-199 各 gt 52 对 port 228）：这些面在 OCCT 里被判失败或无三角，GT 数很小，而端口用 rescue 硬铺了 224 到 336 个。④ 处置：删掉或停用这条路径是实质行为改动（会牵动 a3n00 网格与面积比，且可能是当年为对齐某基线刻意保留），按门禁第 4 条我没有动手，把它作为待人工决定的事项写进卡片 next 与 nextAction；下一轮若动手，第一件事是改掉那段错误注释并把 §9.306 的 no triggering case 更正为 a3n00 上触发 16 次。⑤ 本轮在 discret_root.rs 留了两处 env-gated 临时插桩（TRIRES 与 FALLBACK，默认不输出、零逻辑改变），下一轮用完须撤除。cargo check 0 error，--lib 1255/26（环境性），a3n00 基线未变。",
      ref: "specs/_a3n00_gap_analysis.md §9.351",
    },
    {
      id: "a52",
      at: "2026-09-30",
      title: "节点插入/三角化三处逐行等价；撤回「狭长面」线索；下一步收敛到两个整数",
      tone: "success",
      detail: "本轮转向卡片定的唯一落点。① collectTriangles 逐行等价：OCCT BaseMeshAlgo.cxx:245-278 与端口 node_insertion.rs:414-460 都是每个域内元素输出一个三角形、用 used 映射做紧凑重编号、输出节点数取映射最大值。② 域标记逐行等价：OCCT DataStructureOfDelaun.cxx:170-183 的 AddElement 无条件加入 myElementsOfDomain，:197 的 RemoveElement 移除；端口 delaun_data.rs:433-442 的 add_element 同样无条件插入 elements_of_domain，:477 的 remove_element 移除。③ SetCellSize 与 SetTolerance 逐行等价：OCCT NodeInsertionMeshAlgo.hxx:86-93 与端口 node_insertion.rs:503-508 都是 14.0 乘 tolUV 再除以 delta。④ 撤回狭长面线索：用 UVSUM 的 lenU 与 lenV 算长宽比，与三角数比值并排（47 个可用配对），长宽比大于等于 3 的面全部比值 1.000 且只有 2 到 3 个三角形、两侧完全相同；差异集中在长宽比约等于 1 的面（0.5 到 1.28）；两个桶均值几乎相同（1.107 对 1.000）⇒ 差异不随长宽比缩放。所以 §9.342.3 的异常聚集在狭长面不再成立，该线索撤回，也不得用作启发式（那是 CLAUDE.md 禁止的长度滤边式规则）。⑤ 剩余唯一去向是 Delaunay 三角化本身产出的域内元素数不同，它由 BRepMesh_Delaun 的算法过程决定（插入顺序、ClassifyTriangle、边翻转、erase_free_links），不是能逐行对照就定论的公式。下一步改用比较算法状态：每面的 elements_of_domain().len() 与 nb_nodes 两个整数。该数据已在手上 —— OCCT 侧是 UVSUM 的 triangles= 字段，端口侧是 TriangulationResult 的 triangles.len() 与 nodes.len()（node_insertion.rs:90-97，均为 pub）⇒ 不必新造仪器，只需在 zz_uv_feed 打出这两个数与 UVSUM 按 bbox 配对比较。⑥ 本轮未提交任何源码改动，node_insertion.rs diff 为空（插桩已在 §9.349 撤除），门禁未变。",
      ref: "specs/_a3n00_gap_analysis.md §9.350",
    },
    {
      id: "a51",
      at: "2026-09-30",
      title: "定论：vertices_nb 两侧完全相同（226/226，无 Deleted 节点）⇒ 范围分割器整条输入路径排除",
      tone: "success",
      detail: "本轮把上一轮提出的最后一处可疑点量掉了：vertices_nb 两侧完全相同，范围分割器整条路径就此排除。① 判定测量：在 node_insertion.rs 的 cells_count(indices.len()) 处加 env-gated 打印，开关为 OCCT_TOPO_TRACE_CELLS，同时输出 nb_nodes、live_indices、boundary_uv 与 cells。结果 226 个面里 indices.len() 等于 nb_nodes 全部成立、差异 0 个，说明没有 Deleted 状态的节点，端口传给 getCellsCount 的数与 OCCT 调用点传的 aStructure->NbNodes() 完全相同。② 这同时撤清了 §9.347.3 与 §9.348.3 两处疑虑。前者说端口数未删除节点、OCCT 数全部，属口径差异，实测差异为 0；后者依据我插桩里的 vnb_live 等于 boundary_uv.len() 加 structure.nb_nodes()，推测可能重复计数把格数放大，但那只是我诊断行的写法、不是真实传参，真实传参 indices.len() 等于 nb_nodes()，不存在重复计数，该线索作废。这是一次自我更正：两个可疑点都是我提出的，本轮用一次直接测量同时排除。③ 阶段性定论：范围分割器整条路径的所有量都已核对完毕且与 OCCT 一致。公式 7 处逐行等价，含 AdjustRange、computeLengthU/V、computeTolerance、computeDelta、格数三支、AdjustCellsCounts、ComputeErrFactors；输入方面 delta 在 47/47 与 56/56 两次配对下都完全相同、最坏 2e-15，range_u 与 range_v 最多差 0.68% 且差异是周期区间起点而非跨度，tolerance 最坏 3.3e-2，vertices_nb 226/226 完全相同；基类 GenerateSurfaceNodes 是空实现。所以端口在某些面铺得远多于或远少于 OCCT 不可能由范围分割器产生，剩余唯一去向是节点插入与三角化。④ 下一轮落点唯一且明确：node_insertion.rs 的 insert_nodes 与 finish_mesh 对 OCCT 的 NodeInsertionMeshAlgo::insertNodes 与 BaseMeshAlgo::generateMesh；meshing/delaun 目录的 Delaun::new_with_data_cells、add_triangle、erase_free_links 对 BRepMesh_Delaun；collect_triangles 对 BaseMeshAlgo::collectTriangles 即 cxx 245 到 277 行。做法是在两侧同一管线点打印插入前后的节点数与 collect_triangles 得到的三角形数，对 §9.340 亏损榜的 f=206（gt 68 对 port 336）与盈利榜的 f=171（gt 1569 对 port 349）逐值对照。⑤ 本轮动作：端口的两处临时插桩即 SPLITTER 与 CELLS 打印已完成使命，已在本轮结束时撤除，crates/occt-topo/src 回到只有本会话既有 4 处功能改动的状态，node_insertion.rs 的 diff 为空；OCCT 探针侧的 SPLITTER 与 cells 输出保留，因为探针本就是仪器集合、且下一轮还要用它做同口径输出。撤除前后 --lib 仍是 1255 passed 与 26 failed（环境性），a3n00 基线 stats=226、unmatched=0、mesh 11101/11121 均未变。",
      ref: "specs/_a3n00_gap_analysis.md §9.347 §9.348 §9.349",
    },
    {
      id: "a50",
      at: "2026-09-30",
      title: "九量逐值对照完成：delta 47/47 完全相同、范围最多差 0.68% ⇒ 分歧不在范围分割器，在下游",
      tone: "success",
      detail: "本轮把卡片定的九量对照真正做出来了：结论是范围分割器的输入与公式两侧一致，分歧不在范围分割器而下在下游。① 两侧仪器均建好。端口侧在 meshing/node_insertion.rs 的 list_surface_nodes 两个分支都加了 env-gated eprintln，开关为 OCCT_TOPO_TRACE_SPLITTER，取值点在 adjust_range 之后、generate_surface_nodes 之前，打印 range_u、range_v、len_u、len_v、delta、tol、vnb、defl 与 bbox。OCCT 侧在 --uvsum 循环里按同一管线点重建 splitter，即 Reset(aDF, params)、逐 pcurve 采样 AddPoint、AdjustRange，然后打印同样的量加自身 bbox。两侧都带 bbox，因此按 bbox 六坐标即 D21 配对，因为两侧的面枚举顺序不同：端口走 faces_of，OCCT 的 --uvsum 走 IMeshData_Model::GetFace，实测索引对齐仅 3/226。一个实现细节：OCCT 的 computeLengthU 与 computeLengthV 在发行版头文件里是 private，其符号也按 private 修饰，无法从外部调用，试过 define private public 与派生类暴露两种办法，都因符号名按访问级别修饰而链接失败；所幸 --uvsum 本来就打印 lenU 与 lenV，故该项不重复输出。② 对照结果，47 个精确配对：delta 47/47 完全相同，最坏相对差 2e-15；range_v 45/47，最坏 0.081%；range_u 41/47，最坏 0.68%；tolerance 39/47，最坏 0.87%。最坏两例都是周期性闭合面上的同一条周期，例如 port f=152 的 ru=[0.4058,16.5942] 对 gt f=137 的 ru=[0.2924,16.7076]，两者跨度都恰好 16.1884、相等，只是起点相差 0.1134，来自周期域上取范围的起点约定；rv 那一例也是以 0 为中心的同一区间。所以不是范围算错，而是同一周期区间的起点与对中度不同。③ 关键否定：至此范围分割器整条路径的公式都已逐行等价核对完毕，包括 AdjustRange、computeLengthU/V、computeTolerance、computeDelta、格数三支、AdjustCellsCounts、ComputeErrFactors，见 §9.343 与 §9.344；输入也已核对，delta 完全一致、范围最多差 0.68% 且差异不是跨度；基类 GenerateSurfaceNodes 又是空实现。因此在同一 delta 下格数也应一致，端口某些面过密、某些面过疏既不来自范围分割器也不来自其输入，只能来自下游的节点插入与三角化，或 vertices_nb 这一项与 OCCT 不同。④ 下一轮落点：meshing/delaun 目录、mesh_algo.rs 的 finish_mesh 与 collect_triangles、node_insertion.rs 的 insert_nodes，对照 OCCT 的 BRepMesh_DelaunayBaseMeshAlgo、NodeInsertionMeshAlgo::insertNodes、BaseMeshAlgo::collectTriangles；做法同本轮，在两侧同一管线点打印同一批量，如 structure.nb_nodes()、插入前后的节点数、collect_triangles 得到的三角形数。同时要核清 vertices_nb 的口径：端口用 boundary_uv.len() 加 structure.nb_nodes()，OCCT 用 NodeInsertionMeshAlgo 里结构的节点数，这两者是否同一个量需要确认，因为它进 AdjustCellsCounts，会直接改变格数。⑤ 注意：本轮的端口插桩是临时的、当前保留，因为下一轮要用同一仪器，用完后应撤除；crates/occt-topo/src 因此有未提交改动，仅插桩、逻辑零改变，--lib 仍是 1255 passed 与 26 failed（环境性），a3n00 基线 stats=226、unmatched=0、mesh 11101/11121 未变。",
      ref: "specs/_a3n00_gap_analysis.md §9.346",
    },
    {
      id: "a49",
      at: "2026-09-30",
      title: "格数/误差因子逐行等价（第 16 条排除）；基类 GenerateSurfaceNodes 是空实现，落点收敛到三个输入量",
      tone: "success",
      detail: "本轮继续逐式对照，把格数与误差因子也逐行核完，并发现一个结构性事实，落点因此收敛到三个输入量。① computeTolerance（BRepMesh_DefaultRangeSplitter.cxx:111-127）与端口 splitter.rs:152-174 逐行等价，含 UResolution 乘 1.1 与 max(min(1e-5,res), 1e-7*diff)；computeDelta（.cxx:131-138）与 splitter.rs:177-182 同式。② 格数：OCCT BRepMesh_GeomTool.cxx:485-508 有 Torus、Cylinder 特例、else 三支，与端口 node_insertion.rs:935-946 逐字相同，都是 2 的 ceil(log10(...)) 次方，其中 Cylinder 的 U 除 deltaU 再除 dV；AdjustCellsCounts（GeomTool.cxx:86-144）与 node_insertion.rs:947-986 逐支相同，含末尾 max(2)；ComputeErrFactors（GeomTool.cxx:32-83）与 node_insertion.rs:990-1033 逐支相同，含 Extrusion 与 Revolution 的 errV 除以 Degree 乘 NbKnots、Bezier 按 U 与 V 的 degree、BSpline 按 Degree 乘 NbKnots、Plane 与 default 置 1。③ 两处 Knots 口径也核过：OCCT 的 NbUKnots() 是不同节点数，端口写 nb_u_intervals(0) 加 1，即区间数加一等于节点数，正确。④ 结构性发现：基类 BRepMesh_DefaultRangeSplitter::GenerateSurfaceNodes（.cxx:103-107）直接 return Handle(IMeshData::ListOfPnt2d)()，返回空句柄、不加任何内部节点。所以端口在某些面过密、某些面过疏不可能来自基类的节点生成，那里没有逻辑；只能是 delta 与面 pcurve 节点这两个输入量不同，或下游三角化步骤不同。⑤ 因此落点收敛到三个输入量：range_u 与 range_v 即 AdjustRange 之后的离散范围、面的 pcurve 节点即 NodeInsertionMeshAlgo 收集的 wire 点、vertices_nb 即进 AdjustCellsCounts 的那个数。下一轮在端口与 OCCT 各加一个同口径探针，对 f=171（过疏，port 349 对 gt 1569）与 f=206（过密，port 336 对 gt 68）打印 range_u、range_v、len_u、len_v、delta、tol、cells_u、cells_v、vertices_nb 九个数逐值对照；OCCT 侧可在 --uvsum 循环里顺手打印（那里已有 aDF，GetRangeU、GetRangeV、GetDelta、GetTolerance 都在 RangeSplitter 上），端口侧在 param_set.rs 的 adjust_range_base 末尾打印。这是打印同一批中间量的做法，不是数值拟合。本轮未提交任何源码改动，Rust 库零改动。",
      ref: "specs/_a3n00_gap_analysis.md §9.344",
    },
    {
      id: "a48",
      at: "2026-09-30",
      title: "逐式对照：范围/步长前端与 OCCT 逐行等价（第 15 条排除），落点转向节点生成",
      tone: "success",
      detail: "本轮按卡片把端口与 OCCT 的步长/范围路径逐式对照，结论是前端两侧逐行等价（记为第 15 条被排除的假设），落点转向节点生成。逐式对照（OCCT 源在 D:\source\OCCT-src\src\ModelingAlgorithms\TKMesh\BRepMesh\）：① AdjustRange 在 BRepMesh_DefaultRangeSplitter.cxx:45-81，用 updateRange(aSurface->FirstUParameter(), LastUParameter(), IsUPeriodic(), ...) 再算 length、tolerance、delta；端口在 meshing/range_splitter/param_set.rs:791-832，用 surf.u_range() 与 surf.v_range() 加 update_range(gu0, gu1, is_u_periodic, ...)，逐行等价。② computeLengthU 在 .cxx:142-168，端口在 splitter.rs:104-125，两边都是 du = 0.05*(u1-u0)、3 条参数线（V.first、V.ave、V.second）、20 段、除以 3；computeLengthV 同理。所以步长数学与范围调整都不是分歧点。③ 端口自己已写明这一点：param_set.rs:801-812 的注释说明 GetSurface() 是无限制的 BRepAdaptor_Surface（引 IMeshData_Face.hxx:69），其 First/Last U|VParameter 来自 Geom_Surface::Bounds（引 BRepAdaptor_Surface.cxx:79）而非 BRepTools::UVBounds，并解释若改用 Restriction=true 的 adaptor 会让 computeLengthV 为 0、IsValid 为 false；这是一次有依据的忠实移植说明。④ 附带澄清一个口径陷阱：本轮先做的范围对比（range_compare.py）一度显示 160/167 个面范围不同，但该对比无效——端口 zz_uv_feed --model 打印的 urange 来自 surface().u_range() 即曲面域，而 OCCT --uvsum 打印的来自 pcurve 采样点范围（collectWirePoints 口径，即 AdjustRange 之后的 myRangeU），两者是不同量；例如 port f=223 打 [0, 6.283185] 而 UVSUM f=94 打 [4.712389, 6.283185]，同为整周期只是原点不同。⑤ 因此落点转到节点生成：既然 AdjustRange、computeLengthU/V、updateRange、computeDelta 全部等价，分歧只能在把范围铺成节点那一步，可对照的两处是 BRepMesh_DefaultRangeSplitter::generateSurfaceNodes（端口 param_set.rs:525-564 的 get_undefined_interval、init_params_from_intervals、compute_grain_and_filter）与 BRepMesh_NodeInsertionMeshAlgo::generateSurfaceNodes（端口 meshing/node_insertion.rs）。做法：对 f=171（严重过疏，port 349 对 gt 1569）与 f=206（严重过密，port 336 对 gt 68）各取一个，打印节点生成这一步的 u_params 与 v_params 长度和首尾值、delta、tol_u 与 tol_v、最终插入节点数，与 OCCT 同口径逐值对照。本轮未提交任何源码改动，Rust 库零改动。",
      ref: "specs/_a3n00_gap_analysis.md §9.343",
    },
    {
      id: "a47",
      at: "2026-09-30",
      title: "检查点答案：那两对是真实严重少铺且面定义与 GT 相同；步长假设被否证（第 14 条）",
      tone: "success",
      detail: "本轮执行卡片检查点：核实那两对配对正确、确认面定义与 GT 相同、并把成因从「面定义」排除、指向「三角化密度」；顺带否证了步长假设。① 先核实配对：port f=171 对 UVSUM f=1 的 bbox 逐值相同到 1e-6，port f=129 对 UVSUM f=23 同样相同 ⇒ §9.340 的配对不是错配。② 结构对比显示面定义与 GT 一致：port f=171 的 wires=4 与 GT 的 wires=4 相同、bbox 相同、vrange 相同；f=129 的 wires 是 2 对 GT 的 1。而三角数是 port 349 对 gt 1569（比值 0.222）与 port 18 对 gt 1012（比值 0.018）⇒ 不是面建错了，而是端口在同一面上生成的三角远少于 OCCT，故 reader 建面一侧这条线可以关闭（bbox 与 wires 都相同、只有三角数差）。③ 步长假设的检验被否证（记为第 14 条被否证假设）：把 167 个配对面的 GT lenU/lenV（UVSUM 打印的 DefaultRangeSplitter 弧长近似）与三角数做 log-log 相关，得 corr(log(lenU*lenV), log(gt_t)) = +0.5249、corr(..., log(port_t)) = +0.4155，两侧都只是中等且彼此接近；若端口步长口径与 OCCT 有系统性差别应当看到明显不同的相关强度或斜率，实测没有。④ 异常的两侧特征：端口过度细分（比值 2.2 到 4.9）与严重不足（比值 0.014 到 0.29）两类都聚集在狭长面上（GT lenU/lenV 一个极大一个很小），167 个配对面里 25 个长宽比超过 10，但方向不统一（有的过密有的过疏）⇒ 不得把狭长面一定过密或过疏当规则，那是 CLAUDE.md 明禁的长度滤边式启发式。⑤ 下一轮回到门禁第 2 条的正规做法：在端口三角化管线里对 f=171（严重过疏）与 f=206（严重过密）各取一个，打印该面的细分步长与预计格点数，与 OCCT 的 BRepMesh_DefaultRangeSplitter::computeLengthU/V 逐式对照，不要再在数值上拟合。本轮未提交任何源码改动，Rust 库零改动。",
      ref: "specs/_a3n00_gap_analysis.md §9.341 §9.342",
    },
    {
      id: "a46",
      at: "2026-09-30",
      title: "验证推翻上一轮链条：--mesh 与 --uvsum 索引只对 1/226；已重建合法链条并产出修正表",
      tone: "warning",
      detail: "本轮先做卡片定的验证，结果推翻了 §9.338 的做法。实测取 --mesh 与 --uvsum 同一次遍历的同一索引比 bbox，只对 1/226（--mesh i=43 的 bbox 在 x 属于 [-132.57,-112.43]，而 --uvsum f=43 在 x 属于 [112.50,114.50]，位置相反）⇒ 两套顺序完全不同：--mesh 用 TopExp_Explorer(aShape)，--uvsum 用 IMeshData_Model::GetFace(f)。故 §9.338 用 UVSUM 面号去取 --mesh 同号三角数的链条无效，其逐面表（净 -154、亏损 5584、盈余 5738、f=116 由 1012 变 18、f=63 由 559 变 11 等）作废。索引对齐实测：--mesh 对 --uvsum 为 1/226、--mesh 对 PORTID 为 0/226、PORTID 对 --uvsum 为 3/226；按 bbox 集合配对：PORTID 对 --uvsum 为 167/226、PORTID 对 --mesh 为 47/226、--mesh 对 --uvsum 为 55/226（--mesh 的 bbox 取自三角化后的面，与另三者的几何面 bbox 口径不同）。随后重建合法链条：给 --uvsum 加 triangles 字段（在同一个 GetFace(f) 循环里取 BRep_Tool::Triangulation 的 NbTriangles），并补一次 BRepMesh_IncrementalMesh —— 该分支原本只做 ModelBuilder 与 EdgeDiscret 不三角化，新字段初值实测为 0；补上后 --uvsum 0.215 的逐面三角数之和 = 12462，与 --mesh 0.215 0.5 的 TOTAL triangles 完全相等，自洽校验通过。修正后：167 面上 port 8647 对 gt 9002（净 -355，3.9%），最大亏损 f=206（gt 68 对 port 336）与 f=192 至 f=199 八面（各 gt 52 对 port 228），最大盈利 f=171（gt 1569 对 port 349）与 f=129（gt 1012 对 port 18）。方向性结论未被推翻，但 §9.338 的具体面号多数错位，本节的表取代它；port f=171 对 gt f=1 与 port f=129 对 gt f=23 在两次不同配对下都居盈利榜首，不是配对假象。本轮改动仅在 specs/occt_probe/occt_probe.cpp（--mesh 加 tshape 诊断、--uvsum 加 triangles 并补三角化），build.bat 重建成功；crates/occt-topo/src 零改动",
      ref: "specs/_a3n00_gap_analysis.md §9.339 §9.340",
    },
    {
      id: "a45",
      at: "2026-09-30",
      title: "配对改进后：167 面上总量净差仅 1.8%，但逐面短亏 5584 / 盈余 5738 几乎整批互换",
      tone: "success",
      detail: "按卡片定法改进配对：不用 --mesh 的 bbox（取自离散模型、与 OBJ 口径不同，精确命中仅 47），改用 UVSUM 的 bbox（与 PORTID 同为进口形状口径，精确命中 167），三角数由 UVSUM 面号直接取 --mesh 的 FACE i triangles=。新脚本 .target-gate/face_tris_diff2.py。结果：port faces=226 uvsum faces=226 mesh faces=226、精确配对 167、未配 59；port triangles(paired)=8647 vs gt=8801，净 deficit 154（1.8%）；shortfalls 合计 5584、surpluses 合计 5738。端口几乎没铺的面：f=114(gt 2 vs port 228)、f=191(3 vs 376)、f=115(5 vs 243)、f=171(6 vs 349)、f=200(8 vs 376)；端口铺得远多的面：f=116(gt 1012 vs port 18)、f=63(559 vs 11)、f=95(265 vs 11)、f=100(326 vs 11)。三条结论：(1) 精确配对那批面上两侧网格密度总体相当（净差 1.8%），与 §9.337 的「端口不比 GT 粗」一致，把量级钉住——缺口不在这些面铺得太少；(2) 逐面分配差异极大而总量相抵，这解释了为何总量/面积比这类单一指标长期无法定位问题（它是大额双向误差的净值）；(3) T-99 的正确问题最终表述为「端口把三角放到了与 OCCT 不同的面上」——某面过度细分、另一面几乎放弃，成对出现。三条局限已写入卡片：① --mesh 面号与 UVSUM 面号的一一对应未经独立验证（UVSUM 的 f 来自 IMeshData_Model::GetFace(f)，而离散模型会跳面，FaceMeshStat 的文档即言 the mesh model skips faces）⇒ 若不同序则本节配对有一层错位、会放大表观双向偏差，须先验证才能定案；② GT 侧出现个位数三角读数（uvsum f=3→3、f=151→5、f=156→2）而端口同面 228-376，也可能是错位症状；③ 59 个面未配上（占端口 22% 三角）未纳入统计。本轮未提交任何源码改动",
      ref: "specs/_a3n00_gap_analysis.md §9.338",
    },
    {
      id: "a44",
      at: "2026-09-30",
      title: "关键结论：缺口不是参数假象 —— 端口网格在每个可比档位上都比 OCCT 更细",
      tone: "success",
      detail: "针对 §9.336 提出的假说「面积差可由端口网格较粗（lin=1.076）而 GT 较细（d≈0.215）解释」做决定性对照。用上一轮新加的 --lin 覆写在同名义 deflection 下逐档比较端口 zz_uv_feed --ids（mt 之和，与导出 OBJ 三角数逐字一致）与 OCCT probe --mesh：d≈1.0/1.076 → 端口 11121 对 GT 8124（端口更细 1.37×）；d=0.5 → 22616 对 20046（1.13×）；d=0.3 → 42217 对 23652（1.79×）。⇒ 每个可比档位上端口三角数都多于 OCCT ⇒ 端口的网格不比 GT 粗。这否证了那个假说：若端口更粗则三角数应少于 GT，实际多于。而且 occ-a3n00.obj（12466 三角）与 output/a3n00.obj（11121 三角）是同口径产物，所以三角数差与面积差（195950.71 对 227130.30，86.3%）都是真实信号 ⇒ 端口以更细的网格却只有 86.3% 的面积，少铺是真实的。这很重要：它把 T-99 从「等口径决定」的停滞中释放出来（不必再等 GT OBJ 的确切参数即可判定缺口不是端口偏粗）。顺带否证：用三角数拟合 GT 参数不可靠 —— 细扫非单调（0.200→12844、0.205→12560、0.210→12574、0.215→12462、0.220→12368、0.225→12204），系 BRepMesh 内部角度判据与插入顺序所致，不确定度≥±0.005，不足以支撑逐面对拍。保留 §9.336 的事实：两侧 deflection 语义确实不同（1.076 对 8124 三角 vs 0.5 对 20046 三角，量纲不同），且门禁容差注释 tests/common/mod.rs:44 写的 a few percent 与实测 14% 不一致，但都不改变本节结论。下一轮改进配对：改用 UVSUM 的 bbox（GT 进口形状口径，与 PORTID 精确命中 167）替代 --mesh 的 bbox（只得 47），再重做 §9.334 的逐面三角差。本轮未提交任何源码改动",
      ref: "specs/_a3n00_gap_analysis.md §9.337",
    },
    {
      id: "a43",
      at: "2026-09-30",
      title: "前提问题：端口与 OCCT 的 deflection 语义不同（同数值差 4.6 倍），逐面对拍须先定口径",
      tone: "warning",
      detail: "给 zz_uv_feed 加 --lin <d> 覆写（探针本来就有 override_lin，只是把同一能力补到网格侧；该 example 是 untracked 的工作区仪器），然后在同一 deflection 数值下比三角数：OCCT probe --mesh 0.215 0.5 得 12462，端口 zz_uv_feed --ids --lin 0.215 得 57055 ⇒ 差 4.6 倍。⇒ 两者 deflection 不同义，不能在相同数值下逐面对比。两侧口径核查：端口驱动网格的是 Prs3d::GetDeflection = maxComp(bbox)*0.001*4（export_data_obj 的注释明写 changing 0.1 does not，Prs3d.hxx:82-103），实测 a3n00 computed_lin=1.076007；GT 的 occ-a3n00.obj 生成参数不在本仓库（grep 全仓无生成器），只能由三角数反推 d≈0.215。此外角度单位也不一致：端口用 20.0_f64.to_radians()，而 OCCT --mesh 的 angle 参数是弧度（0.5 即约 28.6°）⇒ 两侧两个参数都不同义。这给面积比 0.863（偏离 0.1373）提供了一个「参数不同」的解释（较粗网格在曲面上覆盖略少），但门禁自己的容差注释与此不一致：tests/common/mod.rs:44 写 Curved shapes: a few percent of chord/deflection difference is expected，而 14% 远超 a few percent ⇒ 要么注释已过时，要么 0.8627 里确有一部分不是 chord/deflection 差异；本轮不能判定是哪一种，因为需要 occ-a3n00.obj 的确切口径。因此 T-99 若要继续需要一个口径决定：把 occ-*.obj 的生成参数记录入库，或让端口支持绝对 deflection 与同单位角度。在此之前不做逐面对拍，否则会把参数差当成缺口（§9.334 的精确配对只有 47/226 与此有关）。--lin 的界限：可用于端口内部的收敛性检查，不能用于与 GT 对拍。本轮 crates/occt-topo/src 零改动",
      ref: "specs/_a3n00_gap_analysis.md §9.336",
    },
    {
      id: "a42",
      at: "2026-09-30",
      title: "首次可靠的逐面三角差：缺口是「三角分配」而非总量，且「曲面铺粗」不成立",
      tone: "success",
      detail: "落实 §9.332 的更正：给 specs/occt_probe/occt_probe.cpp 的 --mesh 模式加 bbox 输出（复用同文件 UVSUM 路径的 BRepBndLib::Add + Bnd_Box::Get），build.bat 重建成功（BUILD OK）—— 证明 OCCT 侧确实可扩展。口径扫描：angle 从 0.1 到 π/4 三角数 89226→19404，最少也高于 GT OBJ 的 12466 ⇒ 仅调 angle 无法复现；改扫 deflection，0.20→12844、0.21→12574、0.22→12368 ⇒ d≈0.215 得 12462（GT OBJ 12466，差 0.03%）；注意 nodes 对不上（探针约 10621 对 OBJ 的 11145）⇒ OBJ 的 v 与 Poly_Triangulation 的 NbNodes 不同义，只有三角数 f 与 triangles 同义。逐面对减结果：port 11121 vs gt 12462、deficit 1341，但是亏 3737 / 盈 5078 的大额相抵，不是单向铺少。三条结论：(1)「曲面铺得偏粗」不成立 —— 若成立端口会在所有曲面面上都少于 GT，实际大额双向偏差、同一面型两侧都有；(2) 我此前所有基于面积排序的判断与此可靠口径不符 —— §9.325 判 f=171 为最大亏损面（−20874），按三角数它是 +1220 的盈面（gt 1569、port 349）；f=129（gt 1012、port 18）、f=216（gt 559、port 1）、f=201（gt 469、port 1）几乎被完全放弃，而 f=206（gt 68、port 336）、f=192-199（各 gt 52、port 228）被过度细分；(3) 缺口是三角被放到别的面上的分配问题，不是总量不足。配对质量：精确 bbox 配对只有 47 个（PORTID 对 UVSUM 是 167），因 --mesh 0.215 的 bbox 由另一套参数重建，其余最近邻兜底可能错配 ⇒ 具体面号谨慎，总量与双向偏差结论不依赖个别配对。另查清 occt_probe.cpp 工作区比 HEAD 多 620 行（--wireinv/--typefaces/--adv/--uvsum），是既有未提交改动、非本轮所写（本轮仅 +10/−1 的 bbox 段）；这些仪器有被 git checkout/stash 抹掉的风险，建议由人决定是否提交",
      ref: "specs/_a3n00_gap_analysis.md §9.334 §9.335",
    },
    {
      id: "a41",
      at: "2026-09-30",
      title: "重要更正：GT 侧本来就有逐面三角数据（探针 --mesh），我此前说「拿不到」是错的",
      tone: "warning",
      detail: "本轮核查「能否给 OCCT 探针加字段」时发现我此前多次写下的结论是错的。实测：① D:\\source\\occt-8.0.0\\inc 下有 7084 个 .hxx 头文件（发行版带头文件）；② C:\\Program Files\\Microsoft Visual Studio\\2022\\Professional\\VC\\Auxiliary\\Build\\vcvars64.bat 存在（编译器可用）；③ specs/occt_probe 自带 occt_probe.cpp（1893 行）与 build.bat（vcvars64 + cl + OCCT inc/lib）⇒ 本就具备给探针加字段并重建的条件；④ 而且这字段早就有 —— 探针的 --mesh 模式（occt_probe.cpp:142-188，注释标为 T-54 oracle）逐面输出 FACE <idx> nodes=<n> triangles=<t>，实测 FACE lines=226、TOTAL faces=226 nodes=15538 triangles=20046 meshvol=1.89747e+06 neg_triangles=4609 deflection=0.1 angle=0.5。但 --mesh 的角度是 0.5 rad，与生成 occ-a3n00.obj（12466 三角）和端口 OBJ（11121 三角）的参数不同 ⇒ 20046 与 12466 都是 GT 却不相等，直接比会把参数差当成缺口。因此下一轮先查清 occ-a3n00.obj 的生成参数，用同一组重跑 --mesh，再与 PORTID 的逐面 mt 按 bbox 六坐标配对逐面对比。另：§9.331 的两条否定结果仍成立且独立 —— 两版 OBJ 都没有面积 < 1e-6 的塌缩三角（area<1e-3 各只有 1 个），最大三角完全相同（1514.86），所以缺口不是塌缩三角造成的。本轮未提交任何源码改动",
      ref: "specs/_a3n00_gap_analysis.md §9.331 §9.332",
    },
    {
      id: "a40",
      at: "2026-09-30",
      title: "定位了一个真实缺陷（4 面只铺部分环），但量级不足，已撤回「59 面是主线」的期望",
      tone: "success",
      detail: "按卡片检查点把「面建小了」与「只铺了一小块」两种解释区分开。① 差异分布：59 个不匹配面里 51 个差异 ≤ 0.1（34 个 ≤ 0.001 属舍入量级），>1.0 的只有 4 个（208↔10、210↔13、212↔15、214↔17，全部同值 4.0013）。② 这 4 个的成因已定：zz_uv_feed --model 取端口结构后与 GT UVSUM 逐值比 —— urange 与 vrange 与 GT 逐值相同（到 1e-14）⇒ 面的几何定义完全相同，§9.329 的「面建小了」被排除；端口 2 条 wire 对 GT 1 条（与既有「29 个面 2 vs 1」统计同族），两条 edge 都是整圆 dpar=2π（参数域 [π,3π] 与 [0,2π]）；而 vrange 到 59.69（2.5 圈）却只铺出 y 跨度 18.93 ⇒ 这 4 个面只被铺了「部分环」，即 §9.329.3 里标为无法区分的第二种解释被坐实。③ 量级评估：这 4 个面 mt=451 三角、承载 GT 面积约 2035（bbox 粗估），而全模型缺口 31179.6 ⇒ 至多几个百分点，不足以解释缺口，故撤回「59 个面是主要落点」的期望。缺口主因仍未定位，需要在配对成功的 167 个面里找，但那需要逐面面积口径（给 MeshFace 加三角访问器；OCCT 的 BRepMeshData_Face 持有 Triangulation，有同等分支、非自造），属改库、须完整走五道门禁。新增的这个可复现缺陷（4 个面 + 2 vs 1 wire）是真实且已完全定位的，可作独立小卡，但按五道门禁第 4 条本轮不动手。本轮未提交任何源码改动",
      ref: "specs/_a3n00_gap_analysis.md §9.330",
    },
    {
      id: "a39",
      at: "2026-09-30",
      title: "排除索引假说、量出 59 面差异幅度很小；聚焦到「中心相同、尺寸 0.7045 倍」的 Extrusion 面",
      tone: "success",
      detail: "三件事。① 排除卡片列的第三个候选：PORTID 的 f 连续 0-225、mt 查找失败 0 次，且 FaceMeshStat 本就按 shape_key 设计以防遍历顺序错位（discret_root.rs:77-80）⇒ 索引不一致不是成因。② 量出 59 个不匹配面的 bbox 差异幅度：最坏 4.0、中位约 0.05、最小 0.0033，每坐标都在 1e-9 内的为 0 个 ⇒ §9.327.3 说它们几何真的不同是过重表述，准确说法是同一几何的轻微边界差异；它能解释为何不能用精确 bbox 配对，但不足以单独解释 31179.6 的缺口（§9.327 的 22% 是配对失败面所承载的三角量，不是 22% 的三角几何错误 —— 这个区别我上一轮没讲清楚，已更正）。③ 聚焦线索：差异成组同值（208/210/212/214 同为 4.00133、209/211/213/215 同为 0.073154、109-112 同为 0.0462028）。取 port 208 对 GT 10 按绝对 bbox 量：x 完全相同 [112.5,132.5]、y/z 中心相同（差 0.032），而 y/z 尺寸只有 GT 的 0.7045 倍（extent 18.931 对 26.870）；面型 Extrusion，成组出现说明是同一段构造逻辑在 4 个象限重复。已写明反证：bbox 受该面被裁剪到的 UV 范围影响，port 208 的 mt=1 意味着只铺了一小块也能让 bbox 偏小，两种解释本轮无法区分，故不断言成因。下一轮用 --model/--sph 取 208/210/212/214 的 urange/vrange/wire 结构与 GT 逐值对比。本轮未提交任何源码改动",
      ref: "specs/_a3n00_gap_analysis.md §9.328 §9.329",
    },
    {
      id: "a38",
      at: "2026-09-30",
      title: "正确口径两个都不可行；新落点：59 个面与 GT 的 bbox 精确配不上（占 22% 三角）",
      tone: "success",
      detail: "按 §9.326.4 列的替代口径逐个试：① 给探针加 --facearea 在端口侧不可行 —— 查 MeshFace 的方法只有 face/surface/wires/points/point(i)/deflection/status 等，没有三角容器（逐面三角只在模型级 MeshModel 之上），ShapeMesh 也是整形状的；② 按 GT 面型归因三角数不可行 —— UVSUM 行里没有三角计数（只有 wires/bbox/urange/vrange/dU/dV/lenU/lenV/valid），GT 三角数只能从 OBJ 反推、又回到不可靠归属。两个原因都记下来以免下轮重试。新落点（不依赖任何面积反推，只用两侧各自 dump 的 bbox）：port faces=226 gt faces=226，按 bbox 六坐标精确配对 167 个、59 个配不上（占端口三角 2474/11121 ≈ 22%）。这不是舍入假象 —— GT 的 bbox 确实带 ±1e-10 epsilon（UVSUM f=1 的 -57.360500611 对 PORTID f=171 的 -57.360501，舍入到 6 位相等，这正是 f=171 能配上的原因），舍入到 6 位仍不相等的 59 个说明 bbox 真的不同；可疑面含 f=190、f=138、f=40、f=176、f=201–217 等。下一轮先查这 59 个的成因（范围建错 / 拆分合并不同 / faces_of 顺序与 PORTID f=i 不一致）。若需逐面面积须先给 MeshFace 加三角访问器（OCCT 的 BRepMeshData_Face 持有 Triangulation，有同等分支、非自造），但那是改库、须完整走五道门禁。本轮未提交任何源码改动",
      ref: "specs/_a3n00_gap_analysis.md §9.327",
    },
    {
      id: "a37",
      at: "2026-09-30",
      title: "更正：§9.325 的逐面面积归因无效（全部 226 个面的 bbox 互相重叠），已撤回",
      tone: "warning",
      detail: "做 §9.325 的三个检查点时读到矛盾数字：f=9 两侧分别是 226.99 / 18049.08（同一个面不可能两侧都大）；f=1 一侧为 0 而它的 bbox 与 port f=171 逐字相同。根因实测：port 侧 bbox 互相重叠的面对数 = 967，涉及 226/226 个面 —— 每个三角质心都落在多个 bbox 里，而脚本是「第一个匹配就 break」，归属基本任意。所以 §9.325 的「净缺口可归因到少数面」「主因不是 F113 而是 f=171/f=200/f=191」「f=171 单面占 67%」「亏损约 −54000、盈利约 +66000」全部作废（脚本里唯一可信的是总量，因为每块面积终被派给某个面、与派给谁无关）。仍然可信的是结构字段（各自 dump 的 wires 数与 bbox，不经分桶）：f=171/200/191 的 wires 数与 GT 完全一致（4v4、10v10、10v10）⇒ 它们不太可能是 107 边差集的来源；真正 wires 不一致的是 173(2v1)、113(4v2)、169(2v1)、189(2v1)，与 pair_a3n00.py 报的「端口在 29 个面上有 2 条 wire 而 GT 只有 1 条」一致。下一轮改用正确口径：推荐给 zz_uv_feed / zz_probe_a3n00 加 --facearea 模式，直接读每个 MeshFace 的三角求和（唯一不依赖反推的口径），再按 bbox 与 GT 配对；落点重新从这 29 个 wire 数不一致的面入手。本轮只做核查与更正，未提交任何源码改动",
      ref: "specs/_a3n00_gap_analysis.md §9.326",
    },
    {
      id: "a36",
      at: "2026-09-30",
      title: "T-99 缺口归因：净缺口 30919.7，主因是 f=171（占 67%）而非 F113",
      tone: "success",
      detail: "在按 §9.324 修正的新基线上做第一次逐面面积归因。配对：按 bbox 六坐标（D21）226/226 全配上、两侧 unpaired 皆 0、port unmeshed=0。工具：新脚本 .target-gate/area_deficit.py 把两个 OBJ 的三角按质心落进 PORTID/UVSUM 的 bbox 桶再逐面对减；口径校验 —— port 侧已分配 195950.71、未分配 0.00，与 pair_area.py 的总面积逐字一致。结果：净缺口 30919.69，是「亏 + 盈」的净额。亏损大项：f=171 −20874.45（对 GT f=1，单面占净缺口 67%）、f=200 −9315.27（端口面积 0.00）、f=191 −9244.65（端口 0.00）、f=173 −4505.07、f=180 −2226.20、f=206 −1832.08。盈利大项：f=190 +14778.90（GT 无对应）、f=189 +11641.24、f=169 +9315.89、f=113 +7338.31、f=40 +6290.31、f=138 +5968.43。⇒ 旧结论「缺口主因是 F113 一个面」不成立（F113 排第 4 且方向是盈）。下一轮三个有界检查点：f=171 对 GT f=1 查 wire/边结构；f=200/f=191 端口面积 0 的原因（三角退化 vs 配对有误）；f=190 是否为重复面。附带线索（未验证）：pair_a3n00.py 报端口在 29 个面上有 2 条 wire 而 GT 只有 1 条，与「盈」的大项面方向一致。本轮未提交任何源码改动",
      ref: "specs/_a3n00_gap_analysis.md §9.325",
    },
    {
      id: "a35",
      at: "2026-09-30",
      title: "漂移定案：旧参照来自本会话丢失的未提交改动（我的操作失误）；T-99 目标已按新基线改写",
      tone: "warning",
      detail: "从 git 历史查清根因。决定性证据：f1e56776 的提交信息自记「已知代价与未结：a3n00 面积比 0.8824→0.8627（F113 出网格 mv 0→223 但整体偏低）；T0M 未网格 6→7、wires 1931→2121、多 wire 面 120→309」—— 与当前 HEAD 实测的 226/0/11101-11121/0.8627/F113 mt=267 逐条吻合。而 git log -S '0.8996' 只命中 f1e56776 的父提交 35d16da9，且该提交记的改前值是 0.8824 而非 0.8996 ⇒ 0.8996 与 edges=1083 / 10863-11941 / F113 不出网格 不对应任何已提交状态。它们来自本会话早期工作区里 read_topology.rs 的 T-93/T-99 reader 侧 seam 接线（该文件当时 2036 行且含该块，现在与 HEAD blob 4f79fa0c 逐字相同、1975 行、不含该块）—— 那 61 行从未提交，stash@{0} 不含它（内容是 bop_builder/fillet_curved），git fsck 的 dangling blob 为 0（已 gc），无法恢复。最可能是我在会话中段的某次 git checkout/stash 序列造成的，属我的操作失误并已记录备查。据此按新参照改写 T-99 的 goal（从「让 F113 出网格」改为「把 a3n00 面积比从 0.8627 推向 1、缺口 31179.6」）与 acceptance（加「T0M 未网格 7 不增加」与「不得改基线或 area_tol 绕过」），并重述了三条受影响的旧结论。本轮未提交任何源码改动",
      ref: "specs/_a3n00_gap_analysis.md §9.324",
    },
    {
      id: "a34",
      at: "2026-09-30",
      title: "基线重建完成：当前干净树已经在「F113 出网格 + 面积比 0.8627」这一端，T-99 的前提不成立",
      tone: "warning",
      detail: "在干净树上重测（HEAD 9036d8ab；read_topology.rs blob 4f79fa0c == HEAD，无未提交改动；tgeometry.rs 无改动）：a3n00 stats=226 unmatched=0 mesh_v/mesh_t=11101/11121 面积比=0.862724 abs diff=31179.5867 edges=976，F113 wires=4 mt=267（出网格），连跑 3 次一致；T0M stats=1765 unmatched=7 mesh=49999/47616。对照本会话旧参照（225 / 1 / 10863 / 11941 / 0.8996 / edges 1083 / F113 不出网格）⇒ 旧参照不可复现，而当前树已在本会话「实验组」的那一端。因此「消除共享键写回让 F113 出网格但面积比掉到 0.8627 故不能落地」对当前树没有意义 —— 实验组与当前干净基线是同一个状态，该二分不存在；追了十四轮的「F113 不出网格」在当前代码上不存在。漂移的代码依据：读者的 resolve_edge 自己就把 STEP 的 seam pcurve 存进 registry（read_topology.rs:938-958），seam 信息来自 STEP 本身、不依赖已在 f1e56776 删除的 T-93 reader 侧补偿；F113 能否出网格主要由 pcurve 存到哪个键决定（repr_key 语义）。本轮未提交任何源码改动、未改卡片状态，只重建基线并记录，并在 T-99 的 next 与 nextAction 里写明建议的两种处置（改写 goal 为「把 a3n00 面积比从 0.8627 推向 1」，或按旧参照不可复现暂停本卡），等人工确认",
      ref: "specs/_a3n00_gap_analysis.md §9.323",
    },
    {
      id: "a33",
      at: "2026-09-30",
      title: "停止信号：基线出现不可解释的漂移，本会话 T-99 所依赖的参照系无法在干净工作树上复现",
      tone: "warning",
      detail: "本轮为「记录被采纳面」插桩时发现读数与卡片不符，遂彻底核查。事实：① 工作区的 read_topology.rs 与 HEAD 的 blob 哈希相同（4f79fa0c），即无未提交改动；② 但该文件里没有 T-93 的 reader 侧 seam 块（grep fix_missing_seam / T-99 均为 0），而 git show --stat f1e56776 显示该提交正是删除这段的；③ 当前干净工作树连跑 3 次完全一致：stats=226、matched=226 unmatched=0、mesh_v/mesh_t=11101/11121、PORTID f=113 mt=267（F113 出网格）、面积比 0.862724。⇒ 当前状态就在本会话一节用来把 T-99 推向「已解决」的那组数字上，而不是会话里反复引作基线的 225/1/10863/11941/0.8996。T-99 与 T-94 的全部结论都建立在那条基线之上（F113 病根是共享键写回、实验掉到 0.8627 故不能落地、43 面各少 2–3 条边），因此参照系不成立时那些数字不能继续作为依据；而 read_topology.rs 与 HEAD 相同又与「本会话多次读到该 seam 块」相互矛盾，说明工作树在会话期间发生过未记录的状态变化。本轮刻意未提交任何源码改动、未试新方案、未改卡片状态，只做核查与记录，并把 nextAction 改为「先重建基线」",
      ref: "specs/_a3n00_gap_analysis.md §9.322",
    },
    {
      id: "a32",
      at: "2026-09-30",
      title: "T-99 插桩收窄范围：swap_seam 读取路径 0 次调用（推翻 §9.296.3），check_wire 早退两版都是 0",
      tone: "success",
      detail: "两处插桩各否决一条假设。① wire_data.rs 的 swap_seam 入口（OCCT_TOPO_TRACE_SWAPSEAM）：基线 0 次调用、实验 0 次调用 ⇒ swap_seam 在 read_step_file 路径上从未被调用，故 §9.296.3 的「D2 的共享 registry 改写其执行者就是 swap_seam」不成立；§9.296.1 的静态结论没错，但那处代码在读取路径不可达，而 §9.296.9/§9.314 的观测是探针单独构造 fix_missing_seam() 调用得来的，与读取路径不是同一条。② check_wire 的 curve_on_surface_oriented == None 早退处（OCCT_TOPO_TRACE_CHECKWIRE）：基线 0、实验 0 ⇒ 我据此提出的「共享键的第二个作用是让 check_wire 找到 pcurve」被否决。合起来：共享键写回确实影响结果（stats 225→226、unmatched 1→0、edges 1083→976，已由 §9.315 与 §9.320 两次独立实现证实），但既不经 swap_seam、也不经 check_wire 早退。最强线索是 unmatched 1→0 —— 它意味着 F113 在两版里是不同的面对象、即 resolve_face 返回值不同。下一轮据此量「哪些面被采纳与否」的差集。三处插桩与实验全部 git checkout 回退，三文件 diff 为空、grep 无残留，基线复验 stats=225 unmatched=1，库代码本轮零改动",
      ref: "specs/_a3n00_gap_analysis.md §9.321",
    },
    {
      id: "a31",
      at: "2026-09-30",
      title: "T-99 核心改动否决：把 repr_key 改为面身份与摘注册实验逐字同结果，单纯改键归属无效",
      tone: "success",
      detail: "按 §9.319.4 把 tgeometry.rs 的 repr_key 从「面的曲面数据指针」改为「面自身的键」（fn repr_key(&self, face_key: usize) -> usize { face_key }），cargo check 通过。a3n00 实测：stats 225 到 226、unmatched 1 到 0、mesh_v/mesh_t 10863/11941 到 11101/11121、面积比 0.899606 到 0.862724、绝对差 22802.58 到 31179.59 —— 三个指标与摘注册实验（§9.315）逐字相同。zz_pcurve_key_probe 进一步确认 edges 976（基线 1083）、edges_without_pcurve 448 两版相同、同样是 43 个面各少 2–3 条边（face=41 4→2、face=42 4→2）。所以改用面键等价于摘掉共享键，两者都阻止 seam 写入抵达原面槽位，得到同一组后果。这个否决把问题问得更准：写回同时承担两件事，而这 43 个面的 seam 边之所以存在，依赖原面槽位里有那个 seam pcurve —— 真修复不能只动键的归属，得看 fix_missing_seam 为什么需要那份 pcurve 才能选中 wire 对并加上 seam 边。另记一条自我更正：我曾推断「变体的 976 里不含 F113 的 22 边」，实测两版 F113 都是 22 边，该推断作废。改动已 git checkout 回退，基线复验 stats=225 unmatched=1，库代码本轮零改动",
      ref: "specs/_a3n00_gap_analysis.md §9.320",
    },
    {
      id: "a30",
      at: "2026-09-30",
      title: "T-99 因果链打通：抑制写回不是解耦，而是切断 seam 步骤自身的读写通路（43 面各少 2–3 条 seam 边）",
      tone: "success",
      detail: "把 zz_pcurve_key_probe 加上 OCCT_TOPO_KEY_DUMP 逐面边键 dump 后对比基线与重建实验。形状键是指针值、跨进程因 ASLR 不稳定，故改比逐面边数：基线 faces=226 total_edges=1083，实验 976，delta −107；差异落在 43 个面上，每个少 2–3 条（4→2 或 5→2，face=138 是 15→13）。这与 §9.316.1 的「实验版少铺 8377」方向一致，且少掉的正是 fix_missing_seam 本该加上的 seam 边（对照 §9.290/§9.292：该步骤会把单闭合边面从 1 条扩到 4 条）。原因：fix_missing_seam 构造 seam 时要读回自己刚写下的 pcurve（check_pcurves_and_shift / curve_on_surface_range 都按 repr_key 取），而摘掉 tmp_f 的曲面注册后写落在 tmp_f 自己的键上、读取侧却仍按 repr_key(tmp_f)（=原面曲面指针）去找 ⇒ 读不到 ⇒ seam 边加不上 ⇒ 面退化。所以症结被精确到一处：tgeometry.rs:524-531 的 repr_key 把面身份与曲面身份混为一谈，写则污染原面、不给共享键则 seam 读不到自己的写入，两难；而 OCCT 把 pcurve 挂在边的 BRep_TEdge 上、按 (edge, face) 定位、与曲面无关，本可并存。修复方向因此是让 pcurve 的键跟随面身份（即 §9.296.7 的 A 方案，现有必要性论证）。实验已 git checkout 回退，基线复验 stats=225 unmatched=1 edges=1083，库代码本轮零改动",
      ref: "specs/_a3n00_gap_analysis.md §9.319",
    },
    {
      id: "a29",
      at: "2026-09-30",
      title: "T-99 不变量检查：判据过强（基线即不满足），但把落点收敛到 107 条边",
      tone: "success",
      detail: "新建探针 zz_pcurve_key_probe 检查「每个面的每条边是否都在 repr_key(F) 下有 pcurve」：基线 a3n00 得 faces=226 edges=1083 edges_without_pcurve=448 faces_with_missing=95（T0M 8630/3674/719）⇒ 该不变量在基线上就不成立，我上一轮的表述过强，作为否定结果记下。但对照实验版得到有信息的差异：两版 edges_without_pcurve 同为 448、faces_with_missing 同为 95，而 edges 为 1083 对 976 ⇒ 基线比实验版多 107 条边且这 107 条全部带 pcurve。所以落点不是「缺 pcurve」而是「这 107 条边的 pcurve 被挂到错误的键下」。下一轮据此只量这个有界差集，并检查这批边所属的面是否被采纳。实验已 git checkout 回退，基线复验 stats=225 unmatched=1 edges=1083，库代码本轮零改动",
      ref: "specs/_a3n00_gap_analysis.md §9.318",
    },
    {
      id: "a28",
      at: "2026-09-30",
      title: "T-99 方向③交付：实验版面积差是「少铺 8377」而非多铺错误面，成因指向「面的表示键不匹配」",
      tone: "success",
      detail: "把基线版与实验版各导出一份 OBJ 逐三角对比：基线 tris=11941 area=204327.7166（对 GT 227130.2998 缺口 22802.58）、实验 tris=11121 area=195950.7131（缺口 31179.59）。两版最大的三角完全同名同值（1514.8633、1326.7472、1322.0572、四个 1069.1455、903.1153、727.9473），说明不存在「实验版多出一张巨大错误面」。实验版的 v 更多（11101 对 10863）却 f 更少（11121 对 11941）、总面积少 8377 ⇒ 它在整体上铺得更粗、覆盖更少，是 under-cover 而非多铺。由此把 T-99 的剩余工作重述为一个可判定的不变量：对每个被采纳的面 F，F 的每条边是否都在 repr_key(F) 下有一条 pcurve？若不成立就是「写入的键与面的表示键不匹配」，修复走方向 ②（形状级替换表语义）。两个对照 OBJ 留档在 .target-gate/，不必重跑实验。实验已 git checkout 回退，工作区干净、基线回到 stats=225 unmatched=1，本轮库代码零改动",
      ref: "specs/_a3n00_gap_analysis.md §9.316",
    },
    {
      id: "a27",
      at: "2026-09-30",
      title: "边界实验定案 T-99 病根：共享键写回污染 F113；但消除它会让面积比 0.8996 掉到 0.8627",
      tone: "success",
      detail: "按 A 方案最小实现做边界实验：新增 tgeometry.rs 的 unregister_face_surface（只摘临时面自己的曲面注册，不动 edge/face 几何），并在 shape_fix_face.rs:453 建 tmp_f 后调用一次，使 repr_key(tmp_f) 回退为 tmp_f 自身键。三条效果：① 写回被完全消除（zz_share_probe 从 T0M 21/2/8、a3n00 3/0/0 变为两边都 0/0/0），§9.296.6 的构造级证明获实验级确认；② a3n00 unmatched 1 到 0、F113 出网格 mt=267、全 226 面 mt=0 为 0 个、mesh_v/mesh_t 10863/11941 到 11101/11121 —— 追了九轮的 F113 出网格了；③ 面积比 0.899606 掉到 0.862724、绝对差 22802.6 到 31179.6，判据 |ratio-1|<0.15 下裕度从 0.0496 缩到 0.0127 ⇒ 否决性，不能落地。结论：该写回是承重的 —— 消除它就等于让 fix_missing_seam 的结果真正生效、改变哪些面被返回，所以「中性解耦」不存在；C 方案（只改 swap_seam）的前提也因此存疑。已 git stash 回退，工作区与实验前逐字一致。T-94 结项，T-99 从「搁置待外部证据」重开为 P1（问题重新表述为「找既消除写回又不变差面积比的机制」）",
      ref: "specs/_a3n00_gap_analysis.md §9.315 · T-94 / T-99",
    },
    {
      id: "a26",
      at: "2026-09-30",
      title: "T-94 步骤①完成：实测 2 面/8 槽位、两种机制；计数从 1/24 更正为 2/24",
      tone: "success",
      detail: "升级 zz_share_probe 快照（从只比 n=/rng= 改为逐槽位打印 pcurve 指纹：域 + d0 端/中点 + is_line）后重跑 T0M，可复现地得到 faces_with_fix=21 faces_changed=2 slots_changed=8 —— 此前记的「1 面/4 槽位」是因为旧快照对「pcRange 相同但端点被平移」的槽位漏检。两种机制被区分开：face 1752 是真交换（slot1↔slot4、slot2↔slot3 两两互换，swap_seam 的特征，与 §9.296.6 构造级证明一致）；face 1537 是四槽位 u 各恰好 +2π（delta/2π=1.000000，点序未换），来自 shape_fix_face.rs 的 adjust_by_period，与 swap_seam 无关。网格状态：zz_uv_feed --ids 显示 f=1537 mt=10、f=1752 mt=15，全量 1772 面 mt=0 为 0 个 ⇒ 两面都出网格，改动不影响网格产出。对方案的影响：只改 wire_data.rs 的 swap_seam（C 方案）只能消掉 face 1752 一半，C 不足以完全解耦。定级 P3 维持",
      ref: "specs/_a3n00_gap_analysis.md §9.314",
    },
    {
      id: "a25",
      at: "2026-09-30",
      title: "T-98 结项：锥面 pcurve 报 (-inf,inf) 与 OCCT 的 Geom2d_Line 逐字一致，端口正确",
      tone: "success",
      detail: "按卡面要求做 OCCT 语义对照，三条原文定案：① Geom_ConicalSurface.cxx:207-215 的 Bounds 就是 V1=-Infinite/V2=Infinite；② Geom2d_Line.cxx:142-150 的 First/LastParameter 返回 ∓Precision::Infinite()（hxx:97-100 文档亦明写 RealFirst/RealLast）；③ GeomProjLib.cxx:81 在 GeomAbs_Line 分支原样返回未修剪的 new Geom2d_Line(Proj.Line())，只有输入是 Geom_TrimmedCurve 才修剪（:118-128）。锥面母线 V 恒定 ⇒ 投影是常 V 直线 ⇒ OCCT 同样报 (-inf,inf)。端口 line.rs:37-38 逐字一致 ⇒ 不修，走卡面第一条验收路径（保持现状），未改 mock/门禁。真正教训在探针侧：上轮把 d0 读成「全 NaN」是用自报的 -inf 端点求值造成（D14 已撤回）；端口早有 is_finite 守卫，正确口径是用 3D 边的参数域。顺带修掉一个真实文档错位：make_pcurve_full 的文档被粘进了 select_forward_seam，导致前者反而没文档，已拆开并把「parameterized over the edge's range [a,b]」这句不准确的描述改准。门禁：cargo check 0 error、occt-topo 1255/26（与改动前相同）、a3n00 与 T0M 基线未变",
      ref: "specs/_a3n00_gap_analysis.md §9.313",
    },
    {
      id: "a24",
      at: "2026-09-30",
      title: "T-51 结项：gcpnts 的两条 arm 早已实现；「uniform_abscissa 外形不同」被证不成立",
      tone: "success",
      detail: "核对 5 处 UNPORTED 注释：2 处真实（UniformDeflection 未移植、adv 版 AbsComposite 未移植）、2 处过期（gcpnts.rs 声称 LengthParametrized/AbsComposite 未移植，而同函数的 match 三分支全在；另一处内联注释说用非-adv walk，而紧接的代码正是调 compute_abs_composite）、1 处指向性引用。关键证据：GCPnts_UniformAbscissa.cxx:510 用 anAbscissa = aL / (theNbPoints - 1)，即 OCCT 以「点数-1」为段数，而端口 uniform_abscissa(c,n) 产 n+1 个 ⇒ n = NbPoints-1 时逐字等价，只是入参命名口径不同（段数 vs 点数）。真实余量因此收敛为 1 个缺口：adv 版 AbsComposite（CPnts_AbscissaPoint.cxx:187-295），它需要 CPnts_MyRootFunction 的容差重载。改动为 3 处纯注释（零行为变更）。门禁：occt-geom 143/0、occt-topo 1255/26（与改动前相同）、a3n00 基线逐字未变",
      ref: "specs/_a3n00_gap_analysis.md §9.312",
    },
    {
      id: "a23",
      at: "2026-09-30",
      title: "T-97 结项：不删死代码，改为给 A 类模块标注；死代码警告 390 到 111",
      tone: "success",
      detail: "用户追问「BOPAlgo 族的忠实移植在 C++ 里哪里被使用」后，两次误删（第 1 批 31 个函数、第 2 批整个 bopalgo_tools_wires 503 行）全部逐字回退，净删除 0 行。核对 OCCT 侧确认：BOPAlgo_Tools::EdgesToWires / WiresToFaces 是 hxx:151/:170 的 Standard_EXPORT，被 BRepFill_AdvancedEvolved.cxx:561,641,642,1499,1500 等调用，只是该消费方在端口里尚未移植 ⇒ 定性为「忠实移植、消费方未移植」，不是死代码。落地：对 60 个模块头自报 OCCT 来源的文件各插入统一 6 行块（5 行待接线说明 + #![allow(dead_code)]），每个文件 diff 恰好 +6/-0。效果：死代码 390 到 111、总警告 513 到 220（-57%），未用 import 47 与局部赋值未读 35 不变（说明只压死代码、未掩盖别的类别）。口径更正：此前「561 条」是 rtk 过滤代理混合两个 crate 的产物，改用直接 cargo 后 occt-topo 自身为 513 条。门禁：cargo check 0 error；a3n00 与 T0M 基线逐字未变；--lib 1255/26 与改动前完全相同（26 条为 %TEMP% 环境性失败）",
      ref: "specs/_a3n00_gap_analysis.md §9.309 §9.310 §9.311",
    },
        { id: "a22", at: "2026-09-30", title: "T-96 结项：结论是「不重做 TrimmedCurve 参数表示」—— 两条 gcpnts arm 早已移植并接线",
       tone: "success",
      detail: "否证了 T-51/T-96 的原始前提。① `compute_type`(`gcpnts.rs:66`, cxx:26-65) 与 `compute_abs_composite`(`:102`, cxx:96-158) 都已移植，并经 `:491-501`/`:579-585` 的 `match` 三分支接线，`LengthParametrized` 与 `AbsComposite` 都可到达（模块头却写着 UNPORTED，已按实测改写）；② `[0,1]` 是**自洽重参数化**：`d1/d2/d3/eval_dn` 的 `s^N`(`trimmed.rs:45-79`)、`parameter_intervals` 结点重映射(`:130-140`)、`resolution` 的 `/s`(`:147-153`) 用同一个 `s`，而这三处正是三条 arm 唯一读取的量，故数值等价。③ 若真改成 basis 参数：12 处构造 + ~8 处显式 `[0,1]` 换算，且消费方按 `[0,1]` 喂参数 ⇒ 改后**静默给出错误几何**（不编译失败），判为高风险无收益 ⇒ 不做。T-51 因此重开（其豁免理由不成立），改为核对真实余量（NbPoints 语义 / uniform_abscissa 外形）。门禁：occt-geom --lib **143/0**、cargo check 0 error；库源码改动仅 gcpnts.rs 一处注释",
      
      ref: "specs/_a3n00_gap_analysis.md §9.295 · occt-geom/src/gcpnts.rs:66 :102 :491 :579 · occt-geom/src/trimmed.rs:45,130,147 · T-51",
    },
        { id: "a21", at: "2026-09-30", title: "T-59 结项：a3n00 配对表出齐，面积缺口定到单面 F113；D18 的一行修复被独立否证",
       tone: "warning",
      detail: "226/226 面按 bbox 六坐标全部配上（D21 口径）；端口只剩 1 个未网格面 F113（wire **4 条 [6,14,1,1]** vs GT **2 条 [22,6]**），该面区域占 +26952 面积缺口（pad 0→2 只变 <2%），而全模型净缺口仅 22802.58 ⇒ F113 是唯一主因，已拆为 **T-99**。wire 数不匹配全模型只有 3 组 ⇒ 缺陷是外科式的。**D18 的 `n1→n2` 修复在 D27 之后重测仍为净负面**：a3n00 面积 204327.72→202232.18（比值 0.8996→0.8904）、未网格仍 1、F113 没被救活 ⇒ 已回退，与 D19 一致。a3n00 新比值 **0.8996**（判据 |ratio-1|<0.15，余量 0.0496），卡片里 0.8627 的旧记已过期",
      
      ref: "specs/_a3n00_gap_analysis.md §9.291 · .target-gate/pair_a3n00.py · T-99",
    },
        { id: "a20", at: "2026-09-30",
      title: "T-93 结项：导入期 seam 修复接线（未网格 7→6、wires 2121→1931、WIREHIST 向 GT 收敛），全门禁保持绿",
       tone: "success",
      detail: "先证同构：`zz_seam_fix` 对 T0M f=1691 产出 4 条边（底圆 [π,3π] + 两条 seam [1.080839000541168,1.570796326794897] + deg 顶点边）、bbox z 由 [2.744,2.744] 变 [2.744,3.244]，与 OCCT 成形后逐项相同。再证缺口是调用点：临时计数器量到 `read_step_file` 期间 `SEAMCALLS = 0`（源码里「UNPORTED seam 构造」注释是过期的）。接线后 `--lib` **1281/0**、`step_obj_gates` **5/5**（588.92s / 704.02s 两次），**未改任何断言或基线**。副作用：a3n00 unmatched 0→1、mesh_t 11121→11941（仍在门禁内）。已知有界缺口：f=1758 是圆锥，`FixMissingSeam` 四个分支都不匹配（cxx:1972-1975 return false）",
      
      ref: "specs/_a3n00_gap_analysis.md §9.290 · crates/occt-topo/src/step/read_topology.rs · D27",
    },
        { id: "a1f", at: "2026-09-30",
      title: "T-93 主线找到：导入期 seam 修复（FixMissingSeam）——STEP 里该球带只声明 1 条边，OCCT 成形后有 4 条",
       tone: "warning",
      detail: "文本级取证：T0M.stp #29920(球 R=4.25) → #4223 → EDGE_LOOP #6320 = (#25460) ⇒ 该面**只声明 1 条闭合边**（圆 R=2@z=2.744）。OCCT 成形后同面 4 条边（底圆 + 两条 seam + 退化顶点边），bbox z 由 [-5.256,3.244] 收到 [2.744,3.244]。端口同面 = 1 条边，与未修复的环边数一致 ⇒ 缺的正是 ShapeFix_Shape::Perform → ShapeFix_Face::Perform(FixMissingSeam)。端口 fix_missing_seam 单独作用该面**能** 1→4，但整形状跑一遍 212/1772 面变动、其中只有 f=1691 这一个是真的（126 个 2→4 属 `shape_fix_face.rs:147-148` 的三处 UNPORTED 造成的伪增量）",
      
      ref: "data/occ/T0M.stp #29920/#6320 · ShapeFix_Face.cxx:492-498 :1722-2330 · D25",
    },
        { id: "a1d", at: "2026-09-30", title: "T-69 收口：未网格面成因是「边界环缺边」，不是 RangeSplitter 判据；D15 配对表被撤回",
       tone: "danger",
      detail: "新增 OCCT 探针 `--uv/--uvsum`（走 BRepMesh_ModelBuilder → BRepMesh_EdgeDiscret 后读 IPCurveHandle，即 collectWirePoints 的真实输入）与端口探针 `zz_uv_feed`。按 bbox 六坐标配对得真孪生：P1691→G409（端口 1 wire/1 edge/常 V/bbox Z 零厚 vs GT 1 wire/**4 edge**、2.744–3.244）、P1758→G566、P405-408→G73-76（GT 各 4 边环）。用同一参数域条件可证 OCCT 拿到端口这批点也会拒 ⇒ 判据无罪。D15 的 P405→G407/P1691→G1697/P1758→G1764 经查 bbox 相差 60–180，是错配（中心最近邻未校验尺寸）",
      
      ref: "specs/_a3n00_gap_analysis.md §9.288 · D21 · D22",
    },
        { id: "a1e", at: "2026-09-30", title: "T-93 主线重定向为「补回缺边」；库内插桩全部摘除、行为逐字不变",
       tone: "warning",
      detail: "§9.287.4 那行修复（未网格 7→3）降为旁支：它换的是 notch 边对、不补缺边，且要付 a3n00 −0.0329（D19）。摘除 node_insertion.rs / param_set.rs 的 env-gated dump 后复跑 `zz_uv_feed --ids`，1772 面与同一 7 个被拒面（index+bbox）逐字相同，仅指针型 shape_key 串不同 ⇒ 插桩惰性、结论可归因到未插桩的二进制。本轮库源码零改动、未新写测试、未改断言",
      
      ref: "crates/occt-topo/src/meshing/node_insertion.rs · crates/occt-topo/src/meshing/range_splitter/param_set.rs · D23 · D24",
    },
        { id: "a1a", at: "2026-09-29", title: "看板自审 + 修死锁：T-93 的 dependsOn 指向不存在的 T-59",
       tone: "danger",
      detail: "T-59 只写在散文里、从未立卡 ⇒ P0 的 T-93 挂在一个永不解除的依赖上。已补立 T-59（mesh/P1，a3n00 面积缺口），并给看板加了「悬空依赖」独立风险与 startTurn 的「前置」字段",
      
      ref: "specs/board.canvas.tsx · D20",
    },
        { id: "a1b", at: "2026-09-29", title: "T-93 收口：fix_dummy_seam 传参错位（真 bug），一行修复因 a3n00 面积门禁回退，落 §9.287",
       tone: "warning",
      detail: "fix_notched_edges 传 n1 而 cxx:4019 传 i ⇒ FixDummySeam 合并/重指错误边对，造出 F171 的中点顶点。修后 T0M 未网格面 7→3，但 a3n00 面积比 0.8627→0.8298 越界（23 模型只有它变差）⇒ 按红线回退、不重定断言",
      
      ref: "crates/occt-topo/src/shhealing/wire_fix.rs:3084 · specs/_a3n00_gap_analysis.md §9.287",
    },
        { id: "a1c", at: "2026-09-29", title: "T-94 与 T-98 各得一个否证结论（两轮均未提交）",
       tone: "info",
      detail: "T-94：把 seam-shift 的 set_edge_pcurve 推迟到 result 确定后提交 —— 实测与基线逐字相同，惰性，已回退；T-98：F1758 的「d0=NaN」是探针用 -inf 端点求值的假象，pcurve 实测误差 3.0e-12，卡由 P1 降 P3",
      
      ref: "crates/occt-topo/src/shhealing/shape_fix_face.rs:480-487 · crates/occt-topo/src/pcurve_full/surface_projector.rs:2876-2891",
    },
        { id: "a0", at: "2026-09-29", title: "T-69 本轮：量到 RangeSplitter 零宽点集，6 个丢面分成两个独立缺陷",
       tone: "warning",
      detail: "① F1691/F1758：所有 pcurve 点的 V 是同一个常数（surfV 非退化）；② F405-408：所有点的 U 都等于 surfU 上界。§9.270 的 pcurve 归一化假设同期被否证（未网格仍 7、WIREHIST 逐字不变）",
      
      ref: "crates/occt-topo/src/meshing/range_splitter/param_set.rs:830-833",
    },
        { id: "a00", at: "2026-09-29", title: "T-69 几何键配对完成：端口拒的 7 个面在 GT 全部出网格（真差异已证实）",
       tone: "danger",
      detail: "P405→G407 tri=185、P406-408→G408-410 各 185、P1691→G1697 tri=73、P1758→G1764 tri=79、P116→G116 tri=226；dist 0.0–0.30，旧序号法偏差 4–6。分歧定位到「喂进 RangeSplitter 的离散点集」",
      
      ref: "specs/occt_probe/occt_probe.cpp --pair（本轮临时加入，已回退）",
    },
        { id: "a0b", at: "2026-09-29",
      title: "T-69 单面复现：F1691/F405 的 pcurve 几何正确、退化来自边界本身；F1758 的「NaN」是探针假象（D14 已撤回）",
       tone: "warning",
      detail: "F1691 单闭合边（V1==V2）、pcurve 常 V=1.080839001、误差 7.1e-15、pcRange==edgeRange；F405 两条边 pcurve 常 U=-1.141432027、误差 8.1e-8/8.2e-4；F1758 pcRange 报 (-inf,inf) 但 d0 有限且误差 3.0e-12。另确认 T0M 的 STEP 无 PCURVE 记录 ⇒ 这些 pcurve 全是端口自造",
      
      ref: "crates/occt-topo/src/pcurve_full/make_pcurve.rs:128-139",
    },
        { id: "a0c", at: "2026-09-29", title: "本轮门禁实测（纯诊断、探针已全部回退）", tone: "success",
      detail: "五 crate check 0 error；occt-topo --lib 1281/0（8.58s）；step_obj_gates 5/5（271.24s）—— 与基线一致，无劣化",
      
      ref: "crates/occt-topo/tests/step_obj_gates.rs",
    },
        { id: "a1", at: "2026-09-28", title: "T-93 (a) 落地：移除 reader 侧 FixMissingSeam，F113 出网格",
       tone: "info\", detail: \"代价量化：T0M 未网格面 6→7、wires 1931→2121（§9.269）",
      
      ref: "crates/occt-topo/src/shhealing/shape_fix_face.rs",
    },
        { id: "a2", at: "2026-09-28", title: "T-69 矛盾解开：§9.265 的「190 个面丢 wire」是补偿副作用",
       tone: "info\", detail: \"补偿即使丢弃 Shell(5) 也改写共享 GeometryRegistry，使解析期探针与最终模型分叉",
      
      ref: "specs/_a3n00_gap_analysis.md §9.268",
    },
        { id: "a3", at: "2026-09-28", title: "§9.270 给出有界路径：只加回 pcurve 归一化（不含 fix_missing_seam 本体）",
       tone: "success\", detail: \"把 T-93/T-69 的张力从「需决策」降为「一次实测」",
      
      ref: "specs/_a3n00_gap_analysis.md §9.270",
    },
        { id: "a4", at: "2026-09-28", title: "T-93 根因定案：亏损来自单面 F171 的第 4 条 wire 少一条边",
       tone: "warning\", detail: \"端口用的端点 P 不在 STEP 里，且是边 #7405 曲线的中点 ⇒ reader 另建顶点",
      
      ref: "specs/_a3n00_gap_analysis.md §9.276–§9.279",
    },
        { id: "a5", at: "2026-09-28", title: "T-41 (a) 完成：repair/report 改派 + 夹具按合法布尔重设计",
       tone: "success\", detail: \"改派后 --lib 1281/0、step_obj_gates 5/5、phase10/19/3/4 = 8/5/4/9",
      
      ref: "crates/occt-topo/tests/bop_builder2_boss.rs",
    },
        { id: "a6", at: "2026-09-28", title: "T-28 结项（步 1–4 完成），步 5 经实测判为有害不做",
       tone: "success\", detail: \"统一入口 PatchIntersection 会让 general_bspline_sphere 曲线离球面 ≈0.67",
      
      ref: "specs/_design_architecture_t25_t28.md · §9.250",
    },
        { id: "a7", at: "2026-09-28", title: "T-25 结项、T-11 / T-51 豁免、T-29 收尾", tone: "success",
      detail: "四个清理/架构项同批关闭，各自证据见 tasks 的 evidence",
      
      ref: "specs/_a3n00_gap_analysis.md §9.249",
    },
        { id: "a8", at: "2026-09-28", title: "STEP→OBJ 三门禁合并为 step_obj_gates", tone: "info",
      detail: "23 模型共享 OnceLock 缓存，2168s → 660s",
      ref: "commit 08c4edf2",
    },
  ],
} as const;

// ---- 派生用的常量与纯函数（只在这里，渲染期不手抄任何结论）----------------
const STATUS_TONE: Record<string, Tone> = {
  pending: "neutral",
  in_progress: "warning",
  blocked: "danger",
  completed: "success",
  cancelled: "neutral",
};

const STATUS_LABEL: Record<string, string> = {
  pending: "待办",
  in_progress: "进行中",
  blocked: "阻塞",
  completed: "完成",
  cancelled: "取消",
};

const PRIORITY_TONE: Record<string, Tone> = { P0: "danger", P1: "warning", P2: "info", P3: "neutral" };
const LANE_TONE: Record<string, Tone> = { bop: "danger", mesh: "warning", "port-gap": "info", arch: "neutral", hygiene: "neutral" };

const DAY_MS = 86400000;

type TaskRow = {
  id: string;
  lane: string;
  title: string;
  status: Status;
  priority: string;
  owner: string;
  progress: number;
  estimate: number;
  actual: number;
  startedAt: string;
  updatedAt: string;
  completedAt: string;
  blocker: string;
  goal: string;
  next: string;
  acceptance: string;
  evidence: string;
  write: string;
  ref: string;
  note: string;
  dependsOn: readonly string[];
};

function statusOf(row: { status?: string }): Status {
  const s = row.status;
  return s === "pending" || s === "in_progress" || s === "blocked" || s === "completed" || s === "cancelled" ? s : "pending";
}

function isOpen(task: TaskRow): boolean {
  return task.status === "pending" || task.status === "in_progress" || task.status === "blocked";
}

function progressOf(task: TaskRow): number {
  return task.status === "completed" ? 100 : task.status === "cancelled" ? 0 : task.progress;
}

function priorityRank(priority: string): number {
  return priority === "P0" ? 0 : priority === "P1" ? 1 : priority === "P2" ? 2 : 3;
}

function daysBetween(from: string, to: string): number | null {
  const a = Date.parse(from);
  const b = Date.parse(to);
  if (Number.isNaN(a) || Number.isNaN(b)) return null;
  return Math.round((b - a) / DAY_MS);
}

/// `dependsOn` 里那些**不在板上**的 id：依赖不存在就永远不会自己解除，属死锁。
function danglingIds(task: TaskRow, all: TaskRow[]): string[] {
  return (task.dependsOn ?? []).filter((id) => !all.some((candidate) => candidate.id === id));
}

export default function ProjectBoard() {
  const dispatch = useCanvasAction();
  // 筛选器与当前选中项是纯 UI 态：agent 不需要知道，刷新后保留即可。
  const [filter, setFilter] = useCanvasState<string>("filter", "open");
  const [activeId, setActiveId] = useCanvasState<string>("active", "T-94");
  // 人的状态改动必须走 overlay：它落在 sidecar 里，agent 下一轮 canvas_read 就能看到。
  const overlay = useCanvasOverlay("tasks", DATA.tasks);

  const tasks = useMemo<TaskRow[]>(
    () =>
      overlay.items.map((row) => ({
        ...row,
        status: statusOf(row),
        progress: row.progress ?? 0,
        blocker: row.blocker ?? "",
        dependsOn: row.dependsOn ?? [],
      })),
    [overlay.items],
  );

  // —— 派生分析：所有数字都从 DATA 现算，不手抄结论 ——————————————————
  const open = tasks.filter(isOpen);
  const closed = tasks.filter((task) => !isOpen(task));
  const wip = tasks.filter((task) => task.status === "in_progress");
  const blocked = tasks.filter((task) => task.status === "blocked");
  const done = tasks.filter((task) => task.status === "completed");
  const counted = tasks.filter((task) => task.status !== "cancelled");
  const estTotal = counted.reduce((sum, task) => sum + task.estimate, 0);
  const weighted = counted.reduce((sum, task) => sum + task.estimate * progressOf(task), 0);
  const progressPct = estTotal === 0 ? 0 : Math.round(weighted / estTotal);
  const overrun = counted.filter((task) => task.actual > task.estimate);
  const stale = open.filter((task) => {
    const age = daysBetween(task.updatedAt, DATA.asOf);
    return age !== null && age > DATA.staleDays;
  });
  const recentDone = done.filter((task) => {
    const age = task.completedAt === "" ? null : daysBetween(task.completedAt, DATA.asOf);
    return age !== null && age <= 7;
  });
  const gated = DATA.gates.filter((gate) => gate.status !== "green" && !gate.acked);
  const gateReds = DATA.gates.filter((gate) => gate.status === "red");

  const risks: string[] = [];
  if (blocked.length > 0) risks.push(blocked.length + " 条阻塞：" + blocked.map((task) => task.id).join("、") + "（在等 " + blocked.map((task) => task.blocker.split("，")[0]).join("；") + "）");
  if (wip.length > DATA.wipLimit) risks.push("在办 " + wip.length + " 条，超过 WIP 上限 " + DATA.wipLimit + "（" + wip.map((task) => task.id).join("、") + "）");
  if (stale.length > 0) risks.push(stale.length + " 条超过 " + DATA.staleDays + " 天未更新：" + stale.map((task) => task.id).join("、"));
  if (overrun.length > 0) risks.push(overrun.length + " 条实际已超出估算：" + overrun.map((task) => task.id + "（" + task.actual + "/" + task.estimate + "）").join("、"));
  if (gateReds.length > 0) risks.push(gateReds.length + " 条门禁红（" + gateReds.map((gate) => gate.id).join("、") + "）；其中 " + gateReds.filter((gate) => gate.acked).length + " 条已被人接受");
  // 锚点自检：截断后 agent 会先信这三块，所以它们指向的东西必须存在
  if (!tasks.some((task) => task.id === DATA.nextAction.taskId)) {
    risks.push("nextAction 指向不存在的任务：" + DATA.nextAction.taskId);
  }
  const premature = wip.filter((task) =>
    (task.dependsOn ?? []).some((id) => {
      const dependency = tasks.find((candidate) => candidate.id === id);
      return dependency === undefined || dependency.status !== "completed";
    }),
  );
  if (premature.length > 0) risks.push(premature.length + " 条在依赖未完成时已开工：" + premature.map((task) => task.id).join("、"));
  // 悬空依赖是**比未完成更严重**的一类错误：依赖的任务根本不在板上，于是
  // 「依赖未完成」永远不会自己解除，卡会静默挂死。必须单独报，不能混进上面那条。
  const dangling = tasks
    .map((task) => ({ id: task.id, missing: (task.dependsOn ?? []).filter((d) => !tasks.some((t) => t.id === d)) }))
    .filter((row) => row.missing.length > 0);
  if (dangling.length > 0) {
    risks.push(
      "**悬空依赖（死锁）** " + dangling.length + " 条：" +
      dangling.map((row) => row.id + " 依赖不存在的 " + row.missing.join("、")).join("；") +
      " —— 依赖不在板上就永远不会解除，必须立卡或删掉该依赖",
    );
  }
  // 同族张力：两条在办的卡共用一个文件/路径时必须一起解，否则会互相打补丁
  const sharedPath = wip.filter((task) => wip.some((other) => other.id !== task.id && other.write === task.write));
  if (sharedPath.length > 0) {
    risks.push("在办共用同一改动位置（只能一起解，禁止各自打补丁）：" + sharedPath.map((task) => task.id + "→" + task.write).join("、"));
  }

  const nextUp = open
    .slice()
    .sort((a, b) => priorityRank(a.priority) - priorityRank(b.priority) || (a.updatedAt < b.updatedAt ? 1 : -1))
    .slice(0, 3);

  const visible =
    filter === "all" ? tasks
      : filter === "closed" ? closed
        : filter === "blocked" ? blocked
          : filter === "open" ? open
            : tasks.filter((task) => task.lane === filter);

  const byLane = DATA.lanes.map((lane) => {
    const rows = tasks.filter((task) => task.lane === lane);
    const laneCounted = rows.filter((task) => task.status !== "cancelled");
    const laneEst = laneCounted.reduce((sum, task) => sum + task.estimate, 0);
    const laneWeighted = laneCounted.reduce((sum, task) => sum + task.estimate * progressOf(task), 0);
    return {
      lane,
      total: rows.length,
      open: rows.filter(isOpen).length,
      blocked: rows.filter((task) => task.status === "blocked").length,
      done: rows.filter((task) => task.status === "completed").length,
      pct: laneEst === 0 ? 0 : Math.round(laneWeighted / laneEst),
    };
  });

  const active = tasks.find((task) => task.id === activeId) ?? nextUp[0] ?? tasks[0];
  const activeAge = active === undefined ? null : daysBetween(active.updatedAt, DATA.asOf);

  const startTurn = (task: TaskRow) => dispatch({
    type: "startTurn",
    prompt:
      "处理 " + task.id + "：" + task.title +
      "\n目标：" + task.goal +
      "\n验收：" + task.acceptance +
      "\n改动位置：" + task.write +
      "\n证据：" + task.evidence +
      "\nOCCT 参考：" + task.ref +
      "\n当前：" + STATUS_LABEL[task.status] + "，进度 " + progressOf(task) + "%，估算 " + task.estimate + " 人日 / 已用 " + task.actual + " 人日" +
      "\n阻塞：" + (task.blocker === "" ? "无" : task.blocker) +
      "\n前置：" + ((task.dependsOn ?? []).length === 0 ? "无" : (task.dependsOn ?? []).join("、") + (danglingIds(task, tasks).length > 0 ? "（⚠ 其中 " + danglingIds(task, tasks).join("、") + " 不在板上，属悬空依赖）" : "")) +
      "\n下一步：" + task.next +
      "\n遵守本板「约束与红线」：无同等 OCCT 分支就标 UNPORTED（写 文件+行号），不加自造谓词/阈值，不为对齐新写或改测试；" +
      "\n只改与本任务相关的代码；完成后更新这张画布的 DATA（进度、更新日期、证据）。",
  });

  const nextAction = DATA.nextAction;

  return (
    <Stack gap={24}>
      <Stack gap={10}>
        <H1>项目看板 · OCCT→Rust 移植</H1>
        <Text tone="secondary">{DATA.goal}</Text>
        <Text size="caption" tone="tertiary">
          快照 {DATA.asOf} · 修订 {DATA.revision} · 跟踪 {tasks.length} 条（在办 {open.length} / 关闭 {closed.length}）· 快照前 7 天完成 {recentDone.length} 条 · WIP 上限 {DATA.wipLimit}
        </Text>
        <Progress
          value={progressPct}
          showValue
          label={"加权进度（按估算人日，计入 " + counted.length + " 条未取消任务）"}
          tone={progressPct >= 85 ? "success" : progressPct >= 60 ? "info" : "warning"}
        />
        <Text size="caption" tone="tertiary">
          口径：跟踪的是「尚未对齐的缺口」，不是交付总量 —— 已对齐的历史条目见文末「已归档」，不参与加权。
        </Text>
      </Stack>

      {/* 概览卡只有「一个数字 + 一句标签」，允许并排；auto-fit 让窄栏自动折成 2 列 */}
      <Grid columns="repeat(auto-fit, minmax(120px, 1fr))" gap={12}>
        <Stat value={open.length} label="在办 / 待办" hint={"进行中 " + wip.length + " · 待办 " + tasks.filter((task) => task.status === "pending").length} tone={open.length > 0 ? "info" : "success"} />
        <Stat value={wip.length + "/" + DATA.wipLimit} label="在办 / WIP 上限" tone={wip.length > DATA.wipLimit ? "danger" : "info"} />
        <Stat value={blocked.length} label="阻塞" tone={blocked.length > 0 ? "danger" : "neutral"} hint={blocked.length > 0 ? blocked.map((task) => task.id).join("、") : "无"} />
        <Stat value={progressPct + "%"} label="加权进度" tone={progressPct >= 85 ? "success" : "warning"} hint={Math.round(weighted / 100) + " / " + estTotal + " 人日"} />
        <Stat value={gated.length} label="门禁需处理" tone={gated.length > 0 ? "danger" : "success"} hint={gated.length > 0 ? gated.map((gate) => gate.id).join("、") : "全部绿或已接受"} />
      </Grid>

      {risks.length === 0 ? (
        <Callout tone="success" title="没有需要立即处理的风险">
          <Text size="small">无阻塞、无超期未更新、WIP 未超限、无超估算、门禁无未处理红项。</Text>
        </Callout>
      ) : (
        <Callout tone={blocked.length > 0 || gated.length > 0 ? "danger" : "warning"} title={"派生风险 " + risks.length + " 项"}>
          <Stack gap={4}>
            {risks.map((risk) => <Text key={risk} size="small">{risk}</Text>)}
          </Stack>
        </Callout>
      )}

      {/* 全局唯一的下一个动作：截断后第一眼要看到的就是它 */}
      <Callout tone="info" title={"现在该做：" + nextAction.action}>
        <Text size="small">
          对应 <Code>{nextAction.taskId}</Code>：{nextAction.why}
        </Text>
      </Callout>

      <Row gap={8} wrap>
        <Pill active={filter === "open"} onClick={() => setFilter("open")}>在办 {open.length}</Pill>
        <Pill active={filter === "blocked"} onClick={() => setFilter("blocked")}>阻塞 {blocked.length}</Pill>
        <Pill active={filter === "closed"} onClick={() => setFilter("closed")}>已关闭 {closed.length}</Pill>
        <Pill active={filter === "all"} onClick={() => setFilter("all")}>全部 {tasks.length}</Pill>
        <Divider />
        {DATA.lanes.map((lane) => (
          <Pill key={lane} active={filter === lane} onClick={() => setFilter(lane)}>
            {lane} {tasks.filter((task) => task.lane === lane && isOpen(task)).length}
          </Pill>
        ))}
      </Row>

      <Stack gap={20}>
        <Stack gap={12}>
          <H2>待办</H2>
          {DATA.lanes.map((lane) => {
            const laneRows = visible.filter((task) => task.lane === lane);
            if (laneRows.length === 0) return null;
            const laneOpen = laneRows.filter(isOpen);
            const summary = byLane.find((item) => item.lane === lane);
            return (
              <CollapsibleSection
                key={lane}
                title={lane}
                count={laneRows.length}
                defaultOpen={laneOpen.length > 0}
                trailing={<Text size="caption" tone="tertiary">{laneOpen.length} 在办 · {summary === undefined ? 0 : summary.pct}%</Text>}
              >
                <Stack gap={8}>
                  <Progress value={summary === undefined ? 0 : summary.pct} size="sm" tone={laneRows.some((task) => task.status === "blocked") ? "danger" : LANE_TONE[lane]} />
                  <TodoList
                    todos={laneRows.map((task) => ({ id: task.id, status: task.status, content: task.priority + " · " + task.title }))}
                    onTodoClick={(todo) => setActiveId(todo.id)}
                  />
                </Stack>
              </CollapsibleSection>
            );
          })}
          <CollapsibleSection
            title="已归档（历史条目，不进当前窗口）"
            count={15}
            trailing={<Text size="caption" tone="tertiary">完整证据在 git 与 specs/</Text>}
          >
            <Stack gap={6}>
              <Text size="small" tone="secondary">
                同批完成的还有 T-80（BOP 根因三处，commit da61a30）、T-05（Offset 解析体积对齐 GT，2610.501436 vs 2610.501440，commit 59fefad）、T-04（phase19 断言按 OCCT 口径订正，GT 不相交 68 / 相交 80）、T-01（groove_cuts_cylinder 转绿，commit 113b4d1）。
              </Text>
              <Text size="small" tone="tertiary">
                这些条目在上一版画布里内联为卡片；按 0.2.0 口径只内联当前窗口，证据留在 git 提交与 specs/_a3n00_gap_analysis.md。
              </Text>
            </Stack>
          </CollapsibleSection>
        </Stack>

        {active ? (
          <Card>
            <CardHeader trailing={<Pill size="sm" tone={STATUS_TONE[active.status]}>{STATUS_LABEL[active.status]}</Pill>}>
              {active.id + " · " + active.priority + " · " + active.lane}
            </CardHeader>
            <CardBody>
              <Stack gap={12}>
                <Text weight="semibold">{active.title}</Text>
                <Progress
                  value={progressOf(active)}
                  size="sm"
                  showValue
                  label={"进度 · 已用 " + active.actual + " / 估算 " + active.estimate + " 人日"}
                  tone={STATUS_TONE[active.status] === "neutral" ? "info" : STATUS_TONE[active.status]}
                />
                <KeyValue
                  dense
                  items={[
                    { label: "负责人", value: active.owner === "unassigned" ? "未认领" : active.owner },
                    { label: "开始", value: active.startedAt === "" ? "—" : active.startedAt },
                    { label: "最近更新", value: active.updatedAt + (activeAge === null ? "" : "（" + activeAge + " 天前）") },
                    { label: "目标", value: active.goal },
                    { label: "验收", value: active.acceptance },
                    { label: "下一步", value: active.next, tone: "info" },
                    { label: "阻塞", value: active.blocker === "" ? "无" : active.blocker, tone: active.blocker === "" ? "neutral" : "danger" },
                    { label: "证据", value: <Code>{active.evidence}</Code> },
                    { label: "改动", value: <Code>{active.write}</Code> },
                    { label: "OCCT", value: <Code>{active.ref}</Code> },
                    { label: "前置", value: (active.dependsOn ?? []).length === 0 ? "无" : (active.dependsOn ?? []).join("、"), tone: "warning" },
                    { label: "备注", value: active.note },
                  ]}
                />
                <Divider />
                <Row gap={8} wrap>
                  <Button variant="primary" onClick={() => startTurn(active)}>交给 agent</Button>
                  <Button disabled={active.status === "in_progress"} onClick={() => overlay.set(active.id, { status: "in_progress" })}>标记进行中</Button>
                  <Button disabled={active.status === "completed"} onClick={() => overlay.set(active.id, { status: "completed" })}>标记完成</Button>
                  <Button disabled={active.status === "blocked"} onClick={() => overlay.set(active.id, { status: "blocked" })}>标记阻塞</Button>
                  <Button onClick={() => overlay.set(active.id, { status: "cancelled" })}>豁免</Button>
                  <Button variant="ghost" onClick={() => overlay.set(active.id, { owner: "human" })}>认领</Button>
                  <Button variant="ghost" onClick={() => overlay.clear(active.id)}>还原源数据</Button>
                  <Button variant="ghost" onClick={() => dispatch({ type: "openFile", path: active.write })}>打开文件</Button>
                </Row>
                <Text size="caption" tone="tertiary">
                  人的改动落进 sidecar（<Code>.canvas/board.state.json</Code>），agent 下一轮 <Code>canvas_read</Code> 就能看到。
                </Text>
              </Stack>
            </CardBody>
          </Card>
        ) : null}
      </Stack>

      <H2>下一步（按优先级取前 3 条）</H2>
      <Grid columns="repeat(auto-fit, minmax(240px, 1fr))" gap={12}>
        {nextUp.map((task) => (
          <Card key={task.id}>
            <CardHeader trailing={<Pill size="sm" tone={PRIORITY_TONE[task.priority]}>{task.priority}</Pill>}>
              {task.id + " · " + task.lane}
            </CardHeader>
            <CardBody>
              <Stack gap={8}>
                <Text size="small" weight="semibold">{task.title}</Text>
                <Text size="caption" tone="secondary">{task.next}</Text>
                <Row gap={6} wrap>
                  <Button size="sm" variant="primary" onClick={() => startTurn(task)}>开始</Button>
                  <Button size="sm" variant="ghost" onClick={() => setActiveId(task.id)}>看详情</Button>
                  {task.status === "blocked" ? <Pill size="sm" tone="danger">阻塞</Pill> : null}
                  {(task.dependsOn ?? []).length === 0 ? null : <Text size="caption" tone="warning">{"前置 " + (task.dependsOn ?? []).join("、")}</Text>}
                </Row>
              </Stack>
            </CardBody>
          </Card>
        ))}
      </Grid>

      <H2>分组成度（派生）</H2>
      <BarChart
        categories={byLane.map((item) => item.lane)}
        series={[{ name: "加权进度 %", data: byLane.map((item) => item.pct), tone: "info" }]}
        beginAtZero
        height={160}
      />
      <Table
        headers={["lane", "条目", "在办", "阻塞", "已完成", "进度%"]}
        columnAlign={["left", "right", "right", "right", "right", "right"]}
        striped
        rows={byLane.map((item) => [item.lane, item.total, item.open, item.blocked, item.done, item.pct + "%"])}
        rowTone={byLane.map((item) => (item.blocked > 0 ? "danger" : item.open === 0 ? "success" : "neutral"))}
      />

      <H2>门禁基线</H2>
      <Text size="caption" tone="tertiary">
        「基线 / 最近」逐字来自 DATA；点红项可标「已接受」（走 overlay，agent 可见）。本轮未重跑门禁 —— 见文末维护说明。
      </Text>
      <Table
        headers={["门禁", "命令", "基线", "最近", "Δ", "状态", "归属", "说明"]}
        columnAlign={["left", "left", "left", "left", "right", "left", "left", "left"]}
        striped
        stickyHeader
        rows={DATA.gates.map((gate) => {
          const acked = overlay.get(gate.id)?.acked === true;
          const ok = gate.status === "green";
          return [
            gate.name,
            <Code>{gate.cmd}</Code>,
            gate.baseline,
            gate.last,
            gate.baseline === gate.last ? "0" : "changed",
            <Pill
              size="sm"
              tone={acked ? "neutral" : ok ? "success" : "danger"}
              onClick={() => overlay.set(gate.id, { acked: true })}
            >
              {acked ? "已接受" : ok ? "pass" : "fail"}
            </Pill>,
            gate.owner === "" ? <Text size="small" tone="tertiary">—</Text> : <Code>{gate.owner}</Code>,
            <Text size="small" tone="tertiary">{gate.note}</Text>,
          ];
        })}
        rowTone={DATA.gates.map((gate) => (gate.status === "green" ? "success" : gate.acked ? "warning" : "danger"))}
      />

      <H2>明细</H2>
      <Table
        headers={["ID", "优先", "分组", "标题", "状态", "负责人", "进度", "估/实", "更新", "证据"]}
        columnAlign={["left", "left", "left", "left", "left", "left", "right", "right", "right", "left"]}
        striped
        stickyHeader
        emptyText="没有符合当前筛选的条目"
        onRowClick={(index) => {
          const row = visible[index];
          if (row !== undefined) setActiveId(row.id);
        }}
        rows={visible.map((task) => [
          <Code>{task.id}</Code>,
          <Pill size="sm" tone={PRIORITY_TONE[task.priority]}>{task.priority}</Pill>,
          task.lane,
          task.title,
          <Pill size="sm" tone={STATUS_TONE[task.status]}>{STATUS_LABEL[task.status]}</Pill>,
          task.owner === "unassigned" ? "未认领" : task.owner,
          <Progress value={progressOf(task)} size="sm" showValue tone={STATUS_TONE[task.status] === "neutral" ? "info" : STATUS_TONE[task.status]} />,
          task.actual + "/" + task.estimate,
          task.updatedAt,
          <Text size="caption" tone="tertiary">{task.evidence}</Text>,
        ])}
        rowTone={visible.map((task) => (task.status === "blocked" ? "danger" : task.status === "completed" ? "success" : task.status === "cancelled" ? "neutral" : "warning"))}
      />

      {/* 上下文锚点：截断或换人接手时先读这三节，避免重开已经定过的事 */}
      <CollapsibleSection title="约束与红线" count={DATA.constraints.length} defaultOpen>
        <Table
          headers={["规则", "为什么", "违反了会怎样"]}
          columnAlign={["left", "left", "left"]}
          rows={DATA.constraints.map((item) => [item.rule, item.because, item.violation])}
          emptyText="没有登记红线"
        />
      </CollapsibleSection>

      <CollapsibleSection
        title="已定决策（含被否方案）"
        count={DATA.decisions.length}
        trailing={<Text size="caption" tone="tertiary">防止截断后重新讨论</Text>}
      >
        <Table
          headers={["决策", "日期", "选了", "已否决", "理由", "参考"]}
          columnAlign={["left", "left", "left", "left", "left", "left"]}
          rows={DATA.decisions.map((item) => [
            <Code>{item.id}</Code>,
            item.at,
            item.chose,
            item.rejected.join("；"),
            item.why,
            <Code>{item.ref}</Code>,
          ])}
          emptyText="没有登记决策"
        />
      </CollapsibleSection>

      <CollapsibleSection
        title="活动"
        count={DATA.activity.length}
        defaultOpen
        trailing={<Text size="caption" tone="tertiary">只留最近 {DATA.activity.length} 条；更早的看提交历史</Text>}
      >
        <Timeline
          events={DATA.activity.map((event) => ({
            id: event.id,
            at: event.at,
            title: event.title,
            tone: event.tone,
            detail: event.detail,
            ref: event.ref,
          }))}
        />
      </CollapsibleSection>

      <CollapsibleSection title="口径与节奏（intake 决议）" count={7}>
        <Stack gap={6}>
          <Text size="small">· 决定：下一批做哪条 + 门禁能否提交 + 哪些红是已知非缺陷。</Text>
          <Text size="small">· 行单位：任务 <Code>T-xx</Code>（lane = 模块/族，共 {DATA.lanes.length} 组）；门禁是另一段，行单位不同。</Text>
          <Text size="small">· 状态：pending / in_progress / blocked / completed / cancelled；blocked 必须写明在等谁（本板只有 T-54）。</Text>
          <Text size="small">· 字段：状态 · 优先级 · owner · progress/estimate/actual · startedAt/updatedAt/completedAt · blocker · next · acceptance · evidence · write · ref · dependsOn。</Text>
          <Text size="small">· 派生（全部现算，不手抄）：加权进度 · WIP 超限 · 阻塞 · 陈旧 · 超估算 · 门禁待处理 · 依赖未完成即开工 · 在办共用改动位置 · 分组完成度 · 下一步 3 条。</Text>
          <Text size="small">· 人的动作：标记进行中/完成/阻塞、豁免、认领、还原、门禁标「已接受」、打开文件、交给 agent —— 全部走 overlay。</Text>
          <Text size="small">· 更新节奏：agent 每完成一批就更新 <Code>DATA</Code>（进度、更新日、证据）；人只动状态/认领/接受。</Text>
          <Text size="small" tone="tertiary">· 归档：completed / cancelled 超过一个窗口就移出当前窗口（见「已归档」折叠区），完整证据留在 git 与 specs/_a3n00_gap_analysis.md。</Text>
          <Text size="small" tone="tertiary">· 本轮门禁（2026-09-29，无代码改动，纯诊断）：五 crate check 0 error、--lib 1281/0、step_obj_gates 5/5 —— 与基线一致；诊断探针已全部回退，工作树只有本画布文件的改动。</Text>
        </Stack>
      </CollapsibleSection>

      <Text size="small" tone="tertiary">
        筛选器与当前选中项只存在本机（<Code>useCanvasState</Code>），agent 看不到；状态与门禁接受落 sidecar，agent 可见。明细表的「证据」列指向 文件:行 / 测试名 / 提交号，结论性描述一律用引用。
      </Text>
    </Stack>
  );
}
