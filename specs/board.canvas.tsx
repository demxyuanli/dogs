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
  asOf: "2026-09-30",
  revision: "r22",
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
  ],
  nextAction:   { taskId: "T-101",
    action: "**【§9.369 之后：a3n00 只剩 2 个面没做 seam 合并（端口 `{1:206,2:10,4:2}` vs OCCT `{1:208,2:9,4:1}`），且法兰孔面密度 72 vs GT 52】** ① 起点：§9.369 把导入期 `ShapeFix_Face::FixMissingSeam` 接回 reader，a3n00 面积比 **0.8627→0.8996**、T0M 未网格 **7→6**、UNPORTED UV-grid rescue **16→0**，带倒角法兰的 16 个面（10 圆柱 + 3 圆锥 + 3 BSpline 型）从 2 wires 变成 **1 wire / 5 条边** 并与 OCCT 模型孪生 f=85 逐项相同。② **第一件事（先做）**：按 bbox 把「端口 10 个 2-wire 面 + 2 个 4-wire 面」与「OCCT 9 个 2-wire + 1 个 4-wire」逐面配出来，锁死那 2 个 OCCT 合并而端口没合并的面，再对着 `ShapeFix_Face.cxx:1899-2330` 的分支看它们为何不进取合路径 —— **不许按面加特例，只按 .cxx 的分支走**。③ **第二件事**：`zz_seam_fix --all` 在 a3n00 上 `faces=226 changed=46`，其中 42 个是 `2→4 / 2→5` 边（seam 合并），另 4 个是 `13→15 / 13→16 / 22→28 / 5→8`；查这 4 个在 OCCT 里是否也改（对应 `FixReorder` `cxx:2029-2032` 与 post-seam 面循环那两处 UNPORTED）。④ **第三件事**：法兰孔/凸台面端口 **72** 个三角 vs GT **52**（rescue 已不参与）⇒ 在 `GenerateSurfaceNodes` / deflection 一侧按 `.cxx` 查密度差的来源，不要动阈值。⑤ 上一轮的仪器与数据可直接复用：探针 `--delaunstruct` / `--boundary`、`.target-gate/delaun_in/`、`.target-gate/pair368b.py`、`zz_uv_feed --ids` / `--model <f>`、`zz_seam_fix <file> <f> [--all]`。门禁：a3n00 面积比从 **0.8996** 上升且不劣化、T0M 未网格 **6** 不增加、`--lib` 1281/0、`step_obj_gates` 5/5。",
    why: "用户要求「a3n00 一个一个处理，先处理带有倒角的法兰盘」，本轮就把那一步做完并达成 T-99 的 accept。① 做法是把 f1e56776 删掉的导入期 seam 步骤接回 reader（`ShapeFixFace::fix_missing_seam` + 结果面 wire 的 `check_pcurves_and_shift`），这不是新造规则：`ShapeProcess_OperLibrary.cxx:785-899` 的 FixShape 算子 → `ShapeFix_Face::Perform` → `FixMissingSeam`（`cxx:482-498`，构造在 `cxx:1722-2330`），`STEPControl_Controller.cxx:201/:221` 默认就开。② 移除理由（T-93(a)/§9.219「OCCT 在这条路径上不跑 Perform」）被本轮实测推翻：探针 `--faceids` 开/关 ShapeProcess 得 `{1:208,2:9,4:1,6:4,10:4}` vs `{1:163,2:53,4:2,6:4,10:4}`，后者与 STEP 的 FACE_BOUND 直方图、与端口逐字相同 ⇒ OCCT 的导入确实合并了 45 个 2-bound 面。③ 结果方向与量级都对：法兰面结构变成 OCCT 的 1 wire / 5 边（边序、参数区间逐项相同），a3n00 面积比 0.8627→0.8996、T0M 未网格 7→6、T0M 面积比 0.9786→0.9987、acs10 0.9025→0.9846、rescue 16→0，且**没有改任何断言或 area_tol**。④ 剩余问题都已量化成有界的下一步（还差 1 个 2-bound + 1 个 4-bound 面没合并；全形状多改 4 个面待查；法兰面密度 72 vs GT 52），所以 T-99 按 accept 结项，另立 T-100 按「面」继续推进，而不是再调全局参数。",
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
      updatedAt: "2026-09-30",
      completedAt: "2026-09-30",
      blocker: "",
      goal: "把 a3n00 的面积比从当前真实值 **0.8627** 推向 1（缺口 31179.6）。**注意**：F113 在当前 HEAD **已经出网格**（`wires=4 mt=267`），所以「让 F113 出网格」不再是本卡目标；旧表述里的 0.8996 / 缺口 22802 / F113 不出网格来自本会话丢失的未提交改动，已按 §9.324 更正",
      next: "**【§9.368 把端口的完整结构（72 个结构节点 + 72 条约束链 + 容差/格/顶点序）喂给 OCCT 的 `BRepMesh_Delaun`：18/18 逐字段完全相同 ⇒ Delaunay 层忠实；病灶在输入端的面拓扑】** ① 做法：端口侧 TEMP 落盘 `delaun_in_f<f>.txt`（`TOL/CELLS/N/L/V` + `AFTER_CTOR`/`FINAL`）与 `delaun_uv_f<f>.txt`（登记序原始 UV）；探针新增 `--delaunstruct <dir>` —— `AddNode` 逐个核对返回值等于声明序号、`AddLink`×72、`SetTolerance/SetCellSize`，再用 public 构造 `BRepMesh_Delaun(aStruct, indices, cellsU, cellsV)` 跑完整条路，**绕开** §9.366 的 LNK2019。② **结果**：18 个面（含 f=0、f=208 两个健康对照）OCCT 的 `NbNodes/NbLinks/ElementsOfDomain` 与四种链状态计数**与端口逐字段相同**（`idx_mismatch=0`、`links_dumped == links_in_struct == 72`；f=0 是 `nodes=71 links=207 domain=66 frontier=68 free=65 deleted=74`，f=192 是 `nodes=75 links=245 domain=0 frontier=72 deleted=173`）⇒ **这一层的实现是忠实的**，§9.359 的「同一批点、顺序不同就归零」不是这层的代码差异。③ **更正 §9.367**：只喂点（0 条链）时的 `domain=0` 是 `cleanupMesh` 造成的 —— `BRepMesh_Delaun.cxx:1028` 无条件调用它，而它遍历 `FreeEdges()` 时只跳过 Frontier 链（`cxx:832-835`），结构里一条 Frontier 都没有时所有三角都被判为外部删掉（`cxx:908-911`）⇒ 那条的「能推出：这 74 个点本身不足以三角化」**作废**。④ **更正 §9.365.2**：a3n00 两侧面序号**只有 3/226 对应**；按六坐标 bbox 配对（201/226 精确 1:1）后 port 192 的真孪生是 model 85、port 169 是 model 87、port 174 是 model 42……那 16 个面的 GT 三角数是 **52–80**（不是 2–15），端口是 204–336 ⇒ 真实倍数约 **4 倍**且**方向双向**（189/204 是端口少铺：36 vs 54）。⑤ **真正的输入差**：这 16 个面在 STEP 里声明 **2 个 `FACE_BOUND`**，端口忠实建成 2 条 wire；OCCT 读入期 `ShapeProcess` 把它们并成 **1 条 wire**（`--faceids` 直方图：默认 `{1:208,2:9,4:1,6:4,10:4}`，`--nofix` `{1:163,2:53,4:2,6:4,10:4}` 与 STEP、与端口**逐字相同**）⇒ 整形合并了 45 个 2-bound 面。逐边对照 port 192 ↔ model 85：OCCT 是 **1 条 wire / 5 条边 / 7 条 pcurve / 83 点**，由**两条 2 点 seam 边**把上下两圆连成一个环、圆 A 被切成 19+19；端口是 **2 条互不相连的闭合 pcurve**（v=-267.805924 与 v=-287.805924 各 37 点）。同时更正 §9.364/§9.365.6：`UVSUM` 的 `wires` 就是拓扑 wire 数，**不是第三个量**。⑥ **下一轮**：实现读入期把 2 个 bound 并成 1 条 seam 闭合 wire 的那段控制流（候选位置 `crates/occt-topo/src/shhealing/shape_fix_face.rs:147` 的 `w2 != null` 合并、`:406` 的 `FixReorder`；对照 `ShapeFix_Face.cxx:492-498` `:1722-2330`）。**可证伪的预测**：接上后 port 192 的结构应变成 1 wire / 83 点 / 含两条 2 点 seam pcurve，`AFTER_CTOR domain` 由 0 变 >0、`mt` 由 228 降到约 52，a3n00 面积比从 0.8627 上升。⑦ 本轮改动：**Rust 零改动**（插桩已按 §9.353 规程用 `edit` 反向撤除，`crates/` diff 为空），探针新增 `--delaunstruct` 与 `--boundary` 保留作仪器。门禁：`cargo check` 0 error、`--lib` 1281 passed / 0 failed、a3n00 `stats=226 unmatched=0 mesh=11101/11121`。",
      acceptance: "a3n00 面积比从 **0.8627** 上升且不劣化；T0M 未网格 7 不增加（`f1e56776` 已记 T0M 6→7 的代价）；--lib 与 step_obj_gates 5/5 不劣化；**不得改基线或 area_tol 来绕过**；不新写测试、不改断言",
      evidence: "**2026-09-30，十三轮**。**(1) 面已 100% 锁定**（§9.303.1）：STEP `#5375` = `CYLINDRICAL_SURFACE` **R=34**、高 68，解析面积 `2πrh = 29053.4` 与该区域缺口逐字吻合。**(2) 下游逐级实测忠实**：reader `FACE_BOUND 到 wire` **1:1**；读入路径**无 healing 阶段**；4 条 wire **不共享顶点 TShape**；4 条 wire 在 3D 上**各自闭合**（gap 精确 0）⇒ `CheckWire` 拒绝**正确**；范围守卫逐行忠实 `cxx:1781-1802`；`ComposeShell` 的 `load_wires` 入口即 `[6,14,1,1]`。**(3) 病根由实验确认**（§9.315）：`shape_fix_face.rs:453` 的 `surf.clone()` 让 `tmp_f` 与原面共用 `repr_key`，seam 写入落到原面槽位；摘掉 `tmp_f` 的曲面注册后 ⇒ `unmatched` **1 到 0**、F113 出网格 `mt=267`、全 226 面出网格，但面积比 **0.8996 到 0.8627** ⇒ 不能落地。**(4) 面积差成因（§9.316 + §9.319）**：逐三角对比显示两版**最大的三角完全同名同值**（不存在多出一张巨大错误面）；实验版 `v` 更多却 `f` 更少、总面积少 8377 ⇒ **少铺**；逐面边数差集显示差异落在 **43 个面上、每个少 2–3 条 seam 边、合计 −107**。原因是 `fix_missing_seam` 要**读回自己刚写下的 pcurve**（`check_pcurves_and_shift`/`curve_on_surface_range` 都按 `repr_key` 取），摘注册后写与读不在同一个键上 ⇒ 读不到 ⇒ seam 边加不上。**(5) 核心改动试并否决（§9.320，本轮）**：把 `repr_key` 改为 `face_key` ⇒ 三个指标与摘注册实验**逐字相同**（`stats` 226、`unmatched` 0、`mesh` 11101/11121、面积比 0.862724、`edges` 976 对基线 1083、同样 43 面退化）⇒ **单纯改键的归属无效**。**(6) 已排除 11 条假设**：CheckWire 口径、`v_range`、ComposeShell 内部拆分、返回 Shell/compound、顶点身份/邻接口径、两个退化边同母线、`heal_shape` 接在读入路径、D18 的 n1 到 n2、该写回可中性解耦、「每面的每条边都要有 pcurve」、以及**让 pcurve 的键跟随面身份**。**(7) 量化**：F113 单面区域占 +26952 面积缺口，全模型净缺口仅 22802.58",
      write: "crates/occt-topo/src/shhealing/shape_fix_face.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.303 §9.314 §9.315 §9.316 §9.318 · crates/occt-topo/examples/zz_pcurve_key_probe.rs · crates/occt-topo/examples/zz_share_probe.rs",
      note: "【§9.321】13 条假设已否决。本轮新增两条：第 12 条「共享键的危害经 swap_seam 生效」被插桩否决（读取路径 0 次调用，并因此**推翻 §9.296.3 的判断**）；第 13 条「共享键的第二个作用是让 check_wire 找到 pcurve」也被否决（两版早退都是 0）。范围已收窄到：差异不在这两处，而 unmatched 1→0 提示 F113 在两版里是不同的面对象 ⇒ 下一步量「哪些面被采纳与否」的差集，而不是再猜哪个函数",
      dependsOn: [],
    },
        { id: "T-59", lane: "mesh", title: "a3n00 面积缺口：配对表已出，缺口定到单面 F113（已拆为 T-99）",
       status: "completed", priority: "P1", owner: "agent", progress: 100, estimate: 5, actual: 2, startedAt: "2026-09-30",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
      goal: "按 bbox 六坐标把 a3n00 的 226 个端口面与 GT 面配对，逐面定位构成面积亏损的 face class",
      
      next: "本卡收口：配对表已出（226/226 全配上），缺口已定到**单个面 F113**，后续修复拆到 T-99。留下两件已量化的事实：① D27 seam 接线后 a3n00 比值 **0.8996**（判据 `|ratio-1| < area_tol=0.15`），比本卡起点 0.8627 高 0.037 ⇒ **§9.287.4 那行修复（当时 −0.0329）已不再是「应用后越界」的阻塞项**，但它本身仍是净负面（见 T-99 (d) 重测：−0.0092、且救不活 F113）⇒ **不要**再把它当主线。② 面类分布（确定配对的 face 数）：Plane 86 / Cylinder 62 / Cone 30 / Extrusion 15 / Torus 7；端口只有 1 面未网格（F113）。工具：`.target-gate/pair_a3n00.py`（配对）、`.target-gate/deficit_a3n00.py`（按 bbox 归因面积）",
       acceptance: "产出可复跑的配对表并定位亏损 face class（不要求修好）；不新写测试、不改断言",
      
      evidence: "**2026-09-30 实测**：`.target-gate/pair_a3n00.py .target-gate/a3n00_port.txt .target-gate/a3n00_gt.txt` ⇒ `port faces=226 gt faces=226`、`paired=226 port unpaired=0 gt unpaired=0`、`pairs with bbox diff > 1e-3: 26`。**wire 数不匹配只有 3 组**：`port wires=4 vs gt wires=2`（1 面，F113）、`port wires=2 vs gt wires=1`（2 面：F140→gt f=39 Cone、F170→gt f=69 Cylinder）⇒ a3n00 的 wire 结构基本已对齐，缺陷是**外科式**的而非系统性。**面积**：ours 204327.72 / occ 227130.30 = **0.8996**（净缺口 22802.58）；按 bbox 把缺口归因到各面后，**F113 单面所在区域就占 +26952.18**（`pad=0.01`；pad 0→2 变动 <2%）⇒ 超过全部净缺口，其余面合计是净盈余。**起点数字（2026-09-29，已被本行取代）**：本卡原记 a3n00 比值 0.8627、应用 §9.287.4 后 0.8298（−0.0329）；D27 seam 接线后这两个数字都已过期，**以 0.8996 为准**。宽 tol 与面类的登记出处：`crates/occt-topo/tests/common/mod.rs:91-93`",
       write: "", ref: "specs/_a3n00_gap_analysis.md §9.291 · §9.287 · tests/common/mod.rs:93 · .target-gate/pair_a3n00.py",
       note: "本卡原为 T-93 的前置（D20 补立）。分析类卡，**库源码零改动**；发现的 F113 缺陷已拆为 T-99，避免把两类缺陷混在一张卡里",
      
      dependsOn: [],
    },
        { id: "T-93", lane: "port-gap", title: "reader 缺带 Context 的 ShapeFix_Face::Perform（导入期 seam 修复）",
       status: "completed", priority: "P0", owner: "agent", progress: 100, estimate: 4, actual: 6, startedAt: "2026-09-28",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
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
       status: "completed", priority: "P1", owner: "agent", progress: 100, estimate: 3, actual: 5, startedAt: "2026-09-24",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
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
       status: "completed", priority: "P2", owner: "agent", progress: 100, estimate: 4, actual: 2, startedAt: "2026-09-20",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
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
       status: "completed", priority: "P1", owner: "agent", progress: 100, estimate: 8, actual: 2, startedAt: "2026-09-28",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
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
       status: "completed", priority: "P2", owner: "agent", progress: 100, estimate: 4, actual: 1, startedAt: "2026-09-30",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
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
       status: "completed", priority: "P2", owner: "agent", progress: 100, estimate: 4, actual: 1, startedAt: "2026-09-30",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
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
       next: "—", acceptance: "--lib 1281/0、step_obj_parity 14/14、export 逐位一致",
       evidence: "a3n00 wire 直方图 2-wire 面 24→10；mesh 11675/13168",
       write: "crates/occt-topo/src/shape_fix_compose_shell",
       ref: "ShapeFix_ComposeShell.cxx · ShapeFix_Face.cxx:1722-2330",
       note: "主路径上仍挂 T-95 的 UNPORTED",
      dependsOn: [],
    },
        { id: "T-82", lane: "bop", title: "BOPAlgo_BOP::BuildSolid 曲面布尔精确化", status: "completed", priority: "P0", owner: "agent",
      progress: 100, estimate: 3, actual: 3, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "曲面布尔走精确 BOPAlgo 路径而非采样近似", next: "—",
      acceptance: "--lib 不劣化；box∪cyl 单实体、体积≈8.5；phase10 / boss 不劣化",
      
      evidence: "HEAD 08c4edf2 实测 box∪cyl FUSE = 1 solid / 8 面 / volume 8.502655，逐面面积与 GT（occt_probe --fuse）逐一相同；卡片的 8.737 与 loops=0 是更早快照",
       write: "crates/occt-topo/src/bop_bop.rs", ref: "BOPAlgo_BOP.cxx:583-711 · BOPAlgo_Builder_3.cxx::BuildSplitSolids",
       note: "AC 未覆盖的其他曲面布尔组合未复核",
      dependsOn: [],
    },
        { id: "T-88", lane: "bop", title: "phase10 curved_face_fillet_sphere_plane", status: "completed", priority: "P0", owner: "agent",
      progress: 100, estimate: 1, actual: 1, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "phase10 的球面/平面倒圆用例转绿", next: "—", acceptance: "--test phase10_integration 8/8",
      evidence: "HEAD 08c4edf2 实测 phase10 8/8；旧板记的 7/8 是 2026-09-21 快照，原诊断未再复现",
       write: "crates/occt-topo/src/bop_bop.rs", ref: "BOPAlgo_BuilderSolid · BOPTools_AlgoTools", note: "—",
      dependsOn: [],
    },
        { id: "T-41", lane: "bop", title: "摘除 bop_curved 的体素/网格布尔（OCCT 无对应）",
       status: "completed", priority: "P1", owner: "agent", progress: 100, estimate: 2, actual: 3, startedAt: "2026-09-20",
      updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "", goal: "布尔入口改派到忠实路径，去掉 OCCT 里不存在的体素/采样布尔",
       next: "—", acceptance: "摘除后 --lib 与门禁不劣化",
      evidence: "改派后实测 --lib 1281/0、step_obj_gates 5/5、phase10/19/3/4 = 8/5/4/9；调用图 curved_boolean_full ← boolean_dispatch:227 / bop_builder_repair:239,405 / bop_builder_report:388，按 (a) 夹具重设计",
       write: "crates/occt-topo/src/bop_curved/region_trim.rs",
       ref: "—", note: "classify_face_general 本身仍是 8×8 采样，不是 BRepClass3d_SolidClassifier",
      
      dependsOn: [],
    },
        { id: "T-67", lane: "port-gap", title: "Extrema_ExtPExtS / ExtPRevS 两臂", status: "completed", priority: "P2", owner: "agent",
      progress: 100, estimate: 3, actual: 3, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "补上点-曲面极值搜索的两条 OCCT 分支，摘掉替换实现",
       next: "—", acceptance: "occt-geom --lib 143/0 + 门禁不劣化",
      evidence: "两臂已移植并分派（point_surface_extrema.rs:60-61/:88-91、:450 perform_ext_ps、:485 perform_rev_ps，对应 Extrema_ExtPS.cxx:292-343）；A15 替身已摘除。实测 occt-geom 143/0、occt-topo 1281/0、step_obj_gates 5/5",
       write: "crates/occt-geom/src/extrema_surf", ref: "Extrema_ExtPS.cxx:292-343", note: "—",
      dependsOn: [],
    },
        { id: "T-25", lane: "arch", title: "GeometryRegistry 侧表 → 几何进 TShape", status: "completed", priority: "P2", owner: "agent",
      progress: 100, estimate: 6, actual: 6, startedAt: "2026-09-20", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "几何不再走全局侧表，而是挂在 TShape 上", next: "—", acceptance: "--lib + 门禁不劣化",
      evidence: "三张几何表已删除（tgeometry.rs:163-181 注释写明 all live on their own TShape now），35 处读写转发到 TShape 槽；设计文档第 5 步（global() 兼容外壳）未做，属文档明示可停点。实测 occt-topo 1281/0、step_obj_gates 5/5",
       write: "specs/_design_architecture_t25_t28.md", ref: "specs/_design_architecture_t25_t28.md",
      note: "保留 global() 兼容外壳是既定停点，不是未完成项",
      
      dependsOn: [],
    },
        { id: "T-28", lane: "arch", title: "ImpPrm HVertex 合并 + intana/intpatch 重叠合并",
       status: "completed", priority: "P2", owner: "agent", progress: 100, estimate: 8, actual: 8, startedAt: "2026-09-20",
      updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "", goal: "ImpPrm 的 HVertex 合并与 intana 的 closed form 单一化",
       next: "—", acceptance: "--lib + 门禁不劣化",
      evidence: "步 1–3（ImpPrm HVertex 合并）与步 4（intana/intpatch closed form 合并，五副本清零）均完成；步 5 经实测判为有害不做（D5）。实测 occt-topo 1281/0、step_obj_gates 5/5",
       write: "specs/_design_architecture_t25_t28.md", ref: "IntPatch_ImpPrmIntersection.cxx:221-469",
      note: "若日后重开步 5，前提是先把 PatchIntersection 对「一般 B 样条 × 球」修正到与 tracer 等价",
      
      dependsOn: [],
    },
        { id: "T-29", lane: "hygiene", title: "过期规格刷新", status: "completed", priority: "P3", owner: "agent", progress: 100, estimate: 1,
      actual: 1, startedAt: "2026-09-21", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "规格文件不再指向已删除的 _board.md，且过期内容有横幅",
       next: "—", acceptance: "—",
      evidence: "两个目标文件已带「仅历史（T-29 登记）」横幅（_coverage.md:1、_brepmesh_align_review.md:3）；指向已改为 specs/board.canvas.tsx，并同步 4 处悬空引用。occt-core check 通过、occt-topo 1281/0",
       write: "specs/_coverage.md", ref: "—", note: "—",
      dependsOn: [],
    },
        { id: "T-11", lane: "hygiene", title: "occt-topo 编译警告清理", status: "cancelled", priority: "P3", owner: "agent", progress: 0,
      estimate: 1, actual: 4, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "", blocker: "",
      goal: "清理 unused import 类警告", next: "已豁免：剩余项经实验判定不可安全清除；死代码部分改由 T-97 承接",
       acceptance: "无行为变化",
      evidence: "四批共清 11 处真实 unused import，每批后 --lib 1281/0；剩余 ~60 条多为「只被同文件 #[cfg(test)] 使用」的假阳性，实测删除后 test 目标报 cannot find type ⇒ 已回退；两条自动化路径亦实测失败",
       write: "crates/occt-topo/src", ref: "—", note: "豁免结论见 D7",
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
       status: "cancelled", priority: "P3", owner: "human", progress: 0, estimate: 1, actual: 1, startedAt: "2026-09-19",
      updatedAt: "2026-09-21", completedAt: "", blocker: "", goal: "该夹具是否应通过",
      next: "维持 1/2 作为已知，不再投入（夹具不合法，非引擎缺口）",
       acceptance: "—", evidence: "夹具期望的「1 solid」在 OCCT 的 BOPAlgo 下也不成立；门禁 G-boss 记为已接受的红",
       write: "crates/occt-topo/tests/bop_builder2_boss.rs",
       ref: "BOPAlgo_BOP.cxx:583-711", note: "改动已按 (a) 夹具重设计落地（见 T-41）",
      dependsOn: [],
    },
        { id: "T-23", lane: "hygiene", title: "ATU01038 顶点/面密度 −1.8%", status: "cancelled", priority: "P3", owner: "human", progress: 0,
      estimate: 1, actual: 0, startedAt: "2026-09-21", updatedAt: "2026-09-21", completedAt: "", blocker: "", goal: "密度差异是否要追",
      next: "不改，作为已知记录", acceptance: "—", evidence: "UV-grid 与 deflection-adaptive 的采样差异；parity 已明确不断言密度",
       write: "crates/occt-topo/tests/common/mod.rs", ref: "—", note: "—",
      dependsOn: [],
    },
        { id: "T-101",
      lane: "mesh",
      title: "带倒角的法兰盘已正确处理；剩下 2 个面未合并 + 法兰面密度 72 vs GT 52",
      status: "pending",
      priority: "P1",
      owner: "agent",
      progress: 25,
      estimate: 8,
      actual: 1,
      startedAt: "2026-09-30",
      updatedAt: "2026-09-30",
      completedAt: "",
      blocker: "",
      goal: "按「面」继续推进 a3n00 到面积比 1：先把剩下的 seam 合并缺口补齐（端口 {1:206,2:10,4:2} vs OCCT {1:208,2:9,4:1}），再查法兰面密度差（端口 72 vs GT 52）",
      next: "**【§9.370 迭代（1）：只剩 3 个面的 seam 缺口（113/140/170），坏在「结果面 vs 结果 shell」】** ① 已定：这 3 个面 `fix_missing_seam` 返回 true 但结果是 **Shell**（5 / 2 / 2 个面），而**合并本身与 OCCT 逐一相同**（后置边数 113: 28=28、140: 16=16、170: 8=8，口径用 `W k edges=N` 不用 `E` 行）⇒ 分歧只能落在 `crates/occt-topo/src/shhealing/shape_fix_face.rs:590-673` 这段尾巴：假想 grid（`cxx:2236-2245`）→ `ComposeShell`（`cxx:2246-2261`）→ `myResult = CompShell.Result()`（`cxx:2263-2268`）→ 两轮剪枝（`FixSmall` / `FixSmallAreaWire`，`cxx:2270-2322`；判据 `crate::shhealing::check_small_area` ↔ `ShapeAnalysis_Wire::CheckSmallArea` `cxx:2004`）。reader 现在只接受 `ShapeType::Face`，所以这 3 个面保持 4/2/2 wires。**修的方向是让那段尾巴给出 OCCT 的单面结果** —— 不是改 seam 位置选择、更不是按下标/bbox 加特例；该段的分歧定位已交给一个并行子任务。② 同批的 138 已核：`--fixms` 对它返回 false，不属 Shell 类（`--all` 里它的 13→15 与 `zz_seam_fix` 的单面口径不同源）。③ 密度（法兰孔面端口约 72 vs GT 52）**先要解决测量**：`zz_probe_a3n00 --fstats` 与 `zz_uv_feed --ids` 的逐面三角数在 **89/226** 个面上不同（而逐面 bbox 226/226 相同）⇒ 逐面读数必须在**同一次运行内**取；下一轮先把逐面统计改成与模型序同源，再逐面比密度。④ 复跑命令见 `specs/_a3n00_gap_analysis.md` §9.370.4。",
      acceptance: "a3n00 面积比从 **0.8996** 上升且不劣化（step_obj_gates 的 step_obj_area 实测口径）；T0M 未网格 **6** 不增加；--lib 1281/0 与 step_obj_gates 5/5 不劣化；不得改基线或 area_tol；不新写测试、不改断言",
      write: "crates/occt-topo/src/shhealing/shape_fix_face.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.368 §9.369 · ShapeFix_Face.cxx:1722-2330",
      dependsOn: [],
      evidence: "§9.370：`zz_probe_a3n00 --fixms` 实测只有 113/140/170 返回 Shell（5/2/2 个面）；`--all` 里“多改”的 4 个面就是 113/138/140/170，且后置边数与 OCCT 孪生逐一相同（113:28=28、140:16=16、170:8=8，`W k edges=N` 口径）",
      note: "① 的尾巴分歧已交给并行子任务对照 ShapeFix_Face.cxx + shape_fix_compose_shell 深挖；不要在 reader 里按面号/bbox 特例地接受 Shell；③ 的逐面统计口径已按 §9.371 修正（FaceMeshStat.index）；§9.372 把密度差按曲面类型归并后 Torus 2.14× 最突出，§9.373 **自我否定**：ParamSet 升序**不是**分歧（OCCT 的 FUN_CalcAverageDUV 形参是非 const 引用、会就地排序 aParamArray，两侧抽稀都在升序上做）⇒ 不要动 ParamSet；③ 的正确下一步是先修配对（199 个配对面的 GT 三角和只有 6392/12324，大面系统性失配），再做 Torus 面的九量对拍（range/delta/r/R/ArcAngularStep/oldDv/nbV/Du/nbU/节点数）",
    },
  ],
  gates: [
        { id: "G-check", name: "五 crate 编译", cmd: "cargo check --all-targets", baseline: "exit 0", last: "exit 0", status: "green",
      acked: false, owner: "",
      note: "五 crate 逐 crate 跑（无根 workspace）：core · math · geom · geom2d · topo 各 0 error（2026-09-29 本轮实测）",
    },
        { id: "G-lib", name: "topo 单测", cmd: "--lib", baseline: "1281 / 0", last: "1281 / 0", status: "green", acked: false, owner: "",
      note: "2026-09-30 本轮实测 1281 passed / 0 failed（12.91s，T-93 seam 修复接线后）；摘除临时计数器后复跑仍 1281/0",
    },
        { id: "G-gates", name: "STEP→OBJ 门禁", cmd: "--test step_obj_gates", baseline: "5 / 5", last: "5 / 5", status: "green", acked: false,
      owner: "",
      note: "2026-09-30 本轮实测 5 passed（588.92s 接线后 / 704.02s 清理后复跑），**T-93 seam 修复接线未越出任何 area_tol**。历史：T-93 (a) 后两处基线重定（tests/common/mod.rs 写明来源）：occ/T0M.stp 面积比 0.9995→0.9786 ⇒ area_tol 0.01→0.025、occ/acs10.stp 0.9846→0.9025 ⇒ 0.05→0.12",
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
      note: "内部回归哨兵（旧板记 16/16 为过期）；保真度见 specs/_occt_mesh_gt.md",
    },
        { id: "G-iges", name: "IGES 结构自洽", cmd: "--example iges_check -- <18 模型>", baseline: "18 / 18", last: "18 / 18", status: "green",
      acked: false, owner: "",
      note: "15 个 data/*.step + data/occ/{ATU01038,bottom,top}.step；全部 unreferenced=1（只剩根）",
    },
        { id: "G-core", name: "core / math / geom / geom2d", cmd: "--lib", baseline: "290·215·143·72", last: "290·215·143·72", status: "green",
      acked: false, owner: "",
      note: "逐项与基线相同",
    },
  ],
  // 新条目插到数组**开头**（brief 取前 3 条当最近活动；早于 a29 的看提交历史）
  activity: [
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
       tone: "info", detail: "代价量化：T0M 未网格面 6→7、wires 1931→2121（§9.269）",
      
      ref: "crates/occt-topo/src/shhealing/shape_fix_face.rs",
    },
        { id: "a2", at: "2026-09-28", title: "T-69 矛盾解开：§9.265 的「190 个面丢 wire」是补偿副作用",
       tone: "info", detail: "补偿即使丢弃 Shell(5) 也改写共享 GeometryRegistry，使解析期探针与最终模型分叉",
      
      ref: "specs/_a3n00_gap_analysis.md §9.268",
    },
        { id: "a3", at: "2026-09-28", title: "§9.270 给出有界路径：只加回 pcurve 归一化（不含 fix_missing_seam 本体）",
       tone: "success", detail: "把 T-93/T-69 的张力从「需决策」降为「一次实测」",
      
      ref: "specs/_a3n00_gap_analysis.md §9.270",
    },
        { id: "a4", at: "2026-09-28", title: "T-93 根因定案：亏损来自单面 F171 的第 4 条 wire 少一条边",
       tone: "warning", detail: "端口用的端点 P 不在 STEP 里，且是边 #7405 曲线的中点 ⇒ reader 另建顶点",
      
      ref: "specs/_a3n00_gap_analysis.md §9.276–§9.279",
    },
        { id: "a5", at: "2026-09-28", title: "T-41 (a) 完成：repair/report 改派 + 夹具按合法布尔重设计",
       tone: "success", detail: "改派后 --lib 1281/0、step_obj_gates 5/5、phase10/19/3/4 = 8/5/4/9",
      
      ref: "crates/occt-topo/tests/bop_builder2_boss.rs",
    },
        { id: "a6", at: "2026-09-28", title: "T-28 结项（步 1–4 完成），步 5 经实测判为有害不做",
       tone: "success", detail: "统一入口 PatchIntersection 会让 general_bspline_sphere 曲线离球面 ≈0.67",
      
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
