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
  goal: "把 STEP→OBJ 几何/网格管线对齐 OCCT 8.0.0（源码树 D:\\source\\OCCT-src @ V…",
  asOf: "2026-10-07",
  revision: "r40",
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
       why: "两条老路互斥：恢复补偿则 F113 又不出网格，保留现状则 T0M wire 结构更远离 GT…",
      
      ref: "specs/_a3n00_gap_analysis.md §9.268 / §9.269",
    },
        { id: "D2", at: "2026-09-28", chose: "移除 reader 侧 ShapeFix_Face::FixMissingSeam 补偿（T-93 (a) 落地）",
       rejected: ["继续用补偿换取 T0M 的 wire 数"],
      why: "补偿即使丢弃 Shell(5) 也会改写共享 GeometryRegistry，使解析期探针与最终模型分叉，取证链被污染",
      
      ref: "specs/_a3n00_gap_analysis.md §9.266 / §9.268",
    },
        { id: "D3", at: "2026-09-28", chose: "boolean_dispatch 改派忠实路径；摘除 bop_curved 体素/网格布尔",
       rejected: ["保留采样计票 classify_face_general 作为曲面分类"],
       why: "体素/采样计票在 OCCT 无对应物（OCCT 走 BRepClass3d_SolidClassifier / IntT…",
      
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
       why: "两条自动路径实测失败（cargo fix 删掉只被 #[cfg(test)] 使用的 import 后报 109 个构建…",
      
      ref: "crates/occt-topo/src · §9.249",
    },
        { id: "D8", at: "2026-09-28", chose: "T-51 豁免：不重做 GeomTrimmedCurve 的 [0,1] 参数表示",
       rejected: ["为移植 gcpnts 两条 Compute arm 先改 TrimmedCurve 参数表示"],
       why: "属架构级改动且无消费方触发那两条 arm；改用 T-96 单独立卡评估",
      
      ref: "crates/occt-geom/src/gcpnts.rs:87-158 · GeomAdaptor_Curve.cxx:239-255",
    },
        { id: "D9", at: "2026-09-29", chose: "§9.270 的假设否证：不再试图用 pcurve 归一化 pass 挽回 T0M 未网格面",
       rejected: ["把被删块的 check_pcurves_and_shift 第二次 pass 加回 reader"],
      
      why: "实测加回后未网格面仍为 7、WIREHIST 逐字不变（1:1463 2:294 3:9 5:1 6:1 7:2 8:1…",
      
      ref: "crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs:441-462 · node_insertion.rs:536",
    },
        { id: "D10", at: "2026-09-29", chose: "按 GT 实测修正「已知收敛差异」的口径：端口丢的 6 个面在 OCCT 都出网格，属真分歧",
       rejected: ["把 T0M 未网格面记成 deflection/启发式导致的容忍差异"],
      
      why: "occt_probe --mesh 0.1 实测：GT 1778 面里 triangles=0 的只有 FACE 132…",
      
      ref: "specs/occt_probe/occt_probe.exe --mesh 0.1 data/occ/T0M.stp · 撤回见 D12/D13",
    },
        { id: "D11", at: "2026-09-29", chose: "T0M 6 个丢面按两个独立缺陷处理，不做「放宽容差/阈值」式合并修法",
       rejected: ["把 IsValid 判据放宽（改 PConfusion 或加最小宽度阈值）让这些面通过"],
      
      why: "端口 PCONFUSION=1e-9 与 OCCT Precision::PConfusion() 完全一致…",
      
      ref: "crates/occt-core/src/precision.rs:12 · BRepMesh_NodeInsertionMeshAlgo.hxx:80-84 · Precision.hxx:334",
    },
        { id: "D12", at: "2026-09-29",
      chose: "撤回 D10：T0M 的 6 个未网格面**不能**判定为「端口丢面」…",
       rejected: ["按 --mesh 的面序号逐面断言端口丢面", "把「7 vs 1」当成已证的缺陷差"],
      
      why: "① 单面复现证明这些 pcurve 几何正确且参数对齐（F1691 误差 7.1e-15 且 pcRange==edge…",
      
      ref: "BRepMesh_DefaultRangeSplitter.cxx:23-88 :202-236 · Precision.hxx:334 · crates/occt-topo/src/meshing/range_splitter/param_set.rs:787-838",
    },
        { id: "D13", at: "2026-09-29", chose: "「1778 vs 1772」的面数差与「7 vs 1」的未网格差解耦：不再用逐面序号对拍这两组量",
       rejected: ["用 occt_probe --mesh 的面序号直接对应端口 F 序号来判丢面"],
      
      why: "GT 面 1327（唯一 triangles=0 者）bbox=(27.75947,-172.18146,-34.051…",
      
      ref: "specs/occt_probe/occt_probe.exe --fbbox 1327",
    },
        { id: "D14", at: "2026-09-29", chose: "撤回 T-98 的「F1758 pcurve d0=NaN」结论：该 pcurve 功能正确，只存在自报参数域不一致",
       rejected: ["按「锥面 NaN pcurve」去改 make_pcurve_full / 补无限域投影分支"],
      
      why: "复核实测：pcRange 虽报 (-inf,inf)，但 d0 在 t∈(0,2π) 九点全部有限…",
      
      ref: "crates/occt-topo/src/pcurve_full/surface_projector.rs:2876-2891",
    },
        { id: "D15", at: "2026-09-29",
      chose: "T-69 的重定结论：端口拒的 7 个面在 GT 侧都出网格…",
       rejected: ["继续用面序号对拍（已由 D13 否证）", "把 7 个未网格面当成忠实的 RejectFace 结案"],
      
      why: "按 bbox 中心最近邻配对（临时 GT --pair 模式 + 端口 zz_pair_probe…",
      
      ref: "crates/occt-topo/src/pcurve_full/make_pcurve.rs:128-139 · BRepMesh_DefaultRangeSplitter.cxx:35-88",
    },
        { id: "D16", at: "2026-09-29",
      chose: "新增测量口径警告：端口的「无 stat 面数」在不同调用序列下不稳定（同一探针内 1757/15 与…",
       rejected: ["把某一次探针输出的 no_stat 数当成模型性质"],
      
      why: "同一进程内实测：先 prs3d_get_deflection 再 from_deflection 得 stats=175…",
      
      ref: "crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs:239-274",
    },
        { id: "D17", at: "2026-09-29", chose: "T-94 第一版实现（推迟 seam-shift 的 pcurve 写入到 result 确定后）实测惰性，已回退不提交",
       rejected: ["把「推迟 shift 写入」当成已消除副作用而提交", "以「行为中性」为由保留无实测收益的改动"],
      
      why: "实测与未修改基线逐字相同：T0M wires=2121 / faces_with_2plus=309 / WIREHIS…",
      
      ref: "crates/occt-topo/src/shhealing/shape_fix_face.rs:480-487 :604-608",
    },
        { id: "D18", at: "2026-09-29",
      chose: "撤回「F171 中点顶点是 reader 自造/未移植行为」的推论…",
      
      rejected: ["按「找 reader 里按中点断边建点的路径并删除/标 UNPORTED」去修 combine_vertex / fix_dummy_seam", "沿用「中点顶点是忠实行为、不存在缺陷」这个半对的结论"],
      
      why: "第一轮 BACKTRACE 定位到唯一调用链 resolve_face(read_topology.rs:714) → …",
      
      ref: "ShapeBuild_Vertex.cxx:38-72 · ShapeFix_Wire.cxx:4019 · :4213-4289 · ShapeFix_Face.cxx:365-520 · crates/occt-topo/src/shhealing/wire_fix.rs:3084",
    },
        { id: "D19", at: "2026-09-29",
      chose: "回退 fix_dummy_seam 的传参修复（`(wire, n1)` → `(wire, n2)…",
       rejected: [ "保留修复并把 a3n00 的 area_tol 0.15→0.18 以让它转绿", "以「未网格面 7→3」为由直接改基线而不量其余模型", ],
      
      why: "23 个门禁模型的面积比逐一对拍：只有 a3n00 变差（0.8627→0.8298，−0.0329）…",
      
      ref: "crates/occt-topo/tests/common/mod.rs:93 · crates/occt-topo/src/shhealing/wire_fix.rs:3084",
    },
        { id: "D20", at: "2026-09-29",
      chose: "看板自审定案：T-93 的 `dependsOn: [\「T-59\」]` 是悬空依赖（T-59 从…",
      
      rejected: [ "删掉 T-93 的 dependsOn 而不立 T-59（会让 a3n00 面积缺口再次消失成散文）", "把悬空依赖并进已有的「依赖未完成即开工」那条风险里（量级不同：未完成会自己解除，不存在不会）", ],
      
      why: "自审发现 3 条硬问题并已修：① 悬空依赖使 P0 静默挂死（`dependsOn` 全量校验：本次补立 T-59 **…",
      
      ref: "specs/board.canvas.tsx · crates/occt-topo/tests/common/mod.rs:91",
    },
        { id: "D21", at: "2026-09-30",
      chose: "撤回 D15 的配对表：端口拒的 7 个面在 GT 的**同名孪生面**上出网格…",
       rejected: [ "沿用 bbox 中心最近邻（不校验尺寸，半径差 0.5 的面可互配）", "按 --mesh/--face 序号配对（D13 已否证）", ],
      
      why: "D15 称 P405→G407 tri=185、P1691→G1697 tri=73、P1758→G1764 tri=7…",
      
      ref: "specs/occt_probe/occt_probe.exe --uvsum（本轮新增）· D15 · D13",
    },
        { id: "D22", at: "2026-09-30",
      chose: "T-69 重定：端口 7 面被拒**不是** RangeSplitter 判据缺陷 —— OCCT …",
      
      rejected: [ "继续在 RangeSplitter / IsValid / updateRange 一侧找缺口（三者已逐行等价）", "按「pcurve 参数域/取样」解释零宽（常宽方向是真实的几何退化边界，不是参数域 bug）", ],
      
      why: "① 判据同构不是推论而是实跑：`collectWirePoints`（BRepMesh_NodeInsertionMes…",
      
      ref: "specs/_a3n00_gap_analysis.md §9.288 · BRepMesh_NodeInsertionMeshAlgo.hxx:142-184 · BRepMesh_DefaultRangeSplitter.cxx:35-41 :202-236",
    },
        { id: "D23", at: "2026-09-30",
      chose: "T-93 的方向重定向：从「换 `fix_dummy_seam` 传参换 T0M 未网格 7→3」改…",
      
      rejected: [ "继续按 §9.287.5 的顺序（先 T-59 面积缺口、再应用 §9.287.4 一行改动）作为主线", "把缺边归到 `make_pcurve` / `pcurve_full`（这些 pcurve 是端口自己造的，但造的是**真实存在的那条**闭合圆）", ],
      
      why: "未网格面的成因既然在边界环，§9.287.4 那行改动（7→3）就不再是这条缺口的主线：它改的是 notch 边对…",
      
      ref: "specs/_a3n00_gap_analysis.md §9.288 · crates/occt-topo/src/meshing/model_builder/wire_builder.rs · crates/occt-topo/src/step/read_topology.rs",
    },
        { id: "D24", at: "2026-09-30",
      chose: "把本轮探针固化成可复跑仪器（OCCT `--uv` / `--uvsum`、端口 `zz_uv_fe…",
       rejected: [ "把库内 env-gated dump 留在源码里当长期仪器", "不固化仪器，下一轮再临时写一遍", ],
      
      why: "库内插桩污染被测代码（且读数与被测二进制同源，难判「插桩本身是否改变行为」）…",
      
      ref: "crates/occt-topo/examples/zz_uv_feed.rs · specs/occt_probe/occt_probe.cpp · §9.288.5",
    },
        { id: "D25", at: "2026-09-30",
      chose: "T-93 主线定位到**导入期 seam 修复**：STEP 里这些面声明的环边数本来就少于 OCC…",
      
      rejected: [ "整形状直接调用现有 fix_missing_seam（实测 212 面变动，含 126 个 2→4 的伪增量）", "按 §9.287.4 的 fix_dummy_seam 传参修复当主线（它换 notch 边对，不补缺边）", ],
      
      why: "取证链（全部可复跑）：① STEP 文本实测 `#29920=ADVANCED_FACE('',(#4223),#77,…",
      
      ref: "data/occ/T0M.stp #29920/#4223/#6320/#25460 · ShapeFix_Face.cxx:492-498 :1722-2330 · ShapeProcess_OperLibrary.cxx:801-899 · crates/occt-topo/src/shhealing/shape_fix_face.rs:147-148 · .target-gate/fixall.txt",
    },
        { id: "D26", at: "2026-09-30",
      chose: "定级结论：**reader 没有丢边**。端口 F1691 的 1 条边就是 STEP `EDGE_…",
      
      rejected: [ "继续在 resolve_face / resolve_loop / bind_edge_loop_vertices 里找「哪条边被丢了」——没有可丢的", "把 1 条边当成 reader 的回归（它忠实于文件）", ],
      
      why: "逐级核过，每一级都是 `file → 逐条进 → 无过滤`：① STEP `#29920=ADVANCED_FACE('…",
      
      ref: "crates/occt-topo/src/step/read_topology.rs:307 :369 :413 :564 :578 :628-644 · data/occ/T0M.stp #6320 · D25 · §9.289.7",
    },
        { id: "D27", at: "2026-09-30",
      chose: "T-93 落地：把导入期 seam 修复**接进 reader**（`resolve_face` 末…",
      
      rejected: [ "继续把 `fix_missing_seam` 留在 reader 之外（D2/T-93(a) 的处理）—— 实测 OCCT 成形期确实加了 seam，不接就等于永久少 3 条边", "为「消灭」211 个伪增量去改 `check_wire` / 配对选择（没有 .cxx 依据说明端口逻辑与 OCCT 不同）", ],
      
      why: "① 先证同构：`zz_seam_fix` 对 T0M f=1691 单独调用…",
      
      ref: "specs/_a3n00_gap_analysis.md §9.290 · ShapeFix_Face.cxx:492-498 :1899-2264 · ShapeProcess_OperLibrary.cxx:785-899 · crates/occt-topo/src/step/read_topology.rs",
    },
    {
      id: "D28",
      at: "2026-09-30",
      chose: "T-99 主线定在**读入期的 bound 合并**（把 2 个 FACE_BOUND 并成 1 条…",
      rejected: [ "继续在 Delaunay 侧（frontier_adjust / cleanup_mesh / 点序）找缺口",
        "继续按 §9.365.2 的按序号 GT 表当目标值" ],
      why: "① 等价性已证：把端口的完整结构（点 + 72 条约束链 + 容差/格/顶点序）喂给 OCCT 的 BRepMesh_D…",
      ref: "specs/_a3n00_gap_analysis.md §9.368 · ShapeFix_Face.cxx:492-498 :1722-2330 · crates/occt-topo/src/shhealing/shape_fix_face.rs:147 :406",
    },
    {
      id: "D29",
      at: "2026-09-30",
      chose: "把导入期 `ShapeFix_Face::FixMissingSeam` 步骤**接回 reader…",
      rejected: [ "继续把 seam 合并留在 reader 之外（靠 rescue 补三角形）",
        "用 repr_key 摘注册一类旁路去换同一个面积比（§9.315/§9.320）" ],
      why: "① 移除理由是「OCCT 在这条路径上不跑 `ShapeFix_Face::Perform`」（link map 实验…",
      ref: "specs/_a3n00_gap_analysis.md §9.369 · crates/occt-topo/src/step/read_topology.rs · ShapeFix_Face.cxx:482-498 :1722-2330 · D27 · D28",
    },
    {
      id: "D30",
      at: "2026-10-03",
      chose: "T0M 的薄面簇对拍用**同一次运行内的 OCCT 真值**做（新增测点 = 既有 `--faces…",
      rejected: [ "把 ④ 当成 T0M 0.9931→0.9809 的回归撤掉",
        "继续按 3D bbox 逐面配对来判 ④ 的优劣（OCCT 与端口的 bbox 在第 6 位小数有差，1778 个里 591 个对不上）" ],
      why: "① 测点同源：偏转取门禁同口径 `prs3d_get_deflection` = maxComp(bbox)*0.001…",
      ref: "specs/_a3n00_gap_analysis.md §9.581 · specs/occt_probe/occt_probe.cpp (--facestats/--wires) · crates/occt-topo/examples/zz_probe_a3n00.rs (--fstats/--ecensus/--wirehist) · ShapeBuild_ReShape.cxx:200-324 · D29",
    },
    {
      id: "D31",
      at: "2026-10-03",
      chose: "T0M 的 6 张残差面按 `ShapeFix_Face::FixMissingSeam` 返回多面…",
      rejected: [ "按 §9.580-F 去改 `CollectWires` 的 wire 合并（该假设已证伪：OCCT face 116 是 healing 新建面，不是端口同 bbox 的那张原始面）",
        "为了凑 1778 去加任何 OCCT 里没有的谓词/启发式（例如『若面是周期面就拆两张』）" ],
      why: "① 先量清方向：T0M.stp 自身 1772 个 `ADVANCED_FACE`…",
      ref: "specs/_a3n00_gap_analysis.md §9.582 · crates/occt-topo/src/step/read_topology.rs · ShapeFix_Face.cxx:2266-2268 · ShapeFix_Shape.cxx:200 :257 :294 · ShapeBuild_ReShape.cxx:282-299 · STEPControl_Controller.cxx:201 · D30",
    },
    {
      id: "D32",
      at: "2026-10-03",
      chose: "T-102 剩余缺口按 `FixShape` 子模式 bisect 钉成两条独立入口：**(1) 2…",
      rejected: [ "继续用 `--set FromSTEP.FixShape.*`（或读前的 `SetShapeFixParameters` / `SetShapeProcessFlags`）做 bisect —— actor 是 `ReadFile` 期间才建起来的（`XSControl_Reader.cxx:480-534`），读前调用是 no-op；8.0 又走 `DE_ShapeFixParameters` 显式 map（`XSAlgo_ShapeProcessor.cxx:609-620`），静态量根本不被读。本轮 T0M 上 14 个单点 `--set` 结果全部等于默认值即为证",
        "把 loop-wire 面按『周期面就拆两条 wire』这类 OCCT 里没有的…",
        "把 `FixMissingSeamMode=0` 的 1973 面当成『seam 步造了 201 张面』（实际是 seam 缺席后 `FixAddNaturalBound` 不再兜住、`FixSplitFace` 在 `cxx:711-717` 开火）" ],
      why: "① `FixFaceMode=0` ⇒ T0M 1772 面、`1:1463 2:294 3:9 5:1 6:1 7:2…",
      ref: "specs/_a3n00_gap_analysis.md §9.583 · specs/occt_probe/occt_probe.cpp (--fix/--set) · ShapeAnalysis_Wire.cxx:2228-2340 · ShapeFix_Face.cxx:583-598 :2398 :2478 · ShapeFix_Shape.cxx:196-210 :705-717 · XSControl_Reader.cxx:480-534 · XSAlgo_ShapeProcessor.cxx:609-620 · D31",
    },
    { id: "D33", at: "2026-10-04",
      chose: "T-102 的剩余缺口**不是** `ComposeShell` 少切/少分类（D32 的两条入口因…",
      rejected: [ "按 D32/§9.586-F(1) 去 `ShapeFix_ComposeShell.cxx:206-270` 的 `Perform` 分段里找『端口少切/少分类的那一段』（实测该处不缺：修朝向后面数与 OCCT 逐项相同）",
        "按 §9.585-C(2) 去 `ShapeFix_Face.cxx:2236-…",
        "为凑 1778 面去移植 `FixLoopWire` 或加任何 OCCT 里没有的谓词（`FixLoopWiresMode=0` 那条 bisect 结论仍成立，但它不是本轮残差的入口）" ],
      why: "① 对着 `.cxx` 定因：`ShapeFix_ComposeShell.cxx:629` 的 `TopoDS_Ite…",
      ref: "specs/_a3n00_gap_analysis.md §9.587 · crates/occt-topo/src/shape_fix_compose_shell/load_wires.rs · ShapeFix_ComposeShell.cxx:629-640 :586-607 · TopoDS_Iterator.cxx:26-83 · D32",
    },
    { id: "D34", at: "2026-10-04",
      chose: "`DispatchWires` 里「给新建边补 3D 曲线」改走真实 `BRepLib::Build…",
      rejected: [ "保留 §9.30 的旧体（`CurveOnSurface` 当 3D 曲线）—— 它给不出一条 OCCT 意义上的 3D 曲线，`SameRange`/容差/后续离散都不同源",
        "为凑某个模型的读数去调 `BuildCurve3d` 的连续性/阶数/`MaxSegment`（应走 `BRepLib.hxx:90-94` 的默认 `GeomAbs_C1` / `MaxDegree=14` / `MaxSegment=0` + `evaluateMaxSegment`）" ],
      why: "① 对着 `.cxx`：`ShapeFix_Edge::FixAddCurve3d`（`cxx:618-638`）在 `…",
      ref: "specs/_a3n00_gap_analysis.md §9.588 · BRepLib.cxx:149-183 :187-263 :275-295 :301-455 · ShapeBuild_Edge.cxx:714-775 · ShapeFix_Edge.cxx:335-464 :618-638 · GeomLib.cxx:559-679 :1051-1165 :2991-3077 :3079-3221 · ElCLib.cxx:1339-1422 · ShapeFix_ComposeShell.cxx:3506-3529 · D33",
    },
    { id: "D35", at: "2026-10-04",
      chose: "`ShapeFix_Face::Perform` 的周期锥面退化臂按 `.cxx` 顺序补齐：读入链…",
      rejected: [ "按面级 bbox/面积/密度门把「塌陷薄面」筛掉或跳过（OCCT 里不存在的规则）",
        "为 T0M 的读数去调 `FixPeriodicDegeneratedMode`/`RecadreOnPeriodic` 的默认值或偏转（默认值本身就是规格）" ],
      why: "① §9.590 钉死第一现场：face 1764 的 wire 塌成 `edges=[1] nv=1`、V 跨度 0 …",
      ref: "specs/_a3n00_gap_analysis.md §9.590 · §9.591 · ShapeFix_Face.cxx:144 :482-498 :1737-1741 :3018-3098 :3101-3259 · ShapeFix_Root.lxx:101 · ShapeProcess_OperLibrary.cxx:830 · D32",
    },
    { id: "D36", at: "2026-10-05",
      chose: "删除 `meshing/model_builder/wire_builder.rs::add_wir…",
      rejected: [ "把 `chain_area()` 换成真实 2D `wire_area_2d()`、或加 `|area|` 阈值/退化保护 —— 该分支本身在 OCCT 里就不存在，改成「更好的启发式」仍是非 OCCT 规则",
        "只对 `chain_area == -0.0` 这一处特判（等于给 motoc 调参，换模型即失效）" ],
      why: "① 现象：`data/occ/motoc.step` 大圆柱 / 右侧小圆柱上多出「横向 2 层环带」…",
      ref: "specs/_a3n00_gap_analysis.md §9.598 · §9.601 · BRepMesh_ShapeVisitor.cxx:99-134 :103 · ShapeExtend_WireData.cxx:583-588 · meshing/wire_order.rs · D34",
    },
    { id: "D37", at: "2026-10-06",
      chose: "`geom2d` 的 `Geom2dConvert` 全链（`CurveToBSplineCurve…",
      rejected: [ "继续用 `ReparamCurve2d` 包装器冒充 `GeomLib::SameRange` 的产出（§9.605：包装器不是 OCCT 会持有的类型，后接分派全落错臂）",
        "按某个模型的面积 / 三角数去调这些基元的常量或分支（基元是共享的，且这些臂在当前 23 模型上不改变行为）" ],
      why: "① §9.604-§9.623 是一条**潜伏缺口**链：这些分支在 23 个基准模型上不改变任何读数（四模型 `--w…",
      ref: "specs/_a3n00_gap_analysis.md §9.604-§9.623 · Geom2dConvert.cxx:244-262 :283-301 :323-345 :379-397 · GeomLib.cxx:871-888 :908-921 :960-968 · ShapeBuild_Edge.cxx:643-655 :687-697 · gp_Ax22d.hxx · gp_Circ2d.hxx · gp_Trsf2d.cxx:198-213 · gp_Trsf.cxx:90-99 :103-109 :207-213 · gp_Mat.cxx:242-334 · gp_Quaternion.cxx · gp_EulerSequence.hxx · gp_GTrsf.cxx · gp_Dir.cxx:129-155 · D34",
    },
    { id: "D38", at: "2026-10-08",
      chose: "`MakeFacesOnPatch` 的切向游走按 `cxx` 原文复刻（`while (stPoint == ON || UNKNOWN)` + 跨 `j` 共享的 `ew` 迭代器 + `aCL` 随 `aCW` 走），`CopyNMVertex` 面重载落地、边重载登记 UNPORTED",
      rejected: [ "保留端口原来的 `loop {}`（不复查 `stPoint`、`aCL` 恒为第 k 条边的 last、`eidx` 每个 `j` 从 k 重来）—— 那不是 `.cxx` 控制流，`stPoint` 已是 IN/OUT 时还会继续走边、把 `stPoint` 拖到错值",
        "把 `ew` 迭代位置在 `j` 之间重置（看似更『干净』，但 `.cxx:3034` 的 `TopoDS_Iterator ew(wr)` 是 i 循环内、j 循环外的单个对象，重置即偏离）",
        "为 `CopyNMVertex` 的边重载硬凑一个调用点（`ShapeFix_Wire_1.cxx` 只有 `FixGaps3d`/`FixGaps2d`/`FixGap3d`/`FixGap2d`，只在 `ShapeFix_Wireframe::FixWireGaps` 与 `ShapeUpgrade_WireDivide` 下可达，二者都未移植、也不在 STEP `FromSTEP` 序列里，且其 `BRep_PointRepresentation` 表示链在端口无数据模型）" ],
      why: "① §9.626 触发面实测（`MFP_PROBE`）：T0M 226 次调用 / 13 次多环、a3n00 47/1、acs10 114/4；多环路径（根识别、切向游走、孔分类）真实触发、`reverse` 两个取值都跑到 ⇒ 切向游走不是死代码。② 顶点孔与单环 `invert` 分支在 8 个 occ 模型上 0 触发 ⇒ 面重载属潜伏分支，但 `cxx:3241` 是它唯一的调用点，必须随 `MakeFacesOnPatch` 一起落地。③ 逐面 `--fstats` 与 §9.625 基线逐位相同（a3n00 11138/12440、acs10 37291/46548、T0M 60640/67524）⇒ 复刻无回归。",
      ref: "specs/_a3n00_gap_analysis.md §9.626 · ShapeFix_ComposeShell.cxx:2978-3271 :3034 :3104-3130 :3197-3209 :3238-3245 :3582 · ShapeAnalysis_TransferParametersProj.cxx:562-710 :715-802 · ShapeFix_Wire_1.cxx:76-98 :101-123 :154-862 :893-2005 · ShapeFix_Wireframe.cxx:92 :1051 :1617 · STEPControl_Controller.cxx:201 · ShapeFix_Root.cxx:26-30 · ShapeFix_Face.cxx:120-128 :1639-1641 · D37",
    },
    { id: "D39", at: "2026-10-08",
      chose: "把 `Geom2dInt_GInter::InternalPerform` 的 (Circle,Circle) 行切到专用 `IntCurveIntConicConic` 重载（`IntCurve_IntConicConic_1.cxx:58-151 :154-355 :807-1206`），并顺带补齐 (Circle,Line)/(Ellipse,Line) 的「反向」行（`gxx:346-358`、`:441-455` 的 `SetReversedParameters(true)`）",
      rejected: [ "继续让 (Circle,Circle) 落在通用 `IntConicCurveGen` 臂上 —— 数学等价，但专用重载在重合（`nbsol==3`）与 `TolTang` 切线退化上的分支更完整，而 `ShapeAnalysis_Wire::CheckSelfIntersectingEdge`/`FixSelfIntersectingEdge` 正是要在这些退化上收敛",
        "把 (Circle,Ellipse)/(Ellipse,Ellipse)/(Circle|Ellipse × Parabola|Hyperbola) 同批一起移植 —— 本轮只为解除自交修复链的前置，其余行仍由通用臂正确处理，先做最小必要集",
        "给 `CircleCircleGeometricIntersection` 的区间长度 > PI 分支加「按角度跨度/包围盒」的启发式替代 —— `_1.cxx:250-300` 的归一化 + 取补是原文控制流，不是可替换的启发式" ],
      why: "① §9.627 D.1 登记的第一条缺口，也是 D.3（`ShapeAnalysis_Wire::CheckSelfIntersectingEdge` / `ShapeFix_Wire::FixSelfIntersectingEdge`）的公共前置。② 新文件 `int_conic_conic_circle_circle.rs` 逐段对 `_1.cxx:58-151`（`ProjectOnC2AndIntersectWithC2Domain`）、`:154-355`（`CircleCircleGeometricIntersection`）、`:807-1206`（`Perform(Circ2d,Circ2d)`）；`gp_Circ2d::reversed()` 对 `gp_Circ2d::Reversed()`（`location`/`xdir`/`radius` 不变，只翻局部 `GpAx22d` 的 `vydir`）。③ **可达性实测**：全仓唯一实代码调用点 `split_by_line.rs:204` 的第一参数 `j_c2d` 是 `Geom2dLine`（由 `line` 直接构造，恒为线）⇒ 本批三条臂当前 0 触发、属潜伏分支；可达条件即 D.3 的 `CheckSelfIntersectingEdge`/`CheckIntersectingEdges` 移植。④ 门禁：`--all-targets` 0 error、`step_obj_gates` 4 passed、`export_data_obj` 15 ok / 0 err，读数与 §9.626 相同（潜伏分支，无行为变化）。",
      ref: "specs/_a3n00_gap_analysis.md §9.628 · IntCurve_IntConicConic_1.cxx:58-151 :154-355 :807-1206 :850 :869 · Geom2dInt_Geom2dGInter.gxx:346-358 :360-374 :441-455 · gp_Circ2d.hxx Reversed · D38",
    },
  ],
  nextAction:   { taskId: "T-109",
    action: "【T-109（2026-10-08 §9.628 闭环后）】`IntCurve_IntConicConic` 的 (Circle,Circle) 专用重载已落地并接线，(Circle,Line)/(Ellipse,Line) 反向行一并补齐；`ginter.rs` 里 T-100 登记的「真实缺口」至此清零。下一步候选（约束：缺口必须先在板上立卡）：(1) **建议立 T-110**：`ShapeAnalysis_Wire::CheckSelfIntersectingEdge`（`ShapeAnalysis_Wire.cxx:1269-1344`）+ `CheckIntersectingEdges`（`:1360-`）→ `ShapeFix_Wire::FixSelfIntersectingEdge` / `FixIntersectingEdges` —— 本批是它们的公共前置，前置已清，且移植后本批三条臂才由潜伏转为可达（§9.628-C 可达性实测：当前唯一调用点 `split_by_line.rs:204` 的第一参数恒为 `Geom2dLine`，本批 0 触发）；(2) `IntCurve_IntConicConic` 余下重载（Circle/Ellipse、Ellipse/Ellipse、Circle|Ellipse × Parabola|Hyperbola）仍走通用臂，属特化非新能力；(3) `bottom.step` 面积比 0.97412、`top.step` 导出三角对参考差 6；(4) `shhealing/wire_fix.rs` 3379 行 > 1000，按仓库规则需拆子模块。",
    why: "【T-109】§9.628：(Circle,Circle) 行原先落在通用 `IntConicCurveGen` 臂上，(Circle,Line)/(Ellipse,Line) 反向行落在 `default:` 臂上；现按 `Geom2dInt_Geom2dGInter.gxx:346-374 :441-455` 切到专用重载，`IntCurve_IntConicConic_1.cxx:58-151 :154-355 :807-1206` 三段按原文复刻（含 `TolTang` 切线管、`nbsol==3` 重合、区间 > PI 归一化取补、`nextafter(PIpPI,0)`）；门禁 `--all-targets` 0 error、`step_obj_gates`、`export_data_obj`。",
  },
  tasks: [
        { id: "T-99",
      lane: "mesh",
      title: "同结构喂 OCCT：18/18 逐字段完全相同 ⇒ Delaunay 层忠实…",
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
      goal: "把 a3n00 的面积比从当前真实值 **0.8627** 推向 1（缺口 31179.6）…",
      next: "**【唯一入口：段集合（BREAK-IN）同构性，须在 OCCT 侧取证】** 端口侧能对照的环节已全部判定忠实（§9.577 清单）…",
      acceptance: "a3n00 面积比从 **0.8627** 上升且不劣化；T0M 未网格 7 不增加（`f1e56776` 已记 T0M 6→7 的代价）…",
      evidence: "**2026-09-30，十三轮**。**(1) 面已 100% 锁定**（§9.303.1）：STEP `#5375` = `CYLIND…",
      write: "crates/occt-topo/src/shhealing/shape_fix_face.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.303 §9.314 §9.315 §9.316 §9.318 · crates/occt-topo/examples/zz_pcurve_key_probe.rs · crates/occt-topo/examples/zz_share_probe.rs",
      note: "【§9.321】13 条假设已否决。本轮新增两条：第 12 条「共享键的危害经 swap_seam …",
      dependsOn: [],
    },
        { id: "T-59", lane: "mesh", title: "a3n00 面积缺口：配对表已出，缺口定到单面 F113（已拆为 T-99）",
      status: "completed", priority: "P1", owner: "agent", progress: 100, estimate: 5, actual: 2, startedAt: "2026-09-30",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
      goal: "按 bbox 六坐标把 a3n00 的 226 个端口面与 GT 面配对，逐面定位构成面积亏损的 face class",
      
      next: "本卡收口：配对表已出（226/226 全配上），缺口已定到**单个面 F113**，后续修复拆到 T-99…",
       acceptance: "产出可复跑的配对表并定位亏损 face class（不要求修好）；不新写测试、不改断言",
      
      evidence: "**2026-09-30 实测**：`.target-gate/pair_a3n00.py .target-gate/a3n00_port.…",
      write: "", ref: "specs/_a3n00_gap_analysis.md §9.291 · §9.287 · tests/common/mod.rs:93 · .target-gate/pair_a3n00.py",
       note: "本卡原为 T-93 的前置（D20 补立）。分析类卡，**库源码零改动**…",
      
      dependsOn: [],
    },
        { id: "T-93", lane: "port-gap", title: "reader 缺带 Context 的 ShapeFix_Face::Perform（导入期 seam 修复）",
      status: "completed", priority: "P0", owner: "agent", progress: 100, estimate: 4, actual: 6, startedAt: "2026-09-28",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
      goal: "端口 F1691 一类「单闭合边」面的边界环补回 GT 的边数（GT 孪生面为 4 条），T0M 未网格面下降…",
      
      next: "本卡结项（D27）。留下的有界旁支，按优先级：① **f=1758（圆锥）** —— `FixMissingSeam` 只为球/BSpline/退化环面构造 s…",
       acceptance: "端口 F1691 一类的环边数与 GT 孪生面一致；T0M 未网格面下降；--lib 与 step_obj_gates 不劣化…",
      
      evidence: "**本轮（2026-09-30）落地实测，只切一个 hunk**（`read_topology.rs` 的 `resolve_face` 末…",
       write: "crates/occt-topo/src/step/read_topology.rs",
      
      ref: "specs/_a3n00_gap_analysis.md §9.287 · §9.288 · §9.289 · §9.290 · ShapeFix_Face.cxx:492-498 :1652-1718 :1899-2264 · ShapeProcess_OperLibrary.cxx:785-899 · STEPControl_Controller.cxx:201 :221",
      
      note: "工具（均已固化）：`zz_seam_fix`（单面修复前后 + 整形状过触发量化）、`zz_uv_f…",
      
      dependsOn: [],
    },
        { id: "T-69", lane: "mesh", title: "T0M 未网格面（6→7）：根因已定位为「边界环缺边」（读入层），非 mesh/判据缺陷",
      status: "completed", priority: "P1", owner: "agent", progress: 100, estimate: 3, actual: 5, startedAt: "2026-09-24",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
      goal: "查清 T0M 未网格面（7 个）的成因：判据缺陷、还是喂进判据的点集缺陷",
      
      next: "本卡收口，后续动作转 T-93（D23）：缺边的修补属 reader 的边界环构建，不在本卡的 mesh 范围内…",
       acceptance: "给出「未网格面成因」的可复核结论，并留下能复跑的对拍仪器（不要求修好）",
      
      evidence: "**结论（D22）**：不是 RangeSplitter/IsValid/updateRange 的判据缺陷 —— 三者与 OCCT 逐行同…",
       write: "",
      ref: "specs/_a3n00_gap_analysis.md §9.288 · BRepMesh_NodeInsertionMeshAlgo.hxx:142-184 · BRepMesh_DefaultRangeSplitter.cxx:202-236",
      
      note: "本卡为分析卡，**库源码零改动**（库内临时插桩已摘除，见 D24…",
      
      dependsOn: [],
    },
        { id: "T-54", lane: "mesh", title: "单面回退 wireframe_face_triangulation：6 个门禁模型上从未被触达（保留，注释已更正）",
      status: "completed", priority: "P2", owner: "agent", progress: 100, estimate: 4, actual: 2, startedAt: "2026-09-20",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
      goal: "回退路径不再是 OCCT 里不存在的耳切/桥洞，而是约束 Delaunay（或按 OCCT 直接暴露主线报错）",
      
      next: "结项（§9.306）：**先只做判断**这一步做完了，结论是**两条路都不走**…",
       acceptance: "确认回退有无实测触发例；据实测决定换 Delaunay / 删除 / 保留并标注…",
      
      evidence: "**2026-09-30 实测（插桩 `OCCT_TOPO_TRACE_FALLBACK`…",
       write: "crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs",
      
      ref: "specs/_a3n00_gap_analysis.md §9.306 · BRepMesh_BaseMeshAlgo.cxx:52-62 · crates/occt-topo/src/meshing/incremental_mesh/discret_root.rs:394 :449-454 :461-463 :479-480 :549",
      
      note: "由 T-71 结项时经 D-201 重开；D22 解除了 T-69 的阻塞…",
      
      dependsOn: [],
    },
        { id: "T-94",
      lane: "port-gap",
      title: "共享键写回不是可解耦的副作用：实验证明它是承重的（消除后 F113 出网格…",
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
      next: "结项（§9.314 + §9.315）—— **结论是「解耦」这个前提不成立**…",
      acceptance: "T0M 的 wire 结构在不回写共享 registry 的前提下不劣化；F113 仍出网格…",
      evidence: "**写入点（2026-09-30）**：`ComposeShell` 共 6 处写 registry —— `dispatch_wires.…",
      write: "crates/occt-topo/src/shape_fix_compose_shell/wire_data.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.296 §9.296.6 §9.314 §9.315 · crates/occt-topo/examples/zz_share_probe.rs · crates/occt-topo/src/tgeometry.rs:524-531",
      note: "本卡由 T-93/T-69 的旁支立卡。**净结论：验证性的负面结果** —— 机制被实验证实…",
      dependsOn: [],
    },
        { id: "T-95", lane: "port-gap", title: "Geom2dInt 曲线适配层与 TheProjPCurOfGInter：已实现且已接线（卡面描述过期）",
      status: "completed", priority: "P1", owner: "agent", progress: 100, estimate: 8, actual: 2, startedAt: "2026-09-28",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
      goal: "把 SplitByLine 依赖的 2D 求交内核按 .cxx 补齐…",
      
      next: "结项（§9.297）：卡面点名的全部已实现且已接线，无需动代码。**真实余量已拆为 T-100**（`ginter.rs` 的两条分派臂）…",
       acceptance: "对应 UNPORTED 清除；--lib 1281/0 与 step_obj_gates 5/5 不劣化",
      
      evidence: "**逐项核对（2026-09-30，`crates/occt-geom2d/src/geom2d_int/`）**：`curve_tool.…",
       write: "crates/occt-geom2d/src/geom2d_int/curve_tool.rs",
       ref: "specs/_a3n00_gap_analysis.md §9.297 · Geom2dInt_Geom2dCurveTool.lxx · Geom2dInt_TheProjPCurOfGInter_0.cxx:27-79",
      
      note: "【2026-10-08 补记 §9.628】T-100 名下登记的「真实缺口」（`ginter.rs` 的 (Circle,Line)/(Ellipse,Line) 反向行与 (Circle,Circle) 行）已全部收口（见 T-109）。结项理由：卡面点名的范围（curve_tool + proj_p_cur）已全部完成且已接线…",
      
      dependsOn: [],
    },
        { id: "T-100", lane: "port-gap", title: "ginter 空分支已移除（通用臂本可胜任）；typ1 != Line 臂确认为本仓不可达",
      status: "completed", priority: "P2", owner: "agent", progress: 100, estimate: 4, actual: 1, startedAt: "2026-09-30",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
      goal: "按 .cxx 补齐 Geom2dInt_GInter 分派里剩下的两条臂，去掉这两处 UNPORTED",
      
      next: "结项（§9.305）。**臂一（Line/Conic）**：原代码在 `Circle|Ellipse|Parabola|Hyperbola` 上是个**空分支*…",
       acceptance: "`ginter.rs` 的两处 UNPORTED 清掉或改为有依据的「无触发例」标注…",
      
      evidence: "**2026-09-30 实测**。**臂一**：新探针 `crates/occt-geom2d/examples/zz_lineconic…",
       write: "crates/occt-geom2d/src/geom2d_int/ginter.rs",
      
      ref: "specs/_a3n00_gap_analysis.md §9.305 · crates/occt-geom2d/examples/zz_lineconic_probe.rs · IntCurve_IntConicConic_1.cxx:2236 :2861 · IntCurve_IntCurveCurveGen.gxx:339- · occt-core/src/intcurve/iconic_tool.rs",
      
      note: "由 T-95 结项拆出。**核心发现：那个 `Circle|Ellipse|Parabola|Hyp…",
      
      dependsOn: [],
    },
        { id: "T-96", lane: "arch", title: "GeomTrimmedCurve 参数表示：评估结论 = 不重做（原前提已被实测否证）",
      status: "completed", priority: "P2", owner: "agent", progress: 100, estimate: 4, actual: 1, startedAt: "2026-09-30",
      updatedAt: "2026-09-30", completedAt: "2026-09-30", blocker: "",
      goal: "判定是否重做 TrimmedCurve 参数表示（保留 basis 区间）…",
      
      next: "结项：**不重做**（§9.295）。原前提「两条 arm 无法忠实表达」已被实测否证 —— ① `compute_type`（`gcpnts.rs:66`…",
       acceptance: "产出一页评估（改动面 / 消费方清单 / 回归风险）并在本板登记结论；不做代码改动",
      
      evidence: "**2026-09-30 实测**：`gcpnts.rs` 的 `compute_type` / `compute_abs_composit…",
       write: "crates/occt-geom/src/trimmed.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.295 · Geom_TrimmedCurve.cxx:212-243 :255-267 · GCPnts_AbscissaPoint.cxx:26-65 :87-158 · occt-geom/src/gcpnts.rs:66 :102 :491 :579 · occt-geom/src/trimmed.rs:45,130,147",
       note: "评估类任务，交付即结论：不重做参数表示。T-51 的豁免理由（「余下两条要先改参数表示」）据此不成立…",
      
      dependsOn: [],
    },
        { id: "T-98", lane: "mesh",
      title: "锥面 pcurve 报 (-inf,inf)：OCCT 语义对照完成 —— 与 Geom2…",
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
      next: "结项（§9.313）：走卡面第一条验收路径 —— **给出 OCCT 语义对照的结论并保持现状**，未改 mock/门禁…",
      acceptance: "（可选卡）要么给出 OCCT 语义对照的结论并保持现状，要么让自报参数域与可行域一致；两者都不得改动 mock/门禁",
      evidence: "**2026-09-30 OCCT 语义对照（§9.313）**。**原有实测（2026-09-29…",
      write: "crates/occt-topo/src/pcurve_full/make_pcurve.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.313 · Geom_ConicalSurface.cxx:207-215 · Geom2d_Line.cxx:142-150 · Geom2d_Line.hxx:97-100 · GeomProjLib.cxx:81 :118-128 · crates/occt-geom2d/src/line.rs:37-38",
      note: "本卡已从 P1「NaN 缺陷」降级为 P3「自报参数域不一致」…",
      dependsOn: [],
    },
        { id: "T-97", lane: "hygiene",
      title: "「never used」定性为「移植已做、消费方未移植」；A 类 60 个模块已标注…",
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
      next: "结项（§9.311）。**方向已按用户选择改为「标注而非删除」**，且已落地。**定性**：这份代码库在逐模块移植 OCCT，OCCT 的调用方往往还没移植…",
      acceptance: "产出一页评估（改动面 / 消费方清单 / 回归风险）并在本板登记结论；不做代码改动",
      evidence: "**2026-09-30，两批误删全部回退，改为标注交付**。**口径更正**：此前引用的「561 条」是两个来源混合的产物（`rtk` 是…",
      write: "crates/occt-topo/src/lib.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.307 §9.308 §9.309 §9.310 §9.311 · .target-gate/final_census.py · .target-gate/annotate_clean.py · .target-gate/raw_warn.txt",
      note: "**本卡的产出是「判据 + 降噪」，不是删除量**（净删除 0 行）…",
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
      
      evidence: "HEAD 08c4edf2 实测 box∪cyl FUSE = 1 solid / 8 面 / volume 8.502655…",
      write: "crates/occt-topo/src/bop_bop.rs", ref: "BOPAlgo_BOP.cxx:583-711 · BOPAlgo_Builder_3.cxx::BuildSplitSolids",
       note: "AC 未覆盖的其他曲面布尔组合未复核",
      dependsOn: [],
    },
        { id: "T-88", lane: "bop", title: "phase10 curved_face_fillet_sphere_plane", status: "completed", priority: "P0", owner: "agent",
      progress: 100, estimate: 1, actual: 1, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "phase10 的球面/平面倒圆用例转绿", next: "—", acceptance: "--test phase1…",
      evidence: "HEAD 08c4edf2 实测 phase10 8/8；旧板记的 7/8 是 2026-09-21 快照，原诊断未再复现",
      write: "crates/occt-topo/src/bop_bop.rs", ref: "BOPAlgo_BuilderSolid · BOPTools_AlgoTools", note: "—",
      dependsOn: [],
    },
        { id: "T-41", lane: "bop", title: "摘除 bop_curved 的体素/网格布尔（OCCT 无对应）",
      status: "completed", priority: "P1", owner: "agent", progress: 100, estimate: 2, actual: 3, startedAt: "2026-09-20",
      updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "", goal: "布尔入口改派到忠实路径，去掉 OCCT 里不存在的体素/采样布尔",
      next: "—", acceptance: "摘除后 --lib 与门禁不劣化",
      evidence: "改派后实测 --lib 1281/0、step_obj_gates 5/5、phase10/19/3/4 = 8/5/4/9…",
       write: "crates/occt-topo/src/bop_curved/region_trim.rs",
      ref: "—", note: "classify_face_general 本身仍是 8×8 采样，不是 BRepClass3d_SolidClassifier",
      
      dependsOn: [],
    },
        { id: "T-67", lane: "port-gap", title: "Extrema_ExtPExtS / ExtPRevS 两臂", status: "completed", priority: "P2", owner: "agent",
      progress: 100, estimate: 3, actual: 3, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "补上点-曲面极值搜索的两条 OCCT 分支，摘掉替换实现",
      next: "—", acceptance: "occt-geom --lib 143/0 + 门禁不劣化",
      evidence: "两臂已移植并分派（point_surface_extrema.rs:60-61/:88-91、:450 perform_ext_ps、:48…",
      write: "crates/occt-geom/src/extrema_surf", ref: "Extrema_ExtPS.cxx:292-343", note: "—",
      dependsOn: [],
    },
        { id: "T-25", lane: "arch", title: "GeometryRegistry 侧表 → 几何进 TShape", status: "completed", priority: "P2", owner: "agent",
      progress: 100, estimate: 6, actual: 6, startedAt: "2026-09-20", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "几何不再走全局侧表，而是挂在 TShape 上", next: "—", acceptance: "--lib + 门禁…",
      evidence: "三张几何表已删除（tgeometry.rs:163-181 注释写明 all live on their own TShape now）…",
      write: "specs/_design_architecture_t25_t28.md", ref: "specs/_design_architecture_t25_t28.md",
      note: "保留 global() 兼容外壳是既定停点，不是未完成项",
      
      dependsOn: [],
    },
        { id: "T-28", lane: "arch", title: "ImpPrm HVertex 合并 + intana/intpatch 重叠合并",
      status: "completed", priority: "P2", owner: "agent", progress: 100, estimate: 8, actual: 8, startedAt: "2026-09-20",
      updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "", goal: "ImpPrm 的 HVertex 合并与 intana 的 closed form 单一化",
      next: "—", acceptance: "--lib + 门禁不劣化",
      evidence: "步 1–3（ImpPrm HVertex 合并）与步 4（intana/intpatch closed form 合并，五副本清零）均完成…",
      write: "specs/_design_architecture_t25_t28.md", ref: "IntPatch_ImpPrmIntersection.cxx:221-469",
      note: "若日后重开步 5，前提是先把 PatchIntersection 对「一般 B 样条 × 球」修正到…",
      
      dependsOn: [],
    },
        { id: "T-29", lane: "hygiene", title: "过期规格刷新", status: "completed", priority: "P3", owner: "agent", progress: 100, estimate: 1,
      actual: 1, startedAt: "2026-09-21", updatedAt: "2026-09-28", completedAt: "2026-09-28", blocker: "",
      goal: "规格文件不再指向已删除的 _board.md，且过期内容有横幅",
      next: "—", acceptance: "—",
      evidence: "两个目标文件已带「仅历史（T-29 登记）」横幅（_coverage.md:1、_brepmesh_align_review.md:3）…",
      write: "specs/_coverage.md", ref: "—", note: "—",
      dependsOn: [],
    },
        { id: "T-11", lane: "hygiene", title: "occt-topo 编译警告清理", status: "cancelled", priority: "P3", owner: "agent", progress: 0,
      estimate: 1, actual: 4, startedAt: "2026-09-26", updatedAt: "2026-09-28", completedAt: "", blocker: "",
      goal: "清理 unused import 类警告", next: "已豁免：剩余项经实验判定不可安全清除…",
       acceptance: "无行为变化",
      evidence: "四批共清 11 处真实 unused import，每批后 --lib 1281/0…",
      write: "crates/occt-topo/src", ref: "—", note: "豁免结论见 D7",
      dependsOn: [],
    },
        { id: "T-51", lane: "port-gap",
      title: "gcpnts 余量核对完成：UNPORTED 清单已对齐代码…",
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
      goal: "核对 gcpnts 的真实余量：两条 Compute arm 是否真的已完整…",
      next: "结项（§9.312）。两条结论：① **两条 Compute arm 都已实现** —— `gcpnts.rs` 里 `:488`（旧行号）声称 `GCPnts…",
      acceptance: "gcpnts 的 UNPORTED 清单与代码现状一致（要么清空、要么逐条列出真实缺口）…",
      evidence: "**2026-09-30 核对（§9.312）**。**改动前 5 处 UNPORTED…",
      write: "crates/occt-geom/src/gcpnts.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.295 §9.312 · GCPnts_UniformAbscissa.cxx:494-515 :510 · GCPnts_AbscissaPoint.cxx:26-65 :87-90 :96-158 :187-295 :436-474 · CPnts_MyRootFunction.cxx:32-37",
      note: "由 T-96 的结论重开。**结论是「余量比卡面说的小」**：两条 arm 早已实现、外形差异不成立…",
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
      next: "不改，作为已知记录", acceptance: "—", evidence: "UV-grid 与 deflection-adaptive 的采样差异…",
      write: "crates/occt-topo/tests/common/mod.rs", ref: "—", note: "—",
      dependsOn: [],
    },
        { id: "T-101",
      lane: "mesh",
      title: "F113（螺帽顶部圆锥斜切面）恢复为 OCCT 的 1 面 / 2 wire（22+6）并…",
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
      goal: "让 a3n00 的 F113（螺帽顶部圆锥斜切面）恢复为 OCCT 的 1 面 / 2 wire 并出网格…",
      next: "**【收口（§9.580）：a3n00 ratio 0.9963→**1.0000**、`WIREHIST 1:208 2:9 4:1 6:4 10:4`（与 …",
      acceptance: "a3n00 面积比从 **0.8996** 上升且不劣化（step_obj_gates 的 step_obj_area 实测口径）…",
      write: "crates/occt-topo/src/shhealing/wire_fix.rs · crates/occt-topo/src/shape_fix_compose_shell/make_faces_on_patch.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.579 · ShapeBuild_Edge.cxx:103-114 · TopoDS_Builder.cxx:74-91 · TopExp.cxx:214-253 · TopoDS_Iterator.cxx:26-83 · ShapeFix_ComposeShell.cxx:3100/3110/3223",
      dependsOn: [],
      evidence: "§9.580（本环境实测）：a3n00 `WIREHIST 1:208 2:9 4:1 6:4 10:4`、`our=227126.78 o…",
      note: "① 的根因**已在 ComposeShell 内部定位**（§9.377）：第一处分歧 = OCCT…",
    },
        { id: "T-102",
      lane: "mesh",
      title: "T0M 残差收口：`load_wires.rs::wire_data_edges` 对 R…",
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
      goal: "把 T0M 的面/wire 结构对齐到 OCCT：`--wires` 真值 `TOTAL faces=1778 wire…",
      next: "**已收口（§9.587，D33）** ⇒ 网格侧残差另立 T-103（T0M 逐面密度 +6.0%、薄面簇 1 张不出网格、4 张多 wire 面的 `MUL…",
      acceptance: "**已满足（§9.587）**：T0M `--wirehist` = `TOTAL faces=1778 wires=1921` + `1:…",
      write: "crates/occt-topo/src/shape_fix_compose_shell/load_wires.rs（`wire_data_edges`，§9.587 收口点）· crates/occt-topo/src/shhealing/shape_fix_face.rs（`FixLoopWire` + 调用点 + wire 首段）· crates/occt-topo/src/shhealing/wire_fix.rs（`ShapeFix_Wire::Analyzer`/`FixReorder` 支撑）· crates/occt-topo/src/step/read_topology.rs（resolve_face / resolve_shell）",
      ref: "specs/_a3n00_gap_analysis.md §9.587 · §9.586 · §9.585 · §9.583 · ShapeAnalysis.cxx:268-273 :48-62 · ShapeFix_ComposeShell.cxx:206-270 :586-607 :629-640 · TopoDS_Iterator.cxx:26-83 · ShapeAnalysis_Wire.cxx:2228-2340 · ShapeFix_Face.cxx:583-598 :1996-2127 :2077-2078 :2107 :2396-2404 :2398 :2478 :365-480 :1722-2330 :2236-2325 :2266-2268 · ShapeFix_Shape.cxx:196-210 :705-717 :200 :257 :294 · ShapeBuild_ReShape.cxx:282-299 · XSControl_Reader.cxx:480-534 · XSAlgo_ShapeProcessor.cxx:609-620 · D33 · D32 · D31",
      dependsOn: ["T-101"],
      evidence: "**2026-10-04，§9.587 实测（收口：`load_wires.rs::wire_data_edges` 去掉二次反转…",
      note: "【§9.587 收口】D32 的两条入口（『`ComposeShell` 少切/少分类』与『`Sha…",
    },
        { id: "T-103",
      lane: "mesh",
      title: "T0M 网格侧残差：ON-only 塌陷面 8/8 归位（§9.608：pcurve 求值…",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 6,
      actual: 0,
      startedAt: "2026-10-04",
      updatedAt: "2026-10-06",
      completedAt: "2026-10-06",
      blocker: "",
      goal: "把 T0M 的**网格**读数对齐到 OCCT。拓扑侧已收口（T-102/§9.587：1778 面 / 1921 wi…",
      next: "**【已停止 — 2026-10-06】本线在 §9.608 闭环后停止：ON-only 塌陷第一现场曾是 §9.603 分类器多边形自重叠（`collect_…",
      acceptance: "T0M `--fstats` 的 `matched/unmatched` 与逐面 `mv/mt` 对 OCCT 收敛（`unmatched`…",
      write: "crates/occt-topo/src/meshing/**（UV 覆盖 / 参数域离散 / 三角化，按命中层）· crates/occt-topo/src/shhealing/shape_fix_face.rs（若落在三角化前的面修补）",
      ref: "specs/_a3n00_gap_analysis.md §9.603 · §9.608 · §9.602 · §9.600 · §9.599 · §9.597 · §9.596 · §9.595 · §9.587-F · §9.581 · §9.373 · ShapeAnalysis_Edge.cxx:192-208 · BRep_Tool.cxx:327-361 · ShapeAnalysis_Wire.cxx:593-651 · ShapeBuild_Edge.cxx:206-334 · ShapeFix_Face.cxx:400 :532 :379-380 :2266 :2268 · ShapeFix_ComposeShell.cxx:505-506 :589-606 :629-630 :2344 :2767 :3451-3457 :3348-3357 · ShapeExtend_WireData.cxx:114-121 · ShapeBuild_ReShape.cxx:284 :317 · TopoDS_Iterator.cxx:73-82 · ShapeFix_Shape.cxx:257 · BRep_TFace.cxx:37-44 · BRep_CurveOnSurface.cxx:58-62 · BRepTools_ReShape.hxx:239-240 · ShapeBuild_ReShape.cxx:249-254 · BRepMesh_DefaultRangeSplitter.cxx:35-41 :45-81 :202-236 · BRepMesh_NodeInsertionMeshAlgo.hxx:142-184 · BaseMeshAlgo.cxx:52-62 · D22 · D30",
      dependsOn: ["T-102"],
      evidence: "**【已停止 — 2026-10-06】T-103 在 §9.608 后停止（`pf=812`/`pf=28` 已由 §9.608 归位…",
    },
        { id: "T-104",
      lane: "port-gap",
      title: "`BRepLib::BuildCurve3d` / `SameRange` / `Shap…",
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
      goal: "把 `DispatchWires` 里「给新建边补 3D 曲线」这条链从`CurveOnSurface` 假曲线换成 O…",
      next: "**已收口（§9.588）**：`cargo check --all-targets` 0 errors、`step_obj_gates` 5/5…",
      acceptance: "**已满足**：编译 0 error、`step_obj_gates` 5/5、4 个基准模型 `--wirehist` 不回落且与 OCC…",
      write: "crates/occt-topo/src/shhealing/shape_build_edge.rs · crates/occt-topo/src/shhealing/shape_fix_edge.rs · crates/occt-topo/src/brep_lib_same_range.rs · crates/occt-geom/src/geom_lib.rs · crates/occt-core/src/elib/clib.rs · crates/occt-topo/src/shape_fix_compose_shell/dispatch_wires.rs · crates/occt-geom2d/src/{curve,offset,offset2d,bezier_curve,bspline_curve}.rs · crates/occt-topo/src/tgeometry.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.588 · BRepLib.cxx:149-183 :187-263 :275-295 :301-455 · ShapeBuild_Edge.cxx:714-775 · ShapeFix_Edge.cxx:335-464 :618-638 · GeomLib.cxx:559-679 :1051-1165 :2991-3077 :3079-3221 · ElCLib.cxx:1339-1422 · ShapeFix_ComposeShell.cxx:3506-3529 · D34",
      dependsOn: ["T-102"],
      evidence: "**2026-10-04，§9.588 实测（本机）**：① 源侧：新增 `brep_lib_same_range.rs`（`CheckSa…",
      note: "UNPORTED（文件头已标）：`GeomLib::To3d` 的有理 Bezier / 周期 BS…",
    },
        { id: "T-105",
      lane: "port-gap",
      title: "锥面周期退化环移植：`IsPeriodicConicalLoop` + `ShapeFix…",
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
      goal: "按 `ShapeFix_Face.cxx:3018-3259` 补齐 `Perform` 在 `FixMissingSe…",
      next: "**已收口（§9.591）**。旁支（未动，见 §9.591-E）：(1) face 1764 端口 `edges=[5]` vs OCCT `edges_pe…",
      acceptance: "**已满足**：`cargo check --manifest-path crates/occt-topo/Cargo.toml --off…",
      write: "crates/occt-topo/src/shhealing/shape_fix_face.rs · crates/occt-topo/src/step/read_topology.rs · crates/occt-topo/examples/zz_probe_a3n00.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.591 · §9.590 · ShapeFix_Face.cxx:482-498 :144 :1737-1741 :3018-3098 :3101-3259 :3226-3232 :3239-3257 · ShapeFix_Root.lxx:101 :34 · ShapeProcess_OperLibrary.cxx:830 · D35",
      dependsOn: ["T-104"],
      evidence: "**2026-10-04，§9.591 实测（本机）**：① 缺口：`ShapeFix_Face::Perform`（`cxx:482-49…",
      note: "OCCT 侧对照：`output/t0m_occt_wires.txt` `FACE 1764 ty…",
    },
        { id: "T-106",
      lane: "mesh",
      title: "motoc 圆柱面多出横向环带（§9.598，13944 三角与 OCCT 逐位相同）+ …",
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
      goal: "修用户报告的 motoc 缺陷：`data/occ/motoc.step` 端口渲染截图里大/小圆柱面上各多出一处横向环…",
      next: "**已收口（§9.598 + §9.601 + §9.604-§9.607 + §9.609-§9.623）**…",
      acceptance: "**已满足**：`cargo check --all-targets` 0 errors；`step_obj_gates` 5/5…",
      write: "crates/occt-topo/src/meshing/model_builder/wire_builder.rs · crates/occt-topo/src/meshing/wire_order.rs · crates/occt-geom2d/src/bspline_curve.rs · crates/occt-geom2d/src/trimmed.rs · crates/occt-geom2d/src/geom2d_convert.rs · crates/occt-geom2d/src/convert_approx_curve.rs · crates/occt-geom2d/src/comp_curve_to_bspline.rs · crates/occt-geom2d/src/bezier_curve.rs · crates/occt-core/src/adv_approx/ · crates/occt-topo/src/geom_bnd_lib_bspline2d.rs · crates/occt-topo/src/pcurve_full/surface_projector.rs · crates/occt-topo/src/shhealing/wire_fix.rs · crates/occt-topo/src/shhealing/shape_build_edge.rs · crates/occt-core/src/gp/ax22d.rs · crates/occt-core/src/gp/circ2d.rs · crates/occt-geom2d/src/curve.rs · crates/occt-geom2d/src/circle.rs · crates/occt-core/src/gp/trsf2d.rs · crates/occt-core/src/gp/dir2d.rs · crates/occt-core/src/gp/ax2d.rs · crates/occt-core/src/gp/lin2d.rs · crates/occt-core/src/gp/elips2d.rs · crates/occt-core/src/gp/hypr2d.rs · crates/occt-core/src/gp/parab2d.rs · crates/occt-core/src/gp/dir.rs · crates/occt-core/src/gp/ax1.rs · crates/occt-core/src/gp/ax2.rs · crates/occt-core/src/gp/ax3.rs · crates/occt-core/src/gp/lin.rs · crates/occt-core/src/gp/circ.rs · crates/occt-core/src/gp/elips.rs · crates/occt-core/src/gp/hypr.rs · crates/occt-core/src/gp/parab.rs · crates/occt-core/src/gp/sphere.rs · crates/occt-core/src/gp/cylinder.rs · crates/occt-core/src/gp/cone.rs · crates/occt-core/src/gp/pln.rs · crates/occt-core/src/gp/torus.rs · crates/occt-core/src/gp/trsf.rs · crates/occt-core/src/gp/mat.rs · crates/occt-core/src/gp/quaternion.rs · crates/occt-core/src/gp/euler_sequence.rs · crates/occt-core/src/gp/mod.rs · crates/occt-core/src/gp/mat2d.rs · crates/occt-core/src/gp/trsf2d.rs · crates/occt-core/src/gp/gtrsf.rs · crates/occt-core/src/gp/gtrsf2d.rs · crates/occt-core/src/gp/vec.rs · crates/occt-core/src/gp/vec2d.rs · crates/occt-core/src/gp/dir2d.rs · crates/occt-core/src/gp/dir.rs · crates/occt-topo/src/iges.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.598 · §9.601 · §9.604 · §9.605 · BRepMesh_ShapeVisitor.cxx:99-134 · ShapeExtend_WireData.cxx:583-588 · BRepBndLib.cxx:109 · GeomBndLib_BSplineCurve2d.cxx:34-57 · GeomBndLib_Curve2d.cxx:154-157 · Geom2d_BSplineCurve.cxx:707 :343 :989 · GeomLib.cxx:842-969 :947-958 :961-969 :920-921 · Geom2dConvert.cxx:347-351 :199-209 :420-423 :379-397 :425-440 · ShapeBuild_Edge.cxx:679-686 :688-698 · GeomBndLib_BSplineSurface.cxx:30 :103 :118 · §9.606 · Geom2dConvert.cxx:181-449 :211-226 :228-264 :266-303 :305-312 :314-321 :347-351 · Geom2dConvert.hxx:169-171 · Convert_ParameterisationType.hxx · Convert_CircleToBSplineCurve.cxx:47-172 · §9.607 · Geom2dConvert_CompCurveToBSplineCurve.cxx:29-243 · Geom2dConvert_CompCurveToBSplineCurve.hxx:28-68 · Geom2d_BezierCurve.cxx:356-391 · Geom2d_BSplineCurve.cxx:236-292 :410-478 · Geom2dConvert.cxx:244-262 :283-301 :323-345 · ProjLib_ProjectedCurve.cxx:179 · BRepOffset_Inter2d.cxx:1154 · §9.609 · Geom2dConvert_ApproxCurve.cxx:31-203 · Geom2dConvert_ApproxCurve.hxx:29-92 · Geom2dAdaptor_Curve.cxx:285-288 :409-486 :490-573 :453-480 :536-566 · AdvApprox_ApproxAFunction.cxx:364-599 :603-630 :663-956 · §9.610 · ShapeBuild_Edge.cxx:596-700 :660-678 :679-687 :688-697 · Geom2dConvert_ApproxCurve.cxx:111-118 · Geom2d_TrimmedCurve.cxx:98-154 · §9.611 · GeomLib.cxx:863-921 :871-888 :890-900 :908-921 · gp_Circ2d.hxx:181-184 :335-344 :347-355 · gp_Ax22d.hxx:360-367 · gp_Dir2d.hxx:210-213 · gp_XY.hxx:159-163 · gp_Trsf2d.hxx SetRotation · §9.612 · GeomLib.cxx:908-921 :924-969 :929-958 :960-968 · Geom2dConvert.cxx:187-209 :442-444 · Geom2d_BSplineCurve.hxx SetKnots · BSplCLib.cxx:756-798 · occt-geom/src/bspline_curve.rs:91 · §9.613 · ShapeBuild_Edge.cxx:643-655 :685-687 :687-697 :596-700 · §9.614 · gp_Circ2d.hxx:181-184 :199-207 :309-315 :335-344 · gp_Circ2d.cxx:23-38 · gp_Ax22d.hxx:335-341 :346-354 :360-367 · gp_Ax22d.cxx:29-52 · gp_Dir2d.hxx:434-440 · gp_Dir2d.cxx:67-77 · gp_Trsf2d.cxx:31-46 · gp_Pnt2d.cxx:52-79 · §9.615 · gp_Lin2d.hxx:142-186 · gp_Ax2d.hxx:155-186 · gp_Ax2d.cxx:51-77 · gp_Dir2d.cxx:108-118 · gp_Elips2d.hxx:245 :355-370 :418-428 · gp_Elips2d.cxx:65-80 · gp_Hypr2d.hxx:500-513 · gp_Parab2d.hxx:313-321 · §9.616 · gp_Dir.hxx:496-501 · gp_Dir.cxx:86-155 · gp_Ax1.hxx:163-211 · gp_Ax1.cxx:46-76 · gp_Ax2.hxx:320-354 · gp_Ax2.cxx:83-124 · gp_Ax3.hxx:265-325 · gp_Ax3.cxx:82-115 · gp_Circ.hxx:236-266 · gp_Elips.hxx:341-394 · gp_Cylinder.hxx:203-237 · gp_Sphere.hxx:208-242 · gp_Trsf.cxx:57-86 :159-168 :172-194 :217-240 :387-390 :397-425 · gp_Trsf.hxx:389-396 :420-431 :464-476 · gp_Mat.cxx:30-43 :102-118 :122-160 · gp_Mat.hxx:246-254 · §9.618 · gp_Trsf.cxx:111-156 :243-280 :284-337 :346-385 :564-708 :712-836 :862-938 · gp_Mat.cxx:45-100 :292-334 · gp_Mat.hxx:55-107 :127-131 · Standard_Integer.hxx:39-42 · §9.619 · gp_Quaternion.cxx:81-88 :91-112 :153-187 :226-296 :302-357 :362-412 :416-444 · gp_Quaternion.hxx:106-129 · gp_EulerSequence.hxx:37-73 · gp_Mat2d.cxx:40-89 :144-179 · gp_Mat2d.hxx · gp_Trsf2d.cxx:72-83 :86-117 :120-196 :198-213 :216-219 :384-548 :550-674 :676-710 :721-746 · §9.620 · gp_GTrsf.cxx:28-197 :202 · gp_GTrsf.hxx:63-91 :118-136 :140-174 :187-213 :316-343 :347-386 :390-428 · gp_GTrsf2d.cxx:24-204 · gp_GTrsf2d.hxx:54-59 :83-131 :144 :166 :193 :235-283 :287-305 · §9.621 · gp_Trsf.hxx:328-352 :355 :358 · gp_GTrsf.hxx:260-303 :305 · gp_Trsf.cxx:943-960 :964-1009 · gp_Mat.cxx:338 · gp_Mat.hxx:301 · gp_GTrsf.cxx:202 · §9.622 · gp_Trsf.cxx:90-99 :103-109 :207-213 · §9.623 · gp_Mat.cxx:242-279 :283-289 :292-334 · gp_XYZ.hxx:367-373 · gp_Trsf.hxx:254 :464-482 · gp_Trsf2d.cxx:198-213 · gp_Dir.cxx:129-155 · gp_Dir.hxx:496-501 · gp_Dir2d.cxx:80-105 · gp_Vec.cxx:120-136 · gp_Vec2d.cxx:112-130 · gp_GTrsf.cxx:44-58 · §9.611-§9.623",
      dependsOn: ["T-105"],
      evidence: "**2026-10-05，§9.598 实测（本机）**：① 拓扑一致：端口 `TOTAL faces=223`、OCCT `FSTAT_O…",
      note: "编号说明：用户本轮任务编号为 T-104，但板上 `T-104` 已被 §9.588（`BRepLi…",
    },
        { id: "T-107",
      lane: "port-gap",
      title: "`ShapeFix_Face::FixOrientation` 单环分支落地（§9.625）：top.step 环面面片 112 → 648 与 OCCT 逐位相同，导出三角 20042 → 22794（参考 22800）",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 3,
      actual: 2,
      startedAt: "2026-10-07",
      updatedAt: "2026-10-07",
      completedAt: "2026-10-07",
      blocker: "",
      goal: "`ShapeFix_Face::Perform` 第二部分循环（`ShapeFix_Face.cxx:687-697`）对 `myResult` 每个面在 `NeedFix(myFixOrientationMode)` 时跑 `FixOrientation(MapWires)`，重建面经 `Context()->Replace`（`cxx:1639-1641`）由 `ShapeFix_Shape` 的 `Apply`（`ShapeFix_Shape.cxx:257`）落地。端口整段 UNPORTED，STEP `.F.` bound 的单环面 wire 朝向保持与自然外环相反，网格侧边界环定向为负、Delaunay 前沿饿死。",
      next: "旁支（未动，见 §9.625-D）：(1) top.step 端口面 321 ↔ OCCT 面 177 的 bbox 各差 0.382684（端口更小）⇒ 两块不是同一面片，`407/720` 与 `370/648` 不同口径，该 bbox 差开关两侧逐位相同、属既有裁边范围差；(2) 其余 41/324 面按 `--fstats` 差 ±1..4 节点/三角（累计 +1.2%），属边界采样差；(3) `-51.85…`/`-103.04…` 两处「同 bbox 面密度差一倍」类别未变。",
      acceptance: "**已满足**：`cargo check --manifest-path crates/occt-topo/Cargo.toml --all-targets` = 0 error；`--test step_obj_gates` = 4 passed；`export_data_obj -- data/occ data` = 23 ok / 0 err；开关对拍（只翻 `fix_orientation_mode`）top `15390/20042 → 16902/22794`、bottom `17190/22744 → 17586/23424`、T0M `60644/67532 → 60640/67524`，其余 5 模型与 `data/` 15 个 `*.step` 逐字节不变。",
      write: "crates/occt-topo/src/shhealing/shape_fix_face.rs · crates/occt-topo/src/step/read_topology.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.625 · ShapeFix_Face.cxx:141 :143 :218 :687-697 :1120-1160 :1165-1275 :1610-1646 · ShapeFix_Root.lxx:101 · ShapeProcess_OperLibrary.cxx:828 · ShapeFix_Shape.cxx:257 · ShapeFix_Shell.cxx:108 :202",
      dependsOn: ["T-105"],
      evidence: "**2026-10-07，§9.625 实测（本机）**：① 逐面（`zz_probe_a3n00 data/occ/top.step --fstats` ↔ OCCT `--facestats 0.451774 0.349066`）改动面恰好 5 个：`face=51/52/53/54` `75/112 → 370/648`（与 OCCT 同 bbox 面的 370/648 逐位相同）、`face=321` `75/112 → 407/720`；合计增量 `4×536+608=2752` 与总数 `22794−20042` 相符；按 (bbox, nodes/triangles) 全同的面 278 → 282 / 324。② 表面积（`obj_area` 口径）：top 1.00000（52132.94 / 52132.85）、motoc 1.00000、a3n00 0.99999、ATU01038 0.99999、acs10 1.00004、T0M 0.99998、TDB 0.99978、bottom 0.97412（仍缺 2.6%）。",
      note: "环面面片只占 top 总面积 ~0.7%，修前缺 12% 三角时面积仍在 1% 门限内 ⇒ 抓到该缺陷的是逐面 `--fstats` 计数。",
    },
        { id: "T-108",
      lane: "port-gap",
      title: "`ShapeFix_ComposeShell::MakeFacesOnPatch` 全段 + `CopyNMVertex` 面重载落地（§9.626）：T0M/a3n00/acs10 逐面 `--fstats` 与 §9.625 逐位相同",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 4,
      actual: 3,
      startedAt: "2026-10-08",
      updatedAt: "2026-10-08",
      completedAt: "2026-10-08",
      blocker: "",
      goal: "`MakeFacesOnPatch`（`ShapeFix_ComposeShell.cxx:2978-3271`）是 `DispatchWires`（`cxx:3582`）拼装输出面的最后一步。端口四处偏离 `.cxx`：单环 `myInvertEdgeStatus` 分支未接 `FixOrientation`；切向游走用 `loop {}`（不复查 `stPoint`、`aCL` 恒为第 k 条边 last、`ew` 位置不跨 `j` 持久）；wire 孔取「首条定向边」而非 `TopoDS_Iterator` 首子边；顶点孔所需的 `CopyNMVertex(V, toFace, fromFace)`（`Proj.cxx:715-802`）缺失。",
      next: "旁支（未动，见 §9.626-D）：(1) `CopyNMVertex` 边重载（`Proj.cxx:562-710`）未移植 —— 它只在 `ShapeFix_Wireframe::FixWireGaps` → `ShapeFix_Wire::FixGaps3d/FixGaps2d` → `FixGap3d/FixGap2d`（`ShapeFix_Wire_1.cxx:76-98 :101-123 :154-862 :893-2005`）与 `ShapeUpgrade_WireDivide` 下可达，两条链都未移植、也不在 STEP `FromSTEP` 序列（只跑 `FixShape`，`STEPControl_Controller.cxx:201`）；等该链移植时一并接；(2) `wire_fix.rs` 3379 行 > 1000，应按仓库规则拆 `shhealing/wire_fix/` 子模块；(3) §9.625 旁支 1/2/3（top.step 面 321 bbox 差、±1..4 节点采样差）仍在。",
      acceptance: "**已满足**：`cargo check --all-targets` = 0 error；`--test step_obj_gates` = 4 passed（23 模型表含 8 个 `data/occ/*`）；`export_data_obj` = 15 ok / 0 err；`MFP_PROBE` 触发面实测 T0M 226 调用 / 13 多环、a3n00 47/1、acs10 114/4（顶点孔与单环 invert 均 0 触发）；逐面 `--fstats` 合计与 §9.625 逐位相同：a3n00 11138/12440、acs10 37291/46548、T0M 60640/67524。",
      write: "crates/occt-topo/src/shape_fix_compose_shell/make_faces_on_patch.rs · crates/occt-topo/src/shhealing/transfer_params.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.626 · ShapeFix_ComposeShell.cxx:2978-3271 :3034 :3104-3130 :3197-3209 :3238-3245 :3582 · ShapeAnalysis_TransferParametersProj.cxx:562-710 :715-802 · ShapeFix_Root.cxx:26-30 · ShapeFix_Face.cxx:120-128 :1639-1641",
      dependsOn: ["T-107"],
      evidence: "**2026-10-08，§9.626 实测（本机）**：① 触发面（临时 `MFP_PROBE`，已撤除）：T0M.stp 226 次调用、13 次多环（9×holes=0 rev=false、2×holes=0 rev=true、1×holes=1、1×holes=7）；a3n00.stp 47/1（holes=1）；acs10.stp 114/4（2×holes=0、1×holes=1、1×holes=7）；top/bottom/TDB/motoc/ATU01038 0 次。② 逐面 `--fstats`（`zz_probe_a3n00 --fstats`）合计：a3n00 11138/12440、acs10 37291/46548、T0M 60640/67524（面数 226/787/1778），与 §9.625 基线逐位相同。③ `export_data_obj` 15 模型读数与 §9.625 全同。",
      note: "多环路径（根识别/切向游走/孔分类）在 T0M/a3n00/acs10 上真实触发 ⇒ 本轮不是纯潜伏改动；顶点孔与单环 `invert` 分支在 8 个 occ 模型上 0 触发（潜伏分支）。",
    },
        { id: "T-109",
      lane: "port-gap",
      title: "`IntCurve_IntConicConic` (Circle,Circle) 专用重载 + (Circle,Line)/(Ellipse,Line) 反向行落地（§9.628）：`ginter.rs` 里 T-100 登记的「真实缺口」清零",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 3,
      actual: 2,
      startedAt: "2026-10-08",
      updatedAt: "2026-10-08",
      completedAt: "2026-10-08",
      blocker: "",
      goal: "`Geom2dInt_GInter::InternalPerform` 的 (Circle,Circle) 行（`gxx:360-374`）应走 `IntConicConic` 专用重载，(Circle,Line)/(Ellipse,Line) 行（`gxx:346-358`、`:441-455`）应 `SetReversedParameters(true)` 后以 `Perform(Line(C2), D2, conic(C1), D1)` 调用；端口此前分别落在通用 `IntConicCurveGen` 臂与 `default:` 臂，且 `IntCurve_IntConicConic_1.cxx:58-151 :154-355 :807-1206`（`ProjectOnC2AndIntersectWithC2Domain` / `CircleCircleGeometricIntersection` / `Perform(Circ2d,Circ2d)`）整段 UNPORTED。这是 `ShapeAnalysis_Wire::CheckSelfIntersectingEdge` / `ShapeFix_Wire::FixSelfIntersectingEdge` 的公共前置。",
      next: "旁支（未动，见 §9.628-C/D）：(1) 仍走通用臂的行 —— Circle/Ellipse(`gxx:377-388`)、Circle/Parabola(`:390-404`)、Circle/Hyperbola(`:406-419`)、Ellipse/Circle 反向(`:458-471`)、Ellipse/Ellipse(`:473-486`)、Ellipse/Parabola(`:488-502`)、Ellipse/Hyperbola(`:504-518`) 及两行 `default:`(`:421-433`、`:517-529`)；属特化非新能力（通用臂的 `IntCurveIConicTool` 已能表示所有圆锥曲线），仅在重合/切线退化上不如专用重载完整；(2) `Extrema_ExtElC2d(gp_Lin2d, gp_Elips2d)` 仍未移植 ⇒ `line_ellipse_geometric_intersection` 的 `D < 0` 切线支仍见 §9.627-B 末；(3) `ShapeAnalysis_Wire::CheckSelfIntersectingEdge` / `CheckIntersectingEdges`（进而 `FixSelfIntersectingEdge` / `FixIntersectingEdges`）仍待移植，本批已解除其前置。",
      acceptance: "**已满足**：`cargo check --manifest-path crates/occt-geom2d/Cargo.toml --all-targets` = 0 error；`cargo check --manifest-path crates/occt-topo/Cargo.toml` = 0 error；`--test step_obj_gates` 4 passed；`export_data_obj` 15 ok / 0 err（与 §9.626 基线相同）。",
      write: "crates/occt-geom2d/src/geom2d_int/int_conic_conic_circle_circle.rs · crates/occt-geom2d/src/geom2d_int/ginter.rs · crates/occt-geom2d/src/geom2d_int/mod.rs · crates/occt-core/src/gp/circ2d.rs",
      ref: "specs/_a3n00_gap_analysis.md §9.628 · IntCurve_IntConicConic_1.cxx:58-151 :154-355 :807-1206 :850 :869 · Geom2dInt_Geom2dGInter.gxx:346-358 :360-374 :441-455 · gp_Circ2d.hxx Reversed",
      dependsOn: ["T-100"],
      evidence: "**2026-10-08，§9.628 实测（本机）**：① 逐段对 `_1.cxx`：`project_on_c2_and_intersect_with_c2_domain`（`:58-151`）、`circle_circle_geometric_intersection`（`:154-355`，`Tol` 管 + `TolTang` 切线管、`nbsol`=0/1/2/3、区间 `> PI` 归一化取补）、`perform_circle_circle`（`:807-1206`，反向/间接圆转正向并重映射域界、周期归一化、端点裁剪、`D2` 导数与切/法向量、`determine_transition_lc`）、`next_after_pi_p_pi`（`:850 :869`）。② `gp_Circ2d::reversed()` 对 `gp_Circ2d::Reversed()`：`location`/`xdir`/`radius` 不变，只翻局部 `GpAx22d` 的 `vydir`。③ 分派臂：`typ1 ∈ {Circle, Ellipse} × typ2 == Line` → `set_reversed_parameters(true)` + `perform_line_circle`/`perform_line_ellipse`；`Circle × Circle` → `set_reversed_parameters(false)` + `perform_circle_circle`。④ **触发面（可达性）**：全仓唯一实代码调用点是 `split_by_line.rs:204`（`ShapeFix_ComposeShell::SplitByLine`，`cxx:1605-1656`），其第一参数 `j_c2d` 是 `Geom2dLine::new(*line.position())`（`cxx:1443-1444`）⇒ **`typ1` 恒为 Line**，本批三条臂当前不可达、属潜伏分支；可达条件是 `ShapeAnalysis_Wire::CheckSelfIntersectingEdge` / `CheckIntersectingEdges`（经 `Geom2dAPI_InterCurveCurve` / `Geom2dInt_GInter::Perform(C1, C2)` 传两条边自身曲线）移植。⑤ 门禁：`--all-targets` 0 error；`step_obj_gates` 4 passed（571.45s）；`export_data_obj` 15 ok / 0 err（读数与 §9.626 相同）。",
      note: "本批只做最小必要集（(Circle,Circle) + 两条反向 Line 行），因为目标是解除自交修复链前置而非补齐全部专用重载；其余行数学等价地由通用臂处理，见 §9.628-D。**与 T-108 同类：新增臂在本仓当前是潜伏分支（0 触发）**，故行为与基线逐位一致 —— 记录以免把「接线完成」误读为「路径已激活」。",
    },
        { id: "T-111",
      lane: "port-gap",
      title: "IGES 写侧 124 变换矩阵未走 `gp_Trsf::SetTransformation`：`IGESConvGeom_GeomBuilder` 的符号零归一化缺失 + 球/环母线 frame 由字面量 `(0,-1,0)` 构造",
      status: "completed",
      priority: "P1",
      owner: "agent",
      progress: 100,
      estimate: 4,
      actual: 4,
      startedAt: "2026-10-08",
      updatedAt: "2026-10-08",
      completedAt: "2026-10-08",
      blocker: "",
      goal: "初等面（柱/锥/球/环）在本写侧已走 `GeomToIGES_GeomSurface.cxx:696-995` 的 **120 回转面 + 124** 路径（`emit_local_revolution_surface`），但 124 的数值与 oracle 不同：`emit_transformation_matrix` 直接抄 frame 的三个方向分量，而 OCCT 经 `IGESConvGeom_GeomBuilder::SetPosition(gp_Ax3)` → `gp_Trsf::SetTransformation(pos, gp::XOY())` → `MakeTransformation`；`gp_Trsf::SetTransformation` 末尾的 `matrix.Multiply(MA1)`（单位阵 × 由 frame 方向作列构成的 `gp_Mat`）会按 IEEE 求和把多数 `-0.` 归一成 `+0.`。此外球/环母线 frame 由字面量 `GpDir::new(0,-1,0)` 构造，丢了 `-gp::DY()` 经 `gp_Dir::Reverse` 得到的 `-0.` 零符号。",
      next: "旁支（未动，见本轮量化）：逐实体多重集（去 `-0.`）仍有 **357** 条真实差，按 OCCT 类型分：126×81（T-117：Bezier 分段/有理臂）、110×42、120×41（面 UV 上界差，如 Cone `2π-U2=-4.14E-13` 端口取 `0.`，属上游 `UVBounds` 缺口）、102×36、100×29、108×26、128×15、124×11、142×8、140×3、402×2。这些不在本批范围内，另开任务。",
      acceptance: "**已满足**（现有基线，不新写测试）：① `cargo check --manifest-path crates/occt-topo/Cargo.toml --all-targets` = 0 error；② `--example iges_check -- <23 模型>`（15 个 `data/*.step` + 8 个 `data/occ/*`）= **23/23 `ok` / 0 problem**；③ `Cylinder/Cone/Sphere/Torus` 的 DE 卡数与 oracle **逐位相同**（28 / 20 / 14 / 20）；④ `--test phase5_integration` = **7 passed**（含 `real_brepmesh_and_iges`）；⑤ `_igesdump/oracle` vs `_igesdump/port` 的 14 模型逐实体多重集差：**486 → 374**（去 `-0.` 归一后 **378 → 357**），其中 **124 类型差 123 → 11**。",
      write: "crates/occt-topo/src/iges.rs",
      ref: "GeomToIGES_GeomSurface.cxx:696-768 :796-860 :879-925 :947-995 · IGESConvGeom_GeomBuilder.cxx:137-143 :218-237 · gp_Trsf.cxx:172-194 · gp_Mat.hxx:380-409 · gp_XYZ.hxx:550-565 :580-594 · gp_Ax2.hxx:73-80 · gp_Dir.hxx:440-455 · gp_Ax3.hxx:477-482",
      dependsOn: [],
      evidence: "**2026-10-08 实测（本机，`crates/occt-topo/src/iges.rs`）**：① `emit_transformation_matrix` 改为 `GpTrsf::set_transformation_from_to(frame, xoy_frame())` 后读 `Value(i,j)`（i=1..3, j=1..4，第 4 列除 `unit`），与 `IGESConvGeom_GeomBuilder::SetPosition/MakeTransformation` 逐行对应；form 用 `GpTrsf::is_negative()`（= `thepos.IsNegative()`）。② `revolution_generatrix_frame` 改为 `GpDir::from_axis(Y).reversed()`（= `-gp::DY()`，零带负号）+ `n ^ (vx ^ n)` 正交化（= `gp_Ax2` 3 参 ctor 的 `CrossCross`），再交给 `GpAx3::new` 算 `vydir = n ^ vxdir`。③ 复跑 `IGES_DUMP` 重导出 14 模型，逐实体多重集（按 `类型|P 文本` 计数）差：总 **486 → 374**；去 `-0.` 后 **378 → 357**；124 类 **123 → 11**。④ `iges_check` 23 模型全 `ok`（T0M 85023 DE / TDB 105732 DE / acs10 36723 DE），`phase5_integration` 7 passed。",
      note: "关键机理（都在 `.cxx` 有对应控制流，非调参）：`gp_Mat::Multiply`（`gp_Mat.hxx:380-409`）对每个分量按 `a*b + c*d + e*f` 三次求和，故单位阵乘 `MA1` 时 `0.0*(-0.0) = -0.0` 与 `+0.0 + (-0.0) = +0.0` 会把多数 `-0.` 归一；`-gp::DY()` 的 `(-0.,-1.,-0.)` 与字面量 `(0.,-1.,0.)` 差别只在两处零符号，但足以让 124 逐字符不同。本批**未**改任何面/边 UV 取值，故 120/110/100 的真实差（面 bounds 差）保留在旁支。",
    },
        { id: "T-120",
      lane: "port-gap",
      title: "IGES 写侧 DE 状态字段（`SubordinateStatus`/`UseFlag`）未移植：`IGESControl_Writer::ComputeModel` 的 `ComputeStatus` + `AutoCorrectModel` 两步缺失",
      status: "completed",
      priority: "P2",
      owner: "agent",
      progress: 100,
      estimate: 3,
      actual: 3,
      startedAt: "2026-10-08",
      updatedAt: "2026-10-08",
      completedAt: "2026-10-08",
      blocker: "",
      goal: "`IGESControl_Writer::ComputeModel`（`IGESControl_Writer.cxx:256-262`）在写文件前跑 `myEditor.ComputeStatus()` 与 `myEditor.AutoCorrectModel()`，端口 `write_shape_iges` 无这两步：`Ent::directory_lines` 把 DE card 1 第 9..12 字段（8 列）写成 `00`×4。oracle（`IGESControl_Controller::Init` + `IGESControl_Writer(\"MM\",0)` + `AddShape` + `ComputeModel` + `Write`，`specs/occt_probe --iges`）实测该 8 列只有 4 种取值：`00|00|00|00` 361 张、`00|01|00|00` 2090、`00|02|00|00` 203、`00|01|05|00` 218（Blank 与 Hierarchy 恒 0）。",
      next: "旁支（未动）：(1) DE 逐号对齐仍受实体编号顺序差影响 —— 本轮 163 张状态不同的 DE **全部**同时类型不同（Offset 147/147、Shape-1 4/4、HoledPlate 8/8、Cylinder 2/2、Torus 2/2），编号顺序（`AddWithRefs` DFS 序）本身是独立缺口；(2) `ComputeStatus` 的 UseFlag 传播臂（2xx/134/116/132）与 `Interface_Graph` 状态数组未移植（见 note）；(3) 其余 P 文本差 1250 条属几何差异（T-111/T-113/T-114 等）。",
      acceptance: "**已满足**（现有基线，不新写测试）：① `--example iges_check -- <23 模型>`（15 个 `data/*.step` + 8 个 `data/occ/*`）= **23/23 `ok` / 0 problem**；② `_igesdump/entside.ps1` 口径的逐实体 P 文本差集 = **1250**（与上一轮记录逐位相同，本批只动 DE 列）；③ DE card 1 第 9..12 字段：`_igesdump/port/*.iges` 与 `_igesdump/oracle/*.iges` 的**全局分布逐位相同**（361 / 2090 / 218 / 203），逐号比对 2872 张中仅 163 张不同且**全部**落在类型也不同的 DE 上（编号顺序差），类型相同的模型（Cube、Extrusion、Sphere、rev、screw、Shape、Shape-2、OffsetPlaneHoleEdge、Cone 等）差异为 **0**（上一轮对照：Sphere 10/14、Cube 48/49、Shape-1 819/1034）。",
      write: "crates/occt-topo/src/iges.rs",
      ref: "IGESControl_Writer.cxx:256-262 · IGESData_BasicEditor.cxx:216-328 :330-433 · IGESData_DirChecker.cxx:33-36 :361-453 · IGESGeom_ToolCurveOnSurface.cxx:144-162 :182-211 · IGESGeom_ToolTrimmedSurface.cxx:264-276 · IGESGeom_SpecificModule.cxx:313-379 · IGESData_IGESWriter.cxx:343-360 :836-840 · Interface_Graph.cxx:248-255 :369-398",
      dependsOn: [],
      evidence: "**2026-10-08 实测（本机，`crates/occt-topo/src/iges.rs` 新增 `compute_status`）**：① 全局 DE card 1 第 9..12 字段分布 oracle = port = `00|00|00|00` 361 / `00|01|00|00` 2090 / `00|01|05|00` 218 / `00|02|00|00` 203。② 逐号（按 DE 号配对）比对 14 个 oracle 模型：matched 2872，状态不同 163 = 类型也不同的 163（Shape-1 4、HoledPlate 8、Cylinder 2、Torus 2、Offset 147），故差异全部来自既有实体编号顺序缺口、非状态计算。③ `iges_check` 23/23 `ok`（含 T0M 85023 DE / TDB 105753 DE / acs10 36723 DE），0 problem。④ P 文本逐实体差集 = 1250（Cube 0、Extrusion 0、Cone 3、Sphere 3、rev 8、Cylinder 12、Torus 15、OffsetPlaneHoleEdge 40、screw 45、Shape 49、Shape-2 128、HoledPlate 192、Offset 353、Shape-1 402），与上一轮逐位相同 ⇒ 本批未触碰 P 段。⑤ `--test phase5_integration` = 7 passed（含 `real_brepmesh_and_iges`）。",
      note: "实现要点（都在 `.cxx` 有对应控制流）：① Subordinate —— 逐实体对 `Ent::refs`（= `OwnSharedCase`）的每个后代 `|= 1`，引用方为 402/404 时 `|= 2`（`BasicEditor.cxx:257-281`），随后 `InitStatus(bl, subs, uf, hy)` 写回（`:315-327`）；DE 部分引用（`Ent::trsf`，field 7）不属 `OwnShared`，故 124 恒 `sub=0`，与 oracle 逐条相符。② UseFlag 传播臂（`:283-309`）**未移植**：只对 2xx/134/116/132 触发，本写侧（100/102/108/110/120/124/126/128/140/142/144/402）不产生这些类型，且需 `Interface_Graph` 的状态数组（其 shared 集含 DE 指针，本写侧不建模）⇒ `G.Status(i)` 恒 0、`uf==0` 时保持 0。③ AutoCorrect 只移植到 `DirChecker::Correct`：本写侧类型中仅 142 声明 `UseFlagRequired(5)`（`IGESGeom_ToolCurveOnSurface.cxx:209`）、144 声明 `UseFlagRequired(0)`（`IGESGeom_ToolTrimmedSurface.cxx:272`），其余字段留在 `-100` 默认「不检查」（`IGESData_DirChecker.cxx:33-36`，`thegraphier = -100` 同时使 DE 图形清理段失效，`:400-420`）。**修正上一轮的理解**：oracle 的 `00|01|05|00` 218 张来自 `DirChecker::Correct` 对 142 自身置 uf=5，**不是** `OwnCorrect` 对 `CurveUV` 置 5 —— OCCT 8.0.0 中 142 的 `OwnCorrect` 注册在 case 9 但 `DeclareAndCast` 成了 `IGESGeom_Boundary`（`IGESGeom_SpecificModule.cxx:337-345`），`anent.IsNull()` 直接返回，故 Sphere 里 142(#6) 的 `CurveUV` = 102(#7) 状态仍是 `00|01|00|00`；端口按此行为照抄，不实现该死分支。",
    },
  ],
  gates: [
        { id: "G-check", name: "五 crate 编译", cmd: "cargo check --all-targets", baseline: "exit 0", last: "exit 0", status: "green",
      acked: false, owner: "",
      note: "五 crate 逐 crate 跑（无根 workspace…",
    },
        { id: "G-lib", name: "topo 单测", cmd: "--lib", baseline: "1281 / 0", last: "1281 / 0", status: "green", acked: false, owner: "",
      note: "2026-10-07 三轮实测：`--lib` **1281 passed / 0 failed**…",
    },
        { id: "G-gates", name: "STEP→OBJ 门禁", cmd: "--test step_obj_gates", baseline: "4 / 4", last: "4 / 4", status: "green", acked: false,
      owner: "",
      note: "2026-10-08 §9.626 复跑 4 passed（375.81s 冷跑）；2026-10-07 三轮实测 4 passed（336.72s 热跑…",
    },
        { id: "G-fuse", name: "fuse 闭环实体", cmd: "fuse_box_cylinder_is_closed_solid", baseline: "1 / 1", last: "1 / 1", status: "green",
      acked: false, owner: "",
      note: "2026-10-06 实测 1 passed（1280 filtered out）…",
    },
        { id: "G-phase", name: "phase3-10 / 19 / 20", cmd: "--test phaseN_integration", baseline: "全绿", last: "全绿", status: "green",
      acked: false, owner: "",
      note: "4·9·7·5·5·5·8·8·5·5 = 61 passed（2026-10-07 二轮重跑…",
    },
        { id: "G-boss", name: "boss 合并", cmd: "--test bop_builder2_boss", baseline: "1 / 2", last: "1 / 2", status: "red", acked: false,
      owner: "T-03",
      note: "夹具不合法（定案 b），非引擎缺口 —— 已接受的红",
    },
        { id: "G-export", name: "导出总检", cmd: "--example export_data_obj", baseline: "15 / 15", last: "15 / 15", status: "green",
      acked: false, owner: "",
      note: "2026-10-08 §9.626 复跑 15 ok / 0 err；2026-10-06 收尾轮仍 15/15（§9.617-§9.623 每轮复跑）；2026-10-07 §9.625 复跑 `-- data/occ data` = 23 ok / 0 err…",
    },
        { id: "G-iges", name: "IGES 结构自洽", cmd: "--example iges_check -- <23 模型>", baseline: "23 / 23", last: "23 / 23", status: "green",
      acked: false, owner: "",
      note: "2026-10-08 复跑 23 / 23（三处写侧忠实化：P 卡按 `Interface_LineBuffer` 字段边界换卡、DE card 1 状态字段按 `%2.2d` 零填充、`ComputeStatus`+`AutoCorrectModel` 填状态 8 列；P 卡总数随之变化，如 Shape-2 7762 → 8521、T0M 302322、TDB 305438；DE 状态 8 列分布现与 oracle 逐位相同 361/2090/218/203；2026-10-07 二轮把集合扩到全部 23 个模型（15 个 data/*.step + 8 个 data/occ/*）；2026-10-08 T-111 把 124 改走 `GpTrsf::set_transformation_from_to`（符号零归一）后复跑 23/23，TDB DE 105753 → 105732）…",
    },
        { id: "G-core", name: "core / math / geom / geom2d", cmd: "--lib", baseline: "290·215·143·72", last: "290·215·143·72", status: "green",
      acked: false, owner: "",
      note: "逐项与基线相同（2026-10-06 §9.619/§9.622 复跑 occt-core / oc…",
    },
  ],
  // 新条目插到数组**开头**（brief 取前 3 条当最近活动；早于 a29 的看提交历史）
  activity: [
    {
      id: "a89",
      at: "2026-10-08",
      title: "T-109 闭环：`IntCurve_IntConicConic` (Circle,Circle) 专用重载落地（§9.628）—— `ginter.rs` 里 T-100 登记的「真实缺口」清零，(Circle,Line)/(Ellipse,Line) 反向行一并补齐",
      tone: "success",
      detail: "§9.628。缺口对着 `.cxx`：`Geom2dInt_Geom2dGInter.gxx:360-374` 的 (Circle,Circle) 行应走 `IntConicConic`，端口此前落在通用 `IntConicCurveGen` 臂；`:346-358`(Circle/Line) 与 `:441-455`(Ellipse/Line) 两行要 `SetReversedParameters(true)` 后以 `Perform(Line(C2), D2, conic(C1), D1)` 调用，端口此前落在 `default:` 臂。本批新增 `crates/occt-geom2d/src/geom2d_int/int_conic_conic_circle_circle.rs`：`project_on_c2_and_intersect_with_c2_domain`（`_1.cxx:58-151`，把 C1 解经圆参数映射搬到 Circle2 与 `DomainC2` 求交、再把存活区间投影回 Circle1）、`circle_circle_geometric_intersection`（`_1.cxx:154-355`，`Tol` 管 + `TolTang` 切线管，`nbsol` = 0/1/2/3，3 为重合；区间 `Length() > PI` 时归一化并取补）、`IntCurveIntConicConic::perform_circle_circle`（`_1.cxx:807-1206`，反向/间接圆转正向并重映射域界、周期域归一化、端点裁剪、`D2` 导数与切/法向量、`determine_transition_lc`、填 `IntRes2dIntersectionPoint`/`Segment`）、`next_after_pi_p_pi`（`_1.cxx:850 :869` 的 `std::nextafter(PIpPI, 0.)`）；`crates/occt-core/src/gp/circ2d.rs` 补 `gp_Circ2d::reversed()`（对 `gp_Circ2d::Reversed()`：`location`/`xdir`/`radius` 不变，翻局部 `GpAx22d` 的 `vydir`）；`ginter.rs::internal_perform` 加三条分派臂。这是 `ShapeAnalysis_Wire::CheckSelfIntersectingEdge` / `ShapeFix_Wire::FixSelfIntersectingEdge` 的公共前置（§9.627 D.1/D.3）。**可达性实测：全仓唯一实代码调用点 `split_by_line.rs:204` 的第一参数是 `Geom2dLine`（`cxx:1443-1444`）⇒ 本批三条臂当前 0 触发、属潜伏分支**，行为与基线逐位一致（`step_obj_gates` 4 passed、`export_data_obj` 15 ok / 0 err）。" ,
      ref: "specs/_a3n00_gap_analysis.md §9.628 · IntCurve_IntConicConic_1.cxx:58-151 :154-355 :807-1206 :850 :869 · Geom2dInt_Geom2dGInter.gxx:346-358 :360-374 :441-455 · gp_Circ2d.hxx Reversed"
    },
    {
      id: "a88",
      at: "2026-10-08",
      title: "T-108 闭环：`ShapeFix_ComposeShell::MakeFacesOnPatch` 全段 + `CopyNMVertex` 面重载 —— 切向游走按 `.cxx` 复刻（跨 j 共享 `ew`），T0M/a3n00/acs10 逐面 `--fstats` 与 §9.625 逐位相同",
      tone: "success",
      detail: "§9.626。缺口对着 `.cxx`：单环 `myInvertEdgeStatus` 分支要 `new ShapeFix_Face` 跑 `FixOrientation`（`cxx:2997-3006`）；切向游走是 `while (stPoint == ON || UNKNOWN)` 沿 i 循环内、j 循环外的共享 `TopoDS_Iterator ew(wr)` 前进、`aCL` 随 `aCW` 更新（`cxx:3104-3130`）；wire 孔取 `TopoDS_Iterator(bw)` 首子边、不筛朝向（`cxx:3197-3209`）；顶点孔要 `CopyNMVertex(V, newFace, myFace)` + `Context()->Replace` + `B.Add`（`cxx:3238-3245`）。端口此前用 `loop {}`（不复查 `stPoint`、`aCL` 恒为第 k 条边 last、`ew` 每个 j 从 k 重来）、孔取首条定向边、`CopyNMVertex` 缺失。移植 `copy_nm_vertex_face`（`Proj.cxx:715-802`，`ValueOfUV` + `Gap` 容差放宽），边重载（`Proj.cxx:562-710`）因其唯一可达链（`ShapeFix_Wireframe::FixWireGaps` → `ShapeFix_Wire::FixGaps3d/2d` → `FixGap3d/2d`、`ShapeUpgrade_WireDivide`）全未移植且不在 STEP `FromSTEP` 序列、表示链亦无数据模型而登记 UNPORTED。触发面实测（`MFP_PROBE`，已撤除）：T0M 226 调用/13 多环、a3n00 47/1、acs10 114/4 ⇒ 多环路径真实触发。门禁：`--all-targets` 0 error、`step_obj_gates` 4/4、`export_data_obj` 15/15、逐面 `--fstats` 与 §9.625 逐位相同（a3n00 11138/12440、acs10 37291/46548、T0M 60640/67524）。",
      ref: "specs/_a3n00_gap_analysis.md §9.626 · ShapeFix_ComposeShell.cxx:2978-3271 :3034 :3104-3130 :3197-3209 :3238-3245 :3582 · ShapeAnalysis_TransferParametersProj.cxx:562-710 :715-802 · ShapeFix_Root.cxx:26-30 · ShapeFix_Face.cxx:120-128 :1639-1641"
    },
    {
      id: "a87",
      at: "2026-10-07",
      title: "T-107 闭环：`ShapeFix_Face::FixOrientation` 单环分支 —— top.step 4 个环面面片 112 → 648 与 OCCT 逐位相同，导出三角 20042 → 22794（参考 22800，−6）",
      tone: "success",
      detail: "§9.625。缺口对着 `.cxx`：`Perform` 第二部分循环（`ShapeFix_Face.cxx:687-697`）对 `myResult` 每个面在 `NeedFix(myFixOrientationMode)` 时跑 `FixOrientation(MapWires)`，重建面经 `Context()->Replace`（`cxx:1639-1641`）由 `ShapeFix_Shape` 的 `Apply`（`ShapeFix_Shape.cxx:257`）落地；端口整段 UNPORTED ⇒ STEP `.F.` bound 的单环面 wire 朝向未翻，网格边界环定向为负、Delaunay 前沿饿死。移植 `is_need_add_natural_bound`（`cxx:1120-1160`）+ `fix_orientation`（`cxx:1165-1275`、重建段 `cxx:1610-1646`），接入 `resolve_face` 两条返回路径。开关对拍（只翻 `fix_orientation_mode`）：top `15390/20042 → 16902/22794`、bottom `17190/22744 → 17586/23424`、T0M `60644/67532 → 60640/67524`，其余 5 模型与 `data/` 15 个 `*.step` 逐字节不变；逐面改动恰好 `face=51/52/53/54`（`75/112 → 370/648`，与 OCCT 同 bbox 面逐位相同）与 `face=321`（→ `407/720`）。",
      ref: "specs/_a3n00_gap_analysis.md §9.625 · ShapeFix_Face.cxx:141 :143 :218 :687-697 :1120-1160 :1165-1275 :1610-1646 · ShapeFix_Root.lxx:101 · ShapeProcess_OperLibrary.cxx:828 · ShapeFix_Shape.cxx:257"
    },
    {
      id: "a86",
      at: "2026-10-07",
      title: "A26 STEP 写侧忠实化：曲面 6×6 采样重拟换成面自身的极点/结点/权重导出（`-…",
      tone: "success",
      detail: "审计 A26 剩的最后一格（写侧脱离 OCCT）闭环。`step/write_context.rs` 里的 `fit_bspline_sur…",
      ref: "specs/_audit/_index.md A26 · specs/_audit/mesh-exchange.md 发现 15 · GeomToStep_MakeBoundedSurface.cxx:41-88 · GeomToStep_MakeBSplineSurfaceWithKnots.cxx:150-163 · data/occ/T0M.stp:901 · data/occ/ATU01038.step:31351-31372"
    },
    {
      id: "a85",
      at: "2026-10-06",
      title: "T-106 收口：`Geom2dConvert` 全链 + `Gp*` 基元变换正本清源（…",
      tone: "info",
      detail: "§9.604-§9.623 一串**潜伏缺口**闭环 —— 这些分支在 23 个基准模型上**不改变任何读数**（四模型 `--wirehi…",
      ref: "specs/_a3n00_gap_analysis.md §9.604 · §9.605 · §9.606 · §9.607 · §9.609-§9.623"
    },
    {
      id: "a84",
      at: "2026-10-06",
      title: "T-103 网格塌陷入口闭环：ON-only 塌陷面 8/8 归位 —— 2D pcurv…",
      tone: "success",
      detail: "§9.608。根因对着 `.cxx`：2D pcurve 求值用了**边的 3D range**…",
      ref: "specs/_a3n00_gap_analysis.md §9.603 · §9.608"
    },
    {
      id: "a83",
      at: "2026-10-05",
      title: "T-104 闭环：motoc 圆柱面「横向 2 层环带」= `add_wire` 里 OC…",
      tone: "success",
      detail: "§9.598。用户截图的蓝圈环带：按 bbox 配对 `TOTAL faces=223` 两侧一致 ⇒ 多的是**网格**不是面…",
      ref: "specs/_a3n00_gap_analysis.md §9.598 · §9.601"
    },
    {
      id: "a82",
      at: "2026-10-04",
      title: "T-105 闭环：移植 `FixPeriodicDegenerated` + `IsPer…",
      tone: "success",
      detail: "§9.590 把「空面」第一现场钉到 T0M face 1764：wire 塌成 `edges=[1] nv=1`、V 跨度 0 ⇒ `in…",
      ref: "specs/_a3n00_gap_analysis.md §9.590 · §9.591"
    },
    {
      id: "a81",
      at: "2026-10-04",
      title: "T-104：`BRepLib::BuildCurve3d` / `CheckSameRan…",
      tone: "success",
      detail: "§9.588。旧端口在 `FixAddCurve3d` 里用 `meshing::edge_discret::CurveOnSurface`…",
      ref: "specs/_a3n00_gap_analysis.md §9.588"
    },
    {
      id: "a80",
      at: "2026-10-04",
      title: "T-102 收口：T0M 面/wire 结构与 OCCT 逐项相同（1778 面 / 19…",
      tone: "success",
      detail: "§9.587。根因对着 cxx：`ShapeFix_ComposeShell.cxx:629` 的 `TopoDS_Iterator aIt…",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.587"
    },
    {
      id: "a79",
      at: "2026-10-02",
      title: "T-101 收口：a3n00 面积比 0.9963 → 1.0000、WIREHIST 与…",
      tone: "success",
      detail: "§9.580。两处新落地的 .cxx 分支：③ `wire_fix.rs::fix_dummy_seam` 改用 `ShapeAnalysi…",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.580"
    },
    {
      id: "a78",
      at: "2026-10-02",
      title: "T-101：F113 恢复为 1 面 / 2 wire（22+6）…",
      tone: "success",
      detail: "§9.579。两处 .cxx 分歧同时在场才成立：① `wire_fix.rs::copy_replace_vertices_with` 把…",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.579"
    },
    {
      id: "a77",
      at: "2026-10-02",
      title: "T-101：goal（160 轮）终态 —— 全链判定忠实、差异锁到「段集合（BREAK-…",
      tone: "info",
      detail: "§9.493-§9.578 逐段对拍，把差异逼到一点：93 段「首末顶点同一、逐边 LastVertex 是另一 TShape」（edge_…",
      ref: "specs/_a3n00_gap_analysis.md §9.578"
    },
    {
      id: "a76",
      at: "2026-10-02",
      title: "T-101：切分侧与 WireSegment 访问器逐条对照完成（全部忠实）",
      tone: "info",
      detail: "§9.534-§9.536 SplitByGrid（UVBounds/Bounds/TOLINT/closed 段位移/U 线循环）逐条一致…",
      ref: "specs/_a3n00_gap_analysis.md §9.577"
    },
    {
      id: "a75",
      at: "2026-10-02",
      title: "T-101：IsShortSegment 两处朝向语义修正落地（门禁 5/5、基线逐字不动…",
      tone: "info",
      detail: "§9.489-§9.492：对拍 OCCT IsShortSegment(ShapeFix_ComposeShell.cxx:2394-24…",
      ref: "specs/_a3n00_gap_analysis.md \u00a79.492",
    },
  ],
};

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

const num = (v: unknown): number => (typeof v === "number" && Number.isFinite(v) ? v : 0);

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
  return task.status === "completed" ? 100 : task.status === "cancelled" ? 0 : num(task.progress);
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
  const [activeId, setActiveId] = useCanvasState<string>("active", "T-103");
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
  const estTotal = counted.reduce((sum, task) => sum + num(task.estimate), 0);
  const weighted = counted.reduce((sum, task) => sum + num(task.estimate) * progressOf(task), 0);
  const progressPct = estTotal === 0 ? 0 : Math.round(weighted / estTotal);
  const overrun = counted.filter((task) => num(task.actual) > num(task.estimate));
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
  if (overrun.length > 0) risks.push(overrun.length + " 条实际已超出估算：" + overrun.map((task) => task.id + "（" + num(task.actual) + "/" + num(task.estimate) + "）").join("、"));
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
            : tasks.filter((task) => task.lane === filter && isOpen(task));

  const byLane = DATA.lanes.map((lane) => {
    const rows = tasks.filter((task) => task.lane === lane);
    const laneCounted = rows.filter((task) => task.status !== "cancelled");
    const laneEst = laneCounted.reduce((sum, task) => sum + num(task.estimate), 0);
    const laneWeighted = laneCounted.reduce((sum, task) => sum + num(task.estimate) * progressOf(task), 0);
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
      "\n当前：" + STATUS_LABEL[task.status] + "，进度 " + progressOf(task) + "%，估算 " + num(task.estimate) + " 人日 / 已用 " + num(task.actual) + " 人日" +
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
                  label={"进度 · 已用 " + num(active.actual) + " / 估算 " + num(active.estimate) + " 人日"}
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
          num(task.actual) + "/" + num(task.estimate),
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

      <CollapsibleSection title="口径与节奏（intake 决议）" count={9}>
        <Stack gap={6}>
          <Text size="small">· 决定：下一批做哪条 + 门禁能否提交 + 哪些红是已知非缺陷。</Text>
          <Text size="small">· 行单位：任务 <Code>T-xx</Code>（lane = 模块/族，共 {DATA.lanes.length} 组）；门禁是另一段，行单位不同。</Text>
          <Text size="small">· 状态：pending / in_progress / blocked / completed / cancelled；blocked 必须写明在等谁（本板当前无 blocked，历史只有 T-54）。</Text>
          <Text size="small">· 字段：状态 · 优先级 · owner · progress/estimate/actual · startedAt/updatedAt/completedAt · blocker · next · acceptance · evidence · write · ref · dependsOn。</Text>
          <Text size="small">· 派生（全部现算，不手抄）：加权进度 · WIP 超限 · 阻塞 · 陈旧 · 超估算 · 门禁待处理 · 依赖未完成即开工 · 在办共用改动位置 · 分组完成度 · 下一步 3 条。</Text>
          <Text size="small">· 人的动作：标记进行中/完成/阻塞、豁免、认领、还原、门禁标「已接受」、打开文件、交给 agent —— 全部走 overlay。</Text>
          <Text size="small">· 更新节奏：agent 每完成一批就更新 <Code>DATA</Code>（进度、更新日、证据）；人只动状态/认领/接受。</Text>
          <Text size="small" tone="tertiary">· 归档：completed / cancelled 超过一个窗口就移出当前窗口（见「已归档」折叠区），完整证据留在 git 与 specs/_a3n00_gap_analysis.md。</Text>
          <Text size="small" tone="tertiary">· 本轮门禁（2026-10-07 三轮，A26 STEP 写侧忠实化）：五 crate <Code>cargo check --all-targets</Code> 0 error、occt-topo <Code>--lib</Code> 1281/0、<Code>step_obj_gates</Code> 4/4（336.72s）；<Code>fuse_box_cylinder_is_closed_solid</Code> 1/1、<Code>export_data_obj</Code> 15/15、<Code>iges_check</Code> 23/23 为上一轮（2026-10-06/07 二轮）实测值、本轮未复跑。上一轮回执（2026-10-06 收尾轮，§9.617-§9.623）：四模型 <Code>--wirehist</Code> 逐项同基线；T-103 网格线临时仪器已全部撤除。</Text>
        </Stack>
      </CollapsibleSection>

      <Text size="small" tone="tertiary">
        筛选器与当前选中项只存在本机（<Code>useCanvasState</Code>），agent 看不到；状态与门禁接受落 sidecar，agent 可见。明细表的「证据」列指向 文件:行 / 测试名 / 提交号，结论性描述一律用引用。
      </Text>
    </Stack>
  );
}
