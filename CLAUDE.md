## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).

## OCCT 移植门禁

移植 OCCT 时按这五道门禁推进。对齐靠对拍现有基线，不靠堆断言。

1. **改完先编译。** 先 `cargo check`（或该 crate 的编译）。编译不过就修编译。
2. **对着 `.cxx` 审控制流再改代码。** 有缺口时先定位缺了哪段 OCCT 控制流，修复走同一套 `.cxx` 分支，不为某一个 STEP / 特征 / 特例调参（网格管线是共享的，特例补丁会带偏已对齐模型）。有同等分支才改；没有就标成未移植（注释写清 OCCT 文件与行号），不加谓词、启发式、面积比、长度滤边、体积门这类 OCCT 里不存在的规则。以前某次 Rust 跑出来的数字不是规格，OCCT 才是。
3. **对齐验证 = 跑当前实例的真实几何，对照现有基线。** 本轮修哪条几何就跑哪条（如 `fuse_box_cylinder_is_closed_solid`、`step_obj_parity` 对 `data/occ-*.obj`、差分 harness）。失配即停，报告差异。不要为对齐新写单元测试，不要用 `StubBopBuilder` 一类假输入证明算法已对齐。全量 `cargo test` 只在提交前跑一次（改单一 crate 时只跑该 crate）。整体目视/网格复查：`rtk cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj`（`data/*.step` → 仓库根 `output/<stem>.obj`，偏转 0.1）。
4. **完成即停。** 当前任务门禁满足 → 报告结果与发现的旁支问题 → 停，等下一步指令。旁支问题写进报告，不要动手修，不要自创下一步。
5. **任务卡的门禁只能引用现有基线**（既有测试名、harness/parity 输出），不得要求新写测试。源码行数可用作批次范围（翻译多少 OCCT），不得用作产出目标或成果汇报（产出多少 Rust、多少测试）。

```text
OK:  ElCLib::LineD1 在 U 处求值（cxx 有同等调用）
OK:  TotCross2D 用定向边 pcurve（BRep_Tool::CurveOnSurface 有同等分支）
NO:  按辐条长度当 unused；包围盒面积比 0.8 冒充 same-domain
NO:  为凑某一 STEP 的顶点数在 addWire 里按面积翻链；用特例 walk 代替 Edge.Reverse()
```
