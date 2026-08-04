# graphify
- **graphify** (`.claude/skills/graphify/SKILL.md`) - any input to knowledge graph. Trigger: `/graphify`
When the user types `/graphify`, use the installed graphify skill or instructions before doing anything else.

# 测试纪律（迭代 vs 门禁）
- **迭代阶段只跑相关模块测试**：`cd crates/occt-topo && cargo test --lib <module>::`（或按需 `--test <integration>`），不要反复跑全量 `cargo test`——全量约 60s+，迭代中重复跑是浪费。
- **全量回归只在提交前跑一次**（`cargo test` 全 crate），作为提交门禁。修改只在单一 crate 时，全量只需跑该 crate + 确认其它 crate 未改（`git status` 检查）。
- 调试/网格化等探索性改动，先用临时 `tests/_dbg.rs` 单测验证，确认后再并入正式测试。
