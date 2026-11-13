# dogs — OCCT → Rust 忠实移植

本仓库把 **Open CASCADE Technology 8.0.1**（OCCT）按控制流逐步搬迁到 Rust：以官方 C++ 行为为准绳，用真实几何对拍既有基线，而不是另写一套「看起来差不多」的几何库。本地对照源码：`D:\source\OCCT-src`（tag `V8_0_1`）。

| | |
|---|---|
| 标准答案 | OCCT 8.0.1 源码行为；缺分支就标注文件:行号 |
| 验收方式 | 对拍既有基线（网格 v/f、体积、拓扑）；不为对齐新写测试 |
| 读者入口 | [`specs/occt_port_interface_review_2026-11-12.pdf`](specs/occt_port_interface_review_2026-11-12.pdf)（图文综合说明） |

**五句话读懂现状（2026-11-12）**

1. 造型 / 布尔 / 偏置 / 特征 / 网格导出等主路径已可用，并用标准模型对拍。
2. 表上「类覆盖未满 100%」多半是宏、别名、调试代码，不是功能没写完。
3. 本轮已收口产品侧肋、拔模、布尔收尾与曲线反向等缺口；特征检查 **17/17**。
4. 不值得为刷覆盖率堆空壳；下一优先是有业务约定时的产品封装与显示边序残差。
5. 工作纪律见 `CLAUDE.md` 与下方五道门禁。

---

## 体量与分层（L0→L8）

2026-11-12 实测：Rust `src` 非空行 **851,773**（3,369 个 `.rs`）；对照 OCCT `.cxx+.hxx` 约 224 万行；lib 单测 1,960；`UNPORTED` 台账 562 处（带行号）。

```
L7/L8  visual / occt 门面          可视化与总导出；occt-core 为兼容门面
L6     hlr / exchange              隐藏线；STEP/IGES/网格交换
L5     occt-topo  ★主战场          布尔 / 特征 / 网格 / 偏置 / 求交（约 51 万行）
L4     topology                    拓扑形状与 BRep 基础
L3     geom / geom2d               曲线曲面与 2D 几何
L0–L2  foundation / gp / geombase / math
```

仓库**没有根 workspace**，构建必须带 `--manifest-path`。

| 层 | crate | 文件 | 非空行 | UNPORTED |
|---|---|---:|---:|---:|
| L0 | occt-foundation | 16 | 1,178 | 1 |
| L1 | occt-gp | 163 | 33,998 | 13 |
| L2 | occt-geombase / occt-math | 74 / 84 | 14,387 / 20,672 | 13 / 9 |
| L3 | occt-geom2d / occt-geom | 221 / 492 | 54,646 / 122,104 | 58 / 88 |
| L4 | occt-topology | 53 | 10,579 | 25 |
| L5 | occt-topo | 1,956 | 512,040 | 281 |
| L6 | occt-hlr / occt-exchange | 115 / 175 | 30,641 / 44,588 | 26 / 45 |
| L7/L8 | occt-visual / occt | 13 / 6 | 4,338 / 2,579 | 2 / 1 |
| — | occt-core（门面） | 1 | 23 | 0 |

---

## 功能域就绪度

| 功能域 | 状态 | 含义 |
|---|---|---|
| 布尔运算 | 已收口 | 曲面布尔走正式建造器；旧体素引擎已删除 |
| 网格化 | 已接线 | 增量网格 + Delaunay；15 模型顶/面数与 OCCT 一致 |
| STEP / IGES | 已接线 | 可读可写；工业模型有已登记小残差 |
| 偏置 / 抽壳 | 已收口 | 主链落地；进度权重表不移植 |
| 圆角 | 已接线 | 角点族齐；少量转发外壳未单独建名 |
| 特征（孔槽肋拔模） | 已接线 | 产品入口已接；部分高级类尚无一键 API |
| 隐藏线 | 已接线 | 精确链通；立方体多边形对拍 0 差 |
| 可视化 / Draw lite | Lite | 占位，非生产级完整工具包 |

### 波次覆盖（类名命中率，勿当完成度）

| 波次 | 覆盖 | 白话 |
|---|---:|---|
| S1 填充/法则 | 100% | 72/72 收口 |
| S2 布尔底座 | 90.2% | 算法满；表差来自宏/别名/调试 |
| S2i 求交 | 95% | typedef / 枚举命名差异 |
| S3 圆角 | 94.6% | 需 S1+S2；余薄转发外壳 |
| S4 特征 | 98.6% | 72/73；余 DRAW 诊断 |
| S5 2D 相切 | 100% | 类全；多数缺生产调用 |
| S6 偏置/拔模 | 100% | 38/38；进度表不移植 |
| S7 隐藏线 | 80% | 算法通；显示命名/边序残差 |

计划内合计约 **645/692** 类（约 93%）。**不要为刷 100% 堆空壳。**

---

## 产品接口

| 产品入口 | 底层 | 可靠依据 |
|---|---|---|
| 开槽 / 颈部 | MakeRevol | 体积增减方向对 |
| 通孔 / 盲孔 | MakeCylindricalHole | 与 OCCT 探针五入口对拍 |
| 凸台 | MakePrism | 体积 / 贯穿 |
| 加强肋 | MakeLinearForm | 草图约定 + rib 检查 |
| 拔模 | DraftAngle | Direction = **中性面法向** |
| 通用布尔 | algo 门面 | 与 Form 同路径 |

使用契约：肋走 OCCT 草图线约定；部分 `tol` 仅为签名兼容；失败返回错误。类已落、产品未包：`MakeDPrism` / `MakePipe` / `MakeRevolutionForm`（有真实几何约定后再加一键入口）。

---

## 对拍与门禁

| 门禁 / 探针 | 验证什么 | 结果 |
|---|---|---|
| 编译检查 | 库可构建 | 相关 crate 0 error |
| `export_data_obj` | 15 模型网格 v/f vs OCCT OBJ | **15/15 逐位一致** |
| `step_obj_gates` | 包围盒 / 面积 | 4/4；工业模型残差 ≤0.15% 已登记 |
| 特征套件 | 孔槽肋拔模等 | **17/17** |
| `fuse_box_cylinder` | 融合后闭合实体 | 通过 |
| 曲面布尔探针 | Fuse/Cut/Common 体积 | 对齐 OCCT 真值 |
| 质量属性 | 原始体体积/面积 | 相对误差约 1e-9 |
| 多边形隐藏线 | 可见边 vs OCCT | Cube 0 差；曲面模型有边序差 |

抽样（本仓与 OCCT v/f 相同）：Cube 24/12 · Cylinder 146/140 · Sphere 642/1244 · Torus 1369/2592 · screw 600/790 · Shape 6150/11372 · …共 15 个。

### 五道移植门禁

1. **改完先编译。** `cargo check`（或该 crate 编译）不过就先修编译。
2. **对着 `.cxx` 审控制流再改。** 同等分支才改；没有就标 `UNPORTED`（OCCT 文件与行号）。不加启发式凑某一模型数字。
3. **对齐 = 跑当前实例真实几何，对照既有基线。** 失配即停。不为对齐新写测试；全量 `cargo test` 只在提交前跑一次。
4. **完成即停。** 门禁满足 → 报告结果与旁支 → 停。
5. **任务卡门禁只能引用现有基线**，不得要求新写测试。

---

## 开放缺口与下一轮

| 缺口 | 建议 |
|---|---|
| 部分高级特征无一键 API | 有业务几何约定再接线 |
| 隐藏线模板类名 / 多边形边序 | 查三角化输入与边序，不重写求交 |
| 调试转储 / 进度表 / DRAW 诊断 | **不移植**，保持登记 |
| UNPORTED 562 处 | 按行号消化；不是回归信号 |
| 2D 相切多数缺调用方 | 有生产链再做数值对齐 |

**优先做**有约定的产品封装；**按需做**显示边序残差；**不要做**为空壳刷 100%。

---

## 环境 · 构建 · 常用命令

- **Rust** stable（无 `rust-toolchain.toml`）
- 文档里的 `rtk cargo …` 可换成普通 `cargo`；日常建议 `--offline`
- 仅跑 `specs/occt_probe/` 时需要本机 OCCT 构建

```bash
# 编译（第一条门禁）
cargo check --manifest-path crates/occt-topo/Cargo.toml --offline

# 单元测试（提交前；改单一 crate 时只跑该 crate）
cargo test  --manifest-path crates/occt-topo/Cargo.toml --offline --lib

# 逐模型网格门禁
cargo test  --manifest-path crates/occt-topo/Cargo.toml --offline --test step_obj_gates -- --nocapture

# 导出 OBJ（data/*.step → data/output/<stem>.obj，偏转 0.1）
cargo run --manifest-path crates/occt-topo/Cargo.toml --offline --example export_data_obj
```

端口探针：`export_data_obj`、`model_probe`、`mesh_uv_feed_probe`、`seam_fix_probe`、`iges_check` 等（`crates/occt-topo/examples/`）。  
OCCT 侧探针：`specs/occt_probe/`（见该目录 README）。对拍跨侧必须按几何配对，不得按面号。

---

## 数据与产物

- 输入：`data/occ/*.stp` / `*.step`，以及 `data/*.step`
- 参考网格（GT）：`data/occ/occ-*.obj`
- 导出：`data/output/<stem>.obj`（gitignore）；门禁写入 `data/output/gates/`（与导出目录分离）
- 临时产物：`.target-gate/`、`tmp_*`、`out*.txt` 等一律 gitignore，不入库

---

## 文档索引

| 文档 | 用途 |
|---|---|
| [`specs/occt_port_interface_review_2026-11-12.pdf`](specs/occt_port_interface_review_2026-11-12.pdf) | **对外图文说明**（忠实度审计 · 状态 · 路线图）；同源 HTML 同名 |
| `specs/occt_migration_guide.md` | 移植指导：忠实度口径、lite 缺口、波次路线图、禁令 |
| `specs/_tasks.md` | 任务台账（含 T-537 等） |
| `specs/_a3n00_gap_analysis.md` | a3n00 缺口与实测数字 |
| `CLAUDE.md` | 工作纪律（改代码前必读） |
| Cursor 画板 | `occt-migration-audit-and-100k-roadmap` · `occt-port-status-overview`（IDE 内，不在仓内） |

本轮故事摘要：肋改走线性成型；拔模改走官方 `DraftAngle`（方向=中性面法向）；`SimplifyResult` 接线；裁剪曲线真正 `Reverse`。故意不搬调试转储、进度表、DRAW 诊断。
