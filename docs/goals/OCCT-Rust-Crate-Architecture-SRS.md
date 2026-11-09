# OCCT Rust 移植 — Crate 架构需求规格说明

| 项目 | 内容 |
|------|------|
| 文档标题 | Open CASCADE Technology (OCCT) Rust 移植 Crate 拆分与架构需求规格说明 |
| 版本 | 0.1.0 |
| 状态 | 草案 |
| 创建日期 | 2026-11-09 |
| 适用范围 | 将 OCCT 核心建模能力移植 / 重写为纯 Rust 实现时的 crate 划分 |

---

## 1. 概述

### 1.1 目的

本文档定义将 Open CASCADE Technology（OCCT）移植到 Rust 时，**crate 的拆分原则、层级结构、职责边界与依赖关系**。

目标是在保证正确性的前提下，最大化：

- 编译并行度
- 增量编译效率
- 模块可维护性
- 后期算法独立演进能力

### 1.2 背景

OCCT 是一个超大型 C++ CAD 内核（约 200 万+ 行代码），采用严格的分层架构。直接将其作为一个巨型 crate 移植会导致：

- 编译时间极长
- 修改任意底层代码引发大规模重编译
- 团队协作冲突严重
- 二进制体积与链接时间失控

因此必须按照 OCCT 原有分层思想，结合 Rust/Cargo 的编译模型进行 crate 级拆分。

### 1.3 范围

**包含：**

- Foundation、Math、Geometry、Topology、BRep
- 核心建模算法（布尔、圆角、偏移、基本体、修复、网格）
- 数据交换（STEP/IGES/STL 等）的基础框架
- 门面（Facade）crate 的组织方式

**不包含（本期）：**

- 完整可视化模块（TKV3d / OpenGl 等）的详细拆分
- OCAF（Application Framework）完整移植
- 与具体 GUI 框架的集成
- 性能基准与数值精度验证规范（另文规定）

---

## 2. 总体架构要求

### 2.1 分层原则

必须严格遵循自底向上的依赖方向，禁止循环依赖：

```
Foundation
    ↑
  Math (gp)
    ↑
 Geometry (Geom / Geom2d)
    ↑
 Topology (TopoDS)
    ↑
   BRep
    ↑
 Algorithm Core
    ↑
┌────┼────┬────────┬─────────┐
│    │    │        │         │
Bool Fillet Offset Prim   Healing / Mesh
    ↑
 Exchange（可选依赖上层算法）
    ↑
  Facade（occt）
```

### 2.2 设计目标

| 编号 | 目标 | 说明 |
|------|------|------|
| G1 | 编译隔离 | 修改高频算法（布尔、圆角）不应导致几何/拓扑层重编译 |
| G2 | 并行编译 | 多个算法 crate 应尽可能可并行编译 |
| G3 | 稳定地基 | Foundation 与 Math 变更频率极低，作为稳定依赖 |
| G4 | 职责单一 | 每个 crate 只负责一个清晰的领域边界 |
| G5 | 可裁剪 | 通过 feature 可只编译需要的算法子集 |
| G6 | 可测试 | 每个 crate 应能独立进行单元测试与集成测试 |

---

## 3. Crate 列表与职责

### 3.1 核心基础层

#### 3.1.1 `occt-foundation`

- **对应 OCCT**：TKernel
- **职责**：
  - 基本类型、句柄语义（或 Rust 所有权替代方案）
  - 集合与容器抽象
  - 错误/异常模型映射
  - 内存与工具函数
- **依赖**：无
- **稳定性要求**：极高，接口变更需严格评审

#### 3.1.2 `occt-math`

- **对应 OCCT**：TKMath（gp 包）
- **职责**：
  - 点、向量、矩阵、变换（Trsf）
  - 坐标系（Ax1、Ax2、Ax3）
  - 基本几何计算
- **依赖**：`occt-foundation`
- **稳定性要求**：极高

### 3.2 建模数据层

#### 3.2.1 `occt-geometry`

- **对应 OCCT**：TKG2d、TKG3d、TKGeomBase
- **职责**：
  - 2D/3D 曲线与曲面定义（解析曲线、BSpline、Bezier 等）
  - GeomAdaptor 相关适配
  - 几何求值与基本属性
- **依赖**：`occt-math`
- **关键约束**：不得依赖任何拓扑概念

#### 3.2.2 `occt-topology`

- **对应 OCCT**：TopoDS
- **职责**：
  - 拓扑实体：Vertex、Edge、Wire、Face、Shell、Solid、Compound 等
  - 方向、位置、形状枚举
  - 拓扑遍历基础接口
- **依赖**：`occt-geometry`
- **关键约束**：只描述“连接关系”，不拥有具体几何实现细节以外的算法

#### 3.2.3 `occt-brep`

- **对应 OCCT**：TKBRep
- **职责**：
  - BRep_Builder、BRep_Tool
  - 几何与拓扑的绑定
  - 位置与方向的实际应用
- **依赖**：`occt-topology`

### 3.3 算法层

#### 3.3.1 `occt-algo-core`

- **职责**：
  - 被多个算法共用的基础能力（相交、投影、极值、分类等）
  - 算法上下文、进度、中断接口
- **依赖**：`occt-brep`

#### 3.3.2 `occt-bool`

- **对应 OCCT**：TKBO、TKBool
- **职责**：布尔运算（Fuse、Cut、Common）及历史记录
- **依赖**：`occt-algo-core`
- **拆分理由**：复杂度高、编译重、变更频繁，必须独立

#### 3.3.3 `occt-fillet`

- **对应 OCCT**：TKFillet
- **职责**：圆角、倒角
- **依赖**：`occt-algo-core`
- **拆分理由**：与布尔同等重量级，独立演进

#### 3.3.4 `occt-offset`

- **对应 OCCT**：TKOffset
- **职责**：偏移、加厚、抽壳
- **依赖**：`occt-algo-core`

#### 3.3.5 `occt-prim`

- **对应 OCCT**：TKPrim 及部分 TopAlgo
- **职责**：
  - 基本体（Box、Cylinder、Sphere、Cone、Torus…）
  - 拉伸、旋转、扫掠、放样、管道
- **依赖**：`occt-algo-core`

#### 3.3.6 `occt-healing`

- **对应 OCCT**：TKShHealing
- **职责**：形状修复、统一域、缝合等
- **依赖**：`occt-brep`（可按需依赖 algo-core）

#### 3.3.7 `occt-mesh`

- **对应 OCCT**：TKMesh
- **职责**：BRep 三角化、网格数据结构
- **依赖**：`occt-brep`

### 3.4 数据交换与门面

#### 3.4.1 `occt-exchange`

- **对应 OCCT**：Data Exchange 相关 Toolkit
- **职责**：STEP、IGES、STL、BREP 等格式的读写框架
- **依赖**：按需依赖建模相关 crate（通过 feature 控制）

#### 3.4.2 `occt`（门面 crate）

- **职责**：
  - 统一 re-export 常用类型与函数
  - 提供更符合 Rust 习惯的高级 API（可选）
  - 通过 feature 聚合功能子集
- **依赖**：上述所有 crate（均为 optional + feature）

---

## 4. 依赖关系约束

### 4.1 强制规则

1. **禁止循环依赖**。任何两个 crate 不得直接或间接互相依赖。
2. **依赖方向只能向上**（从底层指向高层）。
3. 算法 crate（bool、fillet、offset、prim 等）之间**默认不得互相依赖**。如有共用逻辑，必须下沉到 `occt-algo-core`。
4. `occt-geometry` 不得依赖 `occt-topology` 或更高层。
5. 可视化相关功能（若后续增加）不得反向污染建模核心。

### 4.2 推荐依赖图（简化）

```
occt-foundation
└── occt-math
    └── occt-geometry
        └── occt-topology
            └── occt-brep
                └── occt-algo-core
                    ├── occt-bool
                    ├── occt-fillet
                    ├── occt-offset
                    ├── occt-prim
                    ├── occt-healing
                    └── occt-mesh
                        └── occt-exchange
                            └── occt（facade）
```

### 4.3 修订（基于实测依赖图的修正）

本节记录对 §4.1 与 §3.3 的三处修正，依据是对现有代码的模块级依赖统计
（`crates/occt-topo` 的模块级 SCC 计算结果与 `crates/*/Cargo.toml` 实测）。修订前的规则在该仓库的现有代码上不可实现。

1. **§4.1 第 3 条修正为：算法 crate 之间允许有向无环依赖，禁止环。**
   实测 `brep_offset`（TKOffset）依赖 `chfi3d`（TKFillet）、`bop_builder2`/`topopebrep`（TKBool）、`brep_fill_sweep`、`meshing`；
   `chfi3d`（TKFillet）依赖 `topopebrep`、`blend`/`brep_blend`、`intpatch`。这与 OCCT 自身 toolkit 的依赖方向一致，
   因此可实现的顺序是 `brep -> algo-core -> 经典布尔(TopOpeBRep) -> prim -> fillet -> offset`。
2. **§3.3.7 的 `occt-mesh` 位置修正为两层：** 三角化**数据结构**下沉到 brep 层（`Poly_Triangulation` 已在 `occt-core::poly`），
   只有 `BRepMesh_IncrementalMesh` 这类**算法**留在上层。依据是 `fclass2d`、`brep_extrema`、`brep_lib_build_curves3d`、
   `brep_lib_continuity`、`brep_lib_same_parameter`、`brep_offset`、`brep_tools`、`shhealing` 均直接依赖 `meshing`。
3. **§3.4.1 的 `occt-exchange` 只包含文件格式读写，不包含 BRep 交换对象。**
   `step`、`iges`、`iges_read`、`xmlcaf`、`bincaf`、`rwmesh`、`xcaf*`、`model`、`brep_exchange` 位于所有建模算法之上，
   建模层不反向依赖它们；但 `BRepTools_ShapeSet` 一类 BRep 交换对象属于 TKBRep，不得随格式读写一起上移。
   依据是 `step` 依赖 `shhealing`、`iges_read` 依赖 `brep_fill`、`brep_offset_api_make_pipe` 依赖 `brep_fill_pipe`。

另外记录两条实测约束，用于排期：

- `occt-topo` 当前存在一个 **204 个模块的强连通分量**，同时跨越 TKBRep、TKTopAlgo 与全部算法层。
  Rust 的 crate 依赖必须无环，所以在打断 `tshape -> brep_tool`、`brep_extrema -> shape_mesh`、
  `fclass2d -> meshing`、`brep_tools -> meshing` 这类回边之前，同一个 SCC 内的分层无法拆开。
  实测的候选边界与回边代价见 `docs/goals/occt-crate-split/T007-scc-analysis.md`。
- 拆分必须**自底向上**：被拆出的下层 crate 由原 crate 用 `pub use` 重新导出，既有 `occt_core::gp::...` 一类路径才能保持不变。
  反向操作（把上层先搬走）会让原 crate 无法再导出依赖自己的 crate，导致 13 个门禁测试与 27 个 example 的路径全部失效。

### 4.4 当前落地状态（2026-11-10）

| crate | 来源 | 行数 | 状态 |
|-------|------|------|------|
| `occt-foundation` | `occt-core` 的 kernel / precision / message / quantity | 1,178 | 已拆出 |
| `occt-gp` | `occt-core` 的 gp / bnd / elib / cslib / bspl / poly / gprop / toploc / convert / approx / int* / math_* | 34,021 | 已拆出 |
| `occt-geombase` | `occt-core` 的 geom / intcurve / intf / intimpargen / adv_approx / bvh / io | 14,307 | 已拆出 |
| `occt-core` | 上述三个 crate 的 `pub use` 门面 | 23 | 已退化为门面 |
| `occt-math` | 原有 TKMath 求解器端口 | 20,646 | 未动（与 `occt-gp` 职责重叠，待合并） |
| `occt-geom2d` / `occt-geom` | 原有 | 54,672 / 117,428 | 未动 |
| `occt-topology` | `occt-topo` 的 TKBRep 数据层（abs / shape / tshape / iterator / transform / topexp / tools / validate / debug_flags / builder / brep_tool / brep_point_representation / brep_builder_update_vertex / brep_lprop / brep_tools_substitution / tgeometry / topo_tools / adaptor / brep_tools{,_uv,_modifier,_wire_explorer} / brep_uv_bounds / brep_bnd_lib / brep_surface） | 9,761 | 已拆出（批 B1b + B2） |
| `occt-topo` | 剩余的 TKTopAlgo / TKBO / TKBool / TKFillet / TKOffset / TKPrim / TKShHealing / TKMesh | 508,031 | 已缩减，仍含 204 模块 SCC |
| `occt-hlr` | `occt-topo` 的 HLR / Contap | 30,654 | 已拆出 |
| `occt-exchange` | `occt-topo` 的 STEP / IGES / XCAF / RWMesh / BRepExchange / 装配 | 44,301 | 已拆出 |
| `occt-visual` | `occt-topo` 的 viz_scene / render_svg / brep_scene | 4,338 | 已拆出 |
| `occt` | `draw` 门面 | 2,579 | 已拆出 |

当前正常依赖链（`cargo tree -e normal` 实测，无环）：

```text
occt-base  occt-core / occt-math / occt-geom2d / occt-geom
occt-topology  -> occt-core, occt-geom, occt-geom2d, occt-math
occt-topo       -> occt-topology
occt-hlr        -> occt-topo
occt-exchange   -> occt-topo
occt-visual     -> occt-hlr, occt-exchange, occt-topo
occt            -> occt-topo, occt-hlr, occt-exchange, occt-visual
```

`occt-topo` 用 `pub use occt_topology::{abs, shape, tshape, builder, brep_tool, ...}` 保持既有
`occt_topo::*` 路径不变，因此 13 个门禁测试、`examples/*` 与所有下游 crate 无需改路径。

B2 批次已落地（`docs/goals/occt-crate-split/state.yaml` T010）：`cut-B2.txt` 的 17 条回边中 7 条可逆的
已按“下沉共用 helper”打断——`write_triangulation` 留在 `occt-topo`、`GeomBndLib_*` 下沉到 `occt-geom`、
`BRep_Tool::CurveOnSurface` 下沉到 `occt_topology::brep_tool_curve_on_surface`、`SurfaceAdaptor` 下沉到
`occt_topology::adaptor`、`empty_copied_*` 与面-边连续性存储下沉到 `occt_topology::{shape_empty_copied, tshape_continuity}`。
`brep_faces` / `brep_measure` **没有**下沉：`brep_measure::edge_length` 是 `BRepGProp::LinearProperties`
（`BRepGProp.cxx`）的高斯积分，属 TKTopAlgo 算法，OCCT 本身也把 `BRepGProp` 放在 TKTopAlgo，
把它切开不是移动 helper 而是改控制流。
`brep_tools -> mesh`、`brep_tools -> meshing` 属端口缺陷（OCCT 的 `BRep_Tool::Triangulation`
只碰 `Poly_Triangulation`），已按此拆分。

`occt-math` 与本文 §3.1.2 定义的 `occt-math` 语义不一致：现 `occt-math` 只承载求解器（SVD、Gauss、BFGS 等），
TKMath 的 `gp` / `bspl` 部分已归入 `occt-gp`。两者是否合并需单独评审。

---

## 5. 非功能性需求

### 5.1 编译与构建

- 修改 `occt-bool` 或 `occt-fillet` 时，不应触发 `occt-geometry` / `occt-topology` 的重编译。
- 支持 Cargo workspace 统一管理版本与依赖。
- 支持通过 feature 裁剪，实现“最小建模核心”与“完整内核”两种构建配置。

### 5.2 可维护性

- 每个 crate 必须有清晰的模块边界与文档。
- 公共 API 变更需考虑语义化版本。
- 鼓励在 crate 内部再按 OCCT 原 package 进行 module 划分（保持可追溯性）。

### 5.3 可测试性

- 每个 crate 应能独立运行 `cargo test`。
- 算法 crate 应支持对标 OCCT 的回归测试（测试数据与精度要求另文规定）。

### 5.4 二进制与分发（指导性）

- 最终应用可通过 feature 只链接需要的算法，控制最终二进制体积。
- 不强制静态/动态链接策略，由具体实现决定。

---

## 6. 实施优先级

| 阶段 | 内容 | 目标 |
|------|------|------|
| P0 | `occt-foundation` + `occt-math` | 建立稳定地基 |
| P1 | `occt-geometry` → `occt-topology` → `occt-brep` | 完成建模数据层 |
| P2 | `occt-algo-core` + `occt-prim` | 具备基本造型能力 |
| P3 | `occt-bool` + `occt-fillet` | 核心高级建模能力 |
| P4 | `occt-offset` + `occt-healing` + `occt-mesh` | 完善算法集 |
| P5 | `occt-exchange` + `occt` 门面 | 可用性与集成 |

---

## 7. 目录结构示例

```text
occt-rs/
├── Cargo.toml                          # workspace 根
├── crates/
│   ├── foundation/
│   ├── math/
│   ├── geometry/
│   ├── topology/
│   ├── brep/
│   ├── algo-core/
│   ├── bool/
│   ├── fillet/
│   ├── offset/
│   ├── prim/
│   ├── healing/
│   ├── mesh/
│   ├── exchange/
│   └── occt/                           # 门面
└── ...
```

根 `Cargo.toml` 示例关键部分：

```toml
[workspace]
resolver = "2"
members = [
    "crates/foundation",
    "crates/math",
    "crates/geometry",
    "crates/topology",
    "crates/brep",
    "crates/algo-core",
    "crates/bool",
    "crates/fillet",
    "crates/offset",
    "crates/prim",
    "crates/healing",
    "crates/mesh",
    "crates/exchange",
    "crates/occt",
]

[workspace.package]
version = "0.1.0"
edition = "2021"
license = "LGPL-2.1"          # 需根据实际许可策略调整
```

---

## 8. 风险与开放问题

| 编号 | 问题 | 说明 |
|------|------|------|
| R1 | 句柄与生命周期 | OCCT 大量使用 Handle，Rust 中需重新设计所有权与共享模型 |
| R2 | 数值精度与容差 | 移植后必须建立与 OCCT 对齐的容差与回归测试体系 |
| R3 | 算法正确性 | 布尔与圆角是高风险区域，需要大量对标测试 |
| R4 | 编译时间 vs 拆分粒度 | 过细拆分会增加 Cargo 开销，需在实践中平衡 |
| R5 | 许可兼容 | 需明确最终 crate 的许可与 OCCT 原许可的合规性 |

---

## 9. 文档修订记录

| 版本 | 日期 | 说明 |
|------|------|------|
| 0.1.0 | 2026-11-09 | 初始草案，基于 OCCT 官方分层与 Rust 编译特性整理 |
| 0.2.0 | 2026-11-10 | 依据实测模块依赖图修正 §4.1 第 3 条、§3.3.7 与 §3.4.1，新增 §4.3 修订与 §4.4 落地状态 |

---

**文档结束**
