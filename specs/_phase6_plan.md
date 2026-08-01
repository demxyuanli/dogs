# Phase 6 Plan — IntPatch 曲面布尔 + Fillet 深度 + Offset + 交换补全

> 依据 `specs/_coverage.md` 缺口优先级排序。目标：本轮新增 Rust ≥ 10,000 行（最低 5,000，不足则扩阶段）。
>
> 进度（2026-08-01）：P6-C/D/E 完成；P6-A/B 并行编写中；本人并行模块 polygon_boolean + geom extrema 完成，并修复 gp 变换方向向量 bug。

## 缺口优先级（影响 × 可行性）

| # | 缺口 | 内容 | 预估 | 依赖 |
|---|---|---|---|---|
| **1** | IntPatch 曲面求交 | 平面/球/柱/锥/环 精确相交曲线 + 一般 trace 回退 | 1,600 | occt-geom Surface |
| **2** | 曲面布尔 (BRepAlgoAPI) | bop_curved：带曲面面的 Fuse/Cut/Common（球∩盒/球∩球） | 1,400 | intpatch |
| **3** | 精确边倒圆 (TKFillet) | fillet_edge：沿边链滚动球圆柱/球面混合 | 900 | brep_connect |
| **4** | BRepOffset 深度 | 多边形偏移/面偏移/壳偏移（盒子±d 精确） | 900 | plane_plane_plane |
| **5** | BinXCAF + glTF | 二进制装配容器 + glTF 2.0 网格导出 | 1,000 | shape_mesh |
| **6** | math 深度 | B样条曲面插值/QR 特征/最小二乘 NURBS/分布族 | 1,600 | matrix/qr |

## 里程碑

- **M1 (1+2)**: 球/柱/锥相交曲线解析求解 + 曲面实体精确布尔 —— 突破"平面多面体"限制
- **M2 (3+4)**: 实体边倒圆（圆柱混合面）+ 盒子精确偏移
- **M3 (5)**: 双交换格式（BinXCAF 二进制 + glTF JSON/BIN）
- **M4 (6)**: 数值深度（QR 一般矩阵特征、曲面插值、分布族）

## 并行

- 五个 Agent 完全独立，各自只写自己的模块文件（lib.rs 已预注册）。
- intpatch → bop_curved 为同一 Agent 内部依赖。
- math 与 topo 无关。
