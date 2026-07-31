# Phase 5 Plan — 精确建模算法 + 可视化/交换深度

> 依据 `specs/_coverage.md` 缺口优先级排序。目标：本轮新增 Rust ≥ 10,000 行。

## 缺口优先级（影响 × 可行性）

| # | 缺口 | 内容 | 预估 | 依赖 |
|---|---|---|---|---|
| **1** | 精确布尔 (BOPAlgo) | inttools 精确求交 + bop_builder 布尔构建（平面多面体） | 3,000 | edge_split/polygon_ops |
| **2** | 真实 BRepMesh | Deflection 自适应离散（曲线/曲面细分） | 1,300 | gcpnts |
| **3** | IGES 交换 | ANSI Y14.26M 写器 | 1,000 | step 模式 |
| **4** | 可视化-lite | SVG 线框渲染器（投影 + 画家算法） | 700 | hlr |
| **5** | TKFeat-lite | 凸台/凹槽/拔模特征 | 800 | boolean |
| **6** | math 深度 | 分布采样/二维插值/稀疏矩阵 CG/FFT2D | 1,300 | — |
| **7** | XCAF-lite | STEP 装配元数据（名称/颜色/层级） | 600 | step |

## 里程碑

- **M1 (1)**: 平面多面体的精确 Fuse/Cut/Common —— 与体素布尔的可验证替代，首个"精确"建模算法
- **M2 (2+3)**: 任意 BRep 自适应网格化 + IGES 导出 → 双格式交换（STEP+IGES）
- **M3 (4+5)**: 内置线框渲染 + 特征建模（凸台/凹槽）
- **M4 (6+7)**: 数学/元数据补全

## 并行
- inttools 与 brepmesh/iges/feature 相互独立；bop_builder 依赖 inttools（固定 API 契约）。
- math/geom2d 与 topo 无关。
