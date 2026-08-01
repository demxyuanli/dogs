# Phase 7 Plan — NURBS 曲面布尔深度 + BRepFeat + 变半径倒圆 + 形状修复 + XML 装配

> 依据 `specs/_coverage.md` 缺口优先级排序。目标：本轮新增 Rust ≥ 10,000 行（最低 5,000，不足则扩阶段）。

## 缺口优先级（影响 × 可行性）

| # | 缺口 | 内容 | 预估 | 依赖 |
|---|---|---|---|---|
| **1** | NURBS 曲面求交/布尔深度 | intpatch 扩展：B样条曲面自适应 trace → NURBS 曲线；bop_curved 一般曲面面 | 2,500 | spline_surface/bspline_surface |
| **2** | TKFeat 特征深度 (BRepFeat) | 拔模/凹槽/颈部/加强筋（revolve 布尔） | 1,300 | sweep_revolve + bop_builder |
| **3** | TKFillet 变半径/链 | 沿边半径线性/样条变化 + 多边链式 | 1,100 | fillet_edge |
| **4** | TKShHealing 形状修复 | 自由边焊接/小边移除/未闭合线框修复 | 1,100 | topo_tools_full |
| **5** | XmlXCAF | XML 装配容器（补 bincaf 的文本格式） | 900 | bincaf 模式 |
| **6** | math/geom 深度 | RK45 自适应 ODE + Gauss-Newton 非线性 LS + B样条曲面求值 | 1,400 | matrix/qr |

## 里程碑

- **M1 (1)**: 一般曲面（B样条）求交曲线 + 曲面布尔 —— 突破解析曲面限制
- **M2 (2+3)**: BRepFeat 特征族 + 变半径倒圆
- **M3 (4+5)**: 形状修复 + 文本交换（XML）
- **M4 (6)**: 数值深度

## 并行

- 六个 Agent 完全独立，各自只写自己的模块文件（lib.rs 已预注册）。
- P7-A 内部：intpatch 扩展 → bop_curved 扩展（同一 Agent）。
- math/geom 与 topo 无关。
