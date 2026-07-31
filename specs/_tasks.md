# _tasks — OCCT → Rust migration map (current)

> Last updated: 2026-07-31. 5 crates · 179 modules · 18,959 lines · 215 tests
> Phase 3 detail: [specs/_phase3_plan.md](_phase3_plan.md)

## 已完整 (✅)

### occt-core (131 tests)
| 域 | 模块 |
|---|---|
| gp 几何 (37 类型) | gp/ 全部 — 点/线/圆/椭圆/双曲/抛物/平面/柱/锥/球/环/2D/变换 |
| 求值 | elib/ — ElCLib + ElSLib + intersect + surface_eval + measure |
| B样条/Bezier | bspl/ — knots/eval/bezier/plib/surface/poles/rational/曲线工具 |
| 包围体 | bnd/ — Box/Box2d/Sphere/OBB(PCA)/SortBox |
| 空间结构 | bvh/ — BVH + 遍历 + 三角BVH |
| 网格 | poly/ — Triangulation/Polygon2D/Polygon3D |
| 交换 | io/ — OBJ + STL + PLY 读写 |
| 建模基础 | csg/ (体素布尔) · hull/ · geom/ (采样/三角剖分/拟合/折线) |
| 质量属性 | gprop/ — 质心/惯性/面积/体积 |
| 内核 | kernel/ — OCCError/Handle/容器/TCollection/TColStd/GeomAbs/Units/OSD/Resource/字符串 |
| 其他 | numeric/ · quantity/ · message/ · toploc/ · convert/ · cslib/ · gcpnts/ · precision |

### occt-math (36 tests)
vector/matrix/intvec/status + SVD/Crout/Jacobi/Householder/Newton/BFGS/Powell/
Gauss-Legendre/Kronrod/Eigen/Trig/LeastSquares/MultiInt + stats/fft/rng

### occt-geom (TKG3d) + occt-geom2d (TKG2d)
Curve/Surface traits + Line/Circle/Ellipse/Hyperbola/Parabola/Plane/Cylinder/Cone/
Sphere/Torus + B样条/Bezier/Trimmed/Offset/Revolved + convert/transform

### occt-topo (TKBRep, 46 tests)
ShapeType/Orientation/TShape/TopoShape/typed wrappers/TopoBuilder/Explorer/
BRepTool(占位)/mesh(基本体)/primitives/transform/validate/model

## Phase 3 — BRep 几何内核 (进行中, 最高优先)

| # | 任务 | 状态 | 依赖 |
|---|---|---|---|
| A1 | TShape 挂几何 (vertex point/edge curve/face surface) | ⚪ | geom |
| A2 | BRepBuilder 真实顶点/边/面构建 | ⚪ | A1 |
| A3 | BRepTool 完整几何访问 | ⚪ | A1 |
| A4 | 真实 children 树存储 | ⚪ | A1 |
| A5 | 边→折线 / 面→网格逼近 | ⚪ | A4+B |
| A6 | 任意形状网格化 | ⚪ | A5 |
| A7 | 从几何算真实包围盒 | ⚪ | A5 |
| B1-B5 | B样条→Bezier→折线 / 曲面→UV网格 / 转换 | ⚪ | geom+bspl |
| C1-C5 | 面求交/边分裂/壳检查/网格布尔 | ⚪ | A6 |
| D1-D3 | 网格→BRep / BRep→交换 / 场景 | ⚪ | A6 |

## 远期 (⬜)
- STEP/IGES 数据交换 (需 A1 完成几何存取)
- ModelingAlgorithms 精确布尔/倒角 (需 A6 + 求交)
- 可视化渲染器 (外部工具已可看 OBJ/STL)
