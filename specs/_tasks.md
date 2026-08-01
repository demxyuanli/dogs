# _tasks — OCCT → Rust migration map (current)

> Last updated: 2026-08-01. 5 crates · 74,250 lines · 1,069 tests
> Phase 3 detail: [specs/_phase3_plan.md](_phase3_plan.md)

## 已完整 (✅)

### occt-core (208 tests)
| 域 | 模块 |
|---|---|
| gp 几何 (37 类型) | gp/ 全部 — 点/线/圆/椭圆/双曲/抛物/平面/柱/锥/球/环/2D/变换 |
| 求值 | elib/ — ElCLib + ElSLib + intersect + surface_eval + measure |
| B样条/Bezier | bspl/ — knots/eval/bezier/plib/surface/poles/rational/曲线工具 (rational eval 修复) |
| 包围体 | bnd/ — Box/Box2d/Sphere/OBB(PCA)/SortBox |
| 空间结构 | bvh/ — BVH + 遍历 + 三角BVH |
| 网格 | poly/ — Triangulation/Polygon2D/Polygon3D |
| 交换 | io/ — OBJ + STL + PLY 读写 |
| 建模基础 | csg/ (体素布尔) · hull/ · geom/ (采样/三角剖分/拟合/折线/多边形布尔) |
| 质量属性 | gprop/ — 质心/惯性/面积/体积 |
| 内核 | kernel/ — OCCError/Handle/容器/TCollection/TColStd/GeomAbs/Units/OSD/Resource/字符串 |
| 其他 | numeric/ · quantity/ · message/ · toploc/ · convert/ · cslib/ · gcpnts/ · precision |

### occt-math (160 tests)
vector/matrix/intvec/status + SVD/Crout/Jacobi/Householder/Newton/BFGS/Powell/
Gauss-Legendre/Kronrod/Eigen/Trig/LeastSquares/MultiInt + Cholesky/LU/QR-full/
cubic-spline/polyfit/Levenberg/RK4/bisection/secant/golden/brent/polynomial
roots + matrix_ext + stats/fft/rng + **P6**: spline_surface (张量积 B样条曲面插值) ·
eig_qr (Hessenberg+QR 一般矩阵特征) · nurbs_fit (最小二乘 B样条拟合) ·
distrib_extra (χ²/t/F/Γ/Beta/Poisson 分布族)

### occt-geom (TKG3d, 52 tests)
Curve/Surface traits + Line/Circle/Ellipse/Hyperbola/Parabola/Plane/Cylinder/Cone/
Sphere/Torus + B样条/Bezier/Trimmed/Offset/Revolved + **B1-B5 转换** +
**P6 extrema**: 点-曲线/曲线-曲线/曲线-曲面/曲面-曲面最近距离
(Extrema_ExtPC/CC/PS/CS/SS)

### occt-geom2d (TKG2d, 38 tests)
Curve2d trait + Line/Circle/Ellipse/Hyperbola/Parabola/Trimmed/Offset/BSpline +
Geom2dBezierCurve · curve_ops(长度/交点/最近点) · bspline2d_to_bezier

### occt-topo (TKBRep, 381 tests)
- **几何内核**: tgeometry(GeometryRegistry 侧表) · 真实 children 树 ·
  builder(全部 make_*) · brep_tool(BRep_Tool) · primitives(真几何 box/cyl/sphere/cone/torus)
- **网格化**: wireframe · shape_mesh(任意 BRep 网格化+法线+焊接) · brepmesh(Deflection 自适应)
- **拓扑分析**: face_face · edge_split · shell_check · solid_union ·
  boolean_ops · brep_extrema · shape_analysis · brep_gprop
- **精确建模 (P5+P6)**: inttools + bop_builder (平面精确布尔) ·
  **intpatch (曲面解析求交) + bop_curved (曲面实体布尔)** ·
  **fillet_edge (滚动球边/角倒圆)** · **brep_offset (多边形/面/壳偏移)**
- **交换**: STEP(ISO 10303-21) · IGES(ANSI Y14.26M) · OBJ/STL/PLY/VRML ·
  **bincaf (二进制装配) · gltf (glTF 2.0)**
- **其他**: sweep · loft · pipe · sweep_revolve · fillet · hlr · render_svg ·
  feature(凸台/凹槽) · patterns · shape_metrics · brep_projection · xcaf

## 里程碑状态

| 里程碑 | 状态 |
|---|---|
| M1 (A1-A4): BRep 可存几何、可查询 | ✅ make_box 真几何, BRepTool 取回坐标/曲面 |
| M2 (A5-A7 + B): 任意形状可网格化/可视化/导出 | ✅ 任意 BRep 三角化 → OBJ/STL/PLY/IGES/glTF |
| M3 (C+D): 网格导入 + 基础布尔 | ✅ 体素布尔 + mesh→BRep + 交换 |
| M4 (E, P5): 平面精确布尔 + Deflection 网格 + SVG | ✅ 精确平面 Fuse/Cut/Common |
| M5 (F, P6): 曲面求交/布尔 + 边倒圆 + 偏移 | ✅ 球/柱面精确布尔 + 滚动球倒圆 |

## 远期 (⬜)
- 一般 NURBS 曲面布尔的完整拓扑（当前 B样条 trace + 网格重建）
- TKFillet 链式共享顶点角块、变半径链式
- RWMesh 纹理/VRML 读、BinXCAF 完整 schema
- 整个 Visualization（TKOpenGl 硬件、TKV3d 交互选择、TKService 字体/纹理）
- Draw / TEST 交互
- BRepGProp / Extrema 解析深度
