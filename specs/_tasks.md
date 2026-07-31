# _tasks — OCCT → Rust migration map (current)

> Last updated: 2026-07-31. 5 crates · 29,022 lines · 386 tests
> Phase 3 detail: [specs/_phase3_plan.md](_phase3_plan.md)

## 已完整 (✅)

### occt-core (131 tests)
| 域 | 模块 |
|---|---|
| gp 几何 (37 类型) | gp/ 全部 — 点/线/圆/椭圆/双曲/抛物/平面/柱/锥/球/环/2D/变换 |
| 求值 | elib/ — ElCLib + ElSLib + intersect + surface_eval + measure |
| B样条/Bezier | bspl/ — knots/eval/bezier/plib/surface/poles/rational/曲线工具 (rational eval 修复) |
| 包围体 | bnd/ — Box/Box2d/Sphere/OBB(PCA)/SortBox |
| 空间结构 | bvh/ — BVH + 遍历 + 三角BVH |
| 网格 | poly/ — Triangulation/Polygon2D/Polygon3D |
| 交换 | io/ — OBJ + STL + PLY 读写 |
| 建模基础 | csg/ (体素布尔) · hull/ · geom/ (采样/三角剖分/拟合/折线) |
| 质量属性 | gprop/ — 质心/惯性/面积/体积 |
| 内核 | kernel/ — OCCError/Handle/容器/TCollection/TColStd/GeomAbs/Units/OSD/Resource/字符串 |
| 其他 | numeric/ · quantity/ · message/ · toploc/ · convert/ · cslib/ · gcpnts/ · precision |

### occt-math (70 tests)
vector/matrix/intvec/status + SVD/Crout/Jacobi/Householder/Newton/BFGS/Powell/
Gauss-Legendre/Kronrod/Eigen/Trig/LeastSquares/MultiInt + Cholesky/LU/QR-full/
cubic-spline/polyfit/Levenberg/RK4/bisection/secant/golden/brent/polynomial
roots + matrix_ext + stats/fft/rng

### occt-geom (TKG3d, 19 tests)
Curve/Surface traits + Line/Circle/Ellipse/Hyperbola/Parabola/Plane/Cylinder/Cone/
Sphere/Torus + B样条/Bezier/Trimmed/Offset/Revolved + **B1-B5 转换**:
bspline_to_bezier(正确全局算法) · bezier_to_polyline · curve_approx ·
surface_to_grid · convert_geom

### occt-geom2d (TKG2d, 21 tests)
Curve2d trait + Line/Circle/Ellipse/Hyperbola/Parabola/Trimmed/Offset/BSpline +
Geom2dBezierCurve · curve_ops(长度/交点/最近点) · bspline2d_to_bezier

### occt-topo (TKBRep, 145 tests)
- **几何内核**: tgeometry(GeometryRegistry 侧表) · 真实 children 树 ·
  builder(全部 make_*) · brep_tool(BRep_Tool) · primitives(真几何 box/cyl/sphere/cone/torus)
- **网格化**: wireframe(边→折线/面→三角) · shape_mesh(任意 BRep 网格化+法线+焊接) ·
  bbox_from_geometry(几何真实包围盒)
- **拓扑分析**: face_face(面求交) · edge_split · shell_check(闭合壳) ·
  solid_union(体素并集) · boolean_ops(体素交/差) · brep_extrema(距离/inside) ·
  shape_analysis(有效性) · brep_gprop(质量属性)
- **交换**: mesh_to_brep(三角→BRep) · brep_exchange(BRep→OBJ/STL/PLY) · brep_scene(多形状场景)
- **其他**: sweep(棱柱拉伸) · brep_surface(曲面分类/pcurve) · topo_tools_full(TopExp) ·
  brep_measure(解析测量) · mesh · transform · validate · model · topexp

## 里程碑状态

| 里程碑 | 状态 |
|---|---|
| M1 (A1-A4): BRep 可存几何、可查询 | ✅ make_box 真几何, BRepTool 取回坐标/曲面 |
| M2 (A5-A7 + B): 任意形状可网格化/可视化/导出 | ✅ 任意 BRep 三角化 → OBJ/STL/PLY |
| M3 (C+D): 网格导入 + 基础布尔 | ✅ 体素布尔 + mesh→BRep + 交换 |

## 远期 (⬜)
- STEP/IGES 数据交换 (需 location 应用到几何存取)
- 精确 NURBS 布尔/倒角 (当前为体素近似)
- 可视化渲染器 (外部工具已可看 OBJ/STL)
