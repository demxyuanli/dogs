# Phase 13 Plan — TKGeomAlgo 深度：Extrema 72 类 + IntAna + GCPnts + GeomConvert（≥8,500 行）

> 覆盖 `_coverage.md` 最大缺口：Extrema <10%（采样近似）→ 解析/牛顿；GCPnts 45%→90%；IntAna ~40%→90%；GeomConvert 补齐。
> 参考源：`D:\source\occt-src\src\ModelingData\TKGeomBase\{Extrema,GCPnts,IntAna,GeomConvert}` + `ModelingAlgorithms\TKTopAlgo\IntCurvesFace`。
> 基线：5 crate · 114,016 行 · 1,485 测试（topo 912）。纯数值几何 → 复用 C++ 差分 harness（numeric_tol oracle）验证。

## 缺口 → 方案

| # | 缺口 | 方案 | 模块（新建，各 agent 独立文件） | 预估 |
|---|---|---|---|---|
| 1 | 点-曲线极值（3D+2D） | ExtPC/ExtPC2d（网格+牛顿）、ExtPElC/ExtPElC2d（直线/圆/椭/双曲/抛物解析解）、LocateExtPC/GLocateExtPC、POnCurv/POnCurv2d | `occt-geom/extrema_pc.rs` + `occt-geom2d/extrema2d.rs` 扩展 | 2,000 |
| 2 | 曲线-曲线极值（3D+2D） | ExtCC/ExtCC2d、ExtElC/ExtElC2d（解析配对）、LocateExtCC/ECC、GenLocateExtCC | `occt-geom/extrema_cc.rs` + `occt-geom2d/extrema2d.rs` 扩展 | 1,800 |
| 3 | 点/曲线-曲面极值 | ExtPS、ExtPElS（平面/柱/锥/球/torus 解析）、ExtPRevS、ExtPExtS、ExtCS/ExtElCS、FuncExtCS | `occt-geom/extrema_surf.rs` | 1,800 |
| 4 | 曲面-曲面极值 | ExtSS/ExtElSS、GenExtSS、FuncExtSS、FuncPSDist/FuncPSNorm | `occt-geom/extrema_ss.rs` | 1,200 |
| 5 | 解析求交 | IntAna 全（Int3Pln/IntLinTorus/IntConicQuad/IntQuadQuad/QuadQuadGeo）+ IntCurvesFace 2 类 | `occt-geom/intana.rs` + `occt-topo/brep_face_intersect.rs` | 1,200 |
| 6 | 离散点 / B样条转换 | GCPnts（AbscissaPoint/Uniform/QuasiUniform/TangentialDeflection）+ GeomConvert（KnotSplitting/BSpline↔Bezier/CompCurve/ApproxCurve） | `occt-geom/gcpnts.rs` + `occt-geom/convert_bspl.rs` | 1,800 |

## OCCT 类 → Rust 映射（Extrema 72 类主类）

- **枚举/工具**：Extrema_ExtAlgo/ExtFlag/ElementType、CurveTool/Curve2dTool、POnCurv/POnCurv2d、POnSurf/POnSurfParams
- **点-曲线**：ExtPC/2d、ExtPElC/2d、LocateExtPC/2d、GLocateExtPC、GenLocateExtPC
- **曲线-曲线**：ExtCC/2d、ExtElC/2d、LocateExtCC/2d、ECC/2d、GenLocateExtCC、GenExtCC
- **曲面侧**：ExtPS、ExtPElS、ExtPRevS、ExtPExtS、GenExtPS/GenLocateExtPS、ExtCS/ExtElCS、GenExtCS/GenLocateExtCS、ExtSS/ExtElSS、GenExtSS/GenLocateExtSS、FuncExtCS/FuncExtSS/FuncPSDist/FuncPSNorm
- **内部模板（不逐类移植，行为内联）**：`*OF*`/`PCF*`/`PCLocF*`/`LocE*` 等 .lxx 辅助 → 并入对应主类实现
- **全局优化函数**（GlobOptFunc*/GFunc*/GGen*/GGE*）→ 可选，第一波用网格+牛顿替代，标记 `ponytail: global opt 后补`

## 交付波次

- **波 1（点-曲线，~2,000 行）**：枚举/POnCurv + ExtPC（网格初始化 + 牛顿收敛）+ ExtPElC 解析（直线/圆/椭/双曲/抛物）+ LocateExtPC/GLocateExtPC（种子局部）。替换 `point_curve_extrema` 采样路径。
- **波 2（曲线-曲线，~1,800 行）**：ExtCC + ExtElC 解析配对 + LocateExtCC/ECC。替换 `curve_curve_extrema` 采样。
- **波 3（曲面侧，~3,000 行）**：ExtPS/ExtPElS/ExtPRevS/ExtPExtS + ExtCS/ExtElCS/FuncExtCS + ExtSS/ExtElSS/FuncExtSS。替换 `point_surface_extrema`/`curve_surface_extrema`/`surface_surface_extrema`。
- **波 4（求交+离散+转换+接入，~3,000 行）**：IntAna 全 + IntCurvesFace + GCPnts + GeomConvert；接入 `brep_extrema.rs`（point_shape_distance/closest_point_on_edge/face 改走解析），`lib.rs` 统一注册新模块。

## 差分验证（每波必做）

- 对解析曲线/曲面（直线/圆/椭/球/柱/锥/torus）：随机 1e4 组，C++ OCCT 程序输出 `NbExt+参值+距离`，Rust 侧同输入，距离差 < 1e-9、极值数一致、参值差 < 1e-8。
- 对一般 B样条：采样粗网格 + 牛顿收敛 vs OCCT 结果，距离差 < 1e-6。
- 回归：现有 1,485 测试保持绿；`point_*_extrema` 等 API 签名不变（内部换解析实现）。

## 依赖与纪律

- 每个 Agent 只写自己的新模块文件（各不同文件，无冲突）；`lib.rs` 由 orchestrator 统一注册。
- 匹配现有风格（`Result<T,String>`、`mod tests`、GpPnt/GpVec/dyn Curve/Surface）；复用 occt-math 的 newton.rs/roots.rs/optimize.rs。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target` 避免并行文件锁。
- 不碰 occt-topo 的 bop_builder/fillet_curved 等既有模块。
