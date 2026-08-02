# Phase 14 Plan — 基础补全：Poly + GProp + BRepGProp（≥5,000 行）

> 覆盖 `_coverage.md` 缺口：Poly 22%→~85%、GProp/BRepGProp 40%/23%→~80%（网格采样 → 精确 Gauss 积分）。
> 参考源：`D:\source\occt-src\src\FoundationClasses\TKMath\Poly\`（16 cxx 4,514 行）、
> `ModelingData\TKGeomBase\GProp\`（8 cxx 1,751 行）、`ModelingAlgorithms\TKTopAlgo\BRepGProp\`（13 cxx 5,435 行）、
> `FoundationClasses\TKMath\math\`（GlobOptMin + FRPR 补全）。
> 基线：5 crate · 121,337 行 · 1,578 测试（core 247 / math 202 / geom 139 / geom2d 72 / topo 918）。纯增量新模块，不动既有签名。

## 缺口 → 方案（各 agent 独立文件，无冲突）

| # | 缺口 | 方案 | 模块（新建） | 预估 |
|---|---|---|---|---|
| 1 | Poly 网格数据结构主干 | `Poly_Triangulation`+`PolygonOnTriangulation`+`ArrayOfNodes/UVNodes`（triangulation_full.rs）；`Poly_Coherent*` 族（coherent.rs）；`Poly_Connect`（connect.rs）；`Poly_MakeLoops`（make_loops.rs）；`Poly_MergeNodesTool`（merge_nodes.rs） | `occt-core/src/poly/*.rs`（5 文件） | 3,500 |
| 2 | GProp 全局属性框架 | `GProp_GProps` 累加器（Mass/CentreOfMass/MatrixOfInertia/StaticMoments/MomentOfInertia/RadiusOfGyration/PrincipalProperties）+ `PrincipalProps` + `PGProps/SelGProps/VelGProps/CelGProps/PEquation` + `GProp` 工具 | `occt-core/src/gprop/gprops.rs` | 1,400 |
| 3 | BRepGProp 精确属性 | `BRepGProp` 入口（Linear/Surface/Volume/VolumeGK）+ `Gauss`（自适应 Gauss）+ `Face` + `VinertGK/Vinert` + `MeshProps/MeshCinert` + `UFunction/TFunction/Cinert/Sinert` + `EdgeTool/Domain` | `occt-topo/src/brep_gprop_full.rs` | 4,300 |
| 4 | math 残余补全 | `math_GlobOptMin`（确定性全局优化 669 行）+ `math_FRPR`（共轭梯度 224 行） | `occt-math/src/globoptmin.rs` + `frpr.rs` | 700 |

**合计 ≈ 9,900 行**，远超 5,000 目标。

## OCCT 类 → Rust 映射要点

- **GProp_GProps**：`Add(item, density)` 合成（Density<=Resolution 抛错）；内部 `(g, loc, dim, inertia 3×3)`。CentreOfMass 用相对原点合成；MatrixOfInertia 平移到质心。PrincipalProperties 依赖 `math_Jacobi`（Rust `occt_math::jacobi`/`eigen_ext` 已有）。
- **GProp_PGProps/SelGProps/VelGProps/CelGProps**：点集 / 曲面 / 实体 / 曲线 的 GProps 子类，各自 `Add(点/三角形/曲面采样)` 规则。
- **BRepGProp 入口**：`LinearProperties(S)` 沿边积分；`SurfaceProperties(S)` 面/曲面参数域积分（d0/d1，二阶数值差分——Surface trait 无 d2，沿用 Phase 13 惯例）；`VolumeProperties(S)`/`VolumePropertiesGK(S)` 实体积分（GK=Gauss-Kronrod 版本）。入口签名与现有 `brep_gprop.rs` 不同——**新增** `brep_gprop_full.rs`，不替换旧模块。
- **BRepGProp_Gauss**：面自适应 Gauss 积分（4/8 点）+ 误差控制。
- **Poly_Triangulation**：`nodes: Vec<GpPnt>` + `uv_nodes` + `triangles` + `deflection`；`PolygonOnTriangulation` 是边→三角索引序列；`Coherent*` 纯数据结构（neighbor/opposite 链接）；`Poly_Connect` 顶点邻接；`MakeLoops` 边环闭合（图回路）；`MergeNodesTool` 焊点。
- **测试基线**：对拍 `gprop_analytic.rs` 的盒（面积 52/体积 24）/球/圆柱/圆锥/圆环精确结果；Poly 用 BRepMesh 输出自洽性；math GlobOptMin 用 Rosenbrock 已知最小值。

## 交付波次

- **波 1（Poly，~3,500 行，1 agent）**：`occt-core/src/poly/` 下 5 个新文件 + 每文件 `mod tests`。
- **波 2（GProp，~1,400 行，1 agent）**：`occt-core/src/gprop/gprops.rs`。
- **波 3（BRepGProp，~4,300 行，1 agent）**：`occt-topo/src/brep_gprop_full.rs`。
- **波 4（math 残余，~700 行，1 agent）**：`occt-math/src/globoptmin.rs` + `frpr.rs`。

## 差分验证（每波必做）

- GProp/BRepGProp：与 `gprop_analytic.rs` 的盒/圆柱/球/圆锥/圆环精确结果对拍（面积/体积/质心容差 < 1e-9 解析面）；非解析面与旧 `brep_gprop` 网格结果 < 1% 一致性。
- Poly：三角化节点/三角计数/包围盒；Coherent 邻接自洽；MakeLoops 环闭合（单位盒 6 面各自闭合）。
- math：`math_GlobOptMin` 在 Rosenbrock/Rastrigin 已知最小值附近收敛；`FRPR` 与 `bfgs.rs` 结果一致。

## 依赖与纪律

- 每个 Agent 只写自己的新模块文件（各不同文件，无冲突）；`lib.rs`/`poly/mod.rs`/`gprop/mod.rs` 由 orchestrator 统一注册。
- 匹配现有风格（`Result<T,String>`、`mod tests`、GpPnt/GpVec/GpMat、dyn Curve/Surface）；GProp 用 `occt_math::jacobi`/`eigen_ext`。
- 新增 `#[test]`；隔离 `CARGO_TARGET_DIR=<模块目录>/target` 避免并行文件锁。
- 不碰 bop_builder/fillet_curved/meshing 等既有模块；不重写 `brep_gprop.rs`/`gprop_analytic.rs`。
