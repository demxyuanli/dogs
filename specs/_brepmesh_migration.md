# BRepMesh 59 类深度迁移计划

> 覆盖 `_coverage.md` ⬜「BRepMesh 59 类深度」。参考源：`D:\source\occt-src\src\ModelingAlgorithms\TKMesh\`。
> 当前 Rust：occt-topo 有简化 meshing（shape_mesh/brepmesh/wireframe，~1,071 行）——本迁移按 OCCT 现代
> 增量网格架构补足深度：框架 → 数据模型 → UV Delaunay 三角化 → deflection 细化 → 愈合。

## OCCT 架构 → Rust 映射

| OCCT 包 | 类 | Rust 目标 | 角色 |
|---|---|---|---|
| IMeshTools | Parameters, Context, MeshAlgo, ModelAlgo, ShapeExplorer, ShapeVisitor, MeshBuilder, CurveTessellator, MeshAlgoType | `meshing/mod.rs` + traits | 框架参数/上下文/算法接口 |
| IMeshData | Model, Shape, TessellatedShape, Edge, Face, Wire, Curve, PCurve, ParametersList, Status, StatusOwner, Types | `meshing/model.rs` | 网格数据模型（含状态旗标） |
| BRepMeshData | Model, Edge, Face, Wire, Curve, PCurve | 并入 `model.rs` | 具体数据 |
| BRepMesh | IncrementalMesh, DiscretRoot | `meshing/incremental.rs` | 入口（Perform 管线） |
| BRepMesh | GeomTool, Deflection | `meshing/geom_tool.rs` | 点/切矢/法向 + 线性/角度偏差 |
| BRepMesh | EdgeDiscret, EdgeParameterProvider, CurveTessellator | `meshing/edge_discret.rs` | 边→多段线（deflection 自适应） |
| BRepMesh | FaceDiscret, FaceChecker, Classifier | `meshing/face_discret.rs` | 面→UV 顶点 + 面内分类 |
| BRepMesh | Delaun, DataStructureOfDelaun, SelectorOfDataStructureOfDelaun | `meshing/delaun.rs` | 2D Delaunay 三角化（核心，89KB） |
| BRepMesh | Vertex, Triangle, Circle, Edge, OrientedEdge, PairOfIndex | `meshing/delaun_types.rs` | Delaunay 数据结构 |
| BRepMesh | VertexTool/CircleTool + Vertex/CircleInspector | `meshing/delaun_index.rs` | 空间索引（包围盒网格） |
| BRepMesh | DelaunayBaseMeshAlgo, NodeInsertionMeshAlgo, DelaunayDeflectionControlMeshAlgo | `meshing/mesh_algo.rs` | 细化算法（线性/角度偏差分裂） |
| BRepMesh | ModelBuilder, ModelPreProcessor, ModelPostProcessor, ModelHealer | `meshing/model_pipeline.rs` | 构建/预处理/后处理/愈合 |
| BRepMesh | RangeSplitter 族（Default/UVParam/Undefined/BoundaryParams + 解析 Cylinder/Cone/Sphere/Torus/NURBS/Extrusion） | `meshing/range_splitter.rs` | UV 范围分割 |
| BRepMesh | MeshTool, ShapeTool, ShapeVisitor, DiscretFactory, 各 AlgoFactory | `meshing/tools.rs` | 工具/工厂（可选，插件架构锅炉板） |

## 交付波次

- **波 1（框架+数据+几何，~2,500 行）**：MeshParameters 全字段、MeshModel/Edge/Face/Wire/Curve/PCurve、
  Context traits、GeomTool、Deflection、EdgeDiscret、FaceDiscret、Classifier、IncrementalMesh+ModelBuilder+
  PreProcessor —— 端到端：形状→模型→边/面离散化（deflection 自适应，UV 顶点+边界）
- **波 2（Delaunay 核心，~2,500 行）**：Delaun（UV 空圆准则）、DataStructureOfDelaun、Selector、
  Vertex/Triangle/Circle/Edge/OrientedEdge/PairOfIndex、VertexTool/CircleTool+Inspector、
  MeshTool、ShapeTool —— 高质量三角化替代现有 UV 网格
- **波 3（细化+愈合+集成，~2,000 行）**：DelaunayBaseMeshAlgo、NodeInsertionMeshAlgo、
  DeflectionControlMeshAlgo（线性+角度偏差分裂细化）、ModelHealer、ModelPostProcessor、
  RangeSplitter 族、接入现有 shape_mesh/brepmesh（ShapeMesh 输出改由 OCCT 风格管线产生）

## 纪律
- 每个 agent 只写自己的 `meshing/` 新模块；不碰 lib.rs（orchestrator 统一注册）。
- 匹配 occt-topo 现有风格（Result<T,String>、`mod tests`、BRepTool/TopoShape）。
- 新增 `#[test]`，保持现有 1,350 测试绿。隔离 `CARGO_TARGET_DIR`。
- Delaun 移植重「行为等价」：UV 空圆插入、超三角形、顶点/圆索引，输出质量断言（最小角/共面）。
