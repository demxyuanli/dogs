# OCCT 面网格 GT 表（`BRepMesh_IncrementalMesh`）

> ⚠️ **必读（round 138 重大订正）**：本文开头那张「deflection 0.1」表**参数用错了**，不能用来判断端口保真度。导出路径 `export_data_obj` 实际用的是 **`Prs3d` 相对偏转**（`lin = maxComp(bbox) × 0.001 × 4`，`brep_exchange.rs:56-77`）+ **角偏 20°**（`:95`），且走的是 **Delaunay 管线**（`:88-94` 注明 legacy 四叉树/UV 网格回退已按 audit A17 删除）。**按端口真实参数重跑后，15 个文件里 14 个与 OCCT 逐位一致**（见下方「端口真实参数下的 GT 表」）。

> 用途：T-54+T-90（网格保真度）的**验收基准**，以及 §2 `export_data_obj` 期望值重订的依据。
> 生成命令（探针 `--mesh` 模式，round 129 加入）：
> ```
> $env:THIRDPARTY_DIR='D:\source\occt-8.0.0\3rdparty-vc14-64'
> cmd /c "cd /d D:\source\occt-8.0.0 && call env.bat vc14 64 && D:\source\repos\dogs\specs\occt_probe\occt_probe.exe D:\source\repos\dogs\data\<f>.step --mesh 0.1"
> ```
> 逐面明细也在探针输出里（`FACE i nodes=.. triangles=..`）。

## 全表（OCCT vs 端口 `export_data_obj` 基线）

| `data/*.step` | OCCT faces | OCCT nodes | OCCT tris | 端口 v | 端口 f | 判定 |
|---|---|---|---|---|---|---|
| Cone | 2 | 103 | 141 | 195 | 301 | 端口偏多 |
| **Cube** | 6 | **24** | **12** | **24** | **12** | ✅ **逐位一致** |
| Cylinder | 3 | 106 | 100 | 146 | 140 | 端口偏多 |
| **Extrusion** | 6 | **24** | **12** | **24** | **12** | ✅ **逐位一致** |
| **HoledPlate** | 32 | **180** | **128** | **180** | **128** | ✅ **逐位一致** |
| linkrods | 37 | 2184 | 2928 | 3494 | 5078 | 端口偏多 |
| Offset | 26 | 496 | 556 | 712 | 892 | 端口偏多 |
| **OffsetPlaneHoleEdge** | 10 | **48** | **32** | **48** | **32** | ✅ **逐位一致** |
| rev | 6 | 76 | 64 | 104 | 92 | 端口偏多 |
| screw | 10 | 652 | 944 | 600 | 790 | 端口**偏少** |
| Shape-1 | 60 | 2892 | 3724 | 3343 | 4336 | 端口偏多 |
| Shape-2 | 31 | 11212 | 20286 | 3105 | 4792 | 端口**仅 ~1/3.6** |
| Shape | 11 | 18733 | 36444 | 6150 | 11372 | 端口**仅 ~1/3** |
| Sphere | 1 | 273 | 516 | 642 | 1244 | 端口 ~2.4× |
| Torus | 1 | 810 | 1508 | 1369 | 2592 | 端口 ~1.7× |

## 从表里能读出的结论（round 135）

1. **4 个模型已逐位一致**（Cube / Extrusion / HoledPlate / OffsetPlaneHoleEdge）——它们都**以平面面为主**（含带孔平面，如 HoledPlate 32 面）⇒ 端口的**平面路径**（`planar_polygon_triangulate`：「质心径向排序 + 扇形 + 耳切带孔桥接」）在这些用例上**恰好**与 OCCT 的「边离散化 + 约束 Delaunay」给出相同片数。
2. **曲面面系统性偏多**：Cone/Cylinder/rev（~1.3×）、Sphere（2.4×）、Torus（1.7×）、linkrods（1.6×）、Offset（1.4×）——与「端口用**递归 UV 四叉树**顶替 OCCT 的**边离散化 + 约束 Delaunay**」一致（四叉树在中纬度处铺得比 Delaunay 密）。
3. **`Shape` / `Shape-2` 反向、且差 3 倍以上**（端口**偏少**）：这两例面数多（11 / 31）且含大量**带内环/非凸**的平面面——端口对这类面会**回落到 UV 网格**（`planar_polygon_triangulate` 自述「非凸即 None ⇒ 回落 `face_to_triangles`」），且不生成 OCCT 的**内部边**（`BRepMesh` 的 `AddInternalEdges`），所以片数远低于 OCCT。**这是当前最大的单点缺口**。
4. `screw` 偏少（652/944 vs 600/790，约 0.95×）——它多为平直/柱面混合，偏离最小，量级上更像「点数公式差一点」而非算法差。

⇒ 给 T-54+T-90 的优先级建议：**先把 `Shape`/`Shape-2` 这一类（带内环/非凸平面 + 内部边）补上**（收益最大），再统一曲面面的 Deflection→点数公式与约束 Delaunay。

## 逐边离散化 GT（`--edges` 模式，round 136 新增）

> 用途：T-54+T-90 的**第 1 步**（把已移植未接线的 `compute_nb_samples*` 族接进网格管线）的**逐边验收基准** —— 它要复现的正是「每条边放几个点」。
> 命令同上，把 `--mesh` 换成 `--edges`（可选 deflection 参数）。输出 `EDGE face=<面序> i=<边序> nodes=<边上节点数> degenerated=<0/1>` 逐条 + `TOTAL faces/edge_occurrences/edge_nodes/without_polygon`。注意：同一几何边会被它相邻的**每个面各列一次**（`Poly_PolygonOnTriangulation` 是 (edge, face) 对的），所以 `edge_nodes` 是「面×边」的总和。

**示例（`Cylinder.step`，deflection 0.1）**：
```
EDGE face=0 i=2 nodes=27 degenerated=0      <- 底圆（每圈 27 点）
EDGE face=0 i=3 nodes= 2 degenerated=0      <- 缝边 occurrence A
EDGE face=1 i=4 nodes=27 degenerated=0      <- 顶圆
EDGE face=2 i=5 nodes=27 degenerated=0      <- 侧柱面里的圆
TOTAL faces=3 edge_occurrences=6 edge_nodes=112 without_polygon=0
```
⇒ 半径 4、高 10 的圆柱在 deflection 0.1 下 OCCT 给**每圆 27 点**；端口的 `compute_nb_samples*`（对应 `BRepMesh_EdgeDiscret`/`BRepMesh_Deflection` 的 nb-points 公式）必须给出同一数字，之后约束 Delaunay 才能给出同样的 100 三角形。

其它文件的逐边表可用同一条命令自取（`--edges 0.1`），必要时逐面对齐。

## ✅ 端口真实参数下的 GT 表（round 138 订正，**这是判断保真度的唯一正确表**）

参数：对每个文件 `lin = maxComp(bbox) × 0.004`（= `Prs3d::GetDeflection`，`brep_exchange.rs:56-77`）、角偏 **20° = 0.349066 rad**（`brep_exchange.rs:95`）。命令：`occt_probe.exe <file> --mesh <lin> 0.349066`。

| `data/*.step` | bbox maxComp | `lin` | OCCT nodes/tris | 端口 v/f | 判定 |
|---|---|---|---|---|---|
| Cone | 10 | 0.04 | 195 / 301 | 195 / 301 | ✅ 逐位一致 |
| Cube | 10 | 0.04 | 24 / 12 | 24 / 12 | ✅ |
| Cylinder | 10 | 0.04 | 146 / 140 | 146 / 140 | ✅ |
| Extrusion | 20 | 0.08 | 24 / 12 | 24 / 12 | ✅ |
| HoledPlate | 102.144 | 0.408576 | 180 / 128 | 180 / 128 | ✅ |
| linkrods | 5.115 | 0.02046 | 3494 / 5078 | 3494 / 5078 | ✅ |
| Offset | 14 | 0.056 | 712 / 892 | 712 / 892 | ✅ |
| OffsetPlaneHoleEdge | 14 | 0.056 | 48 / 32 | 48 / 32 | ✅ |
| rev | 10 | 0.04 | 104 / 92 | 104 / 92 | ✅ |
| screw | 42.3 | 0.1692 | 600 / 790 | 600 / 790 | ✅ |
| Shape-1 | 77 | 0.308 | 3343 / 4336 | 3343 / 4336 | ✅ |
| **Shape-2** | 150 | 0.6 | 3099 / 4780 | **3105 / 4792** | ❌ 差 6 节点 / 12 三角（**0.2%**） |
| Shape | 112.106 | 0.448424 | 6150 / 11372 | 6150 / 11372 | ✅ |
| Sphere | 10 | 0.04 | 642 / 1244 | 642 / 1244 | ✅ |
| Torus | 25.978 | 0.103912 | 1369 / 2592 | 1369 / 2592 | ✅ |

**结论（取代 round 129/135 的判断）**：
1. **端口的网格管线在相同参数下与 OCCT 逐位一致（14/15）** ⇒ round 129/135 说的「系统性不符」是**我用错参数**（拿 0.1/0.5rad 的 GT 去比 0.004·bbox/20° 的端口输出）造成的假象；**T-90 不成立**（应关闭）。
2. `Shape-2` 的 0.2% 差**已由源码自注给出方向**：`brep_exchange.rs:50-52` 写明「B-spline 面的 STEP `SURFACE_CURVE` pcurve 被丢 ⇒ 边界错 ⇒ Delaunay bbox 漂移」，即 **STEP 导入侧缺口**，不是 BRepMesh 侧。
3. **T-54 也必须重新评估**：导出路径**不再**用 `brepmesh.rs` 的四叉树、也不用 `wireframe.rs` 的平面耳切（`brep_exchange.rs:88-94`：这两个非 OCCT 回退已按 audit A17 删除）——它们只在**失败回退**里存在。因此「移植约束 Delaunay 以提高保真度」这个前提**不成立**；`brepmesh.rs`/`wireframe.rs` 的处置应改为「确认回退是否可达 + 可达则接线/不可达则摘除」，而不是大规模替换。

## 附：旧表（参数 0.1 / 0.5rad，**已作废，仅留档对照**）
