# OCCT 面网格 GT 表（`BRepMesh_IncrementalMesh`，deflection 0.1）

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
