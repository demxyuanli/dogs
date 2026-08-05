# Phase 22 Plan — STEP→Mesh→OBJ 面积/几何对齐（修复 Cone/Cylinder/rev）

> **状态更新（2026-08-05）**：缺口 2（Cylinder u 域塌缩）已修——`wireframe::face_uv_bounds` 对整周期边连续解包 + 全周期检测（wireframe.rs）。incremental（brep_to_obj）面积 133.9→149.35（解析 150.8）；Cylinder bbox 对拍过。**剩余 P2**：`face_to_triangles` 对平面 cap 铺外接正方形（面积 16.0 vs 12.57）→ 复用 `planar_face_mesh` 边界多边形逻辑（见下方波 A2'）。

> 审查结论（已核实 OCCT 源码）：
> - **面积错是 STEP 几何读入错**，不是网格——Cone 解析面积（brep_gprop_full）51.4（应 185），mesh 4.7。OCCT 面积用解析 `BRepGProp_Face`（BRepGProp.cxx L225-230），非网格。
> - 波 A（face_uv_bounds 用 pcurve）**OCCT 证实**：`BRepTools::AddUVBounds`（L179）用 `BRep_Tool::CurveOnSurface` 取边 pcurve + `BndLib_Add2dCurve` 加包围盒。
> - Cone 根因已定位：**make_pcurve_full 对锥面投影 v 域只 0.9（应 10）**——母线 (0,0,0)→(4,0,10) 在锥面的 pcurve 投影错。
> - 波 C 修正：**不要改 surface_area 用网格**——OCCT 用解析 `BRepGProp_Face`（Phase 14 已移植 brep_gprop_full）。应修 STEP 读入使解析面积正确。
> 参考源：OCCT `BRepTools::UVBounds`/`AddUVBounds` + `BRepGProp_Face` + `StepToTopoDS`。

## 缺口清单（审查定位）

| # | 缺口 | 证据 | 根因 | 波及 |
|---|---|---|---|---|
| 1 | Cone 侧壁 v 域错 | 侧壁面积 1.0（应 135） | face_uv_bounds 对锥面（1 母线+顶圆）边投影不完整 | Cone |
| 2 | Cylinder 侧壁 mesh_shape 面积 16 | incremental 109 对 | face_to_triangles 对圆柱侧壁 UV 网格覆盖错 | Cylinder/surface_area |
| 3 | rev 1/4 圆环柱锯齿 | z 层 r=5/10 交替 | 6 面体（2 圆柱+4 平面）面方向/曲面解析错 | rev |
| 4 | Shape 复合实体 0 形状 | 解析到中段 | parse_records 复合实体跨行 type_name 空 | Shape/1/2 |
| 5 | 曲面侧壁 v 细分过密 | Cylinder 侧壁 16 z 层 | OCCT v 直线不细分，我们弦偏差细分 | Cylinder/Cone 密度 |
| 6 | Cone u 段少 | ratio 0.23 | cone_radius_at/arc_angular_step 对锥标定 | Cone 密度 |
| 7 | 面积 API 语义 | OCCT 用解析 BRepGProp_Face | 验证用 brep_gprop_full 解析面积（非 mesh） | 面积/体积 API |
| 8 | 面积对拍门禁缺失 | 仅 bbox | 无面积对比测试 | 质量门禁 |

## 波次方案

### 波 A：STEP 曲面边界修复（缺口 1/2/3，~3,000 行）
- **A1 face_uv_bounds 对锥/柱**：修复边投影对单一曲面（锥 1 母线+顶圆）的 v 域推导。OCCT `BRepAdaptor_Surface::UVBounds` 用**所有边界边 pcurve 的 UV 包围盒**，需补锥面母线+圆的 pcurve 完整覆盖。
- **A2 face_to_triangles 对曲面侧壁**：Cylinder 侧壁 mesh_shape 面积 16 修复——face_uv_bounds（波 A 已修 pcurve）后 UV 网格覆盖完整。（u 域塌缩已修，见状态更新。）
- **A2'（P2，新增）face_to_triangles 平面 cap 铺正方形**：平面面应走边界多边形（复用 `brepmesh::planar_face_mesh` 的收集+径向排序+凸扇形逻辑，上提为共享函数），非 UV 网格。cap 面积 16.0→12.57。mesh_shape 总面积 157.6→~150.8。
- **A3 rev 1/4 圆环柱**：诊断 6 面体曲面拼接（2 圆柱 90° + 4 平面），修面方向/法向。
- 门禁：面积对拍（Cone ~175、Cylinder ~151、rev r=10 圆柱面积 ~628）。

### 波 B：曲面网格标定（缺口 5/6，~1,500 行）
- **B1 解析曲面 v 方向**：OCCT 对圆柱/锥 v（直线母线）不细分——adaptive_face_mesh 对解析曲面 v 方向用 1-2 段（非弦偏差）。
- **B2 Cone u 段**：cone_radius_at 采样修正，u 段按底半径标定。
- 门禁：Cylinder/Cone 密度 ratio → ~1（面数接近 OCCT）。

### 波 C：解析面积验证 + 面积门禁（缺口 7/8，~800 行）
- **C1 解析面积对齐**：OCCT 面积用解析 `BRepGProp_Face`（非网格）。用 `brep_gprop_full::surface_properties`（Phase 14 已移植）作为面积验证——修 STEP 读入后解析面积应正确（Cone 51.4→185、Cylinder→151）。若 brep_gprop_full 对曲面有缺口则修它。
- **C2 面积对拍门禁**：step_obj_parity 加解析面积断言（Cube 600、Sphere 314、Cylinder 151、Cone 175）。
- 门禁：解析面积 ratio 0.9-1.1（比 mesh 面积更准）。

## 验证

- 每波：面积对拍（新增面积断言）作为验收。
- 全量回归：提交前一次（测试纪律）。
- 保留 bbox 门禁 + 密度记录。

## 依赖与纪律

- 波 A 是根因（Cone/Cylinder 面积错），优先。
- 波 B 依赖波 A（v 域修正后标定才准）。
- 波 C 独立可并行。
- 匹配现有风格；对拍门禁是每波验收。
