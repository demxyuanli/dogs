# Phase 22 Plan — STEP→Mesh→OBJ 面积/几何对齐（修复 Cone/Cylinder/rev）

> 审查结论：STEP→OBJ 管线读入/网格化存在几何缺口，面积对比暴露：
> Cube/Sphere 正确，**Cone 面积 4.7（应 175）错 37×、Cylinder 16（应 151）错 9×、rev 锯齿错**、Torus 偏 18%。
> 根因在 STEP 曲面边界解析（face_uv_bounds）与网格化（face_to_triangles 对曲面侧壁）。
> 参考源：OCCT `BRep_Tool::UVBounds`/`BRepAdaptor_Surface` + `StepToTopoDS` 曲面构造。

## 缺口清单（审查定位）

| # | 缺口 | 证据 | 根因 | 波及 |
|---|---|---|---|---|
| 1 | Cone 侧壁 v 域错 | 侧壁面积 1.0（应 135） | face_uv_bounds 对锥面（1 母线+顶圆）边投影不完整 | Cone |
| 2 | Cylinder 侧壁 mesh_shape 面积 16 | incremental 109 对 | face_to_triangles 对圆柱侧壁 UV 网格覆盖错 | Cylinder/surface_area |
| 3 | rev 1/4 圆环柱锯齿 | z 层 r=5/10 交替 | 6 面体（2 圆柱+4 平面）面方向/曲面解析错 | rev |
| 4 | Shape 复合实体 0 形状 | 解析到中段 | parse_records 复合实体跨行 type_name 空 | Shape/1/2 |
| 5 | 曲面侧壁 v 细分过密 | Cylinder 侧壁 16 z 层 | OCCT v 直线不细分，我们弦偏差细分 | Cylinder/Cone 密度 |
| 6 | Cone u 段少 | ratio 0.23 | cone_radius_at/arc_angular_step 对锥标定 | Cone 密度 |
| 7 | surface_area 用 face_to_triangles | Cylinder 16 vs incremental 109 | brep_gprop::surface_area 未用 incremental | 面积/体积 API |
| 8 | 面积对拍门禁缺失 | 仅 bbox | 无面积对比测试 | 质量门禁 |

## 波次方案

### 波 A：STEP 曲面边界修复（缺口 1/2/3，~3,000 行）
- **A1 face_uv_bounds 对锥/柱**：修复边投影对单一曲面（锥 1 母线+顶圆）的 v 域推导。OCCT `BRepAdaptor_Surface::UVBounds` 用**所有边界边 pcurve 的 UV 包围盒**，需补锥面母线+圆的 pcurve 完整覆盖。
- **A2 face_to_triangles 对曲面侧壁**：Cylinder 侧壁 mesh_shape 面积 16 修复——face_uv_bounds（波 A 已修 pcurve）后 UV 网格覆盖完整。
- **A3 rev 1/4 圆环柱**：诊断 6 面体曲面拼接（2 圆柱 90° + 4 平面），修面方向/法向。
- 门禁：面积对拍（Cone ~175、Cylinder ~151、rev r=10 圆柱面积 ~628）。

### 波 B：曲面网格标定（缺口 5/6，~1,500 行）
- **B1 解析曲面 v 方向**：OCCT 对圆柱/锥 v（直线母线）不细分——adaptive_face_mesh 对解析曲面 v 方向用 1-2 段（非弦偏差）。
- **B2 Cone u 段**：cone_radius_at 采样修正，u 段按底半径标定。
- 门禁：Cylinder/Cone 密度 ratio → ~1（面数接近 OCCT）。

### 波 C：输出 API + 面积门禁（缺口 7/8，~800 行）
- **C1 brep_gprop::surface_area/volume** 用 incremental_mesh（或 export_mesh），对齐面积。
- **C2 面积对拍门禁**：step_obj_parity 加面积断言（Cube 600、Sphere 314、Cylinder 151、Cone 175）。
- 门禁：面积 ratio 0.9-1.1。

## 验证

- 每波：面积对拍（新增面积断言）作为验收。
- 全量回归：提交前一次（测试纪律）。
- 保留 bbox 门禁 + 密度记录。

## 依赖与纪律

- 波 A 是根因（Cone/Cylinder 面积错），优先。
- 波 B 依赖波 A（v 域修正后标定才准）。
- 波 C 独立可并行。
- 匹配现有风格；对拍门禁是每波验收。
