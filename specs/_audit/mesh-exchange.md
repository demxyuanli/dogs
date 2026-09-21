# 网格化 / 数据交换 移植忠实度审查（2026-09-20）

## 结论摘要

范围：`crates/occt-topo/src/meshing/**`、`brepmesh.rs`、`wireframe.rs`、`shape_mesh.rs`、`mesh*.rs`、`step/**`、`iges.rs`、`brep_exchange.rs`、`gltf.rs`、`vrml.rs`、`brep_bnd_lib.rs`、`geom_bnd_lib_*.rs`、`brep_uv_bounds.rs`、`occt-core/src/io/{obj,ply,stl}.rs`。**只读，未改任何源码。**

**live 路径**：四道门禁全走 `brep_exchange::brep_to_obj` → `export_mesh` → `IncrementalMesh::from_deflection`（`brep_exchange.rs:88-97`）。实测本轮 `export_data_obj` 与 OCCT 基线（`data/_occ_ref_export.tcl`：`incmesh s lin -angular 20`，`lin = maxComp(bbox)*0.001*4`）**顶点/三角数逐模型一致**：Cube 24/12、Cone 195/301、Cylinder 146/140、Sphere 642/1244、Torus 1369/2592、HoledPlate 180/128、Shape-1 3343/4336、rev 104/92、linkrods 3494/5078、screw 600/790；仅 Shape 差 1/2、Shape-2 差 7/14、ATU01038 差 335/426 ⇒ **活的 Delaunay 管线已高度对齐**，高危集中在**回退路径 / 类型判定 / 边参数推导 / UV 域推导**。

- **自创** 15 ｜ **等价替换未登记** 8 ｜ **已登记（备案）** 6
- **最高危**：`step/format.rs:114` `classify_curve` —— OCCT 的 `StepToTopoDS_TranslateEdge` **没有任何曲线族判定**，族由 STEP 实体类型经 `StepToGeom::MakeCurve` 固定；本端口用 6 点二阶差分猜族，读+写两侧共用。
- **次危（live 未登记）**：`meshing/model_healer.rs:279-290` `adjustSamePoints` 退化支路左右端接反。

## 发现

**1｜`step/format.rs:114-144` 采样分类器决定曲线族（读+写共用）｜自创**
`:121` `let n = 6;` ／ `:132` `if (max_d2 - min_d2) / max_d2 < 0.02 {`。live **是**（读 `step/read_geometry.rs:573`；写 `step/format.rs:670`、`step/write_context.rs:157`）。OCCT `StepToGeom/StepToGeom.cxx:1335` `if (SC->IsKind(STANDARD_TYPE(StepGeom_Line)))`…`:1343`；写侧 `GeomToStep/GeomToStep_MakeCurve.cxx:54` —— `6`/`1e-9`/`0.02` 在任何 `.cxx` 中都不存在。影响**是**（四道门禁）；注释已自承近圆椭圆误判为圆（`step/read_geometry.rs:552-556`，ATU01038 端点错位 7.65e-2）。建议：改用 `occt_geom::Curve` 真实类型标签（`StepToGeom.cxx:1335-1349`），或标「未移植」。

**2｜`brep_exchange.rs:56-98` `deflection` 被 bbox 相对 Prs3d 偏转顶掉 + 三级静默回退｜未登记**
`:87` `let lin = prs3d_get_deflection(shape, deflection);`；`:96-97` `.or_else(|| brepmesh::incremental_mesh(...))` / `.unwrap_or_else(|| shape_mesh::mesh_shape(...))`。live **是**（门禁出口）。数值对得上 `Prs3d.hxx:71` + `Prs3d_Drawer.cxx:95-96`，但**所引 RWObj 链路不存在**：`RWObj_CafWriter` 不网格化，`RWMesh_FaceIterator.cxx:87` `myPolyTriang = BRep_Tool::Triangulation(...)`；相对偏转属显示管线 `StdPrs_ToolTriangulatedShape.cxx:170`；空三角化是**跳过**（`RWMesh_FaceIterator.cxx:89`），不换算法。影响**是**，`brep_to_obj(shape, 0.1)` 的 `0.1` **从不生效**（仅 void-bbox 兜底），密度由包围盒推出。建议：删 `:96-97` 两级回退；形参改名 `maximal_chordial_deviation`，同步 `export_data_obj.rs:88`、`tests/step_obj_parity.rs:91` 的注释。

**3｜`wireframe.rs:392-448` 曲面面片用「偏转取整 UV 栅格 + 64 硬上限」替代 Delaunay｜自创**
`:407-408` `let nu = ((du / def).ceil() as usize + 1).clamp(3, 64);`（`nv` 同）。live **是**（`brep_tools.rs:75`、`brep_exchange.rs:97`、`shape_mesh.rs:29`；`mesh_shape` 的 gltf/vrml/hlr/render_svg/viz_scene/geometry_query/shape_metrics/brep_scene 调用者无条件吃它）。OCCT **未找到**：曲面只走 `MeshAlgoFactory::GetAlgo`（`BRepMesh_MeshAlgoFactory.cxx:60-117`）；`GeomTool::CellsCount`（`BRepMesh_GeomTool.cxx:465-512`）算 Delaunay **内部栅格密度**，无 64 上限。影响**是**，触发即整形状换三角化器（`step_obj_area` 的 4~5% 容差即为此设）。建议：删该分支（只留平面精确 + Delaunay），或按 `BRepMesh_GeomTool.cxx:465-512` 移植。

**4｜`wireframe.rs:257-380` 平面面片用「质心角度排序 + 耳切 + 最近点对桥洞」｜自创**
`:297-304` `let order = |l: &mut Vec<GpPnt>| { l.sort_by(|p,q| { /* atan2 around centroid */ }); };`；另 `:317` `ear_clip`、`:352` `bridge_holes`。live **是**（平面面片主路径）。OCCT **未找到**：「按质心角度排序」只在凸边界等价，凹边界会跨面自交；OCCT 用**约束 Delaunay**（`BRepMesh_DelaunayBaseMeshAlgo.cxx` + `BRepMesh_Delaun.cxx`）。影响**是**（`step_obj_area`、HoledPlate/OffsetPlaneHoleEdge）。建议：把边界链作为约束边交给 `meshing::delaun::Delaun`。同类（已登记，非 live）：`face_discret.rs:573-609` `interior_grid`（`nu=(du/def).ceil()+1` + `4096` 预算 + `:588` `scale=(budget/(nu*nv)).sqrt()`）、`triangulator.rs:65` `MAX_REFINE_PASSES = 4`（OCCT `BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx:85` `aIterationsNb = 11`）、`fast_discret.rs:107` `weld_tol = deflection*0.01`（OCCT 8.0 `BRepMesh_FastDiscret.hxx:19-27` 只剩 deprecated typedef）。

**5｜`incremental_mesh/discret_root.rs:369,460-467` `WIREFRAME_FALLBACK_RATIO_MAX = 0.10` 整形状换算法｜自创**
`:369` `const WIREFRAME_FALLBACK_RATIO_MAX: f64 = 0.10;` → `:462-466` 超限 `return Err`，`perform:210` 整体转 UV 栅格（发现 3）。live **是**。OCCT **未找到**：`BRepMesh_IncrementalMesh::Perform` 无失败率阈值；单面失败只置 `IMeshData_Failure`，其余面照常（`BRepMesh_BaseMeshAlgo.cxx:52-59`）。影响**是**。建议：删阈值，逐面标记失败。
**处置**：阈值已删（T-55，2026-09-20）。**逐面 UV 栅格回退（`wireframe_face_triangulation`，`incremental_mesh/discret_root.rs:539`）仍在**：T-69 第 1 轮（第 41 轮）实测 T0M 有 166 面走它；其中"单闭合边 loop 的管面"在 OCCT 参考 `data/occ-ref/T0M.obj` 里也没有网格（参考逐面写顶点不去重：20958/51160 条重复 `v`；失败管顶圆按半径 9.75±0.01 数只出现一次且只被一个圆盘盖面使用），但同日删该回退的实测被 `step_obj_parity` **13/14** 挡下（T0M `min[2]` 短 0.33，参考另覆盖 `(−0,−11.823,−424.742)`）⇒ 还有一类失败面 OCCT 会网格化；先补该类控制流，再删回退。

**5b｜`delaun/polygon_meshing.rs:346-375` `decomposeSimplePolygon` 自造「先删三角形再 AddElement」｜自创（未登记）**
`:348` `if self.mesh_data.elements_connected_to(id).extent() < 2 { continue; }`；`:362-365` `if element.link_at(k).abs() == id && (element.link_at(k) > 0) == is_forward { ... self.delete_triangle(elem_id, &mut loop_edges);`。live **是**（`Delaun` 内部，所有面）。OCCT `BRepMesh_Delaun.cxx:2259-2274` 直接 `AddLink(...)`+`addTriangle(...)`，**无任何删除**；OCCT 在此退化为第三条连接时 `BRepMesh_PairOfIndex.hxx:41` 抛 `Standard_OutOfRange`，被 `BRepMesh_BaseMeshAlgo.cxx:62` 空 catch 吞掉 ⇒ OCCT 该面**没有三角形**。影响**是**，OCCT 不生成三角网的面本端口删邻三角形后继续建网（为避免 `delaun_types.rs:470` 的 append panic 打的补丁）。建议：让该面失败置 `IMeshData_Failure`（不要改网格）。

**6｜`step/read_geometry.rs:508-604` 边参数用 1e-3 顶点距离替代 `V1.IsSame(V2)` 拓扑判据｜自创**
`:516` `const PRECI: f64 = 1e-3;` ／ `:522` `if p1.distance(p2) < PRECI {`；另 `:618-622` `CurveKind::Other => { ... (0.0, 1.0) }`。live **是**。OCCT `StepToTopoDS_TranslateEdge.cxx:438` `if (V1.IsSame(V2))`；`:443` `temp1 = sac.Project(C1, pnt1, preci, pproj, U1, false);`；`ShapeAnalysis_Curve.cxx:376-377` `case GeomAbs_Hyperbola:` → `ElCLib::Parameter(...)`。影响**是**：①任何两端点距离 <1e-3 的有界非周期曲线（B-spline/Bezier/trimmed）整条结点域被当边域并跳过 Project；②`HYPERBOLA`（`occt-geom/src/hyperbola.rs:17` first/last = ±inf）边域落到 `(0,1)`。建议：改判顶点同一性，补 `ElCLib` 精确臂（`ShapeAnalysis_Curve.cxx:376-400`）。

**6b｜`meshing/model_healer.rs:279-290` `adjustSamePoints` 退化支路左右端接反｜自创（未登记）**
`:279` `if curr_prev_side == curr_next_side {`；`:285-286` `set_pcurve_end(..., curr_prev_side, prev_val)` / `set_pcurve_end(..., other_curr_side, other_next_val)`。live **是**（`incremental_mesh/discret_root.rs:348,587`）。OCCT `BRepMesh_ModelHealer.cxx:491-512` `if (aPrevSqDist - aNextSqDist > gp::Resolution()) { adjustSamePoints(aCurrNextUV, aNextUV, aCurrPrevUV, aPrevUV, ...) } else { adjustSamePoints(aCurrPrevUV, aPrevUV, aCurrNextUV, aNextUV, ...) }` + `hxx:143-151`（两支都是「共享端 ← next 值，另一端 ← prev 值」）；本端口连 `aPrevSqDist - aNextSqDist` 判定都丢了。影响**是**，pcurve 端点吸附到不同邻居 ⇒ `AddNode` 焊接点与闭合边不同。建议：逐行对齐 `BRepMesh_ModelHealer.cxx:491-512` + `hxx:143-151`。

**7｜`step/read_geometry.rs:740-760` 采样圆回退（3 点外心 + 点积 atan2）取代 `Project`｜自创**
`:748` `let pa = curve.d0(lo);` ／ `:751` `let center = circle_center3(&pa, &pb, &pc).unwrap_or_else(GpPnt::zero);`。live **是**（`step/read_geometry.rs:583` `CurveKind::Circle`；在发现 1 把周期 B-spline/近圆椭圆判成 Circle 时可达）。OCCT `ShapeAnalysis_Curve.cxx:160` `GeomAdaptor_Curve GAC(C3D, uMin, uMax);` / `:200` `return ProjectAct(GAC, P3D, preci, proj, param);`。影响**是**。建议：全走 `shape_analysis_project`（已有），或移植 `ElCLib::CircleParameter` + `AdjustPeriodic`。

**8｜`wireframe.rs:462-599` `face_uv_bounds` 用 9 点 pcurve 采样代替精确 2D 包围盒｜自创**
`:502-503` `for k in 0..=8 { let t = a0 + (a1 - a0) * k as f64 / 8.0;`；另 `:553` `if cu_span > 1e-9 {`、`:567` `if cu_max - cu_min >= u_period - 1e-6 {`。live **是**（`brepmesh.rs:115`、`wireframe.rs:399`）。OCCT `ModelingData/TKBRep/BRepTools/BRepTools.cxx:185` `BndLib_Add2dCurve::Add(aC2D, aT1, aT2, 0., aBoxC);`（**精确**包围盒，B-spline 走控制多边形）；周期核实 `:210-268` 用 `aS->Value()` 精确判别；无采样、无 `1e-9`/`1e-6`。影响**是**，漏掉 pcurve UV 极值 ⇒ 窗口偏斜 ⇒ 采样点跑出 trim（board 记 `a3n00` y 溢出 +60.04）。建议：逐行移植 `BRepTools::AddUVBounds`（`BRepTools.cxx:172-330`），复用 `geom_bnd_lib_curve2d`。

**9｜`meshing/range_splitter/param_set.rs:86-133` 面类型由「周期标志 + 半径采样」猜｜自创**
`:112` `(true, false) if v_fin => SurfaceType::Sphere,`；`:132` `(r0 - r1).abs() <= 1e-7 * r0.abs().max(r1.abs()).max(1e-7)`。live **是**（`node_insertion.rs:889`、`range_splitter/nodes.rs:178,203`、`model_builder/preprocessor.rs:170`、`edge_discret.rs:912`）。OCCT `BRepMesh_FaceDiscret.cxx:112` `myAlgoFactory->GetAlgo(aDFace->GetSurface()->GetType(), myParameters)`；`BRepMesh_MeshAlgoFactory.cxx:64` `switch (theSurfaceType)`（11 个 case）；`GeomAdaptor_Surface.cxx:422-529` 亦为 `DynamicType()` 精确比较。影响**是（静默）**，选错 splitter 即换掉整张内部节点网格（Sphere 的 0.7 交错格 vs Cylinder 的 du 列），并决定 `factory_uses_deflection_control`；**任何 V 有限的柱/锥会被判成 Sphere**。建议：`classify_surface` 直接映射真实类型标签，删采样臂。

**10｜`brepmesh.rs:38,167-203` 旧四叉树网格器的深度上限、`/96`、`/4` 全为自造｜自创**
`:38` `const MAX_FACE_DEPTH: usize = 9;` ／ `:174-175` `Some(du.max((u1 - u0) / 96.0))` ／ `:259` `(cv1 - cv0) <= (v1 - v0) / 4.0`。live **是（二级兜底）**（`brep_exchange.rs:96`、`brep_tools.rs:74`）。OCCT **未找到**：`BRepMesh_TorusRangeSplitter.cxx:34-69`、`BRepMesh_SphereRangeSplitter.cxx:25-36`、`GCPnts_TangentialDeflection::ArcAngularStep`；子分深度是 Delaunay 自适应，无 `9`、无 `/96`、无 `/4`。影响**是**。建议：标「未移植（legacy 四叉树）」并删两个 `.or_else` 兜底。

**11｜`step/read_geometry.rs:348-370` `ProjectAct` 缺 Ellipse/Parabola/Hyperbola 精确臂｜未登记**
`:349-366` `if curve.is_line() { ... } else if let Some(c) = curve.gp_circ() { ... } else { project_on_segments(..., 25, ...) }`。live **是**。OCCT `ShapeAnalysis_Curve.cxx:382` `case GeomAbs_Parabola:` / `:383` `theProjParam = ElCLib::Parameter(theCurve.Parabola(), thePoint);` / `:394` `case GeomAbs_Ellipse:`。影响**是**。建议：补 `ElCLib` 精确臂（`ShapeAnalysis_Curve.cxx:376-400`）。

**12｜`brep_exchange.rs:118,125` + `occt-core/src/io/ply.rs:95-102` PLY 焊接 + 头/属性类型偏离 RWPly｜自创/未登记**
`brep_exchange.rs:118` `weld_vertices(&mut mesh, 1e-9);`；`ply.rs:99` `"property list uchar int vertex_indices\nend_header\n"`。live **是**。OCCT `RWPly_PlyWriterContext.cxx:276` `NCollection_Vec3<int>(myVertOffset) + theTri`（**每面重复自身 node+偏移，从不合并**）；`:214` `"property list uchar uint vertex_indices\n"`；`:156` comment 行。影响**是**，PLY 顶点数（盒体 24→8）与面索引全不同。建议：删 `weld_vertices`，改 per-face 顶点池。

**13｜`occt-core/src/io/stl.rs:92,112,214` STL 退化法向阈值 / 头 / 格式嗅探自创｜自创**
`:214` `if len < 1e-12 { [0.0, 0.0, 0.0] }`；另 `:112` `let name = b"solid generated by occt-core";`、`:92` `starts_with("solid") && contains("facet")`。live **是**。OCCT `RWStl.cxx:325` `if (aVNorm.SquareMagnitude() > gp::Resolution())`，`gp.hxx:60` `Resolution() { return RealSmall(); }` = `DBL_MIN`；`:374` 的 80 字节头；`RWStl_Reader.cxx:218` 按前 134 字节非 ASCII 字节嗅探。影响**是**，退化三角形法向被清零（<1e-12 vs OCCT 1e-308 量级）；二进制头字节必然不同。建议：用 `RealSmall()` 等价常量，移植 `RWStl.cxx:374` 头与 `RWStl_Reader.cxx:218` 嗅探。

**14｜`iges.rs:7-11,148-173,389-438` IGES 曲线族采样判定（头已登记，阈值仍自创）｜已登记 + 自创阈值**
`:168` `if c.is_periodic() || (max_d2 - min_d2) / max_d2 < 0.02 { CurveKind::Circle }`。live **否**（无门禁引用 `write_iges`）。OCCT `GeomToIGES_GeomCurve.cxx:504` `if (start->IsKind(STANDARD_TYPE(Geom_Circle)))`；完整球走 `GeomToIGES_GeomSurface.cxx:576-585` → `:1381` `IGESSolid_SphericalSurface`，而 `iges.rs:409-433` 自造 120 回转面 + 两条经线弧。影响：否（门禁外）；纳入门禁时 `<2%`、`:390` 降弦线、`:435-438` 降 108 平面会直接改实体集。建议：已登记充分；纳入门禁须按 `GeomToIGES_*` 重写。

**15｜`step/write_context.rs:16-17,43` + `mesh_to_brep.rs:102,109,184` 写侧/网格回灌的自创规则｜自创**（均非门禁路径）
`step/write_context.rs:43` `let n = 8;`、`:16-17` 曲面 6×6 采样重拟 vs OCCT `GeomToStep_MakeCurve.cxx:94-99` → `MakeBoundedCurve` 写**真实极点**（门禁不走写侧）。`mesh_to_brep.rs:102` `make_wire(&[e_ab, e_bc, e_ca])`；顶点按 1e-9 归并、Vertex 容差 0.0；`edge_use==2` → `Solid`；OCCT `BRepBuilderAPI_MakeShapeOnMesh`：`:109` node 1:1 + 容差 `BRepLib::Precision()`=1e-7、`:127/:138` 跳过退化三角形、`:184` 只返回 `Compound`。`mesh_to_brep` 影响**是（BOP 侧）**（`_board.md` T-01 已定位到此）；建议按 `MakeShapeOnMesh` 移植并成套移 `MakeWire::Add` + `MakeFace` 朝向规则。

**16｜`vrml.rs:92` + `gltf.rs:1-13,56` + `mesh.rs:50-107` 展示类导出（已登记）｜已登记**
`vrml.rs:92`（`solid TRUE`）；`gltf.rs:56` `let mesh = crate::shape_mesh::mesh_shape(shape, deflection);`。live **否**。OCCT `VrmlData_IndexedFaceSet.cxx:493`（`IsSolid=false`，见 `VrmlData_ShapeConvert.cxx:360`）；`RWGltf_CafWriter.cxx:258,356`（Y-up/单位）、`:378`（无 NORMAL 不写）、`:426`（TEXCOORD_0）；`mesh_sphere`/`mesh_cylinder` 连 deflection 参数都没有（OCCT `BRepMesh_IncrementalMesh.cxx:50`）；`compute_vertex_normals` 退化回零矢量 vs `Poly_Triangulation.cxx:468` 的 (0,0,1)。影响：否。建议 `vrml.rs:92` 改 `solid FALSE`。

**17｜已登记的 UNPORTED / PARK 清单（无需动作）｜已登记**
`brep_bnd_lib.rs:13-24`（三角化臂 PARK）、`:58`（pcurve-only 边 PARK）；`model_healer.rs:40-53`（`DEGENERATE_AREA_EPS` 自承非 OCCT 派生）；`delaun_index.rs:612-616`（DUPLICATE PORT 非 live）+`:627-631`（`UNPORTED_DETERMINANT_GUARD`）；`face_discret.rs:471-497` `is_degenerate`（非 live 分支）；`step/read_topology.rs:207,1028,1197`（`TranslateEdgeLoop`/`MakeTrimmedCurve(2d)` 未移植）、`:250-254`（ADVANCED_FACE 参数顺序自适应）；`step/read_geometry.rs:189`、`:856`。影响：`p04/p05` 未移植项是 `Shape`/`Shape-2` 残差（`step_obj_parity` 容差 0.1/0.05）的根因，与 `brep_exchange.rs:50-52` 自承的「SURFACE_CURVE pcurve 被丢」同源。建议：补 `StepToGeom::MakeTrimmedCurve`（`StepToGeom.cxx:2376-2394`）/`MakeTrimmedCurve2d`（`:2517-2561`）。

**18｜低危：参数/常量级偏离（详见文末「低危偏离」节）｜自创为主**
`param_resolution` 的 `h=1e-4` 差分（锥面静默错臂）、`angle_interior` 回退、`shape_size` 下限、Parabola/Hyperbola min-points、Bezier `AdjustRange` 永不触发、`impl DeflectionMesh for Delaun` 未 `Scale`、`deflection.rs` 采样偏转、`ply.rs` 读侧。门禁层面均未报。

## 附：抽样逐行核对确认忠实、不要动的部分

- `range_splitter/splitter.rs:534-627` Torus ↔ `BRepMesh_TorusRangeSplitter.cxx:21-116`（`oldDv`、`Du *= min(oldDv,Du)/aa`、`ru > 1e-16`、`aa < gp::Resolution()`、`nbU = max(nbU, nbV*diffU*R/(diffV*r)/5.)`、`fillParams` 的 `theScale` 0.5/2/3）；`:386-419` Sphere ↔ `BRepMesh_SphereRangeSplitter.cxx:21-56`（`0.7*`、`PCONFUSION` 上界、交错 `aHalfDu`、`computeStep`）；Cylinder 的 `nbV=0`（OCCT 注释块）、Cone、UVParam 亦对齐
- `range_splitter/param_set.rs` 的 `update_range`/`computeLengthU,V`/`computeTolerance`/`computeDelta`/`filterParameters` 忠实；`node_insertion.rs:881-945` ↔ `BRepMesh_GeomTool.cxx:465-512` + `ComputeErrFactors:32-84` + `AdjustCellsCounts:86-144`
- `delaun/constants.rs:150-183` `MakeCircle` ↔ `BRepMesh_CircleTool`；`delaun/{constants,triangulation,frontier}` 的 `createTriangles`/`createTrianglesOnNewVertices`（含 `cxx:703` 无条件 `ProcessConstraints`）/`isBoundToFrontier`/`cleanupMesh`/`frontierAdjust`/`meshLeftPolygonOf`/`findNextPolygonLink`/`checkIntersection`/`cleanupPolygon` 逐行对上
- `face_discret.rs:370-439` FaceChecker ↔ `BRepMesh_FaceChecker.cxx:150-186` + `Accept`（5° 切向阈值、`pi*(2*defl)^2` 小环判据）；`Classifier::RegisterWire/Perform` 忠实
- `edge_discret.rs` 的 `splitSegment`/`splitByDeflection2d`/`AddPoint`/`PerformCircular`（`min(ceil(diff/Du),1e6)`）/`ArcAngularStep`/`PConfusion()/10` 均对齐；`incremental_mesh/discret_root.rs:570` `AMP_ITERS = 5` ↔ `BRepMesh_ModelHealer.cxx:179` `aIterNb = 5`；`deflection_control/debug_flags.rs:203` `MAX_PASSES = 11` ↔ OCCT `aIterationsNb`
- `step/read_geometry.rs:476-504` `shift_displaced_line` ↔ `TranslateEdge.cxx:461-469`；`:605-700` `update_param3d` ↔ `StepToTopoDS_GeometricTool.cxx:238-406`；`shape_analysis_project` ↔ `ShapeAnalysis_Curve.cxx:147-201`（含 `:194` `delta = min(GAC.Resolution(preci), (uMax-uMin)*0.1)`）；`ProjectAct` closed/old-solution ↔ `:480-495`；`[40,20,25,40]` ↔ `:449`；`(w1, w1+2π)` ↔ `ElCLib.cxx:139-148`（可疑但仅 `p1≈p2` 全周期椭圆可达：`step/read_geometry.rs:584-607` 椭圆臂 `atan2(-v·ydir/b, v·xdir/a)` 与 `ElCLib.cxx:1234-1248` 符号相反）
- `brep_bnd_lib.rs:56-141` ↔ `BRepBndLib.cxx:81-215`（`useTriangulation=false`）；`prs3d_get_deflection`（`brep_exchange.rs:56-77`）与 `Prs3d.hxx:83-105` 逐行等同
- `occt-core/src/io/obj.rs:1-9` 已声明「不是 `RWObj_CafWriter` 的移植」：**读侧**可用；**写侧**缺 `vn`（`RWObj_CafWriter.cxx:317` `if (theWriter.HasNormals() && !writeNormals(...))`、`RWObj_ObjWriterContext.cxx:281`），`brep_exchange.rs:15,18` 恒写空 normals，`data/occ-cube.obj` 含 `vn`。消费者只读 `v `/`f `，现有断言不失败 —— 备案。

## 附：低危偏离（行号已核）

`range_splitter/param_set.rs:205-215` `param_resolution` 的 `let h = 1e-4;` 差分分支 vs `GeomAdaptor_Surface.cxx:1828-1888` 解析分支（Plane `R3d`、Cylinder `R3d/(2R)`→`2*asin`、Cone `R3d/R`、Sphere、Torus、Bezier/BSpline）——差分在柱/球/环 ≈`tol/R`，**但锥用 V 中值半径，OCCT 用 V 两端最大半径** ⇒ 静默错臂。`deflection_control/debug_flags.rs:595-600` `angle_interior` 回退成 `angle`，OCCT `BRepMesh_IncrementalMesh.hxx:104-107` 是 `2.0 * Angle`（live 经 `incremental_mesh/discret_root.rs:25-26` 已归一化）。`deflection.rs:69` + `edge_discret.rs:616` `shape_size = box_max_dimension.max(relative_deflection)` vs `BRepMesh_Deflection.cxx:42-43` + `BRepMesh_ShapeTool.cxx:93`（`BoxMaxDimension` **覆盖**初值；live 版本 `model_builder/preprocessor.rs:22` 写法正确）。`edge_discret.rs:949-957` Parabola/Hyperbola min-points 门槛退到 2，OCCT `BRepMesh_CurveTessellator.cxx:100-112` 对 Circle/Ellipse/Parabola/Hyperbola 都是 4（已登记）。`range_splitter/splitter.rs:656-665` OCCT 的 Bezier `AdjustRange` 有效性检查（`myIsValid=false`，`BRepMesh_NURBSRangeSplitter.cxx:298-311`）**永不触发**，因 `classify_surface`（`range_splitter/param_set.rs:86-124`）从不返回 `BezierSurface`。`deflection_control/debug_flags.rs:43-62,466-468` `impl DeflectionMesh for Delaun` 不做 `getRangeSplitter().Scale(...)`，把 face-basis 坐标当 UV 求值（OCCT `BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx:197,375,427`；live 走正确的 `node_insertion.rs:838-861` `ParametricDelaun` ⇒ 仅测试受影响）。`deflection.rs:102-190` `curve_deflection`/`surface_deflection`（按 `samples` 采样 + 重心网格）：OCCT 无对应函数，文件头**未登记**这两支。`occt-core/src/io/ply.rs` **读侧** `parse_ply`/`to_triangulation`：OCCT 8.0.0 中**无对应物（TKDEPLY 只有写出类）**——纯自创但不进门禁。
