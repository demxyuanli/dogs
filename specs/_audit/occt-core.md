# occt-core 移植忠实度审查（2026-09-20）

> 范围：`crates/occt-core/src/**`（196 个 `.rs` / 33,919 行）。只读审查，未改任何源码。
> 方法：① 关键词扫描（`not ported|UNPORTED|PARKED|simplif|fallback|heuristic|sampled|empiric|for now|workaround` → **40 命中**；`Source:` → **198 命中**；`ponytail` → **7 命中**）；② 对 12 个候选逐条读上下文并回查 `D:\source\OCCT-src`；③ 抽样逐行对读 7 个"声称移植"的函数（见 §4）。

## 结论摘要

1. **自创 15 项**（其中 5 项在 live 几何/网格路径上，会静默改变数值结果）：`gcpnts.rs` 的 `UniformDeflection`/`TangentialDeflection`、`elib/intersect.rs::circle_plane_intersection`、`elib/surface_eval.rs::surface_d2`、`geom/polyline_simplify.rs`（RDP）、`elib/measure.rs`（整套自造积分/采样）。
2. **最高危**：`crates/occt-core/src/gcpnts.rs:87-137` —— 用"递归中点二分 + `MAX_DEPTH = 16` 截断"冒充 `GCPnts_UniformDeflection`，缺 OCCT 的 `Linear/Circular/Curved/Composite` 分派与末尾 `Controle` 修正；该函数经 `occt-topo/src/wireframe.rs:48`（`edge_to_polyline`）进入 `brepmesh`/`shape_mesh`/`geometry_query` 活路径。
3. **同一 crate 内已有忠实件**：`gcpnts_perform.rs`（`GCPnts_TangentialDeflection::PerformCurve` 逐行移植，`1.5/0.75*dusave` 已核对 `TangentialDeflection.cxx:836-857`）、`gcpnts_estim.rs`、`intana2d/`、`intres2d/`、`bspl/`、`math_*`。`gcpnts.rs` 与它们并存 ⇒ 属"该用忠实件却另造一套"。
4. **4 项"假出处"**（文件头写 `Source: Xxx/`，但 OCCT 的该包根本没有这些函数）：`convert/`、`cslib/mod.rs`、`gprop/mod.rs`、`bnd/{intersect,obb_pca}.rs`。这比实现差异更危险：会让后续审查误判为已对齐。
5. **7 项已登记**（`io/obj.rs`、`io/stl.rs`、`io/ply.rs`、`kernel/osd.rs::sha1_hex`、`bvh/bvh_query/mesh_queries.rs` 全扫描回退、`geom/csg.rs`（消费侧登记）、`elib/clib.rs` 的两处例外返回）——见 §3，本身可接受但仍列出。
6. **未发现**规则明令禁止的"面积比/长度滤边/体积门"式自创；`bspl/`、`gp/`、`intana2d/`、`intres2d/`、`intf/`、`intcurve/`、`intimpargen/`、`math_*`、`gprop/gprops/`、`poly/{connect,merge_nodes,triangulation_full}` 抽样比对未发现自创。

## 发现

### 1. `UniformDeflection` 用二分+深度截断冒充 `GCPnts_UniformDeflection`（live 网格路径）
- 判定：**自创**
- 证据：`crates/occt-core/src/gcpnts.rs:30` `const MAX_DEPTH: usize = 16;`；`crates/occt-core/src/gcpnts.rs:124`
  `if dev > tol && depth < MAX_DEPTH && params.len() < max_pts { subdivide(...) }`
  `crates/occt-core/src/gcpnts.rs:122` `let dev = point_segment_dist(&pm, &pa, &pb);`（只看中点，不是真最大偏差）
- OCCT 对应：`GCPnts_UniformDeflection.cxx:311-352`（`initialize`：`GetDefType` 分派 + `Controle` 末尾点修正）、`:192-206`（Circular 解析点数）、`CPnts_UniformDeflection.cxx:190-249`（Curved 推进器）。**OCCT 没有深度截断**。
- 影响：**会改变结果**。直线/圆走 OCCT 的解析分支（直线上恰好 2 点）；此实现永远走二分，点数与参数分布都不同。live：`occt-topo/src/wireframe.rs:48` → `edge_to_polyline` → `brepmesh.rs`/`shape_mesh.rs`/`geometry_query.rs`；另 `occt-geom/src/curve_approx.rs:26`、`occt-geom/src/gcpnts.rs:288`。
- 建议：删掉 `subdivide`，改调 `GCPnts_UniformDeflection.cxx:311-352` 的四分支 + `Controle`（`:142-160`）；Curved 分支按 `CPnts_UniformDeflection.cxx` 移植。当前实现按 `:30` 标"未移植"。

### 2. `QuasiUniformDeflection`/`TangentialDeflection` 同为二分近似，且与同 crate 忠实件冲突
- 判定：**自创**
- 证据：`crates/occt-core/src/gcpnts.rs:139-140` `/// count is always odd (matches `GCPnts_QuasiUniformDeflection`)`；`:149-150`
  `let base = UniformDeflection::from_curve(c, a, b, tol, max_pts); ... 插入中点`
  `crates/occt-core/src/gcpnts.rs:214` `if (dev > tol || dangle > angle_tol) && depth < MAX_DEPTH {`
- OCCT 对应：`GCPnts_QuasiUniformDeflection.cxx`（独立迭代器，非"二分+插中点"）；`GCPnts_TangentialDeflection.cxx:522-916`（已由 `gcpnts_perform.rs:52` 逐行移植）。
- 影响：会改变点数与参数。live：`occt-topo/src/meshing/geom_tool.rs:421,443`（`GeomTool::discretize_curve`/`discretize_iso_curve`）、`occt-geom/src/gcpnts.rs:298`。
- 建议：`from_curve_with_deriv` 改为委托 `perform_tangential_curve`（`gcpnts_perform.rs:52`，即 `TangentialDeflection.cxx:522-916`）；`QuasiUniformDeflection` 按 `GCPnts_QuasiUniformDeflection.cxx` 重写，否则标"未移植"。

### 3. `Color::from_name` 取了 OCCT 的 sRGB 列并四舍五入到 0.01
- 判定：**自创 / 不等价**
- 证据：`crates/occt-core/src/quantity/mod.rs:33` `NameOfColor::Brown=>Self{r:0.65,g:0.16,b:0.16},`
- OCCT 对应：`Quantity_Color.hxx:49-51` `Quantity_Color(const Quantity_NameOfColor theName) : myRgb(valuesOf(theName, Quantity_TOC_RGB))`；`Quantity_ColorTable.pxx:185`
  `RawColor(BROWN, sRGB, 0xA52A2A, 0.647059, 0.164706, 0.164706, RGB, 0.376262, 0.023153, 0.023153)` —— 命名构造函数取的是 **RGB（线性）列 0.376262/0.023153/0.023153**，Rust 给的是 sRGB 列再截成两位小数。
- 影响：数值不同（BROWN 差 ~1.7 倍）。当前 `from_name` 只被自身测试 `quantity/mod.rs:108` 调用 ⇒ 目前仅潜在；一旦用于 XCAF/STEP 颜色写回即写出错误颜色。另 OCCT 有数百个名字，此处仅 23 个。
- 建议：按 `Quantity_ColorTable.pxx` 建全表并取 `RGB` 列，entry 对 `Quantity_Color.hxx:49`。

### 4. `circle_plane_intersection` 的"共面"分支凭空造点
- 判定：**自创**
- 证据：`crates/occt-core/src/elib/intersect.rs:58-60`
  `if d.abs() < 1e-12 { // Coplanar — full circle` → `return vec![center, clib::circle_value(c, FRAC_PI_2)];`
  另 `:57` `if d.abs() > c.radius + 1e-12`、`:63-66` 自造 `cos_phi` 公式。
- OCCT 对应：**未找到对应分支**。`ElCLib` 无求交函数；`GeomAPI_IntCS`/`IntAna` 的圆-平面走解析求交，共面情形不返回"圆心 + π/2 点"这两点。文件头 `:7` 只写 `Source: inspired by GeomAPI_IntCS`。
- 影响：**返回几何上错误的点集**（共面时 OCCT 语义是整圆/无孤立交点）。若被调用即产出错误交线。
- 建议：改为移植 `IntAna_Quadric`/`GeomAPI_IntCS` 的对应分支；在写出等价分支前标"未移植"。

### 5. `surface_d2` 对球/环面静默返回零二阶导
- 判定：**自创（静默默认值掩盖失败）**
- 证据：`crates/occt-core/src/elib/surface_eval.rs:82` `/// ... — simplified for planes/cylinders.`；`:100` `SurfaceRef::Sphere(s) => { ...; (p, du, dv, GpVec::zero(), GpVec::zero(), GpVec::zero()) }`（Torus 同 `:104`）
- OCCT 对应：`ElSLib::SphereD2` / `TorusD2` 返回非零 `duu/dvv/duv`；本文件头 `:2` 却写 `Source: ElSLib.cxx D1/D2/Norm functions`。
- 影响：静默给出错误曲率/法向变化；调用方无法察觉（无 `Option`、无 panic）。
- 建议：按 `ElSLib.cxx` 的 `SphereD2`/`TorusD2` 补齐；未补齐前改为返回 `Option`/`Result` 或标"未移植"。

### 6. `Bnd_Sphere::Add` 缺"被包含"分支，`Distance` 语义被改
- 判定：**自创 / 不等价**
- 证据：`crates/occt-core/src/bnd/bsphere.rs:23-25`
  `let new_r = (d + self.radius + other.radius) * 0.5; if new_r <= self.radius { return; } ...`
  `crates/occt-core/src/bnd/bsphere.rs:38` `let d = p.coord.subtracted(&self.center.coord).modulus() - self.radius;`（返回**到球面**距离）
- OCCT 对应：`Bnd_Sphere.cxx:73-96`：`if (myRadius + aDist <= theOther.myRadius) { *this = theOther; return; }`（被包含 → 整体替换，圆心也要换）；`Bnd_Sphere.cxx:64` `Distance(theNode)` 返回**到球心**的 `Modulus()`（到球面的最小距离在 `:56-62` 的 `Distances`）。
- 影响：被包含时圆心不更新 ⇒ 包围球错误（偏心/偏小），会漏检；`distance` 与 OCCT 同名函数返回不同量。当前 `BndSphere` 无仓内调用者 ⇒ 潜在。
- 建议：按 `Bnd_Sphere.cxx:73-96` 补三个分支；`distance` 改名或按 `:64` 实现。

### 7. RDP 抽稀冒充 `ShapeAnalysis_FreeBoundData`
- 判定：**自创（假出处）**
- 证据：`crates/occt-core/src/geom/polyline_simplify.rs:2` `//! Source: ShapeAnalysis_FreeBoundData, math decimation (RDP).`；`:8` `pub fn rdp_simplify(points: &[GpPnt], tol: f64) -> Vec<usize>`
- OCCT 对应：**未找到对应分支**。`ShapeAnalysis_FreeBoundData` 只有自由边界的数据容器，无 Douglas–Peucker 抽稀；OCCT 无 RDP 实现。
- 影响：**会改变交线几何**。live：`occt-topo/src/intpatch_trace.rs:260-262` 对 `poly.len() > 8` 的交线按 `tol` 抽稀后才建 `IntersectionCurve` —— 抽稀阈值 OCCT 中不存在，抽稀后顶点即交线顶点。
- 建议：`intpatch_trace` 侧改按 `IntPatch` 真实交点序列建线；`rdp_simplify` 若仅诊断用则标"未移植 + 非 OCCT"。

### 8. `elib/measure.rs` 整套自造积分/采样，无 `Source:`、无未移植标记
- 判定：**自创（未登记）**
- 证据：`crates/occt-core/src/elib/measure.rs:52` `pub fn simpson<F: Fn(f64) -> f64>(...)`；`:64-76` `curve_arc_length`（累加弦长）；`:111` `point_at_arc_length`（"Walk samples"）；`:133` `surface_patch_area`（中点矩形和）
- OCCT 对应：**未找到对应分支**。`ElCLib`/`ElSLib` 无弧长/面积函数；等价物是 `GCPnts_AbscissaPoint`（对弧长积分做 Newton 迭代，带容差）与 `GProp_*`。
- 影响：**会改变数值**。`curve_arc_length` 被 `occt-topo/src/inttools_range.rs:28,175` 调用（该处 `:13,26` 已登记为近似）。同 crate 另有一份重复实现：`geom/curve_frenet.rs:156`（`curve_arc_length`）与 `:116` `arc_length_parameters`（弦长参数化，声称用于 `GCPnts_UniformAbscissa`，而 OCCT 该函数基于精确弧长）。
- 建议：统一换 `GCPnts_AbscissaPoint`（`CPnts_AbscissaPoint.cxx` + `CPnts_MyRootFunction.cxx`）；`measure.rs` 保留部分在文件头写 `UNPORTED`。

### 9. 体素 CSG 被当作布尔实现
- 判定：**自创**（理由写在注释，但无 OCCT 出处、无 UNPORTED）
- 证据：`crates/occt-core/src/geom/csg.rs:2` `//! Simplified voxel-approximation CSG — robust for prototyping.`；`:63-66` `volume() = count_inside * voxel_volume`
- OCCT 对应：**未找到对应分支**。OCCT 无体素布尔（`BRepAlgoAPI_*` + `BOPAlgo_Builder`）。
- 影响：**几何结果随分辨率变化**，非精确布尔。live：`occt-topo/src/solid_union.rs:13,207-229`（消费侧 `:7-10` 已登记为 "deliberate approximations"）。
- 建议：occt-core 侧文件头补 `UNPORTED`；布尔侧指向 `BOPAlgo_Builder` 缺口（与 `_index.md` §1-A5 同族）。

### 10. `Poly_MakeLoops::chooseLeftWay` 丢掉了 OCCT 的最小夹角选择
- 判定：**自创（缺分支）**
- 证据：`crates/occt-core/src/poly/make_loops.rs:234-237`
  `pub fn choose_left_way(&self, _node: usize, _seg_index: isize, lst_ind_s: &[isize]) -> isize { /* ponytail: first candidate */ lst_ind_s[0] }`
- OCCT 对应：`Poly_MakeLoops.hxx:223` 基类为纯虚 `= 0`；两个实现都在做夹角选择：`Poly_MakeLoops.cxx:611-676`（3D，取 `aAngleMin`）、`:688-700`（2D，另受 `myRightWay` 控制），仅在 `GetNormal`/`GetLastTangent` 失败时 `return theLstIndS.First()`。Rust 注释"geometric subclasses would pick..."把 OCCT 的**自带分支**说成了子类扩展。
- 影响：环路定向错误 ⇒ 三角化/回路朝向不同。当前 `PolyMakeLoops` 仅被自身测试调用（`make_loops.rs:448,471,489,519`）⇒ 目前死路径。
- 建议：补 `Poly_MakeLoops.cxx:611-676` / `:688-700` 两个实现与 `myRightWay`；在此之前标注"未移植"。

### 11. `int/curve_curve.rs` 是 Ericson 线段算法，却标为 `IntCurveCurve`/`IntTools` 移植
- 判定：**自创（假出处）**
- 证据：`crates/occt-core/src/int/curve_curve.rs:2` `//! Source: IntCurveCurve_IntImpCurveCurve, IntTools.`；`:94` `/// Closest points on two 3D segments (Ericson, Real-Time Collision Detection).`；`:29` 自造阈值 `cross.square_magnitude() <= tol * tol * l1 * l1 * l2 * l2`
- OCCT 对应：**未找到对应分支**。`IntCurve_IntImpConicParConic` / `IntCurveCurve_IntImpCurveCurve` 建在参数曲线 + `math_FunctionAllRoots`/`Extrema` 上（见 `intcurve/mod.rs:20-73` 的未移植清单），不做"线段夹紧最近点"。
- 影响：与 OCCT 交点/容差判定不同（`EPS=1e-14`、`tol²l1²l2²` 门限 OCCT 中不存在）。
- 建议：文件头改为"非 OCCT；线段近似"，真实交线走 `intana2d/` 与 `intimpargen/intersector.rs`（已移植）。

### 12. `convert/` 包假借 OCCT `Convert/` 之名，实为 GIS/坐标换算自造件
- 判定：**自创（假出处）**
- 证据：`crates/occt-core/src/convert/mod.rs:1` `//! Coordinate system conversions. Source: `Convert/``；`crates/occt-core/src/convert/mapproj.rs:1`（Mercator/等距/Lambert/UTM）；`crates/occt-core/src/convert/advanced.rs:52`（Rodrigues 旋转）
- OCCT 对应：`src/FoundationClasses/TKMath/Convert/` 只有 28 个 B 样条转换类（`Convert_CircleToBSplineCurve`、`Convert_CompPolynomialToPoles`…），**没有**极坐标/柱坐标/球坐标/地图投影的任何一个。
- 影响：非 OCCT 语义；`mapproj.rs`（`:1-59`）整块属 GIS 领域，OCCT 无对应物。仓内无生产调用者。
- 建议：拆包：B 样条转换按 `Convert/` 移植；其余移到非 `Source: Convert/` 的辅助模块并标注"非 OCCT"。

### 13. `cslib/mod.rs` 的 `classify_point`/`sphere_normal` 在 OCCT `CSLib` 中不存在
- 判定：**自创（假出处）**
- 证据：`crates/occt-core/src/cslib/mod.rs:1` `//! Classifier for parametric surfaces. Source: `CSLib/``；`:24-33` `classify_point(...)` 用"平面投影符号"分类；`:36-41` `sphere_normal`
- OCCT 对应：OCCT `CSLib` 只有 `CSLib::Normal` / `DNNormal` / `Class2d` / `NormalPolyDef`（本 crate 的 `cslib/{normal,dn_normal,class2d,poly_def}.rs` 才是移植件）。**没有** `classify_point`。
- 影响：分类结果与 `BRepClass3d_SolidClassifier` 不同（后者是本仓已有忠实件）。仓内无外部调用者。
- 建议：删除或标注"非 OCCT 辅助函数"，文件头去掉 `Source: CSLib/`。

### 14. `gprop/mod.rs` 的网格 `GProperties` 与真实移植件 `gprop/gprops/` 并存
- 判定：**自创（假出处）**
- 证据：`crates/occt-core/src/gprop/mod.rs:1` `//! ... Source: `GProp/``；`:27-42` `add_triangle(...)` 用 `idx(a)+idx(b)+idx(c)` 的 4 顶点对称和算惯量（`:37`）；`:71-78` `mesh_volume` 四面体和
- OCCT 对应：**未找到对应分支**。`GProp_GProps`/`GProp_SelGProps`/`GProp_VelGProps` 是解析积分（本 crate 已移植于 `gprop/gprops/props.rs:479,539,626`，且那些函数头**明确**写 "approximated by planar triangular facets"）。
- 影响：质量/惯量随三角化变化。live 仅 `centroid_of_points`（`occt-topo/src/brep_gprop.rs:17`）。
- 建议：`gprop/mod.rs` 改为薄层委托 `gprop/gprops/`；`GProperties` 若保留则标注"三角面近似，非 `GProp_GProps`"。

### 15. `bnd/`：两处 `Source:` 指向并不含该功能的 OCCT 头
- 判定：**自创（假出处）**
- 证据：`crates/occt-core/src/bnd/obb_pca.rs:1` `//! PCA-based OBB computation. Source: Bnd_OBB + Jacobi eigenvalue decomposition.`（`:7` `compute_obb_pca`）；`crates/occt-core/src/bnd/intersect.rs:1` `//! ... Source: Bnd_Tools`（`:6` `ray_box_intersect`、`:25` `box_box_distance`）
- OCCT 对应：`Bnd_OBB.hxx:54,61,85` 只有 3 个构造函数（空、显式参数、`Bnd_Box`），**无 PCA/点云拟合**；`Bnd_Tools.hxx:24-42` 只有 `Bnd2BVH` 两个转换函数，**无** ray-box / box-box。射线求交在 OCCT 属 `BVH_Tools::RayBoxIntersection`。
- 影响：无（两文件均无生产调用者，`compute_obb_pca` 只有自身测试）。风险在于出处误导后续审查。
- 建议：改标真实出处（`Bnd_OBB.hxx` 的 3 个 ctor 用 `Bnd_OBB`，PCA 部分写"非 OCCT"；ray-box 指向 `BVH_Tools.hxx`）。

## 3. 已登记的等价替换（可接受，仍列出）

| 位置 | 登记内容 | 备注 |
|---|---|---|
| `io/obj.rs:1-9` | "not a port of an OCCT class" + 列出 `RWObj_*` 真实类名 | 登记完整 |
| `io/stl.rs:1` / `io/ply.rs:1` | `Source: RWStl/RWPLY (simplified)` | 有理由、无行号 |
| `kernel/osd.rs:114-126` | `sha1_hex` = FNV-1a，`ponytail: placeholder ... NOT cryptographic` | **不等价**（非 SHA-1），已登记故列此 |
| `bvh/bvh_query/mesh_queries.rs:365-408` | BVH 探测→全扫描回退，`ponytail:` 说明理由 | 结果仍正确，仅性能 |
| `poly/merge_nodes.rs:176` / `gp/ax1.rs:55` / `gp/quaternion.rs:140` / `gp/trsf.rs:298` | `ponytail:` 实现选择说明 | 不改变结果 |
| `elib/clib.rs:177-196` | `circle2d_parameter` 退化时返回 `0.0` 而非 OCCT 的 `Standard_ConstructionError` | 出口语义不同，已写明 |
| `elib/clib.rs:46-55` | `ellipse_value` 符号修正 + ATU01038 实证 | 修正**向 OCCT 对齐**，非自创 |
| `geom/curve_interp3d.rs:1-6` | Catmull–Rom 明示为"plus ... helpers"、非 `GeomAPI` | 非 OCCT 但已声明；仅测试调用 |

## 4. 抽样逐行对读结果（3–5 处要求的扩展版，共 7 处）

| Rust | OCCT | 结论 |
|---|---|---|
| `elib/clib.rs:213-229` `circle_parameter` | `ElCLib.cxx:1199-1222` | **忠实**（含两个 `gp::Resolution()` 门与 `normalizeAngle`） |
| `elib/clib.rs:200-210` `normalize_angle` | `ElCLib.cxx:56-72` | **忠实**（双 `while` + `RESOLUTION` 收边） |
| `math_direct_poly_roots.rs:47-75` `refine_root` | `math_DirectPolynomialRoots.cxx:92-120` | **忠实**（迭代界、双重 break、更优解回退一致） |
| `gcpnts_perform.rs:317-333` `1.5/0.75*dusave` | `GCPnts_TangentialDeflection.cxx:836-857` | **忠实**（解决 `_index.md:109` 的待核项） |
| `intana2d/ana_intersection.rs:27-80` 状态器 | `IntAna2d_AnaIntersection.lxx:17-83` | **忠实**（`done` 断言、1-based `Point`） |
| `bnd/bsphere.rs:19-29` `add_sphere` | `Bnd_Sphere.cxx:73-96` | **不等价** → 发现 6 |
| `poly/make_loops.rs:234` `choose_left_way` | `Poly_MakeLoops.cxx:611-676,688-700` | **缺分支** → 发现 10 |

## 5. 未发现自创的区域（附证据）

- `bspl/`（30 文件）、`gp/`（26 文件）、`math_*`（8 文件）、`intana2d/`、`intres2d/`、`intf/`、`intcurve/`、`intimpargen/`、`elib/{clib,clib2d,slib}.rs`：198 条 `Source:` 中绝大多数集中于此，且抽查的 4 个函数逐行一致；`intcurve/mod.rs:20-73` 与 `intf/mod.rs:6-14` 的未移植清单**逐文件给了字节数**，属高质量登记。
- `cslib/{normal,dn_normal,class2d,poly_def}.rs`：`CSLib.cxx:84-151/183-387/391-568` 行号可核。
- `gprop/gprops/`：解析件与近似件在注释里分开标注（`gprops/props.rs:479,539,626`）。
- `poly/{connect,merge_nodes,triangulation_full}.rs`：控制流有 `Poly_Connect.cxx`/`Poly_MergeNodesTool.cxx` 出处；`merge_nodes.rs:176` 的 `ponytail` 只是预分配提示。
- 规则禁令专项（面积比/长度滤边/体积门）：本轮 40 条关键词命中中未出现用于生产判定的此类谓词（与 `_index.md:108` 一致）。
- **未逐项核对**（建议下轮）：`bnd/sortbox.rs`（`Bnd_BoundSortBox` 的网格划分方式与本实现 `new(elements,nx,ny,nz)` 不同，需回读 `.cxx`）、`geom/` 其余文件（`delaunay.rs`/`triangulate.rs`/`fit*.rs`/`polygon_*.rs`/`mesh_analysis.rs`/`curve_ops3d.rs`：均为无行号的自造几何辅助，属"整个 `geom/` 包不是移植件"这一系统性事实，非单点缺陷）。
