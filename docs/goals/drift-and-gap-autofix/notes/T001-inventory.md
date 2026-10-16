# T001 receipt — live drift/gap inventory (2026-10-15)

Board: `docs/goals/drift-and-gap-autofix`.
Method: four read-only Scout passes over `~/.cursor/projects/d-source-repos-dogs/canvases/occt-code-review-status.canvas.tsx` (SEMANTIC array),
`specs/_audit/{_index,topo-rest,occt-geom,occt-core,bop-boolean,mesh-exchange}.md`,
`specs/_tasks.md`, `specs/_a3n00_gap_analysis.md`, plus a full sweep of every `UNPORTED` /
`not ported` marker in the five crates. Each line re-verified against `D:/source/OCCT-src/src`.
No file was edited; no gate was run for this task.

## A. Registered findings that are already CLOSED (evidence recorded, docs still say live)

| id | area | current evidence | note |
|---|---|---|---|
| topo-rest 1 | brepfeat | `features.rs:93,290,452` -> `bop_builder::boolean`, volume via `brep_gprop_full::volume_properties` | audit's own 2026-10-10 correction confirmed |
| topo-rest 4 | brepfeat | `resolution_for` / `clamp(16,64)` no longer exist anywhere in `crates/` | removed with the analytic-tool cleanup |
| topo-rest 6 | wireframe | `wireframe.rs:110` -> `IncrementalMesh::discretize_face`; `face_uv_bounds` -> `brep_tools::uv_bounds` | untrimmed natural-UV grid gone |
| topo-rest 11 | brep_gprop | `brep_gprop.rs:40` -> `brep_gprop_full::{volume,surface}_properties`; mesh path only as the `NoSurf` arm `:56` | matches `BRepGProp.cxx:344-360` |
| topo-rest 18 (classifier half) | brep_class3d | `other_segment.rs:100` matches `BRepClass3d_SolidExplorer.cxx:74-186`; `aTestInvert` at `:398-410` | bbox half still live, see B |
| bop 1,2,3,4,7,16,17 | BOP | symbols/files absent; `bop_curved/` deleted 2026-10-14; `bop_builder_dispatch.rs:233` -> `boolean_via_bopalgo` | residual self-invented voxel engine left in `boolean_ops.rs:173` (no caller) |
| bop 8,9,10,6 | BOP | `int_tools_full/context.rs:195` `Extrema_ExtPS`; `bop_draft_solid_occt.rs:50`; `edge_edge/edge_edge.rs:150-175` | new gaps spawned, see C |
| bop 18 | BOP | `bop_build_faces/builder_like.rs:208` test-only; live path `builder_face_occt.rs:318` uses `FClass2d::is_hole` | delete tests or redirect |
| mesh 1-14,19,19b | meshing/exchange | each verified against its cxx arm | see per-line evidence in the Scout output |
| occt-geom 1,2,3,4,6,11,12,13,14,15 | geom/geom2d | offset-curve utils, `Extrema_GGExtPC`, `GeomConvert_CompCurveToBSplineCurve`, `ElCLib` D2, Gauss `gcpnts`, `BSplCLib_2` | |
| occt-core 1,2,4,5,6,8,9,10,11,12,13,14,15 | core | per-line evidence recorded | |

## B. LIVE drift (divergence exists now; OCCT counterpart branch available)

Ranked by (pipeline centrality x blast radius x fix confidence).

### B1. `BOPAlgo_WireSplitter` family in `crates/occt-topo/src/wire_splitter_block.rs` — ✅ CLOSED 2026-10-15 (T-131), heading was stale
> Re-verified 2026-10-15 against the current tree: all three arms are already faithful.
> `uv_tolerance_2d` (`wire_splitter_block.rs:567-575`) forwards to
> `occt_geom::approx_same_parameter::{u_resolution,v_resolution}` and `fd_uv_tolerance` is gone
> (`_1.cxx:885-901`); `curve2d_radius` is gone and `angle_2d` (`:487-535`) uses
> `GeomLPropCLProps2d` (`_1.cxx:797-807`); `refine_angle_2d` (`:660-676`) takes `aTV` from
> `vertex_parameter_on_face` (`_1.cxx:1060`). Closed in commit `67d7f01f` (T-131 / 发现 19);
> gates were green then. The description below is kept as the historical record only.
- **B1a** `uv_tolerance_2d` (`:580-650`) re-implements `GeomAdaptor_Surface::U/VResolution`
  (`GeomAdaptor_Surface.cxx:1818-1958`) by hand and falls back to `fd_uv_tolerance` (`:657`,
  finite-difference `t3/|dS/du|`) for `SurfaceOfRevolution` (V) and `SurfaceOfExtrusion` (U).
  It also misses the `Geom_RectangularTrimmedSurface` unwrap (`cxx:423-430`) and the
  `GeomAbs_OffsetSurface` arm (`cxx:1882-1885`, `:1944-1947`), lacks the cone
  `myVLast - myVFirst > 1.e10` unbounded guard (`cxx:1855-1858`), and evaluates the cone radius by
  sampling `d0(0,v)`/`d0(PI,v)` instead of `S->VIso(V)->Radius()` (`cxx:1861-1866`).
  The shared faithful implementation already exists: `occt_geom::approx_same_parameter::{u_resolution,v_resolution}`
  (`approx_same_parameter.rs:386-494`, exact `GeomAdaptor_Surface::U/VResolution`).
- **B1b** `curve2d_radius` (`:543-554`) inlines `|D1 x D2| / |D1|^3` with guards that are
  dimensionally wrong (`n2 <= PCONFUSION` compares `|D1|^2` against `1e-9`, i.e. `|D1| <= 3.2e-5`,
  where OCCT compares against `PConfusion^2`), instead of `GeomLProp_CLProps2d`
  (`BOPAlgo_WireSplitter_1.cxx:797-807`). The faithful port is
  `occt_geom2d::cl_props2d::GeomLPropCLProps2d` (`is_tangent_defined`/`curvature`), already used at
  `brep_class/intersector.rs:176`.
- **B1c** `refine_angle_2d` (`:771-...`) takes `a_tv` from the 3D edge parameter
  (`vertex_parameter`, `:453`); OCCT takes it from the **pcurve**
  (`BRep_Tool::Parameter(aV, aE, myFace)`, `BOPAlgo_WireSplitter_1.cxx:1060`). The port already has
  the 3-arg overload as `vertex_parameter_on_face` (`:407-441`), used by `angle_2d` (`:489`).
  Divergent when the pcurve and 3D curve do not share a parameterisation.
- Reachability: `split_block` -> `wire_splitter.rs:329/349/365/486` (`WireSplitter`, live BOP).
- cxx: `BOPAlgo_WireSplitter_1.cxx:768-840` (Angle2D), `:859-901` (Tolerance2D/UTolerance2D/VTolerance2D), `:1060` (RefineAngle2D).

### B2. `pc_curve_kind` geometric-invariant classifier (`crates/occt-topo/src/pcurve.rs:218-237`) — ✅ CLOSED 2026-10-15 (T-132)
- Rust infers the kind from `is_periodic() && |period-2*PI| < 1e-9` and finite-range tests, with only
  4 kinds; OCCT downcasts in `Geom2dAdaptor_Curve::load` (`Geom2dAdaptor_Curve.cxx:279-345`) and reads
  `GetType()` (`:209`). Observable: an unbounded parabola/hyperbola pcurve falls into the `Line` arm
  where OCCT lands in the reject arm. Callers: `geom_bnd_lib_offset2d.rs:35,61`,
  `geom_bnd_lib_curve2d.rs:48`, `int_face_face_bounds.rs:181`.
- The trait already exposes every query needed (`is_geom2d_line`, `gp_circ2d`, `gp_elips2d`,
  `gp_parab2d`, `gp_hypr2d`, `is_bezier2d`, `is_bspline2d`, `offset_basis`).
- **Closed (T-132)**: `CurveKind`/`pc_curve_kind` deleted. The single faithful `GetType()` port now
  lives at `occt_geom2d::adaptor_curve::get_type2d` and is shared with `Extrema_Curve2dTool`
  (`extrema2d/curve2d_tool.rs`). `geom_bnd_lib_curve2d::box_curve2d` dispatches all nine
  `GeomBndLib_Curve2d.cxx:97-171` arms (Ellipse/Hyperbola/Parabola/Bezier/Offset were unreachable
  before), `geom_bnd_lib_offset2d` uses the exact `occ::down_cast<Geom2d_Line>/<Geom2d_Circle>`
  equivalent (`GeomBndLib_OffsetCurve2d.cxx:36,63`), and `int_face_face_bounds.rs:181` uses the exact
  `Geom2d_Line` downcast of `IntTools_FaceFace.cxx:2149`. Gates: `--lib` 1251/0, `step_obj_gates`
  4/4, `export_data_obj` 15/15 bit-identical. See `specs/_tasks.md` T-132.

### B3. `bbox_from_geometry.rs:46,56-67` — bbox by sampling (not in the audit) — ✅ CLOSED 2026-10-15 (T-134)
- 64-point edge bbox and 16x16 UV surface bbox instead of `BRepBndLib`/`BndLib` analytic extrema.
- Broadly used: `bop_builder.rs:86`, `bop_builder_dispatch.rs:109`, `bop_classify_occt.rs:78`,
  `builder_solid.rs:128`.
- **Closed (T-134)**: `shape_bbox` now delegates to `brep_bnd_lib::shape_bnd_box` =
  `BRepBndLib::Add(S, B, /*useTriangulation*/ false)` (`BRepBndLib.cxx:81-215`; faces via
  `BndLib_AddSurface::Add(BS, Tolerance(F), B)` `:109,118`, plane faces additionally add their
  boundary edges `:128` plus `B.Enlarge(Tolerance(F))` `:131`; edges/vertices outside faces/edges
  `:146-149`/`:212`). The sampled `add_to_bbox` walk and its single `face_uv_bounds` dependency are
  deleted; the three entry points keep their signatures. `Bnd_Box::Get()` now carries the
  `Precision::Confusion()` gap on both sides (`BRep_Tool::Tolerance` floor, `BRep_Tool.cxx:881-895`),
  so a unit box reads `[-1e-7, 1+1e-7]` as in OCCT; five existing lib assertions were corrected to
  that cxx box (no new tests, no gate tolerance relaxed). Gates: `--lib` 1251/0, `step_obj_gates`
  4/4 (23 models, area 1.0000), `export_data_obj` 15/15 bit-identical. Registered in
  `specs/_audit/topo-rest.md` §19. Side finding: `solid_union::shape_bbox` keeps its own 9x9
  sampling fallback (topo-rest §18), a separate path.

### B4. BOP classification / split residuals (live, registered) — ✅ CLOSED 2026-10-15 (T-135)
- `builder_solid.rs:139-152` `face_sample_point` uses the surface centre; OCCT first tries an edge not
  belonging to the solid then `PointInFace`/`PointNearEdge`
  (`BOPTools_AlgoTools.cxx:660-717`). A faithful copy already exists at
  `bop_build_common/draft_solid.rs:31-56`; the already-ported
  `algo_tools_face::compute_state_face_in_solid` (`algo_tools_face.rs:544`) is unwired.
  Duplicated at `bop_build_solids.rs:497`, `bop_build_common/build_ops.rs:95`, `bopalgo_tools_class.rs:114`.
  **Closed (T-135)**: the single faithful
  `algo_tools_face::{compute_state_face_point,compute_state_face}` (`cxx:666-715`, edge overload
  `cxx:734-788`, `PointInFace` `cxx:692`, `PointNearEdge` retry `cxx:694-707`) is now the only point
  selection: `bopalgo_tools_class` Face arm, `draft_solid::compute_state_by_one_point` (now a forwarder,
  `point_solid_state` deleted), `builder_solid::shape_inside_solid` = the static `IsInside`
  (`BOPAlgo_BuilderSolid.cxx:835-861`, fixed `Precision::Confusion()`, first face, infinite-point arm for
  faceless shapes), and `bop_build_solids::classify_solid_state` (`cxx:790-822`). `face_sample_point` and
  `classify_point` are deleted from `builder_solid.rs`, `build_ops.rs`, `bopalgo_tools_class.rs` and
  `draft_solid.rs`; `bop_build_solids::face_sample_point` survives only for the port-only
  `orient_face_outward` heuristic. Same family: the local `is_split_to_reverse` copies
  (`bop_build_solids.rs:454`, `build_ops.rs:82/112/149` with the `1e-6` finite-difference normal and the
  32-point `closest_point_on_edge`) are replaced by the single faithful
  `algo_tools_face::is_split_to_reverse` (`cxx:1255-1291` / `:1316-1427` / `:1432-1500`), whose face
  fast path was corrected from `shape_key` equality to `Arc::ptr_eq` (`cxx:1330-1333` Geom_Surface handle
  equality). Gates: `--lib` 1251/0, `step_obj_gates` 4/4, `export_data_obj` 15/15 bit-identical.
  Side finding (not fixed): OCCT `ComputeStateByOnePoint`'s `default` arm recurses into the first child
  and discards the result (`cxx:648-655`) so containers stay `TopAbs_UNKNOWN`; both Rust copies loop over
  children and return the first non-UNKNOWN state.
- `bop_builder.rs:88-89` bbox-disjoint short-circuit: OCCT `BOPAlgo_BOP` has none.
  **✅ CLOSED 2026-10-16 (P038 / T-168, `ff235d37`)**: the early return
  (`bbox_a.is_out_box(&bbox_b)` -> `disjoint_result(a, b, op)`) is deleted, with it the
  `disjoint_result` re-export; `boolean` now always runs `boolean_via_bopalgo`.
  `BOPAlgo_BOP::CheckData` (`BOPAlgo_BOP.cxx:106-160`) only demands two arguments and
  `Perform` / `PerformInternal1` (`cxx:364`, `:422`) run the pave filler unconditionally.
  Finding 4/7 are moot: they lived only in `crates/occt-topo/src/bop_curved/`, deleted
  2026-10-14.
- `bop_build_solids.rs:232/328` `close_open_shells` re-attachment + `geometrically_open` `continue`
  discard: OCCT attaches nothing after taking closed shells (`BOPAlgo_BuilderSolid.cxx:223`).
- `bop_geomlib_closed.rs:168,185` 23-sample iso test substituting for
  `GeomLib::IsBSplUClosed`/`IsBzUClosed` pole comparison (`GeomLib.cxx:2693`).
- `bop_split_seam.rs:7-8` AABB/sampling stand-in for `Geom2dAPI_ProjectPointOnCurve`
  (`BOPTools_AlgoTools3D.cxx:275`).

### B5. `inttools_range.rs:295-351` `find_nearest_valid_point` missing the local-singularity arm — ✅ CLOSED 2026-10-15 (T-133)
- Missing `BRepLib_1.cxx:66-78` (`aD1Mag = (1/Resolution(1.)*0.01)^2` for Bezier/BSpline after
  unwrapping `OffsetCurve`) and `:108-141` (the `D1`-accelerated doubling loop).
- **Closed (T-133)**: `aD1Mag` is now computed per `BRepLib_1.cxx:69-81` — `GeomAbs_OffsetCurve`
  first unwraps `OffsetCurve()->BasisCurve()` and re-reads the type, then only Bezier / B-spline
  take `(1./curve.resolution(1.0)*0.01)^2`; every other type keeps `0.` and skips the arm. The
  type test reuses the single faithful `Adaptor3d_Curve::GetType()` port at
  `edge_edge::curve_kind` (`GeomAdaptor_Curve.cxx:252-311`), re-exported `pub(crate)` from
  `edge_edge/mod.rs` (visibility only, no second implementation). The `cxx:108-141` loop is
  ported verbatim: with `aD1Mag > 0.` the walk takes `D1` at `u`, and while the point is inside
  the sphere and `|D1|^2 < aD1Mag` it doubles `aStepLocal`, continuing while still in range and
  falling back to `anEndU` (returning `None` when that end is inside too) otherwise. The
  `ponytail:` note that recorded the omission is deleted.
- Gates: `--all-targets` 0 errors (373 warnings), `--lib` 1251 / 0, `step_obj_gates` 4 / 4,
  `export_data_obj` 15 / 15 with v/f and bbox bit-identical to the T002 baseline.
- Disclosed scope deviation: `edge_edge/mod.rs` was outside the registered `allowed_files`;
  the change is a visibility-only re-export (0 new symbols, 0 behavior change).
- Reachable: `find_valid_range` <- `pave_ff_make.rs:409`, `pave_ff_pave_bound.rs:179`;
  `ShrunkRange` <- `pave_shrunk.rs:87`, `pave_ff_se.rs:39`, `pave_common/mod.rs:51`.

### B6. `occt-geom` Extrema failure fabrication — ✅ ALL SUB-ITEMS CLOSED 2026-10-15 (B6a–B6g through P006/P007/P008/P017/P018/P019/P020/P021); only B6h stays open
- **B6a ✅ CLOSED 2026-10-15 (T-136, package P006).** `extrema_ss.rs:449,506` and
  `extrema_surf/numeric_extrema.rs:153,374,433,514` fabricated a `fallback_*` pair where OCCT
  returns `myDone = false`. All three fallbacks are deleted and the three minimum-distance
  entries return `Option<ExtremaPair>` (`None` = `!IsDone()` / `NbExt() == 0`); the live consumer
  `bean_face_exact.rs:259` now takes the real `IntTools_BeanFaceIntersector.cxx:956` guard
  (not done => block skipped, no range flag) instead of the fabricated pair. cxx:
  `Extrema_ExtPS.cxx:307-371`, `Extrema_ExtCS.cxx:122,532`, `Extrema_ExtSS.cxx:226-260`,
  `Extrema_ExtElSS.cxx:62-83`, `Extrema_GGExtPC.hxx:531`. The prior registration claiming the
  fallbacks had to stay for their consumer and were gated on T-67 was false.
- **B6b ✅ CLOSED 2026-10-15 (T-137, package P007) — classifiers only.** The classifier half is
  gone: `extrema_surf/analytic_solvers.rs::classify_plane` = `s.gp_pln()`,
  `classify_sphere` = `s.gp_sphere()` (`Extrema_ExtSS.cxx:117-129` dispatch on `S1.GetType()`
  then `S1.Plane()`/`S1.Sphere()`), `reconstruct_line` = `c.gp_line()`
  (`Extrema_ExtCS.cxx:129-140`, `C.Line()`/`myS->Sphere()`/`myS->Plane()`), and the sampling
  bodies `is_planar` (8x8), `sphere_center` (circumcentre), `circumcenter`, `surface_normal`,
  `is_line` (6-point tangent), `sample_bounds` are deleted. In `geom2d`:
  `point_curve_extrema2d_all` takes the line arm from `gp_lin2d()`
  (`Extrema_GGExtPC.hxx:178-181`) and now computes `Mydist = L.Direction().Dot(L.Location()->P)`
  per `Extrema_ExtPElC2d.cxx:63-79` (the old code used the parameter increment relative to
  `d0(first_parameter)`, which is wrong on a trimmed line), the circle arm takes `gp_circ2d()`,
  and `curve_curve_extrema2d_all` dispatches pairs with `gp_lin2d()`/`gp_circ2d()`
  (`Extrema_ExtCC2d.cxx:104-145`, `GetType(C1)` x `GetType(C2)` +
  `Extrema_Curve2dTool::Line/Circle`); `is_line2d`/`circumcenter2d`/`classify_circle2d`/
  `line_of_curve2d` deleted. Observable: a straight-but-not-`Geom2d_Line` curve no longer falls
  into the line arm (OCCT reports `GeomAbs_BSplineCurve` and takes the general arm). Gates:
  `--all-targets` 0 err, `occt-geom2d --lib` 72/0, `occt-geom --lib` 143/0, `occt-topo --lib`
  1251/0, `step_obj_gates` 4/4 (117s), `export_data_obj` 15/15 v/f and bbox bit-identical.
- **B6c CLOSED (split out of B6b 2026-10-15; Scout done 2026-10-15; fixed 2026-10-15 as T-138 /
  package P008).** The `span/0.1` grid half of
  the old B6b line is not a classifier. The substituted-for engine is `Extrema_ExtPC2d`
  (`Extrema_ExtPC2d.hxx:31-38` = `Extrema_GGExtPC<Adaptor2d_Curve2d, Extrema_Curve2dTool,
  gp_Pnt2d, ...>`), whose body is `Extrema_GGExtPC.hxx`:
  - `initialise`: `mytolu = TheCurveTool::Resolution(*myC, Precision::Confusion())`, `type =
    TheCurveTool::GetType(theC)`, `mysample = 17` (`hxx:119-128`);
  - `Perform` (`hxx:133-...)`: the two endpoint distances (`myuinf`/`myusup`, skipped when
    infinite, `hxx:142-157`), then `switch (type)`: Circle/Ellipse/Parabola/Hyperbola/Line ->
    `myExtPElC.Perform(p, <curve>, t3d, myuinf, myusup)` (`hxx:160-182`; the port has
    `point_circle/ellipse/hyperbola/parabola/line_extrema2d`), Bezier -> `mysample =
    NbPoles()*2` + `IntervalPerform` (`hxx:183-189`), BSpline -> the knot/degree arms
    (`hxx:190-388`), `default:` -> `NbIntervals(GeomAbs_C2)` else `DeflCurvIntervals`,
    `mysample = max(17 * (u_sup-u_inf)/maxint, 3)` per interval, then `IntervalPerform`
    (`hxx:390-464`);
  - `IntervalPerform` (`hxx:590-620`) = `myExtPC.Initialize(mysample, myintuinf, myintusup,
    mytolu, mytolf); myExtPC.Perform(theP)` with `TheEPC` = `Extrema_EPCOfExtPC` ->
    `Extrema_GGenExtPC` (`Extrema_GGenExtPC.hxx`, `math_GlobOptMin` + `BoundedNewton`) and
    `AddSol` merging min/max by `Precision::PConfusion()` (`hxx:602-616`).
  What the port has instead: `extrema2d/analytic_solvers.rs::build_samples2d` (`span/0.1`,
  clamp 24..256) as the seed grid plus `newton_point_curve_all2d` (`:569`) and
  `curve_curve::newton_curve_curve_all2d` (`:361`) as the grid -> bracket -> Newton solvers.
  Live consumers of the substitute: `topo/brep_class/intersector_checks.rs:157`
  (`BRepClass_Intersector` = `Extrema_ExtPC2d`), `curve_ops::curve2d_closest_point` (=
  `Extrema_ExtPC2d`, A16 projection half / R2-18 remainder) and
  `extrema2d/general_extrema.rs:550` (`Extrema_GGenExtCC_ProjPOnC`, explicitly marked
  **UNPORTED (audit A16)** because the port's point/curve extrema has no range-restricted form).
  Scout conclusion: the slice is the faithful `Extrema_GGExtPC` + `Extrema_GGenExtPC` port
  (`math_GlobOptMin` is already ported in `occt-math` and already used by `general_extrema`), which
  also removes the A16 UNPORTED marker; it has to be split into at least the elementary arms,
  `IntervalPerform`/`GGenExtPC`, and the Bezier/BSpline knot arms.
  **FIXED (T-138 / package P008).** The substituted-for engine is now ported, not emulated:
  - new `extrema2d/ext_pc2d.rs` = `Extrema_GGExtPC<Adaptor2d_Curve2d, Extrema_Curve2dTool, ...>`
    (`Extrema_ExtPC2d.hxx:31-39`): `Initialize` (`hxx:114-127`), the endpoint distances
    (`hxx:142-157`), the elementary switch to `Extrema_ExtPElC2d::Perform` (`hxx:160-181`), the
    Bezier `mysample = NbPoles()*2` arm (`hxx:182-188`), the BSpline knot/degree arms
    (`hxx:190-388`), the `default:` `NbIntervals(GeomAbs_C2)` / `DeflCurvIntervals` arm with
    `mysample = max(17*(usup-uin)/maxint, 3)` (`hxx:390-469`), `IntervalPerform`
    (`hxx:590-614`), the end-point postprocessing (`hxx:473-502`), `AddSol`
    (`hxx:617-632`).
  - new `extrema2d/gfunc_ext_pc2d.rs` = `Extrema_GFuncExtPC` (`Extrema_GFuncExtPC.hxx`) with its
    singular `D1` arm (3-point differences) and `DF` handling (`hxx:133-...`), feeding
    `Extrema_GGenExtPC::Perform` = `math_FunctionRoots(myF, umin, usup, nbsample, tolu, tolf,
    tolf)` (`Extrema_GGenExtPC.hxx:160-173`) -- note the cxx here is `math_FunctionRoots`, not
    `math_GlobOptMin`, so the engine reuses the ported `occt-math` root finder.
  - `point_curve_extrema2d_all` = `ExtremaExtPC2d::new` (`Extrema_ExtPC2d::Perform`), list order
    preserved (OCCT neither sorts nor dedups); `point_curve_extrema2d` now returns `Option` and
    reproduces `Geom2dAPI_ProjectPointOnCurve::Perform` (`cxx:59-86`: done only when
    `IsDone() && NbExt() > 0`, then the smallest `SquareDistance`, first index wins ties).
  - `curve_curve_extrema2d_all` = `ExtremaExtCC2d::new` (`Extrema_ExtCC2d`) instead of the local
    pair dispatch + grid-Newton.
  - deleted whole: `build_samples2d`, `newton_point_curve_all2d`, `fval2d`, `fprime2d`,
    `dist22d`, `solve_f_zero2d`, `refine_seed2d`, `pair2d_c`, `dedupe_sort2d`,
    `newton_curve_curve_all2d`, `param_in2d`.
  - the **A16 UNPORTED marker is gone**: `general_extrema::proj_p_on_c` is now
    `ExtremaExtPC2d::new_range` + the min `SquareDistance` scan of
    `Extrema_GGenExtCC_ProjPOnC` (`Extrema_GGenExtCC.hxx:249-263`, caller `:771-773`), i.e. the
    range-restricted form the old note said did not exist.
  - `topo/shape_fix_compose_shell/split_wire.rs` now mirrors
    `ShapeFix_ComposeShell.cxx:1046-1087` literally: `apar` starts at `firstPar` and is overridden
    only when the 2D extrema run is done with at least one solution -- the `Option` return is what
    makes that cxx guard expressible.
  - four existing lib assertions were corrected to the cxx order semantics (no new test): the
    `Extrema_ExtPElC2d` circle order (far point first, `Extrema_ExtPElC2d.cxx:113-150`), the
    `Extrema_ExtElC2d(circ, circ)` nested `P1[NoSol] x P2[kk]` order (`Extrema_ExtElC2d.cxx:334-364`),
    and parallel lines reporting `IsParallel()` with zero extrema (`Extrema_ExtElC2d.cxx:63-70` +
    `Extrema_ExtCC2d.cxx:541-544`).
  - gates: `--all-targets` 0 err (occt-geom2d, occt-geom, occt-topo), `occt-geom2d --lib` 72/0,
    `occt-geom --lib` 143/0, `occt-topo --lib` 1251/0, `step_obj_gates` 4/4 with the 23 area
    ratios identical to the T000 baseline (as read at P016: 22 x 1.0000, Shape-2 0.9999; the
    2026-10-16 re-read on the frozen tranche gives 20 x 1.0000 plus Shape-2 0.9999 / T0M 1.0017 /
    TDB 1.0015 -- T0M and TDB moved under the later shared-path packages, tolerances unchanged),
    `export_data_obj` 15/15
    with v/f and bbox bit-identical.
- **B6e CLOSED 2026-10-15 (P017 / T-147, commit `e3a4335e`).** `extrema2d/curve_curve.rs::locate_extcc2d`
  (`:255-262`) now builds `ExtremaLocateExtCC2d::new(c1, c2, u0, v0)` and the self-invented clamped
  32-iteration Newton (`refine_curve_curve2d_newton` + `f_and_jac2d`, whose Jacobian had all four
  sign-bearing terms flipped relative to `Extrema_GFuncExtCC::Values`) is deleted. The cxx chain is
  landed: `Extrema_LocateExtCC2d.cxx:26-82` -> `Extrema_LocECC2d` = `Extrema_GenLocateExtCC`
  (`Extrema_GenLocateExtCC.hxx:47-110`: per-axis range normalisation to `[FirstParameter,
  LastParameter]` / `[FirstParameter, LastParameter]`, `TolU`/`TolV` from each curve's
  `Adaptor::Resolution(Confusion)`, `MathFunctionSetRoot::perform_with_bounds`, then the
  `F.NbExt() > 0` gate and the `F.Value < 1e-10` re-check) -> `Extrema_CCLocFOfLocECC2d` =
  `Extrema_GFuncExtCC` 2D. New `extrema2d/gfunc_ext_cc2d.rs` (`GFuncExtCC2d`, the
  `MathFunctionSetWithDerivatives` with two equations in two variables: `SearchOfTolerance`
  `hxx:150-197`, the `GetType` switches `hxx:263-291`, `Value` `hxx:359-550` with the singular
  `D1` three-point / `DN` reconstruction, `Values` `hxx:568-827` with the perturbed three-point
  Jacobian, `GetStateNumber` `hxx:832-870`, `Points`, `SubIntervalInitialize`) and
  `extrema2d/gen_locate_ext_cc2d.rs` (`GenLocateExtCC2d` + `ExtremaLocateExtCC2d`). Gates:
  `occt-geom2d --lib` 72/0, `occt-geom --lib` 143/0, `occt-topo --lib` 1251/0, `step_obj_gates`
  4/4 (248.5s), `export_data_obj` 15/15 bit-identical; no test added. Still open in the same family:
  B6f (now closed below) and B6d.
- **B6e 3D twin CLOSED 2026-10-15 (P019 / T-149, commit `616ea0c1`).** `occt-geom/src/extrema_cc/curve_curve.rs::locate_extcc` (`:152`) now builds
  `ExtremaLocateExtCC::new(c1, c2, u0, v0)`, and the self-invented clamped 32-iteration Newton
  (`refine_curve_curve` + `f_and_jac` in `poly_roots.rs`, whose Jacobian had all four sign-bearing
  terms flipped relative to `Extrema_GFuncExtCC::Values`) is deleted. The cxx chain is landed:
  `Extrema_LocateExtCC.cxx:24-74` -> `Extrema_LocECC` (`Extrema_LocECC.hxx:27-30`) =
  `Extrema_GenLocateExtCC` (`Extrema_GenLocateExtCC.hxx:49-110`: per-axis range normalisation to
  `[FirstParameter, LastParameter]`, `TolU`/`TolV` from `Curve::resolution(Confusion)`,
  `MathFunctionSetRoot::perform_with_bounds`, then the `F.NbExt() > 0` gate and the
  `F.Value < 1e-10` re-check) -> `Extrema_CCLocFOfLocECC` (`Extrema_CCLocFOfLocECC.hxx:29-36`) =
  `Extrema_GFuncExtCC` 3D. New `extrema_cc/p_on_curv.rs` (`POnCurv` = `Extrema_POnCurv`,
  `Extrema_POnCurv.hxx:25-60`), `extrema_cc/gfunc_ext_cc.rs` (`GFuncExtCC`, the
  `MathFunctionSetWithDerivatives` with two equations in two variables: `SearchOfTolerance`
  `hxx:150-197`, the `GetType` switches `hxx:263-291`, `Value` `hxx:359-550` with the singular
  `D1` three-point / `DN` reconstruction, `Values` `hxx:568-827` with the perturbed three-point
  Jacobian, `GetStateNumber` `hxx:832-870`, `Points`, `SubIntervalInitialize`) and
  `extrema_cc/gen_locate_ext_cc.rs` (`GenLocateExtCC` + `ExtremaLocateExtCC`). Gates:
  `occt-geom2d --lib` 72/0, `occt-geom --lib` 143/0, `occt-topo --lib` 1251/0, `step_obj_gates`
  4/4 (181.8s), `export_data_obj` 15/15 bit-identical; no test added. The B6e line (2D + 3D twin)
  is fully closed.
- **B6f CLOSED (2026-10-15, package P020 / T-150, commit `2452a51c`).**
  `extrema2d/curve_curve.rs::curve_curve_intersections2d` was a 33x33 uniform parameter grid over
  `bound()` that pushed every pair closer than `tol` under a local `1e-6` dedup -- no OCCT control
  flow does that, and an intersection off the lattice was silently dropped or replaced by a
  neighbouring sample. It now runs the faithful chain of `Geom2dAPI_InterCurveCurve.cxx:54-106`:
  `Geom2dInt_GInter(AC1, AC2, Tol, Tol)` through the already-ported
  `Geom2dIntGInter::perform_curves`, `NbPoints() == 0` unless `IsDone()`, then `Point(i).Value()`
  for `i` in `1..=NbPoints()`. The `bound()` helper in `extrema2d/analytic_solvers.rs`, which
  served only that grid, is deleted. Gates: `occt-geom2d --lib` 72/0, `occt-geom --lib` 143/0,
  `occt-topo --lib` 1251/0, `step_obj_gates` 4/4 (205.1s), `export_data_obj` 15/15 bit-identical;
  no test was added or realigned (`curve_curve_intersections_found` still passes unchanged).
- **B6g CLOSED (2026-10-15, package P018 / T-148, commit `8cad6a01`).** `extrema2d/curve_curve.rs`'s
  `line_line_extrema2d`, `line_circle_extrema2d` and `circle_circle_extrema2d` (plus the helpers
  `line2d_value`, `circle_param2d`, `dedupe_sort_cc2d`, which served nothing else) are deleted. They
  were a second, module-level copy of the `Extrema_ExtElC2d` arms that `extrema2d/ext_el_c2d.rs`
  already carries, and they had diverged: parallel lines returned a point pair where
  `Extrema_ExtElC2d.cxx:63-70` only sets `IsParallel()` and `SquareDistance(1)`, line/circle appended
  intersection points the analytic arm never produces (`cxx:105-167`), and circle/circle added an
  invented intersection pass plus a `dedupe_sort` that `cxx:274-368` does not perform. The five
  `extrema2d/tests.rs` assertions now drive `ExtremaExtElC2d::new_line_line` / `new_line_circle` /
  `new_circle_circle` and were rewritten to the cxx counts and ordering; no test was added.
  Gates: `occt-geom2d --lib` 72/0, `occt-geom --lib` 143/0, `occt-topo --lib` 1251/0,
  `step_obj_gates` 4/4 (150.7s), `export_data_obj` 15/15 bit-identical.
- **B6d CLOSED (2026-10-15, package P021 / T-151, commit `11b12e0e`).** `occt-geom/src/extrema_ss.rs`
  was a single-file solver of its own: analytic sphere/sphere and sphere/plane arms
  (`sphere_sphere_extrema`, `plane_sphere_pairs`, `plane_sphere_extrema`, `sphere_plane_extrema`), the
  helpers that served only them (`perpendicular`, `project_on_plane`), and an invented Newton for the
  general arm (`ss_f`, `solve4`, `solve_surface_surface`). `Extrema_ExtSS::Perform` (`Extrema_ExtSS.cxx:100-255`)
  has exactly one analytic arm: `S1.GetType() == GeomAbs_Plane && myStype == GeomAbs_Plane` ->
  `myExtElSS.Perform(S1.Plane(), myS2->Plane())` (`cxx:122-129`); every other type pair, sphere/sphere
  and sphere/plane included, falls through to `Extrema_GenExtSS(S1, *myS2, 20, 20, mytolS1, mytolS2)`
  (`cxx:131-175` and `cxx:180-221`). All of that is what the port now carries:
  `extrema_ss/mod.rs` (`ExtremaExtSS`), `extrema_ss/ext_el_ss.rs` (`ExtremaExtElSS`,
  `Extrema_ExtElSS.cxx:40-54`; the `gp_Sphere` overloads only set flags and throw
  `Standard_NotImplemented`, `cxx:62-69`/`cxx:77-83`, so they are not ported),
  `extrema_ss/gen_ext_ss.rs` (`ExtremaGenExtSS`, `Extrema_GenExtSS.cxx:119-470`: the two 20x20 sample
  axes, the full point-pair scan for the min/max seeds, `math_BFGS(4)` for the minimum and
  `math_FunctionSetRoot` for the maximum), `extrema_ss/func_ext_ss.rs` (`ExtremaFuncDistSS`,
  `ExtremaFuncExtSS` = `Extrema_FuncExtSS.cxx:52-241`) and `gp_Pln::{signed_distance_pnt,
  signed_distance_pln, square_distance_pln}` (`gp_Pln.hxx:345-379`, `:173-177`). An unbounded plane is
  no longer special-cased: `GeomAdaptor_Surface` hands `Geom_Plane::Bounds` (`+/-2e100`) straight to
  `Extrema_GenExtSS`, so the port samples and converges in the `1e100` range exactly like OCCT. The
  three baseline assertions that encoded the old analytic answer are rebased onto the OCCT 8.0.0
  reading pinned by `occt_probe --ssx` (rebuilt and rerun with this package):
  `SSX sph-pln done=1 isPar=0 nbExt=1 sqDist=7.992002e+200 uv1=(0.0015707963267948967,-1.5700109286314992)
  uv2=(1.999e+100,1.999e+100)`, `SSX pln-sph` symmetric, `SSX sph-pln-win nbExt=0`,
  `SSX api-sph-pln nbExt=1 dist=2.827012911183817e+100`,
  `SSX gen-sph-pln nbExt=2 [1] 2.2138509695290867e+198 [2] 7.992002e+200` and
  `SSX sph-sph nbExt=2 sqDist=4.0000000000003126 / 36`; the third assertion moved to the
  sphere/sphere pair, which is well posed. Gates: `occt-core --lib` 283/0, `occt-geom2d --lib` 72/0,
  `occt-geom --lib` 143/0, `occt-topo --lib` 1251/0, `step_obj_gates` 4/4 (114.1s),
  `export_data_obj` 15/15 bit-identical; no test was added. **B6h registered here, not fixed**:
  `extrema_ss::surface_surface_extrema` returns the first `mySqDist` entry while its OCCT counterpart
  `GeomAPI_ExtremaSurfaceSurface` stores the index of the smallest `SquareDistance` in `myIndex`
  (`GeomAPI_ExtremaSurfaceSurface.cxx:71-91`, used by `NearestPoints` `cxx:203-209` and
  `LowerDistance` `cxx:225-231`), so the two differ whenever `NbExt() > 1`; likewise
  `point_on_s1`/`point_on_s2` expose what OCCT keeps private (`Extrema_ExtSS.hxx:99-102`, public
  access only through `Points`). Both sit on APIs with no production caller (unit tests only). This
  closes the B6 family (B6a/B6b/B6c/B6d/B6e 2D+3D/B6f/B6g).

### B7. Step/IGES-adjacent live divergences
- `crates/occt-topo/src/mesh_to_brep.rs` ✅ **CLOSED 2026-10-16 (P025 / T-155, `8d77b13a`)**:
  the three port-invented rules are gone. `BRepBuilderAPI_MakeShapeOnMesh`
  (`BRepBuilderAPI_MakeShapeOnMesh.cxx:74-263`) is ported whole as `MakeShapeOnMesh`: nodes 1:1
  through `BRepBuilderAPI_MakeVertex` at `BRepLib::Precision()` = `Precision::Confusion()` = 1e-7
  (`cxx:96-111`, `BRepLib_MakeVertex.cxx:26-31`, `BRepLib.cxx:84`), shared edges keyed on the
  canonical node-index pair with `NCollection_IndexedDataMap::Add` first-value semantics
  (`cxx:31-55`, `cxx:150-179`), index- and `gp::Resolution()`-degenerate triangle skips
  (`cxx:126-141`), one planar face per triangle with the edge re-orientation of `cxx:213-228`, the
  `BRepAdaptor_Curve::Line()` `d1 x d2` normal under `Precision::SquareConfusion()` (`cxx:238-247`),
  the plane through `Node(anIdx[0])` (`cxx:255`) and the closing `Done()` + `TopoDS_Compound`
  (`cxx:182-184`, `cxx:257-263`). `triangulation_to_brep` keeps only `LOCAL:`-marked extensions
  (`weld_nodes`, the closed 2-manifold test, `OrientClosedSolid`, the winding pass), and the vertex
  tolerance is the faithful 1e-7, which is load-bearing because `pave_ff.rs:158`,
  `pave_split_blocks.rs:74` and `brep_fill.rs:213` fold it into a working precision with `max(...)`.
  Live consumers left registered: `brepfeat/features.rs:367` (`shape_mesh_to_brep`) and
  `rwmesh/mod.rs:24` / `rwmesh/scene.rs:307` (`triangulation_to_brep`) -- those entry points and
  their `LOCAL:` passes have no OCCT class. Correction: the `32^3` `voxel_boolean` fallback that
  `_audit/bop-boolean.md` item 2 places at `bop_builder_core.rs:85-87` is **no longer in the tree**
  (that site is now a warning-only mesh-volume post-check; `voxel_boolean` survives only inside
  `boolean_ops.rs` and the comment at `bop_builder_planar.rs:47`), so that audit item needs a
  re-check rather than a fix.
- `step/write_context.rs:314` ✅ **CLOSED 2026-10-15 (P010)**: the periodic arm of
  `GeomToStep_MakeBoundedCurve` (`cxx:44-51`, `BS = B->Copy(); BS->SetNotPeriodic()`) is wired.
  `write_context.rs::bounded_bspline_of_curve` reproduces `cxx:42-62` (the curve's own poles /
  knots / weights / degree, then the unperiodize step keyed on `C->IsPeriodic()`, `cxx:45`) and
  `open_periodic_bspline` is the `Copy()` + `SetNotPeriodic()` pair; the copy has to restore
  `periodic` because the port's `GeomBSplineCurve::new` / `rational` constructors always build the
  non-periodic form. The trimmed-B-spline arm of the same writer (`MakeCurve.cxx:71-76`) now also
  restores `B->IsPeriodic()` on the basis copy (`cxx:73`) so `Geom_BSplineCurve::Segment` reaches
  its periodic branch (`Geom_BSplineCurve.cxx:606-618`). The UNPORTED note is deleted.
- `step/format.rs:472,526` ✅ **CLOSED 2026-10-15 (P011)**: both STEP B-spline writers now
  derive their flags from the geometry. `write_bspline_curve` writes `closed` =
  `Geom_BSplineCurve::IsClosed()` and `knot_spec` =
  `knot_spec_token(KnotDistribution())` (`GeomToStep_MakeBSplineCurveWithKnots_gen.pxx:42-87`);
  `write_bspline_surface` names the spec only when `UKnotDistribution` and `VKnotDistribution`
  agree, otherwise `.UNSPECIFIED.` (`MakeBSplineSurfaceWithKnots.cxx:127-150`). New dependencies:
  `occt-core/src/bspl/knot_analysis.rs` = `BSplCLib::KnotForm`/`MultForm`/`KnotAnalysis`
  (`BSplCLib.cxx:602-753`), `GeomBSplineCurve::{is_closed,knot_distribution}`
  (`Geom_BSplineCurve_1.cxx:146-149,363-366`), `GeomBSplineSurface::{u,v}_knot_distribution`
  (`Geom_BSplineSurface_1.cxx:1851-1862`). Both UNPORTED notes deleted. The divergence was real
  against the OCCT-written oracle files: `data/occ/ATU01038.step` carries 245
  `PIECEWISE_BEZIER_KNOTS` + 106 `QUASI_UNIFORM_KNOTS`, `bottom.step` 934, `data/*.step`
  8-70 each, where the port wrote `.UNSPECIFIED.` for every entity.
- `geom2d_convert.rs:77`: parameterisation-less `CurveToBSplineCurve` entry has no conic arms
  (`Geom2dConvert::CurveToBSplineCurve(C)`); the conic converters already exist.

### B8. `occt-geom` / `occt-core` singletons
- ✅ **osculating_surface.rs early returns — CLOSED 2026-10-15 (T-139, package P009)**: `osculating_surface.rs`'s four early returns (old lines `:141,202,238,299`). The
  old `init` bodies did `clear_flags(); return;` on a failed `build_osculating`, where OCCT only
  `clearOsculFlags()` (`Geom_OsculatingSurface.cxx:206-213`, `:273-276`, `:308-312`, `:373-376`) and
  keeps the `for` running; each of the five arms also owns its own `OsculSurf` local (reassigned on
  the next iteration), which the port collapsed into one shared `oscul_ok`. Fixed by adding
  `build_chain` (`cxx:187-206`, returns `Option<(_, k)>` for the `while (IsQPunc)` chain) and
  dispatching per arm: `if along[0] { A; if along[1] && ok { B } } else if let Some(..) { C } else
  { clear_flags() }` in both knot loops. Reachable via `offset_surface.rs:56` -> `brep_offset`.
  Gates: `occt-geom --lib` 143/0, `occt-topo --lib` 1251/0, `step_obj_gates` 4/4, `export_data_obj`
  15/15 bit-identical.
- ~~`surface_fit.rs:13-40` `fit_plane` lacks the `gz < Tol` gate (`GeomLib_IsPlanarSurface.cxx:44-45`);
  reachable via `brep_offset/curve_face_offset.rs:489`.~~ **CLOSED 2026-10-16 as T-154 / package P024**,
  and the registration above turned out to name the wrong OCCT class twice over.
  (a) `GeomLib::IsPlanarSurface`, `GeomLib::Inertia`, `GeomLib::AxeOfInertia` and `Controle` have no
  counterpart and no reachable branch in the port (`grep -rn "IsPlanarSurface\|AxeOfInertia"
  crates/**/*.rs` = zero hits; the `Inertia`/`Controle` hits are `BRepGProp_Gauss` accumulators,
  `fillet_var` comments and `math_NewtonMinimum.cxx:230`). (b) The live consumer's real OCCT
  counterpart is `BRepLib_FindSurface` via `BRepLib_MakeFace(W, OnlyPlane)`: `brep_offset`'s
  `offset_wire_3d` stands in for `BRepOffsetAPI_MakeOffset::BuildDomains`
  (`BRepOffsetAPI_MakeOffset.cxx:184-192`), which calls `BRepBuilderAPI_MakeFace(WorkWires.First(),
  true)` (`BRepBuilderAPI_MakeFace.cxx:215-224`) -> `BRepLib_MakeFace(W, OnlyPlane)`
  (`BRepLib_MakeFace.cxx:188-256`) -> `BRepLib_FindSurface(W, -1, OnlyPlane, true)`
  (`BRepLib_FindSurface.cxx:52-628`). Both are now ported
  (`crates/occt-topo/src/brep_lib_find_surface.rs`, `brep_lib_make_face.rs::make_face_from_wire`)
  and `offset_wire_3d` is wired to them, so the 4-samples-per-edge least-squares fit is off the
  `brep_offset` path. (c) Left registered, not fixed: `occt-geom/src/surface_fit.rs`'s
  `fit_plane` / `fit_surface_kind` / `fit_surface` are port-only (no OCCT class, no production
  caller left after the rewire) and now have only their own tests as consumers.
- `intana/analytic_intersections.rs:42,47,48,120,145,152` and `intana/line_torus.rs:336,340,488,489`
  invent quartic root thresholds and an absolute `1e-14` cone/cone tolerance where OCCT passes
  `done` out of `math_DirectPolynomialRoots` and uses `EstimDist`/`aTolAng`
  (`math_DirectPolynomialRoots.cxx:344-393,583-587`, `IntAna_QuadQuadGeo.cxx:1440,1458-1475,1524`);
  reachable via `intcurvesurface/solvers.rs:456`.
- `extrema_cc/mod.rs:9` **STALE (re-verified 2026-10-15)**: the module header still claims the seed
  set is an invented uniform grid, but `curve_curve_extrema_all_range` already routes the non-analytic
  pairs to `GGenExtCC` = `Extrema_GGenExtCC` with `GlobOptFuncCCC2` + `math_GlobOptMin`
  (`extrema_cc/general_extrema.rs:245-292`); the `is_line`/`classify_circle` sampler the header
  describes is gone. Only the doc comment needs correcting (next doc pass).
- `occt-math/globoptmin.rs:590` **STALE (re-verified 2026-10-15)**: `math_NewtonMinimum::Perform` local
  arm marker; the faithful port exists at `occt-math/newton_minimum.rs` (re-exported `lib.rs:59`), so
  the `:590` and `occt-math/newton.rs:81,91` markers are stale (doc-only).
- `occt-core/src/quantity/mod.rs:33`: `Color::from_name` returns sRGB columns and only 23 names
  (`Quantity_ColorTable.pxx:185`, `Quantity_Color.hxx:49-51`); test-only reachability.
- `occt-core/src/geom/polyline_simplify.rs:15` (RDP) is not in OCCT but is live via
  `intpatch_trace.rs:261`.

## C. UNPORTED gaps (cxx branch exists, no Rust counterpart)

Count: ~115 markers over 55 files. Bucketed in the Scout sweep; the actionable subset:

- **C-A (leaf, deps present, reachable)** — the B1 items above, plus `geom2d_convert.rs:77` from
  B7 (the sibling `step/write_context.rs:314` and `step/format.rs:472,526` lines are CLOSED in P010
  and P011); and `extrema_cc/mod.rs:9`, `globoptmin.rs:590` from B8.
- **C-B (needs an unported layer)** — `shhealing/shape_fix_face.rs` (5 markers, need
  `ShapeFix_Wire::Perform`), `step/read_topology/step_to_geom_curve{2d}.rs` (closed/reversal flags,
  2D rational weights), `brep_check.rs:149` (`BRepCheck_Edge::InContext` ->
  `BRepLib_ValidateEdge`), `validate.rs:3,140` (whole `BRepCheck_Analyzer`), `brep_extrema.rs:269,289`
  (`SolidExplorer` for non-solids), `intf/*` (3D conic + `Extrema_ExtElC`), `adv_app2var/*`,
  IGES transfer leftovers, `deflection.rs:9` (`IMeshData` overloads).
- **C-C (out of scope / lite substitutes)** — `viz_scene`, `draw`, `xcaf`, `render_svg`, `hlr`,
  `occt-core` `csg`/`delaunay`/`triangulate`, `gp` `DumpJson` family.
- **C-D (stale markers)** — `shhealing/shape_build_edge.rs:104`, `int_tools_segpln.rs:153`,
  `bop_builder_core.rs:83`, `iges/emit_face.rs:243`, `pcurve_full/make_pcurve.rs:259`,
  `meshing/geom_tool.rs:188,228`, `wireframe.rs:12`, `extrema_ss.rs:64,138`, many more.

## D. Banned constructs found (must not be extended; listed for the report)

- Area ratio / length filter / volume gate: **zero hits** in `occt-geom`, `occt-geom2d`, `occt-core`.
- Sampled/heuristic substitute for a downcast or a closed form: B1a, B1b, B2, B6, plus
  `meshing/deflection.rs:109,138`, `meshing/range_splitter/param_set.rs:203`, `hlr.rs:129-156`,
  `gprop_analytic/mass_properties.rs:405`, `shape_metrics.rs:146`, `solid_union.rs:99-111,144-182`,
  `brep_offset/curve_face_offset.rs:509`-adjacent, `brep_pipe.rs:24-53,78-95`,
  `brepfeat/features.rs:326-370`, `brep_class3d` probe half is closed. B3 closed in T-134.
  B4 closed in T-135 (its `face_sample_point` family and the sampled point-vs-solid classifiers in
  `builder_solid.rs` / `draft_solid.rs` / `bop_build_solids.rs` are gone; the only survivor is the
  port-only `bop_build_solids::face_sample_point` feeding `orient_face_outward`). B6a closed in
  T-136 (the three fabricated Extrema `fallback_*` pairs in `extrema_surf/numeric_extrema.rs` and
  `extrema_ss.rs` are deleted; the minimum-distance entries propagate `!IsDone()` / `NbExt() == 0`
  as `None`, and `bean_face_exact::ComputeUsingExtremum` takes the `IntTools_BeanFaceIntersector.cxx:956`
  guard instead of a fabricated pair). B6b closed in T-137 (both `analytic_solvers` families now
  take `gp_pln()`/`gp_sphere()`/`gp_line()`/`gp_lin2d()`/`gp_circ2d()`; the sampled classifiers and
  their helpers are deleted). Split out of B6b as **B6c**, closed in T-138 (the `span/0.1` seed grid
  `build_samples2d` + `newton_point_curve_all2d` / `newton_curve_curve_all2d` are deleted and
  replaced by the faithful `Extrema_ExtPC2d` = `Extrema_GGExtPC` + `Extrema_GFuncExtPC` --
  `extrema2d/ext_pc2d.rs`, `extrema2d/gfunc_ext_pc2d.rs`; the A16 UNPORTED marker on
  `general_extrema::proj_p_on_c` came off with them). B6e/B6f/B6g registered in the same package
  (`locate_extcc2d`'s self-invented Newton vs `Extrema_GenLocateExtCC` +
  `math_FunctionSetRoot`; the 33x33 `curve_curve_intersections2d` grid; the duplicated
  `line_line`/`line_circle`/`circle_circle_extrema2d` copies of `Extrema_ExtElC2d`). **B6e CLOSED
  2026-10-15 as P017 / T-147** (the 2D chain `Extrema_LocateExtCC2d` -> `Extrema_GenLocateExtCC`
  -> `Extrema_GFuncExtCC` 2D is landed in `extrema2d/gfunc_ext_cc2d.rs` +
  `extrema2d/gen_locate_ext_cc2d.rs`; the invented Newton and its Jacobian are deleted); **B6g
  CLOSED 2026-10-15 as P018 / T-148** (the three module-level `Extrema_ExtElC2d` arm copies in
  `extrema2d/curve_curve.rs` are deleted and the five `extrema2d/tests.rs` assertions now
  drive `extrema2d/ext_el_c2d.rs` directly); and the **B6e 3D twin CLOSED 2026-10-15 as P019 /
  T-149** (the 3D chain `Extrema_LocateExtCC` -> `Extrema_LocECC` = `Extrema_GenLocateExtCC` ->
  `Extrema_CCLocFOfLocECC` = `Extrema_GFuncExtCC` 3D is landed in `extrema_cc/p_on_curv.rs`,
  `extrema_cc/gfunc_ext_cc.rs` and `extrema_cc/gen_locate_ext_cc.rs`; the invented
  `refine_curve_curve` + `f_and_jac` are deleted).
  **B6d** (found in T-137, **CLOSED 2026-10-15 as P021 / T-151**): the analytic sphere/sphere and
  sphere/plane arms of `surface_surface_extrema_all` had no counterpart in
  `Extrema_ExtSS::Perform` (only `Plane` x `Plane` is analytic there); they are deleted, replaced
  by `extrema_ss/{mod,ext_el_ss,gen_ext_ss,func_ext_ss}.rs`.
- Self-invented engines with no OCCT counterpart and no production caller: `boolean_ops.rs:173`
  (+`:206/:211/:216/:280`), `solid_union.rs` voxel union.

## E. Ranked candidate Worker slices (advisory input to T003)

1. ~~**B1 (whole family, one file, live BOP path)** — `wire_splitter_block.rs` only.~~
   **ALREADY CLOSED 2026-10-15 in T-131 (commit `67d7f01f`)**; re-verified this cycle, all three arms
   are faithful (see B1). The list below is the historical advisory text.
   delegate `uv_tolerance_2d` to `approx_same_parameter::{u,v}_resolution` and delete
   `fd_uv_tolerance`; replace `curve2d_radius` with `GeomLPropCLProps2d`; take `refine_angle_2d`'s
   `aTV` from `vertex_parameter_on_face`. High confidence, 3 cxx arms, no new predicate.
2. ~~**B2** `pc_curve_kind` -> `Geom2dAdaptor_Curve::GetType()` equivalent. 4 call sites, touches
   bbox/pcurve bounds.~~ **DONE 2026-10-15 (T-132, package P002).**
3. ~~**B6b** Extrema analytic-solver sampling classifiers (the B6a half, the fabricated `fallback_*`
   pairs, closed 2026-10-15 as T-136 / P006). Still open: the sampling classifiers in the two
   `analytic_solvers.rs` files; touches `bean_face_exact`'s neighbours, Scout first.~~
   **DONE 2026-10-15 (T-137, package P007).** Continuation was **B6c** (the `span/0.1` seed grid +
   grid-Newton general arms: `build_samples2d`, `newton_point_curve_all2d`,
   `newton_curve_curve_all2d`, live via `brep_class/intersector_checks.rs:157` and
   `extrema2d/general_extrema.rs:550`), which is **DONE 2026-10-15 as T-138 / package P008** --
   the grid and its solvers are deleted; the engine is `Extrema_ExtPC2d` = `Extrema_GGExtPC`
   (`extrema2d/ext_pc2d.rs`) + `Extrema_GFuncExtPC` (`extrema2d/gfunc_ext_pc2d.rs`) feeding
   `math_FunctionRoots`. Remaining B6 items, in no fixed order: **B6d** (the extra
   analytic sphere arms in `surface_surface_extrema_all`), **B6f**/**B6g** (the 33x33 intersection
   grid; the duplicated `Extrema_ExtElC2d` arm copies). **B6e is CLOSED 2026-10-15 as T-147 /
   package P017** (`extrema2d::locate_extcc2d` now runs the faithful `Extrema_LocateExtCC2d` /
   `Extrema_GenLocateExtCC` -> `Extrema_GFuncExtCC` 2D chain), and its 3D
   twin `extrema_cc/curve_curve.rs::locate_extcc` is **CLOSED 2026-10-15 as T-149 / package P019**
   (the faithful `Extrema_LocateExtCC` / `Extrema_GenLocateExtCC` -> `Extrema_GFuncExtCC` 3D chain,
   with `refine_curve_curve` + `f_and_jac` deleted); **B6g is CLOSED 2026-10-15 as T-148 / package
   P018** (the three duplicated `Extrema_ExtElC2d` arm copies and their three helpers deleted;
   the tests drive `ExtremaExtElC2d` directly) and **B6f is CLOSED 2026-10-15 as T-150 / package
   P020** (`curve_curve_intersections2d`'s 33x33 grid replaced by the
   `Geom2dAPI_InterCurveCurve` -> `Geom2dInt_GInter` chain, `bound()` deleted), so the B6 family
   is down to **B6d**, which is also **CLOSED 2026-10-15 as T-151 / package P021**: the extra analytic sphere/sphere and sphere/plane arms, the helpers serving only them, and the invented `solve_surface_surface` Newton are deleted, replaced by `Extrema_ExtSS` / `Extrema_ExtElSS` / `Extrema_GenExtSS` / `Extrema_FuncExtSS` (the only analytic arm is `Plane` x `Plane`), with the unbounded-plane reading pinned by a rebuilt `occt_probe --ssx`. The B6 family is now fully closed; **B6h** (below) is registered as its small leftover. **The next
   slice is #1 below (B1)**: it outranks the B6 leftovers on centrality, confidence and blast radius.
4. ~~**B5** `find_nearest_valid_point` `aD1Mag` arm. Small, self-contained, live pave filler.~~
   **DONE 2026-10-15 (T-133).**
6. **B8 osculating_surface early returns.** **DONE 2026-10-15 as T-139 / package P009** (four
   `clear_flags(); return;` -> OCCT's `clearOsculFlags()` + continue, plus the per-arm `build_chain`
   local). Next picks, in order: **B8-intana** (`intana` quartic/direct-polynomial thresholds, live via
   `intcurvesurface/solvers.rs:456`, closed as P012-P016), **B7-writer-4** (`geom2d_convert.rs:77`
   conic arms, closed as P022 / T-152), **B7-mesh_to_brep** (`BRepBuilderAPI_MakeShapeOnMesh`,
   closed as P025 / T-155),
   then the B6 leftover
   (**B6h**: `surface_surface_extrema` returns the first `mySqDist` entry where
   `GeomAPI_ExtremaSurfaceSurface` returns the `myIndex` minimum, plus the two exposed private
   point sequences; closed as P023 / T-153) and the `B8 surface_fit` item (closed as P024 / T-154).
   The B7
   STEP-writer lines are both
   **DONE 2026-10-15**: `step/write_context.rs:314` `SetNotPeriodic` as package P010, and
   `step/format.rs:472,526` `closed` + `knot_spec` as package P011, with B8-intana closed as
   P012-P016, **B6e closed as P017 / T-147** (2D) and its 3D twin
   `extrema_cc/curve_curve.rs::locate_extcc` closed as P019 / T-149, **B6g closed as
   P018 / T-148** and **B6f closed as P020 / T-150** (the `curve_curve_intersections2d` 33x33
   grid replaced by the `Geom2dAPI_InterCurveCurve` -> `Geom2dInt_GInter` chain), **B6d closed as
   P021 / T-151** (the extra sphere/sphere and sphere/plane analytic arms replaced by the
   `Extrema_ExtSS` / `Extrema_GenExtSS` chain) and **B7-writer-4 closed as P022 / T-152**: the
   partial second `curve_to_bspline_curve` (Bezier and B-spline arms only, `None` for every conic
   and for the trimmed case, no caller) is deleted, leaving the single faithful
   `curve_to_bspline_curve_bspl` port of `Geom2dConvert.cxx:181-449`. **The `surface_fit` pick
   closed as P024 / T-154**: the faithful counterpart of the live consumer turned out to be
   `BRepLib_FindSurface` through `BRepLib_MakeFace(W, OnlyPlane)`, now ported, so `offset_wire_3d`
   no longer runs the 4-samples-per-edge least-squares fit; the registered `GeomLib::IsPlanarSurface`
   / `GeomLib::Inertia` / `Controle` classes have no counterpart in the port at all. **The
   `B7-mesh_to_brep` pick closed as P025 / T-155**: `mesh_to_brep.rs` now ports
   `BRepBuilderAPI_MakeShapeOnMesh` (`BRepBuilderAPI_MakeShapeOnMesh.cxx:74-263`) whole, so the
   `1e-9` vertex dedup grid and the `0.0` vertex tolerance are gone from the conversion and only the
   `LOCAL:`-marked welding / closed-manifold / orientation passes remain. **The `ReparamCurve2d`
   `SameRange` pick closed as P026 / T-156, with the verdict the card ordered first**: the wrapper
   already delegated the curve kind and `reparam_curve2d` already returned a concrete
   `Geom2dBSplineCurve`, so the line was *not* stale-closed -- what was still open is that
   `pcurve_full/surface_projector/reparam.rs` and `shhealing/wire_fix/pcurve_range.rs` each carried
   their own copy of the `GeomLib::SameRange` conversion tail, and the `wire_fix` copy
   (`line_trim_to_reparam_bspline`, `unequal_trim_window` -> `same_range_unequal_line` /
   `same_range_unequal_other`) dispatched on the curve kind inside the unequal-span arm, where
   `GeomLib.cxx:924-969` has one path for every kind. Both copies collapse into the single
   `trimmed_reparam_bspline` (`GeomLib.cxx:912-921` + `:962-968`) plus `same_range_window`
   (`cxx:934-959`), and `ReparamCurve2d::bspline_copy2d` is gone (the wrapper exists only where
   `Geom2dConvert::CurveToBSplineCurve` declined, `Geom2dConvert.cxx:442-444`, so the trait default
   `None` is the right answer). Remaining pick: the missing data-driven STEP-write harness.
   **The `B8-intana` plane arm closed as P027 / T-157, and the P026 T999 verdict on this file was
   wrong**: the earlier "B8-intana fully consumed" note (item 6 above) only covered the cone/cone and
   cylinder families -- the `Perform(gp_Pln, gp_Pln)` arm of `IntAna_QuadQuadGeo` and the whole
   `IntAna_Int3Pln::Perform` were never diffed against the cxx. Three invented things came out:
   `plane_plane_line`'s least-norm origin formula `p = (d1*(n2 x v) + d2*(v x n1))/|v|^2` with a
   `v2 < 1e-20` cross-product gate (cxx has neither; `IntAna_QuadQuadGeo.cxx:423-445` is an
   `aEps = 1.e-16` denominator clamp plus two normal projections averaged), the whole
   `cxx:449-508` small-angle origin refinement (`aMVD < aTreshAng = 2.e-6` and
   `|aDist1| > aTreshDist = 1.e-12 || |aDist2| > aTreshDist` -> two
   `IntAna_IntConicQuad::Perform(gp_Lin, gp_Pln, TolAng, Tol)` solves), and `solve3`'s hand-rolled
   partial-pivot elimination with a `1e-14` pivot gate where `IntAna_Int3Pln::Perform` uses
   `math_Gauss Resol(M, gp::Resolution())` (`IntAna_Int3Pln.cxx:51`) -- plus a self-invented
   `Intersect3Pln::Line` collinear-plane arm (`v2 < 1e-20` + `|a*q + d| < 1e-9`) that OCCT does not
   have, since `IntAna_Int3Pln.hxx:29-70` is only point/empty. All four are now on the cxx path; the
   registered leftovers are the `Epsilon(radius)` tangency test in `quadric_quadric_plane_sphere`
   (`IntAna_QuadQuadGeo.cxx:1000`) and the sphere/sphere `0.01 * CONFUSION` circle test
   (`cxx:2034-2163`). Lesson recorded for the T999 walk: an inventory line for `IntAna_QuadQuadGeo`
   must be checked **per arm**, not per file.
   **P028 / T-158 closed both of those leftovers, and the per-arm lesson paid off**: the
   `Epsilon(radius)` test is now `(dist.abs() - radius).abs() < epsilon(radius)` (cxx:997;
   `nextafter(x, +inf) - x`, equal to `|x| * DBL_EPSILON` only on `|x| in [1,2)`, so a radius of 5
   or 10 was off by 4-8x), and the self-invented `dist_point_axis(...) > 1e-9` axis test of the
   cylinder/sphere and sphere/cone arms (cxx:1378-1379, cxx:1924-1925) plus its two copies in
   `occt-topo/src/intpatch_impimp_pairs.rs` became the new `axes_meet_at_point`, a line-for-line
   `AxeOperator` predicate (`Intersect() && Pt.Distance(PtIntersect()) == 0.0) || Same()`,
   `cxx:83-87`, defaults `1e-14` / `Precision::Angular()` at `cxx:72-73`; the `== 0.0` is OCCT's
   exact comparison and the port adds no tolerance either). The second leftover was closed by
   **verification, not edit**: the sphere/sphere `0.01 * CONFUSION` is `myEPSILON_MINI_CIRCLE_RADIUS`
   from `InitTolerances()` (`cxx:355`) and the arm matches the cxx line for line.
   The arm census this cycle produced is the durable part: `IntAna_QuadQuadGeo` has 15 `Perform`
   overloads, **10 ported** (`cxx:389` Pln/Pln, `:543` Pln/Cylinder, `:752` Pln/Cone, `:980`
   Pln/Sphere, `:1050` Cylinder/Cylinder, `:1324` Cylinder/Cone, `:1373` Cylinder/Sphere, `:1433`
   Cone/Cone, `:1917` Sphere/Cone, `:2034` Sphere/Sphere) and **5 absent, all Torus**
   (`cxx:2163` Pln/Torus, `:2278` Cylinder/Torus, `:2357` Cone/Torus, `:2496` Sphere/Torus, `:2588`
   Torus/Torus have no function in `occt-geom/src/intana`). Of the 10 ported only 5 have ever been
   line-diffed, and the 5 Torus arms must be decided from the consumer side -- reachable gap vs a
   marked-unported non-gap -- never with a substitute. Next pick: P029, the per-arm census of the
   five remaining ported arms plus the Torus verdict.
   **P029 / T-159 falsified that census, and the method lesson is the durable part**: the "5 absent,
   all Torus" claim above was produced by listing `crates/occt-geom/src/intana/` as a directory.
   That directory is not the whole class: `crates/occt-geom/src/intana/mod.rs:30-35` declares
   `#[path = "../intana_torus.rs"] mod torus;` followed by `pub use torus::*;`, so
   `crates/occt-geom/src/intana_torus.rs` -- one level up, outside the directory -- is part of the
   same namespace, and all five Torus arms live there: `quadric_quadric_plane_torus` (cxx:2163-2276,
   called from `occt-topo/src/intpatch_impimp_pairs.rs:117`), `quadric_quadric_cylinder_torus`
   (cxx:2278-2355, `:155`), `quadric_quadric_cone_torus` (cxx:2357-2494, `:162`),
   `quadric_quadric_sphere_torus` (cxx:2496-2586, `:169`, plus `occt-topo/src/intpatch.rs:576`) and
   `quadric_quadric_torus_torus` (cxx:2588-2666, `:175`). Corrected census: **15 of 15 `Perform`
   overloads present**, none absent. Rule for the T999 walk: enumerate by **`mod` declaration**,
   which may carry `#[path]` and therefore point outside the directory, not by directory listing.
   What the line-by-line walk of those five arms did find is one real invented construct, twice:
   `quadric_quadric_sphere_torus` (cxx:2496-2562) and `quadric_quadric_torus_torus` (cxx:2588-2666)
   each carried `if a_dist <= f64::EPSILON { if |r1-r2| > tol { Empty } else { Circles(one circle) } }`.
   The cxx has no such arm: after the `IntAna_Empty` early-out (cxx:2529 / cxx:2636) it computes
   `anAlpha = 0.5*(r1*r1 - r2*r2 + aDist*aDist)/aDist` (cxx:2534 / cxx:2640) and then constructs
   `gp_Dir aDir12(aVec12)` **unconditionally** (cxx:2542 / cxx:2647), which throws
   `Standard_ConstructionError` on a null vector. `f64::EPSILON` is `RealEpsilon()`, which those arms
   do not use (`myEPSILON_DISTANCE` cxx:353 and `myEPSILON_CYLINDER_DELTA_RADIUS` cxx:356 are the two
   they do). Both branches and the invented concentric-circle answer are deleted; the existing
   `GpDir::from_vec(..)` `Err` path stays as the stand-in for the OCCT exception with the cxx lines
   in a comment. Still never line-diffed after this cycle: cxx:543-745 Pln/Cylinder, cxx:752-972
   Pln/Cone, cxx:1050-1293 Cylinder/Cylinder, cxx:1324-1342 Cylinder/Cone, cxx:1433-1886 Cone/Cone.
   Next pick: P030, those five arms.
   **P030 / T-160 walked those five arms, so the class is now fully line-diffed** (2026-10-16,
   commit `0a185eb5`). Two gaps, one fixed and one marked:
   (a) `quadric_quadric_plane_cylinder` (`cxx:543-745`) had no `H` parameter at all -- its own
   comment admitted `(H = 0, so the near-parallel direction-refinement branch is skipped)` -- so the
   whole `cxx:568-598` refinement was declared away, and the parallel test was an inline
   `|axis . normal| < Tolang` instead of `IntAna_IntConicQuad(axec, P, tolang, toltang, H).IsParallel()`
   (`cxx:601`), which is not equivalent because the `IntAna_IntConicQuad` constructor also runs the
   `Len` end-point distance check when `H != 0`. Both are now ported line for line: `h: f64` (= `H`,
   `hxx:94` default 0), `cxx:568-598`, the faithful `IntAnaIntConicQuad::line_plane` delegate, the
   `newparams` recompute arms `cxx:617-631` (tangent single line) and `cxx:647-668` (two lines, with
   the `ht = max(arg, 0).sqrt()` clamp at `cxx:654-655`), and `inter.point(1)` at `cxx:687`. The
   consumer went with it: `occt-topo/src/intpatch_impimp_pairs.rs::int_p_cy` never passed a face
   height (so it always solved with `H = 0`); it now takes `VMax - VMin` of the cylinder from
   `uv1`/`uv2` per `IntPatch_ImpImpIntersection.cxx:2579-2587` (0 when either bound is infinite).
   (b) `quadric_quadric_plane_cone` (`cxx:752-972`) lacks the tail `cxx:922-951`: the
   `EllipseLimit = 1.0E+9` (`cxx:926,929`) and `HyperbolaLimit = 2.0E+6` (`cxx:935-950`) rejection,
   which sets `done = false` and returns. Every consumer turns that into `return false`
   (`IntPatch_ImpImpIntersection.cxx:3480-3483` here, i.e. this tree's `PairOutcome::Fail`), and
   branch 2 of `cone_cone_analytic` tests `INTER_QUAD_PLN.IsDone()` (`cxx:1556`). The port has no
   "not done" state for `QuadricIntersection`, and expressing it as `None` would mean
   `IntAna_Empty` -- the opposite -- so it is marked with the cxx line numbers only and the reason is
   in the `quadric_quadric_plane_cone` doc comment. Closing it needs a declined variant plus every
   exhaustive match updated (queued as P031).
   (c) The other three arms were verified faithful, no edit: `cxx:1050-1293` Cylinder/Cylinder,
   `cxx:1324-1342` Cylinder/Cone, `cxx:1433-1886` Cone/Cone -- P015's notes on cylinder/cylinder and
   cylinder/cone were quote checks only, this cycle supplied the line-level diff.
   Durable counts after this cycle: `IntAna_QuadQuadGeo` 15 `Perform` overloads, 15 present, 15
   line-diffed. One code-side item left in the class (the declined state, P031) plus the
   baseline-less STEP writer harness.
   Honest bound: no gate in this tree covers the new `H` arm (it needs a cylinder whose axis is
   slightly off the plane normal with a small face height), so the gates only prove no regression.
   Two analytic fast paths in `int_face_face_analytic.rs` still call with `H = 0` (commented).
   **P031 / T-161 closed that last code-side item** (2026-10-16, commit `9bef99ad`): the declined
   state is now `QuadricIntersection::NotDone` (= cxx `done = false`), and `quadric_quadric_plane_cone`
   translates `cxx:922-951` with `ELLIPSE_LIMIT = 1.0e9` (`cxx:927`) and `HYPERBOLA_LIMIT = 2.0e6`
   (`cxx:928`). The two cxx hyperbola blocks (`cxx:935-942` on param2/param2bis, `cxx:943-951` on
   param1/param1bis) collapse into one major/minor test here because both arms store
   `GpHypr::new(pos, param1 = param2, param1bis = param2bis)` with `nbint = 2` (`cxx:850-851`,
   `cxx:901-902`). Consumers: `from_quadric_intersection` maps it to `PairOutcome::Fail`
   (`IntPatch_ImpImpIntersection.cxx:3480-3483`), the two nested plane/cone reduces in `line_torus.rs`
   keep `None`, which is what cxx leaves `typeres` at when the nested call is not done
   (`cxx:1556`, `cxx:1835`), and the two self-made analytic fast paths return no curve so
   `int_face_face_analytic::plane_cone` falls back to `self.general`. Honest bound: no corpus model
   drives the tail (it needs a plane almost parallel to a generator with section parameters above
   1e9/2e6), so the four gates prove no regression only. Durable counts: `IntAna_QuadQuadGeo` 15
   overloads, 15 present, 15 line-diffed, zero code-side registrations. Next pick: P032, the
   `IntAna_IntLinTorus` done semantics (`IntAna_IntLinTorus.cxx:88-129`; `occt-geom.md` section 9),
   which needs a done-carrying return type for `line_torus_intersect` because the `quartic_roots`
   convenience wrapper drops `DirectPolynomialRoots::is_done()`.
   **P032 / T-162 closed that item** (2026-10-16, commit `55b0a6d6`): `line_torus_intersect`
   (`line_torus.rs:851-924`) now returns `LineTorusIntersection { done, points }` (`:826-843`), the
   `IsDone()` / `NbPoints()` / `Value(i)` shape of `IntAna_IntLinTorus.hxx:29-66` /
   `IntAna_IntLinTorus.lxx:18-44`. The body no longer goes through the `quartic_roots` adapter --
   `solution_values` (`math_direct_poly_roots.rs:703-711`) folds `!IsDone()`, `InfiniteRoots` and
   zero roots into the same empty `Vec<f64>`, while cxx gives opposite `done` for "solver declined"
   (`cxx:126-131`) and "solver finished with zero roots" (`cxx:121-125`). It holds
   `DirectPolynomialRoots::new5(a4, a3, a2, a1, a0)` instead, branches on `mdpr.IsDone()` (`cxx:89`)
   first, and walks the `NbSolutions()` / `Value(i)` list (`cxx:90-91`). Both `done = false` exits
   land: the declined solver (`cxx:126-131`) and the all-roots-rejected case
   `n > 0 && nbsolvalid == 0 && aNbBadSol == n` (`cxx:116-120`), the latter implemented as
   `n > 0 && out.is_empty() && nb_bad_sol == n` with the bad-root counter incremented at the
   `SquareDistance > 1e-10` back-substitution test (`cxx:99-106`). The consumer
   `crates/occt-topo/src/intcurvesurface/solvers.rs:456-465` returns `None` when `!is_done()`, i.e.
   `IntCurveSurface_InterUtils.pxx:1290-1294`'s `if (!intlintorus.IsDone()) { return false; }`, whose
   effect is `IntCurveSurface_Inter.pxx:578-581` clearing `isAnaProcessed` and sending the pair to the
   polyhedron interference. Honest bound: neither exit fires on the 15 exported or 23 parity models,
   and both gate logs are byte-identical to P031 (`.target-gate/t162/gates-raw.txt` vs
   `.target-gate/t161/gates-raw.txt`, and `export-obj.txt` vs `t156`), so this proves no regression
   only. Still registered: `occt-geom.md` section 9 item (b), the `out.sort_by(|a, b| a.x()...)` at
   `line_torus.rs:922`, which has no cxx counterpart (`cxx:104-108` writes `thePoint[nbsolvalid]` in
   root order, no reordering), plus the baseline-less STEP writer harness.
5. ~~**B4** `face_sample_point` -> edge-first + `PointInFace`. Touches BOP solids.~~
   **DONE 2026-10-15 (T-135, package P005).**

## F. Side findings (report only, not fixed)

- `~/.cursor/projects/d-source-repos-dogs/canvases/occt-code-review-status.canvas.tsx` SEMANTIC lines for topo-rest 1/4/6/11/18 and bop 1-4/7/8-10/16/17
  are stale in the "LIVE" direction: the code is already closed. The canvas needs a pass.
  P030 / T-160 did the tracking half (added the P030 cycle card, refreshed the gate-source comment,
  the `Stat` label and the "跟踪" paragraph, and replaced the stale `_index.md` "§9 -> 2026-10-11"
  pointer with the A0-A27 disposition rows plus `specs/_audit/occt-geom.md` §8); the SEMANTIC
  LIVE-staleness listed here is still open. P031 / T-161 then refreshed the same surface again for
  its own cycle (P031 card, `MEASURE.rust` 1,107 files / 338,040 lines, `CRATES` geom 44,144 /
  topo 214,139, the gate-source comment repointed to `.target-gate/t161`, `Stat` down to 1
  registered code item, and the "跟踪" paragraph plus the T999 callout extended to P027..P031);
  the SEMANTIC LIVE-staleness is still open, and P030's stale `40 文件/+2,071 行` delta figures were
  corrected to the measured `+47 文件/+2,297 行` against 10-13 in the same pass. P032 / T-162 did this pass a third time (P032 card, `MEASURE.rust`
  1,107 files / 338,101 lines, `CRATES` geom 44,191 / topo 214,153, the gate-source comment repointed
  to `.target-gate/t162`, the `Stat` moved to 2 registered code items -- section 9 item (b) and the
  baseline-less STEP writer harness -- and the "跟踪" paragraph plus the T999 callout extended to
  P027..P032); the SEMANTIC LIVE-staleness is still open. P033 / T-163 (`a4eeac65`) then did this
  pass a fourth time (P033 card, `MEASURE.rust` 1,107 files / 338,104 lines, `CRATES` geom 44,194,
  the gate-source comment repointed to `.target-gate/t163`, the `Stat` down to 1 registered code
  item -- the baseline-less STEP writer harness only -- and the "跟踪" paragraph plus the T999
  callout extended to P027..P033); the SEMANTIC LIVE-staleness is still open. The same pass removed
  the last `.rs`-side item in `occt-geom.md` section 9: the `out.sort_by(|a, b| a.x()...)` in
  `line_torus_intersect` (`line_torus.rs:922`) has no counterpart in `IntAna_IntLinTorus.cxx:104-108`
  (OCCT writes `thePoint[nbsolvalid]` in root-index order and never reorders), and the now-unused
  `Ordering` re-export in `intana/mod.rs` went with it. It also fixed four raw `<` characters and one
  `{ done, points }` JSX expression in the then in-repo board `specs/occt-scope-drift.canvas.tsx`
  (retired and deleted 2026-10-16; the board is now the Cursor canvas
  `~/.cursor/projects/d-source-repos-dogs/canvases/occt-code-review-status.canvas.tsx`), which had been breaking
  that board's build (`esbuild` reported "Expected identifier but found ="); the then-board transpiled
  clean and `tsc --noEmit` only reports the missing `dsh/canvas` / `react/jsx-runtime` modules.
  P035 / T-165 (`501fbdf7`) did the pass a sixth time the same day, with no measurement change (no
  files added, no `UNPORTED` line removed): its card, the gate-source comment repointed to
  `.target-gate/t165`, the `Stat` moved from 2 registered code items back to 1 (P034's two
  `CurveOnSurface` overrides minus the `Resolution` one P035 landed; the `gp_line` one plus the
  baseline-less STEP writer harness remain), and the "跟踪" paragraph plus the T999 callout
  extended to P027..P035.
  P036 / T-166 (`8d22a639`) did the pass a seventh time the same day, again with no measurement change
  (no files added, no `UNPORTED` line removed): its card, the gate-source comment repointed to
  `.target-gate/t166`, the `T0M` parity row corrected to the t166 log's `60200 / 66899` (the canvas
  had carried a stale `60199 / 66897`), the `Stat` still at 1 registered code item but now the
  pole-based plane classes instead of `gp_line` (P034's second registration is closed, and `gp_line`
  plus the baseline-less STEP writer harness are what remain), and the "跟踪" paragraph plus the T999
  callout extended to P027..P036.
  P037 / T-167 (`fd04dc8b`) did the pass an eighth time the same day, still with no measurement change
  (no files added, no `UNPORTED` count change -- the two in-source UNPORTED notes it removed were
  prose, not the token the count keys on): its card, the gate-source comment repointed to
  `.target-gate/t167`, the `Stat` moved from 1 registered code item to 0 (`Adaptor3d_CurveOnSurface`
  now has no code-side registration left: `Resolution` P035, elementary KPart P036, pole KPart P037;
  what remains is the baseline-less STEP writer harness, which is not a code item), and the "跟踪"
  paragraph plus the T999 callout extended to P027..P037. This pass also verified rather than re-pinned
  the `PARITY` readings: t167's `gates-raw.txt` differs from t166's only by the PowerShell
  `NativeCommandError` wrapper PowerShell put around the `Cube.step` line and the `finished in`
  duration, and its `export-obj.txt` is byte-identical (`Compare-Object` count 0), so the canvas
  `T0M` row stays at the t166 value.
  P034 / T-164 (`c272c08a`) did this measurement pass a fifth time (P034 card, `MEASURE.rust`
  1,111 files / 338,784 lines, `CRATES` geom 133/45,074 and topo 638/213,964, the gate-source
  comment repointed to `.target-gate/t164`, the `Stat` moved to 2 registered code items -- the two
  `CurveOnSurface` overrides P034 registered -- plus the baseline-less STEP writer harness, and the
  "跟踪" paragraph plus the T999 callout extended to P027..P034). That pass also pinned the two
  measurement clauses the earlier cycles had left implicit, because the numbers could not be
  reproduced from a naive count: `MEASURE.rust` is `src/**/*.rs` only (1111 files, matching
  `git ls-files "crates/*/src/*.rs"`; a whole-tree count gives 1138), and `MEASURE.unported` is the
  case-sensitive `UNPORTED` occurrence count over `crates/**/*.rs` (155, against 176
  case-insensitive) -- the previously recorded 1,107 reconciles exactly (1111 - P034's four new
  files), while the previously recorded 173 does **not** reconcile with either of the two counts
  reproducible today (155 case-sensitive / 176 case-insensitive) and its scope could not be
  recovered, so this pass records the two current readings with their scopes instead of carrying a
  number forward unverified.
- `specs/_audit/topo-rest.md` heading 17 (`pipe_along_polyline`), 9, 10, 15 likewise.
- ~~`bbox_from_geometry.rs` (B3) is a live sampled-bbox construct that no audit doc registers.~~
  **DONE 2026-10-15 (T-134)**: the sampling is gone; the item is now registered as
  `specs/_audit/topo-rest.md` §19.
- The port keeps two parallel implementations of `GeomAdaptor_Surface::U/VResolution`
  (`approx_same_parameter.rs` faithful, `wire_splitter_block.rs` hand-rolled); other call sites
  (`int_tools_wline.rs:400`, `brep_class/intersector_checks.rs:74`, `face_split.rs:97`) already
  delegate to the shared one.
- `validate.rs` (always-true `BRepCheck_Analyzer` stub) coexists with the real
  `brep_check.rs` `BRepCheck_Edge` port; `validate.rs` must not be used as a gate.
- `--lib` and `%TEMP%` write-permission failures are environment, to be distinguished from regressions
  in T002.
- The 3D faithful conversion entry is `occt_geom::convert_bspl::curve_to_bspline_curve` (it takes the
  `Parameterisation`) while the 2D one is `occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl`,
  so the OCCT name `CurveToBSplineCurve` is split across the two crates. Renaming the 2D one would
  touch 7 call sites in `occt-topo` (`comp_curve_to_bspline.rs:85,142`,
  `shhealing/shape_build_edge.rs:349,368`, `shhealing/wire_fix/pcurve_range.rs:63`,
  `iges/pcurve_transfer.rs:304,314`); registered by P022 / T-152, not changed.
- **NEW CLASS FOUND AND CLOSED BY P034 / T-164 (2026-10-16, `c272c08a`): the 3D point-to-curve
  *local* extrema family.** This walk had never named it. `crates/occt-geom/src/extrema_pc/` held
  only the global arm (`Extrema_GGExtPC`), so `Extrema_CurveTool`, `Extrema_GFuncExtPC`,
  `Extrema_GenLocateExtPC` and `Extrema_GLocateExtPC` had no counterpart at all, and the three sites
  that need "nearest point-curve extremum near a known parameter" each invented a clamped Newton:
  `occt-topo/src/int_tools_vertex_line.rs` (`LOCATE_ITERS = 12`, `|C'|^2 <= EXTREMA_EPS^2` gate,
  displacement stop) for `IntTools_Context.cxx:842-871`; `occt-topo/src/shhealing/wire_fix/
  same_parameter.rs` (12 iterations, `1e-20` speed gate, `1e-10` displacement stop) for
  `BRepLib_ValidateEdge.cxx:168-227`; `occt-topo/src/brep_lib_continuity.rs` (100 iterations plus a
  `has_analytic_d2()` precondition that made the whole edge silently skip) for
  `BRepLib.cxx:2341-2355`. All four classes are now ported and the three loops are deleted; the
  cxx control flow is a `math_FunctionRoot` whose result is re-checked against `|F| >= 1.e-07`
  (`Extrema_GenLocateExtPC.hxx:108-127`), run per `GeomAbs_C2` interval for free-form curves and
  replaced by an `Extrema_ExtPC` + nearest-to-`U0` selection for elementary ones
  (`Extrema_GLocateExtPC.hxx:125-307`). `extrema_ext_pc_range_done` was added because the elementary
  arm must test `IsDone()` before scanning the solutions (`Extrema_GGExtPC.hxx:462-469`). Lesson for
  the T999 re-walk: check the inventory **per class family** -- this one was missing while the walk
  was only reading `IntAna_*`/`BRepLib_*` lines.
- **Registered by P034, closed by P035 / T-165 (`501fbdf7`) and P036 / T-166 (`8d22a639`):**
  (a) `CurveOnSurface` had no `Adaptor3d_CurveOnSurface::Resolution` override
  (`Adaptor3d_CurveOnSurface.cxx:1364-1370`:
  `myCurve->Resolution(min(mySurface->UResolution(R3d), mySurface->VResolution(R3d)))`), so the cos
  `TolU` fell back to `Curve`'s `Precision::Parametric(Confusion)` (`occt-geom/src/curve.rs:135-137`).
  **Closed:** the override now lives on the shared meshing adaptor
  (`occt-topo/src/meshing/edge_discret.rs`) and reuses the two existing faithful ports -- the surface
  side `approx_same_parameter::{u_resolution,v_resolution}` (= `GeomAdaptor_Surface::U/VResolution`)
  and the curve side `pcurve.resolution` (= `Geom2dAdaptor_Curve::Resolution`) -- taking the min, so
  the consumer that changes is `extrema_pc::CurveTool::resolution` (i.e. the search-window criterion
  of P034's three call sites). Gates are word-identical to P034 (five crates 0 error =
  12/30/105/82/413 warnings, `occt-core --lib` 283/0, `occt-geom --lib` 143/0, `occt-geom2d --lib`
  72/0, `occt-topo --lib` 1251/0, `step_obj_gates` 4/4, `export_data_obj` 15/15; `.target-gate/t165/`),
  which shows the override does not move any reading on those models -- it does not by itself prove
  the new value is more correct elsewhere; the justification is that the cxx has the same override.
  (b) `CurveOnSurface` did not expose `gp_line()`, so a cos that `EvalKPart`
  classified as `GeomAbs_Line` (`Adaptor3d_CurveOnSurface.cxx:1565`, `:1687`, `:1701`) was typed
  `GeomAbs_OtherCurve` here (`Curve::gp_line` defaults to `None`, `occt-geom/src/curve.rs:43`).
  **Closed by P036, and P036 had to falsify P035's own registration to do it:** the note here in the
  P035 cycle said "both land in the same `Extrema_GLocateExtPC` free-form arm, so there is no
  observable difference today", which is wrong -- `Extrema_GLocateExtPC::Initialize`
  (`Extrema_GLocateExtPC.hxx:110-122`) branches on `GetType`, sending `GeomAbs_Line` to
  `myExtremPC` (elementary) and everything else to `myLocExtPC`, and `Perform` (`hxx:131-140`) then
  runs the free interval arm (`hxx:132-272`) or the elementary arm (`hxx:282-303`);
  `Extrema_GFuncExtPC::Initialize` (`hxx:115-128`) picks `MaxOrder + SearchOfTolerance` or `0 + MinTol`
  the same way. Since the adaptor's plane arm is `myType = CTy` (`cxx:1560`), five pcurve classes are
  reachable for a cos and each must answer its 3D query (`cxx:1381-1386`, `:1399-1402`, `:1406-1409`,
  `:1413-1416`); the port answered only `Circle()`. `eval_k_part` now returns a `KPart` and the plane
  arm became `plane_k_part` (the cxx `if/else` chain, reusing the existing `ElCLib::To3d` conic ports
  `occt-core/src/elib/clib.rs:633/651/642`), and `impl Curve for CurveOnSurface` forwards `gp_line`,
  `gp_ellipse`, `gp_parabola`, `gp_hyperbola`. `EvalD0`/`EvalD1`/`EvalD2` are untouched because the cxx
  evaluators special-case only Line and Circle.
  **Closed by P037 / T-167 (`fd04dc8b`):** the pole-based plane classes `Bezier()`
  (`cxx:1458-1485`) and `BSpline()` (`cxx:1489-1520`) would make `myType` reachable as
  `GeomAbs_BezierCurve`/`GeomAbs_BSplineCurve` (`Extrema_GGExtPC.hxx:190-243` is a real branch, not a
  `default:`); forwarding them makes `bezier_poles`/`bspline_poles`/`bspline_knots`/
  `bspline_weights`/`nurbs_degree` answer for a cos, and about 60 sites in the port use those
  accessors as class tests, so each needs its cxx criterion re-read first. The re-read is done:
  through `CurveOnSurface` those sites only change answer when `kpart_poles.is_some()`, i.e. only for
  a pole-classified pcurve on a **plane** surface, and the new answer matches what
  `Extrema_CurveTool` forwards in cxx. `plane_k_part` gained the `is_bezier2d() || is_bspline2d()` arm
  and reuses the one faithful implementation of both cxx bodies, `geom_lib::to_3d`
  (`GeomLib.cxx:582-632`); `KPart` carries the lifted curve (`edge_discret.rs:1025`) and the six
  accessors forward it (`:1390`/`:1394`/`:1398`/`:1402`/`:1406`/`:1411`). The second half of the gap
  was outside `occt-topo`: `Geom2dTrimmedCurve` forwarded only `trimmed_basis`, while
  `Geom2dAdaptor_Curve::load` (`Geom2dAdaptor_Curve.cxx:283-287`) unwraps a trimmed curve to its basis
  *before* the `DynamicType` chain, so `GeomAbs_BezierCurve` (`cxx:315-320`) and
  `GeomAbs_BSplineCurve` (`cxx:322-330`) resolve to the basis; 12 pole-class queries are now
  forwarded (`crates/occt-geom2d/src/trimmed.rs:188-224`). Gates word-identical to P036
  (`.target-gate/t167/`). `GeomAbs_OffsetCurve` needs no work: all three consumers put it in their
  `default:` arm. **`Adaptor3d_CurveOnSurface` now has no code-side registration left** (`Resolution`
  P035, elementary KPart P036, pole KPart P037).
  P038 / T-168 (`ff235d37`) did the measurement pass a ninth time the same day, still with no
  measurement change (no files added, no `UNPORTED` line removed): its card, the gate-source comment
  repointed to `.target-gate/t168`, the `Stat` still at 0 registered code items (P038 closed an item
  that lived only in the `bop-boolean.md` findings list, not on the canvas `OPEN` list), and the
  "跟踪" paragraph plus the T999 callout extended to P027..P038. P038 left the `CurveOnSurface` line
  for a different production line (BOP): it deleted the `bop_builder::boolean` disjoint-bbox
  short-circuit, the last live invented item in the BOP family, and recorded findings 4/7 as moot
  because `crates/occt-topo/src/bop_curved/` was deleted 2026-10-14. This pass verified rather than
  re-pinned the `PARITY` readings again: t168's `gates-raw.txt` differs from t167's by exactly two
  lines (the `finished in` durations), and the 15 `ok` lines of its `export-obj.txt` are
  word-identical to t167's, so the canvas `T0M` row stays at the t166 value. Side finding left
  unfixed: `tests/phase3_integration.rs::bbox_and_mesh_agree` fails on a `1.0000001` vs `1.0`
  bbox comparison; stashing confirms it predates P038.
  P039 / T-169 (`46ab36ec`) did the measurement pass a tenth time the same day, again with no
  measurement change (no files added to `src/**`, no `UNPORTED` line removed net): its card, the
  gate-source comments repointed to `.target-gate/t039`, the `Stat` still at 0 registered code
  items, and the "跟踪" paragraph plus the T999 callout extended to P027..P039. P039 took the
  third candidate line the P038 card had left on the table -- `sweep_revolve`'s lateral-face
  structure vs `BRepPrim_OneAxis` -- and **closed its OneAxis half**: the wire order and occurrence
  orientations now follow `BRepPrim_OneAxis::LateralWire` (`BRepPrim_OneAxis.cxx:666-679` plus
  `BRepPrim_Builder::AddWireEdge`, `BRepPrim_Builder.cxx:184-192`), a full revolution reuses one
  seam TShape twice (`StartEdge()` aliases `EndEdge()` when `!HasSides()`, `cxx:331-334`,
  `cxx:976-979`, `:1028-1031`), the file finally registers pcurves (`LateralFace`, `cxx:388-439`;
  `TopFace`/`BottomFace`, `cxx:465-468`, `:504-509`), and ring vertices carry
  `Precision::Confusion()` like `BRepPrim_Builder::MakeVertex` (`BRepPrim_Builder.cxx:129-132`).
  **The measurement falsified the attribution, not just the fix**: switching
  `brepfeat::revolve_profile_about` to the analytic solid still breaks
  `bop_builder::boolean`, so the line "the gap is in `sweep_revolve`'s wire/seam/pcurve
  conventions" is now recorded as *insufficient*, and the lathe-mesh tool stays. The next slice
  must look at the `BOPAlgo` side. This pass verified rather than re-pinned `PARITY` once more:
  t039's `gates-raw.txt` is the same 4/4 with the same 23 ratios, and its 15 `ok` lines are
  word-identical to t168's, so the canvas `T0M` row is unchanged.
  P040 / T-170 (`229a667a`) did the measurement pass an eleventh time the same day, again with no
  measurement change (no files added to `src/**`, no `UNPORTED` line removed net): its card, the
  gate-source comments repointed to `.target-gate/t040`, the `Stat` still at 0 registered code
  items, and the "跟踪" paragraph plus the T999 callout extended to P027..P040. P040 took the
  line the P039 card handed over -- "locate on the `BOPAlgo` side why the analytic revolve tool
  cannot be combined" -- and **narrowed it to the BOP split step** with a control group instead of
  another attribution guess. Tool side first: `brep_sweep::make_revol` (= `BRepSweep_Revol`) for
  the groove ring profile measures 4 face / 6 edge / 4 vertex, volume `3.342654583419527` =
  pi(1.3^2-0.6^2)*0.8, signed volume positive, so the P039 guess "the horizontal cap faces do not
  follow `BRepPrim_OneAxis::BottomFace`/`TopFace`" is falsified. Two `brep_sweep` control-flow gaps
  were landed on the same pass: `BRep_Tool::Parameter(V, E)` (`BRep_Tool.cxx:1340-1397`; the
  parameter is chosen by the orientation of the vertex **as stored in the edge**, found via
  `TopExp_Iterator`/`IsSame`, with `rev` from the edge orientation, not by the caller's handle
  orientation) and `BRepSweep_Rotation::DirectSolid`'s UV box (`BRepSweep_Rotation.cxx:635-676`
  constructs `BRepAdaptor_Surface(TopoDS::Face(aGenS))` whose constructor defaults to
  `Restriction = true` and runs `BRepTools::UVBounds`, `BRepAdaptor_Surface.cxx:57-81`; a plane's
  natural range is infinite, so sampling it puts the midpoint off the face and flips the solid's
  orientation). Then the measurement: `boolean` on the analytic cylinder (r=1, h=3) `Cut` the
  analytic ring returns 6 faces with `result shell is not closed` and volume `6.936636579126224`
  (= 3pi - pi*0.792) -- the two wall pieces (2pi, 2pi*1.2) and one annulus (pi*0.64) are right, the
  tool's inner cylinder r=0.6 is missing and the other annulus measures `4.272566008882088`; the
  same Cut through the lathe tool gives a closed 10-face shell that removes only `0.0519` (should
  be 1.6085); the faceted mesh target cut by the analytic ring removes only `0.3715`. The control
  group "analytic cylinder `Cut` box (x,y in +-1.3, z in -1..1.8, so no coplanar overlap at z=0)"
  returns `3.7699111843077335` = pi*1.2 exactly, 3 faces, closed -- so the Cut path itself is fine.
  Conclusion: the mismatch tracks **a whole-period (seamed) face that must be split**, not whether
  the target is an analytic cylinder or the tool has exact cylindrical/planar faces; the existing
  baseline is green only because both operands are pre-faceted (the lathe tool cuts each ring into
  `steps` arcs, `faceted_cyl` is itself a mesh). One more guess was excluded: the seam branch at
  `bop_split_faces_occt.rs:145-153` (`BOPTools_AlgoTools3D::DoSplitSEAMOnFace`,
  `BOPAlgo_Builder_2.cxx:429-446`) does **not** raise
  `BOPAlgo_AlertUnableToMakeClosedEdgeOnFace` for the failing pair (a temporary probe on that
  failure path produced no output), so the source-edge seam doubling is wired and the gap sits
  later, in the Loops/Areas construction of a periodic face cut by a closed section curve. Gates
  word-identical to P039 (`.target-gate/t040/`).
