# T001 receipt — live drift/gap inventory (2026-10-15)

Board: `docs/goals/drift-and-gap-autofix`.
Method: four read-only Scout passes over `specs/occt-scope-drift.canvas.tsx` (SEMANTIC array),
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

### B6. `occt-geom` Extrema failure fabrication (live)
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
    ratios identical to the T000 baseline (22 x 1.0000, Shape-2 0.9999), `export_data_obj` 15/15
    with v/f and bbox bit-identical.
- **B6e OPEN (found while closing B6c, 2026-10-15).** `extrema2d/curve_curve.rs::locate_extcc2d`
  (and its 3D twin `occt-geom/src/extrema_cc/curve_curve.rs::locate_extcc`) use a self-invented
  clamped 32-iteration Newton (`refine_curve_curve2d_newton` + `f_and_jac2d`) with a 1e-12
  convergence test. The cxx is `Extrema_LocateExtCC2d` (`Extrema_LocateExtCC2d.cxx:28-46`) ->
  `Extrema_LocECC2d` = `Extrema_GenLocateExtCC` (`Extrema_GenLocateExtCC.hxx:44-105`): build
  `Extrema_CCLocFOfLocECC2d` = `Extrema_GFuncExtCC` (`Extrema_CCLocFOfLocECC2d.hxx:27-34`), run
  `math_FunctionSetRoot` with `Tol = (TolU, TolV)` = each curve's
  `Extrema_Curve2dTool::Resolution(Confusion)` (`Extrema_LocateExtCC2d.cxx:32-36`) inside
  `[First,Last]^2`, then require `S.IsDone() && F.NbExt() > 0` and `Ff(i) < 1e-10` on both
  components. None of 32 / 1e-12 / clamping exist in the cxx; `math_FunctionSetRoot` and the
  `Extrema_GFuncExtCC` shape are already available (`occt-math`, and the same pattern is landed in
  `extrema_surf/extrusion_point_extrema.rs`). Consumers: `extrema2d/tests.rs`,
  `extrema_cc/tests.rs` only.
- **B6f OPEN (found while closing B6c, 2026-10-15).**
  `extrema2d/curve_curve.rs::curve_curve_intersections2d` is a 33x33 parameter grid plus a `tol`
  argument. OCCT has no such sampler: the zero-distance case is the engine already wired
  (`Extrema_ExtCC2d`, used by `curve_ops::curve2d_intersections`) or
  `Geom2dAPI_InterCurveCurve` -> `Geom2dInt_GInter` (ported in `occt-geom2d/geom2d_int`).
  Consumer: `extrema2d/tests.rs:47` only.
- **B6g OPEN (duplication found while closing B6c, 2026-10-15).** `extrema2d/curve_curve.rs`'s
  `line_line_extrema2d`, `line_circle_extrema2d` and `circle_circle_extrema2d` are a second,
  module-level copy of the `Extrema_ExtElC2d` arms that `extrema2d/ext_el_c2d.rs` already carries
  (and that the engine now uses); nothing outside `extrema2d/tests.rs` calls them.
- **B6d OPEN (found 2026-10-15, T-137).** `extrema_ss::surface_surface_extrema_all`'s analytic
  sphere/sphere and sphere/plane arms have no counterpart in `Extrema_ExtSS::Perform`, whose only
  analytic arm is `Plane` x `Plane` -> `Extrema_ExtElSS::Perform(S1.Plane(), myS2->Plane())`
  (`Extrema_ExtSS.cxx:122-129`); every other type pair goes to `Extrema_GenExtSS` (`cxx:131-180`).
  This is an extra branch, not a sampled substitute; Scout must decide between wiring the
  non-plane pairs to the `GenExtSS` equivalent (`surface_surface_newton_all`) and keeping them.

### B7. Step/IGES-adjacent live divergences
- `mesh_to_brep.rs:39-41,118-119,157-158`: 1e-9 vertex dedup grid, vertex tolerance 0.0,
  per-triangle planar faces; OCCT `BRepBuilderAPI_MakeShapeOnMesh` is node 1:1, `BRepLib::Precision()`
  = 1e-7, and returns a Compound. Live via `brepfeat/features.rs:367`.
- `step/write_context.rs:314`: periodic B-spline not `SetNotPeriodic()` before writing
  (`GeomToStep_MakeBoundedCurve.cxx:46-51`).
- `step/format.rs:500`: KnotSpec always `.UNSPECIFIED.`
  (`MakeBSplineSurfaceWithKnots.cxx:127-150`).
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
- `surface_fit.rs:13-40` `fit_plane` lacks the `gz < Tol` gate (`GeomLib_IsPlanarSurface.cxx:44-45`);
  reachable via `brep_offset/curve_face_offset.rs:489`. (Re-verified 2026-10-15: `fit_plane` is a
  least-squares fit, not a port of `GeomLib_IsPlanarSurface`, so the faithful fix is to port
  `Controle`/`GeomLib::Inertia` and wire the caller, not to bolt a gate onto the existing fit.)
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

- **C-A (leaf, deps present, reachable)** — the B1 items above, plus `geom2d_convert.rs:77`,
  `step/write_context.rs:314`, `step/format.rs:500` from B7; and `extrema_cc/mod.rs:9`,
  `globoptmin.rs:590` from B8.
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
  `line_line`/`line_circle`/`circle_circle_extrema2d` copies of `Extrema_ExtElC2d`).
  **B6d** (still open, found in T-137): the analytic sphere/sphere and sphere/plane arms of
  `surface_surface_extrema_all` have no counterpart in `Extrema_ExtSS::Perform` (only
  `Plane` x `Plane` is analytic there).
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
   `math_FunctionRoots`. Remaining B6 items, in no fixed order: **B6e** (`locate_extcc2d` Newton
   vs `Extrema_GenLocateExtCC` + `math_FunctionSetRoot`; test-only consumers), **B6d** (the extra
   analytic sphere arms in `surface_surface_extrema_all`), **B6f**/**B6g** (the 33x33 intersection
   grid; the duplicated `Extrema_ExtElC2d` arm copies). **The next slice is #1 below (B1)**: it
   outranks the B6 leftovers on centrality, confidence and blast radius.
4. ~~**B5** `find_nearest_valid_point` `aD1Mag` arm. Small, self-contained, live pave filler.~~
   **DONE 2026-10-15 (T-133).**
6. **B8 osculating_surface early returns.** **DONE 2026-10-15 as T-139 / package P009** (four
   `clear_flags(); return;` -> OCCT's `clearOsculFlags()` + continue, plus the per-arm `build_chain`
   local). Next picks, in order: **B8-intana** (`intana` quartic/direct-polynomial thresholds, live via
   `intcurvesurface/solvers.rs:456`), **B7-writer** (`step/write_context.rs:314` `SetNotPeriodic`,
   `step/format.rs:500` `KnotSpec`, `geom2d_convert.rs:77`), **B7-mesh_to_brep**
   (`BRepBuilderAPI_MakeShapeOnMesh`), then the B6 leftovers (**B6e**/**B6f**/**B6g**/**B6d**) and the
   `B8 surface_fit` / `B6e` Newton items.
5. ~~**B4** `face_sample_point` -> edge-first + `PointInFace`. Touches BOP solids.~~
   **DONE 2026-10-15 (T-135, package P005).**

## F. Side findings (report only, not fixed)

- `specs/occt-scope-drift.canvas.tsx` SEMANTIC lines for topo-rest 1/4/6/11/18 and bop 1-4/7/8-10/16/17
  are stale in the "LIVE" direction: the code is already closed. The canvas needs a pass.
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
