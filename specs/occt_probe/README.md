# OCCT ground-truth probe (T-82 / T-67 / T-05 acceptance)

The port's acceptance for every OCCT-alignment task in this session is a numeric
comparison against **OCCT 8.0.0 itself**. This directory keeps the small C++
probe used for that, so the instrument survives a fresh checkout (the working
copy lives under `.target-gate/`, which is git-ignored).

## Build & run

```
specs\occt_probe\build_run.bat <step-file> [mode]
```

* `build.bat` — `vcvars64.bat` + `cl /std:c++17 /EHsc /MD /MP` against
  `D:\source\occt-8.0.0\inc` and the `win64\vc14\lib` import libraries.
  `TKGeomAlgo.lib` is required for `GeomAPI_ProjectPointOnSurf`.

## Source layout

The probe is one executable built from several translation units (the monolithic
`occt_probe.cpp` exceeded the line budget and was split by command group in T-193):

| file | contents |
|---|---|
| `occt_probe.cpp` | `main`: argument parsing, file read, the `run*` dispatch loop |
| `occt_probe_common.h` | the OCCT/STL include set, `ProbeCtx`, cross-file declarations |
| `modes_*.cpp` | the `argv[2]` command groups (`runXxxModes`, returning `std::optional<int>`) |
| `delaun_probe.cpp` | `--delauncheck` / `--delaunfeed` / `--delaunstruct` / `--delaunsteps` / `--delaunwhy` |
| `boundary_probe.cpp` | `--boundary` |
| `structdump_probe.cpp` | `--lens`, `--structdump`, `--dcstruct`, `--srfnodes` (`ZZ*Algo` subclasses) |
| `srfparams_probe.cpp` | `--srparams`, `--facestats`, `--facedump`, `--edgedu` |
| `meshdbg_probe.cpp` | `--torusdbg`, `--healerdbg`, `--ctxdbg`, `--rangereplay`, `--algospy` |
| `wiredbg_probe.cpp` | `--orderdbg`, `--wiredbg`, `--wiredbgi` |
| `iges_probe.cpp` | `--iges`, `--igespc` |
| `extrema_probe.cpp` | `--offsetsurf`, `--ssx` |
| `modes_default.cpp` | the fallback `BRepCheck` / `BRepGProp` / per-face report |

Splitting changed no control flow: the dispatch blocks were moved verbatim into
their group function and are tried in the original order.
* `run_probe.bat` — the probe is `/MD`, so the exe is copied next to the OCCT
  DLLs before it is started.
* `probe.bat` — `THIRDPARTY_DIR` **must be absolute** (`env.bat`'s default is
  relative, and a missing `tbb12.dll`/`jemalloc.dll` shows up as `0xC0000135`).

## Port-side counterparts

T-195 removed the port's ad-hoc diagnostic writers (`DELAUNDUMP`,
`DELAUNSTEPS`/`DELONESTEP`, `DELSTAGE`/`DELTRIS`/`DCWHY`, `CTRLDUMP`, `PORTDUMP`,
`SORTDUMP`, `CWDBG`/`CBWIN`, `PC3D`/`SURFD2`/`PCUV2`, the model-specific
`PCBTRACE`/`RANGEW`/`GDBG` probes) so that library code no longer writes dump
files. Any probe mode whose input came from one of those writers is **archival**:
it still runs and still prints the OCCT half of the comparison, but the
`delaun_f*.txt` / `delaun_in_f*.txt` / per-step files it reads (or the port output
it is diffed against) now have to come from a run archive, not from the current
port. Live port diagnostics are the flags exported by
`crates/occt-topo/src/debug_flags.rs` (`PCW`, `PCDEF`, `PCRAW`, `PC2D`, `NURBSGRID`,
`BPUV`, `DELREGSKEL`, `IGES_TRACE*`, ...); everything else in the tables above
describes OCCT behavior only.

## Modes

| mode | prints |
|---|---|
| (none) | read a BRep/STEP and report per-solid `BRepGProp` volume |
| `--cylall` | every face of `BRepPrimAPI_MakeCylinder(0.4, 2.0)` with its `BRepAdaptor_Surface` type, face orientation, and per-edge `CurveOnSurface` orientation + UV mid |
| `--faceoff` | `BOPAlgo_Builder` fuse of `box[-1,1]^3` and the same cylinder; for the section circle it prints every adjacent face and the face `BOPTools_AlgoTools::GetFaceOff` picks from each |
| `--dir` | the same GF: per adjacent face, `BOPTools_AlgoTools3D::GetNormalToFaceOnEdge` (`aDN`), `BOPTools_AlgoTools2D::EdgeTangent` (`aDTgt`) and their cross product (`aDB`) |
| `--dir2` | faithful re-implementation of the file-local `FindPointInFace` (`BOPTools_AlgoTools.cxx:2160-2231`) and `GetFaceDir` (`:2110-2152`, including the `GetApproxNormalToFaceOnEdge` fallback at `:2139-2149`) and `MinStep3D` (`:2235-2346`), so the *final* `aDB`/angle of every candidate can be read |
| `--ssx` | `Extrema_ExtSS` / `Extrema_GenExtSS` / `GeomAPI_ExtremaSurfaceSurface` on `sphere(r=1)@origin` x `sphere(r=1)@(4,0,0)`, `sphere` x `plane z=5`, and `plane z=5` x `sphere` (no input file is read; pass any placeholder as argv[1]). Prints the adaptor windows, then `IsDone`/`IsParallel`/`NbExt` and every `SquareDistance`/`PointOnS1`/`PointOnS2` parameter and point, the 12-argument windowed `Extrema_ExtSS`, the `GeomAPI` `LowerDistance`/`NearestPoints`, and the raw `Extrema_GenExtSS(..., 20, 20, ...)` list. This is the ground truth for B6d (`occt_geom::extrema_ss`) |
| `--uv <faceIdx> [--uvdef <defl> <angle>]` | build the discrete model (`BRepMesh_ModelBuilder`) and discretize edges (`BRepMesh_EdgeDiscret`) at the given deflection/angle (default `0.1` / `0.349066`), then print the face's surface range/periodicity and, per wire edge, the `IPCurve` orientation / `IsForward` / point count / first-last UV / parameter range plus every sample as `UV   p<k> uv=(u,v) t=<param>`. This is the GT a port pcurve dump (`PCW=1`) is diffed against. `--uvsum` alone prints the per-face pcurve summary |
| `--pcraw <first> [<tol>]` | `StepToGeom::MakeCurve2d` on every `StepGeom_Pcurve` entity of the STEP model whose built curve starts at `<first>`, i.e. the **constructor** output before `XSAlgo_ShapeProcessor::CheckPCurve` can replace it: entity numbers, `pf`/`pl`, degree, poles (with weights when rational), distinct knots + multiplicities, and the multiplicity-expanded `flat` knot list. Port counterpart: `PCRAW=<first>` (`step_to_geom_curve2d.rs::dump_pcraw`) |
| `--pcdef <first> [<tol>]` | the same definition read off the **imported shape** (`BRep_Tool::CurveOnSurface` over every face/wire/edge), i.e. after the `CheckPCurve` reprojection, plus the `[0]..[4]` and fixed-`t` sample UVs. Port counterpart: `PCDEF=1` |
| `--pcf <faceIdx>` | per-edge pcurve summary of one face in `TopExp_Explorer` order (`c3f`/`c3l`/`pf`/`pl`/type/poles/knots/degree + 5 sample UVs). Port counterpart: `PCDEF` face blocks |
| `--spcheck <first> [<tol>] [<preci>]` | replay of `XSAlgo_ShapeProcessor::CheckPCurve` (`cxx:344-505`): temporary-edge tolerance after `ShapeFix_Edge::FixSameParameter`, whether OCCT re-projects, the reprojected type/knots/degree, `uv_dev` against the source pcurve, and the `take` verdict. Port counterpart: `shhealing::pcurve_ranges::check_pcurve` |
| `--fbox <x1> <y1> <z1> <x2> <y2> <z2>` | find the face whose bbox matches (0.5 tol) and print its surface type / `ShapeAnalysis_Surface` singularities / U,V closure, then every wire edge's vertices, orientation, degeneracy and tolerances |
| `--edges [defl] [angle]` | `BRepMesh_IncrementalMesh` then, per face occurrence, every edge's `PolygonOnTriangulation` node count and degeneracy flag |
| `--edgerep` | per-edge `BRep_TEdge` representation census, the exact input of `BRepCheck_Edge::Minimum` (`BRepCheck_Edge.cxx:88-148`): `EDGEREP n= deg= spar= srange= c3drep= c3dnull= cos=`. `deg` is `BRep_Tool::Degenerated`, `c3drep`/`c3dnull` say whether a `BRep_Curve3D` representation is present and whether its curve handle is null, `cos` counts the `CurveOnSurface` representations. Port counterpart: `record_dump --edgerep` (same format, plus the edge range); the failing `BRepCheck_Status` list is `record_dump --statuses` |
| `--statusdump` | per-sub-shape `BRepCheck_Analyzer` verdict: `VALID 0|1` then one `STATUS type= n= min=<list> [ctx=<type>:<list> ...]` line per sub-shape, where `min` is `BRepCheck_Result::Status()` and every `ctx=` is `StatusOnShape()` for one contextual shape. This is the OCCT half of `record_dump --statuses` and the only way to read *which* `BRepCheck_Status` and *which* sub-shape make a shape invalid |
| `--record <out.txt> [defl] [angle]` | RECORD v1 dump of the shape `main` read (`modes_record.cpp`): root type, `TopExp::MapShapes` census, `BRepCheck_Analyzer` verdict, `BRepGProp` volume/area/centroid, exact bbox, one sorted `FACE` line per face (surface type, area, wires, edges, bbox), and the `BRepMesh_IncrementalMesh` node/triangle totals and 3D triangle area |
| `--record-rt <iges\|step> <mid> <out.txt> [defl] [angle]` | same record after a full write/read round trip with OCCT's own writer and reader (`IGESControl_Writer`/`IGESControl_Reader`, `STEPControl_Writer`/`STEPControl_Reader`) |

### RECORD parity (source-level, no external models)

`--record` / `--record-rt` exist for the source-level parity gate: both sides run
the same pipeline from the same STEP source file, so nothing is compared against
a model produced by an external application.

```
# OCCT half (this probe, linked against the OCCT 8.0.0 sources)
specs\occt_probe\probe.bat <abs>\data\Cube.step --record <out>\cube.occ.txt 0.1 0.5

# port half (same pipeline, port reader/writer)
cargo build --manifest-path crates\occt-topo\Cargo.toml --offline --example record_dump
crates\occt-topo\target\debug\examples\record_dump.exe data\Cube.step <out>\cube.port.txt --defl 0.1 --angle 0.5
```

Both halves for every `data/*.step`, in `read` and `rt-iges` modes, plus the
field diff (structural vs numeric) and the per-record diff files:

```
powershell -NoProfile -File specs\record_parity\run_record_parity.ps1 [-Only Cube,rev]
```

Records and `*.diff.txt` land in `.target-gate/record/`. Fields are graded in four
bands, because the two ways a record can differ are different work:

| band | meaning | fails the gate |
|---|---|---|
| `STRUCT` | a decision differs: `TYPE`, `VALID`, `COUNTS`, face type/wire/edge counts, volume sign | yes |
| `REAL` | geometry differs above the noise floor (> 1e-6 relative and > 1e-4 absolute) | yes |
| `THRESH` | only `MESH.nv`/`MESH.nf` differ | no |
| `TINY` | floating fields below the noise floor | no |

`--record`'s exactness band is `1e-7` absolute plus `1e-9` relative; `TINY` and
`REAL` split everything above it. `THRESH` is its own band because node/triangle
totals are the output of `BRepMesh_DelaunayDeflectionControlMeshAlgo`'s
`deviation > deflection` split: a link that lands on the threshold flips on
summation order alone, and each flipped node adds exactly two triangles (measured
on `data/Shape.step`: `dnf == 2*dnv` on every affected face). Integer equality
there describes the float path, not the port, so it does not fail the gate. The
exit code is the number of `STRUCT` + `REAL` records.

The `--faceoff` / `--dir` / `--dir2` outputs are the acceptance assertions for
task T-82 (see `specs/_tasks.md`).
