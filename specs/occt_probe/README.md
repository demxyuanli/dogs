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
| `--fbox <x1> <y1> <z1> <x2> <y2> <z2>` | find the face whose bbox matches (0.5 tol) and print its surface type / `ShapeAnalysis_Surface` singularities / U,V closure, then every wire edge's vertices, orientation, degeneracy and tolerances |
| `--edges [defl] [angle]` | `BRepMesh_IncrementalMesh` then, per face occurrence, every edge's `PolygonOnTriangulation` node count and degeneracy flag |

The `--faceoff` / `--dir` / `--dir2` outputs are the acceptance assertions for
task T-82 (see `specs/_tasks.md`).
