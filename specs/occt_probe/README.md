# OCCT ground-truth probe (T-82 / T-67 / T-05 acceptance)

The port's acceptance for every OCCT-alignment task in this session is a numeric
comparison against **OCCT 8.0.0 itself**. This directory keeps the small C++
probe used for that, so the instrument survives a fresh checkout (the working
copy lives under `.target-gate/`, which is git-ignored).

## Build & run

```
specs\occt_probe\build_run.bat <step-file> [mode]
```

* `build.bat` — `vcvars64.bat` + `cl /std:c++17 /EHsc /MD` against
  `D:\source\occt-8.0.0\inc` and the `win64\vc14\lib` import libraries.
  `TKGeomAlgo.lib` is required for `GeomAPI_ProjectPointOnSurf`.
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

The `--faceoff` / `--dir` / `--dir2` outputs are the acceptance assertions for
task T-82 (see `specs/board.canvas.tsx`).
