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

The `--faceoff` / `--dir` / `--dir2` outputs are the acceptance assertions for
task T-82 (`specs/_board.md` §3.1b/§3.2).
