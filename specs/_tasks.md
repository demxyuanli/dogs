# _tasks — OCCT → Rust systematic completion map

> Last updated: 2026-07-31. In-progress loop starting from TKernel upwards.

## Legend
```
🟢 done    🟡 in-progress    ⚪ pending    ⬜ deferred
```

## Phase 1: FoundationClasses/TKernel — root infrastructure

| # | Unit | Status | Effort | Notes |
|---|---|---|---|---|
| TK-01 | `TCollection` | 🟡 | S | Typed collection aliases (Vec<T>, HashMap<K,V>) |
| TK-02 | `TColStd` | 🟡 | S | Standard collection aliases (TColStd_Array1OfReal→Vec<f64>) |
| TK-03 | `GeomAbs` | 🟡 | S | Geometry enums (CurveType, SurfaceType, Shape, JoinType) |
| TK-04 | `Quantity` | ⚪ | S | Physical quantity types (Color, Length, Angle, etc.) |
| TK-05 | `Units/UnitsAPI` | ⚪ | M | Unit conversion (mm→m, rad→deg, etc.) |
| TK-06 | `Message` | ⚪ | M | Progress indicator + messaging (use Rust tracing/log) |
| TK-07 | `OSD` | ⚪ | L | OS-dependent layer (files, threads, signals → stdlib) |
| TK-08 | `Resource` | ⚪ | M | Resource file parser (.res format) |
| TK-09 | `Storage` | ⚪ | L | Binary persistence (OCC format reader/writer) |
| TK-10 | `FSD` | ⚪ | M | File system driver abstraction |
| TK-11 | `StdFail` | ⚪ | S | Standard exception subclasses |
| TK-12 | `Plugin` | ⚪ | M | Dynamic library loading |

## Phase 2: FoundationClasses/TKMath — algorithms + data structures

| # | Unit | Status | Effort | Notes |
|---|---|---|---|---|
| TM-01 | `gp` (all 37 types) | 🟢 done | — | Complete |
| TM-02 | `precision` | 🟢 done | — | Complete |
| TM-03 | `elib` (ElCLib+ElSLib) | 🟢 done | — | Complete |
| TM-04 | `bnd` (4/7) | 🟢 done | S-M | BoundSortBox, B2, B3, Tools remaining |
| TM-05 | `math` (4/60+) | 🟡 | L | SVD, Crout, Jacobi done. Need BFGS, Brent, Newton, PSO, Gauss, Kronrod |
| TM-06 | `Poly` | ⚪ | M | Polygon2D, Polygon3D, Triangulation data + algorithms |
| TM-07 | `BVH` | ⚪ | M | Bounding volume hierarchy builder |
| TM-08 | `TopLoc` | ⚪ | M | Topology location (nested datum transforms) |
| TM-09 | `GeomAbs` | 🟡 | S | (see TK-03) |
| TM-10 | `Convert` | ⚪ | M | Coordinate conversion (polar→cartesian, etc.) |
| TM-11 | `CSLib` | ⚪ | M | Classifier for parametric surfaces |
| TM-12 | `PLib` | ⚪ | L | Polynomial evaluation (Horner, Lagrange, etc.) |
| TM-13 | `BSplCLib` | ⚪ | XL | B-spline curve core (knot insertion, evaluation, etc.) |
| TM-14 | `BSplSLib` | ⚪ | XL | B-spline surface core |

## Phase 3: ModelingData — geometry + topology

| # | Unit | Status | Effort | Notes |
|---|---|---|---|---|
| MD-01 | `TKG2d` curves | 🟢 done | — | 5 curve types + Curve2d trait |
| MD-02 | `TKG2d` remaining | ⚪ | L | TrimmedCurve, OffsetCurve, BSplineCurve, BezierCurve |
| MD-03 | `TKG3d` curves | 🟢 done | — | GeomLine/Circle/Ellipse/Hyperbola/Parabola |
| MD-04 | `TKG3d` surfaces | 🟢 done | — | GeomPlane/Cylinder/Cone/Sphere/Torus |
| MD-05 | `TKG3d` remaining | ⚪ | L | TrimmedCurve, OffsetCurve, BSplineCurve, SurfaceOfRevolution, etc. |
| MD-06 | `TKGeomBase` | ⚪ | XL | Geometric utilities (projections, extrema, intersections) |
| MD-07 | `TKBRep` | ⚪ | XL | Boundary representation data structures |

## Phase 4+: Deferred
| # | Layer | Status | Notes |
|---|---|---|---|
| MA-* | ModelingAlgorithms | ⬜ | 13 toolkits, gated by G1 (invariant oracle) |
| DE-* | DataExchange | ⬜ | 13 toolkits, gated by complete ModelingData |
| VI-* | Visualization | ⬜ | 6 toolkits |
| AF-* | ApplicationFramework | ⬜ | 12 toolkits |

## Execution order
```
🟡 TK-01→TK-02→TK-03   (collections + enums, 1 session)
   ↓
🟡 TM-04→TM-05→TM-06   (bnd completions + math solvers + poly)
   ↓
⚪ TM-07→TM-08→TM-10   (BVH + TopLoc + Convert)
   ↓
⚪ TK-04→TK-05→TK-06   (Quantity + Units + Message)
   ↓
⚪ TM-11→TM-12          (CSLib + PLib)
   ↓
⬜ TM-13→TM-14          (BSplCLib + BSplSLib, XL effort, defer if possible)
   ↓
⬜ MD-02→MD-05→MD-06→MD-07 (ModelingData completion)
```
