# _tasks — OCCT → Rust migration atomic breakdown

> Updated: 2026-07-31. 28 units indexed. 14 done/verified.
> Rules: [specs/_rules.md](_rules.md). Machine index: `migration-index.json`.

## Layer 0 — roots (no deps)

| # | Atom | Status | Source | Region | Effort |
|---|---|---|---|---|---|
| 1 | `gp_xyz` | ✅ done | `TKMath/gp/gp_XYZ.cxx` | wide | S |
| 2 | `gp_pnt` | ✅ done | `TKMath/gp/gp_Pnt.cxx` | wide | S |
| 3 | `gp_vec` | ✅ done | `TKMath/gp/gp_Vec.cxx` | wide | S |
| 4 | `gp_dir` | ✅ done | `TKMath/gp/gp_Dir.cxx` | wide | M |
| 5 | `gp_ax1` | ✅ done | `TKMath/gp/gp_Ax1.cxx` | wide | S |
| 6 | `gp_ax2` | ✅ done | `TKMath/gp/gp_Ax2.cxx` | wide | M |
| 7 | `gp_ax3` | ✅ done | `TKMath/gp/gp_Ax3.cxx` | wide | M |
| 8 | `gp_mat` | ✅ done | `TKMath/gp/gp_Mat.cxx` | wide | S |
| 9 | `gp_trsf` | ✅ done | `TKMath/gp/gp_Trsf.cxx` | wide | M |
| 10 | `gp_quaternion` | ✅ done | `TKMath/gp/gp_Quaternion.cxx` | wide | M |
| 11 | `gp_circ` | ✅ done | `TKMath/gp/gp_Circ.cxx` | wide | S |
| 12 | `gp_pln` | ✅ done | `TKMath/gp/gp_Pln.cxx` | wide | S |
| 13 | `gp_cylinder` | ✅ done | `TKMath/gp/gp_Cylinder.cxx` | wide | S |
| 14 | `gp_cone` | ✅ done | `TKMath/gp/gp_Cone.cxx` | wide | S |
| 15 | `gp_sphere` | ✅ done | `TKMath/gp/gp_Sphere.cxx` | wide | S |
| 16 | `gp_torus` | ✅ done | `TKMath/gp/gp_Torus.cxx` | wide | M |
| 17 | `gp_trsf_form` | ✅ done | `TKMath/gp/gp_TrsfForm.hxx` | repl | S |
| 18 | `precision` | ✅ done | `TKernel/Precision/Precision.hxx` | wide | S |
| 19 | `standard_types` | ✅ done | `TKernel/Standard/Standard_TypeDef.hxx` | repl | S |
| 20 | `gp_xy` | ✅ done | `TKMath/gp/gp_XY.hxx` | wide | S |
| 21 | `gp_mat2d` | ✅ done | `TKMath/gp/gp_Mat2d.hxx` | wide | S |

## Layer 1 — depends on layer 0

| # | Atom | Status | Source | Deps | Region | Effort |
|---|---|---|---|---|---|---|
| 22 | `gp_pnt2d` | ✅ done | `TKMath/gp/gp_Pnt2d.cxx` | gp_xy | wide | S |
| 23 | `gp_vec2d` | ✅ done | `TKMath/gp/gp_Vec2d.cxx` | gp_xy | wide | S |
| 24 | `gp_dir2d` | ✅ done | `TKMath/gp/gp_Dir2d.cxx` | gp_xy | wide | S |
| 25 | `gp_ax2d` | ✅ done | `TKMath/gp/gp_Ax2d.cxx` | gp_pnt2d,gp_dir2d | wide | S |
| 26 | `gp_ax22d` | ✅ done | `TKMath/gp/gp_Ax22d.cxx` | gp_pnt2d,gp_dir2d | wide | S |
| 27 | `gp_trsf2d` | ✅ done | `TKMath/gp/gp_Trsf2d.cxx` | gp_xy,gp_mat2d | wide | M |
| 28 | `math_vector` | ✅ done | `TKMath/math/math_VectorBase.hxx` | — | trans | S |
| 29 | `math_matrix` | ✅ done | `TKMath/math/math_Matrix.hxx` | math_vector | trans | S |
| 30 | `math_intvec` | ✅ done | `TKMath/math/math_IntegerVector.hxx` | — | repl | S |
| 31 | `math_status` | ✅ done | `TKMath/math/math_Status.hxx` | — | repl | S |

## Layer 1 — ready to port (no deps on incomplete units)

| # | Atom | Source | Region | Effort | Notes |
|---|---|---|---|---|---|
| 32 | `standard_transient` | `TKernel/Standard/Standard_Transient.*` | repl | M | Arc<T> wrapper |
| 33 | `standard_handle` | `TKernel/Standard/Standard_Handle.hxx` | repl | M | Handle protocol, gate G2 |
| 34 | `standard_failure` | `TKernel/Standard/Standard_Failure.*` | repl | M | Exception hierarchy |
| 35 | `ncollection` | `TKernel/NCollection/NCollection_*` | repl | L | Container templates |

## Layer 2 — ready to port (math solvers, elementary curves)

| # | Atom | Source | Deps | Effort | Notes |
|---|---|---|---|---|---|
| 36 | `math_gauss` | `TKMath/math/math_Gauss.*` | math_vector,math_matrix | M | Linear system solver, already in MathMatrix::solve |
| 37 | `math_svd` | `TKMath/math/math_SVD.*` | math_vector,math_matrix | M | Singular value decomposition |
| 38 | `math_crout` | `TKMath/math/math_Crout.*` | math_vector,math_matrix | M | LU decomposition |
| 39 | `math_householder` | `TKMath/math/math_Householder.*` | math_vector,math_matrix | M | QR decomposition |
| 40 | `math_jacobi` | `TKMath/math/math_Jacobi.*` | math_vector,math_matrix | M | Eigenvalue solver |
| 41 | `el_clib` | `TKMath/ElCLib/ElCLib.cxx` | gp_* | L | Curve evaluation (lines,circles,ellipses,etc.) |
| 42 | `el_slib` | `TKMath/ElSLib/ElSLib.cxx` | gp_* | L | Surface evaluation (planes,cylinders,cones,etc.) |

## Layers 3-7 — deferred

| Layer | Module | Toolkits | Status |
|---|---|---|---|
| 3 | ModelingData | TKG2d,TKG3d,TKGeomBase,TKBRep | ○ blocked by handle protocol (G2) |
| 4 | ModelingAlgorithms | TKBO,TKBool,TKFeat,TKFillet,TKOffset,TKMesh,TKShHealing,TKTopAlgo,TKXMesh | ○ blocked by invariant oracle (G1) |
| 5 | DataExchange | STEP,IGES,STL,glTF,VRML,OBJ,PLY | ○ deferred |
| 6 | Visualization | OpenGL,3D viewer,selection | ○ deferred |
| 7 | ApplicationFramework | OCAF document framework | ○ deferred |
