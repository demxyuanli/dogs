# Loop 6 Plan — 5000 line existing module coverage expansion

BSplCLib completion:
  - knots_advanced.rs: knot removal, degree reduction, periodization (200 lines)
  - rational.rs: rational curve ops, weight management (180 lines)
  - curve_tools.rs: parameterization, reparameterization (150 lines)

BSplSLib completion:
  - surface_advanced.rs: surface weight ops, trim, Coons patches (250 lines)
  - surface_knots.rs: u/v knot management, seam handling (180 lines)

Math completion:
  - eigenvalue.rs: eigenvalue searcher, power iteration (250 lines)
  - trig_roots.rs: trigonometric polynomial roots (200 lines)
  - multi_integrate.rs: Gauss multiple integration (200 lines)
  - linear_utils.rs: linear system utilities, condition number (150 lines)

Bnd completion:
  - obb_pca.rs: PCA-based OBB computation (150 lines)
  - sortbox2d.rs: 2D sorted box (120 lines)
  - intersect.rs: box intersection utilities (100 lines)

Poly completion:
  - polygon3d_full.rs: full 3D polygon with normals (150 lines)
  - coherent.rs: coherent triangulation (120 lines)
  - array_tools.rs: polygon array utilities (100 lines)

Kernel additions:
  - units.rs: unit conversion system (180 lines)
  - stdfail.rs: standard exception types (80 lines)
  - resource.rs: resource manager stub (120 lines)

BVH completion:
  - traversal.rs: BVH traversal + intersection (180 lines)
  - builder_triangle.rs: triangle-based BVH builder (120 lines)

TopLoc completion:  
  - datum.rs: multi-level datum chain (150 lines)
  - compound.rs: compound location (100 lines)

Convert additions:
  - advanced.rs: additional coordinate transforms (120 lines)
  - projection.rs: map projections (100 lines)

TOTAL: ~4,100 lines estimated
