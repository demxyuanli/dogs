//! General geometric utilities — polyline sampling, curve discretization.
//!
//! **UNPORTED (audit A14)**: this package is **not** an OCCT translation. It
//! holds port-local helpers (polyline sampling/discretization, fitting, planar
//! polygon utilities, a Bowyer–Watson 2D triangulation and a voxel CSG) whose
//! OCCT counterparts live elsewhere: `GCPnts_*`/`BRepMesh_*` (discretization and
//! face meshing), `GeomAPI_PointsToBSpline`/`math_*` (fitting),
//! `IntTools`/`BOPTools_AlgoTools` (planar predicates) and `BRepAlgoAPI_*` +
//! `BOPAlgo_Builder` (booleans — OCCT has **no** voxel boolean, see `csg.rs`).
pub mod polyline;
pub mod triangulate;
pub mod fit;
pub mod csg;
pub mod fit2;
pub mod polygon_ops;
pub mod delaunay;
pub mod curve_ops3d;
pub mod polyline_simplify;
pub mod mesh_analysis;
pub mod polygon_boolean;
pub mod curve_frenet;
pub mod curve_interp3d;
