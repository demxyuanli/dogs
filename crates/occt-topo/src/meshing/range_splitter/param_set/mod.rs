//! Parametric parameter set, surface classification, NURBS grid generation and
//! the base range-splitter trait.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_DefaultRangeSplitter.{hxx,cxx}` — base splitter.
//! - `BRepMesh_NURBSRangeSplitter.{hxx,cxx}` — NURBS grid generation
//!   (`getUndefinedInterval`, `computeGrainAndFilterParameters`, `filterParameters`,
//!   `AnalyticalFilter`, `GenerateSurfaceNodes`).
//! - `GeomAdaptor_Surface` / `GeomAdaptor_Curve` interval queries.

mod analytical;
mod intervals;
mod nurbs_grid;
mod set;
mod splitter;
mod surface_type;

pub use self::analytical::*;
pub use self::intervals::*;
pub use self::nurbs_grid::*;
pub use self::set::*;
pub use self::splitter::*;
pub use self::surface_type::*;
