//! Port of OCCT BRepMesh range splitters — Wave 4 BRepMesh.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_DefaultRangeSplitter.{hxx,cxx}` — base splitter: computes the
//!   discrete UV range of a face from its geometric range and the registered
//!   boundary points, plus a face-basis scale (`delta`).
//! - `BRepMesh_{Cylinder,Cone,Sphere,Torus,NURBS}RangeSplitter.{hxx,cxx}` —
//!   analytical splitters generating interior (surface) nodes.
//! - `BRepMesh_UVParamRangeSplitter.hxx`, `BRepMesh_UndefinedRangeSplitter`,
//!   `BRepMesh_BoundaryParamsRangeSplitter`, `BRepMesh_ExtrusionRangeSplitter`.
//!
//! The OCCT classes form an inheritance tree rooted at `BRepMesh_DefaultRangeSplitter`.
//! Rust models that with a [`RangeSplitter`] trait; every concrete splitter embeds
//! its parent as an `inner` field and exposes it through [`RangeSplitter::base`] /
//! [`RangeSplitter::base_mut`], while the virtual hooks (`compute_delta`,
//! `generate_surface_nodes`, `get_undefined_interval_nb`, …) are trait methods the
//! derived types override.
//!
//! `ponytail:` surfaces arrive as `Arc<dyn Surface>` which cannot be downcast, so
//! the analytic type is classified from the periodic flags + `d0` sampling, and
//! radii (`gp_Cylinder::Radius()`, `gp_Torus::MajorRadius()`, …) are measured by
//! sampling diametrically opposite iso-parameter points instead of reading the
//! concrete `gp_*` handles. The NURBS interval machinery reads knots/poles only
//! through `dyn Surface` (`NbUIntervals` / `NbUPoles`).
mod prelude {


pub(crate) use std::f64::consts::PI;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
pub(crate) use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL, SQUARE_CONFUSION};
pub(crate) use occt_geom::Surface;

pub(crate) use super::super::data_model::*;
pub(crate) use super::super::parameters::MeshParameters;

pub(crate) use crate::brep_tool::BRepTool;

}


mod param_set;
mod splitter;
mod nodes;
pub use param_set::*;
pub use splitter::*;
pub use nodes::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
