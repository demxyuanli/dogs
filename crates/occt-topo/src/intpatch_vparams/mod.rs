//! `ComputeVertexParameters` of `IntPatch_WLine`, `IntPatch_RLine`,
//! `IntPatch_GLine` and `IntPatch_ALine`.
//!
//! Source:
//! - `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_WLine.cxx:88-321`
//!   (statics), `:323-1211` (`IntPatch_WLine::ComputeVertexParameters`).
//! - `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_RLine.cxx:143-433`.
//! - `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_GLine.cxx:414-1102`.
//! - `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_ALine.cxx:77-676`.

pub mod aline;
pub mod common;
pub mod gline;
pub mod rline;
pub mod wline;

pub use aline::compute_vertex_parameters_aline;
pub use gline::compute_vertex_parameters_gline;
pub use rline::compute_vertex_parameters_rline;
pub use wline::compute_vertex_parameters as compute_vertex_parameters_wline;
