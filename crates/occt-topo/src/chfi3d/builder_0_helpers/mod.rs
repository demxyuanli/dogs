//! Helpers of `ChFi3d_Builder_0.cxx` (OCCT 8.0.0).
//!
//! `ChFi3d_Builder_0.cxx` carries the file statics the `ChFi3d_Builder`
//! members call; this module hosts the ones the `StartSol` chain and the
//! `PerformSetOfKPart` / `PerformSetOfKGen` walk need.

pub mod inter_plane_edge;
pub mod k_particular;
pub mod perform_el_spine;
pub mod solid_index;

pub use inter_plane_edge::inter_plane_edge;
pub use k_particular::ch_fi3d_k_particular;
pub use perform_el_spine::ch_fi3d_perform_el_spine;
pub use solid_index::ch_fi3d_solid_index;
