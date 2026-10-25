//! Helpers of `ChFi3d_Builder_0.cxx` (OCCT 8.0.0).
//!
//! `ChFi3d_Builder_0.cxx` carries the file statics the `ChFi3d_Builder`
//! members call; this module hosts the ones the `StartSol` chain and the
//! `PerformSetOfKPart` / `PerformSetOfKGen` walk need.

pub mod boite;
pub mod bound_surf;
pub mod compute_arete;
pub mod compute_curves;
pub mod compute_curves_walking;
pub mod compute_pcurv;
pub mod curve_cleaner;
pub mod enlarge_box;
pub mod eval_tol_reached;
pub mod int_traces;
pub mod inter_plane_edge;
pub mod k_particular;
pub mod mkbound;
pub mod perform_el_spine;
pub mod point_in_ds;
pub mod project_pcurv;
pub mod reparam_pcurv;
pub mod solid_index;

pub use boite::{ch_fi3d_boite, ch_fi3d_boite4};
pub use bound_surf::ch_fi3d_bound_surf;
pub use compute_arete::{ch_fi3d_compute_arete, ComputeArete};
pub use compute_curves::{ch_fi3d_compute_curves, ComputeCurves};
pub use compute_pcurv::{ch_fi3d_compute_pcurv, ch_fi3d_compute_pcurv_uv};
pub use curve_cleaner::curve_cleaner;
pub use enlarge_box::{
    ch_fi3d_enlarge_box_curve, ch_fi3d_enlarge_box_edge, ch_fi3d_enlarge_box_pcurve,
    ch_fi3d_enlarge_box_surf_data,
};
pub use eval_tol_reached::ch_fi3d_eval_tol_reached;
pub use int_traces::ch_fi3d_int_traces;
pub use inter_plane_edge::inter_plane_edge;
pub use k_particular::ch_fi3d_k_particular;
pub use mkbound::{
    build_p_curve_dirs, build_p_curve_vec2d, build_p_curve_vec3d, mkbound_curve,
    mkbound_curve_segment, mkbound_pfac2d, mkbound_pfac3d, mkbound_segment,
    mkbound_segment_adaptor,
};
pub use perform_el_spine::ch_fi3d_perform_el_spine;
pub use point_in_ds::{
    ch_fi3d_index_of_surf_data, ch_fi3d_index_point_in_ds, ch_fi3d_set_point_tolerance,
};
pub use project_pcurv::ch_fi3d_project_pcurv;
pub use reparam_pcurv::ch_fi3d_reparam_pcurv;
pub use solid_index::ch_fi3d_solid_index;
