//! OCCT hidden-line removal - pure Rust port of Open CASCADE Technology TKHlr.
//!
//! Sits directly above the modelling core: every module here consumes
//! `occt-topo` shapes and produces projected 2D output.
//!
//! # Architecture
//! ```text
//! hlr_algo/          (projector, interference and intersection records)
//! hlrbrep*/          (HLRBRep data, curves, surfaces, hider, algorithms)
//! hlrtopobrep/       (HLRTopoBRep stage that prepares topological data)
//! contap/            (Contap_Contour - contour/outline extraction used by HLR)
//! ```

pub mod contap;
pub mod hlr;
pub mod hlr_algo;
pub mod hlrbrep;
pub mod hlrbrep_algo;
pub mod hlrbrep_area_limit;
pub mod hlrbrep_c_inter;
pub mod hlrbrep_curve;
pub mod hlrbrep_curve_2d;
pub mod hlrbrep_curve_tool;
pub mod hlrbrep_data;
pub mod hlrbrep_edge_builder;
pub mod hlrbrep_edge_data;
pub mod hlrbrep_edge_face_tool;
pub mod hlrbrep_edge_i_list;
pub mod hlrbrep_edge_interference_tool;
pub mod hlrbrep_face_data;
pub mod hlrbrep_face_iterator;
pub mod hlrbrep_hider;
pub mod hlrbrep_hlr_to_shape;
pub mod hlrbrep_inter_csurf;
pub mod hlrbrep_internal_algo;
pub mod hlrbrep_intersector;
pub mod hlrbrep_line_tool;
pub mod hlrbrep_poly_algo;
pub mod hlrbrep_poly_hlr_to_shape;
pub mod hlrbrep_shape_bounds;
pub mod hlrbrep_shape_to_hlr;
pub mod hlrbrep_sl_props;
pub mod hlrbrep_surface;
pub mod hlrbrep_surface_tool;
pub mod hlrbrep_type_of_resulting_edge;
pub mod hlrbrep_vertex_list;
pub mod hlrtopobrep;
