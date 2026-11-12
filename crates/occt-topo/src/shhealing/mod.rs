//! Shape healing — repair free edges, remove small edges, close non-closed
//! wires, weld coincident vertices and relocate vertices.
//!
//! Port of the OCCT `ShapeHealing` toolkit: `ShapeFix_Wire` (FixSmall,
//! FixClosed, FixVertexTolerance), `ShapeFix_Shape`, `ShapeFix_Edge`,
//! `ShapeFix_Vertex`, and `ShapeAnalysis_FreeBounds`.
//!
//! Every fix follows the same pipeline: collect the real topology from the
//! `TShape` child tree, build maps from the old sub-shapes to their healed
//! replacements, then rebuild a fresh `TopoShape` tree that reuses the shared
//! `TShape` handles (preserving vertex/edge sharing between faces) and
//! re-registers geometry in the `GeometryRegistry` side-table for any newly
//! created shapes. Healing never mutates the input shape — it always returns a
//! new shape, so the caller can keep the original intact.
//!
//! Geometry (`BRep_TVertex`/`BRep_TEdge`/`BRep_TFace`) lives in the
//! process-wide `GeometryRegistry` keyed by `TShape` address. Rebuilding an
//! edge or vertex therefore creates a *new* `TShape` and registers its geometry
//! afresh; edges and vertices that survive unchanged are reused by handle, so
//! their existing side-table entries remain valid.
mod prelude {

    pub(crate) use std::collections::{HashMap, HashSet};
    pub(crate) use std::sync::Arc;

    pub(crate) use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
    pub(crate) use occt_geom::GeomLine;

    pub(crate) use crate::abs::ShapeType;
    #[cfg(test)]
    pub(crate) use crate::brep_measure::edge_length;
    pub(crate) use crate::builder::TopoBuilder;
    pub(crate) use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
    pub(crate) use crate::tgeometry::GeometryRegistry;
    pub(crate) use crate::topo_tools::{
        edge_vertices, edges_of, edges_of_wire, faces_of, is_same, vertex_position, vertices_of,
        wires_of, wires_of_face,
    };
}

mod box_bnd_tree;
mod common;
mod face_geom_helpers;
mod face_natural_bound;
mod face_split;
mod free_bounds;
mod intersection_tool;
mod loop_wire;
mod pcurve_ranges;
mod shape_build_edge;
mod shape_fix_edge;
mod shape_fix_face;
mod shape_fix_face_wires;
mod shape_fix_parameters;
mod shape_fix_same_parameter;
mod shape_fix_shape;
mod shape_fix_shape_modes;
mod shape_fix_shell;
mod shape_fix_vertex_position;
mod shape_upgrade_unify_same_domain;
mod small_area;
mod split_tool;
mod wire_fix;
mod wire_heal;
mod wire_self_inter;
mod wire_self_inter_driver;
mod wire_self_inter_fix;
mod xsalgo_check_pcurve;
pub use common::*;
pub use face_split::split_face_after_removed_segment;
pub use free_bounds::*;
pub use intersection_tool::*;
pub use loop_wire::*;
pub use pcurve_ranges::*;
pub use shape_build_edge::*;
pub use shape_fix_edge::*;
pub use shape_fix_face::*;
pub use shape_fix_parameters::*;
pub use shape_fix_same_parameter::*;
pub use shape_fix_shape::*;
pub use shape_fix_shape_modes::*;
pub use shape_fix_shell::*;
pub use shape_fix_vertex_position::*;
pub(crate) use shape_upgrade_unify_same_domain::UnifySameDomain;
pub use small_area::*;
pub use split_tool::*;
pub use wire_fix::*;
pub use wire_heal::*;
pub use wire_self_inter::*;
pub use wire_self_inter_driver::*;
pub use wire_self_inter_fix::*;
pub use xsalgo_check_pcurve::*;

/// `ShapeAnalysis_Curve` projection helpers (`ShapeAnalysis_Curve.cxx`).
pub mod shape_analysis_curve;
/// `ShapeAnalysis_TransferParameters(Proj)` (`ShapeAnalysis_TransferParameters.cxx`,
/// `ShapeAnalysis_TransferParametersProj.cxx`).
pub mod transfer_params;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
