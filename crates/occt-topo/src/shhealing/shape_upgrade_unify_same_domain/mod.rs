//! `ShapeUpgrade_UnifySameDomain`.
//!
//! Source: OCCT 8.0.0
//! `ModelingAlgorithms/TKShHealing/ShapeUpgrade/ShapeUpgrade_UnifySameDomain.cxx`.
//!
#![allow(dead_code)]

mod boundaries;
mod boundary_rebuild;
mod context;
mod curve_helpers;
mod edge_chain;
mod edge_merge;
mod face_helpers;
mod history;
mod pcurves;
mod periodic_faces;
mod same_domain;
mod same_surface;
mod unify_edges;
mod unify_faces;
mod union_pcurves;
mod uv_bounds;
mod wire_helpers;

pub(crate) use unify_faces::UnifySameDomain;
