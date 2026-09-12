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
pub(crate) use crate::brep_measure::edge_length;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topo_tools_full::{
    edges_of, edges_of_wire, edge_vertices, is_same, vertex_position, vertices_of, wires_of,
};

}


mod p01;
mod p02;
mod p03;
mod p04;
pub use p01::*;
pub use p02::*;
pub use p03::*;
pub use p04::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
