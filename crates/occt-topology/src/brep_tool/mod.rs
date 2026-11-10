//! Geometric access to topology — a port of `BRep_Tool`.
//!
//! Reads real geometry (curves, surfaces, tolerances) back from the
//! `GeometryRegistry` side-table that `TopoBuilder` populates. Unregistered
//! shapes fall back to OCCT-like defaults (origin point, unbounded ranges,
//! no curve/surface).

mod accessors;
mod closedness;
mod edge_face_queries;
mod helpers;
mod integrity;
mod uv_parameter;
mod vertex_parameter;
mod world;

/// Namespace for BRep_Tool-style static queries.
pub struct BRepTool;

#[cfg(test)]
mod tests;
