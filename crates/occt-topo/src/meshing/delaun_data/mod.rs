//! Port of OCCT `BRepMesh_DataStructureOfDelaun` + `BRepMesh_SelectorOfDataStructureOfDelaun`.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_DataStructureOfDelaun.{hxx,cxx}` — indexed storage for vertices,
//!   links (edges) and triangles plus a vertex-cell filter for coincidence
//!   queries (the OCCT class holds this as a `BRepMesh_VertexTool`).
//! - `BRepMesh_SelectorOfDataStructureOfDelaun.{hxx,cxx}` — neighbor selector.
//!
//! ponytail: the OCCT value types carry movability inside the struct
//! (`BRepMesh_Vertex/Edge/Triangle`). The Rust contract types in `delaun_types`
//! (`DelaunVertex`/`DelaunLink`/`DelaunTriangle`) do not, so the data structure
//! keeps parallel per-element state arrays here. When the contract is extended,
//! these arrays can be folded into the value types.
//!
//! Indices follow OCCT: vertices/links/elements are **1-based** in the public
//! API; storage is 0-based `Vec`s.

mod auxiliary;
mod elements;
mod links;
mod nodes;
mod selector;
mod structure;
mod vertex_cells;

pub use selector::DelaunSelector;
pub use structure::DelaunDataStructure;

#[cfg(test)]
mod tests;
