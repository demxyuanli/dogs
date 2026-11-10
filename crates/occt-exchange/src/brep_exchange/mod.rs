//! Module 2: brep_exchange — export BRep shapes to OBJ / STL / PLY.
//! Source: `RWObj`, `RWStl`, `RWPLY` (simplified OCCT writers).

mod caf_nodes;
mod gltf;
mod mesh_export;
mod model_export;
mod style;

pub use caf_nodes::*;
pub use gltf::*;
pub use mesh_export::*;
pub use model_export::*;

#[cfg(test)]
mod tests;
