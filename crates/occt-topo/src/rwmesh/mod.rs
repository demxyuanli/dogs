//! RWMesh — mesh format import (OBJ/PLY/STL/glTF/VRML) and scene assembly.
//! Source: `RWMesh`.
//!
//! Reads triangle meshes from Wavefront OBJ (with MTL materials and UVs), PLY,
//! STL (binary/ASCII auto-detected), glTF 2.0 (external or embedded base64
//! buffers, node TRS/matrix transforms, PBR base-color materials) and VRML 2.0
//! (Transform nesting, IndexedFaceSet, diffuse materials), assembles them into a
//! `MeshScene`, and converts scenes back to B-Rep shapes via `mesh_to_brep`. A
//! tiny format router (`convert_mesh_format`) round-trips between the formats by
//! extension.
mod prelude {

pub(crate) use std::collections::HashMap;
pub(crate) use std::path::Path;

pub(crate) use occt_core::gp::{GpMat, GpPnt, GpPnt2d, GpQuaternion, GpTrsf, GpXyz, TrsfForm};
pub(crate) use occt_core::io::ply::PlyMesh;
pub(crate) use occt_core::io::stl::StlMesh;
pub(crate) use occt_core::poly::Triangulation;
pub(crate) use occt_core::poly::triangulation::Triangle;

pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::gltf::GltfOptions;
pub(crate) use crate::mesh_to_brep::triangulation_to_brep;
pub(crate) use crate::shape::{Compound, TopoShape};

}


mod scene;
mod vrml;
mod gltf;
pub use scene::*;
pub use vrml::*;
pub use gltf::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
