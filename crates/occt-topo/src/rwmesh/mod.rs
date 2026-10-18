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


mod coordinate_system;
mod mesh_name;
mod scene;
mod vrml;
mod gltf;
mod gltf_caf_writer;
mod gltf_document;
mod gltf_extras;
mod gltf_file;
mod gltf_late_data;
mod gltf_material;
mod gltf_prim_array;
mod gltf_shape_reader;
mod gltf_triangulation_reader;
mod gltf_types;
mod obj_caf_writer;
mod ply_caf_writer;
pub use coordinate_system::*;
pub use mesh_name::*;
pub use scene::*;
pub use vrml::*;
pub use gltf::*;
pub use gltf_caf_writer::{
    bin_file_path, write_gltf as perform_gltf, GltfCafNode, GltfCafWriterOptions,
    GltfCafWriterOutput, GltfWriterTrsfFormat,
};
pub use gltf_extras::*;
pub use gltf_file::*;
pub use gltf_material::*;
pub use gltf_prim_array::*;
pub use gltf_shape_reader::*;
pub use obj_caf_writer::*;
pub use ply_caf_writer::{
    cxx_g, cxx_g6, write_ply as perform_ply, PlyCafWriterOptions, PlyWriterContext,
};

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
