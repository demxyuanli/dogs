//! Port of OCCT incremental_mesh — BRepMesh integration.
//!
//! `BRepMesh_IncrementalMesh` is the meshing entry point: it builds the discrete
//! model (`ModelBuilder`), pre-processes it (`ModelPreProcessor`), then runs the
//! OCCT-style per-face pipeline — `EdgeDiscret` collects the boundary UV points
//! from the discretized edges, `Triangulator` triangulates them into 3D
//! triangles with deflection-controlled interior refinement — and assembles the
//! crate's `ShapeMesh` triangle soup.
//!
//! Reference: `BRepMesh_IncrementalMesh.hxx/.cxx`, `BRepMesh_DiscretRoot.hxx`.
//!
//! # Pipeline
//! `ModelBuilder::build_model` → `ModelPreProcessor::perform` →
//! per-face [`IncrementalMesh::triangulate_model_faces`]
//! (`EdgeDiscret` boundary UV → `Triangulator::triangulate`) → assembled
//! `ShapeMesh`.
//!
//! # Fallback
//! If any pipeline stage fails for a shape (unusual/degenerate geometry), the
//! pre-pipeline wireframe UV-grid tessellator
//! ([`IncrementalMesh::build_shape_mesh_wireframe`]) is used so `perform` never
//! panics.
mod prelude {

pub(crate) use std::collections::HashMap;

pub(crate) use occt_core::gp::{GpPnt, GpPnt2d};
pub(crate) use occt_core::poly::triangulation::Triangle;
pub(crate) use occt_geom::Surface;

pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::intpatch::refine_point_on_surface;
pub(crate) use crate::mesh::{mesh_surface_area, ShapeMesh};
pub(crate) use crate::shape::{Edge, Face, TopoShape};
pub(crate) use crate::wireframe;

pub(crate) use super::super::data_model::{MeshEdge, MeshModel, MeshStatus};
pub(crate) use super::super::delaun_types::VertexState;
pub(crate) use super::super::edge_discret::{
    curve_tessellator_value_ok, split_by_deflection2d, tessellator_min_points, CurveOnSurface,
    CurveTessellator, EdgeDiscret, EdgeParameterProvider,
    MeshFace as UvFace,
};
pub(crate) use super::super::shape_tool::ShapeTool;
pub(crate) use super::super::face_discret::FaceChecker;
pub(crate) use super::super::mesh_tool::MeshTool;
pub(crate) use super::super::model_builder::{ModelBuilder, ModelPreProcessor};
pub(crate) use super::super::model_healer::ModelHealer;
pub(crate) use super::super::node_insertion::{DelaunayNodeInsertionMeshAlgo, TriangulationResult};
pub(crate) use super::super::parameters::MeshParameters;
pub(crate) use super::super::range_splitter::{classify_surface, SurfaceType};
pub(crate) use super::super::triangulator::{FaceTriangulation, Triangulator};

}


mod discret_root;
mod stitching;
pub use discret_root::*;
pub use stitching::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
