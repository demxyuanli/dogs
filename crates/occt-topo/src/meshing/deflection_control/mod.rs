//! Port of OCCT `BRepMesh_DelaunayDeflectionControlMeshAlgo` — deflection
//! controlled refinement of a 2D Delaunay triangulation (Wave 3).
//!
//! Source:
//! `src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx`
//!
//! The OCCT class is a template extending `BRepMesh_DelaunayNodeInsertionMeshAlgo`
//! (the sibling `node_insertion` stub). It walks the generated UV-space mesh and
//! splits every triangle whose linear deflection (chord sagitta) or angular
//! deflection (surface normal change) exceeds the requested tolerances:
//!
//! * `compute_triangle_geometry` evaluates the triangle's 3D normal and 2D area,
//!   rejecting degenerate triangles.
//! * `split_links` checks each link's midpoint: the `LineDeviation` functor
//!   measures how far the surface point sits off the 3D chord, and
//!   `check_link_ends_for_angular_deviation` measures the surface-normal angle
//!   between the two endpoints. Control points are collected and handed to the
//!   caller, which re-triangulates, until convergence or `MinSize` rejection.
//!
//! `reject_by_min_size` uses `CircleTool::Select` on the face-basis point
//! (`Scale(p, true)`), then checks only those shot triangles (`hxx:422-454`).
mod prelude {

pub(crate) use std::collections::HashSet;

pub(crate) use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpXY};
pub(crate) use occt_core::precision::{ANGULAR, CONFUSION, RESOLUTION};
pub(crate) use occt_geom::Surface;

pub(crate) use super::super::delaun::Delaun;
pub(crate) use super::super::delaun_types::{DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
pub(crate) use super::super::geom_tool::GeomTool;
pub(crate) use super::super::parameters::MeshParameters;

}


mod p01;
pub use p01::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
