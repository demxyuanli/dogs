//! Port of OCCT `BRepMesh_Delaun` — incremental 2D Delaunay triangulation
//! (Bowyer–Watson / "algorithm of Watson") over UV points.
//!
//! Source: `src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_Delaun.{hxx,cxx}`
//! (88.9K — the full algorithm, including super-triangle setup, circle-cell
//! acceleration, constraint-edge insertion and frontier adjustment).
//!
//! ponytail: the OCCT circle index (`BRepMesh_CircleTool`/`CircleInspector`,
//! sibling `delaun_index.rs`) is re-implemented here as a `HashMap`-backed grid
//! over `DelaunCircle`s because the sibling file is still under construction.
//! Swap for `delaun_index::CircleTool` when it lands.
mod prelude {

pub(crate) use std::collections::{HashMap, HashSet};
pub(crate) use std::f64::consts::PI;

pub(crate) use occt_core::gp::{GpPnt2d, GpVec2d, GpXY};
pub(crate) use occt_core::precision::{ANGULAR, PCONFUSION, RESOLUTION};

pub(crate) use super::super::delaun_data::{DelaunDataStructure, DelaunSelector};
pub(crate) use super::super::delaun_types::{DelaunCircle, DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
pub(crate) use super::super::geom_tool::{GeomTool, IntFlag};

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
