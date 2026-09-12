//! Discrete mesh data model — port of `BRepMeshData::{Model,Edge,Face,Wire,Curve,PCurve}`
//! and `IMeshData::{Status,StatusOwner,TessellatedShape}`.
//!
//! Mirrors OCCT's `IMeshData`/`BRepMeshData` layer. A [`MeshModel`] holds
//! indexed collections of discrete edges, faces and wires:
//!
//! * [`MeshEdge`] carries the topological [`Edge`], the 3D curve geometry
//!   (`Arc<dyn Curve>`), its parameter range, a [`MeshCurve`] discretization
//!   and per-face [`MeshPCurve`]s.
//! * [`MeshFace`] carries the topological [`Face`], the surface geometry
//!   (`Arc<dyn Surface>`), its wires (referenced by model index) and internal
//!   face points.
//! * [`MeshWire`] is an ordered chain of edge indices + orientations.
//! * Every entity is a [`MeshStatus`] status owner; [`MeshModel::status_mask`]
//!   aggregates the statuses of all contained entities.
mod prelude {

pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpPnt, GpPnt2d};
pub(crate) use occt_geom::{Curve, Surface};

pub(crate) use crate::abs::Orientation;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Wire};
pub(crate) use crate::tgeometry::GeometryRegistry;

}


mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
