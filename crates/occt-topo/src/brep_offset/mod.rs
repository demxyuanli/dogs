//! Parallel offsets of curves, faces, shells and solids.
//! Source: `BRepOffsetAPI_MakeOffsetShape`, `BRepOffsetAPI_MakeOffset`,
//! `BRepOffset_Offset`.
mod prelude {

pub(crate) use std::collections::HashMap;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpVec, GpXyz};
pub(crate) use occt_geom::surface_fit::fit_plane;
pub(crate) use occt_geom::{Curve, GeomCircle, GeomCylinder, GeomPlane, GeomSphere, Surface};

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::brep_builder_api::make_edge_arc;
pub(crate) use crate::brep_surface::{face_is_planar, face_plane, is_planar, sphere_center, surface_normal};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
pub(crate) use crate::shape_ops::translated_copy;
pub(crate) use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, vertices_of, wires_of_face};

}


mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
