//! Parallel offsets of curves, faces, shells and solids.
//! Source: `BRepOffsetAPI_MakeOffsetShape`, `BRepOffsetAPI_MakeOffset`,
//! `BRepOffset_Offset`.
mod prelude {

pub(crate) use std::collections::HashMap;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpCylinder, GpPln, GpPnt, GpSphere, GpVec, GpXyz};
pub(crate) use occt_geom::{Curve, GeomCircle, GeomCylinder, GeomPlane, GeomSphere, Surface};

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::brep_builder_api::make_edge_arc;
pub(crate) use crate::brep_surface::{face_is_plane_type, face_plane};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
pub(crate) use crate::shape_ops::translated_copy;
pub(crate) use crate::topo_tools::{edges_of, edges_of_wire, faces_of, vertices_of, wires_of_face};

}


mod analyse;
mod brep_offset_surface;
mod collapse_singularities;
mod curve_face_offset;
mod inter2d;
mod inter3d;
mod interval;
mod make_loops;
mod make_offset;
mod make_simple_offset;
mod offset;
mod offset_enums;
mod shell_offset;
mod simple_offset;
mod tool;
pub use analyse::*;
pub use brep_offset_surface::*;
pub use collapse_singularities::*;
pub use curve_face_offset::*;
pub use inter2d::*;
pub use inter3d::*;
pub use interval::*;
pub use make_loops::*;
pub use make_offset::*;
pub use make_simple_offset::*;
pub use offset::*;
pub use offset_enums::*;
pub use shell_offset::*;
pub use simple_offset::*;
pub use tool::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
