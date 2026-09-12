//! Rolling-ball fillet along an edge chain — a simplified port of
//! `BRepFilletAPI_MakeFillet` / `ChFi3d` (constant-radius fillet).
//!
//! A sharp convex edge of a solid is replaced by a cylindrical blend surface
//! tangent to both adjacent planar faces at radius `R` (a quarter cylinder for
//! a 90° corner). The two adjacent faces are trimmed to their tangency lines,
//! and the two end faces (perpendicular to the edge) are rebuilt with a
//! quarter-circular arc replacing the corner vertex. Edge instances are shared
//! between faces so the rebuilt shell satisfies `shell_is_closed`.
//!
//! Non-planar adjacent faces and reflex (concave) corners are rejected with an
//! error. The corner patch (`fillet_corner_solid`) is implemented as a
//! spherical octant — see the note there.
mod prelude {

pub(crate) use std::collections::HashMap;
pub(crate) use std::f64::consts::PI;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx2, GpAx3, GpCylinder, GpDir, GpPnt, GpSphere, GpVec};
pub(crate) use occt_geom::{GeomCylinder, GeomSphere, Surface};

pub(crate) use crate::brep_surface::{classify_surface, face_plane, SurfaceKind};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
pub(crate) use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, is_same, vertex_position, wires_of_face,
};

}


mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
