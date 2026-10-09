//! Curved-face rolling-ball fillet — a port of `ChFi3d` blending between a
//! curved surface and a plane, or between two curved surfaces.
//!
//! `fillet_edge` only blends two planar faces (the blend is a cylinder). This
//! module extends the rolling-ball idea to edges where at least one adjacent
//! face is a sphere or a cylinder:
//!
//! * **Plane + Sphere** — the rolling ball stays tangent to the plane and to
//!   the sphere, so its centre sweeps the *intersection of the offset plane*
//!   (the plane pushed out by `r`) *and the offset sphere* (radius `R_s + r`).
//!   That intersection is a circle; the envelope of the ball along it is a
//!   **torus** with major radius = circle radius and minor radius = fillet
//!   radius `r`.
//! * **Plane + Cylinder** (axis perpendicular to the plane) — the ball centre
//!   sweeps the intersection of the offset plane and the offset cylinder
//!   (radius `R_c + r`), again a circle, and the blend is a **torus** band
//!   tangent to the cylinder at the inner equator and to the plane at the
//!   outer bottom.
//! * **Sphere + Sphere** — the ball centre sweeps the intersection of the two
//!   offset spheres (radius `R_i + r`), a circle, and the blend is a torus.
//!
//! In every supported case the blend face is a `GeomTorus` band bounded by the
//! two tangency circles (one per adjacent face). The two adjacent faces are
//! rebuilt: the planar face keeps its outer boundary and gets the enlarged
//! plane-tangency circle as a hole; the curved face gets the smaller
//! surface-tangency circle as its new boundary. A cylinder's lateral face is
//! rebuilt with the contact circle, a re-generated seam sub-edge and the
//! original top rim. All new edges are shared between the blend face and the
//! trimmed adjacent face, so the rebuilt shell satisfies `shell_is_closed`.
//!
//! Unsupported pairs (cone, torus, B-spline, a cylinder whose axis is not
//! perpendicular to the plane, ...) return `Err` — documented, never a panic.
mod prelude {

pub(crate) use std::f64::consts::PI;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx2, GpAx3, GpCylinder, GpDir, GpPnt, GpSphere, GpTorus, GpVec};

pub(crate) use occt_geom::{GeomCylinder, GeomSphere, GeomTorus, Surface};

pub(crate) use crate::brep_surface::{
    is_planar, sphere_center, surface_normal, SurfaceKind,
};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::fillet_edge::classify_surface_full;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Wire};
pub(crate) use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, is_same, vertex_position, vertices_of,
    wires_of_face,
};

}


mod surface_info;
mod meridional;
mod torus_blend;
mod rolling_ball;
mod corner_patch;
pub use surface_info::*;
pub use meridional::*;
pub use torus_blend::*;
pub use rolling_ball::*;
pub use corner_patch::*;

#[cfg(test)]
mod tests;
