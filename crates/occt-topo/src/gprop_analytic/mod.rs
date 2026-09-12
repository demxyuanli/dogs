//! Analytic mass properties and analytic extrema — exact `BRepGProp` + `Extrema`.
//!
//! Three families of routines:
//!
//! * **Mass properties** — exact surface area, volume and centroid for shapes
//!   built from analytic surfaces (planes, spheres, cylinders, cones, toruses).
//!   A box reports `2(dx·dy + dy·dz + dx·dz)` and `dx·dy·dz` exactly, unlike the
//!   tessellation-based [`crate::brep_gprop`] which inflates curved faces. Faces
//!   whose surface is not one of the analytic types fall back to the mesh area.
//!
//! * **Higher-order properties** — the 3×3 inertia tensor about the centroid
//!   ([`inertia_tensor`]), the principal moments and principal axes
//!   ([`principal_inertia`]), the density-weighted mass ([`analytic_mass`]), and
//!   the exact analytic curve length ([`analytic_curve_length`]). Analytic
//!   solids (box, sphere, cylinder, cone, torus) are integrated in closed form;
//!   other solids fall back to a vertex-mass approximation. The full set is
//!   assembled by [`full_analytic_properties`].
//!
//! * **Extrema** — closed-form nearest-point solutions for point/line/plane/
//!   sphere/cylinder/cone/torus (no grid search). Unsupported surfaces fall back
//!   to [`occt_geom::extrema::point_surface_extrema`].
//!
//! `Surface` trait objects cannot be downcast (the port keeps no `Any`-typed
//! geometry), so analytic parameters (radius, axis, …) are reconstructed from
//! sampled geometry invariants rather than read from a concrete `Gp*` struct.
//!
//! The inertia tensor is reported about the shape's centroid. For non-analytic
//! solids the vertex-mass fallback distributes the total mass equally over the
//! vertices and shifts the origin tensor to the centroid by the parallel-axis
//! theorem — a coarse approximation, documented per call site.
mod prelude {

pub(crate) use std::f64::consts::PI;

pub(crate) use occt_core::geom::polygon_ops::polygon_area3d;
pub(crate) use occt_core::gp::{GpDir, GpPln, GpPnt, GpVec};
pub(crate) use occt_geom::extrema::{point_surface_extrema, ExtremaPair};
pub(crate) use occt_geom::{GeomCircle, GeomLine, GeomPlane, GeomSphere, Surface};
pub(crate) use occt_math::eigen_ext::jacobi_eigen_symmetric;
pub(crate) use occt_math::{integrate, MathMatrix};

pub(crate) use crate::brep_surface::{
    classify_surface, face_is_planar, face_plane, sphere_center, surface_normal, SurfaceKind,
};
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::shape::{Edge, Face, TopoShape};
pub(crate) use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, vertex_position, vertices_of, wires_of_face,
};

}


mod p01;
mod p02;
mod p03;
pub use p01::*;
pub use p02::*;
pub use p03::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
