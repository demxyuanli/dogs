//! Exact boolean foundation — edge/face intersection queries.
//!
//! Ports the `IntTools_EdgeEdge`, `IntTools_EdgeFace` and `IntTools_FaceFace`
//! algorithms (TKBO) at the level the exact boolean needs:
//!
//! - **edge–edge**: exact for line–line (segment intersection), line–circle
//!   (quadratic in the circle's plane) and coplanar circle–circle (radical
//!   line); all other combinations fall back to a sampling+refine solver.
//! - **edge–face**: exact for line–plane and circle–plane; general curves
//!   cross a plane by sampled sign-change bisection; non-planar faces use a
//!   curve–surface sampler. Results are filtered to the face's 2D boundary
//!   polygon.
//! - **face–face**: exact for planar faces (plane–plane line clipped to both
//!   boundary polygons); non-planar faces yield an approximate polyline.
//!
//! `bop_builder` consumes these exact signatures to split edges at their
//! intersection parameters and build the boolean result.
mod prelude {

pub(crate) use occt_core::geom::polygon_ops::point_in_polygon2d;
pub(crate) use occt_core::gp::{GpPln, GpPnt, GpPnt2d, GpVec, GpVec2d};
pub(crate) use occt_core::int::curve_curve::segment_segment_intersection_3d;
pub(crate) use occt_geom::geom_api;
pub(crate) use occt_geom::{Curve, Surface};

pub(crate) use crate::brep_surface::face_plane;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::edge_split;
pub(crate) use crate::face_face::{plane_plane_intersection, FaceIntersect, face_face_intersection};
pub(crate) use crate::shape::{Edge, Face, TopoShape};
pub(crate) use crate::topo_tools_full::{edge_vertices, edges_of, edges_of_wire, vertex_position, wires_of_face};

}


mod intersections;
mod face_classification;
pub use intersections::*;
pub use face_classification::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
