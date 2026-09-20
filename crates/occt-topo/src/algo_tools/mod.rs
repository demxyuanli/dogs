//! Build-time geometry helpers for the boolean PaveFiller — a port of the
//! `BOPTools_AlgoTools` / `BOPTools_AlgoTools2D` / `BOPTools_AlgoTools3D` /
//! `BOPTools_Set` utility families (TKBO).
//!
//! These are the low-level construction primitives the pave-filler uses while
//! building a boolean result:
//!
//! * **Vertex construction** — [`AlgoTools::make_new_vertex`] builds a vertex
//!   at a 3D point with a tolerance;
//! * **Edge construction** — [`AlgoTools::make_edge`] builds an edge from an
//!   intersection curve plus two endpoint vertices and parameters, bumping the
//!   vertex tolerances so they cover the curve ends (`BOPTools_AlgoTools::MakeEdge`);
//! * **Vertex coincidence** — [`AlgoTools::compute_vv`] decides whether a vertex
//!   and a point interfere within the summed tolerances (`ComputeVV`);
//! * **Point-on-edge** — [`AlgoTools::point_on_edge`] evaluates the edge's 3D
//!   curve at a parameter; [`AlgoTools::update_vertex`] grows a vertex's
//!   tolerance to cover such a point (`UpdateVertex`);
//! * **P-curves** — [`AlgoTools::make_pcurve`] / [`AlgoTools2D::edge_to_face`]
//!   build the edge→face UV curve by delegating to
//!   [`crate::pcurve_full::make_pcurve_full`];
//!   [`AlgoTools2D::adjust_pcurve_on_surf`] brings a pcurve inside the face UV
//!   bounds via [`crate::pcurve_full::trim_pcurve_to_face`]
//!   (`AdjustPCurveOnSurf`);
//! * **Surface normal** — [`AlgoTools::get_normal_to_surface`] computes the
//!   unit normal of a surface at `(u, v)` (`GetNormalToSurface`);
//! * **Shape sets** — [`BOPToolsSet::shape_list`] dedupes a shape collection
//!   and [`BOPToolsSet::type_count`] counts shapes of a given type.
mod prelude {

pub(crate) use std::collections::{HashMap, HashSet};
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
pub(crate) use occt_core::precision::CONFUSION;
pub(crate) use occt_geom::{Curve, Surface};
pub(crate) use occt_geom2d::curve::Curve2d;

pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::brep_class3d::SolidClassifier;
pub(crate) use crate::brep_extrema::closest_point_on_edge;
pub(crate) use crate::brep_measure::edge_length;
pub(crate) use crate::brep_surface::surface_normal;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::connexity_block::ConnexityBlock;
pub(crate) use crate::fclass2d::{FClass2d, FaceState};
pub(crate) use crate::pcurve_full;
pub(crate) use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, wires_of_face};
pub(crate) use crate::tshape::HandleTShape;

}


mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
