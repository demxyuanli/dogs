//! `BRepFill::Face`, `BRepFill::Shell`, `BRepFill::ComputeACR`,
//! `BRepFill::InsertACR` and `BRepFill::SearchOrigin`.
//!
//! Source: `BRepFill.cxx:211-422` (`Face`), `cxx:426-667` (`Shell`),
//! `cxx:1027-1063` (`ComputeACR`), `cxx:1067-1155` (`InsertACR`) and
//! `cxx:885-1023` (`SearchOrigin`), plus the
//! two file-local statics `MakeWire` (`cxx:81-117`) and `CutEdge`
//! (`cxx:119-142`). `TrimEdge` (`cxx:143-207`) is shared with
//! `BRepFill_CompatibleWires` and lives in
//! [`crate::brep_fill_compatible_wires::statics::trim_edge`].
//!
//! `Face` and `Shell` build the ruled patch with `GeomFill_Generator`
//! (`cxx:295-298`, `cxx:524-527`).
//!
//! Translation notes for [`search_origin`]:
//! - `BRepExtrema_DistShapeShape(V, W)` (`cxx:904`) is re-derived from the two
//!   solution families the class itself enumerates for a vertex / wire pair:
//!   every vertex of `W` (vertex-vertex extrema, exact point distance) and every
//!   edge of `W` (point-curve extrema through
//!   [`crate::brep_extrema::closest_point_on_edge`], the port's analytic
//!   `Extrema_ExtPC`). `SupportTypeShape2` (`cxx:918`) is `IsVertex` exactly
//!   when the vertex family wins; the loop at `cxx:909-916` only keeps the
//!   minimum, so the solution order does not matter except for an exact tie.
//! - `ParOnEdgeS2(isol, theparam)` (`cxx:944`) is the parameter
//!   [`crate::brep_extrema::closest_point_on_edge`] returns for the winning
//!   edge.
//! - `BRepAdaptor_Curve(E).D1(theparam, Pe, Ve)` (`cxx:1005-1008`) is
//!   [`BRepTool::edge_curve_world`] evaluated at `theparam`; `E.Orientation()`
//!   is the edge's stored orientation.

mod prelude {

pub(crate) use std::sync::Arc;
pub(crate) use occt_core::gp::{GpAx3, GpDir, GpDir2d, GpPnt, GpPnt2d, GpVec};
pub(crate) use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};
pub(crate) use occt_core::toploc::TopLocLocation;
pub(crate) use occt_geom::curve::Curve;
pub(crate) use occt_geom::surface::Surface;
pub(crate) use occt_geom2d::line::Geom2dLine;
pub(crate) use occt_geom2d::curve::Curve2d;
pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::brep_lib_find_surface::FindSurface;
pub(crate) use crate::brep_lib_make_face::make_face_from_wire;
pub(crate) use crate::brep_lib_same_parameter::same_parameter_shape;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, Shell, TopoShape, Vertex, Wire};
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topo_tools::{
    edge_vertices, edges_of_wire, explorer_shapes, faces_of, is_same, map_shapes_and_ancestors,
    shape_map_key, wires_of_face,
};

}

mod statics;
mod face_shell;
mod acr;
mod search_origin;
mod axe;

pub use face_shell::face;
pub use face_shell::shell;
pub use acr::compute_acr;
pub use acr::insert_acr;
pub use search_origin::search_origin;
pub use axe::axe;
