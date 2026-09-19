//! 2D point-in-face classifier — a port of `IntTools_FClass2d` (TKBO).
//!
//! Classifies a 2D point `(u, v)` in the parameter domain of a face against
//! the face's UV boundary. The face boundary is sampled edge-by-edge into 2D
//! polygons (pcurves preferred, sampled projection as fallback), chained into
//! closed rings by UV continuity, and decomposed into one outer ring plus hole
//! rings (largest |area| ring = outer). A point is:
//!
//! * `On`  — within `tol` of a boundary ring segment;
//! * `In`  — inside the outer ring and outside every hole ring;
//! * `Out` — otherwise.
//!
//! Periodic surfaces (cylinder, cone, sphere, torus) fold out-of-range `u`/`v`
//! by whole periods (`AdjustPeriodic` semantics) before classifying, and try
//! the periodic images in order so `perform((u + 2π, v))` agrees with
//! `perform((u, v))`.
//!
//! Source: `IntTools_FClass2d.hxx/.cxx` (TKBO). The UV rings are built by
//! sampling each edge pcurve (via [`crate::pcurve_full::make_pcurve_full`],
//! falling back to [`crate::brep_surface::edge_pcurve_on_face`]); edge
//! orientation within a wire is recovered by UV-continuity chaining, since this
//! port's flat `TShape` child tree drops the per-edge `TopAbs_Orientation`.
mod prelude {

pub(crate) use std::f64::consts::PI;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::cslib::{Class2d, Class2dResult};
pub(crate) use occt_core::geom::polygon_ops::{point_in_polygon2d, polygon_area2d};
pub(crate) use occt_core::gp::GpPnt2d;
pub(crate) use occt_core::precision::SQUARE_CONFUSION;
pub(crate) use occt_geom::Surface;
pub(crate) use occt_geom2d::curve::Curve2d;

pub(crate) use crate::abs::Orientation;
pub(crate) use crate::brep_surface::{edge_pcurve_on_face, face_uv_bounds};
pub(crate) use crate::curve_sampling_2d::nb_samples;
pub(crate) use crate::pcurve_full::make_pcurve_full;
pub(crate) use crate::shape::{Edge, Face};
pub(crate) use crate::tgeometry::GeometryRegistry;
pub(crate) use crate::topo_tools_full::{edges_of_wire, wires_of_face};

}


mod p01;
pub use p01::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
