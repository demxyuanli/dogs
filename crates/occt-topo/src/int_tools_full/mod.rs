//! Full port of `IntTools_Context` (TKBO) — a cached intersection context.
//!
//! The context bundles the geometry/topology toolkit used by the boolean
//! intersection pipeline and caches the reusable tools so repeated queries on
//! the same shape do not rebuild them:
//!
//! - a per-face 2D classifier ([`FClass2d`]) cache, so point-in-face queries on
//!   the same face sample the UV boundary only once;
//! - point/vertex projection helpers onto edges and surfaces;
//! - the block/point validity predicates used to decide whether an
//!   intersection sub-range lies inside a face.
//!
//! Source: `IntTools_Context.hxx/.cxx` (TKBO, `ModelingAlgorithms/TKBO/IntTools`).
//! This is the *complete* context; the lighter [`crate::inttools_range::IntContext`]
//! (Phase 16) is a shell that fills the classifiers in later.
//!
//! The `is_valid_block_for_*` predicates are self-contained: with no
//! intersection curve handle in [`IntRange`], the 1D range is interpreted over
//! the face's first (`u`) parameter at the middle of its `v` range (the
//! u-midline of the face's UV domain).
mod prelude {

pub(crate) use std::collections::HashMap;

pub(crate) use occt_core::gp::{GpPnt, GpPnt2d};
pub(crate) use occt_core::precision::{CONFUSION, PCONFUSION};
pub(crate) use occt_geom::Curve;

pub(crate) use crate::abs::{Orientation, ShapeType};
pub(crate) use crate::brep_surface::{face_uv_bounds, is_planar, surface_closest_params};
pub(crate) use crate::edge_face_kind::plane_projection;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::fclass2d::{FaceState, FClass2d};
pub(crate) use crate::inttools_data::IntRange;
pub(crate) use crate::inttools_roots;
pub(crate) use crate::iterator::ShapeIterator;
pub(crate) use crate::shape::{Edge, Face, Vertex};
pub(crate) use crate::shape_naming::ShapeId;

}


mod p01;
pub use p01::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
