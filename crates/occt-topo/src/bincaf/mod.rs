//! Binary XCAF container — a compact binary assembly/attributes document.
//! Source: `BinXCAF` (binary `XCAFDoc_*` document driver).
//!
//! A `BinXcaf` is a tree of entries (`BinXcafEntry`), each carrying an
//! optional `TopoShape`, a list of string attributes (name / color / layer /
//! material …) and child entries. The shape is stored topologically: type,
//! then per-type geometry snapshots (vertex points, edge curves, face
//! surfaces with parameter ranges, and the child sub-shape lists). Curves and
//! surfaces cannot be downcast from `Arc<dyn Curve>` / `Arc<dyn Surface>`, so
//! they are classified by sampling invariants (zero second derivative ⇒ line,
//! periodic ⇒ circle, planar ⇒ plane, equidistant samples ⇒ sphere), exactly
//! like the STEP/IGES writers.
//!
//! Binary layout (all integers little-endian):
//! `BINXCAF` magic (8 bytes) · version `u32` · entry tree.
//! One entry: `u8 has_shape` · shape · `i32` attribute count · `(len-prefixed
//! kind, len-prefixed value)` pairs · `u32` child count · child entries.
mod prelude {

pub(crate) use std::f64::consts::PI;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx1, GpAx2, GpAx3, GpCirc, GpDir, GpPln, GpPnt, GpSphere, GpVec};
pub(crate) use occt_geom::{Curve, GeomCircle, GeomLine, GeomPlane, GeomSphere, Surface};

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::brep_surface::{classify_surface, face_plane, sphere_center, SurfaceKind};
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Wire};
pub(crate) use crate::tgeometry::{FaceGeom, GeometryRegistry};

}


mod format;
pub use format::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
