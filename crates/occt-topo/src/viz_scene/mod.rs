//! Phase 6 module: viz_scene — scene graph, camera, and raster pipeline.
//!
//! **UNPORTED (audit A14)**: port-local scene graph / camera / rasterizer. OCCT's
//! visualization stack is TKV3d (`AIS_InteractiveContext`, `V3d_View`, OpenGl
//! drivers); nothing here is a translation of it.
//!
//! A lightweight port of OCCT's `AIS_Shape` (scene item), `V3d_View`
//! (camera) and `V3d_Viewer` (projection) for offline rendering. A
//! [`VizScene`] holds [`SceneShape`]s — each a `TopoShape` plus a world
//! transform and an optional color. A [`Camera`] provides a right-handed
//! look-at view with perspective or orthographic projection; world points are
//! mapped to screen pixels by [`project_point`].
//!
//! Two renderers are provided:
//! - [`render_scene_svg`] — triangles are projected, depth-sorted
//!   (painter's algorithm) and emitted as filled `<polygon>`s; a wireframe
//!   variant ([`render_scene_svg_wireframe`]) emits the projected edges.
//! - [`render_scene_ppm`] / [`render_scene_raster`] — a ray is cast per
//!   pixel through a BVH over the scene's triangles; hits are shaded with a
//!   Lambert term against a front light and written as binary PPM.
//! - [`render_scene_ppm_shaded`] / [`render_scene_raster_zbuffer`] — a
//!   hardware-style software renderer: per-shape [`Material`]s, point
//!   [`Light`]s, an ambient term and flat / Gouraud / Phong shading, produced
//!   either by per-pixel ray casting (the nearest BVH hit acts as the depth
//!   buffer) or by a screen-space z-buffer triangle rasterizer.
//!
//! Camera interaction ports a subset of `V3d_View`: orbit, pan, zoom, fit-all
//! and a world-space picking ray ([`Camera::camera_orbit`],
//! [`Camera::camera_pan`], [`Camera::camera_zoom`], [`Camera::camera_fit`],
//! [`Camera::camera_ray`]).
mod prelude {

pub(crate) use std::cmp::Ordering;

pub(crate) use occt_core::bnd::BndBox;
pub(crate) use occt_core::bvh::bvh_ops::bvh_ray_cast;
pub(crate) use occt_core::bvh::builder_tri::build_tri_bvh;
pub(crate) use occt_core::gp::{GpMat, GpPnt, GpPnt2d, GpTrsf, GpVec, GpXyz};
pub(crate) use occt_core::poly::triangulation::Triangle;

pub(crate) use crate::mesh::ShapeMesh;
pub(crate) use crate::shape::TopoShape;

}


mod p01;
mod p02;
mod p03;
mod p04;
mod p05;
mod p06;
pub use p01::*;
pub use p02::*;
pub use p03::*;
pub use p04::*;
pub use p05::*;
pub use p06::*;

#[cfg(test)]
mod tests;
