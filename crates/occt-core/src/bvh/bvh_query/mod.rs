//! BVH-accelerated mesh queries and `Poly` mesh topology tools.
//!
//! This module layers the standard OCCT query surface on top of the triangle
//! BVH built by [`crate::bvh::builder_tri::build_tri_bvh`], and adds a set of
//! mesh utility functions that OCCT exposes through the `Poly` package and the
//! `BRepMesh` meshing algorithms:
//!
//! * Ray, segment and box queries against a triangle mesh (pruned by the BVH).
//! * Point containment and closest-point queries.
//! * Integral quantities: signed volume (divergence theorem) and surface area.
//! * Edge/vertex topology: boundary edges, manifold-closed test, connected
//!   components, boundary loops and per-vertex incident triangle lists.
//! * 2D polygon ear-clipping triangulation.
//!
//! Source: `BVH_Tools.hxx`, `BVH_BoxSet.hxx`, `Poly_Connect.hxx`,
//! `Poly_Triangulation.hxx`, `BRepMesh`-style acceleration.
//!
//! All mesh-level queries take a [`TriBvh`] built over the *same* triangle
//! slice they are called with — the leaf ranges stored in the BVH nodes index
//! directly into that slice. The module therefore speaks two triangle
//! representations and bridges them:
//!
//! * `&[(GpPnt, GpPnt, GpPnt)]` — a flat triangle soup; this is what the BVH
//!   builder consumes and what every query here accepts.
//! * `&[(usize, usize, usize)]` — an index mesh referencing a vertex array;
//!   this is the form used by the topology helpers ([`mesh_edge_topology`],
//!   [`mesh_connected_components`], [`mesh_boundary_loops`], ...).
//!
//! [`triangulation_triangles`] converts a [`crate::poly::Triangulation`] into
//! the flat soup form so meshes built by the `Poly` package can be queried
//! directly.
//!
//! ## Relationship to OCCT classes
//!
//! | This module | OCCT |
//! |---|---|
//! | [`ray_cast_mesh`], [`segment_query_mesh`] | `BVH_Tools::RaySegmentIntersection` |
//! | [`box_query_union`] | `BVH_BoxSet::Select` / `BVH_Tree` box query |
//! | [`point_inside_box`] | `BRepMesh` even-odd containment |
//! | [`closest_point_mesh`] | `Poly_Connect` + point-triangle distance |
//! | [`mesh_edge_topology`] | `Poly_Connect` edge counts |
//! | [`mesh_to_polygon_indices`] | `Poly_Connect::TriangleToNodes` |
//!
//! ## Performance notes
//!
//! The ray, segment and box queries prune with the BVH first and only test the
//! triangles that land in candidate leaves, so their cost is logarithmic in
//! the leaf count for meshes with a reasonable spatial distribution. The
//! topology helpers ([`mesh_edge_topology`], [`mesh_connected_components`],
//! [`mesh_boundary_loops`], [`mesh_is_manifold_closed`]) operate directly on
//! index triples and do **not** require a BVH at all — they are linear in the
//! number of edges. [`closest_point_mesh`] combines a BVH candidate lookup
//! with a full scan fallback (see its documentation), so its worst case is
//! linear in the triangle count; the BVH only short-circuits on-surface
//! queries.
mod prelude {

pub(crate) use std::collections::HashMap;

pub(crate) use crate::bnd::BndBox;
pub(crate) use crate::bvh::builder_tri::{closest_triangle_to_point, query_triangles, TriBvh};
pub(crate) use crate::bvh::bvh_ops::{bvh_point_in_mesh, ray_hits_bbox, ray_triangle_t};
pub(crate) use crate::bvh::BvhNode;
pub(crate) use crate::gp::{GpPnt, GpPnt2d, GpVec};

}


mod mesh_queries;
mod triangulate;
pub use mesh_queries::*;
pub use triangulate::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
