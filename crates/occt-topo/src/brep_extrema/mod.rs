//! Shape-to-shape distance and extrema queries.
//!
//! Port of `BRepExtrema_DistShapeShape` (TKBRep) plus a few lightweight
//! helpers (closest vertex, point-in-solid, bbox separation).
//!
//! `BRepExtrema_DistShapeShape` computes exact distances between shapes via
//! extremum search on the underlying geometry. This port delegates the
//! point–curve and point–surface minima to the Phase 13 analytic/Newton
//! solvers in `occt_geom` (`extrema_pc` / `extrema_surf`), then clamps the
//! parameters to the bounded edge / face range and re-checks the face
//! boundary edges (a bounded face's closest point can lie on an edge).


mod bbox_ray;
mod closest;
mod collect;
mod distance;
mod inside;
mod mesh_sampling;
mod uv_projection;

#[cfg(test)]
pub(crate) mod test_box;
#[cfg(test)]
mod tests;

pub(in crate::brep_extrema) use bbox_ray::*;
pub use closest::*;
pub(in crate::brep_extrema) use collect::*;
pub use distance::*;
pub use inside::*;
pub use mesh_sampling::*;
pub(in crate::brep_extrema) use uv_projection::*;
