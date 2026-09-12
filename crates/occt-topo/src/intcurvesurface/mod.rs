//! Exact curve–surface intersection.
//!
//! Port of `IntCurveSurface_HInter` (TKGeomAlgo) at the level the precise
//! boolean needs: given a parametric curve and a surface, compute every point
//! (and, for the coincidence degeneracy, segment) where the curve meets the
//! surface.
//!
//! Dispatch mirrors the OCCT `PerformBounds` → `PerformConicSurf` /
//! `InternalPerform` split:
//!
//! - **Analytic path** for conic curves against quadric surfaces: a line
//!   against a plane/sphere/cylinder/cone (substitute the line's affine
//!   parametrization into the surface's implicit equation → linear/quadratic
//!   polynomial), a line against a torus (`IntAna_IntLinTorus` quartic), and a
//!   circle/ellipse against a plane or sphere (`a·cos + b·sin = c` via
//!   `math_TrigonometricFunctionRoots`). Tangency / coincidence produces a
//!   zero-length tangent point or an `On` segment.
//! - **General path** for everything else (B-spline curves, non-quadric
//!   surfaces): dense curve sampling, signed-distance sign-change bisection for
//!   transverse crossings, and golden-section refinement of near-zero distance
//!   minima for tangencies.
//!
//! Trait objects (`Arc<dyn Curve>` / `Arc<dyn Surface>`) cannot be downcast, so
//! analytic dispatch classifies curves and surfaces by geometric invariants —
//! the established pattern in this port (see `brep_surface::classify_surface`).
mod prelude {

pub(crate) use std::cmp::Ordering;

pub(crate) use occt_core::gp::{GpAx3, GpDir, GpLin, GpPnt, GpTorus, GpVec};
pub(crate) use occt_geom::intana::line_torus_intersect;
pub(crate) use occt_geom::{Curve, Surface};

pub(crate) use crate::brep_surface::{is_planar, sphere_center, surface_closest_params, surface_normal};

}


mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
