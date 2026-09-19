//! Full edge→face UV curve (pcurve) construction — the Phase 16b superset.
//!
//! Extends [`crate::pcurve::make_pcurve_on_face`] from the plane/cylinder pair
//! to every analytic surface (plane, cylinder, cone, sphere, torus) and to
//! general (B-spline) faces:
//!
//! * **Plane** — the Phase 16a implementation is reused verbatim, so
//!   [`make_pcurve_full`] and [`crate::pcurve::make_pcurve_on_face`] agree
//!   exactly: a line edge becomes a `Geom2dLine`, a circle edge a
//!   `Geom2dCircle`, anything else a sampled degree-1 B-spline.
//! * **Cylinder / cone / sphere / torus** — the surface's own parameterization
//!   is recovered from sampled geometry invariants and used to project the edge
//!   exactly.  A generatrix (line edge) or latitude/meridian (circle edge) maps
//!   to a `u`- or `v`-isoline, returned as a `Geom2dLine` (unit-speed) or a
//!   degree-1 B-spline otherwise.  Non-isoparametric edges fall through to the
//!   sampling path.
//! * **Anything else (B-spline faces)** — the edge is sampled at 48 parameters
//!   and the projected `(u, v)` points are fitted with a clamped degree-1
//!   B-spline whose parameter range matches the edge range.
//!
//! Source: `BOPTools_AlgoTools2D::MakePCurveOnFace` / `AdjustPCurveOnSurf`
//! (TKBO) with the `ProjLib` analytic projection rules for the revolution
//! surfaces.
mod prelude {

pub(crate) use std::f64::consts::{FRAC_PI_2, PI};
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx2d, GpAx3, GpDir, GpDir2d, GpPln, GpPnt, GpPnt2d, GpTrsf2d, GpVec, GpVec2d};
pub(crate) use occt_geom::{Curve, Surface};
pub(crate) use occt_geom2d::curve::Curve2d;
pub(crate) use occt_geom2d::{Geom2dBSplineCurve, Geom2dLine};

pub(crate) use crate::brep_surface::{
    classify_surface, edge_pcurve_on_face, face_uv_bounds, sphere_center, SurfaceKind,
};
pub(crate) use crate::shape::{Edge, Face};
pub(crate) use crate::tgeometry::GeometryRegistry;

}


mod p01;
mod p02;
mod p03;
mod p04;
pub use p01::*;
pub use p02::*;
pub use p03::*;
pub use p04::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
