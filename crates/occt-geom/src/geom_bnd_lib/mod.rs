//! Port of `GeomBndLib` (TKGeomBase) and its `BndLib_Add*` entry points:
//! exact bounding boxes of a curve or surface over a parameter range.
//!
//! Every module here only needs `occt_core` plus the [`crate::curve::Curve`] /
//! [`crate::surface::Surface`] / `occt_geom2d::Curve2d` traits, so `GeomBndLib`
//! belongs to this crate rather than to `occt-topo`.
//!
//! - [`curve3d`] - `GeomBndLib_Curve.cxx` dispatch (`box_curve`), the analytic
//!   conic arms, the Bezier/BSpline control-net arms and
//!   `GeomBndLib_OtherCurve`'s sampled envelope.
//! - [`surface3d`] - `GeomBndLib_Surface.cxx` (`add_surface`, `surface_kind`).
//! - [`analytic3d`] - the analytic surface arms of `GeomBndLib_Surface`.
//! - [`add2d`] - `BndLib_Add2dCurve::Add`.
//! - [`curve2d`] - `GeomBndLib_Curve2d.cxx` dispatch and the conic arms.
//! - [`bezier2d`] / [`bspline2d`] / [`offset2d`] / [`other_curve2d`] - the
//!   remaining `GeomBndLib_Curve2d` arms.
//! - [`sample2d`] / [`spline_helpers`] - the sampled envelope and the pole-index
//!   helpers the B-spline arms share.
//! - [`circle3d`] / [`inf3d`] / [`elclib2d`] / [`inf2d`] / [`elclib_to3d`] -
//!   shared box and `ElCLib` evaluator helpers.
pub mod add2d;
pub mod analytic3d;
pub mod bezier2d;
pub mod bspline2d;
pub mod circle2d;
pub mod circle3d;
pub mod curve2d;
pub mod curve3d;
pub mod elclib2d;
pub mod elclib_to3d;
pub mod ellipse2d;
pub mod hyperbola2d;
pub mod inf2d;
pub mod inf3d;
pub mod line2d;
pub mod offset2d;
pub mod other_curve2d;
pub mod parabola2d;
pub mod sample2d;
pub mod spline_helpers;
pub mod surface3d;

pub use curve3d::box_curve;
pub use surface3d::{add_surface, surface_kind, GeomAbsSurfaceType};
