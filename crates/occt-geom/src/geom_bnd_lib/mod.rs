//! Port of `GeomBndLib` (TKGeomBase) and its `BndLib_Add3dCurve::Add` entry
//! point: exact bounding boxes of a 3D curve over a parameter range.
//!
//! These modules used to live in `occt-topo` (`geom_bnd_lib_curve3d.rs` and
//! friends) even though every one of them only needs `occt_core` plus
//! [`crate::curve::Curve`] -- `GeomBndLib` belongs to TKGeomBase, i.e. this
//! crate. `occt-geom` consumers (the GeomFill classes, which read the section
//! box for their `Gabarit`) cannot reach `occt-topo` (that crate depends on
//! this one), so the faithful port moved here; `occt-topo` re-exports it to
//! keep its call sites unchanged.
//!
//! - [`curve3d`] - `GeomBndLib_Curve.cxx` dispatch (`box_curve`), the analytic
//!   conic arms, the Bezier/BSpline control-net arms and
//!   `GeomBndLib_OtherCurve`'s sampled envelope.
//! - [`circle3d`] - `GeomBndLib_Circle.cxx`.
//! - [`inf3d`] - the infinite/open box helpers shared by the 3D arms.
//! - [`elclib2d`] - the `ElCLib` 2D evaluators those arms use.
pub mod circle3d;
pub mod curve3d;
pub mod elclib2d;
pub mod inf3d;

pub use curve3d::box_curve;
