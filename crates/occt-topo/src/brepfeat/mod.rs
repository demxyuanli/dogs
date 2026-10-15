//! TKFeat feature modeling — draft, groove, neck, rib, through-boss.
//! Source: `BRepFeat_MakeDPrism`, `BRepFeat_MakeRevol`, `BRepFeat_MakeLinearForm`.
//!
//! A "feature" modifies an existing solid by adding or removing material in a
//! locally parametrized way. This module ports the common TKFeat classes:
//!
//! - `BRepFeat_MakeDPrism` (tapered draft) → [`draft`];
//! - `BRepFeat_MakeRevol` (revolved feature) → [`groove`] (cut) and [`neck`]
//!   (fuse);
//! - `BRepFeat_MakeLinearForm` (rib / slot) → [`rib`];
//! - a through-boss (a `BRepPrimAPI_MakeCylinder` fused across the solid) →
//!   [`boss_thru_all`].
//!
//! Each operation builds a tool solid and combines it with the target:
//! - `draft` rotates the target face's supporting plane about the hinge (the
//!   intersection line of the face plane and the pivot plane) and swaps the
//!   face's surface geometry in a deep copy of the solid;
//! - `groove` revolves a 2D `(r, z)` profile into an annular tool and cuts it
//!   out of the solid;
//! - `neck` revolves a 2D profile into a solid tool and fuses it on;
//! - `rib` extrudes the profile's in-plane footprint into a thin plate of the
//!   requested wall `thickness` and standing `height`, fused on;
//! - `boss_thru_all` fuses a cylinder tall enough to pass through the solid's
//!   bounding-box height.
//!
//! # Geometry strategy
//!
//! `BRepFeat_*` takes the tool `TopoDS_Shape` as given and hands it to
//! `BOPAlgo_BOP`. The through-boss/neck cylinder tool is the analytic
//! `BRepPrimAPI_MakeCylinder` (`translated_cylinder`); the revolved tool
//! currently goes through the lathe-mesh fixture on [`revolve_profile_about`]
//! (see its note: the analytic `sweep_revolve` solid is not yet accepted by the
//! boolean). Audit A12 / T-125.
//!
//! # Booleans
//!
//! Every feature runs [`crate::bop_builder::boolean`] (`BOPAlgo_BOP`) on the
//! target solid and the tool; volumes are computed with
//! [`solid_volume`] (`BRepGProp::VolumeProperties` / `BRepGProp_Vinert`), the
//! divergence-theorem surface integral over each face of the result.
//!
//! # Example
//! ```ignore
//! let box_s = BRepPrimBox::make_box(2.0, 2.0, 1.0);
//! let profile = [GpPnt2d::new(0.0, 1.0), GpPnt2d::new(0.4, 1.0),
//!                GpPnt2d::new(0.4, 2.0), GpPnt2d::new(0.0, 2.0)];
//! let axis = GpAx1::new(GpPnt::new(1.0, 1.0, 0.0), GpDir::new(0.0, 0.0, 1.0)?);
//! let result = neck(&box_s.solid, &profile, &axis, 16, 1e-4)?;
//! assert!(result.volume > 4.0);
//! # Ok::<(), String>(())
//! ```
mod prelude {

pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx1, GpDir, GpMat, GpPln, GpPnt, GpPnt2d, GpTrsf, GpVec, TrsfForm};
pub(crate) use occt_geom::GeomPlane;

pub(crate) use crate::bop_builder::BoolOp;
pub(crate) use crate::shape::{Solid, TopoShape};
pub(crate) use crate::tgeometry::{FaceGeom, GeometryRegistry};
pub(crate) use crate::topo_tools_full::faces_of;

}


mod features;
pub use features::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
