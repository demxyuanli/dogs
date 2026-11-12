//! TKFeat feature modeling — draft, groove, neck, rib, through-boss, holes.
//! Source: `BRepFeat_MakeDPrism`, `BRepFeat_MakeRevol`, `BRepFeat_MakeLinearForm`,
//! `BRepFeat_MakeCylindricalHole`.
//!
//! A "feature" modifies an existing solid by adding or removing material in a
//! locally parametrized way. This module ports the common TKFeat classes:
//!
//! - `BRepOffsetAPI_DraftAngle` (face draft) → [`draft`];
//! - `BRepFeat_MakeRevol` (revolved feature) → [`groove`] (cut) and [`neck`]
//!   (fuse);
//! - `BRepFeat_MakeLinearForm` (rib / slot) → [`rib`] (`Init` + `Perform`);
//! - `BRepFeat_MakePrism` → [`boss`] (`Perform(Length)`) and
//!   [`boss_thru_all`] (`PerformThruAll`);
//! - `BRepFeat_MakeCylindricalHole` to [`through_hole`] (`Perform(Radius)`) and
//!   [`blind_hole`] (`PerformBlind`).
//!
//! Each operation builds a tool solid and combines it with the target:
//! - `draft` runs `BRepOffsetAPI_DraftAngle::Add` / `Build` on the selected
//!   face (`Direction` = face normal, `NeutralPlane` = pivot);
//! - `groove` / `neck` run `BRepFeat_MakeRevol` over the cross-section face;
//! - `rib` runs `BRepFeat_MakeLinearForm` over a sketch wire, plane, and the
//!   pair of opposite in-plane thickness vectors (`Direc` / `Direc1`);
//! - `boss` / `boss_thru_all` run `BRepFeat_MakePrism` over a circular
//!   generating face;
//! - `through_hole` / `blind_hole` cut a cylinder along an axis through
//!   `BRepFeat_MakeCylindricalHole`.
//!
//! # Geometry strategy
//!
//! Faithful TKFeat entry points keep the Boolean algorithm alive for
//! `Modified` / `IsDeleted` (`brep_feat::algo_api`). [`feature_apply`] uses the
//! same facade when a caller still builds a tool solid by hand.
//!
//! # Booleans
//!
//! Feature Booleans go through [`crate::brep_feat::algo_api`] (`BRepAlgoAPI_*`
//! over `BOPAlgo_BOP` with history). Volumes use [`solid_volume`]
//! (`BRepGProp::VolumeProperties` / `BRepGProp_Vinert`).
//!
//! # Example
//! ```ignore
//! let box_s = BRepPrimBox::make_box(2.0, 2.0, 1.0);
//! let profile = [GpPnt2d::new(0.0, 1.0), GpPnt2d::new(0.4, 1.0),
//!                GpPnt2d::new(0.4, 2.0), GpPnt2d::new(0.0, 2.0)];
//! let axis = GpAx1::new(GpPnt::new(1.0, 1.0, 0.0), GpDir::new(0.0, 0.0, 1.0)?);
//! let result = neck(&box_s.solid, &profile, &axis, 1e-4)?;
//! assert!(result.volume > 4.0);
//! # Ok::<(), String>(())
//! ```
mod prelude {

pub(crate) use occt_core::gp::{GpAx1, GpDir, GpPln, GpPnt, GpPnt2d, GpTrsf, GpVec};

pub(crate) use crate::bop_builder::BoolOp;
pub(crate) use crate::shape::{Solid, TopoShape};
pub(crate) use crate::topo_tools::faces_of;

}


mod features;
pub use features::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
