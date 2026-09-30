//! `Bnd_Tools` — convert `Bnd_Box` / `Bnd_Box2d` to the BVH box type.
//!
//! Source: `Bnd_Tools.hxx`. This port uses [`occt_core::bnd::BndBox`] as the
//! 3-D BVH box (`Bnd2BVH(const Bnd_Box&)` is identity on the stored AABB).
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::{BndBox, BndBox2d};

use crate::bvh_box2d::BvhBox2d;

/// `Bnd_Tools::Bnd2BVH(const Bnd_Box2d&)`.
pub fn bnd2bvh(the_box: &BndBox2d) -> BvhBox2d {
    BvhBox2d::from_bnd(the_box)
}

/// `Bnd_Tools::Bnd2BVH(const Bnd_Box&)`.
pub fn bnd2bvh3d(the_box: &BndBox) -> BndBox {
    *the_box
}
