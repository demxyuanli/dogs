//! `Bnd_Tools` — convert `Bnd_Box` / `Bnd_Box2d` to the BVH box type.
//!
//! Source: `Bnd_Tools.hxx`. This port uses [`occt_core::bnd::BndBox`] as the
//! 3-D BVH box (`Bnd2BVH(const Bnd_Box&)` is identity on the stored AABB).

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
