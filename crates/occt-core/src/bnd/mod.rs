//! Bounding volumes — axis-aligned boxes, oriented boxes, bounding spheres.
pub mod box3d;
pub mod box2d;
pub mod bsphere;
pub mod obb;
pub mod sortbox;
pub mod b2b3;
pub mod obb_pca;
pub mod intersect;

pub use box3d::BndBox;
pub use box2d::BndBox2d;
pub use bsphere::BndSphere;
pub use obb::BndOBB;
pub use sortbox::{BoundSortBox, SortBoxElement};
