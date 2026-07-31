//! Bounding volumes — axis-aligned boxes, oriented boxes, bounding spheres.
pub mod box3d;
pub mod box2d;
pub mod bsphere;
pub mod obb;

pub use box3d::BndBox;
pub use box2d::BndBox2d;
pub use bsphere::BndSphere;
pub use obb::BndOBB;
