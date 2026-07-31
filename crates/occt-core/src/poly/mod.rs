//! Polygon and triangulation data structures. Source: `Poly/`
pub mod triangulation;
pub mod polygon3d;
pub mod polygon3d_full;
pub mod polygon2d;

pub use triangulation::Triangulation;
pub use polygon2d::Polygon2D;
pub use polygon3d_full::Polygon3D;
