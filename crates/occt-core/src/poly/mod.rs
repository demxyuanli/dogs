//! Polygon and triangulation data structures. Source: `Poly/`
pub mod triangulation;
pub mod polygon3d;
pub mod polygon2d;

pub use triangulation::Triangulation;
pub use polygon3d::Polygon3D;
pub use polygon2d::Polygon2D;
