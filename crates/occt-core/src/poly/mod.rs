//! Polygon and triangulation data structures. Source: `Poly/`
pub mod triangulation;
pub mod triangulation_full;
pub mod coherent;
pub mod connect;
pub mod make_loops;
pub mod merge_nodes;
pub mod polygon3d;
pub mod polygon3d_full;
pub mod polygon2d;

pub use triangulation::Triangulation;
pub use polygon2d::Polygon2D;
pub use polygon3d_full::Polygon3D;
