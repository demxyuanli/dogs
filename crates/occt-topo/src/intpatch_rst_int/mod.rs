//! `IntPatch_RstInt` (`IntPatch_RstInt.hxx/.cxx`).
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_RstInt.cxx`.
//!
//! * `geom_helpers.rs` - the file-statics (`cxx:48-429`): `Recadre`, the two
//!   `Tol3d`, `CoincideOnArc`, `VerifyTgline`, `GetLinePoint2d`,
//!   `FindParameter`.
//! * `put_vertex_on_line.rs` - `PutVertexOnLine` (`cxx:433-1355`).
//! * `vertex_write.rs` - its two write-back arms (`cxx:978-1227`).
//!
//! `IntPatch_RstInt` has a single static member in OCCT; nothing is kept on
//! `self`, so the class collapses to the free function
//! [`put_vertex_on_line`].

pub mod geom_helpers;
pub mod put_vertex_on_line;
pub mod vertex_write;

pub use put_vertex_on_line::put_vertex_on_line;
