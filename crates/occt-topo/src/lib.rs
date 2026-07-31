//! OCCT topology (TKBRep) — boundary representation data structures.
//!
//! Maps OCCT's TopoDS shape hierarchy to Rust:
//! ```text
//! Shape (enum: Vertex/Edge/Wire/Face/Shell/Solid/Compound)
//!   ├── TShape: geometric + topological data
//!   ├── Location: optional transform
//!   └── Orientation: FORWARD/REVERSED/INTERNAL/EXTERNAL
//! ```

pub mod abs;
pub mod tshape;
pub mod shape;
pub mod builder;
pub mod iterator;
pub mod tools;

pub use abs::{ShapeType, Orientation};
pub use tshape::{TShape, VertexShape, EdgeShape, WireShape, FaceShape, ShellShape, SolidShape, CompoundShape};
pub use shape::{TopoShape, Vertex, Edge, Wire, Face, Shell, Solid, Compound};
pub use builder::TopoBuilder;
pub use iterator::ShapeIterator;
