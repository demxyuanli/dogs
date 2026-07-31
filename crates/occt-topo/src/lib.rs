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
pub mod mesh;
pub mod transform;
pub mod validate;
pub mod primitives;
pub mod tgeometry;

pub use abs::{ShapeType, Orientation};
pub use tshape::{TShape, VertexShape, EdgeShape, WireShape, FaceShape, ShellShape, SolidShape, CompoundShape};
pub use shape::{TopoShape, Vertex, Edge, Wire, Face, Shell, Solid, Compound};
pub use builder::TopoBuilder;
pub use iterator::ShapeIterator;
pub mod model;
pub mod topexp;
pub mod brep_tool;

// Phase 3 modules (filled in by task agents).
pub mod wireframe;
pub mod shape_mesh;
pub mod bbox_from_geometry;
pub mod face_face;
pub mod edge_split;
pub mod shell_check;
pub mod solid_union;
pub mod mesh_to_brep;
pub mod brep_exchange;
pub mod brep_scene;
pub mod brep_extrema;
pub mod shape_analysis;
pub mod brep_gprop;
pub mod sweep;
pub mod brep_surface;
pub mod topo_tools_full;
pub mod boolean_ops;
pub mod brep_measure;
// Phase 4 modules.
pub mod brep_builder_api;
pub mod fillet;
pub mod loft;
pub mod step;
pub mod hlr;
pub mod vrml;
pub mod sweep_revolve;
pub mod shape_ops;
pub mod shape_naming;
pub mod brep_pipe;
pub mod brep_faces;
pub mod brep_sketch;
pub mod brep_compare;
pub mod brep_assembly;
