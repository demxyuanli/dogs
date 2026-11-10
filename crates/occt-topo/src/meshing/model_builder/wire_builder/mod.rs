//! Model wire builder (BRepMesh_ModelBuilder / ShapeAnalysis_Wire wire walk).
//!
//! Split by concern; the public surface keeps its previous paths through the
//! re-exports below.

mod add_wire;
mod builder;
mod next_param;
mod pcurve;
mod wire_edges;
mod wire_walk;

pub use builder::ModelBuilder;
pub use wire_edges::outer_of_wires;
pub(crate) use pcurve::sample_pcurve;
pub(crate) use wire_edges::wire_edges_explorer;
