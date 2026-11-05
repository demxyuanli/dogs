//! `MAT` package (`TKTopAlgo/MAT`): the medial axis graph that `MAT2d_Mat2d`
//! fills and `BRepMAT2d` / `BRepFill_OffsetAncestors` consume.
//!
//! | OCCT | lines | file |
//! |---|---|---|
//! | `MAT_BasicElt` | `MAT_BasicElt.cxx:27-90` | `elements.rs` |
//! | `MAT_Edge` | `MAT_Edge.cxx:22-81` | `elements.rs` |
//! | `MAT_Bisector` | `MAT_Bisector.cxx:25-206` | `elements.rs` |
//! | `MAT_Node` | `MAT_Node.cxx:24-150` | `elements.rs` |
//! | `MAT_Arc` | `MAT_Arc.cxx:27-271` | `elements.rs` |
//! | `MAT_Graph` | `MAT_Graph.cxx:42-582` | `graph.rs` |
//! | `MAT_Zone` | `MAT_Zone.cxx:26-179` | `zone.rs` |
//!
//! `MAT_ListOfBisector` / `MAT_ListOfEdge` (`MAT_ListOfBisector_0.cxx`,
//! `MAT_ListOfEdge_0.cxx`, 28 lines each) are `MAT_TList` instantiations
//! (`MAT_TList.gxx`), ported as [`MatTList`] over the bisector / edge handles.

mod elements;
mod graph;
mod list;
mod zone;

pub use elements::{
    ArcId, BasicEltId, BisectorId, EdgeId, MatArc, MatBasicElt, MatBisector, MatEdge, MatNode,
    MatSide, NodeId,
};
pub use graph::MatGraph;
pub use list::MatTList;
pub use zone::MatZone;
