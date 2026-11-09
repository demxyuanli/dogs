//! The shape map of `BRepFeat_Builder`.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Builder.hxx:123-125`
//! (`myShapes` / `myRemoved` are
//! `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>`).
//!
//! `TopTools_ShapeMapHasher` keys on `TShape` + `Location` (`IsSame`), which is
//! exactly [`crate::locope::shape_maps::ShapeSet`]; the builder reuses that
//! index instead of re-deriving one, so `Contains` / `Add` answer with the same
//! equality as the reference.

pub use crate::locope::shape_maps::ShapeSet;
