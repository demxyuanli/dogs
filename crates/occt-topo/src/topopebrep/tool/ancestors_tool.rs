//! `TopOpeBRepTool_AncestorsTool` — the ancestor map of one sub-shape type.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_AncestorsTool.cxx:21-31`
//! and `.hxx:24-52` (OCCT 8.0.0).
//!
//! `MakeAncestors` is a one-line forward of `TopExp::MapShapesAndAncestors`
//! (`cxx:30`) that exists so the callers stay inside `TopOpeBRepTool`; the port
//! keeps the same forward onto
//! [`map_shapes_and_ancestors_into`](crate::topo_tools::map_shapes_and_ancestors_into),
//! which is `TopExp.cxx:80-120` with the same append-never-reset in/out map.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::ShapeMapKey;

/// `TopOpeBRepTool_AncestorsTool::MakeAncestors(S, TS, TA, M)` (`cxx:23-31`).
///
/// Appends, for every `TS` sub-shape of `S`, its `TA` ancestors to `m`
/// (`TopExp::MapShapesAndAncestors`, `TopExp.cxx:80-120`); the map is never
/// reset, matching the OCCT in/out contract.
pub fn make_ancestors(
    s: &TopoShape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) {
    crate::topo_tools::map_shapes_and_ancestors_into(s, ts, ta, m);
}
