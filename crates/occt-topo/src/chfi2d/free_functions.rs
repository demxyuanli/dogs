//! The two package-level statics of `ChFi2d` (`ChFi2d.cxx:27-92`).
//!
//! | OCCT | lines |
//! |---|---|
//! | `ChFi2d::CommonVertex(E1, E2, V)` | `cxx:27-45` |
//! | `ChFi2d::FindConnectedEdges(F, V, E1, E2)` | `cxx:49-92` |
//!
//! Translation notes:
//! - `TopExp::Vertices(E1, first, last)` is
//!   [`crate::topo_tools::edge_vertices`], and the `.IsSame` chain of
//!   `cxx:34-44` is exactly `TopExp::CommonVertex` (`TopExp.cxx:323-340`),
//!   which the port already has as [`crate::topo_tools::common_vertex`]; the
//!   body is kept as its own function because `ChFi2d` re-exports it and its
//!   two callers (`ChFi2d_Builder`, `ChFi2d_FilletAlgo`) name it,
//! - `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//!   TopTools_ShapeMapHasher>` keyed by `IsSame` is the
//!   `HashMap<ShapeMapKey, Vec<TopoShape>>` of
//!   [`crate::topo_tools::map_shapes_and_ancestors`],
//! - `TopoDS::Edge(iterator.Value())` (`cxx:66`, `:75`) cannot fail: the map is
//!   built with `TA = TopAbs_EDGE`, so every list entry is an edge. Should a
//!   non-edge ever be stored there the port answers `ConnexionError` instead of
//!   OCCT's `Standard_TypeMismatch` throw; the map cannot produce that.

use crate::abs::ShapeType;
use crate::chfi2d::ConstructionError;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, map_shapes_and_ancestors, shape_map_key};

/// `ChFi2d::CommonVertex(E1, E2, V)` (`cxx:27-45`): the endpoint of `E1` that
/// is also an endpoint of `E2`.
///
/// The vertex handles are compared with `IsSame` in the order `first1` against
/// `first2` / `last2`, then `last1` against the same pair (`cxx:34-43`), which
/// is what the two `TopExp::Vertices` calls followed by the four tests do.
pub fn common_vertex(e1: &Edge, e2: &Edge) -> Option<TopoShape> {
    // `TopExp::Vertices(E1, firstVertex1, lastVertex1)` (`cxx:29-32`).
    let (first1, last1) = edge_vertices(e1);
    let (first2, last2) = edge_vertices(e2);
    let (first1, last1) = (first1?, last1?);
    let (first2, last2) = (first2?, last2?);
    // `if (firstVertex1.IsSame(firstVertex2) || firstVertex1.IsSame(lastVertex2))`
    // (`cxx:34-38`).
    if crate::topo_tools::is_same(&first1.0, &first2.0)
        || crate::topo_tools::is_same(&first1.0, &last2.0)
    {
        return Some(first1.0);
    }
    // `if (lastVertex1.IsSame(firstVertex2) || lastVertex1.IsSame(lastVertex2))`
    // (`cxx:39-43`).
    if crate::topo_tools::is_same(&last1.0, &first2.0)
        || crate::topo_tools::is_same(&last1.0, &last2.0)
    {
        return Some(last1.0);
    }
    // `return false;` (`cxx:44`).
    None
}

/// `ChFi2d::FindConnectedEdges(F, V, E1, E2)` (`cxx:49-92`): the two edges of
/// `F` that `V` is an endpoint of.
///
/// `Ok` is `ChFi2d_IsDone`; every `Err` is the single `ChFi2d_ConnexionError`
/// answer of the body — the vertex is absent from the map (`cxx:86-89`) or the
/// ancestor list does not hold exactly two edges (`cxx:57-84`).
pub fn find_connected_edges(face: &Face, v: &TopoShape) -> Result<(Edge, Edge), ConstructionError> {
    // `TopExp::MapShapesAndAncestors(F, TopAbs_VERTEX, TopAbs_EDGE, vertexMap)`
    // (`cxx:54-55`).
    let vertex_map = map_shapes_and_ancestors(&face.0, ShapeType::Vertex, ShapeType::Edge);
    // `if (vertexMap.Contains(V))` (`cxx:57`).
    let Some(edges) = vertex_map.get(&shape_map_key(v)) else {
        return Err(ConstructionError::ConnexionError);
    };
    // `iterator(vertexMap.FindFromKey(V)); if (iterator.More()) E1 = ...`
    // (`cxx:59-70`).
    let Some(e1) = edges.first() else {
        return Err(ConstructionError::ConnexionError);
    };
    let e1 = Edge::wrap(e1.clone()).ok_or(ConstructionError::ConnexionError)?;
    // `if (iterator.More()) E2 = ...` (`cxx:71-80`).
    let Some(e2) = edges.get(1) else {
        return Err(ConstructionError::ConnexionError);
    };
    let e2 = Edge::wrap(e2.clone()).ok_or(ConstructionError::ConnexionError)?;
    // `if (iterator.More()) return ChFi2d_ConnexionError;` (`cxx:82-85`).
    if edges.len() > 2 {
        return Err(ConstructionError::ConnexionError);
    }
    // `return ChFi2d_IsDone;` (`cxx:91`).
    Ok((e1, e2))
}
