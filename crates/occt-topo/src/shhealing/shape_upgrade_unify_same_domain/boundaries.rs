//! Boundary collection from `ShapeUpgrade_UnifySameDomain.cxx:1301-1352`.

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};

/// `AddOrdinaryEdges`: toggles repeated edge identities, records dropped seam
/// occurrences, removes shared boundaries already present in `edges`, and
/// appends the remaining boundary occurrences in indexed-map order.
pub(super) fn add_ordinary_edges(
    edges: &mut Vec<TopoShape>,
    shape: &TopoShape,
    removed_edges: &mut Vec<TopoShape>,
) -> Option<usize> {
    let mut new_edges: Vec<(ShapeMapKey, TopoShape)> = Vec::new();
    for edge in explorer_shapes(shape, ShapeType::Edge) {
        let key = shape_map_key(&edge);
        if let Some(index) = new_edges
            .iter()
            .position(|(candidate, _)| *candidate == key)
        {
            new_edges.remove(index);
            removed_edges.push(edge);
        } else {
            new_edges.push((key, edge));
        }
    }

    let mut first_dropped = None;
    let mut index = 0;
    while index < edges.len() {
        let key = shape_map_key(&edges[index]);
        if let Some(new_index) = new_edges
            .iter()
            .position(|(candidate, _)| *candidate == key)
        {
            new_edges.remove(new_index);
            let dropped = edges.remove(index);
            removed_edges.push(dropped);
            first_dropped.get_or_insert(index);
        } else {
            index += 1;
        }
    }

    edges.extend(new_edges.into_iter().map(|(_, edge)| edge));
    first_dropped
}
