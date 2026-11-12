//! `ShapeUpgrade_UnifySameDomain::FillHistory`
//! (`ShapeUpgrade_UnifySameDomain.cxx:4475-4558`).

use std::collections::HashSet;

use crate::abs::ShapeType;
use crate::bop_hist::BopHistory;
use crate::shape::TopoShape;
use crate::topo_tools::{map_shapes, shape_map_key, ShapeMapKey};

const HISTORY_TYPES: [ShapeType; 4] = [
    ShapeType::Vertex,
    ShapeType::Edge,
    ShapeType::Face,
    ShapeType::Solid,
];

pub(crate) fn fill_history(
    initial_shape: &TopoShape,
    result_shape: &TopoShape,
    context_history: &BopHistory,
    history: &mut BopHistory,
) -> Result<(), String> {
    let input = map_shapes(initial_shape, &HISTORY_TYPES);
    let result = map_shapes(result_shape, &HISTORY_TYPES);

    let result_keys: HashSet<ShapeMapKey> = HISTORY_TYPES
        .iter()
        .flat_map(|kind| result.get(kind).into_iter().flatten())
        .map(shape_map_key)
        .collect();

    let mut operation = BopHistory::new();
    for kind in HISTORY_TYPES {
        let Some(shapes) = input.get(&kind) else {
            continue;
        };
        for source in shapes {
            if result_keys.contains(&shape_map_key(source)) {
                continue;
            }

            let mut removed = true;
            if let Some(images) = context_history.modified(source) {
                for image in images {
                    if result_keys.contains(&shape_map_key(image)) {
                        if !crate::topo_tools::is_same(source, image) {
                            operation.add_modified(source, image.clone())?;
                        }
                        removed = false;
                    }
                }
            }
            if removed {
                operation.add_removed(source.clone())?;
            }
        }
    }

    history.merge(&operation);
    Ok(())
}
