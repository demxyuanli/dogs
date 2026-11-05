//! Insertion-ordered `TopExp::MapShapesAndAncestors`
//! (`TopExp.cxx:80-120`).
//!
//! `GetSpineAndProfile` (`BRepFill_AdvancedEvolved.cxx:110-149`, `:206-220`,
//! `:287-297`) walks the `NCollection_IndexedDataMap` **in insertion order**
//! and stops at the first corner whose tangents are not parallel, so the map
//! order is part of the answer. [`crate::topo_tools::map_shapes_and_ancestors`]
//! is a `HashMap` and cannot serve there.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, explorer_shapes_avoiding, shape_map_key, ShapeMapKey};

/// One entry of the `IndexedDataMap`: the sub-shape and the ancestors found
/// for it, in the order they were met.
pub(super) struct AncestorEntry {
    pub(super) key: TopoShape,
    pub(super) ancestors: Vec<TopoShape>,
}

/// `TopExp::MapShapesAndAncestors(S, T1, T2, M)` (`TopExp.cxx:80-120`) with
/// the `IndexedDataMap` iteration order: an entry is created the first time a
/// `T1` sub-shape is reached (walking the `T2` ancestors in explorer order),
/// and the trailing `TopExp_Explorer(S, T1)` pass (`cxx:110-119`) only
/// appends the `T1` sub-shapes that sit under no `T2` ancestor.
pub(super) fn ordered_map_shapes_and_ancestors(
    shape: &TopoShape,
    ts: ShapeType,
    ta: ShapeType,
) -> Vec<AncestorEntry> {
    let mut out: Vec<AncestorEntry> = Vec::new();
    let mut index: HashMap<ShapeMapKey, usize> = HashMap::new();
    for anc in explorer_shapes(shape, ta) {
        for sub in explorer_shapes(&anc, ts) {
            let key = shape_map_key(&sub);
            match index.get(&key) {
                Some(&i) => out[i].ancestors.push(anc.clone()),
                None => {
                    index.insert(key, out.len());
                    out.push(AncestorEntry {
                        key: sub,
                        ancestors: vec![anc.clone()],
                    });
                }
            }
        }
    }
    // `TopExp.cxx:110-119`: sub-shapes under no ancestor get an empty list.
    for sub in explorer_shapes_avoiding(shape, ts, ta) {
        let key = shape_map_key(&sub);
        if !index.contains_key(&key) {
            index.insert(key, out.len());
            out.push(AncestorEntry {
                key: sub,
                ancestors: Vec::new(),
            });
        }
    }
    out
}
