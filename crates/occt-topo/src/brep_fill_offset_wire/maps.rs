//! `NCollection_IndexedDataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>`
//! as used by `TrimEdge`'s `MapVV` (`BRepFill_OffsetWire.cxx:2386`).
//!
//! `Add` keeps the value already bound to a key and appends a slot only for a
//! new key (`NCollection_IndexedDataMap.hxx`, `addImpl`); the slot order is the
//! insertion order and `ChangeFromIndex(i)` is 1-based.

use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

/// `NCollection_IndexedDataMap<TopoDS_Shape, TopoDS_Shape>`.
#[derive(Clone, Default)]
pub struct ShapeToShapeMap {
    items: Vec<(ShapeMapKey, TopoShape, TopoShape)>,
}

impl ShapeToShapeMap {
    /// An empty map.
    pub fn new() -> Self {
        Self::default()
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// `Extent()`.
    pub fn extent(&self) -> usize {
        self.items.len()
    }

    /// `Contains(K)`.
    pub fn contains(&self, key: &TopoShape) -> bool {
        let k = shape_map_key(key);
        self.items.iter().any(|(existing, _, _)| *existing == k)
    }

    /// `FindFromKey(K)`; `None` where OCCT raises `Standard_NoSuchObject`.
    pub fn find_from_key(&self, key: &TopoShape) -> Option<&TopoShape> {
        let k = shape_map_key(key);
        self.items
            .iter()
            .find(|(existing, _, _)| *existing == k)
            .map(|(_, _, v)| v)
    }

    /// `ChangeFromKey(K)`; `None` for an unbound key.
    pub fn change_from_key(&mut self, key: &TopoShape) -> Option<&mut TopoShape> {
        let k = shape_map_key(key);
        self.items
            .iter_mut()
            .find(|(existing, _, _)| *existing == k)
            .map(|(_, _, v)| v)
    }

    /// `Add(K, T)`: append `key` unless it is already bound.
    pub fn add(&mut self, key: &TopoShape, value: TopoShape) {
        let k = shape_map_key(key);
        if !self.items.iter().any(|(existing, _, _)| *existing == k) {
            self.items.push((k, key.clone(), value));
        }
    }

    /// `FindKey(Index)` (`index` is 1-based).
    pub fn find_key(&self, index: usize) -> Option<&TopoShape> {
        self.items.get(index.checked_sub(1)?).map(|(_, k, _)| k)
    }

    /// `operator()(Index)` (`index` is 1-based).
    pub fn value(&self, index: usize) -> Option<&TopoShape> {
        self.items.get(index.checked_sub(1)?).map(|(_, _, v)| v)
    }

    /// `ChangeFromIndex(Index)` (`index` is 1-based).
    pub fn value_mut(&mut self, index: usize) -> Option<&mut TopoShape> {
        self.items.get_mut(index.checked_sub(1)?).map(|(_, _, v)| v)
    }
}
