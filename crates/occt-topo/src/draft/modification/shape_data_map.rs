//! `NCollection_IndexedDataMap<TopoDS_Shape, T, TopTools_ShapeMapHasher>`: the
//! insertion-ordered, 1-based indexed map `Draft_Modification` keeps its face,
//! edge and vertex records in (`Draft_Modification.hxx:203-207`).
//!
//! Only the operations the draft code performs are provided: `Add`
//! (replace-or-append), `Contains`, `FindFromKey`, `ChangeFromKey`,
//! `FindKey(theIndex)`, `FindFromIndex`, `ChangeFromIndex`, `RemoveKey`,
//! `Extent` and `Clear`.
//!
//! # Translation notes
//!
//! * Keys are compared the way `TopTools_ShapeMapHasher` compares them -
//!   [`shape_map_key`] - so orientation is ignored and two occurrences of one
//!   `TShape` under different locations are distinct keys
//!   (`TopoDS_Shape::IsSame`).
//! * Indices are 1-based like the reference. `RemoveKey`
//!   (`NCollection_IndexedDataMap::RemoveKey`) fills the hole with the last
//!   element, so the surviving indices shift exactly as they do in OCCT; every
//!   caller re-reads `FindKey(i)` after a removal.

use std::collections::HashMap;

use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

/// The 1-based indexed data map, insertion ordered.
pub struct ShapeDataMap<T> {
    items: Vec<(TopoShape, T)>,
    index: HashMap<ShapeMapKey, usize>,
}

impl<T> Default for ShapeDataMap<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> ShapeDataMap<T> {
    /// An empty map.
    pub fn new() -> Self {
        Self {
            items: Vec::new(),
            index: HashMap::new(),
        }
    }

    /// `M.Clear()`.
    pub fn clear(&mut self) {
        self.items.clear();
        self.index.clear();
    }

    /// `M.Extent()`.
    pub fn extent(&self) -> usize {
        self.items.len()
    }

    /// `M.Contains(theKey)`.
    pub fn contains(&self, s: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(s))
    }

    /// `M.Add(theKey, theItem)`: appends a new binding or replaces the value of
    /// an existing one, the index staying untouched.
    pub fn add(&mut self, s: &TopoShape, v: T) -> usize {
        let key = shape_map_key(s);
        match self.index.get(&key) {
            Some(&i) => {
                self.items[i].1 = v;
                i + 1
            }
            None => {
                self.items.push((s.clone(), v));
                let i = self.items.len() - 1;
                self.index.insert(key, i);
                i + 1
            }
        }
    }

    /// `M.FindFromKey(theKey)`.
    pub fn find_from_key(&self, s: &TopoShape) -> Option<&T> {
        self.index.get(&shape_map_key(s)).map(|&i| &self.items[i].1)
    }

    /// `M.ChangeFromKey(theKey)`.
    pub fn change_from_key(&mut self, s: &TopoShape) -> Option<&mut T> {
        match self.index.get(&shape_map_key(s)) {
            Some(&i) => Some(&mut self.items[i].1),
            None => None,
        }
    }

    /// `M.FindKey(theIndex)`, 1-based.
    pub fn find_key(&self, idx: usize) -> &TopoShape {
        &self.items[idx - 1].0
    }

    /// `M.FindFromIndex(theIndex)`, 1-based.
    pub fn find_from_index(&self, idx: usize) -> &T {
        &self.items[idx - 1].1
    }

    /// `M.ChangeFromIndex(theIndex)`, 1-based.
    pub fn change_from_index(&mut self, idx: usize) -> &mut T {
        &mut self.items[idx - 1].1
    }

    /// `M.RemoveKey(theKey)`.
    pub fn remove_key(&mut self, s: &TopoShape) {
        let Some(i) = self.index.remove(&shape_map_key(s)) else {
            return;
        };
        self.items.swap_remove(i);
        if i < self.items.len() {
            // `swap_remove` moved the last element into slot `i`.
            self.index.insert(shape_map_key(&self.items[i].0), i);
        }
    }

    /// The (key, value) pairs in insertion order.
    ///
    /// `NCollection_IndexedDataMap` exposes the same `begin()` / `end()`
    /// traversal; no caller of the draft port walks it yet.
    #[allow(dead_code)]
    pub fn iter(&self) -> impl Iterator<Item = &(TopoShape, T)> {
        self.items.iter()
    }
}
