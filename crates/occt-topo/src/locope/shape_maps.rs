//! Shape-keyed containers shared by the `LocOpe` ports.
//!
//! `NCollection_DataMap<TopoDS_Shape, V>` and
//! `NCollection_IndexedDataMap<TopoDS_Shape, V>` key their entries on
//! `TopTools_ShapeMapHasher` (`TopTools_ShapeMapHasher.hxx:44-52`): the `TShape`
//! together with the `TopLoc_Location`, orientation ignored. That is exactly
//! [`crate::topo_tools::shape_map_key`].
//!
//! OCCT's containers keep the key shape itself (`Key()`, `FindKey(I)`), so these
//! ports store the key shape next to the value rather than only its hash key.
//!
//! Iteration order: OCCT walks the hash buckets, which
//! `NCollection_DataMap` does not specify; the ports walk insertion order, which
//! is deterministic and stable across runs.

use std::collections::HashSet;

use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

/// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>`: one value per key.
#[derive(Debug, Clone, Default)]
pub struct ShapeDataMap {
    entries: Vec<(TopoShape, TopoShape)>,
    index: std::collections::HashMap<ShapeMapKey, usize>,
}

impl ShapeDataMap {
    /// `NCollection_DataMap()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Bind(K, V)`: append, or replace the value already bound to `K`.
    pub fn bind(&mut self, key: &TopoShape, value: &TopoShape) {
        match self.index.get(&shape_map_key(key)) {
            Some(&i) => self.entries[i].1 = value.clone(),
            None => {
                self.entries.push((key.clone(), value.clone()));
                let i = self.entries.len() - 1;
                self.index.insert(shape_map_key(key), i);
            }
        }
    }

    /// `IsBound(K)`.
    pub fn is_bound(&self, key: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(key))
    }

    /// `Find(K)`; OCCT raises `StdFail_NotDone` when unbound, the port answers
    /// `None`.
    pub fn find(&self, key: &TopoShape) -> Option<&TopoShape> {
        self.index.get(&shape_map_key(key)).map(|&i| &self.entries[i].1)
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.index.clear();
    }

    /// `NCollection_DataMap::Iterator`: `(Key, Value)` in map order.
    pub fn iter(&self) -> impl Iterator<Item = &(TopoShape, TopoShape)> {
        self.entries.iter()
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// `NCollection_IndexedDataMap<TopoDS_Shape, V>`: append-ordered key/value
/// pairs with 1-based indices, as OCCT exposes them.
#[derive(Debug, Clone, Default)]
pub struct ShapeIndexedMap<V> {
    entries: Vec<(TopoShape, V)>,
    index: std::collections::HashMap<ShapeMapKey, usize>,
}

impl<V: Clone> ShapeIndexedMap<V> {
    /// `NCollection_IndexedDataMap()`.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            index: std::collections::HashMap::new(),
        }
    }

    /// `Add(K, V)` (`NCollection_IndexedDataMap::Add`): append when `K` is
    /// absent, otherwise answer the existing index unchanged.
    pub fn add(&mut self, key: &TopoShape, value: V) -> usize {
        if let Some(&i) = self.index.get(&shape_map_key(key)) {
            return i;
        }
        self.entries.push((key.clone(), value));
        let i = self.entries.len() - 1;
        self.index.insert(shape_map_key(key), i);
        i
    }

    /// `Contains(K)`.
    pub fn contains(&self, key: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(key))
    }

    /// `Extent()`.
    pub fn extent(&self) -> usize {
        self.entries.len()
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// `FindKey(I)` (`I` is 1-based, as in OCCT).
    pub fn find_key(&self, i: usize) -> &TopoShape {
        &self.entries[i - 1].0
    }

    /// `FindFromIndex(I)` / `operator()(I)` (`I` is 1-based).
    pub fn find_from_index(&self, i: usize) -> &V {
        &self.entries[i - 1].1
    }

    /// `FindFromKey(K)`.
    pub fn find_from_key(&self, key: &TopoShape) -> Option<&V> {
        self.index.get(&shape_map_key(key)).map(|&i| &self.entries[i].1)
    }

    /// `FindIndex(K)` (`NCollection_IndexedMap` / `IndexedDataMap`): the
    /// 1-based index, or `None` for the reference's `0`.
    pub fn find_index(&self, key: &TopoShape) -> Option<usize> {
        self.index.get(&shape_map_key(key)).map(|&i| i + 1)
    }

    /// `ChangeFromKey(K)` (`NCollection_IndexedDataMap`): mutable access to the
    /// value bound to `K`; `None` stands for the reference's raise on an
    /// unbound key.
    pub fn change_from_key(&mut self, key: &TopoShape) -> Option<&mut V> {
        let &i = self.index.get(&shape_map_key(key))?;
        Some(&mut self.entries[i].1)
    }

    /// `M(theIdx)` from an index (`NCollection_IndexedDataMap::operator()`).
    pub fn change_from_index(&mut self, i: usize) -> &mut V {
        &mut self.entries[i - 1].1
    }

    /// `RemoveKey(K)` (`NCollection_IndexedDataMap::RemoveKey`, which calls
    /// `RemoveFromIndex`): drop the entry and shift every later index down.
    pub fn remove_key(&mut self, key: &TopoShape) {
        let Some(&i) = self.index.get(&shape_map_key(key)) else {
            return;
        };
        self.entries.remove(i);
        self.index.clear();
        for (new_i, (k, _)) in self.entries.iter().enumerate() {
            self.index.insert(shape_map_key(k), new_i);
        }
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.index.clear();
    }
}

/// `NCollection_DataMap<TopoDS_Shape, V>`: one value per key, with the
/// reference's in-place `Bind` / `ChangeFind` semantics.
///
/// `change` is `NCollection_DataMap::operator()` (`hxx:715`), which is
/// `ChangeFind` and therefore *raises* when the key is unbound
/// (`hxx:706-712`) - the `LocOpe_SplitShape` call sites always `Bind` first
/// (`Put`, the `VerSecMap` pre-pass, `SectionsTimes`).
#[derive(Debug, Clone, Default)]
pub struct ShapeValueMap<V> {
    entries: Vec<(TopoShape, V)>,
    index: std::collections::HashMap<ShapeMapKey, usize>,
}

impl<V: Clone> ShapeValueMap<V> {
    /// `NCollection_DataMap()`.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
            index: std::collections::HashMap::new(),
        }
    }

    /// `Bind(K, V)`: replace the bound value, or append. Answers `true` when the
    /// key was newly bound, which is what `Bind`'s `bool` return reports.
    pub fn bind(&mut self, key: &TopoShape, value: V) -> bool {
        match self.index.get(&shape_map_key(key)) {
            Some(&i) => {
                self.entries[i].1 = value;
                false
            }
            None => {
                self.entries.push((key.clone(), value));
                let i = self.entries.len() - 1;
                self.index.insert(shape_map_key(key), i);
                true
            }
        }
    }

    /// `IsBound(K)`.
    pub fn is_bound(&self, key: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(key))
    }

    /// `Find(K)`.
    pub fn find(&self, key: &TopoShape) -> Option<&V> {
        self.index.get(&shape_map_key(key)).map(|&i| &self.entries[i].1)
    }

    /// `Find(K)` mutable, the reference's non-const `operator()` minus the
    /// raise: `None` where OCCT throws `Standard_NoSuchObject`.
    pub fn change(&mut self, key: &TopoShape) -> Option<&mut V> {
        match self.index.get(&shape_map_key(key)) {
            Some(&i) => Some(&mut self.entries[i].1),
            None => None,
        }
    }

    /// `UnBind(K)`.
    pub fn unbind(&mut self, key: &TopoShape) {
        let Some(&i) = self.index.get(&shape_map_key(key)) else {
            return;
        };
        self.entries.remove(i);
        self.index.clear();
        for (new_i, (k, _)) in self.entries.iter().enumerate() {
            self.index.insert(shape_map_key(k), new_i);
        }
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.index.clear();
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// `Extent()`.
    pub fn extent(&self) -> usize {
        self.entries.len()
    }

    /// `NCollection_DataMap::Iterator`: `(Key, Value)` in map order.
    pub fn iter(&self) -> impl Iterator<Item = &(TopoShape, V)> {
        self.entries.iter()
    }
}

/// `NCollection_DataMap<TopoDS_Shape, double>`: one scalar per key, used by
/// `LocOpe_WiresOnShape::BindAll` for the projected vertex parameters.
#[derive(Debug, Clone, Default)]
pub struct ShapeDoubleMap {
    entries: Vec<(TopoShape, f64)>,
    index: std::collections::HashMap<ShapeMapKey, usize>,
}

impl ShapeDoubleMap {
    /// `NCollection_DataMap()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Bind(K, V)`.
    pub fn bind(&mut self, key: &TopoShape, value: f64) {
        match self.index.get(&shape_map_key(key)) {
            Some(&i) => self.entries[i].1 = value,
            None => {
                self.entries.push((key.clone(), value));
                let i = self.entries.len() - 1;
                self.index.insert(shape_map_key(key), i);
            }
        }
    }

    /// `IsBound(K)`.
    pub fn is_bound(&self, key: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(key))
    }

    /// `Find(K)`; `None` where OCCT raises `StdFail_NotDone`.
    pub fn find(&self, key: &TopoShape) -> Option<f64> {
        self.index.get(&shape_map_key(key)).map(|&i| self.entries[i].1)
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.index.clear();
    }
}

/// `NCollection_Map<TopoDS_Shape>`: membership only.
#[derive(Debug, Clone, Default)]
pub struct ShapeSet {
    items: Vec<TopoShape>,
    index: HashSet<ShapeMapKey>,
}

impl ShapeSet {
    /// `NCollection_Map()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Add(S)`: answers `true` when the shape was newly inserted
    /// (`NCollection_Map::Add`'s return value).
    pub fn add(&mut self, shape: &TopoShape) -> bool {
        if self.index.insert(shape_map_key(shape)) {
            self.items.push(shape.clone());
            true
        } else {
            false
        }
    }

    /// `Contains(S)`.
    pub fn contains(&self, shape: &TopoShape) -> bool {
        self.index.contains(&shape_map_key(shape))
    }

    /// `Remove(S)` (`NCollection_Map::Remove`): drop the entry when present.
    pub fn remove(&mut self, shape: &TopoShape) {
        if self.index.remove(&shape_map_key(shape)) {
            self.items.retain(|s| shape_map_key(s) != shape_map_key(shape));
        }
    }

    /// `NCollection_Map::Iterator`: the keys in insertion order.
    pub fn iter(&self) -> &[TopoShape] {
        &self.items
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// `Extent()`.
    pub fn extent(&self) -> usize {
        self.items.len()
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.items.clear();
        self.index.clear();
    }
}
