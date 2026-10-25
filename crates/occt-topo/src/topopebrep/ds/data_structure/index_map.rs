//! The `NCollection_IndexedDataMap<TopoDS_Shape, T, TopTools_ShapeMapHasher>`
//! and `NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>` shapes
//! used by `TopOpeBRepDS_DataStructure`.
//!
//! The key is [`ShapeMapKey`] (`TShape` + `Location`, i.e. `TopoDS_Shape::IsSame`),
//! the same identity `TopTools_ShapeMapHasher` tests. Insertion order is kept
//! because OCCT indexes the map 1..Extent and callers iterate
//! `FindKey(1..Extent)`.

use std::collections::HashMap;

use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

/// `NCollection_IndexedDataMap<TopoDS_Shape, T, TopTools_ShapeMapHasher>`.
#[derive(Debug, Clone)]
pub struct IndexedShapeData<T> {
    /// Insertion-ordered keys (`FindKey(i)`).
    keys: Vec<ShapeMapKey>,
    /// First-seen occurrence of each key.
    shapes: HashMap<ShapeMapKey, TopoShape>,
    /// `FindFromKey` / `ChangeFromKey` values.
    values: HashMap<ShapeMapKey, T>,
}

impl<T> Default for IndexedShapeData<T> {
    fn default() -> Self {
        Self {
            keys: Vec::new(),
            shapes: HashMap::new(),
            values: HashMap::new(),
        }
    }
}

impl<T> IndexedShapeData<T> {
    /// `NCollection_IndexedDataMap()` — empty.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Extent()`.
    pub fn extent(&self) -> i32 {
        self.keys.len() as i32
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// `FindIndex(S)`: 0 when `S` is not in the map.
    pub fn find_index(&self, s: &TopoShape) -> i32 {
        let key = shape_map_key(s);
        match self.keys.iter().position(|k| *k == key) {
            Some(i) => i as i32 + 1,
            None => 0,
        }
    }

    /// `Contains(S)`.
    pub fn contains(&self, s: &TopoShape) -> bool {
        self.values.contains_key(&shape_map_key(s))
    }

    /// `Add(S, V)`: appends when `S` is new, returns the 1-based index.
    pub fn add(&mut self, s: &TopoShape, v: T) -> i32 {
        let key = shape_map_key(s);
        if !self.values.contains_key(&key) {
            self.keys.push(key.clone());
            self.shapes.insert(key.clone(), s.clone());
            self.values.insert(key, v);
        }
        self.find_index(s)
    }

    /// `FindKey(I)`.
    pub fn find_key(&self, i: i32) -> Option<&TopoShape> {
        let key = self.keys.get((i - 1).max(0) as usize)?;
        self.shapes.get(key)
    }

    /// `FindFromIndex(I)`.
    pub fn find_from_index(&self, i: i32) -> Option<&T> {
        let key = self.keys.get((i - 1).max(0) as usize)?;
        self.values.get(key)
    }

    /// `ChangeFromIndex(I)`.
    pub fn change_from_index(&mut self, i: i32) -> Option<&mut T> {
        let key = self.keys.get((i - 1).max(0) as usize)?.clone();
        self.values.get_mut(&key)
    }

    /// `FindFromKey(S)`.
    pub fn find_from_key(&self, s: &TopoShape) -> Option<&T> {
        self.values.get(&shape_map_key(s))
    }

    /// `ChangeFromKey(S)`.
    pub fn change_from_key(&mut self, s: &TopoShape) -> Option<&mut T> {
        self.values.get_mut(&shape_map_key(s))
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.keys.clear();
        self.shapes.clear();
        self.values.clear();
    }
}

/// `NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>`.
#[derive(Debug, Clone, Default)]
pub struct IndexedShapeSet {
    /// Insertion-ordered keys (`FindKey(i)`).
    keys: Vec<ShapeMapKey>,
    /// First-seen occurrence of each key.
    shapes: HashMap<ShapeMapKey, TopoShape>,
}

impl IndexedShapeSet {
    /// `NCollection_IndexedMap()` — empty.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Extent()`.
    pub fn extent(&self) -> i32 {
        self.keys.len() as i32
    }

    /// `FindIndex(S)`: 0 when `S` is not in the map.
    pub fn find_index(&self, s: &TopoShape) -> i32 {
        let key = shape_map_key(s);
        match self.keys.iter().position(|k| *k == key) {
            Some(i) => i as i32 + 1,
            None => 0,
        }
    }

    /// `Contains(S)`.
    pub fn contains(&self, s: &TopoShape) -> bool {
        self.shapes.contains_key(&shape_map_key(s))
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// `Add(S)`: appends when `S` is new, returns the 1-based index.
    pub fn add(&mut self, s: &TopoShape) -> i32 {
        let key = shape_map_key(s);
        if !self.shapes.contains_key(&key) {
            self.keys.push(key.clone());
            self.shapes.insert(key, s.clone());
        }
        self.find_index(s)
    }

    /// `FindKey(I)`.
    pub fn find_key(&self, i: i32) -> Option<&TopoShape> {
        let key = self.keys.get((i - 1).max(0) as usize)?;
        self.shapes.get(key)
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.keys.clear();
        self.shapes.clear();
    }
}

