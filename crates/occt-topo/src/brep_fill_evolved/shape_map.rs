//! `NCollection_DataMap` keyed by `TopoDS_Shape::IsSame`, in insertion order.
//!
//! `BRepFill_Evolved` keeps `myMap` as
//! `NCollection_DataMap<TopoDS_Shape, NCollection_DataMap<TopoDS_Shape,
//! NCollection_List<TopoDS_Shape>>>` (`BRepFill_Evolved.hxx:193-197`) and walks
//! it with `Iterator::Initialize` / `More` / `Next` several times
//! (`BRepFill_Evolved.cxx:1097`, `:1195`, `:1493`, `:1888`, `:2014`, `:2481`).
//!
//! `NCollection_BaseMap` stores its pairs in a flat array indexed by insertion
//! order and skips the freed slots on iteration, so an iteration visits the
//! keys in the order they were first bound. `Vec<(key, value)>` reproduces that
//! order, and `Bind` overwrites in place
//! (`NCollection_DataMap.hxx:397-406`, "overrides value previously bound to
//! the key").

use crate::shape::TopoShape;
use crate::topo_tools::is_same;

/// `NCollection_DataMap<TopoDS_Shape, V, TopTools_ShapeMapHasher>`.
#[derive(Clone, Default)]
pub struct ShapeMap<V> {
    entries: Vec<(TopoShape, V)>,
}

impl<V> ShapeMap<V> {
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.entries.clear();
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    fn index_of(&self, key: &TopoShape) -> Option<usize> {
        self.entries.iter().position(|(k, _)| is_same(k, key))
    }

    /// `IsBound(Key)`.
    pub fn is_bound(&self, key: &TopoShape) -> bool {
        self.index_of(key).is_some()
    }

    /// `Find(Key)`; `None` is the reference's `IsBound` test failing, since the
    /// port has no null reference to answer.
    pub fn find(&self, key: &TopoShape) -> Option<&V> {
        self.index_of(key).map(|i| &self.entries[i].1)
    }

    /// `ChangeFind(Key)`.
    pub fn find_mut(&mut self, key: &TopoShape) -> Option<&mut V> {
        match self.index_of(key) {
            Some(i) => Some(&mut self.entries[i].1),
            None => None,
        }
    }

    /// `Bind(Key, Item)`: insert a fresh slot, or overwrite the bound one.
    pub fn bind(&mut self, key: &TopoShape, value: V) {
        match self.index_of(key) {
            Some(i) => self.entries[i].1 = value,
            None => self.entries.push((key.clone(), value)),
        }
    }

    /// `ChangeFind` after an `if (!IsBound(Key)) Bind(Key, Empty)`:
    /// the value slot for `key`, created from `empty` when the key is new.
    ///
    /// Keeping the pair together keeps the insertion order of the reference
    /// without a second `IsBound` walk at every call site.
    pub fn entry_or_bind(&mut self, key: &TopoShape, empty: impl FnOnce() -> V) -> &mut V {
        if let Some(i) = self.index_of(key) {
            return &mut self.entries[i].1;
        }
        self.entries.push((key.clone(), empty()));
        let last = self.entries.len() - 1;
        &mut self.entries[last].1
    }

    /// `Remove(Key)`.
    pub fn remove_key(&mut self, key: &TopoShape) {
        if let Some(i) = self.index_of(key) {
            self.entries.remove(i);
        }
    }

    /// `Iterator` order: the insertion order of the bound keys.
    pub fn iter(&self) -> impl Iterator<Item = (&TopoShape, &V)> {
        self.entries.iter().map(|(k, v)| (k, v))
    }

    pub fn keys(&self) -> impl Iterator<Item = &TopoShape> {
        self.entries.iter().map(|(k, _)| k)
    }
}

/// The two-level `myMap` of `BRepFill_Evolved.hxx:193-197`:
/// `NCollection_DataMap<Shape, NCollection_DataMap<Shape,
/// NCollection_List<Shape>>>`.
impl ShapeMap<ShapeMap<Vec<TopoShape>>> {
    /// The bound value list of `outer` / `inner`, created on demand. This is
    /// the `myMap(A)(B).Append(...)` shape of `cxx:1231`, `:1945`, `:2325`.
    pub fn list_mut(&mut self, outer: &TopoShape, inner: &TopoShape) -> &mut Vec<TopoShape> {
        let inner_map = self.entry_or_bind(outer, ShapeMap::<Vec<TopoShape>>::new);
        inner_map.entry_or_bind(inner, Vec::new)
    }

    /// `myMap(A)(B)` for a read: `None` when either level is unbound.
    pub fn list(&self, outer: &TopoShape, inner: &TopoShape) -> Option<&Vec<TopoShape>> {
        self.find(outer).and_then(|m| m.find(inner))
    }
}

/// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`.
pub type ShapeListMap = ShapeMap<Vec<TopoShape>>;
