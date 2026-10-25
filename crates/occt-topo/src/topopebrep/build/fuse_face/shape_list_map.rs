//! `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//! TopTools_ShapeMapHasher>` as `TopOpeBRepBuild_FuseFace.cxx` uses it.
//!
//! Every map of the class has that exact type (`cxx:216-219`, `:307-308`,
//! `:595`, `:846-849`, `:976-981`, `:1080-1081`), so one ordered container
//! serves all of them:
//!
//! * `Bind(K, L)` keeps the position of an existing key and replaces its list,
//!   like `NCollection_DataMap::Bind`.
//! * `Find(K)` / `ChangeFind(K)` / `IsBound(K)` hash with
//!   [`shape_map_key`] — `TopTools_ShapeMapHasher` (`TopTools_ShapeMapHasher.hxx:44-52`)
//!   keys on `TShape` + location and **ignores the orientation**, which is what
//!   makes two occurrences of one edge in opposite orientations land in the
//!   same list (the whole point of `GroupShape` / `GroupEdge`).
//! * The value is an `NCollection_List`: appends keep duplicates and order.
//!
//! Deviation: OCCT iterates a `DataMap` in bucket order; the port iterates in
//! insertion order, the way every other ported `NCollection_*Map` does. The
//! class builds every list it later walks in insertion order, so the result
//! does not depend on the map order except through the loops that OCCT itself
//! drives in map order.

use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

/// `DataMap<Shape, List<Shape>>` — insertion-ordered key/value store.
#[derive(Debug, Clone, Default)]
pub struct ShapeListMap {
    entries: Vec<(ShapeMapKey, TopoShape, Vec<TopoShape>)>,
}

impl ShapeListMap {
    /// Empty map (`NCollection_DataMap()`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// `Extent()`.
    pub fn extent(&self) -> usize {
        self.entries.len()
    }

    /// `Clear()`.
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// `IsBound(K)`.
    pub fn is_bound(&self, key: &TopoShape) -> bool {
        let k = shape_map_key(key);
        self.entries.iter().any(|(kk, _, _)| *kk == k)
    }

    /// `Find(K)`; `None` when unbound.
    pub fn find(&self, key: &TopoShape) -> Option<&[TopoShape]> {
        let k = shape_map_key(key);
        self.entries
            .iter()
            .find(|(kk, _, _)| *kk == k)
            .map(|(_, _, v)| v.as_slice())
    }

    /// `Bind(K, L)`.
    pub fn bind(&mut self, key: &TopoShape, list: Vec<TopoShape>) {
        let k = shape_map_key(key);
        if let Some(e) = self.entries.iter_mut().find(|(kk, _, _)| *kk == k) {
            e.2 = list;
        } else {
            self.entries.push((k, key.clone(), list));
        }
    }

    /// `ChangeFind(K)`. OCCT throws `StdNoSuchObject` on an unbound key; the
    /// port panics, which is the same "programming error" contract.
    pub fn change_find(&mut self, key: &TopoShape) -> &mut Vec<TopoShape> {
        let k = shape_map_key(key);
        self.entries
            .iter_mut()
            .find(|(kk, _, _)| *kk == k)
            .map(|e| &mut e.2)
            .expect("ShapeListMap::change_find: unbound key")
    }

    /// `ChangeFind(K).Append(S)`.
    pub fn append(&mut self, key: &TopoShape, s: &TopoShape) {
        self.change_find(key).push(s.clone());
    }

    /// The bound key occurrences, in insertion order (`Iterator::Key()`).
    pub fn keys(&self) -> impl Iterator<Item = &TopoShape> {
        self.entries.iter().map(|(_, k, _)| k)
    }

    /// `(Key(), Find(Key()))` pairs in insertion order.
    pub fn pairs(&self) -> impl Iterator<Item = (&TopoShape, &[TopoShape])> {
        self.entries.iter().map(|(_, k, v)| (k, v.as_slice()))
    }
}
