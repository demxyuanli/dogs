//! `NCollection_DataMap<TopoDS_Shape, BRepOffset_Offset>` `myMapSF`
//! (`BiTgte_Blend.hxx:206`): the offset built for every analysed face, edge and
//! vertex of the blend.
//!
//! The reference walks this map with an `NCollection_DataMap::Iterator`, which
//! yields hash-bucket order, and also asks for the key shape itself
//! (`It.Key()` in `BiTgte_Blend.cxx:2565`). The port's [`ShapeMapKey`] is only
//! the `TShape` pointer plus location, so the shape cannot be recovered from a
//! key; this map therefore keeps the key shapes alongside, in insertion order.
//! `Bind` on an existing key replaces the value in place, exactly like
//! `NCollection_DataMap::Bind`.

use std::collections::HashMap;

use crate::brep_offset::BRepOffsetOffset;
use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

/// `NCollection_DataMap<TopoDS_Shape, BRepOffset_Offset, TopTools_ShapeMapHasher>`.
#[derive(Debug, Clone, Default)]
pub(super) struct OffsetMap {
    /// `(key, value)` in first-insertion order.
    entries: Vec<(TopoShape, BRepOffsetOffset)>,
    /// `key` -> position in `entries`.
    index: HashMap<ShapeMapKey, usize>,
}

impl OffsetMap {
    pub(super) fn new() -> Self {
        Self::default()
    }

    /// `IsBound(S)`.
    pub(super) fn is_bound(&self, s: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(s))
    }

    /// `Find(S)` / `operator()(S)`: the bound value.
    pub(super) fn get(&self, s: &TopoShape) -> Option<&BRepOffsetOffset> {
        self.index.get(&shape_map_key(s)).map(|&i| &self.entries[i].1)
    }

    /// `Bind(S, OF)`: replaces the value of an already bound key in place.
    pub(super) fn bind(&mut self, s: &TopoShape, of: BRepOffsetOffset) {
        let key = shape_map_key(s);
        match self.index.get(&key) {
            Some(&i) => self.entries[i].1 = of,
            None => {
                self.index.insert(key, self.entries.len());
                self.entries.push((s.clone(), of));
            }
        }
    }

    /// `Iterator`: `(key, value)` pairs.
    pub(super) fn iter(&self) -> impl Iterator<Item = (&TopoShape, &BRepOffsetOffset)> {
        self.entries.iter().map(|(k, v)| (k, v))
    }
}
