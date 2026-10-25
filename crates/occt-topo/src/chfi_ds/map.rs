//! `ChFiDS_Map` — an encapsulation of
//! `NCollection_IndexedDataMapOfShapeListOfShape`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_Map.hxx:34-65`,
//! `ChFiDS_Map.cxx:24-51` (OCCT 8.0.0).
//!
//! `Fill` forwards to `TopExp::MapShapesAndAncestors` (`cxx:28-31`), which the
//! port already has as
//! [`crate::topopebrep::tool::fuse_edges::IndexedAncestors`]: it keeps the
//! insertion order of both the keys and each ancestor list, matching the
//! 1-based `FindFromIndex` that `ChFi3d` walks.

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topopebrep::tool::fuse_edges::IndexedAncestors;

/// Source: `ChFiDS_Map.hxx:35-65`.
#[derive(Debug, Clone, Default)]
pub struct ChFiDSMap {
    /// `myMap` — `IndexedDataMapOfShapeListOfShape`.
    my_map: IndexedAncestors,
}

impl ChFiDSMap {
    /// `ChFiDS_Map()` (`cxx:24`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Fill(S, T1, T2)` (`cxx:28-31`): sub-shapes of type `T1` as keys and the
    /// list of their `T2` ancestors as items.
    ///
    /// `TopExp::MapShapesAndAncestors` does **not** clear the map first
    /// (`TopExp.cxx:80-120`), so neither does the port; a caller that reuses a
    /// `ChFiDS_Map` must build a fresh one, exactly as in OCCT.
    pub fn fill(&mut self, s: &TopoShape, t1: ShapeType, t2: ShapeType) {
        let built = IndexedAncestors::map_shapes_and_ancestors(s, t1, t2);
        self.my_map = built;
    }

    /// `Contains(S)` (`cxx:34-37`).
    pub fn contains(&self, s: &TopoShape) -> bool {
        self.my_map.find_from_key(s).is_some()
    }

    /// `FindFromKey(S)` (`cxx:41-44`). OCCT raises `Standard_NoSuchObject` when
    /// the key is absent; the port answers an empty slice.
    pub fn find_from_key(&self, s: &TopoShape) -> &[TopoShape] {
        self.my_map.find_from_key(s).unwrap_or(&[])
    }

    /// `operator()(S)` (`hxx:51-54`).
    pub fn get_from_key(&self, s: &TopoShape) -> &[TopoShape] {
        self.find_from_key(s)
    }

    /// `FindFromIndex(I)` (`cxx:47-50`), `I` 1-based. OCCT raises when `I` is
    /// out of range; the port answers an empty slice.
    pub fn find_from_index(&self, i: usize) -> &[TopoShape] {
        let Some(key) = self.my_map.find_key(i) else {
            return &[];
        };
        self.my_map.find_from_key(key).unwrap_or(&[])
    }

    /// `operator()(I)` (`hxx:57`).
    pub fn get_from_index(&self, i: usize) -> &[TopoShape] {
        self.find_from_index(i)
    }

    /// `myMap.Extent()` — the key count.
    pub fn extent(&self) -> usize {
        self.my_map.extent()
    }

    /// `myMap.FindKey(I)`, `I` 1-based.
    pub fn find_key(&self, i: usize) -> Option<&TopoShape> {
        self.my_map.find_key(i)
    }
}
