//! `ChFiDS_StripeMap` — a vertex to stripes map.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_StripeMap.hxx:33-71`,
//! `ChFiDS_StripeMap.lxx:19-29`, `ChFiDS_StripeMap.cxx:22-56` (OCCT 8.0.0).
//!
//! OCCT wraps `NCollection_IndexedDataMap<TopoDS_Vertex,
//! NCollection_List<handle(ChFiDS_Stripe)>, TopTools_ShapeMapHasher>`: keys keep
//! their insertion order and are addressed 1-based, and the hasher compares
//! with `TopoDS_Shape::IsSame` (`TopTools_ShapeMapHasher.hxx:30-38`). The port
//! keeps two parallel vectors and compares with [`crate::topo_tools::is_same`].
//!
//! `FindFromKey` / `FindFromIndex` raise `StdFail_NoSuchObject` in OCCT on a
//! missing key or an out-of-range index; the port answers `Option`.

use crate::shape::TopoShape;
use crate::topo_tools;

use super::stripe::StripeHandle;

/// Source: `ChFiDS_StripeMap.hxx:34-70`.
#[derive(Debug, Clone, Default)]
pub struct StripeMap {
    /// `mymap`'s keys, in insertion order; `Extent()` is their count.
    keys: Vec<TopoShape>,
    /// `mymap`'s items, parallel to `keys`.
    lists: Vec<Vec<StripeHandle>>,
}

impl StripeMap {
    /// `ChFiDS_StripeMap()` (`cxx:22`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Add(V, F)` (`cxx:26-35`): appends `F` to the stripe list of vertex `V`,
    /// creating the list on first sight of `V`.
    pub fn add(&mut self, v: &TopoShape, f: StripeHandle) {
        let index = match self.keys.iter().position(|k| topo_tools::is_same(k, v)) {
            Some(i) => i,
            None => {
                self.keys.push(v.clone());
                self.lists.push(Vec::new());
                self.keys.len() - 1
            }
        };
        self.lists[index].push(f);
    }

    /// `Extent()` (`lxx:19-22`).
    pub fn extent(&self) -> usize {
        self.keys.len()
    }

    /// `FindFromKey(V)` (`cxx:39-43`).
    pub fn find_from_key(&self, v: &TopoShape) -> Option<&[StripeHandle]> {
        let index = self.keys.iter().position(|k| topo_tools::is_same(k, v))?;
        self.lists.get(index).map(|l| l.as_slice())
    }

    /// `operator()(V)` (`hxx:47-51`).
    pub fn get_from_key(&self, v: &TopoShape) -> Option<&[StripeHandle]> {
        self.find_from_key(v)
    }

    /// `FindFromIndex(I)` (`cxx:46-50`), `I` 1-based.
    pub fn find_from_index(&self, i: usize) -> Option<&[StripeHandle]> {
        if i == 0 {
            return None;
        }
        self.lists.get(i - 1).map(|l| l.as_slice())
    }

    /// `operator()(I)` (`hxx:54-58`).
    pub fn get_from_index(&self, i: usize) -> Option<&[StripeHandle]> {
        self.find_from_index(i)
    }

    /// `FindKey(I)` (`lxx:24-28`), `I` 1-based.
    pub fn find_key(&self, i: usize) -> Option<&TopoShape> {
        if i == 0 {
            return None;
        }
        self.keys.get(i - 1)
    }

    /// `Clear()` (`cxx:52-55`).
    pub fn clear(&mut self) {
        self.keys.clear();
        self.lists.clear();
    }
}
