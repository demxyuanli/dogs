//! `BOPAlgo_Tools::FillMap` / `MakeBlocks` for shape-keyed connexity.
//!
//! Source: `BOPAlgo_Tools.hxx:44-102` (the two templates) as used by
//! `BOPAlgo_Builder::FillSameDomainFaces` (`BOPAlgo_Builder_2.cxx:783`,
//! `:820`, `:826`).
//!
//! OCCT stores a bidirectional adjacency
//! `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`.
//! `FillMap(n1, n2)` appends `n2` to the list of `n1` and `n1` to the list
//! of `n2`. `MakeBlocks` then walks that map in insertion order: each key
//! that is not yet in the fence starts a chain, and the chain grows by
//! appending every neighbour that the fence has not seen.
//!
//! Identity is [`crate::bop_occt_util::shape_key`] (`TopTools_ShapeMapHasher`
//! on the TShape handle). The stored `TopoShape` values keep orientation
//! and location of the first insertion.

use std::collections::{HashMap, HashSet};

use crate::bop_occt_util::shape_key;
use crate::shape::TopoShape;

/// Indexed shape-to-neighbours map used as `aDMSLS` in FillSameDomainFaces.
///
/// Mirrors `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>`. Keys are visited in insertion order so
/// `MakeBlocks` reproduces OCCT chain order (`FindKey(1..Extent)`).
#[derive(Debug, Clone, Default)]
pub struct IndexedShapeAdj {
    /// Insertion-ordered TShape keys (`FindKey(i)`).
    keys: Vec<usize>,
    /// Neighbour keys of each shape key (`FindFromKey`).
    lists: HashMap<usize, Vec<usize>>,
    /// First-seen `TopoShape` for each key (`theMILI.FindKey` value).
    shapes: HashMap<usize, TopoShape>,
}

impl IndexedShapeAdj {
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of keys (`Extent`).
    pub fn extent(&self) -> usize {
        self.keys.len()
    }

    /// True when the map holds `s` (`IsBound` / `ChangeSeek` not null).
    pub fn contains(&self, s: &TopoShape) -> bool {
        self.lists.contains_key(&shape_key(s))
    }

    /// The stored shape for `key`, when present.
    pub fn shape_of(&self, key: usize) -> Option<&TopoShape> {
        self.shapes.get(&key)
    }

    /// Neighbour keys of `key`.
    pub fn neighbours(&self, key: usize) -> &[usize] {
        self.lists.get(&key).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Keys in insertion order (`FindKey(1..Extent)`).
    pub fn keys(&self) -> &[usize] {
        &self.keys
    }

    fn ensure_key(&mut self, s: &TopoShape) -> usize {
        let k = shape_key(s);
        if !self.lists.contains_key(&k) {
            self.keys.push(k);
            self.lists.insert(k, Vec::new());
            self.shapes.insert(k, s.clone());
        }
        k
    }

    /// `BOPAlgo_Tools::FillMap(n1, n2, theMILI, theAllocator)`.
    ///
    /// Appends `n2` onto the list of `n1` and `n1` onto the list of `n2`.
    /// Duplicate neighbour entries are kept, matching `NCollection_List::Append`
    /// (MakeBlocks fences by key so duplicates do not duplicate a chain
    /// member).
    pub fn fill_map(&mut self, n1: &TopoShape, n2: &TopoShape) {
        let k1 = self.ensure_key(n1);
        let k2 = self.ensure_key(n2);
        self.lists.get_mut(&k1).expect("k1").push(k2);
        self.lists.get_mut(&k2).expect("k2").push(k1);
    }

    /// `BOPAlgo_Tools::MakeBlocks(theMILI, theMBlocks, theAllocator)`.
    ///
    /// Each not-yet-seen key starts a chain. The chain is grown by walking
    /// the neighbour lists of every member already in the chain (the OCCT
    /// `aItLChain` iterator that keeps going as `aChain.Append` extends the
    /// list). Isolated keys become singleton blocks; FillSameDomainFaces
    /// only ever inserts pairs, so those blocks have size >= 2.
    pub fn make_blocks(&self) -> Vec<Vec<TopoShape>> {
        let mut fence: HashSet<usize> = HashSet::new();
        let mut blocks: Vec<Vec<TopoShape>> = Vec::new();
        for &n in &self.keys {
            if !fence.insert(n) {
                continue;
            }
            let mut chain_keys: Vec<usize> = vec![n];
            let mut i = 0;
            while i < chain_keys.len() {
                let n1 = chain_keys[i];
                if let Some(li) = self.lists.get(&n1) {
                    for &n2 in li {
                        if fence.insert(n2) {
                            chain_keys.push(n2);
                        }
                    }
                }
                i += 1;
            }
            let chain: Vec<TopoShape> = chain_keys
                .into_iter()
                .filter_map(|k| self.shapes.get(&k).cloned())
                .collect();
            blocks.push(chain);
        }
        blocks
    }

    /// Drop every entry (`aDMSLS.Clear()`).
    pub fn clear(&mut self) {
        self.keys.clear();
        self.lists.clear();
        self.shapes.clear();
    }
}

/// `BOPAlgo_Tools::FillMap` free function used at `_2.cxx:783` and `:820`.
pub fn fill_map_shapes(n1: &TopoShape, n2: &TopoShape, map: &mut IndexedShapeAdj) {
    map.fill_map(n1, n2);
}

/// `BOPAlgo_Tools::MakeBlocks` free function used at `_2.cxx:826`.
pub fn make_blocks_shapes(map: &IndexedShapeAdj) -> Vec<Vec<TopoShape>> {
    map.make_blocks()
}

/// Indexed integer-to-integer adjacency (`NCollection_IndexedDataMap<int, List<int>>`).
///
/// Same templates as the shape map; used by `BOPAlgo_PaveFiller` connexity of
/// vertex/edge indices. Insertion order is preserved.
#[derive(Debug, Clone, Default)]
pub struct IndexedIntAdj {
    keys: Vec<usize>,
    lists: HashMap<usize, Vec<usize>>,
}

impl IndexedIntAdj {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn extent(&self) -> usize {
        self.keys.len()
    }

    pub fn contains(&self, n: usize) -> bool {
        self.lists.contains_key(&n)
    }

    pub fn neighbours(&self, n: usize) -> &[usize] {
        self.lists.get(&n).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn keys(&self) -> &[usize] {
        &self.keys
    }

    fn ensure_key(&mut self, n: usize) {
        if !self.lists.contains_key(&n) {
            self.keys.push(n);
            self.lists.insert(n, Vec::new());
        }
    }

    /// `BOPAlgo_Tools::FillMap<int>(n1, n2, theMILI, theAllocator)`.
    pub fn fill_map(&mut self, n1: usize, n2: usize) {
        self.ensure_key(n1);
        self.ensure_key(n2);
        self.lists.get_mut(&n1).expect("n1").push(n2);
        self.lists.get_mut(&n2).expect("n2").push(n1);
    }

    /// `BOPAlgo_Tools::MakeBlocks` over integer keys.
    pub fn make_blocks(&self) -> Vec<Vec<usize>> {
        let mut fence: HashSet<usize> = HashSet::new();
        let mut blocks: Vec<Vec<usize>> = Vec::new();
        for &n in &self.keys {
            if !fence.insert(n) {
                continue;
            }
            let mut chain: Vec<usize> = vec![n];
            let mut i = 0;
            while i < chain.len() {
                let n1 = chain[i];
                if let Some(li) = self.lists.get(&n1) {
                    for &n2 in li {
                        if fence.insert(n2) {
                            chain.push(n2);
                        }
                    }
                }
                i += 1;
            }
            blocks.push(chain);
        }
        blocks
    }

    pub fn clear(&mut self) {
        self.keys.clear();
        self.lists.clear();
    }
}

/// `BOPAlgo_Tools::FillMap` for integer indices.
pub fn fill_map_int(n1: usize, n2: usize, map: &mut IndexedIntAdj) {
    map.fill_map(n1, n2);
}

/// `BOPAlgo_Tools::MakeBlocks` for integer indices.
pub fn make_blocks_int_indexed(map: &IndexedIntAdj) -> Vec<Vec<usize>> {
    map.make_blocks()
}
