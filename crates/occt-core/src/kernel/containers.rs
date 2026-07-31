//! OCCT NCollection → Rust stdlib type aliases.
//!
//! OCCT's NCollection is a template-based container library (header-only).
//! Rust's stdlib provides equivalent containers with better ergonomics.
//!
//! ## Mapping table
//!
//! | OCCT type | Rust equivalent | Notes |
//! |-----------|----------------|-------|
//! | `NCollection_Array1<T>` | `Vec<T>` | 1-indexed in OCCT, 0-indexed in Rust. Offset with `[i-1]`. |
//! | `NCollection_Array2<T>` | Flat `Vec<T>` + `(rows, cols)` | Column-major in OCCT. |
//! | `NCollection_Sequence<T>` | `Vec<T>` | Indexed sequence with append/prepend. |
//! | `NCollection_List<T>` | `std::collections::LinkedList<T>` or `Vec<T>` | Doubly-linked list. Use Vec unless insertion-order matters. |
//! | `NCollection_Map<T>` | `std::collections::HashSet<T>` | Unordered set. |
//! | `NCollection_DataMap<K,V>` | `std::collections::HashMap<K,V>` | Hash map. OCCT uses DefaultHasher. |
//! | `NCollection_IndexedMap<T>` | `Vec<T>` + `HashMap<T, usize>` | Ordered set with fast lookup. |
//! | `NCollection_IndexedDataMap<K,V>` | `Vec<(K,V)>` + `HashMap<K, usize>` | Ordered map. |
//! | `NCollection_DoubleMap<K,V>` | Two `HashMap` | Bidirectional map. |
//! | `NCollection_FlatMap<K,V>` | `Vec<(K,V)>` (sorted) | Small-map optimization. |
//! | `NCollection_PackedMap` | `bitvec` crate or `Vec<u64>` | Bit-set for integer keys. |
//! | `NCollection_Shared<T>` | `Arc<T>` | Shared ownership. |
//! | `NCollection_Handle<T>` | `Arc<T>` | Intrusive shared ownership. |
//! | `NCollection_LocalArray<T,N>` | `[T; N]` | Stack-allocated fixed array. |
//! | `NCollection_Buffer<T>` | `Box<[T]>` | Owned heap buffer. |
//! | `NCollection_Vector<T>` | `Vec<T>` | Owning vector. |
//! | `NCollection_UBTree<T>` | `std::collections::BTreeSet<T>` | Unbalanced binary tree → balanced B-Tree. |
//! | `NCollection_CellFilter<T>` | Spatial hash or grid | Spatial acceleration structure. |
//! | `NCollection_KDTree<T>` | `kdtree` crate or custom | K-dimensional tree. |
//!
//! ## Migration notes
//!
//! - **1-based → 0-based**: All OCCT collections use 1-based indexing.
//!   Ported code must subtract 1 from all indices, or use `get(i-1)` accessors.
//! - **Growth**: `NCollection_Array1` is contiguous; use `Vec::push` / `Vec::resize`.
//! - **Iteration**: OCCT uses explicit iterators (`NCollection_Iterator`).
//!   Rust uses `for x in &vec` / `iter()` / `iter_mut()`.
//! - **Allocators**: OCCT supports custom allocators (`NCollection_IncAllocator`,
//!   `NCollection_HeapAllocator`). Rust uses the global allocator by default;
//!   custom allocators are opt-in via the `Allocator` API (nightly).
//!
//! Ponteil note: 95% of OCCT collections map to Vec/HashMap with zero code.
//! The remaining 5% (PackedMap, CellFilter, KDTree) only matter in specific
//! algorithm hotspots — port those on demand.

// Re-export stdlib types under OCCT-compatible names
pub type NArray1<T> = Vec<T>;
pub type NSequence<T> = Vec<T>;
pub type NMap<T> = std::collections::HashSet<T>;
pub type NDataMap<K, V> = std::collections::HashMap<K, V>;
pub type NList<T> = std::collections::LinkedList<T>;
pub type NSet<T> = std::collections::BTreeSet<T>;

use std::collections::HashMap;

/// Indexed map: Vec for insertion order + HashMap for O(1) lookup.
/// Replaces `NCollection_IndexedMap<T>`.
#[derive(Debug, Clone)]
pub struct IndexedMap<T: std::hash::Hash + Eq + Clone> {
    items: Vec<T>,
    index: HashMap<T, usize>,
}

impl<T: std::hash::Hash + Eq + Clone> IndexedMap<T> {
    pub fn new() -> Self { Self { items: Vec::new(), index: HashMap::new() } }

    pub fn add(&mut self, item: T) -> usize {
        if let Some(&i) = self.index.get(&item) { return i; }
        let i = self.items.len();
        self.items.push(item.clone());
        self.index.insert(item, i);
        i
    }

    pub fn contains(&self, item: &T) -> bool { self.index.contains_key(item) }
    pub fn find_index(&self, item: &T) -> Option<usize> { self.index.get(item).copied() }
    pub fn value(&self, i: usize) -> &T { &self.items[i] }
    pub fn len(&self) -> usize { self.items.len() }
    pub fn is_empty(&self) -> bool { self.items.is_empty() }
}

/// Indexed data map: Vec for insertion order + HashMap for O(1) key lookup.
/// Replaces `NCollection_IndexedDataMap<K,V>`.
#[derive(Debug, Clone)]
pub struct IndexedDataMap<K: std::hash::Hash + Eq + Clone, V: Clone> {
    entries: Vec<(K, V)>,
    index: HashMap<K, usize>,
}

impl<K: std::hash::Hash + Eq + Clone, V: Clone> IndexedDataMap<K, V> {
    pub fn new() -> Self { Self { entries: Vec::new(), index: HashMap::new() } }

    pub fn add(&mut self, key: K, value: V) -> usize {
        if let Some(&i) = self.index.get(&key) { return i; }
        let i = self.entries.len();
        self.entries.push((key.clone(), value));
        self.index.insert(key, i);
        i
    }

    pub fn bind(&mut self, key: K, value: V) -> usize { self.add(key, value) }

    pub fn find_from_key(&self, key: &K) -> Option<&V> {
        self.index.get(key).map(|&i| &self.entries[i].1)
    }

    pub fn find_key(&self, i: usize) -> &K { &self.entries[i].0 }
    pub fn find_from_index(&self, i: usize) -> &V { &self.entries[i].1 }
    pub fn len(&self) -> usize { self.entries.len() }
    pub fn is_empty(&self) -> bool { self.entries.is_empty() }
}

/// PackedMap — bit-set for integer keys 0..N.
/// Replaces `NCollection_PackedMap` for key-based operations.
#[derive(Debug, Clone)]
pub struct PackedMap {
    bits: Vec<u64>,
    count: usize,
}

impl PackedMap {
    pub fn new() -> Self { Self { bits: Vec::new(), count: 0 } }

    fn ensure(&mut self, key: usize) {
        let block = key / 64;
        if block >= self.bits.len() { self.bits.resize(block + 1, 0); }
    }

    pub fn add(&mut self, key: usize) -> bool {
        self.ensure(key);
        let (block, bit) = (key / 64, key % 64);
        let mask = 1u64 << bit;
        if self.bits[block] & mask != 0 { return false; }
        self.bits[block] |= mask;
        self.count += 1;
        true
    }

    pub fn remove(&mut self, key: usize) -> bool {
        if key / 64 >= self.bits.len() { return false; }
        let (block, bit) = (key / 64, key % 64);
        let mask = 1u64 << bit;
        if self.bits[block] & mask == 0 { return false; }
        self.bits[block] &= !mask;
        self.count -= 1;
        true
    }

    pub fn contains(&self, key: usize) -> bool {
        let (block, bit) = (key / 64, key % 64);
        block < self.bits.len() && (self.bits[block] & (1u64 << bit)) != 0
    }

    pub fn len(&self) -> usize { self.count }
    pub fn is_empty(&self) -> bool { self.count == 0 }
    pub fn clear(&mut self) { self.bits.clear(); self.count = 0; }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexed_map_basics() {
        let mut m = IndexedMap::new();
        assert_eq!(m.add("hello".to_string()), 0);
        assert_eq!(m.add("world".to_string()), 1);
        assert_eq!(m.add("hello".to_string()), 0); // duplicate
        assert_eq!(m.len(), 2);
        assert!(m.contains(&"world".to_string()));
        assert_eq!(m.value(0), "hello");
    }

    #[test]
    fn indexed_data_map_basics() {
        let mut m = IndexedDataMap::new();
        assert_eq!(m.add(1, "one"), 0);
        assert_eq!(m.add(2, "two"), 1);
        assert_eq!(m.find_from_key(&1), Some(&"one"));
        assert_eq!(m.find_key(1), &2);
        assert_eq!(m.find_from_index(0), &"one");
    }

    #[test]
    fn packed_map_basics() {
        let mut m = PackedMap::new();
        assert!(m.add(5));
        assert!(m.add(100));
        assert!(!m.add(5)); // duplicate
        assert!(m.contains(5));
        assert_eq!(m.len(), 2);
        assert!(m.remove(5));
        assert!(!m.contains(5));
        assert_eq!(m.len(), 1);
    }
}
