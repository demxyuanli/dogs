//! `IMeshData::MapOfIntegerInteger` = `NCollection_DataMap<int, bool>`.
//!
//! Iteration order follows OCCT bucket walk (`NCollection_BaseMap::Iterator`):
//! bucket index `std::hash<int>(key) % NbBuckets + 1`, within a bucket newest
//! Bind first (list prepend). MSVC / libstdc++ / libc++ hash `int` as identity.
//! After `ReSize`, OCCT prepends while walking newest-first and therefore
//! reverses each chain (head becomes oldest); `resize` matches that.
//!
//! Cocircular UV lattice diagonals depend on this order through
//! `createTriangles` (`BRepMesh_Delaun.cxx:488`).

use std::collections::VecDeque;

/// Prime steps from `NCollection_Primes.cxx` `THE_PRIME_VECTOR`.
const PRIMES: [usize; 24] = [
    101, 1009, 2003, 5003, 10007, 20011, 37003, 57037, 65003, 100019, 209953, 472393, 995329,
    2359297, 4478977, 9437185, 17915905, 35831809, 71663617, 150994945, 301989889, 573308929,
    1019215873, 2038431745,
];

fn next_prime_for_map(n: usize) -> usize {
    let target = n.saturating_add(1);
    for &p in &PRIMES {
        if p >= target {
            return p;
        }
    }
    n.saturating_add(1)
}

/// `NCollection_DataMap<int, bool>` used as Delaunay cavity / loop-edge maps.
#[derive(Clone, Debug)]
pub struct IntBoolDataMap {
    /// Buckets `1..=nb_buckets`; index 0 unused (OCCT HashCode is 1-based).
    buckets: Vec<VecDeque<(i32, bool)>>,
    nb_buckets: usize,
    size: usize,
}

impl Default for IntBoolDataMap {
    fn default() -> Self {
        Self::new()
    }
}

impl IntBoolDataMap {
    /// `MapOfIntegerInteger(10, allocator)` — initial bucket hint before first Bind.
    pub fn new() -> Self {
        Self {
            buckets: Vec::new(),
            nb_buckets: 10,
            size: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.size
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    pub fn clear(&mut self) {
        for b in &mut self.buckets {
            b.clear();
        }
        self.size = 0;
    }

    pub fn contains_key(&self, key: &i32) -> bool {
        self.get(*key).is_some()
    }

    pub fn get(&self, key: i32) -> Option<bool> {
        if self.buckets.is_empty() {
            return None;
        }
        let b = self.bucket_of(key);
        self.buckets[b].iter().find(|(k, _)| *k == key).map(|(_, v)| *v)
    }

    /// `operator()` — panics if unbound (callers mirror OCCT's bound access).
    pub fn value(&self, key: i32) -> bool {
        self.get(key).expect("IntBoolDataMap: unbound key")
    }

    /// `Bind`. Returns `true` when the key was not already present.
    /// Existing keys keep their list position; the value is overwritten.
    pub fn bind(&mut self, key: i32, value: bool) -> bool {
        if self.resizable() {
            self.resize(self.size);
        }
        let b = self.bucket_of(key);
        if let Some(slot) = self.buckets[b].iter_mut().find(|(k, _)| *k == key) {
            slot.1 = value;
            return false;
        }
        self.buckets[b].push_front((key, value));
        self.size += 1;
        true
    }

    /// `std::map`/`BTreeMap::insert` alias used at call sites that seed the map.
    pub fn insert(&mut self, key: i32, value: bool) -> Option<bool> {
        let had = self.get(key);
        self.bind(key, value);
        had
    }

    /// `UnBind`.
    pub fn unbind(&mut self, key: i32) -> bool {
        if self.buckets.is_empty() {
            return false;
        }
        let b = self.bucket_of(key);
        if let Some(pos) = self.buckets[b].iter().position(|(k, _)| *k == key) {
            self.buckets[b].remove(pos);
            self.size -= 1;
            return true;
        }
        false
    }

    /// `BTreeMap::remove` alias.
    pub fn remove(&mut self, key: &i32) -> Option<bool> {
        let v = self.get(*key)?;
        self.unbind(*key);
        Some(v)
    }

    /// Keys in OCCT `Iterator` order (bucket 1..N, each chain newest-first).
    pub fn keys(&self) -> Vec<i32> {
        let mut out = Vec::with_capacity(self.size);
        if self.buckets.is_empty() {
            return out;
        }
        for b in 1..=self.nb_buckets {
            for (k, _) in &self.buckets[b] {
                out.push(*k);
            }
        }
        out
    }

    fn resizable(&self) -> bool {
        self.size == 0 || self.size > self.nb_buckets
    }

    fn resize(&mut self, extent: usize) {
        let new_buckets = next_prime_for_map(extent);
        if !self.buckets.is_empty() && new_buckets <= self.nb_buckets {
            return;
        }
        let mut fresh = vec![VecDeque::new(); new_buckets + 1];
        if !self.buckets.is_empty() {
            // `NCollection_DataMap::ReSize` (`hxx:371-384`): walk each old
            // chain newest-first and prepend into the new bucket, which
            // reverses the chain (head becomes oldest).
            for b in 1..=self.nb_buckets {
                for &(k, v) in &self.buckets[b] {
                    let nb = ((k as usize) % new_buckets) + 1;
                    fresh[nb].push_front((k, v));
                }
            }
        }
        self.buckets = fresh;
        self.nb_buckets = new_buckets;
    }

    fn bucket_of(&self, key: i32) -> usize {
        // `NCollection_DataMap::HashCode`: `std::hash<int>(key) % NbBuckets + 1`.
        ((key as usize) % self.nb_buckets) + 1
    }
}

/// `deleteTriangle` Bind helper — `true` when the edge was newly bound.
pub fn map_bind(map: &mut IntBoolDataMap, key: i32, value: bool) -> bool {
    map.bind(key, value)
}
