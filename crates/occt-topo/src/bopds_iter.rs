//! BOPDS interference-pair iterators.
use std::collections::HashSet;

use occt_core::bnd::BndBox;

use crate::abs::ShapeType;
use crate::bopds::BopdsDS;
use crate::bopds_tools;

/// Extended interference pairs for the vertices of `the_indices` whose
/// tolerance was increased.
///
/// Port of `BOPDS_Iterator::IntersectExt` (`BOPDS_Iterator.cxx`): every vertex
/// of `the_indices` (using its same-domain representative's box) is tested
/// against *every* interfering shape of the DS (using the representative box
/// for map members, the own box otherwise). The resulting pairs are returned
/// bucketed by interference type (VV..ZZ) — `BOPDS_Iterator::Initialize` then
/// serves these buckets *instead of* the regular pair lists, so the repeated
/// V/V, V/E and V/F stages process only the pairs involving the increased
/// vertices.
///
/// `the_indices` holds *source* vertex indices (as built by
/// `BOPAlgo_PaveFiller::RepeatIntersection`). The box-overlap test is done
/// brute-force (mirroring [`BopdsIterator::prepare`]) rather than with the
/// OCCT BVH tree — the semantics are identical.
pub fn intersect_ext_pairs(ds: &BopdsDS, the_indices: &HashSet<usize>) -> Vec<Vec<(usize, usize)>> {
    let mut buckets: Vec<Vec<(usize, usize)>> = vec![Vec::new(); NB_INTERF_TYPES];
    if the_indices.is_empty() {
        return buckets;
    }
    let n = ds.nb_source_shapes();
    // Resolved box of a shape: the same-domain representative's box for the
    // map members (OCCT `ShapeInfo(nVSD).Box()`), the own box otherwise
    // (OCCT `aSI.Box()` for the tree).
    let box_of = |k: usize, in_map: bool| -> BndBox {
        if in_map {
            ds.boxes[ds.get_same_domain_index(k)].clone()
        } else {
            ds.boxes[k].clone()
        }
    };
    let mut fence: HashSet<(usize, usize)> = HashSet::new();
    for i in 0..n {
        let si = &ds.shape_infos[i];
        // `IntersectExt` skips shapes without a boundary representation and
        // solids; `has_brep` is exactly the effective set.
        if !si.is_interfering() || si.kind == ShapeType::Solid {
            continue;
        }
        if !the_indices.contains(&i) {
            // Only map members are queries (every shape still lives in the tree).
            continue;
        }
        let box_i = box_of(i, true);
        for j in 0..n {
            if i == j {
                continue;
            }
            let sj = &ds.shape_infos[j];
            if !sj.is_interfering() || sj.kind == ShapeType::Solid {
                continue;
            }
            if ds.rank(i) == ds.rank(j) {
                continue;
            }
            let box_j = box_of(j, the_indices.contains(&j));
            if box_i.is_out_box(&box_j) {
                continue;
            }
            // Avoid interfering a shape with its own sub-shapes.
            let t1 = bopds_tools::type_to_integer(si.kind);
            let t2 = bopds_tools::type_to_integer(sj.kind);
            if (t1 < t2 && si.has_subshape(j)) || (t1 > t2 && sj.has_subshape(i)) {
                continue;
            }
            let ix = bopds_tools::type_to_integer2(si.kind, sj.kind);
            if (0..NB_INTERF_TYPES as i32).contains(&ix) {
                let (a, b) = (i.min(j), i.max(j));
                if fence.insert((a, b)) {
                    buckets[ix as usize].push((i, j));
                }
            }
        }
    }
    buckets
}

// ---------------------------------------------------------------------------
// BopdsIterator
// ---------------------------------------------------------------------------

/// Iterator over pairs of sub-shapes of the arguments whose bounding boxes
/// interfere. Source: `BOPDS_Iterator.hxx/.cxx`.
///
/// The regular iterator only reports *inter-argument* pairs (sub-shapes of
/// different arguments). Use [`BopdsIteratorSI`] for self-intersections.
#[derive(Debug)]
pub struct BopdsIterator<'a> {
    ds: Option<&'a BopdsDS>,
    /// Pair lists indexed by interference type (VV..ZZ).
    lists: Vec<Vec<(usize, usize)>>,
    /// The active (initialized) pair list.
    current: Vec<(usize, usize)>,
    /// Position in `current`.
    pos: usize,
    /// When true, self-intersection mode (no inter-argument skip, solids
    /// participate) — used by [`BopdsIteratorSI`].
    self_intersection: bool,
}

const NB_INTERF_TYPES: usize = 10;

impl<'a> BopdsIterator<'a> {
    /// Empty constructor (regular, inter-argument mode).
    pub fn new() -> Self {
        Self {
            ds: None,
            lists: vec![Vec::new(); NB_INTERF_TYPES],
            current: Vec::new(),
            pos: 0,
            self_intersection: false,
        }
    }

    /// Empty constructor in self-intersection mode.
    pub fn new_self_intersection() -> Self {
        Self {
            ds: None,
            lists: vec![Vec::new(); NB_INTERF_TYPES],
            current: Vec::new(),
            pos: 0,
            self_intersection: true,
        }
    }

    /// Sets the data structure to process.
    pub fn set_ds(&mut self, ds: &'a BopdsDS) {
        self.ds = Some(ds);
    }

    /// Returns the data structure.
    pub fn ds(&self) -> &'a BopdsDS {
        self.ds.expect("BopdsIterator: DS not set")
    }

    /// Initializes the iterator for the given pair of shape types.
    pub fn initialize(&mut self, t1: ShapeType, t2: ShapeType) {
        self.pos = 0;
        let ix = bopds_tools::type_to_integer2(t1, t2);
        if ix >= 0 {
            let mut list = self.lists[ix as usize].clone();
            list.sort_unstable();
            self.current = list;
        } else {
            self.current.clear();
        }
    }

    /// True if there are still pairs to iterate.
    pub fn more(&self) -> bool {
        self.pos < self.current.len()
    }

    /// Moves the iteration ahead.
    pub fn next(&mut self) {
        self.pos += 1;
    }

    /// Returns the indices (DS) of the interfering shapes, ordered so the
    /// shape of the higher type comes first.
    pub fn value(&self) -> (usize, usize) {
        let (i, j) = self.current[self.pos];
        let ds = self.ds();
        let t1 = bopds_tools::type_to_integer(ds.shape_infos[i].kind);
        let t2 = bopds_tools::type_to_integer(ds.shape_infos[j].kind);
        if t1 < t2 {
            (j, i)
        } else {
            (i, j)
        }
    }

    /// Performs the interference (bounding-box) computation and prepares the
    /// pair lists to be used.
    pub fn prepare(&mut self) {
        for list in &mut self.lists {
            list.clear();
        }
        self.current.clear();
        self.pos = 0;
        let Some(ds) = self.ds else { return };
        let n = ds.nb_source_shapes();
        for i in 0..n {
            let si = &ds.shape_infos[i];
            let ok = if self.self_intersection {
                bopds_tools::is_interfering(si.kind)
            } else {
                bopds_tools::has_brep(si.kind)
            };
            if !ok {
                continue;
            }
            for j in (i + 1)..n {
                let sj = &ds.shape_infos[j];
                let okj = if self.self_intersection {
                    bopds_tools::is_interfering(sj.kind)
                } else {
                    bopds_tools::has_brep(sj.kind)
                };
                if !okj {
                    continue;
                }
                // Inter-argument mode: skip pairs within the same argument.
                if !self.self_intersection && ds.rank(i) == ds.rank(j) {
                    continue;
                }
                // Bounding-box rejection.
                if ds.boxes[i].is_out_box(&ds.boxes[j]) {
                    continue;
                }
                // Avoid interfering a shape with its own sub-shapes.
                let t1 = bopds_tools::type_to_integer(si.kind);
                let t2 = bopds_tools::type_to_integer(sj.kind);
                if (t1 < t2 && si.has_subshape(j)) || (t1 > t2 && sj.has_subshape(i)) {
                    continue;
                }
                let ix = bopds_tools::type_to_integer2(si.kind, sj.kind);
                if (0..NB_INTERF_TYPES as i32).contains(&ix) {
                    self.lists[ix as usize].push((i, j));
                }
            }
        }
    }

    /// Returns the number of pairs in the currently initialized list.
    pub fn expected_length(&self) -> usize {
        self.current.len()
    }

    /// Returns the block length (half the expected length, at least 1).
    pub fn block_length(&self) -> usize {
        let n = self.expected_length();
        if n <= 1 {
            1
        } else {
            (0.5 * n as f64) as usize
        }
    }

    /// Clears the pair lists for interference types above `level`.
    /// Level semantics: 0 — only V/V; 1 — +V/E; 2 — +E/E; 3 — +V/F;
    /// 4 — +E/F; anything else — all interferences.
    pub fn update_by_level_of_check(&mut self, level: i32) {
        for k in (level + 1).max(0) as usize..self.lists.len() {
            self.lists[k].clear();
        }
    }
}

impl Default for BopdsIterator<'_> {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BopdsSubIterator
// ---------------------------------------------------------------------------

/// Iterator over pairs of interfering sub-shapes of two given sub-sets of DS
/// indices. Source: `BOPDS_SubIterator.hxx/.cxx`.
#[derive(Debug)]
pub struct BopdsSubIterator<'a> {
    ds: Option<&'a BopdsDS>,
    /// First set of indices.
    subset1: Vec<usize>,
    /// Second set of indices.
    subset2: Vec<usize>,
    /// The interfering pairs found by [`BopdsSubIterator::prepare`].
    pairs: Vec<(usize, usize)>,
    /// Position in `pairs`.
    pos: usize,
}

impl<'a> BopdsSubIterator<'a> {
    /// Empty constructor.
    pub fn new() -> Self {
        Self { ds: None, subset1: Vec::new(), subset2: Vec::new(), pairs: Vec::new(), pos: 0 }
    }

    /// Sets the data structure to process.
    pub fn set_ds(&mut self, ds: &'a BopdsDS) {
        self.ds = Some(ds);
    }

    /// Returns the data structure.
    pub fn ds(&self) -> &'a BopdsDS {
        self.ds.expect("BopdsSubIterator: DS not set")
    }

    /// Sets the first set of indices to process.
    pub fn set_subset1(&mut self, subset: Vec<usize>) {
        self.subset1 = subset;
    }

    /// Returns the first set of indices.
    pub fn subset1(&self) -> &[usize] {
        &self.subset1
    }

    /// Sets the second set of indices to process.
    pub fn set_subset2(&mut self, subset: Vec<usize>) {
        self.subset2 = subset;
    }

    /// Returns the second set of indices.
    pub fn subset2(&self) -> &[usize] {
        &self.subset2
    }

    /// Initializes the iterator (sorts the pairs).
    pub fn initialize(&mut self) {
        self.pairs.sort_unstable();
        self.pos = 0;
    }

    /// True if there are more pairs of interfering shapes.
    pub fn more(&self) -> bool {
        self.pos < self.pairs.len()
    }

    /// Moves the iteration ahead.
    pub fn next(&mut self) {
        self.pos += 1;
    }

    /// Returns the indices (DS) of the interfering shapes.
    pub fn value(&self) -> (usize, usize) {
        let (i, j) = self.pairs[self.pos];
        let ds = self.ds();
        let t1 = bopds_tools::type_to_integer(ds.shape_infos[i].kind);
        let t2 = bopds_tools::type_to_integer(ds.shape_infos[j].kind);
        if t1 < t2 {
            (j, i)
        } else {
            (i, j)
        }
    }

    /// Performs the intersection of bounding boxes of the two sub-sets and
    /// prepares the results.
    pub fn prepare(&mut self) {
        self.pairs.clear();
        let Some(ds) = self.ds else { return };
        if self.subset1.is_empty() || self.subset2.is_empty() {
            return;
        }
        let mut fence = HashSet::new();
        for &i in &self.subset1 {
            if i >= ds.shape_infos.len() {
                continue;
            }
            let si = &ds.shape_infos[i];
            for &j in &self.subset2 {
                if j >= ds.shape_infos.len() || i == j {
                    continue;
                }
                let (a, b) = (i.min(j), i.max(j));
                if !fence.insert((a, b)) {
                    continue;
                }
                if ds.boxes[i].is_out_box(&ds.boxes[j]) {
                    continue;
                }
                let sj = &ds.shape_infos[j];
                let t1 = bopds_tools::type_to_integer(si.kind);
                let t2 = bopds_tools::type_to_integer(sj.kind);
                if (t1 < t2 && si.has_subshape(j)) || (t1 > t2 && sj.has_subshape(i)) {
                    continue;
                }
                self.pairs.push((a, b));
            }
        }
        self.pairs.sort_unstable();
        self.pos = 0;
    }

    /// Returns the number of interfering pairs.
    pub fn expected_length(&self) -> usize {
        self.pairs.len()
    }
}

impl Default for BopdsSubIterator<'_> {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BopdsIteratorSI
// ---------------------------------------------------------------------------

/// Iterator computing *self-intersections* between the sub-shapes of each
/// argument. Source: `BOPDS_IteratorSI.hxx/.cxx`.
///
/// Unlike [`BopdsIterator`], it considers pairs inside the same argument and
/// lets solids participate, so it can detect an argument that intersects
/// itself (e.g. two overlapping boxes).
#[derive(Debug)]
pub struct BopdsIteratorSI<'a> {
    inner: BopdsIterator<'a>,
}

impl<'a> BopdsIteratorSI<'a> {
    /// Empty constructor (self-intersection mode).
    pub fn new() -> Self {
        Self { inner: BopdsIterator::new_self_intersection() }
    }

    /// Sets the data structure to process.
    pub fn set_ds(&mut self, ds: &'a BopdsDS) {
        self.inner.set_ds(ds);
    }

    /// Clears the pair lists for interference types above `level`
    /// (see [`BopdsIterator::update_by_level_of_check`]).
    pub fn update_by_level_of_check(&mut self, level: i32) {
        self.inner.update_by_level_of_check(level);
    }
}

impl<'a> std::ops::Deref for BopdsIteratorSI<'a> {
    type Target = BopdsIterator<'a>;
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<'a> std::ops::DerefMut for BopdsIteratorSI<'a> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

impl Default for BopdsIteratorSI<'_> {
    fn default() -> Self {
        Self::new()
    }
}
