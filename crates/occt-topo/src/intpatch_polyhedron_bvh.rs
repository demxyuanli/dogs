//! `IntPatch_PolyhedronBVH` and `IntPatch_BVHTraversal`.
//!
//! Sources (OCCT 8.0.0):
//! * `IntPatch_PolyhedronBVH.hxx/.cxx` - the polyhedron wrapped as a
//!   `BVH_PrimitiveSet<double, 3>`.
//! * `IntPatch_BVHTraversal.hxx/.cxx` - the dual-tree walk that collects the
//!   triangle pairs to feed to `IntPatch_InterferencePolyhedron::Intersect`.
//!
//! The generic machinery below them is ported in `occt_core::bvh`:
//! `box_3d.rs` (`BVH_Box`), `tree_3d.rs` (`BVH_Tree`), `radix_sorter_3d.rs`
//! (`BVH_RadixSorter`), `linear_builder_3d.rs` (`BVH_LinearBuilder`) and
//! `pair_traverse_3d.rs` (`BVH_PairTraverse::Select`).
//!
//! Faithful points worth naming:
//! * `Init` only appends the triangles whose component box is not void, so
//!   degenerate triangles never reach the tree (`cxx:56-75`).
//! * `Box` returns the *component* box of the original triangle, i.e. the one
//!   already enlarged by the deflection (`cxx:103-127`).
//! * `Center` is the centroid of the three vertices, not the box centre
//!   (`cxx:132-163`).
//! * `Accept` records original 1-based triangle indices and, in
//!   self-interference mode, drops pairs with `index1 >= index2`
//!   (`cxx:94-113`).
#![allow(dead_code)]

use occt_core::bvh::box_3d::{BvhBox3d, BvhVec3d};
use occt_core::bvh::linear_builder_3d::BvhLinearBuilder3;
use occt_core::bvh::pair_traverse_3d::{select_pair_tree3, PairTraverse3};
use occt_core::bvh::radix_sorter_3d::BvhSet3;
use occt_core::bvh::tree_3d::BvhTree3d;

use crate::intpatch_polyhedron::IntPatchPolyhedron;

/// `IntPatch_PolyhedronBVH` (`hxx:29-82`).
///
/// OCCT keeps a `const IntPatch_Polyhedron*`; the port carries the borrow so
/// the set cannot outlive the polyhedron it indexes.
#[derive(Clone)]
pub struct IntPatchPolyhedronBVH<'a> {
    /// `myPoly` (`hxx:81`).
    poly: Option<&'a IntPatchPolyhedron>,
    /// `myIndexMap` (`hxx:82`): 0-based slot -> original 1-based triangle.
    index_map: Vec<i32>,
    /// `myBVH` (`BVH_PrimitiveSet.hxx:108`).
    bvh: BvhTree3d,
    /// `myIsDirty` (`BVH_Object.hxx`).
    dirty: bool,
    /// `myBox` (`BVH_PrimitiveSet.hxx:110`), cached `BVH_Set::Box()`.
    box_cache: BvhBox3d,
}

impl Default for IntPatchPolyhedronBVH<'_> {
    /// `IntPatch_PolyhedronBVH()` (`cxx:24-31`).
    fn default() -> Self {
        Self {
            poly: None,
            index_map: Vec::new(),
            bvh: BvhTree3d::new(),
            dirty: true,
            box_cache: BvhBox3d::new(),
        }
    }
}

impl<'a> IntPatchPolyhedronBVH<'a> {
    /// `IntPatch_PolyhedronBVH()` (`cxx:24-31`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `IntPatch_PolyhedronBVH(thePoly)` (`cxx:36-42`).
    pub fn with_poly(poly: &'a IntPatchPolyhedron) -> Self {
        let mut me = Self::default();
        me.init(poly);
        me
    }

    /// `Init(thePoly)` (`cxx:55-76`).
    pub fn init(&mut self, poly: &'a IntPatchPolyhedron) {
        self.poly = Some(poly);

        let nb_triangles = poly.nb_triangles();
        let comp_bnd = poly.components_bounding();

        self.index_map.clear();
        for idx in 1..=nb_triangles {
            // `aCompBnd->Value(anIdx).IsVoid()`
            if comp_bnd[idx].get().is_some() {
                self.index_map.push(idx as i32);
            }
        }

        // `MarkDirty()` (`BVH_Object.hxx`).
        self.dirty = true;
    }

    /// `Clear()` (`cxx:81-86`).
    pub fn clear(&mut self) {
        self.poly = None;
        self.index_map.clear();
        self.dirty = true;
    }

    /// `IsInitialized()` (`hxx:78`).
    pub fn is_initialized(&self) -> bool {
        self.poly.is_some()
    }

    /// `Size()` (`cxx:91-98`).
    pub fn size(&self) -> usize {
        if self.poly.is_none() {
            return 0;
        }
        self.index_map.len()
    }

    /// `OriginalIndex(theIndex)` (`cxx:176-185`): 1-based triangle, or `0`
    /// when out of range.
    pub fn original_index(&self, index: i32) -> i32 {
        if index < 0 || index >= self.size() as i32 {
            return 0;
        }
        self.index_map[index as usize]
    }

    /// `Box(theIndex)` (`cxx:103-127`): the pre-computed component box of the
    /// original triangle, which already includes the deflection.
    pub fn box_at(&self, index: i32) -> BvhBox3d {
        let Some(poly) = self.poly else {
            return BvhBox3d::new();
        };
        if index < 0 || index >= self.size() as i32 {
            return BvhBox3d::new();
        }
        let orig_idx = self.index_map[index as usize];
        let bnd = &poly.components_bounding()[orig_idx as usize];
        BvhBox3d::from_bnd(bnd)
    }

    /// `Center(theIndex, theAxis)` (`cxx:132-163`).
    pub fn center(&self, index: i32, axis: i32) -> f64 {
        let Some(poly) = self.poly else {
            return 0.0;
        };
        if index < 0 || index >= self.size() as i32 {
            return 0.0;
        }
        let orig_idx = self.index_map[index as usize];
        let (p1, p2, p3) = poly.triangle(orig_idx as usize);
        let a = poly.point(p1);
        let b = poly.point(p2);
        let c = poly.point(p3);
        match axis {
            0 => (a.x() + b.x() + c.x()) / 3.0,
            1 => (a.y() + b.y() + c.y()) / 3.0,
            2 => (a.z() + b.z() + c.z()) / 3.0,
            _ => 0.0,
        }
    }

    /// `Swap(theIndex1, theIndex2)` (`cxx:168-173`).
    pub fn swap(&mut self, index1: i32, index2: i32) {
        if index1 != index2 {
            self.index_map.swap(index1 as usize, index2 as usize);
        }
    }

    /// `Box()` (`BVH_PrimitiveSet.hxx:69-76`).
    pub fn box_all(&mut self) -> BvhBox3d {
        if self.dirty {
            self.box_cache = <Self as BvhSet3>::box_all(self);
        }
        self.box_cache
    }

    /// `Update()` (`BVH_PrimitiveSet.hxx:98-105`): build the tree when the set
    /// changed.
    pub fn update(&mut self) {
        if !self.dirty {
            return;
        }
        let box_all = self.box_all();
        let builder = BvhLinearBuilder3::default();
        let mut tree = std::mem::take(&mut self.bvh);
        builder.build(self, &mut tree, &box_all);
        self.bvh = tree;
        self.dirty = false;
    }

    /// `BVH()` (`BVH_PrimitiveSet.hxx:78-85`).
    pub fn bvh(&mut self) -> &BvhTree3d {
        if self.dirty {
            self.update();
        }
        &self.bvh
    }
}

impl BvhSet3 for IntPatchPolyhedronBVH<'_> {
    /// `Size()`.
    fn size(&self) -> usize {
        IntPatchPolyhedronBVH::size(self)
    }

    /// `Box(theIndex)`.
    fn box_at(&self, index: usize) -> BvhBox3d {
        IntPatchPolyhedronBVH::box_at(self, index as i32)
    }

    /// `Swap(theIndex1, theIndex2)`.
    fn swap_elements(&mut self, index1: usize, index2: usize) {
        IntPatchPolyhedronBVH::swap(self, index1 as i32, index2 as i32);
    }

    /// `Center(theIndex, theAxis)` (`cxx:132-163`), the vertex centroid.
    fn center_at(&self, index: usize, axis: usize) -> f64 {
        IntPatchPolyhedronBVH::center(self, index as i32, axis as i32)
    }
}

/// `IntPatch_BVHTraversal::TrianglePair` (`hxx:33-44`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrianglePair {
    /// Triangle index in the first polyhedron, 1-based.
    pub first: i32,
    /// Triangle index in the second polyhedron, 1-based.
    pub second: i32,
}

impl TrianglePair {
    pub fn new(the_first: i32, the_second: i32) -> Self {
        Self {
            first: the_first,
            second: the_second,
        }
    }
}

/// `IntPatch_BVHTraversal` (`hxx:28-94`).
#[derive(Debug, Clone, Default)]
pub struct IntPatchBVHTraversal {
    /// `mySelfInterference` (`hxx:91`).
    self_interference: bool,
    /// `myPairs` (`hxx:92`).
    pairs: Vec<TrianglePair>,
}

impl IntPatchBVHTraversal {
    /// `IntPatch_BVHTraversal()` (`cxx:20-25`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Pairs()` (`hxx:64`).
    pub fn pairs(&self) -> &[TrianglePair] {
        &self.pairs
    }

    /// `Clear()` (`hxx:67`).
    pub fn clear(&mut self) {
        self.pairs.clear();
    }

    /// `Perform(theSet1, theSet2, theSelfInterference)` (`cxx:34-61`).
    pub fn perform(
        &mut self,
        set1: &mut IntPatchPolyhedronBVH<'_>,
        set2: &mut IntPatchPolyhedronBVH<'_>,
        self_interference: bool,
    ) -> i32 {
        self.pairs.clear();
        self.self_interference = self_interference;

        if !set1.is_initialized() || !set2.is_initialized() {
            return 0;
        }
        if set1.size() == 0 || set2.size() == 0 {
            return 0;
        }

        // `theSet1.BVH()` builds the trees; neither is null after `Update`,
        // so the null checks at `cxx:55-60` collapse into the empty checks
        // above.
        set1.update();
        set2.update();

        let tree1 = std::mem::take(&mut set1.bvh);
        let tree2 = std::mem::take(&mut set2.bvh);

        let mut selector = TraversalSelector {
            map1: &set1.index_map,
            map2: &set2.index_map,
            self_interference,
            pairs: &mut self.pairs,
        };
        let nb = select_pair_tree3(&tree1, &tree2, &mut selector);

        set1.bvh = tree1;
        set2.bvh = tree2;

        nb
    }
}

/// `Accept` / `RejectNode` of `IntPatch_BVHTraversal` (`cxx:66-113`).
struct TraversalSelector<'a> {
    map1: &'a [i32],
    map2: &'a [i32],
    self_interference: bool,
    pairs: &'a mut Vec<TrianglePair>,
}

impl PairTraverse3 for TraversalSelector<'_> {
    /// `RejectNode` (`cxx:66-90`): AABB overlap on all three axes.
    fn reject_node(
        &self,
        c_min1: BvhVec3d,
        c_max1: BvhVec3d,
        c_min2: BvhVec3d,
        c_max2: BvhVec3d,
        _metric: &mut f64,
    ) -> bool {
        if c_min1.x > c_max2.x || c_max1.x < c_min2.x {
            return true;
        }
        if c_min1.y > c_max2.y || c_max1.y < c_min2.y {
            return true;
        }
        if c_min1.z > c_max2.z || c_max1.z < c_min2.z {
            return true;
        }
        false
    }

    /// `Accept(theIndex1, theIndex2)` (`cxx:94-113`).
    fn accept(&mut self, index1: i32, index2: i32) -> bool {
        if index1 < 0 || index1 as usize >= self.map1.len() {
            return false;
        }
        if index2 < 0 || index2 as usize >= self.map2.len() {
            return false;
        }
        let orig1 = self.map1[index1 as usize];
        let orig2 = self.map2[index2 as usize];

        if self.self_interference && orig1 >= orig2 {
            return false;
        }

        self.pairs.push(TrianglePair::new(orig1, orig2));
        true
    }
}
