//! `BOPTools_PairSelector<2>` — pairs of overlapping 2D BVH leaves.
//!
//! Source: `BOPTools_PairSelector.hxx`. Not used by PerformAreas (that uses
//! `BoxSelector`); kept because `BOPTools_BoxTree.hxx` typedefs
//! `BOPTools_Box2dPairSelector = BOPTools_PairSelector<2>`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use crate::bvh_box2d::{BvhBox2d, BvhVec2d};
use crate::bvh_traverse2d::{select_pair_tree, PairTraverse2d, BvhTree2d};

/// `BOPTools_PairSelector::PairIDs`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PairIds {
    pub id1: i32,
    pub id2: i32,
}

impl PairIds {
    pub fn new(id1: i32, id2: i32) -> Self {
        Self { id1, id2 }
    }
}

impl PartialOrd for PairIds {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PairIds {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.id1.cmp(&other.id1) {
            std::cmp::Ordering::Equal => self.id2.cmp(&other.id2),
            o => o,
        }
    }
}

/// `BOPTools_PairSelector<2>`.
#[derive(Debug, Clone)]
pub struct Box2dPairSelector {
    pairs: Vec<PairIds>,
    same_bvhs: bool,
    tree1: Option<BvhTree2d>,
    tree2: Option<BvhTree2d>,
    elements1: Vec<i32>,
    elements2: Vec<i32>,
    boxes1: Vec<BvhBox2d>,
    boxes2: Vec<BvhBox2d>,
}

impl Box2dPairSelector {
    pub fn new() -> Self {
        Self {
            pairs: Vec::new(),
            same_bvhs: false,
            tree1: None,
            tree2: None,
            elements1: Vec::new(),
            elements2: Vec::new(),
            boxes1: Vec::new(),
            boxes2: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.pairs.clear();
    }

    /// `Sort`.
    pub fn sort(&mut self) {
        self.pairs.sort();
    }

    /// `SetSame`.
    pub fn set_same(&mut self, the_is_same: bool) {
        self.same_bvhs = the_is_same;
    }

    pub fn pairs(&self) -> &[PairIds] {
        &self.pairs
    }

    pub fn set_bvh_sets(
        &mut self,
        tree1: BvhTree2d,
        elements1: Vec<i32>,
        boxes1: Vec<BvhBox2d>,
        tree2: BvhTree2d,
        elements2: Vec<i32>,
        boxes2: Vec<BvhBox2d>,
    ) {
        self.tree1 = Some(tree1);
        self.elements1 = elements1;
        self.boxes1 = boxes1;
        self.tree2 = Some(tree2);
        self.elements2 = elements2;
        self.boxes2 = boxes2;
    }

    /// `RejectElement(theID1, theID2)`.
    pub fn reject_element(&self, the_id1: i32, the_id2: i32) -> bool {
        if self.same_bvhs && the_id1 >= the_id2 {
            return true;
        }
        let i1 = the_id1 as usize;
        let i2 = the_id2 as usize;
        if i1 >= self.boxes1.len() || i2 >= self.boxes2.len() {
            return true;
        }
        self.boxes1[i1].is_out_box(&self.boxes2[i2])
    }

    pub fn select(&mut self) -> i32 {
        self.pairs.clear();
        let (Some(t1), Some(t2)) = (self.tree1.clone(), self.tree2.clone()) else {
            return 0;
        };
        select_pair_tree(&t1, &t2, self)
    }

    fn element1(&self, i: i32) -> i32 {
        let u = i as usize;
        if u < self.elements1.len() {
            self.elements1[u]
        } else {
            i
        }
    }

    fn element2(&self, i: i32) -> i32 {
        let u = i as usize;
        if u < self.elements2.len() {
            self.elements2[u]
        } else {
            i
        }
    }
}

impl Default for Box2dPairSelector {
    fn default() -> Self {
        Self::new()
    }
}

impl PairTraverse2d for Box2dPairSelector {
    /// `RejectNode` — `BVH_Box(theCMin1, theCMax1).IsOut(theCMin2, theCMax2)`.
    ///
    /// The OCCT header constructs a 3D box from 2D corners for Dimension=2 as
    /// well (`BOPTools_PairSelector.hxx:84` uses `BVH_Box<double, 3>`). The
    /// extra Z axis is unused; 2D `IsOut` is the equivalent test.
    fn reject_node(
        &self,
        the_c_min1: BvhVec2d,
        the_c_max1: BvhVec2d,
        the_c_min2: BvhVec2d,
        the_c_max2: BvhVec2d,
        _metric: &mut f64,
    ) -> bool {
        let a = BvhBox2d::from_corners(the_c_min1, the_c_max1);
        a.is_out_corners(the_c_min2, the_c_max2)
    }

    fn accept(&mut self, the_id1: i32, the_id2: i32) -> bool {
        if !self.reject_element(the_id1, the_id2) {
            self.pairs
                .push(PairIds::new(self.element1(the_id1), self.element2(the_id2)));
            true
        } else {
            false
        }
    }
}
