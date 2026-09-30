//! `BOPTools_BoxSelector<2>` — UV-box selection from a 2D BVH.
//!
//! Source: `BOPTools_BoxSelector.hxx`. PerformAreas (`BuilderFace.cxx:501-515`)
//! does `aSelector.SetBVHSet(&aBoxTree); aSelector.SetBox(Bnd2BVH(aBox));
//! aSelector.Select();` then iterates `aSelector.Indices()`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox2d;

use crate::bvh_box2d::{BvhBox2d, BvhVec2d};
use crate::bvh_traverse2d::{select_tree, BoxTraverse2d, BvhTree2d};

/// `BOPTools_BoxSelector<2>`.
#[derive(Debug, Clone)]
pub struct Box2dSelector {
    box_: BvhBox2d,
    indices: Vec<i32>,
    tree: Option<BvhTree2d>,
    elements: Vec<i32>,
    boxes: Vec<BvhBox2d>,
}

impl Box2dSelector {
    pub fn new() -> Self {
        Self {
            box_: BvhBox2d::new(),
            indices: Vec::new(),
            tree: None,
            elements: Vec::new(),
            boxes: Vec::new(),
        }
    }

    /// `Clear`.
    pub fn clear(&mut self) {
        self.indices.clear();
    }

    /// `SetBox`.
    pub fn set_box(&mut self, box_: BvhBox2d) {
        self.box_ = box_;
    }

    pub fn set_box_from_bnd(&mut self, box_: &BndBox2d) {
        self.box_ = BvhBox2d::from_bnd(box_);
    }

    /// Bind the BVH set (tree + per-leaf element id + per-leaf BVH box).
    pub fn set_bvh_set(&mut self, tree: BvhTree2d, elements: Vec<i32>, boxes: Vec<BvhBox2d>) {
        self.tree = Some(tree);
        self.elements = elements;
        self.boxes = boxes;
    }

    /// `Indices`.
    pub fn indices(&self) -> &[i32] {
        &self.indices
    }

    /// `Select` — returns the number of accepted elements.
    pub fn select(&mut self) -> i32 {
        self.indices.clear();
        let Some(tree) = self.tree.clone() else {
            return 0;
        };
        select_tree(&tree, self)
    }

    /// `RejectElement` — `myBox.IsOut(myBVHSet->Box(theIndex))`.
    pub fn reject_element(&self, the_index: i32) -> bool {
        let i = the_index as usize;
        if i >= self.boxes.len() {
            return true;
        }
        self.box_.is_out_box(&self.boxes[i])
    }

    fn element(&self, the_index: i32) -> i32 {
        let i = the_index as usize;
        if i < self.elements.len() {
            self.elements[i]
        } else {
            the_index
        }
    }
}

impl Default for Box2dSelector {
    fn default() -> Self {
        Self::new()
    }
}

impl BoxTraverse2d for Box2dSelector {
    /// `AcceptMetric(theIsInside)` — accept the whole branch when the query
    /// box fully contains the node.
    fn accept_metric(&self, the_is_inside: bool) -> bool {
        the_is_inside
    }

    /// `RejectNode(theCMin, theCMax, theIsInside)`.
    ///
    /// `theIsInside = myBox.Contains(theCMin, theCMax, hasOverlap)`;
    /// return `!hasOverlap`.
    fn reject_node(&self, the_c_min: BvhVec2d, the_c_max: BvhVec2d, the_is_inside: &mut bool) -> bool {
        let (inside, has_overlap) = self.box_.contains_corners(the_c_min, the_c_max);
        *the_is_inside = inside;
        !has_overlap
    }

    /// `Accept(theIndex, theIsInside)`.
    fn accept(&mut self, the_index: i32, the_is_inside: bool) -> bool {
        if the_is_inside || !self.reject_element(the_index) {
            self.indices.push(self.element(the_index));
            true
        } else {
            false
        }
    }
}
