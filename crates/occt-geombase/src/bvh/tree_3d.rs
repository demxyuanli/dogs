//! 3D `BVH_Tree<double, 3, BVH_BinaryTree>`.
//!
//! Source: `FoundationClasses/TKMath/BVH/BVH_Tree.hxx` (`BVH_TreeBase`) and
//! `BVH_BinaryTree.hxx`:
//! * buffers `myNodeInfoBuffer` / `myMinPointBuffer` / `myMaxPointBuffer`
//!   (`Tree.hxx:204-211`), `myDepth` (`Tree.hxx:213`)
//! * `IsOuter` (`Tree.hxx:148-151`), `BegPrimitive` / `EndPrimitive` /
//!   `NbPrimitives` (`Tree.hxx:107-133`), `Level` (`Tree.hxx:137-146`)
//! * `Clear` (`BinaryTree.hxx:78-84`)
//! * `AddLeafNode(Beg, End)` with uninitialized bounds
//!   (`BinaryTree.hxx:131-135`)
//! * `AddInnerNode(Lft, Rgh)` with uninitialized bounds
//!   (`BinaryTree.hxx:138-142`)
//! * `SetInner` / `SetOuter` (`BinaryTree.hxx:40-49`)
//!
//! Node `0` is the root, children are appended in creation order, and a leaf
//! stores the inclusive element range `y..=z` - the exact layout
//! `BVH_Traverse.lxx:135-...` walks.
#![allow(dead_code)]

use super::box_3d::BvhVec3d;

/// `BVH_Vec4i` (`BVH_Types.hxx:96-120`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BvhVec4i {
    pub x: i32,
    pub y: i32,
    pub z: i32,
    pub w: i32,
}

impl BvhVec4i {
    pub fn new(x: i32, y: i32, z: i32, w: i32) -> Self {
        Self { x, y, z, w }
    }

    pub fn coord(&self, i: usize) -> i32 {
        match i {
            0 => self.x,
            1 => self.y,
            2 => self.z,
            _ => self.w,
        }
    }
}

/// `BVH_Constants_MaxTreeDepth` (`BVH_Constants.hxx:33`).
pub const BVH_MAX_TREE_DEPTH: usize = 32;

/// `BVH_Constants_LeafNodeSizeDefault` (`BVH_Constants.hxx:26`).
pub const BVH_LEAF_NODE_SIZE_DEFAULT: usize = 4;

/// `BVH_Tree<double, 3, BVH_BinaryTree>`.
#[derive(Debug, Clone, Default)]
pub struct BvhTree3d {
    node_info: Vec<BvhVec4i>,
    min: Vec<BvhVec3d>,
    max: Vec<BvhVec3d>,
    depth: usize,
}

impl BvhTree3d {
    /// `BVH_Tree()` (`BinaryTree.hxx:32-36`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Depth()` (`Tree.hxx:75`).
    pub fn depth(&self) -> usize {
        self.depth
    }

    /// `updateDepth` (`BVH_Builder.hxx:76-83`).
    pub fn update_depth(&mut self, level: usize) {
        if level > self.depth {
            self.depth = level;
        }
    }

    /// `Length()` (`Tree.hxx:78`).
    pub fn length(&self) -> usize {
        self.node_info.len()
    }

    /// `Clear()` (`BinaryTree.hxx:78-84`).
    pub fn clear(&mut self) {
        self.depth = 0;
        self.min.clear();
        self.max.clear();
        self.node_info.clear();
    }

    /// `NodeInfoBuffer()` (`Tree.hxx:156-161`).
    pub fn node_info(&self) -> &[BvhVec4i] {
        &self.node_info
    }

    /// `MinPointBuffer()` (`Tree.hxx:163-168`).
    pub fn min_buffer(&self) -> &[BvhVec3d] {
        &self.min
    }

    /// `MaxPointBuffer()` (`Tree.hxx:163-168`).
    pub fn max_buffer(&self) -> &[BvhVec3d] {
        &self.max
    }

    /// `MinPointBuffer().Resize(Length())` / `MaxPointBuffer().Resize(...)`
    /// (`BVH_LinearBuilder.hxx:324-325`).
    pub fn resize_bounds(&mut self, size: usize) {
        self.min.resize(size, BvhVec3d::default());
        self.max.resize(size, BvhVec3d::default());
    }

    /// `MinPoint(theNodeIndex)` (`Tree.hxx:84-98`).
    pub fn min_point(&self, node: usize) -> BvhVec3d {
        self.min[node]
    }

    /// `MaxPoint(theNodeIndex)` (`Tree.hxx:89-103`).
    pub fn max_point(&self, node: usize) -> BvhVec3d {
        self.max[node]
    }

    /// `SetMinPoint` / `SetMaxPoint` on the mutable buffers.
    pub fn set_bounds(&mut self, node: usize, min: BvhVec3d, max: BvhVec3d) {
        self.min[node] = min;
        self.max[node] = max;
    }

    /// `IsOuter(theNodeIndex)` (`Tree.hxx:148-151`).
    pub fn is_outer(&self, node: usize) -> bool {
        self.node_info[node].x != 0
    }

    /// `SetOuter` (`BinaryTree.hxx:40-43`).
    pub fn set_outer(&mut self, node: usize) {
        self.node_info[node].x = 1;
    }

    /// `SetInner` (`BinaryTree.hxx:46-49`).
    pub fn set_inner(&mut self, node: usize) {
        self.node_info[node].x = 0;
    }

    /// `BegPrimitive(theNodeIndex)` (`Tree.hxx:116-120`).
    pub fn beg_primitive(&self, node: usize) -> i32 {
        self.node_info[node].y
    }

    /// `EndPrimitive(theNodeIndex)` (`Tree.hxx:123-127`).
    pub fn end_primitive(&self, node: usize) -> i32 {
        self.node_info[node].z
    }

    /// `NbPrimitives(theNodeIndex)` (`Tree.hxx:129-133`).
    pub fn nb_primitives(&self, node: usize) -> i32 {
        self.end_primitive(node) - self.beg_primitive(node) + 1
    }

    /// `AddLeafNode(BegElem, EndElem)` with uninitialized bounds
    /// (`BinaryTree.hxx:131-135`).
    pub fn add_leaf_node(&mut self, beg_elem: i32, end_elem: i32) -> usize {
        self.node_info.push(BvhVec4i::new(1, beg_elem, end_elem, 0));
        self.node_info.len() - 1
    }

    /// `AddInnerNode(LftChild, RghChild)` with uninitialized bounds
    /// (`BinaryTree.hxx:138-142`).
    pub fn add_inner_node(&mut self, lft_child: i32, rgh_child: i32) -> usize {
        self.node_info.push(BvhVec4i::new(0, lft_child, rgh_child, 0));
        self.node_info.len() - 1
    }

    /// `NodeInfoBuffer()[aNode].y() = aLftChild` /
    /// `.z() = aRghChild` (`BVH_LinearBuilder.hxx:144-145`).
    pub fn set_children(&mut self, node: usize, lft_child: i32, rgh_child: i32) {
        self.node_info[node].y = lft_child;
        self.node_info[node].z = rgh_child;
    }

    /// `NodeInfoBuffer()[theNode].y()` / `.z()` - the child pair of an inner
    /// node (`BVH_Traverse.lxx:204-221`).
    pub fn children(&self, node: usize) -> (i32, i32) {
        (self.node_info[node].y, self.node_info[node].z)
    }
}
