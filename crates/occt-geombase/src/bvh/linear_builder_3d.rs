//! 3D `BVH_LinearBuilder<double, 3>` - LBVH built from Morton codes.
//!
//! Source: `FoundationClasses/TKMath/BVH/BVH_LinearBuilder.hxx`:
//! * `lowerBound` (`hxx:88-111`)
//! * `emitHierachy` (`hxx:116-153`)
//! * `BVH::UpdateBounds` (`hxx:156-203`) - the recursive form; OCCT runs it
//!   through `UpdateBoundTask` (`hxx:216-291`), which only splits the first two
//!   levels across threads and then recurses into the same function, so the
//!   sequential port produces the same node bounds.
//! * `Build` (`hxx:298-331`)
//!
//! Everything the builder touches goes through [`BvhSet3`], [`BvhTree3d`] and
//! [`BvhRadixSorter3`], so a set only has to provide `Size` / `Box` / `Swap`.
#![allow(dead_code)]

use super::box_3d::{BvhBox3d, BvhVec3d};
use super::radix_sorter_3d::{BvhEncodedLink, BvhRadixSorter3, BvhSet3};
use super::tree_3d::{BvhTree3d, BVH_LEAF_NODE_SIZE_DEFAULT, BVH_MAX_TREE_DEPTH};

/// `BVH_LinearBuilder<double, 3>` (`hxx:35-69`).
#[derive(Debug, Clone)]
pub struct BvhLinearBuilder3 {
    /// `myLeafNodeSize` (`BVH_Builder.hxx:52`).
    leaf_node_size: usize,
    /// `myMaxTreeDepth` (`BVH_Builder.hxx:50`).
    max_tree_depth: usize,
}

impl Default for BvhLinearBuilder3 {
    fn default() -> Self {
        Self::new(BVH_LEAF_NODE_SIZE_DEFAULT, BVH_MAX_TREE_DEPTH)
    }
}

impl BvhLinearBuilder3 {
    /// `BVH_LinearBuilder(theLeafNodeSize, theMaxTreeDepth)` (`hxx:73-77`).
    pub fn new(leaf_node_size: usize, max_tree_depth: usize) -> Self {
        Self {
            leaf_node_size,
            max_tree_depth,
        }
    }

    /// `MaxTreeDepth()` (`BVH_Builder.hxx:29`).
    pub fn max_tree_depth(&self) -> usize {
        self.max_tree_depth
    }

    /// `lowerBound` (`hxx:88-111`): index of the first link whose `theDigit`
    /// bit is set.
    fn lower_bound(&self, links: &[BvhEncodedLink], start: i32, final_: i32, digit: i32) -> i32 {
        let mut nb_prims = final_ - start;
        let bit = 1u32 << digit;
        let mut start = start;
        while nb_prims > 0 {
            let step = nb_prims / 2;
            if links[(start + step) as usize].code & bit != 0 {
                nb_prims = step;
            } else {
                start += step + 1;
                nb_prims -= step + 1;
            }
        }
        start
    }

    /// `emitHierachy` (`hxx:116-153`).
    fn emit_hierarchy(
        &self,
        tree: &mut BvhTree3d,
        links: &[BvhEncodedLink],
        digit: i32,
        shift: i32,
        start: i32,
        final_: i32,
    ) -> i32 {
        if final_ - start > self.leaf_node_size as i32 {
            let position = if digit < 0 {
                (start + final_) / 2
            } else {
                self.lower_bound(links, start, final_, digit)
            };
            if position == start || position == final_ {
                return self.emit_hierarchy(tree, links, digit - 1, shift, start, final_);
            }

            // Build inner node
            let node = tree.add_inner_node(0, 0) as i32;
            let rgh_node = shift + position - start;

            let lft_child = self.emit_hierarchy(tree, links, digit - 1, shift, start, position);
            let rgh_child = self.emit_hierarchy(tree, links, digit - 1, rgh_node, position, final_);

            tree.set_children(node as usize, lft_child, rgh_child);
            node
        } else {
            // Build leaf node
            tree.add_leaf_node(shift, shift + final_ - start - 1) as i32
        }
    }

    /// `BVH::UpdateBounds` (`hxx:156-203`), returning the node height.
    fn update_bounds<S: BvhSet3>(&self, set: &S, tree: &mut BvhTree3d, node: usize) -> i32 {
        let data = tree.node_info()[node];
        if data.x == 0 {
            let (lft_child, rgh_child) = tree.children(node);

            let lft_depth = self.update_bounds(set, tree, lft_child as usize);
            let rgh_depth = self.update_bounds(set, tree, rgh_child as usize);

            let lft_min = tree.min_point(lft_child as usize);
            let lft_max = tree.max_point(lft_child as usize);
            let rgh_min = tree.min_point(rgh_child as usize);
            let rgh_max = tree.max_point(rgh_child as usize);

            let new_min = BvhVec3d::new(
                lft_min.x.min(rgh_min.x),
                lft_min.y.min(rgh_min.y),
                lft_min.z.min(rgh_min.z),
            );
            let new_max = BvhVec3d::new(
                lft_max.x.max(rgh_max.x),
                lft_max.y.max(rgh_max.y),
                lft_max.z.max(rgh_max.z),
            );
            tree.set_bounds(node, new_min, new_max);
            lft_depth.max(rgh_depth) + 1
        } else {
            let mut min_point = BvhVec3d::default();
            let mut max_point = BvhVec3d::default();
            for prim in data.y..=data.z {
                let a_box = set.box_at(prim as usize);
                let a_min = a_box.corner_min();
                let a_max = a_box.corner_max();
                if prim == data.y {
                    min_point = a_min;
                    max_point = a_max;
                } else {
                    min_point = BvhVec3d::new(
                        min_point.x.min(a_min.x),
                        min_point.y.min(a_min.y),
                        min_point.z.min(a_min.z),
                    );
                    max_point = BvhVec3d::new(
                        max_point.x.max(a_max.x),
                        max_point.y.max(a_max.y),
                        max_point.z.max(a_max.z),
                    );
                }
            }
            tree.set_bounds(node, min_point, max_point);
            0
        }
    }

    /// `Build(theSet, theBVH, theBox)` (`hxx:298-331`).
    pub fn build<S: BvhSet3>(&self, set: &mut S, tree: &mut BvhTree3d, the_box: &BvhBox3d) {
        let set_size = set.size();
        if set_size == 0 {
            return;
        }

        tree.clear();

        // Step 0 / 1 -- Morton codes of the virtual grid, then the radix sort.
        let mut sorter = BvhRadixSorter3::default();
        sorter.perform_range(set, 0, set_size - 1, the_box);

        // Step 2 -- Emitting BVH hierarchy from the sorted Morton codes.
        let height = emit_and_bound(self, set, tree, sorter.encoded_links(), set_size);

        if height > tree.depth() {
            tree.update_depth(height);
        }
    }
}

/// Steps 2 and 3 of `Build` (`hxx:319-329`): emit the hierarchy from the
/// sorted links, size the bound buffers, then fill them bottom-up.
fn emit_and_bound<S: BvhSet3>(
    builder: &BvhLinearBuilder3,
    set: &S,
    tree: &mut BvhTree3d,
    links: &[BvhEncodedLink],
    set_size: usize,
) -> usize {
    builder.emit_hierarchy(tree, links, 29, 0, 0, set_size as i32);

    // `MinPointBuffer().Resize(NodeInfoBuffer().Size())` and the same for the
    // max buffer (`hxx:324-325`).
    tree.resize_bounds(tree.length());

    let height = builder.update_bounds(set, tree, 0);
    height as usize
}
