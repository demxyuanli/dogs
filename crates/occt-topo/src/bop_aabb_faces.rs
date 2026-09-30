//! AABB selector standing in for `BOPTools_BoxTree` / `BOPTools_BoxTreeSelector`.
//!
//! Source: `BOPTools_BoxTree.hxx` (BVH of `Bnd_Box`) and the selector used by
//! `BOPAlgo_FillIn3DParts::Perform` (`BOPAlgo_Tools.cxx:1345-1354`):
//!
//! ```text
//! aSelector.SetBox(Bnd_Tools::Bnd2BVH(myBoxS));
//! aSelector.SetBVHSet(myBBTree);
//! aSelector.Select();
//! const NCollection_List<int>& aLIFP = aSelector.Indices();
//! ```
//!
//! The OCCT tree is a linear BVH over face bounding boxes. This port keeps
//! the same *select* contract (indices of boxes that are not out of a query
//! box) with a binary AABB tree so ClassifyFaces does not scan every face
//! against every solid when the face count grows. For a handful of boolean
//! arguments a linear scan is equivalent; the tree is the translation of the
//! BVH structure, not a heuristic filter.
//!
//! A box is a hit when `!query.is_out_box(&leaf.box)` — the same predicate
//! `Bnd_Box::IsOut` uses in the existing sequential port.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox;

use crate::bopalgo_tools_class::ShapeBox;

/// One leaf of the AABB tree: the original index into the ShapeBox vector
/// plus the stored box.
#[derive(Debug, Clone)]
pub struct AabbLeaf {
    pub index: usize,
    pub box_: BndBox,
}

/// Internal node: children cover `box_`.
#[derive(Debug, Clone)]
struct AabbNode {
    box_: BndBox,
    /// Leaf indices into `leaves`, or children of this node.
    left: isize,
    right: isize,
}

/// Binary AABB tree over a `BOPAlgo_VectorOfShapeBox`.
///
/// `BOPTools_BoxTree::SetSize` / `Add` / `Build` are `set_size` / `add` /
/// `build`. `BOPTools_BoxTreeSelector::Select` is [`AabbTree::select`].
#[derive(Debug, Clone, Default)]
pub struct AabbTree {
    leaves: Vec<AabbLeaf>,
    nodes: Vec<AabbNode>,
    root: isize,
}

impl AabbTree {
    pub fn new() -> Self {
        Self::default()
    }

    /// `SetSize` — reserve leaf slots.
    pub fn set_size(&mut self, n: usize) {
        self.leaves.clear();
        self.leaves.reserve(n);
        self.nodes.clear();
        self.root = -1;
    }

    /// `Add(index, box)`.
    pub fn add(&mut self, index: usize, box_: BndBox) {
        self.leaves.push(AabbLeaf { index, box_ });
    }

    /// Build from a ShapeBox vector (the ClassifyFaces filler).
    pub fn from_shape_boxes(v_sb: &[ShapeBox]) -> Self {
        let mut t = Self::new();
        t.set_size(v_sb.len());
        for (i, sb) in v_sb.iter().enumerate() {
            t.add(i, sb.box_);
        }
        t.build();
        t
    }

    /// `Build` — recursive median split on the longest axis.
    pub fn build(&mut self) {
        self.nodes.clear();
        let n = self.leaves.len();
        if n == 0 {
            self.root = -1;
            return;
        }
        let mut order: Vec<usize> = (0..n).collect();
        self.root = self.build_range(&mut order, 0, n) as isize;
    }

    fn build_range(&mut self, order: &mut [usize], begin: usize, end: usize) -> usize {
        let count = end - begin;
        if count == 1 {
            let i = order[begin];
            self.nodes.push(AabbNode {
                box_: self.leaves[i].box_,
                left: -(i as isize + 1),
                right: -1,
            });
            return self.nodes.len() - 1;
        }
        let union = self.union_range(order, begin, end);
        let axis = longest_axis(&union);
        order[begin..end].sort_by(|a, b| {
            let ca = centroid_axis(&self.leaves[*a].box_, axis);
            let cb = centroid_axis(&self.leaves[*b].box_, axis);
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = begin + count / 2;
        let left = self.build_range(order, begin, mid);
        let right = self.build_range(order, mid, end);
        let mut box_ = self.nodes[left].box_;
        box_ = union_boxes(&box_, &self.nodes[right].box_);
        self.nodes.push(AabbNode {
            box_,
            left: left as isize,
            right: right as isize,
        });
        self.nodes.len() - 1
    }

    fn union_range(&self, order: &[usize], begin: usize, end: usize) -> BndBox {
        let mut u = self.leaves[order[begin]].box_;
        for k in (begin + 1)..end {
            u = union_boxes(&u, &self.leaves[order[k]].box_);
        }
        u
    }

    /// `Select` — indices of leaves whose boxes are not out of `query`.
    pub fn select(&self, query: &BndBox) -> Vec<usize> {
        let mut out = Vec::new();
        if self.root < 0 {
            return out;
        }
        self.select_node(self.root as usize, query, &mut out);
        out.sort_unstable();
        out
    }

    fn select_node(&self, node: usize, query: &BndBox, out: &mut Vec<usize>) {
        let n = &self.nodes[node];
        if query.is_out_box(&n.box_) {
            return;
        }
        if n.left < 0 {
            let leaf = (-n.left - 1) as usize;
            if !query.is_out_box(&self.leaves[leaf].box_) {
                out.push(self.leaves[leaf].index);
            }
            return;
        }
        self.select_node(n.left as usize, query, out);
        if n.right >= 0 {
            self.select_node(n.right as usize, query, out);
        }
    }

    pub fn leaf_count(&self) -> usize {
        self.leaves.len()
    }

    pub fn is_empty(&self) -> bool {
        self.leaves.is_empty()
    }
}

fn longest_axis(b: &BndBox) -> usize {
    let (xmin, ymin, zmin, xmax, ymax, zmax) = box_corners(b);
    let dx = (xmax - xmin).abs();
    let dy = (ymax - ymin).abs();
    let dz = (zmax - zmin).abs();
    if dx >= dy && dx >= dz {
        0
    } else if dy >= dz {
        1
    } else {
        2
    }
}

fn centroid_axis(b: &BndBox, axis: usize) -> f64 {
    let (xmin, ymin, zmin, xmax, ymax, zmax) = box_corners(b);
    match axis {
        0 => 0.5 * (xmin + xmax),
        1 => 0.5 * (ymin + ymax),
        _ => 0.5 * (zmin + zmax),
    }
}

fn box_corners(b: &BndBox) -> (f64, f64, f64, f64, f64, f64) {
    match b.get() {
        Some((xmin, xmax, ymin, ymax, zmin, zmax)) => (xmin, ymin, zmin, xmax, ymax, zmax),
        None => (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    }
}

fn union_boxes(a: &BndBox, b: &BndBox) -> BndBox {
    let mut u = *a;
    u.add_box(b);
    u
}

/// Selector object matching `BOPTools_BoxTreeSelector`.
#[derive(Debug, Clone)]
pub struct AabbSelector<'a> {
    tree: &'a AabbTree,
    query: BndBox,
}

impl<'a> AabbSelector<'a> {
    pub fn new(tree: &'a AabbTree) -> Self {
        Self {
            tree,
            query: BndBox::new(),
        }
    }

    pub fn set_box(&mut self, box_: BndBox) {
        self.query = box_;
    }

    pub fn select(&self) -> Vec<usize> {
        self.tree.select(&self.query)
    }
}
