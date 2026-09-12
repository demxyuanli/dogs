//! 2D AABB tree standing in for `BOPTools_Box2dTree` / `BOPTools_Box2dTreeSelector`.
//!
//! Source: `BOPTools_BoxTree.hxx` (`BOPTools_BoxSet<double, 2, int>` with a
//! linear BVH builder) and the selector used by
//! `BOPAlgo_BuilderFace::PerformAreas` (`BOPAlgo_BuilderFace.cxx:470-509`):
//!
//! ```text
//! aBoxTree.SetSize(aNbH);
//! aBoxTree.Add(i, Bnd_Tools::Bnd2BVH(aBox));
//! aBoxTree.Build();
//! aSelector.SetBVHSet(&aBoxTree);
//! aSelector.SetBox(Bnd_Tools::Bnd2BVH(aBox));
//! aSelector.Select();
//! const NCollection_List<int>& aLI = aSelector.Indices();
//! ```
//!
//! A leaf is a hit when `!query.is_out_box(&leaf)` — the same predicate
//! `Bnd_Box2d::IsOut` uses. The 3D sibling is [`crate::bop_aabb_faces`].

use occt_core::bnd::BndBox2d;

/// One leaf: original 1-based OCCT index plus the stored 2D box.
#[derive(Debug, Clone)]
pub struct Box2dLeaf {
    pub index: i32,
    pub box_: BndBox2d,
}

#[derive(Debug, Clone)]
struct Box2dNode {
    box_: BndBox2d,
    left: isize,
    right: isize,
}

/// Binary AABB tree over 2D boxes (`BOPTools_Box2dTree`).
#[derive(Debug, Clone, Default)]
pub struct Box2dTree {
    leaves: Vec<Box2dLeaf>,
    nodes: Vec<Box2dNode>,
    root: isize,
}

impl Box2dTree {
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

    /// `Add(index, box)`. `index` is the OCCT 1-based map index.
    pub fn add(&mut self, index: i32, box_: BndBox2d) {
        self.leaves.push(Box2dLeaf { index, box_ });
    }

    /// `Build` — recursive median split on the longest UV axis.
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
            self.nodes.push(Box2dNode {
                box_: self.leaves[i].box_,
                left: -(i as isize + 1),
                right: -1,
            });
            return self.nodes.len() - 1;
        }
        let union = self.union_range(order, begin, end);
        let axis = longest_axis_2d(&union);
        order[begin..end].sort_by(|a, b| {
            let ca = centroid_axis_2d(&self.leaves[*a].box_, axis);
            let cb = centroid_axis_2d(&self.leaves[*b].box_, axis);
            ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
        });
        let mid = begin + count / 2;
        let left = self.build_range(order, begin, mid);
        let right = self.build_range(order, mid, end);
        let mut box_ = self.nodes[left].box_;
        box_.add_box(&self.nodes[right].box_);
        self.nodes.push(Box2dNode {
            box_,
            left: left as isize,
            right: right as isize,
        });
        self.nodes.len() - 1
    }

    fn union_range(&self, order: &[usize], begin: usize, end: usize) -> BndBox2d {
        let mut u = self.leaves[order[begin]].box_;
        for k in (begin + 1)..end {
            u.add_box(&self.leaves[order[k]].box_);
        }
        u
    }

    /// `Select` — 1-based indices of leaves whose boxes are not out of `query`.
    pub fn select(&self, query: &BndBox2d) -> Vec<i32> {
        let mut out = Vec::new();
        if self.root < 0 {
            return out;
        }
        self.select_node(self.root as usize, query, &mut out);
        out.sort_unstable();
        out
    }

    fn select_node(&self, node: usize, query: &BndBox2d, out: &mut Vec<i32>) {
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

fn longest_axis_2d(b: &BndBox2d) -> usize {
    match b.get() {
        Some((xmin, ymin, xmax, ymax)) => {
            let dx = (xmax - xmin).abs();
            let dy = (ymax - ymin).abs();
            if dx >= dy {
                0
            } else {
                1
            }
        }
        None => 0,
    }
}

fn centroid_axis_2d(b: &BndBox2d, axis: usize) -> f64 {
    match b.get() {
        Some((xmin, ymin, xmax, ymax)) => {
            if axis == 0 {
                0.5 * (xmin + xmax)
            } else {
                0.5 * (ymin + ymax)
            }
        }
        None => 0.0,
    }
}

/// Selector matching `BOPTools_Box2dTreeSelector`.
#[derive(Debug, Clone)]
pub struct Box2dTreeSelector<'a> {
    tree: Option<&'a Box2dTree>,
    query: BndBox2d,
}

impl<'a> Box2dTreeSelector<'a> {
    pub fn new() -> Self {
        Self {
            tree: None,
            query: BndBox2d::new(),
        }
    }

    /// `SetBVHSet`.
    pub fn set_bvh_set(&mut self, tree: &'a Box2dTree) {
        self.tree = Some(tree);
    }

    /// `SetBox`.
    pub fn set_box(&mut self, box_: BndBox2d) {
        self.query = box_;
    }

    /// `Clear` — drop the query box, keep the tree.
    pub fn clear(&mut self) {
        self.query = BndBox2d::new();
    }

    /// `Select` — returns the number of hits; indices via [`indices`].
    pub fn select(&self) -> Vec<i32> {
        match self.tree {
            Some(t) => t.select(&self.query),
            None => Vec::new(),
        }
    }
}

impl Default for Box2dTreeSelector<'_> {
    fn default() -> Self {
        Self::new()
    }
}
