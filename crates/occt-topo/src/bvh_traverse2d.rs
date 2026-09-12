//! 2D `BVH_Traverse` / `BVH_PairTraverse` selection walk.
//!
//! Source: `BVH_Traverse.lxx` (`Select` at 21, pair `Select` at 135).
//! `BOPTools_BoxSelector<2>` is a `BVH_Traverse<double, 2, BVH_BoxSet, bool>`
//! whose metric is `bool` (`theIsInside`). PerformAreas only needs the
//! single-tree walk; the pair walk is the same algorithm used by
//! `BOPTools_PairSelector<2>`.

use crate::bvh_box2d::{BvhBox2d, BvhVec2d};

/// `BVH_Constants_MaxTreeDepth` — OCCT default 32.
pub const MAX_TREE_DEPTH: usize = 32;

/// Packed node: `x == 0` inner (`y` left, `z` right); `x != 0` leaf
/// (`y..=z` inclusive element indices). Matches `BVH_Vec4i` in `NodeInfoBuffer`.
#[derive(Debug, Clone, Copy)]
pub struct BvhNodeInfo {
    pub kind: i32,
    pub y: i32,
    pub z: i32,
}

impl BvhNodeInfo {
    pub fn inner(left: i32, right: i32) -> Self {
        Self {
            kind: 0,
            y: left,
            z: right,
        }
    }

    pub fn leaf(begin: i32, end: i32) -> Self {
        Self {
            kind: 1,
            y: begin,
            z: end,
        }
    }

    pub fn is_inner(&self) -> bool {
        self.kind == 0
    }
}

/// Binary BVH used by the traverse (min/max per node + NodeInfo).
#[derive(Debug, Clone, Default)]
pub struct BvhTree2d {
    pub nodes: Vec<BvhNodeInfo>,
    pub min: Vec<BvhVec2d>,
    pub max: Vec<BvhVec2d>,
}

impl BvhTree2d {
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    pub fn min_point(&self, i: i32) -> BvhVec2d {
        self.min[i as usize]
    }

    pub fn max_point(&self, i: i32) -> BvhVec2d {
        self.max[i as usize]
    }
}

/// `BVH_NodeInStack`.
#[derive(Debug, Clone, Copy)]
struct NodeInStack {
    node_id: i32,
    metric: bool,
}

/// `BVH_PairNodesInStack`.
#[derive(Debug, Clone, Copy)]
struct PairInStack {
    node_id1: i32,
    node_id2: i32,
    metric: f64,
}

/// Single-tree selector callbacks matching `BVH_Traverse` with `MetricType = bool`.
pub trait BoxTraverse2d {
    fn accept_metric(&self, metric: bool) -> bool {
        let _ = metric;
        false
    }

    fn reject_node(&self, cmin: BvhVec2d, cmax: BvhVec2d, metric: &mut bool) -> bool;

    fn accept(&mut self, index: i32, metric: bool) -> bool;

    fn stop(&self) -> bool {
        false
    }

    fn is_metric_better(&self, left: bool, right: bool) -> bool {
        left && !right
    }

    fn reject_metric(&self, _metric: bool) -> bool {
        false
    }
}

/// Pair-tree selector callbacks matching `BVH_PairTraverse` with `MetricType = double`.
pub trait PairTraverse2d {
    fn reject_node(
        &self,
        cmin1: BvhVec2d,
        cmax1: BvhVec2d,
        cmin2: BvhVec2d,
        cmax2: BvhVec2d,
        metric: &mut f64,
    ) -> bool;

    fn accept(&mut self, index1: i32, index2: i32) -> bool;

    fn stop(&self) -> bool {
        false
    }

    fn is_metric_better(&self, left: f64, right: f64) -> bool {
        left < right
    }

    fn reject_metric(&self, _metric: f64) -> bool {
        false
    }
}

/// `BVH_Traverse::Select(theBVH)` (`BVH_Traverse.lxx:21`).
pub fn select_tree<S: BoxTraverse2d>(the_bvh: &BvhTree2d, sel: &mut S) -> i32 {
    if the_bvh.is_empty() {
        return 0;
    }
    let mut stack = [NodeInStack {
        node_id: 0,
        metric: false,
    }; MAX_TREE_DEPTH];
    let mut a_node = NodeInStack {
        node_id: 0,
        metric: false,
    };
    let mut a_prev = a_node;
    let mut a_head: i32 = -1;
    let mut a_nb_accepted: i32 = 0;
    loop {
        let a_data = the_bvh.nodes[a_node.node_id as usize];
        if a_data.is_inner() {
            if !sel.accept_metric(a_node.metric) {
                let mut a_metric_lft = false;
                let is_good_lft = !sel.reject_node(
                    the_bvh.min_point(a_data.y),
                    the_bvh.max_point(a_data.y),
                    &mut a_metric_lft,
                );
                if sel.stop() {
                    return a_nb_accepted;
                }
                let mut a_metric_rgh = false;
                let is_good_rgh = !sel.reject_node(
                    the_bvh.min_point(a_data.z),
                    the_bvh.max_point(a_data.z),
                    &mut a_metric_rgh,
                );
                if sel.stop() {
                    return a_nb_accepted;
                }
                if is_good_lft && is_good_rgh {
                    if sel.is_metric_better(a_metric_lft, a_metric_rgh) {
                        a_node = NodeInStack {
                            node_id: a_data.y,
                            metric: a_metric_lft,
                        };
                        a_head += 1;
                        stack[a_head as usize] = NodeInStack {
                            node_id: a_data.z,
                            metric: a_metric_rgh,
                        };
                    } else {
                        a_node = NodeInStack {
                            node_id: a_data.z,
                            metric: a_metric_rgh,
                        };
                        a_head += 1;
                        stack[a_head as usize] = NodeInStack {
                            node_id: a_data.y,
                            metric: a_metric_lft,
                        };
                    }
                } else if is_good_lft || is_good_rgh {
                    a_node = if is_good_lft {
                        NodeInStack {
                            node_id: a_data.y,
                            metric: a_metric_lft,
                        }
                    } else {
                        NodeInStack {
                            node_id: a_data.z,
                            metric: a_metric_rgh,
                        }
                    };
                }
            } else {
                a_node = NodeInStack {
                    node_id: a_data.y,
                    metric: a_node.metric,
                };
                a_head += 1;
                stack[a_head as usize] = NodeInStack {
                    node_id: a_data.z,
                    metric: a_node.metric,
                };
            }
        } else {
            let mut i_n = a_data.y;
            while i_n <= a_data.z {
                if sel.accept(i_n, a_node.metric) {
                    a_nb_accepted += 1;
                }
                if sel.stop() {
                    return a_nb_accepted;
                }
                i_n += 1;
            }
        }
        if a_node.node_id == a_prev.node_id {
            if a_head < 0 {
                return a_nb_accepted;
            }
            a_node = stack[a_head as usize];
            a_head -= 1;
            while sel.reject_metric(a_node.metric) {
                if a_head < 0 {
                    return a_nb_accepted;
                }
                a_node = stack[a_head as usize];
                a_head -= 1;
            }
        }
        a_prev = a_node;
    }
}

/// `BVH_PairTraverse::Select` (`BVH_Traverse.lxx:135`).
pub fn select_pair_tree<S: PairTraverse2d>(
    the_bvh1: &BvhTree2d,
    the_bvh2: &BvhTree2d,
    sel: &mut S,
) -> i32 {
    if the_bvh1.is_empty() || the_bvh2.is_empty() {
        return 0;
    }
    let a_max = 3 * MAX_TREE_DEPTH;
    let mut stack = vec![
        PairInStack {
            node_id1: 0,
            node_id2: 0,
            metric: 0.0,
        };
        a_max
    ];
    let mut a_node = PairInStack {
        node_id1: 0,
        node_id2: 0,
        metric: 0.0,
    };
    let mut a_prev = a_node;
    let mut a_head: i32 = -1;
    let mut a_nb_accepted: i32 = 0;
    loop {
        let a_data1 = the_bvh1.nodes[a_node.node_id1 as usize];
        let a_data2 = the_bvh2.nodes[a_node.node_id2 as usize];
        if !a_data1.is_inner() && !a_data2.is_inner() {
            let mut a_metric = 0.0;
            let is_rejected = sel.reject_node(
                the_bvh1.min_point(a_node.node_id1),
                the_bvh1.max_point(a_node.node_id1),
                the_bvh2.min_point(a_node.node_id2),
                the_bvh2.max_point(a_node.node_id2),
                &mut a_metric,
            );
            if !is_rejected {
                let mut i_n1 = a_data1.y;
                while i_n1 <= a_data1.z {
                    let mut i_n2 = a_data2.y;
                    while i_n2 <= a_data2.z {
                        if sel.accept(i_n1, i_n2) {
                            a_nb_accepted += 1;
                        }
                        if sel.stop() {
                            return a_nb_accepted;
                        }
                        i_n2 += 1;
                    }
                    i_n1 += 1;
                }
            }
        } else {
            let mut a_pairs = [PairInStack {
                node_id1: 0,
                node_id2: 0,
                metric: 0.0,
            }; 4];
            let mut a_nb_pairs = 0;
            if a_data1.is_inner() && a_data2.is_inner() {
                a_pairs[a_nb_pairs] = PairInStack {
                    node_id1: a_data1.y,
                    node_id2: a_data2.y,
                    metric: 0.0,
                };
                a_nb_pairs += 1;
                a_pairs[a_nb_pairs] = PairInStack {
                    node_id1: a_data1.y,
                    node_id2: a_data2.z,
                    metric: 0.0,
                };
                a_nb_pairs += 1;
                a_pairs[a_nb_pairs] = PairInStack {
                    node_id1: a_data1.z,
                    node_id2: a_data2.y,
                    metric: 0.0,
                };
                a_nb_pairs += 1;
                a_pairs[a_nb_pairs] = PairInStack {
                    node_id1: a_data1.z,
                    node_id2: a_data2.z,
                    metric: 0.0,
                };
                a_nb_pairs += 1;
            } else if a_data1.is_inner() {
                a_pairs[a_nb_pairs] = PairInStack {
                    node_id1: a_data1.y,
                    node_id2: a_node.node_id2,
                    metric: 0.0,
                };
                a_nb_pairs += 1;
                a_pairs[a_nb_pairs] = PairInStack {
                    node_id1: a_data1.z,
                    node_id2: a_node.node_id2,
                    metric: 0.0,
                };
                a_nb_pairs += 1;
            } else if a_data2.is_inner() {
                a_pairs[a_nb_pairs] = PairInStack {
                    node_id1: a_node.node_id1,
                    node_id2: a_data2.y,
                    metric: 0.0,
                };
                a_nb_pairs += 1;
                a_pairs[a_nb_pairs] = PairInStack {
                    node_id1: a_node.node_id1,
                    node_id2: a_data2.z,
                    metric: 0.0,
                };
                a_nb_pairs += 1;
            }
            let mut a_kept = [PairInStack {
                node_id1: 0,
                node_id2: 0,
                metric: 0.0,
            }; 4];
            let mut a_nb_kept = 0;
            for i_pair in 0..a_nb_pairs {
                let is_pair_rejected = sel.reject_node(
                    the_bvh1.min_point(a_pairs[i_pair].node_id1),
                    the_bvh1.max_point(a_pairs[i_pair].node_id1),
                    the_bvh2.min_point(a_pairs[i_pair].node_id2),
                    the_bvh2.max_point(a_pairs[i_pair].node_id2),
                    &mut a_pairs[i_pair].metric,
                );
                if !is_pair_rejected {
                    let mut i_sort = a_nb_kept;
                    while i_sort > 0
                        && sel.is_metric_better(a_pairs[i_pair].metric, a_kept[i_sort - 1].metric)
                    {
                        a_kept[i_sort] = a_kept[i_sort - 1];
                        i_sort -= 1;
                    }
                    a_kept[i_sort] = a_pairs[i_pair];
                    a_nb_kept += 1;
                }
            }
            if a_nb_kept > 0 {
                a_node = a_kept[0];
                for i_pair in 1..a_nb_kept {
                    a_head += 1;
                    stack[a_head as usize] = a_kept[i_pair];
                }
            }
        }
        if a_node.node_id1 == a_prev.node_id1 && a_node.node_id2 == a_prev.node_id2 {
            if a_head < 0 {
                return a_nb_accepted;
            }
            a_node = stack[a_head as usize];
            a_head -= 1;
            while sel.reject_metric(a_node.metric) {
                if a_head < 0 {
                    return a_nb_accepted;
                }
                a_node = stack[a_head as usize];
                a_head -= 1;
            }
        }
        a_prev = a_node;
    }
}

/// Convert a `BvhBox2d` to the (min, max) pair `RejectNode` expects.
pub fn box_corners(b: &BvhBox2d) -> (BvhVec2d, BvhVec2d) {
    (b.corner_min(), b.corner_max())
}
