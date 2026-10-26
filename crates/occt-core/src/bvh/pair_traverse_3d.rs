//! 3D `BVH_PairTraverse<double, 3>` selection walk.
//!
//! Source: `FoundationClasses/TKMath/BVH/BVH_Traverse.lxx:135-...`
//! (`BVH_PairTraverse::Select(theBVH1, theBVH2)`), with the callback set of
//! `BVH_Traverse.hxx:276-355` (`RejectNode` on a node box pair, `Accept` on a
//! leaf element pair, `Stop` / `IsMetricBetter` / `RejectMetric` from
//! `BVH_BaseTraverse`).
//!
//! The walk keeps at most `3 * BVH_Constants_MaxTreeDepth` pairs queued: each
//! step pushes at most four candidate child pairs and pops one.
#![allow(dead_code)]

use super::box_3d::BvhVec3d;
use super::tree_3d::{BvhTree3d, BVH_MAX_TREE_DEPTH};

/// `BVH_BaseTraverse` hooks plus the pair callbacks
/// (`BVH_Traverse.hxx:290-329`).
pub trait PairTraverse3 {
    /// `RejectNode(cmin1, cmax1, cmin2, cmax2, metric)` (`hxx:297-304`).
    fn reject_node(
        &self,
        c_min1: BvhVec3d,
        c_max1: BvhVec3d,
        c_min2: BvhVec3d,
        c_max2: BvhVec3d,
        metric: &mut f64,
    ) -> bool;

    /// `Accept(theIndex1, theIndex2)` (`hxx:307`).
    fn accept(&mut self, index1: i32, index2: i32) -> bool;

    /// `Stop()` (`BVH_BaseTraverse`, default `false`).
    fn stop(&self) -> bool {
        false
    }

    /// `IsMetricBetter` (`BVH_BaseTraverse`, default `<` for `double`).
    fn is_metric_better(&self, left: f64, right: f64) -> bool {
        left < right
    }

    /// `RejectMetric` (`BVH_BaseTraverse`, default `false`).
    fn reject_metric(&self, _metric: f64) -> bool {
        false
    }
}

/// `BVH_PairNodesInStack` (`BVH_Traverse.hxx:331-347`).
#[derive(Debug, Clone, Copy)]
struct PairInStack {
    node_id1: i32,
    node_id2: i32,
    metric: f64,
}

/// `BVH_PairTraverse::Select(theBVH1, theBVH2)`
/// (`BVH_Traverse.lxx:135-...`).
pub fn select_pair_tree3<S: PairTraverse3>(
    bvh1: &BvhTree3d,
    bvh2: &BvhTree3d,
    sel: &mut S,
) -> i32 {
    if bvh1.length() == 0 || bvh2.length() == 0 {
        return 0;
    }

    // `aMaxNbPairsInStack = 3 * BVH_Constants_MaxTreeDepth`.
    let max_in_stack = 3 * BVH_MAX_TREE_DEPTH;
    let mut stack = vec![
        PairInStack {
            node_id1: 0,
            node_id2: 0,
            metric: 0.0,
        };
        max_in_stack
    ];

    let mut node = PairInStack {
        node_id1: 0,
        node_id2: 0,
        metric: 0.0,
    };
    let mut prev = node;
    let mut head: i32 = -1;
    let mut nb_accepted: i32 = 0;

    loop {
        let data1 = bvh1.node_info()[node.node_id1 as usize];
        let data2 = bvh2.node_info()[node.node_id2 as usize];

        if data1.x != 0 && data2.x != 0 {
            // Outer/Outer - both nodes are leaves
            let mut metric = 0.0;
            let is_rejected = sel.reject_node(
                bvh1.min_point(node.node_id1 as usize),
                bvh1.max_point(node.node_id1 as usize),
                bvh2.min_point(node.node_id2 as usize),
                bvh2.max_point(node.node_id2 as usize),
                &mut metric,
            );

            if !is_rejected {
                let mut i_n1 = data1.y;
                while i_n1 <= data1.z {
                    let mut i_n2 = data2.y;
                    while i_n2 <= data2.z {
                        if sel.accept(i_n1, i_n2) {
                            nb_accepted += 1;
                        }
                        if sel.stop() {
                            return nb_accepted;
                        }
                        i_n2 += 1;
                    }
                    i_n1 += 1;
                }
            }
        } else {
            let mut pairs = [
                PairInStack {
                    node_id1: 0,
                    node_id2: 0,
                    metric: 0.0,
                };
                4
            ];
            let mut nb_pairs = 0;

            if data1.x == 0 && data2.x == 0 {
                // Inner/Inner
                pairs[nb_pairs] = PairInStack {
                    node_id1: data1.y,
                    node_id2: data2.y,
                    metric: 0.0,
                };
                nb_pairs += 1;
                pairs[nb_pairs] = PairInStack {
                    node_id1: data1.y,
                    node_id2: data2.z,
                    metric: 0.0,
                };
                nb_pairs += 1;
                pairs[nb_pairs] = PairInStack {
                    node_id1: data1.z,
                    node_id2: data2.y,
                    metric: 0.0,
                };
                nb_pairs += 1;
                pairs[nb_pairs] = PairInStack {
                    node_id1: data1.z,
                    node_id2: data2.z,
                    metric: 0.0,
                };
                nb_pairs += 1;
            } else if data1.x == 0 {
                // Inner/Outer
                pairs[nb_pairs] = PairInStack {
                    node_id1: data1.y,
                    node_id2: node.node_id2,
                    metric: 0.0,
                };
                nb_pairs += 1;
                pairs[nb_pairs] = PairInStack {
                    node_id1: data1.z,
                    node_id2: node.node_id2,
                    metric: 0.0,
                };
                nb_pairs += 1;
            } else if data2.x == 0 {
                // Outer/Inner
                pairs[nb_pairs] = PairInStack {
                    node_id1: node.node_id1,
                    node_id2: data2.y,
                    metric: 0.0,
                };
                nb_pairs += 1;
                pairs[nb_pairs] = PairInStack {
                    node_id1: node.node_id1,
                    node_id2: data2.z,
                    metric: 0.0,
                };
                nb_pairs += 1;
            }

            let mut kept_pairs = [
                PairInStack {
                    node_id1: 0,
                    node_id2: 0,
                    metric: 0.0,
                };
                4
            ];
            let mut nb_kept = 0;
            // Compute metrics for the nodes
            for i_pair in 0..nb_pairs {
                let is_pair_rejected = sel.reject_node(
                    bvh1.min_point(pairs[i_pair].node_id1 as usize),
                    bvh1.max_point(pairs[i_pair].node_id1 as usize),
                    bvh2.min_point(pairs[i_pair].node_id2 as usize),
                    bvh2.max_point(pairs[i_pair].node_id2 as usize),
                    &mut pairs[i_pair].metric,
                );
                if !is_pair_rejected {
                    // Put the item into the sorted array of pairs
                    let mut i_sort = nb_kept;
                    while i_sort > 0
                        && sel.is_metric_better(pairs[i_pair].metric, kept_pairs[i_sort - 1].metric)
                    {
                        kept_pairs[i_sort] = kept_pairs[i_sort - 1];
                        i_sort -= 1;
                    }
                    kept_pairs[i_sort] = pairs[i_pair];
                    nb_kept += 1;
                }
            }

            if nb_kept > 0 {
                node = kept_pairs[0];

                for i_pair in 1..nb_kept {
                    head += 1;
                    stack[head as usize] = kept_pairs[i_pair];
                }
            }
        }

        if node.node_id1 == prev.node_id1 && node.node_id2 == prev.node_id2 {
            if head < 0 {
                return nb_accepted;
            }
            node = stack[head as usize];
            head -= 1;
            while sel.reject_metric(node.metric) {
                if head < 0 {
                    return nb_accepted;
                }
                node = stack[head as usize];
                head -= 1;
            }
        }

        prev = node;
    }
}
