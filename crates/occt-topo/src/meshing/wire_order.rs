//! Port of `ShapeAnalysis_WireOrder` (TKShHealing, `ShapeAnalysis_WireOrder.cxx`).
//!
//! Orders the edges of a wire into one connected chain, deciding for each edge
//! whether it must be reversed, and reports a status (`Same` / `Reordered` /
//! `Reversed` / `Shifted`). Used by `BRepMesh_ShapeVisitor::addWire` in 2D
//! pcurve mode: every edge contributes its pcurve start/end `(u, v)`, and the
//! chain is grown by snapping the sequence head/tail to the closest free edge.
//! A wire whose pcurves split into several closed loops (a seam traversed
//! twice) is handled by the loop-joining tail of `Perform`.
//!
//! Only the 2D mode is needed (`addWire` calls `CheckOrder(..., isClosed=true,
//! mode3d=false)`), so 3D/both modes are not carried over.

use occt_core::gp::GpPnt2d;
use occt_core::precision::SQUARE_CONFUSION;

/// OCCT `RealSmall()` — below this squared distance the join is considered exact.
const REAL_SMALL: f64 = 1e-12;

/// Result status of the order analysis (mirrors `ShapeAnalysis_WireOrder::Status`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WireOrderStatus {
    /// 0 — order unchanged.
    Same,
    /// 1 — some edges reordered (but none reversed).
    Reordered,
    /// -1 — some edges reversed (`ShapeExtend_DONE3`).
    Reversed,
    /// 3 — edges only shifted (a forward/backward rotation).
    Shifted,
}

fn sq_dist(a: &GpPnt2d, b: &GpPnt2d) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    dx * dx + dy * dy
}

/// Ordered chain of wire edges in 2D space.
pub struct WireOrder {
    begins: Vec<GpPnt2d>,
    ends: Vec<GpPnt2d>,
    /// `ord[i]` holds the signed edge number at position `i` (1-based `i`):
    /// positive = forward, negative = reversed.
    ord: Vec<i32>,
    status: WireOrderStatus,
}

impl WireOrder {
    pub fn new() -> Self {
        Self {
            begins: Vec::new(),
            ends: Vec::new(),
            ord: Vec::new(),
            status: WireOrderStatus::Same,
        }
    }

    /// Append an edge by its (start, end) 2D points, in natural edge direction.
    pub fn add_edge(&mut self, begin: GpPnt2d, end: GpPnt2d) {
        self.begins.push(begin);
        self.ends.push(end);
    }

    /// Number of loaded edges.
    pub fn nb_edges(&self) -> usize {
        self.begins.len()
    }

    /// Status of the last `perform()`.
    pub fn status(&self) -> WireOrderStatus {
        self.status
    }

    /// Signed edge number at chain position `idx` (1-based): positive = forward,
    /// negative = reversed. Identity when `perform()` has not run yet.
    pub fn ordered(&self, idx: usize) -> i32 {
        match self.ord.get(idx.wrapping_sub(1)) {
            Some(&v) if v != 0 => v,
            _ => idx as i32,
        }
    }

    /// Traversal start (`want_start == true`) or end (`false`) point of the edge
    /// whose signed chain number is `idx` (positive = forward, negative = reversed).
    fn point_at(&self, idx: i32, want_start: bool) -> GpPnt2d {
        let e = idx.unsigned_abs() as usize - 1;
        let (b, en) = (self.begins[e], self.ends[e]);
        if idx > 0 {
            if want_start { b } else { en }
        } else if want_start {
            en
        } else {
            b
        }
    }

    /// Make the order analysis (`ShapeAnalysis_WireOrder::Perform`).
    pub fn perform(&mut self) {
        self.status = WireOrderStatus::Same;
        let n = self.begins.len();
        self.ord = vec![0; n];
        if n == 0 {
            return;
        }

        let tol2 = SQUARE_CONFUSION;
        let mut is_used = vec![false; n];

        // Chain under construction; signed edge numbers (negative = reversed).
        let mut edge_seq: Vec<i32> = vec![1];
        is_used[0] = true;
        let mut first_pnt = self.begins[0];
        let mut last_pnt = self.ends[0];
        let mut loops: Vec<Vec<i32>> = Vec::new();

        loop {
            let mut best_joint_type = 3i32;
            let mut best_min = f64::MAX;
            let mut best_edge_num = usize::MAX;
            let mut is_found = false;

            for i in 0..n {
                if is_used[i] {
                    continue;
                }
                let seq_tail_edge_head = sq_dist(&last_pnt, &self.begins[i]);
                let seq_tail_edge_tail = sq_dist(&last_pnt, &self.ends[i]);
                let seq_head_edge_tail = sq_dist(&first_pnt, &self.ends[i]);
                let seq_head_edge_head = sq_dist(&first_pnt, &self.begins[i]);

                let (tail_join_type, min_to_tail) = if seq_tail_edge_head <= seq_tail_edge_tail {
                    (0i32, seq_tail_edge_head)
                } else {
                    (2i32, seq_tail_edge_tail)
                };
                let (head_joint_type, min_to_head) = if seq_head_edge_tail <= seq_head_edge_head {
                    (1i32, seq_head_edge_tail)
                } else {
                    (3i32, seq_head_edge_head)
                };

                let (cur_joint_type, cur_min) = if (min_to_tail - min_to_head).abs() < tol2 {
                    if tail_join_type < head_joint_type {
                        (tail_join_type, min_to_tail)
                    } else {
                        (head_joint_type, min_to_head)
                    }
                } else if min_to_tail <= min_to_head {
                    (tail_join_type, min_to_tail)
                } else {
                    (head_joint_type, min_to_head)
                };

                if best_min > tol2 || cur_joint_type < best_joint_type {
                    if cur_min < best_min
                        || ((cur_min == best_min || cur_min < tol2)
                            && cur_joint_type < best_joint_type)
                    {
                        is_found = true;
                        best_min = cur_min;
                        best_joint_type = cur_joint_type;
                        best_edge_num = i;
                    }
                }
            }

            if !is_found {
                break;
            }

            let close_dist = sq_dist(&first_pnt, &last_pnt);
            if best_min <= REAL_SMALL || best_min < close_dist {
                match best_joint_type {
                    0 => {
                        edge_seq.push(best_edge_num as i32 + 1);
                        last_pnt = self.ends[best_edge_num];
                    }
                    1 => {
                        edge_seq.insert(0, best_edge_num as i32 + 1);
                        first_pnt = self.begins[best_edge_num];
                    }
                    2 => {
                        edge_seq.push(-(best_edge_num as i32 + 1));
                        last_pnt = self.begins[best_edge_num];
                    }
                    _ => {
                        edge_seq.insert(0, -(best_edge_num as i32 + 1));
                        first_pnt = self.ends[best_edge_num];
                    }
                }
            } else {
                // Better to close the current loop and start a new one.
                loops.push(std::mem::take(&mut edge_seq));
                edge_seq = vec![best_edge_num as i32 + 1];
                first_pnt = self.begins[best_edge_num];
                last_pnt = self.ends[best_edge_num];
            }
            is_used[best_edge_num] = true;
        }
        loops.push(edge_seq);

        // Connect the loops (`myKeepLoops == false` path).
        let mut main_loop = loops.remove(0);
        while !loops.is_empty() {
            let mut min_dist1 = f64::MAX;
            let mut loop_num1 = usize::MAX;
            let mut cur_loop_it1 = 0usize;
            let mut direct1 = false;
            let mut main_loop_it1 = 0usize;

            for (loop_it, cur_loop) in loops.iter().enumerate() {
                let mut cur_loop_it2 = 0usize;
                let mut main_loop_it2 = 0usize;
                let mut direct2 = false;
                let mut min_dist2 = f64::MAX;
                let cur_loop_len = cur_loop.len();

                for cur_edge_it in 0..cur_loop_len {
                    let prev_edge_it = if cur_edge_it == 0 { cur_loop_len - 1 } else { cur_edge_it - 1 };
                    let cur_edge_idx = cur_loop[cur_edge_it];
                    let prev_edge_idx = cur_loop[prev_edge_it];
                    let cur_loop_first = self.point_at(cur_edge_idx, true);
                    let cur_loop_last = self.point_at(prev_edge_idx, false);

                    let mut min_dist3 = f64::MAX;
                    let mut main_loop_it3 = 0usize;
                    let mut direct3 = false;
                    let main_loop_len = main_loop.len();
                    for cur_edge_it2 in 0..main_loop_len {
                        if min_dist3 == 0.0 {
                            break;
                        }
                        let next_edge_it2 =
                            if cur_edge_it2 == main_loop_len - 1 { 0 } else { cur_edge_it2 + 1 };
                        let cur_edge_idx2 = main_loop[cur_edge_it2];
                        let next_edge_idx2 = main_loop[next_edge_it2];
                        let main_loop_first = self.point_at(cur_edge_idx2, false);
                        let main_loop_last = self.point_at(next_edge_idx2, true);

                        let direct_dist = sq_dist(&cur_loop_first, &main_loop_first)
                            + sq_dist(&cur_loop_last, &main_loop_last);
                        let reverse_dist = sq_dist(&cur_loop_first, &main_loop_last)
                            + sq_dist(&cur_loop_last, &main_loop_first);
                        let join_dist = if direct_dist < tol2 || direct_dist < 2.0 * reverse_dist {
                            direct_dist
                        } else {
                            reverse_dist
                        };
                        if join_dist < min_dist3 && (min_dist3 - join_dist).abs() > tol2 {
                            min_dist3 = join_dist;
                            direct3 = direct_dist <= reverse_dist;
                            main_loop_it3 = cur_edge_it2;
                        }
                    }
                    if min_dist3 < min_dist2 && (min_dist2 - min_dist3).abs() > tol2 {
                        min_dist2 = min_dist3;
                        direct2 = direct3;
                        main_loop_it2 = main_loop_it3;
                        cur_loop_it2 = cur_edge_it;
                    }
                }
                if min_dist2 < min_dist1 && (min_dist1 - min_dist2).abs() > tol2 {
                    min_dist1 = min_dist2;
                    loop_num1 = loop_it;
                    direct1 = direct2;
                    main_loop_it1 = main_loop_it2;
                    cur_loop_it1 = cur_loop_it2;
                }
            }

            let loop_to_insert = loops[loop_num1].clone();
            let factor = if direct1 { 1 } else { -1 };
            let insert_len = loop_to_insert.len();
            for i in 0..insert_len {
                let value_idx = (cur_loop_it1 + i) % insert_len;
                main_loop.insert(main_loop_it1 + 1 + i, loop_to_insert[value_idx] * factor);
            }
            loops.remove(loop_num1);
        }

        // Assign the result order and derive the status.
        let mut temp_status = 0i32;
        let len = main_loop.len();
        for i in 0..len {
            let val = main_loop[i];
            if (i as i32 + 1) != val && temp_status >= 0 {
                temp_status = if val > 0 { 1 } else { -1 };
            }
            self.ord[i] = val;
        }

        if temp_status == 0 {
            self.status = WireOrderStatus::Same;
            return;
        }

        let l = len as i32;
        let mut is_shift_forward = true;
        let mut is_shift_reverse = true;
        for i in 0..len.saturating_sub(1) {
            let first = main_loop[i];
            let second = main_loop[i + 1];
            if second - first != 1 && !(first == l && second == 1) {
                is_shift_forward = false;
            }
            if first - second != 1 && !(second == l && first == 1) {
                is_shift_reverse = false;
            }
        }
        let first = main_loop[len - 1];
        let second = main_loop[0];
        if second - first != 1 && !(first == l && second == 1) {
            is_shift_forward = false;
        }
        if first - second != 1 && !(second == l && first == 1) {
            is_shift_reverse = false;
        }
        if is_shift_forward || is_shift_reverse {
            temp_status = 3;
        }
        self.status = match temp_status {
            1 => WireOrderStatus::Reordered,
            -1 => WireOrderStatus::Reversed,
            3 => WireOrderStatus::Shifted,
            _ => WireOrderStatus::Same,
        };
    }
}

impl Default for WireOrder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(x: f64, y: f64) -> GpPnt2d {
        GpPnt2d::new(x, y)
    }

    /// An already-ordered unit square: identity order, no reversal.
    #[test]
    fn ordered_square_is_unchanged() {
        let mut o = WireOrder::new();
        o.add_edge(pt(0.0, 0.0), pt(1.0, 0.0));
        o.add_edge(pt(1.0, 0.0), pt(1.0, 1.0));
        o.add_edge(pt(1.0, 1.0), pt(0.0, 1.0));
        o.add_edge(pt(0.0, 1.0), pt(0.0, 0.0));
        o.perform();
        assert_eq!(o.status(), WireOrderStatus::Same);
        assert_eq!(o.ordered(1), 1);
        assert_eq!(o.ordered(2), 2);
        assert_eq!(o.ordered(3), 3);
        assert_eq!(o.ordered(4), 4);
    }

    /// A shuffled square must be reordered into the connected chain 1→4→2→3.
    #[test]
    fn shuffled_square_is_reordered() {
        // Store the four edges of a unit square out of order.
        // e1: (0,0)-(1,0); e2: (1,0)-(1,1); e3: (1,1)-(0,1); e4: (0,1)-(0,0)
        let mut o = WireOrder::new();
        o.add_edge(pt(0.0, 0.0), pt(1.0, 0.0)); // 1
        o.add_edge(pt(1.0, 1.0), pt(0.0, 1.0)); // 3
        o.add_edge(pt(1.0, 0.0), pt(1.0, 1.0)); // 2
        o.add_edge(pt(0.0, 1.0), pt(0.0, 0.0)); // 4
        o.perform();
        // Reordered (status 1) — the chain 1,2,3,4 maps back to edge numbers
        // 1,3,2,4 in stored order.
        assert_eq!(o.status(), WireOrderStatus::Reordered);
        let ordered: Vec<i32> = (1..=4).map(|i| o.ordered(i)).collect();
        assert_eq!(ordered, vec![1, 3, 2, 4]);
    }

    /// A single edge is its own (trivial) order.
    #[test]
    fn single_edge() {
        let mut o = WireOrder::new();
        o.add_edge(pt(0.0, 0.0), pt(1.0, 1.0));
        o.perform();
        assert_eq!(o.status(), WireOrderStatus::Same);
        assert_eq!(o.ordered(1), 1);
    }

    /// Empty order is a no-op.
    #[test]
    fn empty_order() {
        let mut o = WireOrder::new();
        o.perform();
        assert_eq!(o.nb_edges(), 0);
        assert_eq!(o.status(), WireOrderStatus::Same);
    }
}
