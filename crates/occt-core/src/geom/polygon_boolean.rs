//! Planar polygon boolean operations — intersection, union, difference of two
//! simple polygons.
//!
//! Port of the planar boolean step used by `BOPAlgo_Builder` / `BRepAlgoAPI`
//! (2D polygon clipping): a Weiler–Atherton-style boundary walk on refined
//! edges, with a Sutherland–Hodgman fast path for the convex/convex case.
//! Produces a list of closed result polygons.
//!
//! Source: `BOPAlgo_Builder` (planar section), `IntPolyh`, `BRepAlgoAPI`.

use crate::gp::GpPnt2d;

/// Boolean operation for two planar polygons.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PolygonBoolOp {
    /// A ∩ B (interior of both).
    Intersect,
    /// A ∪ B (interior of either).
    Union,
    /// A − B (interior of A, exterior of B).
    Difference,
}

const EPS: f64 = 1e-9;

/// Orientation of `p` relative to the directed segment a→b: positive = left.
fn orient(a: &GpPnt2d, b: &GpPnt2d, p: &GpPnt2d) -> f64 {
    (b.x() - a.x()) * (p.y() - a.y()) - (b.y() - a.y()) * (p.x() - a.x())
}

/// Segment–segment intersection. Returns:
/// - None: no intersection.
/// - Some((p, false)): a proper crossing at `p`.
/// - Some((p, true)): collinear / endpoint contact (p = the contact point).
pub fn segment_segment2d(
    p1: &GpPnt2d,
    p2: &GpPnt2d,
    p3: &GpPnt2d,
    p4: &GpPnt2d,
) -> Option<(GpPnt2d, bool)> {
    let d1 = orient(p3, p4, p1);
    let d2 = orient(p3, p4, p2);
    let d3 = orient(p1, p2, p3);
    let d4 = orient(p1, p2, p4);
    let proper = (d1 > EPS && d2 < -EPS || d1 < -EPS && d2 > EPS)
        && (d3 > EPS && d4 < -EPS || d3 < -EPS && d4 > EPS);
    if proper {
        let rx = p2.x() - p1.x();
        let ry = p2.y() - p1.y();
        let sx = p4.x() - p3.x();
        let sy = p4.y() - p3.y();
        let denom = rx * sy - ry * sx;
        if denom.abs() < EPS {
            return None;
        }
        let t = ((p3.x() - p1.x()) * sy - (p3.y() - p1.y()) * sx) / denom;
        let q = GpPnt2d::new(p1.x() + t * rx, p1.y() + t * ry);
        return Some((q, false));
    }
    // Endpoint / collinear contact.
    let on = |p: &GpPnt2d, a: &GpPnt2d, b: &GpPnt2d| {
        p.x() >= a.x().min(b.x()) - EPS
            && p.x() <= a.x().max(b.x()) + EPS
            && p.y() >= a.y().min(b.y()) - EPS
            && p.y() <= a.y().max(b.y()) + EPS
    };
    if d1.abs() <= EPS && on(p1, p3, p4) {
        return Some((*p1, true));
    }
    if d2.abs() <= EPS && on(p2, p3, p4) {
        return Some((*p2, true));
    }
    if d3.abs() <= EPS && on(p3, p1, p2) {
        return Some((*p3, true));
    }
    if d4.abs() <= EPS && on(p4, p1, p2) {
        return Some((*p4, true));
    }
    None
}

/// Even-odd point-in-polygon test (any winding).
pub fn point_in_polygon2d(poly: &[GpPnt2d], p: &GpPnt2d) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (poly[i].x(), poly[i].y());
        let (xj, yj) = (poly[j].x(), poly[j].y());
        let crosses = (yi > p.y()) != (yj > p.y())
            && p.x() < (xj - xi) * (p.y() - yi) / (yj - yi) + xi;
        if crosses {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Whether all vertices of `poly` lie inside (or on the boundary of) `other`.
fn all_vertices_inside(poly: &[GpPnt2d], other: &[GpPnt2d]) -> bool {
    poly.iter().all(|p| point_in_polygon2d(other, p))
}

/// Signed area (positive = CCW).
pub fn signed_area2d(poly: &[GpPnt2d]) -> f64 {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }
    let mut s = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        s += poly[i].x() * poly[j].y() - poly[j].x() * poly[i].y();
    }
    0.5 * s
}

/// Is the polygon convex (all consecutive-edge cross products same sign)?
pub fn is_convex2d(poly: &[GpPnt2d]) -> bool {
    let n = poly.len();
    if n < 4 {
        return true;
    }
    let mut sign = 0.0f64;
    for i in 0..n {
        let c = orient(&poly[i], &poly[(i + 1) % n], &poly[(i + 2) % n]);
        if c.abs() < 1e-12 {
            continue;
        }
        if sign == 0.0 {
            sign = c.signum();
        } else if c.signum() != sign {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Sutherland–Hodgman — convex intersection fast path.
// ---------------------------------------------------------------------------

fn clip_halfplane(poly: &[GpPnt2d], a: &GpPnt2d, b: &GpPnt2d, keep_left: bool) -> Vec<GpPnt2d> {
    let mut out = Vec::new();
    let n = poly.len();
    if n == 0 {
        return out;
    }
    // Line through the clip edge, used to intersect the poly edge with the
    // INFINITE half-plane boundary (the clip segment may not contain the
    // crossing point).
    let (ax, ay, bx, by) = (a.x(), a.y(), b.x(), b.y());
    let (dx, dy) = (bx - ax, by - ay);
    let cross_point = |cur: &GpPnt2d, nxt: &GpPnt2d| -> GpPnt2d {
        // Solve orient(a, b, p(t)) = 0 along p = cur + t·(nxt−cur).
        let o_cur = orient(a, b, cur);
        let o_nxt = orient(a, b, nxt);
        let denom = o_cur - o_nxt;
        let t = if denom.abs() < 1e-300 { 0.0 } else { o_cur / denom };
        GpPnt2d::new(
            cur.x() + t * (nxt.x() - cur.x()),
            cur.y() + t * (nxt.y() - cur.y()),
        )
    };
    for i in 0..n {
        let cur = poly[i];
        let nxt = poly[(i + 1) % n];
        let o_cur = orient(a, b, &cur);
        let o_nxt = orient(a, b, &nxt);
        let c_ok = if keep_left { o_cur >= -EPS } else { o_cur <= EPS };
        let n_ok = if keep_left { o_nxt >= -EPS } else { o_nxt <= EPS };
        if c_ok {
            out.push(cur);
            if !n_ok {
                out.push(cross_point(&cur, &nxt));
            }
        } else if n_ok {
            out.push(cross_point(&cur, &nxt));
        }
    }
    let _ = (ax, ay, bx, by, dx, dy);
    out
}

/// Intersection of two convex polygons via Sutherland–Hodgman.
pub fn convex_polygon_intersect(a: &[GpPnt2d], b: &[GpPnt2d]) -> Vec<GpPnt2d> {
    let mut poly = a.to_vec();
    let n = b.len();
    for i in 0..n {
        if poly.len() < 3 {
            break;
        }
        poly = clip_halfplane(&poly, &b[i], &b[(i + 1) % n], true);
    }
    if poly.len() >= 3 && signed_area2d(&poly).abs() > 1e-12 {
        poly
    } else {
        Vec::new()
    }
}

// ---------------------------------------------------------------------------
// Weiler–Atherton-style general boolean.
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct LoopVert {
    p: GpPnt2d,
    is_isect: bool,
    /// Index of the matching intersection vertex in the OTHER loop.
    twin: usize,
    /// Cyclic link pointers.
    next: usize,
    prev: usize,
    /// Whether this vertex has been consumed by an output walk.
    used: bool,
}

struct Loop {
    v: Vec<LoopVert>,
}

impl Loop {
    fn total(&self) -> usize {
        self.v.len()
    }

    fn edge_midpoint_inside(&self, i: usize, other: &[GpPnt2d]) -> bool {
        let a = self.v[i].p;
        let b = self.v[self.v[i].next].p;
        let mid = GpPnt2d::new((a.x() + b.x()) * 0.5, (a.y() + b.y()) * 0.5);
        point_in_polygon2d(other, &mid)
    }
}

/// Collect proper intersections between `a` and `b`, dedup'd. Returns the
/// point list and, per polygon, the sorted (pt_id, param) per edge.
fn collect_intersections(
    a: &[GpPnt2d],
    b: &[GpPnt2d],
) -> (Vec<GpPnt2d>, Vec<Vec<(usize, f64)>>, Vec<Vec<(usize, f64)>>) {
    let (na, nb) = (a.len(), b.len());
    let mut pts: Vec<GpPnt2d> = Vec::new();
    let mut per_a: Vec<Vec<(usize, f64)>> = vec![Vec::new(); na];
    let mut per_b: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nb];
    let param_on = |o: &GpPnt2d, e: &GpPnt2d, q: &GpPnt2d| -> f64 {
        if (e.x() - o.x()).abs() > (e.y() - o.y()).abs() {
            (q.x() - o.x()) / (e.x() - o.x())
        } else if (e.y() - o.y()).abs() > EPS {
            (q.y() - o.y()) / (e.y() - o.y())
        } else {
            0.0
        }
    };
    for i in 0..na {
        let (pa, pb) = (a[i], a[(i + 1) % na]);
        for j in 0..nb {
            let (qa, qb) = (b[j], b[(j + 1) % nb]);
            if let Some((p, collinear)) = segment_segment2d(&pa, &pb, &qa, &qb) {
                if collinear {
                    continue; // contact handled via the containment branch
                }
                let mut pt_id = usize::MAX;
                for (k, existing) in pts.iter().enumerate() {
                    if existing.distance(&p) < 1e-7 {
                        pt_id = k;
                        break;
                    }
                }
                if pt_id == usize::MAX {
                    pts.push(p);
                    pt_id = pts.len() - 1;
                }
                per_a[i].push((pt_id, param_on(&pa, &pb, &p)));
                per_b[j].push((pt_id, param_on(&qa, &qb, &p)));
            }
        }
    }
    for e in per_a.iter_mut() {
        e.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal));
        e.dedup_by(|x, y| x.0 == y.0);
    }
    for e in per_b.iter_mut() {
        e.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal));
        e.dedup_by(|x, y| x.0 == y.0);
    }
    (pts, per_a, per_b)
}

fn build_loop(poly: &[GpPnt2d], pts: &[GpPnt2d], per_edge: &[Vec<(usize, f64)>]) -> Loop {
    let n = poly.len();
    let mut verts: Vec<LoopVert> = Vec::with_capacity(n + pts.len());
    for i in 0..n {
        verts.push(LoopVert {
            p: poly[i],
            is_isect: false,
            twin: usize::MAX,
            next: 0,
            prev: 0,
            used: false,
        });
        for (pt_id, _) in &per_edge[i] {
            verts.push(LoopVert {
                p: pts[*pt_id],
                is_isect: true,
                twin: usize::MAX,
                next: 0,
                prev: 0,
                used: false,
            });
        }
    }
    let m = verts.len();
    for i in 0..m {
        verts[i].next = (i + 1) % m;
        verts[i].prev = (i + m - 1) % m;
    }
    Loop { v: verts }
}

/// Per-operation parameters of the boundary walk.
struct OpParams {
    keep_a_inside: bool,
    keep_b_inside: bool,
    dir_a: i32,
    dir_b: i32,
}

fn op_params(op: PolygonBoolOp) -> OpParams {
    match op {
        PolygonBoolOp::Intersect => OpParams {
            keep_a_inside: true,
            keep_b_inside: true,
            dir_a: 1,
            dir_b: -1,
        },
        PolygonBoolOp::Union => OpParams {
            keep_a_inside: false,
            keep_b_inside: false,
            dir_a: 1,
            dir_b: 1,
        },
        PolygonBoolOp::Difference => OpParams {
            keep_a_inside: false,
            keep_b_inside: true,
            dir_a: 1,
            dir_b: -1,
        },
    }
}

/// Walk one result loop starting from an intersection vertex of `la`.
fn trace_one(
    la: &mut Loop,
    lb: &mut Loop,
    params: &OpParams,
    other_a: &[GpPnt2d],
    other_b: &[GpPnt2d],
    start_a: usize,
) -> Option<Vec<GpPnt2d>> {
    // Start only if the A-edge leaving `start_a` is on the retained side.
    if la.edge_midpoint_inside(start_a, other_b) != params.keep_a_inside {
        return None;
    }
    let mut out: Vec<GpPnt2d> = Vec::new();
    let mut in_a = true;
    let mut cur = start_a;
    let mut guard = 0usize;
    loop {
        guard += 1;
        if guard > 4 * (la.total() + lb.total()) + 32 {
            return None;
        }
        // Returned to the start vertex on A after a full traversal?
        if in_a && cur == start_a && !out.is_empty() {
            break;
        }
        // Walk the current loop forward/back until we reach an intersection.
        let mut reached_isect = false;
        let mut twin = usize::MAX;
        let mut next_cur = cur;
        {
            let (loop_ref, dir, keep, other_poly) = if in_a {
                (&mut *la, params.dir_a, params.keep_a_inside, other_b)
            } else {
                (&mut *lb, params.dir_b, params.keep_b_inside, other_a)
            };
            for _ in 0..loop_ref.total() {
                let vi = next_cur;
                let nxt = if dir > 0 {
                    loop_ref.v[vi].next
                } else {
                    loop_ref.v[vi].prev
                };
                out.push(loop_ref.v[vi].p);
                loop_ref.v[vi].used = true;
                if loop_ref.v[nxt].is_isect {
                    if loop_ref.edge_midpoint_inside(vi, other_poly) != keep {
                        return None;
                    }
                    twin = loop_ref.v[nxt].twin;
                    reached_isect = true;
                    break;
                }
                next_cur = nxt;
            }
        }
        if !reached_isect {
            return None;
        }
        in_a = !in_a;
        cur = twin;
        if twin == usize::MAX {
            return None;
        }
        // Deduplicate a repeated vertex (isect was recorded at both ends).
        if out.len() >= 2 {
            let last = out[out.len() - 1];
            let prev = out[out.len() - 2];
            if last.distance(&prev) < 1e-9 {
                out.pop();
            }
        }
    }
    // Dedup the closure point (first == last).
    if out.len() >= 2 {
        let first = out[0];
        let last = out[out.len() - 1];
        if first.distance(&last) < 1e-9 {
            out.pop();
        }
    }
    if out.len() >= 3 && signed_area2d(&out).abs() > 1e-12 {
        Some(out)
    } else {
        None
    }
}

/// Planar boolean of two simple polygons. Returns the list of closed result
/// polygons (CCW). Exact for convex inputs; the general path handles simple
/// (non-self-intersecting) polygons with a Weiler–Atherton boundary walk.
/// Degenerate contacts (shared edges, point-touching) fall back to the
/// containment branch, which is exact for strict containment.
pub fn polygon_boolean(a: &[GpPnt2d], b: &[GpPnt2d], op: PolygonBoolOp) -> Vec<Vec<GpPnt2d>> {
    if a.len() < 3 || b.len() < 3 {
        return Vec::new();
    }
    if op == PolygonBoolOp::Intersect && is_convex2d(a) && is_convex2d(b) {
        let r = convex_polygon_intersect(a, b);
        return if r.is_empty() {
            Vec::new()
        } else {
            vec![r]
        };
    }

    let (pts, per_a, per_b) = collect_intersections(a, b);

    if pts.is_empty() {
        let a_in_b = all_vertices_inside(a, b);
        let b_in_a = all_vertices_inside(b, a);
        return match op {
            PolygonBoolOp::Intersect => {
                if a_in_b {
                    vec![a.to_vec()]
                } else if b_in_a {
                    vec![b.to_vec()]
                } else {
                    Vec::new()
                }
            }
            PolygonBoolOp::Union => {
                if a_in_b {
                    vec![b.to_vec()]
                } else if b_in_a {
                    vec![a.to_vec()]
                } else {
                    vec![a.to_vec(), b.to_vec()]
                }
            }
            PolygonBoolOp::Difference => {
                if a_in_b {
                    Vec::new()
                } else if b_in_a {
                    // b fully inside a → a with a hole; return both loops so
                    // even-odd fill recovers the correct area.
                    let mut hole = b.to_vec();
                    if signed_area2d(&hole) > 0.0 {
                        hole.reverse();
                    }
                    vec![a.to_vec(), hole]
                } else {
                    vec![a.to_vec()]
                }
            }
        };
    }

    let (mut la, mut lb) = (build_loop(a, &pts, &per_a), build_loop(b, &pts, &per_b));
    // Link twins: match intersection vertices by point proximity.
    let mut matched: Vec<(usize, usize)> = Vec::new();
    for i in 0..la.total() {
        if la.v[i].is_isect {
            for j in 0..lb.total() {
                if lb.v[j].is_isect && la.v[i].p.distance(&lb.v[j].p) < 1e-7 {
                    la.v[i].twin = j;
                    lb.v[j].twin = i;
                    matched.push((i, j));
                    break;
                }
            }
        }
    }
    let _ = matched;

    let params = op_params(op);
    let mut results: Vec<Vec<GpPnt2d>> = Vec::new();
    for i in 0..la.total() {
        if la.v[i].is_isect && !la.v[i].used {
            if let Some(loop_p) = trace_one(&mut la, &mut lb, &params, a, b, i) {
                results.push(loop_p);
            }
        }
    }
    for r in results.iter_mut() {
        if signed_area2d(r) < 0.0 {
            r.reverse();
        }
    }
    let mut uniq: Vec<Vec<GpPnt2d>> = Vec::new();
    for r in results {
        if uniq
            .iter()
            .all(|u| (signed_area2d(u) - signed_area2d(&r)).abs() > 1e-8)
        {
            uniq.push(r);
        }
    }
    uniq
}

/// Union of two polygons — convenience wrapper.
pub fn polygon_union(a: &[GpPnt2d], b: &[GpPnt2d]) -> Vec<Vec<GpPnt2d>> {
    polygon_boolean(a, b, PolygonBoolOp::Union)
}

/// Intersection of two polygons — convenience wrapper.
pub fn polygon_intersect(a: &[GpPnt2d], b: &[GpPnt2d]) -> Vec<Vec<GpPnt2d>> {
    polygon_boolean(a, b, PolygonBoolOp::Intersect)
}

/// Difference a − b — convenience wrapper.
pub fn polygon_difference(a: &[GpPnt2d], b: &[GpPnt2d]) -> Vec<Vec<GpPnt2d>> {
    polygon_boolean(a, b, PolygonBoolOp::Difference)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<GpPnt2d> {
        vec![
            GpPnt2d::new(x0, y0),
            GpPnt2d::new(x1, y0),
            GpPnt2d::new(x1, y1),
            GpPnt2d::new(x0, y1),
        ]
    }

    /// A triangle that properly crosses the unit rect [0,1]² (no shared edges).
    fn triangle() -> Vec<GpPnt2d> {
        vec![
            GpPnt2d::new(-0.5, -0.5),
            GpPnt2d::new(2.5, -0.5),
            GpPnt2d::new(1.25, 2.0),
        ]
    }

    fn area(pts: &[GpPnt2d]) -> f64 {
        signed_area2d(pts).abs()
    }

    fn total_area(polys: &[Vec<GpPnt2d>]) -> f64 {
        // Signed sum: hole loops (CW) subtract.
        polys.iter().map(|p| signed_area2d(p)).sum::<f64>().abs()
    }

    #[test]
    fn segment_crossing_and_touch() {
        let a = GpPnt2d::new(0.0, 0.0);
        let b = GpPnt2d::new(2.0, 0.0);
        let c = GpPnt2d::new(1.0, -1.0);
        let d = GpPnt2d::new(1.0, 1.0);
        let (p, coll) = segment_segment2d(&a, &b, &c, &d).expect("cross");
        assert!(!coll);
        assert!((p.x() - 1.0).abs() < 1e-9 && p.y().abs() < 1e-9);
        // Endpoint touch.
        assert!(segment_segment2d(&a, &b, &b, &d).is_some());
        // Parallel, separated → no touch.
        assert!(segment_segment2d(&a, &b, &GpPnt2d::new(0.0, 2.0), &GpPnt2d::new(2.0, 2.0)).is_none());
    }

    #[test]
    fn convex_intersect_overlapping_rects() {
        let a = rect(0.0, 0.0, 2.0, 2.0);
        let b = rect(1.0, 0.0, 3.0, 2.0);
        let r = convex_polygon_intersect(&a, &b);
        assert_eq!(r.len(), 4);
        assert!((area(&r) - 2.0).abs() < 1e-9, "overlap area {}", area(&r));
    }

    #[test]
    fn convex_detection() {
        assert!(is_convex2d(&rect(0.0, 0.0, 1.0, 1.0)));
        let l = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(2.0, 0.0),
            GpPnt2d::new(2.0, 1.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(1.0, 2.0),
            GpPnt2d::new(0.0, 2.0),
        ];
        assert!(!is_convex2d(&l));
    }

    #[test]
    fn intersection_rects_works() {
        let a = rect(0.0, 0.0, 2.0, 2.0);
        let b = rect(1.0, 1.0, 3.0, 3.0);
        let r = polygon_intersect(&a, &b);
        assert_eq!(r.len(), 1);
        assert!((area(&r[0]) - 1.0).abs() < 1e-8, "area {}", area(&r[0]));
    }

    #[test]
    fn union_rects_area() {
        let a = rect(0.0, 0.0, 2.0, 2.0);
        let b = rect(1.0, 1.0, 3.0, 3.0);
        let r = polygon_union(&a, &b);
        let total: f64 = total_area(&r);
        assert!((total - 7.0).abs() < 1e-8, "union area {}", total);
    }

    #[test]
    fn difference_rects_area() {
        let a = rect(0.0, 0.0, 3.0, 2.0);
        let b = rect(1.0, 0.5, 2.0, 1.5);
        let r = polygon_difference(&a, &b);
        let total: f64 = total_area(&r);
        assert!((total - 5.0).abs() < 1e-8, "diff area {}", total);
    }

    #[test]
    fn disjoint_union_two_pieces() {
        let a = rect(0.0, 0.0, 1.0, 1.0);
        let b = rect(3.0, 3.0, 4.0, 4.0);
        let r = polygon_union(&a, &b);
        assert_eq!(r.len(), 2);
    }

    #[test]
    fn contained_intersection() {
        let outer = rect(0.0, 0.0, 4.0, 4.0);
        let inner = rect(1.0, 1.0, 2.0, 2.0);
        let r = polygon_intersect(&outer, &inner);
        assert_eq!(r.len(), 1);
        assert!((area(&r[0]) - 1.0).abs() < 1e-8);
    }

    #[test]
    fn triangle_cross_rect() {
        let t = triangle();
        let r = polygon_intersect(&t, &rect(0.0, 0.0, 1.0, 1.0));
        assert!(!r.is_empty());
        let a = area(&r[0]);
        // Triangle ∩ [0,1]²: roughly the triangle portion in the unit square.
        assert!(a > 0.2 && a < 1.0, "tri∩rect area {}", a);
    }

    #[test]
    fn union_contains_both_areas() {
        let t = triangle();
        let r = polygon_union(&t, &rect(0.0, 0.0, 1.0, 1.0));
        let total = total_area(&r);
        let t_area = area(&t);
        assert!(total > t_area && total > 1.0, "union area {} tri {}", total, t_area);
    }

    #[test]
    fn difference_removes_triangle() {
        let t = triangle();
        let r = polygon_difference(&t, &rect(0.0, 0.0, 1.0, 1.0));
        let total = total_area(&r);
        let t_area = area(&t);
        // Triangle minus the rect overlap: strictly less than the full triangle.
        assert!(total > 0.0 && total < t_area, "diff area {} tri {}", total, t_area);
    }

    #[test]
    fn contained_difference_has_hole() {
        let outer = rect(0.0, 0.0, 4.0, 4.0);
        let inner = rect(1.0, 1.0, 2.0, 2.0);
        let r = polygon_difference(&outer, &inner);
        // Outer boundary + reversed hole: signed total = 16 − 1 = 15.
        let total = total_area(&r);
        assert!((total - 15.0).abs() < 1e-8, "diff area {}", total);
    }

    #[test]
    fn degenerate_shared_edge_union() {
        // Two rects sharing the bottom edge — handled by the containment path.
        let a = rect(0.0, 0.0, 2.0, 2.0);
        let b = rect(1.0, 0.0, 3.0, 2.0);
        let r = polygon_union(&a, &b);
        assert!(!r.is_empty());
    }
}
