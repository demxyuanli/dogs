//! 2-D region tracing for the planar boolean arrangement.

use std::collections::{HashMap, HashSet};


use occt_core::geom::polygon_ops::{point_in_polygon2d, polygon_area2d};
use occt_core::gp::{GpPnt2d};



use crate::bop_builder_core::{
    disjoint_result, empty_result, validate, BoolOp, BooleanResult,
};







use crate::topo_tools_full::{
    edge_vertices, edges_of, edges_of_wire, faces_of, shapes_of, vertex_position, vertices_of,
    wires_of_face,
};


// ---------------------------------------------------------------------------
// 2-D planar arrangement (BOPAlgo_BuilderFace::PerformAreas-style region
// tracing)
// ---------------------------------------------------------------------------

/// A bounded region of the planar arrangement on one face: an outer loop
/// (CCW) plus hole loops (CW), with a point guaranteed inside the region.
pub(crate) struct Region2d {
    pub(crate) outer: Vec<GpPnt2d>,
    pub(crate) holes: Vec<Vec<GpPnt2d>>,
    pub(crate) interior: GpPnt2d,
}

/// 2-D spatial-hash point welder (near-coincident points → one index).
pub(crate) struct Weld2d {
    cell: f64,
    pts: Vec<GpPnt2d>,
    grid: HashMap<(i64, i64), Vec<usize>>,
}

impl Weld2d {
    fn new(cell: f64) -> Self {
        Self { cell: cell.max(1e-9), pts: Vec::new(), grid: HashMap::new() }
    }

    fn weld(&mut self, p: GpPnt2d) -> usize {
        let k = (f64::floor(p.x() / self.cell) as i64, f64::floor(p.y() / self.cell) as i64);
        let mut best: Option<(usize, f64)> = None;
        for dx in -1i64..=1 {
            for dy in -1i64..=1 {
                if let Some(bucket) = self.grid.get(&(k.0 + dx, k.1 + dy)) {
                    for &i in bucket {
                        let d = self.pts[i].distance(&p);
                        if d <= self.cell * 1.01 && best.map_or(true, |(_, bd)| d < bd) {
                            best = Some((i, d));
                        }
                    }
                }
            }
        }
        match best {
            Some((i, _)) => i,
            None => {
                let i = self.pts.len();
                self.pts.push(p);
                self.grid.entry(k).or_default().push(i);
                i
            }
        }
    }
}

/// Every intersection of segments `a` and `b` as `(t_a, t_b, point)`, where
/// `t_a`/`t_b` are the on-segment parameters (0..1). Handles proper crossings,
/// endpoint touches and collinear overlaps (whose two overlap endpoints are
/// reported).
pub(crate) fn seg_intersections(
    a0: GpPnt2d,
    a1: GpPnt2d,
    b0: GpPnt2d,
    b1: GpPnt2d,
) -> Vec<(f64, f64, GpPnt2d)> {
    let (dax, day) = (a1.x() - a0.x(), a1.y() - a0.y());
    let (dbx, dby) = (b1.x() - b0.x(), b1.y() - b0.y());
    let denom = dax * dby - day * dbx;
    let len2_a = dax * dax + day * day;
    let len2_b = dbx * dbx + dby * dby;
    if len2_a < 1e-24 || len2_b < 1e-24 {
        return Vec::new();
    }
    let eps = 1e-9 * len2_a.sqrt().max(len2_b.sqrt()).max(1.0);
    let proj_a = |p: GpPnt2d| ((p.x() - a0.x()) * dax + (p.y() - a0.y()) * day) / len2_a;
    let proj_b = |p: GpPnt2d| ((p.x() - b0.x()) * dbx + (p.y() - b0.y()) * dby) / len2_b;
    let (wax, way) = (b0.x() - a0.x(), b0.y() - a0.y());
    if denom.abs() > eps {
        let t_a = (wax * dby - way * dbx) / denom;
        let t_b = (wax * day - way * dax) / denom;
        let (lo, hi) = (-1e-7, 1.0 + 1e-7);
        if t_a >= lo && t_a <= hi && t_b >= lo && t_b <= hi {
            let ta = t_a.clamp(0.0, 1.0);
            let tb = t_b.clamp(0.0, 1.0);
            let p = GpPnt2d::new(a0.x() + dax * ta, a0.y() + day * ta);
            return vec![(ta, tb, p)];
        }
        return Vec::new();
    }
    // Parallel: report the overlap endpoints when collinear.
    if (dax * way - day * wax).abs() > eps {
        return Vec::new();
    }
    let a_lo = proj_a(b0).min(proj_a(b1)).max(0.0);
    let a_hi = proj_a(b0).max(proj_a(b1)).min(1.0);
    if a_hi - a_lo < 1e-9 {
        return Vec::new();
    }
    let p_lo = GpPnt2d::new(a0.x() + dax * a_lo, a0.y() + day * a_lo);
    let p_hi = GpPnt2d::new(a0.x() + dax * a_hi, a0.y() + day * a_hi);
    vec![
        (a_lo, proj_b(p_lo), p_lo),
        (a_hi, proj_b(p_hi), p_hi),
    ]
}

/// Split every segment of `segs` at its intersections with all the others,
/// returning the resulting sub-segments (collinear overlaps deduplicated by
/// construction).
pub(crate) fn split_segments_2d(segs: &[(GpPnt2d, GpPnt2d)]) -> Vec<(GpPnt2d, GpPnt2d)> {
    let n = segs.len();
    let mut cuts: Vec<Vec<(f64, GpPnt2d)>> = vec![Vec::new(); n];
    for i in 0..n {
        cuts[i].push((0.0, segs[i].0));
        cuts[i].push((1.0, segs[i].1));
    }
    for i in 0..n {
        let (a0, a1) = segs[i];
        for j in (i + 1)..n {
            let (b0, b1) = segs[j];
            let hits = seg_intersections(a0, a1, b0, b1);
            for (ta, tb, p) in hits {
                cuts[i].push((ta, p));
                cuts[j].push((tb, p));
            }
        }
    }
    let mut out: Vec<(GpPnt2d, GpPnt2d)> = Vec::new();
    for k in 0..n {
        let mut cs = cuts[k].clone();
        cs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut prev_t: Option<f64> = None;
        let mut last: Option<GpPnt2d> = None;
        for (t, p) in cs {
            if prev_t.map_or(false, |q| (t - q).abs() < 1e-9) {
                continue;
            }
            prev_t = Some(t);
            if let Some(lp) = last {
                if lp.distance(&p) > 1e-9 {
                    out.push((lp, p));
                }
            }
            last = Some(p);
        }
    }
    out
}

/// Whether `p` is strictly inside `poly` — the even-odd test AND not within
/// `tol` of any boundary edge. A point exactly on the boundary (e.g. a shared
/// loop corner) is not "inside".
pub(crate) fn point_strictly_inside(poly: &[GpPnt2d], p: &GpPnt2d, tol: f64) -> bool {
    if !point_in_polygon2d(poly, p) {
        return false;
    }
    let n = poly.len();
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let (abx, aby) = (b.x() - a.x(), b.y() - a.y());
        let (apx, apy) = (p.x() - a.x(), p.y() - a.y());
        let len2 = abx * abx + aby * aby;
        let t = if len2 < 1e-24 { 0.0 } else { ((apx * abx + apy * aby) / len2).clamp(0.0, 1.0) };
        let px = a.x() + abx * t;
        let py = a.y() + aby * t;
        if (p.x() - px).hypot(p.y() - py) < tol {
            return false;
        }
    }
    true
}

/// Interior point of a 2-D region (outer loop + holes): the outer centroid
/// when it lies in the region, else a point nudged off the outer loop's first
/// edge toward the region interior (left of a CCW loop).
pub(crate) fn region_interior2(outer: &[GpPnt2d], holes: &[Vec<GpPnt2d>]) -> GpPnt2d {
    let c = {
        let mut sx = 0.0;
        let mut sy = 0.0;
        for p in outer {
            sx += p.x();
            sy += p.y();
        }
        let inv = 1.0 / outer.len().max(1) as f64;
        GpPnt2d::new(sx * inv, sy * inv)
    };
    let inside = point_in_polygon2d(outer, &c) && holes.iter().all(|h| !point_in_polygon2d(h, &c));
    if inside {
        return c;
    }
    // Normalize the outer loop to CCW, then nudge the first edge's midpoint
    // toward the interior (left of a CCW edge).
    let mut o = outer.to_vec();
    if polygon_area2d(&o) < 0.0 {
        o.reverse();
    }
    let (a, b) = (o[0], o[1 % o.len()]);
    let (dx, dy) = (b.x() - a.x(), b.y() - a.y());
    let len = dx.hypot(dy);
    if len < 1e-12 {
        return c;
    }
    let m = GpPnt2d::new(0.5 * (a.x() + b.x()), 0.5 * (a.y() + b.y()));
    GpPnt2d::new(m.x() - dy / len * 1e-6, m.y() + dx / len * 1e-6)
}

/// Trace the bounded regions of the planar arrangement built from `boundary`
/// (the face's outer polygon) plus `segs2d` (the on-face section segments).
///
/// Mirrors `BOPAlgo_BuilderFace`'s `PerformLoops`/`PerformAreas`: the boundary
/// and section edges are split at every intersection and the resulting planar
/// subdivision's bounded faces are recovered by half-edge traversal. A region
/// that is a ring (outer boundary plus interior loops) keeps its holes — the
/// box-top annulus around a boss cylinder — while open section lines simply
/// partition the face into several simple regions.
pub(crate) fn trace_planar_regions(boundary: &[GpPnt2d], segs2d: &[(GpPnt2d, GpPnt2d)]) -> Vec<Region2d> {
    let bn = boundary.len();
    if bn < 3 {
        return Vec::new();
    }
    // Normalize the boundary to CCW (positive signed area).
    let mut boundary = boundary.to_vec();
    if polygon_area2d(&boundary) < 0.0 {
        boundary.reverse();
    }

    let mut raw: Vec<(GpPnt2d, GpPnt2d)> = Vec::new();
    for i in 0..bn {
        let (a, b) = (boundary[i], boundary[(i + 1) % bn]);
        if a.distance(&b) > 1e-9 {
            raw.push((a, b));
        }
    }
    for s in segs2d {
        if s.0.distance(&s.1) > 1e-9 {
            raw.push(*s);
        }
    }
    if raw.is_empty() {
        return Vec::new();
    }
    let sub = split_segments_2d(&raw);

    // Weld the split endpoints into shared vertex indices.
    let mut w2 = Weld2d::new(1e-8);
    let mut edges: Vec<(usize, usize)> = Vec::new();
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    for (a, b) in &sub {
        let ia = w2.weld(*a);
        let ib = w2.weld(*b);
        if ia == ib {
            continue;
        }
        let key = (ia.min(ib), ia.max(ib));
        if seen.insert(key) {
            edges.push((ia, ib));
        }
    }
    if edges.is_empty() {
        return Vec::new();
    }
    // OCCT `BOPAlgo_BuilderFace::PerformShapesToAvoid`: repeatedly strip edges
    // whose endpoint vertex is touched by a single edge — a dangling end that
    // can never close a loop. Such a segment cannot partition a region anyway
    // (it is a slit), so stripping it leaves the arrangement's regions intact
    // while keeping the half-edge traversal well-formed (even vertex degree).
    let mut keep = vec![true; edges.len()];
    loop {
        let mut deg: HashMap<usize, usize> = HashMap::new();
        for (i, &(a, b)) in edges.iter().enumerate() {
            if !keep[i] {
                continue;
            }
            *deg.entry(a).or_default() += 1;
            *deg.entry(b).or_default() += 1;
        }
        let mut stripped = false;
        for (i, &(a, b)) in edges.iter().enumerate() {
            if !keep[i] {
                continue;
            }
            if deg.get(&a).copied().unwrap_or(0) <= 1 || deg.get(&b).copied().unwrap_or(0) <= 1 {
                keep[i] = false;
                stripped = true;
            }
        }
        if !stripped {
            break;
        }
    }
    let edges: Vec<(usize, usize)> = edges
        .iter()
        .enumerate()
        .filter(|(i, _)| keep[*i])
        .map(|(_, e)| *e)
        .collect();
    if edges.is_empty() {
        return Vec::new();
    }
    let verts = &w2.pts;

    // Half-edges: index 2e is (u→v), 2e+1 is (v→u).
    let he = edges.len() * 2;
    let tail = |h: usize| if h % 2 == 0 { edges[h / 2].0 } else { edges[h / 2].1 };
    let head = |h: usize| if h % 2 == 0 { edges[h / 2].1 } else { edges[h / 2].0 };
    let mut out_at: HashMap<usize, Vec<usize>> = HashMap::new();
    for h in 0..he {
        out_at.entry(tail(h)).or_default().push(h);
    }
    let ang = |h: usize| {
        let (u, v) = (tail(h), head(h));
        (verts[v].y() - verts[u].y()).atan2(verts[v].x() - verts[u].x())
    };
    // next(h) = the outgoing half-edge at the head whose direction is the next
    // one CLOCKWISE from the reverse direction, tracing the face on the LEFT
    // of h (the incoming direction reversed, rotated CW to the first outgoing
    // edge). This keeps a straight boundary turn, a T-junction onto an on-face
    // line and a hole loop on the correct face.
    let next = |h: usize| -> usize {
        let v = head(h);
        let rev = ang(h ^ 1);
        let mut best: Option<usize> = None;
        let mut best_delta = std::f64::consts::TAU;
        for &h2 in out_at.get(&v).map(|l| l.as_slice()).unwrap_or(&[]) {
            if h2 == (h ^ 1) {
                continue;
            }
            let mut delta = rev - ang(h2);
            if delta <= 0.0 {
                delta += std::f64::consts::TAU;
            }
            if delta < best_delta {
                best_delta = delta;
                best = Some(h2);
            }
        }
        best.unwrap_or(h ^ 1)
    };

    // Trace all directed-edge cycles (each is one boundary loop).
    let mut used = vec![false; he];
    let mut cycles: Vec<(f64, Vec<GpPnt2d>)> = Vec::new(); // (signed area, pts)
    for h0 in 0..he {
        if used[h0] {
            continue;
        }
        let mut pts: Vec<GpPnt2d> = Vec::new();
        let mut h = h0;
        loop {
            if used[h] {
                break;
            }
            used[h] = true;
            pts.push(verts[tail(h)]);
            h = next(h);
            if h == h0 {
                break;
            }
        }
        if pts.len() >= 3 {
            cycles.push((polygon_area2d(&pts), pts));
        }
    }

    // Split into outer (CCW) and hole (CW) cycles, dropping the unbounded
    // face's cycle (the CW boundary with the largest |area|).
    let mut unbounded: Option<usize> = None;
    let mut unbounded_area = 0.0f64;
    for (i, (a, _)) in cycles.iter().enumerate() {
        if *a < 0.0 && a.abs() > unbounded_area {
            unbounded = Some(i);
            unbounded_area = a.abs();
        }
    }
    let mut regions: Vec<Region2d> = Vec::new();
    let mut outers: Vec<usize> = Vec::new(); // indices of CCW cycles
    for (i, (a, _)) in cycles.iter().enumerate() {
        if Some(i) == unbounded {
            continue;
        }
        if *a > 0.0 {
            outers.push(i);
        }
    }
    // Assign each CW (hole) cycle to the innermost CCW cycle STRICTLY
    // containing it. A point of the CW cycle that lies exactly on a CCW
    // cycle's boundary means the two are the same geometric loop traced with
    // opposite orientation (a ring's inner boundary vs its disk's outer
    // boundary) — the CW cycle is then NOT a hole of that disk.
    let mut holes_of: Vec<Vec<usize>> = vec![Vec::new(); outers.len()];
    for (i, (a, pts)) in cycles.iter().enumerate() {
        if Some(i) == unbounded || *a > 0.0 {
            continue;
        }
        let probe = pts[0];
        let mut best: Option<(usize, f64)> = None;
        for (oi, &o) in outers.iter().enumerate() {
            if point_strictly_inside(&cycles[o].1, &probe, 1e-9) {
                let aa = cycles[o].0.abs();
                if best.map_or(true, |(_, ba)| aa < ba) {
                    best = Some((oi, aa));
                }
            }
        }
        if let Some((oi, _)) = best {
            holes_of[oi].push(i);
        }
    }
    for (oi, &o) in outers.iter().enumerate() {
        let outer = cycles[o].1.clone();
        let mut holes: Vec<Vec<GpPnt2d>> = Vec::new();
        for &hi in &holes_of[oi] {
            holes.push(cycles[hi].1.clone());
        }
        let interior = region_interior2(&outer, &holes);
        regions.push(Region2d { outer, holes, interior });
    }
    regions
}
