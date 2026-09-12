//! Marching-squares / polyline chaining for general surface pairs.
//! Source: tracer previously in `intpatch.rs`.

use std::sync::Arc;

use occt_core::geom::polyline_simplify::rdp_simplify;
use occt_core::gp::{GpPnt, GpTrsf, GpVec};
use occt_geom::{Curve, Surface};

use super::geom::{distance_to_surface, project_params, sample_bounds};
use super::{IntersectionCurve, SurfaceIntersection};

// ---------------------------------------------------------------------------
// General tracer
// ---------------------------------------------------------------------------

/// A piecewise-linear `Curve` through sampled points (parameter `[0, 1]`).
/// Used to carry grid-traced intersections.
#[derive(Debug, Clone)]
pub struct PolylineCurve {
    pub pts: Vec<GpPnt>,
}

impl Curve for PolylineCurve {
    fn d0(&self, u: f64) -> GpPnt {
        let n = self.pts.len();
        if n == 0 {
            return GpPnt::zero();
        }
        if n == 1 {
            return self.pts[0];
        }
        let t = u.clamp(0.0, 1.0) * (n - 1) as f64;
        let i = (t.floor() as usize).min(n - 2);
        let f = t - i as f64;
        let a = self.pts[i];
        let b = self.pts[i + 1];
        GpPnt::new(a.x() + f * (b.x() - a.x()), a.y() + f * (b.y() - a.y()), a.z() + f * (b.z() - a.z()))
    }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let p0 = self.d0(u);
        let h = 1e-6;
        (p0, GpVec::from_pnts(&p0, &self.d0(u + h)).divided(h))
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (p, d1) = self.d1(u);
        (p, d1, GpVec::zero())
    }
    fn first_parameter(&self) -> f64 { 0.0 }
    fn last_parameter(&self) -> f64 { 1.0 }
    fn continuity(&self) -> u8 { 0 }
    fn transform(&mut self, t: &GpTrsf) {
        for p in &mut self.pts {
            *p = p.transformed(t);
        }
    }
    fn reverse(&mut self) { self.pts.reverse(); }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}

/// Crossing of a segment `a → b` (with signed field values `va, vb`) at the
/// zero level set.
fn edge_crossing(pa: &GpPnt, va: f64, pb: &GpPnt, vb: f64) -> Option<GpPnt> {
    if va * vb > 0.0 {
        return None;
    }
    if va.abs() < 1e-12 && vb.abs() < 1e-12 {
        return None;
    }
    if va.abs() < 1e-12 {
        return Some(*pa);
    }
    if vb.abs() < 1e-12 {
        return Some(*pb);
    }
    let t = va / (va - vb);
    Some(GpPnt::new(
        pa.x() + t * (pb.x() - pa.x()),
        pa.y() + t * (pb.y() - pa.y()),
        pa.z() + t * (pb.z() - pa.z()),
    ))
}

/// Chain a set of unordered 3D segments into polylines, merging endpoints that
/// are within `tol`, and return the flattened point stream.
fn chain_segments(segs: Vec<(GpPnt, GpPnt)>, tol: f64) -> Vec<GpPnt> {
    let mut remaining = segs;
    let mut out: Vec<GpPnt> = Vec::new();
    while !remaining.is_empty() {
        let (a, b) = remaining.remove(0);
        let mut poly = vec![a, b];
        let mut changed = true;
        while changed {
            changed = false;
            for i in (0..remaining.len()).rev() {
                let (c, d) = remaining[i];
                let first = poly[0];
                let last = *poly.last().unwrap();
                if c.distance(&last) <= tol {
                    poly.push(d);
                    remaining.swap_remove(i);
                    changed = true;
                } else if d.distance(&last) <= tol {
                    poly.push(c);
                    remaining.swap_remove(i);
                    changed = true;
                } else if c.distance(&first) <= tol {
                    poly.insert(0, d);
                    remaining.swap_remove(i);
                    changed = true;
                } else if d.distance(&first) <= tol {
                    poly.insert(0, c);
                    remaining.swap_remove(i);
                    changed = true;
                }
            }
        }
        out.extend(poly);
    }
    out
}

/// General fallback intersection: grid-march `a`'s parameter space and keep the
/// contour where the distance from `a(u, v)` to surface `b` equals `tol`.
///
/// Returns the traced 3D points (concatenated polylines). The level-set
/// contour is extracted with a marching-squares pass over the `(u, v)` grid.
pub fn trace_surface_curve(a: &dyn Surface, b: &dyn Surface, tol: f64) -> Vec<GpPnt> {
    trace_surface_curve_n(a, b, tol, 48, 48)
}

/// `trace_surface_curve` with an explicit marching-squares grid size.
fn trace_surface_curve_n(a: &dyn Surface, b: &dyn Surface, tol: f64, nu: usize, nv: usize) -> Vec<GpPnt> {
    trace_zero_set(a, |p| distance_to_surface(p, b) - tol, nu, nv, tol * 2.0)
}

/// March the zero contour of `field` over `a`'s UV box (`IntPatch_TheSurfFunction`
/// zeros when `field` is signed quadric distance).
pub(crate) fn trace_zero_set(
    a: &dyn Surface,
    field: impl Fn(&GpPnt) -> f64,
    nu: usize,
    nv: usize,
    chain_tol: f64,
) -> Vec<GpPnt> {
    let (u0, u1, v0, v1) = sample_bounds(a);
    let (nu, nv) = (nu.max(4), nv.max(4));
    let mut grid = vec![vec![0.0f64; nv + 1]; nu + 1];
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let p = a.d0(u, v);
            grid[i][j] = field(&p);
        }
    }

    let mut segments: Vec<(GpPnt, GpPnt)> = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            let p = [
                a.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64),
                a.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64),
                a.d0(u0 + (u1 - u0) * (i + 1) as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64),
                a.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * (j + 1) as f64 / nv as f64),
            ];
            let f = [grid[i][j], grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1]];
            let e01 = edge_crossing(&p[0], f[0], &p[1], f[1]);
            let e12 = edge_crossing(&p[1], f[1], &p[2], f[2]);
            let e23 = edge_crossing(&p[2], f[2], &p[3], f[3]);
            let e30 = edge_crossing(&p[3], f[3], &p[0], f[0]);
            let mut pts = Vec::new();
            for e in [e01, e12, e23, e30] {
                if let Some(q) = e {
                    pts.push(q);
                }
            }
            if pts.len() == 2 {
                segments.push((pts[0], pts[1]));
            } else if pts.len() >= 4 {
                segments.push((pts[0], pts[2]));
                segments.push((pts[1], pts[3]));
            }
        }
    }
    chain_segments(segments, chain_tol)
}

/// Chain an unordered set of 3D points into ordered polylines by greedy
/// nearest-neighbour walking from both ends of each chain. Points separated by
/// more than a (diagonal-relative) link distance start a new chain; junctions
/// are left as separate chains rather than forcing a single path through them.
///
/// The link tolerance is `max(tol, 2% of the bounding-box diagonal)`, so chains
/// whose overlapping endpoints lie within `tol` merge into a single polyline
/// while well-separated disjoint curves stay separate.
pub fn chain_intersection_points(pts: &[GpPnt], tol: f64) -> Vec<Vec<GpPnt>> {
    if pts.is_empty() {
        return Vec::new();
    }
    let mut lo = pts[0];
    let mut hi = pts[0];
    for p in pts {
        lo = GpPnt::new(lo.x().min(p.x()), lo.y().min(p.y()), lo.z().min(p.z()));
        hi = GpPnt::new(hi.x().max(p.x()), hi.y().max(p.y()), hi.z().max(p.z()));
    }
    let link_tol = tol.max(lo.distance(&hi) * 0.02);

    let mut used = vec![false; pts.len()];
    let mut chains: Vec<Vec<GpPnt>> = Vec::new();
    loop {
        let Some(start) = (0..pts.len()).find(|&i| !used[i]) else { break };
        used[start] = true;
        let mut chain: Vec<GpPnt> = vec![pts[start]];
        loop {
            let tail = *chain.last().unwrap();
            let head = chain[0];
            let mut best: Option<(f64, usize, bool)> = None; // (dist, idx, at_tail)
            for (i, p) in pts.iter().enumerate() {
                if used[i] {
                    continue;
                }
                let dt = p.distance(&tail);
                let dh = p.distance(&head);
                if best.map_or(true, |(bd, _, _)| dt < bd || dh < bd) {
                    if dt <= dh {
                        best = Some((dt, i, true));
                    } else {
                        best = Some((dh, i, false));
                    }
                }
            }
            match best {
                Some((d, i, at_tail)) if d <= link_tol => {
                    used[i] = true;
                    if at_tail {
                        chain.push(pts[i]);
                    } else {
                        chain.insert(0, pts[i]);
                    }
                }
                _ => break,
            }
        }
        chains.push(chain);
    }
    chains
}

/// Convert a chained intersection polyline into an [`IntersectionCurve`].
///
/// Long polylines (> 8 points) are first simplified with the
/// Ramer–Douglas–Peucker algorithm at `tol`; the resulting vertices become a
/// [`PolylineCurve`], and each vertex is projected onto both surfaces to fill
/// `on_a` / `on_b`.
pub fn polyline_to_curve(poly: &[GpPnt], a: &dyn Surface, b: &dyn Surface, tol: f64) -> Option<IntersectionCurve> {
    if poly.len() < 2 {
        return None;
    }
    let pts: Vec<GpPnt> = if poly.len() > 8 {
        let keep = rdp_simplify(poly, tol);
        keep.iter().map(|&i| poly[i]).collect()
    } else {
        poly.to_vec()
    };
    if pts.len() < 2 {
        return None;
    }
    let curve: Arc<dyn Curve> = Arc::new(PolylineCurve { pts: pts.clone() });
    let on_a: Vec<(f64, f64)> = pts.iter().map(|p| project_params(a, p)).collect();
    let on_b: Vec<(f64, f64)> = pts.iter().map(|p| project_params(b, p)).collect();
    Some(IntersectionCurve { curve, points: pts, on_a, on_b })
}

/// Trace the intersection curves of two general surfaces and return them as
/// ordered polylines (each `Vec<GpPnt>` is one chained curve). Exposed for the
/// curved boolean (`bop_curved`) which needs per-face intersection polylines.
pub fn intersection_curve_points(a: &dyn Surface, b: &dyn Surface, tol: f64, samples: usize) -> Vec<Vec<GpPnt>> {
    let grid = samples.clamp(24, 96);
    let pts = trace_surface_curve_n(a, b, tol, grid, grid);
    if pts.is_empty() {
        return Vec::new();
    }
    chain_intersection_points(&pts, tol)
}

/// Intersect two general (possibly non-analytic) surfaces using the grid
/// tracer, chaining the traced points into polylines and wrapping each in an
/// [`IntersectionCurve`].
pub fn intersect_general_surfaces(a: &dyn Surface, b: &dyn Surface, tol: f64) -> SurfaceIntersection {
    let pts = trace_surface_curve(a, b, tol);
    if pts.is_empty() {
        return SurfaceIntersection::None;
    }
    let chains = chain_intersection_points(&pts, tol);
    let mut curves = Vec::new();
    for chain in &chains {
        if let Some(ic) = polyline_to_curve(chain, a, b, tol) {
            curves.push(ic);
        }
    }
    if curves.is_empty() {
        SurfaceIntersection::None
    } else {
        SurfaceIntersection::Curves(curves)
    }
}

/// Quick overlap test: `true` when the grid tracer finds any intersection
/// point between the two surfaces.
pub fn surfaces_intersect_general(a: &dyn Surface, b: &dyn Surface, tol: f64) -> bool {
    !trace_surface_curve(a, b, tol).is_empty()
}
