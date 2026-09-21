use super::prelude::*;
use super::*;

/// Append an edge's UV points to a boundary chain, reversing the edge when
/// needed so the chain stays a continuous closed polygon.
pub(super) fn stitch_chain(chain: &mut Vec<GpPnt2d>, mut edge_uv: Vec<GpPnt2d>) {
    if edge_uv.is_empty() {
        return;
    }
    if let Some(&last) = chain.last() {
        let first = edge_uv[0];
        let last_edge = edge_uv[edge_uv.len() - 1];
        // The edge's own last point joining the chain's last means the edge is
        // stored in the traversal direction already; otherwise flip it.
        if last.distance(&last_edge) < last.distance(&first) {
            edge_uv.reverse();
        }
    }
    for p in edge_uv {
        if chain.last().map_or(true, |q: &GpPnt2d| q.distance(&p) > 1e-9) {
            chain.push(p);
        }
    }
}

/// Whether a boundary chain is a closed loop (its last point equals its first).
pub(super) fn chain_closed(chain: &[GpPnt2d]) -> bool {
    chain.len() >= 2 && chain[0].distance(&chain[chain.len() - 1]) < 1e-6
}

/// Invert a surface point to its `(u, v)` parameters via Newton iteration,
/// seeded from the surface range center (or the origin for unbounded ranges).
///
/// For analytic surfaces (planes, cylinders, …) the residual is near-linear, so
/// the iterate converges in a couple of steps.
pub(super) fn project_uv(surface: &dyn Surface, p: &GpPnt) -> (f64, f64) {
    let (u0, u1) = surface.u_range();
    let (v0, v1) = surface.v_range();
    let su = if u0.is_finite() && u1.is_finite() { 0.5 * (u0 + u1) } else { 0.0 };
    let sv = if v0.is_finite() && v1.is_finite() { 0.5 * (v0 + v1) } else { 0.0 };
    let (u, v, _) = refine_point_on_surface(surface, *p, su, sv, 12);
    (u, v)
}
