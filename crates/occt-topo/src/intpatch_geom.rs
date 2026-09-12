//! UV sampling, projection, and analytic-surface extraction for IntPatch.
//! Source: helpers previously in `intpatch.rs`.

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt, GpVec};
use occt_geom::Surface;

use crate::brep_surface::sphere_center;

// ---------------------------------------------------------------------------
// Parameter helpers
// ---------------------------------------------------------------------------

/// Finite sampling bounds for a surface; unbounded ranges clamp to ±1
/// (matching `brep_surface::sample_bounds`).
pub fn sample_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// Newton refinement of `(u, v)` minimizing `|s(u, v) − p|`.
///
/// Uses central finite differences of `d0` for the surface partials (many
/// ported analytic surfaces only implement `d0` and return zero `d1` vectors).
/// Returns `(u, v, s(u, v))`.
pub fn refine_point_on_surface(s: &dyn Surface, p: GpPnt, u0: f64, v0: f64, iters: usize) -> (f64, f64, GpPnt) {
    let (umin, umax) = s.u_range();
    let (vmin, vmax) = s.v_range();
    let clamp_u = |x: f64| if umin.is_finite() && umax.is_finite() { x.clamp(umin, umax) } else { x };
    let clamp_v = |x: f64| if vmin.is_finite() && vmax.is_finite() { x.clamp(vmin, vmax) } else { x };
    let (mut u, mut v) = (clamp_u(u0), clamp_v(v0));
    let h = 1e-6;
    for _ in 0..iters {
        let p0 = s.d0(u, v);
        let r = GpVec::from_pnts(&p0, &p);
        let pu = GpVec::from_pnts(&p0, &s.d0(u + h, v)).divided(h);
        let pv = GpVec::from_pnts(&p0, &s.d0(u, v + h)).divided(h);
        let (g11, g12, g22) = (pu.dot(&pu), pu.dot(&pv), pv.dot(&pv));
        let (b1, b2) = (pu.dot(&r), pv.dot(&r));
        let det = g11 * g22 - g12 * g12;
        if det.abs() < 1e-30 {
            break;
        }
        let du = (b1 * g22 - b2 * g12) / det;
        let dv = (b2 * g11 - b1 * g12) / det;
        u = clamp_u(u + du);
        v = clamp_v(v + dv);
        if du * du + dv * dv < 1e-24 {
            break;
        }
    }
    (u, v, s.d0(u, v))
}

/// Coarse grid search + refinement for the `(u, v)` parameters of the surface
/// point nearest `p`.
pub fn project_params(s: &dyn Surface, p: &GpPnt) -> (f64, f64) {
    let (u0, u1, v0, v1) = sample_bounds(s);
    let (nu, nv) = (24usize, 24usize);
    let mut bu = u0;
    let mut bv = v0;
    let mut bd = f64::INFINITY;
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d = p.square_distance(&s.d0(u, v));
            if d < bd {
                bd = d;
                bu = u;
                bv = v;
            }
        }
    }
    let (u, v, _) = refine_point_on_surface(s, *p, bu, bv, 8);
    (u, v)
}

/// Approximate minimum distance from `p` to the surface `s`.
pub fn distance_to_surface(p: &GpPnt, s: &dyn Surface) -> f64 {
    let (u, v) = project_params(s, p);
    s.d0(u, v).distance(p)
}

/// The surface params of `p` when it lies within `tol` of the surface.
pub fn point_on_surface(s: &dyn Surface, p: &GpPnt, tol: f64) -> Option<(f64, f64)> {
    let (u, v) = project_params(s, p);
    if s.d0(u, v).distance(p) <= tol {
        Some((u, v))
    } else {
        None
    }
}

/// Extract a `GpPln` from a plane-like surface by sampling its geometry.
pub fn plane_from_surface(s: &dyn Surface) -> Option<GpPln> {
    if !crate::face_face::is_plane_like(s) {
        return None;
    }
    let (p, n) = crate::face_face::plane_geometry(s);
    let d = GpDir::from_vec(&n).ok()?;
    let z = GpDir::new(0.0, 0.0, 1.0).ok()?;
    let x = if d.is_normal(&z) { z } else { GpDir::new(1.0, 0.0, 0.0).ok()? };
    Some(GpPln::new(GpAx3::new(p, d, &x).ok()?))
}

/// Extract sphere center + radius from a surface by sampling.
pub fn sphere_params(s: &dyn Surface) -> Option<(GpPnt, f64)> {
    let c = sphere_center(s)?;
    let (u0, _, v0, _) = sample_bounds(s);
    let r = s.d0(u0, v0).distance(&c);
    Some((c, r))
}
