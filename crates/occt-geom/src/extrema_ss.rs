//! Surface-surface extrema. Port of `Extrema_ExtSS`, `Extrema_ExtElSS`,
//! `Extrema_FuncExtSS`, `Extrema_FuncPSDist`, `Extrema_FuncPSNorm`
//! (TKGeomBase).
//!
//! Planes and spheres are classified by geometric invariants (see
//! `extrema_surf`) and solved analytically (`Extrema_ExtElSS`: plane-plane,
//! and the closed-form sphere-sphere / plane-sphere distances that OCCT's
//! `ExtElSS` leaves as `Standard_NotImplemented`). Everything else goes
//! through the general Newton path on the 4 orthogonality equations
//! `F = ((S1-S2)·Su1, (S1-S2)·Sv1, (S1-S2)·Su2, (S1-S2)·Sv2) = 0`.
//!
//! The `Surface` trait exposes only `d0`/`d1`, so the 4×4 Jacobian is computed
//! by central finite differences of `d1` (`eps ≈ 1e-6` relative to the
//! parameter range). // ponytail: numeric Jacobian, Surface trait lacks d2

use std::cmp::Ordering;

use occt_core::gp::{GpPln, GpPnt, GpSphere, GpVec};

use crate::extrema::ExtremaPair;
use crate::extrema_surf::{classify_plane, classify_sphere, plane_parameters, sphere_parameters, surf_bound_u, surf_bound_v};
use crate::surface::Surface;

/// Surface-surface pair: `p1` on `s1` (params `u1`,`v1`), `p2` on `s2`
/// (params `u2`,`v2`).
fn ss_pair(p1: GpPnt, u1: f64, v1: f64, p2: GpPnt, u2: f64, v2: f64, d: f64) -> ExtremaPair {
    ExtremaPair {
        p1,
        p2,
        distance: d,
        u1,
        v1: Some(v1),
        u2,
        v2: Some(v2),
    }
}

/// Unit vector perpendicular to `v`.
fn perpendicular(v: &GpVec) -> GpVec {
    let a = if v.x().abs() < v.y().abs() {
        GpVec::new(1.0, 0.0, 0.0)
    } else {
        GpVec::new(0.0, 1.0, 0.0)
    };
    v.crossed(&a).normalized()
}

/// Orthogonal projection of `q` onto the plane.
fn project_on_plane(pl: &GpPln, q: &GpPnt) -> GpPnt {
    let n = GpVec::from_xyz(pl.axis().direction().xyz());
    let v = GpVec::from_pnts(&pl.location(), q).dot(&n);
    q.translated_vec(&n.multiplied_scalar(-v))
}

// ---------------------------------------------------------------------------
// Analytic solvers.
// ---------------------------------------------------------------------------

/// Two spheres: exact min and max. The min is `max(0, |r1-r2|-d, d-r1-r2)`
/// (0 when the surfaces intersect), the max is `d+r1+r2`; the closest /
/// farthest point pairs lie on the line through the centers (or on an
/// intersection circle when they intersect).
///
/// **UNPORTED / 本地扩展 (audit A15)**: OCCT's
/// `Extrema_ExtElSS::Perform(const gp_Sphere&, const gp_Sphere&)`
/// (`Extrema_ExtElSS.cxx:77-83`) sets `myDone`/`myNbExt` and then throws
/// `Standard_NotImplemented` — there is **no** closed-form sphere–sphere extrema
/// in OCCT, so these values cannot be compared against OCCT behaviour.
pub fn sphere_sphere_extrema(sp1: &GpSphere, sp2: &GpSphere) -> Vec<ExtremaPair> {
    let c1 = sp1.location();
    let c2 = sp2.location();
    let r1 = sp1.radius();
    let r2 = sp2.radius();
    let d = c1.distance(&c2);
    let min = (0.0f64).max((r1 - r2).abs() - d).max(d - r1 - r2);
    let max = d + r1 + r2;
    let dir = if d > 1e-12 {
        GpVec::from_pnts(&c1, &c2).divided(d)
    } else {
        GpVec::new(1.0, 0.0, 0.0)
    };
    let (n1, n2) = if d >= r1 + r2 {
        // disjoint: facing points on the line through the centers
        (c1.translated_vec(&dir.multiplied_scalar(r1)), c2.translated_vec(&dir.multiplied_scalar(-r2)))
    } else if d <= (r1 - r2).abs() {
        // one strictly inside the other: closest points on the same ray
        if r2 >= r1 {
            (c1.translated_vec(&dir.multiplied_scalar(r1)), c2.translated_vec(&dir.multiplied_scalar(r2)))
        } else {
            (c1.translated_vec(&dir.multiplied_scalar(-r1)), c2.translated_vec(&dir.multiplied_scalar(-r2)))
        }
    } else if d > 1e-12 {
        // intersecting/tangent: a common point on the intersection circle
        let h = ((r1 * r1 - r2 * r2 + d * d) / (2.0 * d)).clamp(-r1, r1);
        let rho = (r1 * r1 - h * h).max(0.0).sqrt();
        let perp = perpendicular(&dir);
        let p = c1
            .translated_vec(&dir.multiplied_scalar(h))
            .translated_vec(&perp.multiplied_scalar(rho));
        (p, p)
    } else {
        // concentric equal spheres: coincide everywhere
        let p = c1.translated_vec(&dir.multiplied_scalar(r1));
        (p, p)
    };
    let f1 = c1.translated_vec(&dir.multiplied_scalar(-r1));
    let f2 = c2.translated_vec(&dir.multiplied_scalar(r2));
    let (u1, v1) = sphere_parameters(sp1, &n1);
    let (u2, v2) = sphere_parameters(sp2, &n2);
    let (w1, x1) = sphere_parameters(sp1, &f1);
    let (w2, x2) = sphere_parameters(sp2, &f2);
    vec![
        ss_pair(n1, u1, v1, n2, u2, v2, min),
        ss_pair(f1, w1, x1, f2, w2, x2, max),
    ]
}

/// Plane vs sphere exact distances (returns `(plane point, sphere point, d)`).
fn plane_sphere_pairs(pl: &GpPln, sp: &GpSphere) -> Vec<(GpPnt, GpPnt, f64)> {
    let n = GpVec::from_xyz(pl.axis().direction().xyz());
    let c = sp.location();
    let r = sp.radius();
    let delta = GpVec::from_pnts(&pl.location(), &c).dot(&n);
    let ad = delta.abs();
    let min = if ad >= r { ad - r } else { 0.0 };
    let max = ad + r;
    let s = if delta >= 0.0 { -1.0 } else { 1.0 };
    let near = c.translated_vec(&n.multiplied_scalar(s * r));
    let far = c.translated_vec(&n.multiplied_scalar(-s * r));
    vec![
        (project_on_plane(pl, &near), near, min),
        (project_on_plane(pl, &far), far, max),
    ]
}

/// Plane vs sphere (s1 = plane, s2 = sphere).
///
/// **UNPORTED / 本地扩展 (audit A15)**: `Extrema_ExtElSS::Perform(const gp_Pln&,
/// const gp_Sphere&)` (`Extrema_ExtElSS.cxx:62-69`) sets `myDone`/`myNbExt` and
/// throws `Standard_NotImplemented`; OCCT has no closed-form plane–sphere
/// extrema, so this is a port-local extension.
pub fn plane_sphere_extrema(pl: &GpPln, sp: &GpSphere) -> Vec<ExtremaPair> {
    plane_sphere_pairs(pl, sp)
        .into_iter()
        .map(|(pp, ps, d)| {
            let (up, vp) = plane_parameters(pl, &pp);
            let (us, vs) = sphere_parameters(sp, &ps);
            ss_pair(pp, up, vp, ps, us, vs, d)
        })
        .collect()
}

/// Sphere vs plane (s1 = sphere, s2 = plane).
pub fn sphere_plane_extrema(sp: &GpSphere, pl: &GpPln) -> Vec<ExtremaPair> {
    plane_sphere_pairs(pl, sp)
        .into_iter()
        .map(|(pp, ps, d)| {
            let (up, vp) = plane_parameters(pl, &pp);
            let (us, vs) = sphere_parameters(sp, &ps);
            ss_pair(ps, us, vs, pp, up, vp, d)
        })
        .collect()
}

/// Plane vs plane: parallel planes have a constant distance, non-parallel
/// planes meet (distance 0).
fn plane_plane_extrema(p1: &GpPln, p2: &GpPln) -> Vec<ExtremaPair> {
    let n1 = GpVec::from_xyz(p1.axis().direction().xyz());
    let n2 = GpVec::from_xyz(p2.axis().direction().xyz());
    if n1.cross_magnitude(&n2) < 1e-9 {
        let d0 = GpVec::from_pnts(&p2.location(), &p1.location()).dot(&n1).abs();
        let p = p1.location();
        let q = project_on_plane(p2, &p);
        let (u1, v1) = plane_parameters(p1, &p);
        let (u2, v2) = plane_parameters(p2, &q);
        vec![ss_pair(p, u1, v1, q, u2, v2, d0)]
    } else {
        // Find a point on the intersection line: parametrize p1 and hit p2.
        let x1 = GpVec::from_xyz(p1.position().x_direction().xyz());
        let y1 = GpVec::from_xyz(p1.position().y_direction().xyz());
        let d = GpVec::from_pnts(&p2.location(), &p1.location()).dot(&n2);
        let a = x1.dot(&n2);
        let b = y1.dot(&n2);
        let (u, v) = if b.abs() > 1e-12 {
            (0.0, d / b)
        } else if a.abs() > 1e-12 {
            (d / a, 0.0)
        } else {
            (0.0, 0.0)
        };
        let p = p1
            .location()
            .translated_vec(&x1.multiplied_scalar(u))
            .translated_vec(&y1.multiplied_scalar(v));
        let (u2, v2) = plane_parameters(p2, &p);
        vec![ss_pair(p, u, v, p, u2, v2, 0.0)]
    }
}

// ---------------------------------------------------------------------------
// General Newton path (port of `Extrema_GenExtSS` + `Extrema_FuncExtSS`,
// numeric 4×4 Jacobian).
// ---------------------------------------------------------------------------

/// F(u1,v1,u2,v2) = ((S1-S2)·Su1, (S1-S2)·Sv1, (S1-S2)·Su2, (S1-S2)·Sv2).
fn ss_f(s1: &dyn Surface, s2: &dyn Surface, u1: f64, v1: f64, u2: f64, v2: f64) -> Option<[f64; 4]> {
    let (p1, su1, sv1) = s1.d1(u1, v1);
    let (p2, su2, sv2) = s2.d1(u2, v2);
    if !(p1.x().is_finite() && p2.x().is_finite()) {
        return None;
    }
    let w = GpVec::from_pnts(&p2, &p1); // S1 - S2
    Some([w.dot(&su1), w.dot(&sv1), w.dot(&su2), w.dot(&sv2)])
}

/// Solve a 4×4 linear system via Gaussian elimination with partial pivoting.
fn solve4(a: &[[f64; 4]; 4], b: &[f64; 4]) -> Option<[f64; 4]> {
    let mut m = *a;
    let mut rhs = *b;
    for col in 0..4 {
        let mut piv = col;
        for row in (col + 1)..4 {
            if m[row][col].abs() > m[piv][col].abs() {
                piv = row;
            }
        }
        if m[piv][col].abs() < 1e-300 {
            return None;
        }
        m.swap(col, piv);
        rhs.swap(col, piv);
        let dd = m[col][col];
        for row in (col + 1)..4 {
            let f = m[row][col] / dd;
            for c in col..4 {
                m[row][c] -= f * m[col][c];
            }
            rhs[row] -= f * rhs[col];
        }
    }
    let mut x = [0.0; 4];
    for i in (0..4).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..4 {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

/// Newton on F = 0 (4×4) from a seed, clamped to the parameter box.
fn solve_surface_surface(
    s1: &dyn Surface,
    s2: &dyn Surface,
    x0: [f64; 4],
    ua1: f64,
    ub1: f64,
    va1: f64,
    vb1: f64,
    ua2: f64,
    ub2: f64,
    va2: f64,
    vb2: f64,
) -> Option<[f64; 4]> {
    let cl = |x: f64, a: f64, b: f64| x.clamp(a.min(b), a.max(b));
    let lo = [ua1.min(ub1), va1.min(vb1), ua2.min(ub2), va2.min(vb2)];
    let hi = [ua1.max(ub1), va1.max(vb1), ua2.max(ub2), va2.max(vb2)];
    let mut x = [
        cl(x0[0], lo[0], hi[0]),
        cl(x0[1], lo[1], hi[1]),
        cl(x0[2], lo[2], hi[2]),
        cl(x0[3], lo[3], hi[3]),
    ];
    let hs = [
        ((ub1 - ua1).abs() * 1e-6).max(1e-9),
        ((vb1 - va1).abs() * 1e-6).max(1e-9),
        ((ub2 - ua2).abs() * 1e-6).max(1e-9),
        ((vb2 - va2).abs() * 1e-6).max(1e-9),
    ];
    let base = ss_f(s1, s2, x[0], x[1], x[2], x[3])?.iter().map(|x| x.abs()).sum::<f64>();
    let mut converged = false;
    for _ in 0..40 {
        let f = ss_f(s1, s2, x[0], x[1], x[2], x[3])?;
        let mut j = [[0.0; 4]; 4];
        for k in 0..4 {
            let h = hs[k];
            let mut xp = x;
            xp[k] += h;
            let mut xm = x;
            xm[k] -= h;
            let fp = ss_f(s1, s2, xp[0], xp[1], xp[2], xp[3])?;
            let fm = ss_f(s1, s2, xm[0], xm[1], xm[2], xm[3])?;
            for r in 0..4 {
                j[r][k] = (fp[r] - fm[r]) / (2.0 * h);
            }
        }
        let d = match solve4(&j, &f) {
            Some(d) => d,
            None => break,
        };
        let xn = [
            cl(x[0] - d[0], lo[0], hi[0]),
            cl(x[1] - d[1], lo[1], hi[1]),
            cl(x[2] - d[2], lo[2], hi[2]),
            cl(x[3] - d[3], lo[3], hi[3]),
        ];
        let step = (xn[0] - x[0]).abs()
            + (xn[1] - x[1]).abs()
            + (xn[2] - x[2]).abs()
            + (xn[3] - x[3]).abs();
        if (xn[0] - x[0]).abs() < 1e-10 * (1.0 + x[0].abs())
            && (xn[1] - x[1]).abs() < 1e-10 * (1.0 + x[1].abs())
            && (xn[2] - x[2]).abs() < 1e-10 * (1.0 + x[2].abs())
            && (xn[3] - x[3]).abs() < 1e-10 * (1.0 + x[3].abs())
        {
            x = xn;
            converged = true;
            break;
        }
        if step > 1e6 * (1.0 + x.iter().map(|v| v.abs()).sum::<f64>()) {
            return None;
        }
        x = xn;
    }
    // A converged step locates a stationary point even when the residual stays
    // non-zero at a domain boundary (one-sided derivatives); accept it.
    if converged {
        return Some(x);
    }
    let resid = ss_f(s1, s2, x[0], x[1], x[2], x[3])?.iter().map(|x| x.abs()).sum::<f64>();
    if resid < base * 1e-2 + 1e-6 {
        Some(x)
    } else {
        None
    }
}

/// All local extrema of the surface-surface distance via grid seeding +
/// Newton, deduplicated and sorted.
pub(crate) fn surface_surface_newton_all(s1: &dyn Surface, s2: &dyn Surface) -> Vec<ExtremaPair> {
    let (u10, u11) = surf_bound_u(s1);
    let (v10, v11) = surf_bound_v(s1);
    let (u20, u21) = surf_bound_u(s2);
    let (v20, v21) = surf_bound_v(s2);
    let (n1, n2, n3, n4) = (8, 8, 8, 8);
    let (a1, b1, c1, d1) = (n1 + 1, n2 + 1, n3 + 1, n4 + 1);
    let len = a1 * b1 * c1 * d1;
    let mut d2 = vec![0.0; len];
    let idx = |i: usize, j: usize, k: usize, l: usize| ((i * b1 + j) * c1 + k) * d1 + l;
    let at = |i: usize, j: usize, k: usize, l: usize| {
        (
            u10 + (u11 - u10) * i as f64 / n1 as f64,
            v10 + (v11 - v10) * j as f64 / n2 as f64,
            u20 + (u21 - u20) * k as f64 / n3 as f64,
            v20 + (v21 - v20) * l as f64 / n4 as f64,
        )
    };
    for i in 0..=n1 {
        for j in 0..=n2 {
            for k in 0..=n3 {
                for l in 0..=n4 {
                    let (u1, v1, u2, v2) = at(i, j, k, l);
                    let d = s1.d0(u1, v1).square_distance(&s2.d0(u2, v2));
                    d2[idx(i, j, k, l)] = if d.is_finite() { d } else { f64::INFINITY };
                }
            }
        }
    }
    let mut seeds: Vec<(usize, usize, usize, usize)> = Vec::new();
    for i in 1..n1 {
        for j in 1..n2 {
            for k in 1..n3 {
                for l in 1..n4 {
                    let d = d2[idx(i, j, k, l)];
                    let lo = d2[idx(i - 1, j, k, l)] >= d
                        && d2[idx(i + 1, j, k, l)] >= d
                        && d2[idx(i, j - 1, k, l)] >= d
                        && d2[idx(i, j + 1, k, l)] >= d
                        && d2[idx(i, j, k - 1, l)] >= d
                        && d2[idx(i, j, k + 1, l)] >= d
                        && d2[idx(i, j, k, l - 1)] >= d
                        && d2[idx(i, j, k, l + 1)] >= d;
                    let hi = d2[idx(i - 1, j, k, l)] <= d
                        && d2[idx(i + 1, j, k, l)] <= d
                        && d2[idx(i, j - 1, k, l)] <= d
                        && d2[idx(i, j + 1, k, l)] <= d
                        && d2[idx(i, j, k - 1, l)] <= d
                        && d2[idx(i, j, k + 1, l)] <= d
                        && d2[idx(i, j, k, l - 1)] <= d
                        && d2[idx(i, j, k, l + 1)] <= d;
                    if lo || hi {
                        seeds.push((i, j, k, l));
                    }
                }
            }
        }
    }
    // Global grid min/max (boundary extrema) + the 4D corners.
    let (mut gmin, mut gmax) = ((0usize, 0usize, 0usize, 0usize), (0usize, 0usize, 0usize, 0usize));
    for i in 0..=n1 {
        for j in 0..=n2 {
            for k in 0..=n3 {
                for l in 0..=n4 {
                    let d = d2[idx(i, j, k, l)];
                    if d < d2[idx(gmin.0, gmin.1, gmin.2, gmin.3)] {
                        gmin = (i, j, k, l);
                    }
                    if d > d2[idx(gmax.0, gmax.1, gmax.2, gmax.3)] {
                        gmax = (i, j, k, l);
                    }
                }
            }
        }
    }
    seeds.push(gmin);
    seeds.push(gmax);
    for mask in 0..16u8 {
        let i = if mask & 1 != 0 { n1 } else { 0 };
        let j = if mask & 2 != 0 { n2 } else { 0 };
        let k = if mask & 4 != 0 { n3 } else { 0 };
        let l = if mask & 8 != 0 { n4 } else { 0 };
        seeds.push((i, j, k, l));
    }

    let mut out: Vec<ExtremaPair> = Vec::new();
    for (i, j, k, l) in seeds {
        let (u1, v1, u2, v2) = at(i, j, k, l);
        if let Some(x) = solve_surface_surface(s1, s2, [u1, v1, u2, v2], u10, u11, v10, v11, u20, u21, v20, v21) {
            let p1 = s1.d0(x[0], x[1]);
            let p2 = s2.d0(x[2], x[3]);
            if !(p1.x().is_finite() && p2.x().is_finite()) {
                continue;
            }
            if !out
                .iter()
                .any(|e: &ExtremaPair| e.p1.distance(&p1) < 1e-6 && e.p2.distance(&p2) < 1e-6)
            {
                out.push(ss_pair(p1, x[0], x[1], p2, x[2], x[3], p1.distance(&p2)));
            }
        }
    }
    out.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Ordering::Equal));
    out
}

/// Fallback: coarse grid min + coordinate descent (degenerate pairs).
fn fallback_ss(s1: &dyn Surface, s2: &dyn Surface) -> ExtremaPair {
    let (u10, u11) = surf_bound_u(s1);
    let (v10, v11) = surf_bound_v(s1);
    let (u20, u21) = surf_bound_u(s2);
    let (v20, v21) = surf_bound_v(s2);
    let n = 12;
    let mut best = (0.0, 0.0, 0.0, 0.0, f64::INFINITY);
    for i in 0..=n {
        for j in 0..=n {
            let u1 = u10 + (u11 - u10) * i as f64 / n as f64;
            let v1 = v10 + (v11 - v10) * j as f64 / n as f64;
            let p = s1.d0(u1, v1);
            for k in 0..=n {
                for l in 0..=n {
                    let u2 = u20 + (u21 - u20) * k as f64 / n as f64;
                    let v2 = v20 + (v21 - v20) * l as f64 / n as f64;
                    let d = p.square_distance(&s2.d0(u2, v2));
                    if d < best.4 {
                        best = (u1, v1, u2, v2, d);
                    }
                }
            }
        }
    }
    let p1 = s1.d0(best.0, best.1);
    let p2 = s2.d0(best.2, best.3);
    ss_pair(p1, best.0, best.1, p2, best.2, best.3, p1.distance(&p2))
}

// ---------------------------------------------------------------------------
// Public dispatch.
// ---------------------------------------------------------------------------

/// All local extrema of the surface-surface distance, deduplicated and sorted.
/// Spheres and planes are classified and solved analytically (exact); all
/// other pairs go through the grid + Newton path.
pub fn surface_surface_extrema_all(s1: &dyn Surface, s2: &dyn Surface) -> Vec<ExtremaPair> {
    match (classify_sphere(s1), classify_sphere(s2)) {
        (Some(a), Some(b)) => return sphere_sphere_extrema(&a, &b),
        _ => {}
    }
    if let (Some(a), Some(p)) = (classify_sphere(s1), classify_plane(s2)) {
        return sphere_plane_extrema(&a, &p);
    }
    if let (Some(p), Some(a)) = (classify_plane(s1), classify_sphere(s2)) {
        return plane_sphere_extrema(&p, &a);
    }
    if let (Some(p1), Some(p2)) = (classify_plane(s1), classify_plane(s2)) {
        return plane_plane_extrema(&p1, &p2);
    }
    surface_surface_newton_all(s1, s2)
}

/// Minimum distance between surfaces `s1` and `s2`.
pub fn surface_surface_extrema(s1: &dyn Surface, s2: &dyn Surface) -> ExtremaPair {
    match surface_surface_extrema_all(s1, s2).into_iter().next() {
        Some(e) => e,
        None => fallback_ss(s1, s2),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bspline_surface::{bspline_surface_uniform_knots, fit_surface_grid, GeomBSplineSurface};
    use crate::{GeomPlane, GeomSphere};
    use occt_core::gp::{GpAx3, GpPnt, GpSphere as GpSphereT};

    const PI: f64 = std::f64::consts::PI;

    fn unit_sphere_at(p: GpPnt) -> GeomSphere {
        let mut ax = GpAx3::standard();
        ax.set_location(p);
        GeomSphere::new(GpSphereT::new(ax, 1.0).unwrap())
    }

    fn plane_z(z: f64) -> GeomPlane {
        let mut pl = occt_core::gp::GpPln::new(GpAx3::standard());
        pl.set_location(&GpPnt::new(0.0, 0.0, z));
        GeomPlane::new(pl)
    }

    #[test]
    fn sphere_sphere_extrema_min_and_max() {
        let s1 = unit_sphere_at(GpPnt::new(0.0, 0.0, 0.0));
        let s2 = unit_sphere_at(GpPnt::new(4.0, 0.0, 0.0));
        let all = surface_surface_extrema_all(&s1, &s2);
        assert_eq!(all.len(), 2, "sphere-sphere extrema {all:?}");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        let max = all.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 2.0).abs() < 1e-9, "min {}", min.distance);
        assert!((min.p1.x() - 1.0).abs() < 1e-9 && (min.p2.x() - 3.0).abs() < 1e-9, "closest {min:?}");
        assert!((max.distance - 6.0).abs() < 1e-9, "max {}", max.distance);
        assert!((max.p1.x() + 1.0).abs() < 1e-9 && (max.p2.x() - 5.0).abs() < 1e-9, "farthest {max:?}");
    }

    #[test]
    fn surface_surface_sphere_plane_min() {
        let s1 = unit_sphere_at(GpPnt::new(0.0, 0.0, 0.0));
        let s2 = plane_z(5.0);
        let e = surface_surface_extrema(&s1, &s2);
        assert!((e.distance - 4.0).abs() < 1e-9, "min {}", e.distance);
        assert!((e.p1.z() - 1.0).abs() < 1e-9 && (e.p2.z() - 5.0).abs() < 1e-9, "closest {e:?}");
    }

    #[test]
    fn surface_surface_planes_parallel() {
        let p1 = GeomPlane::new(occt_core::gp::GpPln::new(GpAx3::standard()));
        let p2 = plane_z(3.0);
        let e = surface_surface_extrema(&p1, &p2);
        assert!((e.distance - 3.0).abs() < 1e-9, "parallel planes {}", e.distance);
    }

    #[test]
    fn surface_surface_planes_coincident() {
        let p1 = GeomPlane::new(occt_core::gp::GpPln::new(GpAx3::standard()));
        let p2 = GeomPlane::new(occt_core::gp::GpPln::new(GpAx3::standard()));
        let e = surface_surface_extrema(&p1, &p2);
        assert!(e.distance.abs() < 1e-9, "coincident planes {}", e.distance);
    }

    fn paraboloid_plus_one() -> GeomBSplineSurface {
        let (nu, nv) = (3, 3);
        let points: Vec<Vec<GpPnt>> = (0..nu)
            .map(|i| {
                (0..nv)
                    .map(|j| {
                        let u = i as f64 / (nu - 1) as f64;
                        let v = j as f64 / (nv - 1) as f64;
                        GpPnt::new(u, v, u * u + v * v + 1.0)
                    })
                    .collect()
            })
            .collect();
        fit_surface_grid(&points, 2, 2).unwrap()
    }

    #[test]
    fn surface_surface_newton_paraboloid_plane_min() {
        // S1(u,v) = (u, v, u²+v²+1) (exact degree-2 B-spline), S2 = plane
        // z = 0.5 over [0,1]². Minimum distance 0.5 at (u1,v1)=(0,0).
        let s1 = paraboloid_plus_one();
        let s2 = GeomBSplineSurface::new(
            vec![
                vec![GpPnt::new(0.0, 0.0, 0.5), GpPnt::new(0.0, 1.0, 0.5)],
                vec![GpPnt::new(1.0, 0.0, 0.5), GpPnt::new(1.0, 1.0, 0.5)],
            ],
            bspline_surface_uniform_knots(2, 2, 1, 1).0,
            bspline_surface_uniform_knots(2, 2, 1, 1).1,
            1,
            1,
        )
        .unwrap();
        let all = surface_surface_newton_all(&s1, &s2);
        assert!(!all.is_empty(), "no surface-surface extrema");
        let min = all.iter().min_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap();
        assert!((min.distance - 0.5).abs() < 1e-6, "min {}", min.distance);
        // Orthogonality holds (residual is at the boundary FD error, ~1e-5).
        let f = ss_f(&s1, &s2, min.u1, min.v1.unwrap(), min.u2, min.v2.unwrap()).unwrap();
        assert!(f.iter().map(|x| x.abs()).sum::<f64>() < 1e-3, "orthogonality F={f:?}");
    }

    #[test]
    fn sphere_sphere_direct_analytic() {
        let s1 = GpSphereT::new(GpAx3::standard(), 1.0).unwrap();
        let mut ax = GpAx3::standard();
        ax.set_location(GpPnt::new(4.0, 0.0, 0.0));
        let s2 = GpSphereT::new(ax, 1.0).unwrap();
        let all = sphere_sphere_extrema(&s1, &s2);
        assert_eq!(all.len(), 2);
        assert!((all[0].distance - 2.0).abs() < 1e-9, "min {}", all[0].distance);
        assert!((all[1].distance - 6.0).abs() < 1e-9, "max {}", all[1].distance);
    }
}
