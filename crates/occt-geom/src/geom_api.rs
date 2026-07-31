//! OCCT port: GeomAPI projection, intersection and extremum helpers.
//! Source: `GeomAPI_ProjectPointOnCurve.hxx`, `GeomAPI_ProjectPointOnSurf.hxx`,
//! `GeomAPI_IntCS.hxx`, `GeomAPI_ExtremaCurveCurve.hxx`.
//!
//! ponytail: sampler-based approximations throughout — not root-exact like
//! OCCT's Extrema/IntCurve algorithms. Adequate for geometry tooling; replace
//! with subdivision/Newton solvers if exact roots on arbitrary curves are needed.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use crate::curve::Curve;
use crate::surface::Surface;

/// Orthogonal projection of a point onto a curve (nearest solution).
pub struct PointOnCurve {
    pub parameter: f64,
    pub point: GpPnt,
    pub distance: f64,
}

/// Orthogonal projection of a point onto a surface (nearest solution).
pub struct PointOnSurface {
    pub u: f64,
    pub v: f64,
    pub point: GpPnt,
    pub distance: f64,
}

/// Nearest point on `c` to `p` (coarse scan + golden-section refinement).
/// Unbounded curves (lines) are handled by probing an expanding window around
/// `u = 0` until the minimum is interior.
pub fn project_point_on_curve(c: &dyn Curve, p: &GpPnt, _tol: f64) -> Option<PointOnCurve> {
    let mk = |u: f64, q: GpPnt| PointOnCurve { parameter: u, point: q, distance: q.distance(p) };
    let a = c.first_parameter();
    let b = c.last_parameter();
    if !a.is_finite() || !b.is_finite() {
        let base = c.d0(0.0);
        let d = base.distance(p);
        let mut w = (d + 1.0).max(1.0) * 8.0;
        for _ in 0..4 {
            if let Some((u, q)) = refine_closest_curve(c, p, -w, w) {
                if (u + w).abs() > 1e-9 && (u - w).abs() > 1e-9 {
                    return Some(mk(u, q));
                }
            }
            w *= 8.0;
        }
        return refine_closest_curve(c, p, -w, w).map(|(u, q)| mk(u, q));
    }
    refine_closest_curve(c, p, a, b).map(|(u, q)| mk(u, q))
}

/// Nearest (u, v) parameters of `p` on `s` (grid search + hill-climb).
/// Unbounded parameter directions use an expanding window.
pub fn project_point_on_surface(s: &dyn Surface, p: &GpPnt, _tol: f64) -> Option<PointOnSurface> {
    let (u, v) = surface_closest_params(s, p);
    let point = s.d0(u, v);
    Some(PointOnSurface { u, v, point, distance: point.distance(p) })
}

/// Intersection points of curve `c` with surface `s`.
///
/// Samples the curve, detects parameter intervals where the distance to the
/// surface dips through a near-zero local minimum, refines each minimum by
/// golden-section, and returns the 3D intersection points.
pub fn curve_surface_intersections(c: &dyn Curve, s: &dyn Surface, tol: f64, samples: usize) -> Vec<GpPnt> {
    let tol = tol.max(1e-12);
    let bbox = surface_bbox(s, 8, 8);
    let (lo, hi) = curve_window_for_bbox(c, &bbox);
    let n = samples.max(3);
    let step = (hi - lo) / n as f64;

    let mut ds: Vec<(f64, f64)> = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let u = lo + step * i as f64;
        ds.push((u, dist_curve_surface(c, s, u)));
    }

    let mut found: Vec<GpPnt> = Vec::new();
    for i in 1..ds.len() - 1 {
        // A near-zero local minimum of the unsigned distance marks a crossing.
        if ds[i].1 <= ds[i - 1].1 && ds[i].1 <= ds[i + 1].1 {
            let (um, dm) = minimize_1d(&|u| dist_curve_surface(c, s, u), ds[i - 1].0, ds[i + 1].0);
            if dm <= tol {
                let q = c.d0(um);
                let dup = found.iter().any(|p| p.distance(&q) <= tol.max(1e-6));
                if !dup {
                    found.push(q);
                }
            }
        }
    }
    found
}

/// Approximate intersection points of two 3D curves.
///
/// ponytail: approximate sampler-based intersection — not root-exact. Samples
/// both curves on coarse grids, refines candidate pairs by alternating 1-D
/// minimization, and dedupes points within `tol`.
pub fn curve_curve_intersections(c1: &dyn Curve, c2: &dyn Curve, tol: f64) -> Vec<GpPnt> {
    let tol = tol.max(1e-12);
    let na = 256;
    let nb = 256;
    let sa = sample_curve(c1, c2, na);
    let sb = sample_curve(c2, c1, nb);
    if sa.is_empty() || sb.is_empty() {
        return Vec::new();
    }
    let ua0 = sa[0].0;
    let ua1 = sa[sa.len() - 1].0;
    let vb0 = sb[0].0;
    let vb1 = sb[sb.len() - 1].0;
    let du = ((ua1 - ua0) / na as f64).max(1e-12);
    let dv = ((vb1 - vb0) / nb as f64).max(1e-12);

    // Coarse filter sized to the sampling resolution so every crossing's
    // nearest grid pair is caught.
    let coarse = (tol * 100.0 + 1e-9).max(du.max(dv) * 4.0);
    let mut found: Vec<(f64, f64)> = Vec::new();
    for (ui, ai) in &sa {
        for (vj, bj) in &sb {
            if ai.square_distance(bj) < coarse * coarse {
                // Refine the candidate pair by alternating 1-D minimization.
                let mut u = *ui;
                let mut v = *vj;
                for _ in 0..6 {
                    let (nu, _) = minimize_1d(
                        &|t| c1.d0(t).square_distance(&c2.d0(v)),
                        (u - du).max(ua0),
                        (u + du).min(ua1),
                    );
                    u = nu;
                    let (nv, _) = minimize_1d(
                        &|t| c1.d0(u).square_distance(&c2.d0(t)),
                        (v - dv).max(vb0),
                        (v + dv).min(vb1),
                    );
                    v = nv;
                }
                if c1.d0(u).distance(&c2.d0(v)) <= tol {
                    found.push((u, v));
                }
            }
        }
    }

    // Dedupe candidates that converge to the same geometric point.
    let mut out: Vec<GpPnt> = Vec::new();
    for (u, v) in found {
        let pa = c1.d0(u);
        let pb = c2.d0(v);
        let mid = GpPnt::new((pa.x() + pb.x()) * 0.5, (pa.y() + pb.y()) * 0.5, (pa.z() + pb.z()) * 0.5);
        let dup = out.iter().any(|p| p.distance(&mid) <= tol.max(1e-6));
        if !dup {
            out.push(mid);
        }
    }
    out
}

/// Minimum distance between two curves (coarse sampling + coordinate descent).
pub fn curve_curve_distance(c1: &dyn Curve, c2: &dyn Curve, _tol: f64) -> f64 {
    let sa = sample_curve(c1, c2, 128);
    let sb = sample_curve(c2, c1, 128);
    if sa.is_empty() || sb.is_empty() {
        return f64::INFINITY;
    }
    let mut bu = sa[0].0;
    let mut bv = sb[0].0;
    let mut best = f64::INFINITY;
    for (u, pu) in &sa {
        for (v, pv) in &sb {
            let d = pu.square_distance(pv);
            if d < best {
                best = d;
                bu = *u;
                bv = *v;
            }
        }
    }
    let du = ((sa[sa.len() - 1].0 - sa[0].0) / sa.len() as f64).max(1e-12);
    let dv = ((sb[sb.len() - 1].0 - sb[0].0) / sb.len() as f64).max(1e-12);
    for _ in 0..6 {
        let (nu, _) = minimize_1d(&|t| c1.d0(t).square_distance(&c2.d0(bv)), bu - du, bu + du);
        bu = nu;
        let (nv, nd) = minimize_1d(&|t| c1.d0(bu).square_distance(&c2.d0(t)), bv - dv, bv + dv);
        bv = nv;
        best = nd;
    }
    best.sqrt()
}

/// Distance from point `p` to the nearest point on curve `c`.
pub fn distance_point_curve(c: &dyn Curve, p: &GpPnt) -> f64 {
    project_point_on_curve(c, p, 1e-7).map_or(f64::INFINITY, |pc| pc.distance)
}

/// Distance from point `p` to the nearest point on surface `s`.
pub fn distance_point_surface(s: &dyn Surface, p: &GpPnt) -> f64 {
    project_point_on_surface(s, p, 1e-7).map_or(f64::INFINITY, |ps| ps.distance)
}

/// Unit tangent vector of `c` at `u` (from `d1`, normalized).
pub fn tangent_at(c: &dyn Curve, u: f64) -> GpVec {
    let (_, v) = c.d1(u);
    v.normalized()
}

/// Approximate bounding box of a curve by sampling `(pmin, pmax)`.
/// Unbounded curves are sampled over a clamped `[-1, 1]` proxy window.
pub fn curve_bbox(c: &dyn Curve, samples: usize) -> (GpPnt, GpPnt) {
    let (a, b) = finite_range(c.first_parameter(), c.last_parameter());
    let n = samples.max(1);
    let mut min = GpPnt::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = GpPnt::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for i in 0..=n {
        let p = c.d0(a + (b - a) * i as f64 / n as f64);
        min.coord.x = min.coord.x.min(p.x());
        min.coord.y = min.coord.y.min(p.y());
        min.coord.z = min.coord.z.min(p.z());
        max.coord.x = max.coord.x.max(p.x());
        max.coord.y = max.coord.y.max(p.y());
        max.coord.z = max.coord.z.max(p.z());
    }
    (min, max)
}

/// Approximate bounding box of a surface on a `nu × nv` grid `(pmin, pmax)`.
/// Unbounded parameter directions use a clamped `[-1, 1]` proxy window.
pub fn surface_bbox(s: &dyn Surface, nu: usize, nv: usize) -> (GpPnt, GpPnt) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let (u0, u1) = finite_range(u0, u1);
    let (v0, v1) = finite_range(v0, v1);
    let nu = nu.max(1);
    let nv = nv.max(1);
    let mut min = GpPnt::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = GpPnt::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for i in 0..=nu {
        for j in 0..=nv {
            let p = s.d0(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64);
            min.coord.x = min.coord.x.min(p.x());
            min.coord.y = min.coord.y.min(p.y());
            min.coord.z = min.coord.z.min(p.z());
            max.coord.x = max.coord.x.max(p.x());
            max.coord.y = max.coord.y.max(p.y());
            max.coord.z = max.coord.z.max(p.z());
        }
    }
    (min, max)
}

/// Approximate p-curve of `c` on `s`: projection of each curve sample onto the
/// surface, sampled at `samples` parameter values.
pub fn pcurve_of_curve_on_surface(c: &dyn Curve, s: &dyn Surface, samples: usize) -> Vec<GpPnt2d> {
    let (a, b) = finite_range(c.first_parameter(), c.last_parameter());
    let n = samples.max(2);
    (0..n)
        .map(|i| {
            let u = a + (b - a) * i as f64 / (n - 1) as f64;
            let (su, sv) = surface_closest_params(s, &c.d0(u));
            GpPnt2d::new(su, sv)
        })
        .collect()
}

// --- internals -------------------------------------------------------------

/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
fn minimize_1d<F: Fn(f64) -> f64>(f: &F, lo: f64, hi: f64) -> (f64, f64) {
    const GOLD: f64 = 0.6180339887498949;
    let mut a = lo;
    let mut b = hi;
    let mut c = b - GOLD * (b - a);
    let mut d = a + GOLD * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    while (b - a) > 1e-12 {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - GOLD * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + GOLD * (b - a);
            fd = f(d);
        }
    }
    let x = 0.5 * (a + b);
    (x, f(x))
}

/// Coarse-to-fine closest point over a finite `[a, b]`.
fn refine_closest_curve(c: &dyn Curve, p: &GpPnt, a: f64, b: f64) -> Option<(f64, GpPnt)> {
    if !(b > a) {
        return None;
    }
    let n = 64;
    let mut best = a;
    let mut best_d = f64::INFINITY;
    for i in 0..=n {
        let u = a + (b - a) * i as f64 / n as f64;
        let d = c.d0(u).square_distance(p);
        if d < best_d {
            best_d = d;
            best = u;
        }
    }
    let span = (b - a) / n as f64;
    let lo = (best - span).max(a);
    let hi = (best + span).min(b);
    let (u, _) = minimize_1d(&|u| c.d0(u).square_distance(p), lo, hi);
    Some((u, c.d0(u)))
}

/// Closest (u, v) on a surface within a finite window: grid search + hill-climb.
fn closest_params_in_window(s: &dyn Surface, p: &GpPnt, u0: f64, u1: f64, v0: f64, v1: f64) -> (f64, f64) {
    if !(u1 > u0) || !(v1 > v0) {
        return (u0, v0);
    }
    let nu = 16;
    let nv = 16;
    let mut best = (u0, v0);
    let mut best_d = f64::INFINITY;
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d = s.d0(u, v).square_distance(p);
            if d < best_d {
                best_d = d;
                best = (u, v);
            }
        }
    }
    // Coarse-to-fine refinement (6 rounds of local hill-climb).
    let (mut u, mut v) = best;
    let (mut hu, mut hv) = ((u1 - u0) / nu as f64, (v1 - v0) / nv as f64);
    for _ in 0..6 {
        hu *= 0.5;
        hv *= 0.5;
        let mut du = u;
        let mut dv = v;
        let mut dd = s.d0(u, v).square_distance(p);
        for &(su, sv) in &[(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0)] {
            let nu2 = (u + su * hu).clamp(u0, u1);
            let nv2 = (v + sv * hv).clamp(v0, v1);
            let d = s.d0(nu2, nv2).square_distance(p);
            if d < dd {
                dd = d;
                du = nu2;
                dv = nv2;
            }
        }
        u = du;
        v = dv;
    }
    // Golden-section polish per axis over the local neighborhood. The hill-climb
    // steps halve 6 times, so its accuracy is `window/16/2^6`; the polish closes
    // the gap to the true minimum (needed for far/unbounded projections).
    for _ in 0..3 {
        let (nu, _) = minimize_1d(&|t| s.d0(t, v).square_distance(p), (u - hu).max(u0), (u + hu).min(u1));
        u = nu;
        let (nv, _) = minimize_1d(&|t| s.d0(u, t).square_distance(p), (v - hv).max(v0), (v + hv).min(v1));
        v = nv;
    }
    (u, v)
}

/// Closest (u, v) of `p` on `s`, expanding the window for unbounded directions
/// until the solution is interior.
fn surface_closest_params(s: &dyn Surface, p: &GpPnt) -> (f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let unb_u = !u0.is_finite() || !u1.is_finite();
    let unb_v = !v0.is_finite() || !v1.is_finite();
    let (mut lo_u, mut hi_u) = if unb_u { (-1.0, 1.0) } else { (u0, u1) };
    let (mut lo_v, mut hi_v) = if unb_v { (-1.0, 1.0) } else { (v0, v1) };
    let (mut u, mut v) = closest_params_in_window(s, p, lo_u, hi_u, lo_v, hi_v);
    if unb_u || unb_v {
        let mut w = 1.0;
        for _ in 0..8 {
            let on_u = unb_u && ((u - lo_u).abs() < 1e-9 || (u - hi_u).abs() < 1e-9);
            let on_v = unb_v && ((v - lo_v).abs() < 1e-9 || (v - hi_v).abs() < 1e-9);
            if !on_u && !on_v {
                break;
            }
            w *= 8.0;
            if unb_u {
                lo_u = -w;
                hi_u = w;
            }
            if unb_v {
                lo_v = -w;
                hi_v = w;
            }
            (u, v) = closest_params_in_window(s, p, lo_u, hi_u, lo_v, hi_v);
        }
    }
    (u, v)
}

/// Distance from `c(u)` to the nearest point on surface `s`.
fn dist_curve_surface(c: &dyn Curve, s: &dyn Surface, u: f64) -> f64 {
    let p = c.d0(u);
    let (su, sv) = surface_closest_params(s, &p);
    p.distance(&s.d0(su, sv))
}

/// Sample window for `c` that covers the given bounding box. Unbounded curves
/// get a window proportional to the bbox center distance plus size.
fn curve_window_for_bbox(c: &dyn Curve, bbox: &(GpPnt, GpPnt)) -> (f64, f64) {
    let a = c.first_parameter();
    let b = c.last_parameter();
    if a.is_finite() && b.is_finite() {
        return (a, b);
    }
    let (pmin, pmax) = bbox;
    let center = GpPnt::new(0.5 * (pmin.x() + pmax.x()), 0.5 * (pmin.y() + pmax.y()), 0.5 * (pmin.z() + pmax.z()));
    let half_diag = 0.5 * pmin.distance(pmax);
    let d0 = c.d0(0.0).distance(&center);
    let span = (d0 + half_diag + 1.0).max(1.0) * 4.0;
    (-span, span)
}

/// Finite, sane sampling bounds; unbounded ranges clamp to `[-1, 1]`.
fn finite_range(a: f64, b: f64) -> (f64, f64) {
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-1.0, 1.0)
    }
}

/// Sample a curve over a finite parameter window. Unbounded curves get a
/// window that comfortably covers the other curve's bounding box.
fn sample_curve(c: &dyn Curve, other: &dyn Curve, n: usize) -> Vec<(f64, GpPnt)> {
    let a = c.first_parameter();
    let b = c.last_parameter();
    let (lo, hi) = if a.is_finite() && b.is_finite() {
        (a, b)
    } else {
        let (pmin, pmax) = curve_bbox(other, 32);
        let size = pmax.distance(&pmin).max(1e-6);
        let speed = c.d1(0.0).1.magnitude().max(1e-30);
        let span = (4.0 * size / speed).max(1.0);
        (-span, span)
    };
    (0..=n)
        .map(|i| {
            let u = lo + (hi - lo) * i as f64 / n as f64;
            (u, c.d0(u))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circle::GeomCircle;
    use crate::line::GeomLine;
    use crate::plane::GeomPlane;
    use crate::sphere::GeomSphere;
    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpDir, GpPln, GpSphere};

    #[test]
    fn project_point_onto_circle() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let pc = project_point_on_curve(&circle as &dyn Curve, &GpPnt::new(0.0, 2.0, 0.0), 1e-6).unwrap();
        assert!((pc.distance - 1.0).abs() < 1e-6, "distance={}", pc.distance);
        assert!((pc.parameter - std::f64::consts::FRAC_PI_2).abs() < 1e-4, "parameter={}", pc.parameter);
        assert!(pc.point.distance(&GpPnt::new(0.0, 1.0, 0.0)) < 1e-6);
    }

    #[test]
    fn project_point_onto_plane() {
        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let ps = project_point_on_surface(&plane as &dyn Surface, &GpPnt::new(0.5, 0.5, 1.0), 1e-6).unwrap();
        assert!((ps.distance - 1.0).abs() < 1e-6, "distance={}", ps.distance);
        assert!(ps.point.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-6);
        // Far point on the unbounded plane forces window expansion.
        let ps2 = project_point_on_surface(&plane as &dyn Surface, &GpPnt::new(10.0, 10.0, 1.0), 1e-6).unwrap();
        assert!((ps2.distance - 1.0).abs() < 1e-6, "distance={}", ps2.distance);
        assert!(ps2.point.distance(&GpPnt::new(10.0, 10.0, 0.0)) < 1e-6);
    }

    #[test]
    fn line_sphere_intersections() {
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.0, 0.0, 0.5), GpDir::new(0.0, 1.0, 0.0).unwrap());
        let sphere = GeomSphere::new(GpSphere::new(GpAx3::standard(), 1.0).unwrap());
        let hits = curve_surface_intersections(&line as &dyn Curve, &sphere as &dyn Surface, 1e-4, 200);
        assert_eq!(hits.len(), 2, "hits={hits:?}");
        for p in &hits {
            assert!((p.coord.modulus() - 1.0).abs() < 1e-3, "|p|={}", p.coord.modulus());
        }
    }

    #[test]
    fn crossing_lines_intersection() {
        let l1 = GeomLine::from_pnt_dir(GpPnt::zero(), GpDir::new(1.0, 1.0, 0.0).unwrap());
        let l2 = GeomLine::from_pnt_dir(GpPnt::new(2.0, 0.0, 0.0), GpDir::new(-1.0, 1.0, 0.0).unwrap());
        let hits = curve_curve_intersections(&l1 as &dyn Curve, &l2 as &dyn Curve, 1e-6);
        assert_eq!(hits.len(), 1, "hits={hits:?}");
        assert!(hits[0].distance(&GpPnt::new(1.0, 1.0, 0.0)) < 1e-4, "p={:?}", hits[0]);
    }

    #[test]
    fn parallel_lines_distance() {
        let l1 = GeomLine::from_pnt_dir(GpPnt::zero(), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let l2 = GeomLine::from_pnt_dir(GpPnt::new(0.0, 3.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let d = curve_curve_distance(&l1 as &dyn Curve, &l2 as &dyn Curve, 1e-6);
        assert!((d - 3.0).abs() < 1e-4, "d={d}");
    }

    #[test]
    fn distance_and_tangent() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let d = distance_point_curve(&circle as &dyn Curve, &GpPnt::new(3.0, 0.0, 0.0));
        assert!((d - 2.0).abs() < 1e-6, "d={d}");
        let t = tangent_at(&circle as &dyn Curve, 0.0);
        assert!(t.subtracted(&GpVec::new(0.0, 1.0, 0.0)).magnitude() < 1e-6, "tangent={t:?}");
    }

    #[test]
    fn bbox_and_pcurve() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let (pmin, pmax) = curve_bbox(&circle as &dyn Curve, 64);
        assert!((pmin.x() + 1.0).abs() < 1e-6 && (pmax.x() - 1.0).abs() < 1e-6);

        let plane = GeomPlane::new(GpPln::new(GpAx3::standard()));
        let line = GeomLine::from_pnt_dir(GpPnt::new(0.5, 0.5, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let pc = pcurve_of_curve_on_surface(&line as &dyn Curve, &plane as &dyn Surface, 4);
        assert_eq!(pc.len(), 4);
        // Unbounded line sampled over the [-1, 1] proxy window: first/last
        // samples land at (0.5, 0.5, 0) ± (1, 0, 0) → (u, v) = (x, y).
        assert!((pc[0].x() - -0.5).abs() < 1e-6 && (pc[0].y() - 0.5).abs() < 1e-6);
        assert!((pc[3].x() - 1.5).abs() < 1e-6 && (pc[3].y() - 0.5).abs() < 1e-6);
    }
}
