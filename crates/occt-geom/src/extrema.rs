//! Extrema between geometric objects — point/curve/surface minimum and
//! maximum distances.
//!
//! Port of `Extrema_ExtPC`, `Extrema_ExtCC`, `Extrema_ExtPS`,
//! `Extrema_ExtCS`, `Extrema_ExtSS` (multidimensional extrema via
//! discretization + local refinement). Source: `Extrema` (TKGeomBase).

use occt_core::gp::{GpPnt, GpVec};

use crate::{Curve, Surface};

/// A solved extremum between two objects.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtremaPair {
    pub p1: GpPnt,
    pub p2: GpPnt,
    pub distance: f64,
    /// Parameters on each object (curve u, or surface u/v).
    pub u1: f64,
    pub v1: Option<f64>,
    pub u2: f64,
    pub v2: Option<f64>,
}

/// Refine a 1-D closest-point parameter by golden-section / parabolic search
/// around the initial guess. Handles unbounded curves via window expansion.
pub fn refine_curve_point(c: &dyn Curve, p: &GpPnt, a: f64, b: f64) -> (f64, GpPnt) {
    let f = |u: f64| {
        let q = c.d0(u);
        (q.x() - p.x()).powi(2) + (q.y() - p.y()).powi(2) + (q.z() - p.z()).powi(2)
    };
    let (mut lo, mut hi) = if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        // Unbounded: probe expanding windows until a strict minimum is
        // interior (both endpoints larger than the center).
        let mut w = 1.0;
        let mut lo = -w;
        let mut hi = w;
        let mut fm = f(0.0);
        for _ in 0..8 {
            let fl = f(lo);
            let fh = f(hi);
            if fl >= fm && fh >= fm {
                break; // minimum interior to [lo, hi]
            }
            w *= 8.0;
            lo = -w;
            hi = w;
            fm = f(0.0);
        }
        (lo, hi)
    };
    // Golden-section refinement (golden ratio conjugate 0.618...).
    let phi = (5.0f64.sqrt() - 1.0) * 0.5;
    for _ in 0..64 {
        let x1 = hi - phi * (hi - lo);
        let x2 = lo + phi * (hi - lo);
        if f(x1) < f(x2) {
            hi = x2;
        } else {
            lo = x1;
        }
    }
    let u = 0.5 * (lo + hi);
    (u, c.d0(u))
}

/// Minimum distance from point `p` to curve `c` (with the closest point).
pub fn point_curve_extrema(c: &dyn Curve, p: &GpPnt) -> ExtremaPair {
    let (u, q) = refine_curve_point(c, p, c.first_parameter(), c.last_parameter());
    ExtremaPair {
        p1: *p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: None,
        u2: u,
        v2: None,
    }
}

/// Maximum distance from point `p` to curve `c` over its parameter range
/// (sampling-based; exact for bounded ranges).
pub fn point_curve_max_extrema(c: &dyn Curve, p: &GpPnt, samples: usize) -> ExtremaPair {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let (a, b) = if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-1.0, 1.0)
    };
    let mut best: Option<ExtremaPair> = None;
    let n = samples.max(2);
    for i in 0..=n {
        let u = a + (b - a) * i as f64 / n as f64;
        let q = c.d0(u);
        let d = p.distance(&q);
        if best.as_ref().map_or(true, |bp: &ExtremaPair| d > bp.distance) {
            best = Some(ExtremaPair {
                p1: *p,
                p2: q,
                distance: d,
                u1: u,
                v1: None,
                u2: u,
                v2: None,
            });
        }
    }
    best.unwrap()
}

fn bound(c: &dyn Curve) -> (f64, f64) {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-1.0, 1.0)
    }
}

fn surface_bound(s: &dyn Surface) -> (f64, f64) {
    let (a, b) = s.u_range();
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-1.0, 1.0)
    }
}

fn w_of_j(j: usize, n: usize, s: &dyn Surface) -> f64 {
    let (a, b) = s.v_range();
    let (a, b) = if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    a + (b - a) * j as f64 / n as f64
}

/// Minimum distance between two curves (discretize + refine each local dip).
pub fn curve_curve_extrema(c1: &dyn Curve, c2: &dyn Curve, samples: usize) -> Vec<ExtremaPair> {
    let (a1, b1) = bound(c1);
    let (a2, b2) = bound(c2);
    let n = samples.max(4);
    let mut bests: Vec<ExtremaPair> = Vec::new();
    for i in 0..=n {
        let u = a1 + (b1 - a1) * i as f64 / n as f64;
        for j in 0..=n {
            let v = a2 + (b2 - a2) * j as f64 / n as f64;
            let p = c1.d0(u);
            let q = c2.d0(v);
            let d = p.distance(&q);
            bests.push(ExtremaPair {
                p1: p,
                p2: q,
                distance: d,
                u1: u,
                v1: None,
                u2: v,
                v2: None,
            });
        }
    }
    // Non-maximum suppression over the coarse grid, then refine each survivor.
    let mut seeds: Vec<(f64, f64)> = Vec::new();
    for i in 1..n {
        for j in 1..n {
            let d = bests[i * (n + 1) + j].distance;
            let neigh = [
                bests[(i - 1) * (n + 1) + j].distance,
                bests[(i + 1) * (n + 1) + j].distance,
                bests[i * (n + 1) + j - 1].distance,
                bests[i * (n + 1) + j + 1].distance,
            ];
            if d <= neigh[0] && d <= neigh[1] && d <= neigh[2] && d <= neigh[3] {
                seeds.push((a1 + (b1 - a1) * i as f64 / n as f64, a2 + (b2 - a2) * j as f64 / n as f64));
            }
        }
    }
    if seeds.is_empty() {
        let mut gmin = &bests[0];
        for e in bests.iter().skip(1) {
            if e.distance < gmin.distance {
                gmin = e;
            }
        }
        seeds.push((gmin.u1, gmin.u2));
    }
    let mut out: Vec<ExtremaPair> = Vec::new();
    for (u0, v0) in seeds {
        let (u, v) = refine_curve_curve(c1, c2, u0, v0, a1, b1, a2, b2);
        let p = c1.d0(u);
        let q = c2.d0(v);
        let d = p.distance(&q);
        if out.iter().any(|e| (e.u1 - u).abs() < 1e-5 && (e.u2 - v).abs() < 1e-5) {
            continue;
        }
        out.push(ExtremaPair {
            p1: p,
            p2: q,
            distance: d,
            u1: u,
            v1: None,
            u2: v,
            v2: None,
        });
    }
    out.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(std::cmp::Ordering::Equal));
    out
}

/// 2-D local refinement (coordinate descent on the squared distance).
fn refine_curve_curve(
    c1: &dyn Curve,
    c2: &dyn Curve,
    u0: f64,
    v0: f64,
    a1: f64,
    b1: f64,
    a2: f64,
    b2: f64,
) -> (f64, f64) {
    let dist2 = |u: f64, v: f64| {
        let p = c1.d0(u);
        let q = c2.d0(v);
        (p.x() - q.x()).powi(2) + (p.y() - q.y()).powi(2) + (p.z() - q.z()).powi(2)
    };
    let (mut u, mut v) = (u0, v0);
    let mut step = (b1 - a1).abs().max(1.0) / 4.0;
    for _ in 0..128 {
        let cur = dist2(u, v);
        let best = [
            dist2(u + step, v),
            dist2(u - step, v),
            dist2(u, v + step),
            dist2(u, v - step),
        ];
        let mut bi = 0;
        let mut bb = best[0];
        for (k, b) in best.iter().enumerate() {
            if *b < bb {
                bi = k;
                bb = *b;
            }
        }
        if bb < cur {
            match bi {
                0 => u += step,
                1 => u -= step,
                2 => v += step,
                _ => v -= step,
            }
            u = u.clamp(a1.min(b1), a1.max(b1));
            v = v.clamp(a2.min(b2), a2.max(b2));
        } else {
            step *= 0.5;
            if step < 1e-12 {
                break;
            }
        }
    }
    (u, v)
}

/// Minimum distance between a curve and a surface (grid + coordinate descent).
pub fn curve_surface_extrema(c: &dyn Curve, s: &dyn Surface, samples: usize) -> ExtremaPair {
    let (a1, b1) = bound(c);
    let (a2, b2) = surface_bound(s);
    let n = samples.max(4);
    let mut gmin: Option<ExtremaPair> = None;
    for i in 0..=n {
        let u = a1 + (b1 - a1) * i as f64 / n as f64;
        for j in 0..=n {
            let v = a2 + (b2 - a2) * j as f64 / n as f64;
            let p = c.d0(u);
            let q = s.d0(v, w_of_j(j, n, s));
            let d = p.distance(&q);
            if gmin.as_ref().map_or(true, |e: &ExtremaPair| d < e.distance) {
                gmin = Some(ExtremaPair {
                    p1: p,
                    p2: q,
                    distance: d,
                    u1: u,
                    v1: None,
                    u2: v,
                    v2: Some(w_of_j(j, n, s)),
                });
            }
        }
    }
    let g = gmin.unwrap();
    let (u, vv, ww) = refine_curve_surface(c, s, g.u1, g.u2, g.v2.unwrap(), a1, b1, a2, b2);
    let p = c.d0(u);
    let q = s.d0(vv, ww);
    ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: None,
        u2: vv,
        v2: Some(ww),
    }
}

fn refine_curve_surface(
    c: &dyn Curve,
    s: &dyn Surface,
    u0: f64,
    v0: f64,
    w0: f64,
    a1: f64,
    b1: f64,
    a2: f64,
    b2: f64,
) -> (f64, f64, f64) {
    let dist2 = |u: f64, v: f64, w: f64| {
        let p = c.d0(u);
        let q = s.d0(v, w);
        (p.x() - q.x()).powi(2) + (p.y() - q.y()).powi(2) + (p.z() - q.z()).powi(2)
    };
    let (mut u, mut v, mut w) = (u0, v0, w0);
    let mut step = (b1 - a1).abs().max(1.0) / 4.0;
    for _ in 0..256 {
        let cur = dist2(u, v, w);
        let cand = [
            (u + step, v, w),
            (u - step, v, w),
            (u, v + step, w),
            (u, v - step, w),
            (u, v, w + step),
            (u, v, w - step),
        ];
        let mut bi = usize::MAX;
        let mut bb = cur;
        for (k, (cu, cv, cw)) in cand.iter().enumerate() {
            let d = dist2(*cu, *cv, *cw);
            if d < bb {
                bi = k;
                bb = d;
            }
        }
        if bi == usize::MAX {
            step *= 0.5;
            if step < 1e-12 {
                break;
            }
        } else {
            match bi {
                0 => u += step,
                1 => u -= step,
                2 => v += step,
                3 => v -= step,
                4 => w += step,
                _ => w -= step,
            }
            u = u.clamp(a1.min(b1), a1.max(b1));
            v = v.clamp(a2.min(b2), a2.max(b2));
        }
    }
    (u, v, w)
}

/// Minimum distance between two surfaces (grid + coordinate descent).
pub fn surface_surface_extrema(s1: &dyn Surface, s2: &dyn Surface, samples: usize) -> ExtremaPair {
    let (a1, b1) = surface_bound(s1);
    let (a2, b2) = surface_bound(s2);
    let n = samples.max(3);
    let mut gmin: Option<ExtremaPair> = None;
    for i in 0..=n {
        for j in 0..=n {
            let u1 = a1 + (b1 - a1) * i as f64 / n as f64;
            let v1 = w_of_j(j, n, s1);
            let p = s1.d0(u1, v1);
            for k in 0..=n {
                for l in 0..=n {
                    let u2 = a2 + (b2 - a2) * k as f64 / n as f64;
                    let v2 = w_of_j(l, n, s2);
                    let q = s2.d0(u2, v2);
                    let d = p.distance(&q);
                    if gmin.as_ref().map_or(true, |e: &ExtremaPair| d < e.distance) {
                        gmin = Some(ExtremaPair {
                            p1: p,
                            p2: q,
                            distance: d,
                            u1,
                            v1: Some(v1),
                            u2,
                            v2: Some(v2),
                        });
                    }
                }
            }
        }
    }
    let g = gmin.unwrap();
    let (u1, v1, u2, v2) = refine_surface_surface(s1, s2, g.u1, g.v1.unwrap(), g.u2, g.v2.unwrap(), a1, b1, a2, b2);
    let p = s1.d0(u1, v1);
    let q = s2.d0(u2, v2);
    ExtremaPair {
        p1: p,
        p2: q,
        distance: p.distance(&q),
        u1,
        v1: Some(v1),
        u2,
        v2: Some(v2),
    }
}

fn refine_surface_surface(
    s1: &dyn Surface,
    s2: &dyn Surface,
    u10: f64,
    v10: f64,
    u20: f64,
    v20: f64,
    a1: f64,
    b1: f64,
    a2: f64,
    b2: f64,
) -> (f64, f64, f64, f64) {
    let dist2 = |x: f64, y: f64, z: f64, w: f64| {
        let p = s1.d0(x, y);
        let q = s2.d0(z, w);
        (p.x() - q.x()).powi(2) + (p.y() - q.y()).powi(2) + (p.z() - q.z()).powi(2)
    };
    let (mut a, mut b, mut c, mut d) = (u10, v10, u20, v20);
    let mut step = (b1 - a1).abs().max(1.0) / 4.0;
    for _ in 0..512 {
        let cur = dist2(a, b, c, d);
        let cand = [
            (a + step, b, c, d),
            (a - step, b, c, d),
            (a, b + step, c, d),
            (a, b - step, c, d),
            (a, b, c + step, d),
            (a, b, c - step, d),
            (a, b, c, d + step),
            (a, b, c, d - step),
        ];
        let mut bi = usize::MAX;
        let mut bb = cur;
        for (k, (x, y, z, w)) in cand.iter().enumerate() {
            let dd = dist2(*x, *y, *z, *w);
            if dd < bb {
                bi = k;
                bb = dd;
            }
        }
        if bi == usize::MAX {
            step *= 0.5;
            if step < 1e-12 {
                break;
            }
        } else {
            match bi {
                0 => a += step,
                1 => a -= step,
                2 => b += step,
                3 => b -= step,
                4 => c += step,
                5 => c -= step,
                6 => d += step,
                _ => d -= step,
            }
            a = a.clamp(a1.min(b1), a1.max(b1));
            b = b.clamp(a1.min(b1), a1.max(b1));
            c = c.clamp(a2.min(b2), a2.max(b2));
            d = d.clamp(a2.min(b2), a2.max(b2));
        }
    }
    (a, b, c, d)
}

/// Distance from a point to a surface (minimum), reusing the projection path.
pub fn point_surface_extrema(s: &dyn Surface, p: &GpPnt) -> ExtremaPair {
    let (a, b) = surface_bound(s);
    let (u, v) = refine_surface_point(s, p, a, b);
    let q = s.d0(u, v);
    ExtremaPair {
        p1: *p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        v1: None,
        u2: u,
        v2: Some(v),
    }
}

fn refine_surface_point(s: &dyn Surface, p: &GpPnt, a0: f64, b0: f64) -> (f64, f64) {
    let dist2 = |u: f64, v: f64| {
        let q = s.d0(u, v);
        (q.x() - p.x()).powi(2) + (q.y() - p.y()).powi(2) + (q.z() - p.z()).powi(2)
    };
    let (mut u, mut v) = (0.0, 0.0);
    let mut step = (b0 - a0).abs().max(1.0) / 4.0;
    for _ in 0..256 {
        let cur = dist2(u, v);
        let cand = [(u + step, v), (u - step, v), (u, v + step), (u, v - step)];
        let mut bi = usize::MAX;
        let mut bb = cur;
        for (k, (cu, cv)) in cand.iter().enumerate() {
            let dd = dist2(*cu, *cv);
            if dd < bb {
                bi = k;
                bb = dd;
            }
        }
        if bi == usize::MAX {
            step *= 0.5;
            if step < 1e-12 {
                break;
            }
        } else {
            match bi {
                0 => u += step,
                1 => u -= step,
                2 => v += step,
                _ => v -= step,
            }
            u = u.clamp(a0.min(b0), a0.max(b0));
            v = v.clamp(a0.min(b0), a0.max(b0));
        }
    }
    (u, v)
}

/// Build an extrema pair from a point and a tangent vector on a curve —
/// convenience for callers that already have the closest parameter.
pub fn tangent_at(c: &dyn Curve, u: f64) -> GpVec {
    let (_, d) = c.d1(u);
    d
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GeomCircle, GeomLine, GeomPlane, GeomSphere};
    use occt_core::gp::{GpAx2, GpAx3, GpCirc, GpDir, GpLin, GpPln, GpPnt, GpSphere as _GpSphere, GpTrsf, GpVec};

    fn line(p: GpPnt, d: GpDir) -> GeomLine {
        GeomLine::new(GpLin::new(occt_core::gp::GpAx1::new(p, d)))
    }

    fn sphere(r: f64) -> GeomSphere {
        GeomSphere::new(_GpSphere::new(GpAx3::standard(), r).unwrap())
    }

    #[test]
    fn point_line_min_distance() {
        let line = line(GpPnt::new(0.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let e = point_curve_extrema(&line, &GpPnt::new(3.0, 4.0, 0.0));
        assert!((e.distance - 4.0).abs() < 1e-7, "dist {}", e.distance);
        assert!((e.p2.x() - 3.0).abs() < 1e-6, "closest x {}", e.p2.x());
    }

    #[test]
    fn point_circle_min_and_max() {
        let circle = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let p = GpPnt::new(3.0, 0.0, 0.0);
        let e = point_curve_extrema(&circle, &p);
        assert!((e.distance - 2.0).abs() < 1e-6, "min dist {}", e.distance);
        let m = point_curve_max_extrema(&circle, &p, 64);
        assert!((m.distance - 4.0).abs() < 1e-6, "max dist {}", m.distance);
    }

    #[test]
    fn curve_curve_circles_min() {
        let c1 = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let mut c2 = GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0));
        let mut t = GpTrsf::identity();
        t.set_translation_vec(&GpVec::new(3.0, 0.0, 0.0));
        c2.transform(&t);
        let es = curve_curve_extrema(&c1, &c2, 24);
        assert!(!es.is_empty());
        assert!((es[0].distance - 1.0).abs() < 1e-5, "min {}", es[0].distance);
    }

    #[test]
    fn curve_surface_sphere_line() {
        // Line at y=3, sphere radius 1 at origin → min distance 2.
        let line = line(GpPnt::new(0.0, 3.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let sphere = sphere(1.0);
        let e = curve_surface_extrema(&line, &sphere, 16);
        assert!((e.distance - 2.0).abs() < 1e-4, "min {}", e.distance);
    }

    #[test]
    fn point_surface_sphere() {
        let sphere = sphere(1.0);
        let e = point_surface_extrema(&sphere, &GpPnt::new(3.0, 0.0, 0.0));
        assert!((e.distance - 2.0).abs() < 1e-5, "min {}", e.distance);
    }

    #[test]
    fn surface_surface_planes_parallel() {
        let ax = GpAx3::standard();
        let p1 = GeomPlane::new(GpPln::new(ax.clone()));
        let p2 = GeomPlane::new(GpPln::new(ax));
        let e = surface_surface_extrema(&p1, &p2, 3);
        assert!(e.distance.abs() < 1e-6, "coincident planes {}", e.distance);
    }

    #[test]
    fn surface_surface_sphere_plane() {
        let sphere = sphere(1.0);
        let ax = GpAx3::standard();
        let mut pl = GpPln::new(ax);
        pl.set_location(&GpPnt::new(0.0, 0.0, 5.0));
        let plane = GeomPlane::new(pl);
        let e = surface_surface_extrema(&sphere, &plane, 4);
        assert!((e.distance - 4.0).abs() < 1e-4, "min {}", e.distance);
    }

    #[test]
    fn curve_curve_skew_lines() {
        let l1 = line(GpPnt::new(0.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
        let l2 = line(GpPnt::new(0.0, 0.0, 3.0), GpDir::new(0.0, 1.0, 0.0).unwrap());
        let es = curve_curve_extrema(&l1, &l2, 8);
        assert!(!es.is_empty());
        assert!((es[0].distance - 3.0).abs() < 1e-5, "min {}", es[0].distance);
    }
}
