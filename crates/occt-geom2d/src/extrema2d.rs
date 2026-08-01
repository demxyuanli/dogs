//! 2D extrema — point/curve minimum and maximum distances in the plane.
//!
//! Port of `Extrema_ExtPC2d`, `Extrema_ExtCC2d` (discretization + local
//! refinement). Source: `Extrema` (TKGeomBase).

use occt_core::gp::{GpPnt2d, GpVec2d};

use crate::curve::Curve2d;

/// A solved extremum between two 2D objects.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Extrema2d {
    pub p1: GpPnt2d,
    pub p2: GpPnt2d,
    pub distance: f64,
    pub u1: f64,
    pub u2: f64,
}

fn bound(c: &dyn Curve2d) -> (f64, f64) {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-1.0, 1.0)
    }
}

/// Minimum distance from point `p` to curve `c` (with closest point).
/// Golden-section refinement over the parameter range; unbounded curves use
/// an expanding window until a strict interior minimum is bracketed.
pub fn point_curve_extrema2d(c: &dyn Curve2d, p: &GpPnt2d) -> Extrema2d {
    let f = |u: f64| {
        let q = c.d0(u);
        (q.x() - p.x()).powi(2) + (q.y() - p.y()).powi(2)
    };
    let (a0, b0) = (c.first_parameter(), c.last_parameter());
    let (a, b) = if a0.is_finite() && b0.is_finite() && b0 > a0 {
        (a0, b0)
    } else {
        // Unbounded: expand until the minimum is interior (both endpoints ≥ center).
        let mut w = 1.0;
        let mut fm = f(0.0);
        for _ in 0..10 {
            if f(-w) >= fm && f(w) >= fm {
                break;
            }
            w *= 8.0;
            fm = f(0.0);
        }
        (-w, w)
    };
    // Golden-section refinement (golden ratio conjugate 0.618...).
    let phi = (5.0f64.sqrt() - 1.0) * 0.5;
    let (mut lo, mut hi) = (a, b);
    for _ in 0..80 {
        let x1 = hi - phi * (hi - lo);
        let x2 = lo + phi * (hi - lo);
        if f(x1) < f(x2) {
            hi = x2;
        } else {
            lo = x1;
        }
    }
    let u = 0.5 * (lo + hi);
    let q = c.d0(u);
    Extrema2d {
        p1: *p,
        p2: q,
        distance: p.distance(&q),
        u1: u,
        u2: u,
    }
}

/// Maximum distance from point `p` to curve `c` (sampling-based).
pub fn point_curve_max_extrema2d(c: &dyn Curve2d, p: &GpPnt2d, samples: usize) -> Extrema2d {
    let (a, b) = bound(c);
    let n = samples.max(2);
    let mut best: Option<Extrema2d> = None;
    for i in 0..=n {
        let u = a + (b - a) * i as f64 / n as f64;
        let q = c.d0(u);
        let d = p.distance(&q);
        if best.as_ref().map_or(true, |e: &Extrema2d| d > e.distance) {
            best = Some(Extrema2d { p1: *p, p2: q, distance: d, u1: u, u2: u });
        }
    }
    best.unwrap()
}

/// Minimum distance between two curves (grid + coordinate descent).
pub fn curve_curve_extrema2d(c1: &dyn Curve2d, c2: &dyn Curve2d, samples: usize) -> Vec<Extrema2d> {
    let (a1, b1) = bound(c1);
    let (a2, b2) = bound(c2);
    let n = samples.max(4);
    let mut gmin: Option<Extrema2d> = None;
    for i in 0..=n {
        let u = a1 + (b1 - a1) * i as f64 / n as f64;
        for j in 0..=n {
            let v = a2 + (b2 - a2) * j as f64 / n as f64;
            let p = c1.d0(u);
            let q = c2.d0(v);
            let d = p.distance(&q);
            if gmin.as_ref().map_or(true, |e: &Extrema2d| d < e.distance) {
                gmin = Some(Extrema2d { p1: p, p2: q, distance: d, u1: u, u2: v });
            }
        }
    }
    let g = gmin.unwrap();
    let (u, v) = refine_curve_curve(c1, c2, g.u1, g.u2, a1, b1, a2, b2);
    let p = c1.d0(u);
    let q = c2.d0(v);
    vec![Extrema2d { p1: p, p2: q, distance: p.distance(&q), u1: u, u2: v }]
}

fn refine_curve_curve(
    c1: &dyn Curve2d,
    c2: &dyn Curve2d,
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
        (p.x() - q.x()).powi(2) + (p.y() - q.y()).powi(2)
    };
    let (mut u, mut v) = (u0, v0);
    let mut step = (b1 - a1).abs().max(1.0) / 4.0;
    for _ in 0..128 {
        let cur = dist2(u, v);
        let cand = [(u + step, v), (u - step, v), (u, v + step), (u, v - step)];
        let mut bi = usize::MAX;
        let mut bb = cur;
        for (k, (cu, cv)) in cand.iter().enumerate() {
            let d = dist2(*cu, *cv);
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
                _ => v -= step,
            }
            u = u.clamp(a1.min(b1), a1.max(b1));
            v = v.clamp(a2.min(b2), a2.max(b2));
        }
    }
    (u, v)
}

/// All curve–curve intersection points (zero-distance extrema).
pub fn curve_curve_intersections2d(c1: &dyn Curve2d, c2: &dyn Curve2d, tol: f64) -> Vec<GpPnt2d> {
    let (a1, b1) = bound(c1);
    let (a2, b2) = bound(c2);
    let n = 32;
    let mut pts: Vec<GpPnt2d> = Vec::new();
    for i in 0..=n {
        let u = a1 + (b1 - a1) * i as f64 / n as f64;
        for j in 0..=n {
            let v = a2 + (b2 - a2) * j as f64 / n as f64;
            let p = c1.d0(u);
            let q = c2.d0(v);
            if p.distance(&q) < tol {
                if pts.iter().all(|x| x.distance(&p) > 1e-6) {
                    pts.push(p);
                }
            }
        }
    }
    pts
}

/// Distance from a point to the closest point on a polyline.
pub fn point_polyline_extrema2d(poly: &[GpPnt2d], p: &GpPnt2d) -> Extrema2d {
    let mut best: Option<Extrema2d> = None;
    for w in poly.windows(2) {
        let (a, b) = (w[0], w[1]);
        let abx = b.x() - a.x();
        let aby = b.y() - a.y();
        let len2 = abx * abx + aby * aby;
        let t = if len2 < 1e-30 {
            0.0
        } else {
            (((p.x() - a.x()) * abx + (p.y() - a.y()) * aby) / len2).clamp(0.0, 1.0)
        };
        let q = GpPnt2d::new(a.x() + t * abx, a.y() + t * aby);
        let d = p.distance(&q);
        if best.as_ref().map_or(true, |e: &Extrema2d| d < e.distance) {
            best = Some(Extrema2d { p1: *p, p2: q, distance: d, u1: t, u2: t });
        }
    }
    best.unwrap()
}

/// Tangent direction of a 2D curve at `u` (normalized).
pub fn tangent2d(c: &dyn Curve2d, u: f64) -> GpVec2d {
    let (_, d) = c.d1(u);
    let m = d.magnitude();
    if m > 1e-12 {
        GpVec2d::new(d.x() / m, d.y() / m)
    } else {
        GpVec2d::new(1.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Geom2dLine, Geom2dCircle};
    use occt_core::gp::{GpAx2d, GpCirc2d, GpDir2d, GpPnt2d, GpVec2d};

    fn line2d() -> Geom2dLine {
        Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap()))
    }

    fn circle2d() -> Geom2dCircle {
        Geom2dCircle::new(GpCirc2d::new(occt_core::gp::GpAx22d::standard(), 1.0))
    }

    #[test]
    fn point_line_min_distance() {
        let c = line2d();
        let e = point_curve_extrema2d(&c, &GpPnt2d::new(3.0, 4.0));
        assert!((e.distance - 4.0).abs() < 1e-6, "dist {}", e.distance);
        assert!((e.p2.x() - 3.0).abs() < 1e-5, "closest x {}", e.p2.x());
    }

    #[test]
    fn point_circle_min_and_max() {
        let c = circle2d();
        let p = GpPnt2d::new(3.0, 0.0);
        let e = point_curve_extrema2d(&c, &p);
        assert!((e.distance - 2.0).abs() < 1e-6, "min {}", e.distance);
        let m = point_curve_max_extrema2d(&c, &p, 64);
        assert!((m.distance - 4.0).abs() < 1e-6, "max {}", m.distance);
    }

    #[test]
    fn curve_curve_circles_min() {
        // Two unit circles centers 3 apart → min distance 1.
        let c1 = circle2d();
        let c2 = Geom2dCircle::new(GpCirc2d::new(occt_core::gp::GpAx22d::standard(), 1.0).translated_vec(&GpVec2d::new(3.0, 0.0)));
        let es = curve_curve_extrema2d(&c1, &c2, 24);
        assert!(!es.is_empty());
        assert!((es[0].distance - 1.0).abs() < 1e-5, "min {}", es[0].distance);
    }

    #[test]
    fn curve_curve_intersections_found() {
        // Circle at origin and horizontal line y=0 → 2 intersection points.
        let c1 = circle2d();
        let l = line2d();
        let pts = curve_curve_intersections2d(&c1, &l, 1e-3);
        assert_eq!(pts.len(), 2, "circle∩x-axis points {pts:?}");
        for p in &pts {
            assert!((p.x().abs() - 1.0).abs() < 0.05, "on circle x {}", p.x());
        }
    }

    #[test]
    fn point_polyline_closest() {
        let poly = vec![GpPnt2d::new(0.0, 0.0), GpPnt2d::new(2.0, 0.0)];
        let e = point_polyline_extrema2d(&poly, &GpPnt2d::new(1.0, 3.0));
        assert!((e.distance - 3.0).abs() < 1e-9, "dist {}", e.distance);
        assert!((e.p2.x() - 1.0).abs() < 1e-9 && e.p2.y().abs() < 1e-9);
    }

    #[test]
    fn point_polyline_endpoint_t() {
        let poly = vec![GpPnt2d::new(0.0, 0.0), GpPnt2d::new(2.0, 0.0)];
        // Point beyond the segment end clamps to the endpoint.
        let e = point_polyline_extrema2d(&poly, &GpPnt2d::new(5.0, 1.0));
        assert!((e.p2.x() - 2.0).abs() < 1e-9, "clamp x {}", e.p2.x());
    }

    #[test]
    fn tangent_horizontal_line() {
        let c = line2d();
        let t = tangent2d(&c, 0.5);
        assert!((t.x() - 1.0).abs() < 1e-9 && t.y().abs() < 1e-9, "tangent {t:?}");
    }

    #[test]
    fn skew_lines_min_distance() {
        // Horizontal line and vertical line offset — skew in the plane means
        // they cross; min distance 0. Use parallel lines instead: y=0 and y=3.
        let l1 = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 0.0), GpDir2d::new(1.0, 0.0).unwrap()));
        let l2 = Geom2dLine::new(GpAx2d::new(GpPnt2d::new(0.0, 3.0), GpDir2d::new(1.0, 0.0).unwrap()));
        let es = curve_curve_extrema2d(&l1, &l2, 8);
        assert!((es[0].distance - 3.0).abs() < 1e-5, "min {}", es[0].distance);
    }
}
