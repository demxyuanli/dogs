//! Phase 4 module: geom2d_api — high-level 2D geometry algorithms.
//! Ports of `Geom2dAPI_ProjectPointOnCurve` and `Geom2dAPI_InterCurveCurve`.

use occt_core::gp::{GpCirc2d, GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use crate::curve::Curve2d;
use crate::curve_ops::{curve2d_closest_point, curve2d_intersections, lin2d_intersection};

/// Result of projecting a point onto a curve.
#[derive(Debug, Clone)]
pub struct Project2d {
    pub parameter: f64,
    pub point: GpPnt2d,
    pub distance: f64,
}

/// Project `p` onto `c`. Wraps [`curve2d_closest_point`] with a bounded Newton
/// refinement on the tangency condition `(c(u) - p)·c'(u) = 0`.
pub fn project_point_on_curve(c: &dyn Curve2d, p: &GpPnt2d, tol: f64) -> Option<Project2d> {
    let (mut u, mut q) = curve2d_closest_point(c, p, tol)?;
    let a = c.first_parameter();
    let b = c.last_parameter();
    if a.is_finite() && b.is_finite() {
        for _ in 0..8 {
            let (p0, d1) = c.d1(u);
            let (_, _, d2) = c.d2(u);
            let dx = p0.x() - p.x();
            let dy = p0.y() - p.y();
            let f = dx * d1.x() + dy * d1.y();
            let df = d1.x() * d1.x() + d1.y() * d1.y() + dx * d2.x() + dy * d2.y();
            if df.abs() < 1e-15 {
                break;
            }
            let nu = u - f / df;
            if !nu.is_finite() || nu < a || nu > b {
                break;
            }
            if (nu - u).abs() < 1e-14 {
                u = nu;
                break;
            }
            u = nu;
        }
        q = c.d0(u);
    }
    Some(Project2d {
        parameter: u,
        point: q,
        distance: q.distance(p),
    })
}

/// Unit tangent vector at parameter `u`.
pub fn tangent2d(c: &dyn Curve2d, u: f64) -> GpVec2d {
    let (_, v) = c.d1(u);
    v.normalized().unwrap_or_else(|_| GpVec2d::new(1.0, 0.0))
}

/// Unit normal vector at parameter `u` (tangent rotated +90°).
pub fn normal2d(c: &dyn Curve2d, u: f64) -> GpVec2d {
    let t = tangent2d(c, u);
    GpVec2d::new(-t.y(), t.x())
}

/// A point shared by two curves, with the parameter on each.
#[derive(Debug, Clone)]
pub struct CurveIntersection2d {
    pub u1: f64,
    pub u2: f64,
    pub point: GpPnt2d,
}

/// Intersections of two curves. Exact for line/line (via [`intersect_line_line`]);
/// otherwise falls back to the sampler-based [`curve2d_intersections`].
pub fn intersect_curves(a: &dyn Curve2d, b: &dyn Curve2d, tol: f64) -> Vec<CurveIntersection2d> {
    if let (Some(l1), Some(l2)) = (curve_as_line(a), curve_as_line(b)) {
        return match intersect_line_line(&l1, &l2) {
            Some(p) => vec![CurveIntersection2d {
                u1: param_on_line(&l1, &p),
                u2: param_on_line(&l2, &p),
                point: p,
            }],
            None => Vec::new(), // parallel or coincident
        };
    }
    curve2d_intersections(a, b, tol)
        .into_iter()
        .map(|(u1, u2, point)| CurveIntersection2d { u1, u2, point })
        .collect()
}

/// Exact intersection of two infinite lines. `None` when parallel or coincident.
pub fn intersect_line_line(l1: &GpLin2d, l2: &GpLin2d) -> Option<GpPnt2d> {
    lin2d_intersection(l1, l2)
}

/// Exact analytic circle/line intersections. Solve `|w + t·d|² = r²` with `d`
/// unit; returns 0, 1, or 2 points ordered by signed distance along the line
/// from its origin.
pub fn intersect_circle_line(c: &GpCirc2d, l: &GpLin2d) -> Vec<GpPnt2d> {
    let (cx, cy) = (c.location().x(), c.location().y());
    let r = c.radius();
    let (lx, ly) = (l.pos.loc.x(), l.pos.loc.y());
    let (dx, dy) = (l.pos.vdir.x, l.pos.vdir.y);
    // w = line origin − center; line points are L + t·d.
    let (wx, wy) = (lx - cx, ly - cy);
    let wd = wx * dx + wy * dy;
    let disc = r * r - (wx * wx + wy * wy) + wd * wd;
    if disc < -1e-12 {
        return Vec::new();
    }
    let s = disc.max(0.0).sqrt();
    if s < 1e-12 {
        let t = -wd;
        vec![GpPnt2d::new(lx + t * dx, ly + t * dy)]
    } else {
        let (t0, t1) = (-wd - s, -wd + s);
        vec![
            GpPnt2d::new(lx + t0 * dx, ly + t0 * dy),
            GpPnt2d::new(lx + t1 * dx, ly + t1 * dy),
        ]
    }
}

/// Exact analytic circle/circle intersections via the radical line. Concentric
/// circles give an empty result (0 points); otherwise 0, 1, or 2 points.
pub fn intersect_circle_circle(c1: &GpCirc2d, c2: &GpCirc2d) -> Vec<GpPnt2d> {
    let (x1, y1) = (c1.location().x(), c1.location().y());
    let (x2, y2) = (c2.location().x(), c2.location().y());
    let r1 = c1.radius().abs();
    let r2 = c2.radius().abs();
    let (dx, dy) = (x2 - x1, y2 - y1);
    let d2 = dx * dx + dy * dy;
    let d = d2.sqrt();
    if d < 1e-15 {
        return Vec::new(); // concentric
    }
    if d > r1 + r2 + 1e-12 || d < (r1 - r2).abs() - 1e-12 {
        return Vec::new(); // separated, or one circle inside the other
    }
    // Foot of the radical axis along C1→C2, then perpendicular offset h.
    let a = (r1 * r1 - r2 * r2 + d2) / (2.0 * d);
    let h = (r1 * r1 - a * a).max(0.0).sqrt();
    let (ux, uy) = (dx / d, dy / d);
    let (px, py) = (x1 + a * ux, y1 + a * uy);
    if h < 1e-12 {
        vec![GpPnt2d::new(px, py)]
    } else {
        let (nx, ny) = (-uy, ux);
        let mut v = vec![
            GpPnt2d::new(px + h * nx, py + h * ny),
            GpPnt2d::new(px - h * nx, py - h * ny),
        ];
        v.sort_by(|p, q| p.x().partial_cmp(&q.x()).unwrap().then(p.y().partial_cmp(&q.y()).unwrap()));
        v
    }
}

/// Closest pair of points between two bounded curves: min `|c1(u) - c2(v)|` by
/// a coarse grid over both parameters, then alternating golden-section refine.
///
/// ponytail: unbounded curves (lines) return `None` — expand the window like
/// `curve_ops::curve2d_closest_point` if needed.
pub fn closest_points(c1: &dyn Curve2d, c2: &dyn Curve2d, _tol: f64) -> Option<((f64, f64), f64)> {
    let (a1, b1) = (c1.first_parameter(), c1.last_parameter());
    let (a2, b2) = (c2.first_parameter(), c2.last_parameter());
    if !a1.is_finite() || !b1.is_finite() || !a2.is_finite() || !b2.is_finite() {
        return None;
    }
    let n = 64;
    let mut best: Option<(f64, f64, f64)> = None;
    for i in 0..=n {
        let u = a1 + (b1 - a1) * i as f64 / n as f64;
        for j in 0..=n {
            let v = a2 + (b2 - a2) * j as f64 / n as f64;
            let d = c1.d0(u).square_distance(&c2.d0(v));
            if best.as_ref().map_or(true, |(_, _, bd)| d < *bd) {
                best = Some((u, v, d));
            }
        }
    }
    let (mut u, mut v, _) = best?;
    let du = (b1 - a1) / n as f64;
    let dv = (b2 - a2) / n as f64;
    for _ in 0..10 {
        let (nu, _) = minimize_1d(
            &|t| c1.d0(t).square_distance(&c2.d0(v)),
            (u - du).max(a1),
            (u + du).min(b1),
        );
        u = nu;
        let (nv, _) = minimize_1d(
            &|t| c1.d0(u).square_distance(&c2.d0(t)),
            (v - dv).max(a2),
            (v + dv).min(b2),
        );
        v = nv;
    }
    Some(((u, v), c1.d0(u).distance(&c2.d0(v))))
}

/// Signed area swept by the curve over `[a, b]`: `∫(x·dy − y·dx)/2` by the
/// trapezoid rule on the integrand `(x·y' − y·x')/2`.
pub fn curve2d_signed_area(c: &dyn Curve2d, a: f64, b: f64, samples: usize) -> f64 {
    let n = samples.max(2);
    let h = (b - a) / n as f64;
    let integrand = |u: f64| {
        let p = c.d0(u);
        let (_, d) = c.d1(u);
        0.5 * (p.x() * d.y() - p.y() * d.x())
    };
    let mut sum = 0.5 * (integrand(a) + integrand(b));
    for i in 1..n {
        sum += integrand(a + i as f64 * h);
    }
    sum * h
}

/// Axis-aligned bounding box of a bounded curve by sampling: `(min, max)`.
/// Unbounded curves return an empty/infinite box.
pub fn curve2d_bbox(c: &dyn Curve2d, samples: usize) -> (GpPnt2d, GpPnt2d) {
    let a = c.first_parameter();
    let b = c.last_parameter();
    if !a.is_finite() || !b.is_finite() {
        return (
            GpPnt2d::new(f64::INFINITY, f64::INFINITY),
            GpPnt2d::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
    }
    let n = samples.max(2);
    let mut minx = f64::INFINITY;
    let mut maxx = f64::NEG_INFINITY;
    let mut miny = f64::INFINITY;
    let mut maxy = f64::NEG_INFINITY;
    for i in 0..=n {
        let q = c.d0(a + (b - a) * i as f64 / n as f64);
        minx = minx.min(q.x());
        maxx = maxx.max(q.x());
        miny = miny.min(q.y());
        maxy = maxy.max(q.y());
    }
    (GpPnt2d::new(minx, miny), GpPnt2d::new(maxx, maxy))
}

// --- internals -------------------------------------------------------------

/// Extract a `GpLin2d` from a curve that is geometrically a straight line
/// (second derivative identically zero, e.g. `Geom2dLine` or a line offset).
fn curve_as_line(c: &dyn Curve2d) -> Option<GpLin2d> {
    let (_, _, d2) = c.d2(0.0);
    if d2.square_magnitude() > 1e-20 {
        return None;
    }
    let p = c.d0(0.0);
    let (_, d1) = c.d1(0.0);
    let dir = GpDir2d::new(d1.x(), d1.y()).ok()?;
    Some(GpLin2d::from_pnt_dir(p, dir))
}

/// Parameter `u` on a unit-direction line (signed distance from its origin).
fn param_on_line(line: &GpLin2d, p: &GpPnt2d) -> f64 {
    (p.x() - line.pos.loc.x()) * line.pos.vdir.x + (p.y() - line.pos.loc.y()) * line.pos.vdir.y
}

/// Golden-section minimization of `f` over `[lo, hi]`. Returns `(argmin, min)`.
/// Local copy of `curve_ops::minimize_1d`, which is private to that module.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circle::Geom2dCircle;
    use crate::line::Geom2dLine;
    use occt_core::gp::GpAx22d;

    fn circle(center: (f64, f64), r: f64) -> Geom2dCircle {
        let ax = GpAx22d::from_xdir(GpPnt2d::new(center.0, center.1), GpDir2d::default());
        Geom2dCircle::new(GpCirc2d::new(ax, r))
    }
    fn line(p: (f64, f64), d: (f64, f64)) -> Geom2dLine {
        Geom2dLine::from_pnt_dir(GpPnt2d::new(p.0, p.1), GpDir2d::new(d.0, d.1).unwrap())
    }

    #[test]
    fn circle_line_radius2_exact() {
        let c = GpCirc2d::new(GpAx22d::standard(), 2.0);
        let l = GpLin2d::from_pnt_dir(GpPnt2d::zero(), GpDir2d::new(1.0, 0.0).unwrap()); // y = 0
        let hits = intersect_circle_line(&c, &l);
        assert_eq!(hits.len(), 2, "hits={hits:?}");
        let mut xs: Vec<f64> = hits.iter().map(|p| p.x()).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((xs[0] + 2.0).abs() < 1e-9, "x0={}", xs[0]);
        assert!((xs[1] - 2.0).abs() < 1e-9, "x1={}", xs[1]);
        for p in &hits {
            assert!(p.y().abs() < 1e-9);
        }
    }

    #[test]
    fn circle_line_radius2_through_intersect_curves() {
        let c = circle((0.0, 0.0), 2.0);
        let l = line((0.0, 0.0), (1.0, 0.0));
        let hits = intersect_curves(&c as &dyn Curve2d, &l as &dyn Curve2d, 1e-9);
        assert_eq!(hits.len(), 2, "hits={hits:?}");
        let mut xs: Vec<f64> = hits.iter().map(|h| h.point.x()).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((xs[0] + 2.0).abs() < 1e-9, "x0={}", xs[0]);
        assert!((xs[1] - 2.0).abs() < 1e-9, "x1={}", xs[1]);
    }

    #[test]
    fn circle_circle_two_points() {
        let c1 = GpCirc2d::new(GpAx22d::standard(), 2.0);
        let c2 = GpCirc2d::new(
            GpAx22d::from_xdir(GpPnt2d::new(3.0, 0.0), GpDir2d::default()),
            2.0,
        );
        let hits = intersect_circle_circle(&c1, &c2);
        assert_eq!(hits.len(), 2, "hits={hits:?}");
        for p in &hits {
            assert!((p.x() - 1.5).abs() < 1e-9, "x={}", p.x());
        }
        let ys: Vec<f64> = hits.iter().map(|p| p.y()).collect();
        assert!((ys[0].abs() - 1.3228756555).abs() < 1e-6, "ys={ys:?}");
    }

    #[test]
    fn circle_circle_concentric_and_separated() {
        let c1 = GpCirc2d::new(GpAx22d::standard(), 2.0);
        let c2 = GpCirc2d::new(GpAx22d::standard(), 1.0);
        assert!(intersect_circle_circle(&c1, &c2).is_empty(), "concentric");
        let c3 = GpCirc2d::new(
            GpAx22d::from_xdir(GpPnt2d::new(5.0, 0.0), GpDir2d::default()),
            1.0,
        );
        assert!(intersect_circle_circle(&c1, &c3).is_empty(), "separated");
    }

    #[test]
    fn circle_line_no_hit() {
        let c = GpCirc2d::new(GpAx22d::standard(), 1.0);
        let l = GpLin2d::from_pnt_dir(GpPnt2d::new(2.0, 0.0), GpDir2d::new(0.0, 1.0).unwrap()); // x = 2
        assert!(intersect_circle_line(&c, &l).is_empty());

        let gc = circle((0.0, 0.0), 1.0);
        let gl = line((2.0, 0.0), (0.0, 1.0));
        assert!(intersect_curves(&gc as &dyn Curve2d, &gl as &dyn Curve2d, 1e-9).is_empty());
    }

    #[test]
    fn line_line_exact_and_parallel() {
        let l1 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 1.0), GpDir2d::new(1.0, 0.0).unwrap()); // y = 1
        let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(1.0, 0.0), GpDir2d::new(0.0, 1.0).unwrap()); // x = 1
        let p = intersect_line_line(&l1, &l2).expect("lines cross");
        assert!((p.x() - 1.0).abs() < 1e-12 && (p.y() - 1.0).abs() < 1e-12);

        let par = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 5.0), GpDir2d::new(1.0, 0.0).unwrap());
        assert!(intersect_line_line(&l1, &par).is_none(), "parallel");

        let g1 = line((0.0, 1.0), (1.0, 0.0));
        let g2 = line((1.0, 0.0), (0.0, 1.0));
        let hits = intersect_curves(&g1 as &dyn Curve2d, &g2 as &dyn Curve2d, 1e-9);
        assert_eq!(hits.len(), 1, "hits={hits:?}");
        assert!((hits[0].point.x() - 1.0).abs() < 1e-12 && (hits[0].point.y() - 1.0).abs() < 1e-12);
        assert!((hits[0].u1 - 1.0).abs() < 1e-9, "u1={}", hits[0].u1);
        assert!((hits[0].u2 - 1.0).abs() < 1e-9, "u2={}", hits[0].u2);

        let gpar = line((0.0, 5.0), (1.0, 0.0));
        assert!(intersect_curves(&g1 as &dyn Curve2d, &gpar as &dyn Curve2d, 1e-9).is_empty());
    }

    #[test]
    fn project_point_on_circle() {
        let c = circle((0.0, 0.0), 1.0);
        let r = project_point_on_curve(&c as &dyn Curve2d, &GpPnt2d::new(0.0, 2.0), 1e-9).unwrap();
        assert!((r.distance - 1.0).abs() < 1e-9, "distance={}", r.distance);
        assert!((r.point.x() - 0.0).abs() < 1e-9 && (r.point.y() - 1.0).abs() < 1e-9, "point={:?}", r.point);
    }

    #[test]
    fn signed_area_unit_circle() {
        let c = circle((0.0, 0.0), 1.0);
        let area = curve2d_signed_area(&c as &dyn Curve2d, 0.0, 2.0 * std::f64::consts::PI, 1000);
        assert!((area - std::f64::consts::PI).abs() < 1e-9, "area={area}");
    }

    #[test]
    fn tangent_normal_and_bbox() {
        let c = circle((0.0, 0.0), 1.0);
        let t = tangent2d(&c as &dyn Curve2d, 0.0);
        assert!((t.x() - 0.0).abs() < 1e-12 && (t.y() - 1.0).abs() < 1e-12, "t={t:?}");
        let n = normal2d(&c as &dyn Curve2d, 0.0);
        assert!((n.x() + 1.0).abs() < 1e-12 && (n.y() - 0.0).abs() < 1e-12, "n={n:?}");

        let (lo, hi) = curve2d_bbox(&c as &dyn Curve2d, 256);
        assert!((lo.x() + 1.0).abs() < 1e-6 && (lo.y() + 1.0).abs() < 1e-6, "lo={lo:?}");
        assert!((hi.x() - 1.0).abs() < 1e-6 && (hi.y() - 1.0).abs() < 1e-6, "hi={hi:?}");
    }
}
