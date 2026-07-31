//! 2D curve operations: arc length, closest point, intersections.
//! Analytic ports of the line/segment cases (gp_Lin2d, IntAna2d) plus a
//! sampler-based fallback for arbitrary curves (Geom2dAPI_ProjectPointOnCurve).

use occt_core::gp::{GpLin2d, GpPnt2d, GpXY};
use crate::curve::Curve2d;

/// Arc length of a bounded curve by Simpson quadrature of `|d1|`.
/// Returns `NaN` for unbounded curves.
pub fn curve2d_length(c: &dyn Curve2d, n: usize) -> f64 {
    let a = c.first_parameter();
    let b = c.last_parameter();
    if !a.is_finite() || !b.is_finite() {
        return f64::NAN;
    }
    let integrand = |u: f64| c.d1(u).1.magnitude();
    let n = (n / 2 * 2).max(2);
    let h = (b - a) / n as f64;
    let mut sum = integrand(a) + integrand(b);
    for i in 1..n {
        let x = a + i as f64 * h;
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(x);
    }
    sum * h / 3.0
}

/// Point on the curve at parameter `u`.
pub fn curve2d_point(c: &dyn Curve2d, u: f64) -> GpPnt2d {
    c.d0(u)
}

/// Closest point on the curve to `p`. Returns `(parameter, point)`.
/// Works on bounded curves; unbounded curves (lines) are handled by probing an
/// expanding window around `u = 0` until the minimum is interior.
pub fn curve2d_closest_point(c: &dyn Curve2d, p: &GpPnt2d, _tol: f64) -> Option<(f64, GpPnt2d)> {
    let a = c.first_parameter();
    let b = c.last_parameter();
    if !a.is_finite() || !b.is_finite() {
        let base = c.d0(0.0);
        let d = base.distance(p);
        let mut w = (d + 1.0).max(1.0) * 8.0;
        for _ in 0..4 {
            if let Some((u, q)) = refine_closest(c, p, -w, w) {
                if (u + w).abs() > 1e-9 && (u - w).abs() > 1e-9 {
                    return Some((u, q));
                }
            }
            w *= 8.0;
        }
        return refine_closest(c, p, -w, w);
    }
    refine_closest(c, p, a, b)
}

/// Distance from `p` to the nearest point on the curve.
pub fn curve2d_distance_to_point(c: &dyn Curve2d, p: &GpPnt2d, tol: f64) -> f64 {
    curve2d_closest_point(c, p, tol).map_or(f64::INFINITY, |(_, q)| q.distance(p))
}

/// Approximate intersection points of two curves.
/// Returns `(u, v, point)` where `u`/`v` are the parameters on `a`/`b`.
///
/// ponytail: approximate sampler-based intersection — not root-exact. Adequate
/// for geometry tooling; replace with a subdivision/Newton solver (IntAna2d)
/// if exact roots on arbitrary curves are required.
pub fn curve2d_intersections(a: &dyn Curve2d, b: &dyn Curve2d, tol: f64) -> Vec<(f64, f64, GpPnt2d)> {
    let tol = tol.max(1e-12);
    let na = 256;
    let nb = 256;
    let sa = sample_curve(a, b, na);
    let sb = sample_curve(b, a, nb);
    if sa.is_empty() || sb.is_empty() {
        return Vec::new();
    }
    let ua0 = sa[0].0;
    let ua1 = sa[sa.len() - 1].0;
    let vb0 = sb[0].0;
    let vb1 = sb[sb.len() - 1].0;
    let du = ((ua1 - ua0) / na as f64).max(1e-12);
    let dv = ((vb1 - vb0) / nb as f64).max(1e-12);

    let coarse = tol * 100.0 + 1e-9;
    let mut found: Vec<(f64, f64)> = Vec::new();
    for (ui, ai) in &sa {
        for (vj, bj) in &sb {
            if ai.square_distance(bj) < coarse * coarse {
                // Refine the candidate pair by alternating 1-D minimization.
                let mut u = *ui;
                let mut v = *vj;
                for _ in 0..6 {
                    let (nu, _) = minimize_1d(
                        &|t| a.d0(t).square_distance(&b.d0(v)),
                        (u - du).max(ua0),
                        (u + du).min(ua1),
                    );
                    u = nu;
                    let (nv, _) = minimize_1d(
                        &|t| a.d0(u).square_distance(&b.d0(t)),
                        (v - dv).max(vb0),
                        (v + dv).min(vb1),
                    );
                    v = nv;
                }
                if a.d0(u).distance(&b.d0(v)) <= tol {
                    found.push((u, v));
                }
            }
        }
    }

    // Dedupe candidates that converge to the same geometric point.
    let mut out: Vec<(f64, f64, GpPnt2d)> = Vec::new();
    for (u, v) in found {
        let pa = a.d0(u);
        let pb = b.d0(v);
        let mid = GpPnt2d::new((pa.x() + pb.x()) * 0.5, (pa.y() + pb.y()) * 0.5);
        let dup = out.iter().any(|(_, _, p)| p.distance(&mid) <= tol.max(1e-6));
        if !dup {
            out.push((u, v, mid));
        }
    }
    out
}

/// Exact intersection of two infinite lines (cross-product formula).
/// `None` when the lines are parallel.
pub fn lin2d_intersection(l1: &GpLin2d, l2: &GpLin2d) -> Option<GpPnt2d> {
    let d1 = l1.pos.vdir;
    let d2 = l2.pos.vdir;
    let denom = d1.crossed(&d2);
    if denom.abs() < 1e-15 {
        return None;
    }
    let w = l2.pos.loc.xy().subtracted(&l1.pos.loc.xy());
    let d2v = GpXY::new(d2.x, d2.y);
    let t = w.crossed(&d2v) / denom;
    Some(GpPnt2d::new(l1.pos.loc.x() + t * d1.x, l1.pos.loc.y() + t * d1.y))
}

/// Exact intersection of segments `p1-p2` and `p3-p4`.
/// `None` when parallel or the intersection falls outside either segment.
pub fn segment_intersection(p1: &GpPnt2d, p2: &GpPnt2d, p3: &GpPnt2d, p4: &GpPnt2d) -> Option<GpPnt2d> {
    let d1 = GpXY::new(p2.x() - p1.x(), p2.y() - p1.y());
    let d2 = GpXY::new(p4.x() - p3.x(), p4.y() - p3.y());
    let denom = d1.crossed(&d2);
    if denom.abs() < 1e-15 {
        return None;
    }
    let w = GpXY::new(p3.x() - p1.x(), p3.y() - p1.y());
    let t = w.crossed(&d2) / denom;
    let s = w.crossed(&d1) / denom;
    let eps = 1e-12;
    if t >= -eps && t <= 1.0 + eps && s >= -eps && s <= 1.0 + eps {
        Some(GpPnt2d::new(p1.x() + t * d1.x, p1.y() + t * d1.y))
    } else {
        None
    }
}

/// Total length of an open polyline.
pub fn polyline2d_length(pts: &[GpPnt2d]) -> f64 {
    pts.windows(2).map(|w| w[0].distance(&w[1])).sum()
}

/// Project `p` onto segment `a-b`. Returns the closest point and its
/// parameter `t` in `[0, 1]`.
pub fn project_point_on_segment(p: &GpPnt2d, a: &GpPnt2d, b: &GpPnt2d) -> (GpPnt2d, f64) {
    let abx = b.x() - a.x();
    let aby = b.y() - a.y();
    let len2 = abx * abx + aby * aby;
    let t = if len2 < 1e-30 {
        0.0
    } else {
        (((p.x() - a.x()) * abx + (p.y() - a.y()) * aby) / len2).clamp(0.0, 1.0)
    };
    (GpPnt2d::new(a.x() + t * abx, a.y() + t * aby), t)
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
fn refine_closest(c: &dyn Curve2d, p: &GpPnt2d, a: f64, b: f64) -> Option<(f64, GpPnt2d)> {
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

/// Approximate bounding box of a bounded curve by sampling: `(xmin, xmax, ymin, ymax)`.
fn curve_bbox(c: &dyn Curve2d, n: usize) -> Option<(f64, f64, f64, f64)> {
    let a = c.first_parameter();
    let b = c.last_parameter();
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
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
    Some((minx, maxx, miny, maxy))
}

/// Sample a curve over a finite parameter window. Unbounded curves get a
/// window that comfortably covers the other curve's bounding box.
fn sample_curve(c: &dyn Curve2d, other: &dyn Curve2d, n: usize) -> Vec<(f64, GpPnt2d)> {
    let a = c.first_parameter();
    let b = c.last_parameter();
    let (lo, hi) = if a.is_finite() && b.is_finite() {
        (a, b)
    } else {
        let span = match curve_bbox(other, 32) {
            Some((x0, x1, y0, y1)) => {
                let size = (x1 - x0).max(y1 - y0).max(1e-6);
                let speed = c.d1(0.0).1.magnitude().max(1e-30);
                (4.0 * size / speed).max(1.0)
            }
            None => 100.0,
        };
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
    use crate::circle::Geom2dCircle;
    use crate::line::Geom2dLine;
    use occt_core::gp::{GpAx22d, GpCirc2d, GpDir2d, GpVec2d};

    #[test]
    fn circle_length() {
        let circle = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let len = curve2d_length(&circle, 100);
        assert!((len - 2.0 * std::f64::consts::PI).abs() < 1e-6, "len={len}");
    }

    #[test]
    fn circle_line_two_intersections() {
        let circle = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let line = Geom2dLine::from_pnt_dir(GpPnt2d::zero(), GpDir2d::new(1.0, 0.0).unwrap());
        let hits = curve2d_intersections(&circle as &dyn Curve2d, &line as &dyn Curve2d, 1e-6);
        assert_eq!(hits.len(), 2, "hits={hits:?}");
        let mut xs: Vec<f64> = hits.iter().map(|(_, _, p)| p.x()).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((xs[0] + 1.0).abs() < 1e-4, "x0={}", xs[0]);
        assert!((xs[1] - 1.0).abs() < 1e-4, "x1={}", xs[1]);
    }

    #[test]
    fn crossing_segments() {
        let a = GpPnt2d::new(0.0, 0.0);
        let b = GpPnt2d::new(2.0, 2.0);
        let c = GpPnt2d::new(0.0, 2.0);
        let d = GpPnt2d::new(2.0, 0.0);
        let p = segment_intersection(&a, &b, &c, &d).unwrap();
        assert!((p.x() - 1.0).abs() < 1e-12 && (p.y() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn parallel_lines_none() {
        let l1 = GpLin2d::from_pnt_dir(GpPnt2d::zero(), GpDir2d::new(1.0, 0.0).unwrap());
        let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 1.0), GpDir2d::new(1.0, 0.0).unwrap());
        assert!(lin2d_intersection(&l1, &l2).is_none());
    }

    #[test]
    fn line_intersection_point() {
        let l1 = GpLin2d::from_pnt_dir(GpPnt2d::zero(), GpDir2d::new(1.0, 0.0).unwrap());
        let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, 1.0), GpDir2d::new(0.0, 1.0).unwrap());
        let p = lin2d_intersection(&l1, &l2).unwrap();
        assert!((p.x() - 0.0).abs() < 1e-12 && (p.y() - 0.0).abs() < 1e-12);
    }

    #[test]
    fn closest_point_on_line_is_foot() {
        let line = Geom2dLine::from_pnt_dir(GpPnt2d::zero(), GpDir2d::new(1.0, 0.0).unwrap());
        let p = GpPnt2d::new(3.0, 4.0);
        let (u, foot) = curve2d_closest_point(&line as &dyn Curve2d, &p, 1e-6).unwrap();
        assert!((u - 3.0).abs() < 1e-6, "u={u}");
        assert!((foot.x() - 3.0).abs() < 1e-6 && (foot.y() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn project_on_segment_clamps() {
        let p = GpPnt2d::new(3.0, 4.0);
        let a = GpPnt2d::new(0.0, 0.0);
        let b = GpPnt2d::new(10.0, 0.0);
        let (foot, t) = project_point_on_segment(&p, &a, &b);
        assert!((t - 0.3).abs() < 1e-12, "t={t}");
        assert!((foot.x() - 3.0).abs() < 1e-12 && (foot.y() - 0.0).abs() < 1e-12);
        // Point beyond the segment clamps to the end.
        let (foot2, t2) = project_point_on_segment(&GpPnt2d::new(20.0, 0.0), &a, &b);
        assert!((t2 - 1.0).abs() < 1e-12);
        assert!(foot2.distance(&b) < 1e-12);
    }

    #[test]
    fn polyline_length() {
        let pts = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(3.0, 4.0),
            GpPnt2d::new(3.0, 6.0),
        ];
        assert!((polyline2d_length(&pts) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn distance_to_point_on_circle() {
        let circle = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let d = curve2d_distance_to_point(&circle as &dyn Curve2d, &GpPnt2d::new(3.0, 0.0), 1e-6);
        assert!((d - 2.0).abs() < 1e-6, "d={d}");
    }

    #[test]
    fn closest_point_on_circle() {
        let circle = Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), 1.0));
        let (_, q) = curve2d_closest_point(&circle as &dyn Curve2d, &GpPnt2d::new(3.0, 0.0), 1e-6).unwrap();
        assert!((q.x() - 1.0).abs() < 1e-6 && (q.y() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn polyline_length_empty() {
        assert_eq!(polyline2d_length(&[]), 0.0);
        assert_eq!(polyline2d_length(&[GpPnt2d::zero()]), 0.0);
    }
}
