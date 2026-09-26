//! 2D curve construction and analysis tools.
//!
//! Port of `Geom2dAPI_PointsToBSpline` (interpolation), `GC_MakeArcOfCircle2d`
//! (arc through points), `GCPnts`-style uniform points, and 2D offset/round
//! helpers. Source: `Geom2dAPI`, `GC` (TKGeomBase).

use occt_core::gp::{GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpPnt2d, GpVec2d};

use crate::curve::Curve2d;
use crate::{Geom2dBSplineCurve, Geom2dCircle, Geom2dLine};

/// Arc of a 2D circle passing through three points (GC_MakeArcOfCircle2d).
/// Returns (center, radius, start_angle, end_angle).
pub fn arc_from_three_points(p1: GpPnt2d, p2: GpPnt2d, p3: GpPnt2d) -> Option<(GpPnt2d, f64, f64, f64)> {
    let (cx, cy) = circumcenter(p1, p2, p3)?;
    let center = GpPnt2d::new(cx, cy);
    let r = center.distance(&p1);
    if r < 1e-12 {
        return None;
    }
    // P2 enters only through the circumcentre: OCCT's
    // `GC_MakeArcOfCircle(P1, P2, P3)` (`GC_MakeArcOfCircle.cxx:33-61`) sets
    // `Alpha1 = 0.` and `Alpha3 = ElCLib::Parameter(C, P3)` with `sense = true`
    // and leaves the `Alpha2` lines **commented out** (`cxx:41-51`).
    let a1 = (p1.y() - cy).atan2(p1.x() - cx);
    let a3 = (p3.y() - cy).atan2(p3.x() - cx);
    Some((center, r, a1, a3))
}

/// Circumcenter of three 2D points.
pub fn circumcenter(p1: GpPnt2d, p2: GpPnt2d, p3: GpPnt2d) -> Option<(f64, f64)> {
    let d = 2.0 * (p1.x() * (p2.y() - p3.y()) + p2.x() * (p3.y() - p1.y()) + p3.x() * (p1.y() - p2.y()));
    if d.abs() < 1e-12 {
        return None; // collinear
    }
    let x1 = p1.x() * p1.x() + p1.y() * p1.y();
    let x2 = p2.x() * p2.x() + p2.y() * p2.y();
    let x3 = p3.x() * p3.x() + p3.y() * p3.y();
    let ux = (x1 * (p2.y() - p3.y()) + x2 * (p3.y() - p1.y()) + x3 * (p1.y() - p2.y())) / d;
    let uy = (x1 * (p3.x() - p2.x()) + x2 * (p1.x() - p3.x()) + x3 * (p2.x() - p1.x())) / d;
    Some((ux, uy))
}

/// Interpolate a 2D B-spline through `points` (Geom2dAPI_PointsToBSpline).
/// Degree 1 → polyline through points; degree ≥ 2 → global interpolation with
/// clamped uniform knots.
pub fn interpolate_points_2d(points: &[GpPnt2d], degree: usize) -> Result<Geom2dBSplineCurve, String> {
    if points.len() < 2 {
        return Err("interpolate_points_2d: need at least 2 points".into());
    }
    if degree == 0 {
        return Err("interpolate_points_2d: degree must be ≥ 1".into());
    }
    if degree >= points.len() {
        return Err("interpolate_points_2d: degree must be < point count".into());
    }
    // Chord-length parameterization.
    let n = points.len();
    let mut param = vec![0.0f64; n];
    for i in 1..n {
        param[i] = param[i - 1] + points[i - 1].distance(&points[i]);
    }
    let total = param[n - 1];
    if total < 1e-12 {
        return Err("interpolate_points_2d: duplicate points".into());
    }
    for p in param.iter_mut() {
        *p /= total;
    }
    // Clamped uniform knots: poles = n, knots = n + degree + 1.
    let m = n + degree + 1;
    let mut knots = vec![0.0f64; m];
    for i in 0..=degree {
        knots[i] = 0.0;
        knots[m - 1 - i] = 1.0;
    }
    let interior = m - 2 * (degree + 1);
    for j in 0..interior {
        knots[degree + 1 + j] = (j + 1) as f64 / (interior + 1) as f64;
    }
    if degree == 1 {
        // Control points = data points (piecewise linear).
        let xs: Vec<f64> = points.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = points.iter().map(|p| p.y()).collect();
        return Geom2dBSplineCurve::new(xs, ys, knots, 1).map_err(|e| e.to_string());
    }
    // Global interpolation: build collocation matrix B (n×n) of B-spline
    // basis at each parameter; solve B·poles = points per coordinate.
    let (poles_x, poles_y) = solve_collocation(&param, &knots, degree, points)?;
    Geom2dBSplineCurve::new(poles_x, poles_y, knots, degree).map_err(|e| e.to_string())
}

/// 1-D B-spline basis value N_i^p(u) with clamped uniform knots.
fn basis_value(knots: &[f64], degree: usize, i: usize, u: f64) -> f64 {
    let n = knots.len() - degree - 1;
    if i >= n {
        return 0.0;
    }
    if u < knots[0] || u > knots[n] {
        return 0.0;
    }
    if degree == 0 {
        return if u >= knots[i] && u < knots[i + 1] { 1.0 } else { 0.0 };
    }
    // Cox–de Boor recursion: N_i^p uses only knots[i ..= i+p+1], so keep the
    // local window of the basis for pole i. At the right boundary u == knots[n]
    // (the last clamped knot), the final non-zero basis function must equal 1.
    let n_last = knots.len() - degree - 1;
    let at_last = u >= knots[n_last] && (u - knots[n_last]).abs() < 1e-12;
    let mut n0 = vec![0.0f64; degree + 1];
    for j in 0..=degree {
        let hi = if at_last && i + j + 1 == n_last {
            knots[i + j + 1] + 1.0 // include the boundary
        } else {
            knots[i + j + 1]
        };
        n0[j] = if u >= knots[i + j] && u < hi { 1.0 } else { 0.0 };
    }
    let mut cur = n0;
    for p in 1..=degree {
        let mut next = vec![0.0f64; degree + 1 - p];
        for j in 0..(degree + 1 - p) {
            let den1 = knots[i + j + p] - knots[i + j];
            let den2 = knots[i + j + p + 1] - knots[i + j + 1];
            let mut a = 0.0;
            if den1.abs() > 1e-12 {
                a = (u - knots[i + j]) / den1 * cur[j];
            }
            let mut b = 0.0;
            if den2.abs() > 1e-12 {
                b = (knots[i + j + p + 1] - u) / den2 * cur[j + 1];
            }
            next[j] = a + b;
        }
        cur = next;
    }
    cur[0]
}

fn solve_collocation(
    param: &[f64],
    knots: &[f64],
    degree: usize,
    points: &[GpPnt2d],
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let n = points.len();
    // B[i][j] = N_j^p(param[i]) — basis of pole j at sample i.
    let mut b = vec![vec![0.0f64; n]; n];
    for (i, &u) in param.iter().enumerate() {
        for j in 0..n {
            b[i][j] = basis_value(knots, degree, j, u);
        }
    }
    let mut rhs_x = vec![0.0f64; n];
    let mut rhs_y = vec![0.0f64; n];
    for (i, p) in points.iter().enumerate() {
        rhs_x[i] = p.x();
        rhs_y[i] = p.y();
    }
    // Gaussian elimination with partial pivoting.
    let x = solve_gauss(&b, &rhs_x)?;
    let y = solve_gauss(&b, &rhs_y)?;
    Ok((x, y))
}

fn solve_gauss(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, String> {
    let n = b.len();
    let mut m = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        // Partial pivot.
        let mut piv = col;
        for r in (col + 1)..n {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-12 {
            return Err("solve_gauss: singular collocation matrix".into());
        }
        m.swap(col, piv);
        rhs.swap(col, piv);
        let d = m[col][col];
        for c in col..n {
            m[col][c] /= d;
        }
        rhs[col] /= d;
        for r in 0..n {
            if r != col {
                let f = m[r][col];
                if f.abs() > 1e-15 {
                    for c in col..n {
                        m[r][c] -= f * m[col][c];
                    }
                    rhs[r] -= f * rhs[col];
                }
            }
        }
    }
    Ok(rhs)
}

/// Uniformly spaced points along a curve (GCPnts_UniformAbscissa-lite).
pub fn uniform_abscissa_points(c: &dyn Curve2d, n: usize) -> Vec<(f64, GpPnt2d)> {
    if n < 2 {
        return vec![(c.first_parameter(), c.d0(c.first_parameter()))];
    }
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let (a, b) = if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    // Sample arc length via chord lengths on a fine grid.
    let samples = 256;
    let mut cum = vec![0.0f64; samples + 1];
    let mut pts = vec![GpPnt2d::zero(); samples + 1];
    for i in 0..=samples {
        let u = a + (b - a) * i as f64 / samples as f64;
        pts[i] = c.d0(u);
        if i > 0 {
            cum[i] = cum[i - 1] + pts[i - 1].distance(&pts[i]);
        }
    }
    let total = cum[samples];
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let target = total * k as f64 / (n - 1) as f64;
        let mut idx = 0;
        while idx < samples && cum[idx + 1] < target {
            idx += 1;
        }
        let u = a + (b - a) * idx as f64 / samples as f64;
        out.push((u, c.d0(u)));
    }
    out
}

/// Offset a point set outward (inflate a convex 2D shape boundary) — helper
/// for round-corner construction. Positive `dist` moves each point away from
/// the centroid.
pub fn inflate_points(points: &[GpPnt2d], dist: f64) -> Vec<GpPnt2d> {
    if points.is_empty() {
        return Vec::new();
    }
    let cx = points.iter().map(|p| p.x()).sum::<f64>() / points.len() as f64;
    let cy = points.iter().map(|p| p.y()).sum::<f64>() / points.len() as f64;
    points
        .iter()
        .map(|p| {
            let dx = p.x() - cx;
            let dy = p.y() - cy;
            let m = (dx * dx + dy * dy).sqrt();
            if m < 1e-12 {
                *p
            } else {
                GpPnt2d::new(p.x() + dx / m * dist, p.y() + dy / m * dist)
            }
        })
        .collect()
}

/// Round a sharp 2D corner (three points forming a vertex) into a circular
/// arc fillet of radius `r` (2D fillet). Returns the fillet arc curve.
pub fn fillet_corner_2d(p1: GpPnt2d, vertex: GpPnt2d, p3: GpPnt2d, r: f64) -> Option<Geom2dCircle> {
    if r <= 0.0 {
        return None;
    }
    let d1 = GpVec2d::new(p1.x() - vertex.x(), p1.y() - vertex.y());
    let d2 = GpVec2d::new(p3.x() - vertex.x(), p3.y() - vertex.y());
    let m1 = d1.magnitude();
    let m2 = d2.magnitude();
    if m1 < 1e-12 || m2 < 1e-12 {
        return None;
    }
    let u1 = GpVec2d::new(d1.x() / m1, d1.y() / m1);
    let u2 = GpVec2d::new(d2.x() / m2, d2.y() / m2);
    let cross = u1.x() * u2.y() - u1.y() * u2.x();
    if cross.abs() < 1e-12 {
        return None; // collinear
    }
    let cos = (u1.x() * u2.x() + u1.y() * u2.y()).clamp(-1.0, 1.0);
    let half = (1.0 - cos).max(1e-12).sqrt() * 0.5; // sin(θ/2)
    let dist = r / half.max(1e-12); // distance from vertex to tangency points
    let tan1 = GpPnt2d::new(vertex.x() + u1.x() * dist, vertex.y() + u1.y() * dist);
    let tan2 = GpPnt2d::new(vertex.x() + u2.x() * dist, vertex.y() + u2.y() * dist);
    // Center = vertex + direction bisecting u1,u2 · (dist·cos(θ/2)/sin(θ/2))...
    // use the perpendicular bisector of tan1,tan2 at distance r from both.
    let mid = GpPnt2d::new((tan1.x() + tan2.x()) * 0.5, (tan1.y() + tan2.y()) * 0.5);
    let bisect = GpVec2d::new(tan2.x() - tan1.x(), tan2.y() - tan1.y());
    let bm = bisect.magnitude();
    if bm < 1e-12 {
        return None;
    }
    let perp = GpVec2d::new(-bisect.y() / bm, bisect.x() / bm);
    let half_chord = bm * 0.5;
    let center_dist = (r * r - half_chord * half_chord).max(0.0).sqrt();
    // Choose the side toward the vertex.
    let c0 = GpPnt2d::new(mid.x() + perp.x() * center_dist, mid.y() + perp.y() * center_dist);
    let c1 = GpPnt2d::new(mid.x() - perp.x() * center_dist, mid.y() - perp.y() * center_dist);
    let center = if c0.distance(&vertex) < c1.distance(&vertex) { c0 } else { c1 };
    Some(Geom2dCircle::new(GpCirc2d::new(GpAx22d::standard(), r).translated_vec(&GpVec2d::new(center.x(), center.y()))))
}

/// Line through two points (2D).
pub fn line_through(p1: GpPnt2d, p2: GpPnt2d) -> Result<Geom2dLine, String> {
    let dx = p2.x() - p1.x();
    let dy = p2.y() - p1.y();
    let m = (dx * dx + dy * dy).sqrt();
    if m < 1e-12 {
        return Err("line_through: coincident points".into());
    }
    Ok(Geom2dLine::new(GpAx2d::new(p1, GpDir2d::new(dx / m, dy / m).unwrap())))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Curve2d;

    #[test]
    fn arc_three_points_radius() {
        // Points on a unit circle centered at origin.
        let p1 = GpPnt2d::new(1.0, 0.0);
        let p2 = GpPnt2d::new(0.0, 1.0);
        let p3 = GpPnt2d::new(-1.0, 0.0);
        let (c, r, a1, a3) = arc_from_three_points(p1, p2, p3).expect("arc");
        assert!((c.x().abs() < 1e-9) && (c.y().abs() < 1e-9), "center {c:?}");
        assert!((r - 1.0).abs() < 1e-9, "radius {r}");
        assert!((a3 - a1).abs() > 3.0, "spans ~π");
    }

    #[test]
    fn arc_collinear_returns_none() {
        let r = arc_from_three_points(
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(2.0, 0.0),
        );
        assert!(r.is_none());
    }

    #[test]
    fn interp_degree1_polyline() {
        let pts = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(2.0, 0.0),
        ];
        let c = interpolate_points_2d(&pts, 1).expect("interp");
        for (i, p) in pts.iter().enumerate() {
            let u = c.first_parameter() + (c.last_parameter() - c.first_parameter()) * i as f64 / 2.0;
            let q = c.d0(u);
            assert!(p.distance(&q) < 1e-6, "point {i} {p:?} vs {q:?}");
        }
    }

    #[test]
    fn interp_degree2_passes_points() {
        let pts = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(0.5, 1.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.5, 1.0),
            GpPnt2d::new(2.0, 0.0),
        ];
        let c = interpolate_points_2d(&pts, 3).expect("interp");
        for (i, p) in pts.iter().enumerate() {
            let u = c.first_parameter() + (c.last_parameter() - c.first_parameter()) * i as f64 / 4.0;
            let q = c.d0(u);
            assert!(p.distance(&q) < 1e-5, "point {i} {p:?} vs {q:?}");
        }
    }

    #[test]
    fn interp_errors() {
        assert!(interpolate_points_2d(&[GpPnt2d::new(0.0, 0.0)], 1).is_err());
        assert!(interpolate_points_2d(&[GpPnt2d::new(0.0, 0.0), GpPnt2d::new(1.0, 0.0)], 2).is_err());
    }

    #[test]
    fn uniform_points_on_line() {
        // Use a bounded degree-1 B-spline through collinear points.
        let pts = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(2.0, 0.0),
            GpPnt2d::new(4.0, 0.0),
        ];
        let l = interpolate_points_2d(&pts, 1).expect("interp");
        let uniform = uniform_abscissa_points(&l, 5);
        assert_eq!(uniform.len(), 5);
        assert!((uniform[0].1.x() - 0.0).abs() < 1e-6);
        assert!((uniform[4].1.x() - 4.0).abs() < 0.05, "end {}", uniform[4].1.x());
        assert!((uniform[2].1.x() - 2.0).abs() < 0.05, "mid {}", uniform[2].1.x());
    }

    #[test]
    fn inflate_convex_expands() {
        let pts = vec![
            GpPnt2d::new(-1.0, -1.0),
            GpPnt2d::new(1.0, -1.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(-1.0, 1.0),
        ];
        let out = inflate_points(&pts, 0.5);
        // Corner (-1,-1) moves radially by 0.5 → (-1.354, -1.354).
        assert!(out[0].x() < -1.3, "corner moved out {}", out[0].x());
    }

    #[test]
    fn fillet_corner_2d_arc() {
        // Right angle at (1,0) between (0,0) and (1,1) — fillet r=0.3.
        let c = fillet_corner_2d(
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            0.3,
        )
        .expect("fillet");
        // The fillet circle's radius is r.
        assert!((c.circ().radius() - 0.3).abs() < 1e-9, "radius {}", c.circ().radius());
        // The tangency points are at distance r from the fillet center and the
        // arc lies inside the corner.
        let center = c.circ().location();
        assert!(center.distance(&GpPnt2d::new(1.0, 0.0)) < 1.0);
    }

    #[test]
    fn line_through_errors() {
        assert!(line_through(GpPnt2d::new(0.0, 0.0), GpPnt2d::new(0.0, 0.0)).is_err());
        let l = line_through(GpPnt2d::new(0.0, 0.0), GpPnt2d::new(3.0, 4.0)).expect("line");
        // Geom2dLine uses unit-direction parameterization: d0(1) = the unit
        // direction (0.6, 0.8); d0(0) = the origin.
        let p0 = l.d0(0.0);
        assert!((p0.x() - 0.0).abs() < 1e-9 && (p0.y() - 0.0).abs() < 1e-9, "line start {p0:?}");
        let p1 = l.d0(1.0);
        assert!((p1.x() - 0.6).abs() < 1e-9 && (p1.y() - 0.8).abs() < 1e-9, "line dir {p1:?}");
    }

    #[test]
    fn fillet_corner_collinear_none() {
        let c = fillet_corner_2d(
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(2.0, 0.0),
            0.3,
        );
        assert!(c.is_none());
    }
}
