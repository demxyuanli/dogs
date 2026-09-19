//! 2D offset curves and polygon offsetting.
//!
//! `OffsetCurve2d` mirrors OCCT's `Geom2d_OffsetCurve`; the polygon helpers
//! implement 2D parallel polygons (offsetting a polygon boundary outward or
//! inward by a constant distance), including a self-intersection-corrected
//! variant and Minkowski-sum helpers.

use std::sync::Arc;
use occt_core::geom::polygon_ops::{convex_hull2d, polygon_area2d};
use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d, GpXY};
use crate::curve::Curve2d;
use crate::geom2d_api::normal2d;
use crate::HandleCurve2d;

/// A curve at a constant normal distance from a basis curve.
///
/// `d0(u) = basis.d0(u) + offset * direction * normal2d(basis, u)`, where
/// `normal2d` is the left unit normal (tangent rotated +90°). `direction` is a
/// ±1 sign convention: `+1` (the default) offsets along that left normal, `-1`
/// along the opposite (right) normal. For a counter-clockwise circle the left
/// normal points toward the centre, so a positive offset shrinks its radius;
/// offset by a negative amount (or set `direction = -1`) to grow it.
///
/// Reversing the curve reverses the basis and flips `direction`.
#[derive(Clone)]
pub struct OffsetCurve2d {
    pub basis: HandleCurve2d,
    pub offset: f64,
    pub direction: f64,
}

impl OffsetCurve2d {
    pub fn new(basis: HandleCurve2d, offset: f64) -> Self {
        Self { basis, offset, direction: 1.0 }
    }
    pub fn with_direction(basis: HandleCurve2d, offset: f64, direction: f64) -> Self {
        Self { basis, offset, direction }
    }

    /// Derivative of the unit normal w.r.t. the parameter, by central finite
    /// differences. The step scales with `|u|` to keep relative accuracy
    /// across parameter ranges. ponytail: FD is fine for tooling; swap in the
    /// analytic `d1 · d2`-based formula if exact derivatives are ever needed.
    fn normal_derivative(&self, u: f64) -> GpVec2d {
        let h = 1e-6 * (1.0 + u.abs());
        let np = normal2d(&*self.basis, u + h);
        let nm = normal2d(&*self.basis, u - h);
        GpVec2d::new((np.x() - nm.x()) / (2.0 * h), (np.y() - nm.y()) / (2.0 * h))
    }

    fn normal_second_derivative(&self, u: f64) -> GpVec2d {
        let h = 1e-6 * (1.0 + u.abs());
        let np = normal2d(&*self.basis, u + h);
        let n0 = normal2d(&*self.basis, u);
        let nm = normal2d(&*self.basis, u - h);
        GpVec2d::new(
            (np.x() - 2.0 * n0.x() + nm.x()) / (h * h),
            (np.y() - 2.0 * n0.y() + nm.y()) / (h * h),
        )
    }
}

impl Curve2d for OffsetCurve2d {
    fn d0(&self, u: f64) -> GpPnt2d {
        let p = self.basis.d0(u);
        let n = normal2d(&*self.basis, u);
        GpPnt2d::new(
            p.x() + self.offset * self.direction * n.x(),
            p.y() + self.offset * self.direction * n.y(),
        )
    }
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let (_, bd1) = self.basis.d1(u);
        let dn = self.normal_derivative(u);
        let v = GpVec2d::new(
            bd1.x() + self.offset * self.direction * dn.x(),
            bd1.y() + self.offset * self.direction * dn.y(),
        );
        (self.d0(u), v)
    }
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (_, _, bd2) = self.basis.d2(u);
        let d2n = self.normal_second_derivative(u);
        let (p, v1) = self.d1(u);
        let v2 = GpVec2d::new(
            bd2.x() + self.offset * self.direction * d2n.x(),
            bd2.y() + self.offset * self.direction * d2n.y(),
        );
        (p, v1, v2)
    }
    fn first_parameter(&self) -> f64 {
        self.basis.first_parameter()
    }
    fn last_parameter(&self) -> f64 {
        self.basis.last_parameter()
    }
    fn is_periodic(&self) -> bool {
        self.basis.is_periodic()
    }
    fn period(&self) -> f64 {
        self.basis.period()
    }
    fn continuity(&self) -> u8 {
        self.basis.continuity()
    }
    fn transform(&mut self, t: &GpTrsf2d) {
        // Offset distance and sign convention are invariant; only the basis
        // geometry moves.
        let mut c = self.basis.clone_dyn();
        c.transform(t);
        self.basis = Arc::from(c);
    }
    fn reverse(&mut self) {
        let mut c = self.basis.clone_dyn();
        c.reverse();
        self.basis = Arc::from(c);
        self.direction = -self.direction;
    }
    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.clone())
    }

    /// `Geom2d_OffsetCurve::BasisCurve()` (`Geom2d_OffsetCurve.cxx:174-177`).
    fn offset_basis(&self) -> Option<&dyn Curve2d> {
        Some(&*self.basis)
    }
}

/// Wrap a basis curve in an [`OffsetCurve2d`] with the default sign convention
/// (`direction = +1`, offset along the left normal).
pub fn offset_curve(basis: &dyn Curve2d, offset: f64) -> HandleCurve2d {
    Arc::new(OffsetCurve2d {
        basis: Arc::from(basis.clone_dyn()),
        offset,
        direction: 1.0,
    })
}

/// Offset a closed polygon by `offset`. Positive offset moves every edge
/// outward (away from the polygon interior); negative moves inward. For each
/// original vertex the two incident offset edges are intersected to give the
/// new vertex (the exact "parallel polygon"). Degenerate intersections
/// (parallel/coincident incident edges) fall back to the midpoint of the two
/// offset edge points. Self-intersections that arise from large offsets are
/// NOT resolved — use [`offset_polygon_corrected`] for that.
pub fn offset_polygon(points: &[GpPnt2d], offset: f64) -> Result<Vec<GpPnt2d>, String> {
    let n = points.len();
    if n < 3 {
        return Err("offset_polygon: need at least 3 points".into());
    }
    // Winding sign: +1 for CCW, -1 for CW. The outward unit normal of an edge
    // (e.x, e.y) is `sign * (e.y, -e.x) / |e|`.
    let sign = if polygon_area2d(points) >= 0.0 { 1.0 } else { -1.0 };

    // Offset line of each edge: direction e, through a + offset * outward normal.
    let mut edges: Vec<(GpXY, GpPnt2d)> = Vec::with_capacity(n);
    for i in 0..n {
        let a = points[i];
        let b = points[(i + 1) % n];
        let e = GpXY::new(b.x() - a.x(), b.y() - a.y());
        let len = e.modulus();
        if len < 1e-30 {
            return Err("offset_polygon: degenerate edge".into());
        }
        let q = GpPnt2d::new(
            a.x() + offset * sign * e.y / len,
            a.y() - offset * sign * e.x / len,
        );
        edges.push((e, q));
    }

    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + n - 1) % n;
        let (ej, qj) = edges[j];
        let (ei, qi) = edges[i];
        let denom = ej.crossed(&ei);
        let p = if denom.abs() < 1e-15 {
            // Parallel/coincident incident edges: lerp between the offset points.
            GpPnt2d::new((qj.x() + qi.x()) * 0.5, (qj.y() + qi.y()) * 0.5)
        } else {
            let w = qi.xy().subtracted(qj.xy());
            let t = w.crossed(&ei) / denom;
            GpPnt2d::new(qj.x() + t * ej.x, qj.y() + t * ej.y)
        };
        out.push(p);
    }
    Ok(out)
}

/// Like [`offset_polygon`], but removes the inverted loops that large offsets
/// (especially inward offsets) introduce. After computing the raw offset it
/// walks the polygon and drops any vertex whose local turn sign has flipped
/// relative to the corresponding corner of the original polygon, iterating
/// until stable. Robust for convex and mildly concave inputs.
pub fn offset_polygon_corrected(points: &[GpPnt2d], offset: f64) -> Result<Vec<GpPnt2d>, String> {
    let raw = offset_polygon(points, offset)?;
    if raw.len() < 3 {
        return Ok(raw);
    }
    let n = points.len();
    let orig_turn: Vec<f64> = (0..n)
        .map(|i| turn_sign(&points[(i + n - 1) % n], &points[i], &points[(i + 1) % n]))
        .collect();

    // Ring of (original corner index, offset vertex) — the index survives
    // removals so later passes still compare against the right original corner.
    let mut ring: Vec<(usize, GpPnt2d)> = raw.iter().cloned().enumerate().collect();
    loop {
        let m = ring.len();
        if m < 3 {
            break;
        }
        let mut remove = vec![false; m];
        let mut changed = false;
        for i in 0..m {
            let (ci, p) = ring[i];
            let (_, a) = ring[(i + m - 1) % m];
            let (_, c) = ring[(i + 1) % m];
            let o = orig_turn[ci % n];
            if o != 0.0 && turn_sign(&a, &p, &c) != o {
                remove[i] = true;
                changed = true;
            }
        }
        if !changed {
            break;
        }
        ring = ring
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !remove[*i])
            .map(|(_, v)| v)
            .collect();
    }
    let out: Vec<GpPnt2d> = ring.into_iter().map(|(_, p)| p).collect();
    Ok(if out.len() < 3 { Vec::new() } else { out })
}

/// Sign of the turn (cross product of the two incident edges) at vertex `b`:
/// +1 left turn, -1 right turn, 0 collinear.
fn turn_sign(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d) -> f64 {
    let cross = (b.x() - a.x()) * (c.y() - b.y()) - (b.y() - a.y()) * (c.x() - b.x());
    if cross.abs() < 1e-12 {
        0.0
    } else if cross > 0.0 {
        1.0
    } else {
        -1.0
    }
}

/// Signed area of the offset polygon. Returns 0 when the offset is degenerate.
pub fn polygon_offset_area(points: &[GpPnt2d], offset: f64) -> f64 {
    match offset_polygon(points, offset) {
        Ok(p) => polygon_area2d(&p),
        Err(_) => 0.0,
    }
}

/// Minkowski sum of two convex polygons as the convex hull of all pairwise
/// point sums.
pub fn minkowski_sum_convex(a: &[GpPnt2d], b: &[GpPnt2d]) -> Vec<GpPnt2d> {
    let mut sums = Vec::with_capacity(a.len() * b.len());
    for &pa in a {
        for &pb in b {
            sums.push(GpPnt2d::new(pa.x() + pb.x(), pa.y() + pb.y()));
        }
    }
    convex_hull2d(&sums)
}

/// Approximate the Minkowski sum of a point set with a disc of `radius` as a
/// rounded polygon: the convex hull of the input offset by `radius` along each
/// hull edge's outward normal (K samples per edge), with the corners rounded
/// by sampling the arc of radius `radius` centred on each hull vertex between
/// the two incident edge normal directions.
pub fn inflate_point_set(points: &[GpPnt2d], radius: f64) -> Vec<GpPnt2d> {
    let hull = convex_hull2d(points);
    let n = hull.len();
    if n < 3 || radius <= 0.0 {
        return hull;
    }
    const K: usize = 8;
    let mut out: Vec<GpPnt2d> = Vec::with_capacity(n * (K + 1) * 2);
    // K samples per hull edge, offset by `radius` along the outward normal.
    for i in 0..n {
        let a = hull[i];
        let b = hull[(i + 1) % n];
        let e = GpXY::new(b.x() - a.x(), b.y() - a.y());
        let len = e.modulus();
        if len < 1e-30 {
            continue;
        }
        let (nx, ny) = (e.y / len, -e.x / len);
        for k in 0..=K {
            let t = k as f64 / K as f64;
            out.push(GpPnt2d::new(
                a.x() + t * e.x + radius * nx,
                a.y() + t * e.y + radius * ny,
            ));
        }
    }
    // Arc rounding: interpolate the outward-normal direction between the two
    // incident edges at each hull vertex (the normals rotate counter-clockwise
    // around a CCW convex vertex).
    let tau = 2.0 * std::f64::consts::PI;
    for i in 0..n {
        let v = hull[i];
        let prev = hull[(i + n - 1) % n];
        let next = hull[(i + 1) % n];
        let e1 = GpXY::new(v.x() - prev.x(), v.y() - prev.y());
        let e2 = GpXY::new(next.x() - v.x(), next.y() - v.y());
        let l1 = e1.modulus();
        let l2 = e2.modulus();
        if l1 < 1e-30 || l2 < 1e-30 {
            continue;
        }
        let a1 = (-e1.x / l1).atan2(e1.y / l1);
        let a2 = (-e2.x / l2).atan2(e2.y / l2);
        let mut da = a2 - a1;
        while da < 0.0 {
            da += tau;
        }
        if da > tau * 0.999 {
            continue;
        }
        for k in 1..K {
            let ang = a1 + da * k as f64 / K as f64;
            out.push(GpPnt2d::new(v.x() + radius * ang.cos(), v.y() + radius * ang.sin()));
        }
    }
    convex_hull2d(&out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circle::Geom2dCircle;
    use crate::curve_ops::segment_intersection;
    use crate::line::Geom2dLine;
    use occt_core::gp::{GpAx22d, GpCirc2d, GpDir2d};

    fn p2(x: f64, y: f64) -> GpPnt2d {
        GpPnt2d::new(x, y)
    }

    fn square() -> Vec<GpPnt2d> {
        vec![p2(0., 0.), p2(1., 0.), p2(1., 1.), p2(0., 1.)]
    }

    /// True when no two non-adjacent edges of the closed polygon cross.
    fn is_simple_polygon(pts: &[GpPnt2d]) -> bool {
        let n = pts.len();
        if n < 3 {
            return false;
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if j == i + 1 || (i == 0 && j == n - 1) {
                    continue;
                }
                if segment_intersection(
                    &pts[i],
                    &pts[(i + 1) % n],
                    &pts[j],
                    &pts[(j + 1) % n],
                )
                .is_some()
                {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn offset_square_positive_area() {
        let area = polygon_offset_area(&square(), 0.1);
        assert!((area - 1.44).abs() < 1e-6, "area={area}");
    }

    #[test]
    fn offset_square_negative_area() {
        let area = polygon_offset_area(&square(), -0.1);
        assert!((area - 0.64).abs() < 1e-6, "area={area}");
    }

    #[test]
    fn offset_curve_line_is_parallel() {
        let line = Geom2dLine::from_pnt_dir(p2(0., 0.), GpDir2d::new(1.0, 0.0).unwrap());
        let off = offset_curve(&line as &dyn Curve2d, 0.5);
        let p = off.d0(2.0);
        // Parallel line at distance |offset| (normal of the horizontal line is
        // vertical, so the y-coordinate is exactly the offset).
        assert!((p.y() - 0.5).abs() < 1e-12, "p={p:?}");
        // Constant normal ⇒ derivative of the normal is 0, so d1 == basis d1.
        let (_, v) = off.d1(0.0);
        assert!(v.y().abs() < 1e-9 && (v.x() - 1.0).abs() < 1e-9, "v={v:?}");
    }

    #[test]
    fn offset_curve_circle_radius() {
        let ax = GpAx22d::from_xdir(p2(0., 0.), GpDir2d::default());
        let circ = Geom2dCircle::new(GpCirc2d::new(ax, 1.0));
        // normal2d of a CCW circle points inward, so the offset curve is a
        // concentric circle whose radius differs from r by |offset|.
        let off = offset_curve(&circ as &dyn Curve2d, 0.3);
        let r = off.d0(0.0).distance(&p2(0., 0.));
        assert!((r - 0.7).abs() < 1e-9, "r={r}");
        let off2 = offset_curve(&circ as &dyn Curve2d, -0.3);
        let r2 = off2.d0(0.0).distance(&p2(0., 0.));
        assert!((r2 - 1.3).abs() < 1e-9, "r2={r2}");
    }

    #[test]
    fn offset_polygon_triangle_larger() {
        let tri = vec![p2(0., 0.), p2(2., 0.), p2(0., 2.)];
        let off = offset_polygon(&tri, 0.2).unwrap();
        assert_eq!(off.len(), 3);
        let a = polygon_area2d(&tri);
        let ao = polygon_area2d(&off);
        assert!(ao > a, "ao={ao} a={a}");
        // The (0,0) corner has a 90° interior angle, so the new vertex lies at
        // offset / sin(45°) = offset * sqrt(2) from the corner.
        let d = off[0].distance(&tri[0]);
        assert!((d - 0.2 * std::f64::consts::SQRT_2).abs() < 1e-9, "d={d}");
    }

    #[test]
    fn minkowski_sum_two_squares() {
        let sq = square();
        let m = minkowski_sum_convex(&sq, &sq);
        assert_eq!(m.len(), 4);
        assert!((polygon_area2d(&m) - 4.0).abs() < 1e-9, "area={}", polygon_area2d(&m));
        let mut sides: Vec<f64> = (0..m.len()).map(|i| m[i].distance(&m[(i + 1) % m.len()])).collect();
        sides.sort_by(|a, b| a.total_cmp(b));
        assert!((sides[0] - 2.0).abs() < 1e-9, "sides={sides:?}");
    }

    #[test]
    fn corrected_concave_polygon_is_simple() {
        // CCW house-with-notch: (1.5,1) is a concave dent. A 0.3 outward offset
        // flips the turn at the two adjacent corners in the raw offset; the
        // corrected walk drops them and bridges the notch.
        let l = vec![p2(0., 0.), p2(3., 0.), p2(3., 3.), p2(1.5, 1.), p2(0., 3.)];
        let corrected = offset_polygon_corrected(&l, 0.3).unwrap();
        assert!(!corrected.is_empty(), "corrected empty");
        let a = polygon_area2d(&corrected);
        assert!(a.is_finite() && a > 0.0, "area={a}");
        assert!(is_simple_polygon(&corrected), "corrected self-intersects: {corrected:?}");
    }

    #[test]
    fn inflate_point_set_rounds() {
        let inflated = inflate_point_set(&square(), 0.5);
        assert!(inflated.len() >= 8, "len={}", inflated.len());
        let a = polygon_area2d(&inflated);
        // Rounded 1x1 square with radius 0.5: above the chamfered polygon (3.5),
        // approaching the true Minkowski area (1 + 4*0.5 + pi*0.25 ~ 3.785).
        assert!(a > 3.5 && a < 3.8, "area={a}");
    }
}
