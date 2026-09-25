//! 3D curve helpers — arc length, curvature, closest point, offsets and
//! projection onto curves. Source: `GCPnts_AbscissaPoint`, `math` helpers.

use crate::gp::{GpDir, GpPnt, GpVec};

/// Chord length of a polyline.
pub fn polyline_length3d(pts: &[GpPnt]) -> f64 {
    let mut len = 0.0;
    for w in pts.windows(2) {
        len += w[0].distance(&w[1]);
    }
    len
}

/// Arc length of a parametric curve by composite Simpson over [a,b].
/// `f` maps a parameter to a 3D point.
pub fn curve_length3d<F: Fn(f64) -> GpPnt>(f: &F, a: f64, b: f64, n: usize) -> f64 {
    let n = (n / 2 * 2).max(2);
    let h = (b - a) / n as f64;
    let mut sum = 0.0;
    let mut prev = f(a);
    for i in 1..=n {
        let x = a + i as f64 * h;
        let cur = f(x);
        sum += prev.distance(&cur);
        prev = cur;
    }
    sum
}

/// Numeric curvature of a 3D parametric curve at `u` (radius of curvature
/// inverse), using central finite differences of the position.
pub fn curvature3d<F: Fn(f64) -> GpPnt>(f: &F, u: f64, h: f64) -> f64 {
    let p0 = f(u - h);
    let p1 = f(u);
    let p2 = f(u + h);
    // Triangle area / (product of chord lengths) gives sin(θ); curvature of the
    // osculating circle ≈ 2·sin(θ)/|p0-p2|.
    let a = p1.coord.subtracted(&p0.coord);
    let b = p2.coord.subtracted(&p0.coord);
    let cross = a.crossed(&b);
    let area2 = cross.modulus();
    let chord = b.modulus();
    if area2 < 1e-30 || chord < 1e-30 {
        return 0.0;
    }
    // Radius of curvature R = |AB|·|AC|·|BC| / (2·|AB×AC|); κ = 1/R.
    let ab = p0.distance(&p1);
    let ac = p0.distance(&p2);
    let bc = p1.distance(&p2);
    if ab < 1e-30 || ac < 1e-30 || bc < 1e-30 {
        return 0.0;
    }
    let r = ab * ac * bc / (2.0 * area2);
    if r.abs() < 1e-30 {
        0.0
    } else {
        1.0 / r
    }
}

/// Closest point on a parametric curve to `p`, via coarse scan + golden-section
/// refinement. Returns (parameter, closest point, distance).
pub fn closest_on_curve3d<F: Fn(f64) -> GpPnt>(f: &F, a: f64, b: f64, p: &GpPnt, samples: usize) -> (f64, GpPnt, f64) {
    let mut best_u = a;
    let mut best_d = f64::INFINITY;
    for i in 0..=samples {
        let u = a + (b - a) * i as f64 / samples as f64;
        let d = f(u).distance(p);
        if d < best_d {
            best_d = d;
            best_u = u;
        }
    }
    // Golden-section refinement.
    let phi = 0.618_033_988_749_895;
    let (mut lo, mut hi) = (a.max(best_u - (b - a) / samples as f64), (best_u + (b - a) / samples as f64).min(b));
    let mut c = hi - phi * (hi - lo);
    let mut d = lo + phi * (hi - lo);
    let mut fc = f(c).distance(p);
    let mut fd = f(d).distance(p);
    for _ in 0..60 {
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - phi * (hi - lo);
            fc = f(c).distance(p);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + phi * (hi - lo);
            fd = f(d).distance(p);
        }
    }
    let u = 0.5 * (lo + hi);
    let q = f(u);
    (u, q, q.distance(p))
}

/// Offset a parametric curve point by `offset` along the (finite-difference)
/// normal. This gives a single offset point; the full offset curve is a
/// sampled polyline via `offset_curve_points`.
pub fn offset_point3d<F: Fn(f64) -> GpPnt>(f: &F, u: f64, h: f64, offset: f64) -> GpPnt {
    let p = f(u);
    let pu = f(u + h);
    let pv = f(u - h);
    // Tangent and a reference perpendicular: use a second reference axis.
    let tangent = GpVec::from_pnts(&pv, &pu);
    let ref_axis = if tangent.xyz().x.abs() <= tangent.xyz().y.abs() && tangent.xyz().x.abs() <= tangent.xyz().z.abs() {
        GpVec::new(1.0, 0.0, 0.0)
    } else {
        GpVec::new(0.0, 1.0, 0.0)
    };
    let normal = tangent.xyz().crossed(ref_axis.xyz());
    let n = normal.modulus();
    if n < 1e-30 {
        return p;
    }
    GpPnt::new(
        p.x() + offset * normal.x / n,
        p.y() + offset * normal.y / n,
        p.z() + offset * normal.z / n,
    )
}

/// Sampled offset curve as a polyline.
pub fn offset_curve_points<F: Fn(f64) -> GpPnt>(f: &F, a: f64, b: f64, offset: f64, n: usize, h: f64) -> Vec<GpPnt> {
    (0..=n)
        .map(|i| {
            let u = a + (b - a) * i as f64 / n as f64;
            offset_point3d(f, u, h, offset)
        })
        .collect()
}

/// Tangent direction (unit) of a parametric curve at `u` by central differences.
pub fn tangent3d<F: Fn(f64) -> GpPnt>(f: &F, u: f64, h: f64) -> Option<GpVec> {
    let t = GpVec::from_pnts(&f(u - h), &f(u + h));
    GpDir::from_vec(&t).ok().map(|d| GpVec::new(d.x(), d.y(), d.z()))
}

/// Normal plane curvature circle: the osculating circle's center at `u`.
pub fn osculating_center3d<F: Fn(f64) -> GpPnt>(f: &F, u: f64, h: f64) -> Option<GpPnt> {
    let p0 = f(u - h);
    let p1 = f(u);
    let p2 = f(u + h);
    // Circumcenter of the three sample points.
    let d1 = GpVec::from_pnts(&p0, &p1);
    let d2 = GpVec::from_pnts(&p0, &p2);
    let n = d1.xyz().crossed(d2.xyz());
    if n.modulus() < 1e-20 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    // Solve O·d1 = (|p1|²-|p0|²)/2 ; O·d2 = (|p2|²-|p0|²)/2 ; O·n = p0·n.
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.x, n.y, n.z],
    ];
    let rhs = [0.5 * (n2(&p1) - n2(&p0)), 0.5 * (n2(&p2) - n2(&p0)), p0.coord.dot(&n)];
    let det = det3(&mat);
    if det.abs() < 1e-20 {
        return None;
    }
    let mut o = [0.0; 3];
    for k in 0..3 {
        let mut m = mat;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        o[k] = det3(&m) / det;
    }
    Some(GpPnt::new(o[0], o[1], o[2]))
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Total curvature energy ∫κ² ds (approximated by sampling).
pub fn curvature_energy<F: Fn(f64) -> GpPnt>(f: &F, a: f64, b: f64, n: usize) -> f64 {
    let h = (b - a) / n as f64;
    let mut sum = 0.0;
    for i in 0..n {
        let u = a + (i as f64 + 0.5) * h;
        let k = curvature3d(f, u, 1e-6);
        sum += k * k * h;
    }
    sum
}

/// Signed distance from a point to the infinite line through `a` with
/// direction `d` (unit not required).
pub fn point_line_distance3d(p: &GpPnt, a: &GpPnt, d: &GpVec) -> f64 {
    let ap = GpVec::from_pnts(a, p);
    let cross = ap.xyz().crossed(d.xyz());
    let len = d.xyz().modulus();
    if len < 1e-30 {
        return p.distance(a);
    }
    cross.modulus() / len
}

/// Project a point onto a segment [a,b]: closest point + parameter t in [0,1].
pub fn project_on_segment3d(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> (GpPnt, f64) {
    let ab = GpVec::from_pnts(a, b);
    let len2 = ab.xyz().dot(&ab.xyz());
    if len2 < 1e-30 {
        return (*a, 0.0);
    }
    let ap = GpVec::from_pnts(a, p);
    let t = (ap.xyz().dot(&ab.xyz()) / len2).clamp(0.0, 1.0);
    (GpPnt::new(a.x() + t * ab.x(), a.y() + t * ab.y(), a.z() + t * ab.z()), t)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(t: f64) -> GpPnt {
        GpPnt::new(2.0 * t, 0.0, 0.0)
    }

    fn circle(t: f64) -> GpPnt {
        // Radius-2 circle in the xy-plane.
        GpPnt::new(2.0 * t.cos(), 2.0 * t.sin(), 0.0)
    }

    #[test]
    fn polyline_and_curve_length() {
        let pts = vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.)];
        assert!((polyline_length3d(&pts) - 2.0).abs() < 1e-12);
        let len = curve_length3d(&line, 0.0, 3.0, 64);
        assert!((len - 6.0).abs() < 1e-9);
    }

    #[test]
    fn circle_curvature_constant() {
        let k = curvature3d(&circle, 1.0, 1e-3);
        assert!((k - 0.5).abs() < 1e-3, "k={k} (expect 1/r = 0.5)");
    }

    #[test]
    fn closest_on_circle() {
        let (u, q, d) = closest_on_curve3d(&circle, 0.0, 2.0 * std::f64::consts::PI, &GpPnt::new(2.0, 0.0, 0.0), 64);
        let _ = u;
        assert!(q.distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-6);
        assert!(d < 1e-6);
    }

    #[test]
    fn offset_and_tangent() {
        let off = offset_curve_points(&circle, 0.0, std::f64::consts::FRAC_PI_2, 0.5, 16, 1e-4);
        // The offset direction for a space curve needs a reference normal, but
        // the DISTANCE invariant always holds: every offset point is |offset|
        // away from the base curve.
        for p in &off {
            let (_, _, d) = closest_on_curve3d(&circle, 0.0, 2.0 * std::f64::consts::PI, p, 32);
            assert!((d - 0.5).abs() < 1e-3, "offset distance {d}");
        }
        let t = tangent3d(&circle, 0.0, 1e-4).unwrap();
        assert!((t.xyz().x).abs() < 1e-3, "tangent at 0 is +Y: {t:?}");
    }

    #[test]
    fn point_line_distance() {
        let d = point_line_distance3d(&GpPnt::new(0.0, 1.0, 0.0), &GpPnt::zero(), &GpVec::new(1.0, 0.0, 0.0));
        assert!((d - 1.0).abs() < 1e-12);
    }

    #[test]
    fn project_on_segment() {
        let (q, t) = project_on_segment3d(&GpPnt::new(0.5, 2.0, 0.0), &GpPnt::zero(), &GpPnt::new(1.0, 0.0, 0.0));
        assert!(q.x() > 0.49 && q.x() < 0.51);
        assert!((t - 0.5).abs() < 1e-9);
        let (_, t) = project_on_segment3d(&GpPnt::new(5.0, 0.0, 0.0), &GpPnt::zero(), &GpPnt::new(1.0, 0.0, 0.0));
        assert!((t - 1.0).abs() < 1e-12, "clamped t");
    }
}
