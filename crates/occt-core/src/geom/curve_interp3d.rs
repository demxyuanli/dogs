//! 3D curve interpolation — construct smooth parametric curves through a
//! point set.
//!
//! Port of `GeomAPI_PointsToBSpline` (global cubic interpolation) and
//! `GeomAPI_Interpolate`, plus Catmull–Rom and arc-length reparameterization
//! helpers. Source: `GeomAPI` (TKGeomBase), `Geom_BSplineCurve`.

use crate::gp::GpPnt;

/// A parametric 3D curve evaluated via a closure. Used to abstract the
/// interpolated result (Catmull-Rom, B-spline, or polyline) behind one type.
#[derive(Debug, Clone)]
pub struct InterpCurve3d {
    /// Evaluate the curve at parameter t (typically 0..=1 or knot range).
    eval: fn(f64, &InterpData) -> GpPnt,
    data: InterpData,
}

#[derive(Debug, Clone)]
enum InterpData {
    CatmullRom { points: Vec<GpPnt> },
    Polyline { points: Vec<GpPnt> },
    CubicSpline { knots: Vec<f64>, coeffs: Vec<[[f64; 3]; 4]> },
}

impl InterpCurve3d {
    fn new(data: InterpData, kind: InterpKind) -> Self {
        let eval: fn(f64, &InterpData) -> GpPnt = match kind {
            InterpKind::CatmullRom => eval_catmull,
            InterpKind::Polyline => eval_polyline,
            InterpKind::CubicSpline => eval_cubic,
        };
        Self { eval, data }
    }

    /// Evaluate at parameter `t` (0..=1 for Catmull-Rom/polyline, knot range
    /// for cubic spline).
    pub fn point(&self, t: f64) -> GpPnt {
        (self.eval)(t, &self.data)
    }

    /// Sample the curve at `n` evenly spaced parameters.
    pub fn sample(&self, n: usize) -> Vec<GpPnt> {
        let n = n.max(2);
        let (a, b) = self.range();
        (0..n)
            .map(|i| {
                let t = a + (b - a) * i as f64 / (n - 1) as f64;
                self.point(t)
            })
            .collect()
    }

    /// Parameter range.
    pub fn range(&self) -> (f64, f64) {
        match &self.data {
            InterpData::CubicSpline { knots, .. } => (knots[0], knots[knots.len() - 1]),
            _ => (0.0, 1.0),
        }
    }

    /// The knot parameters at which the data points are interpolated (for a
    /// cubic spline: the chord-length parameters of each input point).
    pub fn data_knots(&self) -> Vec<f64> {
        match &self.data {
            InterpData::CubicSpline { knots, .. } => knots.clone(),
            _ => (0..2).map(|i| i as f64).collect(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum InterpKind {
    CatmullRom,
    Polyline,
    CubicSpline,
}

/// Catmull–Rom evaluation with centripetal-style parameterization over [0,1].
/// Each segment uses the four surrounding control points (clamped at ends).
fn eval_catmull(t: f64, data: &InterpData) -> GpPnt {
    let InterpData::CatmullRom { points } = data else { return GpPnt::zero() };
    let n = points.len();
    if n == 0 {
        return GpPnt::zero();
    }
    if n == 1 {
        return points[0];
    }
    // Map t ∈ [0,1] to a segment index + local u ∈ [0,1].
    let seg_count = (n - 1) as f64;
    let tt = t.clamp(0.0, 1.0);
    let seg = ((tt * seg_count) as usize).min(n - 2);
    let u = tt * seg_count - seg as f64;

    let p0 = if seg == 0 { points[0] } else { points[seg - 1] };
    let p1 = points[seg];
    let p2 = points[seg + 1];
    let p3 = if seg + 2 < n { points[seg + 2] } else { points[n - 1] };

    // Standard Catmull–Rom basis on the segment p1→p2 (neighbors p0, p3):
    //   p(u) = 0.5·[ 2p1 + (−p0+p2)u + (2p0−5p1+4p2−p3)u² + (−p0+3p1−3p2+p3)u³ ]
    // At u=0 → p1, at u=1 → p2 (exactly through the control points).
    let u2 = u * u;
    let u3 = u2 * u;
    let x = 0.5
        * (2.0 * p1.x()
            + (-p0.x() + p2.x()) * u
            + (2.0 * p0.x() - 5.0 * p1.x() + 4.0 * p2.x() - p3.x()) * u2
            + (-p0.x() + 3.0 * p1.x() - 3.0 * p2.x() + p3.x()) * u3);
    let y = 0.5
        * (2.0 * p1.y()
            + (-p0.y() + p2.y()) * u
            + (2.0 * p0.y() - 5.0 * p1.y() + 4.0 * p2.y() - p3.y()) * u2
            + (-p0.y() + 3.0 * p1.y() - 3.0 * p2.y() + p3.y()) * u3);
    let z = 0.5
        * (2.0 * p1.z()
            + (-p0.z() + p2.z()) * u
            + (2.0 * p0.z() - 5.0 * p1.z() + 4.0 * p2.z() - p3.z()) * u2
            + (-p0.z() + 3.0 * p1.z() - 3.0 * p2.z() + p3.z()) * u3);
    GpPnt::new(x, y, z)
}

/// Piecewise-linear evaluation (degree-1 interpolation).
fn eval_polyline(t: f64, data: &InterpData) -> GpPnt {
    let InterpData::Polyline { points } = data else { return GpPnt::zero() };
    let n = points.len();
    if n == 0 {
        return GpPnt::zero();
    }
    if n == 1 {
        return points[0];
    }
    let tt = t.clamp(0.0, 1.0) * (n - 1) as f64;
    let i = (tt.floor() as usize).min(n - 2);
    let u = tt - i as f64;
    let a = points[i].coord;
    let b = points[i + 1].coord;
    GpPnt::from_xyz(&a.multiplied(1.0 - u).added(&b.multiplied(u)))
}

/// Cubic B-spline segment evaluation (uniform, per-coordinate coefficients).
fn eval_cubic(t: f64, data: &InterpData) -> GpPnt {
    let InterpData::CubicSpline { knots, coeffs } = data else { return GpPnt::zero() };
    let n = coeffs.len();
    if n == 0 {
        return GpPnt::zero();
    }
    // Find the segment containing t.
    let mut seg = 0;
    for (i, k) in knots.iter().enumerate() {
        if *k <= t + 1e-12 {
            seg = i;
        }
    }
    seg = seg.min(n - 1);
    let k0 = knots[seg];
    let k1 = knots[seg + 1];
    let u = if (k1 - k0).abs() > 1e-30 {
        ((t - k0) / (k1 - k0)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let c = &coeffs[seg];
    // p(u) = a + b·u + c·u² + d·u³ (per coordinate).
    let u2 = u * u;
    let u3 = u2 * u;
    GpPnt::new(
        c[0][0] + c[1][0] * u + c[2][0] * u2 + c[3][0] * u3,
        c[0][1] + c[1][1] * u + c[2][1] * u2 + c[3][1] * u3,
        c[0][2] + c[1][2] * u + c[2][2] * u2 + c[3][2] * u3,
    )
}

/// Build a Catmull–Rom curve through `points`.
pub fn catmull_rom(points: &[GpPnt]) -> InterpCurve3d {
    InterpCurve3d::new(InterpData::CatmullRom { points: points.to_vec() }, InterpKind::CatmullRom)
}

/// Build a degree-1 polyline through `points`.
pub fn polyline_curve(points: &[GpPnt]) -> InterpCurve3d {
    InterpCurve3d::new(InterpData::Polyline { points: points.to_vec() }, InterpKind::Polyline)
}

/// Global cubic B-spline interpolation through `points` (GeomAPI_PointsToBSpline).
/// Uses chord-length parameterization + natural (zero second derivative) end
/// conditions, solving the tridiagonal system per coordinate.
pub fn cubic_spline_interp(points: &[GpPnt]) -> Result<InterpCurve3d, String> {
    if points.len() < 2 {
        return Err("cubic_spline_interp: need at least 2 points".into());
    }
    let n = points.len();
    // Chord-length parameters (normalized to [0,1]).
    let mut knots = vec![0.0f64; n];
    let mut total = 0.0;
    for i in 1..n {
        total += points[i - 1].distance(&points[i]);
    }
    if total < 1e-12 {
        return Err("cubic_spline_interp: coincident points".into());
    }
    let mut acc = 0.0;
    for i in 1..n {
        acc += points[i - 1].distance(&points[i]);
        knots[i] = acc / total;
    }
    // Natural cubic spline second-derivative coefficients (tridiagonal solve).
    // h[i] = knots[i+1] − knots[i] for i in 0..n−1.
    let mut h = vec![0.0f64; n - 1];
    for i in 0..n - 1 {
        h[i] = knots[i + 1] - knots[i];
    }
    // Solve for the interior second derivatives M[1..n−2], with M[0]=M[n−1]=0
    // (natural end conditions). m = n−2 interior unknowns; equation i (0-based
    // over M[1..n−2]) involves M[i], M[i+1], M[i+2]:
    //   (h[i]/6)·M[i] + ((h[i]+h[i+1])/3)·M[i+1] + (h[i+1]/6)·M[i+2]
    //     = (p[i+2]−p[i+1])/h[i+1] − (p[i+1]−p[i])/h[i]
    let m = n - 2;
    let mut a = vec![0.0f64; m];
    let mut b = vec![0.0f64; m];
    let mut c = vec![0.0f64; m];
    let mut dx = vec![0.0f64; m];
    let mut dy = vec![0.0f64; m];
    let mut dz = vec![0.0f64; m];
    for i in 0..m {
        a[i] = h[i] / 6.0;
        b[i] = (h[i] + h[i + 1]) / 3.0;
        c[i] = h[i + 1] / 6.0;
        dx[i] = (points[i + 2].x() - points[i + 1].x()) / h[i + 1]
            - (points[i + 1].x() - points[i].x()) / h[i];
        dy[i] = (points[i + 2].y() - points[i + 1].y()) / h[i + 1]
            - (points[i + 1].y() - points[i].y()) / h[i];
        dz[i] = (points[i + 2].z() - points[i + 1].z()) / h[i + 1]
            - (points[i + 1].z() - points[i].z()) / h[i];
    }
    let mx = solve_tridiagonal(&a, &b, &c, &dx)?;
    let my = solve_tridiagonal(&a, &b, &c, &dy)?;
    let mz = solve_tridiagonal(&a, &b, &c, &dz)?;

    // Per-segment cubic coefficients [c0, c1, c2, c3] per coordinate, in
    // LOCAL segment parameter u ∈ [0,1] (segment i maps to u ∈ [0,1] regardless
    // of the actual knot spacing, since the natural spline segment is written
    // in terms of the chord-length parameter — we use the standard Hermite
    // form with second derivatives M).
    let _m_len = mx.len();
    let mut coeffs: Vec<[[f64; 3]; 4]> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        // Second derivative at node i (M_i) and node i+1 (M_{i+1}); natural
        // ends (node 0 and node n−1) are 0. Interior node j ∈ [1, n−2] has
        // M_j = arr[j−1].
        let m_at = |j: usize, arr: &[f64]| -> f64 {
            if j == 0 || j == n - 1 {
                0.0
            } else if j - 1 < arr.len() {
                arr[j - 1]
            } else {
                0.0
            }
        };
        let mi = |arr: &[f64]| m_at(i, arr);
        let mi1 = |arr: &[f64]| m_at(i + 1, arr);
        let hi = h[i];
        let h2 = hi * hi;
        let p0 = &points[i];
        let p1 = &points[i + 1];
        let (mx_i, mx_i1) = (mi(&mx), mi1(&mx));
        let (my_i, my_i1) = (mi(&my), mi1(&my));
        let (mz_i, mz_i1) = (mi(&mz), mi1(&mz));
        // Natural cubic spline segment in local u = t/hi ∈ [0,1] (Hi = hi²):
        //   p(u) = (1−u)p0 + u·p1 + Hi/6·(1−u)³·M_i + Hi/6·u³·M_{i+1}
        //         − Hi/6·(1−u)·M_i − Hi/6·u·M_{i+1}
        // (The (1−u)³ and (1−u) terms are the Hermite basis; the −Hi/6 terms
        // ensure p(0)=p0, p(1)=p1, p''(0)=M_i, p''(1)=M_{i+1}.)
        // Expanding:
        //   c0 = p0
        //   c1 = (p1−p0) + Hi/6·(−2·M_i − M_{i+1})
        //   c2 = Hi/2·M_i
        //   c3 = Hi/6·(M_{i+1} − M_i)
        let s6 = |m: f64| h2 / 6.0 * m;
        let c0 = [p0.x(), p0.y(), p0.z()];
        let c1 = [
            (p1.x() - p0.x()) + s6(-2.0 * mx_i - mx_i1),
            (p1.y() - p0.y()) + s6(-2.0 * my_i - my_i1),
            (p1.z() - p0.z()) + s6(-2.0 * mz_i - mz_i1),
        ];
        let c2 = [
            h2 / 2.0 * mx_i,
            h2 / 2.0 * my_i,
            h2 / 2.0 * mz_i,
        ];
        let c3 = [
            s6(mx_i1 - mx_i),
            s6(my_i1 - my_i),
            s6(mz_i1 - mz_i),
        ];
        coeffs.push([c0, c1, c2, c3]);
    }
    Ok(InterpCurve3d::new(InterpData::CubicSpline { knots, coeffs }, InterpKind::CubicSpline))
}

/// Thomas algorithm for a tridiagonal system.
fn solve_tridiagonal(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Result<Vec<f64>, String> {
    let n = d.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut cp = vec![0.0f64; n];
    let mut dp = vec![0.0f64; n];
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for i in 1..n {
        let denom = b[i] - a[i] * cp[i - 1];
        if denom.abs() < 1e-30 {
            return Err("solve_tridiagonal: singular".into());
        }
        cp[i] = if i == n - 1 { 0.0 } else { c[i] / denom };
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }
    let mut x = vec![0.0f64; n];
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
    Ok(x)
}

/// Arc length of an interpolated curve via sampling.
pub fn interp_arc_length(c: &InterpCurve3d, samples: usize) -> f64 {
    let (a, b) = c.range();
    let n = samples.max(2);
    let mut len = 0.0;
    let mut prev = c.point(a);
    for i in 1..=n {
        let t = a + (b - a) * i as f64 / n as f64;
        let cur = c.point(t);
        len += prev.distance(&cur);
        prev = cur;
    }
    len
}

/// Reparameterize by arc length: returns (s, t) pairs at n uniform arc-length
/// positions (GCPnts_UniformAbscissa over the interpolated curve).
pub fn arc_length_reparameterize(c: &InterpCurve3d, n: usize, samples: usize) -> Vec<(f64, f64)> {
    if n < 2 {
        return vec![(0.0, c.range().0)];
    }
    let (a, b) = c.range();
    let samples = samples.max(8);
    let mut cum = vec![0.0f64; samples + 1];
    let mut ts = vec![0.0f64; samples + 1];
    let mut prev = c.point(a);
    ts[0] = a;
    for i in 1..=samples {
        let t = a + (b - a) * i as f64 / samples as f64;
        let p = c.point(t);
        cum[i] = cum[i - 1] + prev.distance(&p);
        prev = p;
        ts[i] = t;
    }
    let total = cum[samples];
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let target = total * k as f64 / (n - 1) as f64;
        let mut idx = 0;
        while idx < samples && cum[idx + 1] < target {
            idx += 1;
        }
        let seg = cum[idx + 1] - cum[idx];
        let u = if seg > 1e-30 { (target - cum[idx]) / seg } else { 0.0 };
        let t = ts[idx] + u * (ts[idx + 1] - ts[idx]);
        out.push((target, t));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catmull_through_points() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        let c = catmull_rom(&pts);
        // Catmull-Rom passes through interior points exactly.
        let mid = c.point(0.5);
        assert!((mid.x() - 1.0).abs() < 0.05 && (mid.y() - 1.0).abs() < 0.05, "mid {mid:?}");
    }

    #[test]
    fn polyline_corners() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        let c = polyline_curve(&pts);
        assert!((c.point(0.0).x() - 0.0).abs() < 1e-9);
        assert!((c.point(0.5).x() - 1.0).abs() < 1e-9);
        assert!((c.point(1.0).x() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn cubic_spline_passes_points() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 2.0, 0.0),
            GpPnt::new(2.0, -1.0, 0.0),
            GpPnt::new(3.0, 3.0, 0.0),
        ];
        let c = cubic_spline_interp(&pts).expect("spline");
        // The curve passes through each point at its chord-length knot.
        let knots = c.data_knots();
        assert_eq!(knots.len(), 4);
        for (i, p) in pts.iter().enumerate() {
            let q = c.point(knots[i]);
            assert!(p.distance(&q) < 1e-6, "point {i}: {p:?} vs {q:?}");
        }
    }

    #[test]
    fn cubic_spline_errors() {
        assert!(cubic_spline_interp(&[GpPnt::new(0.0, 0.0, 0.0)]).is_err());
        // Duplicate points → error.
        let dup = vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(0.0, 0.0, 0.0)];
        assert!(cubic_spline_interp(&dup).is_err());
    }

    #[test]
    fn spline_line_is_line() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
            GpPnt::new(3.0, 0.0, 0.0),
        ];
        let c = cubic_spline_interp(&pts).expect("spline");
        let (a, b) = c.range();
        let q = c.point((a + b) / 2.0);
        assert!((q.y() - 0.0).abs() < 1e-9 && (q.z() - 0.0).abs() < 1e-9, "mid {q:?}");
    }

    #[test]
    fn arc_length_polyline() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(3.0, 4.0, 0.0),
        ];
        let c = polyline_curve(&pts);
        let len = interp_arc_length(&c, 32);
        assert!((len - 5.0).abs() < 1e-6, "len {len}");
    }

    #[test]
    fn arc_reparameterize_uniform() {
        let pts = vec![
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        let c = polyline_curve(&pts);
        let p = arc_length_reparameterize(&c, 5, 64);
        assert_eq!(p.len(), 5);
        let step = p[1].0 - p[0].0;
        for k in 1..5 {
            assert!((p[k].0 - p[k - 1].0 - step).abs() < 1e-6, "uniform spacing");
        }
    }

    #[test]
    fn sample_count() {
        let pts = vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(1.0, 0.0, 0.0)];
        let c = catmull_rom(&pts);
        let s = c.sample(10);
        assert_eq!(s.len(), 10);
    }
}
