//! Curve reparameterization and reversal helpers.
//! Source: `Geom_Curve`, `Geom_BSplineCurve`, `GCPnts_AbscissaPoint`.

use std::sync::Arc;
use crate::bspline_curve::GeomBSplineCurve;
use crate::curve::Curve;
use occt_core::bspl::knots::{build_uniform_knots, hunt};
use occt_core::bspl::poles::greville_abscissae;
use occt_core::gp::{GpPnt, GpTrsf, GpVec};

/// Linear reparameterization adapter: maps a new parameter `t ∈ [new_a, new_b]`
/// onto the original curve's parameter range `[first, last]` linearly.
///
/// When the basis is a BSpline / Bezier, knots (remapped) and degree are kept
/// so `Extrema_ExtPC` still takes the `GeomAbs_BSplineCurve` arm
/// (`Extrema_GGExtPC.hxx:190+`) after SameRange / Project reparam.
#[derive(Clone)]
pub struct ReparamCurve {
    curve: Arc<dyn Curve>,
    first: f64,
    last: f64,
    new_a: f64,
    new_b: f64,
    /// Flat knot sequence mapped into `[new_a, new_b]`; `None` if basis is not
    /// a BSpline.
    remapped_knots: Option<Vec<f64>>,
    nurbs_deg: Option<usize>,
}

impl ReparamCurve {
    fn orig(&self, t: f64) -> f64 {
        self.first + (t - self.new_a) * (self.last - self.first) / (self.new_b - self.new_a)
    }

    fn map_param(&self, u: f64) -> f64 {
        let den = self.last - self.first;
        if den.abs() <= 0.0 {
            return self.new_a;
        }
        self.new_a + (u - self.first) * (self.new_b - self.new_a) / den
    }
}

impl Curve for ReparamCurve {
    fn d0(&self, t: f64) -> GpPnt { self.curve.d0(self.orig(t)) }
    fn d1(&self, t: f64) -> (GpPnt, GpVec) {
        let s = (self.last - self.first) / (self.new_b - self.new_a);
        let (p, d) = self.curve.d1(self.orig(t));
        (p, d.multiplied_scalar(s))
    }
    fn d2(&self, t: f64) -> (GpPnt, GpVec, GpVec) {
        let s = (self.last - self.first) / (self.new_b - self.new_a);
        let (p, d1, d2) = self.curve.d2(self.orig(t));
        (p, d1.multiplied_scalar(s), d2.multiplied_scalar(s * s))
    }
    fn first_parameter(&self) -> f64 { self.new_a }
    fn last_parameter(&self) -> f64 { self.new_b }
    fn is_periodic(&self) -> bool { self.curve.is_periodic() }
    fn period(&self) -> f64 {
        let den = self.last - self.first;
        if den.abs() <= 0.0 {
            return self.curve.period();
        }
        self.curve.period() * (self.new_b - self.new_a).abs() / den.abs()
    }
    fn continuity(&self) -> u8 { self.curve.continuity() }
    fn is_line(&self) -> bool { self.curve.is_line() }
    fn gp_circ(&self) -> Option<occt_core::gp::GpCirc> { self.curve.gp_circ() }
    fn circle_radius(&self) -> Option<f64> { self.curve.circle_radius() }
    fn nurbs_degree(&self) -> Option<usize> { self.nurbs_deg.or_else(|| self.curve.nurbs_degree()) }
    fn bspline_knots(&self) -> Option<&[f64]> {
        self.remapped_knots.as_deref()
    }
    fn bspline_poles(&self) -> Option<&[GpPnt]> { self.curve.bspline_poles() }
    fn bspline_weights(&self) -> Option<&[f64]> { self.curve.bspline_weights() }
    fn bezier_poles(&self) -> Option<&[GpPnt]> { self.curve.bezier_poles() }
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        self.curve
            .parameter_intervals(continuity)
            .into_iter()
            .map(|u| self.map_param(u))
            .collect()
    }
    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity)
            .len()
            .saturating_sub(1)
            .max(1) as i32
    }
    fn resolution(&self, r3d: f64) -> f64 {
        let den = (self.last - self.first).abs();
        let span = (self.new_b - self.new_a).abs();
        if den <= 0.0 {
            return self.curve.resolution(r3d);
        }
        self.curve.resolution(r3d) * span / den
    }
    // ponytail: adapter over an immutable Arc; transforms are no-ops.
    fn transform(&mut self, _t: &GpTrsf) {}
    fn reverse(&mut self) {
        // Swap the new range only. Remapped knots stay in the geometric
        // image of the basis knot vector; ExtPC normalizes [uinf,usup].
        std::mem::swap(&mut self.new_a, &mut self.new_b);
    }
    fn clone_dyn(&self) -> Box<dyn Curve> { Box::new(self.clone()) }
}

/// Reparameterize `c` so the returned curve's range is `[new_a, new_b]` while
/// preserving geometry: new `t` maps linearly onto `[first, last]`. Requires a
/// bounded original range (unbounded lines yield a degenerate map).
pub fn reparameterize_curve(c: &dyn Curve, new_a: f64, new_b: f64) -> Arc<dyn Curve> {
    let first = c.first_parameter();
    let last = c.last_parameter();
    let den = last - first;
    let remapped_knots = c.bspline_knots().map(|knots| {
        if den.abs() <= 0.0 {
            knots.to_vec()
        } else {
            knots
                .iter()
                .map(|&u| new_a + (u - first) * (new_b - new_a) / den)
                .collect()
        }
    });
    Arc::new(ReparamCurve {
        curve: Arc::from(c.clone_dyn()),
        first,
        last,
        new_a,
        new_b,
        remapped_knots,
        nurbs_deg: c.nurbs_degree(),
    })
}

/// Reparameterize then trim conceptually: returns the reparameterized curve
/// restricted to `[u0, u1]` (i.e. the reparam adapter with that new range).
pub fn compose_reparam(c: &dyn Curve, u0: f64, u1: f64) -> Arc<dyn Curve> {
    reparameterize_curve(c, u0, u1)
}

/// Reversed copy of `c` via `Curve::reversed`.
pub fn curve_reverse_copy(c: &dyn Curve) -> Arc<dyn Curve> {
    Arc::from(c.reversed())
}

/// Maximum deviation between `c` and its linear reparameterization (the chord
/// from `c(first)` to `c(last)`, sampled at `samples` parameters). Zero for
/// straight curves; diagnostic for how far the curve is from linear.
pub fn arc_reparam_check(c: &dyn Curve, samples: usize) -> f64 {
    let a = c.first_parameter();
    let b = c.last_parameter();
    if !a.is_finite() || !b.is_finite() || b <= a {
        return f64::INFINITY;
    }
    let pa = c.d0(a);
    let pb = c.d0(b);
    let span = b - a;
    let n = samples.max(1);
    let mut max_dev: f64 = 0.0;
    for i in 0..=n {
        let t = a + span * i as f64 / n as f64;
        let p = c.d0(t);
        let f = (t - a) / span;
        let q = GpPnt::new(
            pa.x() + f * (pb.x() - pa.x()),
            pa.y() + f * (pb.y() - pa.y()),
            pa.z() + f * (pb.z() - pa.z()),
        );
        max_dev = max_dev.max(p.distance(&q));
    }
    max_dev
}

/// B-spline basis values at `u`, indexed by pole (Algorithm A2.2 from
/// The NURBS Book), for a clamped knot vector.
///
/// **非 OCCT 出处**：OCCT 的等价件是 `BSplCLib::BasisFuns`/`bspl`（按结点区间
/// 与 `Span` 求基函数），本函数取自 The NURBS Book A2.2 教科书算法，故文件头
/// 只把 `Geom_Curve`/`Geom_BSplineCurve`/`GCPnts_AbscissaPoint` 记为*调用侧*来源。
fn basis_values(knots: &[f64], degree: usize, u: f64) -> Vec<f64> {
    let n_poles = knots.len() - degree - 1;
    let s = hunt(knots, u).max(degree).min(n_poles - 1);
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    let mut n = vec![0.0; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = u - knots[s + 1 - j];
        right[j] = knots[s + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() > 1e-15 { n[r] / denom } else { 0.0 };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    let mut basis = vec![0.0; n_poles];
    for r in 0..=degree {
        let gi = s - degree + r;
        if gi < n_poles {
            basis[gi] = n[r];
        }
    }
    basis
}

/// Solve `A·x = b` for a square `A` by partial-pivoted Gaussian elimination;
/// `None` if singular.
///
/// The singularity test is OCCT's `math_Gauss` pivot threshold
/// `MinPivot = 1.0e-20` (`math_Gauss.hxx:45-49`: "If the largest pivot found is
/// less than MinPivot the matrix A is considered singular"). The previous
/// constant was the invented `1e-14` (audit A15).
fn gauss_solve(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let mut m: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = a[i].clone();
            row.push(b[i]);
            row
        })
        .collect();
    for col in 0..n {
        let mut piv = col;
        let mut best = m[col][col].abs();
        for r in (col + 1)..n {
            if m[r][col].abs() > best {
                best = m[r][col].abs();
                piv = r;
            }
        }
        if best < 1.0e-20 {
            return None;
        }
        if piv != col {
            m.swap(piv, col);
        }
        let pv = m[col][col];
        for r in (col + 1)..n {
            let f = m[r][col] / pv;
            if f == 0.0 {
                continue;
            }
            for c in col..n + 1 {
                m[r][c] -= f * m[col][c];
            }
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = m[r][n];
        for c in (r + 1)..n {
            s -= m[r][c] * x[c];
        }
        x[r] = s / m[r][r];
    }
    Some(x)
}

/// Approximate/interpolate `points` with a degree-`degree` clamped B-spline.
///
/// Degree 1: control points are the data points. Degree ≥ 2: uniform clamped
/// knots with a collocation solve at the Greville abscissae (global
/// interpolation, as in `BSplCLib::Interpolate`).
pub fn resample_bspline(points: &[GpPnt], degree: usize) -> Result<GeomBSplineCurve, String> {
    let n = points.len();
    if n == 0 {
        return Err("resample_bspline: no points".into());
    }
    if degree == 0 {
        return Err("resample_bspline: degree must be >= 1".into());
    }
    if n < degree + 1 {
        return Err(format!("resample_bspline: need at least {} points for degree {degree}", degree + 1));
    }
    let knots = build_uniform_knots(n, degree);
    if degree == 1 {
        return GeomBSplineCurve::new(points.to_vec(), knots, 1).map_err(|e| e.to_string());
    }
    let params = greville_abscissae(&knots, degree, n);
    let a: Vec<Vec<f64>> = params.iter().map(|&u| basis_values(&knots, degree, u)).collect();
    let px = gauss_solve(&a, &points.iter().map(|p| p.x()).collect::<Vec<_>>())
        .ok_or("resample_bspline: singular collocation matrix")?;
    let py = gauss_solve(&a, &points.iter().map(|p| p.y()).collect::<Vec<_>>())
        .ok_or("resample_bspline: singular collocation matrix")?;
    let pz = gauss_solve(&a, &points.iter().map(|p| p.z()).collect::<Vec<_>>())
        .ok_or("resample_bspline: singular collocation matrix")?;
    let poles: Vec<GpPnt> = (0..n).map(|i| GpPnt::new(px[i], py[i], pz[i])).collect();
    GeomBSplineCurve::new(poles, knots, degree).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A straight B-spline line from (0,0,0) to (10,0,0), range [0, 10].
    fn line_0_10() -> GeomBSplineCurve {
        GeomBSplineCurve::new(
            vec![GpPnt::new(0.0, 0.0, 0.0), GpPnt::new(10.0, 0.0, 0.0)],
            vec![0.0, 0.0, 10.0, 10.0],
            1,
        )
        .unwrap()
    }

    #[test]
    fn reparameterize_line_0_10_to_0_1() {
        let c = line_0_10();
        let r = reparameterize_curve(&c, 0.0, 1.0);
        assert!((r.first_parameter() - 0.0).abs() < 1e-12);
        assert!((r.last_parameter() - 1.0).abs() < 1e-12);
        // New t=0.5 maps to original u=5.
        assert!(r.d0(0.5).distance(&c.d0(5.0)) < 1e-9);
        // Derivative scales by orig'(t) = (10-0)/(1-0) = 10.
        let (_, d) = r.d1(0.5);
        assert!((d.x() - 10.0).abs() < 1e-9 && d.y().abs() < 1e-9 && d.z().abs() < 1e-9, "d={d:?}");
    }

    #[test]
    fn resample_bspline_degree1_through_points() {
        let pts = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(2.0, 0.0, 0.0),
        ];
        let c = resample_bspline(&pts, 1).expect("resample");
        assert!(c.d0(0.0).distance(&pts[0]) < 1e-9);
        assert!(c.d0(0.5).distance(&pts[1]) < 1e-9);
        assert!(c.d0(1.0).distance(&pts[2]) < 1e-9);
    }

    #[test]
    fn resample_bspline_cubic_interpolates() {
        let pts = [
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(-1.0, 0.0, 0.0),
            GpPnt::new(0.0, -1.0, 0.0),
        ];
        let c = resample_bspline(&pts, 3).expect("resample");
        let params = greville_abscissae(&c.knots, 3, 4);
        for (i, &u) in params.iter().enumerate() {
            assert!(c.d0(u).distance(&pts[i]) < 1e-6, "point {i} at u={u}: got {:?}, expected {:?}", c.d0(u), pts[i]);
        }
    }

    #[test]
    fn reverse_copy_flips_endpoints() {
        let c = line_0_10();
        let orig0 = c.d0(0.0);
        let orig10 = c.d0(10.0);
        let r = curve_reverse_copy(&c);
        assert!(r.d0(0.0).distance(&orig10) < 1e-9);
        assert!(r.d0(10.0).distance(&orig0) < 1e-9);
    }

    #[test]
    fn compose_reparam_restricts_range() {
        let c = line_0_10();
        let r = compose_reparam(&c, 2.0, 4.0);
        assert!((r.first_parameter() - 2.0).abs() < 1e-12);
        assert!((r.last_parameter() - 4.0).abs() < 1e-12);
        assert!(r.d0(2.0).distance(&c.d0(0.0)) < 1e-9);
        assert!(r.d0(4.0).distance(&c.d0(10.0)) < 1e-9);
    }

    #[test]
    fn arc_reparam_check_line_is_zero() {
        let c = line_0_10();
        let dev = arc_reparam_check(&c, 32);
        assert!(dev < 1e-9, "dev={dev}");
    }
}
