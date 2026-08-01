//! Least-squares B-spline curve fit to a 3-D point cloud.
//! Source: `AppParCurves` / `math_GaussLeastSquare`-adjacent.
//!
//! Points are parameterized by chord length, a clamped uniform knot vector is
//! built for the requested pole count, and the control poles are obtained from
//! the regularized normal equations `(BᵀB + λI) p = Bᵀ y` solved per
//! coordinate. The same Cox–de Boor machinery as `spline_surface` is used for
//! evaluation.

use crate::{MathMatrix, MathVector};

/// Result of a least-squares B-spline curve fit.
///
/// `knots` is the clamped uniform knot vector for `poles.len()` poles of degree
/// `deg`; `weights` is `None` for the polynomial fits produced here (kept for
/// API symmetry with NURBS). `error` is the maximum distance between the fitted
/// curve and the input points at their chord-length parameters.
#[derive(Debug, Clone)]
pub struct BSplineFit {
    pub knots: Vec<f64>,
    pub deg: usize,
    pub poles: Vec<[f64; 3]>,
    pub weights: Option<Vec<f64>>,
    pub error: f64,
}

/// Fit a degree-`deg` B-spline curve with `num_poles` control poles through
/// `points` in the least-squares sense.
///
/// The algorithm parameterizes the points by cumulative chord length onto
/// `[0, 1]`, builds the clamped uniform knot vector for `num_poles` poles, then
/// solves the ridge-regularized normal equations `(BᵀB + λI) p = Bᵀ y` for each
/// coordinate (x, y, z). The `error` field of the result holds the maximum
/// `‖fitted(t_i) − P_i‖`.
///
/// # Errors
/// Returns an error when there are fewer than two points, the degree is zero,
/// `num_poles` does not exceed the degree, or `num_poles` exceeds the point
/// count (an underdetermined system).
pub fn fit_bspline_curve(points: &[[f64; 3]], deg: usize, num_poles: usize) -> Result<BSplineFit, String> {
    let params = chord_length_parameters(points);
    fit_with_parameters(points, deg, num_poles, &params)
}

/// Fit a B-spline curve using centripetal parameterization
/// (`t_{i+1} − t_i = √‖P_{i+1} − P_i‖`).
///
/// Centripetal parameterization distributes the knot spacing more evenly for
/// strongly curved data and often reduces the fit error on loops and arcs
/// compared to plain chord-length parameterization.
pub fn fit_bspline_curve_centripetal(
    points: &[[f64; 3]],
    deg: usize,
    num_poles: usize,
) -> Result<BSplineFit, String> {
    let params = centripetal_parameters(points);
    fit_with_parameters(points, deg, num_poles, &params)
}

/// Knot count sanity helper: a degree-`deg` B-spline curve with `num_poles`
/// poles has `num_poles + deg + 1` knots (`#poles = #knots − deg − 1`).
pub fn bspline_curve_knot_count(num_poles: usize, deg: usize) -> usize {
    num_poles + deg + 1
}

/// Build a [`BSplineFit`] from explicit knot vector, degree, poles and
/// optional weights, validating that the knot count matches the pole count and
/// that the knots are non-decreasing.
///
/// The `error` field is initialized to 0 (the caller should compute residuals
/// against the fitted data with [`fit_curve_residuals`]).
pub fn build_bspline_curve(
    knots: Vec<f64>,
    deg: usize,
    poles: Vec<[f64; 3]>,
    weights: Option<Vec<f64>>,
) -> Result<BSplineFit, String> {
    let num_poles = poles.len();
    if num_poles == 0 {
        return Err("build_bspline_curve: no poles".to_string());
    }
    if knots.len() != bspline_curve_knot_count(num_poles, deg) {
        return Err(format!(
            "build_bspline_curve: knot count {} != poles {} + deg {} + 1",
            knots.len(),
            num_poles,
            deg
        ));
    }
    if knots.windows(2).any(|w| w[1] < w[0]) {
        return Err("build_bspline_curve: knot vector not non-decreasing".to_string());
    }
    if let Some(ws) = &weights {
        if ws.len() != num_poles {
            return Err("build_bspline_curve: weight count does not match poles".to_string());
        }
    }
    Ok(BSplineFit {
        knots,
        deg,
        poles,
        weights,
        error: 0.0,
    })
}

/// Central finite-difference derivative `dS/dt` of the fitted curve at `t`,
/// with step `h` (default `1e-6` when `h <= 0`).
pub fn eval_fit_curve_derivative(f: &BSplineFit, t: f64, h: f64) -> [f64; 3] {
    let h = if h > 0.0 { h } else { 1e-6 };
    let a = eval_curve_with(f, t + h);
    let b = eval_curve_with(f, t - h);
    [
        (a[0] - b[0]) / (2.0 * h),
        (a[1] - b[1]) / (2.0 * h),
        (a[2] - b[2]) / (2.0 * h),
    ]
}

/// Shared implementation: build the collocation matrix for `params`, solve the
/// regularized normal equations per coordinate, and compute the max residual.
fn fit_with_parameters(
    points: &[[f64; 3]],
    deg: usize,
    num_poles: usize,
    params: &[f64],
) -> Result<BSplineFit, String> {
    let n = points.len();
    if n < 2 {
        return Err("fit_bspline_curve: need at least 2 points".to_string());
    }
    if deg < 1 {
        return Err("fit_bspline_curve: degree must be >= 1".to_string());
    }
    if num_poles <= deg {
        return Err("fit_bspline_curve: num_poles must be > degree".to_string());
    }
    if num_poles > n {
        return Err("fit_bspline_curve: num_poles cannot exceed point count".to_string());
    }

    let knots = uniform_clamped_knots(num_poles, deg);
    let b = collocation(&knots, deg, params, num_poles);
    let poles = solve_normal(&b, points, num_poles)?;

    let mut f = BSplineFit {
        knots,
        deg,
        poles,
        weights: None,
        error: 0.0,
    };
    let mut error = 0.0;
    for (i, p) in points.iter().enumerate() {
        let fitted = eval_curve_with(&f, params[i]);
        let d = dist(fitted, *p);
        if d > error {
            error = d;
        }
    }
    f.error = error;
    Ok(f)
}

/// Evaluate the fitted B-spline curve at parameter `t` in `[0, 1]`.
///
/// The curve is clamped, so `eval_fit_curve(f, 0.0)` returns the first pole and
/// `eval_fit_curve(f, 1.0)` returns the last pole. When `weights` is `Some`,
/// the evaluation is rational (NURBS).
pub fn eval_fit_curve(f: &BSplineFit, t: f64) -> [f64; 3] {
    eval_curve_with(f, t)
}

/// Per-point residuals `‖fitted(t_i) − P_i‖` of the fit against `points`,
/// using the same chord-length parameterization the fit used.
pub fn fit_curve_residuals(f: &BSplineFit, points: &[[f64; 3]]) -> Vec<f64> {
    let params = chord_length_parameters(points);
    points
        .iter()
        .zip(&params)
        .map(|(p, &t)| dist(eval_curve_with(f, t), *p))
        .collect()
}

/// Maximum residual `‖fitted(t_i) − P_i‖` over `points`.
pub fn fit_curve_max_error(f: &BSplineFit, points: &[[f64; 3]]) -> f64 {
    fit_curve_residuals(f, points)
        .into_iter()
        .fold(0.0, f64::max)
}

/// Chord-length parameterization of `points` onto `[0, 1]`.
///
/// `t_0 = 0` and `t_i = t_{i−1} + ‖P_i − P_{i−1}‖`, normalized by the total
/// length. When all points coincide the result is all zeros.
pub fn chord_length_parameters(points: &[[f64; 3]]) -> Vec<f64> {
    let n = points.len();
    let mut t = vec![0.0; n];
    let mut acc = 0.0;
    for i in 1..n {
        acc += dist(points[i], points[i - 1]);
        t[i] = acc;
    }
    if acc > 0.0 {
        for x in t.iter_mut() {
            *x /= acc;
        }
    }
    t
}

/// Centripetal parameterization: `t_{i+1} − t_i = √‖P_{i+1} − P_i‖`,
/// normalized onto `[0, 1]`. This reduces parameter bunching on strongly
/// curved data (Lee, "Choosing nodes in parametric curve interpolation").
pub fn centripetal_parameters(points: &[[f64; 3]]) -> Vec<f64> {
    let n = points.len();
    let mut t = vec![0.0; n];
    let mut acc = 0.0;
    for i in 1..n {
        acc += dist(points[i], points[i - 1]).sqrt();
        t[i] = acc;
    }
    if acc > 0.0 {
        for x in t.iter_mut() {
            *x /= acc;
        }
    }
    t
}

/// Full basis vector `(N_{0,deg}(t), …, N_{np−1,deg}(t))` for the clamped knot
/// vector `knots`. Useful for partition-of-unity checks and collocation.
pub fn bspline_basis_1d(knots: &[f64], deg: usize, t: f64) -> Vec<f64> {
    let np = knots.len() - deg - 1;
    let mut out = vec![0.0; np];
    let span = find_span(knots, deg, t);
    let vals = basis_vals(span, deg, t, knots);
    for r in 0..=deg {
        let idx = span - deg + r;
        if idx < np {
            out[idx] = vals[r];
        }
    }
    out
}

/// Evaluate a (possibly rational) B-spline curve with explicit poles.
pub fn eval_bspline_curve(
    poles: &[[f64; 3]],
    weights: Option<&[f64]>,
    knots: &[f64],
    deg: usize,
    t: f64,
) -> [f64; 3] {
    let np = poles.len();
    if np == 0 {
        return [0.0; 3];
    }
    let span = find_span(knots, deg, t);
    let vals = basis_vals(span, deg, t, knots);
    let mut num = [0.0; 3];
    let mut denom = 0.0;
    for r in 0..=deg {
        let idx = span - deg + r;
        if idx >= np {
            continue;
        }
        let w = weights.map_or(1.0, |ws| ws[idx]);
        let b = vals[r] * w;
        num[0] += b * poles[idx][0];
        num[1] += b * poles[idx][1];
        num[2] += b * poles[idx][2];
        denom += b;
    }
    if denom.abs() < 1e-300 {
        return [0.0; 3];
    }
    [num[0] / denom, num[1] / denom, num[2] / denom]
}

/// Evaluation of a [`BSplineFit`] curve at `t`.
fn eval_curve_with(f: &BSplineFit, t: f64) -> [f64; 3] {
    eval_bspline_curve(&f.poles, f.weights.as_deref(), &f.knots, f.deg, t)
}

/// Clamped uniform knot vector for `n` poles of degree `deg`.
fn uniform_clamped_knots(n: usize, deg: usize) -> Vec<f64> {
    let interior = n.saturating_sub(deg + 1);
    let mut k = Vec::with_capacity(n + deg + 1);
    for _ in 0..=deg {
        k.push(0.0);
    }
    for i in 1..=interior {
        k.push(i as f64 / (interior + 1) as f64);
    }
    for _ in 0..=deg {
        k.push(1.0);
    }
    k
}

/// Find the knot span `[knots[s], knots[s+1])` containing `t`.
fn find_span(knots: &[f64], deg: usize, t: f64) -> usize {
    let n = knots.len() - deg - 2;
    if t <= knots[deg] {
        return deg;
    }
    if t >= knots[n + 1] {
        return n;
    }
    let mut lo = deg;
    let mut hi = n + 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if t < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo
}

/// Non-zero basis values on the span containing `t`; `out[r] = N_{span−deg+r}`.
fn basis_vals(span: usize, deg: usize, t: f64, knots: &[f64]) -> Vec<f64> {
    let mut left = vec![0.0; deg + 1];
    let mut right = vec![0.0; deg + 1];
    let mut n = vec![0.0; deg + 1];
    n[0] = 1.0;
    for j in 1..=deg {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let den = right[r + 1] + left[j - r];
            let temp = if den.abs() < 1e-300 { 0.0 } else { n[r] / den };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Build the collocation matrix `B[i][j] = N_j(t_i)` for parameters `params`.
fn collocation(knots: &[f64], deg: usize, params: &[f64], num_poles: usize) -> MathMatrix {
    let n = params.len();
    let mut b = MathMatrix::new(1, n, 1, num_poles);
    for (i, &t) in params.iter().enumerate() {
        let span = find_span(knots, deg, t);
        let vals = basis_vals(span, deg, t, knots);
        for r in 0..=deg {
            let j = span - deg + r;
            b.set_value(i + 1, j + 1, vals[r]);
        }
    }
    b
}

/// Solve the ridge-regularized normal equations per coordinate.
fn solve_normal(b: &MathMatrix, points: &[[f64; 3]], num_poles: usize) -> Result<Vec<[f64; 3]>, String> {
    let n = b.row_count();
    let lambda = 1e-10;
    let mut ata = MathMatrix::new(1, num_poles, 1, num_poles);
    for j in 1..=num_poles {
        for k in 1..=num_poles {
            let mut s = 0.0;
            for i in 1..=n {
                s += b.value(i, j) * b.value(i, k);
            }
            ata.set_value(j, k, s + if j == k { lambda } else { 0.0 });
        }
    }
    let mut poles = vec![[0.0; 3]; num_poles];
    for coord in 0..3 {
        let rhs = MathVector::from_slice(&points.iter().map(|p| p[coord]).collect::<Vec<_>>());
        let mut atb = MathVector::new(1, num_poles);
        for j in 1..=num_poles {
            let mut s = 0.0;
            for i in 1..=n {
                s += b.value(i, j) * rhs.value(i);
            }
            atb.set_value(j, s);
        }
        let sol = ata.solve(&atb).map_err(|e| e.to_string())?;
        for j in 0..num_poles {
            poles[j][coord] = sol.value(j + 1);
        }
    }
    Ok(poles)
}

/// Euclidean distance between two 3-D points.
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// `n` points on the line `(x, 2x+1, x)` for `x` in `[0, 1]`.
    fn line_points(n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let x = i as f64 / (n - 1) as f64;
                [x, 2.0 * x + 1.0, x]
            })
            .collect()
    }

    /// `n` points on the parabola `(x, x², 0)` for `x` in `[-1, 1]`.
    fn parabola_points(n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let x = -1.0 + 2.0 * i as f64 / (n - 1) as f64;
                [x, x * x, 0.0]
            })
            .collect()
    }

    /// `n` points on the unit circle in the xy-plane.
    fn circle_points(n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let th = 2.0 * PI * i as f64 / (n - 1) as f64;
                [th.cos(), th.sin(), 0.0]
            })
            .collect()
    }

    #[test]
    fn fit_line_is_exact() {
        let pts = line_points(20);
        let f = fit_bspline_curve(&pts, 2, 3).unwrap();
        assert!(f.error < 1e-8, "line residual {}", f.error);
        // Poles lie on the line y = 2x + 1, z = x.
        for p in &f.poles {
            assert!((p[1] - (2.0 * p[0] + 1.0)).abs() < 1e-8, "pole {p:?}");
            assert!((p[2] - p[0]).abs() < 1e-8);
        }
    }

    #[test]
    fn fit_parabola() {
        let pts = parabola_points(21);
        let f = fit_bspline_curve(&pts, 3, 13).unwrap();
        assert!(f.error < 1e-3, "parabola residual {}", f.error);
        // The fitted z-coordinate at any sample is close to the input z.
        for r in fit_curve_residuals(&f, &pts) {
            assert!(r < 1e-3);
        }
    }

    #[test]
    fn fit_circle() {
        let pts = circle_points(40);
        let f = fit_bspline_curve(&pts, 3, 12).unwrap();
        assert!(f.error < 0.05, "circle residual {}", f.error);
    }

    #[test]
    fn deg_ge_poles_is_err() {
        let pts = line_points(10);
        assert!(fit_bspline_curve(&pts, 3, 3).is_err(), "num_poles == deg");
        assert!(fit_bspline_curve(&pts, 4, 3).is_err(), "num_poles < deg");
        assert!(fit_bspline_curve(&pts, 0, 3).is_err(), "zero degree");
    }

    #[test]
    fn eval_endpoints_clamped() {
        let pts = line_points(20);
        let f = fit_bspline_curve(&pts, 2, 3).unwrap();
        let a = eval_fit_curve(&f, 0.0);
        let b = eval_fit_curve(&f, 1.0);
        for c in 0..3 {
            assert!((a[c] - pts[0][c]).abs() < 1e-6, "start c{c}");
            assert!((b[c] - pts[19][c]).abs() < 1e-6, "end c{c}");
        }
    }

    #[test]
    fn fit_arc_midpoint() {
        // Quarter circle (1,0,0) → (0,1,0), uniformly sampled in angle. The
        // chord length between consecutive samples is constant, so the
        // chord-length parameter t is linear in the sample index and t = 0.5
        // lands at the θ = 45° sample.
        let n = 25;
        let pts: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let th = PI * 0.5 * i as f64 / (n - 1) as f64;
                [th.cos(), th.sin(), 0.0]
            })
            .collect();
        let f = fit_bspline_curve(&pts, 3, 9).unwrap();
        let mid = eval_fit_curve(&f, 0.5);
        let want = [
            std::f64::consts::FRAC_1_SQRT_2,
            std::f64::consts::FRAC_1_SQRT_2,
            0.0,
        ];
        let d = dist(mid, want);
        assert!(d < 0.05, "arc midpoint error {d}");
    }

    #[test]
    fn partition_of_unity() {
        // All poles equal → the curve is constant; Σ B_i(t) = 1.
        let knots = uniform_clamped_knots(6, 3);
        let f = BSplineFit {
            knots,
            deg: 3,
            poles: vec![[1.0, 1.0, 1.0]; 6],
            weights: None,
            error: 0.0,
        };
        for t in [0.0, 0.13, 0.37, 0.5, 0.83, 1.0] {
            let p = eval_fit_curve(&f, t);
            for c in 0..3 {
                assert!((p[c] - 1.0).abs() < 1e-9, "partition of unity at t={t}: {}", p[c]);
            }
        }
    }

    #[test]
    fn error_field_matches_residuals() {
        let pts = circle_points(25);
        let f = fit_bspline_curve(&pts, 3, 8).unwrap();
        let max_res = fit_curve_max_error(&f, &pts);
        assert!((max_res - f.error).abs() < 1e-12, "{} vs {}", max_res, f.error);
    }

    #[test]
    fn centripetal_fit_works() {
        let pts = circle_points(40);
        let f = fit_bspline_curve_centripetal(&pts, 3, 12).unwrap();
        assert!(f.error < 0.05, "centripetal circle residual {}", f.error);
        assert_eq!(bspline_curve_knot_count(f.poles.len(), f.deg), f.knots.len());
    }

    #[test]
    fn builder_validates_and_derivative_finite() {
        let pts = line_points(10);
        let c = fit_bspline_curve(&pts, 2, 4).unwrap();
        // Rebuild through the public constructor and check the derivative.
        let f = build_bspline_curve(c.knots.clone(), c.deg, c.poles.clone(), None).unwrap();
        assert_eq!(f.knots, c.knots);
        let d = eval_fit_curve_derivative(&f, 0.4, 1e-6);
        assert!(d.iter().all(|x| x.is_finite()), "derivative {d:?}");
        // Validation errors.
        assert!(build_bspline_curve(vec![0.0, 0.0, 1.0], 2, c.poles.clone(), None).is_err());
        assert!(build_bspline_curve(c.knots.clone(), 2, vec![], None).is_err());
    }
}
