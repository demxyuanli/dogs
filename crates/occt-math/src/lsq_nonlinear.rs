//! Gauss-Newton nonlinear least squares (overdetermined systems).
//! Source: `math_GaussLeastSquare` / `math_FunctionSetRoot`.
//!
//! The core iterate is
//! `x_{k+1} = x_k - (JᵀJ + λI)⁻¹ Jᵀ r(x_k)`,
//! where `r` is the residual vector and `J` its Jacobian. The damping term
//! `λI` regularizes the normal matrix when the system is rank deficient and is
//! shrunk on successful steps so that the method behaves like Newton near a
//! solution. An analytic Jacobian can be supplied through [`GaussNewton`];
//! [`GaussNewtonSimple`] and [`gauss_newton_solve_fd`] compute it internally
//! by finite differences, matching the split in OCCT's `math_FunctionSetRoot`
//! (analytic gradient) vs `math_GaussLeastSquare` (numerical).

use crate::matrix_ext::solve_linear;
use crate::{MathMatrix, MathVector};

/// Gauss-Newton driver with an analytic Jacobian.
///
/// `residual(x, r)` fills `r[0..m]` with the residuals; `jacobian(x, j)` fills
/// `j[i][j]` for `i` in `0..m` residuals and `j` in `0..n` parameters.
pub struct GaussNewton<'a> {
    pub residual: &'a dyn Fn(&[f64], &mut [f64]),
    pub jacobian: &'a dyn Fn(&[f64], &mut Vec<Vec<f64>>),
    pub max_iter: usize,
    pub tol: f64,
}

impl Default for GaussNewton<'_> {
    fn default() -> Self {
        Self { residual: &|_, _| {}, jacobian: &|_, _| {}, max_iter: 200, tol: 1e-10 }
    }
}

/// Solve the overdetermined nonlinear least-squares problem
/// `min ||r(x)||²` via damped Gauss-Newton:
/// `x_{k+1} = x_k - (JᵀJ + λI)⁻¹ Jᵀ r`, with `λ` decreased on success and
/// increased on failure (Levenberg-like regularization).
///
/// Returns `(solution, final sum of squared residuals, iterations used)`.
/// Fails with `Err` when the step count is exhausted or no improving step is
/// found (singular Jacobian / stuck iterate).
pub fn gauss_newton_solve(
    gn: &GaussNewton,
    x0: &[f64],
    m: usize,
    n: usize,
) -> Result<(Vec<f64>, f64, usize), String> {
    if m == 0 || n == 0 {
        return Err("gauss_newton_solve: m and n must be positive".to_string());
    }
    if x0.len() != n {
        return Err(format!(
            "gauss_newton_solve: x0 length {} != n {n}",
            x0.len()
        ));
    }
    let mut x = x0.to_vec();
    let mut r = vec![0.0; m];
    (gn.residual)(&x, &mut r);
    let mut cost = r.iter().map(|v| v * v).sum::<f64>();
    if gn.max_iter == 0 {
        return Ok((x, cost, 0));
    }
    let mut lambda = 1e-3;

    for iter in 0..gn.max_iter {
        // Jacobian J (m×n).
        let mut j = vec![vec![0.0; n]; m];
        (gn.jacobian)(&x, &mut j);
        if j.len() != m || j.first().map_or(false, |row| row.len() != n) {
            return Err("gauss_newton_solve: jacobian returned wrong dimensions".to_string());
        }

        // Normal matrix JᵀJ and gradient Jᵀr.
        let mut jtj = vec![vec![0.0; n]; n];
        let mut jtr = vec![0.0; n];
        for a in 0..n {
            for b in 0..n {
                let mut s = 0.0;
                for i in 0..m {
                    s += j[i][a] * j[i][b];
                }
                jtj[a][b] = s;
            }
            for i in 0..m {
                jtr[a] += j[i][a] * r[i];
            }
        }
        let grad_max = jtr.iter().fold(0.0_f64, |acc, &v| acc.max(v.abs()));
        if grad_max < gn.tol {
            return Ok((x, cost, iter));
        }

        // Damped Gauss-Newton step: solve (JᵀJ + λI) δ = -Jᵀr.
        let mut accepted = false;
        for _ in 0..40 {
            let mut a = MathMatrix::new(1, n, 1, n);
            for rr in 1..=n {
                for c in 1..=n {
                    let mut v = jtj[rr - 1][c - 1];
                    if rr == c {
                        v += lambda;
                    }
                    a.set_value(rr, c, v);
                }
            }
            let mut bv = MathVector::new(1, n);
            for k in 1..=n {
                bv.set_value(k, -jtr[k - 1]);
            }
            let delta = solve_linear(&a, &bv)
                .map_err(|e| format!("gauss_newton_solve: singular system: {e}"))?;

            let mut xt = x.clone();
            for k in 0..n {
                xt[k] += delta.value(k + 1);
            }
            let mut rt = vec![0.0; m];
            (gn.residual)(&xt, &mut rt);
            let cost_t = rt.iter().map(|v| v * v).sum::<f64>();
            if cost_t < cost {
                x = xt;
                r = rt;
                cost = cost_t;
                lambda = (lambda / 10.0).max(1e-12);
                accepted = true;
                let max_step = (1..=n).map(|i| delta.value(i).abs()).fold(0.0_f64, f64::max);
                if max_step < gn.tol {
                    return Ok((x, cost, iter + 1));
                }
                break;
            }
            lambda *= 10.0;
            if lambda > 1e15 {
                break;
            }
        }
        if !accepted {
            return Err(
                "gauss_newton_solve: no improving step found (singular Jacobian or stuck)"
                    .to_string(),
            );
        }
    }
    Err(format!(
        "gauss_newton_solve: max iterations ({}) reached",
        gn.max_iter
    ))
}

/// Central-difference Jacobian of `f(x)` (returns the residual vector) at `x`.
/// `m` residual components, `n` parameters, step `h`.
pub fn jacobian_finite_difference(
    f: &dyn Fn(&[f64], &mut [f64]),
    x: &[f64],
    m: usize,
    n: usize,
    h: f64,
) -> Vec<Vec<f64>> {
    let h = if h == 0.0 { 1e-6 } else { h };
    let mut j = vec![vec![0.0; n]; m];
    let mut xp = x.to_vec();
    let mut xm = x.to_vec();
    let mut fp = vec![0.0; m];
    let mut fm = vec![0.0; m];
    for col in 0..n {
        xp[col] = x[col] + h;
        xm[col] = x[col] - h;
        f(&xp, &mut fp);
        f(&xm, &mut fm);
        for i in 0..m {
            j[i][col] = (fp[i] - fm[i]) / (2.0 * h);
        }
        xp[col] = x[col];
        xm[col] = x[col];
    }
    j
}

/// Forward-difference Jacobian of `f(x)` at `x`: uses `n` residual evaluations
/// (plus one base evaluation) instead of `2n` for central differences.
/// Less accurate than [`jacobian_finite_difference`] but cheaper; prefer it
/// when residual evaluations dominate the cost.
pub fn jacobian_forward_difference(
    f: &dyn Fn(&[f64], &mut [f64]),
    x: &[f64],
    m: usize,
    n: usize,
    h: f64,
) -> Vec<Vec<f64>> {
    let h = if h == 0.0 { 1e-7 } else { h };
    let mut j = vec![vec![0.0; n]; m];
    let mut xp = x.to_vec();
    let mut f0 = vec![0.0; m];
    let mut fp = vec![0.0; m];
    f(x, &mut f0);
    for col in 0..n {
        xp[col] = x[col] + h;
        f(&xp, &mut fp);
        for i in 0..m {
            j[i][col] = (fp[i] - f0[i]) / h;
        }
        xp[col] = x[col];
    }
    j
}

/// Approximate parameter standard errors from the covariance matrix
/// `(JᵀJ)⁻¹ · cost/(m - n)` at a converged fit. Requires `m > n` residuals.
///
/// The returned vector holds `sqrt(diag((JᵀJ)⁻¹) · s²)` where
/// `s² = ||r(x)||²/(m - n)` is the residual variance estimate, matching the
/// standard nonlinear least-squares uncertainty formula.
pub fn estimate_parameter_uncertainties(
    gn: &GaussNewton,
    x: &[f64],
    m: usize,
    n: usize,
) -> Result<Vec<f64>, String> {
    if m <= n {
        return Err("estimate_parameter_uncertainties: need m > n (overdetermined)".to_string());
    }
    let mut j = vec![vec![0.0; n]; m];
    (gn.jacobian)(x, &mut j);
    let mut jtj = vec![vec![0.0; n]; n];
    for a in 0..n {
        for b in 0..n {
            let mut s = 0.0;
            for i in 0..m {
                s += j[i][a] * j[i][b];
            }
            jtj[a][b] = s;
        }
    }
    let mut r = vec![0.0; m];
    (gn.residual)(x, &mut r);
    let cost = r.iter().map(|v| v * v).sum::<f64>();
    let sigma2 = cost / (m - n) as f64;

    let mut a = MathMatrix::new(1, n, 1, n);
    for i in 1..=n {
        for k in 1..=n {
            a.set_value(i, k, jtj[i - 1][k - 1]);
        }
    }
    let inv = crate::matrix_ext::inverse(&a)
        .map_err(|e| format!("estimate_parameter_uncertainties: {e}"))?;
    let mut std_errs = vec![0.0; n];
    for k in 1..=n {
        std_errs[k - 1] = (sigma2 * inv.value(k, k)).max(0.0).sqrt();
    }
    Ok(std_errs)
}

/// Sum of squared residuals `Σ r_i(x)²` for an `m`-component residual at `x`.
///
/// This is the objective minimized by [`gauss_newton_solve`]; exposing it as a
/// helper lets callers evaluate the cost independently of the solver loop
/// (e.g. for reporting or for comparing solver runs).
pub fn residual_sse(
    residual: &dyn Fn(&[f64], &mut [f64]),
    x: &[f64],
    m: usize,
) -> f64 {
    let mut r = vec![0.0; m];
    residual(x, &mut r);
    r.iter().map(|v| v * v).sum()
}

/// Gauss-Newton driver that computes the Jacobian internally by central
/// finite differences (`n = x0.len()`).
pub struct GaussNewtonSimple<'a> {
    pub residual: &'a dyn Fn(&[f64], &mut [f64]),
    pub max_iter: usize,
    pub tol: f64,
}

impl Default for GaussNewtonSimple<'_> {
    fn default() -> Self {
        Self { residual: &|_, _| {}, max_iter: 200, tol: 1e-10 }
    }
}

/// Solve `min ||r(x)||²` where the Jacobian is approximated by finite
/// differences. Returns `(solution, final sum of squared residuals, iterations)`.
pub fn gauss_newton_solve_fd(
    gn: &GaussNewtonSimple,
    x0: &[f64],
    m: usize,
) -> Result<(Vec<f64>, f64, usize), String> {
    let n = x0.len();
    let jac = |x: &[f64], j: &mut Vec<Vec<f64>>| {
        *j = jacobian_finite_difference(gn.residual, x, m, n, 1e-6);
    };
    let full = GaussNewton {
        residual: gn.residual,
        jacobian: &jac,
        max_iter: gn.max_iter,
        tol: gn.tol,
    };
    gauss_newton_solve(&full, x0, m, n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qr_full::qr_least_squares;

    const TAU: f64 = std::f64::consts::TAU;

    fn build_matrix(rows: &[Vec<f64>]) -> MathMatrix {
        let m = rows.len();
        let n = rows[0].len();
        let mut a = MathMatrix::new(1, m, 1, n);
        for i in 0..m {
            for j in 0..n {
                a.set_value(i + 1, j + 1, rows[i][j]);
            }
        }
        a
    }

    #[test]
    fn fit_circle_center_and_radius() {
        let (cx, cy, rad) = (0.5, 1.0, 2.0);
        let mut pts = Vec::new();
        for k in 0..16 {
            let a = TAU * k as f64 / 16.0;
            pts.push((cx + rad * a.cos(), cy + rad * a.sin()));
        }
        let residual = |p: &[f64], r: &mut [f64]| {
            for (i, (x, y)) in pts.iter().enumerate() {
                r[i] = ((x - p[0]).powi(2) + (y - p[1]).powi(2)).sqrt() - p[2];
            }
        };
        let gn = GaussNewtonSimple { residual: &residual, max_iter: 200, tol: 1e-10 };
        let (sol, cost, _) = gauss_newton_solve_fd(&gn, &[0.0, 0.0, 1.0], pts.len()).unwrap();
        assert!((sol[0] - cx).abs() < 1e-3, "cx = {}", sol[0]);
        assert!((sol[1] - cy).abs() < 1e-3, "cy = {}", sol[1]);
        assert!((sol[2] - rad).abs() < 1e-3, "r = {}", sol[2]);
        assert!(cost < 1e-6, "residual cost {cost}");
    }

    #[test]
    fn fit_exponential_curve() {
        let xs: Vec<f64> = (0..24).map(|i| 0.2 * i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 * (0.7 * x).exp()).collect();
        let residual = |p: &[f64], r: &mut [f64]| {
            for i in 0..xs.len() {
                r[i] = p[0] * (p[1] * xs[i]).exp() - ys[i];
            }
        };
        let gn = GaussNewtonSimple { residual: &residual, max_iter: 200, tol: 1e-10 };
        let (sol, cost, _) = gauss_newton_solve_fd(&gn, &[1.0, 1.0], xs.len()).unwrap();
        assert!((sol[0] - 2.0).abs() < 1e-3, "a = {}", sol[0]);
        assert!((sol[1] - 0.7).abs() < 1e-3, "b = {}", sol[1]);
        assert!(cost < 1e-8, "cost {cost}");
    }

    #[test]
    fn rosenbrock_residual_norm_shrinks() {
        // r = [10(x2 - x1²), 1 - x1]  ->  zero at (1, 1).
        let residual = |p: &[f64], r: &mut [f64]| {
            r[0] = 10.0 * (p[1] - p[0] * p[0]);
            r[1] = 1.0 - p[0];
        };
        let gn = GaussNewtonSimple { residual: &residual, max_iter: 400, tol: 1e-12 };
        let (sol, cost, _) = gauss_newton_solve_fd(&gn, &[-1.2, 1.0], 2).unwrap();
        let initial = 100.0 * (1.0 - 1.44_f64).powi(2) + 1.44_f64;
        let _ = initial;
        assert!(cost < 1e-8, "cost {cost}");
        assert!((sol[0] - 1.0).abs() < 1e-4, "x1 = {}", sol[0]);
        assert!((sol[1] - 1.0).abs() < 1e-4, "x2 = {}", sol[1]);
    }

    #[test]
    fn overdetermined_linear_matches_ls() {
        let a_data = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 10.0],
            vec![1.0, 1.0, 1.0],
            vec![2.0, 0.0, 1.0],
        ];
        let b_data = vec![6.0, 15.0, 25.0, 3.0, 3.0]; // A·[1,1,1]
        let a = build_matrix(&a_data);
        let b = crate::MathVector::from_slice(&b_data);
        let ls = qr_least_squares(&a, &b);

        let residual = |p: &[f64], r: &mut [f64]| {
            for i in 0..5 {
                let mut s = 0.0;
                for j in 0..3 {
                    s += a_data[i][j] * p[j];
                }
                r[i] = s - b_data[i];
            }
        };
        let jac = |_p: &[f64], j: &mut Vec<Vec<f64>>| {
            *j = a_data.clone();
        };
        let gn = GaussNewton { residual: &residual, jacobian: &jac, max_iter: 50, tol: 1e-12 };
        let (sol, cost, _) = gauss_newton_solve(&gn, &[0.0, 0.0, 0.0], 5, 3).unwrap();
        for j in 0..3 {
            assert!((sol[j] - ls.value(j + 1)).abs() < 1e-6, "x{j} = {}", sol[j]);
        }
        assert!(cost < 1e-10, "cost {cost}");
    }

    #[test]
    fn fd_jacobian_matches_analytic() {
        let f = |p: &[f64], r: &mut [f64]| {
            r[0] = p[0] * p[0] + p[1];
            r[1] = (p[0] * p[1]).sin();
            r[2] = p[0] * p[1] * p[2];
        };
        let x = [0.7, 1.2, -0.5];
        let jfd = jacobian_finite_difference(&f, &x, 3, 3, 1e-6);
        let c = (x[0] * x[1]).cos();
        let jana = [
            [2.0 * x[0], 1.0, 0.0],
            [x[1] * c, x[0] * c, 0.0],
            [x[1] * x[2], x[0] * x[2], x[0] * x[1]],
        ];
        for i in 0..3 {
            for j in 0..3 {
                assert!((jfd[i][j] - jana[i][j]).abs() < 1e-4, "J[{i}][{j}]");
            }
        }
    }

    #[test]
    fn fd_jacobian_forward_matches_analytic() {
        let f = |p: &[f64], r: &mut [f64]| {
            r[0] = p[0] * p[0] + p[1];
            r[1] = (p[0] * p[1]).sin();
            r[2] = p[0] * p[1] * p[2];
        };
        let x = [0.7, 1.2, -0.5];
        let jfd = jacobian_forward_difference(&f, &x, 3, 3, 1e-7);
        let c = (x[0] * x[1]).cos();
        let jana = [
            [2.0 * x[0], 1.0, 0.0],
            [x[1] * c, x[0] * c, 0.0],
            [x[1] * x[2], x[0] * x[2], x[0] * x[1]],
        ];
        for i in 0..3 {
            for j in 0..3 {
                assert!((jfd[i][j] - jana[i][j]).abs() < 1e-4, "J[{i}][{j}]");
            }
        }
    }

    #[test]
    fn parameter_uncertainties_positive() {
        // Fit y = 1 + 2x to slightly noisy points: standard errors > 0.
        let xs: Vec<f64> = (0..10).map(|i| 0.5 * i as f64).collect();
        let ys: Vec<f64> = xs
            .iter()
            .enumerate()
            .map(|(i, &x)| 1.0 + 2.0 * x + 0.02 * (i % 3) as f64)
            .collect();
        let residual = |p: &[f64], r: &mut [f64]| {
            for i in 0..xs.len() {
                r[i] = p[0] + p[1] * xs[i] - ys[i];
            }
        };
        let jac = |_p: &[f64], j: &mut Vec<Vec<f64>>| {
            for i in 0..xs.len() {
                j[i][0] = 1.0;
                j[i][1] = xs[i];
            }
        };
        let gn = GaussNewton { residual: &residual, jacobian: &jac, max_iter: 50, tol: 1e-12 };
        let (sol, _, _) = gauss_newton_solve(&gn, &[0.0, 0.0], xs.len(), 2).unwrap();
        assert!((sol[0] - 1.0).abs() < 0.1 && (sol[1] - 2.0).abs() < 0.1);
        let unc = estimate_parameter_uncertainties(&gn, &sol, xs.len(), 2).unwrap();
        assert!(unc[0] > 0.0 && unc[1] > 0.0, "uncertainties {unc:?}");
        assert!(unc[0] < 1.0 && unc[1] < 1.0);
    }

    #[test]
    fn zero_residual_for_satisfiable_system() {
        // 3 exact points on a parabola -> residual zero achievable (m = n = 3).
        let xs = [0.0, 1.0, 2.0];
        let ys = [1.0, 3.0, 7.0]; // y = 1 + x + x²
        let residual = |p: &[f64], r: &mut [f64]| {
            for i in 0..3 {
                r[i] = p[0] + p[1] * xs[i] + p[2] * xs[i] * xs[i] - ys[i];
            }
        };
        let gn = GaussNewtonSimple { residual: &residual, max_iter: 50, tol: 1e-12 };
        let (sol, cost, _) = gauss_newton_solve_fd(&gn, &[0.0, 0.0, 0.0], 3).unwrap();
        assert!(cost < 1e-16, "cost {cost}");
        assert!((sol[0] - 1.0).abs() < 1e-8 && (sol[1] - 1.0).abs() < 1e-8 && (sol[2] - 1.0).abs() < 1e-8);
    }

    #[test]
    fn max_iter_exceeded_is_error() {
        let residual = |p: &[f64], r: &mut [f64]| {
            r[0] = 10.0 * (p[1] - p[0] * p[0]);
            r[1] = 1.0 - p[0];
        };
        let gn = GaussNewtonSimple { residual: &residual, max_iter: 2, tol: 1e-12 };
        let res = gauss_newton_solve_fd(&gn, &[-1.2, 1.0], 2);
        assert!(res.is_err(), "two iterations cannot solve Rosenbrock");
    }

    #[test]
    fn square_system_reduces_to_newton() {
        // f(x,y) = [x²+y²-4, x-y] -> roots at (±√2, ±√2).
        let residual = |p: &[f64], r: &mut [f64]| {
            r[0] = p[0] * p[0] + p[1] * p[1] - 4.0;
            r[1] = p[0] - p[1];
        };
        let gn = GaussNewtonSimple { residual: &residual, max_iter: 100, tol: 1e-12 };
        let (sol, cost, _) = gauss_newton_solve_fd(&gn, &[1.0, 2.0], 2).unwrap();
        let root = 2.0_f64.sqrt();
        assert!((sol[0] - root).abs() < 1e-6, "x = {}", sol[0]);
        assert!((sol[1] - root).abs() < 1e-6, "y = {}", sol[1]);
        assert!(cost < 1e-12, "cost {cost}");
    }

    #[test]
    fn loose_tolerance_stops_early() {
        let xs: Vec<f64> = (0..24).map(|i| 0.2 * i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 * (0.7 * x).exp()).collect();
        let residual = |p: &[f64], r: &mut [f64]| {
            for i in 0..xs.len() {
                r[i] = p[0] * (p[1] * xs[i]).exp() - ys[i];
            }
        };
        let loose = GaussNewtonSimple { residual: &residual, max_iter: 200, tol: 1e-2 };
        let tight = GaussNewtonSimple { residual: &residual, max_iter: 200, tol: 1e-10 };
        let (_, cost_loose, iters_loose) = gauss_newton_solve_fd(&loose, &[1.0, 1.0], xs.len()).unwrap();
        let (_, cost_tight, iters_tight) = gauss_newton_solve_fd(&tight, &[1.0, 1.0], xs.len()).unwrap();
        assert!(iters_loose <= iters_tight, "{iters_loose} > {iters_tight}");
        assert!(cost_loose >= cost_tight, "loose tol should stop with >= cost");
    }

    #[test]
    fn scale_invariance_sanity() {
        let xs: Vec<f64> = (0..24).map(|i| 0.2 * i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 * (0.7 * x).exp()).collect();
        let base = |p: &[f64], r: &mut [f64]| {
            for i in 0..xs.len() {
                r[i] = p[0] * (p[1] * xs[i]).exp() - ys[i];
            }
        };
        let scaled = |p: &[f64], r: &mut [f64]| {
            for i in 0..xs.len() {
                r[i] = 100.0 * (p[0] * (p[1] * xs[i]).exp() - ys[i]);
            }
        };
        let g1 = GaussNewtonSimple { residual: &base, max_iter: 200, tol: 1e-10 };
        let g2 = GaussNewtonSimple { residual: &scaled, max_iter: 200, tol: 1e-10 };
        let (s1, _, _) = gauss_newton_solve_fd(&g1, &[1.0, 1.0], xs.len()).unwrap();
        let (s2, _, _) = gauss_newton_solve_fd(&g2, &[1.0, 1.0], xs.len()).unwrap();
        assert!((s1[0] - s2[0]).abs() < 1e-3, "a: {} vs {}", s1[0], s2[0]);
        assert!((s1[1] - s2[1]).abs() < 1e-3, "b: {} vs {}", s1[1], s2[1]);
    }
}
