//! Constrained optimization via the penalty and augmented-Lagrangian methods.
//! Source: `math_...` / classical penalty methods (Bertsekas, "Constrained
//! Optimization and Lagrange Multiplier Methods"; Nocedal & Wright, ch. 17).
//!
//! # Problem form
//!
//! The modules solve the general smooth constrained minimization problem
//!
//! ```text
//! min f(x)   s.t.   g_i(x) ≤ 0   (i in I, inequality)
//!                   g_j(x) = 0   (j in E, equality)
//! ```
//!
//! by replacing it with a sequence of *unconstrained* problems.  Each
//! constraint contributes a term that is zero in the feasible region and grows
//! smoothly once the point leaves it, so an unconstrained minimizer is pulled
//! back toward feasibility.
//!
//! # Quadratic penalty method
//!
//! The penalty function for a factor μ is
//!
//! ```text
//! F_μ(x) = f(x) + μ·Σ_{i∈I} max(0, g_i(x))² + μ·Σ_{j∈E} g_j(x)²
//! ```
//!
//! and [`penalty_optimize`] minimizes `F_μ` for the increasing schedule
//! `μ = 1, 10, 100, 1000`, warming each inner solve with the previous iterate
//! (a continuation in μ).  The unconstrained minimizers `x*(μ)` converge to
//! the true constrained solution as `μ → ∞`; the residual infeasibility is
//! `O(1/μ)` and the objective error `O(1/μ)` as well.
//!
//! # Augmented Lagrangian method
//!
//! The Hestenes–Powell–Rockafellar [`lagrange_penalty_optimize`] adds a
//! Lagrange-multiplier estimate λ for every constraint:
//!
//! ```text
//! L(x, λ, μ) = f(x) + Σ_{j∈E} [ λ_j·g_j(x) + (μ/2)·g_j(x)² ]
//!                  + Σ_{i∈I} [ max(0, λ_i + μ·g_i(x))² − λ_i² ] / (2μ)
//! ```
//!
//! After each inner solve the multipliers are updated with the current
//! constraint values and μ grows by a factor of ten.  Because the multiplier
//! estimates absorb the "force" needed to satisfy the constraints, a moderate μ
//! reaches much tighter feasibility than the plain quadratic penalty — an
//! equality constraint converges to machine-visible feasibility within a
//! handful of outer iterations.
//!
//! # Constraint convention
//!
//! A [`ConstraintFn`] stores a scalar `g` and a flag.  For an inequality the
//! feasible set is `{ x : g(x) ≤ 0 }`; for an equality it is `{ x : g(x) = 0 }`.
//! [`constraint_violation`] reports the largest violation and [`feasible`]
//! tests it against a tolerance.

use crate::optimize::nelder_mead;

/// A single scalar constraint `g(x)` with its type tag.
///
/// The convention matches OCCT's `math_...` minimizers: an *inequality*
/// constraint is satisfied when `g(x) ≤ 0`; an *equality* constraint is
/// satisfied when `g(x) = 0`.
pub struct ConstraintFn<'a> {
    pub g: &'a dyn Fn(&[f64]) -> f64,
    pub is_inequality: bool,
}

/// A constrained minimization problem: scalar objective plus a constraint set.
///
/// `max_iter` budgets the inner unconstrained solves (spread across the outer
/// schedule) and `tol` is used both as the inner-solver convergence tolerance
/// and as the feasibility threshold for the augmented-Lagrangian stopping test.
pub struct ConstraintOpt<'a> {
    pub objective: &'a dyn Fn(&[f64]) -> f64,
    pub constraints: &'a [ConstraintFn<'a>],
    pub max_iter: usize,
    pub tol: f64,
}

/// Maximum signed violation of one constraint at `x`.
///
/// Inequality constraints violate by `max(0, g(x))` (the amount by which
/// `g(x)` exceeds its 0 upper bound); equality constraints violate by
/// `|g(x)|`.  The constraint set is violated when the maximum over all
/// constraints is positive.
fn single_violation(c: &ConstraintFn, x: &[f64]) -> f64 {
    let g = (c.g)(x);
    if c.is_inequality {
        g.max(0.0)
    } else {
        g.abs()
    }
}

/// Largest constraint violation over the whole set (0 when feasible).
///
/// The violation of a single constraint is `max(0, g(x))` for an inequality
/// and `|g(x)|` for an equality; the returned value is the maximum over all
/// constraints.  It is exactly zero on the feasible set and grows smoothly
/// outside it, which makes it a convenient scalar measure of infeasibility
/// (e.g. for a convergence monitor or for comparing candidate solutions).
pub fn constraint_violation(opt: &ConstraintOpt, x: &[f64]) -> f64 {
    let mut maxv: f64 = 0.0;
    for c in opt.constraints {
        maxv = maxv.max(single_violation(c, x));
    }
    maxv
}

/// Whether `x` satisfies every constraint within tolerance `tol`.
///
/// Equivalently `constraint_violation(opt, x) <= tol`.  The tolerance absorbs
/// the residual `O(1/μ)` infeasibility left by the penalty method so callers
/// can test near-feasible iterates without demanding exact satisfaction.
pub fn feasible(opt: &ConstraintOpt, x: &[f64], tol: f64) -> bool {
    constraint_violation(opt, x) <= tol
}

/// Shared argument validation for the penalty drivers.
fn validate(opt: &ConstraintOpt, x0: &[f64]) -> Result<(), String> {
    if x0.is_empty() {
        return Err("penalty: empty initial point".to_string());
    }
    if !x0.iter().all(|v| v.is_finite()) {
        return Err("penalty: non-finite initial point".to_string());
    }
    if opt.max_iter == 0 {
        return Err("penalty: max_iter must be at least 1".to_string());
    }
    if !opt.tol.is_finite() || opt.tol < 0.0 {
        return Err("penalty: tolerance must be a non-negative finite value".to_string());
    }
    Ok(())
}

/// Default penalty factor schedule.
const DEFAULT_MU: [f64; 4] = [1.0, 10.0, 100.0, 1000.0];

/// The quadratic penalty function `F_μ(x)` minimized by the penalty method.
///
/// ```text
/// F_μ(x) = f(x) + μ·Σ_ineq max(0, g_i(x))² + μ·Σ_eq g_j(x)²
/// ```
///
/// This is a useful diagnostic: for a feasible point the penalty terms vanish
/// and `F_μ(x) = f(x)`; the gap `F_μ(x) − f(x)` is exactly the weighted
/// squared violation.  Callers can inspect it to monitor how far an iterate is
/// from feasibility, or to compare the penalty strength across μ values.
pub fn penalty_value(opt: &ConstraintOpt, x: &[f64], mu: f64) -> f64 {
    let mut v = (opt.objective)(x);
    for c in opt.constraints {
        let g = (c.g)(x);
        let viol = if c.is_inequality { g.max(0.0) } else { g };
        v += mu * viol * viol;
    }
    v
}

/// Penalty method with a caller-supplied μ schedule.
///
/// Each outer iteration minimizes the penalized objective [`penalty_value`]
///
/// ```text
/// F_μ(x) = f(x) + μ·Σ_ineq max(0, g_i(x))² + μ·Σ_eq g_j(x)²
/// ```
///
/// starting from the previous iterate (a continuation / homotopy in μ).  As μ
/// grows the unconstrained minimizers of `F_μ` approach the feasible region,
/// so the final iterate satisfies the constraints to within `O(1/μ)`.  The
/// inner unconstrained problems are solved with [`nelder_mead`]; each one is
/// budgeted `max(50, max_iter / len(schedule))` iterations.
///
/// The caller controls the schedule, which is useful for experimenting with
/// the μ–accuracy trade-off: a short schedule is cheaper but leaves a larger
/// residual violation, while a long schedule (e.g. up to 1e6) drives the
/// solution arbitrarily close to feasibility at the cost of a steeper inner
/// objective.  [`penalty_optimize`] wraps this with the default schedule.
///
/// # Returns
/// A tuple `(solution, objective at solution, total inner iterations)`.  The
/// reported objective is the *original* `f(x)`, not the penalized `F_μ(x)`.
///
/// # Errors
/// See [`validate`]; an empty `mu_schedule` is also an error.
pub fn penalty_optimize_with_mu(
    opt: &ConstraintOpt,
    x0: &[f64],
    mu_schedule: &[f64],
) -> Result<(Vec<f64>, f64, usize), String> {
    validate(opt, x0)?;
    if mu_schedule.is_empty() {
        return Err("penalty: empty mu schedule".to_string());
    }
    let inner_budget = (opt.max_iter / mu_schedule.len()).max(50);
    let mut x = x0.to_vec();
    let mut iterations: usize = 0;
    for &mu in mu_schedule {
        let penalty = |p: &[f64]| penalty_value(opt, p, mu);
        x = nelder_mead(&penalty, &x, inner_budget, opt.tol);
        iterations += inner_budget;
    }
    let fx = (opt.objective)(&x);
    Ok((x, fx, iterations))
}

/// Quadratic penalty method for constrained minimization.
///
/// Minimizes `f(x) + μ·Σ(max(0,g_i))² + μ·Σ(g_j)²` for the default μ schedule
/// `[1, 10, 100, 1000]` (see [`penalty_optimize_with_mu`] for a custom one).
///
/// # Returns
/// A tuple `(solution, objective at solution, total inner iterations)`.  The
/// solution is the final unconstrained minimizer of the penalty function for
/// the largest μ; the objective is the original `f(x)` evaluated there.
///
/// # Errors
/// Empty or non-finite `x0`, `max_iter == 0`, or a negative tolerance produce
/// an `Err`.
pub fn penalty_optimize(opt: &ConstraintOpt, x0: &[f64]) -> Result<(Vec<f64>, f64, usize), String> {
    penalty_optimize_with_mu(opt, x0, &DEFAULT_MU)
}

/// Augmented-Lagrangian (Hestenes–Powell–Rockafellar) constrained minimization.
///
/// The outer iteration minimizes
///
/// ```text
/// L(x, λ, μ) = f(x)
///            + Σ_eq [ λ_j·g_j(x) + (μ/2)·g_j(x)² ]
///            + Σ_ineq [ max(0, λ_j + μ·g_i(x))² − λ_j² ] / (2μ)
/// ```
///
/// then updates the multipliers
///
/// ```text
/// λ_j ← λ_j + μ·g_j(x)            (equality)
/// λ_i ← max(0, λ_i + μ·g_i(x))    (inequality)
/// ```
///
/// and grows μ tenfold each outer pass.  Convergence is declared when the
/// constraint violation drops below `opt.tol`.  Because the multiplier
/// estimates λ track the constraint gradients, a fixed, moderate μ converges
/// much faster than the plain penalty method (the equality test converges to
/// machine-visible feasibility within a handful of outer iterations).
///
/// # Returns
/// A tuple `(solution, objective at solution, total inner iterations)`.  The
/// objective is the original `f(x)` evaluated at the returned solution.
///
/// # Errors
/// Same validation as [`penalty_optimize`].
pub fn lagrange_penalty_optimize(
    opt: &ConstraintOpt,
    x0: &[f64],
) -> Result<(Vec<f64>, f64, usize), String> {
    validate(opt, x0)?;
    let nc = opt.constraints.len();
    let mut lambda = vec![0.0; nc];
    let mut mu: f64 = 1.0;
    let mut x = x0.to_vec();
    let mut iterations: usize = 0;
    let max_outer = 20;
    let inner_budget = (opt.max_iter / max_outer).max(50);
    for _ in 0..max_outer {
        let lam = &lambda;
        let lag = |p: &[f64]| -> f64 {
            let mut v = (opt.objective)(p);
            for (j, c) in opt.constraints.iter().enumerate() {
                let g = (c.g)(p);
                if c.is_inequality {
                    let s = lam[j] + mu * g;
                    v += (s.max(0.0).powi(2) - lam[j] * lam[j]) / (2.0 * mu);
                } else {
                    v += lam[j] * g + 0.5 * mu * g * g;
                }
            }
            v
        };
        x = nelder_mead(&lag, &x, inner_budget, opt.tol);
        iterations += inner_budget;
        // Multiplier update + convergence check.
        let mut max_viol: f64 = 0.0;
        for (j, c) in opt.constraints.iter().enumerate() {
            let g = (c.g)(&x);
            if c.is_inequality {
                lambda[j] = (lambda[j] + mu * g).max(0.0);
                max_viol = max_viol.max(g.max(0.0));
            } else {
                lambda[j] += mu * g;
                max_viol = max_viol.max(g.abs());
            }
        }
        mu *= 10.0;
        if max_viol < opt.tol {
            break;
        }
    }
    let fx = (opt.objective)(&x);
    Ok((x, fx, iterations))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn eq_opt<'a>(
        obj: &'a dyn Fn(&[f64]) -> f64,
        cons: &'a [ConstraintFn<'a>],
    ) -> ConstraintOpt<'a> {
        ConstraintOpt { objective: obj, constraints: cons, max_iter: 5000, tol: 1e-10 }
    }

    #[test]
    fn equality_line_minimizes_quadratic() {
        // min x²+y²  s.t.  x + y − 1 = 0 → (0.5, 0.5), f = 0.5.
        let obj = |x: &[f64]| x[0] * x[0] + x[1] * x[1];
        let g = |x: &[f64]| x[0] + x[1] - 1.0;
        let cons = [ConstraintFn { g: &g, is_inequality: false }];
        let opt = eq_opt(&obj, &cons);
        let (x, fx, _) = penalty_optimize(&opt, &[0.2, 0.2]).unwrap();
        assert!((x[0] - 0.5).abs() < 1e-2, "x0 = {}", x[0]);
        assert!((x[1] - 0.5).abs() < 1e-2, "x1 = {}", x[1]);
        assert!((fx - 0.5).abs() < 2e-2, "objective = {fx}");
        assert!(feasible(&opt, &x, 1e-3));
    }

    #[test]
    fn inequality_lower_bound() {
        // min x²  s.t.  x ≥ 2  →  x ≈ 2, f ≈ 4.  Constraint written as g = 2−x ≤ 0.
        let obj = |x: &[f64]| x[0] * x[0];
        let g = |x: &[f64]| 2.0 - x[0];
        let cons = [ConstraintFn { g: &g, is_inequality: true }];
        let opt = eq_opt(&obj, &cons);
        let (x, fx, _) = penalty_optimize(&opt, &[0.0]).unwrap();
        assert!(x[0] > 1.98 && x[0] < 2.02, "x = {}", x[0]);
        assert!((fx - 4.0).abs() < 0.08, "objective = {fx}");
    }

    #[test]
    fn boxed_inequalities_origin() {
        // min x²+y²  s.t.  x ≥ 0, y ≥ 0, x + y ≤ 1  →  (0, 0).
        let obj = |x: &[f64]| x[0] * x[0] + x[1] * x[1];
        let g1 = |x: &[f64]| -x[0];
        let g2 = |x: &[f64]| -x[1];
        let g3 = |x: &[f64]| x[0] + x[1] - 1.0;
        let cons = [
            ConstraintFn { g: &g1, is_inequality: true },
            ConstraintFn { g: &g2, is_inequality: true },
            ConstraintFn { g: &g3, is_inequality: true },
        ];
        let opt = eq_opt(&obj, &cons);
        let (x, fx, _) = penalty_optimize(&opt, &[0.5, 0.5]).unwrap();
        assert!(x[0].abs() < 1e-3, "x0 = {}", x[0]);
        assert!(x[1].abs() < 1e-3, "x1 = {}", x[1]);
        assert!(fx < 1e-6, "objective = {fx}");
    }

    #[test]
    fn circle_equality_minimizes_x() {
        // min x  s.t.  x² + y² = 1  →  (−1, 0).
        let obj = |x: &[f64]| x[0];
        let g = |x: &[f64]| x[0] * x[0] + x[1] * x[1] - 1.0;
        let cons = [ConstraintFn { g: &g, is_inequality: false }];
        let opt = eq_opt(&obj, &cons);
        let (x, fx, _) = penalty_optimize(&opt, &[-0.5, 0.5]).unwrap();
        assert!((x[0] + 1.0).abs() < 0.05, "x0 = {}", x[0]);
        assert!(x[1].abs() < 0.05, "x1 = {}", x[1]);
        assert!((fx + 1.0).abs() < 0.05, "objective = {fx}");
    }

    #[test]
    fn feasible_true_and_false() {
        let g = |x: &[f64]| 2.0 - x[0];
        let cons = [ConstraintFn { g: &g, is_inequality: true }];
        let obj = |x: &[f64]| x[0] * x[0];
        let opt = eq_opt(&obj, &cons);
        assert!(feasible(&opt, &[3.0], 1e-6), "x=3 satisfies x ≥ 2");
        assert!(!feasible(&opt, &[1.5], 1e-6), "x=1.5 violates x ≥ 2");
        // Equality version: |g| ≤ tol.
        let ge = |x: &[f64]| x[0] + x[1] - 1.0;
        let cons_e = [ConstraintFn { g: &ge, is_inequality: false }];
        let opt_e = eq_opt(&obj, &cons_e);
        assert!(feasible(&opt_e, &[0.5, 0.5], 1e-9));
        assert!(!feasible(&opt_e, &[0.5, 0.6], 1e-9));
    }

    #[test]
    fn violation_measures_infeasibility() {
        let g = |x: &[f64]| 2.0 - x[0];
        let cons = [ConstraintFn { g: &g, is_inequality: true }];
        let obj = |x: &[f64]| x[0] * x[0];
        let opt = eq_opt(&obj, &cons);
        assert!(constraint_violation(&opt, &[1.5]) > 0.0);
        assert!(constraint_violation(&opt, &[1.5]) - 0.5 < 1e-12);
        assert!(constraint_violation(&opt, &[3.0]) < 1e-12);
    }

    #[test]
    fn penalty_solution_improves_with_mu() {
        // Stronger μ schedule → tighter feasibility and objective closer to
        // the constrained optimum (min x² s.t. x ≥ 2 → f* = 4).
        let obj = |x: &[f64]| x[0] * x[0];
        let g = |x: &[f64]| 2.0 - x[0];
        let cons = [ConstraintFn { g: &g, is_inequality: true }];
        let opt = eq_opt(&obj, &cons);
        let (x_lo, f_lo, _) = penalty_optimize_with_mu(&opt, &[0.0], &[1.0, 10.0]).unwrap();
        let (x_hi, f_hi, _) =
            penalty_optimize_with_mu(&opt, &[0.0], &[1.0, 10.0, 100.0, 1000.0]).unwrap();
        let v_lo = constraint_violation(&opt, &x_lo);
        let v_hi = constraint_violation(&opt, &x_hi);
        assert!(v_hi < v_lo, "violation should shrink with μ: {v_lo} -> {v_hi}");
        assert!(
            (f_hi - 4.0).abs() < (f_lo - 4.0).abs(),
            "objective should approach 4: {f_lo} -> {f_hi}"
        );
    }

    #[test]
    fn lagrange_equality_converges() {
        // min x²+y²  s.t.  x + y = 1.  Augmented Lagrangian reaches near-exact
        // feasibility within a few multiplier updates.
        let obj = |x: &[f64]| x[0] * x[0] + x[1] * x[1];
        let g = |x: &[f64]| x[0] + x[1] - 1.0;
        let cons = [ConstraintFn { g: &g, is_inequality: false }];
        let opt = eq_opt(&obj, &cons);
        let (x, fx, _) = lagrange_penalty_optimize(&opt, &[0.2, 0.2]).unwrap();
        assert!((x[0] - 0.5).abs() < 1e-3, "x0 = {}", x[0]);
        assert!((x[1] - 0.5).abs() < 1e-3, "x1 = {}", x[1]);
        assert!(constraint_violation(&opt, &x) < 1e-4, "violation = {}", constraint_violation(&opt, &x));
        assert!((fx - 0.5).abs() < 2e-3, "objective = {fx}");
    }

    #[test]
    fn non_finite_initial_is_error() {
        let obj = |x: &[f64]| x[0] * x[0];
        let g = |x: &[f64]| 2.0 - x[0];
        let cons = [ConstraintFn { g: &g, is_inequality: true }];
        let opt = eq_opt(&obj, &cons);
        assert!(penalty_optimize(&opt, &[f64::NAN]).is_err());
        assert!(penalty_optimize(&opt, &[f64::INFINITY]).is_err());
        assert!(lagrange_penalty_optimize(&opt, &[f64::NAN]).is_err());
        // Empty initial point.
        assert!(penalty_optimize(&opt, &[]).is_err());
        // Empty mu schedule.
        assert!(penalty_optimize_with_mu(&opt, &[1.0], &[]).is_err());
    }

    #[test]
    fn scaled_problem_sanity() {
        // min 10(x²+y²)  s.t.  x + y = 4  →  (2, 2), f = 80.
        let obj = |x: &[f64]| 10.0 * (x[0] * x[0] + x[1] * x[1]);
        let g = |x: &[f64]| x[0] + x[1] - 4.0;
        let cons = [ConstraintFn { g: &g, is_inequality: false }];
        let opt = eq_opt(&obj, &cons);
        let (x, fx, _) = penalty_optimize(&opt, &[0.0, 0.0]).unwrap();
        assert!((x[0] - 2.0).abs() < 2e-2, "x0 = {}", x[0]);
        assert!((x[1] - 2.0).abs() < 2e-2, "x1 = {}", x[1]);
        assert!((fx - 80.0).abs() < 2.0, "objective = {fx}");
    }
}
