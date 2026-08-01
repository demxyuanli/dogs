//! Adams–Bashforth / Adams–Moulton multistep ODE integrators.
//! Source: `math_...` / classical linear multistep methods (Hairer, Nørsett &
//! Wanner, "Solving Ordinary Differential Equations I").
//!
//! # Background
//!
//! A `k`-step linear multistep method for `dy/dt = f(t, y)` is a recurrence
//!
//! ```text
//! Σ_{j=0}^{k} α_j · y_{n+j} = h · Σ_{j=0}^{k} β_j · f(t_{n+j}, y_{n+j})
//! ```
//!
//! with `α_k = 1`.  When `β_k = 0` the method is *explicit* (the new value is
//! isolated on the left); otherwise it is *implicit* and must be solved for
//! `y_{n+k}`.  The Adams family chooses the coefficients so that the method is
//! exact for all polynomials of degree `k` (or `k+1`), which maximises the
//! *order* for a given number of stored slopes.  The explicit Adams–Bashforth
//! (AB) methods use past slopes only; the implicit Adams–Moulton (AM) methods
//! additionally involve the slope at the unknown new point and are applied
//! here in predictor–corrector form.
//!
//! # Accuracy and stability
//!
//! The global error of an order-`p` method is `O(hᵖ)`; AB2 is second order,
//! AB3 third.  The implicit trapezoidal rule (AM2) is A-stable: its
//! amplification factor satisfies `|R(z)| ≤ 1` for every `Re(z) ≤ 0`, so it
//! cannot blow up on stiff or oscillatory problems regardless of step size.
//! The explicit AB methods have only a finite region of absolute stability —
//! on the imaginary axis AB2 is only weakly stable, which shows up as a tiny
//! amplitude drift over very long runs.  In practice, for a smooth non-stiff
//! problem and a modest step size, AB3 is a cheap and very accurate choice.
//!
//! # Startup
//!
//! A `k`-step method needs `k` past solution values before the recurrence can
//! begin.  The high-order embedded Dormand–Prince 5(4) step ([`rk45_step`]) is
//! used to generate these startup values; its one-step error `O(h⁵)` is far
//! below the global order of the multistep methods, so the contamination from
//! the startup phase is negligible in the asymptotic error balance.

use crate::ode_rk45::rk45_step;

/// Selector for the fixed-step multistep integrator used by
/// [`ode_system_solve`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OdeMultistepMethod {
    /// 2nd-order explicit Adams–Bashforth, 2-step:
    /// `y_{n+1} = y_n + h·(3·f_n − f_{n−1})/2`.
    AdamsBashforth2,
    /// 3rd-order explicit Adams–Bashforth, 3-step:
    /// `y_{n+1} = y_n + h·(23·f_n − 16·f_{n−1} + 5·f_{n−2})/12`.
    AdamsBashforth3,
    /// 2nd-order implicit Adams–Moulton (trapezoidal) in predictor–corrector
    /// form: AB2 predicts, then the AM2 corrector
    /// `y_{n+1} = y_n + h·(f_{n+1} + f_n)/2` is iterated to a fixed point.
    AdamsMoulton2,
}

/// Shared argument validation for every multistep driver.
///
/// All methods integrate a non-empty system with a finite, non-zero step `h`
/// over at least one step.  Returning an error for these degenerate inputs
/// keeps the drivers total (they never panic on malformed arguments).
fn validate_args(name: &str, t0: f64, y0: &[f64], h: f64, steps: usize) -> Result<(), String> {
    if y0.is_empty() {
        return Err(format!("{name}: empty state vector"));
    }
    if !t0.is_finite() || y0.iter().any(|v| !v.is_finite()) {
        return Err(format!("{name}: non-finite initial data"));
    }
    if !h.is_finite() || h == 0.0 {
        return Err(format!("{name}: step size must be finite and non-zero"));
    }
    if steps == 0 {
        return Err(format!("{name}: steps must be at least 1"));
    }
    Ok(())
}

/// Advance `k` fixed steps of the embedded Dormand–Prince 5(4) pair.
///
/// A multistep method needs the solution on several past nodes before its
/// recurrence can start.  The high-order RK45 step is used only to generate
/// these startup values; its one-step error (O(h⁵)) is far below the global
/// order of the multistep methods, so the startup contamination is negligible.
fn rk_seed(
    f: &dyn Fn(f64, &[f64], &mut [f64]),
    t0: f64,
    y0: &[f64],
    h: f64,
    k: usize,
) -> Vec<(f64, Vec<f64>)> {
    let n = y0.len();
    let mut traj = Vec::with_capacity(k + 1);
    traj.push((t0, y0.to_vec()));
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut y_out = vec![0.0; n];
    let mut err = vec![0.0; n];
    for _ in 0..k {
        rk45_step(f, t, &y, h, &mut y_out, &mut err);
        t += h;
        y.copy_from_slice(&y_out);
        traj.push((t, y.clone()));
    }
    traj
}

/// Integrate `dy/dt = f(t, y)` with the 2-step Adams–Bashforth method.
///
/// The recurrence
///
/// ```text
/// y_{n+1} = y_n + h·(3·f_n − f_{n−1})/2
/// ```
///
/// is the 2nd-order explicit member of the Adams family.  Its local
/// truncation error is `5/12 · h³ · y‴(ξ)`, and the method is exact for
/// linear `y`.  The single required startup step is taken with [`rk45_step`]
/// so that `f_{n−1}` is available from the first regular step.
///
/// # Arguments
/// * `f` — right-hand side `f(t, y) → dy/dt`; the system dimension is taken
///   from `y0`.
/// * `t0` — initial time.
/// * `y0` — initial state vector.
/// * `h` — fixed step size (positive or negative, but never zero).
/// * `steps` — number of steps to take.
///
/// # Returns
/// The trajectory of `steps + 1` nodes `(t0 + k·h, y_k)` for `k = 0..=steps`.
///
/// # Errors
/// Empty `y0`, a zero or non-finite step, or `steps == 0` all produce an
/// `Err` (see [`validate_args`]).
pub fn adams_bashforth_2(
    f: &dyn Fn(f64, &[f64], &mut [f64]),
    t0: f64,
    y0: &[f64],
    h: f64,
    steps: usize,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    validate_args("adams_bashforth_2", t0, y0, h, steps)?;
    let n = y0.len();
    let mut traj = rk_seed(f, t0, y0, h, 1);
    let mut f_prev = vec![0.0; n];
    let mut f_curr = vec![0.0; n];
    f(t0, y0, &mut f_prev);
    f(t0 + h, &traj[1].1, &mut f_curr);
    for _ in 1..steps {
        let (t_last, y_last) = (traj.last().unwrap().0, traj.last().unwrap().1.clone());
        let t_next = t_last + h;
        let mut y_next = vec![0.0; n];
        for i in 0..n {
            y_next[i] = y_last[i] + h * (1.5 * f_curr[i] - 0.5 * f_prev[i]);
        }
        f_prev.copy_from_slice(&f_curr);
        f(t_next, &y_next, &mut f_curr);
        traj.push((t_next, y_next));
    }
    Ok(traj)
}

/// Integrate `dy/dt = f(t, y)` with the 3-step Adams–Bashforth method.
///
/// The recurrence
///
/// ```text
/// y_{n+1} = y_n + h·(23·f_n − 16·f_{n−1} + 5·f_{n−2})/12
/// ```
///
/// is the 3rd-order explicit member of the Adams family.  Its local
/// truncation error is `3/8 · h⁴ · y⁗(ξ)`, so at a fixed step size it is
/// roughly a factor of `h` more accurate than AB2.  Two startup steps are
/// taken with [`rk45_step`].
///
/// # Arguments
/// * `f` — right-hand side `f(t, y) → dy/dt`.
/// * `t0` — initial time.
/// * `y0` — initial state vector.
/// * `h` — fixed step size (never zero).
/// * `steps` — number of steps to take.
///
/// # Returns
/// The trajectory of `steps + 1` nodes `(t0 + k·h, y_k)` for `k = 0..=steps`.
///
/// # Errors
/// Same validation as [`adams_bashforth_2`]; additionally a single-step run is
/// handled by seeding one RK45 step (so it never over-seeds and a `steps == 1`
/// call still returns exactly two nodes).
pub fn adams_bashforth_3(
    f: &dyn Fn(f64, &[f64], &mut [f64]),
    t0: f64,
    y0: &[f64],
    h: f64,
    steps: usize,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    validate_args("adams_bashforth_3", t0, y0, h, steps)?;
    if steps == 1 {
        return Ok(rk_seed(f, t0, y0, h, 1));
    }
    let n = y0.len();
    let mut traj = rk_seed(f, t0, y0, h, 2);
    let mut fs = vec![vec![0.0; n]; 3];
    f(t0, y0, &mut fs[0]);
    f(t0 + h, &traj[1].1, &mut fs[1]);
    f(t0 + 2.0 * h, &traj[2].1, &mut fs[2]);
    for _ in 2..steps {
        let (t_last, y_last) = (traj.last().unwrap().0, traj.last().unwrap().1.clone());
        let t_next = t_last + h;
        let mut y_next = vec![0.0; n];
        for i in 0..n {
            y_next[i] = y_last[i] + h * (23.0 * fs[2][i] - 16.0 * fs[1][i] + 5.0 * fs[0][i]) / 12.0;
        }
        fs[0] = fs[1].clone();
        fs[1] = fs[2].clone();
        f(t_next, &y_next, &mut fs[2]);
        traj.push((t_next, y_next));
    }
    Ok(traj)
}

/// Integrate `dy/dt = f(t, y)` with the Adams–Moulton predictor–corrector.
///
/// The AB2 formula is used as the predictor, then the 2nd-order Adams–Moulton
/// corrector (the trapezoidal rule)
///
/// ```text
/// y_{n+1} = y_n + h·(f(t_{n+1}, y_{n+1}) + f_n)/2
/// ```
///
/// is iterated up to three times, warming the RHS with the newly corrected
/// value so the final iterate is a good fixed-point approximation of the
/// implicit method.
///
/// The trapezoidal rule is A-stable: for `y' = λy` with `Re(λ) ≤ 0` its
/// amplification factor `R(z) = (1 + z/2)/(1 − z/2)` satisfies `|R(z)| ≤ 1`
/// for every step size.  For a harmonic oscillator the amplitude is therefore
/// preserved exactly and only the phase drifts — a qualitative improvement
/// over the explicit AB2, which slowly gains or loses amplitude on the
/// imaginary axis.  One startup step comes from [`rk45_step`].
///
/// # Arguments
/// * `f` — right-hand side `f(t, y) → dy/dt`.
/// * `t0` — initial time.
/// * `y0` — initial state vector.
/// * `h` — fixed step size (never zero).
/// * `steps` — number of steps to take.
///
/// # Returns
/// The trajectory of `steps + 1` nodes `(t0 + k·h, y_k)` for `k = 0..=steps`.
///
/// # Errors
/// Same validation as [`adams_bashforth_2`].
pub fn adams_moulton_2(
    f: &dyn Fn(f64, &[f64], &mut [f64]),
    t0: f64,
    y0: &[f64],
    h: f64,
    steps: usize,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    validate_args("adams_moulton_2", t0, y0, h, steps)?;
    let n = y0.len();
    let mut traj = rk_seed(f, t0, y0, h, 1);
    let mut f_prev = vec![0.0; n];
    let mut f_curr = vec![0.0; n];
    f(t0, y0, &mut f_prev);
    f(t0 + h, &traj[1].1, &mut f_curr);
    for _ in 1..steps {
        let (t_last, y_last) = (traj.last().unwrap().0, traj.last().unwrap().1.clone());
        let t_next = t_last + h;
        let mut y_pred = vec![0.0; n];
        for i in 0..n {
            y_pred[i] = y_last[i] + h * (1.5 * f_curr[i] - 0.5 * f_prev[i]);
        }
        // Predictor–corrector: iterate the implicit trapezoidal corrector.
        let mut y_corr = y_pred;
        let mut f_next = vec![0.0; n];
        for _ in 0..3 {
            f(t_next, &y_corr, &mut f_next);
            for i in 0..n {
                y_corr[i] = y_last[i] + 0.5 * h * (f_next[i] + f_curr[i]);
            }
        }
        f_prev.copy_from_slice(&f_curr);
        f(t_next, &y_corr, &mut f_curr);
        traj.push((t_next, y_corr));
    }
    Ok(traj)
}

/// Integrate `dy/dt = f(t, y)` from `t0` to `t1` with a fixed step `h`.
///
/// The number of steps is `round(|t1 − t0| / |h|)`; the direction of `h` must
/// agree with the integration interval (both advance toward `t1`).  The
/// trajectory is produced by the selected [`OdeMultistepMethod`] and contains
/// `steps + 1` nodes, the last of which lies at `t0 + steps·h ≈ t1`.
///
/// This is the convenient entry point for a complete integration over an
/// interval: pick a method from [`OdeMultistepMethod`], a step size `h`, and
/// the driver computes the node count for you.  For finer control over the
/// number of steps, call the individual `adams_*` functions directly.
///
/// # Arguments
/// * `f` — right-hand side `f(t, y) → dy/dt`.
/// * `t0` — initial time.
/// * `y0` — initial state vector.
/// * `t1` — final time.
/// * `h` — fixed step size; its sign must point from `t0` toward `t1`.
/// * `method` — the multistep scheme to use.
///
/// # Errors
/// Empty `y0`, a zero/non-finite step, a direction mismatch between `h` and
/// the interval, or a span shorter than one step all produce an `Err`.
pub fn ode_system_solve(
    f: &dyn Fn(f64, &[f64], &mut [f64]),
    t0: f64,
    y0: &[f64],
    t1: f64,
    h: f64,
    method: OdeMultistepMethod,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    if y0.is_empty() {
        return Err("ode_system_solve: empty state vector".to_string());
    }
    if !h.is_finite() || h == 0.0 {
        return Err("ode_system_solve: step size must be finite and non-zero".to_string());
    }
    let span = t1 - t0;
    if span == 0.0 {
        return Ok(vec![(t0, y0.to_vec())]);
    }
    if h.signum() != span.signum() {
        return Err("ode_system_solve: step direction must match the integration interval".to_string());
    }
    let steps = (span / h).abs().round();
    if steps < 1.0 {
        return Err("ode_system_solve: interval shorter than one step".to_string());
    }
    let steps = steps as usize;
    match method {
        OdeMultistepMethod::AdamsBashforth2 => adams_bashforth_2(f, t0, y0, h, steps),
        OdeMultistepMethod::AdamsBashforth3 => adams_bashforth_3(f, t0, y0, h, steps),
        OdeMultistepMethod::AdamsMoulton2 => adams_moulton_2(f, t0, y0, h, steps),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Harmonic oscillator `y'' = −y` as a first-order system with
    /// `y(0) = 0, y'(0) = 1`, whose exact solution is `(sin t, cos t)`.
    fn harmonic() -> impl Fn(f64, &[f64], &mut [f64]) {
        |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[1];
            dy[1] = -y[0];
        }
    }

    #[test]
    fn ab2_harmonic_tracks_sin() {
        let f = harmonic();
        let traj = adams_bashforth_2(&f, 0.0, &[0.0, 1.0], 0.01, 1000).unwrap();
        let (t, y) = traj.last().unwrap();
        let exact = t.sin();
        assert!(
            (y[0] - exact).abs() < 5e-3,
            "AB2 y({t}) = {} vs sin = {exact}",
            y[0]
        );
    }

    #[test]
    fn ab3_more_accurate_than_ab2() {
        let f = harmonic();
        let ab2_err = {
            let traj = adams_bashforth_2(&f, 0.0, &[0.0, 1.0], 0.01, 1000).unwrap();
            let (t, y) = traj.last().unwrap();
            (y[0] - t.sin()).abs()
        };
        let ab3_err = {
            let traj = adams_bashforth_3(&f, 0.0, &[0.0, 1.0], 0.01, 1000).unwrap();
            let (t, y) = traj.last().unwrap();
            (y[0] - t.sin()).abs()
        };
        assert!(
            ab3_err < ab2_err,
            "AB3 error {ab3_err} should beat AB2 error {ab2_err}"
        );
    }

    #[test]
    fn am2_more_accurate_than_ab2() {
        let f = harmonic();
        let ab2_err = {
            let traj = adams_bashforth_2(&f, 0.0, &[0.0, 1.0], 0.01, 1000).unwrap();
            let (t, y) = traj.last().unwrap();
            (y[0] - t.sin()).abs()
        };
        let am2_err = {
            let traj = adams_moulton_2(&f, 0.0, &[0.0, 1.0], 0.01, 1000).unwrap();
            let (t, y) = traj.last().unwrap();
            (y[0] - t.sin()).abs()
        };
        assert!(
            am2_err < ab2_err,
            "AM2 error {am2_err} should beat AB2 error {ab2_err}"
        );
    }

    #[test]
    fn ab2_exponential_decay() {
        // y' = −y, y(0) = 1 → y = e^(−t). The AB2 principal root matches the
        // exact advancing factor to high accuracy on real eigenvalues.
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[0];
        };
        let traj = adams_bashforth_2(&f, 0.0, &[1.0], 0.01, 500).unwrap();
        let (t, y) = traj.last().unwrap();
        let exact = (-t).exp();
        assert!((y[0] - exact).abs() < 1e-4, "y({t}) = {} vs e^-t = {exact}", y[0]);
    }

    #[test]
    fn two_by_two_linear_matches_analytic() {
        // y1' = −y2, y2' = y1, y(0) = (1, 0) → (cos t, sin t).
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[1];
            dy[1] = y[0];
        };
        let traj = adams_bashforth_2(&f, 0.0, &[1.0, 0.0], 0.01, 314).unwrap();
        let (t, y) = traj.last().unwrap();
        assert!((y[0] - t.cos()).abs() < 1e-3, "y1 = {} vs cos = {}", y[0], t.cos());
        assert!((y[1] - t.sin()).abs() < 1e-3, "y2 = {} vs sin = {}", y[1], t.sin());
    }

    #[test]
    fn ode_system_solve_reaches_endpoint() {
        let f = harmonic();
        let traj = ode_system_solve(
            &f,
            0.0,
            &[0.0, 1.0],
            10.0,
            0.01,
            OdeMultistepMethod::AdamsBashforth2,
        )
        .unwrap();
        let (t, _) = traj.last().unwrap();
        assert!((t - 10.0).abs() < 1e-9, "endpoint t = {t}");
        assert_eq!(traj.len(), 1001);
        // All three methods agree on the node count.
        for m in [
            OdeMultistepMethod::AdamsBashforth3,
            OdeMultistepMethod::AdamsMoulton2,
        ] {
            let tr = ode_system_solve(&f, 0.0, &[0.0, 1.0], 5.0, 0.01, m).unwrap();
            assert_eq!(tr.len(), 501);
        }
    }

    #[test]
    fn empty_state_is_error() {
        let f = harmonic();
        assert!(adams_bashforth_2(&f, 0.0, &[], 0.1, 10).is_err());
        assert!(adams_bashforth_3(&f, 0.0, &[], 0.1, 10).is_err());
        assert!(adams_moulton_2(&f, 0.0, &[], 0.1, 10).is_err());
        assert!(
            ode_system_solve(&f, 0.0, &[], 1.0, 0.1, OdeMultistepMethod::AdamsBashforth2)
                .is_err()
        );
    }

    #[test]
    fn step_size_validation() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[0];
        };
        // Zero step.
        assert!(adams_bashforth_2(&f, 0.0, &[1.0], 0.0, 10).is_err());
        assert!(adams_bashforth_3(&f, 0.0, &[1.0], 0.0, 10).is_err());
        assert!(adams_moulton_2(&f, 0.0, &[1.0], 0.0, 10).is_err());
        // Zero steps.
        assert!(adams_bashforth_2(&f, 0.0, &[1.0], 0.1, 0).is_err());
        // Direction mismatch in ode_system_solve.
        assert!(
            ode_system_solve(&f, 0.0, &[1.0], 1.0, -0.1, OdeMultistepMethod::AdamsBashforth2)
                .is_err()
        );
    }

    #[test]
    fn single_step_returns_two_nodes() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[0];
        };
        let tr = adams_bashforth_2(&f, 0.0, &[1.0], 0.1, 1).unwrap();
        assert_eq!(tr.len(), 2);
        assert!((tr[0].0 - 0.0).abs() < 1e-12 && (tr[1].0 - 0.1).abs() < 1e-12);
        // A one-step AB3 run also yields two nodes.
        let tr3 = adams_bashforth_3(&f, 0.0, &[1.0], 0.1, 1).unwrap();
        assert_eq!(tr3.len(), 2);
    }
}
