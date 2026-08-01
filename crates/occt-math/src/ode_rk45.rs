//! Adaptive Runge-Kutta-Fehlberg 5(4) ODE integrator (Dormand-Prince).
//! Source: `math_RK45` / classical RK45 (DOPRI5), Shampine & Gordon.

/// Coefficients of the Dormand-Prince 5(4) embedded pair.
///
/// Butcher tableau (7 stages):
/// ```text
///  0   |
///  1/5 | 1/5
///  3/10| 3/40     9/40
///  4/5 | 44/45    -56/15   32/9
///  8/9 | 19372/6561  -25360/2187  64448/6561  -212/729
///  1   | 9017/3168 -355/33 46732/5247 49/176 -5103/18656
///  1   | 35/384    0       500/1113   125/192 -2187/6784 11/84
///  ---------------------------------------------------------
///  b   | 35/384 0 500/1113 125/192 -2187/6784 11/84 0
///  b*  | 5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40
/// ```
/// The error estimate is the difference of the two embedded solutions.
mod dopri5 {
    pub const C2: f64 = 1.0 / 5.0;
    pub const C3: f64 = 3.0 / 10.0;
    pub const C4: f64 = 4.0 / 5.0;
    pub const C5: f64 = 8.0 / 9.0;

    pub const A21: f64 = 1.0 / 5.0;
    pub const A31: f64 = 3.0 / 40.0;
    pub const A32: f64 = 9.0 / 40.0;
    pub const A41: f64 = 44.0 / 45.0;
    pub const A42: f64 = -56.0 / 15.0;
    pub const A43: f64 = 32.0 / 9.0;
    pub const A51: f64 = 19372.0 / 6561.0;
    pub const A52: f64 = -25360.0 / 2187.0;
    pub const A53: f64 = 64448.0 / 6561.0;
    pub const A54: f64 = -212.0 / 729.0;
    pub const A61: f64 = 9017.0 / 3168.0;
    pub const A62: f64 = -355.0 / 33.0;
    pub const A63: f64 = 46732.0 / 5247.0;
    pub const A64: f64 = 49.0 / 176.0;
    pub const A65: f64 = -5103.0 / 18656.0;
    pub const A71: f64 = 35.0 / 384.0;
    pub const A73: f64 = 500.0 / 1113.0;
    pub const A74: f64 = 125.0 / 192.0;
    pub const A75: f64 = -2187.0 / 6784.0;
    pub const A76: f64 = 11.0 / 84.0;

    /// 5th-order weights.
    pub const B1: f64 = 35.0 / 384.0;
    pub const B3: f64 = 500.0 / 1113.0;
    pub const B4: f64 = 125.0 / 192.0;
    pub const B5: f64 = -2187.0 / 6784.0;
    pub const B6: f64 = 11.0 / 84.0;

    /// Error estimate weights e = b - b* (4th-order embedded solution subtracted).
    pub const E1: f64 = 71.0 / 57600.0;
    pub const E3: f64 = -71.0 / 16695.0;
    pub const E4: f64 = 71.0 / 1920.0;
    pub const E5: f64 = -17253.0 / 339200.0;
    pub const E6: f64 = 22.0 / 525.0;
    pub const E7: f64 = -1.0 / 40.0;
}

/// Fixed-step adaptive driver state: RHS `f(t, y, dy)` plus tolerance controls.
pub struct OdeRk45<'a> {
    pub f: &'a dyn Fn(f64, &[f64], &mut [f64]),
    pub atol: f64,
    pub rtol: f64,
    pub max_steps: usize,
}

impl<'a> OdeRk45<'a> {
    pub fn new(
        f: &'a dyn Fn(f64, &[f64], &mut [f64]),
        atol: f64,
        rtol: f64,
        max_steps: usize,
    ) -> Self {
        Self { f, atol, rtol, max_steps }
    }

    /// Integrate and return the full accepted trajectory.
    pub fn solve(&self, t0: f64, y0: &[f64], t1: f64) -> Result<Vec<(f64, Vec<f64>)>, String> {
        rk45_adaptive(self.f, t0, y0, t1, self.atol, self.rtol, self.max_steps)
    }

    /// Integrate and return only the final state.
    pub fn final_state(&self, t0: f64, y0: &[f64], t1: f64) -> Result<(f64, Vec<f64>), String> {
        rk45_solve(self.f, t0, y0, t1, self.atol, self.rtol)
    }
}

/// Stateful adaptive RK45 driver, mirroring OCCT's `math_RK45` class.
///
/// The solver keeps the current time, state and step size internally so a
/// caller can drive the integration step-by-step (e.g. to detect events or
/// sample dense output between accepted nodes), in addition to the one-shot
/// free functions above.
pub struct Rk45Integrator<'a> {
    f: &'a dyn Fn(f64, &[f64], &mut [f64]),
    atol: f64,
    rtol: f64,
    t: f64,
    y: Vec<f64>,
    dt: f64,
    y_out: Vec<f64>,
    err: Vec<f64>,
}

impl<'a> Rk45Integrator<'a> {
    pub fn new(f: &'a dyn Fn(f64, &[f64], &mut [f64]), atol: f64, rtol: f64) -> Self {
        Self { f, atol, rtol, t: 0.0, y: Vec::new(), dt: 0.1, y_out: Vec::new(), err: Vec::new() }
    }

    /// Reset to `y(t0) = y0` with an initial step of 0.1.
    pub fn init(&mut self, t0: f64, y0: &[f64]) {
        self.t = t0;
        self.y = y0.to_vec();
        self.y_out = vec![0.0; y0.len()];
        self.err = vec![0.0; y0.len()];
        self.dt = 0.1;
    }

    /// Current integration time.
    pub fn t(&self) -> f64 { self.t }

    /// Current state `y(t)`.
    pub fn y(&self) -> &[f64] { &self.y }

    /// Step size used by the last accepted step.
    pub fn last_step(&self) -> f64 { self.dt }

    /// Advance one accepted step and return the new time.
    ///
    /// The candidate step is retried at half size until the weighted error
    /// estimate passes; if no acceptable step is found within 10 retries the
    /// integration is aborted with an error.
    pub fn step(&mut self) -> Result<f64, String> {
        if self.y.is_empty() {
            return Err("Rk45Integrator::step: init() has not been called".to_string());
        }
        let mut dt = self.dt;
        for _ in 0..10 {
            rk45_step(self.f, self.t, &self.y, dt, &mut self.y_out, &mut self.err);
            let en = error_norm(&self.y, &self.y_out, &self.err, self.atol, self.rtol);
            if en <= 1.0 {
                self.t += dt;
                self.y.copy_from_slice(&self.y_out);
                self.dt = dt;
                let factor = 0.9 * (1.0 / en.max(1e-300)).powf(0.2);
                self.dt *= factor.min(5.0).max(0.2);
                return Ok(self.t);
            }
            dt *= 0.5;
        }
        Err("Rk45Integrator::step: no acceptable step found".to_string())
    }

    /// Integrate forward to `t1`, leaving the internal state at `t1`.
    pub fn integrate(&mut self, t1: f64) -> Result<f64, String> {
        while (self.t - t1).abs() > 1e-12 {
            let rem = t1 - self.t;
            if self.dt.abs() > rem.abs() {
                self.dt = rem;
            }
            if self.dt == 0.0 {
                break;
            }
            self.step()?;
        }
        Ok(self.t)
    }
}

/// One Dormand-Prince 5(4) step: fills `y_out` with the 5th-order solution and
/// `err` with the local error estimate (difference of the embedded pair).
pub fn rk45_step(
    f: &dyn Fn(f64, &[f64], &mut [f64]),
    t: f64,
    y: &[f64],
    dt: f64,
    y_out: &mut [f64],
    err: &mut [f64],
) {
    use dopri5::*;
    let n = y.len();
    let mut k = vec![vec![0.0; n]; 7];
    let mut yt = vec![0.0; n];

    // Stage 1
    f(t, y, &mut k[0]);
    // Stage 2
    for i in 0..n { yt[i] = y[i] + dt * A21 * k[0][i]; }
    f(t + dt * C2, &yt, &mut k[1]);
    // Stage 3
    for i in 0..n { yt[i] = y[i] + dt * (A31 * k[0][i] + A32 * k[1][i]); }
    f(t + dt * C3, &yt, &mut k[2]);
    // Stage 4
    for i in 0..n { yt[i] = y[i] + dt * (A41 * k[0][i] + A42 * k[1][i] + A43 * k[2][i]); }
    f(t + dt * C4, &yt, &mut k[3]);
    // Stage 5
    for i in 0..n { yt[i] = y[i] + dt * (A51 * k[0][i] + A52 * k[1][i] + A53 * k[2][i] + A54 * k[3][i]); }
    f(t + dt * C5, &yt, &mut k[4]);
    // Stage 6
    for i in 0..n { yt[i] = y[i] + dt * (A61 * k[0][i] + A62 * k[1][i] + A63 * k[2][i] + A64 * k[3][i] + A65 * k[4][i]); }
    f(t + dt, &yt, &mut k[5]);
    // Stage 7 (used for the error estimate)
    for i in 0..n { yt[i] = y[i] + dt * (A71 * k[0][i] + A73 * k[2][i] + A74 * k[3][i] + A75 * k[4][i] + A76 * k[5][i]); }
    f(t + dt, &yt, &mut k[6]);

    // 5th-order solution.
    for i in 0..n {
        y_out[i] = y[i] + dt * (B1 * k[0][i] + B3 * k[2][i] + B4 * k[3][i] + B5 * k[4][i] + B6 * k[5][i]);
    }
    // Error estimate = y5 - y4 = dt * Σ (b - b*)·k.
    for i in 0..n {
        err[i] = dt * (E1 * k[0][i] + E3 * k[2][i] + E4 * k[3][i] + E5 * k[4][i] + E6 * k[5][i] + E7 * k[6][i]);
    }
}

/// Weighted RMS error norm used by the adaptive controller.
fn error_norm(y: &[f64], y_out: &[f64], err: &[f64], atol: f64, rtol: f64) -> f64 {
    let n = y.len();
    let mut s = 0.0;
    for i in 0..n {
        let scale = atol + rtol * y[i].abs().max(y_out[i].abs());
        let e = err[i] / scale;
        s += e * e;
    }
    (s / n as f64).sqrt()
}

/// Integrate `dy/dt = f(t, y)` from `t0` to `t1` with adaptive step size.
///
/// The step is accepted when the weighted RMS error is <= 1; otherwise it is
/// retried with a smaller step. The next step is
/// `dt = dt * min(5, max(0.2, 0.9 * (tol/err)^(1/5)))`.
/// Returns the accepted `(t, y)` trajectory, or an `Err` when the step count is
/// exhausted or the step underflows.
pub fn rk45_adaptive(
    f: &dyn Fn(f64, &[f64], &mut [f64]),
    t0: f64,
    y0: &[f64],
    t1: f64,
    atol: f64,
    rtol: f64,
    max_steps: usize,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    let n = y0.len();
    if n == 0 {
        return Err("rk45_adaptive: zero-dimensional system".to_string());
    }
    let span = t1 - t0;
    if span == 0.0 {
        return Ok(vec![(t0, y0.to_vec())]);
    }
    let dir = span.signum();
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut traj = vec![(t0, y0.to_vec())];
    let mut dt = dir * span.abs() * 0.1;
    let mut y_out = vec![0.0; n];
    let mut err = vec![0.0; n];

    for _ in 0..max_steps {
        let remaining = (t1 - t) * dir;
        if remaining <= 1e-15 * span.abs().max(1.0) {
            break;
        }
        if dt.abs() > remaining {
            dt = dir * remaining;
        }
        if dt.abs() < 1e-300 {
            return Err("rk45_adaptive: step size underflow".to_string());
        }
        rk45_step(f, t, &y, dt, &mut y_out, &mut err);
        let en = error_norm(&y, &y_out, &err, atol, rtol);
        if en <= 1.0 {
            t += dt;
            y.copy_from_slice(&y_out);
            traj.push((t, y.clone()));
            if (t - t1).abs() <= 1e-14 * span.abs().max(1.0) {
                break;
            }
        }
        let factor = 0.9 * (1.0 / en.max(1e-300)).powf(0.2);
        let factor = factor.min(5.0).max(0.2);
        dt *= factor;
    }
    if (t - t1).abs() > 1e-9 * span.abs().max(1.0) {
        return Err(format!(
            "rk45_adaptive: max_steps ({max_steps}) exceeded, reached t = {t}"
        ));
    }
    Ok(traj)
}

/// Integrate and return only the final state `(t1, y(t1))`.
pub fn rk45_solve(
    f: &dyn Fn(f64, &[f64], &mut [f64]),
    t0: f64,
    y0: &[f64],
    t1: f64,
    atol: f64,
    rtol: f64,
) -> Result<(f64, Vec<f64>), String> {
    let traj = rk45_adaptive(f, t0, y0, t1, atol, rtol, 100_000)?;
    let (t, y) = traj.last().expect("trajectory non-empty");
    Ok((*t, y.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn harmonic_oscillator_tracks_sin() {
        // y0' = y1, y1' = -y0  ->  y0(0)=0, y1(0)=1  =>  y0 = sin(t).
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[1];
            dy[1] = -y[0];
        };
        let traj = rk45_adaptive(&f, 0.0, &[0.0, 1.0], 10.0, 1e-6, 1e-6, 100_000).unwrap();
        let (t, y) = traj.last().unwrap();
        let exact = t.sin();
        assert!(
            (y[0] - exact).abs() < 1e-4,
            "y({t}) = {} vs sin = {exact}",
            y[0]
        );
    }

    #[test]
    fn exponential_decay_tracks_e() {
        // y' = -y, y(0) = 1  ->  y = e^-t.
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[0];
        };
        let (t, y) = rk45_solve(&f, 0.0, &[1.0], 1.0, 1e-8, 1e-8).unwrap();
        let exact = (-t).exp();
        assert!((y[0] - exact).abs() < 1e-6, "y({t}) = {} vs e^-t = {exact}", y[0]);
    }

    #[test]
    fn error_estimate_shrinks_with_tol() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[1];
            dy[1] = -y[0];
        };
        let err_lo = {
            let (_, y) = rk45_solve(&f, 0.0, &[0.0, 1.0], 10.0, 1e-6, 1e-6).unwrap();
            (y[0] - 10.0_f64.sin()).abs()
        };
        let err_hi = {
            let (_, y) = rk45_solve(&f, 0.0, &[0.0, 1.0], 10.0, 1e-8, 1e-8).unwrap();
            (y[0] - 10.0_f64.sin()).abs()
        };
        assert!(err_hi < err_lo, "tight tol should reduce error: {err_lo} -> {err_hi}");
    }

    #[test]
    fn logistic_stays_bounded() {
        // y' = y(1 - y), y(0) = 0.5  ->  lim y = 1 (no blow-up).
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[0] * (1.0 - y[0]);
        };
        let traj = rk45_adaptive(&f, 0.0, &[0.5], 10.0, 1e-8, 1e-8, 100_000).unwrap();
        for (_, y) in &traj {
            assert!((0.0..=1.0).contains(&y[0]), "logistic escaped bound: {}", y[0]);
        }
        assert!((traj.last().unwrap().1[0] - 1.0).abs() < 1e-3);
    }

    #[test]
    fn rk45_solve_matches_trajectory_endpoint() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[1];
            dy[1] = -y[0];
        };
        let traj = rk45_adaptive(&f, 0.0, &[1.0, 0.0], 3.0, 1e-8, 1e-8, 100_000).unwrap();
        let last = traj.last().unwrap().1.clone();
        let (_, sol) = rk45_solve(&f, 0.0, &[1.0, 0.0], 3.0, 1e-8, 1e-8).unwrap();
        for i in 0..2 {
            assert!((sol[i] - last[i]).abs() < 1e-12, "endpoint mismatch at {i}");
        }
    }

    #[test]
    fn zero_dimension_is_error() {
        let f = |_t: f64, _y: &[f64], _dy: &mut [f64]| {};
        assert!(rk45_adaptive(&f, 0.0, &[], 1.0, 1e-6, 1e-6, 10).is_err());
    }

    #[test]
    fn decay_trajectory_is_monotone() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[0];
        };
        let traj = rk45_adaptive(&f, 0.0, &[1.0], 2.0, 1e-8, 1e-8, 100_000).unwrap();
        for w in traj.windows(2) {
            assert!(w[1].1[0] < w[0].1[0], "decay must be strictly decreasing");
        }
    }

    #[test]
    fn first_step_matches_fixed_step_reference() {
        // One RK45 step of dt = 0.1 on y' = -y must match e^-0.1 closely.
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[0];
        };
        let mut out = [0.0];
        let mut err = [0.0];
        rk45_step(&f, 0.0, &[1.0], 0.1, &mut out, &mut err);
        let exact = (-0.1_f64).exp();
        assert!(
            (out[0] - exact).abs() < 1e-6,
            "RK45 step {} vs e^-0.1 = {exact}",
            out[0]
        );
        assert!(err[0].abs() < 1e-5, "error estimate should be small");
    }

    #[test]
    fn max_steps_exceeded_is_error() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = y[1];
            dy[1] = -y[0];
        };
        let res = rk45_adaptive(&f, 0.0, &[0.0, 1.0], 10.0, 1e-6, 1e-6, 5);
        assert!(res.is_err(), "5 steps cannot cover [0,10]");
    }

    #[test]
    fn ode_rk45_struct_driver() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[0];
        };
        let solver = OdeRk45::new(&f, 1e-8, 1e-8, 100_000);
        let (t, y) = solver.final_state(0.0, &[2.0], 2.0).unwrap();
        let exact = 2.0 * (-t).exp();
        assert!((y[0] - exact).abs() < 1e-6);
    }

    #[test]
    fn rk45_integrator_steps_to_endpoint() {
        let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
            dy[0] = -y[0];
        };
        let mut solver = Rk45Integrator::new(&f, 1e-8, 1e-8);
        solver.init(0.0, &[1.0]);
        // Drive a few individual steps, then integrate the rest.
        assert!((solver.step().unwrap() - 0.0).abs() > 0.0);
        assert!(solver.t() > 0.0);
        let t_end = solver.integrate(2.0).unwrap();
        assert!((t_end - 2.0).abs() < 1e-9);
        let exact = (-2.0_f64).exp();
        assert!((solver.y()[0] - exact).abs() < 1e-6, "y = {}", solver.y()[0]);
        assert!(solver.last_step() > 0.0);
    }
}
