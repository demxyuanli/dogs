//! Runge-Kutta 4th-order ODE integrator.
//! Source: classical RK4 (standard textbook algorithm).

/// Advance `y` by one RK4 step of size `dt` for `dy/dt = f(t, y)`.
pub fn rk4_step<F>(f: &F, t: f64, y: &[f64], dt: f64) -> Vec<f64>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let n = y.len();
    let k1 = f(t, y);
    let mut y2 = vec![0.0; n];
    for i in 0..n {
        y2[i] = y[i] + 0.5 * dt * k1[i];
    }
    let k2 = f(t + 0.5 * dt, &y2);
    for i in 0..n {
        y2[i] = y[i] + 0.5 * dt * k2[i];
    }
    let k3 = f(t + 0.5 * dt, &y2);
    for i in 0..n {
        y2[i] = y[i] + dt * k3[i];
    }
    let k4 = f(t + dt, &y2);
    let mut out = vec![0.0; n];
    for i in 0..n {
        out[i] = y[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    out
}

/// Integrate `dy/dt = f(t, y)` from `t0` to `t1` in `steps` equal steps.
///
/// Returns `steps + 1` states; the first is `y0`, the last is the solution at `t1`.
pub fn rk4<F>(f: &F, t0: f64, t1: f64, y0: &[f64], steps: usize) -> Vec<Vec<f64>>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    assert!(steps > 0, "rk4: steps must be > 0");
    let dt = (t1 - t0) / steps as f64;
    let mut out = Vec::with_capacity(steps + 1);
    out.push(y0.to_vec());
    let mut t = t0;
    let mut y = y0.to_vec();
    for _ in 0..steps {
        y = rk4_step(f, t, &y, dt);
        t += dt;
        out.push(y.clone());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exponential_growth() {
        // y' = y, y(0) = 1 -> y(1) = e.
        let traj = rk4(&|_: f64, y: &[f64]| vec![y[0]], 0.0, 1.0, &[1.0], 100);
        let last = traj.last().unwrap();
        assert!((last[0] - std::f64::consts::E).abs() < 1e-9, "y(1) = {}", last[0]);
        assert_eq!(traj.len(), 101);
        assert_eq!(traj[0][0], 1.0);
    }

    #[test]
    fn step_size_consistency() {
        // Logistic y' = y(1 - y), y(0) = 0.5: stable, bounded on [0, 2].
        let f = |_: f64, y: &[f64]| vec![y[0] * (1.0 - y[0])];
        let y0 = [0.5];
        let a = rk4(&f, 0.0, 2.0, &y0, 50);
        let b = rk4(&f, 0.0, 2.0, &y0, 400);
        let diff = (a.last().unwrap()[0] - b.last().unwrap()[0]).abs();
        assert!(diff < 1e-4, "refinement drift {diff}");
    }
}
