//! Levenberg-Marquardt nonlinear least squares.
//! Source: standard algorithm (Numerical Recipes `mrqmin`).

use crate::{MathMatrix, MathVector};

/// Configuration for [`levenberg_marquardt`].
#[derive(Debug, Clone, Copy)]
pub struct LMConfig {
    pub max_iter: usize,
    pub tol: f64,
    pub lambda: f64,
}

impl Default for LMConfig {
    fn default() -> Self {
        Self { max_iter: 200, tol: 1e-10, lambda: 1e-3 }
    }
}

/// Fit parameters `x` minimizing `||f(x)||²` where `f` returns the residual vector.
///
/// The Jacobian is approximated by central differences. Errors when the
/// damped normal system is singular or convergence is not reached.
pub fn levenberg_marquardt<F>(f: &F, x0: &[f64], params: &LMConfig) -> Result<Vec<f64>, String>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut lambda = params.lambda;
    let mut r = f(&x);
    let mut cost = 0.5 * r.iter().map(|v| v * v).sum::<f64>();
    let eps = 1e-7;

    for _ in 0..params.max_iter {
        let m = r.len();
        // Central-difference Jacobian, m×n.
        let mut j = vec![vec![0.0; n]; m];
        for col in 0..n {
            let h = eps * (1.0 + x[col].abs());
            let mut xp = x.clone();
            xp[col] += h;
            let mut xm = x.clone();
            xm[col] -= h;
            let fp = f(&xp);
            let fm = f(&xm);
            for i in 0..m {
                j[i][col] = (fp[i] - fm[i]) / (2.0 * h);
            }
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
        if jtr.iter().fold(0.0_f64, |acc, &v| acc.max(v.abs())) < params.tol {
            return Ok(x);
        }

        // Solve (JᵀJ + λ·diag(JᵀJ)) δ = -Jᵀr.
        let mut system = MathMatrix::new(1, n, 1, n);
        let mut bvec = MathVector::new(1, n);
        for a in 0..n {
            for b in 0..n {
                let mut v = jtj[a][b];
                if a == b {
                    v += lambda * (jtj[a][b] + 1e-12);
                }
                system.set_value(a + 1, b + 1, v);
            }
            bvec.set_value(a + 1, -jtr[a]);
        }
        let delta = system.solve(&bvec).map_err(|_| "levenberg_marquardt: singular system".to_string())?;

        // Trial step: accept if it reduces the cost, otherwise increase damping.
        let mut xt = x.clone();
        for a in 0..n {
            xt[a] += delta.value(a + 1);
        }
        let rt = f(&xt);
        let cost_t = 0.5 * rt.iter().map(|v| v * v).sum::<f64>();
        if cost_t < cost {
            x = xt;
            r = rt;
            cost = cost_t;
            lambda = (lambda / 10.0).max(1e-15);
        } else {
            lambda = (lambda * 10.0).min(1e15);
        }
        let step = (1..=n).map(|i| delta.value(i).abs()).sum::<f64>();
        if step < params.tol {
            return Ok(x);
        }
    }
    Err("levenberg_marquardt: max iterations reached".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_sin_amplitude_and_frequency() {
        // Recover (a, b) from noisy samples of y = a·sin(b·x).
        let n = 60;
        let xs: Vec<f64> = (0..n).map(|i| 3.0 * i as f64 / (n - 1) as f64).collect();
        let ys: Vec<f64> = xs
            .iter()
            .enumerate()
            .map(|(i, &x)| 2.0 * (3.0 * x).sin() + 0.02 * (37.0 * x + i as f64).sin())
            .collect();

        let resid = |p: &[f64]| -> Vec<f64> {
            xs.iter()
                .zip(&ys)
                .map(|(&x, &y)| p[0] * (p[1] * x).sin() - y)
                .collect()
        };

        let cfg = LMConfig::default();
        let p = levenberg_marquardt(&resid, &[1.0, 1.0], &cfg).unwrap();
        assert!((p[0] - 2.0).abs() < 0.1, "a = {}", p[0]);
        assert!((p[1] - 3.0).abs() < 0.1, "b = {}", p[1]);
    }
}
