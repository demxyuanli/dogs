//! Newton's method for multi-dimensional minimization and root-finding.
//! Source: `math_NewtonMinimum.cxx`, `math_NewtonFunctionSetRoot.cxx`
use crate::{MathVector, MathMatrix, MathStatus};

/// Newton's method for solving F(X) = 0 where F: R^n → R^n.
/// Uses finite-difference Jacobian if no analytical derivative provided.
pub struct NewtonSolver {
    pub max_iter: usize,
    pub tolerance: f64,
}

impl NewtonSolver {
    pub fn new() -> Self { Self { max_iter: 100, tolerance: 1e-10 } }

    /// Solve F(x) = 0 starting from initial guess x0.
    /// F: function returning residual vector.
    /// J: optional Jacobian matrix function J(x).
    pub fn solve<F>(&self, f: F, x0: &MathVector, jacobian: Option<&dyn Fn(&MathVector) -> MathMatrix>) -> Result<MathVector, MathStatus>
    where F: Fn(&MathVector) -> MathVector
    {
        let n = x0.len();
        let mut x = x0.clone();
        let mut fx = f(&x);
        let mut iter = 0usize;

        while iter < self.max_iter {
            // Check convergence
            let err = fx.norm();
            if err < self.tolerance { return Ok(x); }

            // Build Jacobian
            let j = if let Some(ref jfn) = jacobian {
                jfn(&x)
            } else {
                finite_difference_jacobian(&f, &x, &fx)
            };

            // Solve J * dx = -fx
            let neg_fx = fx.opposite();
            let dx = j.solve(&neg_fx).map_err(|_| MathStatus::FunctionError)?;

            // Line search with backtracking
            let mut alpha = 1.0;
            let mut new_x = MathVector::new(1, n);
            for _ in 0..20 {
                for i in 1..=n { new_x.set_value(i, x.value(i) + alpha * dx.value(i)); }
                let new_f = f(&new_x);
                if new_f.norm() < fx.norm() { break; }
                alpha *= 0.5;
            }

            x = new_x;
            fx = f(&x);
            iter += 1;
        }
        if fx.norm() < self.tolerance * 100.0 { Ok(x) }
        else { Err(MathStatus::TooManyIterations) }
    }
}

/// Compute Jacobian via finite differences.
fn finite_difference_jacobian<F: Fn(&MathVector) -> MathVector>(f: &F, x: &MathVector, fx: &MathVector) -> MathMatrix {
    let n = x.len();
    let mut j = MathMatrix::new(1, n, 1, n);
    let eps = 1e-8;
    for col in 1..=n {
        let mut xp = x.clone();
        let h = eps * (1.0 + x.value(col).abs());
        xp.set_value(col, x.value(col) + h);
        let fp = f(&xp);
        for row in 1..=n {
            j.set_value(row, col, (fp.value(row) - fx.value(row)) / h);
        }
    }
    j
}

/// Newton's method for unconstrained minimization of f: R^n → R.
/// Uses gradient and approximate Hessian via BFGS updates.
pub struct NewtonMinimum {
    pub max_iter: usize,
    pub tolerance: f64,
}

impl NewtonMinimum {
    pub fn new() -> Self { Self { max_iter: 200, tolerance: 1e-8 } }

    /// Minimize f(x) starting from x0. grad: gradient function.
    pub fn minimize<F, G>(&self, f: F, grad: G, x0: &MathVector) -> Result<MathVector, MathStatus>
    where F: Fn(&MathVector) -> f64,
          G: Fn(&MathVector) -> MathVector
    {
        let n = x0.len();
        let mut x = x0.clone();
        let mut h = MathMatrix::new(1, n, 1, n); // approximate inverse Hessian (identity)
        for i in 1..=n { h.set_value(i, i, 1.0); }

        let mut g = grad(&x);
        for _iter in 0..self.max_iter {
            if g.norm() < self.tolerance { return Ok(x); }

            // Search direction: d = -H * g (line search direction)
            let mut d = MathVector::new(1, n);
            for i in 1..=n {
                let mut s = 0.0;
                for j in 1..=n { s += h.value(i, j) * g.value(j); }
                d.set_value(i, -s);
            }

            // Line search
            let mut alpha = 1.0;
            let fx = f(&x);
            for _ in 0..30 {
                let mut xn = MathVector::new(1, n);
                for i in 1..=n { xn.set_value(i, x.value(i) + alpha * d.value(i)); }
                if f(&xn) < fx + 1e-4 * alpha * g.dot(&d) { break; }
                alpha *= 0.5;
            }

            // Update x
            let mut s = MathVector::new(1, n);
            for i in 1..=n { s.set_value(i, alpha * d.value(i)); }
            let x_old = x.clone();
            for i in 1..=n { x.set_value(i, x.value(i) + s.value(i)); }

            let g_new = grad(&x);
            let mut y = MathVector::new(1, n);
            for i in 1..=n { y.set_value(i, g_new.value(i) - g.value(i)); }
            g = g_new;

            let sy = s.dot(&y);
            if sy > 1e-15 {
                // BFGS update: H = (I - s*y^T/sy) * H * (I - y*s^T/sy) + s*s^T/sy
                let mut hy = MathVector::new(1, n);
                for i in 1..=n {
                    let mut sum = 0.0;
                    for j in 1..=n { sum += h.value(i, j) * y.value(j); }
                    hy.set_value(i, sum);
                }
                let yhy = y.dot(&hy);
                let factor = 1.0 + yhy / sy;
                for i in 1..=n {
                    for j in 1..=n {
                        let term1 = s.value(i) * s.value(j) * factor / sy;
                        let term2 = (s.value(i) * hy.value(j) + hy.value(i) * s.value(j)) / sy;
                        h.set_value(i, j, h.value(i, j) + term1 - term2);
                    }
                }
            } else {
                // Reset to identity if update is degenerate
                for i in 1..=n { for j in 1..=n { h.set_value(i, j, if i==j {1.0} else {0.0}); } }
            }
        }
        Err(MathStatus::TooManyIterations)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newton_quadratic() {
        let f = |x: &MathVector| -> f64 {
            let a = x.value(1); let b = x.value(2);
            (a - 3.0).powi(2) + (b + 2.0).powi(2)
        };
        let grad = |x: &MathVector| -> MathVector {
            MathVector::from_slice(&[2.0*(x.value(1)-3.0), 2.0*(x.value(2)+2.0)])
        };
        let x0 = MathVector::from_slice(&[0.0, 0.0]);
        let nm = NewtonMinimum::new();
        let x = nm.minimize(f, grad, &x0).unwrap();
        assert!((x.value(1) - 3.0).abs() < 1e-6);
        assert!((x.value(2) + 2.0).abs() < 1e-6);
    }

    #[test]
    fn newton_solve_linear() {
        let f = |x: &MathVector| -> MathVector {
            MathVector::from_slice(&[
                2.0*x.value(1) + x.value(2) - 5.0,
                x.value(1) + 3.0*x.value(2) - 6.0,
            ])
        };
        let x0 = MathVector::from_slice(&[0.0, 0.0]);
        let solver = NewtonSolver::new();
        let x = solver.solve(f, &x0, None).unwrap();
        assert!((x.value(1) - 1.8).abs() < 1e-8);
        assert!((x.value(2) - 1.4).abs() < 1e-8);
    }
}
