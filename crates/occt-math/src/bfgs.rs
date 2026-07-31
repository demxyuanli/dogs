//! BFGS quasi-Newton optimization. Source: `math_BFGS.cxx`
use crate::{MathVector, MathMatrix, MathStatus};

/// BFGS optimizer for unconstrained minimization.
pub struct BFGS { pub max_iter: usize, pub tolerance: f64 }
impl BFGS {
    pub fn new() -> Self { Self { max_iter: 200, tolerance: 1e-8 } }

    /// Minimize f(x) starting from x0. grad must return gradient at x.
    pub fn minimize<F, G>(&self, f: F, grad: G, x0: &MathVector) -> Result<MathVector, MathStatus>
    where F: Fn(&MathVector) -> f64, G: Fn(&MathVector) -> MathVector
    {
        let n = x0.len();
        let mut x = x0.clone();
        let mut g = grad(&x);
        let mut h = MathMatrix::new(1, n, 1, n); // inverse Hessian approx
        for i in 1..=n { h.set_value(i, i, 1.0); }

        for _iter in 0..self.max_iter {
            if g.norm() < self.tolerance { return Ok(x); }
            // Search direction d = -H * g
            let mut d = MathVector::new(1, n);
            for i in 1..=n {
                let mut s = 0.0f64;
                for j in 1..=n { s += h.value(i, j) * g.value(j); }
                d.set_value(i, -s);
            }
            // Line search with backtracking
            let mut alpha = 1.0f64;
            let f0 = f(&x);
            let dg = g.dot(&d);
            for _ls in 0..30 {
                let mut xn = MathVector::new(1, n);
                for i in 1..=n { xn.set_value(i, x.value(i) + alpha * d.value(i)); }
                if f(&xn) < f0 + 1e-4 * alpha * dg { break; }
                alpha *= 0.5;
            }
            // Update x
            let mut s = MathVector::new(1, n);
            for i in 1..=n { s.set_value(i, alpha * d.value(i)); }
            let x_old = x.clone();
            for i in 1..=n { x.set_value(i, x.value(i) + s.value(i)); }
            let g_new = grad(&x);
            // y = g_new - g
            let mut y = MathVector::new(1, n);
            for i in 1..=n { y.set_value(i, g_new.value(i) - g.value(i)); }
            g = g_new;
            let sy = s.dot(&y);
            if sy > 1e-15 {
                let mut hy = MathVector::new(1, n);
                for i in 1..=n { let mut sum=0.0f64; for j in 1..=n { sum+=h.value(i,j)*y.value(j); } hy.set_value(i,sum); }
                let yhy = y.dot(&hy);
                let coef = 1.0 + yhy / sy;
                for i in 1..=n { for j in 1..=n {
                    h.set_value(i, j, h.value(i,j) + (coef * s.value(i) * s.value(j) - s.value(i) * hy.value(j) - hy.value(i) * s.value(j)) / sy);
                }}
            } else {
                for i in 1..=n { for j in 1..=n { h.set_value(i, j, if i==j {1.0} else {0.0}); } }
            }
        }
        Err(MathStatus::TooManyIterations)
    }
}

##[cfg(test)]
#mod tests {
#    use super::*;
#    #[test]
#    fn bfgs_quadratic() {
#        let f = |x: &MathVector| -> f64 { (x.value(1)-3.).powi(2)+(x.value(2)+2.).powi(2) };
#        let grad = |x: &MathVector| -> MathVector { MathVector::from_slice(&[2.*(x.value(1)-3.),2.*(x.value(2)+2.)]) };
#        let b = BFGS::new();
#        let x = b.minimize(f, grad, &MathVector::from_slice(&[0.,0.])).unwrap();
#        assert!((x.value(1)-3.).abs()<1e-6); assert!((x.value(2)+2.).abs()<1e-6);
#    }
#}
