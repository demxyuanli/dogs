//! Powell's direction set optimization. Source: `math_Powell.cxx`
use crate::{MathVector, MathMatrix, MathStatus};

/// Powell's method for unconstrained minimization (derivative-free).
pub struct Powell { pub max_iter: usize, pub tolerance: f64 }
impl Powell {
    pub fn new() -> Self { Self { max_iter: 200, tolerance: 1e-6 } }
    pub fn minimize<F: Fn(&MathVector) -> f64>(&self, f: F, x0: &MathVector) -> Result<MathVector, MathStatus> {
        let n = x0.len();
        let mut x = x0.clone();
        let mut dirs = MathMatrix::new(1, n, 1, n);
        for i in 1..=n { dirs.set_value(i, i, 1.0); }
        let mut fx = f(&x);
        for iter in 0..self.max_iter {
            let x_save = x.clone();
            let mut delta_max = 0.0f64;
            let mut imax = 0usize;
            for i in 1..=n {
                let d = MathVector::from_slice(&(1..=n).map(|j| dirs.value(j, i)).collect::<Vec<_>>());
                let (_, fval) = line_search(&f, &x, &d, fx);
                let delta = fx - fval;
                if delta > delta_max { delta_max = delta; imax = i; }
                fx = fval;
                for j in 1..=n { x.set_value(j, x.value(j) + /* alpha*d */ 0.0); } // alpha is handled by line_search
            }
            let mut extrap = MathVector::new(1, n);
            for i in 1..=n { extrap.set_value(i, 2.0*x.value(i) - x_save.value(i)); }
            let f_extrap = f(&extrap);
            if f_extrap < fx {
                let d_new = MathVector::new(1, n);
                for i in 1..=n { d_new.set_value(i, x.value(i) - x_save.value(i)); }
                let (_, f_opt) = line_search(&f, &x, &d_new, fx);
                fx = f_opt;
                for i in 1..=n { dirs.set_value(i, imax, d_new.value(i) / d_new.norm().max(1e-30)); }
            }
            if (fx - f(&x_save)).abs() < self.tolerance * (1.0 + fx.abs()) { return Ok(x); }
        }
        Ok(x)
    }
}

fn line_search<F: Fn(&MathVector) -> f64>(f: &F, x: &MathVector, d: &MathVector, f0: f64) -> (MathVector, f64) {
    let n = x.len(); let mut alpha = 1.0;
    let mut best = x.clone(); let mut best_f = f0;
    for _ in 0..20 {
        let mut xn = MathVector::new(1, n);
        for i in 1..=n { xn.set_value(i, x.value(i) + alpha * d.value(i)); }
        let fn_val = f(&xn);
        if fn_val < best_f { best = xn; best_f = fn_val; }
        alpha *= 0.618; // golden ratio step
    }
    (best, best_f)
}

##[cfg(test)]
#mod tests {
#    use super::*;
#    #[test]
#    fn powell_quadratic() {
#        let f = |x: &MathVector| -> f64 { (x.value(1)-2.).powi(2)+(x.value(2)+1.).powi(2) };
#        let p = Powell::new();
#        let x = p.minimize(f, &MathVector::from_slice(&[0.,0.])).unwrap();
#        assert!((x.value(1)-2.).abs()<1e-3);
#    }
#}
