//! Eigenvalue searcher (power iteration + QR for symmetric matrices).
//! Source: `math_EigenValuesSearcher.hxx`
use crate::{MathVector, MathMatrix};

/// Compute the largest eigenvalue and eigenvector via power iteration.
pub fn largest_eigen(a: &MathMatrix, x0: &MathVector, max_iter: usize, tol: f64) -> (f64, MathVector) {
    let n = a.row_count();
    let mut x = x0.clone();
    let mut lambda = 0.0f64;

    for _ in 0..max_iter {
        let norm = x.norm();
        if norm < 1e-30 { break; }
        for i in 1..=n { x.set_value(i, x.value(i) / norm); }

        // y = A * x
        let mut y = MathVector::new(1, n);
        for i in 1..=n {
            let mut s = 0.0f64;
            for j in 1..=n { s += a.value(i, j) * x.value(j); }
            y.set_value(i, s);
        }
        // Rayleigh quotient
        let new_lambda = x.dot(&y);
        if (new_lambda - lambda).abs() < tol * (1.0 + new_lambda.abs()) {
            lambda = new_lambda;
            x = y;
            break;
        }
        lambda = new_lambda;
        x = y;
    }
    let norm = x.norm();
    for i in 1..=n { x.set_value(i, x.value(i) / norm.max(1e-30)); }
    (lambda, x)
}

/// Compute smallest eigenvalue (via power iteration on inverse — shifted).
pub fn smallest_eigen(a: &MathMatrix, x0: &MathVector, max_iter: usize, tol: f64) -> (f64, MathVector) {
    let inv = a.inverted().ok();
    let n = a.row_count();
    match inv {
        Some(ai) => {
            let (lambda_inv, v) = largest_eigen(&ai, x0, max_iter, tol);
            (1.0 / lambda_inv.max(1e-30), v)
        }
        None => (0.0, MathVector::new(1, n)),
    }
}

/// Compute all eigenvalues of a symmetric matrix via Jacobi iteration.
pub fn all_eigenvalues(a: &MathMatrix) -> Vec<f64> {
    let jac = crate::Jacobi::new(a);
    (1..=jac.values.len()).map(|i| jac.values.value(i + jac.values.lower() - 1)).collect()
}

/// Condition number = |largest| / |smallest| eigenvalue.
pub fn condition_number(a: &MathMatrix) -> f64 {
    let n = a.row_count();
    let x0 = MathVector::with_init(1, n, 1.0);
    let (lmax, _) = largest_eigen(a, &x0, 100, 1e-10);
    let (lmin, _) = smallest_eigen(a, &x0, 100, 1e-10);
    (lmax.abs() / lmin.abs().max(1e-30)).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn largest_eigen_diag() {
        let mut a = MathMatrix::new(1, 2, 1, 2);
        a.set_value(1,1,5.0); a.set_value(2,2,3.0);
        let x0 = MathVector::from_slice(&[1.0, 1.0]);
        let (lam, v) = largest_eigen(&a, &x0, 50, 1e-10);
        assert!((lam - 5.0).abs() < 0.1, "got {lam}");
        // eigenvector should be dominated by the x-component (eigenvalue 5 > 3)
        assert!(v.value(2).abs() < 0.2);
    }
}
