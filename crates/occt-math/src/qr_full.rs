//! Full QR decomposition via Gram-Schmidt and least-squares solve.
//! Source: analogue of `math_Householder` using (modified) Gram-Schmidt.

use crate::{MathMatrix, MathVector};

/// Thin QR decomposition `A = Q·R` via modified Gram-Schmidt.
///
/// `A` is `m×n`; `Q` is `m×p` and `R` is `p×n` with `p = min(m, n)`.
/// For full-rank `A` the reconstruction `Q·R` reproduces `A` exactly.
pub fn gram_schmidt(a: &MathMatrix) -> (MathMatrix, MathMatrix) {
    let m = a.row_count();
    let n = a.col_count();
    let p = m.min(n);
    let mut q = MathMatrix::new(1, m, 1, p);
    let mut r = MathMatrix::new(1, p, 1, n);
    for j in 1..=n {
        // Work on column j of A.
        let mut v: Vec<f64> = (1..=m).map(|i| a.value(i, j)).collect();
        for i in 1..=p.min(j - 1) {
            let mut dot = 0.0;
            for k in 0..m {
                dot += q.value(k + 1, i) * v[k];
            }
            r.set_value(i, j, dot);
            for k in 0..m {
                v[k] -= dot * q.value(k + 1, i);
            }
        }
        if j <= p {
            let nrm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            r.set_value(j, j, nrm);
            if nrm > 1e-12 {
                for k in 0..m {
                    q.set_value(k + 1, j, v[k] / nrm);
                }
            }
            // rank-deficient column: leave Q column zero, R[j][j] stays ~0
        }
    }
    (q, r)
}

/// Solve the square system `A x = b` via QR back-substitution.
/// Returns NaN components if `A` is numerically rank-deficient.
pub fn qr_solve(a: &MathMatrix, b: &MathVector) -> MathVector {
    let n = a.row_count();
    let (q, r) = gram_schmidt(a);
    // y = Q^T b
    let mut y = MathVector::new(1, n);
    for i in 1..=n {
        let mut s = 0.0;
        for k in 1..=n {
            s += q.value(k, i) * b.value(k);
        }
        y.set_value(i, s);
    }
    // R x = y
    let mut x = MathVector::new(1, n);
    for i in (1..=n).rev() {
        let mut s = y.value(i);
        for j in (i + 1)..=n {
            s -= r.value(i, j) * x.value(j);
        }
        let diag = r.value(i, i);
        x.set_value(i, if diag.abs() < 1e-12 { f64::NAN } else { s / diag });
    }
    x
}

/// Least-squares minimizer of `||A x - b||` for an `m×n` system (usually `m ≥ n`).
/// Returns NaN components if `A` is numerically rank-deficient.
pub fn qr_least_squares(a: &MathMatrix, b: &MathVector) -> MathVector {
    let m = a.row_count();
    let n = a.col_count();
    let p = m.min(n);
    let (q, r) = gram_schmidt(a);
    // y = Q^T b (only the p leading components matter)
    let mut y = MathVector::new(1, p);
    for i in 1..=p {
        let mut s = 0.0;
        for k in 1..=m {
            s += q.value(k, i) * b.value(k);
        }
        y.set_value(i, s);
    }
    // Back-substitute the upper-triangular R (p×n, take the first p rows).
    let mut x = MathVector::new(1, n);
    for i in (1..=p).rev() {
        let mut s = y.value(i);
        for j in (i + 1)..=n {
            s -= r.value(i, j) * x.value(j);
        }
        let diag = r.value(i, i);
        x.set_value(i, if diag.abs() < 1e-12 { f64::NAN } else { s / diag });
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qr_reconstruct() {
        let mut a = MathMatrix::new(1, 3, 1, 2);
        a.set_value(1, 1, 1.0);
        a.set_value(1, 2, 1.0);
        a.set_value(2, 1, 0.0);
        a.set_value(2, 2, 1.0);
        a.set_value(3, 1, 0.0);
        a.set_value(3, 2, 1.0);
        let (q, r) = gram_schmidt(&a);
        let qr = q.multiplied_mat(&r);
        assert!((qr.value(1, 1) - a.value(1, 1)).abs() < 1e-12);
        assert!((qr.value(2, 2) - a.value(2, 2)).abs() < 1e-12);
        assert!((qr.value(3, 2) - a.value(3, 2)).abs() < 1e-12);
    }

    #[test]
    fn qr_solve_square() {
        let mut a = MathMatrix::new(1, 2, 1, 2);
        a.set_value(1, 1, 3.0);
        a.set_value(1, 2, 1.0);
        a.set_value(2, 1, 1.0);
        a.set_value(2, 2, 2.0);
        let b = MathVector::from_slice(&[5.0, 5.0]); // x = [1,2]
        let x = qr_solve(&a, &b);
        assert!((x.value(1) - 1.0).abs() < 1e-12);
        assert!((x.value(2) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn least_squares_line() {
        // Fit y = 2x + 1 from three exact points.
        let mut a = MathMatrix::new(1, 3, 1, 2);
        a.set_value(1, 1, 0.0);
        a.set_value(1, 2, 1.0);
        a.set_value(2, 1, 1.0);
        a.set_value(2, 2, 1.0);
        a.set_value(3, 1, 2.0);
        a.set_value(3, 2, 1.0);
        let b = MathVector::from_slice(&[1.0, 3.0, 5.0]);
        let x = qr_least_squares(&a, &b);
        assert!((x.value(1) - 2.0).abs() < 1e-8, "slope {}", x.value(1));
        assert!((x.value(2) - 1.0).abs() < 1e-8, "intercept {}", x.value(2));
    }
}
