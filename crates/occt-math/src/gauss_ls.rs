//! Gauss least-squares fitting. Source: `math_GaussLeastSquare`.

use crate::{MathMatrix, MathStatus, MathVector};

/// Least-squares polynomial fit result holder.
pub struct GaussLeastSquare {
    pub is_done: bool,
    coeffs: MathVector,
}

impl GaussLeastSquare {
    /// Access to the fitted coefficients (1-indexed, `coeffs[1]` is the constant term).
    pub fn coeffs(&self) -> &MathVector {
        &self.coeffs
    }
}

/// Fit polynomial `y = c0 + c1*x + ... + cd*x^d` to the sample points by solving
/// the normal equations `(A^T A) c = A^T b`.
///
/// Returns the coefficient vector `c` (1-indexed, `c[1] = c0`). Errors with
/// [`MathStatus::FunctionError`] when there are fewer samples than coefficients
/// or the system is singular.
pub fn least_squares_fit(
    x_vals: &[f64],
    y_vals: &[f64],
    degree: usize,
) -> Result<MathVector, MathStatus> {
    let n = x_vals.len();
    let p = degree + 1;
    if y_vals.len() != n || n < p {
        return Err(MathStatus::FunctionError);
    }

    // Design matrix A (n x p): A[i][j] = x_i^(j-1).
    let mut a = MathMatrix::new(1, n, 1, p);
    for i in 1..=n {
        let mut xp = 1.0;
        for j in 1..=p {
            a.set_value(i, j, xp);
            xp *= x_vals[i - 1];
        }
    }

    // Normal equations: (A^T A) c = A^T b.
    let mut ata = MathMatrix::new(1, p, 1, p);
    let mut atb = MathVector::new(1, p);
    for j in 1..=p {
        for k in 1..=p {
            let mut s = 0.0;
            for i in 1..=n {
                s += a.value(i, j) * a.value(i, k);
            }
            ata.set_value(j, k, s);
        }
        let mut s = 0.0;
        for i in 1..=n {
            s += a.value(i, j) * y_vals[i - 1];
        }
        atb.set_value(j, s);
    }

    ata.solve(&atb).map_err(|_| MathStatus::FunctionError)
}

/// Horner evaluation of the fitted polynomial at `x`.
pub fn evaluate_fit(coeffs: &MathVector, x: f64) -> f64 {
    let n = coeffs.len();
    let mut acc = coeffs.value(n);
    for i in (1..n).rev() {
        acc = acc * x + coeffs.value(i);
    }
    acc
}

/// RMS residual of the fit over the sample points.
pub fn residual_error(coeffs: &MathVector, x_vals: &[f64], y_vals: &[f64]) -> f64 {
    let n = x_vals.len().min(y_vals.len());
    if n == 0 {
        return 0.0;
    }
    let mut s = 0.0;
    for i in 0..n {
        let d = evaluate_fit(coeffs, x_vals[i]) - y_vals[i];
        s += d * d;
    }
    (s / n as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_line() {
        let xs = [0.0, 1.0, 2.0];
        let ys = [1.0, 3.0, 5.0];
        let c = least_squares_fit(&xs, &ys, 1).unwrap();
        assert!((c.value(1) - 1.0).abs() < 1e-9, "intercept {}", c.value(1));
        assert!((c.value(2) - 2.0).abs() < 1e-9, "slope {}", c.value(2));
        assert!(residual_error(&c, &xs, &ys) < 1e-9, "residual");
    }
}
