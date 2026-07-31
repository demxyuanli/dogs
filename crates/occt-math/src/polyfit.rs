//! Polynomial least-squares fit via normal equations.
//! Source: analogue of `math_GaussLeastSquare`.

use crate::{MathMatrix, MathVector};

/// Fit a degree-`degree` polynomial to the points, returning coefficients
/// highest-first (e.g. `[a, b, c]` for `a·x² + b·x + c`).
///
/// Errors when there are fewer samples than coefficients or the normal system
/// is singular.
pub fn polyfit(xs: &[f64], ys: &[f64], degree: usize) -> Result<Vec<f64>, String> {
    let n = xs.len();
    let p = degree + 1;
    if ys.len() != n || n < p {
        return Err("polyfit: need at least degree+1 samples and equal-length xs/ys".to_string());
    }
    // Vandermonde A[i][j] = x_i^(degree+1-j) (highest-first).
    let mut a = MathMatrix::new(1, n, 1, p);
    for i in 1..=n {
        for j in 1..=p {
            let power = degree + 1 - j;
            a.set_value(i, j, xs[i - 1].powi(power as i32));
        }
    }
    // Normal equations: (Aᵀ A) c = Aᵀ b.
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
            s += a.value(i, j) * ys[i - 1];
        }
        atb.set_value(j, s);
    }
    let c = ata.solve(&atb).map_err(|e| e.to_string())?;
    Ok((1..=p).map(|j| c.value(j)).collect())
}

/// Horner evaluation of a polynomial with highest-first coefficients.
pub fn polyval(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().fold(0.0, |acc, &c| acc * x + c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_known_quadratic() {
        // y = 3x² + 2x + 1
        let xs = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        let ys: Vec<f64> = xs.iter().map(|&x| 3.0 * x * x + 2.0 * x + 1.0).collect();
        let c = polyfit(&xs, &ys, 2).unwrap();
        assert_eq!(c.len(), 3);
        assert!((c[0] - 3.0).abs() < 1e-9, "a {}", c[0]);
        assert!((c[1] - 2.0).abs() < 1e-9, "b {}", c[1]);
        assert!((c[2] - 1.0).abs() < 1e-9, "c {}", c[2]);
        assert!((polyval(&c, 1.7) - (3.0 * 1.7 * 1.7 + 2.0 * 1.7 + 1.0)).abs() < 1e-9);
    }

    #[test]
    fn insufficient_samples() {
        assert!(polyfit(&[0.0, 1.0], &[0.0, 1.0], 2).is_err());
    }
}
