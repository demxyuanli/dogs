//! Polynomial arithmetic helpers (Horner evaluation and derivative).
//!
//! Root finding for degrees 1 to 4 is not implemented here: it lives in the
//! ported `math_DirectPolynomialRoots`
//! (`occt_core::math_direct_poly_roots`), which `quadratic_roots` and
//! `cubic_roots` below delegate to.

/// Horner evaluation of a polynomial with highest-first coefficients.
pub fn poly_eval(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().fold(0.0, |acc, &c| acc * x + c)
}

/// Derivative coefficients (highest-first); the constant term drops out.
pub fn poly_derivative(coeffs: &[f64]) -> Vec<f64> {
    let n = coeffs.len();
    let mut out = Vec::with_capacity(n.saturating_sub(1));
    for i in 0..n.saturating_sub(1) {
        out.push((n - 1 - i) as f64 * coeffs[i]);
    }
    out
}

/// Real roots of `a*x^2 + b*x + c = 0` as `(Value(1), Value(2))`.
/// `None` when there is no real root. Delegates to the ported
/// `math_DirectPolynomialRoots` (`occt_core::math_direct_poly_roots`); the
/// closed form with its `1e-14` discriminant comparisons used to live here.
pub fn quadratic_roots(a: f64, b: f64, c: f64) -> Option<(f64, f64)> {
    let r = occt_core::math_direct_poly_roots::direct_quadratic_roots(a, b, c);
    if r.len() == 2 {
        Some((r[0], r[1]))
    } else {
        None
    }
}

/// Real roots of `a*x^3 + b*x^2 + c*x + d = 0` in `Value` order, 1 to 3
/// roots. Delegates to the ported `math_DirectPolynomialRoots`.
pub fn cubic_roots(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    occt_core::math_direct_poly_roots::direct_cubic_roots(a, b, c, d)
}

// ---------------------------------------------------------------------------
// The `polynomial_roots` / `durand_kerner` pair that used to live here was an
// invented general-degree solver: `math_DirectPolynomialRoots` stops at degree
// 4 and OCCT has no counterpart for it, so it is gone rather than kept as a
// second root-finding path.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eval_and_derivative() {
        // p(x) = 3x³ + x² - 2x + 5
        let c = [3.0, 1.0, -2.0, 5.0];
        assert!((poly_eval(&c, 2.0) - (24.0 + 4.0 - 4.0 + 5.0)).abs() < 1e-12);
        let d = poly_derivative(&c);
        assert_eq!(d, vec![9.0, 2.0, -2.0]);
    }

    #[test]
    fn quadratic_roots_stable() {
        // x² - 5x + 6 = 0 -> 2, 3
        let (r1, r2) = quadratic_roots(1.0, -5.0, 6.0).unwrap();
        assert!((r1.min(r2) - 2.0).abs() < 1e-12);
        assert!((r1.max(r2) - 3.0).abs() < 1e-12);
        assert!(quadratic_roots(1.0, 0.0, 1.0).is_none()); // no real roots
    }

    #[test]
    fn cubic_roots_closed_form() {
        // (x-1)(x-2)(x-3) = x³ - 6x² + 11x - 6
        let mut r = cubic_roots(1.0, -6.0, 11.0, -6.0);
        r.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(r.len(), 3);
        assert!((r[0] - 1.0).abs() < 1e-10);
        assert!((r[1] - 2.0).abs() < 1e-10);
        assert!((r[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn cubic_roots_single() {
        // x³ - 8 = 0 -> 2 (one real root)
        let r = cubic_roots(1.0, 0.0, 0.0, -8.0);
        assert_eq!(r.len(), 1);
        assert!((r[0] - 2.0).abs() < 1e-10);
    }
}
