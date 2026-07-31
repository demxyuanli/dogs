//! Polynomial root finding: Laguerre's method, Newton grid search, deflation.
//! Source: `math_FunctionRoot`, `math_FunctionRoots`.

use crate::polynomial::{poly_derivative, poly_eval};

/// Find one real root of a polynomial (highest-first coefficients) with
/// Laguerre's method, starting from `initial`.
pub fn laguerre(coeffs: &[f64], initial: f64, tol: f64, max_iter: usize) -> Result<f64, String> {
    if coeffs.is_empty() {
        return Err("laguerre: empty polynomial".to_string());
    }
    let n = coeffs.len() - 1;
    if n == 0 {
        return Err("laguerre: constant polynomial has no root".to_string());
    }
    if coeffs[0].abs() < 1e-300 {
        return Err("laguerre: leading coefficient is zero".to_string());
    }
    let nf = n as f64;
    let deriv = poly_derivative(coeffs);
    let deriv2 = poly_derivative(&deriv);
    let mut x = initial;

    for _ in 0..max_iter {
        let p = poly_eval(coeffs, x);
        if p.abs() < tol {
            return Ok(x);
        }
        let dp = poly_eval(&deriv, x);
        if dp.abs() < 1e-300 {
            // Flat spot: nudge and continue rather than dividing by zero.
            x += 1e-6 * (1.0 + x.abs());
            continue;
        }
        let g = dp / p;
        let ddp = poly_eval(&deriv2, x);
        let h = g * g - ddp / p;
        let rad = nf * h - g * g;
        // For complex conjugate pairs rad < 0; fall back to a scaled step.
        let sq = if rad >= 0.0 { ((nf - 1.0) * rad).sqrt() } else { 0.0 };
        let d1 = g + sq;
        let d2 = g - sq;
        let denom = if d1.abs() > d2.abs() { d1 } else { d2 };
        if denom.abs() < 1e-300 {
            x += 1e-6 * (1.0 + x.abs());
            continue;
        }
        let step = nf / denom;
        x -= step;
        if step.abs() < tol {
            return Ok(x);
        }
    }
    Err("laguerre: max iterations reached".to_string())
}

/// All real roots of a polynomial (highest-first coefficients) by running
/// Newton's method from a grid of starting points over a Cauchy-bounded
/// interval, then de-duplicating within `tol`.
pub fn poly_roots_newton(coeffs: &[f64], tol: f64) -> Vec<f64> {
    let n = coeffs.len().saturating_sub(1);
    if n == 0 || coeffs[0].abs() < 1e-300 {
        return Vec::new();
    }
    // Cauchy's bound: every root lies in |x| <= 1 + max(|a_i|)/|a_n|.
    let max_c = coeffs[1..].iter().fold(0.0f64, |m, &c| m.max(c.abs()));
    let r = 1.0 + max_c / coeffs[0].abs();
    let deriv = poly_derivative(coeffs);

    let mut roots: Vec<f64> = Vec::new();
    let grid = 64usize;
    for i in 0..grid {
        let t = (i as f64 + 0.5) / grid as f64;
        let x0 = -r + 2.0 * r * t;
        if let Some(x) = newton_from(coeffs, &deriv, x0, tol) {
            if poly_eval(coeffs, x).abs() < tol && !roots.iter().any(|&rt| (rt - x).abs() <= tol) {
                roots.push(x);
            }
        }
    }
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots
}

/// Newton iteration from `x0`; returns the converged point or `None` on
/// divergence / non-finite values.
fn newton_from(coeffs: &[f64], deriv: &[f64], x0: f64, tol: f64) -> Option<f64> {
    let mut x = x0;
    for _ in 0..100 {
        let fx = poly_eval(coeffs, x);
        if fx.abs() < tol {
            return Some(x);
        }
        let dfx = poly_eval(deriv, x);
        if dfx.abs() < 1e-300 {
            return None;
        }
        let dx = fx / dfx;
        x -= dx;
        if !x.is_finite() {
            return None;
        }
        if dx.abs() < tol {
            return Some(x);
        }
    }
    None
}

/// Synthetic division of `coeffs` (highest-first) by `(x - root)`; returns the
/// quotient coefficients (highest-first).
pub fn deflate(coeffs: &[f64], root: f64) -> Vec<f64> {
    let n = coeffs.len();
    if n <= 1 {
        return Vec::new();
    }
    let mut b = vec![0.0f64; n - 1];
    let mut carry = 0.0;
    for (i, &c) in coeffs[..n - 1].iter().enumerate() {
        carry = c + root * carry;
        b[i] = carry;
    }
    b
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::poly_eval;

    #[test]
    fn laguerre_sqrt2() {
        // x² - 2 = 0
        let r = laguerre(&[1.0, 0.0, -2.0], 1.5, 1e-12, 100).unwrap();
        assert!((r - 2.0_f64.sqrt()).abs() < 1e-10, "r = {r}");
        let r = laguerre(&[1.0, 0.0, -2.0], -1.5, 1e-12, 100).unwrap();
        assert!((r + 2.0_f64.sqrt()).abs() < 1e-10, "r = {r}");
    }

    #[test]
    fn laguerre_multiple_root() {
        // x³ - 3x² + 3x - 1 = (x-1)³ has a flat derivative at the root.
        let r = laguerre(&[1.0, -3.0, 3.0, -1.0], 0.0, 1e-12, 200).unwrap();
        assert!((r - 1.0).abs() < 1e-10, "r = {r}");
    }

    #[test]
    fn roots_newton_square() {
        // x² - 4 = 0 -> {-2, 2}
        let mut r = poly_roots_newton(&[1.0, 0.0, -4.0], 1e-9);
        r.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(r.len(), 2, "roots: {r:?}");
        assert!((r[0] + 2.0).abs() < 1e-6, "r[0] = {}", r[0]);
        assert!((r[1] - 2.0).abs() < 1e-6, "r[1] = {}", r[1]);
    }

    #[test]
    fn deflate_divides() {
        // (x²-4)/(x-2) = x+2
        let q = deflate(&[1.0, 0.0, -4.0], 2.0);
        assert_eq!(q, vec![1.0, 2.0]);
        // x²-3x+2 = (x-1)(x-2): dividing by (x-1) leaves the other root.
        let q = deflate(&[1.0, -3.0, 2.0], 1.0);
        assert_eq!(q, vec![1.0, -2.0]);
        assert!((poly_eval(&q, 2.0)).abs() < 1e-15);
    }
}
