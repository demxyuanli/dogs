//! Polynomial arithmetic, closed-form quadratic/cubic roots, and
//! Durand-Kerner root finding for general polynomials.
//! Source: `math_DirectPolynomialRoots`.

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

/// Real roots of `a·x² + b·x + c = 0`, `None` if no real roots or not quadratic.
pub fn quadratic_roots(a: f64, b: f64, c: f64) -> Option<(f64, f64)> {
    if a.abs() < 1e-300 {
        return None;
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    // Numerically stable form.
    let q = -0.5 * (b + b.signum() * sq);
    let r1 = q / a;
    let r2 = if q.abs() > 1e-300 { c / q } else { (-b - sq) / (2.0 * a) };
    Some((r1, r2))
}

/// Real roots of `a·x³ + b·x² + c·x + d = 0` (closed form), 1 to 3 roots.
pub fn cubic_roots(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a.abs() < 1e-300 {
        return quadratic_roots(b, c, d).map(|(r1, r2)| vec![r1, r2]).unwrap_or_default();
    }
    let b = b / a;
    let c = c / a;
    let d = d / a;
    // Depressed cubic: t³ + p·t + q with x = t - b/3.
    let p = c - b * b / 3.0;
    let q = 2.0 * b * b * b / 27.0 - b * c / 3.0 + d;
    let disc = q * q / 4.0 + p * p * p / 27.0;
    let shift = b / 3.0;
    let mut roots = Vec::new();
    if disc > 1e-14 {
        // One real root.
        let sq = disc.sqrt();
        let u = (-q / 2.0 + sq).cbrt();
        let v = (-q / 2.0 - sq).cbrt();
        roots.push(u + v - shift);
    } else if disc >= -1e-14 {
        // Coincident real roots.
        if q.abs() < 1e-14 {
            roots.push(-shift); // triple root
        } else {
            let u = (-q / 2.0).cbrt();
            roots.push(2.0 * u - shift);
            roots.push(-u - shift);
        }
    } else {
        // Three distinct real roots (trigonometric form).
        let r = 2.0 * (-p / 3.0).sqrt();
        let theta = ((3.0 * q / (2.0 * p)) * (-3.0 / p).sqrt()).acos() / 3.0;
        let two_pi = 2.0 * std::f64::consts::PI;
        roots.push(r * theta.cos() - shift);
        roots.push(r * (theta - two_pi / 3.0).cos() - shift);
        roots.push(r * (theta + two_pi / 3.0).cos() - shift);
    }
    // De-duplicate near-equal roots and sort ascending.
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut out: Vec<f64> = Vec::new();
    for x in roots {
        if out.last().map_or(true, |&last| (x - last).abs() > 1e-9 * (1.0 + x.abs())) {
            out.push(x);
        }
    }
    out
}

/// Roots of a real polynomial with highest-first coefficients.
///
/// Degrees 1-3 use the closed forms; higher degrees use Durand-Kerner iteration.
/// Returns the real parts of all roots (complex-conjugate pairs included twice).
pub fn polynomial_roots(coeffs: &[f64]) -> Vec<f64> {
    match coeffs.len() {
        0 | 1 => Vec::new(),
        2 => {
            if coeffs[0].abs() < 1e-300 {
                Vec::new()
            } else {
                vec![-coeffs[1] / coeffs[0]]
            }
        }
        3 => quadratic_roots(coeffs[0], coeffs[1], coeffs[2])
            .map(|(r1, r2)| vec![r1, r2])
            .unwrap_or_default(),
        4 => cubic_roots(coeffs[0], coeffs[1], coeffs[2], coeffs[3]),
        _ => durand_kerner(coeffs),
    }
}

/// Durand-Kerner iteration for a real polynomial of degree >= 4.
fn durand_kerner(coeffs: &[f64]) -> Vec<f64> {
    let n = coeffs.len() - 1;
    let lead = coeffs[0];
    let monic: Vec<f64> = coeffs.iter().map(|&c| c / lead).collect();

    // Initial guesses (0.4 + 0.9j)^k, k = 0..n-1.
    let seed = C { re: 0.4, im: 0.9 };
    let mut roots: Vec<C> = Vec::with_capacity(n);
    let mut z = C { re: 1.0, im: 0.0 };
    for _ in 0..n {
        roots.push(z);
        z = z.mul(seed);
    }

    for _ in 0..1000 {
        let mut max_delta = 0.0_f64;
        for i in 0..n {
            let num = poly_eval_c(&monic, roots[i]);
            let mut den = C { re: 1.0, im: 0.0 };
            for j in 0..n {
                if i != j {
                    den = den.mul(roots[i].sub(roots[j]));
                }
            }
            let delta = num.div(den);
            roots[i] = roots[i].sub(delta);
            max_delta = max_delta.max(delta.norm());
        }
        if max_delta < 1e-12 {
            break;
        }
    }
    roots.into_iter().map(|r| r.re).collect()
}

/// Minimal complex arithmetic used by [`durand_kerner`].
#[derive(Debug, Clone, Copy)]
struct C {
    re: f64,
    im: f64,
}

impl C {
    fn norm(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
    fn add(self, o: C) -> C {
        C { re: self.re + o.re, im: self.im + o.im }
    }
    fn sub(self, o: C) -> C {
        C { re: self.re - o.re, im: self.im - o.im }
    }
    fn mul(self, o: C) -> C {
        C { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
    fn div(self, o: C) -> C {
        let den = o.re * o.re + o.im * o.im;
        C { re: (self.re * o.re + self.im * o.im) / den, im: (self.im * o.re - self.re * o.im) / den }
    }
}

/// Horner evaluation at a complex point.
fn poly_eval_c(coeffs: &[f64], x: C) -> C {
    coeffs
        .iter()
        .fold(C { re: 0.0, im: 0.0 }, |acc, &c| acc.mul(x).add(C { re: c, im: 0.0 }))
}

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

    #[test]
    fn durand_kerner_quartic() {
        // (x²-1)(x²-4) = x⁴ - 5x² + 4 -> ±1, ±2
        let mut r = polynomial_roots(&[1.0, 0.0, -5.0, 0.0, 4.0]);
        r.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(r.len(), 4);
        for (got, want) in r.iter().zip([-2.0, -1.0, 1.0, 2.0]) {
            assert!((got - want).abs() < 1e-6, "root {got} vs {want}");
        }
    }
}
