//! Adaptive and composite 1-D quadrature.
//! Source: `math_Kronrod`, `math_Gauss`.

/// Recursive adaptive Simpson quadrature with a depth limit of 30.
pub fn adaptive_simpson<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64) -> f64 {
    adaptive_simpson_rec(f, a, b, tol, 0)
}

fn adaptive_simpson_rec<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64, depth: usize) -> f64 {
    let c = 0.5 * (a + b);
    let (fa, fc, fb) = (f(a), f(c), f(b));
    let whole = (b - a) / 6.0 * (fa + 4.0 * fc + fb);

    let d = 0.5 * (a + c);
    let e = 0.5 * (c + b);
    let fd = f(d);
    let fe = f(e);
    let left = (c - a) / 6.0 * (fa + 4.0 * fd + fc);
    let right = (b - c) / 6.0 * (fc + 4.0 * fe + fb);

    let delta = left + right - whole;
    if depth >= 30 || delta.abs() <= 15.0 * tol {
        // Richardson-extrapolated value (classic adaptive Simpson).
        left + right + delta / 15.0
    } else {
        adaptive_simpson_rec(f, a, c, tol / 2.0, depth + 1)
            + adaptive_simpson_rec(f, c, b, tol / 2.0, depth + 1)
    }
}

/// Romberg integration: Richardson extrapolation over `n` refinement levels
/// of the trapezoid rule. Returns the `[n][n]` entry of the Romberg table.
pub fn romberg<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut r = vec![vec![0.0f64; n + 1]; n + 1];
    let mut h = b - a;
    r[0][0] = 0.5 * h * (f(a) + f(b));
    for i in 1..=n {
        h *= 0.5;
        let mut sum = 0.0;
        let mut x = a + h;
        while x < b {
            sum += f(x);
            x += 2.0 * h;
        }
        r[i][0] = 0.5 * r[i - 1][0] + h * sum;
        for j in 1..=i {
            let factor = 4.0_f64.powi(j as i32);
            r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (factor - 1.0);
        }
    }
    r[n][n]
}

/// Composite trapezoid rule with `n` sub-intervals.
pub fn composite_trapezoid<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let h = (b - a) / n as f64;
    let mut sum = 0.5 * (f(a) + f(b));
    for i in 1..n {
        sum += f(a + i as f64 * h);
    }
    sum * h
}

/// Composite midpoint rule with `n` sub-intervals.
pub fn midpoint_rule<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let h = (b - a) / n as f64;
    (0..n).map(|i| f(a + (i as f64 + 0.5) * h)).sum::<f64>() * h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adaptive_simpson_sin() {
        // ∫₀^π sin(x) dx = 2
        let v = adaptive_simpson(&|x| x.sin(), 0.0, std::f64::consts::PI, 1e-8);
        assert!((v - 2.0).abs() < 1e-8, "v = {v}");
    }

    #[test]
    fn romberg_sin() {
        // ∫₀^π sin(x) dx = 2
        let v = romberg(&|x| x.sin(), 0.0, std::f64::consts::PI, 10);
        assert!((v - 2.0).abs() < 1e-10, "v = {v}");
    }

    #[test]
    fn trapezoid_x2() {
        // ∫₀¹ x² dx = 1/3
        let v = composite_trapezoid(&|x| x * x, 0.0, 1.0, 100);
        assert!((v - 1.0 / 3.0).abs() < 1e-4, "v = {v}");
    }

    #[test]
    fn midpoint_x2() {
        let v = midpoint_rule(&|x| x * x, 0.0, 1.0, 100);
        assert!((v - 1.0 / 3.0).abs() < 1e-4, "v = {v}");
    }
}
