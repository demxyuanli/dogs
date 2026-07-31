//! Trigonometric polynomial root finding. Source: `math_TrigonometricFunctionRoots.hxx`
//! Solves a*cos(u) + b*sin(u) + c = 0 for u in [inf, sup].

/// Find roots of a*cos(u) + b*sin(u) + c = 0 in interval [inf, sup].
/// Returns sorted roots (0..2 roots).
pub fn trig_roots(a: f64, b: f64, c: f64, inf: f64, sup: f64) -> Vec<f64> {
    let mut roots = Vec::new();
    // Special case: a = b = 0 → c = 0 is degenerate
    if a.abs() < 1e-30 && b.abs() < 1e-30 {
        if c.abs() < 1e-30 && inf <= sup { return vec![inf]; }
        return roots;
    }
    // R * cos(u - phi) = -c, where R = sqrt(a²+b²), phi = atan2(b, a)
    let r = (a*a + b*b).sqrt();
    let phi = b.atan2(a);
    let rhs = -c / r;
    if rhs < -1.0 || rhs > 1.0 { return roots; }

    let alpha = rhs.acos();
    // u - phi = ± alpha + 2kπ
    let base = phi + alpha;
    let base2 = phi - alpha;
    let two_pi = 2.0 * std::f64::consts::PI;

    // Generate candidates within [inf, sup]
    let mut k_start = ((inf - base) / two_pi).floor() as i32 - 1;
    let mut candidates = Vec::new();
    for k in k_start..k_start + 4 {
        let u = base + k as f64 * two_pi;
        if u >= inf && u <= sup { candidates.push(u); }
    }
    k_start = ((inf - base2) / two_pi).floor() as i32 - 1;
    for k in k_start..k_start + 4 {
        let u = base2 + k as f64 * two_pi;
        if u >= inf && u <= sup { candidates.push(u); }
    }
    candidates.sort_by(|a, b| a.partial_cmp(b).unwrap());
    candidates.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    roots.extend(candidates);
    roots
}

/// Check if any root exists in interval (for classification).
pub fn has_root(a: f64, b: f64, c: f64, inf: f64, sup: f64) -> bool {
    !trig_roots(a, b, c, inf, sup).is_empty()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sine_roots() {
        // sin(u) = 0 → roots at 0, π, 2π in [0, 2π]
        let roots = trig_roots(0.0, 1.0, 0.0, 0.0, 2.0 * std::f64::consts::PI);
        assert_eq!(roots.len(), 3, "got {roots:?}");
        assert!((roots[0]).abs() < 1e-10);
        assert!((roots[1] - std::f64::consts::PI).abs() < 1e-10);
        assert!((roots[2] - 2.0*std::f64::consts::PI).abs() < 1e-10);
    }

    #[test]
    fn cos_shifted() {
        // cos(u) = 0.5 in [0, 2π] → u = π/3, 5π/3
        let roots = trig_roots(1.0, 0.0, -0.5, 0.0, 2.0 * std::f64::consts::PI);
        assert_eq!(roots.len(), 2);
        assert!((roots[0] - std::f64::consts::FRAC_PI_3).abs() < 1e-10);
    }
}
