//! 1-D interpolation: linear, nearest-neighbor, and Newton divided-difference.
//! Source: `math_Interpolation`, `math_NewtonFunctionRoot`.

/// Interpolation strategy selected when building an [`Interpolant`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterpKind {
    Linear,
    Nearest,
    Newton,
}

/// A ready-to-evaluate interpolant over a sample set.
#[derive(Debug, Clone)]
pub struct Interpolant {
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub kind: InterpKind,
    /// Newton divided-difference coefficients when `kind == Newton`.
    pub coefficients: Option<Vec<f64>>,
}

/// Validate sample arrays: same non-empty length, strictly increasing xs.
fn validate(xs: &[f64], ys: &[f64]) -> Result<(), String> {
    if xs.is_empty() || xs.len() != ys.len() {
        return Err("interp: xs and ys must be non-empty and equal length".to_string());
    }
    if xs.windows(2).any(|w| w[1] <= w[0]) {
        return Err("interp: xs must be strictly increasing".to_string());
    }
    Ok(())
}

/// Piecewise-linear interpolation. Errors if `xs` is unsorted or `x` is outside
/// the sample range.
pub fn linear_interp(xs: &[f64], ys: &[f64], x: f64) -> Result<f64, String> {
    validate(xs, ys)?;
    let (lo, hi) = (xs[0], xs[xs.len() - 1]);
    if x < lo || x > hi {
        return Err("linear_interp: x out of range".to_string());
    }
    // Find the bracketing interval [xs[i], xs[i+1]].
    let mut i = 0;
    while i + 1 < xs.len() && x > xs[i + 1] {
        i += 1;
    }
    if xs[i] == x {
        return Ok(ys[i]);
    }
    let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    Ok(ys[i] + t * (ys[i + 1] - ys[i]))
}

/// Interpolate by selecting the sample whose `xs` is closest to `x`.
pub fn nearest_interp(xs: &[f64], ys: &[f64], x: f64) -> Result<f64, String> {
    validate(xs, ys)?;
    let mut best = 0;
    let mut best_d = (x - xs[0]).abs();
    for i in 1..xs.len() {
        let d = (x - xs[i]).abs();
        if d < best_d {
            best_d = d;
            best = i;
        }
    }
    Ok(ys[best])
}

/// Newton divided-difference coefficients of the interpolating polynomial.
/// `coeffs[k]` multiplies the basis term `∏_{i<k} (x - xs[i])`.
pub fn newton_coefficients(xs: &[f64], ys: &[f64]) -> Result<Vec<f64>, String> {
    validate(xs, ys)?;
    let n = xs.len();
    let mut f: Vec<f64> = ys.to_vec();
    let mut coeffs = vec![f[0]];
    for k in 1..n {
        for i in (k..n).rev() {
            f[i] = (f[i] - f[i - 1]) / (xs[i] - xs[i - k]);
        }
        coeffs.push(f[k]);
    }
    Ok(coeffs)
}

/// Evaluate the Newton interpolating polynomial (nested form) at `x`.
pub fn newton_interp(xs: &[f64], ys: &[f64], x: f64) -> Result<f64, String> {
    let coeffs = newton_coefficients(xs, ys)?;
    let mut acc = *coeffs.last().unwrap();
    for i in (0..coeffs.len() - 1).rev() {
        acc = coeffs[i] + (x - xs[i]) * acc;
    }
    Ok(acc)
}

/// Build an interpolant over the sample set, precomputing coefficients for
/// `InterpKind::Newton`.
pub fn build_interpolant(xs: &[f64], ys: &[f64], kind: InterpKind) -> Result<Interpolant, String> {
    validate(xs, ys)?;
    let coefficients = match kind {
        InterpKind::Newton => Some(newton_coefficients(xs, ys)?),
        _ => None,
    };
    Ok(Interpolant {
        xs: xs.to_vec(),
        ys: ys.to_vec(),
        kind,
        coefficients,
    })
}

impl Interpolant {
    /// Evaluate the interpolant at `x` using the stored strategy.
    pub fn evaluate(&self, x: f64) -> Result<f64, String> {
        match self.kind {
            InterpKind::Linear => linear_interp(&self.xs, &self.ys, x),
            InterpKind::Nearest => nearest_interp(&self.xs, &self.ys, x),
            InterpKind::Newton => {
                let coeffs = self
                    .coefficients
                    .as_deref()
                    .ok_or_else(|| "interp: Newton coefficients missing".to_string())?;
                let mut acc = *coeffs.last().unwrap();
                for i in (0..coeffs.len() - 1).rev() {
                    acc = coeffs[i] + (x - self.xs[i]) * acc;
                }
                Ok(acc)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_midpoint() {
        // xs = [0, 2, 4], ys = [0, 4, 8]: y = 2x.
        let xs = [0.0, 2.0, 4.0];
        let ys = [0.0, 4.0, 8.0];
        assert!((linear_interp(&xs, &ys, 1.0).unwrap() - 2.0).abs() < 1e-12);
        assert!((linear_interp(&xs, &ys, 3.0).unwrap() - 6.0).abs() < 1e-12);
        // Out of range errors.
        assert!(linear_interp(&xs, &ys, 5.0).is_err());
        assert!(linear_interp(&xs, &ys, -1.0).is_err());
    }

    #[test]
    fn newton_reproduces_cubic() {
        // y = x³ sampled at 0,1,2,3.
        let xs = [0.0, 1.0, 2.0, 3.0];
        let ys = [0.0, 1.0, 8.0, 27.0];
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            assert!((newton_interp(&xs, &ys, x).unwrap() - y).abs() < 1e-12);
        }
        // Interior point also lands on the cubic.
        assert!((newton_interp(&xs, &ys, 1.5).unwrap() - 3.375).abs() < 1e-12);
    }

    #[test]
    fn nearest_picks_closest() {
        let xs = [0.0, 1.0, 10.0];
        let ys = [0.0, 1.0, 100.0];
        assert_eq!(nearest_interp(&xs, &ys, 0.4).unwrap(), 0.0);
        assert_eq!(nearest_interp(&xs, &ys, 7.0).unwrap(), 100.0);
    }

    #[test]
    fn interpolant_kinds() {
        let xs = [0.0, 1.0, 2.0];
        let ys = [0.0, 1.0, 4.0];
        let lin = build_interpolant(&xs, &ys, InterpKind::Linear).unwrap();
        assert!((lin.evaluate(0.5).unwrap() - 0.5).abs() < 1e-12);
        let newt = build_interpolant(&xs, &ys, InterpKind::Newton).unwrap();
        assert!((newt.evaluate(1.0).unwrap() - 1.0).abs() < 1e-12);
    }
}
