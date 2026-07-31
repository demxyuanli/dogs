//! Natural cubic spline interpolation.
//! Source: standard algorithm (Numerical Recipes `spline`/`splint`).

/// Natural cubic spline over strictly increasing data points.
///
/// On segment `i` (between `xs[i]` and `xs[i+1]`) the spline is
/// `a[i] + b[i]·t + c[i]·t² + d[i]·t³` with `t = x - xs[i]`.
#[derive(Debug, Clone)]
pub struct CubicSpline {
    xs: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

/// Build the natural cubic spline (zero second derivative at both ends).
///
/// Errors when fewer than two points are given, lengths differ, or `xs` is
/// not strictly increasing.
pub fn natural_cubic_spline(xs: &[f64], ys: &[f64]) -> Result<CubicSpline, String> {
    let n = xs.len();
    if n != ys.len() || n < 2 {
        return Err("CubicSpline: need at least two points and equal-length xs/ys".to_string());
    }
    for i in 0..n - 1 {
        if xs[i + 1] <= xs[i] {
            return Err("CubicSpline: xs must be strictly increasing".to_string());
        }
    }
    let m = n - 1; // number of segments
    let mut h = vec![0.0; m];
    for i in 0..m {
        h[i] = xs[i + 1] - xs[i];
    }

    // Second derivatives M[0..n-1] with natural boundary M[0] = M[n-1] = 0.
    // Interior equation (i = 1..n-2):
    //   h[i-1] M[i-1] + 2(h[i-1]+h[i]) M[i] + h[i] M[i+1]
    //       = 6( (y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1] )
    let mut mm = vec![0.0; n];
    let n_mid = n - 2;
    if n_mid > 0 {
        let mut sub = vec![0.0; n_mid];
        let mut diag = vec![0.0; n_mid];
        let mut sup = vec![0.0; n_mid];
        let mut rhs = vec![0.0; n_mid];
        for i in 0..n_mid {
            let idx = i + 1;
            sub[i] = h[i];
            diag[i] = 2.0 * (h[i] + h[i + 1]);
            sup[i] = h[i + 1];
            rhs[i] = 6.0 * ((ys[idx + 1] - ys[idx]) / h[idx] - (ys[idx] - ys[idx - 1]) / h[idx - 1]);
        }
        // Thomas algorithm (forward elimination).
        for i in 1..n_mid {
            let w = sub[i] / diag[i - 1];
            diag[i] -= w * sup[i - 1];
            rhs[i] -= w * rhs[i - 1];
        }
        // Back substitution. Unknown M[i] is stored at mm[i], i = 1..=n-2.
        mm[n - 2] = rhs[n_mid - 1] / diag[n_mid - 1];
        for i in (1..n_mid).rev() {
            mm[i] = (rhs[i - 1] - sup[i - 1] * mm[i + 1]) / diag[i - 1];
        }
    }

    let mut a = Vec::with_capacity(m);
    let mut b = Vec::with_capacity(m);
    let mut c = Vec::with_capacity(m);
    let mut d = Vec::with_capacity(m);
    for i in 0..m {
        let mi = mm[i];
        let mi1 = mm[i + 1];
        let hi = h[i];
        a.push(ys[i]);
        b.push((ys[i + 1] - ys[i]) / hi - hi * (2.0 * mi + mi1) / 6.0);
        c.push(mi / 2.0);
        d.push((mi1 - mi) / (6.0 * hi));
    }
    Ok(CubicSpline { xs: xs.to_vec(), a, b, c, d })
}

impl CubicSpline {
    fn segment_index(&self, x: f64) -> usize {
        let m = self.a.len();
        if x <= self.xs[0] {
            return 0;
        }
        for i in 0..m {
            if x < self.xs[i + 1] {
                return i;
            }
        }
        m - 1
    }

    /// Evaluate the spline at `x` (endpoint cubics extrapolate outside the data range).
    pub fn evaluate(&self, x: f64) -> f64 {
        let i = self.segment_index(x);
        let t = x - self.xs[i];
        self.a[i] + t * (self.b[i] + t * (self.c[i] + t * self.d[i]))
    }

    /// First derivative of the spline at `x`.
    pub fn derivative(&self, x: f64) -> f64 {
        let i = self.segment_index(x);
        let t = x - self.xs[i];
        self.b[i] + t * (2.0 * self.c[i] + 3.0 * t * self.d[i])
    }

    /// Definite integral of the spline over `[a, b]`.
    pub fn integrate(&self, a: f64, b: f64) -> f64 {
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        let sign = if a <= b { 1.0 } else { -1.0 };
        sign * (self.antiderivative(hi) - self.antiderivative(lo))
    }

    /// Integral from `xs[0]` to `x` (with linear extrapolation outside the range).
    fn antiderivative(&self, x: f64) -> f64 {
        let m = self.a.len();
        let x0 = self.xs[0];
        let xn = self.xs[m];
        if x <= x0 {
            return seg_int(&self.a[0], &self.b[0], &self.c[0], &self.d[0], 0.0, x - x0);
        }
        if x >= xn {
            let mut total = 0.0;
            for i in 0..m {
                let h = self.xs[i + 1] - self.xs[i];
                total += seg_int(&self.a[i], &self.b[i], &self.c[i], &self.d[i], 0.0, h);
            }
            let i = m - 1;
            let t_base = xn - self.xs[i];
            total += seg_int(&self.a[i], &self.b[i], &self.c[i], &self.d[i], t_base, x - self.xs[i]);
            return total;
        }
        let mut total = 0.0;
        for i in 0..m {
            let h = self.xs[i + 1] - self.xs[i];
            if x < self.xs[i + 1] {
                total += seg_int(&self.a[i], &self.b[i], &self.c[i], &self.d[i], 0.0, x - self.xs[i]);
                break;
            }
            total += seg_int(&self.a[i], &self.b[i], &self.c[i], &self.d[i], 0.0, h);
        }
        total
    }
}

/// ∫_{t0}^{t1} (a + b·t + c·t² + d·t³) dt.
fn seg_int(a: &f64, b: &f64, c: &f64, d: &f64, t0: f64, t1: f64) -> f64 {
    let prim = |t: f64| t * (a + t * (b / 2.0 + t * (c / 3.0 + t * (d / 4.0))));
    prim(t1) - prim(t0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolates_input_points() {
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = [0.0, 0.8, 1.1, 1.4, 1.2];
        let s = natural_cubic_spline(&xs, &ys).unwrap();
        for i in 0..xs.len() {
            assert!((s.evaluate(xs[i]) - ys[i]).abs() < 1e-12, "knot {i}");
        }
    }

    #[test]
    fn linear_data_is_exact() {
        let xs = [0.0, 1.0, 2.0, 3.0];
        let ys = [1.0, 2.0, 3.0, 4.0];
        let s = natural_cubic_spline(&xs, &ys).unwrap();
        assert!((s.evaluate(1.5) - 2.5).abs() < 1e-12);
        assert!((s.derivative(1.5) - 1.0).abs() < 1e-12);
        assert!((s.integrate(0.0, 3.0) - 7.5).abs() < 1e-12); // ∫(x+1) = 7.5
    }

    #[test]
    fn smooth_at_interior_knots() {
        let xs: [f64; 6] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let ys: Vec<f64> = xs.iter().map(|&x| x.sin()).collect();
        let s = natural_cubic_spline(&xs, &ys).unwrap();
        for i in 1..xs.len() - 1 {
            let left = s.derivative(xs[i] - 1e-9);
            let right = s.derivative(xs[i] + 1e-9);
            assert!((left - right).abs() < 1e-6, "C1 discontinuity at knot {i}: {left} vs {right}");
        }
    }

    #[test]
    fn validates_input() {
        assert!(natural_cubic_spline(&[0.0, 1.0], &[0.0]).is_err());
        assert!(natural_cubic_spline(&[0.0, 1.0, 1.0], &[0.0, 1.0, 2.0]).is_err());
    }
}
