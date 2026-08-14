//! Curve fitting and interpolation utilities.
use crate::gp::GpPnt;
use crate::bspl::plib;

/// Fit a degree-d polynomial to points (x, y) via least squares.
/// Returns coefficient vector [c0, c1, ..., cd].
pub fn fit_polynomial(xs: &[f64], ys: &[f64], degree: usize) -> Result<Vec<f64>, &'static str> {
    if xs.len() != ys.len() { return Err("fit: length mismatch"); }
    let n = xs.len();
    if n < degree + 1 { return Err("fit: insufficient points"); }

    // Normal equations: sum over samples of x^k for k in 0..=2d
    let m = degree + 1;
    let mut a = vec![vec![0.0f64; m]; m];
    let mut b = vec![0.0f64; m];

    for k in 0..m {
        for j in 0..m {
            let mut sum = 0.0;
            for i in 0..n {
                sum += xs[i].powi((k + j) as i32);
            }
            a[k][j] = sum;
        }
        let mut sum = 0.0;
        for i in 0..n { sum += ys[i] * xs[i].powi(k as i32); }
        b[k] = sum;
    }
    solve_gauss(&a, &b)
}

/// Gaussian elimination with partial pivoting for n×n system. Returns solution.
fn solve_gauss(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, &'static str> {
    let n = a.len();
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    for k in 0..n {
        // Pivot
        let mut max_i = k;
        let mut max_v = a[k][k].abs();
        for i in (k+1)..n {
            if a[i][k].abs() > max_v { max_v = a[i][k].abs(); max_i = i; }
        }
        if max_v < 1e-300 { return Err("fit: singular system"); }
        if max_i != k {
            a.swap(k, max_i);
            b.swap(k, max_i);
        }
        for i in (k+1)..n {
            let f = a[i][k] / a[k][k];
            for j in k..n { a[i][j] -= f * a[k][j]; }
            b[i] -= f * b[k];
        }
    }
    // Back substitution
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i+1)..n { s -= a[i][j] * x[j]; }
        x[i] = s / a[i][i];
    }
    Ok(x)
}

/// Evaluate polynomial fit at x (Horner).
pub fn eval_fit(coeffs: &[f64], x: f64) -> f64 { plib::eval_polynomial(coeffs, x) }

/// Fit a line y = m*x + c to points. Returns (m, c).
pub fn fit_line(xs: &[f64], ys: &[f64]) -> Option<(f64, f64)> {
    let coeffs = fit_polynomial(xs, ys, 1).ok()?;
    Some((coeffs[1], coeffs[0]))
}

/// Total least squares line fit (minimizes perpendicular distance).
pub fn fit_line_orthogonal(pts: &[GpPnt]) -> Option<(GpPnt, crate::gp::GpVec)> {
    let n = pts.len();
    if n < 2 { return None; }
    // Centroid
    let mut cx = 0.0; let mut cy = 0.0;
    for p in pts { cx += p.x(); cy += p.y(); }
    cx /= n as f64; cy /= n as f64;
    // Covariance
    let mut sxx = 0.0; let mut sxy = 0.0; let mut syy = 0.0;
    for p in pts {
        let dx = p.x() - cx; let dy = p.y() - cy;
        sxx += dx*dx; sxy += dx*dy; syy += dy*dy;
    }
    // Principal direction = eigenvector of [[sxx,sxy],[sxy,syy]]
    let theta = 0.5 * (2.0*sxy).atan2(sxx - syy);
    let dir = crate::gp::GpVec::new(theta.cos(), theta.sin(), 0.0);
    Some((GpPnt::new(cx, cy, 0.0), dir))
}

/// Lagrange interpolation through points.
pub fn interpolate_lagrange(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    plib::lagrange_value(xs, ys, x)
}

/// Simple cubic spline interpolation through points (natural boundary).
/// Returns coefficients per segment: for segment i, y = a + b*t + c*t² + d*t³, t∈[0,1].
pub fn cubic_spline_segments(xs: &[f64], ys: &[f64]) -> Result<Vec<[f64; 4]>, &'static str> {
    let n = xs.len();
    if n < 2 || xs.len() != ys.len() { return Err("spline: bad input"); }
    if n == 2 {
        let m = (ys[1] - ys[0]) / (xs[1] - xs[0]);
        return Ok(vec![[ys[0], m * (xs[1]-xs[0]), 0.0, 0.0]]);
    }
    // Solve tridiagonal system for second derivatives (natural splines)
    let mut h = vec![0.0f64; n - 1];
    for i in 0..n-1 { h[i] = xs[i+1] - xs[i]; }
    let mut alpha = vec![0.0f64; n];
    for i in 1..n-1 {
        alpha[i] = 3.0/h[i]*(ys[i+1]-ys[i]) - 3.0/h[i-1]*(ys[i]-ys[i-1]);
    }
    let mut l = vec![0.0f64; n];
    let mut mu = vec![0.0f64; n];
    let mut z = vec![0.0f64; n];
    l[0] = 1.0; mu[0] = 0.0; z[0] = 0.0;
    for i in 1..n-1 {
        l[i] = 2.0*(xs[i+1]-xs[i-1]) - h[i-1]*mu[i-1];
        if l[i].abs() < 1e-300 { return Err("spline: degenerate"); }
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
    }
    l[n-1] = 1.0; z[n-1] = 0.0;
    let mut c = vec![0.0f64; n];
    let mut b = vec![0.0f64; n];
    let mut d = vec![0.0f64; n];
    for j in (0..n-1).rev() {
        c[j] = z[j] - mu[j]*c[j+1];
        b[j] = (ys[j+1]-ys[j])/h[j] - h[j]*(c[j+1]+2.0*c[j])/3.0;
        d[j] = (c[j+1]-c[j])/(3.0*h[j]);
    }
    // Segment coefficients: y(t) = a + b*t + c*t² + d*t³ on [0, h_i]
    let mut segs = Vec::with_capacity(n - 1);
    for i in 0..n-1 {
        segs.push([ys[i], b[i], c[i], d[i]]);
    }
    Ok(segs)
}

/// Evaluate a cubic spline at x given segments + knots.
pub fn eval_spline(xs: &[f64], segs: &[[f64; 4]], x: f64) -> f64 {
    let n = xs.len();
    if n < 2 { return 0.0; }
    // Find segment
    let mut i = 0;
    for k in 0..n-1 { if x >= xs[k] { i = k; } }
    let t = x - xs[i];
    let [a, b, c, d] = segs[i];
    a + b*t + c*t*t + d*t*t*t
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_linear() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![1.0, 3.0, 5.0, 7.0]; // y = 2x + 1
        let coeffs = fit_polynomial(&xs, &ys, 1).unwrap();
        assert!((coeffs[0] - 1.0).abs() < 1e-10);
        assert!((coeffs[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn spline_endpoints() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 0.0];
        let segs = cubic_spline_segments(&xs, &ys).unwrap();
        assert_eq!(segs.len(), 2);
        assert!((eval_spline(&xs, &segs, 0.0)).abs() < 1e-10);
        assert!((eval_spline(&xs, &segs, 1.0) - 1.0).abs() < 1e-10);
        assert!((eval_spline(&xs, &segs, 2.0)).abs() < 1e-10);
    }
}
