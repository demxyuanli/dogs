//! Multiple integration over a box via tensor-product Gauss-Legendre.
//! Source: `math_GaussMultipleIntegration`.

use crate::{gauss_legendre, MathStatus, MathVector};

fn check_box(lo: &[f64], hi: &[f64]) -> Result<usize, MathStatus> {
    if lo.is_empty() || lo.len() != hi.len() {
        return Err(MathStatus::FunctionError);
    }
    Ok(lo.len())
}

/// Recursive tensor-product integration. `axis` is the current 0-based dimension;
/// `point` accumulates the evaluation point across dimensions.
fn rec_integrate<F: Fn(&MathVector) -> f64>(
    f: &F,
    lo: &[f64],
    hi: &[f64],
    n_pts: usize,
    axis: usize,
    point: &mut MathVector,
) -> f64 {
    if axis == lo.len() {
        return f(point);
    }
    let (nodes, weights) = gauss_legendre(lo[axis], hi[axis], n_pts);
    let mut sum = 0.0;
    for (x, w) in nodes.iter().zip(weights.iter()) {
        point.set_value(axis + 1, *x);
        sum += w * rec_integrate(f, lo, hi, n_pts, axis + 1, point);
    }
    sum
}

/// Integrate `f` over the n-dimensional box `[lo, hi]` with `n_pts` Gauss-Legendre
/// points per axis. Returns 0.0 for an invalid or empty box.
pub fn integrate_multidim<F: Fn(&MathVector) -> f64>(
    f: &F,
    lo: &[f64],
    hi: &[f64],
    n_pts: usize,
) -> f64 {
    let dim = match check_box(lo, hi) {
        Ok(d) => d,
        Err(_) => return 0.0,
    };
    if n_pts == 0 {
        return 0.0;
    }
    let mut point = MathVector::new(1, dim);
    rec_integrate(f, lo, hi, n_pts, 0, &mut point)
}

/// Double integral of `f(x, y)` over `[x0,x1] x [y0,y1]` via nested Gauss-Legendre.
pub fn integrate_2d<F: Fn(f64, f64) -> f64>(
    f: &F,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    n: usize,
) -> f64 {
    let (xn, wx) = gauss_legendre(x0, x1, n);
    let (yn, wy) = gauss_legendre(y0, y1, n);
    let mut sum = 0.0;
    for i in 0..n {
        for j in 0..n {
            sum += wx[i] * wy[j] * f(xn[i], yn[j]);
        }
    }
    sum
}

/// Triple integral of `f(x, y, z)` over the box.
pub fn integrate_3d<F: Fn(f64, f64, f64) -> f64>(
    f: &F,
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
    z0: f64,
    z1: f64,
    n: usize,
) -> f64 {
    let (xn, wx) = gauss_legendre(x0, x1, n);
    let (yn, wy) = gauss_legendre(y0, y1, n);
    let (zn, wz) = gauss_legendre(z0, z1, n);
    let mut sum = 0.0;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                sum += wx[i] * wy[j] * wz[k] * f(xn[i], yn[j], zn[k]);
            }
        }
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integral_2d_xy() {
        let f = |x: f64, y: f64| x * y;
        let v = integrate_2d(&f, 0.0, 1.0, 0.0, 1.0, 4);
        assert!((v - 0.25).abs() < 1e-12, "got {v}");
    }

    #[test]
    fn integral_3d_x() {
        let f = |x: f64, y: f64, z: f64| x;
        let v = integrate_3d(&f, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 4);
        assert!((v - 0.5).abs() < 1e-12, "got {v}");
    }

    #[test]
    fn integral_multidim_2d() {
        let f = |p: &MathVector| p.value(1) * p.value(2);
        let v = integrate_multidim(&f, &[0.0, 0.0], &[1.0, 1.0], 4);
        assert!((v - 0.25).abs() < 1e-12, "got {v}");
    }
}
