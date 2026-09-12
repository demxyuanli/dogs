//! Polynomial evaluation and interpolation. Source: `PLib.hxx` + `PLib.cxx`
//! Horner scheme, Lagrange interpolation, derivative evaluation.

/// Evaluate polynomial at x using Horner's scheme.
/// coeffs[0] + coeffs[1]*x + coeffs[2]*x² + ... + coeffs[n]*xⁿ
pub fn eval_polynomial(coeffs: &[f64], x: f64) -> f64 {
    if coeffs.is_empty() { return 0.0; }
    let mut result = coeffs[coeffs.len() - 1];
    for i in (0..coeffs.len() - 1).rev() {
        result = result * x + coeffs[i];
    }
    result
}

/// Evaluate polynomial and its first n_deriv derivatives.
/// Returns vector [f(x), f'(x), f''(x), ...].
pub fn eval_polynomial_derivs(coeffs: &[f64], x: f64, n_deriv: usize) -> Vec<f64> {
    let n = coeffs.len();
    let mut derivs = vec![0.0f64; n_deriv + 1];

    // Horner with derivative accumulation
    for k in 0..=n_deriv {
        let mut val = 0.0;
        let mut fact = 1.0f64;
        for i in k..n {
            val = val * x + coeffs[n - 1 - i + k] * fact;
            fact *= (i + 1) as f64;
        }
        derivs[k] = val;
    }
    derivs
}

/// Evaluate 2-variable polynomial: sum(c[i][j] * x^i * y^j)
pub fn eval_poly2var(coeffs: &[Vec<f64>], x: f64, y: f64) -> f64 {
    let nx = coeffs.len();
    if nx == 0 { return 0.0; }
    let mut row_vals = Vec::with_capacity(nx);
    for row in coeffs {
        row_vals.push(eval_polynomial(row, y));
    }
    eval_polynomial(&row_vals, x)
}

/// Lagrange interpolation: find polynomial passing through (x[i], y[i]).
/// Returns coefficients [c0, c1, ..., cn] where P(x) = c0 + c1*x + c2*x² + ...
pub fn lagrange_coefficients(x_vals: &[f64], y_vals: &[f64]) -> Vec<f64> {
    let n = x_vals.len();
    if n == 0 { return vec![]; }
    let mut coeffs = vec![0.0f64; n];
    let mut denom = vec![1.0f64; n];

    // Compute denominators
    for i in 0..n {
        for j in 0..n {
            if i != j { denom[i] *= x_vals[i] - x_vals[j]; }
        }
    }

    // Build Lagrange basis via product of (x - x_j)
    for i in 0..n {
        let mut poly = vec![1.0f64]; // Start with 1
        for j in 0..n {
            if i != j {
                // Multiply poly by (x - x_j)
                let mut new_poly = vec![0.0f64; poly.len() + 1];
                for k in 0..poly.len() {
                    new_poly[k] += poly[k] * (-x_vals[j]);
                    new_poly[k + 1] += poly[k];
                }
                poly = new_poly;
            }
        }
        let factor = y_vals[i] / denom[i];
        for k in 0..poly.len() {
            if k < coeffs.len() { coeffs[k] += poly[k] * factor; }
        }
    }
    coeffs
}

/// Compute the value of a Lagrange interpolation at given x.
pub fn lagrange_value(x_vals: &[f64], y_vals: &[f64], x: f64) -> f64 {
    let coeffs = lagrange_coefficients(x_vals, y_vals);
    eval_polynomial(&coeffs, x)
}

/// Convert power basis to matrix form for Bezier representation.
pub fn power_to_bezier(degree: usize) -> Vec<f64> {
    // Binomial coefficients for conversion
    let n = degree + 1;
    let mut matrix = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let idx = i * n + j;
            let sign = if (i - j) % 2 == 0 { 1.0 } else { -1.0 };
            let binom = binomial(degree, i) as f64 * binomial(i, j) as f64;
            matrix[idx] = sign * binom / binomial(degree, j) as f64;
        }
    }
    matrix
}

/// `PLib::Trimming` dimension arm (`PLib.cxx:1642-1716`).
/// `coefs` is degree-major, `dim`-minor, length `(degree+1)*dim`.
pub fn trimming(u1: f64, u2: f64, dim: usize, coefs: &mut [f64]) {
    if dim == 0 || coefs.len() < dim {
        return;
    }
    let lsp = u2 - u1;
    let mut len = coefs.len() / dim;
    if len == 0 {
        return;
    }
    len -= 1;
    let upc = coefs.len() - dim;
    for _i in 1..=len {
        let mut indc = upc - dim * (_i - 1);
        for j in 0..dim {
            let lo = indc - dim + j;
            if lo < coefs.len() && indc + j < coefs.len() {
                coefs[lo] += u1 * coefs[indc + j];
            }
        }
        while indc < upc {
            indc += dim;
            for k in 0..dim {
                let lo = indc - dim + k;
                if lo < coefs.len() && indc + k < coefs.len() {
                    coefs[lo] = u1 * coefs[indc + k] + lsp * coefs[lo];
                }
            }
        }
        for j in 0..dim {
            if upc + j < coefs.len() {
                coefs[upc + j] *= lsp;
            }
        }
    }
}

/// `PLib::VTrimming` non-rational (`PLib.cxx:1885-1918`).
pub fn v_trimming(v1: f64, v2: f64, coeffs: &mut [Vec<crate::gp::GpPnt>]) {
    if coeffs.is_empty() {
        return;
    }
    let cols = coeffs[0].len();
    for row in coeffs.iter_mut() {
        let mut temp = vec![0.0; cols * 3];
        for (icol, p) in row.iter().enumerate() {
            temp[icol * 3] = p.x();
            temp[icol * 3 + 1] = p.y();
            temp[icol * 3 + 2] = p.z();
        }
        trimming(v1, v2, 3, &mut temp);
        for icol in 0..cols {
            row[icol] = crate::gp::GpPnt::new(
                temp[icol * 3],
                temp[icol * 3 + 1],
                temp[icol * 3 + 2],
            );
        }
    }
}

/// `PLib::UTrimming` non-rational (`PLib.cxx:1839-1881`).
pub fn u_trimming(u1: f64, u2: f64, coeffs: &mut [Vec<crate::gp::GpPnt>]) {
    if coeffs.is_empty() {
        return;
    }
    let rows = coeffs.len();
    let cols = coeffs[0].len();
    for icol in 0..cols {
        let mut temp = vec![0.0; rows * 3];
        for (irow, row) in coeffs.iter().enumerate() {
            let p = row[icol];
            temp[irow * 3] = p.x();
            temp[irow * 3 + 1] = p.y();
            temp[irow * 3 + 2] = p.z();
        }
        trimming(u1, u2, 3, &mut temp);
        for irow in 0..rows {
            coeffs[irow][icol] = crate::gp::GpPnt::new(
                temp[irow * 3],
                temp[irow * 3 + 1],
                temp[irow * 3 + 2],
            );
        }
    }
}

/// `PLib::Bin` (`PLib.cxx:260-270`). N in `[0, 25]`, P in `[0, N]`.
pub fn bin(n: i32, p: i32) -> f64 {
    if n < 0 || n > 25 || p < 0 || p > n {
        return 0.0;
    }
    binomial(n as usize, p as usize) as f64
}

/// Binomial coefficient C(n, k).
fn binomial(n: usize, k: usize) -> usize {
    if k > n { return 0; }
    let k = k.min(n - k);
    let mut c = 1usize;
    for i in 0..k { c = c * (n - i) / (i + 1); }
    c
}

/// Evaluate multiple polynomials (array of coefficient sets) at x.
pub fn eval_polynomials(coeffs: &[Vec<f64>], x: f64) -> Vec<f64> {
    coeffs.iter().map(|c| eval_polynomial(c, x)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn horner_linear() { assert!((eval_polynomial(&[2.0, 3.0], 4.0) - 14.0).abs() < 1e-14); }

    #[test]
    fn horner_quadratic() { assert!((eval_polynomial(&[1.0, 2.0, 3.0], 2.0) - 17.0).abs() < 1e-14); }

    #[test]
    fn lagrange_linear() {
        let x = vec![0.0, 1.0];
        let y = vec![2.0, 5.0]; // y = 3x + 2
        let v = lagrange_value(&x, &y, 0.5);
        assert!((v - 3.5).abs() < 1e-14);
    }

    #[test]
    fn binomial_test() { assert_eq!(binomial(5, 2), 10); assert_eq!(binomial(4, 2), 6); }
}
