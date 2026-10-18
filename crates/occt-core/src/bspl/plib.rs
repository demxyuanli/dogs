//! Polynomial evaluation and interpolation. Source: `PLib.hxx` + `PLib.cxx`
//! Horner scheme, Lagrange interpolation, derivative evaluation.

/// Evaluate polynomial at x using Horner's scheme.
/// coeffs[0] + coeffs[1]*x + coeffs[2]*x² + ... + coeffs[n]*xⁿ
pub fn eval_polynomial_horner(coeffs: &[f64], x: f64) -> f64 {
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
pub fn eval_poly2var_grid(coeffs: &[Vec<f64>], x: f64, y: f64) -> f64 {
    let nx = coeffs.len();
    if nx == 0 { return 0.0; }
    let mut row_vals = Vec::with_capacity(nx);
    for row in coeffs {
        row_vals.push(eval_polynomial_horner(row, y));
    }
    eval_polynomial_horner(&row_vals, x)
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
    eval_polynomial_horner(&coeffs, x)
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
    trimming_w(u1, u2, dim, coefs, None);
}

/// `PLib::Trimming` with `WCoefs` (`PLib.cxx:1660-1715`).
/// `weights` is 1-based `Lower == 1`, stored here 0-based.
pub fn trimming_w(u1: f64, u2: f64, dim: usize, coefs: &mut [f64], mut weights: Option<&mut [f64]>) {
    if dim == 0 || coefs.len() < dim {
        return;
    }
    let lsp = u2 - u1;
    let mut len = coefs.len() / dim;
    if len == 0 {
        return;
    }
    let rat = match &weights {
        Some(w) => w.len() == len,
        None => false,
    };
    len -= 1;
    let upc = coefs.len() - dim;
    let upw = weights.as_ref().map(|w| w.len() as i32).unwrap_or(0);
    for i in 1..=len {
        let mut indc = upc - dim * (i - 1);
        for j in 0..dim {
            let lo = indc - dim + j;
            if lo < coefs.len() && indc + j < coefs.len() {
                coefs[lo] += u1 * coefs[indc + j];
            }
        }
        let mut indw = upw - i as i32 + 1;
        if rat {
            weight_add(&mut weights, indw - 1, u1, indw);
        }
        while indc < upc {
            indc += dim;
            for k in 0..dim {
                let lo = indc - dim + k;
                if lo < coefs.len() && indc + k < coefs.len() {
                    coefs[lo] = u1 * coefs[indc + k] + lsp * coefs[lo];
                }
            }
            if rat {
                indw += 1;
                weight_mix(&mut weights, indw - 1, indw, u1, lsp);
            }
        }
        for j in 0..dim {
            if upc + j < coefs.len() {
                coefs[upc + j] *= lsp;
            }
        }
        if rat {
            if let Some(w) = weights.as_mut() {
                let at = (upw - 1) as usize;
                if at < w.len() {
                    w[at] *= lsp;
                }
            }
        }
    }
}

fn weight_add(weights: &mut Option<&mut [f64]>, dst: i32, scale: f64, src: i32) {
    let Some(w) = weights.as_mut() else {
        return;
    };
    let d = (dst - 1) as usize;
    let s = (src - 1) as usize;
    if d < w.len() && s < w.len() {
        let value = w[s];
        w[d] += scale * value;
    }
}

fn weight_mix(weights: &mut Option<&mut [f64]>, dst: i32, src: i32, u1: f64, lsp: f64) {
    let Some(w) = weights.as_mut() else {
        return;
    };
    let d = (dst - 1) as usize;
    let s = (src - 1) as usize;
    if d < w.len() && s < w.len() {
        w[d] = u1 * w[s] + lsp * w[d];
    }
}

/// `PLib::VTrimming` (`PLib.cxx:1885-1918`), weights from `WCoeffs`.
pub fn v_trimming_w(
    v1: f64,
    v2: f64,
    coeffs: &mut [Vec<crate::gp::GpPnt>],
    weights: &mut [Vec<f64>],
) {
    if coeffs.is_empty() {
        return;
    }
    let cols = coeffs[0].len();
    for (row, wrow) in coeffs.iter_mut().zip(weights.iter_mut()) {
        let mut temp = vec![0.0; cols * 3];
        for (icol, p) in row.iter().enumerate() {
            temp[icol * 3] = p.x();
            temp[icol * 3 + 1] = p.y();
            temp[icol * 3 + 2] = p.z();
        }
        trimming_w(v1, v2, 3, &mut temp, Some(wrow));
        for icol in 0..cols {
            row[icol] = crate::gp::GpPnt::new(
                temp[icol * 3],
                temp[icol * 3 + 1],
                temp[icol * 3 + 2],
            );
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

/// `PLib::UTrimming` (`PLib.cxx:1839-1881`), weights from `WCoeffs`.
pub fn u_trimming_w(
    u1: f64,
    u2: f64,
    coeffs: &mut [Vec<crate::gp::GpPnt>],
    weights: &mut [Vec<f64>],
) {
    if coeffs.is_empty() {
        return;
    }
    let rows = coeffs.len();
    let cols = coeffs[0].len();
    for icol in 0..cols {
        let mut temp = vec![0.0; rows * 3];
        let mut wcol = Vec::with_capacity(rows);
        for (irow, row) in coeffs.iter().enumerate() {
            let p = row[icol];
            temp[irow * 3] = p.x();
            temp[irow * 3 + 1] = p.y();
            temp[irow * 3 + 2] = p.z();
            wcol.push(weights.get(irow).and_then(|r| r.get(icol).copied()).unwrap_or(0.0));
        }
        trimming_w(u1, u2, 3, &mut temp, Some(&mut wcol));
        for irow in 0..rows {
            coeffs[irow][icol] = crate::gp::GpPnt::new(
                temp[irow * 3],
                temp[irow * 3 + 1],
                temp[irow * 3 + 2],
            );
            if let Some(slot) = weights.get_mut(irow).and_then(|r| r.get_mut(icol)) {
                *slot = wcol[irow];
            }
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

/// `PLib::CoefficientsPoles` non-rational `gp_Pnt` arm
/// (`PLib.cxx:1482-1493` and `1522-1589`). `NoWeights()` skips the rational
/// divide. Point `0` and point `n-1` are copied; interior points are divided
/// by `Bin(n-1, i)` and then accumulated from the low end.
pub fn coefficients_poles(coefs: &[crate::gp::GpPnt]) -> Vec<crate::gp::GpPnt> {
    let reflen = coefs.len();
    if reflen == 0 {
        return Vec::new();
    }
    let mut poles = coefs.to_vec();
    for i in 2..reflen {
        let cnp = bin((reflen as i32) - 1, (i as i32) - 1);
        let p = coefs[i - 1];
        poles[i - 1] = crate::gp::GpPnt::new(p.x() / cnp, p.y() / cnp, p.z() / cnp);
    }
    for i in 1..reflen {
        for j in (i..reflen).rev() {
            let prev = poles[j - 1];
            let cur = poles[j];
            poles[j] = crate::gp::GpPnt::new(cur.x() + prev.x(), cur.y() + prev.y(), cur.z() + prev.z());
        }
    }
    poles
}

/// `PLib::CoefficientsPoles` grid with `WCoefs` (`PLib.cxx:1754-1816`).
pub fn coefficients_poles_grid_w(
    coefs: &[Vec<crate::gp::GpPnt>],
    wcoefs: &[Vec<f64>],
) -> (Vec<Vec<crate::gp::GpPnt>>, Vec<Vec<f64>>) {
    let poles = coefficients_poles_grid(coefs);
    let nrows = coefs.len();
    if nrows == 0 {
        return (poles, Vec::new());
    }
    let ncols = coefs[0].len();
    let mut weights = wcoefs.to_vec();
    for u in 0..nrows {
        let Some(row) = weights.get(u).cloned() else {
            continue;
        };
        if row.len() != ncols {
            continue;
        }
        let mut next = row.clone();
        for col in 1..ncols.saturating_sub(1) {
            let cnp = bin((ncols as i32) - 1, col as i32);
            next[col] = row[col] / cnp;
        }
        for i1 in 1..ncols {
            for i2 in (i1..ncols).rev() {
                next[i2] += next[i2 - 1];
            }
        }
        weights[u] = next;
    }
    for v in 0..ncols {
        for row in 1..nrows.saturating_sub(1) {
            let cnp = bin((nrows as i32) - 1, row as i32);
            if let Some(slot) = weights.get_mut(row).and_then(|r| r.get_mut(v)) {
                *slot /= cnp;
            }
        }
        for i1 in 1..nrows {
            for i2 in (i1..nrows).rev() {
                let below = weights.get(i2 - 1).and_then(|r| r.get(v).copied()).unwrap_or(0.0);
                if let Some(slot) = weights.get_mut(i2).and_then(|r| r.get_mut(v)) {
                    *slot += below;
                }
            }
        }
    }
    (poles, weights)
}

/// `PLib::CoefficientsPoles` non-rational `gp_Pnt` grid
/// (`PLib.cxx:1727-1819`). Rows are U and columns are V. `NoWeights2()`
/// skips the rational divide. Each row is converted first, then each column.
pub fn coefficients_poles_grid(coefs: &[Vec<crate::gp::GpPnt>]) -> Vec<Vec<crate::gp::GpPnt>> {
    let nrows = coefs.len();
    if nrows == 0 {
        return Vec::new();
    }
    let ncols = coefs[0].len();
    let mut poles: Vec<Vec<crate::gp::GpPnt>> = coefs.to_vec();
    for u in 0..nrows {
        for col in 1..ncols.saturating_sub(1) {
            let cnp = bin((ncols as i32) - 1, col as i32);
            let p = coefs[u][col];
            poles[u][col] = crate::gp::GpPnt::new(p.x() / cnp, p.y() / cnp, p.z() / cnp);
        }
        for i1 in 1..ncols {
            for i2 in (i1..ncols).rev() {
                let a = poles[u][i2];
                let b = poles[u][i2 - 1];
                poles[u][i2] =
                    crate::gp::GpPnt::new(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
            }
        }
    }
    for v in 0..ncols {
        for row in 1..nrows.saturating_sub(1) {
            let cnp = bin((nrows as i32) - 1, row as i32);
            let p = poles[row][v];
            poles[row][v] = crate::gp::GpPnt::new(p.x() / cnp, p.y() / cnp, p.z() / cnp);
        }
        for i1 in 1..nrows {
            for i2 in (i1..nrows).rev() {
                let a = poles[i2][v];
                let b = poles[i2 - 1][v];
                poles[i2][v] =
                    crate::gp::GpPnt::new(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
            }
        }
    }
    poles
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
    coeffs.iter().map(|c| eval_polynomial_horner(c, x)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn horner_linear() { assert!((eval_polynomial_horner(&[2.0, 3.0], 4.0) - 14.0).abs() < 1e-14); }

    #[test]
    fn horner_quadratic() { assert!((eval_polynomial_horner(&[1.0, 2.0, 3.0], 2.0) - 17.0).abs() < 1e-14); }

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
