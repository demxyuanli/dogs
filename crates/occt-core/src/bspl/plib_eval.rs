//! Polynomial evaluation and Hermite endpoint interpolation.
//! Source: `PLib::EvalPolynomial`, `PLib::NoDerivativeEvalPolynomial`,
//! `PLib::EvalLagrange`, `PLib::HermiteInterpolate`
//! (`PLib.cxx:945-2045`).

/// `PLib::EvalLagrange` (`PLib.cxx:1122-1250`) for 2D points.
///
/// `values` / `params` hold `degree + 1` points and their assigned parameters.
/// Returns `min(derivative_request, degree) + 1` entries: result `k` is the
/// `k`-th derivative of the Lagrange polynomial at `parameter`. `Err(1)` is
/// OCCT's `ReturnCode` for a repeated parameter (`cxx:1189`).
pub fn eval_lagrange(
    parameter: f64,
    derivative_request: usize,
    degree: usize,
    values: &[[f64; 2]],
    params: &[f64],
) -> Result<Vec<[f64; 2]>, i32> {
    const DIM: usize = 2;
    if degree == 0 || values.len() < degree + 1 || params.len() < degree + 1 {
        return Err(1);
    }
    let local_request = derivative_request.min(degree);
    // `cxx:1170-1200`: in-place divided differences, one row per difference
    // order (`divided_differences_array[jj * Dimension + kk]`).
    let mut dd: Vec<[f64; 2]> = values[..=degree].to_vec();
    for ii in (0..=degree).rev() {
        for jj in ((degree - ii + 1)..=degree).rev() {
            for kk in 0..DIM {
                dd[jj][kk] -= dd[jj - 1][kk];
            }
            let difference = params[jj] - params[jj + ii - degree - 1];
            if difference.abs() < f64::MIN_POSITIVE {
                return Err(1);
            }
            let inv = 1.0 / difference;
            for kk in 0..DIM {
                dd[jj][kk] *= inv;
            }
        }
    }
    // `cxx:1211-1247`: Horner with the divided differences, derivative orders
    // updated from the highest down so order `jj - 1` is still the previous one.
    let mut result = vec![[0.0f64; 2]; local_request + 1];
    result[0] = dd[degree];
    for ii in (1..=degree).rev() {
        let difference = parameter - params[ii - 1];
        for jj in (1..=local_request).rev() {
            let prev = result[jj - 1];
            for kk in 0..DIM {
                result[jj][kk] = result[jj][kk] * difference + prev[kk] * jj as f64;
            }
        }
        for kk in 0..DIM {
            result[0][kk] = result[0][kk] * difference + dd[ii - 1][kk];
        }
    }
    Ok(result)
}

/// Horner evaluation of a vector-valued polynomial (derivative order 0).
/// Coefficients are stored degree-major, dimension-minor:
/// `c0[0..dim), c1[0..dim), ..., cDegree[0..dim)`.
/// `PLib::EvalPoly2Var` (`PLib.cxx:1055-1114`) with derivative requests 0,0.
pub fn eval_poly2var(
    u: f64,
    v: f64,
    u_degree: i32,
    v_degree: i32,
    dimension: usize,
    coeffs: &[f64],
    out: &mut [f64],
) {
    let udim = ((v_degree + 1) as usize) * dimension;
    let mut curve = vec![0.0; udim];
    eval_poly0(coeffs, u_degree, udim, u, &mut curve);
    eval_poly0(&curve, v_degree, dimension, v, out);
}

pub fn eval_poly0(coeffs: &[f64], degree: i32, dimension: usize, u: f64, out: &mut [f64]) {
    for d in 0..dimension {
        out[d] = 0.0;
    }
    for k in (0..=degree).rev() {
        let base = (k as usize) * dimension;
        for d in 0..dimension {
            out[d] = out[d] * u + coeffs.get(base + d).copied().unwrap_or(0.0);
        }
    }
}

/// `PLib::EvalPolynomial` (`PLib.cxx:945-1026`) for `derivative_request >= 0`.
///
/// Coefficients are stored degree-major, dimension-minor (`PLib.cxx:951-966`).
/// `out` must hold `(1 + derivative_request) * dimension` values; the block
/// `[k * dimension, (k + 1) * dimension)` receives derivative order `k`.
/// The `cxx` dispatches orders 1/2 to optimized helpers and everything else
/// (including 0) to the general loop; this is the general loop.
pub fn eval_polynomial(
    par: f64,
    derivative_request: i32,
    degree: i32,
    dimension: usize,
    coeffs: &[f64],
    out: &mut [f64],
) {
    if dimension == 0 {
        return;
    }
    let deriv = derivative_request.max(0) as usize;
    let degree = degree.max(0) as usize;
    let res_size = (1 + deriv) * dimension;
    for v in out[..res_size].iter_mut() {
        *v = 0.0;
    }
    for deg in 0..=degree {
        let base = (degree - deg) * dimension;
        let mut ptr = deriv * dimension;
        for d in (1..=deriv).rev() {
            let orig = ptr - dimension;
            for i in 0..dimension {
                out[ptr + i] = out[ptr + i] * par + out[orig + i] * d as f64;
            }
            ptr = orig;
        }
        for i in 0..dimension {
            out[ptr + i] = out[ptr + i] * par + coeffs.get(base + i).copied().unwrap_or(0.0);
        }
    }
}

/// `PLib::HermiteInterpolate` (`cxx:1931-2045`).
///
/// `first_constr` / `last_constr` are `[dimension][order+1]` with order 0 = value.
/// Output coefficients use the same degree-major layout as [`eval_poly0`].
pub fn hermite_interpolate(
    dimension: usize,
    first_parameter: f64,
    last_parameter: f64,
    first_order: i32,
    last_order: i32,
    first_constr: &[Vec<f64>],
    last_constr: &[Vec<f64>],
    coeffs: &mut [f64],
) -> bool {
    let n = (first_order + last_order + 2) as usize;
    let mut a = vec![vec![0.0; n]; n];
    let pattern = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        [0.0, 0.0, 2.0, 6.0, 12.0, 20.0],
    ];
    for irow in 0..=first_order as usize {
        let mut first_val = 1.0;
        for icol in 0..n {
            a[irow][icol] = pattern[irow][icol] * first_val;
            if irow <= icol {
                first_val *= first_parameter;
            }
        }
    }
    for irow in 0..=last_order as usize {
        let mut last_val = 1.0;
        for icol in 0..n {
            a[irow + first_order as usize + 1][icol] = pattern[irow][icol] * last_val;
            if irow <= icol {
                last_val *= last_parameter;
            }
        }
    }
    let Some(lu) = gauss_factor(&a) else {
        return false;
    };
    for idim in 0..dimension {
        let mut b = vec![0.0; n];
        for icol in 0..=first_order as usize {
            b[icol] = first_constr[idim][icol];
        }
        for icol in 0..=last_order as usize {
            b[first_order as usize + 1 + icol] = last_constr[idim][icol];
        }
        gauss_solve(&lu, &mut b);
        for icol in 0..n {
            coeffs[dimension * icol + idim] = b[icol];
        }
    }
    true
}

struct Lu {
    a: Vec<Vec<f64>>,
    piv: Vec<usize>,
}

fn gauss_factor(a: &[Vec<f64>]) -> Option<Lu> {
    let n = a.len();
    let mut m = a.to_vec();
    let mut piv: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let mut best = k;
        let mut best_v = m[k][k].abs();
        for i in (k + 1)..n {
            let v = m[i][k].abs();
            if v > best_v {
                best_v = v;
                best = i;
            }
        }
        if best_v <= 0.0 {
            return None;
        }
        if best != k {
            m.swap(k, best);
            piv.swap(k, best);
        }
        for i in (k + 1)..n {
            let f = m[i][k] / m[k][k];
            m[i][k] = f;
            for j in (k + 1)..n {
                m[i][j] -= f * m[k][j];
            }
        }
    }
    Some(Lu { a: m, piv })
}

fn gauss_solve(lu: &Lu, b: &mut [f64]) {
    let n = lu.a.len();
    let mut x = vec![0.0; n];
    for i in 0..n {
        x[i] = b[lu.piv[i]];
    }
    for i in 0..n {
        for j in 0..i {
            x[i] -= lu.a[i][j] * x[j];
        }
    }
    for i in (0..n).rev() {
        for j in (i + 1)..n {
            x[i] -= lu.a[i][j] * x[j];
        }
        x[i] /= lu.a[i][i];
    }
    b.copy_from_slice(&x);
}
