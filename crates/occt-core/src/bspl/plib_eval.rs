//! Polynomial evaluation and Hermite endpoint interpolation.
//! Source: `PLib::EvalPolynomial`, `PLib::NoDerivativeEvalPolynomial`,
//! `PLib::HermiteInterpolate` (`PLib.cxx:945-2045`).

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
