//! Schoenberg interpolation of an N-dimensional spline.
//! Source: `BSplCLib::BuildSchoenbergPoints`, `KnotSequence`, `Interpolate`
//! (`BSplCLib.cxx:478`, `3331`, `3353`; `BSplCLib_2.cxx:327-425`).

use crate::bspl::knots::hunt;

/// `BSplCLib::KnotSequence` non-periodic (`cxx:488-516`).
pub fn knot_sequence(knots: &[f64], mults: &[i32], _degree: i32) -> Vec<f64> {
    let mut seq = Vec::new();
    for (i, &k) in knots.iter().enumerate() {
        let m = *mults.get(i).unwrap_or(&1);
        for _ in 0..m.max(0) {
            seq.push(k);
        }
    }
    seq
}

/// `BSplCLib::BuildSchoenbergPoints` (`cxx:3331-3348`).
pub fn schoenberg_points(degree: usize, flat: &[f64], n_poles: usize) -> Vec<f64> {
    let inv = 1.0 / degree as f64;
    let mut p = vec![0.0; n_poles];
    for i in 0..n_poles {
        let mut s = 0.0;
        for j in 1..=degree {
            s += *flat.get(j + i).unwrap_or(&0.0);
        }
        p[i] = s * inv;
    }
    p
}

/// `BSplCLib::EvalBsplineBasis` for derivative request 0 (`BSplCLib_2.cxx:502-520`).
fn eval_basis0(order: usize, flat: &[f64], parameter: f64) -> Option<(usize, Vec<f64>)> {
    let n_poles = flat.len().saturating_sub(order);
    if n_poles == 0 {
        return None;
    }
    let degree = order - 1;
    let mut span = hunt(flat, parameter).max(degree).min(n_poles.saturating_sub(1) + degree);
    if span < degree {
        span = degree;
    }
    if span > n_poles + degree - 1 {
        span = n_poles + degree - 1;
    }
    // `ii` after subtracting Lower: 0-based knot index of the span.
    let ii = span;
    let first_nz = ii + 1 - order; // 0-based pole index
    let mut basis = vec![0.0; order];
    basis[0] = 1.0;
    const RES: f64 = 1e-30;
    for qq in 2..=order {
        basis[qq - 1] = 0.0;
        for pp in 1..qq {
            let lo = ii + 1 + pp - qq;
            let hi = ii + pp;
            let scale = flat.get(hi).copied().unwrap_or(0.0) - flat.get(lo).copied().unwrap_or(0.0);
            if scale.abs() < RES {
                return None;
            }
            let factor = (parameter - flat.get(lo).copied().unwrap_or(0.0)) / scale;
            let saved = factor * basis[pp - 1];
            basis[pp - 1] *= 1.0 - factor;
            basis[pp - 1] += basis[qq - 1];
            basis[qq - 1] = saved;
        }
    }
    Some((first_nz, basis))
}

/// `BSplCLib::Interpolate` N-D, contact order 0 (`cxx:3353-3398`).
///
/// `poles` is row-major `[n_poles][dimension]`, overwritten with solved poles.
pub fn interpolate(
    degree: usize,
    flat: &[f64],
    parameters: &[f64],
    poles: &mut [f64],
    dimension: usize,
) -> Result<(), i32> {
    let n = parameters.len();
    let bw = 2 * degree + 1;
    let mut matrix = vec![0.0; n * bw];
    for (i, &u) in parameters.iter().enumerate() {
        let (first_nz, basis) = eval_basis0(degree + 1, flat, u).ok_or(2)?;
        // `anIndex = Degree + first_nz_0 - i_0` as 0-based column in the band.
        let col0 = degree as isize + first_nz as isize - i as isize;
        for j in 0..=degree {
            let c = col0 + j as isize;
            if c >= 0 && (c as usize) < bw {
                matrix[i * bw + c as usize] = basis[j];
            }
        }
    }
    let mut pivot = 0i32;
    if factor_banded(&mut matrix, n, degree, &mut pivot) != 0 {
        return Err(pivot);
    }
    if solve_banded(&matrix, n, degree, poles, dimension) != 0 {
        return Err(1);
    }
    Ok(())
}

/// `BSplCLib::FactorBandedMatrix` (`BSplCLib_2.cxx:385-425`), 0-based rows.
fn factor_banded(m: &mut [f64], n: usize, degree: usize, pivot: &mut i32) -> i32 {
    let lower = degree;
    let band = 2 * degree + 1;
    for i in 1..n {
        let min_index = if lower + 1 > i { lower + 1 - i } else { 1 };
        for j in min_index..=lower {
            let an_index = i + 1 - lower + j - 1; // 1-based anIndex, then 0-based row
            let row_idx = an_index - 1;
            let piv = m[row_idx * band + lower];
            if piv.abs() <= f64::MIN_POSITIVE {
                *pivot = an_index as i32;
                return 1;
            }
            let inv = -1.0 / piv;
            m[i * band + (j - 1)] *= inv;
            let max_index = band + an_index - (i + 1);
            for k in (j + 1)..=max_index {
                m[i * band + (k - 1)] +=
                    m[i * band + (j - 1)] * m[row_idx * band + (k + (i + 1) - an_index - 1)];
            }
        }
    }
    0
}

/// `BSplCLib::SolveBandedSystem` (`BSplCLib.cxx:3205-3268`).
fn solve_banded(m: &[f64], n: usize, degree: usize, poles: &mut [f64], dim: usize) -> i32 {
    let lower = degree;
    let upper = degree;
    let band = 2 * degree + 1;
    for ii in 1..n {
        let min_index = if ii >= lower { ii - lower } else { 0 };
        for jj in min_index..ii {
            let coeff = m[ii * band + (jj + lower - ii)];
            for kk in 0..dim {
                poles[ii * dim + kk] += poles[jj * dim + kk] * coeff;
            }
        }
    }
    for ii in (0..n).rev() {
        let max_index = (ii + upper).min(n - 1);
        for jj in ((ii + 1)..=max_index).rev() {
            let coeff = m[ii * band + (jj + lower - ii)];
            for kk in 0..dim {
                poles[ii * dim + kk] -= poles[jj * dim + kk] * coeff;
            }
        }
        let divisor = m[ii * band + lower];
        if divisor.abs() <= 1.0e-16 {
            return 1;
        }
        let inv = 1.0 / divisor;
        for kk in 0..dim {
            poles[ii * dim + kk] *= inv;
        }
    }
    0
}
