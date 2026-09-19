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

/// `BSplCLib::EvalBsplineBasis` (`BSplCLib_2.cxx:429-563`).
///
/// Returns `(first_nz_0based, basis)` where `basis[der * order + j]` is
/// derivative `der` of basis function `j`.
pub fn eval_bspline_basis(
    derivative_request: usize,
    order: usize,
    flat: &[f64],
    parameter: f64,
) -> Result<(usize, Vec<f64>), i32> {
    if order == 0 || flat.len() < order {
        return Err(1);
    }
    let n_poles = flat.len() - order;
    let local_req = derivative_request.min(order.saturating_sub(1));
    // `BSplCLib::EvalBsplineBasis` (`BSplCLib_2.cxx:477-484`):
    // `LocateParameter(Order - 1, FlatKnots, Parameter, isPeriodic = false,
    // Order, aNumPoles + 1, ii, aNewParam)`, then
    // `FirstNonZeroBsplineIndex = ii - Order + 1`. The clamp range is
    // `[Order, aNumPoles]` (1-based), not `[degree, numPoles + degree - 1]`.
    let (ii1, _new_u) = crate::bspl::locate::locate_parameter_range(
        flat,
        parameter,
        false,
        order as i32,
        n_poles as i32 + 1,
        0.0,
        1.0,
    );
    let first_nz = ii1 as usize - order;
    let ii = (ii1 - 1) as isize;
    let ncols = order;
    let nrows = local_req + 1;
    let mut basis = vec![0.0; nrows * ncols];
    basis[0] = 1.0;
    const RES: f64 = 1e-30;
    // Knot indices use signed arithmetic like OCCT `BSplCLib_2.cxx:508-514`
    // (`ii - qq + pp + 1`); usize subtraction underflows for small `ii`.
    let knot = |idx: isize| -> f64 {
        if idx < 0 {
            0.0
        } else {
            flat.get(idx as usize).copied().unwrap_or(0.0)
        }
    };
    for qq in 2..=(order - local_req) {
        basis[qq - 1] = 0.0;
        for pp in 1..qq {
            let lo = ii as isize - qq as isize + pp as isize + 1;
            let hi = ii as isize + pp as isize;
            let scale = knot(hi) - knot(lo);
            if scale.abs() < RES {
                return Err(2);
            }
            let factor = (parameter - knot(lo)) / scale;
            let saved = factor * basis[pp - 1];
            basis[pp - 1] *= 1.0 - factor;
            basis[pp - 1] += basis[qq - 1];
            basis[qq - 1] = saved;
        }
    }
    for qq in (order - local_req + 1)..=order {
        for pp in 1..qq {
            let row = order - qq + 1;
            if row < nrows {
                basis[row * ncols + (pp - 1)] = basis[pp - 1];
            }
        }
        basis[qq - 1] = 0.0;
        for ss in (order - local_req + 1)..=qq {
            let row = order - ss + 1;
            if row < nrows {
                basis[row * ncols + (qq - 1)] = 0.0;
            }
        }
        for pp in 1..qq {
            let lo = ii as isize - qq as isize + pp as isize + 1;
            let hi = ii as isize + pp as isize;
            let scale = knot(hi) - knot(lo);
            if scale.abs() < RES {
                return Err(2);
            }
            let inv = 1.0 / scale;
            let factor = (parameter - knot(lo)) * inv;
            let mut saved = factor * basis[pp - 1];
            basis[pp - 1] *= 1.0 - factor;
            basis[pp - 1] += basis[qq - 1];
            basis[qq - 1] = saved;
            let local_inv = (qq - 1) as f64 * inv;
            for ss in (order - local_req + 1)..=qq {
                let row = order - ss + 1;
                if row >= nrows {
                    continue;
                }
                let base = row * ncols;
                saved = local_inv * basis[base + (pp - 1)];
                basis[base + (pp - 1)] *= -local_inv;
                basis[base + (pp - 1)] += basis[base + (qq - 1)];
                basis[base + (qq - 1)] = saved;
            }
        }
    }
    let _ = local_req;
    Ok((first_nz, basis))
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
    let contact = vec![0i32; parameters.len()];
    interpolate_contact(degree, flat, parameters, &contact, poles, dimension)
}

/// `BSplCLib::Interpolate` with `ContactOrderArray` (`cxx:3353-3398`).
pub fn interpolate_contact(
    degree: usize,
    flat: &[f64],
    parameters: &[f64],
    contact: &[i32],
    poles: &mut [f64],
    dimension: usize,
) -> Result<(), i32> {
    let n = parameters.len();
    if contact.len() != n {
        return Err(1);
    }
    let order = degree + 1;
    let bw = 2 * degree + 1;
    let mut matrix = vec![0.0; n * bw];
    for (i, &u) in parameters.iter().enumerate() {
        let der = contact[i].max(0) as usize;
        let (first_nz, basis) = eval_bspline_basis(der, order, flat, u)?;
        let col0 = degree as isize + first_nz as isize - i as isize;
        let row = der * order;
        for j in 0..=degree {
            let c = col0 + j as isize;
            if c >= 0 && (c as usize) < bw {
                matrix[i * bw + c as usize] = basis[row + j];
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
    let lower = degree as isize;
    let band = 2 * degree + 1;
    // OCCT uses `int` for `i` / `anIndex` (`cxx:397-415`). Keep signed
    // arithmetic so `i + 1 - lower` does not underflow when `lower > i`.
    for i in 1..n {
        let i = i as isize;
        let min_index = if lower + 1 > i { lower + 1 - i } else { 1 };
        for j in min_index..=lower {
            // `cxx:404`: `anIndex = i - LowerBandWidth + j - 1` with OCCT's
            // 1-based `i`; `i` here is 0-based, so the `- 1` must not be added.
            let an_index = i - lower + j;
            if an_index < 1 {
                *pivot = an_index as i32;
                return 1;
            }
            let row_idx = (an_index - 1) as usize;
            let piv = m[row_idx * band + lower as usize];
            if piv.abs() <= f64::MIN_POSITIVE {
                *pivot = an_index as i32;
                return 1;
            }
            let inv = -1.0 / piv;
            m[i as usize * band + (j as usize - 1)] *= inv;
            let max_index = band as isize + an_index - (i + 1);
            for k in (j + 1)..=max_index {
                let k_u = k as usize;
                let j_u = j as usize;
                let col = (k + (i + 1) - an_index - 1) as usize;
                m[i as usize * band + (k_u - 1)] +=
                    m[i as usize * band + (j_u - 1)] * m[row_idx * band + col];
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
