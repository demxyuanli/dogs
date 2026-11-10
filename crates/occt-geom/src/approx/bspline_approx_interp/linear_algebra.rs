use occt_math::{MathMatrix, MathVector};

use super::constants::LU_NONZERO_EPS;


/// The `MathLin::LU` result (`MathLin_Gauss.hxx:31-38`), reduced to what
/// `solveLU` (`cxx:592-625`) reads.
pub(super) struct LuResult {
    /// Combined L (below the diagonal, unit diagonal) and U.
    lu: MathMatrix,
    /// `Pivot` (`MathLin_Gauss.hxx:35`), 1-based, `pivot[0]` unused.
    pivot: Vec<usize>,
}

/// `MathLin::LU` (`MathLin_Gauss.hxx:51-155`): scaled implicit pivoting, then
/// Crout elimination. `None` is OCCT's `Status::NumericalError` /
/// `Status::InvalidInput`.
pub(super) fn math_lin_lu(a: &MathMatrix, min_pivot: f64) -> Option<LuResult> {
    let n = a.row_count();
    if n != a.col_count() {
        return None;
    }
    let mut lu = a.clone();
    let mut pivot = vec![0usize; n + 1];

    // Scaling factors for implicit pivoting (`hxx:77-96`).
    let mut scale = vec![0.0f64; n + 1];
    for i in 1..=n {
        let mut max = 0.0f64;
        for j in 1..=n {
            let abs = lu.value(i, j).abs();
            if abs > max {
                max = abs;
            }
        }
        if max < min_pivot {
            return None;
        }
        scale[i] = 1.0 / max;
    }

    // Crout's algorithm with partial pivoting (`hxx:98-143`).
    for k in 1..=n {
        let mut max_scaled = 0.0f64;
        let mut pivot_row = k;
        for i in k..=n {
            let scaled = scale[i] * lu.value(i, k).abs();
            if scaled > max_scaled {
                max_scaled = scaled;
                pivot_row = i;
            }
        }
        if pivot_row != k {
            for j in 1..=n {
                let t = lu.value(pivot_row, j);
                lu.set_value(pivot_row, j, lu.value(k, j));
                lu.set_value(k, j, t);
            }
            scale.swap(pivot_row, k);
        }
        pivot[k] = pivot_row;

        if lu.value(k, k).abs() < min_pivot {
            return None;
        }

        for i in (k + 1)..=n {
            let m = lu.value(i, k) / lu.value(k, k);
            lu.set_value(i, k, m);
            for j in (k + 1)..=n {
                let v = lu.value(i, j) - m * lu.value(k, j);
                lu.set_value(i, j, v);
            }
        }
    }

    Some(LuResult { lu, pivot })
}

/// `solveLU` (`cxx:592-625`).
pub(super) fn solve_lu(res: &LuResult, rhs: &MathVector) -> MathVector {
    let lower = rhs.lower();
    let upper = rhs.upper();
    let mut x = rhs.clone();
    let mut first_non_zero = 0usize;
    for i in lower..=upper {
        let pivot_idx = res.pivot[i];
        let mut sum = x.value(pivot_idx);
        x.set_value(pivot_idx, x.value(i));
        if first_non_zero != 0 {
            for j in first_non_zero..i {
                sum -= res.lu.value(i, j) * x.value(j);
            }
        } else if sum.abs() > LU_NONZERO_EPS {
            first_non_zero = i;
        }
        x.set_value(i, sum);
    }
    for i in (lower..=upper).rev() {
        let mut sum = x.value(i);
        for j in (i + 1)..=upper {
            sum -= res.lu.value(i, j) * x.value(j);
        }
        x.set_value(i, sum / res.lu.value(i, i));
    }
    x
}

/// `math_Matrix::Multiplied(math_Vector)`.
pub(super) fn matrix_times_vector(m: &MathMatrix, v: &MathVector) -> MathVector {
    let rows = m.row_count();
    let cols = m.col_count();
    let mut out = MathVector::new(1, rows);
    for i in 1..=rows {
        let mut sum = 0.0;
        for j in 1..=cols {
            sum += m.value(i, j) * v.value(v.lower() + j - 1);
        }
        out.set_value(i, sum);
    }
    out
}

/// `math_Matrix::Set(R1, R2, C1, C2, theOther)`.
pub(super) fn set_block(
    dst: &mut MathMatrix,
    r1: usize,
    r2: usize,
    c1: usize,
    c2: usize,
    src: &MathMatrix,
) {
    for i in r1..=r2 {
        for j in c1..=c2 {
            dst.set_value(
                i,
                j,
                src.value(i - r1 + src.lower_row(), j - c1 + src.lower_col()),
            );
        }
    }
}

/// `insertKinkKnot` (`cxx:68-119`): raise the multiplicity of a matching knot
/// to `theDegree`, or insert a new knot (which grows both arrays).
pub(super) fn insert_kink_knot(knot: f64, degree: i32, knots: &mut Vec<f64>, mults: &mut Vec<i32>, tol: f64) {
    for i in 0..knots.len() {
        if (knots[i] - knot).abs() < tol {
            mults[i] = mults[i].max(degree);
            return;
        }
    }

    let mut pos = 0;
    while pos < knots.len() && knots[pos] < knot {
        pos += 1;
    }
    knots.insert(pos, knot);
    mults.insert(pos, degree);
}
