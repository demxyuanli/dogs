//! Linear solvers used by `SearchDirection`: `math_Recipes` LU/SVD,
//! `math_Gauss`, `math_GaussLeastSquare` and `math_SVD`.

use super::*;
use crate::SVD;

// ---------------------------------------------------------------------------
// Linear solvers used by SearchDirection
// ---------------------------------------------------------------------------

/// Result of `LU_Decompose` (`math_Recipes.cxx:202-296`).
pub(super) struct LuDecomp {
    pub(super) lu: MathMatrix,
    pub(super) indx: Vec<usize>,
    pub(super) done: bool,
}

/// `LU_Decompose(a, indx, d, TINY)` (`math_Recipes.cxx:202-296`), the implicit
/// row-scaling variant used by `math_Gauss` and `math_GaussLeastSquare`.
pub(super) fn lu_decompose(a: &MathMatrix, tiny: f64) -> LuDecomp {
    let n = a.row_count();
    let mut lu = a.clone();
    let mut indx = vec![0usize; n + 1];
    let mut vv = vec![0.0f64; n + 1];

    for i in 1..=n {
        let mut big = 0.0;
        for j in 1..=n {
            let t = lu.value(i, j).abs();
            if t > big {
                big = t;
            }
        }
        if big <= tiny {
            return LuDecomp {
                lu,
                indx,
                done: false,
            };
        }
        vv[i] = 1.0 / big;
    }

    for j in 1..=n {
        for i in 1..j {
            let mut sum = lu.value(i, j);
            for k in 1..i {
                sum -= lu.value(i, k) * lu.value(k, j);
            }
            lu.set_value(i, j, sum);
        }
        let mut big = 0.0;
        let mut imax = j;
        for i in j..=n {
            let mut sum = lu.value(i, j);
            for k in 1..j {
                sum -= lu.value(i, k) * lu.value(k, j);
            }
            lu.set_value(i, j, sum);
            let dum = vv[i] * sum.abs();
            if dum < big {
                continue;
            }
            big = dum;
            imax = i;
        }
        if j != imax {
            for k in 1..=n {
                let t = lu.value(imax, k);
                lu.set_value(imax, k, lu.value(j, k));
                lu.set_value(j, k, t);
            }
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if lu.value(j, j).abs() <= tiny {
            return LuDecomp {
                lu,
                indx,
                done: false,
            };
        }
        if j != n {
            let dum = 1.0 / lu.value(j, j);
            for i in (j + 1)..=n {
                lu.set_value(i, j, lu.value(i, j) * dum);
            }
        }
    }
    LuDecomp {
        lu,
        indx,
        done: true,
    }
}

/// `LU_Solve(a, indx, b)` (`math_Recipes.cxx:309-344`).
pub(super) fn lu_solve(a: &MathMatrix, indx: &[usize], b: &mut MathVector) {
    let n = a.row_count();
    let nblow = b.lower() - 1;
    let mut ii = 0usize;
    for i in 1..=n {
        let ip = indx[i];
        let mut sum = b.value(ip + nblow);
        let t = b.value(i + nblow);
        b.set_value(ip + nblow, t);
        if ii != 0 {
            for j in ii..i {
                sum -= a.value(i, j) * b.value(j + nblow);
            }
        } else if sum != 0.0 {
            ii = i;
        }
        b.set_value(i + nblow, sum);
    }
    for i in (1..=n).rev() {
        let mut sum = b.value(i + nblow);
        for j in (i + 1)..=n {
            sum -= a.value(i, j) * b.value(j + nblow);
        }
        b.set_value(i + nblow, sum / a.value(i, i));
    }
}

/// `math_SVD` constructor (`math_SVD.cxx:29-39`): pads A with zero rows to
/// `max(RowNumber, ColNumber)` before calling `SVD_Decompose`.
pub(super) fn svd_padded(a: &MathMatrix) -> SVD {
    let m = a.row_count();
    let n = a.col_count();
    let size = m.max(n);
    let mut u = MathMatrix::new(1, size, 1, n);
    for i in 1..=m {
        for j in 1..=n {
            u.set_value(i, j, a.value(i, j));
        }
    }
    SVD::new(&u)
}

/// `math_SVD::Solve` (`math_SVD.cxx:41-68`) and `SVD_Solve`
/// (`math_Recipes.cxx:712-748`), with the default `Eps = 1.0e-6`
/// (`math_SVD.hxx:53`).
pub(super) fn svd_solve_occt(svd: &SVD, b: &MathVector, x: &mut MathVector) {
    let m = svd.u.row_count();
    let n = svd.v.row_count();
    let mut wmax = 0.0f64;
    for i in 1..=svd.w.len() {
        wmax = wmax.max(svd.w.value(i));
    }
    let wmin = 1.0e-6 * wmax;
    let mut w = svd.w.clone();
    for i in 1..=w.len() {
        if w.value(i) < wmin {
            w.set_value(i, 0.0);
        }
    }
    // BB(1..U.RowNumber()) = B padded with zeros.
    let mut bb = MathVector::new(1, m);
    for i in 1..=b.len() {
        bb.set_value(i, b.value(i));
    }
    let mut tmp = MathVector::new(1, n);
    for j in 1..=n {
        let mut s = 0.0;
        if w.value(j) != 0.0 {
            for i in 1..=m {
                s += svd.u.value(i, j) * bb.value(i);
            }
            s /= w.value(j);
        }
        tmp.set_value(j, s);
    }
    for j in 1..=n {
        let mut s = 0.0;
        for jj in 1..=n {
            s += svd.v.value(j, jj) * tmp.value(jj);
        }
        x.set_value(j, s);
    }
}

/// `math_GaussLeastSquare` (`math_GaussLeastSquare.cxx:28-49`): LU of
/// `A^T A` plus the `A^T B` right-hand side.
pub(super) struct GaussLeastSquare {
    pub(super) lu: MathMatrix,
    pub(super) indx: Vec<usize>,
    pub(super) a2: MathMatrix,
    pub(super) done: bool,
}

impl GaussLeastSquare {
    pub(super) fn new(a: &MathMatrix, min_pivot: f64) -> Self {
        let neq = a.row_count();
        let ninc = a.col_count();
        // A2 = A.Transposed() (ninc x neq).
        let mut a2 = MathMatrix::new(1, ninc, 1, neq);
        for i in 1..=neq {
            for j in 1..=ninc {
                a2.set_value(j, i, a.value(i, j));
            }
        }
        let mut lu = MathMatrix::new(1, ninc, 1, ninc);
        lu.multiply_mat(&a2, a); // LU = A2 * A
        let d = lu_decompose(&lu, min_pivot);
        Self {
            lu: d.lu,
            indx: d.indx,
            a2,
            done: d.done,
        }
    }

    /// `Solve(B, X)` (`math_GaussLeastSquare.cxx:40-49`).
    pub(super) fn solve(&self, b: &MathVector, x: &mut MathVector) {
        x.multiply_mat_vec(&self.a2, b); // X = A2 * B
        lu_solve(&self.lu, &self.indx, x);
    }
}
