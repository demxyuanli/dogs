//! `math_Gauss` — Gauss LU decomposition (Crout algorithm) with partial
//! pivoting, plus the underlying `math_Recipes` LU kernels.
//!
//! Source: `math_Gauss.cxx:23-58` (constructor / `Solve` / `Determinant`),
//! `math_Gauss.hxx:44-80` (`MinPivot` = `1.0e-20` default, `IsDone`,
//! `Solve` / `Determinant` contracts), `math_Recipes.cxx:202-344`
//! (`LU_Decompose` with the implicit row-scale vector, `LU_Solve`).
//!
//! Translation notes:
//! - OCCT's `LU` is a `math_Matrix(1, n, 1, n)` assigned from `A`; the port
//!   copies `A` into a 1-based [`MathMatrix`] the same way, so the 1-based
//!   indexing of `LU_Decompose` / `LU_Solve` is reproduced verbatim.
//! - `math_NotSquare_Raise_if` becomes the `assert!` in [`MathGauss::with_min_pivot`].
//! - `StdFail_NotDone_Raise_if` in `Solve` / `Determinant` becomes the `Err`
//!   return of [`MathGauss::solve`] / [`MathGauss::determinant`].
//! - `LU_Decompose` returns `math_Status_SingularMatrix` when a row's scale or
//!   a pivot falls below `TINY`; the port folds that into `done = false`.
//!   `math_Status_UserAborted` (progress break, `cxx:289-292`) has no
//!   counterpart because the port carries no progress range.

use crate::{MathMatrix, MathVector};

/// `math_Gauss::MinPivot` default (`math_Gauss.hxx:49`).
pub const MIN_PIVOT: f64 = 1.0e-20;

/// `math_Gauss` (`math_Gauss.cxx:23-58`): the LU factors plus the pivot row
/// table and the determinantal sign.
#[derive(Debug, Clone)]
pub struct MathGauss {
    /// `LU` — the in-place `L`/`U` factors (`math_Recipes.cxx:228-286`).
    lu: MathMatrix,
    /// `Index` — `indx(i)` is the row swapped with row `i` (`cxx:277`).
    index: Vec<usize>,
    /// `D` — `+1`/`-1`, flipped on each row interchange (`cxx:212`, `:273`).
    d: f64,
    /// `Done` — `Error == 0` (`math_Gauss.cxx:35-36`).
    done: bool,
}

impl MathGauss {
    /// `math_Gauss(A)` (`math_Gauss.cxx:23-37`) with the default `MinPivot`.
    pub fn new(a: &MathMatrix) -> Self {
        Self::with_min_pivot(a, MIN_PIVOT)
    }

    /// `math_Gauss(A, MinPivot)` (`math_Gauss.cxx:23-37`).
    ///
    /// ```text
    /// LU(1, A.RowNumber(), 1, A.ColNumber());  Index(1, A.RowNumber());
    /// D = 0.0;  Done = false;
    /// math_NotSquare_Raise_if(A.RowNumber() != A.ColNumber(), " ");
    /// LU = A;
    /// int Error = LU_Decompose(LU, Index, D, MinPivot, theProgress);
    /// Done = Error == 0;
    /// ```
    pub fn with_min_pivot(a: &MathMatrix, min_pivot: f64) -> Self {
        assert!(a.row_count() == a.col_count(), "math_Gauss: matrix not square");
        let n = a.row_count();
        let mut lu = MathMatrix::new(1, n, 1, n);
        // `LU = A` is `math_Matrix::Initialized` (`math_Matrix.hxx:363`): a
        // positional copy, so `LU(1, 1)` is `A(A.LowerRow(), A.LowerCol())`
        // whatever the bounds of `A` are.
        for i in 1..=n {
            for j in 1..=n {
                lu.set_value(i, j, a.value(a.lower_row() + i - 1, a.lower_col() + j - 1));
            }
        }
        let mut index = vec![0usize; n + 1];
        let mut d = 0.0;
        let mut vv = vec![0.0f64; n + 1];
        let error = lu_decompose(&mut lu, &mut index, &mut d, &mut vv, min_pivot);
        Self { lu, index, d, done: error == 0 }
    }

    /// `IsDone()` (`math_Gauss.hxx:53`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Solve(const math_Vector& B, math_Vector& X)`
    /// (`math_Gauss.cxx:39-47`): `X = B` then `LU_Solve`.
    pub fn solve(&self, b: &MathVector, x: &mut MathVector) -> Result<(), &'static str> {
        if !self.done {
            return Err("math_Gauss::Solve: not done");
        }
        *x = b.clone();
        lu_solve(&self.lu, &self.index, x);
        Ok(())
    }

    /// `Solve(math_Vector& X)` in place (`math_Gauss.cxx:49-58`), including
    /// the `Standard_DimensionError` on a length mismatch.
    pub fn solve_in_place(&self, x: &mut MathVector) -> Result<(), &'static str> {
        if !self.done {
            return Err("math_Gauss::Solve: not done");
        }
        if x.len() != self.lu.row_count() {
            return Err("math_Gauss::Solve: dimension error");
        }
        lu_solve(&self.lu, &self.index, x);
        Ok(())
    }

    /// `Determinant()` (`math_Gauss.cxx:60-68`).
    pub fn determinant(&self) -> Result<f64, &'static str> {
        if !self.done {
            return Err("math_Gauss::Determinant: not done");
        }
        let mut result = self.d;
        for j in 1..=self.lu.upper_row() {
            result *= self.lu.value(j, j);
        }
        Ok(result)
    }
}

/// `LU_Decompose(a, indx, d, vv, TINY)` (`math_Recipes.cxx:202-296`).
///
/// Returns `0` (`math_Status_OK`) or a non-zero status. The only status the
/// port distinguishes is `math_Status_SingularMatrix`.
fn lu_decompose(
    a: &mut MathMatrix,
    indx: &mut [usize],
    d: &mut f64,
    vv: &mut [f64],
    tiny: f64,
) -> i32 {
    let n = a.row_count();
    *d = 1.0;

    // `cxx:218-234`: the implicit row scaling.
    for i in 1..=n {
        let mut big = 0.0;
        for j in 1..=n {
            let temp = a.value(i, j).abs();
            if temp > big {
                big = temp;
            }
        }
        if big <= tiny {
            return 1; // math_Status_SingularMatrix
        }
        vv[i] = 1.0 / big;
    }

    for j in 1..=n {
        // `cxx:236-246`: upper part of the current column.
        for i in 1..j {
            let mut sum = a.value(i, j);
            for k in 1..i {
                sum -= a.value(i, k) * a.value(k, j);
            }
            a.set_value(i, j, sum);
        }
        // `cxx:247-266`: scaled pivot search. `continue` when `dum < big` so
        // that `imax` is updated even for NaN / Inf (`cxx:256-260`).
        let mut big = 0.0;
        let mut imax = 0usize;
        for i in j..=n {
            let mut sum = a.value(i, j);
            for k in 1..j {
                sum -= a.value(i, k) * a.value(k, j);
            }
            a.set_value(i, j, sum);
            let dum = vv[i] * sum.abs();
            if dum < big {
                continue;
            }
            big = dum;
            imax = i;
        }
        // `cxx:267-275`: row interchange.
        if j != imax {
            for k in 1..=n {
                let dum = a.value(imax, k);
                a.set_value(imax, k, a.value(j, k));
                a.set_value(j, k, dum);
            }
            *d = -*d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        // `cxx:277-280`.
        if a.value(j, j).abs() <= tiny {
            return 1; // math_Status_SingularMatrix
        }
        // `cxx:281-287`.
        if j != n {
            let dum = 1.0 / a.value(j, j);
            for i in (j + 1)..=n {
                a.set_value(i, j, a.value(i, j) * dum);
            }
        }
    }
    0 // math_Status_OK
}

/// `LU_Solve(a, indx, b)` (`math_Recipes.cxx:309-344`): `b` is replaced by the
/// solution. `nblow = b.Lower() - 1` reproduces the OCCT offset arithmetic.
fn lu_solve(a: &MathMatrix, indx: &[usize], b: &mut MathVector) {
    let n = a.row_count();
    // OCCT uses `nblow = b.Lower() - 1` and indexes `b(i + nblow)`; written as
    // `b.lower() + i - 1` it is the same offset without underflowing on a
    // vector whose `Lower()` is 0 (`Plate_Plate.cxx:473-474` builds those).
    let nblow = b.lower();
    let mut ii = 0usize;
    for i in 1..=n {
        let ip = indx[i];
        let mut sum = b.value(nblow + ip - 1);
        b.set_value(nblow + ip - 1, b.value(nblow + i - 1));
        if ii != 0 {
            for j in ii..i {
                sum -= a.value(i, j) * b.value(nblow + j - 1);
            }
        } else if sum != 0.0 {
            ii = i;
        }
        b.set_value(nblow + i - 1, sum);
    }
    for i in (1..=n).rev() {
        let mut sum = b.value(nblow + i - 1);
        for j in (i + 1)..=n {
            sum -= a.value(i, j) * b.value(nblow + j - 1);
        }
        b.set_value(nblow + i - 1, sum / a.value(i, i));
    }
}
