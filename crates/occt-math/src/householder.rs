//! `math_Householder` — Householder transformation of `A` into an orthogonal
//! times an upper-triangular matrix, with the right-hand sides transformed in
//! the same sweep.
//!
//! Source: `FoundationClasses/TKMath/math/math_Householder.cxx:38-179` and
//! `math_Householder.hxx` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `math_Householder(const math_Matrix& A, const math_Vector& B, EPS)`
//!   (`cxx:38-51`) and the row/column-window constructor (`cxx:66-84`) both
//!   funnel into the matrix constructor, so the port has one constructor
//!   [`Householder::new`] plus [`Householder::new_slice`] for the window arm.
//! - `Sol` is `math_Matrix(1, A.ColNumber(), 1, B.ColNumber())`; `Done` is
//!   `false` until the last statement of `Perform` (`cxx:86`, `cxx:177`).
//!   `Perform` bails out at `cxx:122-126` when `|g| <= EPS`, leaving `Sol`
//!   at whatever the zero-initialised matrix holds.
//! - `Q` is `math_Matrix(1, l, 1, n)` and only `n` columns are read, exactly
//!   as in OCCT.

use crate::MathMatrix;

/// `math_Householder` (`math_Householder.hxx`).
#[derive(Debug, Clone)]
pub struct Householder {
    /// `Sol` (`hxx`): the least-squares solution, `A.ColNumber()` rows.
    sol: MathMatrix,
    /// `Q` (`hxx`): the transformed system, upper-triangular in the first
    /// `A.ColNumber()` rows.
    q: MathMatrix,
    /// `Done` (`hxx`).
    done: bool,
}

impl Householder {
    /// `math_Householder(A, B, EPS)` (`cxx:52-64`): `A` is `m x n`, `B` is
    /// `m x p` and holds the right-hand sides column by column.
    pub fn new(a: &MathMatrix, b: &MathMatrix) -> Self {
        Self::new_with_eps(a, b, 1.0e-15)
    }

    /// Same as [`Householder::new`] with the caller's `EPS`
    /// (`math_Householder.hxx` default `1.0e-15`).
    pub fn new_with_eps(a: &MathMatrix, b: &MathMatrix, eps: f64) -> Self {
        let n = a.col_count();
        let l = a.row_count();
        let m = b.col_count();
        let mut this = Self {
            sol: MathMatrix::new(1, n, 1, m),
            q: MathMatrix::new(1, l, 1, n),
            done: false,
        };
        this.perform(a, b, 0, 0, eps);
        this
    }

    /// `math_Householder(A, B, lowerArow, upperArow, lowerAcol, upperAcol, EPS)`
    /// (`cxx:66-84`): the window of `A` actually used is
    /// `[lowerArow, upperArow] x [lowerAcol, upperAcol]`, and `B` is read from
    /// its own `LowerRow()` on.
    pub fn new_slice(
        a: &MathMatrix,
        b: &MathMatrix,
        lower_a_row: usize,
        upper_a_row: usize,
        lower_a_col: usize,
        upper_a_col: usize,
    ) -> Self {
        let n = upper_a_col - lower_a_col + 1;
        let l = upper_a_row - lower_a_row + 1;
        let m = b.col_count();
        let mut this = Self {
            sol: MathMatrix::new(1, n, 1, m),
            q: MathMatrix::new(1, l, 1, n),
            done: false,
        };
        // `mylowerArow` / `mylowerAcol` are the windows' lower bounds
        // (`cxx:78-81`); the copy below reads `A(i + mylowerArow - 1, ...)`.
        let row_offset = lower_a_row - 1;
        let col_offset = lower_a_col - 1;
        this.perform(a, b, row_offset, col_offset, 1.0e-15);
        this
    }

    /// `Perform(A, B, EPS)` (`cxx:86-179`).
    fn perform(&mut self, a: &MathMatrix, b: &MathMatrix, row_offset: usize, col_offset: usize, eps: f64) {
        let n = self.q.col_count();
        let l = self.q.row_count();
        let m = b.col_count();
        let mut b2 = MathMatrix::new(1, l, 1, m);

        let lbrow = b.lower_row();
        for i in 1..=l {
            for j in 1..=n {
                self.q.set_value(i, j, a.value(i + row_offset, j + col_offset));
            }
            for j in 1..=m {
                b2.set_value(i, j, b.value(i + lbrow - 1, j));
            }
        }
        // `Standard_DimensionError_Raise_if(l != B.RowNumber() || n > l, " ")`.
        if l != b.row_count() || n > l {
            panic!("Standard_DimensionError: math_Householder::Perform");
        }

        // Traitement de chaque colonne de A:
        for i in 1..=n {
            let mut h = 0.0;
            for k in i..=l {
                let qki = self.q.value(k, i);
                h += qki * qki;
            }
            let f = self.q.value(i, i);
            let g = if f < 1.0e-15 { h.sqrt() } else { -h.sqrt() };
            if g.abs() <= eps {
                self.done = false;
                return;
            }
            h -= f * g;
            let alfaii = g - f;
            for j in (i + 1)..=n {
                let mut scale = 0.0;
                for k in i..=l {
                    scale += self.q.value(k, i) * self.q.value(k, j);
                }
                let cj = (g * self.q.value(i, j) - scale) / h;
                let qij = self.q.value(i, j);
                self.q.set_value(i, j, qij - alfaii * cj);
                for k in (i + 1)..=l {
                    let v = self.q.value(k, j) + cj * self.q.value(k, i);
                    self.q.set_value(k, j, v);
                }
            }

            // Modification de B:
            for j in 1..=m {
                let mut scale = self.q.value(i, i) * b2.value(i, j);
                for k in (i + 1)..=l {
                    scale += self.q.value(k, i) * b2.value(k, j);
                }
                let cj = (g * b2.value(i, j) - scale) / h;
                let bij = b2.value(i, j);
                b2.set_value(i, j, bij - cj * alfaii);
                for k in (i + 1)..=l {
                    let v = b2.value(k, j) + cj * self.q.value(k, i);
                    b2.set_value(k, j, v);
                }
            }
            self.q.set_value(i, i, g);
        }

        // Remontee:
        for j in 1..=m {
            let v = b2.value(n, j) / self.q.value(n, n);
            self.sol.set_value(n, j, v);
            for i in (1..n).rev() {
                let mut scale = 0.0;
                for k in (i + 1)..=n {
                    scale += self.q.value(i, k) * self.sol.value(k, j);
                }
                let v = (b2.value(i, j) - scale) / self.q.value(i, i);
                self.sol.set_value(i, j, v);
            }
        }
        self.done = true;
    }

    /// `IsDone()` (`math_Householder.hxx`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `AllValues()` (`math_Householder.hxx`): the solution matrix, one column
    /// per right-hand side.
    pub fn all_values(&self) -> &MathMatrix {
        &self.sol
    }

    /// `Value(RowIndex, ColIndex)` (`math_Householder.hxx`).
    pub fn value(&self, row: usize, col: usize) -> f64 {
        self.sol.value(row, col)
    }

    /// The transformed matrix `Q` (`hxx`), exposed for the callers that read
    /// the triangular factor.
    pub fn q(&self) -> &MathMatrix {
        &self.q
    }

}
