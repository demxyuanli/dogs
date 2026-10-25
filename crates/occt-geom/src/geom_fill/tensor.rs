//! `GeomFill_Tensor` — the "gradient of gradient" storage the draft-location
//! law needs for its second derivative.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Tensor.{hxx,lxx}` and
//! `GeomFill_Tensor.cxx:20-52` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `NCollection_Array1<double> Tab` (`hxx:79`) is a flat `Vec<f64>`; the
//!   `Value` / `ChangeValue` index formula (`lxx:17-27`)
//!   `Mat + nbmat*(Col-1) + nbmtcl*(Row-1)` is kept verbatim.
//! - `Tab(Row, Col, Mat)` range-checks and raises `Standard_RangeError`; the
//!   port indexes the `Vec` directly, which panics on the same out-of-range
//!   use.

use occt_math::{MathMatrix, MathVector};

/// `GeomFill_Tensor` (`GeomFill_Tensor.hxx:31-82`).
pub struct Tensor {
    /// `Tab` (`hxx:79`), 1-based through [`Tensor::value`].
    tab: Vec<f64>,
    /// `nbrow` (`hxx:80`).
    nb_row: usize,
    /// `nbcol` (`hxx:81`).
    nb_col: usize,
    /// `nbmat` (`hxx:82`).
    nb_mat: usize,
}

impl Tensor {
    /// `GeomFill_Tensor(NbRow, NbCol, NbMat)` (`cxx:20-26`).
    pub fn new(nb_row: usize, nb_col: usize, nb_mat: usize) -> Self {
        Self {
            tab: vec![0.0; nb_row * nb_col * nb_mat],
            nb_row,
            nb_col,
            nb_mat,
        }
    }

    /// `Init(InitialValue)` (`cxx:28-34`).
    pub fn init(&mut self, initial_value: f64) {
        self.tab.fill(initial_value);
    }

    /// `Value(Row, Col, Mat)` (`lxx:17-21`).
    pub fn value(&self, row: usize, col: usize, mat: usize) -> f64 {
        self.tab[self.index(row, col, mat)]
    }

    /// `ChangeValue(Row, Col, Mat)` (`lxx:23-27`).
    pub fn set_value(&mut self, row: usize, col: usize, mat: usize, value: f64) {
        let i = self.index(row, col, mat);
        self.tab[i] = value;
    }

    /// `Tab(Mat + nbmat*(Col-1) + nbmtcl*(Row-1))` with `nbmtcl = nbmat*nbcol`.
    fn index(&self, row: usize, col: usize, mat: usize) -> usize {
        let nbmtcl = self.nb_mat * self.nb_col;
        mat + self.nb_mat * (col - 1) + nbmtcl * (row - 1)
    }

    /// `Multiply(Right, M)` (`cxx:36-52`): `M(i,j) = sum_k Value(i,j,k)*Right(k)`.
    pub fn multiply(&self, right: &MathVector, m: &mut MathMatrix) {
        for i in 1..=self.nb_row {
            for j in 1..=self.nb_col {
                let mut sum = 0.0;
                for k in 1..=self.nb_mat {
                    sum += self.value(i, j, k) * right.value(k);
                }
                m.set_value(i, j, sum);
            }
        }
    }
}
