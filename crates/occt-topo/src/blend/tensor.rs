//! `BlendFunc_Tensor` (`BlendFunc_Tensor.hxx:29-70`, `.lxx:17-24`,
//! `.cxx:21-53`).
//!
//! The "gradient of gradient" storage of the `BlendFunc` solvers: a
//! `nbrow x nbcol x nbmat` cube of doubles flattened into one array, plus the
//! single `Multiply` that contracts it with a `math_Vector` into a
//! `math_Matrix`.
//!
//! Translation notes:
//! - `NCollection_Array1<double> Tab(1, NbRow * NbMat * NbCol)` (`.cxx:23`)
//!   starts empty in OCCT (the constructor leaves the elements untouched); the
//!   port zero-fills, which is the only value a reader can see before
//!   [`Tensor::init`],
//! - `Value` / `ChangeValue` (`.lxx:17-24`) keep OCCT's 1-based layout
//!   `Mat + nbmat * (Col - 1) + nbmtcl * (Row - 1)`, including the `nbmtcl`
//!   stride that is `NbMat * NbCol`,
//! - `operator()` on a `math_Vector` / `math_Matrix` is 1-based there, so
//!   [`Tensor::multiply`] indexes `right` / `product` the same way.

use occt_math::{MathMatrix, MathVector};

/// `BlendFunc_Tensor` (`hxx:29-70`).
#[derive(Debug, Clone)]
pub struct Tensor {
    /// `Tab`, 1-based in OCCT, 0-based here.
    tab: Vec<f64>,
    /// `nbrow`.
    nbrow: usize,
    /// `nbcol`.
    nbcol: usize,
    /// `nbmat`.
    nbmat: usize,
    /// `nbmtcl` = `nbmat * nbcol`.
    nbmtcl: usize,
}

impl Tensor {
    /// `BlendFunc_Tensor(NbRow, NbCol, NbMat)` (`.cxx:21-28`).
    pub fn new(nb_row: usize, nb_col: usize, nb_mat: usize) -> Self {
        Self {
            tab: vec![0.0; nb_row * nb_mat * nb_col],
            nbrow: nb_row,
            nbcol: nb_col,
            nbmat: nb_mat,
            nbmtcl: nb_mat * nb_col,
        }
    }

    /// `Init(InitialValue)` (`.cxx:31-37`): `Tab.Init(InitialValue)` over the
    /// whole cube.
    pub fn init(&mut self, initial_value: f64) {
        self.tab.fill(initial_value);
    }

    /// `Value(Row, Col, Mat)` (`.lxx:17-20`).
    pub fn value(&self, row: usize, col: usize, mat: usize) -> f64 {
        self.tab[self.index(row, col, mat)]
    }

    /// `ChangeValue(Row, Col, Mat)` (`.lxx:22-24`).
    pub fn change_value(&mut self, row: usize, col: usize, mat: usize) -> &mut f64 {
        let i = self.index(row, col, mat);
        &mut self.tab[i]
    }

    /// Sets `(Row, Col, Mat)`; the store side of `ChangeValue`.
    pub fn set_value(&mut self, row: usize, col: usize, mat: usize, value: f64) {
        let i = self.index(row, col, mat);
        self.tab[i] = value;
    }

    /// `Multiply(Right, M)` (`.cxx:39-53`): `M(i, j) = sum_k Value(i, j, k) *
    /// Right(k)`, with `i` over `nbrow` and `j` over `nbcol`.
    ///
    /// `M` is written in place and keeps whatever it held outside
    /// `1..nbrow x 1..nbcol`, exactly like the C++ loop.
    pub fn multiply(&self, right: &MathVector, m: &mut MathMatrix) {
        for i in 1..=self.nbrow {
            for j in 1..=self.nbcol {
                let mut somme = 0.0;
                for k in 1..=self.nbmat {
                    somme += self.value(i, j, k) * right.value(k);
                }
                m.set_value(i, j, somme);
            }
        }
    }

    /// `Tab(Mat + nbmat * (Col - 1) + nbmtcl * (Row - 1))`, 1-based to 0-based.
    fn index(&self, row: usize, col: usize, mat: usize) -> usize {
        assert!(
            row >= 1 && row <= self.nbrow,
            "BlendFunc_Tensor::Value: Row out of range"
        );
        assert!(
            col >= 1 && col <= self.nbcol,
            "BlendFunc_Tensor::Value: Col out of range"
        );
        assert!(
            mat >= 1 && mat <= self.nbmat,
            "BlendFunc_Tensor::Value: Mat out of range"
        );
        mat + self.nbmat * (col - 1) + self.nbmtcl * (row - 1) - 1
    }

    /// `NbRow`.
    pub fn nb_row(&self) -> usize {
        self.nbrow
    }

    /// `NbCol`.
    pub fn nb_col(&self) -> usize {
        self.nbcol
    }

    /// `NbMat`.
    pub fn nb_mat(&self) -> usize {
        self.nbmat
    }
}
