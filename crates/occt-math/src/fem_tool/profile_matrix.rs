//! `FEmTool_SparseMatrix` / `FEmTool_ProfileMatrix` — the banded symmetric
//! matrix the variational assembly factors.
//!
//! Source: `ModelingData/TKGeomBase/FEmTool/FEmTool_SparseMatrix.hxx:28-64` and
//! `FEmTool_ProfileMatrix.cxx:31-321` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `profile(1, i)` / `profile(2, i)` are stored as two `Vec<i32>` with a
//!   leading unused slot (OCCT 1-based).
//! - `Decompose` is the in-profile Cholesky factorisation of
//!   `S * t(S) X = B`; the raw pointer walks `SMA` / `PM` by integer index, so
//!   the port indexes the backing `Vec` with `[k as usize]` where OCCT used
//!   `SMA[k]` on a pointer shifted by one.
//! - `Prepare` / the iterative `Solve` overload raise `Standard_NotImplemented`
//!   (`cxx:203-217`); the port returns `false` / does nothing.
//! - `OutM` / `OutS` (`cxx:279-321`) are OCCT_DEBUG dumps and are UNPORTED.

use crate::vector::MathVector;

/// `FEmTool_SparseMatrix` (`hxx:28-64`).
pub trait SparseMatrix {
    /// `Init(Value)`.
    fn init(&mut self, value: f64);

    /// `ChangeValue(I, J)`.
    fn change_value(&mut self, i: usize, j: usize) -> &mut f64;

    /// `Decompose()`.
    fn decompose(&mut self) -> bool;

    /// `Solve(B, X)`.
    fn solve(&self, b: &MathVector, x: &mut MathVector);

    /// `Prepare()`.
    fn prepare(&mut self) -> bool;

    /// `Multiplied(X, MX)`.
    fn multiplied(&self, x: &MathVector, mx: &mut MathVector);

    /// `RowNumber()`.
    fn row_number(&self) -> usize;

    /// `ColNumber()`.
    fn col_number(&self) -> usize;
}

/// `FEmTool_ProfileMatrix` (`FEmTool_ProfileMatrix.hxx:30-90`).
pub struct ProfileMatrix {
    /// `profile(1, i)` — half bandwidth of row `i`.
    bandwidth: Vec<i32>,
    /// `profile(2, i)` — address of the diagonal of row `i`.
    diagonal: Vec<i32>,
    /// `NextCoeff` (`1, profile(2, Length)`).
    next_coeff: Vec<i32>,
    /// `ProfileMatrix`.
    profile_matrix: Vec<f64>,
    /// `SMatrix` — the factor `S` of `S * t(S)`.
    s_matrix: Vec<f64>,
    is_decomp: bool,
}

impl ProfileMatrix {
    /// `FEmTool_ProfileMatrix(FirstIndexes)` (`cxx:32-69`).
    pub fn new(first_indexes: &[i32]) -> Self {
        let n = first_indexes.len();
        let mut bandwidth = vec![0i32; n + 1];
        let mut diagonal = vec![0i32; n + 1];
        bandwidth[1] = 0;
        diagonal[1] = 1;
        for i in 2..=n {
            bandwidth[i] = i as i32 - first_indexes[i - 1];
            diagonal[i] = diagonal[i - 1] + bandwidth[i] + 1;
        }

        let total = diagonal[n].max(0) as usize;
        let mut next_coeff = vec![0i32; total + 1];

        let mut k = 1usize;
        for i in 1..=n {
            for j in first_indexes[i - 1]..=(i as i32) {
                let mut l = i + 1;
                while l <= n && j < first_indexes[l - 1] {
                    l += 1;
                }
                next_coeff[k] = if l > n { 0 } else { l as i32 };
                k += 1;
            }
        }

        Self {
            bandwidth,
            diagonal,
            next_coeff,
            profile_matrix: vec![0.0; total + 1],
            s_matrix: vec![0.0; total + 1],
            is_decomp: false,
        }
    }

    /// `IsInProfile(i, j)` (`cxx:264-277`).
    pub fn is_in_profile(&self, i: usize, j: usize) -> bool {
        if j <= i {
            (i - j) as i32 <= self.bandwidth[i]
        } else {
            (j - i) as i32 <= self.bandwidth[j]
        }
    }
}

impl SparseMatrix for ProfileMatrix {
    /// `Init(Value)` (`cxx:73-77`).
    fn init(&mut self, value: f64) {
        self.profile_matrix.fill(value);
        self.is_decomp = false;
    }

    /// `ChangeValue(I, J)` (`cxx:81-97`).
    fn change_value(&mut self, i: usize, j: usize) -> &mut f64 {
        let mut ind = i as i32 - j as i32;
        let index = if ind < 0 {
            ind = -ind;
            assert!(ind <= self.bandwidth[j], "FEmTool_ProfileMatrix::ChangeValue");
            (self.diagonal[j] - ind) as usize
        } else {
            assert!(ind <= self.bandwidth[i], "FEmTool_ProfileMatrix::ChangeValue");
            (self.diagonal[i] - ind) as usize
        };
        &mut self.profile_matrix[index]
    }

    /// `Decompose()` (`cxx:101-150`).
    fn decompose(&mut self) -> bool {
        const EPS: f64 = 1.0e-32;

        self.s_matrix.fill(0.0);
        self.is_decomp = false;

        let row_number = self.row_number();
        for j in 1..=row_number {
            let diag_addr = self.diagonal[j];
            let kj = j as i32 - self.bandwidth[j];
            let mut sum = 0.0;
            for k in (diag_addr - self.bandwidth[j])..diag_addr {
                let v = self.s_matrix[k as usize];
                sum += v * v;
            }

            let mut a = self.profile_matrix[diag_addr as usize] - sum;
            if a < EPS {
                return false; // Matrix is not positive defined.
            }
            a = a.sqrt();
            self.s_matrix[diag_addr as usize] = a;

            let mut curr_addr = diag_addr;
            loop {
                let i = self.next_coeff[curr_addr as usize];
                if i <= 0 {
                    break;
                }
                curr_addr = self.diagonal[i as usize] - (i - j as i32);

                // Computation of Sum of S(i,k) * S(j,k) for k = 1,..,j-1.
                let mut sum = 0.0;
                let k_min = std::cmp::max(i - self.bandwidth[i as usize], kj);
                let mut ik = self.diagonal[i as usize] - i + k_min;
                let mut jk = diag_addr - j as i32 + k_min;
                let mut k = k_min;
                while k < j as i32 {
                    sum += self.s_matrix[ik as usize] * self.s_matrix[jk as usize];
                    k += 1;
                    ik += 1;
                    jk += 1;
                }
                self.s_matrix[curr_addr as usize] = (self.profile_matrix[curr_addr as usize] - sum) / a;
            }
        }
        self.is_decomp = true;
        true
    }

    /// `Solve(B, X)` (`cxx:155-201`): resolution of `S * t(S) X = B`.
    fn solve(&self, b: &MathVector, x: &mut MathVector) {
        if !self.is_decomp {
            return; // `StdFail_NotDone` (`cxx:159-162`).
        }

        let row_number = self.row_number();
        let col_number = self.col_number();

        // Resolution of S w = B.
        for i in 1..=row_number {
            let diag_addr = self.diagonal[i];
            let mut sum = 0.0;
            let mut j = i as i32 - self.bandwidth[i];
            let mut jj = diag_addr - (i as i32 - j);
            while j < i as i32 {
                sum += self.s_matrix[jj as usize] * x.value(x.lower() + j as usize - 1);
                j += 1;
                jj += 1;
            }
            let v = (b.value(b.lower() + i - 1) - sum) / self.s_matrix[diag_addr as usize];
            x.set_value(x.lower() + i - 1, v);
        }

        // Resolution of t(S) X = w.
        for i in (1..=col_number).rev() {
            let diag_addr = self.diagonal[i];
            let mut j = self.next_coeff[diag_addr as usize];
            let mut sum = 0.0;
            while j > 0 {
                let curr_addr = self.diagonal[j as usize] - (j - i as i32);
                sum += self.s_matrix[curr_addr as usize] * x.value(x.lower() + j as usize - 1);
                j = self.next_coeff[curr_addr as usize];
            }
            let v = (x.value(x.lower() + i - 1) - sum) / self.s_matrix[diag_addr as usize];
            x.set_value(x.lower() + i - 1, v);
        }
    }

    /// `Prepare()` (`cxx:203-206`): `Standard_NotImplemented`.
    fn prepare(&mut self) -> bool {
        false
    }

    /// `Multiplied(X, MX)` (`cxx:224-252`): `MX = H * X`.
    fn multiplied(&self, x: &MathVector, mx: &mut MathVector) {
        for i in 1..=self.row_number() {
            let diag_addr = self.diagonal[i];
            let mut m = 0.0;
            let mut j = i as i32 - self.bandwidth[i];
            let mut jj = diag_addr - (i as i32 - j);
            while j <= i as i32 {
                m += self.profile_matrix[jj as usize] * x.value(x.lower() + j as usize - 1);
                j += 1;
                jj += 1;
            }

            let mut curr_addr = diag_addr;
            let mut j = self.next_coeff[curr_addr as usize];
            while j > 0 {
                curr_addr = self.diagonal[j as usize] - (j - i as i32);
                m += self.profile_matrix[curr_addr as usize] * x.value(x.lower() + j as usize - 1);
                j = self.next_coeff[curr_addr as usize];
            }
            mx.set_value(mx.lower() + i - 1, m);
        }
    }

    /// `RowNumber()` (`cxx:254-257`).
    fn row_number(&self) -> usize {
        self.bandwidth.len() - 1
    }

    /// `ColNumber()` (`cxx:259-262`).
    fn col_number(&self) -> usize {
        self.bandwidth.len() - 1
    }
}
