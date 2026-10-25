//! `math_Crout` — inversion of a symmetric matrix by the Crout method,
//! `A = L * D * T(L)` with `L` lower triangular and `D` diagonal.
//!
//! Source: `FoundationClasses/TKMath/math/math_Crout.cxx:31-208`,
//! `math_Crout.hxx` and `math_Crout.lxx` (OCCT 8.0.0).
//!
//! Translation notes:
//! - Only the lower triangle of `A` (plus the diagonal) is read
//!   (`cxx:45-68`), and only the lower triangle of `InvA` is written
//!   (`cxx:89-107`); `Solve` (`cxx:113-135`) re-reads the mirror element
//!   `InvA(j, i)` for `j > i`. The port keeps that storage convention
//!   verbatim instead of mirroring the triangle.
//! - There is no pivoting (unlike `crout.rs`, which ports the Numerical
//!   Recipes `ludcmp` used by `math_Recipes`); the singularity test is
//!   `|Diag(i)| <= MinPivot` with `MinPivot = 1.0e-20` (`hxx`), and it aborts
//!   the decomposition at `cxx:63-67` leaving `InvA` zero-initialised.
//! - `Determinant`, `Inverse` and `Invert` raise `StdFail_NotDone` when
//!   `IsDone()` is false (`math_Crout.lxx:27-44`); the port panics.

use crate::{MathMatrix, MathVector};

/// `math_Crout` (`math_Crout.hxx`).
#[derive(Debug, Clone)]
pub struct MathCrout {
    /// `InvA` (`hxx`): the inverse of `A`, lower triangle only.
    inv_a: MathMatrix,
    /// `Det` (`hxx`): the product of the Crout diagonals.
    det: f64,
    /// `Done` (`hxx`).
    done: bool,
}

impl MathCrout {
    /// `math_Crout(A, MinPivot = 1.0e-20)` (`cxx:31-110`).
    pub fn new(a: &MathMatrix) -> Self {
        Self::new_with_min_pivot(a, 1.0e-20)
    }

    /// Same as [`MathCrout::new`] with the caller's `MinPivot` (`hxx:47`).
    pub fn new_with_min_pivot(a: &MathMatrix, min_pivot: f64) -> Self {
        let nctl = a.row_count();
        let lowr = a.lower_row();
        let lowc = a.lower_col();

        let mut this = Self {
            // `math_Crout` accepts an empty `A`; `MathMatrix` needs
            // `upper >= lower`, so the empty shape keeps a 1x1 stand-in that
            // the (empty) loops below never touch.
            inv_a: MathMatrix::with_init(1, nctl.max(1), 1, a.col_count().max(1), 0.0),
            det: 1.0,
            done: false,
        };

        // `math_NotSquare_Raise_if(Nctl != A.ColNumber(), " ")` (`cxx:42`).
        if nctl != a.col_count() {
            panic!("math_NotSquare: math_Crout");
        }

        // `math_Matrix L(1, Nctl, 1, Nctl); math_Vector Diag(1, Nctl);`
        // (`cxx:39-40`), both zero-initialised.
        let mut l = MathMatrix::with_init(1, nctl.max(1), 1, nctl.max(1), 0.0);
        let mut diag = MathVector::with_init(1, nctl.max(1), 0.0);

        // `Det = 1` (`cxx:44`).
        this.det = 1.0;
        for i in 1..=nctl {
            for j in 1..i {
                let mut scale = 0.0;
                for k in 1..j {
                    scale += l.value(i, k) * l.value(j, k) * diag.value(k);
                }
                // `L(i, j) = (A(i + lowr - 1, j + lowc - 1) - scale) / Diag(j)`
                // (`cxx:57`).
                let aij = a.value(i + lowr - 1, j + lowc - 1);
                l.set_value(i, j, (aij - scale) / diag.value(j));
            }
            let mut scale = 0.0;
            for k in 1..i {
                scale += l.value(i, k) * l.value(i, k) * diag.value(k);
            }
            // `Diag(i) = A(i + lowr - 1, i + lowc - 1) - scale` (`cxx:61`).
            let aii = a.value(i + lowr - 1, i + lowc - 1);
            diag.set_value(i, aii - scale);
            this.det *= diag.value(i);
            // `if (std::abs(Diag(i)) <= MinPivot) { Done = false; return; }`
            // (`cxx:63-67`).
            if diag.value(i).abs() <= min_pivot {
                this.done = false;
                return this;
            }
            l.set_value(i, i, 1.0);
        }

        // Calcul de l inverse de L (`cxx:71-86`):
        l.set_value(1, 1, 1.0 / l.value(1, 1));
        for i in 2..=nctl {
            for k in 1..i {
                let mut scale = 0.0;
                for j in k..i {
                    scale += l.value(i, j) * l.value(j, k);
                }
                l.set_value(i, k, -scale / l.value(i, i));
            }
            l.set_value(i, i, 1.0 / l.value(i, i));
        }

        // Calcul de l inverse de Mat (`cxx:89-108`): only the lower triangle
        // of `InvA` is assigned.
        for j in 1..=nctl {
            let mut scale = l.value(j, j) * l.value(j, j) / diag.value(j);
            for k in (j + 1)..=nctl {
                scale += l.value(k, j) * l.value(k, j) / diag.value(k);
            }
            this.inv_a.set_value(j, j, scale);
            for i in (j + 1)..=nctl {
                let mut scale = l.value(i, j) * l.value(i, i) / diag.value(i);
                for k in (i + 1)..=nctl {
                    scale += l.value(k, j) * l.value(k, i) / diag.value(k);
                }
                this.inv_a.set_value(i, j, scale);
            }
        }
        this.done = true;
        this
    }

    /// `IsDone()` (`math_Crout.lxx:17-20`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Solve(B, X)` (`cxx:113-135`): solution of `A . X = B` using the
    /// stored lower triangle (`InvA(i, j)` for `j <= i`, `InvA(j, i)`
    /// otherwise).
    pub fn solve(&self, b: &MathVector, x: &mut MathVector) {
        if !self.done {
            panic!("StdFail_NotDone: math_Crout::Solve");
        }
        let n = self.inv_a.row_count();
        if b.len() != n || x.len() != b.len() {
            panic!("Standard_DimensionError: math_Crout::Solve");
        }
        let lowb = b.lower();
        let lowx = x.lower();

        for i in 1..=n {
            let mut acc = self.inv_a.value(i, 1) * b.value(1 + lowb - 1);
            for j in 2..=i {
                acc += self.inv_a.value(i, j) * b.value(j + lowb - 1);
            }
            for j in (i + 1)..=n {
                acc += self.inv_a.value(j, i) * b.value(j + lowb - 1);
            }
            x.set_value(i + lowx - 1, acc);
        }
    }

    /// `Determinant()` (`math_Crout.lxx:27-32`).
    pub fn determinant(&self) -> f64 {
        if !self.done {
            panic!("StdFail_NotDone: math_Crout::Determinant");
        }
        self.det
    }

    /// `Inverse()` (`math_Crout.lxx:34-38`): the inverse, lower triangle only.
    pub fn inverse(&self) -> &MathMatrix {
        if !self.done {
            panic!("StdFail_NotDone: math_Crout::Inverse");
        }
        &self.inv_a
    }

    /// `Invert(Inv)` (`math_Crout.lxx:40-44`).
    pub fn invert(&self) -> MathMatrix {
        if !self.done {
            panic!("StdFail_NotDone: math_Crout::Invert");
        }
        self.inv_a.clone()
    }
}
