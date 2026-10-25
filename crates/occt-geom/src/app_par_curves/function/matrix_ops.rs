//! Bound-correct `math_Matrix` operations used by `AppParCurves_Function`.
//!
//! `MathMatrix::multiplied_mat` / `MathMatrix::transposed` in this crate read
//! elements with 1-based indices and only transpose square matrices, while
//! OCCT's `math_Matrix::Multiply` / `math_Matrix::Transposed`
//! (`math_Matrix.cxx`) honour each operand's own lower bounds. The constraint
//! matrices here start at `FirstPoint` rows, so the helpers below reproduce
//! the OCCT bounds instead of the crate's 1-based subset.

use occt_math::{MathMatrix, MathVector};

/// `math_Matrix::Transposed()` (`math_Matrix.lxx`): the result has the
/// operand's column bounds as its row bounds, and vice versa.
pub(super) fn transposed(m: &MathMatrix) -> MathMatrix {
    let mut t = MathMatrix::new(m.lower_col(), m.upper_col(), m.lower_row(), m.upper_row());
    for i in m.lower_row()..=m.upper_row() {
        for j in m.lower_col()..=m.upper_col() {
            t.set_value(j, i, m.value(i, j));
        }
    }
    t
}

/// `math_Matrix::Multiply(const math_Matrix&)` (`math_Matrix.cxx`): `a * b`
/// with `a`'s row bounds and `b`'s column bounds.
pub(super) fn mat_mul(a: &MathMatrix, b: &MathMatrix) -> MathMatrix {
    if a.col_count() != b.row_count() {
        panic!("Standard_DimensionError: math_Matrix::Multiply");
    }
    let mut r = MathMatrix::new(a.lower_row(), a.upper_row(), b.lower_col(), b.upper_col());
    for i in a.lower_row()..=a.upper_row() {
        for j in b.lower_col()..=b.upper_col() {
            let mut sum = 0.0;
            for k in a.lower_col()..=a.upper_col() {
                sum += a.value(i, k) * b.value(k, j);
            }
            r.set_value(i, j, sum);
        }
    }
    r
}

/// `math_Matrix::Multiply(const math_Vector&)` (`math_Matrix.cxx`): `m * v`,
/// a vector over the matrix' row bounds.
pub(super) fn mat_vec(m: &MathMatrix, v: &MathVector) -> MathVector {
    if m.col_count() != v.len() {
        panic!("Standard_DimensionError: math_Matrix::Multiply");
    }
    let mut r = MathVector::new(m.lower_row(), m.upper_row());
    for i in m.lower_row()..=m.upper_row() {
        let mut sum = 0.0;
        for j in m.lower_col()..=m.upper_col() {
            sum += m.value(i, j) * v.value(j);
        }
        r.set_value(i, sum);
    }
    r
}

/// `math_Matrix::Col(const int)` (`math_Matrix.cxx`).
pub(super) fn col(m: &MathMatrix, j: usize) -> MathVector {
    let mut v = MathVector::new(m.lower_row(), m.upper_row());
    for i in m.lower_row()..=m.upper_row() {
        v.set_value(i, m.value(i, j));
    }
    v
}

/// `math_Matrix::Row(const int)` (`math_Matrix.cxx`).
pub(super) fn row(m: &MathMatrix, i: usize) -> MathVector {
    let mut v = MathVector::new(m.lower_col(), m.upper_col());
    for j in m.lower_col()..=m.upper_col() {
        v.set_value(j, m.value(i, j));
    }
    v
}

/// `math_Vector operator*(const math_Vector&, const math_Vector&)` — the dot
/// product (`math_Vector.lxx`), with OCCT's dimension check.
pub(super) fn dot(a: &MathVector, b: &MathVector) -> f64 {
    if a.len() != b.len() {
        panic!("Standard_DimensionError: math_Vector::operator*");
    }
    let mut sum = 0.0;
    for i in a.lower()..=a.upper() {
        sum += a.value(i) * b.value(i);
    }
    sum
}

/// `math_Matrix::Inverse()` (`math_Matrix.cxx`): Gauss inversion; the port
/// panics where OCCT raises `StdFail_NotDone`.
pub(super) fn inverse(m: &MathMatrix) -> MathMatrix {
    match m.clone().inverted() {
        Ok(inv) => inv,
        Err(_) => panic!("StdFail_NotDone: math_Matrix::Inverse"),
    }
}
