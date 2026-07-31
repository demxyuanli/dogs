//! Extended matrix operations built on `MathMatrix` / `MathVector`.
//! Source: `math_Recipes.hxx`, `math_Gauss.hxx` (determinant, inversion, rank).
//!
//! All matrices use the crate's 1-based indexing convention.

use crate::lu::{lu_decompose, lu_solve};
use crate::matrix::MathMatrix;
use crate::vector::MathVector;

/// Identity matrix of size `n × n` (1-based bounds 1..n).
pub fn identity(n: usize) -> MathMatrix {
    let mut m = MathMatrix::new(1, n, 1, n);
    for i in 1..=n {
        m.set_value(i, i, 1.0);
    }
    m
}

/// Transpose of `m`.
pub fn transpose(m: &MathMatrix) -> MathMatrix {
    let rows = m.row_count();
    let cols = m.col_count();
    let mut t = MathMatrix::new(1, cols, 1, rows);
    for r in 1..=rows {
        for c in 1..=cols {
            t.set_value(c, r, m.value(r, c));
        }
    }
    t
}

/// Matrix multiplication `a · b` (column-major storage is hidden; sizes are
/// checked). Source: `math_Matrix::Multiply`.
pub fn matmul(a: &MathMatrix, b: &MathMatrix) -> Result<MathMatrix, String> {
    if a.col_count() != b.row_count() {
        return Err(format!(
            "matmul: dimension mismatch {}x{} · {}x{}",
            a.row_count(), a.col_count(), b.row_count(), b.col_count()
        ));
    }
    let n = a.col_count();
    let mut r = MathMatrix::new(1, a.row_count(), 1, b.col_count());
    for i in 1..=a.row_count() {
        for j in 1..=b.col_count() {
            let mut s = 0.0;
            for k in 1..=n {
                s += a.value(i, k) * b.value(k, j);
            }
            r.set_value(i, j, s);
        }
    }
    Ok(r)
}

/// Trace of a square matrix.
pub fn trace(m: &MathMatrix) -> f64 {
    let n = m.row_count().min(m.col_count());
    let mut s = 0.0;
    for i in 1..=n {
        s += m.value(i, i);
    }
    s
}

/// Determinant via LU decomposition with partial pivoting.
pub fn determinant(m: &MathMatrix) -> Result<f64, String> {
    if !m.is_square() {
        return Err("determinant: matrix must be square".into());
    }
    // The LU module's det_from_lu already exists; recompute here for a
    // self-contained formula in case the caller wants the pivot sign.
    let n = m.row_count();
    let mut a = m.clone();
    let mut det = 1.0;
    let mut sign = 1.0;
    for k in 1..=n {
        // Partial pivot.
        let mut p = k;
        let mut maxv = a.value(k, k).abs();
        for i in (k + 1)..=n {
            let v = a.value(i, k).abs();
            if v > maxv {
                maxv = v;
                p = i;
            }
        }
        if maxv < 1e-300 {
            return Ok(0.0);
        }
        if p != k {
            for c in k..=n {
                let tmp = a.value(k, c);
                a.set_value(k, c, a.value(p, c));
                a.set_value(p, c, tmp);
            }
            sign = -sign;
        }
        let piv = a.value(k, k);
        det *= piv;
        for i in (k + 1)..=n {
            let f = a.value(i, k) / piv;
            for c in (k + 1)..=n {
                a.set_value(i, c, a.value(i, c) - f * a.value(k, c));
            }
        }
    }
    Ok(sign * det)
}

/// Inverse via LU decomposition with partial pivoting.
pub fn inverse(m: &MathMatrix) -> Result<MathMatrix, String> {
    if !m.is_square() {
        return Err("inverse: matrix must be square".into());
    }
    let n = m.row_count();
    let lu = lu_decompose(m)?;
    let mut inv = MathMatrix::new(1, n, 1, n);
    for c in 1..=n {
        let mut e = MathVector::new(1, n);
        e.set_value(c, 1.0);
        let col = lu_solve(&lu, &e)?;
        for r in 1..=n {
            inv.set_value(r, c, col.value(r));
        }
    }
    Ok(inv)
}

/// Frobenius norm: sqrt(Σ a_ij²).
pub fn frobenius_norm(m: &MathMatrix) -> f64 {
    let mut s = 0.0;
    for r in 1..=m.row_count() {
        for c in 1..=m.col_count() {
            s += m.value(r, c) * m.value(r, c);
        }
    }
    s.sqrt()
}

/// Matrix rank via Gaussian elimination with partial pivoting.
pub fn matrix_rank(m: &MathMatrix, tol: f64) -> usize {
    let mut a = m.clone();
    let rows = a.row_count();
    let cols = a.col_count();
    let mut rank = 0;
    let mut r = 1;
    let mut c = 1;
    while r <= rows && c <= cols {
        // Find pivot.
        let mut pr = r;
        let mut maxv = a.value(r, c).abs();
        for i in (r + 1)..=rows {
            let v = a.value(i, c).abs();
            if v > maxv {
                maxv = v;
                pr = i;
            }
        }
        if maxv <= tol {
            c += 1;
            continue;
        }
        // Swap rows.
        for cc in c..=cols {
            let tmp = a.value(r, cc);
            a.set_value(r, cc, a.value(pr, cc));
            a.set_value(pr, cc, tmp);
        }
        // Eliminate below.
        let piv = a.value(r, c);
        for i in (r + 1)..=rows {
            let f = a.value(i, c) / piv;
            for cc in c..=cols {
                a.set_value(i, cc, a.value(i, cc) - f * a.value(r, cc));
            }
        }
        rank += 1;
        r += 1;
        c += 1;
    }
    rank
}

/// Direct 3×3 determinant formula.
pub fn det3(
    a: f64, b: f64, c: f64,
    d: f64, e: f64, f: f64,
    g: f64, h: f64, i: f64,
) -> f64 {
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

/// Cross product in 2D (scalar result): `det[[a, b], [c, d]]`.
pub fn det2(a: f64, b: f64, c: f64, d: f64) -> f64 {
    a * d - b * c
}

/// Solve the linear system `a · x = b` (square) via LU. Convenience wrapper.
pub fn solve_linear(a: &MathMatrix, b: &MathVector) -> Result<MathVector, String> {
    if !a.is_square() {
        return Err("solve_linear: matrix must be square".into());
    }
    let lu = lu_decompose(a)?;
    lu_solve(&lu, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mat2(a: f64, b: f64, c: f64, d: f64) -> MathMatrix {
        let mut m = MathMatrix::new(1, 2, 1, 2);
        m.set_value(1, 1, a);
        m.set_value(1, 2, b);
        m.set_value(2, 1, c);
        m.set_value(2, 2, d);
        m
    }

    #[test]
    fn identity_and_trace() {
        let id = identity(3);
        assert_eq!(trace(&id), 3.0);
        assert_eq!(id.value(2, 2), 1.0);
        assert_eq!(id.value(1, 3), 0.0);
    }

    #[test]
    fn transpose_and_matmul() {
        let a = mat2(1.0, 2.0, 3.0, 4.0);
        let t = transpose(&a);
        assert_eq!(t.value(1, 2), 3.0);
        assert_eq!(t.value(2, 1), 2.0);
        // A · A⁻¹ round trip via matmul of A and its inverse.
        let inv = inverse(&a).unwrap();
        let prod = matmul(&a, &inv).unwrap();
        assert!((prod.value(1, 1) - 1.0).abs() < 1e-12);
        assert!((prod.value(1, 2) - 0.0).abs() < 1e-12);
        assert!((prod.value(2, 2) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn determinant_2x2_and_3x3() {
        let a = mat2(1.0, 2.0, 3.0, 4.0);
        assert!((determinant(&a).unwrap() + 2.0).abs() < 1e-12);
        assert!((det2(1.0, 2.0, 3.0, 4.0) + 2.0).abs() < 1e-12);
        assert!((det3(1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn inverse_roundtrip() {
        let a = mat2(4.0, 7.0, 2.0, 6.0);
        let inv = inverse(&a).unwrap();
        let prod = matmul(&a, &inv).unwrap();
        assert!((prod.value(1, 1) - 1.0).abs() < 1e-12);
        assert!((prod.value(2, 2) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rank_and_solve() {
        // Rank 2 vs rank 1.
        let full = mat2(1.0, 2.0, 3.0, 4.0);
        assert_eq!(matrix_rank(&full, 1e-12), 2);
        let rank1 = mat2(1.0, 2.0, 2.0, 4.0);
        assert_eq!(matrix_rank(&rank1, 1e-12), 1);

        // Solve: x + 2y = 5, 3x + 4y = 11 → x=1, y=2.
        let b = {
            let mut v = MathVector::new(1, 2);
            v.set_value(1, 5.0);
            v.set_value(2, 11.0);
            v
        };
        let x = solve_linear(&full, &b).unwrap();
        assert!((x.value(1) - 1.0).abs() < 1e-10);
        assert!((x.value(2) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn frobenius() {
        let a = mat2(3.0, 0.0, 0.0, 4.0);
        assert!((frobenius_norm(&a) - 5.0).abs() < 1e-12);
    }
}
