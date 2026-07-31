//! Cholesky decomposition of symmetric positive-definite matrices.
//! Source: Numerical Recipes `choldc` (port of a `math_Cholesky`-style solver).

use crate::{MathMatrix, MathVector};

/// Cholesky decomposition: returns the lower-triangular factor `L` with `A = L · Lᵀ`.
///
/// Errors if `a` is not square, not symmetric (within tolerance), or not
/// positive-definite.
pub fn cholesky(a: &MathMatrix) -> Result<MathMatrix, String> {
    let n = a.row_count();
    if n != a.col_count() {
        return Err("Cholesky: matrix not square".to_string());
    }
    for i in 1..=n {
        for j in (i + 1)..=n {
            let tol = 1e-12 * (1.0 + a.value(i, j).abs());
            if (a.value(i, j) - a.value(j, i)).abs() > tol {
                return Err("Cholesky: matrix not symmetric".to_string());
            }
        }
    }
    let mut l = MathMatrix::new(1, n, 1, n);
    for j in 1..=n {
        let mut d = a.value(j, j);
        for k in 1..j {
            d -= l.value(j, k) * l.value(j, k);
        }
        if d <= 0.0 {
            return Err("Cholesky: matrix not positive-definite".to_string());
        }
        let dj = d.sqrt();
        l.set_value(j, j, dj);
        for i in (j + 1)..=n {
            let mut s = a.value(i, j);
            for k in 1..j {
                s -= l.value(i, k) * l.value(j, k);
            }
            l.set_value(i, j, s / dj);
        }
    }
    Ok(l)
}

/// Solve `A x = b` for a symmetric positive-definite `A` via Cholesky factors.
pub fn cholesky_solve(a: &MathMatrix, b: &MathVector) -> Result<MathVector, String> {
    let n = a.row_count();
    if b.len() != n {
        return Err("Cholesky solve: b dimension mismatch".to_string());
    }
    let l = cholesky(a)?;
    // Forward substitution: L y = b
    let mut y = MathVector::new(1, n);
    for i in 1..=n {
        let mut s = b.value(i);
        for k in 1..i {
            s -= l.value(i, k) * y.value(k);
        }
        y.set_value(i, s / l.value(i, i));
    }
    // Back substitution: Lᵀ x = y
    let mut x = MathVector::new(1, n);
    for i in (1..=n).rev() {
        let mut s = y.value(i);
        for k in (i + 1)..=n {
            s -= l.value(k, i) * x.value(k);
        }
        x.set_value(i, s / l.value(i, i));
    }
    Ok(x)
}

/// True if `a` is symmetric positive-definite.
pub fn is_positive_definite(a: &MathMatrix) -> bool {
    cholesky(a).is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mat2(data: &[&[f64]]) -> MathMatrix {
        let n = data.len();
        let mut m = MathMatrix::new(1, n, 1, n);
        for (i, row) in data.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                m.set_value(i + 1, j + 1, v);
            }
        }
        m
    }

    #[test]
    fn factor_roundtrip() {
        let a = mat2(&[&[4.0, 2.0], &[2.0, 3.0]]);
        let l = cholesky(&a).unwrap();
        let prod = l.multiplied_mat(&l.transposed());
        assert!((prod.value(1, 1) - 4.0).abs() < 1e-12);
        assert!((prod.value(1, 2) - 2.0).abs() < 1e-12);
        assert!((prod.value(2, 2) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn solve_2x2() {
        let a = mat2(&[&[4.0, 2.0], &[2.0, 3.0]]);
        let b = MathVector::from_slice(&[8.0, 8.0]); // A*[1,2]
        let x = cholesky_solve(&a, &b).unwrap();
        assert!((x.value(1) - 1.0).abs() < 1e-12);
        assert!((x.value(2) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn positive_definiteness() {
        let spd = mat2(&[&[4.0, 2.0], &[2.0, 3.0]]);
        assert!(is_positive_definite(&spd));
        let not_spd = mat2(&[&[1.0, 2.0], &[2.0, 1.0]]);
        assert!(!is_positive_definite(&not_spd));
    }
}
