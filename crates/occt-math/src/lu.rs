//! LU decomposition with partial pivoting.
//! Source: Numerical Recipes `ludcmp`/`lubksb` (port of a `math_LU`-style solver).

use crate::{MathMatrix, MathVector};

/// Result of [`lu_decompose`]: unit-lower-triangular `L`, upper-triangular `U`,
/// and the 1-based pivot row chosen at each elimination column.
#[derive(Debug, Clone)]
pub struct LU {
    pub l: MathMatrix,
    pub u: MathMatrix,
    pub piv: Vec<usize>,
}

/// Doolittle LU decomposition with partial pivoting: `A = P·L·U`.
///
/// Errors if `a` is not square or is singular (a zero pivot column appears).
pub fn lu_decompose(a: &MathMatrix) -> Result<LU, String> {
    let n = a.row_count();
    if n != a.col_count() {
        return Err("LU: matrix not square".to_string());
    }
    let mut m = a.clone(); // compact L+U during elimination
    let mut piv = vec![0usize; n + 1]; // 1-based

    for k in 1..=n {
        let mut max_val = m.value(k, k).abs();
        let mut piv_row = k;
        for i in (k + 1)..=n {
            let v = m.value(i, k).abs();
            if v > max_val {
                max_val = v;
                piv_row = i;
            }
        }
        if max_val < 1e-30 {
            return Err("LU: singular matrix".to_string());
        }
        if piv_row != k {
            for j in 1..=n {
                let t = m.value(k, j);
                m.set_value(k, j, m.value(piv_row, j));
                m.set_value(piv_row, j, t);
            }
        }
        piv[k] = piv_row;
        for i in (k + 1)..=n {
            let factor = m.value(i, k) / m.value(k, k);
            m.set_value(i, k, factor);
            for j in (k + 1)..=n {
                let v = m.value(i, j) - factor * m.value(k, j);
                m.set_value(i, j, v);
            }
        }
    }

    let mut l = MathMatrix::new(1, n, 1, n);
    let mut u = MathMatrix::new(1, n, 1, n);
    for i in 1..=n {
        for j in 1..=n {
            if i > j {
                l.set_value(i, j, m.value(i, j));
            } else {
                u.set_value(i, j, m.value(i, j));
            }
        }
        l.set_value(i, i, 1.0);
    }
    Ok(LU { l, u, piv })
}

/// Solve `A x = b` using the LU factors (permute `b`, forward/back substitute).
pub fn lu_solve(lu: &LU, b: &MathVector) -> Result<MathVector, String> {
    let n = lu.l.row_count();
    if b.len() != n {
        return Err("LU solve: b dimension mismatch".to_string());
    }
    let mut x = b.clone();
    // Apply the pivot permutation to the right-hand side.
    for k in 1..=n {
        let p = lu.piv[k];
        if p != k {
            let t = x.value(k);
            x.set_value(k, x.value(p));
            x.set_value(p, t);
        }
    }
    // Forward substitution: L y = P b (L unit lower triangular).
    for i in 1..=n {
        let mut s = x.value(i);
        for j in 1..i {
            s -= lu.l.value(i, j) * x.value(j);
        }
        x.set_value(i, s);
    }
    // Back substitution: U x = y.
    for i in (1..=n).rev() {
        let mut s = x.value(i);
        for j in (i + 1)..=n {
            s -= lu.u.value(i, j) * x.value(j);
        }
        x.set_value(i, s / lu.u.value(i, i));
    }
    Ok(x)
}

/// Determinant of the original matrix from the LU factors.
pub fn det_from_lu(lu: &LU) -> f64 {
    let n = lu.l.row_count();
    let mut det = 1.0;
    for k in 1..=n {
        if lu.piv[k] != k {
            det = -det;
        }
        det *= lu.u.value(k, k);
    }
    det
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
    fn lu_solve_3x3() {
        let a = mat2(&[&[2.0, 1.0, 0.0], &[1.0, 3.0, 1.0], &[0.0, 1.0, 2.0]]);
        let lu = lu_decompose(&a).unwrap();
        let b = MathVector::from_slice(&[5.0, 11.0, 5.0]); // x = [1,3,1]
        let x = lu_solve(&lu, &b).unwrap();
        assert!((x.value(1) - 1.0).abs() < 1e-12);
        assert!((x.value(2) - 3.0).abs() < 1e-12);
        assert!((x.value(3) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn lu_requires_pivoting() {
        let a = mat2(&[&[0.0, 1.0], &[1.0, 0.0]]); // zero pivot without row swap
        let lu = lu_decompose(&a).unwrap();
        let b = MathVector::from_slice(&[3.0, 4.0]); // x = [4,3]
        let x = lu_solve(&lu, &b).unwrap();
        assert!((x.value(1) - 4.0).abs() < 1e-12);
        assert!((x.value(2) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn determinant() {
        let a = mat2(&[&[4.0, 7.0], &[2.0, 6.0]]);
        let lu = lu_decompose(&a).unwrap();
        assert!((det_from_lu(&lu) - 10.0).abs() < 1e-12);
        let a2 = mat2(&[&[0.0, 1.0], &[1.0, 0.0]]);
        let lu2 = lu_decompose(&a2).unwrap();
        assert!((det_from_lu(&lu2) + 1.0).abs() < 1e-12);
    }
}
