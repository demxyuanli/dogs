//! Sparse square matrix and iterative solvers.
//! Source: OCCT `math_SparseSingularValues`, `math_SparseMatrix`.

/// Row-wise sparse square matrix (CSR-ish). 0-based row/column indices.
/// Only non-zero entries are stored; each row is kept sorted by column.
#[derive(Debug, Clone)]
pub struct SparseMatrix {
    pub n: usize,
    pub rows: Vec<Vec<(usize, f64)>>,
}

impl SparseMatrix {
    /// Empty `n × n` sparse matrix.
    pub fn new(n: usize) -> Self {
        SparseMatrix { n, rows: vec![Vec::new(); n] }
    }

    /// Set entry `(i, j)` to `v`. A zero value removes the entry.
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        let row = &mut self.rows[i];
        if v == 0.0 {
            if let Some(pos) = row.iter().position(|&(cj, _)| cj == j) {
                row.remove(pos);
            }
            return;
        }
        match row.iter().position(|&(cj, _)| cj >= j) {
            Some(pos) if row[pos].0 == j => row[pos].1 = v,
            Some(pos) => row.insert(pos, (j, v)),
            None => row.push((j, v)),
        }
    }

    /// Value at `(i, j)` (0.0 if not stored).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.rows[i].iter().find(|&&(cj, _)| cj == j).map_or(0.0, |&(_, v)| v)
    }

    /// Number of stored (non-zero) entries.
    pub fn nnz(&self) -> usize {
        self.rows.iter().map(|r| r.len()).sum()
    }

    /// Euclidean (L2) norm of row `i`.
    pub fn row_norm(&self, i: usize) -> f64 {
        self.rows[i].iter().map(|&(_, v)| v * v).sum::<f64>().sqrt()
    }

    /// True if `A(i,j) == A(j,i)` for every stored entry, within `tol`.
    pub fn is_symmetric(&self, tol: f64) -> bool {
        for i in 0..self.n {
            for &(j, v) in &self.rows[i] {
                if (self.get(j, i) - v).abs() > tol {
                    return false;
                }
            }
        }
        true
    }

    /// Build a sparse matrix from a dense `MathMatrix` (non-zeros only).
    pub fn from_dense(m: &crate::matrix::MathMatrix) -> SparseMatrix {
        let rows = m.row_count();
        let cols = m.col_count();
        let mut s = SparseMatrix::new(rows);
        for r in 1..=rows {
            for c in 1..=cols {
                let v = m.value(r, c);
                if v != 0.0 {
                    s.set(r - 1, c - 1, v);
                }
            }
        }
        s
    }

    /// `y = A · x`.
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.n, "matvec: dimension mismatch");
        let mut y = vec![0.0; self.n];
        for i in 0..self.n {
            let mut s = 0.0;
            for &(j, v) in &self.rows[i] {
                s += v * x[j];
            }
            y[i] = s;
        }
        y
    }

    /// Scale columns by `d`: result `(i,j) = A(i,j) * d[j]`.
    pub fn mul_diag(&self, d: &[f64]) -> SparseMatrix {
        let mut out = SparseMatrix::new(self.n);
        for i in 0..self.n {
            for &(j, v) in &self.rows[i] {
                out.set(i, j, v * d[j]);
            }
        }
        out
    }
}

/// Residual norm `|Ax - b|₂`.
pub fn residual(a: &SparseMatrix, b: &[f64], x: &[f64]) -> f64 {
    let ax = a.matvec(x);
    ax.iter().zip(b).map(|(u, v)| (u - v).powi(2)).sum::<f64>().sqrt()
}

/// Conjugate gradient for a symmetric positive-definite `A`.
/// Returns the solution, or `Err` if it fails to converge.
pub fn conjugate_gradient(
    a: &SparseMatrix,
    b: &[f64],
    x0: &[f64],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<f64>, String> {
    if b.len() != a.n || x0.len() != a.n {
        return Err("conjugate_gradient: dimension mismatch".into());
    }
    let bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut x = x0.to_vec();
    let ax0 = a.matvec(&x);
    let mut r: Vec<f64> = b.iter().zip(&ax0).map(|(bi, ai)| bi - ai).collect();
    let mut p = r.clone();
    let mut rsold = r.iter().map(|v| v * v).sum::<f64>();

    for _ in 0..max_iter {
        if rsold.sqrt() <= tol * (1.0 + bnorm) {
            return Ok(x);
        }
        let ap = a.matvec(&p);
        let pap: f64 = p.iter().zip(&ap).map(|(pi, ai)| pi * ai).sum();
        if pap.abs() < 1e-300 {
            return Err("conjugate_gradient: breakdown (matrix not SPD?)".into());
        }
        let alpha = rsold / pap;
        for i in 0..a.n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rsnew = r.iter().map(|v| v * v).sum::<f64>();
        if rsnew.sqrt() <= tol * (1.0 + bnorm) {
            return Ok(x);
        }
        let beta = rsnew / rsold;
        for i in 0..a.n {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }
    Err(format!("conjugate_gradient: did not converge in {max_iter} iterations"))
}

/// Jacobi iteration (works well for diagonally-dominant matrices).
pub fn jacobi_iteration(
    a: &SparseMatrix,
    b: &[f64],
    x0: &[f64],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<f64>, String> {
    if b.len() != a.n || x0.len() != a.n {
        return Err("jacobi_iteration: dimension mismatch".into());
    }
    let bnorm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut x = x0.to_vec();
    for _ in 0..max_iter {
        let mut xnew = vec![0.0; a.n];
        for i in 0..a.n {
            let diag = a.get(i, i);
            if diag.abs() < 1e-300 {
                return Err("jacobi_iteration: zero diagonal".into());
            }
            let mut s = b[i];
            for &(j, v) in &a.rows[i] {
                if j != i {
                    s -= v * x[j];
                }
            }
            xnew[i] = s / diag;
        }
        x = xnew;
        if residual(a, b, &x) <= tol * (1.0 + bnorm) {
            return Ok(x);
        }
    }
    Err(format!("jacobi_iteration: did not converge in {max_iter} iterations"))
}

/// Spectral radius (largest |eigenvalue|) via power iteration.
pub fn spectral_radius(a: &SparseMatrix) -> f64 {
    let mut x = vec![1.0; a.n];
    let mut lambda = 0.0;
    for _ in 0..200 {
        let y = a.matvec(&x);
        // x is max-norm normalized (max |x_i| == 1), so ||y||∞ estimates |λ|.
        let new_lambda = y.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let yinf = new_lambda.max(1e-30);
        for i in 0..a.n {
            x[i] = y[i] / yinf;
        }
        if (new_lambda - lambda).abs() <= 1e-10 * (1.0 + new_lambda.abs()) {
            lambda = new_lambda;
            break;
        }
        lambda = new_lambda;
    }
    lambda
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::MathMatrix;

    /// 3×3 SPD tridiagonal: diag 2, off-diag 1.
    fn tri3() -> SparseMatrix {
        let mut a = SparseMatrix::new(3);
        for i in 0..3 {
            a.set(i, i, 2.0);
        }
        for i in 0..2 {
            a.set(i, i + 1, 1.0);
            a.set(i + 1, i, 1.0);
        }
        a
    }

    #[test]
    fn conjugate_gradient_solves_spd() {
        let a = tri3();
        assert!(a.is_symmetric(1e-12));
        // A * [1,2,3] = [4,8,8]
        let b = vec![4.0, 8.0, 8.0];
        let x = conjugate_gradient(&a, &b, &vec![0.0; 3], 100, 1e-10).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-8, "got {:?}", x);
        assert!((x[1] - 2.0).abs() < 1e-8);
        assert!((x[2] - 3.0).abs() < 1e-8);
        assert!(residual(&a, &b, &x) < 1e-6);
    }

    #[test]
    fn jacobi_solves_diagonally_dominant() {
        // [[4,1,0],[1,5,1],[0,1,4]], x=[1,2,3], b=[6,14,14]
        let mut a = SparseMatrix::new(3);
        a.set(0, 0, 4.0);
        a.set(0, 1, 1.0);
        a.set(1, 0, 1.0);
        a.set(1, 1, 5.0);
        a.set(1, 2, 1.0);
        a.set(2, 1, 1.0);
        a.set(2, 2, 4.0);
        let b = vec![6.0, 14.0, 14.0];
        let x = jacobi_iteration(&a, &b, &vec![0.0; 3], 1000, 1e-9).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-6, "got {:?}", x);
        assert!((x[1] - 2.0).abs() < 1e-6);
        assert!((x[2] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn from_dense_roundtrip_and_nnz() {
        let mut m = MathMatrix::new(1, 3, 1, 3);
        m.set_value(1, 1, 1.0);
        m.set_value(1, 2, 0.0);
        m.set_value(1, 3, 2.0);
        m.set_value(2, 1, 0.0);
        m.set_value(2, 2, 3.0);
        m.set_value(2, 3, 4.0);
        m.set_value(3, 1, 5.0);
        m.set_value(3, 2, 0.0);
        m.set_value(3, 3, 6.0);
        let s = SparseMatrix::from_dense(&m);
        assert_eq!(s.n, 3);
        assert_eq!(s.nnz(), 6);
        assert_eq!(s.get(0, 0), 1.0);
        assert_eq!(s.get(0, 1), 0.0);
        assert_eq!(s.get(2, 2), 6.0);
        assert!(!s.is_symmetric(1e-12));
        assert!((s.row_norm(1) - 5.0).abs() < 1e-12); // sqrt(3² + 4²)
    }

    #[test]
    fn matvec_and_residual_exact() {
        let mut a = SparseMatrix::new(2);
        a.set(0, 0, 2.0);
        a.set(0, 1, 1.0);
        a.set(1, 0, 1.0);
        a.set(1, 1, 3.0);
        let x = vec![1.0, 2.0];
        let y = a.matvec(&x);
        assert!((y[0] - 4.0).abs() < 1e-12);
        assert!((y[1] - 7.0).abs() < 1e-12);
        let b = vec![4.0, 7.0];
        assert!(residual(&a, &b, &x) < 1e-12);
        // mul_diag: scale columns by [2, 3]
        let d = a.mul_diag(&[2.0, 3.0]);
        assert_eq!(d.get(0, 0), 4.0);
        assert_eq!(d.get(0, 1), 3.0);
        assert_eq!(d.get(1, 1), 9.0);
    }

    #[test]
    fn spectral_radius_diag() {
        let mut a = SparseMatrix::new(3);
        a.set(0, 0, 4.0);
        a.set(1, 1, 2.0);
        a.set(2, 2, -1.0);
        let rho = spectral_radius(&a);
        assert!((rho - 4.0).abs() < 1e-6, "got {rho}");
    }
}
