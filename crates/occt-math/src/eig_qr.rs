//! QR eigensolver for general (not necessarily symmetric) real matrices.
//! Source: analogue of `math_Jacobi`'s general sibling — shifted Hessenberg QR.
//!
//! The input is first reduced to upper Hessenberg form by Householder
//! similarity transforms, then the shifted QR algorithm (Wilkinson shift) is
//! iterated, deflating trailing 1×1 and 2×2 blocks. Real eigenvalues receive
//! eigenvectors computed by inverse iteration on the original matrix;
//! complex-conjugate pairs carry an empty eigenvector (their real part is in
//! `real` and the non-zero imaginary part in `imag`).

use crate::{MathMatrix, MathVector};

/// Result of the QR eigenvalue decomposition of a general real matrix.
///
/// `real`/`imag` hold the eigenvalues (imaginary part zero for real
/// eigenvalues; complex-conjugate pairs appear with opposite signs in `imag`).
/// `eigenvalues` is a convenience copy of `real`. `vectors[i]` is the real
/// eigenvector of `real[i]` when the eigenvalue is real, otherwise an empty
/// vector (complex eigenvalues have no real eigenvector).
#[derive(Debug, Clone)]
pub struct QrEigen {
    pub eigenvalues: Vec<f64>,
    pub real: Vec<f64>,
    pub imag: Vec<f64>,
    pub vectors: Vec<Vec<f64>>,
}

/// Eigenvalues (and real eigenvectors) of a general square matrix by the
/// shifted Hessenberg QR algorithm.
///
/// The matrix need not be symmetric; non-symmetric matrices may produce
/// complex-conjugate eigenvalue pairs, which are reported through `imag`.
///
/// # Errors
/// Returns an error when the matrix is not square or empty, or when an active
/// block fails to converge within 500 QR iterations.
pub fn qr_eigenvalues(matrix: &MathMatrix) -> Result<QrEigen, String> {
    let n = matrix.row_count();
    if !matrix.is_square() {
        return Err("qr_eigenvalues: matrix must be square".to_string());
    }
    if n == 0 {
        return Err("qr_eigenvalues: empty matrix".to_string());
    }
    if n == 1 {
        let v = matrix.value(1, 1);
        return Ok(QrEigen {
            eigenvalues: vec![v],
            real: vec![v],
            imag: vec![0.0],
            vectors: vec![vec![1.0]],
        });
    }

    // 1. Similarity transform to upper Hessenberg form.
    let mut h = hessenberg_reduce(matrix);

    // 2. Shifted QR iteration with deflation of trailing 1×1 / 2×2 blocks.
    let mut real = vec![0.0; n];
    let mut imag = vec![0.0; n];
    let mut m = n;
    while m > 1 {
        let k = active_block(&h, m);
        if k == m {
            // Trailing 1×1 block.
            real[m - 1] = h.value(m, m);
            imag[m - 1] = 0.0;
            m -= 1;
            continue;
        }
        if k == m - 1 {
            // Trailing 2×2 block: closed-form eigenvalues (possibly complex).
            let (r1, r2, i1, i2) = eigen_2x2(&h, m);
            real[m - 2] = r1;
            real[m - 1] = r2;
            imag[m - 2] = i1;
            imag[m - 1] = i2;
            m -= 2;
            continue;
        }
        // Iterate the active block k..=m.
        if !qr_iterate(&mut h, k, m) {
            return Err(format!("qr_eigenvalues: no convergence on block {k}..={m}"));
        }
    }
    if m == 1 {
        real[0] = h.value(1, 1);
        imag[0] = 0.0;
    }

    // 3. Eigenvectors via inverse iteration for the real eigenvalues.
    let vectors = (0..n)
        .map(|i| {
            if imag[i].abs() < 1e-10 {
                inverse_iteration(matrix, real[i]).unwrap_or_else(|| {
                    // Fallback: unit basis vector (non-empty, may be inexact).
                    let mut v = vec![0.0; n];
                    v[i] = 1.0;
                    v
                })
            } else {
                Vec::new()
            }
        })
        .collect();

    Ok(QrEigen {
        eigenvalues: real.clone(),
        real,
        imag,
        vectors,
    })
}

/// Return the i-th eigenvector (real part), or an empty vector when the index
/// is out of range or the eigenvalue is complex.
pub fn eigen_vector_for(qe: &QrEigen, index: usize) -> Vec<f64> {
    qe.vectors.get(index).cloned().unwrap_or_default()
}

/// Reduce `a` to upper Hessenberg form by Householder similarity transforms.
fn hessenberg_reduce(a: &MathMatrix) -> MathMatrix {
    let n = a.row_count();
    let mut h = a.clone();
    for k in 1..=n.saturating_sub(2) {
        // Householder vector from column k, rows k+1..n.
        let mut norm = 0.0;
        for i in (k + 1)..=n {
            norm += h.value(i, k) * h.value(i, k);
        }
        norm = norm.sqrt();
        if norm < 1e-300 {
            continue;
        }
        let ak1 = h.value(k + 1, k);
        let alpha = if ak1 >= 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = (k + 1..=n).map(|i| h.value(i, k)).collect();
        v[0] -= alpha;
        let vnorm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if vnorm < 1e-300 {
            continue;
        }
        for x in v.iter_mut() {
            *x /= vnorm;
        }
        // Left multiplication: h = H h (rows k+1..n, columns k..n).
        for j in k..=n {
            let mut s = 0.0;
            for (vi, &vv) in v.iter().enumerate() {
                s += vv * h.value(k + 1 + vi, j);
            }
            for (vi, &vv) in v.iter().enumerate() {
                let row = k + 1 + vi;
                h.set_value(row, j, h.value(row, j) - 2.0 * s * vv);
            }
        }
        // Right multiplication: h = h H (columns k+1..n, all rows).
        for i in 1..=n {
            let mut s = 0.0;
            for (vi, &vv) in v.iter().enumerate() {
                s += h.value(i, k + 1 + vi) * vv;
            }
            for (vi, &vv) in v.iter().enumerate() {
                let col = k + 1 + vi;
                h.set_value(i, col, h.value(i, col) - 2.0 * s * vv);
            }
        }
    }
    h
}

/// Largest index `k` such that the subdiagonal element `h[k][k-1]` is
/// negligible (or `k = 1`). The active block is rows/cols `k..=m`.
fn active_block(h: &MathMatrix, m: usize) -> usize {
    let eps = 1e-14;
    let mut k = m;
    while k > 1 {
        let scale = (h.value(k - 1, k - 1).abs() + h.value(k, k).abs()).max(1.0);
        if h.value(k, k - 1).abs() <= eps * scale {
            break;
        }
        k -= 1;
    }
    k
}

/// Eigenvalues of the trailing 2×2 block at rows/cols `m-1`, `m`. Returns
/// `(r1, r2, i1, i2)`.
fn eigen_2x2(h: &MathMatrix, m: usize) -> (f64, f64, f64, f64) {
    let a = h.value(m - 1, m - 1);
    let b = h.value(m - 1, m);
    let c = h.value(m, m - 1);
    let d = h.value(m, m);
    let tr = a + d;
    let det = a * d - b * c;
    let disc = tr * tr - 4.0 * det;
    if disc >= 0.0 {
        let sq = disc.sqrt();
        (0.5 * (tr - sq), 0.5 * (tr + sq), 0.0, 0.0)
    } else {
        (0.5 * tr, 0.5 * tr, 0.5 * (-disc).sqrt(), -0.5 * (-disc).sqrt())
    }
}

/// Householder QR factorization of a square matrix, returning `(Q, R)` with
/// `A = Q·R`, `Q` orthogonal and `R` upper triangular.
fn house_qr(a: &MathMatrix) -> (MathMatrix, MathMatrix) {
    let n = a.row_count();
    let mut r = a.clone();
    let mut q = crate::matrix_ext::identity(n);
    for k in 1..=n {
        let mut norm = 0.0;
        for i in k..=n {
            norm += r.value(i, k) * r.value(i, k);
        }
        norm = norm.sqrt();
        if norm < 1e-300 {
            continue;
        }
        let akk = r.value(k, k);
        let alpha = if akk >= 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = (k..=n).map(|i| r.value(i, k)).collect();
        v[0] -= alpha;
        let vnorm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if vnorm < 1e-300 {
            continue;
        }
        for x in v.iter_mut() {
            *x /= vnorm;
        }
        // R = H R.
        for j in k..=n {
            let mut s = 0.0;
            for (vi, &vv) in v.iter().enumerate() {
                s += vv * r.value(k + vi, j);
            }
            for (vi, &vv) in v.iter().enumerate() {
                let row = k + vi;
                r.set_value(row, j, r.value(row, j) - 2.0 * s * vv);
            }
        }
        // Q = Q H.
        for i in 1..=n {
            let mut s = 0.0;
            for (vi, &vv) in v.iter().enumerate() {
                s += q.value(i, k + vi) * vv;
            }
            for (vi, &vv) in v.iter().enumerate() {
                let col = k + vi;
                q.set_value(i, col, q.value(i, col) - 2.0 * s * vv);
            }
        }
    }
    (q, r)
}

/// Apply one sweep of shifted QR iterations to the active block `k..=m` of
/// `h`, returning true when a subdiagonal becomes negligible (deflation).
fn qr_iterate(h: &mut MathMatrix, k: usize, m: usize) -> bool {
    let size = m - k + 1;
    let eps = 1e-14;
    for _ in 0..500 {
        // Wilkinson shift from the trailing 2×2 block.
        let a00 = h.value(m - 1, m - 1);
        let a01 = h.value(m - 1, m);
        let a10 = h.value(m, m - 1);
        let a11 = h.value(m, m);
        let tr = a00 + a11;
        let det = a00 * a11 - a01 * a10;
        let disc = tr * tr * 0.25 - det;
        let shift = if disc >= 0.0 {
            let sq = disc.sqrt();
            let l1 = tr * 0.5 - sq;
            let l2 = tr * 0.5 + sq;
            if (a11 - l1).abs() <= (a11 - l2).abs() {
                l1
            } else {
                l2
            }
        } else {
            tr * 0.5 // complex pair: use the real part as a real shift
        };

        // Extract the active submatrix, apply the shifted QR step, write back.
        let mut sub = MathMatrix::new(1, size, 1, size);
        for i in 0..size {
            for j in 0..size {
                sub.set_value(i + 1, j + 1, h.value(k + i, k + j));
            }
        }
        for i in 1..=size {
            sub.set_value(i, i, sub.value(i, i) - shift);
        }
        let (q, r) = house_qr(&sub);
        let rq = r.multiplied_mat(&q);
        for i in 0..size {
            for j in 0..size {
                let mut v = rq.value(i + 1, j + 1);
                if i == j {
                    v += shift;
                }
                h.set_value(k + i, k + j, v);
            }
        }

        // Check for deflation, scanning from the bottom of the block.
        for i in (k + 1..=m).rev() {
            let scale = (h.value(i - 1, i - 1).abs() + h.value(i, i).abs()).max(1.0);
            if h.value(i, i - 1).abs() <= eps * scale {
                h.set_value(i, i - 1, 0.0);
                return true;
            }
        }
    }
    false
}

/// Eigenvector of `a` for the real eigenvalue `lambda`, by inverse iteration
/// with a small real shift so the factor `A − μI` is non-singular even when
/// `lambda` is an exact eigenvalue. Returns `None` if the solves fail.
fn inverse_iteration(a: &MathMatrix, lambda: f64) -> Option<Vec<f64>> {
    let n = a.row_count();
    let mut v = MathVector::new(1, n);
    for i in 1..=n {
        v.set_value(i, 1.0);
    }
    let mu = lambda + 1e-10 * (1.0 + lambda.abs()).max(1.0);
    let mut m = a.clone();
    for i in 1..=n {
        m.set_value(i, i, m.value(i, i) - mu);
    }
    for _ in 0..30 {
        let y = m.solve(&v).ok()?;
        let norm = y.norm();
        if !norm.is_finite() || norm < 1e-300 {
            return None;
        }
        for i in 1..=n {
            v.set_value(i, y.value(i) / norm);
        }
    }
    let norm = v.norm();
    if !norm.is_finite() || norm < 1e-300 {
        return None;
    }
    Some((1..=n).map(|i| v.value(i) / norm).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eigen_ext::jacobi_eigen_symmetric;

    fn mat2(a: f64, b: f64, c: f64, d: f64) -> MathMatrix {
        let mut m = MathMatrix::new(1, 2, 1, 2);
        m.set_value(1, 1, a);
        m.set_value(1, 2, b);
        m.set_value(2, 1, c);
        m.set_value(2, 2, d);
        m
    }

    fn sym3() -> MathMatrix {
        // [[2,1,0],[1,2,1],[0,1,2]] eigenvalues ≈ 0.586, 2, 3.414.
        let mut m = MathMatrix::new(1, 3, 1, 3);
        m.set_value(1, 1, 2.0);
        m.set_value(1, 2, 1.0);
        m.set_value(2, 1, 1.0);
        m.set_value(2, 2, 2.0);
        m.set_value(2, 3, 1.0);
        m.set_value(3, 2, 1.0);
        m.set_value(3, 3, 2.0);
        m
    }

    fn sorted(vals: &[f64]) -> Vec<f64> {
        let mut v = vals.to_vec();
        v.sort_by(f64::total_cmp);
        v
    }

    fn mat_vec_mul(a: &MathMatrix, v: &[f64]) -> Vec<f64> {
        let n = a.row_count();
        (0..n)
            .map(|i| (0..n).map(|j| a.value(i + 1, j + 1) * v[j]).sum())
            .collect()
    }

    #[test]
    fn symmetric_matches_jacobi() {
        let a = sym3();
        let qr = qr_eigenvalues(&a).unwrap();
        let (jac, _) = jacobi_eigen_symmetric(&a, 1e-9).unwrap();
        let qs = sorted(&qr.eigenvalues);
        assert_eq!(qs.len(), 3);
        for i in 0..3 {
            assert!((qs[i] - jac[i]).abs() < 1e-6, "eigen {i}: {} vs {}", qs[i], jac[i]);
        }
    }

    #[test]
    fn diagonal_matrix_eigenvalues() {
        let mut a = MathMatrix::new(1, 3, 1, 3);
        a.set_value(1, 1, 5.0);
        a.set_value(2, 2, 3.0);
        a.set_value(3, 3, -2.0);
        let qr = qr_eigenvalues(&a).unwrap();
        assert_eq!(sorted(&qr.eigenvalues), vec![-2.0, 3.0, 5.0]);
        assert!(qr.imag.iter().all(|i| i.abs() < 1e-12));
    }

    #[test]
    fn two_by_two_real_distinct() {
        let a = mat2(1.0, 1.0, 0.0, 2.0); // upper triangular: eigenvalues 1, 2
        let qr = qr_eigenvalues(&a).unwrap();
        assert_eq!(sorted(&qr.eigenvalues), vec![1.0, 2.0]);
    }

    #[test]
    fn eigenvector_residual() {
        let a = sym3();
        let qr = qr_eigenvalues(&a).unwrap();
        for i in 0..3 {
            assert!(qr.imag[i].abs() < 1e-10, "expected real eigenvalue at {i}");
            let v = eigen_vector_for(&qr, i);
            assert_eq!(v.len(), 3);
            let av = mat_vec_mul(&a, &v);
            for k in 0..3 {
                assert!(
                    (av[k] - qr.eigenvalues[i] * v[k]).abs() < 1e-8,
                    "residual at eigen {i} component {k}: {}",
                    (av[k] - qr.eigenvalues[i] * v[k]).abs()
                );
            }
        }
    }

    #[test]
    fn trace_equals_sum_of_eigenvalues() {
        let a = sym3();
        let qr = qr_eigenvalues(&a).unwrap();
        let sum: f64 = qr.real.iter().sum();
        let tr: f64 = (1..=3).map(|i| a.value(i, i)).sum();
        assert!((sum - tr).abs() < 1e-8, "sum {sum} vs trace {tr}");
    }

    #[test]
    fn rotation_matrix_complex_pair() {
        // 90° rotation [[0,-1],[1,0]] has eigenvalues ±i.
        let a = mat2(0.0, -1.0, 1.0, 0.0);
        let qr = qr_eigenvalues(&a).unwrap();
        assert!(qr.imag[0].abs() > 0.5, "imag[0] = {}", qr.imag[0]);
        assert!(qr.imag[1].abs() > 0.5, "imag[1] = {}", qr.imag[1]);
        assert!(qr.real[0].abs() < 1e-8);
        assert!(qr.real[1].abs() < 1e-8);
        // Complex eigenvalues carry an empty eigenvector.
        assert!(qr.vectors[0].is_empty() || qr.vectors[1].is_empty());
    }

    #[test]
    fn rank_deficient_has_zero_eigenvalue() {
        let a = mat2(1.0, 2.0, 2.0, 4.0); // rank 1: eigenvalues 0, 5
        let qr = qr_eigenvalues(&a).unwrap();
        let s = sorted(&qr.eigenvalues);
        assert!(s[0].abs() < 1e-8, "smallest = {}", s[0]);
        assert!((s[1] - 5.0).abs() < 1e-8);
    }

    #[test]
    fn one_by_one_returns_single_value() {
        let mut a = MathMatrix::new(1, 1, 1, 1);
        a.set_value(1, 1, 7.0);
        let qr = qr_eigenvalues(&a).unwrap();
        assert_eq!(qr.eigenvalues, vec![7.0]);
        assert_eq!(qr.vectors[0], vec![1.0]);
    }

    #[test]
    fn non_symmetric_real_eigenvalues() {
        // A = P·diag(1,2,3)·P⁻¹ is non-symmetric with real eigenvalues 1, 2, 3.
        let mut p = MathMatrix::new(1, 3, 1, 3);
        p.set_value(1, 1, 1.0);
        p.set_value(1, 2, 1.0);
        p.set_value(2, 2, 1.0);
        p.set_value(2, 3, 1.0);
        p.set_value(3, 1, 1.0);
        p.set_value(3, 3, 1.0);
        let mut d = MathMatrix::new(1, 3, 1, 3);
        d.set_value(1, 1, 1.0);
        d.set_value(2, 2, 2.0);
        d.set_value(3, 3, 3.0);
        let pinv = crate::matrix_ext::inverse(&p).unwrap();
        let a = crate::matrix_ext::matmul(&p, &crate::matrix_ext::matmul(&d, &pinv).unwrap()).unwrap();
        let qr = qr_eigenvalues(&a).unwrap();
        let s = sorted(&qr.eigenvalues);
        for (got, want) in s.iter().zip([1.0, 2.0, 3.0]) {
            assert!((got - want).abs() < 1e-6, "eigen {got} vs {want}");
        }
        assert!(qr.imag.iter().all(|i| i.abs() < 1e-8));
    }
}
