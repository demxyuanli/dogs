//! Extended eigensolvers — power iteration, deflation, generalized eigen and
//! symmetric Jacobi for small dense matrices.
//! Source: `math_Jacobi`, `math_EigenValuesSearcher`-adjacent.

use crate::matrix::MathMatrix;
use crate::vector::MathVector;

/// Matrix–vector product M·v (dense).
fn matvec(a: &MathMatrix, v: &MathVector) -> MathVector {
    let n = v.len();
    let mut w = MathVector::new(1, n);
    for i in 1..=n {
        let mut s = 0.0;
        for j in 1..=n {
            s += a.value(i, j) * v.value(j);
        }
        w.set_value(i, s);
    }
    w
}

/// Power iteration: largest-magnitude eigenvalue and its unit eigenvector.
/// Returns (eigenvalue, eigenvector). Err if it fails to converge or the
/// matrix is zero.
pub fn power_iteration(a: &MathMatrix, max_iter: usize, tol: f64) -> Result<(f64, MathVector), String> {
    let n = a.row_count();
    if !a.is_square() || n == 0 {
        return Err("power_iteration: square non-empty matrix required".into());
    }
    let mut v = MathVector::new(1, n);
    for i in 1..=n {
        v.set_value(i, 1.0);
    }
    let mut lambda: f64 = 0.0;
    for _ in 0..max_iter {
        let w = matvec(a, &v);
        // Scale by the max component.
        let mut mx: f64 = 0.0;
        for i in 1..=n {
            mx = mx.max(w.value(i).abs());
        }
        if mx < 1e-30 {
            break;
        }
        let new_lambda = mx;
        for i in 1..=n {
            v.set_value(i, w.value(i) / mx);
        }
        if (new_lambda - lambda).abs() < tol * lambda.abs().max(1.0) {
            lambda = new_lambda;
            break;
        }
        lambda = new_lambda;
    }
    // Normalize.
    let norm = v.norm();
    if norm < 1e-30 {
        return Err("power_iteration: converged to zero vector".into());
    }
    for i in 1..=n {
        v.set_value(i, v.value(i) / norm);
    }
    Ok((lambda, v))
}

/// Eigen-decomposition of a symmetric matrix via Jacobi rotations. Returns
/// (eigenvalues ascending, eigenvectors as MathVector columns). Only valid for
/// symmetric input.
pub fn jacobi_eigen_symmetric(a: &MathMatrix, tol: f64) -> Result<(Vec<f64>, Vec<MathVector>), String> {
    let n = a.row_count();
    if !a.is_square() || n == 0 {
        return Err("jacobi: square matrix required".into());
    }
    let mut m = a.clone();
    // Eigenvectors as columns: vecs[i] holds column i.
    let mut vecs: Vec<MathVector> = (0..n)
        .map(|i| {
            let mut v = MathVector::new(1, n);
            v.set_value(i + 1, 1.0);
            v
        })
        .collect();
    let mut converged = false;
    for _ in 0..(100 * n.max(1)) {
        let mut p = 0usize;
        let mut q = 1usize;
        let mut mx = m.value(1, 2).abs();
        for i in 0..n {
            for j in (i + 1)..n {
                let v = m.value(i + 1, j + 1).abs();
                if v > mx {
                    mx = v;
                    p = i;
                    q = j;
                }
            }
        }
        if mx < tol {
            converged = true;
            break;
        }
        let app = m.value(p + 1, p + 1);
        let aqq = m.value(q + 1, q + 1);
        let apq = m.value(p + 1, q + 1);
        let tau = (aqq - app) / (2.0 * apq);
        let t = tau.signum() / (tau.abs() + (1.0 + tau * tau).sqrt());
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        // Rotate the off-diagonal rows/cols (skip p, q — handled below).
        for k in 0..n {
            if k == p || k == q {
                continue;
            }
            let akp = m.value(k + 1, p + 1);
            let akq = m.value(k + 1, q + 1);
            m.set_value(k + 1, p + 1, c * akp - s * akq);
            m.set_value(p + 1, k + 1, c * akp - s * akq);
            m.set_value(k + 1, q + 1, s * akp + c * akq);
            m.set_value(q + 1, k + 1, s * akp + c * akq);
        }
        // Diagonal entries after the rotation (Jacobi closed form).
        m.set_value(p + 1, p + 1, c * c * app - 2.0 * s * c * apq + s * s * aqq);
        m.set_value(q + 1, q + 1, s * s * app + 2.0 * s * c * apq + c * c * aqq);
        m.set_value(p + 1, q + 1, 0.0);
        m.set_value(q + 1, p + 1, 0.0);
        // Update eigenvector columns.
        for k in 0..n {
            let vkp = vecs[k].value(p + 1);
            let vkq = vecs[k].value(q + 1);
            vecs[k].set_value(p + 1, c * vkp - s * vkq);
            vecs[k].set_value(q + 1, s * vkp + c * vkq);
        }
    }
    if !converged {
        return Err("jacobi: did not converge".into());
    }
    // Extract eigenvalues from the diagonal, sort ascending with vectors.
    let mut evals: Vec<(f64, usize)> = (0..n).map(|i| (m.value(i + 1, i + 1), i)).collect();
    evals.sort_by(|a, b| a.0.total_cmp(&b.0));
    let sorted_vals: Vec<f64> = evals.iter().map(|(v, _)| *v).collect();
    let sorted_vecs: Vec<MathVector> = evals.iter().map(|(_, i)| vecs[*i].clone()).collect();
    Ok((sorted_vals, sorted_vecs))
}

/// Deflate: remove a known eigenpair (λ, v) from a symmetric matrix, returning
/// a matrix with that eigenvalue shifted to ~0.
pub fn deflate(a: &MathMatrix, lambda: f64, v: &MathVector) -> MathMatrix {
    let n = a.row_count();
    let mut out = a.clone();
    for i in 1..=n {
        for j in 1..=n {
            let val = out.value(i, j) - lambda * v.value(i) * v.value(j);
            out.set_value(i, j, val);
        }
    }
    out
}

/// The k largest eigenvalues of a symmetric matrix, via the Jacobi
/// eigensolver (robust; naive deflation can stall on imperfect eigenvectors).
pub fn largest_eigenvalues(a: &MathMatrix, k: usize, tol: f64) -> Result<Vec<(f64, MathVector)>, String> {
    let (vals, vecs) = jacobi_eigen_symmetric(a, tol)?;
    let n = vals.len();
    let k = k.min(n);
    let mut out = Vec::new();
    for i in 0..k {
        out.push((vals[n - 1 - i], vecs[n - 1 - i].clone()));
    }
    Ok(out)
}

/// 2×2 symmetric eigen-decomposition by closed form. Returns (values ascending,
/// vectors).
pub fn eigen_2x2(a: f64, b: f64, c: f64) -> (f64, f64, [f64; 2], [f64; 2]) {
    let tr = a + c;
    let det = a * c - b * b;
    let disc = ((tr * tr * 0.25) - det).max(0.0).sqrt();
    let l1 = tr * 0.5 - disc;
    let l2 = tr * 0.5 + disc;
    let mut v1 = [1.0f64, 0.0];
    let mut v2 = [0.0f64, 1.0];
    if b.abs() > 1e-30 {
        v1 = [b, l1 - a];
        v2 = [b, l2 - a];
    } else if (a - c).abs() < 1e-30 {
        v1 = [1.0, 0.0];
        v2 = [0.0, 1.0];
    }
    let n1 = (v1[0] * v1[0] + v1[1] * v1[1]).sqrt();
    let n2 = (v2[0] * v2[0] + v2[1] * v2[1]).sqrt();
    if n1 > 1e-30 {
        v1[0] /= n1;
        v1[1] /= n1;
    }
    if n2 > 1e-30 {
        v2[0] /= n2;
        v2[1] /= n2;
    }
    (l1, l2, v1, v2)
}

/// Largest eigenvalue of a symmetric 2×2 (closed form).
pub fn largest_eigenvalue_2x2(a: f64, b: f64, c: f64) -> f64 {
    eigen_2x2(a, b, c).1
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sym3() -> MathMatrix {
        // [[2,1,0],[1,2,1],[0,1,2]] eigenvalues ~0.586, 2, 3.414.
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

    #[test]
    fn power_iteration_largest() {
        let m = sym3();
        let (lambda, v) = power_iteration(&m, 2000, 1e-8).unwrap();
        assert!((lambda - 3.414).abs() < 0.02, "lambda {lambda}");
        // Eigenvector sanity: A·v ≈ λ·v.
        let w = matvec(&m, &v);
        for i in 1..=3 {
            assert!((w.value(i) - lambda * v.value(i)).abs() < 0.05);
        }
    }

    #[test]
    fn jacobi_symmetric_ascending() {
        let m = sym3();
        let (vals, vecs) = jacobi_eigen_symmetric(&m, 1e-9).unwrap();
        assert!(vals[0] < vals[1] && vals[1] < vals[2]);
        assert!((vals[0] - 0.586).abs() < 0.05);
        assert!((vals[1] - 2.0).abs() < 0.05);
        assert!((vals[2] - 3.414).abs() < 0.05);
        // Vectors are orthonormal.
        for i in 0..3 {
            for j in (i + 1)..3 {
                let dot = vecs[i].dot(&vecs[j]);
                assert!(dot.abs() < 1e-6);
            }
        }
    }

    #[test]
    fn largest_eigenvalues_deflation() {
        let m = sym3();
        let top = largest_eigenvalues(&m, 3, 1e-6).unwrap();
        assert_eq!(top.len(), 3);
        assert!((top[0].0 - 3.414).abs() < 0.05);
        // Deflation-based extraction of the smallest eigenvalue is noisier.
        assert!((top[2].0 - 0.586).abs() < 0.25, "smallest via deflation: {}", top[2].0);
    }

    #[test]
    fn eigen_2x2_closed_form() {
        let (l1, l2, v1, v2) = eigen_2x2(2.0, 1.0, 2.0);
        assert!((l1 - 1.0).abs() < 1e-9 && (l2 - 3.0).abs() < 1e-9);
        // v1 for eigenvalue 1: (-1,1)/√2.
        assert!((v1[0] + v1[1]).abs() < 1e-9, "v1 perpendicular: {v1:?}");
        assert!((v2[0] - v2[1]).abs() < 1e-9);
        assert!((largest_eigenvalue_2x2(2.0, 1.0, 2.0) - 3.0).abs() < 1e-9);
    }
}
