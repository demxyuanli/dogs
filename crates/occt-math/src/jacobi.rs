//! Jacobi eigenvalue decomposition for symmetric matrices.
//! Source: `math_Jacobi.cxx`
//! Finds eigenvalues and eigenvectors of a real symmetric matrix.
use crate::{MathMatrix, MathVector};

/// Result of Jacobi eigenvalue decomposition.
#[derive(Debug, Clone)]
pub struct Jacobi {
    pub values: MathVector,   // eigenvalues in descending order
    pub vectors: MathMatrix,  // eigenvectors as columns
    pub is_done: bool,
}

impl Jacobi {
    /// Compute eigenvalues and eigenvectors for symmetric matrix A.
    /// Uses Jacobi iteration with Givens rotations.
    pub fn new(a: &MathMatrix) -> Self {
        let n = a.row_count();
        let mut eigv = MathVector::new(1, n);
        let mut a_mat = a.clone();
        let mut v = MathMatrix::new(1, n, 1, n);
        for i in 1..=n { v.set_value(i, i, 1.0); }

        // Initialize eigenvalues to diagonal
        for i in 1..=n {
            eigv.set_value(i, a.value(i, i));
        }
        let mut b = MathVector::new(1, n);
        let mut z = MathVector::new(1, n);
        for i in 1..=n { b.set_value(i, eigv.value(i)); z.set_value(i, 0.0); }

        let max_iter = 50;
        for iter in 1..=max_iter {
            // Compute sum of off-diagonal elements
            let mut sm = 0.0;
            for i in 1..n {
                for j in (i + 1)..=n { sm += a_mat.value(i, j).abs(); }
            }
            if sm == 0.0 { break; }

            let threshold = if iter < 4 { 0.2 * sm / (n * n) as f64 } else { 0.0 };

            for p in 1..n {
                for q in (p + 1)..=n {
                    let g = 100.0 * a_mat.value(p, q).abs();
                    if iter > 4 && eigv.value(p).abs() + g == eigv.value(p).abs()
                        && eigv.value(q).abs() + g == eigv.value(q).abs() {
                        a_mat.set_value(p, q, 0.0);
                    } else if a_mat.value(p, q).abs() > threshold {
                        let h = eigv.value(q) - eigv.value(p);
                        let t = if h.abs() + g == h.abs() {
                            a_mat.value(p, q) / h
                        } else {
                            let theta = 0.5 * h / a_mat.value(p, q);
                            let mut t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                            if theta < 0.0 { t = -t; }
                            t
                        };
                        let c = 1.0 / (1.0 + t * t).sqrt();
                        let s = t * c;
                        let tau = s / (1.0 + c);
                        let h = t * a_mat.value(p, q);
                        z.set_value(p, z.value(p) - h);
                        z.set_value(q, z.value(q) + h);
                        eigv.set_value(p, eigv.value(p) - h);
                        eigv.set_value(q, eigv.value(q) + h);
                        a_mat.set_value(p, q, 0.0);
                        for j in 1..p {
                            let g = a_mat.value(j, p);
                            let h = a_mat.value(j, q);
                            a_mat.set_value(j, p, g - s * (h + g * tau));
                            a_mat.set_value(j, q, h + s * (g - h * tau));
                        }
                        for j in (p + 1)..q {
                            let g = a_mat.value(p, j);
                            let h = a_mat.value(j, q);
                            a_mat.set_value(p, j, g - s * (h + g * tau));
                            a_mat.set_value(j, q, h + s * (g - h * tau));
                        }
                        for j in (q + 1)..=n {
                            let g = a_mat.value(p, j);
                            let h = a_mat.value(q, j);
                            a_mat.set_value(p, j, g - s * (h + g * tau));
                            a_mat.set_value(q, j, h + s * (g - h * tau));
                        }
                        for j in 1..=n {
                            let g = v.value(j, p);
                            let h = v.value(j, q);
                            v.set_value(j, p, g - s * (h + g * tau));
                            v.set_value(j, q, h + s * (g - h * tau));
                        }
                    }
                }
            }
            for i in 1..=n {
                b.set_value(i, b.value(i) + z.value(i));
                eigv.set_value(i, b.value(i));
                z.set_value(i, 0.0);
            }
        }

        // Sort eigenvalues descending, permute eigenvectors
        for i in 1..n {
            let mut k = i;
            let mut p = eigv.value(i);
            for j in (i + 1)..=n {
                if eigv.value(j) > p { k = j; p = eigv.value(j); }
            }
            if k != i {
                eigv.set_value(k, eigv.value(i));
                eigv.set_value(i, p);
                for j in 1..=n {
                    let tmp = v.value(j, i);
                    v.set_value(j, i, v.value(j, k));
                    v.set_value(j, k, tmp);
                }
            }
        }

        Self { values: eigv, vectors: v, is_done: true }
    }

    /// Compute OBB axes from point cloud. Returns (center, x_dir, y_dir, z_dir, hx, hy, hz).
    pub fn compute_obb_axes(points: &[MathVector]) -> (MathVector, MathVector, MathVector, MathVector, f64, f64, f64) {
        let n_pts = points.len();
        let n = points[0].len();
        // Build covariance matrix
        let mut cov = MathMatrix::new(1, n, 1, n);
        // Compute centroid
        let mut centroid = MathVector::new(1, n);
        for p in points { centroid.add(p); }
        for i in 1..=n { centroid.set_value(i, centroid.value(i) / n_pts as f64); }
        // Covariance
        for p in points {
            for i in 1..=n {
                let di = p.value(i) - centroid.value(i);
                for j in 1..=n {
                    cov.set_value(i, j, cov.value(i, j) + di * (p.value(j) - centroid.value(j)));
                }
            }
        }
        let jacobi = Jacobi::new(&cov);
        let mut ax0 = MathVector::new(1,n);
        let mut ax1 = MathVector::new(1,n);
        let mut ax2 = MathVector::new(1,n);
        for k in 1..=n {
            ax0.set_value(k, jacobi.vectors.value(k, 1));
            ax1.set_value(k, jacobi.vectors.value(k, 2));
            ax2.set_value(k, jacobi.vectors.value(k, 3));
        }
        let mut hdims = [0.0f64; 3];
        let axes_arr = [&ax0, &ax1, &ax2];
        for (ai, ax) in axes_arr.iter().enumerate() {
            let mut mn = f64::MAX; let mut mx = f64::MIN;
            for p in points {
                let proj = ax.dot(p);
                mn = f64::min(mn, proj); mx = f64::max(mx, proj);
            }
            hdims[ai] = (mx - mn) * 0.5;
        }
        (centroid, ax0, ax1, ax2, hdims[0], hdims[1], hdims[2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jacobi_diagonal() {
        let mut a = MathMatrix::new(1, 2, 1, 2);
        a.set_value(1,1,5.0); a.set_value(1,2,0.0);
        a.set_value(2,1,0.0); a.set_value(2,2,3.0);
        let j = Jacobi::new(&a);
        assert!(j.is_done);
        let evs = [j.values.value(1), j.values.value(2)];
        assert!((evs.iter().max_by(|a,b| a.partial_cmp(b).unwrap()).unwrap() - 5.0).abs() < 1e-14);
    }
}
