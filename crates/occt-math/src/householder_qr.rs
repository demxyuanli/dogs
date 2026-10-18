//! Complete Householder QR decomposition with an explicit Q matrix.
//! Source: `math_Householder`.

use crate::{MathMatrix, MathStatus, MathVector};

/// Full Householder QR result: `A = Q * R`, Q is m x m orthogonal, R is m x n upper triangular.
pub struct HouseholderFull {
    pub q: MathMatrix,
    pub r: MathMatrix,
    pub is_done: bool,
}

/// Householder QR decomposition of `a` (m x n).
///
/// For each column `k` a reflector `H_k` is built from rows `k..=m`, applied to the
/// remaining columns of `R`, and accumulated into `Q`. Zero-norm columns produce an
/// identity reflector.
pub fn qr_decompose(a: &MathMatrix) -> HouseholderFull {
    let m = a.row_count();
    let n = a.col_count();

    // R starts as a copy of A.
    let mut r = MathMatrix::new(1, m, 1, n);
    for i in 1..=m {
        for j in 1..=n {
            r.set_value(i, j, a.value(i, j));
        }
    }

    // Q accumulates the reflectors, starting from identity.
    let mut q = MathMatrix::new(1, m, 1, m);
    for i in 1..=m {
        q.set_value(i, i, 1.0);
    }

    for k in 1..=n.min(m) {
        let len = m - k + 1;

        // Householder vector v for column k (rows k..=m).
        let mut v = vec![0.0; len];
        let mut norm2 = 0.0;
        for i in 0..len {
            let val = r.value(k + i, k);
            v[i] = val;
            norm2 += val * val;
        }
        if norm2 <= f64::EPSILON {
            continue; // zero column -> identity reflector
        }
        let norm = norm2.sqrt();
        let alpha = if v[0] >= 0.0 { -norm } else { norm };
        v[0] -= alpha;
        let vnorm2: f64 = v.iter().map(|x| x * x).sum();
        if vnorm2 <= f64::EPSILON {
            continue;
        }
        let beta = 2.0 / vnorm2;

        // R = H_k R: apply the reflector to columns k..=n.
        for j in k..=n {
            let mut dot = 0.0;
            for i in 0..len {
                dot += v[i] * r.value(k + i, j);
            }
            let scale = beta * dot;
            for i in 0..len {
                let val = r.value(k + i, j) - scale * v[i];
                r.set_value(k + i, j, val);
            }
        }

        // Q = Q * H_k.
        for i in 1..=m {
            let mut dot = 0.0;
            for t in 0..len {
                dot += q.value(i, k + t) * v[t];
            }
            let scale = beta * dot;
            for t in 0..len {
                let val = q.value(i, k + t) - scale * v[t];
                q.set_value(i, k + t, val);
            }
        }
    }

    HouseholderFull { q, r, is_done: true }
}

/// Minimize `||A x - b||` via QR: `y = Q^T b`, then back-solve `R x = y` (top n rows).
/// Errors with [`MathStatus::FunctionError`] on a near-zero R diagonal (rank deficient).
pub fn solve_least_squares(hf: &HouseholderFull, b: &MathVector) -> Result<MathVector, MathStatus> {
    let m = hf.q.row_count();
    let n = hf.r.col_count();
    if b.len() != m {
        return Err(MathStatus::FunctionError);
    }

    // y = Q^T b.
    let mut y = MathVector::new(1, m);
    for i in 1..=m {
        let mut s = 0.0;
        for j in 1..=m {
            s += hf.q.value(j, i) * b.value(j);
        }
        y.set_value(i, s);
    }

    // Back-substitute the top n rows of R.
    let mut x = MathVector::new(1, n);
    for i in (1..=n).rev() {
        let mut s = y.value(i);
        for j in (i + 1)..=n {
            s -= hf.r.value(i, j) * x.value(j);
        }
        let diag = hf.r.value(i, i);
        if diag.abs() < 1e-12 {
            return Err(MathStatus::FunctionError);
        }
        x.set_value(i, s / diag);
    }
    Ok(x)
}

/// Sign of `det(Q)` (Q is orthogonal, so `±1`), via Gaussian elimination with pivoting.
fn det_sign(q: &MathMatrix) -> f64 {
    let n = q.row_count();
    let mut a = MathMatrix::new(1, n, 1, n);
    for i in 1..=n {
        for j in 1..=n {
            a.set_value(i, j, q.value(i, j));
        }
    }
    let mut sign = 1.0;
    for k in 1..=n {
        let mut piv = k;
        let mut best = a.value(k, k).abs();
        for i in (k + 1)..=n {
            let val = a.value(i, k).abs();
            if val > best {
                best = val;
                piv = i;
            }
        }
        if best < 1e-15 {
            return 0.0;
        }
        if piv != k {
            sign = -sign;
            for j in k..=n {
                let t = a.value(k, j);
                a.set_value(k, j, a.value(piv, j));
                a.set_value(piv, j, t);
            }
        }
        let d = a.value(k, k);
        for i in (k + 1)..=n {
            let f = a.value(i, k) / d;
            for j in (k + 1)..=n {
                let val = a.value(i, j) - f * a.value(k, j);
                a.set_value(i, j, val);
            }
        }
    }
    sign
}

/// Determinant from the QR factors: `det(A) = det(Q) * det(R)`.
pub fn determinant_from_r(hf: &HouseholderFull) -> f64 {
    let n = hf.r.row_count().min(hf.r.col_count());
    let mut det = det_sign(&hf.q);
    for i in 1..=n {
        det *= hf.r.value(i, i);
    }
    det
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn least_squares_line() {
        // Fit y = 2x + 1 from the three exact points.
        let mut a = MathMatrix::new(1, 3, 1, 2);
        a.set_value(1, 1, 0.0);
        a.set_value(1, 2, 1.0);
        a.set_value(2, 1, 1.0);
        a.set_value(2, 2, 1.0);
        a.set_value(3, 1, 2.0);
        a.set_value(3, 2, 1.0);
        let mut b = MathVector::new(1, 3);
        b.set_value(1, 1.0);
        b.set_value(2, 3.0);
        b.set_value(3, 5.0);

        let hf = qr_decompose(&a);
        let x = solve_least_squares(&hf, &b).unwrap();
        assert!((x.value(1) - 2.0).abs() < 1e-8, "slope {}", x.value(1));
        assert!((x.value(2) - 1.0).abs() < 1e-8, "intercept {}", x.value(2));
    }
}
