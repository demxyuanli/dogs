//! Singular Value Decomposition. Source: `math_SVD.cxx`
use crate::{MathMatrix, MathVector, MathStatus};

/// SVD result: A = U * diag(W) * V^T
#[derive(Debug, Clone)]
pub struct SVD {
    pub u: MathMatrix,
    pub v: MathMatrix,
    pub w: MathVector,
    pub is_done: bool,
}

/// Sign function returning 1.0 or -1.0 (0 → 1.0)
fn sign(a: f64, b: f64) -> f64 { if b >= 0.0 { a.abs() } else { -a.abs() } }

impl SVD {
    /// Compute SVD: store U, V, W. A = U diag(W) V^T.
    /// Translated from C Numerical Recipes with 1-indexed array access.
    pub fn new(a: &MathMatrix) -> Self {
        let m = a.row_count();
        let n = a.col_count();
        let mp = m.max(n);
        let np = m.min(n);
        let mut u = a.clone();
        let mut v = MathMatrix::new(1, n, 1, n);
        let mut w = MathVector::new(1, np);
        let mut rv1 = MathVector::new(1, np);

        if mp < 2 { w.set_value(1, u.value(1,1)); return Self { u, v, w, is_done: true }; }

        // Householder reduction to bidiagonal form
        let mut g = 0.0; let mut scale = 0.0; let mut anorm = 0.0;
        for i in 1..=np {
            let l = i + 1;
            rv1.set_value(i, scale * g);
            g = 0.0; scale = 0.0;
            if i <= m {
                for k in i..=m { scale += u.value(k, i).abs(); }
                if scale != 0.0 {
                    for k in i..=m { u.set_value(k, i, u.value(k, i) / scale); let s = u.value(k, i); g += s * s; }
                    let mut f = u.value(i, i);
                    g = -sign(g.sqrt(), f);
                    let h = f * g - g * g;
                    u.set_value(i, i, f - g);
                    for j in l..=n {
                        let mut s = 0.0;
                        for k in i..=m { s += u.value(k, i) * u.value(k, j); }
                        f = s / h;
                        for k in i..=m { u.set_value(k, j, u.value(k, j) + f * u.value(k, i)); }
                    }
                    for k in i..=m { u.set_value(k, i, u.value(k, i) * scale); }
                }
            }
            w.set_value(i, scale * g);
            g = 0.0; scale = 0.0;
            if i <= m && i != n {
                for k in l..=n { scale += u.value(i, k).abs(); }
                if scale != 0.0 {
                    for k in l..=n { u.set_value(i, k, u.value(i, k) / scale); let s = u.value(i, k); g += s * s; }
                    let mut f = u.value(i, l);
                    g = -sign(g.sqrt(), f);
                    let h = f * g - f * f;
                    u.set_value(i, l, f - g);
                    for k in l..=n { rv1.set_value(k, u.value(i, k) / h); }
                    for j in l..=m {
                        let mut s = 0.0;
                        for k in l..=n { s += u.value(j, k) * u.value(i, k); }
                        for k in l..=n { u.set_value(j, k, u.value(j, k) + s * rv1.value(k)); }
                    }
                    for k in l..=n { u.set_value(i, k, u.value(i, k) * scale); }
                }
            }
            anorm = f64::max(anorm, w.value(i).abs() + rv1.value(i).abs());
        }

        // Accumulate V
        for i in (1..=np).rev() {
            let l = i + 1;
            if i < n {
                if g != 0.0 {
                    for j in l..=n { v.set_value(j, i, (u.value(i, j) / u.value(i, l)) / g); }
                    for j in l..=n {
                        let mut s = 0.0;
                        for k in l..=n { s += u.value(i, k) * v.value(k, j); }
                        for k in l..=n { v.set_value(k, j, v.value(k, j) + s * v.value(k, i)); }
                    }
                }
                for j in l..=n { v.set_value(i, j, 0.0); v.set_value(j, i, 0.0); }
            }
            v.set_value(i, i, 1.0);
            g = rv1.value(i);
        }

        // Accumulate U
        for i in (1..=np).rev() {
            let l = i + 1;
            g = w.value(i);
            for j in l..=np { u.set_value(i, j, 0.0); }
            if g != 0.0 {
                g = 1.0 / g;
                for j in l..=n {
                    let mut s = 0.0;
                    for k in l..=m { s += u.value(k, i) * u.value(k, j); }
                    let f = (s / u.value(i, i)) * g;
                    for k in i..=m { u.set_value(k, j, u.value(k, j) + f * u.value(k, i)); }
                }
                for j in i..=m { u.set_value(j, i, u.value(j, i) * g); }
            } else { for j in i..=m { u.set_value(j, i, 0.0); } }
            u.set_value(i, i, u.value(i, i) + 1.0);
        }

        // Diagonalization
        for k in (1..=np).rev() {
            for its in 1..=30 {
                let mut flag = true;
                let mut l = k;
                while l >= 2 {
                    let nm = l - 1;
                    if (rv1.value(l).abs() + anorm) == anorm { flag = false; break; }
                    if (w.value(nm).abs() + anorm) == anorm { break; }
                    l -= 1;
                }
                if flag {
                    let mut c = 0.0;
                    let mut sv = 1.0;
                    for i in l..=k {
                        let mut fs = sv * rv1.value(i);
                        rv1.set_value(i, c * rv1.value(i));
                        if (fs.abs() + anorm) == anorm { break; }
                        g = w.value(i);
                        let mut h = (fs * fs + g * g).sqrt();
                        w.set_value(i, h);
                        h = 1.0 / h;
                        c = g * h;
                        sv = -fs * h;
                        let nm_i = if i >= 2 { i - 1 } else { 1 };
                        for j in 1..=m {
                            let y = u.value(j, nm_i);
                            let z = u.value(j, i);
                            u.set_value(j, nm_i, y * c + z * sv);
                            u.set_value(j, i, z * c - y * sv);
                        }
                    }
                }
                let mut z = w.value(k);
                if l == k {
                    if z < 0.0 { w.set_value(k, -z); for j in 1..=n { v.set_value(j, k, -v.value(j, k)); } }
                    break;
                }
                if its == 30 { return Self { u, v, w, is_done: false }; }
                let mut x = w.value(l);
                let nm = k - 1;
                let mut y = w.value(nm);
                g = rv1.value(nm);
                let mut h = rv1.value(k);
                let mut f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                g = (f * f + 1.0).sqrt();
                f = ((x - z) * (x + z) + h * (y / (f + sign(g, f)) - h)) / x;
                let mut c = 1.0;
                let mut sv = 1.0;
                for j in l..=nm {
                    let i = j + 1;
                    g = rv1.value(i);
                    y = w.value(i);
                    h = sv * g;
                    g = c * g;
                    let mut zj = (f * f + h * h).sqrt();
                    rv1.set_value(j, zj);
                    c = f / zj;
                    sv = h / zj;
                    f = x * c + g * sv;
                    g = g * c - x * sv;
                    h = y * sv;
                    y *= c;
                    for jj in 1..=n {
                        let xv = v.value(jj, j);
                        let wv = v.value(jj, i);
                        v.set_value(jj, j, xv * c + wv * sv);
                        v.set_value(jj, i, wv * c - xv * sv);
                    }
                    zj = (f * f + h * h).sqrt();
                    w.set_value(j, zj);
                    if zj != 0.0 { c = f / zj; sv = h / zj; }
                    f = c * g + sv * y;
                    x = c * y - sv * g;
                    for jj in 1..=m {
                        let yv = u.value(jj, j);
                        let wv = u.value(jj, i);
                        u.set_value(jj, j, yv * c + wv * sv);
                        u.set_value(jj, i, wv * c - yv * sv);
                    }
                }
                rv1.set_value(l, 0.0);
                rv1.set_value(k, f);
                w.set_value(k, x);
            }
        }
        Self { u, v, w, is_done: true }
    }

    /// Solve Ax = b. x = V * diag(1/W) * U^T * b
    pub fn solve(&self, b: &MathVector) -> Result<MathVector, MathStatus> {
        let m = self.u.row_count();
        let n = self.v.row_count();
        let min_n = self.w.len();
        let mut tmp = MathVector::new(1, min_n);
        let wmax = (1..=min_n).fold(0.0f64, |mx, i| mx.max(self.w.value(i)));
        let wmin = wmax * 1e-15;

        for j in 1..=min_n {
            let mut s = 0.0;
            if self.w.value(j) > wmin {
                for i in 1..=m { s += self.u.value(i, j) * b.value(i) / self.w.value(j); }
            }
            tmp.set_value(j, s);
        }
        let mut x = MathVector::new(1, n);
        for j in 1..=n {
            let mut s = 0.0;
            for jj in 1..=min_n { s += self.v.value(j, jj) * tmp.value(jj); }
            x.set_value(j, s);
        }
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn svd_small_matrix() {
        let mut a = MathMatrix::new(1, 2, 1, 2);
        a.set_value(1,1,2.0); a.set_value(1,2,1.0);
        a.set_value(2,1,1.0); a.set_value(2,2,3.0);
        let svd = SVD::new(&a);
        assert!(svd.is_done);
        let b = MathVector::from_slice(&[5.0, 6.0]);
        let x = svd.solve(&b).unwrap();
        assert!((x.value(1) - 1.8).abs() < 1e-10);
        assert!((x.value(2) - 1.4).abs() < 1e-10);
    }
}
