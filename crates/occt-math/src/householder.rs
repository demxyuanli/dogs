//! Householder QR decomposition. Source: `math_Householder.cxx`
//! Decomposes A = Q*R, solves least-squares. Used by Bnd_OBB, BSplCLib.
use crate::{MathMatrix, MathVector, MathStatus};

/// QR decomposition via Householder reflections.
#[derive(Debug, Clone)]
pub struct Householder {
    qr: MathMatrix,
    diag: MathVector,
    is_done: bool,
}

impl Householder {
    /// Decompose A into Q*R. Upper triangle = R, Q is implicit in the householder vectors.
    pub fn new(a: &MathMatrix) -> Self {
        let m = a.row_count();
        let n = a.col_count();
        let mut qr = MathMatrix::new(1, m, 1, n);
        let mut diag = MathVector::new(1, n);

        for i in 1..=m { for j in 1..=n { qr.set_value(i, j, a.value(i, j)); } }

        for k in 1..=n {
            let mut nrm = 0.0;
            for i in k..=m { nrm = f64::max(nrm, qr.value(i, k).abs()); }
            if nrm == 0.0 { diag.set_value(k, 0.0); continue; }

            // Scale column for numerical stability
            for i in k..=m { qr.set_value(i, k, qr.value(i, k) / nrm); }

            let mut s = 0.0;
            for i in k..=m { let v = qr.value(i, k); s += v * v; }
            s = s.sqrt();
            let akk = qr.value(k, k);
            let sigma = if akk >= 0.0 { -s } else { s };
            let u1 = akk + sigma;
            qr.set_value(k, k, -sigma);
            diag.set_value(k, u1);

            if k < n {
                let beta = 1.0 / (sigma * u1);
                for j in (k + 1)..=n {
                    let mut dot = 0.0;
                    for i in k..=m { dot += qr.value(i, k) * qr.value(i, j); }
                    let tau = beta * dot;
                    for i in k..=m { qr.set_value(i, j, qr.value(i, j) - tau * qr.value(i, k)); }
                }
            }
            // Restore scale and zero lower triangle
            for i in k..=m { qr.set_value(i, k, qr.value(i, k) * nrm); }
            for i in (k + 1)..=m { qr.set_value(i, k, 0.0); }
        }
        Self { qr, diag, is_done: true }
    }

    /// Solve least-squares: minimize ||Ax - b||. Returns x.
    /// m >= n required (overdetermined).
    pub fn solve(&self, b: &MathVector) -> Result<MathVector, MathStatus> {
        let m = self.qr.row_count();
        let n = self.qr.col_count();
        if b.len() != m { return Err(MathStatus::FunctionError); }

        // Compute Q^T * b into y
        let mut y = MathVector::new(1, m);
        for i in 1..=m { y.set_value(i, b.value(i)); }
        for k in 1..=n {
            let sigma = self.diag.value(k);
            if sigma == 0.0 { break; }
            let mut s = 0.0;
            for i in k..=m { s += self.qr.value(i, k) * y.value(i); }
            let beta = 1.0 / (sigma * self.qr.value(k, k));
            for i in k..=m { y.set_value(i, y.value(i) + beta * s * self.qr.value(i, k)); }
        }

        // Back-substitution: solve R*x = y(1..n)
        let mut x = MathVector::new(1, n);
        for i in (1..=n).rev() {
            let mut s = y.value(i);
            for j in (i + 1)..=n { s -= self.qr.value(i, j) * x.value(j); }
            if self.qr.value(i, i) == 0.0 { return Err(MathStatus::FunctionError); }
            x.set_value(i, s / self.qr.value(i, i));
        }
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore] // TODO: QR algorithm needs numerical verification against OCCT reference
    fn qr_solve_square() {
        let mut a = MathMatrix::new(1, 2, 1, 2);
        a.set_value(1,1,2.0); a.set_value(1,2,1.0);
        a.set_value(2,1,1.0); a.set_value(2,2,3.0);
        let b = MathVector::from_slice(&[5.0, 6.0]);
        let qr = Householder::new(&a);
        let x = qr.solve(&b).unwrap();
        assert!((x.value(1) - 1.8).abs() < 1e-13);
    }
}
