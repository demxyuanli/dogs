//! LU decomposition (Crout algorithm). Source: `math_Crout.cxx`
//! Decomposes A = L * U where L is lower-triangular, U is upper-triangular.
use crate::{MathMatrix, MathVector, MathStatus};

/// LU decomposition via Crout's method with implicit scaling.
#[derive(Debug, Clone)]
pub struct Crout {
    lu: MathMatrix,
    indx: Vec<usize>,
    det_sign: f64,
    is_done: bool,
}

impl Crout {
    /// Decompose matrix A. Returns Err if singular.
    pub fn new(a: &MathMatrix) -> Result<Self, MathStatus> {
        let n = a.row_count();
        if n != a.col_count() { return Err(MathStatus::FunctionError); }
        let mut lu = a.clone();
        let mut indx = vec![0usize; n + 1]; // 1-indexed
        let mut vv = vec![0.0f64; n + 1];
        let mut det_sign = 1.0;

        // Compute implicit scaling for each row
        for i in 1..=n {
            let mut big = 0.0;
            for j in 1..=n { big = f64::max(big, lu.value(i, j).abs()); }
            if big == 0.0 { return Err(MathStatus::FunctionError); }
            vv[i] = 1.0 / big;
        }

        for j in 1..=n {
            for i in 1..j {
                let mut sum = lu.value(i, j);
                for k in 1..i { sum -= lu.value(i, k) * lu.value(k, j); }
                lu.set_value(i, j, sum);
            }
            let mut big = 0.0;
            let mut imax = j;
            for i in j..=n {
                let mut sum = lu.value(i, j);
                for k in 1..j { sum -= lu.value(i, k) * lu.value(k, j); }
                lu.set_value(i, j, sum);
                let dum = vv[i] * sum.abs();
                if dum >= big { big = dum; imax = i; }
            }
            if j != imax {
                for k in 1..=n {
                    let tmp = lu.value(imax, k);
                    lu.set_value(imax, k, lu.value(j, k));
                    lu.set_value(j, k, tmp);
                }
                det_sign = -det_sign;
                vv[imax] = vv[j];
            }
            indx[j] = imax;
            if lu.value(j, j) == 0.0 { lu.set_value(j, j, 1e-30); }
            if j != n {
                let pivot = 1.0 / lu.value(j, j);
                for i in (j + 1)..=n { lu.set_value(i, j, lu.value(i, j) * pivot); }
            }
        }
        Ok(Self { lu, indx, det_sign, is_done: true })
    }

    /// Solve Ax = b using the stored LU decomposition.
    pub fn solve(&self, b: &MathVector) -> Result<MathVector, MathStatus> {
        let n = self.lu.row_count();
        let mut x = b.clone();

        // Forward substitution: solve L*y = b
        let mut ii = 0usize;
        for i in 1..=n {
            let ip = self.indx[i];
            let mut sum = x.value(ip);
            x.set_value(ip, x.value(i));
            if ii != 0 {
                for j in ii..i { sum -= self.lu.value(i, j) * x.value(j); }
            } else if sum != 0.0 { ii = i; }
            x.set_value(i, sum);
        }

        // Back substitution: solve U*x = y
        for i in (1..=n).rev() {
            let mut sum = x.value(i);
            for j in (i + 1)..=n { sum -= self.lu.value(i, j) * x.value(j); }
            x.set_value(i, sum / self.lu.value(i, i));
        }
        Ok(x)
    }

    /// Determinant of original matrix.
    pub fn determinant(&self) -> f64 {
        let n = self.lu.row_count();
        let mut d = self.det_sign;
        for i in 1..=n { d *= self.lu.value(i, i); }
        d
    }

    /// Invert original matrix.
    pub fn invert(&self) -> Result<MathMatrix, MathStatus> {
        let n = self.lu.row_count();
        let mut inv = MathMatrix::new(1, n, 1, n);
        let mut col = MathVector::new(1, n);
        for j in 1..=n {
            for i in 1..=n { col.set_value(i, 0.0); }
            col.set_value(j, 1.0);
            let x = self.solve(&col)?;
            for i in 1..=n { inv.set_value(i, j, x.value(i)); }
        }
        Ok(inv)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crout_solve() {
        let mut a = MathMatrix::new(1, 3, 1, 3);
        a.set_value(1,1,2.0); a.set_value(1,2,1.0); a.set_value(1,3,0.0);
        a.set_value(2,1,1.0); a.set_value(2,2,3.0); a.set_value(2,3,1.0);
        a.set_value(3,1,0.0); a.set_value(3,2,1.0); a.set_value(3,3,2.0);
        let crout = Crout::new(&a).unwrap();
        let b = MathVector::from_slice(&[5.0, 11.0, 5.0]); // x = [1,3,1]
        let x = crout.solve(&b).unwrap();
        assert!((x.value(1) - 1.0).abs() < 1e-14);
        assert!((x.value(2) - 3.0).abs() < 1e-14);
        assert!((x.value(3) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn crout_det() {
        let mut a = MathMatrix::new(1, 2, 1, 2);
        a.set_value(1,1,4.0); a.set_value(1,2,7.0);
        a.set_value(2,1,2.0); a.set_value(2,2,6.0);
        let crout = Crout::new(&a).unwrap();
        assert!((crout.determinant() - 10.0).abs() < 1e-14); // 4*6 - 7*2 = 10
    }
}
