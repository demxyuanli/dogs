//! Dynamic f64 matrix with 1-based indexing and Gauss elimination.
//! Source: `math_Matrix.hxx` + `math_Matrix.cxx`
//!
//! Storage: column-major flat Vec<f64> matching OCCT's NCollection_Array2 layout.

use crate::MathVector;

/// Dynamic matrix with configurable row/column bounds.
/// Elements stored in column-major order.
#[derive(Debug, Clone)]
pub struct MathMatrix {
    data: Vec<f64>,
    lower_row: usize,
    upper_row: usize,
    lower_col: usize,
    upper_col: usize,
}

impl MathMatrix {
    // ---- constructors ----

    /// Uninitialized matrix with given bounds. Source: `math_Matrix.hxx:86`
    pub fn new(lr: usize, ur: usize, lc: usize, uc: usize) -> Self {
        assert!(ur >= lr && uc >= lc, "MathMatrix: invalid bounds");
        let rows = ur - lr + 1;
        let cols = uc - lc + 1;
        Self { data: vec![0.0; rows * cols], lower_row: lr, upper_row: ur, lower_col: lc, upper_col: uc }
    }

    /// Matrix with initial value. Source: `math_Matrix.hxx:94`
    pub fn with_init(lr: usize, ur: usize, lc: usize, uc: usize, init: f64) -> Self {
        let rows = ur - lr + 1;
        let cols = uc - lc + 1;
        Self { data: vec![init; rows * cols], lower_row: lr, upper_row: ur, lower_col: lc, upper_col: uc }
    }

    // ---- sizing ----

    #[inline] pub fn row_count(&self) -> usize { self.upper_row - self.lower_row + 1 }
    #[inline] pub fn col_count(&self) -> usize { self.upper_col - self.lower_col + 1 }
    #[inline] pub fn lower_row(&self) -> usize { self.lower_row }
    #[inline] pub fn upper_row(&self) -> usize { self.upper_row }
    #[inline] pub fn lower_col(&self) -> usize { self.lower_col }
    #[inline] pub fn upper_col(&self) -> usize { self.upper_col }
    #[inline] pub fn is_square(&self) -> bool { self.row_count() == self.col_count() }

    /// Initialize all elements. Source: `math_Matrix.hxx:117`
    pub fn init(&mut self, val: f64) { self.data.fill(val); }

    // ---- element access ----

    /// Get element at (row, col), 1-based. Source: `math_Matrix.hxx` (Value)
    #[inline]
    pub fn value(&self, row: usize, col: usize) -> f64 {
        let ridx = row - self.lower_row;
        let cidx = col - self.lower_col;
        let nr = self.row_count();
        self.data[cidx * nr + ridx]
    }

    /// Set element at (row, col), 1-based. Source: `math_Matrix.hxx` (ChangeValue)
    #[inline]
    pub fn set_value(&mut self, row: usize, col: usize, val: f64) {
        let ridx = row - self.lower_row;
        let cidx = col - self.lower_col;
        let nr = self.row_count();
        self.data[cidx * nr + ridx] = val;
    }

    /// Set entire row. Source: `math_Matrix.hxx` (SetRow)
    pub fn set_row(&mut self, row: usize, v: &MathVector) {
        assert_eq!(v.len(), self.col_count(), "SetRow: dimension mismatch");
        for c in 1..=self.col_count() { self.set_value(row, c, v.value(c)); }
    }

    /// Set entire column. Source: `math_Matrix.hxx` (SetCol)
    pub fn set_col(&mut self, col: usize, v: &MathVector) {
        assert_eq!(v.len(), self.row_count(), "SetCol: dimension mismatch");
        for r in 1..=self.row_count() { self.set_value(r, col, v.value(r)); }
    }

    /// Set diagonal. Source: `math_Matrix.hxx` (SetDiag)
    pub fn set_diag(&mut self, val: f64) {
        let n = self.row_count().min(self.col_count());
        for i in 1..=n { self.set_value(i, i, val); }
    }

    // ---- scalar ops ----

    /// Multiply all elements by scalar in-place. Source: `math_Matrix.hxx:183`
    pub fn multiply_scalar(&mut self, scalar: f64) {
        for x in &mut self.data { *x *= scalar; }
    }

    /// New matrix = self * scalar. Source: `math_Matrix.hxx:189`
    pub fn multiplied_scalar(&self, scalar: f64) -> Self {
        let mut result = self.clone();
        result.multiply_scalar(scalar);
        result
    }

    // ---- matrix arithmetic ----

    /// Check dimension compatibility for binary ops.
    fn assert_same_dims(&self, other: &Self, op: &str) {
        assert_eq!(self.row_count(), other.row_count(), "{op}: row count mismatch");
        assert_eq!(self.col_count(), other.col_count(), "{op}: column count mismatch");
    }

    /// `self += other`. Source: `math_Matrix.hxx` (Add)
    pub fn add(&mut self, other: &Self) {
        self.assert_same_dims(other, "Add");
        for (a, &b) in self.data.iter_mut().zip(&other.data) { *a += b; }
    }

    /// `self + other`.
    pub fn added(&self, other: &Self) -> Self { let mut r = self.clone(); r.add(other); r }

    /// `self -= other`. Source: `math_Matrix.hxx` (Subtract)
    pub fn subtract(&mut self, other: &Self) {
        self.assert_same_dims(other, "Subtract");
        for (a, &b) in self.data.iter_mut().zip(&other.data) { *a -= b; }
    }

    /// `self - other`.
    pub fn subtracted(&self, other: &Self) -> Self { let mut r = self.clone(); r.subtract(other); r }

    // ---- matrix multiply ----

    /// `self = left * right`. Source: `math_Matrix.hxx` (Multiply)
    pub fn multiply_mat(&mut self, left: &Self, right: &Self) {
        let m = left.row_count();
        let n = left.col_count();
        let p = right.col_count();
        assert_eq!(n, right.row_count(), "Multiply: inner dimensions mismatch");
        assert_eq!(m, self.row_count(), "Multiply: result row count mismatch");
        assert_eq!(p, self.col_count(), "Multiply: result col count mismatch");

        for i in 1..=m {
            for j in 1..=p {
                let mut sum = 0.0;
                for k in 1..=n { sum += left.value(i, k) * right.value(k, j); }
                self.set_value(i, j, sum);
            }
        }
    }

    /// `self * other`.
    pub fn multiplied_mat(&self, other: &Self) -> Self {
        let mut result = Self::new(self.lower_row, self.upper_row, other.lower_col, other.upper_col);
        result.multiply_mat(self, other);
        result
    }

    /// Multiply in-place: `self = self * other`. Source: `math_Matrix.hxx` (Multiply in-place)
    pub fn multiply_assign(&mut self, other: &Self) {
        let tmp = self.multiplied_mat(other);
        *self = tmp;
    }

    // ---- transpose ----

    /// Transpose in-place. Source: `math_Matrix.hxx:159`
    pub fn transpose(&mut self) {
        assert!(self.is_square(), "Transpose: not square");
        for i in 1..self.row_count() {
            for j in (i + 1)..=self.col_count() {
                let a = self.value(i, j);
                let b = self.value(j, i);
                self.set_value(i, j, b);
                self.set_value(j, i, a);
            }
        }
    }

    /// Transposed copy.
    pub fn transposed(&self) -> Self {
        let mut result = self.clone();
        result.transpose();
        result
    }

    // ---- determinant ----

    /// Determinant via LU decomposition (Doolittle).
    /// Source: `math_Matrix.cxx` (Determinant)
    pub fn determinant(&self) -> Result<f64, &'static str> {
        if !self.is_square() { return Err("Determinant: not square"); }
        let n = self.row_count();
        if n == 1 { return Ok(self.value(1, 1)); }
        if n == 2 {
            return Ok(self.value(1,1) * self.value(2,2) - self.value(1,2) * self.value(2,1));
        }
        if n == 3 {
            let a = self.value(1,1); let b = self.value(1,2); let c = self.value(1,3);
            let d = self.value(2,1); let e = self.value(2,2); let f = self.value(2,3);
            let g = self.value(3,1); let h = self.value(3,2); let i = self.value(3,3);
            return Ok(a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g));
        }
        // LU decomposition for n >= 4
        let mut lu = self.clone();
        let mut det = 1.0;
        let mut sign = 1.0;
        for k in 1..=n {
            // Partial pivoting
            let mut max_val = lu.value(k, k).abs();
            let mut max_row = k;
            for i in (k + 1)..=n {
                let v = lu.value(i, k).abs();
                if v > max_val { max_val = v; max_row = i; }
            }
            if max_row != k {
                for j in 1..=n {
                    let tmp = lu.value(k, j);
                    lu.set_value(k, j, lu.value(max_row, j));
                    lu.set_value(max_row, j, tmp);
                }
                sign = -sign;
            }
            let pivot = lu.value(k, k);
            if pivot.abs() < 1e-30 { return Ok(0.0); }
            det *= pivot;
            for i in (k + 1)..=n {
                let factor = lu.value(i, k) / pivot;
                for j in k..=n {
                    let v = lu.value(i, j) - factor * lu.value(k, j);
                    lu.set_value(i, j, v);
                }
            }
        }
        Ok(sign * det)
    }

    // ---- invert (Gauss-Jordan) ----

    /// Invert in-place via Gauss-Jordan elimination.
    /// Source: `math_Matrix.cxx` (Invert)
    pub fn invert(&mut self) -> Result<(), &'static str> {
        if !self.is_square() { return Err("Invert: not square"); }
        let n = self.row_count();
        let mut inv = Self::with_init(1, n, 1, n, 0.0);
        for i in 1..=n { inv.set_value(i, i, 1.0); }

        // Augmented matrix [self | inv], apply Gauss-Jordan
        for k in 1..=n {
            // Pivot
            let mut max_val = self.value(k, k).abs();
            let mut max_row = k;
            for i in (k + 1)..=n {
                let v = self.value(i, k).abs();
                if v > max_val { max_val = v; max_row = i; }
            }
            if max_row != k {
                for j in 1..=n {
                    let tmp = self.value(k, j);
                    self.set_value(k, j, self.value(max_row, j));
                    self.set_value(max_row, j, tmp);
                    let ti = inv.value(k, j);
                    inv.set_value(k, j, inv.value(max_row, j));
                    inv.set_value(max_row, j, ti);
                }
            }
            let pivot = self.value(k, k);
            if pivot.abs() < 1e-30 { return Err("Invert: singular matrix"); }
            for j in 1..=n {
                self.set_value(k, j, self.value(k, j) / pivot);
                inv.set_value(k, j, inv.value(k, j) / pivot);
            }
            for i in 1..=n {
                if i == k { continue; }
                let factor = self.value(i, k);
                for j in 1..=n {
                    self.set_value(i, j, self.value(i, j) - factor * self.value(k, j));
                    inv.set_value(i, j, inv.value(i, j) - factor * inv.value(k, j));
                }
            }
        }
        *self = inv;
        Ok(())
    }

    /// Inverted copy.
    pub fn inverted(&self) -> Result<Self, &'static str> { let mut r = self.clone(); r.invert()?; Ok(r) }

    // ---- Gauss elimination (solve Ax = b) ----

    /// Solve linear system `Ax = b` via Gauss elimination with partial pivoting.
    /// Source: `math_Gauss.cxx`
    pub fn solve(&self, b: &MathVector) -> Result<MathVector, &'static str> {
        if !self.is_square() { return Err("Solve: not square"); }
        let n = self.row_count();
        assert_eq!(b.len(), n, "Solve: b dimension mismatch");

        let mut a = self.clone();
        let mut x = b.clone();

        // Forward elimination with partial pivoting
        for k in 1..=n {
            let mut max_val = a.value(k, k).abs();
            let mut max_row = k;
            for i in (k + 1)..=n {
                let v = a.value(i, k).abs();
                if v > max_val { max_val = v; max_row = i; }
            }
            if max_row != k {
                for j in k..=n {
                    let tmp = a.value(k, j);
                    a.set_value(k, j, a.value(max_row, j));
                    a.set_value(max_row, j, tmp);
                }
                let tb = x.value(k);
                x.set_value(k, x.value(max_row));
                x.set_value(max_row, tb);
            }
            let pivot = a.value(k, k);
            // `math_Gauss.hxx:45-49`: `MinPivot = 1.0e-20` ("If the largest pivot
            // found is less than MinPivot the matrix A is considered singular").
            if pivot.abs() < 1.0e-20 { return Err("Solve: singular matrix"); }
            for i in (k + 1)..=n {
                let factor = a.value(i, k) / pivot;
                a.set_value(i, k, 0.0);
                for j in (k + 1)..=n {
                    let v = a.value(i, j) - factor * a.value(k, j);
                    a.set_value(i, j, v);
                }
                x.set_value(i, x.value(i) - factor * x.value(k));
            }
        }

        // Back substitution
        for i in (1..=n).rev() {
            let mut sum = x.value(i);
            for j in (i + 1)..=n { sum -= a.value(i, j) * x.value(j); }
            x.set_value(i, sum / a.value(i, i));
        }
        Ok(x)
    }

    // ---- row/col operations ----

    /// Swap rows. Source: `math_Matrix.hxx` (ExchangeRow)
    pub fn exchange_row(&mut self, r1: usize, r2: usize) {
        for c in 1..=self.col_count() {
            let tmp = self.value(r1, c);
            self.set_value(r1, c, self.value(r2, c));
            self.set_value(r2, c, tmp);
        }
    }

    /// Swap columns. Source: `math_Matrix.hxx` (ExchangeCol)
    pub fn exchange_col(&mut self, c1: usize, c2: usize) {
        for r in 1..=self.row_count() {
            let tmp = self.value(r, c1);
            self.set_value(r, c1, self.value(r, c2));
            self.set_value(r, c2, tmp);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_and_access() {
        let m = MathMatrix::with_init(1, 3, 1, 3, 1.0);
        assert_eq!(m.row_count(), 3);
        assert_eq!(m.col_count(), 3);
        assert_eq!(m.value(2, 3), 1.0);
    }

    #[test]
    fn determinant_3x3() {
        let mut m = MathMatrix::new(1, 3, 1, 3);
        // identity → det = 1
        m.set_value(1,1,1.0); m.set_value(2,2,1.0); m.set_value(3,3,1.0);
        assert!((m.determinant().unwrap() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn transpose() {
        let mut m = MathMatrix::new(1, 2, 1, 2);
        m.set_value(1,1,1.0); m.set_value(1,2,2.0);
        m.set_value(2,1,3.0); m.set_value(2,2,4.0);
        let t = m.transposed();
        assert_eq!(t.value(1, 2), 3.0);
        assert_eq!(t.value(2, 1), 2.0);
    }

    #[test]
    fn solve_2x2() {
        let mut m = MathMatrix::new(1, 2, 1, 2);
        m.set_value(1,1,2.0); m.set_value(1,2,1.0);
        m.set_value(2,1,1.0); m.set_value(2,2,3.0);
        let b = MathVector::from_slice(&[5.0, 6.0]); // 2x+y=5, x+3y=6 → x=1.8, y=1.4
        let x = m.solve(&b).unwrap();
        assert!((x.value(1) - 1.8).abs() < 1e-14);
        assert!((x.value(2) - 1.4).abs() < 1e-14);
    }

    #[test]
    fn invert_roundtrip() {
        let mut m = MathMatrix::new(1, 2, 1, 2);
        m.set_value(1,1,4.0); m.set_value(1,2,7.0);
        m.set_value(2,1,2.0); m.set_value(2,2,6.0);
        let inv = m.inverted().unwrap();
        let prod = m.multiplied_mat(&inv);
        assert!((prod.value(1,1)-1.0).abs() < 1e-13);
        assert!(prod.value(1,2).abs() < 1e-13);
        assert!((prod.value(2,2)-1.0).abs() < 1e-13);
    }
}
