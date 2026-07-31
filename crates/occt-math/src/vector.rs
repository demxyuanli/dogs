//! Dynamic-length real vector with 1-based indexing.
//! Source: `math_VectorBase.hxx` (template instantiation for double)
//! OCCT: `math_Vector = math_VectorBase<double>`

use crate::MathMatrix;

/// Dynamic f64 vector. Uses 0-based Vec internally, exposes 1-based indexing
/// to match OCCT's `math_Vector`.
///
/// Size is fixed at construction; resizing is explicit via `resize()`.
#[derive(Debug, Clone, PartialEq)]
pub struct MathVector {
    data: Vec<f64>,
    lower: usize,
}

impl MathVector {
    // ---- constructors ----

    /// Vector with indices `lower..=upper` (inclusive). Source: `math_VectorBase.hxx:74`
    pub fn new(lower: usize, upper: usize) -> Self {
        assert!(upper >= lower, "MathVector: upper < lower");
        let len = upper - lower + 1;
        Self { data: vec![0.0; len], lower }
    }

    /// Vector with values initialized. Source: `math_VectorBase.hxx:78`
    pub fn with_init(lower: usize, upper: usize, init: f64) -> Self {
        let mut v = Self::new(lower, upper);
        v.init(init);
        v
    }

    /// From slice with 1-based indexing (lower=1). Source: implicit
    pub fn from_slice(values: &[f64]) -> Self {
        Self { data: values.to_vec(), lower: 1 }
    }

    /// From GP XYZ. Source: `math_VectorBase.hxx:88`
    pub fn from_xyz(xyz: &occt_core::GpXyz) -> Self {
        let mut v = Self::new(1, 3);
        v.set_value(1, xyz.x);
        v.set_value(2, xyz.y);
        v.set_value(3, xyz.z);
        v
    }

    /// From GP XY. Source: `math_VectorBase.hxx:85`
    pub fn from_xy(xy: &occt_core::GpXY) -> Self {
        let mut v = Self::new(1, 2);
        v.set_value(1, xy.x);
        v.set_value(2, xy.y);
        v
    }

    // ---- sizing ----

    /// Number of elements. Source: `math_VectorBase.hxx:100`
    #[inline] pub fn len(&self) -> usize { self.data.len() }

    /// Lower index. Source: `math_VectorBase.hxx:103`
    #[inline] pub fn lower(&self) -> usize { self.lower }

    /// Upper index. Source: `math_VectorBase.hxx:106`
    #[inline] pub fn upper(&self) -> usize { self.lower + self.data.len() - 1 }

    /// True if empty. Source: implicit
    #[inline] pub fn is_empty(&self) -> bool { self.data.is_empty() }

    /// Resize, keeping existing data. Source: `math_VectorBase.hxx:314`
    pub fn resize(&mut self, new_size: usize) {
        self.data.resize(new_size, 0.0);
    }

    // ---- initialization ----

    /// Fill all elements. Source: `math_VectorBase.hxx:91`
    pub fn init(&mut self, val: f64) {
        self.data.fill(val);
    }

    // ---- element access ----

    /// Get element at 1-based index. Source: `math_VectorBase.hxx:240`
    #[inline]
    pub fn value(&self, i: usize) -> f64 {
        assert!(i >= self.lower && i <= self.upper(), "MathVector::value index out of range");
        self.data[i - self.lower]
    }

    /// Set element at 1-based index. Source: `math_VectorBase.hxx:243`
    #[inline]
    pub fn set_value(&mut self, i: usize, val: f64) {
        assert!(i >= self.lower && i <= self.upper(), "MathVector::set_value index out of range");
        self.data[i - self.lower] = val;
    }

    /// Operator() access, immutable. Source: `math_VectorBase.hxx:245`
    #[inline]
    pub fn at(&self, i: usize) -> f64 { self.value(i) }

    /// Operator() access, mutable. Source: `math_VectorBase.hxx:247`
    #[inline]
    pub fn at_mut(&mut self, i: usize) -> &mut f64 {
        let lo = self.lower;
        assert!(i >= lo && i <= self.upper(), "MathVector::at_mut index out of range");
        &mut self.data[i - lo]
    }

    // ---- slice access ----

    /// Set sub-range from another vector. Source: `math_VectorBase.hxx:144`
    pub fn set_slice(&mut self, i1: usize, i2: usize, v: &Self) {
        assert!(i1 <= i2, "set_slice: i1 > i2");
        assert!(i2 - i1 + 1 == v.len(), "set_slice: length mismatch");
        self.copy_from(v, i1);
    }

    /// Copy values from another vector starting at offset. Source: OCCT internal
    pub fn copy_from(&mut self, src: &Self, dst_offset: usize) {
        let lo = self.lower;
        for (i, &val) in src.data.iter().enumerate() {
            self.data[dst_offset - lo + i] = val;
        }
    }

    /// Slice with inversion. Source: `math_VectorBase.hxx:151`
    pub fn slice(&self, i1: usize, i2: usize) -> Self {
        let mut result = Self::new(self.lower, self.upper());
        result.data.copy_from_slice(&self.data);
        let lo = self.lower;
        let mut a = i1 - lo;
        let mut b = i2 - lo;
        while a < b {
            result.data.swap(a, b);
            a += 1;
            b -= 1;
        }
        result
    }

    /// Raw data slice for interop. Source: implicit
    #[inline] pub fn as_slice(&self) -> &[f64] { &self.data }

    /// Mutable raw data slice. Source: implicit
    #[inline] pub fn as_mut_slice(&mut self) -> &mut [f64] { &mut self.data }

    // ---- norms ----

    /// Euclidean norm. Source: `math_VectorBase.hxx:109`
    #[inline] pub fn norm(&self) -> f64 { self.data.iter().map(|x| x * x).sum::<f64>().sqrt() }

    /// Squared Euclidean norm. Source: `math_VectorBase.hxx:112`
    #[inline] pub fn norm2(&self) -> f64 { self.data.iter().map(|x| x * x).sum() }

    /// Index of max element. Source: `math_VectorBase.hxx:116`
    pub fn max_index(&self) -> usize {
        let (i, _) = self.data.iter().enumerate().fold((0, f64::NEG_INFINITY), |(mi, mv), (i, &v)| if v > mv { (i, v) } else { (mi, mv) });
        self.lower + i
    }

    /// Index of min element. Source: `math_VectorBase.hxx:119`
    pub fn min_index(&self) -> usize {
        let (i, _) = self.data.iter().enumerate().fold((0, f64::INFINITY), |(mi, mv), (i, &v)| if v < mv { (i, v) } else { (mi, mv) });
        self.lower + i
    }

    // ---- normalize ----

    /// Normalize in-place. Source: `math_VectorBase.hxx:125`
    pub fn normalize(&mut self) {
        let n = self.norm();
        assert!(n > f64::EPSILON, "MathVector::normalize: zero vector");
        for x in &mut self.data { *x /= n; }
    }

    /// Normalized copy. Source: `math_VectorBase.hxx:132`
    pub fn normalized(&self) -> Self {
        let mut result = self.clone();
        result.normalize();
        result
    }

    // ---- invert ----

    /// Invert (1/v for each element). Source: `math_VectorBase.hxx:135`
    pub fn invert(&mut self) {
        for x in &mut self.data { *x = 1.0 / *x; }
    }

    /// Inverted copy. Source: `math_VectorBase.hxx:138`
    pub fn inverse(&self) -> Self {
        let mut result = self.clone();
        result.invert();
        result
    }

    // ---- scalar ops ----

    /// `self *= scalar`. Source: `math_VectorBase.hxx:154`
    pub fn multiply_scalar(&mut self, scalar: f64) {
        for x in &mut self.data { *x *= scalar; }
    }

    /// `self * scalar`. Source: `math_VectorBase.hxx:159`
    pub fn multiplied_scalar(&self, scalar: f64) -> Self {
        let mut result = self.clone();
        result.multiply_scalar(scalar);
        result
    }

    /// `self /= scalar`. Source: `math_VectorBase.hxx:177`
    pub fn divide_scalar(&mut self, scalar: f64) {
        for x in &mut self.data { *x /= scalar; }
    }

    /// `self / scalar`. Source: `math_VectorBase.hxx:183`
    pub fn divided_scalar(&self, scalar: f64) -> Self {
        let mut result = self.clone();
        result.divide_scalar(scalar);
        result
    }

    // ---- vector ops ----

    /// `self += other`. Source: `math_VectorBase.hxx:195`
    pub fn add(&mut self, other: &Self) {
        assert_eq!(self.len(), other.len(), "MathVector::add: dimension mismatch");
        for (a, &b) in self.data.iter_mut().zip(&other.data) { *a += b; }
    }

    /// `self + other`. Source: `math_VectorBase.hxx:202`
    pub fn added(&self, other: &Self) -> Self {
        let mut result = self.clone();
        result.add(other);
        result
    }

    /// `self -= other`. Source: `math_VectorBase.hxx:282`
    pub fn subtract(&mut self, other: &Self) {
        assert_eq!(self.len(), other.len(), "MathVector::subtract: dimension mismatch");
        for (a, &b) in self.data.iter_mut().zip(&other.data) { *a -= b; }
    }

    /// `self - other`. Source: `math_VectorBase.hxx:288`
    pub fn subtracted(&self, other: &Self) -> Self {
        let mut result = self.clone();
        result.subtract(other);
        result
    }

    /// Negate: `-self`. Source: `math_VectorBase.hxx:277`
    pub fn opposite(&self) -> Self {
        let mut result = self.clone();
        for x in &mut result.data { *x = -*x; }
        result
    }

    /// Dot product. Source: `math_VectorBase.hxx:260`
    pub fn dot(&self, other: &Self) -> f64 {
        assert_eq!(self.len(), other.len(), "MathVector::dot: dimension mismatch");
        self.data.iter().zip(&other.data).map(|(a, b)| a * b).sum()
    }

    // ---- matrix-vector ops ----

    /// `self = mat * self`. Source: `math_VectorBase.hxx:215`
    pub fn multiply_mat_vec(&mut self, mat: &MathMatrix, vec: &Self) {
        assert_eq!(mat.col_count(), vec.len(), "multiply_mat_vec: dimension mismatch");
        self.resize(mat.row_count());
        for r in 1..=mat.row_count() {
            let mut sum = 0.0;
            for c in 1..=mat.col_count() { sum += mat.value(r, c) * vec.value(c); }
            self.set_value(r, sum);
        }
    }

    /// `self = vec * mat`. Source: `math_VectorBase.hxx:211`
    pub fn multiply_vec_mat(&mut self, vec: &Self, mat: &MathMatrix) {
        assert_eq!(vec.len(), mat.row_count(), "multiply_vec_mat: dimension mismatch");
        self.resize(mat.col_count());
        for c in 1..=mat.col_count() {
            let mut sum = 0.0;
            for r in 1..=mat.row_count() { sum += vec.value(r) * mat.value(r, c); }
            self.set_value(c, sum);
        }
    }

    /// `vec * mat`. Source: `math_VectorBase.hxx:268`
    pub fn multiplied_mat(&self, mat: &MathMatrix) -> Self {
        let mut result = Self::new(1, mat.col_count());
        result.multiply_vec_mat(self, mat);
        result
    }

    // ---- initialization from another vector ----

    /// Copy from another vector (dimensions must match). Source: `math_VectorBase.hxx:251`
    pub fn initialized(&mut self, other: &Self) {
        self.data.copy_from_slice(&other.data);
        self.lower = other.lower;
    }
}

// ---- std::ops ----

impl std::ops::Add<&MathVector> for &MathVector {
    type Output = MathVector;
    fn add(self, rhs: &MathVector) -> MathVector { self.added(rhs) }
}

impl std::ops::Sub<&MathVector> for &MathVector {
    type Output = MathVector;
    fn sub(self, rhs: &MathVector) -> MathVector { self.subtracted(rhs) }
}

impl std::ops::Mul<&MathVector> for &MathVector {
    type Output = f64;
    fn mul(self, rhs: &MathVector) -> f64 { self.dot(rhs) }
}

impl std::ops::Mul<f64> for &MathVector {
    type Output = MathVector;
    fn mul(self, s: f64) -> MathVector { self.multiplied_scalar(s) }
}

impl std::ops::Div<f64> for &MathVector {
    type Output = MathVector;
    fn div(self, s: f64) -> MathVector { self.divided_scalar(s) }
}

impl std::ops::Neg for &MathVector {
    type Output = MathVector;
    fn neg(self) -> MathVector { self.opposite() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_and_access() {
        let v = MathVector::with_init(1, 5, 2.0);
        assert_eq!(v.len(), 5);
        assert_eq!(v.value(1), 2.0);
        assert_eq!(v.value(5), 2.0);
    }

    #[test]
    fn norm() {
        let v = MathVector::from_slice(&[3.0, 4.0]);
        assert!((v.norm() - 5.0).abs() < 1e-15);
        assert!((v.norm2() - 25.0).abs() < 1e-15);
    }

    #[test]
    fn arithmetic() {
        let a = MathVector::from_slice(&[1.0, 2.0, 3.0]);
        let b = MathVector::from_slice(&[4.0, 5.0, 6.0]);
        let sum = a.added(&b);
        assert_eq!(sum.value(1), 5.0);
        assert_eq!(sum.value(3), 9.0);
        let dot = a.dot(&b);
        assert!((dot - 32.0).abs() < 1e-15);
    }

    #[test]
    fn scalar_ops() {
        let v = MathVector::from_slice(&[1.0, 2.0]);
        let doubled = &v * 2.0;
        assert_eq!(doubled.value(1), 2.0);
        assert_eq!(doubled.value(2), 4.0);
        let halved = &v / 2.0;
        assert_eq!(halved.value(1), 0.5);
    }

    #[test]
    fn normalize() {
        let v = MathVector::from_slice(&[3.0, 4.0]).normalized();
        assert!((v.norm() - 1.0).abs() < 1e-15);
    }
}
