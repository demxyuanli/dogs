//! 3×3 matrix. Source: `gp_Mat.hxx`

use std::ops::{Add, Div, Mul, Sub};
use crate::precision::RESOLUTION;
use crate::gp::xyz::GpXyz;

/// 3×3 matrix stored row-major.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpMat {
    pub m: [[f64; 3]; 3],
}

impl Default for GpMat {
    fn default() -> Self { Self::identity() }
}

impl GpMat {
    pub const fn zero() -> Self {
        Self { m: [[0.0; 3]; 3] }
    }

    pub const fn identity() -> Self {
        Self {
            m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    pub const fn new(a11: f64, a12: f64, a13: f64, a21: f64, a22: f64, a23: f64, a31: f64, a32: f64, a33: f64) -> Self {
        Self { m: [[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]] }
    }

    pub fn from_cols(c1: &GpXyz, c2: &GpXyz, c3: &GpXyz) -> Self {
        Self { m: [[c1.x, c2.x, c3.x], [c1.y, c2.y, c3.y], [c1.z, c2.z, c3.z]] }
    }

    #[inline] pub fn value(&self, row: usize, col: usize) -> f64 { self.m[row - 1][col - 1] }

    #[inline] pub fn row(&self, r: usize) -> GpXyz {
        let r = r - 1;
        GpXyz::new(self.m[r][0], self.m[r][1], self.m[r][2])
    }

    #[inline] pub fn column(&self, c: usize) -> GpXyz {
        let c = c - 1;
        GpXyz::new(self.m[0][c], self.m[1][c], self.m[2][c])
    }

    #[inline] pub fn diagonal(&self) -> GpXyz {
        GpXyz::new(self.m[0][0], self.m[1][1], self.m[2][2])
    }

    pub fn set_identity(&mut self) { *self = Self::identity(); }

    pub fn set_diagonal(&mut self, d: &GpXyz) {
        *self = Self::zero();
        self.m[0][0] = d.x;
        self.m[1][1] = d.y;
        self.m[2][2] = d.z;
    }

    pub fn set_scale(&mut self, s: f64) {
        *self = Self::zero();
        self.m[0][0] = s;
        self.m[1][1] = s;
        self.m[2][2] = s;
    }

    /// Rodrigues' rotation formula: R = I + sin(a)*K + (1-cos(a))*K²
    pub fn set_rotation(&mut self, axis: &GpXyz, angle: f64) {
        let x = axis.x;
        let y = axis.y;
        let z = axis.z;
        let s = angle.sin();
        let c = angle.cos();
        let t = 1.0 - c;

        self.m[0][0] = t * x * x + c;
        self.m[0][1] = t * x * y - s * z;
        self.m[0][2] = t * x * z + s * y;
        self.m[1][0] = t * x * y + s * z;
        self.m[1][1] = t * y * y + c;
        self.m[1][2] = t * y * z - s * x;
        self.m[2][0] = t * x * z - s * y;
        self.m[2][1] = t * y * z + s * x;
        self.m[2][2] = t * z * z + c;
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut r = *self;
        for i in 0..3 { for j in 0..3 { r.m[i][j] += other.m[i][j]; } }
        r
    }

    pub fn subtract(&self, other: &Self) -> Self {
        let mut r = *self;
        for i in 0..3 { for j in 0..3 { r.m[i][j] -= other.m[i][j]; } }
        r
    }

    pub fn multiply_scalar(&self, s: f64) -> Self {
        let mut r = *self;
        for i in 0..3 { for j in 0..3 { r.m[i][j] *= s; } }
        r
    }

    /// this * other
    pub fn multiply(&self, other: &Self) -> Self {
        let mut r = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                r.m[i][j] = self.m[i][0] * other.m[0][j]
                          + self.m[i][1] * other.m[1][j]
                          + self.m[i][2] * other.m[2][j];
            }
        }
        r
    }

    /// other * this
    pub fn pre_multiply(&self, other: &Self) -> Self { other.multiply(self) }

    /// Multiply matrix by column vector.
    pub fn multiplied(&self, v: &GpXyz) -> GpXyz {
        GpXyz::new(
            self.m[0][0] * v.x + self.m[0][1] * v.y + self.m[0][2] * v.z,
            self.m[1][0] * v.x + self.m[1][1] * v.y + self.m[1][2] * v.z,
            self.m[2][0] * v.x + self.m[2][1] * v.y + self.m[2][2] * v.z,
        )
    }

    pub fn divide(&self, s: f64) -> Self { self.multiply_scalar(1.0 / s) }

    pub fn determinant(&self) -> f64 {
        self.m[0][0] * (self.m[1][1] * self.m[2][2] - self.m[1][2] * self.m[2][1])
            - self.m[0][1] * (self.m[1][0] * self.m[2][2] - self.m[1][2] * self.m[2][0])
            + self.m[0][2] * (self.m[1][0] * self.m[2][1] - self.m[1][1] * self.m[2][0])
    }

    pub fn transpose(&self) -> Self {
        Self { m: [
            [self.m[0][0], self.m[1][0], self.m[2][0]],
            [self.m[0][1], self.m[1][1], self.m[2][1]],
            [self.m[0][2], self.m[1][2], self.m[2][2]],
        ]}
    }

    pub fn invert(&self) -> Self {
        let a = self.m;
        let det = self.determinant();
        let inv_det = 1.0 / det;
        Self { m: [
            [
                (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_det,
                -(a[0][1] * a[2][2] - a[0][2] * a[2][1]) * inv_det,
                (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_det,
            ],
            [
                -(a[1][0] * a[2][2] - a[1][2] * a[2][0]) * inv_det,
                (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_det,
                -(a[0][0] * a[1][2] - a[0][2] * a[1][0]) * inv_det,
            ],
            [
                (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_det,
                -(a[0][0] * a[2][1] - a[0][1] * a[2][0]) * inv_det,
                (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_det,
            ],
        ]}
    }

    pub fn is_singular(&self) -> bool { self.determinant().abs() <= RESOLUTION }
}

// --- std::ops ---

impl Add<&GpMat> for &GpMat {
    type Output = GpMat;
    fn add(self, rhs: &GpMat) -> GpMat { self.add(rhs) }
}

impl Sub<&GpMat> for &GpMat {
    type Output = GpMat;
    fn sub(self, rhs: &GpMat) -> GpMat { self.subtract(rhs) }
}

impl Mul<&GpMat> for &GpMat {
    type Output = GpMat;
    fn mul(self, rhs: &GpMat) -> GpMat { self.multiply(rhs) }
}

impl Mul<f64> for &GpMat {
    type Output = GpMat;
    fn mul(self, rhs: f64) -> GpMat { self.multiply_scalar(rhs) }
}

impl Mul<&GpMat> for f64 {
    type Output = GpMat;
    fn mul(self, rhs: &GpMat) -> GpMat { rhs.multiply_scalar(self) }
}

impl Div<f64> for &GpMat {
    type Output = GpMat;
    fn div(self, rhs: f64) -> GpMat { self.divide(rhs) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn identity_det_1() {
        assert!((GpMat::identity().determinant() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rotation_90_z() {
        let mut m = GpMat::zero();
        m.set_rotation(&GpXyz::new(0.0, 0.0, 1.0), FRAC_PI_2);
        // Rotate (1,0,0) by 90° around Z -> (0,1,0)
        let v = GpXyz::new(1.0, 0.0, 0.0);
        let r = m.multiplied(&v);
        assert!((r.x - 0.0).abs() < 1e-12);
        assert!((r.y - 1.0).abs() < 1e-12);
        assert!((r.z - 0.0).abs() < 1e-12);
    }

    #[test]
    fn multiply_associative() {
        let a = GpMat::new(1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0);
        let b = GpMat::new(4.0, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0, 3.0, 2.0);
        let c = GpMat::new(1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0);
        let ab = &a * &b;
        let ab_c = &ab * &c;
        let a_bc = &a * &(&b * &c);
        assert_eq!(ab_c, a_bc);
    }

    #[test]
    fn invert_roundtrip() {
        let m = GpMat::new(1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0);
        assert!(!m.is_singular());
        let inv = m.invert();
        let prod = &m * &inv;
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((prod.m[i][j] - expected).abs() < 1e-12,
                    "prod[{i}][{j}] = {}, expected {expected}", prod.m[i][j]);
            }
        }
    }
}
