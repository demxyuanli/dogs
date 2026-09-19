//! 3D Cartesian coordinate. Source: `gp_XYZ.hxx`

use std::ops::{Add, BitXor, Div, Mul, Neg, Sub};
use crate::precision::RESOLUTION;
use crate::gp::mat::GpMat;

/// Three-dimensional Cartesian coordinate (x, y, z).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpXyz {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Default for GpXyz {
    fn default() -> Self { Self::zero() }
}

impl GpXyz {
    pub const fn zero() -> Self { Self { x: 0.0, y: 0.0, z: 0.0 } }

    pub const fn new(x: f64, y: f64, z: f64) -> Self { Self { x, y, z } }

    #[inline] pub fn set_coord(&mut self, x: f64, y: f64, z: f64) { self.x = x; self.y = y; self.z = z; }

    #[inline] pub fn set_x(&mut self, v: f64) { self.x = v; }
    #[inline] pub fn set_y(&mut self, v: f64) { self.y = v; }
    #[inline] pub fn set_z(&mut self, v: f64) { self.z = v; }

    #[inline] pub fn x(&self) -> f64 { self.x }
    #[inline] pub fn y(&self) -> f64 { self.y }
    #[inline] pub fn z(&self) -> f64 { self.z }

    /// `gp_XYZ::Coord(theIndex)` with a 0-based index (0..3).
    #[inline] pub fn coord(&self, the_index: usize) -> f64 {
        match the_index {
            0 => self.x,
            1 => self.y,
            _ => self.z,
        }
    }

    pub fn modulus(&self) -> f64 { (self.x * self.x + self.y * self.y + self.z * self.z).sqrt() }
    pub fn square_modulus(&self) -> f64 { self.x * self.x + self.y * self.y + self.z * self.z }

    pub fn is_equal(&self, other: &Self) -> bool {
        (self.x - other.x).abs() <= RESOLUTION
            && (self.y - other.y).abs() <= RESOLUTION
            && (self.z - other.z).abs() <= RESOLUTION
    }

    pub fn add(&self, other: &Self) -> Self { Self::new(self.x + other.x, self.y + other.y, self.z + other.z) }
    pub fn added(&self, other: &Self) -> Self { self.add(other) }
    pub fn subtract(&self, other: &Self) -> Self { Self::new(self.x - other.x, self.y - other.y, self.z - other.z) }
    pub fn subtracted(&self, other: &Self) -> Self { self.subtract(other) }

    pub fn multiply_scalar(&self, s: f64) -> Self { Self::new(self.x * s, self.y * s, self.z * s) }
    pub fn multiplied(&self, s: f64) -> Self { self.multiply_scalar(s) }
    pub fn multiply(&self, s: f64) -> Self { self.multiply_scalar(s) }
    pub fn divide(&self, s: f64) -> Self { Self::new(self.x / s, self.y / s, self.z / s) }
    pub fn divided(&self, s: f64) -> Self { self.divide(s) }

    /// Component-wise multiplication: (x1*x2, y1*y2, z1*z2)
    pub fn multiply_xyz(&self, other: &Self) -> Self {
        Self::new(self.x * other.x, self.y * other.y, self.z * other.z)
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// In-place cross product: self = left x right. Aliasing-safe via caching.
    pub fn set_cross(&mut self, left: &Self, right: &Self) {
        let x = left.y * right.z - left.z * right.y;
        let y = left.z * right.x - left.x * right.z;
        let z = left.x * right.y - left.y * right.x;
        self.x = x;
        self.y = y;
        self.z = z;
    }

    /// Return self x other as a new GpXyz.
    pub fn crossed(&self, other: &Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn cross_magnitude(&self, other: &Self) -> f64 { self.crossed(other).modulus() }
    pub fn cross_square_magnitude(&self, other: &Self) -> f64 { self.crossed(other).square_modulus() }

    /// In-place cross-cross: self = (a x b) x c
    pub fn set_cross_cross(&mut self, a: &Self, b: &Self, c: &Self) {
        let ab = a.crossed(b);
        self.set_cross(&ab, c);
    }

    /// Returns (a x b) x c
    pub fn cross_crossed(&self, b: &Self, c: &Self) -> Self {
        let ab = self.crossed(b);
        ab.crossed(c)
    }

    /// Triple product: self · (a x b)
    pub fn dot_cross(&self, a: &Self, b: &Self) -> f64 { self.dot(&a.crossed(b)) }

    /// Normalize in-place. Returns true on success, false if zero-length.
    pub fn normalize(&mut self) -> bool {
        let m = self.modulus();
        if m <= RESOLUTION { return false; }
        self.x /= m;
        self.y /= m;
        self.z /= m;
        true
    }

    pub fn normalized(&self) -> Self {
        let mut r = *self;
        r.normalize();
        r
    }

    pub fn reverse(&mut self) { self.x = -self.x; self.y = -self.y; self.z = -self.z; }
    pub fn reversed(&self) -> Self { Self::new(-self.x, -self.y, -self.z) }

    /// In-place matrix * column vector: self = M * self
    pub fn multiply_mat(&mut self, m: &GpMat) {
        let x = m.m[0][0] * self.x + m.m[0][1] * self.y + m.m[0][2] * self.z;
        let y = m.m[1][0] * self.x + m.m[1][1] * self.y + m.m[1][2] * self.z;
        let z = m.m[2][0] * self.x + m.m[2][1] * self.y + m.m[2][2] * self.z;
        self.x = x;
        self.y = y;
        self.z = z;
    }

    /// Return M * self as a new GpXyz.
    pub fn multiplied_mat(&self, m: &GpMat) -> Self {
        Self::new(
            m.m[0][0] * self.x + m.m[0][1] * self.y + m.m[0][2] * self.z,
            m.m[1][0] * self.x + m.m[1][1] * self.y + m.m[1][2] * self.z,
            m.m[2][0] * self.x + m.m[2][1] * self.y + m.m[2][2] * self.z,
        )
    }

    /// self = a1*v1 + a2*v2
    pub fn set_linear_form(&mut self, a1: f64, v1: &Self, a2: f64, v2: &Self) {
        self.x = a1 * v1.x + a2 * v2.x;
        self.y = a1 * v1.y + a2 * v2.y;
        self.z = a1 * v1.z + a2 * v2.z;
    }

    /// self = a1*v1 + a2*v2 + a3*v3
    pub fn set_linear_form_3(&mut self, a1: f64, v1: &Self, a2: f64, v2: &Self, a3: f64, v3: &Self) {
        self.x = a1 * v1.x + a2 * v2.x + a3 * v3.x;
        self.y = a1 * v1.y + a2 * v2.y + a3 * v3.y;
        self.z = a1 * v1.z + a2 * v2.z + a3 * v3.z;
    }

    /// self = a1*v1 + a2*v2 + a3*v3 + a4*v4
    pub fn set_linear_form_4(&mut self, a1: f64, v1: &Self, a2: f64, v2: &Self, a3: f64, v3: &Self, a4: f64, v4: &Self) {
        self.x = a1 * v1.x + a2 * v2.x + a3 * v3.x + a4 * v4.x;
        self.y = a1 * v1.y + a2 * v2.y + a3 * v3.y + a4 * v4.y;
        self.z = a1 * v1.z + a2 * v2.z + a3 * v3.z + a4 * v4.z;
    }

    #[cfg(feature = "serde_json")]
    pub fn dump_json(&self) -> serde_json::Value {
        serde_json::json!({"x": self.x, "y": self.y, "z": self.z})
    }
}

// --- std::ops ---

impl Add<&GpXyz> for &GpXyz {
    type Output = GpXyz;
    fn add(self, rhs: &GpXyz) -> GpXyz { self.add(rhs) }
}

impl Sub<&GpXyz> for &GpXyz {
    type Output = GpXyz;
    fn sub(self, rhs: &GpXyz) -> GpXyz { self.subtract(rhs) }
}

impl Mul<f64> for &GpXyz {
    type Output = GpXyz;
    fn mul(self, rhs: f64) -> GpXyz { self.multiply_scalar(rhs) }
}

impl Mul<&GpXyz> for f64 {
    type Output = GpXyz;
    fn mul(self, rhs: &GpXyz) -> GpXyz { rhs.multiply_scalar(self) }
}

impl Div<f64> for &GpXyz {
    type Output = GpXyz;
    fn div(self, rhs: f64) -> GpXyz { self.divide(rhs) }
}

impl Neg for &GpXyz {
    type Output = GpXyz;
    fn neg(self) -> GpXyz { self.reversed() }
}

impl BitXor<&GpXyz> for &GpXyz {
    type Output = GpXyz;
    fn bitxor(self, rhs: &GpXyz) -> GpXyz { self.crossed(rhs) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_unit_vectors() {
        let i = GpXyz::new(1.0, 0.0, 0.0);
        let j = GpXyz::new(0.0, 1.0, 0.0);
        let k = &i ^ &j;
        assert!((k.x - 0.0).abs() < 1e-12);
        assert!((k.y - 0.0).abs() < 1e-12);
        assert!((k.z - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cross_anticommutative() {
        let a = GpXyz::new(1.0, 2.0, 3.0);
        let b = GpXyz::new(4.0, 5.0, 6.0);
        let ab = &a ^ &b;
        let ba = &b ^ &a;
        assert_eq!(ab, -&ba);
    }

    #[test]
    fn normalize_unit() {
        let mut v = GpXyz::new(3.0, 0.0, 0.0);
        assert!(v.normalize());
        assert!((v.modulus() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn normalize_zero_fails() {
        let mut v = GpXyz::zero();
        assert!(!v.normalize());
    }

    #[test]
    fn dot_cross_is_triple() {
        let a = GpXyz::new(1.0, 0.0, 0.0);
        let b = GpXyz::new(0.0, 2.0, 0.0);
        let c = GpXyz::new(0.0, 0.0, 3.0);
        assert!((a.dot_cross(&b, &c) - 6.0).abs() < 1e-12);
    }
}
