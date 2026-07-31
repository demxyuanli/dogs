use crate::gp::mat::GpMat;
use crate::gp::vec::GpVec;
use crate::precision::RESOLUTION;

/// Rotation quaternion (x, y, z, w).
#[derive(Debug, Clone, Copy)]
pub struct GpQuaternion {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Default for GpQuaternion {
    fn default() -> Self {
        Self::identity()
    }
}

impl GpQuaternion {
    /// Identity quaternion (0,0,0,1).
    pub fn identity() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 1.0,
        }
    }

    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self { x, y, z, w }
    }

    pub fn x(&self) -> f64 {
        self.x
    }

    pub fn y(&self) -> f64 {
        self.y
    }

    pub fn z(&self) -> f64 {
        self.z
    }

    pub fn w(&self) -> f64 {
        self.w
    }

    /// Set rotation that rotates `from` into `to`.
    pub fn set_rotation_vecs(&mut self, from: &GpVec, to: &GpVec) {
        let f = from.normalized();
        let t = to.normalized();
        let dot = f.dot(&t);
        if dot <= -1.0 + RESOLUTION {
            // 180-degree rotation: pick perpendicular axis
            let axis = if f.x().abs() <= RESOLUTION && f.y().abs() <= RESOLUTION {
                GpVec::new(1.0, 0.0, 0.0).crossed(&f)
            } else {
                GpVec::new(0.0, 0.0, 1.0).crossed(&f)
            };
            let half = std::f64::consts::FRAC_PI_2;
            let s = half.sin();
            let axis_n = axis.normalized();
            let axis = GpVec::new(axis_n.x(), axis_n.y(), axis_n.z());
            self.x = axis.x() * s;
            self.y = axis.y() * s;
            self.z = axis.z() * s;
            self.w = half.cos();
        } else {
            let v = f.crossed(&t);
            let s = ((1.0 + dot) * 2.0).sqrt();
            let inv = 1.0 / s;
            self.x = v.x() * inv;
            self.y = v.y() * inv;
            self.z = v.z() * inv;
            self.w = s * 0.5;
        }
    }

    /// Set from axis-angle rotation (axis normalized internally).
    pub fn set_axis_angle(&mut self, axis: &GpVec, angle: f64) {
        let a = axis.normalized();
        let half = angle * 0.5;
        let s = half.sin();
        self.x = a.x() * s;
        self.y = a.y() * s;
        self.z = a.z() * s;
        self.w = half.cos();
    }

    /// Set from a 3x3 rotation matrix stored in GpMat.
    pub fn set_matrix(&mut self, mat: &GpMat) {
        let trace = mat.value(1, 1) + mat.value(2, 2) + mat.value(3, 3);
        if trace > 0.0 {
            let s = (trace + 1.0).sqrt();
            let inv = 0.5 / s;
            self.w = s * 0.5;
            self.x = (mat.value(3, 2) - mat.value(2, 3)) * inv;
            self.y = (mat.value(1, 3) - mat.value(3, 1)) * inv;
            self.z = (mat.value(2, 1) - mat.value(1, 2)) * inv;
        } else if mat.value(1, 1) > mat.value(2, 2) && mat.value(1, 1) > mat.value(3, 3) {
            let s = (mat.value(1, 1) - mat.value(2, 2) - mat.value(3, 3) + 1.0).sqrt();
            let inv = 0.5 / s;
            self.x = s * 0.5;
            self.y = (mat.value(1, 2) + mat.value(2, 1)) * inv;
            self.z = (mat.value(1, 3) + mat.value(3, 1)) * inv;
            self.w = (mat.value(3, 2) - mat.value(2, 3)) * inv;
        } else if mat.value(2, 2) > mat.value(3, 3) {
            let s = (mat.value(2, 2) - mat.value(1, 1) - mat.value(3, 3) + 1.0).sqrt();
            let inv = 0.5 / s;
            self.x = (mat.value(1, 2) + mat.value(2, 1)) * inv;
            self.y = s * 0.5;
            self.z = (mat.value(2, 3) + mat.value(3, 2)) * inv;
            self.w = (mat.value(1, 3) - mat.value(3, 1)) * inv;
        } else {
            let s = (mat.value(3, 3) - mat.value(1, 1) - mat.value(2, 2) + 1.0).sqrt();
            let inv = 0.5 / s;
            self.x = (mat.value(1, 3) + mat.value(3, 1)) * inv;
            self.y = (mat.value(2, 3) + mat.value(3, 2)) * inv;
            self.z = s * 0.5;
            self.w = (mat.value(2, 1) - mat.value(1, 2)) * inv;
        }
    }

    /// Convert quaternion to 3x3 rotation matrix.
    pub fn get_matrix(&self) -> GpMat {
        let xx = self.x * self.x;
        let yy = self.y * self.y;
        let zz = self.z * self.z;
        let xy = self.x * self.y;
        let xz = self.x * self.z;
        let yz = self.y * self.z;
        let wx = self.w * self.x;
        let wy = self.w * self.y;
        let wz = self.w * self.z;

        let mut m = GpMat::zero();
        // ponytail: build 3x3 rotation matrix directly
        m.m[0][0] = 1.0 - 2.0 * (yy + zz);
        m.m[0][1] = 2.0 * (xy - wz);
        m.m[0][2] = 2.0 * (xz + wy);
        m.m[1][0] = 2.0 * (xy + wz);
        m.m[1][1] = 1.0 - 2.0 * (xx + zz);
        m.m[1][2] = 2.0 * (yz - wx);
        m.m[2][0] = 2.0 * (xz - wy);
        m.m[2][1] = 2.0 * (yz + wx);
        m.m[2][2] = 1.0 - 2.0 * (xx + yy);
        m
    }

    pub fn square_norm(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w
    }

    pub fn norm(&self) -> f64 {
        self.square_norm().sqrt()
    }

    /// Normalize in-place. Returns Err on zero norm.
    pub fn normalize(&mut self) -> Result<(), &'static str> {
        let n = self.norm();
        if n <= RESOLUTION {
            return Err("GpQuaternion::normalize: zero norm");
        }
        let inv = 1.0 / n;
        self.x *= inv;
        self.y *= inv;
        self.z *= inv;
        self.w *= inv;
        Ok(())
    }

    pub fn normalized(&self) -> Result<Self, &'static str> {
        let n = self.norm();
        if n <= RESOLUTION {
            return Err("GpQuaternion::normalized: zero norm");
        }
        let inv = 1.0 / n;
        Ok(Self {
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
            w: self.w * inv,
        })
    }

    pub fn set_identity(&mut self) {
        self.x = 0.0;
        self.y = 0.0;
        self.z = 0.0;
        self.w = 1.0;
    }

    pub fn is_identity(&self) -> bool {
        self.x.abs() <= RESOLUTION
            && self.y.abs() <= RESOLUTION
            && self.z.abs() <= RESOLUTION
            && (self.w - 1.0).abs() <= RESOLUTION
    }

    pub fn reverse(&mut self) {
        self.invert();
    }

    pub fn reversed(&self) -> Self {
        self.inverted()
    }

    pub fn invert(&mut self) {
        let n = self.square_norm();
        if n > RESOLUTION * RESOLUTION {
            let inv = 1.0 / n;
            self.x = -self.x * inv;
            self.y = -self.y * inv;
            self.z = -self.z * inv;
            self.w = self.w * inv;
        }
    }

    pub fn inverted(&self) -> Self {
        let n = self.square_norm();
        if n > RESOLUTION * RESOLUTION {
            let inv = 1.0 / n;
            Self {
                x: -self.x * inv,
                y: -self.y * inv,
                z: -self.z * inv,
                w: self.w * inv,
            }
        } else {
            *self
        }
    }

    pub fn negated(&self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        }
    }

    /// self *= other (consuming self)
    pub fn multiply(mut self, other: &GpQuaternion) -> Self {
        let qx = self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y;
        let qy = self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x;
        let qz = self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w;
        let qw = self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z;
        self.x = qx;
        self.y = qy;
        self.z = qz;
        self.w = qw;
        self
    }

    /// self * other (returns new)
    pub fn multiplied(&self, other: &GpQuaternion) -> Self {
        Self {
            x: self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            y: self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
            z: self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w,
            w: self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
        }
    }

    /// self += other (consuming self)
    pub fn add(mut self, other: &GpQuaternion) -> Self {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
        self.w += other.w;
        self
    }

    /// self + other (returns new)
    pub fn added(&self, other: &GpQuaternion) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            w: self.w + other.w,
        }
    }

    /// self -= other (consuming self)
    pub fn subtract(mut self, other: &GpQuaternion) -> Self {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
        self.w -= other.w;
        self
    }

    /// self - other (returns new)
    pub fn subtracted(&self, other: &GpQuaternion) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            w: self.w - other.w,
        }
    }

    /// self *= scalar (consuming self)
    pub fn scale(mut self, s: f64) -> Self {
        self.x *= s;
        self.y *= s;
        self.z *= s;
        self.w *= s;
        self
    }

    /// self * scalar (returns new)
    pub fn scaled(&self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
            w: self.w * s,
        }
    }

    pub fn is_equal(&self, other: &GpQuaternion) -> bool {
        let d = self.subtracted(other);
        d.square_norm() <= RESOLUTION * RESOLUTION
    }
}

impl std::ops::Mul<&GpQuaternion> for GpQuaternion {
    type Output = GpQuaternion;

    fn mul(self, rhs: &GpQuaternion) -> GpQuaternion {
        self.multiplied(rhs)
    }
}

impl std::ops::Mul<&GpQuaternion> for &GpQuaternion {
    type Output = GpQuaternion;

    fn mul(self, rhs: &GpQuaternion) -> GpQuaternion {
        self.multiplied(rhs)
    }
}
