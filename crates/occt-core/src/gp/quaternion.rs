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

    /// `gp_Quaternion::SetVectorAndAngle(const gp_Vec&, const double)`
    /// (`gp_Quaternion.cxx:81-88`): normalize the axis, then half-angle.
    pub fn set_vector_and_angle(&mut self, axis: &GpVec, angle: f64) {
        let a = axis.normalized();
        let half = 0.5 * angle;
        let sin_a = half.sin();
        self.x = a.x() * sin_a;
        self.y = a.y() * sin_a;
        self.z = a.z() * sin_a;
        self.w = half.cos();
    }

    /// `gp_Quaternion::GetVectorAndAngle(gp_Vec&, double&)`
    /// (`gp_Quaternion.cxx:91-112`). Returns `(axis, angle)` with the angle in
    /// `[-PI, PI]`; a zero vector part yields axis `+Z` and angle `0`.
    pub fn get_vector_and_angle(&self) -> (GpVec, f64) {
        let vl = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if vl > RESOLUTION {
            let ivl = 1.0 / vl;
            let axis = GpVec::new(self.x * ivl, self.y * ivl, self.z * ivl);
            let angle = if self.w < 0.0 {
                2.0 * (-vl).atan2(-self.w)
            } else {
                2.0 * vl.atan2(self.w)
            };
            (axis, angle)
        } else {
            (GpVec::new(0.0, 0.0, 1.0), 0.0)
        }
    }

    /// `gp_Quaternion::SetEulerAngles(...)` (`gp_Quaternion.cxx:302-357`).
    pub fn set_euler_angles(
        &mut self,
        order: crate::gp::euler_sequence::GpEulerSequence,
        alpha: f64,
        beta: f64,
        gamma: f64,
    ) {
        let o = crate::gp::euler_sequence::translate_euler_sequence(order);

        let mut a = alpha;
        let mut b = beta;
        let mut c = gamma;
        if !o.is_extrinsic {
            a = gamma;
            c = alpha;
        }
        if o.is_odd {
            b = -b;
        }

        let ti = 0.5 * a;
        let tj = 0.5 * b;
        let th = 0.5 * c;
        let ci = ti.cos();
        let cj = tj.cos();
        let ch = th.cos();
        let si = ti.sin();
        let sj = tj.sin();
        let sh = th.sin();
        let cc = ci * ch;
        let cs = ci * sh;
        let sc = si * ch;
        let ss = si * sh;

        // values[0] = w, values[1..3] = x, y, z (1-based axis indices).
        let mut values = [0.0f64; 4];
        if o.is_two_axes {
            values[o.i] = cj * (cs + sc);
            values[o.j] = sj * (cc + ss);
            values[o.k] = sj * (cs - sc);
            values[0] = cj * (cc - ss);
        } else {
            values[o.i] = cj * sc - sj * cs;
            values[o.j] = cj * ss + sj * cc;
            values[o.k] = cj * cs - sj * sc;
            values[0] = cj * cc + sj * ss;
        }
        if o.is_odd {
            values[o.j] = -values[o.j];
        }

        self.x = values[1];
        self.y = values[2];
        self.z = values[3];
        self.w = values[0];
    }

    /// `gp_Quaternion::GetEulerAngles(...)` (`gp_Quaternion.cxx:362-412`).
    pub fn get_euler_angles(
        &self,
        order: crate::gp::euler_sequence::GpEulerSequence,
    ) -> (f64, f64, f64) {
        let m = self.get_matrix();
        let o = crate::gp::euler_sequence::translate_euler_sequence(order);

        let (mut alpha, mut beta, mut gamma);
        if o.is_two_axes {
            let sy = (m.value(o.i, o.j) * m.value(o.i, o.j)
                + m.value(o.i, o.k) * m.value(o.i, o.k))
            .sqrt();
            if sy > 16.0 * f64::EPSILON {
                alpha = m.value(o.i, o.j).atan2(m.value(o.i, o.k));
                gamma = m.value(o.j, o.i).atan2(-m.value(o.k, o.i));
            } else {
                alpha = (-m.value(o.j, o.k)).atan2(m.value(o.j, o.j));
                gamma = 0.0;
            }
            beta = sy.atan2(m.value(o.i, o.i));
        } else {
            let cy = (m.value(o.i, o.i) * m.value(o.i, o.i)
                + m.value(o.j, o.i) * m.value(o.j, o.i))
            .sqrt();
            if cy > 16.0 * f64::EPSILON {
                alpha = m.value(o.k, o.j).atan2(m.value(o.k, o.k));
                gamma = m.value(o.j, o.i).atan2(m.value(o.i, o.i));
            } else {
                alpha = (-m.value(o.j, o.k)).atan2(m.value(o.j, o.j));
                gamma = 0.0;
            }
            beta = (-m.value(o.k, o.i)).atan2(cy);
        }
        if o.is_odd {
            alpha = -alpha;
            beta = -beta;
            gamma = -gamma;
        }
        if !o.is_extrinsic {
            let a_first = alpha;
            alpha = gamma;
            gamma = a_first;
        }
        (alpha, beta, gamma)
    }

    /// `gp_Quaternion::StabilizeLength()` (`gp_Quaternion.cxx:416-430`).
    pub fn stabilize_length(&mut self) {
        let cs = self.x.abs() + self.y.abs() + self.z.abs() + self.w.abs();
        if cs > 0.0 {
            self.x /= cs;
            self.y /= cs;
            self.z /= cs;
            self.w /= cs;
        } else {
            self.set_identity();
        }
    }

    /// `gp_Quaternion::Normalize()` (`gp_Quaternion.cxx:434-444`): divide by
    /// `Norm()`, or stabilize the length when degenerate.
    pub fn normalize_or_stabilize(&mut self) {
        let magn = self.norm();
        if magn < RESOLUTION {
            self.stabilize_length();
        } else {
            let inv = 1.0 / magn;
            self.x *= inv;
            self.y *= inv;
            self.z *= inv;
            self.w *= inv;
        }
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

    /// `gp_Quaternion::GetMatrix()` (`gp_Quaternion.cxx:153-187`): the factor
    /// `2/SquareNorm()` means a non-unit quaternion yields a scaled matrix.
    pub fn get_matrix(&self) -> GpMat {
        let s = 2.0 / self.square_norm();
        let x2 = self.x * s;
        let y2 = self.y * s;
        let z2 = self.z * s;
        let xx = self.x * x2;
        let xy = self.x * y2;
        let xz = self.x * z2;
        let yy = self.y * y2;
        let yz = self.y * z2;
        let zz = self.z * z2;
        let wx = self.w * x2;
        let wy = self.w * y2;
        let wz = self.w * z2;

        GpMat::new(
            1.0 - (yy + zz), xy - wz, xz + wy,
            xy + wz, 1.0 - (xx + zz), yz - wx,
            xz - wy, yz + wx, 1.0 - (xx + yy),
        )
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
