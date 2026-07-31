use crate::gp::mat::GpMat;
use crate::gp::pnt::GpPnt;
use crate::gp::xyz::GpXyz;
use crate::precision::RESOLUTION;

/// 3D vector newtype over GpXyz.
#[derive(Debug, Clone, Copy)]
pub struct GpVec {
    pub coord: GpXyz,
}

impl Default for GpVec {
    fn default() -> Self {
        Self::zero()
    }
}

impl GpVec {
    /// Zero vector (0,0,0).
    pub fn zero() -> Self {
        Self {
            coord: GpXyz::new(0.0, 0.0, 0.0),
        }
    }

    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self {
            coord: GpXyz::new(x, y, z),
        }
    }

    pub fn from_xyz(xyz: &GpXyz) -> Self {
        Self { coord: *xyz }
    }

    pub fn from_pnts(p1: &GpPnt, p2: &GpPnt) -> Self {
        Self {
            coord: GpXyz::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z()),
        }
    }

    pub fn x(&self) -> f64 {
        self.coord.x()
    }

    pub fn y(&self) -> f64 {
        self.coord.y()
    }

    pub fn z(&self) -> f64 {
        self.coord.z()
    }

    pub fn xyz(&self) -> &GpXyz {
        &self.coord
    }

    pub fn magnitude(&self) -> f64 {
        self.coord.modulus()
    }

    pub fn square_magnitude(&self) -> f64 {
        self.coord.square_modulus()
    }

    pub fn is_equal(&self, other: &GpVec) -> bool {
        let diff = self.subtracted(other);
        diff.square_magnitude() <= RESOLUTION * RESOLUTION
    }

    /// self + other (consuming self)
    pub fn add(mut self, other: &GpVec) -> Self {
        self.coord = self.coord.added(&other.coord);
        self
    }

    /// self + other (returns new)
    pub fn added(&self, other: &GpVec) -> Self {
        Self {
            coord: self.coord.added(&other.coord),
        }
    }

    /// self - other (consuming self)
    pub fn subtract(mut self, other: &GpVec) -> Self {
        self.coord = self.coord.subtracted(&other.coord);
        self
    }

    /// self - other (returns new)
    pub fn subtracted(&self, other: &GpVec) -> Self {
        Self {
            coord: self.coord.subtracted(&other.coord),
        }
    }

    /// self *= scalar (consuming self)
    pub fn multiply_scalar(mut self, scalar: f64) -> Self {
        self.coord = self.coord.multiplied(scalar);
        self
    }

    /// self * scalar (returns new)
    pub fn multiplied_scalar(&self, scalar: f64) -> Self {
        Self {
            coord: self.coord.multiplied(scalar),
        }
    }

    /// self /= scalar (consuming self)
    pub fn divide(mut self, scalar: f64) -> Self {
        self.coord = self.coord.divided(scalar);
        self
    }

    /// self / scalar (returns new)
    pub fn divided(&self, scalar: f64) -> Self {
        Self {
            coord: self.coord.divided(scalar),
        }
    }

    pub fn dot(&self, other: &GpVec) -> f64 {
        self.coord.dot(&other.coord)
    }

    /// self = self ^ other (consuming self)
    pub fn cross(mut self, other: &GpVec) -> Self {
        self.coord = self.coord.crossed(&other.coord);
        self
    }

    /// self ^ other (returns new)
    pub fn crossed(&self, other: &GpVec) -> Self {
        Self {
            coord: self.coord.crossed(&other.coord),
        }
    }

    pub fn cross_magnitude(&self, other: &GpVec) -> f64 {
        let c = self.crossed(other);
        c.magnitude()
    }

    pub fn cross_square_magnitude(&self, other: &GpVec) -> f64 {
        let c = self.crossed(other);
        c.square_magnitude()
    }

    /// (self ^ other) ^ self (consuming self)
    pub fn cross_cross(mut self, other: &GpVec) -> Self {
        let crossed = self.coord.crossed(&other.coord);
        self.coord = crossed.crossed(&self.coord);
        self
    }

    /// (self ^ other) ^ self (returns new)
    pub fn cross_crossed(&self, other: &GpVec) -> Self {
        let crossed = self.coord.crossed(&other.coord);
        Self {
            coord: crossed.crossed(&self.coord),
        }
    }

    pub fn dot_cross(&self, v1: &GpVec, v2: &GpVec) -> f64 {
        let c = v1.crossed(v2);
        self.dot(&c)
    }

    /// Normalize in-place. Returns magnitude, or 0 if zero vector.
    pub fn normalize(&mut self) -> f64 {
        let m = self.magnitude();
        if m > RESOLUTION {
            self.coord = self.coord.divided(m);
        }
        m
    }

    pub fn normalized(&self) -> Self {
        let m = self.magnitude();
        if m > RESOLUTION {
            Self {
                coord: self.coord.divided(m),
            }
        } else {
            *self
        }
    }

    /// Reverse in-place.
    pub fn reverse(&mut self) {
        self.coord = self.coord.multiplied(-1.0);
    }

    pub fn reversed(&self) -> Self {
        Self {
            coord: self.coord.multiplied(-1.0),
        }
    }

    /// Angle between self and other (radians). Returns 0 if either is zero.
    pub fn angle(&self, other: &GpVec) -> f64 {
        let m1 = self.magnitude();
        let m2 = other.magnitude();
        if m1 <= RESOLUTION || m2 <= RESOLUTION {
            return 0.0;
        }
        let cos = self.dot(other) / (m1 * m2);
        // clamp for numerical stability
        let cos = cos.clamp(-1.0, 1.0);
        cos.acos()
    }

    pub fn is_normal(&self, other: &GpVec) -> bool {
        self.dot(other).abs() <= RESOLUTION
    }

    pub fn is_parallel(&self, other: &GpVec) -> bool {
        let m = self.cross_magnitude(other);
        m <= RESOLUTION
    }

    pub fn is_opposite(&self, other: &GpVec) -> bool {
        let m = self.cross_magnitude(other);
        if m > RESOLUTION {
            return false;
        }
        self.dot(other) < 0.0
    }

    pub fn multiply_mat(mut self, mat: &GpMat) -> Self {
        self.coord = self.coord.multiplied_mat(mat);
        self
    }

    pub fn multiplied_mat(&self, mat: &GpMat) -> Self {
        Self {
            coord: self.coord.multiplied_mat(mat),
        }
    }

    /// self = a1*v1 + a2*v2 (consuming self)
    pub fn set_linear_form_2(mut self, a1: f64, v1: &GpVec, a2: f64, v2: &GpVec) -> Self {
        self.coord = v1.coord.multiplied(a1).added(&v2.coord.multiplied(a2));
        self
    }

    /// self = a1*v1 + a2*v2 + a3*v3 (consuming self)
    pub fn set_linear_form_3(
        mut self,
        a1: f64,
        v1: &GpVec,
        a2: f64,
        v2: &GpVec,
        a3: f64,
        v3: &GpVec,
    ) -> Self {
        self.coord = v1
            .coord
            .multiplied(a1)
            .added(&v2.coord.multiplied(a2))
            .added(&v3.coord.multiplied(a3));
        self
    }

    /// self = a1*v1 + a2*v2 + v3 (consuming self)
    pub fn set_linear_form_2plus(
        mut self,
        a1: f64,
        v1: &GpVec,
        a2: f64,
        v2: &GpVec,
        v3: &GpVec,
    ) -> Self {
        self.coord = v1
            .coord
            .multiplied(a1)
            .added(&v2.coord.multiplied(a2))
            .added(&v3.coord);
        self
    }
}
