use crate::gp::vec::GpVec;
use crate::gp::xyz::GpXyz;
use crate::precision::RESOLUTION;

/// Axis direction constants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirAxis {
    X,
    Y,
    Z,
    NX,
    NY,
    NZ,
}

/// Unit direction (always normalized).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpDir {
    coord: GpXyz,
}

impl Default for GpDir {
    fn default() -> Self {
        Self::default_dir()
    }
}

impl GpDir {
    /// Default direction: X axis.
    pub fn default_dir() -> Self {
        Self {
            coord: GpXyz::new(1.0, 0.0, 0.0),
        }
    }

    /// Create from a named axis.
    pub fn from_axis(axis: DirAxis) -> Self {
        match axis {
            DirAxis::X => Self {
                coord: GpXyz::new(1.0, 0.0, 0.0),
            },
            DirAxis::Y => Self {
                coord: GpXyz::new(0.0, 1.0, 0.0),
            },
            DirAxis::Z => Self {
                coord: GpXyz::new(0.0, 0.0, 1.0),
            },
            DirAxis::NX => Self {
                coord: GpXyz::new(-1.0, 0.0, 0.0),
            },
            DirAxis::NY => Self {
                coord: GpXyz::new(0.0, -1.0, 0.0),
            },
            DirAxis::NZ => Self {
                coord: GpXyz::new(0.0, 0.0, -1.0),
            },
        }
    }

    /// Create from components. Normalizes iff needed. Returns Err on zero.
    pub fn new(x: f64, y: f64, z: f64) -> Result<Self, &'static str> {
        let m = (x * x + y * y + z * z).sqrt();
        if m <= RESOLUTION {
            return Err("GpDir: zero vector");
        }
        // fast-path: already unit within tolerance
        if (m - 1.0).abs() <= RESOLUTION {
            return Ok(Self {
                coord: GpXyz::new(x, y, z),
            });
        }
        let inv = 1.0 / m;
        Ok(Self {
            coord: GpXyz::new(x * inv, y * inv, z * inv),
        })
    }

    pub fn from_xyz(xyz: &GpXyz) -> Result<Self, &'static str> {
        Self::new(xyz.x(), xyz.y(), xyz.z())
    }

    /// Create from vector (normalizes). Returns Err on zero.
    pub fn from_vec(v: &GpVec) -> Result<Self, &'static str> {
        Self::new(v.x(), v.y(), v.z())
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

    pub fn is_equal(&self, other: &GpDir) -> bool {
        let diff = GpXyz::new(
            self.coord.x() - other.coord.x(),
            self.coord.y() - other.coord.y(),
            self.coord.z() - other.coord.z(),
        );
        diff.square_modulus() <= RESOLUTION * RESOLUTION
    }

    pub fn dot(&self, other: &GpDir) -> f64 {
        self.coord.dot(&other.coord)
    }

    /// self = self ^ other (re-normalizes result). Err on zero cross.
    pub fn cross(mut self, other: &GpDir) -> Result<Self, &'static str> {
        let crossed = self.coord.crossed(&other.coord);
        let m = crossed.modulus();
        if m <= RESOLUTION {
            return Err("GpDir::cross: zero cross product");
        }
        self.coord = crossed.divided(m);
        Ok(self)
    }

    /// self ^ other (re-normalizes result). Err on zero cross.
    pub fn crossed(&self, other: &GpDir) -> Result<Self, &'static str> {
        let crossed = self.coord.crossed(&other.coord);
        let m = crossed.modulus();
        if m <= RESOLUTION {
            return Err("GpDir::crossed: zero cross product");
        }
        Ok(Self {
            coord: crossed.divided(m),
        })
    }

    /// (self ^ other) ^ self (re-normalizes result). Err on zero cross.
    pub fn cross_cross(&self, other: &GpDir) -> Result<Self, &'static str> {
        let co = self.crossed(other)?;
        let cc = self.coord.crossed(&co.coord);
        let m = cc.modulus();
        if m <= RESOLUTION {
            return Err("GpDir::cross_cross: zero result");
        }
        Ok(Self {
            coord: cc.divided(m),
        })
    }

    pub fn dot_cross(&self, v1: &GpDir, v2: &GpDir) -> f64 {
        let c = v1.coord.crossed(&v2.coord);
        self.coord.dot(&c)
    }

    pub fn reverse(&mut self) {
        self.coord = self.coord.multiplied(-1.0);
    }

    pub fn reversed(&self) -> Self {
        Self {
            coord: self.coord.multiplied(-1.0),
        }
    }

    /// Angle between self and other (radians). Returns 0 if either is zero.
    pub fn angle(&self, other: &GpDir) -> f64 {
        let cos = self.dot(other).clamp(-1.0, 1.0);
        cos.acos()
    }

    /// Angle between self and other, with reference direction Vref.
    /// Positive if cross product aligns with Vref, negative otherwise.
    pub fn angle_with_ref(&self, other: &GpDir, vref: &GpDir) -> f64 {
        if let Ok(c) = self.crossed(other) {
            let a = self.angle(other);
            if c.dot(vref) < 0.0 {
                return -a;
            }
            a
        } else {
            0.0
        }
    }

    pub fn is_normal(&self, other: &GpDir) -> bool {
        self.dot(other).abs() <= RESOLUTION
    }

    pub fn is_opposite(&self, other: &GpDir) -> bool {
        self.dot(other) + 1.0 <= RESOLUTION
    }

    pub fn is_parallel(&self, other: &GpDir) -> bool {
        self.dot(other).abs() - 1.0 <= RESOLUTION
    }

    #[cfg(feature = "serde_json")]
    pub fn dump_json(&self) -> serde_json::Value {
        serde_json::json!({
            "x": self.coord.x(),
            "y": self.coord.y(),
            "z": self.coord.z(),
        })
    }
}
