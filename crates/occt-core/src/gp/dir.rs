use crate::gp::vec::GpVec;
use crate::gp::xyz::GpXyz;
use crate::precision::{REAL_SMALL, RESOLUTION};

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

    /// `gp_Dir::gp_Dir(double, double, double)` (`gp_Dir.hxx:352-369`).
    ///
    /// The predicate is on the **square modulus**: OCCT's fast path compares
    /// `aSqMod` with `1.0 +- gp::Resolution()`, and `gp::Resolution()` is
    /// `RealSmall()` == `DBL_MIN` (`gp.hxx:60`), so both bounds round to `1.0`
    /// and only an already-unit square modulus skips the division. The zero
    /// guard is `aSqMod <= gp::Resolution() * gp::Resolution()`, and
    /// `DBL_MIN * DBL_MIN` underflows to `0.0`, so only an exactly zero vector
    /// raises. The previous body compared `sqrt(aSqMod)` against
    /// `Precision::Resolution()` (1e-12) and kept any vector within 1e-12 of
    /// unit length unnormalized.
    pub fn new(x: f64, y: f64, z: f64) -> Result<Self, &'static str> {
        let a_sq_mod = x * x + y * y + z * z;
        // Fast path: already normalized (`gp_Dir.hxx:357-361`).
        if a_sq_mod >= 1.0 - REAL_SMALL && a_sq_mod <= 1.0 + REAL_SMALL {
            return Ok(Self {
                coord: GpXyz::new(x, y, z),
            });
        }
        if a_sq_mod <= REAL_SMALL * REAL_SMALL {
            return Err("GpDir: zero vector");
        }
        let inv = 1.0 / a_sq_mod.sqrt();
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
    ///
    /// OCCT counterpart: `gp_Dir::Crossed(const gp_Dir&)`
    /// (`gp_Dir.cxx:299-302`) delegates to the `gp_Dir(gp_XYZ)` constructor,
    /// whose guard is `aSqMod <= gp::Resolution() * gp::Resolution()`
    /// (`gp_Dir.hxx:364`).
    pub fn cross(mut self, other: &GpDir) -> Result<Self, &'static str> {
        let crossed = self.coord.crossed(&other.coord);
        let m = crossed.modulus();
        if m <= REAL_SMALL {
            return Err("GpDir::cross: zero cross product");
        }
        self.coord = crossed.divided(m);
        Ok(self)
    }

    /// self ^ other (re-normalizes result). Err on zero cross.
    pub fn crossed(&self, other: &GpDir) -> Result<Self, &'static str> {
        let crossed = self.coord.crossed(&other.coord);
        let m = crossed.modulus();
        if m <= REAL_SMALL {
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
        if m <= REAL_SMALL {
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

    /// `gp_Dir::Rotate(const gp_Ax1&, angle)` (`gp_Dir.hxx:496-501`): multiply
    /// by the rotation's h-vectorial part (no renormalization).
    pub fn rotate(&mut self, a1: &crate::gp::ax1::GpAx1, ang: f64) {
        let mut t = crate::gp::trsf::GpTrsf::identity();
        if t.set_rotation_ax1(a1, ang).is_err() {
            return;
        }
        let mut xyz = self.coord;
        xyz.multiply_mat(&t.h_vectorial_part());
        self.coord = xyz;
    }
    /// `gp_Dir::Mirror(const gp_Dir&)` (`gp_Dir.cxx:86-102`): reflect the
    /// direction across the line carried by `v`.
    pub fn mirror_dir(&mut self, v: &Self) {
        let (a, b, c) = (v.coord.x(), v.coord.y(), v.coord.z());
        let (x, y, z) = (self.coord.x(), self.coord.y(), self.coord.z());
        let m1 = 2.0 * a * b;
        let m2 = 2.0 * a * c;
        let m3 = 2.0 * b * c;
        self.coord = GpXyz::new(
            (2.0 * a * a - 1.0) * x + m1 * y + m2 * z,
            m1 * x + (2.0 * b * b - 1.0) * y + m3 * z,
            m2 * x + m3 * y + (2.0 * c * c - 1.0) * z,
        );
    }
    /// `gp_Dir::Mirror(const gp_Ax1&)` (`gp_Dir.cxx:104-120`). Ported verbatim,
    /// including OCCT's line 109 reading `C = XYZ.Y()` rather than `Z`.
    pub fn mirror_ax1(&mut self, a1: &crate::gp::ax1::GpAx1) {
        let d = a1.direction();
        let (a, b) = (d.x(), d.y());
        let c = d.y();
        let (x, y, z) = (self.coord.x(), self.coord.y(), self.coord.z());
        let m1 = 2.0 * a * b;
        let m2 = 2.0 * a * c;
        let m3 = 2.0 * b * c;
        self.coord = GpXyz::new(
            (2.0 * a * a - 1.0) * x + m1 * y + m2 * z,
            m1 * x + (2.0 * b * b - 1.0) * y + m3 * z,
            m2 * x + m3 * y + (2.0 * c * c - 1.0) * z,
        );
    }
    /// `gp_Dir::Mirror(const gp_Ax2&)` (`gp_Dir.cxx:122-127`): mirror across
    /// the axis placement's main direction, then reverse.
    pub fn mirror_ax2(&mut self, a2: &crate::gp::ax2::GpAx2) {
        let vz = a2.direction();
        self.mirror_dir(&vz);
        self.reverse();
    }
    /// `gp_Dir::Transform(const gp_Trsf&)` (`gp_Dir.cxx:129-155`).
    pub fn transform(&mut self, t: &crate::gp::trsf::GpTrsf) {
        use crate::gp::trsf_form::TrsfForm;
        match t.form() {
            TrsfForm::Identity | TrsfForm::Translation => {}
            TrsfForm::PntMirror => self.reverse(),
            TrsfForm::Scale => {
                if t.scale_factor() < 0.0 {
                    self.reverse();
                }
            }
            _ => {
                let mut xyz = self.coord;
                xyz.multiply_mat(&t.h_vectorial_part());
                let d = xyz.modulus();
                // `gp_Dir::Transform` (`gp_Dir.cxx:54-70`): normalizes only
                // when `aMagnitude > gp::Resolution()`.
                if d > REAL_SMALL {
                    xyz = xyz.divided(d);
                }
                if t.scale_factor() < 0.0 {
                    xyz.reverse();
                }
                self.coord = xyz;
            }
        }
    }

    /// Angle between self and other (radians). Returns 0 if either is zero.
    pub fn angle(&self, other: &GpDir) -> f64 {
        let cos = self.dot(other).clamp(-1.0, 1.0);
        cos.acos()
    }

    /// `gp_Dir::Angle` (`gp_Dir.cxx:27-52`) - accurate near 0 and PI.
    pub fn angle_tol(&self, other: &GpDir) -> f64 {
        let cosinus = self.dot(other);
        if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
            cosinus.clamp(-1.0, 1.0).acos()
        } else {
            let sinus = self.coord.crossed(&other.coord).modulus();
            if cosinus < 0.0 {
                std::f64::consts::PI - sinus.asin()
            } else {
                sinus.asin()
            }
        }
    }

    /// `gp_Dir::IsParallel(other, theAngularTolerance)` (`gp_Dir.hxx`).
    pub fn is_parallel_tol(&self, other: &GpDir, ang_tol: f64) -> bool {
        let an_ang = self.angle_tol(other);
        an_ang <= ang_tol || (std::f64::consts::PI - an_ang) <= ang_tol
    }

    /// `gp_Dir::AngleWithRef` (`gp_Dir.cxx:55-84`). Opposite directions
    /// use `PI - asin(|cross|)` so the result is `PI`, not 0.
    pub fn angle_with_ref(&self, other: &GpDir, vref: &GpDir) -> f64 {
        let xyz = self.coord.crossed(&other.coord);
        let cosinus = self.coord.dot(&other.coord);
        let sinus = xyz.modulus();
        let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
            cosinus.acos()
        } else if cosinus < 0.0 {
            std::f64::consts::PI - sinus.asin()
        } else {
            sinus.asin()
        };
        if xyz.dot(vref.xyz()) >= 0.0 {
            ang
        } else {
            -ang
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
