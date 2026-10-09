use std::ops::Mul;

use crate::gp::ax1::GpAx1;
use crate::gp::ax2::GpAx2;
use crate::gp::ax3::GpAx3;
use crate::gp::mat::GpMat;
use crate::gp::pnt::GpPnt;
use crate::gp::quaternion::GpQuaternion;
use crate::gp::trsf_form::TrsfForm;
use crate::gp::vec::GpVec;
use crate::gp::xyz::GpXyz;
use crate::precision::RESOLUTION;

#[derive(Debug, Clone)]
pub struct GpTrsf {
    pub scale: f64,
    pub shape: TrsfForm,
    pub matrix: GpMat,
    pub loc: GpXyz,
}

impl Default for GpTrsf {
    fn default() -> Self {
        Self::identity()
    }
}

impl GpTrsf {
    pub fn identity() -> Self {
        GpTrsf {
            scale: 1.0,
            shape: TrsfForm::Identity,
            matrix: GpMat::identity(),
            loc: GpXyz::zero(),
        }
    }

    pub fn scale_factor(&self) -> f64 {
        self.scale
    }

    pub fn form(&self) -> TrsfForm {
        self.shape
    }

    pub fn translation_part(&self) -> GpXyz {
        self.loc.clone()
    }

    /// `gp_Trsf::HVectorialPart()` (`gp_Trsf.hxx:254`): the 3x3 matrix without
    /// the scale factor.
    pub fn h_vectorial_part(&self) -> GpMat {
        self.matrix
    }

    /// `gp_Trsf::VectorialPart()` (`gp_Trsf.hxx:464-482`): the 3x3 matrix with
    /// the scale factor folded in. For `Scale`/`PntMirror` only the diagonal is
    /// scaled; otherwise the whole matrix is multiplied.
    pub fn vectorial_part(&self) -> GpMat {
        if self.scale == 1.0 {
            return self.matrix;
        }
        let mut m = self.matrix;
        if self.shape == TrsfForm::Scale || self.shape == TrsfForm::PntMirror {
            m.set_diagonal(&GpXyz::new(
                self.scale * self.matrix.m[0][0],
                self.scale * self.matrix.m[1][1],
                self.scale * self.matrix.m[2][2],
            ));
        } else {
            m = m.multiply_scalar(self.scale);
        }
        m
    }

    pub fn is_negative(&self) -> bool {
        self.matrix.determinant() < 0.0
    }

    pub fn set_identity(&mut self) {
        self.scale = 1.0;
        self.shape = TrsfForm::Identity;
        self.matrix = GpMat::identity();
        self.loc = GpXyz::zero();
    }

    /// `gp_Trsf::SetTransformation(const gp_Ax3&)` — world to local frame of `a3`.
    pub fn set_transformation(&mut self, a3: &GpAx3) {
        self.shape = TrsfForm::CompoundTrsf;
        self.scale = 1.0;
        self.matrix.set_rows(
            a3.x_direction().xyz(),
            a3.y_direction().xyz(),
            a3.direction().xyz(),
        );
        let mut loc = a3.location().coord;
        loc.multiply_mat(&self.matrix);
        loc.reverse();
        self.loc = loc;
    }

    /// `gp_Trsf::SetTransformation(const gp_Ax3& theFromSystem1, const gp_Ax3& theToSystem2)`
    /// (`gp_Trsf.cxx:172-194`): the change of basis that maps a point expressed
    /// in `from_system1` to its coordinates in `to_system2`. The 1-arg overload
    /// (`gp_Trsf.cxx:196-204`) is [`GpTrsf::set_transformation`], which is the
    /// special case `to_system2 == gp::XOY()`.
    pub fn set_transformation_from_to(&mut self, from_system1: &GpAx3, to_system2: &GpAx3) {
        self.shape = TrsfForm::CompoundTrsf;
        self.scale = 1.0;
        // matrix from XOY to `to_system2` (`cxx:176-180`):
        self.matrix.set_rows(
            to_system2.x_direction().xyz(),
            to_system2.y_direction().xyz(),
            to_system2.direction().xyz(),
        );
        let mut loc = to_system2.location().coord;
        loc.multiply_mat(&self.matrix);
        loc.reverse();
        // matrix from `from_system1` to XOY (`cxx:182-193`). OCCT builds it with
        // `gp_Mat MA1(xDir, yDir, zDir)`, whose constructor stores the three
        // arguments as *columns* (`gp_Mat.cxx:30-43`), not rows.
        let ma1 = GpMat::from_cols(
            from_system1.x_direction().xyz(),
            from_system1.y_direction().xyz(),
            from_system1.direction().xyz(),
        );
        let mut ma1_loc = from_system1.location().coord;
        ma1_loc.multiply_mat(&self.matrix);
        loc = loc.add(&ma1_loc);
        self.matrix = self.matrix.multiply(&ma1);
        self.loc = loc;
    }

    /// `gp_Trsf::SetDisplacement(const gp_Ax3& theFromSystem1, const gp_Ax3& theToSystem2)`
    /// (`gp_Trsf.cxx:217-240`): the placement that carries geometry from
    /// `from_system1` into `to_system2`. Unlike [`GpTrsf::set_transformation_from_to`],
    /// the placement axes enter as *columns* (`SetCol`) and `MA1` is built with
    /// the column constructor then transposed.
    pub fn set_displacement(&mut self, from_system1: &GpAx3, to_system2: &GpAx3) {
        self.shape = TrsfForm::CompoundTrsf;
        self.scale = 1.0;
        // matrix from `to_system2` to XOY (`cxx:221-225`):
        self.matrix = GpMat::from_cols(
            to_system2.x_direction().xyz(),
            to_system2.y_direction().xyz(),
            to_system2.direction().xyz(),
        );
        self.loc = to_system2.location().coord;
        // matrix XOY to `from_system1` (`cxx:226-239`):
        let ma1 = GpMat::from_cols(
            from_system1.x_direction().xyz(),
            from_system1.y_direction().xyz(),
            from_system1.direction().xyz(),
        )
        .transpose();
        let mut ma1_loc = from_system1.location().coord;
        ma1_loc.multiply_mat(&ma1);
        ma1_loc.reverse();
        ma1_loc.multiply_mat(&self.matrix);
        self.loc = self.loc.add(&ma1_loc);
        self.matrix = self.matrix.multiply(&ma1);
    }

    /// `gp_Trsf::GetRotation()` (`gp_Trsf.cxx:387-390`): the rotation part as a
    /// quaternion, built from the h-vectorial matrix (scale excluded).
    pub fn get_rotation(&self) -> GpQuaternion {
        let mut q = GpQuaternion::identity();
        q.set_matrix(&self.matrix);
        q
    }

    /// `gp_Trsf::Value(row, col)` — 1-based, column 4 is the translation.
    pub fn value(&self, row: i32, col: i32) -> f64 {
        if col < 4 {
            self.scale * self.matrix.value(row as usize, col as usize)
        } else {
            match row {
                1 => self.loc.x,
                2 => self.loc.y,
                _ => self.loc.z,
            }
        }
    }

    /// `gp_Trsf::GetMat4(NCollection_Mat4<T>&)` (`gp_Trsf.hxx:328-352`):
    /// row-major 4x4 with rows 1-3 taken from `Value(row, col)` and last row
    /// `[0, 0, 0, 1]`. `NCollection_Mat4` has no port equivalent, so the 4x4 is
    /// mirrored as a plain `[[f64; 4]; 4]` with `SetValue(row, col)` semantics.
    pub fn get_mat4(&self) -> [[f64; 4]; 4] {
        if self.shape == TrsfForm::Identity {
            return [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ];
        }
        let mut m = [[0.0f64; 4]; 4];
        for r in 1..=3usize {
            for c in 1..=4usize {
                m[r - 1][c - 1] = self.value(r as i32, c as i32);
            }
        }
        m[3][3] = 1.0;
        m
    }

    /// UNPORTED: `gp_Trsf::DumpJson` (`gp_Trsf.cxx:943-960`) and
    /// `gp_Trsf::InitFromJson` (`gp_Trsf.cxx:964-1009`) need the OCCT
    /// `Standard_Dump` / `Standard_SStream` / `TCollection_AsciiString`
    /// serialization framework and the `OCCT_DUMP_*` / `OCCT_INIT_*` macros,
    /// none of which exist in this port. No equivalent branch, so not added.
    pub fn set_translation_vec(&mut self, v: &GpVec) {
        self.set_identity();
        self.shape = TrsfForm::Translation;
        self.loc = v.xyz().clone();
    }

    pub fn set_translation_pnts(&mut self, p1: &GpPnt, p2: &GpPnt) {
        let v = GpVec::from_pnts(p1, p2);
        self.set_translation_vec(&v);
    }

    pub fn set_mirror_pnt(&mut self, p: &GpPnt) {
        self.set_identity();
        self.scale = -1.0;
        self.shape = TrsfForm::PntMirror;
        self.loc = p.coord.multiply(2.0);
    }

    /// `gp_Trsf::SetMirror(const gp_Ax1&)` (`gp_Trsf.cxx:57-71`). The stored
    /// matrix is the reflection `2*v*v^T - I` (`scale = 1`); the location is
    /// obtained by applying the pre-negation matrix `I - 2*v*v^T` to the axis
    /// origin and adding the origin back (`loc = 2*o - 2*(v.o)*v`).
    pub fn set_mirror_ax1(&mut self, ax1: &GpAx1) {
        self.shape = TrsfForm::Ax1Mirror;
        self.scale = 1.0;
        let o = ax1.location().coord;
        let d = ax1.direction();
        let vx = d.x();
        let vy = d.y();
        let vz = d.z();

        // matrix = I - 2*v*v^T (SetDot * -2, then diagonal +1):
        let mut m = GpMat::new(
            1.0 - 2.0 * vx * vx, -2.0 * vx * vy,      -2.0 * vx * vz,
            -2.0 * vy * vx,      1.0 - 2.0 * vy * vy, -2.0 * vy * vz,
            -2.0 * vz * vx,      -2.0 * vz * vy,      1.0 - 2.0 * vz * vz,
        );
        // loc = matrix * loc + A1.Location():
        let loc = m.multiplied(&o).add(&o);
        // matrix.Multiply(-1)
        m = m.multiply_scalar(-1.0);

        self.matrix = m;
        self.loc = loc;
    }

    /// `gp_Trsf::SetMirror(const gp_Ax2&)` (`gp_Trsf.cxx:73-86`). The stored
    /// matrix is `2*n*n^T - I` with `scale = -1` (the plane reflection is the
    /// product `scale * matrix`); the location is `2*(n.o)*n`.
    pub fn set_mirror_ax2(&mut self, ax2: &GpAx2) {
        self.shape = TrsfForm::Ax2Mirror;
        self.scale = -1.0;
        let o = ax2.location().coord;
        let n = ax2.direction();
        let nx = n.x();
        let ny = n.y();
        let nz = n.z();

        // matrix = 2*n*n^T - I (SetDot * 2, then diagonal -1):
        let m = GpMat::new(
            2.0 * nx * nx - 1.0, 2.0 * nx * ny,       2.0 * nx * nz,
            2.0 * ny * nx,       2.0 * ny * ny - 1.0, 2.0 * ny * nz,
            2.0 * nz * nx,       2.0 * nz * ny,       2.0 * nz * nz - 1.0,
        );
        // loc = matrix * loc + A2.Location():
        let loc = m.multiplied(&o).add(&o);

        self.matrix = m;
        self.loc = loc;
    }

    pub fn set_mirror_ax3(&mut self, ax3: &GpAx3) {
        let ax2 = ax3.ax2();
        self.set_mirror_ax2(&ax2);
    }

    /// `gp_Trsf::SetRotation(const gp_Ax1&, const double)`
    /// (`gp_Trsf.cxx:90-99`): direct construction, `shape = Rotation` and
    /// `loc = -R*o + o`. `_angle` is kept named `angle` below.
    pub fn set_rotation_ax1(&mut self, ax1: &GpAx1, angle: f64) -> Result<(), &'static str> {
        self.shape = TrsfForm::Rotation;
        self.scale = 1.0;
        let o = ax1.location().coord;
        self.loc = o;
        self.matrix.set_rotation(ax1.direction().xyz(), angle)?;
        self.loc.reverse();
        self.loc.multiply_mat(&self.matrix);
        self.loc = self.loc.add(&o);
        Ok(())
    }

    /// `gp_Trsf::SetTransformation(const gp_Quaternion& R, const gp_Vec& T)`
    /// (`gp_Trsf.cxx:207-213`).
    pub fn set_transformation_quat_vec(&mut self, r: &GpQuaternion, t: &GpVec) {
        self.shape = TrsfForm::CompoundTrsf;
        self.scale = 1.0;
        self.loc = *t.xyz();
        self.matrix = r.get_matrix();
    }

    /// `gp_Trsf::SetRotation(const gp_Quaternion&)` (`gp_Trsf.cxx:103-109`).
    pub fn set_rotation_quat(&mut self, q: &GpQuaternion) {
        self.shape = TrsfForm::Rotation;
        self.scale = 1.0;
        self.loc = GpXyz::zero();
        self.matrix = q.get_matrix();
    }

    /// `gp_Trsf::SetRotationPart(const gp_Quaternion&)` (`gp_Trsf.cxx:111-156`):
    /// replaces the rotation part and reclassifies `shape` accordingly.
    pub fn set_rotation_part(&mut self, q: &GpQuaternion) {
        let has_rotation = !q.is_identity();
        if has_rotation {
            self.matrix = q.get_matrix();
        } else {
            self.matrix.set_identity();
        }

        match self.shape {
            TrsfForm::Identity => {
                if has_rotation {
                    self.shape = TrsfForm::Rotation;
                }
            }
            TrsfForm::Rotation => {
                if !has_rotation {
                    self.shape = TrsfForm::Identity;
                }
            }
            TrsfForm::Translation
            | TrsfForm::PntMirror
            | TrsfForm::Ax1Mirror
            | TrsfForm::Ax2Mirror
            | TrsfForm::Scale
            | TrsfForm::CompoundTrsf
            | TrsfForm::Other => {
                if has_rotation {
                    self.shape = TrsfForm::CompoundTrsf;
                }
            }
        }
    }

    /// `gp_Trsf::SetTranslationPart(const gp_Vec&)` (`gp_Trsf.cxx:243-280`).
    pub fn set_translation_part(&mut self, v: &GpVec) {
        self.loc = v.xyz().clone();
        let loc_null = self.loc.square_modulus() < RESOLUTION;

        match self.shape {
            TrsfForm::Identity => {
                if !loc_null {
                    self.shape = TrsfForm::Translation;
                }
            }
            TrsfForm::Translation => {
                if loc_null {
                    self.shape = TrsfForm::Identity;
                }
            }
            TrsfForm::Rotation
            | TrsfForm::PntMirror
            | TrsfForm::Ax1Mirror
            | TrsfForm::Ax2Mirror
            | TrsfForm::Scale
            | TrsfForm::CompoundTrsf
            | TrsfForm::Other => {
                if !loc_null {
                    self.shape = TrsfForm::CompoundTrsf;
                }
            }
        }
    }

    /// `gp_Trsf::SetScaleFactor(const double)` (`gp_Trsf.cxx:284-337`).
    pub fn set_scale_factor(&mut self, s: f64) -> Result<(), &'static str> {
        if s.abs() <= RESOLUTION {
            return Err("gp_Trsf::SetScaleFactor: scale factor is too small");
        }
        self.scale = s;
        let unit = (self.scale - 1.0).abs() <= RESOLUTION;
        let munit = (self.scale + 1.0).abs() <= RESOLUTION;

        match self.shape {
            TrsfForm::Identity | TrsfForm::Translation => {
                if !unit {
                    self.shape = TrsfForm::Scale;
                }
                if munit {
                    self.shape = TrsfForm::PntMirror;
                }
            }
            TrsfForm::Rotation => {
                if !unit {
                    self.shape = TrsfForm::CompoundTrsf;
                }
            }
            TrsfForm::PntMirror | TrsfForm::Ax1Mirror | TrsfForm::Ax2Mirror => {
                if !munit {
                    self.shape = TrsfForm::Scale;
                }
                if unit {
                    self.shape = TrsfForm::Identity;
                }
            }
            TrsfForm::Scale => {
                if unit {
                    self.shape = TrsfForm::Identity;
                }
                if munit {
                    self.shape = TrsfForm::PntMirror;
                }
            }
            TrsfForm::CompoundTrsf | TrsfForm::Other => {}
        }
        Ok(())
    }

    /// `gp_Trsf::SetValues(...)` (`gp_Trsf.cxx:346-385`): build a compound
    /// transform from the 12 coefficients (3 columns with translation last).
    /// `scale` is the cube root of the determinant, then `Orthogonalize`.
    #[allow(clippy::too_many_arguments)]
    pub fn set_values(
        &mut self,
        a11: f64, a12: f64, a13: f64, a14: f64,
        a21: f64, a22: f64, a23: f64, a24: f64,
        a31: f64, a32: f64, a33: f64, a34: f64,
    ) -> Result<(), &'static str> {
        let matrix = GpMat::new(a11, a12, a13, a21, a22, a23, a31, a32, a33);
        let det = matrix.determinant();
        if det.abs() < RESOLUTION {
            return Err("gp_Trsf::SetValues, null determinant");
        }
        let s = if det > 0.0 {
            det.powf(1.0 / 3.0)
        } else {
            -(-det).powf(1.0 / 3.0)
        };

        self.scale = s;
        self.shape = TrsfForm::CompoundTrsf;
        self.matrix = matrix.divide(s)?;
        self.orthogonalize();
        self.loc = GpXyz::new(a14, a24, a34);
        Ok(())
    }

    /// `gp_Trsf::Orthogonalize()` (`gp_Trsf.cxx:862-938`): Gram-Schmidt on the
    /// columns then on the rows of the h-vectorial matrix.
    pub fn orthogonalize(&mut self) {
        let mut tm = self.matrix;
        let mut v1 = tm.column(1);
        let mut v2 = tm.column(2);
        let mut v3 = tm.column(3);

        v1.normalize();
        v2 = v2.subtract(&v1.multiply(v2.dot(&v1)));
        v2.normalize();
        v3 = v3
            .subtract(&v1.multiply(v3.dot(&v1)))
            .subtract(&v2.multiply(v3.dot(&v2)));
        v3.normalize();
        tm.set_cols(&v1, &v2, &v3);

        let mut r1 = tm.row(1);
        let mut r2 = tm.row(2);
        let mut r3 = tm.row(3);

        r1.normalize();
        r2 = r2.subtract(&r1.multiply(r2.dot(&r1)));
        r2.normalize();
        r3 = r3
            .subtract(&r1.multiply(r3.dot(&r1)))
            .subtract(&r2.multiply(r3.dot(&r2)));
        r3.normalize();
        tm.set_rows(&r1, &r2, &r3);

        self.matrix = tm;
    }

    /// `gp_Trsf::SetScale(const gp_Pnt& theP, const double theS)`
    /// (`gp_Trsf.cxx:159-168`): homothety of factor `theS` centered at `theP`.
    /// The matrix stays identity and the factor is carried by `scale`, so the
    /// stored form is always `Scale` (there is no identity shortcut in OCCT).
    pub fn set_scale(&mut self, p: &GpPnt, s: f64) -> Result<(), &'static str> {
        if s.abs() <= RESOLUTION {
            return Err("gp_Trsf::SetScaleFactor: scale factor is too small");
        }

        self.shape = TrsfForm::Scale;
        self.scale = s;
        self.matrix = GpMat::identity();
        // loc = P * (1 - S)
        self.loc = p.coord.multiply(1.0 - s);
        Ok(())
    }

    pub fn transforms_xyz(&self, xyz: &mut GpXyz) {
        let mut tmp = self.matrix.multiplied(&*xyz);
        if self.scale != 1.0 {
            tmp = tmp.multiply(self.scale);
        }
        *xyz = tmp.add(&self.loc);
    }

    /// Transform a direction (rotation + scale only; no translation).
    /// Source: `gp_Trsf::Transforms` on a direction is the linear part.
    pub fn transforms_xyz_dir(&self, xyz: &mut GpXyz) {
        let mut tmp = self.matrix.multiplied(&*xyz);
        if self.scale != 1.0 {
            tmp = tmp.multiply(self.scale);
        }
        *xyz = tmp;
    }

    pub fn multiply(&mut self, other: &GpTrsf) {
        if self.form() == TrsfForm::Identity {
            *self = other.clone();
            return;
        }
        if other.form() == TrsfForm::Identity {
            return;
        }

        let new_scale = self.scale * other.scale;

        // new_loc = self.scale * self.matrix * other.loc + self.loc
        let tmp = other.loc.multiply(self.scale);
        let new_loc = self.matrix.multiplied(&tmp).add(&self.loc);

        let new_matrix = self.matrix.multiply(&other.matrix);

        let new_shape = if self.shape == TrsfForm::Identity {
            other.shape
        } else if other.shape == TrsfForm::Identity {
            self.shape
        } else {
            TrsfForm::CompoundTrsf
        };

        self.scale = new_scale;
        self.matrix = new_matrix;
        self.loc = new_loc;
        self.shape = new_shape;
    }

    pub fn multiplied(&self, other: &GpTrsf) -> GpTrsf {
        let mut result = self.clone();
        result.multiply(other);
        result
    }

    /// `gp_Trsf::PreMultiply(const gp_Trsf&)` (`gp_Trsf.cxx:712-836`).
    pub fn pre_multiply(&mut self, t: &GpTrsf) {
        use TrsfForm::*;
        if t.shape == Identity {
            return;
        } else if self.shape == Identity {
            self.shape = t.shape;
            self.scale = t.scale;
            self.loc = t.loc;
            self.matrix = t.matrix;
        } else if self.shape == Rotation && t.shape == Rotation {
            self.loc = self.loc.multiplied_mat(&t.matrix).add(&t.loc);
            self.matrix = self.matrix.pre_multiply(&t.matrix);
        } else if self.shape == Translation && t.shape == Translation {
            self.loc = self.loc.add(&t.loc);
        } else if self.shape == Scale && t.shape == Scale {
            self.loc = self.loc.multiply(t.scale).add(&t.loc);
            self.scale = self.scale * t.scale;
        } else if self.shape == PntMirror && t.shape == PntMirror {
            self.scale = 1.0;
            self.shape = Translation;
            self.loc.reverse();
            self.loc = self.loc.add(&t.loc);
        } else if self.shape == Ax1Mirror && t.shape == Ax1Mirror {
            self.shape = Rotation;
            self.loc = self.loc.multiplied_mat(&t.matrix).add(&t.loc);
            self.matrix = self.matrix.pre_multiply(&t.matrix);
        } else if matches!(self.shape, CompoundTrsf | Rotation | Ax1Mirror | Ax2Mirror)
            && t.shape == Translation
        {
            self.loc = self.loc.add(&t.loc);
        } else if matches!(self.shape, Scale | PntMirror) && t.shape == Translation {
            self.loc = self.loc.add(&t.loc);
        } else if self.shape == Translation
            && matches!(t.shape, CompoundTrsf | Rotation | Ax1Mirror | Ax2Mirror)
        {
            self.shape = CompoundTrsf;
            self.matrix = t.matrix;
            if t.scale == 1.0 {
                self.loc = self.loc.multiplied_mat(&t.matrix);
            } else {
                self.scale = t.scale;
                self.loc = self.loc.multiplied_mat(&self.matrix).multiply(self.scale);
            }
            self.loc = self.loc.add(&t.loc);
        } else if matches!(t.shape, Scale | PntMirror) && self.shape == Translation {
            self.loc = self.loc.multiply(t.scale).add(&t.loc);
            self.scale = t.scale;
            self.shape = t.shape;
        } else if matches!(self.shape, PntMirror | Scale) && matches!(t.shape, PntMirror | Scale) {
            self.shape = CompoundTrsf;
            self.loc = self.loc.multiply(t.scale).add(&t.loc);
            self.scale = self.scale * t.scale;
        } else if matches!(self.shape, CompoundTrsf | Rotation | Ax1Mirror | Ax2Mirror)
            && matches!(t.shape, Scale | PntMirror)
        {
            self.shape = CompoundTrsf;
            self.loc = self.loc.multiply(t.scale).add(&t.loc);
            self.scale = self.scale * t.scale;
        } else if matches!(t.shape, CompoundTrsf | Rotation | Ax1Mirror | Ax2Mirror)
            && matches!(self.shape, Scale | PntMirror)
        {
            self.shape = CompoundTrsf;
            self.matrix = t.matrix;
            if t.scale == 1.0 {
                self.loc = self.loc.multiplied_mat(&t.matrix);
            } else {
                self.loc = self.loc.multiplied_mat(&self.matrix).multiply(t.scale);
                self.scale = t.scale * self.scale;
            }
            self.loc = self.loc.add(&t.loc);
        } else {
            self.shape = CompoundTrsf;
            self.loc = self.loc.multiplied_mat(&t.matrix);
            if t.scale != 1.0 {
                self.loc = self.loc.multiply(t.scale);
                self.scale = self.scale * t.scale;
            }
            self.loc = self.loc.add(&t.loc);
            self.matrix = self.matrix.pre_multiply(&t.matrix);
        }
    }

    /// `gp_Trsf::Power(const int)` (`gp_Trsf.cxx:564-708`). `IsOdd(n)` is
    /// `n % 2 == 1` (`Standard_Integer.hxx:39-42`).
    pub fn power(&mut self, n: i32) -> Result<(), &'static str> {
        use TrsfForm::*;
        let is_odd = |v: i32| v % 2 == 1;
        if self.shape == Identity {
            return Ok(());
        }
        if n == 0 {
            self.scale = 1.0;
            self.shape = Identity;
            self.matrix.set_identity();
            self.loc = GpXyz::zero();
        } else if n == 1 {
            // already this^1
        } else if n == -1 {
            self.invert()?;
        } else {
            if n < 0 {
                self.invert()?;
            }
            if self.shape == Translation {
                let mut npower = n.abs() - 1;
                let mut temploc = self.loc;
                loop {
                    if is_odd(npower) {
                        self.loc = self.loc.add(&temploc);
                    }
                    if npower == 1 {
                        break;
                    }
                    temploc = temploc.add(&temploc);
                    npower /= 2;
                }
            } else if self.shape == Scale {
                let mut npower = n.abs() - 1;
                let mut temploc = self.loc;
                let mut tempscale = self.scale;
                loop {
                    if is_odd(npower) {
                        self.loc = self.loc.add(&temploc.multiplied(self.scale));
                        self.scale = self.scale * tempscale;
                    }
                    if npower == 1 {
                        break;
                    }
                    temploc = temploc.add(&temploc.multiplied(tempscale));
                    tempscale = tempscale * tempscale;
                    npower /= 2;
                }
            } else if self.shape == Rotation {
                let mut npower = n.abs() - 1;
                let mut tempmatrix = self.matrix;
                if self.loc.x == 0.0 && self.loc.y == 0.0 && self.loc.z == 0.0 {
                    loop {
                        if is_odd(npower) {
                            self.matrix = self.matrix.multiply(&tempmatrix);
                        }
                        if npower == 1 {
                            break;
                        }
                        tempmatrix = tempmatrix.multiply(&tempmatrix);
                        npower /= 2;
                    }
                } else {
                    let mut temploc = self.loc;
                    loop {
                        if is_odd(npower) {
                            self.loc = self.loc.add(&temploc.multiplied_mat(&self.matrix));
                            self.matrix = self.matrix.multiply(&tempmatrix);
                        }
                        if npower == 1 {
                            break;
                        }
                        temploc = temploc.add(&temploc.multiplied_mat(&tempmatrix));
                        tempmatrix = tempmatrix.multiply(&tempmatrix);
                        npower /= 2;
                    }
                }
            } else if matches!(self.shape, PntMirror | Ax1Mirror | Ax2Mirror) {
                if n % 2 == 0 {
                    self.shape = Identity;
                    self.scale = 1.0;
                    self.matrix.set_identity();
                    self.loc = GpXyz::zero();
                }
            } else {
                self.shape = CompoundTrsf;
                let mut npower = n.abs() - 1;
                let mut temploc = self.loc;
                let mut tempscale = self.scale;
                let mut tempmatrix = self.matrix;
                loop {
                    if is_odd(npower) {
                        self.loc = self
                            .loc
                            .add(&temploc.multiplied_mat(&self.matrix).multiply(self.scale));
                        self.scale = self.scale * tempscale;
                        self.matrix = self.matrix.multiply(&tempmatrix);
                    }
                    if npower == 1 {
                        break;
                    }
                    tempscale = tempscale * tempscale;
                    temploc = temploc.add(&temploc.multiplied_mat(&tempmatrix).multiply(tempscale));
                    tempmatrix = tempmatrix.multiply(&tempmatrix);
                    npower /= 2;
                }
            }
        }
        Ok(())
    }

    pub fn invert(&mut self) -> Result<(), &'static str> {
        if self.scale.abs() < RESOLUTION {
            return Err("transform is singular (zero scale)");
        }

        let inv_scale = 1.0 / self.scale;

        // Full matrix inverse (handles compound transforms with non-orthogonal matrices)
        let inv_matrix = self.matrix.inverted()?;

        // inv_loc = -inv_scale * inv_matrix * loc
        let inv_loc = self.loc.multiply(-1.0);
        let inv_loc = inv_matrix.multiplied(&inv_loc);
        let inv_loc = inv_loc.multiply(inv_scale);

        self.scale = inv_scale;
        self.matrix = inv_matrix;
        self.loc = inv_loc;
        Ok(())
    }

    pub fn inverted(&self) -> Result<GpTrsf, &'static str> {
        let mut result = self.clone();
        result.invert()?;
        Ok(result)
    }
}

// ponytail: only &GpTrsf * &GpTrsf, single Deref-free impl
impl Mul<&GpTrsf> for &GpTrsf {
    type Output = GpTrsf;

    fn mul(self, rhs: &GpTrsf) -> GpTrsf {
        self.multiplied(rhs)
    }
}
