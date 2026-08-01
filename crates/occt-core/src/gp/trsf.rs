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

    pub fn vectorial_part(&self) -> GpMat {
        self.matrix.clone()
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

    pub fn set_mirror_ax1(&mut self, ax1: &GpAx1) {
        // Line reflection through axis: M = 2 * v * v^T - I
        let d = ax1.direction();
        let vx = d.x();
        let vy = d.y();
        let vz = d.z();

        let m = GpMat::new(
            2.0 * vx * vx - 1.0, 2.0 * vx * vy,       2.0 * vx * vz,
            2.0 * vy * vx,       2.0 * vy * vy - 1.0, 2.0 * vy * vz,
            2.0 * vz * vx,       2.0 * vz * vy,       2.0 * vz * vz - 1.0,
        );

        let o = ax1.location().coord;
        let loc_neg = o.multiply(-1.0);
        let loc = m.multiplied(&o).add(&loc_neg);

        self.scale = 1.0;
        self.shape = TrsfForm::Ax1Mirror;
        self.matrix = m;
        self.loc = loc;
    }

    pub fn set_mirror_ax2(&mut self, ax2: &GpAx2) {
        // Plane reflection through ax2 plane (XY plane of the coordinate system)
        let n = ax2.direction(); // normal of the plane
        let nx = n.x();
        let ny = n.y();
        let nz = n.z();

        let m = GpMat::new(
            1.0 - 2.0 * nx * nx, -2.0 * nx * ny, -2.0 * nx * nz,
            -2.0 * ny * nx, 1.0 - 2.0 * ny * ny, -2.0 * ny * nz,
            -2.0 * nz * nx, -2.0 * nz * ny, 1.0 - 2.0 * nz * nz,
        );

        let o = ax2.location().coord;
        let o2 = o.multiply(2.0);
        let loc = m.multiplied(&o).subtracted(&o2);
        let loc = loc.multiply(-1.0);

        self.scale = 1.0;
        self.shape = TrsfForm::Ax2Mirror;
        self.matrix = m;
        self.loc = loc;
    }

    pub fn set_mirror_ax3(&mut self, ax3: &GpAx3) {
        let ax2 = ax3.ax2();
        self.set_mirror_ax2(&ax2);
    }

    pub fn set_rotation_ax1(&mut self, ax1: &GpAx1, angle: f64) -> Result<(), &'static str> {
        // Rodrigues rotation: translate to origin, rotate, translate back
        let o = ax1.location().coord;
        let d = ax1.direction();

        // Build rotation-to-origin translation
        let mut tr1 = GpTrsf::identity();
        let loc_neg = o.multiply(-1.0);
        tr1.loc = loc_neg;
        tr1.shape = TrsfForm::Translation;

        // Build rotation around axis through origin
        let mut rot = GpTrsf::identity();
        rot.matrix.set_rotation(d.xyz(), angle);
        rot.shape = TrsfForm::Rotation;

        // Build translation back
        let mut tr2 = GpTrsf::identity();
        tr2.loc = o.clone();
        tr2.shape = TrsfForm::Translation;

        *self = tr2.multiplied(&rot).multiplied(&tr1);
        Ok(())
    }

    pub fn set_rotation_quat(&mut self, q: &GpQuaternion) {
        self.set_identity();
        self.shape = TrsfForm::Rotation;
        self.matrix = q.get_matrix();
    }

    pub fn set_scale(&mut self, p: &GpPnt, s: f64) -> Result<(), &'static str> {
        if (s - 1.0).abs() < RESOLUTION {
            *self = Self::identity();
            return Ok(());
        }

        if s.abs() < RESOLUTION {
            return Err("scale factor is too small");
        }

        let o = p.coord;
        // translation: p - s * p = (1 - s) * p
        let loc = o.multiply(1.0 - s);

        self.shape = TrsfForm::Scale;
        self.scale = s;
        self.matrix = GpMat::identity();
        self.loc = loc;
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

    pub fn invert(&mut self) -> Result<(), &'static str> {
        if self.scale.abs() < RESOLUTION {
            return Err("transform is singular (zero scale)");
        }

        let inv_scale = 1.0 / self.scale;

        // Full matrix inverse (handles compound transforms with non-orthogonal matrices)
        let inv_matrix = self.matrix.invert();

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
