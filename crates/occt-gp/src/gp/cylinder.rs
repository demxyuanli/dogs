use crate::gp::ax1::GpAx1;
use crate::gp::ax2::GpAx2;
use crate::gp::ax3::GpAx3;
use crate::gp::pnt::GpPnt;
use crate::gp::trsf::GpTrsf;
use crate::gp::vec::GpVec;

#[derive(Debug, Clone)]
pub struct GpCylinder {
    pub pos: GpAx3,
    pub radius: f64,
}

impl Default for GpCylinder {
    fn default() -> Self {
        GpCylinder {
            pos: GpAx3::default(),
            radius: f64::MAX,
        }
    }
}

impl GpCylinder {
    pub fn new(pos: GpAx3, radius: f64) -> Result<Self, &'static str> {
        if radius < 0.0 {
            return Err("radius must be non-negative");
        }
        Ok(GpCylinder { pos, radius })
    }

    pub fn location(&self) -> GpPnt {
        self.pos.location()
    }

    pub fn position(&self) -> GpAx3 {
        self.pos.clone()
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    pub fn set_location(&mut self, p: &GpPnt) {
        self.pos.set_location(*p);
    }

    pub fn set_position(&mut self, ax3: &GpAx3) {
        self.pos = ax3.clone();
    }

    pub fn set_radius(&mut self, r: f64) -> Result<(), &'static str> {
        if r < 0.0 {
            return Err("radius must be non-negative");
        }
        self.radius = r;
        Ok(())
    }

    pub fn axis(&self) -> GpAx1 {
        *self.pos.axis()
    }

    pub fn x_axis(&self) -> GpAx1 {
        self.pos.x_axis()
    }

    pub fn y_axis(&self) -> GpAx1 {
        self.pos.y_axis()
    }

    pub fn is_direct(&self) -> bool {
        self.pos.is_direct()
    }

    /// `gp_Cylinder::Coefficients` — local X^2 + Y^2 - R^2 = 0 in world coords.
    pub fn coefficients(&self) -> [f64; 10] {
        let mut t = GpTrsf::identity();
        t.set_transformation(&self.pos);
        let t11 = t.value(1, 1);
        let t12 = t.value(1, 2);
        let t13 = t.value(1, 3);
        let t14 = t.value(1, 4);
        let t21 = t.value(2, 1);
        let t22 = t.value(2, 2);
        let t23 = t.value(2, 3);
        let t24 = t.value(2, 4);
        [
            t11 * t11 + t21 * t21,
            t12 * t12 + t22 * t22,
            t13 * t13 + t23 * t23,
            t11 * t12 + t21 * t22,
            t11 * t13 + t21 * t23,
            t12 * t13 + t22 * t23,
            t11 * t14 + t21 * t24,
            t12 * t14 + t22 * t24,
            t13 * t14 + t23 * t24,
            t14 * t14 + t24 * t24 - self.radius * self.radius,
        ]
    }

    /// `gp_Cylinder::UReverse` (`gp_Cylinder.hxx:86`): `pos.YReverse()`.
    pub fn u_reverse(&mut self) {
        self.pos.y_reverse();
    }

    /// `gp_Cylinder::VReverse` (`gp_Cylinder.hxx:90`): `pos.ZReverse()`.
    pub fn v_reverse(&mut self) {
        self.pos.z_reverse();
    }

    pub fn mirror_pnt(&mut self, p: &GpPnt) {
        self.pos.mirror_pnt(p);
    }

    pub fn mirrored_pnt(&self, p: &GpPnt) -> GpCylinder {
        let mut result = self.clone();
        result.mirror_pnt(p);
        result
    }

    pub fn mirror_ax1(&mut self, ax1: &GpAx1) {
        self.pos.mirror_ax1(ax1);
    }

    pub fn mirrored_ax1(&self, ax1: &GpAx1) -> GpCylinder {
        let mut result = self.clone();
        result.mirror_ax1(ax1);
        result
    }

    pub fn mirror_ax2(&mut self, ax2: &GpAx2) {
        self.pos.mirror_ax2(ax2);
    }

    pub fn mirrored_ax2(&self, ax2: &GpAx2) -> GpCylinder {
        let mut result = self.clone();
        result.mirror_ax2(ax2);
        result
    }

    pub fn rotate(&mut self, ax1: &GpAx1, angle: f64) {
        self.pos.rotate(ax1, angle);
    }

    pub fn rotated(&self, ax1: &GpAx1, angle: f64) -> GpCylinder {
        let mut result = self.clone();
        result.rotate(ax1, angle);
        result
    }

    pub fn scale(&mut self, p: &GpPnt, s: f64) {
        self.pos.scale(p, s);
        self.radius *= s;
        if self.radius < 0.0 {
            self.radius = -self.radius;
        }
    }

    pub fn scaled(&self, p: &GpPnt, s: f64) -> GpCylinder {
        let mut result = self.clone();
        result.scale(p, s);
        result
    }

    /// `gp_Cylinder::Transform` (`gp_Cylinder.hxx:229-237`).
    pub fn transform(&mut self, t: &GpTrsf) {
        self.pos.transform(t);
        self.radius *= t.scale_factor();
        if self.radius < 0.0 {
            self.radius = -self.radius;
        }
    }

    pub fn transformed(&self, t: &GpTrsf) -> GpCylinder {
        let mut result = self.clone();
        result.transform(t);
        result
    }

    pub fn translate(&mut self, v: &GpVec) {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(v);
        self.transform(&t);
    }

    pub fn translated(&self, v: &GpVec) -> GpCylinder {
        let mut result = self.clone();
        result.translate(v);
        result
    }
}
