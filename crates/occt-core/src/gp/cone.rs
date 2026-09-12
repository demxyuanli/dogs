use crate::gp::ax1::GpAx1;
use crate::gp::ax2::GpAx2;
use crate::gp::ax3::GpAx3;
use crate::gp::dir::GpDir;
use crate::gp::pnt::GpPnt;
use crate::gp::trsf::GpTrsf;
use crate::gp::vec::GpVec;

#[derive(Debug, Clone)]
pub struct GpCone {
    pub pos: GpAx3,
    pub radius: f64,
    pub semi_angle: f64,
}

impl Default for GpCone {
    fn default() -> Self {
        GpCone {
            pos: GpAx3::default(),
            radius: f64::MAX,
            semi_angle: 0.01,
        }
    }
}

impl GpCone {
    pub fn new(pos: GpAx3, radius: f64, semi_angle: f64) -> Result<Self, &'static str> {
        if radius < 0.0 {
            return Err("radius must be non-negative");
        }
        if semi_angle <= 0.0 || semi_angle >= std::f64::consts::FRAC_PI_2 {
            return Err("semi_angle must be in (0, PI/2)");
        }
        Ok(GpCone { pos, radius, semi_angle })
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

    pub fn semi_angle(&self) -> f64 {
        self.semi_angle
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

    pub fn set_semi_angle(&mut self, angle: f64) -> Result<(), &'static str> {
        if angle <= 0.0 || angle >= std::f64::consts::FRAC_PI_2 {
            return Err("semi_angle must be in (0, PI/2)");
        }
        self.semi_angle = angle;
        Ok(())
    }

    /// `gp_Cone::Coefficients`.
    pub fn coefficients(&self) -> [f64; 10] {
        let mut t = GpTrsf::identity();
        t.set_transformation(&self.pos);
        let k_ang = self.semi_angle.tan();
        let t11 = t.value(1, 1);
        let t12 = t.value(1, 2);
        let t13 = t.value(1, 3);
        let t14 = t.value(1, 4);
        let t21 = t.value(2, 1);
        let t22 = t.value(2, 2);
        let t23 = t.value(2, 3);
        let t24 = t.value(2, 4);
        let t31 = t.value(3, 1) * k_ang;
        let t32 = t.value(3, 2) * k_ang;
        let t33 = t.value(3, 3) * k_ang;
        let t34 = t.value(3, 4) * k_ang;
        let r = self.radius;
        [
            t11 * t11 + t21 * t21 - t31 * t31,
            t12 * t12 + t22 * t22 - t32 * t32,
            t13 * t13 + t23 * t23 - t33 * t33,
            t11 * t12 + t21 * t22 - t31 * t32,
            t11 * t13 + t21 * t23 - t31 * t33,
            t12 * t13 + t22 * t23 - t32 * t33,
            t11 * t14 + t21 * t24 - t31 * (r + t34),
            t12 * t14 + t22 * t24 - t32 * (r + t34),
            t13 * t14 + t23 * t24 - t33 * (r + t34),
            t14 * t14 + t24 * t24 - r * r - t34 * t34 - 2.0 * r * t34,
        ]
    }

    pub fn apex(&self) -> GpPnt {
        let d = self.radius / self.semi_angle.tan();
        let direction = self.pos.direction();
        let dir = direction.xyz();
        let offset = dir.multiply(d);
        let apex_xyz = self.pos.location().coord.add(&offset);
        GpPnt::from_xyz(&apex_xyz)
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

    pub fn u_reverse(&mut self) {
        self.pos.x_reverse();
    }

    pub fn v_reverse(&mut self) {
        self.pos.y_reverse();
    }

    pub fn mirror_pnt(&mut self, p: &GpPnt) {
        let mut t = GpTrsf::identity();
        t.set_mirror_pnt(p);
        self.transform(&t);
    }

    pub fn mirrored_pnt(&self, p: &GpPnt) -> GpCone {
        let mut result = self.clone();
        result.mirror_pnt(p);
        result
    }

    pub fn mirror_ax1(&mut self, ax1: &GpAx1) {
        let mut t = GpTrsf::identity();
        t.set_mirror_ax1(ax1);
        self.transform(&t);
    }

    pub fn mirrored_ax1(&self, ax1: &GpAx1) -> GpCone {
        let mut result = self.clone();
        result.mirror_ax1(ax1);
        result
    }

    pub fn mirror_ax2(&mut self, ax2: &GpAx2) {
        let mut t = GpTrsf::identity();
        t.set_mirror_ax2(ax2);
        self.transform(&t);
    }

    pub fn mirrored_ax2(&self, ax2: &GpAx2) -> GpCone {
        let mut result = self.clone();
        result.mirror_ax2(ax2);
        result
    }

    pub fn rotate(&mut self, ax1: &GpAx1, angle: f64) {
        let mut t = GpTrsf::identity();
        let _ = t.set_rotation_ax1(ax1, angle);
        self.transform(&t);
    }

    pub fn rotated(&self, ax1: &GpAx1, angle: f64) -> GpCone {
        let mut result = self.clone();
        result.rotate(ax1, angle);
        result
    }

    pub fn scale(&mut self, p: &GpPnt, s: f64) {
        let mut t = GpTrsf::identity();
        let _ = t.set_scale(p, s);
        self.radius *= s.abs();
        self.transform(&t);
    }

    pub fn scaled(&self, p: &GpPnt, s: f64) -> GpCone {
        let mut result = self.clone();
        result.scale(p, s);
        result
    }

    pub fn transform(&mut self, t: &GpTrsf) {
        let mut loc = self.pos.location().coord;
        t.transforms_xyz(&mut loc);
        self.pos.set_location(GpPnt::from_xyz(&loc));

        let mut x_xyz = *self.pos.x_direction().xyz();
        t.transforms_xyz_dir(&mut x_xyz);
        let mut y_xyz = *self.pos.y_direction().xyz();
        t.transforms_xyz_dir(&mut y_xyz);

        if let (Ok(xd), Ok(yd)) = (GpDir::from_xyz(&x_xyz), GpDir::from_xyz(&y_xyz)) {
            self.pos.set_x_direction(&xd);
            self.pos.set_y_direction(&yd);
        }

        self.radius *= t.scale_factor();
    }

    pub fn transformed(&self, t: &GpTrsf) -> GpCone {
        let mut result = self.clone();
        result.transform(t);
        result
    }

    pub fn translate_vec(&mut self, v: &GpVec) {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(v);
        self.transform(&t);
    }

    pub fn translated_vec(&self, v: &GpVec) -> GpCone {
        let mut result = self.clone();
        result.translate_vec(v);
        result
    }

    pub fn translate_pnts(&mut self, p1: &GpPnt, p2: &GpPnt) {
        let v = GpVec::from_pnts(p1, p2);
        self.translate_vec(&v);
    }

    pub fn translated_pnts(&self, p1: &GpPnt, p2: &GpPnt) -> GpCone {
        let mut result = self.clone();
        result.translate_pnts(p1, p2);
        result
    }
}
