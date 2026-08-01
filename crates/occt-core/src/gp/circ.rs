use crate::gp::ax1::GpAx1;
use crate::gp::ax2::GpAx2;
use crate::gp::dir::GpDir;
use crate::gp::pnt::GpPnt;
use crate::gp::trsf::GpTrsf;
use crate::gp::vec::GpVec;

#[derive(Debug, Clone)]
pub struct GpCirc {
    pub pos: GpAx2,
    pub radius: f64,
}

impl Default for GpCirc {
    fn default() -> Self {
        GpCirc {
            pos: GpAx2::default(),
            radius: f64::MAX,
        }
    }
}

impl GpCirc {
    pub fn new(pos: GpAx2, radius: f64) -> Self {
        GpCirc { pos, radius }
    }

    pub fn location(&self) -> GpPnt {
        self.pos.location()
    }

    pub fn position(&self) -> GpAx2 {
        self.pos.clone()
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    pub fn set_location(&mut self, p: &GpPnt) {
        self.pos.set_location(*p);
    }

    pub fn set_position(&mut self, ax2: &GpAx2) {
        self.pos = ax2.clone();
    }

    pub fn set_radius(&mut self, r: f64) -> Result<(), &'static str> {
        if r < 0.0 {
            return Err("radius must be non-negative");
        }
        self.radius = r;
        Ok(())
    }

    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.radius * self.radius
    }

    pub fn length(&self) -> f64 {
        2.0 * std::f64::consts::PI * self.radius
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

    pub fn mirror_pnt(&mut self, p: &GpPnt) {
        let mut t = GpTrsf::identity();
        t.set_mirror_pnt(p);
        self.transform(&t);
    }

    pub fn mirrored_pnt(&self, p: &GpPnt) -> GpCirc {
        let mut result = self.clone();
        result.mirror_pnt(p);
        result
    }

    pub fn mirror_ax1(&mut self, ax1: &GpAx1) {
        let mut t = GpTrsf::identity();
        t.set_mirror_ax1(ax1);
        self.transform(&t);
    }

    pub fn mirrored_ax1(&self, ax1: &GpAx1) -> GpCirc {
        let mut result = self.clone();
        result.mirror_ax1(ax1);
        result
    }

    pub fn mirror_ax2(&mut self, ax2: &GpAx2) {
        let mut t = GpTrsf::identity();
        t.set_mirror_ax2(ax2);
        self.transform(&t);
    }

    pub fn mirrored_ax2(&self, ax2: &GpAx2) -> GpCirc {
        let mut result = self.clone();
        result.mirror_ax2(ax2);
        result
    }

    pub fn rotate(&mut self, ax1: &GpAx1, angle: f64) {
        let mut t = GpTrsf::identity();
        let _ = t.set_rotation_ax1(ax1, angle);
        self.transform(&t);
    }

    pub fn rotated(&self, ax1: &GpAx1, angle: f64) -> GpCirc {
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

    pub fn scaled(&self, p: &GpPnt, s: f64) -> GpCirc {
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
            self.pos.vxdir = xd;
            self.pos.vydir = yd;
        }

        self.radius *= t.scale_factor();
    }

    pub fn transformed(&self, t: &GpTrsf) -> GpCirc {
        let mut result = self.clone();
        result.transform(t);
        result
    }

    pub fn translate_vec(&mut self, v: &GpVec) {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(v);
        self.transform(&t);
    }

    pub fn translated_vec(&self, v: &GpVec) -> GpCirc {
        let mut result = self.clone();
        result.translate_vec(v);
        result
    }

    pub fn translate_pnts(&mut self, p1: &GpPnt, p2: &GpPnt) {
        let v = GpVec::from_pnts(p1, p2);
        self.translate_vec(&v);
    }

    pub fn translated_pnts(&self, p1: &GpPnt, p2: &GpPnt) -> GpCirc {
        let mut result = self.clone();
        result.translate_pnts(p1, p2);
        result
    }
}
