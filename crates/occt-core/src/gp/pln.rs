use crate::gp::ax1::GpAx1;
use crate::gp::ax2::GpAx2;
use crate::gp::ax3::GpAx3;
use crate::gp::dir::GpDir;
use crate::gp::pnt::GpPnt;
use crate::gp::trsf::GpTrsf;
use crate::gp::vec::GpVec;

#[derive(Debug, Clone)]
pub struct GpPln {
    pub pos: GpAx3,
}

impl Default for GpPln {
    fn default() -> Self {
        GpPln {
            pos: GpAx3::default(),
        }
    }
}

impl GpPln {
    pub fn new(pos: GpAx3) -> Self {
        GpPln { pos }
    }

    pub fn location(&self) -> GpPnt {
        self.pos.location()
    }

    pub fn position(&self) -> GpAx3 {
        self.pos.clone()
    }

    pub fn set_location(&mut self, p: &GpPnt) {
        self.pos.set_location(*p);
    }

    pub fn set_position(&mut self, ax3: &GpAx3) {
        self.pos = ax3.clone();
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

    /// `gp_Pln::Coefficients` — Ax + By + Cz + D = 0.
    pub fn coefficients(&self) -> (f64, f64, f64, f64) {
        let d = self.pos.direction();
        let (a, b, c) = if self.pos.is_direct() {
            (d.x(), d.y(), d.z())
        } else {
            (-d.x(), -d.y(), -d.z())
        };
        let p = self.pos.location();
        let cte = -(a * p.x() + b * p.y() + c * p.z());
        (a, b, c, cte)
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

    pub fn mirrored_pnt(&self, p: &GpPnt) -> GpPln {
        let mut result = self.clone();
        result.mirror_pnt(p);
        result
    }

    pub fn mirror_ax1(&mut self, ax1: &GpAx1) {
        let mut t = GpTrsf::identity();
        t.set_mirror_ax1(ax1);
        self.transform(&t);
    }

    pub fn mirrored_ax1(&self, ax1: &GpAx1) -> GpPln {
        let mut result = self.clone();
        result.mirror_ax1(ax1);
        result
    }

    pub fn mirror_ax2(&mut self, ax2: &GpAx2) {
        let mut t = GpTrsf::identity();
        t.set_mirror_ax2(ax2);
        self.transform(&t);
    }

    pub fn mirrored_ax2(&self, ax2: &GpAx2) -> GpPln {
        let mut result = self.clone();
        result.mirror_ax2(ax2);
        result
    }

    pub fn rotate(&mut self, ax1: &GpAx1, angle: f64) {
        let mut t = GpTrsf::identity();
        let _ = t.set_rotation_ax1(ax1, angle);
        self.transform(&t);
    }

    pub fn rotated(&self, ax1: &GpAx1, angle: f64) -> GpPln {
        let mut result = self.clone();
        result.rotate(ax1, angle);
        result
    }

    pub fn scale(&mut self, p: &GpPnt, s: f64) {
        let mut t = GpTrsf::identity();
        let _ = t.set_scale(p, s);
        self.transform(&t);
    }

    pub fn scaled(&self, p: &GpPnt, s: f64) -> GpPln {
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
    }

    pub fn transformed(&self, t: &GpTrsf) -> GpPln {
        let mut result = self.clone();
        result.transform(t);
        result
    }

    pub fn translate_vec(&mut self, v: &GpVec) {
        let mut t = GpTrsf::identity();
        t.set_translation_vec(v);
        self.transform(&t);
    }

    pub fn translated_vec(&self, v: &GpVec) -> GpPln {
        let mut result = self.clone();
        result.translate_vec(v);
        result
    }

    pub fn translate_pnts(&mut self, p1: &GpPnt, p2: &GpPnt) {
        let v = GpVec::from_pnts(p1, p2);
        self.translate_vec(&v);
    }

    pub fn translated_pnts(&self, p1: &GpPnt, p2: &GpPnt) -> GpPln {
        let mut result = self.clone();
        result.translate_pnts(p1, p2);
        result
    }
}
