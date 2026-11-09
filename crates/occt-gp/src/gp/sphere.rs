//! Sphere. Source: `gp_Sphere.hxx`
use crate::gp::{ax3::GpAx3, ax1::GpAx1, ax2::GpAx2, pnt::GpPnt, vec::GpVec, trsf::GpTrsf};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpSphere { pub pos: GpAx3, pub radius: f64 }

impl GpSphere {
    pub fn new(pos: GpAx3, r: f64) -> Result<Self, &'static str> { if r<0.0 { Err("negative") } else { Ok(Self{pos,radius:r}) } }
    #[inline] pub fn location(&self) -> GpPnt { self.pos.location() }
    #[inline] pub const fn position(&self) -> &GpAx3 { &self.pos }
    #[inline] pub const fn radius(&self) -> f64 { self.radius }
    pub fn set_location(&mut self, l: GpPnt) { self.pos.set_location(l); }
    pub fn set_position(&mut self, p: GpAx3) { self.pos=p; }
    pub fn set_radius(&mut self, r: f64) -> Result<(),&'static str> { if r<0.0 { Err("negative") } else { self.radius=r; Ok(()) } }
    pub fn area(&self) -> f64 { 4.0*std::f64::consts::PI*self.radius*self.radius }
    pub fn volume(&self) -> f64 { (4.0*std::f64::consts::PI*self.radius.powi(3))/3.0 }
    #[inline] pub fn x_axis(&self) -> GpAx1 { GpAx1::new(self.pos.location(), *self.pos.x_direction()) }
    #[inline] pub fn y_axis(&self) -> GpAx1 { GpAx1::new(self.pos.location(), *self.pos.y_direction()) }
    #[inline] pub fn is_direct(&self) -> bool { self.pos.is_direct() }

    /// `gp_Sphere::Coefficients` — local X^2+Y^2+Z^2-R^2 = 0 in world coords.
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
        let t31 = t.value(3, 1);
        let t32 = t.value(3, 2);
        let t33 = t.value(3, 3);
        let t34 = t.value(3, 4);
        [
            t11 * t11 + t21 * t21 + t31 * t31,
            t12 * t12 + t22 * t22 + t32 * t32,
            t13 * t13 + t23 * t23 + t33 * t33,
            t11 * t12 + t21 * t22 + t31 * t32,
            t11 * t13 + t21 * t23 + t31 * t33,
            t12 * t13 + t22 * t23 + t32 * t33,
            t11 * t14 + t21 * t24 + t31 * t34,
            t12 * t14 + t22 * t24 + t32 * t34,
            t13 * t14 + t23 * t24 + t33 * t34,
            t14 * t14 + t24 * t24 + t34 * t34 - self.radius * self.radius,
        ]
    }
    #[inline] pub fn u_reverse(&mut self) { self.pos.y_reverse(); }
    #[inline] pub fn v_reverse(&mut self) { self.pos.z_reverse(); }
    pub fn mirror_pnt(&mut self, p: &GpPnt) { self.pos.mirror_pnt(p); }
    pub fn mirrored_pnt(&self, p: &GpPnt) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax1(&mut self, a1: &GpAx1) { self.pos.mirror_ax1(a1); }
    pub fn mirrored_ax1(&self, a1: &GpAx1) -> Self { let mut r=*self; r.mirror_ax1(a1); r }
    pub fn mirror_ax2(&mut self, a2: &GpAx2) { self.pos.mirror_ax2(a2); }
    pub fn mirrored_ax2(&self, a2: &GpAx2) -> Self { let mut r=*self; r.mirror_ax2(a2); r }
    pub fn rotate(&mut self, a1: &GpAx1, angle: f64) { self.pos.rotate(a1, angle); }
    pub fn rotated(&self, a1: &GpAx1, angle: f64) -> Self { let mut r=*self; r.rotate(a1, angle); r }
    pub fn scale(&mut self, p: &GpPnt, s: f64) { self.pos.scale(p, s); self.radius*=s; if self.radius<0.0 { self.radius=-self.radius; } }
    pub fn scaled(&self, p: &GpPnt, s: f64) -> Self { let mut r=*self; r.scale(p, s); r }
    pub fn transform(&mut self, t: &GpTrsf) { self.pos.transform(t); self.radius*=t.scale_factor(); if self.radius<0.0 { self.radius=-self.radius; } }
    pub fn transformed(&self, t: &GpTrsf) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self, v: &GpVec) { self.pos.translate_vec(v); }
    pub fn translated_vec(&self, v: &GpVec) -> Self { let mut r=*self; r.translate_vec(v); r }
    pub fn translate_pnts(&mut self, p1: &GpPnt, p2: &GpPnt) { self.pos.axis.loc.translate_pnts(p1, p2); }
    pub fn translated_pnts(&self, p1: &GpPnt, p2: &GpPnt) -> Self { let mut r=*self; r.translate_pnts(p1, p2); r }
}

impl Default for GpSphere { fn default() -> Self { Self { pos: GpAx3::standard(), radius: f64::MAX } } }
