//! Ellipse in 3D. Source: `gp_Elips.hxx`
use crate::gp::{ax2::GpAx2, ax1::GpAx1, pnt::GpPnt, vec::GpVec, trsf::GpTrsf};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpElips { pub pos: GpAx2, pub major_radius: f64, pub minor_radius: f64 }

impl GpElips {
    pub fn new(pos: GpAx2, major: f64, minor: f64) -> Self { Self { pos, major_radius: major, minor_radius: minor } }
    #[inline] pub fn area(&self) -> f64 { std::f64::consts::PI * self.major_radius * self.minor_radius }
    #[inline] pub fn focal(&self) -> f64 { 2.0 * (self.major_radius * self.major_radius - self.minor_radius * self.minor_radius).sqrt() }
    pub fn eccentricity(&self) -> f64 { (self.major_radius*self.major_radius - self.minor_radius*self.minor_radius).sqrt() / self.major_radius }
    #[inline] pub fn location(&self) -> GpPnt { self.pos.location() }
    #[inline] pub fn position(&self) -> &GpAx2 { &self.pos }
    #[inline] pub fn axis(&self) -> &GpAx1 { self.pos.axis() }
    #[inline] pub fn major_radius(&self) -> f64 { self.major_radius }
    #[inline] pub fn minor_radius(&self) -> f64 { self.minor_radius }
    pub fn set_location(&mut self, p: GpPnt) { self.pos.set_location(p); }
    pub fn set_position(&mut self, a2: GpAx2) { self.pos = a2; }
    pub fn set_axis(&mut self, a1: GpAx1) { self.pos.set_direction(*a1.direction()); }
    pub fn set_major_radius(&mut self, r: f64) { self.major_radius = r; }
    pub fn set_minor_radius(&mut self, r: f64) { self.minor_radius = r; }
    pub fn x_axis(&self) -> GpAx1 { self.pos.x_axis() }
    pub fn y_axis(&self) -> GpAx1 { self.pos.y_axis() }
    pub fn focus1(&self) -> GpPnt { let f=self.focal()*0.5; GpPnt::from_xyz(&self.pos.x_direction().xyz().multiplied(f).added(&self.location().coord)) }
    pub fn focus2(&self) -> GpPnt { let f=self.focal()*0.5; GpPnt::from_xyz(&self.pos.x_direction().xyz().multiplied(-f).added(&self.location().coord)) }
    pub fn mirror_pnt(&mut self, p: &GpPnt) { let mut t=GpTrsf::identity(); t.set_mirror_pnt(p); self.transform(&t); }
    pub fn mirrored_pnt(&self, p: &GpPnt) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax1(&mut self, a1: &GpAx1) { let mut t=GpTrsf::identity(); t.set_mirror_ax1(a1); self.transform(&t); }
    pub fn mirrored_ax1(&self, a1: &GpAx1) -> Self { let mut r=*self; r.mirror_ax1(a1); r }
    pub fn mirror_ax2(&mut self, a2: &GpAx2) { let mut t=GpTrsf::identity(); t.set_mirror_ax2(a2); self.transform(&t); }
    pub fn mirrored_ax2(&self, a2: &GpAx2) -> Self { let mut r=*self; r.mirror_ax2(a2); r }
    pub fn rotate(&mut self, a1: &GpAx1, angle: f64) { self.pos.axis.loc.rotate(a1, angle); }
    pub fn rotated(&self, a1: &GpAx1, angle: f64) -> Self { let mut r=*self; r.rotate(a1, angle); r }
    pub fn scale(&mut self, p: &GpPnt, s: f64) { self.pos.axis.loc.scale(p, s); self.major_radius*=s.abs(); self.minor_radius*=s.abs(); }
    pub fn scaled(&self, p: &GpPnt, s: f64) -> Self { let mut r=*self; r.scale(p, s); r }
    pub fn transform(&mut self, t: &GpTrsf) { self.pos.axis.loc.transform(t); let s=t.scale_factor().abs(); self.major_radius*=s; self.minor_radius*=s; }
    pub fn transformed(&self, t: &GpTrsf) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self, v: &GpVec) { self.pos.axis.loc.translate_vec(v); }
    pub fn translated_vec(&self, v: &GpVec) -> Self { let mut r=*self; r.translate_vec(v); r }
    pub fn translate_pnts(&mut self, p1: &GpPnt, p2: &GpPnt) { self.pos.axis.loc.translate_pnts(p1, p2); }
    pub fn translated_pnts(&self, p1: &GpPnt, p2: &GpPnt) -> Self { let mut r=*self; r.translate_pnts(p1, p2); r }
}

impl Default for GpElips { fn default() -> Self { Self { pos: GpAx2::standard(), major_radius: 1.0, minor_radius: 1.0 } } }
