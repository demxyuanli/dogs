//! Torus. Source: `gp_Torus.hxx`
use crate::gp::{ax3::GpAx3, ax1::GpAx1, ax2::GpAx2, pnt::GpPnt, vec::GpVec, trsf::GpTrsf};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpTorus { pub pos: GpAx3, pub major_radius: f64, pub minor_radius: f64 }

impl GpTorus {
    pub fn new(pos: GpAx3, major: f64, minor: f64) -> Result<Self, &'static str> { if major<0.0||minor<0.0 { Err("negative") } else { Ok(Self{pos,major_radius:major,minor_radius:minor}) } }
    #[inline] pub fn location(&self) -> GpPnt { self.pos.location() }
    #[inline] pub const fn position(&self) -> &GpAx3 { &self.pos }
    #[inline] pub const fn major_radius(&self) -> f64 { self.major_radius }
    #[inline] pub const fn minor_radius(&self) -> f64 { self.minor_radius }
    pub fn set_location(&mut self, l: GpPnt) { self.pos.set_location(l); }
    pub fn set_position(&mut self, p: GpAx3) { self.pos=p; }
    pub fn set_major_radius(&mut self, r: f64) -> Result<(),&'static str> { if r<0.0 { Err("negative") } else { self.major_radius=r; Ok(()) } }
    pub fn set_minor_radius(&mut self, r: f64) -> Result<(),&'static str> { if r<0.0 { Err("negative") } else { self.minor_radius=r; Ok(()) } }
    pub fn area(&self) -> f64 { 4.0*std::f64::consts::PI.powi(2)*self.major_radius*self.minor_radius }
    pub fn volume(&self) -> f64 { 2.0*std::f64::consts::PI.powi(2)*self.major_radius*self.minor_radius.powi(2) }
    #[inline] pub fn axis(&self) -> &GpAx1 { self.pos.axis() }
    #[inline] pub fn x_axis(&self) -> GpAx1 { GpAx1::new(self.pos.location(), *self.pos.x_direction()) }
    #[inline] pub fn y_axis(&self) -> GpAx1 { GpAx1::new(self.pos.location(), *self.pos.y_direction()) }
    #[inline] pub fn is_direct(&self) -> bool { self.pos.is_direct() }
    #[inline] pub fn u_reverse(&mut self) { self.pos.y_reverse(); }
    #[inline] pub fn v_reverse(&mut self) { self.pos.z_reverse(); }
    fn _mirror(&mut self, t: &GpTrsf) { self.transform(t); }
    pub fn mirror_pnt(&mut self, p: &GpPnt) { let mut t=GpTrsf::identity(); t.set_mirror_pnt(p); self._mirror(&t); }
    pub fn mirrored_pnt(&self, p: &GpPnt) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax1(&mut self, a1: &GpAx1) { let mut t=GpTrsf::identity(); t.set_mirror_ax1(a1); self._mirror(&t); }
    pub fn mirrored_ax1(&self, a1: &GpAx1) -> Self { let mut r=*self; r.mirror_ax1(a1); r }
    pub fn mirror_ax2(&mut self, a2: &GpAx2) { let mut t=GpTrsf::identity(); t.set_mirror_ax2(a2); self._mirror(&t); }
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

impl Default for GpTorus { fn default() -> Self { Self { pos: GpAx3::standard(), major_radius: f64::MAX, minor_radius: f64::MAX } } }
