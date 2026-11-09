//! Parabola in 3D. Source: `gp_Parab.hxx`
use crate::gp::{ax2::GpAx2, ax1::GpAx1, pnt::GpPnt, vec::GpVec, trsf::GpTrsf};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpParab { pub pos: GpAx2, pub focal: f64 }

impl GpParab {
    pub fn new(pos: GpAx2, focal: f64) -> Self { Self { pos, focal } }
    pub fn focus(&self) -> GpPnt { GpPnt::from_xyz(&self.pos.x_direction().xyz().multiplied(self.focal).added(&self.pos.location().coord)) }
    #[inline] pub fn location(&self) -> GpPnt { self.pos.location() }
    #[inline] pub fn position(&self) -> &GpAx2 { &self.pos }
    #[inline] pub fn axis(&self) -> &GpAx1 { self.pos.axis() }
    pub fn set_location(&mut self, p: GpPnt) { self.pos.set_location(p); }
    pub fn set_position(&mut self, a2: GpAx2) { self.pos = a2; }
    pub fn set_axis(&mut self, a1: GpAx1) { self.pos.set_direction(*a1.direction()); }
    pub fn set_focal(&mut self, f: f64) { self.focal = f; }
    pub fn x_axis(&self) -> GpAx1 { self.pos.x_axis() }
    pub fn y_axis(&self) -> GpAx1 { self.pos.y_axis() }
    pub fn directrix(&self) -> GpAx1 { GpAx1::new(self.location(), *self.pos.y_direction()) }
    pub fn mirror_pnt(&mut self, p: &GpPnt) { self.pos.mirror_pnt(p); }
    pub fn mirrored_pnt(&self, p: &GpPnt) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax1(&mut self, a1: &GpAx1) { self.pos.mirror_ax1(a1); }
    pub fn mirrored_ax1(&self, a1: &GpAx1) -> Self { let mut r=*self; r.mirror_ax1(a1); r }
    pub fn mirror_ax2(&mut self, a2: &GpAx2) { self.pos.mirror_ax2(a2); }
    pub fn mirrored_ax2(&self, a2: &GpAx2) -> Self { let mut r=*self; r.mirror_ax2(a2); r }
    pub fn rotate(&mut self, a1: &GpAx1, angle: f64) { self.pos.rotate(a1, angle); }
    pub fn rotated(&self, a1: &GpAx1, angle: f64) -> Self { let mut r=*self; r.rotate(a1, angle); r }
    pub fn scale(&mut self, p: &GpPnt, s: f64) { self.focal*=s; if self.focal<0.0 { self.focal=-self.focal; } self.pos.scale(p, s); }
    pub fn scaled(&self, p: &GpPnt, s: f64) -> Self { let mut r=*self; r.scale(p, s); r }
    pub fn transform(&mut self, t: &GpTrsf) { self.focal*=t.scale_factor(); if self.focal<0.0 { self.focal=-self.focal; } self.pos.transform(t); }
    pub fn transformed(&self, t: &GpTrsf) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self, v: &GpVec) { self.pos.translate_vec(v); }
    pub fn translated_vec(&self, v: &GpVec) -> Self { let mut r=*self; r.translate_vec(v); r }
    pub fn translate_pnts(&mut self, p1: &GpPnt, p2: &GpPnt) { self.pos.axis.loc.translate_pnts(p1, p2); }
    pub fn translated_pnts(&self, p1: &GpPnt, p2: &GpPnt) -> Self { let mut r=*self; r.translate_pnts(p1, p2); r }
}

impl Default for GpParab { fn default() -> Self { Self { pos: GpAx2::standard(), focal: 1.0 } } }
