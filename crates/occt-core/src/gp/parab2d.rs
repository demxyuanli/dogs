//! 2D parabola. Source: `gp_Parab2d.hxx`
use crate::gp::{ax22d::GpAx22d, pnt2d::GpPnt2d, vec2d::GpVec2d, trsf2d::GpTrsf2d, ax2d::GpAx2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpParab2d { pub pos: GpAx22d, pub focal: f64 }

impl GpParab2d {
    pub fn new(pos: GpAx22d, focal: f64) -> Self { Self { pos, focal } }
    /// `gp_Parab2d::Axis()`.
    #[inline] pub const fn axis(&self) -> &GpAx22d { &self.pos }
    /// `gp_Parab2d::Parameter()`: twice the focal length.
    #[inline] pub const fn parameter(&self) -> f64 { 2.0 * self.focal }
    pub fn set_location(&mut self, p: GpPnt2d) { self.pos.set_location(p); }
    pub fn set_axis(&mut self, a: GpAx22d) { self.pos = a; }
    pub fn mirror_pnt(&mut self, p: &GpPnt2d) { let mut t=GpTrsf2d::identity(); t.set_mirror_pnt(p); self.transform(&t); }
    pub fn mirrored_pnt(&self, p: &GpPnt2d) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax2d(&mut self, a: &GpAx2d) { let mut t=GpTrsf2d::identity(); t.set_mirror_ax2d(a); self.transform(&t); }
    pub fn mirrored_ax2d(&self, a: &GpAx2d) -> Self { let mut r=*self; r.mirror_ax2d(a); r }
    pub fn rotate(&mut self, p: &GpPnt2d, angle: f64) { self.pos.point.rotate(p, angle); }
    pub fn rotated(&self, p: &GpPnt2d, angle: f64) -> Self { let mut r=*self; r.rotate(p, angle); r }
    pub fn scale(&mut self, p: &GpPnt2d, s: f64) { self.pos.point.scale(p, s); self.focal*=s.abs(); }
    pub fn scaled(&self, p: &GpPnt2d, s: f64) -> Self { let mut r=*self; r.scale(p, s); r }
    pub fn transform(&mut self, t: &GpTrsf2d) { t.transforms_xy(&mut self.pos.point.coord); self.focal*=t.scale_factor().abs(); }
    pub fn transformed(&self, t: &GpTrsf2d) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self, v: &GpVec2d) { self.pos.point.translate_vec(v); }
    pub fn translated_vec(&self, v: &GpVec2d) -> Self { let mut r=*self; r.translate_vec(v); r }
}

impl Default for GpParab2d { fn default() -> Self { Self { pos: GpAx22d::standard(), focal: 1.0 } } }
