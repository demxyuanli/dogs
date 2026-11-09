//! 2D ellipse. Source: `gp_Elips2d.hxx`
use crate::gp::{ax22d::GpAx22d, pnt2d::GpPnt2d, vec2d::GpVec2d, trsf2d::GpTrsf2d, ax2d::GpAx2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpElips2d { pub pos: GpAx22d, pub major_radius: f64, pub minor_radius: f64 }

impl GpElips2d {
    pub fn new(pos: GpAx22d, major: f64, minor: f64) -> Self { Self { pos, major_radius: major, minor_radius: minor } }
    pub fn area(&self) -> f64 { std::f64::consts::PI * self.major_radius * self.minor_radius }
    /// `gp_Elips2d::Axis()`.
    #[inline] pub const fn axis(&self) -> &GpAx22d { &self.pos }
    /// `gp_Elips2d::XAxis()`.
    #[inline] pub fn x_axis(&self) -> GpAx2d { self.pos.x_axis() }
    /// `gp_Elips2d::YAxis()`.
    #[inline] pub fn y_axis(&self) -> GpAx2d { self.pos.y_axis() }
    pub fn set_location(&mut self, p: GpPnt2d) { self.pos.set_location(p); }
    pub fn set_axis(&mut self, a: GpAx22d) { self.pos = a; }
    /// `gp_Elips2d::SetMajorRadius()`.
    #[inline] pub fn set_major_radius(&mut self, r: f64) { self.major_radius = r; }
    /// `gp_Elips2d::SetMinorRadius()`.
    #[inline] pub fn set_minor_radius(&mut self, r: f64) { self.minor_radius = r; }
    pub fn mirror_pnt(&mut self, p: &GpPnt2d) { self.pos.mirror_pnt(p); }
    pub fn mirrored_pnt(&self, p: &GpPnt2d) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax2d(&mut self, a: &GpAx2d) { self.pos.mirror_ax2d(a); }
    pub fn mirrored_ax2d(&self, a: &GpAx2d) -> Self { let mut r=*self; r.mirror_ax2d(a); r }
    pub fn rotate(&mut self, p: &GpPnt2d, angle: f64) { self.pos.rotate(p, angle); }
    pub fn rotated(&self, p: &GpPnt2d, angle: f64) -> Self { let mut r=*self; r.rotate(p, angle); r }
    pub fn scale(&mut self, p: &GpPnt2d, s: f64) { self.major_radius*=s; if self.major_radius<0.0 { self.major_radius=-self.major_radius; } self.minor_radius*=s; if self.minor_radius<0.0 { self.minor_radius=-self.minor_radius; } self.pos.scale(p, s); }
    pub fn scaled(&self, p: &GpPnt2d, s: f64) -> Self { let mut r=*self; r.scale(p, s); r }
    pub fn transform(&mut self, t: &GpTrsf2d) { let s=t.scale_factor().abs(); self.major_radius*=s; self.minor_radius*=s; self.pos.transform(t); }
    pub fn transformed(&self, t: &GpTrsf2d) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self, v: &GpVec2d) { self.pos.point.translate_vec(v); }
    pub fn translated_vec(&self, v: &GpVec2d) -> Self { let mut r=*self; r.translate_vec(v); r }
}

impl Default for GpElips2d { fn default() -> Self { Self { pos: GpAx22d::standard(), major_radius: 1.0, minor_radius: 1.0 } } }
