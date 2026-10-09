//! 2D circle. Source: `gp_Circ2d.hxx`
//! Not in standard OCCT gp but needed by Geom2d.
use crate::gp::{ax22d::GpAx22d, pnt2d::GpPnt2d, vec2d::GpVec2d, ax2d::GpAx2d, trsf2d::GpTrsf2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpCirc2d { pub pos: GpAx22d, pub radius: f64 }

impl GpCirc2d {
    pub fn new(pos: GpAx22d, radius: f64) -> Self { Self { pos, radius } }
    pub fn area(&self) -> f64 { std::f64::consts::PI * self.radius * self.radius }
    pub fn length(&self) -> f64 { 2.0 * std::f64::consts::PI * self.radius }
    #[inline] pub fn location(&self) -> GpPnt2d { self.pos.point }
    #[inline] pub fn position(&self) -> &GpAx22d { &self.pos }
    #[inline] pub fn radius(&self) -> f64 { self.radius }
    pub fn set_location(&mut self, p: GpPnt2d) { self.pos.set_location(p); }
    pub fn set_radius(&mut self, r: f64) { self.radius = r; }
    pub fn set_axis(&mut self, a: GpAx22d) { self.pos = a; }
    pub fn x_axis(&self) -> GpAx2d { GpAx2d::new(self.pos.point, *self.pos.x_direction()) }
    pub fn y_axis(&self) -> GpAx2d { GpAx2d::new(self.pos.point, *self.pos.y_direction()) }
    /// `gp_Circ2d::IsDirect` (`gp_Circ2d.hxx:181-184`): true when the local
    /// frame is direct, i.e. `XDirection().Crossed(YDirection()) >= 0`.
    #[inline] pub fn is_direct(&self) -> bool { self.pos.x_direction().crossed(self.pos.y_direction()) >= 0.0 }
    /// `gp_Circ2d::Reversed()` (`gp_Circ2d.hxx:288-295`): reverse the Y
    /// direction of the local frame, which flips the implicit orientation of
    /// the circle. `SetYDirection` recomputes `XDirection` as `(vy.Y, -vy.X)`.
    pub fn reversed(&self) -> Self {
        let mut r = *self;
        let mut vy = r.pos.vydir;
        vy.reverse();
        r.pos.set_y_direction(vy);
        r
    }
    pub fn mirror_pnt(&mut self, p: &GpPnt2d) { self.pos.mirror_pnt(p); }
    pub fn mirrored_pnt(&self, p: &GpPnt2d) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax2d(&mut self, a: &GpAx2d) { self.pos.mirror_ax2d(a); }
    pub fn mirrored_ax2d(&self, a: &GpAx2d) -> Self { let mut r=*self; r.mirror_ax2d(a); r }
    pub fn rotate(&mut self, p: &GpPnt2d, angle: f64) { self.pos.rotate(p, angle); }
    pub fn rotated(&self, p: &GpPnt2d, angle: f64) -> Self { let mut r=*self; r.rotate(p, angle); r }
    pub fn scale(&mut self, p: &GpPnt2d, s: f64) { self.radius *= s; if self.radius < 0.0 { self.radius = -self.radius; } self.pos.scale(p, s); }
    pub fn scaled(&self, p: &GpPnt2d, s: f64) -> Self { let mut r=*self; r.scale(p, s); r }
    pub fn transform(&mut self, t: &GpTrsf2d) { self.radius *= t.scale_factor(); if self.radius < 0.0 { self.radius = -self.radius; } self.pos.transform(t); }
    pub fn transformed(&self, t: &GpTrsf2d) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self, v: &GpVec2d) { self.pos.point.translate_vec(v); }
    pub fn translated_vec(&self, v: &GpVec2d) -> Self { let mut r=*self; r.translate_vec(v); r }
}

impl Default for GpCirc2d { fn default() -> Self { Self { pos: GpAx22d::standard(), radius: 1.0 } } }
