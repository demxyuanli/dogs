//! 2D line. Source: `gp_Lin2d.hxx`
use crate::gp::{ax2d::GpAx2d, pnt2d::GpPnt2d, dir2d::GpDir2d, vec2d::GpVec2d, trsf2d::GpTrsf2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpLin2d { pub pos: GpAx2d }

impl GpLin2d {
    pub fn new(a: GpAx2d) -> Self { Self { pos: a } }
    pub fn from_pnt_dir(p: GpPnt2d, d: GpDir2d) -> Self { Self { pos: GpAx2d::new(p, d) } }
    #[inline] pub fn position(&self) -> &GpAx2d { &self.pos }
    #[inline] pub fn location(&self) -> GpPnt2d { self.pos.loc }
    #[inline] pub fn direction(&self) -> &GpDir2d { &self.pos.vdir }
    pub fn set_position(&mut self, a: GpAx2d) { self.pos = a; }
    pub fn angle(&self, other: &Self) -> f64 { self.pos.vdir.angle(&other.pos.vdir) }
    pub fn distance(&self, p: &GpPnt2d) -> f64 {
        let dx = p.x() - self.pos.loc.x(); let dy = p.y() - self.pos.loc.y();
        let dot = dx * self.pos.vdir.x + dy * self.pos.vdir.y;
        (dx*dx + dy*dy - dot*dot).sqrt()
    }
    /// `gp_Lin2d::Coefficients(A, B, C)` (`gp_Lin2d.hxx:91-96`): the normalized
    /// coefficients of `A*X + B*Y + C = 0`.
    pub fn coefficients(&self) -> (f64, f64, f64) {
        let a = self.pos.vdir.y;
        let b = -self.pos.vdir.x;
        let c = -(a * self.pos.loc.x() + b * self.pos.loc.y());
        (a, b, c)
    }
    pub fn contains(&self, p: &GpPnt2d, tol: f64) -> bool { self.distance(p) <= tol }
    pub fn mirror_pnt(&mut self, p: &GpPnt2d) { self.pos.mirror_pnt(p); }
    pub fn mirrored_pnt(&self, p: &GpPnt2d) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax2d(&mut self, a: &GpAx2d) { self.pos.mirror_ax2d(a); }
    pub fn mirrored_ax2d(&self, a: &GpAx2d) -> Self { let mut r=*self; r.mirror_ax2d(a); r }
    pub fn rotate(&mut self, p: &GpPnt2d, angle: f64) { self.pos.rotate(p, angle); }
    pub fn rotated(&self, p: &GpPnt2d, angle: f64) -> Self { let mut r=*self; r.rotate(p, angle); r }
    pub fn scale(&mut self, p: &GpPnt2d, s: f64) { self.pos.scale(p, s); }
    pub fn scaled(&self, p: &GpPnt2d, s: f64) -> Self { let mut r=*self; r.scale(p, s); r }
    pub fn transform(&mut self, t: &GpTrsf2d) { self.pos.transform(t); }
    pub fn transformed(&self, t: &GpTrsf2d) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self, v: &GpVec2d) { self.pos.loc.translate_vec(v); }
    pub fn translated_vec(&self, v: &GpVec2d) -> Self { let mut r=*self; r.translate_vec(v); r }
}

impl Default for GpLin2d { fn default() -> Self { Self { pos: GpAx2d::default() } } }
