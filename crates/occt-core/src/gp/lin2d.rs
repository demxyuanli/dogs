//! 2D line. Source: `gp_Lin2d.hxx`
use crate::gp::{ax2d::GpAx2d, pnt2d::GpPnt2d, dir2d::GpDir2d, vec2d::GpVec2d, trsf2d::GpTrsf2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpLin2d { pub pos: GpAx2d }

impl GpLin2d {
    pub fn new(a: GpAx2d) -> Self { Self { pos: a } }
    pub fn from_pnt_dir(p: GpPnt2d, d: GpDir2d) -> Self { Self { pos: GpAx2d::new(p, d) } }
    #[inline] pub fn position(&self) -> &GpAx2d { &self.pos }
    pub fn set_position(&mut self, a: GpAx2d) { self.pos = a; }
    pub fn angle(&self, other: &Self) -> f64 { self.pos.vdir.angle(&other.pos.vdir) }
    pub fn distance(&self, p: &GpPnt2d) -> f64 {
        let dx = p.x() - self.pos.loc.x(); let dy = p.y() - self.pos.loc.y();
        let dot = dx * self.pos.vdir.x + dy * self.pos.vdir.y;
        (dx*dx + dy*dy - dot*dot).sqrt()
    }
    pub fn contains(&self, p: &GpPnt2d, tol: f64) -> bool { self.distance(p) <= tol }
    pub fn mirror_pnt(&mut self, p: &GpPnt2d) { let mut t=GpTrsf2d::identity(); t.set_mirror_pnt(p); t.transforms_xy(&mut self.pos.loc.coord); }
    pub fn mirrored_pnt(&self, p: &GpPnt2d) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax2d(&mut self, a: &GpAx2d) { let mut t=GpTrsf2d::identity(); t.set_mirror_ax2d(a); t.transforms_xy(&mut self.pos.loc.coord); }
    pub fn mirrored_ax2d(&self, a: &GpAx2d) -> Self { let mut r=*self; r.mirror_ax2d(a); r }
    pub fn rotate(&mut self, p: &GpPnt2d, angle: f64) { self.pos.loc.rotate(p, angle); }
    pub fn rotated(&self, p: &GpPnt2d, angle: f64) -> Self { let mut r=*self; r.rotate(p, angle); r }
    pub fn scale(&mut self, p: &GpPnt2d, s: f64) { self.pos.loc.scale(p, s); }
    pub fn scaled(&self, p: &GpPnt2d, s: f64) -> Self { let mut r=*self; r.scale(p, s); r }
    pub fn transform(&mut self, t: &GpTrsf2d) { t.transforms_xy(&mut self.pos.loc.coord); }
    pub fn transformed(&self, t: &GpTrsf2d) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self, v: &GpVec2d) { self.pos.loc.translate_vec(v); }
    pub fn translated_vec(&self, v: &GpVec2d) -> Self { let mut r=*self; r.translate_vec(v); r }
}

impl Default for GpLin2d { fn default() -> Self { Self { pos: GpAx2d::default() } } }
