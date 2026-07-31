//! Infinite line in 3D. Source: `gp_Lin.hxx`
use crate::gp::{ax1::GpAx1, pnt::GpPnt, dir::GpDir, vec::GpVec, trsf::GpTrsf, ax2::GpAx2};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpLin { pub pos: GpAx1 }

impl GpLin {
    pub fn new(a1: GpAx1) -> Self { Self { pos: a1 } }
    pub fn from_pnt_dir(p: GpPnt, d: GpDir) -> Self { Self { pos: GpAx1::new(p, d) } }
    #[inline] pub fn position(&self) -> &GpAx1 { &self.pos }
    #[inline] pub fn location(&self) -> GpPnt { *self.pos.location() }
    #[inline] pub fn direction(&self) -> GpDir { *self.pos.direction() }
    pub fn set_location(&mut self, p: GpPnt) { self.pos.set_location(p); }
    pub fn set_direction(&mut self, d: GpDir) { self.pos.set_direction(d); }
    pub fn set_position(&mut self, a1: GpAx1) { self.pos = a1; }
    pub fn angle(&self, other: &Self) -> f64 { self.pos.vdir.angle(&other.pos.vdir) }
    pub fn contains(&self, p: &GpPnt, tol: f64) -> bool { self.distance(p) <= tol }
    pub fn distance(&self, p: &GpPnt) -> f64 {
        let v = p.coord.subtracted(&self.pos.loc.coord);
        let dxyz = self.pos.vdir.xyz();
        (v.subtracted(&dxyz.multiplied(v.dot(dxyz)))).modulus()
    }
    pub fn square_distance(&self, p: &GpPnt) -> f64 { let d = self.distance(p); d * d }
    pub fn normal(&self, p: &GpPnt) -> GpLin {
        let v = p.coord.subtracted(&self.pos.loc.coord);
        let dxyz = self.pos.vdir.xyz();
        let proj = dxyz.multiplied(v.dot(dxyz));
        let n = v.subtracted(&proj);
        GpLin::from_pnt_dir(*p, GpDir::from_xyz(&n).unwrap_or(GpDir::default_dir()))
    }
    pub fn mirror_pnt(&mut self, p: &GpPnt) { let mut t=GpTrsf::identity(); t.set_mirror_pnt(p); self.pos.loc.transform(&t); }
    pub fn mirrored_pnt(&self, p: &GpPnt) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax1(&mut self, a1: &GpAx1) { let mut t=GpTrsf::identity(); t.set_mirror_ax1(a1); self.pos.loc.transform(&t); }
    pub fn mirrored_ax1(&self, a1: &GpAx1) -> Self { let mut r=*self; r.mirror_ax1(a1); r }
    pub fn mirror_ax2(&mut self, a2: &GpAx2) { let mut t=GpTrsf::identity(); t.set_mirror_ax2(a2); self.pos.loc.transform(&t); }
    pub fn mirrored_ax2(&self, a2: &GpAx2) -> Self { let mut r=*self; r.mirror_ax2(a2); r }
    pub fn rotate(&mut self, a1: &GpAx1, angle: f64) { self.pos.loc.rotate(a1, angle); }
    pub fn rotated(&self, a1: &GpAx1, angle: f64) -> Self { let mut r=*self; r.rotate(a1, angle); r }
    pub fn scale(&mut self, p: &GpPnt, s: f64) { self.pos.loc.scale(p, s); }
    pub fn scaled(&self, p: &GpPnt, s: f64) -> Self { let mut r=*self; r.scale(p, s); r }
    pub fn transform(&mut self, t: &GpTrsf) { self.pos.loc.transform(t); }
    pub fn transformed(&self, t: &GpTrsf) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self, v: &GpVec) { self.pos.loc.translate_vec(v); }
    pub fn translated_vec(&self, v: &GpVec) -> Self { let mut r=*self; r.translate_vec(v); r }
    pub fn translate_pnts(&mut self, p1: &GpPnt, p2: &GpPnt) { self.pos.loc.translate_pnts(p1, p2); }
    pub fn translated_pnts(&self, p1: &GpPnt, p2: &GpPnt) -> Self { let mut r=*self; r.translate_pnts(p1, p2); r }
}

impl Default for GpLin { fn default() -> Self { Self { pos: GpAx1::default() } } }
