//! 2D axis. Source: `gp_Ax2d.hxx`
use crate::gp::{pnt2d::GpPnt2d, dir2d::GpDir2d, trsf2d::GpTrsf2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpAx2d { pub loc: GpPnt2d, pub vdir: GpDir2d }
impl GpAx2d {
    pub fn new(loc:GpPnt2d,dir:GpDir2d) -> Self { Self{loc,vdir:dir} }
    #[inline] pub const fn location(&self) -> &GpPnt2d { &self.loc }
    #[inline] pub const fn direction(&self) -> &GpDir2d { &self.vdir }
    pub fn set_location(&mut self,l:GpPnt2d) { self.loc=l; }
    pub fn set_direction(&mut self,d:GpDir2d) { self.vdir=d; }
    pub fn is_coaxial(&self,o:&Self,atol:f64,ltol:f64) -> bool { self.vdir.is_parallel(&o.vdir,atol) && self.loc.distance(&o.loc)<=ltol }
    pub fn reverse(&mut self) { self.vdir.reverse(); }
    pub fn reversed(&self) -> Self { Self{loc:self.loc,vdir:self.vdir.reversed()} }
    /// `gp_Ax2d::Rotate` (`gp_Ax2d.hxx:155-159`): rotate the location and the
    /// direction.
    pub fn rotate(&mut self, p: &GpPnt2d, ang: f64) {
        self.loc.rotate(p, ang);
        self.vdir.rotate(ang);
    }
    /// `gp_Ax2d::Scale` (`gp_Ax2d.cxx:51-58`): scale the location and, for a
    /// negative factor, reverse the direction.
    pub fn scale(&mut self, p: &GpPnt2d, s: f64) {
        self.loc.scale(p, s);
        if s < 0.0 { self.vdir.reverse(); }
    }
    /// `gp_Ax2d::Mirror(const gp_Pnt2d&)` (`gp_Ax2d.cxx:60-64`): mirror the
    /// location and reverse the direction.
    pub fn mirror_pnt(&mut self, p: &GpPnt2d) {
        self.loc.mirror_pnt(p);
        self.vdir.reverse();
    }
    /// `gp_Ax2d::Mirror(const gp_Ax2d&)` (`gp_Ax2d.cxx:73-77`): mirror the
    /// location and reflect the direction across `a`'s direction.
    pub fn mirror_ax2d(&mut self, a: &Self) {
        self.loc.mirror_ax2d(a);
        self.vdir.mirror_dir2d(&a.vdir);
    }
    /// `gp_Ax2d::Transform` (`gp_Ax2d.hxx:182-186`): transform the location and
    /// the (unit) direction.
    pub fn transform(&mut self, t: &GpTrsf2d) {
        self.loc.transform(t);
        self.vdir.transform(t);
    }
}
impl Default for GpAx2d { fn default() -> Self { Self{loc:GpPnt2d::zero(),vdir:GpDir2d::default()} } }
