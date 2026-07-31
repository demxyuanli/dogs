//! 2D axis. Source: `gp_Ax2d.hxx`
use crate::gp::{pnt2d::GpPnt2d, dir2d::GpDir2d};

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
}
impl Default for GpAx2d { fn default() -> Self { Self{loc:GpPnt2d::zero(),vdir:GpDir2d::default()} } }
