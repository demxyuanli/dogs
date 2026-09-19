//! 2D coordinate system. Source: `gp_Ax22d.hxx`
use crate::precision::ANGULAR;
use crate::gp::{ax2d::GpAx2d, pnt2d::GpPnt2d, dir2d::GpDir2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpAx22d { pub point: GpPnt2d, pub vxdir: GpDir2d, pub vydir: GpDir2d }
impl GpAx22d {
    pub fn standard() -> Self { Self{point:GpPnt2d::zero(),vxdir:GpDir2d::default(),vydir:GpDir2d{x:0.0,y:1.0}} }
    pub fn from_xdir(p:GpPnt2d,vx:GpDir2d) -> Self { Self{point:p,vxdir:vx,vydir:GpDir2d{x:-vx.y,y:vx.x}} }
    pub fn new(p:GpPnt2d,vx:GpDir2d,vy:GpDir2d) -> Result<Self,&'static str> { if !vx.is_normal(&vy,ANGULAR) { Err("not perpendicular") } else { Ok(Self{point:p,vxdir:vx,vydir:vy}) } }
    #[inline] pub const fn location(&self) -> &GpPnt2d { &self.point }
    #[inline] pub const fn x_direction(&self) -> &GpDir2d { &self.vxdir }
    #[inline] pub const fn y_direction(&self) -> &GpDir2d { &self.vydir }
    /// `gp_Ax22d::XAxis()` (`gp_Ax22d.hxx:151`): the origin plus the X
    /// direction, as a `gp_Ax2d`.
    #[inline] pub fn x_axis(&self) -> GpAx2d { GpAx2d::new(self.point, self.vxdir) }
    /// `gp_Ax22d::YAxis()` (`gp_Ax22d.hxx:157`).
    #[inline] pub fn y_axis(&self) -> GpAx2d { GpAx2d::new(self.point, self.vydir) }
    pub fn set_location(&mut self,p:GpPnt2d) { self.point=p; }
    pub fn set_x_direction(&mut self,vx:GpDir2d) { self.vxdir=vx; self.vydir=GpDir2d{x:-vx.y,y:vx.x}; }
    pub fn set_y_direction(&mut self,vy:GpDir2d) { self.vydir=vy; self.vxdir=GpDir2d{x:vy.y,y:-vy.x}; }
    pub fn x_reverse(&mut self) { self.vxdir.reverse(); }
    pub fn y_reverse(&mut self) { self.vydir.reverse(); }
}
impl Default for GpAx22d { fn default() -> Self { Self::standard() } }
