//! 2D coordinate system. Source: `gp_Ax22d.hxx`
use crate::precision::ANGULAR;
use crate::gp::{ax2d::GpAx2d, pnt2d::GpPnt2d, dir2d::GpDir2d, trsf2d::GpTrsf2d};

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
    /// `gp_Ax22d::Transform` (`gp_Ax22d.hxx:360-367`): transform the location
    /// and both (unit) directions by the trsf.
    pub fn transform(&mut self, t: &GpTrsf2d) {
        self.point.transform(t);
        self.vxdir.transform(t);
        self.vydir.transform(t);
    }
    /// `gp_Ax22d::Rotate` (`gp_Ax22d.hxx:335-341`): rotate the location and
    /// both directions.
    pub fn rotate(&mut self, p: &GpPnt2d, ang: f64) {
        self.point.rotate(p, ang);
        self.vxdir.rotate(ang);
        self.vydir.rotate(ang);
    }
    /// `gp_Ax22d::Scale` (`gp_Ax22d.hxx:346-354`): scale the location and, for
    /// a negative factor, reverse both directions.
    pub fn scale(&mut self, p: &GpPnt2d, s: f64) {
        self.point.scale(p, s);
        if s < 0.0 {
            self.vxdir.reverse();
            self.vydir.reverse();
        }
    }
    /// `gp_Ax22d::Mirror(const gp_Pnt2d&)` (`gp_Ax22d.cxx:29-36`): mirror the
    /// location and reverse both directions.
    pub fn mirror_pnt(&mut self, p: &GpPnt2d) {
        self.point.mirror_pnt(p);
        self.vxdir.reverse();
        self.vydir.reverse();
    }
    /// `gp_Ax22d::Mirror(const gp_Ax2d&)` (`gp_Ax22d.cxx:43-52`): mirror the
    /// location and both directions across the axis.
    pub fn mirror_ax2d(&mut self, a: &GpAx2d) {
        self.vydir.mirror_ax2d(a);
        self.vxdir.mirror_ax2d(a);
        self.point.mirror_ax2d(a);
    }
}
impl Default for GpAx22d { fn default() -> Self { Self::standard() } }
