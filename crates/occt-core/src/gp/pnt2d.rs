//! 2D point. Source: `gp_Pnt2d.hxx`
use crate::gp::{xy::GpXY, vec2d::GpVec2d, trsf2d::GpTrsf2d, ax2d::GpAx2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpPnt2d { pub coord: GpXY }
impl GpPnt2d {
    #[inline] pub const fn zero() -> Self { Self { coord: GpXY::zero() } }
    #[inline] pub const fn new(x:f64,y:f64) -> Self { Self { coord: GpXY::new(x,y) } }
    #[inline] pub const fn from_xy(coord: GpXY) -> Self { Self { coord } }
    #[inline] pub const fn x(&self) -> f64 { self.coord.x } #[inline] pub const fn y(&self) -> f64 { self.coord.y }
    pub fn set_x(&mut self,x:f64) { self.coord.x=x; } pub fn set_y(&mut self,y:f64) { self.coord.y=y; }
    pub fn set_coord(&mut self,x:f64,y:f64) { self.coord.set_coord(x,y); }
    #[inline] pub const fn xy(&self) -> &GpXY { &self.coord }
    pub fn distance(&self,o:&Self) -> f64 { self.coord.subtracted(&o.coord).modulus() }
    pub fn square_distance(&self,o:&Self) -> f64 { self.coord.subtracted(&o.coord).square_modulus() }
    pub fn is_equal(&self,o:&Self,tol:f64) -> bool { self.distance(o)<=tol }
    pub fn mirror_pnt(&mut self,p:&Self) { let mut t=GpTrsf2d::identity(); t.set_mirror_pnt(p); t.transforms_xy(&mut self.coord); }
    pub fn mirrored_pnt(&self,p:&Self) -> Self { let mut r=*self; r.mirror_pnt(p); r }
    pub fn mirror_ax2d(&mut self,a:&GpAx2d) { let mut t=GpTrsf2d::identity(); t.set_mirror_ax2d(a); t.transforms_xy(&mut self.coord); }
    pub fn mirrored_ax2d(&self,a:&GpAx2d) -> Self { let mut r=*self; r.mirror_ax2d(a); r }
    pub fn rotate(&mut self,p:&Self,angle:f64) { let mut t=GpTrsf2d::identity(); t.set_rotation(p,angle); t.transforms_xy(&mut self.coord); }
    pub fn rotated(&self,p:&Self,angle:f64) -> Self { let mut r=*self; r.rotate(p,angle); r }
    pub fn scale(&mut self,p:&Self,s:f64) { let mut xy=p.coord; xy.multiply_scalar(1.0-s); self.coord.multiply_scalar(s); self.coord.add(&xy); }
    pub fn scaled(&self,p:&Self,s:f64) -> Self { let mut r=*self; r.scale(p,s); r }
    pub fn transform(&mut self,t:&GpTrsf2d) { t.transforms_xy(&mut self.coord); }
    pub fn transformed(&self,t:&GpTrsf2d) -> Self { let mut r=*self; r.transform(t); r }
    pub fn translate_vec(&mut self,v:&GpVec2d) { self.coord.add(&v.xy()); }
    pub fn translated_vec(&self,v:&GpVec2d) -> Self { let mut r=*self; r.translate_vec(v); r }
    pub fn translate_pnts(&mut self,p1:&Self,p2:&Self) { self.coord.add(&p2.coord); self.coord.subtract(&p1.coord); }
    pub fn translated_pnts(&self,p1:&Self,p2:&Self) -> Self { let mut r=*self; r.translate_pnts(p1,p2); r }
}
impl Default for GpPnt2d { fn default() -> Self { Self::zero() } }
