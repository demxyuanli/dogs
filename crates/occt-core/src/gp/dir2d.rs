//! 2D unit direction. Source: `gp_Dir2d.hxx`
use crate::precision::RESOLUTION;
use crate::gp::{xy::GpXY, vec2d::GpVec2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpDir2d { pub x: f64, pub y: f64 }
impl GpDir2d {
    pub fn new(x:f64,y:f64) -> Result<Self,&'static str> { let sq=x*x+y*y; if sq>=(1.0-RESOLUTION)&&sq<=(1.0+RESOLUTION) { Ok(Self{x,y}) } else if sq<=RESOLUTION*RESOLUTION { Err("zero") } else { let d=sq.sqrt(); Ok(Self{x:x/d,y:y/d}) } }
    pub fn from_xy(c:&GpXY) -> Result<Self,&'static str> { Self::new(c.x,c.y) }
    pub fn from_vec2d(v:&GpVec2d) -> Result<Self,&'static str> { Self::new(v.x(),v.y()) }
    #[inline] pub const fn x(&self) -> f64 { self.x } #[inline] pub const fn y(&self) -> f64 { self.y }
    #[inline] pub const fn dot(&self,o:&Self) -> f64 { self.x*o.x+self.y*o.y }
    #[inline] pub const fn crossed(&self,o:&Self) -> f64 { self.x*o.y-self.y*o.x }
    pub fn angle(&self,o:&Self) -> f64 { self.dot(o).clamp(-1.0,1.0).acos() }
    pub fn is_equal(&self,o:&Self,tol:f64) -> bool { (self.x-o.x).abs()<tol && (self.y-o.y).abs()<tol }
    pub fn is_normal(&self,o:&Self,tol:f64) -> bool { (std::f64::consts::FRAC_PI_2-self.angle(o)).abs()<=tol }
    pub fn is_parallel(&self,o:&Self,tol:f64) -> bool { let a=self.angle(o); a<=tol||std::f64::consts::PI-a<=tol }
    pub fn is_opposite(&self,o:&Self,tol:f64) -> bool { std::f64::consts::PI-self.angle(o)<=tol }
    pub fn reverse(&mut self) { self.x=-self.x; self.y=-self.y; }
    pub fn reversed(&self) -> Self { Self{x:-self.x,y:-self.y} }
    /// `gp_Dir2d::Transform(gp_Trsf2d)` — vectorial part only (no translation).
    pub fn transform(&mut self, t: &crate::gp::trsf2d::GpTrsf2d) {
        let mut xy = GpXY::new(self.x, self.y);
        xy.multiply_mat2d(t.vectorial_part());
        if t.scale_factor() != 1.0 {
            xy.multiply_scalar(t.scale_factor());
        }
        if let Ok(d) = Self::from_xy(&xy) {
            *self = d;
        }
    }
    pub fn transformed(&self, t: &crate::gp::trsf2d::GpTrsf2d) -> Self {
        let mut r = *self;
        r.transform(t);
        r
    }
}
impl Default for GpDir2d { fn default() -> Self { Self { x: 1.0, y: 0.0 } } }
