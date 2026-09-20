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
    /// `gp_Dir2d::Angle` (`gp_Dir2d.cxx:26-63`): the **signed** angle in
    /// `]-PI; PI]`, using `acos` inside 45 degrees and `asin` outside for
    /// precision. The previous body was `acos(self.dot(o))` — unsigned — which
    /// also made the three predicates below deviate (audit A10 / task T-76).
    pub fn angle(&self, o: &Self) -> f64 {
        const COS_45: f64 = 0.70710678118655;
        let cosinus = self.dot(o);
        let sinus = self.crossed(o);
        if cosinus > -COS_45 && cosinus < COS_45 {
            if sinus > 0.0 {
                cosinus.acos()
            } else {
                -cosinus.acos()
            }
        } else if cosinus > 0.0 {
            sinus.asin()
        } else if sinus > 0.0 {
            std::f64::consts::PI - sinus.asin()
        } else {
            -std::f64::consts::PI - sinus.asin()
        }
    }
    pub fn is_equal(&self,o:&Self,tol:f64) -> bool { (self.x-o.x).abs()<tol && (self.y-o.y).abs()<tol }
    /// `gp_Dir2d::IsNormal` (`gp_Dir2d.hxx:393-406`): `abs(PI/2 - abs(Angle)) <= tol`.
    pub fn is_normal(&self,o:&Self,tol:f64) -> bool { (std::f64::consts::FRAC_PI_2-self.angle(o).abs()).abs()<=tol }
    /// `gp_Dir2d::IsParallel` (`gp_Dir2d.hxx:422-430`): `abs(Angle) <= tol || PI - abs(Angle) <= tol`.
    pub fn is_parallel(&self,o:&Self,tol:f64) -> bool { let a=self.angle(o).abs(); a<=tol||std::f64::consts::PI-a<=tol }
    /// `gp_Dir2d::IsOpposite` (`gp_Dir2d.hxx:410-418`): `PI - abs(Angle) <= tol`.
    pub fn is_opposite(&self,o:&Self,tol:f64) -> bool { std::f64::consts::PI-self.angle(o).abs()<=tol }
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
