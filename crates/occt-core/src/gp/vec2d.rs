//! 2D vector. Source: `gp_Vec2d.hxx`
use crate::precision::RESOLUTION;
use crate::gp::{xy::GpXY, mat2d::GpMat2d, dir2d::GpDir2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpVec2d { pub coord: GpXY }
impl GpVec2d {
    #[inline] pub const fn zero() -> Self { Self { coord: GpXY::zero() } }
    #[inline] pub const fn new(x:f64,y:f64) -> Self { Self { coord: GpXY::new(x,y) } }
    #[inline] pub const fn from_xy(c: GpXY) -> Self { Self { coord: c } }
    pub fn from_dir2d(d: &GpDir2d) -> Self { Self { coord: GpXY::new(d.x,d.y) } }
    #[inline] pub const fn x(&self) -> f64 { self.coord.x } #[inline] pub const fn y(&self) -> f64 { self.coord.y }
    pub fn set_coord(&mut self,x:f64,y:f64) { self.coord.set_coord(x,y); }
    #[inline] pub const fn xy(&self) -> GpXY { self.coord }
    pub fn magnitude(&self) -> f64 { self.coord.modulus() }
    pub const fn square_magnitude(&self) -> f64 { self.coord.square_modulus() }
    pub fn is_equal(&self,o:&Self,tol:f64) -> bool { self.coord.is_equal(&o.coord,tol) }
    pub fn add(&mut self,o:&Self) { self.coord.add(&o.coord); }
    pub fn added(&self,o:&Self) -> Self { Self{coord:self.coord.added(&o.coord)} }
    pub fn subtract(&mut self,o:&Self) { self.coord.subtract(&o.coord); }
    pub fn subtracted(&self,o:&Self) -> Self { Self{coord:self.coord.subtracted(&o.coord)} }
    pub fn multiply_scalar(&mut self,s:f64) { self.coord.multiply_scalar(s); }
    pub fn multiplied_scalar(&self,s:f64) -> Self { Self{coord:self.coord.multiplied_scalar(s)} }
    pub fn divide(&mut self,s:f64) { self.coord.divide(s); }
    pub fn divided(&self,s:f64) -> Self { Self{coord:self.coord.divided(s)} }
    pub const fn dot(&self,o:&Self) -> f64 { self.coord.dot(&o.coord) }
    pub const fn crossed(&self,o:&Self) -> f64 { self.coord.crossed(&o.coord) }
    pub fn normalize(&mut self) -> Result<(),&'static str> { self.coord.normalize() }
    pub fn normalized(&self) -> Result<Self,&'static str> { Ok(Self{coord:self.coord.normalized()?}) }
    pub fn reverse(&mut self) { self.coord.reverse(); }
    pub fn reversed(&self) -> Self { Self{coord:self.coord.reversed()} }
    /// Signed angle from this vector to `o`, in `(-PI, PI]`.
    /// Source: `gp_Vec2d::Angle`.
    pub fn angle(&self, o: &Self) -> f64 {
        let m1 = self.coord.modulus();
        let m2 = o.coord.modulus();
        if m1 <= RESOLUTION || m2 <= RESOLUTION {
            return 0.0;
        }
        let d = m1 * m2;
        let cosinus = (self.coord.dot(&o.coord) / d).clamp(-1.0, 1.0);
        let sinus = self.coord.crossed(&o.coord) / d;
        const COS_45: f64 = std::f64::consts::FRAC_1_SQRT_2;
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
    pub fn is_normal(&self, o: &Self, tol: f64) -> bool {
        (std::f64::consts::FRAC_PI_2 - self.angle(o).abs()).abs() <= tol
    }
    pub fn is_parallel(&self, o: &Self, tol: f64) -> bool {
        let a = self.angle(o).abs();
        a <= tol || std::f64::consts::PI - a <= tol
    }
    pub fn is_opposite(&self, o: &Self, tol: f64) -> bool {
        std::f64::consts::PI - self.angle(o).abs() <= tol
    }
    pub fn multiply_mat2d(&mut self,m:&GpMat2d) { self.coord.multiply_mat2d(m); }
    pub fn multiplied_mat2d(&self,m:&GpMat2d) -> Self { Self{coord:self.coord.multiplied_mat2d(m)} }
    pub fn set_linear_form_2(&mut self,a1:f64,v1:&Self,a2:f64,v2:&Self) { self.coord.set_linear_form_2(a1,&v1.coord,a2,&v2.coord); }
    pub fn set_linear_form_add(&mut self,v1:&Self,v2:&Self) { self.coord.set_linear_form_add(&v1.coord,&v2.coord); }
    pub fn set_linear_form_add_scaled(&mut self,a1:f64,v1:&Self,v2:&Self) { self.coord.set_linear_form_add_scaled(a1,&v1.coord,&v2.coord); }
}
impl Default for GpVec2d { fn default() -> Self { Self::zero() } }
