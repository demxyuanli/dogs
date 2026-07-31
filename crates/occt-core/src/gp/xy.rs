//! 2D coordinate {x,y}. Source: `gp_XY.hxx`
use crate::precision::RESOLUTION;
use crate::gp::mat2d::GpMat2d;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpXY { pub x: f64, pub y: f64 }

impl GpXY {
    #[inline] pub const fn zero() -> Self { Self { x: 0.0, y: 0.0 } }
    #[inline] pub const fn new(x: f64, y: f64) -> Self { Self { x, y } }
    #[inline] pub fn set_coord(&mut self, x: f64, y: f64) { self.x = x; self.y = y; }
    #[inline] pub fn set_x(&mut self, x: f64) { self.x = x; }
    #[inline] pub fn set_y(&mut self, y: f64) { self.y = y; }
    #[inline] pub const fn x(&self) -> f64 { self.x }
    #[inline] pub const fn y(&self) -> f64 { self.y }
    #[inline] pub const fn coords(&self) -> (f64, f64) { (self.x, self.y) }
    #[inline] pub fn modulus(&self) -> f64 { self.square_modulus().sqrt() }
    #[inline] pub const fn square_modulus(&self) -> f64 { self.x * self.x + self.y * self.y }
    #[inline] pub fn is_equal(&self, o: &Self, tol: f64) -> bool { (self.x-o.x).abs()<tol && (self.y-o.y).abs()<tol }
    #[inline] pub fn add(&mut self, o: &Self) { self.x+=o.x; self.y+=o.y; }
    #[inline] pub const fn added(&self, o: &Self) -> Self { Self { x: self.x+o.x, y: self.y+o.y } }
    #[inline] pub fn subtract(&mut self, o: &Self) { self.x-=o.x; self.y-=o.y; }
    #[inline] pub const fn subtracted(&self, o: &Self) -> Self { Self { x: self.x-o.x, y: self.y-o.y } }
    #[inline] pub fn multiply_scalar(&mut self, s: f64) { self.x*=s; self.y*=s; }
    #[inline] pub const fn multiplied_scalar(&self, s: f64) -> Self { Self { x: self.x*s, y: self.y*s } }
    #[inline] pub const fn multiplied(&self, s: f64) -> Self { self.multiplied_scalar(s) }
    #[inline] pub const fn multiply(&self, s: f64) -> Self { self.multiplied_scalar(s) }
    #[inline] pub fn divide(&mut self, s: f64) { self.x/=s; self.y/=s; }
    #[inline] pub const fn divided(&self, s: f64) -> Self { Self { x: self.x/s, y: self.y/s } }
    #[inline] pub fn multiply_xy(&mut self, o: &Self) { self.x*=o.x; self.y*=o.y; }
    #[inline] pub const fn multiplied_xy(&self, o: &Self) -> Self { Self { x: self.x*o.x, y: self.y*o.y } }
    #[inline] pub const fn dot(&self, o: &Self) -> f64 { self.x*o.x + self.y*o.y }
    #[inline] pub const fn crossed(&self, o: &Self) -> f64 { self.x*o.y - self.y*o.x }
    #[inline] pub fn multiply_mat2d(&mut self, m: &GpMat2d) { let ox=self.x; let oy=self.y; self.x=m.data[0][0]*ox+m.data[0][1]*oy; self.y=m.data[1][0]*ox+m.data[1][1]*oy; }
    #[inline] pub fn multiplied_mat2d(&self, m: &GpMat2d) -> Self { Self { x: m.data[0][0]*self.x+m.data[0][1]*self.y, y: m.data[1][0]*self.x+m.data[1][1]*self.y } }
    pub fn normalize(&mut self) -> Result<(), &'static str> { let d=self.modulus(); if d<=RESOLUTION { Err("zero") } else { self.divide(d); Ok(()) } }
    pub fn normalized(&self) -> Result<Self, &'static str> { let d=self.modulus(); if d<=RESOLUTION { Err("zero") } else { Ok(self.divided(d)) } }
    #[inline] pub fn reverse(&mut self) { self.x=-self.x; self.y=-self.y; }
    #[inline] pub const fn reversed(&self) -> Self { Self { x:-self.x, y:-self.y } }
    #[inline] pub fn set_linear_form_2(&mut self, a1:f64, xy1:&Self, a2:f64, xy2:&Self) { self.x=a1*xy1.x+a2*xy2.x; self.y=a1*xy1.y+a2*xy2.y; }
    #[inline] pub fn set_linear_form_add(&mut self, xy1:&Self, xy2:&Self) { self.x=xy1.x+xy2.x; self.y=xy1.y+xy2.y; }
    #[inline] pub fn set_linear_form_add_scaled(&mut self, a1:f64, xy1:&Self, xy2:&Self) { self.x=a1*xy1.x+xy2.x; self.y=a1*xy1.y+xy2.y; }
}

impl Default for GpXY { #[inline] fn default() -> Self { Self::zero() } }

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn dot_unit() { assert!((GpXY::new(1.0,0.0).dot(&GpXY::new(0.0,1.0))-0.0).abs()<1e-15); }
    #[test] fn cross_scalar() { assert!((GpXY::new(1.0,0.0).crossed(&GpXY::new(0.0,1.0))-1.0).abs()<1e-15); }
    #[test] fn cross_anti() { let a=GpXY::new(3.0,-2.0); let b=GpXY::new(1.0,4.0); assert!((a.crossed(&b)+b.crossed(&a)).abs()<1e-15); }
    #[test] fn norm() { let v=GpXY::new(3.0,4.0).normalized().unwrap(); assert!((v.modulus()-1.0).abs()<1e-15); }
}
