//! 2×2 matrix. Source: `gp_Mat2d.hxx`
use crate::precision::RESOLUTION;
use crate::gp::xy::GpXY;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpMat2d { pub data: [[f64; 2]; 2] }

impl GpMat2d {
    #[inline] pub const fn zero() -> Self { Self { data: [[0.0;2];2] } }
    #[inline] pub const fn identity() -> Self { Self { data: [[1.0,0.0],[0.0,1.0]] } }
    #[inline] pub const fn new(a11:f64,a12:f64,a21:f64,a22:f64) -> Self { Self { data: [[a11,a12],[a21,a22]] } }
    pub fn from_cols(c1:&GpXY, c2:&GpXY) -> Self { Self { data: [[c1.x,c2.x],[c1.y,c2.y]] } }
    #[inline] pub const fn value(&self, r:usize, c:usize) -> f64 { self.data[r][c] }
    #[inline] pub fn row(&self, r:usize) -> GpXY { GpXY::new(self.data[r][0], self.data[r][1]) }
    #[inline] pub fn column(&self, c:usize) -> GpXY { GpXY::new(self.data[0][c], self.data[1][c]) }
    #[inline] pub fn diagonal(&self) -> GpXY { GpXY::new(self.data[0][0], self.data[1][1]) }
    #[inline] pub fn set_identity(&mut self) { *self=Self::identity(); }
    #[inline] pub fn set_diagonal(&mut self, x1:f64, x2:f64) { self.data[0][0]=x1; self.data[1][1]=x2; }
    #[inline] pub fn set_scale(&mut self, s:f64) { self.set_diagonal(s,s); }
    pub fn set_rotation(&mut self, angle:f64) { let(s,c)=(angle.sin(),angle.cos()); self.data=[[c,-s],[s,c]]; }
    #[inline] pub const fn determinant(&self) -> f64 { self.data[0][0]*self.data[1][1]-self.data[1][0]*self.data[0][1] }
    #[inline] pub fn add(&mut self, o:&Self) { for r in 0..2 { for c in 0..2 { self.data[r][c]+=o.data[r][c]; } } }
    #[inline] pub fn added(&self, o:&Self) -> Self { let mut m=*self; m.add(o); m }
    #[inline] pub fn subtract(&mut self, o:&Self) { for r in 0..2 { for c in 0..2 { self.data[r][c]-=o.data[r][c]; } } }
    #[inline] pub fn subtracted(&self, o:&Self) -> Self { let mut m=*self; m.subtract(o); m }
    #[inline] pub fn multiply_scalar(&mut self, s:f64) { for r in 0..2 { for c in 0..2 { self.data[r][c]*=s; } } }
    #[inline] pub fn multiplied_scalar(&self, s:f64) -> Self { let mut m=*self; m.multiply_scalar(s); m }
    pub fn multiply(&mut self, o:&Self) { let a=self.data; let b=o.data; let t00=a[0][0]*b[0][0]+a[0][1]*b[1][0]; let t10=a[1][0]*b[0][0]+a[1][1]*b[1][0]; self.data[0][1]=a[0][0]*b[0][1]+a[0][1]*b[1][1]; self.data[1][1]=a[1][0]*b[0][1]+a[1][1]*b[1][1]; self.data[0][0]=t00; self.data[1][0]=t10; }
    #[inline] pub fn multiplied(&self, o:&Self) -> Self { let mut m=*self; m.multiply(o); m }
    pub fn divide(&mut self, s:f64) -> Result<(),&'static str> { if s.abs()<=RESOLUTION { Err("singular") } else { self.multiply_scalar(1.0/s); Ok(()) } }
    pub fn divided(&self, s:f64) -> Result<Self,&'static str> { let mut m=*self; m.divide(s)?; Ok(m) }
    pub fn invert(&mut self) -> Result<(),&'static str> { let det=self.determinant(); if det.abs()<=RESOLUTION { Err("singular") } else { let inv=1.0/det; let a=self.data; self.data=[[a[1][1]*inv,-a[0][1]*inv],[-a[1][0]*inv,a[0][0]*inv]]; Ok(()) } }
    pub fn inverted(&self) -> Result<Self,&'static str> { let mut m=*self; m.invert()?; Ok(m) }
    #[inline] pub fn transpose(&mut self) { let t=self.data[0][1]; self.data[0][1]=self.data[1][0]; self.data[1][0]=t; }
    #[inline] pub fn transposed(&self) -> Self { let mut m=*self; m.transpose(); m }
}

impl Default for GpMat2d { #[inline] fn default() -> Self { Self::identity() } }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::xy::GpXY;
    #[test] fn identity_det() { assert!((GpMat2d::identity().determinant()-1.0).abs()<1e-15); }
    #[test] fn rot90() { let mut m=GpMat2d::identity(); m.set_rotation(std::f64::consts::FRAC_PI_2); let v=GpXY::new(1.0,0.0).multiplied_mat2d(&m); assert!((v.x-0.0).abs()<1e-15); assert!((v.y-1.0).abs()<1e-15); }
    #[test] fn invert_roundtrip() { let m=GpMat2d::new(3.0,1.0,2.0,4.0); let inv=m.inverted().unwrap(); let p=m.multiplied(&inv); assert!((p.data[0][0]-1.0).abs()<1e-14); assert!(p.data[0][1].abs()<1e-14); assert!((p.data[1][1]-1.0).abs()<1e-14); }
}
