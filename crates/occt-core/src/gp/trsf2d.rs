//! 2D affine transform. Source: `gp_Trsf2d.hxx`
use crate::precision::RESOLUTION;
use crate::gp::{mat2d::GpMat2d, xy::GpXY, trsf_form::TrsfForm, pnt2d::GpPnt2d, vec2d::GpVec2d, ax2d::GpAx2d};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpTrsf2d { pub scale: f64, pub shape: TrsfForm, pub matrix: GpMat2d, pub loc: GpXY }
impl GpTrsf2d {
    pub fn identity() -> Self { Self { scale:1.0, shape:TrsfForm::Identity, matrix:GpMat2d::identity(), loc:GpXY::zero() } }
    #[inline] pub const fn scale_factor(&self) -> f64 { self.scale }
    #[inline] pub const fn form(&self) -> TrsfForm { self.shape }
    #[inline] pub const fn translation_part(&self) -> &GpXY { &self.loc }
    #[inline] pub const fn vectorial_part(&self) -> &GpMat2d { &self.matrix }
    #[inline] pub fn is_negative(&self) -> bool { self.matrix.determinant()<0.0 }
    pub fn set_identity(&mut self) { *self=Self::identity(); }
    pub fn set_translation_vec(&mut self,v:&GpVec2d) { self.shape=TrsfForm::Translation; self.scale=1.0; self.matrix.set_identity(); self.loc=v.xy(); }
    pub fn set_translation_pnts(&mut self,p1:&GpPnt2d,p2:&GpPnt2d) { self.shape=TrsfForm::Translation; self.scale=1.0; self.matrix.set_identity(); self.loc=p2.xy().subtracted(p1.xy()); }
    pub fn set_mirror_pnt(&mut self,p:&GpPnt2d) { self.shape=TrsfForm::PntMirror; self.scale=-1.0; self.matrix.set_identity(); self.loc=p.xy().multiplied_scalar(2.0); }
    pub fn set_mirror_ax2d(&mut self,ax:&GpAx2d) { self.shape=TrsfForm::Ax1Mirror; self.scale=-1.0; let d=ax.direction(); let p=ax.location().xy(); let(dx,dy)=(d.x,d.y); self.matrix=GpMat2d::new(2.0*dx*dx-1.0,2.0*dx*dy,2.0*dy*dx,2.0*dy*dy-1.0); self.loc=*p; self.loc.multiply_mat2d(&self.matrix); self.loc.reverse(); self.loc.add(p); }
    pub fn set_rotation(&mut self,p:&GpPnt2d,angle:f64) { self.shape=TrsfForm::Rotation; self.scale=1.0; self.loc=p.xy().reversed(); self.matrix.set_rotation(angle); self.loc.multiply_mat2d(&self.matrix); self.loc.add(p.xy()); }
    pub fn set_scale(&mut self,p:&GpPnt2d,s:f64) -> Result<(),&'static str> { if s.abs()<=RESOLUTION { Err("too small") } else { self.shape=TrsfForm::Scale; self.scale=s; self.matrix.set_identity(); self.loc=p.xy().multiplied_scalar(1.0-s); Ok(()) } }
    pub fn transforms_xy(&self,coord:&mut GpXY) { coord.multiply_mat2d(&self.matrix); if self.scale!=1.0 { coord.multiply_scalar(self.scale); } coord.add(&self.loc); }
    pub fn multiply(&mut self,o:&Self) { let mut nl=o.loc; nl.multiply_mat2d(&self.matrix); if self.scale!=1.0 { nl.multiply_scalar(self.scale); } nl.add(&self.loc); self.loc=nl; self.matrix.multiply(&o.matrix); self.scale*=o.scale; if self.shape==TrsfForm::Identity { self.shape=o.shape; } else if o.shape!=TrsfForm::Identity { self.shape=TrsfForm::CompoundTrsf; } }
    pub fn multiplied(&self,o:&Self) -> Self { let mut r=*self; r.multiply(o); r }
    pub fn invert(&mut self) -> Result<(),&'static str> { self.matrix.invert()?; self.loc.reverse(); self.loc.multiply_mat2d(&self.matrix); if self.scale!=1.0 { self.loc.divide(self.scale); } self.scale=1.0/self.scale; Ok(()) }
    pub fn inverted(&self) -> Result<Self,&'static str> { let mut r=*self; r.invert()?; Ok(r) }
}
impl Default for GpTrsf2d { fn default() -> Self { Self::identity() } }
