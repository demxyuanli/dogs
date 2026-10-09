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
    /// `gp_Trsf2d::HVectorialPart()` (`gp_Trsf2d.hxx`): the 2x2 matrix without
    /// the scale factor.
    #[inline] pub const fn h_vectorial_part(&self) -> &GpMat2d { &self.matrix }
    #[inline] pub fn is_negative(&self) -> bool { self.matrix.determinant()<0.0 }
    pub fn set_identity(&mut self) { *self=Self::identity(); }
    pub fn set_translation_vec(&mut self,v:&GpVec2d) { self.shape=TrsfForm::Translation; self.scale=1.0; self.matrix.set_identity(); self.loc=v.xy(); }
    pub fn set_translation_pnts(&mut self,p1:&GpPnt2d,p2:&GpPnt2d) { self.shape=TrsfForm::Translation; self.scale=1.0; self.matrix.set_identity(); self.loc=p2.xy().subtracted(p1.xy()); }
    pub fn set_mirror_pnt(&mut self,p:&GpPnt2d) { self.shape=TrsfForm::PntMirror; self.scale=-1.0; self.matrix.set_identity(); self.loc=p.xy().multiplied_scalar(2.0); }
    /// `gp_Trsf2d::SetMirror(const gp_Ax2d&)` (`gp_Trsf2d.cxx:31-46`): the
    /// vectorial matrix alone is the negated reflection, the sign being carried
    /// by `scale = -1` (so `VectorialPart = matrix*scale` is the reflection).
    pub fn set_mirror_ax2d(&mut self,ax:&GpAx2d) { self.shape=TrsfForm::Ax1Mirror; self.scale=-1.0; let d=ax.direction(); let (dx,dy)=(d.x,d.y); let(x0,y0)=(ax.location().x(),ax.location().y()); self.matrix=GpMat2d::new(1.0-2.0*dx*dx,-2.0*dx*dy,-2.0*dy*dx,1.0-2.0*dy*dy); self.loc=GpXY::new(-2.0*((dx*dx-1.0)*x0+dx*dy*y0),-2.0*(dx*dy*x0+(dy*dy-1.0)*y0)); }
    pub fn set_rotation(&mut self,p:&GpPnt2d,angle:f64) { self.shape=TrsfForm::Rotation; self.scale=1.0; self.loc=p.xy().reversed(); self.matrix.set_rotation(angle); self.loc.multiply_mat2d(&self.matrix); self.loc.add(p.xy()); }
    /// `gp_Trsf2d::SetTransformation(const gp_Ax2d&, const gp_Ax2d&)` (`gp_Trsf2d.cxx:48-70`):
    /// change of basis from `from_system1` to `to_system2`.
    pub fn set_transformation(&mut self, from_system1: &GpAx2d, to_system2: &GpAx2d) {
        self.shape = TrsfForm::CompoundTrsf;
        self.scale = 1.0;
        // matrix from XOY to theSystem2
        let v1 = GpXY::new(to_system2.vdir.x, to_system2.vdir.y);
        let v2 = GpXY::new(-v1.y, v1.x);
        self.matrix = GpMat2d::from_cols(&v1, &v2);
        self.loc = *to_system2.loc.xy();
        self.matrix.transpose();
        self.loc.multiply_mat2d(&self.matrix);
        self.loc.reverse();
        // matrix from theSystem1 to XOY
        let v3 = GpXY::new(from_system1.vdir.x, from_system1.vdir.y);
        let v4 = GpXY::new(-v3.y, v3.x);
        let ma1 = GpMat2d::from_cols(&v3, &v4);
        let mut ma1loc = *from_system1.loc.xy();
        // matrix * MA1 => fromSystem1 -> toSystem2
        ma1loc.multiply_mat2d(&self.matrix);
        self.loc.add(&ma1loc);
        self.matrix.multiply(&ma1);
    }
    pub fn set_scale(&mut self,p:&GpPnt2d,s:f64) -> Result<(),&'static str> { if s.abs()<=RESOLUTION { Err("too small") } else { self.shape=TrsfForm::Scale; self.scale=s; self.matrix.set_identity(); self.loc=p.xy().multiplied_scalar(1.0-s); Ok(()) } }
    /// `gp_Trsf2d::SetTransformation(const gp_Ax2d&)` (`gp_Trsf2d.cxx:72-83`):
    /// world to the local frame of `a`.
    pub fn set_transformation_ax2d(&mut self, a: &GpAx2d) {
        self.shape = TrsfForm::CompoundTrsf;
        self.scale = 1.0;
        let v1 = GpXY::new(a.vdir.x, a.vdir.y);
        let v2 = GpXY::new(-v1.y, v1.x);
        self.matrix.set_col(1, &v1);
        self.matrix.set_col(2, &v2);
        self.loc = *a.loc.xy();
        self.matrix.transpose();
        self.loc.multiply_mat2d(&self.matrix);
        self.loc.reverse();
    }

    /// `gp_Trsf2d::SetTranslationPart(const gp_Vec2d&)` (`gp_Trsf2d.cxx:86-117`).
    pub fn set_translation_part(&mut self, v: &GpVec2d) {
        self.loc = v.xy();
        let loc_null = self.loc.x().abs() <= RESOLUTION && self.loc.y().abs() <= RESOLUTION;
        if loc_null {
            match self.shape {
                TrsfForm::Identity
                | TrsfForm::PntMirror
                | TrsfForm::Scale
                | TrsfForm::Rotation
                | TrsfForm::Ax1Mirror => {}
                TrsfForm::Translation => self.shape = TrsfForm::Identity,
                _ => self.shape = TrsfForm::CompoundTrsf,
            }
        } else {
            match self.shape {
                TrsfForm::Translation | TrsfForm::Scale | TrsfForm::PntMirror => {}
                TrsfForm::Identity => self.shape = TrsfForm::Translation,
                _ => self.shape = TrsfForm::CompoundTrsf,
            }
        }
    }

    /// `gp_Trsf2d::SetScaleFactor(const double)` (`gp_Trsf2d.cxx:120-196`).
    pub fn set_scale_factor(&mut self, s: f64) {
        if s == 1.0 {
            let x = self.loc.x().abs();
            let y = self.loc.y().abs();
            if x <= RESOLUTION && y <= RESOLUTION {
                match self.shape {
                    TrsfForm::Identity | TrsfForm::Rotation => {}
                    TrsfForm::Scale => self.shape = TrsfForm::Identity,
                    TrsfForm::PntMirror => self.shape = TrsfForm::Translation,
                    _ => self.shape = TrsfForm::CompoundTrsf,
                }
            } else {
                match self.shape {
                    TrsfForm::Identity | TrsfForm::Rotation | TrsfForm::Scale => {}
                    TrsfForm::PntMirror => self.shape = TrsfForm::Translation,
                    _ => self.shape = TrsfForm::CompoundTrsf,
                }
            }
        } else if s == -1.0 {
            match self.shape {
                TrsfForm::PntMirror | TrsfForm::Ax1Mirror => {}
                TrsfForm::Identity | TrsfForm::Scale => self.shape = TrsfForm::PntMirror,
                _ => self.shape = TrsfForm::CompoundTrsf,
            }
        } else {
            match self.shape {
                TrsfForm::Scale => {}
                TrsfForm::Identity | TrsfForm::Translation | TrsfForm::PntMirror => {
                    self.shape = TrsfForm::Scale;
                }
                _ => self.shape = TrsfForm::CompoundTrsf,
            }
        }
        self.scale = s;
    }

    /// `gp_Trsf2d::VectorialPart()` (`gp_Trsf2d.cxx:198-213`): the matrix with
    /// the scale factor folded in.
    pub fn vectorial_part(&self) -> GpMat2d {
        if self.scale == 1.0 { return self.matrix; }
        let mut m = self.matrix;
        if self.shape == TrsfForm::Scale || self.shape == TrsfForm::PntMirror {
            m.set_diagonal(self.matrix.value(0,0) * self.scale, self.matrix.value(1,1) * self.scale);
        } else {
            m.multiply_scalar(self.scale);
        }
        m
    }

    /// `gp_Trsf2d::RotationPart()` (`gp_Trsf2d.cxx:216-219`).
    pub fn rotation_part(&self) -> f64 {
        self.matrix.value(1,0).atan2(self.matrix.value(0,0))
    }

    /// `gp_Trsf2d::Power(const int)` (`gp_Trsf2d.cxx:384-548`). `IsOdd(n)` is
    /// `n % 2 == 1` (`Standard_Integer.hxx:39-42`).
    pub fn power(&mut self, n: i32) -> Result<(),&'static str> {
        use TrsfForm::*;
        let is_odd = |v: i32| v % 2 == 1;
        if self.shape == Identity { return Ok(()); }
        if n == 0 {
            self.scale = 1.0;
            self.shape = Identity;
            self.matrix.set_identity();
            self.loc = GpXY::zero();
        } else if n == 1 {
            // already this^1
        } else if n == -1 {
            self.invert()?;
        } else {
            if n < 0 { self.invert()?; }
            if self.shape == Translation {
                let mut npower = n.abs() - 1;
                let mut temploc = self.loc;
                loop {
                    if is_odd(npower) { self.loc.add(&temploc); }
                    if npower == 1 { break; }
                    temploc = temploc.added(&temploc);
                    npower /= 2;
                }
            } else if self.shape == Scale {
                let mut npower = n.abs() - 1;
                let mut temploc = self.loc;
                let mut tempscale = self.scale;
                loop {
                    if is_odd(npower) {
                        self.loc.add(&temploc.multiplied(self.scale));
                        self.scale = self.scale * tempscale;
                    }
                    if npower == 1 { break; }
                    temploc = temploc.added(&temploc.multiplied(tempscale));
                    tempscale = tempscale * tempscale;
                    npower /= 2;
                }
            } else if self.shape == Rotation {
                let mut npower = n.abs() - 1;
                let mut tempmatrix = self.matrix;
                if self.loc.x() == 0.0 && self.loc.y() == 0.0 {
                    loop {
                        if is_odd(npower) { self.matrix.multiply(&tempmatrix); }
                        if npower == 1 { break; }
                        tempmatrix = tempmatrix.multiplied(&tempmatrix);
                        npower /= 2;
                    }
                } else {
                    let mut temploc = self.loc;
                    loop {
                        if is_odd(npower) {
                            self.loc.add(&temploc.multiplied_mat2d(&self.matrix));
                            self.matrix.multiply(&tempmatrix);
                        }
                        if npower == 1 { break; }
                        temploc = temploc.added(&temploc.multiplied_mat2d(&tempmatrix));
                        tempmatrix = tempmatrix.multiplied(&tempmatrix);
                        npower /= 2;
                    }
                }
            } else if matches!(self.shape, PntMirror | Ax1Mirror) {
                if n % 2 == 0 {
                    self.shape = Identity;
                    self.scale = 1.0;
                    self.matrix.set_identity();
                    self.loc = GpXY::zero();
                }
            } else {
                self.shape = CompoundTrsf;
                let mut npower = n.abs() - 1;
                let s = self.scale;
                self.matrix.set_diagonal(
                    s * self.matrix.value(0, 0),
                    s * self.matrix.value(1, 1),
                );
                let mut temploc = self.loc;
                let mut tempscale = self.scale;
                let mut tempmatrix = self.matrix;
                loop {
                    if is_odd(npower) {
                        self.loc.add(&temploc.multiplied_mat2d(&self.matrix).multiplied(self.scale));
                        self.scale = self.scale * tempscale;
                        self.matrix.multiply(&tempmatrix);
                    }
                    if npower == 1 { break; }
                    tempscale = tempscale * tempscale;
                    temploc = temploc.added(&temploc.multiplied_mat2d(&tempmatrix).multiplied(tempscale));
                    tempmatrix = tempmatrix.multiplied(&tempmatrix);
                    npower /= 2;
                }
            }
        }
        Ok(())
    }

    /// `gp_Trsf2d::PreMultiply(const gp_Trsf2d&)` (`gp_Trsf2d.cxx:550-674`).
    pub fn pre_multiply(&mut self, t: &Self) {
        use TrsfForm::*;
        if t.shape == Identity {
            return;
        } else if self.shape == Identity {
            self.shape = t.shape;
            self.scale = t.scale;
            self.loc = t.loc;
            self.matrix = t.matrix;
        } else if self.shape == Rotation && t.shape == Rotation {
            self.loc.multiply_mat2d(&t.matrix);
            self.loc.add(&t.loc);
            self.matrix.pre_multiply(&t.matrix);
        } else if self.shape == Translation && t.shape == Translation {
            self.loc.add(&t.loc);
        } else if self.shape == Scale && t.shape == Scale {
            self.loc.multiply_scalar(t.scale);
            self.loc.add(&t.loc);
            self.scale = self.scale * t.scale;
        } else if self.shape == PntMirror && t.shape == PntMirror {
            self.scale = 1.0;
            self.shape = Translation;
            self.loc.reverse();
            self.loc.add(&t.loc);
        } else if self.shape == Ax1Mirror && t.shape == Ax1Mirror {
            self.shape = Rotation;
            self.loc.multiply_mat2d(&t.matrix);
            self.loc.multiply_scalar(t.scale);
            self.scale = self.scale * t.scale;
            self.loc.add(&t.loc);
            self.matrix.pre_multiply(&t.matrix);
        } else if matches!(self.shape, CompoundTrsf | Rotation | Ax1Mirror | Ax2Mirror)
            && t.shape == Translation
        {
            self.loc.add(&t.loc);
        } else if matches!(self.shape, Scale | PntMirror) && t.shape == Translation {
            self.loc.add(&t.loc);
        } else if self.shape == Translation
            && matches!(t.shape, CompoundTrsf | Rotation | Ax1Mirror | Ax2Mirror)
        {
            self.shape = CompoundTrsf;
            self.matrix = t.matrix;
            if t.scale == 1.0 {
                self.loc.multiply_mat2d(&t.matrix);
            } else {
                self.scale = t.scale;
                self.loc.multiply_mat2d(&self.matrix);
                self.loc.multiply_scalar(self.scale);
            }
            self.loc.add(&t.loc);
        } else if matches!(t.shape, Scale | PntMirror) && self.shape == Translation {
            self.loc.multiply_scalar(t.scale);
            self.loc.add(&t.loc);
            self.scale = t.scale;
            self.shape = t.shape;
        } else if matches!(self.shape, PntMirror | Scale) && matches!(t.shape, PntMirror | Scale) {
            self.shape = CompoundTrsf;
            self.loc.multiply_scalar(t.scale);
            self.loc.add(&t.loc);
            self.scale = self.scale * t.scale;
        } else if matches!(self.shape, CompoundTrsf | Rotation | Ax1Mirror | Ax2Mirror)
            && matches!(t.shape, Scale | PntMirror)
        {
            self.shape = CompoundTrsf;
            self.loc.multiply_scalar(t.scale);
            self.loc.add(&t.loc);
            self.scale = self.scale * t.scale;
        } else if matches!(t.shape, CompoundTrsf | Rotation | Ax1Mirror | Ax2Mirror)
            && matches!(self.shape, Scale | PntMirror)
        {
            self.shape = CompoundTrsf;
            self.matrix = t.matrix;
            if t.scale == 1.0 {
                self.loc.multiply_mat2d(&t.matrix);
            } else {
                self.loc.multiply_mat2d(&self.matrix);
                self.loc.multiply_scalar(t.scale);
                self.scale = t.scale * self.scale;
            }
            self.loc.add(&t.loc);
        } else {
            self.shape = CompoundTrsf;
            self.loc.multiply_mat2d(&t.matrix);
            if t.scale != 1.0 {
                self.loc.multiply_scalar(t.scale);
                self.scale = self.scale * t.scale;
            }
            self.loc.add(&t.loc);
            self.matrix.pre_multiply(&t.matrix);
        }
    }

    /// `gp_Trsf2d::SetValues(...)` (`gp_Trsf2d.cxx:676-710`): 6 coefficients
    /// (two columns then translation). `scale` is the square root of the
    /// determinant, then `Orthogonalize`.
    pub fn set_values(
        &mut self, a11: f64, a12: f64, a13: f64,
        a21: f64, a22: f64, a23: f64,
    ) -> Result<(),&'static str> {
        let col1 = GpXY::new(a11, a21);
        let col2 = GpXY::new(a12, a22);
        let col3 = GpXY::new(a13, a23);
        let mut m = GpMat2d::from_cols(&col1, &col2);
        let det = m.determinant();
        if det.abs() < RESOLUTION {
            return Err("gp_Trsf2d::SetValues, null determinant");
        }
        let s = if det > 0.0 { det.sqrt() } else { (-det).sqrt() };
        m.divide(s);

        self.scale = s;
        self.shape = TrsfForm::CompoundTrsf;
        self.matrix = m;
        self.orthogonalize();
        self.loc = col3;
        Ok(())
    }

    /// `gp_Trsf2d::Orthogonalize()` (`gp_Trsf2d.cxx:721-746`): Gram-Schmidt on
    /// the columns then on the rows.
    pub fn orthogonalize(&mut self) {
        let mut tm = self.matrix;
        let mut v1 = tm.column(0);
        let mut v2 = tm.column(1);

        let _ = v1.normalize();
        v2 = v2.subtracted(&v1.multiplied(v2.dot(&v1)));
        let _ = v2.normalize();
        tm.set_cols(&v1, &v2);

        let mut r1 = tm.row(0);
        let mut r2 = tm.row(1);

        let _ = r1.normalize();
        r2 = r2.subtracted(&r1.multiplied(r2.dot(&r1)));
        let _ = r2.normalize();
        tm.set_rows(&r1, &r2);

        self.matrix = tm;
    }

    pub fn transforms_xy(&self,coord:&mut GpXY) { coord.multiply_mat2d(&self.matrix); if self.scale!=1.0 { coord.multiply_scalar(self.scale); } coord.add(&self.loc); }
    pub fn multiply(&mut self,o:&Self) { let mut nl=o.loc; nl.multiply_mat2d(&self.matrix); if self.scale!=1.0 { nl.multiply_scalar(self.scale); } nl.add(&self.loc); self.loc=nl; self.matrix.multiply(&o.matrix); self.scale*=o.scale; if self.shape==TrsfForm::Identity { self.shape=o.shape; } else if o.shape!=TrsfForm::Identity { self.shape=TrsfForm::CompoundTrsf; } }
    pub fn multiplied(&self,o:&Self) -> Self { let mut r=*self; r.multiply(o); r }
    pub fn invert(&mut self) -> Result<(),&'static str> { self.matrix.invert()?; self.loc.reverse(); self.loc.multiply_mat2d(&self.matrix); if self.scale!=1.0 { self.loc.divide(self.scale); } self.scale=1.0/self.scale; Ok(()) }
    pub fn inverted(&self) -> Result<Self,&'static str> { let mut r=*self; r.invert()?; Ok(r) }
}
impl Default for GpTrsf2d { fn default() -> Self { Self::identity() } }
