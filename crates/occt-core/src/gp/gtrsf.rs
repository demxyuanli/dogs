//! General (non-uniform) affine 3D transform. Source: `gp_GTrsf.hxx`
//! Unlike GpTrsf, this can represent non-uniform scaling (affinity).
//! Only applicable to coordinates, not geometric objects (may change nature).
use crate::gp::{mat::GpMat, xyz::GpXyz, trsf::GpTrsf, trsf_form::TrsfForm};

/// 3×4 affine matrix: `x' = scale * matrix * x + loc`.
/// Supports non-uniform scaling. Source: `gp_GTrsf.hxx:57`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpGTrsf {
    pub matrix: GpMat,
    pub loc: GpXyz,
    pub scale: f64,
    pub shape: TrsfForm,
}

impl GpGTrsf {
    /// Identity transform. Source: `gp_GTrsf.hxx:63`
    #[inline]
    pub fn identity() -> Self {
        Self { matrix: GpMat::identity(), loc: GpXyz::zero(), scale: 1.0, shape: TrsfForm::Identity }
    }

    /// From uniform transform. Source: `gp_GTrsf.hxx:74`
    pub fn from_trsf(t: &GpTrsf) -> Self {
        Self {
            matrix: t.matrix,
            loc: t.translation_part(),
            scale: t.scale_factor(),
            shape: t.form(),
        }
    }

    /// Set vectorial part. Source: `gp_GTrsf.hxx` (SetVectorialPart)
    pub fn set_vectorial_part(&mut self, m: &GpMat) { self.matrix = *m; }

    /// Set translation part. Source: `gp_GTrsf.hxx` (SetTranslationPart)
    pub fn set_translation_part(&mut self, v: &GpXyz) { self.loc = *v; }

    /// Set form. Source: `gp_GTrsf.hxx` (SetForm)
    pub fn set_form(&mut self, f: TrsfForm) { self.shape = f; }

    /// Apply transform to a point. Source: `gp_GTrsf.hxx` (Transforms)
    pub fn transforms(&self, xyz: &mut GpXyz) {
        let mut r = xyz.multiplied_mat(&self.matrix);
        if self.scale != 1.0 { r = r.multiplied(self.scale); }
        *xyz = r.added(&self.loc);
    }

    /// Transformed copy.
    pub fn transformed(&self, xyz: &GpXyz) -> GpXyz {
        let mut r = *xyz;
        self.transforms(&mut r);
        r
    }

    /// Compose: `self = self * other`. Source: `gp_GTrsf.hxx` (Multiply)
    pub fn multiply(&mut self, other: &Self) {
        let mut new_loc = other.loc;
        new_loc.multiply_mat(&self.matrix);
        if self.scale != 1.0 { new_loc = new_loc.multiplied(self.scale); }
        self.loc = new_loc.added(&self.loc);
        self.matrix = self.matrix.multiply(&other.matrix);
        self.scale *= other.scale;
        if self.shape != TrsfForm::Identity && other.shape != TrsfForm::Identity {
            self.shape = TrsfForm::CompoundTrsf;
        }
    }

    /// `self * other`.
    pub fn multiplied(&self, other: &Self) -> Self { let mut r = *self; r.multiply(other); r }

    /// Invert. Source: `gp_GTrsf.hxx` (Invert)
    pub fn invert(&mut self) -> Result<(), &'static str> {
        if self.scale.abs() < 1e-30 { return Err("singular scale"); }
        let inv = self.matrix.invert();
        let inv_scale = 1.0 / self.scale;
        let mut new_loc = self.loc.reversed();
        new_loc.multiply_mat(&inv);
        new_loc = new_loc.multiplied(inv_scale);
        self.matrix = inv;
        self.loc = new_loc;
        self.scale = inv_scale;
        Ok(())
    }

    /// Inverted copy.
    pub fn inverted(&self) -> Result<Self, &'static str> { let mut r = *self; r.invert()?; Ok(r) }

    /// Pre-multiply: `self = other * self`. Source: `gp_GTrsf.hxx` (PreMultiply)
    pub fn pre_multiply(&mut self, other: &Self) {
        let mut m = other.clone();
        m.multiply(self);
        *self = m;
    }

    /// Set uniform affine transform. Source: `gp_GTrsf.hxx` (SetAffinity)
    pub fn set_affinity(&mut self, a11:f64,a12:f64,a13:f64,a14:f64,
                        a21:f64,a22:f64,a23:f64,a24:f64,
                        a31:f64,a32:f64,a33:f64,a34:f64) {
        self.matrix = GpMat::new(a11,a12,a13, a21,a22,a23, a31,a32,a33);
        self.loc = GpXyz::new(a14, a24, a34);
        self.scale = 1.0;
        self.shape = TrsfForm::Other;
    }
}

impl Default for GpGTrsf { fn default() -> Self { Self::identity() } }
