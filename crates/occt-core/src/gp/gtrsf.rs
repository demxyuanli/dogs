//! General (non-uniform) affine 3D transform. Source: `gp_GTrsf.hxx`/`gp_GTrsf.cxx`.
//! Unlike `GpTrsf`, a `gp_GTrsf` may be non-orthogonal (an affinity) and can
//! change the nature of a geometric object; use `GpTrsf` when the nature must
//! be preserved. When `shape == Other`, `scale` is a `0.0` placeholder and the
//! stored `matrix` is the final vectorial part (see `Transforms`).
use crate::gp::ax1::GpAx1;
use crate::gp::ax2::GpAx2;
use crate::gp::mat::GpMat;
use crate::gp::trsf::GpTrsf;
use crate::gp::trsf_form::TrsfForm;
use crate::gp::xyz::GpXyz;
use crate::precision::RESOLUTION;

/// 3x4 affine matrix. Source: `gp_GTrsf.hxx:57`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpGTrsf {
    pub matrix: GpMat,
    pub loc: GpXyz,
    pub scale: f64,
    pub shape: TrsfForm,
}

impl GpGTrsf {
    /// Identity transformation. Source: `gp_GTrsf.hxx:63-70`
    #[inline]
    pub fn identity() -> Self {
        Self { matrix: GpMat::identity(), loc: GpXyz::zero(), scale: 1.0, shape: TrsfForm::Identity }
    }

    /// `gp_GTrsf(const gp_Trsf&)` (`gp_GTrsf.hxx:73-81`).
    pub fn from_trsf(t: &GpTrsf) -> Self {
        Self { matrix: t.matrix, loc: t.loc, scale: t.scale, shape: t.shape }
    }

    /// `gp_GTrsf(const gp_Mat&, const gp_XYZ&)` (`gp_GTrsf.hxx:85-91`).
    pub fn from_matrix_vector(m: &GpMat, v: &GpXyz) -> Self {
        Self { matrix: *m, loc: *v, scale: 0.0, shape: TrsfForm::Other }
    }

    /// `gp_GTrsf::SetAffinity(const gp_Ax1&, double)` (`gp_GTrsf.hxx:316-328`).
    pub fn set_affinity_ax1(&mut self, a1: &GpAx1, ratio: f64) {
        self.shape = TrsfForm::Other;
        self.scale = 0.0;
        self.matrix.set_dot(a1.direction().xyz());
        self.matrix = self.matrix.multiply_scalar(1.0 - ratio);
        let (v11, v22, v33) =
            (self.matrix.value(1, 1), self.matrix.value(2, 2), self.matrix.value(3, 3));
        self.matrix.set_diagonal(&GpXyz::new(v11 + ratio, v22 + ratio, v33 + ratio));
        let o = *a1.location().xyz();
        self.loc = o;
        self.loc.reverse();
        self.loc.multiply_mat(&self.matrix);
        self.loc = self.loc.add(&o);
    }

    /// `gp_GTrsf::SetAffinity(const gp_Ax2&, double)` (`gp_GTrsf.hxx:333-343`).
    /// The translation is computed with the matrix *before* the diagonal `+1`.
    pub fn set_affinity_ax2(&mut self, a2: &GpAx2, ratio: f64) {
        self.shape = TrsfForm::Other;
        self.scale = 0.0;
        let d = a2.direction();
        self.matrix.set_dot(d.xyz());
        self.matrix = self.matrix.multiply_scalar(ratio - 1.0);
        let o = *a2.location().xyz();
        self.loc = o;
        self.loc.reverse();
        self.loc.multiply_mat(&self.matrix);
        let (v11, v22, v33) =
            (self.matrix.value(1, 1), self.matrix.value(2, 2), self.matrix.value(3, 3));
        self.matrix.set_diagonal(&GpXyz::new(v11 + 1.0, v22 + 1.0, v33 + 1.0));
    }

    /// `gp_GTrsf::SetValue(int, int, double)` (`gp_GTrsf.hxx:347-370`).
    pub fn set_value(&mut self, row: usize, col: usize, value: f64) {
        if col == 4 {
            match row {
                1 => self.loc.set_x(value),
                2 => self.loc.set_y(value),
                _ => self.loc.set_z(value),
            }
            if self.shape == TrsfForm::Identity {
                self.shape = TrsfForm::Translation;
            }
            return;
        }
        if self.shape != TrsfForm::Other && self.scale != 1.0 {
            self.matrix = self.matrix.multiply_scalar(self.scale);
        }
        self.matrix.set_value(row, col, value);
        self.shape = TrsfForm::Other;
        self.scale = 0.0;
    }

    /// `gp_GTrsf::SetVectorialPart` (`gp_GTrsf.hxx:118-123`).
    pub fn set_vectorial_part(&mut self, m: &GpMat) {
        self.matrix = *m;
        self.shape = TrsfForm::Other;
        self.scale = 0.0;
    }

    /// `gp_GTrsf::SetTranslationPart` (`gp_GTrsf.cxx:28-43`).
    pub fn set_translation_part(&mut self, coord: &GpXyz) {
        self.loc = *coord;
        match self.shape {
            TrsfForm::CompoundTrsf | TrsfForm::Other | TrsfForm::Translation => {}
            TrsfForm::Identity => self.shape = TrsfForm::Translation,
            _ => self.shape = TrsfForm::CompoundTrsf,
        }
    }

    /// `gp_GTrsf::SetTrsf` (`gp_GTrsf.hxx:130-136`).
    pub fn set_trsf(&mut self, t: &GpTrsf) {
        self.shape = t.shape;
        self.matrix = t.matrix;
        self.loc = t.loc;
        self.scale = t.scale;
    }

    /// `gp_GTrsf::IsNegative` (`gp_GTrsf.hxx:140`).
    #[inline]
    pub fn is_negative(&self) -> bool {
        self.matrix.determinant() < 0.0
    }

    /// `gp_GTrsf::IsSingular` (`gp_GTrsf.hxx:150`).
    #[inline]
    pub fn is_singular(&self) -> bool {
        self.matrix.is_singular()
    }

    /// `gp_GTrsf::Form` (`gp_GTrsf.hxx:157`).
    #[inline]
    pub fn form(&self) -> TrsfForm {
        self.shape
    }

    /// `gp_GTrsf::SetForm` (`gp_GTrsf.cxx:150-197`): recompute `shape` from the
    /// matrix (uniform => `CompoundTrsf`, otherwise `Other`).
    pub fn set_form(&mut self) -> Result<(), &'static str> {
        const TOL: f64 = 1.0e-12;
        let mut m = self.matrix;
        let mut s = m.determinant();
        if s.abs() < RESOLUTION {
            return Err("gp_GTrsf::SetForm, null determinant");
        }
        s = if s > 0.0 { s.powf(1.0 / 3.0) } else { -(-s).powf(1.0 / 3.0) };
        m = m.divide(s)?;

        let mut tm = m.transpose().multiply(&m);
        tm = tm.subtract(&GpMat::identity());
        if self.shape == TrsfForm::Other {
            self.shape = TrsfForm::CompoundTrsf;
        }
        for i in 1..=3 {
            for j in 1..=3 {
                if tm.value(i, j).abs() > TOL {
                    self.shape = TrsfForm::Other;
                    return Ok(());
                }
            }
        }
        Ok(())
    }

    /// `gp_GTrsf::TranslationPart` (`gp_GTrsf.hxx:170`).
    #[inline]
    pub fn translation_part(&self) -> GpXyz {
        self.loc
    }

    /// `gp_GTrsf::VectorialPart` (`gp_GTrsf.hxx:174`).
    #[inline]
    pub fn vectorial_part(&self) -> GpMat {
        self.matrix
    }

    /// `gp_GTrsf::Value(int, int)` (`gp_GTrsf.hxx:374-386`), 1-based indices.
    pub fn value(&self, row: usize, col: usize) -> f64 {
        if col == 4 {
            return self.loc.coord(row - 1);
        }
        if self.shape == TrsfForm::Other {
            self.matrix.value(row, col)
        } else {
            self.scale * self.matrix.value(row, col)
        }
    }

    /// `gp_GTrsf::Invert` (`gp_GTrsf.cxx:44-58`).
    pub fn invert(&mut self) -> Result<(), &'static str> {
        if self.shape == TrsfForm::Other {
            self.matrix.invert()?;
            self.loc.multiply_mat(&self.matrix);
            self.loc.reverse();
        } else {
            let mut t = self.to_trsf();
            t.invert()?;
            self.set_trsf(&t);
        }
        Ok(())
    }

    /// `gp_GTrsf::Inverted` (`gp_GTrsf.hxx:187-192`).
    pub fn inverted(&self) -> Result<Self, &'static str> {
        let mut r = *self;
        r.invert()?;
        Ok(r)
    }

    /// `gp_GTrsf::Multiply` (`gp_GTrsf.cxx:60-78`): `self = self * t`.
    pub fn multiply(&mut self, t: &Self) {
        if self.shape == TrsfForm::Other || t.shape == TrsfForm::Other {
            self.shape = TrsfForm::Other;
            // Translation uses the current (pre-multiply) matrix.
            self.loc = self.loc.add(&t.loc.multiplied_mat(&self.matrix));
            self.matrix = self.matrix.multiply(&t.matrix);
        } else {
            let mut t1 = self.to_trsf();
            let t2 = t.to_trsf();
            t1.multiply(&t2);
            self.matrix = t1.matrix;
            self.loc = t1.loc;
            self.scale = t1.scale;
            self.shape = t1.shape;
        }
    }

    /// `gp_GTrsf::Multiplied` (`gp_GTrsf.hxx:209-213`).
    pub fn multiplied(&self, t: &Self) -> Self {
        let mut r = *self;
        r.multiply(t);
        r
    }

    /// `gp_GTrsf::PreMultiply` (`gp_GTrsf.cxx:129-148`): `self = t * self`.
    pub fn pre_multiply(&mut self, t: &Self) {
        if self.shape == TrsfForm::Other || t.shape == TrsfForm::Other {
            self.shape = TrsfForm::Other;
            self.loc.multiply_mat(&t.matrix);
            self.loc = self.loc.add(&t.loc);
            self.matrix = self.matrix.pre_multiply(&t.matrix);
        } else {
            let mut t1 = self.to_trsf();
            let t2 = t.to_trsf();
            t1.pre_multiply(&t2);
            self.matrix = t1.matrix;
            self.loc = t1.loc;
            self.scale = t1.scale;
            self.shape = t1.shape;
        }
    }

    /// `gp_GTrsf::Power` (`gp_GTrsf.cxx:80-127`).
    pub fn power(&mut self, n: i32) -> Result<(), &'static str> {
        if n == 0 {
            self.scale = 1.0;
            self.shape = TrsfForm::Identity;
            self.matrix = GpMat::identity();
            self.loc = GpXyz::zero();
        } else if n == 1 {
            // no-op
        } else if n == -1 {
            self.invert()?;
        } else if self.shape == TrsfForm::Other {
            let mut npower = n.abs() - 1;
            let mut temploc = self.loc;
            let mut tempmatrix = self.matrix;
            loop {
                if npower % 2 == 1 {
                    self.loc = self.loc.add(&temploc.multiplied_mat(&self.matrix));
                    self.matrix = self.matrix.multiply(&tempmatrix);
                }
                if npower == 1 {
                    break;
                }
                temploc = temploc.add(&temploc.multiplied_mat(&tempmatrix));
                tempmatrix = tempmatrix.multiply(&tempmatrix);
                npower /= 2;
            }
        } else {
            let mut t = self.to_trsf();
            t.power(n)?;
            self.set_trsf(&t);
        }
        Ok(())
    }

    /// `gp_GTrsf::Powered` (`gp_GTrsf.hxx:244-248`).
    pub fn powered(&self, n: i32) -> Result<Self, &'static str> {
        let mut r = *self;
        r.power(n)?;
        Ok(r)
    }

    /// `gp_GTrsf::Transforms(gp_XYZ&)` (`gp_GTrsf.hxx:390-396`). Note the
    /// factor `scale` is applied only when `shape != Other` (`Other` stores the
    /// final matrix with `scale == 0.0`).
    pub fn transforms(&self, xyz: &mut GpXyz) {
        xyz.multiply_mat(&self.matrix);
        if self.shape != TrsfForm::Other && self.scale != 1.0 {
            *xyz = xyz.multiplied(self.scale);
        }
        *xyz = xyz.add(&self.loc);
    }

    /// Transformed copy.
    pub fn transformed(&self, xyz: &GpXyz) -> GpXyz {
        let mut r = *xyz;
        self.transforms(&mut r);
        r
    }

    /// `gp_GTrsf::GetMat4(NCollection_Mat4<T>&)` (`gp_GTrsf.hxx:260-286`):
    /// row-major 4x4 with rows 1-3 from `Value(row, col)` and last row
    /// `[0, 0, 0, 1]`. `NCollection_Mat4` has no port equivalent, so the 4x4 is
    /// mirrored as a plain `[[f64; 4]; 4]` with `SetValue(row, col)` semantics.
    pub fn get_mat4(&self) -> [[f64; 4]; 4] {
        if self.shape == TrsfForm::Identity {
            return [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ];
        }
        let mut m = [[0.0f64; 4]; 4];
        for r in 1..=3usize {
            for c in 1..=4usize {
                m[r - 1][c - 1] = self.value(r, c);
            }
        }
        m[3][3] = 1.0;
        m
    }

    /// `gp_GTrsf::SetMat4(const NCollection_Mat4<T>&)` (`gp_GTrsf.hxx:288-303`):
    /// takes only the upper-left 3x3 and the fourth column, and sets
    /// `shape = Other`, `scale = 0.0`.
    pub fn set_mat4(&mut self, m: &[[f64; 4]; 4]) {
        self.shape = TrsfForm::Other;
        self.scale = 0.0;
        for r in 1..=3usize {
            for c in 1..=3usize {
                self.matrix.set_value(r, c, m[r - 1][c - 1]);
            }
        }
        self.loc = GpXyz::new(m[0][3], m[1][3], m[2][3]);
    }

    /// UNPORTED: `gp_GTrsf::DumpJson` (`gp_GTrsf.cxx:202-...`) needs the OCCT
    /// `Standard_Dump` / `Standard_OStream` serialization framework, absent from
    /// this port. No equivalent branch, so not added.
    /// `gp_GTrsf::Trsf` (`gp_GTrsf.hxx:416-428`); errors on a non-orthogonal
    /// transformation (`Form() == Other`).
    pub fn trsf(&self) -> Result<GpTrsf, &'static str> {
        if self.shape == TrsfForm::Other {
            return Err("gp_GTrsf::Trsf() - non-orthogonal GTrsf");
        }
        Ok(self.to_trsf())
    }

    /// Build a `GpTrsf` from the stored parts without the `Other` check
    /// (internal helper for the `.cxx` `Trsf()` calls that already exclude it).
    fn to_trsf(&self) -> GpTrsf {
        GpTrsf { scale: self.scale, shape: self.shape, matrix: self.matrix, loc: self.loc }
    }
}

impl Default for GpGTrsf {
    fn default() -> Self {
        Self::identity()
    }
}
