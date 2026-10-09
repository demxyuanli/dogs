//! 2D general (affine) transform. Source: `gp_GTrsf2d.hxx`/`gp_GTrsf2d.cxx`.
//!
//! Unlike `gp_Trsf2d` this may be non-orthogonal (an affinity); it is the type
//! `ShapeBuild_Edge::TransformPCurve` uses to scale the U parameter of a pcurve
//! (`ShapeBuild_Edge.cxx:620-698`). When `shape == Other`, `scale` is a `0.0`
//! placeholder and the stored `matrix` is the final vectorial part.
use crate::gp::{GpAx2d, GpMat2d, GpPnt2d, GpTrsf2d, GpXY, TrsfForm};
use crate::precision::{ANGULAR, RESOLUTION};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpGTrsf2d {
    pub matrix: GpMat2d,
    pub loc: GpXY,
    pub scale: f64,
    pub shape: TrsfForm,
}

impl GpGTrsf2d {
    /// `gp_GTrsf2d::gp_GTrsf2d()` (`gp_GTrsf2d.hxx:54-59`).
    pub fn identity() -> Self {
        Self {
            matrix: GpMat2d::identity(),
            loc: GpXY::zero(),
            scale: 1.0,
            shape: TrsfForm::Identity,
        }
    }

    /// `gp_GTrsf2d(const gp_Trsf2d&)` (`gp_GTrsf2d.hxx:244-251`).
    pub fn from_trsf2d(t: &GpTrsf2d) -> Self {
        Self { matrix: t.matrix, loc: t.loc, scale: t.scale, shape: t.shape }
    }

    /// `gp_GTrsf2d::SetAffinity(A, Ratio)` (`gp_GTrsf2d.cxx:24-38`).
    pub fn set_affinity(&mut self, a: &GpAx2d, ratio: f64) {
        self.shape = TrsfForm::Other;
        self.scale = 0.0;
        let (ax, ay) = (a.direction().x, a.direction().y);
        self.matrix.data[0][0] = (1.0 - ratio) * ax * ax + ratio;
        self.matrix.data[1][1] = (1.0 - ratio) * ay * ay + ratio;
        self.matrix.data[0][1] = (1.0 - ratio) * ax * ay;
        self.matrix.data[1][0] = self.matrix.data[0][1];
        let base = *a.location().xy();
        let mut l = base;
        l.reverse();
        l.multiply_mat2d(&self.matrix);
        l.add(&base);
        self.loc = l;
    }

    /// `gp_GTrsf2d::SetTranslationPart` (`gp_GTrsf2d.cxx:40-54`).
    pub fn set_translation_part(&mut self, coord: &GpXY) {
        self.loc = *coord;
        match self.shape {
            TrsfForm::CompoundTrsf | TrsfForm::Other | TrsfForm::Translation => {}
            TrsfForm::Identity => self.shape = TrsfForm::Translation,
            _ => self.shape = TrsfForm::CompoundTrsf,
        }
    }

    /// `gp_GTrsf2d::SetTrsf2d` (`gp_GTrsf2d.hxx:235-241`).
    pub fn set_trsf2d(&mut self, t: &GpTrsf2d) {
        self.shape = t.shape;
        self.matrix = t.matrix;
        self.loc = t.loc;
        self.scale = t.scale;
    }

    /// `gp_GTrsf2d::SetValue(int, int, double)` (`gp_GTrsf2d.hxx:255-267`).
    pub fn set_value(&mut self, row: usize, col: usize, value: f64) {
        if col == 3 {
            match row {
                1 => self.loc.x = value,
                _ => self.loc.y = value,
            }
        } else {
            self.matrix.data[row - 1][col - 1] = value;
        }
        self.shape = TrsfForm::Other;
    }

    /// `gp_GTrsf2d::SetVectorialPart` (`gp_GTrsf2d.hxx:98-103`).
    pub fn set_vectorial_part(&mut self, m: &GpMat2d) {
        self.matrix = *m;
        self.shape = TrsfForm::Other;
        self.scale = 0.0;
    }

    /// `gp_GTrsf2d::IsNegative` (`gp_GTrsf2d.hxx:107`).
    #[inline]
    pub fn is_negative(&self) -> bool {
        self.matrix.determinant() < 0.0
    }

    /// `gp_GTrsf2d::IsSingular` (`gp_GTrsf2d.hxx:117`).
    #[inline]
    pub fn is_singular(&self) -> bool {
        self.matrix.determinant().abs() <= RESOLUTION
    }

    /// `gp_GTrsf2d::Form` (`gp_GTrsf2d.hxx:124`).
    #[inline]
    pub fn form(&self) -> TrsfForm {
        self.shape
    }

    /// `gp_GTrsf2d::TranslationPart` (`gp_GTrsf2d.hxx:131` translation part).
    #[inline]
    pub fn translation_part(&self) -> GpXY {
        self.loc
    }

    /// `gp_GTrsf2d::VectorialPart` (`gp_GTrsf2d.hxx:131`).
    #[inline]
    pub fn vectorial_part(&self) -> GpMat2d {
        self.matrix
    }

    /// `gp_GTrsf2d::Value(int, int)` (`gp_GTrsf2d.hxx:271-283`), 1-based.
    pub fn value(&self, row: usize, col: usize) -> f64 {
        if col == 3 {
            return if row == 1 { self.loc.x } else { self.loc.y };
        }
        if self.shape == TrsfForm::Other {
            self.matrix.data[row - 1][col - 1]
        } else {
            self.scale * self.matrix.data[row - 1][col - 1]
        }
    }

    /// `gp_GTrsf2d::Invert` (`gp_GTrsf2d.cxx:56-70`).
    pub fn invert(&mut self) -> Result<(), &'static str> {
        if self.shape == TrsfForm::Other {
            if self.is_singular() {
                return Err("gp_GTrsf2d::Invert() - singular matrix");
            }
            self.matrix.invert()?;
            self.loc.multiply_mat2d(&self.matrix);
            self.loc.reverse();
        } else {
            let mut t = self.to_trsf2d();
            t.invert()?;
            self.set_trsf2d(&t);
        }
        Ok(())
    }

    /// `gp_GTrsf2d::Inverted` (`gp_GTrsf2d.hxx:144-148`).
    pub fn inverted(&self) -> Result<Self, &'static str> {
        let mut r = *self;
        r.invert()?;
        Ok(r)
    }

    /// `gp_GTrsf2d::Multiply` (`gp_GTrsf2d.cxx:72-90`): `self = self * t`.
    pub fn multiply(&mut self, t: &Self) {
        if self.shape == TrsfForm::Other || t.shape == TrsfForm::Other {
            self.shape = TrsfForm::Other;
            self.loc = self.loc.added(&t.loc.multiplied_mat2d(&self.matrix));
            self.matrix.multiply(&t.matrix);
        } else {
            let mut t1 = self.to_trsf2d();
            let t2 = t.to_trsf2d();
            t1.multiply(&t2);
            self.matrix = t1.matrix;
            self.loc = t1.loc;
            self.scale = t1.scale;
            self.shape = t1.shape;
        }
    }

    /// `gp_GTrsf2d::Multiplied` (`gp_GTrsf2d.hxx:166-170`).
    pub fn multiplied(&self, t: &Self) -> Self {
        let mut r = *self;
        r.multiply(t);
        r
    }

    /// `gp_GTrsf2d::Power` (`gp_GTrsf2d.cxx:92-143`).
    pub fn power(&mut self, n: i32) -> Result<(), &'static str> {
        if n == 0 {
            self.scale = 1.0;
            self.shape = TrsfForm::Identity;
            self.matrix.set_identity();
            self.loc = GpXY::zero();
        } else if n == 1 {
            // no-op
        } else if n == -1 {
            self.invert()?;
        } else {
            if n < 0 {
                self.invert()?;
            }
            if self.shape == TrsfForm::Other {
                let mut npower = n.abs() - 1;
                let mut temploc = self.loc;
                let mut tempmatrix = self.matrix;
                loop {
                    if npower % 2 == 1 {
                        self.loc = self.loc.added(&temploc.multiplied_mat2d(&self.matrix));
                        self.matrix.multiply(&tempmatrix);
                    }
                    if npower == 1 {
                        break;
                    }
                    temploc = temploc.added(&temploc.multiplied_mat2d(&tempmatrix));
                    tempmatrix = tempmatrix.multiplied(&tempmatrix);
                    npower /= 2;
                }
            } else {
                let mut t = self.to_trsf2d();
                t.power(n)?;
                self.set_trsf2d(&t);
            }
        }
        Ok(())
    }

    /// `gp_GTrsf2d::Powered` (`gp_GTrsf2d.hxx:193-197`).
    pub fn powered(&self, n: i32) -> Result<Self, &'static str> {
        let mut r = *self;
        r.power(n)?;
        Ok(r)
    }

    /// `gp_GTrsf2d::PreMultiply` (`gp_GTrsf2d.cxx:145-163`): `self = t * self`.
    pub fn pre_multiply(&mut self, t: &Self) {
        if self.shape == TrsfForm::Other || t.shape == TrsfForm::Other {
            self.shape = TrsfForm::Other;
            self.loc.multiply_mat2d(&t.matrix);
            self.loc = self.loc.added(&t.loc);
            self.matrix.pre_multiply(&t.matrix);
        } else {
            let mut t1 = self.to_trsf2d();
            let t2 = t.to_trsf2d();
            t1.pre_multiply(&t2);
            self.matrix = t1.matrix;
            self.loc = t1.loc;
            self.scale = t1.scale;
            self.shape = t1.shape;
        }
    }

    /// `gp_GTrsf2d::Transforms(gp_XY&)` (`gp_GTrsf2d.hxx:287-295`).
    pub fn transforms_xy(&self, coord: &mut GpXY) {
        coord.multiply_mat2d(&self.matrix);
        if self.shape != TrsfForm::Other && self.scale != 1.0 {
            coord.multiply_scalar(self.scale);
        }
        coord.add(&self.loc);
    }

    /// Transformed copy (`gp_GTrsf2d::Transformed`, `gp_GTrsf2d.hxx:200-206`).
    pub fn transformed_xy(&self, coord: &GpXY) -> GpXY {
        let mut c = *coord;
        self.transforms_xy(&mut c);
        c
    }

    /// Convenience point wrapper around `transforms_xy` (existing caller API).
    pub fn transforms(&self, p: &GpPnt2d) -> GpPnt2d {
        let mut xy = *p.xy();
        self.transforms_xy(&mut xy);
        GpPnt2d::from_xy(xy)
    }

    /// `gp_GTrsf2d::Trsf2d` (`gp_GTrsf2d.cxx:165-204`); errors on a
    /// non-orthogonal transformation (`Form() == Other`).
    pub fn trsf2d(&self) -> Result<GpTrsf2d, &'static str> {
        let tol = ANGULAR;
        let tol2 = 2.0 * tol;

        if self.shape == TrsfForm::Other {
            return Err("gp_GTrsf2d::Trsf2d() - non-orthogonal GTrsf2d(0)");
        }
        let m = |r: usize, c: usize| self.matrix.data[r - 1][c - 1];

        let mut value = m(1, 1) * m(1, 1) + m(2, 1) * m(2, 1);
        if (value - 1.0).abs() > tol2 {
            return Err("gp_GTrsf2d::Trsf2d() - non-orthogonal GTrsf2d(1)");
        }
        value = m(1, 2) * m(1, 2) + m(2, 2) * m(2, 2);
        if (value - 1.0).abs() > tol2 {
            return Err("gp_GTrsf2d::Trsf2d() - non-orthogonal GTrsf2d(2)");
        }
        value = m(1, 1) * m(1, 2) + m(2, 1) * m(2, 2);
        if value.abs() > tol {
            return Err("gp_GTrsf2d::Trsf2d() - non-orthogonal GTrsf2d(3)");
        }
        Ok(self.to_trsf2d())
    }

    /// Build a `GpTrsf2d` from the stored parts without the `Other` check
    /// (internal helper for the `.cxx` `Trsf2d()` calls that already exclude it).
    fn to_trsf2d(&self) -> GpTrsf2d {
        GpTrsf2d { scale: self.scale, shape: self.shape, matrix: self.matrix, loc: self.loc }
    }
}

impl Default for GpGTrsf2d {
    fn default() -> Self {
        Self::identity()
    }
}
