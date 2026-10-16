//! `Extrema_GlobOptFuncCCC0/1/2` — the two-variable curve/curve distance
//! function that `Extrema_GGenExtCC::Perform` hands to `math_GlobOptMin`.
//!
//! Source: `Extrema_GlobOptFuncCC.cxx` (statics `_NbVariables`/`_Value`/
//! `_Gradient`/`_Hessian` at `:24-185`, the three class shells at `:191-388`) and
//! `Extrema_GlobOptFuncCC.hxx:25-99`.
//!
//! OCCT declares three classes only so that `math_GlobOptMin::computeLocalExtremum`
//! can `dynamic_cast` its way to the best local engine (Hessian → Newton, gradient
//! → BFGS, value → Powell; `math_GlobOptMin.cxx:266-339`). The maths is the same
//! in all three, so one struct implements all three ported objective traits and
//! [`MathMultipleVarFunction::as_hessian`] / `as_gradient` answer for it, exactly
//! as the two `dynamic_cast`s succeed on `Extrema_GlobOptFuncCCC2`.
//!
//! **OCCT quirk kept verbatim**: `_Value` returns the *distance*
//! (`C2.Value(v).Distance(C1.Value(u))`, `cxx:44`) while `_Gradient`/`_Hessian`
//! are the first/second derivatives of the *squared* distance (`cxx:87-91`,
//! `:140-151`, both scaled by 2). The minimizer is the same (distance is monotone
//! in squared distance), the scaling only affects step sizes.

use occt_math::{MathMatrix, MathMultipleVarFunction, MathMultipleVarFunctionWithGradient, MathMultipleVarFunctionWithHessian, MathVector};

use super::prelude::*;

/// `Extrema_GlobOptFuncCCC2(C1, C2)` (`Extrema_GlobOptFuncCC.cxx:306-314`).
pub struct GlobOptFuncCCC2<'a> {
    c1: &'a dyn Curve,
    c2: &'a dyn Curve,
}

impl<'a> GlobOptFuncCCC2<'a> {
    /// `Extrema_GlobOptFuncCCC2(const Adaptor3d_Curve&, const Adaptor3d_Curve&)`.
    pub fn new(c1: &'a dyn Curve, c2: &'a dyn Curve) -> Self {
        Self { c1, c2 }
    }

    /// `_NbVariables()` (`cxx:24-27`).
    pub fn nb_variables(&self) -> i32 {
        2
    }

    /// `_Value` (`cxx:29-46`): the distance `|C2(v) - C1(u)|`; `None` when the
    /// parameters leave the curves' ranges (`cxx:38-42`), which is how OCCT
    /// reports `false` to `math_GlobOptMin`.
    pub fn value_uv(&self, u: f64, v: f64) -> Option<f64> {
        if u < self.c1.first_parameter()
            || u > self.c1.last_parameter()
            || v < self.c2.first_parameter()
            || v > self.c2.last_parameter()
        {
            return None;
        }
        Some(self.c2.d0(v).distance(&self.c1.d0(u)))
    }

    /// `_Gradient` (`cxx:69-93`), the gradient of the squared distance:
    /// `G1 = -2 (C2(u) - C1(u))·C1'`, `G2 = +2 (C2 - C1)·C2'`.
    pub fn gradient_uv(&self, u: f64, v: f64) -> Option<(f64, f64)> {
        if u < self.c1.first_parameter()
            || u > self.c1.last_parameter()
            || v < self.c2.first_parameter()
            || v > self.c2.last_parameter()
        {
            return None;
        }
        let (p1, d1) = self.c1.d1(u);
        let (p2, d2) = self.c2.d1(v);
        let d = GpVec::from_pnts(&p1, &p2);
        let g1 = -d.dot(&d1);
        let g2 = d.dot(&d2);
        Some((2.0 * g1, 2.0 * g2))
    }

    /// `Values(X, F, G)` = `Value && Gradient` (`cxx:365-371`).
    pub fn values_uv(&self, u: f64, v: f64) -> Option<(f64, (f64, f64))> {
        Some((self.value_uv(u, v)?, self.gradient_uv(u, v)?))
    }

    /// `_Hessian` (`cxx:121-153`), the Hessian of the squared distance:
    /// `H11 = 2(|C1'|² - (C2-C1)·C1'')`, `H12 = H21 = -2 C2'·C1'`,
    /// `H22 = 2(|C2'|² + (C2-C1)·C2'')`.
    pub fn hessian_uv(&self, u: f64, v: f64) -> Option<[[f64; 2]; 2]> {
        if u < self.c1.first_parameter()
            || u > self.c1.last_parameter()
            || v < self.c2.first_parameter()
            || v > self.c2.last_parameter()
        {
            return None;
        }
        let (p1, d1, d1b) = self.c1.d2(u);
        let (p2, d2, d2b) = self.c2.d2(v);
        let d = GpVec::from_pnts(&p1, &p2);
        let h11 = 2.0 * (d1.square_magnitude() - d.dot(&d1b));
        let h12 = -2.0 * d2.dot(&d1);
        let h22 = 2.0 * (d2.square_magnitude() + d.dot(&d2b));
        Some([[h11, h12], [h12, h22]])
    }
}

impl MathMultipleVarFunction for GlobOptFuncCCC2<'_> {
    /// `NbVariables()` (`cxx:330-333`).
    fn nb_variables(&self) -> usize {
        2
    }

    /// `Value(X, F)` (`cxx:337-347`).
    fn value(&self, x: &MathVector, f: &mut f64) -> bool {
        match GlobOptFuncCCC2::value_uv(self, x.value(1), x.value(2)) {
            Some(v) => {
                *f = v;
                true
            }
            None => false,
        }
    }

    /// `dynamic_cast<math_MultipleVarFunctionWithGradient*>` succeeds
    /// (`Extrema_GlobOptFuncCC.hxx:73`).
    fn as_gradient(&self) -> Option<&dyn MathMultipleVarFunctionWithGradient> {
        Some(self)
    }

    /// `dynamic_cast<math_MultipleVarFunctionWithHessian*>` succeeds
    /// (`Extrema_GlobOptFuncCC.hxx:73`).
    fn as_hessian(&self) -> Option<&dyn MathMultipleVarFunctionWithHessian> {
        Some(self)
    }
}

impl MathMultipleVarFunctionWithGradient for GlobOptFuncCCC2<'_> {
    /// `Gradient(X, G)` (`cxx:351-361`).
    fn gradient(&self, x: &MathVector, g: &mut MathVector) -> bool {
        match GlobOptFuncCCC2::gradient_uv(self, x.value(1), x.value(2)) {
            Some((g1, g2)) => {
                g.set_value(1, g1);
                g.set_value(2, g2);
                true
            }
            None => false,
        }
    }

    /// `Values(X, F, G)` (`cxx:365-371`).
    fn values(&self, x: &MathVector, f: &mut f64, g: &mut MathVector) -> bool {
        MathMultipleVarFunction::value(self, x, f)
            && MathMultipleVarFunctionWithGradient::gradient(self, x, g)
    }
}

impl MathMultipleVarFunctionWithHessian for GlobOptFuncCCC2<'_> {
    /// `Values(X, F, G, H)` (`cxx:373-388`).
    fn values_hessian(
        &self,
        x: &MathVector,
        f: &mut f64,
        g: &mut MathVector,
        h: &mut MathMatrix,
    ) -> bool {
        let is_hessian_computed = match GlobOptFuncCCC2::hessian_uv(self, x.value(1), x.value(2)) {
            Some(m) => {
                h.set_value(1, 1, m[0][0]);
                h.set_value(1, 2, m[0][1]);
                h.set_value(2, 1, m[1][0]);
                h.set_value(2, 2, m[1][1]);
                true
            }
            None => false,
        };

        // return (Value(X, F) && Gradient(X, G) && isHessianComputed);
        MathMultipleVarFunction::value(self, x, f)
            && MathMultipleVarFunctionWithGradient::gradient(self, x, g)
            && is_hessian_computed
    }
}
