//! `Extrema_GlobOptFuncCS` — the squared point-to-point distance between a
//! curve and a surface, in the three parameters `(cu, su, sv)`.
//! Source: `Extrema_GlobOptFuncCS.cxx:28-174`, `.hxx:27-58`.
//!
//! ```text
//! F(cu,su,sv) = (C(cu) - S(su,sv))^2
//! ```
//!
//! `Value` / `Gradient` / `Values` all start with `checkInputData`
//! (`cxx:90-100`): an out-of-range parameter answers `false`, which is what
//! `math_PSO` reads as "no value at this particle".

use occt_math::{
    MathMatrix, MathMultipleVarFunction, MathMultipleVarFunctionWithGradient,
    MathMultipleVarFunctionWithHessian, MathVector,
};

use super::prelude::*;

/// `Extrema_GlobOptFuncCS` (`Extrema_GlobOptFuncCS.hxx:27-58`).
pub struct GlobOptFuncCS<'a> {
    /// `myC`.
    c: &'a dyn Curve,
    /// `myS`.
    s: &'a dyn Surface,
}

impl<'a> GlobOptFuncCS<'a> {
    /// `Extrema_GlobOptFuncCS(C, S)` (`cxx:104-108`).
    pub fn new(c: &'a dyn Curve, s: &'a dyn Surface) -> Self {
        GlobOptFuncCS { c, s }
    }

    /// `checkInputData(X, cu, su, sv)` (`cxx:90-100`): `None` is the `false`
    /// return.
    fn check_input_data(&self, x: &MathVector) -> Option<(f64, f64, f64)> {
        let start = 1;
        let cu = x.value(start);
        let su = x.value(start + 1);
        let sv = x.value(start + 2);
        let inside = cu >= self.c.first_parameter()
            && cu <= self.c.last_parameter()
            && su >= self.s.u_range().0
            && su <= self.s.u_range().1
            && sv >= self.s.v_range().0
            && sv <= self.s.v_range().1;
        if inside {
            Some((cu, su, sv))
        } else {
            None
        }
    }

    /// `value(cu, su, sv, F)` (`cxx:28-31`).
    fn value_at(&self, cu: f64, su: f64, sv: f64) -> f64 {
        self.c.d0(cu).square_distance(&self.s.d0(su, sv))
    }

    /// `gradient(cu, su, sv, G)` (`cxx:35-49`).
    fn gradient_at(&self, cu: f64, su: f64, sv: f64, g: &mut MathVector) {
        let (cd0, cd1) = self.c.d1(cu);
        let (sd0, sd1u, sd1v) = self.s.d1(su, sv);
        let cd0_sd0 = GpVec::from_pnts(&sd0, &cd0);
        g.set_value(1, cd0_sd0.dot(&cd1));
        g.set_value(2, -cd0_sd0.dot(&sd1u));
        g.set_value(3, -cd0_sd0.dot(&sd1v));
    }

    /// `hessian(cu, su, sv, H)` (`cxx:54-87`).
    fn hessian_at(&self, cu: f64, su: f64, sv: f64, h: &mut MathMatrix) {
        let (cd0, cd1, cd2) = self.c.d2(cu);
        let (sd0, sd1u, sd1v, sd2uu, sd2vv, sd2uv) = self.s.d2(su, sv);
        let cd0_sd0 = GpVec::from_pnts(&sd0, &cd0);

        h.set_value(
            1,
            1,
            cd1.x() * cd1.x() + cd1.y() * cd1.y() + cd1.z() * cd1.z() + cd0_sd0.dot(&cd2),
        );
        h.set_value(1, 2, -cd1.dot(&sd1u));
        h.set_value(1, 3, -cd1.dot(&sd1v));

        h.set_value(2, 1, h.value(1, 2));
        h.set_value(
            2,
            2,
            sd1u.x() * sd1u.x() + sd1u.y() * sd1u.y() + sd1u.z() * sd1u.z()
                - cd0_sd0.dot(&sd2uu),
        );
        h.set_value(
            2,
            3,
            sd1u.dot(&sd1v) - cd0_sd0.dot(&sd2uv),
        );

        h.set_value(3, 1, h.value(1, 3));
        h.set_value(3, 2, h.value(2, 3));
        h.set_value(
            3,
            3,
            sd1v.x() * sd1v.x() + sd1v.y() * sd1v.y() + sd1v.z() * sd1v.z()
                - cd0_sd0.dot(&sd2vv),
        );
    }
}

impl MathMultipleVarFunction for GlobOptFuncCS<'_> {
    /// `NbVariables()` (`cxx:112-115`): 3.
    fn nb_variables(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:119-127`).
    fn value(&self, x: &MathVector, f: &mut f64) -> bool {
        match self.check_input_data(x) {
            Some((cu, su, sv)) => {
                *f = self.value_at(cu, su, sv);
                true
            }
            None => false,
        }
    }
}

impl MathMultipleVarFunctionWithGradient for GlobOptFuncCS<'_> {
    /// `Gradient(X, G)` (`cxx:135-141`).
    fn gradient(&self, x: &MathVector, g: &mut MathVector) -> bool {
        match self.check_input_data(x) {
            Some((cu, su, sv)) => {
                self.gradient_at(cu, su, sv, g);
                true
            }
            None => false,
        }
    }

    /// `Values(X, F, G)` (`cxx:147-157`).
    fn values(&self, x: &MathVector, f: &mut f64, g: &mut MathVector) -> bool {
        match self.check_input_data(x) {
            Some((cu, su, sv)) => {
                *f = self.value_at(cu, su, sv);
                self.gradient_at(cu, su, sv, g);
                true
            }
            None => false,
        }
    }
}

impl MathMultipleVarFunctionWithHessian for GlobOptFuncCS<'_> {
    /// `Values(X, F, G, H)` (`cxx:162-174`).
    fn values_hessian(
        &self,
        x: &MathVector,
        f: &mut f64,
        g: &mut MathVector,
        h: &mut MathMatrix,
    ) -> bool {
        match self.check_input_data(x) {
            Some((cu, su, sv)) => {
                *f = self.value_at(cu, su, sv);
                self.gradient_at(cu, su, sv, g);
                self.hessian_at(cu, su, sv, h);
                true
            }
            None => false,
        }
    }
}
