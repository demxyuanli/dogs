//! `Geom2dGcc_FunctionTanObl`
//! (`Geom2dGcc_FunctionTanObl.{hxx:31-60,cxx:24-64}`): the one-variable system
//! of `Geom2dGcc_Lin2dTanOblIter`.
//!
//! For `T = Curve(u)` the system looks for the parameter where the tangent is
//! parallel to `TheDirection`:
//!
//! ```text
//! F(u) = TheDirection x T / |T|
//! ```
//!
//! Translation notes:
//! - `Value` (`cxx:30-38`) is `D1` reduced to the tangent,
//! - `Derivative` (`cxx:40-50`) and `Values` (`cxx:52-63`) share the same
//!   expression, which OCCT writes as
//!   `dir x T'' / |T'| - (T' . T'' * dir x T') / |T'|`; the classical
//!   derivative would divide the second term by `|T'|^3`, so the written form
//!   is kept term by term - `(T' . T'') * (dir x T')` is multiplied first and
//!   the product divided by `|T'|`,
//! - all three methods return `true` unconditionally; a curve that cannot be
//!   evaluated is a null-handle dereference in OCCT, so the port panics.
//!
//! UNPORTED: none.

use occt_core::gp::{GpDir2d, GpXY};
use occt_math::FunctionWithDerivative;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;

/// `Geom2dGcc_FunctionTanObl` (`Geom2dGcc_FunctionTanObl.hxx:31-60`).
pub struct Geom2dGccFunctionTanObl {
    /// `TheCurv` (`hxx:58`).
    the_curv: Geom2dAdaptorCurve,
    /// `TheDirection` (`hxx:59`).
    the_direction: GpDir2d,
}

impl Geom2dGccFunctionTanObl {
    /// `Geom2dGcc_FunctionTanObl(C, Dir)` (`cxx:24-28`).
    pub fn new(c: &Geom2dAdaptorCurve, dir: &GpDir2d) -> Self {
        Self {
            the_curv: c.clone(),
            the_direction: *dir,
        }
    }

    /// `TheDirection.XY()` (`cxx:36`, `:47`, `:59`).
    fn dir_xy(&self) -> GpXY {
        GpXY::new(self.the_direction.x(), self.the_direction.y())
    }

    /// The shared body of `Derivative` (`cxx:40-50`) and `Values`
    /// (`cxx:52-63`).
    fn deriv_at(&self, vec1: &occt_core::gp::GpVec2d, vec2: &occt_core::gp::GpVec2d) -> f64 {
        let dir = self.dir_xy();
        let norme_d1 = vec1.magnitude();
        dir.crossed(&vec2.xy()) / norme_d1
            - vec1.xy().dot(&vec2.xy()) * dir.crossed(&vec1.xy()) / norme_d1
    }
}

impl FunctionWithDerivative for Geom2dGccFunctionTanObl {
    /// `Value(X, Fval)` (`cxx:30-38`).
    fn value(&mut self, x: f64, fval: &mut f64) -> bool {
        let (_, vect) = curve_tool::d1(&self.the_curv, x).expect("Geom2dAdaptor_Curve::D1");
        let norme_d1 = vect.magnitude();
        *fval = self.dir_xy().crossed(&vect.xy()) / norme_d1;
        true
    }

    /// `Derivative(X, Deriv)` (`cxx:40-50`).
    fn derivative(&mut self, x: f64, deriv: &mut f64) -> bool {
        let (_, vec1, vec2) =
            curve_tool::d2(&self.the_curv, x).expect("Geom2dAdaptor_Curve::D2");
        *deriv = self.deriv_at(&vec1, &vec2);
        true
    }

    /// `Values(X, Fval, Deriv)` (`cxx:52-63`).
    fn values(&mut self, x: f64, fval: &mut f64, deriv: &mut f64) -> bool {
        let (_, vec1, vec2) =
            curve_tool::d2(&self.the_curv, x).expect("Geom2dAdaptor_Curve::D2");
        let norme_d1 = vec1.magnitude();
        *fval = self.dir_xy().crossed(&vec1.xy()) / norme_d1;
        *deriv = self.deriv_at(&vec1, &vec2);
        true
    }
}
