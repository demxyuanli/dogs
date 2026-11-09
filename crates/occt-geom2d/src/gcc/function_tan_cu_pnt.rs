//! `Geom2dGcc_FunctionTanCuPnt`
//! (`Geom2dGcc_FunctionTanCuPnt.{hxx:32-57,cxx:38-90}`): the one-variable
//! system of the curve/point branch of `Geom2dGcc_Lin2d2TanIter`.
//!
//! For `P = Curve(u)` and `C = ThePoint` the system looks for the segment `CP`
//! normal to the curve:
//!
//! ```text
//! F(u) = CP x T / (|T| * |CP|)
//! ```
//!
//! Translation notes:
//! - `TheDirection` is `gp_Vec2d(ThePoint, Point)` (`cxx:50`, `:62`) i.e.
//!   `Point - ThePoint`, and `|CP|` is its magnitude, not a squared one,
//! - `Derivative` (`cxx:57-70`) and `Values` (`cxx:72-90`) share the same
//!   expression; the port keeps one helper for both,
//! - the `std::cout` line (`cxx:87`) is commented out in OCCT and has no port,
//! - all three methods return `true` unconditionally; a curve that cannot be
//!   evaluated is a null-handle dereference in OCCT, so the port panics.
//!
//! UNPORTED: none.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_math::FunctionWithDerivative;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.cxx:83-90`): `P2 - P1`.
fn vec_between(from: &GpPnt2d, to: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(to.x() - from.x(), to.y() - from.y())
}

/// `Geom2dGcc_FunctionTanCuPnt` (`Geom2dGcc_FunctionTanCuPnt.hxx:32-57`).
pub struct Geom2dGccFunctionTanCuPnt {
    /// `TheCurv` (`hxx:56`).
    the_curv: Geom2dAdaptorCurve,
    /// `ThePoint` (`hxx:57`).
    the_point: GpPnt2d,
}

impl Geom2dGccFunctionTanCuPnt {
    /// `Geom2dGcc_FunctionTanCuPnt(C, Point)` (`cxx:38-42`).
    pub fn new(c: &Geom2dAdaptorCurve, point: &GpPnt2d) -> Self {
        Self {
            the_curv: c.clone(),
            the_point: *point,
        }
    }

    /// The shared body of `Derivative` (`cxx:57-70`) and `Values`
    /// (`cxx:72-90`).
    fn deriv_at(&self, point: &GpPnt2d, vec1: &GpVec2d, vec2: &GpVec2d) -> f64 {
        let direction = vec_between(&self.the_point, point);
        let norme_d1 = vec1.magnitude();
        let norme_dir = direction.magnitude();
        direction.crossed(vec2) / (norme_d1 * norme_dir)
            - (direction.crossed(vec1) / (norme_d1 * norme_dir))
                * (vec1.dot(vec2) / (norme_d1 * norme_d1)
                    + vec1.dot(&direction) / (norme_dir * norme_dir))
    }
}

impl FunctionWithDerivative for Geom2dGccFunctionTanCuPnt {
    /// `Value(X, Fval)` (`cxx:44-55`).
    fn value(&mut self, x: f64, fval: &mut f64) -> bool {
        let (point, vect) = curve_tool::d1(&self.the_curv, x).expect("Geom2dAdaptor_Curve::D1");
        let norme_d1 = vect.magnitude();
        let direction = vec_between(&self.the_point, &point);
        let norme_dir = direction.magnitude();
        *fval = direction.crossed(&vect) / (norme_d1 * norme_dir);
        true
    }

    /// `Derivative(X, Deriv)` (`cxx:57-70`).
    fn derivative(&mut self, x: f64, deriv: &mut f64) -> bool {
        let (point, vec1, vec2) =
            curve_tool::d2(&self.the_curv, x).expect("Geom2dAdaptor_Curve::D2");
        *deriv = self.deriv_at(&point, &vec1, &vec2);
        true
    }

    /// `Values(X, Fval, Deriv)` (`cxx:72-90`).
    fn values(&mut self, x: f64, fval: &mut f64, deriv: &mut f64) -> bool {
        let (point, vec1, vec2) =
            curve_tool::d2(&self.the_curv, x).expect("Geom2dAdaptor_Curve::D2");
        let direction = vec_between(&self.the_point, &point);
        let norme_d1 = vec1.magnitude();
        let norme_dir = direction.magnitude();
        *fval = direction.crossed(&vec1) / (norme_d1 * norme_dir);
        *deriv = self.deriv_at(&point, &vec1, &vec2);
        true
    }
}
