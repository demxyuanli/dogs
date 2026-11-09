//! `Geom2dGcc_FunctionTanCirCu`
//! (`Geom2dGcc_FunctionTanCirCu.{hxx:30-57,cxx:50-127}`): the one-variable
//! system of the circle/curve branch of `Geom2dGcc_Circ2d3TanIter`.
//!
//! For `P1 = Curve(X)` and `C` the centre of `TheCirc` the system looks for the
//! point `P2` on the circle with `P1P2 . CP2 = 0` and `|CP2| = R`, i.e. the
//! zeros of
//!
//! ```text
//! F(X) = |CP1|^2 - (CP1.T)^2 / |T|^2 - R^2
//! ```
//!
//! divided by `myWeight`.
//!
//! Translation notes:
//! - the normalization weight (`cxx:56-70`) is the squared distance between the
//!   circle centre and the average of `NbSamples` curve points taken at the
//!   step midpoints, floored at `Radius`; the port keeps the `<=` loop guard and
//!   the `std::max` clamp,
//! - `Value` divides by `NormeD1 * NormeD1` with `NormeD1 = Vect1.Magnitude()`
//!   (`cxx:78-84`) while `Derivative` / `Values` use the already squared
//!   `NormeD1 = Vect1.SquareMagnitude()` (`cxx:95`, `:112`); the port keeps
//!   both spellings because they are the same quantity only up to rounding,
//! - `myWeight` divides the residual and the derivative but not the squared
//!   direction (`cxx:86`, `:101`, `:119`, `:127`),
//! - all three methods return `true` unconditionally; a failure to evaluate the
//!   curve is a null-handle dereference in OCCT, so the port panics.
//!
//! UNPORTED: none.

use occt_core::gp::{GpCirc2d, GpPnt2d, GpVec2d, GpXY};
use occt_math::FunctionWithDerivative;

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.cxx:83-90`): `P2 - P1`.
fn vec_between(from: &GpPnt2d, to: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(to.x() - from.x(), to.y() - from.y())
}

/// `Geom2dGcc_FunctionTanCirCu` (`Geom2dGcc_FunctionTanCirCu.hxx:30-57`).
pub struct Geom2dGccFunctionTanCirCu {
    /// `TheCirc` (`hxx:55`).
    the_circ: GpCirc2d,
    /// `Curve` (`hxx:56`).
    curve: Geom2dAdaptorCurve,
    /// `myWeight` (`hxx:57`).
    my_weight: f64,
}

impl Geom2dGccFunctionTanCirCu {
    /// `Geom2dGcc_FunctionTanCirCu(Circ, Curv)` (`cxx:50-71`).
    pub fn new(circ: &GpCirc2d, curv: &Geom2dAdaptorCurve) -> Self {
        let mut s = Self {
            the_circ: *circ,
            curve: curv.clone(),
            my_weight: 0.0,
        };
        // `Modified by Sergey KHROMOV - Thu Apr  5 09:51:21 2001` (`cxx:56-70`).
        let a_nb_samp = curve_tool::nb_samples(&s.curve) as f64;
        let a_first = curve_tool::first_parameter(&s.curve);
        let a_last = curve_tool::last_parameter(&s.curve);
        let a_step = (a_last - a_first) / a_nb_samp;
        let mut an_x = a_first + a_step / 2.0;
        let mut a_loc = GpXY::zero();
        while an_x <= a_last {
            let point = curve_tool::value(&s.curve, an_x).expect("Geom2dAdaptor_Curve::Value");
            a_loc.add(point.xy());
            an_x += a_step;
        }
        let center = s.the_circ.location();
        let delta = a_loc.subtracted(&GpXY::new(center.x(), center.y()));
        s.my_weight = delta.square_modulus().max(s.the_circ.radius());
        s
    }

    /// `myWeight` (`hxx:57`): the normalization factor.
    pub fn my_weight(&self) -> f64 {
        self.my_weight
    }

    /// `Derivative(X, Deriv)` (`cxx:90-105`) and the derivative half of
    /// `Values` (`cxx:122-127`): the same expression in both.
    fn deriv_at(&self, x: f64, vect1: &GpVec2d, vect2: &GpVec2d, point: &GpPnt2d) -> f64 {
        let norme_d1 = vect1.square_magnitude();
        let direction = vec_between(&self.the_circ.location(), point);
        let cp1dott = direction.dot(vect1);
        let deriv = -2.0
            * (cp1dott / norme_d1)
            * (direction.dot(vect2) - cp1dott * vect1.dot(vect2) / norme_d1);
        let _ = x;
        deriv / self.my_weight
    }
}

impl FunctionWithDerivative for Geom2dGccFunctionTanCirCu {
    /// `Value(X, Fval)` (`cxx:73-88`).
    fn value(&mut self, x: f64, fval: &mut f64) -> bool {
        let (point, vect1) = curve_tool::d1(&self.curve, x).expect("Geom2dAdaptor_Curve::D1");
        let norme_d1 = vect1.magnitude();
        let direction = vec_between(&self.the_circ.location(), &point);
        let squaredir = direction.dot(&direction);
        let r = self.the_circ.radius();
        let dot = direction.dot(&vect1);
        *fval = squaredir - r * r - dot * dot / (norme_d1 * norme_d1);
        *fval /= self.my_weight;
        true
    }

    /// `Derivative(X, Deriv)` (`cxx:90-105`).
    fn derivative(&mut self, x: f64, deriv: &mut f64) -> bool {
        let (point, vect1, vect2) =
            curve_tool::d2(&self.curve, x).expect("Geom2dAdaptor_Curve::D2");
        *deriv = self.deriv_at(x, &vect1, &vect2, &point);
        true
    }

    /// `Values(X, Fval, Deriv)` (`cxx:107-127`).
    fn values(&mut self, x: f64, fval: &mut f64, deriv: &mut f64) -> bool {
        let (point, vect1, vect2) =
            curve_tool::d2(&self.curve, x).expect("Geom2dAdaptor_Curve::D2");
        let norme_d1 = vect1.square_magnitude();
        let direction = vec_between(&self.the_circ.location(), &point);
        let squaredir = direction.square_magnitude();
        let cp1dott = direction.dot(&vect1);
        let r = self.the_circ.radius();
        *fval = squaredir - r * r - cp1dott * cp1dott / norme_d1;
        *fval /= self.my_weight;
        *deriv = self.deriv_at(x, &vect1, &vect2, &point);
        true
    }
}
