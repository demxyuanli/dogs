//! `Geom2dGcc_FunctionTanCuCu`
//! (`Geom2dGcc_FunctionTanCuCu.{hxx:33-80,cxx:49-213}`): the two-variable
//! system of the curve/curve branch of `Geom2dGcc_Lin2d2TanIter`.
//!
//! For `P1 = Curve1(X1)`, `P2 = Curve2(X2)` the system looks for the segment
//! `P1P2` normal to `Curve1` at `P1` and for parallel tangents:
//!
//! ```text
//! F1 = P1P2 x T1 / (|T1| * |P1P2|^2)
//! F2 = T1 x T2 / (|T1| * |T2|)
//! ```
//!
//! Translation notes:
//! - `InitDerivative` (`cxx:22-48`) has a `default:` arm that leaves the six
//!   output points perpendicular to the input; `TheType` is only ever
//!   `Geom2dGcc_CuCu` or `Geom2dGcc_CiCu` (`cxx:56`, `:63`), so the port models
//!   [`Geom2dGccType3`] with exactly those two variants and drops the
//!   unreachable arm,
//! - the `CiCu` arm evaluates the circle through `ElCLib::D2(X(1), TheCirc1,
//!   ...)` (`cxx:43`), i.e. the parameter is the circle's own angular
//!   parameter, not an arc length,
//! - the second constructor (`cxx:58-64`) leaves `TheCurve1` at its default
//!   (null) adaptor; the port keeps the field but only reads it in the `CuCu`
//!   arm,
//! - `F1` divides by `NormeD11 * squaredir` where `squaredir` is the
//!   **squared** length of `P1P2` (`cxx:140-142`), while the derivative
//!   formulas keep OCCT's redundant `NormeD11 * NormeD11` numerator and
//!   `NormeD11^3 * squaredir^3` denominator (`cxx:166-177`); the port keeps
//!   both spellings verbatim instead of cancelling them,
//! - `Vect11.XY();` / `Vect21.XY();` under `#else` (`cxx:160-161`, `:194-195`)
//!   are no-ops, kept as unused reads of the vectors,
//! - all three methods return `true` unconditionally; a curve that cannot be
//!   evaluated is a null-handle dereference in OCCT, so the port panics.
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::circle_d2_ax22d;
use occt_core::gp::{GpCirc2d, GpPnt2d, GpVec2d};
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::adaptor_curve::Geom2dAdaptorCurve;
use crate::gcc::curve_tool;
use crate::gcc::geom2d_gcc_type::Geom2dGccType3;

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.cxx:83-90`): `P2 - P1`.
fn vec_between(from: &GpPnt2d, to: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(to.x() - from.x(), to.y() - from.y())
}

/// `Geom2dGcc_FunctionTanCuCu` (`Geom2dGcc_FunctionTanCuCu.hxx:33-80`).
pub struct Geom2dGccFunctionTanCuCu {
    /// `TheCurve1` (`hxx:76`).
    the_curve1: Geom2dAdaptorCurve,
    /// `TheCurve2` (`hxx:77`).
    the_curve2: Geom2dAdaptorCurve,
    /// `TheCirc1` (`hxx:78`).
    the_circ1: GpCirc2d,
    /// `TheType` (`hxx:79`).
    the_type: Geom2dGccType3,
}

/// The six outputs of `InitDerivative` (`cxx:22-48`).
struct InitValues {
    point1: GpPnt2d,
    point2: GpPnt2d,
    tan1: GpVec2d,
    tan2: GpVec2d,
    d21: GpVec2d,
    d22: GpVec2d,
}

impl Geom2dGccFunctionTanCuCu {
    /// `Geom2dGcc_FunctionTanCuCu(C1, C2)` (`cxx:49-56`).
    pub fn new_two_curves(c1: &Geom2dAdaptorCurve, c2: &Geom2dAdaptorCurve) -> Self {
        Self {
            the_curve1: c1.clone(),
            the_curve2: c2.clone(),
            the_circ1: GpCirc2d::default(),
            the_type: Geom2dGccType3::CuCu,
        }
    }

    /// `Geom2dGcc_FunctionTanCuCu(C1, C2)` with a circle as first operand
    /// (`cxx:58-64`).
    pub fn new_circle_curve(c1: &GpCirc2d, c2: &Geom2dAdaptorCurve) -> Self {
        Self {
            the_curve1: Geom2dAdaptorCurve::default(),
            the_curve2: c2.clone(),
            the_circ1: *c1,
            the_type: Geom2dGccType3::CiCu,
        }
    }

    /// `TheType` (`hxx:79`).
    pub fn curve_type(&self) -> Geom2dGccType3 {
        self.the_type
    }

    /// `InitDerivative(X, Point1, Point2, Tan1, Tan2, D21, D22)`
    /// (`cxx:22-48`).
    fn init_derivative(&self, x: &MathVector) -> InitValues {
        let (point1, tan1, d21) = match self.the_type {
            Geom2dGccType3::CuCu => curve_tool::d2(&self.the_curve1, x.value(1))
                .expect("Geom2dAdaptor_Curve::D2"),
            Geom2dGccType3::CiCu => circle_d2_ax22d(
                x.value(1),
                self.the_circ1.position(),
                self.the_circ1.radius(),
            ),
        };
        let (point2, tan2, d22) = curve_tool::d2(&self.the_curve2, x.value(2))
            .expect("Geom2dAdaptor_Curve::D2");
        InitValues {
            point1,
            point2,
            tan1,
            tan2,
            d21,
            d22,
        }
    }

    /// `Value` (`cxx:128-145`) on an already initialised set.
    fn fill_value(v: &InitValues, fval: &mut MathVector) {
        let norme_d11 = v.tan1.magnitude();
        let norme_d21 = v.tan2.magnitude();
        let direction = vec_between(&v.point1, &v.point2);
        let squaredir = direction.dot(&direction);
        fval.set_value(1, direction.crossed(&v.tan1) / (norme_d11 * squaredir));
        fval.set_value(2, v.tan1.crossed(&v.tan2) / (norme_d11 * norme_d21));
    }

    /// `Derivatives` (`cxx:147-178`) on an already initialised set.
    fn fill_derivatives(v: &InitValues, deriv: &mut MathMatrix) {
        let norme_d11 = v.tan1.magnitude();
        let norme_d21 = v.tan2.magnitude();
        let direction = vec_between(&v.point1, &v.point2);
        let squaredir = direction.dot(&direction);
        deriv.set_value(
            1,
            1,
            direction.crossed(&v.d21) / (norme_d11 * squaredir)
                + (direction.crossed(&v.tan1)
                    * norme_d11
                    * norme_d11
                    * v.tan1.dot(&direction))
                    / (norme_d11
                        * norme_d11
                        * norme_d11
                        * squaredir
                        * squaredir
                        * squaredir),
        );
        deriv.set_value(
            1,
            2,
            v.tan2.crossed(&v.tan1) / (norme_d11 * squaredir)
                - (direction.crossed(&v.tan1)
                    * norme_d11
                    * norme_d11
                    * v.tan2.dot(&direction))
                    / (norme_d11
                        * norme_d11
                        * norme_d11
                        * squaredir
                        * squaredir
                        * squaredir),
        );
        deriv.set_value(
            2,
            1,
            (v.d21.crossed(&v.tan2)) / (norme_d11 * norme_d21)
                - (v.tan1.crossed(&v.tan2))
                    * (v.d21.dot(&v.tan1))
                    * norme_d21
                    * norme_d21
                    / (norme_d11
                        * norme_d11
                        * norme_d11
                        * norme_d21
                        * norme_d21
                        * norme_d21),
        );
        deriv.set_value(
            2,
            2,
            (v.tan1.crossed(&v.d22)) / (norme_d11 * norme_d21)
                - (v.tan1.crossed(&v.tan2))
                    * (v.d22.dot(&v.tan2))
                    * norme_d11
                    * norme_d11
                    / (norme_d11
                        * norme_d11
                        * norme_d11
                        * norme_d21
                        * norme_d21
                        * norme_d21),
        );
    }
}

impl MathFunctionSetWithDerivatives for Geom2dGccFunctionTanCuCu {
    /// `NbVariables()` (`cxx:117-120`).
    fn nb_variables(&self) -> usize {
        2
    }

    /// `NbEquations()` (`cxx:122-126`).
    fn nb_equations(&self) -> usize {
        2
    }

    /// `Value(X, Fval)` (`cxx:128-145`).
    fn value(&mut self, x: &MathVector, fval: &mut MathVector) -> bool {
        let v = self.init_derivative(x);
        // `Vect11.XY(); Vect21.XY();` under `#else` (`cxx:160-161`).
        let _ = (v.tan1.xy(), v.tan2.xy());
        Self::fill_value(&v, fval);
        true
    }

    /// `Derivatives(X, Deriv)` (`cxx:147-178`).
    fn derivatives(&mut self, x: &MathVector, deriv: &mut MathMatrix) -> bool {
        let v = self.init_derivative(x);
        let _ = (v.tan1.xy(), v.tan2.xy());
        Self::fill_derivatives(&v, deriv);
        true
    }

    /// `Values(X, Fval, Deriv)` (`cxx:180-213`).
    fn values(&mut self, x: &MathVector, fval: &mut MathVector, deriv: &mut MathMatrix) -> bool {
        let v = self.init_derivative(x);
        let _ = (v.tan1.xy(), v.tan2.xy());
        Self::fill_value(&v, fval);
        Self::fill_derivatives(&v, deriv);
        true
    }
}
