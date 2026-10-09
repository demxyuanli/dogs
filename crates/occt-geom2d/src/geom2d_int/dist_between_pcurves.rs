//! Port of `IntCurve_DistBetweenPCurvesGen` instantiated as
//! `Geom2dInt_TheDistBetweenPCurvesOfTheIntPCurvePCurveOfGInter`
//! (`IntCurve_DistBetweenPCurvesGen.gxx:33-104`).
//!
//! The 2-equation / 2-variable `math_FunctionSetWithDerivatives` solved by
//! `IntCurve_ExactIntersectionPoint`: `F(u, v) = C1(u) - C2(v)` with the exact
//! Jacobian.

use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::curve::Curve2d;
use super::curve_tool;

/// `Geom2dInt_TheDistBetweenPCurvesOfTheIntPCurvePCurveOfGInter`.
#[derive(Clone, Copy)]
pub struct DistBetweenPCurves<'a> {
    /// `thecurve1` (`gxx:36`).
    pub curve1: &'a dyn Curve2d,
    /// `thecurve2` (`gxx:37`).
    pub curve2: &'a dyn Curve2d,
}

impl<'a> DistBetweenPCurves<'a> {
    /// `IntCurve_DistBetweenPCurvesGen(C1, C2)` (`gxx:34-38`).
    pub fn new(c1: &'a dyn Curve2d, c2: &'a dyn Curve2d) -> Self {
        Self { curve1: c1, curve2: c2 }
    }
}

impl MathFunctionSetWithDerivatives for DistBetweenPCurves<'_> {
    /// `NbVariables()` (`gxx:42-45`).
    fn nb_variables(&self) -> usize {
        2
    }

    /// `NbEquations()` (`gxx:48-51`).
    fn nb_equations(&self) -> usize {
        2
    }

    /// `Value(X, F)` (`gxx:55-63`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let p1 = curve_tool::d0(self.curve1, x.value(1));
        let p2 = curve_tool::d0(self.curve2, x.value(2));
        f.set_value(1, p1.x() - p2.x());
        f.set_value(2, p1.y() - p2.y());
        true
    }

    /// `Derivatives(X, D)` (`gxx:67-81`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (_, t) = curve_tool::d1(self.curve1, x.value(1));
        d.set_value(1, 1, t.x());
        d.set_value(2, 1, t.y());

        let (_, t) = curve_tool::d1(self.curve2, x.value(2));
        d.set_value(1, 2, -t.x());
        d.set_value(2, 2, -t.y());
        true
    }

    /// `Values(X, F, D)` (`gxx:85-104`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (p1, t) = curve_tool::d1(self.curve1, x.value(1));
        d.set_value(1, 1, t.x());
        d.set_value(2, 1, t.y());

        let (p2, t) = curve_tool::d1(self.curve2, x.value(2));
        d.set_value(1, 2, -t.x());
        d.set_value(2, 2, -t.y());

        f.set_value(1, p1.x() - p2.x());
        f.set_value(2, p1.y() - p2.y());
        true
    }
}
