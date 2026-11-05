//! `Bisector_FunctionInter` (`Bisector_FunctionInter.cxx:19-106`,
//! `.hxx:30-58`).
//!
//! The scalar function whose zero gives the parameter at which two bisectors
//! that share a support curve meet:
//!
//! ```text
//! F(X) = ||PC(X) - PB1(X)|| - ||PC(X) - PB2(X)||
//! ```
//!
//! where `PC` is the common (guide) curve and `PB1` / `PB2` are the two
//! bisectors, all three parameterised by `X`. `Bisector_Inter::NeighbourPerform`
//! feeds it to `math_BissecNewton` (`Bisector_Inter.cxx:373-377`).
//!
//! The OCCT class stores the three `occ::handle`s; the port borrows them, which
//! keeps the same aliasing (no copy of the bisector state is made).

use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};
use occt_core::precision::{INFINITE, RESOLUTION};

use crate::curve::Curve2d;

/// `class Bisector_FunctionInter` (`Bisector_FunctionInter.hxx:30-58`).
#[derive(Clone, Copy)]
pub struct BisectorFunctionInter<'a> {
    /// `curve` (`.hxx:53`), the guide shared by the two bisectors.
    curve: Option<&'a dyn Curve2d>,
    /// `bisector1` (`.hxx:54`).
    bisector1: Option<&'a dyn Curve2d>,
    /// `bisector2` (`.hxx:55`).
    bisector2: Option<&'a dyn Curve2d>,
}

impl Default for BisectorFunctionInter<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> BisectorFunctionInter<'a> {
    /// `Bisector_FunctionInter()` (`cxx:29`): defaulted, so all three handles
    /// stay null. `Value` / `Values` on such an instance dereference a null
    /// handle in OCCT; the port reports the null curve instead.
    pub fn new() -> Self {
        Self {
            curve: None,
            bisector1: None,
            bisector2: None,
        }
    }

    /// `Bisector_FunctionInter(C, B1, B2)` (`cxx:33-41`).
    pub fn with_curves(
        curve: &'a dyn Curve2d,
        bisector1: &'a dyn Curve2d,
        bisector2: &'a dyn Curve2d,
    ) -> Self {
        Self {
            curve: Some(curve),
            bisector1: Some(bisector1),
            bisector2: Some(bisector2),
        }
    }

    /// `Perform(C, B1, B2)` (`cxx:45-53`).
    pub fn perform(
        &mut self,
        curve: &'a dyn Curve2d,
        bisector1: &'a dyn Curve2d,
        bisector2: &'a dyn Curve2d,
    ) {
        self.curve = Some(curve);
        self.bisector1 = Some(bisector1);
        self.bisector2 = Some(bisector2);
    }

    fn guide(&self) -> &dyn Curve2d {
        self.curve
            .expect("Bisector_FunctionInter: null curve")
    }

    fn first_bisector(&self) -> &dyn Curve2d {
        self.bisector1
            .expect("Bisector_FunctionInter: null bisector1")
    }

    fn second_bisector(&self) -> &dyn Curve2d {
        self.bisector2
            .expect("Bisector_FunctionInter: null bisector2")
    }

    /// `Value(X, F)` body (`cxx:54-66`).
    fn eval_value(&self, x: f64) -> f64 {
        let pc = self.guide().d0(x);
        let pb1 = self.first_bisector().d0(x);
        let pb2 = self.second_bisector().d0(x);
        pc.distance(&pb1) - pc.distance(&pb2)
    }

    /// `Values(X, F, D)` body (`cxx:73-105`).
    fn eval_values(&self, x: f64) -> (f64, f64) {
        let (pc, tc) = self.guide().d1(x);
        let (pb1, tb1) = self.first_bisector().d1(x);
        let (pb2, tb2) = self.second_bisector().d1(x);

        let f1 = pc.distance(&pb1);
        let f2 = pc.distance(&pb2);
        let f = f1 - f2;

        // `cxx:88-97`: `|F1| < gp::Resolution()` answers `Precision::Infinite()`.
        let df1 = if f1.abs() < RESOLUTION {
            INFINITE
        } else {
            ((pc.x() - pb1.x()) * (tc.x() - tb1.x())
                + (pc.y() - pb1.y()) * (tc.y() - tb1.y()))
                / f1
        };
        let df2 = if f2.abs() < RESOLUTION {
            INFINITE
        } else {
            ((pc.x() - pb2.x()) * (tc.x() - tb2.x())
                + (pc.y() - pb2.y()) * (tc.y() - tb2.y()))
                / f2
        };
        (f, df1 - df2)
    }
}

impl MathFunction for BisectorFunctionInter<'_> {
    /// `Bisector_FunctionInter::Value` (`cxx:54-66`): always returns `true`.
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        *f = self.eval_value(x);
        true
    }
}

impl MathFunctionWithDerivative for BisectorFunctionInter<'_> {
    /// `Bisector_FunctionInter::Derivative` (`cxx:69-71`): forwards to `Values`.
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let (_, dv) = self.eval_values(x);
        *d = dv;
        true
    }

    /// `Bisector_FunctionInter::Values` (`cxx:73-105`).
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let (fv, dv) = self.eval_values(x);
        *f = fv;
        *d = dv;
        true
    }
}

impl occt_math::FunctionWithDerivative for BisectorFunctionInter<'_> {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        *f = self.eval_value(x);
        true
    }

    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let (_, dv) = self.eval_values(x);
        *d = dv;
        true
    }

    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let (fv, dv) = self.eval_values(x);
        *f = fv;
        *d = dv;
        true
    }
}
