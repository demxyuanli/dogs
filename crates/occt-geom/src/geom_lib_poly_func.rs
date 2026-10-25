//! `GeomLib_PolyFunc` — a `math_FunctionWithDerivative` over the derivative of
//! a polynomial given by its coefficients.
//! Source: `ModelingData/TKGeomBase/GeomLib/GeomLib_PolyFunc.cxx:21-58`
//! (OCCT 8.0.0), declaration `GeomLib_PolyFunc.hxx:26-48`.
//!
//! The constructor keeps the *derivative* coefficients (`cxx:26-34`:
//! `myCoeffs(ii) = ii * Coeffs(ii + 1)`), and `Value` / `Derivative` /
//! `Values` evaluate that polynomial with `PLib::EvalPolynomial(..., Degree =
//! myCoeffs.Length() - 1, Dimension = 1, ...)` (`cxx:36-58`).
//!
//! `Value` / `Derivative` / `Values` always answer `true` in OCCT
//! (`cxx:44`, `:53`, `cxx:66`), so the trait methods never fail.

use occt_core::bspl::coefs::eval_polynomial;
use occt_core::math_fn::{MathFunction, MathFunctionWithDerivative};

/// `GeomLib_PolyFunc` (`hxx:26-48`).
#[derive(Debug, Clone)]
pub struct PolyFunc {
    /// `myCoeffs` (`hxx:44`), the derivative coefficients, 0-based: entry `i`
    /// multiplies `X^i`.
    coeffs: Vec<f64>,
}

impl PolyFunc {
    /// `GeomLib_PolyFunc(Coeffs)` (`cxx:21-27`): `myCoeffs(1, Length-1)` with
    /// `myCoeffs(ii) = ii * Coeffs(ii + 1)`.
    pub fn new(coeffs: &[f64]) -> Self {
        let n = coeffs.len();
        let mut derivative = Vec::with_capacity(n.saturating_sub(1));
        // `for (ii = 1; ii <= myCoeffs.Length(); ii++)` with
        // `myCoeffs.Length() = Coeffs.Length() - 1` (`cxx:24-33`): 0-based this
        // is `derivative[ii - 1] = ii * coeffs[ii]` for `ii = 1..=n - 1`.
        for ii in 1..n {
            derivative.push(ii as f64 * coeffs[ii]);
        }
        Self { coeffs: derivative }
    }

    /// `PLib::EvalPolynomial(X, Request, Length - 1, 1, coeff, ff)` with the
    /// 1-based `Length - 1` clamped the way the OCCT call would read it.
    fn eval(&self, x: f64, derivative_request: i32) -> Vec<f64> {
        let degree = self.coeffs.len() as i32 - 1;
        eval_polynomial(x, derivative_request, degree, 1, &self.coeffs)
    }
}

impl MathFunction for PolyFunc {
    /// `Value` (`cxx:29-36`).
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        *f = self.eval(x, 0)[0];
        true
    }
}

impl MathFunctionWithDerivative for PolyFunc {
    /// `Derivative` (`cxx:38-46`): `Aux(2)` of a two-slot evaluation.
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        *d = self.eval(x, 1)[1];
        true
    }

    /// `Values` (`cxx:48-58`): `F = Aux(1)`, `D = Aux(2)`.
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let aux = self.eval(x, 1);
        *f = aux[0];
        *d = aux[1];
        true
    }
}
