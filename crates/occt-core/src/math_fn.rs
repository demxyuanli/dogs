//! `math_Function` / `math_FunctionWithDerivative`.
//! Source: `math_Function.hxx`, `math_FunctionWithDerivative.hxx`.

/// `math_Function` (`math_Function.hxx:29-58`).
pub trait MathFunction {
    /// `Value(X, F)`. Writes `F` and returns whether the evaluation succeeded.
    fn value(&mut self, x: f64, f: &mut f64) -> bool;

    /// `GetStateNumber` (`math_Function.cxx:17-20`). Default 0.
    fn get_state_number(&mut self) -> i32 {
        0
    }
}

/// `math_FunctionWithDerivative` (`math_FunctionWithDerivative.hxx:31-52`).
pub trait MathFunctionWithDerivative: MathFunction {
    /// `Derivative(X, D)`.
    fn derivative(&mut self, x: f64, d: &mut f64) -> bool;

    /// `Values(X, F, D)`.
    fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool;
}
