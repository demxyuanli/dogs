//! `math_MultipleVarFunction` (`math_MultipleVarFunction.hxx:29-52`),
//! `math_MultipleVarFunctionWithGradient`
//! (`math_MultipleVarFunctionWithGradient.hxx`) and
//! `math_MultipleVarFunctionWithHessian`
//! (`math_MultipleVarFunctionWithHessian.hxx:26-57`).
//!
//! These are the scalar-objective interfaces consumed by `math_PSO` and
//! `math_NewtonMinimum`.
//!
//! Rust note: OCCT's `Value`/`Gradient`/`Values` are non-const and are invoked
//! through mutable references. The objectives ported in this crate are pure
//! functions of the parameter vector, so the trait methods take `&self`; the
//! algorithms below (and their arithmetic) are unaffected.

use crate::{MathMatrix, MathVector};

/// `math_MultipleVarFunction` (`math_MultipleVarFunction.hxx:29`).
pub trait MathMultipleVarFunction {
    /// `NbVariables()` (`math_MultipleVarFunction.hxx:38`).
    fn nb_variables(&self) -> usize;

    /// `Value(X, F)` (`math_MultipleVarFunction.hxx:45`): returns false when the
    /// evaluation failed (OCCT leaves `F` untouched in that case).
    fn value(&self, x: &MathVector, f: &mut f64) -> bool;

    /// `GetStateNumber()` (`math_MultipleVarFunction.hxx:51`): default 0.
    fn state_number(&self) -> i32 {
        0
    }

    /// Rust stand-in for OCCT's
    /// `dynamic_cast<math_MultipleVarFunctionWithGradient*>(myFunc)`
    /// (`math_GlobOptMin.cxx:294`): the ported hierarchy is static, so the
    /// object answers for itself and the default `None` is the failing cast.
    fn as_gradient(&self) -> Option<&dyn MathMultipleVarFunctionWithGradient> {
        None
    }

    /// Rust stand-in for OCCT's
    /// `dynamic_cast<math_MultipleVarFunctionWithHessian*>(myFunc)`
    /// (`math_GlobOptMin.cxx:273`).
    fn as_hessian(&self) -> Option<&dyn MathMultipleVarFunctionWithHessian> {
        None
    }
}

/// `math_MultipleVarFunctionWithGradient`.
pub trait MathMultipleVarFunctionWithGradient: MathMultipleVarFunction {
    /// `Gradient(X, G)`.
    fn gradient(&self, x: &MathVector, g: &mut MathVector) -> bool;

    /// `Values(X, F, G)`.
    fn values(&self, x: &MathVector, f: &mut f64, g: &mut MathVector) -> bool;
}

/// `math_MultipleVarFunctionWithHessian`
/// (`math_MultipleVarFunctionWithHessian.hxx:26-57`).
pub trait MathMultipleVarFunctionWithHessian: MathMultipleVarFunctionWithGradient {
    /// `Values(X, F, G, H)` (`math_MultipleVarFunctionWithHessian.hxx:49`).
    fn values_hessian(
        &self,
        x: &MathVector,
        f: &mut f64,
        g: &mut MathVector,
        h: &mut MathMatrix,
    ) -> bool;
}
