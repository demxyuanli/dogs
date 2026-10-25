//! The three "point / restriction" inverse families of the `Blend` package.
//!
//! Source (OCCT 8.0.0), `ModelingAlgorithms/TKFillet/Blend/`:
//! - `Blend_SurfPointFuncInv` (`Blend_SurfPointFuncInv.hxx:33-88`,
//!   `.cxx:21-24`), `X = (w, U, V)`,
//! - `Blend_SurfCurvFuncInv` (`Blend_SurfCurvFuncInv.hxx:33-88`, `.cxx:22-25`),
//!   `X = (wguide, wcurv, wrst)`,
//! - `Blend_CurvPointFuncInv` (`Blend_CurvPointFuncInv.hxx:31-87`, `.cxx:21-24`),
//!   `X = (w, U, V)`.
//!
//! The three have the same shape — `NbVariables`, `NbEquations`, `Value`,
//! `Derivatives`, `Values`, `Set`, `GetTolerance`, `GetBounds`, `IsSolution` —
//! and differ only in the `Set` argument and in the `NbVariables` value (all
//! three return `3` in OCCT except `Blend_CurvPointFuncInv`, which returns `2`).
//! They stay three traits because `BRepBlend_*Inv` classes implement exactly
//! one of them and the walking builders dispatch on that type.
//!
//! Translation notes:
//! - `NbVariables` is a body at this level in each `.cxx`, but it redefines a
//!   required base method, so an implementor writes it: `3`, `3` and `2`
//!   respectively. See [`super::app_function`] for the convention.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathVector};

/// `class Blend_SurfPointFuncInv` (`hxx:33-88`): a solution on a given point of
/// the curve.
pub trait BlendSurfPointFuncInv: MathFunctionSetWithDerivatives {
    /// `Set(P)` (`hxx:70`): the point on which a solution has to be found.
    fn set(&mut self, p: GpPnt);

    /// `GetTolerance(Tolerance, Tol)` (`hxx:75`): the parametric tolerance of
    /// each of the three variables.
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64);

    /// `GetBounds(InfBound, SupBound)` (`hxx:81`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector);

    /// `IsSolution(Sol, Tol)` (`hxx:85`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool;
}

/// `class Blend_SurfCurvFuncInv` (`hxx:33-88`): a solution on a given
/// restriction of the surface.
pub trait BlendSurfCurvFuncInv: MathFunctionSetWithDerivatives {
    /// `Set(Rst)` (`hxx:70`): the restriction on which a solution has to be
    /// found. OCCT does not guard the handle; the port keeps the null case as
    /// `Option`.
    fn set(&mut self, rst: Option<Arc<dyn Curve2d>>);

    /// `GetTolerance(Tolerance, Tol)` (`hxx:75`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64);

    /// `GetBounds(InfBound, SupBound)` (`hxx:81`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector);

    /// `IsSolution(Sol, Tol)` (`hxx:85`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool;
}

/// `class Blend_CurvPointFuncInv` (`hxx:31-87`): a solution on a given point of
/// the curve; `NbVariables` is `2` (`cxx:18-21`).
pub trait BlendCurvPointFuncInv: MathFunctionSetWithDerivatives {
    /// `Set(P)` (`hxx:69`).
    fn set(&mut self, p: GpPnt);

    /// `GetTolerance(Tolerance, Tol)` (`hxx:75`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64);

    /// `GetBounds(InfBound, SupBound)` (`hxx:81`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector);

    /// `IsSolution(Sol, Tol)` (`hxx:85`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool;
}
