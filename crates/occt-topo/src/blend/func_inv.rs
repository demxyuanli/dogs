//! `Blend_FuncInv` — the deferred inverse function that finds a solution on a
//! restriction (a curve on surface) of one of the two supports.
//!
//! Source: `ModelingAlgorithms/TKFillet/Blend/Blend_FuncInv.hxx:37-94`,
//! `Blend_FuncInv.cxx:18-21` (OCCT 8.0.0).
//!
//! The vector `X` is `(t, w, U, V)`: `t` the parameter on the curve on surface,
//! `w` the parameter on the guide line, `U, V` the parametric coordinates on
//! the partner surface (`hxx:43-51`).
//!
//! Translation notes:
//! - `NbVariables()` returns `4` in OCCT (`cxx:18-21`) but redefines a required
//!   base method, so the port cannot default it in this sub-trait; an
//!   implementor writes `fn nb_variables(&self) -> usize { 4 }`. See
//!   [`super::app_function`] for the convention.
//! - `Set(OnFirst, COnSurf)` (`hxx:75-77`) takes an `Adaptor2d_Curve2d` handle
//!   that OCCT does not guard; the port keeps the null case as `Option`.

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathVector};

/// `class Blend_FuncInv : public math_FunctionSetWithDerivatives`
/// (`hxx:37-94`).
pub trait BlendFuncInv: MathFunctionSetWithDerivatives {
    /// `Set(OnFirst, COnSurf)` (`hxx:77`): the curve on surface on which a
    /// solution has to be found, on the first surface when `on_first`, on the
    /// second one otherwise.
    fn set(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>);

    /// `GetTolerance(Tolerance, Tol)` (`hxx:82`): the parametric tolerance of
    /// each of the four variables, for the 3d tolerance `tol`.
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64);

    /// `GetBounds(InfBound, SupBound)` (`hxx:88`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector);

    /// `IsSolution(Sol, Tol)` (`hxx:92`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool;
}
