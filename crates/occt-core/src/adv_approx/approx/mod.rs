//! `AdvApprox_ApproxAFunction` for a 2D or 3D subspace.
//! Source: `AdvApprox_ApproxAFunction.cxx` (constructor `cxx:603-630`,
//! `Perform` `cxx:663-956`, `Approximation` `cxx:364-599`).

/// Crate imports shared by `mod.rs` and every submodule.
mod prelude {
    pub(super) use crate::bspl::comp_poly::CompPolynomialToPoles;
    pub(super) use crate::bspl::plib_eval::eval_polynomial;
    pub(super) use crate::bspl::plib_jacobi::{jacobi_parameters, JacobiPolynomial};
    pub(super) use crate::gp::{GpPnt, GpPnt2d};
    pub(super) use crate::kernel::geomabs::Shape;
    pub(super) use crate::precision::PCONFUSION;
    pub(super) use crate::adv_approx::cutting::{Cutting, DichoCutting};
    pub(super) use crate::adv_approx::simple::{EvalFn, SimpleApprox};
}

use self::prelude::*;

/// Result of a 3D `ApproxAFunction` run (`Perform` `cxx:663-956`).
pub struct ApproxAFunction3d {
    pub done: bool,
    pub has_result: bool,
    pub poles: Vec<GpPnt>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub degree: i32,
    /// `MaxError(3, 1)` (`cxx:1004-1035`): maximum over the produced spans of
    /// the error bound accumulated by `PrepareConvert` (the per-span
    /// `SimpleApprox::MaxError` plus its `Prec`/`Suivant` continuity terms).
    pub max_error: f64,
    /// `AverageError(3, 1)` (`cxx:1039-1070`): mean over the produced spans of
    /// `SimpleApprox::AverageError`.
    pub average_error: f64,
}

/// Result of a 2D `ApproxAFunction` run (`Perform` `cxx:663-956`): the
/// `Num1DSS = 0`, `Num2DSS = 1`, `Num3DSS = 0` layout used by
/// `Geom2dConvert_ApproxCurve` (`Geom2dConvert_ApproxCurve.cxx:137-183`).
pub struct ApproxAFunction2d {
    pub done: bool,
    pub has_result: bool,
    /// `AdvApprox_ApproxAFunction::Poles2d(1, Poles)` (`cxx:823-860`).
    pub poles: Vec<GpPnt2d>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub degree: i32,
    /// `MaxError(2, 1)` (`cxx:1004-1035`): maximum over the produced spans of
    /// the error bound accumulated by `PrepareConvert`.
    pub max_error: f64,
    /// `AverageError(2, 1)` (`cxx:1039-1070`): mean over the produced spans of
    /// `SimpleApprox::AverageError`.
    pub average_error: f64,
}

/// Result of `AdvApprox_ApproxAFunction` with `Num1DSS=2` (`Approx_SameParameter`).
pub struct ApproxAFunction1dPair {
    pub done: bool,
    pub has_result: bool,
    pub poles_u: Vec<f64>,
    pub poles_v: Vec<f64>,
    pub knots: Vec<f64>,
    pub mults: Vec<i32>,
    pub degree: i32,
}
// ---------------------------------------------------------------------------
// Submodules
// ---------------------------------------------------------------------------

mod function_1d_3d;
mod function_1d_pair;
mod function_2d;
mod function_3d;
mod function_3d_multi;
pub(super) mod prepare_convert;

pub use self::function_1d_3d::ApproxAFunction1d3d;
pub use self::function_3d_multi::ApproxAFunction3dMulti;
