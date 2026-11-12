//! `AdvApp2Var_ApproxAFunc2Var` (`AdvApp2Var_ApproxAFunc2Var.cxx`).
//!
//! `Init` builds `AdvApp2Var_Context` through `mma2jmx_`, `mma2roo_`, and
//! `mmapptt_`. A finished square is stored. `Compute3DErrors` runs, then
//! `ConvertBS` returns the B-spline.

use crate::bspline_surface::GeomBSplineSurface;

mod border_error;
mod compute_patches;
mod constraints;
mod context;
mod convert_bs;
mod criterion;
mod degree_error;
mod degree_reduce;
mod discretise;
mod eval_surface;
mod framework;
mod gauss;
mod grid_polynomial;
mod iso_fnc;
mod jacobi_canon;
mod jacobi_max;
mod legendre_roots;
mod make_approx;
mod network;
mod node;
mod param_normalize;
mod square_disc;
mod subspace_error;
mod table_expand;
mod xmaxj4;

pub use context::{context_until_roots, App2VarContext};
pub use criterion::{
    patch_coefficients, Criterion, CriterionBase, CriterionPatch, CriterionRepartition,
    CriterionType,
};
pub use iso_fnc::{approximate_until_cdi, approximate_until_cdi_with_criterion};

/// Result of `AdvApp2Var_ApproxAFunc2Var` (`AdvApp2Var_ApproxAFunc2Var.hxx:86-192`).
///
/// `HasResult()` (`:148`) is `true` whenever this value exists: `Perform`
/// only returns after `ConvertBS` has stored a surface on every subspace
/// (`cxx:330-990`), and the two paths that leave `myHasResult` false
/// (`cxx:379`, `:446`, `:469`, `:503`, `:601`, `:624`, `:712`, `:807`)
/// throw `Standard_ConstructionError` instead of returning.
pub struct ApproxAFunc2VarResult {
    /// `Surface(1)` (`cxx:991-996`).
    pub surface: GeomBSplineSurface,
    /// All `Surface(Index)` results in 3D subspace order.
    pub surfaces: Vec<GeomBSplineSurface>,
    /// `IsDone()` (`cxx:998-1001`): `true` when every square met its
    /// tolerance and `Convert_GridPolynomialToPoles` succeeded.
    pub is_done: bool,
    /// `MaxError(3, 1)` (`cxx:1098-1106`).
    pub max_error: f64,
    /// `AverageError(3, 1)` (`cxx:1110-1118`).
    pub average_error: f64,
    /// `UFrontError(3, 1)` (`cxx:1122-1130`).
    pub u_front_error: f64,
    /// `VFrontError(3, 1)` (`cxx:1134-1142`).
    pub v_front_error: f64,
    /// `CritError(3, 1)` (`cxx:1146-1152`): `myCriterionError`, the maximum
    /// `CritValue` over the surviving squares (`ComputeCritError`,
    /// `cxx:868-887`). Stays `0` when no criterion ran, because
    /// `myNumSubSpaces[2] == 0` in that case.
    pub criterion_error: f64,
}
