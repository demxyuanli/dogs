//! `AdvApp2Var_ApproxAFunc2Var` (`AdvApp2Var_ApproxAFunc2Var.cxx`).
//!
//! `Init` builds `AdvApp2Var_Context` through `mma2jmx_`, `mma2roo_`, and
//! `mmapptt_`. A finished square is stored. `Compute3DErrors` runs, then
//! `ConvertBS` returns the B-spline.

mod border_error;
mod compute_patches;
mod constraints;
mod convert_bs;
mod context;
mod degree_error;
mod degree_reduce;
mod discretise;
mod eval_surface;
mod framework;
mod grid_polynomial;
mod gauss;
mod iso_fnc;
mod jacobi_canon;
mod jacobi_max;
mod make_approx;
mod network;
mod legendre_roots;
mod node;
mod param_normalize;
mod square_disc;
mod table_expand;
mod xmaxj4;

pub use context::{context_until_roots, App2VarContext};
pub use iso_fnc::approximate_until_cdi;
