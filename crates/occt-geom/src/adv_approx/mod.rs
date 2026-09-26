//! Adaptive approximation of a function by a B-spline.
//! Source: `AdvApprox_ApproxAFunction.cxx`, `AdvApprox_SimpleApprox.cxx`,
//! `AdvApprox_Cutting.hxx`, `AdvApprox_DichoCutting.cxx`,
//! `AdvApprox_PrefCutting.cxx`, `AdvApprox_PrefAndRec.cxx`.

mod cutting;
mod simple;
mod approx;

pub use approx::{ApproxAFunction1dPair, ApproxAFunction3d};
pub use cutting::{Cutting, DichoCutting, PrefAndRec, PrefCutting};
