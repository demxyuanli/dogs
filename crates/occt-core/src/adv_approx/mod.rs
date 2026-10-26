//! Adaptive approximation of a function by a B-spline.
//! Source: `AdvApprox_ApproxAFunction.cxx`, `AdvApprox_SimpleApprox.cxx`,
//! `AdvApprox_Cutting.hxx`, `AdvApprox_DichoCutting.cxx`,
//! `AdvApprox_PrefCutting.cxx`, `AdvApprox_PrefAndRec.cxx`.

mod cutting;
mod simple;
mod approx;
mod approx_1d_pair_3d;
mod mixed_dimension;

pub use approx::{
    ApproxAFunction1d3d, ApproxAFunction1dPair, ApproxAFunction2d, ApproxAFunction3d,
    ApproxAFunction3dMulti,
};
pub use approx_1d_pair_3d::ApproxAFunction1dPair3d;
pub use mixed_dimension::ApproxAFunction;
pub use cutting::{Cutting, DichoCutting, PrefAndRec, PrefCutting};
