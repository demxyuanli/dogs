//! Adaptive approximation of a function by a B-spline.
//! Source: `AdvApprox_ApproxAFunction.cxx`, `AdvApprox_SimpleApprox.cxx`.

mod simple;
mod approx;

pub use approx::{ApproxAFunction1dPair, ApproxAFunction3d};
