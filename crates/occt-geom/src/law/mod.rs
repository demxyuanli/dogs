//! 1D evolution laws — the OCCT `Law` package.
//!
//! Source: `TKGeomAlgo/Law` (OCCT 8.0.0):
//! `Law.hxx`/`.cxx`, `Law_Function.hxx`/`.cxx`, `Law_Constant.hxx`/`.cxx`,
//! `Law_Linear.hxx`/`.cxx`, `Law_S.hxx`/`.cxx`, `Law_BSpFunc.hxx`/`.cxx`,
//! `Law_BSpline.hxx`/`.cxx`, `Law_BSplineKnotSplitting.hxx`/`.cxx`,
//! `Law_Interpolate.hxx`/`.cxx`, `Law_Interpol.hxx`/`.cxx`,
//! `Law_Composite.hxx`/`.cxx`.
//!
//! A law maps a parameter to a scalar value $\lambda(t)$ together with its
//! derivatives D0 / D1 / D2, and reports the parametric intervals on which it
//! reaches a requested continuity.
//!
//! Module layout:
//! ```text
//! function   Law_Function (the root trait) and Law_Constant / Law_Linear
//! bspline    Law_BSpline (the 1D B-spline) and Law_BSpFunc
//! knot_split Law_BSplineKnotSplitting
//! interpolate Law_Interpolate and Law_Interpol
//! composite  Law_Composite
//! handle     the shared `handle(Law_Composite)` adapter
//! services   the Law static helpers (MixBnd, MixTgt, Reparametrize, Scale)
//! ```

pub mod bspline;
pub mod composite;
pub mod function;
pub mod handle;
pub mod interpolate;
pub mod knot_split;
pub mod services;

pub use bspline::{LawBSpFunc, LawBSpline};
pub use composite::LawComposite;
pub use function::{LawConstant, LawFunction, LawLinear, LawS};
pub use handle::LawHandle;
pub use interpolate::{LawInterpol, LawInterpolate};
pub use knot_split::LawBSplineKnotSplitting;
pub use services::{
    law_mix_bnd, law_mix_bnd_poles, law_mix_tgt, law_reparametrize, law_scale, law_scale_cub,
};
