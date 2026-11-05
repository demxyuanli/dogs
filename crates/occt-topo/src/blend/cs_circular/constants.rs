//! Constants of `BlendFunc_CSCircular`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_CSCircular.cxx:39`
//! (`#define Eps 1.e-15`), `math_Gauss.hxx:49` (the default `MinPivot`)
//! (OCCT 8.0.0).

/// `Eps` (`cxx:39`): below it the projected normal is not divided by
/// `|nplan x ns|` (`if (norm < Eps) norm = 1.;`).
pub(super) const EPS: f64 = 1.0e-15;

/// `math_Gauss`'s default `MinPivot` (`math_Gauss.hxx:49`): `IsSolution`,
/// `GetSection` and both `Section(P, ...)` build the solver without a pivot
/// argument.
pub(super) const GAUSS_MIN_PIVOT: f64 = 1.0e-20;
