//! Constants of `BlendFunc_CSConstRad`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_CSConstRad.cxx:37-989`
//! (OCCT 8.0.0).

/// `math_Gauss`'s default `MinPivot` (`math_Gauss.hxx:49`), the pivot threshold
/// of the linear system of `IsSolution` (`cxx:179`) and `GetSection`
/// (`cxx:588`) / `Section` (`cxx:870`).
pub(super) const GAUSS_MIN_PIVOT: f64 = 1.0e-20;
