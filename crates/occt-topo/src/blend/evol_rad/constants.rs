//! Numeric literals of `BlendFunc_EvolRad` / `BlendFunc_EvolRadInv`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_EvolRad.cxx:37`
//! (`#define Eps 1.e-15`), `:919` / `:924` / `:1155` / `:1161` (the `Eps`
//! guards), `:1486` (`math_Gauss Resol(DEDX, 1.e-9)`), `:1497`
//! (`Poles/Points ... IsEqual(pts2, 1.e-4)`), `:1507` (`math_SVD ... 1.e-6`),
//! `:1818` (`math_Gauss Resol(DEDX, 1.e-9)`), `:1834` / `:1841` (`math_SVD ...
//! 1.e-6`), `:1826` / `:1834` note the `1.e-14` pivot of `IsSolution`
//! (`cxx:900`) (OCCT 8.0.0).

/// `#define Eps 1.e-15` (`cxx:37`): the degeneracy guard of the section plane
/// terms.
pub const EPS: f64 = 1e-15;

/// `math_Gauss Resol(DEDX, 1.e-9)` in both `Section` overloads
/// (`cxx:1486`, `:1818`).
pub const GAUSS_SOLUTION_EPS: f64 = 1e-9;

/// `double maxpiv = 1.e-14; math_Gauss Resol(DEDX, maxpiv)` of `IsSolution`
/// (`cxx:899-900`). `BlendFunc_ConstRad` uses `1.e-9` there; `EvolRad` does
/// not, and it has no `math_SVD` fallback (`cxx:934-959`).
pub const GAUSS_IS_SOLUTION_MAXPIV: f64 = 1e-14;

/// `pts1.IsEqual(pts2, 1.e-4)` (`cxx:1592`): when the two extremities
/// coincide the derivative system is solved by the SVD fallback instead.
pub const POINTS_EQUAL_TOL: f64 = 1e-4;

/// `math_SVD SingRS(DEDX)` + `SingRS.Solve(..., 1.e-6)` (`cxx:1507`,
/// `:1841`).
pub const SVD_EPS: f64 = 1e-6;

/// `tval = -9.876e100; xval.Init(-9.876e100)` (`cxx:129-130`): the cache
/// sentinel of `ComputeValues`.
pub const CACHE_EMPTY: f64 = -9.876e100;
