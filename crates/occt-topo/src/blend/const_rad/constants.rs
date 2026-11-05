//! Constants of `BlendFunc_ConstRad` and `BlendFunc_ConstRadInv`
//! (OCCT `TKFillet/BlendFunc/BlendFunc_ConstRad.cxx`,
//! `BlendFunc_ConstRadInv.cxx`).

/// `#define Eps 1.e-15` (`BlendFunc_ConstRad.cxx:41`,
/// `BlendFunc_ConstRadInv.cxx:23`): the normal magnitude below which a surface
/// is treated as singular.
pub const EPS: f64 = 1.0e-15;

/// `math_Gauss Resol(DEDX, 1.e-9)` of `IsSolution` (`cxx:864`) and of the two
/// derivative `Section`s (`cxx:1392`, `cxx:1708`).
pub const GAUSS_SOLUTION_EPS: f64 = 1.0e-9;

/// `Eps` of the `math_SVD::Solve` fallback (`cxx:887`, `cxx:1405`,
/// `cxx:1729`): `1.e-6`.
pub const SVD_EPS: f64 = 1.0e-6;

/// `pts1.IsEqual(pts2, 1.e-4)` (`cxx:1387`, `cxx:1704`): below this distance
/// the two extremities are considered merged and the derivative solve is
/// skipped.
pub const POINTS_EQUAL_TOL: f64 = 1.0e-4;

/// Initial value of the cache fields `tval` and `xval`
/// (`BlendFunc_ConstRad.cxx:70-72`): `-9.876e100`.
pub const CACHE_EMPTY: f64 = -9.876e100;
