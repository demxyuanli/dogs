//! Constants of `BRepBlend_SurfRstConstRad`
//! (OCCT `TKFillet/BRepBlend/BRepBlend_SurfRstConstRad.cxx`).

/// `#define Eps 1.e-15` (`cxx:35`): the jacobian magnitude below which `Section`
/// treats the surface as singular (`cxx:853`).
pub const EPS: f64 = 1.0e-15;

/// `NT` floor of `Decroch` (`cxx:566`): `1.e-7`.
pub const DECROCH_NT_EPS: f64 = 1.0e-7;

/// `dot` threshold of `Decroch` (`cxx:573`): `1.e-10`.
pub const DECROCH_DOT_EPS: f64 = 1.0e-10;

/// `Eps` passed to the `math_SVD::Solve` fallback of `IsSolution`
/// (`cxx:373`) and `Section` (`cxx:889`).
pub const SVD_EPS: f64 = 1.0e-6;

/// `math_Gauss` default `MinPivot` (`math_Gauss.hxx:49`), used by the
/// two-argument `math_Gauss Resol(gradsol)` in `IsSolution` (`cxx:359`).
pub const GAUSS_DEFAULT_EPS: f64 = 1.0e-20;

/// `Eps` passed to `math_Gauss(gradsol, 1.e-9)` in `Section`
/// (`cxx:876`).
pub const GAUSS_SECTION_EPS: f64 = 1.0e-9;
