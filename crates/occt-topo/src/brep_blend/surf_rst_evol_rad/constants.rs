//! Constants of `BRepBlend_SurfRstEvolRad`
//! (OCCT `TKFillet/BRepBlend/BRepBlend_SurfRstEvolRad.cxx`).

/// `#define Eps 1.e-15` (`cxx:37`): the jacobian magnitude below which
/// `Section(P, ...)` replaces `norm` by one (`cxx:894-899`).
pub const EPS: f64 = 1.0e-15;

/// `math_Gauss` default `MinPivot` (`math_Gauss.hxx:49`), used by the
/// two-argument calls `math_Gauss Resol(gradsol)` in `IsSolution`
/// (`cxx:420`) and `Section` (`cxx:917`).
pub const GAUSS_DEFAULT_EPS: f64 = 1.0e-20;

/// `Eps` passed to the `math_SVD::Solve` fallback of `IsSolution`
/// (`cxx:434`) and `Section` (`cxx:930`).
pub const SVD_EPS: f64 = 1.0e-6;

/// `NT` floor of `Decroch` (`cxx:634`): `1.e-7`.
pub const DECROCH_NT_EPS: f64 = 1.0e-7;

/// `dot` threshold of `Decroch` (`cxx:639`): `1.e-10`.
pub const DECROCH_DOT_EPS: f64 = 1.0e-10;
