//! Constants of `BRepBlend_RstRstEvolRad`
//! (OCCT `TKFillet/BRepBlend/BRepBlend_RstRstEvolRad.cxx`).

/// `#define Eps 1.e-15` (`cxx:40`). OCCT defines the macro but never uses it in
/// this translation unit; the port keeps it so the source line is not lost.
pub const EPS: f64 = 1.0e-15;

/// `math_Gauss` default `MinPivot` (`math_Gauss.hxx:49`), used by the
/// two-argument `math_Gauss Resol(gradsol)` in `IsSolution` (`cxx:288`).
pub const GAUSS_DEFAULT_EPS: f64 = 1.0e-20;

/// `Eps` passed to `math_Gauss(gradsol, 1.e-9)` in `Section` (`cxx:886`).
pub const GAUSS_SECTION_EPS: f64 = 1.0e-9;

/// `Eps` passed to the `math_SVD::Solve` fallback of `IsSolution` (`cxx:301`)
/// and `Section` (`cxx:900`).
pub const SVD_EPS: f64 = 1.0e-6;

/// `1.E-07` of `CenterCircleRst1Rst2` (`cxx:592`, `cxx:599`): the dead band on
/// `ray^2 - |P1P2|^2 / 4` below which the centre collapses onto the chord
/// midpoint, and below `-1e-7` there is no circle at all. The same bound
/// guards the derivative branch of `Section` (`cxx:935`).
pub const DIST_EPS: f64 = 1.0e-7;

/// `dot` threshold of `Decroch` (`cxx:527`, `cxx:531`, `cxx:541`): `1.e-10`.
pub const DECROCH_DOT_EPS: f64 = 1.0e-10;

/// `Precision::PConfusion() * 0.99` of `FusionneIntervalles` (`cxx:58`).
pub const FUSIONNE_EPS_FACTOR: f64 = 0.99;
