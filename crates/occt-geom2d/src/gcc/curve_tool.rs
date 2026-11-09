//! `Geom2dGcc_CurveTool` (`Geom2dGcc_CurveTool.{hxx:24-62,cxx:21-71}`): the
//! `TheCurveTool` of the `Geom2dGcc_*` iterative solvers.
//!
//! Every entry is a single delegation to the `Geom2dAdaptor_Curve` argument,
//! so the port is a set of free functions over [`Geom2dAdaptorCurve`] rather
//! than a type: there is no state and no overload to dispatch on.
//!
//! Translation notes:
//! - `D1` / `D2` / `D3` (`cxx:47-71`) forward to `Geom2dAdaptor_Curve::D1` /
//!   `D2` / `D3`, which are `D0` plus `DN(1)` / `DN(1,2)` / `DN(1,2,3)`; the
//!   port composes them from [`Geom2dAdaptorCurve::value`] and
//!   [`Geom2dAdaptorCurve::dn`] and returns `None` when the adaptor holds no
//!   curve, which is the state OCCT answers with a null handle dereference,
//! - `NbSamples` (`cxx:26-29`) ignores its argument and returns the constant
//!   `20`, so the port keeps it as a constant plus a function with the OCCT
//!   signature,
//! - `EpsX(C, Tol)` (`cxx:22-25`) is `C.Resolution(Tol)`, i.e.
//!   [`Geom2dAdaptorCurve::resolution`].
//!
//! UNPORTED: none.

use occt_core::gp::{GpPnt2d, GpVec2d};

use crate::adaptor_curve::Geom2dAdaptorCurve;

/// `NbSamples`'s return value (`cxx:28`): the argument is unused.
pub const NB_SAMPLES: i32 = 20;

/// `FirstParameter(C)` (`cxx:37-40`).
pub fn first_parameter(c: &Geom2dAdaptorCurve) -> f64 {
    c.first_parameter()
}

/// `LastParameter(C)` (`cxx:42-45`).
pub fn last_parameter(c: &Geom2dAdaptorCurve) -> f64 {
    c.last_parameter()
}

/// `EpsX(C, Tol)` (`cxx:22-25`): `C.Resolution(Tol)`.
pub fn eps_x(c: &Geom2dAdaptorCurve, tol: f64) -> f64 {
    c.resolution(tol)
}

/// `NbSamples(C)` (`cxx:26-29`): constant `20`.
pub fn nb_samples(_c: &Geom2dAdaptorCurve) -> i32 {
    NB_SAMPLES
}

/// `Value(C, X)` (`cxx:31-34`).
pub fn value(c: &Geom2dAdaptorCurve, x: f64) -> Option<GpPnt2d> {
    c.value(x)
}

/// `D1(C, U, P, T)` (`cxx:47-50`).
pub fn d1(c: &Geom2dAdaptorCurve, u: f64) -> Option<(GpPnt2d, GpVec2d)> {
    Some((c.value(u)?, c.dn(u, 1)?))
}

/// `D2(C, U, P, T, N)` (`cxx:52-60`).
pub fn d2(c: &Geom2dAdaptorCurve, u: f64) -> Option<(GpPnt2d, GpVec2d, GpVec2d)> {
    Some((c.value(u)?, c.dn(u, 1)?, c.dn(u, 2)?))
}

/// `D3(C, U, P, T, N, dN)` (`cxx:62-71`).
pub fn d3(c: &Geom2dAdaptorCurve, u: f64) -> Option<(GpPnt2d, GpVec2d, GpVec2d, GpVec2d)> {
    Some((c.value(u)?, c.dn(u, 1)?, c.dn(u, 2)?, c.dn(u, 3)?))
}
