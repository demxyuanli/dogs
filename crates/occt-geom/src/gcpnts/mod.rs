//! Curve discretization points. Port of the `GCPnts`/`CPnts` packages
//! (TKGeomBase): `CPnts_AbscissaPoint`, `GCPnts_AbscissaPoint`,
//! `GCPnts_UniformAbscissa`, `GCPnts_QuasiUniformAbscissa`, plus a parameter
//! accessor for the tangential-deflection sampler.
//!
//! Arc length is the faithful `CPnts_AbscissaPoint::Length` path: the speed
//! `|C'(u)|` is integrated with `math_GaussSingleIntegration` at OCCT's
//! per-type order (`CPnts_AbscissaPoint.cxx:57-77`) and, for the tolerance
//! overload, with the 13-iteration interval-doubling loop
//! (`math_GaussSingleIntegration.cxx:64-98`). OCCT tabulates the Gauss nodes and
//! weights (`math.cxx` `Point[]`/`Weight[]`, `GaussPointsMax() = 61`); the port
//! computes the same values with [`occt_math::gauss::gauss_legendre`] and clamps
//! the order to 61 the same way.
//!
//! The abscissa *inversion* (`GCPnts_AbscissaPoint::Parameter`) is the faithful
//! `CPnts_AbscissaPoint::Init`/`Perform` path
//! (`CPnts_AbscissaPoint.cxx:270-313`, `:374-432`) driven by
//! `math_FunctionRoot` on `CPnts_MyRootFunction` / `CPnts_MyGaussFunction`
//! (`CPnts_MyRootFunction.cxx:19-90`, `CPnts_MyGaussFunction.cxx:17-27`);
//! `math_FunctionRoot` itself runs `math_FunctionSetRoot`
//! (`math_FunctionRoot.cxx:73-119`).
//!
//! **`GCPnts_AbscissaPoint::Compute` is ported.** Its type dispatch
//! (`GCPnts_AbscissaPoint.cxx:26-65`, `:67-161`) is [`compute_type`], and all
//! three arms are reachable: `GCPnts_LengthParametrized` (`:87-90`, when
//! `NbIntervals(GeomAbs_CN) == 1` and the type is Line/Circle or a 2-pole
//! non-rational Bezier/BSpline), `GCPnts_Parametrized` (`:91-95`, the iterative
//! root search) and `GCPnts_AbsComposite` (`:96-158`, [`compute_abs_composite`]).
//!
//! **`GCPnts_AbscissaPoint::AdvCompute` (`:166-295`) is ported too**: the same
//! three arms with the tolerance-aware integration
//! (`CPnts_AbscissaPoint::Length(..., Tol)`), the `GCPnts_AbsComposite` arm
//! including its `anIndex == 0` preamble (`:202-229`) and the `aSign`-shifted,
//! first/last-parameter-clamped tail (`:273-293`) in
//! [`compute_abs_composite_adv`]. `theEPSILON` there is the caller's `theTol`
//! (`:486`, `:540`), not `Resolution(Precision::Confusion())`.
//!
//! T-96 (2026-09-30) checked the earlier claim that these arms "cannot be
//! expressed faithfully" because this port's `GeomTrimmedCurve` remaps a trim
//! onto `[0, 1]` while OCCT's `GeomAdaptor_Curve` unwraps to the basis. **The
//! claim is stale**: the `[0, 1]` convention is applied consistently to *every*
//! quantity the three arms read —
//! `d1`/`d2`/`d3`/`eval_dn` carry the `s = last - first` chain rule
//! (`trimmed.rs:45-79`), `parameter_intervals` remaps the basis knots by the same
//! `s` (`trimmed.rs:130-140`), and `resolution` divides by the same `s`
//! (`trimmed.rs:147-153`). That is a plain reparameterisation: the numerical
//! values the arms consume keep their meaning, so the arms are expressible and
//! are expressed. See `specs/_a3n00_gap_analysis.md` §9.295 / §9.296.
//!
//! `uniform_abscissa` / `quasi_uniform_abscissa` take an **interval count** and
//! return `n + 1` parameters. That is equivalent to OCCT's `NbPoints` overload,
//! not a shape difference: `GCPnts_UniformAbscissa::initialize` sets
//! `anAbscissa = aL / (theNbPoints - 1)` (`GCPnts_UniformAbscissa.cxx:510`), i.e.
//! it uses `NbPoints - 1` **intervals** for `NbPoints` points (endpoints
//! included). So `uniform_abscissa(c, n)` == OCCT with `NbPoints = n + 1`; only
//! the argument naming differs (intervals here, points there).
//!
//! `GCPnts_UniformDeflection` / `GCPnts_QuasiUniformDeflection` are **UNPORTED**
//! (see `occt-core/src/gcpnts.rs` module docs) — this module no longer exposes a
//! fake "uniform deflection" accessor.

mod abscissa;
mod length;
mod sampling;

#[allow(unused_imports)]
use self::length::{compute_abs_composite, compute_abs_composite_adv, compute_type};

pub use self::abscissa::{
    abscissa_point, abscissa_point_from_guess, abscissa_point_with_tolerance,
};
pub use self::length::{curve_length, curve_length_between, curve_length_range};
pub use self::sampling::{quasi_uniform_abscissa, tangential_deflection, uniform_abscissa};

#[cfg(test)]
mod tests;
