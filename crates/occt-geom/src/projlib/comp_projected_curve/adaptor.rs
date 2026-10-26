//! `ProjLib_CompProjectedCurve` used as an `Adaptor2d_Curve2d`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_CompProjectedCurve.cxx`
//! `Perform` `cxx:1370-1374` builds `HPCur` from `ShallowCopy()` of the
//! projection itself, which is possible because the class derives from
//! `Adaptor2d_Curve2d` (`hxx:36`). The Rust class is not a curve, so this
//! wrapper supplies that interface.
//!
//! The parameter space is the **3D curve's**: `D0(U, P)` is
//! `ProjLib_CompProjectedCurve::D0`, `FirstParameter`/`LastParameter` are the
//! projected curve's (`cxx:1766-1776`), and `NbIntervals`/`Intervals` go
//! through `BuildIntervals` (`cxx:1811-2093`), whose table is what
//! `Adaptor2d_Curve2d` reports (`cxx:1822-1828`).
//!
//! Translation notes:
//! - `Standard_DomainError` from `D0` (`cxx:1578`) is an `Err` in the Rust
//!   evaluator; the trait cannot report it, so it is raised as the same panic
//!   OCCT's uncaught exception would produce (`expect`), matching the rest of
//!   this port's `PrjResolve` call sites.
//! - `Transform` / `Reverse` do not exist on `Adaptor2d_Curve2d`; both are
//!   no-ops, as on the other adaptor wrappers of this port.

use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};
use occt_geom2d::curve::Curve2d;

use super::CompProjectedCurve;

/// `Adaptor2d_Curve2d` view of a `ProjLib_CompProjectedCurve`.
pub struct ProjectedCurve2d {
    curve: CompProjectedCurve,
}

impl ProjectedCurve2d {
    pub fn new(curve: CompProjectedCurve) -> Self {
        Self { curve }
    }
}

impl Curve2d for ProjectedCurve2d {
    /// `D0(U, P)` (`cxx:1563-1704`), the `Adaptor2d_Curve2d::Value`.
    fn d0(&self, u: f64) -> GpPnt2d {
        self.curve.d0(u).expect("ProjLib_CompProjectedCurve::D0")
    }

    /// `D1(U, P, V)` (`cxx:1708-1717`).
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        self.curve.d1(u).expect("ProjLib_CompProjectedCurve::D1")
    }

    /// `D2(U, P, V1, V2)` (`cxx:1720-1727`).
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        self.curve.d2(u).expect("ProjLib_CompProjectedCurve::D2")
    }

    /// `FirstParameter()` (`cxx:1766-1769`).
    fn first_parameter(&self) -> f64 {
        self.curve.first_parameter()
    }

    /// `LastParameter()` (`cxx:1773-1776`).
    fn last_parameter(&self) -> f64 {
        self.curve.last_parameter()
    }

    /// `Continuity()` (`cxx:1781-1796`).
    fn continuity(&self) -> u8 {
        self.curve.continuity()
    }

    /// `Intervals(T, S)` (`cxx:1822-1828`): the fused partition of
    /// `BuildIntervals`. `BuildIntervals` rejects an unknown `GeomAbs_Shape`
    /// with `Standard_OutOfRange` (`cxx:1852`); the trait cannot report it, so
    /// the single whole span is returned for that (unreachable from the
    /// approximation, which passes `C1..C3`) case.
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        self.curve
            .build_intervals(continuity)
            .unwrap_or_else(|_| vec![self.curve.first_parameter(), self.curve.last_parameter()])
    }

    /// `Transform` / `Reverse` are not part of `Adaptor2d_Curve2d`.
    fn transform(&mut self, _t: &GpTrsf2d) {}

    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(Self::new(self.curve.clone()))
    }
}
