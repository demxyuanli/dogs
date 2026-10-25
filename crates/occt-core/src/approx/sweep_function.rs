//! `Approx_SweepFunction` — the callbacks `Approx_SweepApproximation` drives
//! for a general sweep.
//!
//! Source: `TKGeomBase/Approx/Approx_SweepFunction.hxx:34-156` and
//! `Approx_SweepFunction.cxx:25-87` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The `NCollection_Array1` out-parameters become returned vectors; the
//!   caller knows their lengths from `SectionShape` / `Nb2dCurves`, exactly the
//!   OCCT contract.
//! - `Standard_NotImplemented` from the base `D1` / `D2` / `Resolution` /
//!   `BarycentreOfSurf` / `MaximalSection` / `GetMinimalWeight`
//!   (`cxx:24-87`) becomes `Err`.
//! - No `Send` / `Sync` bound: OCCT's `Approx_SweepFunction` is plain
//!   inheritance, and `BRepBlend_AppFuncRoot` (`BRepBlend_AppFuncRoot.cxx:26-46`)
//!   aliases the caller's `Blend_AppFunction&` through an untyped pointer,
//!   which is a single-threaded `Rc<RefCell<..>>` in the port. Nothing in the
//!   crate moves an `Approx_SweepFunction` across a thread boundary.

use crate::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

/// `D1` result of a sweep function (`Approx_SweepFunction.hxx:51-59`).
#[derive(Debug, Clone, Default)]
pub struct SweepD1 {
    pub poles: Vec<GpPnt>,
    pub dpoles: Vec<GpVec>,
    pub poles2d: Vec<GpPnt2d>,
    pub dpoles2d: Vec<GpVec2d>,
    pub weights: Vec<f64>,
    pub dweights: Vec<f64>,
}

/// `D2` result of a sweep function (`Approx_SweepFunction.hxx:64-76`).
#[derive(Debug, Clone, Default)]
pub struct SweepD2 {
    pub poles: Vec<GpPnt>,
    pub dpoles: Vec<GpVec>,
    pub d2poles: Vec<GpVec>,
    pub poles2d: Vec<GpPnt2d>,
    pub dpoles2d: Vec<GpVec2d>,
    pub d2poles2d: Vec<GpVec2d>,
    pub weights: Vec<f64>,
    pub dweights: Vec<f64>,
    pub d2weights: Vec<f64>,
}

/// `Approx_SweepFunction` (`Approx_SweepFunction.hxx:34-156`).
pub trait ApproxSweepFunction {
    /// `Approx_SweepFunction::D0` (`hxx:40-46`, pure virtual).
    fn d0(
        &self,
        param: f64,
        first: f64,
        last: f64,
    ) -> Result<(Vec<GpPnt>, Vec<GpPnt2d>, Vec<f64>), String>;

    /// `Approx_SweepFunction::D1` (`cxx:24-45`): `Standard_NotImplemented`.
    fn d1(&self, _param: f64, _first: f64, _last: f64) -> Result<SweepD1, String> {
        Err("Standard_NotImplemented: Approx_SweepFunction::D1".into())
    }

    /// `Approx_SweepFunction::D2` (`cxx:47-72`): `Standard_NotImplemented`.
    fn d2(&self, _param: f64, _first: f64, _last: f64) -> Result<SweepD2, String> {
        Err("Standard_NotImplemented: Approx_SweepFunction::D2".into())
    }

    /// `Approx_SweepFunction::Nb2dCurves` (`hxx:79`, pure virtual).
    fn nb_2d_curves(&self) -> i32;

    /// `Approx_SweepFunction::SectionShape` (`hxx:83`, pure virtual):
    /// `(nb_poles, nb_knots, degree)`.
    fn section_shape(&self) -> (i32, i32, i32);

    /// `Approx_SweepFunction::Knots` (`hxx:86`, pure virtual).
    fn knots(&self) -> Vec<f64>;

    /// `Approx_SweepFunction::Mults` (`hxx:89`, pure virtual).
    fn mults(&self) -> Vec<i32>;

    /// `Approx_SweepFunction::IsRational` (`hxx:89`, pure virtual).
    fn is_rational(&self) -> bool;

    /// `Approx_SweepFunction::NbIntervals` (`hxx:95`, pure virtual).
    fn nb_intervals(&self, s: u8) -> i32;

    /// `Approx_SweepFunction::Intervals` (`hxx:102-103`, pure virtual).
    fn intervals(&self, s: u8) -> Vec<f64>;

    /// `Approx_SweepFunction::SetInterval` (`hxx:109`, pure virtual).
    fn set_interval(&self, first: f64, last: f64);

    /// `Approx_SweepFunction::Resolution` (`cxx:74-79`):
    /// `Standard_NotImplemented`.
    fn resolution(&self, _index: i32, _tol: f64) -> Result<(f64, f64), String> {
        Err("Standard_NotImplemented: Approx_SweepFunction::Resolution".into())
    }

    /// `Approx_SweepFunction::GetTolerance` (`hxx:125-128`, pure virtual).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64) -> Vec<f64>;

    /// `Approx_SweepFunction::SetTolerance` (`hxx:129`, pure virtual).
    fn set_tolerance(&self, tol3d: f64, tol2d: f64);

    /// `Approx_SweepFunction::BarycentreOfSurf` (`cxx:81-84`):
    /// `Standard_NotImplemented`.
    fn barycentre_of_surf(&self) -> Result<GpPnt, String> {
        Err("Standard_NotImplemented: Approx_SweepFunction::BarycentreOfSurf".into())
    }

    /// `Approx_SweepFunction::MaximalSection` (`cxx:86-89`):
    /// `Standard_NotImplemented`.
    fn maximal_section(&self) -> Result<f64, String> {
        Err("Standard_NotImplemented: Approx_SweepFunction::MaximalSection()".into())
    }

    /// `Approx_SweepFunction::GetMinimalWeight` (`cxx:91-94`):
    /// `Standard_NotImplemented`.
    fn get_minimal_weight(&self) -> Result<Vec<f64>, String> {
        Err("Standard_NotImplemented: Approx_SweepFunction::GetMinimalWeight".into())
    }
}
