//! `GeomFill_TrihedronLaw` — the abstract law that lays a moving frame on the
//! sweep spine, plus the `myCurve` / `myTrimmed` bookkeeping the OCCT base
//! class owns.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_TrihedronLaw.hxx:29-113` and
//! `GeomFill_TrihedronLaw.cxx:24-76` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `occ::handle<Adaptor3d_Curve>` is `Arc<dyn Curve>`.
//! - OCCT shares one law object behind several handles and calls the non-const
//!   `SetCurve` / `SetInterval` through any of them (`GeomFill_Sweep` keeps a
//!   `myTlaw` handle that `GeomFill_Pipe` also owns). The port therefore keeps
//!   the mutated state behind a `Mutex` and takes `&self` in those methods,
//!   which reproduces the handle semantics.
//! - `myTrimmed = myCurve->Trim(First, Last, 0)` (`cxx:57-60`) keeps the
//!   parameter values of `myCurve` (an `Adaptor3d_Curve` trim does not
//!   reparametrize), so [`TrihedronLawCurve`] stores the trimmed range and
//!   keeps evaluating through `myCurve`, exactly as the OCCT laws do.
//! - `Standard_NotImplemented` from the base `D1` / `D2` (`cxx:37-55`) becomes
//!   `Err("Standard_NotImplemented: ...")`.

use std::sync::{Arc, Mutex};

use occt_core::gp::GpVec;

use crate::curve::Curve;
use crate::geom_fill::pipe_error::GeomFillPipeError;

/// The `myCurve` / `myTrimmed` pair of `GeomFill_TrihedronLaw`
/// (`hxx:111-112`).
pub struct TrihedronLawCurve {
    inner: Mutex<CurveState>,
}

#[derive(Clone)]
struct CurveState {
    curve: Option<Arc<dyn Curve>>,
    trimmed: Option<(f64, f64)>,
}

impl Default for TrihedronLawCurve {
    fn default() -> Self {
        Self {
            inner: Mutex::new(CurveState {
                curve: None,
                trimmed: None,
            }),
        }
    }
}

impl Clone for TrihedronLawCurve {
    fn clone(&self) -> Self {
        let guard = self.inner.lock().unwrap_or_else(|e| e.into_inner());
        Self {
            inner: Mutex::new(guard.clone()),
        }
    }
}

impl TrihedronLawCurve {
    /// New empty law data (OCCT leaves both handles null).
    pub fn new() -> Self {
        Self::default()
    }

    /// `GeomFill_TrihedronLaw::SetCurve` (`cxx:25-30`): `myCurve = C`,
    /// `myTrimmed = C`, then `true`.
    pub fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        let mut state = self.inner.lock().unwrap_or_else(|e| e.into_inner());
        state.trimmed = Some((curve.first_parameter(), curve.last_parameter()));
        state.curve = Some(curve);
        true
    }

    /// `GeomFill_TrihedronLaw::SetInterval` (`cxx:57-60`):
    /// `myTrimmed = myCurve->Trim(First, Last, 0)`.
    pub fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        let mut state = self.inner.lock().unwrap_or_else(|e| e.into_inner());
        if state.curve.is_none() {
            return Err("Standard_NullObject: GeomFill_TrihedronLaw::SetInterval".into());
        }
        state.trimmed = Some((first, last));
        Ok(())
    }

    /// `GeomFill_TrihedronLaw::GetInterval` (`cxx:62-67`).
    pub fn interval(&self) -> (f64, f64) {
        let state = self.inner.lock().unwrap_or_else(|e| e.into_inner());
        state.trimmed.unwrap_or((0.0, 0.0))
    }

    /// `myCurve` (`hxx:111`).
    pub fn curve(&self) -> Option<Arc<dyn Curve>> {
        self.inner
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .curve
            .clone()
    }

    /// `myCurve` for the laws that dereference it unconditionally.
    pub fn curve_or_err(&self, who: &str) -> Result<Arc<dyn Curve>, String> {
        self.curve()
            .ok_or_else(|| format!("Standard_NullObject: {who}"))
    }
}

/// `D1` result: the frame and its first derivative
/// (`GeomFill_TrihedronLaw.hxx:51-57`).
#[derive(Debug, Clone, Copy)]
pub struct TrihedronD1 {
    pub tangent: GpVec,
    pub dtangent: GpVec,
    pub normal: GpVec,
    pub dnormal: GpVec,
    pub binormal: GpVec,
    pub dbinormal: GpVec,
}

/// `D2` result: the frame and its first two derivatives
/// (`GeomFill_TrihedronLaw.hxx:62-72`).
#[derive(Debug, Clone, Copy)]
pub struct TrihedronD2 {
    pub tangent: GpVec,
    pub dtangent: GpVec,
    pub d2tangent: GpVec,
    pub normal: GpVec,
    pub dnormal: GpVec,
    pub d2normal: GpVec,
    pub binormal: GpVec,
    pub dbinormal: GpVec,
    pub d2binormal: GpVec,
}

/// `GeomFill_TrihedronLaw` (`GeomFill_TrihedronLaw.hxx:29-113`).
pub trait TrihedronLaw: Send + Sync {
    /// `GeomFill_TrihedronLaw::SetCurve` (`cxx:25-30`).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool;

    /// `GeomFill_TrihedronLaw::Copy` (`hxx:36`, pure virtual).
    fn copy(&self) -> Arc<dyn TrihedronLaw>;

    /// `GeomFill_TrihedronLaw::ErrorStatus` (`cxx:34-37`): `PipeOk`.
    fn error_status(&self) -> GeomFillPipeError {
        GeomFillPipeError::PipeOk
    }

    /// `GeomFill_TrihedronLaw::D0` (`hxx:42-46`, pure virtual).
    fn d0(&self, param: f64) -> Result<(GpVec, GpVec, GpVec), String>;

    /// `GeomFill_TrihedronLaw::D1` (`cxx:37-44`): `Standard_NotImplemented`.
    fn d1(&self, _param: f64) -> Result<TrihedronD1, String> {
        Err("Standard_NotImplemented: GeomFill_TrihedronLaw::D1".into())
    }

    /// `GeomFill_TrihedronLaw::D2` (`cxx:46-55`): `Standard_NotImplemented`.
    fn d2(&self, _param: f64) -> Result<TrihedronD2, String> {
        Err("Standard_NotImplemented: GeomFill_TrihedronLaw::D2".into())
    }

    /// `GeomFill_TrihedronLaw::NbIntervals` (`hxx:77`, pure virtual).
    fn nb_intervals(&self, s: u8) -> i32;

    /// `GeomFill_TrihedronLaw::Intervals` (`hxx:83-85`, pure virtual): the
    /// parameters bounding the intervals of continuity `s`.
    fn intervals(&self, s: u8) -> Vec<f64>;

    /// `GeomFill_TrihedronLaw::SetInterval` (`cxx:57-60`).
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String>;

    /// `GeomFill_TrihedronLaw::GetInterval` (`cxx:62-67`).
    fn get_interval(&self) -> (f64, f64);

    /// `GeomFill_TrihedronLaw::GetAverageLaw` (`hxx:97-100`, pure virtual).
    fn get_average_law(&self) -> (GpVec, GpVec, GpVec);

    /// `GeomFill_TrihedronLaw::IsConstant` (`cxx:69-72`): false.
    fn is_constant(&self) -> bool {
        false
    }

    /// `GeomFill_TrihedronLaw::IsOnlyBy3dCurve` (`cxx:74-77`): false.
    fn is_only_by_3d_curve(&self) -> bool {
        false
    }
}
