//! `GeomFill_LocationLaw` — the abstract law that places the section frame on
//! the sweep spine.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_LocationLaw.hxx:33-164` and
//! `GeomFill_LocationLaw.cxx:24-114` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `occ::handle<Adaptor3d_Curve>` is `Arc<dyn Curve>`.
//! - `NCollection_Array1<gp_Pnt2d>` / `<gp_Vec2d>` become slices; the caller
//!   sizes them from [`LocationLaw::nb_2d_curves`], exactly the OCCT contract.
//! - The base implementations that throw `Standard_NotImplemented`
//!   (`cxx:28-55`, `cxx:92-96`, `cxx:112-115`) become `Err`.
//! - `SetTolerance` (`cxx:98-101`) and the `IsTranslation` / `IsRotation`
//!   defaults (`cxx:103-110`) keep the OCCT "do nothing / return false"
//!   behaviour.

use std::sync::Arc;

use occt_core::gp::{GpMat, GpPnt, GpPnt2d, GpVec, GpVec2d};

use crate::curve::Curve;
use crate::geom_fill::pipe_error::GeomFillPipeError;

/// `D1` result of a location law (`GeomFill_LocationLaw.hxx:61-68`).
#[derive(Debug, Clone)]
pub struct LocationD1 {
    pub m: GpMat,
    pub v: GpVec,
    pub dm: GpMat,
    pub dv: GpVec,
    pub poles2d: Vec<GpPnt2d>,
    pub dpoles2d: Vec<GpVec2d>,
}

/// `D2` result of a location law (`GeomFill_LocationLaw.hxx:72-82`).
#[derive(Debug, Clone)]
pub struct LocationD2 {
    pub m: GpMat,
    pub v: GpVec,
    pub dm: GpMat,
    pub dv: GpVec,
    pub d2m: GpMat,
    pub d2v: GpVec,
    pub poles2d: Vec<GpPnt2d>,
    pub dpoles2d: Vec<GpVec2d>,
    pub d2poles2d: Vec<GpVec2d>,
}

/// `GeomFill_LocationLaw` (`GeomFill_LocationLaw.hxx:33-164`).
pub trait LocationLaw: Send + Sync {
    /// `GeomFill_LocationLaw::SetCurve` (`hxx:39`, pure virtual).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool;

    /// `GeomFill_LocationLaw::GetCurve` (`hxx:41`, pure virtual). OCCT returns
    /// the stored handle; the port clones the `Arc`.
    fn get_curve(&self) -> Option<Arc<dyn Curve>>;

    /// `GeomFill_LocationLaw::SetTrsf` (`hxx:46`): `M(t)` becomes
    /// `Mat * M(t)`.
    fn set_trsf(&self, transfo: &GpMat);

    /// `GeomFill_LocationLaw::Copy` (`hxx:48`, pure virtual).
    fn copy(&self) -> Arc<dyn LocationLaw>;

    /// `GeomFill_LocationLaw::D0(Param, M, V)` (`hxx:51`, pure virtual).
    fn d0(&self, param: f64) -> Result<(GpMat, GpVec), String>;

    /// `GeomFill_LocationLaw::D0(Param, M, V, Poles2d)` (`hxx:54-57`, pure
    /// virtual).
    fn d0_with_poles(&self, param: f64) -> Result<(GpMat, GpVec, Vec<GpPnt2d>), String>;

    /// `GeomFill_LocationLaw::D1` (`cxx:28-37`): `Standard_NotImplemented`.
    fn d1(&self, _param: f64) -> Result<LocationD1, String> {
        Err("Standard_NotImplemented: GeomFill_LocationLaw::D1".into())
    }

    /// `GeomFill_LocationLaw::D2` (`cxx:39-55`): `Standard_NotImplemented`.
    fn d2(&self, _param: f64) -> Result<LocationD2, String> {
        Err("Standard_NotImplemented: GeomFill_LocationLaw::D2".into())
    }

    /// `GeomFill_LocationLaw::Nb2dCurves` (`cxx:57-66`): traces plus the
    /// restrictions this law provides.
    fn nb_2d_curves(&self) -> i32 {
        let mut n = self.trace_number();
        if self.has_first_restriction() {
            n += 1;
        }
        if self.has_last_restriction() {
            n += 1;
        }
        n
    }

    /// `GeomFill_LocationLaw::HasFirstRestriction` (`cxx:68-71`): false.
    fn has_first_restriction(&self) -> bool {
        false
    }

    /// `GeomFill_LocationLaw::HasLastRestriction` (`cxx:73-76`): false.
    fn has_last_restriction(&self) -> bool {
        false
    }

    /// `GeomFill_LocationLaw::TraceNumber` (`cxx:78-81`): 0.
    fn trace_number(&self) -> i32 {
        0
    }

    /// `GeomFill_LocationLaw::ErrorStatus` (`cxx:85-88`): `PipeOk`.
    fn error_status(&self) -> GeomFillPipeError {
        GeomFillPipeError::PipeOk
    }

    /// `GeomFill_LocationLaw::NbIntervals` (`hxx:109`, pure virtual).
    fn nb_intervals(&self, s: u8) -> i32;

    /// `GeomFill_LocationLaw::Intervals` (`hxx:117-118`, pure virtual).
    fn intervals(&self, s: u8) -> Vec<f64>;

    /// `GeomFill_LocationLaw::SetInterval` (`hxx:123`, pure virtual).
    fn set_interval(&self, first: f64, last: f64);

    /// `GeomFill_LocationLaw::GetInterval` (`hxx:127`, pure virtual).
    fn get_interval(&self) -> (f64, f64);

    /// `GeomFill_LocationLaw::GetDomain` (`hxx:133`, pure virtual).
    fn get_domain(&self) -> (f64, f64);

    /// `GeomFill_LocationLaw::Resolution` (`cxx:92-96`):
    /// `Standard_NotImplemented`.
    fn resolution(&self, _index: i32, _tol: f64) -> Result<(f64, f64), String> {
        Err("Standard_NotImplemented: GeomFill_LocationLaw::Resolution".into())
    }

    /// `GeomFill_LocationLaw::SetTolerance` (`cxx:98-101`): does nothing.
    fn set_tolerance(&self, _tol3d: f64, _tol2d: f64) {}

    /// `GeomFill_LocationLaw::GetMaximalNorm` (`hxx:149`, pure virtual).
    fn get_maximal_norm(&self) -> f64;

    /// `GeomFill_LocationLaw::GetAverageLaw` (`hxx:152-153`, pure virtual).
    fn get_average_law(&self) -> (GpMat, GpVec);

    /// `GeomFill_LocationLaw::IsTranslation` (`cxx:103-106`): false, `error`
    /// untouched.
    fn is_translation(&self, _error: &mut f64) -> bool {
        false
    }

    /// `GeomFill_LocationLaw::IsRotation` (`cxx:108-111`): false, `error`
    /// untouched.
    fn is_rotation(&self, _error: &mut f64) -> bool {
        false
    }

    /// `GeomFill_LocationLaw::Rotation` (`cxx:113-116`):
    /// `Standard_NotImplemented`.
    fn rotation(&self) -> Result<GpPnt, String> {
        Err("Standard_NotImplemented: GeomFill_SectionLaw::Rotation".into())
    }

    /// `occ::down_cast<GeomFill_LocationGuide>` on a location-law handle
    /// (`BRepFill_PipeShell.cxx:494`, `:1215`, `:1257`). The port has no RTTI on
    /// trait objects, so every law answers `None` unless it *is* a
    /// [`LocationGuide`](crate::geom_fill::LocationGuide).
    fn as_location_guide(&self) -> Option<&crate::geom_fill::LocationGuide> {
        None
    }
}
