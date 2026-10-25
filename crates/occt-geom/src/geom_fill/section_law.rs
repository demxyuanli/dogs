//! `GeomFill_SectionLaw` — the abstract law that gives the swept section at a
//! parameter.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_SectionLaw.hxx:33-164` and
//! `GeomFill_SectionLaw.cxx:26-91` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `NCollection_Array1<gp_Pnt>` / `<gp_Vec>` / `<double>` become slices; the
//!   caller sizes them from `SectionShape` (`hxx:69`) plus `IsRational`.
//! - `occ::handle<Geom_BSplineSurface>` is `Arc<GeomBSplineSurface>`;
//!   `BSplineSurface` returning a null handle (`cxx:47-52`) becomes `None`.
//! - The base `D1` / `D2` throw `Standard_NotImplemented` (`cxx:27-44`);
//!   `BarycentreOfSurf` / `GetMinimalWeight` throw it too (`cxx:58-67`), while
//!   `ConstantSection` / `CirclSection` throw `Standard_DomainError`
//!   (`cxx:74-91`). All become `Err` with the OCCT class name.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};

use crate::bspline_surface::GeomBSplineSurface;
use crate::curve::Curve;

/// `GeomFill_SectionLaw` (`GeomFill_SectionLaw.hxx:33-164`).
pub trait SectionLaw: Send + Sync {
    /// `GeomFill_SectionLaw::D0` (`hxx:38-41`, pure virtual): poles and
    /// weights at parameter `param`.
    fn d0(&self, param: f64) -> Result<(Vec<GpPnt>, Vec<f64>), String>;

    /// `GeomFill_SectionLaw::D1` (`cxx:27-34`): `Standard_NotImplemented`.
    fn d1(&self, _param: f64) -> Result<(Vec<GpPnt>, Vec<GpVec>, Vec<f64>, Vec<f64>), String> {
        Err("Standard_NotImplemented: GeomFill_SectionLaw::D1".into())
    }

    /// `GeomFill_SectionLaw::D2` (`cxx:36-44`): `Standard_NotImplemented`.
    fn d2(
        &self,
        _param: f64,
    ) -> Result<
        (
            Vec<GpPnt>,
            Vec<GpVec>,
            Vec<GpVec>,
            Vec<f64>,
            Vec<f64>,
            Vec<f64>,
        ),
        String,
    > {
        Err("Standard_NotImplemented: GeomFill_SectionLaw::D2".into())
    }

    /// `GeomFill_SectionLaw::BSplineSurface` (`cxx:46-52`): null by default.
    fn bspline_surface(&self) -> Option<Arc<GeomBSplineSurface>> {
        None
    }

    /// `GeomFill_SectionLaw::SectionShape` (`hxx:69`, pure virtual):
    /// `(nb_poles, nb_knots, degree)`.
    fn section_shape(&self) -> (i32, i32, i32);

    /// `GeomFill_SectionLaw::Knots` (`hxx:72`, pure virtual).
    fn knots(&self) -> Vec<f64>;

    /// `GeomFill_SectionLaw::Mults` (`hxx:75`, pure virtual).
    fn mults(&self) -> Vec<i32>;

    /// `GeomFill_SectionLaw::IsRational` (`hxx:78`, pure virtual).
    fn is_rational(&self) -> bool;

    /// `GeomFill_SectionLaw::IsUPeriodic` (`hxx:81`, pure virtual).
    fn is_u_periodic(&self) -> bool;

    /// `GeomFill_SectionLaw::IsVPeriodic` (`hxx:84`, pure virtual).
    fn is_v_periodic(&self) -> bool;

    /// `GeomFill_SectionLaw::NbIntervals` (`hxx:89`, pure virtual).
    fn nb_intervals(&self, s: u8) -> i32;

    /// `GeomFill_SectionLaw::Intervals` (`hxx:96-97`, pure virtual).
    fn intervals(&self, s: u8) -> Vec<f64>;

    /// `GeomFill_SectionLaw::SetInterval` (`hxx:102`, pure virtual).
    fn set_interval(&self, first: f64, last: f64);

    /// `GeomFill_SectionLaw::GetInterval` (`hxx:106`, pure virtual).
    fn get_interval(&self) -> (f64, f64);

    /// `GeomFill_SectionLaw::GetDomain` (`hxx:112`, pure virtual).
    fn get_domain(&self) -> (f64, f64);

    /// `GeomFill_SectionLaw::GetTolerance` (`hxx:119-122`, pure virtual).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64) -> Vec<f64>;

    /// `GeomFill_SectionLaw::SetTolerance` (`cxx:50-53`): does nothing.
    fn set_tolerance(&self, _tol3d: f64, _tol2d: f64) {}

    /// `GeomFill_SectionLaw::BarycentreOfSurf` (`cxx:55-61`):
    /// `Standard_NotImplemented`.
    fn barycentre_of_surf(&self) -> Result<GpPnt, String> {
        Err("Standard_NotImplemented: GeomFill_SectionLaw::BarycentreOfSurf".into())
    }

    /// `GeomFill_SectionLaw::MaximalSection` (`hxx:139`, pure virtual).
    fn maximal_section(&self) -> f64;

    /// `GeomFill_SectionLaw::GetMinimalWeight` (`cxx:63-67`):
    /// `Standard_NotImplemented`.
    fn get_minimal_weight(&self) -> Result<Vec<f64>, String> {
        Err("Standard_NotImplemented: GeomFill_SectionLaw::GetMinimalWeight".into())
    }

    /// `GeomFill_SectionLaw::IsConstant` (`cxx:69-73`): sets `error = 0` and
    /// returns false.
    fn is_constant(&self, error: &mut f64) -> bool {
        *error = 0.0;
        false
    }

    /// `GeomFill_SectionLaw::ConstantSection` (`cxx:75-80`):
    /// `Standard_DomainError`.
    fn constant_section(&self) -> Result<Arc<dyn Curve>, String> {
        Err("Standard_DomainError: GeomFill_SectionLaw::ConstantSection".into())
    }

    /// `GeomFill_SectionLaw::IsConicalLaw` (`cxx:82-86`): sets `error = 0` and
    /// returns false.
    fn is_conical_law(&self, error: &mut f64) -> bool {
        *error = 0.0;
        false
    }

    /// `GeomFill_SectionLaw::CirclSection` (`cxx:88-91`):
    /// `Standard_DomainError`.
    fn circl_section(&self, _param: f64) -> Result<Arc<dyn Curve>, String> {
        Err("Standard_DomainError: GeomFill_SectionLaw::CirclSection".into())
    }
}
