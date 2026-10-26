//! `ProjLib_CompProjectedCurve` - the projection of a curve onto a surface,
//! cut into continuous parts.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_CompProjectedCurve.cxx`
//! (2390 lines) + `ProjLib_CompProjectedCurve.hxx` (280 lines), OCCT 8.0.0.
//!
//! The class is split by responsibility, one file per OCCT section:
//!
//! - [`jets`]: `d1`, `d2`, `d2CurvOnSurf` (`cxx:137-305`),
//! - [`bounds`]: `ExactBound`, `DichExactBound`, `InitialPoint`
//!   (`cxx:307-544`),
//! - [`splits`]: `BuildCurveSplits`, `SplitOnDirection`, `FindSplitPoint`
//!   (`cxx:2243-2390`),
//! - [`init`]: `Init` and its sequence post-processing (`cxx:657-1227`),
//! - [`eval`]: `D0` / `D1` / `D2` / `DN` / `Value` (`cxx:1551-1755`),
//! - [`intervals`]: `NbIntervals` / `Intervals` / `BuildIntervals`
//!   (`cxx:1810-2092`).
//!
//! `mySequence` holds one `NCollection_HSequence<gp_Pnt>` per continuous part
//! (`hxx:255`); the port stores `Vec<Vec<GpPnt>>`, with OCCT's 1-based
//! `Value(i)` at `[i - 1]`.
//!
//! UNPORTED sections, with the OCCT lines they occupy:
//! - `Perform` (`cxx:1229-1408`): the approximation branch needs
//!   `Approx_CurveOnSurface`, which this port does not have. The isoparametric
//!   and single-point arms are reachable only through `Perform`, so the whole
//!   method is left out until that class lands.
//! - `Trim` (`cxx:2095-2103`): needs `Adaptor2d_Curve2d`, which the port does
//!   not model.
//! - `ShallowCopy` (`cxx:630-656`): only used by `Perform`'s `UTrim` / `VTrim`
//!   path, see above.
//! - `GetSequence` (`cxx:1758-1762`) has no Rust analogue; [`CompProjectedCurve::part`]
//!   exposes one part at a time.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::kernel::geomabs::{CurveType, SurfaceType};
use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;

use crate::approx_same_parameter::{u_resolution, v_resolution};
use crate::curve::Curve;
use crate::projlib::projected_curve::surface_type;
use crate::surface::Surface;

mod bounds;
mod eval;
mod init;
mod intervals;
mod jets;
mod splits;

/// `#define FuncTol 1.e-10` (`cxx:53`): the residual tolerance handed to every
/// `ProjLib_PrjResolve::Perform` call.
pub(crate) const FUNC_TOL: f64 = 1.0e-10;

/// `ProjLib_CompProjectedCurve` (`hxx:36-280`).
pub struct CompProjectedCurve {
    /// `mySurface` (`hxx:251`).
    surface: Option<Arc<dyn Surface>>,
    /// `myCurve` (`hxx:252`).
    curve: Option<Arc<dyn Curve>>,
    /// `myNbCurves` (`hxx:253`).
    nb_curves: i32,
    /// `mySequence` (`hxx:255`).
    sequence: Vec<Vec<GpPnt>>,
    /// `myUIso` (`hxx:256`).
    ui_so: Vec<bool>,
    /// `myVIso` (`hxx:257`).
    v_iso: Vec<bool>,
    /// `mySnglPnts` (`hxx:258`).
    sngl_pnts: Vec<bool>,
    /// `myMaxDistance` (`hxx:259`).
    max_distance: Vec<f64>,
    /// `myTabInt` (`hxx:260`).
    tab_int: Option<Vec<f64>>,
    /// `myTol3d` (`hxx:261`).
    tol3d: f64,
    /// `myContinuity` (`hxx:262`), a `GeomAbs_Shape` code.
    continuity: u8,
    /// `myMaxDegree` (`hxx:263`).
    max_degree: i32,
    /// `myMaxSeg` (`hxx:264`).
    max_seg: i32,
    /// `myProj2d` (`hxx:265`).
    proj2d: bool,
    /// `myProj3d` (`hxx:266`).
    proj3d: bool,
    /// `myMaxDist` (`hxx:267`).
    max_dist: f64,
    /// `myTolU` (`hxx:268`).
    tol_u: f64,
    /// `myTolV` (`hxx:269`).
    tol_v: f64,
    /// `myResultIsPoint` (`hxx:271`).
    result_is_point: Vec<bool>,
    /// `myResult2dUApproxError` (`hxx:272`).
    result_2d_u_approx_error: Vec<f64>,
    /// `myResult2dVApproxError` (`hxx:273`).
    result_2d_v_approx_error: Vec<f64>,
    /// `myResult3dApproxError` (`hxx:274`).
    result_3d_approx_error: Vec<f64>,
    /// `myResult3dPoint` (`hxx:275`).
    result_3d_point: Vec<GpPnt>,
    /// `myResult2dPoint` (`hxx:276`).
    result_2d_point: Vec<GpPnt2d>,
    /// `myResult3dCurve` (`hxx:277`).
    result_3d_curve: Vec<Option<Arc<dyn Curve>>>,
    /// `myResult2dCurve` (`hxx:278`).
    result_2d_curve: Vec<Option<Arc<dyn Curve2d>>>,
}

impl Default for CompProjectedCurve {
    fn default() -> Self {
        Self::new()
    }
}

impl CompProjectedCurve {
    /// `ProjLib_CompProjectedCurve()` (`cxx:546-553`).
    pub fn new() -> Self {
        Self {
            surface: None,
            curve: None,
            nb_curves: 0,
            sequence: Vec::new(),
            ui_so: Vec::new(),
            v_iso: Vec::new(),
            sngl_pnts: Vec::new(),
            max_distance: Vec::new(),
            tab_int: None,
            tol3d: 1.0e-6,
            continuity: 4,
            max_degree: 14,
            max_seg: 16,
            proj2d: true,
            proj3d: false,
            max_dist: 0.0,
            tol_u: 0.0,
            tol_v: 0.0,
            result_is_point: Vec::new(),
            result_2d_u_approx_error: Vec::new(),
            result_2d_v_approx_error: Vec::new(),
            result_3d_approx_error: Vec::new(),
            result_3d_point: Vec::new(),
            result_2d_point: Vec::new(),
            result_3d_curve: Vec::new(),
            result_2d_curve: Vec::new(),
        }
    }

    /// `ProjLib_CompProjectedCurve(S, C, TolU, TolV)` (`cxx:556-577`).
    pub fn with_tolerances(
        surface: Arc<dyn Surface>,
        curve: Arc<dyn Curve>,
        tol_u: f64,
        tol_v: f64,
    ) -> Self {
        let mut this = Self::new();
        this.surface = Some(surface);
        this.curve = Some(curve);
        this.max_dist = -1.0;
        this.tol_u = tol_u;
        this.tol_v = tol_v;
        this.init();
        this
    }

    /// `ProjLib_CompProjectedCurve(S, C, TolU, TolV, MaxDist)` (`cxx:582-602`).
    pub fn with_tolerances_and_max_dist(
        surface: Arc<dyn Surface>,
        curve: Arc<dyn Curve>,
        tol_u: f64,
        tol_v: f64,
        max_dist: f64,
    ) -> Self {
        let mut this = Self::new();
        this.surface = Some(surface);
        this.curve = Some(curve);
        this.max_dist = max_dist;
        this.tol_u = tol_u;
        this.tol_v = tol_v;
        this.init();
        this
    }

    /// `ProjLib_CompProjectedCurve(Tol3d, S, C, MaxDist)` (`cxx:606-626`):
    /// the parameter tolerances come from the surface resolutions.
    pub fn with_tol3d(
        tol3d: f64,
        surface: Arc<dyn Surface>,
        curve: Arc<dyn Curve>,
        max_dist: f64,
    ) -> Self {
        let mut this = Self::new();
        this.surface = Some(surface);
        this.curve = Some(curve);
        this.tol3d = tol3d;
        this.max_dist = max_dist;
        let s = this.surface.as_ref().expect("surface");
        this.tol_u = PCONFUSION.max(u_resolution(s.as_ref(), tol3d));
        this.tol_v = PCONFUSION.max(v_resolution(s.as_ref(), tol3d));
        this.init();
        this
    }

    // -- Load / GetSurface / GetCurve ---------------------------------------

    /// `Load(const occ::handle<Adaptor3d_Surface>&)` (`cxx:1462-1465`).
    pub fn load_surface(&mut self, s: Arc<dyn Surface>) {
        self.surface = Some(s);
    }

    /// `Load(const occ::handle<Adaptor3d_Curve>&)` (`cxx:1468-1471`).
    pub fn load_curve(&mut self, c: Arc<dyn Curve>) {
        self.curve = Some(c);
    }

    /// `GetSurface()` (`cxx:1475-1478`).
    pub fn get_surface(&self) -> Option<&Arc<dyn Surface>> {
        self.surface.as_ref()
    }

    /// `GetCurve()` (`cxx:1483-1486`).
    pub fn get_curve(&self) -> Option<&Arc<dyn Curve>> {
        self.curve.as_ref()
    }

    // -- Set* ---------------------------------------------------------------

    /// `SetTol3d` (`cxx:1412-1415`).
    pub fn set_tol3d(&mut self, tol3d: f64) {
        self.tol3d = tol3d;
    }

    /// `SetContinuity` (`cxx:1418-1421`).
    pub fn set_continuity(&mut self, continuity: u8) {
        self.continuity = continuity;
    }

    /// `SetMaxDegree` (`cxx:1424-1436`): values below 1 are ignored.
    pub fn set_max_degree(&mut self, max_degree: i32) {
        if max_degree < 1 {
            return;
        }
        self.max_degree = max_degree;
    }

    /// `SetMaxSeg` (`cxx:1439-1450`): values below 1 are ignored.
    pub fn set_max_seg(&mut self, max_seg: i32) {
        if max_seg < 1 {
            return;
        }
        self.max_seg = max_seg;
    }

    /// `SetProj3d` (`cxx:1455-1458`).
    pub fn set_proj3d(&mut self, proj3d: bool) {
        self.proj3d = proj3d;
    }

    /// `SetProj2d` (`cxx:1461-1464`).
    pub fn set_proj2d(&mut self, proj2d: bool) {
        self.proj2d = proj2d;
    }

    /// `GetProj2d()` (`hxx:239`).
    pub fn get_proj2d(&self) -> bool {
        self.proj2d
    }

    /// `GetProj3d()` (`hxx:242`).
    pub fn get_proj3d(&self) -> bool {
        self.proj3d
    }

    /// `GetTolerance` (`cxx:1490-1494`).
    pub fn get_tolerance(&self) -> (f64, f64) {
        (self.tol_u, self.tol_v)
    }

    // -- Part queries -------------------------------------------------------

    /// `NbCurves()` (`cxx:1497-1500`).
    pub fn nb_curves(&self) -> i32 {
        self.nb_curves
    }

    /// `Bounds(Index, Udeb, Ufin)` (`cxx:1503-1513`). `None` is
    /// `Standard_NoSuchObject` for an out-of-range index.
    pub fn bounds(&self, index: i32) -> Option<(f64, f64)> {
        let part = self.part(index)?;
        Some((part[0].x(), part[part.len() - 1].x()))
    }

    /// `IsSinglePnt(Index, P)` (`cxx:1516-1524`): `(point, is_single)`.
    pub fn is_single_pnt(&self, index: i32) -> Option<(GpPnt2d, bool)> {
        let part = self.part(index)?;
        Some((
            GpPnt2d::new(part[0].y(), part[0].z()),
            self.sngl_pnts[index as usize - 1],
        ))
    }

    /// `IsUIso(Index, U)` (`cxx:1527-1535`): `(u, is_uiso)`.
    pub fn is_u_iso(&self, index: i32) -> Option<(f64, bool)> {
        let part = self.part(index)?;
        Some((part[0].y(), self.ui_so[index as usize - 1]))
    }

    /// `IsVIso(Index, V)` (`cxx:1541-1549`): `(v, is_viso)`.
    pub fn is_v_iso(&self, index: i32) -> Option<(f64, bool)> {
        let part = self.part(index)?;
        Some((part[0].z(), self.v_iso[index as usize - 1]))
    }

    /// `GetSequence()->Value(Index)`): one continuous part as `(t, u, v)`
    /// triples. `None` is `Standard_NoSuchObject`.
    pub fn part(&self, index: i32) -> Option<&[GpPnt]> {
        if index < 1 || index > self.nb_curves {
            return None;
        }
        self.sequence.get(index as usize - 1).map(|p| p.as_slice())
    }

    /// `MaxDistance(Index)` (`cxx:1801-1807`).
    pub fn max_distance(&self, index: i32) -> Option<f64> {
        if index < 1 || index > self.nb_curves {
            return None;
        }
        self.max_distance.get(index as usize - 1).copied()
    }

    /// `FirstParameter()` (`cxx:1766-1769`).
    pub fn first_parameter(&self) -> f64 {
        self.curve.as_ref().map(|c| c.first_parameter()).unwrap_or(0.0)
    }

    /// `LastParameter()` (`cxx:1773-1776`).
    pub fn last_parameter(&self) -> f64 {
        self.curve.as_ref().map(|c| c.last_parameter()).unwrap_or(0.0)
    }

    /// `Continuity()` (`cxx:1781-1796`): the minimum of the curve's and the
    /// surface's two continuities.
    pub fn continuity(&self) -> u8 {
        let (Some(c), Some(s)) = (self.curve.as_ref(), self.surface.as_ref()) else {
            return self.continuity;
        };
        let mut cont = c.continuity();
        cont = cont.min(s.continuity());
        cont
    }

    /// `GetType()` (`cxx:2107-2110`): always `GeomAbs_OtherCurve`.
    pub fn get_type(&self) -> CurveType {
        CurveType::OtherCurve
    }

    /// `ResultIsPoint(theIndex)` (`cxx:2113-2116`).
    pub fn result_is_point(&self, index: i32) -> Option<bool> {
        self.result_is_point.get(index as usize - 1).copied()
    }

    /// `GetResult2dUApproxError(theIndex)` (`cxx:2120-2123`).
    pub fn get_result_2d_u_approx_error(&self, index: i32) -> Option<f64> {
        self.result_2d_u_approx_error.get(index as usize - 1).copied()
    }

    /// `GetResult2dVApproxError(theIndex)` (`cxx:2128-2131`).
    pub fn get_result_2d_v_approx_error(&self, index: i32) -> Option<f64> {
        self.result_2d_v_approx_error.get(index as usize - 1).copied()
    }

    /// `GetResult3dApproxError(theIndex)` (`cxx:2136-2139`).
    pub fn get_result_3d_approx_error(&self, index: i32) -> Option<f64> {
        self.result_3d_approx_error.get(index as usize - 1).copied()
    }

    /// `GetResult2dC(theIndex)` (`cxx:2144-2147`).
    pub fn get_result_2d_c(&self, index: i32) -> Option<&Option<Arc<dyn Curve2d>>> {
        self.result_2d_curve.get(index as usize - 1)
    }

    /// `GetResult3dC(theIndex)` (`cxx:2151-2154`).
    pub fn get_result_3d_c(&self, index: i32) -> Option<&Option<Arc<dyn Curve>>> {
        self.result_3d_curve.get(index as usize - 1)
    }

    /// `GetResult2dP(theIndex)` (`cxx:2158-2163`). `Err` is OCCT's
    /// `Standard_TypeMismatch` when the part is not punctual.
    pub fn get_result_2d_p(&self, index: i32) -> Result<GpPnt2d, String> {
        if !self.result_is_point(index).unwrap_or(false) {
            return Err("ProjLib_CompProjectedCurve: result is not a point 2d".into());
        }
        self.result_2d_point
            .get(index as usize - 1)
            .copied()
            .ok_or_else(|| "ProjLib_CompProjectedCurve: index out of range".to_string())
    }

    /// `GetResult3dP(theIndex)` (`cxx:2166-2171`). `Err` is OCCT's
    /// `Standard_TypeMismatch` when the part is not punctual.
    pub fn get_result_3d_p(&self, index: i32) -> Result<GpPnt, String> {
        if !self.result_is_point(index).unwrap_or(false) {
            return Err("ProjLib_CompProjectedCurve: result is not a point 3d".into());
        }
        self.result_3d_point
            .get(index as usize - 1)
            .copied()
            .ok_or_else(|| "ProjLib_CompProjectedCurve: index out of range".to_string())
    }

    /// `UpdateTripleByTrapCriteria(gp_Pnt&)` (`cxx:2176-2237`).
    ///
    /// Escapes two known traps: a projection stuck on the `v` border of a
    /// surface of revolution (bug 25892, `cxx:2182-2192`) and one stuck on a
    /// degenerated edge of a sphere (bug 27135, `cxx:2195-2202`).
    pub(crate) fn update_triple_by_trap_criteria(&self, point: &mut GpPnt) {
        let (Some(surface), Some(curve)) = (self.surface.as_ref(), self.curve.as_ref()) else {
            return;
        };
        let s = surface.as_ref();
        let mut is_problems_possible = false;

        if s.is_surface_of_revolution() {
            // Compute maximal deviation from 3D and choose the biggest one.
            let a_v_res = v_resolution(s, occt_core::precision::CONFUSION);
            let a_max_tol = PCONFUSION.max(a_v_res);
            if (point.z() - s.v_range().0).abs() < a_max_tol
                || (point.z() - s.v_range().1).abs() < a_max_tol
            {
                is_problems_possible = true;
            }
        }

        if surface_type(s) == SurfaceType::Sphere
            && ((point.z() - s.v_range().0).abs() < PCONFUSION
                || (point.z() - s.v_range().1).abs() < PCONFUSION
                || (point.y() - s.u_range().0).abs() < PCONFUSION
                || (point.y() - s.u_range().1).abs() < PCONFUSION)
        {
            is_problems_possible = true;
        }

        if !is_problems_possible {
            return;
        }

        let Some((mut u, mut v)) = bounds::initial_point(
            &curve.d0(point.x()),
            point.x(),
            curve.as_ref(),
            s,
            PCONFUSION,
            PCONFUSION,
            self.max_dist,
        ) else {
            return;
        };

        // Restore the original position in case of a period jump.
        if s.is_u_periodic() && ((u - point.y()).abs() - s.u_period()).abs() < PCONFUSION {
            u = point.y();
        }
        if s.is_v_periodic() && ((v - point.z()).abs() - s.v_period()).abs() < PCONFUSION {
            v = point.z();
        }
        point.set_y(u);
        point.set_z(v);
    }
}
